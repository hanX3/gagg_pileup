"""为共享项目选择本机 Python 环境；不在共享目录中安装虚拟环境。"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import venv


PROJECT = Path(__file__).resolve().parent.parent
KERNEL = "gagg-pileup"
DISPLAY_NAME = "GAGG pileup (本机环境)"
REQUIREMENTS = PROJECT / "requirements-analysis.txt"
CONSTRAINTS = PROJECT / "requirements-analysis.lock.txt"
MODULES = ("numpy", "matplotlib", "pandas", "scipy", "uproot", "awkward",
           "IPython", "ipykernel", "nbformat", "nbclient", "nbclassic", "jupyter_server")


def environment_dir() -> Path:
    """各发行版、Python 次版本及项目分别使用本地目录。"""
    override = os.environ.get("GAGG_ENV_DIR")
    if override:
        return Path(override).expanduser().resolve()
    data = Path(os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local/share")).expanduser()
    if not data.is_absolute():
        data = Path.home() / ".local/share"
    os_release = platform.freedesktop_os_release()
    system = "-".join((os_release.get("ID", "linux"), os_release.get("VERSION_ID", "unknown"),
                       platform.machine(), f"py{sys.version_info.major}.{sys.version_info.minor}"))
    system = re.sub(r"[^a-zA-Z0-9_.-]", "_", system)
    project_id = hashlib.sha256(str(PROJECT).encode()).hexdigest()[:12]
    return data / "gagg-pileup/envs" / project_id / system


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runtime_environment(directory: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONHOME", None)
    env["VIRTUAL_ENV"] = str(directory)
    env["GAGG_ENV_DIR"] = str(directory)
    env["GAGG_PROJECT_ROOT"] = str(PROJECT)
    env["PATH"] = str(directory / "bin") + os.pathsep + env.get("PATH", "")
    data = str(directory / "share/jupyter")
    env["JUPYTER_DATA_DIR"] = data
    env["JUPYTER_PATH"] = data + (os.pathsep + env["JUPYTER_PATH"] if env.get("JUPYTER_PATH") else "")
    for name, subdir in (("JUPYTER_CONFIG_DIR", "jupyter_config"),
                         ("JUPYTER_RUNTIME_DIR", "jupyter_runtime"),
                         ("IPYTHONDIR", "ipython"), ("MPLCONFIGDIR", "matplotlib")):
        env.setdefault(name, str(directory / "var" / subdir))
    return env


def check(directory: Path) -> tuple[bool, str]:
    python = directory / "bin/python"
    if not python.is_file():
        return False, "本机环境尚未安装"
    code = "import sys; " + "; ".join(f"import {name}" for name in MODULES)
    try:
        result = subprocess.run([str(python), "-c", code], env=runtime_environment(directory),
                                capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as error:
        return False, str(error)
    if result.returncode:
        return False, result.stderr.strip()
    marker = directory / "gagg-environment.json"
    if not marker.is_file():
        return False, "缺少本项目环境记录"
    saved = json.loads(marker.read_text())
    if saved.get("project_root") != str(PROJECT):
        return False, "本机环境属于其他项目路径"
    if saved.get("requirements_sha256") != digest(REQUIREMENTS):
        return False, "项目依赖清单已更新"
    if CONSTRAINTS.is_file() and saved.get("constraints_sha256") != digest(CONSTRAINTS):
        return False, "项目依赖版本已更新"
    kernel = directory / "share/jupyter/kernels" / KERNEL / "kernel.json"
    if not kernel.is_file() or json.loads(kernel.read_text())["argv"][0] != str(python):
        return False, "项目内核未注册或解释器路径失效"
    return True, "ready"


def setup(directory: Path, record: Path | None) -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    record = (record or PROJECT / "jupyter/run_records" /
              f"local_environment_{stamp}_{directory.name}").resolve()
    record.mkdir(parents=True, exist_ok=True)
    # 安装失败的环境可以继续安装；已失效的解释器先保留独立备份。
    if (directory / "bin/python").exists():
        probe = subprocess.run([str(directory / "bin/python"), "-c", "import encodings"],
                               capture_output=True, env=runtime_environment(directory))
        if probe.returncode:
            backup = directory.with_name(directory.name + ".backup-" + stamp)
            directory.rename(backup)
            print(f"旧本机环境已备份：{backup}", flush=True)
    if not (directory / "bin/python").exists():
        directory.parent.mkdir(parents=True, exist_ok=True)
        venv.EnvBuilder(with_pip=True, symlinks=False).create(directory)
    env = runtime_environment(directory)
    for name in ("JUPYTER_CONFIG_DIR", "JUPYTER_RUNTIME_DIR", "IPYTHONDIR", "MPLCONFIGDIR"):
        Path(env[name]).mkdir(parents=True, exist_ok=True, mode=0o700)
    python = directory / "bin/python"
    command = [str(python), "-m", "pip", "install", "--disable-pip-version-check",
               "--log", str(record / "pip-install-local.log"), "-r", str(REQUIREMENTS)]
    if CONSTRAINTS.is_file():
        command += ["-c", str(CONSTRAINTS)]
    print(f"本机环境：{directory}", flush=True)
    subprocess.run(command, env=env, check=True)
    subprocess.run([str(python), "-m", "ipykernel", "install", "--prefix", str(directory),
                    "--name", KERNEL, "--display-name", DISPLAY_NAME], env=env, check=True)
    kernel = directory / "share/jupyter/kernels" / KERNEL / "kernel.json"
    spec = json.loads(kernel.read_text())
    spec["env"] = {"GAGG_ENV_DIR": str(directory), "GAGG_PROJECT_ROOT": str(PROJECT)}
    kernel.write_text(json.dumps(spec, ensure_ascii=False, indent=2) + "\n")
    versions = subprocess.check_output([str(python), "-m", "pip", "freeze"], env=env, text=True)
    (record / "requirements.freeze.txt").write_text(versions)
    metadata = {"created_utc": datetime.now(timezone.utc).isoformat(), "project_root": str(PROJECT),
                "environment": str(directory), "base_python": sys.executable,
                "python_version": sys.version, "system": platform.freedesktop_os_release(),
                "kernel_name": KERNEL, "requirements_sha256": digest(REQUIREMENTS),
                "constraints_sha256": digest(CONSTRAINTS) if CONSTRAINTS.is_file() else None,
                "record_directory": str(record.relative_to(PROJECT)) if record.is_relative_to(PROJECT) else str(record)}
    (directory / "gagg-environment.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
    (record / "environment_after.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
    ready, reason = check(directory)
    if not ready:
        raise RuntimeError(reason)
    print(f"环境与 {KERNEL} 内核检查通过。", flush=True)


def register_kernel(directory: Path, prefix: Path) -> Path:
    """让已有的本机 Jupyter 服务器发现项目内核。"""
    source = directory / "share/jupyter/kernels" / KERNEL / "kernel.json"
    target = prefix / "share/jupyter/kernels" / KERNEL / "kernel.json"
    spec = json.loads(source.read_text())
    spec.setdefault("env", {}).update({
        "VIRTUAL_ENV": str(directory),
        "PATH": str(directory / "bin") + os.pathsep + "${PATH}",
    })
    content = json.dumps(spec, ensure_ascii=False, indent=2) + "\n"
    if not target.is_file() or target.read_text() != content:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("setup", "info", "register", "run", "start"))
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    directory = environment_dir()
    if args.action == "info":
        print(json.dumps({"project": str(PROJECT), "environment": str(directory),
                          "python": str(directory / "bin/python"), "kernel": KERNEL}, ensure_ascii=False, indent=2))
        return 0
    if args.action == "setup":
        options = argparse.ArgumentParser()
        options.add_argument("--record", type=Path)
        setup(directory, options.parse_args(args.arguments).record)
        return 0
    if args.action == "register":
        options = argparse.ArgumentParser(description="向已有本机 Jupyter 环境注册项目内核")
        options.add_argument("--prefix", type=Path, required=True)
        prefix = options.parse_args(args.arguments).prefix.expanduser().resolve()
    ready, reason = check(directory)
    if not ready:
        print(f"正在自动准备本机 GAGG 环境：{reason}\n环境目录：{directory}",
              file=sys.stderr, flush=True)
        try:
            setup(directory, None)
        except (OSError, subprocess.SubprocessError, RuntimeError) as error:
            print(f"环境准备失败：{error}\n处理上述错误后，重新运行同一条启动命令即可继续。",
                  file=sys.stderr, flush=True)
            return 1
    if args.action == "register":
        target = register_kernel(directory, prefix)
        print(f"GAGG 项目内核已就绪：{target}", flush=True)
        return 0
    env = runtime_environment(directory)
    python = str(directory / "bin/python")
    if args.action == "run":
        os.execve(python, [python, *args.arguments], env)
    os.chdir(PROJECT)
    print(f"使用本机环境：{directory}", file=sys.stderr)
    config = PROJECT / "jupyter/classic_notebook_config.json"
    if not config.is_file():
        config = PROJECT / "tools/classic_notebook_config.json"
    notebook = PROJECT / "jupyter/codex_0001_SiPM_Digitization.ipynb"
    default_url = "/notebooks/jupyter/codex_0001_SiPM_Digitization.ipynb" if notebook.is_file() else "/tree"
    command = [python, "-m", "nbclassic", f"--config={config}",
               f"--ServerApp.root_dir={PROJECT}", "--ServerApp.use_redirect_file=False",
               f"--MultiKernelManager.default_kernel_name={KERNEL}",
               f"--ServerApp.default_url={default_url}", *args.arguments]
    os.execve(python, command, env)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
