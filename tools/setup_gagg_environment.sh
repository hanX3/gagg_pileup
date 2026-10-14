#!/usr/bin/env bash
# 可选的手动安装入口；日常启动脚本会自动准备本机环境。
set -euo pipefail
gagg_project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${GAGG_BASE_PYTHON:-/usr/bin/python3}" "$gagg_project_dir/tools/gagg_environment.py" setup "$@"
