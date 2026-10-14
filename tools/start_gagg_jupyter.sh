#!/usr/bin/env bash
# 从任意目录启动经典 Notebook；自动准备并复用当前系统的本机环境。
set -euo pipefail

gagg_project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${GAGG_BASE_PYTHON:-/usr/bin/python3}" "$gagg_project_dir/tools/gagg_environment.py" start "$@"
