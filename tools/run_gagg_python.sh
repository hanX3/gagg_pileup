#!/usr/bin/env bash
# 自动准备并调用当前系统的项目 Python，参数原样转交。
set -euo pipefail
gagg_project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${GAGG_BASE_PYTHON:-/usr/bin/python3}" "$gagg_project_dir/tools/gagg_environment.py" run "$@"
