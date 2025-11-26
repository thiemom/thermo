#!/usr/bin/env bash
set -euo pipefail

# Simple helper to run clang-tidy only on library / core code.
# - An existing CMake build directory with compile_commands.json is required.
# - Tests under tests/ and Python bindings under python/ are intentionally skipped
#   to keep noise low and highlight issues in src/ and include/.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
COMPILE_COMMANDS="${BUILD_DIR}/compile_commands.json"

if [[ ! -f "${COMPILE_COMMANDS}" ]]; then
  echo "error: ${COMPILE_COMMANDS} not found." >&2
  echo "       Configure the project with CMake using -DCMAKE_EXPORT_COMPILE_COMMANDS=ON" >&2
  echo "       e.g.:" >&2
  echo "         mkdir -p build && cd build" >&2
  echo "         cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON .." >&2
  exit 1
fi

cd "${ROOT_DIR}"

# Collect sources we care about: .cpp in src/ and headers in include/.
# tests/ and python/ are intentionally not included.
mapfile -t FILES < <(git ls-files \
  "src/*.cpp" "src/**/*.cpp" \
  "include/*.h" "include/**/*.h" 2>/dev/null)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No matching source files found under src/ or include/." >&2
  exit 0
fi

# Run clang-tidy. You can pass extra arguments via CLANG_TIDY_FLAGS, e.g.:
#   CLANG_TIDY_FLAGS="-p ${BUILD_DIR}" scripts/run-clang-tidy.sh
# By default we point -p to the build dir with compile_commands.json.

CLANG_TIDY_BIN=${CLANG_TIDY_BIN:-clang-tidy}

"${CLANG_TIDY_BIN}" \
  -p "${BUILD_DIR}" \
  "${FILES[@]}"
