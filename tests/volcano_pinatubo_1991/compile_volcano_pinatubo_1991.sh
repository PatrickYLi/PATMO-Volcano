#!/usr/bin/env bash
set -euo pipefail
CASE_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(CDPATH= cd -- "$CASE_DIR/../.." && pwd)"
BUILD_DIR="$ROOT_DIR/build"

if [[ $# -ne 0 ]]; then
  echo "Usage: $0 (no arguments)" >&2
  exit 1
fi
if ! command -v make >/dev/null 2>&1; then
  echo "ERROR: make is required to compile the model." >&2
  exit 1
fi
for required in test.f90 copylist.pcp; do
  if [[ ! -f "$CASE_DIR/$required" ]]; then
    echo "ERROR: Missing case input: $CASE_DIR/$required" >&2
    exit 1
  fi
done

# Reuse the root conversion/generation workflow and its input manifest.
bash "$ROOT_DIR/compile.sh" "$(basename "$CASE_DIR")"

# The case owns the complete driver; never generate or edit its contents here.
cp "$CASE_DIR/test.f90" "$BUILD_DIR/test.f90"
cmp "$CASE_DIR/test.f90" "$BUILD_DIR/test.f90"
while read -r input || [[ -n "$input" ]]; do
  case "$input" in
    ""|\#*) continue ;;
  esac
  if ! cmp -s "$CASE_DIR/$input" "$BUILD_DIR/$input"; then
    echo "ERROR: Copied input differs from the case: $input" >&2
    exit 1
  fi
done < "$CASE_DIR/copylist.pcp"

echo "[*] Driver and runtime inputs verified. Compiling ..."
make -C "$BUILD_DIR"
echo "Build complete. No simulation was started."
printf 'Next: cd "%s"\n' "$BUILD_DIR"
echo "  ./test_volcano  # optical pre-run"
echo "  ./test          # full chemistry simulation"
