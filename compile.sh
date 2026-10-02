#!/usr/bin/env bash
set -euo pipefail

# Work from the repository root, even when invoked from another directory.
cd "$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
USER_SUBDIR="${1:-}"
while [[ -z "${USER_SUBDIR// }" ]]; do
  read -r -p "[?] Enter the folder name you prepared under ./tests folder: " USER_SUBDIR
  if [[ -n "${USER_SUBDIR// }" ]]; then
    break
  fi
  echo "[-] Folder name cannot be empty. Please try again."
done

BASE_DIR="./tests/${USER_SUBDIR}"
if [[ ! -d "$BASE_DIR" ]]; then
  echo "Folder not found: $BASE_DIR"
  echo "    Please create it and rerun."
  exit 1
fi

INFILE="${BASE_DIR}/reaction_network.xlsx"
OUTFILE="${BASE_DIR}/reaction_network.ntw"

echo "[*] Using base folder : $BASE_DIR"
echo "[*] Input Excel       : $INFILE"
echo "[*] Output NTW        : $OUTFILE"

# -------------------------------
# ------ Dependency check -------
# -------------------------------
need_cmd() { command -v "$1" >/dev/null 2>&1; }

echo "[*] Checking python3 ..."
if ! need_cmd python3; then
  echo "python3 not found. Please install it: sudo apt install -y python3"
  exit 1
fi

# Do not install or change the user's Python environment automatically.
python3 -c "import pandas, openpyxl, numpy" || {
  echo "Missing Python dependency: install pandas openpyxl numpy." >&2
  exit 1
}

# -------------------------------
# Reaction Network Python Converter
# -------------------------------
python3 - "$INFILE" "$OUTFILE" <<'PY'
# -*- coding: utf-8 -*-
import sys
from pathlib import Path
import pandas as pd

EXPECTED_HEADER_LEADING = "@format:idx"
INTERNAL_COLUMNS = ["IDX","R1","R2","R3","P1","P2","P3","RATE"]
EXPORT_HEADER = ["@format:idx","R","R","R","P","P","P","rate"]
HEADER_COMMENT = "#@var:T=Tgas"

def find_columns(df):
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    if EXPECTED_HEADER_LEADING not in df.columns:
        raise ValueError(f"Required column '{EXPECTED_HEADER_LEADING}' not found.")
    idx_col = EXPECTED_HEADER_LEADING
    r_cols = [c for c in df.columns if isinstance(c,str) and (c=="R" or c.startswith("R."))]
    p_cols = [c for c in df.columns if isinstance(c,str) and (c=="P" or c.startswith("P."))]
    if "k" in df.columns:
        rate_col = "k"
    elif "rate" in df.columns:
        rate_col = "rate"
    else:
        raise ValueError("Rate column 'k' or 'rate' not found.")
    r_cols = (r_cols + [None, None, None])[:3]
    p_cols = (p_cols + [None, None, None])[:3]
    return idx_col, r_cols, p_cols, rate_col

def s(series):
    series = series.astype("string").fillna("")
    return series.apply(lambda x: x.strip() if isinstance(x,str) else x)

def build_output_df(df, idx_col, r_cols, p_cols, rate_col):
    out = pd.DataFrame()
    out["IDX"]  = s(df[idx_col])
    out["R1"]   = s(df[r_cols[0]]) if r_cols[0] else ""
    out["R2"]   = s(df[r_cols[1]]) if r_cols[1] else ""
    out["R3"]   = s(df[r_cols[2]]) if r_cols[2] else ""
    out["P1"]   = s(df[p_cols[0]]) if p_cols[0] else ""
    out["P2"]   = s(df[p_cols[1]]) if p_cols[1] else ""
    out["P3"]   = s(df[p_cols[2]]) if p_cols[2] else ""
    out["RATE"] = s(df[rate_col])
    return out[INTERNAL_COLUMNS]

def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: python3 - <in.xlsx> <out.ntw>")
    in_path = Path(sys.argv[1]).expanduser().resolve()
    out_path = Path(sys.argv[2]).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")
    df = pd.read_excel(in_path, dtype="string", engine="openpyxl")
    if df.empty:
        raise SystemExit("Excel file is empty.")
    idx_col, r_cols, p_cols, rate_col = find_columns(df)
    out_df = build_output_df(df, idx_col, r_cols, p_cols, rate_col)
    tmp_path = out_path.with_suffix(".tmp")
    try:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", line_terminator="\n")
    except TypeError:
        out_df.to_csv(tmp_path, index=False, header=EXPORT_HEADER, encoding="utf-8", lineterminator="\n")
    with open(out_path, "w", encoding="utf-8") as f_out, open(tmp_path, "r", encoding="utf-8") as f_in:
        f_out.write(f"{HEADER_COMMENT}\n")
        f_out.writelines(f_in.readlines())
    tmp_path.unlink(missing_ok=True)
    print(f"File generated successfully: {out_path}")

if __name__ == "__main__":
    main()
PY

echo "Reaction network conversion completed successfully."

# Keep the case-owned manifest and driver intact across regeneration.
if [[ ! -f "${BASE_DIR}/copylist.pcp" ]]; then
  echo "Missing copylist.pcp in $BASE_DIR" >&2
  exit 1
fi

# -------------------------------
# Convert settings.xlsx → options.opt
# -------------------------------
SETTINGS_XLSX="${BASE_DIR}/settings.xlsx"
OPTIONS_OPT="${BASE_DIR}/options.opt"

echo "[*] Converting settings.xlsx to options.opt ..."
python3 - "$SETTINGS_XLSX" "$OPTIONS_OPT" "$BASE_DIR" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import os

if len(sys.argv) < 4:
    raise SystemExit("Usage: python3 - <settings.xlsx> <options.opt> <BASE_DIR>")

in_path = Path(sys.argv[1]).expanduser().resolve()
out_path = Path(sys.argv[2]).expanduser().resolve()
base_dir = Path(sys.argv[3]).expanduser().resolve()

if not in_path.exists():
    raise SystemExit(f"settings.xlsx not found: {in_path}")

df = pd.read_excel(in_path, dtype="string", engine="openpyxl")
expected = ["Parameter", "unit", "input"]
if not all(col in df.columns for col in expected):
    raise SystemExit("settings.xlsx must have columns: Parameter, unit, input")
df = df.fillna("")
relative_path = os.path.relpath(base_dir / "reaction_network.ntw", Path.cwd()).replace("./", "")
with open(out_path, "w", encoding="utf-8") as f:
    f.write(f"network = {relative_path}\n")
    for _, row in df.iterrows():
        param = str(row["Parameter"]).strip()
        value = str(row["input"]).strip() if isinstance(row["input"], str) else ""
        f.write(f"{param} = {value}\n")
print(f"options.opt created at: {out_path}")
PY
echo "options.opt created successfully."

# -------------------------------
# Convert profile.xlsx → profile.dat
# -------------------------------
PROFILE_XLSX="${BASE_DIR}/profile.xlsx"
PROFILE_DAT="${BASE_DIR}/profile.dat"

echo "[*] Converting profile.xlsx to profile.dat ..."
python3 - "$PROFILE_XLSX" "$PROFILE_DAT" <<'PY'
import sys
from pathlib import Path
import pandas as pd
import numpy as np

if len(sys.argv) < 3:
    raise SystemExit("Usage: python3 - <profile.xlsx> <profile.dat>")

in_path = Path(sys.argv[1]).expanduser().resolve()
out_path = Path(sys.argv[2]).expanduser().resolve()

if not in_path.exists():
    raise SystemExit(f"profile.xlsx not found: {in_path}")

df = pd.read_excel(in_path, engine="openpyxl")
df = df.where(pd.notnull(df), "")

cols_keep = [0, 1]
for i, col in enumerate(df.columns):
    if i not in cols_keep:
        df[col] = df[col].apply(lambda x: f"{float(x):.4E}" if isinstance(x,(int,float,np.integer,np.floating)) else x)

total_columns = len(df.columns)
env_params = 5
species_count = total_columns - env_params

tmp_path = out_path.with_suffix(".tmp")
try:
    df.to_csv(tmp_path, sep="\t", index=False, header=True, encoding="utf-8", line_terminator="\n")
except TypeError:
    df.to_csv(tmp_path, sep="\t", index=False, header=True, encoding="utf-8", lineterminator="\n")

with open(out_path, "w", encoding="utf-8") as fout:
    fout.write(f"{env_params}\t{species_count}\n")
    with open(tmp_path, "r", encoding="utf-8") as fin:
        fout.writelines(fin.readlines())
tmp_path.unlink(missing_ok=True)
print(f"profile.dat created with header line: {env_params} {species_count}")
PY
echo "profile.dat created successfully."

# Solar conversion is part of the generator and uses its actual bin edges.
echo "[*] Generating PATMO case: ${USER_SUBDIR}"
python3 patmo -test="${USER_SUBDIR}"
echo "Generated build/. Next: cd build && make"
