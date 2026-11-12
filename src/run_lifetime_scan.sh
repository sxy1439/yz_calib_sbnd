#!/bin/bash
# Produce YZ maps with yz_median_planes_sce.C for several electron lifetimes.
# For each lifetime: set the default getELifetime() value in file_handling.cpp,
# rebuild libFileHandling.so (make clean && make), run the macro, and rename the output.
# file_handling.cpp is restored and the library rebuilt when the script exits.
#
# Run from anywhere on the gpvm after setting up ROOT:  ./run_lifetime_scan.sh
# Uses GNU sed (Linux); on macOS change "sed -i" to "sed -i ''".
set -euo pipefail
cd "$(dirname "$0")"

CINTYP="doData_run_2_gen_2_BNBLight"
OUTDIR="../output_files/run_2_sample"
BASE="run_2_gen_2_BNBLight_data_v10_14_02_02"
DEFAULT_OUT="$OUTDIR/$BASE.root"   # where getOutputNameYZ() writes for $CINTYP
TAG="// lifetime scan"

# label                 East [ms]  West [ms]
RUNS=(
  "30usEast30usWest      0.03       0.03"
  "100usEast100usWest    0.1        0.1"
  "200usEast80usWest     0.2        0.08"
  "1000usEast1000usWest  1.0        1.0"
)

if ! grep -q "$TAG" file_handling.cpp; then
  echo "ERROR: no line tagged '$TAG' in file_handling.cpp" >&2
  exit 1
fi
if [ -e "$DEFAULT_OUT" ]; then
  echo "ERROR: $DEFAULT_OUT already exists and would be overwritten. Move it first." >&2
  exit 1
fi

cp file_handling.cpp file_handling.cpp.orig
restore() {
  mv file_handling.cpp.orig file_handling.cpp
  make clean >/dev/null && make >/dev/null
  echo "Restored file_handling.cpp and rebuilt libFileHandling.so"
}
trap restore EXIT

for run in "${RUNS[@]}"; do
  read -r label east west <<< "$run"
  out="$OUTDIR/CORRECTVERSION_${BASE}(${label}).root"
  log="$OUTDIR/CORRECTVERSION_${BASE}(${label}).log"

  if [ -e "$out" ]; then
    echo "Skipping $label: $out already exists"
    continue
  fi

  echo "=== $label: East = $east ms, West = $west ms"
  sed -i -E "\|$TAG|s|\{[^}]*\}|{${east}, ${west}}|" file_handling.cpp
  grep -F "{${east}, ${west}};$TAG" file_handling.cpp >/dev/null || {
    echo "ERROR: failed to set lifetime in file_handling.cpp" >&2
    exit 1
  }

  make clean && make
  root -l -b -q "yz_median_planes_sce.C(\"$CINTYP\")" 2>&1 | tee "$log"

  if [ ! -e "$DEFAULT_OUT" ]; then
    echo "ERROR: macro did not produce $DEFAULT_OUT (see $log)" >&2
    exit 1
  fi
  mv "$DEFAULT_OUT" "$out"
  echo "Wrote $out"
done
