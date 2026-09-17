#!/usr/bin/env bash

# Run four BOMEX-WF advection experiments sequentially on one local GPU,
# then compare them with the archived swadvec=2 baseline.

set -Eeuo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    echo "Usage: $(basename "$0")"
    echo "Runs 2i5, limited 2i5, 2i62, and limited 2i62 sequentially."
    echo "Optional environment variables: CUDA_VISIBLE_DEVICES, MICROHH_EXECUTABLE, PYTHON_EXECUTABLE"
    exit 0
elif (( $# > 0 )); then
    echo "Unexpected argument: $1 (use --help)" >&2
    exit 2
fi

case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${case_dir}/../.." && pwd)"
base_ini="${case_dir}/bomex-wf.ini"
input_script="${case_dir}/bomex-wf_input.py"
analysis_script="${case_dir}/analyze_bomex_advection_tests.py"
microhh_executable="${MICROHH_EXECUTABLE:-${repo_root}/build-develop_combined_wf_sgs/microhh}"

if [[ -n "${PYTHON_EXECUTABLE:-}" ]]; then
    python_executable="${PYTHON_EXECUTABLE}"
elif command -v python >/dev/null 2>&1; then
    python_executable="$(command -v python)"
elif command -v python3 >/dev/null 2>&1; then
    python_executable="$(command -v python3)"
else
    echo "No Python interpreter found; activate the microhh environment first." >&2
    exit 1
fi

: "${CUDA_VISIBLE_DEVICES:=0}"
export CUDA_VISIBLE_DEVICES

test -x "${microhh_executable}" || {
    echo "MicroHH executable is missing or not executable: ${microhh_executable}" >&2
    exit 1
}
test -f "${base_ini}" || {
    echo "Missing base INI: ${base_ini}" >&2
    exit 1
}
test -f "${input_script}" || {
    echo "Missing input generator: ${input_script}" >&2
    exit 1
}
test -f "${analysis_script}" || {
    echo "Missing comparison analysis: ${analysis_script}" >&2
    exit 1
}
"${python_executable}" -c "import matplotlib, netCDF4, numpy" || {
    echo "The selected Python interpreter lacks a required analysis package." >&2
    echo "Activate the microhh environment or set PYTHON_EXECUTABLE." >&2
    exit 1
}

baseline_dir="${case_dir}/old/swadvec2"
for required in \
    bomex-wf.default.0000000.nc \
    bomex-wf.ql.0000000.nc \
    bomex-wf.qlcore.0000000.nc \
    bomex-wf_input.nc
do
    test -s "${baseline_dir}/${required}" || {
        echo "The swadvec=2 baseline is incomplete: ${baseline_dir}/${required}" >&2
        exit 1
    }
done
find "${baseline_dir}" -maxdepth 1 -type f -name 'ql_path.xy.*' -print -quit | grep -q . || {
    echo "The swadvec=2 baseline has no ql_path xy field for Figure 13." >&2
    exit 1
}

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
fi

prepare_ini() {
    local destination=$1
    local scheme=$2
    local limited=$3

    "${python_executable}" - "${base_ini}" "${destination}" "${scheme}" "${limited}" <<'PY'
import configparser
from pathlib import Path
import sys

source, destination, scheme, limited = sys.argv[1:]
config = configparser.ConfigParser(interpolation=None)
config.optionxform = str
with Path(source).open() as stream:
    config.read_file(stream)

config.set("advec", "swadvec", scheme)
config.remove_option("advec", "fluxlimit_list")
if limited == "1":
    config.set("advec", "fluxlimit_list", "thl,qt")

with Path(destination).open("w") as stream:
    config.write(stream, space_around_delimiters=False)
PY
}

run_case() {
    local run_name=$1
    local archive_name=$2
    local scheme=$3
    local limited=$4
    local run_dir="${case_dir}/old/${archive_name}"
    local completion_marker="${run_dir}/.run-complete"

    if [[ -f "${completion_marker}" ]]; then
        echo "Skipping completed run ${run_name}: ${run_dir}"
        return
    fi

    mkdir -p "${run_dir}"
    if find "${run_dir}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        echo "Refusing to overwrite incomplete/nonempty run directory: ${run_dir}" >&2
        echo "Inspect or move that directory before relaunching." >&2
        exit 1
    fi

    prepare_ini "${run_dir}/bomex-wf.ini" "${scheme}" "${limited}"
    cp -- "${input_script}" "${run_dir}/bomex-wf_input.py"
    cp -- "${case_dir}"/BOMEX-IEA10-*.txt "${run_dir}/"

    {
        echo "run_name=${run_name}"
        echo "archive_name=${archive_name}"
        echo "swadvec=${scheme}"
        echo "fluxlimit_list=$([[ "${limited}" == "1" ]] && echo 'thl,qt' || echo 'none')"
        echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
        echo "microhh_executable=${microhh_executable}"
        sha256sum "${microhh_executable}"
        git -C "${repo_root}" rev-parse HEAD
    } > "${run_dir}/run-metadata.txt"

    echo "Starting ${run_name} in ${run_dir}"
    (
        cd "${run_dir}"
        "${python_executable}" bomex-wf_input.py
        "${microhh_executable}" init bomex-wf 2>&1 | tee bomex-wf.init.out
        "${microhh_executable}" run bomex-wf 2>&1 | tee bomex-wf.out
    )

    test -s "${run_dir}/bomex-wf.default.0000000.nc"
    test -s "${run_dir}/bomex-wf.ql.0000000.nc"
    test -s "${run_dir}/bomex-wf.qlcore.0000000.nc"
    find "${run_dir}" -maxdepth 1 -type f -name 'ql_path.xy.*' -print -quit | grep -q .
    touch "${completion_marker}"
    echo "Completed ${run_name}"
}

run_case "swadvec2i5"          "swadvec2i5"    "2i5"  0
run_case "swadvec2i5limited"   "swadvec2i5lim" "2i5"  1
run_case "swadvec2i62"         "swadvec2i62"   "2i62" 0
run_case "swadvec2i62limited"  "swadvec2i62lim" "2i62" 1

echo "All four simulations completed; generating comparison plots."
cd "${case_dir}"
"${python_executable}" "${analysis_script}" --case-dir "${case_dir}"
