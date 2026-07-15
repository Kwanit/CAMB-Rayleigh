#!/usr/bin/env bash

set -euo pipefail

workspace_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
external_dir="${workspace_dir}/external"
cosmorec_dir="${COSMOREC_PATH:-${external_dir}/CosmoRec}"
hyrec_dir="${HYREC_PATH:-${external_dir}/HYREC-2}"
cosmorec_repo="${COSMOREC_REPO:-https://github.com/cmbant/CosmoRec.git}"
hyrec_repo="${HYREC_REPO:-https://github.com/nanoomlee/HYREC-2.git}"

safe_git() {
    env -u GIT_DIR -u GIT_WORK_TREE -u GIT_COMMON_DIR -u GIT_INDEX_FILE git "$@"
}

cosmorec_dir="${cosmorec_dir%/}"
hyrec_dir="${hyrec_dir%/}"

mkdir -p "${external_dir}"

if [[ ! -d "${cosmorec_dir}/.git" ]]; then
    rm -rf "${cosmorec_dir}"
    safe_git clone --depth 1 --filter=blob:none "${cosmorec_repo}" "${cosmorec_dir}"
fi

chmod -R u+rwX "${cosmorec_dir}" >/dev/null 2>&1 || true

if [[ ! -f "${cosmorec_dir}/Makefile" ]]; then
    echo "CosmoRec checkout at ${cosmorec_dir} does not contain a Makefile" >&2
    exit 1
fi

if [[ ! -d "${hyrec_dir}/.git" ]]; then
    rm -rf "${hyrec_dir}"
    safe_git clone --depth 1 --filter=blob:none "${hyrec_repo}" "${hyrec_dir}"
fi

chmod -R u+rwX "${hyrec_dir}" >/dev/null 2>&1 || true

if [[ ! -f "${hyrec_dir}/Makefile" ]]; then
    echo "HYREC-2 checkout at ${hyrec_dir} does not contain a Makefile" >&2
    exit 1
fi
