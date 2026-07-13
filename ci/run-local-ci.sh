#!/usr/bin/env bash

set -Eeuo pipefail

WORKFLOW="${WORKFLOW:-.github/workflows/pyrad_tests_base.yml}"
JOB="${JOB:-unit_tests}"
PYTHON_VERSION="${PYTHON_VERSION:-3.12}"
EVENT="${EVENT:-push}"

# A reasonably complete Ubuntu runner image. The default act image is much

# smaller and may lack tools normally available on GitHub-hosted runners.

RUNNER_IMAGE="${RUNNER_IMAGE:-catthehacker/ubuntu:act-latest}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
TEMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/pyrad-act.XXXXXX")"

cleanup()
{
exit_code=$?

```
echo
echo "Cleaning local CI resources..."

rm -rf "${TEMP_DIR}"

# --rm normally removes the act job container. This catches containers
# left behind when act or Docker is interrupted unexpectedly.
mapfile -t act_containers < <(
    docker ps -aq \
        --filter "name=act-${JOB}" \
        --filter "label=com.github.nektos.act" 2>/dev/null || true
)

if ((${#act_containers[@]} > 0)); then
    docker rm -f "${act_containers[@]}" >/dev/null 2>&1 || true
fi

echo "Cleanup complete."
exit "${exit_code}"
```

}

trap cleanup EXIT INT TERM

usage()
{
cat <<EOF
Usage:
$(basename "$0") [Python version]

Examples:
$(basename "$0")
$(basename "$0") 3.11
PYTHON_VERSION=3.13 $(basename "$0")
WORKFLOW=.github/workflows/test.yml $(basename "$0") 3.12

Environment variables:
WORKFLOW         Workflow file to execute
JOB              GitHub Actions job ID
PYTHON_VERSION   Python matrix version
EVENT            Event to simulate, normally push or workflow_dispatch
RUNNER_IMAGE     Docker runner image
KEEP_ACT_IMAGE   Set to 0 to remove the downloaded runner image afterward

Optional test secrets are read from the current environment:
S3_SECRET_READ
S3_KEY_READ
S3_SECRET_WRITE
S3_KEY_WRITE
GITHUB_TOKEN
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
usage
exit 0
fi

if [[ $# -gt 0 ]]; then
PYTHON_VERSION="$1"
fi

cd "${REPO_ROOT}"

if ! command -v docker >/dev/null 2>&1; then
echo "Error: Docker is not installed or not available in PATH." >&2
exit 1
fi

if ! docker info >/dev/null 2>&1; then
echo "Error: Docker is not running or the current user cannot access it." >&2
exit 1
fi

if ! command -v act >/dev/null 2>&1; then
echo "Error: act is not installed." >&2
echo "See: https://nektosact.com/installation/" >&2
exit 1
fi

if [[ ! -f "${WORKFLOW}" ]]; then
echo "Error: workflow not found: ${WORKFLOW}" >&2
exit 1
fi

# Simulate the branch expected by:

#

# on:

# push:

# branches: [dev]

#

# Otherwise act may use the currently checked-out branch name.

EVENT_FILE="${TEMP_DIR}/event.json"

cat >"${EVENT_FILE}" <<EOF
{
"ref": "refs/heads/dev",
"repository": {
"full_name": "MeteoSwiss/pyrad"
},
"act": true
}
EOF

ACT_ARGS=(
"${EVENT}"
--workflows "${WORKFLOW}"
--job "${JOB}"
--matrix "os:ubuntu-latest"
--matrix "python-version:${PYTHON_VERSION}"
--platform "ubuntu-latest=${RUNNER_IMAGE}"
--eventpath "${EVENT_FILE}"
--rm
)

# Pass only secrets that are currently defined. Using "-s NAME" makes act

# read the value from the environment instead of writing it on the command line.

for secret_name in S3_SECRET_READ S3_KEY_READ S3_SECRET_WRITE S3_KEY_WRITE GITHUB_TOKEN; do
    if [[ -n "${!secret_name:-}" ]]; then
        ACT_ARGS+=(--secret "${secret_name}")
    fi
done

echo "Running local GitHub Actions test"
echo "  Workflow: ${WORKFLOW}"
echo "  Job:      ${JOB}"
echo "  OS:       ubuntu-latest"
echo "  Python:   ${PYTHON_VERSION}"
echo

act "${ACT_ARGS[@]}"

if [[ "${KEEP_ACT_IMAGE:-1}" == "0" ]]; then
echo "Removing runner image ${RUNNER_IMAGE}..."
docker image rm "${RUNNER_IMAGE}" >/dev/null 2>&1 || true
fi
