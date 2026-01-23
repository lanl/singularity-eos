export DEPLOYMENT_VERSION_CURRENT_DEFAULT="2026-01-07"
export DEPLOYMENT_VERSION_DEFAULT="${DEPLOYMENT_VERSION_DEFAULT:-$DEPLOYMENT_VERSION_CURRENT_DEFAULT}"
SCRIPT_PATH=${BASH_SOURCE[0]:-${(%):-%x}}
PARENT_DIR="$( cd "$( dirname "${SCRIPT_PATH}" )" &>/dev/null && pwd )"

if [ -n "$DEPLOYMENT_MR" ]; then
  export DEPLOYMENT_VERSION=${DEPLOYMENT_VERSION:-mr/$DEPLOYMENT_MR/$DEPLOYMENT_VERSION_DEFAULT}
else
  export DEPLOYMENT_VERSION=${DEPLOYMENT_VERSION:-$DEPLOYMENT_VERSION_DEFAULT}
fi

export KESSEL_INIT="source $SCRIPT_PATH"
export SINGULARITY_EOS_CHECKOUT=$(realpath $PARENT_DIR/..)

if command -v jq >/dev/null 2>&1 && command -v sacctmgr >/dev/null 2>&1; then
  SYSTEM_NAME=$(sacctmgr list --json clusters  | jq -r '.clusters[0].name')
elif command -v flux >/dev/null 2>&1; then
  SYSTEM_NAME=$(hostname | sed 's/[0-9]//g')
fi

# During a regular CI run the workflow deployment is always a temporary copy of
# the cluster deployment to allow customization. This can be overwritten by
# setting KESSEL_WORKFLOW_DEPLOYMENT to another persistent location. If the
# folder set in KESSEL_WORKFLOW_DEPLOYMENT doesn't exist, a copy of the cluster
# deployment will be made. If KESSEL_WORKFLOW_DEPLOYMENT set to "upstream", the cluster
# deployment is used, but typically remains read-only.
_KESSEL_WORKFLOW_DEPLOYMENT="$KESSEL_WORKFLOW_DEPLOYMENT"
export KESSEL_WORKFLOW_DEPLOYMENT="${KESSEL_WORKFLOW_DEPLOYMENT:-${TMPDIR:-/tmp}/$USER-ci-envs}"

if [ "$SYSTEM_NAME" == "darwin" ] || [ "$SYSTEM_NAME" == "rocinante" ]; then
  export KESSEL_DEPLOYMENT=${KESSEL_DEPLOYMENT:-/usr/projects/xcap/oss/deployments/$DEPLOYMENT_VERSION/$SYSTEM_NAME}
elif [ "$SYSTEM_NAME" == "rzadams" ] || [ "$SYSTEM_NAME" == "rzvernal" ] || [ "$SYSTEM_NAME" == "elcapitan" ] || [ "$SYSTEM_NAME" == "tuolumne" ]; then
  export KESSEL_DEPLOYMENT=${KESSEL_DEPLOYMENT:-/usr/workspace/xcap/oss/deployments/$DEPLOYMENT_VERSION/$SYSTEM_NAME}
else
  echo "ERROR: Unknown system!" >&2
  return 1
fi

if [ "$KESSEL_WORKFLOW_DEPLOYMENT" = "upstream" ] && [ -d "$KESSEL_DEPLOYMENT" ]; then
  source "$KESSEL_DEPLOYMENT/activate.sh"
else
  if [ -d "$KESSEL_DEPLOYMENT" ] && [ ! -d "$KESSEL_WORKFLOW_DEPLOYMENT" ]; then
    source "$KESSEL_DEPLOYMENT/activate.sh"
    clone-deployment "$KESSEL_WORKFLOW_DEPLOYMENT"
  fi
  if [ "$(uname -s)" = "Darwin" ]; then
    _KESSEL_WORKFLOW_DEPLOYMENT_OWNER=$(stat -f %u "$KESSEL_WORKFLOW_DEPLOYMENT")
  else
    _KESSEL_WORKFLOW_DEPLOYMENT_OWNER=$(stat -c %u "$KESSEL_WORKFLOW_DEPLOYMENT")
  fi
  if [ -z "$_KESSEL_WORKFLOW_DEPLOYMENT" ] && [ "$_KESSEL_WORKFLOW_DEPLOYMENT_OWNER" -ne "$(id -u)" ]; then
    unset _KESSEL_WORKFLOW_DEPLOYMENT_OWNER
    unset _KESSEL_WORKFLOW_DEPLOYMENT
    echo "ERROR: $KESSEL_WORKFLOW_DEPLOYMENT not owned by $USER!" >&2
    return 1
  else
    unset _KESSEL_WORKFLOW_DEPLOYMENT_OWNER
    source "$KESSEL_WORKFLOW_DEPLOYMENT/activate.sh"
  fi
fi

unset _KESSEL_WORKFLOW_DEPLOYMENT
