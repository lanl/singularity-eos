SINGULARITY_GOLDFILES_VERSION="goldfiles-1.8.0"
wget https://github.com/lanl/singularity-eos/releases/download/${SINGULARITY_GOLDFILES_VERSION}/goldfiles.tar.gz

CHECKOUT_ROOT=$(git rev-parse --show-toplevel)

export XCAP_SPACKAGES_CHECKOUT="${CHECKOUT_ROOT}/extern/xcap_spackages"

mkdir -p $(dirname ${XCAP_SPACKAGES_CHECKOUT})
REPO_URL=$(git remote get-url origin)

if [ ! -d "${XCAP_SPACKAGES_CHECKOUT}" ]; then
  git clone "${REPO_URL%/*/*}/spackages.git" "${XCAP_SPACKAGES_CHECKOUT}"
fi

if [ -n "${XCAP_SPACKAGES_MR}" ]; then
  git -C ${XCAP_SPACKAGES_CHECKOUT} fetch origin merge-requests/${XCAP_SPACKAGES_MR}/head:mr-${XCAP_SPACKAGES_MR}
  export XCAP_SPACKAGES_REF="${XCAP_SPACKAGES_REF:-"mr-${XCAP_SPACKAGES_MR}"}"
else
  export XCAP_SPACKAGES_REF="${XCAP_SPACKAGES_REF:-main}"
fi

git -C ${XCAP_SPACKAGES_CHECKOUT} checkout ${XCAP_SPACKAGES_REF}

