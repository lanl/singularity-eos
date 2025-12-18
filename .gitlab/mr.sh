export PROJECT_SPACK_ENV_VERSION_DEFAULT=$(printf '%s\n' "${CI_MERGE_REQUEST_DESCRIPTION}" | grep "^Deployment Version:" | cut -d: -f2)
export PROJECT_SPACK_ENV_MR=$(printf '%s\n' "${CI_MERGE_REQUEST_DESCRIPTION}" | grep "^Requires:" | cut -d: -f2- | grep "xcap/deployment/-/merge_requests" | xargs basename 2>/dev/null)
export XCAP_SPACKAGES_MR=$(printf '%s\n' "${CI_MERGE_REQUEST_DESCRIPTION}" | grep "^Requires:" | cut -d: -f2- | grep "xcap/spackages/-/merge_requests" | xargs basename 2>/dev/null)
