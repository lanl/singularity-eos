#!/bin/bash
PROJECT_NAME=singularity-eos
PROJECT_DEFAULT_BRANCH=main
PROJECT_GROUP=oss
BUILD_DIR=${BUILD_DIR:-build}
SOURCE_DIR=${CI_PROJECT_DIR:-$PWD}
PROJECT_SPACK_ENV=${PROJECT_SPACK_ENV:-/usr/projects/xcap/spack-env/oss/current}
CI_SPACK_ENV=${CI_SPACK_ENV:-$TMPDIR/$USER-ci-env}
UNTIL=${UNTIL:-install}

# colors
COLOR_BLUE='\033[1;34m'
COLOR_MAGENTA='\033[1;35m'
COLOR_CYAN='\033[1;36m'
COLOR_PLAIN='\033[0m'

section() {
  if [ -z "${CI}" ]; then
    echo -e "${3:+${COLOR_CYAN}$3${COLOR_PLAIN}}"
  else
    echo -e $'\e[0K'"section_$1:$(date +%s):$2"$'\r\e[0K'"${3:+${COLOR_CYAN}$3${COLOR_PLAIN}}"
  fi
}

print_usage() {
  echo "usage: build_and_test [-h] [-u {spack,env,cmake,build,test,install}] system environment"
  echo
  echo "This script allows you to active a Spack environment on a given system,"
  echo "configure your project's CMake based on its Spack variants, and then build, test and"
  echo "install it. Use the -u/--until option to stop the script early at a given phase."
  echo
  echo "After the 'env' phase has run, you will be in an active Spack environment."
  echo "You can then use the 'activate_build_env' command to enter a sub-shell"
  echo "with the Spack build environment for your project."
  echo
  echo "You can also make use of the following commands to manually run each phase and customize"
  echo "the environment in between:"
  echo " - prepare_spack"
  echo " - prepare_env"
  echo " - spack_cmake_configure"
  echo " - cmake_build [EXTRA_CMAKE_ARGS]"
  echo " - cmake_test [EXTRA_CTEST_ARGS]"
  echo " - cmake_install"
  echo
  echo "positional arguments:"
  echo "  system         name of the system"
  echo "  environment    name of the Spack environment"
  echo
  echo "options:"
  echo "  -u {spack,env,cmake,build,test,install}, --until {spack,env,cmake,build,test,install}"
  echo "                 run script until the given phase"
  echo "  -h, --help     show this help message and exit"
}

VALID_CMD=false
if [ $# -eq 2 ] && [[ ! "$1" =~ ^-.*  ]] && [[ ! "$2" =~ ^-.*  ]]; then
  VALID_CMD=true
elif [ $# -eq 4 ]; then
  if [ "$1" = "-u" ] || [ "$1" = "--until" ]; then
    if [[ "$2" =~ ^(spack|env|cmake|build|test|install)$ ]]; then
      UNTIL=$2
      VALID_CMD=true
      shift 2
    fi
  fi
fi

if ! $VALID_CMD; then
    print_usage
    return
fi

SYSTEM_NAME=$1
SPACK_ENV_NAME=$2

###############################################################################
# Project-specific steps
###############################################################################

#------------------------------------------------------------------------------
# Configure dependencies that were cloned as submodules
#------------------------------------------------------------------------------
setup_submodule_dependencies() {
  if [ "$(ls -A ${SOURCE_DIR}/utils/spiner)" ]; then
    spack develop --no-clone -p ${SOURCE_DIR}/utils/spiner spiner@main
  fi

  if [ "$(ls -A ${SOURCE_DIR}/utils/ports-of-call)" ]; then
    spack develop --no-clone -p ${SOURCE_DIR}/utils/ports-of-call ports-of-call@main
  fi
}

#------------------------------------------------------------------------------
# Hooks for before and after ctest
#------------------------------------------------------------------------------
pre_test_steps() {
  (
  source ${BUILD_ENV}
  pushd ${BUILD_DIR}
  if [[ -f ./sesame2spiner/sesame2spiner ]]; then
  ./sesame2spiner/sesame2spiner -s materials.sp5 ../sesame2spiner/examples/unit_tests/*.dat;
  ./sesame2spiner/sesame2spiner -s duplicates.sp5 ../sesame2spiner/examples/duplicate-test/*.dat;
  fi
  popd
  )
}

###############################################################################
# Generic steps
###############################################################################

#------------------------------------------------------------------------------
# Prepare Spack instance
#------------------------------------------------------------------------------
prepare_spack() {
  section start "prepare_spack[collapsed=true]" "Prepare Spack"
  umask 0007
  source $PROJECT_SPACK_ENV/replicate.sh $CI_SPACK_ENV
  section end "prepare_spack"
}

#------------------------------------------------------------------------------
# Setup up Spack environment
#------------------------------------------------------------------------------
prepare_env() {
  section start "prepare_env[collapsed=true]" "Create environment"
  echo "Activating ${SPACK_ENV_NAME} environment on ${SYSTEM_NAME}"
  source ${CI_SPACK_ENV}/systems/${SYSTEM_NAME}/activate.sh ${PROJECT_GROUP}/${PROJECT_NAME}/${SPACK_ENV_NAME}
  if [ -d ${SOURCE_DIR}/spack-repo ]; then
    spack repo add ${SOURCE_DIR}/spack-repo
  fi

  setup_submodule_dependencies

  spack develop -b ${BUILD_DIR} -p ${SOURCE_DIR} --no-clone ${PROJECT_NAME}@${PROJECT_DEFAULT_BRANCH}
  spack install --include-build-deps --only dependencies -j $(nproc)
  section end "prepare_env"
}

#------------------------------------------------------------------------------
# Run CMake (configured via Project Spack package and selected Spec)
#------------------------------------------------------------------------------
spack_cmake_configure() {
  section start "spack_cmake_configure[collapsed=true]" "Initial Spack CMake Configure"
  spack install --test root --include-build-deps -u cmake -v ${PROJECT_NAME}
  BUILD_ENV=${BUILD_DIR}/build_env.sh
  spack build-env --dump ${BUILD_ENV} ${PROJECT_NAME}
  section end "spack_cmake_configure"
}

#------------------------------------------------------------------------------
# Build project
#------------------------------------------------------------------------------
cmake_build() {
  section start "cmake_build[collapsed=false]" "CMake Build"
  (
  source ${BUILD_ENV}
  cmake -DCMAKE_VERBOSE_MAKEFILE=off -DCMAKE_INSTALL_PREFIX=$PWD/${BUILD_DIR}/install $@ ${BUILD_DIR}
  cmake --build ${BUILD_DIR} --parallel
  )
  section end "cmake_build"
}

#------------------------------------------------------------------------------
# Run tests
#------------------------------------------------------------------------------
cmake_test() {
  section start "testing[collapsed=false]" "Tests"
  pre_test_steps
  (
  source ${BUILD_ENV}
  export CTEST_OUTPUT_ON_FAILURE=1
  ctest --test-dir ${BUILD_DIR} --output-junit tests.xml $@
  )
  section end "testing"
}

#------------------------------------------------------------------------------
# Install project
#------------------------------------------------------------------------------
cmake_install() {
  section start install "Install"
  (
  source ${BUILD_ENV}
  cmake --build ${BUILD_DIR} --target install
  )
  section end install
}

#------------------------------------------------------------------------------
# Enter build environment shell
#------------------------------------------------------------------------------
activate_build_env() {
  spack build-env ${PROJECT_NAME} bash
}

###############################################################################
# Run script
###############################################################################
echo "Running on $(hostname)"

if [ ${CI} ]; then
  echo " "
  echo -e "${COLOR_BLUE}######################################################################${COLOR_PLAIN}"
  echo " "
  echo -e "${COLOR_BLUE}To recreate this CI run, follow these steps:${COLOR_PLAIN}"
  echo " "
  echo -e "${COLOR_BLUE}ssh ${CLUSTER}${COLOR_PLAIN}"
  echo -e "${COLOR_BLUE}cd /your/${PROJECT_NAME}/checkout${COLOR_PLAIN}"
  echo -e "${COLOR_BLUE}.gitlab/download_prereq.sh${COLOR_PLAIN}"
  echo -e "${COLOR_BLUE}salloc ${SCHEDULER_PARAMETERS}${COLOR_PLAIN}"
  echo -e "${COLOR_BLUE}source .gitlab/build_and_test.sh --until ${UNTIL} ${CLUSTER} ${SPACK_ENV_NAME}${COLOR_PLAIN}"
  echo " "
  echo -e "${COLOR_BLUE}See 'source .gitlab/build_and_test.sh -h' for more options.${COLOR_PLAIN}"
  echo " "
  echo -e "${COLOR_BLUE}######################################################################${COLOR_PLAIN}"
  echo " "
fi

prepare_spack
if [[ "${UNTIL}" == "spack" ]]; then return; fi

prepare_env
if [[ "${UNTIL}" == "env" ]]; then return; fi

spack_cmake_configure
if [[ "${UNTIL}" == "cmake" ]]; then return; fi

cmake_build
if [[ "${UNTIL}" == "build" ]]; then return; fi

cmake_test
if [[ "${UNTIL}" == "test" ]]; then return; fi

cmake_install
if [[ "${UNTIL}" == "install" ]]; then return; fi
