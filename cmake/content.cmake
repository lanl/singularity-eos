include(FetchContent)

###############################################################
# 
# For dependency management, use the `FetchContent` pattern
#
# NOTE: We seek to replicate the implementation of `FetchContent`
###############################################################

# Constructs `FetchContent_Declare` call, and records information
# for population and linking
# 
# ::Overview
#   With `cmake@3.24` and higher, the functionality of `FetchContent_` has
#   been updated to provide the capability to try a `find_package` at the
#   content population stage. This allows for the flexibility provided by 
#   `FetchContent_` to easily coexist with the dependency configuration
#   provided by `find_package`. Some benefits include:
#   - automated fetching of content if it is unavailable
#   - content population at configure. This allows for a package import
#     to utilize `find_package()` or `add_subdirectory()` operations
#   - local clones can be given priority with `FETCHCONTENT_SOURCE_DIR_<uppercaseName>`,
#     which lets developers use dependencies explicitly in-tree
#   The case where i.) the content or dependency is unavailable or ii.) it cannot
#   be fetched (due, for example, to a failure of a download step) is an 
#   unavoidable error.
#
#   As this cmake version is relatively fresh, we also want to provide a bridge
#   to downstream code that still is pinned to an earlier releases of `cmake`. 
#   To that end, these macros provide for replicating the main features introduced 
#   to `FetchContent_` in `cmake@3.24`. Future releases of `singularity` will require 
#   a minimum of `cmake@3.24`.
#
# :: Arguments 
#   pkg_name - the name of the package or content. usually the argument to `find_package(<pkg_name>)`
#   options:
#   - NO_FETCH - disables all download steps. If `find_package(<pkg_name>)` fails, produce an error
#   single_value:
#   - GIT_REPO - the URL of the git repository to clone, if necessary
#   - GIT_TAG  - the tag or commit to use 
#   - PATCH - apply a patch file
#   - NAMESPACE - prefix book-keeping variables with this
#   multi_value:
#   - COMPONENTS - specify required components of <pkg_name>. same as used in to `find_package`
#   - EXPECTED_TARGETS - a list of targets that expected at population.
#                       if these targets are not available after the population stage, an error is produced. 
#                       if not specified, will default to `<pkg_name>::<pkg_name>`
#   - DISABLE/ENABLE_OPTS - a list of cache vars used to configure content that is imported using `add_subdirectory()`
#
macro(singularityeos_content_declare pkg_name)
  set(options
    NO_FETCH
  )
  set(one_value_args
    GIT_REPO
    GIT_TAG
    PATCH
    NAMESPACE
  )
  set(multi_value_args
    COMPONENTS
    EXPECTED_TARGETS
    ENABLE_OPTS
    DISABLE_OPTS
  )

  cmake_parse_arguments(fp "${options}" "${one_value_args}" "${multi_value_args}" "${ARGN}")

  string(TOUPPER ${pkg_name} pkg_CAP)
  string(REPLACE "-" "_" pkg_CAP "${pkg_CAP}")

  message(STATUS
    "[${pkg_name}] FetchContent_Declare wrapper"
  )
  # because the signature is different between versions,
  # we build the cmake call beforehand
  set(_fetch_content_cmd "FetchContent_Declare(${pkg_name}")
  
  if(fp_NO_FETCH)
    message(VERBOSE
      " :: \"${pkg_name}\" is specified not fetchable, will rely on `find_package` for population"
    )
    string(APPEND _fetch_content_cmd " DOWNLOAD_COMMAND \":\"") 
    set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_NOFETCH TRUE)
  else()
    message(VERBOSE 
      " :: \"${pkg_name}\" is fetchable, will fall-back to git clone [${fp_GIT_REPO}] if other population methods fail"
    )
    string(APPEND _fetch_content_cmd " GIT_REPOSITORY ${fp_GIT_REPO} GIT_TAG ${fp_GIT_TAG}")
  endif()

  if(fp_PATCH)
    message(VERBOSE 
      " :: \"${pkg_name}\" has patch [${fp_PATCH}], will apply after update stage"
    )
    #NOTE: the ` || true` at the end of the command is a due to a bug where the apply fails if
    #      it has already been applied, so the command fails and cmake fails.
    #      https://gitlab.kitware.com/cmake/cmake/-/issues/21146
    string(APPEND _fetch_content_cmd " PATCH_COMMAND git apply --ignore-space-change --ignore-whitespace ${fp_PATCH} || true")
  endif()

  # bifurcation on cmake version
  if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    # versions >= 3.24 will do an implicit `find_package`, so pass on
    # requirements to declaration
    string(APPEND _fetch_content_cmd " FIND_PACKAGE_ARGS COMPONENTS ${fp_COMPONENTS}")
  endif()

  # close the command
  string(APPEND _fetch_content_cmd ")")

  message(DEBUG " :: FC cmd: ${_fetch_content_cmd}")
  # execute the built `FetchContent_Declare(...)`
  cmake_language(EVAL CODE "${_fetch_content_cmd}")
 
  # return some info
  list(APPEND ${fp_NAMESPACE}_DECLARED_EXTERNAL_CONTENT ${pkg_name})
  set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_COMPONETS ${fp_COMPONENTS} )
  set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_ENABLEOPTS ${fp_ENABLE_OPTS} )
  set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_DISABLEOPTS ${fp_DISABLE_OPTS} )
  set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_PRIORS ${fp_PRIORS} )

  if(fp_EXPECTED_TARGETS)
      set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_TARGETS ${fp_EXPECTED_TARGETS} )
  else()
      set(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_TARGETS "${pkg_name}::${pkg_name}" )
  endif()

endmacro()

#
# :: Overview
#   This macro wraps a call to `FetchContent_MakeAvailable`.
#   If we are using a `cmake` prior to 3.24, explicitly do a 
#   `find_package()` with the declared options 
# :: Arguments 
#   single_value:
#   - NAMESPACE - the namespace used in the corrisponding `singularity_content_declare` call
#
macro(singularityeos_content_populate)
  set(options)
  set(one_value_args
    NAMESPACE
  )
  set(multi_value_args
  )

  cmake_parse_arguments(fp "${options}" "${one_value_args}" "${multi_value_args}" "${ARGN}")

  message(STATUS 
    "[${fp_NAMESPACE}] Populating declared content"
  )
  # fill lists to populate
  # if cmake@3.24+, these are just the lists prepared in singularity_content_declare
  # otherwise, manually check `find_package` and remove content if found
  foreach(pkg_name ${${fp_NAMESPACE}_DECLARED_EXTERNAL_CONTENT})
    string(TOUPPER ${pkg_name} pkg_CAP)
    string(REPLACE "-" "_" pkg_CAP "${pkg_CAP}")
    # bifurcation on cmake version
    if (NOT CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
      find_package(${pkg_name}
        COMPONENTS ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_COMPONETS}
        QUIET
      )
      if(${pkg_name}_FOUND)
        message(VERBOSE
          "${pkg_name} located with `find_package`"
          "${pkg_name}_DIR: ${${pkg_name}_DIR}"
        )
      else()
        # if no fetching and not found, produce an error
        # conditionally include a custom error msg
        if(${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_NOFETCH)
          message(FATAL_ERROR
            "${pkg_name} is requested, but it was not located and is not declared as fetchable.\n"
            "if ${pkg_name} is installed, set \"-D${pkg_name}_ROOT=<install-dir>\""
          )
        endif() # NOFETCH
        message(VERBOSE
          "${pkg_name} NOT located with `find_package`, appending to populate list."
          " Will attempt to clone repository when content is populated."
        )
        list(APPEND _fetchList ${pkg_name})
        list(APPEND _fetchOptsOn ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_ENABLEOPTS})
        list(APPEND _fetchOptsOff ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_DISABLEOPTS})
        list(APPEND _fetchTars ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_TARGETS})

      endif() # FOUND
    else()
      list(APPEND _fetchList ${pkg_name})
      list(APPEND _fetchOptsOn ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_ENABLEOPTS})
      list(APPEND _fetchOptsOff ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_DISABLEOPTS})
      list(APPEND _fetchTars ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_TARGETS})
    endif() #CMAKE_VERSION
    # collect all targets, reguardless of populated
    list(APPEND _expectedTars ${${fp_NAMESPACE}_DECLARED_EXTERNAL_${pkg_CAP}_TARGETS})
  endforeach()

  # for content to be populated, set some cache options given in singularity_content_declare
  foreach(ext_opt ${_fetchOptsOn})
    message(DEBUG "setting \"${ext_opt}\"=ON")
    set(${ext_opt} ON CACHE INTERNAL "")
  endforeach()
  foreach(ext_opt ${_fetchOptsOff})
    message(DEBUG "setting \"${ext_opt}\"=OFF")
    set(${ext_opt} ON CACHE INTERNAL "")
  endforeach()


  message(VERBOSE "\n"
      " :: Populating dependency targets ${_fetchTars}\n"
      " :: Calling `FetchContent_MakeAvailable` with ${_fetchList}\n"
  )
  message(STATUS
      "FetchContent_MakeAvailable prepared, "
      "this may take a few moments if a download is required...\n"
  )
   # populate
  FetchContent_MakeAvailable(${_fetchList})

  # check that declared targets exist
  foreach(expected_target ${_expectedTars})
    if(NOT TARGET ${expected_target})
      message(FATAL_ERROR
          "target \"${expected_target}\" was expected, but does not exist after population!"
          ) # NOT TARGET
    endif()
  endforeach()

  # return target list
  set(${fp_NAMESPACE}_POPULATED_TARGETS ${_expectedTars} )

  # try to keep scope clean
  unset(${fp_NAMESPACE}_DECLARED_EXTERNAL_CONTENT)
endmacro()


# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
