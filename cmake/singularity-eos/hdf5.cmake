macro(singularity_enable_hdf5 target)
  # TODO not sure how useful this snippet is. I was fine with just leaving it
  # be, but then INCLUDE_DIRS started to go wild and screw up CMake on the CI,
  # so it may be worth blanking.
  if(NOT HDF5_ROOT)
    # changed from HDF5_INCLUDE_DIRS to avoid collision with
    # `find_package(HDF5)`
    find_path(
      _TMP_HDF5_DIR hdf5.h
      HINTS ENV HDF5_ROOT
      PATH_SUFFIXES include)
    # TODO dont set if _TMP_HDF5_DIR-NOTFOUND is set
    get_filename_component(HDF5_ROOT ${_TMP_HDF5_DIR} DIRECTORY)
    unset(_TMP_HDF5_DIR-NOTFOUND)
  endif()

  # mauneyc 20230525 SOOOOO..... HDF5 detection has a lot of unremarked on
  # dependencies and build checkpoints that turn this into an inscrutable
  # failure at times. here's what I've found, though I can't really make sense
  # of when/how/where this applies 0.) TODO one or all of these may be due to
  # modifying `HDF5_` output vars before doing `find_package`, which seems to
  # produce pathological behavior. Test this soon. 1.) if 'C' component is
  # requested, 'enable_language(C)' needs to be done beforehand. otherwise,
  # FindHDF5.cmake _may_ fail to detect anything. 1a.) I say _may_ because this
  # behavior is sporadic across CMake/HDF5 versions. 2.) the 'C' component _may_
  # be required. On some installs, 'HL' is sufficent; on others, both are
  # needed. 3.) MPI behavior: rather than just *tell* us what MPI libraries it
  # needs, we have to just guess. Using both MPI_C and MPI_CXX as a
  # 'kitchen-sink' approach is what works (currently) although this may be an
  # issue in other platforms. Unfortunately, this won't give a configuration
  # error, as the target check isn't performed until export. Will have to think
  # about trying to catch these issues.
  if(NOT HDF5_LIBRARIES)
    find_package(
      HDF5
      COMPONENTS C HL
      REQUIRED)
  endif()

  # I'm bailing out here if these arn't filled in. I don't know if this is
  # correct for every downstream use, but right now we need to enforce some kind
  # of uniformity, otherwise we will need 100s of toolchain files JUST for hdf5

  if(NOT HDF5_LIBRARIES OR NOT HDF5_INCLUDE_DIRS)
    message(SEND_ERROR "Could not locate the required HDF5 libraries!")
    message(SEND_ERROR "-- HDF5_LIBRARIES=\"${HDF5_LIBRARIES}\"")
    message(SEND_ERROR "-- HDF5_INCLUDE_DIRS=\"${HDF5_INCLUDE_DIRS}\"")
    message(FATAL_ERROR "Cannot continue without HDF5 libraries!")
  endif()

  # NOTE: `HDF5_LIBRARIES` should have "all requested bindings" according to
  # docs, but some distributions completely seperate the HL libs to
  # `HDF5_HL_LIBRARIES`. or sometimes they don't. So, check for HL, if it's not
  # set, then check if they are placed in HDF5_LIBRARIES (simple substring
  # search). If, finally, they're not found then bail.
  if(NOT HDF5_HL_LIBRARIES)
    string(FIND ${HDF5_LIBRARIES} "hdf5_hl" _hl_in_libs)
    # if HL libs are not found in HDF5_LIBRARIES, then bail
    if(_hl_in_libs LESS 0)
      message(
        SEND_ERROR
          "Found HDF5 libraries, but could not locate required HDF5_HL libraries!"
      )
      message(SEND_ERROR "-- HDF5_LIBRARIES=\"${HDF5_LIBRARIES}\"")
      message(SEND_ERROR "-- HDF5_HL_LIBRARIES=\"${HDF5_HL_LIBRARIES}\"")
      message(FATAL_ERROR "Cannot continue without HDF5 and HDF5_HL libraries!")
    endif()
  endif()

  target_include_directories(${target} SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${target} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

  if(HDF5_IS_PARALLEL)
    # find_package(MPI COMPONENTS C CXX REQUIRED)
    # target_link_libraries(${target} PUBLIC MPI::MPI_C MPI::MPI_CXX)
    enable_language(C)
    find_package(
      MPI
      COMPONENTS C CXX
      REQUIRED)
    target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
    set(SINGLUARITY_USE_SPINER_WITH_PARALLEL_HDF5 TRUE)
  endif()

  target_compile_definitions(${target} PUBLIC SINGULARITY_USE_HDF5)

endmacro()
