Integrating `singularity-eos` into a CMake project
==================================================

Installing singularity-eos into a prefix, either manually or via Spack gives
you an installation that can  be used in other projects via CMake. Assuming
your installation is visible to CMake, either by installing in a known prefix
or setting ``CMAKE_PREFIX_PATH`` manually, a CMake project can integrate
Singularity-EOS via ``find_package(singularity-eos)`` and using the provided
targets to link to it.

.. code::

   find_package(singularity-eos)
   ...
   target_link_libraries(yourTarget PRIVATE singularity-eos::singularity-eos)
    
The ``singularity-eos`` config module provides the following targets:

``singularity-eos::singularity-eos_Interface``:
  The C++ header-only library, adding the necessary include directories and
  other dependencies such as mpark-variant and ports-of-call.

``singularity-eos::singularity-eos_Library``:
  The static or shared library installed when building with
  ``SINGULARITY_BUILD_CLOSURE=on``. If ``SINGULARITY_USE_FORTRAN=on``, this also includes the
  Fortran bindings and adds include path for its Fortran module.

``singularity-eos::singularity-eos``
  Convenience target that contains either
  ``singularity-eos::singularity-eos_Interface``,
  ``singularity-eos::singularity-eos_Library``, or both, depending on the
  COMPONENTS selection during ``find_package``. By default, if no COMPONENTS
  are specified both are included.

Example: Integrating singularity-eos into a Fortran code
--------------------------------------------------------

Fortran projects do not require the C++ header-only library and its
dependencies, but only the compiled Fortran bindings provided by
Singularity-EOS. To avoid unnecessary dependency checks by CMake, a Fortran
project would integrate Singularity-EOS as follows:

.. code::

   find_package(singularity-eos COMPONENTS Library)
   ...
   target_link_libraries(yourTarget PRIVATE singularity-eos::singularity-eos)
