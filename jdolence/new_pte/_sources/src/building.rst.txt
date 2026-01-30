.. _building:

Building Singularity-EOS
=========================

There are two main ways to build singularity-eos:
1. Through the spack interface
2. Through the cmake system

Installing Via Spack
--------------------

.. warning::
  The spack build is currently experimental. 
  Please report problems you havee as github issues.

Although the spackage has not yet made it to the main `Spack`_
repositories, we provide a spackage for ``singularity-eos`` witin the
the singularity-eos source repository. If you have spack installed,
simply call

.. _Spack: https://spack.io/

.. code-block:: bash

  git clone --recursive git@github.com:lanl/singularity-eos.git
  spack repo add singularity-eos/spack-repo
  spack install singularity-eos

to install ``singularity-eos`` into your spack instance. The spackage
supports a number of relevant variants:

+-----------------------------+-----------------+-----------------------------+
| Variant Name [default]      | Allowed Values  | Description                 |
+=============================+=================+=============================+
| build_extra [none]          | none, sesame,   | Build sesame2spiner         |
|                             | stellarcollapse | or stellarcollapse2spiner   |
+-----------------------------+-----------------+-----------------------------+
| build_type [RelWithDebInfo] | Debug, Release, | Equivalent to               |
|                             | RelWitHDebInfo, | -DCMAKE_BUILD_TYPE          |
|                             | MinSizeRel      | in cmake build              |
+-----------------------------+-----------------+-----------------------------+
| cuda [off]                  | on, off         | Build with cuda             |
+-----------------------------+-----------------+-----------------------------+
| cuda_arch [none]            | see kokkos spec | The target GPU architecture |
+-----------------------------+-----------------+-----------------------------+
| doc [off]                   | on, off         | Build sphinx docs           |
+-----------------------------+-----------------+-----------------------------+
| format [off]                | on, off         | Support for clang-format    |
+-----------------------------+-----------------+-----------------------------+
| fortran [on]                | on, off         | Provide fortran bindings    |
+-----------------------------+-----------------+-----------------------------+
| ipo [off]                   | on, off         | CMake interprocedural       |
|                             |                 | optimization                |
+-----------------------------+-----------------+-----------------------------+
| kokkos [off]                | on, off         | Enable Kokkos backend       |
|                             |                 | Required for cuda support   |
+-----------------------------+-----------------+-----------------------------+
| kokkos-kernels [off]        | on, off         | Use kokkos-kernels for      |
|                             |                 | linear algebra suport,      |
|                             |                 | which is needed with        |
|                             |                 | mixed-cell closures on GPU  |
+-----------------------------+-----------------+-----------------------------+
| mpi [off]                   | on, off         | Build with parallel HDF5    |
|                             |                 | otherwise build with serial |
+-----------------------------+-----------------+-----------------------------+
| openmp [off]                | on, off         | Build Kokkos openmp backend |
+-----------------------------+-----------------+-----------------------------+
| tests [off]                 | on, off         | Build tests                 |
+-----------------------------+-----------------+-----------------------------+

Building Via Cmake
--------------------

The `cmake`_ build offers a few more options. For example, it
supports building without ``hdf5`` and bulding "in-tree" by adding the
``singularity-eos`` directory via cmake's ``add_subdirectory``.

.. _cmake: https://cmake.org/

For example, if ``singularity-eos`` is a submodule in the ``external``
directory of your project, you might call a line like this in your
``CMakeLists.txt``.

.. code-block:: cmake

  add_subdirectory(external/singularity-eos singularity-eos)

At it's simplest, the cmake build process looks like this:

.. code-block:: bash

  git clone --recursive git@github.com:lanl/singularity-eos.git
  cd singularity-eos
  mkdir bin
  cd bin
  cmake ..
  make install

You can set options on the cmake line via, e.g.,

.. code-block:: bash

  cmake -DSINGULARITY_USE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug ..

at the cmake configure line. The following cmake options (in addition
to the standard ones) are supported:

+------------------------------------------+------------+---------+-----------------------------------------------+
| Option [default]                         | Value Type | Default | Description                                   |
+==========================================+============+=========+===============================================+
| SINGULARITY_USE_HDF5                     | boolean    | OFF     | Use HDF5                                      |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_USE_FORTRAN                  | boolean    | ON      | Enable fortran bindings                       |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_USE_KOKKOS                   | boolean    | OFF     | Use Kokkos backend. Required for GPU support  |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_USE_EOSPAC                   | boolean    | OFF     | Enable eospac. Required for sesame2spiner and |
|                                          |            |         | for the eospac equation of state class.       |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_USE_CUDA                     | boolean    | OFF     | Enable cuda                                   |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_USE_KOKKOSKERNELS            | boolean    | OFF     | Use kokkos kernels for linear algebra.        |
|                                          |            |         | Linear algebra is needed for mixed cell       |
|                                          |            |         | closures. And kokkos kernels is required for  |
|                                          |            |         | linear algebra on GPU. If kokkos kernels is   |
|                                          |            |         | disabled, Eigen is used.                      |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_BUILD_TESTS                  | boolean    | OFF     | Turn on testing                               |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_BUILD_EXAMPLES               | boolean    | OFF     | Build code in examples directory              |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_BUILD_SESAME2SPINER          | boolean    | OFF     | Build converter from sesame to sp5 tables     |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER | boolean    | OFF     | Build converter from stellar collapse         |
|                                          |            |         | tables to sp5 format.                         |
|                                          |            |         | This is not required to use the               |
|                                          |            |         | stellar collapse reader, but sp5 files are    |
|                                          |            |         | faster to load.                               |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_BETTER_DEBUG_FLAGS           | boolean    | ON      | Makes for more verbose compiler output        |
|                                          |            |         | but can cause problems for in-tree builds.    |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_HIDE_MORE_WARNINGS           | boolean    | OFF     | Makes for less verbose compiler output        |
|                                          |            |         | but can cause problems for in-tree builds.    |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_SUBMODULE_MODE               | boolean    | OFF     | Set other options for in-tree builds          |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_BUILD_CLOSURE                | boolean    | ON      | Build mixed cell closure models               |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_TEST_SESAME                  | boolean    | OFF     | Test the sesame table readers                 |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_TEST_STELLAR_COLLAPSE        | boolean    | OFF     | Test stellar collapse readers                 |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_USE_SINGLE_LOGS              | boolean    | OFF     | Use single-precision logs. Can harm accuracy. |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_FMATH_USE_ORDER_4            | boolean    | OFF     | Use 4th- or 5th-order accurate fast logs.     |
+------------------------------------------+------------+---------+ This is faster but less accurate.             |
| SINGULARITY_FMATH_USE_ORDER_4            | boolean    | OFF     | The default accuracy is 7th-order.            |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_HDF5_INSTALL_DIR             | string     | NONE    | Location of external library.                 |
|                                          |            |         | Not needed, but a hint for cmake.             |
+------------------------------------------+------------+---------+                                               |
| SINGULARITY_MPI_INSTALL_DIR              | string     | NONE    |                                               |
+------------------------------------------+------------+---------+                                               |
| SINGULARITY_KOKKOS_INSTALL_DIR           | string     | NONE    |                                               |
+------------------------------------------+------------+---------+                                               |
| SINGULARITY_KOKKOSKERNERNELS_INSTALL_DIR | string     | NONE    |                                               |
+------------------------------------------+------------+---------+-----------------------------------------------+
| SINGULARITY_KOKKOSKERNELS_SUB_DIR        | string     | NONE    | Set this to build kokkos-kernels "in-tree"    |
|                                          |            |         | by adding it as a subdirectory.               |
+------------------------------------------+------------+---------+-----------------------------------------------+

Dependencies
------------

``singularity-eos`` has a number of dependencies that are handled in a
number of different ways:

* `spiner`_ is a required dependency, included as a submodule
* `hdf5`_ is an optional dependency. It is needed for the table
  readers. If you want it, it must be installed externally and
  findable by ``cmake``. ``MPI`` is an optional dependency of
  ``hdf5``, but otherwise not needed.
* `eospac`_ is an optional dependency. This is needed if you want to
  use sesame tables. If you want it, it must be installed externally
  and findable by ``cmake``
* `kokkos`_ is an optional dependency. It provides GPU support. If it's
  available externally, ``singularity-eos`` will use the available
  version. If not, ``singularity-eos`` will use its own version,
  packaged as a submodule.
* `Eigen`_ is an optional dependency and is used for linear algebra on
  the CPU when doing mixed-cell closures. If it's available
  externally, ``singularity-eos`` will use the available version. If
  not, ``singularity-eos`` will use its own version, packaged as a
  submodule.
* `kokkos-kernels`_ is an optional dependency. This must be available
  externally if desired, but there are a number of ways to expose
  it. One can set ``SINGULARITY_KOKKOSKERNELS_SUB_DIR`` to tell
  ``cmake`` where to ``add_subdirectory`` to make it available. One
  can also simply let ``cmake`` find a pre-installed version of the
  library.
* A fortran compiler is required if fortran bindings are enabled.

.. _spiner: https://github.com/lanl/spiner

.. _hdf5: https://www.hdfgroup.org/solutions/hdf5/

.. _eospac: https://laws.lanl.gov/projects/data/eos/eospacReleases.php

.. _kokkos: https://github.com/kokkos/kokkos

.. _Eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page

.. _kokkos-kernels: https://github.com/kokkos/kokkos-kernels/

If you use spack, but would like to build ``singularity-eos`` from
source, you can install dependencies via, e.g.,

.. code-block:: bash

  git clone --recursive git@github.com:lanl/singularity-eos.git
  spack repo add singularity-eos/spack-repo
  spack install --only dependencies singularity-eos+cuda cuda_arch=70

which will install all the dependencies for the variant of ``singularity-eos`` you've chosen.

Spack can also be used to generate a cmake configuration file based on the 
package variants, so that your development environment and build configuration
are consistent

.. code-block:: bash

   spack install singularity-eos
   spack load singularity-eos
   cd <to/build/dir>
   cmake -C $SINGULARITY_SPACK_CMAKE_CONFIG <path/to/source/dir>


