Overview
========

The ``singularity-eos`` build system is designed with two goals in mind

1. Portability to a wide range of host codes, system layouts, and
   underlying hardware
2. Ease of code development, and flexibility for developers

These considerations continue to guide development of the tools and
workflows in working with ``singularity-eos``.

Basics
------

The build of ``singularity-eos`` can take two forms:

1. Submodule mode
2. Standalone mode

These will be described in more detail below, but in brief **submodule
mode** is intended for downstream codes that build ``singularity-eos``
source code directly in the build (sometimes referred to as "in-tree"),
while **standalone mode** will build ``singularity-eos`` as an independent
library that can be installed onto the system.

The most important distinction between the modes is how dependencies are
handled. *submodule mode* will use *internal* source clones of key
dependencies (located in ``utils\``), effectively building these
dependencies as part of the overall ``singularity-eos`` build procedure.
It should be noted, however, that there are optional dependencies that
are not provided internally and must be separately available.

In *standalone mode*, **all** dependencies must be available in the
environment, and be discoverable to CMake. While not required, it is
encouraged to use the dependency management tool ``spack`` to help
facilitate constructing a build environment, as well as deploying
``singularity-eos``. Example uses of ``spack`` for these purposes are
provided below.

A CMake configuration option is provided that allows developers to
select a specific mode (``SINGULARITY_FORCE_SUBMODULE_MODE``), however
this is intended for internal development only. The intended workflow is
to let ``singularity-eos`` decide that appropriate mode, which it
decides based on inspecting the project directory that the source
resides in.

Dependencies
------------

``singularity-eos`` has a number of required and optional depdencies. 

====================================== =============================== ===========================================
  Package Name                          Distribution                    Comment 
====================================== =============================== ===========================================
 `spiner`_                              submodule [*]_ / external [*]_   Required
 `ports-of-call`_                       submodule / external             Required
 `mpark_variant`_                       submodule / external             Required
 `hdf5`_                                external only                    Optional; used for table I/O
 `eospac`_                              external only                    Optional; used for sesame tables.
 `kokkos`_                              submodule / external             Optional; enables GPU offloading.
 `Eigen`_                               submodule / external             Optional; used for linear algebra on the CPU when doing mixed-cell closures.
 `kokkos-kernels`_                      submodule / external             Optional; used for linear algebra on the GPU when doing mixed-cell closures.
 `pybind11`_                            external / fetchable [*]_        Optional; 
====================================== =============================== ===========================================

.. [*] availible as a git submodule for in-tree builds
.. [*] located outside the build tree and discoverable by CMake
.. [*] CMake can download and configure this source in-tree

.. _spiner: https://github.com/lanl/spiner

.. _ports-of-call: https://github.com/lanl/spiner

.. _mpark_variant: https://github.com/mpark/variant

.. _hdf5: https://www.hdfgroup.org/solutions/hdf5/

.. _eospac: https://laws.lanl.gov/projects/data/eos/eospacReleases.php

.. _kokkos: https://github.com/kokkos/kokkos

.. _Eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page

.. _kokkos-kernels: https://github.com/kokkos/kokkos-kernels/

.. _pybind11: https://github.com/pybind/pybind11

A FORTRAN compiler is required if fortran bindings are enabled.


Options for configuring the build
---------------------------------

Most configuration options are the same between the two builds.
*standalone* / *submodule* specific options are touched on in the
sections detailing those build modes.

The main CMake options to configure building are in the following table:

====================================== ======= ===========================================
  Option                               Default  Comment 
====================================== ======= ===========================================
 ``SINGULARITY_USE_SPINER``              ON       Enables EOS objects that use ``spiner``.
 ``SINGULARITY_USE_FORTRAN``             ON       Enable Fortran API for equation of state.
 ``SINGULARITY_USE_KOKKOS``              OFF      Uses Kokkos as the portability backend. Currently only Kokkos is supported for GPUs.
 ``SINGULARITY_USE_EOSPAC``              OFF      Link against EOSPAC. Needed for sesame2spiner and some tests.
 ``SINGULARITY_BUILD_TESTS``             OFF      Build test infrastructure.
 ``SINGULARITY_BUILD_PYTHON``            OFF      Build Python bindings.
 ``SINGULARITY_INVERT_AT_SETUP``         OFF      For tests, pre-invert eospac tables.
 ``SINGULARITY_BETTER_DEBUG_FLAGS``      ON       Enables nicer GPU debug flags. May interfere with in-tree builds as a submodule.
 ``SINGULARITY_HIDE_MORE_WARNINGS``      OFF      Makes warnings less verbose. May interfere with in-tree builds as a submodule.
 ``SINGULARITY_FORCE_SUBMODULE_MODE``    OFF      Force build in _submodule_ mode.
 ``SINGULARITY_USE_SINGLE_LOGS``         OFF      Use single precision logarithms (may degrade accuracy).
 ``SINGULARITY_USE_TRUE_LOG_GRIDDING``   OFF      Use grids that conform to logarithmic spacing.
====================================== ======= ===========================================

More options are available to modify only if certain other options or
variables satisfy certain conditions (*dependent options*). *Dependent
options* can only be accessed if their precondition is satisfied.

If the precondition is satisfied, they take on a default value, although
they can be changed. If the precondition is *not* satisfied, then their
value is fixed and cannot be changed. For instance,

.. code:: bash

   # in <top-level>/build
   cmake .. -DSINGULARITY_USE_KOKKOS=OFF -DSINGULARITY_USE_CUDA=ON

will have no effect (i.e. ``SINGULARITY_USE_CUDA`` will be set to
``OFF``), because the precondition of ``SINGULARITY_USE_CUDA`` is for
``SINGULARITY_USE_KOKKOS=ON``.

Generally, *dependent options* should only be used for specific
use-cases where the defaults are not applicable. For most scenarios, the
preconditions and defaults are logically constructed and the most
natural in practice (``SINGULARITY_TEST_*`` are only available if
``SINGLARITY_BUILD_TESTS`` is enabled, for instance).

These options are listed in the following table, along with their
preconditions:

============================================== ================================================================================= ===========================================
  Option                                       Precondition                                                                       Comment 
============================================== ================================================================================= ===========================================
 ``SINGULARITY_USE_SPINER_WITH_HDF5``           ``SINGULARITY_USE_SPINER=ON``                                                     Requests that ``spiner`` be configured for ``HDF5`` support.
 ``SINGULARITY_USE_CUDA``                       ``SINGULARITY_USE_KOKKOS=ON``                                                     Target nvidia GPUs for ``Kokkos`` offloading.
 ``SINGULARITY_USE_KOKKOSKERNELS``              ``SINGULARITY_USE_KOKKOS=ON``                                                     Use Kokkos Kernels for linear algebra. Needed for mixed cell closure models on GPU.
 ``SINGULARITY_BUILD_CLOSURE``                  ``SINGULARITY_USE_KOKKOS=ON`` ``SINGULARITY_USE_KOKKOSKERNELS=ON``                Mixed cell closure.
 ``SINGULARITY_BUILD_SESAME2SPINER``            ``SINGULARITY_USE_SPINER=ON`` ``SINGULARITY_USE_SPINER_WITH_HDF5=ON``             Builds the conversion tool sesame2spiner which makes files readable by SpinerEOS.
 ``SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER``   ``SINGULARITY_USE_SPINER=ON`` ``SINGULARITY_USE_SPINER_WITH_HDF5=ON``             Builds the conversion tool stellarcollapse2spiner which optionally makes stellar collapse files faster to read.
 ``SINGULARITY_TEST_SESAME``                    ``SINGULARITY_BUILD_TESTS=ON`` ``SINGULARITY_BUILD_SESAME2SPINER=ON``             Test the Sesame table readers.
 ``SINGULARITY_TEST_STELLAR_COLLAPSE``          ``SINGULARITY_BUILD_TESTS=ON`` ``SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER=ON``     Test the Stellar Collapse table readers.
 ``SINGULARITY_TEST_PYTHON``                    ``SINGULARITY_BUILD_TESTS=ON`` ``SINGULARITY_BUILD_PYTHON=ON``                    Test the Python bindings.
============================================== ================================================================================= ===========================================

CMake presets
-------------

To further aid the developer, ``singularity-eos`` is distributed with
**Presets**, a list of common build options with naturally named labels
that when used can reduce the need to input and remember the many
options ``singularity-eos`` uses. For a general overview of CMake
presets, see the `cmake documentation on
presets <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`__

.. warning::
  CMake presets are only available if ``singularity-eos`` is the 
  top-level project.

Predefined presets
~~~~~~~~~~~~~~~~~~

Predefined presets are described with a ``json`` schema in the file
``CMakePresets.json``. As an example:

.. code:: bash

   # in <top-level>/build
   $> cmake .. --preset="basic_with_testing"
   Preset CMake variables:

     CMAKE_EXPORT_COMPILE_COMMANDS="ON"
     SINGULARITY_BUILD_TESTS="ON"
     SINGULARITY_USE_EOSPAC="ON"
     SINGULARITY_USE_SPINER="ON"

   # ...

As you can see, CMake reports the configuration variables that the
preset has used, and their values. A list of presets can be easily
examined with:

.. code:: bash

   # in <top-level>/build
   $> cmake .. --list-presets
   Available configure presets:

     "basic"
     "basic_with_testing"
     "kokkos_nogpu"
     "kokkos_nogpu_with_testing"
     "kokkos_gpu"
     "kokkos_gpu_with_testing"

When using presets, additional options may be readily appended to
augment the required build. For example, suppose that the ``basic``
preset is mostly sufficient, but you would like to enable building the
closure models:

.. code:: bash

   # in <top-level>/build
   $> cmake .. --preset="basic_with_testing" -DSINGULARITY_BUILD_CLOSURE=ON
   # ...

User defined presets
~~~~~~~~~~~~~~~~~~~~

The CMake preset functionality includes the ability of developers to
define local presets in ``CMakeUserPresets.json``. ``singularity-eos``
explicitly does not track this file in Git, so developers can construct
their own presets. All presets in the predefined ``CMakePresets.json``
are automatically included by CMake, so developers can build off of
those if needed.

For instance, suppose you have a local checkout of the ``kokkos`` and
``kokkos-kernels`` codes that you're using to debug a GPU build, and you
have these installed in ``~/scratch/``. Your ``CMakeUserPresets.json``
could look like:

.. code:: json

   {
     "version": 1,
     "cmakeMinimumRequired": {
       "major": 3,
       "minor": 19
     },
     "configurePresets": [
       {
         "name": "my_local_build",
         "description": "submodule build using a local scratch install of kokkos",
         "inherits": [
           "kokkos_gpu_with_testing"
         ],
         "cacheVariables": {
           "Kokkos_DIR": "$env{HOME}/scratch/kokkos/lib/cmake/Kokkos",
           "KokkosKernels_DIR": "$env{HOME}/scratch/kokkoskernels/lib/cmake/KokkosKernels",
           "SINGULARITY_BUILD_PYTHON": "ON",
           "SINGULARITY_TEST_PYTHON": "OFF"
         }
       }
     ]
   }

This inherits the predefined ``kokkos_gpu_with_testing`` preset, sets
the ``Kokkos*_DIR`` cache variables to point ``find_package()`` to use
these directories, and finally enables building the python bindings
without including the python tests.

Building in *submodule mode*
----------------------------

For *submodule mode* to activate, a clone of the ``singularity-eos``
source should be placed below the top-level of a host project

.. code:: bash

   # An example directory layout when using singularity-eos in submodule mode
   my_project
   |_CMakeLists.txt
   |_README.md 
   |_src 
   |_include
   |_tpl/singularity-eos

``singularity-eos`` is then imported using the ``add_subdirectory()``
command in CMake

.. code:: cmake

   # In your CMakeLists.txt
   cmake_minimum_required(VERSION 3.19)
   project(my_project)

   add_executable(my_exec src/main.cc)
   target_include_directories(my_exec include)

   add_subdirectory(tpl/singularity-eos)

   target_link_libraries(my_exec singularity-eos::singularity-eos)

This will expose the ``singularity-eos`` interface and library to your
code, along with the interfaces of the internal dependencies

.. code:: c++

   // in source of my_project 

   #include<singularity-eos/eos/eos.hpp>
   // from the internal ports-of-call submodule
   #include<ports-of-call/portability>

   // ...

   using namespace singularity;

``singularity-eos`` will build (along with internal dependencies) and be
linked directly to your executable.

The git submoudles may change during development, either by changing the
pinned hash, addition or removal of submodules. If you have errors that
appear to be the result of incompatible code, make sure you have updated
your submodules with

.. code:: bash

   git submodule update --init --recursive

Building in *standalone mode*
-----------------------------

For *standalone* mode, all required and optional dependencies are
expected to be discoverable by CMake. This can be done several ways

1. (*preferred*) Use Spack to configure and install all the dependencies
   needed to build.
2. Use a system package manager (``apt-get``, ``yum``, &t) to install
   dependencies.
3. Hand-build to a local filesystem, and configure your shell or CMake
   invocation to be aware of these installs

*standalone* mode is the mode used to install ``singularity-eos`` to a
system as a common library. If, for example, you use Spack to to install
packages, ``singularity-eos`` will be built and installed in
*standalone* mode.

Building with Spack
~~~~~~~~~~~~~~~~~~~

Spack is a package management tool that is designed specifically for HPC
environments, but may be used in any compute environment. It is useful
for gathering, configuring and installing software and it's dependencies
self-consistently, and can use existing software installed on the system
or do a "full" install of all required (even system) packages in a local
directory.

Spack remains under active development, and is subject to rapid change
in interface, design, and functionality. Here we will provide an
overview of how to use Spack to develop and deploy ``singularigy-eos``,
but for more in-depth information, please refer to the `official Spack
documentation <spack.readthedocs.io>`__.

Preperation
^^^^^^^^^^^

First, we need to clone the Spack repository. You can place this
anywhere, but note that by default Spack will download and install
software under this directory. This default behavior can be changed,
please refer to the documentation for information of customizing your
Spack instance.

.. code:: bash

   $> cd ~
   $> git clone https://github.com/spack/spack.git

To start using Spack, we use the provided activation script

.. code:: bash

   # equivalent scripts for tcsh, fish are located here as well 
   $> source ~/spack/share/spack/setup-env.sh

You will always need to *activate* spack for each new shell. You may
find it convienant to invoke this Spack setup in your login script,
though be aware that Spack will prepend paths to your environment which
may cause conflicts with other package tools and software.

The first time a Spack command is invoked, it will need to bootstrap
itself to be able to start *concretizing package specs*. This will
download pre-built packages and create a ``${HOME}/.spack`` directory.
This directory is important and is where your *primary* Spack
configuration data will be located. If at any point this configuration
becomes corrupted or too complicated to easily fix, you may safely
remove this directory to restore the default configuration, or just to
try a new approach. Again, refer to the Spack documentaion for more
information.

Setup compilers
^^^^^^^^^^^^^^^

To use Spack effectively, we need to configure it for the HPC
environment we're using. This can be done manually (by editing
``packages.yaml``, ``compilers.yaml``, and perhaps a few others). This
is ideal if you understand how your software environment is installed on
the HPC system, and you are fluent in the Spack configuration schema.

However, Spack has put in a lot of effort to be able to automatically
discover the available tools and software on any given system. While not
perfect, we can get a fairly robust starting point.

Assume we are on an HPC system that has Envionrmental Modules that
provides compilers, MPI implementations, and sundry other common tools.
To help Spack find these, let's load a specific configuration into the
active shell environment.

.. code:: bash

   $> module load cmake/3.19.2 gcc/11.2.0 openmpi/4.1.1 python/3.10
   $> module list

   Currently Loaded Modules:
     1) cmake/3.19.2   2) gcc/11.2.0   3) openmpi/4.1.1   4) python/3.10-anaconda-2023.03

First, let's find the available compilers. (If this is the first Spack
command you've run, it will need to bootstrap)

.. code:: bash

   $> spack compiler find
   ==> Added 2 new compilers to ${HOME}/.spack/linux/compilers.yaml
       gcc@4.8.5  gcc@11.2.0
   ==> Compilers are defined in the following files:
       ${HOME}/.spack/linux/compilers.yaml

Here, we find the default system compiler (``gcc@4.8.5``), along with
the compiler from the module we loaded. Also notice that the
``${HOME}/.spack`` directory has been modified with some new YAML config
files. These are information on the compilers and how Spack will use
them. You are free to modify these files, but for now let's leave them
as is.

*NB*: You can repeat this procedure for other compilers and packages,
though if you need to use many different combinations of
compiler/software, you will find using Spack *environments* `more
convenient <https://spack.readthedocs.io/en/latest/environments.html>`__.

Setup system-provided packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we will try and find system software (e.g.
``ncurses``,\ ``git``,\ ``zlib``) that we can use instead of needing to
build our own. This will also find the module software we loaded
(``cmake``,\ ``openmpi``,\ ``python``). (This command will take a couple
minutes to complete).

.. code:: bash

   $> spack external find --all --not-buildable
   ==> The following specs have been detected on this system and added to ${HOME}/.spack/packages.yaml
   autoconf@2.69       bzip2@1.0.6     coreutils@8.22  dos2unix@6.0.3    gcc@11.2.0        go@1.16.5            hdf5@1.8.12      libfuse@3.6.1         ncurses@6.4.20221231   openssl@1.1.1t     python@3.10.9   sqlite@3.7.17      texlive@20130530
   automake@1.13.4     bzip2@1.0.8     cpio@2.11       doxygen@1.8.5     gettext@0.19.8.1  go@1.18.4            hdf5@1.10.6      libtool@2.4.2         ninja@1.10.2           perl@5.16.3        rdma-core@22.4  sqlite@3.40.1      which@2.20
   bash@4.2.46         ccache@3.7.7    curl@7.29.0     file@5.11         ghostscript@9.25  go-bootstrap@1.16.5  krb5@1.15.1      lustre@2.12.9         opencv@2.4.5           pkg-config@0.27.1  rsync@3.1.2     subversion@1.7.14  xz@5.2.2
   berkeley-db@5.3.21  cmake@2.8.12.2  curl@7.87.0     findutils@4.5.11  git@2.18.4        go-bootstrap@1.18.4  krb5@1.19.4      m4@1.4.16             openjdk@1.8.0_372-b07  python@2.7.5       ruby@2.0.0      swig@2.0.10        xz@5.2.10
   binutils@2.27.44    cmake@3.17.5    cvs@1.11.23     flex@2.5.37       git-lfs@2.10.0    gpgme@1.3.2          libfabric@1.7.2  maven@3.0.5           openssh@7.4p1          python@3.4.10      sed@4.2.2       tar@1.26           zip@3.0
   bison@3.0.4         cmake@3.19.2    diffutils@3.3   gawk@4.0.2        gmake@3.82        groff@1.22.2         libfuse@2.9.2    ncurses@5.9.20130511  openssl@1.0.2k-fips    python@3.6.8       slurm@23.02.1   texinfo@5.1

   -- no arch / gcc@11.2.0 -----------------------------------------
   openmpi@4.1.1

*Generally* you will want to use as much system-provided software as you
can get away with (in Spack speak, these are called **externals**, though
*external packages* are not limited to system provided ones and can
point to, e.g., a manual install). In the above command, we told Spack
to mark any packages it can find as ``not-buildable``, which means that
Spack will never attempt to build that package and will always use the
external one. This *may* cause issues in resolving packages specs when
the external is not compatible with the requirements of an downstream
package.

As a first pass, we will use ``--not-buildable`` for
``spack external find``, but if you have any issues with concretizing
then start this guide over (remove ``${HOME}/.spack`` and go back to
compilers) and do not use ``--not-buildable`` in the previous command.
You may also manually edit the ``packages.yaml`` file to switch the
``buildable`` flag for the troublesome package, but you will need to be
a least familiar with YAML schema.

First install with spack
^^^^^^^^^^^^^^^^^^^^^^^^

Let's walk through a simple Spack workflow for installing. First, we
want to look at the options available for a package. The Spack team and
package developers have worked over the years to provide an impressive
selection of packages. This example will use ``hypre``, a parallel
library for multigrid methods.

.. code:: bash

   $> spack info hypre
   AutotoolsPackage:   hypre

   Description:
       Hypre is a library of high performance preconditioners that features
       parallel multigrid methods for both structured and unstructured grid
       problems.

   Homepage: https://llnl.gov/casc/hypre

   Preferred version:
       2.28.0     https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz

   Safe versions:
       develop    [git] https://github.com/hypre-space/hypre.git on branch master
       2.28.0     https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz

   # ... more versions listed 

   Variants:
       Name [Default]              When       Allowed values          Description
       ========================    =======    ====================    ==============================================

       amdgpu_target [none]        [+rocm]    none, gfx900,           AMD GPU architecture
                                              gfx1030, gfx90c,
                                              gfx90a, gfx1101,
                                              gfx908, gfx1010,
   # ... lots of amd targets listed 
       build_system [autotools]    --         autotools               Build systems supported by the package
       caliper [off]               --         on, off                 Enable Caliper support
       complex [off]               --         on, off                 Use complex values
       cuda [off]                  --         on, off                 Build with CUDA
       cuda_arch [none]            [+cuda]    none, 62, 80, 90,       CUDA architecture
                                              20, 32, 35, 37, 87,
                                              10, 21, 30, 12, 61,
                                              11, 72, 13, 60, 53,
                                              52, 75, 70, 89, 86,
                                              50
       debug [off]                 --         on, off                 Build debug instead of optimized version
       fortran [on]                --         on, off                 Enables fortran bindings
       gptune [off]                --         on, off                 Add the GPTune hookup code
       int64 [off]                 --         on, off                 Use 64bit integers
       internal-superlu [off]      --         on, off                 Use internal SuperLU routines
       mixedint [off]              --         on, off                 Use 64bit integers while reducing memory use
       mpi [on]                    --         on, off                 Enable MPI support
       openmp [off]                --         on, off                 Enable OpenMP support
       rocm [off]                  --         on, off                 Enable ROCm support
       shared [on]                 --         on, off                 Build shared library (disables static library)
       superlu-dist [off]          --         on, off                 Activates support for SuperLU_Dist library
       sycl [off]                  --         on, off                 Enable SYCL support
       umpire [off]                --         on, off                 Enable Umpire support
       unified-memory [off]        --         on, off                 Use unified memory

   Build Dependencies:
       blas  caliper  cuda  gnuconfig  hip  hsa-rocr-dev  lapack  llvm-amdgpu  mpi  rocprim  rocrand  rocsparse  rocthrust  superlu-dist  umpire

   Link Dependencies:
       blas  caliper  cuda  hip  hsa-rocr-dev  lapack  llvm-amdgpu  mpi  rocprim  rocrand  rocsparse  rocthrust  superlu-dist  umpire

   Run Dependencies:
       None

The ``spack info`` commands gives us three important data-points we
need. First, it tells the versions available. If you do not specify a
version, the *preferred* version is default.

Next and most important are the *variants*. These are used to control
how to build the package, i.e. to build with MPI, to build a fortran
interface, and so on. These will have default values, and in practice
you will only need to provide a small number for any particular system.

Finally, we are given the *dependencies* of the package. The
dependencies listed are for *all* configurations, so some dependencies
may not be necessary for your particular install. (For instance, if you
do not build with ``cuda``, then ``cuda`` will not be necessary to
install)

Let's look at what Spack will do when we want to install. We will start
with the default configuration (that is, all variants are left to
default). The ``spack spec`` command will try to use the active Spack
configuration to determine which packages are needed to install
``hypre``, and will print the dependency tree out.

.. code:: bash

   $> spack spec hypre
   Input spec
   --------------------------------
    -   hypre

   Concretized
   --------------------------------
    -   hypre@2.28.0%gcc@11.2.0~caliper~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~sycl~umpire~unified-memory build_system=autotools arch=linux-rhel7-broadwell
    -       ^openblas@0.3.23%gcc@11.2.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-rhel7-broadwell
   [e]          ^perl@5.16.3%gcc@11.2.0+cpanm+opcode+open+shared+threads build_system=generic patches=0eac10e,3bbd7d6 arch=linux-rhel7-broadwell
   [e]      ^openmpi@4.1.1%gcc@11.2.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+pmi+romio+rsh~singularity+static+vt~wrapper-rpath build_system=autotools fabrics=ofi,psm,psm2 schedulers=slurm arch=linux-rhel7-broadwell

Here, we see the full default Spack *spec*, which as a rough guide is
structured as
``<package>@<version>%<compiler>@<compiler_version>{[+/~]variants} <arch_info>``.
The ``+,~`` variant prefixes are used to turn on/off variants with
binary values, while variants with a set of values are given similar to
keyword values (e.g. ``+cuda cuda_arch=70 ~shared``)

If we wanted to install a different configuration, in this case say we
want ``complex`` and ``openmp`` enabled, but we don't need ``fortran``.

.. code:: bash

   $> spack spec hypre+complex+openmp~fortran
   Input spec
   --------------------------------
    -   hypre+complex~fortran+openmp

   Concretized
   --------------------------------
    -   hypre@2.28.0%gcc@11.2.0~caliper+complex~cuda~debug~fortran~gptune~int64~internal-superlu~mixedint+mpi+openmp~rocm+shared~superlu-dist~sycl~umpire~unified-memory build_system=autotools arch=linux-rhel7-broadwell
    -       ^openblas@0.3.23%gcc@11.2.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-rhel7-broadwell
   [e]          ^perl@5.16.3%gcc@11.2.0+cpanm+opcode+open+shared+threads build_system=generic patches=0eac10e,3bbd7d6 arch=linux-rhel7-broadwell
   [e]      ^openmpi@4.1.1%gcc@11.2.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+pmi+romio+rsh~singularity+static+vt~wrapper-rpath build_system=autotools fabrics=ofi,psm,psm2 schedulers=slurm arch=linux-rhel7-broadwell

Here, you can see the full spec has out supplied variants. In general,
variants can control build options and features, and can change which
dependencies are needed.

Notice also the left-aligned string starting each line for a package.
``-`` indicates that Spack isn't aware that this package is installed
(which is expected). ``[+]`` indicates that the package has been
previously installed. ``[e]`` indicates that the package has been marked
as externally installed.

Finally, we can install it. Because ``perl`` and ``openmpi`` are already
present, Spack will not need to download, build, and install these
packages. This can save lots of time! Note, however, that external
packages are loosely constrained and may not be correctly configured for
the requested package.

*NB*: By default, Spack will try to download the package source from the
repository associated with the package. This behavior can be overrided
with Spack *mirrors* , but that is beyond the scope of this doc.

.. code:: bash

Now, we can use Spack similarly to ``module load``,

.. code:: bash

   $> spack load hypre
   $> spack find --loaded

Other options are available for integrating Spack installed packages
into your environment. For more, head over to
https://spack.readthedocs.io

Installing ``singularity-eos`` using Spack
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

. warning::
  The spack build is currently experimental. 
  Please report problems you havee as github issues.

The spackage is available in the main `Spack`_
repositories, and we provide a spackage for ``singularity-eos`` witin the
the singularity-eos source repository. The distributed spackage may be 
more up-to-date than the one in the main `Spack`_ repository. If you 
have spack installed, simply call

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
| hdf5 [off]                  | on, off         | Enable HDF5 I/O for tables  |
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

Developing ``singularigy-eos`` using Spack
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spack is a powerful tool that can help develop ``singularigy-eos`` for a
variety of platforms and hardware.

1. Install the dependencies ``singularigy-eos`` needs using Spack

.. code:: bash

   $> spack install -u cmake singularity-eos@main%gcc@13+hdf5+eospac+mpi+kokkos+kokkos-kernels+openmp^eospac@6.4.0

This command will initiate an install of ``singularity-eos`` using
Spack, but will stop right before ``singularity-eos`` starts to build
(``-u cmake`` means ``until cmake``). This ensures all the necessary
dependencies are installed and visible to Spack

2. Use Spack to construct an *ad-hoc* shell environment

.. code:: bash

   $> spack build-env singularity-eos@main%gcc@13+hdf5+eospac+mpi+kokkos+kokkos-kernels+openmp^eospac@6.4.0 -- bash

This command will construct a shell environment in ``bash`` that has all
the dependency information populated (e.g. ``PREFIX_PATH``,
``CMAKE_PREFIX_PATH``, ``LD_LIBRARY_PATH``, and so on). Even external
packages from a module system will be correctly loaded. Thus, we can
build for a specific combination of dependencies, compilers, and
portability strategies.

.. code:: bash

   $> salloc -p scaling
   # ...
   $> source ~/spack/share/spack/setup-env.sh
   $> spack build-env singularity-eos@main%gcc@12+hdf5+eospac+mpi+kokkos+kokkos-kernels+openmp^eospac@6.4.0 -- bash
   $> mkdir -p build_gpu_mpi ; cd build_gpu_mpi
   $> cmake .. --preset="kokkos_nogpu_with_testing"
