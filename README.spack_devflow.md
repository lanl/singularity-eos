Spack has become a standard package management solution for software development HPC environments. While a powerful tool, even experienced users have found it's use to be occasionally frustrating. Furthermore, much of the downstream use of Spack tends to emphasize it's use as means to dev-ops. This is, of course, Spack's main context, but this can mislead developers looking for a flexible, extensible programming environment.


This walkthrough is intended to provide users with experience in using Spack for a developer workflow, with the goal of showing 

Assume familiarty with basic Spack usage (for example, defining a package's dependencies and installing those dependencies)

## Setup
From my personal experience, the trivial environment and directory layout is the most difficult to adjust to, while all the iterative development that follows is natural and straightforward. That may be a personal failing, but if you see these first steps as strange, I want to assure you that: 1.) you are in good company and 2.) once we begin doing code work, it all makes perfect sense, and 3.) this approach is not required, and with a few obvious options we can work in almost any conceivable setup.

### NOTE STEP ZERO using a fresh spack clone, I recommend doing this.

Make a directory, everything we do will stay in this directory.

```bash
$ mkdir devel
$ cd devel
$ git clone https://github.com/spack/spack.git
$ export SPACK_USER_CACHE_PATH=$(pwd)/.spack
$ source spack/share/spack/setup-env.sh
$ spack spec zlib
```

The last command will initiate Spack bootstraping. Since it has to be done anyway, let's go ahead and get it over with.

Let's start with creating our top-level environment directory. For this exercise, all code checkouts and local installs will go beneath this directory.

```bash
$ ~ > mkdir -p devel/seos-dev
$ ~ > cd devel/seos-dev
$ seos-dev >
```


Now, we want to create our development environment. We will make an _anonymous_ environment; the best way to think about it is that this directory will be synonymous with our development environment - that is, when we work in anything underneth this directory, it is being done in the environment we will construct.



```bash
$ seos-dev > spack env create -d .
==> Created independent environment in: /home/jack/devel/seos-dev
==> Activate with: spack env activate .
$ seos-dev > spack env activate .
```

First, let's point out what we have after the environment is created

```bash
$ seos-dev > ls -la
drwxr-xr-x@ 4 jack  dudes  128 Mar 28 10:29 .
drwxr-xr-x@ 3 jack  dudes   96 Mar 28 10:29 ..
drwxr-xr-x@ 5 jack  dudes  160 Mar 28 10:29 .spack-env
-rw-r--r--@ 1 jack  dudes  230 Mar 28 10:29 spack.yaml
```

`spack.yaml` will be the configuration file for our environment. Right now it's empty, but that will soon change. `.spack-env` will be a local cache, with the main data being a file-system view that collects our environment similar to a (very) lightweight container.

Adding our packages to the environment
```bash
$ seos-dev > spack add singularity-eos@main
```

Let's assume we want `mpi` and `hdf5` on, since we might as well. While we're at it, let's fix the `mpi` provider to `openmpi`

```bash
$ seos-dev > spack config add "packages:all:variants: +hdf5+mpi"
$ seos-dev > spack config add "packages:mpi:require: openmpi"
```

Looking ahead, we want to have access to `spiner` code, so let's explicitly add this to our packages. (NOTE: This step is not required, but shown as a demo. An example of getting a clone of a non-explicit package will be given after)

```bash
$ seos-dev > spack add spiner@main
```

We've got a basic spack environment ready. Let's take a look at the new `spack.yaml`

```yaml
spack:
#  include:
#  - ./external_packages.yaml
#  - ./external_compilers.yaml
  packages:
    all:
      variants: +hdf5+mpi
    mpi:
      require: openmpi
  specs:
  - singularity-eos@main
  - spiner@main
  view: true
  concretizer:
    unify: true
```

(Here, I've added some `include:` sections for external packages and external compilers. These will be platform specific. See - - below for an example to generate these files)

Now, let's recall: 
- I want to start editing and updating `singularity-eos` code in a new clone and checkout.
- `spiner` and `singularity-eos` have a tight coupling, and I also want a new clone and checkout of that package.
- Perhaps another package?
## Getting and updating `singularity-eos`
First, we will just make some updates to the core package `singularity-eos`. 

The regular pipeline of spack is to go get a clone of repo into a temporary directory and execute that package's build procedure. What we want is to tell spack to get and use a non-temporary clone, preferably to a path we can easily `cd` to and start hacking.

This can be done with the `spack develop` command 
```bash
$ seos-dev > spack develop singularity-eos
==> Cloning source code for singularity-eos@=main
```

We can see that spack has done a `git clone` for us!
```bash
$ seos-dev > ls -la
drwxr-xr-x@ 4 jack  dudes  128 Mar 28 10:29 .
drwxr-xr-x@ 3 jack  dudes   96 Mar 28 10:29 ..
drwxr-xr-x@ 5 jack  dudes  160 Mar 28 10:29 .spack-env
drwxrwxr-x@ 2 jack  dudes 4.0K Mar 28 11:11 singularity-eos
-rw-rw-r--@ 1 jack  dudes  24K Mar 28 11:10 spack.lock
-rw-r--r--@ 1 jack  dudes  230 Mar 28 10:29 spack.yaml
```

From here, we can do regular spack things.
```bash
$ seos-dev > spack concretize -f
==> Warning: Unable to resolve the git commit for singularity-eos. An installation of this binary won't have complete binary provenance.
==> Warning: Unable to resolve the git commit for spiner. An installation of this binary won't have complete binary provenance.
==> Concretized 2 specs:
[+]  66ruf4l  singularity-eos@main+closure~cuda+eospac+fortran+hdf5~ipo+kokkos~kokkos-kernels+mpi~openmp~python+spiner build_extra:=none build_system=cmake build_type=Release dev_path=/vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos generator=make arch=linux-rhel8-haswell %c,cxx,fortran=gcc@13.2.0
[+]  cvvxqem      ^cmake@3.31.6~doc+ncurses+ownlibs~qtgui build_system=generic build_type=Release arch=linux-rhel8-haswell %c,cxx=gcc@13.2.0
[e]  ivwo7eo          ^curl@7.61.1+gssapi+ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs:=shared,static tls:=openssl arch=linux-rhel8-haswell
[+]  b4cxilh          ^ncurses@6.5~symlinks+termlib abi=none build_system=autotools patches:=7a351bc arch=linux-rhel8-haswell %c,cxx=gcc@13.2.0
[e]  p3bp3e3          ^zlib@1.2.11+optimize+pic+shared build_system=makefile arch=linux-rhel8-haswell
[+]  mj4ubyb      ^compiler-wrapper@1.0 build_system=generic arch=linux-rhel8-haswell
[+]  withyuf      ^eigen@3.3.8~ipo~nightly~rocm build_system=cmake build_type=RelWithDebInfo generator=make patches:=55daee8,b8877a8 arch=linux-rhel8-haswell %c,cxx=gcc@13.2.0
[+]  il66js2      ^eospac@6.5.12~offload build_system=generic arch=linux-rhel8-haswell %c,cxx,fortran=gcc@13.2.0
[e]  lbz2bqn      ^gcc@13.2.0~binutils+bootstrap~graphite~mold~nvptx~piclibs~profiled~strip build_system=autotools build_type=RelWithDebInfo languages:='c,c++,fortran' arch=linux-rhel8-haswell
[+]  va64afh      ^gcc-runtime@13.2.0 build_system=generic arch=linux-rhel8-haswell
[e]  i5ej3xm      ^glibc@2.28 build_system=autotools arch=linux-rhel8-haswell
[e]  df7ehnf      ^gmake@4.2.1~guile build_system=generic patches:=ca60bd9,fe5b60d arch=linux-rhel8-haswell
[+]  zaqgsqx      ^hdf5@1.14.6~cxx~fortran+hl~ipo~java~map+mpi+shared~subfiling~szip~threadsafe+tools api=default build_system=cmake build_type=Release generator=make arch=linux-rhel8-haswell %c=gcc@13.2.0
[e]  lpanugo          ^openmpi@5.0.2+atomics~cuda~debug+fortran~gpfs~internal-hwloc~internal-libevent~internal-pmix~ipv6~java~lustre~memchecker~openshmem~romio+rsh~static~two_level_namespace+vt~wrapper-rpath build_system=autotools fabrics:=ucx patches:=1a8ff33 romio-filesystem:=none schedulers:=none arch=linux-rhel8-haswell
[e]  ojffucq          ^pkgconf@1.4.2 build_system=autotools arch=linux-rhel8-haswell
[+]  gjjssb7      ^kokkos@4.6.01~aggressive_vectorization~cmake_lang~compiler_warnings+complex_align~cuda~debug~debug_bounds_check~debug_dualview_modify_check~deprecated_code~examples~hip_relocatable_device_code~hpx~hpx_async_dispatch~hwloc~ipo~memkind~numactl~openmp~openmptarget~pic~rocm+serial~shared~sycl~tests~threads~tuning~wrapper build_system=cmake build_type=Release cxxstd=17 generator=make intel_gpu_arch=none arch=linux-rhel8-haswell %c,cxx=gcc@13.2.0
[+]  e5n25ma      ^mpark-variant@1.4.0~ipo build_system=cmake build_type=Release generator=make patches:=21a4f8d,4e173fe,b3501f7 arch=linux-rhel8-haswell %cxx=gcc@13.2.0
[+]  qhzavje      ^ports-of-call@main~ipo~test build_system=cmake build_type=Release commit=9dad5cb335ee92d7bf4a33ca070ef93e745212d0 generator=make arch=linux-rhel8-haswell %c,cxx=gcc@13.2.0
[+]  kq4og7r  spiner@main+hdf5~ipo+kokkos+mpi~python build_system=cmake build_type=Release dev_path=/vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/spiner generator=make arch=linux-rhel8-haswell %c,cxx=gcc@13.2.0
$ seos-dev > spack install
```


Now, let's make an update.
```bash
$ seos-dev > cd singularity-eos
$ seos-dev > git checkout -b new_eos_feature
$ seos-dev > vi singularity-eos/eos/eos_bash.hpp
```

Here, I'll add my amazing new update
```c++
	15   #ifndef _SINGULARITY_EOS_EOS_EOS_BASE_
    16   #define _SINGULARITY_EOS_EOS_EOS_BASE_
    17   
    18   #include <cstring>
    19   #include <limits>
    20   #include <string>
    21   
    22   #include <ports-of-call/portability.hpp>
    23   #include <ports-of-call/portable_errors.hpp>
    24   #include <singularity-eos/base/constants.hpp>
    25   #include <singularity-eos/base/robust_utils.hpp>
    26   #include <singularity-eos/base/root-finding-1d/root_finding.hpp>
    27   #include <singularity-eos/base/variadic_utils.hpp>
    28   
    29   // MY NEW CODE
    30   // THIS WILL FIX EVERYTHING!
    31   static_assert(false);
```

Instead of manually building, let's use spack to check that the build works. (We will move to a build environment later.)
```bash
$ seos-dev > spack install
[+] /usr (external glibc-2.28-i5ej3xmy3uonjbh2yn3lubq7xmue442j)
[+] /projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0 (external openmpi-5.0.2-lpanugo33x6swa22ws2gvnpt4cyjd6p4)
[+] /usr (external pkgconf-1.4.2-ojffucqu32xqeuph3dnufdahuiqzjijl)
[+] /usr (external zlib-1.2.11-p3bp3e3qjussu7ibvpfolmjd36va7rkq)
[+] /usr (external curl-7.61.1-ivwo7eocoj4lubj2xj7t4luqhffzkng4)
[+] /usr (external gmake-4.2.1-df7ehnfeehbynnqobw5vg4ih5aojixew)
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/compiler-wrapper-1.0-mj4ubybewkayq3jjmsuyqdyskwffbre5
[+] /projects/opt/rhel8/x86_64/gcc/13.2.0 (external gcc-13.2.0-lbz2bqno2yi74phtrhypxnc3quk3ylbu)
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/gcc-runtime-13.2.0-va64afhaljvjalhfhh66dudde5ifaib4
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/eigen-3.3.8-withyufrtkjcrt2nng2nauy7jt6qbluw
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/eospac-6.5.12-il66js2nmooobxtng3sykhzmdhobuxv4
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/kokkos-4.6.01-gjjssb76dkwnmwjvu3lgslvtdknk6soq
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/mpark-variant-1.4.0-e5n25mavzyoqywk4w6g3xwbx5f7kerwh
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/ports-of-call-main-qhzavjeh5tdrzm55u7uw74j4po3rd6ra
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/hdf5-1.14.6-zaqgsqxckmewfgc2wdptlmd2yrgpqp54
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/ncurses-6.5-b4cxilhggxjrp3ly2bwkh42rogdbvg2c
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/spiner-main-kq4og7rxelnbfqef45rzntjw3lds2uen
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/cmake-3.31.6-cvvxqemwg7lkzpokbh3mt4vbkma6ohnm
==> No binary for singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc found: installing from source
==> Installing singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc [19/19]
==> No patches needed for singularity-eos
==> singularity-eos: Executing phase: 'cmake'
==> singularity-eos: Executing phase: 'build'
==> Error: ProcessError: Command exited with status 2:
    '/usr/bin/make' '-j16'

14 errors found in build log:
     43     /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/compiler-wrapper-1.0-mj4ubybewkayq3jjmsuyqdyskwffbre5/libexec/spack/gcc/g++ -DH5_BUILT_AS_DYNAMIC_LIB -DKOKKOS_DEPENDENCE -DPORTABILITY_STRATEGY_KOKKOS -DSINGULARITY_BUILD_CLOSURE -DSINGULARITY_USE_EOSPAC -DSINGULARITY_USE_HDF5 -DSINGULARIT
            Y_USE_SPINER -DSINGULARITY_USE_SPINER_WITH_HDF5 -DSINGULARITY_VERSION=\"1.9.2\" -DSINGULARITY_VERSION_MAJOR=1 -DSINGULARITY_VERSION_MINOR=9 -DSINGULARITY_VERSION_PATCH=2 -DSPINER_USE_HDF -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE_SOURCE -D_POSIX_C_SOURCE=200809L -I/vast/home/mauneyc/develop/xcap-de
            v/dev-env-workflow/singularity-eos/eospac-wrapper -I/vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos -I/tmp/mauneyc/spack-stage/spack-stage-singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc/spack-build-66ruf4l/generated -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/l
            inux-haswell/hdf5-1.14.6-zaqgsqxckmewfgc2wdptlmd2yrgpqp54/include -isystem /vast/projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0/include -isystem /projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/eospac-6.5.12-il66js2n
            mooobxtng3sykhzmdhobuxv4/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/kokkos-4.6.01-gjjssb76dkwnmwjvu3lgslvtdknk6soq/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/mpark-variant-1.4.0-e5n25mavzyoqywk4w6g3xwbx5f7kerwh/include -isyste
            m /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/ports-of-call-main-qhzavjeh5tdrzm55u7uw74j4po3rd6ra/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/spiner-main-kq4og7rxelnbfqef45rzntjw3lds2uen/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack
            /opt/spack/linux-haswell/eigen-3.3.8-withyufrtkjcrt2nng2nauy7jt6qbluw/include/eigen3 -O3 -DNDEBUG -pthread -march=core-avx2 -mtune=core-avx2 -Wno-psabi -Wno-class-memaccess -MD -MT CMakeFiles/singularity-eos_Library.dir/singularity-eos/eos/get_sg_eos_rho_p.cpp.o -MF CMakeFiles/singularity-eos_Library.dir/
            singularity-eos/eos/get_sg_eos_rho_p.cpp.o.d -o CMakeFiles/singularity-eos_Library.dir/singularity-eos/eos/get_sg_eos_rho_p.cpp.o -c /vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos/singularity-eos/eos/get_sg_eos_rho_p.cpp
     44     [ 80%] Building CXX object CMakeFiles/singularity-eos_Library.dir/singularity-eos/eos/get_sg_eos_rho_e.cpp.o
     45     /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/compiler-wrapper-1.0-mj4ubybewkayq3jjmsuyqdyskwffbre5/libexec/spack/gcc/g++ -DH5_BUILT_AS_DYNAMIC_LIB -DKOKKOS_DEPENDENCE -DPORTABILITY_STRATEGY_KOKKOS -DSINGULARITY_BUILD_CLOSURE -DSINGULARITY_USE_EOSPAC -DSINGULARITY_USE_HDF5 -DSINGULARIT
            Y_USE_SPINER -DSINGULARITY_USE_SPINER_WITH_HDF5 -DSINGULARITY_VERSION=\"1.9.2\" -DSINGULARITY_VERSION_MAJOR=1 -DSINGULARITY_VERSION_MINOR=9 -DSINGULARITY_VERSION_PATCH=2 -DSPINER_USE_HDF -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE_SOURCE -D_POSIX_C_SOURCE=200809L -I/vast/home/mauneyc/develop/xcap-de
            v/dev-env-workflow/singularity-eos/eospac-wrapper -I/vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos -I/tmp/mauneyc/spack-stage/spack-stage-singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc/spack-build-66ruf4l/generated -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/l
            inux-haswell/hdf5-1.14.6-zaqgsqxckmewfgc2wdptlmd2yrgpqp54/include -isystem /vast/projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0/include -isystem /projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/eospac-6.5.12-il66js2n
            mooobxtng3sykhzmdhobuxv4/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/kokkos-4.6.01-gjjssb76dkwnmwjvu3lgslvtdknk6soq/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/mpark-variant-1.4.0-e5n25mavzyoqywk4w6g3xwbx5f7kerwh/include -isyste
            m /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/ports-of-call-main-qhzavjeh5tdrzm55u7uw74j4po3rd6ra/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/spiner-main-kq4og7rxelnbfqef45rzntjw3lds2uen/include -isystem /vast/home/mauneyc/develop/xcap-dev/spack
            /opt/spack/linux-haswell/eigen-3.3.8-withyufrtkjcrt2nng2nauy7jt6qbluw/include/eigen3 -O3 -DNDEBUG -pthread -march=core-avx2 -mtune=core-avx2 -Wno-psabi -Wno-class-memaccess -MD -MT CMakeFiles/singularity-eos_Library.dir/singularity-eos/eos/get_sg_eos_rho_e.cpp.o -MF CMakeFiles/singularity-eos_Library.dir/
            singularity-eos/eos/get_sg_eos_rho_e.cpp.o.d -o CMakeFiles/singularity-eos_Library.dir/singularity-eos/eos/get_sg_eos_rho_e.cpp.o -c /vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos/singularity-eos/eos/get_sg_eos_rho_e.cpp
     46     In file included from /vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos/singularity-eos/eos/default_variant.hpp:27,
     47                      from /tmp/mauneyc/spack-stage/spack-stage-singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc/spack-build-66ruf4l/generated/singularity-eos/eos/eos.hpp:18,
     48                      from /vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos/singularity-eos/eos/singularity_eos.cpp:17:
  >> 49     /vast/home/mauneyc/develop/xcap-dev/dev-env-workflow/singularity-eos/singularity-eos/eos/eos_base.hpp:31:15: error: static assertion failed
     50        31 | static_assert(false);
     51           |               ^~~~~
```

Oops! Looks like we messed up. Let's fix it
```bash
$ seos-dev > vi singularity-eos/eos/eos_bash.hpp
```

```c++
	15   #ifndef _SINGULARITY_EOS_EOS_EOS_BASE_
    16   #define _SINGULARITY_EOS_EOS_EOS_BASE_
    17   
    18   #include <cstring>
    19   #include <limits>
    20   #include <string>
    21   
    22   #include <ports-of-call/portability.hpp>
    23   #include <ports-of-call/portable_errors.hpp>
    24   #include <singularity-eos/base/constants.hpp>
    25   #include <singularity-eos/base/robust_utils.hpp>
    26   #include <singularity-eos/base/root-finding-1d/root_finding.hpp>
    27   #include <singularity-eos/base/variadic_utils.hpp>
    28   
    29   // MY NEW CODE
    30   // THIS WILL FIX EVERYTHING! NOW TRUE!
    31   static_assert(true);
```

```bash
$ seos-dev > spack install
[+] /usr (external glibc-2.28-i5ej3xmy3uonjbh2yn3lubq7xmue442j)
[+] /projects/opt/rhel8/x86_64/openmpi/5.0.2-gcc_13.2.0 (external openmpi-5.0.2-lpanugo33x6swa22ws2gvnpt4cyjd6p4)
[+] /usr (external pkgconf-1.4.2-ojffucqu32xqeuph3dnufdahuiqzjijl)
[+] /usr (external zlib-1.2.11-p3bp3e3qjussu7ibvpfolmjd36va7rkq)
[+] /usr (external curl-7.61.1-ivwo7eocoj4lubj2xj7t4luqhffzkng4)
[+] /usr (external gmake-4.2.1-df7ehnfeehbynnqobw5vg4ih5aojixew)
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/compiler-wrapper-1.0-mj4ubybewkayq3jjmsuyqdyskwffbre5
[+] /projects/opt/rhel8/x86_64/gcc/13.2.0 (external gcc-13.2.0-lbz2bqno2yi74phtrhypxnc3quk3ylbu)
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/gcc-runtime-13.2.0-va64afhaljvjalhfhh66dudde5ifaib4
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/hdf5-1.14.6-zaqgsqxckmewfgc2wdptlmd2yrgpqp54
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/ports-of-call-main-qhzavjeh5tdrzm55u7uw74j4po3rd6ra
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/ncurses-6.5-b4cxilhggxjrp3ly2bwkh42rogdbvg2c
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/kokkos-4.6.01-gjjssb76dkwnmwjvu3lgslvtdknk6soq
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/eospac-6.5.12-il66js2nmooobxtng3sykhzmdhobuxv4
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/mpark-variant-1.4.0-e5n25mavzyoqywk4w6g3xwbx5f7kerwh
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/eigen-3.3.8-withyufrtkjcrt2nng2nauy7jt6qbluw
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/cmake-3.31.6-cvvxqemwg7lkzpokbh3mt4vbkma6ohnm
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/spiner-main-kq4og7rxelnbfqef45rzntjw3lds2uen
==> No binary for singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc found: installing from source
==> Installing singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc [19/19]
==> No patches needed for singularity-eos
==> singularity-eos: Executing phase: 'cmake'
==> singularity-eos: Executing phase: 'build'
==> singularity-eos: Executing phase: 'install'
==> singularity-eos: Successfully installed singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc
  Stage: 0.00s.  Cmake: 0.00s.  Build: 1m 16.82s.  Install: 2.69s.  Post-install: 0.56s.  Total: 1m 20.29s
[+] /vast/home/mauneyc/develop/xcap-dev/spack/opt/spack/linux-haswell/singularity-eos-main-66ruf4le2mplqe3qiux3bowkck3gyknc
```

Great! We can now commit and request a merge.

And there we go! We've gotten a build with our expected environment, made changes and built all with spack!

## Working with `spiner` 

Now, we want to add a feature to `spiner` to enhance `singularity-eos`, but we need concurrent checkouts to update and test these together.

The procedure here is almost identical.
spack develop spiner

NOTE: Here, we are still using our local clone of `singularity-eos`. At the end of this section we will "remove" these local clones from our environment.

ls

vi #...

spack install

And there we go! We have a shared environment for `spiner` and `singularity-eos`, and `singularity-eos` is using our local, updated `spiner`! Who said spack is hard?

Well, we've done enough work for now. Let's clean up while our merges are getting approved

spack undevelop spiner

spack undevelop singularity-eos

## Working with another dependency
Kokkos

spack config change packages:all:variants: +hdf5+mpi+kokkos

Since we are asking for a package that we don't explicitly make a spack spec for, we need to provide it here. Otherwise, this works as you should expect by now

spack develop kokkos@4

vi #...

## Getting a shell with an environment, manually building

spack build-env singularity-eos -- bash

source $SPACK_ROOT/share/spack/setup-env.sh

spack cd -c singulary-eos
