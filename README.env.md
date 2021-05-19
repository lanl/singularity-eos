# Spack Environment Demo
This document will go over setting up the simple, ideal spack configuration, and how to build and activate the development environment for singularity.

Spack is designed to seperate the abstract dependency graph from the particular depolyment environment. That is, a project will have dependencies that are generally independent of compilers or hardware. Therefore, spack makes use of two configurations: an abstract configuration for a package and/or environment that only represents the dependency graph, and a configuration for the instance where spack is locally installed. The local configuration defines things like what and where compilers are, where packages can be found, ect.

The local configuration has to be maintained by the local user. However, the abstract configuration (in this case, the environment) does not need to bother with these details. In the repo, we define the abstract environment, and (hopefully!) spack will be able to take the local instance configuration and concretize the abstract environment.

NOTE: Moving forward we assume there is no local spack configuration installed. If you do have `~/.spack` directory already present, you should rename it for the purposes of this demonstration. Later, you can move it back to maintain your original settings. Further, we will use a clean spack clone, though in this case you do not need to move your previous instance.

## Preliminaries
### Clone the Spack repository
I clone it to `${HOME}/spack` for convience, but the location does not matter. Moving forward, I will refer to this directory as `${SPACK_TOP}` (this env var is only for reference, and does not gets set in the process)

    cd
    git clone https://github.com/spack/spack.git 

### Do initial spack configuration
Here we will generate a base spack configuration for your local machine. There are ways and there are ways to do this, but we will try to keep it as simple as possible.

##### Step 0: (optional)
In the following steps, Spack will try to look for compilers and packages it can automatically generate configurations for. It can be somewhat clever in where it looks, but if you want it to find some particular package or compiler, you're best bet is to bring it into your `${PATH}`, for instance by doing a `module load X` command.

##### Step 1: activate spack
This step will bring in all the spack machinary into your path.

    . ${SPACK_TOP}/share/spack/setup-env.sh

##### Step 2: Detect compilers and packages
This step will add a local configuration to `${HOME}/.spack`, specific to your workspace machine. 

    spack compiler find
    spack external find --not-buildable

Each command will generate output telling you what was found, and where the configuration was placed.
Note that this process is not perfect. Some packages may be missed, or their specific versions will not be usable as a dependency for later concretization. It is possible for spack to find no compilers or pre-installed packages, spack will just build everything (even compilers). This step, however, hopefully alievates some unnecessary work.

## Using the Dev environment
In the top directory, the is a new file `spack.yaml`:

```yaml
#
# spack config for environment [singularity_eos-env]
#
spack:
  specs:
  - hdf5~mpi+cxx+hl
  - cmake
  - ninja
  concretization: together
  view: true
  mirrors:
    lanl_spack_mirror: https://pe-serve.lanl.gov/spack-mirror/
```
The important parts for now are the entries under `specs:`. These are the minimum required dependencies to build that are currently not in the tree. For the purposes of this demo, we will elide the syntax of these entires for now.

##### Step 0: (optional)
Generally, your compilers are already in path. However, if you want to have spack default to a particular one, run `spack load <compiler spec>`. You can see the compiler specs in `${HOME}/.spack/compilers.yaml`. Example: `spack load gcc@10.2.0`

##### Step 1: Activate development environment
Go to the top directory of the source code. To activate this environment, invoke

`spack env activate .`

The env is active, but currently is just abstract. We need to concretize and install it.

_NOTE: If you have previously done Steps 2-3, go ahead and skip to 4_

##### Step 2: Concretize

This is the point where the abstract dependency graph and the local computer environment get worked.

`spack concretize -f`

(The `-f` option is usually redundent, but forces a concretization even if concretized before)

You should see similar output to

```bash
==> Concretized hdf5+cxx+hl~mpi
[+]  xnchygx  hdf5@1.10.7%gcc@11.1.0+cxx~debug~fortran+hl~java~mpi+pic+shared~szip~threadsafe api=none arch=linux-manjaro21-skylake
[+]  m2ayh7e      ^zlib@1.2.11%gcc@11.1.0+optimize+pic+shared arch=linux-manjaro21-skylake

==> Concretized cmake
[+]  vg5az72  cmake@3.20.2%gcc@11.1.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-manjaro21-skylake
[+]  mgmse5w      ^ncurses@6.2%gcc@11.1.0~symlinks+termlib abi=none arch=linux-manjaro21-skylake
[+]  qaeny5m          ^pkgconf@1.7.4%gcc@11.1.0 arch=linux-manjaro21-skylake
[+]  i4ilzdw      ^openssl@1.1.1k%gcc@11.1.0~docs+systemcerts arch=linux-manjaro21-skylake
[+]  gpy6hwv          ^perl@5.32.1%gcc@11.1.0+cpanm+shared+threads arch=linux-manjaro21-skylake
[+]  rgn2eot              ^berkeley-db@18.1.40%gcc@11.1.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-manjaro21-skylake
[+]  pppq4ti              ^gdbm@1.19%gcc@11.1.0 arch=linux-manjaro21-skylake
[+]  lie5mqn                  ^readline@8.1%gcc@11.1.0 arch=linux-manjaro21-skylake
[+]  m2ayh7e          ^zlib@1.2.11%gcc@11.1.0+optimize+pic+shared arch=linux-manjaro21-skylake

==> Concretized ninja
[+]  lld5sf3  ninja@1.10.2%gcc@11.1.0 arch=linux-manjaro21-skylake
 -   27xp5tl      ^python@3.8.10%gcc@11.1.0+bz2+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93189bc278fbc37a50ed7f183bd8aaf249a8e1670a465f0db6bb4f8cf87 arch=linux-manjaro21-skylake
 -   r26sjzs          ^bzip2@1.0.8%gcc@11.1.0~debug~pic+shared arch=linux-manjaro21-skylake
 -   db32yy6              ^diffutils@3.7%gcc@11.1.0 arch=linux-manjaro21-skylake
 -   caudqfa                  ^libiconv@1.16%gcc@11.1.0 arch=linux-manjaro21-skylake
 -   j5y6q6y          ^expat@2.3.0%gcc@11.1.0+libbsd arch=linux-manjaro21-skylake
 -   dqrm4gc              ^libbsd@0.11.3%gcc@11.1.0 arch=linux-manjaro21-skylake
 -   43bpbme                  ^libmd@1.0.3%gcc@11.1.0 arch=linux-manjaro21-skylake
[+]  pppq4ti          ^gdbm@1.19%gcc@11.1.0 arch=linux-manjaro21-skylake
[+]  lie5mqn              ^readline@8.1%gcc@11.1.0 arch=linux-manjaro21-skylake
[+]  mgmse5w                  ^ncurses@6.2%gcc@11.1.0~symlinks+termlib abi=none arch=linux-manjaro21-skylake
[+]  qaeny5m                      ^pkgconf@1.7.4%gcc@11.1.0 arch=linux-manjaro21-skylake
 -   hgssc66          ^gettext@0.21%gcc@11.1.0+bzip2+curses+git~libunistring+libxml2+tar+xz arch=linux-manjaro21-skylake
 -   t45uko3              ^libxml2@2.9.10%gcc@11.1.0~python arch=linux-manjaro21-skylake
 -   2ef6f7e                  ^xz@5.2.5%gcc@11.1.0~pic arch=linux-manjaro21-skylake
[+]  m2ayh7e                  ^zlib@1.2.11%gcc@11.1.0+optimize+pic+shared arch=linux-manjaro21-skylake
 -   lghvo4m              ^tar@1.34%gcc@11.1.0 arch=linux-manjaro21-skylake
 -   il37yob          ^libffi@3.3%gcc@11.1.0 patches=26f26c6f29a7ce9bf370ad3ab2610f99365b4bdd7b82e7c31df41a3370d685c0 arch=linux-manjaro21-skylake
[+]  i4ilzdw          ^openssl@1.1.1k%gcc@11.1.0~docs+systemcerts arch=linux-manjaro21-skylake
[+]  gpy6hwv              ^perl@5.32.1%gcc@11.1.0+cpanm+shared+threads arch=linux-manjaro21-skylake
[+]  rgn2eot                  ^berkeley-db@18.1.40%gcc@11.1.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-manjaro21-skylake
 -   3gba3da          ^sqlite@3.35.4%gcc@11.1.0+column_metadata+fts~functions~rtree arch=linux-manjaro21-skylake
 -   sjfhl7l          ^util-linux-uuid@2.36.2%gcc@11.1.0 arch=linux-manjaro21-skylake

==> Updating view at /home/mauneyc/scratch/singularity_local/singularity-eos/.spack-env/view
```
The first column indicates if the package spec is already installed (`[+]` if found), the next column is the short hash of the spec, followed by the full string spec. Indented specs with the `^` prefix indicate nodes/leaves of the dependency graph.

##### Step 3: Install

We're ready to install the environment

`spack install`

This may take some time, depending on what you have installed already.
Spack will install these files into the local spack instance at `${SPACK_TOP}`, and by default will generate module files for them. This way, these packages can be used by other environments if needed.

##### Step 4: Ready to build

We are now ready to build using the new environment.

```
cd <build-dir>
cmake -GNinja <path-to-repo>
-- The C compiler identification is GNU 11.1.0
-- The CXX compiler identification is GNU 11.1.0
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /home/mauneyc/spack/opt/spack/linux-manjaro21-skylake/gcc-10.2.0/gcc-11.1.0-ewakllpu44azdc5g2k44bzwc4kynznhp/bin/gcc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /home/mauneyc/spack/opt/spack/linux-manjaro21-skylake/gcc-10.2.0/gcc-11.1.0-ewakllpu44azdc5g2k44bzwc4kynznhp/bin/g++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- The Fortran compiler identification is GNU 11.1.0
-- Detecting Fortran compiler ABI info
-- Detecting Fortran compiler ABI info - done
-- Check for working Fortran compiler: /home/mauneyc/spack/opt/spack/linux-manjaro21-skylake/gcc-10.2.0/gcc-11.1.0-ewakllpu44azdc5g2k44bzwc4kynznhp/bin/gfortran - skipped
-- Checking whether /home/mauneyc/spack/opt/spack/linux-manjaro21-skylake/gcc-10.2.0/gcc-11.1.0-ewakllpu44azdc5g2k44bzwc4kynznhp/bin/gfortran supports Fortran 90
-- Checking whether /home/mauneyc/spack/opt/spack/linux-manjaro21-skylake/gcc-10.2.0/gcc-11.1.0-ewakllpu44azdc5g2k44bzwc4kynznhp/bin/gfortran supports Fortran 90 - yes
-- Patching mpark::variant to support GPUs
Reversed (or previously applied) patch detected!  Skipping patch.
12 out of 12 hunks ignored -- saving rejects to file /home/mauneyc/scratch/singularity_local/singularity-eos/utils/variant/include/mpark/variant.hpp.rej
-- Found HDF5: /home/mauneyc/scratch/singularity_local/singularity-eos/.spack-env/view/lib/libhdf5.so;/home/mauneyc/scratch/singularity_local/singularity-eos/.spack-env/view/lib/libz.so;/usr/lib/libdl.so;/usr/lib/libm.so (found version "1.10.7") found components: C HL 
-- Found Eigen3:   
-- Found PortsofCall:   
-- Configuring done
-- Generating done
-- Build files have been written to: /home/mauneyc/scratch/singularity_local/builds.singularity/b1
```

Notice that we have now found the HDF5 installed by the spack environment (in the above, it also used a `gcc` that I had installed using spack).

## Cleanup

To deactivate the environment, simply

`spack env deactivate`

which will return your shell to a state prior to activation.
