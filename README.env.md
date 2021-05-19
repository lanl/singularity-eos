# Spack Environment Demo
*UPDATE*: The pure environment has been put on the side for now. For now, we use a "dependency package" and repos.

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
In the top directory, the is a folder "spack-repo"

##### Step 0: (optional)
Generally, your compilers are already in path. However, if you want to have spack default to a particular one, run `spack load <compiler spec>`. You can see the compiler specs in `${HOME}/.spack/compilers.yaml`. Example: `spack load gcc@10.2.0`

##### Step 1: Activate development environment
Go to the top directory of the source code. To add the repo,

`spack repo add spack-repo/`

##### Step 2: Choose variants (WIP) and Install

As work continues, environment specificity can be expressed using package variants, for instance `singularity-eos_deps+use_cuda`. This will be more documeneted when there is more done on this aspect.

##### Step 4: Load the package and build

Once installed, loading the package will load all the dependencies with it. This works similiar to `module load <...>`

`spack load singularity-eos_deps`

Note that, once there are more variants, you can install/load different concretizations of the package

`spack load singularity-eos_deps+use_cuda`

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

To unload, simply

`spack unload singularity-eos_deps`

which will return your shell to a state prior to loading.
