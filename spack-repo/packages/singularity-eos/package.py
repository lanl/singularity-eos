# Spackage for Singularity-EOS

import os
from spack import *
import llnl.util.tty as tty

class SingularityEos(CMakePackage, CudaPackage):
    homepage    = "https://lanl.github.io/singularity-eos/main/index.html"
    git         = "http://github.com/lanl/singularity-eos.git"

    version("main", branch="main", submodules=True)
    version("develop", branch="mauneyc/update/cmake-idiomatic-install", submodules=True)

    # build with kokkos, kokkos-kernels for offloading support
    variant("kokkos", default=False, description="Enable kokkos")
    variant("kokkos-kernels", default=False, description="Enable kokkos-kernals for linear algebra")

    variant("openmp", default=False, description="Enable openmp")

    # for compatibility with downstream projects
    variant("mpi", default=False, description="Build with MPI support")

    # build converters for sesame, stellarcollapse eos's
    variant("build_extra", description="Build converters", values=any_combination_of("sesame","stellarcollapse").with_default("none"))

    # build tests (NB: will include tests for selected options in `build_extra`
    variant("tests", default=False, description="Build tests")

    # build the fortran interface
    variant("fortran", default=True, description="Enable building fortran interface")

    # build the documentation
    variant("doc", default=False, description="Sphinx Documentation Support")

    # include depedencies for automatic code formatting (i.e. clang-format)
    variant("format", default=False, description="Clang-Format Support")

    # TODO: do we always depend on eospac? 
    depends_on("eospac")

    # building/testing/docs
    depends_on("cmake@3.14:")
    depends_on("catch2@2.13.7", when="+tests")
#    depends_on("py-h5py", when="+tests build_extra=stellarcollapse")
    depends_on("py-sphinx", when="+doc")
    depends_on("py-sphinx-rtd-theme@0.4.3", when="+doc")
    depends_on("py-sphinx-multiversion", when="+doc")
    # TODO: this can be messy, esp if all we need is clang-format
    depends_on('llvm@12.0.0+clang', when='+format')

    # linear algebra when not using GPUs
    depends_on("eigen@3.3.8", when="~cuda")

    # set up kokkos offloading dependencies
    for _flag in ("~cuda", "+cuda", "~openmp", "+openmp"):
        depends_on("kokkos@3.5.00 ~shared" +_flag, when="+kokkos" + _flag)
        depends_on("kokkos-kernels@3.5.00" + _flag, when="+kokkos-kernels" + _flag)

    # specfic specs when using GPU/cuda offloading
    depends_on("kokkos +wrapper+cuda_lambda+cuda_relocatable_device_code", when="+cuda+kokkos")

    # fix for older spacks
    if spack.version.Version(spack.spack_version) >= spack.version.Version("0.17"):
        depends_on("kokkos-kernels ~shared", when="+kokkos-kernels")

    for _flag in list(CudaPackage.cuda_arch_values):
        depends_on("kokkos cuda_arch=" +_flag, when="+cuda+kokkos cuda_arch=" + _flag)
        depends_on("kokkos-kernels cuda_arch=" +_flag, when="+cuda+kokkos cuda_arch=" + _flag)

    conflicts("cuda_arch=none", when="+cuda",
          msg="CUDA architecture is required")

    # NOTE: we can do depends_on("libfoo cppflags='-fPIC -O2'") for compiler options

    # these are mirrored in the cmake configuration
    conflicts("+cuda", when="~kokkos")
    conflicts("+openmp", when="~kokkos")
    conflicts("+kokkos-kernels", when="~kokkos")

    # NOTE: these are set so that dependencies in downstream projects share common MPI dependence
    for _flag in ("~mpi", "+mpi"):
        depends_on("hdf5~cxx+hl" + _flag, when=_flag)
        depends_on("py-h5py" + _flag, when="+tests build_extra=stellarcollapse "+_flag)
#        depends_on("hdf5+hl" + _flag, when=_flag)
        depends_on("py-h5py" + _flag, when=_flag)
        depends_on("kokkos-nvcc-wrapper" + _flag, when="+cuda+kokkos"+_flag)

    def cmake_args(self):

        args = [
            self.define_from_variant("SINGULARITY_USE_CUDA", "cuda"),
            self.define_from_variant("SINGULARITY_USE_KOKKOS", "kokkos"),
            self.define_from_variant("SINGULARITY_USE_KOKKOSKERNELS", "kokkos-kernels"),
            self.define_from_variant("SINGULARITY_USE_FORTRAN", "fortran"),
            self.define_from_variant("SINGULARITY_BUILD_CLOSURE", "fortran"),
            self.define_from_variant("SINGULARITY_BUILD_TESTS", "tests"),
            self.define("SINGULARITY_BUILD_SESAME2SPINER", "sesame" in self.spec.variants["build_extra"]),
            self.define("SINGULARITY_TEST_SESAME", ("sesame" in self.spec.variants["build_extra"] and "tests" in self.spec)),
            self.define("SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER", "stellarcollapse" in self.spec.variants["build_extra"]),
            self.define("SINGULARITY_TEST_STELLARCOLLAPSE2SPINER", ("stellarcollapse" in self.spec.variants["build_extra"] and "tests" in self.spec)),
            self.define("SINGULARITY_USE_HDF5", "^hdf5" in self.spec),
            self.define("SINGULARITY_USE_EOSPAC", "^eospac" in self.spec)
        ]

        return args

    # specify the name of the auto-generated cmake cache config
    @property
    def cmake_config_fname(self):
        return "singularity-eos_spackconfig.cmake"

    # generate the pre-configured cmake cache file that reflects the spec options
    # NOTE: this file isn't replaced if the same spec is already installed - you may need to uninstall the old spec first
    @run_after('cmake')
    def generate_cmake_configuration(self):
        config_fname = self.cmake_config_fname
        cmake_config = self.cmake_args()

        with working_dir("cmake-gen", create=True):
            with open(config_fname, "w") as cmc:
                for arg in cmake_config:
                    kt, v = arg.replace("-D","").split("=")
                    k, t = kt.split(":")
                    cmc.write("set(" + k + " \"" + v + "\" CACHE " + t + " \"\" FORCE)" + "\n")
            install(config_fname, join_path(prefix, config_fname))
        return       

    # run when loaded
    # NOTE: to use:
    #   cmake -C $SINGULARITY_SPACK_CMAKE_CONFIG ...
    def setup_run_environment(self, env):
        env.set("SINGULARITY_SPACK_CMAKE_CONFIG", os.path.join(self.prefix, self.cmake_config_fname))
