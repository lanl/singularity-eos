# Spackage for Singularity-EOS

from spack import *

class SingularityEos(CMakePackage, CudaPackage):
    homepage    = "https://lanl.github.io/singularity-eos/main/index.html"
    git         = "http://github.com/lanl/singularity-eos.git"

    version("main", branch="main", submodules=True)

    variant("kokkos", default=False, description="Enable kokkos")
    variant("kokkos-kernels", default=False, description="Enable kokkos-kernals for linear algebra")
    variant("openmp", default=False, description="Enable openmp")
    variant("mpi", default=False, description="Build with MPI support")
    variant("build_extra", description="Build converters", values=any_combination_of("sesame","stellarcollapse").with_default("none"))

    variant("tests", default=False, description="Build tests")
    variant("fortran", default=True, description="Enable building fortran interface")

    variant("doc", default=False, description="Sphinx Documentation Support")
    variant("format", default=False, description="Clang-Format Support")

    #depends_on("mpark-variant")
    depends_on("eospac")

    depends_on("cmake@3.12:")
    depends_on("eigen@3.3.9", when="~cuda")
    depends_on("catch2@2.13.4:2.13.6", when="+tests")

    depends_on("py-sphinx", when="+doc")
    depends_on("py-sphinx-rtd-theme@0.4.3", when="+doc")
    depends_on("py-sphinx-multiversion", when="+doc")

    depends_on('llvm@12.0.0+clang', when='+format')

    for _flag in ("~cuda", "+cuda", "~openmp", "+openmp"):
        depends_on("kokkos@3.3:" +_flag, when="+kokkos" + _flag)
        depends_on("kokkos-kernels" + _flag, when="+kokkos-kernels" + _flag)

#   NOTE: we can do depends_on("libfoo cppflags='-fPIC -O2'") for compiler options
    conflicts("+cuda", when="~kokkos")
    conflicts("+openmp", when="~kokkos")
    conflicts("+kokkos-kernels", when="~kokkos")

    depends_on("kokkos@3.3:~shared+wrapper+cuda_lambda+cuda_relocatable_device_code", when="+cuda+kokkos")
#    depends_on("kokkos-nvcc-wrapper~mpi", when="+cuda+kokkos~mpi")
#    depends_on("kokkos-nvcc-wrapper+mpi", when="+cuda+kokkos+mpi")

    for _flag in ("~mpi", "+mpi"):
        depends_on("hdf5+cxx+hl" + _flag, when=_flag)
        depends_on("py-h5py" + _flag, when=_flag)
        depends_on("kokkos-nvcc-wrapper" + _flag, when="+cuda+kokkos"+_flag)


    for _flag in list(CudaPackage.cuda_arch_values):
        depends_on("kokkos@3.3: cuda_arch=" +_flag, when="+cuda+kokkos cuda_arch=" + _flag)
        depends_on("kokkos-kernels cuda_arch=" +_flag, when="+cuda+kokkos cuda_arch=" + _flag)

    conflicts("cuda_arch=none", when="+cuda",
          msg="CUDA architecture is required")

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
            self.define("SINGULARITY_USE_HDF5", "+hdf5" in self.spec),
            self.define("SINGULARITY_USE_EOSPAC", "+eospac" in self.spec)
        ]

        return args

    @run_after('cmake')
    def generate_cmake_configuration(self):
        config_fname = "singularity-eos_spackconfig.cmake"
        cmake_config = self.cmake_args()
        with working_dir("cmake-gen", create=True):
            with open(config_fname, "w") as cmc:
                for arg in cmake_config:
                    cmc.write(arg + "\n")
            install(config_fname, join_path(prefix, config_fname))
        return        
