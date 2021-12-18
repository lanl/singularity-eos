# Spackage for Singularity-EOS

from spack import *

class SingularityEos(CMakePackage, CudaPackage):
    git         = "http://github.com/lanl/singularity-eos.git"

    version("main", branch="main", submodules=True)

    variant("kokkos",
            description="Enable kokkos",
            default=False
    )

    variant("kokkos-kernels",
            description="Enable kokkos-kernals for linear algebra",
            default=False
    )

    variant("openmp",
            description="Enable openmp",
            default=False
    )

    variant("build_extra",
            description="Build converters",
            values=any_combination_of(
                'sesame','stellarcollapse'
            ).with_default('none') 	
    )

    variant("enable_tests",
            default=False,
            description="Build tests"
    )

    variant("enable_fortran",
            default=True,
            description="Enable building fortran interface"
    )

    variant("mpi",
            default=False,
            description="Build with MPI support"
    )

    variant("doc",
            default=False,
            description="Enable Sphinx Documentation Support"
    )

    #depends_on("mpark-variant")
    depends_on("hdf5+cxx+hl~mpi", when="~mpi")
    depends_on("hdf5+cxx+hl+mpi", when="+mpi")
    depends_on("eospac")
    depends_on("py-h5py")

    depends_on("cmake@3.12:")
    depends_on("eigen@3.3.9", when="~kokkos-kernels")
    depends_on("catch2@2.13.4:2.13.6", when="+enable_tests")

    depends_on("py-sphinx", when="+doc")
    depends_on("py-sphinx-rtd-theme@0.4.3", when="+doc")
    depends_on("py-sphinx-multiversion", when="+doc")

    for _flag in ("~cuda", "+cuda", "~openmp", "+openmp"):
        depends_on("kokkos@3.3:" +_flag, when="+kokkos" + _flag)
        depends_on("kokkos-kernels" + _flag, when="+kokkos-kernels" + _flag)

    conflicts("+cuda", when="~kokkos")
    conflicts("+openmp", when="~kokkos")
    conflicts("+kokkos-kernels", when="~kokkos")

    depends_on("kokkos@3.3:~shared+wrapper+cuda_lambda+cuda_relocatable_device_code", when="+cuda+kokkos")
    depends_on("kokkos-nvcc-wrapper~mpi", when="+cuda+kokkos~mpi")
    depends_on("kokkos-nvcc-wrapper+mpi", when="+cuda+kokkos+mpi")
    for _flag in list(CudaPackage.cuda_arch_values):
        depends_on("kokkos@3.3: cuda_arch=" +_flag, when="+cuda+kokkos cuda_arch=" + _flag)
        depends_on("kokkos-kernels cuda_arch=" +_flag, when="+cuda+kokkos cuda_arch=" + _flag)

    conflicts("cuda_arch=none", when="+cuda",
          msg="CUDA architecture is required")

    def cmake_args(self):

        args = [
            self.define_from_variant("SINGULARITY_USE_CUDE", "cuda"),
            self.define_from_variant("SINGULARITY_USE_KOKKOS", "kokkos"),
            self.define_from_variant("SINGULARITY_USE_KOKKOSKERNELS", "kokkos-kernels"),
            self.define_from_variant("SINGULARITY_USE_FORTRAN", "enable_fortran"),
            self.define_from_variant("SINGULARITY_BUILD_CLOSURE", "enable_fortran"),
            self.define_from_variant("SINGULARITY_BUILD_TESTS", "enable_tests"),
            self.define("SINGULARITY_BUILD_SESAME2SPINER", "sesame" in self.spec.variants["build_extra"]),
            self.define("SINGULARITY_TEST_SESAME", ("sesame" in self.spec.variants["build_extra"] and "enable_tests" in self.spec)),
            self.define("SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER", "stellarcollapse" in self.spec.variants["build_extra"]),
            self.define("SINGULARITY_TEST_STELLARCOLLAPSE2SPINER", ("stellarcollapse" in self.spec.variants["build_extra"] and "enable_tests" in self.spec)),
            self.define("SINGULARITY_USE_HDF5", True),
            self.define("SINGULARITY_USE_EOSPAC", True)
        ]

        return args
