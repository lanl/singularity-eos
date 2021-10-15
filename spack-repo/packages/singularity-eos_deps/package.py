# dependency package for singulary-eos

import os
from spack import *

class SingularityEosDeps(BundlePackage, CudaPackage):
    homepage    = "https://github.com/lanl/singularity-eos"
    url         = "https://github.com/lanl/singularity-eos/archive/refs/heads/main.zip"
    git         = "git@github.com:lanl/singularity-eos.git"

    version("main", branch="main")

    variant("kokkos",
            description="Enable kokkos",
            default=True
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

    depends_on("mpark-variant")
    depends_on("hdf5~mpi+cxx+hl")
    depends_on("eospac")

    depends_on("cmake@3.12:")
    depends_on("eigen@3.3.9", when="~kokkos-kernels")
    depends_on("catch2@2.13.4:2.13.6", when="+enable_tests")

    for _flag in ("~cuda", "+cuda", "~openmp", "+openmp"):
        depends_on("kokkos@3.3:" +_flag, when="+kokkos" + _flag)
        depends_on("kokkos-kernels" + _flag, when="+kokkos-kernels" + _flag)

    with when("~kokkos"):
        conflicts("+cuda")
        conflicts("+openmp")
        conflicts("+kokkos-kernels")

    with when("+cuda+kokkos"):
        depends_on("kokkos@3.3:~shared+wrapper+cuda_lambda+cuda_relocatable_device_code")
        depends_on("kokkos-nvcc-wrapper~mpi")
        for _flag in list(CudaPackage.cuda_arch_values):
            depends_on("kokkos@3.3: cuda_arch=" +_flag, when="cuda_arch=" + _flag)
            depends_on("kokkos-kernels cuda_arch=" +_flag, when="cuda_arch=" + _flag)

    conflicts("cuda_arch=none", when="+cuda",
          msg="CUDA architecture is required")

    phases=["install"]

    def setup_run_environment(self, env):
        env.set('HDF5_ROOT', self.spec['hdf5'].prefix)
        env.set('SINGULARITYEOS_TCF', join_path(self.prefix, "singularity_tc.cmake"))
    
    def install(self, spec, prefix):
        cmake_args_map={
            "SINGULARITY_USE_HDF5": "ON",
            "SINGULARITY_USE_FORTRAN": "ON",
            "SINGULARITY_USE_KOKKOS": "OFF",
            "SINGULARITY_USE_EOSPAC": "OFF",
            "SINGULARITY_USE_CUDA": "OFF",
            "SINGULARITY_USE_KOKKOSKERNELS": "OFF",
            "SINGULARITY_BUILD_CLOSURE": "ON",
            "SINGULARITY_BUILD_TESTS": "OFF",
            "SINGULARITY_TEST_SESAME": "OFF",
            "SINGULARITY_TEST_STELLAR_COLLAPSE": "OFF",
            "SINGULARITY_BUILD_SESAME2SPINER": "OFF",	
            "SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER": "OFF",
        }

        if "+cuda" in spec:
            cmake_args_map["SINGULARITY_USE_CUDA"] = "ON"
        if "+kokkos" in spec:
            cmake_args_map["SINGULARITY_USE_KOKKOS"] = "ON"
            if "+kokkos-kernels" in spec:
                cmake_args_map["SINGULARITY_USE_KOKKOSKERNELS"] = "ON"
        
        if "sesame" in spec.variants["build_extra"]:
            cmake_args_map["SINGULARITY_BUILD_SESAME2SPINER"] = "ON"
            if "+enable_tests" in spec:
                cmake_args_map["SINGULARITY_TEST_SESAME"] = "ON"
        if "stellarcollapse" in spec.variants["build_extra"]:
            cmake_args_map["SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER"] = "ON"
            if "+enable_tests" in spec:
                cmake_args_map["SINGULARITY_TEST_STELLAR_COLLAPSE"] = "ON"

        if "+enable_tests" in spec:
            cmake_args_map["SINGULARITY_BUILD_TESTS"] = "ON"

        if "+enable_fortran" in spec:
            cmake_args_map["SINGULARITY_USE_FORTRAN"] = "ON"

        with working_dir('spack-build', create=True):
            with open("singularity_tc.cmake", 'w') as cmtcf:
                for k,v in cmake_args_map.items():
                    cmtcf.write("set(" +k + " " + v + " CACHE BOOL \"\")\n")
            install("singularity_tc.cmake", join_path(prefix, "singularity_tc.cmake"))
            
        

