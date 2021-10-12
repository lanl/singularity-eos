# dependency package for singulary-eos

import os
from spack import *

class SingularityEosDeps(BundlePackage):
    homepage    = "https://github.com/lanl/singularity-eos"
    url         = "https://github.com/lanl/singularity-eos/archive/refs/heads/main.zip"
    git         = "git@github.com:lanl/singularity-eos.git"

    version("main", branch="main")

    variant(
            "cuda",
            default=False,
            description="Use cuda"
    )

    variant(
            "kokkos",
            default=False,
            description="Use kokkos"
    )

    variant(
            "linalg", 
            default="eigen", 
            description="Linear algebra library",
            values=("eigen" ,"kokkos"),
            multi=False
    )
    variant(
            "build_extra",
            values=("sesame","stellarcollapse", "none"), 	
            default="none",
            description="Build converters",
            multi=True
    )

    variant(
            "enable_tests",
            default=False,
            description="Build tests"
    )

    depends_on("hdf5~mpi+cxx+hl")
    depends_on("cuda@11:", when="+cuda")
    depends_on("kokkos@3:", when="+kokkos")
    depends_on("eigen@3.3.9", when="linalg=eigen")
    depends_on("kokkos-kernels", when="linalg=kokkos")
    depends_on("mpark-variant")
    depends_on("eospac")

    depends_on("cmake@3.12:")
    depends_on("catch2@2.13.4:2.13.6")

    conflicts(
            "linalg=eigen",
            when="+cuda"
    )
    

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
            if spec.variants["linalg"].value == "kokkos":
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

        with working_dir('spack-build', create=True):
            with open("singularity_tc.cmake", 'w') as cmtcf:
                for k,v in cmake_args_map.items():
                    cmtcf.write("set(" +k + " " + v + " CACHE BOOL \"\")\n")
            install("singularity_tc.cmake", join_path(prefix, "singularity_tc.cmake"))
            
        

