# dependency package for singulary

import os
from spack import *

class SingularityDeps(BundlePackage):
    homepage    = "https://github.com/lanl/singularity-eos"
    url         = "https://github.com/lanl/singularity-eos/archive/refs/heads/main.zip"
    git         = "git@github.com:lanl/singularity-eos.git"

    version("main", branch="main")

    variant("use_cuda", default=False, description="Enable cuda support")

    depends_on("cmake", type="build")
    depends_on("hdf5~mpi+cxx+hl", type=("build", "run"))
    depends_on("eigen@3.3.9", type="build")
    depends_on("kokkos@3:", when="+use_cuda", type=("build", "run"))
    depends_on("kokkos-kernels", when="+use_cuda", type=("build", "run"))
    depends_on("eospac", type=("build","run"))

    phases=["install"]

    def setup_run_environment(self, env):
        env.set('HDF5_ROOT', self.spec['hdf5'].prefix)
        env.set('SINGULARITY_HELP', os.path.join(self.spec.prefix, f"load_env-{self.spec.full_hash(length=4)}.sh"))

    def install(self, spec, prefix):
        mod_script = os.path.join(spec.prefix, f"load_env-{spec.full_hash(length=4)}.sh")
       
        with open(os.path.join(mod_script), "w") as f:
            f.write(f"# load env {spec.short_spec}")
            f.write("")
            for dep in spec.dependencies(deptype="build"):
                f.write(dep.format())


