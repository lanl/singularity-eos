from kessel.workflows import *
from kessel.workflows.base.spack import BuildEnvironment
from kessel.workflows.base.cmake import CMake as CMakeWorkflow
from pathlib import Path


class Build(BuildEnvironment):
    steps = ["env", "configure", "build", "test", "install"]

    spack_env = environment("singularity_eos_default")
    project_spec = environment("singularity-eos+tests")

    xcap_spackages_checkout = environment(Path.cwd() / "extern" / "xcap_spackages", variable="XCAP_SPACKAGES_CHECKOUT")

    def ci_message(self, args):
        pre_alloc_init = ""
        post_alloc_init = ""

        if "XCAP_SPACKAGES_MR" in self.environ and self.environ["XCAP_SPACKAGES_MR"]:
            pre_alloc_init += f"export XCAP_SPACKAGES_MR={self.environ['XCAP_SPACKAGES_MR']}\n"

        pre_alloc_init += "source .gitlab/download_prereqs.sh"

        if "DEPLOYMENT_VERSION" in self.environ and \
           "DEPLOYMENT_VERSION_CURRENT_DEFAULT" in self.environ and \
           self.environ["DEPLOYMENT_VERSION"] != self.environ["DEPLOYMENT_VERSION_CURRENT_DEFAULT"]:
            post_alloc_init += f"export DEPLOYMENT_VERSION={self.environ['DEPLOYMENT_VERSION']}\n"

        post_alloc_init += "source .gitlab/kessel.sh"
        return super().ci_message(args, pre_alloc_init=pre_alloc_init, post_alloc_init=post_alloc_init)

    @collapsed
    def env(self, args):
        """Prepare Environment"""
        self.xcap_spackages_checkout = self.source_dir / "extern" / "xcap_spackages"

        # Always allow lockfile changes since we use spiner@main and ports-of-call@main
        self.allow_lockfile_changes = True

        if not self.xcap_spackages_checkout.exists():
            self.exec("source .gitlab/download_prereqs.sh")

        super().prepare_env(args)

        if "KESSEL_DEPLOYMENT" not in self.environ:
            self.exec(f"spack repo add {self.xcap_spackages_checkout}/v2/spack_repo/xcap || true")

        # Handle spiner and ports-of-call submodule dependencies
        spiner_path = self.source_dir / "utils" / "spiner"
        if (spiner_path / ".git").exists():
            self.exec(f"spack develop --no-clone -p {spiner_path} spiner@main")

        ports_of_call_path = self.source_dir / "utils" / "ports-of-call"
        if (ports_of_call_path / ".git").exists():
            self.exec(f"spack develop --no-clone -p {ports_of_call_path} ports-of-call@main")

        super().install_env(args)

    def test(self, args):
        """Testing"""
        self.exec(f"""
            pushd {self.build_dir}
            if [[ -f sesame2spiner/sesame2spiner ]]; then
                echo "/usr/projects/data/eos/eos-developmental/Sn2162/v01/sn2162-v01.bin" > sesameFilesDir.txt
                sesame2spiner/sesame2spiner -s materials.sp5 ../sesame2spiner/examples/unit_tests/*.dat
                sesame2spiner/sesame2spiner -s duplicates.sp5 ../sesame2spiner/examples/duplicate-test/*.dat
            fi
            popd""")
        super().test(args)


class Cmake(Build, CMakeWorkflow):
    pass
