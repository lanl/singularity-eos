from kessel.workflows import *
from kessel.workflows.base.spack import BuildEnvironment
from kessel.workflows.base.cmake import CTest as CTestWorkflow
from kessel.workflows.cmake import Build


class Cdash(Build, CTestWorkflow):
    steps = ["env", "configure", "build", "test", "install", "submit"]
    ctest_project_name = environment("singularity-eos", variable="CTEST_PROJECT_NAME")
