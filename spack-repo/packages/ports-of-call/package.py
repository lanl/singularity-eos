# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class PortsOfCall(CMakePackage):
    """Ports of Call: Performance Portability Utilities"""

    homepage    = "https://github.com/lanl/ports-of-call"
    url         = "https://github.com/lanl/ports-of-call/archive/refs/tags/v1.1.0.tar.gz"
    git         = "https://github.com/lanl/ports-of-call.git"

    maintainers = ['rbberger']

    version("main", branch="main")
    version('1.3.0', sha256='54b4a62539c23b1a345dd87c1eac65f4f69db4e50336cd81a15a627ce80ce7d9')
    version('1.2.0', sha256='b802ffa07c5f34ea9839f23841082133d8af191efe5a526cb7e53ec338ac146b')
    version('1.1.0', sha256='c47f7e24c82176b69229a2bcb23a6adcf274dc90ec77a452a36ccae0b12e6e39')

    depends_on("cmake@3.12:")
