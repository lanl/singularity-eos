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
    version("1.4.2", sha256="02e0ae66cda5bec7c2121d16728b9fbaa653b3863d11737e43b303145ddb2cb2")
    version("1.4.1", sha256="82d2c75fcca8bd613273fd4126749df68ccc22fbe4134ba673b4275f9972b78d")
    version("1.4.0", sha256="e08ae556b7c30d14d77147d248d118cf5343a2e8c0847943385c602394bda0fa")
    version('1.3.0', sha256='54b4a62539c23b1a345dd87c1eac65f4f69db4e50336cd81a15a627ce80ce7d9')
    version('1.2.0', sha256='b802ffa07c5f34ea9839f23841082133d8af191efe5a526cb7e53ec338ac146b')
    version('1.1.0', sha256='c47f7e24c82176b69229a2bcb23a6adcf274dc90ec77a452a36ccae0b12e6e39')

    depends_on("cmake@3.12:")
