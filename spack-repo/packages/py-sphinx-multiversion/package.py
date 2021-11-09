# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

### This is a copy of the spackage upstream, but it is not yet available in the release tagged Spack ###

class PySphinxMultiversion(PythonPackage):
    """A Sphinx extension for building self-hosted versioned documentation."""

    homepage = "https://github.com/Holzhaus/sphinx-multiversion"
    url      = "https://github.com/Holzhaus/sphinx-multiversion/archive/refs/tags/v0.2.4.tar.gz"
    #pypi = "sphinx-multiversion/sphinx-multiversion-0.2.4.tar.gz"

    version('0.2.4', sha256='c22abc33160c8ff63b95bca6df7bffea6c418decfa0456b9645e2e93b8b8a99a')
    version('0.2.3', sha256='a9c72cb3e994f61f3e5a5332e633484863a866a25d9fde9df01ed58a996043db')
    version('0.2.2', sha256='a5af32563ba9261e17626d9a839632bec62d6a4597d1736d50432a7eee02b38c')
    version('0.2.1', sha256='468aa861bc606144b1d5212bb7d5e9b83b1de697dffd3d1737708a203ff7ef59')

    depends_on('py-setuptools', type='build')
    depends_on('py-sphinx', type=('build', 'run'))
