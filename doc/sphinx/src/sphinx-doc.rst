.. _sphinx-doc:

.. _Sphinx CheatSheet: https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html

How to Use Sphinx for Writing Docs
===================================

How to Get the Dependencies
---------------------------

Using Docker
^^^^^^^^^^^^

If you are using `Docker`_, then simply pull the docker image specified below:

.. _Docker: https://www.docker.com

.. code-block::

  image: sphinxdoc/sphinx-latexpdf

Then, after running :code:`docker run -it <docker-image-name> /bin/bash`, install the theme we are using with :code:`pip install sphinx_rtd_theme`

Using Spack
^^^^^^^^^^^

If you are using `Spack`_ to provision dependencies, then follow the steps as such:

.. _Spack: https://spack.io

.. literalinclude:: ../../../.gitlab-ci.yml
   :lineno-match:
   :language: yaml
   :lines: 111-120

from :code:`.gitlab-ci.yml`

.. warning::
   If you do not have either Docker or Spack locally, you would need to install one of them first.

   For Docker, refer to their `Get Docker Guide`_.

   For Spack, refer to their `Getting Started Guide`_.

.. _Get Docker Guide: https://docs.docker.com/get-docker

.. _Getting Started Guide: https://spack.readthedocs.io/en/latest/getting_started.html#installation

How to Build .rst into .html
----------------------------

After you have the dependencies in your environment, then simply build your documentation as the following:

.. literalinclude:: ../../../.gitlab-ci.yml
   :lineno-match:
   :language: yaml
   :lines: 121-122

from :code:`.gitlab-ci.yml`

.. note:: 
   You can view the documentation webpage locally on your web browser by passing in the URL as :code:`file:///path/to/singularity-eos/doc/sphinx/_build/html/index.html`

How to Deploy
-------------

#. Submit a PR with your .rst changes for documentation on `Github Singularity-EOS`_
#. Get your PR reviewed and merged into main
#. Make sure the :code:`pages` CI job passes in the CI pipeline

.. _Github Singularity-EOS: https://github.com/lanl/singularity-eos

As soon as the PR is merged into main, this will trigger the Pages deployment automatically if the :code:`pages` CI job passes.

Documentation is available on `github-pages`_ and on `re-git`_

.. _github-pages: https://lanl.github.io/singularity-eos/

.. _re-git: http://xcap.re-pages.lanl.gov/oss/singularity-eos

More Info.
----------

* `Sphinx Installation`_

.. _Sphinx Installation: https://www.sphinx-doc.org/en/master/usage/installation.html

* `Sphinx reStructuredText Documentation`_

.. _Sphinx reStructuredText Documentation: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
