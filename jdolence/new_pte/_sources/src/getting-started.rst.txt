.. _getting-started:

Getting Started
===============

At it's most basic, you can download and install ``singularity-eos`` with:

.. code-block:: bash

  git clone --recursive git@github.com:lanl/singularity-eos.git
  cd singularity-eos
  mkdir bin
  cd bin
  cmake ..
  make install

This will downloady singularity-eos with no optional dependencies and
compile the capabilities available in that form.

If the library is in your include and lib paths (or you built it
in-tree), you can include the eos part of the library with

.. code-block:: cpp

  #include<singularity-eos/eos/eos.hpp>

You can then initialize an equation of state by setting the model
class you want to an ``EOS`` object. For example:

.. code-block:: cpp

  EOS ideal = IdealGas(gm1, Cv);

To see which equations of state are available, see :ref:`The Equation of State API`.

Some equations of state store tabulated data.

.. warning::
  If you want to run one of these on an accelerator device like a GPU,
  you must copy the data over. This is provided by the function

  .. cpp:function:: void EOS::GetOnDevice()

  which can be called as, e.g.,

  .. code-block:: cpp

    ideal.GetOnDevice();

  If you don't want to use GPU's you can ignore this.

You can then call the EOS. For example, to get pressure from density
and temperature, just call (for example),

.. code-block:: cpp

  ideal.PressureFromDensityTemperature(rho, T);

The units are all in cgs. You can ignore the lambda for now.

.. warning::

  When you're done with the model, it's good practice to release
  device memory with a call to

  .. cpp:function:: EOS::Finalize();

  If you're not using device memory, you can ignore this.

And that's it!

Going Deeper
--------------

* You can find code examples in the ``examples`` source directory.
* To learn more about the design philosophy, look :ref:`here <philosophy>`.
* To learn about how to build, look at :ref:`our build document <building>`.
* To learn more about the equation of state API, look :ref:`here <using-eos>`.
* To learn about the available equations of state, look :ref:`here <models>`.
* To learn about our mixed-cell closure models, such as pressure-temperature equilibrium, look at :ref:`using-closures`.
* If you're interested in contributing, check out our :ref:`documentation for developers <contributing>`.

