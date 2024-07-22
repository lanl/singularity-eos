.. _examples:

Examples
=========

The ``examples`` directory of ``singularity-eos`` contains several
examples of using the code. You can build the examples by setting
``-DSINGULARITY_BUILD_EXAMPLES=ON`` at ``CMake`` configuration
time. For example:

.. code-block:: bash

   # from singularity-eos repo
   mkdir -p build && cd build
   cmake .. -DSINGULARITY_BUILD_EXAMPLES=ON ..
   make -j

The available examples are listed below.

Get Sound Speed and Pressure
------------------------------

The ``examples/get_sound_speed_press.cpp`` file implements a call go
``singularity-eos`` that computes pressure and sound speed from
density and energy for an ideal gas equation of state. It demonstrates
how to make this call both through individual calls to pressure and
bulk modulus, as well as by calling the in-one ``FillEos`` API. The
former looks something like this:

.. code-block:: cpp

  // Loop through the cells and use the two function calls
  for (int i = 0; i < Ncells; ++i) {
    double sie = robust::ratio(uu[i], rho[i]); // convert to specific internal energy
    P[i] = eos.PressureFromDensityInternalEnergy(rho[i], sie, lambda.data());
    double bmod = eos.BulkModulusFromDensityInternalEnergy(rho[i], sie, lambda.data());
    cs[i] = std::sqrt(robust::ratio(bmod, rho[i]));
  }

The exact same code is implemented via the python bindings in ``get_sound_speed_press.py``.

Get SESAME State
-------------------

If you build with both ``SpinerEOS`` and ``EOSPAC`` backends for
tabulated data, you can compare tabulated interpolations by calling
the ``get_sesame_state`` executable built via the
``get_sesame_state.cpp`` example file. You can call it as

.. code-block:: bash

   get_sesame_state matid sp5_file_name rho T sie

for some ``SESAME`` material index ``matid`` and some tabulated spiner
file ``sp5_file_name``, and a density, temperature and specific
internal energy to evaluate at.

The example demonstrates how to call the pressure, energy, and
thermodynamic derivatives of a table at that point in phase space.
