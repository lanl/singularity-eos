.. _examples:

Examples
=========

The ``example`` directory of ``singularity-eos`` contains several
examples of using the code. You can build the examples by setting
``-DSINGULARITY_BUILD_EXAMPLES=ON`` at ``CMake`` configuration
time. For example:

.. code-block:: bash

   # from singularity-eos repo
   mkdir -p builddir && cd builddir
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

Custom EOS to SpinerEOS
-------------------------

The ``examples/custom_eos_to_spiner.cpp`` example demonstrates how to
tabulate a custom EOS implementation into ``SpinerEOS`` format using
the generic constructor. This is useful when you have your own EOS
physics model and want to improve performance through table
interpolation or enable GPU portability.

The example creates a simple Mie-Gruneisen-like EOS class that provides
the minimal required interface:

.. code-block:: cpp

  class CustomHostEOS {
  public:
    // Minimal required interface
    Real InternalEnergyFromDensityTemperature(Real rho, Real T) const;
    Real TemperatureFromDensityInternalEnergy(Real rho, Real sie) const;
    Real PressureFromDensityTemperature(Real rho, Real T) const;

    // Optional: improves derivative accuracy
    Real GruneisenParamFromDensityTemperature(Real rho, Real T) const;

    // Optional: material properties
    Real MeanAtomicMass() const;
    Real MeanAtomicNumber() const;
  };

It then demonstrates how to tabulate this custom EOS:

.. code-block:: cpp

  CustomHostEOS custom_eos(rho0, C0, s, Gamma0, Cv);

  // Set up grid parameters
  SpinerTableGridParams params;
  params.rhoMin = 1.0;
  params.rhoMax = 20.0;
  params.TMin = 300.0;
  params.TMax = 50000.0;
  params.numRhoPerDecade = 50;

  // Construct SpinerEOS from custom EOS
  SpinerEOSDependsRhoSie spiner_eos(custom_eos, params);

  // Use tabulated EOS
  Real P = spiner_eos.PressureFromDensityTemperature(rho, T);

The example verifies that the tabulated EOS matches the original custom
EOS with interpolation errors typically less than 1%. The resulting
``SpinerEOS`` can be used in production simulations for better
performance, on GPUs via ``GetOnDevice()``, and with mixed-cell closure
models.

This approach is particularly valuable for:

- Converting analytical models to tables for performance
- Integrating custom physics with codes that expect tabulated EOS
- Creating GPU-portable versions of host code EOS implementations
- Re-gridding existing tables to different resolutions

Plugins
----------

The example directory also contains an example plugin that can be
included via the plugin infrastructure, as described in :ref:`our customization section <customization>`.
