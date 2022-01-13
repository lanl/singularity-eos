.. _models:

EOS Models
===========

EOS Primer
----------

An equation of state (EOS) is a constituitive model that generally relates
thermodynamic quantities such as pressure, temperature, density, internal
energy, free energy, and entropy consistent with the constraints of equillibrium
thermodynamics.

``singularity-eos`` contains a number of equations of state that are most useful
in hydrodynamics codes. All of the equations of state presented here are
considered "complete" equations of state in the sense that they contain all the
information needed to derive a given thermodynamic quantity. However, not all
variables are exposed since the emphasis is on only those that are needed for
hydrodynamics codes. An incomplete equation of state is often sufficient for
pure hydrodynamics since it can relate pressure to a given energy-density state,
but it is missing information that relates the energy to temperature (which in
turn provides a means to calculate entropy).

The Mie-Gruneisen form
``````````````````````

Many of the following equations of state are considered to be of
the "Mie-Gruneisen form", which has many implications but for our purposes
means that the general form of the EOS is

.. math::

    P - P_\mathrm{ref} = \rho \Gamma(\rho) (E - E_\mathrm{ref})

where 'ref' denotes quantities along some reference curve, :math:`P` is the
pressure, :math:`\rho` is the density, :math:`\Gamma` is the Gruneisen
parameter, and :math:`E` is the specific internal energy. In this sense, an EOS
of this form uses the Gruneisen parameter to describe the pressure behavior of
the EOS away from the reference curve. Coupled with a relationship between
energy and temperature (sometimes as simple as a constant heat capacity), the
complete equation of state can be constructed.

To some degree it is the complexity of the reference state and the heat
capacity that will determine an EOS's ability to capture the complex behavior of
a material. At the simplest level, the ideal gas EOS uses a reference state at
zero pressure and energy, while more complex equations of state such as the
Davis EOS use the material's isentrope. In ths way, the reference curve
indicates the conditions under which you can expect the EOS to represent the
intended behavior.

EOS API
```````

The mathematical descriptions of these models are presented below while the
details of using them is presented in the description of the 
:ref:`EOS API <using-eos>`.


Implemented EOS models
----------------------


Ideal Gas
`````````

The ideal gas (aka perfect or gamma-law gas) model in ``singularity-eos`` takes
the form

.. math::

    P = \Gamma \rho \hat{e}

    \hat{e} = \hat{C_\mathrm{V}} T

where :math:`P`, :math:`\Gamma`, :math:`\rho`, :math:`e`, :math:`T`, and 
:math:`\hat{C_\mathrm{V}}` are the pressure, Gruneisen parameter, density,
specifc internal energy, temperature, and specific heat capacity at constant
volume respetively. In this context, "specific" implies per unit mass.

The settable parameters are the Gruneisen parameter and specific heat capacity.
Although this differs from the traditional representation of the ideal gas law
as :math:`P\hat{V} = RT`, the forms are equivalent by recognizing that
:math:`\Gamma = \frac{R}{\tilde{C_\mathrm{V}}}` where :math:`R` is the ideal gas
constant and :math:`\tilde{C_\mathrm{V}}` is the *molar* heat capacity,
relatable to the *specific* heat capacity through the molecular weight of the
gas. Since :math:`\tilde{C_\mathrm{V}} = \frac{5}{2} R` for a diatomic ideal
gas, the corresponding Gruneisen parameter should be 0.4.


Davis EOS
`````````

The Davis reactants and products EOS are both in Mie-Gruneisen forms that use
isentropes for the reference curves. The equations of state are typically used
to represent high explosives and their detonation products and the reference
curves are calibrated to several sets of experimental data.
