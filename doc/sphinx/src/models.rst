.. _models:

EOS Models
===========

An equation of state (EOS) is a constituitive model that generally relates
thermodynamic quantities such as pressure, temperature, density, internal
energy, free energy, and entropy consistent with the constraints of equillibrium
thermodynamics.

``singularity-eos`` contains a number of equations of state that are most useful
in hydrodynamics codes. All of the equations of state here are considered
"complete" equations of state in the sense that they contain all the information
needed to derive a given thermodynamic quantity. However, not all variables are
exposed since the emphasis is on only those that are needed for hydrodynamics
codes.

The mathematical descriptions of these models are presented below while the
details of using them is presented in the description of the 
:ref:`EOS API <using-eos>`.


Ideal Gas
---------

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
---------

The Davis reactants and products EOS are more generally of a Mie-Gruneisen form
that uses an isentrope for the reference curve and a non-constant heat capacity.
