.. _modifiers:

Equation of State Modifiers
============================

``EOS`` models can be *modified* by templated classes we call
*modifiers*. A modifier has exactly the same API as an ``EOS``, but
provides some internal transformation on inputs and outputs. For
example the ``ShiftedEOS`` modifier changes the zero point energy of a
given EOS model by shifting all energies up or down. Modifiers can be
used to, for example, production-harden a model. Only certain
combinations of ``EOS`` and ``modifier`` are permitted by the defualt
``Variant``. For example, only ``IdealGas``, ``SpinerEOS``, and
``StellarCollapse`` support the ``RelativisticEOS`` and ``UnitSystem``
modifiers. All models support the ``ShiftedEOS`` and ``ScaledEOS``
modifiers. However, note that modifiers do not commute, and only one
order is supported. The ordering, inside-out, is ``UnitSystem`` or
``RelativisticEOS``, then ``ScaledEOS``, then ``ShiftedEOS``.

We list below the available modifiers and their constructors.

The Shifted EOS
-----------------

The shifted equation of state modifies zero point energy of an
underlying model by some shift. So for example, it transforms

.. math::

  P(\rho, \varepsilon) \to P(\rho, \varepsilon - \varepsilon_0)

for some shift :math:`\varepsilon_0`. This is a permitted,
thermodynamically consistent operation, the energy that corresponds to
"zero" is a free gauge parameter.

The constructor for the ``ShiftedEOS`` takes the underlying model and
the shift parameter. For example, a shifted ideal gas might be
initialized as:

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = ShiftedEOS<IdealGas>(IdealGas(gm1, Cv), shift);

where the first two parameters are the Gruneisen parameter and
specific heat required by the ideal gas constructor and the latter is
the energy shift.

The Scaled EOS
---------------

To understand the scaled EOS, consider the pressure for an ideal gas:

.. math::

  P = \Gamma \rho \varepsilon

where here :math:`\Gamma` is the Gruneien parameter, :math:`\rho` is
the density, and :math:`\varepsilon` is the specific internal
energy. The pressure is unchanged under the operation

.. math::

  \rho \to s\rho,\ \varepsilon\to \varepsilon/s

for some scale parameter :math:`s`.  The ``ScaledOES`` applies this
transformation to any equation of state, not just an ideal gas, where
the pressure may change.

The constructor for ``ScaledEOS`` accepts the underlying model, and
the scale parameter. For example, a shifted ideal gas might be
initialized as:

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = ScaledEOS<IdealGas>(IdealGas(gm1, Cv), scale);

where the first two parameters are the Gruneisen parameter and
specific heat required by the ideal gas constructor and the latter is
the scale.

The Relativistic EOS
---------------------

The relativistic modifier modifies the bulk modulus to enforce that
the sound speed, defined as

.. math::

  c_s = \sqrt{B_S/\rho}

is always less than the speed of light. It does so by applying the
transformation

.. math::

  B_S \to B_S/h

for the specific enthalpy :math:`h`. This brings the sound speed formula into alignment with the relativistic version,

.. math::

  c_s = \sqrt{B_S/w}

for enthalpy by volume :math:`w`. The ``RelativisticEOS`` constructor
doesn't take any additional parameters, simply construct one with the
parameters for the underlying EOS:

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = RelativisticEOS<IdealGas>(IdealGas(gm1, Cv));

EOS Unit System
-----------------

By default, the ``singularity-eos`` models all use cgs units. However,
it is often desirable to modify the units used to interact with the
library. The ``UnitSystem`` modifier partially implements this
functionality.

In particular, when constructing an EOS modified by the
``UnitSystem``, the user may specify a new unit system either by
thermal units, specific internal energy, and temperature, or by
length, mass, and time units. Then all calls of the modified EOS will
expect values in the new units and return values in the new units.

The way units are specified is via tag dispatch. For example

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = UnitSystem<IdealGas>(IdealGas(gm1, Cv),
    eos_units_init::ThermalUnitsInit(),
    rho_unit, sie_unit, temp_unit);

specifies the unit system by specifying units for density, specific
internal energy, and temperature. On the other hand,

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = UnitSystem<IdealGas>(IdealGas(gm1, Cv),
    eos_units_init::LengthTimeUnitsInit(),
    time_unit, mass_unit, length_unit, temp_unit);

specifies the unit system by specifying units for time, mass, length,
and temperature.

Composing Modifiers
--------------------

Modifiers can be composed. For example:

.. code-block:: cpp

  using namespace singularity;
  auto my_eos = ShiftedEOS<ScaledEOS<IdealGas>(ScaledEOS(IdealGas(gm1, Cv), scale), shift);
