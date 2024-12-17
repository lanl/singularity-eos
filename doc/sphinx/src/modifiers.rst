.. _modifiers:

Equation of State Modifiers
============================

``EOS`` models can be *modified* by templated classes we call
*modifiers*. A modifier has exactly the same API as an ``EOS``, but
provides some internal transformation on inputs and outputs. For
example the ``ShiftedEOS`` modifier changes the zero point energy of a
given EOS model by shifting all energies up or down. Modifiers can be
used to, for example, production-harden a model.

Only certain combinations of ``EOS`` and ``modifier`` are permitted by
the defualt ``Variant``. For example, only ``IdealGas``,
``SpinerEOS``, and ``StellarCollapse`` support the ``RelativisticEOS``
and ``UnitSystem`` modifiers. All models support the ``ShiftedEOS``
and ``ScaledEOS`` modifiers. However, note that modifiers do not
commute, and only one order is supported. The ordering, inside-out, is
``UnitSystem`` or ``RelativisticEOS``, then ``ScaledEOS``, then
``ShiftedEOS``.

We list below the available modifiers and their constructors.

The Shifted EOS
-----------------

.. _modifiers shifted EOS:

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

for some scale parameter :math:`s`.  The ``ScaledEOS`` applies this
transformation to any equation of state, not just an ideal gas, where
the pressure may change for different scaling ratios.

Another way of understanding scaling ratios is that the pressure can be written
as

.. math::

  P = \left(\frac{\partial F}{\partial V} \right)_T

where :math:`F` is the Helmholtz free energy. For a given scaling such that
:math:`\rho_\mathrm{eos} = s\rho_\mathrm{in}`, the volume obeys the inverse
scaling. Since the scaling ratio is constant, it can be substituted into the
above expression so that

.. math::

  P = \left(\frac{\partial F_\mathrm{eos}}{\partial V_\mathrm{eos}} \right)_T
    = \left(\frac{\partial F_\mathrm{in}}{\partial V_\mathrm{in}} \right)_T
    = \left(\frac{\partial F_\mathrm{in}}{s \partial V_\mathrm{eos}} \right)_T
    = \left(\frac{s\partial F_\mathrm{eos}}{s \partial V_\mathrm{eos}} \right)_T

which implies that the Helmholtz free energy must scale in the same way as
volume (inverse to density) in order to preserve the same pressure. Applying
this scaling to the definition of the Helmholtz free energy yields

.. math::

  F_\mathrm{eos} = e_\mathrm{eos} - TS_\mathrm{eos} = \frac{1}{R} F_\mathrm{in}
    = \frac{1}{R}e_\mathrm{in} - T\left(\frac{1}{R}S_\mathrm{in}\right),

where the implicaiton is that this inverse the scaling ratio should also be
applied to energy. The inverse scaling ratio must be applied to the entropy
here in order to ensure that all other thermodynamic potentials
(energy, entropy, and the Gibbs free energy) scale similarly.

where :math:`e` is the internal energy and :math:`S` is the entropy. The
implication is that the same scaling should be applied to the energy and entropy
to maintain thermodynamic consistency.

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

for enthalpy by volume :math:`w`. The ``RelativisticEOS`` constructor accepts
the underlying model, and the speed of light as parameter. For example, a
relativistic ideal gas might be initialized as:

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = RelativisticEOS<IdealGas>(IdealGas(gm1, Cv), cl);

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

Z-Split EOS
-------------

For 3T physics (as described in the models section) it is often
desirable to have a separate equation of state for electrons and a
separate equation of state for ions. The Z-split model takes a total
equation of state and splits it into electron and ion components. The
"Z" here signifies mean ionization state, or the average number of
free electrons contributed per atomic nucleus, which is the mean
atomic number in the case of full ionization, but could be smaller in
the case of partial ionization. (It is zero for an unionized gas.) The
physical model of Z-split can be derived from an ideal gas equation of
state. For an ideal gas made up of electrons and ions, where all
molecular bonds have been broken, the total pressure is given by

.. math::

  P_t = (\left\langle Z\right\rangle + 1) \frac{\rho}{m_p \bar{A}} k_b T

where :math:`\left\langle Z\right\rangle` is the mean ionization
state, :math:`rho` is the ion mass density (the electron ion mass
density is negligible), :math:`m_p` is the proton mass,
:math:`\bar{A}` is the mean atomic mass, :math:`k_b` is Boltzmann's
constant, and :math:`T` is temperature. The contribution from
electrons is proportional to :math:`\left\langle Z\right\rangle`.

The split simply splits the total pressure into normalized
contributions for the electrons:

.. math::

  P_e = \frac{\left\langle Z\right\rangle}{1 + \left\langle Z\right\rangle} P_t(\rho, T_e)

and ions:

.. math::

  P_i = \frac{1}{1 + \left\langle Z\right\rangle} P_t(\rho, T_i)

where here :math:`T_e` is the electron temperature and :math:`T_i` the
ion temperature such that, when the temperatures are equal,

.. math::

  P_e + P_i = P_t

The same split is applied to the specific internal energy:

.. math::

  \varepsilon_e = \frac{\left\langle Z\right\rangle}{1 + \left\langle Z\right\rangle} \varepsilon_t(\rho, T_e)\\
  \varepsilon_i = \frac{1}{1 + \left\langle Z\right\rangle} \varepsilon_t(\rho, T_i)

and the remaining state variables and thermodynamic derivatives can be
derived from these relations.

In ``singularity-eos``, the z-split is implemented as the templated class

.. code-block:: cpp

  template<ZSplitComponent ztype, typename T>
  class ZSplit

where ``ZSPlitComponent`` may either be ``ZSplitComponent::Electrons``
or ``ZSplitComponent::Ions``. As syntactic sugar, ``ZSplitE<T>`` and
``ZSplitI<T>`` are available. The Z-split constructor does not require
any additional parameters, so you may construct one as, e.g.,

.. code-block:: cpp

  using namespace singularity
  auto ion_eos = ZSPlitI<IdealGas>(IdealGas(gm1, Cv);

and similarly for electrons,

.. code-block:: cpp

  auto electron_eos = ZSPlitE<IdealGas>(IdealGas(gm1, Cv);

The Z-split modifier takes the ionization state as an additional
parameter via the lambda. For example:

.. code-block:: cpp

  Real lambda[1] = {Z};
  Real Pe = electron_eos.PressureFromDensityTemperature(rho, temperature, lambda);

.. note::

  For now, the Z-split EOS is not in the default variant provided by
  singularity-eos. If you would like to use it, you must implement
  your own custom variant.

Composing Modifiers
--------------------

Modifiers can be composed. For example:

.. code-block:: cpp

  using namespace singularity;
  auto my_eos = ShiftedEOS<ScaledEOS<IdealGas>>(ScaledEOS(IdealGas(gm1, Cv), scale), shift);

You can build modifiers up iteratively by, for example:

.. code-block:: cpp

  using namespace singularity;
  EOS eos = IdealGas(gm1, cv);
  if (do_shift) {
    eos = eos.template Modify<ShiftedEOS>(shift);
  }
  if (do_scale) {
    eos = eos.template Modify<ScaledEOS>(scale);
  }

Undoing Modifiers
------------------

Modifiers can also be undone, extracting the underlying EOS. Continuing the example above,

.. code-block:: cpp

   auto unmodified = my_eos.GetUnmodifiedObject();

will extract the underlying ``IdealGas`` EOS model out from the scale and shift.
