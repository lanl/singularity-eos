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
the default ``Variant``. All models support the ``ShiftedEOS``,
``ScaledEOS``, and ``BilinearRampEOS`` modifiers, which are the only
modifiers the default variant contains. Modifiers do not commute, and
the default variant provides exactly one nesting order. From the
innermost modifier outward that order is ``ShiftedEOS``, then
``ScaledEOS``, then ``BilinearRampEOS``; a shifted and scaled ideal gas
therefore has the concrete type ``ScaledEOS<ShiftedEOS<IdealGas>>``.

The remaining modifiers documented below --- ``RelativisticEOS``,
``UnitSystem``, ``ScaledFreeEnergy``, ``FlooredEnergy``, and
``ZSplit`` --- are *not* members of the default variant. You may still
use them by naming the concrete modified type directly, or by building
your own variant as described in :ref:`customization`. See
:ref:`modifier-ordering` below for how to choose the nesting order when
you do.

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

where here :math:`\Gamma` is the Gruneisen parameter, :math:`\rho` is
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

where the implication is that this inverse the scaling ratio should also be
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

The scale factor must be strictly positive. A negative scale is not
physically meaningful: it would invert the sign of energy and entropy,
and would turn every minimum bound reported by the bounds introspection
API into a maximum. ``CheckParams`` enforces this.


The Free-Energy Amplitude Scaling EOS
-------------------------------------

The ``ScaledEOS`` modifier changes the density at which an underlying
equation of state is evaluated. ``ScaledFreeEnergy`` performs a
different operation: it scales the Helmholtz free energy at fixed density and
temperature. If the specific Helmholtz free energy of the underlying model is
:math:`f(\rho, T)`, then the amplitude-scaled free energy is

.. math::

   f_\mathrm{scaled}(\rho, T) = a f(\rho, T),

where :math:`a` is the free-energy amplitude scale factor.

The thermodynamic quantities returned by the modifier follow from this free
energy. In particular,

.. math::

   \begin{aligned}
   P_\mathrm{scaled}(\rho, T) &= a P(\rho, T), \\
   e_\mathrm{scaled}(\rho, T) &= a e(\rho, T), \\
   S_\mathrm{scaled}(\rho, T) &= a S(\rho, T), \\
   C_{V,\mathrm{scaled}}(\rho, T) &= a C_V(\rho, T), \\
   K_\mathrm{scaled}(\rho, T) &= a K(\rho, T).
   \end{aligned}

The Gruneisen parameter is unchanged because pressure and specific internal
energy have the same amplitude factor:

.. math::

   \Gamma_\mathrm{scaled} = \frac{1}{\rho}
   \left(\frac{\partial P_\mathrm{scaled}}
   {\partial e_\mathrm{scaled}}\right)_\rho = \Gamma.

For calls with density and specific internal energy as inputs, the underlying
model is evaluated at the inverse-scaled energy:

.. math::

   \begin{aligned}
   T_\mathrm{scaled}(\rho, e) &= T\left(\rho, e/a\right), \\
   P_\mathrm{scaled}(\rho, e) &= a P\left(\rho, e/a\right).
   \end{aligned}

Thus, unlike ``ScaledEOS``, free-energy amplitude scaling does not change the
density coordinate of the underlying EOS. For an ideal gas, ``ScaledEOS``
leaves pressure unchanged, while ``ScaledFreeEnergy`` multiplies it
by :math:`a`.

The constructor for ``ScaledFreeEnergy`` takes the underlying model
and the amplitude. For example, an amplitude-scaled ideal gas might be
initialized as:

.. code-block:: cpp

   using namespace singularity;
   auto my_eos = ScaledFreeEnergy<IdealGas>(IdealGas(gm1, Cv), amplitude);

The amplitude must be strictly positive. A negative amplitude would reverse
the signs of the heat capacity and bulk modulus, violating thermodynamic
stability. Positivity also preserves the ordering of the minimum and maximum
specific-internal-energy bounds, both of which are multiplied by the
amplitude.

.. note::

  ``ScaledFreeEnergy`` is not in the default variant provided
  by ``singularity-eos``. The example above uses ``auto`` for that
  reason: assigning the modified object to the default ``EOS`` variant
  will abort at runtime with ``Modifier must be in variant``. To hold an
  amplitude-scaled model in a variant, you must build your own; see
  :ref:`custom-variant-modifiers`.

.. _scaled-free-energy-ordering:

Where to Place ``ScaledFreeEnergy``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ScaledFreeEnergy`` is multiplicative, so by the rules in
:ref:`modifier-ordering` it commutes exactly with ``ScaledEOS``:
``ScaledFreeEnergy<ScaledEOS<T>>`` and
``ScaledEOS<ScaledFreeEnergy<T>>`` agree on every method. The two are
not the same transformation, though. ``ScaledEOS`` is a mass rescaling,
so it applies no amplitude factor to pressure or bulk modulus, although
it does rescale the density argument. ``ScaledFreeEnergy`` multiplies
both by :math:`a` and leaves the density argument alone.

Its position relative to ``ShiftedEOS`` does matter. For an underlying
energy :math:`e_b`, a shift :math:`q`, a scale :math:`s`, and an
amplitude :math:`a`, the two candidate placements give

.. math::

  \begin{aligned}
  \mathtt{ScaledFreeEnergy}
    &\left[\mathtt{ScaledEOS}\left[\mathtt{ShiftedEOS}\right]\right]
    : \quad e(\rho, T) = a s\, e_b(s \rho, T) + a s q, \\
  \mathtt{ScaledEOS}
    &\left[\mathtt{ShiftedEOS}\left[\mathtt{ScaledFreeEnergy}\right]\right]
    : \quad e(\rho, T) = a s\, e_b(s \rho, T) + s q.
  \end{aligned}

The two differ only in whether the energy offset :math:`q` picks up the
amplitude. Placing ``ScaledFreeEnergy`` outside the shift, as in the
first form, treats the offset the same way ``ScaledEOS`` already does in
the default variant, and is the recommended choice. Placing it
innermost, directly on the base model, leaves a host-supplied reference
energy offset untouched; prefer that form if the amplitude belongs to
the material's free-energy fit rather than to the assembled EOS.

A custom variant holding an amplitude-scaled, shifted, and scaled ideal
gas in the recommended order is then

.. code-block:: cpp

  #include <singularity-eos/eos/eos.hpp>
  #include <singularity-eos/eos/modifiers/scaled_free_energy.hpp>
  using namespace singularity;

  using MyEOS_t =
      Variant<IdealGas, ShiftedEOS<IdealGas>, ScaledEOS<ShiftedEOS<IdealGas>>,
              ScaledFreeEnergy<ScaledEOS<ShiftedEOS<IdealGas>>>>;

Remember that the intermediate types are required as well as the fully
modified one; see :ref:`custom-variant-modifiers`.


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

All of the bounds introspection methods described in the :ref:`EOS API
section<using-eos>`, including ``MinimumInternalEnergy`` and
``MaximumInternalEnergy``, are converted to the new unit system. Because
they are part of the general EOS API, they are also reachable through
the ``singularity::EOS`` variant, so one no longer needs to call
``GetUnmodifiedObject`` to retrieve energy bounds, which would return
them in the unmodified (cgs) unit system:

.. code-block:: cpp

  using namespace singularity;
  EOS my_eos = UnitSystem<SpinerEOSDependsRhoSie>(/* ... */);
  // bounds in the new unit system
  Real sie_min = my_eos.MinimumInternalEnergy();
  Real sie_max = my_eos.MaximumInternalEnergy();

See :ref:`How Modifiers Transform the Energy
Bounds<modifier-energy-bounds>` below.

.. note::

  ``UnitSystem`` does *not* forward the model-specific table bounds
  accessors ``sieMin()``/``sieMax()``. Those are not part of the general
  EOS API, so they cannot be converted for an arbitrary underlying
  model. Use ``MinimumInternalEnergy``/``MaximumInternalEnergy``
  instead, which every model provides.

Z-Split EOS
-------------

For 3T physics (as described in the :ref:`models section<3T-Model>`) it is often
desirable to have a separate equation of state for electrons and a
separate equation of state for ions. The Z-split model takes a total
equation of state and splits it into electron and ion components. The
"Z" here signifies mean ionization state, or the average number of
free electrons contributed per atomic nucleus, which is the mean
atomic number in the case of full ionization, but could be smaller in
the case of partial ionization. (It is zero for an unionized gas.) The
physical model of Z-split can be derived by approximating the
material as an ideal gas. For an ideal gas made up of electrons and ions, where all
molecular bonds have been broken, the total pressure is given by

.. math::

  P_t = (\left\langle Z\right\rangle + 1) \frac{\rho}{m_p \overline{A}} k_b T

where :math:`\left\langle Z\right\rangle` is the mean ionization
state, :math:`\rho` is the ion mass density (the electron ion mass
density is negligible), :math:`m_p` is the proton mass,
:math:`\overline{A}` is the mean atomic mass, :math:`k_b` is Boltzmann's
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

.. note::

  The Z-split modifier uses the mean atomic properties methods
  provided by the underlying equation of state to pick, e.g., the
  total number of nuclei per unit mass. This means you must specify
  ``MeanAtomicProperties`` when constructing the underlying equation
  of state. If you don't specify anything, hydrogen is assumed.

The Z-split modifier takes the ionization state as an additional
parameter via the lambda. For example:

.. code-block:: cpp

  Real lambda[1] = {Z};
  Real Pe = electron_eos.PressureFromDensityTemperature(rho, temperature, lambda);

If you use the indexable types infrastructure, the name is
``IndexableTypes::MeanIonizationState``.

.. warning::

  Several thermodynamic properties are approximated in the z-split
  model due to incomplete information. Notably, the specific heat and
  bulk modulus should depend on the electron chemical potential and
  the ionization model: how much does the ionization state change with
  respect to temperature? However, the ionization model is typically
  decoupled from the equation of state and treated separately. As
  such, this dependence is neglected here. This treatment may be
  extended in the future.

.. note::

  For now, the Z-split EOS is not in the default variant provided by
  singularity-eos. If you would like to use it, you must implement
  your own custom variant.

Energy Floor Modifier
---------------------

The models in ``singularity-eos`` do not guarantee that a lookup will lie on the
EOS surface. The laws of thermodynamics require that volume (or density) and
temperature be positive, but there are no natural lower bounds on other
variables, particularly energy.

In general, this modifer imposes a lower bound on the energy by finding the
energy along the zero kelvin isotherm, i.e. :math:`e(\rho, T=0)`, and then
flooring the energy to that value when performing a lookup. As a result, any
lookups in energy space that would otherwise be below the zero K isotherm return
whatever value would lie upon the isotherm.

Notably though, tabular SESAME tables can result in unexpected behavior. The
``EOSPAC`` and ``SpinerEOSDependsRhoT`` models both rely on the result when
EOSPAC querries the cold energy table (306). This energy *can* differ from the
zero K isotherm for two reasons:

#. The zero K isotherm includes the zero point energy. In other words, it is the
   energy of the lowest energy level whereas the cold curve is the energy of the
   potential well itself. This effect may be insignificant.

#. In the absence of a 306 table for the material, EOSPAC will instead report
   the lowest isotherm's energy in place of a cold curve. If a given table does
   not extend to zero K, then an extrapolaiton is *not* performed. This results
   in an energy floor from the lowest temperature on the table, **not** from
   an extrapolated energy at zero K.

.. note::

  This modifier is only compatible with EOS where the
  ``MinInternalEnergyFromDensity`` function is properly defined. Where it isn't
  defined, an error will occur at runtime.


.. _composing-modifiers:

Composing Modifiers
--------------------

Modifiers can be composed. For example:

.. code-block:: cpp

  using namespace singularity;
  auto my_eos = ScaledEOS<ShiftedEOS<IdealGas>>(ShiftedEOS<IdealGas>(IdealGas(gm1, Cv), shift), scale);

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

.. _modifier-ordering:

Why the Order Matters
^^^^^^^^^^^^^^^^^^^^^^

A modifier wraps the model below it: in ``Mod<T>``, ``Mod`` transforms
the inputs on their way down into ``T`` and the outputs on their way
back up. The *outermost* modifier therefore acts first on every input
and last on every output. Because of this, nesting order is part of the
physics, not a matter of taste.

Whether two modifiers commute depends on the kind of transformation each
applies:

* *Multiplicative* modifiers, such as ``ScaledEOS`` and
  ``ScaledFreeEnergy``, multiply inputs and outputs by constant
  factors. Constant factors commute, so multiplicative modifiers commute
  with one another exactly, and the order you choose among them is free.
  Note that commuting does not make them equivalent: they generally
  scale different quantities.
* ``ShiftedEOS`` is *additive*, so it commutes with neither itself nor a
  multiplicative modifier. A multiplicative modifier placed *outside* a
  shift rescales the shift along with everything else, while one placed
  *inside* leaves the shift at its original magnitude. This is the
  asymmetry to watch: a reference-energy offset does not survive an
  outer rescaling unchanged.
* ``BilinearRampEOS`` should always be outermost, which is where the
  default variant puts it.

The default variant resolves these choices by nesting ``ShiftedEOS``
innermost, then ``ScaledEOS``, then ``BilinearRampEOS``. The shift is
therefore rescaled by the scale factor. When you place an additional
multiplicative modifier in your own variant, deciding whether it sits
inside or outside ``ShiftedEOS`` is the substantive decision; see
:ref:`scaled-free-energy-ordering` for a worked example.

.. _custom-variant-modifiers:

Building a Custom Variant
^^^^^^^^^^^^^^^^^^^^^^^^^^

Any modifier that is not in the default variant must be held in a
variant you declare yourself. For example, to hold a shifted and scaled
ideal gas:

.. code-block:: cpp

  #include <singularity-eos/eos/eos.hpp>
  using namespace singularity;

  using MyEOS_t =
      Variant<IdealGas, ShiftedEOS<IdealGas>, ScaledEOS<ShiftedEOS<IdealGas>>>;

Use ``singularity::Variant``, not ``singularity::eos_variant``. The
latter is a bare alias for ``PortsOfCall::variant`` and carries none of
the EOS API; ``Variant`` is the wrapper that holds it and provides
``Modify``, the lookup methods, and the vector overloads.

Every *intermediate* type must be listed, not only the fully modified
one. ``Modify`` builds the stack one layer at a time, and each
intermediate concrete type has to be a member of the variant:

.. code-block:: cpp

  MyEOS_t eos = IdealGas(gm1, Cv);
  eos = eos.Modify<ShiftedEOS>(shift);   // needs ShiftedEOS<IdealGas>
  eos = eos.Modify<ScaledEOS>(scale);    // needs ScaledEOS<ShiftedEOS<IdealGas>>

Omitting an intermediate type is not a compile-time error. It aborts at
runtime with ``Modifier must be in variant``, because the variant must
be able to hold the result of every ``Modify`` call it might see.

For larger variants, the ``variadic_utils`` helpers that
``singularity-eos`` uses internally apply a modifier across a whole
list of models at once. The pattern, which mirrors
``singularity-eos/eos/eos_type_lists.hpp``, is

.. code-block:: cpp

  #include <singularity-eos/eos/eos_type_lists.hpp>
  using namespace singularity;

  static constexpr const auto my_models = tl<IdealGas, Gruneisen>{};
  static constexpr const auto shifted =
      transform_variadic_list(my_models, al<ShiftedEOS>{});
  static constexpr const auto scaled =
      transform_variadic_list(shifted, al<ScaledEOS>{});
  static constexpr const auto my_list =
      variadic_utils::concat(my_models, shifted, scaled);
  using MyEOS_t = typename decltype(tl_to_Variant(my_list))::vt;

Each call to ``transform_variadic_list`` wraps every model in the list
it is given, so modifiers are layered by chaining the calls, exactly as
in the hand-written form above.

Note that applying a modifier across a list multiplies the number of
types in the variant, and every alternative is instantiated for every
call site. Adding one modifier over the full default list would roughly
double its size, so list only the combinations your code actually uses.

Undoing Modifiers
------------------

Modifiers can also be undone, extracting the underlying EOS. Continuing
the ``my_eos`` example from :ref:`Composing Modifiers<composing-modifiers>`,

.. code-block:: cpp

   auto unmodified = my_eos.GetUnmodifiedObject();

will extract the underlying ``IdealGas`` EOS model out from the scale and shift.

.. _modifier-energy-bounds:

How Modifiers Transform the Energy Bounds
------------------------------------------

Specific internal energy is transformed by several modifiers, so the
energy bounds reported by ``MinimumInternalEnergy`` and
``MaximumInternalEnergy`` (described in the :ref:`EOS API
section<using-eos>`) are transformed to match:

* ``UnitSystem`` divides both bounds by the energy unit.
* ``ShiftedEOS`` adds the shift to both bounds.
* ``ScaledEOS`` multiplies both bounds by the scale factor. Since the
  scale factor is required to be positive, this preserves their ordering.
* ``ScaledFreeEnergy`` multiplies both bounds by the amplitude factor.
* ``RelativisticEOS`` and ``BilinearRampEOS`` leave energy alone, so the
  bounds pass through unchanged.
* ``FlooredEnergy`` clamps energy to the per-density cold curve, which
  lies inside the underlying energy range, so the global bounds pass
  through unchanged.
* ``ZSplit`` scales energy by a factor that depends on the ionization
  state passed through the ``lambda``, but the bounds introspection API
  takes no ``lambda``. The bounds it reports are therefore the
  *un-split* bounds of the underlying EOS.

This file was made in part with generative AI.
