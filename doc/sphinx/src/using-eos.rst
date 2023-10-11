.. _using-eos:

The Equation of State API
=========================

This page describes the equation of state API in detail. For just the
information needed to get started, check out the :ref:`getting started
<getting-started>` page.

The ``Real`` Type
------------------

``singularity-eos`` defines the ``singularity::Real`` type as a proxy
for the ``float`` and ``double`` types. We currently resolve ``Real``
to a double precision number, however we plan to have the option to
select different precisions at compile time in the future.

The parallelism model
----------------------

For the most part, ``singularity-eos`` tries to be agnostic to how you
parallelize your code on-node. (It knows nothing at all about
distributed memory parallelism.) An ``EOS`` object can be copied into
any parallel code block by value (see below) and scalar calls do not
attempt any internal multi-threading, meaning ``EOS`` objects are not
thread-safe, but are compatible with thread safety, assuming the user
calls them appropriately. The main complication is ``lambda`` arrays,
which are discussed below.

The vector ``EOS`` method overloads are a bit different. These are
thread-parallel operations launched by ``singularity-EOS``. They run
in parallel, and ordering between indices of vectors cannot be
assumed. Therefore care must taken in memory layout to avoid race
conditions. The type of parallelism used depends on how
``singularity-eos`` is compiled. If the ``Kokkos`` backend is used,
any parallel dispatch supported by ``Kokkos`` is supported.

.. _variant section:

Variants
---------

The equation of state library is object oriented, and uses a kind of
type erasure called a `Variant`_. (Technically we use a backport of
this C++ feture to C++11, see: `mpark variant`_.) The salient detail
is that a variant is a kind of compile-time polymorphism.

.. _Variant: https://en.cppreference.com/w/cpp/utility/variant

.. _mpark variant: https://en.cppreference.com/w/cpp/utility/variant

The ``singularity::EOS`` class is generic and can be initialized as
any equation of state model listed in :ref:`the models section
<models>`. Unlike with standard polymorphism, you don't need to
initialize your equation of state as a pointer. Rather, just use the
assignment operator. For example:

.. code-block:: cpp

   singularity::EOS my_eos = singularity::IdealGas(gm1, Cv);

To make this machinery work, there's an underlying variatic class,
``singularity::Variant``, defined in
``singularity-eos/eos/eos_variant.hpp``. Only methods defined for the
``singularity::Variant`` class are available for the equation of state
models. Moreover, any new equation of state model must define all
methods defined in the ``singularity::Variant`` class that call the ``visit``
function, or compile errors may occur.

If you wish to extract an underlying EOS model as an independent type,
undoing the type erasure, you can do so with the ``get``
method. ``get`` is templated and type deduction is not possible. You
must specify the type of the class you're pulling out of the
variant. For example:

.. code-block::

   auto my_ideal_gas = my_eos.get<singularity::IdealGas>();

This will give you access to methods and fields which may be unique to
a class but not shared by the ``Variant``.

The EOS model also allows some host-side introspection. The method

.. cpp:function:: static std::string EosType();

returns a string representing the equation of state an ``EOS`` object
currently is. For example:

.. code-block::

  auto tpe_str = my_ideal_gas.EosType();
  // prints "IdealGas"
  std::cout << tpe_str << std::endl;

Similarly the method

.. cpp:function:: void PrintParams() const;

prints relevant parameters that the EOS object was created with, such
as the Gruneisen coefficient and specific heat for an ideal gas model.

If you would like to create your own custom variant with additional
models (or a subset of models), you may do so by using the
``eos_variant`` class. For example,

.. code-block:: cpp

  #include <singularity-eos/eos.hpp>
  using namespace singularity;
  
  using MyEOS_t = eos_variant<IdealGas, Gruneisen>;

This will create a new type, ``MyEOS_t`` which contains only the
``IdealGas`` and ``Gruneisen`` classes. (All of these live under the
``singularity`` namespace.)

Reference Semantics and ``GetOnDevice``
-----------------------------------------

Equation of state objects in ``singularity-eos`` have so-called
*reference-semantics*. This means that when a variable is copied or
assigned, the copy is *shallow*, and underlying data is not moved,
only metadata. For analytic models this is essentially irrelevant, the
only data they contain is metadata, which is copied. For tabulated
models such as ``SpinerEOS``, this matters more.

In a heterogenous environment, e.g., where both a CPU and an GPU are
available, data is allocated on the host by default. It can be copied
to device via

.. cpp:function:: void EOS::GetOnDevice()

which can be called as, e.g.,

.. code-block:: cpp

  eos.GetOnDevice();

Once data is on device, ``EOS`` objects can be trivially copied into
device kernels by value. The copy will be shallow, but the data will
be available on device. In Cuda, this may mean passing the EOS in as a
function parameter into a kernel. In a higher-level abstraction like
Kokkos, simply capture the object into a device lambda by value.

Underlying data is **not** reference-counted, and must be freed by
hand. This can be achieved via the

.. cpp:function:: void EOS::Finalize()
   
method, which can be called as, e.g.,

.. code-block:: cpp

  eos.Finalize();

Vector and Scalar API, Accessors
---------------------------------

Most ``EOS`` methods have both scalar and vector overloads, where the
scalar version returns a value, and the vector version modifies an
array. By default the vector version is called from host on device (if
``singularity-eos`` was compiled for device).

The vector API is templated to accept *accessors*. An accessor is any
object with a square bracket operator. One-dimensional arrays,
pointers, and ``std::vector<double>`` are all examples of what we call
an accessor. However, the value of an accessor is it doesn't have to
be an array. You can create an accessor class that wraps your
preferred memory layout, and ``singularity-eos`` will handle it
appropriately. An accessor that indexes into an array with some stride
might look like this:

.. code-block:: cpp

  struct Indexer {
    Indexer(int stride, double *array) : stride_(stride), A_(array) {}
    double &operator[](int i) {
      return A_[stride*i];
    }
    double *A_;
    int stride_;
  };

We do note, however, that vectorization may suffer if your underlying
data structure is not contiguous in memory.

.. _eospac_vector:

EOSPAC Vector Functions
~~~~~~~~~~~~~~~~~~~~~~~

For performance reasons EOSPAC vector calls only support contiguous memory
buffers as input and output. They also require an additional scratch buffer.

These changes are needed to allow passing buffers directly into EOSPAC, taking
advantage of EOSPAC options, and avoiding unnecessary copies.

The size of the needed scratch buffer depends on which EOS function is called
and the number of elements in the vector. Use the ``scratch_size(func_name, num_elements)``
static member function to determine the size needed for an individual function
or ``max_scratch_size(num_elements)`` to retrieve the maximum needed by any
available member function.

.. code-block:: cpp

   // std::vector<double> density = ...;
   // std::vector<double> energy = ...;
   // std::vector<double> temperature = ...;

   // determine size and allocate needed scratch buffer
   auto sz = EOSPAC::scratch_size("TemperatureFromDensityInternalEnergy", density.size());
   std::vector<double> scratch(sz / sizeof(double));

   // call EOSPAC eos vector function with scratch buffer
   eos.TemperatureFromDensityInternalEnergy(density.data(), energy.data(), temperature.data(),
                                            scratch.data(), density.size());

The Evaluate Method
~~~~~~~~~~~~~~~~~~~

A special call related to the vector calls is the ``Evaluate``
method. The ``Evaluate`` method requests the EOS object to evaluate
almost arbitrary code, but in a way where the type of the underlying
EOS object is resolved *before* this arbitrary code is evaluated. This
means the code required to resolve the type of the variant is only
executed *once* per ``Evaluate`` call. This can enable composite EOS
calls, non-standard vector calls, and vector calls with non-standard
loop structure.

The ``Evaluate`` call has the signature

.. code-block:: cpp

  template<typename Functor_t>
  PORTABLE_INLINE_FUNCTION
  void Evaluate(Functor_t f);

where a ``Functor_t`` is a class that *must* provide a ``void
operator() const`` method templated on EOS type. ``Evaluate`` is
decorated so that it may be evaluated on either host or device,
depending on desired use-case. Alternatively, you may use an anonymous
function with an `auto` argument as the input, e.g.,

.. code-block::

   eos.Evaluate([=](auto eos) { /* my code snippet */ });

.. warning::

  It can be dangerous to use functors with side-effects. Especially
  with GPUs it can produce very unintuitive behaviour. We recommend
  you only make the ``operator()`` non-const if you really know what
  you're doing. And in the anonymous function case, we recommend you
  capture by value, not reference.

To see the utlity of the ``Evaluate`` function, it's probably just
easiest to provide an example. The following code evaluates the EOS on
device and compares against a tabulated pressure. The total difference
is summed using the ``Kokkos::parallel_reduce`` functionality in the
``Kokkos`` performance portability library.

.. code-block:: cpp

  // The functor we use is defined here.
  // This class definition needs to be of appropriately global scope.
  class CheckPofRE {
   public:
    CheckPofRE(Real *P, Real *rho, Real *sie, int N) : P_(P), rho_(rho), sie_(sie), N_(N) {}
    template <typename T>
    // this is a host-only call, but if you wanted to write
    // a function that you wanted to evaluate on device
    // you could add the
    // PORTABLE_INLINE_FUNCTION
    // decorator here.
    void operator()(const T &eos) const {
      // Capturing member functions of a class in a lambda typically causes problems
      // when launching a GPU kernel.
      // Better to pull out new variables to capture before launching a kernel.
      Real *P = P_;
      Real *rho = rho_;
      Real *sie = sie_;
      // reduction target
      Real tot_diff;
      // reduction op
      Kokkos::parallel_reduce(
          "MyCheckPofRE", N_,
          KOKKOS_LAMBDA(const int i, Real &diff) {
            diff += std::abs(P[i] - eos.PressureFromDensityInternalEnergy(rho[i], sie[i]));
          },
          tot_diff);
      std::cout << "Total difference = " << tot_diff << std::endl;
    }
  
   private:
    int N_;
    Real *P_;
    Real *rho_;
    Real *sie_;
  };

  // Here we construct our functor
  // it is assumed the pointers to device memory P, rho, sie, are defined elsewhere.
  CheckPofRE my_op(P, rho, sie, N);

  // Here we call the evaluate function
  eos.Evaluate(my_op);

  // The above two lines could have been called "in-one" with:
  // eos.Evaluate(CheckPofRE(P, rho, sie, N));

Alternatively, you could eliminate the functor and use an anonymous
function with:

.. code-block:: cpp

  eos.Evaluate([=](auto eos) {
    Real tot_diff;
    Kokkos::parallel_reduce(
        "MyCheckPofRE", N_,
        KOKKOS_LAMBDA(const int i, Real &diff) {
          diff += std::abs(P[i] - eos.PressureFromDensityInternalEnergy(rho[i], sie[i]));
        },
        tot_diff);
    std::cout << "Total difference = " << tot_diff << std::endl;
  });

This is not functionality that would be available with the standard
vector calls provided by ``singularity-eos``, at least not without
chaining multiple parallel dispatch calls. Here we can do it in a
single call.

Lambdas and Optional Parameters
--------------------------------

Most methods for ``EOS`` objects accept an optional ``lambda``
parameter, which is a ``Real *``. Unless specified in :ref:`the
models section <models>`, this parameter does nothing. However, some
models require or benefit from additional information. For example
models with internal root finds can leverage initial guesses and
models with composition mixing parameters may need additional input to
return a meaningful state.

``EOS`` models are introspective and can provide the desired/required
size of the lambda array with:

.. cpp:function:: int EOS::nlambda()

which is the desired size of the ``lambda`` array per scalar call. For
vector calls, there should be one such array per grid point. An
accessor for ``lambda`` should return a ``Real *`` pointer at each
index. A trivial example of such an indexer for ``lambda`` might be
the null indexer:

.. code-block:: cpp

  class NullIndexer {
    Real *operator[](int i) { return nullptr; }
  };

As a general rule, to avoid race conditions, you will want at least
one ``lambda`` array (or subview of a larger memory allocation) per
thread. You may want one array per point you are evaluating
on. Ideally these arrays are persistent between ``EOS`` calls, to
minimize latency due to ``malloc`` and ``free``. Several models, such
as ``SpinerEOS`` also use the persistency of these arrays to cache
useful quantities for a performance boost.

EOS Modifiers
--------------

``EOS`` models can be *modified* by templated classes we call
*modifiers*. A modifier has exactly the same API as an ``EOS``, but
provides some internal transformation on inputs and outputs. For
example the ``ShiftedEOS`` modifier changes the reference energy of a
given EOS model by shifting all energies up or down. Modifiers can be
used to, for example, production-harden a model. Only certain
combinations of ``EOS`` and ``modifier`` are permitted by the defualt
``Variant``. For example, only ``IdealGas``, ``SpinerEOS``, and
``StellarCollapse`` support the ``RelativisticEOS`` and ``UnitSystem``
modifiers. All models support the ``ShiftedEOS`` and ``ScaledEOS``
modifiers. However, note that modifiers do not commute, and only one
order is supported. The ordering, inside-out, is ``UnitSystem`` or
``RelativisticEOS``, then ``ScaledEOS``, then ``ShiftedEOS``.

Relevant to the broad ``singularity-eos`` API, EOS models provide
introspection. To check if an EOS is modified, call

.. cpp:function:: bool IsModified() const;

This will return ``true`` for a modified model and ``false``
otherwise. Modifiers can also be undone. To get a completely
unmodified EOS model, call

.. cpp:function:: auto GetUnmodifiedObject();

The return value here will be either the type of the ``EOS`` variant
type or the unmodified model (for example ``IdealGas``) or, depending
on whether this method was callled within a variant or on a standalone
model outside a variant.

If you have chained modifiers, e.g.,
``ShifedEOS<ScaledEOS<IdealGas>``, you can undo only one of the
modifiers with the

.. cpp:function:: auto UnmodifyOnce();

method, which has the same return type pattern as above, but only
undoes one level of modification.

For more details on modifiers, see the :ref:`modifiers<modifiers>`
section. If you need a combination of modifiers not supported by
default, we recommend building a custom variant as described above.

Preferred Inputs
-----------------

Some equations of state, such as those built on tabulated data, are
most performant when quantities, e.g., pressure, are requested in
terms of density and temperature. Others may be most performant for
density and specific internal energy.

Most fluid codes work in terms of density and energy. However, for a
model that prefers density and temperature inputs, it may be better
compute temperature first, then compute other quantities given density
and temperature, rather than computing everything from density and
energy.

``singularity-eos`` offers some introspection to enable users to
determine what the right sequence of calls to make is:

.. cpp:function:: static constexpr unsigned long PreferredInput();

The return value is a bit field, represented as a number, where each
nonzero bit in the field represents some thermodynamic quantity like
density or temperature. You can check whether or not an eos prefers
energy or temperature as an input via code like this:

.. code-block:: cpp

  using namespace singularity;
  auto preferred_input = my_eos.PreferredInput();
  bool en_preferred = preferred_input & thermalqs::specific_internal_energy;
  bool temp_preferred = preferred_input & thermalqs::temperature;

Here the bitwise and operator masks out a specific flag, allowing one
to check whether or not the bitfield contains that flag.

The available flags in the ``singulartiy::thermalqs`` namespace are
currently:
* ``thermalqs::none``
* ``thermalqs::density``
* ``thermalqs::specific_internal_energy``
* ``thermalqs::pressure``
* ``thermalqs::temperature``
* ``thermalqs::specific_heat``
* ``thermalqs::bulk_modulus``
* ``thermalqs::all_values``

however, most EOS models only specify that they prefer density and
temperature or density and specific internal energy.

.. _eos builder section:

EOS Builder
------------

The inclusion of modifiers can make building a desired equation of
state somewhat cumbersome. To handle this, we have implemented the
``EOSBuilder`` machinery. ``EOSBuilder`` is a set of functions that
provides a declarative interface for building an equation of state
object.

The EOS Builder functions and types are defined in the
``singularity::EOSBuilder`` namespace. The key function is

.. cpp:function:: EOS EOSBuilder::buildEOS(EOSBuilder::EOSType t, EOSBuilder::params_t base_params, EOSBuilder::modifiers_t modifiers)

* ``EOSBuilder::EOSType`` is an enum class with names that match the various EOS classes defined in :ref:`the models section <models>`; for example, ``EOSBuilder::EOSType::IdealGas``.
* ``EOSBuilder::params_t`` is a dictionary object with some type erasure, which maps strings to the types ``std::string``, ``int``, or ``Real``. It is used to map parameter names to their values for class constructors.
* ``EOSBuilder::modifiers_t`` is a dictionary from the ``EOSModifier`` enum class, which works identically to the ``EOSType`` enum but for modifiers, to ``params_t`` objects, specifying the constructor values for each modifier.

Putting it all together, initializing an ``IdealGas`` with
``EOSBuilder`` looks something like this:

.. code-block:: cpp

  using namespace singularity;
  EOSBuilder::EOSType type = EOSBuilder::EOSType::IdealGas;
  EOSBuilder::modifiers_t modifiers;
  EOSBuilder::params_t base_params, shifted_params, scaled_params;
  base_params["Cv"].emplace<Real>(Cv);
  base_params["gm1"].emplace<Real>(gm1);
  shifted_params["shift"].emplace<Real>(shift);
  scaled_params["scale"].emplace<Real>(scale);
  modifiers[EOSBuilder::EOSModifier::Shifted] = shifted_params;
  modifiers[EOSBuilder::EOSModifier::Scaled] = scaled_params;
  EOS eos = EOSBuilder::buildEOS(type, base_params, modifiers);


.. _eos methods reference section:

Equation of State Methods Reference
------------------------------------

Below the scalar functions are listed. In general, a vector version of
each of these functions exists, which returns void and takes indexers
of each input followed by each output. All of these functions are
available on both host and device (if compiled for a system with a
discrete accelerator).

Functions are named descriptively, and therefore the method names
should be self explanatory. Unless specified, all units are in
cgs. Unless specified, all functions work on device, if the code is
compiled appropriately. The exceptions are constructors,
``GetOnDevice``, and ``Finalize``, all of which are host-only.

.. cpp:function:: Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie, Rela &lambda = nullptr) const;

Returns temperature in Kelvin. Inputs are density in :math:`g/cm^3`
and specific internal energy in :math:`erg/g`. The vector equivalent
of this function is

.. code-block::

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, const int num,
                                       LambdaIndexer &&lambdas) const;

where ``rhos`` and ``sies`` are input arrays and ``temperatures`` is
an output array. ``num`` is the size of those arrays and ``lambdas``
is an optional array of ``lambda`` arrays. In general, every scalar
function that returns a real number given a thermodynamic state has a
vector function with analogous signature. The optional ``lambda``
parameter is always last in the function signature. As they are all
almost exactly analogous to their scalar counterparts, we will mostly
not list the vector functions here.

.. cpp:function:: Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature, Real *lambda=nullptr) const;

returns specific internal energy in :math:`erg/g` given a density in
:math:`g/cm^3` and a temperature in Kelvin.

.. cpp:function:: Real PressureFromDensityTemperature(const Real rho, const Real temperature, Real *lambda = nullptr) const;

returns pressure in Barye given density in :math:`g/cm^3` and temperature in Kelvin.

.. cpp:function:: Real PressureFromDensityInternalEnergy(const Real rho, const Real temperature, Real *lambda = nullptr) const;

returns pressure in Barye given density in :math:`g/cm^3` and specific
internal energy in :math:`erg/g`.

.. cpp:function:: Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature, Real *lambda = nullptr) const;

returns specific heat capacity at constant volume, in units of
:math:`erg/(g K)` in terms of density in :math:`g/cm^3` and
temperature in Kelvin.

.. cpp:function:: Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie, Real *lambda = nullptr) const;

returns specific heat capacity at constant volume, in units of
:math:`erg/(g K)` in terms of density in :math:`g/cm^3` and specific
internal energy in :math:`erg/g`.

.. cpp:function:: Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature, Real *lambda = nullptr) const;

returns the the bulk modulus

.. math::

  B_s = (\partial P/\partial \rho)_s

in units of :math:`g cm^2/s^2` given density in :math:`g/cm^3` and
temperature in Kelvin. For most material models, the square of the
sound speed is given by

.. math::

   c_s^2 = \frac{B_S}{\rho}

Note that for relativistic models,

.. math::

   c_s^2 = \frac{B_S}{w}

where :math:`w = \rho h` for specific entalpy :math:`h` is the
enthalpy by volume. The sound speed may also differ for, e.g., porous
models, where the pressure is less directly correlated with the
density.

.. cpp:function:: Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie, Real *lambda = nullptr) const;

returns the bulk modulus in units of :math:`g cm^2/s^2` given density
in :math:`g/cm^3` and specific internal energy in :math:`erg/g`.

.. cpp:function:: Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature, Real *lambda = nullptr) const;

returns the unitless Gruneisen parameter

.. math::

  \Gamma = \frac{1}{\rho}\left(\frac{\partial P}{\partial \varepsilon}\right)_\rho

given density in :math:`g/cm^3` and temperature in Kelvin.

.. cpp:function:: Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie, Real *lambda = nullptr) const;

returns the unitless Gruneisen parameter given density in
:math:`g/cm^3` and specific internal energy in :math:`erg/g`.

The function

.. cpp:function:: void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod, Real &dpde, Real &dvdt, Real *lambda = nullptr) const;

fills the density, temperature, specific internal energy, pressure,
and thermodynamic derivatives a specifically chosen characteristic
"reference" state. For terrestrial equations of state, this reference
state is probably close to standard density and pressure. For
astrophysical models, it will be chosen to be close to a
representative energy and density scale.

The function

.. cpp:function:: void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod, const unsigned long output, Real *lambda = nullptr) const;

is a a bit of a special case. ``output`` is a bitfield represented as
an unsigned 64 bit number. Quantities such ``pressure`` and
``specific_internal_energy`` can be represented in the ``output``
field by flipping the appropriate bits. There is one bit per
quantity. ``FillEos`` sets all parameters (passed in by reference)
requested in the ``output`` field utilizing all paramters not
requested in the ``output`` flag, which are assumed to be input.

The ``output`` variable uses the same ``thermalqs`` flags as the
``PreferredInput`` method. If an insufficient number of variables are
passed in as input, or if the input is not a combination supported by
a given model, the function is expected to raise an error. The exact
combinations of inputs and ouptuts supported is model
dependent. However, the user will always be able to use density and
temperature or internal energy as inputs and get all other
quantities as outputs.

Methods Used for Mixed Cell Closures
--------------------------------------

Several methods were developed in support of mixed cell closures. In particular:

.. cpp:function:: Real MinimumDensity() const;

and 

.. cpp:function:: Real MinimumTemperature() const;

provide bounds for valid inputs into a table, which can be used by a
root finder to meaningful bound the root search. Similarly,

.. cpp:function:: Real RhoPmin(const Real temp) const;

returns the density at which pressure is minimized for a given
temperature. This is again useful for root finds.
