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
pointers, and ``std::vector<double>`` are two exmaples of what we call
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

For more details on modifiers, see the :ref:`modifiers<modifiers>`
section. If you need a combination of modifiers not supported by
default, we recommend building a custom variant as described above.

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


