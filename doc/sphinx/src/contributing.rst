.. _contributing-doc:

Contributing
=============

If you have any trouble with the project, or are interested in
participating, please contact us by creating an issue on the github
repository, or submit a pull request!

Pull request protocol
----------------------

There is a pull reuqest template that will be auto-populated when you
submit a pull request. A pull request should have a summary of
changes. You should also add tests for bugs fixed or new features you
add.

We have a changelog file, ``CHANGELOG.md``. After creating your pull
request, add the relevant change and a link to the PR in the
changelog.

Before a pull request will be merged, the code should be formatted. We
use clang-format for this, pinned to version 12. You can automatically
trigger ``clang-format`` in two ways: first you can run the script
``utils/scripts/format.sh``; second you can type ``make
format_singularity`` after configuring the code with ``clang-format``
discoverable by ``cmake``.

Several sets of tests are triggered on a pull request: a static format
check, a docs buld, and unit tests of analytic models and the stellar
collapse model. These are run through github's CPU infrastructure. We
have a second set of tests run on a wider set of architectures that
also access the Sesame library, which we are not able to make public.

The docs are built but not deployed on PRs from forks, and the
internal tests will not be run automatically. So when the code is
ready for merge, you must ask a project maintainer to trigger the
remaining tests for you.

Interwoven Dependencies
------------------------

``singularity-eos`` depends on several other open-source, Los Alamos
maintained, projects. In particular, ``spiner`` and
``ports-of-call``. If you have issues with these projects, ideally
submit issues on the relevant github pages. However, if you can't
figure out where an issue belongs, no big deal. Submit where you can
and we'll engage with you to figure out how to proceed.

Process for adding a new EOS
----------------------------

The basic process for adding a new EOS can be summarized as

#. Create a new header file for your EOS
#. Add the EOS to the ``full_eos_list`` list of EOS in ``eos.hpp``
#. Create tests for your EOS
#. Create a Fortran interface to initialize your EOS into an array of EOS

In addition to these main steps, there are a couple more that are required if
you would like your EOS to work with our fortran interface, which will be
discussed below.

Step 1: Create a new header file for your EOS
`````````````````````````````````````````````

In general, the best practice is to simply copy an existing EOS file and modify
it for the new EOS. However, there are some subtleties here that are important.

- Parameters for the EOS can be initialzed via an initializer list, and
  additional parameter checking can be done in the constructor.
- Any EOS must have a set of member functions that conform to the general
  :ref:`EOS API<eos methods reference section>`. In essence, these functions are
  :ref:`defined by the <variant section>` ``Variant`` :ref:`class <variant
  section>` as a ``visit`` on the underlying member of the specific EOS type
  contained in the ``variant``. If a new EOS doesn't have an appropriate
  member, a compilation error will be thrown when the ``EOS`` type is used to
  instantiate an instance of the new EOS. This will be discussed more in the
  testing section.
- You may find it useful to define other functions that are specific to that EOS
  but won't be available to the general ``EOS`` type. These can be internal
  checking functions or common calculations that need to be performed for
  multiple types of lookups.
- An analytic EOS needs to be "trivially copiable" in order to use the standard
  ``GetOnDevice()`` function that we use for analytic EOS. In general, analytic
  EOS should only need parameters whose size is known at compile time, so this
  should be fairly straight-forward. Any EOS that needs dynamic memory (e.g.
  a tabular EOS) will need more effort in order to ensure that memory is copied
  correctly over to the device.


Step 2: Add the EOS to the ``full_eos_list`` list of EOS in ``eos.hpp``
````````````````````````````````````````````````````````````````````````

As was mentioned previously, we use the ``Variant`` class and a ``visit``
pattern to achieve compile-time polymorphism on a closed set of types. For
convenience, we provide this closed set in the ``eos.hpp`` file through the
type list, ``full_eos_list``.

For most new EOS, you can simply add the EOS to the ``full_eos_list`` and this
will enable all of the modifiers (with certain exceptions) to instantly work
with your EOS. This would effectively look like

.. code-block:: c++

    static constexpr const auto full_eos_list =
        tl<IdealGas, Gruneisen, JWL, DavisReactants, DavisProducts, MyNewEOS

Note the lack of a trailing comma or closing angle bracket.

If your EOS introduces new dependencies to ``singularity-eos``, then you will
need to create a new flag that enables these dependencies. Then you will need to
wrap the inclusion of your EOS in the ``full_eos_list`` with an appropriate
``#ifdef <my_dependency_flag>`` preprocessor directive. For example, the EOSPAC
EOS needs the ``SINGULARITY_USE_EOSPAC`` flag so the inclusion in the list is
wrapped with ``#ifdef SINGULARITY_USE_EOSPAC``. This might look something like

.. code-block:: c++

    static constexpr const auto full_eos_list =
        tl<IdealGas, Gruneisen, JWL, DavisReactants, DavisProducts
    //
    // ...the other EOS that have dependencies
    //
    #ifdef MY_NEW_DEP_FLAG
           ,
           MyNewEOS
    #endif
           >{};

Note the placement of commas and angle brackets. This example excludes 

Step 3: Create tests for your EOS
`````````````````````````````````

**Important:** this is a sublty that highlights the importance of unit tests!
Since our library is header only, the unit tests are often the only place where
a specific EOS may be instantiated when ``singularity-eos`` is compiled. Unit
tests _must_ make use of the ``EOS`` type, i.e.

.. code-block:: c++

    #include <singularity-eos/eos/eos.hpp>
    EOS my_eos = my_new_eos(parameter1, parameter2, ...)

in order to properly test the functionality of a new EOS. Simply using the
new class as the type such as

.. code-block:: c++

    #include <singularity-eos/eos/eos.hpp>
    auto my_eos = my_new_eos(parameter1, parameter2, ...)

won't ensure that the new EOS is working correctly in singularity with the
static polymorphism of the ``EOS`` type.

You may wish to also design tests that operate on member functions or member
data that is particular to the EOS you have developed, and only for those
specific tests should you instantiate an object whose type is your specific
EOS. Otherwise, use the ``EOS`` object.

Step 4: Fortran interface
`````````````````````````

At this point your new EOS should be usable to any host code written in C++. To
allow the EOS to work with Fortran, an initializer wrapper function needs to be
defined and interfaced with Fortran.

First, the C++ intialization function needs to be named soas to avoid namespace
conflicts. We typically name the initialization functions ``init_sg_<EOSName>``.
For example, the function for initialing an ideal gas looks like

.. code-block:: c++

    int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1,
                         const double Cv, int const *const enabled,
                         double *const vals) {
      assert(matindex >= 0);
      EOS eosi = SGAPPLYMODSIMPLE(IdealGas(gm1, Cv));
      if (enabled[3] == 1) {
        singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                               vals[3], vals[4], vals[5]);
      }
      EOS eos_ = SGAPPLYMOD(IdealGas(gm1, Cv));
      eos[matindex] = eos_.GetOnDevice();
      return 0;
    }

Here the ``*eos`` is a pointer to a container of ``EOS`` objects and the
``matindex`` integer indicates the index at which this EOS will reside in that
container. The ``gm1`` and ``Cv`` inputs are all of the required parameters to
initialize the EOS, while the ``enabled`` and ``vals`` variables are used by
the ``SGAPPLYMOD`` and ``SGAPPLYMODSIMPLE`` macros to apply specific modifiers
to the EOS. The return value of the function is an integer error code that may
or may not be relevant to all EOS.

We also overload the initialization function to make the ``enabled`` and
``vals`` variables effectively optional.

.. code-block:: c++

    int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1,
                         const double Cv) {
      return init_sg_IdealGas(matindex, eos, gm1, Cv, def_en, def_v);
    }

Finally the fortran side, we then define a fortran interface to the C++
initialization function,

.. code-block:: fortran

      interface
        integer(kind=c_int) function &
          init_sg_IdealGas(matindex, eos, gm1, Cv, sg_mods_enabled, &
                           sg_mods_values) &
          bind(C, name='init_sg_IdealGas')
          import
          integer(c_int), value, intent(in)      :: matindex
          type(c_ptr), value, intent(in)         :: eos
          real(kind=c_double), value, intent(in) :: gm1, Cv
          type(c_ptr), value, intent(in)         :: sg_mods_enabled, sg_mods_values
        end function init_sg_IdealGas
      end interface

and a fortran wrapper function to call the C++ function:

.. code-block:: fortran

      integer function init_sg_IdealGas_f(matindex, eos, gm1, Cv, &
                                          sg_mods_enabled, sg_mods_values) &
        result(err)
        integer(c_int), value, intent(in) :: matindex
        type(sg_eos_ary_t), intent(in)    :: eos
        real(kind=8), value, intent(in)   :: gm1, Cv
        integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
        real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
        err = init_sg_IdealGas(matindex-1, eos%ptr, gm1, Cv, &
                               c_loc(sg_mods_enabled), c_loc(sg_mods_values))
      end function init_sg_IdealGas_f

Note that the ``eos`` variable of type ``sg_eos_ary_t`` is just a simple wrapper
for the C pointer to the actual EOS object.

A Note on the EOS Builder
`````````````````````````

The :ref:`EOS Builder <eos builder section>` is a tool that eliminates the need
for chaining together an EOS with a series of modifiers by instead specifing
the parameters and modifications in one function. This convenience comes at the
cost of added development complexity though, and so we do not require a new EOS
to be available for the EOS Builder.

At a basic level though, the EOS needs to be declared in the ``EOSType`` enum
and logic needs to be added to initialze the EOS parameters. More effort may be
needed to make the EOS compatible with modifiers and we point the interested
contributor to the existing EOS as examples.


Notes for Contributors on navigating/developing code features
-------------------------------------------------------------

The CRTP class structure and static polymorphism
````````````````````````````````````````````````

Each of the EOS models in ``singularity-eos`` inherits from a base class in
order to centralize default functionality and avoid code duplication. The two
main examples of this are the vector overloads and the ``PTofRE`` scalar lookup
function. In the vector overloads, a simple for loop is used to iterate over
the set of states provided to the function and then call the scalar version on
each state. The ``PTofRE`` function is designed to provide a common method for
getting the needed information for a PTE solve from an EOS. Both of these
features are general to all types of EOS, but are reliant on specific
implementations of the EOS lookups. In both cases, these functions provide a
default behaviour that we might also want to override for a given equation of
state.

As an example, the vector overloads in the base class take the following form
(in pseudocode):

.. code-block:: c++

    template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
    inline void
    TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                         RealIndexer &&temperatures, const int num,
                                         LambdaIndexer &&lambdas) const {
    for (int i = 0; i < num; i++) {
        temperatures[i] = eos.TemperatureFromDensityInternalEnergy(rhos[i],
            sies[i], lambdas[i])

where the base class needs to call the specific implementation of the scalar
lookup for the particular EOS. However, this means that the base class needs to
have knowledge of which class is being derived from it in order to call the
correct EOS implementation.

The standard solution to this problem would be to deduce the type of the EOS at
runtime (often through virtual functions) and then call the apprporiate member
function in the derived class. While this is possible on GPU, it becomes
cumbersome, as the user must be very explicit about class inheritence.
Moreover, run-time inheritence relies on relocatable device code, which is not
as performant on device, thanks to weaker cross-compilation unit optimization.
We note that to obtain full performance on device and to build with compilers
that don't support relocatable device code, the entire library must be made
header-only.

We could have used a similar technique to the modifier classes and
pass the EOS as a template paramter, but then the vector function
calls could only be achieved by creating vector modifiers of all the
implemented EOS, and the user would have to manually specify that they want to
use a vector verison of the class.

Since we wanted to both leverage C++ function overloading while enabling
compile-time polymorphism, we decided to use the "curiously recurring template
pattern" (`CRTP`_). The basic idea is two-fold:

1.  The base class is templated on the derived class to avoid the need for
    vtables.

2.  The ``*this`` pointer for the base class can be statically cast to that of
    the derived class. This is only possible because the base class is inherited
    by the derived class and this is known at compile time.

Through template resolution, the compiler can then know exactly which member
functions need to be called at *compile time*. This allows us to write the EOS
implementation in the derived class and have common functionality that leverages
these implementations in the base class.

The above example modified to take advantage of the CRTP becomes

.. code-block:: c++

    template <typename CRTP>
    class EosBase {
     public:
      template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
      inline void
      TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                           RealIndexer &&temperatures, const int num,
                                           LambdaIndexer &&lambdas) const {
        for (int i = 0; i < num; i++) {
          temperatures[i] = static_cast<CRTP const &>(*this).TemperatureFromDensityInternalEnergy(
            rhos[i], sies[i], lambdas[i]);
      }
    }

The ``EosBase`` class is templated upon the derived class which is passed via the
`CRTP` template parameter. Then the ``EosBase`` class vector implementation
statically casts its own ``*this`` pointer to that of the derived class in order
to call the specific EOS implementation.

The derived class then needs to look something like

.. code-block:: c++

    class EosImplementation : public EosBase<EosImplementation> {
     public:
      static inline Real TemperatureFromDensityInternalEnergy(
          const Real rho, const Real sie, Real *lambda) const {
        // Specific EOS implementation for returning T(rho, e)
        return temperature;
      }
      using EosBase<EosImplementation>::TemperatureFromDensityInternalEnergy
    }

Note that the ``using`` statement needs to be included in order to properly
overload the scalar functionality with the vector functionality. Otherwise the
vector member function is hidden by the derived class method rather than
overloaded.

With several EOS that all inherit from the ``EosBase`` class, we can achieve
static polymorphism in all of the EOS classes without having to duplicate code
in each class.

Note there are several macros to enable the ``using`` statements if
all the functions in the base class can be used freely. Omitting a ``using``
statement allows the developer to provide a custom implementation of a member
function for that particular EOS.

Also note that any new functionality added to the base class needs to be
mirrored in the :ref:`Variant class <variant section>` so that it is accessable
when using the ``EOS`` type.

.. _CRTP: https://www.fluentcpp.com/2017/05/12/curiously-recurring-template-pattern/

Fast Logs and Approximate Log Gridding
```````````````````````````````````````

When spanning many orders of magnitude, Logarithmic grids are a
natural choice. Even spacing in log space corresponds to exponential
spacing in the original linear space. In other words, the grid spacing
is proportional to the value of the independent variable.

One can perform log-linear or log-log interpolation by simply
converting to log space, interpolating as one normally would, and then
converting back out. Unfortunately, logarithms and exponents are
transcendental functions, meaning they are expensive to compute and it
is thus expensive to transform in and out of log space.

To avoid this issue, we construct a space that is *approximately*
logarithmically spaced, but not quite exactly. The requirements for
this space are that the transformation into and out of this space is
fast to compute, continuous, differentiable, analytically invertible,
and close to taking a logarithm or exponentiation (depending on which
way you're going).

To achieve this, we leverage the internal representation of a floating
point number in the IEE standard. In particular, a floating point
number :math:`x` is represented as a mantissa and an exponent in base
2:

.. math::

   x = m 2^e

for mantissa :math:`m` and exponent :math:`e`. The mantiss is
guaranteed to be on the interval :math:`[1/2, 1)`. The standard
library of most low-level languages provides a performant and portable
routine to pick apart this represnetation, ``frexp``, which given a
number :math:`x`, return :math:`m` and :math:`e`.

The log in base 2 ``lg`` of :math:`x` is then given by the logarithm
of the mantissa plus the exponent:

.. math::

   \lg(x) = \lg(m) + e

Therefore, if we can find a fast, invertible approximation to
:math:`\lg(m)`, we will have achieved our goal. It turns out the
expression

.. math::

   2 (x - 1)

works pretty well, so we use that. (To convince yourself of this note
that for :math:`x=1/2` this expression returns -1 and for :math:`x=1`,
it returns 0, which are the correct values of :math:`\lg(x)` at the
bounds of the interval.) Thus our approximate, invertible expression
for :math:`\lg` is just

.. math::

   2 (m - 1) + e

for the mantissa and exponent extracted via ``frexp``. This differs
from :math:`lg` by a maximum of about 0.1, which translates to at most
a 25 percent difference. As discussed above, however, the function
itself is an exact representation of itself and the difference from
:math:`lg` is acceptable.

To invert, we use the built in function that inverts ``frexp``,
``ldexp``, which combines the mantissa and exponent into the original
floating point representation.

This approach is described in more detail in our `short note`_ on the topic.

.. _Short note: https://arxiv.org/abs/2206.08957
