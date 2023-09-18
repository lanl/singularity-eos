.. _contributing-doc:

.. _Catch2: https://github.com/catchorg/Catch2/blob/devel/docs/tutorial.md

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
discoverable by ``cmake``. The former script takes two CLI arguments
that may be useful, ``CFM``, which can be set to the path for your
clang-format binary, and ``VERBOSE``, which if set to ``1`` adds
useful output. For example:

.. code-block:: bash

    CFM=clang-format-12 VERBOSE=1 ./util/scripts/format.sh

Several sets of tests are triggered on a pull request: a static format
check, a docs buld, and unit tests of analytic models and the stellar
collapse model. These are run through github's CPU infrastructure. We
have a second set of tests run on a wider set of architectures that
also access the Sesame library, which we are not able to make public.

The docs are built but not deployed on PRs from forks, and the
internal tests will not be run automatically. So when the code is
ready for merge, you must ask a project maintainer to trigger the
remaining tests for you.

Expectations for code review
-----------------------------

From the perspective of the contributor
````````````````````````````````````````

Code review is an integral part of the development process
for ``singularity-eos``. You can expect at least one, perhaps many,
core developers to read your code and offer suggestions.
You should treat this much like scientific or academic peer review.
You should listen to suggestions but also feel entitled to push back
if you believe the suggestions or comments are incorrect or
are requesting too much effort.

Reviewers may offer conflicting advice, if this is the case, it's an
opportunity to open a discussion and communally arrive at a good
approach. You should feel empowered to argue for which of the
conflicting solutions you prefer or to suggest a compromise. If you
don't feel strongly, that's fine too, but it's best to say so to keep
the lines of communication open.

Big contributions may be difficult to review in one piece and you may
be requested to split your pull request into two or more separate
contributions. You may also receive many "nitpicky" comments about
code style or structure. These comments help keep a broad codebase,
with many contributors uniform in style and maintainable with
consistent expectations accross the code base. While there is no
formal style guide for now, the regular contributors have a sense for
the broad style of the project. You should take these stylistic and
"nitpicky" suggestions seriously, but you should also feel free to
push back.

As with any creative endeavor, we put a lot of ourselves into our
code. It can be painful to receive criticism on your contribution and
easy to take it personally. While you should resist the urge to take
offense, it is also partly code reviewer's responsiblity to create a
constructive environment, as discussed below.

Expectations of code reviewers
````````````````````````````````

A good code review builds a contribution up, rather than tearing it
down. Here are a few rules to keep code reviews constructive and
congenial:

* You should take the time needed to review a contribution and offer
  meaningful advice. Unless a contribution is very small, limit
  the times you simply click "approve" with a "looks good to me."

* You should keep your comments constructive. For example, rather than
  saying "this pattern is bad," try saying "at this point, you may
  want to try this other pattern."

* Avoid language that can be misconstrued, even if it's common
  notation in the commnunity. For example, avoid phrases like "code
  smell."

* Explain why you make a suggestion. In addition to saying "try X
  instead of Y" explain why you like pattern X more than pattern Y.

* A contributor may push back on your suggestion. Be open to the
  possibility that you're either asking too much or are incorrect in
  this instance. Code review is an opportunity for everyone to learn.

* Don't just highlight what you don't like. Also highlight the parts
  of the pull request you do like and thank the contributor for their
  effort.

General principle for everyone
```````````````````````````````

It's hard to convey tone in text correspondance. Try to read what
others write favorably and try to write in such a way that your tone
can't be mis-interpreted as malicious.

Interwoven Dependencies
------------------------

``singularity-eos`` depends on several other open-source, Los Alamos
maintained, projects. In particular, ``spiner`` and
``ports-of-call``. If you have issues with these projects, ideally
submit issues on the relevant github pages. However, if you can't
figure out where an issue belongs, no big deal. Submit where you can
and we'll engage with you to figure out how to proceed.

.. note::
   There are scheduled workflows triggered by GitHub actions that will
   automatically check ``spiner`` and ``ports-of-call`` for Spack updates.  If
   detected, the GitHub action bot will create a PR with the necessary changes.

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

The first step is to create a new ``.cpp`` file for testing your new EOS using
the `Catch2 <Catch2_>`_ framework to design your test. We make use of the
Behavior Driven Development (BDD) style for Catch2, so we recommend you do the
same. In general, we recommend you copy the general structure of one of the
existing EOS-specific unit tests.

After creating your tests, you will need to include the ``.cpp`` for your new
test in the ``CMakeLists.txt`` file, 

.. code-block:: cmake

    add_executable(eos_unit_tests
        catch2_define.cpp
        eos_unit_test_helpers.hpp
        test_eos_unit.cpp
        test_eos_gruneisen.cpp
        test_eos_vinet.cpp
        test_my_new_eos.cpp
    )

in order for the test to be compiled. If your EOS requires any special
dependencies, be sure to block off the test using ``#IFDEF`` blocks.

**Important:** this is a subtlety that highlights the importance of unit tests!
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

If you wish to test error handling in your EOS, you may use the macro
``REQUIRE_MAYBE_THROWS``, which is defined in the ``eos_unit_test_helpers.hpp``
header file. This macro will check if your code throws an exception if
compiled for CPU only and otherwise is a no-op. This is intended to combine with
the ``PORTABLE_THROW_OR_ABORT` macro defined in ``ports-of-call``.

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

Some notes on style and code architecture
``````````````````````````````````````````

* ``singularity-eos`` is primarily designed to provide needed equation
  of state functionality to continuum dynamics codes. It isn't
  supposed to provide the most accurate or complete picture of thermal
  or statistical physics. As such the project tries to limit
  capabilities to this scope.

* A major influence on code style and architecture is the
  `ten rules for developing safety-critical code`_, by Gerard Holzmann.
  Safety critical code is code that exists in a context where failure
  implies serious harm. A flight controler on an airplane or
  spacecraft or the microcontroller in a car are examples of
  safety-critical contexts. ``singularity-eos`` is not safety-critical
  but many of the coding habits advocated for by Holzmann produce
  long-lived, easy to understand, easy to parse, and easy to maintain code.
  And we take many of the rules to heart. Here are a few that are most
  relevant to ``singularity-eos``. They have been adapted slightly to 
  our context.

    #. Avoid complex flow constructs such as gotos.

    #. All loops must have fixed bounds. This prevents runaway
       code. (Note this implies that as a general rule, one should use
       ``for`` loops, not ``while`` loops. It also implies one should
       keep recursion to a minimum.)

    #. Heap memory allocation should only be performed at
       initialization. Heap memory de-allocation should only be
       performed at cleanup.

    #. Restrict the length of functions to a single printed page.

    #. Restrict the scope of data to the smallest possible.

    #. Use the preprocessor sparingly.

    #. Limit pointer use to a single dereference. Avoid pointers of
       pointers when possible.

    #. Be compiler warning aware. Try to address compiler warnings as
       they come up.

.. _ten rules for developing safety-critical code: http://web.eecs.umich.edu/~imarkov/10rules.pdf

* ``singularity-eos`` is a modern C++ code
  and both standard template library capabilities and template
  metaprogramming are leveraged frequently. This can sometimes make
  parsing the code difficult. If you see something you don't
  understand, ask. It may be it can be refactored to be more simple or
  better documented.

* As a general rule, to avoid accidental division by zero, use the
  ``robust::ratio(x, y)`` function provided in
  ``singularity-eos/base/robust_utils.hpp`` instead of writing ``x /
  y``.

Performance portability concerns
`````````````````````````````````

``singularity-eos`` is performance portable, meaning it is designed to
run not only on CPUs, but GPUs from a variety of manufacturers,
powered by a variety of device-side development tools such as Cuda,
OpenMP, and OpenACC. This implies several constraints on code
style. Here we briefly discuss a few things one should be aware of.

* **``ports-of-call`` and portability decorators:** Functions that
  should be run on device needs to be decorated with one of the
  following macros: ``PORTABLE_FUNCTION``,
  ``PORTABLE_INLINE_FUNCTION``,
  ``PORTABLE_FORCEINLINE_FUNCTION``. These macros are imported from
  the `ports-of-call`_ library and resolve to the appropriate
  decorations for a given device-side backend such as cuda so the code
  compiles correctly. Code that doesn't need to run on device, 
  such as EOS class constructors, does not need these decorations.

* **Relocatable device code:** It is common in C++ to split code
  between a header file and an implementation file. Functionality that
  is to be called from within loops run on device should not be split
  in this way. Not all accelerator languages support this and the ones
  that do take a performance hit. Instead implement that functionality
  only in a header file and decorate it with
  ``PORTABLE_INLINE_FUNCTION`` or ``PORTABLE_FORCEINLINE_FUNCTION``.

* **Host and device pointers:** Usually accelerators have different
  memory spaces than the CPU they are attached to. So you need to be
  aware that data needs to be copied to an accelerator device to be
  used. If it is not properly copied, the code will likely crash with
  a segfault. In general scalar data such as a single variable (e.g.,
  ``int x``) can be easily and automatically copied to device and you
  don't need to worry about managing it. Arrays and pointers, however,
  are a different story. If you create an array or point to some
  memory on CPU, then you are pointing to a location in memory on your
  CPU. If you try to access it from your accelerator, your code will
  not behave properly. You need to manually copy data from host to
  device in this case. The libraries `ports-of-call`_ and `spiner`_
  offer some functionality for managing arrays on device.

* **Shallow copies:** As a general rule, large
  amount of data stored within an ``EOS`` object should have
  "reference-semantics." This means that if you copy an EOS object, it
  should always be a shallow copy, not a deep copy, unless a deep copy
  is explicitly requested. This is for performance reasons and also to
  simplify the managment of data on device.

* **Real:** The ``Real`` datatype is either a single precision or
  double precision floating point number, depending on how
  `ports-of-call`_ is configured. For most floating point numbers use
  the ``Real`` type. However, be conscious that sometimes you will
  specifically need a single or double precision number, in which case
  you should specify the type as built into the language.

.. _ports-of-call: https://lanl.github.io/ports-of-call/main/index.html

.. _spiner: https://lanl.github.io/spiner/main/index.html

The CRTP slass structure and static polymorphism
````````````````````````````````````````````````

Each of the EOS models in ``singularity-eos`` inherits from a base class in
order to centralize default functionality and avoid code duplication. The
main example of this are the vector overloads.
In the vector overloads, a simple for loop is used to iterate over
the set of states provided to the function and then call the scalar version on
each state. This feature is
general to all types of EOS, but reliant on specific
implementations of the EOS lookups. These functions provide a
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
