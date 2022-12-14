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

Notes for Contribuors
---------------------------------------

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

    #. Use the preprocessor sparingly. (The same applies for
       non-trivial template metaprogramming.)

    #. Limit pointer use to a single dereference. Avoid pointers of
       pointers when possible.

    #. Be compiler warning aware. Try to address compiler warnings as
       they come up.

.. _ten rules for developing safety-critical code: http://web.eecs.umich.edu/~imarkov/10rules.pdf

* Despite the above rules, ``singularity-eos`` is a modern C++ code
  and both standard template library capabilities and template
  metaprogramming are leveraged frequently. This can sometimes make
  parsing the code difficult. If you see something you don't
  understand, ask. It may be it can be refactored to be more simple or
  better documented.

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
  compiles correctly. Code that doesn't need to run on device does not
  need these decorations.

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

Each of the EOS models in ``singularity-eos`` inherits from a base
class in order to centralize default functionality and avoid code
duplication. The two main examples of this are the vector overloads
and the ``PTofRE`` scalar lookup function. In the vector overloads, a
simple for loop is used to loop over the set of states provided to the
function and then call the scalar version on each state. The
``PTofRE`` function is designed to provide a common method for getting
the needing information for a PTE solve from an EOS. Both of these
features are not dependent on the specific EOS for their definition,
but in the case of the vector overloads, they *do* need to access
methods in the derived class. In both cases, these functions have
default behaviour that may need to be overriden for a given equation
of state.

The vector overloads in the base class take the following form (in pseudocode):

.. code-block:: c++

    template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
    inline void
    TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                         RealIndexer &&temperatures, const int num,
                                         LambdaIndexer &&lambdas) const {
    for (int i = 0; i < num; i++) {
        temperatures[i] = eos.TemperatureFromDensityInternalEnergy(rhos[i],
            sies[i], lambdas[i])

where the base class basically needs to call the implementation of the scalar
lookup in the specific EOS. However, this means that the base class needs to
have knowledge of which class is being derived from it in order to call the
correct EOS implementation.

The standard solution to this problem would be "run-time inheritence,"
where type deduction is performed at run-time. While this is possible
on GPU, it becomes cumbersome, as the user must be very explicit about
class inheritence. Moreover, run-time inheritence relies on
relocatable device code, which is not as performant on device, thanks
to weaker cross-compilation unit optimization. We note that to obtain
full performance on device and to build with compilers that don't
support relocatable device code, the entire library must be made
header-only.

We could have used a similar technique to the modifier classes and
pass the EOS as a template paramter, but then the vector function
calls could only be achieved by creating vector modifiers of all the
implemented EOS.

Instead, the strategy we decided to use in this case was to implement the
polymorphism at compile time through the `CRTP`_ (curiously recurring template
pattern). The basic idea is two-fold:

1.  The base class is templated on the derived class to avoid the need for
    vtables.
2.  The ``*this`` pointer for the base class can be statically cast to that of
    the derived class since the derived class inherits from the base. This is
    only possible because the base class is inherited by the derived class and
    this is known at compile time.

Through template resolution, the compiler can then know exactly which member
functions need to be called at *compile time*. This allows us to write the EOS
implementation in the derived class and have the base class call the appropriate
member function.

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
static polymorphism in all of the EOS classes without having to implement
vector member functions in each class.

Note there are several macros to enable the ``using`` statements if
all the functions in the base class can be used freely.

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
