.. _contributing-doc:

Contributing
=============

Useful code structure information
---------------------------------

The ``Variant`` class and ``mpark::variant``
````````````````````````````````````````````

Work in progress. Things to cover:

*  Type erasure
*  Unified API

The CRTP slass structure and static polymorphism
````````````````````````````````````````````````

Each of the EOS models in ``singularity-eos`` inherits from a base class in
order to centralize default functionality and avoid code duplication. The two
main examples of this are the vector overloads and the ``PTofRE`` scalar lookup
function. In the vector overloads, a simple for loop is used to loop over the
set of states provided to the function and then call the scalar version on each
state. The ``PTofRE`` function is designed to provide a common method for
getting the needing information for a PTE solve from an EOS. Both of these
features are not dependent on the specific EOS for their definition, but in the
case of the vector overloads, they *do* need to access methods in the derived
class.

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

The challenge in making this sort of code performance portable is that on a CPU
this sort of type deduction (i.e. what specific EOS implementaiton to call)
would be done at runtime through vtables. In a heterogenous computing enviroment
though, the vtable may point to member functions that do not exist on the
device, causing a segfault. We could have used a similar technique to the
modifier classes and pass the EOS as a template paramter, but then the vector
function calls could only be achieved by creating vector modifiers of all the
implemented EOS.

Instead, the strategy we decided to use in this case was to implement the
polymorphism at compile time through the `CRTP`_ (curiously recurisve templating
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
