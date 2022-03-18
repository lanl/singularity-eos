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

.. _CRTP: https://www.fluentcpp.com/2017/05/12/curiously-recurring-template-pattern/

