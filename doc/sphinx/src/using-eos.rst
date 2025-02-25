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
parallelize your code on-node. It knows nothing at all about
distributed memory parallelism, with one exception, discussed
below. An ``EOS`` object can be copied into any parallel code block by
value (see below) and scalar calls do not attempt any internal
multi-threading, meaning ``EOS`` objects are not thread-safe, but are
compatible with thread safety, assuming the user calls them
appropriately. The main complication is ``lambda`` arrays, which are
discussed below.

The vector ``EOS`` method overloads are a bit different. These are
thread-parallel operations launched by ``singularity-EOS``. They run
in parallel, and ordering between indices of vectors cannot be
assumed. Therefore care must taken in memory layout to avoid race
conditions. The type of parallelism used depends on how
``singularity-eos`` is compiled. If the ``Kokkos`` backend is used,
any parallel dispatch supported by ``Kokkos`` is supported.

A more generic version of the vector calls exists in the
``EvaluateHost`` and ``EvaluateDevice`` methods, which allow the user
to specify arbitrary parallel dispatch models by writing their own
loops. See the relevant section below.

Serialization and shared memory
--------------------------------

While ``singularity-eos`` makes a best effort to be agnostic to
parallelism, it exposes several methods that are useful in a
distributed memory environment. In particular, there are two use-cases
the library seeks to support:

#. To avoid stressing a filesystem, it may desirable to load a table from one thread (e.g., MPI rank) and broadcast this data to all other ranks.
#. To save memory it may be desirable to place tabulated data, which is read-only after it has been loaded from file, into shared memory on a given node, even if all other data is thread local in a distributed-memory environment. This is possible via, e.g., `MPI Windows`_.

Therefore ``singularity-eos`` exposes several methods that can be used
in this context. The function

.. cpp:function:: std::size_t EOS::SerializedSizeInBytes() const;

returns the amount of memory required in bytes to hold a serialized
EOS object. The return value will depend on the underlying equation of
state model currently contained in the object. The function

.. cpp:function:: std::size_t EOS::SharedMemorySizeInBytes() const;

returns the amount of data (in bytes) that a given object can place into shared memory. Again, the return value depends on the model the object currently represents.

.. note::

  Many models may not be able to utilize shared memory at all. This
  holds for most analytic models, for example. The ``EOSPAC`` backend
  will only utilize shared memory if the ``EOSPAC`` version is sufficiently recent
  to support it and if ``singularity-eos`` is built with serialization
  support for ``EOSPAC`` (enabled with
  ``-DSINGULARITY_EOSPAC_ENABLE_SHMEM=ON``).

The function

.. cpp:function:: std::size_t EOS::Serialize(char *dst);

fills the ``dst`` pointer with the memory required for serialization
and returns the number of bytes written to ``dst``. The function

.. cpp:function:: std::pair<std::size_t, char*> EOS::Serialize();

allocates a ``char*`` pointer to contain serialized data and fills
it.

.. warning::

  Serialization and de-serialization may only be performed on objects
  that live in host memory, before you have called
  ``eos.GetOnDevice()``. Attempting to serialize device-initialized
  objects is undefined behavior, but will likely result in a
  segmentation fault.

The pair is the pointer and its size. The function

.. code-block:: cpp

  std::size_t EOS::DeSerialize(char *src,
                const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS)

Sets an EOS object based on the serialized representation contained in
``src``. It returns the number of bytes read from ``src``. Optionally,
``DeSerialize`` may also write the data that can be shared to a
pointer contained in ``SharedMemSettings``. If you do this, you must
pass this pointer in, but designate only one thread per shared memory
domain (frequently a node or socket) to actually write to this
data. ``SharedMemSettings`` is a struct containing a ``data`` pointer
and a ``is_domain_root`` boolean:

.. code-block:: cpp

  struct SharedMemSettings {
    SharedMemSettings();
    SharedMemSettings(char *data_, bool is_domain_root_)
        : data(data_), is_domain_root(is_domain_root_) {}
    char *data = nullptr; // defaults
    bool is_domain_root = false;
  };

The ``data`` pointer should point to a shared memory allocation. The
``is_domain_root`` boolean should be true for exactly one thread per
shared memory domain.

For example you might call ``DeSerialize`` as

.. code-block:: cpp

  std::size_t read_size = eos.DeSerialize(packed_data,
                singularity::SharedMemSettings(shared_data,
                                               my_rank % NTHREADS == 0));
  assert(read_size == write_size); // for safety

.. warning::

  Note that for equation of state models that have dynamically
  allocated memory, ``singularity-eos`` reserves the right to point
  directly at data in ``src``, so it **cannot** be freed until you
  would call ``eos.Finalize()``. If the ``SharedMemSettings`` are
  utilized to request data be written to a shared memory pointer,
  however, you can free the ``src`` pointer, so long as you don't free
  the shared memory pointer.

Putting everything together, a full sequence with MPI might look like this:

.. code-block:: cpp

  singularity::EOS eos;
  std::size_t packed_size, shared_size;
  char *packed_data;
  if (rank == 0) { // load eos object
    eos = singularity::StellarCollapse(filename);
    packed_size = eos.SerializedSizeInBytes();
    shared_size = eos.SharedMemorySizeInBytes();
  }

  // Send sizes
  MPI_Bcast(&packed_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&spacked_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  // Allocate data needed
  packed_data = (char*)malloc(packed_size);
  if (rank == 0) {
    eos.Serialize(packed_data);
    eos.Finalize(); // Clean up this EOS object so it can be reused.
  }
  MPI_Bcast(packed_data, packed_size, MPI_BYTE, 0, MPI_COMM_WORLD);
  
  // the default doesn't do shared memory.
  // we will change it below if shared memory is enabled.
  singularity::SharedMemSettings settings = singularity::DEFAULT_SHMEM_STNGS;

  char *shared_data;
  char *mpi_base_pointer;
  int mpi_unit;
  MPI_Aint query_size;
  MPI_Win window;
  MPI_Comm shared_memory_comm;
  int island_rank, island_size; // rank in, size of shared memory region
  if (use_mpi_shared_memory) {
    // Generate shared memory comms
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared_memory_comm);
    // rank on a region that shares memory
    MPI_Comm_rank(shared_memory_comm, &island_rank);
    // size on a region that shares memory
    MPI_COMM_size(shared_memory_comm, &island_size);

    // Create the MPI shared memory object and get a pointer to shared data
    // this allocation is a collective and must be called on every rank.
    // the total size of the allocation is the sum over ranks in the shared memory comm
    // of requested memory. So it's valid to request all you want on rank 0 and nothing
    // on the remaining ranks.
    MPI_Win_allocate_shared((island_rank == 0) ? shared_size : 0,
                            1, MPI_INFO_NULL, shared_memory_comm, &mpi_base_pointer,
                            &window);
    // This gets a pointer to the shared memory allocation valid in local address space
    // on every rank
    MPI_Win_shared_query(window, MPI_PROC_NULL, &query_size, &mpi_unit, &shared_data);
    // Mutex for MPI window. Writing to shared memory currently allowed.
    MPI_Win_lock_all(MPI_MODE_NOCHECK, window);
    // Set SharedMemSettings
    settings.data = shared_data;
    settings.is_domain_root = (island_rank == 0);
  }
  eos.DeSerialize(packed_data, settings);
  if (use_mpi_shared_memory) {
    MPI_Win_unlock_all(window); // Writing to shared memory disabled.
    MPI_Barrier(shared_memory_comm);
    free(packed_data);
  }

In the case where many EOS objects may be active at once, you can
combine serialization and comm steps. You may wish to, for example,
have a single pointer containing all serialized EOS's. Same for the
shared memory. ``singularity-eos`` provides machinery to do this in
the ``singularity-eos/base/serialization_utils.hpp`` header. This
provides a helper struct, ``BulkSerializer``:

.. code-block:: cpp

  template<typename Container_t, Resizer_t = MemberResizer>
  singularity::BulkSerializer

which may be initialized by a collection of ``EOS`` objects or by
simply assigning (or constructing) its member field, ``eos_objects``
appropriately. An example ``Container_t`` might be
``std::vector<EOS>``. A specialization for ``vector`` is provided as
``VectorSerializer``. The ``Resizer_t`` is a functor that knows how to
resize a collection. For example, the ``MemberResizor`` functor used
for ``std::vector``s

.. code-block:: cpp

  struct MemberResizer {
    template<typename Collection_t>
    void operator()(Collection_t &collection, std::size_t count) {
      collection.resize(count);
    }
  };

which will work for any ``stl`` container with a ``resize`` method.

The ``BulkSerializer`` provides all the above-described serialization
functions for ``EOS`` objects: ``SerializedSizeInBytes``,
``SharedMemorySizeInBytes``, ``Serialize``, and ``DeSerialize``, but
it operates on all ``EOS`` objects contained in the container it
wraps, not just one. Example usage might look like this:

.. code-block:: cpp

  int packed_size, shared_size;
  singularity::VectorSerializer<EOS> serializer;
  if (rank == 0) { // load eos object
    // Code to initialize a bunch of EOS objects into a std::vector<EOS>
    /*
       Initialization code goes here
     */
    serializer = singularity::VectorSerializer<EOS>(eos_vec);
    packed_size = serializer.SerializedSizeInBytes();
    shared_size = serializer.SharedMemorySizeInBytes();
  }

  // Send sizes
  MPI_Bcast(&packed_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&packed_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  // Allocate data needed
  packed_data = (char*)malloc(packed_size);
  if (rank == 0) {
    serializer.Serialize(packed_data);
    serializer.Finalize(); // Clean up all EOSs owned by the serializer
  }
  MPI_Bcast(packed_data, packed_size, MPI_BYTE, 0, MPI_COMM_WORLD);
  
  singularity::SharedMemSettings settings = singularity::DEFAULT_SHMEM_STNGS;
  // same MPI declarations as above
  if (use_mpi_shared_memory) {
    // same MPI code as above including setting the settings
    settings.data = shared_data;
    settings.is_domain_root = (island_rank == 0);
  }
  singularity::VectorSerializer<EOS> deserializer;
  deserializer.DeSerialize(packed_data, settings);
  if (use_mpi_shared_memory) {
    // same MPI code as above
  }
  // extract each individual EOS and do something with it
  std::vector<EOS> eos_host_vec = deserializer.eos_objects;
  // get on device if you want
  for (auto EOS : eos_host_vec) {
    EOS eos_device = eos.GetOnDevice();
    // ...
  }

It is also possible to (with care) mix serializers... i.e., you might
serialize with a ``VectorSerializer`` and de-serialize with a
different container, as all that is required is that a container have
a ``size``, provide iterators, and be capable of being resized.

.. warning::

  Since EOSPAC is a library, DeSerialization is destructive for EOSPAC
  and may have side-effects.

.. _`MPI Windows`: https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node311.htm

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

.. code-block:: cpp

   auto my_ideal_gas = my_eos.get<singularity::IdealGas>();

This will give you access to methods and fields which may be unique to
a class but not shared by the ``Variant``.

The EOS model also allows some host-side introspection. The method

.. cpp:function:: static std::string EosType();

returns a string representing the equation of state an ``EOS`` object
currently is. For example:

.. code-block:: cpp

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

Accessors and Indexers
-----------------------

Many functions in ``singularity-eos`` accept **accessors**, also
called **indexers**. An accessor is any object with a square bracket
operator. One-dimensional arrays, pointers, and
``std::vector<double>`` are all examples of what we call an
accessor. However, the value of an accessor is it doesn't have to be
an array. You can create an accessor class that wraps your preferred
memory layout, and ``singularity-eos`` will handle it
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

The Vector API and the ``lambda`` optional arguments all use
accessors, as discussed below.

Vector and Scalar API
----------------------

Most ``EOS`` methods have both scalar and vector overloads, where the
scalar version returns a value, and the vector version modifies an
array. By default the vector version is called from host on device (if
``singularity-eos`` was compiled for device).

The vector API is templated to accept accessors. We do note, however,
that vectorization may suffer if your underlying data structure is not
contiguous in memory.

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

The Evaluate Methods
~~~~~~~~~~~~~~~~~~~~~~

A pair of special call related to the vector calls are the
``EvaluateHost`` and ``EvaluateDevice`` methods. These methods request
the EOS object to evaluate almost arbitrary code, but in a way where
the type of the underlying EOS object is resolved *before* this
arbitrary code is evaluated. This means the code required to resolve
the type of the variant is only executed *once* per ``Evaluate``
call. This can enable composite EOS calls, non-standard vector calls,
and vector calls with non-standard loop structure.

The ``EvaluateHost`` call has the signature

.. code-block:: cpp

  template<typename Functor_t>
  void Evaluate(Functor_t f);

and the ``EvaluateDevice`` method has the signature

.. code-block:: cpp

  template<typename Functor_t>
  PORTABLE_INLINE_FUNCTION
  void Evaluate(Functor_t f);


where a ``Functor_t`` is a class that *must* provide a ``void
operator() const`` method templated on EOS type. Alternatively, you
may use an anonymous function with an `auto` argument as the input,
e.g.,

.. code-block:: cpp

   // equivalent to [=], but with device markings
   eos.Evaluate(PORTABLE_LAMBDA(auto eos) { /* my code snippet */ });

.. warning::

  It can be dangerous to use functors with side-effects. Especially
  with GPUs it can produce very unintuitive behaviour. We recommend
  you only make the ``operator()`` non-const if you really know what
  you're doing. And in the anonymous function case, we recommend you
  capture by value, not reference. ``EvaluateDevice`` does not support
  side effects at all and you must pass your functors in by value in
  that case.

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
  eos.EvaluateHost(my_op);

  // The above two lines could have been called "in-one" with:
  // eos.EvaluateHost(CheckPofRE(P, rho, sie, N));

Alternatively, you could eliminate the functor and use an anonymous
function with:

.. code-block:: cpp

  eos.EvaluateHost([=](auto eos) {
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
parameter, which is an accessor as discussed above. ``lambda[i]``
should return a real number unless ``lambda==nullptr``. Unless
specified in :ref:`the models section <models>`, this parameter does
nothing, and the default type is ``Real*`` with a default value of
``nullptr``

However, some models require or benefit from additional
information. For example models with internal root finds can leverage
initial guesses and models with composition mixing parameters may need
additional input to return a meaningful state.

``EOS`` models are introspective and can provide the desired/required
size of the lambda array with:

.. cpp:function:: int EOS::nlambda()

which is the desired size of the ``lambda`` array per scalar call. For
vector calls, there should be one such accessor per grid point. A
vector accessor for ``lambda`` should return an accessor at each
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

Named Arguments to Lambda Indexers
-----------------------------------

Lambdas support a more free-form kind of indexer. In particular, you
may define indexers that take **names** via a type system. For
example, you can write code such as:

.. code-block:: cpp

  Lambda_t lambda;
  lambda[MeanIonizationState()] = zbar;
  Real P = eos.PressureFromDensityTemperature(rho, T, lambda);

The available types currently supported by default are:

.. code-block:: cpp

  struct singularity::IndexableTypes::MeanIonizationState;
  struct singularity::IndexableTypes::LogDensity;
  struct singularity::IndexableTypes::LogTemperature;
  struct singularity::IndexableTypes::MeanAtomicMass;
  struct singularity::IndexableTypes::MeanAtomicNumber;
  struct singularity::IndexableTypes::ElectronFraction;

However if you are not limited to these types. Any type will do and
you can define your own as you like. For example:

.. code-block::

  namespace IndexableTypes {
  struct MyLambdaParameter;
  }

To use an indexable type, you must define an indexer with an overload
of the ``[]`` operator that takes your type. For example:

.. code-block:: cpp

  class MyLambda_t {
  public:
    MyLambda_t() = default;
    PORTABLE_FORCEINLINE_FUNCTION
    Real &operator[](const std::size_t idx) {
      return data_[idx];
    }
    PORTABLE_FORCEINLINE_FUNCTION
    Real &operator[](const MeanIonizationState &zbar) {
      return data_[2];
    }
  private:
    std::array<Real, 3> data_;
  };

which might be used as

.. code-blocK:: cpp

  MyLambda_t lambda;
  lambda[0] = lRho;
  lambda[1] = lT;
  lambda[MeanIonizationState()] = Z;

where ``MeanIonizationState`` is shorthand for index 2, since you
defined that overload. To more easily enable mixing and matching
integer-based indexing with type-based indexing, the function

.. code-block:: cpp

  template<typename Name_t, typename Indexer_t>
  PORTABLE_FORCEINLINE_FUNCTION
  Real &Get(Indexer_t &&lambda, std::size_t idx = 0);

will return a reference to the value at named index ``Name_t()`` if
that overload is defined in ``Indexer_t`` and otherwise return a
reference at index ``idx``.

As a convenience tool, the struct

.. code-block:: cpp

  template<typename... Ts>
  singularity::IndexableTypes::VariadicIndexer;

automatically defines an indexer that accepts all named indices in the
variadic list ``Ts...`` and also integer indexing. It's a fixed-size
array under the hood.

You can check if an equation of state is compatible with a given named
index type by calling

.. code-block:: cpp

   eos.NeedsLambda(NamedIndex())

which returns ``true`` if the EOS is compatible with ``NamedIndex``
and ``false`` otherwise.

All of the above functionality is available in
the header file ``singularity-eos/base/indexable_types.hpp``.

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

A modified equation of state can be built up iteratively. To check if
the equation of state currently stored in the variant can modified,
you may call

.. code-block:: cpp

  bool ModifiedInVariant<Mod>() const;

where ``Mod`` is the type of the modifier you want to apply, for
example ``ShiftedEOS``. If this function returns true, then you can
apply a modifier with the function

.. code-block:: cpp

  Variant Modify<Mod>(Args &&..args) const;

where again ``Mod`` is the modifier you wish to apply, and ``args``
are the arguments to the constructor for that modifier, e.g., the
shift. For example, one might build up a shifted or scaled eos with a
code block like this:

.. code-block:: cpp

  using namespace singularity;
  EOS eos = IdealGas(gm1, cv);
  if (do_shift) {
    eos = eos.template Modify<ShiftedEOS>(shift);
  }
  if (do_scale) {
    eos = eos.template Modify<ScaledEOS>(scale);
  }

Relevant to the broad ``singularity-eos`` API, EOS models provide
introspection. To check if an EOS is modified, call

.. cpp:function:: bool IsModified() const;

This will return ``true`` for a modified model and ``false``
otherwise. Modifiers can also be undone. To get a completely
unmodified EOS model, call

.. cpp:function:: auto GetUnmodifiedObject();

The return value here will be either the type of the ``EOS`` variant
type or the unmodified model (for example ``IdealGas``), depending
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

Modifiers and Lambdas
-----------------------

Modifiers may require lambdas. When this is the case, the lambda
required by the modifier is appended to the end of the lambda
indexer. For example, the ``StellarCollapse`` EOS model requires
``nlambda=2``. The ``ZSplitI`` modifier rquires
``nlambda=1``. Together, ``ZSplitI<StellarCollapse>`` requires a
lambda indexer of length 3, an the ordering is two parameters for
``StellarCollapse`` first, and then the parameter required by
``ZSplitI``.

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
* ``thermalqs::do_lambda``
* ``thermalqs::all_values``

however, most EOS models only specify that they prefer density and
temperature or density and specific internal energy.

.. note::

   The ``thermalqs::do_lambda`` flag is a bit special. It specifies that
   eos-specific operations are to be performed on the additional
   quantities passed in through the ``lambda`` variable.

.. _eos builder section:

EOS Builder
------------

The iterative construction of modifiers described above and in the
:ref:`modifiers<modifiers>` section is object oriented. For
convenience, we also provide a procedural, dispatch-based approach in
the ``EOSBuilder`` namespace and header. The key function is

.. code-block:: cpp

  template <template <class> typename Mod, typename... Ts, typename... Args>
  singularity::Variant<Ts...> Modify(const Variant<Ts...> &eos, Args &&...args);

where ``Mod`` is an EOS modifier, ``Variant`` is either your
user-defined custom EOS variant type, or the pre-defined ``EOS`` type,
the ``eos`` object is an EOS you'd like to modify (stored as a
variant), and ``args`` are the additional arguments to the constructor
of ``Mod`` beyond the object to modify. For example, initializing an
``IdealGas`` equation of state that is optionally shifted and scaled
might look something like this:

.. code-block:: cpp

  using namespace singularity;
  EOS eos = IdealGas(gm1, cv);
  if (do_shift) {
    eos = EOSBuilder::Modify<ShiftedEOS>(eos, shift);
  }
  if (do_scale) {
    eos = EOSBuilder::Modify<ScaledEOS>(eos, scale);
  }

.. _eos methods reference section:

CheckParams
------------

You may check whether or not an equation of state object is
constructed self-consistently and ready for use by calling

.. cpp:function:: void CheckParams() const;

which raise an error and/or print an equation of state specific error
message if something has gone wrong. Most EOS constructors and ways of
building an EOS call ``CheckParams`` by default.

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

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Indexer_t &&lambda = nullptr) const;

Returns temperature in Kelvin. Inputs are density in :math:`g/cm^3`
and specific internal energy in :math:`erg/g`. The vector equivalent
of this function is

.. code-block:: cpp

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

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                             Indexer_t &&lambda = nullptr) const;

returns specific internal energy in :math:`erg/g` given a density in
:math:`g/cm^3` and a temperature in Kelvin.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                       Indexer_t &&lambda = nullptr) const;

returns pressure in Barye given density in :math:`g/cm^3` and temperature in Kelvin.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real PressureFromDensityInternalEnergy(const Real rho, const Real temperature,
                                          Indexer_t &&lambda = nullptr) const;

returns pressure in Barye given density in :math:`g/cm^3` and specific
internal energy in :math:`erg/g`.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                           Indexer_t &&lambda = nullptr) const;

returns specific heat capacity at constant volume, in units of
:math:`erg/(g K)` in terms of density in :math:`g/cm^3` and
temperature in Kelvin.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                              Indexer_t &&lambda = nullptr) const;

returns specific heat capacity at constant volume, in units of
:math:`erg/(g K)` in terms of density in :math:`g/cm^3` and specific
internal energy in :math:`erg/g`.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                          Indexer_t &&lambda = nullptr) const;

returns the the bulk modulus

.. math::

  B_s = \rho (\partial P/\partial \rho)_s

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

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Indexer_t &&lambda = nullptr) const;

returns the bulk modulus in units of :math:`g cm^2/s^2` given density
in :math:`g/cm^3` and specific internal energy in :math:`erg/g`.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                             Indexer_t &&lambda = nullptr) const;

returns the unitless Gruneisen parameter

.. math::

  \Gamma = \frac{1}{\rho}\left(\frac{\partial P}{\partial \varepsilon}\right)_\rho

given density in :math:`g/cm^3` and temperature in Kelvin.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                                Indexer_t &&lambda = nullptr) const;

returns the unitless Gruneisen parameter given density in
:math:`g/cm^3` and specific internal energy in :math:`erg/g`.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                      Indexer_t &&lambda = nullptr) const;

returns the entropy as a function of density in :math:`g/cm^3` and
temperature in Kelvin.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real EntropyFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Indexer_t &&lambda = nullptr) const;

returns the entropy as a function of density in :math:`g/cm^3` and
specific internal energy in :math:`erg/g`.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real GibssFreeEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                              Indexer_t &&lambda = nullptr) const;

returns the Gibbs free energy as a function of density in :math:`g/cm^3` and
temperature in Kelvin.

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   Real GibbsFreeEnergyFromDensityInternalEnergy(const Real rho, const Real sie,
                                                 Indexer_t &&lambda = nullptr) const;


returns the Gibbs free energy as a function of density in :math:`g/cm^3` and
specific internal energy in :math:`erg/g`.

.. warning::

  Not all equations of state provide entropy and Gibbs free
  energy. These are coupled, however, so if one is provided, the other
  will be too. If you call an entropy for a model that does not
  provide it, ``singularity-eos`` will return an error.

The function

.. code-block:: cpp

   template<typename Indexer_t Real*>
   void MeanAtomicMassFromDensityTemperature(const Real rho, const Real T,
                                             Indexer_t &&lambda = nullpter) const;

returns the mean atomic mass (i.e., the number of nucleons) of a
material given density in :math:`g/cm^3` and temperature in
Kelvin. The reason this is allowed to vary with density and
temperature is that some equations of state, such as the Stellar
Collapse and Helmholtz equations of state encapsulate reactive flows
where the average nucleus may depend on thermodynamic variables. For
most materials, however, this is not the case and a convenience
function that drops the dependence is available:

.. code-block:: cpp

   Real MeanAtomicMass() const;

The function

.. code-block:: cpp

   template<typename Indexer_t Real*>
   void MeanAtomicNumberFromDensityTemperature(const Real rho, const Real T,
                                               Indexer_t &&lambda = nullpter) const;

returns the mean atomic number (i.e., the number of protons in the
nucleus) of a material given density in :math:`g/cm^3` and temperature
in Kelvin. The reason this is allowed to vary with density and
temperature is that some equations of state, such as the Stellar
Collapse and Helmholtz equations of state encapsulate reactive flows
where the average nucleus may depend on thermodynamic variables. For
most materials, however, this is not the case and a convenience
function that drops the dependence is available:

.. code-block:: cpp

   Real MeanAtomicNumber() const;



.. warning::

  For materials where the mean atomic mass and number **do** vary with
  density and temperature, the convenience call without this
  dependence will produce an error.

The function

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                               Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                               Indexer_t &&lambda = nullptr) const;

fills the density, temperature, specific internal energy, pressure,
and thermodynamic derivatives a specifically chosen characteristic
"reference" state. For terrestrial equations of state, this reference
state is probably close to standard density and pressure. For
astrophysical models, it will be chosen to be close to a
representative energy and density scale.

The function

.. code-block:: cpp

   template <typename Indexer_t = Real*>
   void FillEos(Real &rho, Real &temp, Real &energy,
                Real &press, Real &cv, Real &bmod,
                const unsigned long output,
                Indexer_t &&lambda = nullptr) const;

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

Several methods were developed in support of mixed cell closures. In particular the function

.. cpp:function:: Real MinimumDensity() const;

the function

.. cpp:function:: Real MinimumTemperature() const;

and the function

.. cpp:function:: Real MaximumDensity() const;

provide bounds for valid inputs into a table, which can be used by a
root finder to meaningful bound the root search.

.. warning::

  For unbounded equations of state, ``MinimumDensity`` and
  ``MinimumTemperature`` will return zero, while ``MaximumDensity``
  will return a very large finite number. Which number you get,
  however, is not guaranteed. You may wish to apply more sensible
  bounds in your own code.

Similarly,

.. cpp:function:: Real RhoPmin(const Real temp) const;

returns the density at which pressure is minimized for a given
temperature. The function

.. cpp:function:: Real MinimumPressure() const;

provides the minimum pressure an equation of state supports, which may
be the most negative tension state. The function

.. cpp:function:: Real MaximumPressureAtTemperature(const Real temp) const;

provides a maximum possible pressure an equation of state supports at
a given temperature. (Most models are unbounded in pressure.) This is
again useful for root finds.

The function

.. code-block:: cpp

  template <typename Indexer_t = Real*>
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Indexer_t &&lambda, Real &rho, Real &sie) const;

is designed for working in Pressure-Temperature space. Given a
pressure ``press`` and temperature ``temp``, it sets a density ``rho``
and specific internal energy ``sie``. The ``lambda`` is optional and
defaults to a ``nullptr``.

Typically this operation requires a root find. You may pass in an
initial guess for the density ``rho`` in-place and most EOS models
will use it.

.. warning::

  Pressure is not necessarily monotone in density and it may be double
  valued. Thus you are not guaranteed to find the correct root and the
  value of your initial guess may determine correctness. The fact that
  ``rho`` may be used as an initial guess means you **must** pass in
  an initialized variable, even if it is zero-initialized. Do not pass
  uninitialized memory into this function!

