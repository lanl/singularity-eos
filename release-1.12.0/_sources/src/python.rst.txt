.. _python:

Python Bindings
===============

Singularity EOS provides Python bindings which can be enabled with the CMake
``SINGULARITY_BUILD_PYTHON`` option. They provide a 1:1 mapping of the C++ EOS
types and give access to both scalar and vector functions.

Where you build or install your python bindings to must be included in
your python module search path. For example, if you do not install,
this may be ``singularity-eos/builddir/python``.

.. note::

   At this time, all Python bindings are host-only.

Example
-------

::

   import singularity_eos as eos

   with eos.context():
      # Parameters for ideal gas
      gm1 = 0.6
      Cv = 2
      ideal_gas = eos.IdealGas(gm1, Cv)
      rho = ...
      sie = ...
      P = ideal_gas.PressureFromDensityInternalEnergy(rho, sie)
      ideal_gas.Finalize()

   # Leaving the context does not finalize the shared runtime.
   eos.finalize()

which should work for both scalar and vector (via numpy array)
calls. A more elaborate example can be found in
``example/get_sound_speed_press.py``.

Runtime Lifecycle
-----------------

The module owns a single runtime context, returned by :func:`context`. Entering
the context calls :func:`initialize`, but leaving the context does **not** call
:func:`finalize`. This lets EOS objects and work created in separate ``with``
blocks share the same runtime safely.

When the bindings are built with Kokkos support, :func:`initialize` first
checks whether Kokkos is already initialized. If another library initialized
Kokkos, Singularity EOS uses that runtime but does not take ownership of it;
in that case, :func:`finalize` does not shut Kokkos down. If Singularity EOS
starts Kokkos, it records ownership and :func:`finalize` shuts it down.

Calls to :func:`initialize` are idempotent while the runtime is active. Kokkos
cannot be initialized again after it has been finalized, so calling
:func:`initialize` after finalization raises an exception. When the bindings
are built without Kokkos, both lifecycle functions are no-ops.

The runtime may also be managed explicitly without a context manager:

::

   import singularity_eos as eos

   eos.initialize()
   try:
      ideal_gas = eos.IdealGas(0.6, 2.0)
      # Use ideal_gas here.
      ideal_gas.Finalize()
   finally:
      eos.finalize()

Call :func:`finalize` only after all Kokkos-dependent work and cleanup are
complete. An EOS object's ``Finalize()`` method releases resources owned by
that EOS object; it is separate from the module-level :func:`finalize`, which
manages the shared runtime.

.. currentmodule:: singularity_eos

Classes
-------
List may not be complete.

 * :class:`SingularityContext`
 * :class:`IdealGas`
 * :class:`Gruneisen`
 * :class:`JWL`
 * :class:`DavisReactants`
 * :class:`DavisProducts`

Modifiers
---------

Similar to what is described in :doc:`modifiers`, the Python bindings allow you to create modified
versions of EOS with modifiers. Beware that the Python variants follow the same
rules as the C++ modifiers, so not all combinations are possible.

 * :func:`Shifted`
 * :func:`Scaled`
 * :func:`BilinearRamp`
 * :func:`Relativistic`
 * :func:`UnitSystem`

To create a modified EOS, simply pass an existing EOS to a modifier function
along with any modifier arguments.

::

   from singularity_eos import IdealGas, Scaled, Shifted

   # Parameters for ideal gas
   gm1 = 0.6
   Cv = 2
   eos = Scaled(Shifted(IdealGas(gm1, Cv), 1), 1)

.. note::
   
   While you are operating with Python types during construction, the final EOS
   object will be backed by a pure C++ type. E.g., the Python expression
   ``Scaled(Shifted(IdealGas(gm1, Cv), shift), scale)`` will return a Python object
   that wraps the ``ScaledEOS<ShiftedEOS<IdealGas>>`` C++ type.

API Reference
-------------
List may not be complete.

.. autofunction:: initialize

.. autofunction:: finalize

.. autofunction:: context

.. autoclass:: SingularityContext

.. autoclass:: EOSState
   :members:

.. autoclass:: IdealGas
   :members:

.. autoclass:: Gruneisen
   :members:

.. autoclass:: JWL
   :members:

.. autoclass:: DavisReactants
   :members:

.. autoclass:: DavisProducts
   :members:

Modifier Reference:
-------------------
List may not be complete.

.. autofunction:: Shifted

.. autofunction:: Scaled

.. autofunction:: BilinearRamp

.. autofunction:: Relativistic

.. autofunction:: UnitSystem

This file was made in part with generative AI.
