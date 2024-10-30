.. _python:

Python Bindings
===============

Singularity EOS provides Python bindings which can be enabled with the CMake
``SINGULARITY_BUILD_PYTHON`` option. They provide a 1:1 mapping of the C++ EOS
types and give access to both scalar and vector functions.

Where you build or install your python bindings to must be included in
your python module search path. For example, if you do not install,
this may be ``singularity-eos/build/python``.

.. note::

   At this time, all Python bindings are host-only.

Example
-------

::

   from singularity_eos import IdealGas

   # Parameters for ideal gas
   gm1 = 0.6
   Cv = 2
   eos = IdealGas(gm1, Cv)
   rho = ...
   sie = ...
   P = eos.PressureFromDensityInternalEnergy(rho, sie)
   eos.Finalize()

A more elaborate example can be found in ``examples/get_sound_speed_press.py``.

.. currentmodule:: singularity_eos

Classes
-------
List may not be complete.

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

Class Reference
---------------
List may not be complete.

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

