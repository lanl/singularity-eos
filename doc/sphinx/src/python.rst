.. _python:

Python Bindings
===============

Singularity EOS provides Python bindings which can be enabled with the CMake
``SINGULARITY_BUILD_PYTHON`` option. They provide a 1:1 mapping of the C++ EOS
types and give access to both scalar and vector functions.

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

Classes
-------

.. currentmodule:: singularity_eos

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
