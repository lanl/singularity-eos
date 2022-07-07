.. _using-closures:

Mixed Cell Closures
====================

In Eulerian multi-material fluid simulations, a single grid point or
finite volume cell may have more than one material in it. In this
situation, one must decide how to compute thermodynamic quantities in
that cell. One choice, though not necessarily the optimal choice, is
to assume pressure-temperature equilibrium (PTE). This implies that there's
a single pressure in the cell and a single temperature, which are the
same for all materials.

``singularity-eos`` provides several methods for finding a PTE
soluution. These methods differ in what they treat as independent
variables, and thus what precise system of equations they solve
for. However they all share some common characteristics.

Governing Equations and Assumptions
------------------------------------

Given a set of :math:`N` materials ranging from :math:`i = 0` to
:math:`N-1`, each material can, in principle, have a distinct pressure
:math:`P_i`, a distinct temperature :math:`T_i`, a distinct material
density :math:`\rho_i`, and a distinct specific material energy
:math:`\varepsilon_i`. For some finite total volume :math:`V`, each
material also occupies some fraction of that volume given by the
*volume fraction* :math:`f_i` such that

.. math::

  \sum_{i=0}^{N - 1} f_i = 1

The average density in a given material-occupied volume is

.. math::

  \bar{\rho} = \rho_i f_i

and thus the density averaged over all materials is

.. math::

  \rho = \sum_{i=0}^{N - 1} \bar{\rho} = \sum_{i=0}^{N-1} \rho_i f_i

Conservation of energy implies that

.. math::

  u = \rho \varepsilon = \sum_{i = 0}^{N - 1} \rho_i \varepsilon_i

where :math:`u = E/V` is the energy density by volume for total energy
:math:`E` and total volume :math:`V`, and :math:`\varepsilon` is the
total specific internal energy within that volume.

The assumption of pressure equilibrium implies that

.. math::

  P = P_0 = P_1 = \ldots = P_{N - 1}

where each pressure is computed as either

.. math::

  P_i = P_i(\rho_i, T_i)

.. math::

  P_i = P_i(\rho_i, \varepsilon_i)

depending on the treatment.

In ``singularity-eos`` the :math:`N` volume fractions are treated as
unknowns, for which one must solve, and the fact that the volume
fractions sum to 1 implies one constraint. At this point, we have 1
constraint and :math:`N` unknowns (the volume fractions). To guarantee
uniqueness (but not existence) of a PTE solution, we must find a way
to have at least as many equations as unknowns. There are several
choices:

The Density-Energy Formulation
---------------------------------

One choice is to treat volume fractions and material energies as
independent quantities. The material energies provide :math:`N` more
unknowns. The equality of pressures provides :math:`N-1` additional
constraints. Additionally, the euqality of material temperatures, evaluated as

.. math::

  T = T_0(\rho_0, \varepsilon_0) = T_1(\rho_1, \varepsilon_1) = \ldots = T_{N-1}(\rho_{N-1},\varepsilon_{N-1})

provides :math:`N-1` additional constraints. Finally, conservation of
energy provides one more constraint. In the end we have :math:`2 N`
constraints and :math:`2 N` unknowns.

In the code this is referred to as the ``PTESolverRhoU``.

The Density-Temperature Formulation
------------------------------------

Another choice is to treat the temperature as an independent
variable. Then the assumption of PTE implies that

.. math::

  T = T_0 = T_1 = \ldots = T_{N - 1}

which leads to a single additional unknown, the temperature
:math:`T`. The equality of pressure, now computed as

.. math::

  P_0(\rho_0, T) = P_1(\rho_1, T) = \ldots = P_{N-1}(\rho_{N-1}, T)

provides an additional :math:`N-1` constraints. Conservation of
energy, now computed as

.. math::

  u = \sum_{i=}^{N-1} \rho_i \varepsilon_i(\rho_i, T)

provides another constraint. This leads to :math:`N+1` constraints and
:math:`N+1` unknowns.

In the code this is referred to as the ``PTESolverRhoT``.

Using the Pressure-Temperature Equilibrium Solver
--------------------------------------------------

The PTE machinery is implemented in the
``singularity-es/closure/mixed_cell_models.hpp`` header. It is
entirely header only.

There are several moving parts. First, one must allocate scratch space
used by the solver. There are helper routines for providing the needed
scratch space, wich will tell you how many bytes per mixed cell are
required. For example:

.. cpp:function:: int PTESolverRhoTRequiredScratch(const int nmat);

and

.. cpp:function:: int PTESolverRhoURequiredScratch(const int nmat);

provide the number of real numbers (i.e., either ``float`` or
``double``) required for a single cell given a number of materials in
equilibriun for either the ``RhoT`` or ``RhoU`` solver. The equivalent
functions

.. cpp:function:: size_t PTESolverRhoTRequiredScratchInBytes(const int nmat);

and

.. cpp:function:: int PTESolverRhoURequiredScratchInBytes(const int nmat);

give the size in bytes needed to be allocated per cell given a number
of materials ``nmat``.

A solver in a given cell is initialized via a ``Solver`` object,
either ``PTESolverRhoT`` or ``PTESolverRhoU``. The constructor takes
the number of materials, some set of total quantities required for the
conservation constraints, and *indexer* objects for the equation of
state, the independent and dependent variables, and the ``lambda``
objects for each equation of state, similar to the vector API for a
given EOS. Here the indexers/vectors are not over cells, but
materials.

The constructor for the ``PTESolverRhoT`` is of the form

.. code-block:: cpp

  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PTESolverRhoT(const int nmat, EOS_t &&eos, const Real vfrac_tot, const Real sie_tot,
                Real_t &&rho, Real_t &&vfrac, Real_t &&sie, Real_t &&temp, Real_t &&press,
                Lambda_t &&lambda, Real *scratch);

where ``nmat`` is the number of materials, ``eos`` is an indexer over
equation of state objects, one per material, and ``vfrac_tot`` is a
number :math:`\in (0,1]` such that the sum over all volume fractions
adds up to ``vfrac_tot``. For a problem in which all materials
participate in PTE, ``vfrac_tot_`` should be 1. ``sie_tot`` is the
total specific internal energy in the problem, ``rho`` is an indexer
over densities, one per material. ``vfract`` is an indexer over volume
fractions, one per material. ``sie`` is an indexer over temperatures,
one per material. ``press`` is an indexer over pressures, one per
material. ``lambda`` is an indexer over lambda arrays, one ``Real *``
object per material. ``scratch`` is a pointer to pre-allocated scratch
memory, as described above. It is assumed enough scratch has been
allocated.

The constructor for the ``PTESolverRhoU`` has the same structure:

.. code-block:: cpp

  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PTESolverRhoU(const int nmat, const EOS_t &&eos, const Real vfrac_tot,
                const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch);

Both constructors are callable on host or device. In gerneral,
densities and internal energies are the required inputs. However, all
indexer quantities are asusmed to be input/output, as the PTE solver
may use unknowns, such as pressure and temperature, as initial guesses
and may reset input quantities, such as material densities, to be
thermodynamically consistent with the equilibrium solution.

Once a PTE solver has been constructed, one performs the solve with
the ``PTESolver`` function, which takes a ``PTESolver`` object as
input and returns a boolean status of either success or failure. For
example:

.. code-block:: cpp

  auto method = PTESolverRhoT<decltype(eos), decltype(rho), decltype(lambda)>(NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda, scratch);
  bool success = PTESolver(method);

For an example of the PTE solver machinery in use, see the
``test_pte.cpp`` file in the tests directory.
