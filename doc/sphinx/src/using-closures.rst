.. _using-closures:

Mixed Cell Closures
====================

In the single-material Euler equations, the mass and energy are typically
evolved and the EOS is called to provide a pressure for the momentum equations.
When transitioning to a multi-material approach, a single velocity is typically
used and the Euler equations are solved with respect to the bulk fluid motion.
In this case, the pressure contribution to momentum isn't well-defined and in
principle each material could have its own pressure contribution to material
motion. Furthermore, the paritioning of volume and energy between the materials
in the flow is not well-defined either.

As a result, a multi-material closure rule is needed to determine both how to
compute the pressure response of the flow and how to partition the volume and
energy between the individual materials. In this situation, one must decide how
to compute thermodynamic quantities in that cell for each material.

Governing Equations and Assumptions
------------------------------------

In a general sense then the mixed
material closure rule takes the form

.. math::

  P_i = F(rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1}) \\
  \rho_i = G(rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1}) \\
  \epsilon_i = H(rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1}) \\
  T_i = J(rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1})

where each material has its own density, :math:`\rho_i`, and specific internal
energy, :math:`e_i`, as well as in principle its own pressure and temperature.
Each will be a function of the bulk density, :math:`\rho`, the bulk specific
internal energy, :math:`\epsilon`, and both the mass fraction, :math:`\mu_j`,
of all materials present.

For convenience of minimizing input arguments though, the solvers in
``singularity-eos`` are typically posed in a different way that can sometimes be
unintuitive.

For some finite total volume :math:`V`, each material occupies some fraction of
that volume given by the *volume fraction*
:math:`f_i` such that

.. math::

  \sum_{i=0}^{N - 1} f_i = f_\mathrm{tot},

where :math:`f_\mathrm{tot}` is the total fraction of the total volume being
considered (in principle, different closure models can be used for different
sets of materials). To consider the entire volume, :math:`f_\mathrm{tot}` can
simply be set to one.

The average density, :math:`\bar{\rho}_i`, (i.e. mass per *total* volume) for a
material in the total volume is

.. math::

  \bar{\rho}_i = \rho_i f_i,

where :math:`\rho_i` is the physical density (i.e. material mass per *material*
volume). The total density (mass of *participating* materials per total volume)
is then

.. math::

  \rho = \sum_{i=0}^{N - 1} \bar{\rho}_i = \sum_{i=0}^{N-1} \rho_i f_i

Similarly the energy can be summed in a similar way so that

.. math::

  u = \rho \epsilon = \sum_{i = 0}^{N - 1} \rho_i \epsilon_i 
  = \sum_{i = 0}^{N - 1} u_i

where :math:`u` is the total internal energy density (internal energy per unit
volume). Similarly, :math:`u_i` is analagous to :math:`\bar{\rho}_i` in that it
is the internal energy for a material averaged over the entire control volume.

Internally, the closer models in ``singularity-eos`` operate on :math:`f_i`,
:math:`\bar{\rho}_i`, and :math:`u_i` as well as their total counterparts. This
is different than the forms stated at the beginning of this section so that in
essence the PTE solver has the form

.. math::

  P_i = F(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1}) \\
  \rho_i = G(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1}) \\
  \epsilon_i = H(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1}) \\
  T_i = J(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1})


The important nuance here is the **the volume fractions are both inputs and
outputs** in the current ``singularity-eos`` formulation of the closure models.
From physics perspective this can be confusing, but from a code perspective this
limits the number of variables that need to be passed to the PTE solver and
provides a convenient way to specify an initial guess for the closure state.

.. note::

  The mass fraction information is encoded in the specifciation of :math:`f_i`
  and :math:`\rho_i`. In order to convert to component densities,
  :math:`\rho_i` from mass fractions, :math:`\mu_i` and total density,
  :math:`\rho`, volume fractions must be ***assumed*** in some consistent way
  so that

  .. math::

    \rho_i = \frac{\mu_i \rho}{f_i}

  .. math::

    \sum\limits_{i=0}^{N-1} f_i = f_\mathrm{tot}.

  One such assumption would be to **equipartition** the volume between all
  materials such that

  .. math::

    f_i = \frac{1}{N}

Pressure-Temperature Equilibirum
================================

At present, ``singularity-eos`` focuses on several methods for finding a PTE
solution, i.e. one where the pressures and temperatures of the individual
materials are all the same. The methods presented differ in what they treat as
independent variables, and thus what precise system of equations they solve.
However they all share the above mathematicaly formulation.

In essence, the PTE equations can be posed as two residual equations:

.. math::

  f_\mathrm{tot} - \sum\limits_{i=0}^{N-1} f_i = 
    \sum\limits_{i=0}^{N-1} f_i^*(x_i*, y_i*) - f_i(x_i, y_i)

.. math::

  u_\mathrm{tot} - \sum\limits_{i=0}^{N-1} u_i = 
    \sum\limits_{i=0}^{N-1} u_i^*(x_i*, y_i*) - u_i(x_i, y_i)

where the superscript :math:`^*` denotes the variables at the PTE state,
:math:`f` corresponds to the volume fractions, and :math:`u` to the energy
density (see the previous section for more information). In these equations,
:math:`x` and :math:`y` represent some choice of independent thermodynamic
variables.

Then the energy and volume fraction constraint equations can be Taylor-expanded
about the equilibrium states :math:`f_i^*(\rho_i, y_i)` and
:math:`u_i^*(\rho_i, y_i)` so that they become

.. math::

  f_\mathrm{tot} - \sum\limits_{i=0}^{N-1} f_i(\rho_i, y_i) \approx
    \sum\limits_{i=0}^{N-1} (\rho_i^* - \rho_i)
      \left(\frac{\partial f_i}{\partial \rho_i}\right)_{y_i}
    + \sum\limits_{i=0}^{N-1} (y_i^* - y_i)
      \left(\frac{\partial f_i}{\partial y_i}\right)_{\rho_i}

.. math::

  u - \sum\limits_{i=0}^{N-1} u_i(\rho_i, y_i) \approx
    \sum\limits_{i=0}^{N-1} (\rho_i^* - \rho_i)
      \left(\frac{\partial u_i}{\partial \rho_i}\right)_{y_i}
    + \sum\limits_{i=0}^{N-1} (y_i^* - y_i)
      \left(\frac{\partial u_i}{\partial y_i}\right)_{\rho_i},

removing the dependence on the unknown equilibrium state. Minor manipulations
are needed to recast the derivatives in terms of accessible thermodynamic
derivatives and then these equations can be written in matrix form to solve for
the unknown distance away from the equilibrum state.

The choice of :math:`x` and :math:`y` is discussed below, but crucially it
determines the number of equations and unknowns needed to specify the system.
For example, if pressure, :math:`P`, and temperature, :math:`T`, are chosen,
then the subscripts are eliminated since we seek a solution where all materials
have the same temperature and pressure. In this formulation, there are two
equations and two unkowns, but due to the difficulty of inverting an
equation of state to be a function of pressure and temperature,
``singularity-eos`` does not have any PTE solvers that are designed to use
pressure and temperature as independent variables.

Instead, all of the current PTE solvers in ``singularity-eos`` are cast in terms
of density and another independent variable. This introduces
:math:`N - 1` additional unknowns since the each material density is independent
except for the last. The assumption of pressure equilibrium naturally leads to
an addition :math:`N - 1` residual equations of the form

.. math::

  P_i(\rho_i, y_j) - P_j(\rho_j, y_j) = 0

The choice of the second independent variable is discussed below:

The Density-Energy Formulation
---------------------------------

One choice is to treat volume fractions and material energies as independent
quantities. However, the material energies provide :math:`N - 1` additional
unknowns. This requires that euqality of material temperatures satisfy the
additional degrees of freedom. As a result, we add :math:`N - 1` residual
equations of the form

.. math::

  T_i(\rho_i, \epsilon_j) - T_j(\rho_j, \epsilon_j) = 0.

In the code this is referred to as the ``PTESolverRhoU``.

The Density-Temperature Formulation
------------------------------------

Another choice is to treat the temperature as an independent
variable. Then the assumption of temperature equilibrium requires no additional
equations, and the energy and volume fraction constraints take the form

.. math::

  f_\mathrm{tot} - \sum\limits_{i=0}^{N-1} f_i(\rho_i, y_i) \approx
    \sum\limits_{i=0}^{N-1} (\rho_i^* - \rho_i)
      \left(\frac{\partial f_i}{\partial \rho_i}\right)_{y_i}
    + (T^* - T)\sum\limits_{i=0}^{N-1}
      \left(\frac{\partial f_i}{\partial T}\right)_{\rho_i}

.. math::

  u - \sum\limits_{i=0}^{N-1} u_i(\rho_i, y_i) \approx
    \sum\limits_{i=0}^{N-1} (\rho_i^* - \rho_i)
      \left(\frac{\partial u_i}{\partial \rho_i}\right)_{y_i}
    + (T^* - T)\sum\limits_{i=0}^{N-1}
      \left(\frac{\partial u_i}{\partial T}\right)_{\rho_i},

where the temperature difference can be factored out of the sum since it doesn't
depend on material index.

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

.. warning::

  It bears repeating: **both the volume fractions and densities act as inputs
  and outputs**. They are used to define the internal :math:`\bar
  {\rho}_i` variables at the beginning of the PTE solve. The volume fractions
  and densities at the end of the PTE solve will represent those for the new
  PTE state.

.. warning::

  The PTE solvers ***require*** that all input densities and volume fractions
  are non-zero. As a result, ``nmat`` refers to the number of *participating*
  materials. The user is encouraged to wrap their data arrays using an
  ``Indexer`` concept where, for example, three paricipating PTE materials
  might be indexed as 5, 7, 20 in the material arrays. This requires overloading
  the square bracket operator to map from PTE idex to material index.

The constructor for the ``PTESolverRhoT`` is of the form

.. code-block:: cpp

  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PTESolverRhoT(const int nmat, EOS_t &&eos, const Real vfrac_tot, const Real sie_tot,
                Real_t &&rho, Real_t &&vfrac, Real_t &&sie, Real_t &&temp, Real_t &&press,
                Lambda_t &&lambda, Real *scratch, const Real Tguess = 0);

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
allocated.  Finally, the optional argument ``Tguess`` allows for host
codes to pass in an initial temperature guess for the solver.  For more
information on initial guesses, see the section below.

The constructor for the ``PTESolverRhoU`` has the same structure:

.. code-block:: cpp

  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PTESolverRhoU(const int nmat, const EOS_t &&eos, const Real vfrac_tot,
                const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
                const Real Tguess = 0);

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

Initial Guesses for PTE Solvers
------------------------------------

As is always the case when solving systems of nonlinear equations, good initial
guesses are important to ensure rapid convergence to the solution.  For the PTE
solvers, this means providing intial guesses for the material densities and the
equilibrium temperature.  For material densities, a good initial guess is often
the previous value obtained from a prior call to the solver. ``singularity-eos``
does not provide any mechanism to cache these values from call to call, so it is
up to the host code to provide these as input to the solvers.  Note that the
input values for the material densities and volume fractions are assumed to be
consistent with the conserved cell-averaged material densities, or in other
words, the produce of the input material densities, volume fractions, and cell
volume should equal the amount of mass of each material in the cell.  This
consistency should be ensured for the input values or else the solvers will not
provide correct answers.

For the temperature initial guess, one can similarly use a previous value for
the cell.  Alternatively, ``singularity-eos`` provides a function that can be
used to provide an initial guess.  This function takes the form

.. code-block:: cpp

  template <typename EOSIndexer, typename RealIndexer>
  PORTABLE_INLINE_FUNCTION Real ApproxTemperatureFromRhoMatU(
    const int nmat, EOSIndexer &&eos, const Real u_tot, RealIndexer &&rho,
    RealIndexer &&vfrac, const Real Tguess = 0.0);

where ``nmat`` is the number of materials, ``eos`` is an indexer over
equation of state objects, ``u_tot`` is the total material internal
energy density (energy per unit volume), ``rho`` is an indexer over
material density, ``vfrac`` is an indexer over material volume fractions,
and the optional argument ``Tguess`` allows for callers to pass in a guess
that could accelerate finding a solution.  This function does a 1-D root find
to find the temperature at which the material internal energies sum to the
total.  The root find does not have a tight tolerance -- instead the
hard-coded tolerance was selected to balance performance with the accuracy
desired for an initial guess in a PTE solve.  If a previous temperature value
is unavailable or some other process may have significantly modified the
temperature since it was last updated, this function can be quite effective.
