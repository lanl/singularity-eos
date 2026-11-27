.. _using-closures:

Mixed Cell Closures
====================

In the single-material Euler equations, mass and energy are typically
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

For an example of how to use the PTE solvers, see
``examples/pte_2mat.cpp`` in the ``examples`` folder.

Governing Equations and Assumptions
------------------------------------

In a general sense the mixed
material closure rule takes the form

.. math::

  P_i = F(\rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1}) \\
  \rho_i = G(\rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1}) \\
  \epsilon_i = H(\rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1}) \\
  T_i = J(\rho, \epsilon, \mu_1, ..., \mu_i, ..., \mu_{N-1})

where each material has its own density, :math:`\rho_i`, and specific internal
energy, :math:`e_i`, as well as in principle its own pressure and temperature.
Each will be a function of the bulk density, :math:`\rho`, the bulk specific
internal energy, :math:`\epsilon`, and the mass fractions, :math:`\mu_j`,
of all other materials present.

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

The average density, :math:`\overline{\rho}_i`, (i.e. mass per *total* volume)
for a material in the total volume is

.. math::

  \overline{\rho}_i = \rho_i f_i,

where :math:`\rho_i` is the physical density (i.e. material mass per *material*
volume). It is important to note here that while the densities, :math:`\rho_i`,
and the volume fractions, :math:`f_i`, will vary as the closure model is
applied, the average densities, :math:`\overline{\rho}_i`, will all remain
constant, motiviating their internal use in the closure solvers. The total
density (mass of *participating* materials per total volume) is then

.. math::

  \rho = \sum_{i=0}^{N - 1} \overline{\rho}_i = \sum_{i=0}^{N-1} \rho_i f_i

Similarly the energy can be summed in a similar way so that

.. math::

  u = \rho \epsilon = \sum_{i = 0}^{N - 1} \rho_i \epsilon_i
  = \sum_{i = 0}^{N - 1} u_i

where :math:`u` is the total internal energy density (internal energy per unit
volume). Similarly, :math:`u_i` is analagous to :math:`\overline{\rho}_i` in
that it is the internal energy for a material averaged over the entire control
volume.

Internally, the closer models in ``singularity-eos`` operate on :math:`f_i`,
:math:`\overline{\rho}_i`, and :math:`u_i` as well as their total counterparts.
This is different than the forms stated at the beginning of this section so
that in essence the PTE solver has the form

.. math::

  P_i = F(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1}) \\
  \rho_i = G(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1}) \\
  \epsilon_i = H(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1}) \\
  T_i = J(\epsilon, f_\mathrm{tot}, f_1, ..., f_i, ..., f_{N-1},
          \rho_1, ..., \rho_i, ..., \rho_{N-1})


The important nuance here is that **the volume fractions and densities are both
inputs and outputs** in the current ``singularity-eos`` formulation of the
closure models. From physics perspective this can be confusing, but from a code
perspective this limits the number of variables that need to be passed to the
PTE solver and provides a convenient way to specify an initial guess for the
closure state.

.. warning::

  The volume fractions and material densities must be consistent so that

  .. math::

    \rho = \sum_{i=0}^{N - 1} \overline{\rho}_i = \sum_{i=0}^{N-1} \rho_i f_i


.. _massandvolumefractions:
.. note::

  Since mass fraction information is encoded in the specification of the
  component densities, :math:`\rho_i`, and component volume fractions,
  :math:`f_i`, they must be consistent so that

  .. math::

    \rho_i = \frac{\mu_i \rho}{f_i}

  .. math::

    \sum\limits_{i=0}^{N-1} f_i = f_\mathrm{tot}.

  For most practical applications, a previous PTE state for the current masses
  (i.e. from a previous timestep in a Lagrangian frame) or an appropriate
  prediction of the new PTE state (i.e. from advected values in an Eulerian
  frame) is usually available. This is usually the preferred input for the
  volume fractions and densities provided that they are consistent with the
  current mass fractions in the control volume.

  When a previous state is not available, an assuption must be made for how volume
  is partitioned between the materials. The simplest (but perhaps not the most
  effective) assumption is that volume is *equipartitioned* such that

  .. math::

    f_i = \frac{1}{N}.

  It is important to note though that this may not be sufficient in *many*
  cases. A better guess just use the mass fractions so that

  .. math::

    f_i = \mu_i = \frac{\overline{\rho}_i}{\rho},

  but this is really only valid when the materials have similar
  compressibilities. A further improvement could be made by weighting the mass
  fractions by the material bulk moduli to reflect the relative
  compressibilities.

Pressure-Temperature Equilibrium
--------------------------------

At present, ``singularity-eos`` focuses on several methods for finding a PTE
solution, i.e. one where the pressures and temperatures of the individual
materials are all the same. The methods presented differ in what they treat as
independent variables, and thus what precise system of equations they solve.
However they all share the above mathematical formulation.

In essence, the PTE equations can be posed as two residual equations:

.. math::

  f_\mathrm{tot} - \sum\limits_{i=0}^{N-1} f_i =
    \sum\limits_{i=0}^{N-1} f_i^*(x_i^*, y_i^*) - f_i(x_i, y_i)

.. math::

  u_\mathrm{tot} - \sum\limits_{i=0}^{N-1} u_i =
    \sum\limits_{i=0}^{N-1} u_i^*(x_i^*, y_i^*) - u_i(x_i, y_i)

where the superscript :math:`^*` denotes the variables at the PTE state,
:math:`f` corresponds to the volume fractions, and :math:`u` to the energy
density (see the previous section for more information). In these equations,
:math:`x` and :math:`y` represent some choice of independent thermodynamic
variables.

These are two non-linear residual equations that will need to be solved. In
``singularity-eos`` a Newton-Raphson method is used that first relies on
Taylor-expanding the equations about the equilibrium state in order to cast the
equations in terms of an update to the unknowns. The expansion about an
equilibrium state described by :math:`f_i^*(x_i, y_i)` and
:math:`u_i^*(x_i, y_i)` is

.. math::

  f_\mathrm{tot} - \sum\limits_{i=0}^{N-1} f_i(x_i, y_i) \approx
    \sum\limits_{i=0}^{N-1} (x_i^* - x_i)
      \left(\frac{\partial f_i}{\partial x_i}\right)_{y_i}
    + \sum\limits_{i=0}^{N-1} (y_i^* - y_i)
      \left(\frac{\partial f_i}{\partial y_i}\right)_{x_i}

.. math::

  u - \sum\limits_{i=0}^{N-1} u_i(x_i, y_i) \approx
    \sum\limits_{i=0}^{N-1} (x_i^* - x_i)
      \left(\frac{\partial u_i}{\partial x_i}\right)_{y_i}
    + \sum\limits_{i=0}^{N-1} (y_i^* - y_i)
      \left(\frac{\partial u_i}{\partial y_i}\right)_{x_i},

providing a means to update the guess for the equilbrium state. Minor
manipulations are needed to recast the derivatives in terms of accessible
thermodynamic derivatives, and then these equations can be written in matrix
form to solve for the unknown distance away from the equilibrum state. At each
iteration of the Newton-Raphson solver, the derivatives are recomputed and a
new update is found until some tolerance is reached. When a good initial guess
is used (such as a previous PTE state), some algorithms may converge in
relatively few iterations.

The choice of :math:`x` and :math:`y` is discussed below, but
crucially it determines the number of equations and unknowns needed to
specify the system.  For example, if pressure, :math:`P`, and
temperature, :math:`T`, are chosen, then the subscripts are eliminated
since we seek a solution where all materials have the same temperature
and pressure. (See :ref:`pressure-temperature-formulation`.) In this
formulation, there are two equations and two unkowns, and one such
solver is described below.

Most of the current PTE solvers in ``singularity-eos`` are cast in terms
of volume fraction and some other independent variable. Using material volume
fractions introduces :math:`N - 1` additional unknowns since all but one of the
volume fractions are independent from each other. The assumption of pressure
equilibrium naturally leads to the addition of :math:`N - 1` residual equations
of the form

.. math::

  P_i(f_i, y_i) - P_j(f_j, y_j) = 0,


These can also be written as a Taylor expansion about the equilibrium state such
that

.. math::

  P_i(f_i, y_j) - P_j(f_j, y_j)
    = (f^*_i - f_i) \left(\frac{\partial P_i}{\partial f_i}\right)_{y_i}
    + (y^*_i - y_i) \left(\frac{\partial P_i}{\partial y_i}\right)_{f_i} \\
    - (f^*_j - f_j) \left(\frac{\partial P_j}{\partial f_j}\right)_{y_j}
    - (y^*_j - y_j) \left(\frac{\partial P_j}{\partial y_j}\right)_{f_j},

where the equations are typically written such that :math:`j = i + 1`. Since the
equlibrium pressure is the same for both materials, the term cancels out and
the material pressures are left.


Formulating the closure equations in terms of volume fractions instead of
densities has the benefit of allowing the volume constraint to be written in
terms of just the volume fractions:

.. math::

  f_\mathrm{tot} - \sum\limits_{i=0}^{N-1} f_i =
    \sum\limits_{i=0}^{N-1} (f_i^* - f_i).

The EOS only returns derivatives in terms of density though, so a the density
derivatives must be transformed to volume fraction derivatives via

.. math::

  \left(\frac{\partial Q}{\partial f_i}\right)_X
    = - \frac{\rho_i^2}{\rho}\left(\frac{\partial Q}{\partial \rho_i}\right)_X,

were :math:`Q` and :math:`X` are arbitrary thermodynamic variables. At this
point, there are :math:`N + 1` equations and unknowns in the PTE sover. The
choice of the second independent variable is discussed below and has
implications for both the number of additional unknowns and the stability of the
method.

.. _pressure-temperature-formulation:
The Pressure-Temperature Formulation
`````````````````````````````````````

An obvious choice is to treat the independent variables as pressure
and temperature. Then one has only two equations and two unknowns. The
residual contains only the volume fraction and energy summmation rules
described above. Taylor expanding these residuals about fixed
temeprature and pressure points leads to two residual equations of the
form

.. math::

  1 - \sum_{i=0}^{N-1} f_i = (T^* - T) \sum_{i = 0}^{N-1} \left(\frac{\partial f_i}{\partial T}\right)_P + (P^* - P) \sum_{i = 0}^{N-1} \left(\frac{\partial f_i}{\partial P}\right)_T\\
  u_{tot} - \sum_{i=0}^{N-1} u_i = (T^* - T) \sum_{i = 0}^{N-1} \left(\frac{\partial u_i}{\partial T}\right)_P + (P^* - P) \sum_{i = 0}^{N-1} \left(\frac{\partial u_i}{\partial P}\right)_T

However, derivatives in the volume fraction are not easily
accessible. To access them, we leverage the fact that

.. math::

  \bar{\rho}_i = \rho_i f_i,

and thus

.. math::

  d f_i = - \frac{\overline{\rho}}{\rho_i^2} d \rho_i.

Thus the residual can be recast as

.. math::

  f_\mathrm{tot} - \sum_{i=0}^{N-1} = -(T^* - T) \sum_{i = 0}^{N-1} \frac{\bar{\rho}_i}{\rho_i^2} \left(\frac{\partial \rho_i}{\partial T}\right)_P - (P^* - P) \sum_{i = 0}^{N-1} \frac{\bar{\rho}_i}{\rho_i^2} \left(\frac{\partial \rho_i}{\partial P}\right)_T\\
  u_\mathrm{tot} - u_i = (T^* - T) \sum_{i = 0}^{N-1} \left(\frac{\partial u_i}{\partial T}\right)_P + (P^* - P) \sum_{i = 0}^{N-1} \left(\frac{\partial u_i}{\partial P}\right)_T

where :math:`\rho_{\mathrm{tot}}` is the sum of densities over all
materials. These residual equations can then be cast as a matrix
equation to solve for pressure and temperature.

The primary advantage of the pressure-temperature space solver is that
it has only two independent variables and two unknowns, meaning the
cost scales only linearly with the number of materials, not
quadratically (or worse). The primary disadvantage, is that most
equations of state are not formulated in terms of pressure and
temperature, meaning additional inversions are required. In the case
where a root-find is required for this inversion, performance may
suffer for a small number of materials compared to a different
formulation.

In the code, this method is referred to as ``PTESolverPT``.

.. _density-energy-formalism:
The Density-Energy Formulation
```````````````````````````````

One choice is to treat volume fractions and material energies as independent
quantities, but the material energies provide :math:`N - 1` additional
unknowns. The additional degrees of freedom are satisfied by requiring that the
material temperatures be equal. As a result, we add :math:`N - 1` residual
equations of the form

.. math::

  T_i(\rho_i, \epsilon_i) - T_j(\rho_j, \epsilon_j) = 0.

Again Taylor expanding about the equilibirum state, this results in a set of
equations of the form

.. math::

  T_i(f_i, \epsilon_i) - T_j(f_j, \epsilon_j)
    = (f^*_i - f_i) \left(\frac{\partial T_i}{\partial f_i}\right)_{\epsilon_i}
    + (\epsilon^*_i - \epsilon_i) \
        \left(\frac{\partial T_i}{\partial \epsilon_i}\right)_{f_i} \\
    - (f^*_j - f_j) \left(\frac{\partial T_j}{\partial f_j}\right)_{\epsilon_j}
    - (\epsilon^*_j - \epsilon_j)
        \left(\frac{\partial T_j}{\partial \epsilon_j}\right)_{f_j}

Here there are a total number of :math:`2N` equations and unknowns, which
results in a fairly large matrix to invert when many materials are present in a
cell. Further, the density-energy derivatives may require inversion of any EOS
with density and temperature as the natural variables. In the case of tabular
EOS, an iterative inversion step may be required to find the density-energy
state by iterating on temperature; there may also be a loss of accuracy in the
derivatives depending on how they are calculated.

In general, the density-temperature formulation of the PTE solver seems to be
more stable and performant and is usually preferrred to this formulation.

In the code this is referred to as the ``PTESolverRhoU``.

The Density-Temperature Formulation
````````````````````````````````````

Another choice is to treat the temperature as an independent variable, requiring
no additional equations. The energy residual equation then takes the form

.. math::

  u - \sum\limits_{i=0}^{N-1} u_i(f_i, T) \approx
    \sum\limits_{i=0}^{N-1} (f_i^* - f_i)
      \left(\frac{\partial u_i}{\partial f_i}\right)_{T}
    + (T^* - T)\sum\limits_{i=0}^{N-1}
      \left(\frac{\partial u_i}{\partial T}\right)_{f_i},

where the temperature difference can be factored out of the sum since it doesn't
depend on material index.

In the code this is referred to as the ``PTESolverRhoT``.

Fixed Pressure or Temperature
"""""""""""""""""""""""""""""

For initialization, the energy of a mixed material region is usually unknown
while the density, mass fractions, and either temperature or pressure *are*
known. To find the energy, a PTE solve is required, but with the added
constraint of the fixed pressure or temperature.

Fixed temperature
^^^^^^^^^^^^^^^^^

When the temperature and total density are known, the equilibrium pressure and
the component densities are unknown. This requires a total of :math:`N`
equations and unknowns. Since the total energy is unknown, it can be dropped
from the contraints leaving just the :math:`N - 1` pressure equality equations
and the volume fraction sum constraint. The pressure residuals can then be
simplified to be

.. math::

  P_i(f_i, T) - P_j(f_j, T)
    = (f^*_i - f_i) \left(\frac{\partial P_i}{\partial f_i}\right)_{T}
    - (f^*_j - f_j) \left(\frac{\partial P_j}{\partial f_j}\right)_{T}

In the code this is referred to as the ``PTESolverFixedT``.

Fixed pressure
^^^^^^^^^^^^^^

When the pressure and total density are known, the procedure is slightly more
complicated. Since the pressure is known but the independent variables are
density and temperature, there are :math:`N + 1` unknowns for the component
densities and the unknown equilibrium temperature.

Once again, the energy constraint is dropped since the energy is unknown, but
since the equilibrium pressure is a *specified* quantity, the pressure residual
equations must be modified to take the form

.. math::

  P_i^*(f^*_i, T) - P_i(f_i, T)
    = (f^*_i - f_i) \left(\frac{\partial P_i}{\partial f_i}\right)_{T}
    - (T^* - T) \left(\frac{\partial P_i}{\partial T}\right)_{f_i}.

Note that this results in :math:`N` equations for each of the individual
material pressures.

In the code this is referred to as the ``PTESolverFixedP``.

Using the Pressure-Temperature Equilibrium Solver
```````````````````````````````````````````````````

The PTE machinery is implemented in the
``singularity-es/closure/mixed_cell_models.hpp`` header. It is
entirely header only.

There are several moving parts. First, one must allocate scratch space
used by the solver. There are helper routines for providing the needed
scratch space, wich will tell you how many bytes per mixed cell are
required. For example:

.. cpp:function:: int PTESolverRhoTRequiredScratch(const size_t nmat);

and

.. cpp:function:: int PTESolverRhoURequiredScratch(const size_t nmat);

provide the number of real numbers (i.e., either ``float`` or
``double``) required for a single cell given a number of materials in
equilibriun for either the ``RhoT`` or ``RhoU`` solver.

The equivalent functions

.. cpp:function:: size_t PTESolverRhoTRequiredScratchInBytes(const size_t nmat);

and

.. cpp:function:: size_t PTESolverRhoURequiredScratchInBytes(const size_t nmat);

give the size in bytes needed to be allocated per cell given a number
of materials ``nmat``. Alternatively, there are a static member functions
for each closure model that provides the same information:

.. cpp:function:: int RequiredScratch(const size_t nmat);

.. cpp:function:: size_t RequiredScratchInBytes(const size_t nmat);

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
  and outputs**. They are used to define the internal :math:`\overline
  {\rho}_i` variables at the beginning of the PTE solve. The volume fractions
  and densities at the end of the PTE solve will represent those for the new
  PTE state. It's important to note that :math:`\overline{\rho}_i` remain
  constant throughout the calculation.

.. warning::

  The PTE solvers **require** that all input densities and volume fractions
  are non-zero. As a result, ``nmat`` refers to the number of *participating*
  materials. The user is encouraged to wrap their data arrays using an
  ``Indexer`` concept where, for example, three paricipating PTE materials
  might be indexed as 5, 7, 20 in the material arrays. This requires overloading
  the square bracket operator to map from PTE idex to material index.

The constructor for the ``PTESolverRhoT`` is of the form

.. code-block:: cpp

  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverRhoT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot, const Real sie_tot,
                Real_t &&rho, Real_t &&vfrac, Real_t &&sie, Real_t &&temp, Real_t &&press,
                Lambda_t &&lambda, Real *scratch, const Real Tnorm = 0.0,
                const MixParams &params = MixParams())

where ``nmat`` is the number of materials, ``eos`` is an indexer over
equation of state objects, one per material, and ``vfrac_tot`` is a
number :math:`\in (0,1]` such that the sum over all volume fractions
adds up to ``vfrac_tot``. For a problem in which all materials
participate in PTE, ``vfrac_tot_`` should be 1. ``sie_tot`` is the
total specific internal energy in the problem, ``rho`` is an indexer
over densities, one per material. ``vfract`` is an indexer over volume
fractions, one per material. ``sie`` is an indexer over temperatures,
one per material. ``press`` is an indexer over pressures, one per
material. ``lambda`` is an indexer over lambda arrays, one per
material. ``scratch`` is a pointer to pre-allocated scratch memory, as
described above. It is assumed enough scratch has been allocated.
The optional argument ``Tnorm`` allows for host codes to pass in a
normalization for the temperature scale. Initial guesses for density
and temperature may be passed in through the ``rho`` and ``temp`` input
parameters.

The optional ``MixParams`` input contains a struct of runtime
parameters that may be used by the various PTE solvers. This struct
contains the following member fields, with default values:

.. code-block:: cpp

  struct MixParams {
    bool verbose = false; // verbosity
    Real derivative_eps = 3.0e-6;
    Real pte_rel_tolerance_p = 1.e-6;
    Real pte_rel_tolerance_e = 1.e-6;
    Real pte_rel_tolerance_t = 1.e-4;
    Real pte_abs_tolerance_p = 0.0;
    Real pte_abs_tolerance_e = 1.e-4;
    Real pte_abs_tolerance_t = 0.0;
    Real pte_residual_tolerance = 1.e-8;
    std::size_t pte_max_iter_per_mat = 128;
    Real line_search_alpha = 1.e-2;
    std::size_t line_search_max_iter = 6;
    Real line_search_fac = 0.5;
    Real vfrac_safety_fac = 0.95;
    Real temperature_limit = 1.0e15;
    Real default_tguess = 300.;
    Real min_dtde = 1.0e-16;
  };

where here ``verbose`` enables verbose output in the PTE solve is,
``derivative_eps`` is the spacing used for finite differences
evaluations of equations of state when building a jacobian. The
``pte_rel_tolerance_p``, ``pte_rel_tolerance_e``, and
``pte_rel_tolerance_t`` variables are relative tolerances for the
error in the pressure, energy, temperature respectively. The
``pte_abs_tolerance_*`` variables are the same but are absolute,
rather than relative tolerances. ``pte_residual_tolerance`` is the
absolute tolerance for the residual.

The maximum number of iterations the solver is allowed to take before
giving up is ``pte_max_iter_per_mat`` multiplied by the number of
materials used. ``line_search_alpha`` is used as a safety factor in
the line search. ``line_search_max_iter`` is the maximum number of
iterations the solver is allowed to take in the line
search. ``line_search_fac`` is the step size in the line
search. ``vfrac_safety_fac`` limites the relative amount the volume
fraction can take in a given iteration. ``temperature_limit`` is the
maximum temperature allowed by the solver. ``default_tguess`` is used
as an initial guess for temperature if a better guess is not passed in
or cannot be inferred. ``min_dtde`` is the minmum that temperature is
allowed to change with respect to energy when computing Jacobians.

.. note::

  If ``MixParams`` are not provided, the default values are used. Not
  all ``MixParams`` are used by every solver.

The constructor for the ``PTESolverRhoU`` has the same structure:

.. code-block:: cpp

  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PTESolverRhoU(const int nmat, const EOS_t &&eos, const Real vfrac_tot,
                const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
                const Real Tnorm = 0, const MixParams &params = MixParams());

Both constructors are callable on host or device. In gerneral,
densities and internal energies are the required inputs. However, all
indexer quantities are asusmed to be input/output, as the PTE solver
may use unknowns, such as pressure and temperature, as initial guesses
and may reset input quantities, such as material densities, to be
thermodynamically consistent with the equilibrium solution.

Once a PTE solver has been constructed, one performs the solve with
the ``PTESolver`` function, which takes a ``PTESolver`` object as
input and returns a ``SolverStatus`` struct:

.. code-block:: cpp

  auto method = PTESolverRhoT<decltype(eos), decltype(rho), decltype(lambda)>(NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda, scratch);
  auto status = PTESolver(method);

The status struct is of the form:

.. code-block:: cpp

  struct SolverStatus {
    bool converged;
    std::size_t max_niter;
    std::size_t max_line_niter;
    Real residual;
  };

where ``converged`` will report whether or not the solver successfully
converged, ``residual`` will report the final value of the residual,
``max_niter`` will report the total number of iterations that the
solver performed and ``max_line_niter`` will report the maximum number
of iterations within a line search that the solver performed. For an
example of the PTE solver machinery in use, see the ``test_pte.cpp``
file in the tests directory.

Initial Guesses for PTE Solvers
`````````````````````````````````

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
material density, ``vfrac`` is an indexer over material volume
fractions, and the optional argument ``Tguess`` allows for callers to
pass in an initial guess that could accelerate finding a solution.
This function does a 1-D root find to find the temperature at which
the material internal energies sum to the total.  The root find does
not have a tight tolerance -- instead the hard-coded tolerance was
selected to balance performance with the accuracy desired for an
initial guess in a PTE solve.  If a previous temperature value is
unavailable or some other process may have significantly modified the
temperature since it was last updated, this function can be quite
effective.
