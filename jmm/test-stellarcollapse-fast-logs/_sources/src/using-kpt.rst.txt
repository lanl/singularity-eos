.. _using-kpt:

.. _WillsKPT: https://www.osti.gov/biblio/2202604

.. _WillsParameterstudy: https://www.osti.gov/biblio/1900445

Kinetic Phase Transition Framework
==================================

A step in adapting to more realistic calculations is to 
get away from that phase transitions are instantanious and the 
system is in, at least, local thermal equilibrium (LTE) all the time.
Allowing for superheated (the material is hotter than the equilibrium state) 
and supercooled (the material is cooler than the equilibrium state) materials is
a large step towards more realistic modelling. Gas bubbles "exploding" when a 
spoon is put into a microwave heated cup of coffee is an example of superheating and
raindrops turning directly to ice when hitting the asphalt is an example
of supercooling.

Phase transitions are governed by the Gibbs free energy, :math:`G(P,T)`. The equilibrium phase is the phase with lowest
Gibbs free energy. If the material is heated or cooled from this equilibrium state and another
phase emerges to have a lower Gibbs free energy, a phase transition occurs and the system eventually transform
to this new phase. However, all phase transitions need some type of *nucleation event* to get started.
For the liquid raindrop cooling down while falling to the ground through very cold air, the hit against the
asphalt is the nucleation event, and the supercooled liquid in the drop instantaniously transform to the equilibrium ice phase. 
For the very hot coffee, the spoon insertion is the "nucleation" event when the very hot liquid coffee 
instantaneously transforms into its equilibrium phase, gas. 

In a solid state phase transition, say between an face-centered cubic (fcc) phase and a body-centered
cubic (bcc) phase, the mere rearrangement of atoms can take some time, also resulting in a phase transition
that doesn't go directly from the previously lowest Gibbs free energy phase to the presently 
lowest Gibbs free energy phase. 

``singularity-eos`` provides the ingredients to build a kinetic phase transition framework in a host code.

A comprehensive guide for implementing KPT into a hydro code is available `here <WillsKPT_>`_.

Mass and volume fractions
-------------------------

The first ingredient needed in a KPT framework is a quantity describing what phase a material currently
is in. *Mass fractions*, :math:`\mu_i`, gives a measure of how much of the total mass of the material in a cell, :math:`M`, 
is in a specific phase

.. math::

 \mu_i = \frac{M_i}{M} \qquad \qquad \sum_i \mu_i = 1

For completeness we also introduce *volume fractions*, :math:`f_i`, that give how much of the volume of 
the material in a cell, :math:`V`, is in a specific phase

.. math::

 f_i = \frac{V_i}{V} \qquad \qquad \sum_i f_i = 1

With these definitions it is clear that both :math:`0 \leq \mu_i \leq 1` and :math:`0 \leq f_i \leq 1`.
Mass and volume fractions are discussed more closely in the :ref:`mixed cell closure <massandvolumefractions>` section.
Note, however, that mass and volume fractions in the KPT framework give how much of each *phase* of *one* material is present,
while in mixed cells they give how much of one *material* is present. Even though we can use some of the same routines
for mixed material cells and KPT, implementing KPT into a mixed cell environment need to use separate mass and volume 
fractions for the mixed cell and within each material in the cell.

Homogenization
--------------

A cell in a simulation can contain several different materials in several different phases. In order to have a well defined state
in this heterogeneous cell we need to have some rules, preferably based in physics, on how to *homogenize* the cell into one average
or mixture "material" with a well defined state. The natural and most common homogenization method for phase transitions is 
*pressure-temperature equilibrium* (PTE), which also is a common method for mixed materials cells.

PTE
'''

It is natural to use PTE for phase transitions since the Gibbs free energy, :math:`G(P,T)`, is uniquely defined if :math:`P` and :math:`T`
are given. PTE means that every phase/material present in a cell has the same pressure and temperature and 
given the specific internal energy, :math:`E`, specific volume, :math:`V`, the set of mass fractions :math:`\{\mu_i\}`, and one Equations of state (EOS) for each phase :math:`i`, we have:

.. math::

 V &=& \sum_i \mu_i V_i \\
 E &=& \sum_i \mu_i E_i \\
 P &=& P_i(V_i,E_i) \\
 T &=& T_i(V_i,E_i)

where :math:`i` denotes the different phases of one material. For :math:`N` phases this gives enough of equations (:math:`2N`) to determine all
:math:`V_i` and :math:`E_i`, and thus :math:`P` and :math:`T` and Gibbs Free energy :math:`G_i(P,T)` for each phase: 

.. math::

 G_i &=& E_i - P V_i - T S_i \\
 G &=& \sum_i \mu_i G_i

Note that we need information about the specific entropy, :math:`S`, for each phase. This means that each phase EOS needs to be a :ref:`complete EOS <Complete eos>` with a full description of the states.

We recommend using the :ref:`density-energy formalism <density-energy-formalism>` (``PTESolverRhoU``) described in the :ref:`mixed cell closure <using-closures>` section, as PTE solver in the KPT framework.

Mass fraction update
--------------------

In a computational cell at a time :math:`t_0`, the set of current mass fractions, :math:`\{\mu_i\}_{t_0}`, the current specific internal energy, :math:`E_{t_0}`, and specific volume, :math:`V_{t_0}`, 
give all the information needed, including the state of the phases, about the current state of the cell via a PTE solve. This state is the ground for advancing the state in the cell to a 
state at time :math:`t = t_0 + dt` via mass, momentum, and energy conservation in a hydro code. If a system is considered to be in equilibrium also through a phase transition, this state also gives
the new mass fractions. However, in order to allow kinetics to influence the phase transition a new material specific model is needed, a mass fraction update model.

Equilibrium phase transitions
'''''''''''''''''''''''''''''

Tables of equilibrium mass fractions as a function of :math:`P` and :math:`T`, can be constructed at the time of constructing complete EOSs for the different phases of a system. 
Phases participation in an equilibrium phase transition all have the same Gibbs free energy even though their specific internal energy, specific volume, and specific entropy, are different
(see equation for Gibbs free energy above). Equilibrium phase transition mass fraction tables will be made available through ``singularity-eos`` in the near future.

Kinetic phase transition models
'''''''''''''''''''''''''''''''

We write the new mass fractions at time :math:`t` as

.. math::

 \mu_i^{t} = \mu_i^{t_0} + \frac{d\mu_i}{dt} dt  

which with discretized time becomes

.. math::

 \mu_i^{t} = \mu_i^{t_0} + \dot{\mu_i} \Delta t

Using mass conservation (:math:`\sum_i \mu_i = 1`), we see that 

.. math::

 0 = \sum_i \dot{\mu_i} = \sum_i \sum_j ( \mu_j R_{ji} - \mu_i R_{ij} )
   
where :math:`R_{ij}` is the mass transportation rate from phase :math:`i` to phase :math:`j`, and we can derive the master equation

.. math::

 \dot{\mu_i} = \sum_j ( \mu_j R_{ji} - \mu_i R_{ij} )

This equation simply states the fact that all mass transforming from one phase ends up in another phase, it is
just mass conservation.

A KPT model is a model for the :math:`R_{ij}`.

Carl Greeff's KPT model
'''''''''''''''''''''''

Carl Greeff formulated an empirical mass fraction update model as

.. math::

 R_{ij} = \nu_{ij} \theta(G_i-G_j) \frac{G_i-G_j}{B_{ij}} \exp\left[ \left( \frac{G_i-G_j}{B_{ij}}\right)^2 \right]

where :math:`\nu_{ij}` and :math:`B_{ij}` are material dependent fitting constants, and :math:`\theta` is the Heaviside step function that is
:math:`0` for a negative argument and :math:`1` for a positive argument. 
Mattsson-Wills performed a `parameter study <WillsParameterstudy_>`_ of this model, 
which can be used as a guide for how to choose these parameters for a specific material and a specific phase transition.

This model is included in ``singularity-eos`` and its signature is

.. code-block:: cpp

  LogRatesCGModel(const Real *w, const Real *b, const int num_phases, const Real *gibbs,
                  const int *gibbsorder, Real *logRij, int *fromto)

where ``w`` is :math:`\nu_{ij}`, ``b`` is :math:`B_{ij}`, ``num_phases`` is
:math:`N`, the number of phases, and ``gibbs`` is :math:`G_i`. 
``gibbsorder`` is an array of length :math:`N`, where the ``gibbs`` phase indices are ordered 
from the largest Gibbs free energy phase to the lowest Gibbs free energy phase (see figure below),
``logRij`` is :math:`\log(R_{ij})`, and ``fromto`` is a map between the phase indices :math:`i` and :math:`j` and the :math:`ij` phase transition indices in :math:`R_{ij}`. 

The ``gibbs`` and ``gibbsorder`` arrays are of length :math:`N` while ``w`` and ``b`` are arrays of length :math:`N^2` representing the :math:`N \times N` matrices, row by row. 
Note that it is *NOT* assumed that the phase transition parameters are the same when going from :math:`i \rightarrow j` and  :math:`j \rightarrow i`, that is, :math:`\{\nu,B\}_{ij} \neq \{\nu,B\}_{ji}`.
Also note that :math:`R_{ii} = 0`.

``logRij`` is an array of length :math:`N (N-1)/2` giving the logarithm of the non-zero mass transportation rates between phases.  
Using ``gibbsorder`` indices, :math:`j` and :math:`k`, we see that all :math:`R_{jk}` with :math:`j \geq k` are :math:`0` because of 
:math:`\theta(G_j - G_k) = 0` (since :math:`G_k > G_j`). Writing :math:`R_{jk}` on matrix form would give that only the upper triangular part is non-zero, 
giving :math:`N (N-1)/2` non-zero elements.

The ``fromto`` array 
gives which two phases in the ``gibbs`` array, each rate in ``logRij`` is associated with: ``logRij[k]`` is the logarithm of the mass 
transportation rate from/to phases ``fromto[k]``, with ``k`` a phase transition index according to the figure below. The integer
in ``fromto``, "ij", is composed from the ``gibbs`` index of the "from" phase, :math:`i`, and the ``gibbs`` index of the "to" phase, 
:math:`j`, as :math:`i*10+j`, and with a single digit integer, "x", interpreted as "0x".

.. image:: ../GibbsOrder.pdf
  :width: 500
  :alt: Figure: How the phase transition index used in several arrays relate to the phase index in the  gibbsorder array.

``gibbsorder`` can be obtained with any sorting algorithm. In ``singularity-eos``, ``SortGibbs`` can be used

.. code-block:: cpp
 
 SortGibbs(const int num_phases, const Real *gibbs, int *gibbsorder) 

where ``num_phases`` is :math:`N`, ``gibbs`` is an array with length :math:`N`, with Gibbs free energy for each phase.
``gibbsorder`` gives the indices of ``gibbs``, the phase indices, in order from highest Gibbs free energy to lowest Gibbs free energy (see figure above). This means
that ``gibbs[gibbsorder[0]]`` is the highest Gibbs free energy and ``gibbs[gibbsorder[N-1]]`` is the lowest, that is, the Gibbs free energy of the equilibrium phase.
                                        
                                        
The time step
'''''''''''''

If a timestep would be truely infinitesimal, :math:`R_{ij} dt \leq 1` would always hold, since however big the 
rate :math:`R_{ij}` is, :math:`dt < \frac{1}{R_{ij}}`. This means that the new
mass fractions would always obey :math:`0 \leq \mu_i \leq 1`. However, with a discretized time step, :math:`R_{ij} \Delta t` can become larger than :math:`1`, and it can be that even
if the master equation holds, it results in some phase mass fractions becomming negative and some being above :math:`1`, which is unphysical. 

One way of dealing with this is to use a time step, :math:`\Delta t`, that is smaller than the inverse of the largest rate from an active phase. A routine 
suggesting a maximum timestep is available in ``singularity-eos``:  

.. code-block:: cpp
 
 Real LogMaxTimeStep(const int num_phases, const Real *mfs, 
                     const int *gibbsorder, const Real *logRij)

where ``num_phases`` is :math:`N`, ``mfs`` is the array containing each phase's (old) mass fraction, :math:`\mu_i`, ``gibbsorder`` contains the ``gibbs`` indices of the phases, 
ordered from the phase with the largest to the phase with the smallest Gibbs free energy (see figure above), and ``logRij`` contains :math:`\log(R_{ij})`. The function
gives out :math:`\log(\Delta t_{max})` since it will be used together with :math:`\log(R_{ij})` and both can be very large numbers but with opposite signs so that the difference is
small and can be safely evaluated inside an exponential.

The update method
'''''''''''''''''

Because of the numerical sensitivity to the size of the time step, several different methods have been developed
for how to perform the update. The first method made available in ``singularity-eos`` is suitable for simulations where a small timestep can be used:

.. code-block:: cpp

 SmallStepMFUpdate(const Real logdt, const int num_phases, const Real *massfractions,
                   const int *gibbsorder, const Real *logRij, Real *dmfs, Real *newmfs) 

where ``logdt`` is :math:`\log(\Delta t)`, ``num_phases`` is :math:`N`, ``massfractions`` is the array containing each phase's (old) mass fraction,
``gibbsorder`` contains the indices of the phases, ordered from the phase with the largest to the phase with 
the smallest Gibbs free energy (see figure above), ``logRij`` is :math:`\log(R_{ij})`, ``dmfs`` is the mass transformed from one phase to another,
:math:`\mu_i R_{ij} \Delta t`, for each phase transition in the order described in the figure above, and ``newmfs`` is containing the new, updated, 
mass fractions.

The advantage of using Gibbs ordered phases in the internal calculations is shown in the figure above. 
All phase transitions will always go from a higher Gibbs free energy phase to a smaller Gibbs free energy phase, 
and by using the indexing scheme in the figure all mass transformed will always go from a phase with a lower index 
to a phase with a larger index. In addition, the rates are usually larger when the Gibbs free energy difference is larger
(even though the material and phase transition fitting constants could reverse the order of the rates), and dealing with the phase transitions
in the order shown in the figure facilitates the calculations. Using the Gibbs order indices, the connection between these indices, :math:`j` and :math:`k`, and 
the phase transition indices :math:`jk` is

.. math::

 jk = (j+1)(N-1) - (j-1)j/2 - k

as can be verified by hand in the figure above.

Note that this method depleats phases in order of mass transfer to the lowest Gibbs free energy first, then the next lowest, and so on (see figure above), 
but stops once the originating phase is depleated. If this is reflecting the physical reality, is up to the user to decide. The size of the time step 
problem is taken care of by never transfering more that the existing mass in a phase to any other phase.


