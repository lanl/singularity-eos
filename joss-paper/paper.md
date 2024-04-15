---
title: 'Singularity-EOS: Performance Portable Equations of State and Mixed Cell Closures'
tags:
  - C++
  - Fortran
  - Python
  - equations of state
  - thermodynamics
  - performance portability
authors:
  - name: Jonah M. Miller
    orcid: 0000-0001-6432-7860
    affiliation: "1, 2"
    corresponding: true
  - name: Daniel A. Holladay
    affiliation: "2, 3"
  - name: Jeffrey H. Peterson
    orcid: 0000-0001-9425-4674
    affiliation: 4
  - name: Christopher M. Mauney
    affiliation: "2, 5"
  - name: Richard Berger
    orcid: 0000-0002-3044-8266
    affiliation: 3
  - name: Anna Pietarila Graham
    affiliation: 5
  - name: Karen Tsai
    affiliation: 3
  - name: Brandon Barker
    affiliation: "1, 2, 6, 7"
  - name: Alexander Holas
    affiliation: "2, 3, 8"
  - name: Ann E. Wills
    affiliation: 9
  - name: Mariam Gogilashvili
    affiliation: "1, 2, 7, 10"
  - name: Joshua C. Dolence
    affiliation: "1, 2"
  - name: Chad D. Meyer
    affiliation: 11
  - name: Sriram Swaminarayan
    affiliation: 3
  - name: Christoph Junghans
    orcid: 0000-0003-0925-1458
    affiliation: 3
affiliations:
  - name: CCS-2, Computational Physics and Methods, Los Alamos National Laboratory, USA
    index: 1
  - name: Center for Theoretical Astrophysics, Los Alamos National Laboratory, Los Alamos, NM
    index: 2
  - name: CCS-7, Applied Computer Scienc, Los Alamos National Laboratory, USA
    index: 3
  - name: XCP-2, Eulerian Codes, Los Alamos National Laboratory, USA
    index: 4
  - name: HPC-ENV, HPC Environments, Los Alamos National Laboratory, USA
    index: 5
  - name: Department of Physics and Astronomy, Michigan State University, USA
    index: 6
  - name: Center for Nonlinear Studies, Los Alamos National Laboratory, USA
    index: 7
  - name: Heidelberg Institute for Theoretical Studies, Germany
    index: 8
  - name: XCP-5, Materials and Physical Data, Los Alamos National Laboratory, USA
    index: 9
  - name: Department of Physics, Florida State University, USA
    index: 10
  - name: XCP-4, Continuum Models and Numerical Methods, Los Alamos National Laboratory, USA
date: 13 April 2024
bibliography: refs.bib

---

# Summary

We present Singularity-EOS, a new performance-portable library for
equations of state and related capabilities. Singularity-EOS provides
a large set of analytic equations of state, such as the Gruneisen
equation of state, and tabulated equation of state data under a
unified interface. It also provides support capabilities around these
equations of state, such as Python wrappers, solvers for finding
pressure-temperature equilibrium between multiple equations of state,
and a unique *modifier* framework, allowing the user to transform a
base equation of state, for example by shifting or scaling the
specific internal energy. All capabilities are performance portable,
meaning they compile and run on both CPU and GPU for a wide variety of
architectures.

# Statement of need

When expressed mathematically for continuous materials, the laws of
conservation of mass, energy, and momentum form the Navier-Stokes
equations of fluid dynamics. In the limit of zero molecular viscosity,
they become the Euler equations. These laws can describe everything
from the flow of air over an airplane wing to the cataclysmic deaths
of stars. However, the fluid equations are not complete, and the
system must be *closed* by a description of the material at a
sub-continuum (e.g., molecular or atomic) scale. This closure is
commonly called the *equation of state* (EOS).

Equations of state vary from the simple ideal gas law, to
sophisticated descriptions multi-phase descriptions of the lattice
structure of ice or wood, to models of quark-gluon plasma and nuclear
pasta at ultra high densities. A common form to write an equation of
state is as a pair of relations:

$$P = P(\rho, T, \vec{\lambda}) \text{ and } \varepsilon = \varepsilon(\rho, T, \vec{\lambda}),$$

which relate the pressure $P$ and specific internal energy
$\varepsilon$ to density $\rho$, temperature $T$, and potentially some
unknown set of additional quantities $\vec{\lambda}$. However, other
representations are possible, and in common parlance an EOS is the
collection of knowledge needed to reconstruct some intrinsic
thermodynamic quantities from others. For example, the speed of sound
through a material or the specific heat capacity, which are
thermodynamic derivatives of the pressure and the specific internal
energy, are both determined by the EOS.

In multi-material fluid dynamics simulations, one often will end up
with a so-called *mixed cell*, where two materials exist within
the same simulation zone. This can be an artifact of the numerical
representation; for example a steel bar and the surrounding air may
end up sharing a finite volume cell if the boundaries of the cell do
not align exactly with the surface of the steel bar. Or it may
represent physical reality; for example, air is a mixture of nitrogen
and oxygen gases, as well as water vapor. Regardless of the nature of
the mixed cell, one must somehow provide to the fluid code what the
material properties of the cell are as a whole. This is called a
*mixed cell closure.* One such closure is
*pressure-temperature equilibrium* (PTE), where all materials
in the cell are assumed to be at the same pressure and temperature.

# State of the Field

Typically fluid dynamics codes each develop an EOS package
individually to meet a given problem's needs. Databases of tabulated
equations of state, such as the Sesame [@sesame] and Stellar
Collapse [@stellarcollapseweb] databases often come with tabulated
data readers, for example, the EOSPAC library
[@PimentelDavidA2021EUMV] and Stellar Collapse library
[@stellarcollapsetables]. However, these libraries typically do
not include analytic equations of state or provide a unified API. They
also don't provide extra equation-of-state capabilities, such as
equilibrium solvers or production hardening. With a few exceptions,
these libraries are also typically not GPU-capable.

We present Singularity-EOS, which aims to be a "one stop shop" for
EOS models for fluid and continuum dynamics codes. It provides a
unified interface for both analytic and tabulated equations of
state. It also provides useful surrounding capabilities, such as
Python wrappers, *modifiers*, which allow the user to transform a
given EOS, and solvers which can find the state in which multiple
EOS's are in PTE. To support usability, the library is extensively
documented and tested and supports builds through both ``cmake``
and ``Spack`` [@spack]. Singularity-EOS leverages the
``Kokkos`` [@Kokkos] library for performance portability,
meaning the code can run on both CPUs and GPUs, as well as other
accelerators. This fills an important need, as modern super computing
capabilities increasingly rely on GPUs for
performance. Singularity-EOS is now used in the ongoing open-source
[Phoebus](https://github.com/lanl/phoebus) project which
has a separate code paper in-prep.

# Design Principles and Feature Highlights

Here we enumerate several design principles underlying
Singularity-EOS, and highlight a few feature of the library.

## Flexibility in loop patterns

Singularity-EOS provides both scalar and vector APIs, allowing the
user to make EOS calls on both single points in thermodynamic space,
and on collections of points. The vector calls may be more performant
(as they may vectorize), however care is made to ensure both APIs
operate at acceptable performance, to accommodate different code
structures downstream.

## Flexibility in memory layout

The vector calls in Singularity-EOS use an *accessor* API and (with a
few exceptions) accept any C++ object that has a ``operator[]``
function defined. This allows users to lay out their memory as they
see fit and use Singularity-EOS even on strided or sparsely allocated
memory.

## Expose APIs to aid performance

Many equations of state are most naturally represented as functions of
density and temperature. However, fluid codes require pressure as a
function of density and internal energy. Extracting this often
requires computing a root find to invert the relation

$$\varepsilon = \varepsilon(\rho, T).$$

In these cases, we expose an initial guess for temperature, which
helps guarantee rapid convergence. Similarly, the performance of a
sequence of EOS calls may depend on the ordering of the calls. For
example, if both temperature and pressure are required from an
equation of state that requires inversion, requesting pressure first
will be less performant than requesting temperature first, as the
former requires two root finds, and the latter requires only one. To
enable this, we expose a function ``FillEos``, in which the user
may request multiple quantities at once, and the code uses ordering
knowledge to compute them as performantly as possible.

## Performance-portable polymorphism

Accelerators provide new challenges to standard object-oriented
programming. In particular, not all compiler stacks (such as Sycl
[@SYCL] or OpenMP Target Offload [@chandra2001parallel])
support relocatable device code, which is required for standard C++
polymorphism. Even in programming models, such as CUDA [@cuda],
which do support relocatable device code, polymorphism can be slower
than naively expected, and the user-level API can be cumbersome,
requiring operations such as ``placement new``. To sidestep these
issues, we use the C++ language feature ``std::variant`` to
implement a polymorphism mechanism that works on device.

## Modifiers

A given code may need to modify an EOS model to make it suitable for a
given application. For example, the zero-point of the energy may need
to be shifted, a porosity model may need to be added, or the unit
system may need to be changed. We implement this with a system of
*modifiers*, which can be applied on top of an EOS in a generic
way. Modifiers may also be chained.

## Fast log-lookups

To span the required orders of magnitude, tabulated equations of state
are often tabulated on log-spaced grids. Logarithms and exponentials
are, however, expensive operations and the performance of lookups can
suffer. We instead use the not-quite-transcendental lookups described
in [@NQT] to significantly enhance performance of log-like lookups.

# Acknowledgements

This work was supported through the Laboratory Directed Research and
Development program, the Center for Space and Earth Sciences, and the
center for Nonlinear Studies under project numbers 20240477CR-SES and
20220564ECR at Los Alamos National Laboratory (LANL). LANL is operated
by Triad National Security, LLC, for the National Nuclear Security
Administration of U.S. Department of Energy (Contract
No. 89233218CNA000001). This research used resources provided by the
Darwin testbed at LANL which is funded by the Computational Systems
and Software Environments subprogram of LANL's Advanced Simulation and
Computing program (NNSA/DOE). This work is approved for unlimited
release with report number LA-UR-24-23364.
