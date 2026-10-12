.. _philosophy:

Design Philosophy
=================

Here we discuss some of the design principles that went in to
``singularity-eos``.

Host-Code First
-----------------

``singularity--eos`` is designed for use in continuum dynamics code
from a "host-code" first perspective. In other words, it serves up
what a fluid dynamics code may need, in the way it needs it. No more
and no less.

This guiding principle motivates several design decisions:

* We are not interested in being a one-stop-shop for material physics,
  we expose only the kind of material physics used in continuum
  dynamics.
* A host code should be able to use the looping and data access
  patterns it wants to. We provide functions that operate on *scalar*
  quantities, rather than vectors of data. If a host-code needs to
  loop over vector quantities, we use the ``indexer`` notation, rather
  than explicitly demanding contiguous arrays of data, so that the
  host code can choose its memory layout if it likes.
* We trust our users to know the memory layout, access patterns, and
  information they need. We thus are not afraid to expose, e.g.,
  initial guesses, root finder caching, etc, in the name of
  performance and flexibility.
* When possible, we expose trades-spaces such as accuracy vs. memory
  usage or computational cost.

Material Physics is Messy, but Reproducibility is Important
-------------------------------------------------------------

We recognize that material physics is complicated and try to bake this
understanding into the code. This motivates, for example, the multiple
kinds of PTE equilibrium solvers, as well as the ``preferredInput``
and ``FillEos`` machinery.

This also motivates the hot-swappable modifiers for equation of state
models. Capturing real-world material, or production-hardening an
equation of state can be challenging and imperfect. The modifiers
allow the user to ensure an EOS meets their needs, in a way that is
reproducible and comparable across host codes.

No Compromise on Performance Portability
------------------------------------------

All pieces of ``singularity-eos`` are performance portable and will
run natively on CPU, GPU, and whatever comes next.

Performance, Flexibility, Usability, and Extendability
--------------------------------------------------------

We recognize that performance, runtime usability and flexibility, and
extendability are all important, and do our best to navigate this
trade-space. We write our code in as modular a way as possible, but we
recognize that sometimes abstraction gets in the way.
