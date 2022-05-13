.. _using-eos:

The Equation of State API
=========================

This page describes the equation of state API in detail. For just the
information needed to get started, check out the :ref:`getting started
<getting-started>` page.

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
<models>. Unlike with standard polymorphism, you don't need to
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

.. code-block::

   auto my_ideal_gas = my_eos.get<singularity::IdealGas>();

This will give you access to methods and fields which may be unique to
a class but not shared by the ``Variant``.

