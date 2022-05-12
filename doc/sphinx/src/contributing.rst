.. _contributing-doc:

Contributing
=============


Fast Logs and Approximate Log Gridding
---------------------------------------

When spanning many orders of magnitude, Logarithmic grids are a
natural choice. Even spacing in log space corresponds to exponential
spacing in the original linear space. In other words, the grid spacing
is proportional to the value of the independent variable.

One can perform log-linear or log-log interpolation by simply
converting to log space, interpolating as one normally would, and then
converting back out. Unfortunately, logarithms and exponents are
transcendental functions, meaning they are expensive to compute and it
is thus expensive to transform in and out of log space.

To avoid this issue, we construct a space that is *approximately*
logarithmically spaced, but not quite exactly. The requirements for
this space are that the transformation into and out of this space is
fast to compute, continuous, differentiable, analytically invertible,
and close to taking a logarithm or exponentiation (depending on which
way you're going).

To achieve this, we leverage the internal representation of a floating
point number in the IEE standard. In particular, a floating point
number :math:`x` is represented as a mantissa and an exponent in base
2:

.. math::

   x = m 2^e

for mantissa :math:`m` and exponent :math:`e`. The mantiss is
guaranteed to be on the interval :math:`[1/2, 1)`. The standard
library of most low-level languages provides a performant and portable
routine to pick apart this represnetation, ``frexp``, which given a
number :math:`x`, return :math:`m` and :math:`e`.

The log in base 2 ``lg`` of :math:`x` is then given by the logarithm
of the mantissa plus the exponent:

.. math::

   \lg(x) = \lg(m) + e

Therefore, if we can find a fast, invertible approximation to
:math:`\lg(m)`, we will have achieved our goal. It turns out the
expression

.. math::

   2 (x - 1)

works pretty well, so we use that. (To convince yourself of this note
that for :math:`x=1/2` this expression returns -1 and for :math:`x=1`,
it returns 0, which are the correct values of :math:`\lg(x)` at the
bounds of the interval.) Thus our approximate, invertible expression
for :math:`\lg` is just

.. math::

   2 (m - 1) + e

for the mantissa and exponent extracted via ``frexp``. This differs
from :math:`lg` by a maximum of about 0.1, which translates to at most
a 25 percent difference. As discussed above, however, the function
itself is an exact representation of itself and the difference from
:math:`lg` is acceptable.

To invert, we use the built in function that inverts ``frexp``,
``ldexp``, which combines the mantissa and exponent into the original
floating point representation.
