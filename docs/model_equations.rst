===============
Model equations
===============

.. toctree::
   :maxdepth: 1


The choice of the model's equations is the first one to be made for an
experiment. It is controlled by

.. code::

   param.modelname = 'euler'

The common core of fluid2d's equations is a vorticity formulation
allowed by the incompressibility assumption. The velocity is deduced
from the streamfunction, itself deduced from the vorticity. The model
thus integrates in time a scalar quantity, the vorticity, instead of a
vector quantity, the velocity. Pressure is totally absent from both
the equations and the code. Of course, pressure exists. For all the
flows you can simulate there is always an underlying pressure
field. But it has to be diagnosed and so far the diagnostic is not
implemented (it was in an older version). A nice property of a
vorticity formulation is that the Kelvin's theorem is enforced at every
grid point and at all time. This is a very desirable property.


Euler equations
===============

This is the model for the dynamics of vorticity

.. math::

   \partial_t \omega + J(\psi,\omega) = 0

   \partial^2_{xx}\psi+ \partial^2_{yy} \psi = \omega

where :math:`\omega` is the vorticity and :math:`\psi` is the
stream-function. The boundary condition for the Poisson equation is
:math:`\psi=0` along solid boundaries. The constant 0 might be replaced
along each island by a given value :math:`\psi_k`.

There is a possibility to add an explicit viscous term on the r.h.s.,
controlled by a diffusion parameter.

This models allows to study the **dynamics of vorticity** in general from
**interaction of isolated coherent structures** to **2D turbulence**.

Transport equations
===================

In this model the stream-function is frozen. Only passive a tracer :math:`\phi` is transported

.. math::

   \partial_t \phi + J(\psi,\phi) = 0

This model is used to study the **properties of transport per se**. It can
also be used to test the **numerical properties of advection
schemes**. There is a practical on this question.

Quasi-geostrophic equations
===========================

This model is similar to the Euler one except that the Coriolis force is now included, yet hidden in the deformation radius,

.. math::

   \partial_t \omega + J(\psi,\omega) = 0

   \partial^2_{xx}\psi+ \partial^2_{yy} \psi -R^{-2}_d \psi - \beta y = \omega

where :math:`\omega` is now the total potential vorticity (PV), :math:`R_d` is
the Rossby deformation radius and :math:`\beta` is the beta parameter
encoding the Earth's curvature. It is worth noting that quite often the quasi-geostrophic equations are written for the PV anomaly :math:`\omega'` such that

.. math::
   \omega = \omega' + \beta y

The evolution equation for the PV anomaly involves the :math:`-\beta
v` source term. To preserve a conservative form for PV, fluid2d uses
the total PV.


This model allows to study **wind-driven circulation, Rossby waves,
turbulence on the beta plane** etc.

Boussinesq equations
====================

This model describes the motion of a stratified fluid in a 2D vertical slice. The equations are the same as Euler, complemented by a transport equation for buoyancy :math:`b`. The vorticity equation is coupled to the buoyancy via the torque of the weight, the :math:`\partial_x b` term

.. math::

   \partial_t \omega + J(\psi,\omega) = \partial_x b

   \partial_t b + J(\psi,b) = 0

   \partial^2_{xx}\psi+ \partial^2_{yy} \psi = \omega

This model allows to study stratified processes such as: **internal
waves, convection, Kelvin-Helmholz instabilities, lee-waves** etc.

Thermal wind equations
======================

This model describes the secondary (ageostrophic) circulation that develops around an axisymmetric jet. The equations couple the Boussinesq-like dynamics in the vertical plane  to the jet speed

.. math::

   \partial_t \omega + J(\psi,\omega) = \partial_x b - f\partial_z v

   \partial_t v + J(\psi,v) = f \partial_z\psi

   \partial_t b + J(\psi,b) = 0

   \partial^2_{xx}\psi+ \partial^2_{yy} \psi = \omega

where :math:`v` is the (horizontal) jet velocity and :math:`f` is the
Coriolis parameter. In the absence of secondary circulation,
i.e. :math:`\psi=0`, the buoyancy and the jet velocity are
in thermal wind-balance (first equation).

This model allows to study **symmetric instabilities, frontogenesis**
etc. This model is not yet fully debuged.


.. _viscous-interaction:
Viscous Interaction
===================

The no-slip boundary condition is handled by added a localized source
term :math:`S` along the boundary. The source term is recomputed at each time
step. It is adjusts such that the tangential velocity along the
boundary vanishes.

Last update: 12/10/2021
