.. fluid2d documentation master file, created by
   sphinx-quickstart on Thu Mar  2 22:46:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2


Welcome to Fluid2d documentation
================================

Fluid2d is a versatile Python-Fortran CFD code that solves a large
class of 2D flows. It is designed to be used easily by Students
learning Fluid Mechanics or Geophysical Fluid Dynamics and by their
professors willing to illustrate their course and to organize
numerical practicals. The idea is to visualize flows on the fly, as
they are computed. The effect of parameter changes can be seen
immediately. The key quantity of fluid2D is the vorticity. If you feel
weak on vorticity dynamics, this code is for you. You should rapidly
become as expert as the experts.

You can learn how basic processes work because of the power of
animations. It is quite easy to go beyond textbooks and to reach
research questions.

Several features are particulary cool

- the code handles many different sets of equations: transport, Euler,
  quasi-geostrophic, Boussinesq and even the thermal wind equations.

- the code handles a mask system that allows to have complicated
  geometries (closed domain with arbitratry shape, reentrant channel, with
  multiple island etc)

- the code tends to have a very low level of dissipation because the
  dissipation is handled implicitely by the numerics.

- a no-slip condition that allows to study boundary layers questions

- the code is parallelized and can be run on cluster if high
  resolution is desired.

- a hopefully not too hard user's interface. You tell me. An experiment
  boils down to one script where you set everything up. For
  forced-dissipated flows or fancy add-ons (like introducing a dye
  tracer), you may need a second script.

Credits:
--------

This code have been developed over many years of teaching numerical
methods for CFD. The first version of this code was in Matlab under
the impulse of the FDSE Summer School. A special credit to Caroline
Muller who makes this code happen in the first place. A full writing
in python has been done with the great help of Clement Vic and Nicolas
Grima (CNRS engineer). Others Students have contributed to either the
coding or a better understanding of how the numerics handles
dissipation. They should be thanked here: Ljuba Novi for the work on
RK3 and up5, Alexander Siteur for his work on the pressure (not yet
implemented in this version), Milan Kloewer for his work on a level
set implementation (not present in this version), Mathieu Morvan for
the dipole wall interaction and Markus Reinert for his many valuable
feedbacks. Finally, the code has benefited from my enthusiastic
collaboration with Professor J. C. McWilliams.


Where to start?
---------------

You may have a quick look at the :ref:`gallery` to see the possibilities. Then you need to  :ref:`install` it and start your first experiment.



Cite the code
-------------

The code has never been published and will likely never be. If you use the
code for your classes or your research I would appreciate your
feedbacks (mailto: roullet AT univ-brest.fr)

Disclaimer
----------

As most codes, this code almost certainly has bugs.  Please report
them and possibly their fix if you find it.

A word of caution
-----------------

The flows simulated by fluid2d are two-dimensional. This has many
impacts. The most important to be aware of  is the systematic tendency for the
inverse cascade of kinetic energy, i.e. the tendency for the vorticity
to form bigger and bigger structures. This holds for all flows. It is
a realistic feature in many cases but not always. Because our world is
three-dimensional, vorticity is actually a vector, not a scalar. In
3D, vorticity form tubes than can get twisted like spaghettis. Flows
look then very differently than in 2D.


Contents:
---------

.. toctree::
   :maxdepth: 2

   gallery/gallery
   docs/model_equations
   docs/model_numerics
   docs/howto
