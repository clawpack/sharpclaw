.. SharpClaw documentation master file, created by
   sphinx-quickstart on Thu Oct  7 12:28:46 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview of SharpClaw
=====================================
SharpClaw is a library of Fortran routines designed to solve
hyperbolic systems of PDEs with arbitrarily high order accuracy. 
The numerical method used to solve these equations is 
described in David Ketcheson's Ph.D. thesis.

To solve a particular hyperbolic system, SharpClaw requires only
that the user provide a Riemann solver.  The Riemann solver should
be written in the same format as that used in Randy LeVeque's 
Clawpack sofware (available at http://www.clawpack.org).  Like Clawpack,
SharpClaw is based on a *wave propagation* approach: the Riemann
solver does not need to provide fluxes, but only the waves and their
speeds.

Contents:

.. toctree::
   :maxdepth: 2

   quickstart
   userguide
   future
   about
   bib

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

