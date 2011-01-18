=========================
Setting up a new problem
=========================

The easiest way to set up a new problem is typically to copy an existing
problem directory that is similar.


Adapting a Clawpack problem setup to SharpClaw
================================================

It is relatively quick and straightforward to adapt a problem that is
set up to run in Clawpack to run in SharpClaw.  The following modifications
need to be made:

    * Modify :file:`setrun.py`
        * Set `time_integrator`, `tfluct_solver`, `char_decomp`, `src_term`,
          `lim_type`, `mbc`
        * Set `mbc` to 3
    * Modify :file:`Makefile`
        * For this, it is best to copy an existing SharpClaw Makefile
    * Slight modifications to:
        * :file:`rpnN.f`
        * :file:`bcN.f`
        * :file:`outN.f`
    * Source terms: :file:`srcN.f` should be converted from full-discrete
      to semi-discrete form.
    * Optionally, provide :file:`tfluctN.f` and/or :file:`evec.f`

It is expected that in Clawpack 5.0, the setup for using SharpClaw or Clawpack
will be even less pronounced, as the interfaces to these latter functions will
be unified.

User-provided functions
========================

Riemann solver
---------------
The Riemann solver is the basic building block of any Godunov-type method,
including the methods in SharpClaw.  The routine required takes the form
of a wave-propagation solver, just like those used in Clawpack.  For many 
systems of hyperbolic equations, a wave propagation Riemann solver already
exists within Clawpack or SharpClaw.  For new systems of equations,
this function must be written by the user.

The inputs to the Riemann solver are the left and right states appearing
in the Riemann problem, as well as any auxiliary variables (such as spatially
varying coefficients) that are required for the solution of the Riemann problem.
The outputs of the Riemann solver routine are the waves, speeds, and fluctuations.

The following routines must be provided by the user, if necessary.

Source terms
-------------
If the problem includes source terms, these should be included in a routine
:file:`src.f`.  Unlike the source term routine for Clawpack, this routine
should not integrate the source term over a time step.  Rather, it should
compute the instantaneous value of the source term in the semi-discrete
system.

Auxiliary variables
--------------------

Eigenvectors
--------------------

Total Fluctuations
--------------------
