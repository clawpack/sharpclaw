=======================
Problem Parameters
=======================

The following values should be set in :file:`setrun.py`.

Problem definition
===================

**ndim** (*integer*)
    Number of dimensions in problem domain (1 or 2).

**mx** (*integer*)
    Number of computational cells.

**meqn** (*integer*)
    Number of equations in the system to be solved.

**mwaves** (*integer*)
    Number of waves with non-zero speed in Riemann solution.

**maux** (*integer*)
    Number of auxiliary variables.

**mcapa** (*0 or 1*)
    * 0: No kappa array used.
    * 1: Kappa array used.

**xlower**, **xupper** (*real*)
    Endpoints of computational domain.

**t0** (*real*)
    Initial simulation time.

**src_term** (*0 or 1*)

**mbc**
    Number of ghost cells.

**mthbc_xlower**, **mthbc_xupper**
    Boundary condition to be used at left and right boundaries.

    * 0: User-specified boundary condition (must modify :file:bcN.f
      to use this option).
    * 1: Zero-order extrapolation (non-reflecting)
    * 2: Periodic (must be specified at both boundaries).
    * 3: Solid wall BC for systems where *q(2)* is the normal velocity.

Output Options
==================

**outstyle** (*1, 2, or 3*)
    * 1: Produce *nout* outputs at equally spaced intervals.
    * 2: Produce output at a provided list of times (*tout*).
    * 3: Produce output every *iout* timesteps with a total of *ntot* time steps.

**verbosity** (*0 or 1*)
    * 0: Output a message to screen only when writing output files.
    * 1: Output a message to screen at every time step.

Time stepping
================

**dt_variable** (*0 or 1*)
    * 0: Use *dt_initial* for the size of all time steps.
    * 1: Adjust the time step to achieve a desired CFL number.

**dt_initial** (*real*)
    Step size for first time step.

**dt_max**
    Maximum allowable time step size; usually set to `10^99`.

**cfl_desired**
    The time step will be adjusted at each step in order to achieve
    approximately this CFL at the next time step.

**cfl_max**
    If this CFL number is exceeded in a step, the step will be rejected
    and retaken.  This should be set somewhat larger than *cfl_desired*
    to avoid having many rejected steps.

**max_steps**
    SharpClaw will abort if more than this number of steps is taken between
    outputs.


Numerical Scheme
====================

**time_integrator** (*1, 2, 3, 4*)
    * 1: 2-stage, 2nd-order SSP Runge-Kutta method
    * 2: 4-stage, 2nd-order SSP Runge-Kutta method
    * 3: 3-stage, 3rd-order SSP Runge-Kutta method
    * 4: 10-stage, 4th-order SSP Runge-Kutta method

**lim_type** (*0, 1, 2*)
    The type of limiter to use:

    * 0: No limiting (piecewise-polynomial interpolation)
    * 1: TVD reconstruction
    * 2: WENO reconstruction

**mthlim** (*integer*)
    Which limiter to use.

