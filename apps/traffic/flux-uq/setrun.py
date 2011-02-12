""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
from pyclaw import data 


#------------------------------
def setrun(claw_pkg='sharpclaw'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "classic" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    assert claw_pkg.lower() == 'sharpclaw',  "Expected claw_pkg = 'sharpclaw'"

    rundata = data.ClawRunData(pkg=claw_pkg, ndim=1)
    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    probdata.add_param('u1-u2',      [2.0,2.0],       'Velocities')
    probdata.add_param('rho1-rho2',  [0.3333,0.1],  'Densities')
    #probdata.add_param('x0-y0',     [0.0,0.0],  'Coordinates of the dam center')
    #probdata.add_param('r0',        0.5,        'Initial radius')
    #probdata.add_param('hin',       2.0,        'Internal height')
    #probdata.add_param('hout',      1.0,        'External height')
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = 1
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = -45
    clawdata.xupper = 45
    

    # Number of grid cells:
    clawdata.mx = 800
    

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 1

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 1
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0
    
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.outstyle = 1

    if clawdata.outstyle==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.nout = 5
        clawdata.tfinal = 20.0

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.5, 1.0]   # used if outstyle == 2
        clawdata.nout = len(clawdata.tout)

    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 5
        clawdata.iout = [iout, ntot]
    


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.05
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.1
    clawdata.cfl_max = 0.9
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 5000

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Time integrator
    clawdata.time_integrator = 2
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 1
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [5]
    
    # User-supplied total fluctuation solver?
    clawdata.tfluct_solver = 0

    # Use characteristic decomposition in reconstruction step?
    clawdata.char_decomp = 0

    # Source term ### NOT USED!!!!
    clawdata.src_term = 0


    # Limiter type: 0=None, 1=TVD, 2=WENO
    clawdata.lim_type = 2

    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = 3
    
    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity
    
    clawdata.mthbc_xlower = 1
    clawdata.mthbc_xupper = 1
    
    return rundata
    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
	rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
    
