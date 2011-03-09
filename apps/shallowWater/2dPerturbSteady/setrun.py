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

    rundata = data.ClawRunData(pkg=claw_pkg, ndim=2)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    probdata.add_param('grav',      9.81,        'Gravitational constant')
    probdata.add_param('eps',       0.00,        'Perturbation')
 
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = 2
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = 0.0
    clawdata.xupper = 2.0
    
    clawdata.ylower = 0.0
    clawdata.yupper = 1.0
        

    # Number of grid cells:
    clawdata.mx = 200
    clawdata.my = 100

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 3

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
        clawdata.nout = 10
        clawdata.tfinal = 0.12

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.12]   # used if outstyle == 2
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
    clawdata.verbosity = 1
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.0001
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.5
    clawdata.cfl_max = 0.6
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 5000

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Time integrator
    clawdata.time_integrator = 4

    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [5, 5, 5]

    #User-supplied total fluctuation solver?
    clawdata.tfluct_solver = 1

    #Use characteristic decomposition in reconstruction step?
    clawdata.char_decomp = 0

    # Source terms? NOT USED!!!!
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
    
    clawdata.mthbc_ylower = 1
    clawdata.mthbc_yupper = 1
    
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
    
