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

    probdata.add_param('ic',   9,  'IC')
    probdata.add_param('a',    1.0,  'Pulse half-width')
    probdata.add_param('x0',   -4.0,  'Pulse center')
    probdata.add_param('ninter',   1,  '# of interfaces')
    probdata.add_param('interface',   0.0,  'interface locations')
    probdata.add_param('rho',   [1.0,4.0],  'Densities')
    probdata.add_param('c',   [1.0,0.5],   'Sound speeds')

    
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
    clawdata.xlower = -10.0
    clawdata.xupper =  10.0

    # Number of grid cells:
    clawdata.mx = 200
    
    

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 2

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 2
    
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
        clawdata.tfinal = 8.0

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
    clawdata.dt_initial = 0.025
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 2.85
    clawdata.cfl_max = 2.9
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 50000

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Time integrator
    clawdata.time_integrator = 4
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 2
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [2,2]
    
    #User-supplied total fluctuation solver?
    clawdata.tfluct_solver = 1

    #Use characteristic decomposition in reconstruction step?
    clawdata.char_decomp = 1

    # Source terms?
    clawdata.src_term = 0
    
    
    # Limiter type: 0=None, 1=TVD, 2=WENO
    clawdata.lim_type = 1

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 3)
    clawdata.mbc = 3
    
    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity
    
    clawdata.mthbc_xlower = 2
    clawdata.mthbc_xupper = 2
    
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
    
