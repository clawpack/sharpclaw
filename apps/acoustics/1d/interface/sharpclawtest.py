#!/usr/bin/env python
#
# Python script to run xsclaw with several different parameter choices.
#

import shelve
from pylab import *
from clawtools import *
import actools as ac


#-------------------------------------
def runtests(testname='interface 1', plotsoln=False):
#-------------------------------------
    '''
    ================================
    Main code to run a set of tests:
    ================================
    '''

    if plotsoln:
        from plotclaw import plotclaw
    
    # set data values:
    data = ClawData()
    parms= ClawData()
    data.tfinal = 8     # final time
    data.ic =  9           #C^5 Newton polynomial initial condition
    data.a  =  1.0         #Pulse half-width
    data.x0 = -4.0         #Pulse half-width
    data.xlower = -10.0
    data.xupper =  10.0
    data.mthbc_xlower = 2
    data.mthbc_xupper = 2
    data.mthbc=[2,2]
    parms.time_integrator = 4
    parms.cfl_max         = 2.9
    parms.cfl_desired     = 2.85

    if plotsoln:
        data.nout = 10          # output at several times for plotting
        data.write('setplot1.data')
    else:
        data.nout = 1          # only output at final time

    if testname=='homogeneous':
        # Homogenous medium
        data.ninter=1
        data.interface=0.0
        data.c = 1.0,1.0
        data.rho = 1.0,1.0
        data.write('setprob.data')
        data.write('sharpclaw.data')
        #Set parameters that are uniform for all runs
        parms.char_decomp     = 0
        parms.lim_type        = 0
        parms.write('sharpclaw.data')
        parms.write('claw.data')

        #Set run-specific parameters
        runnames=['TVD2','WENO5','UC7','ClawPack']
        runparms={}
        for name in runnames: 
            runparms[name]=ClawData()
            runparms=setparms(runparms,name)
        mxvals = array([200,400,800,1600,3200])        # mx values 

    elif testname=='interface 1':
        # Weak interface
        data.ninter=1
        data.interface=0.0
        data.c = 1.0,0.5
        data.rho = 1.0,4.0
        data.write('setprob.data')
        data.write('sharpclaw.data')
        #Set parameters that are uniform for all runs
        parms.char_decomp     = 0
        parms.lim_type        = 0
        parms.write('sharpclaw.data')
        parms.write('claw.data')

        #Set run-specific parameters
        #runnames=['WENO5 Wave','UC3 Wave','UC7 Wave']
        runnames=['TVD2','WENO5','UC7','ClawPack']
        runparms={}
        for name in runnames: 
            runparms[name]=ClawData()
            runparms=setparms(runparms,name)
        mxvals = array([200,400,800,1600,3200])        # mx values 

    elif testname=='interface 2':
        # Strong interface
        data.ninter=1
        data.interface=0.0
        data.c = 1.0,0.5
        data.rho = 1.0,4000.0
        data.write('setprob.data')
        data.write('sharpclaw.data')
        data.write('claw.data')
        #Set parameters that are uniform or default
        parms.char_decomp     = 0
        parms.lim_type        = 0

        #Set run-specific parameters
        runnames=['TVD2','WENO5','WENO5 Char','ClawPack']
        runparms={}
        for name in runnames: 
            runparms[name]=ClawData()
            runparms=setparms(runparms,name)
        mxvals = array([200,400,800,1600,3200])        # mx values 
 
    elif testname=='periodic':
        data.a = 1.0         #Pulse half-width
        data.x0 = -1.0         #Pulse center
        data.xlower = -10.0
        data.xupper =  10.0
        data.tfinal = 8.      # final time
        #data.interface=range(data.xlower+1,data.xupper)
        delta=2.0
        data.interface=list(arange(data.xlower+delta,data.xupper,delta))
        data.ninter=len(data.interface)
        print 'ninter: ',data.ninter
        data.c   = [1.,1.]*((data.ninter+1)/2)
        data.rho = [1.,3.]*((data.ninter+1)/2)
        print len(data.c)
        data.write('setprob.data')
        data.write('sharpclaw.data')
        data.write('claw.data')

        #Set run-specific parameters
        runnames=['TVD2','WENO5','WENO5 Char','ClawPack']
        runparms={}
        for name in runnames: 
            runparms[name]=ClawData()
            runparms=setparms(runparms,name)
        mxvals = array([200,400,800,1600,3200])        # mx values 
 
    else:
        print "Unrecognized input:  itest = ", itest
        sys.exit(1)

    print "------------------------------------------------------------------"

    exact = {}
    
    table = {}  # dictionary of data and results for each test
    
    for run in runnames:
        table[run] = {} # dictionary to hold data
                                          # and results for this test
        this_table = table[run]  # short name
        this_table['mxvals'] = mxvals       # grid resolutions to test
        this_table['ave_cell_area'] = (data.xupper-data.xlower) / mxvals
        this_table['tfinal'] = data.tfinal
        this_table['errors'] = empty(len(mxvals))
        this_table['runtimes'] = empty(len(mxvals))
        
        for igrid in range(len(mxvals)):
            mx = mxvals[igrid]
            parms.mx = mx
            #data.write('setprob.data')
            #os.system('dos2unix claw2ez.data')
            #os.system('dos2unix setprob.data')

            print "-------------------------------------------------------"
            print 'Running code for run %s with mx = %g ...' \
                                  % (run,mx)
        
            odir = '.'
            # or use next line to put each set of output in separate subdir:
            # odir = 'output_lim%ggrid%gmx%g' % (mthlim,igrid,mx)

            # run fortran code:
            t1 = time.time()
            if run=='ClawPack':
                parms.write('claw.data')
                runparms[run].write('claw.data')
                runclaw(outdir=odir, overwrite=True, xclawcmd='xclaw')
            else:
                parms.write('sharpclaw.data')
                runparms[run].write('sharpclaw.data')
                runclaw(outdir=odir, overwrite=True, xclawcmd='xsclaw')
            CPUtime = time.time() - t1
            this_table['runtimes'][igrid] = CPUtime
            print "CPU time = ",CPUtime
        
            # compute exact solution if not already computed
            if not exact.has_key(mx):
                exact[mx]=compute_exact_solution(odir,frame=data.nout)
            # compute errors:
            try:
                errors = compute_errors(exact[mx],odir,frame=data.nout)
            except:
                errors = array([inf])
    
            # max-norm of error:
            errormax = abs(errors).max()   

            # approx 1-norm of error:
            errorsum = abs(errors).sum() 
            error1 = errorsum * this_table['ave_cell_area'][igrid]
    
            this_table['errors'][igrid] = error1
        
            if plotsoln:
                # put each set of plots in a different directory:
                pdir = 'plots_%g_%gmx%g' % (prm1,prm2,data.mx)
                plotclaw(outdir=odir, plotdir=pdir, overwrite=True, \
                     indexentry='%s = %2i, %s = %2i, mx = %4i, \
                     ' % (param1,prm1,param2,prm2,mx))
    
    # save results:
    db = shelve.open('table.db')
    db['tfinal'] = data.tfinal
    db['table'] = table
    db['runnames'] = runnames
    db.close
    
    # output tables:
    #make_text_table(table, runnames, mxvals, fname='errortables.txt')
    make_latex_table(table, runnames, testname, fname='errortables.tex')


#------------------------------------
def compute_errors(qq,odir='.',frame=0):
#------------------------------------
    # change to output directory and read in solution from fort.q000N
    # for frame N.
    os.chdir(odir)
    clawframe = read_clawframe(frame)
    t = clawframe.t
    grid = clawframe.grids[0]
    xc = grid.xcenter
    dx=grid.xedge[1]-grid.xedge[0]

    errors = qq[0,:] - grid.q[:,0]
    ion()
    clf()
    plot(xc,qq[0,:],'k',xc,grid.q[:,0],'+b')
    draw()
    ioff()
    return errors

#------------------------------------
def compute_exact_solution(odir='.',frame=0):
#------------------------------------
    from scipy.integrate import fixed_quad
    ap=ac.AcousticsProblem('setprob.data','sharpclaw.data')

    # change to output directory and read in solution from fort.q000N
    # for frame N.
    os.chdir(odir)
    # read in clawpack solution for this frame:
    clawframe = read_clawframe(frame)
    t = clawframe.t
    grid = clawframe.grids[0]
    xc = grid.xcenter
    dx=grid.xedge[1]-grid.xedge[0]
    
    print "Computing exact solution for mx = ",grid.mx
    qq=zeros([2,size(xc)])
    for ii in range(size(xc)):
	qq[0,ii],dummy = fixed_quad(ap.pvec,grid.xedge[ii],grid.xedge[ii+1],args=(t,))
	qq[1,ii],dummy = fixed_quad(ap.uvec,grid.xedge[ii],grid.xedge[ii+1],args=(t,))
    qq/=dx
    return qq
    
    
#-------------------------------------------------------------------
def make_text_table(table, runnames, mxvals, fname='errortables.txt'):
#-------------------------------------------------------------------
    errfile = open(fname,'w')
    errfile.write('Errors computed by clawtest.py on %s \n\n' \
                   % current_time())
    for iparam1 in range(len(param1vals)):
	prm1=param1vals[iparam1]
        errfile.write('\n\n%s = %g\n===========\n'  % (param1,prm1))
        for iparam2 in range(len(param2vals)):
	    prm2=param2vals[iparam2]
            errfile.write('\n\n%s = %g\n-------------\n'  % (param2,prm2))
            errfile.write('  mx    ave Delta x         error  ')
            errfile.write('      observed order       CPU Time\n\n')

            # table values for this table:
            te = table[(prm1,prm2)]   
            for igrid in range(len(te['mxvals'])):
                mx = te['mxvals'][igrid]
                ave_cell_area = te['ave_cell_area'][igrid]
                ave_dx = ave_cell_area
                errornorm = te['errors'][igrid]
		CPUtime = te['runtimes'][igrid]
                if igrid>0:
                    E1 = te['errors'][igrid-1]
                    E2 = te['errors'][igrid]
                    A1 = te['ave_cell_area'][igrid-1]
                    A2 = te['ave_cell_area'][igrid]
                    porder = log(E1/E2)/log(A1/A2)
                else:
                    porder = nan
                errfile.write('%4i  %14.4e  %16.4e  %16.4e %16.2f\n' \
                            % (mx, ave_dx, errornorm, porder, CPUtime))

    errfile.write('\n\nTotal CPU time: %g seconds\n\n' % time.clock())
    errfile.close()

    print "------------------------------------------------------------------"
    print "Table of errors is in file ",   fname
    print "------------------------------------------------------------------"
    

#-------------------------------------------------------------------
def make_latex_table(table, runnames, testname, fname='errortables.tex'):
#-------------------------------------------------------------------

    errfile = open(fname,'w')
    errfile.write(r'% ')
    errfile.write(' Errors computed by clawtest.py on %s \n\n' \
                   % current_time())
    
    headerline=""
    headerline2="mx "
    for name in runnames:
        headerline+="& \multicolumn{2}{c|}{"+name+"}"
        headerline2+="& Error & Order"
    headerline+="\\\\ \\hline \n"
    headerline2+="\\\\ \\hline \n"
    errfile.write(r'''

      \begin{table}[th]
      \caption{Errors for %s problem \label{tbl:%s}}
      \begin{center}
    '''  % (testname, testname))
    errfile.write("\\begin{tabular}{c|"+"cc|"*len(runnames)+"}"+"\n")
    errfile.write(headerline)
    errfile.write(headerline2)

    mxvals=table[runnames[0]]['mxvals']

    # lines in table for each grid resolution:
    for imx in range(len(mxvals)):
        t0=table[runnames[0]]
        mx = t0['mxvals'][imx]
        ave_cell_area = t0['ave_cell_area'][imx]
        ave_dx = ave_cell_area
        errfile.write('\\shadeRow'*mod(imx,2)) #Shade even rows
        errfile.write('%4i ' % (mx))
        for run in runnames:
            t0=table[run]
            mx = t0['mxvals'][imx]
            ave_cell_area = t0['ave_cell_area'][imx]
            ave_dx = ave_cell_area
            errornorm0 = t0['errors'][imx]
            if imx>0:
                E1 = t0['errors'][imx-1]
                E2 = t0['errors'][imx]
                A1 = t0['ave_cell_area'][imx-1]
                A2 = t0['ave_cell_area'][imx]
                porder0 = log(E1/E2)/log(A1/A2)
            else:
                porder0 = nan
            errfile.write('& %16.2e & %16.2f' % (errornorm0, porder0))
            #errfile.write(' & %16.4f & %16.2f \\\\ \n' \
            #            % (errornorm3, porder3))
        errfile.write('\\\\ \n')


    errfile.write(r'''
      \hline
      \end{tabular}
      \end{center}
      \end{table}

    ''')

    errfile.close()

    print "------------------------------------------------------------------"
    print "Latex table of errors is in file ",  fname
    print "------------------------------------------------------------------"
    
#-------------------------------------------------------------------
def setparms(runparms,name):
#-------------------------------------------------------------------
    if name=='TVD2':
        runparms[name].lim_type=1
        runparms[name].char_decomp=0
        runparms[name].mthlim=[4,4]
    elif name=='UC3':
        runparms[name].lim_type=0
        runparms[name].char_decomp=0
        runparms[name].mthlim=[3,3]
    elif name=='UC3 Wave':
        runparms[name].lim_type=0
        runparms[name].char_decomp=1
        runparms[name].mthlim=[3,3]
    elif name=='UC7':
        runparms[name].lim_type=0
        runparms[name].char_decomp=0
        runparms[name].mthlim=[7,7]
    elif name=='UC7 Wave':
        runparms[name].lim_type=0
        runparms[name].char_decomp=1
        runparms[name].mthlim=[7,7]
    elif name=='WENO5':
        runparms[name].lim_type=2
        runparms[name].char_decomp=0
        runparms[name].mthlim=[5,5]
    elif name=='WENO5 Wave':
        runparms[name].lim_type=2
        runparms[name].char_decomp=1
        runparms[name].mthlim=[5,5]
    elif name=='ClawPack':
        runparms[name].mthlim=[3,3]
        runparms[name].cfl_max=1.0
        runparms[name].cfl_desired=0.9
    elif name=='TVD2 Char':
        runparms[name].lim_type=1
        runparms[name].char_decomp=2
        runparms[name].mthlim=[4,4]
    elif name=='WENO5 Char':
        runparms[name].lim_type=2
        runparms[name].char_decomp=2
        runparms[name].mthlim=[5,5]
    return runparms

if __name__ == '__main__':
    import string
    # execute runtests with any command line arguments:
    cmd = 'runtests(%s)' % string.join(sys.argv[1:])
    exec(cmd)

