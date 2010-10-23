import os,sys,shutil,glob
import string,re
import time

try:
    import subprocess
except:
    print 'Must use a more recent version of Python with module subprocess'
    print 'Python 2.5 is recommended'
    sys.exit(1)  

try:
    from numpy import *
except:
    'NumPy not found'
    sys.exit(1)  

#===========================================================================
#
#

#================
class ClawData:
#================
    """
    Classes and functions for manipulating clawpack input data in files *.data
    """

    #=============================
    def __init__(self,*datafilelist):
    #=============================
        self.datafiles = []
        for datafile in datafilelist:
            self.read(datafile)

    #=============================
    def read(self,*datafilelist):
    #=============================
        """
        Read data in from a clawpack style data file.
        Any lines of the form
           values  =:  name
        is used to set an attribute of self.
        """
        verbose = False

        if len(datafilelist) == 0:
            datafilelist = self.datafiles

        for datafile in datafilelist:
            if datafile not in self.datafiles:
                self.datafiles.append(datafile)
            if verbose:
                print "    Reading data file:  ",datafile,"  <p>"
            file = open(datafile,'r')
            regexp = re.compile(r"(?P<values>.*)=:(?P<name>.*)")
            for line in file:
                result = regexp.search(line)
                if result:
                    name = re.split(r'\s+', result.group('name').strip())[0]
                    values = re.split(r'\s+', result.group('values').strip())
		    if values[0]=='T':
		        values = ['True']
		    if values[0]=='F':
		        values = ['False']
                    if len(values) > 0:
                        setval = "self." + name + " = " +  ','.join(values)
                    else:
                        print "Error -- empty values for name ", name
                        
                    exec(setval)
                
    #=============================
    def write(self,*datafilelist):
    #=============================
        """
        Write data out to a clawpack datafile
        Any lines of the form
           values  =:  name
        are overwritten using the current attribute of self for values.

        If datafilelist is empty, all data files in self.datafiles are
        overwritten.
        """
        verbose = False

        if len(datafilelist) == 0:
            datafilelist = self.datafiles

        if len(datafilelist) == 0:
            print "Error - no data files specified"
            return

        for datafile in datafilelist:
            try:
                file = open(datafile,'r')
            except:
                print "Error -- data file does not exist: ",datafile
                return

            newdatafile = datafile + 'xxxxx'
            try:
                newfile = open(newdatafile,'wb')
            except:
                print "Error -- data file does not exist: ",datafile

            regexp = re.compile(r"(?P<values>.*)=:(?P<name>.*)")
            for line in file:
                result = regexp.search(line)
                if result:
                    name = re.split(r'\s+', result.group('name').strip())[0]
                    values = re.split(r'\s+', result.group('values').strip())
                    if hasattr(self,name):
                        newvalues = getattr(self,name)
                        if isinstance(newvalues,tuple) \
                                   | isinstance(newvalues,list):
                            # remove brackets or parens:
                            newstring = repr(newvalues)[1:-1]
                            # remove commas:
                            newstring = newstring.replace(',','')
                        # Need to fix?
                        #elif isinstance(newvalues,bool):
                        elif 0:
			    if newvalues==True:
                                newstring = 'T'
                            elif newvalues==False:
                                newstring = 'F'
                        else:
                            newstring = repr(newvalues)
                         
                        newstart = string.ljust(newstring,25)
                        line = line.replace(result.group('values')+"=:", \
                                newstart +" =:")
                newfile.write(line)

            file.close()
            newfile.close()
            try:
                os.remove(datafile)
                os.rename(newdatafile,datafile)
                if verbose:
                    print "Wrote new version of data file:  ",datafile,"  <p>"
            except:
                print 'could not rename ',newdatafile, ' as ',datafile
                         

#=========================================
def change_data(datafile,param,newstring):
#=========================================
    verbose = False

    file = open(datafile,'r')
    newdatafile = datafile + 'xxxxx'
    newfile = open(newdatafile,'wb')
    regexp = re.compile(r"(?P<values>.*)=:\s*" + param)
    
    for line in file:
        result = regexp.search(line)
        if result:
            newstart = string.ljust(newstring,25)
            line = line.replace(result.group('values')+"=:", \
                                newstart +" =:")
            if verbose:
                print "Changed data value ",param," to ", newstring
        newfile.write(line)

    file.close()
    newfile.close()
    try:
        os.remove(datafile)
        os.rename(newdatafile,datafile)
    except:
        print 'could not rename ',newdatafile, ' as ',datafile



# ----------------------------------------------------------------------
#
# Classes and functions for reading clawpack solutions from files fort.*
#
# Simplest use:
#     clawframe = read_clawframe(N,outdir)
# reads data from both fort.t000N and fort.q000N 
#
# To read from fort.t but not the bigger fort.q, e.g. if you only need 
# to determine ndim:
#     clawframe = ClawFrame(N,outdir)
#     clawframe.read_fortt()
#     ndim = clawframe.ndim
# ----------------------------------------------------------------------

    
#=========================================================
def read_clawline(inputfile, type='float', filetype='ascii'):
#=========================================================
    """
    Reads one line from an input file fort.t* or fort.q*
    type is an optional type declaration, understands 'int', 'float'

    """
    line = inputfile.readline()
    l = line.split()
    if type=='int':
        return int(l[0])
    else:
        return float(l[0])


#==============
class ClawGrid:
#==============
    """
    
    ClawGrid - Data structure for a single grid.
               One Frame may consist of many grids if AMR is used.
               clawframe.grids[m] is of type ClawGrid for m=1..clawframe.ngrids

    Members:
      Grid Information - 
        gridno - Current grid number
        level - Level of the current grid in AMR structure

      Dimension Information - 
        (mx,my,mz) - Number of grid cells in each dimension
        (xlow,ylow,zlow) - Lower bound on the current grid
        (dx,dy,dz) - Spatial descretization

      Data - 
        q - Data for this frame of dimension (mx,my,mz,meqn)

      Grid Data - Automatically generated grid
        (xedge,yedge,zedge) - grid of size (mx+1,my+1,mz+1) with location 
          of left edge of each cell (including 1 ghost cell on right)
        (xcenter,ycenter,zcenter) - grid of size (mx,my,mz) with location
          of center of each cell

    """
    #==================
    def __init__(self):
    #==================
        # Grid information
        self.gridno = []
        self.level = []

        # Dimension information
        self.mx = []
        self.xlow = []
        self.dx = []
        self.my = []
        self.ylow = []
        self.dy = []
        self.mz = []
        self.zlow = []
        self.dz = []

        # Data for this frame
        self.q = []

        # Grid data
        self.xedge = []
        self.xcenter = []
        self.yedge = []
        self.ycenter = []
        self.zedge = []
        self.zcenter = []
        
        # Viewable members of this class
        self.members = ["mx","xlow","dx","my","ylow","dy","mz","zlow","dz","gridno","level"]
        self.members.sort()

    #==================
    def __repr__(self):
    #==================
        output = ""
        for attr in self.members:
            output = output+"%s = %s, " % (attr,eval("self.%s" % attr))
        return str(output)


#================
class ClawFrame:
#================
    """
    
    ClawFrame - General data structure class for all grids at one output time

    Members:
      General Frame Information - 
        frame            - Frame number N from fort.q000N for this object
        meqn             - Number of equations
        ngrids           - Number of grids
        maux             - Dimension of the auxillary array
        t                - Output time for this frame
        grids            - Array of objects of type ClawGrid of size ngrids

    Functions:
      __init__           - initialize attributes 
      read_fortt         - read t,meqn,ngrids,ndim,maux  from fort.t file
      read_fortq         - read solution from ascii fort.q file

    """

    #=====================================
    def __init__(self,frame,outdir='./'):
    #=====================================
        # General frame information
        self.frame = frame
        self.outdir = outdir
        self.meqn = -1
        self.ngrids = -1
        self.maux = -1
        self.t = -1
        self.grids = []
        
        # Viewable members of this class
        self.members = ["frame","meqn","ngrids","maux","t","grids"]
        self.members.sort()

    #==================
    def __repr__(self):
    #==================
        output = ""
        for attr in self.members:
            output = output+"%s = %s\n" % (attr,eval("self.%s" % attr))
        return str(output)
    
    #=======================================
    def read_fortt(self):
    #=======================================
        """
    
        Reads in clawpack file fort.t000N where N = self.frame
    
        """
    
        frame = self.frame
        if frame < 0:
            # Don't construct file names with negative frame values.
            print " "
            print "Frame " + str(frame) + " does not exist ***"
            clawframe = []
            return clawframe
        
    
        # Open fort.t file
	outdir = self.outdir
        fname = os.path.join(outdir, 'fort.t') + str(frame).zfill(4)
        try:
	    f = open(fname,'r')
	except:
            print " "
            print "File " + fname + " does not exist ***"
	    sys.exit(1)
    
        # Read in values from fort.t file:
        try:
            self.t = read_clawline(f)
            self.meqn = read_clawline(f,'int')
            self.ngrids = read_clawline(f,'int')
            self.ndim = read_clawline(f,'int')
            self.maux = read_clawline(f,'int')
            f.close()
        except:
            print "File " + fname + " should contain t, meqn, ngrids, ndim, maux"
    	    sys.exit(1)

    #=======================================
    def read_fortq(self, filetype='ascii'):
    #=======================================
        """
    
        Reads in clawpack file fort.q000N   where N = self.frame
	(currently only works for ascii output files)
    
        """

	if filetype != 'ascii':
	    print "Only works for filetype = 'ascii'"
	    sys.exit(1)
    
        frame = self.frame
        if frame < 0:
            # Don't construct file names with negative frame values.
            print " "
            print "Frame " + str(frame) + " does not exist ***"
	    sys.exit(1)
        
    
        # Open fort.q file

	outdir = self.outdir
        fname = os.path.join(outdir, 'fort.q') + str(frame).zfill(4)
        try:
	    f = open(fname,'r')
	except:
            print " "
            print "File " + fname + " does not exist ***"
	    sys.exit(1)
    
    
        # Read in solution from fort.q file:
        for ng in range(self.ngrids):
            grid = ClawGrid()
            grid.gridno = read_clawline(f,'int')
            grid.level = read_clawline(f,'int')
    
            grid.mx = read_clawline(f,'int')
            if self.ndim > 1:
                grid.my = read_clawline(f,'int')
            if self.ndim > 2:
                grid.mz = read_clawline(f,'int')
    
            grid.xlow = read_clawline(f)
            if self.ndim > 1:
                grid.ylow = read_clawline(f)
            if self.ndim > 2:
                grid.zlow = read_clawline(f)
    
            grid.dx = read_clawline(f)
            if self.ndim > 1:
                grid.dy = read_clawline(f)
            if self.ndim > 2:
                grid.dz = read_clawline(f)
    
            blank = f.readline()
    
            # Read the data for this grid and create the grid
            if self.ndim == 1:
                grid.xedge = empty(grid.mx+1,'d')
                grid.xcenter = empty(grid.mx,'d')
                grid.q = empty((grid.mx,self.meqn),'d')
                for i in range(0,grid.mx):
                    grid.xedge[i] = grid.xlow + i*grid.dx
                    grid.xcenter[i] = grid.xlow + (i+0.5)*grid.dx
                    line = f.readline()
                    l = line.split()
                    for m in range(0,self.meqn):
                        grid.q[i,m] = float(l[m])
    	    # right edge of last cell:
                grid.xedge[grid.mx] = grid.xlow + grid.mx*grid.dx
    
            elif self.ndim == 2:
                grid.xedge = empty(grid.mx+1,'d')
                grid.xcenter = empty(grid.mx,'d')
                grid.yedge = empty(grid.my+1,'d')
                grid.ycenter = empty(grid.my,'d')
                grid.q = empty((grid.mx,grid.my,self.meqn),'d')
                for j in range(grid.my):
                    grid.yedge[j] = grid.ylow + j*grid.dy
                    grid.ycenter[j] = grid.ylow + (j+0.5)*grid.dy
                    for i in range(grid.mx):
                        grid.xedge[i] = grid.xlow + i*grid.dx
                        grid.xcenter[i] = grid.xlow + (i+0.5)*grid.dx
                        line = f.readline()
                        l = line.split()
                        for m in range(0,self.meqn):
                            grid.q[i,j,m] = float(l[m])
                    blank = f.readline()
                grid.xedge[grid.mx] = grid.xlow + grid.mx*grid.dx
                grid.yedge[grid.my] = grid.ylow + grid.my*grid.dy
    
            else:
                print "ERROR:  3D not implemented!"
                grid.q = emtpy((grid.mx,grid.my,grid.mz,self.meqn),'d')
    
            # Add grid to the clawframe data object
            self.grids.append(grid)
            
    

#==================================================
def read_clawframe(frame, outdir='./', filetype='ascii'):
#==================================================
    """

    read_clawframe:  General function call to read in clawpack data.  
    Returns an object that is of type ClawFrame.
    
    Arguments:
    frame        - Frame number N: read files fort.t000N and fort.q000N
    outdir       - Directory containing fort files, defaults to  './'
    filetype     - Type of file to be read in 
                     (ascii is the only file type currently supported)

    Returns:
    clawframe    - Object of class ClawFrame containing all the information 
                   about the frame requested

    """

    clawframe = ClawFrame(frame,outdir)
    clawframe.read_fortt()
    clawframe.read_fortq(filetype)
    return clawframe


#------------------------------------------------------------

#-----------------------------
def current_time(addtz=False):
#-----------------------------
    # determine current time and reformat:
    time1 = time.asctime()
    year = time1[-5:]
    day = time1[:-14]
    hour = time1[-13:-5]
    current_time = day + year + ' at ' + hour
    if addtz:
        current_time = current_time + ' ' + time.tzname[time.daylight]
    return current_time



#------------------------------------------------------------
def runclaw(xdir='.', rundir='.', outdir='.', overwrite=False, \
             xclawcmd = 'xclaw', xclawout=None, xclawerr=None, \
             savecode=False, runmake=False):
#------------------------------------------------------------
    '''
    Run the command xdir/xclawcmd, directing the output fort.*
    files to outdir, writing unit 6 timestepping info to file xclawout.
    Runtime error messages are written to file xclawerr.
    If xclawout(xclawerr) is None, then output to stdout(stderr).

    If savecode==True, archive a copy of the code into directory outdir.

    This function returns the returncode from the process running xclawcmd,
    which will be nonzero if a runtime error occurs.
    '''

    debug = False
    if debug:
        print "<p>In runclaw... xdir = ",xdir
        print "<p>outdir = ",outdir
        print "<p>xclawout = ",xclawout
        print "<p>"
    
    startdir = os.getcwd()
    xdir = os.path.abspath(xdir)
    xclawcmd = os.path.join(xdir,xclawcmd)

    try:
	os.chdir(xdir)
    except:
        print "Cannot change to directory xdir = ",xdir
	return 999


    if runmake:
        try:
            os.system('make')
        except:
            print 'Warning: no make file in directory xdir = ',xdir

    
    if outdir != '.':
        if os.path.isfile(outdir):
            print "Error: outdir specified is a file"
            return 999
        elif (os.path.isdir(outdir) & overwrite):
            print "Directory ", outdir, " already exists, "
            print "   will be removed and recreated..."
            try:
                shutil.rmtree(outdir)
            except:
                print "Cannot remove directory ",outdir
                return 999
        elif (os.path.isdir(outdir) & (not overwrite)):
            print "Directory ", outdir, " already exists."
            print "Remove directory with 'rm -r ",outdir,"' and try again,"
            print "  or use overwrite=True in call to runclaw"
            return 999

        try:
            os.mkdir(outdir)
        except:
            print "Cannot make directory ",outdir
            return 999
        #print "Created directory ",outdir

    
    try:
        os.chdir(rundir)
    except:
        print "Cannot change to directory ",rundir
        return 999

    datafiles = glob.glob('*.data')
    if datafiles == ():
        print "Warning: no data files found in directory ",rundir
    else:
        if rundir != outdir:
            for file in datafiles:
                shutil.copy(file,os.path.join(outdir,file))

    if xclawout:
	xclawout = open(xclawout,'wb')
    if xclawerr:
	xclawerr = open(xclawerr,'wb')

    os.chdir(outdir)

    #print "\nIn directory outdir = ",outdir,"\n"

    # execute command to run fortran program:

    try:
        print "\nExecuting ",xclawcmd, "  ...  "
        pclaw = subprocess.Popen(xclawcmd,stdout=xclawout,stderr=xclawerr)
        pclaw.wait()   # wait for code to run

        if pclaw.returncode == 0:
            print "\nFinished executing\n"
        else:
            print "\n *** Runtime error: return code = %s\n " % pclaw.returncode
    except:
        print "\nError: could not execute command\n"
        return 999

    os.chdir(startdir)
    return pclaw.returncode

#---------------------------------
def svn_revision(dir="CLAW"):
#---------------------------------
    """
    Determine the svn revision number of the version of code being used.
    If dir="CLAW" it returns the revision number of the claw directory.
    This only checks the top most level and assumes this is accurate.
    """
   
    if dir=="CLAW":
        dir  = os.environ['CLAW']
    svn_entries = dir + '/.svn/entries'
    try:
        f = open(svn_entries,'r')
        lines = f.readlines()
        revision = int(lines[3])  # fourth line of file, converted to int
        f.close()
    except:
        revision = None
    
    return revision
