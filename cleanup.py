 #
# cleanup script for claw
#

import os,sys

print """
    WARNING: This will do 'make veryclean' or 'make clean'
             in all directories with this capacity. 

             This will remove *.o, fort.*, frame*, *.html *.o, etc.  
             files from all subdirectories.

             Any previously computed results in WILL BE LOST! 
              (except in directories with 'SAVE' in their name)

     """
rootdir = os.path.abspath('.')
sys.stdout.write('  OK to cleanup all subdirectories of %s?  ' \
                % rootdir)
answer = raw_input()
if answer not in ['y','Y','yes']:
    print '  Aborting!'
    sys.exit(1)


for (dirpath, subdirs, files) in os.walk(rootdir):
    currentdir = os.path.abspath(os.getcwd())
    cleandir = os.path.abspath(dirpath)
    if (cleandir.find('SAVE') == -1) & (cleandir.find('svn') == -1):
        # skip over directories if 'SAVE' is in the full path name
        # skip over directories if 'svn' is in the full path name
        os.chdir(cleandir)
        if 'Makefile' in files:
            print "Cleaning up ", dirpath
            os.system('make clobber')
#            M = open('Makefile', 'r')
#	    if M.read().find('veryclean') > -1:
#	        os.system('make veryclean')
#	    elif M.read().find('clean') > -1:
#	        os.system('make clean')
        os.chdir(currentdir)
    

