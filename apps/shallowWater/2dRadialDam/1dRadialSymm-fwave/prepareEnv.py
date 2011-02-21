#!/usr/bin/env python


import os 
import sys
import glob

print "\n"
print "Hello, Word! This is my first python script. It calls few times 'make' with some options to compile the 'internal' SharpClaw libraries necessary to run the current test case! \n"


# Get current directory
cwd = os.getcwd()

# Clean current directory
print "CLEAN CURRENT DIRECTORY \n"
os.system("make clean")
os.system("make clobber")

# Create libraries paths needed to clean-up the `internal' SharpClaw subroutines
sharpLibDir = os.path.expandvars('$SCLAW/lib/')
sharpLibDir1D = os.path.expandvars('$SCLAW/lib/1d')
sharpLibDir2D = os.path.expandvars('$SCLAW/lib/2d')

# Clean `internal' SharpClaw subroutine
os.chdir(sharpLibDir1D)
os.system("make clean")

os.chdir(sharpLibDir2D)
os.system("make clean")

# INSERT THE PROBLEM SPACE DIMENSION
print "\n"
print "IS THE TEST CASE 1D or 2D? "
dim = raw_input("Type 1 or 2: ")

if dim.lower() == '1':
    # Compile the 1D subroutines
    os.chdir(sharpLibDir1D)
    os.system("make")
elif dim.lower() == '2':
    # Compile the 2D subroutines
    os.chdir(sharpLibDir2D)
    os.system("make")
else:
    print "Aborting."
    sys.exit()








