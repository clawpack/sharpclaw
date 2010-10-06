import os

rootdir = '.'

for (dirpath, subdirs, files) in os.walk(rootdir):
    currentdir = os.path.abspath(os.getcwd())
    os.chdir(os.path.abspath(dirpath))
    os.system('rm -f *.o')
    os.system('rm -f fort.*')
    os.system('rm -f xwclaw')
    os.system('rm -f xsclaw')
    os.system('rm -f xclaw')

    os.chdir(currentdir)
