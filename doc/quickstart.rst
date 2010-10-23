=================================
Quick Start Guide
=================================

Obtaining SharpClaw
=======================

The current development version of SharpClaw can be obtained via
Mercurial::

    hg clone http://bitbucket.org/ketch/sharpclaw

If you don't have Mercurial, you can download it from
http://mercurial.selenic.com/downloads/.

Alternatively, the most recent release version of SharpClaw can
be downloaded from http://web.kaust.edu.sa/faculty/davidketcheson/Software.html.

Installing SharpClaw
=========================

You must set the environment variable `SCLAW` to the path where you
put the `sharpclaw` directory. In csh the relevant line in your .cshrc file
is::

    setenv SCLAW /path/to/sharpclaw/

Then, to build the 1D code::

    cd $SCLAW/lib/1d
    make

Or to build the 2D code::

    cd $SCLAW/lib/2d
    make

If you have built the code for 1D and later need to rebuild for 2D,
(or vice versa), it is necessary to `make clean` in the `/1d/` directory.

Examples
=========

To run an example, just cd to any subdirectory of `examples/` and type::

    make all

This will compile the executable, run the code, and create plots.
Finally, you will be given a message directing you to open a webpage
containing the plots.
