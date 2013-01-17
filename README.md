README 
=========

#Python 58

This is the README file for Python 58, an old but stable version.

==========

# Getting the radiative transfer code 'Python'

You can download the required structure under the structure branch.

Releases can be found under 'tags'.

Consult the wiki for how to install Python.

===========

# Basic Git Instructions


===========

# Original README file from KSL

Updated 090830 -- ksl

Note that old versions of the software are should be/are largely archived 
on the central store at

/user/long/archive_progs/Python_archive



This is a directory structure for Python, the radiative transfer code
developed by Knox Long, Christian Knigge, and lately Stuart Sim.  This
structure is supposed to allow for development as well as running Python.

To install python on a new machine, one should first copy the entire
directory structure.  This simplest way to do this (assuming you have
an account on the linux cluster at STScI) is with the rsync command:

Go to the directory where you want to create a programming environment
and simply issue the following command.

	rsync -az -e ssh --delete sockeye.stsci.edu:/data/sierr1/long/Python .


Once you have the directory structure there, you can begin to work on
the various python programs (if you are running Redhat linux).  If
you are working on a different system, e.g SunOs, however, you must 
recompile varous libraries, as follows:

Here is an outline of what needs to be done if you want to compile
python on a different operating system.  If you are working on Redhat
you should be able to skip to step 4

1.  cd software/cfitsio and rebuild cfitsio.  The install step does not
    mv the library to the correct location so mv libcfitsio.a to lib
2.  cd progs/kpar and make libkpar.a
3.  cd software.  Unpack and recompile the gsl libraries, if they are
    not on your machine already.  Make sure the binary ends up in the
    lib directory

At this point you should be able to recompile python, py_wind or any other
program you would like to be using, so:

4.  cd progs/python and make python, py_wind and any of the other programs
    you would like to be using.

The philosophy of the directory structure is as follows:
    progs: Location of working programs
    dev:   Location of developmental version of the programs.  If you tar
	a directory in progs and mv it to dev, and untar it, then all of
	your makefiles should work just as they do in progs.
    data: location for all datafiles.  Files that are mainly for reference
	should be gzipped to save space. Such files are not recreated in
       the standalone environments for running python (that you will be
       creating, especially for beowulf type environents).
    data44, etc -- These are versions of the data files needed to run older
	versions of python.
    data66: The datafiles needed to run the current verisons of the program
    bin: The location of the executables.  (It is a good idea to put
	this directory in your path)
    base_environment:  This is a base directory that is used by the routine
	MakePyEnviron in order to create a new set of directories for running
	lots of models, e.g for comparison with your favorite CV spectrum.
    Example:  A directory with a few recent examples of python runs.
    

If you are not (foolishly) trying to improve the program, but actually model
something, one uses the MakePyEnvrion to create a standalone environment
for running python.  This routine copies the binaries and the base_environment
to a directory that you designate.  (The one thing to remember is that to
use the environment the commands are bin/py rather py.  The former uses
the binary that you have trasnfered, whereas the former uses py in this
master directory.


Please send comments on this file to long@stsci.edu
