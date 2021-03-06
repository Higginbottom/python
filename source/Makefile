# This is the makefile for the python related programs
#
# usage      make [CC=...] [D] python
#
# Adding D causes the routine to be run in a way that profiling and ddd can be used.
# Otherwise the run will be optimized to run as fast as possible. CC is an option to choose
# a different compiler other than mpicc.
#
# History
# 05jan	ksl	54f  Modified Makefile so that the Version number is automatically
# 		copied to version.h  Also fixed so that one does not need to put
# 		a new filename in two separate places, and so that the latest
# 		compile version automatically becomve shte 
# 05apr	ksl	55c   Modified in a way that implies that one must have gsl as part
# 		of the Python directory tree.  This is because I was having problems
# 		with gsl after we went to Redhat Enterprise at the Institute, and
# 		so that Stuart and I could standardise on the distribution.
# 08jul	ksl	Removed pfop from routines so no need to complile with g77
# 13jun jm      Added capability to switch to use debugger routne
#
# 13jun jm/ss	SS added parallel flag -DMPI_ON and mpicc wrapper for gcc
#		is now used as C compiler
# 13jul jm	can now compiled with gcc using 'make CC=gcc python'
# 13jul jm	kpar is now integrated into python and compiled here


#MPICC is now default compiler- currently code will not compile with gcc
CC = mpicc
# CC = gcc	can use GCC either from command line or by uncommenting this
FC = g77
# FC = gfortran
# speciify any extra compiler flags here
EXTRA_FLAGS =


#Check a load of compiler options
#this is mostly to address issue $100
ifeq (mpicc, $(CC))		
	# you're using mpicc, so we'll compile with the DMPI_ON flag
	MPI_FLAG = -DMPI_ON

	# check what underlying compiler mpi is using, and the version
	# we'll use this to print to the user and also to do some checks
	MPI_COMPILER := $(shell mpicc --showme:command)
	COMPILER_VERSION := $(shell expr `'$(CC)' -dumpversion`)
	COMPILER_PRINT_STRING = Compiling with mpicc wrapper, for compiler $(MPI_COMPILER) $(COMPILER_VERSION)

	# if it's gcc we want to check if the version is 4.8 or later
	# if it is then we'll disable aggressive loop optimizations, see #100
	ifeq (gcc, $(MPI_COMPILER))
		GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1-2 -d.` \>= 4.8)
		ifeq ("$(GCCVERSIONGTEQ4)", "1")
			EXTRA_FLAGS += -fno-aggressive-loop-optimizations 	# add the flag to EXTRA_FLAGS
			COMPILER_PRINT_STRING += with -fno-aggressive-loop-optimizations
		endif
	endif	

else ifeq (gcc, $(CC))
	# no mpicc,
	MPI_FLAG =

	# check the version
	# we'll use this to print to the user and also to do some checks
	COMPILER_VERSION := $(shell expr `gcc -dumpversion`)

	# if it's gcc we want to check if the version is 4.8 or later
	# if it is then we'll disable aggressive loop optimizations, see #100
	GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1-2 -d.` \>= 4.8)
	COMPILER_PRINT_STRING = Compiling with $(CC) $(COMPILER_VERSION)
	ifeq ("$(GCCVERSIONGTEQ4)", "1")
		EXTRA_FLAGS += -fno-aggressive-loop-optimizations
		COMPILER_PRINT_STRING += with -fno-aggressive-loop-optimizations
	endif

else	# you must be using clang or icc
	MPI_FLAG =

	# check the version we'll use this to print to the user
	COMPILER_VERSION = $(shell expr `$(CC) -dumpversion`)
	COMPILER_PRINT_STRING = Compiling with $(CC) $(COMPILER_VERSION)
endif

# this command finds out how many files with uncommitted changes there are
GIT_DIFF_STATUS := $(shell expr `git status --porcelain 2>/dev/null| grep "^ M" | wc -l`)
GIT_COMMIT_HASH := $(shell expr `git rev-parse HEAD`)

# Test if indent is installed. We use a different version on mac and linux.
OS_TYPE := $(shell uname -s)
ifeq ($(OS_TYPE), Darwin)
	INDENT_PATH := $(shell which gnuindent)
	INDENT_VERSION = $(shell gnuindent -version)
else
	INDENT_PATH := $(shell which indent)
	INDENT_VERSION := $(shell indent -version)
endif

ifeq (, $(INDENT_PATH))
	#If indent is not installed, error
	INDENT_STRING = Gnuindent not installed so not indenting code
	INDENT_CMD = 
else
	ifeq (GNU,$(firstword $(INDENT_VERSION)))
		INDENT_STRING = Indenting code with $(INDENT_VERSION)
		# this command runs indent with GNU styling
		INDENT_CMD = $(INDENT_PATH) -gnu -l140 -bli0 -nut *.c *.h
	else
		# if it is not GNU indent then do not run
		INDENT_STRING = $(INDENT_PATH) is not the GNU version
		INDENT_CMD = 
	endif
endif


INCLUDE = ../include
LIB = ../lib
BIN = ../bin

ifeq (D,$(firstword $(MAKECMDGOALS)))
# use pg when you want to use gprof the profiler
# to use profiler make with arguments "make D python" 
# this can be altered to whatever is best	
	CFLAGS = -g -pg -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
	FFLAGS = -g -pg   
	PRINT_VAR = DEBUGGING, -g -pg -Wall flags
else
# Use this for large runs
	CFLAGS = -O3 -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
	FFLAGS =         
	PRINT_VAR = LARGE RUNS, -03 -Wall flags
endif



# next line for debugging when concerned about memory problems and duma installed in python directory
# LDFLAGS= -L$(LIB)  -lm -lkpar  -lgsl -lgslcblas ../duma_2_5_3/libduma.a -lpthread
# next line if you want to use kpar as a library, rather than as source below
# LDFLAGS= -L$(LIB)  -lm -lkpar -lcfitsio -lgsl -lgslcblas 

LDFLAGS= -L$(LIB) -lm -lgsl -lgslcblas 

#Note that version should be a single string without spaces. 



VERSION = 81c_dev

CHOICE=1             // Compress plasma as much as possible
# CHOICE=0           //  Keep relation between plasma and wind identical

startup:
	@echo $(COMPILER_PRINT_STRING)			# prints out compiler information
	@echo 'YOU ARE COMPILING FOR' $(PRINT_VAR)	# tells user if compiling for optimized or debug
	@echo 'MPI_FLAG=' $(MPI_FLAG)
	echo "#define VERSION " \"$(VERSION)\" > version.h
	echo "#define CHOICE"   $(CHOICE) >> version.h
	echo "#define GIT_COMMIT_HASH" \"$(GIT_COMMIT_HASH)\" >> version.h
	echo "#define GIT_DIFF_STATUS" $(GIT_DIFF_STATUS)\ >> version.h

indent:
	@echo $(INDENT_STRING)
	$(INDENT_CMD)


foo: foo.o signal.o time.o
	$(CC) ${cfllags} foo.o signal.o time.o ${LDFLAGS}  -o foo


# these are the objects required for compiltion of python
# note that the kpar_source is now separate from this
python_objects = bb.o get_atomicdata.o photon2d.o photon_gen.o \
		saha.o spectra.o wind2d.o wind.o  vvector.o debug.o recipes.o \
		trans_phot.o phot_util.o resonate.o radiation.o \
		wind_updates2d.o windsave.o extract.o pdf.o roche.o random.o \
		stellar_wind.o homologous.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o  continuum.o get_models.o emission.o recomb.o diag.o \
		sv.o ionization.o  ispy.o   levels.o gradv.o reposition.o \
		anisowind.o util.o density.o  detail.o bands.o time.o \
		matom.o estimators.o wind_sum.o yso.o elvis.o cylindrical.o rtheta.o spherical.o  \
		cylind_var.o bilinear.o gridwind.o partition.o signal.o auger_ionization.o \
		agn.o shell_wind.o compton.o torus.o zeta.o dielectronic.o \
		spectral_estimators.o variable_temperature.o matom_diag.o \
		log.o lineio.o rdpar.o direct_ion.o pi_rates.o matrix_ion.o para_update.o \
		setup.o photo_gen_matom.o macro_gov.o \
		reverb.o paths.o setup2.o run.o brem.o search_light.o synonyms.o
		


# For reasons that are unclear to me.  get_models.c cannot be included in the sources
# Problems ocurr due to the prototypes that are generated.  ksl 160705
python_source= bb.c get_atomicdata.c python.c photon2d.c photon_gen.c \
		saha.c spectra.c wind2d.c wind.c  vvector.c debug.c recipes.c \
		trans_phot.c phot_util.c resonate.c radiation.c \
		wind_updates2d.c windsave.c extract.c pdf.c roche.c random.c \
		stellar_wind.c homologous.c hydro_import.c corona.c knigge.c  disk.c\
		lines.c  continuum.c emission.c recomb.c diag.c \
		sv.c ionization.c  ispy.c  levels.c gradv.c reposition.c \
		anisowind.c util.c density.c  detail.c bands.c time.c \
		matom.c estimators.c wind_sum.c yso.c elvis.c cylindrical.c rtheta.c spherical.c  \
		cylind_var.c bilinear.c gridwind.c partition.c signal.c auger_ionization.c \
		agn.c shell_wind.c compton.c torus.c zeta.c dielectronic.c \
		spectral_estimators.c variable_temperature.c matom_diag.c \
		direct_ion.c pi_rates.c matrix_ion.c para_update.c setup.c \
		photo_gen_matom.c macro_gov.c \
		reverb.c paths.c setup2.c run.c brem.c search_light.c synonyms.c

# kpar_source is now declared seaprately from python_source so that the file log.h 
# can be made using cproto
kpar_source = lineio.c rdpar.c log.c synonyms.c

additional_py_wind_source = py_wind_sub.c py_wind_ion.c py_wind_write.c py_wind_macro.c py_wind.c 

prototypes: 
	cp templates.h templates.h.old
	cproto -I$(INCLUDE) $(python_source) ${additional_py_wind_source} > foo.h  
	cp foo.h templates.h
	cproto -I$(INCLUDE) $(kpar_source) > log.h 

python: startup  python.o $(python_objects)
	$(CC)  ${CFLAGS} python.o $(python_objects) $(kpar_objects) $(LDFLAGS) -o python
		cp $@ $(BIN)/py
		mv $@ $(BIN)/py$(VERSION)

#This line is jsut so you can use make D python for debugging
D:	
	@echo 'Debugging Mode'

py_wind_objects = py_wind.o get_atomicdata.o py_wind_sub.o windsave.o py_wind_ion.o \
		emission.o recomb.o util.o detail.o \
		pdf.o random.o recipes.o saha.o \
		stellar_wind.o homologous.o sv.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o vvector.o wind2d.o wind.o  ionization.o  py_wind_write.o levels.o \
		radiation.o gradv.o phot_util.o anisowind.o resonate.o density.o \
		matom.o estimators.o yso.o elvis.o photon2d.o cylindrical.o rtheta.o spherical.o  \
		cylind_var.o bilinear.o gridwind.o py_wind_macro.o partition.o auger_ionization.o\
		spectral_estimators.o shell_wind.o compton.o torus.o zeta.o dielectronic.o \
        variable_temperature.o bb.o rdpar.o log.o direct_ion.o diag.o matrix_ion.o \
		pi_rates.o photo_gen_matom.o macro_gov.o \
		time.o reverb.o paths.o synonyms.o



py_wind: startup $(py_wind_objects)
	$(CC) $(CFLAGS) $(py_wind_objects) $(LDFLAGS) -o py_wind
	cp $@ $(BIN)
	mv $@ $(BIN)/py_wind$(VERSION)


summarize_atomic: summarize_atomic.o get_atomicdata.o log.o rdpar.o synonyms.o  
	$(CC) $(CFLAGS) summarize_atomic.o get_atomicdata.o log.o rdpar.o  synonyms.o  -o summarize_atomic
	mv $@ $(BIN)


table_objects = windsave2table.o get_atomicdata.o py_wind_sub.o windsave.o py_wind_ion.o \
		emission.o recomb.o util.o detail.o \
		pdf.o random.o recipes.o saha.o \
		stellar_wind.o homologous.o sv.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o vvector.o wind2d.o wind.o  ionization.o  py_wind_write.o levels.o \
		radiation.o gradv.o phot_util.o anisowind.o resonate.o density.o \
		matom.o estimators.o yso.o elvis.o photon2d.o cylindrical.o rtheta.o spherical.o  \
		cylind_var.o bilinear.o gridwind.o py_wind_macro.o partition.o auger_ionization.o\
		spectral_estimators.o shell_wind.o compton.o torus.o zeta.o dielectronic.o \
        	variable_temperature.o bb.o rdpar.o log.o direct_ion.o diag.o matrix_ion.o \
		pi_rates.o photo_gen_matom.o macro_gov.o reverb.o paths.o time.o synonyms.o    



windsave2table: $(table_objects) 
	$(CC) $(CFLAGS) $(table_objects) $(LDFLAGS) -o windsave2table
	cp $@ $(BIN)
	mv $@ $(BIN)/windsave2table$(VERSION)

py_smooth: py_smooth.o 
	$(CC) $(CFLAGS) py_smooth.o  $(LDFLAGS)  -o py_smooth
		mv $@ $(BIN)

test_bb: bb.o test_bb.o pdf.o recipes.o bilinear.o time.o 
	$(CC)  ${CFLAGS} bb.o pdf.o test_bb.o  recipes.o bilinear.o time.o $(LDFLAGS) -o test_bb
	
test_pow: test_pow.o pdf.o recipes.o bilinear.o time.o 
	$(CC)  ${CFLAGS} pdf.o test_pow.o  recipes.o bilinear.o time.o $(LDFLAGS) -o test_pow

test_matrix: test_matrix.o $(python_objects)
	$(CC) ${CFLAGS} test_matrix.o $(python_objects) $(LDFLAGS) -o test_matrix
		mv $@ $(BIN)

test_dielectronic: test_dielectronic.o $(python_objects)
	$(CC) ${CFLAGS} test_dielectronic.o $(python_objects) $(LDFLAGS) -o test_dielectronic
		mv $@ $(BIN)

t_bilinear:  t_bilinear.o bilinear.o  
	$(CC) $(CFLAGS) t_bilinear.o  bilinear.o  $(LDFLAGS)  -o t_bilinear

	
plot_roche: plot_roche.o roche.o vvector.o phot_util.o recipes.o 
	${CC} ${CFLAGS} plot_roche.o roche.o vvector.o phot_util.o recipes.o  \
		$(LDFLAGS) -o plot_roche
		mv $@ $(BIN)/plot_roche



# This diagnostic routine has not been kept up to date
py_grid: bb.o get_atomicdata.o py_grid.o photon2d.o photon_gen.o \
		saha.o spectra.o wind2d.o wind.o  vvector.o debug.o recipes.o \
		trans_phot.o phot_util.o resonate.o radiation.o \
		wind_updates2d.o windsave.o extract.o pdf.o roche.o random.o \
		stellar_wind.o homologous.o sv.o hydro_import.o corona.o knigge.o  disk.o\
		lines.o continuum.o get_models.o emission.o recomb.o util.o anisowind.o \
		sv.o ionization.o  ispy.o   levels.o gradv.o reposition.o \
		yso.o elvis.o cylindrical.o rtheta.o \
		matom.o estimators.o wind_sum.o \
		density.o bands.o detail.o 
	$(CC)  ${CFLAGS} py_grid.o get_atomicdata.o photon2d.o photon_gen.o \
		bb.o\
		saha.o spectra.o wind2d.o wind.o  vvector.o debug.o recipes.o \
		trans_phot.o phot_util.o resonate.o radiation.o \
		stellar_wind.o homologous.o  corona.o util.o anisowind.o \
		wind_updates2d.o windsave.o extract.o ispy.o bands.o\
		pdf.o roche.o random.o continuum.o get_models.o  \
		lines.o ionization.o emission.o  recomb.o reposition.o \
		sv.o hydro_import.o knigge.o disk.o  levels.o gradv.o density.o \
		yso.o elvis.o cylindrical.o rtheta.o \
		matom.o estimators.o wind_sum.o \
		detail.o \
		$(LDFLAGS) -o py_grid
		cp $@ $(BIN)/py_grid
		mv $@ $(BIN)/py_grid$(VERSION)


FILE = get_atomicdata.o atomic.o

libatomic.a:  get_atomicdata.o atomic.o
	ar ru libatomic.a get_atomicdata.o atomic.o
	ranlib libatomic.a
	mv libatomic.a $(LIB)
	cp atomic.h  $(INCLUDE)


saha_sources = saha_inv.c get_atomicdata.c

startup_saha: startup $(saha_sources)
	cproto -I$(INCLUDE) $(saha_sources) > saha_templates.h

saha_inv: saha_inv.o get_atomicdata.c
	$(CC) ${CFLAGS} saha_inv.o get_atomicdata.o
	cp $@ $(BIN)/saha_inv
	mv $@ $(BIN)/saha_inv$(VERSION)


clean :
	rm -f *.o  *~ 
