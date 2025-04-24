# This is the makefile for diablo.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#
# As the user, set the following definitions:

#**********   User Defined Variables Below *********
# Fortran 90 complier to use:
COMPILER = h5pfc
#COMPILER = mpif90
#COMPILER = ifort

# Faster code
#COMPOPTS = -xAVX -ipo -O3 -no-prec-div -opt-prefetch
COMPOPTS = 

# Any desired options for the compiler (e.g. -O2, -g, etc.)
#USEROPTS = -O3 -ftree-vectorize
#USEROPTS = -align all -O3 -xHost -fpp
#USEROPTS = -g -O3 -fno-align-commons -cpp
USEROPTS = -fno-align-commons -g -fallow-argument-mismatch -cpp -O0 -fcheck=all

# Location where fftw and netcdf (optional) libraries are installed
#LINKDIR = /usr/local/lib
#LINKDIR = /usr/local/include
LINKDIR = /home/teaves001/my_fftw/lib

# Location where the optional netcdf include file (netcdf.inc) is installed
INCLUDEDIR = /usr/local/include

# Option to compile with MPI libraries
PARALLEL = TRUE
HDF5 = TRUE

# Option to enable the LES model (loads the required variables into memory)
LES = FALSE

# Option to compile with the NetCDF libraries
NETCDF = FALSE

# Option to run different flavors (basic, ensemble, etc.)
ENSEM = FALSE
BATCH = FALSE
# **********    END of user definitions ************

# Use the parameters to set flags
ifeq ($(LES),TRUE)
LES_o = les.o
else
LES_o = no_les.o
endif

ifeq ($(PARALLEL),TRUE)
COMPILER = mpif90
MPI = mpi.o
ifeq ($(HDF5),TRUE)
HDF5_o = hdf5.o
COMPILER = h5pfc
HDF5OPTS=-DHDF5
endif
else
MPI = mpi_serial.o
endif

MAIN = diablo.f
HEADER = header
ENSEM_HOOKS = dummy_code/ensem_dummy.f
BATCH_HOOKS = dummy_code/batch_dummy.f
HOOKS = batch_hooks.o ensem_hooks.o
ADJOINT = 

ifeq ($(ENSEM),TRUE)
MAIN = ensemble.f
HEADER = header header_ensem
COMPILER = mpif90
ENSEM_HOOKS = ensem_hooks.f
endif

ifeq ($(BATCH),TRUE)
MAIN = batch.f
HEADER = header header_batch
BATCH_HOOKS = batch_hooks.f
#ADJOINT = adj_chan.o adj_per.o
ADJOINT = adj_per.o
endif

ifeq ($(NETCDF),TRUE)
COMPOPTS = $(USEROPTS) $(HDF5OPTS) -I$(INCLUDEDIR)
LINKOPTS = -L$(LINKDIR) -lrfftw -lfftw -lnetcdf
NETCDF_o = netcdf.o
else
COMPOPTS = $(USEROPTS) $(HDF5OPTS)
LINKOPTS = -L$(LINKDIR) -lrfftw -lfftw -lm
#LINKOPTS = -L$(LINKDIR) -lfftw3 -lm
NETCDF_o = no_netcdf.o
endif


diablo: $(MAIN) diablo_io.o channel.o $(LES_o) $(NETCDF_o) \
	fft.o rand.o $(HOOKS) $(ADJOINT) $(MPI) \
	$(HEADER) grid_def $(HDF5_o)
	$(COMPILER) $(COMPOPTS) $(MAIN) -o diablo \
	diablo_io.o channel.o $(LES_o) $(NETCDF_o) \
	fft.o rand.o $(HOOKS) $(ADJOINT) \
	$(MPI) $(LINKOPTS) $(HDF5_o)

diablo_io.o: diablo_io.f header grid_def
	$(COMPILER) $(COMPOPTS) -c diablo_io.f

channel.o: channel.f fft.o $(MPI) header grid_def
	$(COMPILER) $(COMPOPTS) -c channel.f

ifeq ($(LES),TRUE) 
les.o: les.f fft.o header header_les grid_def
	$(COMPILER) $(COMPOPTS) -c les.f
else
no_les.o: dummy_code/no_les.f
	$(COMPILER) $(COMPOPTS) -c dummy_code/no_les.f
endif

ifeq ($(NETCDF),TRUE)
netcdf.o: netcdf.f header grid_def
	$(COMPILER) $(COMPOPTS) -c netcdf.f
else
no_netcdf.o: dummy_code/no_netcdf.f 
	$(COMPILER) $(COMPOPTS) -c dummy_code/no_netcdf.f
endif

ifeq ($(PARALLEL),TRUE)
mpi.o: mpi.f header header_mpi grid_def
	$(COMPILER) $(COMPOPTS) -c mpi.f
else
mpi_serial.o: dummy_code/mpi_serial.f header header_mpi grid_def
	$(COMPILER) $(COMPOPTS) -c dummy_code/mpi_serial.f
endif

hdf5.o : hdf5.f
	$(COMPILER) $(COMPOPTS) -c hdf5.f    

ensem_hooks.o: $(ENSEM_HOOKS) header header_ensem grid_def
	$(COMPILER) $(COMPOPTS) -c $(ENSEM_HOOKS) -o ensem_hooks.o

batch_hooks.o: $(BATCH_HOOKS) header header_batch grid_def
	$(COMPILER) $(COMPOPTS) -c $(BATCH_HOOKS) -o batch_hooks.o

ifeq ($(BATCH),TRUE)
#adj_chan.o: adj_chan.f header header_batch grid_def
#	$(COMPILER) $(COMPOPTS) -c adj_chan.f

adj_per.o: adj_chan.f header header_batch grid_def
	$(COMPILER) $(COMPOPTS) -c adj_per.f
endif

fft.o:  fft.f header grid_def
	$(COMPILER) $(COMPOPTS) -lrfftw -lfftw -lm -c fft.f

rand.o:  rand.f
	$(COMPILER) $(COMPOPTS) -c rand.f

clean:
	rm -f *.o fort.* *~ diablo core

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
# 
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add 
# trailing underscores to the name of every fftw function where they
# appear throughout the code.

