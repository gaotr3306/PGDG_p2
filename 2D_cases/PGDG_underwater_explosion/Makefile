################################################################################
# Build script for MPM program
################################################################################

#-----------------------------------------------------------------------------
# Self-defined flags
#-----------------------------------------------------------------------------

#*****************************************************#
#                 General setups                      #
#*****************************************************#

PARALLEL = SINGLECPU

PERCISION = DOUBLE
ifeq ($(PERCISION),DOUBLE)
DFLAGS	+= -D_DOUBLE_
else
DFLAGS	+= -D_SINGLE_
endif

# DFLAGS	+= -D_G_


KERNEL = -D_WINDLAND_C4_
# KERNEL = -D_CUBIC_SPLINE_
#KERNEL = -D_GAUSSIAN_
#KERNEL = -D_QUINTIC_SPLINE_

DFLAGS  += $(KERNEL)

#*****************************************************#
#           Physical and partition dimension          #
#*****************************************************#

DIMENSION = DIM2
ifeq ($(DIMENSION),DIM1)
DFLAGS	+= -D_DIM1_
else ifeq ($(DIMENSION),DIM2)
DFLAGS	+= -D_DIM2_
else
DFLAGS	+= -D_DIM3_
endif



#*****************************************************#
#                     Compiler                        #
#*****************************************************#

CXX := mpicxx
EXECUTABLE	:= MPMRUN
# C++ source files (compiled with gcc / c++)

CCFILES	:= mpm.cpp cell_list.cpp main.cpp solver.cpp system.cpp level_infor.cpp particle.cpp 
CCFILES += particle_flow.cpp mpm_ini.cpp 


#-----------------------------------------------------------------------------
# BOOST MPI & TBB COMPILE LIB
#-----------------------------------------------------------------------------
BOOST_PREFIX =/home/tianrun/Downloads/software_gao/boost_1_80_0/install
BOOST_MPI_INCLUDES := -I${BOOST_PREFIX}/include 
BOOST_MPI_FLAGS   := -L${BOOST_PREFIX}/lib -lboost_mpi -lboost_serialization -Wl,-rpath, -L${BOOST_PREFIX}/lib

TBB_PREFIX= /home/tianrun/Downloads/software_gao/oneTBB-2022.0.0/install

TBB_INCLUDE:=-I${TBB_PREFIX}/include 
TBB_MPI_FLAGS:=-L${TBB_PREFIX}/lib 

EIGEN_INCLUDE:=-I/home/tianrun/Downloads/software_gao/eigen-3.4.0

#-----------------------------------------------------------------------------
#  TBB TECPLOT and common libs
#-----------------------------------------------------------------------------

CCFLAGS   += -O3
LIBPARALLEL = -ltbb -ltbbmalloc


#----------------------------------------------------------------------------
# PREPARE AND COMPILE:
#-----------------------------------------------------------------------------

CXXFLAGS += $(CCFLAGS)

# Common includes and LD_flags
INCLUDES      := -I. -I.. $(BOOST_MPI_INCLUDES) ${TBB_INCLUDE} ${EIGEN_INCLUDE} 
LDFLAGS       := $(LIBPARALLEL) -lrt $(LIBS) ${BOOST_MPI_FLAGS} ${TBB_MPI_FLAGS}  

# Location of output and obj files
OUTDATADIR	?= ./outdata
TPOUTDATADIR	?= ./tpout
OBJECTDIR 	?= ./obj
RESTARTDIR	?= ./rstfile

# Set up object files
 OBJECTS 	 =  $(patsubst %.cpp,$(OBJECTDIR)/%.cpp.o,$(notdir $(CCFILES)))

# Target rules
all: directories $(EXECUTABLE)

# Link commands:

$(OBJECTDIR)/%.cpp.o : %.cpp
	$(CXX) $(CXXFLAGS) -O  $(DFLAGS) $(INCLUDES)  -o $@ -c $<
# -fopenmp 
$(EXECUTABLE): directories $(OBJECTS) Makefile
	$(CXX) $(CCFLAGS) $(DFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LDFLAGS)

directories:
	mkdir -p $(OBJECTDIR)
	mkdir -p $(TPOUTDATADIR)
	mkdir -p $(OUTDATADIR)
	mkdir -p $(RESTARTDIR)


run: all
	 ./MPMRUN >& run.log

clean:
	rm -rf $(EXECUTABLE) $(OBJECTDIR)/* *.bin

cleandata:
	rm -rf $(OUTDATADIR)/* $(TPOUTDATADIR)/* 

