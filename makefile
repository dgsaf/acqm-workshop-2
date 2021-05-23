# here the make accepts the compiler, gpu, optimisation level, profiling
COMPILERTYPE ?= GCC
GPUTYPE ?= NVIDIA
OPTLEVEL ?= 2
PROFILER ?= ON

# optmisation flags
OPTFLAGS = -O$(OPTLEVEL)
# profiling flags if desired
ifeq ($(PROFILER), ON)
	OPTFLAGS += -pg -g
endif

# formatting characters for info output
NULL :=
TAB := $(NULL)  $(NULL)

# compilers definitions
GCC = gcc
GCCCXX = g++
GCCFORT = gfortran
CLANG = clang
CLANGCXX = clang++
CLANGFORT = flang
AOMP = aompcc
AOMPCXX = aompcc
CRAYCC = cc
CRAYCXX = CC
CRAYFORT = ftn
INTELCC = icc
INTELCXX = icpc
INTELFORT = ifort

# C flags. Currently only really need to enforce c11 standard
GCCCFLAGS = -std=c11

# fortran flags
GCCFFLAGS = -cpp -ffixed-line-length-none -dM \
	-Wall -Wextra -Wconversion -pedantic
INTELFFLAGS = -cpp -extend-source -D_INTELFTN

# cuda
NCC = nvcc
NCXX = nvcc
CUDA_FLAGS = -DUSECUDA

# hip
HCC = hipcc
HCXX = hipcc
HIP_FLAGS = -DUSEHIP
OCL_FLAGS = -lOpenCL -DUSEOPENCL

# mpi
MPICC = mpicc
MPICXX = mpicxx
MPIFORT = mpif90

# openmp
GCCOMP_FLAGS = -fopenmp -DUSEOPENMP
GCCOMPTARGET_FLAGS = -fopenmp -DUSEOPENMPTARGET
INTELOMP_FLAGS = -qopenmp -DUSEOPENMP
INTELMPTARGET_FLAGS = -qopenmp -DUSEOPENMPTARGET
GCCOACC_FLAGS = -fopenacc -fopt-info-optimized-omp -DUSEOPENACC
INTELOACC_FLAGS = -qopenacc -fopt-info-optimized-omp -DUSEOPENACC

# set default compilers
CC = $(GCC)
CXX = $(GCCCXX)
FORT = $(GCCFORT)
GPUCC = $(NCC)
GPUCXX = $(NCXX)
FFLAGS = $(GCCFFLAGS)
OMP_FLAGS = $(GCCOMP_FLAGS)
OMPTARGET_FLAGS = $(GCCOMPTARGET_FLAGS)
OACC_FLAGS = $(GCCOACC_FLAGS)
CFLAGS = $(GCCCFLAGS)

# change if required
ifeq ($(COMPILERTYPE), CLANG)
	CC = $(CLANG)
	CXX = $(CLANGCXX)
	FORT = $(CLANGFORT)
endif
ifeq ($(COMPILERTYPE), AOMP)
	CC = $(AOMP)
	CXX = $(AOMPCXX)
endif
ifeq ($(COMPILERTYPE), CRAY)
	CC = $(CRAYCC)
	CXX = $(CRAYCXX)
	FORT = $(CRAYFORT)
	FFLAGS = -eZ -ffree
endif
ifeq ($(COMPILERTYPE), INTEL)
	CC = $(INTELCC)
	CXX = $(INTELCXX)
	FORT = $(INTELFORT)
	FFLAGS = $(INTELFFLAGS)
	OMP_FLAGS = $(INTELOMP_FLAGS)
	OMPTARGET_FLAGS = $(INTELOMPTARGET_FLAGS)
endif
ifeq ($(GPUTYPE), AMD)
	GPUCC = $(HCC)
	GPUCCXX = $(HCXX)
endif

# compiler used for openmp
OMPCC = $(CC)
OMPCXX = $(CXX)
OMPFORT = $(FORT)

# compilers used for openacc
OACCCC = $(CC)
OACCCXX = $(CXX)
OACCFORT = $(FORT)

# compilers used for opencl
OCLC = $(CCGPU)
OCLCXX = $(CXXGPU)

COMMONFLAGS = $(OPTFLAGS)

#
.PHONY : dirs allinfo configinfo buildinfo makecommands clean
.PHONY : homonuclear_diatomic

#
all : dirs buildinfo homonuclear_diatomic

#
allinfo : configinfo buildinfo makecommands

# information about current configuration
configinfo :
	$(info Compiler options:)
	$(info > Compiler can be selected with COMPILERTYPE=(GCC|CLANG|AOMP|CRAY))
	$(info > GPU compiler can be selected with GPUTYPE=(NVIDIA|AMD))
	$(info > Optimisation level can be selected with OPTLEVEL=(0|1|2|3|fast))
	$(info > Profiling can be turned on/off with PROFILER=(ON|OFF))
	$(info )

# information about current build given the commands make was passed
buildinfo :
	$(info Current compilers selected:)
	$(info > Compiling with ($(CC)|$(CXX)|$(FORT)) for CPU focused codes)
	$(info > Compiling with ($(MPICC)|$(MPICXX)|$(MPIFORT)) \
	for MPI-CPU focused codes)
	$(info > Compiling with ($(GPUCC)|$(GPUCXX)) for GPU focused codes)
	$(info > Compiling with ($(OMPCC)|$(OMPCXX)|$(OMPFORT)) \
	for OpenMP directive GPU focused codes)
	$(info > Compiling with ($(OACCCC)|$(OACCCXX)|$(OACCFORT)) \
	for OpenACC directive GPU focused codes)
	$(info )

# information about current build given the commands make was passed
makecommands :
	$(info Make commands:)
	$(info > Make is configured so that the following can be compiled \
	if provided this argument:)
	$(info > )
	$(info > homonuclear_diatomic : compiles FORTRAN program \
	homonuclear_diatomic which calculates one-electron homonuclear-diatomic \
	molecule)
	$(info > $(TAB)sources : homonuclear_diatomic.f90)
	$(info )

#
dirs :
	[ -d obj ] || mkdir obj
	[ -d mod ] || mkdir mod
	[ -d bin ] || mkdir bin

#
clean :
	rm obj/*
	rm mod/*
	rm bin/*

#
homonuclear_diatomic : bin/homonuclear_diatomic

#
obj/rsg.o : src/rsg.f
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/rsg.f -o obj/rsg.o

#
obj/intp.o : src/intp.f
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/intp.f -o obj/intp.o

#
obj/wigner.o : src/wigner.f
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/wigner.f -o obj/wigner.o

#
obj/laguerre.o : src/laguerre.f90
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/laguerre.f90 -o obj/laguerre.o \
	-J mod/

#
obj/io.o : src/io.f90
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/io.f90 -o obj/io.o \
	-J mod/

#
bin/homonuclear_diatomic : src/homonuclear_diatomic.f90 \
	obj/rsg.o obj/laguerre.o obj/io.o

	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/homonuclear_diatomic.f90 \
	-o obj/homonuclear_diatomic.o -I mod/

	$(FORT) $(COMMONFLAGS) $(FFLAGS) -o bin/homonuclear_diatomic \
	obj/homonuclear_diatomic.o obj/rsg.o obj/laguerre.o obj/io.o
