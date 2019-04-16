FC = gfortran
CXX = g++

UNAME := $(shell uname -s)

FCFLAGS= 
CCFLAGS=-I $(EIGEN_PATH) -I $(HODLR_PATH)/include -DUSE_DOUBLE
LDFLAGS=-lstdc++

all: test.x
test.o : kernel_mod.o

%.x : %.o kernel_mod.o kernel_capi.o Kernel.o
	${FC} $^ -o $@ ${LDFLAGS}

%.o : %.f90
	${FC} ${FCFLAGS} -c $< -o $@

%.o : %.cpp
	${CXX} ${CCFLAGS} -c $^ -o $@

.PHONY : clean

clean :
	${RM} -rf *.o *.mod test.x
