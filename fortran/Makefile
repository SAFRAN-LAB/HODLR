FC = gfortran
CXX = g++

UNAME := $(shell uname -s)

FCFLAGS= 
CCFLAGS=
LDFLAGS=-lstdc++

all: test.x
test.o : foo_mod.o

%.x : %.o foo_mod.o foo_capi.o Foo.o
	${FC} $^ -o $@ ${LDFLAGS}

%.o : %.f90
	${FC} ${FCFLAGS} -c $< -o $@

%.o : %.cpp
	${CXX} ${CCFLAGS} -c $^ -o $@

.PHONY : clean

clean :
	${RM} -rf *.o *.mod test.x
