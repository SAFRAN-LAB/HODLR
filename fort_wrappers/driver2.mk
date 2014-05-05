CC	=g++
FC 	=gfortran
CFLAGS	=-c -Wall -O2 -ffast-math -ffinite-math-only -I ../header/ -I ~/Dropbox/Eigen/
FFLAGS  =-c -Wall -O3
LDFLAGS	=-lstdc++
CSOURCES	= fort_hodlr_wrappers.cpp Matrix_Entry_Routine.cpp
FSOURCES 	= driver2.f
COBJECTS	=$(CSOURCES:.cpp=.o)
FOBJECTS	=$(FSOURCES:.f=.o)
EXECUTABLE	= DRIVER2
all: $(FSOURCES) $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(FOBJECTS) $(COBJECTS)
	$(FC) $(FOBJECTS) $(COBJECTS) -o $@ $(LDFLAGS) 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.f.o:
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -rf $(COBJECTS) $(FOBJECTS) $(EXECUTABLE)
