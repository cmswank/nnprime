#simple makefile for mirror neutron trajectories
#Use g++ compiler 
CC=g++
#CFLAGS are used to make .o files. 
CFLAGS=-c
#Please export the Boost library location. Or put write it in here.  
IDIRS=-I${BOOST_INCLUDE}
#using library of math, gsl and fftw, (I think the first two are default, but who knows.)
LDINC=-lm -lboost_thread -lboost_system #-lgsl #-lfftw3
#Required source code for objects
SOURCES= particleTrajectory.cpp  mirrorConductor.cpp
#location of main(). 
MAIN=runMirror.cpp
#The executable
EXE=runMirror
OBJECTS:=$(SOURCES:.cpp=.o)

all:	$(EXE) $(OBJECTS)

default: $(EXE) $(OBJECTS)

#maybe you want to only do objects? 
objects: $(OBJECTS)

$(EXE): $(OBJECTS)
	$(CC) $(IDIRS) $(OBJECTS) $(MAIN) -o $@ $(LDINC)
#Objects don't need each other as of yet...
$(OBJECTS):
	$(CC) $(CFLAGS) $(IDIRS) $(SOURCES) $(LDINC)

clean:
	rm *.o $(EXE)