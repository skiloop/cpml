###########################
## compiler
###########################
CC			=gcc
CXX			=g++

###########################
## latex maker
###########################
PDFLATEX=pdflatex

###########################
## base compiler options
###########################
CFLAGS		=-Wall -g -Wno-sign-compare#-fopenmp
CXXFLAGS	=-Wall -g -Wno-sign-compare#-fopenmp

###########################
## link options
###########################
LIB=-lm -fopenmp

###########################
## virtual path
###########################
VPATH=doc:src:test


###########################
## targets
###########################
TAGET=formulations.pdf fdtd-1d-cpml.pdf

###########################
## objects
###########################
OBJECTS=main.o fdtd.o 


###########################
## rules
###########################
.PHONY:all clean

all: doc main
doc:formulations.pdf fdtd-1d-cpml.pdf


main:$(OBJECTS)
	$(CXX) -o $@ $^ $(LIB)



$(OBJECTS):%.o:%.cpp
	$(CXX) -c $< $(CXXFLAGS)


## formulations
%.pdf:%.tex
	$(PDFLATEX) $<

clean:
	-rm -f $(TAGET) $(OBJECTS) *.log *.aux main
