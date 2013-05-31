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
CFLAGS		=-Wall -g
CXXFLAGS	=-Wall -g

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
OBJECTS=cpml.o inputChecker.o main.o fdtd.o


###########################
## rules
###########################
.PHONY:all clean

doc:formulations.pdf fdtd-1d-cpml.pdf

all: doc main

main:$(OBJECTS)
	$(CXX) -o $@ $^



$(OBJECTS):%.o:%.cpp
	$(CXX) -c $< $(CXXFLAGS)


## formulations
%.pdf:%.tex
	$(PDFLATEX) $<

clean:
	-rm -f $(TAGET) $(OBJECTS) *.log *.aux
