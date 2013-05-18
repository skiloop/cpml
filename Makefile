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
VPATH=doc:src


###########################
## targets
###########################
TAGET=formulations.pdf

###########################
## objects
###########################
OBJECTS=cpml.o


###########################
## rules
###########################
.PHONY:all clean

all: formulations.pdf cpml.o



$(OBJECTS):%.o:%.cpp
	$(CXX) -c $< $(CXXFLAGS)


## formulations
formulations.pdf:formulations.tex
	$(PDFLATEX) $<

clean:
	-rm -f $(TAGET) $(OBJECTS) *.log *.aux
