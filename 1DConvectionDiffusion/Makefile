F90=gfortran
LINK=$(F90)
INCLUDE=$(shell pkg-config --cflags plplotd-f95)
LIBS=$(shell pkg-config --libs plplotd-f95)
VPATH=src:bin
BPATH=bin
SPATH=src

SRCS=$(shell ls src)
OBJS=$(SRCS:.f90=.o)

ifeq ($(OS),Windows_NT)
INCLUDE="C:/Program Files (x86)/plplot/include/"
INCLUDE="C:/Program Files (x86)/plplot/lib/fortran/modules/plplot"
LIBS="C:/Program Files (x86)/plplot/lib/fortran/modules/plplot"
endif

#FLAGS=-O3 -march=native -mfpmath=sse -m64 
#FLAGS=-O2 -g -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -fdiagnostics-color=never
#FLAGS=-Ofast -g -mfpmath=sse -march=native -Wall -Wtabs -fbacktrace -fdiagnostics-color=never
#=-Ofast -g -mfpmath=sse -march=native -Wall -Wtabs -fcheck=all -fbacktrace -fdiagnostics-color=never
FLAGS=-O2 -g -mfpmath=sse -march=native -Wall -Wtabs -fcheck=all -fbacktrace
FLAGS=-O1 -g -Wall -Wtabs -fcheck=all -fbacktrace

EXE=main

all: $(EXE)
main.o: plplotlib.o kinds.o case.o util.o list.o Makefile
plplotlib.o: kinds.o Makefile
case.o: kinds.o util.o
kinds.o: Makefile
util.o: kinds.o
list.o: kinds.o

$(EXE): $(OBJS) Makefile
	@echo 'Linking [$(EXE)] from [$(OBJS)] using [$(LINK)]'
	@$(LINK) $(FLAGS) -o $(EXE) $(addprefix $(BPATH)/,$(OBJS)) $(LIBS)

%.o: %.f90 Makefile
	@echo 'Compiling [$@] from [$<] using [$(F90)]'
	@$(F90) $(FLAGS) $(INCLUDE) -J $(BPATH) -c $< -o $(BPATH)/$@

clean:
	@-rm $(BPATH)/* $(EXE)
	@-touch $(BPATH)/.keep
