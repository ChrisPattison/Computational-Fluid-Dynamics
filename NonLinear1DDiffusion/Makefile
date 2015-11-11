F90=gfortran
LINK=$(F90)
INCLUDE=$(shell pkg-config --cflags plplotd-f95)
LIBS=$(shell pkg-config --libs plplotd-f95)
VPATH=src:bin
BPATH=bin
SPATH=src

SRCS=$(shell ls src)
OBJS=$(SRCS:.f90=.o)

FLAGS=-O2 -march=native -mfpmath=sse -m64 
#FLAGS=-g -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -fdiagnostics-color=never

EXE=main

all: $(EXE)
main.o: plplotlib.o kinds.o case.o Makefile
plplotlib.o: kinds.o Makefile
case.o: kinds.o
kinds.o: Makefile

$(EXE): $(OBJS) Makefile
	@echo 'Linking [$(EXE)] from [$(OBJS)] using [$(LINK)]'
	@$(LINK) $(FLAGS) -o $(EXE) $(addprefix $(BPATH)/,$(OBJS)) $(LIBS)

%.o: %.f90 Makefile
	@echo 'Compiling [$@] from [$<] using [$(F90)]'
	@$(F90) $(FLAGS) $(INCLUDE) -J $(BPATH) -c $< -o $(BPATH)/$@

clean:
	@-rm $(BPATH)/* $(EXE)
	@-touch $(BPATH)/.keep
