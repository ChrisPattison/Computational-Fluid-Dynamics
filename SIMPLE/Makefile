F90=gfortran
LINK=$(F90)
VPATH=src:bin
BPATH=bin
SPATH=src
SRCS=$(shell ls src)
OBJS=$(SRCS:.f90=.o)

INCLUDE=
LIBS=

#FLAGS=-O3 -march=native -mfpmath=sse -m64 
#FLAGS=-O3 -g -mfpmath=sse -march=native -m64 -msse4.2 -mavx -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
#FLAGS=-O3 -mfpmath=sse -march=native -Wall -Wtabs
FLAGS=-O3 -mfpmath=sse -march=native -Wall -Wtabs -m64

EXE=main

all: $(EXE)
main.o: kinds.o case.o util.o types.o Makefile
case.o: kinds.o types.o coeffutil.o
kinds.o: Makefile
util.o: kinds.o
types.o: kinds.o

$(EXE): $(OBJS) Makefile
	@echo 'Linking [$(EXE)] from [$(OBJS)] using [$(LINK)]'
	@$(LINK) $(FLAGS) -o $(EXE) $(addprefix $(BPATH)/,$(OBJS)) $(LIBS)

%.o: %.f90 Makefile
	@echo 'Compiling [$@] from [$<] using [$(F90)]'
	@$(F90) $(FLAGS) $(INCLUDE) -J $(BPATH) -c $< -o $(BPATH)/$@

clean:
	@-rm $(BPATH)/* $(EXE)
	@-touch $(BPATH)/.keep
