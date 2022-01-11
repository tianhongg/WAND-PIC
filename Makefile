PROGRAM.O = domain.o steps.o partition.o mesh.o multigrid.o multigridC.o wakefield.o pushtrajs.o saveload.o commute.o particles.o pushparts.o trajectory.o pulse.o cell.o namelist.o wand_PIC.o

M_MAKE = $(MAKE) -j4

FC = mpiicpc
LD = mpiicpc

Optimized = -O3

FCFLAGS = $(Optimized) $(DEFINES)
LDFLAGS = 

INC = -I$/usr/local/include
LIB = -L$/usr/local/lib -lpnetcdf

EXE = WAND

all : program

program: $(EXE)

%.o : %.cpp
		$(FC) $(FCFLAGS) -c -o $@ $< $(INC)

$(EXE): $(PROGRAM.O)
		$(LD) $(LDFLAGS) -o $@ $^ $(LIB)

float: 
	$(M_MAKE) DEFINES="${DEFINE4ALL}"

double: 
	$(M_MAKE) DEFINES="${DEFINE4ALL} -D_WTYPE"

clean:
		rm -r $(PROGRAM.O) $(EXE)
