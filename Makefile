PROGRAM.O = domain.o steps.o partition.o mesh.o multigrid.o multigridC.o wakefield.o pushtrajs.o saveload.o commute.o particles.o pushparts.o trajectory.o pulse.o cell.o namelist.o wand_PIC.o

FC = mpiicpc
LD = mpiicpc

FCFLAGS = -O3
LDFLAGS = -O3

INC = -I$(TACC_PNETCDF_INC)
LIB = -L$(TACC_PNETCDF_LIB) -lpnetcdf

EXE = WAND

all : program

program: $(EXE)

$(EXE): $(PROGRAM.O)
		$(LD) $(LDFLAGS) -o $@ $^ $(LIB)
%.o : %.cpp
		$(FC) $(FCFLAGS) -c -o $@ $< $(INC)
clean:
		rm -r $(PROGRAM.O) $(EXE)
