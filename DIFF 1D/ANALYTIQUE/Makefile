FC = gfortran
# Avec options de debuggage activees
OPT = -g -O0 -fbounds-check 
# Sans option de debuggage
#OPT = -g

OBJ = m_type.o prog.o sousprog.o VTSWriter.o
EXE = prog.exe

prog:	$(OBJ)
	$(FC) $(OPT) $(OBJ) -o $(EXE)

m_type.o:	m_type.f90
	$(FC) $(OPT) m_type.f90 -c	

prog.o :	prog.f90
	$(FC) $(OPT) prog.f90 -c

sousprog.o :	sousprog.f90
	$(FC) $(OPT) sousprog.f90 -c

VTSWriter.o :	VTSWriter.f90
	$(FC) $(OPT) VTSWriter.f90 -c

clean :
	/bin/rm -f $(OBJ) $(EXE) $(MODS) *.mod

