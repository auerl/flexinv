F90=ifort
FLAGS=-assume byterecl -i8 -O2 -mcmodel=medium -openmp 
BINDIR=./bin
SRCDIR=./src

all:	pro
	rm -f *.o *.mod

pro: 	 project_adpx_module.mod $(BINDIR)/project $(BINDIR)/addhtc


$(BINDIR)/project: project_adpx_main.o project_adpx_module.o param_adpx_angle_depth.o project_matrix.o  
	$(F90) $(FLAGS) project_adpx_module.o \
	param_adpx_angle_depth.o \
        project_matrix.o \
        project_adpx_main.o -o $(BINDIR)/project

$(BINDIR)/addhtc: addhtc.o
	$(F90) -o $(BINDIR)/addhtc addhtc.o

project_adpx_module.mod:
	$(F90)  $(FLAGS) -c $(SRCDIR)/project_adpx_module.f90 

param_adpx_angle_depth.o: $(SRCDIR)/param_adpx_angle_depth.f90
	$(F90) $(FLAGS) -c $(SRCDIR)/param_adpx_angle_depth.f90 -o param_adpx_angle_depth.o

project_matrix.o: $(SRCDIR)/project_matrix.f90
	$(F90) $(FLAGS) -c $(SRCDIR)/project_matrix.f90 -o project_matrix.o

project_adpx_main.o: $(SRCDIR)/project_adpx_main.f90
	$(F90) $(FLAGS) -c $(SRCDIR)/project_adpx_main.f90 -o project_adpx_main.o

addhtc.o: $(SRCDIR)/addhtc.f90
	$(F90) -assume byterecl -i8 -O2 -c $(SRCDIR)/addhtc.f90 -o addhtc.o
