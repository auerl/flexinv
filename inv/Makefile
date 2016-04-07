FLAGS= -i8 -O3 -shared-intel -mcmodel=large -assume byterecl -132
FLAGS90= -i8 -O3 -shared-intel -mcmodel=large -assume byterecl -132

BINDIR=./bin
SRCDIR=./src
HRVLIB=../mat/suwa/kernels/hrvlib
LIB=

all: $(BINDIR)/adotm $(BINDIR)/flexinv90 \
	$(BINDIR)/flexinv_module.mod	
	rm -f *.o *.mod

$(BINDIR)/flexinv_module.mod:
	ifort  $(FLAGS90) -c $(SRCDIR)/flexinv_module.f90 

$(BINDIR)/flexinv90: flexinv_module.o flexinv90.o flexinv_subroutines.o
	ifort $(FLAGS90) -o $(BINDIR)/flexinv90 $(LIB) \
	flexinv90.o \
	flexinv_module.o \
	flexinv_subroutines.o

$(BINDIR)/adotm: adotm.o
	ifort -132 -O3 -shared-intel -mcmodel=large -o $(BINDIR)/adotm adotm.o

flexinv90.o: $(SRCDIR)/flexinv90.f90
	ifort $(FLAGS90) -c $(SRCDIR)/flexinv90.f90 -o flexinv90.o
flexinv_module.o: $(SRCDIR)/flexinv_module.f90
	ifort $(FLAGS90) -c $(SRCDIR)/flexinv_module.f90 -o flexinv_module.o
flexinv_subroutines.o: $(SRCDIR)/flexinv_subroutines.f90
	ifort $(FLAGS90) -c $(SRCDIR)/flexinv_subroutines.f90 -o flexinv_subroutines.o
adotm.o: $(SRCDIR)/adotm.f
	ifort -132 -O3 -shared-intel -mcmodel=large -c $(SRCDIR)/adotm.f -o adotm.o
