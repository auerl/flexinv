FFLAGS = -assume byterecl -extend_source
BINDIR = .
LIBDIR = ./gulib
LIB = libkernel

all:	lib matrix_bw_vx
	rm -f *.o

lib:
	cd ./gulib/libkernel/; make -f library_vox.mk

matrix_bw_vx: $(BINDIR)/matrix_bw_vx

$(BINDIR)/matrix_bw_vx:  $(LIBDIR)/libkernel.a \
	matrix_bw_vx.o
	ifort $(FFLAGS) -o $(BINDIR)/matrix_bw_vx matrix_bw_vx.o \
	$(LIBDIR)/libkernel.a

matrix_bw_vx.o: matrix_bw_vx.f
	ifort -c $(FFLAGS) -o matrix_bw_vx.o matrix_bw_vx.f
