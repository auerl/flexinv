# 
# Makefile to make all executables required 
# to compute surface wave kernels
# 
# Ludwig Auer, 18.04.12 - 04:21 AM
# 
# Please ignore all individual .mk files in the
# Source subfolders!
#


FLAGS= -132 -static -i8 -L/usr/lib/x86_64-redhat-linux5E/lib64/
BINDIR=./bin
HRVLIB=./hrvlib
SRCDIR=./src

LIB1=./hrvlib/libio
LIB2=./hrvlib/libkern
LIB3=./hrvlib/libgeo

all: $(BINDIR)/integrate_modes \
	$(BINDIR)/min_kernel \
	$(BINDIR)/ptlyrmod \
	$(BINDIR)/crlyrmod_10 \
	$(BINDIR)/crlyrmod_20 \
	$(BINDIR)/vox2pert
	rm -f *.o

min_kernel: $(BINDIR)/min_kernel
	rm -f *.o

vox2pert: $(BINDIR)/vox2pert
	rm -f *.o

$(BINDIR)/integrate_modes: integrate_modes.o \
	gemodelib.o \
	initmodel.o \
	integrate.o \
	corfac.o \
	initmode.o \
	libio_all.o \
	radial_basis_layers.o \
	ebspl.o \
	rspln.o \
	csubs_linux.o
	ifort $(FLAGS) -o $(BINDIR)/integrate_modes \
        integrate_modes.o \
        initmodel.o \
        initmode.o \
        integrate.o \
        corfac.o \
        gemodelib.o \
	libio_all.o \
	radial_basis_layers.o \
	ebspl.o \
	rspln.o \
	csubs_linux.o

(BINDIR)/min_kernel: min_kernel.o \
	minsubs_ge.o \
	libio_all.o \
	csubs_linux.o
	ifort $(FLAGS) -o $(BINDIR)/min_kernel \
	minsubs_ge.o \
	min_kernel.o \
	libio_all.o \
	csubs_linux.o

$(BINDIR)/ptlyrmod: ptlyrmod.o libio_use.o radial_basis_layers.o ebspl.o
	ifort $(FLAGS) -o $(BINDIR)/ptlyrmod ptlyrmod.o libio_use.o \
	radial_basis_layers.o ebspl.o

$(BINDIR)/crlyrmod_10: crlyrmod_10.o remsubs.o geolib.o libio_use.o
	ifort $(FLAGS) -o $(BINDIR)/crlyrmod_10 crlyrmod_10.o remsubs.o geolib.o \
	libio_use.o

$(BINDIR)/crlyrmod_20: crlyrmod_20.o remsubs.o geolib.o libio_use.o \
	crust_20_subs.o
	ifort $(FLAGS) -o $(BINDIR)/crlyrmod_20 crlyrmod_20.o remsubs.o geolib.o \
	libio_use.o crust_20_subs.o

$(BINDIR)/vox2pert: vox2pert.o 
	ifort $(FLAGS) -o $(BINDIR)/vox2pert vox2pert.o
	rm -f *.o

vox2pert.o: $(SRCDIR)/vox2pert.f
	ifort $(FLAGS) -c -o vox2pert.o $(SRCDIR)/vox2pert.f

ptlyrmod.o: $(SRCDIR)/ptlyrmod.f
	ifort  $(FLAGS) -c -o ptlyrmod.o $(SRCDIR)/ptlyrmod.f

libio_use.o: $(HRVLIB)/libio/libio_use.f
	ifort $(FLAGS) -c -o libio_use.o \
	$(HRVLIB)/libio/libio_use.f 

radial_basis_layers.o: $(SRCDIR)/radial_basis_layers.f
	ifort $(FLAGS) -c -o radial_basis_layers.o $(SRCDIR)/radial_basis_layers.f

crlyrmod_10.o: $(SRCDIR)/crlyrmod_10.f
	ifort  $(FLAGS) -c -o crlyrmod_10.o $(SRCDIR)/crlyrmod_10.f

crlyrmod_20.o: $(SRCDIR)/crlyrmod_20.f
	ifort  $(FLAGS) -c -o crlyrmod_20.o $(SRCDIR)/crlyrmod_20.f

remsubs.o: $(SRCDIR)/remsubs.f
	ifort $(FLAGS) -c -o remsubs.o $(SRCDIR)/remsubs.f

geolib.o: $(SRCDIR)/geolib.f 
	ifort $(FLAGS) -c -o geolib.o $(SRCDIR)/geolib.f


ebspl.o: $(LIB2)/ebspl.f
	ifort $(FLAGS) -c $(LIB2)/ebspl.f -o ebspl.o

integrate_modes.o: $(SRCDIR)/integrate_modes.f $(SRCDIR)/gemodes.h
	ifort $(FLAGS) -c $(SRCDIR)/integrate_modes.f -o integrate_modes.o

initmodel.o: $(SRCDIR)/initmodel.f $(SRCDIR)/gemodes.h $(SRCDIR)/gemodl.h
	ifort $(FLAGS) -c $(SRCDIR)/initmodel.f -o initmodel.o

initmode.o: $(SRCDIR)/initmode.f $(SRCDIR)/gemodes.h $(SRCDIR)/gemodl.h
	ifort $(FLAGS) -c $(SRCDIR)/initmode.f -o initmode.o

integrate.o: $(SRCDIR)/integrate.f $(SRCDIR)/gemodl.h
	ifort $(FLAGS) -c $(SRCDIR)/integrate.f -o integrate.o

corfac.o: $(SRCDIR)/corfac.f $(SRCDIR)/gemodl.h
	ifort $(FLAGS) -c $(SRCDIR)/corfac.f -o corfac.o

gemodelib.o: $(SRCDIR)/gemodelib.f $(SRCDIR)/gemodes.h
	ifort $(FLAGS) -c $(SRCDIR)/gemodelib.f -o gemodelib.o

rspln.o: $(LIB3)/rspln.f
	ifort $(FLAGS) -c $(LIB3)/rspln.f -o rspln.o

min_kernel.o: $(SRCDIR)/min_kernel.f
	ifort $(FLAGS) -c $(SRCDIR)/min_kernel.f -o min_kernel.o

libio_all.o: $(LIB1)/libio_all.f
	ifort $(FLAGS) -c $(LIB1)/libio_all.f -o libio_all.o

crust_20_subs.o: $(SRCDIR)/crust_20_subs.f
	ifort $(FLAGS) -c -o crust_20_subs.o $(SRCDIR)/crust_20_subs.f

csubs_linux.o: $(LIB1)/csubs_linux.c
	cc -Wsystem-headers -H -c $(LIB1)/csubs_linux.c -o csubs_linux.o

minsubs_ge.o: $(SRCDIR)/minsubs_ge.f
	ifort $(FLAGS) -c $(SRCDIR)/minsubs_ge.f -o minsubs_ge.o
