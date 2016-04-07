FLAGS= -132 -static -i8
HRVLIB=../HRVLIB
MODES=../MODES

../BIN/ptlyrmod: ptlyrmod.o ptlyrmod.mk libio_use.o radial_basis_layers.o ebspl.o
	ifort $(FLAGS) -o ../BIN/ptlyrmod ptlyrmod.o libio_use.o \
	radial_basis_layers.o ebspl.o
	rm -f *.o

ptlyrmod.o: ptlyrmod.f
	ifort  $(FLAGS) -c -o ptlyrmod.o ptlyrmod.f

libio_use.o: $(HRVLIB)/libio/libio_use.f
	ifort $(FLAGS) -c -o libio_use.o \
	$(HRVLIB)/libio/libio_use.f 

radial_basis_layers.o: $(MODES)/radial_basis_layers.f
	ifort $(FLAGS) -c -o radial_basis_layers.o $(MODES)/radial_basis_layers.f

ebspl.o: $(HRVLIB)/libkern/ebspl.f
	ifort $(FLAGS) -c -o ebspl.o $(HRVLIB)/libkern/ebspl.f

