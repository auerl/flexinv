FLAGS= -assume byterecl -extend_source -shared-intel -mcmodel=large #-132 -i8 -fPIC -shared-intel -mcmodel=large

SWDIR=./kernels

matrix_sw_vx: Makefile matrix_sw_vx.o premgephvelo.o libio_use.o rsple.o rspln.o
	ifort $(FLAGS) -o matrix_sw_vx \
	matrix_sw_vx.o \
	premgephvelo.o \
	rsple.o \
	rspln.o \
	libio_use.o
	rm -f *.o

matrix_sw_vx.o: matrix_sw_vx.f
	ifort $(FLAGS) -c matrix_sw_vx.f -o matrix_sw_vx.o

premgephvelo.o: $(SWDIR)/src/premgephvelo.f
	ifort $(FLAGS) -c \
	$(SWDIR)/src/premgephvelo.f \
	-o premgephvelo.o

libio_use.o: $(SWDIR)/hrvlib/libio/libio_use.f
	ifort $(FLAGS) -c -o libio_use.o \
	$(SWDIR)/hrvlib/libio/libio_use.f

rsple.o: $(SWDIR)/hrvlib/libgeo/rsple.f
	ifort $(FLAGS) -c -o rsple.o \
	$(SWDIR)/hrvlib/libgeo/rsple.f

rspln.o: $(SWDIR)/hrvlib/libgeo/rspln.f
	ifort $(FLAGS) -c -o rspln.o \
	$(SWDIR)/hrvlib/libgeo/rspln.f

