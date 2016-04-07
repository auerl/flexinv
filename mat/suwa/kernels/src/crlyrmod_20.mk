FLAGS= -132 -static -i8
HRVLIB= ../HRVLIB


../BIN/crlyrmod_20: crlyrmod_20.o crlyrmod_20.mk remsubs.o geolib.o libio_use.o \
	crust_20_subs.o
	ifort $(FLAGS) -o ../BIN/crlyrmod_20 crlyrmod_20.o remsubs.o geolib.o \
	libio_use.o crust_20_subs.o
	rm -f *.o

crlyrmod_20.o: crlyrmod_20.f
	ifort  $(FLAGS) -c -o crlyrmod_20.o crlyrmod_20.f

remsubs.o: remsubs.f
	ifort $(FLAGS) -c -o remsubs.o remsubs.f

geolib.o: geolib.f 
	ifort $(FLAGS) -c -o geolib.o geolib.f

libio_use.o: $(HRVLIB)/libio/libio_use.f
	ifort $(FLAGS) -c -o libio_use.o \
	$(HRVLIB)/libio/libio_use.f 

crust_20_subs.o: crust_20_subs.f
	ifort $(FLAGS) -c -o crust_20_subs.o crust_20_subs.f
