FLAGS= -132 -static -i8

../BIN/vox2pert: vox2pert.o vox2pert.mk param.o
	ifort $(FLAGS) -o ../BIN/vox2pert vox2pert.o param.o
	rm -f *.o

vox2pert.o: vox2pert.f
	ifort $(FLAGS) -c -o vox2pert.o vox2pert.f

param.o: ../../../../LIB/param.f
	ifort $(FLAGS) -c -o param.o ../../../../LIB/param.f

