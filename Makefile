# Linux Intel compiler
#FC       = ifort
#FCOPT    = -O5
##FCOPT    = -O0 -g -traceback -check all -check bounds -debug all 
#FCOPTOMP = -qopenmp
#LIB      = -mkl

#Linux GNU compiler
FC       = gfortran   
FCOPT    = -O3
#FCOPT    = -O0 -g -fbacktrace -fcheck=all
FCOPTOMP = -fopenmp
LIB      = -llapack -lblas

export FC
export FCOPT
export FCOPTOMP
export LIB

all:
	cd src/analysis/         ; make; cd ../../
	cd src/free-surface_v2.4/; make; cd ../../
	cd src/free-warp/        ; make; cd ../../
	cd src/panca/            ; make; cd ../../
	cd src/warp_v1.5omp/     ; make; cd ../../
	cd src/pressure/         ; make; cd ../../
	cd src/tatra/            ; make; cd ../../
	cd src/geomod/		 ; make; cd ../../
	cd src/metamodels/       ; make; cd ../../
	$(FC) $(FCOPT) $(FCOPTOMP) -o bin/L2-Sea \
				src/main.f90 \
				lib/libanalysis.a \
				lib/libfree-warp.a \
				lib/libpanca.a \
				lib/libfree-surface.a \
				lib/libwarp.a \
				lib/libpressure.a \
				lib/libtatra.a \
				lib/libgeomod.a \
				lib/libmetamodels.a \
				$(LIB)

clean:
	rm src/*/*.o
	rm lib/*.a
	rm bin/*


