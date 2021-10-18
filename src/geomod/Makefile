LIB   =geomod
LNAME =lib$(LIB).a
LOBJS =orthopatches3D.o morphing.o regrid.o geodiff.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

orthopatches3D.o: orthopatches3D.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c orthopatches3D.f90

morphing.o: morphing.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c morphing.f90

regrid.o: regrid.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c regrid.f

geodiff.o: geodiff.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c geodiff.f90

