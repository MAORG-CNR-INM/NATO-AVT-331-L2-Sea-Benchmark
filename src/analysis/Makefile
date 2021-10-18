LIB   =analysis
LNAME =lib$(LIB).a
LOBJS =analysis.o 
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

analysis.o: analysis.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c analysis.f90


