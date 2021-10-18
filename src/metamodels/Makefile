LIB   =metamodels
LNAME =lib$(LIB).a
LOBJS =metamodels.o 
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

metamodels.o: metamodels.f90
	$(FC) $(FCOPT) $(FCOPTOMP) -c metamodels.f90



