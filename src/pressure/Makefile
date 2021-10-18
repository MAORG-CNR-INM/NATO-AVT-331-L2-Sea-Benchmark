LIB   =pressure
LNAME =lib$(LIB).a
LOBJS =pressure.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

pressure.o: pressure.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c pressure.f

