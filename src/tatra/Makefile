LIB   =tatra
LNAME =lib$(LIB).a
LOBJS =tatra.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

tatra.o: tatra.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c tatra.f

