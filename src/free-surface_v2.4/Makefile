LIB   =free-surface
LNAME =lib$(LIB).a
LOBJS =free-surface.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

free-surface.o: free-surface.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c free-surface.f

clean: 
	rm *.o *.a 

