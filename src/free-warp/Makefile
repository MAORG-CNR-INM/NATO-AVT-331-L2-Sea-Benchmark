LIB   =free-warp
LNAME =lib$(LIB).a
LOBJS =free-warp.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

free-warp.o: free-warp.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c free-warp.f

clean: 
	rm *.o *.a 

