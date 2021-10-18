LIB   =panca
LNAME =lib$(LIB).a
LOBJS =panca.o asset.o jumbo.o subs.o legge.o dislo.o
LDIR  =../../lib

$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

panca.o: panca.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c panca.f

asset.o: asset.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c asset.f

jumbo.o: jumbo.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c jumbo.f

subs.o: subs.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c subs.f

legge.o: legge.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c legge.f

dislo.o: dislo.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c dislo.f

clean: 
	rm *.o *.a 
