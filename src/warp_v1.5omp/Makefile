#-----------------------------------------------------------------------------#
#                    ####               M A K E F I L E
#                  ###                        per       
#                ##                         W A R P
#               _|__                 WAve Resistance Program
#   ____________|  |____
#----\________________/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#-----------------------------------------------------------------------------#
#FCOPT =

LIB   =warp
LNAME =lib$(LIB).a
LDIR  =../../lib

LOBJS = warp.o caldino.o shippo.o ittc.o hemiqua.o hemitri.o panlib.o \
         spln1.o indomo.o riama.o calve.o hemiquall.o hemitrill.o \
         primaqua.o primatri.o spln2.o fidell.o jellqua.o jelltri.o \
         splynqua.o splyntri.o transom3.o tracut.o amado.o \
         recipe.o tecout.o spline.o simqit.o flotta.o


$(LNAME):$(LOBJS)
	ar -vru $(LDIR)/$(LNAME) $(LOBJS)
	ranlib  $(LDIR)/$(LNAME)

warp.o: warp.f warp.cmn
	$(FC) $(FCOPT) $(FCOPTOMP) -c warp.f

hemiqua.o: hemiqua.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c hemiqua.f 

hemitri.o: hemitri.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c hemitri.f 

hemiquall.o: hemiquall.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c hemiquall.f 

hemitrill.o: hemitrill.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c hemitrill.f 

spline.o: spline.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c spline.f

spln1.o: spln1.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c spln1.f

spln2.o: spln2.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c spln2.f

caldino.o: caldino.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c caldino.f

fidell.o: fidell.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c fidell.f

jellqua.o: jellqua.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c jellqua.f

jelltri.o: jelltri.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c jelltri.f

calve.o: calve.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c calve.f

panlib.o: panlib.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c panlib.f

indomo.o: indomo.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c indomo.f

splynqua.o: splynqua.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c splynqua.f

splyntri.o: splyntri.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c splyntri.f

primaqua.o: primaqua.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c primaqua.f

primatri.o: primatri.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c primatri.f

transom3.o: transom3.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c transom3.f

tracut.o: tracut.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c tracut.f

amado.o: amado.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c amado.f

riama.o: riama.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c riama.f

shippo.o: shippo.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c shippo.f

simqit.o: simqit.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c simqit.f

ittc.o: ittc.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c ittc.f

tecout.o: tecout.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c tecout.f

flotta.o: flotta.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c flotta.f

glotta.o: glotta.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c glotta.f

recipe.o: recipe.f
	$(FC) $(FCOPT) $(FCOPTOMP) -c recipe.f

clean: 
	rm *.o *.a 

