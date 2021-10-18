c
c                    ####  
c                  ###    
c                ##                       S W A R P
c               _|__    INSEAN   -   Ship WAve Resistance Program
c   ____________|  |____
c~~~~\________________/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
c
c
c-----------------------------------------------------------------------
c
c   
c   --> Calcolo della Resistenza d'onda, della Resistenza totale 
c       e del campo ondoso  generato da una carena che avanza con
c       velocita' costante in acque calme.
c       Il problema e' formulato in un sistema di coordinate x,y,z, solidali
c       alla carena.
c   --> La Rw viene calcolata integrando la pressione  sullo scafo
c       immerso secondo progetto. La sub. TRACUT fornisce i rilievi per
c       l'analisi con la tecnica del taglio d'onda trasversale (TATRA).
c   --> La Rt viene calcolata ipotizzando valide, per tutti i pannelli,
c       le eq. di strato limite di lastra piana.
c   --> La formulazione e' linearizzata : se nekel = 0 e' del tipo Dawson, 
c                                             altrimenti e' Neuman-Kelvin.
c   --> La derivata seconda viene calcolata per via analitica.
c   --> Il file fr.inp contiene il numero di simulazioni richieste e la 
c                               sequenza dei valori del numero di Froude
c   --> La carena puo' essere di tipo SWATH
c   --> L'eventuale poppa di tipo TRANSOM deve avere x=cost
c   --> In caso di poppa di tipo TRANSOM viene aggiunto il termine di 
c       resistenza idrostatica
c-----------------------------------------------------------------------
c       PARAMETER da controllare
c   --> ATTENZIONE ! La parameter (nti=100) (Num. MAX di punti transom)
c                    in transom3.f
c   --> ATTENZIONE ! La dimensione max di x_0,y_0 in hydroloads
c   --> ATTENZIONE ! La dimensione max di x_0,y_0,yappo,... in elepre
c----------------------------------------------------------------------*

      subroutine warp(ntraint)

c-----------------------------------------------------------------------

c      implicit none

      include "warp.cmn"

!      integer, parameter :: stdout=6

      integer     ntraint,ntraest,nlo,nlovallec,ntrapoc
      integer     ntra,nplib

c-----------------------------------------------------------------------

      integer itriqua(ntotd)
      integer indx(ntotd)
      integer iwork(ntotd)
      real x1(ntotd),y1(ntotd),z1(ntotd),
     &     x2(ntotd),y2(ntotd),z2(ntotd),
     &     x3(ntotd),y3(ntotd),z3(ntotd),
     &     x4(ntotd),y4(ntotd),z4(ntotd)
      real area(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real xns(ntotd),yns(ntotd),zns(ntotd)
      real amat(ntotd,ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real vx(ntotd),vy(ntotd),vz(ntotd)
      real xh(ntotd),yh(ntotd),zh(ntotd)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd),
     &     fi1x(ntotd),fi1y(ntotd),fi1z(ntotd)
      real vll(ntotd),fi0ll
      real pres0(ncard),pres1(ncard),pres12(ncard),presto(ncard)
      real sigma1(ntotd),sigma0(ntotd),eta(nplibd),r(ntotd)
      real fr
      real appo(ntotd),work(4*ntotd)
      real fextz,mexty

      character*80 carena
      character*20 nfile
      character*1 cata,drycho,trans

c-----------------------------------------------------------------------

      real actsink,acttrim,alfa,alpha,atrim
      real Cf,Cf_ittc,Cw0,Cyd
      real Cmy,Ct,Ct_ittc,Cw,Cfy,Cfz,Cfz0
      real diana,dislo,draft,dryfr,dryfr2,dx,eps
      real frn,fx0
      real ro,nu,gi
      real dislo0,xG0,zG0

      integer iaso,isolv,igeomod,ngrid,igrid,iunit,numfilest,numfilint

      namelist/FLUID_PARAMETERS/ ro,nu,gi
      namelist/MAIN_PARAMETERS/ iaso,isolv,igeomod,ngrid,igrid,
     &                          numfilest,numfilint
      namelist/WARP_PARAMETERS/ diana,kbase,modo,nekel,drycho,dryfr,zG0

      iunit  = 24

      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=MAIN_PARAMETERS)
      close(iunit)
      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=FLUID_PARAMETERS)
      close(iunit)
      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=WARP_PARAMETERS)
      close(iunit)

      fi0x = 0.
      fi0y = 0.
      fi0z = 0.
      fi1x = 0.
      fi1y = 0.
      fi1z = 0.
      vx   = 0.
      vy   = 0.
      vz   = 0.
      r    = 0.
      sigma0=0.
      sigma1=0.

c      write(*,*) drycho,dryfr

      if(drycho.eq.'n') then 
        dryfr = 0.
      end if

c      write(*,*) dryfr

c----------------------------------------------------------------------*
c     Dati di ingresso ed in uscita
c----------------------------------------------------------------------*
      if (ntraint.gt.0.and.ntrapoc.gt.0) cata='y'
      pig    = 2.0*2.0*atan(1.0)
      eps    = 1.0 
 1    eps    = 0.5 * eps
      upe    = 1.0 + eps
      if (upe.ne.1.0) go to 1
      eps    = eps*2.0
c
      open (unit=24,file='panca.aux',status='old')
        read(24,*)
        read(24,*)
        read(24,*)
        read(24,*) xG
      close(24)
c
c      if(modo.eq.(-1)) nekel = 1
c
      open (unit=24,file='fr.inp',status='old')
        read(24,*) nfr
c        do i=1,nfr
c          read(24,*) fr(i),xG0,dislo0,fextz(i),mexty(i)
          read(24,*) fr,xG0,dislo0,fextz,mexty
c        end do
      close(24)

      open(unit=nun,file='warp.out',status='unknown')
        rewind(nun)
      open(unit=nct,file='Ct.out',status='unknown',
     &              form='formatted',recl=64000)
!        write(nct,111)
        write(nct,*) '# Fr Cw Ch Cf Cf_ittc Ct Ct_ittc Sink Trim WSurf'
 111    format(4x,'Fr',10x,'Cw',10x,'Cyd',10x,'Cf',10x,'Cf_ittc',
     &         10x,'Ct',10x,'Ct_ittc',10x,'snkg',10x,'trim')
c----------------------------------------------------------------------*
c    Geometria della carena: legge CONFINE.GRD
c                            SUP = sup. TOTALE carena (anche parte y<0)
c---------------------------------------------------------------------*
      ycamax=0.0
c

c      open(20,file='cosn.dbg',status='unknown',recl=10000)

      open(unit=22,file='confine.grd',status='old',
     &             form='formatted',recl=64000)
      open(unit=33,file='confine.plt',status='unknown',
     &             form='formatted',recl=64000)
      write(33,*) 'VARIABLES = "X" "Y" "Z"'
      write(33,*) 'ZONE'
      rewind(22)
      read(22,100) carena
      read(22,*) ncar,lontra,idum,idum,idum,
     &           sup,shiplen,actsink,acttrim,
     &           draft,dryfr2,swl,wli,dislo
      if (dryfr.eq.0.) dryfr=dryfr2
      !write(*,*) 'DRYFR  - ', dryfr,' - '
      if (ncar.gt.ncard) then
c         write(*,*)' DISACCORDO TRA WARP E CONFINE.GRD'
c         write(*,*)' ncar(confine.grd) =',ncar
c         write(*,*)' ncar(warp.f)      =',ncard
         stop
      end if
      if (ncar.eq.0) then
c         write(*,*)' ERRORE - NON C''E'' LA CARENA!!! '
         stop
      end if
c
c      write(*,*) 'ncar=',ncar
      do i=1,ncar
        read(22,*) x1(i),x2(i),x3(i),x4(i),xn(i),xh(i),cosxn(i)
c        write(*,*) i
c        write(20,*) cosxn(i)
      end do
      do i=1,ncar
        read(22,*) y1(i),y2(i),y3(i),y4(i),yn(i),yh(i),cosyn(i)
        ycamax = max(ycamax,y1(i),y2(i),y3(i),y4(i))
c        write(*,*) i
c        write(20,*) cosyn(i)
      end do
      do i=1,ncar
        read(22,*) z1(i),z2(i),z3(i),z4(i),zn(i),zh(i),coszn(i)
c        write(*,*) i
c        write(20,*) coszn(i)
      end do
      do i=1,ncar
        read(22,*) area(i),itriqua(i)
        write(33,222) xn(i),yn(i),zn(i)
c        write(*,*) i
      end do
      close(33)
      close(22)

c      close(20)

c      write(*,*) 'dryfr=',dryfr

 100  format(1x,a)
c----------------------------------------------------------------------*
c     Geometria della superficie libera: legge FSGRID.GRD
c----------------------------------------------------------------------*

      call panlib(ycamax,
     &            xn,yn,zn,x1,x2,x3,x4,xh,y1,y2,y3,y4,yh,
     &            z1,z2,z3,z4,zh,cosxn,cosyn,coszn,
     &            area,ncar,alfa,ncard,nplib,nun,
     &            ntra,ntraest,ntraint,ntrapoc,nlovallec,nlo,
     &            ntotd,nplibd)
      nulib  = ntra*nlo     !!  Numero di pannelli sulla sup lib
      nlast  = ncar+nulib   !!  L'ultimo pannello prima della coda
      ntot   = ncar+nplib

      if(nplib.gt.nplibd) then
        write(*,*) 'STOP!!!'
        write(*,*) 'Errore nel dimensionamento'
        write(*,*) 'NPLIB  = ',nplib
        write(*,*) 'NPLIBD = ',nplibd
      end if

      if(ncar.gt.ncard) then
        write(*,*) 'STOP!!!'
        write(*,*) 'Errore nel dimensionamento'
        write(*,*) 'NCAR  = ',ncar
        write(*,*) 'NCARD = ',ncard
      end if

c----------------------------------------------------------------------*
c      call system('rm -rf fsgrid.grd')
c----------------------------------------------------------------------*
c     Stampa dati 
c----------------------------------------------------------------------*
      call uscita(nfr,diana,dryfr,nekel,ncar,nplib,ntra,
     &            nlo,ntraint,ntraest,ntrapoc,nlovallec,
     &            modo,lontra,cata,nun)
c-----------------------------------------------------------------------
c     Tutti i pannelli della sup.lib. sono quadrilateri !
c-----------------------------------------------------------------------
      do i=ncar+1,ntot
        itriqua(i) = 4
      end do
c----------------------------------------------------------------------*
c     I punti di controllo della sup.lib. sono traslati in avanti
c----------------------------------------------------------------------*
      call shippo(xn,yn,zn,xns,yns,zns,
     &            ncar,nplib,ntra,nlast,
     &            ncard,ntotd)
c----------------------------------------------------------------------*
c     Plotta la griglia shiftata
c----------------------------------------------------------------------*

c      call glotta(xns,yns,zns,ncar,nulib,nplib,
c     &            ntra,ntraint,ntrapoc,nlo,nlovallec)

c----------------------------------------------------------------------*
c     Calcolo delle sigma0  --- doppio modello
c----------------------------------------------------------------------*

      call amado(x1,x2,x3,x4,xn,xh,y1,y2,y3,y4,yn,yh,
     &           z1,z2,z3,z4,zn,zh,cosxn,cosyn,coszn,
     &           xns,yns,zns,
     &           area,fi0x,fi0y,fi0z,vx,vy,vz,
     &           amat,sigma0,indx,appo,diana,nekel,
     &           ncar,nlast,eps,itriqua,modo,
     &           nplib,ntra)

c----------------------------------------------------------------------*
c     Calcolo pressione e resistenza del doppio modello
c     I risultati si riferiscono alla carena COMPLETA
c     (anche la parte con y<0)
c----------------------------------------------------------------------*

      fx0=0.0
      fz0=0.0
      do i=1,ncar
        pres0(i)=-0.5*(fi0x(i)**2+fi0y(i)**2+fi0z(i)**2+2.0*fi0x(i))
        fx0=fx0-pres0(i)*cosxn(i)*area(i)
        fz0=fz0-pres0(i)*coszn(i)*area(i)
      end do

c      open(unit=99,file='psuni-dm.dat',status='unknown',
c     &             form='formatted',recl=64000)
c      write(99,*) ncar
c      do i=1,ncar
c        vx0 = 1.+fi0x(i)
c        vy0 =    fi0y(i)
c        vz0 =    fi0z(i)
c        write (99,222) xns(i),yns(i),zns(i),vx0,vy0,vz0,pres0(i)
c      end do
c      close(99)
c-----------------------------------------------------------------------
c Fx0 e Fz0 sono le componenti dinamiche secondo x e z.  In particolare,
c Fz0 non e' nulla solo perche' si riferisce a meta' (z<0) carena.
c Vengono raddoppiate per tener conto della parte simmetrica (y<0).
c-----------------------------------------------------------------------
      fx0=2.0*fx0 
      fz0=2.0*fz0
c-----------------------------------------------------------------------
c cw0 e cfz0 sono i coefficienti di resistenza e portanza
c sup e' l'area TOTALE della carena
c-----------------------------------------------------------------------
      cw0 = 2.0*fx0/sup
      cfz0 = 2.0*fz0/sup
c-----------------------------------------------------------------------
c     STOP se si vuole solo il doppio modello.
c-----------------------------------------------------------------------
      if (kbase.eq.1) goto 999
c-----------------------------------------------------------------------
c     CICLO  SU  PIU' FROUDE
c-----------------------------------------------------------------------

c      do 1000 ifr=1,nfr

c-----------------------------------------------------------------------
c    Riempie ... la matrice amat
c-----------------------------------------------------------------------

      call riama(x1,x2,x3,x4,xn,xh,y1,y2,y3,y4,yn,yh,
     &           z1,z2,z3,z4,zn,zh,cosxn,cosyn,coszn,
     &           xns,yns,zns,nun,
     &           area,fi0x,fi0y,fi0z,vx,vy,vz,fi0ll,vll,
     &           amat,sigma0,sigma1,
     &           diana,nekel,fr,itriqua,modo,
     &           ncar,nplib,nlast,ntra,
     &           ncard,ntotd)

c-----------------------------------------------------------------------
c     Condizione TRANSOM (solo se ntrapoc > 0 )
c-----------------------------------------------------------------------

      if(ntrapoc.gt.0) then
        nfile = 'transom.grd'
c        call checkfile(nfile,ierr)
c        if(ierr.eq.1) then
c          zztra = 0.
c          do i=1,ncar
c            if(xn(i).ge.0.4) zztra = min(zztra,zn(i))
c          end do
c          lontra = 2
c          open(44,file='stern.profile',status='unknown')
c          read(44,*) xa,ya,za
c          read(44,*) xb,yb,zb
c          close(44)
c          dzdx=(zb-za)/(xb-xa)
c          open(44,file=nfile,status='unknown')
c          write(44,*)  0.5,0.0
c          write(44,*) -1.0,zztra
c          write(44,*)  1.0,zztra
c          write(44,*)  dzdx
c          write(44,*)  dzdx
c          close(44)
c        end if
        call transom3(x1,x2,x3,x4,xh,xn,cosxn,
     $                y1,y2,y3,y4,yh,yn,cosyn,
     $                z1,z2,z3,z4,zh,zn,coszn,
     $                xns,yns,zns,area,vx,vy,vz,
     $                fi0x,fi0y,fi0z,amat,sigma1,
     $                nlast,ncar,nulib,nekel,
     $                fr,dryfr,lontra,diana,itriqua,cata,
     $                zmaxtra,nplib,ntrapoc,nfile)
      end if

c-----------------------------------------------------------------------
c     Soluzione sistema lineare
c-----------------------------------------------------------------------
      !write(*,*) 'Soluzione sistema'
      nsys = ncar+nplib
      sing = 1.E-5
      dx   = 1.E-5
      alpha = 1./float(nsys)
      do i=1,nsys
        r(i) = sigma1(i)
      end do
c -- Traspone
      do i=1,nsys-1
        do j=i+1,nsys
          rsave = amat(i,j)
          amat(i,j) = amat(j,i)
          amat(j,i) = rsave
        end do
      end do
      call simqit(amat,r,nsys,ntotd,alpha,dx,sigma1,ierr)
      !write(*,*) 'Sistema risolto'

c----------------------------------------------------------------------*
c     Calcolo fi1x, fi1y, fi1z
c----------------------------------------------------------------------*

      !write(*,*)'Calcolo di FI1x, FI1y e FI1z'

      do i=1,ntot
         if(i.le.ncar) then
           ising=i
         else if(i.gt.ncar.and.i.le.ncar+ntra) then
           ising=0
         else if(i.gt.ncar+ntra.and.i.le.nlast) then
           ising=i-ntra
         else if(i.gt.nlast) then
           ising=i
         end if
c         !write(*,fmt='($,A,I4.4)')'\b\b\b\b\b', i
         call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &                y1,y2,y3,y4,yh,yn,cosyn,
     &                z1,z2,z3,z4,zh,zn,coszn,diana,
     &                area,vx,vy,vz,ncar,
     &                xns(i),yns(i),zns(i),ising,1,
     &                ntot,itriqua,ncard,ntotd,nplib)
         call calve(i,ncar,vx,vy,vz,fi1x,fi1y,fi1z,sigma1,
     &              ncard,ntotd,nplib)

      end do

c      open(20,file='f1.dbg',status='unknown',recl=10000)
c      do i=1,ntot
c        write(20,*) fi1x(i),fi1y(i),fi1z(i),sigma1(i)
c      end do
c      close(20)

c----------------------------------------------------------------------*
c     Calcola l'elevazione
c----------------------------------------------------------------------*

      call elepre(nekel,ncar,nplib,
     &            eta,fr,xns,yns,
     &            fi0x,fi0y,fi0z,fi1x,fi1y,fi1z,
     &            pres0,pres1,pres12,presto,
     &            ncard,nplibd,ntotd)

c----------------------------------------------------------------------*
c     Hydrodynamics loads
c----------------------------------------------------------------------*

      cyd = 0.0
      if(ntrapoc.gt.0) then
        call idrosta(dryfr,lontra,fr,cyd)
      end if

      call hydroloads(rw,cw,cyd,cfy,cfz,cmy,snkg,atrim,
     &                cosxn,cosyn,coszn,
     &                area,fr,xn,yn,zn,presto,
     &                xG,xG0,zG0,dislo0,fextz,mexty,
     &                sup,dislo,swl,wli,actsink,acttrim,
     &                ncar,nplib,ntraint,ncard,ntotd)

      ztrasla = draft
      if(ztrasla.eq.0) then
        ztrala = 0.
        do i=1,ncar
          ztrasla = max(ztrasla,-z1(i))
          ztrasla = max(ztrasla,-z2(i))
          ztrasla = max(ztrasla,-z3(i))
          ztrasla = max(ztrasla,-z4(i))
        end do
      end if

      call ittcnl(fr,cw,snkg,atrim,ztrasla,
     &            Cf_ittc,Cf,Ct_ittc,Ct,Cyd,
     &            sup,shiplen,area,
     &            fi0x,fi0y,fi0z,
     &            fi1x,fi1y,fi1z,
     &            presto,
     &            xn,yn,zn,
     &            cosxn,cosyn,coszn,
     &            ncar,nplib,ncard,ntotd)

      call risulta(rw,cw,cyd,cw0,Cf_ittc,Cf,Ct_ittc,Ct,snkg,atrim,
     &             swl,wli,cmy,cfz,sup,fr,carena,eta,
     &             ncar,ncard,nplib,nun,nct,nplibd)

c----------------------------------------------------------------------*
c     Chiama le routine per tatra e per il plottaggio
c----------------------------------------------------------------------*

      call tracut(nlast,nlo,nlovallec,ntraest,ntraint,ntrapoc,
     &            ncar,nplib,fr,sup,eta,xns,yns)
c
c      call elunca(ncar,ntotd,nplibd,nplib,eta,xn,zn,fr(ifr),
c     &            ntraint,ntra,nlo)

      if((iaso.eq.1).or.(iaso.eq.4)) then
        call tecout(fr,fi0x,fi0y,fi0z,fi1x,fi1y,fi1z,
     &              eta,presto,xns,yns,zns,
     &              ncar,nulib,ncard,nplib,
     &              ntra,ntraint,ntrapoc,nlo,
     &              nlovallec,ntotd,nplibd)
      end if
c      frn = fr(ifr)
c      call flotta(frn,cw,ncar,fi0x,fi1x,fi0y,fi1y,
c     &            fi0z,fi1z,eta,presto,
c     &            nplib,ntrapoc,nlovallec,
c     &            ntra,ntraint,nlo,nulib,
c     &            ncard,nplibd,ntotd)

      call uniras(fr,xns,yns,zns,eta,presto,
     &            fi0x,fi0y,fi0z,fi1x,fi1y,fi1z,
     &            sup,itriqua,ncar,ncard,nplibd,ntotd,nplib)

c##################################################################
c 1000 continue
c##################################################################

 999  continue
c      tend = secondi()
c      tend = tend - tstart

c      write( * ,*)
c     &  '-----------------------------------------------------'
c      !write(*,*) 'Conto completato in',tend,'secondi.'
      !write(*,*)'----------========== C L O S E D ==========----------'

c      write(nun,*)
c     &  '-----------------------------------------------------'
c      write(nun,*) 'Conto completato in',tend,'secondi.'
      write(nun,*)
     &          '----------========== C L O S E D ==========----------'

      close(nun)
      close(nct)

!      close(stdout)

  222 format(100e16.8) 
c ----------------------------------------------------------------------
c      stop 'WARP ha finito'
      end
c-----------------------------------------------------------------------

      subroutine risulta(rw,cw,cyd,cw0,Cf_ittc,Cf,Ct_ittc,Ct,snkg,atrim,
     &                   swl,wli,cmy,cfz,sup,fr,carena,eta,
     &                   ncar,ncard,nplib,nun,nct,nplibd)

c-----------------------------------------------------------------------
      implicit real (a-h,o-z)
      real eta(nplibd)
      character*80 carena
c-----------------------------------------------------------------------

      write(nun,*)
     &  '----------========== O U T P U T ==========----------'
      write(nun,103)'Ship                         :  Carena  =',carena
      write(nun,102)'Waterline area               :  Swl     =',swl
      write(nun,102)'Waterline inertia moment     :  Wli     =',wli
      write(nun,102)'Wetted surface area          :  Sup     =',sup
      write(nun,102)'Froude                       :  Fr      =',fr
      write(nun,102)'Wave resistance (linear) WARP:  Cw      =',Cw
      write(nun,102)'Hydrost. resistance          :  Cyd     =',Cyd
      write(nun,102)'Frictional resistance    WARP:  Cf      =',Cf
      write(nun,102)'Frictional resistance    ITTC:  Cf_ittc =',Cf_ittc
      write(nun,102)'Total resistance         WARP:  Ct      =',Ct
      write(nun,102)'Total resistance         ITTC:  Ct_ittc =',Ct_ittc
      write(nun,102)'Vertical Force Coefficient   :  Cfz     =',Cfz
      write(nun,102)'Horizontal Moment Coefficient:  Cmy     =',Cmy
      write(nun,102)'Sinkage                      :  snkg    =',snkg
      write(nun,102)'Trim                         :  trim    =',atrim
c
      !write( * ,*)
      !&  '----------========== O U T P U T ==========----------'
      !write( * ,103)'Ship                         :  Carena  =',carena
      !write( * ,102)'Waterline area               :  Swl     =',swl
      !write( * ,102)'Waterline inertia moment     :  Wli     =',wli
      !write( * ,102)'Wetted surface area          :  Sup     =',sup
      !write( * ,102)'Froude                       :  Fr      =',fr
      !write( * ,102)'Wave resistance (linear) WARP:  Cw      =',Cw
      !write( * ,102)'Hydrost. resistance          :  Cyd     =',Cyd
      !write( * ,102)'Frictional resistance    WARP:  Cf      =',Cf
      !write( * ,102)'Frictional resistance    ITTC:  Cf_ittc =',Cf_ittc
      !write( * ,102)'Total resistance         WARP:  Ct      =',Ct
      !write( * ,102)'Total resistance         ITTC:  Ct_ittc =',Ct_ittc
      !write( * ,102)'Vertical Force Coefficient   :  Cfz     =',Cfz
      !write( * ,102)'Horizontal Moment Coefficient:  Cmy     =',Cmy
      !write( * ,102)'Sinkage                      :  snkg    =',snkg
      !write( * ,102)'Trim                         :  trim    =',atrim
c
c - as
c      etamin = 1.e6
c      etamax =-1.e6
c      do i = 1,nplib
c          if(eta(i).le.etamin) etamin=eta(i)
c          if(eta(i).ge.etamax) etamax=eta(i)
c      end do
       etamin = minval(eta(1:nplib))
       etamax = maxval(eta(1:nplib))
c - as

      write(nun,*)'ETA min. =',etamin,' ETA max. =',etamax 
      write(nct,101) fr,Cw,Cyd,Cf,Cf_ittc,Ct,Ct_ittc,snkg,atrim,sup

c-----------------------------------------------------------------------
 100  format(1x,g13.6,1x,g13.6,1x,g13.6,1x,g13.6,1x,g13.6,
     &       1x,g13.6,1x,g13.6,1x,g13.6,1x,g13.6)
 101  format(1x,e13.6,1x,e13.6,1x,e13.6,1x,e13.6,1x,e13.6,
     &       1x,e13.6,1x,e13.6,1x,e13.6,1x,e13.6,1x,e13.6)
 102  format(1x,a41,e12.5)
 103  format(1x,a41,1x,a11)
c-----------------------------------------------------------------------
      return                                                    
      end
c-----------------------------------------------------------------------
      subroutine uscita(nfr,diana,dryfr,nekel,ncar,nplib,ntra,
     &                  nlo,ntraint,ntraest,ntrapoc,nlovallec,
     &                  modo,lontra,cata,nun)
c-----------------------------------------------------------------------
c
      implicit real (a-h,o-z)
      character*50  cha,cha1  
      character*1   cata  
c-----------------------------------------------------------------------
c
      !write( * ,1101) nfr
      !write( * ,1102) ncar,nplib
      !write( * ,1103) ntra,nlo
      !write( * ,1110) ntraint,ntraest
      !write( * ,1004) ntrapoc,nlovallec
      !write( * ,1109) lontra
      !write( * ,1106) diana,dryfr
      !write( * ,1107) nekel
      !write( * ,1108) modo
c
      write(nun,1101) nfr
      write(nun,1102) ncar,nplib
      write(nun,1103) ntra,nlo
      write(nun,1110) ntraint,ntraest
      write(nun,1004) ntrapoc,nlovallec
      write(nun,1109) lontra
      write(nun,1106) diana,dryfr
      write(nun,1107) nekel
      write(nun,1108) modo

      if(ntrapoc.eq.0) then
        cha=' con poppa chiusa '
      else
        cha=' con poppa transom'
      end if

      if(ntraint.eq.0) then
        cha1=' Monocarena'//cha
      else
        cha1='      SWATH'//cha
      end if
      !write(*,110) cha1

c-----------------------------------------------------------------------
 110  format(a)
 1004 format(4x,' ntrapoc =',i4,4x,' nlovallec =',i3)
 1101 format(/,20x,'nfr =',i3,/)
 1102 format(4x,'    ncar =',i4,4x,'    nplib  =',i6)
 1103 format(4x,'    ntra =',i4,4x,'     nlo   =',i3)
 1106 format(4x,'   diana =',f7.3,4x,'  dryfr =',f5.3)
 1107 format(4x,'   nekel =',i3)
 1108 format(4x,'    modo =',i3)
 1109 format(4x,'  lontra =',i3)
 1110 format(4x,' ntraint =',i4,4x,'   ntraest =',i3)
c-----------------------------------------------------------------------
      return                                                    
      end
c-----------------------------------------------------------------------
      subroutine swarp
c-----------------------------------------------------------------------

      !write(*,*)' '
      !write(*,*)'                    ####'
      !write(*,*)'                  ###'
      !write(*,*)'                ##           W A R P         V 1.5'
      !write(*,*)'               _|__ INSEAN - WAve Resistance Program'
      !write(*,*)'   ____________|*_|__.'
      !write(*,*)'~~~\\._______________/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      !write(*,*)' '

c-----------------------------------------------------------------------
      return 
      end
c-----------------------------------------------------------------------
      subroutine swarp_sa
c-----------------------------------------------------------------------

      character*80 s1

c      call system ('hostname -A > swarp.aux')
c      open (999,file='swarp.aux',status='old')
c      read (999,*) s1
c      close(999)

      write(*,*)' '
      write(*,*)'                    ####'
      write(*,*)'                  ###'
      write(*,*)'                ##                    W A R P'
      write(*,*)
     &'               _|__      CNR-INSEAN, WAve Resistance Program'
      write(*,*)'   ____________|*_|__.'
      write(*,*)
     &'~~~\\._______________/~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*)' '
      write(*,*)'                                Stand-Alone Version'
      write(*,*)'                                  V1.1 (28Mar2014)'
      write(*,*)' '
      write(*,*)'                                Embeds:'
      write(*,*)'                                  Panca'
      write(*,*)'                                  Free-surface V2.4'
      write(*,*)'                                  Warp V1.5 OMP'
      write(*,*)'                                  Pressione'
      write(*,*)'                                  Free-warp'
      write(*,*)' '
c      write(*,*)
c     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c      write(*,*)'Running on'
c      write(*,*)' ',trim(s1)
      write(*,*)
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*) ' '

c-----------------------------------------------------------------------
      return 
      end
c-----------------------------------------------------------------------
      subroutine elunca(ncar,ntotd,nplibd,nplib,eta,xn,zn,fr,
     &                  ntraint,ntra,nlo)
c-----------------------------------------------------------------------

      implicit real (a-h,o-z)

      real  xn(ntotd),zn(ntotd),eta(nplibd)
      character*20 nfile
      character*5  cha  

c-----------------------------------------------------------------------

      ifr = fr*1000.0 + 0.1
      write(cha,99) ifr
      nfile= 'e'//cha
      if (ntraint.gt.0) nfile= 'ee'//cha
c
      open(unit=30,file=nfile,status='unknown')
      rewind(30)
      do i=1,nlo
      icont=ncar+1+ntraint+ntra*(i-1)
         write(30,*)xn(icont),eta(icont-ncar)
      end do
      close(30)
c
      if (ntraint.gt.0) then
c
        nfile= 'ei'//cha
        open(unit=30,file=nfile,status='unknown')
        rewind(30)
        do i=1,nlo
          icont=ncar+ntraint+ntra*(i-1)
          write(30,*)xn(icont),eta(icont-ncar)
        end do
        close(30)
      end if

c-----------------------------------------------------------------------
   99 format(i5.5)
c-----------------------------------------------------------------------
      return 
      end
c-----------------------------------------------------------------------
      subroutine uniras(fr,xns,yns,zns,eta,presto,
     &                  fi0x,fi0y,fi0z,fi1x,fi1y,fi1z,
     &                  sup,itriqua,ncar,ncard,nplibd,ntotd,nplib)
c-----------------------------------------------------------------------

      implicit real (a-h,o-z)

      integer itriqua(ntotd)
      real xns(ntotd),yns(ntotd),zns(ntotd)
      real presto(ncard)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd),
     &     fi1x(ntotd),fi1y(ntotd),fi1z(ntotd)
      real eta(nplibd)
      character*10 fina

c-----------------------------------------------------------------------

      fina(1:6) = 'psuni.'
      nufru = nint(fr*1000)
      write(fina(7:10),'(i4.4)') nufru
      open(unit=99,file=fina,status='unknown',
     &             form='formatted',recl=64000)
      write(99,*) ncar
      do i=1,ncar
        vx = 1.+fi0x(i)+fi1x(i)
        vy =    fi0y(i)+fi1y(i)
        vz =    fi0z(i)+fi1z(i)
        write (99,222) xns(i),yns(i),zns(i),vx,vy,vz,presto(i)
      end do
      close(99)

  222 format(100e16.8) 
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c Hydrodynamics loads
c-----------------------------------------------------------------------
c     Calcolo della resistenza d'onda (Rw) e del coeff. (Cw) con
c             pressioni statiche e dinamiche 
c     Calcolo della forza verticale (Fz) e del momento intorno ad y (Ym)
c     Calcolo di sinkage & trim. 
c      SNKG < 0 ===> Immersioni
c      TRIM > 0 ===> Appoppamenti
c                                                                       
c                                       ^ z                             
c                                       |
c                                       |
c                                 trim<---
c        x <--------------------------- |----------------------
c                  \                    |                     / 
c                   ----\               V snkg               /
c                        \                                  /
c                         \--------------------------------/
c
c     TUTTE LE GRANDEZZE SONO RELATIVE ALLA CARENA COMPLETA (y>0 & y<0)
c-----------------------------------------------------------------------
      subroutine hydroloads(rw,cw,cyd,cfy,cfz,cmy,snkg,atrim,
     &                      cosxn,cosyn,coszn,
     &                      area,fr,xn,yn,zn,presto,
     &                      xG,xG0,zG0,dislo0,fextz,mexty,
     &                      sup,dislo,swl,wli,actsink,acttrim,
     &                      ncar,nplib,ntraint,ncard,ntotd)
c-----------------------------------------------------------------------

      implicit real (a-h,o-z)

      real xn(ntotd),yn(ntotd),zn(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real area(ntotd)
      real presto(ncard)
      real fextz,mexty

c-----------------------------------------------------------------------

      rw     = 0.0
      fz     = 0.0
      fy     = 0.0
      ym     = 0.0
      fxstat = 0.0
      fystat = 0.0
      fzstat = 0.0
      xmh    = 0.0
      ymh    = 0.0
      zmh    = 0.0
      frq    = fr*fr

      xpolo = 0.0
      ypolo = 0.0
      zpolo = 0.0

c      open(888,file='hydro.dbg',status='unknown',recl=10000)

c
c Calcola solo su meta' carena (y>0)
c
      do i=1,ncar
        rw     =     rw - presto(i)*cosxn(i)*area(i)
        fy     =     fy - presto(i)*cosyn(i)*area(i)
        fz     =     fz - presto(i)*coszn(i)*area(i)
        ym     =     ym - presto(i)*(zn(i)*cosxn(i)
     &                              -xn(i)*coszn(i))*area(i)
       
c        write(888,*) fz,xn(i),zn(i),presto(i),
c     &               cosxn(i),cosyn(i),coszn(i),area(i)

        dfxstat = zn(i) * cosxn(i) * area(i) / frq
        dfystat = zn(i) * cosyn(i) * area(i) / frq
        dfzstat = zn(i) * coszn(i) * area(i) / frq
        fxstat  = fxstat + dfxstat
        fystat  = fystat + dfystat
        fzstat  = fzstat + dfzstat
        xb = xn(i)-xpolo
        yb = yn(i)-ypolo
        zb = zn(i)-zpolo
c--md
c        xmh     = xmh - (yb*dfzstat-zb*dfystat)
c        ymh     = ymh - (zb*dfxstat-xb*dfzstat)
c        zmh     = zmh - (xb*dfystat-yb*dfxstat)
        xmh     = xmh + (yb*dfzstat-zb*dfystat)
        ymh     = ymh + (zb*dfxstat-xb*dfzstat)
        zmh     = zmh + (xb*dfystat-yb*dfxstat)
c--md

      end do

c      close(888)

c -- Ora su TUTTA la carena

      rw      = 2.*(rw + fxstat)
      fz      = 2.*fz
      ym      = 2.*ym
      fzstat  = 2.*fzstat
      ymh     = 2.*ymh

c      write(*,*) 'fz,fextz,ym,mexty'
c      write(*,*) fz,fextz,ym,mexty

c --

      fz = fz + fextz
      ym = ym + mexty

c --
      cw  = 2.0*rw/sup 
      cfy = 2.0*fy/sup 
      cfz = 2.0*fz/sup 
      cmy = 2.0*ym/sup 

cc -- Componente peso
c
c      zGr = zG + actsink
c      xGr = xG + zGr*tan(acttrim)
c
cc -- Il dislocamento e' gia' quello totale
c      
c      ym  = ym + (xGr-xpolo)*dislo*0.5
c
cc Sinkage & trim: non vanno inserite le azioni di massa
c
c      snkg  =  fz*frq/swl
c      atrim =  ym*frq/wli


c--md

      xGr =  xG0 * cos(acttrim) + zG0 * sin(acttrim)

      rym = ym                             !momento esterno              
     &         + (xG -xpolo)*(-dislo /frq) !momento idrostatico + idrodinamico
     &         + (xGr-xpolo)*(+dislo0/frq) !momento forza peso
      rfz = fz                             !forza esterna
     &         + dislo /frq                !forza idrostatica + idrodinamica
     &         - dislo0/frq                !forza peso

c      rym = ym + ymh
c     &         + (xGr-xpolo)*(+dislo0/frq)
c      rfz = fz + fzstat 
c     &         - dislo0/frq 

      snkg  =   rfz*frq/swl + actsink
      atrim =   rym*frq/wli + acttrim

c--md


c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
c Calcola elevazione e pressione
c-----------------------------------------------------------------------

      subroutine elepre(nekel,ncar,nplib,
     &                  eta,fr,xns,yns,
     &                  fi0x,fi0y,fi0z,fi1x,fi1y,fi1z,
     &                  pres0,pres1,pres12,presto,
     &                  ncard,nplibd,ntotd)

c-----------------------------------------------------------------------

      implicit  real (a-h,o-z)

      real xns(ntotd),yns(ntotd)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd),
     &     fi1x(ntotd),fi1y(ntotd),fi1z(ntotd)
      real pres0(ncard),pres1(ncard),pres12(ncard),presto(ncard)
      real eta(nplibd)

c-----------------------------------------------------------------------
      !call checknana(xns,ntotd,ierr)
      !print *, ierr
      !call checknana(yns,ntotd,ierr)
      !print *, ierr
      !call checknana(fi0x,ntotd,ierr)
      !print *, ierr
      !call checknana(fi0y,ntotd,ierr)
      !print *, ierr
      !call checknana(fi0z,ntotd,ierr)
      !print *, ierr
      !call checknana(fi1x,ntotd,ierr)
      !print *, ierr
      !call checknana(fi1y,ntotd,ierr)
      !print *, ierr
      !call checknana(fi1z,ntotd,ierr)
      !print *, ierr
      !call checknana(pres0,ncard,ierr)
      !print *, ierr
      !call checknana(pres1,ncard,ierr)
      !print *, ierr
      !call checknana(pres12,ncard,ierr)
      !print *, ierr
      !call checknana(presto,ncard,ierr)
      !print *, ierr
      !call checknana(eta,nplib,ierr)
      !print *, ierr
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c     Calcola l'elevazione
c-----------------------------------------------------------------------

      frq = fr*fr
      open(22,file='eta-raw.plt',form='formatted',status='unknown')
      write(22,*) 'VARIABLES = "X" "Y" "h"'
      write(22,*) 'ZONE'
      do i=1,ncar
        write(22,*) xns(i),yns(i),0.
      end do
      if(nekel.eq.1) then
        do i=1,nplib
          l = i + ncar 
          eta(i) = -frq*(fi0x(l)+fi1x(l)) 
          write(22,*) xns(l),yns(l),eta(i)
        end do
      else
        do i=1,nplib
          l = i + ncar
          eta(i)=-0.5*frq*(fi0x(l)**2+fi0y(l)**2+fi0z(l)**2+
     &           2.0*(fi0x(l)+fi1x(l)+fi0x(l)*fi1x(l)+fi0y(l)*fi1y(l)+
     &           fi0z(l)*fi1z(l)))
          write(22,*) xns(l),yns(l),eta(i)
        end do
      end if
      close(22)

c----------------------------------------------------------------------*
c     Calcola la pressione: 
c     P = p0 + p1 + p_statica + p_nonlin
c----------------------------------------------------------------------*

c      open(887,file='pres.dbg',status='unknown',recl=10000)

      if(nekel.eq.1) then
        do i=1,ncar
          pres1(i) = -fi1x(i)
          presto(i)=pres0(i)+pres1(i)     ! pressione totale
          if(abs(presto(i)).gt.0.5)
     &           presto(i)=0.5*presto(i)/abs(presto(i))
        end do 
      else if(nekel.eq.0) then
        do i=1,ncar
          pres1(i) = -(fi1x(i)+fi0x(i)*fi1x(i)+fi0y(i)*fi1y(i)+
     &                 fi0z(i)*fi1z(i))  
          pres12(i)= -0.5*(fi1x(i)**2+fi1y(i)**2+fi1z(i)**2)
          presto(i)= pres0(i)+pres1(i)+pres12(i)

c          write(887,*) pres0(i),pres1(i),pres12(i),presto(i)

          if(abs(presto(i)).gt.0.5)
     &           presto(i)=0.5*presto(i)/abs(presto(i))
        end do 
      end if

c      close(887)
c-----------------------------------------------------------------------
      return
      end

c-----------------------------------------------------------------------
c  Check NaN
c-----------------------------------------------------------------------

      subroutine checknana(rarray,ndim,ierr)

      real rarray(ndim)

      ierr=0
      icount=0

      do i=1,ndim
        if ((rarray(i).eq.rarray(i)+1.)  .or.
     &      (rarray(i).ne.rarray(i)   )) then
          rarray(i)=0.
          ierr=1
          icount=i+1
        endif
      enddo

      print *, "found NaN in elepre: icount =",icount

      return
      end

c-----------------------------------------------------------------------
c  Termine idrostatico
c-----------------------------------------------------------------------

      subroutine idrosta(dryfr,long,fr,cyd)

c-----------------------------------------------------------------------

      implicit real (a-h,o-z)
      dimension yt(100),zt(100)

c-----------------------------------------------------------------------

      !write(*,*) 'Calcolo del termine idrostatico'

      cyd  = 0.0
      !write(*,*) 'DryFr =',dryfr
      !write(*,*) 'LONG  =',long
      !write(*,*) 'Fr    =',fr
      if(dryfr.eq.0) return
      if(long.le.2) return

      open(unit=20,file='transom.grd',status='old')
      rewind(20)
      read(20,*) 
      do i=1,long
        read(20,*) yt(i),zt(i)
      end do
      close(20)

      cyd1 = 0.0
      do i=1,long
        dz = (zt(i+1)+zt(i))
        dy = (yt(i+1)+yt(i))
        if(dy.ne.0.and.dz.ne.0) then
          cyd1=cyd1+((dz*0.5)**2*
     &             abs(dy)*2.)/fr**2
        end if
      end do

      if(fr.lt.dryfr) then
        do i=1,long
          zt(i) = zt(i)*(fr/dryfr)**2
        end do
      end if

      cyd2 = 0.0
      do i=1,long
        dz = (zt(i+1)+zt(i))
        dy = (yt(i+1)+yt(i))
        if(dy.ne.0.and.dz.ne.0) then
          cyd2=cyd2+((dz*0.5)**2*
     &             abs(dy)*2.)/fr**2
        end if
      end do

      call checknan(cyd1)
      call checknan(cyd2)

      !write(*,*) 'Aggiungo ',cyd1
      !write(*,*) 'Tolgo    ',cyd2

      cyd = 0.5*(cyd1-cyd2)

c ----------------------------------------------------------------------
      return 
      end
c ----------------------------------------------------------------------
c
c      real function secondi()
c
c ----------------------------------------------------------------------

c      real*8 omp_get_wtime
c      external omp_get_wtime

c ----------------------------------------------------------------------


c  For IBM machine
c      second = mclock()/1000.
c      t = etime_(tarray)
c      secondi = tarray(1)

c  For Macintosh machine
c       n = LONG(362)
c       secondi = FLOAT(n)/60.0

c  For Unix machine
c      t = etime(tarray)
c      secondi = tarray(1)
c
c  For OpenMP
c      secondi = sngl(omp_get_wtime())
c      secondi = 1.

c ----------------------------------------------------------------------
c      return
c      end
c ----------------------------------------------------------------------

      subroutine checkfile(filen,ierr)

c ----------------------------------------------------------------------
c
c  Calcola il numero di righe del file 'filen'
c  e lo restituisce in 'nl'
c
c ----------------------------------------------------------------------

      character*80 filen

c ----------------------------------------------------------------------

      ierr = 1
      open(85,file=filen,status='old',err=1000)
      ierr = 0

c ----------------------------------------------------------------------
 1000 close(85)
      end
c ----------------------------------------------------------------------

c      subroutine checknan(x)
c
c ----------------------------------------------------------------------
c
c      if(x.ne.x) x = 0.
c      if(x.eq.x+1) x = 0.
c
c ----------------------------------------------------------------------
c      return
c      end
c ----------------------------------------------------------------------
