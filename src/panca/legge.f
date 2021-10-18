      subroutine LEGGE(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xn,yn,zn,
     &                 cosxn,cosyn,coszn,lontra,carena,npan,itra,
     &                 geofile,swath,yswath,shiplen,xtot,eps,ktra,sym,
     &                 n_est,n_int,numfile,dist_pp,ltot,xtrasla,ztrasla)
c=======================================================================
c  1) Lettura dei file di griglia (uno per ogni dominio)
c  2) Adimensionalizzazione della geometria
c  3) Costruzione dei vettori dei punti di confine e di controllo
c  4) Calcolo delle normali con i punti non complanari
c  5) Scrittura file CARENA_0.DAT
c
c     ATTENZIONE !
c
c       * NPAN e' la somma dei pannelli dei domini
c       * LONTRA e' il long del dominio transom che DEVE ESSERE UNICO !
c
c=======================================================================
c
      parameter(ksez=200,klin=200)
      dimension xappo(ksez,klin),yappo(ksez,klin),zappo(ksez,klin),
     &          appox(ksez,klin),appoy(ksez,klin),appoz(ksez,klin)
      dimension x1(*),x2(*),x3(*),x4(*),
     &          y1(*),y2(*),y3(*),y4(*),
     &          z1(*),z2(*),z3(*),z4(*),
     &          xn(*),yn(*),zn(*),cosxn(*),cosyn(*),coszn(*)
      dimension itra(*)
      dimension x_0(10000),y_0(10000)
      dimension x_1(10000),y_1(10000)
c
      real         dz(20)
      real         ltot
      integer      nlift(20)
      integer      ngrid
c
      integer      iuno,iunit,ndom,lat_est,lontra,ktra
c
      character*80 fina
      character*80 carena
      character*80, dimension(:), allocatable ::filenamei,filenameo
      character*20 namepres
      character*3  cha
      character*1  flare,geofile
      character*1  fsdom(20),trdom(20),jumbo(20),rfloi(20),rfloj(20)
      character*1  ex_ij(20),cutsec(20),lift(20)
      character*1  xmirror,ymirror,zmirror
      character*1  swath,sym

      integer      inode,jnode
      integer      iaso,isolv,igeomod,igrid,numfilest,numfilint
c
      common  /appo/ lat,long,snkg,trim,xmaxfs,ymaxfs,kpan

      namelist/MAIN_PARAMETERS/   iaso,isolv,igeomod,ngrid,igrid,
     &                            numfilest,numfilint
      namelist/GRID_PARAMETERS/   filenamei,filenameo,inode,jnode
      namelist/PANCA_PARAMETERS/  carena,shiplen,flare,nflare,swath,
     &                            yswath,sym,dist_pp,stre_le,xtrasla,
     &                            ztrasla,xmirror,ymirror,zmirror,
     &                            areamin,geofile,fsdom,trdom,
     &                            jumbo,cutsec,lift,nlift,dz,ex_ij,
     &                            rfloi,rfloj

c------------------------------------------------------------------------

      iunit  = 20

      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=MAIN_PARAMETERS)
      close(iunit)

      allocate(filenamei(numfilest+numfilint))
      allocate(filenameo(numfilest+numfilint))

      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=GRID_PARAMETERS)
      close(iunit)

      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=PANCA_PARAMETERS)
      close(iunit)

      open(iunit,file="pancanml.aux",status="unknown",recl=10000)
        write(iunit,*) ztrasla
        write(iunit,*) xmirror
        write(iunit,*) dist_pp
      close(iunit)

      open(iunit,file='sinktrim.aux',status='unknown',recl=10000)
        read(iunit,*) snkg
        read(iunit,*) trim
      close(iunit)

      iuno = 1
c
      lat_0   = 0          ! Contatore punti di intersezione carena-f.s.
      lat_est = 0          ! Contatore punti di intersezione carena-f.s.
      lat_int = 0          ! Contatore punti di intersezione carena-f.s.
      lat_0   = 0          ! Contatore punti di intersezione carena-f.s.
      kpan    = 0          ! Contatore dei pannelli.
      ktra    = 0          ! diventa 1 se e' una carena transom
      xmaxfs  = -1.e10
      xmax    = -1.e10
      xmin    =  1.e10
      ymax    = -1.e10
      ymin    =  1.e10
      zmax    = -1.e10
      zmin    =  1.e10

      numfile = numfilest + numfilint
c
      pig  = acos(-1.0)
      xtot =  1.0 + (stre_le  / dist_pp)  !!! Se non si allunga, XTOT=1
c
c  !write(*,*)'XTOT =',xtot,':  DIST_PP=',dist_pp,': STRE_LE=',stre_le
c
      if(swath.eq.'y') then
         !write(*,*) 'Carena CAT/SWATH: interasse scafi =',(yswath*2.)/xtot
      else
         yswath=0.
      end if
      yswath  = yswath  / dist_pp
      stre_le = stre_le / dist_pp
c-----------------------------------------------------------------------
c Verifica compatibilita' parametri iniziali
c-----------------------------------------------------------------------
         if(abs(snkg).le.eps.and.abs(trim).le.eps)then
              snkg=5.*eps
              trim=5.*eps
         end if
c-----------------------------------------------------------------------
c     CICLO SUI DOMINI ESTERNI
c-----------------------------------------------------------------------
      do 333 ndom=1,numfilest
c 100    format(a)
        !write(*,*) '=---------------------------------------------='
        !write(*,*) '              Dominio =',ndom
        !write(*,*) '=---------------------------------------------='
c-----------------------------------------------------------------------
c  Legge un file di griglia
c  Normalizza le coordinate e yswath  ...........(dist_pp)
c  Centra il sistema di riferimento .....(xtrasla,ztrasla)
c  Trasla i domini se la carena viene allungata  (stre_le)
c    (se la carena viene allungata, sposta tutti i domini)
c-----------------------------------------------------------------------
        CALL INPGRI(filenameo(ndom),ksez,nodi_i,nodi_j,xmirror,ymirror,
     &              zmirror,appox,appoy,appoz,ztrasla,xtrasla,dist_pp,
     &              jumbo(ndom),stre_le)
c-----------------------------------------------------------------------
c--- Scambia i con j
c-----------------------------------------------------------------------
        if(ex_ij(ndom).eq.'y') then
           CALL SCAMBIO(ksez,nodi_i,nodi_j,appox,appoy,appoz,
     &                  xappo,yappo,zappo)
           CALL AZZERA(ksez,klin,appox,appoy,appoz,xappo,yappo,zappo)
        end if
c-----------------------------------------------------------------------
c--- Inverte i
c-----------------------------------------------------------------------
        if(rfloi(ndom).eq.'y') then
           CALL INV_I(ksez,nodi_i,nodi_j,appox,appoy,appoz,
     &                xappo,yappo,zappo)
           CALL AZZERA(ksez,klin,appox,appoy,appoz,xappo,yappo,zappo)
        end if
c-----------------------------------------------------------------------
c--- Inverte j
c-----------------------------------------------------------------------
        if(rfloj(ndom).eq.'y') then
           CALL INV_J(ksez,nodi_i,nodi_j,appox,appoy,appoz,
     &                xappo,yappo,zappo)
           CALL AZZERA(ksez,klin,appox,appoy,appoz,xappo,yappo,zappo)
        end if
c-----------------------------------------------------------------------
c--- Allunga la carena
c-----------------------------------------------------------------------
        if(jumbo(ndom).ne.'n'.and.cutsec(ndom).eq.'y') then
           CALL JUMBOSHIP(appox,appoy,appoz,xappo,yappo,zappo,
     &                    stre_le,jumbo(ndom),nodi_i,nodi_j,ksez)
        end if
c-----------------------------------------------------------------------
c--- Alza la carena
c-----------------------------------------------------------------------
        if(lift(ndom).ne.'n') then
           CALL LIFTING(appox,appoy,appoz,xappo,yappo,zappo,
     &                  nodi_i,nodi_j,dz(ndom),nlift(ndom),ksez)
        end if
c-----------------------------------------------------------------------
c--- Pannelli del dominio
c-----------------------------------------------------------------------
      lat  = nodi_i - 1 
      long = nodi_j - 1 
c-----------------------------------------------------------------------
c    Dimensioni del dominio
c-----------------------------------------------------------------------
      xmax  = -1.e10
      xmin  =  1.e10
      ymax  = -1.e10
      ymin  =  1.e10
      zmax  = -1.e10
      zmin  =  1.e10
      do j=1,nodi_j
      do i=1,nodi_i
        if(xmax.le.appox(i,j)) xmax=appox(i,j)
        if(xmin.ge.appox(i,j)) xmin=appox(i,j)
c
        if(ymax.le.appoy(i,j)) ymax=appoy(i,j)
        if(ymin.ge.appoy(i,j)) ymin=appoy(i,j)
c
        if(zmax.le.appoz(i,j)) zmax=appoz(i,j)
        if(zmin.ge.appoz(i,j)) zmin=appoz(i,j)
      end do
      end do
c
c-----------------------------------------------------------------------
c    Assetto della carena
c    Le coordinate vengono ora scalate sulla nuova lunghezza XTOT
c-----------------------------------------------------------------------
      call ASSETTO(appox,appoy,appoz,fsdom(ndom),trdom(ndom),
     &             x_0,y_0,lat_est,lontra,
     &             itra,ktra,x1,x2,x3,x4,xn,y1,y2,y3,y4,yn,
     &             z1,z2,z3,z4,zn,areamin,xtot,ndom,shiplen,iclose)
c-----------------------------------------------------------------------
c--- File per plot3d
c-----------------------------------------------------------------------

      if(iclose.eq.1) y_0(lat_est) = y_0(1)

      call scrivigriglia(appox,appoy,appoz,nodi_i,nodi_j,
     &                   ndom,xtot,yswath)

c-----
      if(swath.eq.'y') then
c-----
        if(sym.eq.'y')
     &    call scrivigrigliasimm(appox,appoy,appoz,nodi_i,nodi_j,
     &                           ndom+numfile,xtot,yswath)
c-----
      end if
c-----
 333  continue

      n_est=kpan
      call reinterpola(x_0,y_0,lat_est)
      lat_0=lat_est

c-----------------------------------------------------------------------
c     CICLO SUI DOMINI INTERNI
c-----------------------------------------------------------------------
c
      do 444 ndom=numfilest+1,numfile
        !write(*,*) '=---------------------------------------------='
        !write(*,*) '              Dominio =',ndom
        !write(*,*) '=---------------------------------------------='
c-----------------------------------------------------------------------
c  Legge un file di griglia
c  Normalizza le coordinate e yswath  ...........(dist_pp)
c  Centra il sistema di riferimento .....(xtrasla,ztrasla)
c  Trasla i domini se la carena viene allungata  (stre_le)
c    (se la carena viene allungata, sposta tutti i domini)
c-----------------------------------------------------------------------
        CALL INPGRI(filenameo(ndom),ksez,nodi_i,nodi_j,xmirror,ymirror,
     &              zmirror,appox,appoy,appoz,ztrasla,xtrasla,dist_pp,
     &              jumbo(ndom),stre_le)
c-----------------------------------------------------------------------
c--- Scambia i con j
c-----------------------------------------------------------------------
        if(ex_ij(ndom).eq.'y') then
           CALL SCAMBIO(ksez,nodi_i,nodi_j,appox,appoy,appoz,
     &                  xappo,yappo,zappo)
           CALL AZZERA(ksez,klin,appox,appoy,appoz,xappo,yappo,zappo)
        end if
c-----------------------------------------------------------------------
c--- Inverte i
c-----------------------------------------------------------------------
        if(rfloi(ndom).eq.'y') then
           CALL INV_I(ksez,nodi_i,nodi_j,appox,appoy,appoz,
     &                xappo,yappo,zappo)
           CALL AZZERA(ksez,klin,appox,appoy,appoz,xappo,yappo,zappo)
        end if
c-----------------------------------------------------------------------
c--- Inverte j
c-----------------------------------------------------------------------
        if(rfloj(ndom).eq.'y') then
           CALL INV_J(ksez,nodi_i,nodi_j,appox,appoy,appoz,
     &                xappo,yappo,zappo)
           CALL AZZERA(ksez,klin,appox,appoy,appoz,xappo,yappo,zappo)
        end if
c-----------------------------------------------------------------------
c--- Allunga la carena
c-----------------------------------------------------------------------
        if(jumbo(ndom).ne.'n'.and.cutsec(ndom).eq.'y') then
           CALL JUMBOSHIP(appox,appoy,appoz,xappo,yappo,zappo,
     &                    stre_le,jumbo(ndom),nodi_i,nodi_j,ksez)
        end if
c-----------------------------------------------------------------------
c--- Alza la carena
c-----------------------------------------------------------------------
        if(lift(ndom).ne.'n') then
           CALL LIFTING(appox,appoy,appoz,xappo,yappo,zappo,
     &                  nodi_i,nodi_j,dz(ndom),nlift(ndom),ksez)
        end if
c-----------------------------------------------------------------------
c--- Pannelli del dominio
c-----------------------------------------------------------------------
      lat  = nodi_i - 1 
      long = nodi_j - 1 
c-----------------------------------------------------------------------
c    Dimensioni del dominio
c-----------------------------------------------------------------------
      xmax  = -1.e10
      xmin  =  1.e10
      ymax  = -1.e10
      ymin  =  1.e10
      zmax  = -1.e10
      zmin  =  1.e10
      do j=1,nodi_j
      do i=1,nodi_i
        if(xmax.le.appox(i,j)) xmax=appox(i,j)
        if(xmin.ge.appox(i,j)) xmin=appox(i,j)
c
        if(ymax.le.appoy(i,j)) ymax=appoy(i,j)
        if(ymin.ge.appoy(i,j)) ymin=appoy(i,j)
c
        if(zmax.le.appoz(i,j)) zmax=appoz(i,j)
        if(zmin.ge.appoz(i,j)) zmin=appoz(i,j)
      end do
      end do
c
      !write(*,*) 'Xmax =',xmax,': Xmin =',xmin
      !write(*,*) 'Ymax =',ymax,': Ymin =',ymin
      !write(*,*) 'Zmax =',zmax,': Zmin =',zmin
c
c-----------------------------------------------------------------------
c    Assetto della carena
c    Le coordinate vengono ora scalate sulla nuova lunghezza XTOT
c-----------------------------------------------------------------------
      CALL ASSETTO(appox,appoy,appoz,fsdom(ndom),trdom(ndom),x_1,y_1,
     &             lat_int,lontra,itra,ktra,x1,x2,x3,x4,xn,y1,y2,y3,y4,
     &             yn,z1,z2,z3,z4,zn,areamin,xtot,ndom,shiplen,iclose)
c-----------------------------------------------------------------------
c--- File per plot3d
c-----------------------------------------------------------------------

      if(iclose.eq.1) y_1(lat_int) = y_1(1)

      call scrivigriglia(appox,appoy,appoz,nodi_i,nodi_j,
     &                   ndom,xtot,yswath)

c-----
      if(swath.eq.'y') then
c-----
        if(sym.eq.'y')
     &    call scrivigrigliasimm(appox,appoy,appoz,nodi_i,nodi_j,
     &                           ndom+numfile,xtot,yswath)
c-----
      end if
c-----
 444  continue
      close(21)
c
      if(numfilint.gt.0) then
        call reinterpola(x_1,y_1,lat_int)
        lat_0=lat_int+lat_est
      end if
      npan = kpan                               ! Totale pannelli carena
c
c-----------------------------------------------------------------------
c     Calcolo delle normali con i punti non complanari
c-----------------------------------------------------------------------
       do i=1,n_est
 445    if(x1(i).eq.x3(i).and.y1(i).eq.y3(i).and.z1(i).eq.z3(i)) then
          do ii=i,n_est-1
            x1(ii)=x1(ii+1)
            x2(ii)=x2(ii+1)
            x3(ii)=x3(ii+1)
            x4(ii)=x4(ii+1)
            y1(ii)=y1(ii+1)
            y2(ii)=y2(ii+1)
            y3(ii)=y3(ii+1)
            y4(ii)=y4(ii+1)
            z1(ii)=z1(ii+1)
            z2(ii)=z2(ii+1)
            z3(ii)=z3(ii+1)
            z4(ii)=z4(ii+1)
          end do
          n_est = n_est - 1
          npan  = npan - 1
          go to 445
        end if
        sden =                                      ! normale esterna !
     & ((y4(i)-y2(i))*(z3(i)-z1(i))-(y3(i)-y1(i))*(z4(i)-z2(i)))**2 +
     & ((z4(i)-z2(i))*(x3(i)-x1(i))-(z3(i)-z1(i))*(x4(i)-x2(i)))**2 +
     & ((x4(i)-x2(i))*(y3(i)-y1(i))-(x3(i)-x1(i))*(y4(i)-y2(i)))**2
        if(sden.ne.0) then
          sden = 1. / sqrt(sden)
        else
          do ii=i,n_est-1
            x1(ii)=x1(ii+1)
            x2(ii)=x2(ii+1)
            x3(ii)=x3(ii+1)
            x4(ii)=x4(ii+1)
            y1(ii)=y1(ii+1)
            y2(ii)=y2(ii+1)
            y3(ii)=y3(ii+1)
            y4(ii)=y4(ii+1)
            z1(ii)=z1(ii+1)
            z2(ii)=z2(ii+1)
            z3(ii)=z3(ii+1)
            z4(ii)=z4(ii+1)
          end do
          n_est = n_est - 1
          npan  = npan - 1
          if(i.gt.n_est) go to 446
          go to 445
        end if
        cosxn(i) = + sden*
     & ((y4(i)-y2(i))*(z3(i)-z1(i))-(y3(i)-y1(i))*(z4(i)-z2(i)))
        cosyn(i) = + sden*
     & ((z4(i)-z2(i))*(x3(i)-x1(i))-(z3(i)-z1(i))*(x4(i)-x2(i)))
        coszn(i) = + sden*
     & ((x4(i)-x2(i))*(y3(i)-y1(i))-(x3(i)-x1(i))*(y4(i)-y2(i)))
       xl = sqrt(cosxn(i)**2.+cosyn(i)**2.+coszn(i)**2.)
       cosxn(i)=cosxn(i)/xl
       cosyn(i)=cosyn(i)/xl
       coszn(i)=coszn(i)/xl
 446  end do
c
      if (numfilint.gt.0) then
       do i=n_est+1,npan
 447    if(x1(i).eq.x3(i).and.y1(i).eq.y3(i).and.z1(i).eq.z3(i)) then
          do ii=i,npan-1
            x1(ii)=x1(ii+1)
            x2(ii)=x2(ii+1)
            x3(ii)=x3(ii+1)
            x4(ii)=x4(ii+1)
            y1(ii)=y1(ii+1)
            y2(ii)=y2(ii+1)
            y3(ii)=y3(ii+1)
            y4(ii)=y4(ii+1)
            z1(ii)=z1(ii+1)
            z2(ii)=z2(ii+1)
            z3(ii)=z3(ii+1)
            z4(ii)=z4(ii+1)
          end do
          npan = npan - 1
          if(i.gt.npan) go to 448
          go to 447
        end if
        sden =                                       ! normale interna !
     & ((y4(i)-y2(i))*(z3(i)-z1(i))-(y3(i)-y1(i))*(z4(i)-z2(i)))**2 +
     & ((z4(i)-z2(i))*(x3(i)-x1(i))-(z3(i)-z1(i))*(x4(i)-x2(i)))**2 +
     & ((x4(i)-x2(i))*(y3(i)-y1(i))-(x3(i)-x1(i))*(y4(i)-y2(i)))**2 
        if(sden.eq.0) then
          do ii=i,npan-1
            x1(ii)=x1(ii+1)
            x2(ii)=x2(ii+1)
            x3(ii)=x3(ii+1)
            x4(ii)=x4(ii+1)
            y1(ii)=y1(ii+1)
            y2(ii)=y2(ii+1)
            y3(ii)=y3(ii+1)
            y4(ii)=y4(ii+1)
            z1(ii)=z1(ii+1)
            z2(ii)=z2(ii+1)
            z3(ii)=z3(ii+1)
            z4(ii)=z4(ii+1)
          end do
          npan = npan - 1
          go to 447
        end if
        sden = -1. / sqrt(sden)
        cosxn(i) = + sden*
     & ((y4(i)-y2(i))*(z3(i)-z1(i))-(y3(i)-y1(i))*(z4(i)-z2(i)))
        cosyn(i) = + sden*
     & ((z4(i)-z2(i))*(x3(i)-x1(i))-(z3(i)-z1(i))*(x4(i)-x2(i)))
        coszn(i) = + sden*
     & ((x4(i)-x2(i))*(y3(i)-y1(i))-(x3(i)-x1(i))*(y4(i)-y2(i)))
       xl = cosxn(i)**2.+cosyn(i)**2.+coszn(i)**2.
       if(xl.ne.0) xl = sqrt(xl)
       cosxn(i)=cosxn(i)/xl
       cosyn(i)=cosyn(i)/xl
       coszn(i)=coszn(i)/xl
 448  end do
      end if
c
      xmaximm = -1E+33
      do i=1,npan
        xmp    =max(x1(i),x2(i),x3(i),x4(i))
        xmaximm=max(xmaximm,xmp)
      end do

c----------- Verifica ---------
      if (geofile.eq.'y') then
c------------------------------ 
      open (99, file='norm.dat',form='formatted',status='unknown',
     &         recl=64000)
      write(99,*) 'VARIABLES="x","y","z","vx","vy","vz"'
      write(99,*) 'ZONE T="Normali" I=',npan
      do i=1,npan
        write(99,*) xn(i),(yn(i)+yswath/xtot),zn(i),
     &              cosxn(i),cosyn(i),coszn(i)
      end do
      close(99)
      end if
c----------- Fine verifica ----
c
c---- E' un file necessario per LITRA & SWARP
c
c
c---- E' una carena SWATH ?
c
      if (swath.eq.'y') then
c
c Ordina i punti x_0,y_0 da prua a poppa
c
c         CALL PIKSR2(lat_est,x_0,y_0)
         nswath=1
         if (sym.eq.'y') then
         nswath=2
            do i = 1,lat_est
               x_0(i+lat_est) = x_0(i)
               y_0(i+lat_est) = -y_0(i)
            end do
         end if
      else
         nswath=1
c         CALL PIKSR2(lat_est,x_0,y_0)
      end if
c
c   Passa al profilo interno
c
      if (numfilint.gt.0) then
c---- E' una carena SWATH ?
c
         if (swath.eq.'y') then
c
c Ordina i punti x_1,y_1 da prua a poppa
c
            nswath=1
c            CALL PIKSR2(lat_int,x_1,y_1)
            if (sym.eq.'y') then
            nswath=2
               do i = 1,lat_int
                  x_1(i+lat_int) = x_1(i)
                  y_1(i+lat_int) = - y_1(i)
               end do
            end if
         else
            nswath=1
c            CALL PIKSR2(lat_int,x_1,y_1)
         end if
c
      end if

      lat_0=lat_est+lat_int
c
      if (swath.eq.'n'.and.lontra.eq.0) y_0(lat_0) = 0.

c      print *, lat_est, lat_int, lat_0
c      print *, y_0(1), y_1(1)

      open(unit=66,file='carena_0.dat',status='unknown')
      write(66,*) lat_0*nswath
c-- md
c      do i = 1,lat_est*nswath
c        write(66,*) x_0(i),y_0(i)+yswath/xtot
c      end do

      write(66,*) x_0(1),y_0(1)+yswath/xtot
      do i = 2,lat_est-1
        if (y_0(i).lt.y_0(1)) y_0(i)=y_0(1) !+ 1.e-6
        write(66,*) x_0(i),y_0(i)+yswath/xtot
      end do
      write(66,*) x_0(lat_est),y_0(lat_est)+yswath/xtot

      if (sym.eq.'y') then
        write(66,*) x_0(lat_est+1),y_0(lat_est+1)+yswath/xtot
        do i = lat_est+2,lat_est*nswath-1
          if (y_0(i).gt.y_0(1)) y_0(i)=y_0(1) !- 1.e-6
          write(66,*) x_0(i),y_0(i)+yswath/xtot
        end do
        write(66,*) x_0(lat_est*nswath),y_0(lat_est*nswath)+yswath/xtot
      endif
c-- md

      ltot = x_0(lat_est*nswath)-x_0(1)
      if (numfilint.gt.0) then

c--md
c         do i = 1,lat_int*nswath
c           write(66,*) x_1(i),y_1(i)+yswath/xtot
c         end do
      write(66,*) x_1(1),y_1(1)+yswath/xtot
      do i = 2,lat_int-1
        if (y_1(i).gt.y_1(1)) y_1(i)=y_1(1) !+ 1.e-6
        write(66,*) x_1(i),y_1(i)+yswath/xtot
      end do
      write(66,*) x_1(lat_int),y_1(lat_int)+yswath/xtot
c--md
      end if
      close(66)

      deallocate(filenamei,filenameo)

c-----
      return
      end 
c ----------------------------------------------------------------------

      subroutine scrivigriglia(appox,appoy,appoz,nodi_i,nodi_j,ndom,
     &                         xtot,yswath)

c ----------------------------------------------------------------------

      parameter(ksez=200,klin=200)
      real appox(ksez,klin),appoy(ksez,klin),appoz(ksez,klin)
      character*3  cha
      character*80 namepres

c ----------------------------------------------------------------------

      iuno = 1
      write(cha,99) ndom
      namepres='gri'//cha//'.tru'
      open(67,file=namepres,form='formatted',status='unknown')
      write(67,*) iuno
      write(67,*) nodi_i,nodi_j,iuno
      write(67,*) (( appox(i,j)/xtot,i=1,nodi_i),j=1,nodi_j),
     &            (((appoy(i,j)+yswath)/xtot,i=1,nodi_i),j=1,nodi_j),
     &            (( appoz(i,j)/xtot,i=1,nodi_i),j=1,nodi_j)
      close(67)
      write(cha,99) ndom
      namepres='gri'//cha//'.dat'
      open (67,file=namepres,form='formatted',status='unknown')
      write(67,*) 'VARIABLES="x","y","z"'
      write(67,*) 'ZONE T="Griglia" I=',nodi_j,'J=',nodi_i
      do i=1,nodi_i
        do j=1,nodi_j
          write(67,*) (appox(i,j)/xtot),
     &               ((appoy(i,j)+yswath)/xtot),
     &                (appoz(i,j)/xtot)
        end do
      end do
      write(67,*) 'ZONE T="J=1" I=',nodi_i,'J=2'
      do i=1,nodi_i
        write(67,*) (appox(i,1)/xtot),
     &             ((appoy(i,1)+yswath)/xtot),
     &              (appoz(i,1)/xtot)
      end do
      do i=1,nodi_i
        write(67,*) (appox(i,2)/xtot),
     &             ((appoy(i,2)+yswath)/xtot),
     &              (appoz(i,2)/xtot)
      end do
      write(67,*) 'ZONE T="I=1" I=',nodi_j,'J=2'
      do i=1,nodi_j
        write(67,*) (appox(1,i)/xtot),
     &             ((appoy(1,i)+yswath)/xtot),
     &              (appoz(1,i)/xtot)
      end do
      do i=1,nodi_j
        write(67,*) (appox(2,i)/xtot),
     &             ((appoy(2,i)+yswath)/xtot),
     &              (appoz(2,i)/xtot)
      end do
      close(67)

   99 format(i3.3)

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine scrivigrigliasimm(appox,appoy,appoz,nodi_i,nodi_j,
     &                             ndom,xtot,yswath)

c ----------------------------------------------------------------------

      parameter(ksez=200,klin=200)
      real appox(ksez,klin),appoy(ksez,klin),appoz(ksez,klin)
      character*3  cha
      character*80 namepres

c ----------------------------------------------------------------------

      iuno = 1
      write(cha,99) ndom
      namepres='gri'//cha//'.tru'
      open(67,file=namepres,form='formatted',status='unknown')
      write(67,*) iuno
      write(67,*) nodi_i,nodi_j,iuno
      write(67,*) ( ( appox(i,j)/xtot,i=1,nodi_i),j=1,nodi_j),
     &            (((-appoy(i,j)+yswath)/xtot,i=1,nodi_i),j=1,nodi_j),
     &            ( ( appoz(i,j)/xtot,i=1,nodi_i),j=1,nodi_j)
      close(67)
      write(cha,99) ndom
      namepres='gri'//cha//'.dat'
      open (67,file=namepres,form='formatted',status='unknown')
      write(67,*) 'VARIABLES="x","y","z"'
      write(67,*) 'ZONE T="Griglia" I=',nodi_j,'J=',nodi_i
      do i=1,nodi_i
        do j=1,nodi_j
          write(67,*) ( appox(i,j)/xtot),
     &               ((-appoy(i,j)+yswath)/xtot),
     &                ( appoz(i,j)/xtot)
        end do
      end do
      write(67,*) 'ZONE T="J=1" I=',nodi_i,'J=2'
      do i=1,nodi_i
        write(67,*) ( appox(i,1)/xtot),
     &             ((-appoy(i,1)+yswath)/xtot),
     &              ( appoz(i,1)/xtot)
      end do
      do i=1,nodi_i
        write(67,*) ( appox(i,2)/xtot),
     &             ((-appoy(i,2)+yswath)/xtot),
     &              ( appoz(i,2)/xtot)
      end do
      write(67,*) 'ZONE T="I=1" I=',nodi_j,'J=2'
      do i=1,nodi_j
        write(67,*) ( appox(1,i)/xtot),
     &             ((-appoy(1,i)+yswath)/xtot),
     &              ( appoz(1,i)/xtot)
      end do
      do i=1,nodi_j
        write(67,*) ( appox(2,i)/xtot),
     &             ((-appoy(2,i)+yswath)/xtot),
     &              ( appoz(2,i)/xtot)
      end do
      close(67)

   99 format(i3.3)

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine leggigriglia(appox,appoy,appoz,ni,nj,ndom)

c ----------------------------------------------------------------------

      parameter(ksez=200,klin=200)
      real appox(ksez,klin),appoy(ksez,klin),appoz(ksez,klin)
      character*3  cha
      character*80 namepres

c ----------------------------------------------------------------------

      write(cha,99) ndom
      namepres='gri'//cha//'.tru'
      open(67,file=namepres,form='formatted',status='unknown')
      read(67,*) nb
      read(67,*) ni,nj,nk
      read(67,*) ((appox(i,j),i=1,ni),j=1,nj),
     &           ((appoy(i,j),i=1,ni),j=1,nj),
     &           ((appoz(i,j),i=1,ni),j=1,nj)
      close(67)

   99 format(i3.3)

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine reinterpola(x,y,np)

c ----------------------------------------------------------------------

      real s(1000),x(*),y(*)
      real si(100),xi(100),yi(100)

c ----------------------------------------------------------------------

      if(np.eq.0) then
        np = 2
        x(1) = -0.5
        y(1) =  0.0
        x(2) =  0.5
        y(2) =  0.0
      end if

      s(1) = 0.
      do i=2,np
        dx = x(i)-x(i-1)
        dy = y(i)-y(i-1)
        ds = dx*dx+dy*dy
        if(ds.ne.0) ds = sqrt(ds)
        s(i) = s(i-1)+ds
      end do

      ni = 100
      do i=1,ni
        si(i) = s(np)*(i-1.)/(ni-1.)
      end do

      call interpola(s,x,np,si,xi,ni)
      call interpola(s,y,np,si,yi,ni)

      np = ni
      do i=1,ni
        x(i) = xi(i)
        y(i) = yi(i)
      end do

c ----------------------------------------------------------------------
      return
      end 
c ----------------------------------------------------------------------
