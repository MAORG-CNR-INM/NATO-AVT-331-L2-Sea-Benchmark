c-----------------------------------------------------------------------

      subroutine ittcnl(fr,cw,sink,atrim,ztrasla,
     &                  Cf_ittc,Cf,Ct_ittc,Ct,Cyd,
     &                  sup,shiplen,area,
     &                  fi0x,fi0y,fi0z,
     &                  fix,fiy,fiz,presto,
     &                  xn,yn,zn,
     &                  cosxn,cosyn,coszn,
     &                  ncar,nplib,ncard,ntotd)

c-----------------------------------------------------------------------

      implicit real (a-h,o-z)
c      real ro,nu
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd)
      real fix(ntotd),fiy(ntotd),fiz(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real area(ntotd)
      real presto(ncard)
      character*9 fina
      character*11 finanl
      logical simm

c-----------------------------------------------------------------------
c - as
      real cfx,reyx
      real ro,nu,gi

      namelist/FLUID_PARAMETERS/ ro,nu,gi

      open(unit=20,file='SBDF.nml',action='read',recl=10000)
        read(20,NML=FLUID_PARAMETERS)
      close(20)
c - as

      simm = .true.

c-----------------------------------------------------------------------
c   Calcolo di Cf e Ct
c-----------------------------------------------------------------------

      cvp = 0.
      rv  = 0.
c      cfx = 0.
c      vinf= fr*sqrt(shiplen*9.81)
c      vinf= fr*sqrt(shiplen*9.8033)
      vinf = fr*sqrt(shiplen*gi)

      do i=1,ncar
        u     = (1.0+fi0x(i)+fix(i))*vinf
        v     = (    fi0y(i)+fiy(i))*vinf
        w     = (    fi0z(i)+fiz(i))*vinf
        vel   = sqrt(u**2+v**2+w**2) 
        rey   = vel*shiplen/nu
        areal2= area(i)*shiplen**2
        rvi   = 0.455*(log10(rey)**(-2.58)) * ro/2.*vel**2.*areal2
        rv    = rv + 2.*rvi*u/vel
c - as
c        reyx  = vel*area(i)*shiplen/nu
c        cfx   = cfx + ((2*log10(rey)-0.65)**(-2.3))*ro*vel*areal2*u
c - as
      end do

c -- Calcolo Cf secondo formula ITTC

      Cf_ittc = 0.075/(log10(vinf*abs(shiplen)/nu)-2.)**2

c -- Calcolo Cf secondo Schlichting

      Cf      = rv/0.5/ro/(sup*abs(shiplen)**2)/vinf**2

c -- as
c      rey = vinf*shiplen/nu
cc -- Prandtl-Schlichting (1932)
c      Cf_PS   = 0.455*((log10(rey))**(-2.58))
cc -- Schlichting (1979)
c      Cf_S    = cfx/0.5/ro/(sup*abs(shiplen)**2)/vinf**2
cc --
c      open(20,file='friction.out',status='unknown',access='append',
c     &        recl=10000)
c       write(20,*) 'Fr,Cf_ittc,Cf,Cf_PS,Cf_S'
c       write(20,*)  fr,Cf_ittc,Cf,Cf_PS,Cf_S
c      close(20)
c -- as

c -- Calcolo Ct 

      Ct      = Cw + Cf      !- Cyd
      Ct_ittc = Cw + Cf_ittc !- Cyd

      close(1)
c ----------------------------------------------------------------------
c   File di uscita
c ----------------------------------------------------------------------

      fina(1:5) = 'ittc.'
      nufru = nint((fr)*1000)
      write(fina(6:9),'(i4.4)') nufru
c
      open(99,file=fina,form='unformatted',status='unknown')
         write(99) fr,cw,sup,shiplen,ncar,cyd,sink,atrim,ztrasla
         write(99) (fi0x(i),i=1,ncar),
     &             (fi0y(i),i=1,ncar),
     &             (fi0z(i),i=1,ncar),
     &             (fix(i) ,i=1,ncar),
     &             (fiy(i) ,i=1,ncar),
     &             (fiz(i) ,i=1,ncar),
     &             (area(i),i=1,ncar)
      close(99)

      finanl(1:7) = 'ittcnl.'
      nufru = nint((fr)*1000)
      write(finanl(8:11),'(i4.4)') nufru

      cfy = 0.
      cfz = 0.
      cmx = 0.
      cmy = 0.
      cmz = 0.

      open(99,file=finanl,form='unformatted',status='unknown')
      write(99) fr,cw,cyd,cfy,cfz,cmx,cmy,cmz,sink,atrim,ztrasla,simm,
     &          sup,shiplen,ncar,1,ncar
      write(99) (fi0x(i)+fix(i) ,i=1,ncar),
     &          (fi0x(i)+fiy(i) ,i=1,ncar),
     &          (fi0x(i)+fiz(i) ,i=1,ncar),
     &          (0.0            ,i=1,ncar),
     &          (0.0            ,i=1,ncar),
     &          (0.0            ,i=1,ncar),
     &          (xn(i)          ,i=1,ncar),
     &          (yn(i)          ,i=1,ncar),
     &          (zn(i)          ,i=1,ncar),
     &          (cosxn(i)       ,i=1,ncar),
     &          (cosyn(i)       ,i=1,ncar),
     &          (coszn(i)       ,i=1,ncar),
     &          (presto(i)      ,i=1,ncar),
     &          (1.0            ,i=1,ncar),
     &          (4              ,i=1,ncar),
     &          (area(i)        ,i=1,ncar)
      close(99)

c-----------------------------------------------------------------------
      return 
      end
c-----------------------------------------------------------------------
