c-----------------------------------------------------------------
      subroutine freesurface(igrid,ngrid,ntraint)
c-----------------------------------------------------------------

      integer ighost,ntraint,antraest(ngrid),ancarlo(ngrid),
     &                       anmonte(ngrid),anvalle(ngrid)
      integer ntraest,ncarlo,nmonte,nvalle
      integer iwrite,iunit
      real    xstart,xend,ystart,yend,beta,strey

      real xfi(500,500),yfi(500,500)
      real xfs(500,500),yfs(500,500)
      real xtr(500,500),ytr(500,500)
      real x(500,500),y(500,500)

      real xn(1000),yn(1000),xs(1000),ys(1000)
      real xe(1000),ye(1000),xo(1000),yo(1000)

      real x0e(1000),y0e(1000),x0i(1000),y0i(1000)

      real sc(1000),xc(1000),yc(1000)
      real si(1000),xi(1000),yi(1000)

      real xstern(1000),ystern(1000),zstern(1000)
      real xdot(1),zdot(1)

      character*60 carena

      real tens
      common/tensione/tens

      PARAMETER(N=20) 
      COMMON/CONST/CHOICE,NA,DS1,DS2,NI,NJ
      COMMON/ATTR/IAL(N),IAX(N),IAY(N),JAL(N),JAX(N),JAY(N) 
      COMMON/COEF/AI(N),BI(N),CI(N),DI(N),AJ(N),BJ(N) 
      COMMON/COEF/CJ(N),DJ(N) 
      INTEGER CHOICE 

      namelist/FREE_SURFACE/ ighost,ntraint,antraest,ancarlo,anmonte,
     &                       anvalle,iwrite,xstart,xend,ystart,yend,
     &                       beta,strey

c-----------------------------------------------------------------

      iunit=20

      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=FREE_SURFACE)
      close(iunit)

      ntraest = antraest(igrid)
      ncarlo  = ancarlo(igrid)
      nmonte  = anmonte(igrid)
      nvalle  = anvalle(igrid)

      open(72,file='confine.grd',form='formatted',
     &        status='old',recl=64000)
      read(72,200) carena
      read(72,*) npancar,lontra,nvort,nscia,ntip,
     &           x1,x2,x3,x4,x5,x6,x7,x8
      close(72)
      if(npancar.eq.0) stop

      pig = acos(-1.0)

c      open(61,file='free-surface.inp',status='old')
c      read(61,*)
c      read(61,*)
c      read(61,*)
c
c      read(61,*) ntraint
c      read(61,*) ntraest
c      read(61,*) ncarlo
c      read(61,*) nmonte
c      read(61,*) nvalle
c      read(61,*)
c      read(61,*) ighost
c      read(61,*)
c      read(61,*) xstart
c      read(61,*) xend
c      read(61,*) ystart
c      read(61,*) yend
c      read(61,*)
c      read(61,*) beta
      beta = beta*pig/180.
      itrx = 1
c      read(61,*)
c      read(61,*) iwrite
c      read(61,*)
c      read(61,*) strey
c      close(61)

      open(1,file='carena_0.dat',status='old')
      if(ntraint.eq.0) then
        read(1,*) nde
        do i=1,nde
          read(1,*) x0e(i),y0e(i)
        end do
        ymed = 0.
      else
        read(1,*) nde
        nde = nde/2
        ndi = nde
        do i=1,nde
          read(1,*) x0e(i),y0e(i)
        end do
        ymini = 1E+33
        do i=1,ndi
          read(1,*) x0i(i),y0i(i)
          if(x0i(i).lt.0.4.and.x0i(i).gt.-0.4) ymini = min(ymini,y0i(i))
        end do
        nlo = nmonte+ncarlo+nvalle
        dx = (xend-xstart)/float(nlo)
c--md        ntraint = int(ymini/dx/0.95)+1
!        ntraint = floor(ymini/(dx*1.010))+1
        ntraint = max(1,ntraint)
c--md
c        if(ntraint.lt.1) ntraint = 1
        x0e(1) = (x0e(1)+x0i(1))/2.
        y0e(1) = (y0e(1)+y0i(1))/2.
        ymed = y0e(1)
        x0i(1) = x0e(1)
        y0i(1) = y0e(1)
      end if
      close(1)

      if(ntraint.ne.0) then
        do i=nde/2,nde
          if(y0i(i).gt.y0e(i)) then
            appo   = y0i(i)
            y0i(i) = y0e(i)
            y0e(i) = appo
          end if
        end do
      end if

      nvalle  = nvalle  + 1
      ntraest = ntraest + 1
      if(ntraint.ne.0) ntraint = ntraint + 1
      ntra = ntraest + ntraint
      ntralo = nvalle

c -- Verifica se serve il transom

      if(ighost.eq.0) then
        open(44,file='transom.grd',status='old',err=4)
        close(44)
      end if

      if(lontra.eq.0.and.ighost.eq.0) then
        ntrans = 0
        if(ntraint.ne.0) then
          y0e(nde) = (y0e(nde)+y0i(ndi))/2.
          y0i(ndi) = y0e(nde)
          ymedtrans = y0e(nde)
        else
          y0e(nde) = 0.
          ymedtrans = y0e(nde)
        end if
      else
        x0 = x0e(nde)
        x1 = xend
        y0 = 0.0
        y1 = y0e(nde)
c--md        dx = (x1-x0)/float(nvalle)*0.95
        dx = (x1-x0)/float(nvalle)
c--md
c        if(ntraint.ne.0) then
c          y0 = y0i(ndi)
c          dx = y0/float(ntraint)
c        end if
c--md         ntrans = int((y1-y0)/dx)
       ntrans = floor((y1-y0)/(dx*1.005))
       ntrans = max(ntrans,1)
c--md      
        if(ntrans.eq.0) ntrans = ighost
        if(ntrans.eq.0) then
          if(ntraint.ne.0) then
            y0e(nde) = (y0e(nde)+y0i(ndi))/2.
            y0i(ndi) =  y0e(nde)-float(ighost)*dx/2.
            y0e(nde) = y0e(nde)+float(ighost)*dx/2.
            ntrans = ighost
            do i=3*nde/4,nde-1
              if(y0e(i).lt.y0e(nde)) y0e(i) = y0e(nde)
            end do
            do i=3*ndi/4,ndi-1
              if(y0i(i).gt.y0i(nde)) y0i(i) = y0i(nde)
            end do
          else
            if(ighost.eq.0) then
              y0e(nde) = 0.
            else
              y0e(nde) = y0e(nde)+float(ighost)*dx
              ntrans = ighost
              do i=3*nde/4,nde-1
                if(y0e(i).lt.y0e(nde)) y0e(i) = y0e(nde)
              end do
            end if
          end if
        end if
        ymedtrans = y0e(nde)
      end if

c --

   4  nlo = nmonte + ncarlo + nvalle
      ntotp = (ntra-1)*(nlo-1) + (ntralo-1)*(ntrans-1)

      if(ntrans.ne.0) then
        do i=3*nde/4,nde-1
          if(y0e(i).lt.y0e(nde)) y0e(i) = y0e(nde)
        end do
        if(ntraint.ne.0) then
          do i=3*ndi/4,ndi-1
            if(y0i(i).gt.y0i(nde)) y0i(i) = y0i(nde)
          end do
        end if
      end if

      if(ntrans.eq.0) then
        if(ntraint.eq.0) then
          if(ighost.ne.0) then
            ntrans = ighost
            y0e(nde) = dx*float(abs(ighost))
            do i=3*nde/4+1,nde-1
              if(y0e(i).lt.y0e(nde)) y0e(i)= y0e(nde)
            end do
          else
            y0e(nde) = 0.0
          end if
        else
          ymed = (y0e(nde)+y0i(nde))/2.
          if(ighost.ne.0) then
            y0e(nde) = ymed+dx/2.*float(abs(ighost))
            y0i(ndi) = ymed-dx/2.*float(abs(ighost))
          else
            y0e(nde) = ymed
            y0i(ndi) = ymed
          end if
        end if
      end if
      if(ntraint.ne.0) then
        yprua = (y0e(1)+y0i(1))/2.
      else
        yprua = 0.
      end if

c-----------------------------------------------------------------
c -- Modifica un po' di cose per il transom fittizio
c-----------------------------------------------------------------

      if(ighost.ne.0.or.ntrans.ne.0) then

        open(44,file='transom.grd',status='old',err=5)
        close(44)
        go to 7

   5    call buildtrans(x0e,nde,ntrans)

   7    continue

        if(ighost.eq.0) then
          if(ntrans.gt.0) then
            open(44,file='transom.grd',status='old')
            read(44,*) x0e(nde)
            x0e(nde) = x0e(nde) + 0.001
            if(ntraint.ne.0) x0i(ndi) = x0e(nde)
            close(44)
          end if
        else if(ighost.ne.0) then

          call buildtrans(x0e,nde,ntrans)

          open(72,file='confine.grd',form='formatted',
     &            status='old',recl=64000)
          open(73,file='confine.new',form='formatted',
     &            status='unknown',recl=4000)
          read(72,200) carena
          write(73,200) carena

          if(itrx.eq.1) then
            read(72,*)  npancar,lontra,nvort,nscia,ntip,
     &                x1,x2,x3,x4,x5,x6,x7,x8
            write(73,250) npancar,2,nvort,nscia,ntip,
     &                  x1,x2,x3,x4,x5,x6,x7,x8
          else
            read(72,*)  npancar,lontra,
     &                x1,x2,x3,x4,x5,x6,x7,x8
            write(73,260) npancar,2,
     &                  x1,x2,x3,x4,x5,x6,x7,x8
          end if

          dxmin = 1000.
          do i=1,npancar
            read(72,*)  x1,x2,x3,x4,x5,x6,x7
            write(73,300) x1,x2,x3,x4,x5,x6,x7
          end do

          do i=1,npancar
            read(72,*)  x1,x2,x3,x4,x5,x6,x7
            write(73,300) x1,x2,x3,x4,x5,x6,x7
          end do

          do i=1,npancar
            read(72,*)  x1,x2,x3,x4,x5,x6,x7
            write(73,300) x1,x2,x3,x4,x5,x6,x7
          end do

          do i=1,npancar
            read(72,*)  x1,idum
            write(73,*) x1,idum
          end do

          do i=1,nscia+ntip
            read(72,*)  idum1,idum2
            write(73,400) idum1,idum2
            read(72,*)  x1,x2,x3
            write(73,300) x1,x2,x3
            read(72,*)  x1,x2,x3
            write(73,300) x1,x2,x3
          end do

          close(73)
          close(72)
          call system('mv confine.new confine.grd')
        end if

      end if

c-----------------------------------------------------------------
c -- Zona esterna (c'e' sempre)
c-----------------------------------------------------------------

      if(ntraint.eq.0) then
        ntra = ntraest
        call profiloe(ntraest,nlo,ncarlo,nmonte,nvalle,
     &                x0e,y0e,nde,
     &                xstart,xend,ystart,yend,
     &                xn,yn,xs,ys,xe,ye,xo,yo)
      else
        call profiloe(ntraest,nlo,ncarlo,nmonte,nvalle,
     &                x0e,y0e,nde,
     &                xstart,xend,y0e(1),yend,
     &                xn,yn,xs,ys,xe,ye,xo,yo)
      end if

      do i=1,nlo
        x0 = xs(i)
        y0 = ys(i)
        x1 = xn(nlo-i+1)
        y1 = yn(nlo-i+1)
        call congiungi(x0,y0,x1,y1,xc,yc,ntraest)
        do j=1,ntraest
          xfs(i,j) = xc(j)
          yfs(i,j) = yc(j)
        end do
      end do

c--md
      if(x7.eq.0.) then ! la carena è interamente sommersa
        do i=1,nlo
          do j=1,ntraest
            xfs(i,j) = xstart+(i-1.)/(nlo-1.)*(xend-xstart)
            yfs(i,j) = ystart+(j-1.)/(ntraest-1.)*(yend-ystart)
          enddo
        enddo
      endif
c--md

c-----------------------------------------------------------------
c -- Zona interna (se c'e')
c-----------------------------------------------------------------

      if(ntraint.ne.0) then
        call profiloi(ntraint,nlo,ncarlo,nmonte,nvalle,
     &                x0i,y0i,ndi,
     &                xstart,xend,ystart,y0e(1),
     &                xn,yn,xs,ys,xe,ye,xo,yo)

        do i=1,nlo
          x0 = xs(i)
          y0 = ys(i)
          x1 = xn(i)
          y1 = yn(i)
          call congiungi(x0,y0,x1,y1,xc,yc,ntraint)
          do j=1,ntraint
            xfi(i,j) = xc(j)
            yfi(i,j) = yc(j)
          end do
        end do

      end if

c-----------------------------------------------------------------
c -- Transom
c-----------------------------------------------------------------

  8   ntotp = (ntra-1)*(nlo-1) + (ntralo-1)*(ntrans-1)

      if(ntrans.ne.0) ntrans  = ntrans  + 1

  10  if(ntrans.ne.0) then

        k = 0
        do i=nlo-nvalle,nlo
          k = k + 1
          x0 = xfs(i,1)
          x1 = xfs(i,1)
          y0 = yfs(i,1)
          y1 = 0
          if(ntraint.ne.0) y1 = yfi(i,ntraint)
          call congiungi(x0,y0,x1,y1,xc,yc,ntrans)
          do j=1,ntrans
            xtr(k,j) = xc(j)
            ytr(k,j) = yc(j)
          end do
        end do
        ntralo = k

      end if

c -- Stretch

c -- Esterno

      yhmax = 0.
      imax  = 1

      do i=1,nlo
!--md        if(yfs(i,1).gt.yhmax) then
        if(yfs(i,1).ge.yhmax) then
          yhmax = yfs(i,1)
          imax  = i
        end if
      end do

      ymax = yfs(1,ntraest)
      dx   = 0.9*(xfs(imax,1)-xfs(imax-1,1))*cos(beta)
      jj   = 1
      yy   = yhmax+dx
      do while(yy.le.ymax)
        dx = dx*strey
        yy = yy + dx
        jj = jj + 1
      end do
      
      ntraest = jj-1

      do i=1,nlo
        x0 = xfs(i,1)
        y0 = yfs(i,1)
        if(i.eq.1) then
          dx = xfs(2,1)-xfs(1,1)
        else
          dx = xfs(i,1)-xfs(i-1,1)
        end if
        do j=2,ntraest
          xfs(i,j) = xfs(i,1)
          yfs(i,j) = yfs(i,j-1)+dx
          dx = dx*strey
        end do
      end do

c -- Riporta tutti al bordo esterno

      do i=1,nlo
        scalay = (ymax-yfs(i,1))/(yfs(i,ntraest)-yfs(i,1))
        do j=2,ntraest
          yfs(i,j) = yfs(i,1)+(yfs(i,j)-yfs(i,1))*scalay
        end do
      end do

c -- Inclinazione

      do j=1,ntraest
        do i=1,nlo
          xfs(i,j) = xfs(i,j) + (yfs(i,j)-yfs(i,1))*tan(beta)
        end do
      end do

C_DPE
c      call angles(xfs,yfs,nlo,ntraest,1)
      call tm(xfs,yfs,nlo,ntraest)
c      call laplace(xfs,yfs,nlo,ntraest)

c -- Scrittura file

      open(1,file='fs.plt',status='unknown')
      open(2,file='GNUfs.dat',status='unknown')
      write(1,*) 'VARIABLES = "x" "y" "z"'
      write(1,*) 'ZONE I=',ntraest,' J=',nlo
      do i=1,nlo
        do j=1,ntraest
          write(1,*) xfs(i,j),yfs(i,j),0.0
          write(2,*) xfs(i,j),yfs(i,j),0.0
        end do
        write(2,*)
      end do
      write(2,*)
      if(ntraint.gt.1) then
        write(1,*) 'ZONE I=',ntraint,' J=',nlo
        do i=1,nlo
          do j=1,ntraint
            write(1,*) xfi(i,j),yfi(i,j),0.0
            write(2,*) xfi(i,j),yfi(i,j),0.0
          end do
          write(2,*)
        end do
        write(2,*)
      end if
      if(ntrans.gt.1) then
        write(1,*) 'ZONE I=',ntrans,' J=',ntralo
        do i=1,ntralo
          do j=1,ntrans
            write(1,*) xtr(i,j),ytr(i,j),0.0
            write(2,*) xtr(i,j),ytr(i,j),0.0
          end do
          write(2,*)
        end do
        write(2,*)
      end if
      close(2)
      close(1)

      if(ntrans.eq.0) ntrans = 1
      if(ntraint.eq.0) ntraint = 1
      open(50,file='fsgrid.grd',status='unknown')
      if(itrx.eq.1) then
        write(50,100) ntraint-1,ntraest-1,nlo-1,nvalle-1,
     &                ntrans-1,ntralo-1,0,0.,0.
      else
        write(50,110) ntraint-1,ntraest-1,nlo-1,nvalle-1,
     &                ntrans-1,0,0.,0.
      end if
      if(ntraint.lt.1) then
        do i=1,nlo-1
          do j=1,ntraest-1
            write(50,70) xfs(i,j+1),xfs(i,j),xfs(i+1,j),xfs(i+1,j+1)
            write(50,70) yfs(i,j+1),yfs(i,j),yfs(i+1,j),yfs(i+1,j+1)
            write(50,70) 0.0       ,0.0         ,0.0       ,0.0
          end do
        end do
      else
        do i=1,nlo-1
          do j=1,ntraint-1
            write(50,70) xfi(i,j),xfi(i,j+1),xfi(i+1,j+1),xfi(i+1,j)
            write(50,70) yfi(i,j),yfi(i,j+1),yfi(i+1,j+1),yfi(i+1,j)
            write(50,70) 0.0       ,0.0         ,0.0       ,0.0
          end do
          do j=1,ntraest-1
            write(50,70) xfs(i,j+1),xfs(i,j),xfs(i+1,j),xfs(i+1,j+1)
            write(50,70) yfs(i,j+1),yfs(i,j),yfs(i+1,j),yfs(i+1,j+1)
            write(50,70) 0.0       ,0.0         ,0.0       ,0.0
          end do
        end do
      end if
      if(ntrans.gt.1) then
        do i=1,ntralo-1
          do j=1,ntrans-1
            write(50,70) xtr(i,j),xtr(i,j+1),xtr(i+1,j+1),xtr(i+1,j)
            write(50,70) ytr(i,j),ytr(i,j+1),ytr(i+1,j+1),ytr(i+1,j)
            write(50,70) 0.0       ,0.0         ,0.0       ,0.0
          end do
        end do
      end if
      close(50)


c -- Aspect ratio pannelli

      armin =  1E+33
      armax = -1E+33

      do i=1,nlo-1
        do j=1,ntraest-1
          ar = abs(yfs(i,j+1)-yfs(i,j))/abs(xfs(i+1,j)-xfs(i,j))
          armin = min(ar,armin)
          armax = max(ar,armax)
        end do
      end do

      do i=1,nlo-1
        do j=1,ntraint-1
          ar = abs(yfi(i,j+1)-yfi(i,j))/abs(xfi(i+1,j)-xfi(i,j))
          armin = min(ar,armin)
          armax = max(ar,armax)
        end do
      end do

c -- Output
      
      ntraint = 0      

      arg    = abs(xfs(2,1)-xfs(1,1))
      fr = 0.1
      nlam1  = int(2.*pig*fr**2/arg)
      nlamd1 = int(2.*pig*fr**2*cos(0.61365777)/arg)
      fr = 0.2
      nlam2  = int(2.*pig*fr**2/arg)
      nlamd2 = int(2.*pig*fr**2*cos(0.61365777)/arg)
      fr = 0.3
      nlam3  = int(2.*pig*fr**2/arg)
      nlamd3 = int(2.*pig*fr**2*cos(0.61365777)/arg)
      fr = 0.4
      nlam4  = int(2.*pig*fr**2/arg)
      nlamd4 = int(2.*pig*fr**2*cos(0.61365777)/arg)

      !if(iwrite.eq.1) then
        !write(*,*) '--------------------------------'
        !write(*,*) 'NTRAINT   = ',0
        !write(*,*) 'NTRAEST   = ',ntra-1
        !write(*,*) 'NLO       = ',nlo-1
        !write(*,*) 'NLOVALLEC = ',nvalle-1
        !write(*,*) 'NTRAPOC   = ',ntrans-1
        !write(*,*) '--------------------------------'
        !write(*,*) 'NTOT      = ',ntotp
        !write(*,*) '--------------------------------'
        !write(*,*) 'Aspect ratio min. :',armin
        !write(*,*) 'Aspect ratio max. :',armax
        !write(*,*) '--------------------------------'
        !if(ntrans-1.gt.0) then
          !write(*,*) 'Transom'
          !write(*,*) 'Aspect ratio min. :',atmin
          !write(*,*) 'Aspect ratio max. :',atmax
          !write(*,*) '--------------------------------'
        !end if
        !write(*,*) ' Fr = 0.1 '
        !write(*,*) ' Pannelli per lunghezza d''onda (L) = ',nlam1
        !write(*,*) ' Pannelli per lunghezza d''onda (T) = ',nlamd1
        !write(*,*) '--------------------------------'
        !write(*,*) ' Fr = 0.2 '
        !write(*,*) ' Pannelli per lunghezza d''onda (L) = ',nlam2
        !write(*,*) ' Pannelli per lunghezza d''onda (T) = ',nlamd2
        !write(*,*) '--------------------------------'
        !write(*,*) ' Fr = 0.3 '
        !write(*,*) ' Pannelli per lunghezza d''onda (L) = ',nlam3
        !write(*,*) ' Pannelli per lunghezza d''onda (T) = ',nlamd3
        !write(*,*) '--------------------------------'
        !write(*,*) ' Fr = 0.4 '
        !write(*,*) ' Pannelli per lunghezza d''onda (L) = ',nlam4
        !write(*,*) ' Pannelli per lunghezza d''onda (T) = ',nlamd4
        !write(*,*) '--------------------------------'
      !end if

   70 format(4f10.5)
  100 format(7i5,2f10.5)
  110 format(6i5,2f10.5)
  200 format(1x,a)
  250 format(5i5,8f10.5)
  260 format(2i5,8f10.5)
  300 format(8f10.5)
  400 format(i5)

c-----------------------------------------------------------------
      end
c-----------------------------------------------------------------

      subroutine profiloe(ntra,nlo,ncarlo,nmonte,nvalle,
     &                    xc,yc,nc,
     &                    xstart,xend,ystart,yend,
     &                    xn,yn,xs,ys,xe,ye,xo,yo)

c-----------------------------------------------------------------

      real xn(1000),yn(1000),xs(1000),ys(1000)
      real xe(1000),ye(1000),xo(1000),yo(1000)
      real xc(1000),yc(1000)

      real appox(1000),appoy(1000),s(1000),snew(1000)

      real tens
      common/tensione/tens

c-----------------------------------------------------------------
c -- Profilo nord
c-----------------------------------------------------------------

      x0 =  xend
      x1 =  xstart
      y0 =  yend
      y1 =  yend

      dx = x1-x0
      call stre(xn,nlo,dx,0)
      do i=1,nlo
        xn(i) = xn(i) + x0
      end do

      dy = y1-y0
      call stre(yn,nlo,dy,0)
      do i=1,nlo
        yn(i) = yn(i) + y0
      end do

c-----------------------------------------------------------------
c -- Profilo sud
c-----------------------------------------------------------------

c -- Prua

      x0 =  xstart
      x1 =  xc(1)
      y0 =  ystart
      y1 =  yc(1)

      dx = x1-x0
      dx1 = dx/float(nmonte)
      dx2 = (xc(nc)-xc(1))/float(ncarlo+1)/1.1
      call vinokur(xs,nmonte,dx,dx1,dx2)
      do i=1,nmonte
        xs(i) = xs(i) + x0
      end do

      dy = y1-y0
      call stre(ys,nmonte,dy,0)
      do i=1,nmonte
        ys(i) = ys(i) + y0
      end do

c -- Carena

      dx = xc(nc)-xc(1)
      dx1 = dx2
      dx2 = dx/float(ncarlo+1)/1.1
      call vinokur(appox,ncarlo+1,dx,dx1,dx2)
      do i=1,ncarlo+1
        appox(i) = appox(i) + xc(1)
      end do
      call interpola(xc,yc,nc,appox,appoy,ncarlo+1)

      do i=1,ncarlo
        xs(i+nmonte) = appox(i+1)
        ys(i+nmonte) = appoy(i+1)
      end do
      xs(nmonte+ncarlo) = xc(nc)
      ys(nmonte+ncarlo) = yc(nc)

c -- Poppa

      x0 =  appox(ncarlo+1)
      x1 =  xend
      y0 =  yc(nc)
      y1 =  yc(nc)

      dx = x1-x0
      dx1 = dx2
      dx2 = dx/float(nvalle+1)
      call vinokur(appox,nvalle+1,dx,dx1,dx2)
      do i=1,nvalle
        xs(i+nmonte+ncarlo) = appox(i+1) + x0
      end do

      dy = y1-y0
      call stre(appoy,nvalle+1,dy,0)
      do i=1,nvalle
        ys(i+nmonte+ncarlo) = appoy(i+1) + y0
      end do

c-----------------------------------------------------------------
c -- Lato est
c-----------------------------------------------------------------

      x0 =  xend
      x1 =  xend
      y0 =  yc(nc)
      y1 =  yend

      dx = x1-x0
      call stre(xe,ntra,dx,0)
      do i=1,ntra
        xe(i) = xe(i) + x0
      end do

      dy = y1-y0
      call stre(ye,ntra,dy,0)
      do i=1,ntra
        ye(i) = ye(i) + y0
      end do

c-----------------------------------------------------------------
c -- Lato ovest
c-----------------------------------------------------------------

      x0 =  xstart
      x1 =  xstart
      y0 =  yend
      y1 =  yc(1)

      dx = x1-x0
      call stre(xo,ntra,dx,0)
      do i=1,ntra
        xo(i) = xo(i) + x0
      end do

      dy = y1-y0
      call stre(yo,ntra,dy,0)
      do i=1,ntra
        yo(i) = yo(i) + y0
      end do

c-----------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine profiloi(ntra,nlo,ncarlo,nmonte,nvalle,
     &                    xc,yc,nc,
     &                    xstart,xend,ystart,yend,
     &                    xn,yn,xs,ys,xe,ye,xo,yo)

c-----------------------------------------------------------------

      real xn(1000),yn(1000),xs(1000),ys(1000)
      real xe(1000),ye(1000),xo(1000),yo(1000)
      real xc(1000),yc(1000)

      real appox(1000),appoy(1000),s(1000),snew(1000)

      real tens
      common/tensione/tens

c-----------------------------------------------------------------
c -- Profilo nord
c-----------------------------------------------------------------

c -- Prua

      x0 =  xstart
      x1 =  xc(1)
      y0 =  yc(1)
      y1 =  yc(1)

      dx = x1-x0
      dx1 = dx/float(nmonte)
      dx2 = (xc(nc)-xc(1))/float(ncarlo+1)/1.1
      call vinokur(xn,nmonte,dx,dx1,dx2)
      do i=1,nmonte
        xn(i) = xn(i) + x0
      end do

      dy = y1-y0
      call stre(yn,nmonte,dy,0)
      do i=1,nmonte
        yn(i) = yn(i) + y0
      end do

c -- Carena

      dx = xc(nc)-xc(1)
      dx1 = dx2
      dx2 = dx/float(ncarlo+1)/1.1
      call vinokur(appox,ncarlo+1,dx,dx1,dx2)
      do i=1,ncarlo+1
        appox(i) = appox(i) + xc(1)
      end do
      call interpola(xc,yc,nc,appox,appoy,ncarlo+1)

      do i=1,ncarlo
        xn(i+nmonte) = appox(i+1)
        yn(i+nmonte) = appoy(i+1)
      end do
      xn(ncarlo+nmonte) = xc(nc)
      yn(ncarlo+nmonte) = yc(nc)

c -- Poppa

      x0 =  xc(nc)
      x1 =  xend
      y0 =  yc(nc)
      y1 =  yc(nc)

      dx = x1-x0
      dx1 = dx2
      dx2 = dx/float(nvalle+1)
      call vinokur(appox,nvalle+1,dx,dx1,dx2)
      do i=1,nvalle
        xn(i+nmonte+ncarlo) = appox(i+1) + x0
      end do

      dy = y1-y0
      call stre(appoy,nvalle+1,dy,0)
      do i=1,nvalle
        yn(i+nmonte+ncarlo) = appoy(i+1) + y0
      end do

c-----------------------------------------------------------------
c -- Profilo sud
c-----------------------------------------------------------------

      do i=1,nlo
        xs(i) = xn(i)
        ys(i) = 0.
      end do

c-----------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------
      subroutine congiungi(x0,y0,x1,y1,xc,yc,np)
c-----------------------------------------------------------------
      real xc(1000),yc(1000)
c-----------------------------------------------------------------

      dx = x1-x0
      dy = y1-y0
      xsegno = 1.0
      ysegno = 1.0
      if(dx.lt.0) xsegno = -1.0
      if(dy.lt.0) ysegno = -1.0
      if(xsegno.eq.-1) dx = x0-x1
      if(ysegno.eq.-1) dy = y0-y1
      call stre(xc,np,dx,0)
      call stre(yc,np,dy,0)
      do i=1,np
        xc(i) = xsegno*xc(i) + x0
        yc(i) = ysegno*yc(i) + y0
      end do

c-----------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine skew(x,y,ni,nj)

c----------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj,nk

      real x(ki,kj),y(ki,kj)

c----------------------------------------------------------

      pig = acos(-1.0)

      do i=2,ni-1
        do j=1,nj-1

          dx = x(i+1,j+1)-x(i,j+1)
          dy = y(i+1,j+1)-y(i,j+1)

          alpha = -atan2(dy,dx)*180./pig

          if(alpha.gt.15) then
            alpha = 15.*pig/180.
            dy = dx*tan(alpha)
            y(i+1,j+1) = y(i,j+1)-dy
          end if

          if(x(i,j+1).lt.x(i,j)) x(i,j+1) = x(i,j)

        end do
      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine testar(x,y,ni,nj,armin)

c-----------------------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj

      real x(ki,kj),y(ki,kj)

c----------------------------------------------------------

      armin = 10.
      do i=1,ni-1
        do j=1,nj-1
          dx = x(i+1,j)-x(i,j)
          dy = y(i,j+1)-y(i,j)
          if(dx.ne.0) then
            ar = dy/dx
            if(ar.gt.0.) armin = min(ar,armin)
          end if
        end do
      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------

      subroutine bordo(x,y,xb,yb,nb,ni,nup,ndw,nj)

c-----------------------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj

      real x(ki,kj),y(ki,kj)
      real xb(ki),yb(ki)
      real xo(ki),yo(ki)
      real xi(ki),yi(ki)

c----------------------------------------------------------

      do j=1,nj
        y(1,j)  = y(2,j)
        y(ni-3,j) = y(ni-4,j)
        y(ni-2,j) = y(ni-3,j)
        y(ni-1,j) = y(ni-2,j)
        y(ni,j) = y(ni-1,j)
      end do

      y0 = y(1,1)
      y1 = y(1,nj)
      do j=1,nj
        y(1,j)  = y0 + (y1-y0)*(j-1.)/(nj-1.)
        y(ni,j) = y(1,j)
      end do

      do j=1,nj

c -- Monte

        x0 = x(1,j)
        y0 = y(1,j)
        x1 = x(nup,j)
        y1 = y(nup,j)
        do i=1,nup
          x(i,j) = x0 + (x1-x0)*(i-1.)/(nup-1.)
          y(i,j) = y0 + (y1-y0)*(i-1.)/(nup-1.)
        end do

c -- Carena

        x0 = x(nup,j)
        x1 = x(ni-ndw,j)
        k = 0
        do i=nup,ni-ndw
          k = k + 1
          x(i,j) = x0 + (x1-x0)*(k-1.)/(ni-nup-ndw)
        end do

c -- Coda

        x0 = x(ni-ndw,j)
        y0 = y(ni-ndw,j)
        x1 = x(ni,j)
        y1 = y(ni,j)
        k = 0
        do i=ni-ndw,ni
          k = k + 1
          x(i,j) = x0 + (x1-x0)*(k-1.)/float(ndw)
          y(i,j) = y0
        end do
      end do

      k = 0
      do i=nup,ni-ndw
        k = k + 1
        xo(k) = x(i,1)
      end do
      call interpola(xb,yb,nb,xo,yo,k)
      k = 0
      do i=nup,ni-ndw
        k = k + 1
        y(i,1) = yo(k)
      end do

      do j=1,nj
        y(1,j)  = y(2,j)
        y(ni-3,j) = y(ni-4,j)
        y(ni-2,j) = y(ni-3,j)
        y(ni-1,j) = y(ni-2,j)
        y(ni,j) = y(ni-1,j)
      end do


c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------

      subroutine esterno(x,y,ni,nj)

c-----------------------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj,nk

      real x(ki,kj),y(ki,kj)

c----------------------------------------------------------

      x0 = x(1,nj)
      x1 = x(ni,nj)
      do i=1,ni
        x(i,nj)   = x0 + (x1-x0)*(i-1.)/(ni-1.)
      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine thompson(x,y,ki,kj,imax,jmax)

      PARAMETER(NID=500,NJD=500,N=20) 
      real*4  X(ki,kj),Y(ki,kj) 
      COMMON/CONST/CHOICE,NA,DS1,DS2,NI,NJ
      COMMON/ATTR/IAL(N),IAX(N),IAY(N),JAL(N),JAX(N),JAY(N) 
      COMMON/COEF/AI(N),BI(N),CI(N),DI(N),AJ(N),BJ(N) 
      COMMON/COEF/CJ(N),DJ(N) 
      INTEGER CHOICE 
      DIMENSION P(NID,NJD),Q(NID,NJD),XX(0:NID,0:NJD),YY(0:NID,0:NJD) 
      DIMENSION X1(NID,NJD),X2(NID,NJD),Y1(NID,NJD),Y2(NID,NJD) 
C  
C     BOUNDARY INTERPOLATION  
C  
C      X          X ARRAY OF XI-ETA COORDINATE  
C      Y          Y ARRAY OF XI-ETA COORDINATE  
C      IMAX       MAX. NUMBER OF GRID IN XI AXIS  
C      JMAX       MAX. NUMBER OF GRID IN ETA AXIS  
C      NA         MAX. NUMBER OF ATTRACTIONS  
C      DS1        SPECIFIED LENGTH OF INITIAL INTERVAL  
C      DS2        SPECIFIED LENGTH OF FINAL INTERVAL   
C      ATTR       ARRAY OF ATTRACTION TO LINES/POINTS  
C      COEF       ARRAY OF COEFFICIENT FOR ATTRACTION  
C  
C      CHOICE    
C         1       VERTICAL INTERPOLATION     
C         2       HORIZONTAL INTERPOLATION   
C         3       TENSOR PRODUCT INTERPOLATION   
C         4       TRANSFINITE INTERPOLATION      
C         5       HERMITE CUBIC INTERPOLATION    
C         6       HYPERBOLIC TANGENT INTERPOLATION    
C         7       ELLIPTIC GRID GENERATION ( SOR ITERATION )       
C         8       ATTRACTION TO COORDINATES           
C             
       IF(CHOICE.EQ.1) GO TO 100  
       IF(CHOICE.EQ.2) GO TO 200   
       IF(CHOICE.EQ.3) GO TO 300                      
       IF(CHOICE.EQ.4) GO TO 400 
       IF(CHOICE.EQ.5) GO TO 500                            
       IF(CHOICE.EQ.6) GO TO 600  
       IF(CHOICE.EQ.7) GO TO 700                                  
       IF(CHOICE.EQ.8) GO TO 800  
C  
C      **** VERTICAL INTERPOLATION ****   
C             
 100   DO 110 I=1,IMAX        
       DO 110 J=2,JMAX-1
       RJ1=FLOAT(JMAX-J)/FLOAT(JMAX-1)  
       RJ2=FLOAT(J-1)/FLOAT(JMAX-1)     
C      *** ( EQ. 8-1 )  
       X(I,J)=RJ1*X(I,1)+RJ2*X(I,JMAX)            
 110   Y(I,J)=RJ1*Y(I,1)+RJ2*Y(I,JMAX)         
       RETURN                                     
C                                                  
C      **** HORIZONTAL INTERPOLATION ****             
C                                                  
 200   DO 210 I=1,JMAX                              
       DO 210 J=1,IMAX                            
       RI1=FLOAT(IMAX-I)/FLOAT(IMAX-1)            
       RI2=FLOAT(I-1)/FLOAT(IMAX-1)   
C      *** ( EQ. 8-1 )                                
       X(I,J)=RI1*X(1,J)+RI2*X(IMAX,J)            
 210   Y(I,J)=RI1*Y(1,J)+RI2*Y(IMAX,J)              
       RETURN                                     
C                                                  
C      **** TENSOR PRODUCT INTERPOLATION ****         
C                                                   
 300   DO 310 I=1,IMAX                              
       DO 310 J=1,JMAX                            
       RI1=FLOAT(IMAX-I)/FLOAT(IMAX-1)            
       RI2=FLOAT(I-1)/FLOAT(IMAX-1)               
       RJ1=FLOAT(JMAX-J)/FLOAT(JMAX-1)            
       RJ2=FLOAT(J-1)/FLOAT(JMAX-1)               
C      *** ( EQ. 8-69 )                               
       X(I,J)=RI1*RJ1*X(1,1)+RI1*RJ2*X(1,JMAX) 
     & +RI2*RJ1*X(IMAX,1)+RI2*RJ2*X(IMAX,JMAX) 
       Y(I,J)=RI1*RJ1*Y(1,1)+RI1*RJ2*Y(1,JMAX) 
     & +RI2*RJ1*Y(IMAX,1)+RI2*RJ2*Y(IMAX,JMAX) 
 310   CONTINUE                                      
       RETURN                                      
C 
C      **** TRANSFINITE INTERPOLATION ****             
C 
 400   DO 410 I=1,IMAX                               
       DO 410 J=1,JMAX                             
       RI1=FLOAT(I-1)/FLOAT(IMAX-1)                
       RI2=FLOAT(IMAX-I)/FLOAT(IMAX-1)             
       X1(I,J)=RI1*X(IMAX,J)+RI2*X(1,J)            
 410   Y1(I,J)=RI1*Y(IMAX,J)+RI2*Y(1,J)              
       DO 420 I=1,IMAX 
       DO 420 J=1,JMAX 
       RJ1=FLOAT(J-1)/FLOAT(JMAX-1) 
       RJ2=FLOAT(JMAX-J)/FLOAT(JMAX-1) 
       X2(I,J)=RJ1*(X(I,JMAX)-X1(I,JMAX))+RJ2*(X(I,1)-X1(I,1)) 
 420   Y2(I,J)=RJ1*(Y(I,JMAX)-Y1(I,JMAX))+RJ2*(Y(I,1)-Y1(I,1)) 
C      *** ( EQ. 8-73 ) 
       DO 430 I=1,IMAX                                         
       DO 430 J=1,JMAX                                         
       X(I,J)=X1(I,J)+X2(I,J)                                  
 430   Y(I,J)=Y1(I,J)+Y2(I,J)    
       IF(CHOICE.NE.4) GO TO 740                               
       RETURN                                                  
C                                                               
C      ***  HERMITE CUBIC INTERPOLATION (ORTHOGONAL BOUNDARY) *** 
C                                                               
 500   DO 510 I=1,IMAX                                          
       DO 510 J=1,JMAX                                     
       XX(I,J)=X(I,J)                                 
 510   YY(I,J)=Y(I,J)                       
       DO 520 J=1,JMAX                      
       XX(0,J)=XX(IMAX-1,J)            
       YY(0,J)=YY(IMAX-1,J)       
       XX(IMAX+1,J)=XX(2,J)  
 520   YY(IMAX+1,J)=YY(2,J)  
       DO 530 I=1,IMAX
C_DPE
       DS1 = X(I,2)-X(I,1)
       DS2 = X(I,JMAX)-X(I,JMAX-1)
       DO 530 J=1,JMAX
       RJJ=FLOAT(J-1)/FLOAT(JMAX-1)
C      *** ( EQ. 8-6 a and b, n=2 )   
       PHI1=(1.+2.*RJJ)*(1.-RJJ)*(1.-RJJ) 
       PHI2=(3.-2.*RJJ)*RJJ*RJJ           
       PSI1=(1.-RJJ)*(1.-RJJ)*RJJ         
       PSI2=(RJJ-1.)*RJJ*RJJ              
C  
C      ** CAL. NORMAL DERIV. **                
C  
       XXI1=.5*(XX(I+1,1)-XX(I-1,1))      
       XXI2=.5*(XX(I+1,JMAX)-XX(I-1,JMAX))
       YXI1=.5*(YY(I+1,1)-YY(I-1,1))      
       YXI2=.5*(YY(I+1,JMAX)-YY(I-1,JMAX))
       UNIT1=SQRT(XXI1*XXI1+YXI1*YXI1)    
       UNIT2=SQRT(XXI2*XXI2+YXI2*YXI2)  
C      *** ( EQ. 3-108 )                    
       XN1=-YXI1/UNIT1*DS1           
       XN2=-YXI2/UNIT2*DS2      
       YN1=XXI1/UNIT1*DS1  
       YN2=XXI2/UNIT2*DS2  
C      *** ( EQ. 8-5 )       
       XX(I,J)=PHI1*XX(I,1)+PHI2*XX(I,JMAX)+PSI1*XN1+PSI2*XN2 
 530   YY(I,J)=PHI1*YY(I,1)+PHI2*YY(I,JMAX)+PSI1*YN1+PSI2*YN2 
       DO 540 I=1,IMAX  
       DO 540 J=1,JMAX 
       X(I,J)=XX(I,J) 
 540   Y(I,J)=YY(I,J)   
       RETURN                                                
C 
C      **** HYPERBOLIC TANGENT SPACING INTERPOLATION **** 
C   
 600   TOL=1.0E-10    
C      *** ( EQ. 8-49, 50 and 51 
       A=SQRT(DS2/DS1)                                       
       B=1./(FLOAT(JMAX-1)*SQRT(DS1*DS2))                    
C      *** INITIAL GUESS BY SERIES EXPANSION                    
       DELTA=SQRT(6.*(B-1.))                
       DO 610 IT=1,20                    
       RESID=SINH(DELTA)/(DELTA*B)-1.         
       IF(ABS(RESID).LT.TOL) GO TO 630      
C_DPE 610   CALL AITKEN(DELTA,RESID,DELTO,RO,RSO)  
 610   CONTINUE
C_DPE       PRINT 620, RESID,DELTA,IT-1  
 620   FORMAT(//, 5X, 'DELTA IS NOT CONVERGE ?', 5X, 2E15.5, 
     & 5X, I3, //)                   
       GO TO 660                   
 630   CONTINUE                     
C      *** ( EQ. 8-52, 53 and 54 )    
       DO 650 I=1,IMAX             
       DO 650 J=2,JMAX-1           
       RATIO=FLOAT(J-1)/FLOAT(JMAX-1)  
       U=.5*(1.+TANH(DELTA*(RATIO-.5))/TANH(.5*DELTA)) 
       S=U/(A+(1.-A)*U)     
       X(I,J)=X(I,1)+(X(I,JMAX)-X(I,1))*S   
 650   Y(I,J)=X(I,1)+(Y(I,JMAX)-Y(I,1))*S    
 660   RETURN    
C   
C      **** ELLIPTIC GRID GENERATION ( SOR ITERATION ) ****  
C      *** CAL. P AND Q ON THE BOUNDARY **   
C   
 700   DO 710 I=1,IMAX    
       DO 710 J=1,JMAX   
       XXI=.5*(X(I+1,J)-X(I-1,J))       
       XXIXI=X(I+1,J)-2.*X(I,J)+X(I-1,J)  
       XETA=.5*(X(I,J+1)-X(I,J-1))    
       XETA2=X(I,J+1)-2.*X(I,J)+X(I,J-1)  

       YXI=.5*(Y(I+1,J)-Y(I-1,J))         
       YXIXI=Y(I+1,J)-2.*Y(I,J)+Y(I-1,J)  
       YETA=.5*(Y(I,J+1)-Y(I,J-1))        
       YETA2=Y(I,J+1)-2.*Y(I,J)+Y(I,J-1)  

       IF(ABS(XETA2).LT.1E-3) XETA2=0.   
       IF(ABS(YETA2).LT.1E-3) YETA2=0.   
       RXI2=XXI*XXI+YXI*YXI  
       RETA2=XETA*XETA+YETA*YETA                               
C      *** ( EQ. 8-70 )   
       P(I,J)=(XXI*XXIXI+YXI*YXIXI)/RXI2     
       Q(I,J)=(XETA*XETA2+YETA*YETA2)/RETA2  
 710   CONTINUE   
C  
C      **    INTERPOLATE P AND Q BETWEEN BOUNDARY **  
C            P : VERTICAL, Q : HORIZONTAL        
C  
       DO 720 I=1,IMAX       
       DO 720 J=1,JMAX          
       RJ1=FLOAT(JMAX-J)/FLOAT(JMAX-1)             
       RJ2=FLOAT(J-1)/FLOAT(JMAX-1)                   
       RI1=FLOAT(IMAX-I)/FLOAT(IMAX-1)                   
       RI2=FLOAT(I-1)/FLOAT(IMAX-1)                         
       P(I,J)=RJ1*P(I,1)+RJ2*P(I,JMAX)                         
       Q(I,J)=RI1*Q(1,J)+RI2*Q(IMAX,J)                         
 720   CONTINUE 
C  
C      **    INITIAL GUESS WITH TRANSFINITE INTERPOLATION **  
C      
       GO TO 400  
 740   CONTINUE  
C  
C      *** ITERATION ( SOR ) ***  
C  
       ITMAX=200 
       TOL=1.E-5                                              
       W=1.8                                                   
       DO 760 IT=1,ITMAX                                       
       ERRX=0.                                                 
       ERRY=0.                                                 
       DO 750 J=2,JMAX-1                                       
       DO 750 I=2,IMAX-1                                       
       XXI=.5*(X(I+1,J)-X(I-1,J))                              
       YXI=.5*(Y(I+1,J)-Y(I-1,J))                      
       XXIXI=X(I+1,J)+X(I-1,J)                                 
       YXIXI=Y(I+1,J)+Y(I-1,J)                                 
       XETA=.5*(X(I,J+1)-X(I,J-1))                             
       YETA=.5*(Y(I,J+1)-Y(I,J-1))                             
       XXIETA=.25*(X(I+1,J+1)-X(I+1,J-1)-X(I-1,J+1)+X(I-1,J-1))
       YXIETA=.25*(Y(I+1,J+1)-Y(I+1,J-1)-Y(I-1,J+1)+Y(I-1,J-1))
       XETA2=X(I,J+1)+X(I,J-1)                                  
       YETA2=Y(I,J+1)+Y(I,J-1)                                  
C      *** ( EQ. 6-18 and 6-20 )    
       G11=XXI*XXI+YXI*YXI                                      
       G22=XETA*XETA+YETA*YETA 
       G12=XXI*XETA+YXI*YETA 
       XTEMP=.5*(G22*(P(I,J)*XXI+XXIXI)+G11*(Q(I,J)*XETA+XETA2) 
     & -2.*G12*XXIETA)/(G11+G22)                                     
       YTEMP=.5*(G22*(P(I,J)*YXI+YXIXI)+G11*(Q(I,J)*YETA+YETA2) 
     & -2.*G12*YXIETA)/(G11+G22)                                    
       XTEMP=W*XTEMP+(1.-W)*X(I,J)                                  
       YTEMP=W*YTEMP+(1.-W)*Y(I,J)                              
       ERRX=MAX(ERRX,ABS(XTEMP-X(I,J)))                           
       ERRY=MAX(ERRY,ABS(YTEMP-Y(I,J)))                           
       X(I,J)=XTEMP                                                 
       Y(I,J)=YTEMP                                                 
 750   CONTINUE            
       IF(ERRX.LT.TOL.AND.ERRY.LT.TOL) GOTO 780
 760   CONTINUE                 
C_DPE       PRINT 770,ERRX,ERRY,IT-1                                     
 770   FORMAT(//, 5X, 'X AND Y ARE NOT CONVERGE ?', 2E15.5,         
     & 5X, I5, //)                                                  
 780   CONTINUE 
       IF(CHOICE.EQ.8) GO TO 830                                    
       RETURN   
C 
C      **** ATTRACTION TO COORDINATE LINE/POINT **** 
C 
 800   DO 810 I=1,IMAX  
       DO 810 J=1,JMAX                                              
       P(I,J)=0.                                            
 810   Q(I,J)=0.    
       DO 820 NS=1,NA                                      
       DO 820 I=1,IMAX                                  
       DO 820 J=1,JMAX                               

       XL=FLOAT(I-IAL(NS))
       XI=FLOAT(I-IAX(NS))
       XJ=FLOAT(J-IAY(NS))

       YL=FLOAT(J-JAL(NS))
       YI=FLOAT(I-JAX(NS))
       YJ=FLOAT(J-JAY(NS))

C      *** ( EQ. 6-30 )
       XDUM = XI*XI + XJ*XJ
       IF (XDUM.NE.0) XDUM = SQRT(XDUM)
       YDUM = YI*YI + YJ*YJ
       IF (YDUM.NE.0) YDUM = SQRT(YDUM)

C -- AI: Peso della parte dipendente dalla distanza da XL
C -- BI: Peso della parte dipendente dalla distanza da XI
C -- CI: Decadimento dell'attrazione al variare di XL
C -- DI: Decadimento dell'attrazione al variare di XI

       F1 = 0
       IF(XL.NE.0) F1 = AI(NS)*(XL/ABS(XL))
       F2 = EXP(-CI(NS)*ABS(XL))  
       F3 = 0
       IF(XI.NE.0) F3 = BI(NS)*(XI/ABS(XI))
       F4 = EXP(-DI(NS)*XDUM)  
       DP=F1*F2+F3*F4
       P(I,J)=P(I,J)-DP

C       P(I,J)=P(I,J)-AI(NS)*(XL/ABS(XL))*EXP(-CI(NS)*ABS(XL))  
C     &              -BI(NS)*(XI/ABS(XI))*EXP(-DI(NS)*XDUM)     

C -- AJ: Peso della parte dipendente dalla distanza da YL
C -- BJ: Peso della parte dipendente dalla distanza da YI
C -- CJ: Decadimento dell'attrazione al variare di YL
C -- DJ: Decadimento dell'attrazione al variare di YI

       F1 = 0
       IF(YL.NE.0) F1 = AJ(NS)*(YL/ABS(YL))
       F2 = EXP(-CJ(NS)*ABS(YL))  
       F3 = 0
       IF(YJ.NE.0) F3 = BJ(NS)*(YJ/ABS(YJ))
       F4 = EXP(-DJ(NS)*YDUM)  
       DQ=F1*F2+F3*F4
       Q(I,J)=Q(I,J)-DQ

C       Q(I,J)=Q(I,J)-AJ(NS)*(YL/ABS(YL))*EXP(-CJ(NS)*ABS(YL))  
C     &              -BJ(NS)*(YJ/ABS(YJ))*EXP(-DJ(NS)*YDUM)     

 820   CONTINUE  
       GO TO 400 
 830   CONTINUE  
       RETURN  
       END
c ----------------------------------------------------------------------
c     SUBROUTINES
c-----------------------------------------------------------------------

C________________________________________________________________
      subroutine ordina2(x,y,n)
c
c     Ordina due vettori in base al primo
c
      integer n,i,j
      real x(*),y(*)
      real appo
c
      do i=1,n-1
        do j=i+1,n
          if(x(j).lt.x(i)) then
            appo = x(i)
            x(i) = x(j)
            x(j) = appo
            appo = y(i)
            y(i) = y(j)
            y(j) = appo
          end if
        end do
      end do
c
c
      return
      end
c----------------------------------------------------------

      subroutine acurv2(x,y,nl,s)
c
      integer i,nl
      real x(*),y(*),s(*)
      real dx,dy,dz
c
      s(1)  = 0.
c
      do i=2,nl
        dx =  x(i)-x(i-1)
        dy =  y(i)-y(i-1)
        s(i) = s(i-1)+sqrt(dx**2+dy**2)
      end do
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine stre(z,npan,l,istre)
c-----------------------------------------------------------------------
C  Restituisce in z() npan valori spaziati tra 0 ed l
C  secondo l'indice 'istre'
C  istre=0  Equispaziati
C  istre=1  Coseno: addensati agli estremi
C  istre=2  Seno:   addensati verso l'estremo superiore
C  istre=3  Seno:   addensati verso lo zero
c-----------------------------------------------------------------------
      real z(*)
      integer i,npan,istre
      real l,dz,pig
c-----------------------------------------------------------------------
C Genera un vettore, con valori tra 0 ed l
c-----------------------------------------------------------------------
      npan = npan-1
      pig = 4.0*atan(1.0)
C
      if (istre.eq.0) then
        dz = 1./npan
        do i=1,npan+1
          z(i) = l*dz*(i-1)
        end do
      end if
      if (istre.eq.1) then
        dz = pig/npan
        do i=1,npan+1
          z(i) = l*(1-cos(dz*(i-1)))/2.
        end do
      end if
      if (istre.eq.2) then
        dz = pig/(2.*npan)
        do i=1,npan+1
          z(i) = l*(sin(dz*(i-1)))
        end do
      end if
      if (istre.eq.3) then
        dz = pig/(2.*npan)
        do i=1,npan+1
          z(i) = l*(1.-sin(dz*(npan+1-i)))
        end do
      end if
C
      npan = npan+1

      return
      end
c-----------------------------------------------------------------

      subroutine strex(z,npan,elle,istre,peso)

c-----------------------------------------------------------------
c
c  Restituisce in z() npan valori spaziati tra 0 ed l
c  secondo l'indice 'istre'
c  istre=0  Equispaziati
c  istre=1  Coseno: addensati agli estremi
c  istre=2  Seno:   addensati verso l'estremo superiore
c  istre=3  Seno:   addensati verso lo zero
c-----------------------------------------------------------------
      integer i,npan,istre
      real z(*),z1(10000),z2(10000)
      real elle,dz,pig,peso
c-----------------------------------------------------------------
c
c Genera un vettore, con valori tra 0 ed l
c
c-----------------------------------------------------------------

      pig = acos(-1.0)

      dz = 1.0/float(npan-1)
      do i=1,npan
        z1(i) = elle*dz*float(i-1)
        z2(i) = elle*dz*float(i-1)
      end do
      if (istre.eq.1) then
        dz = pig/float(npan-1)
        do i=1,npan
          z2(i) = elle*(1.0-cos(dz*float(i-1)))/2.0
        end do
      end if
      if (istre.eq.2) then
        dz = pig/(2.0*float(npan)-1.0)
        do i=1,npan
          z2(i) = elle*(sin(dz*float(i-1)))
        end do
      end if
      if (istre.eq.3) then
        dz = pig/2.0/float(npan-1)
        do i=1,npan
          z2(i) = elle*(1.0-sin(dz*float(npan-i)))
        end do
      end if

      do i=1,npan
        z(i) = peso*z1(i)+(1.0-peso)*z2(i)
      end do

c-----------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------
C
      subroutine interpolas(x,y,np,xi,yi,ni,tens)
C________________________________________________________________
C
C  Interpola una curva.
C  In ingresso: ascisse ed ordinate, numero di punti.
C  Le ascisse devono essere in senso crescente.
C________________________________________________________________
C
      real x(*),y(*),temp(5000),yp(5000),d1,d2,tens
      real xi(*),yi(*)
      integer np,ni,ierr,ispl
C
      d1 = 0.
      d2 = 0.
      ispl = 3
      itest = 0
      do i=2,np
        if(y(i).ne.y(1)) itest = 1
      end do
      if(itest.eq.0) then
        do i=1,ni
          yi(i)=y(1)
        end do
        return
      end if
      call curv1 (np,x,y,d1,d2,ispl,yp,temp,tens,ierr)
      do i=1,ni
        yi(i) = curv2(xi(i),np,x,y,yp,tens)
      end do
C________________________________________________________________
      return
      end

c-----------------------------------------------------------------------

      subroutine buildtrans(x0e,nde,ntrans)

c-----------------------------------------------------------------------

      real x0e(*),xstern(1000),ystern(1000),zstern(1000)
      real xdot(1),zdot(1)

c-----------------------------------------------------------------------

      iuno = 1
      ierr = 1
      open(20,file='stern.profile',
     &        status='old',
     &        form='formatted',
     &        err=123)
      close(20)
      ierr = 0

      open(20,file='stern.profile',
     &        status='old',
     &        form='formatted')
      read(20,*) nstern
      do i=1,nstern
        read(20,*) xstern(i),ystern(i),zstern(i)
      end do
      close(20)
      iuno = 1
      xdot(1) = x0e(nde)-0.01
      call interpola(xstern,zstern,nstern,xdot,zdot,iuno)
      zp = zdot(1)
      call interpola(xstern,ystern,nstern,xdot,zdot,iuno)
      yp = zdot(1)
      xdot(1) = x0e(nde)-0.02
      call interpola(xstern,zstern,nstern,xdot,zdot,iuno)
      zm = zdot(1)
      dzdx = (zp-zm)/0.01
      dzdx = atan(dzdx)

 123  if(ierr.eq.0) then
        open(44,file='transom.grd',status='unknown')
        write(44,*) x0e(nde),yp
        write(44,*) -1000.,zp
        write(44,*)  1000.,zp
        write(44,*)  dzdx
        write(44,*)  dzdx
        close(44)
      else
        ntrans = 0
      end if

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
      subroutine vinokur(s,lmax,smax,ds1e,ds2e)
c--------------------------------------------------------------------
c
c  stretches points on a surface so that a specified spacing
c  at the boundaries is satisfied. Taken from NASA CR 3313 by
c  Vinokur (1980).
c
c  In this version, 4 distinct iterations are made to better
c  match the resulting delta-s values to the requested values.
c  The four iterations are summarized below:
c
c  1.  delta-s is set equal to the desired value.
c  2.  delta-s from the last iteration is corrected from a Taylor
c      series expansion.
c  3.  delta-s is calculated from a linear fit between the first two
c      guesses.
c  4.  delta-s is calculated from a quadratic fit between the first
c      three guesses, if indeed a quadratic will pass through the
c      desired value.  If it doesn't, it takes the value calculated
c      after three swipes.
c
c  Additionally, this version uses the approximate inverse solution
c  for y=sin(x)/x and y=sinh(x)/x rather than a Newton iteration.  The
c  approximate solution was also taken from NASA CR 3313.
c
c
      dimension s(1000), d1(4,2),d2(4,2)
c
c---------------------------------
c  for an IRIS 2500,
c--------------------------------
c
      dsavg=smax/float(lmax-1)
      if(ds1e.eq.0.0.and.ds2e.eq.0.0)then
	 kase=0
	 ds1e=dsavg
	 ds2e=dsavg
	 nlast=4
      else if(ds1e.eq.0.0)then
	 kase=1
	 nlast=1
 23      write(*,106)
c        call realval(0,1,slop,no,no,*23,*101)
	 read(*,*) slop
	 if(slop.lt.0.0.or.slop.gt.1.0)go to 23
	 ds1e=-slop
      else if(ds2e.eq.0.0)then
	 kase=2
	 nlast=1
 24      write(*,106)
	 read(*,*) slop
c        call realval(0,1,slop,no,no,*24,*101)
	 if(slop.lt.0.0.or.slop.gt.1.0)go to 24
	 ds2e=-slop
      else
	 kase=0
	 nlast=4
      end if
      dss1=0.0
      dss2=0.0
c
      do 6 n=1,nlast-2
      if(n.le.2)then
	 ds1=ds1e-0.5*dss1
	 ds2=ds2e+0.5*dss2
	 d1(n,1)=ds1
	 d2(n,1)=ds2
      else if(n.eq.3)then
	 ds1=-d1(1,2)*(d1(2,1)-d1(1,1))/(d1(2,2)-d1(1,2))+d1(1,1)
	 ds2=-d2(1,2)*(d2(2,1)-d2(1,1))/(d2(2,2)-d2(1,2))+d2(1,1)
	 if(ds1.lt.0.0)ds1=0.5*min(d1(1,1),d1(2,1))
	 if(ds2.lt.0.0)ds2=0.5*min(d2(1,1),d2(2,1))
	 d1(n,1)=ds1
	 d2(n,1)=ds2
      else if(n.eq.4)then
	 denom=-(d1(1,1)-d1(2,1))*(d1(2,1)-d1(3,1))*(d1(3,1)-d1(1,1))
	 a11=d1(2,1)-d1(3,1)
	 a21=d1(3,1)**2-d1(2,1)**2
	 a31=d1(2,1)*d1(3,1)*(d1(2,1)-d1(3,1))
	 a12=d1(3,1)-d1(1,1)
	 a22=d1(1,1)**2-d1(3,1)**2
	 a32=d1(3,1)*d1(1,1)*(d1(3,1)-d1(1,1))
	 a13=d1(1,1)-d1(2,1)
	 a23=d1(2,1)**2-d1(1,1)**2
	 a33=d1(1,1)*d1(2,1)*(d1(1,1)-d1(2,1))
	 b1=(a11*d1(1,2)+a12*d1(2,2)+a13*d1(3,2))/denom
	 b2=(a21*d1(1,2)+a22*d1(2,2)+a23*d1(3,2))/denom
	 b3=(a31*d1(1,2)+a32*d1(2,2)+a33*d1(3,2))/denom
	 disc=(b2*b2-4.*b1*b3)
	 if(disc.lt.0.0)go to 8
	 dd1=(-b2+sqrt(disc))/(2.*b1)
	 dd2=(-b2-sqrt(disc))/(2.*b1)
	 dd3=d1(3,1)
	 if(abs(dd1-dd3).lt.abs(dd2-dd3))then
	    ds1=dd1
	 else
	    ds1=dd2
	 end if
	 denom=-(d2(1,1)-d2(2,1))*(d2(2,1)-d2(3,1))*(d2(3,1)-d2(1,1))
	 a11=d2(2,1)-d2(3,1)
	 a21=d2(3,1)**2-d2(2,1)**2
	 a31=d2(2,1)*d2(3,1)*(d2(2,1)-d2(3,1))
	 a12=d2(3,1)-d2(1,1)
	 a22=d2(1,1)**2-d2(3,1)**2
	 a32=d2(3,1)*d2(1,1)*(d2(3,1)-d2(1,1))
	 a13=d2(1,1)-d2(2,1)
	 a23=d2(2,1)**2-d2(1,1)**2
	 a33=d2(1,1)*d2(2,1)*(d2(1,1)-d2(2,1))
	 b1=(a11*d2(1,2)+a12*d2(2,2)+a13*d2(3,2))/denom
	 b2=(a21*d2(1,2)+a22*d2(2,2)+a23*d2(3,2))/denom
	 b3=(a31*d2(1,2)+a32*d2(2,2)+a33*d2(3,2))/denom
	 disc=(b2*b2-4.*b1*b3)
	 if(disc.le.0.0)go to 8
	 dd1=(-b2+sqrt(disc))/(2.*b1)
	 dd2=(-b2-sqrt(disc))/(2.*b1)
	 dd3=d2(3,1)
	 if(abs(dd1-dd3).lt.abs(dd2-dd3))then
	    ds2=dd1
	 else
	    ds2=dd2
	 end if
	 if(ds1.lt.0.0.or.ds2.lt.0.0)go to 8
      end if
c
c  calculate constants
      s0=smax/float(lmax-1)/ds1
      s1=smax/float(lmax-1)/ds2
      b=sqrt(s0*s1)
      a=sqrt(s0/s1)
      if(kase.eq.1)then
	 b=s1
      else if(kase.eq.2)then
	 b=s0
      end if
c
c  calculate x based on value of B
      if(b-1.)1,2,3
c
c  x is real
 1    if(b.lt.0.26938972)then
	 pi=4.*atan(1.)
	 x=pi*(1.   -b     +          b**2 - (1.+pi**2/6.)*b**3
     *    + 6.794732*b**4  -13.205501*b**5     + 11.726095*b**6)
      else
	 c=1.-b
	 x= sqrt(6.*c)*(1.
     *           +0.15*c    + 0.057321429*c**2 +0.048774238*c**3
     *    -0.053337753*c**4 + 0.075845134*c**5)
      end if
      go to 4
c
c  x is zero
 2    x=0.
      go to 4
c
c  x is imaginary
 3    if(b.lt.2.7829681)then
	 c=b-1.
	 x= sqrt(6.*c)*(1.
     *           -0.15*c    + 0.057321429*c**2 -0.024907295*c**3
     *   +0.0077424461*c**4 -0.0010794123*c**5)
      else
	 v=log(b)
	 w=1./b - 0.028527431
	 x= v + (1.+1./v)*log(2.*v) -0.02041793
     *    + 0.24902722*w     +  1.9496443*w**2
     *    -  2.6294547*w**3  + 8.5679591*w**4
      end if
c
c  distribute points along boundary
c
 4    continue
      if(kase.eq.1.or.kase.eq.2)then
	 s(1   ) = 0.0
	 s(lmax) = smax
	 do 9 i=2,lmax-1
	      j= lmax+1-i
	 xi=float(i-1)/(lmax-1)
	 if(b.gt.1.0001)then
	    u1=1. + tanh(x/2.*(xi-1.))/tanh(x/2.)
	 else if(b.lt.0.9999)then
	    u1=1. + tan (x/2.*(xi-1.))/tan (x/2.)
	 else
	    u1= xi*(1.-.5*(b-1.)*(1.-xi)*(2.-xi))
	 end if
	 u2=sinh(xi*x)/sinh(x)
	 if(kase.eq.1)then
	    fact=abs(ds1e)
	    s(j) = ( (1.-fact)*(1.-u1) + fact*(1.-u2) ) *smax
	 else if(kase.eq.2)then
	    fact=abs(ds2e)
	    s(i) = ( (1.-fact)*    u1  + fact*    u2  ) *smax
	 end if
 9       continue
      else
	 do 5 i=1,lmax
	 xi=float(i-1)/float(lmax-1)
	 cnum=x*(xi-0.5)
	 cden=x/2.
	 if(b.lt.0.9999)then
	    cc=tan(cnum)/tan(cden)
	    u=0.5*(1.+cc)
	 else if(b.ge.0.9999.and.b.le.1.0001)then
	    u=xi*(1.+2.*(b-1.)*(xi-0.5)*(1.-xi))
	 else if(b.gt.1.0001)then
	    cc=tanh(cnum)/tanh(cden)
	    u=0.5*(1.+cc)
	 end if
 5       s(i)=u*smax/(a+(1.-a)*u)
      end if
c
      if(lmax.ge.4)then
        dss1=(  -s(4)   +4.*s(3)     -5.*s(2)     +2.*s(1))
        dss2=(2.*s(lmax)-5.*s(lmax-1)+4.*s(lmax-2)   -s(lmax-3))/2.
      end if
c
      es1=s(2)-s(1)
      es2=s(lmax)-s(lmax-1)
      if(n.ne.4)then
        d1(n,2)=es1-ds1e
        d2(n,2)=es2-ds2e
      end if
 6    continue
c
 8    esmin= 1.0e+08
      esmax=-1.0e+08
      do 7 j=2,lmax
      stmp=s(j)-s(j-1)
      if(stmp.lt.esmin)then
        esmin=stmp
      end if
      if(stmp.gt.esmax)then
        esmax=stmp
      end if
 7    continue
c
 101  return
 103  format(/,6x,'enter delta s at beginning of arclength',
     *  /,6x,'(default = ',g12.5,' 0.= auto-spacing)',t59,'>',$)
 104  format(/,6x,'enter delta s at end of arclength',
     *  /,6x,'(default = ',g12.5,' 0.= auto-spacing)',t59,'>',$)
 105  format(/,6x,'computed spacing at beginning:',g12.5,/,
     *6x,         '                          end:',g12.5,/,
     *6x,         'minimum spacing   (i=',i3,',',i3,'):',g12.5,/,
     *6x,         'maximum spacing   (i=',i3,',',i3,'):',g12.5)
 106  format(6x,'   enter the degree of stretching',/,
     *       6x,'   (between 0. (tanh) and 1. (sinh) )',t59,'>',$)
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      end
c-----------------------------------------------------------------------

      subroutine tm(x,y,ni,nj)

c----------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj,nk

      real x(ki,kj),y(ki,kj)
      real xnew(ki,kj),ynew(ki,kj)
      real so(ki),xo(ki),yo(ki)
      real si(ki),xi(ki),yi(ki)

c----------------------------------------------------------

      domega  = 0.8
      domegaB = 0.8
      ncyc   = max(ni,nj)*4

      do icyc =1,ncyc

        do i=2,ni-1
          do j=2,nj-1

            ds1 = sqrt( (x(i+1,j)-x(i,j))**2 + 
     &                  (y(i+1,j)-y(i,j))**2 )

            ds2 = sqrt( (x(i-1,j)-x(i,j))**2 + 
     &                  (y(i-1,j)-y(i,j))**2 )

            ds3 = sqrt( (x(i,j+1)-x(i,j))**2 + 
     &                  (y(i,j+1)-y(i,j))**2 )

            ds4 = sqrt( (x(i,j-1)-x(i,j))**2 + 
     &                  (y(i,j-1)-y(i,j))**2 )

            sums = ds1 + ds2 + ds3 + ds4

            dmu1 = 1.0-ds1/sums
            dmu2 = 1.0-ds2/sums
            dmu3 = 1.0-ds3/sums
            dmu4 = 1.0-ds4/sums

            xnew(i,j) = (dmu1*x(i+1,j) +
     &                   dmu2*x(i-1,j) +
     &                   dmu3*x(i,j+1) +
     &                   dmu4*x(i,j-1) +
     &                    4.0*x(i,j)   )/7.

            ynew(i,j) = (dmu1*y(i+1,j) +
     &                   dmu2*y(i-1,j) +
     &                   dmu3*y(i,j+1) +
     &                   dmu4*y(i,j-1) +
     &                    4.0*y(i,j)   )/7.

          end do
        end do

        domega = domega*1.01
        if(domega.gt.domegaB) domega = domegaB

        do i=2,ni-1
          do j=2,nj-1

            x(i,j) = xnew(i,j)*(1.-domega) + x(i,j)*domega
            y(i,j) = ynew(i,j)*(1.-domega) + y(i,j)*domega

          end do
        end do

        if(icyc.le.ncyc/2) then 
          do j=1,nj/2
            do i=2,ni-1
              dx1 = x(i  ,j)-x(i-1,j)
              dy1 = y(i  ,j)-y(i-1,j)
              t1  = dy1/dx1
              dx2 = x(i+1,j)-x(i  ,j)
              dy2 = y(i+1,j)-y(i  ,j)
              t2  = dy2/dx2
              tl  = 0.5*(t1+t2)
              dy1 = y(i,j+1)-y(i,j)
              dx1 = x(i,j+1)-x(i,j)
              tt  = dx1/dy1
              dt  = tt - tl
              if(dt.gt.0) then
                dt =  1.
              else if(dt.lt.0) then
                dt = -1.
              else
                dt =  0.
              end if
              x(i,j+1)  = x(i,j+1) + 1.E-4*dt
            end do
          end do
        end if

        do j=1,nj
          y(1 ,j) = y(   2,j)
          y(ni,j) = y(ni-1,j)
        end do

      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine laplace(x,y,ni,nj)

c----------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj,nk

      real x(ki,kj),y(ki,kj)
      real xnew(ki,kj),ynew(ki,kj)
      real v(4,2)

c----------------------------------------------------------

      domegaB = 0.5
      domega  = 0.1
      ncyc   = max(ni,nj)*4

      do icyc =1,ncyc

        do i=2,ni-1
          do j=2,nj-1

            v(1,1) = x(i+1,j)-x(i,j)
            v(1,2) = y(i+1,j)-y(i,j)

            v(2,1) = x(i-1,j)-x(i,j)
            v(2,2) = y(i-1,j)-y(i,j)

            v(3,1) = x(i,j+1)-x(i,j)
            v(3,2) = y(i,j+1)-y(i,j)

            v(4,1) = x(i,j-1)-x(i,j)
            v(4,2) = y(i,j-1)-y(i,j)

            dx = 0.
            dy = 0.
            do k=1,4
              dx = dx + v(k,1)
              dy = dy + v(k,1)
            end do
            dx = dx/4.
            dy = dy/4.

            xnew(i,j) = x(i,j) + dx
            ynew(i,j) = y(i,j) + dy

          end do
        end do

        domega = domega*1.01
        if(domega.gt.domegaB) domega = domegaB

        do i=2,ni-1
          do j=2,nj-1

            x(i,j) = xnew(i,j)*(1.-domega) + x(i,j)*domega
            y(i,j) = ynew(i,j)*(1.-domega) + y(i,j)*domega

          end do
        end do

        do j=1,nj
          y(1 ,j) = y(   2,j)
          y(ni,j) = y(ni-1,j)
        end do

      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine angles(x,y,ni,nj,ilast)

c----------------------------------------------------------

      parameter (ki=500,kj=500)
      integer ni,nj,nk

      integer iv(4,4),jv(4,4)
      real x(ki,kj),y(ki,kj)
      real xnew(ki,kj),ynew(ki,kj)
      real xp(ki),xi(ki)
      real sp(ki),si(ki)
      real v(3,2),dmod(3)

c----------------------------------------------------------

      domegaB = 0.95
      domega  = 0.95
      ncyc    = min(ni,nj)
      ncyc    = 20

      si(1) = 0.
      do j=2,nj
        dx = x(ni/2,j)-x(ni/2,j-1)
        dy = y(ni/2,j)-y(ni/2,j-1)
        ds = sqrt(dx*dx + dy*dy)
        si(j) = si(j-1) + ds
      end do
      do j=2,nj
        si(j) = si(j)/si(nj)
      end do

      do icyc =1,ncyc

        domega = domega*1.25
        if(domega.gt.domegaB) domega = domegaB

        do i=2,ni-1
          do j=2,nj-1

            call vertici(i,j,iv,jv)

            dx = 0.
            dy = 0.

            xb =(x(i-1,j+1)+
     &           x(i+1,j+1)+
     &           x(i+1,j-1)+
     &           x(i-1,j-1))/4.
            yb =(y(i-1,j+1)+
     &           y(i+1,j+1)+
     &           y(i+1,j-1)+
     &           y(i-1,j-1))/4.

            do 100 kv=1,4

              do kk = 1,3
                v(kk,1) = x(iv(kk,kv),jv(kk,kv))-x(iv(4,kv),jv(4,kv))
                v(kk,2) = y(iv(kk,kv),jv(kk,kv))-y(iv(4,kv),jv(4,kv))
              end do

              do kk=1,3
                dmod(kk) = v(kk,1)*v(kk,1)+v(kk,2)*v(kk,2)
                if(dmod(kk).ne.0) then
                  dmod(kk) = sqrt(dmod(kk))
                else
                  go to 100
                end if
              end do

              do kk = 1,3
                v(kk,1) = v(kk,1)/dmod(kk)
                v(kk,2) = v(kk,2)/dmod(kk)
              end do

              a1 = v(1,1)*v(2,1)+v(1,2)*v(2,2)
              a2 = v(2,1)*v(3,1)+v(2,2)*v(3,2)
              if(a1.gt.1)  go to 100
              if(a2.gt.1)  go to 100
              if(a1.lt.-1) go to 100
              if(a2.lt.-1) go to 100

              a1 = acos(a1)
              a2 = acos(a2)

              beta  = 0.5*(a1-a2)
              erre  = dmod(2)*tan(beta)
              dx = dx - erre*cos(beta)
              dy = dy + erre*sin(beta)

  100       continue

c -- Somma

            xnew(i,j) = x(i,j) + dx/4.
            ynew(i,j) = y(i,j) + dy/4.

            xnew(i,j) = domega*xnew(i,j) + (1.-domega)*xb
            ynew(i,j) = domega*ynew(i,j) + (1.-domega)*yb

          end do
        end do

        do i=2,ni-1
          do j=2,nj-1

            x(i,j) = xnew(i,j)*(1.-domega) + x(i,j)*domega
            y(i,j) = ynew(i,j)*(1.-domega) + y(i,j)*domega

          end do
        end do

        do i=2,ni-1

          sp(1) = 0.
          do j=2,nj
            dx = x(1,j)-x(1,j-1)
            dy = y(1,j)-y(1,j-1)
            ds = sqrt(dx*dx + dy*dy)
            sp(j) = sp(j-1) + ds
          end do
          do j=2,nj
            sp(j) = sp(j)/sp(nj)
          end do

          do j=1,nj
            xp(j) = x(i,j)
          end do
          call interpola(sp,xp,nj,si,xi,nj)
          do j=1,nj
            x(i,j) = xi(j)
          end do

          do j=1,nj
            xp(j) = y(i,j)
          end do
          call interpola(sp,xp,nj,si,xi,nj)
          do j=1,nj
            y(i,j) = xi(j)
          end do

        end do

        if(ilast.eq.0) njlast = nj-1
        if(ilast.eq.1) njlast = nj
        do i=1,ni
          do j=2,ilast
            x(i,j) = max(x(i,j-1),x(i,j))
          end do
        end do

        do j=1,nj
          y(1 ,j) = y(   2,j)
          y(ni,j) = y(ni-1,j)
        end do

      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------

      subroutine vertici(i,j,iv,jv)

c-----------------------------------------------------------------

      integer iv(4,4),jv(4,4)

c-----------------------------------------------------------------

c -- Primo vertice

      iv(1,1) = i-1
      jv(1,1) = j

      iv(2,1) = i
      jv(2,1) = j

      iv(3,1) = i
      jv(3,1) = j+1

      iv(4,1) = i-1
      jv(4,1) = j+1

c -- Secondo vertice

      iv(1,2) = i
      jv(1,2) = j+1

      iv(2,2) = i
      jv(2,2) = j

      iv(3,2) = i+1
      jv(3,2) = j

      iv(4,2) = i+1
      jv(4,2) = j+1

c -- Terzo vertice

      iv(1,3) = i+1
      jv(1,3) = j

      iv(2,3) = i
      jv(2,3) = j

      iv(3,3) = i
      jv(3,3) = j-1

      iv(4,3) = i+1
      jv(4,3) = j-1

c -- Quarto vertice

      iv(1,4) = i
      jv(1,4) = j-1

      iv(2,4) = i
      jv(2,4) = j

      iv(3,4) = i-1
      jv(3,4) = j

      iv(4,4) = i-1
      jv(4,4) = j-1

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------
