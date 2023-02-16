c ----------------------------------------------------------------------

      subroutine tracut(nlast,nlo,nlovallec,ntraest,ntraint,ntrapoc,
     &                  ncar,nplib,fr,sup,eta,xns,yns)

c ----------------------------------------------------------------------

      include 'warp.cmn'

      integer nlast,nlo,nlovallec,ntraest,ntraint,ntrapoc,ncar,nplib
      real fr,sup
      real xns(ntotd),yns(ntotd),eta(nplibd)
      real xfs(400,400),yfs(400,400),zfs(400,400),zx(400,400)
      real sp(400),xp(400),yp(400),zp(400)
      real si(400),xi(400),yi(400),zi(400)

      real d(2),cint(4000),vint(5),wint(4000)

      character*80 fina
      character*4  cha

c ----------------------------------------------------------------------

      newtra = 256

      !write(*,*)'----------========== T R A C U T ==========----------'

c -- Riorganizza i dati

      i0 = ncar
      ymax = 0.

      ifr = nint(fr*1000)
      write(cha,99) ifr
      fina='tatra'//cha//'.raw'
      open(22,file=fina)
      do i=ncar+1,ncar+nplib
        write(22,*) xns(i),yns(i),eta(i-ncar)
      end do
      close(22)

c -- Dominio "normale"

      do i=1,nlo
        do j=1,ntraest+ntraint
          i0 = i0 + 1
          xfs(i,j) = xns(i0)
          yfs(i,j) = yns(i0)
          zfs(i,j) = eta(i0-ncar)
        end do
      end do

c -- Transom

      do i=1,nlo-nlovallec
        do j=ntraest+ntraint+1,ntraest+ntraint+ntrapoc
          xfs(i,j) = xfs(i,1)
          yfs(i,j) = yfs(i,1)
          zfs(i,j) = zfs(i,1)
        end do
      end do

      do i=nlo-nlovallec+1,nlo
        do j=ntraest+ntraint+1,ntraest+ntraint+ntrapoc
          i0 = i0 + 1
          xfs(i,j) = xns(i0)
          yfs(i,j) = yns(i0)
          zfs(i,j) = eta(i0-ncar)
        end do
      end do

      ni = nlo
      nj = ntraest+ntraint+ntrapoc
      !write(*,*) 'Pannelli in direzione longitudinale : ',ni
      !write(*,*) 'Pannelli in direzione trasversale   : ',nj

c -- Ordina in x

      do j=1,nj
        do i=1,ni
          xp(i) = xfs(i,j)
          yp(i) = yfs(i,j)
          zp(i) = zfs(i,j)
        end do
        call ordina3(xp,yp,zp,ni)
        do i=1,ni
          xfs(i,j) = xp(i)
          yfs(i,j) = yp(i)
          zfs(i,j) = zp(i)
        end do
      end do

c -- Ordina in y

      do i=1,ni
        do j=1,nj
          xp(j) = xfs(i,j)
          yp(j) = yfs(i,j)
          zp(j) = zfs(i,j)
        end do
        call ordina3(yp,xp,zp,nj)
        do j=1,nj
          xfs(i,j) = xp(j)
          yfs(i,j) = yp(j)
          zfs(i,j) = zp(j)
        end do
      end do

      ifr = nint(fr*1000)
      write(cha,99) ifr
      fina='tatra'//cha//'.ori'
      open(22,file=fina)
      write(22,*) sup,fr,ni,nj
      do i=1,ni
        do j=1,nj
          write(22,*) xfs(i,j),yfs(i,j),zfs(i,j) 
        end do
      end do
      close(22)

c -- Rispazia

      xmin = 0.6
      xmax = xfs(ni,1)
      !write(*,*) 'Estremi X:',xmin,xmax

c -- Passata a j costante

      do j=1,nj

        do i=1,ni
          xp(i) = xfs(i,j)
          yp(i) = yfs(i,j)
          zp(i) = zfs(i,j)
        end do

        d(1) = (yp(2)-yp(1))/(xp(2)-xp(1))
        d(2) = (yp(ni)-yp(ni-1))/(xp(ni)-xp(ni-1))
        call spln1(ni,xp,yp,1,d,cint,wint)
        do i=1,newtra
          xi(i)=xmin+(xmax-xmin)*(i-1.)/(newtra-1.)
          vint(1)=xi(i)
          call spln2(ni,xp,yp,cint,vint)
          yi(i)=vint(2)
        end do

        d(1) = (zp(2)-zp(1))/(xp(2)-xp(1))
        d(2) = (zp(ni)-zp(ni-1))/(xp(ni)-xp(ni-1))
        call spln1(ni,xp,zp,1,d,cint,wint)
        do i=1,newtra
          xi(i)=xmin+(xmax-xmin)*(i-1.)/(newtra-1.)
          vint(1)=xi(i)
          call spln2(ni,xp,zp,cint,vint)
          zi(i)=vint(2)
        end do

        do i=1,newtra
          xfs(i,j) = xi(i)
          yfs(i,j) = yi(i)
          zfs(i,j) = zi(i)
        end do

      end do

      ni = newtra

      ymax = yfs(ni,nj)
      !write(*,*) 'Estremo Y:',ymax

c -- Passata ad i costante

      do i=1,ni

        do j=1,nj
          xp(j) = xfs(i,j)
          yp(j) = yfs(i,j)
          zp(j) = zfs(i,j)
        end do

        do j=1,nj
          xp(nj+j) = xfs(i,j)
          yp(nj+j) =-yfs(i,j)
          zp(nj+j) = zfs(i,j)
        end do

        nj2 = 2*nj

        call ordina3(yp,xp,zp,nj2)

        d(1) = (zp(2)-zp(1))/(yp(2)-yp(1))
        d(2) = (zp(nj2)-zp(nj2-1))/(yp(nj2)-yp(nj2-1))
        call spln1(nj2,yp,zp,1,d,cint,wint)
        do j=1,newtra
          yi(j)=ymax*(j-1.)/(newtra-1.)
          vint(1)=yi(j)
          call spln2(nj2,yp,zp,cint,vint)
          zi(j)=vint(2)
        end do

        do j=1,newtra
          xfs(i,j) = xp(1)
          yfs(i,j) = yi(j)
          zfs(i,j) = zi(j)
        end do

      end do

      nj = newtra

c -- Derivata

      do j=1,nj

        do i=1,ni
          xp(i) = xfs(i,j)
          yp(i) = yfs(i,j)
          zp(i) = zfs(i,j)
        end do

        d(1) = 0.
        d(2) = 0.
        call spln1(ni,xp,zp,2,d,cint,wint)
        do i=1,ni
          vint(1)=xp(i)
          call spln2(ni,xp,zp,cint,vint)
          zx(i,j)=vint(3)
        end do
      end do

c -- Apre il file TATRAxxx.dat

      ifr = nint(fr*1000)
      write(cha,99) ifr
      fina='tatra'//cha//'.dat'
   99 format(i4.4)

      open(unit=20,file=fina,status='unknown')
      write(20,*) sup,fr,nj,ni

c      write(20,*) xfs(1,1),xfs(ni,1)
c--md
      write(20,*) xfs(1,1),xfs(ni,1),1,ni
c--md
      do i=1,ni
        write(20,*) xfs(i,1)
        do j=1,nj
          write(20,*) yfs(i,j),zfs(i,j),zx(i,j) 
        end do
      end do

      close (20)

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
