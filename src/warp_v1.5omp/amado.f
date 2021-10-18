c-----------------------------------------------------------------------

      subroutine amado(x1,x2,x3,x4,xn,xh,y1,y2,y3,y4,yn,yh,
     &                 z1,z2,z3,z4,zn,zh,cosxn,cosyn,coszn,
     &                 xns,yns,zns,
     &                 area,fi0x,fi0y,fi0z,vx,vy,vz,
     &                 amat,sigma0,indx,appo,diana,nekel,
     &                 ncar,nlast,eps,itriqua,modo,
     &                 nplib,ntra)

c-----------------------------------------------------------------------

      implicit real (a-h,o-z)

      include "warp.cmn"

      integer itriqua(ntotd)
      integer indx(ntotd)

      real x1(ntotd),y1(ntotd),z1(ntotd),
     &     x2(ntotd),y2(ntotd),z2(ntotd),
     &     x3(ntotd),y3(ntotd),z3(ntotd),
     &     x4(ntotd),y4(ntotd),z4(ntotd)
      real area(ntotd)
      real xh(ntotd),yh(ntotd),zh(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real vx(ntotd),vy(ntotd),vz(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real xns(ntotd),yns(ntotd),zns(ntotd)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd)
      real amat(ntotd,ntotd)
      real sigma0(ntotd)
      real r(ntotd)
c
      integer iwork(ntotd)
      real appo(ntotd)
      real  work(4*(ntotd))
      character*1 trans

c --

      write(nun,*)
     &  '----------======== B A S E F L O W ========----------'
      !write( * ,*)
      !&  '----------======== B A S E F L O W ========----------'
  10  if (nekel.eq.1) then
        !write(*,*)'                  Neumann Kelvin'
        write(nun,*)'                 Neumann Kelvin'
        do i=1,ncar+nplib
          sigma0(i) = 0.
          fi0x(i)   = 0.
          fi0y(i)   = 0.
          fi0z(i)   = 0.
        end do
        go to 330
      else if (nekel.eq.2) then
        !write(*,*)'                  Sigma = coseno'
        write(nun,*)'                 Sigma = coseno'
        do i=1,ncar+nplib
          sigma0(i) = -cosxn(i)
        end do
        go to 304
      else
        !write(*,*)'                  Double model'
        write(nun,*)'                 Double model'
      end if

      !write(nun,*)
      !&  '----------======== L E T U S G O ! ========----------'
      !write( * ,*)
      !&  '----------======== L E T U S G O ! ========----------'
c ----------------------------------------------------------------------
c     calcolo delle sigma0  --- doppio modello
c ----------------------------------------------------------------------

      do 300 i=1,ncar
c
c-----------------------------------------------------------------------
c     INDOMO calcola per il pannello i le velocita' indotte dai
c     pannelli della carena 
c-----------------------------------------------------------------------
c
        call indomo(i,xns(i),yns(i),zns(i),1,ncar,diana,xn,yn,zn,
     $              x1,x2,x3,x4,xh,y1,y2,y3,y4,yh,z1,z2,z3,z4,zh,
     $              area,cosxn,cosyn,coszn,vx,vy,vz,itriqua,modo,
     $              ncard,ntotd,nplib)
c
        do k=1,ncar
          amat(i,k) = cosxn(i)*vx(k)+cosyn(i)*vy(k)+coszn(i)*vz(k)
        end do
c
        sigma0(i) = -cosxn(i)

  300 continue

c      open(20,file='indomo_vel.dbg',status='unknown',recl=10000)
c        do i=1,ncar+nplib
c        write(20,*) vx(i),vy(i),vz(i)
c        end do
c      close(20)

      !write(*,*)   'Precisione macchina    =',eps
      write(nun,*) 'Precisione macchina    =',eps
c-----------------------------------------------------------------------
c            INVERSIONE DELLA MATRICE
c-----------------------------------------------------------------------
      !write(*,*) 'Soluzione sistema'
      nsys = ncar
      sing = 1.E-5
      dx   = 1.E-5
      alpha = 1./float(nsys)
      do i=1,nsys
        r(i) = sigma0(i)
      end do
c -- Traspone
      do i=1,nsys-1
        do j=i+1,nsys
          raux = amat(i,j)
          amat(i,j) = amat(j,i)
          amat(j,i) = raux
        end do
      end do

c-- md
c      call simqit(amat,r,nsys,ntotd,sing,alpha,dx,sigma0,ierr)
      call simqit(amat,r,nsys,ntotd,alpha,dx,sigma0,ierr)
c-- md

c-----------------------------------------------------------------------
c Accuratezza della soluzione
c-----------------------------------------------------------------------
  304 ss0=0.0
      do 305 i=1,ncar
        sigma0(i) = float(modo)*sigma0(i)
        ss0=ss0+sigma0(i)*area(i)
  305 continue
      write(nun,*)  'Integrale delle sorgenti sul corpo = ',ss0
      !write( * ,*)  'Integrale delle sorgenti sul corpo = ',ss0

  310 continue
c
c-----------------------------------------------------------------------
c     calcolo di fi0x, fi0y, fi0z 
c-----------------------------------------------------------------------
c

c        write(*,*) 'amado'
c        open(20,file='f0.dbg',status='unknown',recl=10000)
c        write(*,*) ncar+nplib

      do 320 i=1,ncar+nplib
        fi0x(i)= 0.0
        fi0y(i)= 0.0
        fi0z(i)= 0.0
c------------------------------
c Controlli per la singolarita'
c------------------------------
        if(i.le.ncar) then                          !! Carena senza shift
          ising=i
        else
          ising=0
        end if
c
        call indomo(ising,xns(i),yns(i),zns(i),1,ncar,diana,xn,yn,zn,
     $              x1,x2,x3,x4,xh,y1,y2,y3,y4,yh,z1,z2,z3,z4,zh,
     $              area,cosxn,cosyn,coszn,vx,vy,vz,itriqua,modo,
     $              ncard,ntotd,nplib)
c

        do k=1,ncar
          fi0x(i)= fi0x(i)+sigma0(k)*vx(k)
          fi0y(i)= fi0y(i)+sigma0(k)*vy(k)
          fi0z(i)= fi0z(i)+sigma0(k)*vz(k)
        end do
c

c       write(20,*) fi0x(i),fi0y(i),fi0z(i)

  320 continue

c       close(20) 
c
c ----------------------------------------------------------------------
  330 return
      end
c ----------------------------------------------------------------------
