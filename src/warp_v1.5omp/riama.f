c-----------------------------------------------------------------------
      subroutine riama(x1,x2,x3,x4,xn,xh,y1,y2,y3,y4,yn,yh,
     &                 z1,z2,z3,z4,zn,zh,cosxn,cosyn,coszn,
     &                 xns,yns,zns,nun,
     &                 area,fi0x,fi0y,fi0z,vx,vy,vz,fi0ll,vll,
     &                 amat,sigma0,sigma1,
     &                 diana,nekel,fr,itriqua,modo,
     &                 ncar,nplib,nlast,ntra,
     &                 ncard,ntotd)
c-----------------------------------------------------------------------
      integer ncard,nplib,ntotd
      integer itriqua(ntotd)
      real x1(ntotd),y1(ntotd),z1(ntotd),
     &     x2(ntotd),y2(ntotd),z2(ntotd),
     &     x3(ntotd),y3(ntotd),z3(ntotd),
     &     x4(ntotd),y4(ntotd),z4(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real xh(ntotd),yh(ntotd),zh(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real xns(ntotd),yns(ntotd),zns(ntotd)
      real area(ntotd)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd)
      real vx(ntotd),vy(ntotd),vz(ntotd)
      real vll(ntotd)
      real amat(ntotd,ntotd)
      real sigma1(ntotd),sigma0(ntotd)
c
      integer i,k,nun,ntra,ncar,nlast,modo,ising,nekel
      real fi0l2,fi0l,fi0ll,cosxl,cosyl,coszl,fr,diana
c-----------------------------------------------------------------------
      !write( * ,*)
      !&  '-----------------------------------------------------'
      !write( * ,*)'                   Froude =',fr
      !write( * ,*)
      !&  '-----------------------------------------------------'
      write(nun,*)
     &  '-----------------------------------------------------'
      write(nun,*)'                  Froude =',fr
      write(nun,*)
     &  '-----------------------------------------------------'
c-----------------------------------------------------------------------
c  Matrice AMAT: condizione sul corpo: I = 1, NCAR : K = 1, NCAR+NPLIB
c  La sub. CALDINO restituisce vx,vy,vz, influenza dei pannelli j su i
c-----------------------------------------------------------------------
c
      !write(*,*)'Condizione di impermeabilita'' sul corpo'
c

      do 335 i=1,ncar
        call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &               y1,y2,y3,y4,yh,yn,cosyn,
     &               z1,z2,z3,z4,zh,zn,coszn,diana,
     &               area,vx,vy,vz,ncar,
     &               xns(i),yns(i),zns(i),i,1,ncar+nplib,itriqua,
     &               ncard,ntotd,nplib)
c        !write(*,fmt='($,A,I4.4)')'\b\b\b\b\b', i
        do k=1,ncar+nplib
           amat(i,k)=vx(k)*cosxn(i)+vy(k)*cosyn(i)+vz(k)*coszn(i)
        end do
        if (nekel.eq.1) then
          sigma1(i) = -cosxn(i)  ! Termine noto nk
        else
          sigma1(i) = 0.0   ! Termine noto dawson
        end if

  335 continue

c      open(20,file='cald1_vel.dbg',status='unknown',recl=10000)
c        do i=1,ncar+nplib
c        write(20,*) vx(i),vy(i),vz(i)
c        end do
c      close(20)

c-----------------------------------------------------------------------
c  Matrice AMAT: condizione sulla superficie libera:
c  I = NCAR+1,NCAR+NPLIB : K = 1, NCAR+NPLIB
c  La sub. CALDINO restituisce vx,vy,vz, influenza dei pannelli j su i
c  La sub. FIDELL  restituisce fi0ll e vll     
c-----------------------------------------------------------------------
      !write(*,*)'Condizione di superficie libera'
c

      do 123 i = ncar+1,ncar+nplib  ! con la coda
      if(i.gt.ncar.and.i.le.ncar+ntra) then        !! Prima riga
           ising=0
      else if(i.gt.ncar+ntra.and.i.le.nlast) then  !! Sup.lib
           ising=i-ntra
      else if(i.gt.nlast) then                     !! Coda non shiftata
           ising=i
      end if
c
      call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &             y1,y2,y3,y4,yh,yn,cosyn,
     &             z1,z2,z3,z4,zh,zn,coszn,diana,
     &             area,vx,vy,vz,ncar,
     &             xns(i),yns(i),zns(i),ising,1,ncar+nplib,itriqua,
     &             ncard,ntotd,nplib)
c
c      !write(*,fmt='($,A,I4.4)')'\b\b\b\b\b', i
c
      call fidell(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
     &            xh,yh,zh,xn,yn,zn,area,cosxn,cosyn,coszn,sigma0,
     &            xns,yns,zns,
     &            fi0x,fi0y,fi0z,diana,fi0ll,vll,
     &            i,ncar,nekel,itriqua,modo,ncard,ntotd,nplib)

c
c... Calcolo di cosxl,cosyl per linearizzazione tipo Dawson
c
      fi0l2 = (1.0+fi0x(i))**2+fi0y(i)**2+fi0z(i)**2
      if(fi0l2.gt.0) then
        fi0l = sqrt(fi0l2)
      else
        !write(*,*) 'WARNING! FHI_L ?'
        fi0l = 1.
      end if
      cosxl = (1.0+fi0x(i))/fi0l
      cosyl = fi0y(i)/fi0l
      coszl = fi0z(i)/fi0l
c
c... Riempie la matrice
c
      do k = 1,ncar+nplib     
         if(nekel.eq.1) then
           amat(i,k) = vll(k) + vz(k)/fr**2 
         else 
           amat(i,k)=fi0l2*vll(k) +
     &             2.0*fi0l*fi0ll*(cosxl*vx(k)+cosyl*vy(k)+coszl*vz(k))+
     &               vz(k)/fr**2
         end if
      end do
c
c... Costruisce il termine noto
c
      if(nekel.eq.1) then
        sigma1(i) = 0.
      else
        sigma1(i) = - fi0l2 * fi0ll - fi0z(i) / fr**2
      end if
c
 123  end do

c      open(20,file='cald2_vel.dbg',status='unknown',recl=10000)
c        do i=1,ncar+nplib
c        write(20,*) vx(i),vy(i),vz(i),sigma1(i)
c        end do
c      close(20)


c-----------------------------------------------------------------------
      return 
      end
c-----------------------------------------------------------------------
