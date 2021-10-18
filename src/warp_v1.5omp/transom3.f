c-----------------------------------------------------------------------
      subroutine transom3(x1,x2,x3,x4,xh,xn,cosxn,
     $                    y1,y2,y3,y4,yh,yn,cosyn,
     $                    z1,z2,z3,z4,zh,zn,coszn,
     $                    xns,yns,zns,area,vx,vy,vz,
     $                    fi0x,fi0y,fi0z,amat,sigma1,
     $                    nlast,ncar,nulib,nekel,
     &                    fr,dryfr,lontra,diana,itriqua,cata,
     &                    zmaxtra,nplib,ntrapoc,nfile)
c-----------------------------------------------------------------------

      implicit real (a-h,o-z)

      include "warp.cmn"
      parameter  (nti=200)

c-----------------------------------------------------------------------

      integer itriqua(ntotd)
      real  x1(ntotd),y1(ntotd),z1(ntotd),
     &      x2(ntotd),y2(ntotd),z2(ntotd),
     &      x3(ntotd),y3(ntotd),z3(ntotd),
     &      x4(ntotd),y4(ntotd),z4(ntotd)
      real  area(ntotd)
      real  xh(ntotd),yh(ntotd),zh(ntotd)
      real  xn(ntotd),yn(ntotd),zn(ntotd)
      real  xns(ntotd),yns(ntotd),zns(ntotd)
      real  amat(ntotd,ntotd)
      real  cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real  vx(ntotd),vy(ntotd),vz(ntotd)
      real  vx2(ntotd),vy2(ntotd),vz2(ntotd)
      real  fi0x(ntotd),fi0y(ntotd),fi0z(ntotd)
      real  sigma1(ntotd)
      character*1  cata
      character*20 nfile

c-----------------------------------------------------------------------

      integer icanc(nti)
      real yt(nti),zt(nti),tx(nti),ty(nti),
     &     yint(nti),zint(nti),tint(nti),
     &     der(2),cint(3*nti-3),vint(5),wint(3*nti-3)


c      open(20,file='trans.dbg',status='unknown',recl=10000)
c      do i=1,ntotd
c        write(20,*) fi0x(i),fi0y(i),fi0z(i),sigma1(i)
c      end do
c      close(20)
c      open(777,file='elev.dbg',status='unknown',recl=10000)

c-----------------------------------------------------------------------
c     Interpolo coord. ed angoli del transom con le splyne
c-----------------------------------------------------------------------

      frq = fr*fr
      frlim = 1.

      open(unit=20,file=nfile,status='old')
      read(20,*) poppax,alfalat
      do i=1,lontra
        read(20,*)yint(i),zint(i)
      end do
      do i=1,lontra
        read(20,*)tint(i)
      end do
      close(20)

c -- Elimina duplicati

      do i=1,lontra
        icanc(i) = 0
      end do

      do i=1,lontra-1
        do j=i+1,lontra
          if(yint(i).eq.yint(j)) icanc(j) = 1
        end do
      end do

      k = 0
      do i=1,lontra
        if(icanc(i).eq.0) then
          k = k + 1
          yint(k)=yint(i)
          zint(k)=zint(i)
          tint(k)=tint(i)
        end if
      end do
      lontraR = k

c      write(*,*) k,'lontraR'
c      open(20,file='lontra.dbg',status='unknown',recl=10000)
c      do i=1,lontraR
c        write(20,*) yint(i),zint(i),tint(i)
c      end do
c      close(20)


c-----------------------------------------------------------------------
c     Dead water condition
c-----------------------------------------------------------------------

c      write(*,*) fr,dryfr

      if(fr.lt.dryfr) then
        do i=1,lontraR
          zint(i) = zint(i)*(fr/dryfr)**2.
          tint(i) = tint(i)*(fr/dryfr)**2.
        end do
      end if

c      open(20,file='lontra.dbg',status='unknown',recl=10000)
c      do i=1,lontraR
c        write(20,*) yint(i),zint(i),tint(i)
c      end do
c      close(20)


c-----------------------------------------------------------------------
c     ZTRANSOM:  interpolo i valori delle zn.
c-----------------------------------------------------------------------

c -- Z

      der(1)=(zint(2)-zint(1))/(yint(2)-yint(1))
      der(2)= (zint(lontraR)-zint(lontraR-1))/
     &        (yint(lontraR)-yint(lontraR-1))
      
c      write(*,*) der(1),der(2)

      call spln1(lontraR,yint,zint,1,der,cint,wint)
      do i=1,ntrapoc
        yt(i)=yn(nlast+i)
        vint(1)=yt(i)
        call spln2(lontraR,yint,zint,cint,vint)
        zt(i)=vint(2)
        ty(i)=vint(3)
      end do

c -- tx

      der(1)=(tint(2)-tint(1))/(yint(2)-yint(1))
      der(2)=(tint(lontraR)-tint(lontraR-1))/
     &       (yint(lontraR)-yint(lontraR-1))
      call spln1(lontraR,yint,tint,1,der,cint,wint)
      do i=1,ntrapoc
        vint(1)=yt(i)
        call spln2(lontraR,yint,tint,cint,vint)
        tx(i)=vint(2)
        if(tx(i).gt. 0.1) tx(i) =  0.1
        if(tx(i).lt.-0.1) tx(i) = -0.1
      end do

      tmed = 0.
      do i=1,ntrapoc
        tmed = tmed + tx(i)
      end do
      tmed = tmed/float(ntrapoc)

      open(38,file='transom.int',status='unknown',form='formatted')
      do i=1,ntrapoc
        write(38,*) poppax,yt(i),zt(i),tx(i),tmed
        tx(i) = tmed
      end do
      close(38)

c-----------------------------------------------------------------------
c     Condizione transom: impongo l'elevazione estrapolando con Taylor
c-----------------------------------------------------------------------
      if(nekel.eq.1) then    !! nekel
c-----------------------------------------------------------------------
        do 345 i=nlast+1,nlast+ntrapoc     ! main loop 
          icont = i-nlast                  ! icont va da 1 ---> ntrapoc
          ii    = ncar+nplib-ntrapoc+icont ! contatore per l'ultima fila
          in2   = i+ntrapoc                ! contatore per la fila 2
          elev  = zt(icont)+tx(icont)*(xns(i)-poppax)
          elev2 = zt(icont)+tx(icont)*(xns(in2)-poppax)
c----------------------------------------------------------------------*
c        velocita' nel punto 1
c----------------------------------------------------------------------*

          call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &                 y1,y2,y3,y4,yh,yn,cosyn,
     &                 z1,z2,z3,z4,zh,zn,coszn,diana,
     &                 area,vx,vy,vz,ncar,
     &                 xns(i),yns(i),zns(i),i,1,ncar+nplib,itriqua,
     &                 ncard,ntotd,nplib)
c----------------------------------------------------------------------*
c        velocita' nel punto 2
c----------------------------------------------------------------------*

          call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &                 y1,y2,y3,y4,yh,yn,cosyn,
     &                 z1,z2,z3,z4,zh,zn,coszn,diana,
     &                 area,vx2,vy2,vz2,ncar,
     &                 xns(in2),yns(in2),zns(in2),
     &                 in2,1,ncar+nplib,itriqua,
     &                 ncard,ntotd,nplib)

c----------------------------------------------------------------------*

          dx   = xn(in2)-xn(i) ! dx e' calcolato tra 2 file di pannelli

          do 320 k=1,ncar+nplib

            if (fr.lt.frlim) then

c -- elevazione

              amat(i,k)   = vx(k)
              amat(in2,k) = vx2(k)

c -- Cond. di Sup.Lib. unificata

              tx1 = zt(icont)/frq
              tx2 = 1.0
              ty2 = 0.0
              tz1 = dx/frq
              tz2 = 0.0
              amat(ii,k) = vx(k)*tx1+vx2(k)*tx2+vz(k)*tz1

            else

c -- Elevazione

              amat(i,k)   = vx(k)

            end if

 320      continue 

          if (fr.lt.frlim) then

c -- elevazione

            sigma1(i)   = -elev/frq
            sigma1(in2) = -elev2/frq

c -- Cond. di Sup.Lib. unificata

            sigma1(ii) = -zt(icont)/frq

          else

            sigma1(i)   = -zt(icont)/frq

          end if

  345   continue
c
c----------------------------------------------------------------------*
      else         !! dawson
c----------------------------------------------------------------------*
c     la coda va da [nlast+1---->ncar+nplib]
c----------------------------------------------------------------------*

        do 350 i=nlast+1,nlast+ntrapoc      ! main loop 
          icont = i-nlast                   ! icont va da 1 ---> ntrapoc
          ii    = ncar+nplib-ntrapoc+icont  ! contatore per l'ultima fila
          in2   = i+ntrapoc                 ! contatore per la fila 2

c        USO UNA DERIVATA PRIMA A 2 PUNTI CENTRATA

          elev =zt(icont)+tx(icont)*(xn(i)-poppax)
          elev2=zt(icont)+tx(icont)*(xn(in2)-poppax)

c            write(777,*) zt(icont),tx(icont),xn(i),xn(in2)

          fi021  = fi0x(i)**2  +fi0y(i)**2  +fi0z(i)**2
          fi022  = fi0x(in2)**2+fi0y(in2)**2+fi0z(in2)**2
          fixp1  = fi0x(i)+1.0
          fix2p1 = fi0x(in2)+1.0
          dx     = xn(in2)-xn(i) ! dx e' calcolato tra 2 file di pannelli

c----------------------------------------------------------------------*
c        velocita' nel punto 1
c----------------------------------------------------------------------*

          call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &                 y1,y2,y3,y4,yh,yn,cosyn,
     &                 z1,z2,z3,z4,zh,zn,coszn,diana,
     &                 area,vx,vy,vz,ncar,
     &                 xns(i),yns(i),zns(i),i,1,ncar+nplib,itriqua,
     &                 ncard,ntotd,nplib)

c----------------------------------------------------------------------*
c        velocita' nel punto 2
c----------------------------------------------------------------------*

          call caldino(x1,x2,x3,x4,xh,xn,cosxn,
     &                 y1,y2,y3,y4,yh,yn,cosyn,
     &                 z1,z2,z3,z4,zh,zn,coszn,diana,
     &                 area,vx2,vy2,vz2,ncar,
     &                 xns(in2),yns(in2),zns(in2),
     &                 in2,1,ncar+nplib,itriqua,
     &                 ncard,ntotd,nplib)

c----------------------------------------------------------------------*
c===============================================================
          do 520 k=1,ncar+nplib

            if(fr.lt.frlim) then

c --  elevazione

              amat(i,k)=vx(k)*fixp1+vy(k)*fi0y(i)+vz(k)*fi0z(i)
              amat(in2,k)=vx2(k)*fix2p1+
     &                    vy2(k)*fi0y(in2)+
     &                    vz2(k)*fi0z(in2)

c --  Cond. di Sup.Lib. unificata

              tx1=0.5*fi022+fi0x(in2)+zt(icont)/frq
              tx2=fixp1*(1.0+fi0x(in2))
              ty2=fi0y(in2)*fixp1
              tz1=dx/frq
              tz2=fi0z(in2)*fixp1
              amat(ii,k)=vx(k)*tx1+vx2(k)*tx2+
     &                   vy2(k)*ty2+vz(k)*tz1+vz2(k)*tz2

            else

c -- Elevazione

              amat(i,k)=vx(k)*fixp1+vy(k)*fi0y(i)+vz(k)*fi0z(i)

            end if

 520      continue

          if(fr.lt.frlim) then

c --  elevazione

            sigma1(i)   = -(elev/frq + fi021*0.5 + fi0x(i))
            sigma1(in2) = -(elev2/frq + fi022*0.5 + fi0x(in2))

c --  Cond. di Sup.Lib. unificata

            sigma1(ii) = -fixp1*(0.5*fi022+zt(icont)/frq+fi0x(in2))-
     &                    dx/frq*fi0z(i)

          else

c -- Elevazione

            sigma1(i)   = -(zt(icont)/frq + fi021*0.5 + fi0x(i))

          end if

  350   continue
c===============================================================
      end if

c      close(777)

c-----------------------------------------------------------------------
      return
      end 
c-----------------------------------------------------------------------

      subroutine ordina3(x,y,z,n)

c-----------------------------------------------------------------------
c     Ordina tre vettori in base al primo
c-----------------------------------------------------------------------

      integer n,ni,i,j
      real x(*),y(*),z(*)
      real appo

c-----------------------------------------------------------------------

      if(n.gt.0) then
        do i=1,n-1
          do j=i+1,n
            if(x(j).lt.x(i)) then
              appo = x(i)
              x(i) = x(j)
              x(j) = appo
              appo = y(i)
              y(i) = y(j)
              y(j) = appo
              appo = z(i)
              z(i) = z(j)
              z(j) = appo
            end if
          end do
        end do

      else

        ni = -n
        do i=1,ni-1
          do j=i+1,ni
            if(x(j).gt.x(i)) then
              appo = x(i)
              x(i) = x(j)
              x(j) = appo
              appo = y(i)
              y(i) = y(j)
              y(j) = appo
              appo = z(i)
              z(i) = z(j)
              z(j) = appo
            end if
          end do
        end do
      end if

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------

      subroutine filtra(x,y,np)

c-----------------------------------------------------------------------

      real x(*),y(*)
      real yf(1000)

c-----------------------------------------------------------------------

      c1 = 1./17.
      c2 = 3./17.
      c3 = 9./17.
      c4 = 3./17.
      c5 = 1./17.

      do iter=1,5
        do j=3,np-2
          yf(j) = c1*y(j-2) +
     &            c2*y(j-1) +
     &            c3*y(j  ) +
     &            c4*y(j+1) +
     &            c5*y(j+2)
        end do
        do j=3,np-2
          y(j) = yf(j)
        end do
      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
