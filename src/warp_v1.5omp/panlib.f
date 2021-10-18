c-----------------------------------------------------------
      subroutine panlib(ycamax,
     &                  xn,yn,zn,
     &                  x1,x2,x3,x4,xh,
     &                  y1,y2,y3,y4,yh,
     &                  z1,z2,z3,z4,zh,
     &                  cosxn,cosyn,coszn,
     &                  area,ncar,alfa,ncard,nplib,nun,
     &                  intra,intraest,intraint,intrapoc,invallec,
     &                  inlo,ntotd,nplibd)
c
c----------------------------------------------------------------
c   Griglia per la superficie libera.
c----------------------------------------------------------------
      implicit real (a-h,o-z)
      real x1(ntotd),y1(ntotd),z1(ntotd),
     &     x2(ntotd),y2(ntotd),z2(ntotd),
     &     x3(ntotd),y3(ntotd),z3(ntotd),
     &     x4(ntotd),y4(ntotd),z4(ntotd)
      real area(ntotd)
      real xh(ntotd),yh(ntotd),zh(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
c-----------------------------------------------------------------------

      !write(*,*)'----------========== P A N L I B ==========----------'
      pig      = acos(-1.0)

c----------------------------------------------------------------------*
c     Legge il file di griglia
c----------------------------------------------------------------------*
      open (50,file='fsgrid.grd',form='formatted',status='old')
      read (50,*) intraint,intraest,inlo,invallec,intrapoc,idum,alfa
c
       intra=intraint+intraest
       nplib=(intra*inlo+intrapoc*invallec)
       itota=nplib+ncar
       itest=ntotd-itota
         if(itest.lt.0) then 
         !write(*,*) '-----===== ATTENZIONE !!!!! =====-----'    
         !write(*,*) 'Disaccordo grave tra la parameter di warp '
         !write(*,*) 'ed il file di griglia fsgrid.grd '
         !write(*,*) ' ntotd   = ',ntotd,  ' itotali = ',itota
         !write(*,*) ' nplibd  = ',nplibd, ' nplib   = ',nplib
         stop 
       end if
c --
       xmonte =  1.e+20
       xvalle = -1.e+20
       ylato  = -1.e+20
c --

      do i=ncar+1,ncar+intra*inlo+invallec*intrapoc
c
c  lettura punti della griglia
c
        read(50,*) x1(i),x2(i),x3(i),x4(i)
        read(50,*) y1(i),y2(i),y3(i),y4(i)
        read(50,*) z1(i),z2(i),z3(i),z4(i)
c
c  punto di controllo
c
        xn(i) = 0.25*(x1(i)+x2(i)+x3(i)+x4(i))
        yn(i) = 0.25*(y1(i)+y2(i)+y3(i)+y4(i))
        zn(i) = 0.25*(z1(i)+z2(i)+z3(i)+z4(i))
c
c   punto h
c
        hmod34 = (x3(i)-x4(i))**2 
     &          +(y3(i)-y4(i))**2 
     &          +(z3(i)-z4(i))**2 
        scal34 = (x1(i)-x4(i))*(x3(i)-x4(i)) 
     &          +(y1(i)-y4(i))*(y3(i)-y4(i)) 
     &          +(z1(i)-z4(i))*(z3(i)-z4(i)) 
        scal34 = scal34/hmod34
        xh(i) = x4(i) + scal34*(x3(i)-x4(i))
        yh(i) = y4(i) + scal34*(y3(i)-y4(i))
        zh(i) = z4(i) + scal34*(z3(i)-z4(i))
c
c  normale
c
        cosxn(i) =  0.0
        cosyn(i) =  0.0
        coszn(i) = -1.0 

c
c  area pannello
c
        a13x = x1(i) - x3(i)
        a13y = y1(i) - y3(i)
        a13z = z1(i) - z3(i)
        a24x = x2(i) - x4(i)
        a24y = y2(i) - y4(i)
        a24z = z2(i) - z4(i)
        ax   = a13y*a24z - a13z*a24y
        ay   = a13z*a24x - a13x*a24z
        az   = a13x*a24y - a13y*a24x
        area(i) = ax*ax+ay*ay+az*az
        if(area(i).ne.0) area(i) = 0.5*sqrt(area(i))
c
      end do
      close(50)
c
      return  
      end
c-------------------------------------------------------------------
