c----------------------------------------------------------------------*
c     costruzione delle matrici di influenza 
c----------------------------------------------------------------------*
      subroutine caldino(xa,xb,xc,xd,xx,xn,cosxq,
     &                   ya,yb,yc,yd,yy,yn,cosyq,
     &                   za,zb,zc,zd,zz,zn,coszq,diana,
     &                   arpa,vvx,vvy,vvz,ncar,
     &                   xxx,yyy,zzz,ising,kin,kfin,itriqua,
     &                   ncard,ntotd,nplib) 
c----------------------------------------------------------------------*
c   Calcolo l'influenza dei pannelli k sul pannello i
c   Il pannello sorgente e' (xx,yy,zz) , il punto di collocazione 
c   e' shiftato a monte e si trova in (xxx,yyy,zzz)
c----------------------------------------------------------------------*

      implicit real (a-h,o-z)

      integer itriqua(ntotd)
      real xa(ntotd),ya(ntotd),za(ntotd),
     &     xb(ntotd),yb(ntotd),zb(ntotd),
     &     xc(ntotd),yc(ntotd),zc(ntotd),
     &     xd(ntotd),yd(ntotd),zd(ntotd),
     &     xx(ntotd),yy(ntotd),zz(ntotd)
      real arpa(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real cosxq(ntotd),cosyq(ntotd),coszq(ntotd)
      real vvx(ntotd),vvy(ntotd),vvz(ntotd)

c-----------------------------------------------------------------------

      pig  = acos(-1.0)
c
c=======================================================================
c      caratteristiche geom. del pannello sorgente
c=======================================================================
c
c$omp parallel
c$omp do private(k,w123,w341,dilato,drelpq,iikk,
c$omp&           x12,y12,z12,
c$omp&           x23,y23,z23,
c$omp&           x13,y13,z13,
c$omp&           x34,y34,z34,
c$omp&           x41,y41,z41,
c$omp&           u1,u2,
c$omp&           v1,v2,
c$omp&           w1,w2)
       do 200 k = kin,kfin
c++++++++++++++++++++++++++++++++++++++++
          if (itriqua(k).eq.3) then 
c++++++++++++++++++++++++++++++++++++++++
              call primatri(xa(k),xb(k),xc(k),ya(k),yb(k),yc(k),
     $            za(k),zb(k),zc(k),xx(k),yy(k),zz(k),arpa(k),
     $            w123,dilato,x12,y12,z12,x23,y23,z23,x13,y13,z13)
c++++++++++++++++++++++++++++++++++++++++
          else
c++++++++++++++++++++++++++++++++++++++++
              call primaqua(xa(k),xb(k),xc(k),xd(k),ya(k),yb(k),yc(k),
     $         yd(k),za(k),zb(k),zc(k),zd(k),xx(k),yy(k),zz(k),arpa(k),
     $         w123,w341,dilato,x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $         x34,y34,z34,x41,y41,z41)
c++++++++++++++++++++++++++++++++++++++++
          end if
c++++++++++++++++++++++++++++++++++++++++
c
        drelpq=sqrt((xxx-xn(k))**2+(yyy-yn(k))**2+(zzz-zn(k))**2) 
     $         /dilato
c
c=======================================================================
        if (drelpq.ge.diana) then
c=======================================================================
c++++++++++++++++++++++++++++++++++++++++
        if (itriqua(k).eq.3) then
c++++++++++++++++++++++++++++++++++++++++
           call splyntri(x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $                   xxx,yyy,zzz,w123,u1,v1,w1)
           call splyntri(x12,-y12,z12,x23,-y23,z23,x13,-y13,z13,
     $                   xxx,yyy,zzz,w123,u2,v2,w2)
c++++++++++++++++++++++++++++++++++++++++
        else
c++++++++++++++++++++++++++++++++++++++++
           call splynqua(x12,y12,z12,x23,y23,z23,x13,y13,z13,x34,
     $            y34,z34,x41,y41,z41,xxx,yyy,zzz,w123,w341,u1,v1,w1)
           call splynqua(x12,-y12,z12,x23,-y23,z23,x13,-y13,z13,x34,
     $          -y34,z34,x41,-y41,z41,xxx,yyy,zzz,w123,w341,u2,v2,w2)
c++++++++++++++++++++++++++++++++++++++++
        end if
c++++++++++++++++++++++++++++++++++++++++
        vvx(k) = u1 + u2
        vvy(k) = v1 + v2
        vvz(k) = w1 + w2
c=======================================================================
         else
c=======================================================================
c
	  if(ising.eq.k) then
            iikk=1       
	  else
	    iikk=0
	  end if 
c
c++++++++++++++++++++++++++++++++++++++++
        if (itriqua(k).eq.3) then
c++++++++++++++++++++++++++++++++++++++++
          call hemitri(xn(k),yn(k),zn(k),xxx,yyy,zzz,
     $             xa(k),xb(k),xc(k),xx(k),yy(k),zz(k),
     $             ya(k),yb(k),yc(k),za(k),zb(k),zc(k),
     $             cosxq(k),cosyq(k),coszq(k),u1,v1,w1,pig,iikk)
          iikk=0
          call hemitri(xn(k),-yn(k),zn(k),xxx,yyy,zzz,
     $            xa(k),xb(k),xc(k),xx(k),-yy(k),zz(k),
     $           -ya(k),-yb(k),-yc(k),za(k),zb(k),zc(k),
     $            cosxq(k),-cosyq(k),coszq(k),u2,v2,w2,pig,iikk)
c++++++++++++++++++++++++++++++++++++++++
        else
c++++++++++++++++++++++++++++++++++++++++
          call hemiqua(xn(k),yn(k),zn(k),xxx,yyy,zzz,
     $             xa(k),xb(k),xc(k),xd(k),xx(k),yy(k),zz(k),
     $             ya(k),yb(k),yc(k),yd(k),za(k),zb(k),zc(k),zd(k),
     $             cosxq(k),cosyq(k),coszq(k),u1,v1,w1,pig,iikk)
          iikk=0
          call hemiqua(xn(k),-yn(k),zn(k),xxx,yyy,zzz,
     $            xa(k),xb(k),xc(k),xd(k),xx(k),-yy(k),zz(k),
     $           -ya(k),-yb(k),-yc(k),-yd(k),za(k),zb(k),zc(k),zd(k),
     $            cosxq(k),-cosyq(k),coszq(k),u2,v2,w2,pig,iikk)
c++++++++++++++++++++++++++++++++++++++++
        end if
c++++++++++++++++++++++++++++++++++++++++
        vvx(k) = u1 + u2
        vvy(k) = v1 + v2
        vvz(k) = w1 + w2
c
c=======================================================================
        end if
c=======================================================================
        call checknan(vvx(k))
        call checknan(vvy(k))
        call checknan(vvz(k))
c
  200 continue 

c$omp end do
c$omp end parallel
c 
c ----------------------------------------------------------------------
       return
       end 
c ----------------------------------------------------------------------
