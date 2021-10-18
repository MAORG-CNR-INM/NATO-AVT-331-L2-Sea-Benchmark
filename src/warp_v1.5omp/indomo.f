c
      subroutine indomo(ising,xni,yni,zni,nino,nfino,diana,xn,yn,zn,
     &              x1,x2,x3,x4,xh,y1,y2,y3,y4,yh,z1,z2,z3,z4,zh,
     &              area,cosxn,cosyn,coszn,vx,vy,vz,itriqua,modo,
     &              ncard,ntotd,nplib)
c 
      implicit real (a-h,o-z)
c
c-----------------------------------------------------------------------
      integer itriqua(ntotd)
      real x1(ntotd),y1(ntotd),z1(ntotd),
     &     x2(ntotd),y2(ntotd),z2(ntotd),
     &     x3(ntotd),y3(ntotd),z3(ntotd),
     &     x4(ntotd),y4(ntotd),z4(ntotd)
      real area(ntotd)
      real xh(ntotd),yh(ntotd),zh(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real vx(ntotd),vy(ntotd),vz(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
c-----------------------------------------------------------------------
c
c     ISING e' l'indice del pannello singolare
c
      pig=2.0*2.0*atan(1.0)
c   
c-----------------------------------------------------------------------
c    calcolo degli integrali 'corpo-corpo' vx,vy,vz
c-----------------------------------------------------------------------
c
c$omp parallel
c$omp do private(k,w123,w341,dilato,drelpq,iikk,
c$omp&           x12,y12,z12,
c$omp&           x23,y23,z23,
c$omp&           x13,y13,z13,
c$omp&           x34,y34,z34,
c$omp&           x41,y41,z41,
c$omp&           u1,u2,u3,u4,
c$omp&           v1,v2,v3,v4,
c$omp&           w1,w2,w3,w4)
      do k = nino, nfino               ! pannello sorgente 
c++++++++++++++++++++++++++++++++++++++++
       if(itriqua(k).eq.3) then 
c++++++++++++++++++++++++++++++++++++++++
          call primatri(x1(k),x2(k),x3(k),y1(k),y2(k),y3(k),
     $            z1(k),z2(k),z3(k),xh(k),yh(k),zh(k),area(k),
     $            w123,dilato,x12,y12,z12,x23,y23,z23,x13,y13,z13)
c++++++++++++++++++++++++++++++++++++++++
       else
c++++++++++++++++++++++++++++++++++++++++
         call primaqua(x1(k),x2(k),x3(k),x4(k),y1(k),y2(k),y3(k),y4(k),
     $            z1(k),z2(k),z3(k),z4(k),xh(k),yh(k),zh(k),area(k),
     $            w123,w341,dilato,x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $            x34,y34,z34,x41,y41,z41)
c++++++++++++++++++++++++++++++++++++++++
       end if
c++++++++++++++++++++++++++++++++++++++++
c
        drelpq=sqrt((xni-xn(k))**2+(yni-yn(k))**2+(zni-zn(k))**2) 
     $         / dilato
c
c=======================================================================
         if (drelpq.ge.diana) then
c=======================================================================
c
c++++++++++++++++++++++++++++++++++++++++
       if(itriqua(k).eq.3) then 
c++++++++++++++++++++++++++++++++++++++++
        call splyntri(x12,y12,z12,x23,y23,z23,
     $                x13,y13,z13,xni,yni,zni,
     $                w123,u1,v1,w1)
        call splyntri(x12,-y12,z12,x23,-y23,z23,
     $                x13,-y13,z13,xni,yni,zni,
     $                w123,u2,v2,w2)
        call splyntri(x12,y12,-z12,x23,y23,-z23,
     $                x13,y13,-z13,xni,yni,zni,
     $                w123,u3,v3,w3)
        call splyntri(x12,-y12,-z12,x23,-y23,-z23,
     $                x13,-y13,-z13,xni,yni,zni,
     $                w123,u4,v4,w4)
         vx(k) = u1 + u2 + float(modo)*(u3+u4) 
         vy(k) = v1 + v2 + float(modo)*(v3+v4)
         vz(k) = w1 + w2 + float(modo)*(w3+w4)
c++++++++++++++++++++++++++++++++++++++++
       else
c++++++++++++++++++++++++++++++++++++++++
        call splynqua(x12,y12,z12,x23,y23,z23,
     $             x13,y13,z13,x34,y34,z34,
     $             x41,y41,z41,xni,yni,zni,
     $             w123,w341,u1,v1,w1)
        call splynqua(x12,-y12,z12,x23,-y23,z23,
     $             x13,-y13,z13,x34,-y34,z34,
     $             x41,-y41,z41,xni,yni,zni,
     $             w123,w341,u2,v2,w2)
        call splynqua(x12,y12,-z12,x23,y23,-z23,
     $             x13,y13,-z13,x34,y34,-z34,
     $             x41,y41,-z41,xni,yni,zni,
     $             w123,w341,u3,v3,w3)
        call splynqua(x12,-y12,-z12,x23,-y23,-z23,
     $             x13,-y13,-z13,x34,-y34,-z34,
     $             x41,-y41,-z41,xni,yni,zni,
     $             w123,w341,u4,v4,w4)
         vx(k) = u1 + u2 + float(modo)*(u3+u4) 
         vy(k) = v1 + v2 + float(modo)*(v3+v4)
         vz(k) = w1 + w2 + float(modo)*(w3+w4)
c++++++++++++++++++++++++++++++++++++++++
       end if
c++++++++++++++++++++++++++++++++++++++++
c
c=======================================================================
         else
c=======================================================================
c
         if(ising.eq.k)then
           iikk=1
         else
           iikk=0
         end if
c++++++++++++++++++++++++++++++++++++++++
       if(itriqua(k).eq.3) then 
c++++++++++++++++++++++++++++++++++++++++
         call hemitri(xn(k),yn(k),zn(k),xni,yni,zni,
     $                x1(k),x2(k),x3(k),xh(k),yh(k),zh(k),
     $                y1(k),y2(k),y3(k),z1(k),z2(k),z3(k),
     $                cosxn(k),cosyn(k),coszn(k),u1,v1,w1,pig,iikk)
         iikk=0
         call hemitri(xn(k),-yn(k),zn(k),xni,yni,zni,
     $                x1(k),x2(k),x3(k),xh(k),-yh(k),zh(k),
     $                -y1(k),-y2(k),-y3(k),z1(k),z2(k),z3(k),
     $                cosxn(k),-cosyn(k),coszn(k),u2,v2,w2,pig,iikk)
         iikk=0
         call hemitri(xn(k),yn(k),-zn(k),xni,yni,zni,
     $                x1(k),x2(k),x3(k),xh(k),yh(k),-zh(k),
     $                y1(k),y2(k),y3(k),-z1(k),-z2(k),-z3(k),
     $                cosxn(k),cosyn(k),-coszn(k),u3,v3,w3,pig,iikk)
         iikk=0
         call hemitri(xn(k),-yn(k),-zn(k),xni,yni,zni,
     $                x1(k),x2(k),x3(k),xh(k),-yh(k),-zh(k),
     $                -y1(k),-y2(k),-y3(k),-z1(k),-z2(k),-z3(k),
     $                cosxn(k),-cosyn(k),-coszn(k),u4,v4,w4,pig,iikk)
         vx(k) = u1 + u2 + float(modo)*(u3+u4) 
         vy(k) = v1 + v2 + float(modo)*(v3+v4)
         vz(k) = w1 + w2 + float(modo)*(w3+w4)
c++++++++++++++++++++++++++++++++++++++++
       else
c++++++++++++++++++++++++++++++++++++++++
         call hemiqua(            xn(k),yn(k),zn(k),xni,yni,zni,
     $                      x1(k),x2(k),x3(k),x4(k),xh(k),yh(k),zh(k),
     $                y1(k),y2(k),y3(k),y4(k),z1(k),z2(k),z3(k),z4(k),
     $                  cosxn(k),cosyn(k),coszn(k),u1,v1,w1,pig,iikk)
         iikk=0
         call hemiqua(           xn(k),-yn(k),zn(k),xni,yni,zni,
     $                     x1(k),x2(k),x3(k),x4(k),xh(k),-yh(k),zh(k),
     $            -y1(k),-y2(k),-y3(k),-y4(k),z1(k),z2(k),z3(k),z4(k),
     $                 cosxn(k),-cosyn(k),coszn(k),u2,v2,w2,pig,iikk)
         iikk=0
         call hemiqua(           xn(k),yn(k),-zn(k),xni,yni,zni,
     $                     x1(k),x2(k),x3(k),x4(k),xh(k),yh(k),-zh(k),
     $            y1(k),y2(k),y3(k),y4(k),-z1(k),-z2(k),-z3(k),-z4(k),
     $                 cosxn(k),cosyn(k),-coszn(k),u3,v3,w3,pig,iikk)
         iikk=0
         call hemiqua(          xn(k),-yn(k),-zn(k),xni,yni,zni,
     $                    x1(k),x2(k),x3(k),x4(k),xh(k),-yh(k),-zh(k),
     $        -y1(k),-y2(k),-y3(k),-y4(k),-z1(k),-z2(k),-z3(k),-z4(k),
     $                cosxn(k),-cosyn(k),-coszn(k),u4,v4,w4,pig,iikk)
         vx(k) = u1 + u2 + float(modo)*(u3+u4) 
         vy(k) = v1 + v2 + float(modo)*(v3+v4)
         vz(k) = w1 + w2 + float(modo)*(w3+w4)
c++++++++++++++++++++++++++++++++++++++++
       endif
c++++++++++++++++++++++++++++++++++++++++
c
c=======================================================================
        end if
c=======================================================================
c
      end do
c$omp end do
c$omp end parallel

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
