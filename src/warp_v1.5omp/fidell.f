c
        subroutine fidell(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
     $             xh,yh,zh,xn,yn,zn,area,cosxn,cosyn,coszn,sigma0,
     $             xns,yns,zns,
     $             fi0x,fi0y,fi0z,diana,fi0ll,vll,
     $             ilib,ncar,nekel,itriqua,modo,ncard,ntotd,nplib)

c-----------------------------------------------------------------------

      integer ncard,nplib,ntotd
      integer itriqua(ntotd)

      real x1(ntotd),y1(ntotd),z1(ntotd),
     &     x2(ntotd),y2(ntotd),z2(ntotd),
     &     x3(ntotd),y3(ntotd),z3(ntotd),
     &     x4(ntotd),y4(ntotd),z4(ntotd)
      real xh(ntotd),yh(ntotd),zh(ntotd)
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real area(ntotd)
      real cosxn(ntotd),cosyn(ntotd),coszn(ntotd)
      real sigma0(ntotd)
      real xns(ntotd),yns(ntotd),zns(ntotd)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd)
      real vll(ntotd)

c-----------------------------------------------------------------------

      integer k,modo,ilib,nekel,ncar
      real x12,x13,x23,x34,x41
      real y12,y13,y23,y34,y41
      real z12,z13,z23,z34,z41
      real w123,w341
      real cosxl,cosyl,coszl,dilato,drelpq
      real fil,diana,fi0ll,vll1,vll2

c----------------------------------------------------------------------*
c    Pannello collocante p: xns(i), yns(i), zns(i)
c    la linea di corrente passante per il pann. coll. e' diversa nelle
c    differenti linearizzazioni
c----------------------------------------------------------------------*
        if(nekel.eq.0) then                      !!! dawson
           fil   = sqrt((1.0+fi0x(ilib))**2+fi0y(ilib)**2+fi0z(ilib)**2)
           cosxl = (1.0+fi0x(ilib))/fil
           cosyl = fi0y(ilib)/fil
           coszl = fi0z(ilib)/fil
        else if(nekel.eq.1) then                 !!! kelvin
           cosxl = 1.0 
           cosyl = 0.0
           coszl = 0.0
        end if

c$omp parallel
c$omp do private(k,w123,w341,dilato,drelpq,vll1,vll2,
c$omp&           x12,y12,z12,
c$omp&           x23,y23,z23,
c$omp&           x13,y13,z13,
c$omp&           x34,y34,z34,
c$omp&           x41,y41,z41,
c$omp&           u1,u2,u3,u4,
c$omp&           v1,v2,v3,v4,
c$omp&           w1,w2,w3,w4)

        do 100 k=1,ncar+nplib
c++++++++++++++++++++++++++++++++++++++++
        if(itriqua(k).eq.3) then
c++++++++++++++++++++++++++++++++++++++++
           call primatri(x1(k),x2(k),x3(k),y1(k),y2(k),y3(k),
     $            z1(k),z2(k),z3(k),xh(k),yh(k),zh(k),area(k),
     $         w123,dilato,x12,y12,z12,x23,y23,z23,x13,y13,z13)
c++++++++++++++++++++++++++++++++++++++++
        else
c++++++++++++++++++++++++++++++++++++++++
           call primaqua(x1(k),x2(k),x3(k),x4(k),
     $                   y1(k),y2(k),y3(k),y4(k),
     $            z1(k),z2(k),z3(k),z4(k),xh(k),yh(k),zh(k),area(k),
     $         w123,w341,dilato,x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $                                      x34,y34,z34,x41,y41,z41)
c++++++++++++++++++++++++++++++++++++++++
        end if
c++++++++++++++++++++++++++++++++++++++++
c
        drelpq = sqrt( (xns(ilib)-xn(k))**2 + ( yns(ilib)-yn(k) )**2
     $               + (zns(ilib)-zn(k))**2 ) / dilato
c
c=======================================================================
        if (drelpq.ge.diana) then
c=======================================================================
c++++++++++++++++++++++++++++++++++++++++
        if(itriqua(k).eq.3) then
c++++++++++++++++++++++++++++++++++++++++
        call jelltri(x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $               xns(ilib),yns(ilib),zns(ilib),w123,
     $               cosxl,cosyl,coszl,vll1)
        call jelltri(x12,-y12,z12,x23,-y23,z23,x13,-y13,z13,
     $               xns(ilib),yns(ilib),zns(ilib),w123,
     $               cosxl,cosyl,coszl,vll2)
c++++++++++++++++++++++++++++++++++++++++
        else
c++++++++++++++++++++++++++++++++++++++++
        call jellqua(x12,y12,z12,x23,y23,z23,x13,y13,z13,x34,y34,z34,
     $           x41,y41,z41,xns(ilib),yns(ilib),zns(ilib),w123,w341,
     $           cosxl,cosyl,coszl,vll1)
        call jellqua(x12,-y12,z12,x23,-y23,z23,x13,-y13,z13,
     $           x34,-y34,z34,
     $           x41,-y41,z41,xns(ilib),yns(ilib),zns(ilib),w123,w341,
     $           cosxl,cosyl,coszl,vll2)
c++++++++++++++++++++++++++++++++++++++++
        end if
c++++++++++++++++++++++++++++++++++++++++
c
c=======================================================================
        else
c=======================================================================
c++++++++++++++++++++++++++++++++++++++++
        if(itriqua(k).eq.3) then
c++++++++++++++++++++++++++++++++++++++++
        call hemitrill(xn(k),yn(k),zn(k),xns(ilib),yns(ilib),zns(ilib),
     $           x1(k),x2(k),x3(k),xh(k),yh(k),zh(k),y1(k),y2(k),y3(k),
     $           z1(k),z2(k),z3(k),cosxn(k),cosyn(k),coszn(k),
     $           vll1,cosxl,cosyl,coszl)
       call hemitrill(xn(k),-yn(k),zn(k),xns(ilib),yns(ilib),zns(ilib),
     $       x1(k),x2(k),x3(k),xh(k),-yh(k),zh(k),-y1(k),-y2(k),-y3(k),
     $       z1(k),z2(k),z3(k),cosxn(k),-cosyn(k),coszn(k),
     $       vll2,cosxl,cosyl,coszl)
c++++++++++++++++++++++++++++++++++++++++
        else
c++++++++++++++++++++++++++++++++++++++++
        call hemiquall(xn(k),yn(k),zn(k),xns(ilib),yns(ilib),zns(ilib),
     $           x1(k),x2(k),
     $           x3(k),x4(k),xh(k),yh(k),zh(k),y1(k),y2(k),y3(k),y4(k),
     $              z1(k),z2(k),z3(k),z4(k),cosxn(k),cosyn(k),coszn(k),
     $                                          vll1,cosxl,cosyl,coszl)
       call hemiquall(xn(k),-yn(k),zn(k),xns(ilib),yns(ilib),zns(ilib),
     $      x1(k),x2(k),
     $      x3(k),x4(k),xh(k),-yh(k),zh(k),-y1(k),-y2(k),-y3(k),-y4(k),
     $             z1(k),z2(k),z3(k),z4(k),cosxn(k),-cosyn(k),coszn(k),
     $                                          vll2,cosxl,cosyl,coszl)
c++++++++++++++++++++++++++++++++++++++++
       endif
c++++++++++++++++++++++++++++++++++++++++
c
c=======================================================================
        end if
c=======================================================================
c
        vll(k) = vll1 + vll2

  100   continue 
c$omp end do
c$omp end parallel

c
        fi0ll = 0.0
        do k=1,ncar
          fi0ll = fi0ll + float(1+modo) * vll(k) * sigma0(k)
        end do
c
c
c ----------------------------------------------------------------------
        return
        end 
c ----------------------------------------------------------------------
