c CALCOLA IL DISLOCAMENTO DELLA CARENA
c-----------------------------------------------------------------------
      subroutine DISLO(x1,x2,x3,x4,y1,y2,y3,y4,yn,z1,z2,z3,z4,shiplen,
     &                 disloca,n_est,n_int)
c
      dimension x1(*),x2(*),x3(*),x4(*),y1(*),y2(*),y3(*),y4(*),
     &          z1(*),z2(*),z3(*),z4(*),yn(*)
c-----------------------------------------------------------------------
      disloca = 0.
      do i=1,n_est
         v12x = (x2(i)-x1(i))*shiplen
         v12z = (z2(i)-z1(i))*shiplen
         v14x = (x4(i)-x1(i))*shiplen
         v14z = (z4(i)-z1(i))*shiplen
         area = sqrt(v12x**2+v12z**2) *
     &          sqrt(v14x**2+v14z**2)
         volmin = area*yn(i)*shiplen
         disloca= disloca + volmin
      end do
      do i=n_est+1,n_est+n_int
         v12x = (x2(i)-x1(i))*shiplen
         v12z = (z2(i)-z1(i))*shiplen
         v14x = (x4(i)-x1(i))*shiplen
         v14z = (z4(i)-z1(i))*shiplen
         area = sqrt(v12x**2+v12z**2) *
     &          sqrt(v14x**2+v14z**2)
         volmin = area*yn(i)*shiplen
         disloca= disloca - volmin
      end do
      disloca = disloca * 2.
c-----------------------------------------------------------------------
      return
      end
