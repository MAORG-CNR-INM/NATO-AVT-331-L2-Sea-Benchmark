c=======================================================================
c    splyn
c
c    Formule di quadratura di Lyness & Jespersen
c    INST. MATHS. APPLICS (1975) 15, 19-32
c
c=======================================================================
c
       subroutine splyntri(x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $   xp,yp,zp,w123,u,v,w)
c
       implicit real (a-h,o-z)
c
       xj1 = x23 - xp
       yj1 = y23 - yp
       zj1 = z23 - zp
       xj2 = x12 - xp
       yj2 = y12 - yp
       zj2 = z12 - zp
       xj3 = x13 - xp
       yj3 = y13 - yp
       zj3 = z13 - zp
c
       a1   = sqrt ( xj1*xj1 + yj1*yj1 + zj1*zj1 )
       a2   = sqrt ( xj2*xj2 + yj2*yj2 + zj2*zj2 )
       a3   = sqrt ( xj3*xj3 + yj3*yj3 + zj3*zj3 )
       den1 = 1.0/(a1*a1*a1)
       den2 = 1.0/(a2*a2*a2)
       den3 = 1.0/(a3*a3*a3)
c
       fpx1 = xj1*den1
       fpx2 = xj2*den2
       fpx3 = xj3*den3
       fpy1 = yj1*den1
       fpy2 = yj2*den2
       fpy3 = yj3*den3
       fpz1 = zj1*den1
       fpz2 = zj2*den2
       fpz3 = zj3*den3
c
       u = - w123 * (fpx1+fpx2+fpx3)
       v = - w123 * (fpy1+fpy2+fpy3)
       w = - w123 * (fpz1+fpz2+fpz3)
c
c
       return
       end
