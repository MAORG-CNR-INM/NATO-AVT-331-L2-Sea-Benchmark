c=======================================================================
c    JELL
c
c    Formule di quadratura di Lyness & Jespersen
c    INST. MATHS. APPLICS (1975) 15, 19-32
c
c=======================================================================
c
       subroutine jelltri(x12,y12,z12,x23,y23,z23,x13,y13,z13,
     $                    xp,yp,zp,w123,cosxl,cosyl,coszl,vll)
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
       a1     = sqrt ( xj1*xj1 + yj1*yj1 + zj1*zj1 )
       a2     = sqrt ( xj2*xj2 + yj2*yj2 + zj2*zj2 )
       a3     = sqrt ( xj3*xj3 + yj3*yj3 + zj3*zj3 )
       den1_3 = 1.0/(a1*a1*a1)
       den2_3 = 1.0/(a2*a2*a2)
       den3_3 = 1.0/(a3*a3*a3)
       den1_5 = 1.0/(a1*a1*a1*a1*a1)
       den2_5 = 1.0/(a2*a2*a2*a2*a2)
       den3_5 = 1.0/(a3*a3*a3*a3*a3)
c
       fpx1 = 3.0 * xj1 * xj1 * den1_5 - den1_3
       fpx2 = 3.0 * xj2 * xj2 * den2_5 - den2_3
       fpx3 = 3.0 * xj3 * xj3 * den3_5 - den3_3
       fpy1 = 3.0 * yj1 * yj1 * den1_5 - den1_3
       fpy2 = 3.0 * yj2 * yj2 * den2_5 - den2_3
       fpy3 = 3.0 * yj3 * yj3 * den3_5 - den3_3
       fpz1 = 3.0 * zj1 * zj1 * den1_5 - den1_3
       fpz2 = 3.0 * zj2 * zj2 * den2_5 - den2_3
       fpz3 = 3.0 * zj3 * zj3 * den3_5 - den3_3
       fxy1 = 3.0 * xj1 * yj1 * den1_5 
       fxy2 = 3.0 * xj2 * yj2 * den2_5
       fxy3 = 3.0 * xj3 * yj3 * den3_5
       fxz1 = 3.0 * xj1 * zj1 * den1_5 
       fxz2 = 3.0 * xj2 * zj2 * den2_5
       fxz3 = 3.0 * xj3 * zj3 * den3_5
       fyz1 = 3.0 * yj1 * zj1 * den1_5 
       fyz2 = 3.0 * yj2 * zj2 * den2_5
       fyz3 = 3.0 * yj3 * zj3 * den3_5
c
       u = - w123 * (fpx1+fpx2+fpx3)
       v = - w123 * (fpy1+fpy2+fpy3)
       w = - w123 * (fpz1+fpz2+fpz3)
       uv= - w123 * (fxy1+fxy2+fxy3)
       uw= - w123 * (fxz1+fxz2+fxz3)
       vw= - w123 * (fyz1+fyz2+fyz3)
c
       vll = cosxl**2 * u + 2.0 * cosxl * cosyl * uv + cosyl**2 * v +
     $       coszl**2 * w + 2.0 * cosxl * coszl * uw +
     $                      2.0 * cosyl * coszl * vw
c
       return
       end
