c=======================================================================
c    prima
c=======================================================================
c
       subroutine primatri(x1,x2,x3,y1,y2,y3,z1,z2,z3,
     $      xh,yh,zh,area,w123,dilato,x12,y12,z12,x23,y23,
     $      z23,x13,y13,z13)
       implicit real (a-h,o-z)
c
       dist12  = sqrt ( (x1-x2)**2+(y1-y2)**2+(z1-z2)**2 )
       dist23  = sqrt ( (x2-x3)**2+(y2-y3)**2+(z2-z3)**2 )
       dist31  = sqrt ( (x3-x1)**2+(y3-y1)**2+(z3-z1)**2 )
       dilato  = max (dist12,dist23,dist31)
       altezza = sqrt ( (x1-xh)**2+(y1-yh)**2+(z1-zh)**2 ) 
c
       w123    = area / (2.0+1.0)
c
       x12 = 0.5 * ( x1 + x2 )
       y12 = 0.5 * ( y1 + y2 )
       z12 = 0.5 * ( z1 + z2 )
       x23 = 0.5 * ( x2 + x3 )
       y23 = 0.5 * ( y2 + y3 )
       z23 = 0.5 * ( z2 + z3 )
       x13 = 0.5 * ( x1 + x3 )
       y13 = 0.5 * ( y1 + y3 )
       z13 = 0.5 * ( z1 + z3 )
c
       return
       end
