
c=======================================================================
c    prima
c=======================================================================
c
       subroutine primaqua(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
     $      xh,yh,zh,area,w123,w341,dilato,x12,y12,z12,x23,y23,
     $      z23,x13,y13,z13,x34,y34,z34,x41,y41,z41)
       implicit real (a-h,o-z)
c
       dist12  = sqrt ( (x1-x2)**2+(y1-y2)**2+(z1-z2)**2 )
       dist23  = sqrt ( (x2-x3)**2+(y2-y3)**2+(z2-z3)**2 )
       dist34  = sqrt ( (x3-x4)**2+(y3-y4)**2+(z3-z4)**2 )
       dist41  = sqrt ( (x4-x1)**2+(y4-y1)**2+(z4-z1)**2 )
       dilato  = max (dist12,dist23,dist34,dist41)
       altezza = sqrt ( (x1-xh)**2+(y1-yh)**2+(z1-zh)**2 ) 
c
       area341 = dist34 * altezza * 0.5
       area123 = area - area341
       w123    = area123 / (2.0+1.0)
       w341    = area341 / (2.0+1.0)
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
       x34 = 0.5 * ( x3 + x4 )
       y34 = 0.5 * ( y3 + y4 )
       z34 = 0.5 * ( z3 + z4 )
       x41 = 0.5 * ( x4 + x1 )
       y41 = 0.5 * ( y4 + y1 )
       z41 = 0.5 * ( z4 + z1 )
c
       return
       end
