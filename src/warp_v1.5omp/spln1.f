c-----------------------------------------------------------
      subroutine spln1 (n,x,y,j,d,c,w)                                 
c     this subroutine fit a set of n points exactly with an
c     integrating function made of (n-1) cubics.       
c     the cubic between x(i-1) and x(i) must agree with the cubic      
c               between points x(i) and x(i+1) at the point x(i) in the
c               first and second derivatives.                          
c     spln1  determines the coefficiente of the cubics.               
c                                                                      
c     -----------------------------------------------------------------
c     argument definition                                              
c        n      is the number of data points.fortran integer (input).  
c                                                                      
c        x      is  a one-dimensional array of n independent variable  
c               points.the elements of x must be monotonically
c               increasing .floating-point array (input).
c                                                                      
c        y      is a one dimensional array of n dependent variable      
c               points corresponding to the x array.floating-point
c               array (input).                                               
c                                                                      
c        j      =1 if the first derivatives y:(x) are given;           
c               =2 if the second derivatives y@(x) are given;          
c               fortran integer (input).                               
c                                                                      
c        d      is an array of two elements:                           
c               d(1) is the jth derivative at x(1);                    
c               d(z) is the jth derivative at x(n).                    
c               floating-point array (input).                          
c                                                                      
c        c      is a temporary storage array of (3n-3) celles which hol
c               the coefficients of the interpolating cubics.          
c               floating-point array.                                  
c                                                                      
c        w      is a temporary storage array of (3n-3) celles.         
c               floatting-point array.                                 
c     -----------------------------------------------------------------
c     restrictions                                                     
c      there must be at least two data points.                         
c     -----------------------------------------------------------------
c     other subprograms required                                       
c        none.                                                         
c     -----------------------------------------------------------------
      implicit real (a-h,o-z)
      dimension x(n),y(n),d(2),c(3*n-3),w(3*n-3)                       
c     -----------------------------------------------------------------
c             over the interval x(i) to x(i+1), the interpolating      
c             polynomial                                               
c                  y=y(i)+a(i)*z+b(i)*z**2+e(i)*z**3                   
c             where z=(x-x(i))/(x(i+1)-x(i))                           
c             is used. the coefficients a(i),b(i) and e(i) are computed
c             by spln1 and stored in locations c(3*i-2),c(3*i-1) and   
c             c(3*i) respectively.                                     
c             while working in the ith interval,the variable q will    
c             represent q=x(i+1) - x(i), and y(i) will represent       
c             y(i+1)-y(i)                                              
c     -----------------------------------------------------------------
c                                                                      
      do i=1,3*n-3
        c(i) = 0.0
        w(i) = 0.0
      end do
      q=x(2) - x(1)                                                    
      yi =y(2) - y(1)                                                  
      if (j.eq.2) go to 100                                            
c     -----------------------------------------------------------------
c             if the first derivative at the end points is given,      
c             a(1) is known, and the second equation becomes           
c             merely b(1)+e(1)=yi - q*d(1).                            
c     -----------------------------------------------------------------
      c(1)=q*d(1)                                                      
      c(2)=1.0
      w(2)=yi-c(1)                                                     
      go to 200                                                        
c     -----------------------------------------------------------------
c             if the second derivative at the end points is given      
c             b(1) is known, the second equation becomes               
c             a(1)+e(1)=yi-0.5*q*q*d(1). during the solution of        
c             the 3n-4 equations,a1 will be kept in cell c(2)          
c             instead of c(1) to retain the tridiagonal form of the    
c             coefficient matrix.                                      
c     -----------------------------------------------------------------
100   c(2)=0.0
      w(2)=0.5*q*q*d(1)                                                
200   m=n-2                                                            
      if(m.le.0) go to 350                                             
c     -----------------------------------------------------------------
c             upper triangularization of the tridiagonal system of     
c             equations for the coefficient matrix follows--           
c     -----------------------------------------------------------------
      do 300 i=1,m                                                     
      ai=q                                                             
      q=x(i+2)- x(i+1)                                                 
      h=ai/q                                                           
      c(3*i)=-h/(2.0-c(3*i-1))                                         
      w(3*i)=(-yi-w(3*i-1))/(2.0 - c(3*i-1))                           
      c(3*i+1)=-h*h/(h-c(3*i))                                         
      w(3*i+1)=(yi-w(3*i))/(h-c(3*i))                                  
      yi=y(i+2)- y(i+1)                                                
      c(3*i+2)=1.0/(1.0-c(3*i+1))                                      
300   w(3*i+2)=(yi-w(3*i+1))/(1.0-c(3*i+1))                            
c     -----------------------------------------------------------------
c             e(n-1) is determined directly from the last equation     
c             obtained above, and the first or second derivative       
c             value given at the end point.                            
c     -----------------------------------------------------------------
350   if(j.eq.1) go to 400                                             
      c(3*n-3)=(q*q*d(2)/2.0-w(3*n-4))/(3.0- c(3*n-4))                 
      go to 500                                                        
400   c(3*n-3)=(q*d(2)-yi-w(3*n-4))/(2.0-c(3*n-4))                     
500   m=3*n-6                                                          
      if(m.le.0) go to 700                                             
c     -----------------------------------------------------------------
c             back solution for all coefficents except                 
c             a(1) and b(1) follows--                                  
c     -----------------------------------------------------------------
      do 600 ii=1,m                                                    
      i=m-ii+3                                                         
600   c(i)=w(i)-c(i)*c(i+1)                                            
700   if(j.eq.1) go to 800                                             
c     -----------------------------------------------------------------
c             if the second derivative is given at the end points,     
c             a(1) can now be computed from the known values of        
c             b(1) and e(1). then a(1) and b(1) are put into their     
c             proper places in the c array.                            
c     -----------------------------------------------------------------
      c(1)=y(2) - y(1)-w(2)-c(3)                                       
      c(2)=w(2)                                                        
      return                                                           
800   c(2)=w(2)-c(3)                                                   
      return                                                           
      end                                                              
