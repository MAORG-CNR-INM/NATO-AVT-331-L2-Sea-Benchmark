c-----------------------------------------------------------
      subroutine spln2 (n,x,y,c,v)                                     
c     this subroutine fit a set of n points exactly with an interpolati
c     g function  made up of (n-1) cubics.                             
c        the cubic between x(i)-1 and x(i) must agree with the cubic   
c     between points x(i) and x(i)+(1) at the point x(i) in their first
c     and second derivatives                                           
c                                                                      
c                                                                      
c     -----------------------------------------------------------------
c     argument definition                                              
c        n      it the number of data points.fortran integer (input).  
c                                                                      
c        x      is  a one-dimensional array of n independent variable  
c               points.the elements of x must be monotonically increasi
c               g.floating-point array (input).                        
c                                                                      
c        y      is a one dimensional array of n dependent varianle     
c               points corresponding to the x array.floating-point arra
c               (input).                                               
c        c      is a temporary storage array of (3n-3) celles.         
c                                                                      
c        v      is an array of five elemenets:                         
c               v(1)= point z at which to interpolate;                 
c               v(2)= value of y(z);                                   
c               v(3)= value of y'(z);                                  
c               v(4)= value of y"(z);                                  
c               v(5)= 1,o for  z < x(1)                                  
c                     3.o for  z > x(n)                                 
c                     2.o for  x(1) < x < x(n).                         
c               floating-point array (input/output).                   
c     -----------------------------------------------------------------
c     -------   restrictions                                           
c      there must be at least to data points.                          
c     -----------------------------------------------------------------
c      other subprograms required.                                     
c       none.                                                          
c     -----------------------------------------------------------------
      implicit real (a-h,o-z)
      dimension x(n),y(n),c(3*n-3),v(5)                                
      v(5)=2.0
      lim=n-1                                                          
c                  ----------------------------------------------------
c                  determine in which interval the independent         
c                  variable,v(1),lies.                                 
c                  ----------------------------------------------------
      do 10 i=2,lim                                                    
10    if(v(1).lt.x(i)) go to 20                                        
      i=n                                                              
      if(v(1).gt.x(n)) v(5)=3.0
      go to 30                                                         
20    if(v(1).lt.x(1)) v(5) =1.0
c                  ----------------------------------------------------
c                  q is the size of the interval containing v(1).      
c                  ----------------------------------------------------
c                  z is a linear transformation of the interval        
c                  onto (0,1) and is the variable for which            
c                  the coefficients were computed by spln1.            
c                  ----------------------------------------------------
   30 q=x(i)-x(i-1)                                                    
      z=(v(1)-x(i-1))/q                                                
      v(2)=((z*c(3*i-3)+c(3*i-4))*z+c(3*i-5))*z+y(i-1)                 
      v(3)=((3.0*z*c(3*i-3)+2.0*c(3*i-4))*z+c(3*i-5))/q                 
      v(4)=(6.0*z*c(3*i-3)+2.0*c(3*i-4))/(q*q)                          
      return                                                           
      end                                                              
