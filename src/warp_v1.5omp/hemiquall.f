c-----------------------------------------------------------
      subroutine hemiquall(xn,yn,zn,xnt,ynt,znt,
     &                     xv1,xv2,xv3,xv4,xh,yh,zh,yv1,yv2,yv3,yv4,
     &                     zv1,zv2,zv3,zv4,
     &                     coszitax,coszitay,coszitaz,
     &                     vll,cosxl,cosyl,coszl)
c-----------------------------------------------------------
      integer i
      real hsc(4),hss(4),hsdist(4),rpic(4),um(4),den(4)
      real x1,x2,x3,x4,y1,y2,y3,y4
      real radx,coscsix,coscsiy,coscsiz,rady,cosetax,cosetay,cosetaz
      real xx,yy,zz
      real xn,yn,zn,xnt,ynt,znt,xv1,xv2,xv3,xv4,
     &     xh,yh,zh,yv1,yv2,yv3,yv4,zv1,zv2,zv3,zv4,
     &     coszitax,coszitay,coszitaz,cosxl,cosyl,coszl,
     &     vvxx,vvyy,vvzz,vvxy,vvxz,vvyz,vll,
     &     vxyloc,vxzloc,vyzloc,
     &     vxxloc,vyyloc,vzzloc
      real hsqx1,hsqx2,hsqx3,hsqx4
      real hsqy1,hsqy2,hsqy3,hsqy4
      real hsqz1,hsqz2,hsqz3,hsqz4
c-----------------------------------------------------------
c
c.......... RIFERIMENTO LOCALE ................................
c
c
      c1 = xv3-xv4
      c2 = yv3-yv4
      c3 = zv3-zv4

      radx=sqrt(c1**2+c2**2+c3**2)
      coscsix=c1/radx
      coscsiy=c2/radx
      coscsiz=c3/radx  

      c1 = xv1-xh
      c2 = yv1-yh
      c3 = zv1-zh
      rady=sqrt(c1**2+c2**2+c3**2)
      cosetax=c1/rady
      cosetay=c2/rady
      cosetaz=c3/rady  
c
      xntxn = xnt-xn
      yntyn = ynt-yn
      zntzn = znt-zn

      xx=xntxn*coscsix+yntyn*coscsiy+zntzn*coscsiz
      yy=xntxn*cosetax+yntyn*cosetay+zntzn*cosetaz
      zz=xntxn*coszitax+yntyn*coszitay+zntzn*coszitaz
c
      xv1xn = xv1-xn
      yv1yn = yv1-yn
      zv1zn = zv1-zn
      xv2xn = xv2-xn
      yv2yn = yv2-yn
      zv2zn = zv2-zn
      xv3xn = xv3-xn
      yv3yn = yv3-yn
      zv3zn = zv3-zn
      xv4xn = xv4-xn
      yv4yn = yv4-yn
      zv4zn = zv4-zn

      x1=xv1xn*coscsix+yv1yn*coscsiy+zv1zn*coscsiz
      x2=xv2xn*coscsix+yv2yn*coscsiy+zv2zn*coscsiz
      x3=xv3xn*coscsix+yv3yn*coscsiy+zv3zn*coscsiz
      x4=xv4xn*coscsix+yv4yn*coscsiy+zv4zn*coscsiz
      y1=xv1xn*cosetax+yv1yn*cosetay+zv1zn*cosetaz
      y2=xv2xn*cosetax+yv2yn*cosetay+zv2zn*cosetaz
      y3=xv3xn*cosetax+yv3yn*cosetay+zv3zn*cosetaz
      y4=xv4xn*cosetax+yv4yn*cosetay+zv4zn*cosetaz
c
      x2x1 = x2-x1
      x3x2 = x3-x2
      x4x3 = x4-x3
      x1x4 = x1-x4
      y2y1 = y2-y1
      y3y2 = y3-y2
      y4y3 = y4-y3
      y1y4 = y1-y4

      hsdist(1)=sqrt(x2x1**2+y2y1**2)
      hsdist(2)=sqrt(x3x2**2+y3y2**2)
      hsdist(3)=sqrt(x4x3**2+y4y3**2)
      hsdist(4)=sqrt(x1x4**2+y1y4**2)

      hsc(1)=x2x1
      hsc(2)=x3x2
      hsc(3)=x4x3 
      hsc(4)=x1x4
      hss(1)=y2y1
      hss(2)=y3y2 
      hss(3)=y4y3
      hss(4)=y1y4

      do i = 1,4
        hsc(i)=hsc(i)/hsdist(i)
        hss(i)=hss(i)/hsdist(i)
      end do
c
      x1xx = x1-xx
      x2xx = x2-xx
      x3xx = x3-xx
      x4xx = x4-xx

      y1yy = y1-yy
      y2yy = y2-yy
      y3yy = y3-yy
      y4yy = y4-yy

      rpic(1)=sqrt(x1xx**2+y1yy**2+zz**2)
      rpic(2)=sqrt(x2xx**2+y2yy**2+zz**2)
      rpic(3)=sqrt(x3xx**2+y3yy**2+zz**2)
      rpic(4)=sqrt(x4xx**2+y4yy**2+zz**2)
      um(1) =rpic(1)+rpic(2)+hsdist(1)
      den(1)=rpic(1)+rpic(2)-hsdist(1)
      um(2) =rpic(2)+rpic(3)+hsdist(2)
      den(2)=rpic(2)+rpic(3)-hsdist(2)
      um(3) =rpic(3)+rpic(4)+hsdist(3)
      den(3)=rpic(3)+rpic(4)-hsdist(3)
      um(4) =rpic(4)+rpic(1)+hsdist(4)
      den(4)=rpic(4)+rpic(1)-hsdist(4)
c
c...................... DERIVATE SECONDE .........................
c
      hsqx1=1./(rpic(1)+rpic(2)+hsdist(1))*((xx-x1)/rpic(1)+
     %      (xx-x2)/rpic(2))-1./(rpic(1)+rpic(2)-hsdist(1))*
     %      ((xx-x1)/rpic(1)+(xx-x2)/rpic(2))
      hsqx2=1./(rpic(2)+rpic(3)+hsdist(2))*((xx-x2)/rpic(2)+
     %      (xx-x3)/rpic(3))-1./(rpic(2)+rpic(3)-hsdist(2))*
     %      ((xx-x2)/rpic(2)+(xx-x3)/rpic(3))
      hsqx3=1./(rpic(3)+rpic(4)+hsdist(3))*((xx-x3)/rpic(3)+
     %      (xx-x4)/rpic(4))-1./(rpic(3)+rpic(4)-hsdist(3))*
     %      ((xx-x3)/rpic(3)+(xx-x4)/rpic(4))
      hsqx4=1./(rpic(4)+rpic(1)+hsdist(4))*((xx-x1)/rpic(1)+
     %      (xx-x4)/rpic(4))-1./(rpic(1)+rpic(4)-hsdist(4))*
     %      ((xx-x1)/rpic(1)+(xx-x4)/rpic(4))
      hsqy1=1./(rpic(1)+rpic(2)+hsdist(1))*((yy-y1)/rpic(1)+
     %      (yy-y2)/rpic(2))-1./(rpic(1)+rpic(2)-hsdist(1))*
     %      ((yy-y1)/rpic(1)+(yy-y2)/rpic(2))
      hsqy2=1./(rpic(2)+rpic(3)+hsdist(2))*((yy-y2)/rpic(2)+
     %      (yy-y3)/rpic(3))-1./(rpic(2)+rpic(3)-hsdist(2))*
     %      ((yy-y2)/rpic(2)+(yy-y3)/rpic(3))
      hsqy3=1./(rpic(3)+rpic(4)+hsdist(3))*((yy-y3)/rpic(3)+
     %      (yy-y4)/rpic(4))-1./(rpic(3)+rpic(4)-hsdist(3))*
     %      ((yy-y3)/rpic(3)+(yy-y4)/rpic(4))
      hsqy4=1./(rpic(4)+rpic(1)+hsdist(4))*((yy-y1)/rpic(1)+
     %      (yy-y4)/rpic(4))-1./(rpic(1)+rpic(4)-hsdist(4))*
     %      ((yy-y1)/rpic(1)+(yy-y4)/rpic(4))
      hsqz1=zz*(1./um(1)-1./den(1))*(1./rpic(1)+1./rpic(2))
      hsqz2=zz*(1./um(2)-1./den(2))*(1./rpic(2)+1./rpic(3))
      hsqz3=zz*(1./um(3)-1./den(3))*(1./rpic(3)+1./rpic(4))
      hsqz4=zz*(1./um(4)-1./den(4))*(1./rpic(4)+1./rpic(1))

      vxxloc = - (hss(1)*hsqx1+hss(2)*hsqx2+hss(3)*hsqx3+hss(4)*hsqx4)

      vyyloc = + (hsc(1)*hsqy1+hsc(2)*hsqy2+hsc(3)*hsqy3+hsc(4)*hsqy4)

      vzzloc = - (vxxloc+vyyloc)

      vxyloc = - (hss(1)*hsqy1+hss(2)*hsqy2+hss(3)*hsqy3+hss(4)*hsqy4)

      vxzloc = - (hss(1)*hsqz1+hss(2)*hsqz2+hss(3)*hsqz3+hss(4)*hsqz4)

      vyzloc = + (hsc(1)*hsqz1+hsc(2)*hsqz2+hsc(3)*hsqz3+hsc(4)*hsqz4)

      vvxx = vxxloc*coscsix**2 + vyyloc*cosetax**2 + vzzloc*coszitax**2
     % + 2.*vxyloc*coscsix*cosetax + 2.*vxzloc*coscsix*coszitax +
     % 2.*vyzloc*cosetax*coszitax
c
      vvyy = vxxloc*coscsiy**2 + vyyloc*cosetay**2 + vzzloc*coszitay**2
     % + 2.*vxyloc*coscsiy*cosetay + 2.*vxzloc*coscsiy*coszitay +
     % 2.*vyzloc*cosetay*coszitay
c
      vvzz = - (vvxx+vvyy)
c
      vvxy = vxxloc*coscsix*coscsiy + vyyloc*cosetax*cosetay + 
     % vzzloc*coszitax*coszitay + 
     % vxyloc*(coscsix*cosetay+coscsiy*cosetax) + 
     % vxzloc*(coscsix*coszitay+coszitax*coscsiy) +
     % vyzloc*(cosetax*coszitay+coszitax*cosetay)
c
      vvxz = vxxloc*coscsix*coscsiz + vyyloc*cosetax*cosetaz + 
     % vzzloc*coszitax*coszitaz + 
     % vxyloc*(coscsix*cosetaz+coscsiz*cosetax) + 
     % vxzloc*(coscsix*coszitaz+coszitax*coscsiz) +
     % vyzloc*(cosetax*coszitaz+coszitax*cosetaz)
c
      vvyz = vxxloc*coscsiy*coscsiz + vyyloc*cosetay*cosetaz + 
     % vzzloc*coszitay*coszitaz + 
     % vxyloc*(coscsiy*cosetaz+coscsiz*cosetay) + 
     % vxzloc*(coscsiy*coszitaz+coszitay*coscsiz) +
     % vyzloc*(cosetay*coszitaz+coszitay*cosetaz)
c
c
      vll = cosxl**2*vvxx + cosyl**2*vvyy + coszl**2*vvzz + 
     &       2.*cosxl*cosyl*vvxy + 2.*cosxl*coszl*vvxz + 
     &       2.*cosyl*coszl*vvyz 
c
      return 
      end
