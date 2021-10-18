c-----------------------------------------------------------
      subroutine hemitrill(xn,yn,zn,xnt,ynt,znt,
     %  xv1,xv2,xv3,xh,yh,zh,yv1,yv2,yv3,zv1,
     %  zv2,zv3,coszitax,coszitay,coszitaz,
     %  vll,cosxl,cosyl,coszl)
c
      implicit real (a-h,o-z)
      real hsc(3),hss(3),hsdist(3),rpic(3),um(3),den(3)
c
c
c
c.......... RIFERIMENTO LOCALE ................................
c
      radx=sqrt((xv2-xv3)**2+(yv2-yv3)**2+(zv2-zv3)**2)
      coscsix=(xv2-xv3)/radx
      coscsiy=(yv2-yv3)/radx
      coscsiz=(zv2-zv3)/radx  
      rady=sqrt((xv1-xh)**2+(yv1-yh)**2+(zv1-zh)**2)
      cosetax=(xv1-xh)/rady
      cosetay=(yv1-yh)/rady
      cosetaz=(zv1-zh)/rady
      xx=(xnt-xn)*coscsix+(ynt-yn)*coscsiy+(znt-zn)*coscsiz
      yy=(xnt-xn)*cosetax+(ynt-yn)*cosetay+(znt-zn)*cosetaz
      zz=(xnt-xn)*coszitax+(ynt-yn)*coszitay+(znt-zn)*coszitaz
      x1=(xv1-xn)*coscsix+(yv1-yn)*coscsiy+(zv1-zn)*coscsiz
      x2=(xv2-xn)*coscsix+(yv2-yn)*coscsiy+(zv2-zn)*coscsiz
      x3=(xv3-xn)*coscsix+(yv3-yn)*coscsiy+(zv3-zn)*coscsiz
      y1=(xv1-xn)*cosetax+(yv1-yn)*cosetay+(zv1-zn)*cosetaz
      y2=(xv2-xn)*cosetax+(yv2-yn)*cosetay+(zv2-zn)*cosetaz
      y3=(xv3-xn)*cosetax+(yv3-yn)*cosetay+(zv3-zn)*cosetaz
c
c..............................................
      hsdist(1)=sqrt((x1-x2)**2+(y1-y2)**2)
      hsdist(2)=sqrt((x2-x3)**2+(y2-y3)**2)
      hsdist(3)=sqrt((x3-x1)**2+(y3-y1)**2)
      hsc(1)=(x2-x1) 
      hsc(2)=(x3-x2)	 
      hsc(3)=(x1-x3)	 
      hss(1)=(y2-y1)
      hss(2)=(y3-y2)	 
      hss(3)=(y1-y3)	 
      do i = 1,3
        hsc(i)=hsc(i)/hsdist(i)
        hss(i)=hss(i)/hsdist(i)
      end do
c	 
c---------------------------------------------------
c
      rpic(1)=sqrt((xx-x1)**2+(yy-y1)**2+zz**2)
      rpic(2)=sqrt((xx-x2)**2+(yy-y2)**2+zz**2)
      rpic(3)=sqrt((xx-x3)**2+(yy-y3)**2+zz**2)
      um(1) =rpic(1)+rpic(2)+hsdist(1)
      den(1)=rpic(1)+rpic(2)-hsdist(1)
      um(2) =rpic(2)+rpic(3)+hsdist(2)
      den(2)=rpic(2)+rpic(3)-hsdist(2)
      um(3) =rpic(3)+rpic(1)+hsdist(3)
      den(3)=rpic(3)+rpic(1)-hsdist(3)
c
c
c...................... DERIVATE SECONDE .........................
c
      hsqx1=1./(rpic(1)+rpic(2)+hsdist(1))*((xx-x1)/rpic(1)+
     %      (xx-x2)/rpic(2))-1./(rpic(1)+rpic(2)-hsdist(1))*
     %      ((xx-x1)/rpic(1)+(xx-x2)/rpic(2))
      hsqx2=1./(rpic(2)+rpic(3)+hsdist(2))*((xx-x2)/rpic(2)+
     %      (xx-x3)/rpic(3))-1./(rpic(2)+rpic(3)-hsdist(2))*
     %      ((xx-x2)/rpic(2)+(xx-x3)/rpic(3))
      hsqx3=1./(rpic(3)+rpic(1)+hsdist(3))*((xx-x3)/rpic(3)+
     %      (xx-x1)/rpic(1))-1./(rpic(3)+rpic(1)-hsdist(3))*
     %      ((xx-x3)/rpic(3)+(xx-x1)/rpic(1))
      hsqy1=1./(rpic(1)+rpic(2)+hsdist(1))*((yy-y1)/rpic(1)+
     %      (yy-y2)/rpic(2))-1./(rpic(1)+rpic(2)-hsdist(1))*
     %      ((yy-y1)/rpic(1)+(yy-y2)/rpic(2))
      hsqy2=1./(rpic(2)+rpic(3)+hsdist(2))*((yy-y2)/rpic(2)+
     %      (yy-y3)/rpic(3))-1./(rpic(2)+rpic(3)-hsdist(2))*
     %      ((yy-y2)/rpic(2)+(yy-y3)/rpic(3))
      hsqy3=1./(rpic(3)+rpic(1)+hsdist(3))*((yy-y3)/rpic(3)+
     %      (yy-y1)/rpic(1))-1./(rpic(3)+rpic(1)-hsdist(3))*
     %      ((yy-y3)/rpic(3)+(yy-y1)/rpic(1))
      hsqz1=zz*(1./um(1)-1./den(1))*(1./rpic(1)+1./rpic(2))
      hsqz2=zz*(1./um(2)-1./den(2))*(1./rpic(2)+1./rpic(3))
      hsqz3=zz*(1./um(3)-1./den(3))*(1./rpic(3)+1./rpic(1))
c
c
c
      vxxloc = - (hss(1)*hsqx1+hss(2)*hsqx2+hss(3)*hsqx3)
c
      vyyloc = + (hsc(1)*hsqy1+hsc(2)*hsqy2+hsc(3)*hsqy3)
c
      vzzloc = - (vxxloc+vyyloc)
c
      vxyloc = - (hss(1)*hsqy1+hss(2)*hsqy2+hss(3)*hsqy3)
c
      vxzloc = - (hss(1)*hsqz1+hss(2)*hsqz2+hss(3)*hsqz3)
c
      vyzloc = + (hsc(1)*hsqz1+hsc(2)*hsqz2+hsc(3)*hsqz3)
c
c
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
