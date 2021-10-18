c------------------------------------------------------------
      subroutine hemitri(xn,yn,zn,xnt,ynt,znt,xv1,xv2,xv3,
     $ xh,yh,zh,yv1,yv2,yv3,zv1,zv2,zv3,coszitax,
     $ coszitay,coszitaz,vvx,vvy,vvz,pig,iikk)
c
      implicit real (a-h,o-z)
      real hsc(3),hss(3),hsdist(3),spicuno(3),spicdue(3),
     &     hsq(3),hsj(3),rpic(3),rr(3),arg1(3),arg2(3)
c
c
      radx=sqrt((xv2-xv3)**2+(yv2-yv3)**2+(zv2-zv3)**2)
      coscsix=(xv2-xv3)/radx
      coscsiy=(yv2-yv3)/radx
      coscsiz=(zv2-zv3)/radx  
      rady=sqrt((xv1-xh)**2+(yv1-yh)**2+(zv1-zh)**2)
      cosetax=(xv1-xh)/rady
      cosetay=(yv1-yh)/rady
      cosetaz=(zv1-zh)/rady
c
      xx=(xnt-xn)*coscsix+(ynt-yn)*coscsiy+(znt-zn)*coscsiz
      yy=(xnt-xn)*cosetax+(ynt-yn)*cosetay+(znt-zn)*cosetaz
      zz=(xnt-xn)*coszitax+(ynt-yn)*coszitay+(znt-zn)*coszitaz
c
      x1=(xv1-xn)*coscsix+(yv1-yn)*coscsiy+(zv1-zn)*coscsiz
      x2=(xv2-xn)*coscsix+(yv2-yn)*coscsiy+(zv2-zn)*coscsiz
      x3=(xv3-xn)*coscsix+(yv3-yn)*coscsiy+(zv3-zn)*coscsiz
      y1=(xv1-xn)*cosetax+(yv1-yn)*cosetay+(zv1-zn)*cosetaz
      y2=(xv2-xn)*cosetax+(yv2-yn)*cosetay+(zv2-zn)*cosetaz
      y3=(xv3-xn)*cosetax+(yv3-yn)*cosetay+(zv3-zn)*cosetaz
c
      hsdist(1)=sqrt((x1-x2)**2+(y1-y2)**2)
      hsdist(2)=sqrt((x2-x3)**2+(y2-y3)**2)
      hsdist(3)=sqrt((x3-x1)**2+(y3-y1)**2)
      hsc(1)=(x2-x1) 
      hsc(2)=(x3-x2)	 
      hsc(3)=(x1-x3)	 
      hss(1)=(y2-y1)
      hss(2)=(y3-y2)	 
      hss(3)=(y1-y3)	 
c
      do i = 1,3
        hsc(i)=hsc(i)/hsdist(i)
        hss(i)=hss(i)/hsdist(i)
      end do
c
      spicuno(1)=(x1-xx)*hsc(1)+(y1-yy)*hss(1)	   
      spicdue(1)=(x2-xx)*hsc(1)+(y2-yy)*hss(1)	   
      spicuno(2)=(x2-xx)*hsc(2)+(y2-yy)*hss(2)	
      spicdue(2)=(x3-xx)*hsc(2)+(y3-yy)*hss(2)	
      spicuno(3)=(x3-xx)*hsc(3)+(y3-yy)*hss(3)	
      spicdue(3)=(x1-xx)*hsc(3)+(y1-yy)*hss(3)	
c
      rr(1)=(xx-x1)*hss(1)-(yy-y1)*hsc(1)       
      rr(2)=(xx-x2)*hss(2)-(yy-y2)*hsc(2)
      rr(3)=(xx-x3)*hss(3)-(yy-y3)*hsc(3)
c--------------------------------------------------
      deteta=2.0*pig
      do ii=1,3      
        if(rr(ii).lt.0.0) then
          deteta=0.0
          go to 9191
	end if
      end do
 9191 continue	
c---------------------------------------------------
      rpic(1)=sqrt((xx-x1)**2+(yy-y1)**2+zz**2)
      rpic(2)=sqrt((xx-x2)**2+(yy-y2)**2+zz**2)
      rpic(3)=sqrt((xx-x3)**2+(yy-y3)**2+zz**2)
c
      hsq(1)=log((rpic(1)+rpic(2)+hsdist(1))/(rpic(1)+
     $             rpic(2)-hsdist(1)))
      hsq(2)=log((rpic(2)+rpic(3)+hsdist(2))/(rpic(2)+
     $             rpic(3)-hsdist(2)))
      hsq(3)=log((rpic(3)+rpic(1)+hsdist(3))/(rpic(3)+
     $             rpic(1)-hsdist(3)))
      arg1(1)=rr(1)*abs(zz)*(rpic(1)*spicdue(1)-rpic(2)*spicuno(1))
      arg2(1)=rpic(1)*rpic(2)*(rr(1)**2)+(zz**2)*spicuno(1)*spicdue(1)
      arg1(2)=rr(2)*abs(zz)*(rpic(2)*spicdue(2)-rpic(3)*spicuno(2))
      arg2(2)=rpic(2)*rpic(3)*(rr(2)**2)+(zz**2)*spicuno(2)*spicdue(2)
      arg1(3)=rr(3)*abs(zz)*(rpic(3)*spicdue(3)-rpic(1)*spicuno(3))
      arg2(3)=rpic(3)*rpic(1)*(rr(3)**2)+(zz**2)*spicuno(3)*spicdue(3)
c
      do i = 1,3           
        hsj(i)=atan2(arg1(i),arg2(i))
      end do	
c
      vxloc=-(hss(1)*hsq(1)+hss(2)*hsq(2)+hss(3)*hsq(3))
      vyloc=+(hsc(1)*hsq(1)+hsc(2)*hsq(2)+hsc(3)*hsq(3))
      vzloc=+(deteta-hsj(1)-hsj(2)-hsj(3))
c
      if(iikk.eq.1) zz = 1.0
      if(zz.gt.0.0) then
        vzloc=vzloc
      else
        vzloc=-vzloc
      end if
c 
      vvx = vxloc*coscsix + vyloc*cosetax + vzloc*coszitax
      vvy = vxloc*coscsiy + vyloc*cosetay + vzloc*coszitay
      vvz = vxloc*coscsiz + vyloc*cosetaz + vzloc*coszitaz
c
      return 
      end
