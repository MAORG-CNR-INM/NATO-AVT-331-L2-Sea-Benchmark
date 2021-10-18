c------------------------------------------------------------
      subroutine hemiqua(xn,yn,zn,xnt,ynt,znt,xv1,xv2,xv3,xv4,
     &                   xh,yh,zh,yv1,yv2,yv3,yv4,zv1,zv2,zv3,zv4,
     &                   coszitax,coszitay,coszitaz,
     &                   vvx,vvy,vvz,pig,iikk)
c------------------------------------------------------------
c
      integer i,ii,iikk
      real hsc(4),hss(4),hsdist(4),spicuno(4),spicdue(4),
     &     hsq(4),hsj(4),rpic(4),rr(4),arg1(4),arg2(4)
      real x1,x2,x3,x4,y1,y2,y3,y4
      real radx,coscsix,coscsiy,coscsiz,rady,cosetax,cosetay,cosetaz
      real xx,yy,zz,deteta,vxloc,vyloc,vzloc
      real xn,yn,zn,xnt,ynt,znt,xv1,xv2,xv3,xv4,
     &     xh,yh,zh,yv1,yv2,yv3,yv4,zv1,zv2,zv3,zv4,
     &     coszitax,coszitay,coszitaz,
     &     vvx,vvy,vvz,pig
c------------------------------------------------------------

      if(iikk.eq.1) then
        vvx = 2.*pig*coszitax
        vvy = 2.*pig*coszitay
        vvz = 2.*pig*coszitaz
        return
      end if

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

      spicuno(1)=x1xx*hsc(1)+y1yy*hss(1)	   
      spicdue(1)=x2xx*hsc(1)+y2yy*hss(1)	   
      spicuno(2)=x2xx*hsc(2)+y2yy*hss(2)	
      spicdue(2)=x3xx*hsc(2)+y3yy*hss(2)	
      spicuno(3)=x3xx*hsc(3)+y3yy*hss(3)	
      spicdue(3)=x4xx*hsc(3)+y4yy*hss(3)	
      spicuno(4)=x4xx*hsc(4)+y4yy*hss(4)	
      spicdue(4)=x1xx*hsc(4)+y1yy*hss(4)	
c
      rr(1)=-x1xx*hss(1)+y1yy*hsc(1)       
      rr(2)=-x2xx*hss(2)+y2yy*hsc(2)
      rr(3)=-x3xx*hss(3)+y3yy*hsc(3)
      rr(4)=-x4xx*hss(4)+y4yy*hsc(4)

c--------------------------------------------------
      deteta=2.0*pig
      do ii=1,4      
        if(rr(ii).lt.0.0) then
          deteta=0.0
          go to 9191
	end if
      end do
 9191 continue	
c---------------------------------------------------
      rpic(1)=sqrt(x1xx**2+y1yy**2+zz**2)
      rpic(2)=sqrt(x2xx**2+y2yy**2+zz**2)
      rpic(3)=sqrt(x3xx**2+y3yy**2+zz**2)
      rpic(4)=sqrt(x4xx**2+y4yy**2+zz**2)
c
      hsq(1)=log((rpic(1)+rpic(2)+hsdist(1))/
     &           (rpic(1)+rpic(2)-hsdist(1)))
      hsq(2)=log((rpic(2)+rpic(3)+hsdist(2))/
     &           (rpic(2)+rpic(3)-hsdist(2)))
      hsq(3)=log((rpic(3)+rpic(4)+hsdist(3))/
     &           (rpic(3)+rpic(4)-hsdist(3)))
      hsq(4)=log((rpic(4)+rpic(1)+hsdist(4))/
     &           (rpic(4)+rpic(1)-hsdist(4)))    
      arg1(1)=rr(1)*abs(zz)*(rpic(1)*spicdue(1)-rpic(2)*spicuno(1))
      arg2(1)=rpic(1)*rpic(2)*(rr(1)**2)+(zz**2)*spicuno(1)*spicdue(1)
      arg1(2)=rr(2)*abs(zz)*(rpic(2)*spicdue(2)-rpic(3)*spicuno(2))
      arg2(2)=rpic(2)*rpic(3)*(rr(2)**2)+(zz**2)*spicuno(2)*spicdue(2)
      arg1(3)=rr(3)*abs(zz)*(rpic(3)*spicdue(3)-rpic(4)*spicuno(3))
      arg2(3)=rpic(3)*rpic(4)*(rr(3)**2)+(zz**2)*spicuno(3)*spicdue(3)
      arg1(4)=rr(4)*abs(zz)*(rpic(4)*spicdue(4)-rpic(1)*spicuno(4)) 
      arg2(4)=rpic(4)*rpic(1)*(rr(4)**2)+(zz**2)*spicuno(4)*spicdue(4)
c
      do i = 1,4           
        if(arg1(i).eq.0.and.arg2(i).eq.0) then
          hsj(i) = 0.0
        else
          hsj(i)=atan2(arg1(i),arg2(i))
        end if
      end do	
c
      vxloc=-(hss(1)*hsq(1)+hss(2)*hsq(2)+hss(3)*hsq(3)+hss(4)*hsq(4))
      vyloc=+(hsc(1)*hsq(1)+hsc(2)*hsq(2)+hsc(3)*hsq(3)+hsc(4)*hsq(4))
      vzloc=+(deteta-hsj(1)-hsj(2)-hsj(3)-hsj(4))
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
c------------------------------------------------------------
