c-----------------------------------------------------------
      subroutine calve(i,ncar,vx,vy,vz,
     &                 fi1x,fi1y,fi1z,sigma1,
     &                 ncard,ntotd,nplib)
c-----------------------------------------------------------

      implicit real (a-h,o-z)

      real fi1x(ntotd),fi1y(ntotd),fi1z(ntotd)
      real vx(ntotd),vy(ntotd),vz(ntotd)
      real sigma1(ntotd)

c-----------------------------------------------------------------------

      somx=0.0
      somy=0.0
      somz=0.0
      do k=1,ncar+nplib
        somx=somx+sigma1(k)*vx(k)
        somy=somy+sigma1(k)*vy(k)
        somz=somz+sigma1(k)*vz(k)
      end do
      fi1x(i)=somx
      fi1y(i)=somy
      fi1z(i)=somz

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
