c ----------------------------------------------------------------------

      subroutine tecout(fr,fi0x,fi0y,fi0z,fi1x,fi1y,fi1z,
     &                  eta,presto,xns,yns,zns,
     &                  ncar,nulib,ncard,nplib,
     &                  ntra,ntraint,ntrapoc,nlo,
     &                  nlovallec,ntotd,nplibd)

c ----------------------------------------------------------------------

      integer ncar,nulib,ncard,nplib,ntra,ntraint,nlo,nlovallec
      integer i,j,ifr,ipan,ntrapoc
      real xns(ntotd),yns(ntotd),zns(ntotd)
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd)
      real fi1x(ntotd),fi1y(ntotd),fi1z(ntotd)
      real presto(ncard)
      real eta(nplibd)
      real x,y,z,vx,vy,vz,et,fr
      character*14 fina

      real pred (100*nlo*ntra,3)
      real train(nulib+ntrapoc*nlovallec,3)

c ----------------------------------------------------------------------

      fina(1:6) = 'solfr.'
      ifr = fr*1000
      write(fina(7:10),'(i4.4)') ifr
      fina(11:14) = '.plt'
      open (99,file=fina,form='formatted',status='unknown',recl=64000)
      write(99,*) 'VARIABLES = "x" "y" "z" "vx" "vy" "vz" "eta"'

c -- Zona tra i due scafi

      if(ntraint.gt.0) then
        write(99,*)'ZONE T="Int.Prua"I=',ntraint,'J=',nlo-nlovallec+1
        do i=1,nlo-nlovallec+1
          do j=1,ntraint
            ipan = ntra*(i-1) + j
            x  = xns(ncar+ipan)
            y  = yns(ncar+ipan)
            z  = zns(ncar+ipan)
            vx = 1.+fi0x(ncar+ipan)+fi1x(ncar+ipan)
            vy =    fi0y(ncar+ipan)+fi1y(ncar+ipan)
            vz =    fi0z(ncar+ipan)+fi1z(ncar+ipan)
            et = eta(ipan)
            write(99,*) x,y,z,vx,vy,vz,et
          end do
        end do

        write(99,*)'ZONE T="Int.Poppa"I=',ntraint,'J=',nlovallec
        do i=nlo-nlovallec+1,nlo
          do j=1,ntraint
            ipan = ntra*(i-1) + j
            x  = xns(ncar+ipan)
            y  = yns(ncar+ipan)
            z  = zns(ncar+ipan)
            vx = 1.+fi1x(ncar+ipan)
            vy =    fi1y(ncar+ipan)
            vz =    fi1z(ncar+ipan)
            et = eta(ipan)
            write(99,*) x,y,z,vx,vy,vz,et
          end do
        end do
      end if

c -- Zona esterna

      write(99,*)'ZONE T="Est.Prua"I=',ntra-ntraint,'J=',nlo-nlovallec+1
      do i=1,nlo-nlovallec+1
        do j=1,ntra-ntraint
          ipan = ntra*(i-1) + j + ntraint
          x  = xns(ncar+ipan)
          y  = yns(ncar+ipan)
          z  = zns(ncar+ipan)
          vx = 1.+fi1x(ncar+ipan)
          vy =    fi1y(ncar+ipan)
          vz =    fi1z(ncar+ipan)
          et = eta(ipan)
          write(99,*) x,y,z,vx,vy,vz,et
        end do
      end do
      
      write(99,*)'ZONE T="Est.Poppa"I=',ntra-ntraint,'J=',nlovallec
      do i=nlo-nlovallec+1,nlo
        do j=1,ntra-ntraint
          ipan = ntra*(i-1) + j + ntraint
          x  = xns(ncar+ipan)
          y  = yns(ncar+ipan)
          z  = zns(ncar+ipan)
          vx = 1.+fi1x(ncar+ipan)
          vy =    fi1y(ncar+ipan)
          vz =    fi1z(ncar+ipan)
          et = eta(ipan)
          write(99,*) x,y,z,vx,vy,vz,et
        end do
      end do
      
c -- Transom

      if(ntrapoc.gt.0) then
        if(ntrapoc.gt.1) then
      write(99,*) 'ZONE T="Transom" I=',ntrapoc,' J=',nlovallec
        else
          write(99,*) 'ZONE T="Transom" I=2 J=',nlovallec
        end if
        do i=1,nlovallec
          do j=1,ntrapoc
            ipan = ntrapoc*(i-1) + j + nulib
            x  = xns(ncar+ipan)
            y  = yns(ncar+ipan)
            z  = zns(ncar+ipan)
            vx = 1.+fi1x(ncar+ipan)
            vy =    fi1y(ncar+ipan)
            vz =    fi1z(ncar+ipan)
            et = eta(ipan)
            write(99,*) x,y,z,vx,vy,vz,et
            if(ntrapoc.eq.1) write(99,*) x,y,z,vx,vy,vz,et
          end do
        end do
      end if

      close(99)


c --- Interpolation by Gaussian weighting

      ni=10*nlo
      nj=10*ntra
      npred=ni*nj
      ntrain=nulib+ntrapoc*nlovallec

      xmin=minval(xns(ncar+1:ncar+nulib))
      xmax=maxval(xns(ncar+1:ncar+nulib))
      ymin=0. !minval(yns(ncar+1:ncar+nulib))
      ymax=maxval(yns(ncar+1:ncar+nulib))

      ! free surf.
      do i=1,nulib
         train(i,1)=xns(ncar+i)
         train(i,2)=yns(ncar+i)
         train(i,3)=eta(i)
      enddo
      ! transom
      do i=1,ntrapoc*nlovallec
         train(nulib+i,1)=xns(ncar+nulib+i)
         train(nulib+i,2)=yns(ncar+nulib+i)
         train(nulib+i,3)=eta(nulib+i)
      enddo
 
      k=1
      do i=1,ni
        do j=1,nj
           pred(k,1)=((1.*i-1.)/(ni-1.)*(xmax-xmin)+xmin)
           pred(k,2)=((1.*j-1.)/(nj-1.)*(ymax-ymin)+ymin)
           k=k+1
        enddo
      enddo

      call gauss(2,ntrain,npred,train,pred,0.75*(xmax-xmin)/nlo)


      fina(1:6) = 'intfr.'
      ifr = fr*1000
      write(fina(7:10),'(i4.4)') ifr
      fina(11:14) = '.plt'
      open (99,file=fina,form='formatted',status='unknown',recl=64000)
      write(99,*) 'VARIABLES = "X" "Y" "Z" "eta"'

      write(99,*) 'ZONE T = "Free-surface domain" I=',nj,'J=',ni
      do i=1,npred  
        write(99,*) pred(i,1), pred(i,2), pred(i,3), pred(i,3) 
      enddo

      close(99)

      open(99,file='intfr.aux',status='unknown',recl=10000)
        write(99,*) nj,ni,fr
      close(99)

c      rewind(666)
c      do i=1,ntrain
c        write(666,*) train(i,:)
c      enddo

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
