subroutine geodiff(filenamei,filenameo,filediff,mnode,nnode)

implicit none

character*80 filenamei,filenameo
character*15 filediff,filemode
integer      iunit
integer      knode,mnode,nnode,inode,jnode
real         x(mnode*nnode),y(mnode*nnode),z(mnode*nnode)
real         xnew(mnode*nnode),ynew(mnode*nnode),znew(mnode*nnode)
real         xdiff(mnode*nnode),ydiff(mnode*nnode),zdiff(mnode*nnode)
real         xij(mnode,nnode),yij(mnode,nnode),zij(mnode,nnode)
real         xdiffij(mnode,nnode),ydiffij(mnode,nnode),zdiffij(mnode,nnode)

!--------------------------------------------------------------------------------

!write(*,*) '------------------------------------------------------------------'
!write(*,*) 'GeometryDifference'

x       = 0.
y       = 0.
z       = 0.
xnew    = 0.
ynew    = 0.
znew    = 0.
xdiff   = 0.
ydiff   = 0.
zdiff   = 0.
xij     = 0.
yij     = 0.
zij     = 0.
xdiffij = 0.
ydiffij = 0.
zdiffij = 0.

iunit=20

!---------------------------------------------------------------------------------

open(unit=iunit,file=filenamei,status='old',recl=10000)
  read(iunit,*)
  read(iunit,*) 
  read(iunit,*) (x(knode),knode=1,mnode*nnode),&
                (y(knode),knode=1,mnode*nnode),&
                (z(knode),knode=1,mnode*nnode)
close(iunit)

open(unit=iunit,file=filenameo,status='old',recl=10000)
  read(iunit,*)
  read(iunit,*) 
  read(iunit,*) (xnew(knode),knode=1,mnode*nnode),&
                (ynew(knode),knode=1,mnode*nnode),&
                (znew(knode),knode=1,mnode*nnode)
close(iunit)

do knode=1,mnode*nnode
  xdiff(knode)=xnew(knode)-x(knode)
  ydiff(knode)=ynew(knode)-y(knode)
  zdiff(knode)=znew(knode)-z(knode)
end do

open(unit=iunit,file='griddiff.kle',status='unknown',recl=10000)
  do knode=1,mnode*nnode
    write(iunit,*) xdiff(knode),ydiff(knode),zdiff(knode)
  end do
close(iunit)

do inode=1,mnode
  do jnode=1,nnode
    knode=(jnode-1)*mnode+inode
    xij(inode,jnode)=xnew(knode)
    yij(inode,jnode)=ynew(knode)
    zij(inode,jnode)=znew(knode)
    xdiffij(inode,jnode)=xdiff(knode)
    ydiffij(inode,jnode)=ydiff(knode)
    zdiffij(inode,jnode)=zdiff(knode)
  enddo
enddo

open(iunit,file=filediff,form='formatted',status='unknown',recl=10000)
  write(iunit,*) 'VARIABLES = "X" "Y" "Z" "XDIFF" "YDIFF" "ZDIFF" "MODDIFF"'
  write(iunit,*) 'ZONE I=',mnode,' J=',nnode
  do jnode=1,nnode
    do inode=1,mnode
      write(iunit,100) xij(inode,jnode),yij(inode,jnode),zij(inode,jnode),&
                       xdiffij(inode,jnode),ydiffij(inode,jnode),zdiffij(inode,jnode),&
                       sqrt(xdiffij(inode,jnode)**2+ydiffij(inode,jnode)**2+zdiffij(inode,jnode)**2)
    end do
  end do
close(iunit)

filemode = filediff
write(filemode(13:15),'(A3)') 'mod'

open(iunit,file=filemode,status='unknown',recl=10000)
  write(iunit,*) 'I, J, PHI_X, PHI_Y, PHI_Z, |PHI|'
  do jnode=1,nnode
    do inode=1,mnode
      write(iunit,*) inode,jnode,xdiffij(inode,jnode),ydiffij(inode,jnode),zdiffij(inode,jnode),&
                       sqrt(xdiffij(inode,jnode)**2+ydiffij(inode,jnode)**2+zdiffij(inode,jnode)**2)
    end do
  end do
close(iunit)



100 format(20e16.8)

!---------------------------------------------------------------------------------

return
end

