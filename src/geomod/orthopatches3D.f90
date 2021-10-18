subroutine orthopatches3D(npatches,nbox,filenamei,mnode,nnode,filenameo,xG2,rd2,check,penalty,bofeasible)

character*80 filenamei,filenameo
character*15 filediff,folder
character    lpp,beam,disp,cxmirr
character    bcheck,dcheck,xGcheck,xcfcheck
character    ZC_ycheck,ZC_zcheck
integer      ZC_inode(2),ZC_jnode(2)
integer      iunit,ipatch,icheck
integer      mnode,nnode,knode,inode,jnode
integer      npatches,nbox
integer,dimension(:), allocatable:: n,m,l,e
real         alpha(npatches)
real         A,B,C,BIF_x,BIF_y,BIF_z
real         xori,yori,zori
real         phiA(npatches),phiB(npatches),phiC(npatches)
real         x(mnode*nnode),y(mnode*nnode),z(mnode*nnode)
real         xo(mnode*nnode),yo(mnode*nnode),zo(mnode*nnode)
real         xij(mnode,nnode),yij(mnode,nnode),zij(mnode,nnode)
real         zone_y(mnode*nnode),zone_z(mnode*nnode)
real         eps,pi
real         brange(2),drange(2),xGrange(2),xcfrange(2)
real         bcontrol,dcontrol,xGcontrol,xcfcontrol
real         ZC_ymin,ZC_zmin
real         ZC_ycontrol,ZC_zcontrol
real         check(nbox),penalty,penalty2
real         rlen,rdum,xG,xcf
real         rd0,rl0,rb0,rt0,xG0,xcf0
real         rd1,rl1,rb1,rt1,xG1,xcf1
real         rd2,rl2,rb2,rt2,xG2,xcf2
real         rxmax,rymax,rzmax
real         rxmin,rymin,rzmin
real         rxref,ryref,rdeltaz
real         sink,trim,ymin
logical      bofeasible

namelist/SINK_TRIM/       sink,trim
namelist/ORTHOPATCHES_3D/ BIF_x,BIF_y,BIF_z,&
                          n,m,l,e,&
                          xori,yori,zori
namelist/CONSTRAINTS/     lpp,beam,disp,&
                          bcheck,dcheck,xGcheck,xcfcheck,&
                          brange,drange,xGrange,xcfrange,&
                          ZC_ycheck,ZC_zcheck,&
                          ZC_ymin,ZC_zmin,&
                          ZC_inode,ZC_jnode

! ------------------------------------------------------------

allocate(n(npatches))
allocate(m(npatches))
allocate(l(npatches))
allocate(e(npatches))

iunit=20

open(unit=20,file="SBDF.nml",action="read",recl=10000)
  read(20,NML=SINK_TRIM)
close(20)

open(unit=20,file="SBDF.nml",action="read",recl=10000)
  read(20,NML=ORTHOPATCHES_3D)
close(20)

open(unit=20,file="SBDF.nml",action="read",recl=10000)
  read(20,NML=CONSTRAINTS)
close(20)


open(iunit,file='sinktrim.aux',status='unknown',recl=10000)
  write(iunit,*) sink, 'sinakge'
  write(iunit,*) trim,  'trim'
close(iunit)

write(*,*)
!write(*,*) 'OrthoPatches3D'

x      = 0.
y      = 0.
z      = 0.
xo      = 0.
yo      = 0.
zo      = 0.
xaux   = 0.
yaux   = 0.
zaux   = 0.
xij    = 0.
yij    = 0.
zij    = 0.
check  = 0.
zone_y = 0.
zone_z = 0.
eps    = -1.e-6

iunit  = 20
pi     = acos(-1.)

bofeasible = .true.

!---Read the orthogonal functions coefficients

open(unit=iunit,file='variables.inp',status='old')
do ipatch=1,npatches
  read(iunit,*) alpha(ipatch)
enddo
close(iunit)

!---Read grid input file

open(unit=iunit,file=filenamei,status='old')
  read(iunit,*)
  read(iunit,*) 
  read(iunit,*) (xo(knode),knode=1,mnode*nnode),&
                (yo(knode),knode=1,mnode*nnode),&
                (zo(knode),knode=1,mnode*nnode)
close(iunit)

!---Evaluation of box dimension incresed by BIF

A = abs(maxval(xo)-minval(xo))*(BIF_x/100. + 1.)
B = abs(maxval(yo)-minval(yo))*(BIF_y/100. + 1.)
C = abs(maxval(zo)-minval(zo))*(BIF_z/100. + 1.)

!---Evaluation of phases 

do ipatch=1,npatches
  phiA(ipatch) = -n(ipatch)*pi*xori/A
  phiB(ipatch) = -m(ipatch)*pi*yori/B
  phiC(ipatch) = -l(ipatch)*pi*zori/C
end do

!---Copy grid input file in output file

!filenameo = 'grid_aux.grd'

open(unit=iunit,file=filenameo,status='unknown')
  write(iunit,*) 1
  write(iunit,*) mnode, nnode, 1
  write(iunit,*) (xo(knode),knode=1,mnode*nnode),&
                 (yo(knode),knode=1,mnode*nnode),&
                 (zo(knode),knode=1,mnode*nnode)
close(iunit)

!---Evaluation of minimum value of Y for the overstepping control

do inode=1,mnode
  do jnode=1,nnode
    knode = (jnode-1)*mnode+inode
    yij(inode,jnode) = yo(knode)
  enddo
enddo

ymin = minval(yij(2:mnode,1:nnode-1))

!---Call PANCA to evaluate the hull parameters

call panca

open(999,file='panca.aux',status='old')
  read(999,*) rlen
  read(999,*) rdum, rd0
  read(999,*) rdum, rdum
  read(999,*) xG
  read(999,*) xcf
  read(999,*) rxmin, rxmax
  read(999,*) rymin, rymax
  read(999,*) rzmin, rzmax
close(999)

open(999,file='pancanml.aux',status='old')
 read(999,*) rdeltaz
 read(999,*) cxmirr
 read(999,*) dist_pp
close(999)

rl0  = (rxmax-rxmin)*rlen
!rb0  = (rymax-rymin)*rlen
rb0  = 2*maxval(yo)*rlen/dist_pp
rt0  =        -rzmin*rlen
xG0  =            xG*rlen
xcf0 =           xcf*rlen

write(*,*)   'OBF 3D    Displ       LBP        BOA         T        X CoG      X CoF'
write(*,210) 'Initial',rd0,rl0,rb0,rt0,xG0,xcf0
!write(*,*)   '--------------------------------------------------------------------------'

!---Patches application

x = xo
y = yo
z = zo

do ipatch=1,npatches
  do knode=1,mnode*nnode
    if     (e(ipatch).eq.1) then
      x(knode) = x(knode)+&
                 alpha(ipatch)*&
                 sin((pi*n(ipatch)*xo(knode)/A)+phiA(ipatch))*&
                 sin((pi*m(ipatch)*yo(knode)/B)+phiB(ipatch))*&
                 sin((pi*l(ipatch)*zo(knode)/C)+phiC(ipatch))
    else if(e(ipatch).eq.2) then
      y(knode) = y(knode)+&
                 alpha(ipatch)*&
                 sin((pi*n(ipatch)*xo(knode)/A)+phiA(ipatch))*&
                 sin((pi*m(ipatch)*yo(knode)/B)+phiB(ipatch))*&
                 sin((pi*l(ipatch)*zo(knode)/C)+phiC(ipatch))
    else if(e(ipatch).eq.3) then
      z(knode) = z(knode)+&
                 alpha(ipatch)*&
                 sin((pi*n(ipatch)*xo(knode)/A)+phiA(ipatch))*&
                 sin((pi*m(ipatch)*yo(knode)/B)+phiB(ipatch))*&
                 sin((pi*l(ipatch)*zo(knode)/C)+phiC(ipatch))
    end if
  end do
enddo

!---Fixing possible overstepping of the simmetry y-plane (only buondary TO FIX)

do inode=1,mnode
  do jnode=1,nnode
    knode = (jnode-1)*mnode+inode
    yij(inode,jnode) = y(knode)
  enddo
enddo

penalty2 = 0.

do inode=2,mnode
 do jnode=1,nnode-1
   if(yij(inode,jnode).lt.ymin) then
     bofeasible = .false.
     penalty2 = penalty2 + abs(yij(inode,jnode))
     yij(inode,jnode) = ymin
   end if
 end do
end do

do inode=1,mnode
  do jnode=1,nnode
    knode = (jnode-1)*mnode+inode
    y(knode) = yij(inode,jnode)
  enddo
enddo

!---Writing new geometry on output file

open(unit=iunit,file=filenameo,status='unknown')
  write(iunit,*) 1
  write(iunit,*) mnode, nnode, 1
  write(iunit,*) (x(knode),knode=1,mnode*nnode),&
                 (y(knode),knode=1,mnode*nnode),&
                 (z(knode),knode=1,mnode*nnode)
close(iunit)

!---Scaling

!open(999,file='pancanml.aux',status='old')
! read(999,*) rdeltaz
! read(999,*) cxmirr
!close(999)

call panca()
open(999,file='panca.aux',status='old')
  read(999,*) rlen
  read(999,*) rdum, rd1
  read(999,*) rdum, rdum
  read(999,*) rdum
  read(999,*) rdum
  read(999,*) rxmin, rxmax
  read(999,*) rymin, rymax
  read(999,*) rzmin, rzmax
close(999)

rl1 = (rxmax-rxmin)*rlen
!rb1 = (rymax-rymin)*rlen
rb1 = 2*maxval(y)*rlen/dist_pp
rt1 =        -rzmin*rlen

rdratio = rd0/rd1
rlratio = rl0/rl1
rbratio = rb0/rb1

if (cxmirr=='y') then
  rxref=-rxmax
else
  rxref=rxmin
  ryref=rymin
endif

if (disp.eq.'y') then
  if (lpp.eq.'y') then
    do knode=1,mnode*nnode
      x(knode) = rxref*rlen+(x(knode)/rlen-rxref)*rlratio*rlen
      if(beam.eq.'n') then
        y(knode) =  y(knode)         *sqrt(rdratio/rlratio)
        z(knode) = (z(knode)+rdeltaz)*sqrt(rdratio/rlratio)-rdeltaz
      else
        y(knode) = (y(knode)/rlen)   *rbratio*rlen
        z(knode) = (z(knode)+rdeltaz)*rdratio/(rlratio*rbratio)-rdeltaz
      end if
    end do
  end if
end if


do inode=1,mnode
  do jnode=1,nnode
    knode = (jnode-1)*mnode+inode
    xij(inode,jnode) = x(knode)
    yij(inode,jnode) = y(knode)
    zij(inode,jnode) = z(knode)
    if(((inode.ge.ZC_inode(1)).and.(inode.le.ZC_inode(2))).and.&
       ((jnode.ge.ZC_jnode(1)).and.(jnode.le.ZC_jnode(2)))) then
      zone_y(knode)=yij(inode,jnode)
      zone_z(knode)=zij(inode,jnode)
    end if
  enddo
enddo

!---Writing scaled geometry on output file

open(unit=iunit,file=filenameo,status='unknown')
  write(iunit,*) 1
  write(iunit,*) mnode, nnode, 1
  write(iunit,*) (x(knode),knode=1,mnode*nnode),&
                 (y(knode),knode=1,mnode*nnode),&
                 (z(knode),knode=1,mnode*nnode)
close(iunit)

!---Final check of equality and inequality constraints

call panca()
open(999,file='panca.aux',status='old')
  read(999,*) rlen
  read(999,*) rdum, rd2
  read(999,*) rdum, rdum
  read(999,*) xG
  read(999,*) xcf
  read(999,*) rxmin, rxmax
  read(999,*) rymin, rymax
  read(999,*) rzmin, rzmax
  close(999)

rl2 = (rxmax-rxmin)*rlen
!rb2 = (rymax-rymin)*rlen
rb2 = 2*maxval(y)*rlen/dist_pp
rt2 =        -rzmin*rlen
xG2 =            xG*rlen
xcf2=           xcf*rlen

!write(*,*)   '          displ.  LBP    BOA    T     xG    xcf'
write(*,210) 'Final',rd2,rl2,rb2,rt2,xG2,xcf2
write(*,*)   '--------------------------------------------------------------------------'


icheck = 0

if(bcheck.eq.'y') then
  icheck = icheck + 1
  bcontrol = ((rb2-rb0)/rb0)*100
  if(bcontrol.lt.brange(1)) then 
    check(icheck)   = bcontrol/brange(1)-1.
    write(*,*) ' The new Beam   is over the lower limit =', bcontrol, '%, UNFEASIBLE'
  else if (bcontrol.gt.brange(2)) then
    check(icheck+1) = bcontrol/brange(2)-1.
    write(*,*) ' The new Beam   is over the upper limit =', bcontrol, '%, UNFEASIBLE'
  end if
  icheck = icheck + 1
end if

if(dcheck.eq.'y') then
  icheck = icheck + 1
  dcontrol = ((rt2-rt0)/rt0)*100
  if(dcontrol.lt.drange(1)) then
    check(icheck)   = dcontrol/drange(1)-1.
    write(*,*) ' The new Draft  is over the lower limit =', dcontrol, '%, UNFEASIBLE'
  else if (dcontrol.gt.drange(2)) then
    check(icheck+1) = dcontrol/drange(2)-1.
    write(*,*) ' The new Draft  is over the upper limit =', dcontrol, '%, UNFEASIBLE'
  end if
  icheck = icheck + 1
end if

if(xGcheck.eq.'y') then
  icheck = icheck + 1
  xGcontrol = ((xG2-xG0)/xG0)*100
  if(xGcontrol.lt.xGrange(1)) then
    check(icheck)   = xGcontrol/xGrange(1)-1.
    write(*,*) ' The new xG     is over the lower limit =', xGcontrol, '%, UNFEASIBLE'
  else if (xGcontrol.gt.xGrange(2)) then
    check(icheck+1) = xGcontrol/xGrange(2)-1.
    write(*,*) ' The new xG     is over the upper limit =', xGcontrol, '%, UNFEASIBLE'
  end if
  icheck = icheck + 1
end if

if(xcfcheck.eq.'y') then
  icheck = icheck + 1
  xcfcontrol = ((xcf2-xcf0)/xcf0)*100
  if(xcfcontrol.lt.xcfrange(1)) then
    check(icheck)   = xcfcontrol/xcfrange(1)-1.
    write(*,*) ' The new xcf    is over the lower limit =', xcfcontrol, '%, UNFEASIBLE'
  else if (xcfcontrol.gt.xcfrange(2)) then
    check(icheck+1) = xcfcontrol/xcfrange(2)-1.
    write(*,*) ' The new xcf    is over the upper limit =', xcfcontrol, '%, UNFEASIBLE'
  end if
  icheck = icheck + 1
end if

if(ZC_ycheck.eq.'y') then
   icheck = icheck + 1
   ZC_ycontrol=2*abs(maxval(zone_y))
   if(ZC_ycontrol.lt.ZC_ymin) then
     check(icheck) = ZC_ycontrol/ZC_ymin-1.
     write(*,*) 'The new sonar dome diameter is too small =', ZC_ycontrol, 'UNFEASIBLE'
   end if
end if

if(ZC_zcheck.eq.'y') then
   icheck = icheck + 1
   ZC_zcontrol=abs(maxval(zone_z)-minval(zone_z))
   if(ZC_zcontrol.lt.ZC_zmin) then
     check(icheck) = ZC_zcontrol/ZC_zmin-1.
     write(*,*) 'The new sonar dome height   is too small =', ZC_zcontrol, 'UNFEASIBLE'
   end if
end if

write(*,*)   '--------------------------------------------------------------------------'


penalty = sum(abs(check(:))) + penalty2

!---Genarate geometry file with isodifferences

filediff = 'grid.plt'
!folder   = 'XX.XX.XX'
!
!do ipatch = 1,npatches
!  if (alpha(ipatch).eq.1.) then
!    write(folder(1:2),'(i2.2)') n(ipatch)
!    write(folder(4:5),'(i2.2)') m(ipatch)
!    write(folder(7:8),'(i2.2)') l(ipatch)
!  end if
!end do
!
!call system('mkdir '//folder)

call geodiff(filenamei,filenameo,filediff,mnode,nnode)

210  format(a7,6es11.3)

deallocate(n,m,l,e)

!call system('mv grid.plt '//folder)
!---------------------------------------------------------------------------------------
return
end 

