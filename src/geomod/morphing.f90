subroutine morphing(nvar,nbox,filenamei,filenameo,xG2,rd2,check,penalty,bofeasible)

implicit none

integer ivar,iunit,mnode,nnode,icheck
integer ionode,jonode,konode,tonode
integer inode,jnode,knode,tnode,ttnode
integer znode,zcnode
integer nbox,nvar

real penalty,penalty2
real beta(nvar)
real,dimension(:),allocatable:: x,y,z,xo,yo,zo,xnew,ynew,znew
real,dimension(:),allocatable:: ak,zone_y,zone_z
real,dimension(:,:),allocatable:: nVR,eigx,eigy,eigz
real,dimension(:,:),allocatable:: xij,yij,zij

character*80 filenamei,filenameo
character*15 filediff
character    lpp,beam,disp,cxmirr
character    bcheck,dcheck,xGcheck,xcfcheck
character    ZC_ycheck,ZC_zcheck
integer      ZC_inode(2),ZC_jnode(2)
real         brange(2),drange(2),xGrange(2),xcfrange(2)
real         bcontrol,dcontrol,xGcontrol,xcfcontrol
real         ZC_ymin,ZC_zmin
real         ZC_ycontrol,ZC_zcontrol
real         check(nbox)


real  sink,trim,dist_pp,ymin
real  rlen,rdum,xG,xcf
real  rd0,rl0,rb0,rt0,xG0,xcf0
real  rd1,rl1,rb1,rt1,xG1,xcf1
real  rd2,rl2,rb2,rt2,xG2,xcf2
real  rxmax,rymax,rzmax
real  rxmin,rymin,rzmin
real  rxref,ryref,rdeltaz
real  rdratio,rlratio,rbratio

logical bofeasible

namelist/SINK_TRIM/       sink,trim
namelist/CONSTRAINTS/     lpp,beam,disp,&
                          bcheck,dcheck,xGcheck,xcfcheck,&
                          brange,drange,xGrange,xcfrange,&
                          ZC_ycheck,ZC_zcheck,&
                          ZC_ymin,ZC_zmin,&
                          ZC_inode,ZC_jnode

!-----------------------------------------------------------------------------------------

iunit = 20
check = 0.

open(unit=20,file="SBDF.nml",action="read",recl=10000)
  read(20,NML=SINK_TRIM)
close(20)

open(unit=20,file="SBDF.nml",action="read",recl=10000)
  read(20,NML=CONSTRAINTS)
close(20)

!open(unit=20,file="SBDF.nml",action="read",recl=10000)
!  read(20,NML=KLE_NODES)
!close(20)

open(iunit,file='sinktrim.aux',status='unknown',recl=10000)
  write(iunit,*) sink, 'sinakge'
  write(iunit,*) trim,  'trim'
close(iunit)

!---Read the coefficients

open(unit=iunit,file='variables.inp',status='old')
  do ivar=1,nvar
    read(iunit,*) beta(ivar)
  enddo
close(iunit)

write(*,*) 'Morphing from KLE eigenvector'

bofeasible = .true.

! --- Original grid

open(iunit,file=filenamei,status='old',recl=10000)
  read(iunit,*)
  read(iunit,*) inode,jnode
  tnode = inode*jnode

  allocate(x(tnode))
  allocate(y(tnode))
  allocate(z(tnode))

  read(iunit,*) (x(knode),knode=1,tnode),&
                (y(knode),knode=1,tnode),&
                (z(knode),knode=1,tnode)
close(iunit)

!---Copy grid input file in output file

open(unit=iunit,file=filenameo,status='unknown')
  write(iunit,*) 1
  write(iunit,*) inode, jnode, 1
  write(iunit,*) (x(knode),knode=1,tnode),&
                 (y(knode),knode=1,tnode),&
                 (z(knode),knode=1,tnode)
close(iunit)

!---Evaluation of minimum value of Y for the overstepping control

allocate(xij(inode,jnode))
allocate(yij(inode,jnode))
allocate(zij(inode,jnode))
znode = (ZC_inode(2)-ZC_inode(1)+1)*(ZC_jnode(2)-ZC_jnode(1)+1)
!write(*,*) znode
allocate(zone_y(tnode),zone_z(tnode))
zone_y = 0.
zone_z = 0.

do mnode=1,inode
  do nnode=1,jnode
    knode = (nnode-1)*inode+mnode
    yij(mnode,nnode) = y(knode)
  enddo
enddo

ymin = minval(yij(2:inode,1:jnode-1))

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
rb0  = 2*maxval(y)*rlen/dist_pp
rt0  =        -rzmin*rlen
xG0  =            xG*rlen
xcf0 =           xcf*rlen

!write(*,*) '  Initial displacement, LBP, BOA, T, xG, xcf:'
!write(*,*) rd0,rl0,rb0,rt0,xG0,xcf0
write(*,*)   'KLE EIG   Displ       LBP        BOA         T        X CoG      X CoF'
write(*,210) 'Initial',rd0,rl0,rb0,rt0,xG0,xcf0
!write(*,*)   '--------------------------------------------------------------------------'

allocate(xo(tnode))
allocate(yo(tnode))
allocate(zo(tnode))

ttnode = tnode

allocate(nVR(3*ttnode,nvar))
allocate(ak(nvar))
allocate(eigx(tnode,nvar))
allocate(eigy(tnode,nvar))
allocate(eigz(tnode,nvar))

! --- KLE eigenvector

open(iunit,file='eigenvecto.out',status='unknown',recl=10000)
  do knode=1,3*ttnode
    read(iunit,*) nVR(knode,:)
  end do
close(iunit)

open(iunit,file='alfak.out',status='unknown',recl=10000)
  read(iunit,*) 
  do ivar=1,nvar
    read(iunit,*) rdum,rdum,ak(ivar)
  end do
close(iunit)

do knode=1,tnode
  eigx(knode,:) = sqrt(3.*ak(:))*nVR(         knode,:)
  eigy(knode,:) = sqrt(3.*ak(:))*nVR(  ttnode+knode,:)
  eigz(knode,:) = sqrt(3.*ak(:))*nVR(2*ttnode+knode,:)
end do

do knode=1,tnode
  xo(knode) = x(knode) + sum(beta(:)*eigx(knode,:))
  yo(knode) = y(knode) + sum(beta(:)*eigy(knode,:))
  zo(knode) = z(knode) + sum(beta(:)*eigz(knode,:))
end do

!---Fixing possible overstepping of the simmetry y-plane (only buondary TO FIX)
zcnode = 0 
do mnode=1,inode
  do nnode=1,jnode
    knode = (nnode-1)*inode+mnode
    xij(mnode,nnode) = xo(knode)
    yij(mnode,nnode) = yo(knode)
    zij(mnode,nnode) = zo(knode)
    if(((mnode.ge.ZC_inode(1)).and.(mnode.le.ZC_inode(2))).and.&
       ((nnode.ge.ZC_jnode(1)).and.(nnode.le.ZC_jnode(2)))) then
      zcnode = zcnode+1
      zone_y(zcnode)=yij(mnode,nnode)
      zone_z(zcnode)=zij(mnode,nnode)
    end if
  enddo
enddo
!write(*,*) zcnode

penalty2 = 0.

do mnode=2,inode
 do nnode=1,jnode-1
   if(yij(mnode,nnode).lt.ymin) then
     bofeasible = .false.
     penalty2 = penalty2 + abs(yij(mnode,nnode))
     yij(mnode,nnode) = ymin
   end if
 end do
end do

do mnode=1,inode
  do nnode=1,jnode
    knode = (nnode-1)*inode+mnode
    yo(knode) = yij(mnode,nnode)
  enddo
enddo

!---Writing new geometry on output file

open(unit=iunit,file=filenameo,status='unknown')
  write(iunit,*) 1
  write(iunit,*) inode, jnode, 1
  write(iunit,*) (xo(knode),knode=1,tnode),&
                 (yo(knode),knode=1,tnode),&
                 (zo(knode),knode=1,tnode)
close(iunit)

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
rb1 = 2*maxval(yo)*rlen/dist_pp
rt1 =        -rzmin*rlen
xG1 =            xG*rlen
xcf1=           xcf*rlen

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
    do knode=1,inode*jnode
      xo(knode) = rxref*rlen+(xo(knode)/rlen-rxref)*rlratio*rlen
      if(beam.eq.'n') then
        yo(knode) =  yo(knode)         *sqrt(rdratio/rlratio)
        zo(knode) = (zo(knode)+rdeltaz)*sqrt(rdratio/rlratio)-rdeltaz
      else
        yo(knode) = (yo(knode)/rlen)   *rbratio*rlen
        zo(knode) = (zo(knode)+rdeltaz)*rdratio/(rlratio*rbratio)-rdeltaz
      end if
    end do
  end if
end if

do mnode=1,inode
  do nnode=1,jnode
    knode = (nnode-1)*inode+mnode
    xij(mnode,nnode) = xo(knode)
    yij(mnode,nnode) = yo(knode)
    zij(mnode,nnode) = zo(knode)
    if(((mnode.ge.ZC_inode(1)).and.(mnode.le.ZC_inode(2))).and.&
       ((nnode.ge.ZC_jnode(1)).and.(nnode.le.ZC_jnode(2)))) then
      zone_y(knode)=yij(mnode,nnode)
      zone_z(knode)=zij(mnode,nnode)
    end if
  enddo
enddo

!---Writing scaled geometry on output file

open(unit=iunit,file=filenameo,status='unknown')
  write(iunit,*) 1
  write(iunit,*) inode, jnode, 1
  write(iunit,*) (xo(knode),knode=1,inode*jnode),&
                 (yo(knode),knode=1,inode*jnode),&
                 (zo(knode),knode=1,inode*jnode)
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
rb2 = 2*maxval(yo)*rlen/dist_pp
rt2 =        -rzmin*rlen
xG2 =            xG*rlen
xcf2=           xcf*rlen

write(*,210) 'Final',rd2,rl2,rb2,rt2,xG2,xcf2
write(*,*)   '--------------------------------------------------------------------------'


! --- Box constraints check zone

icheck = 0

if(bcheck.eq.'y') then
  icheck = icheck + 1
  bcontrol = ((rb2-rb0)/rb0)*100
  if(bcontrol.lt.brange(1)) then
  !  check(icheck)   = bcontrol/brange(1)-1.
    write(*,*) ' The new Beam is over the lower limit =', bcontrol, '%, UNFEASIBLE'
  else if (bcontrol.gt.brange(2)) then
  !  check(icheck+1) = bcontrol/brange(2)-1.
    write(*,*) ' The new Beam is over the upper limit =', bcontrol, '%, UNFEASIBLE'
  end if
  !icheck = icheck + 1
  check(icheck) = abs(bcontrol)/abs(brange(2))-1.
end if

if(dcheck.eq.'y') then
  icheck = icheck + 1
  dcontrol = ((rt2-rt0)/rt0)*100
  if(dcontrol.lt.drange(1)) then
  !  check(icheck)   = dcontrol/drange(1)-1.
    write(*,*) ' The new Draft is over the lower limit =', dcontrol, '%, UNFEASIBLE'
  else if (dcontrol.gt.drange(2)) then
  !  check(icheck+1) = dcontrol/drange(2)-1.
    write(*,*) ' The new Draft is over the upper limit =', dcontrol, '%, UNFEASIBLE'
  end if
  !icheck = icheck + 1
  check(icheck) = abs(dcontrol)/abs(drange(2))-1.
end if

if(xGcheck.eq.'y') then
  icheck = icheck + 1
  xGcontrol = ((xG2-xG0)/xG0)*100
  if(xGcontrol.lt.xGrange(1)) then
    check(icheck)   = xGcontrol/xGrange(1)-1.
    write(*,*) ' The new xG is over the lower limit =', xGcontrol, '%, UNFEASIBLE'
  else if (xGcontrol.gt.xGrange(2)) then
    check(icheck+1) = xGcontrol/xGrange(2)-1.
    write(*,*) ' The new xG is over the upper limit =', xGcontrol, '%, UNFEASIBLE'
  end if
  icheck = icheck + 1
end if

if(xcfcheck.eq.'y') then
  icheck = icheck + 1
  xcfcontrol = ((xcf2-xcf0)/xcf0)*100
  if(xcfcontrol.lt.xcfrange(1)) then
    check(icheck)   = xcfcontrol/xcfrange(1)-1.
    write(*,*) ' The new xcf is over the lower limit =', xcfcontrol, '%, UNFEASIBLE'
  else if (xcfcontrol.gt.xcfrange(2)) then
    check(icheck+1) = xcfcontrol/xcfrange(2)-1.
    write(*,*) ' The new xcf is over the upper limit =', xcfcontrol, '%, UNFEASIBLE'
  end if
  icheck = icheck + 1
end if

if(ZC_ycheck.eq.'y') then
   icheck = icheck + 1
   ZC_ycontrol=2*abs(maxval(zone_y))
   if(ZC_ycontrol.lt.ZC_ymin) then
   !  check(icheck) = ZC_ycontrol/ZC_ymin-1.
     write(*,*) 'The new sonar dome diameter is too small =', ZC_ycontrol, 'UNFEASIBLE'
   end if
!   check(icheck) = ZC_ycontrol/ZC_ymin-1.
   check(icheck) = 1.-ZC_ycontrol/ZC_ymin
end if

if(ZC_zcheck.eq.'y') then
   icheck = icheck + 1
   ZC_zcontrol=abs(maxval(zone_z)-minval(zone_z))
   if(ZC_zcontrol.lt.ZC_zmin) then
   !  check(icheck) = ZC_zcontrol/ZC_zmin-1.
     write(*,*) 'The new sonar dome height is too small =', ZC_zcontrol, 'UNFEASIBLE'
   end if
!   check(icheck) = ZC_zcontrol/ZC_zmin-1.
   check(icheck) = 1.-ZC_zcontrol/ZC_zmin
end if

penalty = sum(abs(check(:))) + penalty2

!---Genarate geometry file with isodifferences

filediff = 'grid.plt'
call geodiff(filenamei,filenameo,filediff,inode,jnode)

210  format(a7,6es11.3)

!---------------------------------------------------------------------------------------

deallocate(x,y,z,xo,yo,zo)!,xnew,ynew,znew)
deallocate(ak,zone_y,zone_z)
deallocate(nVR,eigx,eigy,eigz)
deallocate(xij,yij,zij)


return
end 
