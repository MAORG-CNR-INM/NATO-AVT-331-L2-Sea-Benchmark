subroutine analysis(nfunc,nbox)

implicit none

character*80, dimension(:), allocatable:: filenamei,filenameo

integer     numfilest,numfilint
integer     iunit
integer     iaso,isolv,igeomod,igrid
integer     ndv,nfunc,nfcw,nfsk,nbox,nfr
integer     inode,jnode
integer     nobj,ngrid

real        weight(nfunc),constr(nbox)
real        fobj(nfunc),fdum(1)
real        xG,displ,penalty
logical     bofeasible
real, dimension(:), allocatable:: faux,U


namelist/MAIN_PARAMETERS/    iaso,isolv,igeomod,ngrid,igrid,numfilest,numfilint
namelist/PROBLEM_PARAMETERS/ ndv,nfunc,nfcw,nfsk,nbox
namelist/GRID_PARAMETERS/    filenamei,filenameo,inode,jnode
namelist/NUMBER_FROUDE/      nfr

!-------------------------------------------------------------------------------

weight = 0.
constr = 0.
fobj   = 0.
fdum   = 0.
bofeasible = .true.

iunit  = 20

open(iunit,file="SBDF.nml",action="read",recl=10000)
  read(iunit,NML=MAIN_PARAMETERS)
close(iunit)

open(iunit,file="SBDF.nml",action="read",recl=10000)
  read(iunit,NML=PROBLEM_PARAMETERS)
close(iunit)

open(iunit,file="SBDF.nml",action="read",recl=10000)
  read(iunit,NML=NUMBER_FROUDE)
close(iunit)

allocate(faux(nfr))
allocate(U(ngrid))

allocate(filenamei(numfilest+numfilint))
allocate(filenameo(numfilest+numfilint))

open(iunit,file="SBDF.nml",action="read",recl=10000)
  read(iunit,NML=GRID_PARAMETERS)
close(iunit)

if(igeomod.eq.1) then
  call orthopatches3D(ndv,nbox,filenamei(1),inode,jnode,filenameo(1),xG,displ,constr,penalty,bofeasible)
else if(igeomod.eq.2) then
  call morphing(ndv,nbox,filenamei(1),filenameo(1),xG,displ,constr,penalty,bofeasible)
end if 

call regrid(igrid,ngrid,filenameo(1),filenameo(1))   

U(7) = 12.2435756853397
U(6) = 8.94650774731822
U(5) = 5.78367103694875
U(4) = 3.79630512514897
U(3) = 2.68553039332537
U(2) = 1.74724672228843
U(1) = 0.55530393325387

if(isolv.eq.1) then
  call freewarp(igrid,ngrid,faux,nfr)
  write(*,*) '============================================'
  write(*,*) 'Objective function:                         '
  write(*,*) '============================================'
  write(*,*) 'Total resistance [N]            =',faux(1:nfr)
  write(*,*)
  write(*,*) '============================================'
  write(*,*) 'Constraints (non-dimensional, <0 satisifed) '
  write(*,*) '============================================'
  write(*,*) 'Beam                            =',constr(1)
  write(*,*) 'Draft                           =',constr(2)
  write(*,*) 'Sonar dome diameter             =',constr(3)
  write(*,*) 'Sonar dome height               =',constr(4)
  write(*,*)
 
  open(unit=33,file='objective.out',status='unknown',recl=10000)
    write(33,*) "#Grid level,U grid,objective,constraints(1,2,3,4)"
    write(33,*) igrid,U(igrid),faux(1:nfr),constr
  close(33)
end if

deallocate(faux,U,filenamei,filenameo)

return
end

