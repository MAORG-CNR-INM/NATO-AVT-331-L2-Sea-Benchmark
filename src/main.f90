Program main

implicit none

real t_start,t_end,t_exec

character*6 dire
character*12 command1
character*18 command2

integer iunit,ierrmpi,clock
integer iaso,isolv,igeomod,igrid,numfilest,numfilint
integer ndv,nfunc,nfcw,nfsk,nbox
integer nsample,ngrid
integer ialgo

namelist/MAIN_PARAMETERS/         iaso,isolv,igeomod,ngrid,igrid,numfilest,numfilint
namelist/PROBLEM_PARAMETERS/      ndv,nfunc,nfcw,nfsk,nbox
namelist/OPTIMIZATION_PARAMETERS/ ialgo

! **********************************************************************

call system_clock(count=clock)
t_start = clock

  write(*,*)
  write(*,*) '------------------------------------------------------------------------- '
  write(*,*) '|         _____   _      _   _______       _____   _____    ___          |'
  write(*,*) '|        /  _  \ | \    / | |__   __|     |___  \ |___  \  /   \         |'
  write(*,*) '|        | /_\ | \ \    / /    | |    ___   __| |   __| | /_/| |         |'
  write(*,*) '|        |  _  |  \ \  / /     | |   |___| |__  |  |__  |    | |         |'
  write(*,*) '|        | | | |   \ \/ /      | |         ___| |  ___| |   _| |_        |'
  write(*,*) '|        |_| |_|    \__/       |_|        |_____/ |_____/  |_____|       |'
  write(*,*) '|                                                                        |'
  write(*,*) '|                 NATO-AVT-331, L2 Sea benchmark problem                 |'
  write(*,*) '|                                  v1.0                                  |'
  write(*,*) '|                       22 Sep. 2020 ... release 1.0                     |'
  write(*,*) '|                           CNR-INM, Rome, Italy                         |'
  write(*,*) '|                          Serani A. and  Diez M.                        |'
  write(*,*) '|                                                                        |'
  write(*,*) '--------------------------------------------------------------------------'
  write(*,*)

! ----------------------------------------------------------------------
! -- Initialize Dirs
! ----------------------------------------------------------------------

command1 = 'mkdir CPU000'
call system(command1)

command2 = 'cp -fr *.* CPU000/'
call system(command2)

dire = 'CPU000'
call chdir(dire)

! ----------------------------------------------------------------------

iunit = 20

open(iunit,file="SBDF.nml",action="read",recl=10000)
  read(iunit,NML=MAIN_PARAMETERS)
close(iunit)

open(iunit,file="SBDF.nml",action="read",recl=10000)
  read(iunit,NML=PROBLEM_PARAMETERS)
close(iunit)

if(iaso.eq.1) then
        call analysis(nfunc,nbox)
end if

call system_clock(count=clock)
t_end = clock
t_exec = (t_end-t_start)/10000

open(33,file='elapsed_CPUtime.out',status='unknown',recl=10000)
  write(33,*) "Elapsed CPU time = ",t_exec,"seconds"
close(33)

call system('sh ../../../scripts/clean.sh')

end program


