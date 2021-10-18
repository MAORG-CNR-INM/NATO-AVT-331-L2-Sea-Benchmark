c----------------------------------------------------------------------*

      subroutine simqit(a,r,n,nmax,alpha,dx,x,ierr)

c ----------------------------------------------------------------------
c  Iterative solution of a non-singular, real linear system of equations
c  Ax = b, where A, X and R are (n*n), (n*1) and (n*1) matrices,
c  respectively. (n*c) means: n rows, c columns.
c
c  Parameter A (input; changed): Two index array of size 1:NMAX >= N
c         containing the elements of the matrix A. First index counts
c         columns, last index counts rows.
c  R (input; changed): One-index array containing the elements of the
c         matrix B
c  N:     Number of unknowns.
c  NMAX:  Maximum value of the first index of A in dimension stateament of
c         calling the program.
c  ALPHA: Reccomended value: sqrt(1/N)
c  DX:    Iteration stops when the maximum absolute change of any unknown
c         within one iteration step is smaller than DX
c  X (input and result): array of size >N containing an estimated solution
c         when calling the subroutine (is no estimation is available, 0
c         is recommended), and the solution vector X after return. 
c ----------------------------------------------------------------------

      implicit none

      integer iz,is,nmax,n,ierr,iter,iz1,izpiv
      real rsave,alpha,dx,aa,alimit,averh,bmax,dxmax, dxmaxv
      real a(nmax,nmax),r(nmax),x(nmax)

c ----------------------------------------------------------------------

      ierr = 0
      iter = 0

c -- Initial values for R and (n+2)th column of A (equation numbers)

c$omp parallel
c$omp do private(is)
      do is=1,n
        a(is,n+2) = is + 0.5
      end do
c$omp end do
c$omp end parallel

c$omp parallel
c$omp do private(iz,is)
      do iz=1,n
        do is=1,n
          r(iz) = r(iz)-x(is)*a(is,iz)
        end do
      end do
c$omp end do
c$omp end parallel

      averh = alpha*2.
    3 averh = averh/2.

c -- For all rows

      do 50 iz=1,n

c -- Serach for pivot element

      bmax = 0.
      do 5 iz1=iz,n
      if(abs(a(iz,iz1)).gt.bmax) then
        izpiv = iz1
        bmax = abs(a(iz,iz1))
      end if
    5 continue

c -- Exchange of pivot and IZ rows; normalization of pivot row

      alimit = bmax*averh
      bmax = sign(bmax,a(iz,izpiv))
      if(izpiv.eq.iz) then
c$omp parallel
c$omp do private(is)
        do 9 is=1,n
    9   a(is,iz)=a(is,iz)/bmax
c$omp end do
c$omp end parallel
      else
c$omp parallel
c$omp do private(is,rsave)
        do 10 is=1,n
        rsave=a(is,izpiv)
        a(is,izpiv)=a(is,iz)
   10   a(is,iz)=rsave/bmax
c$omp end do
c$omp end parallel
      end if
      a(izpiv,n+2)=a(iz,n+2)
      rsave=r(izpiv)
      r(izpiv)=r(iz)
      r(iz)=rsave/bmax

c -- Elimination of large coefficients

c$omp parallel
c$omp do private(iz1,is,aa)
      do 20 iz1=iz+1,n
      aa=-a(iz,iz1)
      if(abs(aa).gt.alimit) then
        do 15 is=1,n
   15   a(is,iz1)=a(is,iz1)+aa*a(is,iz)
        r(iz1)   =r(iz1)   +aa*r(iz)
      end if
   20 continue
c$omp end do
c$omp end parallel
   50 continue

c -- Determination of changes of X stored in A(IZ,N+1)

   60 dxmaxv=1.E+30
      dxmax = 0.

cc$omp parallel
cc$omp do private(iz) reduction(max:dxmax)
      do 80 iz=n,1,-1
      a(iz,n+1)=r(iz)
      do 75 is=iz+1,n
   75 a(iz,n+1)=a(iz,n+1)-a(is,iz)*a(is,n+1)
      dxmax=max(dxmax,abs(a(iz,n+1)))
   80 x(iz)=x(iz)+a(iz,n+1)
cc$omp end do
cc$omp end parallel

      if(dxmax.lt.dx) go to 200

c -- Update of absolute value vector R

c$omp parallel
c$omp do private(iz,is)
      do 100 iz=1,n
      r(iz)=0.
      do 100 is=1,iz-1
  100 r(iz)=r(iz)-a(is,n+1)*a(is,iz)
c$omp end do
c$omp end parallel

      if(dxmax.gt.0.8*dxmaxv) go to 3
      iter = iter + 1
      if(iter.gt.15) go to 200
      dxmaxv=dxmax
      go to 60

c ----------------------------------------------------------------------
  200 return
      end
c ----------------------------------------------------------------------
