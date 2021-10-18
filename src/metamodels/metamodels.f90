!-------------------------------------------------------------------------------
!-- subroutine for Gauss wighting interpolation
!-------------------------------------------------------------------------------
subroutine gauss(ndim,ntrain,npred,train,pred,sigma)
! Input:
!   ndim      number of dimensions
!   ntrain    training set size
!   npred     number of predictions
!   train     training set
!   pred      columns 1:ndim, prediction points
!   p         tuning parameter
! Output:
!   pred      column ndim+1, predictions

integer ndim,ntrain,npred
real    train(ntrain,ndim+1),pred(npred,ndim+1),w(ntrain)
real    sigma

!$OMP PARALLEL SHARED(pred,train) PRIVATE(i,j,k,wsum,d,w,f)
!$OMP DO
do i=1,npred
  wsum=0.
  do j=1,ntrain
    d=0.
    do k=1,ndim
      d=d+(pred(i,k)-train(j,k))**2.
    enddo
    d=sqrt(d)
    w(j)=exp(-0.5*(d/sigma)**2)
    wsum=wsum+w(j)
  enddo
  !print *, 'wsum', wsum

  f=0.
  do j=1,ntrain
    f=f+w(j)*train(j,ndim+1)
  enddo
  f=f/wsum
  pred(i,ndim+1)=f
  !print *, 'f',f
enddo
!$OMP END DO
!$OMP END PARALLEL

return
end

!-------------------------------------------------------------------------------
!-- subroutine for Shepard's method (inverse distance weighting, IDW)
!-------------------------------------------------------------------------------
subroutine shepard(ndim,ntrain,npred,train,pred,p)
! Input:
!   ndim      number of dimensions
!   ntrain    training set size
!   npred     number of predictions
!   train     training set
!   pred      columns 1:ndim, prediction points
!   p         tuning parameter
! Output:
!   pred      column ndim+1, predictions

integer ndim,ntrain,npred,p
real    train(ntrain,ndim+1),pred(npred,ndim+1),w(ntrain)

do i=1,npred
  wsum=0.
  do j=1,ntrain
    d=0.
    do k=1,ndim
      d=d+(pred(i,k)-train(j,k))**2.
    enddo
    d=sqrt(d)
    !print *, 'd', d
    w(j)=1./((d**p)+1.e-16)
    wsum=wsum+w(j)
  enddo
  !print *, 'wsum', wsum

  f=0.
  do j=1,ntrain
    f=f+w(j)*train(j,ndim+1)
  enddo
  f=f/wsum
  pred(i,ndim+1)=f
  !print *, 'f',f
enddo

return
end


!-------------------------------------------------------------------------------
!-- subroutine for polyharmonic spline
!-------------------------------------------------------------------------------
subroutine phspline(ndim,ntrain,npred,train,pred,p)
! Input:
!   ndim      number of dimensions
!   ntrain    training set size
!   npred     number of predictions
!   train     training set
!   pred      columns 1:ndim, prediction points
!   p         tuning parameter
! Output:
!   pred      column ndim+1, predictions

integer ndim,ntrain,npred,p
real    train(ntrain,ndim+1),pred(npred,ndim+1)
integer ipiv(ntrain)
real    a(ntrain,ntrain),c(ntrain)

caux =1.-abs(mod(p,2))

do i=1,ntrain
  do j=1,ntrain
    d=0.
    do k=1,ndim
      d=d+(train(i,k)-train(j,k))**2.
    enddo
    d=sqrt(d)+1.e-16
    if (caux==0.)   a(i,j)=d**p
    if (caux.ne.0.) a(i,j)=d**p*log(d)
  enddo
enddo
do i=1,ntrain
  c(i)=train(i,ndim+1)
enddo

call sgesv(ntrain,1,a,ntrain,ipiv,c,ntrain,info)
if (info.ne.0) print *, info

do i=1,npred
  f=0.
  do j=1,ntrain
    d=0.
    do k=1,ndim
      d=d+(pred(i,k)-train(j,k))**2.
    enddo
    d=sqrt(d)+1.e-16
    if (caux==0.)   f=f+c(j)*(d**p)
    if (caux.ne.0.) f=f+c(j)*(d**p*log(d))
  enddo
  pred(i,ndim+1)=f
enddo

return
end


!-------------------------------------------------------------------------------
!-- subroutine for radial basis functions network
!-------------------------------------------------------------------------------
subroutine rbf(ndim,ntrain,npred,train,pred,irbf,eps)
! Input:
!   ndim      number of dimensions
!   ntrain    training set size
!   npred     number of predictions
!   train     training set
!   pred      columns 1:ndim, prediction points
!   irbf      rbf type, 1:exponential, 2:multiquadric, 3:inverse multiquadric
!   eps       tuning parameter
! Output:
!   pred      column ndim+1, predictions

integer ndim,ntrain,npred,irbf
real    train(ntrain,ndim+1),pred(npred,ndim+1)
integer ipiv(ntrain)
real    a(ntrain,ntrain),c(ntrain)

do i=1,ntrain
  do j=1,ntrain
    d=0.
    do k=1,ndim
      d=d+(train(i,k)-train(j,k))**2.
    enddo
    d=sqrt(d)
    if (irbf==1) a(i,j)=exp(  -(eps*d)**2.)    !exp
    if (irbf==2) a(i,j)=sqrt(1+(eps*d)**2.)    !multiquadric
    if (irbf==3) a(i,j)=1./sqrt(1+(eps*d)**2.) !inverse multiquadric
  enddo
enddo
do i=1,ntrain
  c(i)=train(i,ndim+1)
enddo

call sgesv(ntrain,1,a,ntrain,ipiv,c,ntrain,info)
if (info.ne.0) print *, info

do i=1,npred
  f=0.
  do j=1,ntrain
    d=0.
    do k=1,ndim
      d=d+(pred(i,k)-train(j,k))**2.
    enddo
    d=sqrt(d)
    if (irbf==1) f=f+c(j)*exp(  -(eps*d)**2.)    !exp
    if (irbf==2) f=f+c(j)*sqrt(1+(eps*d)**2.)    !multiquadric
    if (irbf==3) f=f+c(j)*1./sqrt(1+(eps*d)**2.) !inverse multiquadric
  enddo
  pred(i,ndim+1)=f
enddo

return
end


!-------------------------------------------------------------------------------
!--subroutine for least square support vector machine (LS-SVM)
!-------------------------------------------------------------------------------
subroutine lssvm(ndim,ntrain,npred,train,pred,irbf,eps)
! Input:
!   ndim      number of dimensions
!   ntrain    training set size
!   npred     number of predictions
!   train     training set
!   pred      columns 1:ndim, prediction points
!   irbf      rbf type, 1:exponential, 2:multiquadric, 3:inverse multiquadric
!   eps       tuning parameter
! Output:
!   pred      column ndim+1, predictions

integer ndim,ntrain,npred,irbf
real    train(ntrain,ndim+1),pred(npred,ndim+1)
real    xi(ndim),xj(ndim)
real    a(ntrain+1,ntrain+1),kernel_rbf,gam
real    alpha(ntrain),b(ntrain+1)
integer ipiv(ntrain+1)

gam=1.e3

naux =ntrain+1

a(1,1)     =0.0
a(1,2:naux)=1.0
a(2:naux,1)=1.0
b(1)    =0.0
do i=2,naux
   xi(1:ndim)=train(i-1,1:ndim)
   b(i)=train(i-1,ndim+1)
   do j=2,naux
      xj(1:ndim)=train(j-1,1:ndim)
      a(i,j)=kernel_rbf(ndim,xi,xj,irbf,eps)
      if (i.eq.j) then
         a(i,j)=a(i,j)+1.0/gam
      end if
   end do
end do

call sgesv(naux,1,a,naux,ipiv,b,naux,info)
if (info.ne.0) print *, info

beta=b(1)
alpha(1:naux-1)=b(2:naux)    

do i=1,npred
   xi(1:ndim)=pred(i,1:ndim)
   pred(i,ndim+1)=beta
   do j=1,ntrain
      xj(1:ndim) = train(j,1:ndim)
      pred(i,ndim+1)=pred(i,ndim+1)+alpha(j)*kernel_rbf(ndim,xi,xj,irbf,eps)
   end do
end do

return
end

!-------------------------------------------------------------------------------
function kernel_rbf(ndim,xi,xj,irbf,eps)

real    kernel_RBF
real    xi(ndim),xj(ndim)

d=0.
do k=1,ndim
  d=d+(xi(k)-xj(k))**2.
enddo
d=sqrt(d)

if (irbf==1) kernel_rbf=exp(  -(eps*d)**2.)    !exp
if (irbf==2) kernel_rbf=sqrt(1+(eps*d)**2.)    !multiquadric
if (irbf==3) kernel_rbf=1./sqrt(1+(eps*d)**2.) !inverse multiquadric

return
end



!---------------------------------------------------------------------------
! Stochastic Radial Basis Function Network (SRBFN)
!---------------------------------------------------------------------------

subroutine srbfn(nmm,ndim,ntrain,npred,train,pred,irbf,pmin,pmax,bosolve,crbf)

   integer:: nmm,npred,ntrain,i,j,k,imm,irbfn
   real::    pred(npred,ndim+1),train(ntrain,ndim+1),crbf(ntrain,nmm)
   real::    f,p,pmin,pmax
   integer:: ipiv(ntrain)
   integer,allocatable:: err(:)
   real::    A(ntrain,ntrain),c(ntrain),faux(npred,2)
   real,allocatable:: fout(:,:)
   real::    ran(ndim)
   logical boscale, bosolve

!   nmm=100                                         !Number of MC metamodels
!   boscale=.true. !.false.

   allocate(fout(npred,nmm),err(nmm))
!   call srand(0)   

!   if (boscale) then
!     ranf=maxval(train(:,ndim+1))-minval(train(:,ndim+1))
!     minf=minval(train(:,ndim+1))
!     do k=1,ndim
!       ran(k)=maxval(train(:,k))-minval(train(:,k))
!     enddo
!   else
!     ranf=1.
!     minf=0.
!     do k=1,ndim
!       ran(k)=1.
!     enddo
!   endif
    
   !Building metamodels
!---------------------------------------------------------------------------

 if (bosolve) then

   do imm=1,nmm

     write(*,*) imm
 
!     p=pmin+rand(0)*(pmax-pmin)                    !random value [pmin;pmax]   
     p=pmin+(imm-1)*(pmax-pmin)/(nmm-1)                     !LHS-MC sampling
     do i=1,ntrain
       do j=1,ntrain
         d=0.
         do k=1,ndim
!           d=d+((train(i,k)-train(j,k))/ran(k))**2.
           d=d+(train(i,k)-train(j,k))**2.
        enddo
         d=sqrt(d)
         if (irbf==1) A(i,j)=d**p                                     !power
         if (irbf==2) A(i,j)=sqrt(1+(p*d)**2.)                 !multiquadric
         if (irbf==3) A(i,j)=1./sqrt(1+(p*d)**2.)      !inverse multiquadric 
       enddo
     enddo
     do i=1,ntrain
!       c(i)=(train(i,ndim+1)-minf)/ranf
       c(i)=train(i,ndim+1)
     enddo
     
!     call SGESV((ntrain),1,A,(ntrain),IPIV,C,(ntrain),INFO) !Solving system
     call SGESV(ntrain,1,A,ntrain,IPIV,C,ntrain,INFO)        !Solving system
     if (info.ne.0) then
       err(imm)=1    
     else
       err(imm)=0
     endif

     do i=1,ntrain
       crbf(i,imm)=c(i)
!       print *, 'c_t', c(i)
     enddo

     do j=1,npred
       f=0.
       do i=1,ntrain
         d=0.
         do k=1,ndim
!           d=d+((pred(j,k)-train(i,k))/ran(k))**2.
           d=d+(pred(j,k)-train(i,k))**2.
         enddo
         d=sqrt(d)
         if (irbf==1) f=f+c(i)*(d**p)                                 !power
         if (irbf==2) f=f+c(i)*(sqrt(1+(p*d)**2.))             !multiquadric
         if (irbf==3) f=f+c(i)*(1./sqrt(1+(p*d)**2.))  !inverse multiquadric
       enddo
!       fout(j,imm)=f*ranf+minf
       fout(j,imm)=f
     enddo

   enddo
  
!--------------------------------------------------------------------------

   if (sum(err).gt.0) print *, 'SRBF: Error solving the linear system', sum(err)

 else

   do imm=1,nmm

!     p=pmin+rand(0)*(pmax-pmin)                    !random value [pmin;pmax]   
     p=pmin+(imm-1)*(pmax-pmin)/(nmm-1)                     !LHS-MC sampling
     do i=1,ntrain
       c(i)=crbf(i,imm)
 !      print *, 'c_f', c(i)
     enddo

     do j=1,npred
       f=0.
       do i=1,ntrain
         d=0.
         do k=1,ndim
!           d=d+((pred(j,k)-train(i,k))/ran(k))**2.
           d=d+(pred(j,k)-train(i,k))**2.
         enddo
         d=sqrt(d)
         if (irbf==1) f=f+c(i)*(d**p)                                 !power
         if (irbf==2) f=f+c(i)*(sqrt(1+(p*d)**2.))             !multiquadric
         if (irbf==3) f=f+c(i)*(1./sqrt(1+(p*d)**2.))  !inverse multiquadric
       enddo
!       fout(j,imm)=f*ranf+minf
       fout(j,imm)=f
     enddo
 
   enddo

  endif

  do j=1,npred
    faux(j,1)=sum(fout(j,:))/nmm               !mean value of the function 
    pred(j,ndim+1)=faux(j,1)                              !function output
    !faux(j,2)=maxval(abs(fout(j,:)-faux(j,1))) !max gap between metamodels
  enddo

return 
end subroutine
