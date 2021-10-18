c--------------------------------------------------------------
      subroutine LU(AV,LDA,LMAT,LICN,N,IRN,ICN,IKEEP,IDISP)
c--------------------------------------------------------------
      integer  LDA,N,LICN
      integer  IKEEP(5*LDA),IRN(LICN),ICN(LICN),
     &         IDISP(10),IW(16000)
      real*8   W(16000)
      integer  NZ,IFAIL
      real*8   AV(LMAT)
      real*8   PIVOT
      logical  LBLOCK,GROW,ABORT(4)
      external F01BRF
 
c--------------------------------------------------------------

      NZ = N*N

      k = 0
      do j=1,N
        do i=1,N
          k = k + 1
          irn(k) = i
          icn(k) = j
          !write(*,*) av(k),irn(k),icn(k)
        end do
        do i=1,LMAT-k-(LDA-N)
          av(k+i) = av(k+i+(LDA-N))
        end do
      end do

      PIVOT = 0.1
      LBLOCK = .FALSE.
      GROW = .TRUE.
      ABORT(1) = .TRUE.
      ABORT(2) = .TRUE.
      ABORT(3) = .FALSE.
      ABORT(4) = .TRUE.
      IFAIL = 110

      !write(*,*) 'Inizio F01BRF'
      call F01BRF(N,NZ,AV,LICN,IRN,LICN,ICN,PIVOT,
     &            IKEEP,IW,W,LBLOCK,GROW,
     &            ABORT,IDISP,IFAIL)
      !write(*,*) 'Fine F01BRF',IFAIL

c--------------------------------------------------------------
      return
      end
 
c--------------------------------------------------------------
      subroutine SOL(AV,LDA,LICN,LMAT,N,ICN,IKEEP,IDISP,B)
c--------------------------------------------------------------
      integer LDA,N,LICN
      integer IKEEP(5*LDA),IRN(LICN),ICN(LICN),
     &        IDISP(10),IW(8*LDA)
      integer MTYPE
      real*8  B(LDA),W(16000)
      real*8  AV(LMAT),RESID
      external F04AXF
c--------------------------------------------------------------
      MTYPE = 1
      do i=1,N
        !write(*,*) B(i)
      end do
      !write(*,*) 'Inizio F04AXF'
      call F04AXF(N,AV,LICN,ICN,IKEEP,B,W,MTYPE,
     &            IDISP, RESID)
      !write(*,*) 'Fine F04AXF'
c--------------------------------------------------------------
      return
      end
c--------------------------------------------------------------
