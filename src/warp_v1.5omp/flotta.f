      subroutine flotta(fr,cw,ncar,fi0x,fi1x,fi0y,fi1y,
     &                  fi0z,fi1z,eta,presto,
     &                  nplib,ntrapoc,nlovallec,
     &                  ntra,ntraint,nlo,nulib,
     &                  ncard,nplibd,ntotd)
c-----------------------------------------------------------------------
      integer ncard,nplib,nplibd,ntotd
      real fi0x(ntotd),fi0y(ntotd),fi0z(ntotd),
     &     fi1x(ntotd),fi1y(ntotd),fi1z(ntotd)
      real presto(ncard)
      real eta(nplibd)
c-----------------------------------------------------------------------
      real uappo(ntotd+1000)
      real vappo(ntotd+1000)
      real wappo(ntotd+1000)
      real etappo(nplibd+1000)
      character*10 fina
c-----------------------------------------------------------------------
      real fr,cw
      integer ncar,ntrapoc,nlovallec,ntra,ntraint,nlo
c-----------------------------------------------------------------------
      integer i,j,k,ncodac,nufru,nstop,nstopnew,nco,idummy
c-----------------------------------------------------------------------

      ncodac=ntrapoc*nlovallec
      nco  = (ntrapoc+1)*(nlovallec-1)
c
      fina(1:6) = 'solfr.'
      nufru = nint((fr)*1000)
      write(fina(7:10),'(i4.4)') nufru
c
      open(99,file=fina,form='unformatted',status='unknown')
         rewind(99)
         write(99) 2
         write(99) ncar,1,1,ntra+ntrapoc,nlo-1,1
c
c     CARENA
c
      write(99) 1.,1.,fr,cw
      write(99) (1.                ,i=1,ncar),
     &          (1.+fi0x(i)+fi1x(i),i=1,ncar),
     &          (   fi0y(i)+fi1y(i),i=1,ncar),
     &          (   fi0z(i)+fi1z(i),i=1,ncar),
     &          (         presto(i),i=1,ncar)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c (i) e' il VECCHIO indice di scor. globale (compresa la carena)
c (j) e' il NUOVO   indice di scor. globale (compresa la carena)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nstop    = ncar + ntra * (nlo-nlovallec+1)
      nstopnew = ncar + (ntra+ntrapoc) * (nlo-nlovallec+1)
c----------------------------------------------------------------------
c ATTENZIONE : il (+1) in nstop e nstopnew e' dovuto allo shift a monte
c----------------------------------------------------------------------
      j = ncar
      do i=ncar+1,nstop,ntra 
         do idummy=1,ntraint
            j=j+1
            uappo(j)      = 1.+fi0x(i+idummy-1)+fi1x(i+idummy-1)
            vappo(j)      =    fi0y(i+idummy-1)+fi1y(i+idummy-1)
            wappo(j)      =    fi0z(i+idummy-1)+fi1z(i+idummy-1)
            etappo(j-ncar)=    eta(i+idummy-1-ncar) 
         end do
c----
         do idummy=1,ntrapoc
            j=j+1
            uappo(j)      = 1.+fi0x(i+ntraint)+fi1x(i+ntraint)
            vappo(j)      = fi0y(i+ntraint)+fi1y(i+ntraint)
            wappo(j)      = fi0z(i+ntraint)+fi1z(i+ntraint)
            etappo(j-ncar)= eta(i+ntraint-ncar) 
         end do
c----
         do idummy=1,ntra-ntraint
            j=j+1
            uappo(j)      = 1.+fi0x(i+idummy+ntraint-1)+
     &                      fi1x(i+idummy+ntraint-1)
            vappo(j)      = fi0y(i+idummy+ntraint-1)+
     &                      fi1y(i+idummy+ntraint-1)
            wappo(j)      = fi0z(i+idummy+ntraint-1)+
     &                      fi1z(i+idummy+ntraint-1)
            etappo(j-ncar)= eta(i+idummy+ntraint-1-ncar) 
         end do
      end do
c----
      k=0
      do i=nstop+1,nplib+ncar,ntra
         do idummy=1,ntraint
            j=j+1
            uappo(j)      = 1.+fi0x(i+idummy-1)+fi1x(i+idummy-1)
            vappo(j)      = fi0y(i+idummy-1)+fi1y(i+idummy-1)
            wappo(j)      = fi0z(i+idummy-1)+fi1z(i+idummy-1)
            etappo(j-ncar)= eta(i+idummy-1-ncar) 
         end do
         do idummy=1,ntrapoc
            j=j+1
            k=k+1
            uappo(j)      = 1.+fi0x(nulib+ncar+k)+fi1x(nulib+ncar+k)
            vappo(j)      = fi0y(nulib+ncar+k)+fi1y(nulib+ncar+k)
            wappo(j)      = fi0z(nulib+ncar+k)+fi1z(nulib+ncar+k)
            etappo(j-ncar)= eta(nulib+k) 
         end do
         do idummy=1,ntra-ntraint
            j=j+1
            uappo(j)      = 1.+fi0x(i+idummy+ntraint-1)+
     &                      fi1x(i+idummy+ntraint-1)
            vappo(j)      = fi0y(i+idummy+ntraint-1)+
     &                      fi1y(i+idummy+ntraint-1)
            wappo(j)      = fi0z(i+idummy+ntraint-1)+
     &                      fi1z(i+idummy+ntraint-1)
            etappo(j-ncar)= eta(i+idummy+ntraint-1-ncar) 
         end do
      end do
c
c     SUP.LIB.
c
      write(99) 1.,1.,fr,cw
      write(99) (1.        ,i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1)),
     &          (uappo(i)  ,i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1)),
     &          (vappo(i)  ,i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1)),
     &          (wappo(i)  ,i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1)),
     &          (etappo(i) ,i=     1,     (ntra+ntrapoc)*(nlo-1))
c
      close(99)
c
c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
