c-----------------------------------------------------------------------

      subroutine glotta(xns,yns,zns,ncar,nulib,nplib,
     &                  ntra,ntraint,ntrapoc,nlo,nlovallec)

c-----------------------------------------------------------------------

      include "warp.cmn"

      integer ncar,nulib,nplib,ntra,ntraint,nlo,nlovallec
      integer i,j,nstop,idummy
      real xappo(ntotd+1000),yappo(ntotd+1000),zappo(ntotd+1000)
      real xns(ntotd),yns(ntotd),zns(ntotd)

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     TRUCCO PER PLOT3D
c-----------------------------------------------------------------------
c
      do i=1,ncar
        xappo(i)=xns(i)
        yappo(i)=yns(i)
        zappo(i)=zns(i)
      end do
c-----------------------------------------------------------------------
c (i) e' il VECCHIO indice di scor. globale (compresa la carena)
c (j) e' il NUOVO   indice di scor. globale (compresa la carena)
c-----------------------------------------------------------------------
c
      nstop    = ncar + ntra * (nlo-nlovallec+1) 
c
c-----------------------------------------------------------------------
c ATTENZIONE : il (+1) in nstop e' dovuto allo shift a monte
c ==> l'ultima fila a valle non si plotta!!     
c-----------------------------------------------------------------------
c........................... griglia della f.s. fino al transom
c
      j = ncar
      do i=ncar+1,nstop,ntra
c
         do idummy=1,ntraint
            j=j+1
            xappo(j)=xns(i+idummy-1)
            yappo(j)=yns(i+idummy-1)
            zappo(j)=zns(i+idummy-1)
         end do
c
         do idummy=1,ntrapoc !!! File fantasma a monte e sulla carena
            j=j+1
            xappo(j)=xns(i+ntraint)
            yappo(j)=yns(i+ntraint)
            zappo(j)=zns(i+ntraint)
         end do
c
         do idummy=1,ntra-ntraint
            j=j+1
            xappo(j)=xns(i+idummy+ntraint-1)
            yappo(j)=yns(i+idummy+ntraint-1)
            zappo(j)=zns(i+idummy+ntraint-1)
         end do
      end do
c........................... griglia della f.s. dopo il transom
      k=0
      do i=nstop+1,nulib+ncar,ntra
         do idummy=1,ntraint
            j=j+1
            xappo(j)=xns(i+idummy-1)
            yappo(j)=yns(i+idummy-1)
            zappo(j)=zns(i+idummy-1)
         end do
         do idummy=1,ntrapoc
            j=j+1
            k=k+1
            xappo(j)=xns(nulib+ncar+k)
            yappo(j)=yns(nulib+ncar+k)
            zappo(j)=zns(nulib+ncar+k)
         end do
         do idummy=1,ntra-ntraint
            j=j+1
            xappo(j)=xns(i+ntraint+idummy-1)
            yappo(j)=yns(i+ntraint+idummy-1)
            zappo(j)=zns(i+ntraint+idummy-1)
         end do
      end do
c-----------------------------------------------------------------------
      open (99, file='solfr.grd',form='unformatted',status='unknown')
      write(99) 2
      write(99) ncar,1,1,ntra+ntrapoc,nlo-1,1
      write(99) (xappo(i),i=1,ncar) ,
     &          (yappo(i),i=1,ncar) ,
     &          (zappo(i),i=1,ncar)  
      write(99) (xappo(i),i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1)),
     &          (yappo(i),i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1)),
     &          (zappo(i),i=ncar+1,ncar+(ntra+ntrapoc)*(nlo-1))
      close(99)
c-----------------------------------------------------------------------
      return
      end
