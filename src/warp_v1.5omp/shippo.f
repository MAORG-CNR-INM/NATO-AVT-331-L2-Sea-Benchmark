c-----------------------------------------------------------------------
      subroutine shippo(xn,yn,zn,xns,yns,zns,
     &                  ncar,nplib,ntra,nlast,
     &                  ncard,ntotd)
c-----------------------------------------------------------------------
      integer ncar,ncard,nlast,nplib,i
      real xn(ntotd),yn(ntotd),zn(ntotd)
      real xns(ntotd),yns(ntotd),zns(ntotd)
c-----------------------------------------------------------------------
c Carena - Ovviamente non c'e' shift per i punti carena
c----------------------------------------------------------------------

      do i=1,ncar
         xns(i)=xn(i)
         yns(i)=yn(i)
         zns(i)=zn(i)
      end do

c----------------------------------------------------------------------
c Prima fila sup.lib.
c----------------------------------------------------------------------

      do i = ncar+1,ncar+ntra
         xns(i)= 2.0*xn(i)-xn(i+ntra)
         yns(i)= 2.0*yn(i)-yn(i+ntra)
         zns(i)= 2.0*zn(i)-zn(i+ntra)
      end do

c----------------------------------------------------------------------
c Tutta la sup.lib. (tranne la coda e la prima fila)
c----------------------------------------------------------------------

      do i=ncar+ntra+1,nlast
         xns(i)=xn(i-ntra)
         yns(i)=yn(i-ntra)
         zns(i)=zn(i-ntra)
      end do

c----------------------------------------------------------------------
c Coda non shiftata
c----------------------------------------------------------------------

      do i=nlast+1,ncar+nplib
         xns(i)=xn(i)
         yns(i)=yn(i)
         zns(i)=zn(i)
      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
