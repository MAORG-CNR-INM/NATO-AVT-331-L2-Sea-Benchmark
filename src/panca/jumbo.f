c-----------------------------------------------------------------------

      subroutine JUMBOSHIP(appox,appoy,appoz,xappo,yappo,zappo,
     &                     stre_le,jumbo,nodi_i,nodi_j,ksez)

c-----------------------------------------------------------------------
c  delta_x    : passo della griglia nella sezione di taglio PRIMA 
c               dell'allungamento.
c  xjumbo     : ascissa di taglio
c  nodi_nuovi : nuovi nodi nella parte allungata
c  delta_x_new: passo della griglia DOPO l'allungamento
c  Il flag "Cut" va messo 'y' per tutti e due i blocchi.
c-----------------------------------------------------------------------

      dimension     appox(ksez,*),appoy(ksez,*),appoz(ksez,*)
      dimension     xappo(ksez,*),yappo(ksez,*),zappo(ksez,*)
      character*1   jumbo

c-----------------------------------------------------------------------

      if(stre_le.eq.0) return

c -- memorizzazione appox,y,z

      do j = 1,nodi_j
         do i=1,nodi_i
            xappo(i,j) = appox(i,j)
            yappo(i,j) = appoy(i,j)
            zappo(i,j) = appoz(i,j)
         end do
      end do

c-----------------------------------------------------------------------
c     Calcolo il delta_x dei vecchi pannelli vicini al taglio.
c     Calcolo quanti nodi_nuovi devono essere aggiunti
c     Trasla l'ascissa di taglio.
c     Ricavo il numero di sezioni (i) da inserire.
c-----------------------------------------------------------------------

      if (jumbo.eq.'d')then
         delta_x  = abs(appox(1,1)-appox(2,1)) 
         xjumbo   = appox(1,1)-0.5*stre_le
      else
         delta_x  = abs(appox(nodi_i,1)-appox(nodi_i-1,1))
         xjumbo   = appox(nodi_i,1)+0.5*stre_le
      end if

      nodi_nuovi  = int((0.5*stre_le/delta_x)+0.5)

      if (nodi_nuovi.eq.0) then
        delta_x_new = 0.0
        if (jumbo.eq.'d')then
          delta_x_new = stre_le
          nodi_nuovi  = 1
          xjumbo   = appox(1,1)-stre_le
        end if
      else
        delta_x_new = (stre_le*0.5)/float(nodi_nuovi)
      end if
c
c-----------------------------------------------------------------------
c     Inserisce i nuovi nodi
c-----------------------------------------------------------------------

c..........................
      if (jumbo.eq.'d')then  !!!!!!!! allunga dominio downstream
c..........................

        do j = 1,nodi_j
          do i = nodi_nuovi+1 , nodi_nuovi + nodi_i !!! Vecchi nodi
            xappo(i,j)=appox(i-nodi_nuovi,j)
            yappo(i,j)=appoy(i-nodi_nuovi,j)
            zappo(i,j)=appoz(i-nodi_nuovi,j)
          end do
          do i = 1 , nodi_nuovi                     !!! Nuovi nodi   
            xappo(i,j)=xjumbo+float(i-1)*delta_x_new 
            yappo(i,j)=appoy(1,j)
            zappo(i,j)=appoz(1,j)
          end do
        end do
c..........................
      else                   !!!!!!!! allunga dominio upstream
c..........................
        do j = 1,nodi_j
          inodi=0
          do i = nodi_i+nodi_nuovi,nodi_i+1,-1
            xappo(i,j)=xjumbo-float(inodi)*delta_x_new 
            yappo(i,j)=appoy(nodi_i,j)
            zappo(i,j)=appoz(nodi_i,j)
            inodi=inodi+1
          end do
        end do
c..........................
      end if
c..........................

c-----------------------------------------------------------------------
c Aggiorna il numero dei nodi del dominio
c-----------------------------------------------------------------------

      nodi_i = nodi_i + nodi_nuovi

      do j = 1,nodi_j
         do i=1,nodi_i
            appox(i,j) = xappo(i,j)
            appoy(i,j) = yappo(i,j)
            appoz(i,j) = zappo(i,j)
         end do
      end do

c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
