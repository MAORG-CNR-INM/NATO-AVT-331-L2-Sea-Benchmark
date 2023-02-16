       subroutine tatra(frd,cw_tatra)
c===============================================================
c                                           (v. 1.0 del 25/1/95)
c                                           (v. 1.1 del 2/6/95)
c                                           (v. 1.2 del 5/6/95)
c
c          T A T R A -  TAglio d'onda TRAsversale
c
c  Programma per il calcolo di Rw e di Cw tramite il metodo
c  del TAGLIO D'ONDA TRASVERSALE.
c  (cfr. M.Mandarino e Eggers,Sharma & Ward). 
c
c  I valori in ingresso si suppongono gia' adimensionalizzati con Lpp. 
c  N.B. _____ LA GRIGLIA DEVE ESSERE A PASSO Y=COST _____
c  Il programma legge il  file TATRAFILE che contiene l'elenco
c  dei file TATRAxxx.dat da analizzare (uno per ogni Fr).
c  Per ciascuno dei file TATRAxxx.dat il risultato e' la media dei
c  valori di Rw e Cw ottenuti con i tagli contenuti nella zona
c  a valle della carena. L'estensione della zona valida e' 
c  richiesta in input. L'output e' un file TATOUT0xxx con il
c  Cw per ogni ascissa e un file tarta.dat con il Fr ed il Cw.
c
c  Il file di input TATRAxxx.dat deve essere cosi' strutturato:
c
c ********** due righe di dati generali ... ********************
c
c 1)  Sup. carena; Fr; Num. valori per ogni taglio; Nlo
c 2)  Ascissa primo taglio, 
c           ascissa ultimo taglio, 
c               indice primo taglio, 
c                    indice ultimo taglio
c
c *** ... e per ogni taglio trasversale ... ********************
c
c *)  Ascissa del taglio
c *)  Y(i),ETA(i),ETA_x(i).
c
c ATTENZIONE !!!!!!!!!!!!!
c *) Il campo non e' simmetrico in y (eta non e' fz. pari)
c *) I tagli devono avere y CRESCENTI.
c *) Tutti i Froude devono essere relativi alla stessa griglia
c
c===============================================================
c ----------------------------------------------------------------------
c      program tatra
c ----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)                                        
      real*8 konda,kondas,ker,k0,k
      real*4 cw_tatra,frd
      parameter (nmax=2048)                                    

      dimension y_tg(0:nmax),eta_tg(0:nmax),eta_xtg(0:nmax),
     &          y_tgold(0:nmax),eta_tgold(0:nmax),
     &          yy(nmax),eta(nmax),eta_x(nmax),
     &          cc(nmax),ss(nmax),work(nmax),
     &          cc_x(nmax),ss_x(nmax)

      real*4 yt(nmax),et(nmax),etx(nmax)
      real*4 yi(nmax),ei(nmax),eix(nmax)
      real*4 ymax,ymin,dy,froude
      complex*16 HWS(nmax),CE
      character*80 inpname
      character*10 outname
      logical scrivi
      integer iunit

      namelist/TATRA_INPUT/ xini,xfin
c ----------------------------------------------------------------------

      iunit = 37

      open(iunit,file="SBDF.nml",action="read",recl=10000)
        read(iunit,NML=TATRA_INPUT)
      close(iunit)

      pig = dacos(-1.0d0)      
      cwmax = 0.0
      fr = dble(frd)

c===============================================================
c INIZIO CICLO SUI FILE TATRAxxx.dat
c===============================================================

c 999  continue

      open (30,file='tatra.dat',status='unknown',form='formatted')
      call tinput(inpname,outname,xini,xfin,ntra,nlo,fr,sup,ktot,iend)
c      print *, iend
c      if (iend.eq.1) then 
c        return
c      endif
      frq = fr*fr
      fr4 = frq*frq

c      print *, outname
      open(21,file=outname,status='unknown')

      outname(1:10) = 'cs0000.dat'
      nufru = nint(fr*1000)
      write(outname(3:6),'(i4.4)') nufru

c      print *, outname
      open(23,file=outname,status='unknown',form='formatted')
      write(23,*) 'VARIABLES = "u" "F" "G" "mod(H)"'

      kres = 0
      cwtot= 0.0d0
      rwtot= 0.0d0

c=======================================================================
c CALCOLO PER UN TAGLIO - kta indice del taglio di ascissa xt
c I valori di ETA ed ETA_x sono scalati con Fr**2
c=======================================================================

      ktot = ktot - 1
c      do 9999 kta=1,ktot
      do kta=1,ktot

c -- Legge il taglio

        call leggetaglio(ntra,nf,frq,xt,y_tg,eta_tg,eta_xtg,nmax)

        scrivi = .false.
        if (abs(xt).ge.abs(xini).and.
     &      abs(xt).le.abs(xfin))
     &      scrivi = .true.
        if(scrivi) write (23,*) 'ZONE'

c -- Fissa il passo y ed il du

        ydelta = dabs((y_tg(nf)-y_tg(1)))/dfloat(nf-1)
        du     = 2.*dfloat(1)/(dfloat(nf-1)*ydelta)

c -------------------------------------------------------------------- 
c     Chiamo la  D.F.T.
c -------------------------------------------------------------------- 

        call dft(eta_tg ,nf,ydelta,cc  ,ss  ,y_tg,nmax)
        call dft(eta_xtg,nf,ydelta,cc_x,ss_x,y_tg,nmax)

c -------------------------------------------------------------------- 
c     Resistenza (rw) e spettro onde libere (HWS)
c -------------------------------------------------------------------- 

        rw  = 0.0d0

        do i=1,nf

          u =-ydelta+dfloat(i-1)*du

          konda = (0.5d0*(1.0d0+dsqrt(1.0d0+4.0d0*u*u)))
          kondas = dsqrt(konda)

c------------------------------------------------------ TELSTE REED ----
c -- 
c -- w = sec(theta) = 1./cos(theta) = kondas
c -- cos(theta) = 1./w = 1/kondas
c -- sin(theta)**2 = 1. - cos(theta)**2 = 1. - 1./konda
c -- 
         esse   = dsqrt((1.+dsqrt(1.+4.*u*u))/2.)
         HWS(i) = dcmplx(cc(i),cc_x(i)/esse)
         CE     = dcmplx(dcos(esse*xt),dsin(esse*xt))
         HWS(i) = 2.*HWS(i)*CE
         F   = dimag(HWS(i))
         G   = dreal(HWS(i))
         ker = du*dsqrt(1.0d0+4.0d0*u*u)/(1.d0+dsqrt(1.0d0+4.0d0*u*u))*
     &         cdabs(HWS(i))**2
c
c------------------------------------------------------ MANDARINO   ----
c  w^2 = k0*k
c  k^2 = u^2+w^2
c  w   = kondas
c  k   = sqrt(u^2+kondas^2)
c  k0  = konda/(sqrt(u^2+konda)
c
c          k   = sqrt(u*u+konda)
c          k0  = konda/k
c          ker = 2.*du*(cc(i)*cc(i)+cc_x(i)*cc_x(i)*(2.-k0/k)/konda)
c

          if (scrivi) write(23,200) u/frq,F,G,cdabs(HWS(i))

          rw = rw + ker

        end do
c-------------------------------------------------------------------------

c-- Fourier transform ??
        rw = rw/2./pig

c-- Simmetry ??
        rw = 2.*rw
c--
        cw = fr4*rw/sup
c--
        write(21,*) xt,cw
        if(xt.gt.xini) cwmax = max(cw,cwmax)
        if (abs(xt).ge.abs(xini).and.
     &      abs(xt).le.abs(xfin)) then
            kres=kres+1
            cwtot=cwtot+cw
        end if 
        close(98)
  198 format(i4,3f12.6)

c=========================
c 9999 continue
      end do
c=========================

      cwm=cwtot/dfloat(kres)
c      write(*,*) 'Xini=',xini,'  :  Xfin=',xfin
c      write(*,*) 'Cw =',cwm
      write(30,*) fr, cwm
     
      cw_tatra = real(cwm)

      close(20)
      close(21)
      close(23)

c      go to 999

c      close(20)
c      close(21)
c      close(23)
c      close(37)                                                        

  200 format(200e16.8)

c 1000 continue
      return
c ----------------------------------------------------------------------
      end                                                              
c ----------------------------------------------------------------------

      subroutine dft(vect,np,ydelta,c,s,y,ndim)

c ----------------------------------------------------------------------

      implicit real*8(a-h,o-z)                                         
      dimension vect(0:ndim),c(ndim),s(ndim),y(0:ndim)

c ----------------------------------------------------------------------

      do j=1,np                                                    
        u = (j-1.0d0)/(np*ydelta)               
        somcc = 0.0d0                                                  
        somss = 0.0d0                                                  
        do i=2,np-1                                                  
          somcc = somcc + vect(i) * dcos(u*y(i)) * ydelta             
          somss = somss + vect(i) * dsin(u*y(i)) * ydelta             
        end do                                                         
        c(j) =   somcc +  0.5d0 * ydelta * 
     .                 (vect(1) * dcos(u*y(1))+ 
     .                 vect(np) * dcos(u*y(np)))
        s(j) = - somss -  0.5d0 * ydelta * 
     .                 (vect(1) * dsin(u*y(1))+ 
     .                 vect(np) * dsin(u*y(np)))
      end do                                                           

c-----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine tinput(inpname,outname,xini,xfin,ntra,nlo,fr,sup,ktot,
     &   iend)

c ----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)                                        
      real*4 froude
      character*80 inpname
      character*10 outname
      character*4  cha

c ----------------------------------------------------------------------

c      iend=0
c   25 read(37,*,end=999,err=999) froude
      froude = fr
c      print *, froude
      ifr = nint(froude*1000)
      write(cha,199) ifr
      inpname='tatra'//cha//'.dat'
  199 format(i4.4)
c      open(20,file=inpname,status='old',form='formatted',err=25)
c      print *, inpname
      open(20,file=inpname,status='old',form='formatted')
      read (20,*) sup,fr,ntra,nlo
      read (20,*) 
      ktot = nlo

c      write(*,*) '================================================='
c      write(*,*) '              FR =',fr
c      write(*,*) '================================================='
c      write(*,*) 'NTRA =',ntra,':    NLO =',nlo
      outname(1:6) = 'tatout'
      nufru = nint(fr*1000)
      write(outname(7:10),'(i4.4)') nufru
c      write(*,*) 'File di out : ',outname

c      write(*,*) '----------------------------------------'
c      write(*,*) 'Intervallo selezionato:'
c      write(*,*) 'Ascissa iniz=',xini,':  Ascissa fin=',xfin
c      write(*,*) '----------------------------------------'

      return
c ----------------------------------------------------------------------
c999   iend=1
      end
c ----------------------------------------------------------------------

      subroutine leggetaglio(ntra,nf,frq,xt,y,eta,etax,nmax)

c ----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      parameter(mdim=2048)

      dimension y(0:nmax),eta(0:nmax),etax(0:nmax)

      real*4 yt(mdim),et(mdim),etx(mdim)
      real*4 yi(mdim),ei(mdim),exi(mdim)

c ----------------------------------------------------------------------

      izero = 0
      read (20,*) xt

      do i=1,ntra
        read (20,*) yt(i),et(i),etx(i)
      end do

      call ordina3(yt,et,etx,ntra)

      nf = ntra

      do i=1,nf
        y(i)     = dble(yt(i)/frq)
        eta(i)   = dble(et(i)/frq)
        etax(i)  = dble(etx(i))
      end do

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
