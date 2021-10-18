c ----------------------------------------------------------------------
c
c   16 Giugno 1997 - Variazione del formato di scrittura di sinkage
c                    e trim nel file panca.inp (aumentate le cifre
c                    significative).
c
c    4 Luglio 1997 - Limitazione alle iterazioni su ogni Froude.
c
c
c   29 Ottobre 1997 - Modifica dinamica del Dry-Froude.
c
c   11 Settembre 2006 - Estrapolazione valori asintotici.
c
c ----------------------------------------------------------------------

      subroutine freewarp(igrid,ngrid,fobj,nfr)

c ----------------------------------------------------------------------

      character*1  trimoff,sinkoff,drycho
      integer      nfilest,nfilint,nfr,nitermax,k,i,conta
      integer      iunit,iobj,nobj_cw,igrid,ngrid
      
      real         fobj(nfr)
      real         sup_ref,p_din
      real         ct_all(nfr),cw_all(nfr),cf_all(nfr)
      real         ch_all(nfr),ct_ittc_all(nfr),cf_ittc_all(nfr)
      real         ct_tatra_all(nfr),cw_tatra_all(nfr) 
      real         Rt_all(nfr),Rw_all(nfr),Rf_all(nfr)
      real         Rh_all(nfr),Rt_tatra_all(nfr),Rw_tatra_all(nfr)
      real         Rt_ittc_all(nfr),Rf_ittc_all(nfr),obj_ori(nfr)
      real         sinkhyst(3),trimhyst(3)
      real, dimension(:), allocatable :: fr,fz,my,sinki,trimi

      integer      kbase,modo,nekel
      real         diana,dryfr,zG0
      real         ro,nu,gi

      
      namelist/FREE_WARP/        trimoff,sinkoff,conv,nitermax,
     &                           fr,fz,my,sinki,trimi,iobj
      namelist/WARP_PARAMETERS/  diana,kbase,modo,nekel,drycho,dryfr,zG0

      namelist/FLUID_PARAMETERS/ ro,nu,gi

c ----------------------------------------------------------------------

      allocate(fr(nfr),fz(nfr),my(nfr),sinki(nfr),trimi(nfr))

      iunit = 20

      open(iunit,file='SBDF.nml',action='read',recl=10000)
        read(iunit,NML=FREE_WARP)
      close(iunit)

      open(iunit,file='SBDF.nml',action='read',recl=10000)
        read(iunit,NML=WARP_PARAMETERS)
      close(iunit)

      open(unit=20,file='SBDF.nml',action='read',recl=10000)
        read(20,NML=FLUID_PARAMETERS)
      close(20)

      if(drycho.eq.'n') then 
        dryfr = 0.
      end if

      ct_all       = 0.
      ct_ittc_all  = 0.
      ct_tatra_all = 0.
      cw_all       = 0.
      cw_tatra_all = 0.
      ch_all       = 0.
      cf_all       = 0.
      cf_ittc_all  = 0.

      Rt_all       = 0.
      Rt_ittc_all  = 0.
      Rt_tatra_all = 0.
      Rw_all       = 0.
      Rw_tatra_all = 0.
      Rh_all       = 0.
      Rf_all       = 0.
      Rf_ittc_all  = 0.
      obj_ori       = 0.

      sinkhyst     = 0.
      trimhyst     = 0.

c ---
c      open(77,file='objective-ori.inp',status='old',recl=10000)
c        do k=1,nfr
c          read(77,*) obj_ori(k)
c        end do
c      close(77)

      obj_ori = 1.

      open(77,file='free-warp.log',form='formatted',status='unknown',
     &               recl=64000)
      write(77,*) '# It Cw Ct WSurf SinkG Trim DSink% DTrim%'

      open(78,file='free-warp.out',form='formatted',status='unknown')
      write(78,*) '# Fr Cw Ch Cf Cf_ittc Ct Ct_ittc SinkG Trim WSurf
     & Cw_tatra Ct_tatra'

      open(79,file='coef-stat.out',form='formatted',status='unknown')
      write(79,*) '# Fr Cw_pi Cw_wc Ch Cf Cf_ittc Ct_pi Ct_wc Ct_ittc 
     & SinkG Trim WSurf_s WSurf_d'

      open(80,file='force-dyn.out',form='formatted',status='unknown')
      write(80,*) '# Fr Rw_pi Rw_wc Rh Rf Rf_ittc Rt_pi Rt_wc Rt_ittc'

      call panca

      open(999,file='panca.aux',status='old')
        read(999,*) rlen
        read(999,*) rd0, rd1
        read(999,*) rs0, rs1
        read(999,*) xG
        read(999,*) xcf,swl
        read(999,*) xmin,xmax
        read(999,*) ymin,ymax
        read(999,*) zmin,zmax
      close(999) 

      write(*,*) '(Non-dimensional values)'
      write(*,*) ' '

      write(*,*) 'Static condition (Fr=0, SinkG=0, Trim=0)'
      write(*,'(A,ES14.5)') '     L0                       :', rlen
      write(*,'(A,ES14.5)') '     Displ/L0^3               :', rd0
      write(*,'(A,ES14.5)') '     WSurf/L0^2               :', rs0
      write(*,'(A,ES14.5)') '     X CoG/L0  (X CoB/L0)     :', xG
      if (swl.gt.0.) then
        write(*,'(A,ES14.5)') '     X CoF/L0                 :', xcf
      else
        write(*,'(A)') '     X CoF/L0                 :   Submerged'
      endif
      write(*,'(A,ES14.5)') '     B/L0                     :', ymax-ymin
      write(*,'(A,ES14.5)') '     T/L0                     :', -zmin
      write(*,*) ' '

c --
c Ciclo per diversi numeri di Froude
c --

      do 2000 k=1,nfr

        write(*,90) fr(k)
        write(*,*)
     &' It     Cw         Ct      WSurf    SinkG      Trim     DSink%  D
     &Trim%'
        
        conta = 0
        iconv = 0

        open(66,file='fr.inp',form='formatted',status='unknown')
          write(66,*) 1
          write(66,*) fr(k),xG,rd0,fz(k),my(k)
        close(66)

c --
c   Azzera sink and trim
c --

      sink_old = sinki(k)
      trim_old = trimi(k)
      sink     = sink_old
      atrim    = trim_old

      open(iunit,file='sinktrim.aux',status='unknown',recl=10000)
        write(iunit,*) sink, 'sinakge'
        write(iunit,*) atrim,  'trim'
      close(iunit)


c--md 

      deltam = 1.e32

c --
c   Inizia il ciclo interno
c --

      do while(conta.lt.nitermax.and.deltam.gt.conv)

        conta = conta + 1

        call panca
        call freesurface(igrid,ngrid,ntraint)
        call warp(ntraint)
        call pressure

c --
c   Verifica il nuovo assetto
c --

        call readCt(cw,ch,cf,cf_ittc,ct,ct_ittc,sink,atrim,wsurf)

c --

        if(atrim.gt.0.1)  atrim =  0.1
        if(atrim.lt.-0.1) atrim = -0.1
        if(sink.gt.0.1)   sink  =  0.1
        if(sink.lt.-0.1)  sink  = -0.1

c --

        alfa = 0.9
        beta = 0.9

        if ((sinkoff.eq.'n').or.(sinkoff.eq.'N')) then
          sink = alfa*sink  + (1.-alfa)*sink_old
        else
          sink = sink_old
        end if

        if ((trimoff.eq.'n').or.(trimoff.eq.'N')) then
          atrim = beta*atrim + (1.-beta)*trim_old
        else
          atrim = trim_old
        end if

        call checknan(sink)
        call checknan(atrim)

c - md        if (conta.eq.1) then
c          deltas = 1.e32
c        else
          deltas = 100.*abs((sink_old-sink)/(sink+1.e-12))
c        end if 

c        if (conta.eq.1) then
c          deltat = 1.e32
c        else
          deltat = 100.*abs((trim_old-atrim)/(atrim+1.e-12))
c - md       end if

        deltam = max(deltas,deltat)

        trimG = trim_old
c -as        sinkG = sink_old+(xG*sin(-trimG)+zG*cos(trimG)-zG)
        sinkG = sink_old+(xG*sin(-trimG)+zG0*cos(trimG)-zG0)

        write(*,210)  conta,cw,ct,wsurf,sinkG,trimG,deltas,deltat
        write(77,220) conta,cw,ct,wsurf,sinkG,trimG,deltas,deltat

c - as
        if((cw.lt.0.).or.(ct.lt.0.)) then
          write(*,*) '*****************************************'
          write(*,*) 'WARNING: unfeasible physic'
          write(*,*) '         the actual geometry is discarded'
          write(*,*) '*****************************************'
          fobj       = 1.e+32
          return
        elseif(isnan(cw).or.isnan(ct)) then
          write(*,*) '*****************************************'
          write(*,*) 'WARNING: unfeasible physic'
          write(*,*) '         the actual geometry is discarded'
          write(*,*) '*****************************************'
          fobj       = 1.e+32
          return
        end if
c - as

        sink_old = sink
        trim_old = atrim

c --
c   Riscrive panca.inp
c --

        open(iunit,file='sinktrim.aux',status='unknown',recl=10000)
          write(iunit,*) sink, 'sinakge'
          write(iunit,*) atrim,  'trim'
        close(iunit)

      end do

c -- Tatra

      call tatra(fr(k),cw_tatra)
      ct_tatra = cw_tatra + cf
c      write(*,310) 'TWC',cw_tatra,ct_tatra


c -- Sinkage & trim a convergenza

      write(78,100) fr(k),cw,ch,cf,cf_ittc,ct,ct_ittc,
     &              sinkG,trimG,wsurf,cw_tatra,ct_tatra

c -- Coef static

      cw_all(k)       =       cw*wsurf/rs0
      cw_tatra_all(k) = cw_tatra*wsurf/rs0
      ch_all(k)       =       ch*wsurf/rs0
      cf_all(k)       =       cf*wsurf/rs0
      cf_ittc_all(k)  =  cf_ittc*wsurf/rs0
      ct_all(k)       =       ct*wsurf/rs0
      ct_ittc_all(k)  =  ct_ittc*wsurf/rs0
      ct_tatra_all(k) = ct_tatra*wsurf/rs0

      write(79,100) fr(k),cw_all(k),cw_tatra_all(k),ch_all(k),
     &              cf_all(k),cf_ittc_all(k),
     &              ct_all(k),ct_tatra_all(k),ct_ittc_all(k),
     &              sinkG,trimG,rs0,wsurf

c--as
c -- Calcolo resistenze
c --------------------------------------------------------
      p_din   = 0.5*ro*(fr(k)**2)*gi*rlen
      sup_ref = wsurf*(rlen**2)

      Rw_all(k)       = p_din*sup_ref*cw
      Rw_tatra_all(k) = p_din*sup_ref*cw_tatra
      Rh_all(k)       = p_din*sup_ref*ch
      Rf_all(k)       = p_din*sup_ref*cf
      Rf_ittc_all(k)  = p_din*sup_ref*cf_ittc
      Rt_all(k)       = p_din*sup_ref*ct
      Rt_ittc_all(k)  = p_din*sup_ref*ct_ittc
      Rt_tatra_all(k) = p_din*sup_ref*ct_tatra

      write(80,100) fr(k),Rw_all(k),Rw_tatra_all(k),
     &              Rh_all(k),Rf_all(k),Rf_ittc_all(k),
     &              Rt_all(k),Rt_tatra_all(k),Rt_ittc_all(k)
c
c--as

      write(77,*)
      write(77,*)

      write(*,*)

 2000 continue

c--as
c -- Objective function
c ---------------------------------------------------------
      if(iobj.eq.1) then
        do k=1,nfr
           fobj(k) = (Rt_all(k)/obj_ori(k))
        end do
      end if

      if(iobj.eq.2) then
        do k=1,nfr
           fobj(k) = (Rt_tatra_all(k)/obj_ori(k))
        end do
      end if

      if(iobj.eq.3) then
        do k=1,nfr
           fobj(k) = Rt_all(k) 
        end do
      end if

      if(iobj.eq.4) then
c        fobj(:) = sum(Rt_all(:))/float(nfr)
        do k=1,nfr
           fobj(k) = Rt_all(k)
        end do
      end if

c--as

c --

      close(77)
      close(78)
      close(79)
      close(80)

  90  format('  Fr',f6.3)
 100  format(100es18.8)
 110  format(a16,100es18.8)
 200  format(i3,100es18.8)
 210  format(i4,2es11.3,f7.3,2es11.3,2(2X,f6.2))
 220  format(i3,5es18.8,2(1X,f10.6))
 300  format(a1)
 310  format(a4,2es11.3)

      deallocate(fr,fz,my,sinki,trimi)

c ----------------------------------------------------------------------
      end
c ----------------------------------------------------------------------
c***********************************************************
      subroutine stampa
c
      write(*,*)'------------------------------------------ '
      write(*,*)'       FFFFF  RRRRR   EEEEEE  EEEEEE' 
      write(*,*)'       FF     RR  RR  EE      EE' 
      write(*,*)'       FFFFE  RRRRR   EEEEEE  EEEEEE' 
      write(*,*)'       FF     RR RR   EE      EE' 
      write(*,*)'       FF     RR  RR  EEEEEE  EEEEEE' 
      write(*,*)'------------------------------------------ '
      write(*,*)'      Free model REsistance Evaluation'
c ----------------------------------------------------------------------
      return 
      end
c ----------------------------------------------------------------------

      subroutine readCt(cw,ch,cf,cf_ittc,ct,ct_ittc,sink,atrim,wsurf)

c ----------------------------------------------------------------------

      real dum,cw,b,c,d,e,f,sink,atrim
      character*120 testa

c ----------------------------------------------------------------------

      open(22,file='Ct.out',form='formatted',status='old',err=100)
      read(22,*) testa
      read(22,*,iostat=ios) dum,cw,ch,cf,cf_ittc,ct,ct_ittc,sink,atrim,
     &                          wsurf
      if(ios.ne.0) then
        cw = 0.
        ch  = 0.
        cf  = 0.
        ct  = 0.
        cf_ittc  = 0.
        ct_ittc  = 0.
c--md
        sink  = 0.
        atrim = 0.
c--md
        wsurf = 0.
      end if
 100  close(22)

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine ValutaAsintoti(nlast,vsink,vtrim,vcw,
     &                          fs,fs1,ft,ft1,fw,fw1,
     &                          dxw1,dxw2,dxs1,dxs2,dxt1,dxt2)

c ----------------------------------------------------------------------

      integer nlast
      real    fs,fs1,ft,ft1,fw,fw1,dxw1,dxw2,dxs1,dxs2,dxt1,dxt2
      real    vsink(50),vtrim(50),vcw(50)
      logical conv

c ----------------------------------------------------------------------

      itry = nlast
      conv = .false.

      do while((itry.ge.4).or.(.not.conv))

        fs  = vsink(itry)
        fs1 = vsink(itry-1)
        fs2 = vsink(itry-2)
        fs3 = vsink(itry-3)

        ft  = vtrim(itry)
        ft1 = vtrim(itry-1)
        ft2 = vtrim(itry-2)
        ft3 = vtrim(itry-3)

        fw  = vcw(itry)
        fw1 = vcw(itry-1)
        fw2 = vcw(itry-2)
        fw3 = vcw(itry-3)

        dxs0 = abs(fs3-fs2)
        dxs1 = abs(fs2-fs1)
        dxs2 = abs(fs1-fs )

        dxt0 = abs(ft3-ft2)
        dxt1 = abs(ft2-ft1)
        dxt2 = abs(ft1-ft )

        dxw0 = abs(fw3-fw2)
        dxw1 = abs(fw2-fw1)
        dxw2 = abs(fw1-fw )

        if((dxs2.gt.dxs1)
     &     .and.
     &    (dxt2.gt.dxt1)) conv = .true.

        itry = itry - 1

      end do

      if(conv) then
        c1    = 4./3.
        c2    = 1./3.
        sink1 = c1*fs - c2*fs1
        trim1 = c1*ft - c2*ft1

        sink2 = c1*fs1 - c2*fs2
        trim2 = c1*ft1 - c2*ft2

        sink3 = c1*fs2 - c2*fs3
        trim3 = c1*ft2 - c2*ft3

        deltas1 = sink2-sink1
        deltas2 = sink3-sink2

        deltat1 = trim2-trim1
        deltat2 = trim3-trim2

        conv = .false.

        if((deltas2.gt.deltas1)
     &     .and.
     &    (deltat2.gt.deltat1)) conv = .true.

      end if

      if(.not.conv) then

        fs  = vsink(nlast)
        fs1 = vsink(nlast)
        fs2 = vsink(nlast)

        ft  = vtrim(nlast)
        ft1 = vtrim(nlast)
        ft2 = vtrim(nlast)

        fw  = vcw(nlast)
        fw1 = vcw(nlast)
        fw2 = vcw(nlast)

        dxs1 = abs(fs2-fs1)
        dxs2 = abs(fs1-fs )

        dxt1 = abs(ft2-ft1)
        dxt2 = abs(ft1-ft )

        dxw1 = abs(fw2-fw1)
        dxw2 = abs(fw1-fw )

        write(*,*) 'Not converged'

      else
        write(*,*) 'Convergence at ',itry+1
      end if

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------

      subroutine checknan(x)

c ----------------------------------------------------------------------

      if(x.ne.x) x = 0.
      if(x.eq.x+1) x = 0.

c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
