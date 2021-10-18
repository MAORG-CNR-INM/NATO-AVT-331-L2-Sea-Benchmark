c-----------------------------------------------------------------------
      subroutine pressure
c-----------------------------------------------------------------------
c
c Appiccica la pressione sul fine di griglia
c
c-----------------------------------------------------------------------

      integer ni(20),nj(20)
      integer ndot
      integer i,j,k,nk,nf,nfile,nufru

      real x(20,200,200),y(20,200,200),
     &     z(20,200,200),p(20,200,200),d(20,200,200),
     &     vx(20,200,200),vy(20,200,200),vz(20,200,200)
      real xc(40000),yc(40000),zc(40000),pc(40000),
     &     vcx(40000),vcy(40000),vcz(40000),sc(40000)
      real s(40000),b(40000),u(40000),v(40000),w(40000)
      real h(40000),t(40000)

      character*3 cha
      character*4 cha2
      character*10 fina
      character*11 namepres
      character*15 fina15
      character*80 namefile
c-----------------------------------------------------------------------

      eps = 1.E-6
      call system('ls gri*tru > a0d4dfies')
      namefile ='a0d4dfies'
      call nuli(namefile,nfile)
      call system('rm a0d4dfies')
      ndot = 5

c -- Legge le pressioni

      open(71,file='fr.inp',status='old')
        read(71,*) nufru
c      do 55 inufru=1,nufru
        read(71,*) fr
      close(71)

        fina(1:6) = 'psuni.'
        nufru = nint(fr*1000)
        write(fina(7:10),'(i4.4)') nufru
c        open(67,file=fina,form='formatted', 
c     &         status='old',recl=64000,err=55)
        open(67,file=fina,form='formatted',
     &         status='old',recl=64000)
        rewind(67)
        read(67,*) ncar
        do i=1,ncar
          read(67,*) xc(i),yc(i),zc(i),
     &              vcx(i),vcy(i),vcz(i),
     &              pc(i)
        end do
        close(67)
        fina15 = 'preRAW.0000.plt'
        nufru  = nint(fr*1000)
        write(fina15(8:11),'(i4.4)') nufru
        open(2,file=fina15,form='formatted', 
     &         status='unknown')
        write(2,*) 'VARIABLES = "X" "Y" "Z" "Vx" "Vy" "Vz" "P"'
        write(2,*) 'ZONE'
        do i=1,ncar
          write(2,*) xc(i),yc(i),zc(i),
     &               vcx(i),vcy(i),vcz(i),
     &               pc(i)
        end do
        close(2)

c -- Inizializza il file per Tecplot

        i = nint(fr*1000)
        write(cha2,999) i
        namepres='pre'//cha2//'.plt'
        open(68,file=namepres,form='formatted',
     &           status='unknown')
        write(68,*) 'VARIABLES =  "X" "Y" "Z" "Vx" "Vy" "Vz" "P"'

        do nf=1,nfile

          write(cha,99) nf
          namepres='gri'//cha//'.tru'
          open(1,file=namepres,form='formatted',
     &           status='old',recl=64000)

          read(1,*) 
          read(1,*) ni(nf),nj(nf)
          nk = 1
          read(1,*) ((x(nf,i,j),i=1,ni(nf)),j=1,nj(nf)),
     &              ((y(nf,i,j),i=1,ni(nf)),j=1,nj(nf)),
     &              ((z(nf,i,j),i=1,ni(nf)),j=1,nj(nf))
          close(1)

          do i=1,ni(nf)
            do j=1,nj(nf)
  
               sminx0 = 0.03
               sminy0 = 0.03
               sminz0 = 0.03

               if(i.eq.1) then
                 sminx=3.5*abs(x(nf,2,j)-x(nf,1,j))
               else
                 sminx=2.5*abs(x(nf,i,j)-x(nf,i-1,j))
               end if
  
               if(j.eq.1) then
                 sminy=2.5*abs(y(nf,i,2)-y(nf,i,1))
                 sminz=2.5*abs(z(nf,i,2)-z(nf,i,1))
               else
                 sminy=2.5*abs(y(nf,i,j)-y(nf,i,j-1))
                 sminz=2.5*abs(z(nf,i,j)-z(nf,i,j-1))
               end if
  
               sminx = max(sminx,sminx0)
               sminy = max(sminy,sminy0)
               sminz = max(sminy,sminz)
               sminz = max(sminz,sminz0)
  
               smin = 2.0*sqrt(sminx**2+sminz**2)
               smin = 1000.
  
               kk = 0
               kz = 0
               do k=1,ncar
                 dsx=abs(x(nf,i,j)-xc(k))
                 dsy=abs(y(nf,i,j)-yc(k))
                 dsz=abs(z(nf,i,j)-zc(k))
                 if(dsx.le.sminx
     &                 .and.
     &              dsz.lt.sminz) then
                   kk = kk + 1
                   s(kk) = sqrt(eps+dsx**2+dsy**2+dsz**2)
                   u(kk) = vcx(k)
                   v(kk) = vcy(k)
                   w(kk) = vcz(k)
                   b(kk) = pc(k)
                   h(kk) = z(nf,i,j)
                   if(abs(b(kk)).gt.1) s(kk) = 1000.
                 end if
               end do

               if(kk.lt.1.or.kz.eq.kk) then
                  d(nf,i,j) = 0.0
                  p(nf,i,j) = 0.0
                 vx(nf,i,j) = 0.0
                 vy(nf,i,j) = 0.0
                 vz(nf,i,j) = 0.0
                 go to 10
               end if
  
               npoints = kk
  
               call ordina7(s,u,v,w,b,t,h,npoints)
 
               peso       = 0.0
                d(nf,i,j) = 0.0
                p(nf,i,j) = 0.0
               vx(nf,i,j) = 0.0
               vy(nf,i,j) = 0.0
               vz(nf,i,j) = 0.0
  
               if (s(1).eq.0) then
                  d(nf,i,j) = t(k)
                  p(nf,i,j) = b(k)
                 vx(nf,i,j) = u(k)
                 vy(nf,i,j) = v(k)
                 vz(nf,i,j) = w(k)
                 go to 10
               end if
  
               if (s(1).gt.smin) then
                  d(nf,i,j) = 0.0
                  p(nf,i,j) = 0.0
                 vx(nf,i,j) = 0.0
                 vy(nf,i,j) = 0.0
                 vz(nf,i,j) = 0.0
                 go to 10
               end if
  
               npoints = min(npoints,ndot)
               npoints = min(npoints,3)
  
               dhmax = -1000.
               do k=1,npoints
                 dhmax = max(dhmax,h(k))
               end do
               if (dhmax.gt.0) then
                  p(nf,i,j) = 0.0
                  d(nf,i,j) = 0.0
                 vx(nf,i,j) = 0.0
                 vy(nf,i,j) = 0.0
                 vz(nf,i,j) = 0.0
                 go to 10
               end if
  
               do k=1,npoints
                 peso = peso + 1.0/s(k)
                  d(nf,i,j) =  d(nf,i,j)+t(k)/s(k)
                  p(nf,i,j) =  p(nf,i,j)+b(k)/s(k)
                 vx(nf,i,j) = vx(nf,i,j)+u(k)/s(k)
                 vy(nf,i,j) = vy(nf,i,j)+v(k)/s(k)
                 vz(nf,i,j) = vz(nf,i,j)+w(k)/s(k)
               end do

               d(nf,i,j) =  d(nf,i,j)/peso
               p(nf,i,j) =  p(nf,i,j)/peso
              vx(nf,i,j) = vx(nf,i,j)/peso
              vy(nf,i,j) = vy(nf,i,j)/peso
              vz(nf,i,j) = vz(nf,i,j)/peso

  10           continue
             end do
           end do

c -- Scrive

           write(68,*) 'ZONE I=',ni(nf),' J=',nj(nf)
           do j=1,nj(nf)
             do i=1,ni(nf)
               write(68,100) x(nf,i,j),
     &                       y(nf,i,j),
     &                       z(nf,i,j),
     &                      vx(nf,i,j),
     &                      vy(nf,i,j),
     &                      vz(nf,i,j),
     &                       p(nf,i,j)
            end do
          end do

          write(cha,99) nf
          namepres='gri'//cha//'.pre'
          open(1,file=namepres,form='formatted',
     &           status='unknown')
          write(1,*) 1
          write(1,*) ni(nf),nj(nf)
          write(1,*) ((x(nf,i,j),i=1,ni(nf)),j=1,nj(nf)),
     &               ((y(nf,i,j),i=1,ni(nf)),j=1,nj(nf)),
     &               ((z(nf,i,j),i=1,ni(nf)),j=1,nj(nf)),
     &               ((p(nf,i,j),i=1,ni(nf)),j=1,nj(nf))
          close(1)

        end do
        close(68)

c  55  continue
c      close(71)

   99 format(i3.3)
  100 format(20e16.8)
  999 format(i4.4)
c-----------------------------------------------------------------------
c      stop
      end
c-----------------------------------------------------------------------

      subroutine nuli(filen,nl)

c ----------------------------------------------------------------------
c
c  Calcola il numero di righe del file 'filen'
c  e lo restituisce in 'nl'
c
c ----------------------------------------------------------------------

      integer ios,nl
      character*80 filen

c ----------------------------------------------------------------------

      nl = 0
      open(85,file=filen,status='old',err=1000)
  11  read(85,*,iostat=ios)
      if (ios.ge.0) then
        nl=nl+1
        goto 11
      end if
      close(85)

c ----------------------------------------------------------------------
 1000 return
      end
c ----------------------------------------------------------------------

      subroutine ordina7(s,x,y,z,p,q,r,n)

c ----------------------------------------------------------------------

      integer n,i,j
      real s(*),x(*),y(*),z(*),p(*),q(*),r(*)
      real appo

c ----------------------------------------------------------------------

      do i=1,n-1
        do j=i+1,n
          if(s(j).lt.s(i)) then
            appo = s(i)
            s(i) = s(j)
            s(j) = appo
            appo = x(i)
            x(i) = x(j)
            x(j) = appo
            appo = y(i)
            y(i) = y(j)
            y(j) = appo
            appo = z(i)
            z(i) = z(j)
            z(j) = appo
            appo = p(i)
            p(i) = p(j)
            p(j) = appo
            appo = q(i)
            q(i) = q(j)
            q(j) = appo
            appo = r(i)
            r(i) = r(j)
            r(j) = appo
          end if
        end do
      end do
c ----------------------------------------------------------------------
      return
      end
c ----------------------------------------------------------------------
