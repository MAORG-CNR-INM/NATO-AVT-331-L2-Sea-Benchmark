c=======================================================================
      subroutine panca
c=======================================================================
c
c    Versione 1.0 (Dicembre '94)  -  PANnellazione della CArena
c    Versione 1.1 (Ottobre '95)   -  Anche per carene SWATH
c    Versione 1.2 (Giugno '96)    -  Con assetto imposto
c    Versione 1.3 (Settembre '96) -  Allungamento nella sez.max.
c    Versione 1.4 (Ottobre '96)   -  Elimina i nodi ripetuti
c    Versione 1.5 (Gennaio '97)   -  Parte emersa della carena "wall side"
c    Versione 1.6 (Febbraio '97)  -  Carene SWATH asimmetriche (non transom) 
c    Versione 1.7 (Giugno '97)    -  Correzione sull'assetto (ora coerente
c                                    con WARP)
c    Versione 1.8 (Agosto '97)    -  Correzione sull'allungamento (caso di
c                                    un solo pannello da aggiungere)
c    Versione 1.9 (Febbraio '98)  -  Calcolo superficie transom e DryFroude
c                                    (Vanden-Broeck)
c
c       INPUT:
c       ~~~~~~
c       La griglia puo' essere composta piu' file, FORMATTATI, uno per ogni
c       dominio. Il numero di domini e le caratteristiche di
c       ciascuno devono essere specificate in panca.inp, ed
c       in particolare il verso di lettura dei domini, tenendo 
c       presente che alla fine l'ordinamento deve essere il seguente:
c       (visto dalla normale positiva uscente dalla carena)
c
c                                                        ^ z
c                                                        |
c      i  <-----------------------------------           |
c               |    |    |    |3,1 |2,1 |1,1 | x <------O     <----- U
c               ------------------------------          /
c               |    |    |    |    |    |1,2 |        / 
c               ------------------------------        V  y
c               |    |    |    |    |    |1,3 |      
c               ------------------------------  
c           (lat,long)                        V
c                                               j
c
c       A questo scopo (i) e (j) possono essere scambiate e/o il verso
c       di lettura puo' essere invertito.
c       Il sistema di riferimento e' il seguente: 
c       (x) centrata a meta' carena, (z) > 0 verso l'alto
c
c
c       File in uscita:
c       ~~~~~~~~~~~~~~~
c       confine.grd  : coord. punti di confine della carena:    --> SWARP
c       transom.grd  : punti controllo ed inclinazione transom: --> SWARP
c       carena_0.dat : linea di galleggiamento:                 --> LITRA
c       gri---.dat   : geometria in domini strutturati (puo' essere
c                      visualizzata con GRI.COM. Con VERIFICA.COM si
c                      puo' controllare anche il verso di lettura)
c
c
c       Assetto libero
c       ~~~~~~~~~~~~~~
c       SNKG < 0 ===> Immersioni       
c       TRIM > 0 ===> Appoppamenti
c
c       Trasformazione    |  x'= x * cos(trim) + z * sin(trim)
c            di           |  y'= y
c        coordinate       |  z'= x * sin(trim) + z * cos(trim) + snkg
c
c                                       ^ z
c                                       |
c                                       |
c                                      ---|
c        x <---------------------------------------------------/
c                 |                     | V trim > 0          /
c                 -----\                V                    (_
c                       \              snkg < 0                 \
c                        \______________________________________/
c
c       
c       File di controllo della geometria (solo se geofile='y'):
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       quapa.----            : punti confine pann. quadrilatero degenere
c       tripa.----            : punti confine pannello triangolare
c       trapa.----            : punti di confine pannelli transom
c       nor_g.dat & nor_s.dat : coord. punti di controllo e componenti 
c                               delle normali ai pannelli  (pu essere 
c                               visualizzata con VEC.COM)
c=======================================================================

c      integer, parameter :: stdout=6

      parameter(ksez=200,klin=200,ktot=ksez*klin)
      character*80 fina
      character*80 carena
      character*1  geofile,swath,sym
      integer idel(ktot)
      dimension xappo(klin,klin)
      dimension yappo(klin,klin)
      dimension zappo(klin,klin)
      dimension x1(ktot),x2(ktot),x3(ktot),x4(ktot),
     .          y1(ktot),y2(ktot),y3(ktot),y4(ktot),
     .          z1(ktot),z2(ktot),z3(ktot),z4(ktot),
     .          xn(ktot),yn(ktot),zn(ktot)
      dimension x_0(10000),y_0(10000)
      dimension xh(ktot),yh(ktot),zh(ktot),
     $          cosxn(ktot),cosyn(ktot),coszn(ktot),area(ktot),
     $          appo(ktot),xtra(ktot),ytra(ktot),ztra(ktot),
     $          atr(3,3),itriqua(ktot),itra(klin),appox(200),
     $          appoy(200)
c
      real ltot,Cb,Volume
c
      common  /appo/ lat,long,snkg,trim,xmaxfs,ymaxfs,kpan

      open(unit=101,file="panca.out",form="formatted")

c-----------------------------------------------------------------------
c     Serve per cercare i triangoli 
c-----------------------------------------------------------------------

      disloca = 0.0
      xG      = 0.0
      eps     = 1.0
 1    eps     = 0.5 * eps
      upe     = 1.0 + eps
      if (upe.ne.1.0) go to 1
      eps     = eps*5.0
c
      CALL SCRIVI

!      CALL rimuovi()
      open(47,file='realgrid.plt',
     &             status='unknown',
     &             form='formatted',
     &             recl=64000)
      write(47,*) 'VARIABLES = "X" "Y" "Z"'

      open(48,file='untouchedgrid.plt',
     &             status='unknown',
     &             form='formatted',
     &             recl=64000)
      write(48,*) 'VARIABLES = "X" "Y" "Z"'

c-----------------------------------------------------------------------
c     Lettura dei file di dominio
c-----------------------------------------------------------------------
      lontra = 0
      CALL LEGGE (x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xn,yn,zn,
     &            cosxn,cosyn,coszn,lontra,carena,npan,itra,
     &            geofile,swath,yswath,shiplen,xtot,eps,ktra,sym,
     &            n_est,n_int,numfile,dist_pp,ltot,xtrasla,ztrasla)
      if(sym.eq.'y') numfile = 2*numfile
c-----------------------------------------------------------------------
c     Inizia il ciclo sui pannelli cercando prima i triangoli
c-----------------------------------------------------------------------
      npantri = 0
      npaneff = 0
      do 1000 i = 1,npan                                                    !
        area(i) = 0.0
        CALL CETRIOLO(x1(i),x2(i),x3(i),x4(i),y1(i),y2(i),y3(i),y4(i),
     &          z1(i),z2(i),z3(i),z4(i),itriqua(i),i,eps,geofile)
        if(itriqua(i).eq.3) then ! Ha trovato un quadrilatero degenere
          npantri = npantri + 1 
          npaneff = npaneff + 1 
          idel(i) = 1
          CALL CANCOSH(xn(i),yn(i),zn(i),xh(i),yh(i),zh(i),x1(i),x2(i),
     &                 x3(i),y1(i),y2(i),y3(i),z1(i),z2(i),z3(i),
     &                 cosxn(i),cosyn(i),coszn(i),area(i))
          write(47,*) 'ZONE'
          write(47,333) x1(i),y1(i),z1(i)
          write(47,333) x2(i),y2(i),z2(i)
          write(47,333) x3(i),y3(i),z3(i)
          write(47,333) x1(i),y1(i),z1(i)
          go to 1000
        end if
c-----------------------------------------------------------------------
c     minimi quadrati: piano medio                                     !
c-----------------------------------------------------------------------

        atr(1,1) = (x1(i)-xn(i))**2+(x2(i)-xn(i))**2+
     $             (x3(i)-xn(i))**2+(x4(i)-xn(i))**2
        atr(2,2) = (y1(i)-yn(i))**2+(y2(i)-yn(i))**2+
     $             (y3(i)-yn(i))**2+(y4(i)-yn(i))**2
        atr(3,3) = (z1(i)-zn(i))**2+(z2(i)-zn(i))**2+
     $             (z3(i)-zn(i))**2+(z4(i)-zn(i))**2
        atr(1,2) = (x1(i)-xn(i))*(y1(i)-yn(i))+
     $             (x2(i)-xn(i))*(y2(i)-yn(i))+
     $             (x3(i)-xn(i))*(y3(i)-yn(i))+
     $             (x4(i)-xn(i))*(y4(i)-yn(i))
        atr(1,3) = (x1(i)-xn(i))*(z1(i)-zn(i))+
     $             (x2(i)-xn(i))*(z2(i)-zn(i))+
     $             (x3(i)-xn(i))*(z3(i)-zn(i))+
     $             (x4(i)-xn(i))*(z4(i)-zn(i))
        atr(2,3) = (y1(i)-yn(i))*(z1(i)-zn(i))+
     $             (y2(i)-yn(i))*(z2(i)-zn(i))+
     $             (y3(i)-yn(i))*(z3(i)-zn(i))+
     $             (y4(i)-yn(i))*(z4(i)-zn(i))
        atr(2,1) = atr(1,2)
        atr(3,1) = atr(1,3)
        atr(3,2) = atr(2,3)
        tre = atr(1,1)*atr(2,2)*atr(3,3) + atr(1,2)*atr(2,3)*atr(3,1) +
     $        atr(1,3)*atr(2,1)*atr(3,2) - atr(3,1)*atr(2,2)*atr(1,3) -
     $        atr(3,2)*atr(2,3)*atr(1,1) - atr(3,3)*atr(2,1)*atr(1,2)
c        if(abs(tre).ge.eps) then
c        !write(*,*)'****************************************************'
c        !write(*,*) ' pannello ',i,'         det 3x3 =',tre
c        end if
c
        det_a = atr(2,2)*atr(3,3) - atr(3,2)*atr(2,3)
        det_b = atr(1,1)*atr(3,3) - atr(3,1)*atr(1,3)
        det_c = atr(1,1)*atr(2,2) - atr(2,1)*atr(1,2)
        absa = abs(det_a)
        absb = abs(det_b)
        absc = abs(det_c)
        dem = max(absa,absb,absc)
c
        if(dem.eq.absa) then
          sud = 1./det_a 
          a = 1. 
          b = (-atr(2,1)*atr(3,3) + atr(3,1)*atr(2,3))*sud   
          c = (-atr(2,2)*atr(3,1) + atr(3,2)*atr(2,1))*sud   
        else if(dem.eq.absb) then
          sud = 1./det_b 
          a = (-atr(1,2)*atr(3,3) + atr(3,2)*atr(1,3))*sud  
          b = 1.
          c = (-atr(1,1)*atr(3,2) + atr(3,1)*atr(1,2))*sud  
        else if(dem.eq.absc) then
          sud = 1.0/det_c 
          a = (-atr(1,3)*atr(2,2) + atr(2,3)*atr(1,2))*sud  
          b = (-atr(1,1)*atr(2,3) + atr(2,1)*atr(1,3))*sud  
          c = 1.
        end if
        d = -(a*xn(i)+b*yn(i)+c*zn(i))
c-----------------------------------------------------------------------
c     costruzione del pannello                                         !
c-----------------------------------------------------------------------
        deta = a*a*a + a*b*b + a*c*c 
        detb = a*a*b + b*b*b + b*c*c
        detc = - (a*a*c + b*b*c + c*c*c)
        absa = abs(deta)
        absb = abs(detb)
        absc = abs(detc)
        dem = max(absa,absb,absc)
        if(absa.eq.dem) then
          sud = 1./deta
c         
          xx1 = (-a*a*d + a*c*(c*x1(i)-a*z1(i)) 
     %              + a*b*(b*x1(i)-a*y1(i))) * sud
          yy1 = (-a*a*(b*x1(i)-a*y1(i)) + b*c*(c*x1(i)-a*z1(i))
     $          - c*c*(b*x1(i)-a*y1(i)) - a*b*d) * sud
          zz1 = (-a*a*(c*x1(i)-a*z1(i)) + b*c*(b*x1(i)-a*y1(i)) 
     %          - b*b*(c*x1(i)-a*z1(i)) - a*c*d) * sud
          xx2 = (-a*a*d + a*c*(c*x2(i)-a*z2(i)) 
     %              + a*b*(b*x2(i)-a*y2(i))) * sud
          yy2 = (-a*a*(b*x2(i)-a*y2(i)) + b*c*(c*x2(i)-a*z2(i))
     $          - c*c*(b*x2(i)-a*y2(i)) - a*b*d) * sud
          zz2 = (-a*a*(c*x2(i)-a*z2(i)) + b*c*(b*x2(i)-a*y2(i)) 
     %          - b*b*(c*x2(i)-a*z2(i)) - a*c*d) * sud
          xx3 = (-a*a*d + a*c*(c*x3(i)-a*z3(i)) 
     %              + a*b*(b*x3(i)-a*y3(i))) * sud
          yy3 = (-a*a*(b*x3(i)-a*y3(i)) + b*c*(c*x3(i)-a*z3(i))
     $          - c*c*(b*x3(i)-a*y3(i)) - a*b*d) * sud
          zz3 = (-a*a*(c*x3(i)-a*z3(i)) + b*c*(b*x3(i)-a*y3(i)) 
     %          - b*b*(c*x3(i)-a*z3(i)) - a*c*d) * sud
          xx4 = (-a*a*d + a*c*(c*x4(i)-a*z4(i)) 
     %              + a*b*(b*x4(i)-a*y4(i))) * sud
          yy4 = (-a*a*(b*x4(i)-a*y4(i)) + b*c*(c*x4(i)-a*z4(i))
     $          - c*c*(b*x4(i)-a*y4(i)) - a*b*d) * sud
          zz4 = (-a*a*(c*x4(i)-a*z4(i)) + b*c*(b*x4(i)-a*y4(i)) 
     %          - b*b*(c*x4(i)-a*z4(i)) - a*c*d) * sud
        else if(absb.eq.dem) then
          sud = 1.0/detb
          xx1 = (-d*a*b + c*(b*x1(i)-a*y1(i))*c + 
     %          (c*y1(i)-b*z1(i))*a*c + b*(b*x1(i)-a*y1(i))*b) * sud
          yy1 = (-a*(b*x1(i)-a*y1(i))*b + c*b*(c*y1(i)-b*z1(i))
     %          - b*b*d) * sud
          zz1 = (-a*a*(c*y1(i)-b*z1(i)) - d*b*c - 
     %          c*(b*x1(i)-a*y1(i))*a - (c*y1(i)-b*z1(i))*b*b) * sud
          xx2 = (-d*a*b + c*(b*x2(i)-a*y2(i))*c + 
     %          (c*y2(i)-b*z2(i))*a*c + b*(b*x2(i)-a*y2(i))*b) * sud
          yy2 = (-a*(b*x2(i)-a*y2(i))*b + c*b*(c*y2(i)-b*z2(i))
     %          - b*b*d) * sud
          zz2 = (-a*a*(c*y2(i)-b*z2(i)) - d*b*c - 
     %          c*(b*x2(i)-a*y2(i))*a - (c*y2(i)-b*z2(i))*b*b) * sud
          xx3 = (-d*a*b + c*(b*x3(i)-a*y3(i))*c + 
     %          (c*y3(i)-b*z3(i))*a*c + b*(b*x3(i)-a*y3(i))*b) * sud
          yy3 = (-a*(b*x3(i)-a*y3(i))*b + c*b*(c*y3(i)-b*z3(i))
     %          - b*b*d) * sud
          zz3 = (-a*a*(c*y3(i)-b*z3(i)) - d*b*c - 
     %          c*(b*x3(i)-a*y3(i))*a - (c*y3(i)-b*z3(i))*b*b) * sud
          xx4 = (-d*a*b + c*(b*x4(i)-a*y4(i))*c + 
     %          (c*y4(i)-b*z4(i))*a*c + b*(b*x4(i)-a*y4(i))*b) * sud
          yy4 = (-a*(b*x4(i)-a*y4(i))*b + c*b*(c*y4(i)-b*z4(i))
     %          - b*b*d) * sud
          zz4 = (-a*a*(c*y4(i)-b*z4(i)) - d*b*c - 
     %          c*(b*x4(i)-a*y4(i))*a - (c*y4(i)-b*z4(i))*b*b) * sud
        else if(absc.eq.dem) then
          sud = 1.0/detc
          xx1 = (d*c*a - b*b*(c*x1(i)-a*z1(i)) - 
     $          (c*x1(i)-a*z1(i))*c*c + a*(c*y1(i)-b*z1(i))*b) * sud
          yy1 = (- a*(c*y1(i)-b*z1(i))*a + d*b*c - 
     %          c*(c*y1(i)-b*z1(i))*c + (c*x1(i)-a*z1(i))*b*a) * sud
          zz1 = (a*c*(c*x1(i)-a*z1(i)) + 
     %          b*(c*y1(i)-b*z1(i))*c + c*c*d) * sud
          xx2 = (d*c*a - b*b*(c*x2(i)-a*z2(i)) - 
     $          (c*x2(i)-a*z2(i))*c*c + a*(c*y2(i)-b*z2(i))*b) * sud
          yy2 = (- a*(c*y2(i)-b*z2(i))*a + d*b*c - 
     %          c*(c*y2(i)-b*z2(i))*c + (c*x2(i)-a*z2(i))*b*a) * sud
          zz2 = (a*c*(c*x2(i)-a*z2(i)) + 
     %          b*(c*y2(i)-b*z2(i))*c + c*c*d) * sud
          xx3 = (d*c*a - b*b*(c*x3(i)-a*z3(i)) - 
     $          (c*x3(i)-a*z3(i))*c*c + a*(c*y3(i)-b*z3(i))*b) * sud
          yy3 = (- a*(c*y3(i)-b*z3(i))*a + d*b*c - 
     %          c*(c*y3(i)-b*z3(i))*c + (c*x3(i)-a*z3(i))*b*a) * sud
          zz3 = (a*c*(c*x3(i)-a*z3(i)) + 
     %          b*(c*y3(i)-b*z3(i))*c + c*c*d) * sud
          xx4 = (d*c*a - b*b*(c*x4(i)-a*z4(i)) - 
     $          (c*x4(i)-a*z4(i))*c*c + a*(c*y4(i)-b*z4(i))*b) * sud
          yy4 = (- a*(c*y4(i)-b*z4(i))*a + d*b*c - 
     %          c*(c*y4(i)-b*z4(i))*c + (c*x4(i)-a*z4(i))*b*a) * sud
          zz4 = (a*c*(c*x4(i)-a*z4(i)) + 
     %          b*(c*y4(i)-b*z4(i))*c + c*c*d) * sud
        end if
        x1(i) = xx1
        y1(i) = yy1
        z1(i) = zz1
        x2(i) = xx2
        y2(i) = yy2
        z2(i) = zz2
        x3(i) = xx3
        y3(i) = yy3
        z3(i) = zz3
        x4(i) = xx4
        y4(i) = yy4
        z4(i) = zz4
c
c-----------------------------------------------------------------------
c    calcolo di xh,yh,zh, area e coseni direttori                      !
c-----------------------------------------------------------------------

        a1 = x4(i)-x3(i) 
        b1 = y4(i)-y3(i) 
        c1 = z4(i)-z3(i) 
        d1 = x1(i)*(x4(i)-x3(i)) + y1(i)*(y4(i)-y3(i)) +
     $       z1(i)*(z4(i)-z3(i))
        a2 = y4(i)-y3(i) 
        b2 = x3(i)-x4(i) 
        c2 = 0.
        d2 = x3(i)*(y4(i)-y3(i)) + y3(i)*(x3(i)-x4(i))
        a3 = 0.0
        b3 = z4(i)-z3(i) 
        c3 = y3(i)-y4(i) 
        d3 = y3(i)*(z4(i)-z3(i)) + z3(i)*(y3(i)-y4(i))
        a4 = z4(i)-z3(i)
        b4 = 0.
        c4 = x3(i)-x4(i)
        d4 = x3(i)*(z4(i)-z3(i)) + z3(i)*(x3(i)-x4(i))
        d123 = a1*b2*c3 + b1*c2*a3 + c1*a2*b3 -a3*b2*c1 -
     $       b3*c2*a1 - c3*a2*b1
        d124 = a1*b2*c4 + b1*c2*a4 + c1*a2*b4 -a4*b2*c1 -
     $       b4*c2*a1 - c4*a2*b1
        d134 = a1*b3*c4 + b1*c3*a4 + c1*a3*b4 -a4*b3*c1 -
     $       b4*c3*a1 - c4*a3*b1
        abs123 = abs(d123)
        abs124 = abs(d124)
        abs134 = abs(d134)
        dem = max(abs123,abs124,abs134)
        if(abs123.eq.dem) then
          sud = 1./d123
          xh(i) = (d1*b2*c3 + b1*c2*d3 + c1*d2*b3 -d3*b2*c1 -
     $            b3*c2*d1 - c3*d2*b1) * sud 
          yh(i) = (a1*d2*c3 + d1*c2*a3 + c1*a2*d3 -a3*d2*c1 -
     $            d3*c2*a1 - c3*a2*d1) * sud 
          zh(i) = (a1*b2*d3 + b1*d2*a3 + d1*a2*b3 -a3*b2*d1 -
     $            b3*d2*a1 - d3*a2*b1) * sud
        else if(abs124.eq.dem) then
          sud = 1./d124
          xh(i) = (d1*b2*c4 + b1*c2*d4 + c1*d2*b4 -d4*b2*c1 -
     $            b4*c2*d1 - c4*d2*b1) * sud 
          yh(i) = (a1*d2*c4 + d1*c2*a4 + c1*a2*d4 -a4*d2*c1 -
     $            d4*c2*a1 - c4*a2*d1) * sud 
          zh(i) = (a1*b2*d4 + b1*d2*a4 + d1*a2*b4 -a4*b2*d1 -
     $            b4*d2*a1 - d4*a2*b1) * sud
        else if(abs134.eq.dem) then
          sud = 1./d134
          xh(i) = (d1*b3*c4 + b1*c3*d4 + c1*d3*b4 -d4*b3*c1 -
     $            b4*c3*d1 - c4*d3*b1) * sud 
          yh(i) = (a1*d3*c4 + d1*c3*a4 + c1*a3*d4 -a4*d3*c1 -
     $            d4*c3*a1 - c4*a3*d1) * sud 
          zh(i) = (a1*b3*d4 + b1*d3*a4 + d1*a3*b4 -a4*b3*d1 -
     $            b4*d3*a1 - d4*a3*b1) * sud
        end if

c  Primo triangolo

          u1  = x4(i)-x1(i)
          u2  = y4(i)-y1(i)
          u3  = z4(i)-z1(i)
          v1  = x2(i)-x1(i)
          v2  = y2(i)-y1(i)
          v3  = z2(i)-z1(i)
          w1  = u2*v3-u3*v2
          w2  =-u1*v3+u3*v1
          w3  = u1*v2-u2*v1
          yht = (y1(i)+y2(i)+y4(i))/3.!+yswath/xtot
          aaa = sqrt(w1**2+w2**2+w3**2)
          vvv = yht*abs(w2)

c  Secondo triangolo

          u1  = x4(i)-x3(i)
          u2  = y4(i)-y3(i)
          u3  = z4(i)-z3(i)
          v1  = x2(i)-x3(i)
          v2  = y2(i)-y3(i)
          v3  = z2(i)-z3(i)
          w1  = u2*v3-u3*v2
          w2  =-u1*v3+u3*v1
          w3  = u1*v2-u2*v1
          yht = (y2(i)+y3(i)+y4(i))/3.!+yswath/xtot
          aaa = (aaa + sqrt(w1**2+w2**2+w3**2))/2.
          vvv = (vvv + yht*abs(w2))/2.

          if(cosyn(i).ge.0) then
            segno =  1.
          else
            segno = -1.
          end if

          u1  = x2(i)-x1(i)
          u2  = y2(i)-y1(i)
          u3  = z2(i)-z1(i)
          elle1 = u1*u1+u2*u2+u3*u3
          u1  = x3(i)-x2(i)
          u2  = y3(i)-y2(i)
          u3  = z3(i)-z2(i)
          elle2 = u1*u1+u2*u2+u3*u3
          u1  = x4(i)-x3(i)
          u2  = y4(i)-y3(i)
          u3  = z4(i)-z3(i)
          elle3 = u1*u1+u2*u2+u3*u3
          u1  = x1(i)-x4(i)
          u2  = y1(i)-y4(i)
          u3  = z1(i)-z4(i)
          elle4 = u1*u1+u2*u2+u3*u3
          elle = elle1+elle2+elle3+elle4
          if(elle.lt.0.1) then
            npaneff = npaneff + 1
            idel(i) = 0
            disloca = disloca + vvv*segno
            xG      = xG + vvv*segno*xn(i)
            area(i) = aaa
          else
            idel(i) = 1
          end if

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       La versione attuale calcola la normale utilizzando i 
c       vertici appartenenti al corpo (NON COMPLANARI).
c       Per calcolare la normale con i vertici complanari 
c       scommentare le seguenti righe
c
c       sden = 1. / sqrt(                           ! normale esterna !
c    $ ((y4(i)-y2(i))*(z3(i)-z1(i))-(y3(i)-y1(i))*(z4(i)-z2(i)))**2 +
c    $ ((z4(i)-z2(i))*(x3(i)-x1(i))-(z3(i)-z1(i))*(x4(i)-x2(i)))**2 +
c    $ ((x4(i)-x2(i))*(y3(i)-y1(i))-(x3(i)-x1(i))*(y4(i)-y2(i)))**2 )
c       cosxn(i) = + sden*
c    $ ((y4(i)-y2(i))*(z3(i)-z1(i))-(y3(i)-y1(i))*(z4(i)-z2(i)))
c       cosyn(i) = + sden*
c    $ ((z4(i)-z2(i))*(x3(i)-x1(i))-(z3(i)-z1(i))*(x4(i)-x2(i)))
c       coszn(i) = + sden*
c    $ ((x4(i)-x2(i))*(y3(i)-y1(i))-(x3(i)-x1(i))*(y4(i)-y2(i)))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c-----------------------------------------------------------------------
 1000 continue                                                         ! 
c-----------------------------------------------------------------------
c E' una carena SWATH ?
c-------------------------------------------------------------------

      if(swath.eq.'y'.and.sym.eq.'y') then
c
        nswath  = 2
c
c--- Raddoppia i vettori della geometria, cambiando segno
c    alle y ed al cosyn
c
        do i=1,npan
          x1(i+npan) = x1(i)
          x2(i+npan) = x2(i)
          x3(i+npan) = x3(i)
          x4(i+npan) = x4(i)
          xn(i+npan) = xn(i)
          xh(i+npan) = xh(i)
c
          y1(i+npan) = -y1(i)
          y2(i+npan) = -y2(i)
          y3(i+npan) = -y3(i)
          y4(i+npan) = -y4(i)
          yn(i+npan) = -yn(i)
          yh(i+npan) = -yh(i)
c
          z1(i+npan) = z1(i)
          z2(i+npan) = z2(i)
          z3(i+npan) = z3(i)
          z4(i+npan) = z4(i)
          zn(i+npan) = zn(i)
          zh(i+npan) = zh(i)
c
          cosxn(i+npan)   =  cosxn(i)
          cosyn(i+npan)   = -cosyn(i)
          coszn(i+npan)   =  coszn(i)
          itriqua(i+npan) =  itriqua(i)
          area(i+npan)    =  area(i)

cc  Primo triangolo
c
c          u1  = x4(i)-x1(i)
c          u3  = z4(i)-z1(i)
c          v1  = x2(i)-x1(i)
c          v3  = z2(i)-z1(i)
c          w2  =-u1*v3+u3*v1
c          yht = (y1(i)+y2(i)+y4(i))/3.!+yswath/xtot
c          vvv = yht*abs(w2)
c
cc  Secondo triangolo
c
c          u1  = x4(i)-x3(i)
c          u3  = z4(i)-z3(i)
c          v1  = x2(i)-x3(i)
c          v3  = z2(i)-z3(i)
c          w2  =-u1*v3+u3*v1
c          yht = (y2(i)+y3(i)+y4(i))/3.!+yswath/xtot
c          vvv = (vvv + yht*w2)/2.
c
c          if(cosyn(i+npan).ge.0) then
c            segno =  1.
c          else
c            segno = -1.
c          end if
c
c          disloca = disloca + vvv*segno
c          xG      = xG + vvv*segno*xn(i)

        end do

        disloca = 2.*disloca

        do i=1,lontra
          idum           = itra(i)
          itra(i+lontra) = itra(i)
          itra(i)        = idum + npan
        end do
        lontra  = 2*lontra

      else  
        nswath = 1
      end if

c -- Traslazione in corrispondenza del baricentro

      xG = xG/disloca

cc -- Riscrive le griglie
c
      !write(*,*) 'Posizione longitudinale del baricentro: ',xG
c      zero = 0.0
c      uno  = 1.0
c      do ib=1,numfile
c        call leggigriglia(xappo,yappo,zappo,ni,nj,ib)
c        do i=1,ni
c          do j=1,nj
c            xappo(i,j) = xappo(i,j) - xG
c          end do
c        end do
c        call scrivigriglia(xappo,yappo,zappo,ni,nj,
c     &                     ib,uno,zero)
c      end do
c
cc -- Sposta i pannelli
c
      do i=1,npan*nswath
c
c        x1(i) = x1(i) - xG
c        x2(i) = x2(i) - xG
c        x3(i) = x3(i) - xG
c        x4(i) = x4(i) - xG
c        xn(i) = xn(i) - xG
c        xh(i) = xh(i) - xG
c
        if(idel(i).eq.0) then
          write(47,*) 'ZONE'
          write(47,333) x1(i),y1(i),z1(i)
          write(47,333) x2(i),y2(i),z2(i)
          write(47,333) x3(i),y3(i),z3(i)
          write(47,333) x4(i),y4(i),z4(i)
          write(47,333) x1(i),y1(i),z1(i)
        end if

      end do
c
c-----------------------------------------------------------------------
c  Superficie
c-----------------------------------------------------------------------

      sup = 0.
      dryfr = 0.0
      do i=1,npan*nswath
        sup = sup + area(i)
      end do
      if(ktra.eq.1) then
        supt = 0.
        if(lontra.ne.0) then
          do i=1,lontra
            appox(i) = yn(itra(i))+yswath/xtot
            appoy(i) = zn(itra(i))
          end do

          tramin = 1000.
          do i=1,lontra
            if (zn(itra(i)).lt.tramin) tramin=zn(itra(i))
          end do

          xmax  = -1.e10
          xmin  =  1.e10
          do i=1,npan*nswath
            xmax_i=max(x1(i),x2(i),x3(i),x4(i))
            xmin_i=min(x1(i),x2(i),x3(i),x4(i))
            if(xmax_i.gt.xmax) xmax=xmax_i
            if(xmin_i.lt.xmin) xmin=xmin_i
          end do
c
          call carea(appox,appoy,lontra,supt)
          !write(*,*) 'Superficie transom =',supt*shiplen**2
c
          dryfr = 4.0*(abs(tramin)**0.5)
          !write(*,*) 'DryFroude =',dryfr
c
          if(supt.gt.0) sup = sup + supt
        end if
      end if
      disloca = 2.*disloca
      sup     = 2.*sup     !!!!!!! Area TOTALE della carena (anche y<0)
c                        !!!!!!! (si riferisce ancora al singolo scafo)
c-----------------------------------------------------------------------
c File di griglia della carena
c-----------------------------------------------------------------------

c -- Area al galleggiamento - Momento d'inerzia

      open(unit=22,file='carena_0.dat',status='old')
      rewind(22)
      read(22,*) lat
      do i=1,lat
        read(22,*) x_0(i),y_0(i)
c        x_0(i) = x_0(i) - xG
      end do
      close(22)

c      open(unit=22,file='carena_0.dat',status='old')
c      rewind(22)
c      write(22,*) lat
c      do i=1,lat
c        write(22,*) x_0(i),y_0(i)
c      end do
c      close(22)

      swl = 0.0
      wli = 0.0
      yswl = 0.
      lcf = 0.
      if(lat.ne.2) then
        if(n_int.ne.0) then
          !write(*,*) 'CATAMARANO'
          yswl = y_0(1)    ! Carena SWATH
          do i=1,lat/2-1
            xm = 0.5*(x_0(i)+x_0(i+1))
            ym = 0.5*(y_0(i)+y_0(i+1))-yswl
            dx = abs(x_0(i+1)-x_0(i))
            swl_i = abs(ym*dx)
            swl   = swl + swl_i
            wli   = wli + swl_i*xm*xm
            xcf   = lcf + swl_i*xm

            write(888,*) swl,xcf

          end do
          do i=lat/2+1,lat-1
            xm = 0.5*(x_0(i)+x_0(i+1))
            ym = 0.5*(y_0(i)+y_0(i+1))-yswl
            dx = abs(x_0(i+1)-x_0(i))
            swl_i = abs(ym*dx)
            swl   = swl + swl_i
            wli   = wli + swl_i*xm*xm
            xcf   = lcf + swl_i*xm

            write(888,*) swl,xcf

          end do
        else                  ! Monocarena
          !write(*,*) 'MONOHULL'
c - cl
          if(swath.eq.'y') yswl = y_0(1)
c - cl
          do i=1,lat-1
            xm = 0.5*(x_0(i)+x_0(i+1))
            ym = 0.5*(y_0(i)+y_0(i+1))-yswl
            dx = abs(x_0(i+1)-x_0(i))
            swl_i = abs(ym*dx)
            swl   = swl + swl_i
            wli   = wli + swl_i*xm*xm
            xcf   = lcf + swl_i*xm
          end do
        end if 
      end if 

      xcf = xcf/swl
!      print *, 'xcf=',xcf
      swl = 2.*swl
      wli = 2.*wli

      xmax  = -1.e10
      xmin  =  1.e10
      ymax  = -1.e10
      ymin  =  1.e10
      zmax  = -1.e10
      zmin  =  1.e10
c

c - as
c      write(*,*) 'npan=',npan,'npaneff=',npaneff
c - as
      open(22,file='confine.grd',status='unknown',
     &        form='formatted',recl=64000)
        rewind(22)
        write(22,100) carena
        write(22,200) (npaneff-npantri)*nswath,lontra,
     &                0,0,0,
     &                sup,shiplen,
     &                snkg,trim,ztrasla,dryfr,swl,wli,disloca
        do i=1,npan*nswath
          if(idel(i).eq.0) then
            write(22,300) x1(i),x2(i),x3(i),x4(i),
     &                    xn(i),xh(i),cosxn(i)
            xmax_i=max(x1(i),x2(i),x3(i),x4(i))
            xmin_i=min(x1(i),x2(i),x3(i),x4(i))
            if(xmax_i.gt.xmax) xmax=xmax_i
            if(xmin_i.lt.xmin) xmin=xmin_i
          end if
        end do
        do i=1,npan*nswath
          if(idel(i).eq.0) then
            write(22,300) (y1(i)+yswath/xtot),
     &                    (y2(i)+yswath/xtot),
     &                    (y3(i)+yswath/xtot),
     &                    (y4(i)+yswath/xtot),
     &                    (yn(i)+yswath/xtot),
     &                    (yh(i)+yswath/xtot),
     &                     cosyn(i)
            ymax_i=max(y1(i),y2(i),y3(i),y4(i))
            ymin_i=min(y1(i),y2(i),y3(i),y4(i))
            if(ymax_i.gt.ymax) ymax=ymax_i
            if(ymin_i.lt.ymin) ymin=ymin_i
          end if
        end do
        do i=1,npan*nswath
          if(idel(i).eq.0) then
            write(22,300) z1(i),z2(i),z3(i),z4(i),
     &                    zn(i),zh(i),coszn(i)
            zmax_i=max(z1(i),z2(i),z3(i),z4(i))
            zmin_i=min(z1(i),z2(i),z3(i),z4(i))
            if(zmax_i.gt.zmax) zmax=zmax_i
            if(zmin_i.lt.zmin) zmin=zmin_i
          end if
        end do
        do i=1,npan*nswath
          if(idel(i).eq.0) write(22,400) area(i),itriqua(i)
        end do
      close(22)

      ymax = ymax + yswath/xtot
      ymin = ymin + yswath/xtot
      supT     = sup*shiplen*shiplen
      dislocaT = disloca*shiplen**3
      Volume   = 2.*(ymax-ymin)*zmin*(xmax-xmin)
      Volume   = abs(Volume)
      Volume   = Volume*shiplen**3
      Cb       = dislocaT/Volume

      write(101,*) '=--------------------------------------------='
      write(101,*) '               Final dimensions               '
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Tot. num. of panels          =',npan*nswath
      write(101,*) 'Tot. num. of effect. panels  =',
     &              (npaneff-npantri)*nswath
      write(101,*) 'Tot. num. of triang. panles  =',npantri*nswath
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Displacement (m^3)           =',dislocaT
      write(101,*) 'Hull surface (m^2)           =',supT
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Xmax =',xmax*shiplen,': Xmin =',xmin*shiplen
      write(101,*) 'Ymax =',ymax*shiplen,': Ymin =',ymin*shiplen
      write(101,*) 'Zmax =',zmax*shiplen,': Zmin =',zmin*shiplen
      write(101,*) '=--------------------------------------------='
      write(101,*) '=--        Non-dimensional values          --='
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Volum  (non-dimensional)     =',disloca
      write(101,*) 'Surf.  (non-dimensional)     =',sup
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Figure area buoyancy         =',swl
      write(101,*) 'Bouyancy inertia moment      =',wli
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Cb                           =',Cb
      write(101,*) '=--------------------------------------------='
      write(101,*) 'Xmax =',xmax,': Xmin =',xmin
      write(101,*) 'Ymax =',ymax,': Ymin =',ymin
      write(101,*) 'Zmax =',zmax,': Zmin =',zmin
      write(101,*) '=--------------------------------------------='
c
      write(101,*) 'xG                           =',xG
      write(101,*) '=--------------------------------------------='

c
      open(unit=999,file="panca.aux",status="unknown",recl=500)
      write(999,*) shiplen
      write(999,*) disloca, dislocaT
      write(999,*) sup, supT
c attenzione chi e' sup rif???
      write(999,*) xg !equivalent to xb
      write(999,*) xcf, swl
      write(999,*) xmin, xmax
      write(999,*) ymin, ymax
      write(999,*) zmin, zmax
      close(999)

      
c-------------------------------------------------------------------
c   Se KTRA=1       =====> CARENA TRANSOM <=====
c   File di definizione del transom: 
c   i pannelli della sezione transom devono essere ordinati in senso
c   crescente verso l'esterno (y>0).
c
c itra(i) e' un vettore i cui elementi sono puntatori ai pan. transom
c itra(1) o itra(lontra) punta al pannello piu' vicino all'asse x
c
c-------------------------------------------------------------------
c.......................
      if (lontra.eq.0) ktra = 0
      if(ktra.eq.1) then
c.......................
c        Scrive i punti transom in ordine crescente verso l'esterno
c        Viene riordinato il vettore ITRA in funzione di YN
c.......................
c
        if (lontra.eq.1) then
          lontra = 2
          appo(1)=yn(itra(1))
          appo(2)=yn(itra(1))
          itra(2)=itra(1)
        else
          do i=1,lontra
            appo(i)=yn(itra(i))
          end do
          CALL IPIKSR2(lontra,appo,itra)
        end if

        open(44,file='transom.grd',status='unknown')

        write(44,*) xn(itra(1)),cosyn(itra(1))

        do i=1,lontra
          write (44,*) yn(itra(i))+yswath/xtot,zn(itra(i))
        end do

        do i=1,lontra
          dum = -cosxn(itra(i))/coszn(itra(i))
          if(abs(dum).gt.1.0) dum = 0.
          if(dum.ne.dum) dum = 0.
          if(dum.eq.dum+1) dum = 0.
          write (44,*) dum
        end do

      close(44)
c.......................
      end if
c.......................
c
      close(47)
      close(48)
 100  format(1x,a)
 200  format(5i5,100e13.5)
 300  format(7e13.5)
 333  format(3f12.8)
 400  format(1e13.5,1i5)

 999  continue

c      close(stdout)
      close(101)
      end
C________________________________________________________________
C
      subroutine carea(x,y,nl,sup)
C________________________________________________________________
C
C  Calcola l'area sottesa ad una curva. Calcola i coefficienti
C  di una spline cubica e quindi integra analiticamente le
C  curve risultanti.
C  In ingresso: ascisse ed ordinate, numero di punti.
C  Le ascisse devono essere in senso crescente. Occhio
C  al segno: se y<0, area negativa...
C________________________________________________________________
C
      real      x(200),y(200),c(603),w(603),d(2),sup
      integer   nl,i
c
      if(nl.le.1) then
        sup = 0.
        return
      end if

      sup=0.
      do i=1,nl-1
        dx = abs(x(i+1)-x(i))
        dy = abs(0.5*(y(i+1)+y(i)))
        sup=sup+dx*dy
      end do 

C________________________________________________________________
      return
      end
C________________________________________________________________
c
c-----------------------------------------------------------
c
      subroutine spln1 (n,x,y,j,d,c,w)                                 
c     this subroutine fit a set of n points exactly with an
c     integrating function made of (n-1) cubics.       
c     the cubic between x(i-1) and x(i) must agree with the cubic      
c               between points x(i) and x(i+1) at the point x(i) in the
c               first and second derivatives.                          
c     spln1  determines the coefficiente of the cubics.               
c                                                                      
c     -----------------------------------------------------------------
c     argument definition                                              
c        n      is the number of data points.fortran integer (input).  
c                                                                      
c        x      is  a one-dimensional array of n independent variable  
c               points.the elements of x must be monotonically
c               increasing .floating-point array (input).
c                                                                      
c        y      is a one dimensional array of n dependent variable      
c               points corresponding to the x array.floating-point
c               array (input).                                               
c                                                                      
c        j      =1 if the first derivatives y:(x) are given;           
c               =2 if the second derivatives y@(x) are given;          
c               fortran integer (input).                               
c                                                                      
c        d      is an array of two elements:                           
c               d(1) is the jth derivative at x(1);                    
c               d(z) is the jth derivative at x(n).                    
c               floating-point array (input).                          
c                                                                      
c        c      is a temporary storage array of (3n-3) celles which hol
c               the coefficients of the interpolating cubics.          
c               floating-point array.                                  
c                                                                      
c        w      is a temporary storage array of (3n-3) celles.         
c               floatting-point array.                                 
c     -----------------------------------------------------------------
c     restrictions                                                     
c      there must be at least two data points.                         
c     -----------------------------------------------------------------
c     other subprograms required                                       
c        none.                                                         
c     -----------------------------------------------------------------
      implicit real*4 (a-h,o-z)
      real*4 x(n),y(n),d(2),c(3*n-3),w(3*n-3)                       
c     -----------------------------------------------------------------
c             over the interval x(i) to x(i+1), the interpolating      
c             polynomial                                               
c                  y=y(i)+a(i)*z+b(i)*z**2+e(i)*z**3                   
c             where z=(x-x(i))/(x(i+1)-x(i))                           
c             is used. the coefficients a(i),b(i) and e(i) are computed
c             by spln1 and stored in locations c(3*i-2),c(3*i-1) and   
c             c(3*i) respectively.                                     
c             while working in the ith interval,the variable q will    
c             represent q=x(i+1) - x(i), and y(i) will represent       
c             y(i+1)-y(i)                                              
c     -----------------------------------------------------------------
c                                                                      
      q=x(2) - x(1)                                                    
      yi =y(2) - y(1)                                                  
      if (j.eq.2) go to 100                                            
c     -----------------------------------------------------------------
c             if the first derivative at the end points is given,      
c             a(1) is known, and the second equation becomes           
c             merely b(1)+e(1)=yi - q*d(1).                            
c     -----------------------------------------------------------------
      c(1)=q*d(1)                                                      
      c(2)=1.0
      w(2)=yi-c(1)                                                     
      go to 200                                                        
c     -----------------------------------------------------------------
c             if the second derivative at the end points is given      
c             b(1) is known, the second equation becomes               
c             a(1)+e(1)=yi-0.5*q*q*d(1). during the solution of        
c             the 3n-4 equations,a1 will be kept in cell c(2)          
c             instead of c(1) to retain the tridiagonal form of the    
c             coefficient matrix.                                      
c     -----------------------------------------------------------------
100   c(2)=0.0
      w(2)=0.5*q*q*d(1)                                                
200   m=n-2                                                            
      if(m.le.0) go to 350                                             
c     -----------------------------------------------------------------
c             upper triangularization of the tridiagonal system of     
c             equations for the coefficient matrix follows--           
c     -----------------------------------------------------------------
      do 300 i=1,m                                                     
      ai=q                                                             
      q=x(i+2)- x(i+1)                                                 
      h=ai/q                                                           
      c(3*i)=-h/(2.0-c(3*i-1))                                         
      w(3*i)=(-yi-w(3*i-1))/(2.0 - c(3*i-1))                           
      c(3*i+1)=-h*h/(h-c(3*i))                                         
      w(3*i+1)=(yi-w(3*i))/(h-c(3*i))                                  
      yi=y(i+2)- y(i+1)                                                
      c(3*i+2)=1.0/(1.0-c(3*i+1))                                      
300   w(3*i+2)=(yi-w(3*i+1))/(1.0-c(3*i+1))                            
c     -----------------------------------------------------------------
c             e(n-1) is determined directly from the last equation     
c             obtained above, and the first or second derivative       
c             value given at the end point.                            
c     -----------------------------------------------------------------
350   if(j.eq.1) go to 400                                             
      c(3*n-3)=(q*q*d(2)/2.0-w(3*n-4))/(3.0- c(3*n-4))                 
      go to 500                                                        
400   c(3*n-3)=(q*d(2)-yi-w(3*n-4))/(2.0-c(3*n-4))                     
500   m=3*n-6                                                          
      if(m.le.0) go to 700                                             
c     -----------------------------------------------------------------
c             back solution for all coefficents except                 
c             a(1) and b(1) follows--                                  
c     -----------------------------------------------------------------
      do 600 ii=1,m                                                    
      i=m-ii+3                                                         
600   c(i)=w(i)-c(i)*c(i+1)                                            
700   if(j.eq.1) go to 800                                             
c     -----------------------------------------------------------------
c             if the second derivative is given at the end points,     
c             a(1) can now be computed from the known values of        
c             b(1) and e(1). then a(1) and b(1) are put into their     
c             proper places in the c array.                            
c     -----------------------------------------------------------------
      c(1)=y(2) - y(1)-w(2)-c(3)                                       
      c(2)=w(2)                                                        
      return                                                           
800   c(2)=w(2)-c(3)    
      return                                                           
      end                                                              
