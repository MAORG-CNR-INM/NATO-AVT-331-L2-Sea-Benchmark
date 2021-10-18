c-----------------------------------------------------------------------
c CALCOLA LE COORDINATE DELLA CARENA CON ASSETTO IMPOSTO (SNKG,TRIM)
c ATTENZIONE: le coordinate sono gia' adimens. con stre_le, e la carena
c e' centrata  nel sistema di riferimento.
c-----------------------------------------------------------------------
      subroutine ASSETTO(appox,appoy,appoz,fsdom,trdom,x_0,y_0,lat_0,
     &                   lontra,itra,ktra,x1,x2,x3,x4,xn,
     &                   y1,y2,y3,y4,yn,z1,z2,z3,z4,zn,
     &                   areamin,xtot,ndom,shiplen,iclose)
c-----------------------------------------------------------------------

      parameter(ksez=200,klin=200)
      dimension appox(ksez,klin),appoy(ksez,klin),appoz(ksez,klin),
     &          jzero(ksez)
      dimension x1(*),x2(*),x3(*),x4(*),y1(*),y2(*),y3(*),y4(*),itra(*),
     &          z1(*),z2(*),z3(*),z4(*),xn(*),yn(*),zn(*),x_0(*),y_0(*) 
      real xdot(1),ydot(1),zdot(1)
      real profilox(200),profiloz(200)
      integer itriqua,ktra
      logical trovato,secco,bagnato,prua
      character*1 trdom,fsdom
      common  /appo/ lat,long,snkg,trim,xmaxfs,ymaxfs,kpan

      itriqua = 3
c-----------------------------------------------------------------------
c Inizializza il vettore jzero   
c-----------------------------------------------------------------------
      do i=1,ksez
        jzero(i)=0
      end do      
c-----------------------------------------------------------------------
c ROTO-TRASLAZIONE DEL DOMINIO
c SNKG < 0 ====> Affonda // al piano z=0
c TRIM > 0 ====> Emerge la prua, affonda la poppa
c-----------------------------------------------------------------------
      do j = 1,long+1
        do i = 1,lat+1
          ciccio=appox(i,j)
          appox(i,j)= ciccio*cos(trim)+appoz(i,j)*sin(trim) 
          appoz(i,j)=-ciccio*sin(trim)+appoz(i,j)*cos(trim)+snkg
        end do
      end do
c-----------------------------------------------------------------------
c VERIFICA DOMINI COMPLETAMENTE ASCIUTTI
c Se tutti i piu' profondi sono sopra a z=0 ==> il dominio e' emerso
c-----------------------------------------------------------------------
      secco=.true.
      do i = 1,lat+1
         if(appoz(i,long+1).lt.0.0) secco=.false.
      end do
      if(secco) then
         !write(*,*)' ------======    W A R N I N G    ======------'
         !write(*,*)'        Dominio completamente EMERSO' 
         go to 999
      end if
c-----------------------------------------------------------------------
c VERIFICA DOMINI COMPLETAMENTE IMMERSI 
c Se tutti i piu' superficiali sono sopra a z=0 ==> il dominio e' immerso
c-----------------------------------------------------------------------
      bagnato=.true.
      do i = 1,lat+1
        if(appoz(i,1).gt.0.0) bagnato=.false.
        jzero(i)=long+1      !!! Tutti i nodi sono da scartare
      end do
      if(bagnato.and.fsdom.eq.'n') then !!! Uscita regolare
        !write(*,*)' ====>>>> Dominio completamente IMMERSO' 
        do i=1,lat+1
          jzero(i)=0      !!! Non ci sono nodi da scartare
        end do 
        go to 996 
      end if
      if(bagnato.and.fsdom.eq.'y') then !!! Uscita in errore
         !write(*,*)' ------======    E R R O R E      ======------'
         !write(*,*)' Intersez. non trovata per dominio di F.S.'
         fsdom = 'n'
         do i=1,lat+1
           jzero(i)=0      !!! Non ci sono nodi da scartare
         end do 
      end if
c-----------------------------------------------------------------------
c RICERCA DELLA INTERSEZIONE CON IL PIANO Z=0 PER ASSETTO VARIABILE
c-----------------------------------------------------------------------
      iclose = 0
      do 997 i=1,lat+1
         trovato=.false.
         if(appoz(i,long+1).gt.0.0) then
           !write(*,*)'  Warning: fila',i,'  asciutta' 
           jzero(i) = long+2 !!! Per questa fila tutti i nodi vengono scartati
           if(i.eq.lat+1) iclose = 1
           go to 997
         end if
c
c-AS
         if(appoz(i,1).lt.0.0) then
c            write(*,*)'  Warning: fila',i,'  completamente immersa' 
c-AS         jzero(i) = long+2 !!! Per questa fila tutti i nodi vengono scartati
           jzero(i) = 1
c            write(*,*) long
c-AS         go to 997
         end if
c-AS
c
         j=long
         do while (.not.trovato.and.j.ge.1)  !!! Spazza tutta la sezione
           if(appoz(i,j).ge.0.0) then
             trovato=.true.          !!! Intersez. con la f.s. trovata
             jzero(i)=j                               !!! Nodi esclusi 
             xtmp =(appoz(i,j+1) * appox(i,j) 
     &            - appoz(i,j)   * appox(i,j+1))/
     &             (appoz(i,j+1) - appoz(i,j))
             ytmp =(appoz(i,j+1) * appoy(i,j) 
     &            - appoz(i,j)   * appoy(i,j+1))/
     &             (appoz(i,j+1) - appoz(i,j))
             appox(i,jzero(i)) = xtmp             !!! Nuove coordinate
             appoy(i,jzero(i)) = ytmp
             appoz(i,jzero(i)) = 0.0 
           end if
           j=j-1
         end do
c
 997  continue
c...............................
 996  continue
c-----------------------------------------------------------------------
c  E' un dominio con un lato sulla f.s.: costruisce i vettori (x_0,y_0).
c  (XMAXFS,YMAXFS) sara' l'ultimo punto della linea di galleggiamento
c-----------------------------------------------------------------------
      if(iclose.eq.1) trdom = 'n'
      if(fsdom.eq.'y') then
         !write(*,*) '-----------> E''un dominio con SUPERFICIE LIBERA'
         kk = lat_0
         do i=1,lat+1
           if(jzero(i).lt.long+1) then
             kk = kk + 1
             x_0(kk) = appox(i,jzero(i))/xtot
             y_0(kk) = appoy(i,jzero(i))/xtot
             if(xmaxfs.le.x_0(kk)) then 
               xmaxfs = x_0(kk)
               ymaxfs = y_0(kk)
             end if
           end if
         end do
         lat_0=kk
      end if
c-----------------------------------------------------------------------
c    Costruzione dei vettori dei punti di confine e di controllo
c    di una sezione, dal fondo alla sup.lib.
c-----------------------------------------------------------------------
      do 998 i = 1,lat
c-----------------------------------------------------------------------
         jsafe = max( jzero(i),jzero(i+1) ) + 1  !!! Ultimo nodo sicuro
         if(jsafe.eq.0) jsafe = 1
         if(jzero(i).gt.long+1) go to 998

         do 777 j = long,jsafe,-1

           x1a =  appox(  i,  j)/xtot
           x2a =  appox(  i,j+1)/xtot
           x3a =  appox(i+1,j+1)/xtot
           x4a =  appox(i+1,  j)/xtot

           y1a =  appoy(  i,  j)/xtot
           y2a =  appoy(  i,j+1)/xtot
           y3a =  appoy(i+1,j+1)/xtot
           y4a =  appoy(i+1,  j)/xtot

           z1a =  appoz(  i,  j)/xtot
           z2a =  appoz(  i,j+1)/xtot
           z3a =  appoz(i+1,j+1)/xtot
           z4a =  appoz(i+1,  j)/xtot

           ds1 = (x1a-x2a)**2 +
     &           (y1a-y2a)**2 +
     &           (z1a-z2a)**2
           ds2 = (x2a-x3a)**2 +
     &           (y2a-y3a)**2 +
     &           (z2a-z3a)**2
           ds3 = (x3a-x4a)**2 +
     &           (y3a-y4a)**2 +
     &           (z3a-z4a)**2
           ds4 = (x4a-x1a)**2 +
     &           (y4a-y1a)**2 +
     &           (z4a-z1a)**2
           if(ds1.ne.0) ds1 = sqrt(ds1)
           if(ds2.ne.0) ds2 = sqrt(ds2)
           if(ds3.ne.0) ds3 = sqrt(ds3)
           if(ds4.ne.0) ds4 = sqrt(ds4)

           if(ds1.lt.1.0E-4) then
             x2a = x3a
             x3a = x4a
             x4a = x1a
             y2a = y3a
             y3a = y4a
             y4a = y1a
             z2a = z3a
             z3a = z4a
             z4a = z1a
             itriqua = 3
           else if(ds2.lt.1.0E-4) then
             x3a = x4a
             x4a = x1a
             y3a = y4a
             y4a = y1a
             z3a = z4a
             z4a = z1a
             itriqua = 3
           else if(ds3.lt.1.0E-4) then
             x4a = x1a
             y4a = y1a
             z4a = z1a
             itriqua = 3
           else if(ds4.lt.1.0E-4) then
             itriqua = 3
           end if
           den = 1.0
           if(itriqua.eq.3) then
             den = ((y2a-y1a)*(z3a-z2a) -
     &              (z2a-z1a)*(y3a-y2a))**2 +
     &             ((z2a-z1a)*(x3a-x2a) -
     &              (x2a-x1a)*(z3a-z2a))**2 +
     &             ((x2a-x1a)*(y3a-y2a) -
     &              (y2a-y1a)*(x3a-x2a))**2
             den = abs(den)
             if(den.ne.0) den = sqrt(den)
           end if

c  Primo triangolo

           u1 = x4a-x1a
           u2 = y4a-y1a
           u3 = z4a-z1a
           v1 = x2a-x1a
           v2 = y2a-y1a
           v3 = z2a-z1a
           w1 = u2*v3-u3*v2
           w2 =-u1*v3+u3*v1
           w3 = u1*v2-u2*v1
           aaa = w1**2+w2**2+w3**2
           if(aaa.ne.0) aaa = sqrt(aaa)

c  Secondo triangolo

           u1  = x2a-x3a
           u2  = y2a-y3a
           u3  = z2a-z3a
           v1  = x4a-x3a
           v2  = y4a-y3a
           v3  = z4a-z3a
           w1  = u2*v3-u3*v2
           w2  =-u1*v3+u3*v1
           w3  = u1*v2-u2*v1
           bbb = w1**2+w2**2+w3**2
           if(bbb.ne.0) bbb = sqrt(bbb)
           aaa = (aaa + bbb)/2.

           if (den.gt.1.0E-7.and.aaa.gt.1.0E-7) then 

             kpan = kpan + 1

             x1(kpan) = x1a
             x2(kpan) = x2a
             x3(kpan) = x3a
             x4(kpan) = x4a

             y1(kpan) = y1a
             y2(kpan) = y2a
             y3(kpan) = y3a
             y4(kpan) = y4a

             z1(kpan) = z1a
             z2(kpan) = z2a
             z3(kpan) = z3a
             z4(kpan) = z4a

             xn(kpan) =  .25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
             yn(kpan) =  .25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
             zn(kpan) =  .25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))

             if(trdom.eq.'y'.and.i.eq.lat) then   
               lontra = lontra + 1
               itra(lontra) = kpan
             end if

             if(itriqua.ne.3) itriqua = 4

             if(itriqua.ne.3) call spezza(x1,x2,x3,x4,
     &                                    y1,y2,y3,y4,
     &                                    z1,z2,z3,z4,
     &                                    xn,yn,zn,kpan)
           end if
c
 777     continue

c-----------------------------------------------------------------------
c Se l'area dell'ultimo pannello e' piccola rispetto a quella del pannello
c precedente (sulla stessa fila) allora compatta (AREAMIN)
c-----------------------------------------------------------------------

         if(fsdom.eq.'y'.and.kpan.ne.0) then

            a1j = appox(  i,jsafe  )
            a2j = appox(  i,jsafe+1)
            a3j = appox(i+1,jsafe+1)
            a4j = appox(i+1,jsafe  )

            b1j = appoy(  i,jsafe  )
            b2j = appoy(  i,jsafe+1)
            b3j = appoy(i+1,jsafe+1)
            b4j = appoy(i+1,jsafe  )

            c1j = appoz(  i,jsafe  )
            c2j = appoz(  i,jsafe+1)
            c3j = appoz(i+1,jsafe+1)
            c4j = appoz(i+1,jsafe  )
   
            CALL COMPAREA(a1j,a2j,a3j,a4j,b1j,b2j,b3j,b4j,
     &                    c1j,c2j,c3j,c4j,area_j)

            a1j = appox(  i,jsafe-1)
            a2j = appox(  i,jsafe  )
            a3j = appox(i+1,jsafe  )
            a4j = appox(i+1,jsafe-1)

            b1j = appoy(  i,jsafe-1)
            b2j = appoy(  i,jsafe  )
            b3j = appoy(i+1,jsafe  )
            b4j = appoy(i+1,jsafe-1)

            c1j = appoz(  i,jsafe-1)
            c2j = appoz(  i,jsafe  )
            c3j = appoz(i+1,jsafe  )
            c4j = appoz(i+1,jsafe-1)

            CALL COMPAREA(a1j,a2j,a3j,a4j,b1j,b2j,b3j,b4j,
     &                    c1j,c2j,c3j,c4j,area_jm1)

            if(area_j.ne.0) then
              rappo = area_jm1/area_j
            else
              rappo = areamin
            end if

            if (rappo.gt.areamin) then
               kpan     =  kpan + 1    !! AVANZA KPAN: E'UN NUOVO PANNELLO !!
               x1(kpan) =  appox(  i,jzero(i))/xtot
               x2(kpan) =  appox(  i,jsafe)/xtot
               x3(kpan) =  appox(i+1,jsafe)/xtot
               x4(kpan) =  appox(i+1,jzero(i+1))/xtot

               y1(kpan) =  appoy(  i,jzero(i))/xtot
               y2(kpan) =  appoy(  i,jsafe)/xtot
               y3(kpan) =  appoy(i+1,jsafe)/xtot
               y4(kpan) =  appoy(i+1,jzero(i+1))/xtot

               z1(kpan) =  appoz(  i,jzero(i))/xtot
               z2(kpan) =  appoz(  i,jsafe)/xtot
               z3(kpan) =  appoz(i+1,jsafe)/xtot
               z4(kpan) =  appoz(i+1,jzero(i+1))/xtot

               xn(kpan) =  .25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
               yn(kpan) =  .25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
               zn(kpan) =  .25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))
               if(i.eq.lat.and.trdom.eq.'y') then   
                 lontra = lontra + 1
                 itra(lontra) = kpan
               end if

               write(48,*) 'ZONE'
               write(48,333) x1(kpan),y1(kpan),z1(kpan)
               write(48,333) x2(kpan),y2(kpan),z2(kpan)
               write(48,333) x3(kpan),y3(kpan),z3(kpan)
               write(48,333) x1(kpan),y1(kpan),z1(kpan)

            else                                 !!!! NON AVANZA KPAN !!!!!
c              !write(*,*) 'COMPATTA: rappo=',rappo,
c     &                   ': i=',i,': long(i)=',long-jsafe+1
              if((long-jsafe+1).lt.1) then
                kpan = kpan - 1
              else
                x1(kpan) =  appox(  i,jzero(i))/xtot
                x4(kpan) =  appox(i+1,jzero(i+1))/xtot
   
                y1(kpan) =  appoy(  i,jzero(i))/xtot
                y4(kpan) =  appoy(i+1,jzero(i+1))/xtot

                z1(kpan) =  appoz(  i,jzero(i))/xtot
                z4(kpan) =  appoz(i+1,jzero(i+1))/xtot

                xn(kpan) =  .25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
                yn(kpan) =  .25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
                zn(kpan) =  .25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))
              end if
            end if
ciiiiiiiiiiiiiiiiiiiiiiiiiiiii
         end if
ciiiiiiiiiiiiiiiiiiiiiiiiiiiii
 998  continue
c-----------------------------------------------------------------------
c    E' un dominio transom
c-----------------------------------------------------------------------
      if(trdom.eq.'y') then
        open(20,file='stern.profile',status='unknown',form='formatted')
        write(20,*) lat+1
        ndot = lat+1
        do i=1,lat+1
          write(20,*) appox(i,long+1),appoy(i,long+1),appoz(i,long+1)
        end do
        if(lontra.ge.1) ktra = 1 
      end if
c-----------------------------------------------------------------------
 999  continue
 333  format(3f12.8)
c-----------------------------------------------------------------------
      return 
      end
c-----------------------------------------------------------------------
      subroutine LIFTING(appox,appoy,appoz,xappo,yappo,zappo,
     &                nodi_i,nodi_j,dz,nlift,ksez)
c-----------------------------------------------------------------------
c  dz= altezza del pannello
c  nlift=numero di pannelli
c-----------------------------------------------------------------------
      dimension     appox(ksez,*),appoy(ksez,*),appoz(ksez,*)
      dimension     xappo(ksez,*),yappo(ksez,*),zappo(ksez,*)
c
c--- sposta i nodi per far spazio alla nuova fila
c
      do j = nodi_j,1,-1
         do i=1,nodi_i
            appox(i,j+nlift) = appox(i,j)
            appoy(i,j+nlift) = appoy(i,j)
            appoz(i,j+nlift) = appoz(i,j)
         end do
      end do
c
c--- piazza le nuove file
c
      do i=1,nodi_i
         do j=nlift,1,-1
            appox(i,j) = appox(i,j+1)
            appoy(i,j) = appoy(i,j+1)
            appoz(i,j) = appoz(i,j+1)+dz
         end do
      end do
c
      nodi_j=nodi_j+nlift
c
c--- memorizzazione xappo,y,z
c
      do i = 1,nodi_i
         do j=1,nodi_j
            xappo(i,j) = appox(i,j)
            yappo(i,j) = appoy(i,j)
            zappo(i,j) = appoz(i,j)
         end do
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine spezza(x1,x2,x3,x4,
     &                  y1,y2,y3,y4,
     &                  z1,z2,z3,z4,
     &                  xn,yn,zn,kpan)
c-----------------------------------------------------------------------
      real x1(*),x2(*),x3(*),x4(*)
      real y1(*),y2(*),y3(*),y4(*)
      real z1(*),z2(*),z3(*),z4(*)
      real xn(*),yn(*),zn(*)
c-----------------------------------------------------------------------
      v1x = x1(kpan)-x2(kpan)
      v1y = y1(kpan)-y2(kpan)
      v1z = z1(kpan)-z2(kpan)
      r = sqrt(v1x**2+v1y**2+v1z**2)
      v1x = v1x/r
      v1y = v1y/r
      v1z = v1z/r

      v2x = x3(kpan)-x2(kpan)
      v2y = y3(kpan)-y2(kpan)
      v2z = z3(kpan)-z2(kpan)
      r = sqrt(v2x**2+v2y**2+v2z**2)
      v2x = v2x/r
      v2y = v2y/r
      v2z = v2z/r

      cosxn = abs(v1y*v2z-v1z*v2y)

      if(cosxn.gt.0.25.and.z1(kpan).lt.-0.005
     &                 .or.
     &   x1(kpan).lt.-0.45.and.z1(kpan).lt.-0.005) then

        xa = x1(kpan)
        xb = x2(kpan)
        xc = x3(kpan)
        xd = x4(kpan)
        ya = y1(kpan)
        yb = y2(kpan)
        yc = y3(kpan)
        yd = y4(kpan)
        za = z1(kpan)
        zb = z2(kpan)
        zc = z3(kpan)
        zd = z4(kpan)

        x1(kpan) =  xa
        x2(kpan) = (xa+xb)/2.
        x3(kpan) = (xa+xb+xc+xd)/4.
        x4(kpan) = (xa+xd)/2.
        y1(kpan) =  ya
        y2(kpan) = (ya+yb)/2.
        y3(kpan) = (ya+yb+yc+yd)/4.
        y4(kpan) = (ya+yd)/2.
        z1(kpan) =  za
        z2(kpan) = (za+zb)/2.
        z3(kpan) = (za+zb+zc+zd)/4.
        z4(kpan) = (za+zd)/2.
        xn(kpan) =  0.25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
        yn(kpan) =  0.25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
        zn(kpan) =  0.25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))

        kpan = kpan + 1
        x1(kpan) = (xa+xd)/2.
        x2(kpan) = (xa+xb+xc+xd)/4.
        x3(kpan) = (xc+xd)/2.
        x4(kpan) =  xd
        y1(kpan) = (ya+yd)/2.
        y2(kpan) = (ya+yb+yc+yd)/4.
        y3(kpan) = (yc+yd)/2.
        y4(kpan) =  yd
        z1(kpan) = (za+zd)/2.
        z2(kpan) = (za+zb+zc+zd)/4.
        z3(kpan) = (zc+zd)/2.
        z4(kpan) =  zd
        xn(kpan) =  0.25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
        yn(kpan) =  0.25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
        zn(kpan) =  0.25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))

        kpan = kpan + 1
        x1(kpan) = (xa+xb)/2.
        x2(kpan) =  xb
        x3(kpan) = (xb+xc)/2.
        x4(kpan) = (xa+xb+xc+xd)/4.
        y1(kpan) = (ya+yb)/2.
        y2(kpan) =  yb
        y3(kpan) = (yb+yc)/2.
        y4(kpan) = (ya+yb+yc+yd)/4.
        z1(kpan) = (za+zb)/2.
        z2(kpan) =  zb
        z3(kpan) = (zb+zc)/2.
        z4(kpan) = (za+zb+zc+zd)/4.
        xn(kpan) =  0.25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
        yn(kpan) =  0.25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
        zn(kpan) =  0.25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))

        kpan = kpan + 1
        x1(kpan) = (xa+xb+xc+xd)/4.
        x2(kpan) = (xb+xc)/2.
        x3(kpan) =  xc
        x4(kpan) = (xc+xd)/2.
        y1(kpan) = (ya+yb+yc+yd)/4.
        y2(kpan) = (yb+yc)/2.
        y3(kpan) =  yc
        y4(kpan) = (yc+yd)/2.
        z1(kpan) = (za+zb+zc+zd)/4.
        z2(kpan) = (zb+zc)/2.
        z3(kpan) =  zc
        z4(kpan) = (zc+zd)/2.
        xn(kpan) =  0.25*(x1(kpan)+x2(kpan)+x3(kpan)+x4(kpan))
        yn(kpan) =  0.25*(y1(kpan)+y2(kpan)+y3(kpan)+y4(kpan))
        zn(kpan) =  0.25*(z1(kpan)+z2(kpan)+z3(kpan)+z4(kpan))

      end if
c-----------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------------
