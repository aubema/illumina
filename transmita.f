c-----------------------------------------------------------------------
c
c=======================================================================
c Routine transmita (Andre Morin, Alex Neron, Etienne Rousseau 2004)
c
c Determine la transmittance des aerosols sur un parcours entre les 
c cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
c Est fonction de l'epaisseur optique des aerosols
c Recoit des longueurs d'ondes en nanometre et les transforme en microns.
c La pression doit etre en KPa
c Retourne la transmittance transa
c
c  *** J'ai valide le calcul zenith tout atm avec modtran et 
c      cela concorde M. Aubé mars 2010
c
c pour utilisation avec Illumina 
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2009  Martin Aube
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c    Contact: martin.aube@cegepsherbrooke.qc.ca
c
c
      subroutine transmita(angle,anaz,x_i,y_i,z_i,x_f,y_f,z_f,
     + dx,dy,taua,transa)
      real angle,deltam,e,transa,pi                                       ! Declaration des variables.
      real dist1,dist2,dist1m,dist2m                                      ! angle is the zenith angle
      real z_i,z_f,dx,dy,dist,taua
      integer x_i,y_i,x_f,y_f,k
      real cell_h(50),cell_th(50),anaz    
      integer zinf,zsup
      data cell_th /0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,           ! Epaisseur des niveaux.
     a 2.68,3.24,3.92,4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,
     b 25.94,31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
     c 171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,778.66,
     d 940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,3530.84,
     e 4265.26,5152.44/
                                                                          ! Matrice de la hauteur du centre de chaque niveau (metre).
      data cell_h /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,            ! Hauteur du centre de chaque niveau.
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/
      pi=3.1415926                                                        ! Attribution d'une valeur a la constante pi.
      e=2.718281828     
       dist1m=3000000.      
       dist2m=3000000.                                                    ! Attribution d'une valeur a la constante e.
       do k=1,50                                                          ! Trouver le niveau initial.
         dist1=abs(z_i-cell_h(k))/cell_th(k)
         if (dist1.lt.dist1m) then
            dist1m=dist1
            zinf=k
         endif                                                            ! Trouver le niveau final.
         dist2=abs(z_f-cell_h(k))/cell_th(k)
         if (dist2.lt.dist2m) then
            dist2m=dist2
            zsup=k
         endif  
      enddo 
         if ((zinf.lt.1).or.(zinf.gt.50)) then
            print*,'ERREUR zinf hors limite! (2)'
            stop
         endif   
         if ((zsup.lt.1).or.(zsup.gt.50)) then
           print*,'ERREUR zsup hors limite! (2)'
           stop
         endif    
         dist=sqrt((real(x_f-x_i)*dx)**2.+(real(y_f-y_i)*dy)**2.)  
         if (dist.lt.dx) dist=dx       
         if (abs(angle-pi/2.).lt.abs(atan(cell_th(zsup)/dist)))           ! angle under which the cell is crossed horizontally
     +   then 



            anaz=abs(anaz-real(nint(anaz/(pi/2.)))*(pi/2.))                   ! angle equivalent de projection sur l'axe x premier quadrant, nécessaire car on calcule toujours la transmittance avec deux cellules voisines sur l'axe des x


            deltam=(exp(-1.*cell_h(zsup)/2000.)*dist)/2000./
     +      sin(angle)/abs(cos(anaz))
             if (sin(angle).eq.0.) then
               print*,'ERREUR sin(angle)=0 (1b), angle=',angle
               print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
               stop
             endif   
             if (cos(anaz).eq.0.) then
               print*,'ERREUR cos(anaz)=0 (1b), anaz=',anaz
               print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
               stop
             endif 
          else                                                            ! usual case where the cell is crossed vertically
            if (angle.ge.pi)  angle=pi               
            deltam=abs((((exp(-1.*z_i/2000.))-(exp(-1.*z_f/2000.))))/
     +      cos(angle))
              if (cos(angle).eq.0.) then
               print*,'ERREUR cos(angle)=0 (1b), angle=',angle
               stop
              endif
         endif         
       transa=exp(-1.*deltam*taua)
       if (transa.gt.1.) then 
          print*,'ERREUR avec transa',transa,deltam,taua
          stop
       endif
      return
      end
