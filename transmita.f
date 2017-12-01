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
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=100)
      real angle,deltam,e,transa,pi                                       ! Declaration des variables.
      real dist1,dist2,dist1m,dist2m                                      ! angle is the zenith angle
      real z_i,z_f,dx,dy,dist,taua
      integer x_i,y_i,x_f,y_f,k
      real cell_h(height),cell_th(height),anaz    
      integer zinf,zsup
      call verticalscale(cell_th,cell_h)                                  ! define vertical scale
      pi=3.1415926                                                        ! Attribution d'une valeur a la constante pi.
      e=2.718281828     
       dist1m=3000000.      
       dist2m=3000000.                                                    ! Attribution d'une valeur a la constante e.
       do k=1,height                                                          ! Trouver le niveau initial.
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
         if ((zinf.lt.1).or.(zinf.gt.height)) then
            print*,'ERREUR zinf hors limite! (2)'
            stop
         endif   
         if ((zsup.lt.1).or.(zsup.gt.height)) then
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
