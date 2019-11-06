c-----------------------------------------------------------------------
c
c=======================================================================
c Routine transmitm (Andre Morin, Alex Neron, Etienne Rousseau 2004)
c debuggee par Martin Aube 2004
c Determine la transmittance des molecules atmospheriques sur un parcours
c entre les cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
c Recoit des longueurs d'ondes en nanometre et les transforme en microns.
c La pressi doit etre en KPa
c Retourne la transmittance transm
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
      subroutine transmitm(angle,anaz,x_i,y_i,z_i,x_f,y_f,z_f,
     + dx,dy,transm,distd,tranam)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=1024)
      real angle,e,transm                                                 ! Declaration des variables.
      real lamdm,dist1,dist2,dist1m,dist2m
      real tranae                                                         ! vertical transmittance of the complete atmosphere (molecules)
      real z_i,z_f,dx,dy,distd,pi,z1,z2
      integer x_i,y_i,x_f,y_f,k
      real cell_h(height),cell_t(height),anaz    
      integer zinf,zsup
      call verticalscale(dx,cell_t,cell_h)                                ! defining vertical scale
      pi=3.1415926                                                        ! Attribution d'une valeur a la constante pi.
      e=2.71828182844                                                     ! Attribution d'une valeur a la constante e.
      lamdm=lambda/1000.                                                  ! Les equations suivantes utilisent des lambdas en microns.
       dist1m=3000000.      
       dist2m=3000000.
       do k=1,height 
         dist1=abs(z_i-cell_h(k))/cell_t(k)
         if (dist1.lt.dist1m) then
            dist1m=dist1
            zinf=k
         endif                                                            ! Trouver le niveau initial.
         dist2=abs(z_f-cell_h(k))/cell_t(k)
         if (dist2.lt.dist2m) then
            dist2m=dist2
            zsup=k
         endif  
      enddo
         if ((zinf.lt.1).or.(zinf.gt.height)) then
            print*,'ERREUR zinf hors limite! (1)',z_i,dist1,dist1m
            stop
         endif
         if ((zsup.lt.1).or.(zsup.gt.height)) then
            print*,'ERREUR zsup hors limite! (1)'
            stop
         endif   
         distd=sqrt((real(x_f-x_i)*dx)**2.+(real(y_f-y_i)*dy)**2.)  
         if (distd.lt.dx) distd=dx       
         if (abs(angle-pi/2.).lt.abs(atan(cell_t(zsup)/distd)))            ! angle under which the cell is crossed horizontally
     +   then 
            anaz=abs(anaz-real(nint(anaz/(pi/2.)))*(pi/2.))                   ! angle equivalent de projection sur l'axe x premier quadrant, nécessaire car on calcule toujours la transmittance avec deux cellules voisines sur l'axe des x
            distd=distd/sin(angle)/abs(cos(anaz))
             if (sin(angle).eq.0.) then
               print*,'ERREUR sin(angle)=0 (1a), angle=',angle
               print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
               stop
             endif 
             if (cos(anaz).eq.0.) then
               print*,'ERREUR cos(anaz)=0 (1a), anaz=',anaz
               print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
               stop
             endif 
          else                                                            ! usual case where the cell is crossed vertically 
            if (angle.ge.pi)  angle=pi               
            distd=abs((z_i-z_f)/cos(angle))
              if (cos(angle).eq.0.) then
               print*,'ERREUR cos(angle)=0 (1a), angle=',angle
               stop
              endif
         endif
         if (z_i.gt.z_f) then
            z2=z_i
            z1=z_f
         else
            z1=z_i
            z2=z_f
         endif     
      transm=1.-((tranam-1.)*(exp(-1.*z2/8000.)-exp(-1.*z1/8000.)))
      if ((tranam.lt.0.).or.(tranam.gt.1.)) then
         print*,'ERREUR avec transm',tranam,z_f,z_i,angle,
     +   lamdm,lambda
         stop
      endif
      return
      end
