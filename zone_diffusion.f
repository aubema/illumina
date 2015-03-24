c=======================================================================
c  Routine zone diffusion 2010
c
c  Determine les cellules se trouvant dans la zone de diffusion (effet)
c  des cellules (x1,y1,z1) et (x2,y2,z2) 
c  Retourne la matrice des cellules diffusantes (diffusion) ainsi que le nombre
c  de cellules diffusantes (ncellule)
c
c  
c------------------------------------------------------------------------
c   
c    Copyright (C) 2010  Martin Aube
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
       subroutine zone_diffusion(x_1,y_1,z1,x_2,y_2,z_2,dx,dy,effet,
     +         nbx,nby,alt_sol,zondif,ncell)
       integer x_1,y_1,x_2,y_2,z_2,nbx,nby,i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax
       integer zondif(3000000,4)
       real x1,y1,z1,x2,y2,z2,x0,y0,z0,alt_sol(1024,1024)
       real dx,dy,effet,dmin,aire,a,b,c,s,delta,d,deltmx
       real cell_h(50) 
      data cell_h /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,            ! Hauteur du centre de chaque niveau
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/

 10    ncell=0
       neffet=nint(effet/(dx/2.+dy/2.))+2
c calcul de position en metre
       x1=real(x_1)*dx
       y1=real(y_1)*dy
       x2=real(x_2)*dx
       y2=real(y_2)*dy
       z2=cell_h(z_2)
       if (x_1.lt.x_2) then
         imin=x_1-neffet
       else
         imin=x_2-neffet
       endif
       if (imin.lt.1) imin=1      
c  
       if (y_1.lt.y_2) then
         jmin=y_1-neffet
       else
         jmin=y_2-neffet
       endif
       if (jmin.lt.1) jmin=1 
c       
       if (x_1.gt.x_2) then
         imax=x_1+neffet
       else
         imax=x_2+neffet
       endif
       if (imax.gt.nbx) imax=nbx   
c       
       if (y_1.gt.y_2) then
         jmax=y_1+neffet
       else
         jmax=y_2+neffet
       endif
       if (jmax.gt.nby) jmax=nby     
c
       do i=imin,imax
        do j=jmin,jmax
         do k=1,50
          x0=real(i)*dx
          y0=real(j)*dy
          z0=cell_h(k)
          if (z0.gt.alt_sol(i,j)) then
           a=sqrt((x1-x0)**2.+(y1-y0)**2.+(z1-z0)**2.)                     ! voir http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
           b=sqrt((x1-x2)**2.+(y1-y2)**2.+(z1-z2)**2.)                     ! et http://mathworld.wolfram.com/TriangleArea.html
           c=sqrt((x2-x0)**2.+(y2-y0)**2.+(z2-z0)**2.)
           delta=abs(a-c)
           s=(a+b+c)/2.
           aire=sqrt(s*(s-a)*(s-b)*(s-c))
           dmin=2.*aire/b                                                  ! dmin entre la droite definie par les points 1 et 2 et le point 0
           if ((a.gt.c)) then                                              ! Cas ou le dmin pointe hors du segment 1-2 
            d=sqrt(b**2.+a**2.)                                           ! alors l'un des angles touchant b est superieur a 90deg dans ce cas
            deltmx=abs(a-d)
            if (delta.gt.deltmx) dmin=c                                 ! on prendra le plus petit cote entre c et a
           else 
            d=sqrt(b**2.+c**2.)
            deltmx=abs(c-d)
            if (delta.gt.deltmx) dmin=a                                           
           endif                                                             
           if (dmin.le.effet) then      
            ncell=ncell+1
            if (ncell.gt.3000000) then
             effet=effet*0.9
             print*,'Reducing 2nd order scat radius:',effet
             goto 10
            endif
            zondif(ncell,1)=i                                   
            zondif(ncell,2)=j 
            zondif(ncell,3)=k
c           print*,i,j,k
           endif
          endif                                                           ! fin condition au-dessus du sol
         enddo
        enddo
       enddo
c       print*,x_1,y_1,z1,x_2,y_2,z_2
c       print*,'Nb cell diffusantes=',ncell
c       stop
       
       
       return
       end 
