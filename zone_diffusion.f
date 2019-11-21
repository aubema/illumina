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
c    Copyright (C) 2019  Martin Aube
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
     +         nbx,nby,alt_sol,zondif,ncell,stepdi)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=1024)
       integer x_1,y_1,x_2,y_2,z_2,nbx,nby,i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax,kmax
       integer zondif(3000000,4),keep,stepdi
       real x1,y1,z1,x2,y2,z2,x0,y0,z0,alt_sol(width,width)
       real dx,dy,effet,dmin,aire,a,b,c,s,delta,d,deltmx
       real d1,d2
       real cell_t(height),cell_h(height)
       call verticalscale(dx,cell_t,cell_h)                                  ! define the vertical scale
 10    ncell=0
       neffet=nint(effet/(dx/2.+dy/2.))+1
       keep=0
c calcul de position en metre
       x1=real(x_1)*dx
       y1=real(y_1)*dy
       x2=real(x_2)*dx
       y2=real(y_2)*dy
       z2=cell_h(z_2)
       if (x_1.le.x_2) then
         imin=x_1-neffet
         imax=x_2+neffet
       else
         imin=x_2-neffet
         imax=x_1+neffet
       endif
       if (imin.lt.1) imin=1 
       if (imax.gt.nbx) imax=nbx
       if (y_1.le.y_2) then
         jmin=y_1-neffet
         jmax=y_2+neffet
       else
         jmin=y_2-neffet
         jmax=y_1+neffet
       endif
       if (jmin.lt.1) jmin=1
       if (jmax.gt.nby) jmax=nby
       kmax=z_2+neffet 
       do i=imin,imax
        do j=jmin,jmax
         do k=1,kmax
          x0=real(i)*dx
          y0=real(j)*dy
          z0=cell_h(k)
             d1=sqrt((x1-x0)**2.+(y1-y0)**2.+(1.-z0)**2.)
             d2=sqrt((x2-x0)**2.+(y2-y0)**2.+(z2-z0)**2.)
             if ((d1.le.effet).or.(d2.le.effet)) then 
               if (keep.eq.0) then
                  ncell=ncell+1
                  zondif(ncell,1)=i                                   
                  zondif(ncell,2)=j 
                  zondif(ncell,3)=k
               endif
               keep=keep+1
               if (keep.eq.stepdi) keep=0
             endif

         enddo
        enddo
       enddo
       return
       end 
