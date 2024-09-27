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
      subroutine zone_diffusion(effet,zondif,ncell,siz)
       implicit none
       integer i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax,kmin,kmax
       real*8 x0,y0,z0
       real*8 effet,dmax,d
       real*8 zondif(3000000,3),siz
       real*8 pi
       pi=3.141592654D0
       neffet=idnint(effet/siz)
       dmax=effet
c limits of the calculations loops
       imin=-neffet
       imax=+neffet
       jmin=-neffet
       jmax=+neffet
       kmin=-neffet
       kmax=neffet
       ncell=0
       do i=imin,imax
         x0=dble(i)*siz
         do j=jmin,jmax
           y0=dble(j)*siz
           do k=kmin,kmax                                                     
             z0=dble(k)*siz
             d=dsqrt(x0**2.+y0**2.+z0**2.)
             if (d.le.dmax) then                                                ! ensure spherical zone r<dmax
               ncell=ncell+1
               zondif(ncell,1)=x0                                   
               zondif(ncell,2)=y0 
               zondif(ncell,3)=z0
             endif
           enddo
         enddo
       enddo
       return
       end 
