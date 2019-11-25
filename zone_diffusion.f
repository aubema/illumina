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
     +         nbx,nby,alt_sol,cloudz,zondif,ncell,stepdi,n2nd)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=1024)
       integer x_1,y_1,x_2,y_2,z_2,nbx,nby,i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax,kmax
       integer zondif(3000,4),keep,stepdi,cloudz,n2nd,dstep,nexp
       real x1,y1,z1,x2,y2,z2,x0,y0,z0,alt_sol(width,width)
       real dx,dy,effet,dmin,aire,a,b,c,s,delta,d,deltmx,d1,d2
       real cell_t(height),cell_h(height)
       call verticalscale(dx,cell_t,cell_h)                                ! define the vertical scale
       neffet=nint(effet/(cell_t(1)))

       
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
       if (kmax.gt.cloudz) then
          kmax=cloudz
       endif
       nexp=((imax-imin+1)*(jmax-jmin+1))*(kmax-1)      
       dstep=1
 10    ncell=0
       keep=0

       do i=imin,imax
        do j=jmin,jmax
         do k=2,kmax                                                      ! forbid an artefact coming from source too close to the voxel (2 is the minimum)
          x0=real(i)*dx
          y0=real(j)*dy
          z0=cell_h(k)
          d1=sqrt((x1-x0)**2.+(y1-y0)**2.+(z1-z0)**2.)
          d2=sqrt((x2-x0)**2.+(y2-y0)**2.+(z2-z0)**2.)
          d=d1+d2
          dmin=sqrt((x1-x2)**2.+(y1-y2)**2.+(z1-z2)**2.)
          if (z0.gt.alt_sol(i,j)) then
            if (d.le.dmin+effet) then
               if (ncell.gt.n2nd) then
                 if (nexp/(stepdi*2).gt.n2nd) then
                   dstep=dstep*2
                 else
                   dstep=dstep-1
                   if (dstep.lt.1) dstep=1
                 endif
                 stepdi=stepdi+dstep
c        print*,' Set step of scattering to:',stepdi,ncell,nexp,dstep,
c     +  ncell,neffet
                 goto 10
               endif
               if (keep.eq.0) then
                  zondif(ncell,1)=i                                   
                  zondif(ncell,2)=j 
                  zondif(ncell,3)=k
                  ncell=ncell+1
               endif
               keep=keep+1
               if (keep.eq.stepdi) keep=0
            endif
          endif                                                           ! fin condition au-dessus du sol
         enddo
        enddo
       enddo
       ncell=ncell-1
       return
       end 
