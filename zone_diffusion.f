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
      subroutine zone_diffusion(x1,y1,z1,x2,y2,z2,
     +effet,alts,cloudbase,zondif,ncell,stepdi,dstep,n2nd,siz)
       integer height                                                      ! Matrix dimension in Length/width and height
       parameter (height=1024)
       integer x_1,y_1,z_1,x_2,y_2,z_2,nbx,nby,i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax,kmax
       integer keep,stepdi,cloudz,n2nd,dstep
       real x1,y1,z1,x2,y2,z2,x0,y0,z0,alts
       real dx,dy,effet,dmin,aire,a,b,c,s,delta,d,deltmx,d1,d2
       real zondif(3000,3),siz
       neffet=nint(effet/siz)

       
c limits of the calculations loops
       x_1=nint(x1/siz)
       y_1=nint(y1/siz)
       x_2=nint(x2/siz)
       y_2=nint(y2/siz)
       z_2=nint(z2/siz)+1
       if (x_1.le.x_2) then
         imin=x_1-neffet
         imax=x_2+neffet
       else
         imin=x_2-neffet
         imax=x_1+neffet
       endif
       if (imin.lt.1) imin=1 
       if (y_1.le.y_2) then
         jmin=y_1-neffet
         jmax=y_2+neffet
       else
         jmin=y_2-neffet
         jmax=y_1+neffet
       endif
       if (jmin.lt.1) jmin=1
       kmax=z_2+neffet
       if (z2.gt.cloudbase) then
          kmax=nint(cloudbase/siz)
       endif
       if (stepdi.eq.1) stepdi=(kmax-2)*(imax-imin)*(jmax-jmin)/n2nd
       if (dstep.ne.1) dstep=dstep-1
 10    ncell=0
       keep=0
c       print*,imin,imax,jmin,jmax,kmax,x1,x2,y1,y2,x_1,x_2,y_1,y_2
c       print*,neffet,effet
c       stop
       do i=imin,imax
        do j=jmin,jmax
         do k=2,kmax                                                      ! forbid an artefact coming from source too close to the voxel (2 is the minimum)
          x0=real(i)*siz
          y0=real(j)*siz
          z0=real(k)*siz
          d1=sqrt((x1-x0)**2.+(y1-y0)**2.+(z1-z0)**2.)
          d2=sqrt((x2-x0)**2.+(y2-y0)**2.+(z2-z0)**2.)
          d=d1+d2
          dmin=sqrt((x1-x2)**2.+(y1-y2)**2.+(z1-z2)**2.)
          if (z0.gt.alts) then
            if (d.le.dmin+effet) then
               if (ncell.gt.n2nd) then
                 stepdi=stepdi+dstep
                 dstep=dstep*2
c        print*,' Set step of scattering to:',stepdi,dstep
                 goto 10
                endif
               if (keep.eq.0) then
                  zondif(ncell,1)=x0                                   
                  zondif(ncell,2)=y0 
                  zondif(ncell,3)=z0
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
