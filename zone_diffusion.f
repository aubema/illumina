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

       integer x_1,y_1,z_1,x_2,y_2,z_2,i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax,kmax
       integer keep,stepdi,cloudz,n2nd,step
       real x1,y1,z1,x2,y2,z2,x0,y0,z0,alts
       real effet,dmin,d,d1,d2
       real zondif(1000000,3),siz
       if (siz.lt.40.) siz=40.
       neffet=nint(effet/siz)
       dmin=sqrt((x1-x2)**2.+(y1-y2)**2.+(z1-z2)**2.)
c find an approximate value to stepdi
       stepdi=nint((dmin+2.*effet)*3.14159/siz)*neffet*neffet/n2nd
       if (stepdi.eq.0) stepdi=1

       step=nint(real(stepdi)**(1./3.))
       if (step.le.0) step=1
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
       if (y_1.le.y_2) then
         jmin=y_1-neffet
         jmax=y_2+neffet
       else
         jmin=y_2-neffet
         jmax=y_1+neffet
       endif
       kmax=z_2+neffet
       if (z2.gt.cloudbase) then
          kmax=nint(cloudbase/siz)
       endif
 10    ncell=0
       keep=0
       do i=imin,imax,step
        do j=jmin,jmax,step
         do k=1,kmax,step                                                      ! forbid an artefact coming from source too close to the voxel (2 is the minimum)
          x0=real(i)*siz
          y0=real(j)*siz
          z0=real(k)*siz
          d1=sqrt((x1-x0)**2.+(y1-y0)**2.+(z1-z0)**2.)
          d2=sqrt((x2-x0)**2.+(y2-y0)**2.+(z2-z0)**2.)
          d=d1+d2
          if ((z0.gt.alts).and.(z0.lt.35000.)) then
            if (d.le.dmin+2.*effet) then
                  ncell=ncell+1
                  zondif(ncell,1)=x0                                   
                  zondif(ncell,2)=y0 
                  zondif(ncell,3)=z0
            endif
          endif                                                           ! fin condition au-dessus du sol
         enddo
        enddo
       enddo
       stepdi=step**3
       return
       end 
