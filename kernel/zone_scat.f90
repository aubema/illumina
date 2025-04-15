!=======================================================================
!  Routine zone scattering 2024
!
!  Determine les cellules se trouvant dans la zone de diffusion (effet)
!  des cellules (x1,y1,z1) et (x2,y2,z2) 
!  Retourne la matrice des cellules diffusantes (diffusion) ainsi que le nombre
!  de cellules diffusantes (ncellule)
!
!  
!------------------------------------------------------------------------
!   
!    Copyright (C) 2024  Martin Aube
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Publi! License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Publi! License for more details.
!
!    You should have received a copy of the GNU General Publi! License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    Contact: martin.aube@cegepsherbrooke.qc.ca
!
!
      subroutine zone_scat(xi,yi,zi,xf,yf,zf,effet,zondif,ncell,siz,siz_0)
       implicit none
       integer i,j,k
       integer ncell,neffet,imin,imax,jmin,jmax,kmin,kmax,stepdi
       real*8 x0,y0,z0,xi,yi,zi,xf,yf,zf
       real*8 effet,dmax,d
       real*8 xmin,xmax,ymin,ymax,zmin,zmax
       real*8 zondif(3000000,3),siz,siz_0
       real*8 pi
       pi=3.141592654D0
       dmax=dsqrt((xf-xi)**2.+(yf-yi)**2.+(zf-zi)**2.)
       siz=siz_0*(((4.*effet**3.)/3.+dmax*effet**2.)/((4.*effet**3.)/3.))**(1./3.)  ! explain that calculation
       ! I think this is to keep the number of cells approximately constant but the possible problem is that 
       ! we do not use the right resolution to perform the extrapolation of the radiance à infinite resolution. 
       ! I assume this error is not so important but that should be investigated. MA
       dmax=dmax+effet
       xmin=3000000.
       ymin=3000000.
       zmin=3000000.
       xmax=0.
       ymax=0.
       zmax=0.
       if (xf.lt.xmin) xmin=xf
       if (xi.lt.xmin) xmin=xi
       if (yf.lt.ymin) ymin=yf
       if (yi.lt.ymin) ymin=yi
       if (zf.lt.zmin) zmin=zf
       if (zi.lt.zmin) zmin=zi   
       if (xf.gt.xmax) xmax=xf
       if (xi.gt.xmax) xmax=xi
       if (yf.gt.ymax) ymax=yf
       if (yi.gt.ymax) ymax=yi
       if (zf.gt.zmax) zmax=zf
       if (zi.gt.zmax) zmax=zi          
       
                
       neffet=idnint(1.2*effet/siz)

! limits of the calculations loops
       imin=idnint(xmin/siz)-neffet
       imax=idnint(xmax/siz)+neffet
       jmin=idnint(ymin/siz)-neffet
       jmax=idnint(ymax/siz)+neffet
       kmin=idnint(zmin/siz)-neffet
       kmax=idnint(zmax/siz)+neffet
       ncell=0
       do i=imin,imax
         x0=dble(i)*siz
         do j=jmin,jmax
           y0=dble(j)*siz
           do k=kmin,kmax                                                    
             z0=dble(k)*siz
             d=dsqrt((x0-xi)**2.+(y0-yi)**2.+(z0-zi)**2.)+dsqrt((x0-xf)**2.+(y0-yf)**2.+(z0-zf)**2.)
             if (d.le.dmax) then                                                ! ensure spherical zone r<dmax
               ncell=ncell+1
               zondif(ncell,1)=x0                                   
               zondif(ncell,2)=y0 
               zondif(ncell,3)=z0
             endif
           enddo
         enddo
       enddo

       
!       open(unit=1,file='scatzon.txt',status='unknown')
!         do i=1,ncell
!            write(1,*) zondif(i,1),zondif(i,2),zondif(i,3)
!         enddo
!       close(1)
!       open(unit=1,file='points.txt',status='unknown')
!          write(1,*) xi,yi,zi
!          write(1,*) xf,yf,zf
!       close(1)
       return
       end 
