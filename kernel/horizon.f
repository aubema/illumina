c=======================================================================
c  Routine horizon (Martin Aube 2017)
c
c
c  Determine si la lumiere est bloquee par l'horizon
c
c  pour utilisation avec Illumina
c-----------------------------------------------------------------------
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
      subroutine horizon(x,y,z,dx,dy,altsol,anga,zhoriz,d)
      integer width                                                       ! Matrix dimension in Length/width and height
      parameter (width=512)
      integer x,y,nx,ny
      real*8 dx,dy,altsol(width,width),anga,zout,pi,angaz1,ix,iy
      real*8 hcur,distc                                                           ! Earth curvature terrain
      real*8 posx,posy,scalef,zhoriz,z,d,dout
      pi=3.141592654
      angaz1=anga
      ix = (dcos(angaz1))                                                  ! viewing vector components
      iy = (dsin(angaz1))
      scalef=dx/3.
      posx=dble(x)*dx
      posy=dble(y)*dy
      zhoriz=pi
      d=(dble(width))*dsqrt(1.+dtan(angaz1)**2.)
      do while (((posx.le.dble(width)*dx).and.(posx.ge.1.*dx)).and.
     +((posy.le.dble(width)*dy).and.(posy.ge.1.*dy)))
        posx=posx+ix*scalef
        posy=posy+iy*scalef
        nx=idnint(posx/dx)
        ny=idnint(posy/dy)
c earth curvature (first order correction)
        distc=dsqrt((dx*dble(nx-x))**2.+(dy*dble(ny-y))**2.)
        call curvature(distc,hcur)
        if ((nx.eq.x).and.(ny.eq.y)) then                                  ! to forbid division by zero
           if (z.gt.altsol(nx,ny)-hcur) then                               ! reverse curvature to limit horizontal distance. curv is negative. This is a hack
              zout=pi
              d=0.
           else
              zout=0.
              d=0.
           endif
        else
        dout=distc
        zout=pi/2.-datan((altsol(nx,ny)-hcur-z)/dout)
        if (altsol(nx,ny)-hcur.eq.z) then
           zout=pi/2.-0.0001*pi/180.                                      ! bug for zhoriz=pi, anyway in the real world pi is almost impossible
        endif
        if (zout.lt.zhoriz) then
           zhoriz=zout
           d=dout
        endif
        endif
      enddo
      end
