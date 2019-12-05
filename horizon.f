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
      subroutine horizon(x,y,z,dx,dy,anga,zhoriz,d) 
      integer width                                                       ! Matrix dimension in Length/width and height
      parameter (width=512)      
      integer x,y,nx,ny
      real dx,dy,altsol(width,width),anga,zout,pi,angaz1,ix,iy,dist
      real posx,posy,scalef,zhoriz,z,d
      pi=3.1415926
c      angaz1 = (pi*anga)/180.
       andaz1=anga
      ix = (cos(angaz1))                                                  ! viewing vector components
      iy = (sin(angaz1))
      dist=(real(width))*sqrt(1.+tan(angaz1)**2.)
      scalef=dx/3.     
      posx=real(x)*dx
      posy=real(y)*dy
      zhoriz=pi
      do while (((posx.le.real(width)*dx).and.(posx.gt.1.*dx)).and.
     +((posy.le.real(width)*dy).and.(posy.gt.1.*dy)))
        posx=posx+ix*scalef
        posy=posy+iy*scalef
        nx=nint(posx/dx)
        ny=nint(posy/dy)
c        print*,nx,ny,x,y
        if (altsol(nx,ny).gt.z) then        
          zout=pi/2.-atan((altsol(nx,ny)-z)/sqrt(dx**
     +    2.*real((nx-x))**2.+dy**2.*real((ny-y))**2.))
          d=sqrt(dx**2.*real((nx-x))**2.+dy**2.*real((ny-y))**2.)
        else
          zout=pi/2.-0.5*pi/180.                                          ! bug for zhoriz=pi, anyway in the real world pi is almost impossible 
          d=real(width)*dx
        endif        
        if (zout.lt.zhoriz) then
           zhoriz=zout
        endif
      enddo
      return
      end 
