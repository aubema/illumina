c=======================================================================
c  Routine horizon (Martin Aube 2017)
c
c
c  Determine si la lumiere est bloquee par l'horizon
c
c  pour utilisation avec Illumina 
c-----------------------------------------------------------------------
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
      subroutine horizon(x_c,y_c,z_c,dx,dy,nbx,nby,altsol,
     +latitu,angz,anga,zhoriz) 
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=1024)      
      integer i,x_c,y_c,nbx,nby,cloudz,lcib(width,3),ncib
      integer ii,jj
      real z_c,altsol(width,width),dx,dy,pi,zz
      real latitu,lat,a,b,rterre,angvis,anga,zhoriz,angz,hormin
      real scalef
      pi=3.1415926 
      scalef=1.01  
      cloudz=height
      zhoriz=pi/2.
      lat=latitu*pi/180.
      a=6378137.0                                                         ! earth equatorial radius
      b=6356752.3                                                         ! earth polar radius
      rterre=a**2.*b**2./((a*cos(lat))**2.+(b*sin(lat))**2.)
     + **(1.5)
      angvis=90.-angz
      hormin=pi/2.
      call lignevisee(x_c,y_c,z_c,dx,dy,angvis,anga,nbx,nby,              ! Determination of the viewing line (line of sight voxels).
     +cloudz,lcib,ncib,scalef) 
      do i=1,ncib                           
         ii=lcib(i,1)
         jj=lcib(i,2)
         zz=altsol(ii,jj)
         if (altsol(ii,jj).gt.altsol(x_c,y_c)) then
            zhoriz=pi/2.-atan((altsol(ii,jj)-altsol(x_c,y_c))/sqrt(dx**
     +      2.*real((ii-x_c))**2.+dy**2.*real((jj-y_c))**2.))
         else
            zhoriz=pi/2.-0.5*pi/180.                                      ! bug for zhoriz=pi, anyway in the real world pi is almost impossible                                                                                       
         endif
         if (zhoriz.lt.hormin) hormin=zhoriz
      enddo
      zhoriz=hormin
      return
      end 
