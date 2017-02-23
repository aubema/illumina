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
      subroutine horizon(x_c,y_c,z_c,dx,dy,nbx,nby,alt_sol,
     +latitu,angz,anga,zhoriz) 
      integer i,x_c,y_c,nbx,nby,vistep,cloudz,lcib(1024,3),ncib
      integer ii,jj
      real z_c,cell_h(50),alt_sol(1024,1024),dx,dy,pi,zz
      real latitu,lat,a,b,rterre,angvis,anga,zhoriz,angz
      data cell_h /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,            ! Hauteur du centre de chaque niveau
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/
      pi=3.1415926   
      vistep=1
      cloudz=50
      zhoriz=pi/2.
      lat=latitu*pi/180.
      a=6378137.0                                                         ! earth equatorial radius
      b=6356752.3                                                         ! earth polar radius
      rterre=a**2.*b**2./((a*cos(lat))**2.+(b*sin(lat))**2.)
     + **(1.5)
      angvis=90.-angz
      call lignevisee(x_c,y_c,z_c,dx,dy,angvis,anga,nbx,nby,              ! Determination of the viewing line (line of sight voxels).
     + vistep,cloudz,lcib,ncib) 
      do i=1,ncib 
         ii=lcib(i,1)
         jj=lcib(i,2)
         zz=cell_h(lcib(i,3))
         if (zz.le.alt_sol(ii,jj)) then
            call anglezenithal(x_c,y_c,z_c,ii,jj,zz,dx,dy,zhoriz)         ! calculate and return the zenithal angle in radian: zhoriz
         endif
      enddo
      return
      end 
