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
      subroutine horizon(d2,x_c,y_c,z_c,dx,dy,nbx,nby,altsol,
     +latitu,angz,anga,zhoriz) 
      integer i,x_c,y_c,nbx,nby,vistep,cloudz,lcib(1024,3),ncib
      integer ii,jj,d2,dsqr
      real z_c,cell_h(50),altsol(1024,1024),dx,dy,pi,zz
      real latitu,lat,a,b,rterre,angel,angaz,anga,zhoriz,angz
      real minzhor,celthi(50)
      data cell_h /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,            ! Hauteur du centre de chaque niveau
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/
      data celthi /0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,            ! Epaisseur des niveaux
     a 2.68,3.24,3.92,4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,
     b 25.94,31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
     c 171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,778.66,
     d 940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,3530.84,
     e 4265.26,5152.44/
      pi=3.1415926   
      vistep=1
      cloudz=50
      zhoriz=pi                                                           ! minimum near horizon angle allowed
      lat=latitu*pi/180.
      a=6378137.0                                                         ! earth equatorial radius
      b=6356752.3                                                         ! earth polar radius
      rterre=a**2.*b**2./((a*cos(lat))**2.+(b*sin(lat))**2.)
     + **(1.5)
      angel=(pi/2.-angz)*180./pi                                          ! lignevisee require angles in degrees and elevation angle
      angaz=anga*180./pi
      call lignevisee(x_c,y_c,z_c,dx,dy,angel,angaz,nbx,nby,              ! Determination of the viewing line (line of sight voxels).
     +vistep,cloudz,lcib,ncib)
      minzhor=pi
      do i=1,ncib-1 
         ii=lcib(i,1)
         jj=lcib(i,2)
         zz=cell_h(lcib(i,3))
         dsqr=(x_c-ii)**2+(y_c-jj)**2
         if (dsqr.le.d2) then
            if (zz.le.altsol(ii,jj)) then
c               print*,zz,altsol(ii,jj)
               call anglezenithal(x_c,y_c,z_c,ii,jj,altsol(ii,jj),dx,dy
     +         ,zhoriz)                                                   ! calculate and return the zenithal angle in radian: zhoriz
c               print*,pi/2.-zhoriz
               if (zhoriz.lt.minzhor) minzhor=zhoriz
            endif
         endif
      enddo
      zhoriz=minzhor
c              print*,zhoriz,pi,pi/2.,dsqr
      if (abs(zhoriz-pi/2.).le.(3.*pi/180.)) then                        ! semble y avoir un probleme avec un angle pres de pi/2
         if (zhoriz.le.pi/2.) zhoriz=zhoriz-3.*pi/180.
      endif
      if (zhoriz.ge.(pi/2.-1.5*pi/180.)) zhoriz=pi
c         print*,zhoriz
      return
      end 
