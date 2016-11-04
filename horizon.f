c=======================================================================
c  Routine horizon (Martin Aube 2010)
c
c
c  Determine l'horizon avec une resolution de 1 deg.
c  Retourne une matrice de 360 valeurs
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
      subroutine horizon(x_c,y_c,z_c,d2,alt_sol,nbx,nby,dx,
     +dy,zen_horiz,latitu,anga) 
      integer nbx,nby,ii,jj,x_c,y_c,az,daz,azin,n,nn
      integer bx1,bx2,by1,by2
      real dx,dy,pi,angazi,angzen,alt_sol(1024,1024)
      real zen_horiz(360),z_c,d2,d2p,altitu
      real latitu,corcur,a,b,rterre,anga,lat
      pi=3.1415926   
      lat=latitu*pi/180.
      corcur=0.
      a=6378137.0                                                         ! earth equatorial radius
      b=6356752.3                                                         ! earth polar radius
      rterre=a**2.*b**2./((a*cos(lat))**2.+(b*sin(lat))**2.)
     + **(1.5)
c      print*,'===',x_c,y_c,z_c,d2,'==='
      do jj=1,360
       zen_horiz(jj)=5.
      enddo 
c      do ii=1,nbx
c       do jj=1,nby
        if ((anga.ge.0.).and.(anga.lt.pi/2.)) then
           bx1=x_c
           bx2=x_c+nint(d2/dy)
           by1=y_c
           by2=y_c+nint(d2/dy)
        elseif ((anga.ge.pi/2.).and.(anga.lt.pi)) then
           bx1=x_c-nint(d2/dy)
           bx2=x_c
           by1=y_c
           by2=y_c+nint(d2/dy)
        elseif ((anga.ge.pi).and.(anga.lt.3.*pi/2.)) then
           bx1=x_c-nint(d2/dy)
           bx2=x_c
           by1=y_c-nint(d2/dy)
           by2=y_c
        else
           bx1=x_c
           bx2=x_c+nint(d2/dy)
           by1=y_c-nint(d2/dy)
           by2=y_c
        endif

      do ii=bx1,bx2
       do jj=by1,by2
        d2p=sqrt((real(ii-x_c)*dx)**2.+(real(jj-y_c)*dy)**2.)

c this is a firts attempt to make a simplistic correction for the local 
c curvature of earth the basic idea is to consider that the curvature 
c limit the horizon angle as a mountain located between the observer 
c and observed target. This is not exact because in that approch
c the line of sight toward a remote point at the same elevation of 
c the observer is considered to be at z=90 deg. En fact it should be 
c larger than that because the point should fall under the horizon. 
c Anyway the correction will result in blocking the light ray if the 
c curvature should block the ray.
c 
c To fall back to non correction, please simply comment 2 lines below
c        corcur=sqrt(rterre**2.-(d2/2.-d2p)**2.)-sqrt(rterre**2.-
c     +  (d2/2.)**2)


        if (d2p.le.d2*1.001) then
         if ((ii.eq.x_c).and.(jj.eq.y_c)) then
c calcul de l'angle azimutal de la cell cible à la ligne d'horizon
         else
          call angleazimutal(x_c,y_c,ii,jj,dx,dy,angazi)
          az=nint(angazi*180./pi)+1
          daz=nint(atan(1./sqrt(real(x_c-ii)**2.+real(y_c-jj)**2.))*
     +    180./pi/2.)+2
          altitu=alt_sol(ii,jj)+corcur
          call anglezenithal
     +    (x_c,y_c,z_c,ii,jj,altitu,dx,dy,
     +    angzen) 
c            if (angzen.lt.zen_horiz(az)) then
          azin=az-daz
          if (azin.lt.1) azin=azin+360
          n=azin
          do nn=1,2*daz+1
           if (angzen.lt.zen_horiz(n)) then              
            zen_horiz(n)=angzen                                           ! horizon est l'angle azmutal maximal en rad
           endif
           n=n+1
           if (n.gt.360) n=1
          enddo
c            endif
         endif
        endif
       enddo
      enddo
c      print*,zen_horiz
      return
      end 
