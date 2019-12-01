c======================================================================
c  Routine lignevisee (Martin Aube 2005)
c
c  Determine les cellules se trouvant dans la ligne de visee (cellules cibles)
c  entre les cellules (x1,y1,z1) et (x2,y2,z2) 
c  Retourne la matrice des cellules cibles (visee) ainsi que le nombre
c  de cellules cibles (ncell)
c
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2009  Martin Aube
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

      subroutine lignevisee (x1,y1,z1,dx,dy,angvis,angazi,
     + nbx,nby,cloudz,visfin,ncfin,scalef)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=256,height=2048)
      integer x1,y1,cx,cy,visee(width,3)
      integer visfin(width,3)
      integer viseef(width,3),ncellf,cxm,cym,czm,cloudz
      integer ncell,nbx,nby,a,cxp,cyp,czp,cz,ncfin,i,ii,j
      real z1,xn,yn,zn,dx,dy
      real celthi(height),cell_height(height),pi
      real angvis,angazi,ix,iy,iz,da,angaz1,angvi1
      real r,dr
      real dseuil,dp,dm,dn,scalef
      integer cxn,cyn,czn
      parameter (pi=3.1415926)
      dseuil=sqrt(0.75)                                                   ! no more than one cell away
c
c determining the vertical scale
c
      call verticalscale(dx,celthi,cell_height)
c
      cz=1
      if ((cloudz.ne.height).and.(verbose.eq.1)) then
         print*,'Cloud base vertical level:',cloudz,'/',height
         print*,'Cloud base height (m):',cell_height(cloudz)
      endif
      ncell=0
      xn=int((real(x1-1)*dx))                                            ! Transfer cell coordinates to meter
      yn=int((real(y1-1)*dy))                                            ! reference system starting at the upper left corner to the observer position
      zn=z1                                                               
      angvi1 = (pi*angvis)/180.
      angaz1 = (pi*angazi)/180.
      ix = ( sin((pi/2.)-angvi1) ) * (cos(angaz1))                        ! viewing vector components
      iy = ( sin((pi/2.)-angvi1) ) * (sin(angaz1))
      iz = (sin(angvi1))
      do k = 1,height                                                     ! find initial level
         if ((z1 .lt. cell_height(k)+celthi(k)/2.) .and. 
     +   (z1 .ge. cell_height(k)-celthi(k)/2.)) then
            cz= k
         endif
      enddo
      if (celthi(cz).gt.dx) then
         dr=(dx+dy)/2./100.
      else
         dr=celthi(cz)/100.
      endif                                                               ! linear increment along the line of sight
      r=0.
      cx=x1
      cy=y1
      cxp=0
      cyp=0
      czp=0
      do while ((cx.lt.nbx).and.(cx.ge.1).and.(cy.lt.nby).and.(cy.ge.1)   ! inside the domain?
     +      .and.(cz.lt.height).and.(cz.ge.1))
            r=r+dr
            dr=dr*scalef                        
            cx = x1 + int(ix*r/dx)
            cy = y1 + int(iy*r/dy)
            z = z1 + iz*r
            do k = 1,height
               if ((z .lt.cell_height(k)+celthi(k)/2.).and. 
     +         (z .ge. cell_height(k)-celthi(k)/2.)) then
                  cz= k
               endif
            enddo
            if ((cx.eq.cxp).and.(cy.eq.cyp).and.(cz.eq.czp)) then
            else
               ncell=ncell+1
               visee(ncell,1)=cx
               visee(ncell,2)=cy
               visee(ncell,3)=cz
               cxp=cx
               cyp=cy
               czp=cz
            endif
      enddo  
      do jj=1,1                                                           ! allow multiple pass

      i=0
      ncellf=1
      viseef(ncellf,1)=visee(1,1)
      viseef(ncellf,2)=visee(1,2)
      viseef(ncellf,3)=visee(1,3)
      do while (i.le.ncell-2) 
         i=i+1
         cx=visee(i,1)
         cy=visee(i,2)
         cz=visee(i,3)
         cxp=visee(i+1,1)
         cyp=visee(i+1,2)
         czp=visee(i+1,3)
         cxm=visee(i+2,1)
         cym=visee(i+2,2)
         czm=visee(i+2,3)
         cxn=visee(i+3,1)
         cyn=visee(i+3,2)
         czn=visee(i+3,3)
         dp=(real(cx-cxp)**2.+real(cy-cyp)**2.+real(cz-czp)**2.)**0.5
         dm=(real(cx-cxm)**2.+real(cy-cym)**2.+real(cz-czm)**2.)**0.5
         dn=(real(cx-cxn)**2.+real(cy-cyn)**2.+real(cz-czn)**2.)**0.5
         if ((dp.lt.2.1*dseuil).and.(dm.lt.2.1*dseuil).and.                       ! if the 3 next voxels are within a cell width distance retain le farthest
     +   (dn.lt.2.1*dseuil)) then
            ncellf=ncellf+1
            viseef(ncellf,1)=cxn
            viseef(ncellf,2)=cyn
            viseef(ncellf,3)=czn
            i=i+2
         else
            if ((dp.lt.2.1*dseuil).and.(dm.lt.2.1*dseuil)) 
     +      then
               ncellf=ncellf+1
               viseef(ncellf,1)=cxm
               viseef(ncellf,2)=cym
               viseef(ncellf,3)=czm
               i=i+1
            else
               ncellf=ncellf+1
               viseef(ncellf,1)=cxp
               viseef(ncellf,2)=cyp
               viseef(ncellf,3)=czp
            endif
         endif
      enddo
      do i=1,ncellf
         visee(i,1)=viseef(i,1)
         visee(i,2)=viseef(i,2)
         visee(i,3)=viseef(i,3)
      enddo
      ncell=ncellf
      enddo
c
c stop line of sight to the cloud base and forbid cells outside the domain
c
      ncfin=1
      do ii=1,ncellf
         if (viseef(ii,3).le.cloudz) then
            if ((viseef(ncfin,1).le.nbx).and.(viseef(ncfin,1).ge.1)       ! still inside the domain?
     +      .and.(viseef(ncfin,2).le.nby).and.(viseef(ncfin,2).ge.1)
     +      .and.(viseef(ncfin,3).le.height).and.(viseef(ncfin,3).
     +      ge.1)) then   
               visfin(ii,1)=viseef(ncfin,1)
               visfin(ii,2)=viseef(ncfin,2)
               visfin(ii,3)=viseef(ncfin,3)
               ncfin=ncfin+1
            endif
         endif
      enddo
      vistep=1
      ncfin=ncfin-1 
      return
      end
