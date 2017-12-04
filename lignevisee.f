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
     + nbx,nby,vistep,cloudz,visfin,ncfin)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=100)
      integer x1,y1,cx,cy,visee(width,3),alim,vistep,visfin(width,3)
      integer viseef(width,3),ncellf,cxm,cym,czm,elimf,cloudz
      integer ncell,nbx,nby,a,cxp,cyp,czp,cz,ncfin,domain
      real z1,xn,yn,zn,dx,dy,distance
      real celthi(height),cell_height(height),pi
      real angvis,angazi,ix,iy,iz,amax,da,angaz1,angvi1
      real dminx,dminy,dminz,r,dr,rmax

      parameter (pi=3.1415926)
c
c determining the vertical scale
c
      call verticalscale(celthi,cell_height)
c
      cz=1
      if (cloudz.ne.height) then
         print*,'Cloud base vertical level:',cloudz,'/',height
         print*,'Cloud base height (m):',cell_height(cloudz)
      endif
      rmax=sqrt((real(nbx)*dx)**2.+(real(nby)*dy)**2.+30000.**2.)
      ncell=0
      xn=nint((real(x1-1)*dx))                                            ! Transfert des coordonnees des cellules en mtre
      yn=nint((real(y1-1)*dy))                                            ! Le zero du systeme de reference est au coin superieur gauche du pixel observateur
      zn=z1                                                               !
      angvi1 = (pi*angvis)/180.
      angaz1 = (pi*angazi)/180.
      ix = ( sin((pi/2.)-angvi1) ) * (cos(angaz1))                        ! determiner les projections du vecteur de visee selon chaque axe
      iy = ( sin((pi/2.)-angvi1) ) * (sin(angaz1))
      iz = (sin(angvi1))
c      print*,'iz=',iz,ix,iy,angvis,angvi1
      cxp=0
      cyp=0
      czp=0
      do k = 1,height                                                         ! trouver le niveau initial
         if ((z1 .lt. cell_height(k)+celthi(k)/2.) .and. 
     +   (z1 .ge. cell_height(k)-celthi(k)/2.)) then
            cz= k
         endif
      enddo
      if (celthi(cz).gt.dx) then
         da=(dx+dy)/2./100.
      else
         da=celthi(cz)/100.
      endif                                                               ! determiner l'increment lineaire pour le calcul de la ligne de visee
      if ((real(nbx)*dx.gt.real(nby)*dy).and.(real(nbx)*dx.gt.
     +cell_height(height)+celthi(height)/2.)) then                        ! determiner la dimension maximale a parcourir
         amax=real(nbx)*dx
      elseif ((real(nby)*dy.gt.real(nbx)*dx).and.(real(nby)*dy.gt.
     +cell_height(height)+celthi(height)/2.)) then
         amax=real(nby)*dy
      else
         amax=cell_height(height)+celthi(height)/2.
      endif
      if (amax.lt.cell_height(height)) amax=cell_height(height)
      if (abs(iz).gt.0.017) then 
           alim=2*nint(amax/da/iz)
c           print*,'oblique or vertical'
      else
           alim=2*nint(amax/da)
c           print*,'horizontal'
      endif
      domain=1
      a=0
      r=0.
      dr=da
      cx=x1
      cy=y1
      cz=25                                                               ! somewhere in the middle vertically - in order to begin inside the domain
      do while ((cx.le.nbx).and.(cx.ge.1).and.(cy.le.nby).and.(cy.ge.1)   ! verifier si nous sommes dans le domaine
     +      .and.(cz.le.height).and.(cz.ge.1).and.(r.lt.rmax))
            r=r+dr
            dr=dr*1.0005
            cx = x1 + nint(ix*r/dx)
            cy = y1 + nint(iy*r/dy)
            z = z1 + iz*r
            do k = 1,height
               if ((z .lt.cell_height(k)+celthi(k)/2.).and. 
     +         (z .ge. cell_height(k)-celthi(k)/2.)) then
                  cz= k
               endif
            enddo
               domain=1
               if (z.le.cell_height(height)) then
                  dminx=abs((ix*real(a)*da/dx-real(nint(ix*real(a)*       ! calcul de la distance entre le centre de la cellule et la position du vecteur en unite de largeur de cellule
     +            da/dx))))
                  dminy=abs((iy*real(a)*da/dy-real(nint(iy*real(a)*
     +            da/dy))))
                  dminz=abs((z-cell_height(cz))/celthi(cz))
                  distance=sqrt(dminx**2.+dminy**2.+dminz**2.)
                  if (distance.lt.0.99) then                               ! ne retenir que les positions s'approchant a moins de la demi d'une cellule
                     if ((cx.eq.cxp).and.(cy.eq.cyp).and.(cz.eq.czp))     ! s'assurer de ne pas compter plus d'une fois la meme cellule
     +                then
                     else   
                        ncell=ncell+1
                        visee(ncell,1)=cx
                        visee(ncell,2)=cy
                        visee(ncell,3)=cz  
                        cxp=cx
                        cyp=cy
                        czp=cz 
                     endif
                  endif
               endif
         a=a+1
c      print*,'a=',a,r,cx,cy,cz
      enddo



c      print*,'a=',a



c
c  eviter les angles droits successifs
c
      elimf=0
      ncellf=1
      viseef(ncellf,1)=visee(1,1)
      viseef(ncellf,2)=visee(1,2)
      viseef(ncellf,3)=visee(1,3)
      do i=2,ncell-1
         cx=visee(i,1)
         cy=visee(i,2)
         cz=visee(i,3)
         cxp=visee(i+1,1)
         cyp=visee(i+1,2)
         czp=visee(i+1,3)
         cxm=visee(i-1,1)
         cym=visee(i-1,2)
         czm=visee(i-1,3)
         if ((((cxp.eq.cx).and.(cyp.eq.cy).and.(czp.eq.cz+1).and.
     +   ((cxm.ne.cxp).or.(cym.ne.cyp))).or.((cxm.eq.cx).and.
     +   (cym.eq.cy).and.(czm.eq.cz-1).and.((cxm.ne.cxp).or.
     +   (cym.ne.cyp)))).and.(elimf.ne.1)) then
c    un cas a eliminer
            elimf=1
         else      
            ncellf=ncellf+1
            viseef(ncellf,1)=cx
            viseef(ncellf,2)=cy
            viseef(ncellf,3)=cz  
            elimf=0 
         endif
      enddo
      ncellf=ncellf+1
      viseef(ncellf,1)=visee(ncell,1)
      viseef(ncellf,2)=visee(ncell,2)
      viseef(ncellf,3)=visee(ncell,3)
c
c
c arreter le ligne de visee au nuage and forbid cells outside the domain
c
          ncfin=1
          do ii=1,ncellf
             if (viseef(ii,3).le.cloudz) then
                if ((viseef(ncfin,1).le.nbx).and.(viseef(ncfin,1).ge.1)   ! verifier si nous sommes dans le domaine
     +          .and.(viseef(ncfin,2).le.nby).and.(viseef(ncfin,2).ge.1)
     +          .and.(viseef(ncfin,3).le.height).and.(viseef(ncfin,3).
     +          ge.1)) then   
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
