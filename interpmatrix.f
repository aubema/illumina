c    programme pour interpoler une carte d epaisseur optique
c    sous echantillonnee a partir de quelques mesures ponctuelles 
c    effectuees par des photometres solaires ou par inversion d images 
c    de teledetection.
c
c    L entree du programme est une carte pgm d epaisseur optique
c    qui peut etre incomplete.
c      
c    L interpolation est effectuee a l aide d une fonction poids 
c    ou statistiques soit:
c
c    nearest neighbour ... 0
c    linear .............. 1
c    2nd order ........... 2
c    mean ................ 3
c    minimum ............. 4
c    maximum ............. 5
c    
c    les axes sont definis comme suit
c    x = direction nord = degres
c    y = direction est = degres
c
c

c    Copyright (C) 2011  Martin Aube
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
c -----------------
c   identification des variables 
c
c
c --------------------
c
c   programme principal
c
      subroutine interpmatrix(dat,x1,x2,y1,y2,intype,window,
     +autom,defval)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=100)
c
c ----------
c
c   declaration des variables
c
      real pi2,pi,dmin,maxim,minim,npix
      real dat(width,width),datout(width,width)
      real distan,dtotal,window,ndat
      real defval
      integer i,j,intype,autom
      integer dens(10,width,width)
      integer x1,x2,y1,y2,nbx,nby
       

c   
c ----------
c
c   initialisation des variables
c
       pi=3.141592654
       pi2=pi*pi
       do i=1,width
        do j=1,width
          datout(i,j)=0.
        enddo
       enddo
c
c   calcul de la taille maximale d interpolation
c    
      nbx=x2-x1
      nby=y2-y1
      if (nbx.gt.nby) then
         maxside=2*nbx+1
      else
         maxside=2*nby+1
      endif   

c
c    Calcul de la matrice de voisinage 
c
c    definition des orientations de theta
c               
c    4 3 2   
c    5 X 1
c    6 7 8
c 
c    le dixieme element est le nombre de voisin pour l instant, seul
c    cette valeur est utilisee
c             
              do i=x1,x2
              do j=y1,y2
                do n=1,10
                   dens(n,i,j)=0
                enddo 
               if (dat(i,j).ne.defval) dens(9,i,j)=1
            if (j.ne.nby) then 
               if (dat(i,j+1).ne.defval) dens(1,i,j)=1
            elseif ((i.ne.nbx).and.(j.ne.nby)) then
               if (dat(i+1,j+1).ne.defval) dens(2,i,j)=1
            elseif (i.ne.nby) then
               if (dat(i+1,j).ne.defval) dens(3,i,j)=1 
            elseif ((i.ne.nbx).and.(j.ne.1)) then               
               if (dat(i+1,j-1).ne.defval) dens(4,i,j)=1
            elseif (j.ne.1) then            
               if (dat(i,j-1).ne.defval) dens(5,i,j)=1 
            elseif ((i.ne.1).and.(j.ne.1)) then 
               if (dat(i-1,j-1).ne.defval) dens(6,i,j)=1 
            elseif (i.ne.1) then               
               if (dat(i-1,j).ne.defval) dens(7,i,j)=1
            elseif ((i.ne.1).and.(j.ne.nby)) then                
                if (dat(i-1,j+1).ne.defval) dens(8,i,j)=1  
            endif
                dens(10,i,j)=dens(1,i,j)+ dens(2,i,j)+dens(3,i,j)+ 
     +          dens(4,i,j)+dens(5,i,j)+dens(6,i,j)+dens(7,i,j)+
     +          dens(8,i,j)+dens(9,i,j)                                        
              enddo            
              enddo   
c
c    Taille de la fenetre d interpolation
c                         
c       print*,'Interpolation box size:',window      
c
c ---------
c
c   Interpoler 
c
 225    do 113 nx=x1,x2
           do 114 ny=y1,y2
              dtotal=0.
              dmin=2000000.
              npix=0.
              minim=2000000.
              maxim=-2000000.
              ndat=0.
c
c    corriger pour la repartition inegale des donnees
c    
           
             do i=nx-nint(window),nx+nint(window)
              do j=ny-nint(window),ny+nint(window)
                if (dat(i,j).ne.defval) then
                 distan=sqrt(real(ny-j)**2.+real(nx-i)**2.)
                  
              if (distan.le.0.71*window) then
                 if (distan.lt.0.01) distan=0.01
                 
                          ndat=ndat+1.
                          if (intype.eq.0) then
                             if (distan.lt.dmin) then
                                dmin=distan
                                datout(nx,ny)=dat(i,j)
                             endif                                        
                          elseif (intype.eq.1) then
                             datout(nx,ny)=dat(i,j)/(distan*
     + (real(dens(10,i,j))))+datout(nx,ny) 
                  dtotal=1./(distan*(real(dens(10,i,j))))
     + +dtotal  
                          elseif (intype.eq.2) then
                             datout(nx,ny)=dat(i,j)/(distan**2.*
     + (real(dens(10,i,j))))+datout(nx,ny)
                             dtotal=1./(distan**2.*(real(dens(10,i,j))
     + ))+dtotal
                          elseif (intype.eq.3) then
                             npix=npix+1.
                             datout(nx,ny)=dat(i,j)+datout(nx,ny)
                          elseif (intype.eq.4) then
                             if (dat(i,j).lt.minim) then
                                minim=dat(i,j)
                             endif
                          elseif (intype.eq.5) then
                             if (dat(i,j).gt.maxim) then
                                maxim=dat(i,j)
                             endif
                          endif

                    endif
               endif
              enddo
             enddo

              if (dtotal.eq.0.) dtotal=0.000001   
              if ((intype.eq.1).or.(intype.eq.2)) then
                 datout(nx,ny)=datout(nx,ny)/dtotal   
                 if (ndat.lt.3.) then
  
                    datout(nx,ny)=defval
                    if (autom.eq.1) then
                       window=window*2.
                       goto 225
                    endif
                 endif 
              elseif (intype.eq.3) then
                 if (npix.eq.0.) then
                    datout(nx,ny)=defval
                    if (autom.eq.1) then
                       window=window+8.
                       goto 225
                    endif
                 else
                    datout(nx,ny)=datout(nx,ny)/npix
                 endif
              elseif (intype.eq.4) then
                 datout(nx,ny)=minim
              elseif (intype.eq.5) then
                 datout(nx,ny)=maxim                 
              elseif (intype.eq.0) then
                 if (ndat.lt.1.) then
                    datout(nx,ny)=defval
                 endif                 
              endif
                 if (autom.eq.1) then
                    window=window*2.
                    goto 225
                 endif
 114       continue
 113    continue
       do nx=x1,x2
        do ny=y1,y2
         dat(nx,ny)=datout(nx,ny)
        enddo
       enddo
       return
       end

