c algorithme pour varier la resolution en fct de la distance
c
       program varres
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
       integer width
       parameter (width=1024)
       real maxi,rhomax
       real pixsiz,gain,offset,lumlp(width,width),refle(width,width)
       real xcell0,ycell0,lumlpo(width,width),reflo(width,width)
       real lum,rho
       integer imin,jmin,imax,jmax
       integer nx9ma,nx9mi,ny9ma,ny9mi,nx3ma,nx3mi,ny3ma,ny3mi
       integer l3,l1,nx,ny,i,j,valmax,n,px,py
       integer k,l,lx,ly,xi,xf,yi,yf,point(width*width,3),m,p
       character*72 lin,rin,lout,rout
       character*12 nom
       xcell0=0.
       ycell0=0.
       pixsiz=1.
       nom='varres'
       l1=27
c width of the 1x1 resolution window
c l1 must be an odd multiple of 9 (e.g. 9, 27, 45, 63, 81, ...) 
       l3=81     
c width of the  3x3 resolution window
c l3 must be an odd multiple of 9 and l3>l1
       open(unit=11,file='varres.in',status='unknown')
          read(11,*) lin
c initial lumlp file
          read(11,*) rin
c initial reflectance file
          read(11,*) lout
c final lumlp file
          read(11,*) rout
c final reflectance file
          read(11,*) px
c observer x position
          read(11,*) py
c observer y position
          read(11,*) l1, l3
c limit distance of full resolution and 3x3 resolution
       close(unit=11)
       l1=nint(real(l1)/2.-1.)
       l3=nint(real(l3)/2.-1.)
       lx=0
       ly=0
       do i=1,width
       do j=1,width
          lumlpo(i,j)=0.
          reflo(i,j)=0.
       enddo
       point(i,1)=0
       point(i,2)=0
       point(i,3)=0
       enddo
       max=0.

c load lumlp file
          call intrants2d(lin,lumlp,xcell0,ycell0,pixsiz,
     +    nx,ny)
          maxi=0.
          n=0
c create new lumlp file using a variable mesh grid 
          nx9ma=(nx-px)/9
          nx9mi=(px-1)/9
          ny9ma=(ny-py)/9
          ny9mi=(py-1)/9
          imin=-nx9mi+1
          if (imin.lt.1) imin=1
          jmin=-ny9mi+1
          if (jmin.lt.1) jmin=1
          do i=imin,nx9ma-1
            do j=jmin,ny9ma-1
              lx=i*9
              ly=j*9
              k=px+lx
              l=py+ly
              if ((abs(lx).gt.l3).or.(abs(ly).gt.l3)) then
                do m=-4,4
                  do p=-4,4
                     lumlpo(k,l)=lumlpo(k,l)+lumlp(k+m,l+p)
                  enddo
                enddo
                n=n+1
                point(n,1)=k
                point(n,2)=l
                point(n,3)=9
                if (lumlpo(k,l).gt.maxi) maxi=lumlpo(k,l)
              endif
            enddo
          enddo
          nx3ma=(nx-px)/3
          nx3mi=(px-1)/3
          ny3ma=(ny-py)/3
          ny3mi=(py-1)/3
          imin=-nx3mi
          if (imin.lt.2) imin=2
          jmin=-ny3mi
          if (jmin.lt.2) jmin=2
          do i=imin,nx3ma
            do j=jmin,ny3ma
              lx=i*3
              ly=j*3
              k=px+lx
              l=py+ly
              if ((abs(lx).le.l3).and.(abs(ly)
     +        .le.l3)) then
                do m=-1,1
                  do p=-1,1 
                    lumlpo(k,l)=lumlpo(k,l)+lumlp(k+m,l+p)
                  enddo
                enddo
                n=n+1
                point(n,1)=k
                point(n,2)=l
                point(n,3)=3
                if (lumlpo(k,l).gt.maxi) maxi=lumlpo(k,l)
              endif
            enddo
          enddo
          imin=px-l1
          if (imin.lt.1) imin=1
          jmin=py-l1
          if (jmin.lt.1) jmin=1
          imax=px+l1
          if (imax.gt.nx) imax=nx
          jmax=py+l1
          if (jmax.gt.ny) jmax=ny
          do i=imin,imax
            do j=jmin,jmax
              lumlpo(i,j)=lumlp(i,j)
              n=n+1
              point(n,1)=i
              point(n,2)=j
              point(n,3)=1
            enddo
          enddo
          print*,'Maximum acceleration factor: ',nx*ny/n,'X'
c load reflectance
          call intrants2d(rin,refle,xcell0,ycell0,pixsiz,
     +    nx,ny) 
c creating the reflectance file weighted by the lumlp
          rhomax=0.
          do i=1,n
            xi=point(i,1)-(point(i,3)-1)/2
            xf=point(i,1)+(point(i,3)-1)/2
            yi=point(i,2)-(point(i,3)-1)/2
            yf=point(i,2)+(point(i,3)-1)/2
            lum=0.
            rho=0.
            if (xi.lt.1) xi=1
            if (yi.lt.1) yi=1
            do k=xi,xf
              do l=yi,yf
                 rho=rho+refle(k,l)*lumlp(k,l)
                 lum=lum+lumlp(k,l)
              enddo
            enddo
            if (lum.ne.0.) then
               rho=rho/lum
            else
               rho=refle(k,l)
            endif
            if (rho.gt.rhomax) rhomax=rho
            do k=xi,xf
              do l=yi,yf
                 reflo(k,l)=rho
              enddo
            enddo
          enddo
c writing average reflectance pgm
       gain=rhomax/65535.
       offset=0.
       valmax=65535
       nom='avreflectanc'
       call extrants2d(rout,reflo,nom,xcell0,ycell0,pixsiz,
     + gain,offset,nx,ny,valmax)
c writing lumlp output file with variable resolution
       gain=maxi/250000.
       offset=0.
       valmax=65535
       nom='lumlp-newres'
       call extrants2d(lout,lumlpo,nom,xcell0,ycell0,pixsiz,
     + gain,offset,nx,ny,valmax)
c writing grid points informations
       open(unit=1,file='grid.txt',status='unknown')
         write(1,*) n
         do i=1,n
           write(1,*) i,point(i,1),point(i,2),point(i,3)
         enddo
       close(unit=1)
       stop
       end
