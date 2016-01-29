c recreate the full coverage contribution and sensitivity maps
c with the grid.txt file and the lumlp file
c
       program retres
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
       real maxi
       real pixsiz,gain,offset,lumlp(width,width),donnee(width,width)
       real xcell0,ycell0,output(width,width)
       real lum
       integer nx,ny,i,j,valmax,n
       integer k,l,xi,xf,yi,yf,point(width*width,3)
       character*72 lin,fin,fout
       character*12 nom
       nom='sensi_contri'
       do i=1,width
       do j=1,width
          donnee(i,j)=0.
          lumlp(i,j)=0.
          output(i,j)=0.
       enddo
       point(i,1)=0
       point(i,2)=0
       point(i,3)=0
       enddo
       max=0.
       open(unit=1,file='retres.in',status='unknown')
          read(1,*) lin
c initial lumlp file
          read(1,*) fin
c variable resolution file to convert to full resolution
          read(1,*) fout
c full resolution file
       close(unit=1)
c load lumlp file
          call intrants2d(lin,lumlp,xcell0,ycell0,pixsiz,
     +    nx,ny)
          maxi=0.
          n=0
c load variable resolution file
          call intrants2d(fin,donnee,xcell0,ycell0,pixsiz,
     +    nx,ny) 
c loading grid points informations
       open(unit=1,file='grid.txt',status='unknown')
         read(1,*) n
         do i=1,n
           read(1,*) toto,point(i,1),point(i,2),point(i,3)
         enddo
       close(unit=1)
c creating the output file weighted by the lumlp
          maxi=0.
          do i=1,n
            xi=point(i,1)-(point(i,3)-1)/2
            xf=point(i,1)+(point(i,3)-1)/2
            yi=point(i,2)-(point(i,3)-1)/2
            yf=point(i,2)+(point(i,3)-1)/2
            lum=0.
            rho=0.
            do k=xi,xf
              do l=yi,yf
                 lum=lum+lumlp(k,l)
              enddo
            enddo
            do k=xi,xf
              do l=yi,yf
                if (lum.gt.0.) then
                  output(k,l)=lumlp(k,l)/lum*donnee(point(i,1),
     +            point(i,2))
                  if (output(k,l).gt.maxi) maxi=output(k,l)
                else
                  output(k,l)=0.
                endif
              enddo
            enddo
          enddo
c writing output pgm
       gain=maxi/65535.
       offset=0.
       valmax=65535
       nom='finalsenscon'
       call extrants2d(fout,output,nom,xcell0,ycell0,pixsiz,
     + gain,offset,nx,ny,valmax)
       stop
       end
