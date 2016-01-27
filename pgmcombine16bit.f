c programme pour additionner des pgm 16 bit
c
       program pgmcombine16bit
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
       real total,max
       real pixsiz,gain,offset,dat1(width,width)
       real xcell0,ycell0,dat3(width,width)
       integer nbx,nby,i,j,valmax,nimg,n
       character*72 name1,name3
       character*12 nom
       do i=1,width
       do j=1,width
          dat3(i,j)=0.
       enddo
       enddo
       max=0.
       print*,'Number of pgm to add?'
       read*,nimg
       print*,'Output pgm file name?'
       read*,name3
       do n=1,nimg
          print*,'pgm file name ',n
          read*,name1
          call intrants2d(name1,dat1,xcell0,ycell0,pixsiz,nbx,nby)
          do i=1,nbx
             do j=1,nby
                dat3(i,j)=dat3(i,j)+dat1(i,j)
                if (dat3(i,j).gt.max) max=dat3(i,j)
             enddo
          enddo
       enddo
       gain=max/65535.
       offset=0.
       valmax=65535
       nom='sum'
       call extrants2d (name3,dat3,nom,xcell0,ycell0,pixsiz,
     + gain,offset,nbx,nby,valmax)
       stop
       end
