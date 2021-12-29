c programm to convert bin file to pgm 16 bit
c
       program bintopgm16bit
c   
c    Copyright (C) 2018  Martin Aube
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
       parameter (width=512)
       real total,max
       real pixsiz,gain,offset,dat1(width,width)
       real xcell0,ycell0,dat3(width,width)
       integer nbx,nby,i,j,valmax
       character*72 name1,name3
       character*12 nom
       print*,'Enter bin file name, pixel size in meter and output pgm f
     +ile name:'
       read*,name1,pixsiz,name3
       call twodin(nbx,nby,name1,dat1)
       print*,'Domain size:',nbx,'x',nby
       max=0.
       do i=1,nbx
          do j=1,nby
             if (dat1(i,j).gt.max) max=dat1(i,j)
             if (dat1(i,j).lt.0.) dat1(i,j)=0.
          enddo
       enddo
       gain=max/65535.
       print*,'Max val=',max,'Gain=',gain
       offset=0.
       valmax=65535
       nom='converted'
       xcell0=0.
       ycell0=0.
       call extrants2d(name3,dat1,nom,xcell0,ycell0,pixsiz,
     + gain,offset,nbx,nby,valmax)
       stop
       end
