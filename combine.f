c programm to add 2 bin illumina files
c
       program combine
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
       parameter (width=256)
       real total
       real dat1(width,width)
       real xcell0,ycell0,dat3(width,width)
       integer nbx,nby,i,j,nimg,n
       character*72 name1,name3
       do i=1,width
       do j=1,width
          dat3(i,j)=0.
       enddo
       enddo
       max=0.
       print*,'Number of bin to add?'
       read*,nimg
       print*,'Output bin file name?'
       read*,name3
       do n=1,nimg
          print*,'bin file name ',n
          read*,name1
          call twodin(nbx,nby,name1,dat1)
          do i=1,nbx
             do j=1,nby
                dat3(i,j)=dat3(i,j)+dat1(i,j)
             enddo
          enddo
       enddo
       call twodout(nbx,nby,name3,dat3)       
       stop
       end
