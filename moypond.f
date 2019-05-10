c programm to average two bin illumina files
c
       program moypond
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
       parameter (width=1024)
       real total
       real dat1(width,width),weight1,weight2
       real dat2(width,width),dat3(width,width)
       integer nbx,nby,i,j
       character*72 name1,name2,name3
       print*,'1 st bin file name and averaging weight?'
       read*,name1,weight1
       print*,'2 nd bin file name and averaging weight?'
       read*,name2,weight2       
       print*,'Output bin file name?'
       read*,name3
       call twodin(nbx,nby,name1,dat1)
       call twodin(nbx,nby,name2,dat2)
       total=weight1+weight2
       do i=1,nbx
          do j=1,nby
             dat3(i,j)=(weight1*dat1(i,j)+weight2*dat2(i,j))/total
          enddo
       enddo
       call twodout(nbx,nby,name3,dat3)
       stop
       end
