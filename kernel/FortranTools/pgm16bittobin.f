c programm to convert pgm 16 bit to bin file
c
       program pgm16bittobin
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
       real pixsiz,gain,offset,dat1(width,width)
       real xcell0,ycell0
       integer nbx,nby,i,j,valmax
       character*72 name1,name3
       print*,'Enter pgm file name and output bin file name:'
       read*,name1,name3
       call intrants2d(name1,dat1,xcell0,ycell0,pixsiz,nbx,nby)
       call twodout(nbx,nby,name3,dat1)
       stop
       end
