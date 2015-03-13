c programme pour faire une moyenne ponderee de deux pgm 16 bit
c
       program moypondpgm16bit
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
       real pixsiz,gain,offset,dat1(width,width),weight1,weight2
       real xcell0,ycell0,dat2(width,width),dat3(width,width)
       integer nbx,nby,i,j,valmax
       character*72 name1,name2,name3
       character*12 nom
       print*,'1 st pgm file name and averaging weight?'
       read*,name1,weight1
       print*,'2 nd pgm file name and averaging weight?'
       read*,name2,weight2       
       print*,'Output pgm file name?'
       read*,name3
     
       nom='data1 '
       call intrants2d(name1,dat1,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)
       nom='data2 '
        call intrants2d(name2,dat2,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)

           max=0.
           total=weight1+weight2
            do i=1,nbx
             do j=1,nby
                dat3(i,j)=(weight1*dat1(i,j)+weight2*dat2(i,j))/total
                if (dat3(i,j).gt.max) max=dat3(i,j)
             enddo
            enddo
            gain=max/65535.
            offset=0.
            valmax=65535
        nom='average'
        call extrants2d (name3,dat3,nom,xcell0,ycell0,pixsiz,
     +  gain,offset,nbx,nby,valmax)
       stop
       end
