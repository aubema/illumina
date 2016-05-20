c programme pour integrer les signal sur une surface circulaire de position x,y
c en conformite avec les coordonnees de imagemagick (display) avec un rayon r
c
       program circularpgminteg
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
       real data(width,width),total
       real pixsiz,r,x,y
       real xcell0,ycell0
       integer nbx,nby,i,j
       character*72 basename
       character*12 nom
       print*,'pgm file name?'
       read*,basename
       print*,'Radius (pixels)?'
       read*,radius
       print*,'x y coordinates (according to imagemagick definition)?'
       read*,x,y
     
       nom='data '
       call intrants2d(basename,data,xcell0,ycell0,pixsiz,nbx,nby)
            total=0.
            do i=1,nbx
             do j=1,nby
                r=sqrt(real(i-x-1)**2.+real(j-(nby-y))**2.)
                if (r.le.radius) then
                  total=total+data(i,j)
              endif
             enddo
            enddo
         print*,'Integrated value over a circular area of radius=',
     +   radius,' centered at x=',x,' and y=',y,': ',total
       stop
       end
