c programm to integrate values on a circle centered at position x,y
c with a radius r from a bin file
c 
c
       program circular-integ
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
       real data(width,width),total
       real pixsiz,r,x,y
       integer nbx,nby,i,j
       character*72 basename
       print*,'bin file name?'
       read*,basename
       print*,'Radius (pixels)?'
       read*,radius
       print*,'x y coordinates (pixel)?'
       read*,x,y
       call 2din(nbx,nby,basename,data)
       total=0.
       do i=1,nbx
          do j=1,nby
             r=sqrt(real(i-x)**2.+real(j-y)**2.)
             if (r.le.radius) then
                total=total+data(i,j)
             endif
          enddo
       enddo
       print*,'Integrated value over a circular area of radius=',
     + radius,' centered at x=',x,' and y=',y,': ',total
       stop
       end
