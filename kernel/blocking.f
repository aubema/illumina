c=======================================================================
c  Routine blocking
c
c  
c------------------------------------------------------------------------
c   
c    Copyright (C) 2024  Martin Aube
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
c    Contact: aubema@gmail.com
c
c
      subroutine blocking(x_1,y_1,z_1,x_2,y_2,z_2,dx,dy,nbx,nby,
     +altsol,drefle,ofill,obsH,hh,ff1,ff2)
      integer width                                                  ! Matrix dimension in Length/width and height
      parameter (width=512)
      real*8 hh,ff1,ff2,ofill(width,width),drefle(width,width)
      real*8 altsol(width,width),dx,dy,z_1,z_2,rx_1,rx_2,ry_1,ry_2,dho
      real*8 angazi,zhoriz,dh,obsH(width,width),azen,pi
      integer nbx,nby,x_1,x_2,y_1,y_2
      pi=3.141592653589793238462
      ff1=0.
      ff2=0.
      hh=1.
      rx_1=dble(x_1)*dx
      ry_1=dble(y_1)*dy
      rx_2=dble(x_1)*dx
      ry_2=dble(y_1)*dy
      dho=sqrt((rx_1-rx_2)**2.+(ry_1-ry_2)**2.)
      if (dho.gt.0.) then
        call anglezenithal(rx_1,ry_1,z_1,rx_2,ry_2,z_2,azen)
        call angleazimutal(rx_1,ry_1,rx_2,ry_2,angazi)
        if (azen.gt.pi/4.) then                                          ! 45deg. it is unlikely to have a 1km high mountain less than 1
          call horizon(x_1,y_1,z_1,dx,dy,altsol,angazi,zhoriz,dh)
          if ((dho.le.dh).or.((dho.gt.dh).and.
     +    (azen-zhoriz.lt.0.00001))) then
            hh=1.
          else
            hh=0.
          endif
        endif   
        if ((x_1.lt.1).or.(x_1.gt.nbx).or.(y_1.lt.1).or.
     +  (y_1.gt.nby)) then
          ff1=0.
        else
c sub-grid obstacles
          if (dho.gt.drefle(x_1,y_1)) then               ! light path to observer larger than the mean free path -> subgrid obstacles
            angmin=pi/2.-datan2((altsol(x_1,y_1)+
     +      obsH(x_1,y_1)-z_1),drefle(x_1,y_1))
       
            if (azen.ge.angmin) then                                     ! condition sub-grid obstacles direct.
              ff1=ofill(x_1,y_1)
            endif
          endif
        endif                                                            ! end light path to the observer larger than mean free path
        call anglezenithal(rx_2,ry_2,z_2,rx_1,ry_1,z_1,azen)
        if ((x_2.lt.1).or.(x_2.gt.nbx).or.(y_2.lt.1).or.
     +  (y_2.gt.nby)) then
          ff2=0.
        else
          if (dho.gt.drefle(x_2,y_2)) then                               ! light path from source larger than the mean free path -> subgrid obstacles
            angmin=pi/2.-datan2((altsol(x_2,y_2)+
     +      obsH(x_2,y_2)-z_2),drefle(x_2,y_2))
            if (azen.ge.angmin) then                                     ! condition sub-grid obstacles direct.
              ff2=ofill(x_2,y_2)
            endif
          endif  
        endif
      else
        hh=1.
        ff1=0.
        ff2=0.
      endif                                                            ! end light path to the observer larger than mean free path  
      end
