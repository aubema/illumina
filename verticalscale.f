c=======================================================================
c  Routine verticalscale (Martin Aube 2017)
c
c
c  Determine the cell height and thickness for the illumina model
c
c-----------------------------------------------------------------------
c
c    Copyright (C) 2017  Martin Aube
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
      subroutine verticalscale(dx,cthick,cellh)
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=1024)
      integer nz
      real cthick(height)                                                 ! voxel thickness array (meter)
      real cellh(height)                                                  ! voxel height array (meter)
      real expo                                                           ! multiplicative factor to increase the cell thinkness Thick=thick_0*expo**(zcell_c-1)
c      expo=1.111                                                          ! the magic number of Mont-Megantic observatory ;-) Cheers Bernard Malenfant!
      expo=1.
c      cthick(1)=50                                                       ! thickness of the first level
      cthick(1)=dx
      cellh(1)=cthick(1)/2.
      do nz=2,height
         cthick(nz)=cthick(nz-1)*expo
         cellh(nz)=cellh(nz-1)+cthick(nz-1)/2.+cthick(nz)/2.
      enddo
      return
      end
