c=======================================================================
c  Routine angleazimutal (Martin Aube 2010)
c
c
c  Determine l'angle azimutal entre les points (x1,y1,z1) et (x2,y2,z2)
c  Retourne l'angle anglezen en radians
c
c  pour utilisation avec Illumina 
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2010  Martin Aube
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
      subroutine angleazimutal(x1,y1,x2,y2,dx,dy,angleazi) 
      integer x1,y1,x2,y2
      real dx,dy,pi,angleazi                                                   
      parameter (pi=3.1415926)   
      angleazi=abs(atan(real(y2-y1)*dy/(real(x2-x1)*dx)))
      if (((x2-x1).eq.0).and.((y2-y1).eq.0)) then
         angleazi=0.
      else
        if (x2-x1.ge.0.) then
         if (y2-y1.lt.0.) then
           angleazi=2.*pi-angleazi 
         endif 
        else
         if (y2-y1.lt.0.) then
           angleazi=angleazi+pi 
         else
           angleazi=pi-angleazi
         endif       
        endif
        if ((angleazi.lt.0.).or.(angleazi.gt.2.*pi)) then
          print*,'ERREUR angleazi=',angleazi,x1,y2,x2,y2
          stop
        endif
c        if ((x1.eq.x2).and.(y1.eq.y2)) then
c          print*,'ERREUR cant compute angle between identical points!'
c          stop
c        endif
      endif
      return
        end 
