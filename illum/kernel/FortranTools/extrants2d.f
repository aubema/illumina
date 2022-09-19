c routine d ecriture des extrants 2d
c
c 
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
         subroutine extrants2d(outfile,valeur,nom,xcell0,ycell0,pixsiz,
     +   gain,offset,nbx,nby,valmax)
         integer width
         parameter (width=512)
         real valeur(width,width),xcell0,ycell0,gain,offset
         integer i,j,nbx,nby,valmax
         character*12 nom
         character*72 outfile
         print*,'Ecriture du fichier de ',nom,': ',outfile
         open(unit=1,file=outfile,status='unknown')
         write(1,1000)                                                    ! Tag P2).
         write(1,1001) xcell0
         write(1,1002) ycell0
         write(1,1003) pixsiz
         write(1,1004) gain
         write(1,1005) offset
         write(1,*) nbx,nby,valmax
         if (gain.eq.0.) gain=1.                                          ! case when all the image is 0 (eg. empty zone)
         do i=1,nbx                                                       ! Debut de la boucle sur toutes les cases en x.
           do j=1,nby                                                     ! Debut de la boucle sur toutes les cases en y.
              valeur(i,j)=(valeur(i,j)-offset)/gain                       ! Transformation des donnees avec le gain et l'offset et recherche 
           enddo                                                          ! Fin de la boucle sur toutes les cases en y.
         enddo          
         write(1,*) ((nint(valeur(i,j)),i=1,nbx),j=nby,1,-1)              ! Ecriture de toutes les donnees qui sont ensuite inscrites
c                                                                         ! dans la matricel. Ce sont des boucles imbriquees dans la
c                                                                         ! fonction "write" qui couvrent tout le domaine delimite par nbx et nby.
c                                                                         ! L'increment de la boucle sur les ranges (latitute) est de -1 car 
c                                                                         ! on considere que dans les fichiers,
c                                                                         ! le nombre en haut a gauche des donnees est a la
c                                                                         ! coordonnee (1,nby).
         close(1)                                                         ! Fermeture du fichier relief.pgm.

 1000 format('P2')
 1001 format('# lon0 ',F12.5)
 1002 format('# lat0 ',F12.5)
 1003 format('# pixsiz ',F12.3)
 1004 format('# gain ',E15.9)
 1005 format('# offset ',F17.11)
      return
      end
