c routine de lecture des intrants 2d
c
c 
c   
c    Copyright (C) 2009  Martin Aube
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
         subroutine intrants2d (infile,valeur,xcell0,ycell0,
     +   pixsiz,nbx,nby)
         integer width
         parameter (width=512)
         real valeur(width,width),val,xcell0,ycell0,gain,offset
         real pixsiz
         integer hcnt,i,j,nbx,nby,valmax
         character*72 infile,tag,bidon
         xcell0=0.
         ycell0=0.
         print*,'Loading file : ',infile
         bidon='#'                                                        ! Le caractere # indique que la ligne contient une constante.
         hcnt=0                                                           ! Le nombre de lignes de constantes au debut du fichier (hcnt) est 
c                                                                         ! fixe a 0.
         offset=0.
         gain=0.
         open(unit=1,file=infile,status='old')
         read(1,*)                                                        ! Passe la premiere ligne (il y a toujours un p2).
         do i=1,12                                                        ! Debut de la boucle sur les 12 premieres lignes du fichier.
            read(1,*,end=3,err=1) bidon,tag,val                           ! Lecture du fichier relief.pgm.
            goto 2
 1         backspace 1
            read(1,*,end=3) bidon
 2         if (bidon(1:1).eq.'#') then                                    ! Si le caractere # est present sur une ligne, il continue a 
c                                                                         ! chercher une constante.
               hcnt=hcnt+1                                                ! Permet de calculer le nombre de lignes de constantes au 
c                                                                         ! debut du fichier.
               if (tag(1:4).eq.'lat0') ycell0=val                         ! Lecture de la latitude inferieure du domaine.
               if (tag(1:4).eq.'lon0') xcell0=val                         ! Lecture de la longitude inferieure du domaine.
               if (tag(1:4).eq.'gain') gain=val                           ! Lecture du gain, qui permet de multiplier les donnees par une 
c                                                                         ! constante.
               if (tag(1:6).eq.'offset') then                             ! Lecture de l'offset, qui permet d'additionner une constante aux donnees.
                 offset=val
               else
                 offset=0.
               endif
               if (tag(1:6).eq.'pixsiz') pixsiz=val
            else
               goto 3
            endif  
         enddo                                                            ! Fin de la boucle sur les 12 premieres lignes du fichier.    

 3      rewind 1
         read(1,*)                                                        ! Passe la premiere ligne (il y a toujours un p2).
         do i=1,hcnt                                                      ! Debut de la boucle sur les lignes de constantes au debut du fichier.
             read(1,*)                                                    ! Lecture de hcnt lignes pour passer les lignes de constantes.
         enddo                                                            ! Fin de la boucle sur les lignes de constantes au debut du fichier.
         read(1,*) nbx,nby,valmax                                         ! Lecture de la taille du fichier (en cases) et la valeur maximale 
c                                                                         ! des donnees.
          read(1,*) ((valeur(i,j),i=1,nbx),j=nby,1,-1)                    ! Lecture de toutes les donnees qui sont ensuite inscrites
c                                                                         ! dans la matrice alt_sol. Ce sont des boucles imbriquees dans la
c                                                                         ! fonction "read" qui couvrent tout le domaine delimite par nbx et nby.
c                                                                         ! L'increment de la boucle sur les ranges (latitute) est de -1 car 
c                                                                         ! on considere que dans les fichiers,
c                                                                         ! le nombre en haut a gauche des donnees est a la
c                                                                         ! coordonnee (1,nby).
      close(1)                                                            ! Fermeture du fichier relief.pgm.
      do i=1,nbx                                                          ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
           valeur(i,j)=valeur(i,j)*gain+offset                            ! Transformation des donnees avec le gain et l'offset et recherche 
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
      enddo 
      return
      end
