c routine d ecriture des intrants 2d                                                                
c                                                                                                  
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
         subroutine write2d (outfile,valeur,lat0,lon0,pixsiz,nbx,nby,                              
     +   valmax,gain,offset)  
         integer width                                                                      
         parameter (width=1024)                                                                    
         real valeur(width,width),lat0,lon0,gain,offset                                            
         integer i,j,nbx,nby,valmax                                                                
         character*40 outfile                                                                      
c definition des variables                                                                          
c width=taille de la matrice reservee pour le calcul                                                
c lat0=latitude du pixel du coin inferieur gauche                                                  
c lon0=longitude du pixel du coin inferieur gauche                                                  
c gain= pente de la fonction lineaire qui permet de passer de la                                    
c       valeur numerique a la valeur physique                                                      
c offset= ordonnee a l origine de la fonction lineaire qui                                          
c         permet de passer de la valeur numerique a la valeur physique                              
c nbx= nombre de pixel horizontaux                                                                  
c nby= nombre de pixels verticaux                                                                  
c valmax= valeur numerique maximale (normalement on utilises soit                                  
c         65535 ou 255                                                                              
c outfile= nom du fichier de sortie
c
c         print*,'Ecriture du fichier de ',outfile
 1000    format('P2')
 1001    format('# lat0 ',F10.4)
 1002    format('# lon0 ',F10.4)
 1003    format('# pixsiz ',F10.4)
 1004    format('# gain ',F15.4)
 1005    format('# offset ',F15.4)
         open(unit=1,file=outfile,status='unknown')
           write(1,1000)
           write(1,1001) lat0
           write(1,1002) lon0
           write(1,1003) pixsiz
           write(1,1004) gain
           write(1,1005) offset
           write(1,*) nbx,nby,valmax
           write(1,*) ((nint(valeur(i,j)),i=1,nbx),j=nby,1,-1)
           return
           end
