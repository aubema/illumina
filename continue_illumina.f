c programme pour faire la somme des flux dans un fichier .out de illumina
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

         character*60 name
	 integer i
	 real value, flux
         open (unit=1,file='continue_illumina.in',status='old')
            read(1,*) name
         close (unit=1)
         name=name//'.out'
         open(unit=2,file=name,status='old')
	   flux=0.
           do i=1,1500
              read(2,*,end=10) bidon,bidon,bidon,bidon,value
	      flux=flux+value
           enddo
  10     close (unit=2)
          open (unit=1,file='continue_illumina.out',status='unknown')
            write(1,*) flux, 'Total flux for experiment ',name
         close (unit=1)     
	 stop
	 end   
