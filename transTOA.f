c calculation of the total vertical transmittance of the atmosphere
c (aerosols and molecules separately
c
c pour utilisation avec Illumina 
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2021  Martin Aube
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

      subroutine transtoa(lambm,bandw,taua,layaod,pressi,tranam,tranaa
     +,tranal)
      real m,lambm,taua,pressi,tranam,tranaa,tranal,layaod,tabs,bandw
      m=(pressi/101.3)
c  transmittance tiree de Kneizys et al. (1980)   
      tranam=exp(-1.*m/(((lambm/1000.)**4.)*(115.6406-(1.335/((lambm/
     +1000.)**2.)))))
      call molabs(lambm,bandw,tabs)
      tranam=tranam*tabs
      tranaa=exp(-1.*taua)
      tranal=exp(-1.*layaod)
      return
      end
