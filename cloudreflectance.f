c **      http://cegepsherbrooke.qc.ca/~aubema/index.php/Prof/IllumEn?action=download&upname=intensite_lumineuse.pdf  **
c **                                                                                                                  **
c **********************************************************************************************************************
c   
c    Copyright (C) 2015 Martin Aube
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
c2345678901234567890123456789012345678901234567890123456789012345678901234
c
      subroutine cloudreflectance(angzen,cloudt,rcloud)
c
c=======================================================================
c     Variables declaration
c=======================================================================
      integer cloudt
c  fitted parameters for the cloud reflectance as a function of the incident zenith angle
c  rho(z)=a0+a1*cos z + a2 * cos^2 z + a3 * cos^3 z according to Shapiro 1982 Table 10
c cloud no 1 are thin cirrus & cirrostratus
c cloud no 2 are thick cirrus & cirrostratus
c cloud no 3 are altostratus & altocumulus
c cloud no 4 are stratocumulus & stratus
c cloud no 5 are cumulus & cumulonimbus
      real rhocld(5,4)
      real angzen,rcloud                                               

      rhocld(1,1)=0.25674
      rhocld(1,2)=-0.18077
      rhocld(1,3)=-0.21961
      rhocld(1,4)=0.25272
      rhocld(2,1)=0.60540
      rhocld(2,2)=-0.55142
      rhocld(2,3)=-0.23389
      rhocld(2,4)=0.43648
      rhocld(3,1)=0.66152
      rhocld(3,2)=-0.14863
      rhocld(3,3)=-0.08193
      rhocld(3,4)=0.13442
      rhocld(4,1)=0.71214
      rhocld(4,2)=-0.15033
      rhocld(4,3)=0.00696
      rhocld(4,4)=0.03904
      rhocld(5,1)=0.67072
      rhocld(5,2)=-0.13805
      rhocld(5,3)=-0.10895
      rhocld(5,4)=0.09460
       rcloud=rhocld(cloudt,1)+rhocld(cloudt,2)*cos(angzen)
     + +rhocld(cloudt,3)*(cos(angzen))**2.+rhocld(cloudt,4)*
     + (cos(angzen))**3.
c        print*,rcloud,angzen,cos(angzen)
      return
      end

