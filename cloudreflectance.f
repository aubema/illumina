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
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=50)

      integer cloudt
c  fitted parameters for the cloud reflectance as a function of the incident zenith angle
c  rho(z)=A1+A2*exp(-1.*(cos(z)/A3)**2.) 
      real rhocld(5,3)
      real angzen,rcloud                                               

      rhocld(1,1)=0.1
      rhocld(1,2)=0.146
      rhocld(1,3)=0.357
      rhocld(2,1)=0.24
      rhocld(2,2)=0.3327
      rhocld(2,3)=0.344
      rhocld(3,1)=0.56
      rhocld(3,2)=0.0932
      rhocld(3,3)=0.34
      rhocld(4,1)=0.515
      rhocld(4,2)=0.144
      rhocld(4,3)=0.49
      rhocld(5,1)=0.61
      rhocld(5,2)=0.0907
      rhocld(5,3)=0.4087
      rcloud=rhocld(cloudt,1)+rhocld(cloudt,2)*exp(-1.*(cos(angzen)/
     +rhocld(cloudt,3))**2.)
c        print*,rcloud,angzen,cos(angzen)
      return
      end

