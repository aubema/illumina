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
      subroutine cloudtransmitance(angzen,cloudt,tcloud)
c
c=======================================================================
c     Variables declaration
c=======================================================================
      integer cloudt
c  fitted parameters for the cloud reflectance as a function of the incident zenith angle
c  rho(z)=b0+b1*cos z + b2 * cos^2 z + b3 * cos^3 z according to Shapiro 1982 Table 11
      real thocld(5,4)
      real angzen,tcloud                                               

      thocld(1,1)=0.63547
      thocld(1,2)=0.35229
      thocld(1,3)=0.08709
      thocld(1,4)=-0.22902
      thocld(2,1)=0.26458
      thocld(2,2)=0.66829
      thocld(2,3)=0.24228
      thocld(2,4)=-0.49357
      thocld(3,1)=0.19085
      thocld(3,2)=0.32817
      thocld(3,3)=-0.08613
      thocld(3,4)=-0.08197
      thocld(4,1)=0.13610
      thocld(4,2)=0.29964
      thocld(4,3)=-0.14041
      thocld(4,4)=0.00952
      thocld(5,1)=0.17960
      thocld(5,2)=0.34855
      thocld(5,3)=-0.14875
      thocld(5,4)=0.01962
       tcloud=thocld(cloudt,1)+thocld(cloudt,2)*cos(angzen)
     + +thocld(cloudt,3)*(cos(angzen))**2.+thocld(cloudt,4)*
     + (cos(angzen))**3.
c        print*,rcloud,angzen,cos(angzen)
      return
      end

