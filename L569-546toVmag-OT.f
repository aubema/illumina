c program to transform the 569 and 546 radiances into magnitude V
c for OT
c    Copyright (C) 2011  Martin Aube
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
       real m,mn,mz,mnz,la(25),la1(25),la2(25),laz,el1(25),el2(25)
       real az1(25),az2(25),laz1,laz2,az(25),el(25)
       real slope569,slope546,wt546,wt569,FWHM,ele,azi,l
       integer i,j
c zenith after midnight V magnitude measured by Javier
       mz=21.40
c natural magnitude     
       mnz=21.9
c model calibration constant for each line       
       slope546=1.26E10
       slope569=7.51E9
c FWHM of V filter
       FWHM=88.
c response of the V filter
       wt546=exp(-1.*(551.-546.)**2./2./(FWHM/2.35)**2.)
       wt569=exp(-1.*(569.-551.)**2./2./(FWHM/2.35)**2.)
c we assume that mnz is constant all over the sky
       mn=mnz
       open(unit=1,file='input546.txt',status='unknown')
       do i=1,25
         read(1,*) el1(i),az1(i),la1(i)
       enddo
       close(unit=1)
       open(unit=1,file='after546.txt',status='unknown')
       do i=1,25
         read(1,*) ele,azi,l
         if (ele.eq.90.) laz1=l
       enddo
       close(unit=1)       
       open(unit=1,file='input569.txt',status='unknown')
       do i=1,25
         read(1,*) el2(i),az2(i),la2(i)
       enddo
       close(unit=1)
       open(unit=1,file='after569.txt',status='unknown')
       do i=1,25
         read(1,*) ele,azi,l
         if (ele.eq.90.) laz2=l
       enddo
       close(unit=1)
       laz=laz1*slope546*wt546+laz2*slope569*wt569
       open(unit=2,file='output.txt',status='unknown')
       do i=1,25
        do j=1,25
         if ((az2(j).eq.az1(i)).and.(el2(j).eq.el1(i))) then
           az(i)=az1(i)
           el(i)=el1(i)
           la(i)=la1(i)*slope546*wt546+la2(j)*slope569*wt569
c a fit to the elevation angle vs airglow relationship from Benn and 
c Ellison, 1998
           mn=mnz-(1.37029-0.0344429*el1(i)+0.000221429*el1(i)**2.)           
           m=-2.5*log10(la(i)/laz*(10.**(-mz/2.5)-10.**(-mnz/2.5))+
     +     10.**(-mn/2.5))
           write(2,*) el1(i),az1(i),m
         endif
        enddo
       enddo
       close(unit=2)
       stop
       end

