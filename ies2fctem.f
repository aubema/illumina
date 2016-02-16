c
c  programme pour moyenner le ies tabulee toutes les directions
c  azimutales apres avoir incline la fonction (inclinaison du luminaire)
c
c  la matrice fctpp(i,j) est la fonction de photometrie inclinee
c  l'indice i est relie a l'angle zenithal de sorte que i=1 correspond
c  au nadir, i=90 au zenith
c  l'indice j est l'azimuth de sorte que j=1 correspond a la direction Est
c  et j=90 a la direction Ouest (ici nous supposons que la direction de
c  l'eclairage est l'est
c----------------------------------------------------------------------

c
c    Copyright (C) 2005  Martin Aube
c
c    This program is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program; if not, write to the Free Software
c    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c----------------------------------------------------------------------
      integer i,headend,nligne,ncol,ntheta,nphi
      integer ii,iip,jj,jjp,flag,ilimi,ilims,jlimi,jlims
      real tilt,fct(90,90),dthetas,thetas0,ies(150,150),phi(150)
      real theta(150),moye(150),moyo(150),fctp(90,90),fctem(180)
      real uplight,downlight,pi,totlight,gap,thetas,azim,azim0
      real thetasp,azimp,fctpp(90,90),angle(90),iesmax
      real epsi,epsilon,m,b
      character bidon,bidon1
      character*64 nom,outfile
      thetas0=1.
      azim0=1.
      dthetas=2.
      thetas=0.
      pi=3.14159265359
      do i=1,90
      do j=1,90
         fct(i,j)=-1.
         fctp(i,j)=-2.
         fctpp(i,j)=-3.
      enddo
      enddo
      do i=1,150 
         phi(i)=0.
         theta(i)=0.
         moyo(i)=0.
         moye(i)=0.
         do j=1,150
         ies(i,j)=0.
         enddo
      enddo
       print*,'ies2fctem - An application to resample-tilt-and-integrate
     + a IES photometric file'
       print*,'CopyRight Martin Aube, MEMO Environnement, 2005'
       print*,'======================================================='
       open (unit=1,file='ies2fctem.in',status='old')
       print*,'Reading parameter file ies2fctem.in...'
          read(1,*) nom
          read(1,*) tilt
          read(1,*) outfile
       close (unit=1)
       tilt=-tilt*pi/180.
       open (unit=2,file=nom,status='old')
         print*,'Reading ies tab file...'
         read(2,*) ntheta, nphi
         read(2,*) (theta(i),i=1,ntheta)
         read(2,*) (phi(j),j=1,nphi)
         read(2,*) ((ies(i,j),i=1,ntheta),j=1,nphi)
       close (unit=2)   
c
c  determining the maximum value
c 
       print*,'Searching the maximum...'    
       iesmax=0.
       do i=1,ntheta
         do j=1,nphi
            if (ies(i,j).gt.iesmax) iesmax=ies(i,j)
         enddo
       enddo
       print*,'Maximum value=',iesmax
       print*,'Angle vertical maximal du fichier ies=',theta(ntheta)
c
       gap=(theta(ntheta)-theta(1))/real(ntheta)
       print*,'Resampling function to 2 deg resolution...'
       do jj=1,90
        azim=real(jj)*dthetas-azim0
        do ii=1,90
         thetas=real(ii)*dthetas-thetas0
        if (thetas.le.theta(ntheta)) then
         dmin=100000.         
         do j=1,nphi
          do i=1,ntheta
             dist=sqrt((theta(i)-thetas)**2.+(phi(j)-azim)**2.)
             if (dist.lt.2.4*gap) then
               if (dist.lt.dmin) then
                  dmin=dist
                  fct(ii,jj)=ies(i,j)
                  if (dist.lt.dthetas) goto 400
               endif
             endif
          enddo
         enddo
        endif
 400    enddo
       enddo
c
c  inclinaison de la fonction
c      
      if (tilt.ne.0.) then
        print*,'Tilting function by ',-tilt,' rad...'
        do jj=1,90
          azim=real(jj)*dthetas-azim0
          azim=azim*pi/180.
          do ii=1,90
            thetas=-91.+real(ii)*dthetas+thetas0
            thetas=thetas*pi/180.
            azimp=atan(-1.*sin(-azim)/(sin(tilt+pi/2.)*cos(azim)-
     +      cos(tilt+pi/2.)*tan(thetas)))
            azimp=azimp*180./pi
                 
            if (azimp.lt.0.) azimp=azimp+180.
            

            thetasp=asin(cos(tilt+pi/2.)*cos(azim)*cos(thetas)+
     +      sin(tilt+pi/2.)*sin(thetas))
            thetasp=thetasp*180./pi
            if (thetasp.gt.90.) print*,'thetasp gt 90',thetasp
            if (thetasp.lt.-90.) print*,'thetap lt 90',thetasp
           
            jjp=nint((azimp+azim0)/dthetas)
            iip=nint((thetasp+91.-thetas0)/dthetas) 
            if ((iip.le.90).and.(iip.ge.1).and.(jjp.le.90).and.
     +      (jjp.ge.1)) then
               fctp(iip,jjp)=fct(ii,jj)
            endif
          enddo
        enddo  

       print*,'Resampling the tilted function...'
        do jj=1,90
        azim=real(jj)*dthetas-azim0
        do ii=1,90
            flag=0       
         thetas=real(ii)*dthetas-thetas0
         if (fctp(ii,jj).eq.-1.) fctp(ii,jj)=0.
         if (fctp(ii,jj).lt.-1.) then  
 
           ilimi=ii-int(gap/2.)-1
           if (ilimi.lt.1) ilimi=1
           ilims=ii+int(gap/2.)+1 
           if (ilims.gt.90) ilims=90                  
           jlimi=jj-int(gap/2.)-1
           if (jlimi.lt.1) jlimi=1           
           jlims=jj+int(gap/2.)+1
           if (jlims.gt.90) jlims=90 
               dmin=100000.                          
           do j=jlimi,jlims
             do i=ilimi,ilims

               dist=sqrt((theta(i)-thetas)**2.+(phi(j)-azim)**2.)
               if (fctp(i,j).gt.0.) then
                   if (dist.lt.dmin) then
                     flag=1
                     dmin=dist
                     fctpp(ii,jj)=fctp(i,j)
                   endif
               endif   
             enddo
           enddo
           if (flag.eq.0) fctpp(ii,jj)=0.
         else
           fctpp(ii,jj)=fctp(ii,jj)
         endif
         
       enddo
       enddo  
      else
         do i=1,90
            do j=1,90
               fctpp(i,j)=fct(i,j)
               if (fctpp(i,j).lt.0.) fctpp(i,j)=0.
            enddo
         enddo          
      endif                                                              ! fin condition angle pas egal a 0          
c
c  moyenner horizontalement
c
       print*,'Horizontal averaging...'
       do i=1,90
       do j=1,90
          fctem(i)=(fctpp(i,j))+fctem(i)
       enddo
       enddo
c
c calculer le uplight et downlight
c
       print*,'Integrating uplight fraction...'
       uplight=0.
       downlight=0.
       do i=1,90
          thetas=180.-real(i)*dthetas+thetas0
          angle(i)=thetas
          if (i.le.45) then
            downlight=downlight+fctem(i)*2.*pi*sin(thetas*pi/180.)
     +      *dthetas
          else
            uplight=uplight+fctem(i)*2.*pi*sin(thetas*pi/180.)
     +      *dthetas
          endif
       enddo
       totlight=uplight+downlight
       uplight=uplight/totlight
       downlight=downlight/totlight   
       open(unit=3,file=outfile,status='unknown')
          iii=89
c reechantillonner a tous les deg.
          do i=1,181
              epsi=real(i-1)*0.5
              epsilon=real(i-1)*1.

              ii=iii-int(epsi)+1
              if (angle(ii).le.epsilon) then
                m=(fctem(ii+1)-fctem(ii))/(angle(ii+1)-angle(ii))
                b=fctem(ii+1)-m*angle(ii+1)         
              else
                m=(fctem(ii)-fctem(ii+1))/(angle(ii)-angle(ii+1))
                b=fctem(ii)-m*angle(ii) 
              endif 
              write(3,*) m*epsilon+b,epsilon
          enddo
       close(unit=3)
       write(*,100) uplight
       write(*,101) downlight
       open (unit=2,file='intflux.tmp',status='unknown')
          write(2,102) uplight
          write(2,102) downlight
       close(2)
c      
c  create a pgm image of the ies tilted function
c 
c   normaliser la matrice
c
       do i=1,90
         do j=1,90
              fctpp(i,j)=fctpp(i,j)*255./iesmax
         enddo
       enddo
       open(unit=5,file='iestilt.pgm',status='unknown')      
          write (5,103)
          write (5,*) '90 90'
          write (5,*) '255'
          write (5,*) ((nint(fctpp(i,j)),j=90,1,-1),i=90,1,-1)
       close(5)
  100  format('Upward light fraction  =',1x,F5.3)
  101  format('Downward light fraction=',1x,F5.3)
  102  format(1x,F5.3)
  103  format('P2')
       stop
       end
