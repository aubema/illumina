c programme interpole les reflectances de modis pour chaque longueur
c d'onde nominale de la modelisation 
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
c    along with this program.  If not, see <http://www.gnu.org/licensce
c
      subroutine interp_modis()
      integer width
      parameter (width=1024)
      real modband(3),refl(width,width,3),avgwav,gain,offset,xcell0
      real ycell0,dist,dat(width,width)
      real dist1,dist2,wav1,wav2,m,b
      integer n_modis,i,n_bands,nb,ii,jj,valmax,pixsiz,n1,n2
      REAL, DIMENSION(:), ALLOCATABLE :: modis_wav
      REAL, DIMENSION(:,:), ALLOCATABLE :: bands
      character*72 mod_file,outfile,bands_file,reflex_file
      character*12 nom
      character*3 lambda
      pixsiz=1000.
      xcell0=0.
      ycell0=0.
      valmax=65535
      gain = 1./real(valmax)
      offset=0.
c modis reflectance file list file name
c format of the file
c Line 1: nfile
c Following lines: wavelength_nm(integer) modis_file_name
c
      reflex_file='modis.dat' 
c modelling spectral bands
      bands_file='integration_limits.dat'
c
c reading modis reflectance files and wavelength
c

      OPEN(UNIT=42,FILE=reflex_file,STATUS='OLD')
        READ(42,*) n_modis
        ALLOCATE(modis_wav(n_modis))
        DO i=1,n_modis
          READ(42,*) modis_wav(i),mod_file
          nom='reflexion '
          CALL intrants2d(mod_file,refl(:,:,i),xcell0,ycell0,
     +    pixsiz,nbx,nby)
        ENDDO
      CLOSE(42)
c
c reading wavelength bands list
c
      OPEN(UNIT=42,FILE=bands_file,STATUS='OLD')
        READ(42,*) n_bands
        ALLOCATE(bands(n_bands+1,2))
        DO nb=1,n_bands+1
          READ(42,*) bands(nb,1)
        enddo
        do nb=1,n_bands+1
          if (nb.le.n_bands) then
            bands(nb,2)=bands(nb+1,1)          
            avgwav=(bands(nb,1)+bands(nb,2))/2.
          else
            avgwav=700.
          endif
          write(lambda, '(I3.3)' ) int(avgwav)
          outfile='modis_'//lambda//'.pgm'
c interpolate 
          dist1=1000000.
          dist2=1000000.
          do n=1,n_modis 
c find 2 nearest modis wavelengths
            dist=abs(modis_wav(n)-avgwav)
            if (dist.lt.dist1) then
               dist1=dist
               wav1=modis_wav(n)
               n1=n
            endif
          enddo
          do n=1,n_modis
            dist=abs(modis_wav(n)-avgwav)
            if ((dist.lt.dist2).and.(dist.ge.dist1)) then
               dist2=dist
               wav2=modis_wav(n)
               n2=n
            endif
          enddo
          do ii=1,nbx
            do jj=1,nby
              m=(refl(ii,jj,n1)-refl(ii,jj,n2))/(wav1-wav2)
              b=refl(ii,jj,n1)-m*wav1
              dat(ii,jj)=m*avgwav+b
              if (dat(ii,jj).gt.1.) dat(ii,jj)=1.
              if (dat(ii,jj).lt.0.) dat(ii,jj)=0.
            enddo
          enddo
            nom='reflectance'
            call extrants2d (outfile,dat,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)
        ENDDO
      CLOSE(42)
      DEALLOCATE(modis_wav)
      DEALLOCATE(bands)
c ajouter une bande fixe a 700nm to filter water surface



      return
      end
