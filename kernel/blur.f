c    programm to blur a bin file

c
c -----------------
c   identification des variables 
c
c
c --------------------
c
c   programme principal
c
      program blur
c
c ----------
c
c   declaration des variables
c
      real npix
      real tau(512,512),tauout(512,512)
      character*72 nomo,nomi
      integer i,j,nx,ny,nbx,nby,side
c
c -----------
c
c   fichier de parametres
c
c
c   choix du nom de la racine de fichiers
c
      print*,'Name of the bin file to blur ?'  
      read*,nomi
      print*,'Bluring radiux (pixels) ?'
      read*,side  
      print*,'Name of the output bin file ?'  
      read*,nomo     
      call twodin(nbx,nby,nomi,tau)   
c    
c ---------
c
c   Interpoler la carte d epaisseur optique
c
      print*,'Interpolating map...'
      do 113 nx=1,nbx
         do 114 ny=1,nby
            if (((nx.gt.side).and.(nx.lt.nbx-side)).and.
     +((ny.gt.side).and.(ny.lt.nby-side))) then
               tauout(nx,ny)=0.
               npix=0.
               do 115 i=nx-side,nx+side
               do 116 j=ny-side,ny+side
                  npix=npix+1.
                  tauout(nx,ny)=tau(i,j)+tauout(nx,ny)
 116           continue
 115           continue
               tauout(nx,ny)=tauout(nx,ny)/npix
             else
               tauout(nx,ny)=tau(nx,ny)
             endif
 114       continue
 113    continue
c
c ----------
c
c   fabrication d'un nouveau fichier bin
c
        call twodout(nbx,nby,nomo,tauout)
       stop    
       end
