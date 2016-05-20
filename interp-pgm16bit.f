c    programme pour degrader la resolution d'un pgm

c
c -----------------
c   identification des variables 
c
c
c --------------------
c
c   programme principal
c
      program interp
c
c ----------
c
c   declaration des variables
c
      real xcell0,ycell0,npix
      real pixsiz,tau(1024,1024),tauout(1024,1024)
      real taumax
      real gain,offset
      character*72 nomo,nomi
      character*12 nom
      integer i,j,factor,nx,ny,nbx,nby,side,valmax
c   
c ----------
c
c   initialisation des variables
c
       taumax=0.
       valmax=65535
c
c -----------
c
c   fichier de parametres
c
c
c   choix du nom de la racine de fichiers
c
      print*,'Name of the pgm file ?'  
      read*,nomi
      print*,'Output resolution (pixels)?'
      read*,factor  
      print*,'Name of the output pgm file ?'  
      read*,nomo     
c
       call intrants2d(nomi,tau,xcell0,ycell0,pixsiz,nbx,nby)      
       side=factor/2
       print*,'Averaging radius=',side
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
             if (tauout(nx,ny).gt.taumax) then
            
                taumax=tauout(nx,ny)
             endif
 114       continue
 113    continue
c
c ----------
c
c   fabrication d'un nouveau fichier pgm
c
        gain=taumax/real(valmax)
        offset=0.
        nom='average'
        call extrants2d (nomo,tauout,nom,xcell0,ycell0,pixsiz,
     +  gain,offset,nbx,nby,valmax)
       stop    
       end
