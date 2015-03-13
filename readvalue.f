c    programme lire une valeur dans un pgm

c
c -----------------
c   identification des variables 
c
c
c --------------------
c
c   programme principal
c
      program readvalue
c
c ----------
c
c   declaration des variables
c
      real xcell0,ycell0
      real pixsiz,dat(1024,1024)
      character*72 nomi
      character*12 nom
      integer nx,ny,nbx,nby
c   
c ----------
c
c   initialisation des variables
c
c
c -----------
c
c   fichier de parametres
c
c
c   choix du nom de la racine de fichiers
c
      open(unit=1,file='readvalue.in',status='old')
        read(1,*) nomi
        read(1,*) nx
        read(1,*) ny
      close(unit=1)
c
 
       nom='input data '
       call intrants2d(nomi,dat,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)      
       print*,'Value at x=',nx,' and y=',ny,' : ',dat(nx,ny)

       open(unit=2,file='readvalue.out',status='unknown')
          write(2,100) dat(nx,ny)
       close(unit=2)
  100  format(F10.2)
       stop    
       end

