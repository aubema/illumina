c    programm to read a pixel value un a .bin file
c
c --------------------
c
c   Main programm
c
      program readvalue
c
c ----------
c
c   declare variables
c
      real dat(1024,1024)
      character*72 nomi
      integer nx,ny,nbx,nby
c
c -----------
c
c   input parameter file
c
c
      open(unit=1,file='readvalue.in',status='old')
        read(1,*) nomi
        read(1,*) nx
        read(1,*) ny
      close(unit=1)
      call twodin(nbx,nby,nomi,dat)
      print*,'Value at x=',nx,' and y=',ny,' : ',dat(nx,ny)
      open(unit=2,file='readvalue.out',status='unknown')
         write(2,100) dat(nx,ny)
      close(unit=2)
  100 format(F10.2)
      stop    
      end

