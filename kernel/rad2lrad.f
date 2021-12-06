c conversion des fichier all sky el, az, radiance en x y log(radiance)
       integer n,nl
       real lrad,rad,x,y,el,az
       character*72 outname
       open(unit=3,file='nligne.tmp',status='old')
            read(3,*) nl
            read(3,*) outname
       close(unit=3)
       open (file='rad2lograd.tmp',unit=1,status='unknown')
       open(unit=2,file=outname,status='unknown')
       do n=1,nl
          read(1,*) el,az,rad
          lrad=log10(rad)
          write(2,*) el,az,lrad
       enddo
       close(unit=1)
       close(unit=2)
       stop
       end

