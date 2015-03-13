c conversion des fichier all sky el, az, radiance en coord polaires
       integer n,nl
       real rad,x,y,el,az,radmax,radmin
       open(unit=3,file='nligne.tmp',status='old')
            read(3,*) nl
       close(unit=3)
       open (file='rad2polar.tmp',unit=1,status='unknown')
c       open (file='rad2polarscale.tmp',unit=4,status='unknown')
       open(unit=2,file='data.polar',status='unknown')
       radmax=-1.E19
       radmin=1.E19
       do n=1,nl
          read(1,*) el,az,rad
          if (rad.gt.radmax) radmax=rad
          if (rad.lt.radmin) radmin=rad
          x=(90.-el)*COS(3.14159*(90.-az)/180.)
          y=(90.-el)*sin(3.14159*(90.-az)/180.)
          write(2,*) x,y,rad
       enddo
       
c       write(4,*) nint(radmin-1), nint(radmax+1), (nint(radmax+1)-
c     + nint(radmin-1))/12.
c       close(unit=4)
       close(unit=1)
       close(unit=2)
       stop
       end

