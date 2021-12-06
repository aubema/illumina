c      program CompareToMiro
c
      integer i,n,cnt,lennom
      real az,z,el,radmiro,radaube,cf,nano
      character*60 mironame,aubename,outname,newmiro,mot,newaube
      nano=1e-9
      print*,'Miro file name?'
      read*,mironame
      lennom=index(mironame,' ')-1
      print*,'Aube file name?'
      read*,aubename
      print*,'Calibration factor?'
      read*,cf
      newmiro=mironame(1:lennom)//'miro.txt'
      newaube=mironame(1:lennom)//'aube.txt'
      outname=mironame(1:lennom)//'diff.txt'
      print*,newmiro
      open(unit=1,file=mironame,status='unknown')
      open(unit=2,file=aubename,status='old')
      open(unit=8,file=newmiro,status='unknown')
      open(unit=9,file=outname,status='unknown')
      open(unit=10,file=newaube,status='unknown')
      
      do nn=1,200
      read(2,*,end=100) el, aza,radaube
      radaube=radaube*cf
      
      mot='bidon'
      n=0
      do while (mot.ne.'data')
      n=n+1
      read(1,*) mot
      enddo
      cnt=0
      do i=1,1387
     
        read(1,*) z,az,radmiro
        if (nint(z).eq.nint(90.-el)) then
          if (nint(abs(az-aza)).lt.3) then
            print*,el,aza,radaube,radmiro,radaube-radmiro
            write(8,*) el, aza, radmiro*nano/2.
            write(9,*) el,aza,radaube-radmiro*nano/2.
            write(10,*) el,aza,radaube
          endif
        endif

      enddo
      rewind(1)
      enddo


 100  print*,n
      close(unit=10)
      close(unit=9)
      close(unit=8)
      close(unit=2)
      close(unit=1)
      stop
      end
      
