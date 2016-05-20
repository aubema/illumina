c      program AllSkyRatio
c
      integer i,n,cnt,lennom1,lennom2,nb1,nb2
      real az1,az2,el1,el2,rad2,rad1,cf,nano
      character*60 file1,file2,outname
      nano=1e-9
c First radiance txt file name nb_lines
      read*,file1,nb1
      lennom1=index(file1,' ')-1
c Second radiance txt file name nb_lines
      read*,file2,nb2
      lennom2=index(file2,' ')-1
c Calibration factor
      read*,cf
c Output txt file name
      read*,outname
      open(unit=1,file=file1,status='old')
      open(unit=2,file=file2,status='old')
      open(unit=9,file=outname,status='unknown')
      
      do nn=1,nb1
c         print*,nn
         read(1,*) el1,az1,rad1
         rad1=rad1*cf
         do n=1,nb2
            read(2,*) el2,az2,rad2
            if (nint(el1).eq.nint(el2)) then
               if (nint(az1).eq.nint(az2)) then
                  rad2=rad2*cf
c                  print*,el1,az1,rad1,rad2,rad1/rad2
                  write(9,*) el1,az1,rad1/rad2
               endif
            endif
         enddo
         rewind(2)
        enddo
      print*,'End of computations'
      close(unit=9)
      close(unit=2)
      close(unit=1)
      stop
      end
