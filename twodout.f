       subroutine twodout(nbx,nby,filename,bindata)
c write double precision array in binary
       integer nbx,nby,i,j
       real bindata(1024,1024)
       character*72 filename
       open(unit=1,form='unformatted',file=filename,action='write')
         write(1) nbx,nby
         do j=nby,1,-1
            do i=1,nbx
               write(1) bindata(i,j)
            enddo
         enddo
       close(unit=1)
       stop
       end
