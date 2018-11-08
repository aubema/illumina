       subroutine doubleout(nbx,nby,filename,bindata)
c write double precision array in binary
       integer nbx,nby,i,j
       real*8 bindata(1024,1024)
       character*72 filename
       open(unit=1,form='unformatted',file=filename,action='write')
         do j=1,nby
            do i=1,nbx
               write(1) bindata(i,j)
            enddo
         enddo
       close(unit=1)
       stop
       end
