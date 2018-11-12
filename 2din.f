       subroutine 2din(nbx,nby,filename,bindata)
c read double precision array in binary
       integer nbx,nby,i,j
       real bindata(1024,1024)
       character*72 filename
       open(unit=1,form='unformatted',file=filename,action='read')
         read(1) nbx,nby
         do j=nby,1,-1
            do i=1,nbx
               read(1) bindata(i,j)
            enddo
         enddo
       close(unit=1)
       stop
       end
