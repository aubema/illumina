       subroutine twodin(nbx,nby,filename,bindata)
c read double precision array in binary
       integer nbx,nby,i,j
       real bindata(256,256)
       character*72 filename
       open(unit=1,form='unformatted',file=filename,action='read')
         read(1) nbx,nby
         if ((nbx.gt.width).or.(nby.gt.width)) then
          print*,'You try to use a domain larger than the maximum'
          print*,'allowed. Please restrict it to no more that 256 x 256'
          print*,'Computation aborted'
          stop
         endif
         do j=nby,1,-1
            do i=1,nbx
               read(1) bindata(i,j)
            enddo
         enddo
       close(unit=1)
       return
       end
