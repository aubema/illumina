       subroutine twodin(nbx,nby,filename,bindata,width,height)
       integer width,height            ! Matrix dimension in Length/width and height
       CHARACTER*72 filename
! read double precision array in binary
       integer nbx,nby,i,j
       real, dimension(:,:), allocatable ::  bindata
       allocate ( bindata(width,height) )
       open(unit=1,form='unformatted',file=filename,action='read')
         read(1) nbx,nby
         if ((nbx.gt.width).or.(nby.gt.height)) then
          print*,'You try to use a domain larger than the maximum'
          print*,'allowed. Please restrict it to no more that 512 x 512'
          print*,'Your domain size is: ',nbx,'x',nby
          print*,'Computation aborted'
          stop
         endif
         do j=nby,1,-1
            do i=1,nbx
               read(1) bindata(i,j)
            enddo
         enddo
       close(unit=1)
       deallocate ( bindata )
       return
       end
