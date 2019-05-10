c     find the max height
      integer height
      parameter (height=100) 
      real cthick(height)                                                 ! voxel thickness array (meter)
      real cellh(height)                                                  ! voxel height array (meter)
      call verticalscale(cthick,cellh)
      open(unit=1,file='maxheight.tmp',status='unknown')
      write(1,*) cellh(height)
      close(unit=1)
      stop
      end
