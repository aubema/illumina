c     find the max height
      integer height
      parameter (height=1024) 
      real cthick(height)                                                 ! voxel thickness array (meter)
      real cellh(height)                                                  ! voxel height array (meter)
      call verticalscale(dx,cthick,cellh)
      open(unit=1,file='maxheight.tmp',status='unknown')
      write(1,*) cellh(height)
      close(unit=1)
      stop
      end
