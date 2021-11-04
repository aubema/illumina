c  adding earth curvature to the topography
c
      subroutine molabs(lambda,bandw,tabs)
      real lambda,bandw,tabs,wmax,wmin,wla,ta,bwa,bwabs
      wmin=lambda-bandw/2.
      wmax=lambda+bandw/2.
      tabs=0.
      wlabs=0.
      wla=0.
      open(unit=1,file='MolecularAbs.txt',status='old')
        read(1,*) 
        do while (wla.lt.wmax) 
           read(1,*) wla,ta,bwa
           if (wla.gt.wmin) then
              tabs=tabs+tabs*bwa
              bwabs=bwabs+bwa
           endif
        enddo
        tabs=tabs/bwabs
      close(unit=1)
      return
      end
