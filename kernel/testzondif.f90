     ! test zone_diffusion

     real*8 effet,zondif(30000000,3),siz
!     do i=1,3000000
!       do j=1,3
!         zondif(i,j)=0.
!       enddo
!     enddo
     effet=100000.
     siz=1000.
     call zone_diffusion(effet,zondif,ncell,siz)
     print*,ncell
     stop
     end
     
