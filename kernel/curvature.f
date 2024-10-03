c  adding earth curvature to the topography
c
      subroutine curvature(distc,hcur)
      real*8 Rearth,distc,hcur
      Rearth=6371000.D0  
      hcur=Rearth-dsqrt(Rearth**2.-distc**2.)
      hcur=-1.*hcur
      return
      end
