c  adding earth curvature to the topography
c
      subroutine curvature(distc,hcur)
      real Rearth,distc,hcur
      Rearth=6371000.  
      hcur=Rearth-sqrt(Rearth**2.-distc**2.)
      hcur=-1.*hcur
      return
      end
