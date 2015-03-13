      subroutine GaussLisse(x,x0,y,y0,a1,a2,b1,b2,
     + e1,e2,p1max,p2max,sig1,sig2,cte,p)

c déclaration des variables
      real a1,a2,b1,b2,e1,e2,p,p1,p2,p1max,p2max,sig1,sig2,cte
      real angle,d2,r2,Pi,sigma
      integer x,y,x0,y0

c détermination de d2 (distance par rapport au centre)
      d2=(x-x0)**2.+(y-y0)**2.
c détermination de l'angle
      if(d2.ne.0.)then
	angle=atan(real(y-y0)/real(x-x0))
      else
	angle=0.
      endif

c centre-ville

c détermination du rayon
      r2=b1**2./(1.-(e1**2.)*(cos(angle)**2.))

c détermination de sigma
      sigma=sig1*sqrt(r2)/a1	  

c plug dans la gaussienne
      p1=p1max*exp(-d2/(2.*sigma**2.))

c ville
c détermination du rayon
      r2=b2**2./(1.-(e2**2.)*(cos(angle)**2.))

c détermination de sigma
      sigma=sig2*sqrt(real(r2))/a2	  

c plug dans la gaussienne
      p2=p2max*exp(-d2/(2.*sigma**2.))

      p=p1+p2+cte

      return
      end       