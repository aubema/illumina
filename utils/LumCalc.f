      SUBROUTINE LumCalc(x,x0,y,y0,a1,a2,b1,b2,e1,e2,sig1,sig2,
     + p,b,p1max,p2max,cte,tire,points)


c déclaration des variables
      REAL a1,a2,b1,b2,e1,e2,p,b,p1,p2,p1max,p2max,sig1,sig2,cte,sigma
      REAL angle,d2,r2,Pi,points(10 000,4)
      INTEGER x,y,xb,yb,x0,y0,tire
      
      
c définiton des variables

c réels:

c a1,b1		: demis grand et petit axes du centre-ville
c a2,b2		: demis grand et petit axese de la ville
c e1		: excentricité du centre-ville
c e2		: excentricité de la ville
c p		: p.lum. totale
c b		: bruit
c p1 		: p. lum. de la gaussienne 1 (centre-ville)
c p2 		: p. lum. de la gaussienne 2 (ville)
c p1max 	: p. lum. max. du centre-ville
c p2max 	: p. lum. max. de la ville
c sig1		: écart type gauss 1
c sig2		: écart type gauss 2
c cte		: valeur de la constante
c Pi		: Duh (3.1415926535)
c sigma 	: écart-type variable
c angle		: coordonnée angulaire du point
c d2		: distance entre le point et le centre de l'ellipse au carré
c r2		: carré du rayon de la ville
c points()	: matrice des coordonnées des centres sinus-bruit et du bruit associé

c entiers:

c x,y		: coordonnées des valeurs dans la matrice
c x0,y0		: coordonnées du centre de la ville (origine)
c xb,yb		: coordonnées des centres sinus-bruit
c tire		: nombre de centres sinus-bruit

c initialisation des variables
      p=0.
      b=0.
      p1=0.
      p2=0.
      p3=0.
      Pi=3.1415926535
      sigma=0.
      angle=0.
      d2=0.
      r2=0.
      i=0
      j=0

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



c calcul du bruit
      if(tire.ne.0)then
	do i=1,tire
c extraction des coordonnéées à partir de la matrice
	  xb=int(points(i,1)/1000.)
	  yb=points(i,1)-1000.*xb

c calcul du rayon par rapport au point-bruit
	  r2=(x-xb)**2.+(y-yb)**2.
	

c calcul de la puissance selon la fonction gaussienne
	    points(i,4)= points(i,2)*exp(-r2/(2.*points(i,3)**2.))
	enddo

c calcul de la puissance lumineuse
	p=p1+p2+cte

c ajout de la valeur du bruit   
	do i=1,tire
	  b=points(i,4)+b
	enddo
	p=b+p
	if(p.gt.255.)p=255.
	if(p.lt.0.)p=0.
      else
	p=p1+p2+cte
      endif

      return
      end