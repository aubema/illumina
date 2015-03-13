      program PuissanceLumineuse

c déclaration des variables
      parameter(nmax=1024,mmax=1024)
      real sx,sy,pixsize,pixmin,a1,a2,b1,b2
      real e1,e2,p,p1max,p2max,hwhm1,hwhm2,sig1,sig2,cte
      real ampbruit,rmin,rmax,r,xx,yy,Pi
      integer x,y,x0,y0,n,m,s1,s2,tire,k
      integer valmax
      real puissance(nmax,mmax),points(10 000,4),bruit(nmax,mmax)
      real lat0,lon0,offset,gain
      character*40 outfile

c ouverture du fichier de données sortantes
      open(unit=1, file='PLum.in', status='unknown')

c définiton des variables

c réels:

c sx,sy 	: dimension en x/y du terrain étudié
c pixsize 	: dimension du coté d'une case carré
c pixmin	: dimension de case minimale
c a1,b1		: demi grand et petit axes du centre-ville
c a2,b2		: demi grand et petit axese de la ville
c e1		: excentricité du centre-ville
c e2		: excentricité de la ville
c p		: puissance lumineuse totale
c p1max 	: p. lum. max. du centre-ville
c p2max 	: p. lum. max. de la ville
c hwhm1		: demie-largeur à mi hauteur gauss 1
c hwhm2		: demie-largeur à mi hauteur gauss 2
c sig1		: écart type gauss 1
c sig2		: écart type gauss 2
c cte		: valeur de la constante
c ampbruit	: fonctions gaussiennes radiales des points-bruit
c rmin, rmax	: rayon max des taches sombres ou clair
c r		: rayon entre point-bruit et centre
c Pi		: Duh (Que J'Aime À Faire Apprendre Un Nombre Utile Aux....)
c puissance()	: matrice
c points()	: coordonnées des régions sombres/claires
c bruit()	: matrice-bruit

c entiers:

c x,y		: coordonnées des cases dans la matrice
c x0,y0		: coordonnées du centre de la ville
c n,m,		: entiers de boucles
c s1		: switch poser une géométrie elliptique du centre-ville
c s2		: switch poser une géométrie elliptique de la ville
c tire		: nombre de points de bruit à tirer

c initialisation des variables2
      sx=0.
      sy=0.
      pixsize=0.
      pixmin=0.
      a1=0.
      a2=0.
      b1=0.
      b2=0.
      e1=0.
      e2=0.
      p=0.
      p1max=0.
      p2max=0.
      whm1=0.
      whm2=0.
      sig1=0.2
      sig2=0.
      cte=0.
      ampbruit=0.
      rmin=0.
      rmax=0.
      r=0.
      Pi=3.1415926535
      x=0
      y=0
      x0=0
      y0=0
      n=0
      m=0
      s1=0
      s2=0
      tire=0

c lecture des données
      read(1,*) sx, sy
      read(1,*) x0,y0
      read(1,*) pixsize
      read(1,*) s1, b1
      read(1,*) s2, b2
      read(1,*) p1max,hwhm1
      read(1,*) p2max,hwhm2
      read(1,*) cte
      read(1,*) rmin,rmax
      read(1,*) ampbruit
      read(1,*) tire

c détermine la dimension minimale de case
      if(sx.gt.sy) then
	pixmin=Sx/1024.
      else
	pixmin=Sy/1024.
      endif

c erreur si taille de case trop petite
      if(pixsize.lt.pixmin)then
	WRITE(*,*) 'pixsize doit etre plus grande que ', pixmin
	stop
      endif

c détermination du nombre de cases de la matrice
      n=int(sx/pixsize)+1
      m=int(sy/pixsize)+1

c initialisation de la matrice
      do x=1, n
	do y=1, m
	  puissance(x,y)=-1.
	  bruit(x,y)=0.
	enddo
      enddo

c calcul des sigma
      sig1=2.*hwhm1/2.35482
      sig2=2.*hwhm2/2.35482

c détermination de a1 et b1
      a1=1.
      if(s1.eq.0)then
	b1=a1
      endif

c détermination de a2 et b2
      a2=1.
      if(s1.eq.0)then
	b2=a2
      endif

c séléctionne l'axe le plus grand pour le calcul de e1
      if (a1.gt.b1) then
        e1=(sqrt(a1**2.-b1**2.))/a1 
      else
        e1=(sqrt(b1**2.-a1**2.))/b1
      endif

c séléctionne l'axe le plus grand pour le calcul de e2
      if (a2.gt.b2) then
        e2=(sqrt(a2**2.-b2**2.))/a2 
      else
        e2=(sqrt(b2**2.-a2**2.))/b2
      endif

c purge de la routine aléatoire
      open(unit=3, file='purgatoire.rdm', status='unknown')
      read(3,*)x
      if(x.gt.666)then
	x=1
	rewind(unit=3)
	write(3,*)x
      endif
      do i=1,x
	call random_number(xx)
      enddo
      rewind(unit=3)
      write(3,*)x+1
      close(unit=3)

c tire un certain nombre de points dans la ville
      if(tire.ne.0)then
	i=1
	do i=1,tire
	  call random_number(xx)
	  x=int(xx*real(n))+1
	  call random_number(yy)
	  y=int(yy*real(m))+1
	    points(i,1)=1000.*x+y
	    points(i,4)=0.
	enddo

c maximum des taches aléatoire
	do i=1,tire
	  call random_number(xx)
	  points(i,2)=2.*(xx-0.5)*ampbruit
	enddo

c sigma des taches aléatoires
	do i=1,tire
	  call random_number(xx)
	  points(i,3)=2.*(xx*(rmax-rmin)+rmin)/2.35482
	enddo     
      endif

      do x=1,n
	do y=1,m
c appel de la sous-routine de calcul de la luminosité
	  call LumCalc(x,x0,y,y0,a1,a2,b1,b2,e1,e2,sig1,sig2,
     + p,b,p1max,p2max,cte,tire,points)
	  puissance(x,y)=p
	  bruit(x,y)=b
	enddo
      enddo

c appel de la routine d'écriture de l'image
      outfile='PLum.pgm'
      lat0=45.
      lon0=72.
      valmax=255.
      gain=64.
      offset=0.
         call write2d (outfile,puissance,lat0,lon0,pixsize,n-1,m-1,
     +   valmax,gain,offset)

c appel de la routine d'écriture de l'image pour le bruit
      outfile='Bruit.pgm'
      lat0=45.
      lon0=72.
      valmax=255.
      gain=64.
      offset=0.
         call write2d (outfile,bruit,lat0,lon0,pixsize,n-1,m-1,
     +   valmax,gain,offset)

      close(unit=1)

      stop
      end       
