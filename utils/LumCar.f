        program FullCar

      real LumTot,LumMoy,ET,Somme,LumTotBan,LumMoyBan,ETBan,SommeBan
      real RBan,BruitTot,BruitMoy,ETBruit,SommeBruit
      integer Lum(1024,1024),Diagr(256),bruit(1024,1024)
      integer i,j,n,m,x,y,x0,y0,xmax,ymax,LumMax,LumMin

      real k,l,a1,a2,b1,b2,p,p1max,p2max,hwhm1,hwhm2,sig1,sig2,cte
      integer s1,s2

      real lat0,lon0,offset,gain
      integer valmax
      character*40 outfile


      k=0.
      l=0.
      xmax=0.
      ymax=0.
      LumMax=0
      LumMin=255
      Somme=0.
      ET=0.
      x0=31
      y0=24
      RBan=30.


c ============================================
c lecture des données pour le bruit
      open(unit=1,file='PLum.in',status='unknown')
      read(1,*)
      read(1,*) x0,y0
      read(1,*) 
      read(1,*) s1, b1
      read(1,*) s2, b2
      read(1,*) p1max,hwhm1
      read(1,*) p2max,hwhm2
      read(1,*) cte
      close(unit=1)

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
c ==========================================
c transforme l'image en matrice
	open (unit=2,file='PLum.pgm',status='unknown')
	do i=1,6
	 read(2,*) 
	enddo
	read(2,*) n,m	
	read(2,*) ((Lum(i,j),i=1,n),j=m,1,-1)
	close(unit=2)

      do i=1,256
	Diagr(i)=0.
      enddo

      do x=1,n
        do y=1,m
c calcul la lum totale et la moyenne
          LumTot=Lum(x,y)+LumTot
          k=k+1.
c calcul la lum tot et moy en banlieue
          if(sqrt(real((x-x0)**2.+(y-y0)**2.)).gt.(RBan))then
	    LumTotBan=Lum(x,y)+LumTotBan
	    l=l+1
	  endif  
C diagramme à bandes
          Diagr(Lum(x,y)+1)=Diagr(Lum(x,y)+1)+1.
c max et min
	  if(Lum(x,y).gt.LumMax)LumMax=(Lum(x,y))
	  xmax=x
	  ymax=y
	  if(Lum(x,y).lt.LumMin)LumMin=(Lum(x,y))
c extraction du bruit
	  call GaussLisse(x,x0,y,y0,a1,a2,b1,b2,
     + e1,e2,p1max,p2max,sig1,sig2,cte,p)
	  bruit(x,y)=Lum(x,y)-p
        enddo
      enddo

      LumMoy=LumTot/k
      LumMoyBan=LumTotBan/l

      k=0

      do x=1,n
        do y=1,m
c calcul la lum totale et la moyenne de la matrice bruit
          BruitTot=bruit(x,y)+BruitTot
          k=k+1.
        enddo
      enddo
c moyenne du bruit
      BruitMoy=BruitTot/k

c écarts-type
      k=0.
      do x=1,n
	do y=1,m
	  Somme=((Lum(x,y)-LumMoy)**2.+Somme)
	  SommeBruit=((bruit(x,y)-BruitMoy)**2.+SommeBruit)
	  k=k+1.
          if(sqrt(real((x-x0)**2.+(y-y0)**2.)).gt.(30.))then
	    SommeBan=((Lum(x,y)-LumMoyBan)**2.+SommeBan)
	    l=l+1
	  endif
	enddo
      enddo
      ET=sqrt(Somme/k)
      ETBruit=sqrt(SommeBruit/k)
      ETBan=sqrt(SommeBan/l)

c appel de la routine d'écriture de l'image pour le bruit
      outfile='Bruit.pgm'
      lat0=45.
      lon0=72.
      pixsize=1.
      valmax=255.
      gain=64.
      offset=0.
         call write2d (outfile,real(bruit),lat0,lon0,pixsize,n-1,m-1,
     +   valmax,gain,offset)

c écriture de la caractérisation
      open(unit=3,file='PLum.car',status='unknown')
      write(3,*) 'caractérisation de la ville + banlieu'
      write(3,*) ''
      write(3,*) 'Lum Totale :'
      write(3,*) LumTot
      write(3,*) 'Lum Moyenne :'
      write(3,*) LumMoy
      write(3,*) 'Écart-Type :'  
      write(3,*) ET  
      write(3,*) 'Max:   x=   y='
      write(3,*) LumMax,xmax,ymax
      write(3,*) 'Min :'
      write(3,*) LumMin
      write(3,*) ''


      write(3,*) 'caractérisation de la banlieu'
      write(3,*) ''
      write(3,*) 'Lum Totale :'
      write(3,*) LumTotBan
      write(3,*) 'Lum Moyenne :'
      write(3,*) LumMoyBan
      write(3,*) 'Écart-Type :'
      write(3,*) ETBan    
      write(3,*) ''

      write(3,*) 'caractérisation de la matrice bruit'
      write(3,*) ''
      write(3,*) 'Bruit Total :'
      write(3,*) BruitTot
      write(3,*) 'Bruit Moyen :'
      write(3,*) BruitMoy
      write(3,*) 'Écart-Type :'
      write(3,*) ETBruit
      write(3,*) ''

      write(3,*)'diagramme des fréquences de la luminosité générale'
      write(3,*)''

c diagramme à bandes   
      write(3,*) 'luminosité, effectif'
      do i=1,256
	write(3,*)i-1,Diagr(i)
      enddo
      close(unit=3)

      stop
      end