c   programme permettant de combiner deux psd
c   telles que generees par MakePSD.f
c   Dans un premier temps, les psd sont normalisees (int de tous les 
c   bins =1.  Il seront ensuite additionnes en proportion definies par
c   l utilisateur
c
       program cmbpsd
c
c   declaration des variables
c
       integer nbns1,nbns2,nbns,l,m,n,npsd,psdf1l,psdf2l,psdfol,inttyp
       real*8 x(5005), x1(5005),x2(5005),psd1(5005),psd2(5005)
       real*8 psd(5005), norm1, norm2,poids1, poids2
       real*8 xmin,xmax,dx,dx1,dx2
       character*60 psdf1,psdf2,psdoutf, bidon
c
c   entree
c
          print*,'First psd file root name ?' 
          read*,psdf1
 600      print*,'First psd weigth (0-1) ?'
          read*,poids1
          if (poids1.gt.1.) goto 600
          print*,'Second psd file root name ?'
          read*,psdf2
          poids2=1.-poids1
          print*,'Second psd weigth =',poids2
          print*,'Output psd file root name ?'
          read*,psdoutf        
c
c   lecture des fichiers psd
c
       psdf1l=INDEX(psdf1,' ')-1
       psdf2l=INDEX(psdf2,' ')-1
       psdfol=INDEX(psdoutf,' ')-1      
       psdf1=psdf1(1:psdf1l)//'.mie.psd'
       psdf2=psdf2(1:psdf2l)//'.mie.psd'
       psdoutf=psdoutf(1:psdfol)//'.mie.psd'
       xmin=1000000.
       xmax=0.
       nbns=0
          open(unit=1,file=psdf1,status='old')
	  print*,'Reading first psd file...'
             read(1,*) bidon
             read(1,*) dx1
             read(1,*) nbns1
             if (nbns1.gt.nbns) nbns=nbns1
             read(1,*) inttyp
             do 200 n=1,nbns1
             read(1,*) x1(n), psd1(n)
 200         continue 
             read(1,*) x1(nbns1+1)
             if (x1(1).lt.xmin) xmin=x1(1)
             if (x1(nbns1+1).gt.xmax) xmax=x1(nbns1+1)
          close(unit=1)
          open(unit=2,file=psdf2,status='old')
	  print*,'Reading second psd file...'
             read(2,*) bidon
             read(2,*) dx2
             read(2,*) nbns2
             if (nbns2.gt.nbns) nbns=nbns2
             read(2,*) bidon
             do 300 n=1,nbns2
             read(2,*) x2(n), psd2(n)
 300         continue 
             read(2,*) x2(nbns2+1)
             if (x2(1).lt.xmin) xmin=x2(1)
             if (x2(nbns2+1).gt.xmax) xmax=x2(nbns2+1)
          close(unit=2)
        nbns=10*nbns
	if (dx1.ne.dx2) then
	   print*,'Error: PSD scales do not match.'
	   stop
	else 
	   dx=dx1
	endif
	if (nbns1.ne.nbns2) then
	   print*,'Error: PSD lengths do not match.'
	   stop
	else 
	   nbns=nbns1
	endif	
	if (x1(1).ne.x2(1)) then
	   print*,'Error: PSD origins do not match.'
	   stop
	endif	
	
c        dx=(xmax-xmin)/(DBLE(nbns))
c
c   integration et normalisation des psd
c
        print*,'Normalisation des psd...'
        norm1=0.
        do 800 n=1,nbns1
           norm1=psd1(n)*(x1(n+1)-x1(n))+norm1   
 800    continue
        do 900 n=1,nbns1
           psd1(n)=poids1*psd1(n)/norm1
 900    continue
        norm2=0.
        do 1000 n=1,nbns2
           norm2=psd2(n)*(x2(n+1)-x2(n))+norm2   
 1000   continue
        do 1100 n=1,nbns2
           psd2(n)=poids2*psd2(n)/norm2
 1100   continue
          print*,'Constantes de normalisation=', norm1, norm2
c
c   Creation de la psd combinee
c   
        print*,'Combinaison des psd...'
        x(1)=xmin
           do 1515 n=1,nbns
                 psd(n)=psd1(n)+psd2(n)
 1515      continue      
c
c Normalisation de la psd combinee
c
        norm1=0.
        do 1800 n=1,nbns
           norm1=psd(n)*(x1(n+1)-x1(n))+norm1   
 1800    continue
        do 1900 n=1,nbns1
           psd(n)=psd(n)/norm1
 1900    continue
         print*,'Constante de normalisation finale=',norm1
c
c   ecriture du fichier de psd combinee
c
        open(unit=3,file=psdoutf,status='unknown')
	   print*,'Writing output psd...'
           write(3,*) '7  ** type de distribution de taille **'
           write(3,*) dx,'  ** Increment dX constant (99999=variable)
     +  **'
           write(3,*) nbns,' ** nombre de bins de la distribution **'
	  write(3,*) inttyp,' **interpolation type'
           do 2000 n=1,nbns
              write(3,*) x1(n),psd(n), ' ** X, PSD **'
 2000      continue
           write(3,*) x1(nbns)+dx,' ** X final **'
        close(unit=3)
        end
