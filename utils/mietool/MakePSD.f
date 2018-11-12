C    PROGRAMME DE CONSTRUCTION DE LA PSD POUR UNE POPULATION DE 
C    PARTICULES SPHERIQUES LE FICHIER DE SORTIE DE CE PROGRAMME 
C    SERT D ENTREE AU PROGRAMME MieS.f
C    LE FICHIER DE SORTIE EST EXPRIME EN TERME DE PARAMETRE DE TAILLE
C    X=2 PI / LAMBDA
C
C -----------------
C     Identification des variables (Martin Aube 1998)
C
C        ITYPE = type de particule
C        PSD(*) = distribution de taille pour chaque dX (dN/dX)
C        M = indice de refraction complexe
C        IPSD = type de distribution de taille
C        RINC = increment de la distribution de taille (dX) X=size 
C               param
C               N.B. le size param = 2 pi r / lambda
C        ISIZE = nombre d'intervalles dX de la distribution de taille
C                par defaut  5001
C        RLOW = Valeur minimale du size parameter
C        RUP = Valeur maximale du size parameter
C        IUP = entier representant le bin maximal 
C        ILOW = entier representant le bin minimal
C        RINCINV = 1/dX
C        XM = taille moyenne de la distribution log-normale
C        SG = sigma de la distribution log-normale
C        ANGLEL = angle minimal pour la fonction de phase
C        ANGLEU = angle maximal pour la fonction de phase
C        DELTA = increment de l'angle pour la fonction de phase
C        ANS = reponse a une question
C        NANG = nombre d'angles pour la fonction de phase
C        RIMP = RELATIVE IMPORTANCE D UN MODE ICI COMME ON NE TRAITE
C               QUE LE CAS UNI-MODAL, RIMP=1.
C        INFILE = NOM DU FICHIER D ENTREE
C        OUTFILE = NOM DU FICHIER DE SORTIE
C        NBINS = nombre de classe de taille de la distribution par bins
C        BEGINBIN = TAILLE DU DEBUT DE CHAQUE BINS DE LA DISTRIBUTION
C        LARBIN = LARGEUR DE CHAQUE BIN (EN PARAMETRE DE TAILLE)
C        DBINMIN = PLUS PETITE DECIMALE DE TOUTES LES LARBIN
C        BINDENSITY = NUMBER DENSITY POUR CHAQUE BIN (SUPPOSEE
C                     CONSTANTE A L INTERIEUR D UN MEME BIN)
C
C --------------------C
C   Programme principal
C
       PROGRAM makepsd
C
C ----------
C
C   DECLARATION DES VARIABLES
C
       REAL*8 PSD(5002),RN,RIMP,PI
       REAL*8 BEGINBIN(5002),LARBIN(5002)
       REAL*8 DBINMIN, BINDENSITY(5002), RINC, RUP, RLOW,A,R2PI
       REAL*8 GEORAD, SIGLN, DS, DU, DA, DB, DG, DN,X
       CHARACTER*60 INFILE, OUTFILE
       INTEGER INFILELENGTH, NBINS, II, IPSD, interp,outflen,ILOW,IUP
       integer I,ISIZE
C   
C ----------
C
C   INITIALISATION DE VARIABLES
C
C
C     ** mod ** initialization mod. because original code was only 
C     testedon compilers with automatic zero initialization (for PSD
C     this is a real problem because a zero value is assumed in it's
C     definition below.
C
       DO I=1,5002
         PSD(I)=0.
       enddo
       PI=3.14159265359
       RN = 0.
       RINC = 0.
       ILOW=5001
       IUP=0
       R2PI=SQRT(2.*PI)
       ISIZE=5001
       RIMP=1.
       NBINS=0
       DO I=1,5002
         BEGINBIN(I)=0.
         LARBIN(I)=0.
       enddo
       II=0
       DBINMIN=100000000.

C ----------
C
C   FORMAT D AFFICHAGES
C
 1008 FORMAT(1X,'WHICH PARTICLE SIZE DISTRIBUTION ON A SIZE 
     1PARAMETER SCALE:')
 1009 FORMAT(1X,'1) MONODISPERSED          X')
 1010 FORMAT(1X,'2) GATES-GAUDIN-SCHUMANN  X**(-A), A>1')
 1011 FORMAT(1X,'3) LOG-NORMAL             1/(SG*X*SQRT[2*pi])*EXP(-(LOG
     +(X)-LOG(XM))**2/(2*SG**2))')
 1012 FORMAT(1X,'4) GAMMA                  U**(U+1)/U!*(X**U/S**(U+1))*E
     +XP(-U*X/S)')
 1023 FORMAT(1X,'5) MODIFIED GAMMA         X**A*EXP(-B*X**G)')
 1013 FORMAT(1X,'6) ROSIN-RAMMLER          X**(N-1)*EXP(-B*X**N)')
 1015 FORMAT(1X,'7) INDEPENDANT SIZE-BINS')
 1041 FORMAT(1X,'INCREMENT MUST BE .01 OR GREATER')
 1042 FORMAT(1X,'INCREMENT TOO SMALL! END OF PROGRAM')
 1038 FORMAT(1X,'WHICH POLARIZATION STATE ? 1) PARALLEL')
 1035 FORMAT(1X,'                           2) PERPENDICULAR')
 1032 FORMAT(1X,'                        OR 3) RANDOM')
 1030 FORMAT(1X,'WITH RESPECT TO SCATTERING PLANE.',/)
 1069 FORMAT(1X,'ENTER PARTICLE SIZE OR LENGTH PARAMETER',/)
 1075 FORMAT(1X,F6.2,1X,G12.5)
 1076 FORMAT(1X,'QEXT=',G15.8,1X,'QSCA=',G15.8,1X,'QABS=',G15.8)
 1078 FORMAT(1X,'MASS EXTINCTION COEF.=',G15.8,'/LAMBDA/DENSITY')
 1077 FORMAT(1X,'LIDAR RATIO=',G17.8)
 1027 FORMAT(1X,'INPUT INTEGRATION dX [d(SIZE PARAMETER)] OVER SIZE 
     +DIST.',/)
 1128 FORMAT(1X,'SIZE PARAMETER IS DEFINED AS X=2*PI*R/WAVELENGTH',/)
 1057 FORMAT(1X,'LOWER AND UPPER LIMITS OF PSD DIST. IN SIZE PARAMETER 
     +SCALE (BETWEEN X=0 AND X=',F9.2,')',/)
 1022 FORMAT(1X,'TOO MANY SIZES, PLEASE REDUCE LOWER AND/OR UPPER
     +LIMITS',/)
 1021 FORMAT(1X,'INPUT mean SIZE PARAMETER, sigma',/)
 7001 FORMAT(1X,'DO YOU WANT TO SEE PSD ? (0=YES)',/)
 9991 FORMAT(1x, /'Press return when ready to proceed'/)
 7000 FORMAT(4(1X,G9.3,1X,G9.4))
 7002 FORMAT(1X,'NUMBER OF SIZE BINS (max 500)')
 7004 FORMAT(1X,'SIZE PARAMETER OF THE BEGINNING OF THE BIN #',I3,
     +'(X=2*PI*R/WAVELENGTH)?')
 7009 FORMAT(1X,'NUMBER DENSITY/dX FOR THE BEGINNING OF THE BIN #',I3
     +,'( dN/dX )?')
 7007 FORMAT(1X,'MAXIMUM SIZE PARAMETER THE LAST BIN?')
C 7010 FORMAT(1X,'CENTER WAVELENGTH OF THE SENSOR (MICRONS)?')
 1028 FORMAT(1X,'VALUE FOR A ?',/)
 1029 FORMAT(1X,'INPUT S-MOST PROB. RADIUS & U WHERE 1/(U+1)=
     +HALF WIDTH',/)
 1031 FORMAT(1X,'INPUT A,B,G',/)
 1033 FORMAT(1X,'INPUT B=1/MODE OF DISTRIBUTION & N',/)
C
C ----------
C
C   LONGUEUR DU NOM DU FICHIER DE SORTIE
C
 2001  PRINT*,'Root name of the PSD output file ?'
       READ(*,2005,ERR=2001) OUTFILE
       outflen=INDEX(OUTFILE,' ')-1
       OUTFILE=OUTFILE(1:outflen)//'.mie.psd'
 2005  FORMAT(A)
C
C ----------
C
C   OUVERTURE DU FICHIER DE SORTIE 
C
C
       OPEN (UNIT=7,FILE=OUTFILE,STATUS='UNKNOWN')
C
C   FORMAT DU FICHIER DE SORTIE:
C 
C   TYPE DE DISTRIBUTION (IPSD)
C   NOMBRE DE BINS (NBINS)  
C   TAILLE DEBUT DU 1ER BIN (MICRON), POPULATION DU 1ER BIN (RELATIVE) 
C   TAILLE DEBUT DU 2E BIN (MICRON), POPULATION DU 2E BIN (RELATIVE)
C    .
C    .
C    .
C   TAILLE DE LA FIN DU [NBINS]E BIN (MICRON)
C
C   TYPE DE DISTRIBUTION DE TAILLE (PSD) POSSIBLES
C
 2010      WRITE(*,1008)
           WRITE(*,1009)
           WRITE(*,1010)
           WRITE(*,1011)
           WRITE(*,1012)
           WRITE(*,1023)
           WRITE(*,1013)
           WRITE(*,1015)
C
C --------------------
C
C   Choix de la distribution de taille
C
C
           READ(*,*,ERR=2010) IPSD
           IF ((IPSD.GT.7).OR.(IPSD.LT.1)) GOTO 2010
C
C   Type d'interpolation entre les valeurs
C
           if (IPSD.ne.1) then
 2015       print*,'Interpolation type betwen PSD s values ?'
            print*,'0 = Constant value inside the bin'
            print*,'1 = Linear interpolation'
            print*,'9 = Discret mode (non zero values for the begin. of
     +the bin)'            
            read(*,*,ERR=2015) interp 
            if ((interp.NE.1).and.(interp.NE.0).and.(interp.NE.9)) then
              goto 2015
            endif            
           else
            interp=0
           endif 
 10        format(I1)
C
C
           IF (IPSD.EQ.1) THEN
C
C   MONODISPERSED
C
                WRITE(8,1009)
                NBINS=1.
 2030           WRITE(*,1069)
                READ(*,*,ERR=2030) X
                IF (X.LE.0) GOTO 2030
                BEGINBIN(1)=X
                BEGINBIN(2)=X+X/100.
                PSD(1)=1.
                RINC=X/100.
           ELSE
C
C   POLYDISPERSED
C 
                IF (IPSD.NE.7) THEN
 2040           WRITE(*,1027)
                WRITE(*,1128)
                READ(*,*,ERR=2040) RINC
                IF (RINC.LE.0) GOTO 2040
 2050           WRITE(*,1057) FLOAT(ISIZE-1)*RINC
                READ(*,*,ERR=2050) RLOW,RUP
                IF ((RLOW.LT.0).OR.(RUP.LE.0)) GOTO 2050
                IF (RLOW.GE.RUP) GOTO 2050
                IF(RLOW.LT.RINC) RLOW=RINC
                RINCINV=1./RINC
                ILOW1=RLOW*RINCINV+0.5
                IUP1=RUP*RINCINV+0.5
                IF(IUP1*RINC.GT.ISIZE) THEN
                     WRITE(*,1022)
                     GOTO 2050
                ELSE
                     IF(ILOW.GT.ILOW1) ILOW=ILOW1
                     IF(IUP.LT.IUP1) IUP=IUP1
                ENDIF
                ENDIF
                IF (IPSD.EQ.2) THEN
C
C  DISTRIBUTION GATES-GAUDIN-SCHUMANN
C
                     WRITE(*,1010)
 2060                WRITE(*,1028)
                     READ(*,*,ERR=2060) A
                     IF (A.LE.1.) GOTO 2060
                     NBINS=IUP1-ILOW1+1
                     I=ILOW1
                     DO 6200 II=1,NBINS
                          XI=I*RINC
                          BEGINBIN(II)=XI
                          PSD(II)=(I*RINC)**(-A)*(A-1.)*RLOW**
     +                    (A-1.)*RIMP
                     I=I+1
 6200                CONTINUE
                     BEGINBIN(NBINS+1)=XI+RINC
                ELSEIF (IPSD.EQ.3) THEN
C
C   DISTRIBUTION LOG-NORMALE
C
                     WRITE(*,1011)
 2070                WRITE(*,1021)
                     READ(*,*,ERR=2070) XM,SIGLN                           !** mod **

                     IF (XM.LE.0.) then
                       print*,'Please enter a mean size parameter > 0'
                       GOTO 2070
                     elseIF (SIGLN.LE.0.) then
                       print*,'Please enter a sigma > 1'
                       GOTO 2070
                     endif
c                     SG = LOG(SIGLN) !** mod **
                      SG=SIGLN
		    
                     NBINS=IUP1-ILOW1+1
                     I=ILOW1
                     DO 612 II=1,NBINS
                          XI=I*RINC
                          BEGINBIN(II)=XI
                          PSD(II)=(1./SG/XI/R2PI)
     +                    *EXP(-((LOG(XI)-LOG(XM))**2./2./SG**2.))*RIMP
c                          PSD(II)=PSD(II)+(0.434294482/SG/XI/R2PI)
c     +                    *EXP(-((LOG(XI)-LOG(XM))**2./2./SG**2.))*RIMP
                     I=I+1
  612                CONTINUE
                     BEGINBIN(NBINS+1)=XI+RINC
                ELSEIF (IPSD.EQ.4) THEN
C
C   DISTRIBUTION GAMMA
C
                     WRITE(*,1012)
 2080                WRITE(*,1029)
                     READ(*,*,ERR=2080) DS,DU
                     IF ((DU.LE.0.).OR.(DS.LE.0.)) GOTO 2080
                     IU=INT(DU)
                     FU=DU-IU
                     GAMMA=1.-.5748646*FU+.9512363*FU**2.-.6998588*FU**
     +               3.+.4245549*FU**4.-.1010678*FU**5.
                     IF (IU.NE.0) THEN
                          DO 698 II=1,IU
                               GAMMA=(FLOAT(II)+FU)*GAMMA
 698                      CONTINUE
                          NBINS=IUP1-ILOW1+1
                          I=ILOW1
                     ENDIF
 699                 DO 701 II=1,NBINS
                          XI=I*RINC
                          BEGINBIN(II)=XI      
                          PSD(II)=1.D0/GAMMA*DU**(DU+1.)*(XI**
     +                    DU/DS**(DU+1.))*EXP(-DU*XI/DS)*RIMP
                          I=I+1
 701                 CONTINUE
                     BEGINBIN(NBINS+1)=XI+RINC
                ELSEIF (IPSD.EQ.5) THEN
C
C   DISTRIBUTION MODIFIED GAMMA
C
                     WRITE(*,1023)
 2090                WRITE(*,1031)
                     READ(*,*,ERR=2090) DA,DB,DG
                     IF (DB.LE.0.) GOTO 2090
                     DU=(DA+1.)/DG
                     IU=INT(DU)
                     FU=DU-IU
                     GAMMA=1-.5748646*FU+.9512363*FU**2.-.6998588*FU**3.
     +               +.4245549*FU**4.-.1010678*FU**5.
                     IF (IU.NE.0) THEN
                          DO 753 II=1,IU
                          GAMMA=(FLOAT(II)+FU)*GAMMA
 753                      CONTINUE
                     ENDIF
 752                 RN=(1./DG)*DB**(-DU)*GAMMA
                     NBINS=IUP1-ILOW1+1
                     I=ILOW1
                     DO 751 II=1,NBINS
                          XI=I*RINC
                          BEGINBIN(II)=XI
                          PSDARG = XI
                          PSD(II)=PSDARG**DA*EXP(-DB*PSDARG**DG)
     +                    /RN*RIMP!** mod **
                          I=I+1
 751                 CONTINUE
                     BEGINBIN(NBINS+1)=XI+RINC
                ELSEIF (IPSD.EQ.6) THEN
C
C   DISTRIBUTION ROSIN-RAMMLER
C
                     WRITE(*,1013)
 2100                WRITE(*,1033)
                     READ (*,*,ERR=2100) DB,DN
                     IF (DB.LE.0.) GOTO 2100
                     NBINS=IUP1-ILOW1+1
                     I=ILOW1
                     DO 776 II=1,NBINS
                          XI=I*RINC
                          BEGINBIN(II)=XI
                          PSD(II)=XI**DN/XI*EXP(-DB*XI**DN)*DB*
     +                    DN*RIMP
                          I=I+1
 776                 CONTINUE
                     BEGINBIN(NBINS+1)=XI+RINC
                ELSEIF (IPSD.EQ.7) THEN
C
C -----------------
C
C   DISTRIBUTION VARIABLE DEPENDANTE D UN NOMBRE DE GAMMES DE TAILLE
C   DEFINIE PAR L USAGER (DANS LE FICHIER D ENTREE)
C
                RINC=99999
                WRITE(*,1015)
 2110           WRITE(*,7002)
                READ(*,*,ERR=2110) NBINS
                IF (NBINS.LT.1) GOTO 2110
                DO 7005 II=1,NBINS
 2120                WRITE(*,7004) II
                     READ(*,*,ERR=2120) BEGINBIN(II)
                     IF ((BEGINBIN(II).LT.0.).OR.(BEGINBIN(II).LE.
     +               BEGINBIN(II-1))) GOTO 2120
 2130                WRITE(*,7009) II
                     READ(*,*,ERR=2130) BINDENSITY(II)
                     IF (BINDENSITY(II).LT.0.) GOTO 2130
 7005           CONTINUE
 2140           WRITE(*,7007) 
                READ(*,*,ERR=2140) BEGINBIN(NBINS+1)
                IF ((BEGINBIN(NBINS+1).LE.0.).OR.(BEGINBIN(NBINS+1).LE.
     +          BEGINBIN(NBINS))) GOTO 2140
C
C   DEFINITION DE LA DISTRIBUTION DE TAILLE
C
                          DO 7058 II=1,NBINS
                                    PSD(II)=BINDENSITY(II)
 7058                     CONTINUE
                ENDIF
           ENDIF
C
C ------------------------------------------------
C
C   ECRITURE DANS LE FICHIER
C
C
       WRITE(7,*) IPSD,' ** size distribution type **'
       WRITE(7,*) RINC,' ** dX constant increment (99999=variable) **'
       WRITE(7,*) NBINS,' ** Number of bins **'
       WRITE(7,*) interp,' ** Interpolation (0=terrasse, 1=linear, 9=dis
     +cret **'
C       WRITE(7,*) WAVLEN,' ** Wavelength in microns **'
       DO II=1, NBINS
            WRITE(7,*) BEGINBIN(II), PSD(II),' ** X, Psd **'
       enddo
            WRITE(7,*) BEGINBIN(NBINS+1),' ** last X **'
C
C ---------------
C
C   FERMETURE DU FICHIER
C
       CLOSE(UNIT=7)
C
C ---------------
C
C   FIN DU MAIN
C
       STOP
       END
