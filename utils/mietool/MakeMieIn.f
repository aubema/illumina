C   LOGICIEL PERMETTANT DE CREER LE FICHIER D ENTREE DE DONNEES POUR 
C   LE LOGICIEL MieS. 
C
C   Martin Aube 1998
C
C -----------------
C
C     Identification des variables
C
C        ITYPE = type de particule
C        REM = partie reelle de l incide de refraction
C        IMM = partie imaginaire de l indice de refraction
C        ANGLEL = angle minimal pour la fonction de phase
C        ANGLEU = angle maximal pour la fonction de phase
C        DELTA = increment de l'angle pour la fonction de phase
C        ANS = reponse a une question
C
C --------------------
C
C   Programme principal
C
       PROGRAM MakeMieIn
C
C ----------
C
C   DECLARATION DES VARIABLES
C
      REAL*8 REM,IMM,MABS, ANGLEL, ANGLEU, DELTA,REM1,IMM1,Y,X
     1       ,WAVELEN
      CHARACTER*60 , INFILE, PSDFILE
      INTEGER INFILELENGTH, ANS, ITYPE, IMPRIMER
      COMMON /DATA/ ANGLEL,ANGLEU,DELTA
C
C ----------
C
C   NOM DU FICHIER D ENTREE ET OUVERTURE
C
 100   PRINT*,'Root name of the imies input file ?'
       READ(*,110 ,ERR=100 ) INFILE
       INFILELENGTH=INDEX(INFILE,' ')-1
       INFILE=INFILE(1:INFILELENGTH)//'.mie' 
       OPEN(UNIT=1,FILE=INFILE,STATUS='UNKNOWN',ERR=100)
 110   FORMAT(A)
C
C ----------
C
C   INITIALISATION DE VARIABLES
C
       REM = 0.
       IMM = 0.
       ANGLEL = 0.
       ANGLEU = 0.
       DELTA = 0.
       IMPRIMER=0
C
C
C ----------
C
C   DEBUT
C
C ----------
C
C   FORMAT D AFFICHAGES
C
 2000 FORMAT(1X,'WHICH PARTICLE TYPE:')
 2010 FORMAT(1X,'1)  SPHERE                    EXACT')
 2020 FORMAT(1X,'2)  COATED SPHERE             EXACT')
 2021 FORMAT(1X,'ENTER SENSOR WAVELENGTH (IN MICRONS):')
 2120 FORMAT(1X,'INDEX OF REFRACTION m & k (of m - ik) ->',/)
 2130 FORMAT(1X,'INCREMENT MUST BE .01 OR GREATER')
 2140 FORMAT(1X,'INCREMENT TOO SMALL! END OF PROGRAM')
 2150 FORMAT(1X,'ENTER LOWEST,HIGHEST AND INCREMENT OF ANGLES IN PHASE
     $ FUNCTION') 
 2160 FORMAT(1X,'WHICH POLARIZATION STATE ? 1) PARALLEL')
 2170 FORMAT(1X,'                           2) PERPENDICULAR')
 2180 FORMAT(1X,'                        OR 3) RANDOM')
 2190 FORMAT(1X,'WITH RESPECT TO SCATTERING PLANE.',/)
 2200  FORMAT(I3,1X,' ** Aerosol shape**')
 2210  FORMAT(F10.5,1X,E10.5,' ** Refractive index (a-ib) **')
 2220  FORMAT(F6.2,1X,F6.2,1X,F6.2,' ** i angle , f angle , delta **')
 2230  FORMAT(I2,1X,' ** Polarisation state (1=//,2=_|_,3=random) **')
 2240 FORMAT(1X,'INDEX OF REFRACTION FOR COATING m1 & k1 ?',/)
 2250 FORMAT(1X,'INDEX OF REFRACTION FOR CORE m2 & k2 ?',/)
 2260 FORMAT(1X,'(CORE RADIUS)/(COATING RADIUS)= ? (0=CORE RADIUS IS
     $ CONSTANT)',/)
 2270 FORMAT(1X,'INPUT SIZE PARAMETER OF CORE RADIUS',/)
 2280 FORMAT(1X,'*** WARNING: CALCULATION MAY BE INACCURATE DUE TO
     $ BESSEL FCNS ***',/)
C
C   OPTIONS D AFFICHAGE A L ECRAN
C
 2295  PRINT*,'Verbose imies mode? (1 = yes)'
       READ(*,*,ERR=2295) IMPRIMER
       WRITE(1,*) IMPRIMER,' ** 1 = VERBOSE **'

C
C
C ----------
C
C   CHOIX DU TYPE DE FORME DES AEROSOL (ITYPE)
C
 2300  WRITE(*,2000)
       WRITE(*,2010)
       WRITE(*,2020)
       READ(*,*,ERR=2300) ITYPE
C
C------------
C WAVELENGTH
       WRITE(*,2021)
       READ(*,*) WAVELEN
       WRITE(1,*) WAVELEN,' ** SENSOR WAVELENGTH **'
       IF ((ITYPE.GT.2).OR.(ITYPE.LT.1)) GOTO 2300
       WRITE(1,2200) ITYPE
C
C   PROPRIETES OPTIQUES
C
 2310  WRITE(*,2150)
       WRITE(*,2130)
       READ(*,*) ANGLEL,ANGLEU,DELTA
       IF (DELTA.LT.0.01) GOTO 2310
       IF (DELTA.GT.(ANGLEU-ANGLEL)) GOTO 2310
       IF (ANGLEL.GE.ANGLEU) GOTO 2310
       IF ((ANGLEL.LT.0).OR.(ANGLEL.GT.360)) GOTO 2310
       IF ((ANGLEU.LT.0).OR.(ANGLEU.GT.360)) GOTO 2310
       WRITE(1,2220) ANGLEL,ANGLEU,DELTA
 2320  WRITE(*,2160)
       WRITE(*,2170)
       WRITE(*,2180)
       WRITE(*,2190)
       READ(*,*,ERR=2320) ANS
       IF ((ANS.GT.3).OR.(ANS.LT.1)) GOTO 2320
       WRITE(1,2230) ANS
C
C   SPHERE UNIFORME
C
       IF (ITYPE.EQ.1) THEN
 2330       WRITE(*,2120)
            READ(*,*,ERR=2330) REM,IMM
            IF ((REM.LT.1).OR.(IMM.LT.0)) GOTO 2330
            WRITE(1,2210) REM, IMM
C
C   HOMOGENEOUS COATED SPHERE
C
       ELSEIF (ITYPE.EQ.2) THEN
 2340       WRITE(*,2240)
            READ(*,*,ERR=2340) REM1,IMM1
            IF ((REM1.LT.1).OR.(IMM1.LT.0)) GOTO 2340
            WRITE(1,*) REM1,IMM1,' ** REM1,IMM1 **'
 2350       WRITE(*,2250)
            READ(*,*) REM,IMM
            IF ((REM.LT.1).OR.(IMM.LT.0)) GOTO 2350
            WRITE(1,*) REM,IMM,' ** REM, IMM **'
 2360       WRITE(*,2260)
            READ(*,*) Y
            IF(Y.GE.1.AND.Y.LT.0.) GOTO 2360
            IF (Y.EQ.0) THEN
                 WRITE(*,2270)
                 READ(*,*) X
                 WRITE(1,*) Y,' ** Y **'
                 WRITE(1,*) X,' ** X **'
                 IF(X*IMM.GT.30..OR.X*IMM1.GT.30.) WRITE(*,2280)
            ENDIF
            WRITE(1,*) Y,' ** Y **'
       ENDIF
       CLOSE(UNIT=1)
       STOP
       END
