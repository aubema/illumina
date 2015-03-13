C   UTILISATION INTERACTIVE DU CODE MIES
C
C   PROGRAMME PRINCIPAL
C
      PROGRAM intermie
C
C   DECLARATION DES VARIABLES
C
      INTEGER BINTERAC
      REAL*8 BNBINTER
      CHARACTER*40 BPSDFILE
C
C   INITIALISATION DES VARIABLES
C
      BNBINTER=5000.
      BINTERAC=1
C
C   APPEL DE LA ROUTINE MIES
C
      CALL MIES(BPSDFILE,BNBINTER,BINTERAC)
      STOP
      END
C
C ----------------------------------------------------------------------
C
C   ROUTINE MIES
C
C   CODE MIE UTILISANT UN INPUT PSD PROVENANT DE MakePSD.f
C   ET LES AUTRES DONNEES MORPHOLOGIQUES ET OPTIQUES PROVENANT
C   DU LOGICIEL MakePSD.f
C
C   Martin Aube 1998 
C
C**********************************************************************
C PHASE FUNCTION AND EFFICIENCY FACTORS FOR MOST PARTICLE SHAPES AND
C PARTICLE SIZE DISTRIBUTIONS
C BY BLAIR EVANS/DREV
C JANUARY 1988
C
C   Programme restructure et modifie par Martin Aube aout 1998
C   
C   -Le cas nonodispersed ne fonctionnait pas dans le code original
C    j ai generalise le programme ce qui a permis de traiter le 
C    monodispersed comme un cas particulier.  En fait maintenant toutes
C    les PSD sont traitees par le meme algorithme fonctionnel.  
C   -Les sorties sont uniformisees et un peu plus documentees.
C   -Ajout de commentaires et identification des variables.
C   -Refonte de l'algorithme pour passer du programme spaguetti a
C    un programme un peu plus structure.
C   -Mode d entree des donnees par fichier batch deux applications
C    complementaires (MakeMieIn.f et MakePSD.f) permettent d editer les 
C    fichiers d entree.
C   -Ajout d'une PSD variable ou les limites des size-bins sont choisies
C    par l usager de meme que la valeur de la PSD pour chaque bin.
C   -J ai enleve l option de multi-mode.  Cette option sera integree
C    par la conception d'un logiciel de service qui combinera les 
C    fichiers de sortie (extension .mie.out) pour en fabriquer un 
C    fichier de parametres optiques integres.
C   -Seul les sous-routines n ont pas ete restructurees
C   -Plusieurs test de comparaison des sorties de mie.f et de MieS.f 
C    ont ete realises avec succes.
C   -Autre validation, j ai verifie si la fonction de phase est 
C    normalisee sur toute la sphere (int ff domeg = 4 pi).
C   -Dans la formulation de la PSD lognormale (programme MakePSD.f)
C    la variable SIGLN est entree en microns mais etait consideree 
C    comme si elle etait en unites de parametre de taille (X).  Une 
C    ligne de conversion micron -> X est ajoutee: SIGLN=2*PI*SIGLN/
C    WAVLEN
C   
C   
C   Mods (Norm O'Neill, 06/07/93):
C   1. replace RANDOMF & RANDOMC % RANDOMCC BY RANDF & RANDC % RANDCC
C   since this Fortran sees only 6 characters.
C   2. replace subroutine name MIE by MIEQ to avoid duplication with
C   program name
C   3. Replace COATCYL & COATCYLPHASE by COCYL & COCYLPHASE
C   since this Fortran sees only 6 characters (it will see COCYL
C   and COCYLP).
C   4. replace C:\IPHASE\IPH.DAT by IPH.DAT
C   5. replace READ(5 & WRITE(6 by READ(* and WRITE(*
C   6. add initialization of all mainline variables (known that vector
C   PSD was a problem) to avoid intialization problems
C   7. In MIEPHASE replace ETA1(NL2)=DCMPLX(0.D0,0.D0) by
C   ETA1(NL2 + 1)=DCMPLX(0.D0,0.D0) and ETA2(NL3) = 0.
C   BY ETA2(NL3 + 1) = 0. for reasons which are obvious in the DO
C   loops
C   which follow these initializations (this omission before being
C   fixed caused all kinds of weird problems).
C   8. In SPHEREPHASE added clearer WRITE statements
C   9. In SPHEREPHASE add proper optical weighting (QSCA) to PSD
C   integration of
C   phase function. Also added a WRITE to indicate where processing
C   was
C   in terms of size parameter integration (26 April, 1998). A factor
C   of 4*pi
C   had to be added to the phase function expression to yield
C   acceptable validation results (see mie.tst).
C   10. The modified log-normal PSD expression was incorrectly
C   expressed. The
C   x (RR) argument of this expression was replaced by x*lamda/2*pi
C   and the the user is now asked for a value of lamda (June 7, 1998)
C
C   Below are some reasonable musings which seem to jive with the
C   results
C
C   LOG
C   For a log-normal PSD a logical definition is (L = lambda, x 
C   = 2*pi*r/L, s = sigma):
C   <Q> = cross section / pi*r2**2 (r2 = rN exp[ln^2 s]; see MEO 
C   report)
C       = int{Q*pi*{L*x/(2*pi)}**2*(1/N)*dN/dx}dx / pi*r2**2
C       = int{Q*pi*{L*x/(2*pi)}**2*(1/N)*dN/dr*L/(2pi)}dx / pi*r2**2
C   but:
C   dN/dr = log(e)/r*dN/dlog(r)
C   = log(e)/r/sqrt(2*pi)/log(s)*
C   exp(-(log(r) - log(rN))**2/(2*log(s)**2)
C   = 1/r/sqrt(2*pi)/ln(s)*
C   exp(-(ln(r) - ln(rN))**2/(2*ln(s)**2)
C   = [(2pi)/L]/*1/x/sqrt(2*pi)/ln(s)*
C   exp(-(ln(x) - ln(xN))**2/(2*ln(s)**2) where xN = 2*pi*rN/L
C   Let "PSD" = dN/dr*L/(2pi) = 1/x/sqrt(2*pi)/ln(s)*
C   exp(-(log(x) - ln(xN))**2/(2*ln(s)**2) as per the
C   defining Fortran statement below:
C
C   612 PSD(I)=PSD(I)+(1.D0/SG/XI/R2PI)*EXP(-((LOG(XI)-LOG(XM))**2/
C      +2/SG**2))
C   
C   A little algebra then yields:
C   <Q> = 1/x2**2/sqrt(2*pi)/ln(s)*
C   int{Q*x*exp(-(log(x) - ln(xN))**2/(2*ln(s)**2)}dx
C   a result which is only dependent on size parameter 
C   (x2 = 2*pi*r2/L)
C
C   LOG
C
C   If XXX = EXT, SCA or ABS then one can infer from the code that:
C   SUMN = int{dN/dx}dx = N (total number of all particles) 
C   SUMNS = int{dN/dx*x^2}dx, 
C   so that, area of all part. = A = int{dN/dr*pi*r**2}dr
C   = [lambda^2/(4*pi)] * SUMNS
C   QXXX1 = int{Q*x^2*dN/dx}dx / (4*pi*N)
C   QXXX2 = int{Q*x^2*dN/dx}dx / SUMNS = 4*pi*N*QXXX1 / SUMNS
C
C   then the efficiency factor of <Q>XXX must = QXXX2
C   and the cross sect. (XXX) = <Q>XXX*(A/N) (<Q> * avg. area of 1
C   part.) 
C   = QXXX2*( [lambda^2/(4*pi)] * SUMNS / N )
C   = (4*pi*N*QXXX1 / SUMNS) ( [lambda^2/(4*pi)] * SUMNS / N )
C   = QXXX1*lambda^2 with the units of lambda^2
C
C   Summary:
C   <Q>XXX = QXXX2
C   cross section (XXX) = QXXX1*lambda^2
C   QXXX2 = 4*pi*SUMN*QXXX1 / SUMNS
C
C   Phase function considerations
C
C   Below, in SPHEREPHASE, it is noted that the integrated phase
C   function is: 
C   (4*pi) * int { RI * QSCA * PSD * x**2 } dx / int { QSCA * PSD 
C   * x**2 } dx
C   but clearly PSD(i) must = f( [x(i)*L / 2*pi] ) for uniformly
C   incrementing 
C   values of x(i) (i.e. the quadrature expressions over size
C   distribution).
C   This presents no problem for the log-normal case as long as xN is
C   properly defined (as shown above), but is less simple for the
C   modified gamma.
C
C   MGD
C   modified gamma distribution:
C   From CMGD section below (near statement 751): 
C   PSD = x**DA*EXP(-DB*x**DG)/RN
C   But the proper expression is PSD = r**DA*EXP(-DB*r**DG)/RN
C   = (x*L/2*pi)**DA*EXP(-DB*(x*L/2*pi)**DG)/RN
C   which clearly cannot be separated out into a multiplicative factor
C   times thef(x) expression of statement 751
C   Therefore the original PSD expression cannot be correct for any
C   lambda and x (RR) must be replaced by x*L/2*pi.
C   
C
C   Tests: see mie.tst in Documentation (sub folder to this folder)
C
C     References:
C
C     1. Deirmendjian, D., (1969), Electromagnetic Scattering on
C     Spherical Poly-Dispersions, American Elsevier, New York.
C
C -----------------
C     Identification des variables (Martin Aube 1998)
C
C        ITYPE = type de particule
C        PSD(*) = distribution de taille pour chaque dX (dN/dX)
C        M = indice de refraction complexe
C        REM = partie reelle de l incide de refraction
C        IMM = partie imaginaire de l indice de refraction
C        MABS = Valeur absolue = (RE^2+IM^2)^.5
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
C        POL = ETAT DE POLARISATION
C        NANG = nombre d'angles pour la fonction de phase
C        INFILE = NOM DU FICHIER D ENTREE
C        OUTFILE = NOM DU FICHIER DE SORTIE
C        NBINS = nombre de classe de taille de la distribution par bins
C        BEGINBIN = TAILLE DU DEBUT DE CHAQUE BINS DE LA DISTRIBUTION
C        BINDENSITY = NUMBER DENSITY POUR CHAQUE BIN (SUPPOSEE
C                     CONSTANTE A L INTERIEUR D UN MEME BIN)
C        SUMN = NOMBRE DE PARTICULES PAR UNITE DE VOLUME = N
C        SUMNS= INTEGRALE DE dN/dX*r**2 dX
C        PROD1 = X**2 * dN/dX
C
C --------------------
C       
C     Subroutines Called:
C   
C     SPHEREPHASE PSD quadrature of optical parameters for spheres
C     MIEPHASE    phase function and optical parameters for
C                 homogeneous sphere (i.e. the classical Mie case)
C     RANDCC      Phase function for randomn orientation of coated
C                 perfectly reflecting infinitie cylinders
C     COAT        Mie calculations (including phase function) for
C                 coated spheres
C
C     PROGRAM (FROM RUCK 1970 WITH CORRECTIONS)
C
C --------------------
C
C   FORMAT DU FICHIER DE SORTIE (*.out)
C
C   En premier lieu, il y a deux colonnes de chiffres representant
C   la fonction de phase (angle de diffusion,fonction de phase).  La
C   fonction de phase est definie de telle sorte que si elle est 
C   divisee par 4 pi, on trouve la probabilite de diffusion de la 
C   radiation incidente par unite d angle solide.  La fonction de 
C   phase est definie de telle sorte, que l integrale de la fonction 
C   sur tous les angles solides divise par 4 pi est normalise (=1).
C   
C --------------------
C
C   Programme principal
C
      SUBROUTINE MIES(FILE,NBINTER,INTERAC)
C
C ----------
C
C   DECLARATION DES VARIABLES
C
      REAL*4 RI(0:18000),SUMISS(0:18000),SUMIDL(0:18000)
      REAL*8 QSCA1(5002),QABS(5002),QEXT2,QSCA2,X1,QALPHAD
      REAL*8 PSD(5002),REM,IMM,REM1,IMM1,MABS,RC,RIC,RT,RIT,RN,RIN
      REAL*8 DUMMYR8(10500), BEGINBIN(5002),WAVLEN
      REAL*8 PI, BINDENSITY(5002),DX,DY,NBINTER
      COMPLEX*8 Z,M1,DUMC8(3003)
      COMPLEX*16 M,CM,DUMC16(17500)
      CHARACTER*40 , INFILE, OUTFILE, PSDFILE, FILE
      INTEGER INFILELENGTH, NBINS, II, III, IMPRIMER, PSDFILELENGTH
      INTEGER bncount(5002), bncmax, interp, INTERAC, POL
      COMMON /PHASE/ RI
      COMMON /DATA/ ANGLEL,ANGLEU,DELTA,NANG,ILOW,IUP,RINC,ANGLEU1
      COMMON /INIT/ Z,Q,EP,EDP,PHI
      COMMON /POLAR/ IPARA,IPERP
      COMMON /GLOBALR4/ SUMISS
      COMMON /GLOBALR8/ PSD,QSCA1,QABS,DUMMYR8,WAVLEN
      COMMON /GLOBALC8/ DUMC8
      COMMON /GLOBALC16/ DUMC16
C
C ----------
C
C   INITIALISATION DE VARIABLES
C
      DO I = 1,5002
            PSD(I) = 0.
            QSCA1(I) = 0.
            QABS(I) = 0.
      END DO
      DO I = 0,18000
            RI(I) = 0.
            SUMISS(I) = 0.
            SUMIDL(I) = 0.
      END DO
      QEXT2 = 0.
      QSCA2 = 0.
      X1 = 0.
      QALPHAD = 0.
      REM = 0.
      IMM = 0.
      REM1 = 0.
      IMM1 = 0.
      MABS = 0.
      RC = 0.
      RIC = 0.
      RT = 0.
      RIT = 0.
      RN = 0.
      RIN = 0.
      DX = 0.
      DY = 0.
      ANGLEL = 0.
      ANGLEU = 0.
      DELTA = 0.
      NANG = 0
      RINC = 0.
      ANGLEU1 = 0.
      M1 = CMPLX(0.,0.)
      Z = CMPLX(0.,0.)
      M = DCMPLX(0.,0.)
      CM = DCMPLX(0.,0.)
      EP = 0.
      EDP = 0.
      PHI = 0.
      Q = 0.
      DO I = 1,3003
            DUMC8(I) = CMPLX(0.,0.)
      END DO
      DO I = 1,3003
            DUMC8(I) = CMPLX(0.,0.)
      END DO
      DO I = 1,17500
            DUMC16(I) = DCMPLX(0.,0.)
      END DO
      IPARA=0
      IPERP=0
      IFL=0
      ILOW=5001
      IUP=0
      PI=3.141592653
      R2PI=SQRT(2.*PI)
      ISIZE=5001
      IMODE=0
      IJ=0
      RIMP=1.
      NBINS=0
      DO I = 1,5002
           BEGINBIN(I) = 0.
      END DO
      II=0
      bncmax=5000
      interp=0
C
C ----------
C
C   FORMAT D AFFICHAGES
C
 1002 FORMAT(1X,'WHICH PARTICLE TYPE:')
 3004 FORMAT(1X,'1)  SPHERE                    EXACT')
 3001 FORMAT(1X,'2)  COATED SPHERE             EXACT')
 3006 FORMAT(1X,'3)  ANISOTROPIC COATED SPHERE EXACT')
 3002 FORMAT(1X,'4)  INFINITE CYLINDER         EXACT')
 3003 FORMAT(1X,'5)  FINITE CYLINDER           1ST ORDER VARIATIONAL')
 3005 FORMAT(1X,'6)  COATED INFINITE CYLINDER  EXACT')
 1003 FORMAT(1X,'7)  CUBE                      SEMI-EMPIRICAL')
 1004 FORMAT(1X,'8)  OCTAHEDRA                 SEMI-EMPIRICAL')
 1005 FORMAT(1X,'9)  FLAKE                     SEMI-EMPIRICAL')
 1006 FORMAT(1X,'10) CONVEX-CONCAVE            SEMI-EMPIRICAL')
 1036 FORMAT(1X,'11) OTHER IRREGULAR           SEMI-EMPIRICAL',/)
 1000 FORMAT(1X,'INDEX OF REFRACTION m & k (of m - ik) ->',/)
 1041 FORMAT(1X,'INCREMENT MUST BE .01 OR GREATER')
 1042 FORMAT(1X,'INCREMENT TOO SMALL! END OF PROGRAM')
 1039 FORMAT(1X,'ENTER LOWEST,HIGHEST AND INCREMENT OF ANGLES IN PHASE
     + FUNCTION') 
 1038 FORMAT(1X,'WHICH POLARIZATION STATE ? 1) PARALLEL')
 1035 FORMAT(1X,'                           2) PERPENDICULAR')
 1032 FORMAT(1X,'                        OR 3) RANDOM')
 1030 FORMAT(1X,'WITH RESPECT TO SCATTERING PLANE.',/)
 1075 FORMAT(1X,F6.2,1X,G12.5)
 1076 FORMAT(1X,'QEXT=',G15.8,1X,'QSCA=',G15.8,1X,'QABS=',G15.8)
 1078 FORMAT(1X,'MASS EXTINCTION COEF.=',G15.8,'/LAMBDA/DENSITY')
 1077 FORMAT(1X,'LIDAR RATIO=',G17.8)
 1022 FORMAT(1X,'Too many sizes, please reduce lower and/or upper
     +limits',/)
 7001 FORMAT(1X,'PARTICLE SIZE DISTRIBUTION (PSD)')
 7002 FORMAT(1X,'SIZE PARAM         PSD')
 9991 FORMAT(1x, /'Press return when ready to proceed'/)
 7000 FORMAT(4(1X,G9.3,1X,G9.4))
 1102 FORMAT(1X,'*** warning: calculation may be inaccurate due to
     + bessel fcns ***',/)
 1101 FORMAT(1X,'*** warning: calculation for larger sizes may be
     + inaccurate due to bessel fcns ***')
 5012   FORMAT(1X,'|    ',F9.3,'| ',F13.6,'|')
 5008  FORMAT(1X,'Indice de refraction= ',F5.3,' -i ',E10.5)
C
C --------------
C
C   MESSAGE D INTRO
C
      IF (INTERAC.EQ.1) THEN
c      PRINT*,' '
c      PRINT*,'                  LOGICIEL MieS'
c      PRINT*,' '
c      PRINT*,'    Ce logiciel permet de calculer les parametres '
c      PRINT*,'optiques tels que la fonction de phase et les' 
c      PRINT*,'efficacites' 
c      PRINT*,'(diffusion, absorption, extinction) pour une forme'
c      PRINT*,'d aerosol donnee et pour une distribution de taille'
c      PRINT*,'quelconque.'
c      PRINT*,' '
c      PRINT*,'  - La distribution de taille doit etre construite '
c      PRINT*,'    en executant le logiciel MakePSD .'
c      PRINT*,' '
c      PRINT*,'  - Les autres parametres d entree sont inseres '
c      PRINT*,'    en executant le logiciel MakeMieIn .'
c      PRINT*,' '
c      PRINT*,' '
c      PRINT*,'Martin Aube 1998'
c      PRINT*,' '
c      PRINT*,'Code Mie utilise: Blair Evans/DREV  Janvier 1988'
c      PRINT*,'Modifications du code Mie par Norm O Neill 1993 '
c      PRINT*,'Modifications du code Mie par Martin Aube 1998 '
c      PRINT*,' '
C
C ----------
C
C   NOM DU FICHIER D ENTREE ET OUVERTURE
C
 2001  PRINT*,'Root name for input-output files ?(MAX 40 CHAR):'
       PRINT*,'Extensions .mie (input file), .mie.out (output file) and
     + .mie.psd (psd file)'
       PRINT*,'will be add to the root name.'
       READ(*,2005,ERR=2001) FILE
       ENDIF
       INFILELENGTH=INDEX(FILE,' ')-1
       INFILE=FILE(1:INFILELENGTH)//'.mie'
       PRINT*,'Reading input file:',INFILE(1:INFILELENGTH+4),'.
     +..'
       OPEN(UNIT=1,FILE=INFILE,STATUS='OLD')
 2005  FORMAT(A)
C
C   AFFICHAGE DES MESSAGES
C
       READ(1,*) IMPRIMER
C
C  WAVELENGTH
C
       READ(1,*) WAVLEN
C
C
C
C   LONGUEUR DU NOM DU FICHIER D ENTREE ET NOM DU FICHIER DE SORTIE
C
       OUTFILE=FILE(1:INFILELENGTH)//'.mie.out' 
C
C ----------
C
C   OUVERTURE DU FICHIER DE SORTIE 
C
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'Opening output file:',OUTFILE(1:
     +INFILELENGTH+8),'...'
       ENDIF
       OPEN (UNIT=7,FILE=OUTFILE,STATUS='UNKNOWN')
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'... Termine.'
       ENDIF
C
C ----------
C
C   NOM DU FICHIER PSD ET OUVERTURE
C
       PSDFILE=FILE(1:INFILELENGTH)//'.mie.psd'
       PSDFILELENGTH=INDEX(PSDFILE,' ')-1
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'Ouverture du fichier psd:',PSDFILE(1:PSDFILELENGTH),
     +'...'
       ENDIF
       OPEN(UNIT=2,FILE=PSDFILE,STATUS='OLD')
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'...Done.'
       ENDIF
C
C
C ----------
C
C   LECTURE DU FICHIER contenant la PSD genere par MakePSD
C
C 
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'Reading psd file:',PSDFILE(1:PSDFILELENGTH),'.
     +..'
       ENDIF     
       READ(2,*) IPSD
       READ(2,*) RINC
       READ(2,*) NBINS
       READ(2,*) interp
C       READ(2,*) WAVLEN
       DO II=1, NBINS
            READ(2,*) BEGINBIN(II), BINDENSITY(II)
       END DO
            READ(2,*) BEGINBIN(NBINS+1)
C
C   FRONTIERE SUPERIEURE POUR L INTERPOLATION LINEAIRE
C
            BINDENSITY(NBINS+1)=BINDENSITY(NBINS)
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'...Done.'
       ENDIF
       if (interp.eq.9) then
             bncmax=1
       endif
       IF (INTERAC.EQ.1) THEN
        if (NBINS.ne.1) then
          print*
          print*
 34       PRINT*,'In how many Sub-bins do you wish to split PSD? (MAX.
     + 5000., MIN',NBINS+1,')'
          READ*,NBINTER
          IF (NBINTER.LT.(NBINS))  GOTO 34
          IF (NBINTER.gt.5000.)  GOTO 34
        else 
          NBINTER=2
        endif  
       ENDIF
C
C ----------
C   
C   FERMETURE DU FICHIER PSD
C
       IF (IMPRIMER.EQ.1) THEN
            PRINT*,'Closing psd file.'
       ENDIF     
       CLOSE(UNIT=2)
C ----------
C
C   CHOIX DU TYPE DE FORME DES AEROSOL (ITYPE)
C  
      READ(1,*) ITYPE
      IF (IMPRIMER.EQ.1) THEN
           PRINT*,'Particle shape:'
           IF (ITYPE.EQ.1) THEN
                WRITE(*,3004)
           ELSEIF (ITYPE.EQ.2) THEN
                WRITE(*,3001)
           ELSE
                PRINT*,'Invalid aerosol shape!'
                stop
           ENDIF
      ENDIF
C
C -----------------
C
C   DISTRIBUTION VARIABLE DEPENDANTE D UN NOMBRE DE GAMMES DE TAILLE
C   DEFINIE PAR L USAGER (DANS LE FICHIER D ENTREE)
C
C   TOUTES LES DISTRIBUTION DE TAILLES SONT TRAITEES PAR CET PARTIE
C   DE PROGRAMME.  LES RESULTATS OBTENUS. 
C
C   DETERMINATION DE LINTERVALLE DU PARAMETRE DE TAILLE (DX) 
C
      IF (IMPRIMER.EQ.1) THEN
            PRINT*,'Defining size increment...'
      ENDIF
c   
c   pour que le programme traite toutes les distributions variables
c   en sous ›chantillonnant les bins, il faut mettre en commentaire
c   le IF et le ENDIF qui suit.  Ainsi, la division de l'intervalle 
c   par NBINTER sera toujours effectu›e.
c
           RINC=1.000001*(BEGINBIN(NBINS+1)-BEGINBIN(1))/NBINTER
           RLOW=BEGINBIN(1)
           RUP=BEGINBIN(NBINS+1)-RINC
           IF(RLOW.LT.RINC) RLOW=RINC
           RINCINV=1./RINC
           ILOW1=RLOW*RINCINV+0.5
           IUP1=RUP*RINCINV+0.5
      IF (IMPRIMER.EQ.1) THEN
            PRINT*,'...Done'
      ENDIF
           IF(IUP1*RINC.GT.ISIZE) THEN
                IF (IMPRIMER.EQ.1) THEN
                     WRITE(*,1022)
                ENDIF
                STOP
           ELSE
                IF(ILOW.GT.ILOW1) ILOW=ILOW1
                IF(IUP.LT.IUP1) IUP=IUP1
C
C   DEFINITION DE LA DISTRIBUTION DE TAILLE
C
                IF (IMPRIMER.EQ.1) THEN
                     PRINT*,'Defining size distribution...'
                ENDIF
                DO 7058 I=ILOW1,IUP1
                     XI=I*RINC
                     DO 7058 II=1,NBINS
                          IF ((XI.GE.BEGINBIN(II)).AND.(XI.LT.
     +                    BEGINBIN(II+1))) THEN
c
c   bncount est introduit dans le but de pouvoir traiter le cas 
c   discret.  I.e. lorsqu il y a des particules a certaines valeurs
c   de x mais il n'y en a pas entre ces valeurs.  On fixe alors la 
c   condition (bncount(II).gt.1).  Autrement, (bncount(II).gt.1000)
c   tous les sous bins seront consideres. Nous nous rapprochons
c   alors du cas continu.  Le choix entre 1 et 5000 se fait par 
c   l entree de donnee a l aide de MakePSD.f (variable: interp=9)
c
                               bncount(II)=bncount(II)+1
                               if (bncount(II).gt.bncmax) then
                                    PSD(I)=0.
                               else
                                    PSD(I)=BINDENSITY(II)
                               endif 
                               if (interp.eq.0) then
                                    PSD(I)=BINDENSITY(II)
                               endif
                               if (interp.eq.1) then
                                        PSD(I)=BINDENSITY(II)+(XI-
     +  BEGINBIN(II))*(BINDENSITY(II+1)-BINDENSITY(II))/(BEGINBIN(II+1)-
     +  BEGINBIN(II))
                               endif   
                          ENDIF
 7058           CONTINUE
                IF (IMPRIMER.EQ.1) THEN
                     PRINT*,'|----------------------------|'
                     PRINT*,'|      SIZE DISTRIBUTION     |'
                     PRINT*,'|-------------|--------------|'
                     PRINT*,'| SIZE PARAM. | PSD          |'
                     DO 7059 I=ILOW1,IUP1
                          XI=I*RINC
                          if (PSD(I).ne.0.) then
                               WRITE(*,5012) XI, PSD(I)
                          endif
7059                 CONTINUE
                     PRINT*,'|----------------------------|'
                     IF (IMPRIMER.EQ.1) THEN
                          PRINT*,'...Done.'
                     ENDIF
                ENDIF
C
C   CALCULS DES PROPRIETES OPTIQUES
C
                IF (IMPRIMER.EQ.1) THEN
                     PRINT*,'Optical computations...'
                ENDIF
                READ(1,*) ANGLEL,ANGLEU,DELTA
                IF(DELTA.LT.0.01) THEN
                     IF (IMPRIMER.EQ.1) THEN
                          WRITE(*,1042)
                          PRINT*,'Angular increment to small (<0.01)!'
                     ENDIF
                     STOP
                ENDIF
                READ(1,*) POL
                IF (IMPRIMER.EQ.1) THEN
                     IF(POL.EQ.1) THEN
                          PRINT*,'Parallel prolarization state'
                     ELSEIF(POL.EQ.2) THEN
                          PRINT*,'Perpendicular prolarization state'
                     ELSEIF(POL.EQ.3) THEN
                          PRINT*,'Random polarization state'
                     ELSE
                          PRINT*,'Bad polarisation state !'
                     ENDIF
                ENDIF
                IF(POL.NE.2.) IPARA=1
                IF(POL.NE.1.) IPERP=1
                NANG=(ANGLEU-ANGLEL)/DELTA
                ANGLEU1=NANG*DELTA+ANGLEL
C
C -----------------------------------
C
C   SPHERES UNIFORMES
C
                IF (ITYPE.EQ.1) THEN
                     READ(1,*) REM,IMM
                     IF (IMPRIMER.EQ.1) THEN
                          WRITE(*,5008) REM, IMM
                     ENDIF
                     M=DCMPLX(REM,-IMM)
                     MABS=CDABS(M)
                     CALL SPHEREPHASE(M,IMPRIMER) 
c               ENDIF
C
C -------------------------------------
C
C   HOMOGENEOUS COATED SPHERES 
C
                ELSEIF (ITYPE.EQ.2) THEN
                     READ(1,*) REM1,IMM1
                     M1=CMPLX(REM1,IMM1)
                     READ(1,*) REM,IMM
                     M=DCMPLX(REM,IMM)
                     READ(1,*) Y
                     IF(Y.GE.1.AND.Y.LE.0.) THEN
                          IF (Y.EQ.0) THEN
                               READ(1,*) X
                               WRITE(1,*) X,' ** X **'
                               IF (IMPRIMER.EQ.1) THEN
                                    IF(X*IMM.GT.30..OR.X*IMM1.GT.30.)
     +WRITE(*,1102)
                               ENDIF
                          ENDIF
                     ENDIF
                     CM=DCMPLX(-1.0001D0,0.0D0)
                     IF (IMPRIMER.EQ.1) THEN
                          IF(IUP*RINC*IMM1.GT.30.) WRITE(*,1101)
                     ENDIF
                     CALL COATPHASE(X,Y,M,M1,CM)
                ENDIF
                IF (IMPRIMER.EQ.1) THEN
                     PRINT*,'...Done.'
                ENDIF
                IF (IMPRIMER.EQ.1) THEN
                     PRINT*,'Closing input-output files...'
                ENDIF
C
C ----------
C   
C   FERMETURE DU FICHIER D ENTREE
C
                CLOSE(UNIT=1)
C
C   FERMETURE DU FICHIER DE SORTIE
C
                CLOSE(UNIT=7)
                IF (IMPRIMER.EQ.1) THEN
                     PRINT*,'...Done.'
                ENDIF
           ENDIF
           IF (IMPRIMER.EQ.1) THEN
                PRINT*,' '
                PRINT*,'Computation result in file: ',OUT
     +FILE(1:INFILELENGTH+8)
                PRINT*,'...end of mies !'
           ENDIF
           return
           end
C
C ---------------------------------------------------
C
C PHASE FUNCTION/MIE CODE FOR SPHERICAL PARTICLES
C BY BLAIR EVANS/DREV
C JANUARY 1988
C
C -------------
C     Identification des variables (Martin Aube 1998)
C
C        M = indice de refraction complexe
C        RINC = increment de la distribution de taille (dX) 
c               X=size param
C               N.B. le size param = 2 pi r / lambda
C        IUP = entier representant le bin maximal 
C        ILOW = entier representant le bin minimal
C        ANGLEL = angle minimal pour la fonction de phase
C        ANGLEU = angle maximal pour la fonction de phase
C        DELTA = increment de l'angle pour la fonction de phase
C        NANG = nombre d'angles pour la fonction de phase
C        IDIFF = nombre d'intervalles de taille (bins)
C        X = valeur du param de Mie 
C        
C -------------
      SUBROUTINE MIEPHASE(X,M,QEXT,QSCA,QALPHA)
      REAL*4 RI(0:18000)
      REAL*8 PSI(3500),CHI(3500),ETA2(3500),WAVLEN
      REAL*8 X,P0,SCALE,MABS,C0,YL,T2,QEXT,QSCA
      REAL*8 XL,PI1,PI2,PTMP,PP1,PP2,RAD,TAU,QPH1,QPH2,QPH3,FACT
      REAL*8 DUMMY(15006), somme
      COMPLEX*16 ETA3(3500),ETA1(7000),A(3500),B(3500),M,ZETA
      COMPLEX*16 CDI,Z,QPHASE1,QPHASE2
      COMMON /PHASE/ RI
      COMMON /DATA/ ANGLEL,ANGLEU,DELTA,NANG,ILOW,IUP,RINC,ANGLEU1
      COMMON /POLAR/ IPARA,IPERP
      COMMON /GLOBALR8/ DUMMY,PSI,CHI,ETA2,WAVLEN
      COMMON /GLOBALC16/ ETA3,ETA1,A,B
      SAVE
      DO 111 I=1,3500
  111 CHI(I)=1.E32
      MABS=CDABS(M)
      IF(IPARA.EQ.0.AND.IPERP.EQ.0) GOTO9
      FACT=1./DBLE(IPARA+IPERP)
    9 NL3=1.5*X+10
C CALCULATE PSI(X) FOR ALL ORDERS BY BACKWARD RECURRENCE
C (ABRAMOWITZ & STEGUN P 452)
      IF(X.GT.100) GOTO2
      NL1=2*X+10
      GOTO1
    2 NL1=MAX0(NL3,IDINT(X+4.*X**.533+2))
    1 PSI(NL1+2)=0.D0
      PSI(NL1+1)=1.D0
      KK=0
      DO 20 I=NL1,1,-1
      IF (PSI(I+2).LT.1.E25) GOTO20
      K1=I+2
      KK=KK+1
      PSI(I+1)=PSI(I+1)*1.E-30
      PSI(I+2)=PSI(I+2)*1.E-30
   20 PSI(I)=DBLE(2*I+3)/X*PSI(I+1)-PSI(I+2)
      P0=3.D0/X*PSI(1)-PSI(2)
      SCALE=P0/DSIN(X)
      DO 22 I=1,NL1
   22 PSI(I)=PSI(I)/SCALE
      IF(KK.EQ.0) GOTO999
      DO 23 I=K1+1,NL1
   23 PSI(I)=PSI(I)*1.E-30
  999 P0=DSIN(X)
C CALCULATE ETA1(MX) FOR ALL ORDER BY BACKWARD RECURRENCE
C (KATTAWAR & PLASS APPL. OPT. 6 1377-1382)
      NL2=2.0*MABS*X+10
      Z=M*DCMPLX(X,0D0)
      ETA1(NL2 + 1)=DCMPLX(0.D0,0.D0)! ** mod. ** NL2 + 1 for NL2
      DO 40 I=NL2,1,-1
      CDI=DBLE(I+1)/Z
   40 ETA1(I)=CDI-1.D0/(ETA1(I+1)+CDI)
C CALCULATE ETA2(X) FOR ALL ORDERS BY BACKWARD RECURRENCE
C (KERKER P 67-68)
      ETA2(NL3 + 1)=0.! ** mod. ** NL3 + 1 for NL3
      DO 50 I=NL3,1,-1
   50 ETA2(I)=DBLE(I+1)/X-1./(ETA2(I+1)+DBLE(I+1)/X)
C CALCULATE ETA3(X) AND CHI(X) FOR ALL ORDERS BY FORWARD RECURRENCE
C AND THE MIE SCATTERING COEFFICIENTS
C (KERKER P 67-68)
      C0=DCOS(X)
      CHI(1)=DCOS(X)/X+DSIN(X)
      ETA3(1)=DCMPLX(P0,C0)/DCMPLX(PSI(1),CHI(1))-1.D0/X
      CHI(2)=DBLE(3)/X*CHI(1)-C0
      ETA3(2)=1.D0/(2.D0/X-ETA3(1))-2.D0/X
      ZETA=DCMPLX(PSI(1),CHI(1))
      A(1)=PSI(1)*((ETA1(1)-M*ETA2(1))/(ETA1(1)-M*ETA3(1)))/ZETA
      B(1)=PSI(1)*((ETA2(1)-M*ETA1(1))/(ETA3(1)-M*ETA1(1)))/ZETA
      ZETA=DCMPLX(PSI(2),CHI(2))
      A(2)=PSI(2)*((ETA1(2)-M*ETA2(2))/(ETA1(2)-M*ETA3(2)))/ZETA
      B(2)=PSI(2)*((ETA2(2)-M*ETA1(2))/(ETA3(2)-M*ETA1(2)))/ZETA
      K2=0
      DO 60 I=3,NL3
      ETA3(I)=1.D0/(DBLE(I)/X-ETA3(I-1))-DBLE(I)/X
      IF(K2.EQ.1) GOTO65
      CHI(I)=(2.D0*DBLE(I)-1.D0)/X*CHI(I-1)-CHI(I-2)
      IF(CHI(I).GT.1.E32) K2=1
   65 ZETA=DCMPLX(PSI(I),CHI(I))
      A(I)=PSI(I)*((ETA1(I)-M*ETA2(I))/(ETA1(I)-M*ETA3(I)))/ZETA
      B(I)=PSI(I)*((ETA2(I)-M*ETA1(I))/(ETA3(I)-M*ETA1(I)))/ZETA
   60 CONTINUE
C CALCULATE LEGENDRE FUNCTIONS PI1,PI2,PP1,PP2 AND TAU AT ANGLE THETA
C FOR THETA=0 TO 180 DEGREES
      QEXT=0.D0
      QSCA=0.D0
      I=0
    5 IF (I.LT.X) GOTO6
      IF (I.GE.NL3) GOTO74
      IF (DBLE(A(I)).LE.1.E-17) GOTO74
    6 I=I+1
      QEXT=QEXT+DBLE(2*I+1)*(DBLE(A(I)+B(I)))
      QSCA=QSCA+DBLE(2*I+1)*((CDABS(A(I)))**2+(CDABS(B(I)))**2)
      GOTO5
   74 QEXT=2.D0/X**2*QEXT
      QSCA=2.D0/X**2*QSCA
      QALPHA=4.71238898037/X*QEXT
C
C     ____________________________________________________________
      DO 100 J=0,NANG
      RAD=3.1415926535897939/180.D0*(ANGLEL+DBLE(J)*DELTA)
      XL=DCOS(RAD)
      YL=DSIN(RAD)
      PI1=1.D0
      PI2=3.D0*XL
      PP1=0.D0
      PP2=3.D0
      T2=XL*PI2-YL*YL*3.D0
C AND CALCULATE PHASE FUNCTION (QPHASE1 AND QPHASE2) FOR BOTH
C POLARIZATION STATES AND THE AVERAGE QPH3
      QPHASE1=1.5D0*(A(1)+B(1)*XL)+5.D0/6.D0*(A(2)*PI2+B(2)*T2)
      QPHASE2=1.5D0*(B(1)+A(1)*XL)+5.D0/6.D0*(B(2)*PI2+A(2)*T2)
      I=2
    3 IF (I.LT.X) GOTO4
      IF (I.GE.NL3) GOTO75
      IF (DBLE(A(I)).LE.1.E-17) GOTO75
    4 I=I+1
      PI1=DBLE(2*I-1)/DBLE(I-1)*XL*PI2-DBLE(I)/DBLE(I-1)*PI1
      PP1=DBLE(2*I-1)*PI2+PP1
      TAU=XL*PI1-YL*YL*PP1
      QPHASE1=QPHASE1+DBLE(2*I+1)/DBLE(I*(I+1))*(A(I)*PI1+B(I)*TAU)
      QPHASE2=QPHASE2+DBLE(2*I+1)/DBLE(I*(I+1))*(A(I)*TAU+B(I)*PI1)
      PTMP=PI1
      PI1=PI2
      PI2=PTMP
      PTMP=PP1
      PP1=PP2
      PP2=PTMP
      GOTO3
   75 QPH1=4.D0/X**2*CDABS(QPHASE1)**2
      QPH2=4.D0/X**2*CDABS(QPHASE2)**2
      QPH3=(QPH1*IPARA+QPH2*IPERP)*FACT
  100 RI(J)=QPH3/12.566371/QSCA
      RETURN
      END
C
C ------------------------------------------------------------
C
C INTEGRATE THE PHASE FUNCTION FOR SPHERES OVER THE PARTICLE SIZE
C DISTRIBUTION USING SIMPSON'S INTEGRATION
C BY BLAIR EVANS/DREV 
C JANUARY 1988
C -------------
C     Identification des variables (Martin Aube 1998)
C
C        PSD(*) = distribution de taille pour chaque dX
C        M = indice de refraction complexe
C        RINC = increment de la distribution de taille (dX) X=size param
C               N.B. le size param = 2 pi r / lambda
C        IUP = entier representant le bin maximal 
C        ILOW = entier representant le bin minimal
C        ANGLEL = angle minimal pour la fonction de phase
C        ANGLEU = angle maximal pour la fonction de phase
C        DELTA = increment de l'angle pour la fonction de phase
C        NANG = nombre d'angles pour la fonction de phase
C        IDIFF = nombre d'intervalles de taille (bins)
C        
C        
C -------------
      SUBROUTINE SPHEREPHASE(M,IMPRIMER)
      REAL*4 RI(0:18000),SUMISS(0:18000)
      REAL*8 PSD(5002),DUMR8(20504),QEXT,QSCA,QEXT1,QSCA1,RR,R1,WAVLEN
      COMPLEX*16 M
      INTEGER IMPRIMER
      COMMON /PHASE/ RI
      COMMON /DATA/ ANGLEL,ANGLEU,DELTA,NANG,ILOW,IUP,RINC,ANGLEU1
      COMMON /GLOBALR4/ SUMISS
      COMMON /GLOBALR8/ PSD,DUMR8,WAVLEN
      SAVE
      PI=3.141592653
      IDIFF=IUP-ILOW
      IF(IDIFF/2*2.EQ.IDIFF) GOTO10 !Simpson's "1/3" rule only works for an
C                  even number of intervals (see A&S p. 886) ** mod **
      IUP=IUP+1
C
C     initialisation of PSD quadratures
   10 SUMNS=PSD(ILOW)*ILOW**2*RINC**2
      SUMN=PSD(ILOW)
      SUMA=SUMNS*ILOW*RINC
      RR=ILOW*RINC
      CALL MIEPHASE(RR,M,QEXT,QSCA,QALPHA)
      QEXT1=QEXT*SUMNS
      QSCA1=QSCA*SUMNS
      DO 2001 I=0,NANG
 2001 SUMISS(I)=RI(I)*QSCA*SUMNS!** mod ** (originally RI(I)*SUMNS)
C
C     Actual quadrature.
C     For RI (which by comparison with standard results must = phase function / (4*pi) )
C     this means a Simpson's rule approximation to: int { RI * QSCA * x**2 } dx.
C     Below the quadrature loop the program simply calculates:
C     (4*pi) * int { RI * QSCA * PSD * x**2 } dx / int { QSCA * PSD * x**2 } dx
C     (note that dx is never multiplied but it clearly doesnt matter)
C










      print*,'X		Qext		Qsca'












      DO 2000 I=ILOW+1,IUP,2
        RR=I*RINC
        R1=RR+RINC
        PROD1=PSD(I)*RR**2
        PROD2=PSD(I+1)*R1**2
      	SUMNS=SUMNS+4.*PROD1+2.*PROD2
      	SUMN=SUMN+4.*PSD(I)+2.*PSD(I+1)
      	SUMA=SUMA+4.*PROD1*RR+2.*PROD2*R1
c  cette routine retourne le QEXT pour la valeur de X
      	CALL MIEPHASE(RR,M,QEXT,QSCA,QALPHA)









        print*,RR,QEXT,QSCA










c  convertir integrale ponderee de QEXT a partir de la PSD*X**2
      	QEXT1=QEXT1+4.*QEXT*PROD1
      	QSCA1=QSCA1+4.*QSCA*PROD1
      	DO 2002 II=0,NANG
 2002 	   SUMISS(II)=SUMISS(II)+4.*RI(II)*QSCA*PROD1!** mod ** (originally RI(II)*PROD1)
      	CALL MIEPHASE(R1,M,QEXT,QSCA,QALPHA)
      	QEXT1=QEXT1+2.*QEXT*PROD2
      	QSCA1=QSCA1+2.*QSCA*PROD2
      	DO 2003 II=0,NANG
 2003 	   SUMISS(II)=SUMISS(II)+2.*RI(II)*QSCA*PROD2!** mod ** (originally RI(II)*PROD2)
C
C     	** mod **
      	IF((50*(I/50) - I).EQ.0)THEN
           IF (IMPRIMER.EQ.1) THEN
      	   WRITE(*, 4013) R1
 4013 	   FORMAT(1X, 'finished Mie computations for size parameter = ',
     + f8.2)
           ENDIF
      	END IF
C
 2000 CONTINUE
      SUMNS=SUMNS-PROD2 !subtract out the last even term so that rather
C		than 2 last even terms there is only one (A&S p. 886) ** mod **
      SUMN=SUMN-PSD(IUP)
      SUMA=SUMA+PROD2*R1
      QEXT1=QEXT1-QEXT*PROD2
      QSCA1=QSCA1-QSCA*PROD2
      IF(SUMNS.GT.0.) GOTO14
      WRITE(*,1306)
      STOP
c  diviser la pondration par l integrele du coeff de ponder.
c  i.e. QEXT2= QEXT moyen
   14 QEXT2=QEXT1/SUMNS
      QSCA2=QSCA1/SUMNS
      SUMA=QEXT1*3.*PI/2./SUMA
      QEXT1=.25*QEXT1/SUMN/PI
      QSCA1=.25*QSCA1/SUMN/PI 
      WRITE(7,*) '|--------|------------|'
      WRITE(7,*) '| ANGLE  |  PHASE FCT |'
      WRITE(7,*) '|--------|------------|'
      DO 2004 II=0,NANG
      SUMISS(II)=4.*PI*(SUMISS(II)-RI(II)*QSCA*PROD2)/(QSCA2*SUMNS)!** mod ** 
 2004 WRITE(7,1000) ANGLEL+FLOAT(II)*DELTA,SUMISS(II)
 1000 FORMAT(3X,F6.2,1X,G12.5)
      WRITE(7,*) '|--------|------------|'
      WRITE(7,1021) WAVLEN
      WRITE(7,1023) QEXT2
      WRITE(7,1024) QSCA2
      WRITE(7,1025) QEXT2-QSCA2
 1021 FORMAT(G12.6,' micron = Wavelength')
 1023 FORMAT(G12.6,' = Extinction efficiency (Q)')
 1024 FORMAT(G12.6,' = Scattering efficiency (Qs)')
 1025 FORMAT(G12.6,' = Absorption efficiency (Qa)')
      IF(ANGLEU1.EQ.180.) WRITE(7,1006) 4.*PI*QEXT2/(QSCA2*SUMISS(NANG))
      write(7,1026) SUMNS*(WAVLEN**2./(4.*PI))/SUMN
 1026 format(G12.6,' micron^2 = <b> = Average surface area per particle' 
     +)
C
C     ** mod. **
c      WRITE(7,10010) SUMNS
c      WRITE(7,10011) SUMN
c      WRITE(7,10012) SUMNS*WAVLEN**2./(4*PI)
c10010 FORMAT(G12.6,' (*) = a = int{dN/dx*x^2}dx ')
c10011 FORMAT(G12.6,' particles / micron^3 (*) = N = int{dN/dx}dx ')
c10012 FORMAT(G12.6,' micron^2 (*) = A = a * [lambda^2/(4*pi)] = Area of
c     + all particles')
      WRITE(7,1001) QEXT1*WAVLEN**2.
      WRITE(7,1003) QSCA1*WAVLEN**2.
      WRITE(7,1004) (QEXT1-QSCA1)*WAVLEN**2.
      WRITE(7,1011)
      WRITE(7,*)'| COMMENTAIRES:                                      |'
      WRITE(7,1011)
      WRITE(7,*)'|   -Cross section / <b> should = efficiency         |'
      WRITE(7,*)'|    is independant of the normalisation of dN/dx    |'
      WRITE(7,*)'|   -The integral of the phase fct / 4 pi over all   |' 
      WRITE(7,*)'|    directions (solid angle) = 1.                   |'
      WRITE(7,*)'|   -The volume coefficient (k) is equal to C * N    |'
      WRITE(7,*)'|      = cross section * Number density              |'
      WRITE(7,*)'|      typical units: [um^2]*[particles/cm^3] =      |'
      WRITE(7,*)'|      [10^-3 km-1]                                  |'
      WRITE(7,*)'|      N = M/(d*4/3*pi*r^3) where d is the density   |'
      WRITE(7,*)'|      of particle, M is the mass per unit volume    |'
      WRITE(7,*)'|      and 4/3*pi*r^3 is the average volume  of a    |'
      WRITE(7,*)'|      particle.                                     |'
      WRITE(7,*)'|      typical units:  M -> [ug/m^3], d -> [g/cm^3]  |'
      WRITE(7,*)'|                      r -> [um]                     |'
      WRITE(7,*)'|      Lidar ratio= 4 pi Qe / phase_fct(180)*Qs      |'  
      WRITE(7,1011)
 1001 FORMAT(G12.6,' microns^2 = Extinction cross section (C)')
 1003 FORMAT(G12.6,' microns^2 = Scattering cross section (Cs)')
 1004 FORMAT(G12.6,' microns^2 = Absorption cross section (Ca)')
 1011 FORMAT(' |----------------------------------------------------|')
C     ** mod. **
C
 1006 FORMAT(G12.6,' = Lidar ratio')
 1306 FORMAT(1X,'THERE ARE NO PARTICLES IN THIS SIZE RANGE')
      RETURN
      END
C
C
C
C INTEGRATE COATED SPHERE PHASE FUNCTION OVER PARTICLE SIZE
C DISTRIBUTION BY SIMPSON'S METHOD
C BY BLAIR EVANS/DREV
C JANUARY 1988
      SUBROUTINE COATPHASE(X,Y,M,M2,CM)
      REAL*4 RI(0:18000),SUMISS(0:18000)
      REAL*8 PSD(5002),DUMR8(20504),RC,RIC,RT,RIT,RN,RIN,QEXT3
      REAL*8 QSCA3,QALPHA3,X3,Y3,WAVLEN
      COMPLEX*8 M2
      COMPLEX*16 M,CM
      COMMON /PHASE/ RI
      COMMON /DATA/ ANGLEL,ANGLEU,DELTA,NANG,ILOW,IUP,RINC,ANGLEU1
      COMMON /GLOBALR4/ SUMISS
      COMMON /GLOBALR8/ PSD,DUMR8,WAVLEN
      SAVE
      PI=3.141592653
      IFL=DINT(DBLE(CM))
      IDIFF=IUP-ILOW
      IF(IDIFF/2*2.EQ.IDIFF) GOTO10
      IUP=IUP+1
   10 C=X
      F=Y
      IF(Y.NE.0.) C=0.
      SUMNS=PSD(ILOW)*ILOW**2*RINC**2
      SUMN=PSD(ILOW)
      SUMA=SUMNS*ILOW*RINC
      Y=ILOW*RINC
      X=F*Y+C
      IF(IFL.EQ.-1) THEN
      CALL COAT(X,Y,M,M2,QEXT,QSCA,QALPHA)
      ELSE
      RC=DBLE(M)
      RIC=DIMAG(M)
      RT=REAL(M2)
      RIT=AIMAG(M2)
      RN=DBLE(CM)
      RIN=DIMAG(CM)
      X3=X
      Y3=Y
      CALL ANISO(RC,RIC,RT,RIT,RN,RIN,X3,Y3,QEXT3,QSCA3,QALPHA3)
      QEXT=QEXT3
      QSCA=QSCA3
      ENDIF
      QEXT1=QEXT*SUMNS
      QSCA1=QSCA*SUMNS
      DO 2001 I=0,NANG
 2001 SUMISS(I)=RI(I)*SUMNS
      DO 2000 I=ILOW+1,IUP,2
      RR=I*RINC
      R1=RR+RINC
      PROD1=PSD(I)*RR**2
      PROD2=PSD(I+1)*R1**2
      SUMNS=SUMNS+4.*PROD1+2.*PROD2
      SUMN=SUMN+4.*PSD(I)+2.*PSD(I+1)
      SUMA=SUMA+4.*PROD1*RR+2.*PROD2*R1
      Y=RR
      X=F*Y+C
      IF(IFL.EQ.-1) THEN
      CALL COAT(X,Y,M,M2,QEXT,QSCA,QALPHA)
      ELSE
      X3=X
      Y3=Y
      CALL ANISO(RC,RIC,RT,RIT,RN,RIN,X3,Y3,QEXT3,QSCA3,QALPHA3)
      QEXT=QEXT3
      QSCA=QSCA3
      ENDIF
      QEXT1=QEXT1+4.*QEXT*PROD1
      QSCA1=QSCA1+4.*QSCA*PROD1
      DO 2002 II=0,NANG
 2002 SUMISS(II)=SUMISS(II)+4.*RI(II)*PROD1
      Y=R1
      X=F*Y+C
      IF(IFL.EQ.-1) THEN
      CALL COAT(X,Y,M,M2,QEXT,QSCA,QALPHA)
      ELSE
      X3=X
      Y3=Y
      CALL ANISO(RC,RIC,RT,RIT,RN,RIN,X3,Y3,QEXT3,QSCA3,QALPHA3)
      QEXT=QEXT3
      QSCA=QSCA3
      ENDIF
      QEXT1=QEXT1+2.*QEXT*PROD2
      QSCA1=QSCA1+2.*QSCA*PROD2
      DO 2003 II=0,NANG
 2003 SUMISS(II)=SUMISS(II)+2.*RI(II)*PROD2
 2000 CONTINUE
      SUMNS=SUMNS-PROD2
      SUMN=SUMN-PSD(IUP)
      SUMA=SUMA-PROD2*RR
      QEXT1=QEXT1-QEXT*PROD2
      QSCA1=QSCA1-QSCA*PROD2
      IF(SUMNS.GT.0.) GOTO14
      WRITE(*,1306)
      STOP
   14 QEXT2=QEXT1/SUMNS
      QSCA2=QSCA1/SUMNS
      SUMA=QEXT1*3.*PI/2./SUMA
      QEXT1=.25*QEXT1/SUMN/PI
      QSCA1=.25*QSCA1/SUMN/PI
      WRITE(7,*) '|--------|------------|'
      WRITE(7,*) '| ANGLE  |  PHASE FCT |'
      WRITE(7,*) '|--------|------------|'
      DO 2004 II=0,NANG
      SUMISS(II)=(SUMISS(II)-RI(II)*PROD2)/SUMNS
 2004 WRITE(7,1000) ANGLEL+FLOAT(II)*DELTA,SUMISS(II)
 1000 FORMAT(3X,F6.2,1X,G12.5)
      WRITE(7,*) '|--------|------------|'
      WRITE(7,1021) WAVLEN
      WRITE(7,1023) QEXT2
      WRITE(7,1024) QSCA2
      WRITE(7,1025) QEXT2-QSCA2
 1021 FORMAT(G12.6,' micron = Wavelength')
 1023 FORMAT(G12.6,' = Extinction efficiency (Q)')
 1024 FORMAT(G12.6,' = Scattering efficiency (Qs)')
 1025 FORMAT(G12.6,' = Absorption efficiency (Qa)')
      IF(ANGLEU1.EQ.180.) WRITE(7,1006) 4.*PI*QEXT2/(QSCA2*SUMISS(NANG))
      write(7,1026) SUMNS*(WAVLEN**2./(4*PI))/SUMN
 1026 format(G12.6,' micron^2 = <b> = Average surface area per particle' 
     +)
C
C     ** mod. **
c      WRITE(7,10010) SUMNS
c      WRITE(7,10011) SUMN
c      WRITE(7,10012) SUMNS*WAVLEN**2./(4*PI)
c10010 FORMAT(G12.6,' (*) = a = int{dN/dx*x^2}dx ')
c10011 FORMAT(G12.6,' particles / micron^3 (*) = N = int{dN/dx}dx ')
c10012 FORMAT(G12.6,' micron^2 (*) = A = a * [lambda^2/(4*pi)] = Area of
c     + all particles')
      WRITE(7,1001) QEXT1*WAVLEN**2.
      WRITE(7,1003) QSCA1*WAVLEN**2.
      WRITE(7,1004) (QEXT1-QSCA1)*WAVLEN**2.
      WRITE(7,1011)
      WRITE(7,*)'| COMMENTAIRES:                                      |'
      WRITE(7,1011)
      WRITE(7,*)'|   -Cross section / <b> should = efficiency         |'
      WRITE(7,*)'|    is independant of the normalisation of dN/dx    |'
      WRITE(7,*)'|   -The integral of the phase fct / 4 pi over all   |' 
      WRITE(7,*)'|    directions (solid angle) = 1.                   |'
      WRITE(7,*)'|   -The volume coefficient (k) is equal to C * N    |'
      WRITE(7,*)'|      = cross section * Number density              |'
      WRITE(7,*)'|      typical units: [um^2]*[particles/cm^3] =      |'
      WRITE(7,*)'|      [10^-3 km-1]                                  |'
      WRITE(7,*)'|      N = M/(d*4/3*pi*r^3) where d is the density   |'
      WRITE(7,*)'|      of particle, M is the mass per unit volume    |'
      WRITE(7,*)'|      and 4/3*pi*r^3 is the average volume  of a    |'
      WRITE(7,*)'|      particle.                                     |'
      WRITE(7,*)'|      typical units:  M -> [ug/m^3], d -> [g/cm^3]  |'
      WRITE(7,*)'|                      r -> [um]                     |'
      WRITE(7,*)'|      Lidar ratio= 4 pi Qe / phase_fct(180)*Qs      |'
      WRITE(7,1011)
 1001 FORMAT(G12.6,' microns^2 = Extinction cross section (C)')
 1003 FORMAT(G12.6,' microns^2 = Scattering cross section (Cs)')
 1004 FORMAT(G12.6,' microns^2 = Absorption cross section (Ca)')
 1011 FORMAT(' |----------------------------------------------------|')
C     ** mod. **
C
 1006 FORMAT(G12.6,' = Lidar ratio')
 1306 FORMAT(1X,'THERE ARE NO PARTICLES IN THIS SIZE RANGE')
      RETURN
      END
C COATED SPHERE MIE CALCULATION (FROM BOHREN AND HUFFMAN 1983)
C MODIFIED FOR CALCULATING PHASE FUNCTION ADDITIONALLY
C BY BLAIR EVANS/DREV
C JANUARY 1988
      SUBROUTINE COAT(X,Y,M,M2,QEXT,QSCA,QALPHA)
      REAL*4 RI(0:18000)
      COMPLEX*8 M2
      COMPLEX*16 A(8750),B(8750),X1,X2,Y2,REFREL
      COMPLEX*16 D1X1,D0X1,D1X2,D0X2,D1Y2,D0Y2,QPHASE1,QPHASE2
      COMPLEX*16 XI0Y,XI1Y,XIY,CHI0Y2,CHI1Y2,CHIY2,CHI0X2,CHI1X2,CHIX2
      COMPLEX*16 CHIPX2,CHIPY2,ANCAP,BNCAP,DNBAR,GNBAR,AN,BN,CRACK,BRACK
      COMPLEX*16 AMESS1,AMESS2,AMESS3,AMESS4
      COMPLEX*16 M
      COMMON /PHASE/ RI
      COMMON /DATA/ ANGLEL,ANGLEU,DELTA,NANG,ILOW,IUP,RINC,ANGLEU1
      COMMON /POLAR/ IPARA,IPERP
      COMMON /GLOBALC16/ A,B
      SAVE
      IF(IPARA.EQ.0.AND.IPERP.EQ.0) GOTO1
      FACT=1./FLOAT(IPARA+IPERP)
    1 DEL=1.0E-8
      X1=M*X
      X2=M2*X
      Y2=M2*Y
      YSTOP=Y+4.*Y**.3333+2.
      REFREL=M2/M
      NSTOP=YSTOP
      D0X1=CDCOS(X1)/CDSIN(X1)
      D0X2=CDCOS(X2)/CDSIN(X2)
      D0Y2=CDCOS(Y2)/CDSIN(Y2)
      PSI0Y=COS(Y)
      PSI1Y=SIN(Y)
      CHI0Y=-SIN(Y)
      CHI1Y=COS(Y)
      XI0Y=CMPLX(PSI0Y,-CHI0Y)
      XI1Y=CMPLX(PSI1Y,-CHI1Y)
      CHI0Y2=-CDSIN(Y2)
      CHI1Y2=CDCOS(Y2)
      CHI0X2=-CDSIN(X2)
      CHI1X2=CDCOS(X2)
      QSCA=0.
      QEXT=0.
      N=1
      IFLAG=0
  200 RN=N
      PSIY=(2.*RN-1.)*PSI1Y/Y-PSI0Y
      CHIY=(2.*RN-1.)*CHI1Y/Y-CHI0Y
      XIY=CMPLX(PSIY,-CHIY)
      D1Y2=1./(RN/Y2-D0Y2)-RN/Y2
      IF(IFLAG.EQ.1) GOTO999
      D1X1=1./(RN/X1-D0X1)-RN/X1
      D1X2=1./(RN/X2-D0X2)-RN/X2
      CHIX2=(2.*RN-1.)*CHI1X2/X2-CHI0X2
      CHIY2=(2.*RN-1.)*CHI1Y2/Y2-CHI0Y2
      CHIPX2=CHI1X2-RN*CHIX2/X2
      CHIPY2=CHI1Y2-RN*CHIY2/Y2
      ANCAP=REFREL*D1X1-D1X2
      QPHASE2=REFREL*D1X1*CHIX2-CHIPX2
      IF(CDABS(QPHASE2).LT.1.E-25) IFLAG=1
      IF(QPHASE2.EQ.CMPLX(0.,0.)) QPHASE2=CMPLX(1.E-25,1.E-25)
      ANCAP=ANCAP/QPHASE2
      QPHASE1=CHIX2*D1X2-CHIPX2
      IF(CDABS(QPHASE1).LT.1.E-25) IFLAG=1
      IF(QPHASE1.EQ.CMPLX(0.,0.)) QPHASE1=CMPLX(1.E-25,1.E-25)
      ANCAP=ANCAP/QPHASE1
      BRACK=ANCAP*(CHIY2*D1Y2-CHIPY2)
      BNCAP=REFREL*D1X2-D1X1
      QPHASE2=REFREL*CHIPX2-D1X1*CHIX2
      IF(CDABS(QPHASE2).LT.1.E-25) IFLAG=1
      IF(QPHASE2.EQ.CMPLX(0.,0.)) QPHASE2=CMPLX(1.E-25,1.E-25)
      BNCAP=BNCAP/QPHASE2
      BNCAP=BNCAP/QPHASE1
      CRACK=BNCAP*(CHIY2*D1Y2-CHIPY2)
      AMESS1=BRACK*CHIPY2
      AMESS2=BRACK*CHIY2
      AMESS3=CRACK*CHIPY2
      AMESS4=CRACK*CHIY2
      IF(CDABS(AMESS1).GT.DEL*CDABS(D1Y2)) GOTO999
      IF(CDABS(AMESS2).GT.DEL) GOTO999
      IF(CDABS(AMESS3).GT.DEL*CDABS(D1Y2)) GOTO999
      IF(CDABS(AMESS4).GT.DEL) GOTO999
      BRACK=CMPLX(0.,0.)
      CRACK=CMPLX(0.,0.)
      IFLAG=1
  999 DNBAR=D1Y2-BRACK*CHIPY2
      DNBAR=DNBAR/(1.-BRACK*CHIY2)
      GNBAR=D1Y2-CRACK*CHIPY2
      GNBAR=GNBAR/(1.-CRACK*CHIY2)
      AN=(DNBAR/M2+RN/Y)*PSIY-PSI1Y
      AN=AN/((DNBAR/M2+RN/Y)*XIY-XI1Y)
      A(N)=AN
      BN=(M2*GNBAR+RN/Y)*PSIY-PSI1Y
      BN=BN/((M2*GNBAR+RN/Y)*XIY-XI1Y)
      B(N)=BN
      QSCA=QSCA+(2.*RN+1.)*(CDABS(AN)*CDABS(AN)+CDABS(BN)*CDABS(BN))
      QEXT=QEXT+(2.*RN+1.)*(DBLE(AN)+DBLE(BN))
      PSI0Y=PSI1Y
      PSI1Y=PSIY
      CHI0Y=CHI1Y
      CHI1Y=CHIY
      XI1Y=CMPLX(PSI1Y,-CHI1Y)
      CHI0X2=CHI1X2
      CHI1X2=CHIX2
      CHI0Y2=CHI1Y2
      CHI1Y2=CHIY2
      D0X1=D1X1
      D0X2=D1X2
      D0Y2=D1Y2
      N=N+1
      IF(N-1-NSTOP) 200,300,300
  300 QSCA=(2./(Y*Y))*QSCA
      QEXT=(2./(Y*Y))*QEXT
      QALPHA=4.71238898037*QEXT/Y
      DO 100 J=0,NANG
      RAD=3.14159265/180.*(ANGLEL+FLOAT(J)*DELTA)
      XL=COS(RAD)
      YL=SIN(RAD)
      PI1=1.
      PI2=3.*XL
      PP1=0.
      PP2=3.
      T2=XL*PI2-YL*YL*3.
      QPHASE1=1.5*(A(1)+B(1)*XL)+5./6.*(A(2)*PI2+B(2)*T2)
      QPHASE2=1.5*(B(1)+A(1)*XL)+5./6.*(B(2)*PI2+A(2)*T2)
      I=2
    3 I=I+1
      PI1=FLOAT(2*I-1)/FLOAT(I-1)*XL*PI2-FLOAT(I)/FLOAT(I-1)*PI1
      PP1=FLOAT(2*I-1)*PI2+PP1
      TAU=XL*PI1-YL*YL*PP1
      QPHASE1=QPHASE1+FLOAT(2*I+1)/FLOAT(I*(I+1))*(A(I)*PI1+B(I)*TAU)
      QPHASE2=QPHASE2+FLOAT(2*I+1)/FLOAT(I*(I+1))*(B(I)*PI1+A(I)*TAU)
      PTMP=PI1
      PI1=PI2
      PI2=PTMP
      PTMP=PP1
      PP1=PP2
      PP2=PTMP
      IF(I.LE.NSTOP) GOTO3
      QPH1=4./Y**2*CDABS(QPHASE1)**2
      QPH2=4./Y**2*CDABS(QPHASE2)**2
      QPH3=(QPH1*IPARA+QPH2*IPERP)*FACT
  100 RI(J)=QPH3/12.566371/QSCA
      RETURN
      END
C ANISOTROPIC COATED SPHERE MIE CALCULATION (FROM ROTH AND DIGNAM
C J. OPTICAL SOC. AM. V63 N3 1973) WITH CORRECTIONS.
C BY BLAIR EVANS/DREV
C NOVEMBER 1989
      SUBROUTINE ANISO(RC,RIC,RT,RIT,RN,RIN,X,Y,QEXT,QSCA,QALPHA)
      IMPLICIT REAL*8 (O-Z)
      IMPLICIT COMPLEX*16 (A-H)
      REAL*4 RI(0:18000),RANGLEL,RANGLEU,RDELTA,RINC,ANGLEU1
      COMPLEX*16 M,MT,MN,PSIMX(0:1100),PSIY(0:1100),PSIMTX(0:1100)
      COMPLEX*16 PSIMTY(0:1100),PSIMY(0:1100),AN(1100),BN(1100),LADD
      COMPLEX*16 DUMC16(9795)
      COMMON /PHASE/ RI
      COMMON /DATA/ RANGLEL,RANGLEU,RDELTA,NANG,ILOW,IUP,RINC,ANGLEU1
      COMMON /POLAR/ IPARA,IPERP
      COMMON /GLOBALC16/ PSIMX,PSIY,PSIMTX,PSIMTY,PSIMY,AN,BN,DUMC16
      SAVE
      PI=3.141592653589793238D0
      IF(IPARA.EQ.0.AND.IPERP.EQ.0) GOTO9
      FACT=1./DBLE(IPARA+IPERP)
    9 QSCA=0.D0
      QEXT=0.D0
      MT=DCMPLX(RT,RIT)
      MN=DCMPLX(RN,RIN)
      M=DCMPLX(RC,RIC)
      CMTX=MT*X
      CMX=M*X
      CMTY=MT*Y
      CMY=M*Y
      CKMTX=(PI*CMTX/2.D0)**.5
      CKMTY=(PI*CMTY/2.D0)**.5
      CKMY=(PI*CMY/2.D0)**.5
      CKMX=(PI*CMX/2.D0)**.5
      CKY=(PI*Y/2.D0)**.5
C BACKWARD RECURRENCE FOR HALF INTEGER ORDER BESSEL FUNCTIONS OF
C FIRST KIND + CONVERTING TO RICCATI-BESSEL FUNCTIONS.
      LY=Y+4.05*Y**.333333+2
      LX=LY
C PSI(MT*Y)
      Z1=DBLE(CMTY)
      Z2=DIMAG(CMTY)
      V1=DBLE(LY+1)+.5D0
      V2=0.
      CALL BESSEL(V1,V2,Z1,Z2,PSIMTY(LY+1))
      PSIMTY(LY+1)=CKMTY*CDEXP(PSIMTY(LY+1))
      V1=V1-1.D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMTY(LY))
      PSIMTY(LY)=CKMTY*CDEXP(PSIMTY(LY))
      DO 10 I=LY-1,0,-1
   10 PSIMTY(I)=DBLE(2*I+3)/CMTY*PSIMTY(I+1)-PSIMTY(I+2)
C PSI(Y)
      Z1=Y
      Z2=0.D0
      V1=DBLE(LY+1)+.5D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIY(LY+1))
      PSIY(LY+1)=CKY*CDEXP(PSIY(LY+1))
      V1=V1-1.D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIY(LY))
      PSIY(LY)=CKY*CDEXP(PSIY(LY))
      DO 20 I=LY-1,0,-1
   20 PSIY(I)=DBLE(2*I+3)/Y*PSIY(I+1)-PSIY(I+2)
C PSI(MTX)
      Z1=DBLE(CMTX)
      Z2=DIMAG(CMTX)
      V1=DBLE(LX+1)+.5D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMTX(LX+1))
      PSIMTX(LX+1)=CKMTX*CDEXP(PSIMTX(LX+1))
      V1=V1-1.D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMTX(LX))
      PSIMTX(LX)=CKMTX*CDEXP(PSIMTX(LX))
      DO 30 I=LX-1,0,-1
   30 PSIMTX(I)=DBLE(2*I+3)/CMTX*PSIMTX(I+1)-PSIMTX(I+2)
C PSI(MX)
      Z1=DBLE(CMX)
      Z2=DIMAG(CMX)
      V1=DBLE(LX+1)+.5D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMX(LX+1))
      PSIMX(LX+1)=CKMX*CDEXP(PSIMX(LX+1))
      V1=V1-1.D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMX(LX))
      PSIMX(LX)=CKMX*CDEXP(PSIMX(LX))
      DO 40 I=LX-1,0,-1
   40 PSIMX(I)=DBLE(2*I+3)/CMX*PSIMX(I+1)-PSIMX(I+2)
C PSI(MY)
      Z1=DBLE(CMY)
      Z2=DIMAG(CMY)
      V1=DBLE(LY+1)+.5D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMY(LY+1))
      PSIMY(LY+1)=CKMY*CDEXP(PSIMY(LY+1))
      V1=V1-1.D0
      CALL BESSEL(V1,V2,Z1,Z2,PSIMY(LY))
      PSIMY(LY)=CKMY*CDEXP(PSIMY(LY))
      DO 50 I=LY-1,0,-1
   50 PSIMY(I)=DBLE(2*I+3)/CMY*PSIMY(I+1)-PSIMY(I+2)
C GET INTIAL VALUES FOR NEGATIVE HALF INTEGER ORDER BESSEL FUNCTIONS
C OF FIRST KIND FOR FORWARD RECURRANCE + MAKING THEM RICCATI-BESSEL
C FUNCTIONS.
C CHI(Y)
      CHIY0=-DCMPLX(DCOS(Y),0.D0)
      CHIY1=-DCMPLX(DSIN(Y),0.D0)
C CHI(MTY)
      CMTY0=-CDCOS(CMTY)
      CMTY1=-CDSIN(CMTY)
C CHI(MY)
      CMY0=-CDCOS(CMY)
      CMY1=-CDSIN(CMY)
C CHI(MTX)
      CMTX0=-CDCOS(CMTX)
      CMTX1=-CDSIN(CMTX)
C START OF MAIN LOOP TO CALCULATE MIE COEFFICIENTS. INCLUDES
C CALCULATION OF BESSEL FUNCTIONS OF THE FIRST KIND OF COMPLEX ORDER.
      ZMTX1=DBLE(CMTX)
      ZMTX2=DIMAG(CMTX)
      ZMTY1=DBLE(CMTY)
      ZMTY2=DIMAG(CMTY)
      LL=1
      L=0
      CMTMN2=(MT/MN)**2
  100 L=L+1
      CW=(L*L*CMTMN2+L*CMTMN2+.25)**.5
      V1=DBLE(CW)
      V2=DIMAG(CW)
C CALCULATE NEGATIVE HALF INTEGER ORDER BESSEL FUNCTIONS OF FIRST KIND
C + RICCATI-BESSEL FUNCTIONS.
      CZ=(-2.*L+1)*CHIY0/Y-CHIY1
      CHIY1=CHIY0
      CHIY0=CZ
      CZ=(-2.*L+1)*CMTY0/CMTY-CMTY1
      CMTY1=CMTY0
      CMTY0=CZ
      CZ=(-2.*L+1)*CMY0/CMY-CMY1
      CMY1=CMY0
      CMY0=CZ
      CZ=(-2.*L+1)*CMTX0/CMTX-CMTX1
      CMTX1=CMTX0
      CMTX0=CZ
C CALCULATE POSITIVE AND NEGATIVE COMPLEX ORDER BESSEL FUNCTIONS OF
C FIRST KIND.  THE ORDER BEING CW. MAKE THEM RICCATI-BESSEL LIKE.
C CHI_W(MTX)
      CALL BESSEL(-V1,-V2,ZMTX1,ZMTX2,CWMTX)
      CWMTX=CDLOG(-CKMTX)+CWMTX
C CHI_W(MTY)
      CALL BESSEL(-V1,-V2,ZMTY1,ZMTY2,CWMTY)
      CWMTY=-CKMTY*CDEXP(CWMTY)
C PSI_W(MTX)
      CALL BESSEL(V1,V2,ZMTX1,ZMTX2,BWMTX)
      BWMTX=CKMTX*CDEXP(BWMTX)
C PSI_W(MTY)
      CALL BESSEL(V1,V2,ZMTY1,ZMTY2,BWMTY)
      BWMTY=CKMTY*CDEXP(BWMTY)
C CALCULATE THE DERIVATIVES OF ALL RICCATI-BESSEL FUNCTIONS (INCLUDING
C THOSE OF COMPLEX ORDER.
C PSI'(MTY)
      BMTYP=PSIMTY(L-1)-DBLE(L)*PSIMTY(L)/CMTY
C PSI'(Y)
      BYP=PSIY(L-1)-DBLE(L)*PSIY(L)/Y
C PSI'(MTX)
      BMTXP=PSIMTX(L-1)-DBLE(L)*PSIMTX(L)/CMTX
C PSI'(MX)
      BMXP=PSIMX(L-1)-DBLE(L)*PSIMX(L)/CMX
C PSI'(MY)
      BMYP=PSIMY(L-1)-DBLE(L)*PSIMY(L)/CMY
C CHI'(Y)
      CHIYP=-CHIY1-DBLE(L)*CHIY0/Y
C CHI'(MTY)
      CMTYP=-CMTY1-DBLE(L)*CMTY0/CMTY
C CHI'(MY)
      CMYP=-CMY1-DBLE(L)*CMY0/CMY
C CHI'(MTX)
      CMTXP=-CMTX1-DBLE(L)*CMTX0/CMTX
C PSI_W'(MTX)
      CALL BESSEL(V1-1.,V2,ZMTX1,ZMTX2,BWMTXP)
      BWMTXP=CKMTX*CDEXP(BWMTXP)+(.5D0-CW)*BWMTX/CMTX
C PSI_W'(MTY)
      CALL BESSEL(V1-1.,V2,ZMTY1,ZMTY2,BWMTYP)
      BWMTYP=CKMTY*CDEXP(BWMTYP)+(.5D0-CW)*BWMTY/CMTY
C CHI_W'(MTX)
      CALL BESSEL(-V1-1.,-V2,ZMTX1,ZMTX2,CWMTXP)
      ALX=CDLOG(-CKMTX)+CWMTXP
      ALY=CDLOG((.5D0+CW))+CWMTX-CDLOG(CMTX)
      CWMTXP=LADD(ALX,ALY,1)
C CHI_W'(MTY)
      CALL BESSEL(-V1-1.,-V2,ZMTY1,ZMTY2,CWMTYP)
      CWMTYP=-CKMTY*CDEXP(CWMTYP)+(.5D0+CW)*CWMTY/CMTY
      LL=-LL
C COMPUTE ZETA(Y) AND ZETA(MY) AND THEIR DERIVATIVES.
      CZY=PSIY(L)+CHIY0*DCMPLX(0.D0,1.D0)*LL
C      CZMY=PSIMY(L)+CMY0*DCMPLX(0.D0,1.D0)*LL
      CZYP=BYP+CHIYP*DCMPLX(0.D0,1.D0)*LL
C      CZMYP=BMYP+CMYP*DCMPLX(0.D0,1.D0)*LL
C CALCULATE WRONSKIANS OF PSI, CHI AND ZETA AND THEIR MIXTURES.
C WITH CORRECTIONS (INCLUDING ERRATA J. OPTICAL SOC. AM. V66 N9 P981)
      CPPXP=M*BWMTXP*PSIMX(L)-MT*BWMTX*BMXP
      CCPQXP=CWMTYP*PSIY(L)-MT*CWMTY*BYP
      IF(PSIMX(L).EQ.DCMPLX(0.D0,0.D0).OR.BMXP.EQ.DCMPLX(0.D0,0.D0))
     +THEN
      CCPXP=DCMPLX(0.D0,0.D0)
      ELSE
      ALX=CDLOG(M)+CWMTXP+CDLOG(PSIMX(L))
      ALY=CDLOG(MT)+CWMTX+CDLOG(BMXP)
      CCPXP=CDEXP(LADD(ALX,ALY,-1))
      ENDIF
      CPPQXP=BWMTYP*PSIY(L)-MT*BWMTY*BYP
      CCZQXP=CWMTYP*CZY-MT*CWMTY*CZYP
      CPZQXP=BWMTYP*CZY-MT*BWMTY*CZYP
      CPPXD=MT*BMTXP*PSIMX(L)-M*PSIMTX(L)*BMXP
      CCPQXD=MT*CMTYP*PSIY(L)-CMTY0*BYP
      CCPXD=MT*CMTXP*PSIMX(L)-M*CMTX0*BMXP
      CPPQXD=MT*BMTYP*PSIY(L)-PSIMTY(L)*BYP
      CCZQXD=MT*CMTYP*CZY-CMTY0*CZYP
      CPZQXD=MT*BMTYP*CZY-PSIMTY(L)*CZYP
C CALCULATE THE ANISOTROPIC MIE COEFFICIENTS.
      IF(CPPXP.EQ.DCMPLX(0.D0,0.D0).AND.CCPXP.EQ.DCMPLX(0.D0,0.D0)) THEN
      CRAT=DCMPLX(0.D0,0.D0)
      ELSE
      CRAT=CPPXP/CCPXP
      ENDIF
      IF(CPPXD.EQ.DCMPLX(0.D0,0.D0).AND.CCPXD.EQ.DCMPLX(0.D0,0.D0)) THEN
      CRAT1=DCMPLX(0.D0,0.D0)
      ELSE
      CRAT1=CPPXD/CCPXD
      ENDIF
      AN(L)=(CRAT*CCPQXP-CPPQXP)/(CRAT*CCZQXP-CPZQXP)
      BN(L)=(CRAT1*CCPQXD-CPPQXD)/(CRAT1*CCZQXD-CPZQXD)
C CALCULATE THE EFFICIENCY FACTORS.
      QSCA=QSCA+(2.*DBLE(L)+1.D0)*(CDABS(AN(L))**2+CDABS(BN(L))**2)
      QEXT=QEXT+(2.*DBLE(L)+1.D0)*(DBLE(AN(L))+DBLE(BN(L)))
      IF(L-LY) 100,200,100
  200 QSCA=(2.D0/(Y*Y))*QSCA
      QEXT=(2.D0/(Y*Y))*QEXT
      QALPHA=4.71238898037/Y*QEXT
      DO 101 J=0,NANG
      RAD=3.1415926535897939/180.D0*(RANGLEL+DBLE(J)*RDELTA)
      XL=DCOS(RAD)
      YL=DSIN(RAD)
      PI1=1.D0
      PI2=3.D0*XL
      PP1=0.D0
      PP2=3.D0
      T2=XL*PI2-YL*YL*3.D0
C AND CALCULATE PHASE FUNCTION (CPHA1 AND CPHA2) FOR BOTH
C POLARIZATION STATES AND THE AVERAGE QPH3
      CPHA1=1.5D0*(AN(1)+BN(1)*XL)+5.D0/6.D0*(AN(2)*PI2+BN(2)*T2)
      CPHA2=1.5D0*(BN(1)+AN(1)*XL)+5.D0/6.D0*(BN(2)*PI2+AN(2)*T2)
      I=2
    3 IF (I.LT.Y) GOTO4
      IF (I.GE.LY) GOTO75
      IF (DBLE(AN(I)).LE.1.D-17) GOTO75
    4 I=I+1
      PI1=DBLE(2*I-1)/DBLE(I-1)*XL*PI2-DBLE(I)/DBLE(I-1)*PI1
      PP1=DBLE(2*I-1)*PI2+PP1
      TAU=XL*PI1-YL*YL*PP1
      CPHA1=CPHA1+DBLE(2*I+1)/DBLE(I*(I+1))*(AN(I)*PI1+BN(I)*TAU)
      CPHA2=CPHA2+DBLE(2*I+1)/DBLE(I*(I+1))*(AN(I)*TAU+BN(I)*PI1)
      PTMP=PI1
      PI1=PI2
      PI2=PTMP
      PTMP=PP1
      PP1=PP2
      PP2=PTMP
      GOTO3
   75 QPH1=4.D0/Y**2*CDABS(CPHA1)**2
      QPH2=4.D0/Y**2*CDABS(CPHA2)**2
      QPH3=(QPH1*IPARA+QPH2*IPERP)*FACT
  101 RI(J)=QPH3/12.566371/QSCA
      RETURN
      END
      FUNCTION LADD(ALX,ALY,IS)
      COMPLEX*16 LADD,ALX,ALY
      ILX=DINT(DBLE(ALX))
      ALY=ALY-ILX
      ALX=ALX-ILX
      IF(CDABS(ALY).LT.42.D0) THEN
      LADD=CDLOG(CDEXP(ALX)+IS*CDEXP(ALY))+ILX
      ELSE
      IF(DBLE(ALY).GT.DBLE(ALX)) THEN
      LADD=ALY+ILX
      ELSE
      LADD=ALX+ILX
      ENDIF
      ENDIF
      END
C SUBROUTINE TO CALCULATE BESSEL FUNCTION OF COMPLEX ORDER AND
C ARGUMENT.
      SUBROUTINE BESSEL(V1X,V2,Z3,Z4,LJV)
      REAL*8 V1,V2,Z1,Z2,V3,PI,Z3,Z4,CS2,V1X
      COMPLEX*16 V,Z,JV,RNUM,RDEN,T,RAT,CHI,MU,P,Q,T1,T2,T3,Z64,CS,LJV
      SAVE
      V1=V1X
      LJV=DCMPLX(0.D0,0.D0)
      IF(V2.NE.0..OR.DINT(V1).NE.V1.OR.V1.GE.0.) GOTO3
      Z3=-Z3
      V1=-V1X
    3 Z1=DABS(Z3)
      Z2=DABS(Z4)
      V=CMPLX(V1,V2)
      Z=CMPLX(Z1,Z2)
      PI=3.1415926535897932D0
      IF(Z1.EQ.0..AND.Z2.EQ.0.) GOTO50
      IF (CDABS(V*V/Z).GT.10..OR.CDABS(Z).LE.15.) GOTO60
C BESSEL FUNCTION OF FIRST KIND COMPLEX ORDER AND ARGUMENT
C FROM A+S ASYMPTOTIC EXPANSIONS 9.2.5
      K=INT(CDABS(V*V/Z))+5
      Z64=64*Z**2
      CHI=Z-(.5D0*V+.25D0)*PI
      MU=4.D0*V*V
      P=DCMPLX(1.D0,0.D0)
      Q=(MU-1.D0)/8/Z
      T3=DCMPLX(1.D0,0.D0)
      T2=(MU-1)/8/Z
      DO 70 I=1,K
      T1=(MU-(4*I-3)**2)*(MU-(4*I-1)**2)/(2*I-1)/(2*I)/Z64*(-T3)
      T3=T1
      T2=(MU-(4*I-1)**2)*(MU-(4*I+1)**2)/(2*I)/(2*I+1)/Z64*(-T2)
      P=P+T1
      Q=Q+T2
   70 CONTINUE
      IF(Z2.GT.345.) GOTO30
      JV=(2/PI/Z)**.5*(P*CDCOS(CHI)-Q*CDSIN(CHI))
      LJV=CDLOG(JV)
      GOTO15
   30 CONTINUE
      CS2=(.25+.5*V1)*PI-Z1
      CS2=CS2-2.*PI*DINT((CS2+DSIGN(PI,CS2))/2./PI)
      CS=DCMPLX(Z2-.5*V2*PI-DLOG(2.D0),CS2)
      LJV=.5*CDLOG(2./PI/Z)+CDLOG(P)+CS+CDLOG((1.D0-Q/P*
     +DCMPLX(0.D0,1.D0)))
      GOTO15
C USE BACKWARD RECURRANCE AND NEUMANN EXPANSION FOR NORMALIZATION
C MAKING SURE THAT ROUTINE IS IN A SAFE REGION.
   60 N=3
      IVT=INT(CDABS(V))+10
      V3=V1
      IF(V1.LT.0.) V3=V1-INT(V1)
      IMAG=.01368804*Z2*Z2+22.*(1.-1.826*EXP(-Z2*Z2/1200.))+.5
      IREAL=.02653561*Z2*Z2-17.6*(1.-EXP(-Z2/21.8))+.5
      IF(V2.LT.IMAG.AND.V1.LT.IREAL) V3=IREAL+V3-INT(V3)+2.
      IF(V3.EQ.0.AND.V2.EQ.0.) V3=V3+1.
   20 CALL CBESSEL(V3,V2,Z1,Z2,LJV,ERR)
C	  WRITE(*,*) LJV,ERR
      IF(ERR.GT.10.0) GOTO10
      V3=V3+IVT
      GOTO20
   10 IDIFF=V3-V1+.1
C IF BACKWARD RECURRANCE WAS NOT INTIALLY IN A GOOD REGION THEN
C USE CONTINUED FRACTIONS TO GET THERE.
      IF(IDIFF.EQ.0) GOTO15
      DO 40 I=1,IDIFF
      V=CMPLX(V3,V2)-I+1
      N=3
      RNUM=2.*V/Z
      RDEN=-2.*(V+1.)/Z
      RAT=RNUM
      RNUM=RDEN+1./RNUM
      RAT=RAT*RNUM/RDEN
   18 T=(-1)**(N+1)*2.*(V+N-1)/Z
      RNUM=T+1./RNUM
      RDEN=T+1./RDEN
      IF(CDABS(RNUM-RDEN).LT.1.E-13) GOTO100
      N=N+1
      RAT=RAT*RNUM/RDEN
      GOTO18
  100 CONTINUE
C     WRITE(*,*) LJV,RAT
   40 LJV=LJV+CDLOG(RAT)
   15 IF(Z1.GT.Z3.AND.Z2.GT.Z4) LJV=LJV+DCMPLX(V1,V2)*PI*
     +DCMPLX(0.D0,1.D0)
      IF(Z1.GT.Z3.AND.Z2.EQ.Z4) LJV=DCONJG(LJV+DCMPLX(V1,V2)
     +*PI*DCMPLX(0.D0,1.D0))
      IF(Z1.EQ.Z3.AND.Z2.GT.Z4) LJV=DCONJG(LJV)
C     WRITE(*,*) LJV
C     IF(ABS(DBLE(LJV)).LT.345.) WRITE(*,*) CDEXP(LJV)
      RETURN
   50 IF(V1.NE.0..OR.V2.NE.0) GOTO55
C     WRITE(*,*) (1.D0,0.D0)
      LJV=0.
      RETURN
C  55 WRITE(*,*) (0.D0,0.D0)
   55 LJV=-1.D150
      RETURN
      END
C SUBROUTINE FOR CALCULATING BESSEL FUNCTIONS OF COMPLEX
C ORDER AND ARGUMENT.  BACKWARD RECURSION AND GENERALIZED
C NEUMANN EXPANSIONS (A+S 9.1.87).  ERR RETURNS ACCURACY
C ESTIMATE IF LESS THAN 10 DIGITS.
      SUBROUTINE CBESSEL(V1,V2,Z1,Z2,JV,ERR)
C     REAL*8 GAM1,V1,V2,Z1,Z2,V3,V4,PI,X,Y,YY,SUMM,ERR
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 ERR
      COMPLEX*16 V,Z,JV,JV1,JV2,CN,SUM,GAM,T1,T2
	  SAVE
      ERR=100.
      PI=3.141592653589793238D0
      V=DCMPLX(V1,V2)
      Z=DCMPLX(Z1,Z2)
      Y=CDABS(V)
      X=CDABS(Z)
C L IS NUMBER OF TERMS IN NEUMANN EXPANSION REQUIRED FOR CONVERGENCE
C (DERIVED EMPIRICALLY).
      L=2*DINT((X+5.5*DLOG10(X)+3.*DSQRT(X)+8.5)/2.)
      IF (L.LT.12) L=12
      IF (Y.GT.5) GOTO30
      V3=DABS(V1+L/2)
      V4=DABS(V2)
      CN=DCMPLX(V3,V4)
      JV1=DCMPLX(0.D0,0.D0)
      JV2=DCMPLX(1.D0,0.D0)
      GAM=(CN-.5D0)*CDLOG(CN)-CN+.5D0*DLOG(6.2831853071795865D0)+
     +1.D0/12/CN-1.D0/360/CN**3+1.D0/1260/CN**5-1.D0/1680/CN**7+
     +1.D0/1188/CN**9-.001917526917D0/CN**11+1.D0/156/CN**13
      GAM1=0.D0
      DO 5 I=1,L/2
    5 GAM1=GAM1+DLOG(DBLE(I))
      IF(V3.GT.V1+L/2) GAM=DCMPLX(0.D0,PI)+DLOG(PI)-GAM-CDLOG(
     +CN*CDSIN(PI*CN))
      IF(V3.GT.V1+L/2.AND.V4.EQ.V2) GAM=DCONJG(GAM)
      IF(V3.EQ.V1+L/2.AND.V4.GT.V2) GAM=DCONJG(GAM)
      T2=DCMPLX(1.D0,0.D0)
      SUM=(V+L)*T2
      CN=V+L+1
      DO 10 I=1,L
      T1=2*(CN-I)/Z
      JV=T1*JV2-JV1
      JV1=JV2
      JV2=JV
      IF (((L-I)/2)*2.NE.(L-I)) GOTO 10
      K=(L-I)/2
      T2=T2*(K+1)/(V+K)
      SUM=SUM+(V+L-I)*JV*T2
   10 CONTINUE
      T2=V*CDLOG(.5*Z)-GAM+GAM1
      JV=CDLOG(JV)-CDLOG(SUM)+T2
      GOTO25
   30 V3=DABS(V1)
      V4=DABS(V2)
      CN=DCMPLX(V3,V4)
      JV1=DCMPLX(0.D0,0.D0)
      JV2=DCMPLX(1.D-20,0.D0)
      GAM=(CN-.5D0)*CDLOG(CN)-CN+.5D0*DLOG(6.2831853071795865D0)+
     +1.D0/12/CN-1.D0/360/CN**3+1.D0/1260/CN**5-1.D0/1680/CN**7
      IF(CDABS(CN).GT.100.) GOTO1
      GAM=GAM+1.D0/1188/CN**9-.001917526917D0/CN**11+1.D0/156/CN**13
    1 IF(V3.GT.V1) GAM=DCMPLX(0.D0,PI)+DLOG(PI)-GAM-CDLOG(CN*
     +CDSIN(PI*CN))
      IF (V3.GT.V1.AND.V4.EQ.V2) GAM=DCONJG(GAM)
      IF (V3.EQ.V1.AND.V4.GT.V2) GAM=DCONJG(GAM)
      SUM=(V+L)*1.E-20
      CN=V+L+1
      FAC=0.
      SUMM=-1.D10
      IFAC=0
      DO 15 I=1,L
      T1=2*(CN-I)/Z
      JV=T1*JV2-JV1
      JV1=JV2
      JV2=JV
      IF (((L-I)/2)*2.NE.(L-I)) GOTO15
      IF (CDABS(JV).GT.1.E-32) GOTO2
      FAC=FAC+1.
      JV1=JV1*1.E30
      JV2=JV2*1.E30
      JV=JV*1.E30
    2 K=(L-I)/2
      T2=(K+1)/(V+K)
      JV1=T2*JV1
      JV2=T2*JV2
      IF (FAC.NE.0.) GOTO15
      SUM=SUM+(V+L-I)*JV*T2
      YY=DLOG(CDABS(SUM))+IFAC*69.077553
      IF(YY-SUMM.LT.-16.118) ERR=17.D0+DLOG10(DEXP(YY-SUMM))
      IF(YY.GT.SUMM) SUMM=YY
      IF (CDABS(SUM).LT.1.E30) GOTO15
      IFAC=IFAC+1
      SUM=SUM*1.E-30
      JV1=JV1*1.E-30
      JV2=JV2*1.E-30
      JV=JV*1.E-30
   15 CONTINUE
      T2=V*CDLOG(.5*Z)
      JV=(CDLOG(JV/V)-GAM-CDLOG(SUM)+T2-FAC*69.0775527898213D0)
   25 CONTINUE
      RETURN
      END
