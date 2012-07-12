

********************************************************************
*                                                                  *
*        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                 *
*  D.de Florian, R.Sassot, M.Stratmann   hep-ph/0703242            *
*                                                                  *
*     CALL fDSS (IH,IC,IO, X, Q2, U, UB, D, DB, S, SB, C, B, GL)   *
*                                                                  *	
*  INPUT:                                                          *
*  IH = hadron type    1: PION                                     *
*                      2: KAON                                     *
*                      3: PROTON                                   *
*                      4: CHARGED HADRONS                          *
*                                                                  *
*  IC = Hadron Charge  0: 0 (as average of + and -)                *
*                      1: +                                        *
*                     -1: -                                        *
*                                                                  *
*  IO= Order           0: LO                                       *
*                      1: NLO                                      *
*                                                                  *
*            X                    (between  0.05   and  1.0)       *
*            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      *
*             (for values outside the allowed range the program    *
*              writes a warning and extrapolates to the x and      *
*              Q2 values requested)                                *
*                                                                  *
*   OUTPUT: U, UB, D, DB, S, SB,   C,           B,       GL        *
*           U Ubar D Dbar S Sbar Charm=Cbar Bottom=Bbar Gluon      *
*           Always X times the distribution is returned            *
*                                                                  *
*                                                                  *
*   COMMON:  The main program or the calling routine has to have   *
*            a common block  COMMON / FRAGINI / FINI , and  FINI   *
*            has always to be zero when DSS is called for the      *
*            first time or when the SET has been changed.          *
*                                                                  *
********************************************************************

c 10/23/08
c JGL - modified to calculate and return only requested frag function
c ityp argument follows KKP convention for parton numbers:
c     Parton label:
c     0    1    2    3    4    5    6    7    8     9    10
c     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
c NOTE: cbar and bbar not in DSS, returns c or b instead
c DH is the fragmentation function (X factor divided out)
c

!> DSS fragmentation function routine.
!! 
!!  INPUT:                                                       
!! IH = hadron type    1: PION                                  
!!                     2: KAON                                  
!!                     3: PROTON                                
!!                     4: CHARGED HADRONS                       
!!                                                              
!! IC = Hadron Charge  0: 0 (as average of + and -)             
!!                     1: +                                     
!!                    -1: -                                     
!!                                                              
!! IO= Order           0: LO                                    
!!                     1: NLO                                   
!!                                                              
!!           X                    (between  0.05   and  1.0)    
!!           Q2 = scale in GeV**2 (between  1.0    and  1.D5)   
!!            (for values outside the allowed range the program 
!!             writes a warning and extrapolates to the x and   
!!             Q2 values requested)                             
!!                                                              
!! ityp argument follows KKP convention for parton numbers:
!!
!!     Parton label:
!!
!!     0    1    2    3    4    5    6    7    8     9    10
!!
!!     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
!!
!! NOTE: cbar and bbar not in DSS, returns c or b instead
!!
!! DH is the fragmentation function (X factor divided out)!!                                                              
!!
!!  COMMON:  The main program or the calling routine has to have
!!           a common block  COMMON / FRAGINI / FINI , and  FINI
!!           has always to be zero when DSS is called for the   
!!           first time or when the SET has been changed.       
      SUBROUTINE fDSS (IH,IC,IO, X, Q2, ITYP, DH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPART=9, NX=35, NQ=24, NARG=2)
      DIMENSION PARTON (NPART,NQ,NX-1)
      DIMENSION QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ)

      DIMENSION XUTOTFS(NX,NQ,4,0:1), XDTOTFS(NX,NQ,4,0:1)
      DIMENSION XSTOTFS(NX,NQ,4,0:1)
      DIMENSION XUVALFS(NX,NQ,4,0:1), XDVALFS(NX,NQ,4,0:1) 
      DIMENSION XSVALFS(NX,NQ,4,0:1)
      DIMENSION XCTOTFS(NX,NQ,4,0:1), XBTOTFS(NX,NQ,4,0:1)
      DIMENSION XGFS(NX,NQ,4,0:1)

      COMMON / FRAGINI / FINI

      SAVE XUTOTFS, XDTOTFS, XSTOTFS, XCTOTFS, XBTOTFS, XGFS, NA, ARRF
      SAVE XUVALFS, XDVALFS, XSVALFS

*     Equivalence statements allow quick access to different parts of
*     the fragmentation tables
      
*     LO pions

      DIMENSION XUTOTFPILO(NX,NQ), XDTOTFPILO(NX,NQ), XSTOTFPILO(NX,NQ)
      DIMENSION XUVALFPILO(NX,NQ), XDVALFPILO(NX,NQ), XSVALFPILO(NX,NQ)
      DIMENSION XCTOTFPILO(NX,NQ), XBTOTFPILO(NX,NQ)
      DIMENSION XGFPILO(NX,NQ)
      EQUIVALENCE(XUTOTFPILO,XUTOTFS(1,1,1,0))
      EQUIVALENCE(XDTOTFPILO,XDTOTFS(1,1,1,0))
      EQUIVALENCE(XSTOTFPILO,XSTOTFS(1,1,1,0))
      EQUIVALENCE(XUVALFPILO,XUVALFS(1,1,1,0))
      EQUIVALENCE(XDVALFPILO,XDVALFS(1,1,1,0))
      EQUIVALENCE(XSVALFPILO,XSVALFS(1,1,1,0))
      EQUIVALENCE(XCTOTFPILO,XCTOTFS(1,1,1,0))
      EQUIVALENCE(XBTOTFPILO,XBTOTFS(1,1,1,0))
      EQUIVALENCE(XGFPILO,XGFS(1,1,1,0))

*     LO kaons

      DIMENSION XUTOTFKALO(NX,NQ), XDTOTFKALO(NX,NQ), XSTOTFKALO(NX,NQ)
      DIMENSION XUVALFKALO(NX,NQ), XDVALFKALO(NX,NQ), XSVALFKALO(NX,NQ)
      DIMENSION XCTOTFKALO(NX,NQ), XBTOTFKALO(NX,NQ)
      DIMENSION XGFKALO(NX,NQ)
      EQUIVALENCE(XUTOTFKALO,XUTOTFS(1,1,2,0))
      EQUIVALENCE(XDTOTFKALO,XDTOTFS(1,1,2,0))
      EQUIVALENCE(XSTOTFKALO,XSTOTFS(1,1,2,0))
      EQUIVALENCE(XUVALFKALO,XUVALFS(1,1,2,0))
      EQUIVALENCE(XDVALFKALO,XDVALFS(1,1,2,0))
      EQUIVALENCE(XSVALFKALO,XSVALFS(1,1,2,0))
      EQUIVALENCE(XCTOTFKALO,XCTOTFS(1,1,2,0))
      EQUIVALENCE(XBTOTFKALO,XBTOTFS(1,1,2,0))
      EQUIVALENCE(XGFKALO,XGFS(1,1,2,0))

*     LO protons

      DIMENSION XUTOTFPRLO(NX,NQ), XDTOTFPRLO(NX,NQ), XSTOTFPRLO(NX,NQ)
      DIMENSION XUVALFPRLO(NX,NQ), XDVALFPRLO(NX,NQ), XSVALFPRLO(NX,NQ)
      DIMENSION XCTOTFPRLO(NX,NQ), XBTOTFPRLO(NX,NQ)
      DIMENSION XGFPRLO(NX,NQ)
      EQUIVALENCE(XUTOTFPRLO,XUTOTFS(1,1,3,0))
      EQUIVALENCE(XDTOTFPRLO,XDTOTFS(1,1,3,0))
      EQUIVALENCE(XSTOTFPRLO,XSTOTFS(1,1,3,0))
      EQUIVALENCE(XUVALFPRLO,XUVALFS(1,1,3,0))
      EQUIVALENCE(XDVALFPRLO,XDVALFS(1,1,3,0))
      EQUIVALENCE(XSVALFPRLO,XSVALFS(1,1,3,0))
      EQUIVALENCE(XCTOTFPRLO,XCTOTFS(1,1,3,0))
      EQUIVALENCE(XBTOTFPRLO,XBTOTFS(1,1,3,0))
      EQUIVALENCE(XGFPRLO,XGFS(1,1,3,0))

*     LO charged hadrons

      DIMENSION XUTOTFCHLO(NX,NQ), XDTOTFCHLO(NX,NQ), XSTOTFCHLO(NX,NQ)
      DIMENSION XUVALFCHLO(NX,NQ), XDVALFCHLO(NX,NQ), XSVALFCHLO(NX,NQ)
      DIMENSION XCTOTFCHLO(NX,NQ), XBTOTFCHLO(NX,NQ)
      DIMENSION XGFCHLO(NX,NQ)
      EQUIVALENCE(XUTOTFCHLO,XUTOTFS(1,1,4,0))
      EQUIVALENCE(XDTOTFCHLO,XDTOTFS(1,1,4,0))
      EQUIVALENCE(XSTOTFCHLO,XSTOTFS(1,1,4,0))
      EQUIVALENCE(XUVALFCHLO,XUVALFS(1,1,4,0))
      EQUIVALENCE(XDVALFCHLO,XDVALFS(1,1,4,0))
      EQUIVALENCE(XSVALFCHLO,XSVALFS(1,1,4,0))
      EQUIVALENCE(XCTOTFCHLO,XCTOTFS(1,1,4,0))
      EQUIVALENCE(XBTOTFCHLO,XBTOTFS(1,1,4,0))
      EQUIVALENCE(XGFCHLO,XGFS(1,1,4,0))

*     NLO pions

      DIMENSION XUTOTFPINLO(NX,NQ),XDTOTFPINLO(NX,NQ),XSTOTFPINLO(NX,NQ)
      DIMENSION XUVALFPINLO(NX,NQ),XDVALFPINLO(NX,NQ),XSVALFPINLO(NX,NQ)
      DIMENSION XCTOTFPINLO(NX,NQ), XBTOTFPINLO(NX,NQ)
      DIMENSION XGFPINLO(NX,NQ)
      EQUIVALENCE(XUTOTFPINLO,XUTOTFS(1,1,1,1))
      EQUIVALENCE(XDTOTFPINLO,XDTOTFS(1,1,1,1))
      EQUIVALENCE(XSTOTFPINLO,XSTOTFS(1,1,1,1))
      EQUIVALENCE(XUVALFPINLO,XUVALFS(1,1,1,1))
      EQUIVALENCE(XDVALFPINLO,XDVALFS(1,1,1,1))
      EQUIVALENCE(XSVALFPINLO,XSVALFS(1,1,1,1))
      EQUIVALENCE(XCTOTFPINLO,XCTOTFS(1,1,1,1))
      EQUIVALENCE(XBTOTFPINLO,XBTOTFS(1,1,1,1))
      EQUIVALENCE(XGFPINLO,XGFS(1,1,1,1))

*     NLO kaons

      DIMENSION XUTOTFKANLO(NX,NQ),XDTOTFKANLO(NX,NQ),XSTOTFKANLO(NX,NQ)
      DIMENSION XUVALFKANLO(NX,NQ),XDVALFKANLO(NX,NQ),XSVALFKANLO(NX,NQ)
      DIMENSION XCTOTFKANLO(NX,NQ), XBTOTFKANLO(NX,NQ)
      DIMENSION XGFKANLO(NX,NQ)
      EQUIVALENCE(XUTOTFKANLO,XUTOTFS(1,1,2,1))
      EQUIVALENCE(XDTOTFKANLO,XDTOTFS(1,1,2,1))
      EQUIVALENCE(XSTOTFKANLO,XSTOTFS(1,1,2,1))
      EQUIVALENCE(XUVALFKANLO,XUVALFS(1,1,2,1))
      EQUIVALENCE(XDVALFKANLO,XDVALFS(1,1,2,1))
      EQUIVALENCE(XSVALFKANLO,XSVALFS(1,1,2,1))
      EQUIVALENCE(XCTOTFKANLO,XCTOTFS(1,1,2,1))
      EQUIVALENCE(XBTOTFKANLO,XBTOTFS(1,1,2,1))
      EQUIVALENCE(XGFKANLO,XGFS(1,1,2,1))

*     NLO protons

      DIMENSION XUTOTFPRNLO(NX,NQ),XDTOTFPRNLO(NX,NQ),XSTOTFPRNLO(NX,NQ)
      DIMENSION XUVALFPRNLO(NX,NQ),XDVALFPRNLO(NX,NQ),XSVALFPRNLO(NX,NQ)
      DIMENSION XCTOTFPRNLO(NX,NQ), XBTOTFPRNLO(NX,NQ)
      DIMENSION XGFPRNLO(NX,NQ)
      EQUIVALENCE(XUTOTFPRNLO,XUTOTFS(1,1,3,1))
      EQUIVALENCE(XDTOTFPRNLO,XDTOTFS(1,1,3,1))
      EQUIVALENCE(XSTOTFPRNLO,XSTOTFS(1,1,3,1))
      EQUIVALENCE(XUVALFPRNLO,XUVALFS(1,1,3,1))
      EQUIVALENCE(XDVALFPRNLO,XDVALFS(1,1,3,1))
      EQUIVALENCE(XSVALFPRNLO,XSVALFS(1,1,3,1))
      EQUIVALENCE(XCTOTFPRNLO,XCTOTFS(1,1,3,1))
      EQUIVALENCE(XBTOTFPRNLO,XBTOTFS(1,1,3,1))
      EQUIVALENCE(XGFPRNLO,XGFS(1,1,3,1))

*     NLO charged hadrons

      DIMENSION XUTOTFCHNLO(NX,NQ),XDTOTFCHNLO(NX,NQ),XSTOTFCHNLO(NX,NQ)
      DIMENSION XUVALFCHNLO(NX,NQ),XDVALFCHNLO(NX,NQ),XSVALFCHNLO(NX,NQ)
      DIMENSION XCTOTFCHNLO(NX,NQ), XBTOTFCHNLO(NX,NQ)
      DIMENSION XGFCHNLO(NX,NQ)
      EQUIVALENCE(XUTOTFCHNLO,XUTOTFS(1,1,4,1))
      EQUIVALENCE(XDTOTFCHNLO,XDTOTFS(1,1,4,1))
      EQUIVALENCE(XSTOTFCHNLO,XSTOTFS(1,1,4,1))
      EQUIVALENCE(XUVALFCHNLO,XUVALFS(1,1,4,1))
      EQUIVALENCE(XDVALFCHNLO,XDVALFS(1,1,4,1))
      EQUIVALENCE(XSVALFCHNLO,XSVALFS(1,1,4,1))
      EQUIVALENCE(XCTOTFCHNLO,XCTOTFS(1,1,4,1))
      EQUIVALENCE(XBTOTFCHNLO,XBTOTFS(1,1,4,1))
      EQUIVALENCE(XGFCHNLO,XGFS(1,1,4,1))


*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 1.d0, 1.25D0, 1.5D0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
       DATA XB /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6        0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9 , 0.93, 1.0/
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.0.05D0) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
C          STOP
       ENDIF
       IF ( (Q2.LT.1.D0) .OR. (Q2.GT.1.D5) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
C          STOP
       ENDIF

*...INITIALIZATION :

*    SELECTION AND READING OF THE GRID :
      IF (FINI.NE.0) GOTO 16

*    SAVE IH AND IO, READ EVERYTHING IN, RESTORE IH AND IO
      IHSAVE = IH
      IOSAVE = IO

      DO 35 IH = 1, 4
      DO 36 IO = 0, 1

      IF ((IH.EQ.1).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='fragdata/dss/PINLO.GRID')
      ELSEIF ((IH.EQ.1).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='fragdata/dss/PILO.GRID')
      ELSEIF ((IH.EQ.2).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='fragdata/dss/KANLO.GRID')
       ELSEIF ((IH.EQ.2).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='fragdata/dss/KALO.GRID')
      ELSEIF ((IH.EQ.3).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='fragdata/dss/PRONLO.GRID')       
      ELSEIF ((IH.EQ.3).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='fragdata/dss/PROLO.GRID')       
      ELSEIF ((IH.EQ.4).and.(IO.EQ.1)) THEN
       IIREAD=11       
       OPEN(IIREAD,FILE='fragdata/dss/HNLO.GRID')       
      ELSEIF ((IH.EQ.4).and.(IO.EQ.0)) THEN
       IIREAD=12       
       OPEN(IIREAD,FILE='fragdata/dss/HLO.GRID')       
      ELSE
         WRITE(6,93)
 93      FORMAT (2X,' WRONG SET')
         STOP
      END IF
C
       DO 15 M = 1, NX-1 
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M),
     2                 PARTON(7,N,M), PARTON(8,N,M), PARTON(9,N,M)
  90   FORMAT (9(1PE10.3))
  15   CONTINUE
       CLOSE(IIREAD)

*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUTOTFS(IX,IQ,IH,IO) = PARTON(1,IQ,IX) / (XB1**4 * XB0**0.5)
        XDTOTFS(IX,IQ,IH,IO) = PARTON(2,IQ,IX) / (XB1**4 * XB0**0.5)
        XSTOTFS(IX,IQ,IH,IO) = PARTON(3,IQ,IX) / (XB1**4 * XB0**0.5) 
        XCTOTFS(IX,IQ,IH,IO) = PARTON(4,IQ,IX) / (XB1**7 * XB0**0.3) 
        XBTOTFS(IX,IQ,IH,IO) = PARTON(5,IQ,IX) / (XB1**7 * XB0**0.3)
        XGFS(IX,IQ,IH,IO)    = PARTON(6,IQ,IX) / (XB1**4 * XB0**0.3)
        XUVALFS(IX,IQ,IH,IO) = PARTON(7,IQ,IX) / (XB1**4 * XB0**0.5)
        XDVALFS(IX,IQ,IH,IO) = PARTON(8,IQ,IX) / (XB1**4 * XB0**0.5)
        XSVALFS(IX,IQ,IH,IO) = PARTON(9,IQ,IX) / (XB1**4 * XB0**0.5)
  20  CONTINUE
        XUTOTFS(NX,IQ,IH,IO) = 0.D0
        XDTOTFS(NX,IQ,IH,IO) = 0.D0
        XSTOTFS(NX,IQ,IH,IO) = 0.D0
        XCTOTFS(NX,IQ,IH,IO) = 0.D0
        XBTOTFS(NX,IQ,IH,IO) = 0.D0
        XGFS(NX,IQ,IH,IO)    = 0.D0
        XUVALFS(NX,IQ,IH,IO) = 0.D0
        XDVALFS(NX,IQ,IH,IO) = 0.D0
        XSVALFS(NX,IQ,IH,IO) = 0.D0
  10  CONTINUE  

 36   CONTINUE
 35   CONTINUE

      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE

      FINI = 1
      IH = IHSAVE
      IO = IOSAVE

* Start from here if tables have already been loaded

  16  CONTINUE

*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)

      Up = 0.D0
      UBp = 0.D0
      Dp = 0.D0
      DBp = 0.D0
      Sp = 0.D0
      SBp = 0.D0
      Cp = 0.D0
      Bp = 0.D0
      GL = 0.D0

* LO FF's

      IF(IO.EQ.0) THEN

        IF(IH.EQ.1) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT = FINT(NARG,XT,NA,ARRF,XUTOTFPILO)*(1.D0-X)**4 * X**0.5
            UVAL = FINT(NARG,XT,NA,ARRF,XUVALFPILO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT = FINT(NARG,XT,NA,ARRF,XDTOTFPILO)*(1.D0-X)**4 * X**0.5
            DVAL = FINT(NARG,XT,NA,ARRF,XDVALFPILO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT = FINT(NARG,XT,NA,ARRF,XSTOTFPILO)*(1.D0-X)**4 * X**0.5
            SVAL = FINT(NARG,XT,NA,ARRF,XSVALFPILO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT = FINT(NARG,XT,NA,ARRF,XCTOTFPILO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT = FINT(NARG,XT,NA,ARRF,XBTOTFPILO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFPILO)  * (1.D0-X)**4 * X**0.3
          END IF

        ELSEIF(IH.EQ.2) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT = FINT(NARG,XT,NA,ARRF,XUTOTFKALO)*(1.D0-X)**4 * X**0.5
            UVAL = FINT(NARG,XT,NA,ARRF,XUVALFKALO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT = FINT(NARG,XT,NA,ARRF,XDTOTFKALO)*(1.D0-X)**4 * X**0.5
            DVAL = FINT(NARG,XT,NA,ARRF,XDVALFKALO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT = FINT(NARG,XT,NA,ARRF,XSTOTFKALO)*(1.D0-X)**4 * X**0.5
            SVAL = FINT(NARG,XT,NA,ARRF,XSVALFKALO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT = FINT(NARG,XT,NA,ARRF,XCTOTFKALO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT = FINT(NARG,XT,NA,ARRF,XBTOTFKALO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFKALO)  * (1.D0-X)**4 * X**0.3
          END IF

        ELSEIF(IH.EQ.3) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT = FINT(NARG,XT,NA,ARRF,XUTOTFPRLO)*(1.D0-X)**4 * X**0.5
            UVAL = FINT(NARG,XT,NA,ARRF,XUVALFPRLO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT = FINT(NARG,XT,NA,ARRF,XDTOTFPRLO)*(1.D0-X)**4 * X**0.5
            DVAL = FINT(NARG,XT,NA,ARRF,XDVALFPRLO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT = FINT(NARG,XT,NA,ARRF,XSTOTFPRLO)*(1.D0-X)**4 * X**0.5
            SVAL = FINT(NARG,XT,NA,ARRF,XSVALFPRLO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT = FINT(NARG,XT,NA,ARRF,XCTOTFPRLO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT = FINT(NARG,XT,NA,ARRF,XBTOTFPRLO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFPRLO)  * (1.D0-X)**4 * X**0.3
          END IF
      
        ELSEIF(IH.EQ.4) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT = FINT(NARG,XT,NA,ARRF,XUTOTFCHLO)*(1.D0-X)**4 * X**0.5
            UVAL = FINT(NARG,XT,NA,ARRF,XUVALFCHLO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT = FINT(NARG,XT,NA,ARRF,XDTOTFCHLO)*(1.D0-X)**4 * X**0.5
            DVAL = FINT(NARG,XT,NA,ARRF,XDVALFCHLO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT = FINT(NARG,XT,NA,ARRF,XSTOTFCHLO)*(1.D0-X)**4 * X**0.5
            SVAL = FINT(NARG,XT,NA,ARRF,XSVALFCHLO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT = FINT(NARG,XT,NA,ARRF,XCTOTFCHLO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT = FINT(NARG,XT,NA,ARRF,XBTOTFCHLO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFCHLO)  * (1.D0-X)**4 * X**0.3
          END IF

        END IF

* NLO FF's

      ELSEIF(IO.EQ.1) THEN

        IF(IH.EQ.1) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT =FINT(NARG,XT,NA,ARRF,XUTOTFPINLO)*(1.D0-X)**4 * X**0.5
            UVAL =FINT(NARG,XT,NA,ARRF,XUVALFPINLO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT =FINT(NARG,XT,NA,ARRF,XDTOTFPINLO)*(1.D0-X)**4 * X**0.5
            DVAL =FINT(NARG,XT,NA,ARRF,XDVALFPINLO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT =FINT(NARG,XT,NA,ARRF,XSTOTFPINLO)*(1.D0-X)**4 * X**0.5
            SVAL =FINT(NARG,XT,NA,ARRF,XSVALFPINLO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT =FINT(NARG,XT,NA,ARRF,XCTOTFPINLO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT =FINT(NARG,XT,NA,ARRF,XBTOTFPINLO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFPINLO) * (1.D0-X)**4 * X**0.3
          END IF

        ELSEIF(IH.EQ.2) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT =FINT(NARG,XT,NA,ARRF,XUTOTFKANLO)*(1.D0-X)**4 * X**0.5
            UVAL =FINT(NARG,XT,NA,ARRF,XUVALFKANLO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT =FINT(NARG,XT,NA,ARRF,XDTOTFKANLO)*(1.D0-X)**4 * X**0.5
            DVAL =FINT(NARG,XT,NA,ARRF,XDVALFKANLO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT =FINT(NARG,XT,NA,ARRF,XSTOTFKANLO)*(1.D0-X)**4 * X**0.5
            SVAL =FINT(NARG,XT,NA,ARRF,XSVALFKANLO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT =FINT(NARG,XT,NA,ARRF,XCTOTFKANLO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT =FINT(NARG,XT,NA,ARRF,XBTOTFKANLO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFKANLO) * (1.D0-X)**4 * X**0.3
          END IF

        ELSEIF(IH.EQ.3) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT =FINT(NARG,XT,NA,ARRF,XUTOTFPRNLO)*(1.D0-X)**4 * X**0.5
            UVAL =FINT(NARG,XT,NA,ARRF,XUVALFPRNLO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT =FINT(NARG,XT,NA,ARRF,XDTOTFPRNLO)*(1.D0-X)**4 * X**0.5
            DVAL =FINT(NARG,XT,NA,ARRF,XDVALFPRNLO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT =FINT(NARG,XT,NA,ARRF,XSTOTFPRNLO)*(1.D0-X)**4 * X**0.5
            SVAL =FINT(NARG,XT,NA,ARRF,XSVALFPRNLO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT =FINT(NARG,XT,NA,ARRF,XCTOTFPRNLO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT =FINT(NARG,XT,NA,ARRF,XBTOTFPRNLO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFPRNLO) * (1.D0-X)**4 * X**0.3
          END IF

        ELSEIF(IH.EQ.4) THEN

          IF (ITYP.EQ.1 .OR. ITYP.EQ.2 ) THEN
            UTOT =FINT(NARG,XT,NA,ARRF,XUTOTFCHNLO)*(1.D0-X)**4 * X**0.5
            UVAL =FINT(NARG,XT,NA,ARRF,XUVALFCHNLO)*(1.D0-X)**4 * X**0.5
            Up  = (UTOT+UVAL)/2.
            UBp = (UTOT-UVAL)/2.
          ELSEIF (ITYP.EQ.3 .OR. ITYP.EQ.4 ) THEN
            DTOT =FINT(NARG,XT,NA,ARRF,XDTOTFCHNLO)*(1.D0-X)**4 * X**0.5
            DVAL =FINT(NARG,XT,NA,ARRF,XDVALFCHNLO)*(1.D0-X)**4 * X**0.5 
            Dp  = (DTOT+DVAL)/2.
            DBp = (DTOT-DVAL)/2.
          ELSEIF (ITYP.EQ.5 .OR. ITYP.EQ.6 ) THEN 
            STOT =FINT(NARG,XT,NA,ARRF,XSTOTFCHNLO)*(1.D0-X)**4 * X**0.5
            SVAL =FINT(NARG,XT,NA,ARRF,XSVALFCHNLO)*(1.D0-X)**4 * X**0.5
            Sp  = (STOT+SVAL)/2.
            SBp = (STOT-SVAL)/2.
          ELSEIF (ITYP.EQ.7 .OR. ITYP.EQ.8 ) THEN 
            CTOT =FINT(NARG,XT,NA,ARRF,XCTOTFCHNLO)*(1.D0-X)**7 * X**0.3
            Cp  =  CTOT/2.
          ELSEIF (ITYP.EQ.9 .OR. ITYP.EQ.10 ) THEN 
            BTOT =FINT(NARG,XT,NA,ARRF,XBTOTFCHNLO)*(1.D0-X)**7 * X**0.3
            Bp  =  BTOT/2.
          ELSEIF (ITYP.EQ.0 ) THEN
            GL   = FINT(NARG,XT,NA,ARRF,XGFCHNLO) * (1.D0-X)**4 * X**0.3
          END IF

        END IF
      
      END IF
              
      IF (IC.EQ.1) THEN
       U  = Up
       UB = UBp
       D  = Dp 
       DB = DBp
       S  = Sp
       SB = SBp
       C  = Cp
       B  = Bp
      ELSEIF (IC.EQ.-1) THEN
       U  = UBp
       UB = Up
       D  = DBp 
       DB = Dp
       S  = SBp
       SB = Sp
       C  = Cp
       B  = Bp
      ELSEIF (IC.EQ.0) THEN
       U  = (UBp+Up)/2.
       UB =  U
       D  = (DBp+Dp)/2. 
       DB =  D
       S  = (SBp+Sp)/2.
       SB =  S
       C  =  Cp
       B  =  Bp 
      ELSE
         WRITE(6,94)
 94      FORMAT (2X,' WRONG CHARGE')
         STOP
      END IF
      
      IF(ITYP.EQ.0) THEN
         DH = GL/X
      ELSEIF(ITYP.EQ.1) THEN
         DH = U/X
      ELSEIF(ITYP.EQ.2) THEN
         DH = UB/X
      ELSEIF(ITYP.EQ.3) THEN
         DH = D/X
      ELSEIF(ITYP.EQ.4) THEN
         DH = DB/X
      ELSEIF(ITYP.EQ.5) THEN
         DH = S/X
      ELSEIF(ITYP.EQ.6) THEN
         DH = SB/X
      ELSEIF(ITYP.EQ.7 .OR. ITYP.EQ.8) THEN
         DH = C/X
      ELSEIF(ITYP.EQ.9 .OR. ITYP.EQ.10) THEN
         DH = B/X
      ELSE
         DH = 0.D0
      END IF
 
      RETURN
      END


!> CERN LIBRARY ROUTINE E104 (INTERPOLATION)
!     Thomas Burton 29th Feb 2012
!     Changed:
!        ARG(5) to ARG(2)
!        NENT(5) to NENT(2)
!        ENT(63) to ENT(59)
!        TABLE(882) to TABLE(840)
!     to suppress multiple warnings like:
!        "Warning: Actual argument contains too few
!        elements for dummy argument..."
!     due to arrays being passed to FINT() having fewer elements.
!     The numbers of elements in the arrays passed from fDSS() now matches
!     the numbers in the FINT() arrays.
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(2),NENT(2),ENT(59),TABLE(840)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
      
