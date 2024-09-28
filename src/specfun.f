C       COMPUTATION OF SPECIAL FUNCTIONS
C
C          Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs.
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
C
C       Scipy changes:
C       - Compiled into a single source file and changed REAL To DBLE throughout.
C       - Changed according to ERRATA.
C       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
C       - Made functions return sf_error codes in ISFER variables instead
C         of printing warnings. The codes are
C         - SF_ERROR_OK        = 0: no error
C         - SF_ERROR_SINGULAR  = 1: singularity encountered
C         - SF_ERROR_UNDERFLOW = 2: floating point underflow
C         - SF_ERROR_OVERFLOW  = 3: floating point overflow
C         - SF_ERROR_SLOW      = 4: too many iterations required
C         - SF_ERROR_LOSS      = 5: loss of precision
C         - SF_ERROR_NO_RESULT = 6: no result obtained
C         - SF_ERROR_DOMAIN    = 7: out of domain
C         - SF_ERROR_ARG       = 8: invalid input parameter
C         - SF_ERROR_OTHER     = 9: unclassified error
C
        FUNCTION DNAN()
        DOUBLE PRECISION DNAN
        DNAN = 0.0D0
        DNAN = 0.0D0/DNAN
        END

        FUNCTION DINF()
        DOUBLE PRECISION DINF
        DINF = 1.0D300
        DINF = DINF*DINF
        END

C       **********************************

        SUBROUTINE PSI_SPEC(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=INT(XA)
           DO 10 K=1 ,N-1
            S=S+1.0D0/K
10         CONTINUE
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=INT(XA-.5)
           DO 20 K=1,N
            S=S+1.0/(2.0D0*K-1.0D0)
20         CONTINUE
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
               S=S+1.0D0/(XA+K)
30            CONTINUE
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
       
C       **********************************

        SUBROUTINE CHGUIT(A,B,X,HU,ID)
C
C       ======================================================
C       Purpose: Compute hypergeometric function U(a,b,x) by
C                using Gaussian-Legendre integration (n=60)
C       Input  : a  --- Parameter ( a > 0 )
C                b  --- Parameter
C                x  --- Argument ( x > 0 )
C       Output:  HU --- U(a,b,z)
C                ID --- Estimated number of significant digits
C       Routine called: GAMMA2 for computing Г(x)
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION T(30),W(30)
        DATA T/ .259597723012478D-01, .778093339495366D-01,
     &          .129449135396945D+00, .180739964873425D+00,
     &          .231543551376029D+00, .281722937423262D+00,
     &          .331142848268448D+00, .379670056576798D+00,
     &          .427173741583078D+00, .473525841761707D+00,
     &          .518601400058570D+00, .562278900753945D+00,
     &          .604440597048510D+00, .644972828489477D+00,
     &          .683766327381356D+00, .720716513355730D+00,
     &          .755723775306586D+00, .788693739932264D+00,
     &          .819537526162146D+00, .848171984785930D+00,
     &          .874519922646898D+00, .898510310810046D+00,
     &          .920078476177628D+00, .939166276116423D+00,
     &          .955722255839996D+00, .969701788765053D+00,
     &          .981067201752598D+00, .989787895222222D+00,
     &          .995840525118838D+00, .999210123227436D+00/
        DATA W/ .519078776312206D-01, .517679431749102D-01,
     &          .514884515009810D-01, .510701560698557D-01,
     &          .505141845325094D-01, .498220356905502D-01,
     &          .489955754557568D-01, .480370318199712D-01,
     &          .469489888489122D-01, .457343797161145D-01,
     &          .443964787957872D-01, .429388928359356D-01,
     &          .413655512355848D-01, .396806954523808D-01,
     &          .378888675692434D-01, .359948980510845D-01,
     &          .340038927249464D-01, .319212190192963D-01,
     &          .297524915007890D-01, .275035567499248D-01,
     &          .251804776215213D-01, .227895169439978D-01,
     &          .203371207294572D-01, .178299010142074D-01,
     &          .152746185967848D-01, .126781664768159D-01,
     &          .100475571822880D-01, .738993116334531D-02,
     &          .471272992695363D-02, .202681196887362D-02/
        ID=9
C       DLMF 13.4.4, integration up to C=12/X
        A1=A-1.0D0
        B1=B-A-1.0D0
        C=12.0D0/X
        HU0=0.0D0
        DO 20 M=10,100,5
           HU1=0.0D0
           G=0.5D0*C/M
           D=G
           DO 15 J=1,M
              S=0.0D0
              DO 10 K=1,30
                 T1=D+G*T(K)
                 T2=D-G*T(K)
                 F1=EXP(-X*T1)*T1**A1*(1.0D0+T1)**B1
                 F2=EXP(-X*T2)*T2**A1*(1.0D0+T2)**B1
                 S=S+W(K)*(F1+F2)
10            CONTINUE
              HU1=HU1+S*G
              D=D+2.0D0*G
15         CONTINUE
           IF (DABS(1.0D0-HU0/HU1).LT.1.0D-9) GO TO 25
           HU0=HU1
20      CONTINUE
25      CALL GAMMA2(A,GA)
        HU1=HU1/GA
C       DLMF 13.4.4 with substitution t=C/(1-u)
C       integration u from 0 to 1, i.e. t from C=12/X to infinity
        DO 40 M=2,10,2
           HU2=0.0D0
           G=0.5D0/M
           D=G
           DO 35 J=1,M
              S=0.0D0
              DO 30 K=1,30
                 T1=D+G*T(K)
                 T2=D-G*T(K)
                 T3=C/(1.0D0-T1)
                 T4=C/(1.0D0-T2)
                 F1=T3*T3/C*EXP(-X*T3)*T3**A1*(1.0D0+T3)**B1
                 F2=T4*T4/C*EXP(-X*T4)*T4**A1*(1.0D0+T4)**B1
                 S=S+W(K)*(F1+F2)
30            CONTINUE
              HU2=HU2+S*G
              D=D+2.0D0*G
35         CONTINUE
           IF (DABS(1.0D0-HU0/HU2).LT.1.0D-9) GO TO 45
           HU0=HU2
40      CONTINUE
45      CALL GAMMA2(A,GA)
        HU2=HU2/GA
        HU=HU1+HU2
        RETURN
        END




C       **********************************

        SUBROUTINE GAMMA2(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function Г(x)
C       Input :  x  --- Argument of Г(x)
C                       ( x is not equal to 0,-1,-2,…)
C       Output:  GA --- Г(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION G(26)
       PARAMETER (PI=3.141592653589793D0)
       DATA G/1.0D0,0.5772156649015329D0, -0.6558780715202538D0,
     &          -0.420026350340952D-1, 0.1665386113822915D0,
     &          -.421977345555443D-1, -.96219715278770D-2,
     &          .72189432466630D-2, -.11651675918591D-2, 
     &          -.2152416741149D-3, .1280502823882D-3, 
     &          -.201348547807D-4, -.12504934821D-5, 
     &          .11330272320D-5, -.2056338417D-6, .61160950D-8, 
     &          .50020075D-8, -.11812746D-8, .1043427D-9, 
     &          .77823D-11, -.36968D-11, .51D-12, 
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/     
       IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=INT(X-1)
              DO 10 K=2,M1
               GA=GA*K
10               END DO
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           R=1.0D0
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              DO 15 K=1,M
               R=R*(Z-K)
15            CONTINUE
              Z=Z-M
           ELSE
              Z=X
           ENDIF

           GR=G(26)
           DO 20 K=25,1,-1
            GR=GR*Z+G(K)
20         CONTINUE
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CHGU(A,B,X,HU,MD,ISFER)
C
C       =======================================================
C       Purpose: Compute the confluent hypergeometric function
C                U(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument  ( x > 0 )
C       Output:  HU --- U(a,b,x)
C                MD --- Method code
C                ISFER --- Error flag
C       Routines called:
C            (1) CHGUS for small x ( MD=1 )
C            (2) CHGUL for large x ( MD=2 )
C            (3) CHGUBI for integer b ( MD=3 )
C            (4) CHGUIT for numerical integration ( MD=4 )
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL IL1,IL2,IL3,BL1,BL2,BL3,BN
        AA=A-B+1.0D0
        ISFER=0
        IL1=A.EQ.INT(A).AND.A.LE.0.0
        IL2=AA.EQ.INT(AA).AND.AA.LE.0.0
        IL3=ABS(A*(A-B+1.0))/X.LE.2.0
        BL1=X.LE.5.0.OR.(X.LE.10.0.AND.A.LE.2.0)
        BL2=(X.GT.5.0.AND.X.LE.12.5).AND.(A.GE.1.0.AND.B.GE.A+4.0)
        BL3=X.GT.12.5.AND.A.GE.5.0.AND.B.GE.A+5.0
        BN=B.EQ.INT(B).AND.B.NE.0.0
        ID1=-100
        HU1=0.0D0
        IF (B.NE.INT(B)) THEN
           CALL CHGUS(A,B,X,HU,ID1)
           MD=1
           IF (ID1.GE.9) RETURN
           HU1=HU
        ENDIF
        IF (IL1.OR.IL2.OR.IL3) THEN
           CALL CHGUL(A,B,X,HU,ID)
           MD=2
           IF (ID.GE.9) RETURN
           IF (ID1.GT.ID) THEN
              MD=1
              ID=ID1
              HU=HU1
           ENDIF
        ENDIF
        IF (A.GE.1.0) THEN
           IF (BN.AND.(BL1.OR.BL2.OR.BL3)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
           ELSE
              CALL CHGUIT(A,B,X,HU,ID)
              MD=4
           ENDIF
        ELSE
           IF (B.LE.A) THEN
              A00=A
              B00=B
              A=A-B+1.0D0
              B=2.0D0-B
              CALL CHGUIT(A,B,X,HU,ID)
              HU=X**(1.0D0-B00)*HU
              A=A00
              B=B00
              MD=4
           ELSE IF (BN.AND.(.NOT.IL1)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
           ENDIF
        ENDIF
        IF (ID.LT.6) ISFER=6
        RETURN
        END

C       **********************************

        SUBROUTINE CHGM(A,B,X,HG)
C
C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,... )
C                x  --- Argument
C       Output:  HG --- M(a,b,x)
C       Routine called: CGAMA for computing complex ln[Г(x)]
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z)
        IMPLICIT COMPLEX(KIND=KIND(0.0D0)) (C)
        PI=3.141592653589793D0
        A0=A
        A1=A
        X0=X
        HG=0.0D0
C       DLMF 13.2.39
        IF (X.LT.0.0D0) THEN
           A=B-A
           A0=A
           X=DABS(X)
        ENDIF
        NL=0
        LA=0
        IF (A.GE.2.0D0) THEN
C       preparing terms for DLMF 13.3.1
           NL=1
           LA=INT(A)
           A=A-LA-1.0D0
        ENDIF
        Y0=0.0D0
        Y1=0.0D0
        DO 30 N=0,NL
           IF (A0.GE.2.0D0) A=A+1.0D0
           IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                 RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                 HG=HG+RG
                 IF (HG.NE.0D0.AND.DABS(RG/HG).LT.1.0D-15) THEN
C       DLMF 13.2.39 (cf. above)
                    IF (X0.LT.0.0D0) HG=HG*EXP(X0)
                    GO TO 25
                 ENDIF
15            CONTINUE
           ELSE
C       DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum
              Y=0.0D0
              CALL CGAMA(A,Y,0,TAR,TAI)
              CTA = DCMPLX(TAR, TAI)
              Y=0.0D0
              CALL CGAMA(B,Y,0,TBR,TBI)
              CTB = DCMPLX(TBR, TBI)
              XG=B-A
              Y=0.0D0
              CALL CGAMA(XG,Y,0,TBAR,TBAI)
              CTBA = DCMPLX(TBAR, TBAI)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8
                 R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                 SUM1=SUM1+R1
                 SUM2=SUM2+R2
20            END DO
              IF (X0.GE.0.0D0) THEN
                 HG1=DBLE(EXP(CTB-CTBA))*X**(-A)*DCOS(PI*A)*SUM1
                 HG2=DBLE(EXP(CTB-CTA+X))*X**(A-B)*SUM2
              ELSE
C       DLMF 13.2.39 (cf. above)
                 HG1=DBLE(EXP(CTB-CTBA+X0))*X**(-A)*DCOS(PI*A)*SUM1
                 HG2=DBLE(EXP(CTB-CTA))*X**(A-B)*SUM2
              ENDIF
              HG=HG1+HG2
           ENDIF
25         IF (N.EQ.0) Y0=HG
           IF (N.EQ.1) Y1=HG
30      CONTINUE
        IF (A0.GE.2.0D0) THEN
C       DLMF 13.3.1
           DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
              A=A+1.0D0
35          END DO
        ENDIF
        A=A1
        X=X0
        RETURN
        END

C       **********************************

        SUBROUTINE CCHG(A,B,Z,CHG)
C
C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,z) with real parameters a, b and a
C                complex argument z
C       Input :  a --- Parameter
C                b --- Parameter
C                z --- Complex argument
C       Output:  CHG --- M(a,b,z)
C       Routine called: CGAMA for computing complex ln[Г(x)]
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX(KIND=KIND(0.0D0)) (C,Z)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=A
        A1=A
        Z0=Z
        CY1=0
        CHW=0
        IF (B.EQ.0.0.OR.B.EQ.-INT(ABS(B))) THEN
           CHG=(1.0D+300,0.0D0)
        ELSE IF (A.EQ.0.0D0.OR.Z.EQ.0.0D0) THEN
           CHG=(1.0D0,0.0D0)
        ELSE IF (A.EQ.-1.0D0) THEN
           CHG=1.0D0-Z/B
        ELSE IF (A.EQ.B) THEN
           CHG=EXP(Z)
        ELSE IF (A-B.EQ.1.0D0) THEN
           CHG=(1.0D0+Z/B)*EXP(Z)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           CHG=(EXP(Z)-1.0D0)/Z
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           CR=(1.0D0,0.0D0)
           CHG=(1.0D0,0.0D0)
           DO 10 K=1,M
              CR=CR*(A+K-1.0D0)/K/(B+K-1.0D0)*Z
              CHG=CHG+CR
10         END DO
        ELSE
           X0=DBLE(Z)
           IF (X0.LT.0.0D0) THEN
              A=B-A
              A0=A
              Z=-Z
           ENDIF
           NL=0
           LA=0
           IF (A.GE.2.0D0) THEN
              NL=1
              LA=INT(A)
              A=A-LA-1.0D0
           ENDIF
           NS=0
           DO 30 N=0,NL
              IF (A0.GE.2.0D0) A=A+1.0D0
              IF (ABS(Z).LT.20.0D0+ABS(B).OR.A.LT.0.0D0) THEN
                 CHG=(1.0D0,0.0D0)
                 CRG=(1.0D0,0.0D0)
                 DO 15 J=1,500
                    CRG=CRG*(A+J-1.0D0)/(J*(B+J-1.0D0))*Z
                    CHG=CHG+CRG
                    IF (ABS((CHG-CHW)/CHG).LT.1.D-15) GO TO 25
                    CHW=CHG
15               CONTINUE
              ELSE
                 Y=0.0D0
                 CALL CGAMA(A,Y,0,G1R,G1I)
                 CG1 = DCMPLX(G1R, G1I)
                 Y=0.0D0
                 CALL CGAMA(B,Y,0,G2R,G2I)
                 CG2 = DCMPLX(G2R,G2I)
                 BA=B-A
                 Y=0.0D0
                 CALL CGAMA(BA,Y,0,G3R,G3I)
                 CG3 = DCMPLX(G3R, G3I)
                 CS1=(1.0D0,0.0D0)
                 CS2=(1.0D0,0.0D0)
                 CR1=(1.0D0,0.0D0)
                 CR2=(1.0D0,0.0D0)
                 DO 20 I=1,8
                    CR1=-CR1*(A+I-1.0D0)*(A-B+I)/(Z*I)
                    CR2=CR2*(B-A+I-1.0D0)*(I-A)/(Z*I)
                    CS1=CS1+CR1
                    CS2=CS2+CR2
20               END DO
                 X=DBLE(Z)
                 Y=DIMAG(Z)
                 IF (X.EQ.0.0.AND.Y.GE.0.0) THEN
                    PHI=0.5D0*PI
                 ELSE IF (X.EQ.0.0.AND.Y.LE.0.0) THEN
                    PHI=-0.5D0*PI
                 ELSE
                    PHI=DATAN(Y/X)
                 ENDIF
                 IF (PHI.GT.-0.5*PI.AND.PHI.LT.1.5*PI) NS=1
                 IF (PHI.GT.-1.5*PI.AND.PHI.LE.-0.5*PI) NS=-1
                 CFAC=EXP(NS*CI*PI*A)
                 IF (Y.EQ.0.0D0) CFAC=DCOS(PI*A)
                 CHG1=EXP(CG2-CG3)*Z**(-A)*CFAC*CS1
                 CHG2=EXP(CG2-CG1+Z)*Z**(A-B)*CS2
                 CHG=CHG1+CHG2
              ENDIF
25            IF (N.EQ.0) CY0=CHG
              IF (N.EQ.1) CY1=CHG
30         CONTINUE
           IF (A0.GE.2.0D0) THEN
              DO 35 I=1,LA-1
                 CHG=((2.0D0*A-B+Z)*CY1+(B-A)*CY0)/A
                 CY0=CY1
                 CY1=CHG
                 A=A+1.0D0
35             END DO
           ENDIF
           IF (X0.LT.0.0D0) CHG=CHG*EXP(-Z)
        ENDIF
        A=A1
        Z=Z0
        RETURN
        END

C       **********************************

        SUBROUTINE CGAMA(X,Y,KF,GR,GI)
C
C       =========================================================
C       Purpose: Compute the gamma function Г(z) or ln[Г(z)]
C                for a complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for ln[Г(z)]
C                       KF=1 for Г(z)
C       Output:  GR --- Real part of ln[Г(z)] or Г(z)
C                GI --- Imaginary part of ln[Г(z)] or Г(z)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        PARAMETER (PI=3.141592653589793D0)
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           GR=1.0D+300
           GI=0.0D0
           RETURN
        ELSE IF (X.LT.0.0D0) THEN
           X1=X
           Y1=Y
           X=-X
           Y=-Y
        ELSE
           Y1=0.0D0
           X1=X
        ENDIF
        X0=X
        NA=0
        IF (X.LE.7.0) THEN
           NA=INT(7-X)
           X0=X+NA
        ENDIF
        Z1=DSQRT(X0*X0+Y*Y)
        TH=DATAN(Y/X0)
        GR=(X0-.5D0)*DLOG(Z1)-TH*Y-X0+0.5D0*DLOG(2.0D0*PI)
        GI=TH*(X0-0.5D0)+Y*DLOG(Z1)-Y
        DO 10 K=1,10
           T=Z1**(1-2*K)
           GR=GR+A(K)*T*DCOS((2.0D0*K-1.0D0)*TH)
           GI=GI-A(K)*T*DSIN((2.0D0*K-1.0D0)*TH)
10      END DO
        IF (X.LE.7.0) THEN
           GR1=0.0D0
           GI1=0.0D0
           DO 15 J=0,NA-1
              GR1=GR1+.5D0*DLOG((X+J)**2+Y*Y)
              GI1=GI1+DATAN(Y/(X+J))
15         CONTINUE
           GR=GR-GR1
           GI=GI-GI1
        ENDIF
        IF (X1.LT.0.0D0) THEN
           Z1=DSQRT(X*X+Y*Y)
           TH1=DATAN(Y/X)
           SR=-DSIN(PI*X)*DCOSH(PI*Y)
           SI=-DCOS(PI*X)*DSINH(PI*Y)
           Z2=DSQRT(SR*SR+SI*SI)
           TH2=DATAN(SI/SR)
           IF (SR.LT.0.0D0) TH2=PI+TH2
           GR=DLOG(PI/(Z1*Z2))-GR
           GI=-TH1-TH2-GI
           X=X1
           Y=Y1
        ENDIF
        IF (KF.EQ.1) THEN
           G0=EXP(GR)
           GR=G0*DCOS(GI)
           GI=G0*DSIN(GI)
        ENDIF
        RETURN
        END

     
C       **********************************

        SUBROUTINE CHGUS(A,B,X,HU,ID)
C
C       ======================================================
C       Purpose: Compute confluent hypergeometric function
C                U(a,b,x) for small argument x
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,...)
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       Routine called: GAMMA2 for computing gamma function
C       ======================================================
C
C       DLMF 13.2.42 with prefactors rewritten according to
C       DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ID=-100
        PI=3.141592653589793D0
        CALL GAMMA2(A,GA)
        CALL GAMMA2(B,GB)
        XG1=1.0D0+A-B
        CALL GAMMA2(XG1,GAB)
        XG2=2.0D0-B
        CALL GAMMA2(XG2,GB2)
        HU0=PI/DSIN(PI*B)
        R1=HU0/(GAB*GB)
        R2=HU0*X**(1.0D0-B)/(GA*GB2)
        HU=R1-R2
        HMAX=0.0D0
        HMIN=1.0D+300
        H0=0.0D0
        DO 10 J=1,150
           R1=R1*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
           R2=R2*(A-B+J)/(J*(1.0D0-B+J))*X
           HU=HU+R1-R2
           HUA=DABS(HU)
           IF (HUA.GT.HMAX) HMAX=HUA
           IF (HUA.LT.HMIN) HMIN=HUA
           IF (DABS(HU-H0).LT.DABS(HU)*1.0D-15) GO TO 15
           H0=HU
10      END DO
15      D1=LOG10(HMAX)
        D2=0.0D0
        IF (HMIN.NE.0.0) D2=LOG10(HMIN)
        ID=INT(15-ABS(D1-D2))
        RETURN
        END

C       **********************************

C       **********************************

        SUBROUTINE CHGUL(A,B,X,HU,ID)
C
C       =======================================================
C       Purpose: Compute the confluent hypergeometric function
C                U(a,b,x) for large argument x
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL IL1,IL2
        ID=-100
        AA=A-B+1.0D0
        IL1=A.EQ.INT(A).AND.A.LE.0.0
        IL2=AA.EQ.INT(AA).AND.AA.LE.0.0
        NM=0
        IF (IL1) NM=INT(ABS(A))
        IF (IL2) NM=INT(ABS(AA))
C       IL1: DLMF 13.2.7 with k=-s-a
C       IL2: DLMF 13.2.8
        IF (IL1.OR.IL2) THEN
           HU=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              HU=HU+R
10         CONTINUE
           HU=X**(-A)*HU
           ID=10
        ELSE
C       DLMF 13.7.3
           HU=1.0D0
           R=1.0D0
           DO 15 K=1,25
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              RA=DABS(R)
              IF (K.GT.5.AND.RA.GE.R0.OR.RA.LT.1.0D-15) GO TO 20
              R0=RA
              HU=HU+R
15         END DO
20         ID=INT(ABS(LOG10(RA)))
           HU=X**(-A)*HU
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CHGUBI(A,B,X,HU,ID)
C
C       ======================================================
C       Purpose: Compute confluent hypergeometric function
C                U(a,b,x) with integer b ( b = ±1,±2,... )
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       Routines called:
C            (1) GAMMA2 for computing gamma function Г(x)
C            (2) PSI_SPEC for computing psi function
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ID=-100
        EL=0.5772156649015329D0
        N=INT(ABS(B-1))
        RN1=1.0D0
        RN=1.0D0
        DO 10 J=1,N
           RN=RN*J
           IF (J.EQ.N-1) RN1=RN
10      CONTINUE
        CALL PSI_SPEC(A,PS)
        CALL GAMMA2(A,GA)
        IF (B.GT.0.0) THEN
           A0=A
           A1=A-N
           A2=A1
           CALL GAMMA2(A1,GA1)
           UA=(-1)**(N-1)/(RN*GA1)
           UB=RN1/GA*X**(-N)
        ELSE
           A0=A+N
           A1=A0
           A2=A
           CALL GAMMA2(A1,GA1)
           UA=(-1)**(N-1)/(RN*GA)*X**N
           UB=RN1/GA1
        ENDIF
        HM1=1.0D0
        R=1.0D0
        HMAX=0.0D0
        HMIN=1.0D+300
        H0=0D0
        DO 15 K=1,150
           R=R*(A0+K-1.0D0)*X/((N+K)*K)
           HM1=HM1+R
           HU1=DABS(HM1)
           IF (HU1.GT.HMAX) HMAX=HU1
           IF (HU1.LT.HMIN) HMIN=HU1
           IF (DABS(HM1-H0).LT.DABS(HM1)*1.0D-15) GO TO 20
           H0=HM1
15      END DO
20      DA1=LOG10(HMAX)
        DA2=0.0D0
        IF (HMIN.NE.0.0) DA2=LOG10(HMIN)
        ID=INT(15-ABS(DA1-DA2))
        HM1=HM1*DLOG(X)
        S0=0.0D0
        DO 25 M=1,N
           IF (B.GE.0.0) S0=S0-1.0D0/M
           IF (B.LT.0.0) S0=S0+(1.0D0-A)/(M*(A+M-1.0D0))
25      END DO
        HM2=PS+2.0D0*EL+S0
        R=1.0D0
        HMAX=0.0D0
        HMIN=1.0D+300
        DO 50 K=1,150
           S1=0.0D0
           S2=0.0D0
           IF (B.GT.0.0) THEN
              DO 30 M=1,K
                 S1=S1-(M+2.0D0*A-2.0D0)/(M*(M+A-1.0D0))
30            END DO
              DO 35 M=1,N
                 S2=S2+1.0D0/(K+M)
35            END DO
           ELSE
              DO 40 M=1,K+N
                 S1=S1+(1.0D0-A)/(M*(M+A-1.0D0))
40            END DO
              DO 45 M=1,K
                 S2=S2+1.0D0/M
45            END DO
           ENDIF
           HW=2.0D0*EL+PS+S1-S2
           R=R*(A0+K-1.0D0)*X/((N+K)*K)
           HM2=HM2+R*HW
           HU2=DABS(HM2)
           IF (HU2.GT.HMAX) HMAX=HU2
           IF (HU2.LT.HMIN) HMIN=HU2
           IF (DABS((HM2-H0)/HM2).LT.1.0D-15) GO TO 55
           H0=HM2
50      END DO
55      DB1=LOG10(HMAX)
        DB2=0.0D0
        IF (HMIN.NE.0.0) DB2=LOG10(HMIN)
        ID1=INT(15-ABS(DB1-DB2))
        IF (ID1.LT.ID) ID=ID1
        HM3=1.0D0
        IF (N.EQ.0) HM3=0.0D0
        R=1.0D0
        DO 60 K=1,N-1
           R=R*(A2+K-1.0D0)/((K-N)*K)*X
           HM3=HM3+R
60      END DO
        SA=UA*(HM1+HM2)
        SB=UB*HM3
        HU=SA+SB
        ID2=0.0D0
        IF (SA.NE.0.0) ID1=INT(LOG10(ABS(SA)))
        IF (HU.NE.0.0) ID2=INT(LOG10(ABS(HU)))
        IF (SA*SB.LT.0.0) ID=ID-ABS(ID1-ID2)
        RETURN
        END
C       **********************************

        
