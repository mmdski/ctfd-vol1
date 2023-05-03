
C
C     FIVOL applies the finite volume method to the solution of
C     Laplace's equation in Cartesian coordinates on a polar grid.
C     The discretized equation is solve by SOR

      DIMENSION X(21,21),Y(21,21),QAB(21,21),PAB(21,21),QBC(21,21),
     1PBC(21,21),QCD(21,21),PCD(21,21),QDA(21,21),PDA(21,21),
     2PHI(21,21),PHIX(21,21)
C
      OPEN(1,FILE='FIVOL.DAT')
      OPEN(6,FILE='FIVOL.OUT')
      READ(1,1)JMAX,KMAX,NMAX
      READ(1,2)RW,RX,RY,RZ,THEB,THEN,EPS,OM
    1 FORMAT(8I5)
    2 FORMAT(8E10.3)
C
      WRITE(6,3)
      WRITE(6,4)JMAX,KMAX,NMAX,EPS,OM
      WRITE(6,5)RW,RX,RY,RZ,THEB,THEN
    3 FORMAT(' LAPLACE EQUATION BY FINTE VOLUME METHOD',//)
    4 FORMAT(' JMAX=',I2,' KMAX=',I2,'   NMAX=',I5,
     15X,' EPS=',E10.3,'  OM=',F5.3)
    5 FORMAT('   RW=',F5.3,'  RX=',F5.3,'   RY=',F5.3,'  RZ=',F5.3,
     15X,'  THEB=',F5.1,'  THEN=',F5.1,//)
C
      JMAP = JMAX - 1
      KMAP = KMAX - 1
      AJM = JMAP
      AKM = KMAP
      DRWX = (RX - RW)/AJM
      DRZY = (RY - RZ)/AKM
      DTH = (THEN-THEB)/AKM
      PI = 3.1415927
C
C     SET X, Y, exact and initial PHI
C
      DO 7 K = 1,KMAX
      AK = K - 1
      THK = (THEB + AK*DTH)*PI/180.
      CK = COS(THK)
      SK = SIN(THK)
      DR = DRWX + (DRZY - DRWX)*AK/AKM
      RWZ = RW + (RZ - RW)*AK/AKM
      DO 6 J = 1,JMAX
      AJ = J - 1
      R = RWZ + AJ*DR
      X(J,K) = R*CK
      Y(J,K) = R*SK
      PHIX(J,K) = SK/R
      PHI(J,K) = PHIX(J,K)
    6 CONTINUE
    7 CONTINUE
C
C     SET BOUNDARY VALUES OF PHI
C
      DO 8 J = 1,JMAX
      PHI(J,1) = 0.
      PHI(J,KMAX) = PHIX(J,KMAX)
    8 CONTINUE
      DO 9 K = 1,KMAX
      PHI(1,K) = PHIX(1,K)
      PHI(JMAX,K) = PHIX(JMAX,K)
    9 CONTINUE
C
C    SET GRID RELATED PARAMETERS
C
      DO 11 K = 2,KMAP
      KM = K - 1
      KP = K + 1
      DO 10 J = 2,JMAP
      JM = J - 1
      JP = J + 1
      XA = 0.25*(X(J,K) + X(JM,K) + X(JM,KM) + X(J,KM))
      YA = 0.25*(Y(J,K) + Y(JM,K) + Y(JM,KM) + Y(J,KM))
      XB = 0.25*(X(J,K) + X(J,KM) + X(JP,KM) + X(JP,K))
      YB = 0.25*(Y(J,K) + Y(J,KM) + Y(JP,KM) + Y(JP,K))
      XC = 0.25*(X(J,K) + X(JP,K) + X(JP,KP) + X(J,KP))
      YC = 0.25*(Y(J,K) + Y(JP,K) + Y(JP,KP) + Y(J,KP))
      XD = 0.25*(X(J,K) + X(J,KP) + X(JM,KP) + X(JM,K))
      YD = 0.25*(Y(J,K) + Y(J,KP) + Y(JM,KP) + Y(JM,K))
C
C     SIDE AB
C
      DXA = XB - XA
      DYA = YB - YA
      DXK = X(J,K) - X(J,K-1)
      DYK = Y(J,K) - Y(J,K-1)
      SAB = ABS(DXA*DYK - DXK*DYA)
      QAB(J,K) = (DXA*DXA + DYA*DYA)/SAB
      PAB(J,K) = (DXA*DXK + DYA*DYK)/SAB
C
C     SIDE BC
C
      DXB = XC - XB
      DYB = YC - YB
      DXJ = X(J,K) - X(J+1,K)
      DYJ = Y(J,K) - Y(J+1,K)
      SBC = ABS(DYJ*DXB - DXJ*DYB)
      QBC(J,K) = (DXB*DXB + DYB*DYB)/SBC
      PBC(J,K) = (DXB*DXJ + DYB*DYJ)/SBC
C
C     SIDE CD
C
      DXC = XD - XC
      DYC = YD - YC
      DXK = X(J,K) - X(J,K+1)
      DYK = Y(J,K) - Y(J,K+1)
      SCD = ABS(DXC*DYK - DYC*DXK)
      QCD(J,K) = (DXC*DXC + DYC*DYC)/SCD
      PCD(J,K) = (DXC*DXK + DYC*DYK)/SCD
C
C     SIDE DA
C
      DXD = XA - XD
      DYD = YA - YD
      DXJ = X(J,K) - X(J-1,K)
      DYJ = Y(J,K) - Y(J-1,K)
      SDA = ABS(DXJ*DYD - DYJ*DXD)
      QDA(J,K) = (DXD*DXD + DYD*DYD)/SDA
      PDA(J,K) = (DXD*DXJ + DYD*DYJ)/SDA
   10 CONTINUE
   11 CONTINUE
C
C     Iterate using SOR
C
      DO 14 N = 1,NMAX
      SUM = 0.
      DO 13 K = 2,KMAP
      KM = K - 1
      KP = K + 1
      DO 12 J = 2,JMAP
      JM = J - 1
      JP = J + 1
      PHD = 0.25*(PCD(J,K)-PDA(J,K))*PHI(JM,KP)
      PHD = PHD + (QCD(J,K) + 0.25*(PBC(J,K)-PDA(J,K)))*PHI(J,KP)
      PHD = PHD + 0.25*(PBC(J,K)-PCD(J,K))*PHI(JP,KP)
      PHD = PHD + (QDA(J,K) + 0.25*(PCD(J,K)-PAB(J,K)))*PHI(JM,K)
      PHD = PHD + (QBC(J,K) + 0.25*(PAB(J,K)-PCD(J,K)))*PHI(JP,K)
      PHD = PHD + 0.25*(PDA(J,K) - PAB(J,K))*PHI(JM,KM)
      PHD = PHD + (QAB(J,K) + 0.25*(PDA(J,K)-PBC(J,K)))*PHI(J,KM)
      PHD = PHD + 0.25*(PAB(J,K) - PBC(J,K))*PHI(JP,KM)
      PHD = PHD/(QAB(J,K)+QBC(J,K)+QCD(J,K)+QDA(J,K))
      DIF = PHD - PHI(J,K)
      SUM = SUM + DIF*DIF
      PHI(J,K) = PHI(J,K) + OM*DIF
   12 CONTINUE
   13 CONTINUE
      RMS = SQRT(SUM/(AJM-1.)/(AKM-1.))
      IF(RMS .LT. EPS)GOTO 16
   14 CONTINUE
      WRITE(6,15)NMAX,RMS
   15 FORMAT(' CONVERGENCE NOT ACHIEVED IN',I5,' STEPS',5X,' RMS=',
     1E12.5)
C
C    Compare solution with exact
C
   16 SUM = 0.
      DO 21 K = 1,KMAX
      WRITE(6,17)K
   17 FORMAT(/,'  K=',I2)
      DO 18 J = 1,JMAX
      DIF = PHI(J,K) - PHIX(J,K)
      SUM = SUM + DIF*DIF
   18 CONTINUE
      WRITE(6,19)(PHI(J,K),J=1,JMAX)
      WRITE(6,20)(PHIX(J,K),J=1,JMAX)
   19 FORMAT('  PHI=',10F7.4)
   20 FORMAT('  PHX=',10F7.4)
   21 CONTINUE
      RMS = SQRT(SUM/(AJM-1.)/(AKM-1.))
      WRITE(6,22)N,RMS
   22 FORMAT(/,' CONVERGED AFTER ',I3,'  STEPS',4X,'  RMS=',E12.5)
      STOP
      END
