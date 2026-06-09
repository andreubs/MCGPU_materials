CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  penelope_abridged.f                                                  C
C  Photon-physics subset of PENELOPE 2006 for MCGPU_material_generator. C
C  Contains only the subroutines needed for Rayleigh, Compton,          C
C  photoelectric and pair-production MFP table generation.              C
C
C    EGRID ------------------> PENELOPE log energy grid (50 eV – 1 GeV, 200 pts)
C    PHMFP ------------------> Photon mean free path (Rayleigh/Compton/Photo/Pair)
C    GRAaI, GRAaT, GRAaD, GRAaD1 --> Rayleigh cross section + RITA grid
C    GCOaT, GCOaD -----------> Compton cross section (impulse approximation)
C    GPHa0 ------------------> Initialise photoelectric tables to zero
C    GPPa0 ------------------> Initialise pair-production sampling parameters
C    BLOCK DATA PENDAT ------> Hardcoded element data: ATW, EPX, RA1–RA5, RSCR, ETA
C   SPLINE, FINDI  ----------> Spline fit and binary search
C    SUMGA, RITAI0, SIMPSU --> Numerical integration (Gaussian, RITA, Simpson)  
C                                                                       C
C  PENELOPE/PENGEOM (version 2006). Functions extracted from penelope.f C
C  Copyright (c) 2001-2006 Universitat de Barcelona                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  *********************************************************************
C                       SUBROUTINE EGRID
C  *********************************************************************
      SUBROUTINE EGRID(EMIN,EMAX)
C
C  This subroutine sets the energy grid where transport functions are
C  tabulated. The grid is logarithmically spaced and we assume that it
C  is dense enough to permit accurate linear log-log interpolation of
C  the tabulated functions.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C
C  ****  Consistency of the interval end-points.
C
      IF(EMIN.LT.50.0D0) EMIN=50.0D0
      IF(EMIN.GT.EMAX-1.0D0) THEN
        WRITE(26,2100) EMIN,EMAX
 2100   FORMAT(/3X,'EMIN =',1P,E11.4,' eV,  EMAX =',E11.4,' eV')
        STOP 'EGRID. The energy interval is too narrow.'
      ENDIF
C
C  ****  Energy grid points.
C
      EL=0.99999D0*EMIN
      EU=1.00001D0*EMAX
      DLFC=LOG(EU/EL)/DBLE(NEGP-1)
      DLEMP1=LOG(EL)
      DLEMP(1)=DLEMP1
      ET(1)=EL
      DO I=2,NEGP
        DLEMP(I)=DLEMP(I-1)+DLFC
        ET(I)=EXP(DLEMP(I))
      ENDDO
      DLFC=1.0D0/DLFC
C
C  NOTE: To determine the interval KE where the energy E is located, we
C  do the following,
C     XEL=LOG(E)
C     XE=1.0D0+(XEL-DLEMP1)*DLFC
C     KE=XE
C     XEK=XE-KE  ! 'fractional' part of XE (used for interpolation).
C
      RETURN
      END

C  *********************************************************************
C                  FUNCTION PHMFP
C  *********************************************************************
      FUNCTION PHMFP(E,KPAR,M,ICOL)
C
C  This function computes the mean free path (in cm) of particles of
C  type KPAR and energy E between hard interactions of kind ICOL in
C  material M. If ICOL does not correspond to a hard interaction type,
C  the result is set equal to 1.0D16.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  DW1EL(MAXMAT,NEGP),DW2EL(MAXMAT,NEGP),
     5  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     6  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  DW1PL(MAXMAT,NEGP),DW2PL(MAXMAT,NEGP),
     5  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     6  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
C
      IF(E.LE.EL.OR.E.GE.EU) THEN
        PHMFP=1.0D50
        RETURN
      ENDIF
      XEL=LOG(E)
      XE=1.0D0+(XEL-DLEMP1)*DLFC
      KE=XE
      IF(KE.LT.1) KE=1
      IF(KE.GE.NEGP) KE=NEGP-1
      XEK=XE-KE
C
      HMFP=1.0D-35
      IF(KPAR.EQ.1) THEN
        IF(ICOL.EQ.2) THEN
          HMFP=EXP(SEHEL(M,KE)+(SEHEL(M,KE+1)-SEHEL(M,KE))*XEK)
        ELSE IF(ICOL.EQ.3) THEN
          HMFP=EXP(SEHIN(M,KE)+(SEHIN(M,KE+1)-SEHIN(M,KE))*XEK)
        ELSE IF(ICOL.EQ.4) THEN
          HMFP=EXP(SEHBR(M,KE)+(SEHBR(M,KE+1)-SEHBR(M,KE))*XEK)
        ELSE IF(ICOL.EQ.5) THEN
          HMFP=EXP(SEISI(M,KE)+(SEISI(M,KE+1)-SEISI(M,KE))*XEK)
        ELSE IF(ICOL.EQ.8) THEN
          HMFP=EXP(SEAUX(M,KE)+(SEAUX(M,KE+1)-SEAUX(M,KE))*XEK)
        ENDIF
      ELSE IF(KPAR.EQ.2) THEN
        IF(ICOL.EQ.1) THEN
          HMFP=EXP(SGRA(M,KE)+(SGRA(M,KE+1)-SGRA(M,KE))*XEK)
        ELSE IF(ICOL.EQ.2) THEN
          HMFP=EXP(SGCO(M,KE)+(SGCO(M,KE+1)-SGCO(M,KE))*XEK)
        ELSE IF(ICOL.EQ.3) THEN
          PTOT=0.0D0
          DO IEL=1,NELEM(M)
            IZZ=IZ(M,IEL)
C  ****  Binary search.
            I=IPHF(IZZ)
            IU=IPHL(IZZ)
    1       IT=(I+IU)/2
            IF(XEL.GT.EPH(IT)) THEN
              I=IT
            ELSE
              IU=IT
            ENDIF
            IF(IU-I.GT.1) GO TO 1
C
            DEE=EPH(I+1)-EPH(I)
            IF(DEE.GT.1.0D-15) THEN
              PCSL=XPH(I,1)+(XPH(I+1,1)-XPH(I,1))*(XEL-EPH(I))/DEE
            ELSE
              PCSL=XPH(I,1)
            ENDIF
            PTOT=PTOT+STF(M,IEL)*EXP(PCSL)
          ENDDO
          HMFP=PTOT*VMOL(M)
        ELSE IF(ICOL.EQ.4.AND.E.GT.1.022D6) THEN
          HMFP=EXP(SGPP(M,KE)+(SGPP(M,KE+1)-SGPP(M,KE))*XEK)
        ELSE IF(ICOL.EQ.6) THEN
          HMFP=EXP(SGAUX(M,KE)+(SGAUX(M,KE+1)-SGAUX(M,KE))*XEK)
        ENDIF
      ELSE IF(KPAR.EQ.3) THEN
        IF(ICOL.EQ.2) THEN
          HMFP=EXP(SPHEL(M,KE)+(SPHEL(M,KE+1)-SPHEL(M,KE))*XEK)
        ELSE IF(ICOL.EQ.3) THEN
          HMFP=EXP(SPHIN(M,KE)+(SPHIN(M,KE+1)-SPHIN(M,KE))*XEK)
        ELSE IF(ICOL.EQ.4) THEN
          HMFP=EXP(SPHBR(M,KE)+(SPHBR(M,KE+1)-SPHBR(M,KE))*XEK)
        ELSE IF(ICOL.EQ.5) THEN
          HMFP=EXP(SPISI(M,KE)+(SPISI(M,KE+1)-SPISI(M,KE))*XEK)
        ELSE IF(ICOL.EQ.6) THEN
          HMFP=EXP(SPAN(M,KE)+(SPAN(M,KE+1)-SPAN(M,KE))*XEK)
        ELSE IF(ICOL.EQ.8) THEN
          HMFP=EXP(SPAUX(M,KE)+(SPAUX(M,KE+1)-SPAUX(M,KE))*XEK)
        ENDIF
      ENDIF
C
      PHMFP=1.0D0/MAX(HMFP,1.0D-35)
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE GRAaI
C  *********************************************************************
      SUBROUTINE GRAaI(M)
C
C  Initialisation of random sampling for Rayleigh (coherent) scattering
C  of photons.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
      PARAMETER (RREV=1.0D0/REV)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C
      PARAMETER (NM=512)
      COMMON/CRITA/XTI(NM),PACI(NM),AI(NM),BI(NM),NPI,
     1             ITTLI(NM),ITTUI(NM),NPM1I
C
      PARAMETER (NP=128)
      COMMON/CGRA/XCO(NP,MAXMAT),PCO(NP,MAXMAT),ACO(NP,MAXMAT),
     1  BCO(NP,MAXMAT),PMAX(NEGP,MAXMAT),ITLCO(NP,MAXMAT),
     2  ITUCO(NP,MAXMAT)
C
      PARAMETER (NIP=51)
      DIMENSION XI(NIP),FUN(NIP),SUM(NIP)
      COMMON/CGRA00/FACTE,X2MAX,MM,MOM
      EXTERNAL GRAaD1
C
      IZZ=0
      DO I=1,NELEM(M)
        IZZ=MAX(IZZ,IZ(M,I))
      ENDDO
C
      MM=M
      X2MIN=0.0D0
      X2MAX=4.0D0*20.6074D0**2*(200.0D0*IZZ)**2
      NPT=NP
      NU=NPT/4
      CALL RITAI0(GRAaD1,X2MIN,X2MAX,NPT,NU,ERRM,0)
      IF(NPI.NE.NP) THEN
        WRITE(26,*) 'RITA initialisation error in GRAaI.'
        WRITE(26,*) 'The number of fixed grid points is ',NPI
        WRITE(26,*) 'The required number of grid points was ',NP
        STOP 'RITA initialisation error in GRAaI.'
      ENDIF
      IF(ERRM.GT.1.0D-5) THEN
        WRITE(26,*) 'RITA interpolation error is too large in GRAaI.'
        WRITE(26,*) 'The interpolation error is ',ERRM
        STOP 'RITA interpolation error is too large in GRAaI.'
      ENDIF
C
C  ****  Upper limit of the X2 interval for the PENELOPE grid energies.
C
      DO IE=1,NEGP
        XM=2.0D0*20.6074D0*ET(IE)*RREV
        X2M=XM*XM
        IF(X2M.GT.XTI(1)) THEN
          IF(X2M.LT.XTI(NP)) THEN
            I=1
            J=NPI
    1       IT=(I+J)/2
            IF(X2M.GT.XTI(IT)) THEN
              I=IT
            ELSE
              J=IT
            ENDIF
            IF(J-I.GT.1) GO TO 1
C
            X1=XTI(I)
            X2=X2M
            DX=(X2-X1)/DBLE(NIP-1)
            DO K=1,NIP
              XI(K)=X1+DBLE(K-1)*DX
              TAU=(XI(K)-XTI(I))/(XTI(I+1)-XTI(I))
              CON1=2.0D0*BI(I)*TAU
              CI=1.0D0+AI(I)+BI(I)
              CON2=CI-AI(I)*TAU
              IF(ABS(CON1).GT.1.0D-16*ABS(CON2)) THEN
                ETAP=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))
     1              /CON1
              ELSE
                ETAP=TAU/CON2
              ENDIF
              FUN(K)=(PACI(I+1)-PACI(I))
     1              *(1.0D0+(AI(I)+BI(I)*ETAP)*ETAP)**2
     2              /((1.0D0-BI(I)*ETAP*ETAP)*CI*(XTI(I+1)-XTI(I)))
            ENDDO
            CALL SIMPSU(DX,FUN,SUM,NIP)
            PMAX(IE,M)=PACI(I)+SUM(NIP)
          ELSE
            PMAX(IE,M)=1.0D0
          ENDIF
        ELSE
          PMAX(IE,M)=PACI(1)
        ENDIF
      ENDDO
C
      DO I=1,NP
        XCO(I,M)=XTI(I)
        PCO(I,M)=PACI(I)
        ACO(I,M)=AI(I)
        BCO(I,M)=BI(I)
        ITLCO(I,M)=ITTLI(I)
        ITUCO(I,M)=ITTUI(I)
      ENDDO
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE GRAaT
C  *********************************************************************
      SUBROUTINE GRAaT(E,CS,M)
C
C  Total cross section for Rayleigh (coherent) photon scattering. Born
C  approximation with analytical atomic form factors.
C
C  Input arguments:
C    E ........ photon energy (eV).
C    M ........ material where photons propagate.
C  Output argument:
C    CS ....... coherent total cross section (cm**2/molecule).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
      PARAMETER (RREV=1.0D0/REV)
      PARAMETER (ELRAD=2.817940325D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CGRA00/FACTE,X2MAX,MM,MOM
C
      EXTERNAL GRAaD
C
      MM=M
      MOM=0
      IZZ=0
      DO I=1,NELEM(M)
        IZZ=MAX(IZZ,IZ(M,I))
      ENDDO
      EC=MIN(E,REV*IZZ)
      FACTE=2.0D0*20.6074D0**2*(EC*RREV)**2
      X2MAX=2.0D0*FACTE
      CS=SUMGA(GRAaD,-1.0D0,1.0D0,1.0D-6)
      CS=PIELR2*CS*(EC/E)**2
      RETURN
      END

C  *********************************************************************
C                       FUNCTION GRAaD
C  *********************************************************************
      FUNCTION GRAaD(CDT)
C
C  Differential x-section for Rayleigh scattering.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
C
      COMMON/CGRA00/FACTE,X2MAX,M,MOM
C
      X2=FACTE*(1.0D0-CDT)
      IF(X2.GT.X2MAX) THEN
        GRAaD=0.0D0
        RETURN
      ENDIF
      GRAaD=(1.0D0+CDT*CDT)*GRAaD1(X2)
      IF(MOM.GT.0) GRAaD=GRAaD*((1.0D0-CDT)*0.5D0)**MOM
      RETURN
      END

C  *********************************************************************
C                       FUNCTION GRAaD1
C  *********************************************************************
      FUNCTION GRAaD1(X2)
C
C  Squared molecular form factor (additivity rule).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (SL=137.03599911D0)  ! Speed of light (1/alpha)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
      COMMON/CGRA00/FACTE,X2MAX,M,MOM
C
      GRAaD1=0.0D0
      IF(X2.GT.X2MAX) RETURN
      X=SQRT(X2)
      DO I=1,NELEM(M)
C  ****  Atomic form factors.
        IZZ=IZ(M,I)
        FA=IZZ*(1.0D0+X2*(RA1(IZZ)+X*(RA2(IZZ)+X*RA3(IZZ))))
     1    /(1.0D0+X2*(RA4(IZZ)+X2*RA5(IZZ)))**2
        IF(IZZ.GT.10.AND.FA.LT.2.0D0) THEN
          PA=(IZZ-0.3125D0)/SL
          PG=SQRT(1.0D0-PA*PA)
          PQ=2.426311D-2*X/PA
          FB=SIN(2.0D0*PG*ATAN2(PQ,1.0D0))/(PG*PQ*(1.0D0+PQ*PQ)**PG)
          FA=MAX(FA,FB)
        ENDIF
        GRAaD1=GRAaD1+STF(M,I)*FA**2
      ENDDO
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE GCOaT
C  *********************************************************************
      SUBROUTINE GCOaT(E,CS,M)
C
C  Total cross section for incoherent (Compton) scattering. Relativistic
C  Impulse approximation with analytical Compton profiles.
C
C  Input arguments:
C    E ........ photon energy (eV).
C    M ........ material where photons propagate.
C  Output argument:
C    CS ....... incoherent total cross section (cm**2/molecule).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (MAXMAT=10)
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940325D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C
      COMMON/CGCO00/EE,EP,MM
      EXTERNAL GCOaD
C
      IF(E.LT.5.0D6) THEN
        EE=E
        MM=M
        CS=SUMGA(GCOaD,-1.0D0,1.0D0,1.0D-5)
      ELSE
C  ****  Klein-Nishina total cross section.
        EK=E/REV
        EK3=EK*EK
        EK2=1.0D0+EK+EK
        EK1=EK3-EK2-1.0D0
        T0=1.0D0/(1.0D0+EK+EK)
        CSL=0.5D0*EK3*T0*T0+EK2*T0+EK1*LOG(T0)-1.0D0/T0
        CS=0.0D0
        DO 1 I=1,NOSCCO(M)
          TAU=(E-UICO(M,I))/E
          IF(TAU.LT.T0) GO TO 1
          CSU=0.5D0*EK3*TAU*TAU+EK2*TAU+EK1*LOG(TAU)-1.0D0/TAU
          CS=CS+FCO(M,I)*(CSU-CSL)
    1   CONTINUE
        CS=PIELR2*CS/(EK*EK3)
      ENDIF
      RETURN
      END

C  *********************************************************************
C                        FUNCTION GCOaD
C  *********************************************************************
      FUNCTION GCOaD(CDT)
C
C  Single differential cross section for photon Compton scattering, dif-
C  ferential in the direction of the scattered photon only. Evaluated
C  from the incoherent scattering function.
C
C  The energy E of the primary photon is entered through common CGCO00.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (MAXMAT=10)
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
      PARAMETER (ELRAD=2.817940325D-13)  ! Class. electron radius (cm)
      PARAMETER (PI=3.1415926535897932D0, PIELR2=PI*ELRAD*ELRAD)
      PARAMETER (D2=1.4142135623731D0, D1=1.0D0/D2, D12=0.5D0)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C
      COMMON/CGCO00/E,EP,M
C
      CDT1=1.0D0-CDT
C  ****  Energy of the Compton line.
      EOEC=1.0D0+(E/REV)*CDT1
      ECOE=1.0D0/EOEC
C  ****  Incoherent scattering function (analytical profile).
      SIA=0.0D0
      DO 1 I=1,NOSCCO(M)
        IF(E.LT.UICO(M,I)) GO TO 1
        AUX=E*(E-UICO(M,I))*CDT1
        PZIMAX=(AUX-REV*UICO(M,I))/(REV*SQRT(AUX+AUX+UICO(M,I)**2))
        X=FJ0(M,I)*PZIMAX
        IF(X.GT.0.0D0) THEN
          SIAP=1.0D0-0.5D0*EXP(D12-(D1+D2*X)**2)
        ELSE
          SIAP=0.5D0*EXP(D12-(D1-D2*X)**2)
        ENDIF
        SIA=SIA+FCO(M,I)*SIAP
    1 CONTINUE
C  ****  Klein-Nishina X-factor.
      XKN=EOEC+ECOE-1.0D0+CDT*CDT
C  ****  Differential cross section.
      GCOaD=PIELR2*ECOE**2*XKN*SIA
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE GPHa0
C  *********************************************************************
      SUBROUTINE GPHa0
C
C  This subroutine sets all variables in common /CGPH00/ to zero.
C  It has to be invoked before reading the first material definition
C  file.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),NCUR
C
      DO I=1,99
        NPHS(I)=0
        IPHF(I)=0
        IPHL(I)=0
      ENDDO
C
      DO I=1,NTP
        EPH(I)=0.0D0
        DO J=1,10
          XPH(I,J)=1.0D-35
        ENDDO
      ENDDO
      NCUR=0
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE GPPa0
C  *********************************************************************
      SUBROUTINE GPPa0(M)
C
C  Initialisation of the sampling algorithm for electron-positron pair
C  production by photons in material M. Bethe-Heitler differential cross
C  section.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (REV=5.10998918D5)  ! Electron rest energy (eV)
      PARAMETER (SL=137.03599911D0)  ! Speed of light (1/alpha)
C  ****  Element data.
      CHARACTER*2 LASYMB
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Pair-production cross section parameters.
      COMMON/CGPP00/ZEQPP(MAXMAT),F0(MAXMAT,2),BCB(MAXMAT)
C
C  ***  Effective atomic number.
C
      FACT=0.0D0
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        FACT=FACT+IZZ*ATW(IZZ)*STF(M,I)
      ENDDO
      ZEQPP(M)=FACT/AT(M)
      IZZ=ZEQPP(M)+0.25D0
      IF(IZZ.LE.0) IZZ=1
      IF(IZZ.GT.99) IZZ=99
C  ****  DBM Coulomb correction.
      ALZ=ZEQPP(M)/SL
      A=ALZ*ALZ
      FC=A*(0.202059D0-A*(0.03693D0-A*(0.00835D0-A*(0.00201D0-A*
     1 (0.00049D0-A*(0.00012D0-A*0.00003D0)))))+1.0D0/(A+1.0D0))
C  ****  Screening functions and low-energy correction.
      BCB(M)=2.0D0/RSCR(IZZ)
      F0(M,1)=4.0D0*LOG(RSCR(IZZ))
      F0(M,2)=F0(M,1)-4.0D0*FC
      RETURN
      END

C  *********************************************************************
C                       BLOCK DATA PENDAT
C  *********************************************************************
      BLOCK DATA PENDAT
C
C  Physical data for the elements Z=1-99.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
C
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C
C  ************  Chemical symbols of the elements.
C
      DATA LASYMB/'H ','He','Li','Be','B ','C ','N ','O ','F ',
     1     'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ',
     2     'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
     3     'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ',
     4     'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     5     'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr',
     6     'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',
     7     'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au',
     8     'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     9     'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es'/
C
C  ************  Molar masses of the elements (g/mol).
C
      DATA ATW  /1.0079D0,4.0026D0,6.9410D0,9.0122D0,1.0811D1,
     1  1.2011D1,1.4007D1,1.5999D1,1.8998D1,2.0179D1,2.2990D1,
     2  2.4305D1,2.6982D1,2.8086D1,3.0974D1,3.2066D1,3.5453D1,
     3  3.9948D1,3.9098D1,4.0078D1,4.4956D1,4.7880D1,5.0942D1,
     4  5.1996D1,5.4938D1,5.5847D1,5.8933D1,5.8690D1,6.3546D1,
     5  6.5390D1,6.9723D1,7.2610D1,7.4922D1,7.8960D1,7.9904D1,
     6  8.3800D1,8.5468D1,8.7620D1,8.8906D1,9.1224D1,9.2906D1,
     7  9.5940D1,9.7907D1,1.0107D2,1.0291D2,1.0642D2,1.0787D2,
     8  1.1241D2,1.1482D2,1.1871D2,1.2175D2,1.2760D2,1.2690D2,
     9  1.3129D2,1.3291D2,1.3733D2,1.3891D2,1.4012D2,1.4091D2,
     A  1.4424D2,1.4491D2,1.5036D2,1.5196D2,1.5725D2,1.5893D2,
     B  1.6250D2,1.6493D2,1.6726D2,1.6893D2,1.7304D2,1.7497D2,
     C  1.7849D2,1.8095D2,1.8385D2,1.8621D2,1.9020D2,1.9222D2,
     D  1.9508D2,1.9697D2,2.0059D2,2.0438D2,2.0720D2,2.0898D2,
     E  2.0898D2,2.0999D2,2.2202D2,2.2302D2,2.2603D2,2.2703D2,
     F  2.3204D2,2.3104D2,2.3803D2,2.3705D2,2.3905D2,2.4306D2,
     G  2.4707D2,2.4707D2,2.5108D2,2.5208D2/
C
C  ************  Mean excitation energies of the elements (eV).
C
      DATA EPX / 19.2D0, 41.8D0, 40.0D0, 63.7D0, 76.0D0, 81.0D0,
     1   82.0D0, 95.0D0,115.0D0,137.0D0,149.0D0,156.0D0,166.0D0,
     2  173.0D0,173.0D0,180.0D0,174.0D0,188.0D0,190.0D0,191.0D0,
     3  216.0D0,233.0D0,245.0D0,257.0D0,272.0D0,286.0D0,297.0D0,
     4  311.0D0,322.0D0,330.0D0,334.0D0,350.0D0,347.0D0,348.0D0,
     5  343.0D0,352.0D0,363.0D0,366.0D0,379.0D0,393.0D0,417.0D0,
     6  424.0D0,428.0D0,441.0D0,449.0D0,470.0D0,470.0D0,469.0D0,
     7  488.0D0,488.0D0,487.0D0,485.0D0,491.0D0,482.0D0,488.0D0,
     8  491.0D0,501.0D0,523.0D0,535.0D0,546.0D0,560.0D0,574.0D0,
     9  580.0D0,591.0D0,614.0D0,628.0D0,650.0D0,658.0D0,674.0D0,
     A  684.0D0,694.0D0,705.0D0,718.0D0,727.0D0,736.0D0,746.0D0,
     B  757.0D0,790.0D0,790.0D0,800.0D0,810.0D0,823.0D0,823.0D0,
     C  830.0D0,825.0D0,794.0D0,827.0D0,826.0D0,841.0D0,847.0D0,
     D  878.0D0,890.0D0,902.0D0,921.0D0,934.0D0,939.0D0,952.0D0,
     E  966.0D0,980.0D0/
C
C  ************  Atomic form factor parameters.
C
      DATA RA1/0.0D0, 3.9265D+0, 4.3100D+1, 5.2757D+1, 2.5021D+1,
     1     1.2211D+1, 9.3229D+0, 3.2455D+0, 2.4197D+0, 1.5985D+0,
     2     3.0926D+1, 1.5315D+1, 7.7061D+0, 3.9493D+0, 2.2042D+0,
     3     1.9453D+1, 1.9354D+1, 8.0374D+0, 8.3779D+1, 5.7370D+1,
     4     5.2310D+1, 4.7514D+1, 4.3785D+1, 4.2048D+1, 3.6972D+1,
     5     3.3849D+1, 3.1609D+1, 2.8763D+1, 2.7217D+1, 2.4263D+1,
     6     2.2403D+1, 1.8606D+1, 1.5143D+1, 1.4226D+1, 1.1792D+1,
     7     9.7574D+0, 1.2796D+1, 1.2854D+1, 1.2368D+1, 1.0208D+1,
     8     8.2823D+0, 7.4677D+0, 7.6028D+0, 6.1090D+0, 5.5346D+0,
     9     4.2340D+0, 4.0444D+0, 4.2905D+0, 4.7950D+0, 5.1112D+0,
     A     5.2407D+0, 5.2153D+0, 5.1639D+0, 4.8814D+0, 5.8054D+0,
     B     6.6724D+0, 6.5104D+0, 6.3364D+0, 6.2889D+0, 6.3028D+0,
     C     6.3853D+0, 6.3475D+0, 6.5779D+0, 6.8486D+0, 7.0993D+0,
     D     7.6122D+0, 7.9681D+0, 8.3481D+0, 6.3875D+0, 8.0042D+0,
     E     8.0820D+0, 7.6940D+0, 7.1927D+0, 6.6751D+0, 6.1623D+0,
     F     5.8335D+0, 5.5599D+0, 4.6551D+0, 4.4327D+0, 4.7601D+0,
     G     5.2872D+0, 5.6084D+0, 5.7680D+0, 5.8041D+0, 5.7566D+0,
     H     5.6541D+0, 6.3932D+0, 6.9313D+0, 7.0027D+0, 6.8796D+0,
     I     6.4739D+0, 6.2405D+0, 6.0081D+0, 5.5708D+0, 5.3680D+0,
     J     5.8660D+0, 5.6375D+0, 5.1719D+0, 4.8989D+0/
      DATA RA2/0.0D0, 1.3426D-1, 9.4875D+1,-1.0896D+2,-4.5494D+1,
     1    -1.9572D+1,-1.2382D+1,-3.6827D+0,-2.4542D+0,-1.4453D+0,
     2     1.3401D+2, 7.9717D+1, 6.2164D+1, 4.0300D+1, 3.1682D+1,
     3    -1.3639D+1,-1.5950D+1,-5.1523D+0, 1.8351D+2, 1.2205D+2,
     4     1.0007D+2, 8.5632D+1, 7.9145D+1, 6.3675D+1, 6.2954D+1,
     5     5.6601D+1, 5.4171D+1, 4.8752D+1, 3.8062D+1, 3.9933D+1,
     6     4.8343D+1, 4.2137D+1, 3.4617D+1, 2.9430D+1, 2.4010D+1,
     7     1.9744D+1, 4.0009D+1, 5.1614D+1, 5.0456D+1, 3.9088D+1,
     8     2.6824D+1, 2.2953D+1, 2.4773D+1, 1.6893D+1, 1.4548D+1,
     9     9.7226D+0, 1.0192D+1, 1.1153D+1, 1.3188D+1, 1.4733D+1,
     A     1.5644D+1, 1.5939D+1, 1.5923D+1, 1.5254D+1, 2.0748D+1,
     B     2.6901D+1, 2.7032D+1, 2.4938D+1, 2.1528D+1, 2.0362D+1,
     C     1.9474D+1, 1.8238D+1, 1.7898D+1, 1.9174D+1, 1.9023D+1,
     D     1.8194D+1, 1.8504D+1, 1.8955D+1, 1.4276D+1, 1.7558D+1,
     E     1.8651D+1, 1.7984D+1, 1.6793D+1, 1.5469D+1, 1.4143D+1,
     F     1.3149D+1, 1.2255D+1, 9.2352D+0, 8.6067D+0, 9.7460D+0,
     G     1.1749D+1, 1.3281D+1, 1.4326D+1, 1.4920D+1, 1.5157D+1,
     H     1.5131D+1, 1.9489D+1, 2.3649D+1, 2.4686D+1, 2.4760D+1,
     I     2.1519D+1, 2.0099D+1, 1.8746D+1, 1.5943D+1, 1.4880D+1,
     J     1.6345D+1, 1.5283D+1, 1.3061D+1, 1.1553D+1/
      DATA RA3/0.0D0, 2.2648D+0, 1.0579D+3, 8.6177D+2, 2.4422D+2,
     1     7.8788D+1, 3.8293D+1, 1.2564D+1, 6.9091D+0, 3.7926D+0,
     2     0.0000D+0, 0.0000D+0, 1.6759D-9, 1.3026D+1, 3.0569D+0,
     3     1.5521D+2, 1.2815D+2, 4.7378D+1, 9.2802D+2, 4.7508D+2,
     4     3.6612D+2, 2.7582D+2, 2.1008D+2, 1.5903D+2, 1.2322D+2,
     5     9.2898D+1, 7.1345D+1, 5.1651D+1, 3.8474D+1, 2.7410D+1,
     6     1.9126D+1, 1.0889D+1, 5.3479D+0, 8.2223D+0, 5.0837D+0,
     7     2.8905D+0, 2.7457D+0, 6.7082D-1, 0.0000D+0, 0.0000D+0,
     8     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     9     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     A     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     B     0.0000D+0, 0.0000D+0, 0.0000D+0, 1.7264D-1, 2.7322D-1,
     C     3.9444D-1, 4.5648D-1, 6.2286D-1, 7.2468D-1, 8.4296D-1,
     D     1.1698D+0, 1.2994D+0, 1.4295D+0, 0.0000D+0, 8.1570D-1,
     E     6.9349D-1, 4.9536D-1, 3.1211D-1, 1.5931D-1, 2.9512D-2,
     F     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     G     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     H     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     I     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0,
     J     0.0000D+0, 0.0000D+0, 0.0000D+0, 0.0000D+0/
      DATA RA4/1.1055D1,6.3519D0,4.7367D+1, 3.9402D+1, 2.2896D+1,
     1     1.3979D+1, 1.0766D+1, 6.5252D+0, 5.1631D+0, 4.0524D+0,
     2     2.7145D+1, 1.8724D+1, 1.4782D+1, 1.1608D+1, 9.7750D+0,
     3     1.6170D+1, 1.5249D+1, 9.1916D+0, 5.4499D+1, 4.1381D+1,
     4     3.7395D+1, 3.3815D+1, 3.1135D+1, 2.8273D+1, 2.6140D+1,
     5     2.3948D+1, 2.2406D+1, 2.0484D+1, 1.8453D+1, 1.7386D+1,
     6     1.7301D+1, 1.5388D+1, 1.3411D+1, 1.2668D+1, 1.1133D+1,
     7     9.8081D+0, 1.3031D+1, 1.4143D+1, 1.3815D+1, 1.2077D+1,
     8     1.0033D+1, 9.2549D+0, 9.5338D+0, 7.9076D+0, 7.3263D+0,
     9     5.9996D+0, 6.0087D+0, 6.2660D+0, 6.7914D+0, 7.1501D+0,
     A     7.3367D+0, 7.3729D+0, 7.3508D+0, 7.1465D+0, 8.2731D+0,
     B     9.3745D+0, 9.3508D+0, 8.9897D+0, 8.4566D+0, 8.2690D+0,
     C     8.1398D+0, 7.9183D+0, 7.9123D+0, 8.1677D+0, 8.1871D+0,
     D     8.1766D+0, 8.2881D+0, 8.4227D+0, 7.0273D+0, 8.0002D+0,
     E     8.1440D+0, 7.9104D+0, 7.5685D+0, 7.1970D+0, 6.8184D+0,
     F     6.5469D+0, 6.3056D+0, 5.4844D+0, 5.2832D+0, 5.5889D+0,
     G     6.0919D+0, 6.4340D+0, 6.6426D+0, 6.7428D+0, 6.7636D+0,
     H     6.7281D+0, 7.5729D+0, 8.2808D+0, 8.4400D+0, 8.4220D+0,
     I     7.8662D+0, 7.5993D+0, 7.3353D+0, 6.7829D+0, 6.5520D+0,
     J     6.9181D+0, 6.6794D+0, 6.1735D+0, 5.8332D+0/
      DATA RA5/0.0D0, 4.9828D+0, 5.5674D+1, 3.0902D+1, 1.1496D+1,
     1     4.8936D+0, 2.5506D+0, 1.2236D+0, 7.4698D-1, 4.7042D-1,
     2     4.7809D+0, 4.6315D+0, 4.3677D+0, 4.9269D+0, 2.6033D+0,
     3     9.6229D+0, 7.2592D+0, 4.1634D+0, 1.3999D+1, 8.6975D+0,
     4     6.9630D+0, 5.4681D+0, 4.2653D+0, 3.2848D+0, 2.7354D+0,
     5     2.1617D+0, 1.7030D+0, 1.2826D+0, 9.7080D-1, 7.2227D-1,
     6     5.0874D-1, 3.1402D-1, 1.6360D-1, 3.2918D-1, 2.3570D-1,
     7     1.5868D-1, 1.5146D-1, 9.7662D-2, 7.3151D-2, 6.4206D-2,
     8     4.8945D-2, 4.3189D-2, 4.4368D-2, 3.3976D-2, 3.0466D-2,
     9     2.4477D-2, 3.7202D-2, 3.7093D-2, 3.8161D-2, 3.8576D-2,
     A     3.8403D-2, 3.7806D-2, 3.4958D-2, 3.6029D-2, 4.3087D-2,
     B     4.7069D-2, 4.6452D-2, 4.2486D-2, 4.1517D-2, 4.1691D-2,
     C     4.2813D-2, 4.2294D-2, 4.5287D-2, 4.8462D-2, 4.9726D-2,
     D     5.5097D-2, 5.6568D-2, 5.8069D-2, 1.2270D-2, 3.8006D-2,
     E     3.5048D-2, 3.0050D-2, 2.5069D-2, 2.0485D-2, 1.6151D-2,
     F     1.4631D-2, 1.4034D-2, 1.1978D-2, 1.1522D-2, 1.2375D-2,
     G     1.3805D-2, 1.4954D-2, 1.5832D-2, 1.6467D-2, 1.6896D-2,
     H     1.7166D-2, 1.9954D-2, 2.2497D-2, 2.1942D-2, 2.1965D-2,
     I     2.0005D-2, 1.8927D-2, 1.8167D-2, 1.6314D-2, 1.5522D-2,
     J     1.3141D-2, 1.2578D-2, 1.2238D-2, 1.0427D-2/
C
C  ************  Pair-production cross section parameters.
C
C  ****  Screening parameter (R mc/hbar).
      DATA RSCR  /1.2281D2,7.3167D1,6.9228D1,6.7301D1,6.4696D1,
     1   6.1228D1,5.7524D1,5.4033D1,5.0787D1,4.7851D1,4.6373D1,
     2   4.5401D1,4.4503D1,4.3815D1,4.3074D1,4.2321D1,4.1586D1,
     3   4.0953D1,4.0524D1,4.0256D1,3.9756D1,3.9144D1,3.8462D1,
     4   3.7778D1,3.7174D1,3.6663D1,3.5986D1,3.5317D1,3.4688D1,
     5   3.4197D1,3.3786D1,3.3422D1,3.3068D1,3.2740D1,3.2438D1,
     6   3.2143D1,3.1884D1,3.1622D1,3.1438D1,3.1142D1,3.0950D1,
     7   3.0758D1,3.0561D1,3.0285D1,3.0097D1,2.9832D1,2.9581D1,
     8   2.9411D1,2.9247D1,2.9085D1,2.8930D1,2.8721D1,2.8580D1,
     9   2.8442D1,2.8312D1,2.8139D1,2.7973D1,2.7819D1,2.7675D1,
     A   2.7496D1,2.7285D1,2.7093D1,2.6911D1,2.6705D1,2.6516D1,
     B   2.6304D1,2.6108D1,2.5929D1,2.5730D1,2.5577D1,2.5403D1,
     C   2.5245D1,2.5100D1,2.4941D1,2.4790D1,2.4655D1,2.4506D1,
     D   2.4391D1,2.4262D1,2.4145D1,2.4039D1,2.3922D1,2.3813D1,
     E   2.3712D1,2.3621D1,2.3523D1,2.3430D1,2.3331D1,2.3238D1,
     F   2.3139D1,2.3048D1,2.2967D1,2.2833D1,2.2694D1,2.2624D1,
     G   2.2545D1,2.2446D1,2.2358D1,2.2264D1/
C  ****  Asymptotic triplet contribution (eta).
      DATA ETA   /1.1570D0,1.1690D0,1.2190D0,1.2010D0,1.1890D0,
     1   1.1740D0,1.1760D0,1.1690D0,1.1630D0,1.1570D0,1.1740D0,
     2   1.1830D0,1.1860D0,1.1840D0,1.1800D0,1.1780D0,1.1750D0,
     3   1.1700D0,1.1800D0,1.1870D0,1.1840D0,1.1800D0,1.1770D0,
     4   1.1660D0,1.1690D0,1.1660D0,1.1640D0,1.1620D0,1.1540D0,
     5   1.1560D0,1.1570D0,1.1580D0,1.1570D0,1.1580D0,1.1580D0,
     6   1.1580D0,1.1660D0,1.1730D0,1.1740D0,1.1750D0,1.1700D0,
     7   1.1690D0,1.1720D0,1.1690D0,1.1680D0,1.1640D0,1.1670D0,
     8   1.1700D0,1.1720D0,1.1740D0,1.1750D0,1.1780D0,1.1790D0,
     9   1.1800D0,1.1870D0,1.1940D0,1.1970D0,1.1960D0,1.1940D0,
     A   1.1940D0,1.1940D0,1.1940D0,1.1940D0,1.1960D0,1.1970D0,
     B   1.1960D0,1.1970D0,1.1970D0,1.1980D0,1.1980D0,1.2000D0,
     C   1.2010D0,1.2020D0,1.2040D0,1.2050D0,1.2060D0,1.2080D0,
     D   1.2070D0,1.2080D0,1.2120D0,1.2150D0,1.2180D0,1.2210D0,
     E   1.2240D0,1.2270D0,1.2300D0,1.2370D0,1.2430D0,1.2470D0,
     F   1.2500D0,1.2510D0,1.2520D0,1.2550D0,1.2560D0,1.2570D0,
     G   1.2590D0,1.2620D0,1.2620D0,1.2650D0/
C
      END

C  *********************************************************************
C                       SUBROUTINE SPLINE
C  *********************************************************************
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C
C  Cubic spline interpolation of tabulated data.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     Y(I) (I=1:N) ... corresponding function values.
C     S1,SN .......... second derivatives at X(1) and X(N). The natural
C                      spline corresponds to taking S1=SN=0.
C     N .............. number of grid points.
C  Output:
C     A(1:N),B(1:N),C(1:N),D(1:N) ... spline coefficients.
C
C  The interpolating cubic polynomial in the I-th interval, from X(I) to
C  X(I+1), is
C               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))
C
C  Reference: M.J. Maron, 'Numerical Analysis: a Practical Approach',
C             MacMillan Publ. Co., New York, 1982.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      IF(N.LT.4) THEN
        WRITE(26,10) N
   10   FORMAT(5X,'Spline interpolation cannot be performed with',
     1    I4,' points. Stop.')
        STOP 'SPLINE. N is less than 4.'
      ENDIF
      N1=N-1
      N2=N-2
C  ****  Auxiliary arrays H(=A) and DELTA(=D).
      DO I=1,N1
        IF(X(I+1)-X(I).LT.1.0D-13) THEN
          WRITE(26,11)
   11     FORMAT(5X,'Spline X values not in increasing order. Stop.')
          STOP 'SPLINE. X values not in increasing order.'
        ENDIF
        A(I)=X(I+1)-X(I)
        D(I)=(Y(I+1)-Y(I))/A(I)
      ENDDO
C  ****  Symmetric coefficient matrix (augmented).
      DO I=1,N2
        B(I)=2.0D0*(A(I)+A(I+1))
        K=N1-I+1
        D(K)=6.0D0*(D(K)-D(K-1))
      ENDDO
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
C  ****  Gauss solution of the tridiagonal system.
      DO I=2,N2
        R=A(I)/B(I-1)
        B(I)=B(I)-R*A(I)
        D(I+1)=D(I+1)-R*D(I)
      ENDDO
C  ****  The SIGMA coefficients are stored in array D.
      D(N1)=D(N1)/B(N2)
      DO I=2,N2
        K=N1-I+1
        D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      ENDDO
      D(N)=SN
C  ****  Spline coefficients.
      SI1=S1
      DO I=1,N1
        SI=SI1
        SI1=D(I+1)
        H=A(I)
        HI=1.0D0/H
        A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)
     1      +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))
     2      +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))
        B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)
     1      +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)
        C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))
        D(I)=(HI/6.0D0)*(SI1-SI)
      ENDDO
C  ****  Natural cubic spline for X.GT.X(N).
      FN=Y(N)
      FNP=B(N1)+X(N)*(2.0D0*C(N1)+X(N)*3.0D0*D(N1))
      A(N)=FN-X(N)*FNP
      B(N)=FNP
      C(N)=0.0D0
      D(N)=0.0D0
C
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE FINDI
C  *********************************************************************
      SUBROUTINE FINDI(X,XC,N,I)
C
C  Finds the interval (X(I),X(I+1)) that contains the value XC.
C
C  Input:
C     X(I) (I=1:N) ... grid points (the X values must be in increasing
C                      order).
C     XC ............. point to be located.
C     N  ............. number of grid points.
C  Output:
C     I .............. interval index.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N)
C
      IF(XC.GT.X(N)) THEN
        I=N
        RETURN
      ENDIF
      IF(XC.LT.X(1)) THEN
        I=1
        RETURN
      ENDIF
      I=1
      I1=N
    1 IT=(I+I1)/2
      IF(XC.GT.X(IT)) THEN
        I=IT
      ELSE
        I1=IT
      ENDIF
      IF(I1-I.GT.1) GO TO 1
      RETURN
      END

C  *********************************************************************
C                       FUNCTION SUMGA
C  *********************************************************************
      FUNCTION SUMGA(FCT,XL,XU,TOL)
C
C  This function calculates the value SUMGA of the integral of the
C  (external) function FCT over the interval (XL,XU) using the 20-point
C  Gauss quadrature method with an adaptive bipartition scheme.
C
C  TOL is the tolerance, i.e. maximum allowed relative error; it should
C  not exceed 1.0D-13. A warning message is written in unit 6 when the
C  required accuracy is not attained.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NP=10, NST=256, NCALLS=20000)
      DIMENSION X(NP),W(NP),S(NST),SN(NST),XR(NST),XRN(NST)
C  ****  Gauss 20-point integration formula.
C  Abscissas.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  Weights.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C  ****  Error control.
      CTOL=MIN(MAX(TOL,1.0D-13),1.0D-2)
      PTOL=0.01D0*CTOL
      ERR=1.0D35
C  ****  Gauss integration from XL to XU.
      H=XU-XL
      SUMGA=0.0D0
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO I1=2,NP
        C=A*X(I1)
        D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      ENDDO
      ICALL=NP+NP
      LH=1
      S(1)=D*A
      XR(1)=XL
C  ****  Adaptive bipartition scheme.
    1 CONTINUE
      HO=H
      H=0.5D0*H
      SUMR=0.0D0
      LHN=0
      DO I=1,LH
        SI=S(I)
        XA=XR(I)
        XB=XA+H
        XC=XA+HO
        A=0.5D0*(XB-XA)
        B=0.5D0*(XB+XA)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I2=2,NP
          C=A*X(I2)
          D=D+W(I2)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S1=D*A
        A=0.5D0*(XC-XB)
        B=0.5D0*(XC+XB)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I3=2,NP
          C=A*X(I3)
          D=D+W(I3)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S2=D*A
        ICALL=ICALL+4*NP
        S12=S1+S2
        IF(ABS(S12-SI).LE.MAX(PTOL*ABS(S12),1.0D-35)) THEN
          SUMGA=SUMGA+S12
        ELSE
          SUMR=SUMR+S12
          LHN=LHN+2
          IF(LHN.GE.NST) GO TO 2
          SN(LHN)=S2
          XRN(LHN)=XB
          SN(LHN-1)=S1
          XRN(LHN-1)=XA
        ENDIF
        IF(ICALL.GT.NCALLS) GO TO 2
      ENDDO
      ERR=ABS(SUMR)/MAX(ABS(SUMR+SUMGA),1.0D-35)
      IF(ERR.LT.CTOL.OR.LHN.EQ.0) RETURN
      LH=LHN
      DO I=1,LH
        S(I)=SN(I)
        XR(I)=XRN(I)
      ENDDO
      GO TO 1
C  ****  Warning (low accuracy) message.
    2 CONTINUE
      WRITE(26,11)
   11 FORMAT(/2X,'>>> SUMGA. Gauss adaptive-bipartition quadrature.')
      WRITE(26,12) XL,XU,TOL
   12 FORMAT(2X,'XL =',1P,E19.12,',  XU =',E19.12,',  TOL =',E8.1)
      WRITE(26,13) ICALL,SUMGA,ERR,LHN
   13 FORMAT(2X,'NCALLS = ',I5,',  SUMGA =',1P,E20.13,',  ERR =',E8.1,
     1      /2X,'Number of open subintervals = ',I3)
      WRITE(26,14)
   14 FORMAT(2X,'WARNING: the required accuracy has not been ',
     1  'attained.'/)
      RETURN
      END

C  *********************************************************************
C                        SUBROUTINE RITAI0
C  *********************************************************************
      SUBROUTINE RITAI0(PDF,XMIN,XMAX,N,NU,ERRM,IWR)
C
C     Initialisation of the RITA algorithm for random sampling of a
C  continuous random variable X from a probability distribution function
C  PDF(X) defined in the domain (XMIN,XMAX). N is the number of points
C  in the sampling grid. These points are determined by means of an
C  adaptive strategy that minimises local interpolation errors; the
C  first NU grid points are uniformly spaced in (XMIN,XMAX).
C
C  ****  Interpolation coefficients and PDF tables are printed on
C        separate files (UNIT=8) if IWR=1.
C
C  Other subprograms needed: EXTERNAL function PDF,
C                            subroutine SIMPSU.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10, ZERO=1.0D-35, ZEROT=0.1D0*ZERO)
      PARAMETER (NM=512)
C
C     The information used by the sampling function RITAI is exported
C  through the following common block,
      COMMON/CRITA/X(NM),PAC(NM),A(NM),B(NM),NP,
     1             ITTL(NM),ITTU(NM),NPM1
C  where
C    X(I) ...... grid points, in increasing order.
C    PAC(I) .... value of the cumulative pdf at X(I).
C    A(I), B(I) ... rational inverse cumulative distribution parameters.
C    NP ........ number of grid points (8 .LE. NP .LE. NM).
C
C    ITTL(I) ... largest J for which PAC(J) < (I-1)/(NP-1).
C    ITTU(I) ... smallest K for which PAC(K) > I/(NP-1).
C    NPM1 ...... =NP-1.
C
      PARAMETER (NIP=51)
      DIMENSION AREA(NM),ERR(NM),C(NM),XI(NIP),PDFI(NIP),SUMI(NIP)
      EXTERNAL PDF
C
      IF(N.LE.16) THEN
        WRITE(26,'('' Error in RITAI0: N must be larger than 16.'',
     1    ''  N='',I3)') N
        STOP 'RITAI0: N must be larger than 16.'
      ENDIF
      IF(N.GT.NM) THEN
        WRITE(26,'('' Error in RITAI0: N must be less than NM=512.'')')
        STOP 'RITAI0: N must be less than NM=512.'
      ENDIF
      IF(XMIN.GT.XMAX-EPS) THEN
        WRITE(26,'('' Error in RITAI0: XMIN must be larger than XMAX.'')
     1    ')
        STOP 'RITAI0: XMIN must be larger than XMAX.'
      ENDIF
C
C  ****  We start with a grid of NUNIF points uniformly spaced in the
C        interval (XMIN,XMAX).
C
      NUNIF=MIN(MAX(8,NU),N/2)
      NP=NUNIF
      DX=(XMAX-XMIN)/DBLE(NP-1)
      X(1)=XMIN
      DO I=1,NP-1
        X(I+1)=XMIN+I*DX
      ENDDO
      X(NP)=XMAX
C
      DO I=1,NP-1
        DXI=(X(I+1)-X(I))/DBLE(NIP-1)
        PDFMAX=0.0D0
        DO K=1,NIP
          XI(K)=X(I)+DBLE(K-1)*DXI
          PDFI(K)=MAX(PDF(XI(K)),ZEROT)
          PDFMAX=MAX(PDFMAX,PDFI(K))
        ENDDO
        CALL SIMPSU(DXI,PDFI,SUMI,NIP)
        AREA(I)=SUMI(NIP)
        FACT=1.0D0/AREA(I)
        DO K=1,NIP
          SUMI(K)=FACT*SUMI(K)
        ENDDO
C  ****  When the PDF vanishes at one of the interval end points, its
C        value is modified.
        IF(PDFI(1).LT.ZERO) PDFI(1)=1.0D-5*PDFMAX
        IF(PDFI(NIP).LT.ZERO) PDFI(NIP)=1.0D-5*PDFMAX
C
        PLI=PDFI(1)*FACT
        PUI=PDFI(NIP)*FACT
        B(I)=1.0D0-1.0D0/(PLI*PUI*DX*DX)
        A(I)=(1.0D0/(PLI*DX))-1.0D0-B(I)
        C(I)=1.0D0+A(I)+B(I)
        IF(C(I).LT.ZERO) THEN
          A(I)=0.0D0
          B(I)=0.0D0
          C(I)=1.0D0
        ENDIF
C
C  ****  ERR(I) is defined as the integral of the absolute difference
C        between the rational interpolation and the true PDF, extended
C        over the interval (X(I),X(I+1)).
C
        ICASE=1
  100   CONTINUE
        ERR(I)=0.0D0
        DO K=1,NIP
          RR=SUMI(K)
          PAP=AREA(I)*(1.0D0+(A(I)+B(I)*RR)*RR)**2/
     1       ((1.0D0-B(I)*RR*RR)*C(I)*(X(I+1)-X(I)))
          IF(K.EQ.1.OR.K.EQ.NIP) THEN
            ERR(I)=ERR(I)+0.5D0*ABS(PAP-PDFI(K))
          ELSE
            ERR(I)=ERR(I)+ABS(PAP-PDFI(K))
          ENDIF
        ENDDO
        ERR(I)=ERR(I)*DXI
C  ****  If ERR(I) is too large, the PDF is approximated by a uniform
C        distribution.
        IF(ERR(I).GT.0.10D0*AREA(I).AND.ICASE.EQ.1) THEN
          B(I)=0.0D0
          A(I)=0.0D0
          C(I)=1.0D0
          ICASE=2
          GO TO 100
        ENDIF
      ENDDO
      X(NP)=XMAX
      A(NP)=0.0D0
      B(NP)=0.0D0
      C(NP)=0.0D0
      ERR(NP)=0.0D0
      AREA(NP)=0.0D0
C
C  ****  New grid points are added by halving the subinterval with the
C        largest absolute error.
C
  200 CONTINUE
      ERRM=0.0D0
      LMAX=1
      DO I=1,NP-1
C  ****  ERRM is the largest of the interval errors ERR(I).
        IF(ERR(I).GT.ERRM) THEN
          ERRM=ERR(I)
          LMAX=I
        ENDIF
      ENDDO
C
      NP=NP+1
      DO I=NP,LMAX+1,-1
        X(I)=X(I-1)
        A(I)=A(I-1)
        B(I)=B(I-1)
        C(I)=C(I-1)
        ERR(I)=ERR(I-1)
        AREA(I)=AREA(I-1)
      ENDDO
      X(LMAX+1)=0.5D0*(X(LMAX)+X(LMAX+2))
      DO I=LMAX,LMAX+1
        DX=X(I+1)-X(I)
        DXI=(X(I+1)-X(I))/DBLE(NIP-1)
        PDFMAX=0.0D0
        DO K=1,NIP
          XI(K)=X(I)+DBLE(K-1)*DXI
          PDFI(K)=MAX(PDF(XI(K)),ZEROT)
          PDFMAX=MAX(PDFMAX,PDFI(K))
        ENDDO
        CALL SIMPSU(DXI,PDFI,SUMI,NIP)
        AREA(I)=SUMI(NIP)
        FACT=1.0D0/AREA(I)
        DO K=1,NIP
          SUMI(K)=FACT*SUMI(K)
        ENDDO
C
        IF(PDFI(1).LT.ZERO) PDFI(1)=1.0D-5*PDFMAX
        IF(PDFI(NIP).LT.ZERO) PDFI(NIP)=1.0D-5*PDFMAX
        PLI=PDFI(1)*FACT
        PUI=PDFI(NIP)*FACT
        B(I)=1.0D0-1.0D0/(PLI*PUI*DX*DX)
        A(I)=(1.0D0/(PLI*DX))-1.0D0-B(I)
        C(I)=1.0D0+A(I)+B(I)
        IF(C(I).LT.ZERO) THEN
          A(I)=0.0D0
          B(I)=0.0D0
          C(I)=1.0D0
        ENDIF
C
        ICASE=1
  300   CONTINUE
        ERR(I)=0.0D0
        DO K=1,NIP
          RR=SUMI(K)
          PAP=AREA(I)*(1.0D0+(A(I)+B(I)*RR)*RR)**2/
     1       ((1.0D0-B(I)*RR*RR)*C(I)*(X(I+1)-X(I)))
          IF(K.EQ.1.OR.K.EQ.NIP) THEN
            ERR(I)=ERR(I)+0.5D0*ABS(PAP-PDFI(K))
          ELSE
            ERR(I)=ERR(I)+ABS(PAP-PDFI(K))
          ENDIF
        ENDDO
        ERR(I)=ERR(I)*DXI
C
        IF(ERR(I).GT.0.10D0*AREA(I).AND.ICASE.EQ.1) THEN
          B(I)=0.0D0
          A(I)=0.0D0
          C(I)=1.0D0
          ICASE=2
          GO TO 300
        ENDIF
      ENDDO
C
      IF(NP.LT.N) GO TO 200
      NPM1=NP-1
C
C  ****  Renormalisation.
C
      WS=0.0D0
      DO I=1,NPM1
        WS=WS+AREA(I)
      ENDDO
      WS=1.0D0/WS
      ERRM=0.0D0
      DO I=1,NPM1
        AREA(I)=AREA(I)*WS
        ERR(I)=ERR(I)*WS
        ERRM=MAX(ERRM,ERR(I))
      ENDDO
C
      PAC(1)=0.0D0
      DO I=1,NPM1
        PAC(I+1)=PAC(I)+AREA(I)
      ENDDO
      PAC(NP)=1.0D0
C
C  ****  Pre-calculated limits for the initial binary search in
C        subroutine RITAI.
C
      BIN=1.0D0/DBLE(NPM1)
      ITTL(1)=1
      DO I=2,NPM1
        PTST=(I-1)*BIN
        DO J=ITTL(I-1),NP
          IF(PAC(J).GT.PTST) THEN
            ITTL(I)=J-1
            ITTU(I-1)=J
            GO TO 400
          ENDIF
        ENDDO
  400   CONTINUE
      ENDDO
      ITTU(NPM1)=NP
      ITTL(NP)=NP-1
      ITTU(NP)=NP
C
C  ****  Print interpolation tables (optional, active only if IWR=1).
C
      IF(IWR.EQ.1) THEN
        OPEN(8,FILE='param.dat')
        WRITE(8,1000)
 1000   FORMAT(1x,'#',5X,'X',11X,'PDF(X)',10X,'A',13X,'B',13X,'C',
     1    11X,'error')
        DO I=1,NPM1
          PDFE=MAX(PDF(X(I)),ZEROT)*WS
          WRITE(8,'(1P,7E14.6)') X(I),PDFE,A(I),B(I),C(I),ERR(I)
        ENDDO
        CLOSE(8)
C
        OPEN(8,FILE='table.dat')
        WRITE(8,2000)
 2000   FORMAT(1X,'#',6X,'X',13X,'PDF_ex',10X,'PDF_ap',11X,'err')
        DO I=1,NPM1
          DX=(X(I+1)-X(I))/DBLE(NIP-1)
          DO K=1,NIP
            XT=X(I)+(K-1)*DX
            P1=MAX(PDF(XT),ZEROT)*WS
C  ****  Rational interpolation.
            TAU=(XT-X(I))/(X(I+1)-X(I))
            CON1=2.0D0*B(I)*TAU
            CON2=C(I)-A(I)*TAU
            IF(ABS(CON1).GT.1.0D-10*ABS(CON2)) THEN
              ETA=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))/CON1
            ELSE
              ETA=TAU/CON2
            ENDIF
            P2=AREA(I)*(1.0D0+(A(I)+B(I)*ETA)*ETA)**2
     1        /((1.0D0-B(I)*ETA*ETA)*C(I)*(X(I+1)-X(I)))
            WRITE(8,'(1P,5E16.8)') XT,P1,P2,(P1-P2)/P1
          ENDDO
        ENDDO
        CLOSE(8)
C
        OPEN(8,FILE='limits.dat')
        WRITE(8,3000)
 3000   FORMAT(1X,'#  I',5X,'PAC(ITTL)',7X,'(I-1)/NPM1',9X,'I/NPM1',
     1    10X,'PAC(ITTU)')
        DO I=1,NPM1
          WRITE(8,'(I5,4E17.9)') I,PAC(ITTL(I)),(I-1)*BIN,I*BIN,
     1      PAC(ITTU(I))
          IF(PAC(ITTL(I)).GT.(I-1)*BIN+EPS.OR.
     1       PAC(ITTU(I)).LT.I*BIN-EPS) THEN
            WRITE(8,3001)
 3001       FORMAT(' #  WARNING: The four values should be in in',
     1        'creasing order.')
          ENDIF
        ENDDO
        CLOSE(8)
      ENDIF
C
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE SIMPSU
C  *********************************************************************
      SUBROUTINE SIMPSU(H,F,SUMF,N)
C
C     Simpson's integration of a uniformly tabulated function (coded
C  using the 3-point Lagrange quadrature formula).
C
C     H ...... grid spacing,
C     F ...... array of function values (ordered abscissas),
C     SUMF ... array of integral values defined as
C              SUMF(I)=INTEGRAL(F) from X(1) to X(I)=X(1)+(I-1)*H,
C     N ...... number of points in the table.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (F1O12=1.0D0/12.0D0)
      DIMENSION F(N),SUMF(N)
      CONS=H*F1O12
      HCONS=0.5D0*CONS
      IF(N.LT.4) STOP 'Not enough data points.'
      SUMF(1)=0.0D0
      SUMF(2)=SUMF(1)+(5.0D0*F(1)+8.0D0*F(2)-F(3))*CONS
      DO I=3,N-1
        SUMF(I)=SUMF(I-1)+(13.0D0*(F(I-1)+F(I))-F(I+1)-F(I-2))*HCONS
      ENDDO
      SUMF(N)=SUMF(N-1)+(5.0D0*F(N)+8.0D0*F(N-1)-F(N-2))*CONS
      RETURN
      END

