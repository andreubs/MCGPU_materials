
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   MCGPU_material_generator.f                                         C
C   ============================                                       C
C   Creates material data files (.mcgpu) for use with MC-GPU, the      C
C   GPU-accelerated Monte Carlo X-ray transport code.                  C
C                                                                      C
C   This program combines the material file creation routines from     C
C   PENELOPE's 'material.f' (subroutine PEMATW) with the MC-GPU        C
C   mean-free-path table generator 'MC-GPU_create_material_data.f'     C
C   into a single interactive utility. No PENELOPE material file       C
C   needed or created in the generation of the MCGPU material file.    C
C                                                                      C
C   The PENDBASE_photons/ database folder (101 files) must be present  C
C   in the current working directory. It is created by the companion   C
C   script create_PENDBASE_photons.py from the full PENDBASE.          C
C                                                                      C
C   Physics: photon transport at diagnostic and radiotherapy energies. C
C   Three interactions are simulated in MC-GPU:                        C
C     - Rayleigh (coherent) scattering                                 C
C     - Compton (incoherent) scattering                                C
C     - Photoelectric absorption                                       C
C   Pair production is NOT simulated (no e+/e- generated), but its     C
C   cross section IS included in the total MFP for correct primary     C
C   beam attenuation above 1.022 MeV. Fluorescence is excluded.        C
C                                                                      C
C   Compilation:                                                       C
C   gfortran MCGPU_materials.f penelope_photons.f -o MCGPU_materials.x C
C                                                                      C
C   References:                                                        C
C     - F. Salvat et al., PENELOPE 2006, NEA-OECD (2006)               C
C     - A. Badal & A. Badano, Med. Phys. 36, pp.4878-4880 (2009)       C
C                                                                      C
C   Based on PENELOPE/PENGEOM (version 2006)                           C
C   Copyright (c) 2001-2006 Universitat de Barcelona                   C
C   Adapted for MC-GPU by Andreu Badal, FDA/CDRH/OSEL/DIDSR            C
C   AI assistant: Claude Code Sonnet 4.6                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  *********************************************************************
C                           MAIN PROGRAM
C  *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*62 NAME
      CHARACTER*80 OUTNAME
      DIMENSION FBW(30)
      PARAMETER (AVOG=6.0221415D23)

C  ****  Material composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Simulation parameters (set defaults for photon-only, no e/p).
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Energy grid and interpolation constants.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),
     1  NCUR
C  ****  Rayleigh scattering: RITA sampling data.
      PARAMETER (NP=128,NPM1=NP-1)
      COMMON/CGRA/XCO(NP,MAXMAT),PCO(NP,MAXMAT),ACO(NP,MAXMAT),
     1  BCO(NP,MAXMAT),PMAX(NEGP,MAXMAT),ITLCO(NP,MAXMAT),
     2  ITUCO(NP,MAXMAT)
C  ****  Compton scattering shell data.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C  ****  Auxiliary arrays for MFP table.
      DIMENSION E_MFP(6)
      PARAMETER (MAX_ENERGY_BINS=60005)
      DIMENSION PMAX_linear_energy(MAX_ENERGY_BINS)

C  ****  Material number (always 1 in this program).
      M=1

C  ****  Set default simulation parameters (photon-only, disable e/p).
      DO MM=1,MAXMAT
        EABS(1,MM) = 50.0D0
        EABS(2,MM) = 50.0D0
        EABS(3,MM) = 50.0D0
        C1(MM)     = 0.0D0
        C2(MM)     = 0.0D0
        WCC(MM)    = 0.0D0
        WCR(MM)    =-10.0D0
      ENDDO

      WRITE(6,*) ' '
      WRITE(6,*) ' '
      WRITE(6,*)
     & '   *****************************************************'
      WRITE(6,*)
     & '   ***      MCGPU_material_generator  (2026)        ***'
      WRITE(6,*)
     & '   *****************************************************'
      WRITE(6,*) ' '
      WRITE(6,*)
     & '   Creates .mcgpu material files for MC-GPU simulations.'
      WRITE(6,*)
     & '   Based on PENELOPE 2006 photon interaction physics.'
      WRITE(6,*)
     & '   Database: PENDBASE_photons/ must be in current folder.'
      WRITE(6,*) ' '

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  PART A: Material composition input
C  (adapted from PENELOPE's PEMATW subroutine, photon-relevant only)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IREAD2=0
      WRITE(6,'(a)') '  -- Select one option to define the material:'
      WRITE(6,'(a)') '     1: Enter composition data from keyboard'
      WRITE(6,'(a)') '     2: Read from PENDBASE_photons/pdcompos.p06'
      READ(5,*) IREAD

      IF(IREAD.EQ.1) THEN
C  ****  Keyboard input.
        WRITE(6,'(a)') '  -- Enter material name (up to 60 characters):'
        READ(5,'(A62)') NAME

        WRITE(6,'(a)') '  -- Number of elements in the molecule:'
        READ(5,*) NELEM(M)
        IF(NELEM(M).GT.30.OR.NELEM(M).LT.1) THEN
          WRITE(6,*) ' STOP. NELEM must be between 1 and 30.'
          STOP
        ENDIF

        IF(NELEM(M).EQ.1) THEN
          WRITE(6,'(a)') '  -- Atomic number of the element (Z):'
          READ(5,*) IZZ
          IZ(M,1)=IZZ
          STF(M,1)=1.0D0
        ELSE
          WRITE(6,'(a)') '  -- Select composition format:'
          WRITE(6,'(a)') '     1: Stoichiometric formula (atoms/molec.)'
          WRITE(6,'(a)') '     2: Fractions by weight'
          READ(5,*) IREAD2
          IF(IREAD2.EQ.2) THEN
            WRITE(6,'(a)') '  -- Element Z and weight fraction:'
            DO I=1,NELEM(M)
              WRITE(6,'(''     Element'',I3,'':'')') I
              READ(5,*) IZ(M,I),FBW(I)
            ENDDO
          ELSE
            WRITE(6,'(a)')
     &        '  -- For each element: Z and atoms per molecule'
            DO I=1,NELEM(M)
              WRITE(6,'(''     Element'',I3,'':'')') I
              READ(5,*) IZ(M,I),STF(M,I)
            ENDDO
          ENDIF
        ENDIF

C  ****  Compute atomic weights and stoichiometric fractions.
        WRITE(6,'(a)') '  -- Enter mass density (g/cm**3):'
        READ(5,*) RHO(M)

      ELSE
C  ****  Read composition from pdcompos.p06 database.
        WRITE(6,'(a)') '  -- Enter material identification number'
        WRITE(6,'(a)') '     (1-99 = elements, 100-280 = compounds):'
        READ(5,*) IDNUM
        IF(IDNUM.LT.1.OR.IDNUM.GT.280) THEN
          WRITE(6,*) ' STOP. Material ID must be between 1 and 280.'
          STOP
        ENDIF

        OPEN(3,FILE='PENDBASE_photons/pdcompos.p06')
        DO I=1,15
          READ(3,'(A62)') NAME
        ENDDO
        DO K1=1,300
          READ(3,'(I3,2X,A62)',END=902) IORD,NAME
          READ(3,*) NELEM(M),HOLLOW,EXPOT0,RHO(M)
          IF(NELEM(M).GT.30.OR.NELEM(M).LT.1) THEN
            WRITE(6,*) ' STOP. Corrupt pdcompos.p06 file.'
            STOP
          ENDIF
          DO I=1,NELEM(M)
            READ(3,*) IZ(M,I),HOLLOW,STF(M,I)
          ENDDO
          IF(IORD.EQ.IDNUM) GO TO 903
        ENDDO
  902   CONTINUE
        WRITE(6,*) ' STOP. Material not found in pdcompos.p06.'
        STOP
  903   CONTINUE
        CLOSE(UNIT=3)
        WRITE(6,'(/2X,I3,1X,A62)') IORD,NAME
        WRITE(6,'(2X,''  Density = '',1P,E13.6,'' g/cm**3'')')RHO(M)
      ENDIF

C  ****  Handle weight fractions: convert to stoichiometric indices.
      IF(NELEM(M).GT.1.AND.IREAD.EQ.1.AND.IREAD2.EQ.2) THEN
        DO I=1,NELEM(M)
          CALL PEATWT(IZ(M,I),ATWTMP)
          STF(M,I)=FBW(I)/ATWTMP
        ENDDO
      ENDIF

C  ****  Ask whether material should be treated as insulator.
C        Conductors (metals, semiconductors): auto-plasmon for outer
C        electrons. Insulators (polymers, tissue, glass): Fcb=Wcb=0,
C        keeping electrons at their actual binding energies.
      WRITE(6,'(a)') '  -- Is this material an insulator?'
      WRITE(6,'(a)')
     &  '     1: yes, Fcb=Wcb=0 (polymers, tissue, glass, ...)'
      WRITE(6,'(a)')
     &  '     2: no, keep auto-plasmon (metals, semiconductors)'
      READ(5,*) IFINS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  MC-GPU energy table parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') '  -- Enter energy range: Emin  Emax (eV):'
      WRITE(6,'(a)') '     e.g.:  5000  120000'
      READ(5,*) EMIN, EMAX
      WRITE(6,'(a)') '  -- Enter energy step DE (eV), e.g.: 5'
      WRITE(6,'(a)')
     & '     (one extra bin above Emax will be added for'//
     & ' interpolation)'
      READ(5,*) DE
C  ****  Compute number of bins: covers [Emin, Emax] plus one extra step.
C        NINT handles any floating-point imprecision in (Emax-Emin)/DE.
      NBINS = NINT((EMAX-EMIN)/DE) + 2
      EMAX  = EMIN + DBLE(NBINS-1)*DE   ! actual table maximum
      WRITE(6,'(a,I6,a,1pe12.5,a)')
     &  '     ',NBINS,' bins  --  table Emax = ',EMAX,' eV'
      WRITE(6,'(a)') '  -- Enter output file name (e.g., water.mcgpu):'
      READ(5,'(A80)') OUTNAME

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  PART B: Photon-only material initialization
C  (replaces PEINIT + PEMATR from PENELOPE, photon interactions only)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  ****  Set photon absorption energy = Emin so EGRID lower bound is set.
      EABS(2,M) = EMIN
      EABS(1,M) = EMIN
      EABS(3,M) = EMIN

      WRITE(6,*) ' '
      WRITE(6,*) ' Processing material data. Please wait...'
      CALL PEINIT_PHOTONS(EMIN, EMAX, M, 6, 1, IFINS)
      WRITE(6,*) ' '

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  PART C: MFP table computation and .mcgpu file output
C  (adapted from MC-GPU_create_material_data.f, unchanged logic)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  ****  Compute Rayleigh cumulative probability on linear energy grid.
      CALL GRAaI_linear_energy(M, NBINS, EMIN, DE, PMAX_linear_energy)

C  ****  Open output file and write header.
      OPEN(1, FILE=OUTNAME)
      WRITE(1,'(a)')
     & '#[MATERIAL DEFINITION FOR MC-GPU: interaction'//
     & ' mean free path and sampling data from PENELOPE 2006]'
      WRITE(1,'(a)') '#[MATERIAL NAME]'
      WRITE(1,1001) NAME
 1001 FORMAT('# ',a)
      WRITE(1,'(a)') '#[NOMINAL DENSITY (g/cm^3)]'
      WRITE(1,1002) RHO(M)
 1002 FORMAT('# ',f12.8)
      WRITE(1,'(a)') '#[NUMBER OF DATA VALUES]'
      WRITE(1,1003) NBINS
 1003 FORMAT('# ',I6)
      WRITE(1,'(a)')
     & '#[MEAN FREE PATHS (cm)'//
     & ' (ie, average distance between interactions)]'
      WRITE(1,'(a)')
     & '#[Energy (eV)   | Rayleigh   |'//
     & ' Compton   | Photoelectric |'//
     & ' TOTAL (+pair prod) (cm) |'//
     & ' Rayleigh: max cumul prob F^2]'

C  ****  Compute MFP table: one row per linear energy bin.
C        Pair production is included in the total (ICOL=4) to give
C        correct photon beam attenuation at radiotherapy energies.
C        Below 1022 keV PHMFP returns ~1D35 cm (effectively infinite),
C        so the pair production term contributes zero at diagnostic range.
      KPAR=2    ! photons
      DO I=1,NBINS
        E = EMIN + (I-1)*DE
        E_MFP(1) = E
        E_MFP(2) = PHMFP(E,KPAR,M,1)    ! Rayleigh MFP (cm)
        E_MFP(3) = PHMFP(E,KPAR,M,2)    ! Compton MFP (cm)
        E_MFP(4) = PHMFP(E,KPAR,M,3)    ! Photoelectric MFP (cm)
        E_MFP(5) = PHMFP(E,KPAR,M,4)    ! Pair production MFP (cm)
        E_MFP(6) = 1.0D0/( 1.0D0/E_MFP(2) + 1.0D0/E_MFP(3)
     &                   + 1.0D0/E_MFP(4) + 1.0D0/E_MFP(5) )
        WRITE(1,'(6(1x,1pe12.5))')
     &    E_MFP(1),E_MFP(2),E_MFP(3),E_MFP(4),E_MFP(6),
     &    PMAX_linear_energy(I)
      ENDDO

C  ****  Write Rayleigh interaction data (RITA sampling of form factors).
      WRITE(1,'(a)')
     & '#[RAYLEIGH INTERACTIONS (RITA sampling'//
     & '  of atomic form factor from EPDL database)]'
      WRITE(1,'(a)')
     & '#[DATA VALUES TO SAMPLE SQUARED MOLECULAR FORM FACTOR (F^2)]'
      WRITE(1,1003) NP
      WRITE(1,'(a)')
     & '#[SAMPLING DATA FROM COMMON/CGRA/: X, P, A, B, ITL, ITU]'
      DO I=1,NP
        WRITE(1,5555) XCO(I,M),PCO(I,M),ACO(I,M),
     1    BCO(I,M),ITLCO(I,M),ITUCO(I,M)
      ENDDO
 5555 FORMAT(4(1x,1pe12.5),1x,i4,1x,i4)

C  ****  Write Compton interaction data (impulse approximation).
      WRITE(1,'(a)')
     & '#[COMPTON INTERACTIONS (relativistic impulse model with'//
     & '  approximated one-electron analytical profiles)]'
      WRITE(1,'(a)') '#[NUMBER OF SHELLS]'
      WRITE(1,1003) NOSCCO(M)
      WRITE(1,'(a)')
     & '#[SHELL INFORMATION FROM COMMON/CGCO/:'//
     & ' FCO, UICO, FJ0, KZCO, KSCO]'
      DO I=1,NOSCCO(M)
        WRITE(1,5107) FCO(M,I),UICO(M,I),FJ0(M,I),
     &    KZCO(M,I),KSCO(M,I)
      ENDDO
 5107 FORMAT(3(1X,1pe12.5),2(1X,I4))

      WRITE(1,'(a)') ' '
      CLOSE(1)

      WRITE(6,'(a)')
     & '*** Material file correctly generated. Have a nice simulation!'
      WRITE(6,*) ' '

      END


C  *********************************************************************
C                  SUBROUTINE PEINIT_PHOTONS
C  *********************************************************************
      SUBROUTINE PEINIT_PHOTONS(EMIN, EMAX, M, IWR, INFO, IFINS)
C
C  Photon-only material initialization for MC-GPU.
C  Replaces PEINIT + PEMATR from PENELOPE 2006, keeping only the
C  subroutines needed for photon transport (Rayleigh, Compton,
C  photoelectric). Electron and positron physics are omitted.
C
C  Reads from PENDBASE_photons/:
C    pdatconf.p06  -- atomic shell structure and Compton profiles
C    pdgph##.p06   -- photoelectric cross sections per element
C
C  Populates COMMON blocks used by PHMFP and GRAaI_linear_energy:
C    CEGRID  -- PENELOPE log energy grid (NEGP=200 points)
C    COMPOS  -- material composition (ZT, AT, VMOL)
C    CADATA  -- element shell data (EB, IFI, IKS, NSHT) [via pdatconf]
C    CGCO    -- Compton shell groups (FCO, UICO, FJ0, KZCO, KSCO)
C    CGPH00  -- photoelectric tables (EPH, XPH, IPHF, IPHL)
C    CGIMFP  -- photon cross sections on log grid (SGRA, SGCO)
C    CGRA    -- Rayleigh RITA sampling grid (XCO, PCO, ACO, BCO)
C
C  Based on PEMATW (material file writer) and PEMATR (reader) from
C  PENELOPE 2006 (F. Salvat, J. M. Fernandez-Varea and J. Sempau).
C  Compton shell grouping: PEMATW lines 1539-2000.
C  Photon cross section setup: PEMATR lines 971-1038.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      CHARACTER*5 CH5
      PARAMETER (A0B=5.291772108D-9)    ! Bohr radius (cm)
      PARAMETER (HREV=27.2113845D0)     ! Hartree energy (eV)
      PARAMETER (AVOG=6.0221415D23)     ! Avogadro's number
      PARAMETER (REV=5.10998918D5)      ! Electron rest energy (eV)
      PARAMETER (SL=137.03599911D0)     ! Speed of light (1/alpha)
      PARAMETER (PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Element data (RA1-RA5, RSCR, ATW etc. are set via BLOCK DATA
C        PENDAT in penelope.f -- no file read needed for Rayleigh).
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
C  ****  Simulation parameters.
      COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
     1  WCR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  E/P inelastic collisions (declared for oscillator build step;
C        not used in MC-GPU photon transport).
      PARAMETER (NO=64)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Compton scattering.
      PARAMETER (NOCO=64)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     2  KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),NOSCCO(MAXMAT)
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Local oscillator arrays (from PEMATW, used to build CGCO).
      PARAMETER (NOM=400)
      DIMENSION FF(NOM),UUI(NOM),WWRI(NOM),FFJ0(NOM),KKZ(NOM),KKS(NOM)
      DIMENSION FC(NOM),UIC(NOM),FJ0C(NOM),KZC(NOM),KSC(NOM)
      DIMENSION FFT(NOM),UIT(NOM),WRIT(NOM),KZT(NOM),KST(NOM)
C  ****  Local Compton profile array (read from pdatconf.p06).
      DIMENSION CP(99,30)

C  ****  Step 1: Set up PENELOPE log energy grid.
C        EGRID sets EL=0.99999*EMIN and EU=1.00001*EMAX, so that
C        the first tabulation energy E=EMIN satisfies E > EL.
      CALL EGRID(EMIN, EMAX)

C  ****  Step 2: Initialize photoelectric and pair production tables.
      CALL GPHa0
      CALL GPPa0(M)

C  ****  Step 3: Compute AT, ZT, VMOL from user-supplied IZ, STF, NELEM.
C        (STF, IZ, NELEM, RHO already set by Part A user interaction.)
      ZT(M)=0.0D0
      AT(M)=0.0D0
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        ZT(M)=ZT(M)+IZZ*STF(M,I)
        AT(M)=AT(M)+ATW(IZZ)*STF(M,I)
      ENDDO
      VMOL(M)=AVOG*RHO(M)/AT(M)

C
C  ****  Steps 4+5: Read pdatconf.p06 and build oscillator / Compton data.
C        Code adapted from PEMATW (penelope.f lines 1480-2000).
C        The oscillator table (CEIN) is built as an intermediate step
C        since the Compton shell grouping derives from it.
C

C  ----  Initialise CADATA shell data arrays for this material's elements.
      DO I=1,99
        NSHT(I)=0
        DO J=1,30
          EB(I,J)=0.0D0
          CP(I,J)=0.0D0
          IFI(I,J)=0
          IKS(I,J)=0
        ENDDO
      ENDDO

C  ----  Read pdatconf.p06 for each element (once per unique Z).
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
C  ****  NSHT(IZZ).EQ.0 means element data not yet loaded.
        IF(NSHT(IZZ).EQ.0) THEN
          OPEN(3,FILE='PENDBASE_photons/pdatconf.p06')
          DO J=1,19
            READ(3,'(A5)') CH5
          ENDDO
          NS=0
          IZZT=0
          DO J=1,150000
            READ(3,2005,END=905) IIZ,IS,CH5,IIF,IE,CCP
 2005       FORMAT(I3,1X,I2,1X,A5,1X,I1,1X,I6,E9.2)
            IF(IIZ.EQ.IZZ) THEN
              NS=NS+1
              IF(NS.GT.30) THEN
                WRITE(6,'(/1X,''NS ='',I4)') NS
                WRITE(6,*) ' STOP. Too many shells.'
                STOP
              ENDIF
              IF(IS.LT.1.OR.IS.GT.30) THEN
                WRITE(6,'(/1X,''IS ='',I4)') IS
                WRITE(6,*) ' STOP. Wrong shell number.'
                STOP
              ENDIF
              IZZT=IZZT+IIF
              EB(IZZ,IS)=IE           ! shell binding energy (eV)
              CP(IZZ,IS)=CCP          ! one-electron Compton profile
              IFI(IZZ,IS)=IIF         ! electrons per shell
              IKS(IZZ,NS)=IS          ! shell index mapping
            ENDIF
          ENDDO
  905     CONTINUE
          NSHT(IZZ)=NS
          IF(IZZ.NE.IZZT) THEN
            WRITE(6,*) ' STOP. Unbalanced charges (element Z=',IZZ,')'
            STOP
          ENDIF
          CLOSE(3)
        ENDIF
      ENDDO

C  ----  Build oscillator table from shell data (PEMATW lines 1543-1626).
      DO I=1,NO
        F(M,I)=0.0D0
        UI(M,I)=0.0D0
        WRI(M,I)=0.0D0
        KZ(M,I)=0
        KS(M,I)=0
      ENDDO
      DO I=1,NOM
        FF(I)=0.0D0
        UUI(I)=0.0D0
        WWRI(I)=0.0D0
        FFJ0(I)=0.0D0
        KKZ(I)=0
        KKS(I)=0
      ENDDO
      FT=0.0D0
C  ****  1st oscillator = conduction band (binding energy < 13 eV).
      NOS=1
      FF(1)=0.0D0
      UUI(1)=0.0D0
      FFJ0(1)=0.0D0
      KKZ(1)=0
      KKS(1)=30
      DO I=1,NELEM(M)
        IZZ=IZ(M,I)
        DO K=1,30
          JS=IKS(IZZ,K)
          IF(JS.GT.0) THEN
            IF(IFI(IZZ,JS).GT.0) THEN
              NOS=NOS+1
              IF(NOS.GT.NOM) THEN
                WRITE(6,*) ' STOP. Too many oscillators (NOM).'
                STOP
              ENDIF
              FF(NOS)=IFI(IZZ,JS)*STF(M,I)
              UUI(NOS)=EB(IZZ,JS)
              FFJ0(NOS)=CP(IZZ,JS)
              KKZ(NOS)=IZZ
              IF(IZZ.GT.6.AND.JS.LT.10) THEN
                KKS(NOS)=JS
              ELSE
                KKS(NOS)=30
              ENDIF
              FT=FT+FF(NOS)
              IF(EB(IZZ,K).LT.13.0D0)
     1          FF(1)=FF(1)+IFI(IZZ,JS)*STF(M,I)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF(ABS(FT-ZT(M)).GT.1.0D-10*ZT(M)) THEN
        WRITE(6,*) ' STOP. Unbalanced charges (compound).'
        STOP
      ENDIF
C  ****  Sort oscillators by increasing binding energy.
      DO I=1,NOS-1
        DO J=I+1,NOS
          IF(UUI(I).GE.UUI(J)) THEN
            SAVE=UUI(I)
            UUI(I)=UUI(J)
            UUI(J)=SAVE
            SAVE=FF(I)
            FF(I)=FF(J)
            FF(J)=SAVE
            SAVE=FFJ0(I)
            FFJ0(I)=FFJ0(J)
            FFJ0(J)=SAVE
            ISAVE=KKZ(I)
            KKZ(I)=KKZ(J)
            KKZ(J)=ISAVE
            ISAVE=KKS(I)
            KKS(I)=KKS(J)
            KKS(J)=ISAVE
          ENDIF
        ENDDO
      ENDDO

C  ----  Conduction band / plasmon setup (PEMATW lines 1628-1716).
C        Uses default values: FP=FF(1), EP=EPP (equivalent to IPLOSP=2).
      OP2(M)=FOURPI*ZT(M)*VMOL(M)*A0B**3*HREV**2
      OMEGA=SQRT(OP2(M))
      EPP=OMEGA*SQRT(FF(1)/ZT(M))
      EP=EPP
      FP=FF(1)
C  ****  Apply insulator/conductor choice from user (IFINS argument).
C        If IFINS=1 (insulator) or FF(1) is negligible: set Fcb=Wcb=0.
      IF(IFINS.EQ.1) THEN
        FP=0.0D0
        EP=0.0D0
      ENDIF
      IF(INFO.GE.1) THEN
        IF(EP.LT.1.0D0.OR.FP.LT.0.5D0) THEN
          WRITE(IWR,'(1X,A)')
     1      '  Insulator mode: Fcb=Wcb=0'
        ELSE
          WRITE(IWR,'(1X,A,1PE12.5)')
     1      '  Conductor mode: Fcb=',FP
        ENDIF
      ENDIF

      IF(EP.LT.1.0D0.OR.FP.LT.0.5D0) THEN
C  ****  Insulator: remove conduction band oscillator.
        DO J=1,NOS-1
          FF(J)=FF(J+1)
          UUI(J)=UUI(J+1)
          FFJ0(J)=FFJ0(J+1)
          KKZ(J)=KKZ(J+1)
          KKS(J)=KKS(J+1)
        ENDDO
        NOS=NOS-1
      ELSE
C  ****  Conductor: consolidate outer shells into conduction band.
        IDEAD=0
        FPP=FP
        I=1
  907   I=I+1
        IF(FF(I).LT.FPP) THEN
          FPP=FPP-FF(I)
          FF(I)=0.0D0
          IDEAD=IDEAD+1
          GO TO 907
        ELSE
          FF(I)=FF(I)-FPP
          IF(ABS(FF(I)).LT.1.0D-12) THEN
            FP=FP+FF(I)
            FF(I)=0.0D0
            IDEAD=IDEAD+1
          ENDIF
        ENDIF
        FF(1)=FP
        UUI(1)=0.0D0
        WWRI(1)=EP
        FFJ0(1)=0.75D0/SQRT(3.0D0*PI*PI*VMOL(M)*A0B**3*FP)
        KKZ(1)=0
        KKS(1)=30
        IF(IDEAD.GT.0) THEN
          DO J=2,NOS-IDEAD
            FF(J)=FF(J+IDEAD)
            UUI(J)=UUI(J+IDEAD)
            FFJ0(J)=FFJ0(J+IDEAD)
            KKZ(J)=KKZ(J+IDEAD)
            KKS(J)=KKS(J+IDEAD)
          ENDDO
          NOS=NOS-IDEAD
        ENDIF
      ENDIF
C  ****  Check f-sum rule.
      SUM=0.0D0
      DO J=1,NOS
        SUM=SUM+FF(J)
      ENDDO
      IF(ABS(SUM-ZT(M)).GT.1.0D-6*ZT(M)) THEN
        WRITE(6,*) ' STOP. Inconsistent oscillator strength data.'
        STOP
      ENDIF
      IF(ABS(SUM-ZT(M)).GT.1.0D-12*ZT(M)) THEN
        FACT=ZT(M)/SUM
        DO J=1,NOS
          FF(J)=FACT*FF(J)
        ENDDO
      ENDIF

C  ----  Copy oscillator arrays to Compton working arrays
C        (PEMATW line 1739; done BEFORE any Sternheimer modification
C         of FF/UUI, so Compton parameters are free of e/p corrections).
      NOSTC=NOS
      CSUMT=0.0D0
      DO I=1,NOSTC
        FC(I)=FF(I)
        UIC(I)=UUI(I)
        FJ0C(I)=FFJ0(I)
        KZC(I)=KKZ(I)
        KSC(I)=KKS(I)
        CSUMT=CSUMT+FC(I)*FJ0C(I)
      ENDDO

C  ----  Compton shell grouping (PEMATW lines 1901-1999).
C        Groups outer shells with similar binding energies; inner shells
C        (K, L, M with E > WISCUT) are kept separate for accuracy.
      IZMAX=0
      DO I=1,NELEM(M)
        IZMAX=MAX(IZMAX,IZ(M,I))
      ENDDO
      WISCUT=MAX(200.0D0,EB(IZMAX,10))

      RGROUP=1.50D0
C  ****  Sort by increasing ionisation energy.
      IF(NOSTC.GT.1) THEN
        DO I=1,NOSTC-1
          DO J=I+1,NOSTC
            IF(UIC(I).GT.UIC(J)) THEN
              SAVE=FC(I)
              FC(I)=FC(J)
              FC(J)=SAVE
              SAVE=UIC(I)
              UIC(I)=UIC(J)
              UIC(J)=SAVE
              SAVE=FJ0C(I)
              FJ0C(I)=FJ0C(J)
              FJ0C(J)=SAVE
              ISAVE=KZC(I)
              KZC(I)=KZC(J)
              KZC(J)=ISAVE
              ISAVE=KSC(I)
              KSC(I)=KSC(J)
              KSC(J)=ISAVE
            ENDIF
          ENDDO
        ENDDO
      ENDIF

  920 CONTINUE
      IELIM=0
      IF(NOSTC.GT.2) THEN
        DO 921 I=1,NOSTC-1
        IF(UIC(I).GT.WISCUT) GO TO 921
        IF(UIC(I).LT.1.0D0.OR.UIC(I+1).LT.1.0D0) GO TO 921
        IF(UIC(I+1).GT.RGROUP*UIC(I)) GO TO 921
        UIC(I)=(FC(I)*UIC(I)+FC(I+1)*UIC(I+1))/(FC(I)+FC(I+1))
        FJ0C(I)=(FC(I)*FJ0C(I)+FC(I+1)*FJ0C(I+1))/(FC(I)+FC(I+1))
        FC(I)=FC(I)+FC(I+1)
        IF(KZC(I).NE.KZC(I+1)) KZC(I)=0
        KSC(I)=30
        IF(I.LT.NOSTC-1) THEN
          DO J=I+1,NOSTC-1
            FC(J)=FC(J+1)
            UIC(J)=UIC(J+1)
            FJ0C(J)=FJ0C(J+1)
            KZC(J)=KZC(J+1)
            KSC(J)=KSC(J+1)
          ENDDO
        ENDIF
        IELIM=IELIM+1
        FC(NOSTC)=0.0D0
        UIC(NOSTC)=0.0D0
        FJ0C(NOSTC)=0.0D0
        KZC(NOSTC)=0
        KSC(NOSTC)=0
  921   CONTINUE
      ENDIF
      IF(IELIM.GT.0) THEN
        NOSTC=NOSTC-IELIM
        GO TO 920
      ENDIF

C  ****  Transfer grouped Compton shells to CGCO common block.
      IF(NOSTC.LT.NOCO) THEN
        NOSCCO(M)=NOSTC
        DO I=1,NOSCCO(M)
          FCO(M,I)=FC(I)
          UICO(M,I)=UIC(I)
          FJ0(M,I)=FJ0C(I)
          KZCO(M,I)=KZC(I)
          KSCO(M,I)=KSC(I)
          CSUMT=CSUMT-FCO(M,I)*FJ0(M,I)
        ENDDO
        IF(ABS(CSUMT).GT.1.0D-9) THEN
          WRITE(6,*) ' STOP. Error in grouping the Compton profiles.'
          WRITE(6,'(''  Residual sum ='',1p,e12.5)') ABS(CSUMT)
          STOP
        ENDIF
      ELSE
        RGROUP=RGROUP**2
        GO TO 920
      ENDIF
      WRITE(IWR,'(1X,'' Compton grouping factor = '',1P,E12.5)') RGROUP
      WRITE(IWR,'(1X,'' Number of Compton shells  = '',I3)') NOSCCO(M)

C  ****  Scale FJ0 by SL = 1/alpha (speed of light in a.u.).
C        Both PEMATW (line 2040) and PEMATR (line 440) apply this scaling
C        after transferring FJ0 to the CGCO common block.
      DO I=1,NOSCCO(M)
        FJ0(M,I)=FJ0(M,I)*SL
      ENDDO

C  ****  Step 6: Read pdgph##.p06 and populate photoelectric tables.
      CALL GPHaDirect(M, IWR, INFO)

C  ****  Step 6b: Read pdgpp##.p06 and populate pair production table.
C        SGPP(M,KE) is used by PHMFP(E,2,M,4) to give correct total
C        attenuation including pair production above 1.022 MeV.
      CALL GPPaDirect(M, IWR, INFO)

C  ****  Steps 7+8: Compute Rayleigh and Compton log cross sections
C        on PENELOPE's 200-point log energy grid (SGRA, SGCO arrays).
C        Uses GRAaT and GCOaT directly at each grid point
C        (equivalent to PEMATR's spline-interpolation of pre-computed CS).
      VMOLL=LOG(VMOL(M))
      DO I=1,NEGP
        EE=ET(I)
        CALL GRAaT(EE,CSR,M)
        IF(CSR.LT.1.0D-35) CSR=1.0D-35
        SGRA(M,I)=LOG(CSR)+VMOLL

        CALL GCOaT(EE,CSC,M)
        IF(CSC.LT.1.0D-35) CSC=1.0D-35
        SGCO(M,I)=LOG(CSC)+VMOLL
      ENDDO

C  ****  Step 9: Initialise Rayleigh RITA sampling grid.
      CALL GRAaI(M)

      RETURN
      END


C  *********************************************************************
C                  SUBROUTINE GPHaDirect
C  *********************************************************************
      SUBROUTINE GPHaDirect(M, IWR, INFO)
C
C  Reads photoelectric cross sections directly from the PENDBASE_photons
C  database files (pdgph##.p06) and populates the CGPH00 common block
C  (EPH, XPH, IPHF, IPHL, NPHS, NCUR) used by PHMFP.
C
C  Combines the reading logic of GPHaW (penelope.f lines 6479-6511)
C  with the CGPH00 population logic of GPHaR (lines 6325-6346).
C  Only the total cross section column (XPH(*,1)) is stored because
C  PHMFP's photoelectric branch (ICOL=3) uses only the total.
C  Shell-resolved columns and SGPH are not needed for MC-GPU.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Photoelectric cross sections.
      PARAMETER (NTP=8000)
      COMMON/CGPH00/EPH(NTP),XPH(NTP,10),IPHF(99),IPHL(99),NPHS(99),
     1  NCUR
C  ****  Local work arrays (raw data from database file).
      PARAMETER (NPHM=400)
      DIMENSION XS(10),E0(NPHM),XS0(NPHM)

      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
C  ****  Build file name pdgph##.p06 using element atomic number IZZ.
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdgph'//LDIG2//LDIG1//'.p06'
C  ****  Read directly from PENDBASE_photons/.
        OPEN(3,FILE='PENDBASE_photons/'//FILEN)
        READ(3,*) IZZZ,NSHR
        IF(IZZZ.NE.IZZ) THEN
          WRITE(6,*) ' STOP. GPHaDirect: corrupt file ',FILEN
          STOP
        ENDIF
        IF(NSHR.GT.9) THEN
          WRITE(6,*) ' STOP. GPHaDirect: too many shells in ',FILEN
          STOP
        ENDIF
        NPTAB=0
        DO IE=1,1000
          READ(3,*,END=1) ER,(XS(IS),IS=1,NSHR+1)
          IF(ER.GT.49.9D0.AND.ER.LT.1.01D9) THEN
            NPTAB=NPTAB+1
            IF(NPTAB.GT.NPHM) THEN
              WRITE(6,*) ' STOP. GPHaDirect: NPHM too small.'
              STOP
            ENDIF
            E0(NPTAB)=ER
C  ****  Total cross section is column 1 (XS(1), in barn). Convert to cm^2.
            XS0(NPTAB)=XS(1)*1.0D-24
          ENDIF
        ENDDO
    1   CONTINUE
        CLOSE(3)
        IF(INFO.GE.2) THEN
          WRITE(IWR,2001) IZZ,NSHR,NPTAB
 2001     FORMAT(/1X,'***  Photoelectric cross sections,  IZ =',I3,
     1      ',  NSHELL =',I3,',  NDATA =',I4)
        ENDIF

C  ****  Store in CGPH00 (only once per element, identified by NPHS(IZZ)).
        IF(NPHS(IZZ).EQ.0) THEN
          IPHF(IZZ)=NCUR+1
          IF(NCUR+NPTAB.GT.NTP) THEN
            WRITE(6,*) ' STOP. GPHaDirect: increase NTP to',NCUR+NPTAB
            STOP
          ENDIF
          DO IE=1,NPTAB
            IC=NCUR+IE
            EPH(IC)=LOG(E0(IE))
            XPH(IC,1)=LOG(MAX(XS0(IE),1.0D-35))
          ENDDO
          NCUR=NCUR+NPTAB
          IPHL(IZZ)=NCUR
          NPHS(IZZ)=NSHR
        ENDIF
      ENDDO

      RETURN
      END


C  *********************************************************************
C                  SUBROUTINE GPPaDirect
C  *********************************************************************
      SUBROUTINE GPPaDirect(M, IWR, INFO)
C
C  Reads pair production cross sections directly from PENDBASE_photons
C  (pdgpp##.p06) and populates SGPP(M,KE) on PENELOPE's log energy grid.
C  SGPP is used by PHMFP(E,2,M,4) to include pair production in the
C  total photon attenuation, giving correct primary beam fluence at
C  radiotherapy energies (>1 MeV). Below threshold PHMFP returns ~1D35
C  so the pair production term is negligible at diagnostic energies.
C
C  Logic adapted from PEMATW/GPPaW (penelope.f lines 6677-6725) for
C  reading pdgpp files, and PEMATR (lines 1016-1034) for SGPP setup.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*12 FILEN
      CHARACTER*1 LDIG(10),LDIG1,LDIG2
      DATA LDIG/'0','1','2','3','4','5','6','7','8','9'/
C  ****  Composition data.
      PARAMETER (MAXMAT=10)
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      PARAMETER (NEGP=200)
      COMMON/CEGRID/EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)
C  ****  Local arrays for pair production data.
C        NPPMAX must be >= entries in pdgpp files (typically ~88).
      PARAMETER (NPPMAX=1500)
      DIMENSION EIT(NPPMAX),XG0(NPPMAX)
      DIMENSION F4(NPPMAX),FL(NPPMAX)
      DIMENSION A(NEGP),B(NEGP),C(NEGP),D(NEGP)

C  ****  Initialise combined cross section accumulator.
      DO I=1,NPPMAX
        XG0(I)=0.0D0
      ENDDO

C  ****  Read pdgpp files and accumulate total cross section per molecule.
C        All pdgpp files share the same energy grid (EIT is overwritten
C        each element but values are identical; XG0 accumulates).
      NPTAB=0
      DO IEL=1,NELEM(M)
        IZZ=IZ(M,IEL)
        WGHT=STF(M,IEL)*1.0D-24    ! stoichiometry x (barn -> cm^2)
        NLD=IZZ
        NLD1=NLD-10*(NLD/10)
        NLD2=(NLD-NLD1)/10
        LDIG1=LDIG(NLD1+1)
        LDIG2=LDIG(NLD2+1)
        FILEN='pdgpp'//LDIG2//LDIG1//'.p06'
        OPEN(3,FILE='PENDBASE_photons/'//FILEN)
        READ(3,*) IZZZ
        IF(IZZZ.NE.IZZ) THEN
          WRITE(IWR,*) ' STOP. GPPaDirect: corrupt file ',FILEN
          STOP
        ENDIF
        DO I=1,NPPMAX
          READ(3,*,END=1) EIT(I),XG0P
          XG0(I)=XG0(I)+WGHT*XG0P
          NPTAB=I
          IF(EIT(I).GT.0.999D9) GO TO 1
        ENDDO
    1   CONTINUE
        CLOSE(3)
      ENDDO

C  ****  Build log-log spline of pair production CS above 1.022 MeV.
C        Same logic as PEMATR lines 1016-1034.
      VMOLL=LOG(VMOL(M))
      NP=0
      DO 2 I=1,NPTAB
        IF(EIT(I).LT.1.023D6) GO TO 2
        NP=NP+1
        F4(NP)=LOG(EIT(I))
        FL(NP)=LOG(MAX(XG0(I),1.0D-35))
    2 CONTINUE

      IF(NP.GT.1) THEN
        CALL SPLINE(F4,FL,A,B,C,D,0.0D0,0.0D0,NP)
        DO I=1,NEGP
          IF(ET(I).LT.1.023D6) THEN
            SGPP(M,I)=-80.6D0
          ELSE
            EC=DLEMP(I)
            CALL FINDI(F4,EC,NP,J)
            SGPP(M,I)=A(J)+EC*(B(J)+EC*(C(J)+EC*D(J)))+VMOLL
          ENDIF
        ENDDO
      ELSE
C  ****  No data above threshold (e.g. very low-Z material + low Emax).
        DO I=1,NEGP
          SGPP(M,I)=-80.6D0
        ENDDO
      ENDIF

      IF(INFO.GE.2) THEN
        WRITE(IWR,'(1X,A,I4,A)')
     1    '  Pair production: read ',NPTAB,' energy points'
      ENDIF

      RETURN
      END


C  *********************************************************************
C     SUBROUTINE GRAaI_linear_energy
C     (verbatim from MC-GPU_create_material_data.f, A. Badal 2009)
C     with PENDBASE_photons/ path -- not actually needed since
C     Rayleigh form factors are analytical (BLOCK DATA PENDAT).
C  *********************************************************************
      SUBROUTINE GRAaI_linear_energy(M, nbins, emin, de, PMAX_linear_e)
C
C  Re-init random sampling for Rayleigh scattering using the input
C  linear energy scale. Computes the maximum cumulative probability
C  PMAX_linear_e(IE) for each linear energy bin IE.
C
C  Based on PENELOPE's subroutine GRAaI (coherent Rayleigh scattering
C  initialisation). Adapted for the linear energy grid used by MC-GPU
C  (MC-GPU uses direct linear interpolation, not PENELOPE's log-log).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (REV=5.10998918D5)   ! Electron rest energy (eV)
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
c     COMMON/CGRA/XCO(NP,MAXMAT),PCO(NP,MAXMAT),ACO(NP,MAXMAT),
c    1  BCO(NP,MAXMAT),PMAX(NEGP,MAXMAT),ITLCO(NP,MAXMAT),
c    2  ITUCO(NP,MAXMAT)
C
      PARAMETER (NIP=51)
      DIMENSION XI(NIP),FUN(NIP),SUM(NIP)
      COMMON/CGRA00/FACTE,X2MAX,MM,MOM
      EXTERNAL GRAaD1

C  ****  Output array (linear energy grid, NBINS points).
      PARAMETER (MAX_ENERGY_BINS=50000)
      DIMENSION PMAX_linear_e(MAX_ENERGY_BINS)

      IZZ=0
      DO I=1,NELEM(M)
        IZZ=MAX(IZZ,IZ(M,I))
      ENDDO

      MM=M
      X2MIN=0.0D0
      X2MAX=4.0D0*20.6074D0**2*(200.0D0*IZZ)**2
      NPT=NP
      NU=NPT/4
      CALL RITAI0(GRAaD1,X2MIN,X2MAX,NPT,NU,ERRM,0)

C  ****  Upper limit of the X2 interval for linear energy bins.
      do IE = 1, nbins
        XM=2.0D0*20.6074D0*(emin+(IE-1)*de)*RREV

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
            PMAX_linear_e(IE) = PACI(I)+SUM(NIP)
          ELSE
            PMAX_linear_e(IE) = 1.0D0
          ENDIF
        ELSE
          PMAX_linear_e(IE) = PACI(1)
        ENDIF
      ENDDO

      RETURN
      END


C  *********************************************************************
C     SUBROUTINE PEATWT
C     Returns atomic weight for element IZZ (from BLOCK DATA PENDAT).
C  *********************************************************************
      SUBROUTINE PEATWT(IZZ, ATWT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      COMMON/CADATA/ATW(99),EPX(99),RA1(99),RA2(99),RA3(99),RA4(99),
     1  RA5(99),RSCR(99),ETA(99),EB(99,30),IFI(99,30),IKS(99,30),
     2  NSHT(99),LASYMB(99)
      ATWT=ATW(IZZ)
      RETURN
      END
