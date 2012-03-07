      SUBROUTINE Datain(IN,IOUT,NAtomax,NA,IC,Iopt,IP,IHam,
     1 IPR,IG,ISO1,ISO2,ISO3,IER,istop,IOPD,iupac,Ipent,IPH,
     1 IV1,IV2,IV3,icyl,PS,TolX,R5,R6,Rdist)
      IMPLICIT REAL*8 (A-H,O-Z)
      Character DATEN*80
      Namelist /Coord/ IC,NA,IP,IV1,IV2,IV3,TolR,R5,R6,icyl
      Namelist /Opt/ Iopt
      Namelist /Hamilton/ IHam,iupac
      Namelist /Isomers/ IPR,IPH,IStop
      Namelist /Graph/ IG,ISO1,ISO2,ISO3,PS
C Input send to output   
      WRITE(IOUT,100)
      DO 10 I=1,200
      READ(IN,'(A80)',END=11) DATEN
   10 WRITE(IOUT,60) DATEN
   11 WRITE(IOUT,101)
      REWIND IN

      IOPD=0    !  Switch for producing the nth dual
      iupac=1   !  Switch for producing the Iupac nomenclature
                !  iupac=0 just count Hamiltonian Cycles
      Ipent=0   !  Initial flag for Spriral pentagon input
      IER=0     !  Error flag
      Tol=0.33d0 ! Tolerance
      IP=0      !  Print option
      IPR=0     !  Print Isomers
      IPH=0     !  Printin Hamiltonian cycles for each isomer
      NA=60     !  Number of Atoms
      IC=1      !  Input for fullerene structure
      iupac=0   !  Print IUPAC numbers
      IV1=2     !  Eigenvector option for fullerene construction
      IV2=3     !  Eigenvector option for fullerene construction
      IV3=4     !  Eigenvector option for fullerene construction
      istop=0   !  Option for stopping after isomer list
      icyl=0    !  Option for producing input for ploting program CYLVIEW
      ISO1=0    !  Option for fullerene orientation for Schlegel projection
      ISO2=0    !  Option for fullerene orientation for Schlegel projection
      ISO3=0    !  Option for fullerene orientation for Schlegel projection
      PS=0.d0   !  For Schlegel projection
      TolR=0.d0 !  Tolerance for finding ring connections
      R=1.391d0 !  
      R5=1.455d0 ! Distance in 5-Ring
      R6=R       ! Distance in 6-Ring

C Now process namelist input
C     Old input    
C     Read(IN,*) NA,IC,Iopt,IP,IHam,IPR,IS,ISO2,ISO2,ISO3,PS,TolR
C     New input    
      Read(IN,'(A80)') DATEN
      WRITE(IOUT,60) DATEN
      WRITE(IOUT,101)
      Read(IN,nml=Coord)
      Read(IN,nml=Opt)
      Read(IN,nml=Hamilton)
      Read(IN,nml=Isomers)
      Read(IN,nml=Graph)

C Set Parameters
C  Check on number of atoms (vertices)
      NAtom=IABS(NAtom)
      if(NA.gt.NAtomax) WRITE(IOUT,102) NA
      if(NA.lt.20.or.NA.eq.22) then
      Write(IOUT,103) NA
      IER=1
      return
      endif
      IF (NA/2*2.ne.NA) then
      Write(IOUT,104) NA
      IER=1
      return
      endif

C     Setting minimum distance
      if(R6.ne.R.and.R6.gt.1.d0) then
      Rdist=R6
      WRITE(Iout,106) Rdist
      else
      Rdist=R
      WRITE(Iout,107) Rdist
      endif

C  Output list
      if(IP.gt.0) then
       WRITE(IOUT,105)
       IP=1
      endif
      if(IP.lt.0) then
       IP=0
      endif

C  Tolerance for finding 5- and 6-ring connections
      if(TolR.le.0.d0) then
      TolX=Tol
      else
      TolX=TolR*0.01d0
      endif

      if(IPR.le.0) then
        IPR=-1
      endif
      if(IPR.eq.1) then
        if(NA.lt.60) IPR=0
      endif
      if(IPR.ge.2) then
        IPR=0
      endif
      
   60 FORMAT(1X,A80)
  100 FORMAT(1X,80(1H-),/1X,'I N P U T ',/1X,5H0....,
     161H....1.........2.........3.........4.........5.........6......,
     214H...7.........8,/1X,39H123456789012345678901234567890123456789,
     341H01234567890123456789012345678901234567890,/)
  101 FORMAT(1X,80(1H-))
  102 FORMAT(1X,' Number of Atoms exceed allowed limit of ',I4,
     1 ' Increase Parameter natom')
  103 FORMAT(1x,'Fullerene with requested number of carbon atoms ',
     1 I4,' not possible')
  104 FORMAT(1x,'Fullerene with odd number of carbon atoms ',
     1 I4,' not possible')
  105 FORMAT(1x,'Larger output requested')
  106 Format(1X,'Minimum bond distance set to input value: ',F12.6)
  107 Format(1X,'Minimum bond distance set to default value ',
     1 'taken from C60 bond distance: ',F12.6)
      RETURN
      END
