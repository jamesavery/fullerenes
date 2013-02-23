      SUBROUTINE Datain(IN,IOUT,NAtomax,ICart,Iopt,IP,IHam,
     1 nohueckel,KE,IPR,IPRC,ISchlegel,ISO1,ISO2,ISO3,IER,istop,
     1 leap,IGCtrans,iupac,Ipent,IPH,kGC,lGC,IV1,IV2,IV3,
     1 ixyz,ichk,isonum,loop,mirror,ilp,ISW,IYF,IBF,nzeile,ifs,
     1 ipsphere,ndual,nosort,ispsearch,novolume,ihessian,isearch,
     1 iprinth,ndbconvert,ihamstore,nhamcyc,isomerl,isomerh,
     1 PS,TolX,R5,R6,Rdist,rvdwc,scale,scalePPG,ftol,scaleRad,jumps,
     1 force,forceP,boost,filename,filenameout,DATEN)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      integer NA,iopt
      real(8) force(ffmaxdim),forceP(ffmaxdim) ! user chosen FF (and a backup)
      integer endzeile
      integer jumps(10)
      Character*1 DATEN(nzeile)
      Character filename*50
      Character filenameout*50
      Namelist /General/ NA,IP,TolR,R5,R6,ixyz,ichk,
     1 nohueckel,loop,ndbconvert,
     1 filename,filenameout,ipsphere,nosort,ispsearch,novolume
      Namelist /Coord/ ICart,IV1,IV2,IV3,R5,R6,leap,isonum,IPRC,
     1 kGC,lGC,IGCtrans,ISW,KE,mirror,IYF,IBF,scaleRad,jumps
      Namelist /FFChoice/ Iopt,ftol,ihessian,iprinth
      Namelist /FFParameters/ fCoulomb,WuR5,WuR6,WuA5,WuA6,WufR5,WufR6,
     1 WufA5,WufA6,ExtWuR55,ExtWuR56,ExtWuR66,ExtWuA5,ExtWuA6,ExtWuDppp,
     1 ExtWuDhpp,ExtWuDhhp,ExtWuDhhh,ExtWufR,ExtWufA,ExtWufD
      Namelist /Hamilton/ IHam,iupac,ihamstore
      Namelist /Isomers/ IPR,IPH,IStop,IChk,ISearch,isomerl,isomerh
      Namelist /Graph/ ISchlegel,ISO1,ISO2,ISO3,nhamcyc,ifs,ndual,PS,
     1 scale,scalePPG,boost

C Input send to output
      if(ilp.eq.0) then   
        WRITE(IOUT,100)
        Do I=1,200
          READ(IN,'(132(A1))',END=11) (DATEN(j),j=1,nzeile)
            endzeile=0
            do j=1,nzeile
              if(DATEN(j).ne.' ') endzeile=j
            enddo
            WRITE(IOUT,60) (DATEN(j),j=1,endzeile)
          enddo
   11     WRITE(IOUT,101)
          REWIND IN
          ilp=1
        else
        WRITE(IOUT,108)
      endif

      do i=1,ffmaxdim
        force(i)=0.d0
        forceP(i)=0.d0
      enddo

C Defining the HO force field using Fowler force constants 
C Distances are taken in Angstroems and angles in rad
C Force constants in N/m for distances and N/m A^2/rad^2 for angles (default values)
      WuR5=1.455d0! in angstroem from solid-state
      WuR6=1.391d0
      WuA5=1.08d2! in deg
      WuA6=1.20d2
      WufR5=390.7d0! from Ceulemans, Fowler
      WufR6=499.7d0
      WufA5=47.88d0*1.45d0**2
      WufA6=80.86d0*1.45d0*1.37d0
      fcoulomb=0.d0
C tolerance parameter (to be used in all force fields)
      ftol=1.d-7


C Defining an extension of the Wu force field (default values)
c     three distances: zero values
      ExtWuR55=1.455d0 ! capprox. that of C20
      ExtWuR56=1.455d0
      ExtWuR66=1.391d0
c     two angles: zero values
      ExtWuA5=1.08d2
      ExtWuA6=1.20d2
c     four dihedrals: zero values (all guessed)
      ExtWuDppp=4.0d1
      ExtWuDhpp=3.5d1
      ExtWuDhhp=2.4d1
      ExtWuDhhh=0.0d0
c     three distances: forces (let's assume they are all the same)
      ExtWufR=390.7d0
c     three angles: forces (let's assume they are all the same)
      ExtWufA=160.4d0*1.45d0*1.37d0
c     four dihedrals: forces (let's assume they are all the same)
      ExtWufD=1.0d2

C Default parameters for external files
      filename= 'Fullerene'
      filenameout= 'Fullerene'

C Integers
      isomerl=1 ! Flag for print start of isomer list
      isomerh=Nisoloop ! Flag for print end of isomer list
      ndbconvert=0 ! Flag for conversion of database
      nosort=0  !  Flag for sorting cartesian coordinates
      ispsearch=1 ! Flag for searching for canonical spiral
      ihamstore=0 ! Flag for storing all Hamiltonian cycles
      nhamcyc=0 ! Flag for reading Hamiltonian cycle for 2D graph
      novolume=0 ! Flag for volume calculation
      IGCtrans=0 ! Initial flag for Goldberg-Coxeter transformed fullerene
      ICart=1   !  Input for fullerene structure
      ichk=0    !  Option for restarting the isomer list
      IER=0     !  Error flag
      ifs=0     !  Option for .dat and .tex files
      iham=0    !  Number of Hamiltonian cycles
      iFS=0     !  Option for producing files for 2D fullerene graphs
      nohueckel=0 ! Option for diagonalizing the Hueckel matrix
      iopt=0    !  No (force field) optimization
      ihessian=0 ! No Hessian matrix produced
      iprinth=0 !  No Hessian matrix printed
      Ipent=0   !  Initial flag for Spriral pentagon input
      IP=0      !  Print option
      IPH=0     !  Print Hamiltonian cycles for each isomer
      IPR=-1    !  Print Isomers
      ISearch=0 !  Search for isomers with RSPIs to closest icosahedral fullerene
      IPRC=0    !  Option for isomer list
      IPSphere=0 ! Option for projecting vertices on minimum covering sphere
      isonum=0  !  Isomer number in database
      ISchlegel=0 !  For graph production, option for type of graph
      ISO1=0    !  Option for fullerene orientation for Schlegel projection
      ISO2=0    !  Option for fullerene orientation for Schlegel projection
      ISO3=0    !  Option for fullerene orientation for Schlegel projection
      istop=0   !  Option for stopping after isomer list
      iupac=0   !  Switch for producing the Iupac nomenclature
                !  iupac=0 just counts Hamiltonian Cycles
      ISW=0     !  Option for Stone-Wales transformation
      IBF=0     !  Option for Brinkmann-Fowler transformation
      IYF=0     !  Option for Yoshido-Fowler transformation
      IV1=2     !  Eigenvector option for fullerene construction
      IV2=3     !  Eigenvector option for fullerene construction
      IV3=4     !  Eigenvector option for fullerene construction
      ixyz=0    !  Option for producing input for ploting program CYLVIEW
      KE=0      !  Endo-Kroto C2 insertion
      kGC=0     !  First Goldberg-Coxeter index
      lGC=0     !  second Goldberg-Coxeter index
      leap=0    !  Initial flag for leapfrog fullerene
      loop=0    !  Option for compound job
      mirror=0  !  Invert coordinates
      NA=60     !  Number of Atoms
      ndual=0   !  Option for plotting dual graph as well

C Reals
      PS=0.d0       ! For graph production, angle input for Schlegel diagram
      scale=2.5d0   ! For graph production, scale Tutte graph
      scalePPG=1.d0 ! For graph production exponential factor in Plestenjak alg.
      boost=1.2d0   ! Extra boost in cone Schlegel projection
      R=1.391d0     ! C-C distance 
      R5=1.455d0    ! Distance in 5-Ring
      R6=R          ! Distance in 6-Ring
      Tol=0.33d0    ! Tolerance
      TolR=0.d0     ! Tolerance for finding ring connections
C     Van der Waals radius of carbon, adjusted approximately to the
C     solid-state results of P.A.Heiney et al., Phys. Rev. Lett. 66, 2911 (1991)
      RVdWC=1.415d0

      scaleRad=4    ! scale size of initial tutte sphere by factor.  The more non-spherical the structure is, the larger this factor should be

      do k=1,10
        jumps(k)=0
      enddo

C Now process namelist input
      READ(IN,'(132(A1))') (DATEN(j),j=1,nzeile)
      endzeile=0
      do j=1,nzeile
        if(DATEN(j).ne.' ') endzeile=j
      enddo
      WRITE(IOUT,60) (DATEN(j),j=1,endzeile)
      WRITE(IOUT,101)
C Read Namelist
      Read(IN,nml=General,Err=99,end=99)
      Read(IN,nml=Coord,Err=99,end=99)
      Read(IN,nml=FFChoice,Err=99,end=99)

C Filenames
      if(filenameout.eq.'Fullerene') filenameout=filename
C Set Parameters for force field
c set forceP (default parameters)[needs to be done after iopt and before opt is read]
  99  if(iopt.eq.1 .or. iopt.eq.2)then
C Wu force field
        force(1)=WuR5
        force(2)=WuR6
        force(3)=WuA5
        force(4)=WuA6
        force(5)=WufR5
        force(6)=WufR6
        force(7)=WufA5
        force(8)=WufA6
        force(9)=fCoulomb
      else if(iopt.eq.3 .or. iopt.eq.4)then
C ExtWu force field
        force(1)=ExtWuR55
        force(2)=ExtWuR56
        force(3)=ExtWuR66
        force(4)=ExtWuA5
        force(5)=ExtWuA6
        force(6)=ExtWuDppp
        force(7)=ExtWuDhpp
        force(8)=ExtWuDhhp
        force(9)=ExtWuDhhh
        force(10)=ExtWufR
        force(11)=ExtWufR
        force(12)=ExtWufR
        force(13)=ExtWufA
        force(14)=ExtWufA
        force(15)=ExtWufD
        force(16)=ExtWufD
        force(17)=ExtWufD
        force(18)=ExtWufD
        force(19)=fCoulomb
      endif

      Read(IN,nml=FFParameters,Err=99,end=99)
      Read(IN,nml=Hamilton,Err=99,end=99)
      Read(IN,nml=Isomers,Err=99,end=99)
       if(isearch.ne.0.and.IPR.eq.-1) IPR=1
      Read(IN,nml=Graph,Err=99,end=99)

c set force (custom parameters)
      if(iopt.eq.1 .or. iopt.eq.2)then
C Wu force field
        force(1)=WuR5
        force(2)=WuR6
        force(3)=WuA5
        force(4)=WuA6
        force(5)=WufR5
        force(6)=WufR6
        force(7)=WufA5
        force(8)=WufA6
        force(9)=fCoulomb
      else if(iopt.eq.3 .or. iopt.eq.4)then
C ExtWu force field
        force(1)=ExtWuR55
        force(2)=ExtWuR56
        force(3)=ExtWuR66
        force(4)=ExtWuA5
        force(5)=ExtWuA6
        force(6)=ExtWuDppp
        force(7)=ExtWuDhpp
        force(8)=ExtWuDhhp
        force(9)=ExtWuDhhh
        force(10)=ExtWufR
        force(11)=ExtWufR
        force(12)=ExtWufR
        force(13)=ExtWufA
        force(14)=ExtWufA
        force(15)=ExtWufD
        force(16)=ExtWufD
        force(17)=ExtWufD
        force(18)=ExtWufD
        force(19)=fCoulomb
      endif

      do i=1,ffmaxdim
        forceP(i)=force(i)
      enddo

C Set IC and ichk parameters
      if(ICart.lt.0) ICart=0
      if(ICart.gt.9) ICart=9
      if(ichk.ne.0) istop=1
      if(ihamstore.ne.0.or.nhamcyc.ne.0) then
        nosort=1
        iupac=1
      endif

C  Check on number of atoms (vertices)
      number_vertices=IABS(NA)
      if(number_vertices.gt.NAtomax) WRITE(IOUT,102) number_vertices
      if(number_vertices.lt.20.or.number_vertices.eq.22) then
        Write(IOUT,103) number_vertices
        IER=1
        return
      endif
      IF (number_vertices/2*2.ne.number_vertices) then
        Write(IOUT,104) number_vertices
        IER=1
        return
      endif

C  Setting minimum distance
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

      if(IPRC.lt.0) IPRC=0
      if(IPRC.eq.1) IPRC=1
      if(IPRC.gt.1) IPRC=0
      if(IPR.le.0) then
        IPR=-1
      endif
      if(IPR.eq.1) then
        if(number_vertices.lt.60) IPR=0
      endif
      if(IPR.ge.2) then
        IPR=0
      endif
      
   60 FORMAT(1X,132A1)
  100 FORMAT(1X,80('-'),/1X,'I N P U T ',/1X,5H0....,
     161H....1.........2.........3.........4.........5.........6......,
     214H...7.........8,/1X,39H123456789012345678901234567890123456789,
     341H01234567890123456789012345678901234567890,/)
  101 FORMAT(1X,132('-'))
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
  108 Format(1X,'Start new job',F12.6)
      RETURN
      END
