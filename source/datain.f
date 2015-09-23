      SUBROUTINE Datain(IN,IOUT,NAtomax,ICart,Iopt,IP,IHam,
     1 nohueckel,KE,IPR,IPRC,ISchlegel,ISO1,ISO2,ISO3,IER,istop,
     1 leap,IGCtrans,iupac,Ipent,IPH,kGC,lGC,IV1,IV2,IV3,IPMC,
     1 irext,iwext,ichk,isonum,loop,mirror,ilp,ISW,IYF,IBF,ifs,
     1 ipsphere,ndual,labelvert,nosort,ispsearch,novolume,ihessian,
     1 isearch,iprinth,ndbconvert,ihamstore,ihamstats,nhamcyc,isomerl,
     1 isomerh,ngaudiene,imcs,itop,
     1 PS,TolX,R5,R6,Rdist,rvdwc,scale,scalePPG,ftol,scaleRad,
     1 rspi,jumps,force,forceP,boost,dualdist,
     1 filename,filenameout,DATEN)
C-----------------------------------------------------------------
C  This is the main routine handling the input
C  It is called from the main program
C-----------------------------------------------------------------
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (nlines=9999)
      integer NA,iopt
      real(8) force(ffmaxdim),forceP(ffmaxdim) ! user chosen FF (and a backup)
      integer endzeile
      integer rspi(12), jumps(10)
      Character*1 DATEN(nzeile)
      Character filename*50,filenameout*50,flagsym*3
      Namelist /General/ NA,IP,TolR,R5,R6,irext,iwext,
     1 nohueckel,loop,ndbconvert,iPFcount,IPMC,imcs,itop,
     1 filename,filenameout,ipsphere,nosort,ispsearch,novolume
      Namelist /Coord/ ICart,IV1,IV2,IV3,R5,R6,IPRC,leap,isonum,
     1 kGC,lGC,IGCtrans,ISW,KE,mirror,IYF,IBF,scaleRad,rspi,jumps,
     1 nanotube,dualdist,ngaudiene
      Namelist /FFChoice/ Iopt,ftol,ihessian,iprinth
      Namelist /FFParameters/ fCoulomb,WuR5,WuR6,WuA5,WuA6,WufR5,WufR6,
     1 WufA5,WufA6,ExtWuR55,ExtWuR56,ExtWuR66,ExtWuA5,ExtWuA6,ExtWuD555,
     1 ExtWuD556,ExtWuD566,ExtWuD666,ExtWufR55,ExtWufR56,ExtWufR66,
     1 ExtWufA5,ExtWufA6,ExtWufD555,ExtWufD556,ExtWufD566,ExtWufD666
      Namelist /Hamilton/ IHam,iupac,ihamstore,ihamstats
      Namelist /Isomers/ IPR,IPH,IStop,Ichk,ISearch,isomerl,isomerh
      Namelist /Graph/ ISchlegel,ISO1,ISO2,ISO3,nhamcyc,ifs,ndual,
     1 labelvert,PS,scale,scalePPG,boost

C Input send to output
      if(ilp.eq.0) then   
        WRITE(IOUT,100)
        Do I=1,nlines
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

C tolerance parameter (to be used in all force fields)
      fcoulomb=0.d0
      ftol=1.d-7

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

C     Defining an extension of the Wu force field (default values)
c     three distances: zero values
      ExtWuR55=1.479d0
      ExtWuR56=1.458d0
      ExtWuR66=1.401d0
c     two angles: zero values
      ExtWuA5=1.08d2
      ExtWuA6=1.20d2
c     four dihedrals: zero values (according to ideal_dihedral)
      ExtWuD555=37.38d0
      ExtWuD556=29.20d0
      ExtWuD566=23.49d0
      ExtWuD666=0.0d0
c     three distances: forces (let's assume they are all the same)
      ExtWufR55=260.d0
      ExtWufR56=390.d0
      ExtWufR66=450.d0
c     three angles: forces (let's assume they are all the same)
      ExtWufA5=100.d0
      ExtWufA6=100.d0
c     four dihedrals: forces (let's assume they are all the same)
      ExtWufD555=35.d0
      ExtWufD556=65.d0
      ExtWufD566=85.d0
      ExtWufD666=270.d0

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
      ihamstats=0 ! Flag for Hamiltonian cycle statistics
      nhamcyc=0 ! Flag for reading Hamiltonian cycle for 2D graph
      novolume=0 ! Flag for volume calculation
      nohueckel=0 ! Option for diagonalizing the Hueckel matrix
      nanotube=0!  Flag for creating nanotubes
      IGCtrans=0 ! Initial flag for Goldberg-Coxeter transformed fullerene
      ICart=1   !  Input for fullerene structure
      ichk=0    !  Option for restarting the isomer list
      IER=0     !  Error flag
      ifs=0     !  Option for .dat and .tex files
      iham=0    !  Number of Hamiltonian cycles
      iFS=0     !  Option for producing files for 2D fullerene graphs
      imcs=0    !  Option for cartesian input and only producing minimum covering sphere
      itop=0    !  Option for stopping after topological analysis (and not creating 3D coordinates after GC transformation)
      iPMC=0    !  Option for perfect match count
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
      irext=0   !  Option for reading coordinates/connectivities from external file
      iwext=0   !  Option for writing coordinates/connectivities from to file
      KE=0      !  Endo-Kroto C2 insertion
      kGC=0     !  First Goldberg-Coxeter index
      lGC=0     !  second Goldberg-Coxeter index
      leap=0    !  Initial flag for leapfrog fullerene
      loop=0    !  Option for compound job
      mirror=0  !  Invert coordinates
      NA=60     !  Number of Atoms
      ndual=0   !  Option for plotting dual graph as well
      ngaudiene=0   !  Option for plotting gaudiene structure
      labelvert=0   !  Option labeling vertices in the 2D graph

C Reals
      PS=0.d0       ! For graph production, angle input for Schlegel diagram
      scale=2.5d0   ! For graph production, scale Tutte graph
      scalePPG=1.d0 ! For graph production exponential factor in Plestenjak alg.
      boost=1.2d0   ! Extra boost in cone Schlegel projection
      R=1.391d0     ! C-C distance 
      R5=1.455d0    ! Distance in 5-Ring
      R6=R          ! Distance in 6-Ring
      dualdist=R    ! Scaling dual distances
      Tol=0.33d0    ! Tolerance
      TolR=0.d0     ! Tolerance for finding ring connections
C     Van der Waals radius of carbon, adjusted approximately to the
C     solid-state results of P.A.Heiney et al., Phys. Rev. Lett. 66, 2911 (1991)
      RVdWC=1.415d0

      scaleRad=4    ! scale size of initial tutte sphere by factor.  The more non-spherical the structure is, the larger this factor should be

c init of rspi (always 12)
      do k=1,12
        rspi(k)=0
      enddo
c init of jumps (should be more than 10 (should ... ))
      do k=1,10
        jumps(k)=0
      enddo

C Now process namelist input
      READ(IN,'(132(A1))',Err=98,end=98) (DATEN(j),j=1,nzeile)
   98 endzeile=0
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
      else if(iopt.eq.3.or.iopt.eq.4.or.iopt.eq.5.or.iopt.eq.6)then
C ExtWu force field
        force(1)=ExtWuR55
        force(2)=ExtWuR56
        force(3)=ExtWuR66
        force(4)=ExtWuA5
        force(5)=ExtWuA6
        force(6)=ExtWuD555
        force(7)=ExtWuD556
        force(8)=ExtWuD566
        force(9)=ExtWuD666
        force(10)=ExtWufR55
        force(11)=ExtWufR56
        force(12)=ExtWufR66
        force(13)=ExtWufA5
        force(14)=ExtWufA6
        force(15)=ExtWufD555
        force(16)=ExtWufD556
        force(17)=ExtWufD566
        force(18)=ExtWufD666
        force(19)=fCoulomb
      endif

      Read(IN,nml=FFParameters,Err=99,end=99)
      Read(IN,nml=Hamilton,Err=99,end=99)
      Read(IN,nml=Isomers,Err=99,end=99)
       if(isearch.ne.0.and.IPR.eq.-1) IPR=1
      Read(IN,nml=Graph,Err=99,end=99)

c set ischlegel if ifs is non-zero
      if(ifs.ne.0.and.ischlegel.eq.0) ischlegel=1
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
      else if(iopt.eq.3.or.iopt.eq.4.or.iopt.eq.5.or.iopt.eq.6)then
C ExtWu force field
        force(1)=ExtWuR55
        force(2)=ExtWuR56
        force(3)=ExtWuR66
        force(4)=ExtWuA5
        force(5)=ExtWuA6
        force(6)=ExtWuD555
        force(7)=ExtWuD556
        force(8)=ExtWuD566
        force(9)=ExtWuD666
        force(10)=ExtWufR55
        force(11)=ExtWufR56
        force(12)=ExtWufR66
        force(13)=ExtWufA5
        force(14)=ExtWufA6
        force(15)=ExtWufD555
        force(16)=ExtWufD556
        force(17)=ExtWufD566
        force(18)=ExtWufD666
        force(19)=fCoulomb
      endif

      do i=1,ffmaxdim
        forceP(i)=force(i)
      enddo

C Check imcs parameter
      if(imcs.eq.1.and.ICart.ne.1) then
        write(*,*)"icms and icart parameters don't fit together"
        call exit(1)
      endif
C Set IC and ichk parameters
      if(ICart.lt.0 .or. icart.gt.10) then
        write(*,*)"Invalic value for icart given.  Exiting ..."
        call exit(1)
      endif
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

C Create rspi if nanotube.ne.0
      if(nanotube.ne.0)then
       flagsym='   '
       Call RNanotube(nanotube,number_vertices,rspi,ierrnano,flagsym)
       if(ierrnano.ne.0) then
        Write(IOUT,110) number_vertices
        Stop
       endif
       Write(IOUT,109) flagsym,(rspi,I=1,12)
      endif
C  Setting minimum distance
      if(imcs.eq.0) then
      if(R6.ne.R.and.R6.gt.1.d0) then
        Rdist=R6
        WRITE(Iout,106) Rdist
      else
        Rdist=R
        WRITE(Iout,107) Rdist
      endif
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
  109 Format(1X,'Create ring spiral pentagon indices for smallest ',
     1 'nanotube of ',A3,' symmetry:',/,' rspi= ',11(I6,','),I6,/)
  110 Format(1X,'Number of vertices ',I6,
     1 ' not consistent with a nanotube')
      RETURN
      END

      Subroutine RNanotube(nano,NV,rspi,ierrnano,csym)
      IMPLICIT INTEGER (A-Z)
      Character*3 csym
      integer rspi(12)
C     number of faces
      Nf=NV/2+2
      ierrnano=0
C     Smallest nanotube
      if(nano.eq.1) then
       ncount=NV/10
       if(ncount*10.ne.NV) then
        ierrnano=1
        return
       endif
       neven=ncount/2
       if(neven*2.eq.ncount) then
        csym='D5h'
       else
        csym='D5d'
       endif
       rspi(1)=1
       rspi(2)=2
       rspi(3)=3
       rspi(4)=4
       rspi(5)=5
       rspi(6)=6
       rspi(7)=Nf-5
       rspi(8)=Nf-4
       rspi(9)=Nf-3
       rspi(10)=Nf-2
       rspi(11)=Nf-1
       rspi(12)=Nf
      endif
      if(nano.eq.2) then
       ncount=NV/12
       if(ncount*12.ne.NV) then
        ierrnano=1
        return
       endif
       neven=ncount/2
       if(neven*2.eq.ncount) then
        csym='D6h'
       else
        csym='D6d'
       endif
       rspi(1)=2
       rspi(2)=3
       rspi(3)=4
       rspi(4)=5
       rspi(5)=6
       rspi(6)=7
       rspi(7)=Nf-6
       rspi(8)=Nf-5
       rspi(9)=Nf-4
       rspi(10)=Nf-3
       rspi(11)=Nf-2
       rspi(12)=Nf-1
      endif
      RETURN
      END
      Subroutine ReadFromFile(nchoice,iextfile,iout,iatom,IC3,
     1 extfilename,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C-----------------------------------------------------------------
C  Routine to read cartesian coordinates from external file
C  It is called from the main program
C  Formats: 
C   nchoice=1   .xyz file
C   nchoice=2   .cc1 file
C   nchoice=3   .mol2 file (TRYPOS format)
C  iextfile: unit number for external file
C  iout:     unit number for output
C  iatom:    Field for atom number for each atom (6 for carbon)
C  IC3:      Connectivity field for cubic graph 
C             (short form of adjacency matrix)
C  extfilename: external file name
C  Dist(3,i): Field of (x,y,z) coordinates for each atom i
C-----------------------------------------------------------------
      DIMENSION Dist(3,Nmax),Iatom(Nmax),IC3(Nmax,3)
      CHARACTER*50 extfilename
      Character*1 TEXTINPUT(nzeile)
      CHARACTER*2 element,bondtype
      CHARACTER*5 atomname,atomtype,substname
      CHARACTER*132 Line
      Integer endzeile
      Logical lexist

      inquire(file=extfilename,exist=lexist)
       if(lexist.neqv..True.) then
         Write(Iout,1001) extfilename
         stop
       endif
      Open(unit=iextfile,file=extfilename,form='formatted')
      WRITE(Iout,1000) extfilename,nchoice

C .xyz files
      if(nchoice.eq.1) then
       Read(iextfile,*,end=99) number_vertices
       Write(Iout,1007) number_vertices
       Read(iextfile,1002,end=99) (TEXTINPUT(I),I=1,nzeile)
       endzeile=0
       do j=1,nzeile
         if(TEXTINPUT(j).ne.' ') endzeile=j
       enddo 
       WRITE(Iout,1003) (TEXTINPUT(I),I=1,endzeile)
       Do J=1,number_vertices
        Read(iextfile,*,end=1,err=1) element,(Dist(I,J),I=1,3)
        Write(Iout,1009) element,(Dist(I,J),I=1,3)
        Iatom(j)=6
       enddo
      endif

C .cc1 files
      if(nchoice.eq.2) then
       Read(iextfile,*,end=99) number_vertices
       Write(Iout,1007) number_vertices
       Do J=1,number_vertices
        Do I=1,3
         IC3(J,I)=0
        enddo
        Read(iextfile,'(A132)',err=1,end=1) Line
        Read(Line,*,end=10,err=10) 
     1   element,JJ,(Dist(I,J),I=1,3),ncc1flag,(IC3(J,I),I=1,3)
   10    Write(Iout,1008) element,JJ,(Dist(I,J),I=1,3),ncc1flag,
     1    (IC3(J,I),I=1,3)
        Iatom(j)=6
       enddo
      endif

C .mol2 files: tripos mol2 standard format
      if(nchoice.eq.3) then
C    Read comment section
       do I =1,1000
        Read(iextfile,1002,end=99) (TEXTINPUT(J),J=1,nzeile)
        endzeile=0
        do j=1,nzeile
         if(TEXTINPUT(j).ne.' ') endzeile=j
        enddo 
        WRITE(Iout,1005) (TEXTINPUT(J),J=1,endzeile)
        do J=1,endzeile
         if(TEXTINPUT(J).eq.'@'.and.TEXTINPUT(J+1).eq.'<') then
          if(TEXTINPUT(J+9).eq.'M'.or.TEXTINPUT(J+9).eq.'m') go to 21
         endif
        enddo
       enddo
       go to 99
C    Comment section over, now read molecule section after @<TRIPOS>MOLECULE
  21   Read(iextfile,1002,end=99) (TEXTINPUT(I),I=1,nzeile)
       WRITE(Iout,1006) (TEXTINPUT(I),I=1,endzeile)
        endzeile=0
        do j=1,nzeile
         if(TEXTINPUT(j).ne.' ') endzeile=j
        enddo 
       Read(iextfile,*,end=99,err=99) number_vertices,number_edges,
     1  isubstruct,ifeatures,isets
       WRITE(Iout,1010) number_vertices,number_edges,
     1  isubstruct,ifeatures,isets
       Do I=1,1000
        Read(iextfile,1002,end=99) (TEXTINPUT(J),J=1,nzeile)
         endzeile=0
         do j=1,nzeile
          if(TEXTINPUT(j).ne.' ') endzeile=j
         enddo 
        WRITE(Iout,1006) (TEXTINPUT(J),J=1,endzeile)
         do J=1,endzeile
          if(TEXTINPUT(J).eq.'@'.and.TEXTINPUT(J+1).eq.'<') then
           if(TEXTINPUT(J+9).eq.'A'.or.TEXTINPUT(J+9).eq.'a') go to 22
          endif
         enddo
        enddo
  22   Do I=1,number_vertices
        Read(iextfile,*,end=99) idatom,atomname,(Dist(J,I),J=1,3),
     1   atomtype,idsubst,substname,charge 
        WRITE(Iout,1011) idatom,atomname,(Dist(J,I),J=1,3),
     1   atomtype,idsubst,substname,charge
       enddo
       Do I=1,1000
        Read(iextfile,1002,end=99) (TEXTINPUT(J),J=1,nzeile)
         endzeile=0
         do j=1,nzeile
          if(TEXTINPUT(j).ne.' ') endzeile=j
         enddo 
        WRITE(Iout,1006) (TEXTINPUT(J),J=1,endzeile)
         do J=1,endzeile
          if(TEXTINPUT(J).eq.'@'.and.TEXTINPUT(J+1).eq.'<') then
           if(TEXTINPUT(J+9).eq.'B'.or.TEXTINPUT(J+9).eq.'b') go to 23
          endif
         enddo
        enddo
C    @<TRIPOS>ATOM section done, now read @<TRIPOS>BOND section
       Do J=1,number_vertices
        Do I=1,3
         IC3(J,I)=0
        enddo
       enddo
  23   Do I=1,number_edges
        Read(iextfile,*,end=99) idbond,idoriginatom,idtargetatom,
     1   bondtype
        WRITE(Iout,1012) idbond,idoriginatom,idtargetatom,bondtype
        Call FillIC3(idoriginatom,idtargetatom,IC3)
       enddo
       Do I=1,1000
        Read(iextfile,1002,end=1) (TEXTINPUT(J),J=1,nzeile)
         endzeile=0
         do j=1,nzeile
          if(TEXTINPUT(j).ne.' ') endzeile=j
         enddo 
        WRITE(Iout,1006) (TEXTINPUT(J),J=1,endzeile)
        enddo
      endif

   1     close(unit=7)
         return

  99     WRITE(Iout,1004)
         close(unit=7)
         stop
 1000 FORMAT(/1X,'Read coordinates from external file: ',A60,
     1 /1X,'Choice: ',I1,', Data content:',/)
 1001 Format(1X,'Filename ',A50,' in database not found ==> ABORT')
 1002 FORMAT(132A1)
 1003 FORMAT(1X,132A1)
 1004 FORMAT(1X,'File cannot be read')
 1005 FORMAT(1X,132A1)
 1006 FORMAT(1X,132A1)
 1007 FORMAT(1X,'File content:'/1X,'Number of vertices: ',I10)
 1008 FORMAT(1X,A2,I5,3F12.5,4I5)
 1009 FORMAT(1X,A2,3F12.5)
 1010 Format(5I6)
 1011 Format(I5,1X,A5,1X,3F10.5,1X,A5,1X,I2,1X,A5,1X,F10.5)
 1012 Format(I6,1X,I6,1X,I6,3X,A2)
      END

      Subroutine FillIC3(IA,IB,IC3)
      use config
      DIMENSION IC3(Nmax,3)
      Do I=1,3
       if(IC3(IA,I).eq.0) then
        IC3(IA,I)=IB
        go to 1
       endif
      enddo
  1   Do I=1,3
       if(IC3(IB,I).eq.0) then
        IC3(IB,I)=IA
        return
       endif
      enddo
      return
      END

