!                                                                           !
!                  P R O G R A M     F U L L E R E N E                      !
!                                                                           !
!---------------------------------------------------------------------------!
!  This is an open-source code, see copyright message in the User Manual.   !
!---------------------------------------------------------------------------!
!                                                                           !
! A PROGRAM FOR STRUCTURE GENERATION AND TOPOLOGICAL ANALYSIS OF FULLERENES !
!    The program creates cartesian coordinates for fullerenes isomers       !
!        and performs a topological/graph theoretical analysis.             !
!      The results can be used for plotting 2D/3D fullerene graphs          !
!    (e.g. Schlegel diagrams) and structures, and as a starting point       !
!            for further quantum theoretical treatment.                     !
!        Version 4.5 incorporates C++ routines linked to the                !
!       original Fortran program using much improved algorithms.            !
!---------------------------------------------------------------------------!

!---------------------------------------------------------------------------!
!  This main routine calls all major subroutines                            !
!---------------------------------------------------------------------------!

      PROGRAM Fullerene
      use iso_c_binding
      use config

      IMPLICIT REAL*8 (A-H,O-Z)
C    Set the dimensions for the distance matrix
      real(8) force(ffmaxdim), forceP(ffmaxdim)
      Real(4) TimeX
      DIMENSION CRing5(3,Mmax),CRing6(3,Mmax),cmcs(3),CR(3,Mmax)
      DIMENSION DistMat(NmaxL),Dist(3,Nmax),DistCM(3),Dist2D(2,Nmax)
      DIMENSION DistStore(3,Nmax)
      DIMENSION A(Nmax,Nmax),evec(Nmax),df(Nmax)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6),Iring(Mmax)
      DIMENSION distP(Nmax),IDA(Nmax,Nmax)
      DIMENSION IATOM(Nmax),IC3(Nmax,3),Nring(Mmax)
      integer IVR3(nmax+4,3) ! in ring.f up to four values to many are read
      DIMENSION NringA(Emax),NringB(Emax)
      DIMENSION NringC(Emax),NringD(Emax)
      DIMENSION NringE(Emax),NringF(Emax)
      DIMENSION IDual(Mmax,Mmax),nSW(4,66),nFM(4,66),nYF(6,66),nBF(5,66)
      DIMENSION NEK(3,Nmax)
      DIMENSION Symbol(Mmax)
      CHARACTER CDAT*8,CTIM*10,Zone*5
      CHARACTER*1 Symbol
      CHARACTER*2 El(99)
      CHARACTER*7 Namecc1,Namexyz,Namemol
      CHARACTER*4 Endcc1,Endxyz,Endmol
      CHARACTER*15 routine
      CHARACTER*50 filename,filenameout,extname
      CHARACTER*50 xyzname,cc1name,molname
      Character*1 TEXTINPUT(nzeile)
      CHARACTER*3 GROUP
      Integer Values(8)
      integer istop
      integer mdist(nmax,nmax)
      integer rspi(12),jumps(10)
      integer leapspiral,SWspiral, scaleRad

C Set parameters
C Element Names
      DATA El/' H','HE','LI','BE',' B',' C',' N',' O',' F','NE','NA',
     1 'MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V','CR',
     2 'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',    
     3 'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD',   
     4 'IN','SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',  
     5 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF', 
     6 'TA',' W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',
     7 'AT','RN','FR','RA','AC','TH','PA',' U','NP','PU','AM','CM',   
     8 'BK','CF','ES'/                                               

C External file names
      Namecc1='-3D.cc1'
      Namexyz='-3D.xyz'
      Namemol='-3D.mol'
      Endcc1='.cc1'
      Endxyz='.xyz'
      Endmol='.mol'
      nxyz=0
      ncc1=0
      nmol=0

C Input / Output
      IN=5
      Iout=6
      Iext=7

C Set parameters to zero
      nloop=0
      ilp=0
      iprev=0
      ihalma=0
      VolSphere=0.d0
      ASphere=0.d0
      Do I=1,Nmax
       IAtom(I)=6
        Do J=1,Nmax
          IDA(I,J)=0
        enddo
      enddo
      Do I=1,Nmax
       Dist(1,I)=0.d0
       Dist(2,I)=0.d0
       Dist(3,I)=0.d0
       Dist2D(1,I)=0.d0
       Dist2D(2,I)=0.d0
      enddo

C------------------TIME-AND-DATE-----------------------------------
C Get time and date
      CALL date_and_time(CDAT,CTIM,zone,values)
      TIMEX=0.d0
      CALL Timer(TIMEX)
      WRITE(Iout,1000) Values(3),Values(2),Values(1),Values(5),
     1  Values(6),Values(7),Nmax

C------------------------------------------------------------------
C  S T A R T   O F   I N P U T   S E C T I O N
C------------------------------------------------------------------

C------------------DATAIN------------------------------------------
C  INPUT and setting parameters for running the subroutines
 9    routine='DATAIN         '
      Write(Iout,1008) routine
      Group='   '
      isort=0
      leapspiral=0
      SWspiral=0
C  Next two flags tells us if input has information on Cartesian Coordinates
C   or Adjacency Matrix (through connectivity IC3)
      ncartflag=0
      nadjacencyflag=0
C  Call Datain (main input routine)
      CALL Datain(IN,Iout,Nmax,Icart,Iopt,iprintf,IHam,
     1 nohueckel,KE,IPR,IPRC,ISchlegel,IS1,IS2,IS3,IER,istop,
     1 leap,leapGC,iupac,Ipent,iprintham,IGC1,IGC2,IV1,IV2,IV3,iPMC,
     1 irext,iwext,ichk,isonum,loop,mirror,ilp,ISW,IYF,IBF,ifs,
     1 ipsphere,ndual,labelvert,nosort,ispsearch,novolume,ihessian,
     1 isearch,iprinthessian,ndbconvert,ihamstore,ihamstats,nhamcyc,
     1 isomerl,isomerh,ngaudiene,imcs,itop,
     1 ParamS,TolX,R5,R6,Rdist,rvdwc,scales,scalePPG,
     1 ftolP,scaleRad,rspi,jumps,force,forceP,boost,
     1 dualdist,filename,filenameout,TEXTINPUT)

C  Simple checks
       Rmin5=R5
       Rmin6=R6
C  Stop if error in input
      If(IER.ne.0) go to 99
C  Stop if isomer closest to icosahedral is searched for
      if(isearch.ne.0) then
        istop=1
        go to 98
      endif
C  Only do isomer statistics
      if(istop.ne.0) go to 98

C------------------COMPRESS----------------------------------------
C  Convert printed database into a more compressed file
C  This is a routine used by the programmers to compresse the
C  database files to a reasonable format
      if(ndbconvert.ne.0) then
      routine='COMPRESS       '
      Write(Iout,1008) routine
       Call CompressDatabase(Iout,filename)
       go to 99
      endif

C-------------Coordinates and Connectivities-----------------------
C This controls how fullerene structure is read in

C Input Cartesian coordinates/connectivities from external file
C if irext.ne.0: overwrites completely the Icart option
      if(irext.ne.0) then
       if(irext.eq.1) then

C Read from .xyz file
        xyzname=trim(filename)//".xyz"
        Call ReadFromFile(1,Iext,iout,iatom,IC3,xyzname,Dist)
        ncartflag=1
        xyzname=trim(filename)//'-3D.new.xyz'
        go to 40
       endif

C Read from .cc1 file
       if(irext.eq.2) then
        cc1name=trim(filename)//".cc1"
        Call ReadFromFile(2,Iext,iout,iatom,IC3,cc1name,Dist)
        ncartflag=1
C    This routine complements missing entries in IC3
        Call CheckIC3(IERROR,IC3)
        if(IERROR.eq.1) then
         nadjacencyflag=0
         Write(iout,1012)
        else
         nadjacencyflag=1
        endif
        cc1name=trim(filename)//'-3D.new.cc1'
        go to 40
       endif

C Read from .mol2 file
       if(irext.eq.3) then
        molname=trim(filename)//".mol2"
        Call ReadFromFile(3,Iext,iout,iatom,IC3,molname,Dist)
        ncartflag=1
C    This routine complements missing entries in IC3
        Call CheckIC3(IERROR,IC3)
        if(IERROR.eq.1) then
         nadjacencyflag=0
         Write(iout,1012)
        else
         nadjacencyflag=1
        endif
        molname=trim(filename)//'-3D.new.mol2'
        go to 40
       endif
      endif

C Options for direct Input 
      go to (10,20,30,30,30,30,30,30,30,30,98) Icart+1

C  Cartesian coordinates produced for Ih C20 or C60 using basic geometry
   10 routine='COORDC20/60    '
      Write(Iout,1008) routine
      CALL CoordC20C60(Iout,R5,R6,Dist)
      ncartflag=1
      Do I=1,number_vertices
        IAtom(I)=6
      enddo
      Go to 40

C Read cartesian coordinates directly
   20 Do J=1,number_vertices
        Read(IN,*,end=21) IAtom(J),(Dist(I,J),I=1,3)
       enddo
       ncartflag=1
       if(imcs.ne.0) then
        WRITE(Iout,1002)
        CALL MoveCM(Iout,Iprint,IAtom,mirror,isort,
     1   nosort,SP,Dist,DistCM,El)
        CALL Diameter(Iout,Dist,distp)
        Rmin=1.d15
        Do I=1,number_vertices
        Do J=I+1,number_vertices
         X=Dist(1,I)-Dist(1,J)
         Y=Dist(2,I)-Dist(2,J)
         Z=Dist(3,I)-Dist(3,J)
         R2=X*X+Y*Y+Z*Z
         R=dsqrt(R2)
         if(R.lt.Rmin) Rmin=R
        enddo
        enddo
        go to 97
       endif
       Go to 40
   21  WRITE(Iout,1013)
       Go to 99

   30 Ipent=1
      routine='COORDBUILD     '
      Write(Iout,1008) routine
      CALL CoordBuild(IN,Iout,IDA,IDual,
     1 Icart,IV1,IV2,IV3,IGC1,IGC2,isonum,IPRC,nohueckel,
     1 iprev,ihalma,A,evec,df,Dist,Dist2D,distp,Rdist,scaleRad,
     1 rspi,jumps,GROUP,filename)
      Do I=1,number_vertices
        IAtom(I)=6
      enddo
      ncartflag=1

   40 WRITE(Iout,1001) number_vertices,TolX*100.d0

C------------------------------------------------------------------
C  E N D   O F   I N P U T   S E C T I O N
C------------------------------------------------------------------

C------------------------------------------------------------------
C  S T A R T   O F   T O P O L O G Y   S E C T I O N
C------------------------------------------------------------------
C------------------ISOMERS-----------------------------------------
C Some general infos on isomers and spiral routine
C of Fowler and Manolopoulos. Set parameter IPR for independent
C pentagon rule as full list beyond C60 is computer time 
C intensive
  98  routine='ISOMERS        '
      Write(Iout,1008) routine
      CALL Isomers(IPR,isearch,In,Iout,iprintham,ihamstats,
     1 isomerl,isomerh,ichk,IDA,A,filename)
      if(istop.ne.0) go to 99

C------------------MOVECM------------------------------------------
C Move carbon cage to Atomic Center
  999 routine='MOVECM_1       '
      Write(Iout,1008) routine
      Iprint=iprintf
      CALL MoveCM(Iout,Iprint,IAtom,mirror,isort,
     1 nosort,SP,Dist,DistCM,El)
      mirror=0

C------------------DIAMETER----------------------------------------
C Calculate largest and smallest atom-to-atom diameters
C Also get moment of inertia
       routine='DIAMETER       '
       Write(Iout,1008) routine
       CALL Diameter(Iout,Dist,distp)

C------------------DISTMATRIX--------------------------------------
C Calculate the distance Matrix and print out distance Matrix
      routine='DISTANCEMATRIX '
      Write(Iout,1008) routine
      CALL DistMatrix(Iout,iprintf,0,Iopt,
     1 Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)

C------------------CONNECT-----------------------------------------
C Establish Connectivities
      routine='CONNECT        '
      Write(Iout,1008) routine
      CALL Connect(MCon2,Ipent,Iout,IC3,IDA,TolX,DistMat,Rmin)

C------------------REORDER-----------------------------------------
C Reorder atoms such that distances in internal coordinates are bonds
      if(nosort.eq.0) then
       routine='REORDER        '
       Write(Iout,1008) routine
       CALL Permute(Iout,nperm,IC3,IDA,Dist)
        if(nperm.ne.0) then
         CALL MoveCM(Iout,Iprint,IAtom,mirror,isort,
     1    nosort,SP,Dist,DistCM,El)
         CALL DistMatrix(Iout,0,1,Iopt,
     1    Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)
        endif
      endif

C------------------HUECKEL-----------------------------------------
C Hueckel matrix and eigenvalues
      if(ipent.eq.0) then
        routine='HUECKEL        '
        Write(Iout,1008) routine
        CALL Hueckel(Iout,IC3,nohueckel,IDA,A,evec,df)
      endif

C------------------GOLDBERG-COXETER-------------------------------
C Produce the nth leapfrog of the fullerene
      if(leap.gt.0.or.leapGC.gt.0) then
        routine='GOLDBERGCOXETER'
        Write(Iout,1008) routine
        CALL GoldbergCoxeter(Iout,leap,leapGC,IGC1,IGC2,itop,
     1   nohueckel,LeapErr,IDA,A,evec,df,Dist,Dist2D,distp,Rdist,
     1   scaleRad)
        leap=0
        leapGC=0
        ipent=1
        leapspiral=1
        if(number_vertices.gt.100) IHam=0
C   Write out IC3 on external file
        if(itop.eq.2) then
         Call CubeConnect(Iout,IDA,IC3)
         iext=1
         extname='ic3file'
         Open(unit=Iext,file=extname,form='formatted')
         Write(iext,*) number_vertices
         Do I=1,number_vertices
          Write(iext,*) (IC3(I,J),J=1,3)
         enddo
         Close(unit=Iext)
         Go to 888
        endif
        if(LeapErr.eq.0) go to 999 ! moveCM
      endif

C     Check if only topological analysis is needed 
      if(itop.ne.0) go to 888

C------------------HAMILTON---------------------------------------
C Generate IUPAC name and locate Hamiltonian cycles. 
C Routine written by D. Babic. Note routine
C is called only if IPR>0 as computer time is extensive beyond
C C100 (PN-hard problem). Last routine uses the adjaceny matrix
C to calculate the number of all distinct paths between 
C adjacent vertices
      routine='HAMILTON       '
      Write(Iout,1008) routine
      maxiter=maxit
      if(IHam.gt.1.and.IHam.le.9) then
        maxiter=10**IHam
      endif
      if(IHam.ne.0 .and. 
     1    ke.eq.0 .and. isw.eq.0 .and. iyf.eq.0 .and. ibf.eq.0) then
        if(iupac.ne.0.or.iprintf.ne.0) then
          CALL Hamilton(Iout,iprintf,ihamstore,maxiter,IC3,filename)
          Close(unit=8)
        else
          CALL HamiltonCyc(maxiter,Iout,nbatch,IC3,Nhamilton)
          WRITE(Iout,1010) Nhamilton
          if(nbatch.ne.0) WRITE(Iout,1011)
        endif
      endif
      CALL PathStatistic(Iout,iprintf,IDA,A,evec,df)

C------------------RING-------------------------------------------
C Establish all closed ring systems
      routine='RING           '
      Write(Iout,1008) routine
      CALL Ring(Medges,MCon2,Iout,N5Ring,N6Ring,
     1 IC3,IVR3,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,DistMat)

C------------------RINGC------------------------------------------
C Analyze ring connections
      routine='RINGC          '
      Write(Iout,1008) routine
      CALL RingC(Medges,Iout,iprintf,
     1 N5Ring,N6Ring,NRing,Iring5,Iring6,Iring56,
     1 n565,numberSW,numberFM,numberYF,numberBF,
     1 N5MEM,N6MEM,NringA,NringB,NringC,NringD,NringE,NringF,
     1 IC3,IVR3,nEK,nSW,nFM,nYF,nBF,SmallRingDist,DIST,CRing5,CRing6)
C     Print edge coordinates (barycenter)
      if(iprintf.ne.0) Call EdgeCoord(Iout,DIST,IC3)
      if(iprintf.ne.0) Call RingCoord(Iout,dualdist,R6,
     1 SmallRingDist,DIST,N5Ring,N6Ring,N5MEM,N6MEM)

C------------------STONE-WALES------------------------------------
C Perform Stone-Wales transformation
      if(ISW.ne.0) then
        routine='STONE-WALES    '
        Write(Iout,1008) routine
        CALL StoneWalesTrans(IN,Iout,numbersw,nSW,
     1   nohueckel,IDA,N6MEM,IC3,A,evec,df,Dist,Dist2D,distp,
     1   Rdist,scalerad)
        ISW=0
        ipent=1
        SWspiral=1
        go to 999 ! moveCM
      endif

C------------------ENDO-KROTO-------------------------------------
C Perform Endo-Kroto 2-vertex insertion
      if(KE.ne.0) then
        routine='ENDO-KROTO     '
        Write(Iout,1008) routine
        CALL EndoKrotoTrans(IN,Iout,n565,NEK,nohueckel,
     1   IDA,N5MEM,N6MEM,A,evec,df,Dist,Dist2D,distp,Rdist,scalerad)
        KE=0
        ipent=1
        SWspiral=1
        go to 999 ! moveCM
      endif

C------------------YOSHIDA-FOWLER---------------------------------
C Perform Yoshida-Fowler 4-or 6-vertex insertion
      if(IYF.ne.0) then
        routine='YOSHIDAFOWLER  '
        Write(Iout,1008) routine
        if(IYF.le.2) then
         CALL YoshidaFowler4(IN,Iout,JERR,numberFM,IYF,nFM,
     1    nohueckel,IDA,N5MEM,N6MEM,A,evec,df,Dist,Dist2D,distp,Rdist,
     1    scalerad)
        else
          IYF=IYF-2
          CALL YoshidaFowler6(IN,Iout,JERR,numberYF,IYF,nYF,
     1     nohueckel,IDA,N6MEM,IC3,A,evec,df,Dist,Dist2D,distp,Rdist,
     1     scalerad)
        endif
        IYF=0
        ipent=1
        SWspiral=1
        if(JERR.eq.0) go to 999 ! moveCM
      endif

C------------------BRINKMANN-FOWLER-------------------------------
C Perform Brinkmann-Fowler 6-vertex 6-55-55 insertion
      if(IBF.ne.0) then
        routine='BRINKMANNFOWLER'
        Write(Iout,1008) routine
        CALL BrinkmannFowler(IN,Iout,JERR,numberBF,IBF,
     1   nBF,nohueckel,IDA,N5MEM,N6MEM,IC3,
     1   A,evec,df,Dist,Dist2D,distp,Rdist,scalerad)
        IBF=0
        ipent=1
        SWspiral=1
        if(JERR.eq.0) go to 999 ! moveCM
      endif

C------------------SPIRALSEARCH-----------------------------------
C Now produce clockwise spiral ring pentagon count a la Fowler and Manolopoulos
      if((ipent.eq.0 .or. leapspiral.ne.0.or.SWspiral.ne.0.
     1   or.Icart.eq.6.or.Icart.eq.7.or.ihalma.eq.1).or.
     1   ispsearch.ne.0) then
        routine='SPIRALSEARCH   '
        ispcount=0
        if(ispsearch.gt.1) ispcount=1
        Write(Iout,1008) routine
        CALL SpiralSearch(Iout,Iring5,Iring6,Iring56,NringA,NringB,
     1   NringC,NringD,NringE,NringF,rspi,GROUP,ispcount)
      endif

C--------------TOPOLOGICAL INDICATORS-----------------------------
      routine='TOPOLOINDICATOR'
 888  Write(Iout,1008) routine
C Topological Indicators
      CALL TopIndicators(Iout,iPMC,IDA,mdist)
C Check if vertex number allows for icosahedral fullerenes
      Call IcoFullDetect(Iout)
C Determine if fullerene is chiral
      CALL Chiral(Iout,GROUP)
      if(itop.ne.0) go to 99
C------------------------------------------------------------------
C  E N D   O F   T O P O L O G Y   S E C T I O N
C------------------------------------------------------------------

C------------------------------------------------------------------
C  S T A R T   O F   3D   S T R U C T U R E   S E C T I O N
C------------------------------------------------------------------

C Cartesian coordinates produced ring from spiral pentagon list
C or from adjacency matrix. Uses the Fowler-Manolopoulos algorithm 
C with P-type eigenvectors or the Tutte algorithm to construct 
C the 3D fullerene
  
C  More to come here
      routine='COORDBUILD     '
      Write(Iout,1008) routine

C------------------OPTFF------------------------------------------
C Optimize Geometry through force field method
      If(Iopt.ne.0) then
c       Store distances
        Do I=1,3
          Do J=1,number_vertices 
            DistStore(I,J)=Dist(I,J)
          enddo
        enddo
        routine='OPTFORCEFIELD  '
        ftol=ftolP
        Write(Iout,1008) routine
        call flush(iout)
        if(Iopt.eq.1 .or. Iopt.eq.2) then ! vanilla Wu or Wu + Coulomb
          if(Iopt.eq.2) then ! Wu + Coulomb
            ftol=ftolP*1.d3
            Write(Iout,1003)
            CALL OptFF(Iout,ihessian,iprinthessian,iopt,IDA,
     1        Dist,dist2D,ftol,force)
            do i=1,9
              force(i)=forcep(i)
            enddo
            iopt=1
          endif
          CALL OptFF(Iout,ihessian,iprinthessian,iopt,IDA, ! vanilla Wu
     1      Dist,dist2D,ftolP,force)
          Iprint=0
c       extended Wu, 19 parameters
        else if(iopt.eq.3.or.iopt.eq.4.or.iopt.eq.5.or.iopt.eq.6)then
          CALL OptFF(Iout,ihessian,iprinthessian,iopt,IDA,
     1      Dist,dist2D,ftolP,force)
        endif
c       Compare structures
        CALL CompareStruct(Iout,IC3,Dist,DistStore)
        routine='MOVECM_2       '
        Write(Iout,1008) routine
        CALL MoveCM(Iout,Iprint,IAtom,mirror,isort,
     1   nosort,SP,Dist,DistCM,El)
        routine='DISTANCEMATRIX '
        Write(Iout,1008) routine
        CALL DistMatrix(Iout,Iprintf,0,0,
     1   Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)
        routine='DIAMETER       '
        Write(Iout,1008) routine
        CALL Diameter(Iout,Dist,distp)
        routine='RING           '
        Write(Iout,1008) routine
c  call ring again, this needs some reprogramming as ring duplicates some
c  stuff previously done, but is ok for now, as it takes not much time
        CALL Ring(Medges,MCon2,Iout,N5Ring,N6Ring,
     1   IC3,IVR3,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,DistMat)
      endif
      if(iprintf.ne.0) Call EdgeCoord(Iout,DIST,IC3)
      if(iprintf.ne.0.or.dualdist.ne.R6) Call RingCoord(Iout,
     1 dualdist,R6,SmallRingDist,DIST,N5Ring,N6Ring,N5MEM,N6MEM)
      if(iprintf.ne.0.or.ngaudiene.ne.0) Call Gaudiene(Iout,IC3,DIST)

C------------------XYZ-and-CC1-FILES------------------------------
C Print out Coordinates used as input for CYLview, VMD or other programs

C .xyz format
      if(iwext.ne.0) then
       if(iwext.eq.1) then
        nxyz=nxyz+1
        routine='PRINTCOORD     '
        Write(Iout,1008) routine
        Call WriteToFile(1,Iext,nxyz,ifind,Iout,IERROR1,IAtom,
     1   IC3,El,Dist,filenameout,xyzname,Namexyz,Endxyz,TEXTINPUT)
        if(IERROR1.eq.1) go to 9999
       endif

C .cc1 format
       if(iwext.eq.2) then
        ncc1=ncc1+1
        routine='PRINTCOORD     '
        Write(Iout,1008) routine
        Call WriteToFile(2,Iext,ncc1,ifind,Iout,IERROR1,IAtom,
     1   IC3,El,Dist,filenameout,cc1name,Namecc1,Endcc1,TEXTINPUT)
       endif

C .mol2 format
       if(iwext.eq.3) then
        nmol=nmol+1
        routine='PRINTCOORD     '
        Write(Iout,1008) routine
        Call WriteToFile(3,Iext,nmol,ifind,Iout,IERROR1,IAtom,
     1   IC3,El,Dist,filenameout,molname,Namemol,Endmol,TEXTINPUT)
       endif
      endif
      
C------------------VOLUME-----------------------------------------
C Calculate the volume
  97  if(novolume.eq.0) then
      if(imcs.eq.0) then
       routine='VOLUME         '
       Write(Iout,1008) routine
       CALL Volume(Iout,N5MEM,N6MEM,
     1  IDA,N5Ring,N6Ring,DIST,CRing5,CRing6,VolSphere,ASphere,
     2  Atol,VTol,Rmin5,Rmin6,Rmax5,Rmax6)!,filename)
      endif

C------------------MINCOVSPHERE-----------------------------------
C Calculate the minimum covering sphere and volumes
      routine='MINCOVSPHERE2  '
      Write(Iout,1008) routine
      CALL MinCovSphere2(Iout,imcs,SP,Dist,Rmin,Rmax,
     1 VolSphere,ASphere,Atol,VTol,distP,cmcs,rmcs,RVdWC)

C------------------MINDISTSPHERE----------------------------------
C Calculate the minimum distance sphere
      routine='MINDISTSPHERE  '
      Write(Iout,1008) routine
      CALL MinDistSphere(Iout,Dist,cmcs)

C------------------MAXINSPHERE------------------------------------
C Calculate the maximum inner sphere
      routine='MAXINSPHERE    '
      Write(Iout,1008) routine
      CALL MaxInSphere(Iout,imcs,Dist,cmcs,RVdWC)
      if(imcs.ne.0) go to 99 

C------------------PROJECTSPHERE----------------------------------
C Projecting vertices on minimum covering sphere
C  producing a spherical fullerene
      if(ipsphere.ne.0) then
       call ProjectSphere(ipsphere,Iout,IAtom,
     1 IC3,Dist,cmcs,rmcs,filename,El,TEXTINPUT)
      endif
C-----------------------------------------------------------------
      endif

C------------------------------------------------------------------
C  E N D   O F   3D   S T R U C T U R E   S E C T I O N
C------------------------------------------------------------------
C------------------GRAPH2D----------------------------------------
C Note: In the major restructuring, this routine may be subdivided
C   into whether 3D structure is required or not for Graph2D
C Calculate Schlegel diagram
      if(ISchlegel.ne.0) then
        routine='GRAPH2D        '
        Write(Iout,1008) routine
        if(number_vertices.lt.10000) then
          anglew=45.d0
          if(ISchlegel.eq.2) then
            if(ParamS.le.1.d0.or.ParamS.gt.8.9d1) then
              ParamS=anglew
              WRITE(Iout,1006) ParamS
            endif
          else
            ParamS=dabs(ParamS)
          endif
          CALL Graph2D(Iout,IS1,IS2,IS3,N5MEM,N6MEM,N5Ring,N6Ring,NRing,
     1   Iring,Ischlegel,ifs,ndual,labelvert,IC3,IDA,mdist,nhamcyc,Dist,
     1   ParamS,Rmin,TolX,scales,scalePPG,boost,CR,CRing5,CRing6,
     1   Symbol,filename)
        else
          Write(Iout,1007) number_vertices
        endif
      endif
C------------------END--------------------------------------------
C  E N D   O F   P R O G R A M
C-----------------------------------------------------------------
  99  if(loop-1) 100,101,102
 100  go to 9999
 101  iprev=0
      nloop=nloop+1
      WRITE(Iout,1005) nloop ! line of dashes
      go to 9 ! datain
 102  iprev=1
      nloop=nloop+1
      WRITE(Iout,1005) nloop ! line of dashes
      go to 9 ! datain
9999  call date_and_time(CDAT,CTIM,zone,values)
        WRITE(Iout,1004) Values(3),Values(2),Values(1),Values(5),
     1    Values(6),Values(7)
      CALL Timer(TIMEX)
      Hours=TIMEX/3.6d3
      WRITE(Iout,1009) TIMEX,Hours

C Formats 
C VERSION_NUMBER is set in the Makefile
 1000 FORMAT(
     1  1X,' ________________________________________________________ ',
     1 /1X,'|                                                        |',
     1 /1X,'|          P R O G R A M   F U L L E R E N E             |',
     1 /1X,'|    Fortran/C++ Program for the topological analysis    |',
     1 /1X,'|      of regular fullerenes (pentagons and hexagons)    |',
     1 /1X,'|    Written by Peter Schwerdtfeger, Lukas Wirz          |',
     1 /1X,'|         and James Avery with routines from             |',
     1 /1X,'|            Fowler, Manolopoulos and Babic              |',
     1 /1X,'|    Massey University,  Auckland,  New Zealand          |',
     1 /1X,'|    First version: 1.0                from 08/06/10     |',
     1 /1X,'|    This  version: ',VERSION_NUMBER,
     1                            ', last revision from 03/09/15     |',
     1 /1X,'|________________________________________________________|',
     1 //1X,'Date: ',I2,'/',I2,'/',I4,10X,'Time: ',I2,'h',I2,'m',I2,'s',
     1 /1X,'Limited to ',I6,' Atoms',
     1 /1X,'For citation when running this program use:',/1X,
     1 '1) P. Schwerdtfeger, L. Wirz, J. Avery, Program Fullerene - ',
     1 'A Software Package for Constructing and Analyzing Structures ',
     1 'of Regular Fullerenes (Version 4.4), J. Comput. Chem.,',
     1 'in press.',/1X,'If possible also cite the folowing two ',
     1 'references:',/1X,
     1 '2) P. W. Fowler, D. E. Manolopoulos, An Atlas of Fullerenes',
     1 ' (Dover Publ., New York, 2006).',/1X,
     1 '3) D. Babic, Nomenclature and Coding of Fullerenes,',
     1 ' J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).',/1X,
     1 'See the Manual for further literature and input instructions ',
     1 'concerning this program')
 1001 FORMAT(/1X,'Number of Atoms: ',I4,', and distance tolerance: ',
     1 F12.2,'%')
 1002 FORMAT(/1X,'Perform only minimum covering sphere calculations')
 1003 FORMAT(1X,'Pre-optimization using the Wu force field with ',
     1 'input parameter')
 1004 FORMAT(140(1H-),/1X,'DATE: ',I2,'/',I2,'/',I4,10X,
     1 'TIME: ',I2,'h',I2,'m',I2,'s')
 1005 FORMAT(140('='),/1X,'Loop ',I2)
 1006 FORMAT(/1X,'Angle for Schlegel diagram reset to ',
     1 F10.4,' degrees')
 1007 Format(1X,'2D Graph for such a large fullerene with ',I6,
     1 ' vertices is not meaningful ===> RETURN')
 1008 FORMAT(140('-'),/1x,'--> Enter Subroutine ',A15)
 1009 FORMAT(1x,'CPU Seconds: ',F15.2,', CPU Hours: ',F13.5)
 1010 FORMAT(1X,'Number of Hamiltonian cycles: ',I10)
 1011 FORMAT(3X,'(Add to this batches from previous cycles!)')
 1012 FORMAT(1X,'Connectivity field IC3 in input is errorneous: ',
     1 'taking only cartesian coordinates from input')
 1013 FORMAT(/1X,'End of file reached ==> Stop')
      STOP 
      END
