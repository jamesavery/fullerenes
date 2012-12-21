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
!        Version 4.3 now incorporates C++ routines linked to the            !
!       original Fortran program using much improved algorithms.            !
!---------------------------------------------------------------------------!

      PROGRAM Fullerene
      use iso_c_binding
      use config

      IMPLICIT REAL*8 (A-H,O-Z)
C    Set the dimensions for the distance matrix
      parameter (nzeile=132)
      DIMENSION CRing5(3,Mmax),CRing6(3,Mmax),cmcs(3),CR(3,Mmax)
      DIMENSION DistMat(NmaxL),Dist(3,Nmax),DistCM(3),Dist2D(2,Nmax)
      DIMENSION DistStore(3,Nmax)
      DIMENSION A(Nmax,Nmax),evec(Nmax),df(Nmax)
      real(8) force(ffmaxdim), forceP(ffmaxdim)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6),Iring(Mmax)
      DIMENSION distP(Nmax),IDA(Nmax,Nmax)
      DIMENSION IATOM(Nmax),IC3(Nmax,3),Nring(Mmax),IVR3(Nmax,3)
      DIMENSION NringA(Emax),NringB(Emax)
      DIMENSION NringC(Emax),NringD(Emax)
      DIMENSION NringE(Emax),NringF(Emax)
      DIMENSION IDual(Mmax,Mmax),nSW(4,66),nFM(4,66),nYF(6,66),nBF(5,66)
      DIMENSION NEK(3,Nmax),JP(12)
      DIMENSION Symbol(Mmax)
      Real*4 TimeX
      CHARACTER CDAT*8,CTIM*10,Zone*5
      CHARACTER*1 Symbol
      CHARACTER*2 El(99)
      CHARACTER*7 Namecc1,Namexyz
      CHARACTER*4 Endcc1,Endxyz
      CHARACTER*15 routine
      CHARACTER*50 filename
      CHARACTER*50 filenameout
      CHARACTER*50 xyzname
      CHARACTER*50 cc1name
      CHARACTER*20 element
      Character*1 TEXTINPUT(nzeile)
      CHARACTER*3 GROUP
      Integer endzeile,Values(8)
      integer istop
      Logical lexist
      integer mdist(nmax,nmax)

C Set parameters
      DATA El/' H','HE','LI','BE',' B',' C',' N',' O',' F','NE','NA',
     1 'MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V','CR',
     2 'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',    
     3 'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD',   
     4 'IN','SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',  
     5 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF', 
     6 'TA',' W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',
     7 'AT','RN','FR','RA','AC','TH','PA',' U','NP','PU','AM','CM',   
     8 'BK','CF','ES'/                                               
      Namecc1='-3D.cc1'
      Namexyz='-3D.xyz'
      Endcc1='.cc1'
      Endxyz='.xyz'
      IN=5
      Iout=6
      nloop=0
      nxyz=0
      ncc1=0
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

C Get time and date
      CALL date_and_time(CDAT,CTIM,zone,values)
      TIMEX=0.d0
      CALL Timer(TIMEX)
      WRITE(Iout,1000) Values(3),Values(2),Values(1),Values(5),
     1  Values(6),Values(7),Nmax

C------------------DATAIN------------------------------------------
C  INPUT and setting parameters for running the subroutines
 9    routine='DATAIN         '
      Group='   '
      isort=0
      leapspiral=0
      SWspiral=0
      Write(Iout,1008) routine
        CALL Datain(IN,Iout,Nmax,Icart,Iopt,iprintf,IHam,
     1  nohueckel,KE,IPR,IPRC,ISchlegel,IS1,IS2,IS3,IER,istop,
     1  leap,leapGC,iupac,Ipent,iprintham,IGC1,IGC2,IV1,IV2,IV3,
     1  icyl,ichk,isonum,loop,mirror,ilp,ISW,IYF,IBF,nzeile,ifs,
     1  ipsphere,ndual,nosort,ispsearch,novolume,ihessian,isearch,
     1  iprinthessian,
     1  ParamS,TolX,R5,R6,Rdist,rvdwc,scales,scalePPG,ftolP,scaleRad,
     1  force,forceP,boost,filename,filenameout,TEXTINPUT)

C  Stop if isomer closest to icosahedral is searched for
      if(isearch.ne.0) then
        istop=1
        go to 98
      endif
C  Stop if error in input
      If(IER.ne.0) go to 99
C  Only do isomer statistics
      if(istop.ne.0) go to 98

C------------------Coordinates-------------------------------------
C Options for Input coordinates
      go to (10,20,30,30,30,30,30,98) Icart+1

C  Cartesian coordinates produced for Ih C20 or C60
   10 routine='COORDC20/60    '
      Write(Iout,1008) routine
      CALL CoordC20C60(Iout,R5,R6,Dist)
      Do I=1,60
        IAtom(I)=6
      enddo
      Go to 40

C Input Cartesian coordinates for fullerenes
   20 if(icyl.eq.2.or.icyl.eq.3.or.icyl.eq.5) then
        if(icyl.eq.5) then
         cc1name=trim(filename)//".cc1"
         inquire(file=cc1name,exist=lexist)
          if(lexist.neqv..True.) then
            Write(Iout,1023) cc1name
            stop
          endif
         Open(unit=7,file=cc1name,form='formatted')
         WRITE(Iout,1021) cc1name 
         Read(7,*) number_vertices
         Do J=1,number_vertices
           Read(7,*,end=23) element,JJ,(Dist(I,J),I=1,3)
           Iatom(j)=6
         enddo
   23    close(unit=7)
         cc1name=trim(filename)//'-3D.new.xyz'
        else
         xyzname=trim(filename)//".xyz"
         Open(unit=7,file=xyzname,form='formatted')
         WRITE(Iout,1015) xyzname 
         Read(7,*) number_vertices
         Read(7,1018) (TEXTINPUT(I),I=1,nzeile)
         endzeile=0
         do j=1,nzeile
           if(TEXTINPUT(j).ne.' ') endzeile=j
         enddo 
         WRITE(Iout,1017) number_vertices,(TEXTINPUT(I),I=1,endzeile)
         Do J=1,number_vertices
           Read(7,*,end=22) element,(Dist(I,J),I=1,3)
           Iatom(j)=6
         enddo
   22    close(unit=7)
         xyzname=trim(filename)//'-3D.new.xyz'
        endif

      else

       Do J=1,number_vertices
        Read(IN,*,end=21) IAtom(J),(Dist(I,J),I=1,3)
       enddo
      endif
      Go to 40
   21 WRITE(Iout,1016)
      Go to 99

C Cartesian coordinates produced ring from spiral pentagon list
C or from adjacency matrix. Uses the Fowler-Manolopoulos algorithm 
C with P-type eigenvectors or the Tutte algorithm to construct 
C the 3D fullerene
   30 Ipent=1
      routine='COORDBUILD     '
      Write(Iout,1008) routine
      CALL CoordBuild(IN,Iout,IDA,IDual,
     1 Icart,IV1,IV2,IV3,IGC1,IGC2,isonum,IPRC,nohueckel,JP,
     1 iprev,ihalma,A,evec,df,Dist,Dist2D,distp,Rdist,scaleRad,
     1 GROUP,filename)
      Do I=1,number_vertices
        IAtom(I)=6
      enddo

   40 WRITE(Iout,1001) number_vertices,TolX*100.d0

C------------------ISOMERS-----------------------------------------
C Some general infos on isomers and spiral routine
C of Fowler and Manolopoulos. Set parameter IPR for independent
C pentagon rule as full list beyond C60 is computer time 
C intensive
  98  routine='ISOMERS        '
      Write(Iout,1008) routine
      CALL Isomers(IPR,isearch,In,Iout,iprintham,ichk,IDA,A,filename)
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
        CALL GoldbergCoxeter(Iout,leap,leapGC,IGC1,IGC2,
     1   nohueckel,LeapErr,IDA,A,evec,df,Dist,Dist2D,distp,Rdist,
     1   scaleRad)
        leap=0
        leapGC=0
        ipent=1
        leapspiral=1
        if(number_vertices.gt.100) IHam=0
        if(LeapErr.eq.0) go to 999 ! moveCM
      endif

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
        if(iupac.ne.0) then
          CALL Hamilton(Iout,iprintf,maxiter,IC3)
        else
          CALL HamiltonCyc(maxiter,Iout,nbatch,IDA,Nhamilton)
          WRITE(Iout,1010) Nhamilton
          if(nbatch.ne.0) WRITE(Iout,1014)
        endif
      endif
      CALL Paths(Iout,iprintf,IDA,A,evec,df)

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
     1 IC3,IVR3,nEK,nSW,nFM,nYF,nBF,DIST,CRing5,CRing6)
C     Print edge coordinates (barycenter)
      if(iprintf.ne.0) Call EdgeCoord(Iout,DIST,IC3)

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
      ispcount=0
      if(ispsearch.gt.1) ispcount=1
      if((ipent.eq.0.or.leapspiral.ne.0.or.SWspiral.ne.0.
     1   or.Icart.eq.6.or.Icart.eq.7.or.ihalma.eq.1).or.
     1   ispsearch.ne.0) then
        routine='SPIRALSEARCH   '
        Write(Iout,1008) routine
        CALL SpiralSearch(Iout,Iring5,Iring6,Iring56,NringA,NringB,
     1   NringC,NringD,NringE,NringF,JP,GROUP,ispcount)
      endif

C--------------TOPOLOGICAL INDICATORS-----------------------------
        routine='TOPOLOINDICATOR'
        Write(Iout,1008) routine
C Topological Indicators
      CALL TopIndicators(Iout,IDA,mdist)
C Check if vertex number allows for icosahedral fullerenes
      Call IcoFullDetect(Iout)
C Determine if fullerene is chiral
      CALL Chiral(Iout,GROUP)
C Produce perfect matchings (Kekule structures) and analyze
c      CALL PerfectMatching(Iout,IDA)

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
        else if(Iopt.eq.3 .or. iopt.eq.4) then ! extended Wu, 19 parameters
          CALL OptFF(Iout,ihessian,iprinthessian,iopt,IDA,
     1      Dist,dist2D,ftolP,force)
        endif
c       Compare structures
        CALL CompareStruct(Iout,IDA,Dist,DistStore)
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

C------------------XYZ-and-CC1-FILES------------------------------
C Print out Coordinates used as input for CYLview, VMD or other programs

C xyz format
      if(icyl.le.2) then
      nxyz=nxyz+1
      Call FileMod(filenameout,xyzname,Namexyz,Endxyz,nxyz,ifind)
        if(ifind.ne.0) then
         Write(Iout,1022)
         go to 9999
        endif
        Open(unit=3,file=xyzname,form='formatted')
        routine='PRINTCOORD     '
        Write(Iout,1008) routine
        WRITE(Iout,1002) xyzname 
        endzeile=0
        do j=1,nzeile
          if(TEXTINPUT(j).ne.' ') endzeile=j
        enddo
        if(number_vertices.lt.100) WRITE(3,1011)
     1    number_vertices,number_vertices, (TEXTINPUT(I),I=1,endzeile)
        if(number_vertices.ge.100.and.number_vertices.lt.1000) 
     1    WRITE(3,1012) number_vertices,number_vertices,
     1    (TEXTINPUT(I),I=1,endzeile)
        if(number_vertices.ge.1000.and.number_vertices.lt.10000) 
     1    WRITE(3,1013) number_vertices,number_vertices,
     1    (TEXTINPUT(I),I=1,endzeile)
        if(number_vertices.ge.10000) 
     1    WRITE(3,1020) number_vertices,number_vertices,
     1    (TEXTINPUT(I),I=1,endzeile)
        Do J=1,number_vertices
          IM=IAtom(J)      
          Write(3,1007) El(IM),(Dist(I,J),I=1,3)
        enddo
        Close(unit=3)
      endif

C cc1 format
      if(icyl.ge.4) then
C     Name handling
      ncc1=ncc1+1
      Call FileMod(filenameout,cc1name,Namecc1,Endcc1,ncc1,ifind)
        if(ifind.ne.0) then
         Write(Iout,1022)
         go to 9999
        endif
       Open(unit=3,file=cc1name,form='formatted')
        routine='PRINTCOORD     '
        Write(Iout,1008) routine
        WRITE(Iout,1002) cc1name
        if(number_vertices.lt.100) WRITE(3,1025) number_vertices
        if(number_vertices.ge.100.and.number_vertices.lt.1000)
     1    WRITE(3,1026) number_vertices
        if(number_vertices.ge.1000.and.number_vertices.lt.10000)
     1    WRITE(3,1027) number_vertices
        if(number_vertices.ge.10000) WRITE(3,1028) number_vertices
        Do J=1,number_vertices
          IM=IAtom(J)
          Write(3,1005) El(IM),J,(Dist(I,J),I=1,3),(IC3(J,I),I=1,3)
        enddo
        Close(unit=3)
      endif
      
      if(novolume.eq.0) then
C------------------VOLUME-----------------------------------------
C Calculate the volume
      routine='VOLUME         '
      Write(Iout,1008) routine
      CALL Volume(Iout,N5MEM,N6MEM,
     1 IDA,N5Ring,N6Ring,DIST,CRing5,CRing6,VolSphere,ASphere,
     2 Atol,VTol,Rmin5,Rmin6,Rmax5,Rmax6)!,filename)

C------------------MINCOVSPHERE-----------------------------------
C Calculate the minimum covering sphere and volumes
      routine='MINCOVSPHERE2  '
      Write(Iout,1008) routine
      CALL MinCovSphere2(Iout,SP,Dist,Rmin,Rmax,
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
      CALL MaxInSphere(Iout,Dist,cmcs,RVdWC)

C------------------PROJECTSPHERE----------------------------------
C Projecting vertices on minimum covering sphere
C  producing a spherical fullerene
      if(ipsphere.ne.0) then
       call ProjectSphere(ipsphere,Iout,IAtom,nzeile,
     1 IC3,Dist,cmcs,rmcs,filename,El,TEXTINPUT)
      endif
C-----------------------------------------------------------------
      endif

C------------------GRAPH2D----------------------------------------
C Calculate Schlegel diagram
      if(ISchlegel.ne.0) then
        anglew=45.d0
        routine='GRAPH2D        '
        Write(Iout,1008) routine
        if(ISchlegel.eq.2) then
          if(ParamS.le.1.d0.or.ParamS.gt.8.9d1) then
            ParamS=anglew
            WRITE(Iout,1006) ParamS
          endif
        else
          ParamS=dabs(ParamS)
        endif
        CALL Graph2D(Iout,IS1,IS2,IS3,N5MEM,N6MEM,N5Ring,N6Ring,
     1   NRing,Iring,Ischlegel,ifs,ndual,IC3,IDA,mdist,Dist,ParamS,Rmin,
     1   TolX,scales,scalePPG,boost,CR,CRing5,CRing6,Symbol,filename)
      endif
C------------------END--------------------------------------------
C  E N D   O F   P R O G R A M
  99  if(loop-1) 100,101,102
 100  go to 9999
 101  iprev=0
      nloop=nloop+1
      WRITE(Iout,1019) nloop ! line of dashes
      go to 9 ! datain
 102  iprev=1
      nloop=nloop+1
      WRITE(Iout,1019) nloop ! line of dashes
      go to 9 ! datain
9999  call date_and_time(CDAT,CTIM,zone,values)
        WRITE(Iout,1004) Values(3),Values(2),Values(1),Values(5),
     1    Values(6),Values(7)
      CALL Timer(TIMEX)
      Hours=TIMEX/3.6d3
      WRITE(Iout,1009) TIMEX,Hours

C Formats 
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
     1 /1X,'|    First version: 1.0:               from 08/06/10     |',
     1 /1X,'|    This  version: 4.3.1, last revision from 21/12/12   |',
     1 /1X,'|________________________________________________________|',
     1 //1X,'Date: ',I2,'/',I2,'/',I4,10X,'Time: ',I2,'h',I2,'m',I2,'s',
     1 /1X,'Limited to ',I6,' Atoms',
     1 /1X,'For citation when running this program use:',/1X,
     1 '1) P. Schwerdtfeger, L. Wirz, J. Avery, Program Fullerene - ',
     1 'A Software Package for Constructing and Analyzing Structures ',
     1 'of Regular Fullerenes (Version 4.3), submitted to J. Comput. ',
     1 'Chem.',/1X,
     1 '2) P. W. Fowler, D. E. Manolopoulos, An Atlas of Fullerenes',
     1 ' (Dover Publ., New York, 2006).',/1X,
     1 '3) D. Babic, Nomenclature and Coding of Fullerenes,',
     1 ' J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).',/1X,
     1 'See the Manual for further literature and input instructions',
     1 'concerning this program')
 1001 FORMAT(/1X,'Number of Atoms: ',I4,', and distance tolerance: ',
     1 F12.2,'%')
 1002 FORMAT(/1X,'Input coordinates to be used for plotting program',
     1 ' CYLVIEW, PYMOL or AVOGADRO',/1X,'Output written into ',A31)
 1003 FORMAT(1X,'Pre-optimization using the Wu force field with ',
     1 'input parameter')
 1004 FORMAT(140(1H-),/1X,'DATE: ',I2,'/',I2,'/',I4,10X,
     1 'TIME: ',I2,'h',I2,'m',I2,'s')
 1005 FORMAT(A2,I5,3F12.6,'    2',3I5)
 1006 FORMAT(/1X,'Angle for Schlegel diagram reset to ',
     1 F10.4,' degrees')
 1007 FORMAT(A2,6X,3(F15.6,2X))
 1008 FORMAT(140('-'),/1x,'--> Enter Subroutine ',A15)
 1009 FORMAT(1x,'CPU Seconds: ',F15.2,', CPU Hours: ',F13.5)
 1010 FORMAT(1X,'Number of Hamiltonian cycles: ',I10)
 1011 FORMAT(I5,/,'C',I2,'/  ',132A1)
 1012 FORMAT(I5,/,'C',I3,'/  ',132A1)
 1013 FORMAT(I5,/,'C',I4,'/  ',132A1)
 1014 FORMAT(3X,'(Add to this batches from previous cycles!)')
 1015 FORMAT(/1X,'Read coordinates from xyz file: ',A60)
 1016 FORMAT(/1X,'End of file reached ==> Stop')
 1017 FORMAT(1X,'Number of Atoms: ',I5,/1X,132A1)
 1018 FORMAT(132A1)
 1019 FORMAT(140('='),/1X,'Loop ',I2)
 1020 FORMAT(I8,/,'C',I8,'/  ',132A1)
 1021 FORMAT(/1X,'Read coordinates from cc1 file: ',A60)
 1022 FORMAT(/1X,'You try to write into the database filesystem',
     1 ' which is not allowed  ===>  ABORT')
 1023 Format(1X,'Filename ',A50,' in database not found ==> ABORT')
 1025 FORMAT(I2)
 1026 FORMAT(I3)
 1027 FORMAT(I4)
 1028 FORMAT(I8)
      STOP 
      END
