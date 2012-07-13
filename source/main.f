!                                                                           !
!           P R O G R A M     F U L L E R E N E                             !
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
!        Version 4 now incorporates C++ routines linked to the              !
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
      DIMENSION A(Nmax,Nmax),evec(Nmax),df(Nmax)
      real(8) force(ffmaxdim),forceP(ffmaxdim)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6),Iring(Mmax)
      DIMENSION Icon2(Nmax*Nmax),distP(Nmax),IDA(Nmax,Nmax)
      DIMENSION IATOM(Nmax),IC3(Nmax,3),Nring(Mmax)
      DIMENSION NringA(Emax),NringB(Emax)
      DIMENSION NringC(Emax),NringD(Emax)
      DIMENSION NringE(Emax),NringF(Emax)
      DIMENSION IDual(Mmax,Mmax),nSW(4,66),nFM(4,66),nYF(6,66),
     1 nWS(5,8)
      DIMENSION NEK(3,Nmax),JP(12)
      DIMENSION Symbol(Mmax)
      Real*4 TimeX
CG77  CHARACTER CDAT*9,CTIM*8
      CHARACTER CDAT*8,CTIM*10,Zone*5
      CHARACTER*1 Symbol
      CHARACTER*2 El(99)
      CHARACTER*13 routine
      CHARACTER*20 filename
      CHARACTER*31 xyzname
      CHARACTER*20 element
      Character*1 TEXTINPUT(nzeile)
      CHARACTER*3 GROUP
      Integer endzeile,Values(8)
      Integer MDist(Nmax,Nmax)
      type(c_ptr) :: graph, new_fullerene_graph
      DATA El/' H','HE','LI','BE',' B',' C',' N',' O',' F','NE','NA',
     1 'MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V','CR',
     2 'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',    
     3 'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD',   
     4 'IN','SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',  
     5 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF', 
     6 'TA',' W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',
     7 'AT','RN','FR','RA','AC','TH','PA',' U','NP','PU','AM','CM',   
     8 'BK','CF','ES'/                                               
      DATA Tol,anglew,Rdist/0.33d0,45.d0,1.391d0/
C     Van der Waals radius of carbon, adjusted approximately to the
C     solid-state results of P.A.Heiney et al., Phys. Rev. Lett. 66, 2911 (1991)
      DATA RVdWC/1.415d0/
      IN=5
      IOUT=6
      ilp=0
      iprev=0
      Do I=1,Nmax
       IAtom(I)=6
      Do J=1,Nmax
       IDA(I,J)=0
      enddo
      enddo

C  You might like to comment these 2 lines out 
C  (and same at the end of this routine) or substitute them with your
C  compiler specific option. Next two are g77 options
CG77    CALL Date(CDAT)
CG77    CALL Time(CTIM)
        call date_and_time(CDAT,CTIM,zone,values)
        TIMEX=0.d0
        CALL Timer(TIMEX)
C       WRITE(IOUT,1000) CDAT,CTIM,Nmax
        WRITE(IOUT,1000) Values(3),Values(2),Values(1),Values(5),
     1    Values(6),Values(7),Nmax
C  INPUT and setting parameters for running the subroutines
 9    routine='DATAIN       '
      leapspiral=0
      SWspiral=0
      Write(Iout,1008) routine
        CALL Datain(IN,IOUT,Nmax,MAtom,Icart,Iopt,iprintf,IHam,
     1  Ihueckel,KE,IPR,IPRC,ISchlegel,IS1,IS2,IS3,IER,istop,
     1  leap,leapGC,iupac,Ipent,iprintham,ISW,IGC1,IGC2,IV1,IV2,IV3,
     1  icyl,ichk,isonum,loop,mirror,ilp,IYF,IWS,nzeile,ifs,ndual,
     1  ParamS,TolX,R5,R6,Rdist,scales,scalePPG,ftolP,force,forceP,
     1  filename,TEXTINPUT)

      xyzname=trim(filename)//"-3D.xyz"

C  Stop if error in input
      If(IER.ne.0) go to 99
C  Only do isomer statistics
      if(istop.ne.0) go to 98

C Options for Input coordinates
      go to (10,20,30,30,30,30) Icart+1
C  Cartesian coordinates produced for Ih C60
   10 routine='COORDC20/60  '
      Write(Iout,1008) routine
      CALL CoordC20C60(Iout,MAtom,R5,R6,Dist)
      Do I=1,60
        IAtom(I)=6
      enddo
      Go to 40
C Input Cartesian coordinates for fullerenes

   20 if(icyl.ge.2) then
        Open(unit=7,file=xyzname,form='formatted')
        WRITE(Iout,1015) xyzname 
        Read(7,*) MAtom
        Read(7,1018) (TEXTINPUT(I),I=1,nzeile)
        endzeile=0
        do j=1,nzeile
          if(TEXTINPUT(j).ne.' ') endzeile=j
        enddo 
        WRITE(Iout,1017) MAtom,(TEXTINPUT(I),I=1,endzeile)
        Do J=1,MAtom
          Read(7,*,end=21) element,(Dist(I,J),I=1,3)
          Iatom(j)=6
        enddo
        close(unit=7)
        xyzname=trim(filename)//'-3D.new.xyz'
      else
        Do J=1,MAtom
          Read(IN,*,end=21) IAtom(J),(Dist(I,J),I=1,3)
        enddo
      endif
      Go to 40
   21 WRITE(Iout,1016)
      Go to 99
C Cartesian coordinates produced ring from spiral pentagon list
C currently using the Fowler-Manolopoulos algorithm to
C identify P-type eigenvectors and construct the 3D fullerene
   30 Ipent=1
      routine='COORDBUILD   '
      Write(Iout,1008) routine
      CALL CoordBuild(MAtom,IN,Iout,IDA,IDual,
     1 Icart,IV1,IV2,IV3,IGC1,IGC2,isonum,IPRC,ihueckel,JP,
     1 iprev,A,evec,df,Dist,Dist2D,distp,Rdist,GROUP)
      Do I=1,Matom
       IAtom(I)=6
      enddo
      CALL Chiral(Iout,GROUP)

   40 WRITE(Iout,1001) MAtom,TolX*100.d0

C Some general infos on isomers and spiral routine
C of Fowler and Manolopoulos. Set parameter IPR for independent
C pentagon rule as full list beyond C60 is computer time 
C intensive
  98  routine='ISOMERS      '
      Write(Iout,1008) routine
      CALL Isomers(MAtom,IPR,IOUT,
     1 maxit,iprintham,ichk,IDA,A,trim(filename)//".chkpnt")
      if(istop.ne.0) go to 99

C Move carbon cage to Atomic Center
  999 Group='   '
      routine='MOVECM       '
      Write(Iout,1008) routine
      Iprint=1
      Call MoveCM(Matom,Iout,Iprint,IAtom,mirror,
     1 Dist,DistCM,El)

C Calculate largest and smallest atom-to-atom diameters
C Also get moment of inertia (to be implemented)
      routine='DIAMETER     '
      Write(Iout,1008) routine
      CALL Diameter(MAtom,IOUT,Dist,distp)

C Calculate the distance Matrix and print out distance Matrix
      routine='DISTMATRIX   '
      Write(Iout,1008) routine
      CALL Distmatrix(MAtom,IOUT,iprintf,Iopt,
     1 Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)

C Establish Connectivities
      routine='CONNECT      '
      Write(Iout,1008) routine
      CALL Connect(MCon2,MAtom,Ipent,IOUT,
     1 Icon2,IC3,IDA,TolX,DistMat,Rmin)

C Hueckel matrix and eigenvalues
      if(ipent.eq.0) then
      routine='HUECKEL      '
      Write(Iout,1008) routine
      CALL Hueckel(MAtom,IOUT,IC3,ihueckel,IDA,A,evec,df)
      endif

C Produce the nth leapfrog of the fullerene
      if(leap.gt.0.or.leapGC.gt.0) then
      routine='Leapfrog'
      Write(Iout,1008) routine
      CALL GoldbergCoxeter(MAtom,Iout,leap,leapGC,IGC1,IGC2,
     1 ihueckel,LeapErr,IDA,A,evec,df,Dist,Dist2D,distp,Rdist)
      leap=0
      leapGC=0
      ipent=1
      leapspiral=1
      if(MAtom.gt.100) IHam=0
      if(LeapErr.eq.0) go to 999
      endif

c$$$      graph = new_fullerene_graph(Nmax,MAtom,IDA)
c$$$      call tutte_layout(graph,Dist2D)
c$$$      call set_layout2d(graph,Dist2D)
c$$$
c$$$      call lukas_edges(graph,MAtom)
c$$$      call lukas_corners(graph,MAtom)
c$$$      call lukas_dihedrals(graph,MAtom)

C Generate IUPAC name and locate Hamiltonian cycles. 
C Routine written by D. Babic. Note routine
C is called only if IPR>0 as computer time is extensive beyond
C C100 (PN-hard problem). Last routine uses the adjaceny matrix
C to calculate the number of all distinct paths between 
C adjacent vertices
      routine='HAMILTON     '
      Write(Iout,1008) routine
       maxiter=maxit
      if(IHam.gt.1.and.IHam.le.9) then
       maxiter=10**IHam
      endif
      if(IHam.ne.0.and.ISW.eq.0) then
       if(iupac.ne.0) then
         CALL Hamilton(MAtom,Iout,iprintf,maxiter,IC3)
       else
         CALL HamiltonCyc(MAtom,maxiter,Iout,nbatch,IDA,Nhamilton)
         WRITE(Iout,1010) Nhamilton
         if(nbatch.ne.0) WRITE(Iout,1014)
       endif
      endif
      CALL Paths(MAtom,IOUT,IDA,A,evec,df)

C Establish all closed ring systems
      routine='RING         '
      Write(Iout,1008) routine
      CALL Ring(MCon2,MAtom,IOUT,
     1 N5Ring,N6Ring,IC3,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,
     1 DistMat)

C Optimize Geometry through force field method
c we check for ISW because the coordinates shouldn't be optimized before
c a stone wales (or any other transformation) is done
      If(Iopt.ne.0.and.ISW.eq.0.and.iyf.eq.0.and.iws.eq.0) then
        routine='OPTFF        '
        ftol=ftolP
        Write(Iout,1008) routine
        if(Iopt.eq.1 .or. Iopt.eq.2) then ! vanilla Wu or Wu + Coulomb
          if(Iopt.eq.2) then ! Wu + Coulomb
            ftol=ftolP*1.d3
            Write(Iout,1003)
            CALL OptFF(MAtom,Iout,IDA,N5Ring,N6Ring,
     1        N5MEM,N6MEM,Dist,Rdist,ftol,force,iopt)
          endif
          iopt=1
          CALL OptFF(MAtom,Iout,IDA,N5Ring,N6Ring, ! vanilla Wu
     1      N5MEM,N6MEM,Dist,Rdist,ftolP,forceP,iopt)
          Iprint=0
        else if(Iopt.eq.3) then ! extended Wu, 18 parameters
          CALL OptFF(MAtom,Iout,IDA,N5Ring,N6Ring,
     1      N5MEM,N6MEM,Dist,Rdist,ftolP,forceP,iopt)
        endif
        Call MoveCM(Matom,Iout,Iprint,IAtom,mirror,
     1   Dist,DistCM,El)
        routine='DISTMATRIX   '
        Write(Iout,1008) routine
        CALL Distmatrix(MAtom,IOUT,Iprintf,0,
     1   Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)
        routine='DIAMETER     '
        Write(Iout,1008) routine
        CALL Diameter(MAtom,IOUT,Dist,distp)
      endif

C Print out Coordinates used as input for CYLview
      if(icyl.le.2.and.ISW.eq.0) then
      Open(unit=3,file=xyzname,form='formatted')
      routine='PRINTCOORD   '
      Write(Iout,1008) routine
      WRITE(Iout,1002) xyzname 
      endzeile=0
       do j=1,nzeile
        if(TEXTINPUT(j).ne.' ') endzeile=j
       enddo
      if(MAtom.lt.100) WRITE(3,1011) MAtom,MAtom,
     1 (TEXTINPUT(I),I=1,endzeile)
      if(MAtom.ge.100.and.MAtom.lt.1000) 
     1 WRITE(3,1012) MAtom,MAtom,(TEXTINPUT(I),I=1,endzeile)
      if(MAtom.ge.1000.and.MAtom.lt.10000) 
     1 WRITE(3,1013) MAtom,MAtom,(TEXTINPUT(I),I=1,endzeile)
      Do J=1,MAtom
      IM=IAtom(J)      
      Write(3,1007) El(IM),(Dist(I,J),I=1,3)
      enddo
      Close(unit=3)
      endif

C Rings
      routine='RING         '
      Write(Iout,1008) routine
      CALL Ring(MCon2,MAtom,IOUT,
     1 N5Ring,N6Ring,IC3,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,
     1 DistMat)

C Analyze ring connections
      routine='RINGC        '
      Write(Iout,1008) routine
      CALL RingC(Matom,Iout,iprintf,
     1 N5MEM,N6MEM,N5Ring,N6Ring,NRing,Iring5,Iring6,Iring56,NringA,
     1 NringB,NringC,NringD,NringE,NringF,numbersw,nSW,n565,NEK,
     1 numberFM,nFM,numberYF,nYF,numberWS,nWS,DIST,CRing5,CRing6)

C Perform Stone-Wales transformation
      if(ISW.ne.0) then
      routine='STONE-WALES  '
      Write(Iout,1008) routine
      CALL StoneWalesTrans(Matom,IN,Iout,numbersw,nSW,
     1 ihueckel,IDA,N6MEM,IC3,A,evec,df,Dist,Dist2D,distp,Rdist)
      ISW=0
      ipent=1
      SWspiral=1
      go to 999
      endif

C Perform Endo-Kroto 2-vertex insertion
      if(KE.ne.0) then
      routine='ENDO-KROTO   '
      Write(Iout,1008) routine
      CALL EndoKrotoTrans(Matom,IN,Iout,n565,NEK,ihueckel,
     1 IDA,N5MEM,N6MEM,A,evec,df,Dist,Dist2D,distp,Rdist)
      KE=0
      ipent=1
      SWspiral=1
      go to 999
      endif

C Perform Yoshida-Fowler 4-or 6-vertex insertion
      if(IYF.ne.0) then
      routine='YOSHIDAFOWLER'
      Write(Iout,1008) routine
      if(IYF.le.2) then
      CALL YoshidaFowler4(Matom,IN,Iout,JERR,numberFM,IYF,nFM,
     1 ihueckel,IDA,N5MEM,N6MEM,A,evec,df,Dist,Dist2D,distp,Rdist)
      else
      IYF=IYF-2
      CALL YoshidaFowler6(Matom,IN,Iout,JERR,numberYF,IYF,nYF,
     1 ihueckel,IDA,N6MEM,IC3,A,evec,df,Dist,Dist2D,distp,Rdist)
      endif
      IYF=0
      ipent=1
      SWspiral=1
      if(JERR.eq.0) go to 999
      endif

C Perform Wirz-Schwerdtfeger 6-vertex 6-55-55 insertion
      if(IWS.ne.0) then
      routine='WIRZSCHWERD  '
      Write(Iout,1008) routine
      CALL WirzSchwerd(Matom,IN,Iout,JERR,numberWS,IWS,
     1 nWS,ihueckel,IDA,N5MEM,N6MEM,IC3,
     1 A,evec,df,Dist,Dist2D,distp,Rdist)
      IWS=0
      ipent=1
      SWspiral=1
      if(JERR.eq.0) go to 999
      endif

C Now produce clockwise spiral ring pentagon count a la Fowler and Manolopoulos
      if(ipent.eq.0.or.leapspiral.ne.0.or.SWspiral.ne.0) then
      routine='SPIRALSEARCH '
      Write(Iout,1008) routine
      CALL SpiralSearch(Nspirals,MAtom,Iout,Iring5,
     1 Iring6,Iring56,NringA,NringB,NringC,NringD,NringE,NringF,JP,
     1 GROUP)
C Determine if fullerene is chiral
      CALL Chiral(Iout,GROUP)
      endif
C Topological Indicators
      Call TopIndicators(Matom,Iout,IDA,Mdist)

C Calculate the volume
      routine='VOLUME       '
      Write(Iout,1008) routine

      CALL Volume(Matom,Iout,N5MEM,N6MEM,
     1 IDA,N5Ring,N6Ring,DIST,CRing5,CRing6,VolSphere,ASphere,
     2 Atol,VTol,Rmin5,Rmin6,Rmax5,Rmax6,filename)

C Calculate the minimum covering sphere and volumes
      routine='MINCOVSPHERE2'
      Write(Iout,1008) routine
      CALL MinCovSphere2(MAtom,IOUT,Dist,Rmin,Rmax,
     1 VolSphere,ASphere,Atol,VTol,distP,cmcs,rmcs,RVdWC)

C Calculate the minimum distance sphere
      routine='MINDISTSPHERE'
      Write(Iout,1008) routine
      CALL MinDistSphere(MAtom,IOUT,Dist,cmcs)

C Calculate the maximum inner sphere
      routine='MAXINSPHERE'
      Write(Iout,1008) routine
      CALL MaxInSphere(MAtom,IOUT,Dist,cmcs,RVdWC)

C Calculate Schlegel diagram
      if(ISchlegel.ne.0) then
      routine='GRAPH2D      '
      Write(Iout,1008) routine
      if(ISchlegel.eq.2) then
       if(ParamS.le.1.d0.or.ParamS.gt.8.9d1) then
       ParamS=anglew
       WRITE(IOUT,1006) ParamS
       endif
      else
       ParamS=dabs(ParamS)
      endif
      CALL Graph2D(MAtom,IOUT,IS1,IS2,IS3,N5MEM,N6MEM,N5Ring,N6Ring,
     1 NRing,Iring,Ischlegel,ifs,ndual,IC3,IDA,Mdist,Dist,ParamS,Rmin,
     1 TolX,scales,scalePPG,CR,CRing5,CRing6,Symbol,filename)
      endif

C  E N D   O F   P R O G R A M
  99  if(loop-1) 100,101,102
 100  go to 9999
 101  iprev=0
      WRITE(IOUT,1019) ! line of dashes
      go to 9 ! datain
 102  iprev=1
      WRITE(IOUT,1019) 
      go to 9 
CG77 99  CALL TIME(CTIM)
CG77 99  CALL TIME(CTIM)
9999  call date_and_time(CDAT,CTIM,zone,values)
        WRITE(IOUT,1004) Values(3),Values(2),Values(1),Values(5),
     1    Values(6),Values(7)
      CALL Timer(TIMEX)
      Hours=TIMEX/3.6d3
      WRITE(IOUT,1009) TIMEX,Hours
C Formats 
 1000 FORMAT(
     1  1X,' ________________________________________________________ ',
     1 /1X,'|                                                        |',
     1 /1X,'|          P R O G R A M   F U L L E R E N E             |',
     1 /1X,'|    Fortran/C++ Program for the topological analysis    |',
     1 /1X,'|      of regular fullerenes (pentagons and hexagons)    |',
     1 /1X,'|    Written by Peter Schwerdtfeger and James Avery      |',
     1 /1X,'|      with routines from Fowler, Manolopoulos and Babic |',
     1 /1X,'|    Massey University,  Auckland,  New Zealand          |',
     1 /1X,'|    First version: 1.0:               from 08/06/10     |',
     1 /1X,'|    This  version: 4.0, last revision from 12/07/12     |',
     1 /1X,'|________________________________________________________|',
CG77 1 /1X,'DATE: ',A9,10X,'TIME: ',A8,/1X,'Limited to ',I6,' Atoms',
     1 //1X,'Date: ',I2,'/',I2,'/',I4,10X,'Time: ',I2,'h',I2,'m',I2,'s',
     1 /1X,'Limited to ',I6,' Atoms',
     1 /1X,'For citation when running this program use:',/1X,
     1 '1) P. Schwerdtfeger, J. Avery, Topological Aspects of ',
     1 'Fullerenes - A Fortran Program',/9X,'(Version 4.0, Massey ',
     1 'University Albany, Auckland, New Zealand, 2012).',/1X,
     1 '2) P. W. Fowler, D. E. Manolopoulos, An Atlas of Fullerenes',
     1 ' (Dover Publ., New York, 2006).',/1X,
     1 '3) D. Babic, Nomenclature and Coding of Fullerenes,',
     1 ' J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).',/1X,
     1 'See README file for further literature and input instructions ',
     1 'concerning this program')
 1001 FORMAT(/1X,'Number of Atoms: ',I4,', and distance tolerance: ',
     1 F12.2,'%')
 1002 FORMAT(/1X,'Input coordinates to be used for plotting program',
     1 ' CYLVIEW, PYMOL or AVOGADRO',/1X,'Output written into ',A20)
CG77 1004 FORMAT(1X,140(1H-),/1X,6HTIME: ,A8)
 1003 FORMAT(1X,'Pre-optimization using the Wu force field with ',
     1 'input parameter')
 1004 FORMAT(140(1H-),/1X,'DATE: ',I2,'/',I2,'/',I4,10X,
     1 'TIME: ',I2,'h',I2,'m',I2,'s')
 1006 FORMAT(/1X,'Angle for Schlegel diagram reset to ',
     1 F10.4,' degrees')
 1007 FORMAT(A2,6X,3(F15.6,2X))
 1008 FORMAT(140('-'),/1x,'--> Enter Subroutine ',A13)
 1009 FORMAT(1x,'CPU Seconds: ',F15.2,', CPU Hours: ',F13.5)
 1010 FORMAT(1X,'Number of Hamiltonian cycles: ',I10)
 1011 FORMAT(I5,/,'C',I2,'/  ',132A1)
 1012 FORMAT(I5,/,'C',I3,'/  ',132A1)
 1013 FORMAT(I5,/,'C',I4,'/  ',132A1)
 1014 FORMAT(3X,'(Add to this batches from previous cycles!)')
 1015 FORMAT(/1X,'Read coordinates from xyz file: ',A20)
 1016 FORMAT(/1X,'End of file reached ==> Stop')
 1017 FORMAT(1X,'Number of Atoms: ',I5,/1X,132A1)
 1018 FORMAT(132A1)
 1019 FORMAT(140('='))
      STOP 
      END

      SUBROUTINE TIMER(TIMEX)
      Real TA(2)
      Call DTIME(TA,time)
      TIMEX=TIME
      RETURN
      END
