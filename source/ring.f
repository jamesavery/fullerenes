      SUBROUTINE Ring(NAtom,Nedges,Nfaces,natomL,Natom2,MCon2,MAtom,
     1 IOUT,Ncount5,Ncount6,IC3,Icon2,N5MEM,N6MEM,Rmin5,
     1 Rmin6,Rmax5,Rmax6,DistMat)
C     Get all 6 and 5 ring systems by checking all possible branches (vertices)
C     for each atom
C     I am sure there are better algorithms, but this one is not too bad and fast
C     enough and fast.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IC3(NAtom,3)
      DIMENSION Icon2(natom2)
      DIMENSION N5MEM(Nfaces,5),N5MEMS(Nfaces,5)
      DIMENSION N6MEM(Nfaces,6),N6MEMS(Nfaces,6)
      DIMENSION IPa(6,96)
      DIMENSION DistMat(natomL)
      DIMENSION Rd(6)
      DIMENSION Rmem(Nedges)
      Data Tol/5.d-4/
      Ndif=0
      NForbid=0
      Ncount5=1
      Ncount6=1
      Do I=1,6
      Rd(I)=0.d0
      enddo
      Do IS=1,MAtom
      Do 1 I=1,6
      Do 1 J=1,96
    1 IPa(I,J)=0
      Do 2 I=1,3
    2 IPa(1,I)=IC3(IS,I)     
      Do I2=1,3
      IX1=IPa(1,I2)
      if(IX1.ne.0) CALL Step(NAtom,2,I2,IS,IX1,IPA,IC3)
      Do I3=1,6
      IV=I3/2
      IV1=(I3+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(1,IV)
      IX2=IPa(2,I3)
      if(IX2.ne.0) CALL Step(NAtom,3,I3,IForbid,IX2,IPA,IC3)
      Do I4=1,12
      IV=I4/2
      IV1=(I4+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(2,IV)
      IX3=IPa(3,I4)
      if(IX3.ne.0) CALL Step(NAtom,4,I4,IForbid,IX3,IPA,IC3)
      Do I5=1,24
      IV=I5/2
      IV1=(I5+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(3,IV)
      IX4=IPa(4,I5)
      if(IX4.ne.0) CALL Step(NAtom,5,I5,IForbid,IX4,IPA,IC3)
      Do I6=1,48
      IV=I6/2
      IV1=(I6+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(4,IV)
      IX5=IPa(5,I6)
      if(IX5.ne.0) CALL Step(NAtom,6,I6,IForbid,IX5,IPA,IC3)
      enddo
      enddo
      enddo
      enddo
      enddo
C     Print*,' Tree ',IS
C     Do 3 I=1,6
C     Write(IOUT,1010) (IPa(I,J),J=1,48)
C   3 Write(IOUT,1011) (IPa(I,J),J=49,96)
C     Identify all 5-membered rings
      Do I=1,48
      IN5=IPa(5,I)
      IF(IN5.eq.IS) then
      CALL Ring5(Ncount5,I,IN5,natom,Nfaces,IPA,N5MEM,N5MEMS)
      endif
      enddo
C     Identify all 6-membered rings
      Do I=1,96
      IN6=IPa(6,I)
      IF(IN6.eq.IS) then
      CALL Ring6(Ncount6,I,IN6,natom,Nfaces,IPA,N6MEM,N6MEMS)
      endif
      enddo
      enddo

      Ncount5=Ncount5-1
      Write(IOUT,1000) Ncount5
C     Check bond distances
      Do I=1,Ncount5
      Rsum=0.
      Rmin5=1000.d0
      Rmax5=0.d0
      Do J=1,5
      IAT1=N5MEM(I,J)
      J1=J+1
      IF(J1.eq.6) J1=1
      IAT2=N5MEM(I,J1)
      DM=FunDistMat(IAT1,IAT2,natomL,DistMat)
      Rd(J)=DM
      if(Rd(J).LT.Rmin5) Rmin5=Rd(J)
      if(Rd(J).GT.Rmax5) Rmax5=Rd(J)
      RSum=Rsum+Rd(J)
      enddo
      RAv5=RSum*.2d0
      Rsum=0.
      Do J=1,5
      Rsum=Rsum+(RAv5-Rd(J))**2
      enddo
      Rrmsd=dsqrt(Rsum*.2d0)
      Write(IOUT,1001) (N5MEM(I,J),J=1,5),RAv5,Rrmsd,(Rd(J),J=1,5)
      CALL DifDist(Ndif,5,Tol,Rd,Rmem)
      enddo
      Write(IOUT,1007) Rmin5,Rmax5
C     Check bond angles
      asmall=1.d10
      abig=-1.d10
      Do I=1,Ncount5
      Do J=1,5
      IAT1=N5MEM(I,J)
      J1=J+1
      IF(J1.eq.6) J1=1
      IAT2=N5MEM(I,J1)
      J2=J+2
      IF(J2.eq.6) J2=1
      IF(J2.eq.7) J2=2
      IAT3=N5MEM(I,J2)
C     Deviation from ideal pentagon angle of 108 deg
      AngleM=FunAngleMat(IAT1,IAT2,IAT3,natomL,DistMat)
      if(AngleM.gt.abig) then
      abig=AngleM
      endif
      if(AngleM.lt.asmall) then
      asmall=AngleM
      endif
      enddo
      enddo
      asmalldif=asmall-1.08d2
      abigdif=abig-1.08d2
      Write(IOUT,1015) asmall,asmalldif,abig,abigdif
C
      Ncount6=Ncount6-1
      Write(IOUT,1002) Ncount6
      If(Ncount6.ne.0) Write(IOUT,1006)
C     Check bond distances
      Rmin6=1000.d0
      Rmax6=0.d0
      Do I=1,Ncount6
      Rsum=0.
      Do J=1,6
      IAT1=N6MEM(I,J)
      J1=J+1
      IF(J1.eq.7) J1=1
      IAT2=N6MEM(I,J1)
      DM=FunDistMat(IAT1,IAT2,natomL,DistMat)
      Rd(J)=DM
      if(Rd(J).LT.Rmin6) Rmin6=Rd(J)
      if(Rd(J).GT.Rmax6) Rmax6=Rd(J)
      RSum=Rsum+Rd(J)
      enddo
      RAv6=RSum/6.d0
      Rsum=0.
      Do J=1,6
      Rsum=Rsum+(RAv6-Rd(J))**2
      enddo
      Rrmsd=dsqrt(Rsum*.2d0)
      Write(IOUT,1003) (N6MEM(I,J),J=1,6),RAv6,Rrmsd,(Rd(J),J=1,6)
      CALL DifDist(Ndif,6,Tol,Rd,Rmem)
      enddo
      If(Ncount6.ne.0) Write(IOUT,1008) Rmin6,Rmax6
      If(Rmin5.lt.Rmin6) then
      NMin=5
      RminT=Rmin5
      else
      NMin=6
      RminT=Rmin6
      endif
      If(Rmax5.ge.(Rmax6-1.d-12)) then
      NMax=5
      RmaxT=Rmax5
      else
      NMax=6
      RmaxT=Rmax6
      endif
C     Check bond angles
      asmall=1.d10
      abig=-1.d10
      Do I=1,Ncount5
      Do J=1,6
      IAT1=N6MEM(I,J)
      J1=J+1
      IF(J1.eq.7) J1=1
      IAT2=N6MEM(I,J1)
      J2=J+2
      IF(J2.eq.7) J2=1
      IF(J2.eq.8) J2=2
      IAT3=N6MEM(I,J2)
C     Deviation from ideal hexagon angle of 120 deg
      AngleM=FunAngleMat(IAT1,IAT2,IAT3,natomL,DistMat)
      if(AngleM.gt.abig) then
      abig=AngleM
      endif
      if(AngleM.lt.asmall) then
      asmall=AngleM
      endif
      enddo
      enddo
      asmalldif=asmall-1.2d2
      abigdif=abig-1.2d2
      Write(IOUT,1016) asmall,asmalldif,abig,abigdif

C     Check Euler characteristic
      Ncount56=Ncount5+Ncount6
      MEuler=MAtom-Mcon2+Ncount56
      Mv=(5*Ncount5+6*Ncount6)/3
      Me=(5*Ncount5+6*Ncount6)/2
      Write(IOUT,1004) MAtom,Mcon2,Ncount56,MEuler,Ncount5,Ncount6,Mv,Me
      If(MEuler.ne.2) Write(IOUT,1005)
      If(Ncount5.ne.12) then
      Write(IOUT,1014)
      stop
      endif
      Write(IOUT,1009) NMin,RminT,NMax,RmaxT
      ameas=dfloat(Ndif)/dfloat(Mcon2)
      Write(Iout,1012) Ndif,ameas,Tol
      Write(Iout,1013) (Rmem(I),I=1,Ndif)
 1000 Format(/1X,I3,' five-membered-rings identified',/,
     1 ' Atom numbers in ring, Ni, mean distance, dm, and root mean',
     1 ' square deviation for distances, RMSD,'
     1 ' and distances in the ring:',
     1 /3X,'N1   N2   N3   N4   N5',9X,'dm',11X,'RMSD',
     1 15X,'R1',12X,'R2',12X,'R3',12X,'R4',12X,'R5')
 1001 Format(1X,5(I4,1X),3X,2(d12.6,2X),5X,5(d12.6,2X))
 1002 Format(/1X,I3,' six-membered-rings identified')
 1003 Format(1X,6(I4,1X)3X,2(d12.6,2X),5X,6(d12.6,2X))
 1004 Format(//1X,'Checking the Euler polyhedron formula:',/1X,
     1 'Number of vertices Nv: ',I4,/1X,
     1 'Number of edges Ne:    ',I4,/1X,
     1 'Number of faces Nf:    ',I4,/1X,
     1 'Euler number Nv-Ne+Nf: ',I4,
     1 ' (should be 2 for spherical polyhedra '
     1 'or planar connected graphs)',/1X,
     1 'Number of pentagons:   ',I4,/1X,
     1 'Number of hexagons:    ',I4,/1X,
     1 'Mv=',I4,1X,' Me=',I4)
 1005 Format(1X,' **** Capped Polydron does not fulfill Eulers theorem')
 1006 Format(/1X,' Atom numbers in ring,',
     1 ' Ni, mean distance, dm, and root mean',
     1 ' square deviation for distances, RMSD,'
     1 ' and distances in the ring:',
     1 /3X,'N1   N2   N3   N4   N5   N6',9X,'dm',11X,'RMSD',
     1 15X,'R1',12X,'R2',12X,'R3',12X,'R4',12X,'R5'12X,'R6')
 1007 Format(1X,' 5-ring minimum bond distance: ',d12.6,
     1 3X,' maximum bond distance: ',d12.6) 
 1008 Format(1X,' 6-ring minimum bond distance: ',d12.6,
     1 3X,' maximum bond distance: ',d12.6) 
 1009 Format(1X,' Minimum bond distance in ',I1,'-ring: ',d12.6,
     1 3X,' maximum bond distance in ',I1,'-ring: ',d12.6) 
 1010 Format(1X,96(I3))
 1011 Format(3X,96(I3))
 1012 Format(/1X,'Number of different bond distances Nr=',I4,
     1 ' and Nr/Ne= ',d12.6,
     1 ' (within tolerance ',d6.1,'),  Nonequivalent bond distances:') 
 1013 Format(10(1X,D12.6,1X))
 1014 Format(1X,'**** Severe error: 12 pentagons expected. STOP ****')
 1015 Format(1X,' 5-ring minimum bond angle (deviation from 108 deg): ',
     1 F6.2,' (',F6.2,'), maximum bond angle: ',F6.2,' (',F6.2,')')
 1016 Format(1X,' 6-ring minimum bond angle (deviation from 120 deg): ',
     1 F6.2,' (',F6.2,'), maximum bond angle: ',F6.2,' (',F6.2,')')
      RETURN
      END

      SUBROUTINE RingC(NAtom,Nfaces,Nedges,NAtom2,Matom,nat11,Iout,
     1 iprint,N5MEM,N6MEM,N5Ring,N6Ring,Nring,Iring5,Iring6,Iring56,
     1 NringA,NringB,NringC,NringD,NringE,NringF,DIST,CRing5,CRing6)
      IMPLICIT REAL*8 (A-H,O-Z)
C     Determine the center of each 5-and 6-ring system
      DIMENSION Dist(3,natom),Distac(6)
      DIMENSION CRing5(3,Nfaces),CRing6(3,Nfaces)
      DIMENSION N5MEM(Nfaces,5),N6MEM(Nfaces,6),Nring(Nfaces)
      DIMENSION IedgeA(Nedges),IedgeB(Nedges)
      DIMENSION IedgeC(Nedges),IedgeD(Nedges)
      DIMENSION IedgeE(Nedges),IedgeF(Nedges)
      DIMENSION NringA(Nedges),NringB(Nedges)
      DIMENSION NringC(Nedges),NringD(Nedges)
      DIMENSION NringE(Nedges),NringF(Nedges)
      Integer n3r(3,nat11),n3ra(3,natom),n3rb(3,natom)
      Integer IRhag5(0:5),IRhag6(0:6)
      Character*6,Label
C     Center for 5-rings
      Write(Iout,1000)
      Do I=1,N5Ring
      Nring(I)=I
      Sum5x=0.
      Sum5y=0.
      Sum5z=0.
      Do J=1,5
      IAt=N5MEM(I,J)
      Sum5x=Sum5x+Dist(1,IAt)
      Sum5y=Sum5y+Dist(2,IAt)
      Sum5z=Sum5z+Dist(3,IAt)
      enddo
      CRing5(1,I)=Sum5x*.2d0
      CRing5(2,I)=Sum5y*.2d0
      CRing5(3,I)=Sum5z*.2d0
C     Calculate distance of atoms from the center
      Do J=1,5
      IAt=N5MEM(I,J)
      X1=CRing5(1,I)-Dist(1,IAt)
      Y1=CRing5(2,I)-Dist(2,IAt)
      Z1=CRing5(3,I)-Dist(3,IAt)
      Distac(J)=dsqrt(X1*X1+Y1*Y1+Z1*Z1)
      enddo
      Write(Iout,1002)Nring(I),(N5MEM(I,JT),JT=1,5),
     1 (CRing5(JT1,I),JT1=1,3),(Distac(JT),JT=1,5)
      enddo
C     Center for 6-rings
      If(N6Ring.eq.0) Go to 2000
      Write(Iout,1001)
      Do I=1,N6Ring
      NRing(N5Ring+I)=N5Ring+I
      Sum6x=0.
      Sum6y=0.
      Sum6z=0.
      Do J=1,6
      IAt=N6MEM(I,J)
      Sum6x=Sum6x+Dist(1,IAt)
      Sum6y=Sum6y+Dist(2,IAt)
      Sum6z=Sum6z+Dist(3,IAt)
      enddo
      CRing6(1,I)=Sum6x/6.d0
      CRing6(2,I)=Sum6y/6.d0
      CRing6(3,I)=Sum6z/6.d0
C     Calculate distance of atoms from the center
      Do J=1,6
      IAt=N6MEM(I,J)
      X1=CRing6(1,I)-Dist(1,IAt)
      Y1=CRing6(2,I)-Dist(2,IAt)
      Z1=CRing6(3,I)-Dist(3,IAt)
      Distac(J)=dsqrt(X1*X1+Y1*Y1+Z1*Z1)
      enddo
      Write(Iout,1003)Nring(N5Ring+I),(N6MEM(I,JT),JT=1,6),
     1 (CRing6(JT1,I),JT1=1,3),(Distac(JT),JT=1,6)
      enddo

C     Analyzing the ring fusions
C     All 2-ring fusions
 2000 Write(Iout,1004)
      IR1=5
      IR2=5
      N2ring=0
      IRing5=0
      IRing6=0
      IRing56=0
C     (5-5) 2-ring fusions
      CALL Ring55(natom,Nfaces,Nedges,IRing5,N5Ring,NringA,NringB,Nring,
     1 N5MEM,IedgeA,IedgeB)
      Write(Iout,1005) IR1,IR2,IRing5
      if(IRing5.ne.0) Write(Iout,1006) (NringA(I),NringB(I),I=1,IRing5)
      N2ring=IRing5
      If(N6Ring.eq.0) Go to 3000

C     (5-6) 2-ring fusions
      IR2=6
      CALL Ring56(natom,Nfaces,Nedges,IRing56,N5Ring,N6Ring,
     1 NringC,NringD,Nring,N5MEM,N6MEM)
      Write(Iout,1005) IR1,IR2,IRing56
      if(IRing56.ne.0) Write(Iout,1006)(NringC(I),NringD(I),I=1,IRing56)
      N2ring=N2ring+IRing56

C     (6-6) 2-ring fusions
      IR1=6
      CALL Ring66(natom,Nfaces,Nedges,IRing6,N5Ring,N6Ring,NringE,
     1 NringF,Nring,N6MEM)
      Write(Iout,1005) IR1,IR2,IRing6
      if(IRing6.ne.0) Write(Iout,1006) (NringE(I),NringF(I),I=1,IRing6)
      N2ring=N2ring+IRing6

C     Final 2-ring
 3000 Write(Iout,1007) N2ring

C     Get Rhagavachari-Fowler-Manolopoulos neighboring pentagon and hexagon indices
C     First pentagon indices
      Do I=0,5
       IRhag5(I)=0
      enddo
      If(IRing5.eq.0) then
       IRhag5(0)=12
       go to 111
      endif
      do I=1,12
       IRcount=0
       do J=1,IRing5
        If(NRingA(J).eq.I.or.NRingB(J).eq.I) then
         IRcount=IRcount+1
        endif
       enddo
       IRhag5(IRcount)=IRhag5(IRcount)+1
      enddo
C     Pentagon index
  111 Ifus5=0
      Do I=1,5
      IFus5=IFus5+I*IRhag5(I)
      enddo
      IFus5G=IFus5/2
      Write(Iout,1013) (IRhag5(I),I=0,5),IFus5G
      If(IFus5G.eq.IRing5) then
      Write(Iout,1015) IFus5G
      else
      Write(Iout,1016) IFus5G
      endif
C     Now hexagon indices
      if(N6Ring.eq.0) go to 113
      Do I=0,6
      IRhag6(I)=0
      enddo
      ihk=0
      If(IRing6.eq.0) then
      IRhag6(0)=N6Ring
      go to 112
      endif
      do I=13,12+N6Ring
      IRcount=0
      do J=1,IRing6
      If(NRingE(J).eq.I.or.NRingF(J).eq.I) then
      IRcount=IRcount+1
      endif
      enddo
      IRhag6(IRcount)=IRhag6(IRcount)+1
      enddo
C     Strain Parameter
      khk=0
      k2hk=0
      Do I=3,6
      ihk=ihk+IRhag6(I)
      IIR=I*IRhag6(I)
      khk=khk+IIR
      k2hk=k2hk+I*IIR
      enddo
      aihk=dfloat(ihk)
      akhk2=(dfloat(khk)/aihk)**2
      ak2hk=dfloat(k2hk)/aihk
      sigmah=dsqrt(dabs(ak2hk-akhk2))
  112 Write(Iout,1020) (IRhag6(I),I=0,6),sigmah
      Ifus6=0
      Do I=3,6
      IFus6=IFus6+IRhag6(I)
      enddo
      IFus6G=IFus6*2+20
      If(IFus6G.eq.Matom) then
      Write(Iout,1018) IFus6G
      else
      Write(Iout,1019) IFus6G
      endif

C     All 3-ring fusions
  113 Write(Iout,1014) 
C     (c5-5-5) 3-ring fusions
      Label='closed'
      IR1=5
      IR2=5
      IR3=5
      N3Ring=0
      KRing3=0
      LRing3=0
      if(IRing5.gt.0) CALL Ring555(natom,nfaces,Nedges,nat11,Iout,
     1 IRing5,NringA,NringB,KRing3,LRing3,n3r,n3ra,n3rb)
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C     (o5-5-5) 3-ring fusions
      Label='open  '
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3
      If(N6Ring.eq.0) Go to 4000

C     (c5-5-6) 3-ring fusions with (5-5)
      KRing3=0
      LRing3=0
      if(IRing5.gt.0.and.IRing56.gt.0) then
      CALL Ring556(natom,nfaces,Nedges,nat11,Iout,IRing5,IRing56,
     1 NringA,NringB,NringC,NringD,KRing3,LRing3,n3r,n3ra,n3rb)
      endif
      IR3=6
      Label='closed'
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3

C     (o5-5-6) 3-ring fusions
      Label='open  '
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C     (l5-6-5) 3-ring fusions
      KRing3=0
      LRing3=0
      if(IRing56.gt.0) then
      CALL Ring565(natom,nfaces,Nedges,nat11,Iout,IRing5,IRing56,
     1 NringA,NringB,NringC,NringD,KRing3,LRing3,n3r,n3ra,n3rb,
     1 N5Ring, N6Ring,N5MEM,N6MEM)
      endif
      Label='bent  '
      IR2=6
      IR3=5
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3

C     (b5-6-5) 3-ring fusions
      Label='linear'
      Write(Iout,1009) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C     (c6-5-6) 3-ring fusions
      KRing3=0
      LRing3=0
      if(IRing56.gt.0) then
      CALL Ring656(natom,nfaces,Nedges,nat11,Iout,IRing5,IRing6,IRing56,
     1 NringC,NringD,NringE,NringF,KRing3,LRing3,n3r,n3ra,n3rb)
      endif
      Label='closed'
      IR1=6
      IR2=5
      IR3=6
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3

C     (o6-5-6) 3-ring fusions
      Label='open  '
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C     (l5-6-6) 3-ring fusions
      KRing3=0
      LRing3=0
      Label='bent  '
      IR1=5
      IR2=6
      if(IRing56.gt.0.and.IRing6.gt.0) then
      CALL Ring566(natom,nfaces,Nedges,nat11,Iout,IRing6,IRing56,
     1 NringC,NringD,NringE,NringF,N5MEM,N6MEM,KRing3,LRing3,
     1 n3r,n3ra,n3rb)
      endif
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3

C     (b5-6-6) 3-ring fusions
      Label='linear'
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C     (c6-6-6) 3-ring fusions
      KRing3=0
      LRing3l=0
      LRing3b=0
      IR1=6
      Label='closed'
      if(IRing6.gt.0) then
      CALL Ring666(natom,nfaces,Nedges,nat11,Iout,IRing6,
     1 NringE,NringF,N6MEM,KRing3,LRing3b,LRing3l,n3r,n3ra)
      endif
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3

C     (b6-6-6) 3-ring fusions
      Label='bent  '
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3b
      if(Lring3b.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3r(J,I),J=1,3),i=1,Lring3b)
      N3Ring=N3Ring+LRing3b

C     (l6-6-6) 3-ring fusions
      Label='linear'
      ndimh=nat11/2
      ncount=Lring3l-ndimh
      if(Lring3l.lt.ndimh) ncount=0
      Write(Iout,1008) Label,IR1,IR2,IR3,ncount
      if(ncount.ne.0.and.iprint.eq.1) write(Iout,1011) 
     1  ((n3r(J,I),J=1,3),i=ndimh+1,Lring3l)
      N3Ring=N3Ring+ncount

C     Final 3-ring count
 4000 N3ringexp=11*(Matom/2)-30
      Write(Iout,1010) N3ring,N3ringexp
      Ndif=N3ringexp-N3ring
      if(Ndif.ne.0) then
      Write(Iout,1012)
      endif
C     Check all connected 4-and 5-ring pentagons to create
C     pentagon indices a la Fowler and Manolopoulus

C     Similar for hexagon indices a la Fowler and Manolopoulus

 1000 Format(/1X,'Center for 5-rings',/1X,
     1 ' Ring Number RN, Atom Numbers Ni, Ring centers X,Y,Z and '
     1 'distances di from ring center to atoms',/2X,
     1 'RN     N1  N2  N3  N4  N5',9X,'X',12X,'Y',12X,'Z',
     1 12X,'d1',11X,'d2',11X,'d3',11X,'d4',11X,'d5')
 1001 Format(/1X,'Center for 6-rings',/1X,
     1 ' Ring Number RN, Atom Numbers Ni, Ring centers X,Y,Z and '
     1 'distances di from ring center to atoms',/2X,
     1 'RN     N1  N2  N3  N4  N5  N6',9X,'X',12X,'Y',12X,'Z',
     1 12X,'d1',11X,'d2',11X,'d3',11X,'d4',11X,'d5',11X,'d6')
 1002 Format(1X,I3,3X,5I4,3X,3(D12.6,1X),2X,5(D12.6,1X))
 1003 Format(1X,I3,3X,6I4,3X,3(D12.6,1X),2X,6(D12.6,1X))
 1004 Format(/1X,'Analyzing basic two- and three-ring fusions',
     1 //1X,'2-ring fusions between rings (RNI,RNJ):') 
 1005 Format(2X,'(',I1,'-',I1,') fusions: ',I5,' in total')
 1006 Format(15(1X,'(',I3,',',I3,')'))
 1007 Format(1X,'Total number of distinct two-ring fusions:',I4,
     1 ' (should be identical to the number of edges Ne)',/)
 1008 Format(2X,A6,1X,'(',I1,'-',I1,'-',I1,') fusions: ',I4)
 1009 Format(2X,A6,1X,'(',I1,'-',I1,'-',I1,') fusions: ',I4,
     1 ' (natural constraint of the second kind)')
 1010 Format(1X,'Total number of distinct three-ring fusions:',I5,
     1 ' (expected: ',I5,')')
 1011 Format(10(1X,'(',I3,',',I3,','I3,')'))
 1012 Format(' WARNING: expected 3-ring count does not match ',
     1 'number found')
 1013 Format(1X,'Rhagavachari-Fowler-Manolopoulos neighboring '
     1 'pentagon indices: (',5(I2,','),I2,')',
     1 ' and number of pentagon-pentagon fusions: ',I2)
 1014 Format(//1X,'3-ring fusions between rings (RNI,RNJ,RNK):') 
 1015 Format(1X,'Number of (5,5) fusions matches the ',I2,
     1 ' value obtained from Rhagavachar-/Fowler-Manolopoulos '
     1 'neighboring pentagon indices')
 1016 Format(1X,'Error: Number of (5,5) fusions does not match the ',I2,
     1 ' value obtained from Rhagavachari-Fowler-Manolopoulos '
     1 'neighboring pentagon indices')
 1017 Format(1X,'Rhagavachari/Fowler neighboring hexagon indices: (',
     1 6(I3,','),I3,')')
 1018 Format(1X,'Number vertices matches the ',I5,
     1 ' value obtained from Rhagavachari-Fowler-Manolopoulos '
     1 'neighboring hexagon indices ---> Fullerene is IPR')
 1019 Format(1X,'Number vertices does not match the ',I5,
     1 ' value obtained from Rhagavachari-Fowler-Manolopoulos '
     1 'neighboring hexagon indices ---> Fullerene is not IPR')
 1020 Format(1X,'Rhagavachari/Fowler neighboring hexagon indices: (',
     1 6(I3,','),I3,')  and strain parameter sigma = ',F12.6)
      Return
      END
 
      SUBROUTINE Ring666(natom,nfaces,Nedges,ndim,Iout,I6C,
     1 NrE,NrF,N6MEM,KRing3,LRing3b,LRing3l,n666,c666)
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION N6MEM(Nfaces,6),NrE(Nedges),NrF(Nedges)
      DIMENSION IS(6),JS(6)
      Integer n666(3,ndim),c666(3,natom)
C     Write out all 3-ring 566 connections first then sort them out
      i666=0
      Do I=1,I6C
      IR1=NrE(I)
      IR2=NrF(I)
      Do J=I+1,I6C
      IR3=NRE(J)
      IR4=NRF(J)
C     Check for each ring fusion which ring is connected to IR1.and.IR2
C     and store, there will be duplicates
      If(IR3.eq.IR1) then
       i666=i666+1
       if(i666.gt.ndim) stop
        if(IR2.gt.IR4) then
         n666(1,i666)=IR4 
         n666(2,i666)=IR1 
         n666(3,i666)=IR2
        else
         n666(1,i666)=IR2 
         n666(2,i666)=IR1 
         n666(3,i666)=IR4
       endif
      endif
      If(IR4.eq.IR1) then
       i666=i666+1
       if(i666.gt.ndim) stop
        if(IR2.gt.IR3) then
         n666(1,i666)=IR3 
         n666(2,i666)=IR1 
         n666(3,i666)=IR2
        else
         n666(1,i666)=IR2 
         n666(2,i666)=IR1 
         n666(3,i666)=IR3
       endif
      endif
      If(IR3.eq.IR2) then
       i666=i666+1
       if(i666.gt.ndim) stop
        if(IR1.gt.IR4) then
         n666(1,i666)=IR4 
         n666(2,i666)=IR2 
         n666(3,i666)=IR1
        else
         n666(1,i666)=IR1 
         n666(2,i666)=IR2 
         n666(3,i666)=IR4
       endif
      endif
      If(IR4.eq.IR2) then
       i666=i666+1
       if(i666.gt.ndim) stop
        if(IR1.gt.IR3) then
         n666(1,i666)=IR3 
         n666(2,i666)=IR2 
         n666(3,i666)=IR1
        else
         n666(1,i666)=IR1 
         n666(2,i666)=IR2 
         n666(3,i666)=IR3
       endif
      endif
      enddo
      enddo
      if(i666.eq.0) return

C     Determine closed structures from the list
C      and remove them from the main list
      j666=0
      Kring3=0
      Do I=1,i666
       IS(1)=n666(1,I)
       IS(2)=n666(2,I)
       IS(3)=n666(3,I)
       Do J=1,I6C
        JR1=NRE(J)
        JR2=NRF(J)
        if(IS(1).eq.JR1.and.IS(3).eq.JR2) then
           Kring3=Kring3+1
           CALL SortI(3,IS,JS)
         if(Kring3.gt.1) then
C     Sort out duplicates
          Do k=1,Kring3-1
           if(c666(1,k).eq.JS(1).and.c666(2,k).eq.JS(2).
     1      and.c666(3,k).eq.JS(3)) then
            Kring3=Kring3-1
            go to 2
           endif
          enddo
         endif
         if(Kring3.gt.natom) stop
          c666(1,Kring3)=JS(1)
          c666(2,Kring3)=JS(2)
          c666(3,Kring3)=JS(3)
          go to 2
        endif
       enddo
       j666=j666+1
       n666(1,j666)=n666(1,i)
       n666(2,j666)=n666(2,i)
       n666(3,j666)=n666(3,i)
   2  continue
      enddo
      i666=j666
      if(i666.eq.0) return

C     Now sort out linear and bent 666 structures and dividing the space
C      of aray n666 in half for bent and linear structures
      LRing3b=0
      ndimh=ndim/2
      LRing3l=ndimh
      Do I=1,i666
       IS1=n666(1,i)
       ISM=n666(2,i)
       IS3=n666(3,i)
       ISM6=ISM-12
        do j=1,6
          J1=J+1
          if(J1.eq.7) J1=1
          IP1=N6MEM(ISM6,J)  
          IP2=N6MEM(ISM6,J1)  
          itag1=0
          itag2=0
          do k=1,6
           KP1=N6MEM(IS1-12,k)
           if(KP1.eq.IP1) itag1=itag1+1
           if(KP1.eq.IP2) itag2=itag2+1
          enddo
          do k=1,6
           KP1=N6MEM(IS3-12,k)
           if(KP1.eq.IP2) itag1=itag1+1
           if(KP1.eq.IP1) itag2=itag2+1
          enddo
          if(itag1.eq.2.or.itag2.eq.2) then
           Lring3b=Lring3b+1
           if(Lring3b.gt.ndimh) stop
           n666(1,Lring3b)=IS1
           n666(2,Lring3b)=ISM
           n666(3,Lring3b)=IS3
           go to 3
          endif
        enddo
           Lring3l=Lring3l+1
           if(Lring3l.gt.ndim) stop
           n666(1,Lring3l)=IS1
           n666(2,Lring3l)=ISM
           n666(3,Lring3l)=IS3
   3  continue
      enddo 
      
C     Sort out duplicates for bent structures
      ncount=0
      Do I=1,Lring3b
       IS1=n666(1,I)
       IS2=n666(2,I)
       IS3=n666(3,I)
      Do J=I+1,Lring3b
       JS1=n666(1,J)
       JS2=n666(2,J)
       JS3=n666(3,J)
       if(IS1.eq.JS1.and.IS2.eq.JS2.and.IS3.eq.JS3) n666(1,J)=0
      enddo
      enddo
      Do I=1,Lring3b
       If(n666(1,I).ne.0) then
       ncount=ncount+1
       n666(1,ncount)=n666(1,I)
       n666(2,ncount)=n666(2,I)
       n666(3,ncount)=n666(3,I)
       endif
       enddo
      Lring3b=ncount

C     Sort out duplicates for linear structures
      ncount=ndimh
      Do I=ndimh+1,Lring3l
       IS1=n666(1,I)
       IS2=n666(2,I)
       IS3=n666(3,I)
      Do J=I+1,Lring3l
       JS1=n666(1,J)
       JS2=n666(2,J)
       JS3=n666(3,J)
       if(IS1.eq.JS1.and.IS2.eq.JS2.and.IS3.eq.JS3) n666(1,J)=0
      enddo
      enddo
      Do I=ndimh+1,Lring3l
       If(n666(1,I).ne.0) then
       ncount=ncount+1
       n666(1,ncount)=n666(1,I)
       n666(2,ncount)=n666(2,I)
       n666(3,ncount)=n666(3,I)
       endif
       enddo
      Lring3l=ncount

C     Sort both arrays n666    
C     First the bent structures
      Do I=1,Lring3b
       IS1=n666(1,I)
       IS3=n666(3,I)
      Do J=I+1,Lring3b
       JS1=n666(1,J)
       JS3=n666(3,J)
      IF(IS1-JS1) 10,11,12
  11  IF(JS3.gt.IS3) go to 10
  12  IS1=n666(1,J)
      IS2=n666(2,J)
      IS3=n666(3,J)
      n666(1,J)=n666(1,I)
      n666(2,J)=n666(2,I)
      n666(3,J)=n666(3,I)
      n666(1,I)=IS1
      n666(2,I)=IS2
      n666(3,I)=IS3
  10  continue
      enddo
      enddo

C     Next the linear structures
      Do I=ndimh+1,Lring3l
       IS1=n666(1,I)
       IS3=n666(3,I)
      Do J=I+1,Lring3l
       JS1=n666(1,J)
       JS3=n666(3,J)
      IF(IS1-JS1) 20,21,22
  21  IF(JS3.gt.IS3) go to 20
  22  IS1=n666(1,J)
      IS2=n666(2,J)
      IS3=n666(3,J)
      n666(1,J)=n666(1,I)
      n666(2,J)=n666(2,I)
      n666(3,J)=n666(3,I)
      n666(1,I)=IS1
      n666(2,I)=IS2
      n666(3,I)=IS3
  20  continue
      enddo
      enddo

      Return
      END
 
      SUBROUTINE Ring566(natom,nfaces,Nedges,ndim,Iout,I6C,I56C,
     1 NrC,NrD,NrE,NrF,N5MEM,N6MEM,KRing3,LRing3,n566,b566,l566)
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Nfaces,5),N6MEM(Nfaces,6)
      DIMENSION NrE(Nedges),NrF(Nedges)
      DIMENSION NrC(Nedges),NrD(Nedges)
      Integer n566(3,ndim),l566(3,natom),b566(3,natom)
C     Write out all 3-ring 566 connections first then sort them out
      i566=0
      Do I=1,I56C
       IR1=NrC(I)
       IR2=NrD(I)
      Do J=1,I6C
       IR3=NRE(J)
       IR4=NRF(J)
C     Check for each ring fusion which ring is connected to IR1.and.IR2
C     and store, there will be duplicates
      If(IR2.eq.IR3) then
       i566=i566+1
       n566(1,i566)=IR1 
       n566(2,i566)=IR2 
       n566(3,i566)=IR4
      endif 
      If(IR2.eq.IR4) then
       i566=i566+1
       n566(1,i566)=IR1 
       n566(2,i566)=IR2 
       n566(3,i566)=IR3
      endif
      enddo
      enddo
      if(i566.eq.0) return
      
C     Remove closed structures from the list
      j566=0
      Do I=1,i566
      IR1=n566(1,I)
      IR3=n566(3,I)
      Do J=1,I56C
       JR1=NRC(J)
       JR2=NRD(J)
       if(IR1.eq.JR1.and.IR3.eq.JR2) then
       go to 2
       endif
      enddo
      j566=j566+1
      n566(1,j566)=n566(1,i)
      n566(2,j566)=n566(2,i)
      n566(3,j566)=n566(3,i)
   2  continue
      enddo
      i566=j566
      if(i566.eq.0) return

C     Now sort out linear and bent 566 structures
      Do I=1,i566
      IS5=n566(1,i)
      ISM=n566(2,i)
      IS6=n566(3,i)
      ISM6=ISM-12
        do j=1,6
          J1=J+1
          if(J1.eq.7) J1=1
          IP1=N6MEM(ISM6,J)  
          IP2=N6MEM(ISM6,J1)  
          itag1=0
          itag2=0
          do k=1,5
           KP1=N5MEM(IS5,k)
           if(KP1.eq.IP1) itag1=itag1+1
           if(KP1.eq.IP2) itag2=itag2+1
          enddo
          do k=1,6
           KP1=N6MEM(IS6-12,k)
           if(KP1.eq.IP2) itag1=itag1+1
           if(KP1.eq.IP1) itag2=itag2+1
          enddo
          if(itag1.eq.2.or.itag2.eq.2) then
           Kring3=Kring3+1
           b566(1,Kring3)=IS5
           b566(2,Kring3)=ISM
           b566(3,Kring3)=IS6
           go to 3
          endif
        enddo
           Lring3=Lring3+1
           l566(1,Lring3)=IS5
           l566(2,Lring3)=ISM
           l566(3,Lring3)=IS6
   3  continue   
      enddo      

      Return
      END
 
      SUBROUTINE Ring656(natom,nfaces,Nedges,ndim,Iout,I5C,I6C,
     1 I56C,NrC,NrD,NrE,NrF,KRing3,LRing3,n656,c656,o656)
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrE(Nedges),NrF(Nedges)
      DIMENSION NrC(Nedges),NrD(Nedges)
      Integer n656(3,ndim),c656(3,natom),o656(3,natom)
C     Write out all 3-ring 656 connections first then sort them out
      i656=0
      Do I=1,I56C
       IR1=NrC(I)
       IR2=NrD(I)
      Do J=I+1,I56C
       IR3=NRC(J)
       IR4=NRD(J)
C     Check for each ring fusion which ring is connected to IR1.and.IR2
C     and store, there will be duplicates
       If(IR3.eq.IR1) then
       i656=i656+1
       if(IR2.gt.IR4) then
        mr=IR2
        IR2=IR4
        IR4=mr
       endif
       n656(1,i656)=IR2 
       n656(2,i656)=IR1 
       n656(3,i656)=IR4
       endif 
      enddo
      enddo
      if(i656.eq.0) return
      
C     Now sort out closed and open 656 structures
      Do 1 I=1,i656
      IS1=n656(1,i)
      ISM=n656(2,i)
      IS2=n656(3,i)
        do j=1,I6C
        IS3=NRE(J)
        IS4=NRF(J)
         if(IS1.eq.IS3.and.IS2.eq.IS4) then
          Kring3=Kring3+1
          c656(1,Kring3)=IS1
          c656(2,Kring3)=ISM
          c656(3,Kring3)=IS2
          go to 1
         endif
         enddo      
          Lring3=Lring3+1
          o656(1,Lring3)=IS1
          o656(2,Lring3)=ISM
          o656(3,Lring3)=IS2
   1  Continue

      Return
      END
 
      SUBROUTINE Ring565(natom,nfaces,Nedges,ndim,Iout,I5C,I56C,NrA,NrB,
     1 NrC,NrD,KRing3,LRing3,n565,b565,l565,N5Ring,N6Ring,N5MEM,N6MEM)
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Nfaces,5),N6MEM(Nfaces,6)
      DIMENSION NrA(Nedges),NrB(Nedges)
      DIMENSION NrC(Nedges),NrD(Nedges)
      Integer n565(3,ndim),b565(3,natom),l565(3,natom)
C     Write out all 3-ring 565 connections first then sort them out
      i565=0
      Do I=1,I56C
       IR1=NrC(I)
       IR2=NrD(I)
      Do J=I+1,I56C
       IR3=NRC(J)
       IR4=NRD(J)
C     Check for each ring fusion which ring is connected to IR1.and.IR2
C     and store, there will be duplicates
       If(IR2.eq.IR4) then
        i565=i565+1
       if(IR1.gt.IR3) then
        mr=IR1
        IR1=IR3
        IR3=mr
       endif
        n565(1,i565)=IR1 
        n565(2,i565)=IR2 
        n565(3,i565)=IR3
       endif 
      enddo
      enddo
      if(i565.eq.0) return

C     Remove the (5,5) connections
      do i=1,i565
       I1=n565(1,i)
       I2=n565(3,i)
      do j=1,I5C
       J1=NRA(j)
       J2=NRB(j)
       if(J1.eq.I1.and.J2.eq.I2) then
        n565(1,i)=0
       endif
      enddo
      enddo

C     Now look for bent fusion and add them into b556
C     If found add a 0 for n555(1,I)
      Do 1 I=1,i565
      IS1=n565(1,i)
      If(IS1.ne.0) then
       IS2=n565(2,i)
       IS3=n565(3,i)
C     Determine if 3rings are aligned linear (jring=0) or bent (jring=1)
       Do J=1,6
        i5b=0
        J1=J-1
        if(J1.eq.0) J1=6
          IS12=IS2-12
          I6a=N6MEM(IS12,J1)
          I6b=N6MEM(IS12,J)
          i5b1=0
          i5b2=0
            do k=1,5
             I5x=N5MEM(IS1,k)
             I5y=N5MEM(IS3,k)
             if(I5x.eq.I6a) i5b1=i5b1+1
             if(I5y.eq.I6b) i5b1=i5b1+1
             if(I5x.eq.I6b) i5b2=i5b2+1
             if(I5y.eq.I6a) i5b2=i5b2+1
            enddo
        if(i5b1.eq.2.or.i5b2.eq.2) go to 2
        enddo
        go to 1
   2      Kring3=Kring3+1
C         Sort in ascending order
          b565(1,Kring3)=IS1
          b565(2,Kring3)=IS2
          b565(3,Kring3)=IS3
          n565(1,i)=0
        endif
   1  Continue

C     Now store the linear fusions
      Do I=1,i565
       if(n565(1,I).ne.0) then
        Lring3=Lring3+1
        l565(1,Lring3)=n565(1,i)
        l565(2,Lring3)=n565(2,i)
        l565(3,Lring3)=n565(3,i)
       endif
      enddo

      Return
      END
 
      SUBROUTINE Ring556(natom,nfaces,Nedges,ndim,Iout,I5C,I56C,
     1 NrA,NrB,NrC,NrD,KRing3,LRing3,n556,c556,o556)
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrA(Nedges),NrB(Nedges)
      DIMENSION NrC(Nedges),NrD(Nedges)
      DIMENSION IS(6),JS(6)
      Integer n556(3,ndim),c556(3,natom),o556(3,natom)
C     Write out all 3-ring 556 connections first then sort them out
      i556=0
      Do I=1,I5C
       IR1=NrA(I)
       IR2=NrB(I)
      Do J=1,I56C
       IR3=NRC(J)
       IR4=NRD(J)
C     Check for each ring fusion which ring is connected to IR1.and.IR2
C     and store, there will be duplicates
       If(IR3.eq.IR1) then
        i556=i556+1
        n556(1,i556)=IR2 
        n556(2,i556)=IR1 
        n556(3,i556)=IR4
       endif 
       If(IR3.eq.IR2) then
        i556=i556+1
        n556(1,i556)=IR1 
        n556(2,i556)=IR2 
        n556(3,i556)=IR4
       endif 
      enddo
      enddo
      if(i556.eq.0) return

C     Determine closed 556 structures from the list
C      and remove them from the main list
      Do I=1,i556
       IS(1)=n556(1,I)
       IS(2)=n556(2,I)
       IS(3)=n556(3,I)
       Do J=1,I56C
        JR1=NRC(J)
        JR2=NRD(J)
        if(IS(1).eq.JR1.and.IS(3).eq.JR2) then
           Kring3=Kring3+1
           If(IS(1).gt.IS(2)) then
           IM=IS(1)
           IS(1)=IS(2)
           IS(2)=IM
           endif
         if(Kring3.gt.1) then
C     Sort out duplicates
          Do k=1,Kring3-1
           if(c556(1,k).eq.IS(1).and.c556(2,k).eq.IS(2).
     1      and.c556(3,k).eq.IS(3)) then
            Kring3=Kring3-1
            go to 2
           endif
          enddo
         endif
          c556(1,Kring3)=IS(1)
          c556(2,Kring3)=IS(2)
          c556(3,Kring3)=IS(3)
          go to 2
        endif
       enddo
       Lring3=Lring3+1
       o556(1,Lring3)=IS(1)
       o556(2,Lring3)=IS(2)
       o556(3,Lring3)=IS(3)
   2  continue
      enddo

      Return
      END
 
      SUBROUTINE Ring555(natom,nfaces,Nedges,ndim,Iout,I5C,
     1 NrA,NrB,KRing3,LRing3,n555,n555f,m555)
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrA(Nedges),NrB(Nedges)
      DIMENSION IS(6),JS(6)
      DIMENSION n555(3,ndim),m555(3,natom),n555f(3,natom)
C     Write out all 3-ring 555 connections first then sort them out
      i555=0
      Do I=1,I5C
      IR1=NrA(I)
      IR2=NrB(I)
      Do J=I+1,I5C
      IR3=NrA(J)
      IR4=NrB(J)
C     Check for each ring fusion which ring is connected to IR1.and.IR2
C     and store, there will be duplicates, permute first and last ting
C     if necessary
      If(IR3.eq.IR1) then
       i555=i555+1
       if(IR4.lt.IR2) then
        n555(1,i555)=IR4 
        n555(2,i555)=IR1 
        n555(3,i555)=IR2
       else
        n555(1,i555)=IR2 
        n555(2,i555)=IR1 
        n555(3,i555)=IR4
       endif 
      endif 
      If(IR4.eq.IR1) then
       i555=i555+1
       if(IR3.lt.IR2) then
        n555(1,i555)=IR3 
        n555(2,i555)=IR1 
        n555(3,i555)=IR2
       else
        n555(1,i555)=IR2 
        n555(2,i555)=IR1 
        n555(3,i555)=IR3
       endif 
      endif 
      If(IR3.eq.IR2) then
       i555=i555+1
       if(IR1.lt.IR4) then
        n555(1,i555)=IR1 
        n555(2,i555)=IR2 
        n555(3,i555)=IR4
       else
        n555(1,i555)=IR4 
        n555(2,i555)=IR2 
        n555(3,i555)=IR1
       endif 
      endif 
      If(IR4.eq.IR2) then
       i555=i555+1
       if(IR1.lt.IR3) then
        n555(1,i555)=IR1 
        n555(2,i555)=IR2 
        n555(3,i555)=IR3
       else
        n555(1,i555)=IR3 
        n555(2,i555)=IR2 
        n555(3,i555)=IR1
       endif 
      endif 
      enddo
      enddo
      if(i555.eq.0) return

C     Separate closed from open structures from the list
      Do I=1,i555
       IS(1)=n555(1,I)
       IS(2)=n555(2,I)
       IS(3)=n555(3,I)
       Do J=1,I5C
        JR1=NRA(J)
        JR2=NRB(J)
        if(IS(1).eq.JR1.and.IS(3).eq.JR2) then
           Lring3=Lring3+1
           CALL SortI(3,IS,JS)
         if(Lring3.gt.1) then
C     Sort out duplicates
          Do k=1,Lring3-1
           if(m555(1,k).eq.JS(1).and.m555(2,k).eq.JS(2).
     1      and.m555(3,k).eq.JS(3)) then
            Lring3=Lring3-1
            go to 2
           endif
          enddo
         endif
          m555(1,Lring3)=JS(1)
          m555(2,Lring3)=JS(2)
          m555(3,Lring3)=JS(3)
          go to 2
        endif
       enddo
       Kring3=Kring3+1
       n555f(1,Kring3)=IS(1)
       n555f(2,Kring3)=IS(2)
       n555f(3,Kring3)=IS(3)
   2  continue
      enddo

C     Sort array n555f
      Do I=1,Kring3
       IS1=n555f(1,I)
       IS3=n555f(3,I)
      Do J=I+1,Kring3
       JS1=n555f(1,J)
       JS3=n555f(3,J)
      IF(IS1-JS1) 10,11,12
  11  IF(JS3.gt.IS3) go to 10
  12  IS1=n555f(1,J)
      IS2=n555f(2,J)
      IS3=n555f(3,J)
      n555f(1,J)=n555f(1,I)
      n555f(2,J)=n555f(2,I)
      n555f(3,J)=n555f(3,I)
      n555f(1,I)=IS1
      n555f(2,I)=IS2
      n555f(3,I)=IS3
  10  continue
      enddo
      enddo

      Return
      END
 
      SUBROUTINE Ring56(natom,Nfaces,Nedges,IR56,N5R,N6R,NrA,NrB,Nring,
     1 N5MEM,N6MEM)
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Nfaces,5),N6MEM(Nfaces,6),Nring(Nfaces)
      DIMENSION NrA(Nedges),NrB(Nedges)
      IR56=0
C     (5-6) 2-ring fusions
      Do I=1,N5R
      Do J=1,N6R
       Do K1=1,5
       IN1=N5MEM(I,K1)
        IF(K1.EQ.5) then
        IN2=N5MEM(I,1)
        else
        IN2=N5MEM(I,K1+1)
        endif
       Do K2=1,6
        IM1=N6MEM(J,K2)
        IF(K2.EQ.6) then
        IM2=N6MEM(J,1)
        else
        IM2=N6MEM(J,K2+1)
        endif
        IF(K2.EQ.1) then
        IM3=N6MEM(J,6)
        else
        IM3=N6MEM(J,K2-1)
        endif
         IF(IN2.LT.IN1) then
          IMem=IN2
          IN2=IN1
          IN1=IMem
         endif
         IM3R=IM1
         IF(IM2.LT.IM1) then
          IMem=IM2
          IM2=IM1
          IM1=IMem
         endif
           IF(IN1.EQ.IM1.AND.IN2.EQ.IM2) then
           IR56=IR56+1
           NrA(IR56)=NRing(+I)
           NrB(IR56)=NRing(N5R+J)
           Go to 1000
           endif
         IM1=IM3R
         IF(IM3.LT.IM1) then
          IMem=IM3
          IM3=IM1
          IM1=IMem
         endif
           IF(IN1.EQ.IM1.AND.IN2.EQ.IM3) then
           IR56=IR56+1
           NrA(IR56)=NRing(+I)
           NrB(IR56)=NRing(N5R+J)
           Go to 1000
           endif
       enddo
       enddo
 1000 Continue
      enddo
      enddo
      Return
      END
 
      SUBROUTINE Ring55(natom,Nfaces,Nedges,IR5,N5R,NrA,NrB,Nring,
     1 N5MEM,IedA,IedB)
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Nfaces,5),Nring(Nfaces)
      DIMENSION NrA(Nedges),NrB(Nedges)
      DIMENSION IedA(Nedges),IedB(Nedges)
      IR5=0
C     (5-5) 2-ring fusions
      Do I=1,N5R
      Do J=I+1,N5R
       Do K1=1,5
        IN1=N5MEM(I,K1)
        IF(K1.EQ.5) then
         IN2=N5MEM(I,1)
        else
         IN2=N5MEM(I,K1+1)
        endif
       Do K2=1,5
        IM1=N5MEM(J,K2)
        IF(K2.EQ.5) then
         IM2=N5MEM(J,1)
        else
         IM2=N5MEM(J,K2+1)
        endif
        IF(K2.EQ.1) then
         IM3=N5MEM(J,5)
        else
         IM3=N5MEM(J,K2-1)
        endif
         IF(IN2.LT.IN1) then
          IMem=IN2
          IN2=IN1
          IN1=IMem
         endif
         IM3R=IM1
         IF(IM2.LT.IM1) then
          IMem=IM2
          IM2=IM1
          IM1=IMem
         endif
           IF(IN1.EQ.IM1.AND.IN2.EQ.IM2) then
           IR5=IR5+1
           NrA(IR5)=NRing(I)
           NrB(IR5)=NRing(J)
           IedA(IR5)=IN1
           IedB(IR5)=IN2
           Go to 1000
           endif
         IM1=IM3R
         IF(IM3.LT.IM1) then
          IMem=IM3
          IM3=IM1
          IM1=IMem
         endif
           IF(IN1.EQ.IM1.AND.IN2.EQ.IM3) then
           IR5=IR5+1
           NrA(IR5)=NRing(I)
           NrB(IR5)=NRing(J)
           IedA(IR5)=IN1
           IedB(IR5)=IN2
           Go to 1000
           endif
       enddo
       enddo
 1000 Continue
      enddo
      enddo
      Return
      END
 
      SUBROUTINE Ring66(natom,nfaces,Nedges,IR6,N5R,N6R,NrA,NrB,
     1 Nring,N6MEM)
      IMPLICIT INTEGER (A-Z)
      DIMENSION N6MEM(Nfaces,6),Nring(Nfaces)
      DIMENSION NrA(Nedges),NrB(Nedges)
      IR6=0
C     (6-6) 2-ring fusions
      Do I=1,N6R
      Do J=I+1,N6R
       Do K1=1,6
       IN1=N6MEM(I,K1)
        IF(K1.EQ.6) then
        IN2=N6MEM(I,1)
        else
        IN2=N6MEM(I,K1+1)
        endif
       Do K2=1,6
        IM1=N6MEM(J,K2)
        IF(K2.EQ.6) then
        IM2=N6MEM(J,1)
        else
        IM2=N6MEM(J,K2+1)
        endif
        IF(K2.EQ.1) then
        IM3=N6MEM(J,6)
        else
        IM3=N6MEM(J,K2-1)
        endif
         IF(IN2.LT.IN1) then
          IMem=IN2
          IN2=IN1
          IN1=IMem
         endif
         IM3R=IM1
         IF(IM2.LT.IM1) then
          IMem=IM2
          IM2=IM1
          IM1=IMem
         endif
           IF(IN1.EQ.IM1.AND.IN2.EQ.IM2) then
           IR6=IR6+1
           NrA(IR6)=NRing(N5R+I)
           NrB(IR6)=NRing(N5R+J)
           Go to 1000
           endif
         IM1=IM3R
         IF(IM3.LT.IM1) then
          IMem=IM3
          IM3=IM1
          IM1=IMem
         endif
           IF(IN1.EQ.IM1.AND.IN2.EQ.IM3) then
           IR6=IR6+1
           NrA(IR6)=NRing(N5R+I)
           NrB(IR6)=NRing(N5R+J)
           Go to 1000
           endif
       enddo
       enddo
 1000 Continue
      enddo
      enddo
      Return
      END
 
      SUBROUTINE Distmatrix(NAtom,natomL,MAtom,IOUT,Iprint,Iopt,
     1 Dist,DistMat,Rmin,Rmax,Vol,ASphere)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,natom),DistMat(natomL)
C     Calculate distance matrix between atoms from cartesian coordinates
      DATA API/3.14159265358979d0/
      IMat=0
      if(Iopt.eq.0.and.Iprint.eq.1) Write(IOUT,1000) 
      Do I=2,MAtom
      Do J=1,I-1
      Sum=0.d0
      Do K=1,3
      Sum=Sum+(Dist(K,I)-Dist(K,J))**2
      enddo
      IMat=IMat+1
      DistMat(IMat)=dsqrt(Sum)
      enddo
      enddo
      Do I=2,MAtom
      IZ=((I-1)*(I-2))/2
      if(Iopt.eq.0.and.Iprint.eq.1) 
     1   Write(IOUT,1001) (I,J,DistMat(IZ+J),J=1,I-1)
      enddo
C     Determine minimum and maximum distances
      Rmin=1.d20
      Rmax=1.d0
      RmaxS=1.d0
      Do I=1,MAtom
      Sum=0.d0
      Do K=1,3
      Sum=Sum+Dist(K,I)**2
      enddo
      Sumd=dsqrt(Sum)
      If(Sumd.gt.RmaxS) RmaxS=Sumd
      Do J=I+1,MAtom
      DM=FunDistMat(I,J,natomL,DistMat)
      If(DM.lt.Rmin) then
      Rmin=DM
      Imin=I
      Jmin=J
      endif
      If(DM.gt.Rmax) then
      Rmax=DM
      Imax=I
      Jmax=J
      endif
      enddo
      enddo
      Vol=4.d0/3.d0*API*RmaxS**3
      ASphere=4.d0*API*RmaxS**2
      Write(IOUT,1002) Rmin,Imin,Jmin,Rmax,Imax,Jmax,RmaxS,Vol,
     1 ASphere,ASphere/Vol
C     Volume of ideal capped icosahedron
      If(Matom.eq.60) then
      Redge=Rmin*3.d0
      dsr5=dsqrt(5.d0)
C     Calculate square cos36 and sin36
      coss36=.125d0*(3.d0+dsr5)
      sins36=1.d0-coss36
      fac1=dsqrt(4.d0-1.d0/sins36)
      Rico=Redge/fac1
      VolIco=(5.d0/12.d0)*(3.d0+dsr5)*Redge**3
      VolIcocap=0.25d0*(125.d0+43.d0*dsr5)*(Redge/3.d0)**3
      Write(IOUT,1003) Rico,Redge,VolIco,VolIcocap 
      endif 
C     Calculate the volume for C50 using R5
      If(Matom.eq.20) then
      dsr5=dsqrt(5.d0)
      dsr3=dsqrt(3.d0)
      fac1=.25d0*(15.d0+7.d0*dsr5)
      Rdode=.25d0*dsr3*(1.d0+dsr5)*Rmin
      VolIcocap=fac1*Rmin**3
      Write(IOUT,1004) Rdode,Rmin,VolIcocap 
      endif
 1000 Format(/1X,'Distance Matrix:')
 1001 Format(5(1X,'('I3,',',I3,')',1X,D15.8))
 1002 Format(/1X,'Minimum distance ',D15.8,
     1 ' between atoms ',I3,' and ',I3,
     1      /1X,'Maximum distance ',D15.8,
     1 ' between atoms ',I3,' and ',I3,
     1 /1X,' Radius of covering central sphere ',D15.8,
     1 /1X,' Volume of covering central sphere (upper limit)',
     1 D15.8,' in units cube of distances'
     1 /1X,' Area of covering central sphere (upper limit)  ',D15.8,
     1 ' in units square of distances',
     1 /1X,' Area to volume ratio of covering central sphere (3/R) ',
     1 D15.8,' in units of inverse distance')
 1003 Format(1X,' Radius of (uncapped) icosahedron (obtained from '
     1 'minimum distance): ',D15.8,' and edge length: ',D15.8,/,
     1 '  Volume of ideal icosahedron for C60:',D15.8,
     1 ' and of ideal capped icosahedron:',
     1 D15.8,' (in units cube of distances)')
 1004 Format(1X,' Radius of dodecahedron: ',D15.8,
     1 ' and edge length: ',D15.8,/,
     1 '  Volume of ideal dodecahedron for C20:',D15.8)
      RETURN
      END

      SUBROUTINE Ring5(Ncount5,IN,I5,natom,nfaces,IPA,N5MEM,N5MEMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPa(6,96),N5MEM(Nfaces,5),N5MEMS(Nfaces,5)
      DIMENSION ISM(6),JSM(6)
C     Identify all 5-membered rings
C     Print*,'5-membered ring detected',IN
      ISM(5)=I5
      IT4=IN/2 
      IT4x=(IN+1)/2
      IDIF=IT4x-IT4
      IF(IDIF.EQ.1) IT4=IT4+1
      ISM(4)=IPA(4,IT4)
      IT3=IT4/2 
      IT3x=(IT4+1)/2
      IDIF=IT3x-IT3
      IF(IDIF.EQ.1) IT3=IT3+1
      ISM(3)=IPA(3,IT3)
      IT2=IT3/2 
      IT2x=(IT3+1)/2
      IDIF=IT2x-IT2
      IF(IDIF.EQ.1) IT2=IT2+1
      ISM(2)=IPA(2,IT2)
      IT1=IT2/2 
      IT1x=(IT2+1)/2
      IDIF=IT1x-IT1
      IF(IDIF.EQ.1) IT1=IT1+1
      ISM(1)=IPA(1,IT1)
      IF(Ncount5.eq.1) then
      Do I=1,5
      N5MEM(1,I)=ISM(I)
      enddo
      NCount5=NCount5+1
      CALL SortI(5,ISM,JSM)
      Do I=1,5
      N5MEMS(1,I)=JSM(I)
      enddo
C      Print*,ISM,' Sorted ',JSM
      return
      endif
      CALL SortI(5,ISM,JSM)
      IREP=0
      Do I=1,NCount5-1
      IF(N5MEMS(I,1).EQ.JSM(1).AND.N5MEMS(I,2).EQ.JSM(2).AND.
     1  N5MEMS(I,3).EQ.JSM(3).AND.N5MEMS(I,4).EQ.JSM(4).AND.
     1  N5MEMS(I,5).EQ.JSM(5)) IREP=1
      enddo
      IF(IREP.eq.0) then
      Do I=1,5
      N5MEM(NCount5,I)=ISM(I)
      N5MEMS(NCount5,I)=JSM(I)
      enddo
      NCount5=NCount5+1
      endif
C     Print*,ISM,' Sorted ',JSM
      RETURN
      END

      SUBROUTINE Ring6(Ncount6,IN,I6,natom,nfaces,IPA,N6MEM,N6MEMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPa(6,96),N6MEM(Nfaces,6),N6MEMS(Nfaces,6)
      DIMENSION ISM(6),JSM(6)
C     Identify all 6-membered rings
C     Print*,'6-membered ring detected',IN
      ISM(6)=I6
      IT5=IN/2 
      IT5x=(IN+1)/2
      IDIF=IT5x-IT5
      IF(IDIF.EQ.1) IT5=IT5+1
      ISM(5)=IPA(5,IT5)
      IT4=IT5/2 
      IT4x=(IT5+1)/2
      IDIF=IT4x-IT4
      IF(IDIF.EQ.1) IT4=IT4+1
      ISM(4)=IPA(4,IT4)
      IT3=IT4/2 
      IT3x=(IT4+1)/2
      IDIF=IT3x-IT3
      IF(IDIF.EQ.1) IT3=IT3+1
      ISM(3)=IPA(3,IT3)
      IT2=IT3/2 
      IT2x=(IT3+1)/2
      IDIF=IT2x-IT2
      IF(IDIF.EQ.1) IT2=IT2+1
      ISM(2)=IPA(2,IT2)
      IT1=IT2/2 
      IT1x=(IT2+1)/2
      IDIF=IT1x-IT1
      IF(IDIF.EQ.1) IT1=IT1+1
      ISM(1)=IPA(1,IT1)
      IF(Ncount6.eq.1) then
      Do I=1,6
      N6MEM(1,I)=ISM(I)
      enddo
      NCount6=NCount6+1
      CALL SortI(6,ISM,JSM)
      Do I=1,6
      N6MEMS(1,I)=JSM(I)
      enddo
C     Print*,ISM,' Sorted ',JSM
      return
      endif
      CALL SortI(6,ISM,JSM)
      IREP=0
      Do I=1,NCount6-1
      IF(N6MEMS(I,1).EQ.JSM(1).AND.N6MEMS(I,2).EQ.JSM(2).AND.
     1  N6MEMS(I,3).EQ.JSM(3).AND.N6MEMS(I,4).EQ.JSM(4).AND.
     1  N6MEMS(I,5).EQ.JSM(5).AND.N6MEMS(I,6).EQ.JSM(6)) IREP=1
      enddo
      Do I=1,NCount6-1
      IF(JSM(1).Eq.JSM(2).OR.JSM(2).Eq.JSM(3).OR.JSM(3).Eq.JSM(4)
     1 .OR.JSM(4).Eq.JSM(5).OR.JSM(5).Eq.JSM(6)) IREP=1
      enddo
      IF(IREP.eq.0) then
      Do I=1,6
      N6MEM(NCount6,I)=ISM(I)
      N6MEMS(NCount6,I)=JSM(I)
      enddo
      NCount6=NCount6+1
      endif
C     Print*,ISM,' Sorted ',JSM
      RETURN
      END

      SUBROUTINE Step(Natom,ITree,ILoop,IForbid,INum,IPa,IC3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IC3(Natom,3)
      DIMENSION IPa(6,96)
      DIMENSION IM(2)
C     Filling in Tree structure at level Istep
      I=1
      I1=IC3(INum,1)
      IF(I1.NE.IForbid) then
      IM(1)=I1
      I=2
      endif
      I2=IC3(INum,2)
      IF(I2.NE.IForbid) then
      IM(I)=I2
      endif
      I3=IC3(INum,3)
      IF(I3.NE.IForbid) then
      IM(2)=I3
      endif
      I=2*ILoop-1
      IPA(ITree,I)=IM(1)
      IPA(ITree,I+1)=IM(2)
      RETURN
      END

      SUBROUTINE Connect(NAtom,natomL,Natom2,MCon2,MAtom,
     1 Ipent,IOUT,Icon2,IC3,IDA,Tol,DistMat,Rmin)
C     Get the connectivities between 2 and 3 atoms
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DistMat(natomL)
      DIMENSION Icon2(natom2),IDA(NAtom,NAtom)
      DIMENSION NCI(12),NCJ(12)
      DIMENSION IC3(Natom,3)
      Rtol=Rmin*(1.d0+Tol)
      Mcon2=0
      Do I=1,NAtom
      Do J=1,3
      IC3(I,J)=0
      enddo
      enddo
      if(Ipent.eq.0) then
      Do I=1,MAtom
      Do J=I+1,MAtom
      DM=FunDistMat(I,J,natomL,DistMat)
      If (DM.lt.Rtol) then
      Mcon2=Mcon2+1
      Ncount=I*MAtom+J
      Icon2(Mcon2)=Ncount
      endif
      enddo
      enddo
      else
      Do I=1,MAtom
      Do J=I+1,MAtom
      If (IDA(I,J).eq.1) then
      Mcon2=Mcon2+1
      Ncount=I*MAtom+J
      Icon2(Mcon2)=Ncount
      endif
      enddo
      enddo
      endif
      Write(IOUT,1000) Mcon2
      If(Mcon2.lt.MAtom) then
      Write(IOUT,1004)
      Stop
      endif
      Do I=1,Mcon2,12
      Do J=1,12
      NCI(J)=0
      NCJ(J)=0
      enddo
      M12=12
      Do J=1,12
      IArray=I+J-1
      CALL Num2(MAtom,Icon2(IArray),NCI(J),NCJ(J))
      If(NCI(J).Eq.0) then
      M12=J-1
      Go to 11
      endif
      enddo
   11 Write(IOUT,1001) (NCI(J),NCJ(J),J=1,M12)
      enddo
      Write(IOUT,1002)
C     Get all vertices
      Do I=1,MAtom
      IZ=0
      Do J=1,MCon2
      CALL Num2(MAtom,Icon2(J),IX,IY)
      IF(IX.EQ.I) then
      IZ=IZ+1
      IC3(I,IZ)=IY
      endif
      IF(IY.EQ.I) then
      IZ=IZ+1
      IC3(I,IZ)=IX
      endif
      IF(IZ.EQ.3) Go to 10
      enddo
   10 Continue
      Write(IOUT,1003) I,(IC3(I,J),J=1,3)
      enddo
C     Check if structure is alright at this point
      nexpedge=MAtom*3/2
      if(Mcon2.ne.nexpedge) then
      Write(IOUT,1005) Mcon2,nexpedge
      stop
      endif
 1000 Format(/1X,' Number of connected surface atoms (edges, bonds):',
     1 I4,/1X,' Connectivities (edge set):')
 1001 Format(1X,12('{',I3,',',I3,'} '))
 1002 Format(1X,' Calculate all vertices N and corresponding ',
     1 'adjacencies Ni of 3-connected graph:',
     1 /1X,'   N       N1  N2  N3')
 1003 Format(1X,I4,'    (',3I4,')')
 1004 Format(1X,'**** Error, not enough connected atoms',
     1 ' check coordinates')
 1005 Format(1X,'**** Severe error, number of edges (bonds) not as ',
     1 'expected from number of atoms: ',I4,' (expected: ',I4,')')
      RETURN
      END
