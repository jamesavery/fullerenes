      SUBROUTINE Ring(Me,MCon2,IOUT,Ncount5,Ncount6,
     1 IC3,IVR3,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,DistMat)
      use config
C     Get all 6 and 5 ring systems by checking all possible branches (vertices)
C     for each atom
C     I am sure there are better algorithms, but this one is not too bad and fast
C     enough and fast.
      IMPLICIT REAL*8 (A-H,O-Z)
C     IC3 contains vertex adjacencies and IVR3 ring numbers for a vertex
      DIMENSION IC3(Nmax,3),IVR3(Nmax+4,3)
      DIMENSION N5MEM(Mmax,5),N5MEMS(Mmax,5)
      DIMENSION N6MEM(Mmax,6),N6MEMS(Mmax,6)
      DIMENSION IPa(6,96)
      DIMENSION DistMat(NmaxL)
      DIMENSION Rd(6)
      DIMENSION Rmem(Emax)
      Data Tol/5.d-4/
      Ndif=0
      NForbid=0
      Ncount5=1
      Ncount6=1
      Do I=1,6
      Rd(I)=0.d0
      enddo
      Do IS=1,number_vertices
      Do 1 I=1,6
      Do 1 J=1,96
    1 IPa(I,J)=0
      Do 2 I=1,3
    2 IPa(1,I)=IC3(IS,I)     
      Do I2=1,3
      IX1=IPa(1,I2)
      if(IX1.ne.0) CALL Step(2,I2,IS,IX1,IPA,IC3)
      Do I3=1,6
      IV=I3/2
      IV1=(I3+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(1,IV)
      IX2=IPa(2,I3)
      if(IX2.ne.0) CALL Step(3,I3,IForbid,IX2,IPA,IC3)
      Do I4=1,12
      IV=I4/2
      IV1=(I4+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(2,IV)
      IX3=IPa(3,I4)
      if(IX3.ne.0) CALL Step(4,I4,IForbid,IX3,IPA,IC3)
      Do I5=1,24
      IV=I5/2
      IV1=(I5+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(3,IV)
      IX4=IPa(4,I5)
      if(IX4.ne.0) CALL Step(5,I5,IForbid,IX4,IPA,IC3)
      Do I6=1,48
      IV=I6/2
      IV1=(I6+1)/2
      IDIF=IV1-IV
      IF(IDIF.EQ.1) IV=IV+1
      IForbid=IPa(4,IV)
      IX5=IPa(5,I6)
      if(IX5.ne.0) CALL Step(6,I6,IForbid,IX5,IPA,IC3)
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
      CALL Ring5(Ncount5,I,IN5,Mmax,IPA,N5MEM,N5MEMS)
      endif
      enddo
C     Identify all 6-membered rings
      Do I=1,96
      IN6=IPa(6,I)
      IF(IN6.eq.IS) then
      CALL Ring6(Ncount6,I,IN6,Mmax,IPA,N6MEM,N6MEMS)
      endif
      enddo
      enddo

      Ncount5=Ncount5-1
      Write(IOUT,1000) Ncount5
C     Check bond distances
      Rmin5=1000.d0
      Rmax5=0.d0
      Do I=1,Ncount5
       Rsum=0.
       Do J=1,5
        IAT1=N5MEM(I,J)
        J1=J+1
        IF(J1.eq.6) J1=1
        IAT2=N5MEM(I,J1)
        Rd(J)=FunDistMat(IAT1,IAT2,DistMat)
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
      AngleM=FunAngleMat(IAT1,IAT2,IAT3,DistMat)
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
      Rd(J)=FunDistMat(IAT1,IAT2,DistMat)
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
      MaxN=5
      RmaxT=Rmax5
      else
      MaxN=6
      RmaxT=Rmax6
      endif

C     Check bond angles
      asmall=1.d10
      abig=-1.d10
      Do I=1,Ncount6
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
      AngleM=FunAngleMat(IAT1,IAT2,IAT3,DistMat)
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
      if(Ncount6.ne.0)
     1  Write(IOUT,1016) asmall,asmalldif,abig,abigdif

C     Establish ring numbers for specific vertex
      do I=1,nmax+4
        Do J=1,3
          IVR3(I,J)=0
        enddo
      enddo
C     First pentagons
      Do I=1,Ncount5
      Do J=1,5
       IAT=N5MEM(I,J)
       if(IVR3(IAT,1).eq.0) then
        IVR3(IAT,1)=I
       else
        if(IVR3(IAT,2).eq.0) then
         IVR3(IAT,2)=I
        else
         IVR3(IAT,3)=I
        endif
       endif
      enddo
      enddo
C     Next hexagons
      Do I=1,Ncount6
       II=I+12
      Do J=1,6
       IAT=N6MEM(I,J)
       if(IVR3(IAT,1).eq.0) then
        IVR3(IAT,1)=II
       else
        if(IVR3(IAT,2).eq.0) then
         IVR3(IAT,2)=II
        else
         IVR3(IAT,3)=II
        endif
       endif
      enddo
      enddo
      Write(IOUT,1017)
      Do I=1,number_vertices,5
        Write(IOUT,1018) I,(IVR3(I,J),J=1,3),
     1   I+1,(IVR3(I+1,J),J=1,3),I+2,(IVR3(I+2,J),J=1,3),
     1   I+3,(IVR3(I+3,J),J=1,3),I+4,(IVR3(I+4,J),J=1,3)
      enddo
      
C     Check Euler characteristic
      Ncount56=Ncount5+Ncount6
      MEuler=number_vertices-Mcon2+Ncount56
      Mv=(5*Ncount5+6*Ncount6)/3
      Me=(5*Ncount5+6*Ncount6)/2
      Write(IOUT,1004) number_vertices,Mcon2,Ncount56,MEuler,Ncount5,
     1  Ncount6,Mv,Me
      If(MEuler.ne.2) Write(IOUT,1005)
      If(Ncount5.ne.12) then
      Write(IOUT,1014)
      stop 20
      endif
      Write(IOUT,1009) NMin,RminT,MaxN,RmaxT
      ameas=dfloat(Ndif)/dfloat(Mcon2)
      Write(Iout,1012) Ndif,ameas,Tol
      Write(Iout,1013) (Rmem(I),I=1,Ndif)

      call flush(iout)

 1000 Format(/1X,I3,' five-membered-rings identified',/,
     1 ' Atom numbers in ring, Ni, mean distance, dm, and root mean',
     1 ' square deviation for distances, RMSD,'
     1 ' and distances in the ring:',
     1 /4X,'N1    N2    N3    N4    N5',9X,'dm',11X,'RMSD',
     1 15X,'R1',12X,'R2',12X,'R3',12X,'R4',12X,'R5')
 1001 Format(1X,5(I5,1X),3X,2(d12.6,2X),5X,5(d12.6,2X))
 1002 Format(/1X,I6,' six-membered-rings identified')
 1003 Format(1X,6I6,3X,2(d12.6,2X),5X,6(d12.6,2X))
 1004 Format(//1X,'Checking the Euler polyhedron formula:',/1X,
     1 'Number of vertices Nv: ',I6,/1X,
     1 'Number of edges Ne:    ',I6,/1X,
     1 'Number of faces Nf:    ',I6,/1X,
     1 'Euler number Nv-Ne+Nf: ',I6,
     1 ' (should be 2 for spherical polyhedra '
     1 'or planar connected graphs)',/1X,
     1 'Number of pentagons:   ',I6,/1X,
     1 'Number of hexagons:    ',I6,/1X,
     1 'Mv=',I6,1X,' Me=',I6)
 1005 Format(1X,' **** Capped Polydron does not fulfill Eulers theorem')
 1006 Format(/1X,' Atom numbers in ring,',
     1 ' Ni, mean distance, dm, and root mean',
     1 ' square deviation for distances, RMSD,'
     1 ' and distances in the ring:',
     1 /5X,'N1    N2    N3    N4    N5    N6',8X,'dm',11X,'RMSD',
     1 15X,'R1',12X,'R2',12X,'R3',12X,'R4',12X,'R5'12X,'R6')
 1007 Format(1X,' 5-ring minimum bond distance: ',d12.6,
     1 3X,' maximum bond distance: ',d12.6) 
 1008 Format(1X,' 6-ring minimum bond distance: ',d12.6,
     1 3X,' maximum bond distance: ',d12.6) 
 1009 Format(1X,' Minimum bond distance in ',I1,'-ring: ',d12.6,
     1 3X,' maximum bond distance in ',I1,'-ring: ',d12.6) 
C1010 Format(1X,96(I3))
C1011 Format(3X,96(I3))
 1012 Format(/1X,'Number of different bond distances Nr=',I4,
     1 ' and Nr/Ne= ',d12.6,
     1 ' (within tolerance ',d6.1,'),  Nonequivalent bond distances:') 
 1013 Format(10(1X,D12.6,1X))
 1014 Format(1X,'**** Severe error: 12 pentagons expected. STOP ****')
 1015 Format(1X,' 5-ring minimum bond angle (deviation from 108 deg): ',
     1 F6.2,' (',F8.2,'), maximum bond angle: ',F6.2,' (',F8.2,')')
 1016 Format(1X,' 6-ring minimum bond angle (deviation from 120 deg): ',
     1 F6.2,' (',F8.2,'), maximum bond angle: ',F6.2,' (',F8.2,')')
 1017 Format(/1X,'List of ring numbers RNj containing vertex Ni  ',
     1 '(Ni: RNj, RNk, RNl)',/1X,135('-'))
 1018 Format(5(1X,'(',I5,':',I5,',',I5,',',I5,') '))
      RETURN
      END

      SUBROUTINE RingC(Medges,Iout,iprint,
     1 N5Ring,N6Ring,Nring,Iring5,Iring6,Iring56,
     1 nl565,numbersw,numberFM,numberYF,numberBF,
     1 N5MEM,N6MEM,NringA,NringB,NringC,NringD,NringE,NringF,
     1 IC3,IVR3,n3rc,nSW,nFM,nYF,nBF,SmallRingDist,DIST,CRing5,CRing6)
C     This routine analyzes the pentagons and hexagons
C     The first 12 faces are pentagons (many routines need this order)
C     All other which follow are hexagons
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Determine the center of each 5-and 6-ring system
      DIMENSION Dist(3,Nmax),Distac(6)
      DIMENSION CRing5(3,Mmax),CRing6(3,Mmax)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6),Nring(Mmax)
      DIMENSION IC3(Nmax,3),IVR3(nmax+4,3) ! up to four values past the required ones are read
      DIMENSION IedgeA(Emax),IedgeB(Emax)
      DIMENSION NringA(Emax),NringB(Emax)
      DIMENSION NringC(Emax),NringD(Emax)
      DIMENSION NringE(Emax),NringF(Emax)
      Integer n3r(3,Nmax*11),n3ra(3,Nmax),n3rb(3,Nmax),nSW(4,66),
     1 nFM(4,66),n3rc(3,Nmax),n3rd(3,Nmax),nYF(6,66),nBF(5,66)
      Integer MPatt(30)
      Character*6,Label
      ic566=0
      ib566=0
      il566=0
      ic666=0
      ib666=0
      il666=0
C     Center for 5-rings
      numbersw=0
      numberFM=0
      numberBF=0
      numberYF=0
      sigmah=0.d0
      Write(Iout,1000)
      Do I=1,N5Ring
      Nring(I)=I
      Sum5x=0.d0
      Sum5y=0.d0
      Sum5z=0.d0
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
      Sum6x=0.d0
      Sum6y=0.d0
      Sum6z=0.d0
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

C     Get the largest ring to ring distance
 2000 Rmin5=1.d10
      Rmin6=1.d10
      Rmin56=1.d10
      Rmax5=0.d0
      Rmax6=0.d0
      Rmax56=0.d0
      Do I=1,N5Ring
      Do J=I+1,N5Ring
       X=CRing5(1,I)-CRing5(1,J)
       Y=CRing5(2,I)-CRing5(2,J)
       Z=CRing5(3,I)-CRing5(3,J)
       R=dsqrt(X*X+Y*Y+Z*Z)
       if(R.gt.Rmax5) Rmax5=R
       if(R.lt.Rmin5) Rmin5=R
      enddo
      enddo
      if(N6Ring.eq.0) then
       Write(Iout,1025) RMin5,Rmax5
       go to 2001
      endif
      Do I=1,N6Ring
      Do J=I+1,N6Ring
       X=CRing6(1,I)-CRing6(1,J)
       Y=CRing6(2,I)-CRing6(2,J)
       Z=CRing6(3,I)-CRing6(3,J)
       R=dsqrt(X*X+Y*Y+Z*Z)
       if(R.gt.Rmax6) Rmax6=R
       if(R.lt.Rmin6) Rmin6=R
      enddo
      enddo
      Do I=1,N5Ring
      Do J=1,N6Ring
       X=CRing5(1,I)-CRing6(1,J)
       Y=CRing5(2,I)-CRing6(2,J)
       Z=CRing5(3,I)-CRing6(3,J)
       R=dsqrt(X*X+Y*Y+Z*Z)
       if(R.gt.Rmax56) Rmax56=R
       if(R.lt.Rmin56) Rmin56=R
      enddo
      enddo
      Write(Iout,1026) Rmin5,Rmin6,Rmin56,Rmax5,Rmax6,Rmax56
C     Smallest ring distance
 2001 SmallRingDist=Rmin6
      if(Rmin56.lt.Rmin6) SmallRingDist=Rmin56
      if(Rmin5.lt.SmallRingDist) SmallRingDist=Rmin5

C     Analyzing the ring fusions
C     All 2-ring fusions
      Write(Iout,1004)
      IR1=5
      IR2=5
      N2ring=0
      IRing5=0
      IRing6=0
      IRing56=0
C     (5-5) 2-ring fusions
      CALL Ring55(IRing5,N5Ring,NringA,NringB,Nring,
     1 N5MEM,IedgeA,IedgeB)
      Write(Iout,1005) IR1,IR2,IRing5
      if(IRing5.ne.0) Write(Iout,1006) (NringA(I),NringB(I),I=1,IRing5)
      N2ring=IRing5
      If(N6Ring.eq.0) Go to 3000

C     (5-6) 2-ring fusions
      IR2=6
      CALL Ring56(IRing56,N5Ring,N6Ring,
     1 NringC,NringD,Nring,N5MEM,N6MEM)
      Write(Iout,1005) IR1,IR2,IRing56
      if(IRing56.ne.0) Write(Iout,1006)(NringC(I),NringD(I),I=1,IRing56)
      N2ring=N2ring+IRing56

C     (6-6) 2-ring fusions
      IR1=6
      CALL Ring66(IRing6,N5Ring,N6Ring,NringE,
     1 NringF,Nring,N6MEM)
      Write(Iout,1005) IR1,IR2,IRing6
      if(IRing6.ne.0) Write(Iout,1006) (NringE(I),NringF(I),I=1,IRing6)
      N2ring=N2ring+IRing6

C     Final 2-ring
 3000 Write(Iout,1007) N2ring

C Analysis using pentagon and hexagon indices
      Call PentHexIndex(Iout,IPR,IFus5G,IRing5,IRing6,N6Ring,
     1 NRingA,NRingB,NRingE,NRingF)

C     All 3-ring fusions
      Write(Iout,1014) 
C     (c5-5-5) 3-ring fusions
      Label='closed'
      IR1=5
      IR2=5
      IR3=5
      N3Ring=0
      KRing3=0
      LRing3=0
      if(IRing5.gt.0) CALL Ring555(IRing5,NringA,
     1 NringB,KRing3,LRing3,n3r,n3ra,n3rb)
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

C     (c5-5-6) 3-ring fusions with (5-5)
      If(N6Ring.eq.0) Go to 4000
      KRing3=0
      LRing3=0
      if(IRing5.gt.0.and.IRing56.gt.0) then
      CALL Ring556(IRing5,IRing56,NringA,NringB,NringC,NringD,
     1 KRing3,LRing3,n3r,n3ra,n3rb)
      endif
      IR3=6
      Label='closed'
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      ic556=KRing3
      if(Kring3.ne.0) then
       if(iprint.eq.1) write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
C       Search for Brinkmann-Fowler D2h 55-6-55 patterns
        Call SixvertexinsertWS(Kring3,n3ra,numberBF,nBF)
        N3Ring=N3Ring+KRing3
       endif

C     (o5-5-6) 3-ring fusions
      Label='open  '
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C     (b5-6-5) 3-ring fusions
      KRing3=0
      LRing3=0
      if(IRing56.gt.0) then
      CALL Ring565(IRing5,IRing56,NringA,NringB,NringC,NringD,
     1 KRing3,LRing3,n3r,n3ra,n3rb,N5MEM,N6MEM)
      endif
      Label='bent  '
      IR2=6
      IR3=5
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3
C     Search for Stone-Wales pattern
      Call StoneWales(Kring3,n3ra,numbersw,nSW)
C     Search for Yoshida-Fowler D3h 6555 patterns - C60-like corner patch
      Call Threevertexinsert(Kring3,n3ra,numberFM,nFM)
C     (l5-6-5) 3-ring fusions
      Label='linear'
      Write(Iout,1009) Label,IR1,IR2,IR3,LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3
C     Store linear 5-6-5 for Cioslowski's scheme and Endo-Kroto insertion
       Do J1=1,3
       Do J2=1,Lring3
        n3rc(J1,J2)=n3rb(J1,J2)
       enddo
       enddo
       nl565=Lring3

C     (c6-5-6) 3-ring fusions
      KRing3=0
      LRing3=0
      if(IRing56.gt.0) then
      CALL Ring656(IRing6,IRing56,NringC,NringD,NringE,NringF,
     1 KRing3,LRing3,n3r,n3ra,n3rb)
      endif
      Label='closed'
      IR1=6
      IR2=5
      IR3=6
      ic566=KRing3
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3
C     Store closed 6-5-6 for Yashido-Fowler 9-vertex insertion
       Do J1=1,3
       Do J2=1,Kring3
        n3rd(J1,J2)=n3ra(J1,J2)
       enddo
       enddo
       nc656=Lring3

C     (o6-5-6) 3-ring fusions
      Label='open  '
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3
      io566=LRing3
      if(Lring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3rb(J,I),J=1,3),i=1,Lring3)
      N3Ring=N3Ring+LRing3

C Cioslowski's increment scheme for IPR fullerenes
      if(IPR.eq.1) Call Cioslowski(Kring3,nl565,IRing56,IRing6,
     1 N6Ring,NringC,NringD,NringE,NringF,
     1 n3ra,n3rc,Mpatt,E1,E2)
C     (l5-6-6) 3-ring fusions
      KRing3=0
      LRing3=0
      Label='bent  '
      IR1=5
      IR2=6
      if(IRing56.gt.0.and.IRing6.gt.0) then
      CALL Ring566(IRing6,IRing56,NringC,NringD,NringE,NringF,
     1 N5MEM,N6MEM,KRing3,LRing3,n3r,n3ra,n3rb)
      endif
      il566=KRing3
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3

C     (b5-6-6) 3-ring fusions
      Label='linear'
      ib566=LRing3
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
      CALL Ring666(IRing6,NringE,NringF,
     1 N6MEM,KRing3,LRing3b,LRing3l,n3r,n3ra)
      endif
      ic666=KRing3
      Write(Iout,1008) Label,IR1,IR2,IR3,KRing3
      if(Kring3.ne.0.and.iprint.eq.1) 
     1 write(Iout,1011) ((n3ra(J,I),J=1,3),i=1,Kring3)
      N3Ring=N3Ring+KRing3
C     Search for Yoshida-Fowler D3h 666555 patterns - C80-like corner patch
      Call Ninevertexinsert(Kring3,n3ra,nc656,n3rd,numberYF,nYF)

C     (b6-6-6) 3-ring fusions
      Label='bent  '
      Write(Iout,1008) Label,IR1,IR2,IR3,LRing3b
      ib666=LRing3b
      if(Lring3b.ne.0.and.iprint.eq.1.and.number_vertices.lt.1000) 
     1 write(Iout,1011) ((n3r(J,I),J=1,3),i=1,Lring3b)
      if(Lring3b.ne.0.and.iprint.eq.1.and.number_vertices.ge.1000) 
     1 write(Iout,1021) ((n3r(J,I),J=1,3),i=1,Lring3b)
      N3Ring=N3Ring+LRing3b

C     (l6-6-6) 3-ring fusions
      Label='linear'
      ndimh=(Nmax*11)/2
      ncount=Lring3l-ndimh
      if(Lring3l.lt.ndimh) ncount=0
      il666=ncount
      Write(Iout,1008) Label,IR1,IR2,IR3,ncount
      if(ncount.ne.0.and.iprint.eq.1) write(Iout,1011) 
     1  ((n3r(J,I),J=1,3),i=ndimh+1,Lring3l)
      N3Ring=N3Ring+ncount

C     Final 3-ring count
 4000 N3ringexp=11*(number_vertices/2)-30
      Write(Iout,1010) N3ring,N3ringexp
      Ndif=N3ringexp-N3ring
      if(Ndif.ne.0) then
      Write(Iout,1012)
      endif

C     Get the hexagon signatures from the paper of
C     Y. Ju, H. Liang, J. Zhang, and F. Bai, MATCH Commun. 
C     Math. Comput. Chem. 64, 419424 (2010) and
C     D. Stevanovic, MATCH Commun. Math. Comput. Chem. 66, 285292 (2011).
      ih0=number_vertices/2-10
      ih1=3*(number_vertices-20)-iring56
      ih2=18*(number_vertices-20)-(6*iring56+2*ic566+il566+ib566)
      Write(Iout,1015) ih0,ih1,ih2
      if(IPR.eq.1) then
       ih1=2*iring6
       ih2=9*number_vertices+2*(ib666+il666)-480
       Write(Iout,1013) ih0,ih1,ih2
      endif

C Print Stone-Wales patterns
      if(numbersw.eq.0.and.N6Ring.ne.0) then 
       Write(Iout,1030)
      else
       Write(Iout,1031) numbersw
       Write(Iout,1032) ((nsw(I,J),I=1,4),J=1,numbersw)
      endif
C Reti-Laszlo Analysis using Stone-Wales patterns
      Call ArmsIndices(Iout,IFus5G,numbersw,nsw)

C Print Yoshida-Fowler D3h 6555 patterns
      if(numberFM.eq.0) then 
       Write(Iout,1040)
      else
       Write(Iout,1041) numberFM
       Write(Iout,1042) ((nFM(I,J),I=1,4),J=1,numberFM)
      endif

C Print Yoshida-Fowler D3h 6555 patterns
      if(numberYF.eq.0) then 
       Write(Iout,1043)
      else
       Write(Iout,1044) numberYF
       Write(Iout,1045) ((nYF(I,J),I=1,6),J=1,numberYF)
      endif

C Print Brinkmann-Fowler D2h 55-6-55 patterns
      if(numberBF.eq.0) then 
       Write(Iout,1046)
      else
       Write(Iout,1047) numberBF
       Write(Iout,1048) ((nBF(I,J),I=1,5),J=1,numberBF)
      endif

C Print Alcami's heat of formation from structural patterns
       Write(Iout,1024) Medges
       Call Alcami(Iout,Medges,IC3,IVR3)

C Print Cioslowsky analysis and check of correctness
      if(IPR.eq.1) then 
       Write(Iout,1033) N6Ring
       Write(Iout,1034) (Mpatt(I),I=1,30)
       Call CheckCioslowski(Iout,number_vertices,Mpatt)
        Esum=E1+E2
        EC60=Esum/dfloat(number_vertices)-6.1681359415162888d2/6.d1
       Write(Iout,1035) E1,E2,Esum
       if(number_vertices.lt.100)
     1   Write(Iout,1036) number_vertices,number_vertices,EC60
       if(number_vertices.ge.100.and.number_vertices.lt.1000) 
     1   Write(Iout,1037) number_vertices,number_vertices,EC60
       if(number_vertices.ge.1000)
     1   Write(Iout,1038) number_vertices,number_vertices,EC60
       graphenel=30.336d0*.5d0-6.1681359415162888d2/6.d1
       Write(Iout,1039) graphenel
      endif

 1000 Format(/1X,'Center for 5-rings',/1X,
     1 ' Ring Number RN, Atom Numbers Ni, Ring centers X,Y,Z and '
     1 'distances di from ring center to atoms',/2X,
     1 'RN       N1    N2    N3    N4    N5',9X,'X',12X,'Y',12X,'Z',
     1 12X,'d1',11X,'d2',11X,'d3',11X,'d4',11X,'d5')
 1001 Format(/1X,'Center for 6-rings',/1X,
     1 ' Ring Number RN, Atom Numbers Ni, Ring centers X,Y,Z and '
     1 'distances di from ring center to atoms',/4X,
     1 'RN       N1    N2    N3    N4    N5    N6',
     1  9X,'X',12X,'Y',12X,'Z',
     1 12X,'d1',11X,'d2',11X,'d3',11X,'d4',11X,'d5',11X,'d6')
 1002 Format(I4,3X,5I6,3X,3(D12.6,1X),2X,5(D12.6,1X))
 1003 Format(I6,3X,6I6,3X,3(D12.6,1X),2X,6(D12.6,1X))
 1004 Format(/1X,'Analyzing basic two- and three-ring fusions',
     1 //1X,'2-ring fusions between rings (RNI,RNJ):') 
 1005 Format(2X,'(',I1,'-',I1,') fusions: ',I5,' in total')
 1006 Format(12(1X,'(',I5,',',I5,')'))
 1007 Format(1X,'Total number of distinct two-ring fusions:',I5,
     1 ' (should be identical to the number of edges Ne)',/)
 1008 Format(2X,A6,1X,'(',I1,'-',I1,'-',I1,') fusions: ',I5)
 1009 Format(2X,A6,1X,'(',I1,'-',I1,'-',I1,') fusions: ',I5,
     1 ' (natural constraint of the second kind, Endo-Kroto pattern)')
 1010 Format(1X,'Total number of distinct three-ring fusions:',I6,
     1 ' (expected: ',I6,')')
 1011 Format(10(1X,'(',I3,',',I3,','I3,')'))
 1012 Format(' WARNING: expected 3-ring count does not match ',
     1 'number found')
 1013 Format(1X,'nth moment hexagon signatures for IPR-fullerenes ',
     1 'from Ju et al: H0 = ',I5,', H1 = ',I5,', H2 = ',I5)
 1014 Format(//1X,'3-ring fusions between rings (RNI,RNJ,RNK):') 
 1015 Format(1X,'nth moment hexagon signatures from Stevanovic:',
     1 ' H0 = ',I5,', H1 = ',I5,', H2 = ',I5)
 1021 Format(10(1X,'(',I3,',',I3,','I3,')'))
 1024 Format(/1X,'Calculate Standard Enthalpy for fullerene ',
     1 'from structural motifs (M)',/2X,'M.Alcami, G.Sanchez, ',
     1 'S.Diaz-Tendero, Y.Wang, F.Martin, J. Nanosci. Nanotech. ',
     1 '7, 1329 (2007)',
     1 /1X,'Loop through all ',I5,' edges.')
 1025 Format(/1X,'Pentagon to pentagon distance: smallest ',D12.6,
     1 ', largest ',D12.6)
 1026 Format(/1X,'Smallest distances between faces: ',
     1 D12.6,' (5-5), ',D12.6,' (6-6), 'D12.6,' (5-6)',
     1 /1X,'Largest  distances between faces: ',
     1 D12.6,' (5-5), ',D12.6,' (6-6), 'D12.6,' (5-6)')
 1030 Format(/1X,'No Stone-Wales patterns found')
 1031 Format(/1X,I2,' Stone-Wales patterns found:')
 1032 Format(7(' (',I2,',',I5,',',I5,',',I2,')'))
 1033 Format(/1X,'Calculate Standard Enthalpy for IPR fullerene ',
     1 'from structural motifs (M)',/2X,'J.Cioslowski, ',
     1 'N.Rao, D.Moncrieff, J. Am. Chem. Soc. 122, 8265 (2000)',
     1 /1X,'Loop through all ',I5,' hexagons. Ring patterns:')
 1034 Format  (' M666666: ',I6,', M666665: ',I4,', M666655: ',I4,
     1        ', M666565: ',I4,', M665665: ',I4,', M666555: ',I4,
     1       /,' M665655: ',I6,', M656565: ',I4,', M665555: ',I4,
     1        ', M656555: ',I4,', M655655: ',I4,', M655555: ',I4,
     1       /,' M555555: ',I6,', M6666  : ',I4,', M6665  : ',I4,
     1        ', M6656  : ',I4,', M6655  : ',I4,', M6565  : ',I4,
     1       /,' M6556  : ',I6,', M5665  : ',I4,', M6555  : ',I4,
     1        ', M5655  : ',I4,', M5555  : ',I4,', M13/66 : ',I4,
     1       /,' M13/56 : ',I6,', M13/55 : ',I4,', M14/66 : ',I4,
     1        ', M14/56 : ',I4,', M14/55 : ',I4,', M135   : ',I4)
 1035 Format(1X,'Enthalpy H of formation:',/1X,
     1 'Motif term:              ',F12.3,' kcal/mol',/1X,
     1 'Curvature term:          ',F12.3,' kcal/mol',/1X,
     1 'Total enthalpy:          ',F12.3,' kcal/mol')
 1036 Format(' H(C60)/60 - H(C',I2,')/',I2,4X,F12.3,' kcal/mol')
 1037 Format(' H(C60)/60 - H(C',I3,')/',I3,2X,F12.3,' kcal/mol')
 1038 Format(' H(C60)/60 - H(C',I5,')/',I5,F12.3,' kcal/mol')
 1039 Format('       (graphene limit: ',F12.3,' kcal/mol)')
 1040 Format(/1X,'No Yoshida-Fowler D3h 6555 pattern (C60-like corner',
     1 ' patch) found')
 1041 Format(/1X,I2,' Yoshida-Fowler D3h 6555 patterns (C60-like ', 
     1 'corner  patch) found:')
 1042 Format(8(' (',I5,',',I2,',',I2,',',I2,')'))
 1043 Format(/1X,'No Yoshida-Fowler D3h 666555 pattern (C80-like ',
     1 'corner patch) found')
 1044 Format(/1X,I2,' Yoshida-Fowler D3h 666555 patterns (C80-like ', 
     1 'corner  patch) found:')
 1045 Format(4(' (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,') '))
 1046 Format(/1X,'No B-F D2h 55-6-55 pattern found')
 1047 Format(/1X,I2,' B-F D2h 55-6-55 patterns found:')
 1048 Format(5(' (',I2,',',I2,',',I5,',',I2,',',I2,') '))
      Return
      END

      Subroutine PentHexIndex(Iout,IPR,IFus5G,IRing5,IRing6,N6Ring,
     1 NRingA,NRingB,NRingE,NRingF)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NringA(Emax),NringB(Emax)
      DIMENSION NringE(Emax),NringF(Emax)
      Integer IRhag5(0:5),IRhag6(0:6)
C Get Rhagavachari-Fowler-Manolopoulos neighboring pentagon and hexagon indices
C     First pentagon indices
      IPR=0
      ihk=0
      sigmah = 0
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
  111 IFus5G=IPentInd(IRhag5)
      Write(Iout,1013) (IRhag5(I),I=0,5),IFus5G
      If(IFus5G.eq.IRing5) then
       Write(Iout,1015) IFus5G
      else
       Write(Iout,1016) IFus5G
      endif
C     Now hexagon indices
      if(N6Ring.eq.0) Return
      Do I=0,6
        IRhag6(I)=0
      enddo
      If(IRing6.ne.0) then
        do I=13,12+N6Ring
          IRcount=0
          do J=1,IRing6
            If(NRingE(J).eq.I.or.NRingF(J).eq.I) then
              IRcount=IRcount+1
            endif
          enddo
          IRhag6(IRcount)=IRhag6(IRcount)+1
        enddo
C       Hexagon Neighbor Index
        ih0=0
        ih1=0
        ih2=0
        Do I=0,6
          ih0=ih0+IRhag6(I)
          ih1=ih1+I*IRhag6(I)
          ih2=ih2+I*I*IRhag6(I)
        enddo
C       Strain Parameter
        sigmah=HexInd(IRhag6,ihk)
        Write(Iout,1024) ih0,ih1,ih2 
        if(ihk.eq.0) Write(Iout,1027) 
      else
        IRhag6(0)=N6Ring
      endif
      Write(Iout,1020) (IRhag6(I),I=0,6),sigmah
      Ifus6=0
      Do I=3,6
        IFus6=IFus6+IRhag6(I)
      enddo
      IFus6G=IFus6*2+20
      If(IFus6G.eq.number_vertices) then
        Write(Iout,1018) IFus6G
      else
        Write(Iout,1019) IFus6G
      endif
      If(IRing5.eq.0) then
        IPR=1
        Write(Iout,1022) 
      else
        Write(Iout,1023)
      endif

 1013 Format(1X,'Rhagavachari-Fowler-Manolopoulos neighboring '
     1 'pentagon indices: (',5(I2,','),I2,') and number of fused ',
     1 'pentagon pairs (pentagon signature) Np: ',I2)
 1015 Format(1X,'Number of (5,5) fusions matches the ',I2,
     1 ' value obtained from Rhagavachari-/Fowler-Manolopoulos '
     1 'neighboring pentagon indices')
 1016 Format(1X,'Error: Number of (5,5) fusions does not match the ',I2,
     1 ' value obtained from Rhagavachari-Fowler-Manolopoulos '
     1 'neighboring pentagon indices')
 1018 Format(1X,'Number vertices matches the ',I5,
     1 ' value obtained from Rhagavachari-Fowler-Manolopoulos '
     1 'neighboring hexagon indices h3+h4+h5+h6')
 1019 Format(1X,'Number vertices does not match the ',I5,
     1 ' value obtained from Rhagavachari-Fowler-Manolopoulos '
     1 'neighboring hexagon indices h3+h4+h5+h6')
 1020 Format(1X,'Rhagavachari/Fowler neighboring hexagon indices: (',
     1 6(I3,','),I5,')  and strain parameter sigma = ',F12.6)
 1022 Format(1X,'--> Fullerene is IPR')
 1023 Format(1X,'--> Fullerene is not IPR')
 1024 Format(1X,'nth moment hexagon signatures from neighboring ',
     1 'hexagon indices: H0 = ',I5,', H1 = ',I5,', H2 = ',I5)
 1027 Format(1X,'sum hk is zero -> sigmah set to zero')
      Return
      END

      Subroutine ArmsIndices(Iout,Np,numbersw,nsw)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer nSW(4,66),IArm(0:5),IA(12)
C Calculate Reti-Laszlo Index
      M1=0
      M2=0
      Do I=0,5
       IArm(I)=0
      enddo
      Do I=1,12
       IA(I)=0
      enddo
      if(numbersw.eq.0) then
       IArm(0)=12
      else
       do I=1,numbersw
        I1=nsw(1,I)
        I2=nsw(4,I)
        IA(I1)=IA(I1)+1
        IA(I2)=IA(I2)+1
       enddo
        do I=1,12
         nf=IA(I)
         IArm(nf)=IArm(nf)+1
        enddo
      endif
      ibal=0
      isum=0
      do I=0,5
       isum=isum+IArm(I)
       if(IArm(I).eq.12) ibal=I
      enddo
      Write(Iout,1000) (IArm(I),I=0,5),isum
C Calculate pentagon arm indices
C Still to do
      Do I=1,5
       ic=I*IArm(I)
       M1=M1+ic
       M2=M2+I*ic
      enddo
      NA=M1/12
      AM1=dfloat(M1)/12.d0
      AM2=dfloat(M2)/12.d0
      VAR=AM2-AM1*AM1
      c=dsqrt(1.2d2*AM2/(1+7.d0*AM1))/(1.d0+.9d0*(AM2-AM1**2)**.2d0)
      psi=(3.d1+6.d0*AM1)/(1.d0+4.5d0*dfloat(Np)+c)
      Write(Iout,1001) AM1,AM2,NA,Np+NA
      Write(Iout,1002) ibal
      if(NP+NA.eq.0) Write(Iout,1003) 

 1000 Format(/1X,'Reti-Laszlo topological analysis using Stone-Wales ',
     1 'patterns:',/1X,'Pentagon arm indices: ',
     1 '(',5(I2,','),I2,')',' SUM= ',I2)
 1001 Format(1X,'M1= ',F8.2,', M2= ',F8.2,', NA= ',I2,', NP+NA= ',I2)
 1002 Format(1X,'Fullerene is ',I1,'-balanced')
 1003 Format(1X,'This is a strongly isolated fullerene')
      Return
      END

      SUBROUTINE EdgeCoord(Iout,Dist,IC3)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax)
      DIMENSION IC3(Nmax,3)
C     Print center of edges
      Write(Iout,1000)
      Do I=1,number_vertices
      Do J=1,3
       IAtom=IC3(I,J)
       if(IAtom.gt.I) then
        X=(Dist(1,I)+Dist(1,IAtom))*.5d0
        Y=(Dist(2,I)+Dist(2,IAtom))*.5d0
        Z=(Dist(3,I)+Dist(3,IAtom))*.5d0
        Write(Iout,1001) I,IAtom,X,Y,Z
       endif
      enddo
      enddo

 1000 Format(/1X,'Print center coordinates of edges:',
     1 /1X,'    I    J       X            Y            Z',
     1 /1X,49('-')) 
 1001 Format(1X,2I5,3(1X,F12.8))
      Return
      END

      SUBROUTINE RingCoord(Iout,dualdist,R6,
     1 SmallRingDist,Dist,N5,N6,N5M,N6M)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax),N5M(Mmax,5),N6M(Mmax,6)
C     Print center of rings
      factor=1.d0
      if(dualdist.ne.R6) then
       factor=dualdist/SmallRingDist
       Write(Iout,1002) factor,dualdist
      endif
      Write(Iout,1000)
      IR=5
      Do I=1,N5
        X=0.
        Y=0.
        Z=0.
      Do J=1,5
        X=X+Dist(1,N5M(I,J))
        Y=Y+Dist(2,N5M(I,J))
        Z=Z+Dist(3,N5M(I,J))
      enddo
        X=X*factor/5.
        Y=Y*factor/5.
        Z=Z*factor/5.
        Write(Iout,1001) I,IR,X,Y,Z
      enddo

      IR=6
      Do I=1,N6
        X=0.
        Y=0.
        Z=0.
      Do J=1,6
        X=X+Dist(1,N6M(I,J))
        Y=Y+Dist(2,N6M(I,J))
        Z=Z+Dist(3,N6M(I,J))
      enddo
        X=X*factor/6.
        Y=Y*factor/6.
        Z=Z*factor/6.
        Write(Iout,1001) I,IR,X,Y,Z
      enddo

 1000 Format(/1X,'Print center coordinates of faces:',
     1 /1X,'    I    IR      X            Y            Z',
     1 /1X,49('-')) 
 1001 Format(1X,2I5,3(1X,F12.8))
 1002 Format(1X,'Coordinates multiplied by ',F12.8,
     1 ' to reach distance n dual of ',F12.8)
      Return
      END

      SUBROUTINE Alcami(Iout,Medges,IC3,IVR3)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      integer IC3(Nmax,3),IVR3(number_vertices+4,3)
      integer nring(4),npattern(9)
      real*8 eps(9)
C     Finds different ring patterns of 4 rings connected
C     see M.Alcami, G.Sanchez, S.Diaz-Tendero, Y.Wang, F.Martin, 
C      J. Nanosci. Nanotech. 7, 1329 (2007)
C     Parameter from Table IV in the paper
      Data eps/19.8d0,17.6d0,10.3d0,15.7d0,12.4d0,7.8d0,
     1 6.2d0,4.7d0,1.7d0/

      do I=1,9
       npattern(i)=0
      enddo

C     Go through all edges (IV1,IV2)
      Do IV1=1,number_vertices
      Do J=1,3
       IV2=IC3(IV1,J)
       if(IV1.lt.IV2) then
C     Analyze the rings
C     Find rings left-right to edge
        do nr=1,4
         nring(nr)=0
        enddo
         nrfind=0
        do k1=1,3
        do k2=1,3
         If(IVR3(IV1,k1).eq.IVR3(IV2,k2)) then
          nrfind=nrfind+1
          nring(nrfind)=IVR3(IV1,k1)
         endif
        enddo
        enddo
        do k1=1,3
         if(nring(1).ne.IVR3(IV1,k1).and.nring(2).ne.IVR3(IV1,k1)) 
     1   nring(3)=IVR3(IV1,k1)
        enddo
        do k2=1,3
         if(nring(1).ne.IVR3(IV2,k2).and.nring(2).ne.IVR3(IV2,k2)) 
     1   nring(4)=IVR3(IV2,k2)
        enddo

C      Distinguish between 9 cases of patterns
       if(nring(1).le.12.and.nring(2).le.12) then

C      55 case        
C       55-55
        if(nring(3).le.12.and.nring(4).le.12) then
         npattern(1)=npattern(1)+1
        else
C       55-66
         if(nring(3).gt.12.and.nring(4).gt.12) then
          npattern(3)=npattern(3)+1
         else
C       55-56
          npattern(2)=npattern(2)+1
         endif
        endif

       else

        if(nring(1).gt.12.and.nring(2).gt.12) then

C      66 case
C       66-55
        if(nring(3).le.12.and.nring(4).le.12) then
         npattern(7)=npattern(7)+1
        else
C       66-66
         if(nring(3).gt.12.and.nring(4).gt.12) then
          npattern(9)=npattern(9)+1
         else
C       66-56
          npattern(8)=npattern(8)+1
         endif
        endif

        else

C      56 or 65 case
C       56-55
        if(nring(3).le.12.and.nring(4).le.12) then
         npattern(4)=npattern(4)+1
        else
C       56-66
         if(nring(3).gt.12.and.nring(4).gt.12) then
          npattern(6)=npattern(6)+1
         else
C       56-56
          npattern(5)=npattern(5)+1
         endif
        endif

       endif        
       endif        

      endif
      enddo
      enddo

      Write(Iout,1000) (npattern(i),i=1,9)

C     Calculate the enthalpy of formation per bond
      energy=0.d0
      Do I=1,9
       energy=energy+eps(i)*dfloat(npattern(i))
      enddo
      Ebond=energy/dfloat(Medges)
      Write(Iout,1001) energy,Ebond
       EC60=energy/dfloat(number_vertices)-6.54d2/6.d1
       if(number_vertices.lt.100)
     1   Write(Iout,1002) number_vertices,number_vertices,EC60
       if(number_vertices.ge.100.and.number_vertices.lt.1000)
     1   Write(Iout,1003) number_vertices,number_vertices,EC60
       if(number_vertices.ge.1000)
     1   Write(Iout,1004) number_vertices,number_vertices,EC60
       graphenel=1.7d0*.5d0-6.1681359415162888d2/6.d1
       Write(Iout,1005) graphenel

 1000 Format(1x,'Structural motifs for four connected compact rings:',
     1 /1x,' pp-pp: ',I5,', pp-hp: ',I5,', pp-hh: ',I5,
     1 /1x,' ph-pp: ',I5,', ph-hp: ',I5,', ph-hh: ',I5,
     1 /1x,' hh-pp: ',I5,', hh-hp: ',I5,', hh-hh: ',I5)
 1001 Format(1x,'Enthalpy of formation: ',F12.3,' kcal/mol',
     1 ' (per bond: ',F12.3,' kcal/mol)')
 1002 Format(' H(C60)/60 - H(C',I2,')/',I2,4X,F12.3,' kcal/mol')
 1003 Format(' H(C60)/60 - H(C',I3,')/',I3,2X,F12.3,' kcal/mol')
 1004 Format(' H(C60)/60 - H(C',I5,')/',I5,F12.3,' kcal/mol')
 1005 Format('       (graphene limit: ',F12.3,' kcal/mol)')
      Return
      END

      SUBROUTINE SixvertexinsertWS(Kring3,n3ra,numberBF,nBF)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION n3ra(3,Nmax),nBF(5,66)
C Find Brinkmann-Fowler D2h 55-6-55 patterns
      numberBF=0
      do I=1,5
      do J=1,66
       nBF(I,J)=0
      enddo
      enddo
      ntrans=0
      do I=1,Kring3
      do J=I+1,Kring3
       if(n3ra(3,I).eq.n3ra(3,J)) then
        if(n3ra(1,I).ne.n3ra(1,J).and.n3ra(2,I).ne.n3ra(1,J).and.
     1   n3ra(1,I).ne.n3ra(2,J).and.n3ra(2,I).ne.n3ra(2,J)) then
         ntrans=ntrans+1 
         nBF(3,ntrans)=n3ra(3,I)
         nBF(1,ntrans)=n3ra(1,I)
         nBF(2,ntrans)=n3ra(2,I)
         nBF(4,ntrans)=n3ra(1,J)
         nBF(5,ntrans)=n3ra(2,J)
        endif
       endif
      enddo
      enddo
      numberBF=ntrans
      
      Return
      END

      SUBROUTINE Ninevertexinsert(Kring3,n3ra,nc656,n3rd,numberYF,nYF)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION n3ra(3,Nmax),nFM(4,66)
      DIMENSION n3rd(3,Nmax),nYF(6,66),mem(3)
C Find Yoshida-Fowler D3h 666555 pattern for 9-vertex insertion
      numberYF=0
      if(Kring3.eq.0.or.nc656.eq.0) Return
      do I=1,6
      do J=1,66
       nfm(I,J)=0
      enddo
      enddo
      do I=1,Kring3
       I1=n3ra(1,I)
       I2=n3ra(2,I)
       I3=n3ra(3,I)
       icount=0
       do J=1,nc656
        J1=n3rd(1,J)
        J2=n3rd(3,J)
        if(I1.eq.J1.and.I2.eq.J2) then
         icount=icount+1
         mem(icount)=n3rd(2,J)
        endif
        if(I2.eq.J1.and.I3.eq.J2) then
         icount=icount+1
         mem(icount)=n3rd(2,J)
        endif
        if(I1.eq.J1.and.I3.eq.J2) then
         icount=icount+1
         mem(icount)=n3rd(2,J)
        endif
       enddo
      if(icount.eq.3) then
       numberYF=numberYF+1
       nYF(1,numberYF)=I1
       nYF(2,numberYF)=I2
       nYF(3,numberYF)=I3
       nYF(4,numberYF)=mem(1)
       nYF(5,numberYF)=mem(2)
       nYF(6,numberYF)=mem(3)
      endif
      enddo
      Return
      END

      SUBROUTINE Threevertexinsert(Kring3,n3ra,numberFM,nFM)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION n3ra(3,Nmax),nFM(4,66)
      DIMENSION IC(4)
C Find Yoshida-Fowler D3h 6555 pattern for 3-vertex insertion
      numberfm=0
      if(Kring3.eq.0) Return
      icount=0
      mring=0
      do I=1,4
      do J=1,66
       nfm(I,J)=0
      enddo
      enddo
      do I=1,Kring3
      do J=I+1,Kring3
       if(n3ra(2,I).eq.n3ra(2,J).and.n3ra(2,I).ne.mring) then
        icount=icount+1
        nfm(1,icount)=n3ra(2,I)
        mring=n3ra(2,I)
        IC(1)=n3ra(1,I)
        IC(2)=n3ra(3,I)
        IC(3)=n3ra(1,J)
        IC(4)=n3ra(3,J)
        Do K1=1,4
        Do K2=K1+1,4
         if(IC(K2).lt.IC(K1)) then
          II=IC(K1)
          IC(K1)=IC(K2)
          IC(K2)=II
         endif
        enddo
        enddo
        Do K=2,4
         if(IC(K).eq.IC(K-1)) then
          Do K1=K+1,4
           IC(K1-1)=IC(K1)
          enddo
         endif
        enddo
        nfm(2,icount)=IC(1)
        nfm(3,icount)=IC(2)
        nfm(4,icount)=IC(3)
        go to 9
       endif
      enddo
   9  continue
      enddo

C Sort and check if any duplicates
      do I=1,icount
      do J=I+1,icount
       if(nfm(1,J).lt.nfm(1,I)) then
        I1=nfm(1,I)
        I2=nfm(2,I)
        I3=nfm(3,I)
        I4=nfm(4,I)
        nfm(1,I)=nfm(1,J)
        nfm(2,I)=nfm(2,J)
        nfm(3,I)=nfm(3,J)
        nfm(4,I)=nfm(4,J)
        nfm(1,J)=I1
        nfm(2,J)=I2
        nfm(3,J)=I3
        nfm(4,J)=I4
       endif
      enddo
      enddo

      numberfm=icount
      do I=2,icount
       if(nfm(1,I).eq.nfm(1,I-1)) then
        do J=I+1,icount
         do k=1,4
          nfm(k,J-1)=nfm(k,J)
         enddo
        enddo
        numberfm=numberfm-1
       endif
      enddo
 
      Return
      END

      SUBROUTINE CheckCioslowski(Iout,number_vertices,p)
      IMPLICIT INTEGER (A-Z)
      DIMENSION p(30)
      nhex=number_vertices/2-10
      sum=0
C Check if number of patches is correct
      do I=1,30
       sum=sum+p(I)
      enddo
      if(sum.ne.nhex) Write(Iout,1000) sum,nhex
C Check if eqs.(2a-e) from Cioslowski's paper are fulfilled

      c2a=2*(p(2)+p(3)+p(6)+p(9)+p(12)-p(14))+
     1 4*(p(4)+p(5)+p(7)+p(10)+p(11))+6*p(8)-p(15)+p(18)+p(22)

C There is a typo in Cioslowski's paper concerning N655555. 
C Equation (2b) should be
C    N666655 + 2 N666555 +   N665655 + 3 N665555 + 2 N656555 + 
C  2 N655655 + 4 N655555 + 6 N555555 +   N6556    +  N6555   + 
C    N5555   -   N13/66  = 0

      c2b=2*(p(6)+p(10)+p(11))+3*p(9)+6*p(13)+p(3)+p(7)+4*p(12)
     1 +p(19)+p(21)+p(23)-p(24)

C In Cioslowski's paper 13/65 should be 13/56, same with 14/65
      c2c=2*p(23)+p(17)+p(21)+p(22)-p(25)

C There is a typo in Cioslowski's paper concerning N14/66. 
C Equation (2d) should be
C    N6665 + N6565 + 2 N5665 + N5655 + 2 N13/55 + 2 N13/56 -
C  4 N14/66 - 2 N14/56 = 0
      c2d=p(15)+p(18)+p(22)-4*p(27)+2*(p(20)+p(25)+p(26)-p(28))

      c2e=p(14)+p(15)+p(16)+p(17)+p(18)+p(19)+p(20)+p(21)+p(22)+
     1 p(23)+2*(p(24)+p(25)+p(26)+p(27)+p(28)+p(29))+3*p(30)

      Write(Iout,1001) c2a,c2b,c2c,c2d,c2e
      if(c2a.ne.0.or.c2b.ne.0.or.c2c.ne.0.or.c2d.ne.0.or.c2e.ne.60) 
     1 Write(Iout,1002) 
 1000 Format(' Severe Error: Sum of patches ',I5,' not equal to'
     1 ' number of hexagons ',I5)
 1001 Format(1X,'Check equations (2a)-(2e): Solutions ',5I3) 
 1002 Format(' Severe Error: Cioslowski equations not fulfilled')
      Return
      END

      SUBROUTINE Cioslowski(K656,nl565,IRing56,IRing6,N6Ring,
     1 NringC,NringD,NringE,NringF,n3ra,n3rc,Mpatt,E1,E2)
      use config
      IMPLICIT INTEGER (A-Z)
      Real*8 EMC(30),E1,E2
      DIMENSION MPatt(30),n3ra(3,Nmax),n3rc(3,Nmax)
      DIMENSION NringC(Emax),NringD(Emax)
      DIMENSION NringE(Emax),NringF(Emax)
      DIMENSION nhex(6),ihex(6),ipent(6)
      data EMC/30.336,18.636,13.560,6.351,-0.145,4.162,-2.468,
     1 -5.563,1.559,-0.807,-6.576,-0.133,-0.313,18.498,16.476,
     1 14.792,12.519,14.779,14.255,13.427,8.007,11.087,1.245,
     1 30.422,34.103,26.377,31.455,32.167,29.980,44.281/
      do I=1,30
       Mpatt(I)=0
      enddo

C Loop through all 6-rings and find 4 basic patterns
C depending the number of 5-rings connected to the 6-ring
      Do I=13,N6ring+12
C  Search for 5-rings
       n5count=0
       N135=0
       do J=1,6
        ipent(J)=0
        nhex(J)=0
       enddo
       Do J=1,IRing56
        if(NringD(J).eq.I) then
         n5count=n5count+1
         ipent(n5count)=NringC(J)
       endif
       enddo
       if(n5count.gt.3) stop 1
       go to (10,20,30,40),n5count+1

C  zero pentagons on main hexagon, all hexagons
C  Get hexagon ring numbers
   10  i6=0
       Do J=1,IRing6
        if(NringE(J).eq.I) then
         i6=i6+1
         if(i6.gt.6) stop 2
         nhex(i6)=NringF(J)
        endif
        if(NringF(J).eq.I) then
         i6=i6+1
         if(i6.gt.6) stop 3
         nhex(i6)=NringE(J)
        endif
       enddo
       if(i6.lt.6) stop 4
C Now sort hexagon numbers in ring of hexagons
C according to their adjacencies
       Call sorthex(number_vertices,IRing6,nhex,NringE,NringF)
C Now get the right structure motif. Search through c-6-5-6
       do K=1,6
        I1=K
        I2=K+1
        if(I2.eq.7) I2=1
        IH1=nhex(I1)
        IH2=nhex(I2)
        if(IH2.lt.IH1) then
         imem=IH1
         IH1=IH2
         IH2=imem
        endif
         ihex(K)=0
        do J=1,K656
         if(IH1.eq.n3ra(1,J).and.IH2.eq.n3ra(3,J)) then
          ihex(K)=1
          go to 19
         endif 
        enddo
   19   continue
       enddo
C Now analyze ihex and get M values for patterns 1-13
       Call HexPattern(ihex,M)
       go to 100

C---  One pentagon on main hexagon
   20  i6=0
       nhex(6)=0
       Do J=1,IRing6
        if(NringE(J).eq.I) then
         i6=i6+1
         if(i6.gt.5) stop 5
         nhex(i6)=NringF(J)
        endif
        if(NringF(J).eq.I) then
         i6=i6+1
         if(i6.gt.5) stop 6
         nhex(i6)=NringE(J)
        endif
       enddo
       if(i6.ne.5) stop 7
C Now sort hexagon numbers in ring of adjacent hexagons
       Call sorthex1(number_vertices,IRing6,IRing56,nhex,
     1  NRingC,NRingD,NringE,NringF,ipent)
C Now get the right structure motif. Search throgh c-6-5-6
       do K=1,4
        I1=K
        I2=K+1
        IH1=nhex(I1)
        IH2=nhex(I2)
        if(IH2.lt.IH1) then
         imem=IH1
         IH1=IH2
         IH2=imem
        endif
         ihex(K)=0
        do J=1,K656
          if(IH1.eq.n3ra(1,J).and.IH2.eq.n3ra(3,J)) then
          ihex(K)=1
          go to 29
         endif 
        enddo
   29   continue
       enddo
C Now analyze ihex and get M
       Call HexPattern5(ihex,M)
       go to 100

C---  Two pentagons on main hexagon
   30  i6=0
       nhex(5)=0
       nhex(6)=0
       Do J=1,IRing6
        if(NringE(J).eq.I) then
         i6=i6+1
         nhex(i6)=NringF(J)
        endif
        if(NringF(J).eq.I) then
         i6=i6+1
         nhex(i6)=NringE(J)
        endif
       enddo
       if(i6.ne.4) stop 10
C Now sort hexagon numbers in ascending order
      Do J1=1,4
      Do J2=J1+1,4
       if(nhex(J2).lt.nhex(J1)) then
        imem=nhex(J1)
        nhex(J1)=nhex(J2)
        nhex(J2)=imem
       endif
      enddo
      enddo
C  Sort if main 5-6-5 is linear or bent
       do J=1,nl565
        if(I.eq.n3rc(2,J)) go to 35
       enddo
C  bent 5-6-5 main pattern 
C   Get 3 adjacent hexagons
       Call sorthex3(number_vertices,IRing6,nhex,NringE,NringF)
C   Count 656 rings
       I1=nhex(1)
       I2=nhex(2)
       if(I1.gt.I2) then
        imem=I1
        I1=I2
        I2=imem
       endif
       I3=nhex(2)
       I4=nhex(3)
       if(I3.gt.I4) then
        imem=I3
        I3=I4
        I4=imem
       endif
       icount=0
C    Loop through closed 656
       do J=1,K656
        if(I1.eq.n3ra(1,J).and.I2.eq.n3ra(3,J)) icount=icount+1
        if(I3.eq.n3ra(1,J).and.I4.eq.n3ra(3,J)) icount=icount+1
       enddo
       if(icount.gt.2) stop 11
C 13/66
       if(icount.eq.0) M=24
C 13/65
       if(icount.eq.1) M=25
C 13/55
       if(icount.eq.2) M=26
       go to 100
C  linear 5-6-5 main pattern
C  Find adjacent hexagons
   35  Call sorthex4(number_vertices,IRing6,nhex,NringE,NringF)
C    Loop through closed 656
       I1=nhex(1)
       I2=nhex(2)
       I3=nhex(3)
       I4=nhex(4)
       icount=0
C    Loop through closed 656
       do J=1,K656
        if(I1.eq.n3ra(1,J).and.I2.eq.n3ra(3,J)) icount=icount+1
        if(I3.eq.n3ra(1,J).and.I4.eq.n3ra(3,J)) icount=icount+1
       enddo
       if(icount.gt.2) stop 12
C 14/66
       if(icount.eq.0) M=27
C 14/65
       if(icount.eq.1) M=28
C 14/55
       if(icount.eq.2) M=29
       go to 100

C---  Three pentagons on main hexagon
   40  M=30
       
  100  Mpatt(M)=Mpatt(M)+1
      enddo

C Get energy
      E1=0.d0
      E2=-8.050751d3/(dfloat(number_vertices)-30.050)
      Do I=1,30
       E1=E1+dfloat(Mpatt(I))*EMC(I)
      enddo

      Return
      END

      SUBROUTINE sorthex(number_vertices,IRing6,nhex,NringE,NringF)
      IMPLICIT INTEGER (A-Z)
      DIMENSION NringE(3*number_vertices/2),NringF(3*number_vertices/2)
      DIMENSION nhex(6)
C     Sort hexagon numbers according to their adjacencies
      do I=1,6
       IS=nhex(I)
      do J=I+1,6
       JS=nhex(J)
       do K=1,IRing6
        if(IS.eq.NringE(K).and.JS.eq.NringF(K)) then
         nhex(J)=nhex(I+1)
         nhex(I+1)=JS
         go to 10
        endif
        if(IS.eq.NringF(K).and.JS.eq.NringE(K)) then
         nhex(J)=nhex(I+1)
         nhex(I+1)=JS
         go to 10
        endif
       enddo
      enddo
  10  continue
      enddo
      Return
      END

      SUBROUTINE sorthex1(number_vertices,IRing6,Iring56,nhex,
     1 NRingC,NRingD,NringE,NringF,ipent)
      IMPLICIT INTEGER (A-Z)
C Sort hexagon numbers in ring of adjacent hexagons
C after the pentagon
      DIMENSION NringC(3*number_vertices/2),NringD(3*number_vertices/2)
      DIMENSION NringE(3*number_vertices/2),NringF(3*number_vertices/2)
      DIMENSION nhex(6),ipent(6)
C     Get hexagon adjacent to pentagon
      ifound=0
      do I=1,5
       IH=nhex(I)
       Do J=1,IRing56
        if(IH.eq.NringD(J).and.NRingC(J).eq.ipent(1)) then
         imem=nhex(1)
         nhex(1)=IH
         nhex(I)=imem
         ifound=1
         go to 1
        endif
       enddo
      enddo
   1  if(ifound.eq.0) stop 14
 
C     Now sort in spiral
      do I=1,5
       IS=nhex(I)
      do J=I+1,5
       JS=nhex(J)
       do K=1,IRing6
        if(IS.eq.NringE(K).and.JS.eq.NringF(K)) then
         nhex(J)=nhex(I+1)
         nhex(I+1)=JS
         go to 10
        endif
        if(IS.eq.NringF(K).and.JS.eq.NringE(K)) then
         nhex(J)=nhex(I+1)
         nhex(I+1)=JS
         go to 10
        endif
       enddo
      enddo
  10  continue
      enddo
      Return
      END

      SUBROUTINE sorthex3(number_vertices,IRing6,nhex,NringE,NringF)
      IMPLICIT INTEGER (A-Z)
      DIMENSION NringE(3*number_vertices/2),NringF(3*number_vertices/2),
     1 icount(4),nhex(6)
C     Find ring connected to two others
      do I=1,4
       icount(i)=0
      enddo
      do I=1,4
      do J=I+1,4
       do K=1,IRing6
        if(NringE(K).eq.nhex(I).and.NringF(K).eq.nhex(J)) then
         icount(I)=icount(I)+1
         icount(J)=icount(J)+1
        endif
       enddo
      enddo
      enddo
      idis=0
C     Eliminate disconnected one 
      do I=1,4
       if(icount(I).eq.0) then
        nhex(i)=nhex(4)
        nhex(4)=0
        icount(i)=icount(4)
        icount(4)=0
        idis=1
       endif
      enddo
      if(idis.eq.0) stop 19
C put middle ring at position 2
      do I=1,3
       if(icount(I).eq.2) ipos=I
      enddo
C Now swap
      imem=nhex(2)
      nhex(2)=nhex(ipos)
      nhex(ipos)=imem
      Return
      END

      SUBROUTINE sorthex4(number_vertices,IRing6,nhex,NringE,NringF)
      IMPLICIT INTEGER (A-Z)
      DIMENSION NringE(3*number_vertices/2),NringF(3*number_vertices/2),
     1  nhex(6)
      ifound=0
      do I=2,4
       do K=1,IRing6
        if(NringE(K).eq.nhex(1).and.NringF(K).eq.nhex(I)) then
         ifound=I
         go to 10
        endif
       enddo
      enddo
  10  if(ifound.eq.0) stop 13
       imem=nhex(2)
       nhex(2)=nhex(ifound)
       nhex(ifound)=imem
      if(nhex(3).gt.nhex(4)) then
       imem=nhex(3)
       nhex(3)=nhex(4)
       nhex(4)=imem
      endif
      Return
      END
      
      SUBROUTINE HexPattern(ihex,M)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ihex(6)
       isum=0
       do I=1,6
        isum=isum+ihex(i)
       enddo
       go to (1,2,3,4,5,6,7) isum+1
C 666666
    1  M=1
       Return
C 666665
    2  M=2
       Return
C 666655 or 666565 or 665665
C 665665
    3  M=5
       do I=1,6
        I1=I+1
        if(I.eq.6) I1=1
        if(ihex(I).eq.1.and.ihex(I1).eq.1) then
C 666655
         M=3
         Return
        endif
        I2=I+2
        if(I.eq.5) I2=1
        if(I.eq.6) I2=2
        if(ihex(I).eq.1.and.ihex(I2).eq.1) then
C 666565
         M=4
         Return
        endif
       enddo
       Return
C 666555 or 665655 or 656565
C Care needs to be taken as 665655 is chiral
C so this is taken as the first one
C 665655
    4  M=7
       do I=1,6
        I1=I+1
        I2=I+2
        if(I.eq.5) I2=1
        if(I.eq.6) then
         I1=1
         I2=2
        endif
        if(ihex(I).eq.1.and.ihex(I1).eq.1.and.ihex(I2).eq.1) then
C 666555
         M=6
         Return
        endif
        I4=I2+2
        if(I.eq.3) I4=1
        if(I.eq.4) I4=2
        if(ihex(I).eq.1.and.ihex(I2).eq.1.and.ihex(I4).eq.1) then
C 656565
         M=8
         Return
        endif
       enddo
       Return
C 665555 or 656555 or 655655
C 656555
    5  M=10
       do I=1,6
        I1=I+1
        I2=I+2
        I3=I+3
        if(I.eq.6) then
         I1=1
         I2=2
         I3=3
        endif
        if(I.eq.5) then
         I2=1
         I3=2
        endif
        if(I.eq.4) then
         I3=1
        endif
        if(ihex(I).eq.1.and.ihex(I1).eq.1.and.ihex(I2).eq.1.
     1   and.ihex(I3).eq.1) then
C 665555
         M=9
         Return
        endif
        I4=I2+2
        if(I.eq.3) I4=1
        if(I.eq.4) I4=2
        if(ihex(I).eq.1.and.ihex(I1).eq.1.and.ihex(I3).eq.1.
     1   and.ihex(I4).eq.1) then
C 655655
         M=11
         Return
        endif
       enddo
       Return
C 655555
    6  M=12
       Return
C 555555
    7  M=13
      Return
      END

      SUBROUTINE HexPattern5(ihex,M)
      IMPLICIT INTEGER (A-Z)
C Analyze ihex and get M
      DIMENSION ihex(6)
       isum=0
       do I=1,4
        isum=isum+ihex(i)
       enddo
       go to (1,2,3,4,5) isum+1
C 6666
    1  M=14
       Return
C 6665 or 6656
C 6656
    2  M=16
C 6665 
        if(ihex(1).eq.1.or.ihex(4).eq.1) M=15
       Return
C 6655 or 6565 or 6556 or 5665
C 6655
    3  M=17
       if(ihex(2).eq.1.and.ihex(4).eq.1.or.
     1  ihex(1).eq.1.and.ihex(3).eq.1) then
C 6565 
        M=18
        Return
       endif
       if(ihex(2).eq.1.and.ihex(3).eq.1) then
C 6556
        M=19
        Return
       endif
       if(ihex(1).eq.1.and.ihex(4).eq.1) then
C 5665
        M=20
        Return
       endif
       Return
C 6555 or 5655 
C 5655 
    4  M=22
       if(ihex(2).eq.1.and.ihex(3).eq.1.and.ihex(4).eq.1) then
C 6555
        M=21
        Return
       endif
       if(ihex(1).eq.1.and.ihex(2).eq.1.and.ihex(3).eq.1) then
C 6555
        M=21
        Return
       endif
       Return
C 5555
    5  M=23
      Return
      END

      SUBROUTINE StoneWales(Kring3,n3ra,numbersw,nSW)
      use config
      IMPLICIT INTEGER (A-Z)
      Integer n3ra(3,Nmax),nSW(4,66)
       numbersw=0 
      if(Kring3.eq.0) Return
      do I=1,Kring3
      do J=I+1,Kring3
       nI1=n3ra(1,I)
       nI3=n3ra(3,I)
       nJ1=n3ra(1,J)
       nJ3=n3ra(3,J)
       if(NI1.eq.NJ1.and.NI3.eq.NJ3) then
        numbersw=numbersw+1
        nsw(1,numbersw)=nI1
        nsw(4,numbersw)=nI3
        nsw(2,numbersw)=n3ra(2,I)
        nsw(3,numbersw)=n3ra(2,J)
       endif
      enddo
      enddo
      Return
      END
 
      SUBROUTINE Ring666(I6C,NrE,NrF,N6MEM,
     1 KRing3,LRing3b,LRing3l,n666,c666)
      use config
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION N6MEM(Mmax,6),NrE(Emax),NrF(Emax)
      DIMENSION IS(6),JS(6)
      Integer n666(3,Nmax*11),c666(3,Nmax)
C     Write out all 3-ring 566 connections first then sort them out
      ndim = Nmax*11
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
         if(Kring3.gt.Nmax) stop
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
 
      SUBROUTINE Ring566(I6C,I56C,NrC,NrD,NrE,NrF,
     1 N5MEM,N6MEM,KRing3,LRing3,n566,b566,l566)
      use config
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6)
      DIMENSION NrE(Emax),NrF(Emax)
      DIMENSION NrC(Emax),NrD(Emax)
      Integer n566(3,Nmax*11),l566(3,Nmax),b566(3,Nmax)
C     Write out all 3-ring 566 connections first then sort them out
      ndim = Nmax*11
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
 
      SUBROUTINE Ring656(I6C,I56C,NrC,NrD,NrE,NrF,
     1  KRing3,LRing3,n656,c656,o656)
      use config
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrE(Emax),NrF(Emax)
      DIMENSION NrC(Emax),NrD(Emax)
      Integer n656(3,Nmax*11),c656(3,Nmax),o656(3,Nmax)
C     Write out all 3-ring 656 connections first then sort them out
      ndim = Nmax*11
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
 
      SUBROUTINE Ring565(I5C,I56C,NrA,NrB,NrC,NrD,KRing3,LRing3,
     1 n565,b565,l565,N5MEM,N6MEM)
      use config
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6)
      DIMENSION NrA(Emax),NrB(Emax)
      DIMENSION NrC(Emax),NrD(Emax)
      Integer n565(3,Nmax*11),b565(3,Nmax),l565(3,Nmax)
C     Write out all 3-ring 565 connections first then sort them out
      ndim = Nmax*11
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
 
      SUBROUTINE Ring556(I5C,I56C,NrA,NrB,NrC,NrD,
     1 KRing3,LRing3,n556,c556,o556)
      use config
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrA(Emax),NrB(Emax)
      DIMENSION NrC(Emax),NrD(Emax)
      DIMENSION IS(6)
      Integer n556(3,Nmax*11),c556(3,Nmax),o556(3,Nmax)
C     Write out all 3-ring 556 connections first then sort them out
      ndim = Nmax*11
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
 
      SUBROUTINE Ring555(I5C,NrA,NrB,KRing3,LRing3,n555,n555f,m555)
      use config
C     Search for (o5-5-5) and (c5-5-5) 3-ring fusions
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrA(Emax),NrB(Emax)
      DIMENSION IS(6),JS(6)
      DIMENSION n555(3,Nmax*11),m555(3,Nmax),n555f(3,Nmax)
C     Write out all 3-ring 555 connections first then sort them out
      ndim = Nmax*11
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
 
      SUBROUTINE Ring56(IR56,N5R,N6R,NrA,NrB,Nring,
     1 N5MEM,N6MEM)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6),Nring(Mmax)
      DIMENSION NrA(Emax),NrB(Emax)
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
 
      SUBROUTINE Ring55(IR5,N5R,NrA,NrB,Nring,
     1 N5MEM,IedA,IedB)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION N5MEM(Mmax,5),Nring(Mmax)
      DIMENSION NrA(Emax),NrB(Emax)
      DIMENSION IedA(Emax),IedB(Emax)
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
 
      SUBROUTINE Ring66(IR6,N5R,N6R,NrA,NrB,
     1 Nring,N6MEM)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION N6MEM(Mmax,6),Nring(Mmax)
      DIMENSION NrA(Emax),NrB(Emax)
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
 
      SUBROUTINE DistMatrix(IOUT,Iprint,ireturn,Iopt,
     1 Dist,DistMat,Rmin,Rmax,Vol,ASphere)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax),DistMat(NmaxL)
C     Calculate distance matrix between atoms from cartesian coordinates
      IMat=0
      if(Iopt.eq.0.and.Iprint.eq.1) Write(IOUT,1000) 
      Do I=2,number_vertices
      Do J=1,I-1
       Sum=0.d0
      Do K=1,3
       Sum=Sum+(Dist(K,I)-Dist(K,J))**2
      enddo
       IMat=IMat+1
       DistMat(IMat)=dsqrt(Sum)
      enddo
      enddo
      Do I=2,number_vertices
      IZ=((I-1)*(I-2))/2
      if(Iopt.eq.0.and.Iprint.eq.1) 
     1   Write(IOUT,1001) (I,J,DistMat(IZ+J),J=1,I-1)
      enddo
      if(ireturn.eq.1) return

C     Determine minimum and maximum distances
      Rmin=1.d20
      Rmax=1.d0
      RmaxS=1.d0
      Do I=1,number_vertices
       Sum=0.d0
      Do K=1,3
       Sum=Sum+Dist(K,I)**2
      enddo
       Sumd=dsqrt(Sum)
       If(Sumd.gt.RmaxS) RmaxS=Sumd
      Do J=I+1,number_vertices
       DM=FunDistMat(I,J,DistMat)
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
      Vol=4.d0/3.d0*dpi*RmaxS**3
      ASphere=4.d0*dpi*RmaxS**2
      Write(IOUT,1002) Rmin,Imin,Jmin,Rmax,Imax,Jmax,RmaxS,Vol,
     1 ASphere,ASphere/Vol
C     Volume of ideal capped icosahedron
      If(number_vertices.eq.60) then
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
      If(number_vertices.eq.20) then
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
     1 ' between atoms ',I5,' and ',I5,
     1      /1X,'Maximum distance ',D15.8,
     1 ' between atoms ',I5,' and ',I5,
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

      SUBROUTINE Ring5(Ncount5,IN,I5,Mmax,IPA,N5MEM,N5MEMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPa(6,96),N5MEM(Mmax,5),N5MEMS(Mmax,5)
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

      SUBROUTINE Ring6(Ncount6,IN,I6,Mmax,IPA,N6MEM,N6MEMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPa(6,96),N6MEM(Mmax,6),N6MEMS(Mmax,6)
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

      SUBROUTINE Step(ITree,ILoop,IForbid,INum,IPa,IC3)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IC3(Nmax,3)
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

      SUBROUTINE Connect(MCon2,Ipent,IOUT,
     1  IC3,IDA,Tol,DistMat,Rmin)
      use config
C     Get the connectivities between 2 and 3 atoms
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DistMat(NmaxL)
      integer Icon2(3*number_vertices/2+1)
      integer IDA(Nmax,Nmax)
      integer NCI(12),NCJ(12)
      DIMENSION IC3(Nmax,3)
      Rtol=Rmin*(1.d0+Tol)
      Mcon2=0
c '+1' because num2 relies on finding a '0' after the end of the array.  wtf.
      Do I=1,3*number_vertices/2+1
        Icon2(I)=0
      enddo
      Do I=1,Nmax
        Do J=1,3
          IC3(I,J)=0
        enddo
      enddo
      if(Ipent.eq.0) then
        Do I=1,number_vertices
          Do J=I+1,number_vertices
            DM=FunDistMat(I,J,DistMat)
            If (DM.lt.Rtol) then
              Mcon2=Mcon2+1
              Ncount=I*number_vertices+J
              Icon2(Mcon2)=Ncount
            endif
          enddo
        enddo
      else
        Do I=1,number_vertices
          Do J=I+1,number_vertices
            If (IDA(I,J).eq.1) then
              Mcon2=Mcon2+1
              Ncount=I*number_vertices+J
              Icon2(Mcon2)=Ncount
            endif
          enddo
        enddo
      endif
      Write(IOUT,1000) Mcon2
      If(Mcon2.lt.number_vertices) then
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
          CALL Num2(Icon2(IArray),NCI(J),NCJ(J))
          If(NCI(J).Eq.0) then
            M12=J-1
            Go to 11
          endif
        enddo
   11   if(number_vertices.lt.100)
     1    Write(IOUT,1001) (NCI(J),NCJ(J),J=1,M12)
        if(number_vertices.ge.100.and.number_vertices.lt.1000) 
     1    Write(IOUT,1006) (NCI(J),NCJ(J),J=1,M12)
        if(number_vertices.ge.1000.and.number_vertices.lt.10000) 
     1    Write(IOUT,1007) (NCI(J),NCJ(J),J=1,M12)
        if(number_vertices.ge.10000) 
     1    Write(IOUT,1008) (NCI(J),NCJ(J),J=1,M12)
      enddo
      Write(IOUT,1002)
C     Get all vertices
      Do I=1,number_vertices
        IZ=0
        Do J=1,MCon2
          CALL Num2(Icon2(J),IX,IY)
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
   10   Continue
        Write(IOUT,1003) I,(IC3(I,J),J=1,3)
      enddo
C     Check if structure is alright at this point
      nexpedge=number_vertices*3/2
      if(Mcon2.ne.nexpedge) then
        Write(IOUT,1005) Mcon2,nexpedge
        stop
      endif
 1000 Format(/1X,' Number of connected surface atoms (edges, bonds): ',
     1 I5,/1X,' Connectivities (edge set):')
 1001 Format(1X,12('{',I2,',',I2,'} '))
 1002 Format(1X,' Calculate all vertices N and corresponding ',
     1 'adjacencies Ni of 3-connected graph:',
     1 /1X,'    N        N1   N2   N3')
 1003 Format(1X,I5,'    (',3I6,')')
 1004 Format(1X,'**** Error, not enough connected atoms',
     1 ' check coordinates')
 1005 Format(1X,'**** Severe error, number of edges (bonds) not as ',
     1 'expected from number of atoms: ',I4,' (expected: ',I4,')')
 1006 Format(1X,12('{',I3,',',I3,'} '))
 1007 Format(1X,12('{',I4,',',I4,'} '))
 1008 Format(1X,12('{',I5,',',I5,'} '))
      RETURN
      END
