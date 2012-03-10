      SUBROUTINE MoveCM(Natom,Matom,Iout,Iprint,IAtom,Dist,DistCM,El)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,natom),DistCM(3),Ainert(3,3),evec(3),df(3)
      DIMENSION IATOM(natom)
      CHARACTER*2 El(99)
      Data FNorm,STol/4.97369255d2,1.d-2/
      if(Iprint.ne.0) WRITE(IOUT,1000) 
      AnumX=0.d0
      AnumY=0.d0
      AnumZ=0.d0
      Do J=1,MAtom
        IM=IAtom(J)
        if(Iprint.ne.0) Write(IOUT,1002),J,IM,El(IM),(Dist(I,J),I=1,3)
        AnumX=AnumX+Dist(1,J)
        AnumY=AnumY+Dist(2,J)
        AnumZ=AnumZ+Dist(3,J)
      enddo
      Adenom=dfloat(MAtom)
      DistCM(1)=AnumX/Adenom
      If(dabs(DistCM(1)).lt.1.d-13) DistCM(1)=0.d0
       DistCM(2)=AnumY/Adenom
       If(dabs(DistCM(2)).lt.1.d-13) DistCM(2)=0.d0
         DistCM(3)=AnumZ/Adenom
         If(dabs(DistCM(3)).lt.1.d-13) DistCM(3)=0.d0
      Write(IOUT,1003) (DistCM(I),I=1,3)
      Do J=1,MAtom
       IM=IAtom(J)
         Do I=1,3
          Dist(I,J)=Dist(I,J)-DistCM(I)
         enddo
       Write(IOUT,1002),J,IM,El(IM),(Dist(I,J),I=1,3)
      enddo
C     Calculate distance part of moment of inertia
      Do I=1,3
      Do J=1,3
      Ainert(I,J)=0.d0
      enddo
      enddo
      Do J=1,MAtom
      x=Dist(1,J)
      y=Dist(2,J)
      z=Dist(3,J)
      x2=x*x
      y2=y*y
      z2=z*z
      xy=x*y
      xz=x*z
      yz=y*z
      Ainert(1,1)=Ainert(1,1)+y2+z2
      Ainert(2,2)=Ainert(2,2)+x2+z2
      Ainert(3,3)=Ainert(3,3)+x2+y2
      Ainert(1,2)=Ainert(1,2)+xy
      Ainert(1,3)=Ainert(1,3)+xz
      Ainert(2,3)=Ainert(2,3)+yz
      Ainert(2,1)=Ainert(1,2)
      Ainert(3,1)=Ainert(1,3)
      Ainert(3,2)=Ainert(2,3)
      enddo
      xmin=0.d0
      Do I=1,3
      if(Ainert(I,I).gt.xmin) xmin=Ainert(I,I)
      enddo
      Do I=1,3
      Do J=1,3
      Ainert(I,J)=Ainert(I,J)/xmin
      enddo
      enddo
C Diagonalize without producing eigenvectors
      call tred2l(Ainert,3,3,evec,df)
      call tqlil(evec,df,3,3)
      Do I=1,3
      evec(i)=evec(i)*xmin
      enddo
C Sort eigenvalues
      Do I=1,3
       e0=evec(I)
       jmax=I
        Do J=I+1,3
         e1=evec(J)
          if(e1.gt.e0) then 
           jmax=j
           e0=e1
          endif
        enddo
        if(i.ne.jmax) then
         ex=evec(jmax)
         evec(jmax)=evec(I)
         evec(I)=ex
        endif
      enddo
      Write(Iout,1001) (evec(i),I=1,3)
      do i=1,3
      evec(i)=evec(i)/FNorm
      enddo
      Write(Iout,1004) FNorm,(evec(i),I=1,3)
C Determine shape
      S1=dabs(evec(1)-evec(2))/evec(1)
      S2=dabs(evec(1)-evec(3))/evec(1)
      S3=S1+S2
C Symmetric top
      if(S3.lt.Stol) then
      Write(Iout,1005)
      Return
      endif
      if(S3.lt.5.d0*Stol) then
      Write(Iout,1006)
      Return
      endif
C Prolate
      dif=evec(1)/evec(3)
      if(S1.lt.Stol) then
      Write(Iout,1007)
      Return
      endif
      if(S1.lt.5.d0*Stol) then
      Write(Iout,1008)
      Return
      endif
C Oblate
      S4=dabs(evec(2)-evec(3))/evec(1)
      if(S4.lt.Stol) then
      Write(Iout,1009)
      Return
      endif
      if(S4.lt.5.d0*Stol) then
      Write(Iout,1010)
      Return
      endif
C Asymmetric
      Write(Iout,1011)
 1000 FORMAT(/1x,'Cartesian Input',
     1  /1X,'  I      Z  Element Cartesian Coordinates')
 1001 FORMAT(/1x,'Moment of inertia with setting the masses to unity:',
     1 /1x,'Eigenvalues (principal axis system): ',3F16.6)
 1002 FORMAT(1X,I3,1X,I6,1X,A2,6X,3(D18.12,2X))
 1003 FORMAT(/1X,'Shift Molecule to the centre of points:',
     1 /1X,'Original Centre: ',3(D15.9,1X),
     1 /1X,'New Coordinates:',
     1  /1X,'  I      Z  Element Cartesian Coordinates')
 1004 FORMAT(1x,'Using C60 ideal icosahedron to normalize eigenvalues',
     1 ' (',F6.2,' Angstroem^2)',/,' Eigenvalues (normalized): ',3F16.6)
 1005 FORMAT(1x,'Fullerene is symmetric top')
 1006 FORMAT(1x,'Fullerene is distorted symmetric top')
 1007 FORMAT(1x,'Fullerene is prolate')
 1008 FORMAT(1x,'Fullerene is distorted prolate')
 1009 FORMAT(1x,'Fullerene is oblate')
 1010 FORMAT(1x,'Fullerene is distorted oblate')
 1011 FORMAT(1x,'Fullerene is asymmetric')
      return
      END

      SUBROUTINE Diameter(Natom,M,IOUT,Dist,diam)
C Calculate largest and smallest atom-to-atom diameters
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Natom),diam(Natom),imirror(Natom),jmirror(Natom)
      Data difeps/1.d-6/
      Do i=1,M
C     Search for the closest carbon atom to dmirror
      distmin=1.d10
      Do k=1,M
      if(i.ne.k) then
      x=-Dist(1,i)-Dist(1,k)
      y=-Dist(2,i)-Dist(2,k)
      z=-Dist(3,i)-Dist(3,k)
      distm=dsqrt(x*x+y*y+z*z)
      if(distm.lt.distmin) then
      distmin=distm
      kmin=k
      endif
      endif
      enddo
      x=Dist(1,i)-Dist(1,kmin)
      y=Dist(2,i)-Dist(2,kmin)
      z=Dist(3,i)-Dist(3,kmin)
      diam(i)=dsqrt(x*x+y*y+z*z)
      imirror(i)=kmin
      enddo
C     Sort from the largest to the smallest diameter
C     and delete duplicates
      CALL SortR(Natom,M,Mnew,imirror,jmirror,diam)
      Write(IOUT,1000) Mnew
      If(Mnew.ne.M/2) Write(IOUT,1002) Mnew,M/2
      Write(IOUT,1001) (imirror(i),jmirror(i),diam(i),i=1,MNew)
      dif=diam(1)-diam(MNew)
      if(dif.lt.difeps) then
      Write(IOUT,1003)
      else
C     Calculate the moment of inertia tensor and diagonalize (no mass)
      endif
 1000 Format(/1X,' Atom-to-atom largest diameters obtained from ',
     1 'inversion through center of points sorted largest to '
     2 'smallest: (',I4,' values)')
 1001 Format(5(1X,'(',I4,',',I4,')',2X,D14.8))
 1002 Format(1X,' Error: Check subroutine diameter',2I4)
 1003 Format(1X,' Diameters indicate that points lie on a sphere')
      Return
      END
 
      SUBROUTINE CoordC60(natom,IN,Iout,IAtom,R5,R6,Dist)
      IMPLICIT REAL*8 (A-H,O-Z)
C     This routine constructs coordinates for the ideal capped icosahedron
C     as described by P. Senn, J. Chem. Ed. 72, 302 (1995)
      DIMENSION IATOM(natom),NUM(30),NUN(30)
      DIMENSION Dist(3,natom)
      DIMENSION DIco(3,12)
      Data NUM/1,1,1,1,1,2,2,2,2,3,3,3,4,4,4,5,5,5,
     1   6,6,6,7,7,7,8,8,8,9,9,10/
      Data NUN/2,3,9,10,11,3,4,10,12,4,5,11,5,6,12,
     1   6,7,11,7,8,12,8,9,11,9,10,12,10,11,12/
      DATA API/3.14159265358979d0/
      Write(Iout,1000) R5,R6
      DIF1=R5-R6
      Dif=Dabs(DIF1)
      If(Dif.lt.1.D-9) Write(Iout,1001)
C     Calculate the 12 coordinates of icosahedron
      dsqrt5=dsqrt(5.d0)
      AKSI=R5+.5d0*R6
      Rmax=AKSI*dsqrt((5.d0+dsqrt5)*.5d0)
      Rcov=.5d0*dsqrt(2.d0*(3.d0+dsqrt5)*AKSI**2+R6**2)
      Write(Iout,1002) Rmax,Rcov
      Fac1=Rmax/dsqrt5
      Do I=1,10
C     Angle multiple of 36 deg, or Pi/10., multiplied by 2
      ang=.2d0*dfloat(I)*API
      fcos=dcos(ang)
      fsin=dsin(ang)
      DIco(1,I)=2.d0*Fac1*fcos
      DIco(2,I)=2.d0*Fac1*fsin
      DIco(3,I)=Fac1*(-1.d0)**I
      enddo
      DIco(1,11)=0.d0
      DIco(1,12)=0.d0
      DIco(2,11)=0.d0
      DIco(2,12)=0.d0
      DIco(3,11)=-RMAX
      DIco(3,12)=RMAX
      Write(Iout,1003)
      Do I=1,12
      Write(Iout,1004) (DIco(J,I),J=1,3)
      enddo
C     Construct vertices of capped icosahedron
      fac2=R5/(2.d0*R5+R6)
C     Now calculated the coordinates
      ILoop=1
      Do I=1,30
      I1=NUM(I)
      I2=NUN(I)
      Do J=1,3
      Val1=DIco(J,I1)+fac2*(DIco(J,I2)-DIco(J,I1))
      Dist(J,ILoop)=Val1
      Val2=DIco(J,I2)+fac2*(DIco(J,I1)-DIco(J,I2))
      Dist(J,ILoop+1)=Val2
      enddo
      ILoop=ILoop+2
      enddo
 1000 Format(/1X,' Construct coordinates for 60-Fullerene',
     1 /1X,' Distances: R5= ',D12.6,'  R6= ',D12.6)
 1001 Format(/1X,' Identical distances chosen')
 1002 Format(/1X,' Radius of icosahedron (not capped): ',D12.6,
     1 ' Radius of covering central sphere for capped icosahedron: ',
     2 D12.6)
 1003 FORMAT(/1X,'Coordinates of Icosaeder:',/7X,
     1 'X',12X,'Y',12X,'Z')
 1004 Format(1X,3(D12.6,1X))
      Return
      END 

      FUNCTION FunDistMat(I,J,natomL,DistMat)
C     Unpack distance matrix value from linear vector
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension DistMat(natomL)
      I1=I
      J1=J
      IF(I.lt.J) then
      I1=J
      J1=I
      endif
      IMat=((I1-1)*(I1-2))/2+J1
      FunDistMat=DistMat(IMat)
      Return
      END
 
      FUNCTION FunAngleMat(I,J,K,natomL,DistMat)
C     Unpack distance matrix value from linear vector
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension DistMat(natomL)
      Data dpi/3.14159265358979d0/
      I1=I
      J1=J
      IF(I.lt.J) then
      I1=J
      J1=I
      endif
      IMat=((I1-1)*(I1-2))/2+J1
      R1=DistMat(IMat)
      I1=J
      J1=K
      IF(J.lt.K) then
      I1=K
      J1=J
      endif
      IMat=((I1-1)*(I1-2))/2+J1
      R2=DistMat(IMat)
      I1=I
      J1=K
      IF(I.lt.K) then
      I1=K
      J1=I
      endif
      IMat=((I1-1)*(I1-2))/2+J1
      RM=DistMat(IMat)
      Fac=(R1*R1+R2*R2-RM*RM)/(2.d0*R1*R2)
      FunAngleMat=180.d0*dacos(Fac)/dpi
      Return
      END 
