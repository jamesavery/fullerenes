      SUBROUTINE MoveCM(Matom,Iout,Iprint,IAtom,mirror,isort,
     1 nosort,SP,Dist,DistCM,El)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax),DistCM(3),Ainert(3,3),evec(3),df(3)
      DIMENSION IATOM(Nmax)
      CHARACTER*2 El(99)
      Data FNorm,STol/4.97369255d2,1.d-2/
      if(Iprint.ne.0) WRITE(IOUT,1000) 
      AnumX=0.d0
      AnumY=0.d0
      AnumZ=0.d0
      if(mirror.ne.0) then
       Write(IOUT,1012)
       Do J=1,MAtom
       Do I=1,3
        Dist(I,J)=-Dist(I,J)
       enddo
       enddo
      endif
      Do J=1,MAtom
        IM=IAtom(J)
        if(Iprint.ne.0) Write(IOUT,1002),J,IM,El(IM),(Dist(I,J),I=1,3)
        AnumX=AnumX+Dist(1,J)
        AnumY=AnumY+Dist(2,J)
        AnumZ=AnumZ+Dist(3,J)
      enddo

C     Convert to internal coordinates
      Call CartInt(Dist,Matom,Iout,isort)
      if(isort.ne.0.and.nosort.eq.0) return

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
C Calculate distortion parameter
      a2=(evec(1)-evec(2))**2
      b2=(evec(2)-evec(3))**2
      c2=(evec(3)-evec(1))**2
      SP=dsqrt(a2+b2+c2)/evec(1)
      Write(Iout,1004) FNorm,(evec(i),I=1,3),SP
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
     1  /1X,'  I       Z Element Cartesian Coordinates')
 1001 FORMAT(/1x,'Moment of inertia with setting the masses to unity:',
     1 /1x,'Eigenvalues (principal axis system): ',3(' ',D18.12))
 1002 FORMAT(1X,I4,1X,I6,1X,A2,6X,3(D18.12,2X))
 1003 FORMAT(/1X,'Shift Molecule to the centre of points:',
     1 /1X,'Original Centre: ',3(D15.9,1X),
     1 /1X,'New Coordinates:',
     1  /1X,'  I      Z  Element Cartesian Coordinates')
 1004 FORMAT(1x,'Using C60 ideal icosahedron to normalize eigenvalues',
     1 ' (',F6.2,' Angstroem^2)',/,' Eigenvalues (normalized): ',
     1 3(' ',D15.9),/1X,'Sphericity parameter normed to largest ',
     1 'rotational constant: ',D15.9)
 1005 FORMAT(1x,'Fullerene is symmetric top')
 1006 FORMAT(1x,'Fullerene is distorted symmetric top')
 1007 FORMAT(1x,'Fullerene is prolate')
 1008 FORMAT(1x,'Fullerene is distorted prolate')
 1009 FORMAT(1x,'Fullerene is oblate')
 1010 FORMAT(1x,'Fullerene is distorted oblate')
 1011 FORMAT(1x,'Fullerene is asymmetric')
 1012 FORMAT(1x,'Fullerene coordinates inverted to get mirror image') 
      return
      END

      SUBROUTINE Diameter(M,IOUT,Dist,diam)
      use config
C Calculate largest and smallest atom-to-atom diameters
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax),diam(Nmax),imirror(Nmax),jmirror(Nmax)
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
      CALL SortR(M,Mnew,imirror,jmirror,diam)
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
 
      SUBROUTINE CoordC20C60(Iout,MAtom,R5,R6,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     This routine constructs coordinates for the ideal capped icosahedron
C     as described by P. Senn, J. Chem. Ed. 72, 302 (1995)
C     or for a dodecahedron
      DIMENSION NUM(30),NUN(30)
      DIMENSION Dist(3,Nmax)
      DIMENSION DIco(3,12)
      Data NUM/1,1,1,1,1,2,2,2,2,3,3,3,4,4,4,5,5,5,
     1   6,6,6,7,7,7,8,8,8,9,9,10/
      Data NUN/2,3,9,10,11,3,4,10,12,4,5,11,5,6,12,
     1   6,7,11,7,8,12,8,9,11,9,10,12,10,11,12/

      dsqrt5=dsqrt(5.d0)
C   Capped Isosahedron
      If(Matom.ne.20) then
      Matom=60
      Write(Iout,1000) R5,R6
      DIF1=R5-R6
      Dif=Dabs(DIF1)
      If(Dif.lt.1.D-9) Write(Iout,1001)
C     Calculate the 12 coordinates of icosahedron
      AKSI=R5+.5d0*R6
      Rmax=AKSI*dsqrt((5.d0+dsqrt5)*.5d0)
      Rcov=.5d0*dsqrt(2.d0*(3.d0+dsqrt5)*AKSI**2+R6**2)
      Write(Iout,1002) Rmax,Rcov
      Fac1=Rmax/dsqrt5
      Do I=1,10
C     Angle multiple of 36 deg, or Pi/10., multiplied by 2
      ang=.2d0*dfloat(I)*dpi
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
      endif


C   Dodecahedron
      If(Matom.eq.20) then
      Write(Iout,1005) R5
      phi=(1.d0+dsqrt5)*.5d0
      edge=dsqrt5-1.d0
      phiinv=edge*.5d0
      fac=R5/edge
C     Atom 1
      Dist(1,1)=1.d0
      Dist(2,1)=1.d0
      Dist(3,1)=1.d0
C     Atom 2
      Dist(1,2)=-1.d0
      Dist(2,2)=1.d0
      Dist(3,2)=1.d0
C     Atom 3
      Dist(1,3)=1.d0
      Dist(2,3)=-1.d0
      Dist(3,3)=1.d0
C     Atom 4
      Dist(1,4)=1.d0
      Dist(2,4)=1.d0
      Dist(3,4)=-1.d0
C     Atom 5
      Dist(1,5)=1.d0
      Dist(2,5)=-1.d0
      Dist(3,5)=-1.d0
C     Atom 6
      Dist(1,6)=-1.d0
      Dist(2,6)=1.d0
      Dist(3,6)=-1.d0
C     Atom 7
      Dist(1,7)=-1.d0
      Dist(2,7)=-1.d0
      Dist(3,7)=1.d0
C     Atom 8
      Dist(1,8)=-1.d0
      Dist(2,8)=-1.d0
      Dist(3,8)=-1.d0
C     Atom 9
      Dist(1,9)=0.d0
      Dist(2,9)=phiinv
      Dist(3,9)=phi
C     Atom 10
      Dist(1,10)=0.d0
      Dist(2,10)=-phiinv
      Dist(3,10)=phi
C     Atom 11
      Dist(1,11)=0.d0
      Dist(2,11)=phiinv
      Dist(3,11)=-phi
C     Atom 12
      Dist(1,12)=0.d0
      Dist(2,12)=-phiinv
      Dist(3,12)=-phi
C     Atom 13
      Dist(1,13)=phi
      Dist(2,13)=0.d0
      Dist(3,13)=phiinv
C     Atom 14
      Dist(1,14)=-phi
      Dist(2,14)=0.d0
      Dist(3,14)=phiinv
C     Atom 15
      Dist(1,15)=phi
      Dist(2,15)=0.d0
      Dist(3,15)=-phiinv
C     Atom 16
      Dist(1,16)=-phi
      Dist(2,16)=0.d0
      Dist(3,16)=-phiinv
C     Atom 17
      Dist(1,17)=phiinv
      Dist(2,17)=phi
      Dist(3,17)=0.d0
C     Atom 18
      Dist(1,18)=-phiinv
      Dist(2,18)=phi
      Dist(3,18)=0.d0
C     Atom 19
      Dist(1,19)=phiinv
      Dist(2,19)=-phi
      Dist(3,19)=0.d0
C     Atom 20
      Dist(1,20)=-phiinv
      Dist(2,20)=-phi
      Dist(3,20)=0.d0

C     Scale to get R5 distance
      Do I=1,Matom
      Do J=1,3
       Dist(J,I)=Dist(J,I)*fac
      enddo
      enddo
      RCS=.5d0*dsqrt(3.d0)*phi*R5
      Write(Iout,1006) RCS
      endif
 1000 Format(/1X,'Construct coordinates for capped icosahedron IPR C60',
     1 /1X,' Distances: R5= ',D12.6,'  R6= ',D12.6)
 1001 Format(/1X,'Identical distances chosen')
 1002 Format(/1X,'Radius of icosahedron (not capped): ',D12.6,
     1 ' Radius of covering central sphere for capped icosahedron: ',
     2 D12.6)
 1003 FORMAT(/1X,'Coordinates of Icosaeder:',/7X,
     1 'X',12X,'Y',12X,'Z')
 1004 Format(1X,3(D12.6,1X))
 1005 Format(/1X,'Construct coordinates for dodecahedron C20',
     1 /1X,'Distance: R5= ',D12.6)
 1006 Format(1X,'Radius of minimum covering sphere: ',D12.6)
      Return
      END 

      FUNCTION FunDistMat(I,J,DistMat)
      use config
C     Unpack distance matrix value from linear vector
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension DistMat(NmaxL)
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
 
      FUNCTION FunAngleMat(I,J,K,DistMat)
      use config
C     Unpack distance matrix value from linear vector
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension DistMat(NmaxL)
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

      SUBROUTINE EndoKrotoTrans(Matom,IN,Iout,n565,NEK,
     1 ihueckel,IDA,N5MEM,N6MEM,A,evec,df,Dist,layout2D,distp,
     1 CDist)
      use config
C     use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION IDA(Nmax,Nmax)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6)
      DIMENSION NEK(3,Nmax),IP(2,66),KD(4,66),KSW(3)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
C     type(c_ptr) :: g, new_fullerene_graph
      Write(Iout,1000)
      if(n565.eq.0) then
       Write(Iout,1017)
       return
      endif
      Do J=1,66
      Do I=1,2
       IP(I,J)=0
      enddo
      enddo

C     read pentagon numbers
      Read(IN,*,Err=100,end=100) ((IP(I,J),I=1,2),J=1,66)
C     Check if any is a Endo-Kroto pattern
  100 ntrans=0
      Do J=1,66
       if(IP(1,J).eq.0) go to 10
       ntrans=ntrans+1
      enddo
      if(ntrans.lt.1) then
       Write(Iout,1016) 
       return
      endif
   10 Do J=1,ntrans
       if(IP(2,J).eq.0) then
        J1=IP(1,J)
        if(J1.gt.n565) stop 33
        IP(1,J)=NEK(1,J1)
        IP(2,J)=NEK(3,J1)
       endif
      enddo
      Write(Iout,1001) ntrans,((IP(I,J),I=1,2),J=1,ntrans)
      nswt=0
      Do J=1,ntrans
       I1=IP(1,J)
       I2=IP(2,J)
       If(I1.gt.I2) then
        I3=I1
        I1=I2
        I2=I3
       endif
       Do K=1,n565
        J1=NEK(1,K)
        J2=NEK(3,K)
        if(I1.eq.J1.and.I2.eq.J2) then
         nswt=nswt+1
         if(nswt.ne.k) then
          do L=1,3
          KSW(L)=NEK(L,K)
          NEK(L,k)=NEK(L,nswt)
          NEK(L,nswt)=KSW(L)
          enddo
         endif
        endif
       enddo
      enddo
       If(nswt.ne.ntrans) then
        Write(Iout,1002) nswt,ntrans
        Write(Iout,1003) ((NEK(I,J),I=1,3),J=1,n565)
        Write(Iout,1004) 
        Return
       endif
        Write(Iout,1007) nswt
        Write(Iout,1003) ((NEK(I,J),I=1,3),J=1,nswt)

C Perform Endo-Kroto transformation
C Find bonds between 5- and 6-rings to be deleted 
      Do I=1,nswt
      Do J=1,4
       KD(J,I)=0
      enddo
      enddo
      Do I=1,nswt
       IR1=NEK(1,I)
       IR2=NEK(2,I)-12
       IR3=NEK(3,I)
        ibond=0
        ibond1=0
        Do J1=1,5
         IAtom1=N5MEM(IR1,J1)
         IAtom3=N5MEM(IR3,J1)
        Do J2=1,6
         IAtom2=N6MEM(IR2,J2)
         if(IAtom1.eq.IAtom2) then
          ibond=ibond+1
          KD(ibond,I)=IAtom1
         endif
         if(IAtom3.eq.IAtom2) then
          ibond1=ibond1+1
          KD(ibond1+2,I)=IAtom3
         endif
        enddo
        enddo
      enddo
      ierr=0
      Do I=1,nswt
      Do J=1,4
       if(KD(J,I).eq.0) ierr=1
      enddo
      enddo
      if(ierr.eq.1) then
       Write(Iout,1018) 
       stop
      endif
C Delete edges
      Write(Iout,1019) 
      Write(Iout,1012) ((KD(J,I),J=1,4),I=1,nswt)
      Do I=1,nswt
       I1=KD(1,I)
       I2=KD(2,I)
       IDA(I1,I2)=0
       IDA(I2,I1)=0
       I3=KD(3,I)
       I4=KD(4,I)
       IDA(I3,I4)=0
       IDA(I4,I3)=0
      enddo
C Add edges
      Write(Iout,1020)
      Do I=1,nswt
       II=2*I 
       IAD1=MAtom+II-1
       IAD2=MAtom+II
       if(IAD2.gt.Nmax) then
        Write(Iout,1023) Nmax
        stop
       endif
       I1=KD(1,I)
       I2=KD(2,I)
       I3=KD(3,I)
       I4=KD(4,I)
       Write(Iout,1021) I1,IAD1,I2,IAD1,I3,IAD2,I4,IAD2,IAD1,IAD2
       IDA(IAD1,IAD2)=1
       IDA(IAD2,IAD1)=1
       IDA(I1,IAD1)=1
       IDA(IAD1,I1)=1
       IDA(I2,IAD1)=1
       IDA(IAD1,I2)=1
       IDA(I3,IAD2)=1
       IDA(IAD2,I3)=1
       IDA(I4,IAD2)=1
       IDA(IAD2,I4)=1
      enddo
       MAtom=IAD2
       Write(Iout,1022) MAtom

C Adjacency matrix constructed
C Now analyze the adjacency matrix if it is correct
      nsum=0
      Do I=1,MAtom
      isum=0
      Do J=1,MAtom
      isum=isum+IDA(I,J)
      enddo
      If(isum.ne.3) nsum=nsum+1
      enddo
      if(nsum.ne.0) then
      WRITE(Iout,1005) nsum,isum
      stop
      else
      WRITE(Iout,1006)
      endif
      Call Tutte(Matom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
 1000 Format(/1X,'Endo-Kroto insertion of 2 vertices:',
     1 /1X,'Read pentagon ring numbers (between 1-12)')
 1001 Format(/1X,'Number of Endo-Kroto insertions: ',I2,
     1 /1X,'Endo-Kroto pentagon ring numbers:',/,
     1 10(' (',I2,','I2,') '))
 1002 Format(/1X,I2,' Number of pentagon rings for Endo-Kroto ',
     1 'insertion does not match list, only ',I2,' found:')
 1003 Format(7(' (',I2,',',I5,',',I2,')'))
 1004 Format(/1X,'==> RETURN')
 1005 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1006 FORMAT(1X,'Graph checked, it is cubic')
 1007 Format(1X,'Perform ',I2,' Endo-Kroto transformations ',
     1 /1X,'Modifying adjacency matrix for rings (P,H,P):')
 1012 Format(1X,6(' (',I5,'-',I5,') (',I5,'-',I5,')'))
 1016 Format(/1X,'No input found ==> RETURN')
 1017 Format(/1X,'No Endo-Kroto pattern in this fullerene ==> RETURN')
 1018 Format(/1X,'Found less vertices then anticipated ==> STOP')
 1019 Format(/1X,'Delete the following edges between vertices:')
 1020 Format(/1X,'Add new edges between vertices:')
 1021 Format(1X,5(' (',I5,'-',I5,') '))
 1022 Format(/1X,'Number of vertices in Endo-Kroto transformed ',
     1 'fullerene: ',I5)
 1023 Format(/1X,'Number of vertices in fullerene over the current',
     1 'limit set by Nmaxs=',I5,' ==> STOP')
      Return
      END

      SUBROUTINE WirzSchwerd(Matom,IN,Iout,JERR,numberWS,
     1 IWS,nWS,ihueckel,IDA,N5MEM,N6MEM,IC3, 
     1 A,evec,df,Dist,layout2D,distp,Cdist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION IDA(Nmax,Nmax),IC3(Nmax,3)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6)
      DIMENSION nWS(5,8),IP(20),KWS(5),IBWS(8,8)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
      Write(Iout,1000)
      JERR=0
      Do I=1,10
       IP(I)=0
      enddo
C     read single hexagon numbers
      Read(IN,*,Err=100,end=100) (IP(I),I=1,66)

C     Check if any is a Wirz-Schwerdtfeger D2h 55-6-55 pattern
  100 ntrans=0
      Do I=1,10
       if(IP(I).eq.0) go to 10
       ntrans=ntrans+1
      enddo
      if(ntrans.lt.1) then
       Write(Iout,1016) 
       JERR=1
       return
      endif
   10 if(ntrans.gt.numberWS) then
       Write(Iout,1017) ntrans,numberWS
       JERR=1
       return
      endif
      Write(Iout,1001) ntrans,(IP(J),J=1,ntrans)
 
      if(IWS.eq.1) then
       nfound=0
       Do I=1,ntrans
        I1=IP(I)
        Do J=1,numberWS
         J1=nWS(3,J)
         if(I1.eq.J1) then
          nfound=nfound+1
          IP(nfound)=J
         endif
        enddo
       enddo
       if(nfound.ne.ntrans) then
        Write(Iout,1018)
        JERR=1
        return
       endif
      ntrans=nfound
      endif

C     Sort array nWS
      Do I=1,ntrans
       I1=IP(I)
       Do J=1,5
        KWS(J)=nWS(J,I1)
        nWS(J,I1)=nWS(J,I)
        nWS(J,I)=KWS(J)
       enddo
      enddo

C     Check for shared pentagons or hexagons which is not allowed
      if(ntrans.gt.1) then
       KERR1=0
       KERR2=0
       Do I=1,ntrans
       Do J=I+1,ntrans
        if(nWS(3,I).eq.nWS(3,J)) KERR1=1
        if(nWS(1,I).eq.nWS(1,J)) KERR1=2
        if(nWS(1,I).eq.nWS(2,J)) KERR1=2
        if(nWS(1,I).eq.nWS(4,J)) KERR1=2
        if(nWS(1,I).eq.nWS(5,J)) KERR1=2
        if(nWS(2,I).eq.nWS(1,J)) KERR1=2
        if(nWS(2,I).eq.nWS(2,J)) KERR1=2
        if(nWS(2,I).eq.nWS(4,J)) KERR1=2
        if(nWS(2,I).eq.nWS(5,J)) KERR1=2
        if(nWS(4,I).eq.nWS(1,J)) KERR1=2
        if(nWS(4,I).eq.nWS(2,J)) KERR1=2
        if(nWS(4,I).eq.nWS(4,J)) KERR1=2
        if(nWS(4,I).eq.nWS(5,J)) KERR1=2
        if(nWS(5,I).eq.nWS(1,J)) KERR1=2
        if(nWS(5,I).eq.nWS(2,J)) KERR1=2
        if(nWS(5,I).eq.nWS(4,J)) KERR1=2
        if(nWS(5,I).eq.nWS(5,J)) KERR1=2
        if(KERR1.gt.0) Write(Iout,1024) I,J,
     1   (nWS(J1,I),J1=1,5),(nWS(J1,J),J1=1,5)
        if(KERR2.gt.0) Write(Iout,1020) I,J,
     1   (nWS(J1,I),J1=1,5),(nWS(J1,J),J1=1,5)
       enddo
       enddo
       endif
       if(KERR1.gt.0.or.KERR2.gt.0) return
   
        Write(Iout,1007) ntrans
        Write(Iout,1003) ((nWS(I,J),I=1,5),J=1,ntrans)
    
C Perform Wirz-Schwerdtfeger 6-vertex insertion
C First find common vertices between pentagons and middle hexagon
       Do I=1,ntrans
       IH=nWS(3,I)-12
       icount=0
C Pentagon 1
       I1=nWS(1,I)
        do J=1,5
         IPN=N5MEM(I1,J)
        do K=1,6
         IHN=N6MEM(IH,K)
         if(IPN.eq.IHN) then
          icount=icount+1
          IBWS(I,icount)=IPN
         endif
        enddo
        enddo
C Pentagon 2
       I1=nWS(2,I)
        do J=1,5
         IPN=N5MEM(I1,J)
        do K=1,6
         IHN=N6MEM(IH,K)
         if(IPN.eq.IHN) then
          icount=icount+1
          IBWS(I,icount)=IPN
         endif
        enddo
        enddo
C Pentagon 3
       I1=nWS(4,I)
        do J=1,5
         IPN=N5MEM(I1,J)
        do K=1,6
         IHN=N6MEM(IH,K)
         if(IPN.eq.IHN) then
          icount=icount+1
          IBWS(I,icount)=IPN
         endif
        enddo
        enddo
C Pentagon 4
       I1=nWS(5,I)
        do J=1,5
         IPN=N5MEM(I1,J)
        do K=1,6
         IHN=N6MEM(IH,K)
         if(IPN.eq.IHN) then
          icount=icount+1
          IBWS(I,icount)=IPN
         endif
        enddo
        enddo
       if(icount.ne.8) stop 50
       enddo
C Swap or leave the lest two pentagons
C See if pentagon 2 is next to 3 or not
C Find vertex first
      do I=1,ntrans
       IVP=IBWS(I,4)
       IVN=IBWS(I,3)
       if(IVP.eq.IBWS(I,1).or.IVP.eq.IBWS(I,2)) then
        IVP=IBWS(I,3)
        IVN=IBWS(I,4)
       endif
C Find opposite vertex IVNN
       Do J=1,3
        KWS(J)=IC3(IVP,J)
       enddo
       Do J=1,3
        if(KWS(J).eq.IVN) KWS(J)=KWS(3)
       enddo
       ifound=0
       Do J=1,5
        IP3=N5MEM(IBWS(I,2),J)
        if(IP3.eq.KWS(1)) ifound=1
       enddo
       IVNN=KWS(1)
       if(ifound.eq.1) IVNN=KWS(2)
C Check if IVNN is in pentagon 3
       ifound=0
       Do J=1,5
        IP3=N5MEM(IBWS(I,4),J)
        if(IP3.eq.IVNN) ifound=1
       enddo
C Swap vertices if IVNN not found
       if(ifound.eq.0) then
        IVP1=IBWS(I,5)
        IVP2=IBWS(I,6)
        IBWS(I,5)=IBWS(I,7)
        IBWS(I,6)=IBWS(I,8)
        IBWS(I,7)=IVP1
        IBWS(I,8)=IVP2
       endif
      enddo

C Transform adjacency matrix
C Delete edges
      Nlimit=MAtom+6*ntrans
      if(Nlimit.gt.Nmax) then
       Write(Iout,1019) Nmax
      endif
      Write(Iout,1008)
      do I=1,ntrans
       IE1=IBWS(I,1)
       IE2=IBWS(I,2)
       IE3=IBWS(I,3)
       IE4=IBWS(I,4)
       IE5=IBWS(I,5)
       IE6=IBWS(I,6)
       IE7=IBWS(I,7)
       IE8=IBWS(I,8)
       Write(Iout,1013) IE1,IE2,IE3,IE4,IE5,IE6,IE7,IE8
       IDA(IE1,IE2)=0
       IDA(IE2,IE1)=0
       IDA(IE3,IE4)=0
       IDA(IE4,IE3)=0
       IDA(IE5,IE6)=0
       IDA(IE6,IE5)=0
       IDA(IE7,IE8)=0
       IDA(IE8,IE7)=0
      enddo
      Write(Iout,1009)
C Add edges
      do I=1,ntrans
       idim=MAtom+6*(I-1)
       IE1=IBWS(I,1)
       IE2=IBWS(I,2)
       IE3=IBWS(I,3)
       IE4=IBWS(I,4)
       IE5=IBWS(I,5)
       IE6=IBWS(I,6)
       IE7=IBWS(I,7)
       IE8=IBWS(I,8)
       IV1=idim+1
       IV2=idim+2
       IV3=idim+3
       IV4=idim+4
       IV5=idim+5
       IV6=idim+6
       IDA(IE1,IV1)=1
       IDA(IV1,IE1)=1
       IDA(IE2,IV1)=1
       IDA(IV1,IE2)=1
       IDA(IE3,IV2)=1
       IDA(IV2,IE3)=1
       IDA(IE4,IV2)=1
       IDA(IV2,IE4)=1
       IDA(IE5,IV3)=1
       IDA(IV3,IE5)=1
       IDA(IE6,IV3)=1
       IDA(IV3,IE6)=1
       IDA(IE7,IV4)=1
       IDA(IV4,IE7)=1
       IDA(IE8,IV4)=1
       IDA(IV4,IE8)=1
       IDA(IV1,IV5)=1
       IDA(IV5,IV1)=1
       IDA(IV2,IV6)=1
       IDA(IV6,IV2)=1
       IDA(IV3,IV6)=1
       IDA(IV6,IV3)=1
       IDA(IV4,IV5)=1
       IDA(IV5,IV4)=1
       IDA(IV5,IV6)=1
       IDA(IV6,IV5)=1
       Write(Iout,1010) IE1,IV1,IE2,IV1,IE3,IV2,IE4,IV2,
     1  IE5,IV3,IE6,IV3,IE7,IV4,IE8,IV4,IV1,IV5,IV2,IV6,
     1  IV3,IV6,IV4,IV5,IV5,IV6
      enddo
      MAtom=MAtom+6*ntrans

C Adjacency matrix constructed
C Now analyze the adjacency matrix if it is correct
      nsum=0
      Do I=1,MAtom
      isum=0
      Do J=1,MAtom
      isum=isum+IDA(I,J)
      enddo
      If(isum.ne.3) nsum=nsum+1
      enddo
      if(nsum.ne.0) then
      WRITE(Iout,1005) nsum,isum
      stop
      else
      WRITE(Iout,1006)
      endif
      Call Tutte(Matom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
 1000 Format(/1X,'W-S 6-vertex insertion to D2h 55-6-55',
     1 ' ring pattern:',/1X,'Read hexagon ring numbers or position',
     1 ' in list of patterns')
 1001 Format(/1X,'Input for W-S 6-vertex insertions: ',
     1 I2,' entries with numbers ',15(' ',I5))
 1003 Format(6(' (',I2,',',I2,',',I5,',',I2,',',I2,')'))
 1005 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1006 FORMAT(1X,'Graph checked, it is cubic')
 1007 Format(1X,'Perform ',I2,' W-S 6-vertex insertions',
     1 /1X,'Modifying adjacency matrix for rings (P,P,H,P,P):')
 1008 Format(1X,'Transform adjacency matrix',/1X,
     1 'Bonds to be added and deleted',/1X,
     1 'Delete edges:')
 1009 Format(1X,'Add edges: ')
 1010 Format(1X,7(' (',I5,','I5,') '),/1X,6(' (',I5,','I5,') '))
 1013 Format(1X,4(' (',I5,','I5,') '))
 1016 Format(/1X,'No input found ==> RETURN')
 1017 Format(/1X,'Number of insertions ',I2,' exceeds number of ',
     1 'possible W-S insertions ',I2,' ==> RETURN')
 1018 Format(/1X,'Hexagon numbers do not match list of ',
     1 'W-S list ==> RETURN')
 1019 Format(/1X,'Dimension of new adjacency matrix exceeds the ',
     1 'Nmax limit set at ',I5,' ==> STOP')
 1020 Format(1X,'Shared pentagons found between input pattern ',I2,
     1 ' and ',I2,: ' (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')',
     1 ' / (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')')
 1024 Format(1X,'Shared hexagons found between input pattern ',I2,
     1 ' and ',I2,: ' (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')',
     1 ' / (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')')

      Return
      END

      SUBROUTINE YoshidaFowler6(Matom,IN,Iout,JERR,
     1 numberfm,IYF,nfm,ihueckel,IDA,N6MEM,IC3,
     1 A,evec,df,Dist,layout2D,distp,CDist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION IDA(Nmax,Nmax),IC3(Nmax,3)
      DIMENSION N6MEM(Mmax,6)
      DIMENSION nFM(6,66),IP(66),KFM(6),IBFM(66,10)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
      Write(Iout,1000)
      JERR=0
      Do I=1,66
       IP(I)=0
      enddo

C     read n triple of hexagon numbers
      Read(IN,*,Err=100,end=100) (IP(I),I=1,66)

C     Check if any is a Yoshida-Fowler D3h 666555 pattern
  100 ntrans=0
      Do I=1,66
       if(IP(I).eq.0) go to 10
       ntrans=ntrans+1
      enddo
      if(ntrans.lt.1) then
       Write(Iout,1016) 
       JERR=1
       return
      endif
   10 if(ntrans.gt.numberFM) then
       Write(Iout,1017) ntrans,numberFM
       JERR=1
       return
      endif
      ntr=ntrans/3
      ntr3=ntr*3
      if(ntr3.ne.ntrans) then
       Write(Iout,1022) 
       JERR=1
       return
      endif
      Write(Iout,1001) ntrans,(IP(J),J=1,ntrans)
 
      if(IYF.eq.1) then
       nfound=0
       Do I=1,ntrans,3
        if(IP(I).gt.IP(I+1)) then
         mem=IP(I+1)
         IP(I+1)=IP(I)
         IP(I)=mem
        endif
        if(IP(I).gt.IP(I+2)) then
         mem=IP(I+2)
         IP(I+2)=IP(I)
         IP(I)=mem
        endif
        if(IP(I+1).gt.IP(I+2)) then
         mem=IP(I+2)
         IP(I+2)=IP(I+1)
         IP(I+1)=mem
        endif
        I1=IP(I)
        I2=IP(I+1)
        I3=IP(I+2)
        Do J=1,numberFM
         J1=nFM(1,J)
         J2=nFM(2,J)
         J3=nFM(3,J)
         if(I1.eq.J1.and.I2.eq.J2.and.I3.eq.J3) then
          nfound=nfound+1
          IP(nfound)=J
         endif
        enddo
       enddo
       if(nfound*3.ne.ntrans) then
        Write(Iout,1018)
        JERR=1
        return
       endif
      ntrans=nfound
      endif

C     Sort array nFM
      Do I=1,ntrans
       I1=IP(I)
       Do J=1,6
        KFM(J)=nFM(J,I1)
        nFM(J,I1)=nFM(J,I)
        nFM(J,I)=KFM(J)
       enddo
      enddo

C     Check for shared pentagons or hexagons which is not allowed
      if(ntrans.gt.1) then
       KERR=0
       Do I=1,ntrans
       Do J=I+1,ntrans
       Do K1=1,3
       Do K2=1,3
        if(nFM(K1,I).eq.nFM(K2,J)) KERR=2
        if(nFM(K1+3,I).eq.nFM(K2+3,J)) KERR=1
        if(KERR.gt.0) then
         if(KERR.eq.2) Write(Iout,1024) I,J,
     1    (nFM(J1,I),J1=1,6),(nFM(J1,J),J1=1,6)
         if(KERR.eq.1) Write(Iout,1020) I,J,
     1    (nFM(J1,I),J1=1,6),(nFM(J1,J),J1=1,6)
         JERR=1
         Return
        endif
       enddo
       enddo
       enddo
       enddo
       endif

        Write(Iout,1007) ntrans
        Write(Iout,1003) ((nFM(I,J),I=1,6),J=1,ntrans)
    
C Perform Yoshida-Fowler 6-vertex insertion
C First find center vertex and adjacent vertices
        Do I=1,66
        Do J=1,10 
         IBFM(I,J)=0 
        enddo
        enddo
      Do I=1,ntrans
       IR1=nFM(1,I)-12
       IR2=nFM(2,I)-12
       IR3=nFM(3,I)-12
        Do J1=1,6
         IAtom1=N6MEM(IR1,J1)
        Do J2=1,6
         IAtom2=N6MEM(IR2,J2)
        Do J3=1,6
         IAtom3=N6MEM(IR3,J3)
         if(IAtom1.eq.IAtom2.and.IAtom2.eq.IAtom3) then
          IBFM(I,1)=IAtom1
          go to 20
         endif
        enddo
        enddo
        enddo
  20  if(IAtom1.eq.0) stop 33
       IBFM(I,2)=IC3(IAtom1,1)
       IBFM(I,3)=IC3(IAtom1,2)
       IBFM(I,4)=IC3(IAtom1,3)
C Now find outer 6 vertices between 5- and 6-rings
       icount=4
       do J=1,3
        NC=IBFM(I,J+1)
       do K=1,3
        if(IC3(NC,K).ne.IAtom1)then
         icount=icount+1
         IBFM(I,icount)=IC3(NC,K)
        endif
       enddo
       enddo
       if(icount.ne.10) stop 34
      enddo
      Write(Iout,1023) (IBFM(I,1),I=1,ntrans)

C Transform adjacency matrix
C Delete edges
      Nlimit=MAtom+6*ntrans
      if(Nlimit.gt.Nmax) then
       Write(Iout,1019) Nmax
      endif
      Write(Iout,1008)
      do I=1,ntrans
       IE1=IBFM(I,1)
       IE2=IBFM(I,2)
       IE3=IBFM(I,3)
       IE4=IBFM(I,4)
       IE5=IBFM(I,5)
       IE6=IBFM(I,6)
       IE7=IBFM(I,7)
       IE8=IBFM(I,8)
       IE9=IBFM(I,9)
       IE10=IBFM(I,10)
       Write(Iout,1013)IE2,IE5,IE2,IE6,IE3,IE7,IE3,IE8,IE4,IE9,IE4,IE10
       IDA(IE2,IE5)=0
       IDA(IE5,IE2)=0
       IDA(IE2,IE6)=0
       IDA(IE6,IE2)=0
       IDA(IE3,IE7)=0
       IDA(IE7,IE3)=0
       IDA(IE3,IE8)=0
       IDA(IE8,IE3)=0
       IDA(IE4,IE9)=0
       IDA(IE9,IE4)=0
       IDA(IE4,IE10)=0
       IDA(IE10,IE4)=0
C Add edges
       idim=MAtom+6*(I-1)
       IV1=idim+1
       IV2=idim+2
       IV3=idim+3
       IV4=idim+4
       IV5=idim+5
       IV6=idim+6
       IDA(IE5,IV1)=1
       IDA(IV1,IE5)=1
       IDA(IE6,IV2)=1
       IDA(IV2,IE6)=1
       IDA(IE7,IV3)=1
       IDA(IV3,IE7)=1
       IDA(IE8,IV4)=1
       IDA(IV4,IE8)=1
       IDA(IE9,IV5)=1
       IDA(IV5,IE9)=1
       IDA(IE10,IV6)=1
       IDA(IV6,IE10)=1
       IDA(IV1,IV2)=1
       IDA(IV2,IV1)=1
       IDA(IV3,IV4)=1
       IDA(IV4,IV3)=1
       IDA(IV5,IV6)=1
       IDA(IV6,IV5)=1
C Now start merry go round. Start with vertex IE6
       IDA(IV1,IE2)=1
       IDA(IE2,IV1)=1
       IDA(IV2,IE3)=1
       IDA(IE3,IV2)=1
C Determine to which hexagon IE6 belongs to
       Do J=1,3
        I1=nFM(J,I)-12
        Do K=1,6
         if(IE6.eq.N6MEM(I1,K)) then
          IX=I1
          go to 30
         endif
        enddo
       enddo
   30  ifound=0
       Do K=1,6
        if(IE7.eq.N6MEM(IX,K)) ifound=7
        if(IE8.eq.N6MEM(IX,K)) ifound=8
        if(IE9.eq.N6MEM(IX,K)) ifound=9
        if(IE10.eq.N6MEM(IX,K)) ifound=10
       enddo
       if(ifound.eq.0) stop 35
       if(ifound.eq.7.or.ifound.eq.8) then
         INext=9
        if(ifound.eq.7) then
         IDA(IV3,IE3)=1
         IDA(IE3,IV3)=1
         IDA(IV4,IE4)=1
         IDA(IE4,IV4)=1
         IADD1=IE3
         IADD2=IV3
         IADD3=IE4
         IADD4=IV4
         Ivert=IE8
        else
         IDA(IV4,IE3)=1
         IDA(IE3,IV4)=1
         IDA(IV3,IE4)=1
         IDA(IE4,IV3)=1
         IADD1=IE3
         IADD2=IV4
         IADD3=IE4
         IADD4=IV3
         Ivert=IE7
        endif
       else
         INext=7
        if(ifound.eq.9) then
         IDA(IV5,IE3)=1
         IDA(IE3,IV5)=1
         IDA(IV6,IE4)=1
         IDA(IE4,IV6)=1
         IADD1=IE3
         IADD2=IV5
         IADD3=IE4
         IADD4=IV6
         Ivert=IE10
        else
         IDA(IV6,IE3)=1
         IDA(IE3,IV6)=1
         IDA(IV5,IE4)=1
         IDA(IE4,IV5)=1
         IADD1=IE3
         IADD2=IV6
         IADD3=IE4
         IADD4=IV5
         Ivert=IE9
        endif
       endif
C Determine to which hexagon Ivert belongs to
       Do J=1,3
        I1=nFM(J,I)-12
        Do K=1,6
         if(Ivert.eq.N6MEM(I1,K)) then
          IX=I1
          go to 40
         endif
        enddo
       enddo
   40  ifound=0
       Do K=1,6
        if(IBFM(I,INext).eq.N6MEM(IX,K)) ifound=INext
        if(IBFM(I,INext+1).eq.N6MEM(IX,K)) ifound=INext+1
       enddo
       if(ifound.eq.0) stop 36
       if(ifound.eq.7.or.ifound.eq.8) then
        if(ifound.eq.7) then
         IDA(IV3,IE4)=1
         IDA(IE4,IV3)=1
         IDA(IV4,IE2)=1
         IDA(IE2,IV4)=1
         IADD5=IE4
         IADD6=IV3
         IADD7=IE2
         IADD8=IV4
        else
         IDA(IV4,IE4)=1
         IDA(IE4,IV4)=1
         IDA(IV3,IE2)=1
         IDA(IE2,IV3)=1
         IADD5=IE4
         IADD6=IV4
         IADD7=IE2
         IADD8=IV3
        endif
       else
        if(ifound.eq.9) then
         IDA(IV5,IE4)=1
         IDA(IE4,IV5)=1
         IDA(IV6,IE2)=1
         IDA(IE2,IV6)=1
         IADD5=IE4
         IADD6=IV5
         IADD7=IE2
         IADD8=IV6
        else
         IDA(IV6,IE4)=1
         IDA(IE4,IV6)=1
         IDA(IV5,IE2)=1
         IDA(IE2,IV5)=1
         IADD5=IE4
         IADD6=IV6
         IADD7=IE2
         IADD8=IV5
        endif
       endif
       Write(Iout,1009) IE5,IV1,IE6,IV2,IE7,IV3,IE8,IV4,IE9,IV5,IE10,
     1 IV6,IV1,IV2,IV3,IV4,IV5,IV6,IE2,IV1,IE3,IV2,IADD1,IADD2,IADD3,
     1 IADD4,IADD5,IADD6,IADD7,IADD8
      enddo
      MAtom=MAtom+6*ntrans

C Adjacency matrix constructed
C Now analyze the adjacency matrix if it is correct
      nsum=0
      Do I=1,MAtom
      isum=0
      Do J=1,MAtom
      isum=isum+IDA(I,J)
      enddo
      If(isum.ne.3) nsum=nsum+1
      enddo
      if(nsum.ne.0) then
      WRITE(Iout,1005) nsum,isum
      stop
      else
      WRITE(Iout,1006)
      endif
      Call Tutte(Matom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
 1000 Format(/1X,'Yoshida-Fowler 6-vertex insertion to D3h 666555 ',
     1 'ring pattern:',/1X,'Read hexagon ring numbers or position',
     1 ' in list of patterns')
 1001 Format(/1X,'Input for Yoshida-Fowler 6-vertex insertions: ',I2,
     1 ' entries with numbers ',15(' ',I5))
 1003 Format(6(' (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')'))
 1005 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1006 FORMAT(1X,'Graph checked, it is cubic')
 1007 Format(1X,'Perform ',I2,' Yoshida-Fowler 6-vertex insertions ',
     1 /1X,'Modifying adjacency matrix for rings (H,H,H,P,P,P):')
 1008 Format(1X,'Transform adjacency matrix',/1X,
     1 'Yoshida-Fowler bonds to be added and deleted')
 1009 Format(1X,'Add edges: ',8(' (',I5,'-',I5,') '),
     1 /1X,7(' (',I5,'-',I5,') '))
 1013 Format(1X,'Delete edges: '6(' (',I5,','I5,') '))
 1016 Format(/1X,'No input found ==> RETURN')
 1017 Format(/1X,'Number of insertions ',I2,' exceeds number of ',
     1 'possible Yoshida-Fowler insertions ',I2,' ==> RETURN')
 1018 Format(/1X,'Hexagon numbers do not match list of ',
     1 'Yoshida-Fowler list ==> RETURN')
 1019 Format(/1X,'Dimension of new adjacency matrix exceeds the ',
     1 'Nmax limit set at ',I5,' ==> STOP')
 1020 Format(1X,'Shared pentagons found between input pattern ',I2,
     1 ' and ',I2,: ' (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')',
     1 ' / (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')')
 1022 Format(1X,'Error: Need 3 hexagon numbers per pattern ==> RETURN')
 1023 Format(1X,'Central vertices: ',20I6)
 1024 Format(1X,'Shared hexagons found between input pattern ',I2,
     1 ' and ',I2,: ' (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')',
     1 ' / (',I5,',',I5,',',I5,',',I2,',',I2,',',I2,')')
      Return
      END

      SUBROUTINE YoshidaFowler4(Matom,IN,Iout,JERR,
     1 numberfm,IYF,nfm,ihueckel,IDA,N5MEM,N6MEM,
     1 A,evec,df,Dist,layout2D,distp,CDist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION IDA(Nmax,Nmax)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6)
      DIMENSION nFM(4,66),IP(66),KFM(4),IBFM(66,6)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
      Write(Iout,1000)
      JERR=0
      Do I=1,66
       IP(I)=0
      enddo

C     read hexagon numbers
      Read(IN,*,Err=100,end=100) (IP(I),I=1,66)

C     Check if any is a Yoshida-Fowler D3h 6555 pattern
  100 ntrans=0
      Do I=1,66
       if(IP(I).eq.0) go to 10
       ntrans=ntrans+1
      enddo
      if(ntrans.lt.1) then
       Write(Iout,1016) 
       JERR=1
       return
      endif
   10 if(ntrans.gt.numberFM) then
       Write(Iout,1017) ntrans,numberFM
       JERR=1
       return
      endif 
      if(IYF.eq.1) then
       nfound=0
       Do I=1,ntrans
        I1=IP(I)
        Do J=1,numberFM
         J1=nFM(1,J)
         if(I1.eq.J1) nfound=nfound+1
        enddo
       enddo
       if(nfound.ne.ntrans) then
        Write(Iout,1018)
        JERR=1
        return
       endif
      else
        Do I=1,ntrans
         I1=nFM(1,IP(I))
         IP(I)=I1
        enddo
      endif
      Write(Iout,1001) ntrans,(IP(J),J=1,ntrans)

C     Sort array nFM
      nFMt=0
      Do I=1,ntrans
       I1=IP(I)
       Do J=1,numberFM
        J1=nFM(1,J)
        if(I1.eq.J1) then
         nFMt=nFMt+1
         if(nFMt.ne.J) then
          do L=1,4
          KFM(L)=nFM(L,J)
          nFM(L,J)=nFM(L,nFMt)
          nFM(L,nFMt)=KFM(L)
          enddo
         endif
        endif
       enddo
      enddo

C     Check for shared pentagons which is not allowed
      if(nFMt.gt.1) then
       KERR=0
       Do I=1,nFMt
       Do J=I+1,nFMt
        Do K1=2,4
        Do K2=2,4
         if(nFM(K1,I).eq.nFM(K2,J)) then
          Write(Iout,1020) nFM(1,I),nFM(1,J),nFM(K1,I)
          KERR=1
         endif
        enddo
        enddo
       enddo
       enddo

C Perform Yoshida-Fowler 4-vertex insertion
C Find bonds between the middle 6-ring and the three 5-rings
       If(nFMt.ne.ntrans) then
        Write(Iout,1002) nFMt,ntrans
        Write(Iout,1003) ((nFM(I,J),I=1,4),J=1,nFMt)
        Write(Iout,1004)
        JERR=1 
        Return
       endif
        Write(Iout,1007) nFMt
        Write(Iout,1003) ((nFM(I,J),I=1,4),J=1,nFMt)
       if(KERR.ne.0) then
        Write(Iout,1021)
        JERR=1 
        Return
       endif
      endif
      Do I=1,nFMt
       IR1=nFM(1,I)-12
       IR2=nFM(2,I)
       IR3=nFM(3,I)
       IR4=nFM(4,I)
        ibond2=0
        ibond3=0
        ibond4=0
        Do J1=1,6
         IAtom1=N6MEM(IR1,J1)
        Do J2=1,5
         IAtom2=N5MEM(IR2,J2)
         IAtom3=N5MEM(IR3,J2)
         IAtom4=N5MEM(IR4,J2)
         if(IAtom1.eq.IAtom2) then
          ibond2=ibond2+1
          IBFM(I,ibond2)=IAtom2
         endif
         if(IAtom1.eq.IAtom3) then
          ibond3=ibond3+1
          IBFM(I,ibond3+2)=IAtom3
         endif
         if(IAtom1.eq.IAtom4) then
          ibond4=ibond4+1
          IBFM(I,ibond4+4)=IAtom4
         endif
        enddo
        enddo
        ibond=ibond2+ibond3+ibond4
        If(ibond.ne.6) then
         Write(Iout,1012)
        endif
        if(IBFM(I,1).gt.IBFM(I,2)) then
        IX=IBFM(I,1)
        IBFM(I,1)=IBFM(I,2)
        IBFM(I,2)=IX
        endif
        if(IBFM(I,3).gt.IBFM(I,4)) then
        IX=IBFM(I,3)
        IBFM(I,3)=IBFM(I,4)
        IBFM(I,4)=IX
        endif
        if(IBFM(I,5).gt.IBFM(I,6)) then
        IX=IBFM(I,5)
        IBFM(I,5)=IBFM(I,6)
        IBFM(I,6)=IX
        endif
      enddo

C Transform adjacency matrix
      Write(Iout,1008)
      Write(Iout,1013) ((IBFM(I,J),J=1,6),I=1,nFMt)
      Write(Iout,1015)
      do I=1,nFMt
       IE1=IBFM(I,1)
       IE2=IBFM(I,2)
       IE3=IBFM(I,3)
       IE4=IBFM(I,4)
       IE5=IBFM(I,5)
       IE6=IBFM(I,6)
       IDA(IE1,IE2)=0
       IDA(IE2,IE1)=0
       IDA(IE3,IE4)=0
       IDA(IE4,IE3)=0
       IDA(IE5,IE6)=0
       IDA(IE6,IE5)=0
       idim=MAtom+4*(I-1)
       IV1=idim+1
       IV2=idim+2
       IV3=idim+3
       IV4=idim+4
       if(IV4.gt.Nmax) then
        Write(Iout,1019)
        stop
       endif
       IDA(IE1,IV1)=1
       IDA(IV1,IE1)=1
       IDA(IE2,IV1)=1
       IDA(IV1,IE2)=1
       IDA(IE3,IV2)=1
       IDA(IV2,IE3)=1
       IDA(IE4,IV2)=1
       IDA(IV2,IE4)=1
       IDA(IE5,IV3)=1
       IDA(IV3,IE5)=1
       IDA(IE6,IV3)=1
       IDA(IV3,IE6)=1
       IDA(IV1,IV4)=1
       IDA(IV4,IV1)=1
       IDA(IV2,IV4)=1
       IDA(IV4,IV2)=1
       IDA(IV3,IV4)=1
       IDA(IV4,IV3)=1
       Write(Iout,1009) IE1,IV1,IE2,IV1,IE3,IV2,IE4,IV2,IE5,IV3,IE6,IV3,
     1 IV1,IV4,IV2,IV4,IV3,IV4
      enddo
       MAtom=MAtom+4*nFMt 

C Adjacency matrix constructed
C Now analyze the adjacency matrix if it is correct
      nsum=0
      Do I=1,MAtom
      isum=0
      Do J=1,MAtom
      isum=isum+IDA(I,J)
      enddo
      If(isum.ne.3) nsum=nsum+1
      enddo
      if(nsum.ne.0) then
      WRITE(Iout,1005) nsum,isum
      stop
      else
      WRITE(Iout,1006)
      endif
      Call Tutte(Matom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
 1000 Format(/1X,'Yoshida-Fowler 4-vertex insertion to D3h 6555 ',
     1 'ring pattern:',/1X,'Read hexagon ring numbers or position',
     1 ' in list of patterns')
 1001 Format(/1X,'Number of Yoshida-Fowler 4-vertex insertions: ',I2,
     1 /1X,'Yoshida-Fowler hexagon ring numbers:',20(' ',I5))
 1002 Format(/1X,I2,' Number of hexagon rings for Yoshida-Fowler ',
     1 'insertion does not match list, only ',I2,' found:')
 1003 Format(7(' (',I5,',',I2,',',I2,',',I2,')'))
 1004 Format(/1X,'==> RETURN')
 1005 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1006 FORMAT(1X,'Graph checked, it is cubic')
 1007 Format(1X,'Perform ',I2,' Yoshida-Fowler 4-vertex insertions ',
     1 /1X,'Modifying adjacency matrix for rings (H,P,P,P):')
 1008 Format(1X,'Transform adjacency matrix:',/1X,
     1 'Delete edges:')
 1009 Format(1X,9(' (',I5,'-',I5,') '))
 1012 Format(1X,'Error: did not find the 3 edges in ring connections ')
 1013 Format(1X,'Yoshida-Fowler bonds to be deleted. ',
     1 'Delete edges:',/,6(' (',I5,','I5,') '))
 1015 Format(1X,'Add edges:')
 1016 Format(/1X,'No input found ==> RETURN')
 1017 Format(/1X,'Number of insertions ',I2,' exceeds number of ',
     1 'possible Yoshida-Fowler insertions ',I2,' ==> RETURN')
 1018 Format(/1X,'Hexagon numbers do not match list of ',
     1 'Yoshida-Fowler list ==> RETURN')
 1019 Format(/1X,'Dimension of new adjacency matrix exceeds the ',
     1 'Nmax limit set at ',I5,' ==> STOP')
 1020 Format(1X,'Shared pentagons found between hexagons ',I5,
     1 ' and ',I5,'. Corresponding pentagon numbers: ',I2)
 1021 Format(1X,'Abort Yoshida-Fowler insertion ==> RETURN')
      Return
      END

      SUBROUTINE StoneWalesTrans(Matom,IN,Iout,numbersw,
     1 nSW,ihueckel,IDA,N6MEM,IC3,A,evec,df,Dist,layout2D,distp,
     1 CDist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION IDA(Nmax,Nmax),IC3(Nmax,3)
      DIMENSION N6MEM(Mmax,6)
      DIMENSION nsw(4,66),IP(2,66),KSW(4),IBSW(66,2),JBSW(66,2),IAtSW(4)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
      Write(Iout,1000)
      Do J=1,66
      Do I=1,2
       IP(I,J)=0
      enddo
      enddo

C     read pentagon numbers
      Read(IN,*,Err=100,end=100) ((IP(I,J),I=1,2),J=1,66)

C     Check if any is a Stone-Wales pattern
  100 ntrans=0
      Do J=1,66
       if(IP(1,J).eq.0) go to 10
       ntrans=ntrans+1
      enddo
      if(ntrans.lt.1) then
       Write(Iout,1016) 
       return
      endif
   10 Do J=1,ntrans
       if(IP(2,J).eq.0) then
        J1=IP(1,J)
        IP(1,J)=nsw(1,J1)
        IP(2,J)=nsw(4,J1)
       endif
      enddo
      Write(Iout,1001) ntrans,((IP(I,J),I=1,2),J=1,ntrans)
      nswt=0
      Do J=1,ntrans
       I1=IP(1,J)
       I2=IP(2,J)
       If(I1.gt.I2) then
        I3=I1
        I1=I2
        I2=I3
       endif
       Do K=1,numbersw
        J1=nsw(1,K)
        J2=nsw(4,K)
        if(I1.eq.J1.and.I2.eq.J2) then
         nswt=nswt+1
         if(nswt.ne.k) then
          do L=1,4
          KSW(L)=nsw(L,K)
          nsw(L,k)=nsw(L,nswt)
          nsw(L,nswt)=KSW(L)
          enddo
         endif
        endif
       enddo
      enddo
       If(nswt.ne.ntrans) then
        Write(Iout,1002) nswt,ntrans
        Write(Iout,1003) ((nsw(I,J),I=1,4),J=1,nswt)
        Write(Iout,1004) 
        Return
       endif
        Write(Iout,1007) nswt
        Write(Iout,1003) ((nsw(I,J),I=1,4),J=1,nswt)

C Perform Stone-Wales transformation
C Find bonds between 6-rings for 90 degree rotation
      Do I=1,nswt
       IR1=nsw(2,I)-12
       IR2=nsw(3,I)-12
        ibond=0
        Do J1=1,6
         IAtom1=N6MEM(IR1,J1)
        Do J2=1,6
         IAtom2=N6MEM(IR2,J2)
         if(IAtom1.eq.IAtom2) then
          ibond=ibond+1
          IBSW(I,ibond)=IAtom1
         endif
        enddo
        enddo
        If(IBSW(I,1).gt.IBSW(I,2)) then
         IX=IBSW(I,1)
         IBSW(I,1)=IBSW(I,2)
         IBSW(I,2)=IX
        endif
      enddo
      Write(Iout,1012) ((IBSW(I,J),J=1,2),I=1,nswt)

C Get atoms connected to these for bond rotation
      Write(Iout,1013) 
      Do I=1,nswt
       icount=0
        IC1=IBSW(I,1)
        IC2=IBSW(I,2)
       do J=1,3
        if(IC3(IC1,J).ne.IC2) then
         icount=icount+1
         IAtSW(icount)=IC3(IC1,J)
        endif
       enddo
       do J=1,3
        if(IC3(IC2,J).ne.IC1) then
         icount=icount+1
         IAtSW(icount)=IC3(IC2,J)
        endif
       enddo
       ifirst=IAtSW(1)
       ihit=0
       isix1=nsw(2,I)
       isix2=nsw(3,I)
      do J=1,6
       if(ifirst.eq.N6MEM(isix1,J)) ihit=1
      enddo
       if(ihit.eq.0) then
        isix1=nsw(3,I)
        isix2=nsw(2,I)
       endif
       isecond=IAtSW(3)
       ihit=0
      do J=1,6
       if(isecond.eq.N6MEM(isix2,J)) ihit=1
      enddo
       if(ihit.eq.0) then
        isecond=IAtSW(4)
       endif
       JBSW(I,1)=ifirst
       JBSW(I,2)=isecond
      Write(Iout,1014) JBSW(I,1),IBSW(I,1),IBSW(I,2),JBSW(I,2)
      enddo

C Transform adjacency matrix
      Write(Iout,1008)
      do I=1,nswt
       IE1=JBSW(I,1)
       IE2=IBSW(I,1)
       JE1=IBSW(I,2)
       JE2=JBSW(I,2)
       IDA(IE1,IE2)=0
       IDA(IE2,IE1)=0
       IDA(JE1,JE2)=0
       IDA(JE2,JE1)=0
       Write(Iout,1009) IE1,IE2,JE1,JE2
      enddo 
      Write(Iout,1015)
      do I=1,nswt
       IE1=JBSW(I,1)
       IE2=IBSW(I,2)
       JE1=IBSW(I,1)
       JE2=JBSW(I,2)
       IDA(IE1,IE2)=1
       IDA(IE2,IE1)=1
       IDA(JE1,JE2)=1
       IDA(JE2,JE1)=1
       Write(Iout,1009) IE1,IE2,JE1,JE2
      enddo 

C Adjacency matrix constructed
C Now analyze the adjacency matrix if it is correct
      nsum=0
      Do I=1,MAtom
      isum=0
      Do J=1,MAtom
      isum=isum+IDA(I,J)
      enddo
      If(isum.ne.3) nsum=nsum+1
      enddo
      if(nsum.ne.0) then
      WRITE(Iout,1005) nsum,isum
      stop
      else
      WRITE(Iout,1006)
      endif
      Call Tutte(Matom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
 1000 Format(/1X,'Stone-Wales transformation:',
     1 /1X,'Read pentagon ring numbers (between 1-12)')
 1001 Format(/1X,'Number of Stone-Wales transformations: ',I2,
     1 /1X,'Stone-Wales pentagon ring numbers:',/,
     1 10(' (',I2,','I2,') '))
 1002 Format(/1X,I2,' Number of pentagon rings for Stone-Wales ',
     1 'transformation does not match list, only ',I2,' found:')
 1003 Format(7(' (',I2,',',I5,',',I5,',',I2,')'))
 1004 Format(/1X,'==> RETURN')
 1005 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1006 FORMAT(1X,'Graph checked, it is cubic')
 1007 Format(1X,'Perform ',I2,' Stone-Wales transformations ',
     1 /1X,'Modifying adjacency matrix for rings (P,H,H,P):')
 1008 Format(1X,'Transform adjacency matrix:',/1X,
     1 'Delete edges:')
 1009 Format(1X,' (',I5,'-',I5,') (',I5,'-',I5,')')
 1012 Format(1X,'Stone-Wales bonds to be rotated by 90 degrees. ',
     1 'Edge numbers:',/,10(' (',I5,','I5,') '))
 1013 Format(1X,'Vertices affected in Stone-Wales transformation:')
 1014 Format(1X,' (',I5,',',I5,',',I5,',',I5,')')
 1015 Format(1X,'Add edges:')
 1016 Format(/1X,'No input found ==> RETURN')
      Return
      END

      SUBROUTINE GoldbergCoxeter(Matom,Iout,leap,leapGC,kGC,lGC,
     1 ihueckel,LeapErr,IDA,A,evec,df,Dist,layout2D,distp,CDist) 
C     Construct Leapfrog fullerene through adjacency matrix
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax),IDA(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
      integer graph_is_a_fullerene
      type(c_ptr) :: g, frog, halma, new_fullerene_graph,
     1 leapfrog_fullerene, halma_fullerene
      LeapErr=0 

C Leapfrog fullerene
  10  if(leap.gt.0) then
      MLeap=(3**leap)*MAtom

      if(Mleap.gt.Nmax) then
         Write(Iout,1002) MLeap,Nmax
         LeapErr=1
         return
      endif

      if(leap.eq.1) then
       Write(Iout,1000) MAtom,MLeap
      else
       Write(Iout,1001) leap,MAtom,MLeap
      endif 

      g = new_fullerene_graph(Nmax,MAtom,IDA)
      frog = leapfrog_fullerene(g,leap)

C Test that the created leapfrog graph is a fullerene graph
C (cubic, satisfies Eulers formula, has 12 pentagons, remaining 
C  faces are hexagons)
      isafullerene = graph_is_a_fullerene(frog)

      if(isafullerene .eq. 1) then
         write (Iout,1014) 
      else
         write (Iout,1015) 
      endif
C Produce adjacency matrix
      write (Iout,1005) 
      call adjacency_matrix(frog,Nmax,IDA)
      endif

C Goldberg-Coxeter transform of initial fullerene
C Input: initial graph, and GC indices (kGC,lGC) 
      if(leapGC.gt.0) then
      Write(Iout,1010) kGC,lGC,kGC,lGC
      if(kGC.eq.1.and.lGC.eq.0) then
       Write(Iout,1006) kGC,lGC
       return
      endif
      if(kgc.eq.1.and.lGC.eq.1) then
        write(Iout,1008)
        leapGC=0
        leap=1
        go to 10       
      endif
      if(lGC .ne. 0) then
        write(Iout,1011)
        stop
      endif
      g = new_fullerene_graph(Nmax,MAtom,IDA)
      halma = halma_fullerene(g,kGC-1)
      isafullerene = graph_is_a_fullerene(halma)
      IF (isafullerene.eq.1) then
        write (iout,1013)
      else
        write (iout,1009)
        stop
      endif
C Update fortran structures
      MLeap  = NVertices(halma)
      Medges = NEdges(halma)
        write(Iout,1012)  MLeap,Medges
C Produce adjacency matrix 
      write (Iout,1005) 
      call adjacency_matrix(halma,Nmax,IDA)
      write (Iout,1007) 
      endif
      MAtom = MLeap

      Call Tutte(MAtom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
      if(leap.ne.0) call delete_fullerene_graph(frog)
      if(leapGC.ne.0) call delete_fullerene_graph(halma)
      call delete_fullerene_graph(g)
      write (Iout,1004) 
 1000 Format(/1X,'Creating the adjacency matrix of the next leap-frog',
     1 ' fullerene: ',I4,' --> ',I4)
 1001 Format(/1X,'Creating the adjacency matrix of the ',I2,
     1 'th leap-frog fullerene: ',I4,' --> ',I4)
 1002 Format(1X,'Error: Dimension of leapfrof fullerene is ',I4,
     1 ' greater than dimension of Nmax (',I4,') set in program')
 1004 FORMAT(1X,'Fullerene graph deleted')
 1005 Format(1x,'Produce new adjacency matrix')
 1006 Format(/1x,'Goldberg-Coxeter transformation with indices ',
     1 '(k,l) = (',I2,',',I2,') is the identity ==> Return')
 1007 Format(1x,'Adjacency matrix produced')
 1008 Format(1x,'Goldberg-Coxeter transformation is of leapfrog type')
 1009 Format(1x,'Halma fullerene is not a fullerene')
 1010 Format(/1x,'Goldberg-Coxeter transformation with indices ',
     1 '(k,l) = (',I2,',',I2,') of initial fullerene: GC(',I2,',',I2,
     1 ')[G0]')
 1011 Format(/1x,'Goldberg-Coxeter construction not implemented',
     1 ' for l > 0.')
 1012 Format(1x,'Updating number of vertices (',I5,') and edges (',
     1 I5,')')
 1013 Format(1x,'Halma fullerene is a fullerene')
 1014 Format(1X,'Leapfrog graph satisfies all fullerene conditions')
 1015 Format(1X,'Leapfrog graph does not satisfy all fullerene ',
     1 'conditions')
      Return
      END

      SUBROUTINE Tutte(Matom,Iout,ihueckel,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist)
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D(2,Nmax)
      DIMENSION evec(Nmax),df(Nmax),A(Nmax,Nmax)
      DIMENSION Dist(3,Nmax),distP(Nmax)
      DIMENSION IDA(Nmax,Nmax)
      type(c_ptr) :: g, new_fullerene_graph

C Produce Hueckel matrix and diagonalize
C     Diagonalize
      if(ihueckel.ne.0) then
      Do I=1,MAtom
         Do J=1,MAtom
            A(I,J)=dfloat(IDA(I,J))
         enddo
      enddo
      call tred2(A,MAtom,Nmax,evec,df)
      call tqli(evec,df,MAtom,Nmax,A)
      Write(Iout,1005) MAtom,MAtom
C     Sort eigenvalues evec(i) and eigenvectors A(*,i)
      Do I=1,MAtom
         e0=evec(I)
         jmax=I
      Do J=I+1,MAtom
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
      Do k=1,MAtom
      df(k)=A(k,jmax)
      A(k,jmax)=A(k,I)
      A(k,I)=df(k)
      enddo
      endif
      enddo
C Analyze eigenenergies
      Call HueckelAnalyze(MAtom,NMax,Iout,iocc,df,evec)
      endif

C   Tutte algorithm for the 3D structure (see pentindex.f):
         write (Iout,1011)
         g = new_fullerene_graph(Nmax,MAtom,IDA)
         call tutte_layout(g,layout2d)
         call spherical_layout(g,Dist)
         write (Iout,1013)
         call delete_fullerene_graph(g)
         write (Iout,1014)

C     Obtain smallest distance for further scaling
C     Now this contracts or expands the whole fullerene to set the
C     smallest bond distance to Cdist
c the same functionality is in pentindex.f, twice.
c     corraction: setting the shortest bond to to cdist is not a good idea.  It is beneficial to set the avarage bond length to some value, like e.g. 4*cdist
      R0=1.d10
      Rsum=0.d0
      Do I=1,Matom
        Do J=I+1,Matom
          if (IDA(i,j) .ne. 0) then
            X=Dist(1,I)-Dist(1,J)
            Y=Dist(2,I)-Dist(2,J)
            Z=Dist(3,I)-Dist(3,J)
            R=dsqrt(X*X+Y*Y+Z*Z)
            rsum=rsum+r
c            if(R.lt.R0) R0=R
          endif
        enddo
      enddo
      rsum=rsum/(3*matom/2)
c      fac=CDist/R0
      fac=4.0*CDist/Rsum
      Do I=1,Matom
      Dist(1,I)=Dist(1,I)*fac
      Dist(2,I)=Dist(2,I)*fac
      Dist(3,I)=Dist(3,I)*fac
c      write(*,*)dist(1,i),dist(2,i),dist(3,i),
c     1  dsqrt(dist(1,i)**2+dist(2,i)**2+dist(3,i)**2)
      enddo
C     Check distances
      Write(IOUT,1015) fac
      Do J=1,MAtom
      Write(IOUT,1016) J,(Dist(I,J),I=1,3)
      enddo
      CALL Distan(MAtom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1017) Rmin,Rmax,rms
      ratio=(Rmax/Rmin-1.d0)*1.d2
      iratio=dint(ratio)
      CALL Diameter(MAtom,IOUT,Dist,distp)
      if(iratio.lt.33) then
      Write(IOUT,1018) iratio
      else
      Write(IOUT,1019) iratio
      endif
 1005 FORMAT(/1X,'Construct the (',I4,','I4,') Hueckel ',
     1 ' matrix, diagonalize (E=alpha+x*beta) and get eigenvectors',
     1 /1X,'Eigenvalues are between [-3,+3]')
 1011 FORMAT(/1X,'Using the Tutte embedding algorithm to construct ',
     1 'the fullerene')
 1013 FORMAT(1X,'Projected on sphere')
 1014 FORMAT(1X,'Fullerene graph deleted')
 1015 FORMAT(1X,'Coordinates from Tutte embedding scaled by a factor'
     1 ' of ',D18.12)
 1016 FORMAT(1X,I4,5X,3(D18.12,2X))
 1017 FORMAT(1X,'Minimum distance: ',F12.6,', Maximum distance: ',F12.6,
     1 ', RMS distance: ',F12.6)
 1018 Format(1X,'Maximum bond distance ',I5,'% larger than minimum ',
     1 'distance')
 1019 Format(1X,'Maximum bond distance ',I5,'% larger than minimum ',
     1 'distance. Fullerene strongly distorted!',/1X,
     1 'Optimization of geometry recommended')
      Return
      END

      Subroutine Permute(Matom,Iout,Dist,IC3)
      use config
      implicit double precision (a-h,o-z)
      dimension Dist(3,Nmax),IC3(Nmax,3),Iperm(Nmax,2)
C Sort Cartesian coordinated so atom number i is connected to
C  atoms 1,...,I-1

      nperm=0
      do 10 i=2,MAtom
      do j=i,MAtom-1
       do k=1,3
        if(IC3(j,k).le.i-1) then
         if(i.ne.j) then
C    change distances
          x=Dist(1,i)
          y=Dist(2,i)
          z=Dist(3,i)
          Dist(1,i)=Dist(1,j)
          Dist(2,i)=Dist(2,j)
          Dist(3,i)=Dist(3,j)
          Dist(1,j)=x
          Dist(2,j)=y
          Dist(3,j)=z
C    change IC3
         iz1=IC3(i,1)
         iz2=IC3(i,2)
         iz3=IC3(i,3)
         IC3(i,1)=IC3(j,1)
         IC3(i,2)=IC3(j,2)
         IC3(i,3)=IC3(j,3)
         IC3(j,1)=iz1
         IC3(j,2)=iz2
         IC3(j,3)=iz3
         do k1=1,MAtom
         do 11 k2=1,3
          if(IC3(k1,k2).eq.i) then
           IC3(k1,k2)=j
           go to 11
          endif
          if(IC3(k1,k2).eq.j) then
           IC3(k1,k2)=i
          endif
   11    continue
         enddo

C    record change
          nperm=nperm+1
          Iperm(nperm,1)=i
          Iperm(nperm,2)=j
         endif
        go to 10
        endif
       enddo

      enddo
  10  continue
      Write(Iout,1000) nperm
      Write(Iout,1001) (Iperm(i,1),Iperm(i,2),i=1,nperm)

 1000 Format(1X,'Number of permutations for cartesian coordinates ',
     1  'performed: ',I5,', Permutations:')
 1001 Format(10(1X,'(',I4,',',I4,')'))
      Return
      END

      Subroutine CartInt(Dist,MAtom,Iout,isort)
      use config
      implicit double precision (a-h,o-z)
      dimension Dist(3,Nmax),na(Nmax),nb(Nmax),nc(Nmax),zmatrix(3,Nmax)
C Modified routine xyzint
C Cartint works out the internal coordinates (z-matrix) of a molecule
C        atoms N1,N2,N3,N4 defined with distances, bond angles and dihedrals
C        angles in the range 15 to 165 degrees if possible.
C   on input  Dist    = cartesian array of MAtom atoms
C   on output ZMatrix = Z-Matrix
C        MAtom: number of atoms
C        degree = 360/2Pi = 57.29578..., angles are in degrees
C
      degree=1.8d2/dpi
 
      Write(Iout,1000) 
       do i=1,MAtom
        zmatrix(1,i)=0.d0
        zmatrix(2,i)=0.d0
        zmatrix(3,i)=0.d0
        na(i)=0
        nb(i)=0
        nc(i)=0
       enddo

       do 30 i=1,MAtom
        na(i)=2
        nb(i)=3
        nc(i)=4
        im1=i-1
        if(im1.eq.0) go to 30
        sum=1.d30
        do j=1,im1
          r=(Dist(1,i)-Dist(1,j))**2+
     1     (Dist(2,i)-Dist(2,j))**2+
     2     (Dist(3,i)-Dist(3,j))**2
          if(r.lt.sum.and.na(j).ne.j.and.nb(j).ne.j) then
            sum=r
            k=j
          endif
        enddo
c   atom i is nearest to atom k
      na(i)=k
      if(i.gt.2)nb(i)=na(k)
      if(i.gt.3)nc(i)=nb(k)
c   find any atom to relate to na(i)
   30    continue

      na(1)=0
      nb(1)=0
      nc(1)=0
      nb(2)=0
      nc(2)=0
      nc(3)=0

c   na, nb, nc are determined, now get zmatrix
      call Distgeo(Dist,MAtom,na,nb,nc,degree,zmatrix)

c     print
       do i=1,MAtom
        Write(Iout,1001) i,(zmatrix(J,I),J=1,3),na(i),nb(i),nc(i)
       enddo

c   Check if distances are within certain range
       rmindist=zmatrix(1,2)
       rmaxdist=zmatrix(1,2)
       do i=3,MAtom
        if(zmatrix(1,i).lt.rmindist) rmindist=zmatrix(1,i)
        if(zmatrix(1,i).gt.rmaxdist) rmaxdist=zmatrix(1,i)
       enddo
      ratio=rmaxdist/rmindist
      if(ratio.lt.1.5d0) then
       Write(Iout,1002) rmindist,rmaxdist
       isort=0
      else
       Write(Iout,1003) rmindist,rmaxdist
       isort=1
      endif
      if(rmindist.lt.1.d0) then
       Write(Iout,1004)
      endif
      
 1000 Format(/1X,'Convert cartesian into internal coordinates:',
     1 /4X,'N1',9X,'R',10X,'Angle',6X,'Dihedral',5X,'N2',4X,
     1 'N3',4X,'N4',/1X,64('-'))
 1001 Format(1X,I5,3(1X,F12.6),3(1X,I5))
 1002 Format(1X,'Analysis of distances: All are bond distances',
     1 /1X,'Smallest distance: ',F12.6,', Largest distance: ',F12.6)
 1003 Format(1X,'Analysis of distances: Requires sorting of cartesian',
     1 ' coordinates after getting connectivities',
     1 /1X,'Smallest distance: ',F12.6,', Largest distance: ',F12.6)
 1004 Format('WARNING: Problem with smallest bond distance detected!',
     1 /1X,' Geometry of fullerene may be wrong')
      return
      end

      subroutine Distgeo(Dist,MAtom,na,nb,nc,degree,zmatrix)
      use config
      implicit double precision (a-h,o-z)
      dimension Dist(3,Nmax),na(Nmax),nb(Nmax),nc(Nmax),zmatrix(3,Nmax)
C Distgeo converts coordinates from cartesian to internal.
C  input Dist  = array of cartesian coordinates
C   MAtom= number of atoms
C   na   = numbers of atom to which atoms are related by distance
C   nb   = numbers of atom to which atoms are related by angle
C   nc   = numbers of atom to which atoms are related by dihedral
C  output zmatrix  = internal coordinates in angstroms, radians, and radians

      Data tol/0.2617994d0/

      do 30 i=2,MAtom
         j=na(i)
         k=nb(i)
         l=nc(i)
         if(i.lt.3) go to 30
         ii=i
         call angle(Dist(1,ii),Dist(2,ii),Dist(3,ii),Dist(1,j),
     1           Dist(2,j),Dist(3,j),Dist(1,k),Dist(2,k),Dist(3,k),
     1           zmatrix(2,i))
         zmatrix(2,i)=zmatrix(2,i)*degree
         if(i.lt.4) go to 30
c   make sure dihedral is meaningful
         call angle(Dist(1,j),Dist(2,j),Dist(3,j),Dist(1,k),
     1           Dist(2,k),Dist(3,k),Dist(1,l),Dist(2,l),Dist(3,l),angl)
         if(angl.gt.dpi-tol.or.angl.lt.tol)then
c  angle is unsatisfactory, let's search for another atom for
c  defining the dihedral.
   10       sum=1.d2
            do i1=1,ii-1
               r=(Dist(1,i1)-Dist(1,k))**2+
     1          (Dist(2,i1)-Dist(2,k))**2+
     2          (Dist(3,i1)-Dist(3,k))**2
               if(r.lt.sum.and.i1.ne.j.and.i1.ne.k) then
         call angle(Dist(1,j),Dist(2,j),Dist(3,j),Dist(1,k),
     1           Dist(2,k),Dist(3,k),Dist(1,i1),Dist(2,i1),Dist(3,i1),
     1           angl)
                  if(angl.lt.dpi-tol.and.angl.gt.tol)then
                     sum=r
                     l=i1
                     nc(ii)=l
                  endif
               endif
            enddo
            if(sum.gt.99.d0.and.tol.gt.0.1d0)then
c
c anything within 5 degrees?
c
               tol=0.087266d0
               go to 10
            endif
         endif
         call dihedral(Dist(1,ii),Dist(2,ii),Dist(3,ii),Dist(1,j),
     1           Dist(2,j),Dist(3,j),Dist(1,k),Dist(2,k),Dist(3,k),
     1           Dist(1,l),Dist(2,l),Dist(3,l),zmatrix(3,i))
         zmatrix(3,i)=zmatrix(3,i)*degree
   30 zmatrix(1,i)=dsqrt((Dist(1,i)-Dist(1,j))**2+
     1           (Dist(2,i)-Dist(2,j))**2+
     2           (Dist(3,i)-Dist(3,j))**2)
      zmatrix(1,1)=0.d0
      zmatrix(2,1)=0.d0
      zmatrix(3,1)=0.d0
      zmatrix(2,2)=0.d0
      zmatrix(3,2)=0.d0
      zmatrix(3,3)=0.d0
      return
      end

