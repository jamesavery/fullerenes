      SUBROUTINE MoveCM(Natom,Matom,Iout,Iprint,IAtom,Dist,DistCM,El)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,natom),DistCM(3),Ainert(3,3),evec(3),df(3)
      DIMENSION IATOM(natom),layout2d(2,NAtom)
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
 1002 FORMAT(1X,I4,1X,I6,1X,A2,6X,3(D18.12,2X))
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

      Subroutine Leapfrog(NAtom,MAtom,Iout,leap,LeapErr,IDA,
     1 A,evec,df,Dist,layout2D,distp,CDist) 
C     Construct Leapfrog fullerene through adjacency matrix
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 layout2D
      DIMENSION evec(NAtom),df(NAtom),A(NAtom,NAtom),IDA(NAtom,NAtom)
      DIMENSION IDG(NAtom),Dist(3,NAtom),distP(NAtom),layout2D(2,NAtom)
      Character*10 Symbol
      Data Tol1/.15d0/
      integer graph_is_a_fullerene
      type(c_ptr) :: g, frog, new_fullerene_graph, leapfrog_fullerene,
     1 new_graph
C     Print and test 
      LeapErr=0 

      MLeap=(3**leap)*MAtom

      if(Mleap.gt.NAtom) then
         Write(Iout,1002) MLeap,NAtom
         LeapErr=1
         return
      endif

      if(leap.eq.1) then
       Write(Iout,1000) MAtom,MLeap
      else
       Write(Iout,1001) leap,MAtom,MLeap
      endif 

      g = new_fullerene_graph(NAtom,MAtom,IDA)
      frog = leapfrog_fullerene(g,leap)

C Test that the created leapfrog graph is a fullerene graph
C (cubic, satisfies Eulers formula, has 12 pentagons, remaining 
C  faces are hexagons)
      isafullerene = graph_is_a_fullerene(frog)

      if(isafullerene .eq. 1) then
         write (Iout,1020) 
      else
         write (Iout,1021) 
      endif

      call adjacency_matrix(frog,NAtom,IDA)

C Produce Hueckel matrix and diagonalize
C     Diagonalize
      Do I=1,MLeap
         Do J=1,MLeap
            A(I,J)=dfloat(IDA(I,J))
         enddo
      enddo
      call tred2(A,MLeap,NAtom,evec,df)
      call tqli(evec,df,MLeap,NAtom,A)
      Write(Iout,1005) MLeap,MLeap
C     Sort eigenvalues evec(i) and eigenvectors A(*,i)
      Do I=1,MLeap
         e0=evec(I)
         jmax=I
      Do J=I+1,MLeap
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
      Do k=1,MLeap
      df(k)=A(k,jmax)
      A(k,jmax)=A(k,I)
      A(k,I)=df(k)
      enddo
      endif
      enddo
C     Now sort degeneracies
      df(1)=evec(1)
      ieigv=1
      ideg=1
      IDG(1)=ideg
      Do I=2,MLeap
      diff=dabs(evec(I-1)-evec(I))
      if(diff.lt.Tol) then
      ideg=ideg+1
      IDG(ieigv)=ideg
      else
      ieigv=ieigv+1
      ideg=1
      IDG(ieigv)=ideg
      df(ieigv)=evec(I)
      endif
      enddo

C     Now Print
      ntot=0
      nopen=0
      nflag=0
      iocc=0
      Write(Iout,1006)
      Do I=1,ieigv
      NE=2*idg(i)
      NE1=NE
      ntot=ntot+NE
      Symbol='(occupied)'
      if(ntot.gt.MLeap) then
      if(nflag.eq.0) then
      nflag=1
      bandgap=df(i-1)-df(i)
      endif
      NE=0
      Symbol='(empty)   '
      endif
      if(ntot.gt.MLeap.and.(ntot-NE1).lt.MLeap) then
      NE=MLeap-ntot+NE1
      Symbol='(fractocc)'
      nopen=1
      endif
      if(NE.ne.0.and.NE.eq.idg(i)*2) iocc=iocc+idg(i)
      Write(Iout,1007) df(I),idg(i),NE,Symbol
      enddo
      Write(Iout,1008)
      if(nopen.eq.1) then
      Write(Iout,1009)
      else
      Write(Iout,1010) bandgap
      if(bandgap.lt.Tol1) Write(Iout,1009)
      endif

C   Tutte algorithm for the 3D structure (see pentindex.f):
         write (Iout,1012)
         call tutte_layout(frog,layout2d)
         call spherical_layout(frog,Dist)
         write (Iout,1013)

C     Obtain smallest distance for further scaling
C     Now this contracts or expands the whole fullerene to set the
C     smallest bond distance to Cdist
      R0=1.d10
      Do I=1,MLEAP
      Do J=I+1,MLEAP
      X=Dist(1,I)-Dist(1,J)
      Y=Dist(2,I)-Dist(2,J)
      Z=Dist(3,I)-Dist(3,J)
      R=dsqrt(X*X+Y*Y+Z*Z)
      if(R.lt.R0) R0=R
      enddo
      enddo
      fac=CDist/R0
      Do I=1,MLEAP
      Dist(1,I)=Dist(1,I)*fac
      Dist(2,I)=Dist(2,I)*fac
      Dist(3,I)=Dist(3,I)*fac
      enddo
C     Check distances
      Write(IOUT,1015) fac
      Do J=1,MLeap
      Write(IOUT,1016) J,(Dist(I,J),I=1,3)
      enddo
      CALL Distan(NAtom,MLeap,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1017) Rmin,Rmax,rms
      ratio=(Rmax/Rmin-1.d0)*1.d2
      iratio=dint(ratio)
      CALL Diameter(NAtom,MLeap,IOUT,Dist,distp)
      if(iratio.lt.33) then
      Write(IOUT,1018) iratio
      else
      Write(IOUT,1019) iratio
      endif

      MAtom = MLeap
      call adjacency_matrix(frog,NAtom,IDA)

      call delete_fullerene_graph(frog)
      call delete_fullerene_graph(g)

 1000 Format(/1X,'Creating the adjacency matrix of the next leap-frog',
     1 ' fullerene: ',I4,' --> ',I4)
 1001 Format(/1X,'Creating the adjacency matrix of the ',I2,
     1 'th leap-frog fullerene: ',I4,' --> ',I4)
 1002 Format(1X,'Error: Dimension of leapfrof fullerene is ',I4,
     1 ' greater than dimension of NAtom (',I4,') set in program')
 1003 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1004 FORMAT(1X,'Graph checked, it is cubic')
 1005 FORMAT(/1X,'Using the Tutte embedding algorithm to construct ',
     1 'the fullerene',/1X,'Construct the (',I4,','I4,') Hueckel ',
     1 ' matrix, diagonalize (E=alpha+x*beta) and get eigenvectors',
     1 /1X,'Eigenvalues are between [-3,+3]')
 1006 FORMAT(1X,'       x     deg NE   type    ',/1X,32('-'))
 1007 FORMAT(1X,F12.6,I4,1X,I4,3X,A10)
 1008 FORMAT(1X,32('-'))
 1009 FORMAT(1X,'Fullerene has open-shell character (zero band gap)!')
 1010 FORMAT(1X,'Bandgap delta x = ',F12.6,' (in units of |beta|)')
 1011 FORMAT(/1X,'Using the Tutte-embedding algorithm to construct ',
     1 'the fullerene',/1X,'Construct the Tutte planar graph and ',
     1 'project on sphere')
 1012 Format(1X,'Calculating Tutte-embedding')
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
 1020 Format(1X,'Leapfrog graph satisfies all fullerene conditions')
 1021 Format(1X,'Leapfrog graph does not satisfy all fullerene ',
     1 'conditions')

      Return
      END

      SUBROUTINE GoldbergCoxeter(NMAX,MMAX,LMAX,MAtom,kGC,lGC,IN,Iout,
     1 IDA,A,evec,df,Dist,layout2d,distp,Cdist)
C 
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 layout2d
      DIMENSION layout2d(2,NMAX)
      DIMENSION Dist(3,NMAX),distP(NMAX)
      DIMENSION A(NMAX,NMAX),IDA(NMAX,NMAX)
      DIMENSION evec(NMAX),df(NMAX),dipol(3,3)
      type(c_ptr) :: g, halma, halma_fullerene, new_fullerene_graph, 
     1               new_graph
C Construct Goldberg Coxeter fullerene using the two indices (kGC,lGC)

      Write(Iout,1000) kGC,lGC

      if(lGC .ne. 0) then
        write(Iout,1010)
        stop
      endif

      g = new_graph(NMAX,MAtom,IDA)
      call print_graph("g",g)
      stop

      g = new_fullerene_graph(NMAX,MAtom,IDA)

      call print_graph("g",g)
      halma = halma_fullerene(g,kGC-1)
      call print_graph("halma",halma)

C      call adjacency_matrix(halma,NMAX,IDA)
      Print*,IDA
      stop

 1000 Format(/1x,'Halma fullerene with Coxeter indices (k,l) = (',I2,','
     1 ,I2,')')
 1010 Format(/1x,'Goldberg-Coxeter construction not implemented',
     1 ' for l > 0.')
      Return 
      END
