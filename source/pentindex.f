      SUBROUTINE CoordPent(NMAX,MMAX,LMAX,MAtom,IN,Iout,IDA,D,ICart,
     1 IV1,IV2,IV3,A,evec,df,Dist,layout2d,distp,Cdist,GROUP)
C Cartesian coordinates produced from ring spiral pentagon list
C using either the Fowler-Manopoulus or the Tutte embedding
C algorithm 
C Fowler-Manopoulus algorithm: identify P-type eigenvectors and 
C construct the 3D fullerene
C Tutte embedding algorithm: Tutte embedding and sphere projection
C If nalgorithm=0 use Fowler-Manopoulus algorithm
C If nalgorithm=1 use Fowler-Manopoulus algorithm but Laplacian instead
C If nalgorithm=2 use Tutte algorithm
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 layout2d
      Integer D,S,RT,Spiral
      DIMENSION layout2d(2,NMAX)
      DIMENSION D(MMAX,MMAX),S(MMAX),Dist(3,NMAX),distP(NMAX)
      DIMENSION NMR(6),JP(12),A(NMAX,NMAX),IDA(NMAX,NMAX)
      DIMENSION evec(NMAX),df(NMAX),dipol(3,3)
      DIMENSION IDG(NMax)
      DIMENSION Spiral(12,NMAX)
      Character*10 Symbol
      CHARACTER*3 GROUP
      Data Tol,Tol1,Tol2,ftol/1.d-5,.15d0,1.5d1,1.d-10/
      type(c_ptr) :: g, new_fullerene_graph, new_graph

      nalgorithm=ICart-2
      M=Matom/2+2
      istop=0
      Group='NA '
C Read pentagon list and produce adjacency matrix
      Read(IN,*) (JP(I),I=1,12)
C     Produce the Spiral S using the program WINDUP and UNWIND
      do I=1,MMAX
      do J=1,MMAX
      D(I,J)=0
      enddo
      enddo
      Do I=1,6
      NMR(I)=0
      enddo
      Do I=1,M
      S(I)=6
      enddo
C     Search where the 5-rings are in the spiral
      Do I=1,12
      S(JP(I))=5
      enddo
      IPR=0
      IER=0
      CALL Windup(MMAX,M,IPR,IER,S,D)              ! Wind up spiral into dual
      IF(IER.gt.0) then
      WRITE(Iout,1000) IER
      endif
      IT=1
      Do I=1,12
      Spiral(I,1)=JP(I)
      enddo
      CALL Unwind(NMAX,MMAX,LMAX,M,IER,IT,ispiral,
     1 Spiral,S,D,NMR,Group)                       ! Unwind dual into spirals
      K=0
      DO J=1,6
         IF(NMR(J).EQ.0) GO TO 3
         K=J
      enddo
  3   If(K.le.0) then
      WRITE(Iout,1020) M,Matom,GROUP,(JP(I),I=1,12)
      else
      WRITE(Iout,1001) M,Matom,GROUP,(JP(I),I=1,12),(NMR(J),J=1,K)
      endif
      if(ispiral.ge.2) then
       if(ispiral.eq.2) then
       WRITE(Iout,1023)
       Do II=1,12
       JP(II)=spiral(II,2)
       enddo 
       else
       WRITE(Iout,1019) ispiral-1
       endif
      Do JJ=2,ispiral 
      WRITE(Iout,1021) (spiral(II,JJ),II=1,12)
      enddo
      else
      WRITE(Iout,1022)
      endif
      if(ispiral.gt.2) then
      CALL CanSpiral(NMAX,ispiral,spiral,JP)
      WRITE(Iout,1023)
      WRITE(Iout,1021) (JP(I),I=1,12)
      endif
      Do I=1,M
      S(I)=6
      enddo
      Do I=1,12
      S(JP(I))=5
      enddo
      WRITE(Iout,1024)
      WRITE(Iout,1025) (S(I),I=1,M)
C End of Spiral Program, dual matrix in D(i,j)

C Now produce the Hueckel matrix from the dual matrix
      CALL DUAL(NMAX,D,MMAX,IDA,Matom,IER)
      IF(IER.ne.0) then
      WRITE(Iout,1002) IER
      stop
      endif
      Do I=1,MAtom
      Do J=1,MAtom
      A(I,J)=dfloat(IDA(I,J))
      enddo
      enddo
C Analyze the adjacency matrix if it is correct
      nsum=0
      Do I=1,MAtom
      isum=0
      Do J=1,MAtom
      isum=isum+IDA(I,J)
      enddo
      If(isum.ne.3) nsum=nsum+1
      enddo
      if(nsum.ne.0) then
      WRITE(Iout,1037) nsum,isum
      stop
      else
      WRITE(Iout,1038)
      endif

C Produce Hueckel matrix and diagonalize
C     Diagonalize
      call tred2(A,Matom,NMax,evec,df)
      call tqli(evec,df,Matom,NMax,A)
      Write(Iout,1004) Matom,Matom
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
C     Now sort degeneracies
      df(1)=evec(1)
      ieigv=1
      ideg=1
      IDG(1)=ideg
      Do I=2,MAtom
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
      Write(Iout,1003) 
      Do I=1,ieigv
      NE=2*idg(i)
      NE1=NE
      ntot=ntot+NE
      Symbol='(occupied)'
      if(ntot.gt.Matom) then 
      if(nflag.eq.0) then
      nflag=1
      bandgap=df(i-1)-df(i)
      endif
      NE=0
      Symbol='(empty)   '
      endif
      if(ntot.gt.Matom.and.(ntot-NE1).lt.Matom) then 
      NE=Matom-ntot+NE1
      Symbol='(fractocc)'
      nopen=1
      endif
      if(NE.ne.0.and.NE.eq.idg(i)*2) iocc=iocc+idg(i)
      Write(Iout,1005) df(I),idg(i),NE,Symbol
      enddo
      Write(Iout,1006)
      if(nopen.eq.1) then
      Write(Iout,1007)
      else
      Write(Iout,1008) bandgap
      if(bandgap.lt.Tol1) Write(Iout,1009)
      endif
C     Laplacian matrix is used instead of the adjaceny matrix
C     Add degree of vertex to diagonal
      if(nalgorithm.eq.1) then
      Write(Iout,1039)
      Do I=1,MAtom
      Do J=1,MAtom
      A(I,J)=-dfloat(IDA(I,J))
      enddo
       A(I,I)=3.d0
      enddo
C     Diagonalize
      call tred2(A,Matom,NMax,evec,df)
      call tqli(evec,df,Matom,NMax,A)
      Write(Iout,1004) Matom,Matom
C     Sort eigenvalues evec(i) and eigenvectors A(*,i)
C     Sorting is different to adjacency matrix
C      Here from the lowest to highest eigenvalue
      Do I=1,MAtom
      e0=evec(I)
      jmax=I
      Do J=I+1,MAtom
      e1=evec(J)
      if(e1.lt.e0) then 
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
      endif

C Now produce the 3D image
C   Algorithm 1 or 2 (Fowler-Manopoulus):
C     Now search for lowest energy P-type vectors
C     This needs to be changed
      if(nalgorithm.le.1) then
      icand=0
      Do I=1,iocc
      mneg=0
      mpos=0
      z=0.d0
      Do J=1,MATOM
      if(A(J,I).lt.-1.d-9) mneg=mneg+1
      if(A(J,I).gt.1.d-9) mpos=mpos+1
      enddo
      if(mneg.eq.mpos) then
      icand=icand+1
      idg(icand)=I
      endif
      enddo
C     Analyzing remaining occupied eigenvectors
C     Create cartesian coordinates (in Angstroems) and scale them
C     if required
      Ix1=idg(1)
      Ix2=idg(2)
      Ix3=idg(3)
      if(Ix1.ne.2.or.Ix2.ne.3.or.Ix3.ne.4) then
      Write(Iout,1027) Ix1,Ix2,Ix3
      endif
      I1=IV1
      I2=IV2
      I3=IV3
      Write(Iout,1028) I1,I2,I3
      Do I=1,MATOM
      Dist(1,I)=A(I,I1)
      Dist(2,I)=A(I,I2)
      Dist(3,I)=A(I,I3)
      enddo
      CALL Distan(NMAX,Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      ratiotest=Rminall/Rmax
C     Search for better eigenvectors (not implemented yet)
      if(ratiotest.lt.1.d-6) then
      Write(Iout,1033) ratiotest
      istop=1
      endif
      fac1=1.d0/dsqrt(3.d0-evec(I1))
      fac2=1.d0/dsqrt(3.d0-evec(I2))
      fac3=1.d0/dsqrt(3.d0-evec(I3))
      ratio1=(Rmax/Rmin-1.d0)*1.d2
      Do I=1,MATOM
      Dist(1,I)=A(I,I1)*fac1
      Dist(2,I)=A(I,I2)*fac2
      Dist(3,I)=A(I,I3)*fac3
      enddo
      CALL Distan(NMAX,Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      ratio=(Rmax/Rmin-1.d0)*1.d2
      if(ratio1.lt.ratio) then
      Write(Iout,1026)
      Do I=1,MATOM
      Dist(1,I)=A(I,I1)/fac1
      Dist(2,I)=A(I,I2)/fac2
      Dist(3,I)=A(I,I3)/fac3
      enddo
      fac1=1.d0
      fac2=1.d0
      fac3=1.d0
      CALL Distan(NMAX,Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      endif

C     Obtain smallest distance for further scaling
C     Now this contracts or expandes the whole fullerene to set the
C     smallest bond distance to Cdist
      R0=1.d10
      Do I=1,MATOM
      Do J=I+1,MATOM
      X=Dist(1,I)-Dist(1,J)
      Y=Dist(2,I)-Dist(2,J)
      Z=Dist(3,I)-Dist(3,J)
      R=dsqrt(X*X+Y*Y+Z*Z)
      if(R.lt.R0) R0=R
      enddo
      enddo
      fac=CDist/R0
      Write(Iout,1010) icand,I1,I2,I3,fac1,fac2,fac3,fac,Cdist
      Do k=1,2*IL3
      ICN=0
      ICP=0
      do kk=1,Matom
      if(A(kk,k).gt.0.d0) ICP=ICP+1
      if(A(kk,k).lt.0.d0) ICN=ICN+1
      enddo
      Write(Iout,1011) k,evec(k),ICN,ICP
      Write(Iout,1012) (A(J,k),J=1,Matom)
      enddo
      if(R0.lt.1.d-5.or.istop.eq.1) then
      Write(IOUT,1032) R0,fac
      stop
      endif
      Do I=1,MATOM
      Dist(1,I)=Dist(1,I)*fac
      Dist(2,I)=Dist(2,I)*fac
      Dist(3,I)=Dist(3,I)*fac
      enddo
C     Check distances
      Write(IOUT,1013)     
      Do J=1,MAtom
      Write(IOUT,1014) J,(Dist(I,J),I=1,3)
      enddo
      CALL Distan(NMAX,Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1015) Rmin,Rmax,rms
      ratio=(Rmax/Rmin-1.d0)*1.d2
      iratio=dint(ratio)
      CALL Diameter(NMax,MAtom,IOUT,Dist,distp)
      if(iratio.lt.33) then
      Write(IOUT,1016) iratio
      else
      Write(IOUT,1029) iratio
      endif
C     Calculate P-type dipole moment
       Call Dipole(NMax,MAtom,I1,I2,I3,IOUT,dipol,Dist,A)
       Write(IOUT,1030)
       Do I=1,3
        Write(IOUT,1031) I,(dipol(I,J),J=1,3)
       enddo

      else
C   Algorithm 2 (Tutte):
C     Input: Integer Adjacency Matrix IDA(NMax,NMax)
C     Output: Real*8 Cartesian Coordinates  Dist(3,NMax)
C     NMax: Max Dimension of Matrix
C     MAtom: Working Dimension of Matrix

C     Algorithm 2 (Tutte):
         write (Iout,1017) 
         g = new_fullerene_graph(NMax,MAtom,IDA)
         write (Iout,1018)
         call tutte_layout(g,layout2d)
         call spherical_layout(g,Dist)
         write (Iout,1035)
         call delete_fullerene_graph(g)
         write (Iout,1036)

C     Obtain smallest distance for further scaling
C     Now this contracts or expands the whole fullerene to set the
C     smallest bond distance to Cdist
      R0=1.d10
      Do I=1,MATOM
      Do J=I+1,MATOM
      X=Dist(1,I)-Dist(1,J)
      Y=Dist(2,I)-Dist(2,J)
      Z=Dist(3,I)-Dist(3,J)
      R=dsqrt(X*X+Y*Y+Z*Z)
      if(R.lt.R0) R0=R
      enddo
      enddo
      fac=CDist/R0
      Do I=1,MATOM
      Dist(1,I)=Dist(1,I)*fac
      Dist(2,I)=Dist(2,I)*fac
      Dist(3,I)=Dist(3,I)*fac
      enddo
C     Check distances
      Write(IOUT,1034) fac     
      Do J=1,MAtom
      Write(IOUT,1014) J,(Dist(I,J),I=1,3)
      enddo
      CALL Distan(NMAX,Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1015) Rmin,Rmax,rms
      ratio=(Rmax/Rmin-1.d0)*1.d2
      iratio=dint(ratio)
      CALL Diameter(NMax,MAtom,IOUT,Dist,distp)
      if(iratio.lt.33) then
      Write(IOUT,1016) iratio
      else
      Write(IOUT,1029) iratio
      endif

      endif

 1000 FORMAT(/1X,'Cannot produce dual matrix, error IER= ',I2,
     1 ' Check your input for pentagon locations')
 1001 FORMAT(/1X,'Program to create cartesian coordinates through ',
     1 'pentagon index list producing the dual matrix and finally '
     1 'the Hueckel matrix',/1X,'Number of faces: ',I3,
     1 ', Number of atoms (vertices): ',I4,
     1 ', Point group of fullerene (in ideal symmetry): ',A3,/1X,
     1 'Ring spiral pentagon positions: ',12I4,
     1 /1X,'NMR pattern: ',3(I3,' x',I3,:,','))
 1002 FORMAT(/1X,'D contains IER = ',I6,' separating triangles and is ',
     1 'therefore NOT a fullerene dual')
 1003 FORMAT(1X,'       x     deg NE   type    ',/1X,32('-'))
 1004 FORMAT(/1X,'Using the Fowler-Manopoulus algorithm to construct ',
     1 'the fullerene',/1X,'Construct the (',I3,','I3,') Hueckel ',
     1 ' matrix, diagonalize (E=alpha+x*beta) and get eigenvectors',
     1 /1X,'Eigenvalues are between [-3,+3]')
 1005 FORMAT(1X,F12.6,I3,1X,I3,3X,A10)
 1006 FORMAT(1X,32('-'))
 1007 FORMAT(1X,'Fullerene has open-shell character (zero band gap)!')
 1008 FORMAT(1X,'Bandgap delta x = ',F12.6,' (in units of |beta|)')
 1009 FORMAT(1X,'Caution: Bandgap small, possibility '
     1 'for open-shell character')
 1010 FORMAT(1X,I3,' potential candidates in occupied space discovered',
     1 ' for odd symmetry eigenvectors, and ',3I4,' taken',
     1 /1X,'Create cartesian coordinates and scale a la Fowler ',
     1 /1X,'Scaling factors sx =',D12.6,', sy =',D12.6,', sz =',D12.6,
     1 ', and final scaling factor of ',D16.8,' setting Rmin= ',F12.6,
     1 /1X,'P-type eigenvalues and eigenvectors (path search not yet ',
     1 'implemented and cartesian coordinates may not be correct)',/1X,
     1 'Note: Procedure may give bond lengths which may vary strongly')
 1011 FORMAT(1X,'eigenvalue',I4,': ',F12.6,', eigenvector: (',I3,
     1 ' negative and ',I3,' positive values)')
 1012 FORMAT(10(1X,F12.6))
 1013 FORMAT(1X,'Fowler-Manopoulus Coordinates')
 1014 FORMAT(1X,I3,5X,3(D18.12,2X))
 1015 FORMAT(1X,'Minimum distance: ',F12.6,', Maximum distance: ',F12.6,
     1 ', RMS distance: ',F12.6)
 1016 Format(1X,'Maximum bond distance ',I5,'% larger than minimum ',
     1 'distance')
 1017 FORMAT(/1X,'Using the Tutte-embedding algorithm to construct ',
     1 'the fullerene',/1X,'Construct the Tutte planar graph and ',
     1 'project on sphere')
 1018 Format(1X,'Calculating Tutte-embedding') 
 1019 Format(1X,'Spiral list of pentagon positions with ',
     1 'higher priority: (',I3,' spirals found)')
 1020 FORMAT(/1X,'Program to create cartesian coordinates through ',
     1 'pentagon index list producing the dual matrix and finally '
     1 'the Hueckel matrix',/1X,'Number of faces: ',I3,
     1 ', Number of atoms (vertices): ',I4,
     1 ', Point group of fullerene (in ideal symmetry): ',A3,/1X,
     1 'Ring spiral pentagon positions: ',12I4)
 1021 Format(12(1X,I3))
 1022 Format(1X,'Input spiral is canonical')
 1023 Format(1X,'Canonical spiral list of pentagon positions:')
 1024 Format(1X,'Canonical spiral list of hexagons and pentagons:')
 1025 Format(1X,100I1)
 1026 Format(1X,'No Fowler scaling')
 1027 Format(1X,'Analysis of ',3I3,' eigenvectors finds P-type ',
     1 'eigenvectors not in sequence 2 3 4',
     1 /,' You might analyze eigenvectors and read in numbers for '
     1 'eigenvectors to be used to construct coordinates instead')
 1028 Format(1X,'Take ',3I3,' eigenvectors as required by input')
 1029 Format(1X,'Maximum bond distance ',I5,'% larger than minimum ',
     1 'distance. Fullerene strongly distorted!',/1X,
     1 'Optimization of geometry recommended')
 1030 Format(1X,'P-type dipole moments as indicators:')
 1031 Format(1X,I1,2X,'(',F9.2,',',F9.2,',',F9.2,')')
 1032 Format(1X,'Very small bond distance of ',D18.12,' found',
     1 ' leading to a unreasonable scaling factor of ',D18.12,
     1 /1X,'**** Program stops (choose other set of eigenvectors)')
 1033 Format(1X,'Rmin/Rmax ratio too small: ',D18.12,
     1 /1X,'Program will stop (choose other set of eigenvectors)')
 1034 FORMAT(1X,'Coordinates from Tutte embedding scaled by a factor'
     1 ' of ',D18.12)
 1035 FORMAT(1X,'Projected on sphere')
 1036 FORMAT(1X,'Fullerene graph deleted')
 1037 FORMAT(1X,'Graph is not cubic, ',I4,' vertices detected which ',
     1 'are not of degree 3, last one is of degree ',I4)
 1038 FORMAT(1X,'Graph checked, it is cubic')
 1039 FORMAT(1X,'Laplacian Matrix taken instead of adjacency matrix')
      Return 
      END

      SUBROUTINE Dipole(NMax,MAtom,I1,I2,I3,IOUT,dipol,Dist,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,NMAX),A(NMAX,NMAX),dipol(3,3)
      tol=1.d-7
      do I=1,3
      do j=1,3
      dipol(I,J)=0.d0
      enddo
      enddo
      do I=1,MAtom
      do j=1,3
      dipol(1,j)=dipol(1,j)+Dist(j,I)*A(I,I1)
      dipol(2,j)=dipol(2,j)+Dist(j,I)*A(I,I2)
      dipol(3,j)=dipol(3,j)+Dist(j,I)*A(I,I3)
      enddo
      enddo
      do i=1,3
      do j=1,3
      if(dabs(dipol(i,j)).lt.tol) dipol(i,j)=0.d0
      enddo
      enddo
      Return 
      END
