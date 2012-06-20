      SUBROUTINE CoordBuild(MAtom,IN,Iout,IDA,D,ICart,
     1 IV1,IV2,IV3,kGC,lGC,isonum,IPRC,ihueckel,JP,iprev,
     1 A,evec,df,Dist,layout2d,distp,Cdist,GROUP)
C Cartesian coordinates produced from ring spiral pentagon list
C or Coxeter-Goldberg construction to get the adjacency matrix
C This is followed by using either the Fowler-Manolopoulos matrix
C eigenvector or the Tutte 3D embedding algorithm 
C Fowler-Manolopoulos matrix eigenvector algorithm: identify P-type 
C eigenvectors and construct the 3D fullerene
C Tutte embedding algorithm: Tutte barycentric embedding and sphere
C mapping 
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 layout2d
      Integer D,S,Spiral
      integer :: graph_is_a_fullerene
      DIMENSION layout2d(2,NMAX)
      DIMENSION D(MMAX,MMAX),S(MMAX),Dist(3,NMAX),distP(NMAX)
      DIMENSION NMR(6),JP(12),A(NMAX,NMAX),IDA(NMAX,NMAX)
      DIMENSION evec(NMAX),df(NMAX),dipol(3,3)
      DIMENSION Spiral(12,NMAX)
      DIMENSION idg(NMAX)
      CHARACTER*3 GROUP
      Data Tol,Tol1,Tol2,ftol/1.d-5,.15d0,1.5d1,1.d-10/
      type(c_ptr) :: g, halma, new_fullerene_graph,
     1               new_C20, halma_fullerene
C If nalgorithm=0 use ring-spiral and matrix eigenvector algorithm
C If nalgorithm=1 use ring-spiral and Tutte algorithm
C If nalgorithm=2 use Coxeter Goldberg and matrix eigenvector 
C                 algorithm
C If nalgorithm=3 use Coxeter Goldberg and Tutte algorithm
      nalgorithm=ICart-2

      M=Matom/2+2
      istop=0
      Group='NA '


C Ring spiral first:
C Read pentagon list and produce adjacency matrix
      if(nalgorithm.le.1) then
      if(isonum.eq.0) then
        if(iprev.eq.0) Read(IN,*) (JP(I),I=1,12)
      else
C Read from database
        Call Isomerget(Matom,Iout,Isonum,IPRC,JP)
      endif
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
      IPRS=0
      IER=0
      CALL Windup(M,IPRS,IER,S,D)              ! Wind up spiral into dual
      IF(IER.gt.0) then
      WRITE(Iout,1000) IER
      endif
      IT=1
      Do I=1,12
      Spiral(I,1)=JP(I)
      enddo
      CALL Unwind(M,IER,IT,ispiral,
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
      CALL CanSpiral(ispiral,spiral,JP)
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

C Now produce the adjaceny matrix from the dual matrix
      CALL DUAL(D,MMAX,IDA,Matom,IER)
      IF(IER.ne.0) then
      WRITE(Iout,1002) IER
      stop
      endif

C End ring spiral
      endif


C Start Goldberg-Coxeter
      if(nalgorithm.gt.1) then
      Write(Iout,1040) kGC,lGC,kGC,lGC
      if(lGC .ne. 0) then
        write(Iout,1041)
        stop
      endif
      g = new_C20();
      halma = halma_fullerene(g,kGC-1)
      isafullerene = graph_is_a_fullerene(halma)
      IF (isafullerene .eq. 1) then
        write (iout,1043) 
      else
        write (iout,1044)
        stop
      endif
C Update fortran structures
      MAtom  = NVertices(halma)
      Medges = NEdges(halma)
        write(Iout,1042)  MAtom,Medges
      call adjacency_matrix(halma,NMax,IDA)
C End Goldberg-Coxeter
      endif


C Adjacency matrix constructed
C Now analyze the adjacency matrix if it is correct
      Do I=1,MAtom
      Do J=1,MAtom
       A(I,J)=dfloat(IDA(I,J))
      enddo
      enddo
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
      if(ihueckel.eq.0.or.nalgorithm.eq.0.or.nalgorithm.eq.2) then
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

C Analyze eigenenergies
      Call HueckelAnalyze(MAtom,NMax,Iout,iocc,df,evec)
C     End of Hueckel
      endif


C Now produce the 3D image

C   Fowler-Manolopoulos matrix eigenvector algorithm
C     Now search for lowest energy P-type vectors
C     This needs to be changed
      if(nalgorithm.eq.0.or.nalgorithm.eq.2) then
      Write(Iout,1045)
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
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
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
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
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
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
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
      Do k=1,2*I3
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
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1015) Rmin,Rmax,rms
      ratio=(Rmax/Rmin-1.d0)*1.d2
      iratio=dint(ratio)
      CALL Diameter(MAtom,IOUT,Dist,distp)
      if(iratio.lt.33) then
      Write(IOUT,1016) iratio
      else
      Write(IOUT,1029) iratio
      endif
C     Calculate P-type dipole moment
       Call Dipole(MAtom,I1,I2,I3,dipol,Dist,A)
       Write(IOUT,1030)
       Do I=1,3
        Write(IOUT,1031) I,(dipol(I,J),J=1,3)
       enddo
      endif

      if(nalgorithm.eq.1.or.nalgorithm.eq.3) then
C   Tutte 3D embedding Algorithm:
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
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1015) Rmin,Rmax,rms
      ratio=(Rmax/Rmin-1.d0)*1.d2
      iratio=dint(ratio)
      CALL Diameter(MAtom,IOUT,Dist,distp)
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
     1 'the Hueckel matrix',/1X,'Number of faces: ',I4,
     1 ', Number of atoms (vertices): ',I4,
     1 ', Point group of fullerene (in ideal symmetry): ',A3,/1X,
     1 'Ring spiral pentagon positions: ',12I5,
     1 /1X,'NMR pattern: ',3(I3,' x',I3,:,','))
 1002 FORMAT(/1X,'D contains IER = ',I6,' separating triangles and is ',
     1 'therefore NOT a fullerene dual')
 1004 FORMAT(/1X,'Construct the (',I4,','I4,') Hueckel ',
     1 ' matrix, diagonalize (E=alpha+x*beta) and get eigenvectors',
     1 /1X,'Eigenvalues are between [-3,+3]')
 1010 FORMAT(1X,I4,' potential candidates in occupied space discovered',
     1 ' for odd symmetry eigenvectors, and ',3I4,' taken',
     1 /1X,'Create cartesian coordinates and scale a la Fowler ',
     1 /1X,'Scaling factors sx =',D12.6,', sy =',D12.6,', sz =',D12.6,
     1 ', and final scaling factor of ',D16.8,' setting Rmin= ',F12.6,
     1 /1X,'P-type eigenvalues and eigenvectors (path search not yet ',
     1 'implemented and cartesian coordinates may not be correct)',/1X,
     1 'Note: Procedure may give bond lengths which may vary strongly')
 1011 FORMAT(1X,'eigenvalue',I4,': ',F12.6,', eigenvector: (',I4,
     1 ' negative and ',I4,' positive values)')
 1012 FORMAT(10(1X,F12.6))
 1013 FORMAT(1X,'Fowler-Manolopoulos Coordinates')
 1014 FORMAT(1X,I4,5X,3(D18.12,2X))
 1015 FORMAT(1X,'Minimum distance: ',F12.6,', Maximum distance: ',F12.6,
     1 ', RMS distance: ',F12.6)
 1016 Format(1X,'Maximum bond distance ',I5,'% larger than minimum ',
     1 'distance')
 1017 FORMAT(/1X,'Using the Tutte-embedding algorithm to construct ',
     1 'the fullerene',/1X,'Construct the Tutte planar graph and ',
     1 'project on sphere')
 1018 Format(1X,'Calculating Tutte-embedding') 
 1019 Format(1X,'Spiral list of pentagon positions with ',
     1 'higher priority: (',I5,' spirals found)')
 1020 FORMAT(/1X,'Program to create cartesian coordinates through ',
     1 'pentagon index list producing the dual matrix and finally '
     1 'the Hueckel matrix',/1X,'Number of faces: ',I4,
     1 ', Number of atoms (vertices): ',I4,
     1 ', Point group of fullerene (in ideal symmetry): ',A3,/1X,
     1 'Ring spiral pentagon positions: ',12I4)
 1021 Format(12(1X,I5))
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
 1040 Format(/1x,'Goldberg-Coxeter fullerene with indices (k,l) = (',
     1 I2,',',I2,') taking C20 as the input graph: GC(',I2,',',I2,
     1 ')[G0] with G0=C20')
 1041 Format(/1x,'Goldberg-Coxeter construction not implemented',
     1 ' for l > 0.')
 1042 Format(1x,'Updating number of vertices (',I5,') and edges (',
     1 I5,')')
 1043 Format(1x,'Halma fullerene is a fullerene')
 1044 Format(1x,'Halma fullerene is not a fullerene')
 1045 FORMAT(/1X,'Using the Fowler-Manolopoulos adjacency matrix '
     1 'eigenvector algorithm to construct the fullerene')
      Return 
      END

      SUBROUTINE Dipole(MAtom,I1,I2,I3,dipol,Dist,A)
      use config
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

      SUBROUTINE Isomerget(Matom,Iout,Isonum,IPR,JP)
      IMPLICIT Integer (A-Z)
C Routine to get Isomer number from database
      Real*8 sigmah
      Character*31 databasefile
      Character*13 dbdir
      Character*9 fend
      Character*1 fstart,fnum1,Dummy
      Character*2 fnum2
      Character*3 fnum,fnum3,GROUP,ident
      Dimension JP(12)
      Dimension Iso(119),IsoIPR(123)
      Logical lexist
      Data Iso/1,0,1,1,2,3,6,6,15,17,40,45,89,116,199,
     * 271,437,580,924,1205,1812,2385,3465,4478,6332,8149,
     * 11190,14246,19151,24109,31924,39718,51592,63761,81738,
     * 99918,126409,153493,191839,231017,285914,341658,419013,
     * 497529,604217,713319,860161,1008444,1207119,1408553,
     * 1674171,1942929,2295721,2650866,3114236,3580637,4182071,
     * 4787715,5566949,6344698,7341204,8339033,9604411,10867631,
     * 12469092,14059174,16066025,18060979,20558767,23037594,
     * 26142839,29202543,33022573,36798433,41478344,46088157,
     * 51809031,57417264,64353269,71163452,79538751,87738311,
     * 97841183,107679717,119761075,131561744,145976674,159999462,
     * 177175687,193814658,214127742,233846463,257815889,281006325,
     * 309273526,336500830,369580714,401535955,440216206,477420176,
     * 522599564,565900181,618309598,668662698,729414880,787556069,
     * 857934016,925042498,1006016526,1083451816,1176632247,
     * 1265323971,1372440782,1474111053,1596482232,1712934069,
     * 1852762875,1985250572,2144943655/
      Data IsoIPR/1,0,0,0,0,1,1,1,2,5,7,9,24,19,35,46,86,134,187,
     * 259,450,616,823,1233,1799,2355,3342,4468,6063,8148,10774,13977,
     * 18769,23589,30683,39393,49878,62372,79362,98541,121354,151201,
     * 186611,225245,277930,335569,404667,489646,586264,697720,
     * 836497,989495,1170157,1382953,1628029,1902265,2234133,
     * 2601868,3024383,3516365,4071832,4690880,5424777,6229550,
     * 7144091,8187581,9364975,10659863,12163298,13809901,15655672,
     * 17749388,20070486,22606939,25536557,28700677,32230861,
     * 36173081,40536922,45278722,50651799,56463948,62887775,
     * 69995887,77831323,86238206,95758929,105965373,117166528,
     * 129476607,142960479,157402781,173577766,190809628,209715141,
     * 230272559,252745513,276599787,303235792,331516984,362302637,
     * 395600325,431894257,470256444,512858451,557745670,606668511,
     * 659140287,716217922,776165188,842498881,912274540,987874095,
     * 1068507788,1156161307,1247686189,1348832364,1454359806,
     * 1568768524,1690214836,1821766896,1958581588,2109271290/
      LimitAll=100
      LimitIPR=120
      nhamcycle=0

C     Check if parameters are set correctly
      If(IPR.lt.0) then
        Write(Iout,1008)
        Stop
      endif
      If(IPR.eq.0) then
       if(MAtom.gt.LimitAll) then
        Write(Iout,1000) MAtom,LimitAll
        stop
       endif
       M1=MAtom/2-9
       isoN=Iso(M1)
       if(isonum.gt.isoN) then
        Write(Iout,1003) isonum,isoN
        stop
       endif
       dbdir='database/All/'
       itotal=isoN
       ident='GEN'

      else

       if(MAtom.gt.LimitIPR) then
        Write(Iout,1001) MAtom,LimitIPR
        stop
       endif
       if(MAtom.lt.60.or.(MAtom.gt.60.and.MAtom.lt.70)) then
        Write(Iout,1002) MAtom
        stop
       endif
       M2=MAtom/2-29
         isoIPRN=IsoIPR(M2)
       if(isonum.gt.isoIPRN) then
        Write(Iout,1004) isonum,isoIPRN
        stop
       endif
       dbdir='database/IPR/'
       itotal=isoIPRN
       ident='IPR'

      endif

C     Now produce filename in database
      fend='.database'
      fstart='c'
      if(Matom.lt.100) then
      fnum1='0'
       write(fnum2,'(I2)') MAtom
       fnum=fnum1//fnum2
      else
       write(fnum,'(I3)') MAtom
      endif
      fnum3='all'
      if(IPR.eq.1) fnum3='IPR'
      databasefile=dbdir//fstart//fnum//fnum3//fend
      if(IPR.eq.0) then 
       Write(Iout,1005) MAtom,isonum,databasefile
      else
       Write(Iout,1006) MAtom,isonum,databasefile
      endif

C     Open file
      inquire(file=databasefile,exist=lexist)
      if(lexist.neqv..True.) then
       Write(Iout,1007) databasefile
       stop
      endif
      Open(UNIT=4,FILE=databasefile,STATUS='old',FORM='FORMATTED')
       Read(4,*) IN,IP,IH
      Write(Iout,1009) itotal
      Do I=1,isonum-1
       Read(4,*) Dummy
      enddo
      if(IH.eq.0) then
       Read(4,2000) L,GROUP,(JP(J),J=1,12),IFus5G,sigmah
      if(L.lt.100) Write(Iout,1010) ident,L,GROUP
      if(L.ge.100.and.L.lt.1000) Write(Iout,1011) ident,L,GROUP
      if(L.ge.1000.and.L.lt.10000) Write(Iout,1012) ident,L,GROUP
      if(L.ge.10000.and.L.lt.100000) Write(Iout,1013) ident,L,GROUP
      if(L.ge.100000.and.L.lt.1000000) Write(Iout,1014) ident,L,GROUP
      if(L.ge.1000000.and.L.lt.10000000) Write(Iout,1015) ident,L,GROUP
      if(L.ge.10000000.and.L.lt.100000000) Write(Iout,1016)ident,L,GROUP
      if(L.ge.100000000) Write(Iout,1017) ident,L,GROUP
      Write(Iout,1020) (JP(J),J=1,12),IFus5G,sigmah
      else
       Read(4,2001) L,GROUP,(JP(J),J=1,12),IFus5G,sigmah,nhamcycle
      if(L.lt.100) Write(Iout,1010) ident,L,GROUP
      if(L.ge.100.and.L.lt.1000) Write(Iout,1011) ident,L,GROUP
      if(L.ge.1000.and.L.lt.10000) Write(Iout,1012) ident,L,GROUP
      if(L.ge.10000.and.L.lt.100000) Write(Iout,1013) ident,L,GROUP
      if(L.ge.100000.and.L.lt.1000000) Write(Iout,1014) ident,L,GROUP
      if(L.ge.1000000.and.L.lt.10000000) Write(Iout,1015) ident,L,GROUP
      if(L.ge.10000000.and.L.lt.100000000) Write(Iout,1016)ident,L,GROUP
      if(L.ge.100000000) Write(Iout,1017) ident,L,GROUP
      Write(Iout,1021) (JP(J),J=1,12),IFus5G,sigmah,nhamcycle
      endif
      Close(unit=4)

 1000 Format(1X,'Number of atoms ',I4,' larger than data base limit ',
     1 I3,' for general isomers ==> ABORT')
 1001 Format(1X,'Number of atoms ',I4,' larger than data base limit ',
     1 I3,' for IPR isomers isomers ==> ABORT')
 1002 Format(1X,'Number of atoms ',I4,' does not match IPR fullerene ',
     1 '==> ABORT')
 1003 Format(1X,'Isomer number ',I9,' too large to be in data base:',
     1 ' number of isomers in ALL data base ',I9,' ==> ABORT')
 1004 Format(1X,'Isomer number ',I9,' too large to be in data base:',
     1 ' number of isomers in IPR data base ',I9,' ==> ABORT')
 1005 Format(1X,'Number of atoms ',I5,' and isomer number ',I9,
     1 ' in general isomer list',/1X,'Search in file: ',A29)
 1006 Format(1X,'Number of atoms ',I5,' and isomer number ',I9,
     1 ' in IPR isomer list',/1X,'Search in file: ',A29)
 1007 Format(1X,'Filename ',A29,' in database not found ==> ABORT')
 1008 Format(1X,'IPR parameter not set ==> ABORT')
 1009 Format(1X,'File contains ',I9,' isomers')
 1010 Format(1X,'Read isomer ',A3,I2,' from database (point group: ',
     1 A3,')')
 1011 Format(1X,'Read isomer ',A3,I3,' from database (point group: ',
     1 A3,')')
 1012 Format(1X,'Read isomer ',A3,I4,' from database (point group: ',
     1 A3,')')
 1013 Format(1X,'Read isomer ',A3,I5,' from database (point group: ',
     1 A3,')')
 1014 Format(1X,'Read isomer ',A3,I6,' from database (point group: ',
     1 A3,')')
 1015 Format(1X,'Read isomer ',A3,I7,' from database (point group: ',
     1 A3,')')
 1016 Format(1X,'Read isomer ',A3,I8,' from database (point group: ',
     1 A3,')')
 1017 Format(1X,'Read isomer ',A3,I9,' from database (point group: ',
     1 A3,')')
 1020 Format(/1X,'Ring spiral pentagon indices: ',12I5,
     1 /1X,'Np= ',I2,', sigmah= ',F12.6)
 1021 Format(/1X,'Ring spiral pentagon indices: ',12I5,
     1 /1X,'Np= ',I2,', sigmah= ',F12.6,', Number of distinct ',
     1 'Hamiltonian cycles: ',I10)
 2000 Format(I9,2X,A3,1X,12I4,23X,I2,27X,F8.5)
 2001 Format(I9,2X,A3,1X,12I4,23X,I2,27X,F8.5,25X,I9)

      Return 
      END
