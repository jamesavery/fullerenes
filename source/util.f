c      SUBROUTINE PerfectMatching(Iout,IDA)
c      use config
c      IMPLICIT REAL*8 (A-H,O-Z)
c      DIMENSION IDA(Nmax,Nmax)
c       Write(Iout,1000) number_vertices
c       Write(Iout,1001) 
c       Write(Iout,1002) 
c 1000 Format(1X,'Upper limit for number of perfect matchings',
c     1 ' in cubic graphs: 2**N with N=',I5)
c 1001 Format(1X,'Counting the number of perfect matchings using',
c     1 ' the Fisher-Kasteleyn-Temperley (FKT) algorithm')
c 1002 Format(1X,'Not implemented yet')
c      RETURN
c      END

      SUBROUTINE CompressDatabase(Iout,filename)
C This routine turns a database file from the output into a new compreesed one
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer RSPI(12),PNI(0:5),HNI(0:6),INMR(6)
      CHARACTER*50 filename,dbfilename,dbfilenew
      CHARACTER*26 nmrstring
      CHARACTER*18 dbstring,readstring
      CHARACTER*6  shell
      CHARACTER*3  Group
      CHARACTER*1  dummy
      Logical lexist

      dbfilename=trim(filename)
      inquire(file=dbfilename,exist=lexist)
      if(lexist.neqv..True.) then
        Write(Iout,1000) dbfilename
        close(unit=1)
        Return
      else
        Open(UNIT=1,FILE=dbfilename,STATUS='old',Action='Read',
     1   FORM='FORMATTED')
        Write(Iout,1001) dbfilename
        dbfilenew=trim(filename)//".new"
        Open(UNIT=2,FILE=dbfilenew,STATUS='unknown',FORM='FORMATTED')
        Read(1,*,Iostat=ierr) Nvert,IP,IH

C First line not 3 integers, going to search for it
        if(ierr.ne.0) then
         dbstring=' Isomer List Start'
         Write(Iout,1002) dbstring
         do I=1,1000
          Read(1,FMT='(A18)',ERR=199) readstring
          if(readstring.eq.dbstring) go to 10
         enddo
 10      Read(1,*,Iostat=ierr1) Nvert,IP,IH
         do I=1,5    
          Read(1,FMT='(A1)') dummy
         enddo
        endif

C Now read database
        Write(Iout,1003) Nvert,IP,IH
        Write(2,1005) Nvert,IP,IH
        nlines=0
        Do J=1,Nisoloop
         if(IH.eq.1) then
         Read(1,2000,ERR=199) number,Group,(RSPI(i),I=1,12),
     1    (PNI(I),I=0,5),NP,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,
     1    HLgap,shell,ncycHam,nmrstring
         else
         Read(1,2001,ERR=199) number,Group,(RSPI(i),I=1,12),
     1    (PNI(I),I=0,5),NP,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,
     1    HLgap,shell,nmrstring
         endif
         nlines=nlines+1
         if(nlines.ne.number) then
          Write(Iout,1010) number,nlines
          return
         endif
        lenstring=LEN_TRIM(nmrstring)
        if(lenstring.le.9) then
         read(nmrstring(1:3),'(i3)') INMR(1)
         read(nmrstring(7:9),'(i3)') INMR(2)
         INMR(3)=0
         INMR(4)=0
         INMR(5)=0
         INMR(6)=0
        else
         if(lenstring.gt.20) then
          read(nmrstring(1:3),'(i3)') INMR(1)
          read(nmrstring(6:8),'(i3)') INMR(2)
          read(nmrstring(10:12),'(i3)') INMR(3)
          read(nmrstring(15:17),'(i3)') INMR(4)
          read(nmrstring(19:21),'(i3)') INMR(5)
          read(nmrstring(24:26),'(i3)') INMR(6)
         else
          read(nmrstring(1:3),'(i3)') INMR(1)
          read(nmrstring(6:8),'(i3)') INMR(2)
          read(nmrstring(10:12),'(i3)') INMR(3)
          read(nmrstring(15:17),'(i3)') INMR(4)
          INMR(5)=0
          INMR(6)=0
         endif
        endif
        if(IP.eq.0) then
         if(IH.eq.1) then
Case 1 All isomers with Hamiltonian cycles IP=0 IH=1
          Write(2,1004) Group,(RSPI(i),I=1,12),(PNI(I),I=0,4),
     1     (HNI(I),I=0,5),NeHOMO,NedegHOMO,HLgap,ncycHam,(INMR(I),I=1,6)
         else
Case 2 All isomers without Hamiltonian cycles IP=0 IH=0
          Write(2,1007) Group,(RSPI(i),I=1,12),(PNI(I),I=0,4),
     1     (HNI(I),I=0,5),NeHOMO,NedegHOMO,HLgap,(INMR(I),I=1,6)
         endif
        else
         if(IH.eq.1) then
Case 3 IPR isomers with Hamiltonian cycles IP=1 IH=1
          Write(2,1008) Group,(RSPI(i),I=1,12),(HNI(I),I=3,5),
     1     NeHOMO,NedegHOMO,HLgap,ncycHam,(INMR(I),I=1,6)
         else
Case 4 IPR isomers without Hamiltonian cycles IP=1 IH=0
          Write(2,1009) Group,(RSPI(i),I=1,12),(HNI(I),I=3,5),
     1     NeHOMO,NedegHOMO,HLgap,(INMR(I),I=1,6)
         endif
        endif
        enddo
      endif

  199 close(unit=1)
      close(unit=2)
      Write(Iout,1006) nlines

 1000 Format(1X,'Cannot find database file ',A50,' ===> RETURN')
 1001 Format(1X,'Open database file ',A50)
 1002 Format(1X,'Search for starting point of database',
     1 /1X,'Searching for string',A18)
 1003 Format(1X,'Number of vertices: ',I6,', IPR flag: ',I2,
     1 ', Hamiltonian cycle flag: ',I2)
 1004 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,I7,6I3)
 1005 Format(I3,2I1)
 1006 Format(1X,'Number of isomers written to new database: ',I10)
 1007 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,6I3)
 1008 Format(A3,12I3,3I2,I2,I1,F7.5,I7,6I3)
 1009 Format(A3,12I3,3I2,I2,I1,F7.5,6I3)
 1010 Format(1X,'Isomer number ',I10,' not identical to list number ',
     1 I10)
 2000 Format(I9,2X,A3,1X,12I4,3X,6(I2,1X),2X,I2,3X,6(I2,1X),
     1 I3,1X,F10.5,1X,2I3,F9.5,1X,A6,1X,I10,2X,A26)
 2001 Format(I9,2X,A3,1X,12I4,3X,6(I2,1X),2X,I2,3X,6(I2,1X),
     1 I3,1X,F10.5,1X,2I3,F9.5,1X,A6,2X,A26)
      RETURN
      END


      Function IPentInd(IRhag5)
      Integer IRhag5(0:5)
      Ifus5=0
       Do I=1,5
        IFus5=IFus5+I*IRhag5(I)
       enddo
       IPentInd=IFus5/2
      Return
      End

      Double Precision Function HexInd(IRhag6,ihk)
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer IRhag6(0:6)
      khk=0
      k2hk=0
      ihk=0
       Do I=0,6
        ihk=ihk+IRhag6(I)
        IIR=I*IRhag6(I)
        khk=khk+IIR
        k2hk=k2hk+I*IIR
       enddo
      if(ihk.eq.0) then
       Hexind=0.d0
       Return
      endif
      aihk=dfloat(ihk)
      akhk2=(dfloat(khk)/aihk)**2
      ak2hk=dfloat(k2hk)/aihk
      HexInd=dsqrt(dabs(ak2hk-akhk2))
      Return
      End

      SUBROUTINE Sortr(Mnew,imirror,jmirror,diam)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION imirror(Nmax),jmirror(Nmax),diam(Nmax)
      DIMENSION imirrorw(Nmax),jmirrorw(Nmax),diamw(Nmax)
      ICOUNT=0
      DO I=1,number_vertices
      DO K=I+1,number_vertices
C     Delete duplicates
       IF(I.eq.imirror(k).and.K.eq.imirror(i)) then
        ICOUNT=ICOUNT+1
        imirrorw(ICOUNT)=I
        jmirrorw(ICOUNT)=imirror(I)
        diamw(ICOUNT)=diam(I)
       endif
      enddo
      enddo
      Mnew=Icount
C     Now sort values of diamw, output diam
      DO I=1,MNew
       dMax=diamw(I)
       im=imirrorw(i)
       jm=jmirrorw(i)
       ivec=i
      DO K=I+1,MNew
       IF (dMax.LT.diamw(K)) then
        im=imirrorw(K)
        jm=jmirrorw(K)
        dMax=diamw(K)
        Ivec=K
       endif
      enddo
       imirror(i)=im
       jmirror(i)=jm
       diam(I)=dMax
       if(ivec.ne.i) then
        imirrorw(ivec)=imirrorw(i)
        jmirrorw(ivec)=jmirrorw(i)
        diamw(ivec)=diamw(i)
       endif
      enddo
      RETURN
      END

      SUBROUTINE TopIndicators(Iout,IDA,MDist)
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer MDist(Nmax,Nmax),Edges(2,3*number_vertices/2)
      DIMENSION IDA(Nmax,Nmax),wi(Nmax)
      type(c_ptr) :: graph, new_fullerene_graph
C     This routine calculates the Wiener index, Hyperwiener index,
C     minimal and maximal vertex contribution, rho and rhoE,
C     Schultz index and Balaban index
C     For details see D. Vukicevic,F. Cataldo, O. Ori, A. Graovac,
C     Chem. Phys. Lett. 501, 442–445 (2011).

      Write(Iout,1000) number_vertices

C     Get topological distance matrix
      graph = new_fullerene_graph(Nmax,number_vertices,IDA)
      call all_pairs_shortest_path(graph,number_vertices,Nmax,MDist)
      call edge_list(graph,edges,NE)
c     and finally delete the graph to free the mem
      call delete_fullerene_graph(graph)     

C     Wiener and hyper Wiener index, topological radius and diameter
      Xatom=dfloat(number_vertices)
      wiener1=0.d0
      wiener=0.d0
      wienermin=0
      wienermax=0
      hyperwiener=0.d0
      maxdist=0
      mRadius=100000000
      Do I=1,number_vertices
        wi(I)=0.d0
        maxdistrow=0
       Do J=1,number_vertices
        idist=MDist(I,J)
        adist=dfloat(idist)
        wi(I)=wi(i)+adist
        if(idist.gt.maxdistrow) maxdistrow=idist
         if(J.gt.I) then
           hyperwiener=hyperwiener+adist*(1.d0+adist)/2.d0
         endif
       enddo
       if(maxdistrow.gt.maxdist) maxdist=maxdistrow
       if(maxdistrow.lt.mRadius) mRadius=maxdistrow
       wiener1=wiener1+wi(i)
        if(I.ne.1) then
         if(wi(i).lt.wienermin) then
          wienermin=wi(i)
         endif
         if(wi(i).gt.wienermax) then
          wienermax=wi(i)
         endif
        else
         wienermin=wi(1)
         wienermax=wi(1)
        endif
      enddo

C     Balaban index
      balaban=0.d0
      Do I=1,number_vertices
      Do J=I+1,number_vertices
       if(IDA(I,J).eq.1) then
        wii=wi(I)
        wij=wi(J)
        if(wii.lt.1.d-15.or.wij.lt.1.d-15) then
         Write(Iout,1006)
         Return
        endif
        balaban=balaban+1.d0/(dsqrt(wii*wij))
       endif
      enddo
      enddo
      vertnum=dfloat(number_vertices)
      fac=3.d0*vertnum/(vertnum+4.d0)
      balabanindex=balaban*fac

C     Szeged index
      Call Szeged(Edges,MDist,Sz)

C     Final
      over=1.d-10
      wiener=wiener1/2.d0
      wav=wiener1/vertnum
      rho=wav/wienermin
      rhoE=wienermax/wienermin
      isize=number_vertices*(number_vertices-1)
      Avdist=2.d0*wiener/dfloat(isize)
      izagreb=number_vertices*9
      schultz=wiener*6.d0
      wienerfac=wiener/(9.d0*vertnum**3)
      Wienerbalaban=wienerfac*balabanindex*4.d0*(vertnum+4.d0)
C     ori=wiener/vertnum**2.5
      reversewiener=XATOM*(XATOM-1.d0)*dfloat(maxdist)/2.d0-wiener

      Write(Iout,1001) dint(wiener+over),
     1 dint(hyperwiener+over),
     1 dint(wienermin+over),
     1 dint(wienermax+over),
     1 dint(reversewiener+over),
     1 wav,rho,rhoE,izagreb,
     1 dint(schultz+over),
     1 balabanindex,dint(Sz+over)
      Write(Iout,1002) Wienerbalaban
      Write(Iout,1003) maxdist,mRadius,Avdist
C     Write(Iout,1004) ori

 1000 Format(1X,'Topological Indicators:',/1X,
     1 'For definitions see Vukicevic et al., Chem. Phys. Lett. ',
     1 '501, 442 (2011), and Behtoei et al., Appl. Math. Lett. ',
     1 '22, 1571 (2009)',/1X,'Cn with n=',I5)
 1001 Format(' Wiener index W: ',F20.0,/,' Hyper Wiener index WW: ',
     1 F20.0,/,' Minimal vertex contribution to W: Wmin= ',F20.0,
     1 ' Maximal vertex contribution to W: Wmax= ',F20.0,/,
     1 ' Reverse Wiener index n(n-1)D/2-W: ',F20.0,/,
     1 ' Average vertex contribution (wav): ',D15.9,/,
     1 ' rho: ',D15.9,', rhoE: ',D15.9,/,
     1 ' Zagreb index = nv*3^2 = ',I12,
     1 ' (trivial for regular fullerenes)',/,
     1 ' Schultz index = 6*W = ',F20.0,' (related to Wiener index for ',
     1 'regular fullerenes)',/,' Balaban index = ',D15.9,/,
     1 ' Szeged index = ',F20.0,
     1 /,' For the Estrada index see Subroutine Hueckel output')
 1002 Format(' f*Wiener*Balaban = 4WB(n+4)/(9n^3) = ',D15.9,/,'   ',
     1 ' (should be exactly 1.0 for cubic polyhedra with equal row ',
     1 'sums in distance matrix, i.e. Wmin=Wmax)')
 1003 Format(' Topological distances are between 1 and ',I6,
     1 ' (topological diameter D)',/,
     1 ' Topological radius R: ',I6,
     1 ', and average topological distance: ',F12.6)
C1004 Format(' Ori constant for Wiener index: ',D15.9)
 1006 Format(' Something wrong with Wiener sum')

      RETURN
      END

      SUBROUTINE Szeged(Edges,mdist,Sz)
      use config
C     This routine calculates the Szeged index
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer mdist(Nmax,Nmax),Edges(2,3*number_vertices/2)
      Sz=0.d0
C     Sum over all edges
      Do I=1,3*number_vertices/2
       IE1=Edges(1,I)+1
       IE2=Edges(2,I)+1
C      Get ni and nj for Szeged index
       ni=0
       nj=0
       do J=1,number_vertices
C       if(IE1.ne.J.and.IE2.ne.J) then
         if(mdist(IE1,J).lt.mdist(IE2,J)) ni=ni+1
         if(mdist(IE1,J).gt.mdist(IE2,J)) nj=nj+1
C       endif
       enddo
       Sz=Sz+dfloat(ni*nj)
      enddo
      RETURN
      END

      SUBROUTINE Chiral(Iout,Group)
      Character*3 Group
      Character*3 CPG(12)
      Data CPG/' C1',' C2',' C3',' C5',' C6',
     1 ' D2',' D3',' D5',' D5',' D6','  T','  I'/
      do I=1,12
       if(Group.eq.CPG(I)) then
        Write(Iout,1000) Group
       endif
      enddo
 1000 Format(1X,'Chiral fullerene detected belonging to',
     1 ' point group ',A3)
      RETURN
      END

      SUBROUTINE IcoFullDetect(Iout)
      use config
      IMPLICIT Integer (A-Z)
      N=number_vertices/20
      ntest=20*N
      ico=0
      icoh=0
      if(ntest.eq.number_vertices) then
       loopmax=int(sqrt(float(N)))
       loopmin=int(sqrt(float(N/3)))
       Do I=loopmin,loopmax
       Do J=0,I
        nico=20*(I*I+J*J+I*J)
        if(nico.eq.number_vertices) then
         if(J.eq.0.or.J.eq.I) then
          icoh=icoh+1
         else
          ico=ico+1
         endif
        endif
       enddo
       enddo
      endif
      icot=ico+icoh
      if(icot.ne.0) then
       if(icot.eq.1.and.icoh.eq.1) Write(Iout,1000) 
       if(icot.eq.1.and.ico.eq.1) Write(Iout,1001) 
       if(icot.gt.1) Write(Iout,1002) icot,icoh,ico
      else
       Write(Iout,1003)
      endif
 1000 Format(1X,'For this vertex number we can have 1 ',
     1 ' icoshedral fullerene of Ih-symmetry')
 1001 Format(1X,'For this vertex number we can have 1 ',
     1 ' icoshedral fullerene of I-symmetry')
 1002 Format(1X,'For this vertex number we can have ',I3,
     1 ' icoshedral fullerenes. Out of this ',I3,' are of',
     1 ' Ih-symmetry and ',I3,' of I-symmetry')
 1003 Format(1X,'There are no icosahedral fullerenes for ',
     1 'this vertex number')
      RETURN
      END

      SUBROUTINE Distan(IDA,Dist,Rmin,Rminall,Rmax,rms)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,NMAX)
      DIMENSION IDA(NMAX,NMAX)
C     Calculate minimal, maximal and root mean square distances
C     from adjacancy matrix IDA and cartesian coordinates Dist
      Rmin=1.d10
      Rminall=1.d10
      Rmax=0.d0
      Rrms=0.d0
      mc=0
      Do I=1,number_vertices
      Do J=I+1,number_vertices
      X=Dist(1,I)-Dist(1,J)
      Y=Dist(2,I)-Dist(2,J)
      Z=Dist(3,I)-Dist(3,J)
      R2=X*X+Y*Y+Z*Z
      R=dsqrt(R2)
      mc=mc+1
      if(R.lt.Rminall) Rminall=R
      if(IDA(I,J).ne.0) then
      if(R.lt.Rmin) Rmin=R
      if(R.gt.Rmax) Rmax=R
      Rrms=Rrms+R2
      endif
      enddo
      enddo
      if(Rmax.eq.0.d0.or.Rmin.eq.1.d10) then
      Print*,'**** Error in subroutine Distan'
      stop
      endif
      rms=dsqrt(Rrms/dfloat(mc))
      Return 
      END

      SUBROUTINE SortI(N,IS,JS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IS(6),JS(6)
C     Sort the N integer numbers, input IS, Output JS
      DO I=1,N
      JS(I)=IS(I)
      enddo
      N1=N-1
      DO 15 I=1,N1
      NAMIN=JS(I)
      M=I+1
      DO 25 K=M,N
      IF (NAMIN.LE.JS(K)) GO TO 25
      NZ=NAMIN
      NAMIN=JS(K)
      JS(K)=NZ
 25   CONTINUE
 15   JS(I)=NAMIN
      RETURN
      END

      SUBROUTINE Num2(IArray,I,J)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      I=IArray/number_vertices
      J=IArray-I*number_vertices
      If(J.eq.0) then
      J=number_vertices
      I=I-1
      endif
      If(I.eq.-1) then
      I=0
      J=0
      endif
      RETURN
      END

      SUBROUTINE DifDist(Ndif,N,Tol,Distac,Rmem)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Distac(6),Rmem(100)
      If(Ndif.eq.0) then
       Rmem(1)=Distac(1)
       Ndif=1
      endif
       do I=1,N
       do J=1,Ndif
       difR=dabs(Distac(I)-Rmem(J))
         if(difR.lt.Tol) Go to 1
       enddo
         Ndif=Ndif+1
         Rmem(Ndif)=Distac(I)
    1  Continue
       enddo
      Return
      END

      SUBROUTINE TIMER(TIMEX)
      Real TA(2)
       CALL DTIME(TA,time)
       TIMEX=TIME
      RETURN
      END

      SUBROUTINE FileMod(filenameIn,filenameOut,EndName,End,Ncyc,ifind)
      Implicit Integer (A-Z)
      CHARACTER*4  End
      CHARACTER*20 Number
      CHARACTER*7  EndName
      CHARACTER*50 filenameIn,filenameOut
C     Construct filename for xyz or cc1 files
C     If Ncyc > 1 then add number
        if(Ncyc.le.1) then
         filenameOut=trim(filenameIn)//EndName
        else
         write(Number,*) Ncyc
         filenameOut=trim(filenameIn)//trim(adjustl(Number))//End
        endif
        ichar1=index(filenameOut,'database/ALL')
        ichar2=index(filenameOut,'database/IPR')
        ichar3=index(filenameOut,'database/Yoshida')
        ichar4=index(filenameOut,'database/HOG')
        ifind=ichar1+ichar2+ichar3+ichar4
      RETURN
      END


