      SUBROUTINE SpiralRestart(N,IPR,Iout,
C     Restart version of Subroutine Spiral
     1 Isonum,IsoIPR,iham,IDA,A,chkname)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),IDA(NMAX,NMAX)
      Real*8 A(NMAX,NMAX),gap
      Integer Isonum(119),IsoIPR(123)
      DIMENSION NMR(6),IRhag5(0:5),IRhag6(0:6)
      DIMENSION Spiralx(12,NMAX)
      CHARACTER*1 DUMMY
      CHARACTER*3 GROUP
      CHARACTER*6 Occup
      CHARACTER*50 chkname
      CHARACTER*18 Start,Char,Last
      Real*8 sigmah,sigmahlow,sigmahhigh
      Logical lexist

C Set parameters
      Start=' Isomer List Start'
      Last=' Isomer List Compl'
      Nloop=1000000000
      skip=4
      nhamcycle=0
      IPRdect=0
      maxiter=10000000
      islow=1
      islowIPR=1
      ishigh=1
      ishighIPR=1
      IFusL=1
      IFusH=1
      ISigmaL=1
      ISigmaH=1

C Test if file is the right one
      Write(Iout,1000) chkname
         inquire(file=chkname,exist=lexist)
          if(lexist.neqv..True.) then
            Write(Iout,1023) chkname
            stop
          endif
      Open(UNIT=1,FILE=chkname,ERR=500,STATUS='old',FORM='FORMATTED')
      do I=1,Nloop
      Read(1,2002,end=100) Char
       if(Char.eq.Start) go to 200
      enddo
  100  Write(Iout,1001) I,Start
       Close(unit=1)
       stop
  200 Write(Iout,1002)

C Test if atom number is right
      Read(1,*) IN,IP,IH
      if(IN.ne.N) then
       Write(Iout,1003) N,IN
       Close(unit=1)
       stop
      endif
      if(IH.ne.0) then
       Write(Iout,1004)
        iham=1
       else
        iham=0
      endif
      if(IP.eq.0) then
       Write(Iout,1005)
       IPR=0
      else
       Write(Iout,1006)
       IPR=1
      endif

C Now loop over all data to find the last one
C Check parameters as well
      do I=1,skip
       Read(1,*,end=300) Dummy
      enddo
C Without Hamiltonian cycles
      if(iham.eq.0) then
       Read(1,2000,Err=300) L,GROUP,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,
     1 K11,K12,IFus5G,sigmah
      IFus5Ghigh=IFus5G
      IFus5Glow=IFus5G
      sigmahhigh=sigmah
      sigmahlow=sigmah
      if(L.ne.1) Go to 300
      do I=2,Nloop
       Read(1,2000,Err=400,end=400) L,GROUP,K1,K2,K3,K4,K5,K6,K7,K8,
     1 K9,K10,K11,K12,IFus5G,sigmah
       if(IFus5G.le.IFus5Glow) then
        IFus5Glow=IFus5G
        IFusL=L
       endif
       if(IFus5G.ge.IFus5Ghigh) then
        IFus5Ghigh=IFus5G
        IFusH=L
       endif
       if(sigmah.le.sigmahlow) then
        sigmahlow=sigmah
        ISigmaL=L
       endif
       if(sigmah.ge.sigmahhigh) then
        sigmahhigh=sigmah
        ISigmaH=L
       endif
      enddo
      endif
C With Hamiltonian cycles
      if(iham.ne.0) then
       Read(1,2001,Err=300) L,GROUP,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,
     1 K11,K12,IFus5G,sigmah,nhamcycle
      IFus5Ghigh=IFus5G
      IFus5Glow=IFus5G
      sigmahhigh=sigmah
      sigmahlow=sigmah
      IPRhamlow=nhamcycle
      IPRhamhigh=nhamcycle
      hamlow=nhamcycle
      hamhigh=nhamcycle
      if(L.ne.1) Go to 300
      do I=2,Nloop
       Read(1,2001,Err=400,end=400) L,GROUP,K1,K2,K3,K4,K5,K6,K7,K8,
     1 K9,K10,K11,K12,IFus5G,sigmah,nhamcycle
       if(IFus5G.le.IFus5Glow) then
        IFus5Glow=IFus5G
        IFusL=L
       endif
       if(IFus5G.ge.IFus5Ghigh) then
        IFus5Ghigh=IFus5G
        IFusH=L
       endif
       if(sigmah.le.sigmahlow) then
        sigmahlow=sigmah
        ISigmaL=L
       endif
       if(sigmah.ge.sigmahhigh) then
        sigmahhigh=sigmah
        ISigmaH=L
       endif
        if(nhamcycle.le.hamlow) then
         hamlow=nhamcycle
         islow=L
        endif
        if(nhamcycle.ge.hamhigh) then
         hamhigh=nhamcycle
         ishigh=L
        endif
        if(IPR.eq.0.and.IFus5G.eq.0) then
         IPRdect=1
         if(nhamcycle.lt.IPRhamlow) then
          IPRhamlow=nhamcycle
          islowIPR=L
         endif
         if(nhamcycle.gt.IPRhamhigh) then
          IPRhamhigh=nhamcycle
          ishighIPR=L
         endif
        endif
      enddo
      endif
  400  Write(Iout,1008) L,GROUP,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12
      if(L.eq.0) then
       WRITE (Iout,1015)
       return
      endif
       Close(unit=1)
      if(iham.ne.0) then
       WRITE (Iout,1010) hamlow,islow,hamhigh,ishigh
       if(IPR.eq.0.and.IPRdect.eq.1) then
        WRITE (Iout,1011) IPRhamlow,islowIPR,IPRhamhigh,ishighIPR
       endif
      endif
      WRITE (Iout,1012) IFus5Glow,IFusL,IFus5Ghigh,IFusH,
     1  sigmahlow,ISigmaL,sigmahhigh,ISigmaH

C Now check if list is complete or not
       M1=N/2-9
       M2=N/2-29
       if(IPR.eq.0) then
        if(M1.le.119.and.L.ge.Isonum(M1)) then
         Write(Iout,1009) Isonum(M1)
         return
        endif
       else
        if(M2.le.123.and.M2.gt.0.and.L.ge.IsoIPR(M2)) then
         Write(Iout,1009) IsoIPR(M2)
         return
        endif
       endif

C Now do the calculation for the remainder
      Write(Iout,1013) 
      if(iham.eq.0) then
      if(IPR.EQ.0) then 
         IF(N.lt.100) WRITE(Iout,601) N
         IF(N.ge.100) WRITE(Iout,602) N
      else
         IPR=1 
         IF(N.lt.100) WRITE(Iout,603) N
         IF(N.ge.100) WRITE(Iout,604) N
      endif
      else
      if(IPR.EQ.0) then 
         IF(N.lt.100) WRITE(Iout,701) N
         IF(N.ge.100) WRITE(Iout,702) N
      else
         IPR=1 
         IF(N.lt.100) WRITE(Iout,703) N
         IF(N.ge.100) WRITE(Iout,704) N
      endif
      endif

      do I=1,MMAX
      do J=1,MMAX
      D(I,J)=0
      enddo
      enddo
      IER=0
      IT=0
      JPR=IPR+1
      M=N/2+2
      itest=1
      DO 1  J1= 1    ,M-11*JPR !   Open loop over spiral
       if(J1.lt.K1.and.itest.eq.1) go to 1
      DO 2  J2=J1+JPR,M-10*JPR !   combinations
       if(J2.lt.K2.and.itest.eq.1) go to 2
      DO 3  J3=J2+JPR,M-9*JPR
       if(J3.lt.K3.and.itest.eq.1) go to 3
      DO 4  J4=J3+JPR,M-8*JPR
       if(J4.lt.K4.and.itest.eq.1) go to 4
      DO 5  J5=J4+JPR,M-7*JPR
       if(J5.lt.K5.and.itest.eq.1) go to 5
      DO 6  J6=J5+JPR,M-6*JPR
       if(J6.lt.K6.and.itest.eq.1) go to 6
      DO 7  J7=J6+JPR,M-5*JPR
       if(J7.lt.K7.and.itest.eq.1) go to 7
      DO 8  J8=J7+JPR,M-4*JPR
       if(J8.lt.K8.and.itest.eq.1) go to 8
      DO 9  J9=J8+JPR,M-3*JPR
       if(J9.lt.K9.and.itest.eq.1) go to 9
      DO 10 J10=J9+JPR,M-2*JPR
       if(J10.lt.K10.and.itest.eq.1) go to 10
      DO 11 J11=J10+JPR,M-JPR
       if(J11.lt.K11.and.itest.eq.1) go to 11
      DO 12 J12=J11+JPR,M
       if(J12.lt.K12.and.itest.eq.1) go to 12
       if(J12.eq.K12.and.itest.eq.1) then
        itest=0
        go to 12
       endif
        DO J=1,M               ! Form spiral code in S
         S(J)=6
        enddo
      S(J1)=5
      S(J2)=5
      S(J3)=5
      S(J4)=5
      S(J5)=5
      S(J6)=5
      S(J7)=5
      S(J8)=5
      S(J9)=5
      S(J10)=5
      S(J11)=5
      S(J12)=5
      CALL Windup(M,IPR,IER,S,D)      !      Wind up spiral into dual 
      IF(IER.EQ.12) GO TO 12               !      and check for closure 
      IF(IER.EQ.11) GO TO 11
      IF(IER.EQ.10) GO TO 10
      IF(IER.EQ.9)  GO TO 9
      IF(IER.EQ.8)  GO TO 8
      IF(IER.EQ.7)  GO TO 7
      IF(IER.EQ.6)  GO TO 6
      IF(IER.EQ.5)  GO TO 5
      IF(IER.EQ.4)  GO TO 4
      IF(IER.EQ.3)  GO TO 3
      IF(IER.EQ.2)  GO TO 2
      IF(IER.EQ.1)  GO TO 1
      CALL Unwind(M,IER,IT,ispiral,
     1 Spiralx,S,D,NMR,Group)                            ! Unwind dual into spirals 
      IF(IER.EQ.13) GO TO 13                             ! and check for uniqueness
      K=0
      L=L+1                                              ! Spiral S is canonical      
      DO J=1,6
         IF(NMR(J).EQ.0) GO TO 16
         K=J
      enddo
C     Analyze dual matrix
   16  CALL DualAnalyze(nmax,mmax,N,M,D,IRhag5,IRhag6,
     1 IFus5G,IDA,nelec,ndeg,sigmah,A,gap)
       if(2*ndeg.eq.nelec) then 
        Occup='closed'
       else
        Occup='open  '
        gap=0.d0
       endif
       if(IFus5G.le.IFus5Glow) then
        IFus5Glow=IFus5G
        IFusL=L
       endif
       if(IFus5G.ge.IFus5Ghigh) then
        IFus5Ghigh=IFus5G
        IFusH=L
       endif
       if(sigmah.le.sigmahlow) then
        sigmahlow=sigmah
        ISigmaL=L
       endif
       if(sigmah.ge.sigmahhigh) then
        sigmahhigh=sigmah
        ISigmaH=L
       endif
       if(iham.eq.0) then
       WRITE(Iout,607) L,GROUP,J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,
     1 (IRhag5(J),J=0,5),IFus5G,(IRhag6(J),J=0,6),sigmah,
     2 nelec,ndeg,gap,Occup,(NMR(J),J=1,K)
       else
        Call HamiltonCyc(N,maxiter,Iout,nbatch,IDA,nhamcycle)
        WRITE(Iout,608) L,GROUP,J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,
     1   (IRhag5(J),J=0,5),IFus5G,(IRhag6(J),J=0,6),sigmah,
     2   nelec,ndeg,gap,Occup,nhamcycle,(NMR(J),J=1,K)
        if(nhamcycle.le.hamlow) then
         hamlow=nhamcycle
         islow=L
        endif
        if(nhamcycle.ge.hamhigh) then
         hamhigh=nhamcycle
         ishigh=L
        endif
        if(IPR.eq.0.and.IFus5G.eq.0) then
         IPRdect=1
         if(nhamcycle.lt.IPRhamlow) then
          IPRhamlow=nhamcycle
          islowIPR=L
         endif
         if(nhamcycle.gt.IPRhamhigh) then
          IPRhamhigh=nhamcycle
          ishighIPR=L
         endif
        endif
       endif
       if(IPR.eq.0) then
        if(M1.le.119.and.L.ge.Isonum(M1)) go to 99
       else
        if(M2.le.123.and.M2.gt.0.and.L.ge.IsoIPR(M2)) go to 99
      endif
 13     CONTINUE 
 12     CONTINUE        ! Close loop over spiral 
 11     CONTINUE        ! combinations      
 10     CONTINUE
 9      CONTINUE
 8      CONTINUE
 7      CONTINUE
 6      CONTINUE
 5      CONTINUE
 4      CONTINUE
 3      CONTINUE
 2      CONTINUE
 1      CONTINUE
 99   WRITE (Iout,612)
      if(iham.ne.0) then
       WRITE (Iout,609) hamlow,islow,hamhigh,ishigh
       if(IPR.eq.0.and.IPRdect.eq.1) then
        WRITE (Iout,610) IPRhamlow,islowIPR,IPRhamhigh,ishighIPR
       endif
      endif
      if(IPR.eq.0) then
       WRITE (Iout,611) IFus5Glow,IFusL,IFus5Ghigh,IFusH,
     1  sigmahlow,ISigmaL,sigmahhigh,ISigmaH
      else
       WRITE (Iout,613) sigmahlow,ISigmaL,sigmahhigh,ISigmaH
      endif
      WRITE (Iout,606)
      Return

  300 Write(Iout,1007)
      Close(unit=1)
      Return

  500 Write(Iout,1014) chkname
      stop

 601  FORMAT(1X,'General fullerene isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 602  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 603  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 604  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 606  FORMAT(/1X,170('-'),/1X,'End of subroutine Spiral')
 607  FORMAT(1X,I8,2X,A3,1X,12I4,2X,'(',5(I2,','),I2,')  ',I2,
     1 2X,'(',6(I2,','),I3,')  ',F8.5,2X,I2,1X,I2,1X,F8.5,
     1 1X,A6,2X,3(I3,' x',I3,:,','))
 608  FORMAT(1X,I8,2X,A3,1X,12I4,2X,'(',5(I2,','),I2,')  ',I2,
     1 2X,'(',6(I2,','),I3,')  ',F8.5,2X,I2,1X,I2,1X,F8.5,
     1 1X,A6,2X,I9,2X,3(I3,' x',I3,:,','))
 609  FORMAT(1X,'Lowest number of Hamiltonian cycles     ',I10,
     1 ' for isomer ',I10,
     1      /1X,'Highest number of Hamiltonian cycles    ',I10,
     1 ' for isomer ',I10)
 610  FORMAT(1X,'Lowest number of IPR Hamiltonian cycles ',I10,
     1 ' for isomer ',I10,
     1      /1X,'Highest number of IPR Hamiltonian cycles',I10,
     1 ' for isomer ',I10)
 611  FORMAT(1X,'Lowest  Np= ',I3,' for isomer ',I10,
     1      /1X,'Highest Np= ',I3,' for isomer ',I10,
     1      /1X,'Lowest  Sigmah= ',F8.5,' for isomer ',I10,
     1      /1X,'Highest Sigmah= ',F8.5,' for isomer ',I10)
 612  FORMAT(1X,'Isomer List Complete')
 613  FORMAT(1X,'Lowest  Sigmah= ',F8.5,' for isomer ',I10,
     1      /1X,'Highest Sigmah= ',F8.5,' for isomer ',I10)
 701  FORMAT(1X,'General fullerene isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 702  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 703  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 704  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 1000 Format(1X,'Opening file: ',A20)
 1001 Format(1X,'Lopped over ',I9,' lines. Characterestic card ',
     1 'for isomer list not found: ',A18,' ==> ABORT')
 1002 Format(1X,'Reading input')
 1003 Format(1X,'Atom number ',I5,' requested not identical to ',
     1 I5,' in the data file  ==> ABORT')
 1004 Format(1X,'File contains number of distinct Hamiltonian ',
     1 'cycles')
 1005 Format(1X,'All isomers considered')
 1006 Format(1X,'Only IPR isomers considered')
 1007 Format(1X,'First Data card not found ==> ABORT')
 1008 Format(1X,'Last isomer in list: ',1X,I8,2X,A3,1X,12I4,
     1 /1X,'General statistics from external file data:')
 1009 Format(1X,'List has ',I9,' isomers and is therefore complete',
     1 /1X,'Nothing else to do --> Return')
 1010 FORMAT(1X,'Lowest number of Hamiltonian cycles     ',I10,
     1 ' for isomer ',I10,
     1      /1X,'Highest number of Hamiltonian cycles    ',I10,
     1 ' for isomer ',I10)
 1011 FORMAT(1X,'Lowest number of IPR Hamiltonian cycles ',I10,
     1 ' for isomer ',I10,
     1      /1X,'Highest number of IPR Hamiltonian cycles',I10,
     1 ' for isomer ',I10)
 1012 FORMAT(1X,'Lowest  Np= ',I3,' for isomer ',I10,
     1      /1X,'Highest Np= ',I3,' for isomer ',I10,
     1      /1X,'Lowest  Sigmah= ',F8.5,' for isomer ',I10,
     1      /1X,'Highest Sigmah= ',F8.5,' for isomer ',I10)
 1013 FORMAT(/1X,'Starting with new list:')
 1014 Format(1X,'File: ',A20,' does not exist ==> ABORT')
 1015 Format(1X,'Last entry is zero, so check that you do not have ',
     1 'an empty line at the end of the input')
 1023 Format(1X,'Filename ',A50,' in database not found ==> ABORT')
 2000 Format(I9,2X,A3,1X,12I4,23X,I2,27X,F8.5)
 2001 Format(I9,2X,A3,1X,12I4,23X,I2,27X,F8.5,25X,I9)
 2002 Format(A18)
      END

      SUBROUTINE SpiralIco(N,ivar,In,Iout,IDA,A)
C     This subroutine comes directly from the book of Fowler and 
C     Manolopoulos "An Atlas of Fullerenes" (Dover Publ., New York, 2006),
C     and has been modified to search for ring spirals around an
C     icosahedral fullerene clostest to the vertex number required.           
C     This sub-program catalogues fullerenes with a given number of      
C     vertices using the spiral algorithm and a uniqueness test 
C     based on equivalent spirals. IPR = 1  for isolated-pentagon isomers. 
C     The resulting output is a catalogue of the isomers found containing 
C     their idealized point groups, canonical spirals, and NMR patterns.        
C     N is the nuclearity of the fullerene.
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),SS(MMAX),IDA(NMAX,NMAX)
      DIMENSION NMR(6),IRhag5(0:5),IRhag6(0:6)
      DIMENSION Spiralx(12,NMAX),IVL(12),IVH(12)
      CHARACTER*3 GROUP
      CHARACTER*6 Occup
      Real*8 sigmah,sigmahlow,sigmahhigh
      Real*8 A(NMAX,NMAX),gap,dex

      IPR=1
      ivarlimit=4
      if(ivar.gt.ivarlimit) then
       Write(Iout,805)
       Return
      endif
      if(N.gt.NMAX) Return     ! Increase NMAX 
      if(2*(N/2).ne.N) Return  ! N must be even 
C   Search for icosahedral fullerene Cm closest to Cm with m<n
      dex=dfloat(N)
      icon=0
      loopmax=int(dsqrt(dex/2.d1))
      loopmin=int(dsqrt(dex/6.d1))
      if(loopmin.lt.1) loopmin=1
       Do I1=loopmin,loopmax
       Do J1=0,I1
        nico=20*(I1*I1+J1*J1+I1*J1)
        if(nico.le.N.and.nico.gt.icon) then
         icon=nico
         I=I1
         J=J1
        endif
       enddo
       enddo
       if(icon.eq.0) then
        Write(Iout,800) 
        stop
       endif
       if(icon.eq.N) then
        Write(Iout,801) icon,I,J
C       Return
       else
        Write(Iout,802) icon,N,I,J
       endif

      Read(IN,*,Err=400,end=400) (SS(I),I=1,12)
      Go to 300

  400 I2=I*I
      J2=J*J
C     Getting exponents
      IA=(5*(I+J)**2-5*I-3*J-2)/2
      IB=I+J-1
      IC=(5*I+1)*(I-1)+J*(5*I-3)
      ID=IB
      IE=(5*I2+15*J2-3*I-7*J)/2
      WRITE(Iout,803) icon,IA,IB,IC,ID,IE
C     Construct the ring spiral
      M=N/2+2
      SS(1)=1
      SS(2)=IA+2
      SS(3)=SS(2)+IB+1
      SS(4)=SS(3)+IB+1
      SS(5)=SS(4)+IB+1
      SS(6)=SS(5)+IB+1
      SS(7)=SS(6)+IC+1
      SS(8)=SS(7)+ID+1
      SS(9)=SS(8)+ID+1
      SS(10)=SS(9)+ID+1
      SS(11)=SS(10)+ID+1
      SS(12)=SS(11)+IE+1
      Write(Iout,804) (SS(I),I=1,12),ivar,ivar

C  Set loop limits
  300 IVL(1)=SS(1)
      IVH(1)=SS(1)+1
      do I=2,12
       IVL(I)=SS(I)-ivar
       IVH(I)=SS(I)+ivar
      enddo

C  Set parameters
      IRSPI=0
      maxiter=10000000
      islow=0
      islowIPR=0
      ishigh=0
      ishighIPR=0
      IFus5Glow=1000
      sigmahlow=1.d10
      IFus5Ghigh=0
      sigmahhigh=-1.d0
      L=0
      IER=0
      IT=0
      do I=1,MMAX
      do J=1,MMAX
       D(I,J)=0
      enddo
      enddo

      Write(iout,600)

C  Search for clashes in RSPIs, e.g. SS(j).le.SS(j-1)
C     if(IVH(12).gt.M) IVH(12)=M
C     do I=1,11
C      ndif=IVL(I+1)-IVH(I)
C      if(ndif.le.0) then
C       Write(Iout,806) I,I+1
C       Return
C      endif
C     enddo

      DO 1  J1=IVL(1),IVH(1)
       IVL2=IVL(2)
       IVH2=IVH(2)
       if(J1.ge.IVL2) IVL2=J1+1
       if(IVL2.gt.IVH2) IVH2=IVL2
       if(IVH2.gt.M) Return
      DO 2  J2=IVL2,IVH2         !   combinations
       IVL3=IVL(3)
       IVH3=IVH(3)
       if(J2.ge.IVL3) IVL3=J2+1
       if(IVL3.gt.IVH3) IVH3=IVL3
       if(IVH3.gt.M) Return
      DO 3  J3=IVL3,IVH3
       IVL4=IVL(4)
       IVH4=IVH(4)
       if(J3.ge.IVL4) IVL4=J3+1
       if(IVL4.gt.IVH4) IVH4=IVL4
       if(IVH4.gt.M) Return
      DO 4  J4=IVL4,IVH4
       IVL5=IVL(5)
       IVH5=IVH(5)
       if(J4.ge.IVL5) IVL5=J4+1
       if(IVL5.gt.IVH5) IVH5=IVL5
       if(IVH5.gt.M) Return
      DO 5  J5=IVL5,IVH5
       IVL6=IVL(6)
       IVH6=IVH(6)
       if(J5.ge.IVL6) IVL6=J5+1
       if(IVL6.gt.IVH6) IVH6=IVL6
       if(IVH6.gt.M) Return
      DO 6  J6=IVL6,IVH6
       IVL7=IVL(7)
       IVH7=IVH(7)
       if(J6.ge.IVL7) IVL7=J6+1
       if(IVL7.gt.IVH7) IVH7=IVL7
       if(IVH7.gt.M) Return
      DO 7  J7=IVL7,IVH7
       IVL8=IVL(8)
       IVH8=IVH(8)
       if(J7.ge.IVL8) IVL8=J7+1
       if(IVL8.gt.IVH8) IVH8=IVL8
       if(IVH8.gt.M) Return
      DO 8  J8=IVL8,IVH8
       IVL9=IVL(9)
       IVH9=IVH(9)
       if(J8.ge.IVL9) IVL9=J8+1
       if(IVL9.gt.IVH9) IVH9=IVL9
       if(IVH9.gt.M) Return
      DO 9  J9=IVL9,IVH9
       IVL10=IVL(10)
       IVH10=IVH(10)
       if(J9.ge.IVL10) IVL10=J9+1
       if(IVL10.gt.IVH10) IVH10=IVL10
       if(IVH10.gt.M) Return
      DO 10 J10=IVL10,IVH10
       IVL11=IVL(11)
       IVH11=IVH(11)
       if(J10.ge.IVL11) IVL11=J10+1
       if(IVL11.gt.IVH11) IVH11=IVL11
       if(IVH11.gt.M) Return
      DO 11 J11=IVL11,IVH11
       IVL12=IVL(12)
       IVH12=IVH(12)
       if(J11.ge.IVL12) IVL12=J11+1
       if(IVL12.gt.IVH12) IVH12=IVL12
       if(IVH12.gt.M) Return
      DO 12 J12=IVL12,IVH12

      DO J=1,M   ! Form spiral code in S
       S(J)=6
      enddo
      S(J1)=5
      S(J2)=5
      S(J3)=5
      S(J4)=5
      S(J5)=5
      S(J6)=5
      S(J7)=5
      S(J8)=5
      S(J9)=5
      S(J10)=5
      S(J11)=5
      S(J12)=5
      CALL Windup(M,IPR,IER,S,D)      !      Wind up spiral into dual 
      IF(IER.EQ.12) GO TO 12               !      and check for closure 
      IF(IER.EQ.11) GO TO 11
      IF(IER.EQ.10) GO TO 10
      IF(IER.EQ.9)  GO TO 9
      IF(IER.EQ.8)  GO TO 8
      IF(IER.EQ.7)  GO TO 7
      IF(IER.EQ.6)  GO TO 6
      IF(IER.EQ.5)  GO TO 5
      IF(IER.EQ.4)  GO TO 4
      IF(IER.EQ.3)  GO TO 3
      IF(IER.EQ.2)  GO TO 2
      IF(IER.EQ.1)  GO TO 1
      CALL Unwind(M,IER,IT,ispiral,
     1 Spiralx,S,D,NMR,Group)                            ! Unwind dual into spirals 
      IF(IER.EQ.13) GO TO 13                             ! and check for uniqueness      
      K=0
      L=L+1                                              ! Spiral S is canonical      
      DO J=1,6
         IF(NMR(J).EQ.0) GO TO 16
         K=J
      enddo
C     Analyze dual matrix
   16  CALL DualAnalyze(nmax,mmax,N,M,D,IRhag5,IRhag6,
     1 IFus5G,IDA,nelec,ndeg,sigmah,A,gap)
       if(2*ndeg.eq.nelec) then 
        Occup='closed'
       else
        Occup='open  '
        gap=0.d0
       endif
       if(IFus5G.le.IFus5Glow) then
        IFus5Glow=IFus5G
        IFusL=L
       endif
       if(IFus5G.ge.IFus5Ghigh) then
        IFus5Ghigh=IFus5G
        IFusH=L
       endif
       if(sigmah.le.sigmahlow) then
        sigmahlow=sigmah
        ISigmaL=L
       endif
       if(sigmah.ge.sigmahhigh) then
        sigmahhigh=sigmah
        ISigmaH=L
       endif
       WRITE(Iout,607) L,GROUP,J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,
     1 (IRhag5(J),J=0,5),IFus5G,(IRhag6(J),J=0,6),sigmah,
     2 nelec,ndeg,gap,Occup,(NMR(J),J=1,K)
        IRSPI=IRSPI+1
 13     CONTINUE 
 12     CONTINUE        ! Close loop over spiral 
 11     CONTINUE        ! combinations      
 10     CONTINUE
 9      CONTINUE
 8      CONTINUE
 7      CONTINUE
 6      CONTINUE
 5      CONTINUE
 4      CONTINUE
 3      CONTINUE
 2      CONTINUE
 1      CONTINUE

      if(IRSPI.ne.0) then
       WRITE (Iout,611) IFus5Glow,IFusL,IFus5Ghigh,IFusH,
     1  sigmahlow,ISigmaL,sigmahhigh,ISigmaH
       WRITE (Iout,612)
      else
       WRITE (Iout,613)
      endif
      WRITE (Iout,606)
 600  FORMAT(/1X,'Subroutine Spiral from Fowler and Manolopoulos',
     1 ' (An Atlas of Fullerenes, Dover Publ., New York, 2006)',
     2 /1X,'(Symmetries are given for undistorted fullerenes)',
     1 /1X,' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 606  FORMAT(/1X,170('-'),/1X,'End of subroutine Spiral')
 607  FORMAT(1X,I8,2X,A3,1X,12I4,2X,'(',5(I2,','),I2,')  ',I2,
     1 2X,'(',6(I2,','),I3,')  ',F8.5,2X,I2,1X,I2,1X,F8.5,
     1 1X,A6,2X,3(I3,' x',I3,:,','))
 611  FORMAT(1X,'Lowest  Np= ',I3,' for isomer ',I10,
     1      /1X,'Highest Np= ',I3,' for isomer ',I10,
     1      /1X,'Lowest  Sigmah= ',F8.5,' for isomer ',I10,
     1      /1X,'Highest Sigmah= ',F8.5,' for isomer ',I10)
 612  FORMAT(1X,'Isomer List Complete')
 613  FORMAT(1X,'No isomers found')
 800  FORMAT(1X,'Nothing found ==> STOP')
 801  FORMAT(1X,'Your vertex number matches icoshedral fullerene ',
     1 'with N=',I6,' ,(k,l)=(',I2,',',I2,')')
 802  FORMAT(1X,'Closest vertex number is ',I6,' for icoshedral ',
     1 'fullerene with N=',I6,'; (k,l)=(',I2,',',I2,')')
 803  Format(1X,'General Goldberg-Coxeter transformation of C20 -> Cn',
     1 ' with n=: ',I5,/1X,'Construction of icosahedral fullerenes ',
     1 'using the spiral code representation of Fowler and Rogers',/1X,
     1 ' Lit: P. W. Fowler, K. M. Rogers, J. Chem. Inf. Comput. Sci.',
     1 ' 41, 108-111 (2001)',/1X,'Exponents of spiral series ',
     1 '5(6)^A (5(6)^B)^4 5(6)^C (5(6)^D)^4 5(6)^E 5:',
     1 '  A=',I5,', B=',I5,', C=',I5,', D=',I5,', E=',I5)
 804  Format(1X,'Ring spiral pentagon indices RSPI: ',12I5,
     1 /1X,'Variation of each RSPI from RSPI-',I1,' to ','RSPI+',I1)
 805  Format(1X,'Variation too large, too many loops ==> RETURN')
C806  Format(1X,'Overlaping RSPIs between ',I2,' and ',I2,' ==> RETURN')

C     RETURN
      END

      SUBROUTINE Spiral(N,IPR,Iout,
     1 Isonum,IsoIPR,iham,IDA,A)
C     This subroutine comes directly from the book of Fowler and 
C     Manolopoulos "An Atlas of Fullerenes" (Dover Publ., New York, 2006).           
C     This sub-program catalogues fullerenes with a given number of      
C     vertices using the spiral algorithm and a uniqueness test 
C     based on equivalent spirals. The required input is IPR, 
C     where IPR = 0 for general and 1 for isolated-pentagon isomers. 
C     The resulting output is a catalogue of the isomers found containing 
C     their idealized point groups, canonical spirals, and NMR patterns.        
C     N is the nuclearity of the fullerene.
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),IDA(NMAX,NMAX)
      Real*8 A(NMAX,NMAX),gap
      Integer Isonum(119),IsoIPR(123)
      DIMENSION NMR(6),IRhag5(0:5),IRhag6(0:6)
      DIMENSION Spiralx(12,NMAX)
      CHARACTER*3 GROUP
      CHARACTER*6 Occup
      Real*8 sigmah,sigmahlow,sigmahhigh

      if(N.gt.NMAX) Return            ! Increase NMAX 
      if(2*(N/2).ne.N) Return         ! N must be even 
      M1=N/2-9
      M2=N/2-29
      nhamcycle=0
      IPRdect=0
      maxiter=10000000
      hamlow=2000000000
      IPRhamlow=2000000000
      hamhigh=0
      IPRhamhigh=0
      islow=0
      islowIPR=0
      ishigh=0
      ishighIPR=0
      IPRhamlow=2000000000
      IPRhamhigh=-1
      IFus5Glow=1000
      sigmahlow=1.d10
      IFus5Ghigh=0
      sigmahhigh=-1.d0
      WRITE (Iout,600) N,IPR,iham
      if(iham.eq.0) then
      if(IPR.EQ.0) then 
         IF(N.lt.100) WRITE(Iout,601) N
         IF(N.ge.100) WRITE(Iout,602) N
      else
         IPR=1 
         IF(N.lt.100) WRITE(Iout,603) N
         IF(N.ge.100) WRITE(Iout,604) N
      endif
      else
      if(IPR.EQ.0) then 
         IF(N.lt.100) WRITE(Iout,701) N
         IF(N.ge.100) WRITE(Iout,702) N
      else
         IPR=1 
         IF(N.lt.100) WRITE(Iout,703) N
         IF(N.ge.100) WRITE(Iout,704) N
      endif
      endif
      do I=1,MMAX
      do J=1,MMAX
       D(I,J)=0
      enddo
      enddo
      L=0
      IER=0
      IT=0
      JPR=IPR+1
      M=N/2+2
      DO 1  J1= 1    ,M-11*JPR !   Open loop over spiral
      DO 2  J2=J1+JPR,M-10*JPR !   combinations
      DO 3  J3=J2+JPR,M-9*JPR
      DO 4  J4=J3+JPR,M-8*JPR
      DO 5  J5=J4+JPR,M-7*JPR
      DO 6  J6=J5+JPR,M-6*JPR
      DO 7  J7=J6+JPR,M-5*JPR
      DO 8  J8=J7+JPR,M-4*JPR
      DO 9  J9=J8+JPR,M-3*JPR
      DO 10 J10=J9+JPR,M-2*JPR
      DO 11 J11=J10+JPR,M-JPR
      DO 12 J12=J11+JPR,M
        DO J=1,M               ! Form spiral code in S
         S(J)=6
        enddo
      S(J1)=5
      S(J2)=5
      S(J3)=5
      S(J4)=5
      S(J5)=5
      S(J6)=5
      S(J7)=5
      S(J8)=5
      S(J9)=5
      S(J10)=5
      S(J11)=5
      S(J12)=5
      CALL Windup(M,IPR,IER,S,D)      !      Wind up spiral into dual 
      IF(IER.EQ.12) GO TO 12               !      and check for closure 
      IF(IER.EQ.11) GO TO 11
      IF(IER.EQ.10) GO TO 10
      IF(IER.EQ.9)  GO TO 9
      IF(IER.EQ.8)  GO TO 8
      IF(IER.EQ.7)  GO TO 7
      IF(IER.EQ.6)  GO TO 6
      IF(IER.EQ.5)  GO TO 5
      IF(IER.EQ.4)  GO TO 4
      IF(IER.EQ.3)  GO TO 3
      IF(IER.EQ.2)  GO TO 2
      IF(IER.EQ.1)  GO TO 1
      CALL Unwind(M,IER,IT,ispiral,
     1 Spiralx,S,D,NMR,Group)                            ! Unwind dual into spirals 
      IF(IER.EQ.13) GO TO 13                             ! and check for uniqueness      
      K=0
      L=L+1                                              ! Spiral S is canonical      
      DO J=1,6
         IF(NMR(J).EQ.0) GO TO 16
         K=J
      enddo
C     Analyze dual matrix
   16  CALL DualAnalyze(nmax,mmax,N,M,D,IRhag5,IRhag6,
     1 IFus5G,IDA,nelec,ndeg,sigmah,A,gap)
       if(2*ndeg.eq.nelec) then 
        Occup='closed'
       else
        Occup='open  '
        gap=0.d0
       endif
       if(IFus5G.le.IFus5Glow) then
        IFus5Glow=IFus5G
        IFusL=L
       endif
       if(IFus5G.ge.IFus5Ghigh) then
        IFus5Ghigh=IFus5G
        IFusH=L
       endif
       if(sigmah.le.sigmahlow) then
        sigmahlow=sigmah
        ISigmaL=L
       endif
       if(sigmah.ge.sigmahhigh) then
        sigmahhigh=sigmah
        ISigmaH=L
       endif
       if(iham.eq.0) then
       WRITE(Iout,607) L,GROUP,J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,
     1 (IRhag5(J),J=0,5),IFus5G,(IRhag6(J),J=0,6),sigmah,
     2 nelec,ndeg,gap,Occup,(NMR(J),J=1,K)
       else
        Call HamiltonCyc(N,maxiter,Iout,nbatch,IDA,nhamcycle)
        WRITE(Iout,608) L,GROUP,J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,
     1   (IRhag5(J),J=0,5),IFus5G,(IRhag6(J),J=0,6),sigmah,
     2   nelec,ndeg,gap,Occup,nhamcycle,(NMR(J),J=1,K)
        if(nhamcycle.le.hamlow) then
         hamlow=nhamcycle
         islow=L
        endif
        if(nhamcycle.ge.hamhigh) then
         hamhigh=nhamcycle
         ishigh=L
        endif
        if(IPR.eq.0.and.IFus5G.eq.0) then
         IPRdect=1
         if(nhamcycle.lt.IPRhamlow) then
          IPRhamlow=nhamcycle
          islowIPR=L
         endif
         if(nhamcycle.gt.IPRhamhigh) then
          IPRhamhigh=nhamcycle
          ishighIPR=L
         endif
        endif
       endif
       if(IPR.eq.0) then
        if(M1.le.119.and.L.eq.Isonum(M1)) go to 99
       else
        if(M2.le.123.and.M2.gt.0.and.L.eq.IsoIPR(M2)) go to 99
      endif
 13     CONTINUE 
 12     CONTINUE        ! Close loop over spiral 
 11     CONTINUE        ! combinations      
 10     CONTINUE
 9      CONTINUE
 8      CONTINUE
 7      CONTINUE
 6      CONTINUE
 5      CONTINUE
 4      CONTINUE
 3      CONTINUE
 2      CONTINUE
 1      CONTINUE
 99   WRITE (Iout,612)
      if(iham.ne.0) then
       WRITE (Iout,609) hamlow,islow,hamhigh,ishigh
       if(IPR.eq.0.and.IPRdect.eq.1) then
        WRITE (Iout,610) IPRhamlow,islowIPR,IPRhamhigh,ishighIPR
       endif
      endif
      WRITE (Iout,611) IFus5Glow,IFusL,IFus5Ghigh,IFusH,
     1  sigmahlow,ISigmaL,sigmahhigh,ISigmaH
      WRITE (Iout,606)
 600  FORMAT(/1X,'Subroutine Spiral from Fowler and Manolopoulos',
     1 ' (An Atlas of Fullerenes, Dover Publ., New York, 2006)',
     2 /1X,'(Symmetries are given for undistorted fullerenes)',
     3 /1X,'Isomer List Start ',/,I5,2I2)
 601  FORMAT(1X,'General fullerene isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 602  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 603  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 604  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 606  FORMAT(/1X,170('-'),/1X,'End of subroutine Spiral')
 607  FORMAT(1X,I8,2X,A3,1X,12I4,2X,'(',5(I2,','),I2,')  ',I2,
     1 2X,'(',6(I2,','),I3,')  ',F8.5,2X,I2,1X,I2,1X,F8.5,
     1 1X,A6,2X,3(I3,' x',I3,:,','))
 608  FORMAT(1X,I8,2X,A3,1X,12I4,2X,'(',5(I2,','),I2,')  ',I2,
     1 2X,'(',6(I2,','),I3,')  ',F8.5,2X,I2,1X,I2,1X,F8.5,
     1 1X,A6,2X,I9,2X,3(I3,' x',I3,:,','))
 609  FORMAT(1X,'Lowest number of Hamiltonian cycles     ',I10,
     1 ' for isomer ',I10,
     1      /1X,'Highest number of Hamiltonian cycles    ',I10,
     1 ' for isomer ',I10)
 610  FORMAT(1X,'Lowest number of IPR Hamiltonian cycles ',I10,
     1 ' for isomer ',I10,
     1      /1X,'Highest number of IPR Hamiltonian cycles',I10,
     1 ' for isomer ',I10)
 611  FORMAT(1X,'Lowest  Np= ',I3,' for isomer ',I10,
     1      /1X,'Highest Np= ',I3,' for isomer ',I10,
     1      /1X,'Lowest  Sigmah= ',F8.5,' for isomer ',I10,
     1      /1X,'Highest Sigmah= ',F8.5,' for isomer ',I10)
 612  FORMAT(1X,'Isomer List Complete')
 701  FORMAT(1X,'General fullerene isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 702  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 703  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 704  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'N  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
      Return
      END

      SUBROUTINE spwindup(M,MP,D,S,JP,IER)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX)
      DIMENSION JP(12),IR(12),JR(12),JS(12)
C This routine tries to find the spiral from a preset of
C three connected rings stored in S(1), S(2) and S(3) using
C information on ring fusions from the dual matrix

      IER=0
      iring=0
      mloop=3
C   Big loop from ring 4 to the end
      Do 10 i=4,M
       if(i.gt.5) mloop=5
       ifound=0
       IP=S(i-1)
       if(IP.eq.0) go to 10
C    Find new ring connection to ring IP, but keep previous ones out
C    and store in IR
         IC=0
         Do j=1,M
           if(D(j,IP).eq.1) then
           iflag=0
            do k=1,i-2
            if(j.eq.s(k)) iflag=1
            enddo
            if(iflag.eq.0) then
              IC=IC+1
              IR(IC)=J
            endif
           endif
         enddo
       if(IC.gt.6.or.IC.eq.0) Go to 10
C    Find the ring adjacent to IP and to a previous ring
C    in the spiral and store
           iring=0
          do j=1,IC
           do k=1,i-mloop
            IX=S(k)
            if(D(IR(j),IX).eq.1) then
            iring=iring+1
            JR(iring)=IR(j)
            JS(iring)=k
            endif
           enddo
          enddo
C    Check if spiral is at dead end
      if(iring.eq.0) then
        IER=1
      Return
      endif
C    Now determine neighboring ring
         k=100000
         do j=1,iring
          if(JS(j).lt.k) then
          s(i)=JR(j)
          k=JS(j)
          endif
         enddo
C    Now check if it is a pentagon
            if(s(i).le.12) then
            MP=MP+1
            JP(MP)=i
            endif
  10  Continue
C    Spiral found if loop went through
C     Print*,M,IER,'/',(S(I),I=1,M)
C     Print*,iring,'/',JP
      Return
      END

      SUBROUTINE SpiralSearch(NSP,N,Iout,IRG55,IRG66,
     1 IRG56,NrA,NrB,NrC,NrD,NrE,NrF,JP,GROUP)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrA(EMAX),NrB(EMAX),NrC(EMAX),NrD(EMAX)
      DIMENSION NrE(EMAX),NrF(EMAX),NMR(6),JP(12)
      DIMENSION Spiral(12,NMAX)
      DIMENSION SpiralT(12,NSP),SpiralF(MMAX,NSP)
      DIMENSION D(MMAX,MMAX),S(MMAX)
      CHARACTER*3 GROUP

C     This subroutine has been modified from the original one of Fowler and 
C     Manolopoulos "An Atlas of Fullerenes" (Dover Publ., New York, 2006).           
C     It is used if dual matrix is already known. See subroutine Spiral for details.
C     N is the nuclearity of the fullerene.
      M=N/2+2
      ispiral=0
      WRITE (Iout,600)
         IF(N.lt.100) WRITE(Iout,601) N,M
         IF(N.ge.100.and.N.lt.1000) WRITE(Iout,602) N,M
         IF(N.ge.1000) WRITE(Iout,632) N,M
      do I=1,MMAX
       S(I)=0
      do J=1,MMAX
       D(I,J)=0
      enddo
      enddo
      do I=1,12
      do J=1,NMAX
       Spiral(I,J)=0
      enddo
      do J=1,NSP
       SpiralT(I,J)=0
      enddo
      enddo
C     Set up dual matrix
      do I=1,IRG55
       I1=NrA(I)
       I2=NrB(I)
       D(I1,I2)=1
       D(I2,I1)=1
      enddo
      do I=1,IRG56
       I1=NrC(I)
       I2=NrD(I)
       D(I1,I2)=1
       D(I2,I1)=1
      enddo
      do I=1,IRG66
       I1=NrE(I)
       I2=NrF(I)
       D(I1,I2)=1
       D(I2,I1)=1
      enddo

C     Search for all spirals 
C     Set up first three rings then wind
C     Start ring spiral algorithm. Quit after first successful spiral
C     JP contains the pentagon indices, S the ring numbers
      nspiral=0
      
C     Loop over all (5,5) fusions
      if(IRG55.eq.0) then
       write(Iout,610)
      else
       write(Iout,611) 2*IRG55
       do I=1,2*IRG55
        do j=4,M
         s(j)=0
        enddo
       if(I.le.IRG55) then
        I1=NrA(I)
        I2=NrB(I)
       else
        IR=I-IRG55
        I1=NrB(IR)
        I2=NrA(IR)
       endif
        S(1)=I1
        S(2)=I2
        JP(1)=1
        JP(2)=2
       do J=1,M
        if(D(I1,J).eq.1.and.D(I2,J).eq.1) then
         S(3)=J
         MP=2
         if(J.le.12) then
          JP(3)=3
          MP=3
         endif
          CALL spwindup(M,MP,D,S,JP,IER)
         do K=1,12
          if(JP(K).eq.0.or.JP(K).gt.M) IER=1
         enddo
         if(IER.eq.0) then
          nspiral=nspiral+1
          If(nspiral.gt.NSP) then
           Write(Iout,626) nspiral,nsp
           nspiral=nspiral-1
           Go to 199
          endif
          do k=1,12
           SpiralT(k,nspiral)=JP(k)
          enddo 
          do k=1,M
           SpiralF(k,nspiral)=S(k)
          enddo 
         endif
        endif
       enddo 
       enddo 
      endif

C     Loop over all (5,6) fusions
      if(IRG56.le.0) then
      write(Iout,615)
      else
      write(Iout,612) 2*IRG56
      do I=1,2*IRG56
      do j=4,M
      s(j)=0
      enddo
      if(I.le.IRG56) then
      I1=NrC(I)
      I2=NrD(I)
      else
      IR=I-IRG56
      I1=NrD(IR)
      I2=NrC(IR)
      endif
      S(1)=I1
      S(2)=I2
      if(I1.le.12) then
      JP(1)=1
      else
      JP(2)=1
      endif
      do J=1,M
      if(D(I1,J).eq.1.and.D(I2,J).eq.1) then
      S(3)=J
      MP=1
      if(J.le.12) then
      JP(2)=3
      MP=2
      endif
      CALL spwindup(M,MP,D,S,JP,IER)
         do K=1,12
          if(JP(K).eq.0.or.JP(K).gt.M) IER=1
         enddo
      if(IER.eq.0) then
      nspiral=nspiral+1
         If(nspiral.gt.NSP) then
         Write(Iout,627) nspiral,nsp
         nspiral=nspiral-1
         Go to 199
         endif
      do k=1,12
      SpiralT(k,nspiral)=JP(k)
      enddo 
      do k=1,M
      SpiralF(k,nspiral)=S(k)
      enddo 
      endif
      endif
      enddo 
      enddo 
      endif
      
C     Loop over all (6,6) fusions
      if(IRG66.eq.0) then
      write(Iout,616)
      else
      write(Iout,613) 2*IRG66
      do I=1,2*IRG66
       do j=4,M
        s(j)=0
       enddo
       if(I.le.IRG66) then
        I1=NrE(I)
        I2=NrF(I)
       else
        IR=I-IRG66
        I1=NrF(IR)
        I2=NrE(IR)
       endif
       S(1)=I1
       S(2)=I2
       do J=1,M
        if(D(I1,J).eq.1.and.D(I2,J).eq.1) then
          S(3)=J
          MP=0
         if(J.le.12) then
          JP(1)=3
          MP=1
         endif
       CALL spwindup(M,MP,D,S,JP,IER)
         do K=1,12
          if(JP(K).eq.0.or.JP(K).gt.M) IER=1
         enddo
       if(IER.eq.0) then
        nspiral=nspiral+1
         If(nspiral.gt.NSP) then
         Write(Iout,628) nspiral,nsp
         nspiral=nspiral-1
         Go to 199
         endif
       do k=1,12
        SpiralT(k,nspiral)=JP(k)
       enddo
       do k=1,M
        SpiralF(k,nspiral)=S(k)
       enddo 
       endif
      endif
      enddo 
      enddo 
      endif
     
C     Now loop over with found spiral until success with
C     Fowler algorithm
  199 IT=1
      IPR=0
      nspfound=0
      Do 13 msp=1,nspiral
       Do I=1,M
        S(I)=6
       enddo
       Do I=1,12
        S(SpiralT(I,msp))=5
        JP(I)=SpiralT(I,msp)
       enddo
       IER=0
       CALL Windup(M,IPR,IER,S,D)      !      Wind up spiral into dual 
       IF(IER.ne.0) GO TO 13                !      and check for closure 
       Do I=1,12 
        Spiral(I,1)=JP(I)
       enddo
       CALL Unwind(M,IER,IT,ispiral,
     1  Spiral,S,D,NMR,Group)             ! Unwind dual into spirals 
       K=0
       DO J=1,6
         IF(NMR(J).EQ.0) GO TO 16
         K=J
       enddo
 16    nspfound=nspfound+1
       if(nspfound.eq.1) write(Iout,614) nspiral
       If(K.le.0) then
        WRITE(Iout,603) GROUP,(JP(I),I=1,12)
       else
        WRITE(Iout,605) GROUP,(JP(I),I=1,12),(NMR(J),J=1,K)
       endif
       WRITE(Iout,604) 
       if(M.lt.1000) WRITE(Iout,618) (SpiralF(I,msp),I=1,M)
       if(M.ge.1000) WRITE(Iout,629) (SpiralF(I,msp),I=1,M)
       if(ispiral.ge.2) then
        if(ispiral.eq.2) then
         WRITE(Iout,623)
         Do II=1,12
          JP(II)=spiral(II,2)
         enddo 
        else
         WRITE(Iout,619) ispiral-1
        endif
       Do JJ=2,ispiral 
        WRITE(Iout,607) (spiral(II,JJ),II=1,12)
       enddo
       else
        WRITE(Iout,608)
       endif
       if(ispiral.gt.2) then
        CALL CanSpiral(ispiral,spiral,JP)
        WRITE(Iout,623)
        WRITE(Iout,621) (JP(I),I=1,12)
       endif
       Do I=1,M
        S(I)=6
       enddo
       Do I=1,12
        S(JP(I))=5
       enddo
       WRITE(Iout,624)
       WRITE(Iout,625) (S(I),I=1,M)
       go to 99
 13   CONTINUE 
 99   if(IER.eq.0) then
       if(ispiral.ge.2) then
C     Print ring numbers
        WRITE(Iout,620) nspiral 
       Do msp=1,nspiral
         jpc=0
       Do ipent=1,12
         jpc=jpc+iabs(JP(ipent)-SpiralT(ipent,msp))
       enddo
        if(jpc.eq.0) then
        WRITE(Iout,618) (SpiralF(I,msp),I=1,M)
        return
        endif
       enddo
       endif
      else 
      WRITE(Iout,617)
      endif
 600  FORMAT(/1X,'Modified spiral algorithm Fowler and Manolopoulos',
     1 ' (An Atlas of Fullerenes, Dover Publ., New York, 2006)')
 601  FORMAT(1X,'Spiral for fullerene isomers of C',I2,':',
     1 ' (',I2,' faces)')
 602  FORMAT(1X,'Spiral for fullerene isomers of C',I3,':',
     1 ' (',I3,' faces)')
 603  FORMAT(1X,A3,9X,12I4)
 604  FORMAT(1X,90('-'),/1X,'Corresponding ring numbers:') 
 605  FORMAT(1X,A3,9X,12I4,2X,3(I3,' x',I3,:,','))
 607  Format(12(1X,I4))
 608  Format(1X,'Input spiral is canonical')
 610  Format(1X,'This is an IPR fullerene, no (5,5) fusions to ',
     1 'loop over')
 611  Format(1X,'Loop over (5,5) fusions, ',I5,' max in total')
 612  Format(1X,'Loop over (5,6) fusions, ',I5,' max in total')
 613  Format(1X,'Loop over (6,6) fusions, ',I5,' max in total')
 614  Format(1X,I6,' potential spirals found',/1X,
     1 'Point group   Ring spiral pentagon positions',
     2 19X,'NMR pattern (for fullerene in ideal symmetry)',/1X,90('-')) 
 615  Format(1X,'This is C20, no (5,6) fusions to loop over')
 616  Format(1X,'No (6,6) fusions to loop over')
 617  Format(1X,'Failed to find ring spiral: Fullerene most likely a ',
     1 'non-spiral one')
 618  Format(20(1X,32(I4,'-'),/))
 619  Format(1X,'Spiral list of pentagon positions with ',
     1 'higher priority: (',I4,' spirals found)') 
 620  Format(1X,'Search ',I6,' spirals to produce canonical'
     1 ' list of atoms:')
 621  Format(12(1X,I4))
 623  Format(1X,'Canonical spiral list of pentagon positions:')
 624  Format(1X,'Canonical spiral list of hexagons and pentagons:')
 625  Format(1X,100I1)
 626  Format(1X,'**** Severe Warning: Number of detected spirals is ',
     1 I5,' greater than dimension in field ',I5,' detected in (5,5)',
     1 'list',/1X,'Routine will stop here and tries to work with ',
     1 'existing spirals (otherwise increase NSpScale parameter ',
     1 'in main program')
 627  Format(1X,'**** Severe Warning: Number of detected spirals is ',
     1 I5,' greater than dimension in field ',I5,' detected in (5,6)',
     1 'list',/1X,'Routine will stop here and tries to work with ',
     1 'existing spirals (otherwise increase NSpScale parameter ',
     1 'in main program')
 628  Format(1X,'**** Severe Warning: Number of detected spirals is ',
     1 I5,' greater than dimension in field ',I5,' detected in (6,6)',
     1 'list',/1X,'Routine will stop here and tries to work with ',
     1 'existing spirals (otherwise increase NSpScale parameter ',
     1 'in main program')
 629  Format(20(1X,32(I4,'-'),/))
 632  FORMAT(1X,'Spiral for fullerene isomers of C',I4,':',
     1 ' (',I4,' faces)')
      Return
      END
      
      SUBROUTINE CanSpiral(MS0,S,P)
      use config
      IMPLICIT Integer (A-Z)
      DIMENSION S(12,NMAX),P(12),PI(12),SM(12,NMAX)
C     Find canonical spiral by sorting
      IS=1
      MS=MS0
         Do I=1,12
         PI(I)=0
         enddo
C     Find lowest value
  1   Smax=100000
          Do I=1,MS
          IF(S(IS,I).le.Smax) Smax=S(IS,I)
          enddo
      IZ=0
          Do I=1,MS
           if(S(IS,I).eq.Smax) then
           IZ=IZ+1
            do j=1,12
            SM(j,IZ)=S(j,I)
            enddo
           endif
           enddo
      MS=IZ
       if(MS.eq.1) then
          Do I=1,12
          P(I)=SM(I,1)
          enddo
          return
       else
         Do I=1,MS
         Do J=1,12
          S(J,I)=SM(J,I)
         enddo
         enddo
       IS=IS+1
       go to 1
      endif
      return
      END

      SUBROUTINE DualAnalyze(nmax,mmax,N,M,D,IRhag5,IRhag6,
     1 IFus5G,IDA,nelec,ndeg,sigmah,A,gap)
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer D(MMAX,MMAX),Ddiag(MMAX),NR5(12),IDA(NMAX,NMAX)
      Integer IRhag5(0:5),IRhag6(0:6),IDG(NMAX)
      Dimension A(NMAX,NMAX),evec(NMAX),df(NMAX)
C     Analyze dual matrix and get pentagon and hexagon indices
C     and determine if molecule is open shell
      Data Tol/1.d-5/

      NV5=0
      NV6=0
      IER=0
      Do I=1,M
       Ddiag(I)=0
       Do J=1,M
        Ddiag(I)=Ddiag(I)+D(J,I)
       enddo
       Itest=Ddiag(I)
       If(Itest.eq.5) then
        NV5=NV5+1
        NR5(NV5)=I
       endif
       If(Itest.eq.6) NV6=NV6+1
      enddo
C     NV=NV5+NV6
C     If(NV.ne.M) Write(Iout,1000)
C     Get Rhagavachari-Fowler-Manolopoulos neighboring pentagon and hexagon indices
C     First pentagon indices
      Do I=0,5
       IRhag5(I)=0
      enddo
      do I=1,12
       IRcount=0
       IRing5=NR5(I)
       do J=1,12
        JRing5=NR5(J)
        If(D(JRing5,IRing5).eq.1) then
         IRcount=IRcount+1
        endif
       enddo
       IRhag5(IRcount)=IRhag5(IRcount)+1
      enddo
      Ifus5=0
      Do I=1,5
       IFus5=IFus5+I*IRhag5(I)
      enddo
      IFus5G=IFus5/2

C     Now hexagon indices
      Do I=0,6
       IRhag6(I)=0
      enddo
      do 10 I=1,M
       IRcount=0
       IR5=Ddiag(I)
       if(IR5.eq.5) go to 10
       do J=1,M
        JR5=Ddiag(J)
        If(JR5.ne.5.and.D(I,J).eq.1) then
         IRcount=IRcount+1
        endif
       enddo
       IRhag6(IRcount)=IRhag6(IRcount)+1
   10 Continue
C     Strain Parameter
      khk=0
      k2hk=0
      Do I=3,6
      ihk=ihk+IRhag6(I)
      IIR=I*IRhag6(I)
      khk=khk+IIR
      k2hk=k2hk+I*IIR
      enddo
      if(ihk.eq.0) go to 112
      aihk=dfloat(ihk)
      akhk2=(dfloat(khk)/aihk)**2
      ak2hk=dfloat(k2hk)/aihk
      sigmah=dsqrt(dabs(ak2hk-akhk2))

C     Now produce adjacency matrix
 112  CALL DUAL(D,MMAX,IDA,N,IER)
      Do I=1,N
       df(I)=0.d0
      Do J=I,N
        A(I,J)=0.d0
        if(IDA(I,J).eq.1) A(I,J)=1.d0
        A(J,I)=A(I,J)
      enddo
      enddo
C Diagonalize without producing eigenvectors
      do i=1,nmax
         evec(i) = 0
      enddo
      call tred2l(A,N,NMAX,evec,df)
      call tqlil(evec,df,N,NMAX)
C Sort eigenvalues
      Do I=1,N
       e0=evec(I)
       jmax=I
        Do J=I+1,N
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

C Now sort degeneracies
      df(1)=evec(1)
      ieigv=1
      ideg=1
      IDG(1)=ideg
      Do I=2,N
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
       
C Produce number of electrons in HOMO, degeneracy and gap
      Noc=N/2
      Norb=0
      Do I=1,ieigv
      Iorb=I
      Norb=Norb+IDG(I)
      if(Norb.eq.Noc) then 
      gap=df(Iorb)-df(Iorb+1)
      ndeg=IDG(Iorb)
      nelec=ndeg*2
      go to 111
      endif
      if(Norb.gt.Noc) then 
      gap=df(Iorb)-df(Iorb+1)
      ndeg=IDG(Iorb)
      nelec=(ndeg-Norb+Noc)*2
      go to 111
      endif
      enddo
  111 Continue
      
      Return
      END

C     Inputs: M, IPR, S(M)
C     Output: D(M,M)
      SUBROUTINE Windup2(M,IPR,S,ier,D)
      use config
      integer bigD(Mmax,Mmax), D(M,M), i, j, ier

      call Windup(M, IPR, S, ier, bigD)
      do i=1,M
         do j=1,M
            D(i,j) = bigD(i,j)
         enddo
      enddo   
      END

      SUBROUTINE Windup(M,IPR,IER,S,D)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION S(MMAX),D(MMAX,MMAX)
      DIMENSION R(MMAX),C(MMAX,6)
C       This subroutine attempts to wind up an input spiral S into
C       a fullerene dual (face) adjacency matrix D. It returns with      
C       IER = P if the spiral shorts or is discovered to be open-ended      
C       after P pentagons have been added. Otherwise IER = 0 on return.
      J=1
      C(1,1)=2
      C(2,1)=1
      R(1)=2
      R(2)=2
      E=1
      P=12-S(1)-S(2)
      DO 5 K=3,M-1
         P=P+6-S(K)
         R(K)=1
         I=K-1
 1       IF(IPR.EQ.1.AND.S(I).EQ.5.AND.S(K).EQ.5) GO TO 10      
         IF(R(K).GE.S(K)) GO TO 10      
         C(I,R(I))=K    ! Connect face K to the last open face I
         C(K,R(K))=I    ! in the preceding spiral      
         R(I)=R(I)+1
         R(K)=R(K)+1
         IF(R(I).GT.S(I)) THEN
            L=I-1              ! If this closes face I update I and go
            DO I=L,J+1,-1    ! back to connect face K to the new I      
             IF(R(I).LE.S(I)) GO TO 1
            enddo
            GO TO 10
         ENDIF
 3       IF(IPR.EQ.1.AND.S(J).EQ.5.AND.S(K).EQ.5) GO TO 10
         IF(R(K).GE.S(K)) GO TO 10
         C(J,R(J))=K    ! Connect face K to the first open face J
         C(K,R(K))=J    ! in the preceding spiral 
         R(J)=R(J)+1
         R(K)=R(K)+1
         IF (R(J).GT.S(J)) THEN
            L=J+1              ! If this closes face J update J and go
            DO J=L,I-1,+1      ! back to connect face K to the new J 
             IF(R(J).LE.S(J)) GO TO 3
            enddo
            GO TO 10
         ENDIF
         H=K-P
         E=E+R(K)-1     ! Use Euler's theorem to streamline the 
         V=3+2*P+3*H-E  ! search. F is a lower bound on the # of 
         F=(V+1)/2+1    ! additional faces required for closure 
         IF(F.GT.M-K) GO TO 10
 5    CONTINUE
      P=12
      R(M)=1
      DO 6 K=J,M-1
         IF (R(K).LT.S(K)) GO TO 10
         IF (R(K).GT.S(K)) GO TO 6
         IF (R(M).GT.S(M)) GO TO 10
         IF (IPR.EQ.1.AND.S(K).EQ.5.AND.S(M).EQ.5) GO TO 10
         C(K,R(K))=M    ! Connect face M to all remaining
         C(M,R(M))=K    ! open faces (including face M-1)
         R(K)=R(K)+1
         R(M)=R(M)+1
 6    CONTINUE
      IF(R(M).LE.S(M)) GO TO 10
      P=0               ! Successful spiral
      DO J=1,M        ! Form dual adjacency matrix in D
         DO I=1,M
            D(I,J)=0
         enddo
         DO K=1,S(J)
            I=C(J,K)
            D(I,J)=1
         enddo
      enddo
 10   IER=P
      RETURN
      END

      SUBROUTINE Unwind(M,IER,IT,ispiral,
     1 Spiral,S,D,NMR,GROUP)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),NMR(6)
      DIMENSION P(MMAX),R(MMAX)
      DIMENSION V(3,NMAX),E(2,EMAX)
      DIMENSION VP(NMAX,MMAX),EP(EMAX,MMAX)
      DIMENSION FP(MMAX,MMAX),Spiral(12,NMAX)
      DIMENSION MV(12),ME(12),MF(12),MS(12)
      CHARACTER*3 GROUP
C       This subroutine unwinds a fullerene dual adjacency matrix D 
C       into each of its constituent spirals and checks that none has      
C       a lexicographically smaller code than the input spiral S. The 
C       idealized point group and NMR signature of the fullerene are      
C       also calculated if this test is passed, in which case the      
C       spiral is canonical and IER = 0 on return. Otherwise IER = 13.
C       IF IT=0 this option is supressed as noncanonical spirals are
C       analyzed as well.      
      SO=0
      ispiral=0
      DO 10 I1=1,M      ! Begin multiple loop over all      
         P(1)=I1        ! 6*N possible spiral starts 
         FLAG1=0        ! with initial faces I1,I2,I3      
         IF (S(P(1)).NE.S(1)) THEN
            IF (S(P(1)).GT.S(1)) GO TO 10
            FLAG1=1
         ENDIF
         DO 9 I2=1,M
            IF(D(I1,I2).EQ.0) GO TO 9     
            P(2)=I2
            FLAG2=FLAG1
            IF(FLAG2.EQ.0.AND.S(P(2)).NE.S(2)) THEN
                 IF(S(P(2)).GT.S(2)) GO TO 9
                 FLAG2=2
            ENDIF
            DO 8 I3=1,M
           IF(D(I1,I3).EQ.0.OR.D(I2,I3).EQ.0) GO TO 8
           IF(SO.EQ.0) THEN
              SO=1        !     Store a face permutation for
              DO K=1,M    !     each symmetry operation in FP, 
               FP(K,SO)=K !     with the identity operation 
              enddo       !     (here) in column 1      
              GO TO 8
           ENDIF
           P(3)=I3
           FLAG3=FLAG2
           IF(FLAG3.EQ.0.AND.S(P(3)).NE.S(3)) THEN
              IF(S(P(3)).GT.S(3)) GO TO 8
              FLAG3=3
           ENDIF
           DO J=1,M
              R(J)=0
           enddo
           R(P(1))=2
           R(P(2))=2
           R(P(3))=2
           I=1
         DO 6 J=4,M
 3            IF(R(P(I)).EQ.S(P(I))) THEN
                 I=I+1
                 IF(I.EQ.J-1) GO TO 8
                 GO TO 3
              ENDIF
                 IFX=P(I)    !  These are the first (IF)      
                 IL=P(J-1)   !  and last (IL) open faces      
                 DO 5 IJ=1,M !  in the preceding spiral 
                IF (D(IJ,IFX).EQ.0.OR.D(IJ,IL).EQ.0) GO TO 5
                IF (R(IJ).GT.0) GO TO 5
                P(J)=IJ
                IF (FLAG3.EQ.0.AND.S(P(J)).NE.S(J)) THEN
                   IF (S(P(J)).GT.S(J)) GO TO 8
                   FLAG3=J    ! This spiral has a smaller      
                ENDIF         ! code than S, but it may not      
                DO 4 K=1,J-1  ! close properly. Flag it      
                   IF (D(P(J),P(K)).EQ.1) THEN
                  R(P(J))=R(P(J))+1
                  R(P(K))=R(P(K))+1
                   ENDIF
 4            CONTINUE
                GO TO 6
 5             CONTINUE
                 GO TO 8
 6          CONTINUE
              IF (FLAG3.EQ.0) THEN
                 SO=SO+1      ! Arrive here once for each 
                 DO 7 K=1,M   ! spiral with the same code as 
                FP(K,SO)=P(K) ! S, which is once for each 
 7             CONTINUE       ! symmetry operation SO      
              ELSE
                 If(IT.eq.0) then
                 IER=13       ! The flagged spiral has closed, 
                 RETURN       ! so call it a day  
                 endif    
              ENDIF
C            Store all non-identical spirals
             If(IT.ne.0) then
             ISP=0
             ispiral=ispiral+1
             Do II=1,M
             IF(S(P(II)).eq.5) then
             ISP=ISP+1
             Spiral(ISP,ispiral)=II
             endif
             enddo
             Do II=1,ispiral-1
             nv=0
             Do JJ=1,12
             nv=nv+iabs(Spiral(JJ,II)-Spiral(JJ,ispiral))
             if(Spiral(JJ,ispiral).eq.0) go to 86
             if(JJ.ne.1) then
             ndif=Spiral(JJ-1,ispiral)-Spiral(JJ,ispiral)
             if(ndif.eq.0) go to 86
             endif
             enddo
             if(nv.eq.0) go to 86
             enddo
             go to 87
 86          ispiral=ispiral-1
 87          continue
             endif
 8           CONTINUE
 9            CONTINUE
 10            CONTINUE
         IER=0      !         Spiral S is canonical, and 
         ORDER=SO   !         SO is the point group order.
         N=0        !         Now calculate GROUP and NMR: 
         L=0
         DO 13 K=2,M
            DO 12 J=1,K-1
           IF (D(J,K).EQ.0) GO TO 12
           DO 11 I=1,J-1
              IF (D(I,J).EQ.0.OR.D(I,K).EQ.0) GO TO 11
              N=N+1
              V(1,N)=I    ! Associate the three mutually 
              V(2,N)=J    ! adjacent faces I,J,K 
              V(3,N)=K    ! with vertex N 
 11      CONTINUE
           L=L+1
           E(1,L)=J    ! And the two mutually adjacent 
           E(2,L)=K    ! faces J,K with edge L 
 12     CONTINUE
 13    CONTINUE
       DO 18 SO=1,ORDER
         DO 15 J=1,N
           J1=FP(V(1,J),SO)
           J2=FP(V(2,J),SO)
           J3=FP(V(3,J),SO)
           I1=MIN(J1,J2,J3)
           I3=MAX(J1,J2,J3)
           I2=J1+J2+J3-I1-I3
           DO 14 I=1,N
              IF (V(1,I).EQ.I1.AND.V(2,I).EQ.I2.AND.V(3,I).EQ.I3) THEN 
                 VP(J,SO)=I ! Store a vertex permutation for 
                 GO TO 15   ! each symmetry operation in VP
              ENDIF
 14         CONTINUE
 15    CONTINUE
       DO 17 J=1,L
           J1=FP(E(1,J),SO)
           J2=FP(E(2,J),SO)
           I1=MIN(J1,J2)
           I2=J1+J2-I1
           DO 16 I=1,L
              IF(E(1,I).EQ.I1.AND.E(2,I).EQ.I2) THEN
                 EP(J,SO)=I ! And similarly an edge permutation 
                 GO TO 17   ! in EP
              ENDIF
 16         CONTINUE
 17     CONTINUE
 18    CONTINUE
         DO 19 K=1,12
            MV(K)=0
            ME(K)=0
            MF(K)=0
 19      CONTINUE
         DO 21 J=1,N
            IF (VP(J,1).EQ.0) GO TO 21
            VP(J,1)=0
            K=1
            DO 20 SO=2,ORDER
              I=VP(J,SO)
              IF(VP(I,1).EQ.0) GO TO 20
              VP(I,1)=0
              K=K+1
 20         CONTINUE
            K=ORDER/K        ! Count vertex orbits with 
            MV(K)=MV(K)+1    ! site group order K in MV(K)      
 21      CONTINUE
         DO 22 J=1,N
            VP(J,1)=J
 22         CONTINUE
         DO 24 J=1,L
            IF (EP(J,1).EQ.0) GO TO 24
            EP(J,1)=0
            K=1
         DO 23 SO=2,ORDER
           I=EP(J,SO)
           IF(EP(I,1).EQ.0) GO TO 23
           EP(I,1)=0
           K=K+1
 23      CONTINUE
            K=ORDER/K        ! And edge orbits with 
            ME(K)=ME(K)+1    ! site group order K in ME(K)      
 24         CONTINUE
         DO 25 J=1,L
            EP(J,1)=J
 25         CONTINUE
         DO 27 J=1,M
            IF (FP(J,1).EQ.0) GO TO 27
            FP(J,1)=0
            K=1
         DO 26 SO=2,ORDER
           I=FP(J,SO)
           IF(FP(I,1).EQ.0) GO TO 26
           FP(I,1)=0
           K=K+1
 26      CONTINUE
            K=ORDER/K      ! And face orbits with
            MF(K)=MF(K)+1  ! site group order K in MF(K)
 27         CONTINUE
         DO 28 J=1,M
            FP(1,J)=J
 28         CONTINUE
         DO 29 K=1,12                !  And ALL special point orbits
            MS(K)=MV(K)+ME(K)+MF(K)  !  with site group order K in MS(K)
 29         CONTINUE
         DO 30 J=1,6
            NMR(J)=0
 30         CONTINUE
         J=0
         DO 31 K=6,1,-1              ! Use the vertex orbit counts
            IF (MV(K).EQ.0) GO TO 31 ! to calculate the NMR pattern
            J=J+1
            NMR(J)=MV(K)
            J=J+1
            NMR(J)=ORDER/K
 31         CONTINUE
         GROUP='???'                 !  And, finally, the full
         IF (ORDER.EQ.1) THEN        !  special point orbit counts
            GROUP = ' C1'            !  (in conjunction with the
         ELSE IF (ORDER .EQ. 2) THEN !  point group order) to assign
            IF (MS(2).EQ.0) THEN     !  the point group
           GROUP = ' Ci'
            ELSE IF (MS(2).EQ.2) THEN
           GROUP = ' C2'
            ELSE IF (MS(2).GT.2) THEN
           GROUP = ' Cs'
           ENDIF
         ELSE IF (ORDER.EQ.3) THEN
            GROUP = ' C3'
         ELSE IF (ORDER.EQ.4) THEN
            IF (MS(4).EQ.0) THEN
           IF (MS(2) .EQ. 1) THEN
              GROUP = ' S4'
           ELSE IF (MS(2).EQ.3) THEN
              GROUP = ' D2'
           ELSE IF (MS(2).GT.3) THEN
              GROUP = 'C2h'
           ENDIF
            ELSE IF (MS(4).EQ.2) THEN
           GROUP = 'C2v'
            ENDIF
         ELSE IF (ORDER.EQ.6) THEN
            IF (MS(6).EQ.0) THEN
           IF (MS(2).EQ.0) THEN
              GROUP = ' S6'
           ELSE IF (MS(2).EQ.2) THEN
              GROUP = ' D3'
           ELSE IF (MS(2).GT.2) THEN
              GROUP = 'C3h'
           ENDIF
            ELSE IF (MS(6).EQ.2) THEN
           GROUP = 'C3v'
            ENDIF
         ELSE IF (ORDER.EQ.8) THEN
            IF (MS(4).EQ.1) THEN
           GROUP = 'D2d'
            ELSE IF (MS(4).EQ.3) THEN
           GROUP = 'D2h'
            ENDIF
         ELSE IF (ORDER.EQ.10) THEN
            GROUP = ' D5'
         ELSE IF (ORDER.EQ.12) THEN
            IF (MS(6).EQ.0) THEN
           GROUP ='  T'
            ELSE IF (MS(6).EQ.1) THEN
           IF (MS(4).EQ.0) THEN
              IF (MS(2).EQ.2) THEN
                 GROUP = ' D6'
              ELSE IF (MS(2).GT.2) THEN
                 GROUP = 'D3d'
              ENDIF
           ELSE IF (MS(4).EQ.2) THEN
              GROUP = 'D3h'
           ENDIF
            ENDIF
         ELSE IF (ORDER.EQ.20) THEN
            IF (MS(4).EQ.0) THEN
           GROUP = 'D5d'
            ELSE IF (MS(4).EQ.2) THEN
           GROUP = 'D5h'
            ENDIF
         ELSE IF (ORDER.EQ.24) THEN
            IF (MS(12).EQ.0) THEN
           IF (MS(6).EQ.0) THEN
              GROUP = ' Th'
           ELSE IF (MS(6) .EQ. 2) THEN
              GROUP = ' Td'
          ENDIF
           ELSE IF (MS(12) .EQ. 1) THEN
          IF (MS(4).EQ.0) THEN
             GROUP = 'D6d'
          ELSE IF (MS(4) .EQ. 2) THEN
             GROUP = 'D6h'
          ENDIF
           ENDIF 
        ELSE IF (ORDER.EQ.60) THEN
           GROUP ='  I'
        ELSE IF (ORDER.EQ.120) THEN
           GROUP = ' Ih'
        ENDIF
        RETURN
        END

      SUBROUTINE DUAL(D,M,A,N,IER)
      use config
      IMPLICIT INTEGER (A-Z)
        DIMENSION D(M,M),A(NMAX,NMAX)
        DIMENSION V(3,NMAX)
c       Given a fullerene dual adjacency matrix D, this subroutine 
c       constructs the corresponding fullerene adjacency matrix A. 
c       IER = 0 on return if the construction is successful. 	 
        I=0
        DO 3 L = 1,M 
           DO 2 K=1,L 
              IF (D(K,L).EQ.0) GO TO 2
              DO 1 J = 1,K
                 IF (D(J,K).EQ.0.OR.D(J,L).EQ.0) GO TO 1
                  I = I+1
                 IF (I.GT.N) GO TO 1
                 V(1,I) = J ! Associate the three mutually adjacent 
                 V(2,I) = K ! dual vertices (fullerene faces) J,K,L 
                 V(3,I) = L ! with fullerene vertex I 	 
 1            CONTINUE
 2        CONTINUE
 3      CONTINUE
        IER = I-N
        IF (IER .NE. 0) RETURN ! D contains IER > 0 separating triangles 
        DO 7 J = 1,N           ! and is therefore NOT a fullerene dual 
           DO 6 I = 1,J
              K = 0
              DO 5 JJ = 1,3
                 DO 4 II = 1,3
                  IF(V(II,I).EQ.V(JJ,J)) K = K+1
 4               CONTINUE
 5             CONTINUE
              IF (K.EQ.2) THEN
                 A(I,J)=1   ! Fullerene vertices I and J are adjacent 
                 A(J,I)=1   ! if they have 2 dual vertices in common
               ELSE
                 A(I,J)=0
                 A(J,I)=0
              ENDIF
 6         CONTINUE
 7      CONTINUE
        RETURN
        END
