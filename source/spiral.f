      SUBROUTINE SpiralRestart(IPR,Iout,
     1 Isonum,IsoIPR,iham,IDA,A,chkname)
C     Restart version of Subroutine Spiral
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),IDA(NMAX,NMAX),IC3(NMAX,3)
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
      maxiter=1000000000
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
      if(IN.ne.number_vertices) then
       Write(Iout,1003) number_vertices,IN
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
       WRITE (Iout,1015) I
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
       M1=number_vertices/2-9
       M2=number_vertices/2-29
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
         IF(number_vertices.lt.100) WRITE(Iout,601) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,602) number_vertices
      else
         IPR=1 
         IF(number_vertices.lt.100) WRITE(Iout,603) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,604) number_vertices
      endif
      else
      if(IPR.EQ.0) then 
         IF(number_vertices.lt.100) WRITE(Iout,701) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,702) number_vertices
      else
         IPR=1 
         IF(number_vertices.lt.100) WRITE(Iout,703) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,704) number_vertices
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
      M=number_vertices/2+2
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
   16  CALL DualAnalyze(M,D,IRhag5,IRhag6,
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
         do ia=1,number_vertices
           ka=0
            do ja=1,number_vertices
              if(IDA(Ia,Ja).eq.1) then
               ka=ka+1
               IC3(ia,ka)=ja
               endif
            end do
         end do
        Call HamiltonCyc(maxiter,Iout,nbatch,IC3,nhamcycle)
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
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 602  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 603  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 604  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
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
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 702  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 703  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 704  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 1000 Format(1X,'Opening file: ',A20)
 1001 Format(1X,'Looped over ',I9,' lines. Characterestic card ',
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
     1 'an empty line at the end of the input, line ,'I10)
 1023 Format(1X,'Filename ',A50,' in database not found ==> ABORT')
 2000 Format(I9,2X,A3,1X,12I4,23X,I2,27X,F8.5)
 2001 Format(I9,2X,A3,1X,12I4,23X,I2,27X,F8.5,25X,I9)
 2002 Format(A18)
      END subroutine spiralrestart


      SUBROUTINE SpiralFind(IPR,ivar,In,Iout,IDA,A)
C     This subroutine comes directly from the book of Fowler and 
C     Manolopoulos "An Atlas of Fullerenes" (Dover Publ., New York, 2006),
C     and has been modified to search for ring spirals around an
C     icosahedral fullerene closest to the vertex number required.           
C     This sub-program catalogues fullerenes with a given number of      
C     vertices using the spiral algorithm and a uniqueness test 
C     based on equivalent spirals. IPR = 0  for all isomers. 
C     The resulting output is a catalogue of the isomers found containing 
C     their idealized point groups, canonical spirals, and NMR patterns.        
C     number_vertices is the nuclearity of the fullerene.
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),IDA(NMAX,NMAX)
      DIMENSION NMR(6),IRhag5(0:5),IRhag6(0:6)
      DIMENSION Spiralx(12,NMAX),IVL(12),IVH(12),SS(12)
      Real*8 sigmah,sigmahlow,sigmahhigh
      Real*8 A(NMAX,NMAX),gap,dex
      CHARACTER*3 GROUP
      CHARACTER*6 Occup

      I=0
      J=0
      ivarlimit=4
      if(ivar.gt.ivarlimit) then
       Write(Iout,805)
       Return
      endif
      if(number_vertices.gt.NMAX) Return     ! Increase NMAX 
      if(2*(number_vertices/2).ne.number_vertices) Return  ! number_vertices must be even 

C   Search for icosahedral fullerene Cm closest to Cm with m<n
      dex=dfloat(number_vertices)
      icon=0
      loopmax=int(dsqrt(dex/2.d1))
      loopmin=int(dsqrt(dex/6.d1))
      if(loopmin.lt.1) loopmin=1
       Do I1=loopmin,loopmax
       Do J1=0,I1
        nico=20*(I1*I1+J1*J1+I1*J1)
        if(nico.le.number_vertices.and.nico.gt.icon) then
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
       if(icon.eq.number_vertices) then
        Write(Iout,801) icon,I,J
       else
        Write(Iout,802) icon,number_vertices,I,J
       endif

      M=number_vertices/2+2
      Read(IN,*,Err=400,end=400) (SS(I),I=1,12)
      Write(Iout,807)
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

  300 Write(Iout,804) (SS(I),I=1,12),ivar,ivar
C  Set loop limits
      IVL(1)=1
      IVH(1)=2
      do I=2,12
       IVL(I)=SS(I)-ivar
       IVH(I)=SS(I)+ivar
      enddo
      if(SS(12).gt.M) then
       Write(Iout,*) 'Last RSPI ',SS(12),'>',M,' ==> RETURN'
       return
      endif
      

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

      if(IVH(12).gt.M) IVH(12)=M
      DO 1  J1=IVL(1),IVH(1)
       IVL2=IVL(2)
       IVH2=IVH(2)
       if(J1.ge.IVL2) IVL2=J1+1
       if(IVL2.gt.IVH2) IVH2=IVL2
       if(IVH2.gt.M) then
        Write(Iout,*) 'Loop1 ',IVH2,'>',M,' ==> RETURN'
        Return
       endif 
      DO 2  J2=IVL2,IVH2         !   combinations
       IVL3=IVL(3)
       IVH3=IVH(3)
       if(J2.ge.IVL3) IVL3=J2+1
       if(IVL3.gt.IVH3) IVH3=IVL3
       if(IVH3.gt.M)  then
        Write(Iout,*) 'Loop2 ',IVH3,'>',M,' ==> RETURN'
        Return
       endif
      DO 3  J3=IVL3,IVH3
       IVL4=IVL(4)
       IVH4=IVH(4)
       if(J3.ge.IVL4) IVL4=J3+1
       if(IVL4.gt.IVH4) IVH4=IVL4
       if(IVH4.gt.M)  then
        Write(Iout,*) 'Loop3 ',IVH4,'>',M,' ==> RETURN'
        Return
       endif
      DO 4  J4=IVL4,IVH4
       IVL5=IVL(5)
       IVH5=IVH(5)
       if(J4.ge.IVL5) IVL5=J4+1
       if(IVL5.gt.IVH5) IVH5=IVL5
       if(IVH5.gt.M)  then
        Write(Iout,*) 'Loop4 ',IVH5,'>',M,' ==> RETURN'
        Return
       endif
      DO 5  J5=IVL5,IVH5
       IVL6=IVL(6)
       IVH6=IVH(6)
       if(J5.ge.IVL6) IVL6=J5+1
       if(IVL6.gt.IVH6) IVH6=IVL6
       if(IVH6.gt.M)  then
        Write(Iout,*) 'Loop5 ',IVH6,'>',M,' ==> RETURN'
        Return
       endif
      DO 6  J6=IVL6,IVH6
       IVL7=IVL(7)
       IVH7=IVH(7)
       if(J6.ge.IVL7) IVL7=J6+1
       if(IVL7.gt.IVH7) IVH7=IVL7
       if(IVH7.gt.M)  then
        Write(Iout,*) 'Loop6 ',IVH7,'>',M,' ==> RETURN'
        Return
       endif
      DO 7  J7=IVL7,IVH7
       IVL8=IVL(8)
       IVH8=IVH(8)
       if(J7.ge.IVL8) IVL8=J7+1
       if(IVL8.gt.IVH8) IVH8=IVL8
       if(IVH8.gt.M)  then
        Write(Iout,*) 'Loop7 ',IVH8,'>',M,' ==> RETURN'
        Return
       endif
      DO 8  J8=IVL8,IVH8
       IVL9=IVL(9)
       IVH9=IVH(9)
       if(J8.ge.IVL9) IVL9=J8+1
       if(IVL9.gt.IVH9) IVH9=IVL9
       if(IVH9.gt.M)  then
        Write(Iout,*) 'Loop8 ',IVH9,'>',M,' ==> RETURN'
        Return
       endif
      DO 9  J9=IVL9,IVH9
       IVL10=IVL(10)
       IVH10=IVH(10)
       if(J9.ge.IVL10) IVL10=J9+1
       if(IVL10.gt.IVH10) IVH10=IVL10
       if(IVH10.gt.M)  then
        Write(Iout,*) 'Loop9 ',IVH10,'>',M,' ==> RETURN'
        Return
       endif
      DO 10 J10=IVL10,IVH10
       IVL11=IVL(11)
       IVH11=IVH(11)
       if(J10.ge.IVL11) IVL11=J10+1
       if(IVL11.gt.IVH11) IVH11=IVL11
       if(IVH11.gt.M)  then
        Write(Iout,*) 'Loop10 ',IVH11,'>',M,' ==> RETURN'
        Return
       endif
      DO 11 J11=IVL11,IVH11
       IVL12=IVL(12)
       IVH12=IVH(12)
       if(J11.ge.IVL12) IVL12=J11+1
       if(IVL12.gt.IVH12) IVH12=IVL12
       if(IVH12.gt.M)  then
        Write(Iout,*) 'Loop11 ',IVH12,'>',M,' ==> RETURN'
        Return
       endif
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
      IF(IER.EQ.12) GO TO 12          !      and check for closure 
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
   16  CALL DualAnalyze(M,D,IRhag5,IRhag6,
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
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
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
     1 'with number_vertices=',I6,' ,(k,l)=(',I2,',',I2,')')
 802  FORMAT(1X,'Closest vertex number is ',I6,' for icoshedral ',
     1 'fullerene with number_vertices=',I6,'; (k,l)=(',I2,',',I2,')')
 803  Format(1X,'General Goldberg-Coxeter transformation of C20 -> Cn',
     1 ' with n=: ',I5,/1X,'Construction of icosahedral fullerenes ',
     1 'using the spiral code representation of Fowler and Rogers',/1X,
     1 ' Lit: P. W. Fowler, K. M. Rogers, J. Chem. Inf. Comput. Sci.',
     1 ' 41, 108-111 (2001)',/1X,'Exponents of spiral series ',
     1 '5(6)^A (5(6)^B)^4 5(6)^C (5(6)^D)^4 5(6)^E 5:',
     1 '  A=',I5,', B=',I5,', C=',I5,', D=',I5,', E=',I5)
 804  Format(1X,'Ring spiral pentagon indices RSPI: ',12I6,
     1 /1X,'Variation of each RSPI from RSPI-',I1,' to ','RSPI+',I1)
 805  Format(1X,'Variation too large, too many loops ==> RETURN')
 807  Format(1X,'Pentagon indices taken from input')

C     RETURN
      END subroutine spiralfind

      SUBROUTINE Spiral(IPR,Iout,Isonum,IsoIPR,iham,IDA,A)
C     This subroutine comes directly from the book of Fowler and 
C     Manolopoulos "An Atlas of Fullerenes" (Dover Publ., New York, 2006).           
C     This sub-program catalogues fullerenes with a given number of      
C     vertices using the spiral algorithm and a uniqueness test 
C     based on equivalent spirals. The required input is IPR, 
C     where IPR = 0 for general and 1 for isolated-pentagon isomers. 
C     The resulting output is a catalogue of the isomers found containing 
C     their idealized point groups, canonical spirals, and NMR patterns.        
C     number_vertices is the nuclearity of the fullerene.
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(MMAX,MMAX),S(MMAX),IDA(NMAX,NMAX),IC3(NMAX,3)
      Real*8 A(NMAX,NMAX),gap
      Integer Isonum(119),IsoIPR(123)
      DIMENSION NMR(6),IRhag5(0:5),IRhag6(0:6)
      DIMENSION Spiralx(12,NMAX)
      CHARACTER*3 GROUP
      CHARACTER*6 Occup
      Real*8 sigmah,sigmahlow,sigmahhigh

      if(number_vertices.gt.NMAX) Return            ! Increase NMAX 
      if(2*(number_vertices/2).ne.number_vertices) Return         ! number_vertices must be even 
      M1=number_vertices/2-9
      M2=number_vertices/2-29
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
      WRITE (Iout,600) number_vertices,IPR,iham
      if(iham.eq.0) then
        if(IPR.EQ.0) then 
           IF(number_vertices.lt.100) WRITE(Iout,601) number_vertices
           IF(number_vertices.ge.100) WRITE(Iout,602) number_vertices
        else
           IPR=1 
           IF(number_vertices.lt.100) WRITE(Iout,603) number_vertices
           IF(number_vertices.ge.100) WRITE(Iout,604) number_vertices
        endif
      else
        if(IPR.EQ.0) then 
           IF(number_vertices.lt.100) WRITE(Iout,701) number_vertices
           IF(number_vertices.ge.100) WRITE(Iout,702) number_vertices
        else
           IPR=1 
           IF(number_vertices.lt.100) WRITE(Iout,703) number_vertices
           IF(number_vertices.ge.100) WRITE(Iout,704) number_vertices
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
      M=number_vertices/2+2
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
      CALL Unwind(M,IER,IT,ispiral,Spiralx,S,D,NMR,Group)! Unwind dual into spirals and check for uniqueness
      IF(IER.EQ.13) GO TO 13                
      K=0
      L=L+1                                              ! Spiral S is canonical      
      DO J=1,6
        IF(NMR(J).EQ.0) GO TO 16
        K=J
      enddo
C     Analyze dual matrix
   16  CALL DualAnalyze(M,D,IRhag5,IRhag6,
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
         do ia=1,number_vertices
           ka=0
            do ja=1,number_vertices
              if(IDA(Ia,Ja).eq.1) then
               ka=ka+1
               IC3(ia,ka)=ja
               endif
            end do
         end do
        Call HamiltonCyc(maxiter,Iout,nbatch,IC3,nhamcycle)
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
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 602  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 603  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NMR pattern',
     5 /1X,170('-')) 
 604  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' and gap the HOMO-LUMO gap in units of beta)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
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
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 702  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 703  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
 704  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain paramter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-')) 
      Return
      END subroutine spiral


      SUBROUTINE spwindup(IER,number_faces,MP,D1,S,JP,FreeRing)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D1(6,MMAX),S(MMAX),JP(12),FreeRing(6)
C This routine tries to find the spiral from a preset of
C three connected rings stored in S(1), S(2) and S(3) using
C information on ring fusions from the reduced dual matrix D1
C If pentagon is found MP increases by 1.
C IER=1 implies that no spiral found
C This algorithm is much faster than the one in version 4.2

      IER=0
C   Big loop from ring 4 onwards to the end of the spiral
C---- Start loop
      nsmall=1
      faces: Do i=4,number_faces

C       First collect next possible faces connected to the previous 
C       5- or 6-ring (nloop=5 or 6) numbered IP=S(i-1) 
C       and make sure it is not one of the previous rings
C       in the spiral
        IP=S(i-1)
        nloop=6
        nring=0
        if(D1(6,IP).eq.0) nloop=5
        loops: do j=1,nloop
C         Get all adjacent rings to previous one
          nr=D1(j,IP)
C         Make sure it is not one in the existing spiral
          do j1=I-2,1,-1
            if(nr.eq.S(j1)) cycle loops
          enddo
C       Collect them
        nring=nring+1
        FreeRing(nring)=nr
        enddo loops
          
C       Now it needs to be connected to a previous ring
C       Last 2 are not needed
        do k=nsmall,i-3
          KP=S(k)
          do k1=1,6
            do j2=1,nring
              nr=FreeRing(j2)
              if(D1(k1,KP).eq.nr) then
                nsmall=k
                S(i)=nr
                if(D1(6,nr).eq.0) then
                  MP=MP+1
                  JP(MP)=i
                endif
                cycle faces
              endif
            enddo
          enddo
        enddo

        if(S(i).eq.0) then
          IER=1   ! Spiral has dead end
          Return
        endif
      enddo faces

C  Finally success, spiral found
      Return
      END


      SUBROUTINE SpiralSearch(Iout,IRG55,IRG66,IRG56,
     1 NrA,NrB,NrC,NrD,NrE,NrF,JP,GROUP,spcount)
      use config
      use iso_c_binding
      IMPLICIT INTEGER (A-Z)
      DIMENSION NrA(EMAX),NrB(EMAX),NrC(EMAX),NrD(EMAX)
      DIMENSION NrE(EMAX),NrF(EMAX),NMR(6),JP(12)
      DIMENSION Spiral(12,NMAX),D1(6,MMax)
      DIMENSION SpiralT(12,MaxSpirals),SpiralF(MMAX,MaxSpirals)
      DIMENSION D(MMAX,MMAX),S(MMAX),FreeRing(6)
      CHARACTER*3 GROUP
      type(c_ptr) :: dg, fg, new_graph, dual_graph
      integer gen_rspi(12), gen_jumps(10)

C     Search for all spirals, set up first three rings then wind.
C     Start ring spiral algorithm. Quit after first successful spiral.
C     This subroutine has been modified from the original one of Fowler and 
C     Manolopoulos "An Atlas of Fullerenes" (Dover Publ., New York, 2006).           
C     It is used if dual matrix is already known. See subroutine Spiral for details.
C     number_vertices is the nuclearity of the fullerene.
C     JP contains the pentagon indices, S the ring numbers
C     SpiralT contains the pentagon indices, SpiralF the face numbers
C     Timing O(6 n_v n_f^2)

      number_faces=number_vertices/2+2
      ispiral=0
      maxgen_jumps=10
      WRITE (Iout,600)
      IF(number_vertices.lt.100)   
     1  WRITE(Iout,601) number_vertices,number_faces
      IF(number_vertices.ge.100.and.number_vertices.lt.1000)
     1  WRITE(Iout,602) number_vertices,number_faces
      IF(number_vertices.ge.1000.and.number_vertices.lt.10000)
     1  WRITE(Iout,632) number_vertices,number_faces
      IF(number_vertices.ge.10000)
     1  WRITE(Iout,633) number_vertices,number_faces
      do I=1,MMAX
        do J=1,MMAX
          D(I,J)=0
        enddo
      enddo
      do I=1,12
        do J=1,NMAX
          Spiral(I,J)=0
        enddo
        do J=1,MaxSpirals
          SpiralT(I,J)=0
        enddo
      enddo
      do I=1,6
        FreeRing(i)=0
      enddo

C     Set up dual matrix (adjacency matrix for faces)
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
C     Check if pentagons are in order (dual theorem)
      Do K1=1,12
       isumpent=0
      Do K2=1,number_faces
       isumpent=isumpent+D(K1,K2)
      enddo
       if(isumpent.ne.5) then
        Write(Iout,622) K1
        return
       endif
      enddo
C     Compress dual matrix
      Do K1=1,number_faces
      Do K2=1,6
       D1(K2,K1)=0
      enddo
      enddo
      Do K1=1,number_faces
      ncount=0
      Do K2=1,number_faces
       if(D(K1,K2).eq.1) then
        ncount=ncount+1
        D1(ncount,K1)=K2
       endif
      enddo
      enddo

C     Start searching for spiral from pentagon 1 to 12, then stop
C      and store the pentagon indices. Note this may still be an
C      unsuccessful ring spiral. Also throw duplicates out.
      nspiral=0
      nspiral5=0
      nspiralT=0

C Loop over all (5,5) fusions
      if(IRG55.eq.0) then
       write(Iout,610)
      else
       write(Iout,611) 2*IRG55
C----  outer loop, ensures both left and right spirals are included
       do I=1,2*IRG55
C     Get first two faces
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
C     Get third face, Note that J=I1 or J=I2 is excluded
C---- Inner loop, ensures that both 3rd faces connected to I1 and I2 are included
       do J=1,number_faces
        if(D(I1,J).eq.1.and.D(I2,J).eq.1) then
         MPent=2
         S(3)=J
C       Reset arrays
        do k=4,number_faces
         s(k)=0
        enddo
        do k=3,12
         JP(k)=0
        enddo
C     See if it is a pentagon, first 12 in D are
         if(J.le.12) then
          JP(3)=3
          MPent=3
         endif
C      Now search
          CALL spwindup(IER,number_faces,MPent,D1,S,JP,FreeRing)
C      Check if RSPIs are in order
         if(IER.eq.0) then
          nspiral=nspiral+1
          nspiral5=nspiral5+1
          nspiralT=nspiralT+1
C       Check if enough space
          If(nspiral.gt.MaxSpirals) then
           Write(Iout,626) nspiral,MaxSpirals
           nspiral=nspiral-1
           Go to 199
          endif
C       Now everything works fine
          do k=1,12
           SpiralT(k,nspiral)=JP(k)
          enddo 
          do k=1,number_faces
           SpiralF(k,nspiral)=S(k)
          enddo 
C      Jump count if spcount=0
        if(spcount.eq.0) then
         Write(Iout,635)
         go to 199
        endif
C       Delete identical spirals
          if(nspiral.gt.1) Call SpiralCheck(nspiral,SpiralT)
         endif
        endif
       enddo 
C---- End inner loop
       enddo 
C---- End outer loop
      endif
      nspiral55=nspiral
      nspiralT55=nspiralT

C Loop over all (5,6) fusions
C Dito, see above
      if(IRG56.le.0) then
       write(Iout,615)
      else
       write(Iout,612) 2*IRG56
      do I=1,2*IRG56
C     Reset arrays
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
       n5=1
       JP(1)=1
      else
       JP(1)=2
      endif
      do J=1,number_faces
       if(D(I1,J).eq.1.and.D(I2,J).eq.1) then
        S(3)=J
        MPent=1
        do k=4,number_faces
         s(k)=0
        enddo
        do k=2,12
         JP(k)=0
        enddo
        if(J.le.12) then
         JP(2)=3
         MPent=2
        endif
        CALL spwindup(IER,number_faces,MPent,D1,S,JP,FreeRing)
        if(IER.eq.0) then
         if(JP(1).eq.1) nspiral5=nspiral5+1
         nspiral=nspiral+1
         nspiralT=nspiralT+1
         If(nspiral.gt.MaxSpirals) then
          Write(Iout,627) nspiral,MaxSpirals
          nspiral=nspiral-1
          Go to 199
         endif
         do k=1,12
          SpiralT(k,nspiral)=JP(k)
         enddo 
         do k=1,number_faces
          SpiralF(k,nspiral)=S(k)
         enddo 
C      Jump count if spcount=0
        if(spcount.eq.0) then
         Write(Iout,635)
         go to 199
        endif
        if(nspiral.gt.1) Call SpiralCheck(nspiral,SpiralT)
       endif
       endif
      enddo 
      enddo 
      endif
      nspiral56=nspiral-nspiral55
      nspiralT56=nspiralT-nspiralT55
      
C Loop over all (6,6) fusions
C Dito, see above
      if(IRG66.eq.0) then
       write(Iout,616)
      else
       write(Iout,613) 2*IRG66
      do I=1,2*IRG66
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
       do J=1,number_faces
        if(D(I1,J).eq.1.and.D(I2,J).eq.1) then
          S(3)=J
          MPent=0
       do k=4,number_faces
        s(k)=0
       enddo
        do k=1,12
         JP(k)=0
        enddo
         if(J.le.12) then
          JP(1)=3
          MPent=1
         endif
       CALL spwindup(IER,number_faces,MPent,D1,S,JP,FreeRing)
       if(IER.eq.0) then
        nspiral=nspiral+1
        nspiralT=nspiralT+1
         If(nspiral.gt.MaxSpirals) then
          Write(Iout,628) nspiral,MaxSpirals
          nspiral=nspiral-1
          Go to 199
         endif
        do k=1,12
         SpiralT(k,nspiral)=JP(k)
        enddo
        do k=1,number_faces
         SpiralF(k,nspiral)=S(k)
        enddo 
C      Jump count if spcount=0
        if(spcount.eq.0) then
         Write(Iout,635)
         go to 199
        endif
        if(nspiral.gt.1) Call SpiralCheck(nspiral,SpiralT)
       endif
      endif
      enddo 
      enddo 
      endif
      spcount=1
C---- End of search
 199  if(nspiral.eq.0) then
        WRITE(Iout,630) nspiralT,6*number_vertices

        dg = new_graph(mmax,number_faces,D)
        fg = dual_graph(dg)
        call get_general_spiral(fg, gen_rspi, gen_jumps)
        call delete_graph(dg) 
        call delete_fullerene_graph(fg) 
        write(Iout,636) (gen_rspi(Mrspi),Mrspi=1,12)
        maxgj=2
        Do MJ=1,maxgen_jumps
         if(gen_jumps(MJ).eq.0) then
          maxgj=MJ-1
          go to 640
         endif
        enddo
  640   Write(Iout,637) maxgj/2
        write(Iout,638) (gen_jumps(MJ),MJ=1,maxgj)
        return
      else
        nspiral5sym=0
        do i=1,nspiral
          if(SpiralT(1,i).eq.1) nspiral5sym=nspiral5sym+1
        enddo
        WRITE(Iout,634) nspiral,nspiralT,6*number_vertices,
     1  nspiral5,nspiral5sym
      endif
      if(spcount.ne.0) then
       nspiral66=nspiral-nspiral55-nspiral56
       nspiralT66=nspiralT-nspiralT55-nspiralT56
       WRITE(Iout,631) nspiral55,nspiralT55,nspiral56,nspiralT56,
     1  nspiral66,nspiralT66
      endif
     
C     Now loop over with found spiral until success with
C     Fowler algorithm
      IT=1
      IPR=0
      nspfound=0
      Do 13 msp=1,nspiral
       Do I=1,number_faces
        S(I)=6
       enddo
       Do I=1,12
        S(SpiralT(I,msp))=5
        JP(I)=SpiralT(I,msp)
       enddo
       IER=0
       CALL Windup(number_faces,IPR,IER,S,D)     !      Wind up spiral into dual 
       IF(IER.ne.0) GO TO 13                     !      and check for closure 
       Do I=1,12 
        Spiral(I,1)=JP(I)
       enddo
        CALL Unwind(number_faces,IER,IT,ispiral,
     1   Spiral,S,D,NMR,Group)                    ! Unwind dual into spirals 
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
       if(number_faces.lt.1000) 
     1  WRITE(Iout,618) (SpiralF(I,msp),I=1,number_faces)
       if(number_faces.ge.1000) 
     1  WRITE(Iout,629) (SpiralF(I,msp),I=1,number_faces)
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
       Do I=1,number_faces
        S(I)=6
       enddo
       Do I=1,12
        S(JP(I))=5
       enddo
       WRITE(Iout,624)
       WRITE(Iout,625) (S(I),I=1,number_faces)
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
        WRITE(Iout,618) (SpiralF(I,msp),I=1,number_faces)
        return
        endif
       enddo
       endif
      else 
       WRITE(Iout,617) nspiral
      endif
 600  FORMAT(/1X,'Modified spiral algorithm Fowler and Manolopoulos',
     1 ' (An Atlas of Fullerenes, Dover Publ., New York, 2006)')
 601  FORMAT(1X,'Spiral for fullerene isomers of C',I2,':',
     1 ' (',I2,' faces)')
 602  FORMAT(1X,'Spiral for fullerene isomers of C',I3,':',
     1 ' (',I3,' faces)')
 603  FORMAT(1X,A3,9X,12I6)
 604  FORMAT(1X,131('-'),/1X,'Corresponding ring numbers:') 
 605  FORMAT(1X,A3,9X,12I6,2X,3(I4,' x',I4,:,','))
 607  Format(12(1X,I6))
 608  Format(1X,'Spiral is canonical')
 610  Format(1X,'This is an IPR fullerene, no (5,5) fusions to ',
     1 'loop over')
 611  Format(1X,'Loop over (5,5) fusions, ',I5,' max in total')
 612  Format(1X,'Loop over (5,6) fusions, ',I5,' max in total')
 613  Format(1X,'Loop over (6,6) fusions, ',I5,' max in total')
 614  Format(1X,I7,' RSPIs taken from spiral list to be analyzed:',/1X,
     2 'Point group   Ring spiral pentagon positions',
     3 42X,'NMR pattern (for fullerene in ideal symmetry)',
     4 /1X,131('-')) 
 615  Format(1X,'This is C20, no (5,6) fusions to loop over')
 616  Format(1X,'No (6,6) fusions to loop over')
 617  Format(1X,'--->',I7,' ring spirals found ',
     1 'but RSPIs detected were not valid')
 618  Format(30(1X,20(I5,'-'),/))
 619  Format(1X,'Spiral list of pentagon positions with ',
     1 'higher priority: (',I4,' spirals found)') 
 620  Format(1X,'Search ',I7,' spirals to produce canonical'
     1 ' list of faces:')
 621  Format(84('-'),/,12(1X,I6),/,84('-'))
 622  Format(1X,'Pentagons are not in sequence, stopped at pentagon ',
     1 I6)
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
 629  Format(100(1X,20(I5,'-'),/))
 630  Format(1X,'Failed to find ring spiral: ',I7,
     1 ' detected (maximum possible: ',I7,')',/1X,
     1 'This is a non-spiral fullerene: Entering the general spiral'
     1 ' algorithm',/1X,'Canonical General Spiral:')
 631  FORMAT(1X,I7,' distinct (55)      RSPIs found out of ',I7,
     1 /1X,I7,' distinct (56)/(65) RSPIs found out of ',I7,
     1 /1X,I7,' distinct (66)      RSPIs found out of ',I7)
 632  FORMAT(1X,'Spiral for fullerene isomers of C',I4,':',
     1 ' (',I4,' faces)')
 633  FORMAT(1X,'Spiral for fullerene isomers of C',I5,':',
     1 ' (',I5,' faces)')
 634  FORMAT(1X,I7,' distinct RSPIs found out of total',I7,
     1 ' (maximum possible: ',I7,')',/1X,I7,' spirals with a ',
     1 'pentagon start   (',I7,'  symmetry distinct)')
 635  FORMAT(1X,' Spiral found, avoid counting of spirals ',/2X,
     1 'WARNING: output contains only 1 spiral! Set ispsearch=2 ',
     1 'to get the complete spiral count.')
 636  FORMAT(1X,'RSPI : ',12I6)
 637  FORMAT(1X,'Number of jumps: ',I2)
 638  FORMAT(1X,'Jumps: ',10I6)
      Return
      END
     
      SUBROUTINE SpiralCheck(nspiral,SpiralT) 
      use config
      IMPLICIT Integer (A-Z)
      DIMENSION SpiralT(12,MaxSpirals),dif(12)
      Do I=1,nspiral-1
        sum=0
        Do J=1,12
          dif(J)=SpiralT(J,I)-SpiralT(J,nspiral)
          sum=sum+iabs(dif(J))
        enddo
        if(sum.eq.0) then
          nspiral=nspiral-1
          return
        endif
      enddo
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


      SUBROUTINE DualAnalyze(M,D,IRhag5,IRhag6,
     1 IFus5G,IDA,nelec,ndeg,sigmah,A,gap)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer D(MMAX,MMAX),Ddiag(MMAX),NR5(12),IDA(NMAX,NMAX)
      Integer IRhag5(0:5),IRhag6(0:6),IDG(NMAX)
      Dimension A(NMAX,NMAX),evec(NMAX),df(NMAX)
C     Analyze dual matrix and get pentagon and hexagon indices
C     and determine if molecule is open shell

      Tol=1.d-5

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
      IFus5G=IPentInd(IRhag5)

C     Now hexagon indices
      Do I=0,6
        IRhag6(I)=0
      enddo
      faces: do I=1,M
        IRcount=0
        IR5=Ddiag(I)
        if(IR5.eq.5) cycle faces
        do J=1,M
          JR5=Ddiag(J)
          If(JR5.ne.5.and.D(I,J).eq.1) then
            IRcount=IRcount+1
          endif
        enddo
        IRhag6(IRcount)=IRhag6(IRcount)+1
      enddo faces
C     Strain Parameter
      sigmah=HexInd(IRhag6,ihk)

C     Now produce adjacency matrix
      CALL DUAL(D,MMAX,IDA,IER)
      Do I=1,number_vertices
        df(I)=0.d0
        Do J=I,number_vertices
          A(I,J)=0.d0
          if(IDA(I,J).eq.1) A(I,J)=1.d0
          A(J,I)=A(I,J)
        enddo
      enddo
C Diagonalize without producing eigenvectors
      do i=1,nmax
        evec(i) = 0
      enddo
      call tred2l(A,number_vertices,NMAX,evec,df)
      call tqlil(evec,df,number_vertices,NMAX)
C Sort eigenvalues
      Do I=1,number_vertices
        e0=evec(I)
        jmax=I
        Do J=I+1,number_vertices
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
      Do I=2,number_vertices
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
      Noc=number_vertices/2
      Norb=0
      eigv: Do I=1,ieigv
        Iorb=I
        Norb=Norb+IDG(I)
        if(Norb.eq.Noc) then 
          gap=df(Iorb)-df(Iorb+1)
          ndeg=IDG(Iorb)
          nelec=ndeg*2
          exit eigv
        endif
        if(Norb.gt.Noc) then 
          gap=df(Iorb)-df(Iorb+1)
          ndeg=IDG(Iorb)
          nelec=(ndeg-Norb+Noc)*2
          exit eigv
        endif
      enddo eigv
      
      Return
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

      SUBROUTINE Unwind(M,IER,IT,ispiral,Spiral,S,D,NMR,GROUP)
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
C       Input: M: Number of Faces
C         1 specific spiral S which contains a sequence of numbers 
C             5 and 6 of the spiral
C         D which is the dual adjacency matrix


C-------------------------------------------------------------------
C----    STEP 1
C        In this step one loops over all 6N spirals and compares with
C         input spiral S. If identical because of symmetry, and there
C         should be |G| many (the order of the group), it is used
C         to obtain the face permutation stored in FP. For each 
C         symmetry equivalent spiral the permutation is in field P.
C         Getting all symmetry equivalent spirals results in the
C         order of the group stored in SO.
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
         Print*,' Order of the Group:',SO
C-------------------------------------------------------------------

C----    STEP 2
C         For this part the face permutations FP and dual adjacency
C          matrix are needed.
C         First get all the faces related to vertices (3 of them)
C          and edges (2 of them) and store them in fields V and E
C          respectively. Then get the vertex-, edge- and face-orbits
C          and assign to each orbit the site-symmetry. Count the number
C          of specific site-symmetries for each vertex-, edge- and face-orbit.
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


C       
C-------------------------------------------------------------------
C       Interlude: Get NMP pattern
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
C-------------------------------------------------------------------

C      Step 3: Go through the tree structure to get the point group
C       by using MS (site group order) and ORDER (full group order)
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

      SUBROUTINE DUAL(D,M,A,IER)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D(M,M),A(NMAX,NMAX)
      DIMENSION V(3,NMAX)
c     Given a fullerene dual adjacency matrix D, this subroutine 
c     constructs the corresponding fullerene adjacency matrix A. 
c     IER = 0 on return if the construction is successful. 	 
      I=0
      DO 3 L = 1,M 
        DO 2 K=1,L 
          IF (D(K,L).EQ.0) GO TO 2
          DO 1 J = 1,K
            IF (D(J,K).EQ.0.OR.D(J,L).EQ.0) GO TO 1
            I = I+1
            IF (I.GT.number_vertices) GO TO 1
            V(1,I) = J ! Associate the three mutually adjacent 
            V(2,I) = K ! dual vertices (fullerene faces) J,K,L 
            V(3,I) = L ! with fullerene vertex I 	 
 1        CONTINUE
 2      CONTINUE
 3    CONTINUE
      IER = I-number_vertices
      IF (IER .NE. 0) RETURN ! D contains IER > 0 separating triangles 
      DO 7 J = 1,number_vertices           ! and is therefore NOT a fullerene dual 
        DO 6 I = 1,J
          K = 0
          DO 5 JJ = 1,3
            DO 4 II = 1,3
              IF(V(II,I).EQ.V(JJ,J)) K = K+1
 4          CONTINUE
 5        CONTINUE
          IF (K.EQ.2) THEN
            A(I,J)=1   ! Fullerene vertices I and J are adjacent 
            A(J,I)=1   ! if they have 2 dual vertices in common
          ELSE
            A(I,J)=0
            A(J,I)=0
          ENDIF
 6      CONTINUE
 7    CONTINUE
      RETURN
      END SUBROUTINE DUAL
