      SUBROUTINE Isomers(NAtom,Nfaces,Nedges,N,IPR,IOUT,
     1 maxAiso,iham,ichk,IDA,A,chkname)
C Information on number of isomers with or without fulfilling the
C the IPR rule. The routine also calls SPIRAL using the subroutines
C written by P. W. Fowler and D. E. Manopoulus, "An Atlas of Fullerenes"
C (Dover Publ., New York, 2006), which gives information on the
C isomers point group, pentagon ring spiral indices and NMR pattern. 
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer Isonum(119),IsonumIPR(123),IDA(NAtom,NAtom)
      Character*1 fstart,fnum1
      Character*2 fnum2
      Character*3 fnum,fnum3,GROUP,ident
      Character*9 fend
      Character*11 Isoc(52),IsocIPR(28)
      Character*20 chkname
      Character*13 dbdir
      Character*31 databasefile
      Logical lexist
      Dimension A(NAtom,NAtom)
      Data Isonum/1,0,1,1,2,3,6,6,15,17,40,45,89,116,199,
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
      Data IsonumIPR/1,0,0,0,0,1,1,1,2,5,7,9,24,19,35,46,86,134,187,
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
       Data Isoc/'2295793276','2477017558','2648697036','2854536850',
     *           '3048609900','3282202941','3501931260','3765465341',
     *           '4014007928','4311652376','4591045471','4926987377',
     *           '5241548270','5618445787','5972426835','6395981131',
     *           '6791769082','7267283603','7710782991','8241719706',
     *          '8738236515','9332065811','9884604767','10548218751',
     *       '11164542762','11902015724','12588998862','13410330482',
     *       '14171344797','15085164571','15930619304','16942010457',
     *       '17880232383','19002055537','20037346408','21280571390',
     *       '22426253115','23796620378','25063227406','26577912084',
     *       '27970034826','29642262229','31177474996','33014225318',
     *       '34705254287','36728266430','38580626759','40806395661',
     *       '42842199753','45278616586','47513679057','50189039868'/
      Data IsocIPR/'2266138871','2435848971','2614544391','2808510141',
     *       '3009120113','3229731630','3458148016','3704939275',
     *       '3964153268','4244706701','4533465777','4850870260',
     *       '5178120469','5531727283','5900369830','6299880577',
     *       '6709574675','7158963073','7620446934','8118481242',
     *       '8636262789','9196920285','9768511147','10396040696',
     *       '11037658075','11730538496','12446446419','13221751502'/

      IPRERR=0
C     Number of Isomers
      ISOMAX=256
      If(N.le.256) then
C     Both values fit 32bit signed integer
       M1=N/2-9
       M2=N/2-29
        if(M2.gt.0) then
         isomIPR=IsonumIPR(M2)
        else
         isomIPR=0
        endif
       Write(Iout,1000) Isonum(M1),isomIPR
       AisoNIPR=dfloat(Isonum(M1))
       AisoIPR=dfloat(isomIPR)
      endif

C     IPR value only fits 32bit signed integer
        If(N.le.304.and.N.gt.256) then
           fac1=10.d0**3.502347
         M1=N/2-128
         M2=N/2-29
         Write(Iout,1001) isoc(M1),IsonumIPR(M2)
        endif

C     Both values do not fit 32bit signed integer
C      Output in characters
          If(N.le.360.and.N.gt.304) then
           M1=N/2-128
           M2=N/2-152
         Write(Iout,1002) isoc(M1),IsocIPR(M2)
          endif

C     Both values not known, estimate by polynomial
C      fitted to asymptodic 
          If(N.gt.256) then
           fac1=10.d0**3.502347
           fac2=10.d0**1.503692
           AisoNIPR=fac1*(dfloat(N)/6.d1)**9.250947d0
           AisoIPR=fac2*(dfloat(N)/6.d1)**11.08800d0
           If(N.gt.360) Write(Iout,1003) AisoNIPR,AisoIPR
          endif
 
      if(IPR.eq.0) then
       if(AisoNIPR.gt.dfloat(maxAiso)) then
        Write(Iout,1004) maxAiso
        Return
       endif
      else
       if(AisoIPR.gt.dfloat(maxAiso)) then
        Write(Iout,1004) maxAiso
        Return
       endif
      endif

      if(IPR.eq.1.and.N.lt.60) IPRERR=1
      if(IPR.eq.1.and.(N.gt.60.and.N.lt.70)) IPRERR=1
      if(IPRERR.eq.1) then
       Write(Iout,1007)
       return
      endif

C SPIRAL uses the subroutines written by Fowler and Manopoulus
      If(ichk.gt.0) then
       Write(Iout,1006)
       CALL SpiralRestart(Natom,Nfaces,Nedges,N,IPR,Iout,Isonum,
     1 IsonumIPR,iham,IDA,A,chkname)
       return
      endif

C Check if database can be taken instead
      If((IPR.eq.0.and.N.le.100).or.(IPR.eq.1.and.N.le.120)) then
      dbdir='database/All/'
      if(IPR.eq.1) dbdir='database/IPR/'
      fend='.database'
      fstart='c'
      if(N.lt.100) then
      fnum1='0'
       write(fnum2,'(I2)') N
       fnum=fnum1//fnum2
      else
       write(fnum,'(I3)') N
      endif
      fnum3='all'
      if(IPR.eq.1) fnum3='IPR'
      databasefile=dbdir//fstart//fnum//fnum3//fend
      if(IPR.eq.0) then
       Write(Iout,1008) databasefile
      else
       Write(Iout,1009) databasefile
      endif
      inquire(file=databasefile,exist=lexist)
      if(lexist.neqv..True.) then
       Write(Iout,1010) databasefile
       Go to 99
      else
       Write(Iout,1011) databasefile
       call Printdatabase(N,Iout,databasefile)
      endif
      return
      endif

C Produce list from ring spiral algorithm
  99  If(IPR.ge.0) then
       Write(Iout,1005)
       CALL Spiral(Natom,Nfaces,Nedges,N,IPR,Iout,Isonum,
     1 IsonumIPR,iham,IDA,A)
      endif

 1000 Format(/1X,'Number of possible fullerene isomers: ',I10,
     1 ' (IPR isomers: ',I10,')')
 1001 Format(/1X,'Number of possible fullerene isomers: ',A11,
     1 ' (IPR isomers: ',I10,')')
 1002 Format(/1X,'Number of possible fullerene isomers: ',A11,
     1 ' (IPR isomers: ',A11,')')
 1003 Format(/1X,'Exact number of isomers not known for such a'
     1 ' large fullerene, projected polynomial value:',D12.4,
     1 ' (general), ',D12.4,' (IPR)',
     1 /2X,'obtained from fit of points between Nc=250 and 360)')
 1004 Format(/1X,'Number of Isomers larger than max value of ',I8)
 1005 Format(/1X,'Enter Spiral code for a general list of all ',
     1 'possible isomers (IPR or not depending on input)')
 1006 Format(/1X,'RESTART isomer file from previous run')
 1007 Format(/1X,'Zero IPR isomers -> Return')
 1008 Format(1X,'Search for file: ',A29,' in general isomer list')
 1009 Format(1X,'Search for file: ',A29,' in IPR isomer list')
 1010 Format(1X,'Filename ',A29,' in database not found: ',
     1 'Do it the hard way')
 1011 Format(1X,'Print from file ',A29,' in database')
      Return
      END
 
      SUBROUTINE Printdatabase(N,Iout,databasefile)
      Character*31 databasefile
      Character*1 Text(200),Textind
      Open(UNIT=4,FILE=databasefile,STATUS='old',FORM='FORMATTED')
       Textind=' '
       Nlimit=1000000000
       Read(4,*) IN,IP,IH
       if(IN.ne.N) then
        Write(Iout,1002) IN,N
        return
       endif
       Write(Iout,1000) IN,IP,IH
      if(IH.eq.0) then
      if(IP.EQ.0) then
         IF(N.lt.100) WRITE(Iout,601) N
         IF(N.ge.100) WRITE(Iout,602) N
      else
         IF(N.lt.100) WRITE(Iout,603) N
         IF(N.ge.100) WRITE(Iout,604) N
      endif
      else
      if(IP.EQ.0) then
         IF(N.lt.100) WRITE(Iout,701) N
         IF(N.ge.100) WRITE(Iout,702) N
      else
         IF(N.lt.100) WRITE(Iout,703) N
         IF(N.ge.100) WRITE(Iout,704) N
      endif
      endif

       do I =1,Nlimit
        Read(4,1001,Err=99,end=99) (Text(J),J=1,200)
        do k=1,200
         l=200-k+1
         if(Text(l).ne.Textind) then
          NChar=l
          go to 10
         endif
        enddo
  10    Write(Iout,1001) (Text(J),J=1,NChar)
       enddo
  99  Close(unit=4)
 1000 Format(/1X,I10,2I2)
 1001 Format(200A1)
 1002 Format(/1X,'Atom number ',I5,' not identical to that on file: ',
     1 I5)
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
