      SUBROUTINE Isomers(IPR,isearch,IN,IOUT,iham,ihamstat,
     1 isomerl,isomerh,ichk,IDA,A,filename)
C Information on number of isomers with or without fulfilling the
C the IPR rule. The routine also calls SPIRAL using the subroutines
C written by P. W. Fowler and D. E. Manopoulus, "An Atlas of Fullerenes"
C (Dover Publ., New York, 2006), which gives information on the
C isomers point group, pentagon ring spiral indices and NMR pattern. 
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer Isonum(119),IsonumIPR(123),IDA(Nmax,Nmax)
      Character*1 fstart,fnum1
      Character*2 fnum2
      Character*3 fnum,fnum3
      Character*9 fend
      Character*12 Isoc(72),IsocIPR(48)
      Character*50 filename
      Character*50 chkname
      Character*13 dbdir
      Character*50 databasefile
      Logical lexist
      Dimension A(Nmax,Nmax)
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
     *       '42842199753','45278616586','47513679057','50189039868',
     *       '52628839448','55562506886','58236270451','61437700788',
     *       '64363670678','67868149215','71052718441','74884539987',
     *       '78364039771','82532990559','86329680991','90881152117',
     *      '95001297565','99963147805','104453597992','109837310021',
     *    '114722988623','120585261143','125873325588','132247999328'/

      Data IsocIPR/'2266138871','2435848971','2614544391','2808510141',
     *       '3009120113','3229731630','3458148016','3704939275',
     *       '3964153268','4244706701','4533465777','4850870260',
     *       '5178120469','5531727283','5900369830','6299880577',
     *       '6709574675','7158963073','7620446934','8118481242',
     *       '8636262789','9196920285','9768511147','10396040696',
     *       '11037658075','11730538496','12446446419','13221751502',
     *       '14010515381','14874753568','15754940959','16705334454',
     *       '17683643273','18744292915','19816289281','20992425825',
     *       '22186413139','23475079272','24795898388','26227197453',
     *       '27670862550','29254036711','30852950986','32581366295',
     *       '34345173894','36259212641','38179777473','40286153024'/

      Write(Iout,1012)
      isocount=0
      chkname=trim(filename)//".chkpnt"
      IPRERR=0
C     Number of Isomers
      ISOMAX=256
      If(number_vertices.le.256) then
C     Both values fit 32bit signed integer
       M1=number_vertices/2-9
       M2=number_vertices/2-29
        if(M2.gt.0) then
         isomIPR=IsonumIPR(M2)
        else
         isomIPR=0
        endif
       Write(Iout,1000) Isonum(M1),isomIPR
       isocount=Isonum(M1)
       AisoNIPR=dfloat(Isonum(M1))
       AisoIPR=dfloat(isomIPR)
      endif

C     IPR value only fits 32bit signed integer
        If(number_vertices.le.304.and.number_vertices.gt.256) then
           fac1=10.d0**3.502347
         M1=number_vertices/2-128
         M2=number_vertices/2-29
         Write(Iout,1001) isoc(M1),IsonumIPR(M2)
        endif

C     Both values do not fit 32bit signed integer
C      Output in characters
          If(number_vertices.le.400.and.number_vertices.gt.304) then
           M1=number_vertices/2-128
           M2=number_vertices/2-152
         Write(Iout,1002) isoc(M1),IsocIPR(M2)
          endif

C     Both values not known, estimate by polynomial
C      fitted to asymptodic 
          If(number_vertices.gt.256) then
           fac1=10.d0**3.502347
           fac2=10.d0**1.503692
           AisoNIPR=fac1*(dfloat(number_vertices)/6.d1)**9.250947d0
           AisoIPR=fac2*(dfloat(number_vertices)/6.d1)**11.08800d0
           If(number_vertices.gt.400) Write(Iout,1003) AisoNIPR,AisoIPR
          endif

C     Limit number of cycles
      maxRSI=maxRS
      if(iham.ne.0) maxRSI=maxRS/10
      if(IPR.eq.0) then
       if(AisoNIPR.gt.dfloat(maxRSI)) then
        Write(Iout,1004) maxRSI
        if(isearch.ne.0) then
         Write(Iout,1013) isearch
         go to 99
        endif
        Return
       endif
      else
       if(AisoIPR.gt.dfloat(maxRSI)) then
        Write(Iout,1004) maxRSI
        if(isearch.ne.0) then
         Write(Iout,1013) isearch
         go to 99
        endif
        Return
       endif
      endif

      if(IPR.eq.1.and.number_vertices.lt.60) IPRERR=1
      if(IPR.eq.1.and.(number_vertices.gt.60.and.number_vertices.lt.70))
     1   IPRERR=1
      if(IPRERR.eq.1) then
       Write(Iout,1007)
       return
      endif

C SPIRAL uses the subroutines written by Fowler and Manopoulus
      If(ichk.gt.0) then
       Write(Iout,1006)
       CALL SpiralRestart(IPR,Iout,Isonum,IsonumIPR,iham,IDA,A,chkname)
       return
      endif

C Check if database can be taken instead
      If((IPR.eq.0.and.number_vertices.le.LimitAll).or.
     1 (IPR.eq.1.and.number_vertices.le.LimitIPR)) then
      dbdir='database/All/'
      if(IPR.eq.1) dbdir='database/IPR/'
      fend='.database'
      fstart='c'
      if(number_vertices.lt.100) then
        fnum1='0'
         write(fnum2,'(I2)') number_vertices
         fnum=fnum1//fnum2
      else
        write(fnum,'(I3)') number_vertices
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
       call Printdatabase(Iout,iham,ihamstat,isocount,
     1  isomerl,isomerh,databasefile)
      endif
      return
      endif


C Produce list from ring spiral algorithm
  99  If(IPR.ge.0.or.isearch.ne.0) then
        if(isearch.eq.0) then
          Write(Iout,1005)
          CALL Spiral(IPR,Iout,Isonum,IsonumIPR,iham,IDA,A)
        else
          Write(Iout,1015)
          Call SpiralFind(IPR,isearch,In,Iout,IDA,A)
        endif
      endif

 1000 Format(1X,'Number of possible fullerene isomers: ',I10,
     1 ' (IPR isomers: ',I10,')')
 1001 Format(1X,'Number of possible fullerene isomers: ',A12,
     1 ' (IPR isomers: ',I10,')')
 1002 Format(1X,'Number of possible fullerene isomers: ',A12,
     1 ' (IPR isomers: ',A11,')')
 1003 Format(1X,'Exact number of isomers not known for such a'
     1 ' large fullerene, projected polynomial value:',D12.4,
     1 ' (general), ',D12.4,' (IPR)',
     1 /2X,'obtained from fit of points between Nc=250 and 360)')
 1004 Format(1X,'Number of Isomers larger than max value of ',I8,
     1 ' to produce isomer list')
 1005 Format(1X,'Enter Spiral code for a general list of all ',
     1 'possible isomers (IPR or not depending on input)')
 1006 Format(1X,'RESTART isomer file from previous run')
 1007 Format(1X,'Zero IPR isomers -> Return')
 1008 Format(1X,'Search for file: ',A50,' in general isomer list')
 1009 Format(1X,'Search for file: ',A50,' in IPR isomer list')
 1010 Format(1X,'Filename ',A50,' in database not found: ',
     1 'Do it the hard way')
 1011 Format(1X,'Print from file: ',A50,' in database')
 1012 Format(/1X,'Data for isomer numbers from House of Graphs website:'
     1,' http://hog.grinvin.org/')
 1013 Format(1X,'Search for nearest icosahedral ring spiral indices',
     1' instead with variation V= ',I2)
 1015 Format(1X,'Enter Spiral code for search of possible RSPIs')
      Return
      END
 
      SUBROUTINE Printdatabase(Iout,iham,ihamstat,isocount,
     1  isomerl,isomerh,databasefile)
C---------------------------------------------------------------------
C  This routine reads from the database using a specific format, and prints
C  various parameters for each isomer, that is
C  IN,IP,IH  (IN: number of vertices, IP: 0 (all isomers), 1 (only IPR)
C             IH: 0 (no) 1 (yes) for number of Hamilton cycles
C  Group,(RSPI(i),I=1,12),(PNI(I),I=0,4),(HNI(I),I=0,5),NeHOMO,NedegHOMO,HLgap,
C         ncycHam,(INMR(I),I=1,6)
C  Symmetry group, face spiral pentagon indices RSPI, pentagon indices PNI,
C    hexagon indices HNI, NeHOMO number of electrons in HOMO,
C    NedegHOMO degeneracy of HOMO, HLgap HOMO-LUMO gap, ncycHam number of
C    Hamilton cycles, INMR NMR pattern 
C---------------------------------------------------------------------
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Character*50 databasefile
      CHARACTER*3  Group
      CHARACTER*6  Occup
      Integer hamlow,hamhigh,hamlowIPR,hamhighIPR
      Integer RSPI(12),PNI(0:5),HNI(0:6),INMR(6)
      Integer D(MMAX,MMAX),S(MMAX),IDA(NMAX,NMAX),IC3(NMAX,3)
      Integer IsoExceptl(100),IsoExcepth(100),nbarval(nbardim)

      if(isomerl.ne.1) Write(Iout,1011) isomerl
      if(isomerh.ne.Nisoloop) Write(Iout,1012) isomerh
      nbatch=0
      nexceptionl=0
      nexceptionh=0
      maxiter=1000000000
      Open(UNIT=4,FILE=databasefile,STATUS='old',ACTION='Read',
     1  FORM='FORMATTED')
       Read(4,1003) IN,IP,IH
       if(IN.ne.number_vertices) then
        Write(Iout,1002) IN,number_vertices
        return
       endif
       number_faces=number_vertices/2+2
       Write(Iout,1000) IN,IP,IH
      if(IH.eq.0) then
       if(ihamstat.ne.0) then
        Write(Iout,1013)
        ihamstat=0
       endif
      if(IP.EQ.0) then
         IF(number_vertices.lt.100) WRITE(Iout,601) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,602) number_vertices
      else
         IF(number_vertices.lt.100) WRITE(Iout,603) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,604) number_vertices
      endif
      else
       vertnum=dfloat(number_vertices)
       ahamlow=5.*2.**(vertnum/10.-1.)
       nhamlow=int(ahamlow)
       an=2.*3.**(vertnum/20.-1.)+1.
       ahamhigh=ahamlow*an
       nhamhigh=int(ahamlow*an)
       Write(Iout,619) nhamlow,nhamhigh
       if(ihamstat.ne.0) then
        do i=1,nbardim
         nbarval(i)=0
        enddo
        nhamdif=nhamhigh-nhamlow
        nbars=int(.9*nbardim+1)
        nwidth=nhamdif/nbars
        if(nwidth.lt.10) then
         ihamstat=0
         Write(Iout,1014) nwidth
        else
         Write(Iout,1017) nwidth
        endif
       endif
      if(IP.EQ.0) then
         IF(number_vertices.lt.100) WRITE(Iout,701) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,702) number_vertices
      else
         IF(number_vertices.lt.100) WRITE(Iout,703) number_vertices
         IF(number_vertices.ge.100) WRITE(Iout,704) number_vertices
      endif
      endif

       do I=1,5
        PNI(I)=0
       enddo
       PNI(0)=12
       do I=0,2
        HNI(I)=0
       enddo
       IFus5G=0

       IFus5Glow=100000
       IFus5Ghigh=0
       sigmahlow=1.d10
       sigmahhigh=0.d0
       hamlow=1000000000
       hamhigh=0
       hamlowIPR=1000000000
       hamhighIPR=0
       no5ringstart=0

       do J =1,Nisoloop
        if(IP.eq.0) then
         if(IH.eq.1) then
C Case 1 All isomers with Hamiltonian cycles IP=0 IH=1
          Read(4,1004,ERR=99,end=99) Group,(RSPI(i),I=1,12),
     1     (PNI(I),I=0,4),(HNI(I),I=0,5),NeHOMO,NedegHOMO,HLgap,
     2     ncycHam,(INMR(I),I=1,6)
          if(ncycHam.lt.nhamlow) then
           nexceptionl=nexceptionl+1
           IsoExceptl(nexceptionl)=J  
          endif
          if(ncycHam.gt.nhamhigh) then
           nexceptionh=nexceptionh+1
           IsoExcepth(nexceptionh)=J  
          endif
          if(J.ge.isomerl) then
          if(J.gt.isomerh) go to 99
          if(RSPI(1).ne.1) no5ringstart=no5ringstart+1
          PNI(5)=12-PNI(0)-PNI(1)-PNI(2)-PNI(3)-PNI(4)
          IFus5G=IPentInd(PNI)
          HNI(6)=number_vertices/2-10
     1           -HNI(0)-HNI(1)-HNI(2)-HNI(3)-HNI(4)-HNI(5)
          sigmah=HexInd(HNI,ihk)
          if(2*NedegHOMO.eq.NeHOMO) then
           Occup='closed'
          else
           Occup='open  '
          endif
          nmrloop=2
          if(INMR(3).ne.0) nmrloop=4
          if(INMR(5).ne.0) nmrloop=6
          WRITE(Iout,608) J,GROUP,(RSPI(i),I=1,12),(PNI(I),I=0,5),
     1     IFus5G,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,HLgap,
     2     Occup,ncycHam,(INMR(I),I=1,nmrloop)
          if(ihamstat.ne.0) then
           ifield=(ncycHam-nhamlow)/nwidth+1           
           if(ifield.gt.nbardim) ifield=nbardim
           if(ifield.le.0) ifield=1
           nbarval(ifield)=nbarval(ifield)+1
          endif
          if(IFus5G.le.IFus5Glow) then
           IFus5Glow=IFus5G
           IFusL=J
          endif
          if(IFus5G.ge.IFus5Ghigh) then
           IFus5Ghigh=IFus5G
           IFusH=J
          endif
          if(sigmah.le.sigmahlow) then
           sigmahlow=sigmah
           ISigmaL=J
          endif
          if(sigmah.ge.sigmahhigh) then
           sigmahhigh=sigmah
           ISigmaH=J
          endif
          if(ncycham.le.hamlow) then
           hamlow=ncycham
           islow=J
          endif
          if(ncycham.ge.hamhigh) then
           hamhigh=ncycham
           ishigh=J
          endif
          if((number_vertices.eq.60.or.number_vertices.ge.70).
     1     and.IFus5G.eq.0) then
           if(ncycham.le.hamlowIPR) then
            hamlowIPR=ncycham
            islowIPR=J
           endif
           if(ncycham.ge.hamhighIPR) then
            hamhighIPR=ncycham
            ishighIPR=J
           endif
         endif
         endif

         else
C Case 2 All isomers without Hamiltonian cycles IP=0 IH=0
          ihamstat=0
          Read(4,1007,ERR=99,end=99) Group,(RSPI(i),I=1,12),
     1    (PNI(I),I=0,4),(HNI(I),I=0,5),NeHOMO,NedegHOMO,HLgap,
     2    (INMR(I),I=1,6)
          if(J.ge.isomerl) then
          if(J.gt.isomerh) go to 99
          if(RSPI(1).ne.1) no5ringstart=no5ringstart+1
          PNI(5)=12-PNI(0)-PNI(1)-PNI(2)-PNI(3)-PNI(4)
          IFus5G=IPentInd(PNI)
          HNI(6)=number_vertices/2-10
     1           -HNI(0)-HNI(1)-HNI(2)-HNI(3)-HNI(4)-HNI(5)
          sigmah=HexInd(HNI,ihk)
          if(2*NedegHOMO.eq.NeHOMO) then
           Occup='closed'
          else
           Occup='open  '
          endif
          nmrloop=2
          if(INMR(3).ne.0) nmrloop=4
          if(INMR(5).ne.0) nmrloop=6
          if(iham.eq.0) then
           WRITE(Iout,607) J,GROUP,(RSPI(i),I=1,12),(PNI(I),I=0,5),
     1      IFus5G,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,HLgap,
     2      Occup,(INMR(I),I=1,nmrloop)
          else
C   Calculate Hamiltonian cycles
C   Get dual and then adjacency matrix
           do K1=1,number_faces
           do K2=1,number_faces
            D(K1,K2)=0
           enddo
           enddo
           Do K=1,number_faces
            S(K)=6
           enddo
           Do K=1,12
            S(RSPI(K))=5
           enddo
           IPRS=0
           IER=0
           CALL Windup(number_faces,IPRS,IER,S,D) ! Wind up spiral into dual
           IF(IER.gt.0) then
            WRITE(Iout,1010) IER
            return
           endif
           CALL DUAL(D,MMAX,IDA,IER)
C  Now do Hamiltonian cycles
C create IC3 from IDA        
           do ia=1,number_vertices
             ka=0
             do ja=1,number_vertices
               if(IDA(Ia,Ja).eq.1) then
                 ka=ka+1
                 IC3(ia,ka)=ja
               endif
             end do
           end do
           Call HamiltonCyc(maxiter,Iout,nbatch,IC3,ncycham)
           WRITE(Iout,608) J,GROUP,(RSPI(i),I=1,12),(PNI(I),I=0,5),
     1      IFus5G,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,HLgap,
     2      Occup,ncycHam,(INMR(I),I=1,nmrloop)
           if(ncycham.le.hamlow) then
            hamlow=ncycham
            islow=J
           endif
           if(ncycham.ge.hamhigh) then
            hamhigh=ncycham
            ishigh=J
           endif
          endif
          if(IFus5G.le.IFus5Glow) then
           IFus5Glow=IFus5G
           IFusL=J
          endif
          if(IFus5G.ge.IFus5Ghigh) then
           IFus5Ghigh=IFus5G
           IFusH=J
          endif
          if(sigmah.le.sigmahlow) then
           sigmahlow=sigmah
           ISigmaL=J
          endif
          if(sigmah.ge.sigmahhigh) then
           sigmahhigh=sigmah
           ISigmaH=J
          endif
         endif
         endif

        else
         if(IH.eq.1) then
C Case 3 IPR isomers with Hamiltonian cycles IP=1 IH=1
          ihamstat=0
          Read(4,1008,ERR=99,end=99) Group,(RSPI(i),I=1,12),
     1     (HNI(I),I=3,5),NeHOMO,NedegHOMO,HLgap,ncycHam,
     2     (INMR(I),I=1,6)
          if(ncycHam.lt.nhamlow) then
           nexceptionl=nexceptionl+1
           IsoExceptl(nexceptionl)=J  
          endif
          if(ncycHam.gt.nhamhigh) then
           nexceptionh=nexceptionh+1
           IsoExcepth(nexceptionh)=J  
          endif
          if(J.ge.isomerl) then
          if(J.gt.isomerh) go to 99
          HNI(6)=number_vertices/2-10-HNI(3)-HNI(4)-HNI(5)
          sigmah=HexInd(HNI,ihk)
          if(2*NedegHOMO.eq.NeHOMO) then
           Occup='closed'
          else
           Occup='open  '
          endif
          nmrloop=2
          if(INMR(3).ne.0) nmrloop=4
          if(INMR(5).ne.0) nmrloop=6
          WRITE(Iout,608) J,GROUP,(RSPI(i),I=1,12),(PNI(I),I=0,5),
     1     IFus5G,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,HLgap,
     2     Occup,ncycHam,(INMR(I),I=1,nmrloop)
          if(sigmah.le.sigmahlow) then
           sigmahlow=sigmah
           ISigmaL=J
          endif
          if(sigmah.ge.sigmahhigh) then
           sigmahhigh=sigmah
           ISigmaH=J
          endif
          if(ncycham.le.hamlow) then
           hamlow=ncycham
           islow=J
          endif
          if(ncycham.ge.hamhigh) then
           hamhigh=ncycham
           ishigh=J
          endif
          endif

         else
C Case 4 IPR isomers without Hamiltonian cycles IP=1 IH=0
          ihamstat=0
          Read(4,1009,ERR=99,end=99) Group,(RSPI(i),I=1,12),
     1     (HNI(I),I=3,5),NeHOMO,NedegHOMO,HLgap,(INMR(I),I=1,6)
          HNI(6)=number_vertices/2-10-HNI(3)-HNI(4)-HNI(5)
          if(J.ge.isomerl) then
          if(J.gt.isomerh) go to 99
          sigmah=HexInd(HNI,ihk)
          if(2*NedegHOMO.eq.NeHOMO) then
           Occup='closed'
          else
           Occup='open  '
          endif
          nmrloop=2
          if(INMR(3).ne.0) nmrloop=4
          if(INMR(5).ne.0) nmrloop=6
          if(iham.eq.0) then
           WRITE(Iout,607) J,GROUP,(RSPI(i),I=1,12),(PNI(I),I=0,5),
     1      IFus5G,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,HLgap,
     2      Occup,(INMR(I),I=1,nmrloop)
          else
C   Calculate Hamiltonian cycles
C   Get dual and then adjacency matrix
           do K1=1,number_faces
           do K2=1,number_faces
            D(K1,K2)=0
           enddo
           enddo
           Do K=1,number_faces
            S(K)=6
           enddo
           Do K=1,12
            S(RSPI(K))=5
           enddo
           IPRS=0
           IER=0
           CALL Windup(number_faces,IPRS,IER,S,D) ! Wind up spiral into dual
           IF(IER.gt.0) then
            WRITE(Iout,1010) IER
            return
           endif
           CALL DUAL(D,MMAX,IDA,IER)
C create IC3 from IDA        
           do ia=1,number_vertices
             ka=0
             do ja=1,number_vertices
               if(IDA(Ia,Ja).eq.1) then
                 ka=ka+1
                 IC3(ia,ka)=ja
               endif
             end do
           end do
           Call HamiltonCyc(maxiter,Iout,nbatch,IC3,ncycham)
           WRITE(Iout,608) J,GROUP,(RSPI(i),I=1,12),(PNI(I),I=0,5),
     1      IFus5G,(HNI(I),I=0,6),sigmah,NeHOMO,NedegHOMO,HLgap,
     2      Occup,ncycHam,(INMR(I),I=1,nmrloop)
           if(sigmah.le.sigmahlow) then
            sigmahlow=sigmah
            ISigmaL=J
           endif
           if(sigmah.ge.sigmahhigh) then
            sigmahhigh=sigmah
            ISigmaH=J
           endif
          endif
           if(sigmah.le.sigmahlow) then
           sigmahlow=sigmah
           ISigmaL=J
          endif
          if(sigmah.ge.sigmahhigh) then
           sigmahhigh=sigmah
           ISigmaH=J
          endif
         endif
         endif

        endif
       enddo

C Final statistics
  99  if(IP.eq.0) then
        WRITE(Iout,611) IFus5Glow,IFusL,IFus5Ghigh,IFusH,
     1   sigmahlow,ISigmaL,sigmahhigh,ISigmaH
        if(IH.eq.1) then
         WRITE(Iout,609) hamlow,islow,hamhigh,ishigh
          if(number_vertices.eq.60.or.number_vertices.ge.70) then
           WRITE(Iout,610) hamlowIPR,islowIPR,hamhighIPR,ishighIPR
          endif
           if(nexceptionl.ne.0) then
            Write(iout,616) nexceptionl
            Write(iout,618) (IsoExceptl(J),J=1,nexceptionl)
           endif
           if(nexceptionh.ne.0) then
            Write(iout,617) nexceptionh
            Write(iout,618) (IsoExcepth(J),J=1,nexceptionh)
           endif
        endif
      else
        WRITE(Iout,612) sigmahlow,ISigmaL,sigmahhigh,ISigmaH
       if(IH.eq.1) WRITE(Iout,610) hamlow,islow,hamhigh,ishigh
      endif
      if(no5ringstart.eq.0) then
       Write(Iout,614)
      else
       Write(Iout,615) no5ringstart
      endif
      WRITE(Iout,613) 
      if(ihamstat.ne.0) then
       write(Iout,1015)
        nwidthhalf=nwidth/2
        ibars=0
        mem=1
       do i=1,nbardim
        nhamcount=nbarval(i)
        if(nhamcount.ne.0) then
         jhamcyc=i*nwidth+nhamlow-nwidthhalf
         barnormal=0
         if(isocount.ne.0) then
          if(mem.eq.1) then
           nlowest=jhamcyc
           mem=0
          endif
          ibars=ibars+1
          barnormal=dfloat(nhamcount)/dfloat(isocount)
         endif
         write(Iout,1016) i,jhamcyc,nhamcount,barnormal
        endif
       enddo
         write(Iout,1018) ibars,nwidth,nlowest
      endif
      Close(unit=4)
 1000 Format(/1X,I10,2I2)
 1002 Format(/1X,'Atom number ',I5,' not identical to that on file: ',
     1 I5)
 1003 Format(I3,2I1)
 1004 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,I7,6I3)
 1007 Format(A3,12I3,5I2,6I2,I2,I1,F7.5,6I3)
 1008 Format(A3,12I3,3I2,I2,I1,F7.5,I7,6I3)
 1009 Format(A3,12I3,3I2,I2,I1,F7.5,6I3)
 1010 FORMAT(/1X,'Cannot produce dual matrix, error IER= ',I2,
     1 ' Error in Database ==> Return')
 1011 FORMAT(1X,'Start at isomer ',I10)
 1012 FORMAT(1X,'End   at isomer ',I10)
 1013 FORMAT(1X,'File does not contain Hamilton cycle count,'
     1 ' statistics not performed')
 1014 FORMAT(1X,'Number of Hamilton cycles to small to do ',
     1 'statistics, width of bar would be ',I4)
 1015 Format(/1X,'Frequency of Hamilton cycles',
     1 /,1X,'midpoint gives number of Hamiltonian cycles ',
     1 ' at center of bar, range = midpoint plusminus width/2',
     1 /,1X,'bar',4x,'midpoint',2x,'hamcount',2X,'normalized',
     1 /1X,42('-'))
 1016 FORMAT(1X,I3,1X,I9,1X,I9,3X,E12.6)
 1017 FORMAT(1X,'Performing Hamilton cycle statistics with width ',I7)
 1018 FORMAT(/1X,'Number of bars =',I5,', width =',I9,', lowest',
     1 ' bar sits at Hamilton cycle count of ',I9)
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
 612  FORMAT(1X,'Lowest  Sigmah= ',F8.5,' for isomer ',I10,
     1      /1X,'Highest Sigmah= ',F8.5,' for isomer ',I10)
 613  FORMAT(1X,'Isomer List Complete')
 614  FORMAT(1X,'All isomers have ring spirals starting from a ',
     1 'pentagon')
 615  FORMAT(1X,'Number of isomers with ring spirals without ',
     1 'a pentagon start: ',I5)
 616  Format(1X,I3,' exceptions found for lower bound with isomer',
     1 ' numbers:')
 617  Format(1X,I3,' exceptions found for upper bound with isomer',
     1 ' numbers:')
 618  Format(10I9)
 619  Format(' Semi-tight lower and upper limits for Hamiltonian ',
     1 'cycles:'I9,'/',I9)
  701  FORMAT(1X,'General fullerene isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-'))
 702  FORMAT(1X,'General fullerene isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     2 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-'))
 703  FORMAT(1X,'Isolated-pentagon isomers of C',I2,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-'))
 704  FORMAT(1X,'Isolated-pentagon isomers of C',I3,':',
     1 ' (Np=0 implies IPR isomer, sigmah is the strain parameter, ',
     1 ' Ne the number of HOMO electrons, deg the HOMO degeneracy, ',
     1 /35x,' gap the HOMO-LUMO gap in units of beta, and NHamCyc the ',
     1 ' number of Hamiltonian cycles)',
     1 /8X,'n  PG   Ring spiral pentagon positions',
     3 19X,'Pentagon indices',5x,'Np  Hexagon indices',11x,'Sigmah',
     4 '   Ne  deg  gap    c/o     NHamCyc   NMR pattern',
     5 /1X,170('-'))

      Return
      END 
