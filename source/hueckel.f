      SUBROUTINE Hueckel(MAtom,IOUT,IC3,ihueckel,IDA,A,evec,df)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C Perform Hueckel matrix diagonalization to obtain eigenvalues
C This gives a good hint if the fullerene is closed-shell
      DIMENSION IC3(Nmax,3),A(Nmax,Nmax),evec(Nmax),df(Nmax)
      DIMENSION IDA(Nmax,Nmax)
C Produce adjacency matrix
      WRITE(IOUT,1000) Matom,Matom 
      Do I=1,MAtom
      Do K=1,MAtom
        A(I,K)=0.d0
        IDA(I,K)=0
      enddo
      enddo
      Do I=1,MAtom
      Do J=1,3
        IP=IC3(I,J)
        A(I,IP)=1.d0
        A(IP,I)=1.d0
        IDA(I,IP)=1
        IDA(IP,I)=1
      enddo
      enddo

      if(ihueckel.eq.0) then
       WRITE(IOUT,1001) 
       return
      endif
C Diagonalize without producing eigenvectors
      call tred2l(A,Matom,Nmax,evec,df)
      call tqlil(evec,df,Matom,Nmax)

C Sort eigenvalues
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
        endif
      enddo

C Analyze eigenenergies
      Call HueckelAnalyze(MAtom,NMax,Iout,iocc,df,evec)
      
 1000 FORMAT(/1X,'Construct the (',I5,','I5,') Hueckel matrix')
 1001 FORMAT(/1X,'Skip diagonalization of Hueckel matrix')
      return
      END

      Subroutine HueckelAnalyze(MAtom,NMax,Iout,iocc,df,evec)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION evec(Nmax),df(Nmax),IDG(Nmax)
      Character*10 Symbol
C Parameters alpha and beta are in atomic units and are adjusted 
C  to HOMO DFT orbital energies
      Data Tol,Tol1,alpha,beta/1.d-5,.15d0,-.21d0,-0.111/
C Perform Hueckel matrix diagonalization to obtain eigenvalues
C Now sort degeneracies
      WRITE(IOUT,1007) 
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
      Etot=0.d0
      iocc=0
      estrada=0.d0
      bipartivity=0.d0
      Write(Iout,1000)
      iproper=0 
      Do I=1,ieigv
       NE=2*idg(i)
       NE1=NE
       ntot=ntot+NE
       Symbol='(occupied)'
       if(ntot.gt.Matom) then 
         if(nflag.eq.0) then
          nflag=1
          bandgap=df(i-1)-df(i)
          if(df(i).gt.-1.d-5) iproper=1
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
       epsilon=alpha+df(i)*beta
       Etot=Etot+df(I)*dfloat(NE)
       aidg=dfloat(idg(i))
       bipartivity=bipartivity+aidg*dcosh(df(i))
       estrada=estrada+aidg*dexp(df(i))
       Write(Iout,1002) df(I),epsilon,idg(i),NE,Symbol
      enddo

C     Other useful properties from Hueckel matrix
C     Babic's resonance energy
      TRE=1.024296d0*Etot/dfloat(MAtom)-1.562211d0
      DTRE=TRE-2.82066353359331501d-2
      DTREkcal=DTRE*beta*6.27509541D+02
      Graphene=0.0468d0
      Write(Iout,1003) Etot,TRE,Graphene,DTRE,DTREkcal

C     Estrada and bibartivity index
      Write(Iout,1008) estrada,bipartivity/estrada

C     Bandgap
      if(nopen.eq.1) then
       Write(Iout,1004)
      else
       if(iproper.eq.0) then
         Write(Iout,1009)
        else
         Write(Iout,1010)
       endif
       excite=-bandgap*beta*27.2117
       Write(Iout,1005) bandgap,excite
       if(bandgap.lt.Tol1) Write(Iout,1006)
      endif

      call flush(iout)
 
 1000 FORMAT(8X,'x',13X,'E',4X,'deg NE   type    ',/1X,45('-'))
 1002 FORMAT(2(1X,F12.6),I3,1X,I3,3X,A10)
 1003 FORMAT(1X,45('-'),/1X,'Total pi-energy in units of beta: ',F12.6,
     1 /1X,'Total resonance energy per atom in units of beta ',
     1     'according to Babic: ',F12.6,' (graphene limit: ',
     1     F12.6,')',
     1 /1X,'Difference of total resonance energy per atom ',
     1     'compared to C60 in units of beta: ',F12.6,
     1 /1X,'Difference of total resonance energy per atom ',
     1     'compared to C60 in kcal/mol: ',F12.6)
 1004 FORMAT(1X,'Hueckel theory indicates that fullerene has',
     1 ' open-shell character (zero band gap)!')
 1005 FORMAT(1X,'Bandgap delta x =',F12.6,' in units of |beta| =',
     1 F12.6,' eV')
 1006 FORMAT(1X,'Caution: Bandgap small, possibility '
     1 'for open-shell character')
 1007 FORMAT(/1X,'Diagonalize Hueckel matrix (E=alpha+x*beta; E in au)',
     1 /1X,'Eigenvalues are between [-3,+3] (in units of |beta|)',
     1 /1X,'deg: degeneracy; NE: number of electrons')
 1008 FORMAT(1X,'Estrada index: ',F12.6,
     1 ' and bipartivity index: ',F12.6)
 1009 FORMAT(/1X,'Fullerene has a properly closed shell')
 1010 FORMAT(/1X,'Fullerene has a pseudo closed shell')
      Return
      END 
