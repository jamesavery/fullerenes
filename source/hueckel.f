      SUBROUTINE Hueckel(IOUT,IC3,nohueckel,IDA,A,evec,df)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C Perform Hueckel matrix diagonalization to obtain eigenvalues
C This gives a good hint if the fullerene is closed-shell
      DIMENSION IC3(Nmax,3),A(Nmax,Nmax),evec(Nmax),df(Nmax)
      DIMENSION IDA(Nmax,Nmax)

C Produce adjacency matrix
      Do I=1,number_vertices
      Do K=1,number_vertices
        A(I,K)=0.d0
        IDA(I,K)=0
      enddo
      enddo
      Do I=1,number_vertices
      Do J=1,3
        IP=IC3(I,J)
        A(I,IP)=1.d0
        A(IP,I)=1.d0
        IDA(I,IP)=1
        IDA(IP,I)=1
      enddo
      enddo

      if(nohueckel.ne.0) then
       WRITE(IOUT,1001) 
       return
      endif

C Diagonalize without producing eigenvectors
      WRITE(IOUT,1000) number_vertices,number_vertices 
      call tred2l(A,number_vertices,Nmax,evec,df)
      call tqlil(evec,df,number_vertices,Nmax)

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

C Analyze eigenenergies
      Call HueckelAnalyze(Iout,iocc,df,evec)
      
 1000 FORMAT(/1X,'Construct the (',I5,','I5,') Hueckel matrix')
 1001 FORMAT(/1X,'Skip diagonalization of Hueckel matrix')
      return
      END

      Subroutine HueckelAnalyze(Iout,iocc,df,evec)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION evec(Nmax),df(Nmax),IDG(Nmax),specmom(0:50)
      Character*10 Symbol
C Parameters alpha and beta are in atomic units and are adjusted 
C  to HOMO DFT orbital energies
      Data Tol,Tol1,alpha,beta,scale/1.d-5,.15d0,-.21d0,-0.111,3.1/
C Perform Hueckel matrix diagonalization to obtain eigenvalues
C Now sort degeneracies
      WRITE(IOUT,1007) 
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
       if(ntot.gt.number_vertices) then 
         if(nflag.eq.0) then
          nflag=1
          bandgap=df(i-1)-df(i)
          if(df(i).gt.-1.d-5) iproper=1
         endif
        NE=0
        Symbol='(empty)   '
       endif
       if(ntot.gt.number_vertices.and.
     1       (ntot-NE1).lt.number_vertices) then
        NE=number_vertices-ntot+NE1
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
      TRE=1.024296d0*Etot/dfloat(number_vertices)-1.562211d0
C     Subtracting the C60 value
      DTRE=TRE-2.82066353359331501d-2
      DTREkcal=DTRE*beta*6.27509541D+02
      Graphene=0.0468d0
      Write(Iout,1003) Etot,TRE,Graphene,DTRE,DTREkcal
      Write(Iout,1013) beta
      Write(Iout,1014) scale,DTREkcal*scale
C     Calculate Hueckel resonance energy
      vertnum=dfloat(number_vertices)
      Eres=Etot-vertnum
      EresPA=Eres/vertnum
      Write(Iout,1012) Eres,EresPA,Eres*beta,EresPA*beta

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

C     Spectral moments
      nspec=15
      vertnum=dfloat(number_vertices)
      specmom(0)=vertnum
       do J=1,nspec
        specmom(J)=0.d0
       enddo
      Do I=1,number_vertices
       do J=1,nspec
        specmom(J)=specmom(J)+evec(i)**J
       enddo
      enddo
      Write(Iout,1001) 
      specmom(0)=dfloat(number_vertices)
       do I=0,nspec
        Write(Iout,1011) I,dint(specmom(i)+1.d-5)
       enddo

      call flush(iout)
 
 1000 FORMAT(8X,'x',13X,'E',4X,'deg NE   type    ',/1X,45('-'))
 1001 FORMAT(/1X,'Spectral moments M(I):',/1X,26('-'))
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
 1008 FORMAT(1X,'Estrada index: ',F16.6,
     1 ' and bipartivity index: ',F14.8)
 1009 FORMAT(/1X,'Fullerene has a properly closed shell')
 1010 FORMAT(/1X,'Fullerene has a pseudo closed shell')
 1011 FORMAT(1X,I2,5X,F20.0)
 1012 FORMAT(1X,'Total resonance energy in units of beta: ',F12.6,
     1 ' (per atom: ',F12.6,')',
     2 /1X,'Total resonance energy in atomic units ',F12.6,
     3 ' (per atom: ',F12.6,')')
 1013 FORMAT(1X,'beta used in this program in atomic units: ',F8.3)
 1014 FORMAT(1X,'Scaling of ',F5.2,' applied to represent DFT value: ',
     1 /2X,'Difference of total resonance energy per atom ',
     1 'compared to C60 in kcal/mol: ',F12.6)
      Return
      END 
