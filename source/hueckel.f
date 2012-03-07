      SUBROUTINE Hueckel(NAtom,MAtom,IOUT,IC3,IDA,A,evec,df)
      IMPLICIT REAL*8 (A-H,O-Z)
C Perform Hueckel matrix diagonalization to obtain eigenvalues
C This gives a good hint if the fullerene is closed-shell
      DIMENSION IC3(natom,3),A(NAtom,NAtom),evec(Natom),df(Natom)
      DIMENSION IDG(natom),IDA(Natom,Natom)
      Character*10 Symbol
C Parameters alpha and beta are in atomic units and are adjusted 
C  to HOMO DFT orbital energies
      Data Tol,Tol1,alpha,beta/1.d-5,.15d0,-.21d0,-0.111/
C Produce adjacency matrix
      WRITE(IOUT,1001) Matom,Matom 
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

C Diagonalize without producing eigenvectors
      call tred2l(A,Matom,Natom,evec,df)
      call tqlil(evec,df,Matom,Natom)

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

C Now sort degeneracies
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
      Write(Iout,1000) 
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
       epsilon=alpha+df(i)*beta
       Write(Iout,1002) df(I),epsilon,idg(i),NE,Symbol
      enddo
      Write(Iout,1003)
      if(nopen.eq.1) then
       Write(Iout,1004)
      else
       excite=-bandgap*beta*27.2117
       Write(Iout,1005) bandgap,excite
       if(bandgap.lt.Tol1) Write(Iout,1006)
      endif
 
 1000 FORMAT(8X,'x',13X,'E',4X,'deg NE   type    ',/1X,45('-'))
 1001 FORMAT(/1X,'Construct the (',I3,','I3,') Hueckel matrix '
     1 'and diagonalize (E=alpha+x*beta; E in au)',
     1 /1X,'Eigenvalues are between [-3,+3] (in units of |beta|)',
     1 /1X,'deg: degeneracy; NE: number of electrons')
 1002 FORMAT(2(1X,F12.6),I3,1X,I3,3X,A10)
 1003 FORMAT(1X,45('-'))
 1004 FORMAT(1X,'Hueckel theory indicates that fullerene has',
     1 ' open-shell character (zero band gap)!')
 1005 FORMAT(1X,'Bandgap delta x =',F12.6,' in units of |beta| =',
     1 F12.6,' eV')
 1006 FORMAT(1X,'Caution: Bandgap small, possibility '
     1 'for open-shell character')
      Return
      End 
