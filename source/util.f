      SUBROUTINE PerfectMatching(MAtom,Iout,IDA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDA(Nmax,Nmax)
       Write(Iout,1000) MAtom
       Write(Iout,1001) 
       Write(Iout,1002) 
 1000 Format(/1X,'Upper limit for number of perfect matchings',
     1 ' in cubic graphs: 2**N with N=',I5)
 1001 Format(1X,'Counting the number of perfect matchings using',
     1 ' the Fisher-Kasteleyn-Temperley (FKT) algorithm')
 1002 Format(1X,'Not implemented yet')
      RETURN
      END

      SUBROUTINE Sortr(M,Mnew,imirror,jmirror,diam)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION imirror(Nmax),jmirror(Nmax),diam(Nmax)
      DIMENSION imirrorw(Nmax),jmirrorw(Nmax),diamw(Nmax)
      ICOUNT=0
      DO I=1,M
      DO K=I+1,M
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

      SUBROUTINE TopIndicators(Matom,Iout,IDA,Mdist)
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer MDist(Nmax,Nmax),wi(Nmax),wienermin,wienermax
      DIMENSION IDA(Nmax,Nmax)
      type(c_ptr) :: g, new_fullerene_graph
C     This routine calculates the Wiener index, Hyperwiener index,
C     minimal and maximal vertex contribution, rho and rhoE,
C     Schultz index and Balaban index
C     For details see D. Vukicevic,F. Cataldo, O. Ori, A. Graovac,
C     Chem. Phys. Lett. 501, 442–445 (2011).

      Write(Iout,1000)
      Do I=1,Nmax
      Do J=1,Nmax
       MDist(I,J)=0
      enddo
      enddo

C     Get topological distance matrix
      g = new_fullerene_graph(Nmax,Matom,IDA)
      call all_pairs_shortest_path(g,Matom,Nmax,MDist)

      nlimit=dint(dfloat(intmax8)/1.2d0)
      iwiener1=0
      iwiener=0
      ihyperwiener=0
      ierrorhw=0
      ierrorw=0
      ierrors=0
      maxdist=0
      Do I=1,MAtom
        wi(I)=0
       Do J=1,MAtom
        idist=MDist(I,J)
        wi(I)=wi(i)+idist
         if(J.gt.I) then
          if(idist.gt.maxdist) maxdist=idist
          if(ierrorhw.eq.0) then 
           ihyperwiener=ihyperwiener+(idist*(1+idist))/2
           if(ihyperwiener.gt.nlimit) then 
            ierrorhw=1
            ihyperwiener=0
           endif
          endif
         endif
       enddo
      if(ierrorw.eq.0) then
       iwiener1=iwiener1+wi(i)
       if(iwiener1.gt.nlimit) then
        ierrorw=1
        iwiener1=0
       endif
      endif
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
      Do I=1,MAtom
      Do J=I+1,MAtom
       if(IDA(I,J).eq.1) then
        wii=dfloat(wi(I))
        wij=dfloat(wi(J))
        if(wii.lt.1.d-15.or.wij.lt.1.d-15) then
         Write(Iout,1006)
         Return
        endif
        balaban=balaban+1.d0/(dsqrt(wii)*dsqrt(wij))
       endif
      enddo
      enddo
      vertnum=dfloat(MAtom)
      fac=3.d0*vertnum/(vertnum+4.d0)
      balabanindex=balaban*fac

      iwiener=iwiener1/2
      wav=dfloat(iwiener1)/vertnum
      rho=wav/dfloat(wienermin)
      rhoE=dfloat(wienermax)/dfloat(wienermin)
      isize=Matom*(Matom-1)
      Avdist=2.d0*dfloat(iwiener)/dfloat(isize)
      izagreb=MAtom*9
      if(iwiener.lt.intmax8/10) then
       ischultz=iwiener*6
      else
       ierrors=1
       ischultz=0
      endif
      wienerfac=dfloat(iwiener)/(9.d0*vertnum**3)
      Wienerbalaban=wienerfac*balabanindex*4.d0*(vertnum+4.d0)

      if(ierrorw.eq.1)  Write(Iout,1003)     
      if(ierrorhw.eq.1) Write(Iout,1004)     
      if(ierrors.eq.1)  Write(Iout,1005)     
      Write(Iout,1001) iwiener,ihyperwiener,wienermin,wienermax,
     1 wav,rho,rhoE,izagreb,ischultz,balabanindex
      Write(Iout,1007) Wienerbalaban
      Write(Iout,1002) maxdist,Avdist

 1000 Format(1X,'Topological Indicators:',/1X,
     1 'For definitions see Vukicevic et al., Chem. Phys. Lett. ',
     1 '501, 442 (2011), and Behtoei et al., Appl. Math. Lett. ',
     1 '22, 1571 (2009)')
 1001 Format(' Wiener index W: ',I12,/,' Hyper Wiener index WW: ',I12,
     1 /,' Minimal and maximal vertex contribution to W: ',I12,' and ',
     1 I12,', average vertex contribution (wav): ',D15.9,/,
     1 ' rho: ',D15.9,', rhoE: ',D15.9,/,' Zagreb index = nv*3^2 = ',
     1 I12,' (trivial for regular fullerenes)',/,
     1 ' Schultz index = 6*W = ',I12,' (related to Wiener index for ',
     1 'regular fullerenes)',/,' Balaban index = ',D15.9,
     1 /,' For the Estrada index see Subroutine Hueckel output')
 1002 Format(' Topological distances are between 1 and ',I6,/,
     1 ' Average topological distance: ',F12.6)
 1003 Format(' Integer overflow for Wiener index, W set to zero')
 1004 Format(' Integer overflow for Hyper Wiener index, WW set to zero')
 1005 Format(' Integer overflow for Schultz index,6*W set to zero')
 1006 Format(' Something wrong with Wiener sum')
 1007 Format(' f*Wiener*Balaban = 4WB(n+4)/(9n^3) = ',D15.9,
     1 ' (should be exactly 1.0 for cubic graphs)')

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

      SUBROUTINE Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
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
      Do I=1,Matom
      Do J=I+1,Matom
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

      SUBROUTINE Num2(MAtom,IArray,I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      I=IArray/MAtom
      J=IArray-I*MAtom
      If(J.eq.0) then
      J=MAtom
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
