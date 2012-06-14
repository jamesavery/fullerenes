      SUBROUTINE OptGraph(IOP,MAtom,Iout,IDA,IS,IC3,MDist,maxl,
     1 scalePPG,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C  This subroutine optimizes the fullerene graph using spring embedding
      DIMENSION Dist(2,NMAX),IC3(NMAX,3)
      DIMENSION IDA(NMAX,NMAX),IS(6),MDist(NMAX,NMAX)
      Data Rdist,ftol,dpi,conv/1.d0,.5d-10,3.14159265358979d0,
     1 6.0221367d-3/
      rmin=1.d10
      rmax=0.d0
      rper=0.d0
      maxd=0
      do i=1,MAtom
      do j=i+1,MAtom
       if(IDA(I,J).eq.1) then
        x=Dist(1,I)-Dist(1,J)
        y=Dist(2,I)-Dist(2,J)
        rd=dsqrt(x*x+y*y)
        if(rd.lt.rmin) rmin=rd
        if(rd.gt.rmax) rmax=rd
       endif
      enddo
        rv=dsqrt(Dist(1,I)**2+Dist(2,I)**2)
        if(rv.gt.rper) rper=rv
      enddo
      scale=Rdist/rmin
      If(IOP.eq.1) Write(IOUT,101)
      If(IOP.eq.2) Write(IOUT,102)
      If(IOP.eq.3) then
       do i=1,MAtom
       do j=i+1,MAtom
        if(Mdist(i,j).gt.maxd) maxd=Mdist(i,j)
       enddo
       enddo
       rper=rper*scale
       Write(IOUT,103) maxd,rper,scalePPG
       RAA=scalePPG
      endif
      if(IOP.eq.4) then
       RAA=rmax*scale/dfloat(maxl)
       Write(IOUT,104) maxl,RAA
      endif
      Write(IOUT,1000) rmin,Rdist
      do i=1,MAtom
       Dist(1,i)=Dist(1,i)*scale 
       Dist(2,i)=Dist(2,i)*scale 
       WRITE(IOUT,1001) I,Dist(1,I),Dist(2,I),(IC3(I,J),J=1,3)
      enddo
      CALL frprmng(IOP,MAtom,IDA,Iout,IS,MDist,
     1 maxd,Dist,ftol,iter,fret,E0,RAA)
      if(fret-E0.gt.1.d-2) then
       fretn=(fret-E0)/dfloat(MATOM)
       Write(IOUT,1002) fretn
      endif

  101 Format(/1X,'Optimization of fullerene graph using a ',
     1 'simple spring embedding algorithm for edges',
     1 /1X,'Energy at starting point set to zero')
  102 Format(/1X,'Optimization of fullerene graph using a ',
     1 'spring embedding algorithm for edges plus ',
     1 'Coulomb repulsion between barycenter and vertices',
     1 /1X,'Energy at starting point set to zero')
  103 Format(/1X,'Optimization of fullerene graph using a ',
     1 'scaled spring embedding algorithm for edges ',
     1 '(Pisanski-Plestenjak-Graovac algorithm)',
     1 /1X,'Max graph distance to periphery: ',I5,
     1 /1X,'Distance from barycenter to periphery: ',F15.4,
     1 /1X,'Scaling factor for exponential: ',F15.4,
     1 /1X,'Energy at starting point set to zero')
  104 Format(/1X,'Optimization of fullerene graph using a ',
     1 'Kamada-Kawai embedding algorithm',/1X,
     1 'Max integer graph distance: ',I5,', length of ',
     1 'display square area: ',F15.4,
     1 /1X,'Energy at starting point set to zero')
 1000 Format(/1X,'Fletcher-Reeves-Polak-Ribiere optimization',
     1 /1X,'Smallest distance in Tutte graph: ',F12.6,
     1 /1X,'Smallest distance reset to ',F12.6,
     1 /1X,'Rescaled Tutte graph coordinates:',
     1 /1X,'  Atom       X            Y        N1   N2   N3')
 1001 Format(1X,I4,2(1X,F12.6),1X,3(1X,I4))
 1002 Format(1X,'Energy gain per vertex: ',F12.6)

      Return 
      END

      SUBROUTINE frprmng(IOP,MATOM,AH,Iout,IS,MDist,
     1 maxd,p,ftol,iter,fret,E0,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,EPS=1.d-10)
      Real*8 p(NMAX*2),g(NMAX*2),h(NMAX*2),xi(NMAX*2)
      Real*8 pcom(NMAX*2),xicom(NMAX*2)
      Integer AH(NMAX,NMAX),IS(6),MDist(NMAX,NMAX), N
C     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
C     is performed on a function func, using its gradient as calculated by a routine dfunc.
C     The convergence tolerance on the function value is input as ftol. Returned quantities are
C     p (the location of the minimum), iter (the number of iterations that were performed),
C     and fret (the minimum value of the function). The routine linmin is called to perform
C     line minimizations. AH is the Hueckel adjacency matrix of atoms.
C     Parameters: NMAX is the maximum anticipated value of n; ITMAX is the maximum allowed
C     number of iterations; EPS is a small number to rectify special case of converging to exactly
C     zero function value.
C     USES dfuncg,funcg,linming
C     func input vector p of length n user defined to be optimized
C     IOP=1: spring embedding
C     IOP=2: spring + Coulomb embedding
C     IOP=3: Pisanski-Plestenjak-Graovac algorithm
C     IOP=4: Kamada-Kawai embedding
      N = 2*Matom
      iter=0
      CALL funcg(IOP,N,AH,IS,MDist,maxd,p,fp,RAA)
       E0=fp
      Write(Iout,1003) E0
C     dfunc input vector p of length N, output gradient of length n user defined
      CALL dfuncg(IOP,N,AH,IS,MDist,maxd,p,xi,RAA)
      grad2=0.d0
      do I=1,N
       grad2=grad2+xi(i)*xi(i)
      enddo
      grad=dsqrt(grad2)
      Write(Iout,1001) iter,fp-E0,grad
      if(grad.lt.ftol) return
      do j=1,N
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
      enddo
        fret=0.d0
      do its=1,ITMAX
        iter=its
        call linming(IOP,N,Iout,AH,IS,MDist,maxd,
     1       p,pcom,xi,xicom,fret,RAA)
         grad2=0.d0
         do I=1,n
          grad2=grad2+xi(i)*xi(i)
         enddo
         grad=dsqrt(grad2)
        Write(Iout,1001) iter,fret-E0,grad
        if(2.d0*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+EPS))then
          Write(Iout,1002) fret-E0,fret-fp
          return
        endif
        fp=fret
        CALL dfuncg(IOP,N,AH,IS,MDist,maxd,p,xi,RAA)
        gg=0.d0
        dgg=0.d0
        do j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
        enddo
        if(gg.eq.0.d0)return
        gam=dgg/gg
        do j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
        enddo   
      enddo
      Write(Iout,1000) fret,fret-fp
 1000 Format(' WARNING: Subroutine frprmn: maximum iterations exceeded',
     1 /1X,'energy ',F15.9,', diff= ',D12.3)
 1001 Format(' Iteration ',I4,', energy ',D14.8,', gradient ',D14.8)
 1002 Format(/1X,'Convergence achieved, energy ',F20.7,', diff= ',D12.3)
 1003 Format(/1X,'E0= ',D12.3)
      return
      END

      SUBROUTINE linming(IOP,n,Iout,AH,IS,MDist,
     1 maxd,p,pcom,xi,xicom,fret,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 p(NMAX*2),pcom(NMAX*2),xicom(NMAX*2),xi(NMAX*2)
      Integer AH(NMAX,NMAX),IS(6),MDist(NMAX,NMAX)
      PARAMETER (TOL=1.d-8)
C     USES brent,f1dim,mnbrak
      do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
      enddo
      ax=0.d0
      xx=1.d0
      CALL mnbrakg(IOP,n,Iout,AH,IS,MDist,maxd,
     1 ax,xx,bx,fa,fx,fb,xicom,pcom,RAA)
      CALL brentg(IOP,n,Iout,AH,IS,MDist,maxd,
     1 fret,ax,xx,bx,TOL,xmin,xicom,pcom,RAA)
      do j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
      enddo
      return
      END

      SUBROUTINE f1dimg(IOP,n,A,IS,MDist,maxd,
     1 f1dimf,x,xicom,pcom,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 pcom(NMAX*2),xt(NMAX*2),xicom(NMAX*2)
      Integer A(NMAX,NMAX),IS(6),MDist(NMAX,NMAX)
C     USES funcg
      do j=1,n
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      CALL funcg(IOP,n,A,IS,MDist,maxd,xt,f1dimf,RAA)
      return
      END

      SUBROUTINE mnbrakg(IOP,n,Iout,AH,IS,DD,maxd,
     1 ax,bx,cx,fa,fb,fc,xicom,pcom,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20)
      Integer AH(NMAX,NMAX),IS(6)
      Integer DD(NMAX,NMAX)
      REAL*8 pcom(NMAX*2),xicom(NMAX*2)
      CALL f1dimg(IOP,n,AH,IS,DD,maxd,fa,ax,xicom,pcom,
     1 RAA)
      CALL f1dimg(IOP,n,AH,IS,DD,maxd,fb,bx,xicom,pcom,
     1 RAA)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      CALL f1dimg(IOP,n,AH,IS,DD,maxd,fc,cx,xicom,pcom,
     1 RAA)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        else if((cx-u)*(u-ulim).gt.0.)then
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        else
          u=cx+GOLD*(cx-bx)
        if(u.gt.1.d10) then
        Write(Iout,1000)
        return
        endif
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
 1000 Format('**** Error in Subroutine mnbrakg')
      END

      SUBROUTINE brentg(IOP,n,Iout,AH,IS,DD,maxd,
     1 fx,ax,bx,cx,tol,xmin,xicom,pcom,RAA)
      use config
C BRENT is a FORTRAN library which contains algorithms for finding zeros 
C or minima of a scalar function of a scalar variable, by Richard Brent. 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,CGOLD=.3819660,ZEPS=1.d-10)
      REAL*8 pcom(NMAX*2),xicom(NMAX*2)
      Integer AH(NMAX,NMAX),IS(6)
      Integer DD(NMAX,NMAX)
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      CALL f1dimg(IOP,n,AH,IS,DD,maxd,fx,x,xicom,pcom,
     1 RAA)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.d0*tol1
        if(dabs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(dabs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).
     1      or.p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(dabs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        CALL f1dimg(IOP,n,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1    RAA)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      Write(Iout,1000)
 1000 Format(' WARNING: Subroutine brent: maximum iterations exceeded')
3     xmin=x
      return
      END

      SUBROUTINE OptFF(MAtom,Iout,IDA,N5,N6,
     1 N5MEM,N6MEM,Dist,Rdist,ftol,forceWu)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C  This subroutine optimizes the fullerene 3D structure using the Wu force field:
C   Z. C. Wu, D. A. Jelski, T. F. George, "Vibrational Motions of
C   Buckminsterfullerene", Chem. Phys. Lett. 137, 291-295 (1987).
C   Original force field in forceWu
C   Angstroem and rad is used for bond distances and bond length
C   Data from Table 1 of Wu in dyn/cm = 10**-3 N/m
      DIMENSION Dist(3,NMAX)
      DIMENSION IDA(NMAX,NMAX)
      DIMENSION N5MEM(MMAX,5),N6MEM(MMAX,6),force(9),forceWu(9)
      Data dpi,conv/3.14159265358979d0,6.0221367d-3/
      IOP=1
      force(1)=forceWu(1)
      force(2)=forceWu(2)
C     Conversion of angles in rad
      conrad=dpi/180.d0
      force(3)=forceWu(3)*conrad
      force(4)=forceWu(4)*conrad
C     Conversion of dyn/cm in a.u. / Angstroem**2
      force(5)=.5d0*forceWu(5)*conv*3.80879844d-4
      force(6)=.5d0*forceWu(6)*conv*3.80879844d-4
      force(7)=.5d0*forceWu(7)*conv*3.80879844d-4
      force(8)=.5d0*forceWu(8)*conv*3.80879844d-4
C     Leave parameter for Coulomb force as it is
      force(9)=forceWu(9)
      M=Matom/2+2
C     Optimize
      Write(IOUT,1000) Rdist
      Write(Iout,1003) (force(i),i=1,9),ftol
      if(forceWu(9).gt.0.d0) Write(Iout,1004) forceWu(9)
      CALL frprmn(IOP,MATOM*3,IDA,Iout,N5,N6,N5MEM,N6MEM,
     1 Dist,force,ftol,iter,fret)
      if(fret.gt.1.d-2) then
      fretn=fret/dfloat(MATOM)
      Write(IOUT,1002) fretn
      endif
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      fac=RDist/Rmin
      Do I=1,MATOM
      Dist(1,I)=Dist(1,I)*fac
      Dist(2,I)=Dist(2,I)*fac
      Dist(3,I)=Dist(3,I)*fac
      enddo
      CALL Distan(Matom,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(IOUT,1001) Rmin,Rmax,rms
 1000 Format(1X,'Optimization of geometry using harmonic oscillators',
     1 ' for stretching and bending modes using the force-field of ',
     1 ' Wu et al.',/1X,'Fletcher-Reeves-Polak-Ribiere algorithm used',
     1 /1X,'Smallest bond distance set to ',F12.6)
 1001 FORMAT(1X,'Minimum distance: ',F12.6,', Maximum distance: ',F12.6,
     1 ', RMS distance: ',F12.6)
 1002 FORMAT(1X,'Distances and angles defined in the force field can',
     1 ' not be reached',/1X,'Energy per atom in atomic units: ',F12.6)
 1003 Format(' Force field parameters: ',9F12.6,', Tolerance= ',D9.3,/)
 1004 Format(' Coulomb repulsion from center of origin with force ',
     1 F12.6,/)
     
      Return 
      END

      SUBROUTINE frprmn(IOP,N,AH,Iout,N5,N6,N5M,N6M,
     1 p,force,ftol,iter,fret)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=99999,EPS=1.d-9)
      Real*8 p(NMAX*3),g(NMAX*3),h(NMAX*3),xi(NMAX*3)
      Real*8 pcom(NMAX*3),xicom(NMAX*3),force(9)
      Integer AH(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6)
C     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
C     is performed on a function func, using its gradient as calculated by a routine dfunc.
C     The convergence tolerance on the function value is input as ftol. Returned quantities are
C     p (the location of the minimum), iter (the number of iterations that were performed),
C     and fret (the minimum value of the function). The routine linmin is called to perform
C     line minimizations. AH is the Hueckel adjacency matrix of atoms.
C     Parameters: NMAX is the maximum anticipated value of n; ITMAX is the maximum allowed
C     number of iterations; EPS is a small number to rectify special case of converging to exactly
C     zero function value.
C     USES dfunc,func,linmin
C     func input vector p of length n user defined to be optimized
C     IOP=1: Wu force field optimization
      iter=0
      CALL func(IOP,N,IERR,AH,N5,N6,N5M,N6M,p,fp,force)
      if(IERR.ne.0) then
      Write(Iout,1004)
      return
      endif
C     dfunc input vector p of length N, output gradient of length n user defined
      CALL dfunc(IOP,N,AH,N5,N6,N5M,N6M,p,xi,force)
      grad2=0.d0
      do I=1,N
       grad2=grad2+xi(i)*xi(i)
      enddo
      grad=dsqrt(grad2)
      Write(Iout,1001) iter,fp,grad
      if(grad.lt.ftol) return
      do j=1,N
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
      enddo
        fret=0.d0
      do its=1,ITMAX
        iter=its
        call linmin(IOP,N,AH,N5,N6,N5M,N6M,
     1   p,pcom,xi,xicom,fret,force)
         grad2=0.d0
         do I=1,n
          grad2=grad2+xi(i)*xi(i)
         enddo
         grad=dsqrt(grad2)
        Write(Iout,1001) iter,fret,grad
        if(2.d0*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+EPS))then
          Write(Iout,1002) fret,fret-fp
          return
        endif
        fp=fret
        CALL dfunc(IOP,N,AH,N5,N6,N5M,N6M,p,xi,force)
        gg=0.d0
        dgg=0.d0
        do j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
        enddo
        if(gg.eq.0.d0)return
        gam=dgg/gg
        do j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
        enddo   
      enddo
      Write(Iout,1000) fret,fret-fp
 1000 Format(' WARNING: Subroutine frprmn: maximum iterations exceeded',
     1 /1X,'energy ',F15.9,', diff= ',D12.3)
 1001 Format(' Iteration ',I6,', energy ',D14.8,', gradient ',D14.8)
 1002 Format(/1X,'Convergence achieved, energy ',D14.8,', diff= ',D12.3)
 1004 Format('**** Severe error in angle, check input coordiantes:',
     1 ' One angle either 0 or 180 degrees, ill-alligned structure',
     1 /1X,'Cannot optimize structure, check eigenvector input')
      return
      END

      SUBROUTINE linmin(IOP,n,AH,N5,N6,N5M,N6M,
     1 p,pcom,xi,xicom,fret,c)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 p(NMAX*3),pcom(NMAX*3),xicom(NMAX*3),xi(NMAX*3),c(9)
      Integer AH(NMAX,NMAX)
      Integer N5M(MMAX,5),N6M(MMAX,6)
      PARAMETER (TOL=1.d-5)
C     USES brent,f1dim,mnbrak
      do j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
      enddo
      ax=0.d0
      xx=1.d0
      CALL mnbrak(IOP,n,AH,N5,N6,N5M,N6M,
     1 ax,xx,bx,fa,fx,fb,xicom,pcom,c)
      CALL brent(IOP,n,AH,N5,N6,N5M,N6M,Iout,fret,
     1 ax,xx,bx,TOL,xmin,xicom,pcom,c)
      do j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
      enddo
      return
      END

      SUBROUTINE f1dim(IOP,n,A,N5,N6,N5M,N6M,
     1 f1dimf,x,xicom,pcom,c)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 pcom(NMAX*3),xt(NMAX*3),xicom(NMAX*3),c(9)
      Integer A(NMAX,NMAX)
      Integer N5M(MMAX,5),N6M(MMAX,6)
C     USES func
      do j=1,n
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      CALL func(IOP,n,IERR,A,N5,N6,N5M,N6M,xt,f1dimf,c)
      return
      END

      SUBROUTINE mnbrak(IOP,n,AH,N5,N6,N5M,N6M,
     1 ax,bx,cx,fa,fb,fc,xicom,pcom,c)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20)
      Integer AH(NMAX,NMAX)
      Integer N5M(MMAX,5),N6M(MMAX,6)
      REAL*8 pcom(NMAX*3),xicom(NMAX*3),c(9)
      CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1 fa,ax,xicom,pcom,c)
      CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1 fb,bx,xicom,pcom,c)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1 fc,cx,xicom,pcom,c)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
        else if((cx-u)*(u-ulim).gt.0.)then
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
        else
          u=cx+GOLD*(cx-bx)
        if(u.gt.1.d10) then
        Print*,'**** Error in Subroutine mnbrak'
        return
        endif
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END

      SUBROUTINE brent(IOP,n,AH,N5,N6,N5M,N6M,Iout,
     1 fx,ax,bx,cx,tol,xmin,xicom,pcom,c)
      use config
C BRENT is a FORTRAN library which contains algorithms for finding zeros 
C or minima of a scalar function of a scalar variable, by Richard Brent. 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,CGOLD=.3819660,ZEPS=1.d-10)
      REAL*8 pcom(NMAX*3),xicom(NMAX*3),c(9)
      Integer AH(NMAX,NMAX)
      Integer N5M(MMAX,5),N6M(MMAX,6)
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1 fx,x,xicom,pcom,c)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.d0*tol1
        if(dabs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(dabs(p).ge.abs(.5d0*q*etemp).or.p.le.q*(a-x).
     1      or.p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(dabs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        CALL f1dim(IOP,n,AH,N5,N6,N5M,N6M,
     1   fu,u,xicom,pcom,c)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      Write(Iout,1000)
 1000 Format(' WARNING: Subroutine brent: maximum iterations exceeded')
3     xmin=x
      return
      END

      SUBROUTINE powell(n,iter,Iout,IOP,ier,Matom,ftol,AN,RMDSI,
     1     p,pmax,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=20,TINY=1.D-20)
      REAL*8 p(n),pcom(n),xicom(n),xi(n,n),pt(n),ptt(n),xit(n),step(n)
      REAL*8 Dist(3,Nmax),pmax(n)
      
c numerical recipies cluster to do Powell minimization
c
c   p(n)      ... initial starting point (input); best point (output)
c   xi(n,n)   ... matrix containing the initial set of directions (input)
c   n         ... dimension; number of variables (input)
c   ftol      ... fractional tolerance in the function value (input)
c   iter      ... number of iterations taken
c   fret      ... value of f at p
      iter=0
      ier=0
      TOL=ftol
      If(IOP.eq.0) then
       call MDSnorm(n,Matom,fret,RMDSI,p,Dist)
      else
       call MAInorm(n,Matom,IP,AN,p,Dist)
       fret=-AN
      endif
      WRITE(IOUT,1002)
      WRITE(IOUT,1005) ftol,n
      do j=1,n
       step(j)=0.1d0
       pt(j)=p(j)
      enddo
      WRITE(IOUT,1006) (p(j),j=1,n),fret
        do i=1,n
        do j=1,i
          xi(j,i)=0.d0
          xi(i,j)=0.d0
          if(i.eq.j) xi(i,i)=step(i)
        enddo
        enddo
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.d0
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linminx(n,Matom,IOP,ier,TOL,p,xit,fret,pcom,xicom,Dist)
          if(ier.eq.1) Return
        if(dabs(fptt-fret).gt.del)then
          del=dabs(fptt-fret)
          ibig=i
        endif
13    continue
        WRITE(IOUT,1004) iter,(p(j),j=1,3),fret
      If(IOP.eq.0) then
       if(dabs(p(1)).gt.pmax(1)) ier=2
       if(dabs(p(2)).gt.pmax(2)) ier=2
       if(dabs(p(3)).gt.pmax(3)) ier=2
       if(ier.eq.2) return
      endif
      if(2.*dabs(fp-fret).le.ftol*(dabs(fp)+dabs(fret))+TINY) Return
      if(iter.eq.ITMAX) Go to 2
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      If(IOP.eq.0) then
       Call MDSnorm(n,Matom,fptt,RMDSI,ptt,Dist)
      else
       Call MAInorm(n,Matom,IP,AN,ptt,Dist)
       fptt=-AN
      endif
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linminx(n,Matom,IOP,ier,TOL,p,xit,fret,pcom,xicom,Dist)
      if(ier.eq.1) Return
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
 2    WRITE(IOUT,1007)
 1002 FORMAT(/1x,'Start Powell optimization')
 1004 FORMAT(' Iter: ',I6,' C(S): ',3(D14.8,1X)' Norm: ',D14.8)
 1005 Format(' Fractional tolerance ftol = ',D12.5,
     * ', dimension of problem n = ',i1)
 1006 FORMAT(' Start:       C(S): ',3(D14.8,1X),' Norm: ',D14.8)
 1007 Format(' WARNING: Optimizer Powell exceeding maximum iterations')
      Return
      END

      SUBROUTINE linminx(n,Matom,IOP,ier,
     1 TOL,p,xi,fret,pcom,xicom,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 p(n),xi(n),pcom(n),xicom(n)
      REAL*8 Dist(3,Nmax)
CU    USES brentx,f1dimx,mnbrakx
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call mnbrakx(n,Matom,IOP,ier,
     1 pcom,xicom,ax,xx,bx,fa,fx,fb,Dist)
      if(ier.eq.1) Return
      fret=brentx(n,Matom,IOP,pcom,xicom,ax,xx,bx,TOL,xmin,Dist)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END

      SUBROUTINE mnbrakx(ncom,Matom,IOP,ier,pcom,xicom,ax,bx,cx,
     1 fa,fb,fc,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20,HUGE=1.d10)
      REAL*8 pcom(ncom),xicom(ncom)
      REAL*8 Dist(3,Nmax)
      fa=f1dimx(ncom,Matom,IOP,ax,pcom,xicom,Dist)
      fb=f1dimx(ncom,Matom,IOP,bx,pcom,xicom,Dist)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=f1dimx(ncom,Matom,IOP,cx,pcom,xicom,Dist)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        if(ax.gt.HUGE.or.bx.gt.HUGE.or.CX.gt.HUGE) then
         ier=1
         Return
         endif
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
        else if((cx-u)*(u-ulim).gt.0.d0)then
          fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
        else
          u=cx+GOLD*(cx-bx)
          fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END

      DOUBLE PRECISION FUNCTION brentx(ncom,Matom,IOP,
     1 pcom,xicom,ax,bx,cx,tol,xmin,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Dist(3,Nmax)
      REAL*8 pcom(ncom),xicom(ncom)
      INTEGER ITMAX
      PARAMETER (ITMAX=1000,CGOLD=.3819660,ZEPS=1.D-10)
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f1dimx(ncom,Matom,IOP,x,pcom,xicom,Dist)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=dabs(q)
          etemp=e
          e=d
          if(dabs(p).ge.dabs(.5*q*etemp).or.p.le.q*(a-x).or.
     * p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(dabs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f1dimx(ncom,Matom,IOP,u,pcom,xicom,Dist)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      Print*, 'WARNING: brent exceed maximum iterations'
3     xmin=x
      brentx=fx
      return
      END
