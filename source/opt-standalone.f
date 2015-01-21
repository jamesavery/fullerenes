      SUBROUTINE SA_OptGraph(N, IOP,Iout,IDA,IS,MDist,
     1                       maxl,scalePPG,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C  This subroutine optimizes the fullerene graph using spring embedding
      DIMENSION Dist(2,N)
      DIMENSION IDA(N,N),IS(6),MDist(N,N)
      real*8 Rdist,ftol
      rdist=1.d0
      ftol=0.5d-10
      rmin=1.d10
      rmax=0.d0
      rper=0.d0
      maxd=0
      do i=1,N
        do j=i+1,N
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
      If(IOP.eq.1) Write(*,101)
      If(IOP.eq.2) Write(*,102)
      If(IOP.eq.3) then
        do i=1,N
        do j=i+1,N
          if(Mdist(i,j).gt.maxd) maxd=Mdist(i,j)
        enddo
        enddo
        rper=rper*scale
        Write(*,103) maxd,rper,scalePPG
        RAA=scalePPG
      endif
      if(IOP.eq.4) then
        RAA=rmax*scale/dfloat(maxl)
        Write(*,104) maxl,RAA
      endif
      Write(*,1000) rmin,Rdist
      do i=1,N
        Dist(1,i)=Dist(1,i)*scale
        Dist(2,i)=Dist(2,i)*scale
        WRITE(*,1001) I,Dist(1,I),Dist(2,I)
      enddo
      CALL SA_frprmn2d(N,IOP,IDA,Iout,IS,MDist,
     1 maxd,Dist,ftol,iter,fret,E0,RAA)
      if(fret-E0.gt.1.d-2) then
        fretn=(fret-E0)/dfloat(N)
        Write(*,1002) fretn
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
 1001 Format(1X,I4,2(1X,F12.6),1X)
 1002 Format(1X,'Energy gain per vertex: ',F12.6)

      Return
      END

      SUBROUTINE SA_frprmn2d(N,IOP,AH,Iout,IS,MDist,
     1 maxd,p,ftol,iter,fret,E0,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,EPS=1.d-10)
      Real*8 p(N*2),g(N*2),h(N*2),xi(N*2)
      Real*8 pcom(N*2),xicom(N*2)
      Integer AH(N,N),IS(6),MDist(N,N)
C     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
C     is performed on a function func3d, using its gradient as calculated by a routine dfunc2d.
C     The convergence tolerance on the function value is input as ftol. Returned quantities are
C     p (the location of the minimum), iter (the number of iterations that were performed),
C     and fret (the minimum value of the function). The routine linmin3d is called to perform
C     line minimizations. AH is the Hueckel adjacency matrix of atoms.
C     Parameters: N is the maximum anticipated value of n; ITMAX is the maximum allowed
C     number of iterations; EPS is a small number to rectify special case of converging to exactly
C     zero function value.
C     USES dfunc2d,func2d,linmin2d
C     func3d input vector p of length n user defined to be optimized
C     IOP=1: spring embedding
C     IOP=2: spring + Coulomb embedding
C     IOP=3: Pisanski-Plestenjak-Graovac algorithm
C     IOP=4: Kamada-Kawai embedding

      iter=0
      CALL SA_func2d(N,IOP,AH,IS,MDist,maxd,p,fp,RAA)
       E0=fp
      Write(Iout,1003) E0
C     dfunc2d input vector p of length 2*N, output gradient of length 2*N user defined
      CALL SA_dfunc2d(N,IOP,AH,IS,MDist,maxd,p,xi,RAA)
      grad2=0.d0
      do I=1,2*N
       grad2=grad2+xi(i)*xi(i)
      enddo
      grad=dsqrt(grad2)
      Write(Iout,1001) iter,fp-E0,grad
      if(grad.lt.ftol) return
      do j=1,2*N
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
      enddo
        fret=0.d0
      do its=1,ITMAX
        iter=its
        CALL SA_linmin2d(N,IOP,Iout,AH,IS,MDist,maxd,
     1       p,pcom,xi,xicom,fret,RAA)
         grad2=0.d0
         do I=1,2*N
          grad2=grad2+xi(i)*xi(i)
         enddo
         grad=dsqrt(grad2)
        Write(Iout,1001) iter,fret-E0,grad
        if(2.d0*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+EPS))then
          Write(Iout,1002) fret-E0,fret-fp
          return
        endif
        fp=fret
        CALL SA_dfunc2d(N,IOP,AH,IS,MDist,maxd,p,xi,RAA)
        gg=0.d0
        dgg=0.d0
        do j=1,2*N
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
        enddo
        if(gg.eq.0.d0)return
        gam=dgg/gg
        do j=1,2*N
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
        enddo
      enddo
      Write(Iout,1000) fret,fret-fp
 1000 Format(' WARNING: Subroutine frprmn2d: maximum iterations
     1 exceeded',/1X,'energy ',F15.9,', diff= ',D12.3)
 1001 Format(' Iteration ',I4,', energy ',D14.8,', gradient ',D14.8)
 1002 Format(/1X,'Convergence achieved, energy ',F20.7,', diff= ',D12.3)
 1003 Format(/1X,'E0= ',D12.3)
      return
      END

      SUBROUTINE SA_linmin2d(N,IOP,Iout,AH,IS,MDist,
     1 maxd,p,pcom,xi,xicom,fret,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 p(N*2),pcom(N*2),xicom(N*2),xi(N*2)
      Integer AH(N,N),IS(6),MDist(N,N)
      PARAMETER (TOL=1.d-8)
C     USES brent2d,f1dim2d,mnbrak2d
      do j=1,2*N
        pcom(j)=p(j)
        xicom(j)=xi(j)
      enddo
      ax=0.d0
      xx=1.d0
      CALL SA_mnbrak2d(N,IOP,Iout,AH,IS,MDist,maxd,
     1 ax,xx,bx,fa,fx,fb,xicom,pcom,RAA)
      CALL SA_brent2d(N,IOP,Iout,AH,IS,MDist,maxd,
     1 fret,ax,xx,bx,TOL,xmin,xicom,pcom,RAA)
      do j=1,2*N
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
      enddo
      return
      END

      SUBROUTINE SA_f1dim2d(N,IOP,A,IS,MDist,maxd,
     1 f1dimf,x,xicom,pcom,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 pcom(N*2),xt(N*2),xicom(N*2)
      Integer A(N,N),IS(6),MDist(N,N)
C     USES func2d
      do j=1,2*N
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      CALL SA_func2d(N,IOP,A,IS,MDist,maxd,xt,f1dimf,RAA)
      return
      END

      SUBROUTINE SA_mnbrak2d(N,IOP,Iout,AH,IS,DD,maxd,
     1 ax,bx,cx,fa,fb,fc,xicom,pcom,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20)
      Integer AH(N,N),IS(6)
      Integer DD(N,N)
      REAL*8 pcom(N*2),xicom(N*2)
      CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fa,ax,xicom,pcom,
     1 RAA)
      CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fb,bx,xicom,pcom,
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
      CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fc,cx,xicom,pcom,
     1 RAA)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
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
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        else if((cx-u)*(u-ulim).gt.0.)then
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        else
          u=cx+GOLD*(cx-bx)
        if(u.gt.1.d10) then
        Write(Iout,1000)
        return
        endif
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
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
 1000 Format('**** Error in Subroutine mnbrak2d')
      END

      SUBROUTINE SA_brent2d(N,IOP,Iout,AH,IS,DD,maxd,
     1 fx,ax,bx,cx,tol,xmin,xicom,pcom,RAA)
      use config
C BRENT is a FORTRAN library which contains algorithms for finding zeros
C or minima of a scalar function of a scalar variable, by Richard Brent.
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,CGOLD=.3819660,ZEPS=1.d-10)
      REAL*8 pcom(N*2),xicom(N*2)
      Integer AH(N,N),IS(6)
      Integer DD(N,N)
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      d=0.d0
      CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fx,x,xicom,pcom,
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
        CALL SA_f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
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
 1000 Format('WARNING: Subroutine brent2d: maximum iterations exceeded')
3     xmin=x
      return
      END



      SUBROUTINE SA_func2d(N,IOP,A,IS,MDist,maxd,p,fc,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Embedding algorithms for fullerene graph, energy
      Real*8 p(N*2)
      Integer A(N,N),IS(6),MDist(N,N)
      Data r,f,coulomb/2.0d0,1.d-1,.3d0/
      fc=0.d0
C     simple spring embedding
      if(IOP.le.2) then
      ehook=0.d0
      Do I=1,2*N,2
        I1=(I+1)/2
      Do J=I+2,2*N,2
        J1=(J+1)/2
        if(A(I1,J1).ne.0) then
         px=p(I)-p(J)
         py=p(I+1)-p(J+1)
         ratom=dsqrt(px*px+py*py)
          ehook=ehook+(ratom-r)**2
         endif
      enddo
      enddo
C     total energy
      fc=f*ehook
      endif

C     Add Repulsive Coulomb
      if(IOP.eq.2) then
      fcoul=0.d0
      Do I=1,2*N,2
        rv=dsqrt(p(I)**2+p(I+1)**2)
        fcoul=fcoul+1.d0/rv
      enddo
        fc=fc+coulomb*fcoul
      endif

C     Pisanski-Plestenjak-Graovac algorithm
      if(IOP.eq.3) then
      rpermax=dfloat(maxd)
      af=RAA*(1.4d0+rpermax/24.d0)
        iloop=6
        if(IS(6).eq.0) iloop=5
      ehook=0.d0
      Do I=1,2*N,2
        I1=(I+1)/2
        maxu=1000000
        do K=1,iloop
         K1=IS(K)
         if(Mdist(K1,I1).lt.maxu) maxu=Mdist(K1,I1)
        enddo
      Do J=I+2,2*N,2
        J1=(J+1)/2
        maxp=1000000
        do K=1,iloop
         K2=IS(K)
         if(Mdist(K2,J1).lt.maxp) maxp=Mdist(K2,J1)
        enddo
        if(A(I1,J1).ne.0) then
         px=p(I)-p(J)
         py=p(I+1)-p(J+1)
         ratom2=px*px+py*py
         dup=dfloat(2*maxd-maxu-maxp)
         weight=dexp(af*dup/rpermax)
         ehook=ehook+ratom2*weight
        endif
      enddo
      enddo
C     total energy
      fc=ehook
      endif

C     Kamada-Kawai embedding
      if(IOP.eq.4) then
      ehook=0.d0
      Do I=1,2*N,2
        I1=(I+1)/2
      Do J=I+2,2*N,2
        J1=(J+1)/2
        if(MDist(I1,J1).le.0) then
         Print*,I1,J1,MDist(I1,J1)
         stop
        endif
        DD=dfloat(MDist(I1,J1))
         px=p(I)-p(J)
         py=p(I+1)-p(J+1)
         ratom=dsqrt(px*px+py*py)
          ehook=ehook+((ratom-RAA*DD)/DD)**2
      enddo
      enddo
C     total energy
      fc=f*ehook
      endif

      Return
      END

      SUBROUTINE SA_dfunc2d(N,IOP,A,IS,MDist,maxd,p,x,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Embedding algorithms for fullerene graph, gradient
      Real*8 p(N*2),x(N*2)
      Integer A(N,N),IS(6),MDist(N,N)
      Data r,f,coulomb/2.0d0,1.d-1,.3d0/
C     simple spring embedding
      if(IOP.le.2) then
      Do I=1,2*N,2
        ehookx=0.d0
        ehooky=0.d0
        I1=(I+1)/2
      Do J=1,2*N,2
        J1=(J+1)/2
        if(A(I1,J1).ne.0) then
          px=p(I)-p(J)
          py=p(I+1)-p(J+1)
          ratom=dsqrt(px*px+py*py)
            fac=f*(ratom-r)/ratom
            ehookx=ehookx+fac*px
            ehooky=ehooky+fac*py
        endif
      enddo
C     Fix outer vertices
        if(I1.eq.IS(1).or.I1.eq.IS(2).or.I1.eq.IS(3).
     1   or.I1.eq.IS(4).or.I1.eq.IS(5).or.I1.eq.IS(6)) then
         x(I)=0.d0
         x(I+1)=0.d0
        else
         x(I)  =2.d0*ehookx
         x(I+1)=2.d0*ehooky
        endif
      enddo
      endif

C     Repulsive Coulomb
      if(IOP.eq.2) then
      Do I=1,2*N,2
        I1=(I+1)/2
        if(I1.ne.IS(1).and.I1.ne.IS(2).and.I1.ne.IS(3).
     1   and.I1.ne.IS(4).and.I1.ne.IS(5).and.I1.ne.IS(6)) then
        xyf=(p(i)**2+p(i+1)**2)**(-1.5d0)
          x(I)=x(I)-p(I)*xyf*coulomb
          x(I+1)=x(I+1)-p(I+1)*xyf*coulomb
        endif
      enddo
      endif

C     Pisanski-Plestenjak-Graovac algorithm
      if(IOP.eq.3) then
C     simple spring embedding + exponential weighting
      rpermax=dfloat(maxd)
      af=RAA*(1.4d0+rpermax/24.d0)
        iloop=6
        if(IS(6).eq.0) iloop=5
      Do I=1,2*N,2
        I1=(I+1)/2
        maxu=1000000
        do K=1,iloop
         K1=IS(K)
         if(Mdist(K1,I1).lt.maxu) maxu=Mdist(K1,I1)
        enddo
C     Fix outer vertices
       if(I1.eq.IS(1).or.I1.eq.IS(2).or.I1.eq.IS(3).
     1  or.I1.eq.IS(4).or.I1.eq.IS(5).or.I1.eq.IS(6)) then
        x(I)=0.d0
        x(I+1)=0.d0
       else
        ehookx=0.d0
        ehooky=0.d0
       Do J=1,2*N,2
        J1=(J+1)/2
        if(A(I1,J1).ne.0) then
        maxp=1000000
        do K=1,iloop
         K2=IS(K)
         if(Mdist(K2,J1).lt.maxp) maxp=Mdist(K2,J1)
        enddo
          px=p(I)-p(J)
          py=p(I+1)-p(J+1)
          dup=dfloat(2*maxd-maxu-maxp)
          weight=dexp(af*dup/rpermax)
            ehookx=ehookx+px*weight
            ehooky=ehooky+py*weight
        endif
       enddo
        x(I)  =ehookx
        x(I+1)=ehooky
       endif
      enddo
      endif

C     Kamada-Kawai embedding
      if(IOP.eq.4) then
      Do I=1,2*N,2
        ehookx=0.d0
        ehooky=0.d0
        I1=(I+1)/2
      Do J=1,2*N,2
        if(I.ne.J) then
        J1=(J+1)/2
        DD=dfloat(MDist(I1,J1))
          px=p(I)-p(J)
          py=p(I+1)-p(J+1)
          ratom=dsqrt(px*px+py*py)
            fac=(ratom-RAA*DD)/(ratom*DD*DD)
            ehookx=ehookx+fac*px
            ehooky=ehooky+fac*py
        endif
      enddo
         x(I)  =f*ehookx
         x(I+1)=f*ehooky
      enddo
      endif

      return
      END

C---------------------------------------------
C     OPTIMIZE 3D GEOMETRY
C---------------------------------------------
C
C Input:
C     * graph: C++ graph object
C     * int N: Number of vertices
C     * bool ihessian: If 1, compute hessian matrix
C     * bool iprinthessian: If 1, print hessian matrix
C     * int iopt: Optimization method
C     * real ftol: Convergence tolerance
C     * ? force
C Input/output:
C     * real*8 Dist(3,N): On input, initial coordinate vector; on output, final coordinates.

      SUBROUTINE SA_OptFF(graph,N,ihessian,iprinthessian,iopt,
     1  Dist,ftol,force)
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
C  This subroutine optimizes the fullerene 3D structure using a force field
c  (e.g. the Wu force field):  Z. C. Wu, D. A. Jelski, T. F. George, "Vibrational
c  Motions of Buckminsterfullerene", Chem. Phys. Lett. 137, 291-295 (1987).
C  Angstroem and rad is used for bond distances and bond length
C  Data from Table 1 of Wu in dyn/cm = 10**-3 N/m
      DIMENSION Dist(3,N)
      DIMENSION IDA(N,N)
      real(8) force(ffmaxdim)
      integer iopt,ideg(N*3)
      real(8) hessian(N*3,N*3),
     1 evec(N*3),df(N*3)
      type(c_ptr) :: graph

c edges with 0, 1, 2 pentagons
      integer e_hh(2,3*N/2), e_hp(2,3*N/2), e_pp(2,3*N/2)
      integer a_h(3,3*N-60), a_p(3,60)
      integer d_hhh(4,N), d_hhp(4,N), d_hpp(4,N), d_ppp(4,N)
c counter for edges with 0, 1, 2 pentagons neighbours
      integer ne_hh,ne_hp,ne_pp
      integer nd_hhh,nd_hhp,nd_hpp,nd_ppp
      integer N
      number_vertices=N

      call adjacency_matrix(graph,N,IDA)
      call SA_get_edges(graph,N,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp)
      call SA_get_corners(graph,N,a_h,a_p)
      if(iopt .eq. 3 .or. iopt.eq.4) then
        call SA_get_dihedrals(graph,N,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      endif

      select case(iopt)
      case(1)
        Write(*,1000)
        Write(*,1016) ftol,(force(i),i=1,8)
      case(2)
        Write(*,1000)
        Write(*,1019) ftol,(force(i),i=1,9)
      case(3)
        Write(*,1007)
        Write(*,1020) ftol,(force(i),i=1,18)
      case(4)
        Write(*,1007)
        Write(*,1018) ftol,(force(i),i=1,19)
      end select
C       Conversion to kJ/mol
C       energies in KJ/mol, gradients in kJ/mol/A and hessian kJ/mol/A^2
        unitconv=1.d-20 * 6.02214129d20
      if(iopt.eq.1 .or. iopt.eq.2)then
c        force(1)=force(1)
c        force(2)=force(2)
C       Conversion of angles in rad
        force(3)=force(3)*deg2rad
        force(4)=force(4)*deg2rad
        force(5)=force(5)*unitconv
        force(6)=force(6)*unitconv
        force(7)=force(7)*unitconv
        force(8)=force(8)*unitconv
C       Leave parameter for Coulomb force as it is
c        force(9)=force(9)
      else if (iopt.eq.3 .or. iopt.eq.4) then
c        force(1)=force(1)
c        force(2)=force(2)
c        force(3)=force(3)
C       Conversion of angles and dihedrals in rad
        force(4)=force(4)*deg2rad
        force(5)=force(5)*deg2rad
        force(6)=force(6)*deg2rad
        force(7)=force(7)*deg2rad
        force(8)=force(8)*deg2rad
        force(9)=force(9)*deg2rad
C       Conversion of dyn/cm in a.u. / Angstroem**2
        force(10)=force(10)*unitconv
        force(11)=force(11)*unitconv
        force(12)=force(12)*unitconv
        force(13)=force(13)*unitconv
        force(14)=force(14)*unitconv
        force(15)=force(15)*unitconv
        force(16)=force(16)*unitconv
        force(17)=force(17)*unitconv
        force(18)=force(18)*unitconv
C       Leave parameter for Coulomb force as it is
c        force(19)=force(19)
      end if
      select case(iopt)
      case(1)
        Write(*,1006) (force(i),i=1,8)
      case(2)
        Write(*,1003) (force(i),i=1,9)
      case(3)
        Write(*,1005) (force(i),i=1,18)
      case(4)
        Write(*,1008) (force(i),i=1,19)
      end select
      if(iopt.eq.2 .and. force(9).gt.0.d0) Write(*,1004) force(9)


C OPTIMIZE
      CALL SA_frprmn3d(N,
     1 Dist,force,iopt,ftol,iter,fret,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      if(fret.gt.1.d-2) then
        fretn=fret/dfloat(N)
        Write(*,1002) fretn
      endif
      CALL SA_Distan(N,IDA,Dist,Rmin,Rminall,Rmax,rms)
      Write(*,1001) Rmin,Rmax,rms


C HESSIAN
      if(ihessian.ne.0) then
        call SA_get_hessian(N,dist, force, iopt, hessian,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        if(iprinthessian.gt.1) then
          write(*,1023)
          write(*,1024)
     1      ((hessian(i,j),i=1,3*N),j=1,3*N)
        endif
C Diagonalize without producing eigenvectors
C  Mass of 12-C used
        amassC=12.0d0
        fachess=3.80879844d-4/amassC
        convw=2720.21d0
C       Test if Hessian is symmetric
        symmetric=0.d0
        test=1.d-10
        do i=1,3*N
          do j=1,3*N
            symmetric=symmetric+dabs(hessian(i,j)-hessian(j,i))
          enddo
        enddo
        asym=symmetric*.5d0
        if(asym.gt.test) then
          Write(*,1013) asym
        else
          Write(*,1015) asym
        endif
C       Mass-weight Hessian
        do i=1,3*N
          do j=1,3*N
            hessian(i,j)=hessian(i,j)*fachess
          enddo
        enddo
        call tred2l(hessian,3*N,3*N,evec,df)
        call tqlil(evec,df,3*N,3*N)
C Sort eigenvalues
        negeig=0
        Do I=1,N*3
          e0=evec(I)
          jmax=I
          Do J=I+1,N*3
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
        if(iprinthessian.ne.0) then
          write(*,1009)
          write(*,1010) (evec(i),i=1,3*N)
        endif
        Do I=1,N*3
          if(evec(i).lt.0.d0) then
            negeig=negeig+1
            evec(i)=-dsqrt(-evec(i))
          else
            evec(i)=dsqrt(evec(i))
          endif
        enddo
        write(*,1011) negeig
        Do I=1,N*3
          evec(i)=evec(i)*convw
        enddo
        if(iprinthessian.ne.0) then
          write(*,1012)
          write(*,1010) (evec(i),i=1,N*3)
        endif
C Zero-point vibrational energy
        zerops=0.d0
        Do I=1,N*3-6
          zerops=zerops+evec(i)
        enddo
        zerop=zerops*.5d0
        zeropwn=zerop
        zeropau=zerop/au2wavenumbers
        zeropeV=zeropau*au2eV
        peratom=dfloat(N)
        write(*,1014) zeropau,zeropeV,zeropwn
        write(*,1026) zeropau/peratom,zeropeV/peratom,zeropwn/peratom
C Sort for degeneracies
        tolfreq=1.d-1
        icount=0
        idegc=0
        ndeg=0
        Do I=1,N*3-6
          idegc=idegc+1
          dif=evec(i)-evec(i+1)
          if(dif.gt.tolfreq) then
            icount=icount+1
            evec(icount)=evec(i)
            ideg(icount)=idegc
            ndeg=ndeg+idegc
            idegc=0
          endif
        enddo
        write(*,1021) ndeg,N*3-6
        write(*,1022) (evec(i),ideg(i),i=1,icount)
        write(*,1025)
     1    (evec(i),i=N*3-5,3*N)
      endif

 1000 Format(1X,'Optimization of geometry using harmonic oscillators',
     1 ' for stretching and bending modes using the force-field of',
     1 ' Wu et al.',/1X,'Fletcher-Reeves-Polak-Ribiere algorithm used')
 1001 FORMAT(1X,'Minimum distance: ',F12.6,', Maximum distance: ',F12.6,
     1 ', RMS distance: ',F12.6)
 1002 FORMAT(1X,'Distances and angles defined in the force field can',
     1 ' not be reached',/1X,'Energy per atom in atomic units: ',F12.6)
 1003 Format(' Force field parameters in au/A^2 and au/rad^2:',
     1 /1X,9F12.6,/)
 1004 Format(' Coulomb repulsion from center of origin with force ',
     1 F12.6,/)
 1005 Format(' Force field parameters in au/A^2 and au/rad^2:',
     1 /1X,18F12.6,/)
 1006 Format(' Force field parameters in au/A^2 and au/rad^2:',
     1 /1X,8F12.6,/)
 1007 Format(1X,'Optimization of geometry using harmonic oscillators',
     1 ' for bonds, angles, and dihedrals using an extension of the',
     1 ' force-field of Wu et al.',/1X,'Fletcher-Reeves-Polak-Ribiere',
     1 ' algorithm used')
 1008 Format(' Force field parameters in au/A^2 and au/rad^2:',
     1 /1X,19F12.6,/)
 1009 Format(' Eigenvalues of mass-weighted Hessian:')
 1010 Format(10(1X,D12.6))
 1011 Format(' Number of zero and negative eigenvalues: ',I6)
 1012 Format(' Frequencies (cm-1):')
 1013 Format(' Severe problem. Hessian is not symmetric: asym= ',d12.6)
 1014 Format(' Zero-point vibrational energy: ',d12.6,' a.u. , ',
     1 d12.6,' eV , ',d12.6,' cm-1 , ')
 1015 Format(' Hessian is symmetric: asym= ',d12.6)
 1016 Format(' Tolerance= ',D9.3,', Force field parameters in ',
     1 'A, deg, N/m:',/1X,8F12.3)
 1018 Format(' Tolerance= ',D9.3,', Force field parameters in ',
     1 'A, deg, N/m:'/1X,19F12.2)
 1019 Format(' Tolerance= ',D9.3,', Force field parameters in ',
     1 'A, deg, N/m:'/1X,9F12.2)
 1020 Format(' Tolerance= ',D9.3,', Force field parameters in ',
     1 'A, deg, N/m:'/1X,18F12.2)
 1021 Format(1X,I6,' non-zero frequencies (should be ',I6,').',
     1 ' Frequencies (in cm-1) and (quasi) degeneracies (n):')
 1022 Format(10(' ',f7.1,'(',I2,')'))
 1023 Format(' Hessian matrix:')
 1024 Format(8(d12.6,' '))
 1025 Format(' Zero frequencies for translation and rotation: ',
     1 6(d12.6,' '))
 1026 Format(' Zero-point vibrational energy per atom: ',d12.6,
     1 ' a.u. , ',d12.6,' eV , ',d12.6,' cm-1 , ')

      Return
      END

      SUBROUTINE SA_Distan(N,IDA,Dist,Rmin,Rminall,Rmax,rms)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,N)
      DIMENSION IDA(N,N)
C     Calculate minimal, maximal and root mean square distances
C     from adjacancy matrix IDA and cartesian coordinates Dist
      Rmin=1.d10
      Rminall=1.d10
      Rmax=0.d0
      Rrms=0.d0
      mc=0
      Do I=1,N
         Do J=I+1,N
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

      SUBROUTINE SA_frprmn3d(N,
     1 p,force,iopt,ftol,iter,fret,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=99999,EPS=1.d-9)
      Real*8 p(N*3),g(N*3),h(N*3),xi(N*3)
      Real*8 pcom(N*3),xicom(N*3)
      real*8 force(ffmaxdim)
      integer iopt
      integer e_hh(2,3*number_vertices/2), e_hp(2,3*number_vertices/2),
     1  e_pp(2,3*number_vertices/2)
      integer a_h(3,3*number_vertices-60), a_p(3,60)
      integer d_hhh(4,number_vertices), d_hhp(4,number_vertices),
     1  d_hpp(4,number_vertices), d_ppp(4,number_vertices)

C     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
C     is performed on a function func3d, using its gradient as calculated by a routine dfunc3d.
C     The convergence tolerance on the function value is input as ftol. Returned quantities are
C     p (the location of the minimum), iter (the number of iterations that were performed),
C     and fret (the minimum value of the function). The routine linmin3d is called to perform
C     line minimizations. AH is the Hueckel adjacency matrix of atoms.
C     Parameters: NMAX is the maximum anticipated value of n; ITMAX is the maximum allowed
C     number of iterations; EPS is a small number to rectify special case of converging to exactly
C     zero function value.
C     USES dfunc3d,func3d,linmin3d
C     func3d input vector p of length n user defined to be optimized
C     IOPT=1: Wu force field optimization
      iter=0
      CALL func3d(IERR,p,fp,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      if(IERR.ne.0) then
        Write(*,1004)
        return
      endif
C     dfunc3d input vector p of length N, output gradient of length n user defined
      CALL dfunc3d(p,xi,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      grad2=0.d0
      do I=1,3*N
        grad2=grad2+xi(i)*xi(i)
      enddo
      grad=dsqrt(grad2)
      Write(*,1001) iter,fp,grad
      if(grad.lt.ftol) return
      do j=1,3*N
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
      enddo
      fret=0.d0
      do its=1,ITMAX
c       turn off coulomb pot towards the end (and go to iopt=3 to indicate that coulomb has been shut of)
        if(iopt.eq.4 .and. force(19).gt.0.d0 .and. grad.le.1.d1) then
          force(19)=0.0d0
          iopt=3
          write(*,*)'Switching off coulomb repulsive potential.'
        endif
        iter=its
        call SA_linmin3d(N,p,pcom,xi,xicom,fret,
     1    force,iopt,e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1    a_h,a_p,
     1    d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        grad2=0.d0
        do I=1,3*N
          grad2=grad2+xi(i)*xi(i)
        enddo
        grad=dsqrt(grad2)
c        if(damping.eq.0) then
          write(*,1001) iter,fret,grad
c        else
c          write(*,1002) iter,fret,grad,damping
c        endif
        if(2.d0*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+EPS))then
          fretperatom=3.d0*fret/dfloat(3*N)
          Write(*,1003) fret,fret-fp,fretperatom
          return
        endif
        fp=fret
        CALL dfunc3d(p,xi,force,iopt,
     1    e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1    a_h,a_p,
     1    d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        gg=0.d0
        dgg=0.d0
        do j=1,3*N
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
        enddo
        if(gg.eq.0.d0)return
        gam=dgg/gg
        do j=1,3*N
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
        enddo
      enddo
      Write(*,1000) fret,fret-fp
 1000 Format(' WARNING: Subroutine frprmn3d: maximum iterations
     1 exceeded',/1X,'energy ',F15.9,', diff= ',D12.3)
 1001 Format(' Iteration ',I6,', energy [kJ/mol] ',D14.8,
     1 ', gradient [kJ/mol/A] ',D14.8)
c 1002 Format(' Iteration ',I6,', energy ',D14.8,', gradient ',D14.8,
c     1 ' The displacements of ',I4,' atoms were damped.')
 1003 Format(/1X,'Convergence achieved, energy [kJ/mol] ',D14.8,
     1 ', diff= ',D12.3,/1X,'Energy per atom [kJ/mol]: ',D14.8)
 1004 Format('**** Severe error in angle, check input coordiantes:',
     1 ' One angle either 0 or 180 degrees, ill-alligned structure',
     1 /1X,'Cannot optimize structure, check eigenvector input')
      return
      END

      SUBROUTINE SA_linmin3d(N,p,pcom,xi,xicom,fret,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)!,damping)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 p(N*3),pcom(N*3),xicom(N*3),xi(N*3)
      real*8 force(ffmaxdim)
      PARAMETER (TOL=1.d-5)
      integer e_hh(2,3*number_vertices/2), e_hp(2,3*number_vertices/2),
     1  e_pp(2,3*number_vertices/2)
      integer a_h(3,3*number_vertices-60), a_p(3,60)
      integer d_hhh(4,number_vertices), d_hhp(4,number_vertices),
     1  d_hpp(4,number_vertices), d_ppp(4,number_vertices)
c      real*8 length, cutoff, xi_tmp(nmax*3)
c      integer damping

C     USES brent3d,f1dim3d,mnbrak3d
      do j=1,3*N
        pcom(j)=p(j)
        xicom(j)=xi(j)
      enddo
      ax=0.d0
      xx=1.d0
      CALL SA_mnbrak3d(N,
     1 ax,xx,bx,fa,fx,fb,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      CALL SA_brent3d(N,fret,
     1 ax,xx,bx,TOL,xmin,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
c lets scale all displacements that are longer than a chosen cutoff to that cutoff.
c the direction of the displacement vector is maintained
      do j=1,3*N
        xi(j)=xmin*xi(j)
c        xi_tmp(j)=xi(j)
      enddo
c larger cutoffs result in faster convergence and are less safe ...
c      cutoff=3.0d1
c     count the number of atoms/displacements that were damped
c      damping=0
c      do j=1,n,3
c        length=dsqrt(xi(j)*xi(j) + xi(j+1)*xi(j+1) + xi(j+2)*xi(j+2))
c        if (length .gt. cutoff) then
c          xi_tmp(j)  =xi_tmp(j)  *(cutoff/length)
c          xi_tmp(j+1)=xi_tmp(j+1)*(cutoff/length)
c          xi_tmp(j+2)=xi_tmp(j+2)*(cutoff/length)
c          damping=damping + 1
c        endif
c      enddo
      do j=1,3*N
c        p(j)=p(j)+xi_tmp(j)
        p(j)=p(j)+xi(j)
      enddo
      return
      END


      SUBROUTINE SA_f1dim3d(N,
     1 f1dimf,x,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 pcom(N*3),xt(N*3),xicom(N*3),force(ffmaxdim)
      integer e_hh(2,3*N/2), e_hp(2,3*N/2), e_pp(2,3*N/2)
      integer a_h(3,3*N-60), a_p(3,60)
      integer d_hhh(4,N), d_hhp(4,N), d_hpp(4,N), d_ppp(4,N)

C     USES func3d
      do j=1,3*N
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      CALL func3d(IERR,xt,f1dimf,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      return
      END

      SUBROUTINE SA_mnbrak3d(N,
     1 ax,bx,cx,fa,fb,fc,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20)
      REAL*8 pcom(N*3),xicom(N*3),force(ffmaxdim)
      integer e_hh(2,3*N/2), e_hp(2,3*N/2), e_pp(2,3*N/2)
      integer a_h(3,3*N-60), a_p(3,60)
      integer d_hhh(4,N), d_hhp(4,N), d_hpp(4,N), d_ppp(4,N)

      CALL SA_f1dim3d(N,
     1 fa,ax,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      CALL SA_f1dim3d(N,
     1 fb,bx,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      CALL SA_f1dim3d(N,
     1 fc,cx,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
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
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        else if((cx-u)*(u-ulim).gt.0.)then
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     2   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        else
          u=cx+GOLD*(cx-bx)
        if(u.gt.1.d10) then
        Print*,'**** Error in Subroutine mnbrak3d'
        return
        endif
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
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

      SUBROUTINE SA_brent3d(N,
     1 fx,ax,bx,cx,tol,xmin,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
C BRENT is a FORTRAN library which contains algorithms for finding zeros
C or minima of a scalar function of a scalar variable, by Richard Brent.
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,CGOLD=.3819660,ZEPS=1.d-10)
      REAL*8 pcom(N*3),xicom(N*3),force(ffmaxdim)
      integer e_hh(2,3*N/2), e_hp(2,3*N/2), e_pp(2,3*N/2)
      integer a_h(3,3*N-60), a_p(3,60)
      integer d_hhh(4,N), d_hhp(4,N), d_hpp(4,N), d_ppp(4,N)

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      d=0.d0
      CALL SA_f1dim3d(N,
     1 fx,x,xicom,pcom,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
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
        CALL SA_f1dim3d(N,
     1   fu,u,xicom,pcom,force,iopt,
     1   e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1   a_h,a_p,
     1   d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
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
      Write(*,1000)
 1000 Format(' WARNING: Subroutine brent3d: maximum iterations
     1 exceeded')
3     xmin=x
      return
      END

      SUBROUTINE SA_linminx(n,IOP,ier,
     1 TOL,p,xi,fret,pcom,xicom,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 p(n),xi(n),pcom(n),xicom(n)
      REAL*8 Dist(3,N)
CU    USES brentx,f1dimx,mnbrakx
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call SA_mnbrakx(n,IOP,ier,
     1 pcom,xicom,ax,xx,bx,fa,fx,fb,Dist)
      if(ier.eq.1) Return
      fret=sa_brentx(n,IOP,ier,pcom,xicom,ax,xx,bx,TOL,xmin,Dist)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END

      SUBROUTINE SA_mnbrakx(N,IOP,ier,pcom,xicom,ax,bx,cx,
     1 fa,fb,fc,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20,HUGE=1.d10)
      REAL*8 pcom(N),xicom(N)
      REAL*8 Dist(3,N)
      fa=sa_f1dimx(N,IOP,ier,ax,pcom,xicom,Dist)
      fb=sa_f1dimx(N,IOP,ier,bx,pcom,xicom,Dist)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=sa_f1dimx(N,IOP,ier,cx,pcom,xicom,Dist)
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
          fu=sa_f1dimx(N,IOP,ier,u,pcom,xicom,Dist)
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
          fu=sa_f1dimx(N,IOP,ier,u,pcom,xicom,Dist)
        else if((cx-u)*(u-ulim).gt.0.d0)then
          fu=sa_f1dimx(N,IOP,ier,u,pcom,xicom,Dist)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=sa_f1dimx(N,IOP,ier,u,pcom,xicom,Dist)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          fu=sa_f1dimx(N,IOP,ier,u,pcom,xicom,Dist)
        else
          u=cx+GOLD*(cx-bx)
          fu=sa_f1dimx(ncom,IOP,ier,u,pcom,xicom,Dist)
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

      DOUBLE PRECISION FUNCTION SA_brentx(N,IOP,ier,
     1 pcom,xicom,ax,bx,cx,tol,xmin,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Dist(3,N)
      REAL*8 pcom(N),xicom(N)
      INTEGER ITMAX
      PARAMETER (ITMAX=1000,CGOLD=.3819660,ZEPS=1.D-10)
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      d=0.d0
      fx=sa_f1dimx(N,IOP,ier,x,pcom,xicom,Dist)
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
        fu=sa_f1dimx(N,IOP,ier,u,pcom,xicom,Dist)
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
      sa_brentx=fx
      return
      END


C     Lukas: The following three subroutines iterate over
C      1. edges
C      2. corners
C      3. dihedrals
C     in linear time (i.e. constant time per element), while providing the information you asked for.
C
C     Assumptions:
C     - g is a fullerene graph and has had a call to set_layout2d(g,layout2d) with its Tutte embedding as layout.
C       (or another strictly planar layout)
C     - N = Nvertices(g)
C     ------------------------------------------------------------
C                              EDGES
C     ------------------------------------------------------------
      SUBROUTINE SA_get_edges(graph,N,
     1 edges_hh,edges_hp,edges_pp,na_hh,na_hp,na_pp)
      use iso_c_binding
      type(c_ptr) :: graph
      integer edges(2,3*N/2), faceA(6), faceB(6), NE, np, i, u, v, lA,lB
      integer edges_hh(2,3*N/2), edges_hp(2,3*N/2), edges_pp(2,3*N/2)
c     counter for edges with 0, 1, 2 pentagons neighbours
      integer na_hh,na_hp,na_pp
      na_hh=0
      na_hp=0
      na_pp=0

      do j=1,2
      do i=1,3*N/2
        edges_pp(j,i)=0
        edges_hp(j,i)=0
        edges_hh(j,i)=0
      enddo
      enddo
      call edge_list(graph,edges,NE)

      do i=1,NE
C       Edge u--v
        u = edges(1,i)
        v = edges(2,i)

C       Edge is part of how many pentagons?
        call get_arc_face(graph,u,v,faceA,lA) ! O(1) operation
        call get_arc_face(graph,v,u,faceB,lB) ! O(1) operation
        np = 12-lA-lB

C       Do what needs to be done to u--v here
        select case(np)
        case(0)
          na_hh=na_hh+1
          edges_hh(1,na_hh)=u+1
          edges_hh(2,na_hh)=v+1
        case(1)
          na_hp=na_hp+1
          edges_hp(1,na_hp)=u+1
          edges_hp(2,na_hp)=v+1
        case(2)
          na_pp=na_pp+1
          edges_pp(1,na_pp)=u+1
          edges_pp(2,na_pp)=v+1
        case default
          write(*,*)'Something went horribly wrong: bond not adjacent ',
     1 'to 0, 1 or 2 pentagons'
          exit
        end select

c        write (*,*) "Edge ",(/u,v/)," connects ",np,"pentagons.",lA,lB
      end do

      END SUBROUTINE

C     ------------------------------------------------------------
C                              CORNERS
C     ------------------------------------------------------------
      SUBROUTINE SA_get_corners(graph,N,a_h,a_p)
c     here, n is the number of atoms
      use iso_c_binding
      type(c_ptr) :: graph
      integer pentagons(5,12), hexagons(6,N/2-10), u,v,w,i,j
c     arrays for atoms that are part of angles ...
      integer a_p(3,60), a_h(3,3*n-60)
c     counter for angles around hexagons and pentagons
      integer NH

      NH = N/2-10

      call compute_fullerene_faces(graph,pentagons,hexagons) ! Yields faces with vertices ordered CCW. O(N)
C     Every directed edge u->v is part of two "angles", corresponding to the neighbours to v that aren't u,
C     corresponding to both a CW and a CCW traversal starting in the edge. (Likewise, every edge is part of four)
C     The angles are each counted exactly once if we trace the outline of each face in e.g. CCW order.
      do i=1,12
C     iterate over angles u--v--w
c         write (*,*) "Pentagon number",i,"has corners:"
         do j=1,5
            u = pentagons(j,i)
            v = pentagons(MOD(j,5)+1,i)
            w = pentagons(MOD(j+1,5)+1,i)
C     Do what needs to be done to u--v--w here. Each of these are part of a pentagon, obviously.

            a_p(1,5*(i-1)+j)=u+1
            a_p(2,5*(i-1)+j)=v+1
            a_p(3,5*(i-1)+j)=w+1

c            write (*,*) j,":",(/u,v,w/)
         end do
      end do

      do i=1,NH
c         write (*,*) "Hexagon number",i,"has corners:"
C     iterate over angles u--v--w
         do j=1,6
            u = hexagons(j,i)
            v = hexagons(MOD(j,6)+1,i)
            w = hexagons(MOD(j+1,6)+1,i)
C     Do what needs to be done to u--v--w here. Each of these are part of a hexagon, obviously.

            a_h(1,6*(i-1)+j)=u+1
            a_h(2,6*(i-1)+j)=v+1
            a_h(3,6*(i-1)+j)=w+1

c            write (*,*) j,":",(/u,v,w/)
         end do
      end do
      END SUBROUTINE


C     ------------------------------------------------------------
C                              DIHEDRALS
C     ------------------------------------------------------------
      SUBROUTINE SA_get_dihedrals(graph,N,
     1 d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use iso_c_binding
      integer neighbours(3,N), face(6), lA,lB,lC, u,s,r,t
      type(c_ptr) :: graph
c     arrays for dihedrals. one per atom, starting in the middle
      integer d_hhh(4,n),d_hhp(4,n),d_hpp(4,n),d_ppp(4,n)
c     counter for dihedrals with 0, 1, 2, 3 pentagons neighbours
      integer nd_hhh,nd_hhp,nd_hpp,nd_ppp
      nd_hhh=0
      nd_hhp=0
      nd_hpp=0
      nd_ppp=0

      do j=1,4
      do i=1,n
        d_hhh(j,i)=0
        d_hhp(j,i)=0
        d_hpp(j,i)=0
        d_ppp(j,i)=0
      enddo
      enddo

      call adjacency_list(graph,3,neighbours)

      do u=1,N
C neighbours should be ordered CCW
C      t   B   s
C        \   /
C      C   u   A
C          |
C          r
         r = neighbours(1,u)
         s = neighbours(2,u)
         t = neighbours(3,u)

c         write (*,*) "Dihedral ",u-1,r-1,s-1,t-1
         call get_face(graph,s,u,r,6,face,lA)
         call get_face(graph,t,u,s,6,face,lB)
         call get_face(graph,r,u,t,6,face,lC)

         select case ( lA+lB+lC )
         case ( 15 )            ! (5,5,5) - all pentagons
c            write (*,*) "555"
C     Do stuff here

          nd_ppp=nd_ppp+1
          d_ppp(1,nd_ppp)=u
          d_ppp(2,nd_ppp)=r
          d_ppp(3,nd_ppp)=s
          d_ppp(4,nd_ppp)=t

         case ( 16 )            ! Two pentagons, one hexagon
C     Do stuff common to all three (2,1)-cases here

C     Do case specific stuff here
            select case ( lA*100+lB*10+lC )
            case ( 655 )  ! BC are pentagons, u--t common edge
c               write (*,*) "655"
              nd_hpp=nd_hpp+1
              d_hpp(1,nd_hpp)=u
              d_hpp(2,nd_hpp)=t
              d_hpp(3,nd_hpp)=r
              d_hpp(4,nd_hpp)=s

            case ( 565 )  ! AC are pentagons, u--r common edge
c               write (*,*) "565"
              nd_hpp=nd_hpp+1
              d_hpp(1,nd_hpp)=u
              d_hpp(2,nd_hpp)=r
              d_hpp(3,nd_hpp)=s
              d_hpp(4,nd_hpp)=t

            case ( 556 )  ! AB are pentagons, u--s common edge
c               write (*,*) "556"
              nd_hpp=nd_hpp+1
              d_hpp(1,nd_hpp)=u
              d_hpp(2,nd_hpp)=s
              d_hpp(3,nd_hpp)=t
              d_hpp(4,nd_hpp)=r

            end select

         case ( 17 )            ! One pentagon, two hexagons
C     Do stuff common to all three (1,2)-cases here

C     Do case specific stuff here
            select case ( lA*100+lB*10+lC )
            case ( 566 )  ! BC are hexagons, u--t common edge
c               write (*,*) "566"
              nd_hhp=nd_hhp+1
              d_hhp(1,nd_hhp)=u
              d_hhp(2,nd_hhp)=t
              d_hhp(3,nd_hhp)=r
              d_hhp(4,nd_hhp)=s

            case ( 656 )  ! AC are hexagons, u--r common edge
c               write (*,*) "656"
              nd_hhp=nd_hhp+1
              d_hhp(1,nd_hhp)=u
              d_hhp(2,nd_hhp)=r
              d_hhp(3,nd_hhp)=s
              d_hhp(4,nd_hhp)=t

            case ( 665 )  ! AB are hexagons, u--s common edge
c               write (*,*) "665"
              nd_hhp=nd_hhp+1
              d_hhp(1,nd_hhp)=u
              d_hhp(2,nd_hhp)=s
              d_hhp(3,nd_hhp)=t
              d_hhp(4,nd_hhp)=r

            end select

         case ( 18 )            ! (6,6,6) - all hexagons
C     Do stuff here
c            write (*,*) "666"

          nd_hhh=nd_hhh+1
          d_hhh(1,nd_hhh)=u
          d_hhh(2,nd_hhh)=r
          d_hhh(3,nd_hhh)=s
          d_hhh(4,nd_hhh)=t

         case DEFAULT
            write (*,*) "INVALID: ",(/lA,lB,lC/)
         end select

      end do

      END SUBROUTINE


      SUBROUTINE SA_get_hessian(N,coord, force, iopt, hessian,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p,
     1  d_hhh,d_hhp,d_hpp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
c      use iso_c_binding
      use config
c      type(c_ptr) :: graph
      implicit real*8 (a-h,o-z)
      integer iopt, i, j, m
      integer e_hh(2,3*N/2), e_hp(2,3*N/2), e_pp(2,3*N/2)
      integer ne_hh, ne_hp, ne_pp
      integer a_h(3,3*N-60), a_p(3,60)
      integer d_hhh(4,N), d_hhp(4,N), d_hpp(4,N), d_ppp(4,N)
      integer nd_hhh, nd_hhp, nd_hpp, nd_ppp
      integer a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
      real(8) coord(N*3), force(ffmaxdim),
     1 hessian(3*N,3*N)
      real(8) k

c initialize variables to *something*
      ah=0.0
      ap=0.0
      dhhh=0.0
      dhhp=0.0
      dhpp=0.0
      dppp=0.0
      fah=0.0
      fap=0.0
      fco=0.0
      fdhhh=0.0
      fdhhp=0.0
      fdhpp=0.0
      fdppp=0.0
      frh=0.0
      frhh=0.0
      frhp=0.0
      frp=0.0
      frpp=0.0
      rh=0.0
      rhh=0.0
      rhp=0.0
      rp=0.0
      rpp=0.0

c init
      do i=1,3*N
        do j=1,3*N
          hessian(i,j)=0.0
        enddo
      enddo

c get force constants
      select case(iopt)
        case(1, 2)
          rp=force(1)
          rh=force(2)
          ap=force(3)
          ah=force(4)
          frp=force(5)
          frh=force(6)
          fap=force(7)
          fah=force(8)
          fco=force(9)
        case(3, 4)
          rpp=force(1)
          rhp=force(2)
          rhh=force(3)
          ap=force(4)
          ah=force(5)
          dppp=force(6)
          dhpp=force(7)
          dhhp=force(8)
          dhhh=force(9)
          frpp=force(10)
          frhp=force(11)
          frhh=force(12)
          fap=force(13)
          fah=force(14)
          fdppp=force(15)
          fdhpp=force(16)
          fdhhp=force(17)
          fdhhh=force(18)
          fco=force(19)
        case default
          write(*,*)'Something went horribly wrong: illegal iopt'
      end select

c edges
      edge_types: do i=1,3
        select case(10*iopt + i)
          case(11,21)
            k=frh
            r_naught=rh
            m=ne_hh
          case(12,22)
            k=frp
            r_naught=rp
            m=ne_hp
          case(13,23)
            k=frp
            r_naught=rp
            m=ne_pp
          case(31,41)
            k=frhh
            r_naught=rhh
            m=ne_hh
          case(32,42)
            k=frhp
            r_naught=rhp
            m=ne_hp
          case(33,43)
            k=frpp
            r_naught=rpp
            m=ne_pp
          case default
            write(*,*)'Something went horribly wrong: illegal iopt',
     1 ' or edge type'
            exit
        end select
        if(m.gt.0) then
          edges: do j=1,m
            select case(i)
              case(1)
                a1=3*e_hh(1,j)-2
                a2=3*e_hh(1,j)-1
                a3=3*e_hh(1,j)
                a4=3*e_hh(2,j)-2
                a5=3*e_hh(2,j)-1
                a6=3*e_hh(2,j)
              case(2)
                a1=3*e_hp(1,j)-2
                a2=3*e_hp(1,j)-1
                a3=3*e_hp(1,j)
                a4=3*e_hp(2,j)-2
                a5=3*e_hp(2,j)-1
                a6=3*e_hp(2,j)
              case(3)
                a1=3*e_pp(1,j)-2
                a2=3*e_pp(1,j)-1
                a3=3*e_pp(1,j)
                a4=3*e_pp(2,j)-2
                a5=3*e_pp(2,j)-1
                a6=3*e_pp(2,j)
              case default
                write(*,*)'Something went horribly wrong'
                exit
            end select
            ax=coord(a1)
            ay=coord(a2)
            az=coord(a3)
            bx=coord(a4)
            by=coord(a5)
            bz=coord(a6)
            call dddist(ax, ay, az, bx, by, bz,
     2       dax, day, daz, dbx, dby, dbz,
     3       daxax, daxay, daxaz, daxbx, daxby, daxbz, dayay, dayaz,
     1       daybx, dayby, daybz, dazaz, dazbx, dazby, dazbz, dbxbx,
     1       dbxby, dbxbz, dbyby, dbybz, dbzbz, dist)
c           (partial^2 E/partial x_i partial x_j)
c              = (partial^2 E/partial x_i partial r)(partial r/partial x_j) + (partial^2 r/partial x_i partial x_j)(partial E/partial r)
c              = (partial (k(r - r_0))/partial x_i)(partial r/partial x_j) + (partial^2 r/partial x_i partial x_j)(k(r - r_0))
c              = k * ((partial r/partial x_i)(partial r/partial x_j) + (partial^2 r/partial x_i partial x_j)(r - r_0))
            diff=dist - r_naught
            hessian(a1,a1)=hessian(a1,a1) + k * (dax*dax + daxax*diff)
            hessian(a1,a2)=hessian(a1,a2) + k * (dax*day + daxay*diff)
            hessian(a1,a3)=hessian(a1,a3) + k * (dax*daz + daxaz*diff)
            hessian(a1,a4)=hessian(a1,a4) + k * (dax*dbx + daxbx*diff)
            hessian(a1,a5)=hessian(a1,a5) + k * (dax*dby + daxby*diff)
            hessian(a1,a6)=hessian(a1,a6) + k * (dax*dbz + daxbz*diff)
            hessian(a2,a2)=hessian(a2,a2) + k * (day*day + dayay*diff)
            hessian(a2,a3)=hessian(a2,a3) + k * (day*daz + dayaz*diff)
            hessian(a2,a4)=hessian(a2,a4) + k * (day*dbx + daybx*diff)
            hessian(a2,a5)=hessian(a2,a5) + k * (day*dby + dayby*diff)
            hessian(a2,a6)=hessian(a2,a6) + k * (day*dbz + daybz*diff)
            hessian(a3,a3)=hessian(a3,a3) + k * (daz*daz + dazaz*diff)
            hessian(a3,a4)=hessian(a3,a4) + k * (daz*dbx + dazbx*diff)
            hessian(a3,a5)=hessian(a3,a5) + k * (daz*dby + dazby*diff)
            hessian(a3,a6)=hessian(a3,a6) + k * (daz*dbz + dazbz*diff)
            hessian(a4,a4)=hessian(a4,a4) + k * (dbx*dbx + dbxbx*diff)
            hessian(a4,a5)=hessian(a4,a5) + k * (dbx*dby + dbxby*diff)
            hessian(a4,a6)=hessian(a4,a6) + k * (dbx*dbz + dbxbz*diff)
            hessian(a5,a5)=hessian(a5,a5) + k * (dby*dby + dbyby*diff)
            hessian(a5,a6)=hessian(a5,a6) + k * (dby*dbz + dbybz*diff)
            hessian(a6,a6)=hessian(a6,a6) + k * (dbz*dbz + dbzbz*diff)
          end do edges
        endif
      end do edge_types

c angles
      angle_types: do i=1,2
c       iopt doesn't matter in this case
        select case(i)
          case(1)
            k=fah
            a_naught=ah
            m=3*N-60
          case(2)
            k=fap
            a_naught=ap
            m=60
          case default
            write(*,*)'Something went horribly wrong'
            exit
        end select
        if(m.gt.0)then
          angles: do j=1,m
            select case(i)
              case(1)
                a1=3*a_h(1,j)-2
                a2=3*a_h(1,j)-1
                a3=3*a_h(1,j)
                a4=3*a_h(2,j)-2
                a5=3*a_h(2,j)-1
                a6=3*a_h(2,j)
                a7=3*a_h(3,j)-2
                a8=3*a_h(3,j)-1
                a9=3*a_h(3,j)
              case(2)
                a1=3*a_p(1,j)-2
                a2=3*a_p(1,j)-1
                a3=3*a_p(1,j)
                a4=3*a_p(2,j)-2
                a5=3*a_p(2,j)-1
                a6=3*a_p(2,j)
                a7=3*a_p(3,j)-2
                a8=3*a_p(3,j)-1
                a9=3*a_p(3,j)
              case default
                write(*,*)'Something went horribly wrong'
                exit
            end select
            ax=coord(a1)
            ay=coord(a2)
            az=coord(a3)
            bx=coord(a4)
            by=coord(a5)
            bz=coord(a6)
            cx=coord(a7)
            cy=coord(a8)
            cz=coord(a9)
            call ddangle(ax, ay, az, bx, by, bz, cx, cy, cz,
     1       dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz,
     1       daxax, daxay, daxaz, daxbx, daxby, daxbz, daxcx, daxcy,
     1       daxcz, dayay, dayaz, daybx, dayby, daybz, daycx, daycy,
     1       daycz, dazaz, dazbx, dazby, dazbz, dazcx, dazcy, dazcz,
     1       dbxbx, dbxby, dbxbz, dbxcx, dbxcy, dbxcz, dbyby, dbybz,
     1       dbycx, dbycy, dbycz, dbzbz, dbzcx, dbzcy, dbzcz, dcxcx,
     1       dcxcy, dcxcz, dcycy, dcycz, dczcz,
     1       angle_abc)
            diff=angle_abc - a_naught
            hessian(a1,a1)=hessian(a1,a1) + k * (dax*dax + daxax*diff)
            hessian(a1,a2)=hessian(a1,a2) + k * (dax*day + daxay*diff)
            hessian(a1,a3)=hessian(a1,a3) + k * (dax*daz + daxaz*diff)
            hessian(a1,a4)=hessian(a1,a4) + k * (dax*dbx + daxbx*diff)
            hessian(a1,a5)=hessian(a1,a5) + k * (dax*dby + daxby*diff)
            hessian(a1,a6)=hessian(a1,a6) + k * (dax*dbz + daxbz*diff)
            hessian(a1,a7)=hessian(a1,a7) + k * (dax*dcx + daxcx*diff)
            hessian(a1,a8)=hessian(a1,a8) + k * (dax*dcy + daxcy*diff)
            hessian(a1,a9)=hessian(a1,a9) + k * (dax*dcz + daxcz*diff)
            hessian(a2,a2)=hessian(a2,a2) + k * (day*day + dayay*diff)
            hessian(a2,a3)=hessian(a2,a3) + k * (day*daz + dayaz*diff)
            hessian(a2,a4)=hessian(a2,a4) + k * (day*dbx + daybx*diff)
            hessian(a2,a5)=hessian(a2,a5) + k * (day*dby + dayby*diff)
            hessian(a2,a6)=hessian(a2,a6) + k * (day*dbz + daybz*diff)
            hessian(a2,a7)=hessian(a2,a7) + k * (day*dcx + daycx*diff)
            hessian(a2,a8)=hessian(a2,a8) + k * (day*dcy + daycy*diff)
            hessian(a2,a9)=hessian(a2,a9) + k * (day*dcz + daycz*diff)
            hessian(a3,a3)=hessian(a3,a3) + k * (daz*daz + dazaz*diff)
            hessian(a3,a4)=hessian(a3,a4) + k * (daz*dbx + dazbx*diff)
            hessian(a3,a5)=hessian(a3,a5) + k * (daz*dby + dazby*diff)
            hessian(a3,a6)=hessian(a3,a6) + k * (daz*dbz + dazbz*diff)
            hessian(a3,a7)=hessian(a3,a7) + k * (daz*dcx + dazcx*diff)
            hessian(a3,a8)=hessian(a3,a8) + k * (daz*dcy + dazcy*diff)
            hessian(a3,a9)=hessian(a3,a9) + k * (daz*dcz + dazcz*diff)
            hessian(a4,a4)=hessian(a4,a4) + k * (dbx*dbx + dbxbx*diff)
            hessian(a4,a5)=hessian(a4,a5) + k * (dbx*dby + dbxby*diff)
            hessian(a4,a6)=hessian(a4,a6) + k * (dbx*dbz + dbxbz*diff)
            hessian(a4,a7)=hessian(a4,a7) + k * (dbx*dcx + dbxcx*diff)
            hessian(a4,a8)=hessian(a4,a8) + k * (dbx*dcy + dbxcy*diff)
            hessian(a4,a9)=hessian(a4,a9) + k * (dbx*dcz + dbxcz*diff)
            hessian(a5,a5)=hessian(a5,a5) + k * (dby*dby + dbyby*diff)
            hessian(a5,a6)=hessian(a5,a6) + k * (dby*dbz + dbybz*diff)
            hessian(a5,a7)=hessian(a5,a7) + k * (dby*dcx + dbycx*diff)
            hessian(a5,a8)=hessian(a5,a8) + k * (dby*dcy + dbycy*diff)
            hessian(a5,a9)=hessian(a5,a9) + k * (dby*dcz + dbycz*diff)
            hessian(a6,a6)=hessian(a6,a6) + k * (dbz*dbz + dbzbz*diff)
            hessian(a6,a7)=hessian(a6,a7) + k * (dbz*dcx + dbzcx*diff)
            hessian(a6,a8)=hessian(a6,a8) + k * (dbz*dcy + dbzcy*diff)
            hessian(a6,a9)=hessian(a6,a9) + k * (dbz*dcz + dbzcz*diff)
            hessian(a7,a7)=hessian(a7,a7) + k * (dcx*dcx + dcxcx*diff)
            hessian(a7,a8)=hessian(a7,a8) + k * (dcx*dcy + dcxcy*diff)
            hessian(a7,a9)=hessian(a7,a9) + k * (dcx*dcz + dcxcz*diff)
            hessian(a8,a8)=hessian(a8,a8) + k * (dcy*dcy + dcycy*diff)
            hessian(a8,a9)=hessian(a8,a9) + k * (dcy*dcz + dcycz*diff)
            hessian(a9,a9)=hessian(a9,a9) + k * (dcz*dcz + dczcz*diff)
          end do angles
        endif
      end do angle_types

c dihedrals
      dihedral_types: do i=1,4
        select case(iopt)
          case(1,2)
c           no dihedrals in case of iopt=1,2
            exit dihedral_types
          case(3,4)
        end select
        select case(i)
          case(1)
            k=fdhhh
            d_naught=dhhh
            m=nd_hhh
          case(2)
            k=fdhhp
            d_naught=dhhp
            m=nd_hhp
          case(3)
            k=fdhpp
            d_naught=dhpp
            m=nd_hpp
          case(4)
            k=fdppp
            d_naught=dppp
            m=nd_ppp
          case default
            write(*,*)'Something went horribly wrong'
            exit
        end select
        if(m.gt.0)then
          dihedrals: do j=1,m
            select case(i)
              case(1)
                a1=3*d_hhh(1,j)-2
                a2=3*d_hhh(1,j)-1
                a3=3*d_hhh(1,j)
                a4=3*d_hhh(2,j)-2
                a5=3*d_hhh(2,j)-1
                a6=3*d_hhh(2,j)
                a7=3*d_hhh(3,j)-2
                a8=3*d_hhh(3,j)-1
                a9=3*d_hhh(3,j)
                a10=3*d_hhh(4,j)-2
                a11=3*d_hhh(4,j)-1
                a12=3*d_hhh(4,j)
              case(2)
                a1=3*d_hhp(1,j)-2
                a2=3*d_hhp(1,j)-1
                a3=3*d_hhp(1,j)
                a4=3*d_hhp(2,j)-2
                a5=3*d_hhp(2,j)-1
                a6=3*d_hhp(2,j)
                a7=3*d_hhp(3,j)-2
                a8=3*d_hhp(3,j)-1
                a9=3*d_hhp(3,j)
                a10=3*d_hhp(4,j)-2
                a11=3*d_hhp(4,j)-1
                a12=3*d_hhp(4,j)
              case(3)
                a1=3*d_hpp(1,j)-2
                a2=3*d_hpp(1,j)-1
                a3=3*d_hpp(1,j)
                a4=3*d_hpp(2,j)-2
                a5=3*d_hpp(2,j)-1
                a6=3*d_hpp(2,j)
                a7=3*d_hpp(3,j)-2
                a8=3*d_hpp(3,j)-1
                a9=3*d_hpp(3,j)
                a10=3*d_hpp(4,j)-2
                a11=3*d_hpp(4,j)-1
                a12=3*d_hpp(4,j)
              case(4)
                a1=3*d_ppp(1,j)-2
                a2=3*d_ppp(1,j)-1
                a3=3*d_ppp(1,j)
                a4=3*d_ppp(2,j)-2
                a5=3*d_ppp(2,j)-1
                a6=3*d_ppp(2,j)
                a7=3*d_ppp(3,j)-2
                a8=3*d_ppp(3,j)-1
                a9=3*d_ppp(3,j)
                a10=3*d_ppp(4,j)-2
                a11=3*d_ppp(4,j)-1
                a12=3*d_ppp(4,j)
              case default
                write(*,*)'Something went horribly wrong'
                exit
            end select
            ax=coord(a1)
            ay=coord(a2)
            az=coord(a3)
            bx=coord(a4)
            by=coord(a5)
            bz=coord(a6)
            cx=coord(a7)
            cy=coord(a8)
            cz=coord(a9)
            dx=coord(a10)
            dy=coord(a11)
            dz=coord(a12)
            call dddihedral(
     1       ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz,
     1       dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz, ddx, ddy, ddz,
     1       daxax, daxay, daxaz, daxbx, daxby, daxbz, daxcx, daxcy,
     1       daxcz, daxdx, daxdy, daxdz, dayay, dayaz, daybx, dayby,
     1       daybz, daycx, daycy, daycz, daydx, daydy, daydz,
     1       dazaz, dazbx, dazby, dazbz, dazcx, dazcy, dazcz, dazdx,
     1       dazdy, dazdz, dbxbx, dbxby, dbxbz, dbxcx,
     1       dbxcy, dbxcz, dbxdx, dbxdy, dbxdz,
     1       dbyby, dbybz, dbycx, dbycy, dbycz, dbydx, dbydy, dbydz,
     1       dbzbz, dbzcx, dbzcy, dbzcz, dbzdx, dbzdy, dbzdz,
     1       dcxcx, dcxcy, dcxcz, dcxdx, dcxdy, dcxdz,
     1       dcycy, dcycz, dcydx, dcydy, dcydz,
     1       dczcz, dczdx, dczdy, dczdz, ddxdx, ddxdy, ddxdz,
     1       ddydy, ddydz, ddzdz,
     1       dihedral_abcd)
            diff=dihedral_abcd - d_naught
            hessian(a1,a1) =hessian(a1,a1)   + k *(dax*dax + daxax*diff)
            hessian(a1,a2) =hessian(a1,a2)   + k *(dax*day + daxay*diff)
            hessian(a1,a3) =hessian(a1,a3)   + k *(dax*daz + daxaz*diff)
            hessian(a1,a4) =hessian(a1,a4)   + k *(dax*dbx + daxbx*diff)
            hessian(a1,a5) =hessian(a1,a5)   + k *(dax*dby + daxby*diff)
            hessian(a1,a6) =hessian(a1,a6)   + k *(dax*dbz + daxbz*diff)
            hessian(a1,a7) =hessian(a1,a7)   + k *(dax*dcx + daxcx*diff)
            hessian(a1,a8) =hessian(a1,a8)   + k *(dax*dcy + daxcy*diff)
            hessian(a1,a9) =hessian(a1,a9)   + k *(dax*dcz + daxcz*diff)
            hessian(a1,a10)=hessian(a1,a10)  + k *(dax*ddx + daxdx*diff)
            hessian(a1,a11)=hessian(a1,a11)  + k *(dax*ddy + daxdy*diff)
            hessian(a1,a12)=hessian(a1,a12)  + k *(dax*ddz + daxdz*diff)
            hessian(a2,a2) =hessian(a2,a2)   + k *(day*day + dayay*diff)
            hessian(a2,a3) =hessian(a2,a3)   + k *(day*daz + dayaz*diff)
            hessian(a2,a4) =hessian(a2,a4)   + k *(day*dbx + daybx*diff)
            hessian(a2,a5) =hessian(a2,a5)   + k *(day*dby + dayby*diff)
            hessian(a2,a6) =hessian(a2,a6)   + k *(day*dbz + daybz*diff)
            hessian(a2,a7) =hessian(a2,a7)   + k *(day*dcx + daycx*diff)
            hessian(a2,a8) =hessian(a2,a8)   + k *(day*dcy + daycy*diff)
            hessian(a2,a9) =hessian(a2,a9)   + k *(day*dcz + daycz*diff)
            hessian(a2,a10)=hessian(a2,a10)  + k *(day*ddx + daydx*diff)
            hessian(a2,a11)=hessian(a2,a11)  + k *(day*ddy + daydy*diff)
            hessian(a2,a12)=hessian(a2,a12)  + k *(day*ddz + daydz*diff)
            hessian(a3,a3) =hessian(a3,a3)   + k *(daz*daz + dazaz*diff)
            hessian(a3,a4) =hessian(a3,a4)   + k *(daz*dbx + dazbx*diff)
            hessian(a3,a5) =hessian(a3,a5)   + k *(daz*dby + dazby*diff)
            hessian(a3,a6) =hessian(a3,a6)   + k *(daz*dbz + dazbz*diff)
            hessian(a3,a7) =hessian(a3,a7)   + k *(daz*dcx + dazcx*diff)
            hessian(a3,a8) =hessian(a3,a8)   + k *(daz*dcy + dazcy*diff)
            hessian(a3,a9) =hessian(a3,a9)   + k *(daz*dcz + dazcz*diff)
            hessian(a3,a10)=hessian(a3,a10)  + k *(daz*ddx + dazdx*diff)
            hessian(a3,a11)=hessian(a3,a11)  + k *(daz*ddy + dazdy*diff)
            hessian(a3,a12)=hessian(a3,a12)  + k *(daz*ddz + dazdz*diff)
            hessian(a4,a4) =hessian(a4,a4)   + k *(dbx*dbx + dbxbx*diff)
            hessian(a4,a5) =hessian(a4,a5)   + k *(dbx*dby + dbxby*diff)
            hessian(a4,a6) =hessian(a4,a6)   + k *(dbx*dbz + dbxbz*diff)
            hessian(a4,a7) =hessian(a4,a7)   + k *(dbx*dcx + dbxcx*diff)
            hessian(a4,a8) =hessian(a4,a8)   + k *(dbx*dcy + dbxcy*diff)
            hessian(a4,a9) =hessian(a4,a9)   + k *(dbx*dcz + dbxcz*diff)
            hessian(a4,a10)=hessian(a4,a10)  + k *(dbx*ddx + dbxdx*diff)
            hessian(a4,a11)=hessian(a4,a11)  + k *(dbx*ddy + dbxdy*diff)
            hessian(a4,a12)=hessian(a4,a12)  + k *(dbx*ddz + dbxdz*diff)
            hessian(a5,a5) =hessian(a5,a5)   + k *(dby*dby + dbyby*diff)
            hessian(a5,a6) =hessian(a5,a6)   + k *(dby*dbz + dbybz*diff)
            hessian(a5,a7) =hessian(a5,a7)   + k *(dby*dcx + dbycx*diff)
            hessian(a5,a8) =hessian(a5,a8)   + k *(dby*dcy + dbycy*diff)
            hessian(a5,a9) =hessian(a5,a9)   + k *(dby*dcz + dbycz*diff)
            hessian(a5,a10)=hessian(a5,a10)  + k *(dby*ddx + dbydx*diff)
            hessian(a5,a11)=hessian(a5,a11)  + k *(dby*ddy + dbydy*diff)
            hessian(a5,a12)=hessian(a5,a12)  + k *(dby*ddz + dbydz*diff)
            hessian(a6,a6) =hessian(a6,a6)   + k *(dbz*dbz + dbzbz*diff)
            hessian(a6,a7) =hessian(a6,a7)   + k *(dbz*dcx + dbzcx*diff)
            hessian(a6,a8) =hessian(a6,a8)   + k *(dbz*dcy + dbzcy*diff)
            hessian(a6,a9) =hessian(a6,a9)   + k *(dbz*dcz + dbzcz*diff)
            hessian(a6,a10)=hessian(a6,a10)  + k *(dbz*ddx + dbzdx*diff)
            hessian(a6,a11)=hessian(a6,a11)  + k *(dbz*ddy + dbzdy*diff)
            hessian(a6,a12)=hessian(a6,a12)  + k *(dbz*ddz + dbzdz*diff)
            hessian(a7,a7) =hessian(a7,a7)   + k *(dcx*dcx + dcxcx*diff)
            hessian(a7,a8) =hessian(a7,a8)   + k *(dcx*dcy + dcxcy*diff)
            hessian(a7,a9) =hessian(a7,a9)   + k *(dcx*dcz + dcxcz*diff)
            hessian(a7,a10)=hessian(a7,a10)  + k *(dcx*ddx + dcxdx*diff)
            hessian(a7,a11)=hessian(a7,a11)  + k *(dcx*ddy + dcxdy*diff)
            hessian(a7,a12)=hessian(a7,a12)  + k *(dcx*ddz + dcxdz*diff)
            hessian(a8,a8) =hessian(a8,a8)   + k *(dcy*dcy + dcycy*diff)
            hessian(a8,a9) =hessian(a8,a9)   + k *(dcy*dcz + dcycz*diff)
            hessian(a8,a10)=hessian(a8,a10)  + k *(dcy*ddx + dcydx*diff)
            hessian(a8,a11)=hessian(a8,a11)  + k *(dcy*ddy + dcydy*diff)
            hessian(a8,a12)=hessian(a8,a12)  + k *(dcy*ddz + dcydz*diff)
            hessian(a9,a9) =hessian(a9,a9)   + k *(dcz*dcz + dczcz*diff)
            hessian(a9,a10)=hessian(a9,a10)  + k *(dcz*ddx + dczdx*diff)
            hessian(a9,a11)=hessian(a9,a11)  + k *(dcz*ddy + dczdy*diff)
            hessian(a9,a12)=hessian(a9,a12)  + k *(dcz*ddz + dczdz*diff)
            hessian(a10,a10)=hessian(a10,a10)+ k *(ddx*ddx + ddxdx*diff)
            hessian(a10,a11)=hessian(a10,a11)+ k *(ddx*ddy + ddxdy*diff)
            hessian(a10,a12)=hessian(a10,a12)+ k *(ddx*ddz + ddxdz*diff)
            hessian(a11,a11)=hessian(a11,a11)+ k *(ddy*ddy + ddydy*diff)
            hessian(a11,a12)=hessian(a11,a12)+ k *(ddy*ddz + ddydz*diff)
            hessian(a12,a12)=hessian(a12,a12)+ k *(ddz*ddz + ddzdz*diff)
          end do dihedrals
        endif
      end do dihedral_types

c coulomb
      if(iopt.eq.2 .or. iopt.eq.4) then
        atoms: do j=1,N
          k=fco
          a1=3*N-2
          a2=3*N-1
          a3=3*N
          ax=coord(a1)
          ay=coord(a2)
          az=coord(a3)
          call ddcoulomb(ax, ay, az, dax, day, daz,
     1      daxax, daxay, daxaz, dayay, dayaz, dazaz, c)
          hessian(a1,a1)=hessian(a1,a1) + k * (dax*dax + daxax*c)
          hessian(a1,a2)=hessian(a1,a2) + k * (dax*day + daxay*c)
          hessian(a1,a3)=hessian(a1,a3) + k * (dax*daz + daxaz*c)
          hessian(a2,a2)=hessian(a2,a2) + k * (day*day + dayay*c)
          hessian(a2,a3)=hessian(a2,a3) + k * (day*daz + dayaz*c)
          hessian(a3,a3)=hessian(a3,a3) + k * (daz*daz + dazaz*c)
        end do atoms
      endif

c copy hessian to the other half
      do i=1,3*N
        do j=i+1,3*N
          hessian(i,j)=hessian(j,i)+hessian(i,j)
          hessian(j,i)=hessian(i,j)
        enddo
      enddo

      return
      END SUBROUTINE


      DOUBLE PRECISION FUNCTION SA_f1dimx(N,IOP,ier,
     1 x,pcom,xicom,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION pcom(N),xicom(N),xt(N)
      REAL*8 Dist(3,N)
      HUGE=dfloat(N)*1.d2
      do 11 j=1,N
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      If(IOP.eq.0) then
        Call SA_MDSnorm(N,AN,R,xt,Dist)
        sa_f1dimx=AN
      else
        Call SA_MAInorm(N,IP,AN,xt,Dist)
        sa_f1dimx=-AN
      endif
      if(AN.gt.Huge) then
        Print*,'**** Severe Error, check coordinates'
        ier=1
        stop
      endif
      return
      END


      SUBROUTINE SA_MAInorm(N,IP,dm,c,d)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION d(3,N),c(N)
C     Calculate minimum distance to center c
      dm=1.d8
      Do i=1,N
       dx=d(1,I)-c(1)
       dy=d(2,I)-c(2)
       dz=d(3,I)-c(3)
       dp=dsqrt(dx*dx+dy*dy+dz*dz)
       if(dp.lt.dm) then
        dm=dp
        IP=i
       endif
      enddo
      Return
      End

      SUBROUTINE SA_MDSnorm(N,A,davd,c,d)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION d(3,N),c(N),dp(N)
C     Calculate norm for minimum distance sphere
      XM=dfloat(N)
      dav=0.d0
      Do i=1,N
       dx=d(1,I)-c(1)
       dy=d(2,I)-c(2)
       dz=d(3,I)-c(3)
       dp(i)=dsqrt(dx*dx+dy*dy+dz*dz)
       dav=dav+dp(i)
      enddo
      davd=dav/XM
      ddav=0.d0
      Do i=1,N
       ddav=ddav+dabs(davd-dp(i))
      enddo
      A=ddav/XM
      Return
      End

      subroutine default_force_parameters(iopt,force)
      use config
      implicit none
      integer iopt
      real*8 force(ffmaxdim)
      real*8 fcoulomb, ftol
      real*8 WuR5, WuR6, WuA5, WuA6, WufR5, WufR6, WufA5, WufA6
      real*8 ExtWuR55, ExtWuR56, ExtWuR66,
     1 ExtWuA5, ExtWuA6,
     1 ExtWuD555, ExtWuD556, ExtWuD566, ExtWuD666,
     1 ExtWufR55, ExtWufR56, ExtWufR66,
     1 ExtWufA5, ExtWufA6,
     1 ExtWufD555, ExtWufD556, ExtWufD566, ExtWufD666

C     tolerance parameter (to be used in all force fields)
      fcoulomb=0.d0
      ftol=1.d-7

C     Defining the HO force field using Fowler force constants
C     Distances are taken in Angstroems and angles in rad
C     Force constants in N/m for distances and N/m A^2/rad^2 for angles (default values)
      WuR5=1.455d0              ! in angstroem from solid-state
      WuR6=1.391d0
      WuA5=1.08d2               ! in deg
      WuA6=1.20d2
      WufR5=390.7d0             ! from Ceulemans, Fowler
      WufR6=499.7d0
      WufA5=47.88d0*1.45d0**2
      WufA6=80.86d0*1.45d0*1.37d0

C     Defining an extension of the Wu force field (default values)
c     three distances: zero values
      ExtWuR55=1.479d0
      ExtWuR56=1.458d0
      ExtWuR66=1.401d0
c     two angles: zero values
      ExtWuA5=1.08d2
      ExtWuA6=1.20d2
c     four dihedrals: zero values (according to ideal_dihedral)
      ExtWuD555=37.38d0
      ExtWuD556=29.20d0
      ExtWuD566=23.49d0
      ExtWuD666=0.0d0
c     three distances: forces (let's assume they are all the same)
      ExtWufR55=260.d0
      ExtWufR56=390.d0
      ExtWufR66=450.d0
c     three angles: forces (let's assume they are all the same)
      ExtWufA5=100.d0
      ExtWufA6=100.d0
c     four dihedrals: forces (let's assume they are all the same)
      ExtWufD555=35.d0
      ExtWufD556=65.d0
      ExtWufD566=85.d0
      ExtWufD666=270.d0


      if(iopt.eq.1 .or. iopt.eq.2)then
C     Wu force field
         force(1)=WuR5
         force(2)=WuR6
         force(3)=WuA5
         force(4)=WuA6
         force(5)=WufR5
         force(6)=WufR6
         force(7)=WufA5
         force(8)=WufA6
         force(9)=fCoulomb
      else if(iopt.eq.3 .or. iopt.eq.4)then
C     ExtWu force field
         force(1)=ExtWuR55
         force(2)=ExtWuR56
         force(3)=ExtWuR66
         force(4)=ExtWuA5
         force(5)=ExtWuA6
         force(6)=ExtWuD555
         force(7)=ExtWuD556
         force(8)=ExtWuD566
         force(9)=ExtWuD666
         force(10)=ExtWufR55
         force(11)=ExtWufR56
         force(12)=ExtWufR66
         force(13)=ExtWufA5
         force(14)=ExtWufA6
         force(15)=ExtWufD555
         force(16)=ExtWufD556
         force(17)=ExtWufD566
         force(18)=ExtWufD666
         force(19)=fCoulomb
      endif

      END SUBROUTINE

