      SUBROUTINE OptGraph(N, IOP,Iout,IDA,IS,MDist,maxl,scalePPG,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C  This subroutine optimizes the fullerene graph using spring embedding
      DIMENSION Dist(2,N)
      DIMENSION IDA(N,N),IS(6),MDist(N,N)
      Data Rdist,ftol,conv/1.d0,.5d-10,1 6.0221367d-3/
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
      CALL frprmn2d(N,IOP,IDA,Iout,IS,MDist,
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

      SUBROUTINE frprmn2d(N,IOP,AH,Iout,IS,MDist,
     1 maxd,p,ftol,iter,fret,E0,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=500,EPS=1.d-10)
      Real*8 p(N*2),g(N*2),h(N*2),xi(N*2)
      Real*8 pcom(N*2),xicom(N*2)
      Integer AH(N,N),IS(6),MDist(N,N)
C     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
C     is performed on a function func3d, using its gradient as calculated by a routine dfunc3d.
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
      CALL func2d(N,IOP,AH,IS,MDist,maxd,p,fp,RAA)
       E0=fp
      Write(Iout,1003) E0
C     dfunc3d input vector p of length 2*N, output gradient of length 2*N user defined
      CALL dfunc2d(N,IOP,AH,IS,MDist,maxd,p,xi,RAA)
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
        CALL linmin2d(N,IOP,Iout,AH,IS,MDist,maxd,
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
        CALL dfunc2d(N,IOP,AH,IS,MDist,maxd,p,xi,RAA)
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

      SUBROUTINE linmin2d(N,IOP,Iout,AH,IS,MDist,
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
      CALL mnbrak2d(N,IOP,Iout,AH,IS,MDist,maxd,
     1 ax,xx,bx,fa,fx,fb,xicom,pcom,RAA)
      CALL brent2d(N,IOP,Iout,AH,IS,MDist,maxd,
     1 fret,ax,xx,bx,TOL,xmin,xicom,pcom,RAA)
      do j=1,2*N
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
      enddo
      return
      END

      SUBROUTINE f1dim2d(N,IOP,A,IS,MDist,maxd,
     1 f1dimf,x,xicom,pcom,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 pcom(N*2),xt(N*2),xicom(N*2)
      Integer A(N,N),IS(6),MDist(N,N)
C     USES func2d
      do j=1,2*N
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      CALL func2d(N,IOP,A,IS,MDist,maxd,xt,f1dimf,RAA)
      return
      END

      SUBROUTINE mnbrak2d(N,IOP,Iout,AH,IS,DD,maxd,
     1 ax,bx,cx,fa,fb,fc,xicom,pcom,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (GOLD=1.618034d0,GLIMIT=1.d2,TINY=1.d-20)
      Integer AH(N,N),IS(6)
      Integer DD(N,N)
      REAL*8 pcom(N*2),xicom(N*2)
      CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fa,ax,xicom,pcom,
     1 RAA)
      CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fb,bx,xicom,pcom,
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
      CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fc,cx,xicom,pcom,
     1 RAA)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
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
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        else if((cx-u)*(u-ulim).gt.0.)then
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
     1 RAA)
        else
          u=cx+GOLD*(cx-bx)
        if(u.gt.1.d10) then
        Write(Iout,1000)
        return
        endif
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
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

      SUBROUTINE brent2d(N,IOP,Iout,AH,IS,DD,maxd,
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
      CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fx,x,xicom,pcom,
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
        CALL f1dim2d(N,IOP,AH,IS,DD,maxd,fu,u,xicom,pcom,
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



      SUBROUTINE func2d(N,IOP,A,IS,MDist,maxd,p,fc,RAA)
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

      SUBROUTINE dfunc2d(N,IOP,A,IS,MDist,maxd,p,x,RAA)
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
