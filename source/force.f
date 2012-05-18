      SUBROUTINE func(IOP,NMAX,NMAX3,MMAX,n,IERR,A,N5,N6,N5M,N6M,
     1 p,fc,c)
      IMPLICIT REAL*8 (A-H,O-Z)

C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, energy
      Real*8 p(NMAX3),c(9)
      Integer A(NMAX,NMAX)
      Integer N5M(MMAX,5),N6M(MMAX,6)
      IERR=0
      rp=c(1)
      rh=c(2)
      ap=c(3)
      ah=c(4)
      frp=c(5)
      frh=c(6)
      fap=c(7)
      fah=c(8)
      fco=c(9)
      coulomb=0.d0

C     Stretching
      ehookrp=0.d0
      ehookrh=0.d0
      Do I=1,n,3
        I1=(I+2)/3
      Do J=I+3,n,3
        J1=(J+2)/3
        if(A(I1,J1).ne.0) then
         px=p(I)-p(J)
         py=p(I+1)-p(J+1)
         pz=p(I+2)-p(J+2)
         ratom=dsqrt(px*px+py*py+pz*pz)
         ratominv=1.d0/ratom
C        Check if bond is part of 5-ring
         ibond=6
          do IB=1,12
           ir1=0
           ir2=0
           do JB=1,5
            if(I1.eq.N5M(IB,JB)) ir1=1
            if(J1.eq.N5M(IB,JB)) ir2=1
           enddo
            if(ir1.eq.1.and.ir2.eq.1) then
C           5-ring
             ehookrp=ehookrp+(ratom-rp)**2
             go to 1
            endif
           enddo
C           6-ring
           ehookrh=ehookrh+(ratom-rh)**2
        endif
  1   continue
      enddo
      enddo

C     Bending
C     Loop over 5-rings
      ehookap=0.d0
      Do I=1,N5
      Do J=1,5
        JLX=J-1
        JRX=J+1
        if(JLX.eq.0) JLX=5
        if(JRX.eq.6) JRX=1
        JM=3*N5M(I,J)-2
        JL=3*N5M(I,JLX)-2
        JR=3*N5M(I,JRX)-2
         pxL=p(JM)  -p(JL)
         pyL=p(JM+1)-p(JL+1)
         pzL=p(JM+2)-p(JL+2)
        r2L=pxL*pxL+pyL*pyL+pzL*pzL
        r1L=dsqrt(r2L)
         pxR=p(JM)  -p(JR)
         pyR=p(JM+1)-p(JR+1)
         pzR=p(JM+2)-p(JR+2)
        r2R=pxR*pxR+pyR*pyR+pzR*pzR
        r1R=dsqrt(r2R)
         pxM=p(JL)  -p(JR)
         pyM=p(JL+1)-p(JR+1)
         pzM=p(JL+2)-p(JR+2)
        r2M=pxM*pxM+pyM*pyM+pzM*pzM
         cosarg=.5d0*(r2L+r2R-r2M)/(r1L*r1R)
         if(cosarg.ge.1.d0.or.cosarg.le.-1.d0) then
         IERR=1
         return
         endif
         anglep=dacos(cosarg)
         ehookap=ehookap+(anglep-ap)**2
      enddo
      enddo

C     Loop over 6-rings
      ehookah=0.d0
      if(N6.eq.0) go to 2
        Do I=1,N6
        Do J=1,6
        JLX=J-1
        JRX=J+1
        if(JLX.eq.0) JLX=6
        if(JRX.eq.7) JRX=1
        JM=3*N6M(I,J)  -2
        JL=3*N6M(I,JLX)-2
        JR=3*N6M(I,JRX)-2
         pxL=p(JM)  -p(JL)
         pyL=p(JM+1)-p(JL+1)
         pzL=p(JM+2)-p(JL+2)
        r2L=pxL*pxL+pyL*pyL+pzL*pzL
        r1L=dsqrt(r2L)
         pxR=p(JM)  -p(JR)
         pyR=p(JM+1)-p(JR+1)
         pzR=p(JM+2)-p(JR+2)
        r2R=pxR*pxR+pyR*pyR+pzR*pzR
        r1R=dsqrt(r2R)
         pxM=p(JL)  -p(JR)
         pyM=p(JL+1)-p(JR+1)
         pzM=p(JL+2)-p(JR+2)
        r2M=pxM*pxM+pyM*pyM+pzM*pzM
        cosarg=.5d0*(r2L+r2R-r2M)/(r1L*r1R)
         if(cosarg.ge.1.d0.or.cosarg.le.-1.d0) then
         IERR=1
         return
         endif
        angleh=dacos(cosarg)
        ehookah=ehookah+(angleh-ah)**2
      enddo
      enddo

C     Coulomb repulsion from origin
      if (fco.ne.0.d0)  then
       Do I=1,n,3
        rinv=1.d0/dsqrt(p(I)**2+p(I+1)**2+p(I+2)**2)
        coulomb=coulomb+rinv
       enddo
      endif

C     total energy  
  2   fc=frp*ehookrp+frh*ehookrh+fap*ehookap+fah*ehookah+fco*coulomb
      Return
      END

      SUBROUTINE dfunc(IOP,NMAX,NMAX3,MMAX,n,A,N5,N6,N5M,N6M,p,x,c)
      IMPLICIT REAL*8 (A-H,O-Z)
C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, gradient
      Real*8 p(NMAX3),x(NMAX3),c(9)
      Integer A(NMAX,NMAX)
      Integer N5M(MMAX,5),N6M(MMAX,6)
      rp=c(1)
      rh=c(2)
      ap=c(3)
      ah=c(4)
      frp=c(5)
      frh=c(6)
      fap=c(7)
      fah=c(8)
      fco=c(9)
C     Stretching
      Do I=1,n,3
        ehookx=0.d0
        ehooky=0.d0
        ehookz=0.d0
        I1=(I+2)/3
      Do J=1,n,3
        J1=(J+2)/3
        if(A(I1,J1).ne.0) then
          px=p(I)-p(J)
          py=p(I+1)-p(J+1)
          pz=p(I+2)-p(J+2)
          ratom=dsqrt(px*px+py*py+pz*pz)
          ratominv=1.d0/ratom
C         Check if bond is part of 5-ring
          do IB=1,12
           ir1=0
           ir2=0
           do JB=1,5
            if(I1.eq.N5M(IB,JB)) ir1=1
            if(J1.eq.N5M(IB,JB)) ir2=1
           enddo
           if(ir1.eq.1.and.ir2.eq.1) then
C           5-ring
            fac=frp*ratominv*(ratom-rp)
            ehookx=ehookx+fac*px
            ehooky=ehooky+fac*py
            ehookz=ehookz+fac*pz
            go to 1
           endif
           enddo
C           6-ring
            fac=frh*ratominv*(ratom-rh)
            ehookx=ehookx+fac*px
            ehooky=ehooky+fac*py
            ehookz=ehookz+fac*pz
        endif
  1   continue
      enddo
        x(I)  =2.d0*ehookx
        x(I+1)=2.d0*ehooky
        x(I+2)=2.d0*ehookz
      enddo
        
C     Bending
C     Loop over 5-rings
      Do I=1,N5
      Do J=1,5
        JLX=J-1
        JRX=J+1
        if(JLX.eq.0) JLX=5
        if(JRX.eq.6) JRX=1
        JM=3*N5M(I,J)-2
        JL=3*N5M(I,JLX)-2
        JR=3*N5M(I,JRX)-2
         pxL=p(JM)  -p(JL)
         pyL=p(JM+1)-p(JL+1)
         pzL=p(JM+2)-p(JL+2)
        r2L=pxL*pxL+pyL*pyL+pzL*pzL
        r1L=dsqrt(r2L)
        r3L=r1L*r2L
         pxR=p(JM)  -p(JR)
         pyR=p(JM+1)-p(JR+1)
         pzR=p(JM+2)-p(JR+2)
        r2R=pxR*pxR+pyR*pyR+pzR*pzR
        r1R=dsqrt(r2R)
        r3R=r1R*r2R
         pxM=p(JL)  -p(JR)
         pyM=p(JL+1)-p(JR+1)
         pzM=p(JL+2)-p(JR+2)
        r2M=pxM*pxM+pyM*pyM+pzM*pzM
        r1M=dsqrt(r2M)
        r3M=r1M*r2M
         cosarg=.5d0*(r2L+r2R-r2M)/(r1L*r1R)
         if(cosarg.gt.1.d0) cosarg=1.d0
         if(cosarg.lt.-1.d0) cosarg=-1.d0
         anglep=dacos(cosarg)
         anglesin=dabs(dsin(anglep))
         fac=fap*(anglep-ap)/anglesin
C     Derivative of central atom
         fac1=fac/(r3R*r3L)
         r2RL=r2R-r2L
         r2LR=-r2RL
         fac2=r2RL-r2M
         fac3=r2LR-r2M
         fac4=r2R*fac2
         fac5=r2L*fac3
        x(JM)  =x(JM)  +fac1*(pxL*fac4+pxR*fac5)
        x(JM+1)=x(JM+1)+fac1*(pyL*fac4+pyR*fac5)
        x(JM+2)=x(JM+2)+fac1*(pzL*fac4+pzR*fac5)
C     Derivative of left atom
         fac6=-fac/(r3L*r1R)
        x(JL)  =x(JL)  +fac6*(pxL*fac3-2.d0*pxM*r2L)
        x(JL+1)=x(JL+1)+fac6*(pyL*fac3-2.d0*pyM*r2L)
        x(JL+2)=x(JL+2)+fac6*(pzL*fac3-2.d0*pzM*r2L)
C     Derivative of right atom
         fac7=-fac/(r3R*r1L)
        x(JR)  =x(JR)  +fac7*(pxR*fac2+2.d0*pxM*r2R)
        x(JR+1)=x(JR+1)+fac7*(pyR*fac2+2.d0*pyM*r2R)
        x(JR+2)=x(JR+2)+fac7*(pzR*fac2+2.d0*pzM*r2R)
      enddo
      enddo
      
C     Loop over 6-rings
      if(N6.eq.0) return
      Do I=1,N6
      Do J=1,6
        JLX=J-1
        JRX=J+1
        if(JLX.eq.0) JLX=6
        if(JRX.eq.7) JRX=1
        JM=3*N6M(I,J)-2
        JL=3*N6M(I,JLX)-2
        JR=3*N6M(I,JRX)-2
         pxL=p(JM)  -p(JL)
         pyL=p(JM+1)-p(JL+1)
         pzL=p(JM+2)-p(JL+2)
        r2L=pxL*pxL+pyL*pyL+pzL*pzL
        r1L=dsqrt(r2L)
        r3L=r1L*r2L
         pxR=p(JM)  -p(JR)
         pyR=p(JM+1)-p(JR+1)
         pzR=p(JM+2)-p(JR+2)
        r2R=pxR*pxR+pyR*pyR+pzR*pzR
        r1R=dsqrt(r2R)
        r3R=r1R*r2R
         pxM=p(JL)  -p(JR)
         pyM=p(JL+1)-p(JR+1)
         pzM=p(JL+2)-p(JR+2)
        r2M=pxM*pxM+pyM*pyM+pzM*pzM
        r1M=dsqrt(r2M)
        r3M=r1M*r2M
         cosarg=.5d0*(r2L+r2R-r2M)/(r1L*r1R)
         if(cosarg.gt.1.d0) cosarg=1.d0
         if(cosarg.lt.-1.d0) cosarg=-1.d0
         angleh=dacos(cosarg)
         anglesin=dabs(dsin(angleh))
         fac=fah*(angleh-ah)/anglesin
C     Derivative of central atom
         fac1=fac/(r3R*r3L)
         r2RL=r2R-r2L
         r2LR=-r2RL
        fac2=r2RL-r2M
        fac3=r2LR-r2M
        fac4=r2R*fac2
        fac5=r2L*fac3
        x(JM)  =x(JM)  +fac1*(pxL*fac4+pxR*fac5)
        x(JM+1)=x(JM+1)+fac1*(pyL*fac4+pyR*fac5)
        x(JM+2)=x(JM+2)+fac1*(pzL*fac4+pzR*fac5)
C     Derivative of left atom
         fac6=-fac/(r3L*r1R)
        x(JL)  =x(JL)  +fac6*(pxL*fac3-2.d0*pxM*r2L)
        x(JL+1)=x(JL+1)+fac6*(pyL*fac3-2.d0*pyM*r2L)
        x(JL+2)=x(JL+2)+fac6*(pzL*fac3-2.d0*pzM*r2L)
C     Derivative of right atom
         fac7=-fac/(r3R*r1L)
        x(JR)  =x(JR)  +fac7*(pxR*fac2+2.d0*pxM*r2R)
        x(JR+1)=x(JR+1)+fac7*(pyR*fac2+2.d0*pyM*r2R)
        x(JR+2)=x(JR+2)+fac7*(pzR*fac2+2.d0*pzM*r2R)
      enddo
      enddo

C     Coulomb repulsion from origin
      if (fco.ne.0.d0)  then
       Do I=1,n,3
        rinv=(p(I)**2+p(I+1)**2+p(I+2)**2)**(-1.5d0)
        x(I)=x(I)-fco*rinv*p(I)
        x(I+1)=x(I+1)-fco*rinv*p(I+1)
        x(I+2)=x(I+2)-fco*rinv*p(I+2)
       enddo
      endif

      return
      END
