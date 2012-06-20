C---- Redundant Subroutine -----
      SUBROUTINE Mincovsphere1(ndim,M,IOUT,Dist,Rmin,Rmax,
     1 VCS,ACS,Atol,VTol,c,radius,RVdW)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,ndim)
      DIMENSION u(ndim),c(3),Kappa(ndim)
C     Get the minimum covering sphere using algorithm 1 by E.A.Yildirim
C     Less sufficient than algorithm 1 but shorter code
      DATA API/3.14159265358979d0/
      DATA itermax/10000000/
      DATA eps/1.d-7/
      epsilon=(1.d0+eps)**2-1.d0
      Write(IOUT,1000) M,epsilon
 1000 Format(1X,'Calculate the minimum covering sphere using',
     1 ' the algorithm 1 given by E.A.Yildirim,'
     2 ' SIAM Journal on Optimization 19(3), 1368 (2008)',
     3 //1X,' Number of points in 3-D space: ',I4,
     4 //1X,' lower limit: ',d12.6)
C
C    Initial step
C
      do i=1,ndim
        u(i)=0.d0
      enddo
      xmax=0.d0
      ialpha=0
C    Algorithm changed here by taking the furthest 
C    point from the center of points
      do i=1,M
        dx=Dist(1,i)
        dy=Dist(2,i)
        dz=Dist(3,i)
C       dx=Dist(1,1)-Dist(1,i)
C       dy=Dist(2,1)-Dist(2,i)
C       dz=Dist(3,1)-Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.gt.xmax) then
            ialpha=i
            xmax=Xnorm
          endif
      enddo
      xmax=0.d0
      ibeta=0
      do i=1,M
        dx=Dist(1,ialpha)-Dist(1,i)
        dy=Dist(2,ialpha)-Dist(2,i)
        dz=Dist(3,ialpha)-Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.gt.xmax) then
            ibeta=i
            xmax=Xnorm
          endif
      enddo
      u(ialpha)=.5d0
      u(ibeta)=.5d0
C     Size of set Kappa (number of vectors)
      nkappa=2
C     Kappa(1)=ialpha
C     Kappa(2)=ibeta
      do j=1,3
        c(j)=.5d0*(Dist(j,ialpha)+Dist(j,ibeta))
      enddo
      gammak=fpsi(ndim,M,u,dist)
      xmax=0.d0
      ikappa=0
      do i=1,M
        dx=c(1)-Dist(1,i)
        dy=c(2)-Dist(2,i)
        dz=c(3)-Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.gt.xmax) then
            ikappa=i
            xmax=Xnorm
          endif
      enddo
      dx=c(1)-Dist(1,ikappa)
      dy=c(2)-Dist(2,ikappa)
      dz=c(3)-Dist(3,ikappa)
      delta=(dx*dx+dy*dy+dz*dz)/gammak-1.d0
      k=0
      Write(IOUT,1001) k,nkappa,delta,(c(i),i=1,3)
C     Now loop
C     ----------------
      iter10=1
      do iter=1,itermax
      if(delta.le.epsilon) Go to 10
        alamdak=.5d0*delta/(1.d0+delta)
        k=k+1
        blambdak=1.d0-alamdak
          do i=1,M
            u(i)=blambdak*u(i)
            if(i.eq.ikappa) u(i)=u(i)+alamdak
          enddo
          do j=1,3
            c(j)=blambdak*c(j)+alamdak*dist(j,ikappa)
          enddo
        itest=0
        do i=1,nkappa
        if(ikappa.eq.Kappa(i)) itest=1
        enddo
        if(itest.eq.0) then
        nkappa=nkappa+1
        Kappa(nkappa)=ikappa
        endif
C     New parameters
      gammak=fpsi(ndim,M,u,dist)
      xmax=0.d0
      ikappa=0
      do i=1,M
        dx=c(1)-Dist(1,i)
        dy=c(2)-Dist(2,i)
        dz=c(3)-Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.gt.xmax) then
            ikappa=i
            xmax=Xnorm
          endif
      enddo
      dx=c(1)-Dist(1,ikappa)
      dy=c(2)-Dist(2,ikappa)
      dz=c(3)-Dist(3,ikappa)
      delta=(dx*dx+dy*dy+dz*dz)/gammak-1.d0
      if(nkappa.eq.0.or.nkappa.gt.M) go to 10
      enddo
C     ----------------
C     End of loop
C     
   10 if(iter.ge.itermax) Write(IOUT,1006) itermax
      Write(IOUT,1005) k,nkappa,delta,(c(i),i=1,3)
      rmcs=0.d0
      Do i=1,M
      RI=dsqrt((Dist(1,I)-c(1))**2+(Dist(2,I)-c(2))**2
     1 +(Dist(3,I)-c(3))**2)
      if (ri.gt.rmcs) rmcs=ri
      enddo
      dist0=dsqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
      Write(IOUT,1002) (c(i),i=1,3),dist0
C     Check if there is no point outside the sphere
C     Finally calculate the surface and volume and compare to previous results
      VMCS=4.d0/3.d0*Api*RMCS**3
      AMCS=4.d0*Api*RMCS**2
      RatioMCS=AMCS/VMCS
      RatioCS=ACS/VCS
      RatioT=Atol/Vtol
      Write(IOUT,1004) VMCS,VCS,Vtol,AMCS,ACS,Atol,
     1 RatioMCS,RatioCS,RatioT
C     Do statistics
      keq=0
      kgt=0
      klo=0
      rmcsl=-10.d0*epsilon*rmcs
      rmcsu=10.d0*epsilon*rmcs
      Do i=1,M
      RZ=dsqrt((Dist(1,I)-c(1))**2+(Dist(2,I)-c(2))**2
     1 +(Dist(3,I)-c(3))**2)-rmcs
      if (rz.lt.rmcsl) klo=klo+1
      if (rz.gt.rmcsu) kgt=kgt+1
      if (rz.le.rmcsu.and.rz.ge.rmcsl) then
      keq=keq+1
      RI=dsqrt((Dist(1,I)-c(1))**2+(Dist(2,I)-c(2))**2
     1 +(Dist(3,I)-c(3))**2)
      if(keq.le.4) Write(IOUT,1003) i,(Dist(j,i),j=1,3),RI
      endif
      enddo
      radius=dsqrt((1.d0+delta)*gammak)
      Write(IOUT,1008) radius
      Write(IOUT,1007) M,keq,klo,kgt
 1001 Format(1x,' Cycle ',I4,', nr. points: ',I4,', zero-condition: ',
     1 F14.8,', Center of Sphere (X,Y,Z): ',3(D12.6,2X))
 1002 Format(1x,' End of iteration',
     1 /1X,' Center of minimum covering sphere (x,y,z):',3(D14.8,2X),
     2 ' at distance ',D14.8,2X,'from origin',
     3 /1X,' Points defining the sphere (x,y,z) ',
     4 'and distances from center:')
 1003 Format(1x,' P',I4,'= ',3(D14.8,2X),'  RCP = ',D14.8)
 1004 Format(1x,' Final result and comparison (in units of input):',
     1 /1x,' Volume of minimum covering sphere: ',D14.8,
     2 /1x,' Volume of central covering sphere: ',D14.8,
     3 /1x,' Volume of cage molecule          : ',D14.8,
     4 /1x,' Area of minimum covering sphere: ',D14.8,
     5 /1x,' Area of central covering sphere: ',D14.8,
     6 /1x,' Area of cage molecule          : ',D14.8,
     4 /1x,' Ratio area/volume of minimum covering sphere: ',D14.8,
     5 /1x,' Ratio area/volume of central covering sphere: ',D14.8,
     6 /1x,' Ratio area/volume of cage molecule          : ',D14.8)
 1005 Format(1x,' Number of cycles ',I8,', nr. points: ',I4,
     1 ', convergence: ',F14.8,', Center of Sphere (X,Y,Z): ',
     2 3(D12.6,2X))
 1006 Format(1x,' Maximum number of allowed iterations ',I8,' reached')
 1007 Format(1x,' Final statistics:',/3X,
     1 ' Number of points: ',I4,/3X,
     2 ' Number of points on the sphere: ',I4,/3X,
     3 ' Number of points below the sphere: ',I4,/3X,
     4 ' Number of points above the sphere: ',I4)
 1008 Format(1x,' Radius of minimum covering sphere: ',D14.8)
      return
      END
