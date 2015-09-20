      SUBROUTINE MaxInSphere(IOUT,imcs,Dist,cmcs,RVdWC)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax),c(3),cmax(3),cmcs(3)
      DATA itermax/10000000/
      DATA eps,epsc/1.d-8,2.d0/
      IOP=1
C     Get the maximum inner sphere
      Write(IOUT,1000) 
C    Initial step
      do i=1,3
       c(i)=cmcs(i)
      enddo
      cnorm=dsqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
      nchoice=0
      if(cnorm.gt.epsc) then
       c(1)=0.d0
       c(2)=0.d0
       c(3)=0.d0
      nchoice=1
      endif
      Call MAInorm(3,IP,RMIS,c,Dist)
      if(nchoice.eq.0) then
      Write(IOUT,1001) RMIS,IP,(c(i),i=1,3)
      else
      Write(IOUT,1006) RMIS,IP,(c(i),i=1,3)
      endif
      VMIS=4.d0/3.d0*dpi*RMIS**3
      AMIS=4.d0*dpi*RMIS**2
      Write(IOUT,1002) RMIS,VMIS,AMIS
C     Start Iteration
      CALL powell(3,iter,Iout,IOP,ier,eps,AN,RMIS,c,cmax,Dist)
C     End Iteration
      if(ier.eq.1) then
       Write(IOUT,1010)
       Return
      endif
      Call MAInorm(3,IP,RMIS,c,Dist)
      Write(IOUT,1003) RMIS,IP,(c(i),i=1,3)
      VMIS=4.d0/3.d0*dpi*RMIS**3
      AMIS=4.d0*dpi*RMIS**2
      Write(IOUT,1004) RMIS,VMIS,AMIS
      if(imcs.eq.0) then
       RealMIS=RMIS-RVdWC
       VVdWC=4.d0/3.d0*dpi*RealMIS**3
       Write(IOUT,1005) RVdWC,RealMIS,VVdWC
      endif
 1000 Format(/1X,'Calculate the maximum inner sphere')
 1001 Format(1X,'Initial inner radius: ',d12.6,' to point ',I5,
     1 ' taken from center of MCS at (X,Y,Z): ',3(D14.8,2X))
 1002 Format(/1x,'Initial values (in units of input):',
     1 /1x,' Radius of maximum inner sphere: ',D14.8,
     1 /1x,' Volume of maximum inner sphere: ',D14.8,
     1 /1x,' Area   of maximum inner sphere: ',D14.8)
 1003 Format(1X,'Inner radius: ',d12.6,' to point ',I5,
     1 ' taken from center of MIS at (X,Y,Z): ',3(D14.8,2X))
 1004 Format(/1x,'Final values (in units of input):',
     1 /1x,' Radius of maximum inner sphere: ',D14.8,
     1 /1x,' Volume of maximum inner sphere: ',D14.8,
     1 /1x,' Area   of maximum inner sphere: ',D14.8)
 1005 Format(/1x,'Subtracting Van der Waals radius of carbon ',
     1 F8.4,' gives VdW inner sphere radius of ',D14.8,
     1 ' and volume of ',D14.8)
 1006 Format(1X,'Initial inner radius: ',d12.6,' to point ',I5,
     1 ' taken from barycenter at (X,Y,Z): ',3(D14.8,2X))
 1010 Format(/1x,'**** Error in MaxInSphere, Problem ill-defined')
      RETURN
      END

      SUBROUTINE MinDistSphere(IOUT,Dist,cMCS)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax),c(3),cmax(3),cMCS(3)
C     Get the minimum distance sphere
      
      DATA itermax/10000000/
      DATA eps/1.d-8/

      Write(IOUT,1000) 

C    Initial step
C    Center of MDS needs to be confined in convex hull
C    Here weaker test is carried out if it is within a cube
C    enclosing the fullerine of diameter xmax*2, ymax*2 and zmax*2
C    contained in the vector cmax
C    Get maximum values
      Do J=1,3
       cmax(J)=0.d0
       c(J)=0.d0
      enddo
      Do I=1,Nmax
      Do J=1,3
       if(dabs(Dist(J,I)).gt.cmax(J)) cmax(J)=dabs(Dist(J,I))
      enddo
      enddo
C    Calculate the MDS norm

      IOP=0
      Call MDSnorm(3,AN,RMDSI,c,Dist)
      Write(IOUT,1001) RMDSI,(c(i),i=1,3)
      VMDSI=4.d0/3.d0*dpi*RMDSI**3
      AMDSI=4.d0*dpi*RMDSI**2
      VMDS=VMDSI
      AMDS=AMDSI
      RMDS=RMDSI
      Write(IOUT,1002) RMDSI,VMDS,AMDS
      Write(IOUT,1003) AN
C     Start Iteration
      CALL powell(3,iter,Iout,IOP,ier,eps,AN,RMDSI,c,cmax,Dist)
C     End Iteration
      if(ier.eq.2) then
       Write(IOUT,1006)
       do I=1,3
        c(i)=cMCS(i)
       enddo
      endif
      if(ier.eq.1) then
       Write(IOUT,1007)
       Return
      endif
      Call MDSnorm(3,AN,RMDSI,c,Dist)
      Write(IOUT,1004) iter,(c(i),i=1,3),AN
      VMDSI=4.d0/3.d0*dpi*RMDSI**3
      AMDSI=4.d0*dpi*RMDSI**2
      VMDS=VMDSI
      AMDS=AMDSI
      RMDS=RMDSI
      rminMDS=1.d8
      Do i=1,number_vertices
       dx=dist(1,I)-c(1)
       dy=dist(2,I)-c(2)
       dz=dist(3,I)-c(3)
       dp=dsqrt(dx*dx+dy*dy+dz*dz)
       if(dp.lt.rminMDS) rminMDS=dp
      enddo
      distortion=100.d0*AN/rminMDS
      Write(IOUT,1005) RMDSI,VMDS,AMDS,rminMDS,distortion
 1000 Format(/1X,'Calculate the minimum distance sphere')
 1001 Format(1X,'Initial average radius: ',d12.6,
     1 ' taken from barycenter at (X,Y,Z): ',3(D14.8,2X))
 1002 Format(/1x,'Initial values (in units of input):',
     1 /1x,' Radius of minimum distance sphere: ',D14.8,
     1 /1x,' Volume of minimum distance sphere: ',D14.8,
     1 /1x,' Area of  minimum distance  sphere: ',D14.8)
 1003 Format(1x,' Norm: ',D14.8)
 1004 Format(/1x,'Final result after ',I4,' iterations',
     1 ' (in units of input): '
     1 /1X,'Center (x,y,z) of MDS: ',3(D14.8,2X),
     1 3X,'Norm: ',D14.8)
 1005 Format(/1x,'Final values (in units of input):',
     1 /1x,' Radius of minimum distance sphere:   ',D14.8,
     1 /1x,' Volume of minimum distance sphere:   ',D14.8,
     1 /1x,' Area of  minimum distance  sphere:   ',D14.8,
     1 /1x,' Smallest point distance from center: ',D14.8,
     1 /1x,' MDS distortion parameter in percent: ',D14.8)
 1006 Format(/1x,'Center of MDS tries to move out of convex hull.',
     1 ' MDS problem ill-defined due to large distortion.',
     1 /1x,'Taking center of MCS instead:')
 1007 Format(/1x,'**** Error in MaxInSphere, Problem ill-defined')
      Return
      End

      SUBROUTINE MAInorm(ncom,IP,dm,c,d)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION d(3,Nmax),c(ncom)
C     Calculate minimum distance to center c
      dm=1.d8
      Do i=1,number_vertices
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

      SUBROUTINE MDSnorm(ncom,A,davd,c,d)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION d(3,Nmax),c(ncom),dp(Nmax)
C     Calculate norm for minimum distance sphere
      XM=dfloat(number_vertices)
      dav=0.d0
      Do i=1,number_vertices
       dx=d(1,I)-c(1)
       dy=d(2,I)-c(2)
       dz=d(3,I)-c(3)
       dp(i)=dsqrt(dx*dx+dy*dy+dz*dz)
       dav=dav+dp(i)
      enddo
      davd=dav/XM
      ddav=0.d0
      Do i=1,number_vertices
       ddav=ddav+dabs(davd-dp(i)) 
      enddo
      A=ddav/XM
      Return
      End

      SUBROUTINE ProjectSphere(ipsphere,Iout,IAtom,
     1 IC3,Dist,c,rmcs,filename,El,TEXTINPUT)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     This routine projects vertices on the minimum covering sphere
C      producing a spherical fullerene
      DIMENSION Dist(3,Nmax),c(3),IATOM(Nmax),IC3(Nmax,3)
      CHARACTER*2 El(99)
      Character*1 TEXTINPUT(nzeile)
      Character*50 filename
      CHARACTER*50 xyzname
      Integer endzeile

      Write(Iout,1000)
      Do I=1,number_vertices
       r2=0.d0
      Do J=1,3
       Dist(J,I)=Dist(J,I)-c(J)
       r2=r2+Dist(J,I)**2
      enddo
       rdist=dsqrt(r2)
       rscale=rmcs/rdist
      Do J=1,3
       Dist(J,I)=Dist(J,I)*rscale
      enddo
      enddo

      Do J=1,number_vertices
       IM=IAtom(J)
       Write(IOUT,1001),J,IM,El(IM),(Dist(I,J),I=1,3)
      enddo

C Write out xyz file
      if(ipsphere.eq.2) then
       xyzname=trim(filename)//"-3DMCS.xyz"
       Open(unit=9,file=xyzname,form='formatted')
        WRITE(Iout,1002) xyzname
        endzeile=0
        do j=1,nzeile
          if(TEXTINPUT(j).ne.' ') endzeile=j
        enddo
        if(number_vertices.lt.100) WRITE(9,1011) number_vertices,
     1    number_vertices,(TEXTINPUT(I),I=1,endzeile)
        if(number_vertices.ge.100.and.number_vertices.lt.1000)
     1    WRITE(9,1012) number_vertices,number_vertices,
     1    (TEXTINPUT(I),I=1,endzeile)
        if(number_vertices.ge.1000.and.number_vertices.lt.10000)
     1    WRITE(9,1013) number_vertices,number_vertices,
     1    (TEXTINPUT(I),I=1,endzeile)
        if(number_vertices.ge.10000)
     1    WRITE(9,1014) number_vertices,number_vertices,
     1    (TEXTINPUT(I),I=1,endzeile)
        Do J=1,number_vertices
          IM=IAtom(J)
          Write(9,1003) El(IM),(Dist(I,J),I=1,3)
        enddo
        Close(unit=9)
      endif

C Write out cc1 file
      if(ipsphere.ge.3) then
       xyzname=trim(filename)//"-3DMCS.cc1"
       Open(unit=9,file=xyzname,form='formatted')
        WRITE(Iout,1002) xyzname
        if(number_vertices.lt.100) WRITE(9,1005) number_vertices
        if(number_vertices.ge.100.and.number_vertices.lt.1000)
     1    WRITE(9,1006) number_vertices
        if(number_vertices.ge.1000.and.number_vertices.lt.10000)
     1    WRITE(9,1007) number_vertices
        if(number_vertices.ge.10000) WRITE(9,1008) number_vertices
        Do J=1,number_vertices
          IM=IAtom(J)
          Write(9,1004) El(IM),J,(Dist(I,J),I=1,3),(IC3(J,I),I=1,3)
        enddo
        Close(unit=9)
      endif

 1000 Format(/1X,'Projecting vertices on minimum covering sphere',
     1 /1x,'Cartesian Input',
     1  /1X,'  I      Z  Element Cartesian Coordinates')
 1001 FORMAT(1X,I4,1X,I6,1X,A2,6X,3(D18.12,2X))
 1002 FORMAT(/1X,'Input coordinates for spherical fullerene to be ',
     1 'used for plotting program CYLVIEW, PYMOL or AVOGADRO',
     1 /1X,'Output written into ',A60)
 1003 FORMAT(A2,6X,3(F15.6,2X))
 1004 FORMAT(A2,I5,3F12.6,'    2',3I5)
 1005 FORMAT(I2)
 1006 FORMAT(I3)
 1007 FORMAT(I4)
 1008 FORMAT(I8)
 1011 FORMAT(I5,/,'C',I2,'/  ',132A1)
 1012 FORMAT(I5,/,'C',I3,'/  ',132A1)
 1013 FORMAT(I5,/,'C',I4,'/  ',132A1)
 1014 FORMAT(I8,/,'C',I8,'/  ',132A1)
      Return
      End

      SUBROUTINE MinCovSphere2(IOUT,imcs,SP,Dist,Rmin,Rmax,
     1 VCS,ACS,Atol,VTol,u,c,radius,RVdWC)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(3,Nmax)
      DIMENSION u(Nmax),c(3),Kappa(Nmax)
      Real*8 meandist
C     Get the minimum covering sphere using algorithm 2 by E.A.Yildirim
      DATA itermax/10000000/
      DATA eps,disteps/1.d-8,1.d-5/
      epsilon=(1.d0+eps)**2-1.d0
      Write(IOUT,1000) number_vertices,epsilon
 1000 Format(/1X,'Calculate the minimum covering sphere using',
     1 ' the algorithm 2 given by E.A.Yildirim,'
     2 ' SIAM Journal on Optimization 19(3), 1368 (2008):',
     3 /1X,' Number of points in 3-D space: ',I4,
     4 ', epsilon for convergence: ',d12.6,/1X,
     5 ' Start search using the center of points and radius to ',
     6 'furthest point from the center as initial guess')
C   
C    Initial step
C 
      xmax=0.d0
      ialpha=0
      do i=1,number_vertices
        u(i)=0.d0
      enddo
C    Algorithm changed here by taking the furthest 
C    point from the center of points
      rav=0.d0
      do i=1,number_vertices
        dx=Dist(1,i)
        dy=Dist(2,i)
        dz=Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
        rav=rav+dsqrt(Xnorm)
          if(Xnorm.gt.xmax) then
            ialpha=i
            xmax=Xnorm
          endif
      enddo
      meandist=rav/dfloat(number_vertices)
      asym=0.d0
C     Fowler asymmetry parameter
      do i=1,number_vertices
        dx=Dist(1,i)
        dy=Dist(2,i)
        dz=Dist(3,i)
        ri=dsqrt(dx*dx+dy*dy+dz*dz)
        asym=asym+(ri-meandist)**2
      enddo
        asym=asym/meandist**2
      xmax=0.d0
      ibeta=0
      do i=1,number_vertices
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
      Kappa(1)=ialpha
      Kappa(2)=ibeta
      do j=1,3
        c(j)=.5d0*(Dist(j,ialpha)+Dist(j,ibeta))
      enddo 
      gammak=fpsi(u,dist)
      xmax=0.d0
      ikappa=0
      do i=1,number_vertices
        dx=c(1)-Dist(1,i)
        dy=c(2)-Dist(2,i)
        dz=c(3)-Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.gt.xmax) then
            ikappa=i
            xmax=Xnorm
          endif
      enddo
      xmin=1.d12
      iksi=0
      do i=1,nkappa
        ik=Kappa(i)
        dx=c(1)-Dist(1,ik)
        dy=c(2)-Dist(2,ik)
        dz=c(3)-Dist(3,ik)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.lt.xmin) then
            iksi=Kappa(i)
            xmin=Xnorm
          endif
      enddo
      dx=c(1)-Dist(1,ikappa)
      dy=c(2)-Dist(2,ikappa)
      dz=c(3)-Dist(3,ikappa)
      deltp=(dx*dx+dy*dy+dz*dz)/gammak-1.d0
      dx=c(1)-Dist(1,iksi)
      dy=c(2)-Dist(2,iksi)
      dz=c(3)-Dist(3,iksi)
      deltm=1.d0-(dx*dx+dy*dy+dz*dz)/gammak
      deltak=dmax1(deltp,deltm)
      k=0
      Write(IOUT,1001) k,nkappa,deltak,(c(i),i=1,3)
C     Now loop
C     ----------------
      do iter=1,itermax
      if(deltak.le.epsilon) Go to 10
      if(deltak.gt.deltm) then
        alamdak=.5d0*deltak/(1.d0+deltak)
        k=k+1
        blambdak=1.d0-alamdak
          do i=1,number_vertices
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
      else
        arg1=.5d0*deltm/(1.d0-deltm)
        arg2=u(iksi)/(1.d0-u(iksi))
        alamdak=arg1
        signt=arg2-arg1
        if(signt.lt.0.d0) alamdak=arg2
          if(signt.lt.0.d0.and.nkappa.gt.1) then
            km=1
            isetk=0
            do i=1,nkappa
            imem=Kappa(i)
            if(imem.ne.iksi) then
            Kappa(km)=imem
            km=km+1
            else
            isetk=isetk-1
            endif
            enddo
            nkappa=nkappa+isetk
          endif
        k=k+1
          do i=1,number_vertices
            u(i)=(1.d0+alamdak)*u(i)
            if(i.eq.iksi) u(i)=u(i)-alamdak
          enddo
          do j=1,3
            c(j)=(1.d0+alamdak)*c(j)-alamdak*dist(j,iksi)
          enddo
      endif
C     New parameters
      gammak=fpsi(u,dist)
      xmax=0.d0
      ikappa=0
      do i=1,number_vertices
        dx=c(1)-Dist(1,i)
        dy=c(2)-Dist(2,i)
        dz=c(3)-Dist(3,i)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.gt.xmax) then
            ikappa=i
            xmax=Xnorm
          endif
      enddo
      xmin=1.d12
      iksi=0
      do i=1,nkappa
        ik=Kappa(i)
        dx=c(1)-Dist(1,ik)
        dy=c(2)-Dist(2,ik)
        dz=c(3)-Dist(3,ik)
        Xnorm=dx*dx+dy*dy+dz*dz
          if(Xnorm.lt.xmin) then
            iksi=ik
            xmin=Xnorm
          endif
      enddo
      dx=c(1)-Dist(1,ikappa)
      dy=c(2)-Dist(2,ikappa)
      dz=c(3)-Dist(3,ikappa)
      deltp=(dx*dx+dy*dy+dz*dz)/gammak-1.d0
      dx=c(1)-Dist(1,iksi)
      dy=c(2)-Dist(2,iksi)
      dz=c(3)-Dist(3,iksi)
      deltm=1.d0-(dx*dx+dy*dy+dz*dz)/gammak
      deltak=dmax1(deltp,deltm)
      enddo
C     ----------------
   10 if(iter.ge.itermax) Write(IOUT,1006) itermax
      Write(IOUT,1005) k,nkappa,deltak,(c(i),i=1,3)
      rmcs=0.d0
      Do i=1,number_vertices
      RI=dsqrt((Dist(1,I)-c(1))**2+(Dist(2,I)-c(2))**2
     1 +(Dist(3,I)-c(3))**2)
      if (ri.gt.rmcs) rmcs=ri
      enddo
      dist0=dsqrt(c(1)*c(1)+c(2)*c(2)+c(3)*c(3))
      Write(IOUT,1002) (c(i),i=1,3),dist0
C     Check if there is no point outside the sphere
C     Finally calculate the surface and volume and compare to previous results
      VMCS=4.d0/3.d0*dpi*RMCS**3
      AMCS=4.d0*dpi*RMCS**2
      if(imcs.ne.0) then
       Write(IOUT,1016) VMCS,AMCS
       Go to 11
      endif
      RatioMCS=AMCS/VMCS
      RatioCS=ACS/VCS
      RatioT=Atol/Vtol
      RatioV=VMCS/Vtol
      AIPQ=36.d0*dpi*Vtol**2/Atol**3
      DIPQ=(1.d0-AIPQ)*100.d0
      Write(IOUT,1004) VMCS,VCS,Vtol,AMCS,ACS,Atol,
     1 RatioMCS,RatioCS,RatioT,RatioV,AIPQ,DIPQ,
     1 SP,SP*1.d2,asym
C     Do statistics
  11  Write(IOUT,1011)
      keq=0
      kgt=0
      klo=0
      rmcsl=-10.d0*epsilon*rmcs
      rmcsu=10.d0*epsilon*rmcs
      Do i=1,number_vertices
      RZ=dsqrt((Dist(1,I)-c(1))**2+(Dist(2,I)-c(2))**2
     1 +(Dist(3,I)-c(3))**2)-rmcs
      u(i)=RZ
      if (rz.lt.rmcsl) klo=klo+1
      if (rz.gt.rmcsu) kgt=kgt+1
      if (rz.le.rmcsu.and.rz.ge.rmcsl) then
      keq=keq+1
      RI=dsqrt((Dist(1,I)-c(1))**2+(Dist(2,I)-c(2))**2
     1 +(Dist(3,I)-c(3))**2)
      if(keq.le.4) Write(IOUT,1003) i,(Dist(j,i),j=1,3),RI
      endif
      enddo
      radius=dsqrt((1.d0+deltak)*gammak)
      Write(IOUT,1008) radius
C     Add the carbon Van der Waals radius
      if(imcs.eq.0) then
       RVdWF=radius+RVdWC
       VVdWF=4.d0/3.d0*dpi*RVdWF**3
       VFCC=VVdWF*3.d0*dsqrt(2.d0)/dpi
       Write(IOUT,1009) RVdWC,RVdWF,VVdWF
       ALC=2.d0*RVdWF*dsqrt(2.d0)
       Write(IOUT,1010) ALC,VFCC,VFCC*0.60221367d0
      endif
C     For each vertex (atom) use the distance to the sphere
C     to calculate the root mean square as a measure for distortion
      rsum=0.d0
      rmax=0.d0
      Do I=1,number_vertices
      rsum=rsum+u(i)
      if(u(i).lt.rmax) rmax=u(i)
      enddo
      rmean=rsum/dfloat(number_vertices)
      distortion=-rmean/rmin*100.d0
      Write(IOUT,1014) 
      Write(IOUT,1013) (i,u(i),i=1,number_vertices)
      Write(IOUT,1007) number_vertices,keq,klo,kgt
      Write(IOUT,1012) rmax,rmean,distortion,number_vertices
      if(distortion.lt.disteps) then
      Write(IOUT,1015)
      endif
 1001 Format(1x,' Cycle ',I4,', nr. points: ',I4,', convergence: ',
     1 D14.8,', Center of Sphere (X,Y,Z): ',3(D12.6,2X))
 1002 Format(1x,' End of iteration',
     1 /1X,' Center of minimum covering sphere (x,y,z):',3(D14.8,2X),
     2 ' at distance ',D14.8,2X,'from origin')
 1003 Format(1x,' P',I4,'= ',3(D14.8,2X),'  RMCS = ',D14.8)
 1004 Format(/1x,' Final result and comparison (in units of input):',
     1 /1x,' Volume of minimum covering sphere: ',D14.8,
     1 /1x,' Volume of central covering sphere: ',D14.8,
     1 /1x,' Volume of fullerene              : ',D14.8,
     1 /1x,' Area of minimum covering sphere  : ',D14.8,
     1 /1x,' Area of central covering sphere  : ',D14.8,
     1 /1x,' Area of fullerene                : ',D14.8,
     1 /1x,' Ratio area/volume of minimum covering sphere: ',D14.8,
     1 /1x,' Ratio area/volume of central covering sphere: ',D14.8,
     1 /1x,' Ratio area/volume of fullerene              : ',D14.8,
     1 /1x,' Ratio V(MCS)/V(cage)                        : ',D14.8,
     1 /1x,' Isoperimetric quotient (q_IP=36*PI*V^2/A^3) : ',D14.8,
     1 /1x,' Deformation parameter from q_IP (in %)      : ',D14.8,
     1 /1x,' Diaz-Tendero sphericity parameter divided by ',
     1     'rotational constant A: ',D14.8,' (',D14.8,'%)',
     1 /1x,' Fowler asymmetry parameter       : ',D14.8)
 1005 Format(1x,' Number of cycles ',I8,', nr. points: ',I4,
     1 ', convergence: ',F14.8,', Center of Sphere (X,Y,Z): ',
     2 3(D12.6,2X))
 1006 Format(1x,' Maximum number of allowed iterations ',I8,' reached')
 1007 Format(/1x,' Final statistics:',/3X,
     1 ' Number of points: ',I4,/3X,
     2 ' Number of points on the sphere: ',I4,/3X,
     3 ' Number of points below the sphere: ',I4,/3X,
     4 ' Number of points above the sphere: ',I4)
 1008 Format(1x,' Radius of minimum covering sphere: ',D14.8)
 1009 Format(/1x,' Van der Waals radius of carbon: ',D14.8,' Angstroem,'
     1 1x,' Van der Waals radius of fullerene: ',D14.8,
     2 ' Angstroem',/1x,' Van der Waals volume: ',D14.8,
     3 ' Angstroem cube ',' (NB: Input coordinates need ',
     4 'to be in Angstroem)')
 1010 Format(1x,' FCC Lattice constant: ',D14.8,' Angstroem,'
     1 ' FCC volume: ',D14.8,' Angstroem cube,'
     2 ' FCC volume: ',D14.8,' cm3/mol')
 1011 Format(/1x,' Points on sphere (minimum 2, at most 4 points for
     1 the algorithm) coordinates and distance:')
 1012 Format(/3x,' Largest vertex distance from sphere: ',D12.6,
     1 /3x,' Mean distance from sphere: ',D12.6, 
     1 /3x,' Distortion (in percent) from the MCS relative to the'
     1 '  minimum bond distance: ',D12.6,' for ',I4,' atoms')
 1013 Format(8(1X,'(',I5,')',1X,D10.4))
 1014 Format(/1x,' Atom and distance from sphere:')
 1015 Format(1x,' Atoms lie on a sphere!')
 1016 Format(/1x,' Final result (in units of input):',
     1 /1x,' Volume of minimum covering sphere: ',D14.8,
     1 /1x,' Area of minimum covering sphere  : ',D14.8)
      return
      END

      DOUBLE PRECISION FUNCTION f1dimx(ncom,IOP,ier,
     1 x,pcom,xicom,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION pcom(ncom),xicom(ncom),xt(ncom)
      REAL*8 Dist(3,Nmax)
      HUGE=dfloat(number_vertices)*1.d2
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      If(IOP.eq.0) then
      Call MDSnorm(ncom,AN,R,xt,Dist)
      f1dimx=AN
      else
      Call MAInorm(ncom,IP,AN,xt,Dist)
      f1dimx=-AN
      endif
      if(AN.gt.Huge) then
       Print*,'**** Severe Error, check coordinates'
       ier=1
       stop
      endif
      return
      END

      FUNCTION fpsi(u,dist)
      use config
C     Calculate the Lagrangian dual for the minimum sphere problem
      IMPLICIT REAL*8 (A-H,O-Z)
      Dimension u(Nmax),dist(3,Nmax),x(3)
      sum1=0.d0
      x(1)=0.d0
      x(2)=0.d0
      x(3)=0.d0
      Do i=1,number_vertices
      sum1=sum1+u(i)*(dist(1,i)**2+dist(2,i)**2+dist(3,i)**2)
      Do j=1,3
      x(j)=x(j)+u(i)*dist(j,i)
      enddo
      enddo
      sum2=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      fpsi=sum1-sum2
      Return
      END

