      SUBROUTINE Schlegel(NAtom,Nfaces,Nedges,M,msrs,IOUT,IS1,IS2,IS3,
     1 N5M,N6M,N5R,N6R,NRing,Iring,ISchlegel,IC3,IDA,Dist,angle,Rmin,
     1 Tol,CR,CR5,CR6,Symbol)
      use iso_c_binding
C Produce points in 2D-space for Schlegel diagrams using the cone-
C projection method and the perspective projection, or Tutte
C embedding method and optimization to make it distant transitive. 
C The fullerene is rotated first such that the desired point, edge 
C or face is at the top. Euler angles are used for rotation.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 layout2d
      DIMENSION IDA(Natom,Natom)
      DIMENSION layout2d(2,NAtom)
      DIMENSION Dist(3,NAtom),IAtom(NAtom),NRing(Nfaces),distw(3),c(3)
      DIMENSION CR5(3,Nfaces),CR6(3,Nfaces),vec1(3),Iring(Nfaces)
      DIMENSION N5M(Nfaces,5),N6M(Nfaces,6),Rot(3,3),CR(3,Nfaces)
      DIMENSION Rotz(3,3),Symbol(Nfaces),DistS(2,NAtom),RingS(2,Nfaces)
      DIMENSION IC3(natom,3)
      Character*1  Symbol,SRS(msrs,2*msrs),satom,sring,s5ring,s6ring
      Character*12 Symbol1
      type(c_ptr) :: g, new_fullerene_graph, read_fullerene_graph

      Open(unit=2,file='qmga.dat',form='formatted')
      Data epsf,dpi/.12d0,3.14159265358979d0/
C     Parameter set for Program QMGA
      Data DPoint,Dedge/0.5d0,0.1d0/


      satom='o'
      s5ring='^'
      s6ring='*'
      do I=1,msrs
      do J=1,2*msrs
      SRS(I,J)=' '
      enddo
      enddo
      eps=Rmin*epsf
      iorig=0
      Do I=1,3
      Do J=1,3
      Rot(I,J)=0.d0
      enddo
      enddo
      Rot(1,1)=1.d0
      Rot(2,2)=1.d0
      Rot(3,3)=1.d0
      if(ISchlegel.eq.2) then
      WRITE(IOUT,1013) angle
      endif

C Align fullerene for Schlegel diagram
C First check input and determine the coordinates for alignment
C We allow for vertices, edges and ring centers
      If(IS1.NE.0) then
      WRITE(IOUT,1006) IS1,IS2,IS3

C Check of errors in input
      if(is1.lt.0.or.is2.lt.0.or.is3.lt.0) then
       WRITE(IOUT,1007)
       Return
      endif
      if(is2.gt.0.and.is1.eq.is2) then
       WRITE(IOUT,1007)
       Return
      endif
      if(is3.gt.0.and.(is3.eq.is2.or.is3.eq.is1)) then
       WRITE(IOUT,1007)
       Return
      endif

C Vertex
      If(is2.eq.0) then
        Do I=1,3
         c(I)=DIST(I,is1)
        enddo
       WRITE(IOUT,1008) is1,(c(i),i=1,3)
      endif

C Edge
      If(is2.ne.0.and.is3.eq.0) then
       rd=0.d0
        Do I=1,3
         rd=rd+(DIST(I,is1)-DIST(I,is2))**2
         c(I)=.5d0*(DIST(I,is1)+DIST(I,is2))
        enddo
       rdist=dsqrt(rd)
       rtol=rmin*(1.d0+Tol)
        If(rdist.gt.rtol) then
         WRITE(IOUT,1012) rdist,rtol
         Return
        endif
       WRITE(IOUT,1009) is1,is2,(c(i),i=1,3)
      endif

C Ring center
      If(is3.ne.0.and.is2.ne.0.and.is1.ne.0) then

C   Find ring
C   Search in 5-ring
      n6hit=0
      Do I=1,N5R
       n5hit=0
        Do J=1,5
         if(is1.eq.N5M(I,J)) n5hit=n5hit+1
         if(is2.eq.N5M(I,J)) n5hit=n5hit+1
         if(is3.eq.N5M(I,J)) n5hit=n5hit+1
          if(n5hit.eq.3) then
           c(1)=cr5(1,I) 
           c(2)=cr5(2,I) 
           c(3)=cr5(3,I) 
           IR=I
           go to 150
          endif
        enddo
      enddo
C   Search in 6-ring
      Do I=1,N6R
      n6hit=0
       Do J=1,6
        if(is1.eq.N6M(I,J)) n6hit=n6hit+1
        if(is2.eq.N6M(I,J)) n6hit=n6hit+1
        if(is3.eq.N6M(I,J)) n6hit=n6hit+1
        if(n6hit.eq.3) then 
         c(1)=cr6(1,I) 
         c(2)=cr6(2,I) 
         c(3)=cr6(3,I)
         IR=N5R+I 
         go to 150
        endif
       enddo
      enddo
 150  if(n5hit.eq.3) WRITE(IOUT,1010) IR,is1,is2,is3,(c(i),i=1,3)
      if(n6hit.eq.3) WRITE(IOUT,1016) IR,is1,is2,is3,(c(i),i=1,3)
      nhit=n5hit+n6hit
      if(nhit.lt.3) then
       WRITE(IOUT,1017) nhit
       Return
       endif
      endif

C  Now define the vector of the projection axis and rotate molecule
C  such that it coincides with the z-axis
C  **** This still needs some work for better alignement
      R=1.d0/dsqrt(c(1)**2+c(2)**2+c(3)**2)
      Do I=1,3
      distw(i)=c(i)*R
      enddo
      WRITE(IOUT,1014) (distw(i),i=1,3)
C   Construct rotation matrix
      CALL Rotmat(IOUT,Rot,distw)

C   Now rotate all vertices
      DO I=1,M
      distw(1)=Dist(1,I)
      distw(2)=Dist(2,I)
      distw(3)=Dist(3,I)
      CALL Rotate(Rot,distw,vec1)
      Dist(1,I)=vec1(1)
      Dist(2,I)=vec1(2)
      Dist(3,I)=vec1(3)
      enddo

C   Rotate all ring centers
      Do I=1,N5R
      distw(1)=CR5(1,I)
      distw(2)=CR5(2,I)
      distw(3)=CR5(3,I)
      CALL Rotate(Rot,distw,vec1)
      CR5(1,I)=vec1(1)
      CR5(2,I)=vec1(2)
      CR5(3,I)=vec1(3)
      enddo
      Do I=1,N6R
      distw(1)=CR6(1,I)
      distw(2)=CR6(2,I)
      distw(3)=CR6(3,I)
      CALL Rotate(Rot,distw,vec1)
      CR6(1,I)=vec1(1)
      CR6(2,I)=vec1(2)
      CR6(3,I)=vec1(3)
      enddo

      else
       WRITE(IOUT,1011)
       iorig=1
      endif

C   Sort distances according to z-values
      DO I=1,M
      IAtom(I)=I
      enddo
      DO I=1,M
      dMaxx=Dist(3,I)
      iMax=I
      DO K=I+1,M
      IF (dMaxx.lt.Dist(3,K)) then
      iMax=K
      dMaxx=Dist(3,K)
      endif
      enddo

C   Swap
      Do ii=1,3
      iw=IAtom(iMax)
      IAtom(iMax)=IAtom(I)
      IAtom(I)=iw
      distw(ii)=Dist(ii,imax)
      Dist(ii,imax)=Dist(ii,i)
      Dist(ii,i)=distw(ii)
      enddo
      enddo

C   Now do the same with ring centers
      NR=N5R+N6R
      DO I=1,NR
      IRing(I)=NRing(I)
      if(I.LE.N5R) then
      CR(1,I)=CR5(1,I)
      CR(2,I)=CR5(2,I)
      CR(3,I)=CR5(3,I)
      else
      CR(1,I)=CR6(1,I-N5R)
      CR(2,I)=CR6(2,I-N5R)
      CR(3,I)=CR6(3,I-N5R)
      endif
      enddo
      DO I=1,NR
      dMaxx=CR(3,I)
      iMax=I
      DO K=I+1,NR
      dMax1=CR(3,K)
      IF (dMaxx.lt.dmax1) then
      iMax=K
      dMaxx=dMax1
      endif
      enddo

C   Swap 
      ival=IRing(imax)
      IRing(imax)=IRing(I)
      IRing(I)=ival
      Do ii=1,3
      dval=CR(ii,imax)
      CR(ii,imax)=CR(ii,I)
      CR(ii,I)=dval
      enddo
      enddo

C Rotate around the z-axis such that in Schlegel projection
C an edge of the last ring is on the bottom
C Find first the edge and the rotation angle and build the
C rotation matrix
      if(iorig.eq.0) then
      ILast=IRing(NR)
C   Get shortest distance
      if(ILast.le.N5R) then 
      Symbol1='pentagon (P)'
      rminL=1.d10
      Do I=1,5
      Do J=I+1,5
      IA=N5M(ILast,I)
      IB=N5M(ILast,J)
      X=Dist(1,IA)-Dist(1,IB)
      Y=Dist(2,IA)-Dist(2,IB)
      Z=Dist(3,IA)-Dist(3,IB)
      RL=dsqrt(X**2+Y**2+Z**2)
      if(RL.lt.rminL) then
      rminL=RL
      IAL=IA
      IBL=IB
      endif
      enddo
      enddo
      else
      Symbol1='hexagon  (H)'
      rminL=1.d10
      ILast6=ILast-N5R
      Do I=1,6
      Do J=I+1,6
      IA=N6M(ILast6,I)
      IB=N6M(ILast6,J)
      X=Dist(1,IA)-Dist(1,IB)
      Y=Dist(2,IA)-Dist(2,IB)
      Z=Dist(3,IA)-Dist(3,IB)
      RL=dsqrt(X**2+Y**2+Z**2)
      if(RL.lt.rminL) then
      rminL=RL
      IAL=IA
      IBL=IB
      endif
      enddo
      enddo
      endif

C   Determine the rotation angle
C   Coordinates of vector to rotate
C   Get atom label
      Do j=1,M
      if(IAL.eq.IATOM(j)) IAN=J
      if(IBL.eq.IATOM(j)) IBN=J
      enddo
      X=.5d0*(Dist(1,IAN)+Dist(1,IBN))
      Y=.5d0*(Dist(2,IAN)+Dist(2,IBN))
      Z=.5d0*(Dist(3,IAN)+Dist(3,IBN))
      RL=dsqrt(X**2+Y**2)
      If(dabs(Y).gt.1.d-10) then
      anglev=datan(X/Y)
      else
      anglev=.5d0*dpi
      endif
      anglev=anglev+dpi
      anglevd=anglev*180.d0/dpi
      WRITE(IOUT,1020) Symbol1,anglevd,IAL,IBL,RL
C   Build the rotation matrix
      rw1=dcos(anglev)
      rw2=dsin(anglev)
      Rotz(1,1)=rw1
      Rotz(1,2)=-rw2
      Rotz(1,3)=0.d0
      Rotz(2,1)=rw2
      Rotz(2,2)=rw1
      Rotz(2,3)=0.d0
      Rotz(3,1)=0.d0
      Rotz(3,2)=0.d0
      Rotz(3,3)=1.d0
      WRITE(IOUT,1021)
C   Rotate the vertices
      DO I=1,M
      distw(1)=Dist(1,I)
      distw(2)=Dist(2,I)
      distw(3)=Dist(3,I)
      CALL Rotate(Rotz,distw,vec1)
      Dist(1,I)=vec1(1)
      Dist(2,I)=vec1(2)
      Dist(3,I)=vec1(3)
      enddo
C   Rotate the ring centers
      Do I=1,NR
      distw(1)=CR(1,I)
      distw(2)=CR(2,I)
      distw(3)=CR(3,I)
      CALL Rotate(Rotz,distw,vec1)
      CR(1,I)=vec1(1)
      CR(2,I)=vec1(2)
      CR(3,I)=vec1(3)
      enddo
C   Print Complete rotation matrix
      CALL Rotprint(Iout,Rot,Rotz)
      endif
C   Print the sorted vertices
      WRITE(IOUT,1000) eps
      Natomcirc=0
      Ncircle=1
      DO I=1,M
      Natomcirc=Natomcirc+1
      If(I.gt.1) then
      dif=dabs(Dist(3,I)-Dist(3,I-1))
      if(dif.gt.eps) then
      Ncircle=Ncircle+1
      Natomcirc=1
      endif
      endif
      WRITE(IOUT,1001) IAtom(I),(Dist(ii,I),II=1,3),Ncircle,Natomcirc
      enddo
      if(Natomcirc.ne.5.and.Natomcirc.ne.6) then
      WRITE(IOUT,1005)
      endif
C   Print the sorted ring centers
      WRITE(IOUT,1015) eps
      Natomcirc=0
      Ncircle=1
      DO I=1,NR
      Natomcirc=Natomcirc+1
      If(I.gt.1) then
      dif=dabs(CR(3,I)-CR(3,I-1))
      if(dif.gt.eps) then
      Ncircle=Ncircle+1
      Natomcirc=1
      endif
      endif
      If(IRing(I).le.N5R) then
      Symbol(I)='P'
      iboost=5
      else
      Symbol(I)='H'
      iboost=6
      endif
      WRITE(IOUT,1019) IRing(I),(CR(ii,I),II=1,3),
     1 Ncircle,Natomcirc,Symbol(I)
      enddo
      if(Natomcirc.ne.1) then
      WRITE(IOUT,1005)
      endif

C     Prepare for Program QMGA
      Write(2,901) M,DPoint,Dedge

      If(ISchlegel-2) 10,20,30
 
C     Cone projection using the input angle
C     Calculate distance of vertices from z-axis for projection
  20  app=rmin+Dist(3,1)
      WRITE(IOUT,1002)
C     Write out on file unit=2 for qmga
C     See http://qmga.sourceforge.net/
      sfac=dtan(angle*dpi/180.d0)
      Do I=1,M
      X=Dist(1,I)
      Y=Dist(2,I)
      Z=Dist(3,I)
      R=dsqrt(X*X+Y*Y)
      If(R.gt.eps) then
      Fac=(app-Z)/R*sfac
      else
      Fac=1.d0
      If(I.ne.1) Write(IOUT,1003)
      endif
C   Extra boost for the last ring points
      IVert=M-Iboost
      If(I.gt.IVert) then
      Fac=Fac*1.2d0
      endif
      DistS(1,I)=Dist(1,I)*Fac
      DistS(2,I)=Dist(2,I)*Fac

C   Print
      IAT=IAtom(i)
      WRITE(IOUT,1004) IAT,DistS(1,I),DistS(2,I),
     1 IC3(IAT,1),IC3(IAT,2),IC3(IAT,3),Fac
      Write(2,902) IAT,DistS(1,I),DistS(2,I),
     1 IC3(IAT,1),IC3(IAT,2),IC3(IAT,3)
      enddo
      WRITE(IOUT,1032)

C   Calculate distance of ring centers from z-axis for projection
      WRITE(IOUT,1018)
C   Choose point for angle determination in Schlegel diagram
      X0=CR(1,1)
      Y0=CR(2,1)
      XMAX=0.d0
      YMAX=1.d0
      Do I=1,NR
      X=CR(1,I)
      Y=CR(2,I)
      Z=CR(3,I)
      R=dsqrt(X**2+Y**2)
      If(R.gt.eps) then
      Fac=(app-Z)/R*sfac
      else
      Fac=1.d0
      endif
C   Angle between the two vectors cos^-1(a.b/|a||b|)*(180/pi)
      angle=0.d0
      if(I.ne.1) then
      dx=X-X0
      dy=Y-Y0
      dxm=Xmax-X0
      dym=Ymax-Y0
      r1=dsqrt(dx**2+dy**2)
      r2=dsqrt(dxm**2+dym**2)
      anglef=R1*R2
      if(dabs(anglef).gt.1.d-3) then
      farg=(dx*dxm+dy*dym)/anglef
      angle=dacos(farg)*180.d0/dpi
      else
      angle=0.d0
      endif
      if(dx.lt.0.d0) angle=360.d0-angle
      if(dabs(angle-360.d0).lt.1.d-2) angle=0.d0
      endif
      RingS(1,I)=CR(1,I)*Fac
      RingS(2,I)=CR(2,I)*Fac
C   Print
      WRITE(IOUT,1023) IRing(i),Symbol(i),RingS(1,I),RingS(2,I),
     1 Angle,Fac
      enddo
      go to 999

C   Perspective projection using the input distance
C   Algorithm to produce symmetric Schlegel diagrams
C   Smallest and largest ring z-coordinate
 10   Z=CR(3,1)
      scale=.45d0
      zmin=CR(3,NR)
C   Setting point of projection
      if(angle.eq.0.d0) then
      app=scale*rmin+Z
      else
      app=angle+Z
      WRITE(IOUT,1031) angle
      endif
      WRITE(IOUT,1022) app
      Zproj=app-zmin
      WRITE(IOUT,1030)
C   Atoms
      Do I=1,M
      X=Dist(1,I)
      Y=Dist(2,I)
      Z=Dist(3,I)
      Fac=Zproj/(app-Z)
      DistS(1,I)=Dist(1,I)*Fac
      DistS(2,I)=Dist(2,I)*Fac
      IAT=IAtom(i)
      WRITE(IOUT,1028) IAT,DistS(1,I),DistS(2,I),
     1 IC3(IAT,1),IC3(IAT,2),IC3(IAT,3)
      Write(2,902) IAT,DistS(1,I),DistS(2,I),
     1 IC3(IAT,1),IC3(IAT,2),IC3(IAT,3)
      enddo
      WRITE(IOUT,1032)
      WRITE(IOUT,1029)
C   Rings
      Do I=1,NR
      X=CR(1,I)
      Y=CR(2,I)
      Z=CR(3,I)
C     Extra boost for rings
      Fac=1.1d0*Zproj/(app-Z)
      RingS(1,I)=CR(1,I)*Fac
      RingS(2,I)=CR(2,I)*Fac
      WRITE(IOUT,1027) IRing(i),Symbol(i),RingS(1,I),RingS(2,I)
      enddo
      go to 999

C     Use Tutte method and optimize graph
C     J. Avery
C     Implement Tutte method here
C     Note that fullerene is oriented already according to input for Schlegel projection
C     Take top ring for outer rim

C     Algorithm 3 (Tutte):
 30   write (Iout,1033)
      g = new_fullerene_graph(NAtom,MAtom,IDA)
      write (Iout,1034)
      call tutte_layout_b(g,is1,is2,is3,layout2d)
      call delete_fullerene_graph(g)
      write (Iout,1035)


C   Produce Schlegel picture
C   Atoms
 999   span=0.d0
       iflaga=0
      do i=1,M
       xabs=dabs(DistS(1,I))
       yabs=dabs(DistS(2,I))
       if(yabs.gt.xabs) xabs=yabs
       if(xabs.gt.span) span=xabs
      enddo
       msrs2=msrs/2
       grid=dfloat(msrs2-1)
       grids=grid/span
      do i=1,M
      ix=int(DistS(1,I)*grids)
      ixpos=ix+msrs2
      iy=int(DistS(2,I)*grids)
      iypos=iy+msrs2
      if(ixpos.le.0) ixpos=1
      if(iypos.le.0) iypos=1 
      if(ixpos.gt.msrs) ixpos=msrs
      if(iypos.gt.msrs) iypos=msrs
      if(SRS(ixpos,2*iypos).eq.' ') then
      SRS(ixpos,2*iypos)=satom
       else
        iflaga=iflaga+1
       endif
      enddo
C   Rings
      iflagr=0
      do i=1,NR
      ix=int(RingS(1,I)*grids)
      ixpos=ix+msrs2
      iy=int(RingS(2,I)*grids)
      iypos=iy+msrs2
      if(ixpos.le.0) ixpos=1
      if(iypos.le.0) iypos=1 
      if(ixpos.gt.msrs) ixpos=msrs
      if(iypos.gt.msrs) iypos=msrs
      sring=s6ring
      if(IRing(i).le.12) sring=s5ring
      if(SRS(ixpos,2*iypos).eq.' ') then
        SRS(ixpos,2*iypos)=sring
       else
        iflagr=iflagr+1
       endif
      enddo
C   Print Schlegel picture
      WRITE(IOUT,1024) msrs,2*msrs 
      if(iflaga.ne.0.or.iflagr.gt.1) Write(6,1026) iflaga,iflagr
      do I=1,msrs
        WRITE(IOUT,1025) (SRS(I,J),J=1,2*msrs)
      enddo

  901 Format(I6,2F12.6)
  902 Format(I6,2(1X,F12.6),1X,3(1X,I6))
 1000 Format(/1X,'Schlegel diagram for selected projection point',
     1 /1X,' First sort z-distances and produce NC circles'
     2  ' with NA atoms on it: Tolerance= ',F12.6,/1X,
     3  '  Atom       X            Y            Z     NC  NA')
 1001 Format(1X,I4,3(1X,F12.6),2(1X,I3))
 1002 Format(/1X,'Now project the vertices and print adjacent '
     1 'vertices:',/1X,'  Atom       X            Y     '
     1 '  N1   N2   N3   Scale factor')
 1003 Format(1X,'Warning: Fullerene probably misaligned for'
     1 ' Schlegel diagram')
 1004 Format(1X,I4,2(1X,F12.6),1X,3(1X,I4),' / ',F12.6)
 1005 Format(1X,'Warning: Fullerene probably misaligned for'
     1 ' Schlegel diagram (last circumfence should be a'
     2 ' pentagon or hexagon with ring center in the middle)')
 1006 Format(1X,'Aligning Fullerene, Input: ',3I4)
 1007 Format(1X,'ERROR in input')
 1008 Format(1X,'Vertex of atom ',I4,' chosen for Schlegel alignment.',
     1 ' Coordinates of vertex: X= ',F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1009 Format(1X,'Edge between atoms ',I4,' and ',I4,
     1 ' chosen for Schlegel alignment.',
     2 ' Coordinates of edge center: X= ',
     3 F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1010 Format(1X,'5-Ring center 'I4,' defined by atoms ',I4,' , ',
     1 I4,' and ',I4,' chosen for Schlegel alignment.',
     2 /1X,' Coordinates of ring center: X= ',
     3 F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1011 Format(1X,'Original input coordinates chosen'
     1 ' for Schlegel alignment')
 1012 Format(1X,'ERROR: distance for edge ',F12.6,' greater than',
     1 ' tolerance distance of ',F12.6)
 1013 Format(/1X,'Algorithm for producing the Schlegel diagram'
     1 //1X,'Cone Projection using projection angle of ',
     1 d10.4,' degrees')
 1014 Format(1X,'Now do the rotation mapping vector (',
     1 F8.4,',',F8.4,',',F8.4,') into (0,0,1)')
 1015 Format(/1X,'Same for ring centers producing NC circles'
     2  ' with NA centers on it: Tolerance= ',F12.6,/1X,
     3  '  Ring       X            Y            Z     NC  NA  PG',
     4  ' (P: pentagon, H: hexagon)')
 1016 Format(1X,'6-Ring center ',I4,' defined by atoms ',I4,' , ',
     1 I4,' and ',I4,' chosen for Schlegel alignment.',
     2 /1X,' Coordinates of ring center: X= ',
     3 F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1017 Format(1X,'ERROR: No ring found with given atom numbers',
     1 ' (indices found: ',I3,')')
 1018 Format(/1X,'Finally project the ring centers:',/1X,
     2  '  Ring         X            Y      Spiral Angle  Scale factor')
 1019 Format(1X,I4,3(1X,F12.6),2(1X,I3),3X,A1)
 1020 Format(1X,'Last ring is a ',A12,/1X,
     1 'Rotate around the z-axis by ',F12.4,' degrees',/1X,
     2 ' This moves the edge of the last ring between atoms ',
     3 I4,' and ',I4,' with distance ',F12.4,' to the bottom',
     4 ' the Schlegel diagram')
 1021 Format(1X,'')
 1022 Format(/1X,'Perspective Projection using Z= ',
     1 d10.4,' on z-axis')
 1023 Format(1X,I4,1X,A1,2(1X,F12.6),4X,F7.2,' / ',F12.6)
 1024 Format(/1X,'Schlegel diagram of points (o) and rings ',
     1 '(^ for 5 and * for 6) on a (',I3,','I3,') matrix:')
 1025 Format(2X,200A1)
 1026 Format(1X,'Fullerene too large for Schlegel picture, ',
     1 'could not print ',I3,' atom symbols and ',I3,' ring symbols: ',
     1 /2X,'Either increase parameter  msrs  or use proper ',
     1 'plotting ptogram')
 1027 Format(1X,I4,1X,A1,2(1X,F12.6))
 1028 Format(1X,I4,2(1X,F12.6),1X,3(1X,I4))
 1029 Format(/1X,'Finally project the ring centers:',/1X,
     2  '  Ring         X            Y')
 1030 Format(/1X,'Now project the vertices and print adjacent vertices:'
     2 /1X,'  Atom       X            Y       N1   N2   N3')
 1031 Format(/1X,'Reset focal point distance to nearest ring to ',F12.6)
 1032 Format(/1X,'File qmga.dat written out for input into',
     1 ' program QMGA')
 1033 FORMAT(/1X,'Using the Tutte-embedding algorithm to construct ',
     1 'the fullerne',/1X,'Construct the Tutte planar graph and ',
     1 'make it distant transparant')
 1034 Format(1X,'Calculating Tutte-embedding',/1X,
     1 'Coordinates of Tutte graph:')
 1035 FORMAT(1X,'Fullerene graph deleted') 
      Close(unit=2)
      Return
      END

      SUBROUTINE Rotprint(Iout,Rot,Rotz)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Rot(3,3),Rotz(3,3),RotC(3,3)
      DO I =1,3
      DO K =1,3
      RotC(I,K)=0.d0
      DO J =1,3
      RotC(I,K)=RotC(I,K)+Rotz(I,J)*Rot(J,K)
      enddo
      enddo
      enddo
      Write(Iout,1001)
      Write(Iout,1002) ((RotC(I,J),J=1,3),I=1,3)
 1001 Format(1X,'Complete rotation matrix:')
 1002 Format(3(' |',F12.6,1X),' |')
      Return
      END
 
      SUBROUTINE Rotmat(Iout,Rot,vec)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Rot(3,3),vec(3),vec1(3)
      data dpi,eps/3.14159265358979d0,1.d-8/
C     Construct rotation matrix for turning vec into (0,0,1)
C     First rotate around y-axis to produce (0,y1,z1)
      beta=datan(-vec(1)/vec(3))
      facw=180.d0/dpi
      betaw=beta*facw
      cosb=dcos(beta)
      sinb=dsin(beta)
      xi=cosb*vec(1)+sinb*vec(3)
      yi=vec(2)
      zi=cosb*vec(3)-sinb*vec(1)
C     Then  rotate around x-axis to produce (0, 0,z2)
      alpha=datan(yi/zi)
      alphaw=alpha*facw
      cosa=dcos(alpha)
      sina=dsin(alpha)
C     Produce matrix
      Rot(1,1)=cosb
      Rot(1,2)=0.d0
      Rot(1,3)=sinb
      Rot(2,1)=sina*sinb
      Rot(2,2)=cosa
      Rot(2,3)=-sina*cosb
      Rot(3,1)=-cosa*sinb
      Rot(3,2)=sina
      Rot(3,3)=cosa*cosb
C     Now checking on the original vector
      CALL Rotate(Rot,vec,vec1)
      anorm=dsqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
      dif=dabs(anorm-1.d0)
      If(dif.gt.eps) then
      Write(Iout,1000) dif
      endif
      vec1(1)=0.d0
      vec1(2)=0.d0
      if(vec1(3).lt.0.d0) then
      vec1(3)=-vec1(3)
      do i=1,3
      do j=1,3
      Rot(I,J)=-Rot(I,J)
      enddo
      enddo
      endif
      Write(Iout,1001) betaw,alphaw,anorm,vec1(3)
      Write(Iout,1002) 
      Write(Iout,1003) ((Rot(I,J),J=1,3),I=1,3)
 1000 Format(1X,'WARNING: Vector should be unit vector but is ',F12.6)
 1001 Format(1X,'Rotate around y-axis by ',F12.6,' degrees, then'
     1          ' rotate around x-axis by ',F12.6,' degrees.',
     2          ' Norm: ',F12.6,', Z= ',F12.6)
 1002 Format(1X,'Rotation matrix:')
 1003 Format(3(' |',F12.6,1X),' |')
      Return
      END
 
      SUBROUTINE Rotate(Rot,vec,vec1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Rot(3,3),vec(3),vec1(3)
      vec1(1)=Rot(1,1)*vec(1)+Rot(1,2)*vec(2)+Rot(1,3)*vec(3)
      vec1(2)=Rot(2,1)*vec(1)+Rot(2,2)*vec(2)+Rot(2,3)*vec(3)
      vec1(3)=Rot(3,1)*vec(1)+Rot(3,2)*vec(2)+Rot(3,3)*vec(3)
      Return
      END
