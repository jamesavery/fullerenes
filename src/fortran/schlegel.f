      SUBROUTINE Graph2D(IOUT,IS1,IS2,IS3,N5M,N6M,N5R,N6R,NRing,
     1 Iring,ISchlegel,ifs,ndual,nlabel_vertices,IC3,IDA,Mdist,
     1 nhamcyc,Dist,angle,Rmin,Tol,fscale,scalePPG,boost,CR,CR5,
     1 CR6,Symbol,filename)
      use config
      use iso_c_binding
C Produce points in 2D-space for Schlegel diagrams using the cone-
C projection method and the perspective projection, or Tutte
C embedding method and optimization to make it distant transitive. 
C The fullerene is rotated first such that the desired point, edge 
C or face is at the top. Euler angles are used for rotation.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 layout2d
      DIMENSION IDA(Nmax,Nmax)
      DIMENSION layout2d(2,Nmax)
      DIMENSION Dist(3,Nmax),IAtom(Nmax),NRing(Mmax),distw(3),c(3)
      DIMENSION CR5(3,Mmax),CR6(3,Mmax),vec1(3),Iring(Mmax)
      DIMENSION N5M(Mmax,5),N6M(Mmax,6),Rot(3,3),CR(3,Mmax)
      DIMENSION Rotz(3,3),Symbol(Mmax),RingS(2,Mmax)
      DIMENSION IC3(Nmax,3),IS(6)
      Integer MDist(Nmax,Nmax),mhamfield(Nmax)
      Character*1  Symbol,SRS(msrs,2*msrs),satom,sring,s5ring,s6ring
      Character*1  Dummy
      Character*12 Symbol1
      Character*50 filename
      Character*50 hamname
      type(c_ptr) :: g, new_fullerene_graph
      integer layout_is_spherical

      Data epsf/.12d0/
C     Parameter set for Program QMGA
      Data DPoint,Dedge/0.5d0,0.1d0/

C     Prepare for Program QMGA
      if(ifs.ge.2) then 
       Open(unit=2,file=trim(filename)//"-2D.dat",form='formatted')
       Write(2,901) number_vertices,DPoint,Dedge
      endif

C     Construct a graph object from the adjacency matrix
      g = new_fullerene_graph(Nmax,number_vertices,IDA)
     
      dummy=' '
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
      If(ISchlegel.gt.2) Go to 30
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
      n5hit=0
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
 150  nhit=0
      if(n5hit.eq.3) then
        WRITE(IOUT,1010) IR,is1,is2,is3,(c(i),i=1,3)
        nhit=1
      endif
      if(n6hit.eq.3) then
        WRITE(IOUT,1016) IR,is1,is2,is3,(c(i),i=1,3)
        nhit=1
      endif
        if(nhit.eq.0) then
          WRITE(IOUT,1017)
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
      DO I=1,number_vertices
       distw(1)=Dist(1,I)
       distw(2)=Dist(2,I)
       distw(3)=Dist(3,I)
       CALL Rotate(Rot,distw,vec1)
       Dist(1,I)=vec1(1)
       Dist(2,I)=vec1(2)
       Dist(3,I)=vec1(3)
      enddo

C   Rotate all ring centers
C   Pentagons
      Do I=1,N5R
       distw(1)=CR5(1,I)
       distw(2)=CR5(2,I)
       distw(3)=CR5(3,I)
       CALL Rotate(Rot,distw,vec1)
       CR5(1,I)=vec1(1)
       CR5(2,I)=vec1(2)
       CR5(3,I)=vec1(3)
      enddo
C   Hexagons
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
C   Original input chosen
       WRITE(IOUT,1011)
       iorig=1
      endif

C   Sort distances according to z-values
      DO I=1,number_vertices
       IAtom(I)=I
      enddo
      DO I=1,number_vertices
       dMaxx=Dist(3,I)
       iMax=I
      DO K=I+1,number_vertices
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
      IF(dMaxx.lt.dmax1) then
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
      ian=0
      ibn=0
      Do j=1,number_vertices
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
      anglevd=anglev*deg2rad
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
      DO I=1,number_vertices
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
      DO I=1,number_vertices
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
      iboost = 0 ! only to remove compiler warning
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
     1   Ncircle,Natomcirc,Symbol(I)
      enddo
      if(Natomcirc.ne.1) then
        WRITE(IOUT,1005)
      endif

C   Choice between Schlegel projection or Tutte embedding
      If(ISchlegel-2) 10,20,30

C   Algorithm 2: 
C   Cone projection using the input angle
C   Calculate distance of vertices from z-axis for projection
  20  app=rmin+Dist(3,1)
      WRITE(IOUT,1002)
C     Write out on file unit=2 for schlegel.dat
      sfac=dtan(angle*deg2rad)
      Do I=1,number_vertices
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
C     Extra boost for the last ring points
        IVert=number_vertices-Iboost
        If(I.gt.IVert) then
          Fac=Fac*boost
        endif
        IAT=IAtom(I)
        layout2d(1,IAT)=Dist(1,I)*Fac
        layout2d(2,IAT)=Dist(2,I)*Fac

C     Print
        WRITE(IOUT,1004) IAT,layout2d(1,IAT),layout2d(2,IAT),
     1   IC3(IAT,1),IC3(IAT,2),IC3(IAT,3),Fac
        if(ifs.ge.2) then
          Write(2,902) IAT,layout2d(1,IAT),layout2d(2,IAT),
     1     IC3(IAT,1),IC3(IAT,2),IC3(IAT,3)
        endif
      enddo
      if(ifs.ge.2) WRITE(IOUT,1032) trim(filename)//"-2D.dat"

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
      angle=dacos(farg)*rad2deg
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

C   Algorithm 1: 
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
      Do I=1,number_vertices
      X=Dist(1,I)
      Y=Dist(2,I)
      Z=Dist(3,I)
      Fac=Zproj/(app-Z)
      IAT=IAtom(I)
      layout2d(1,IAT)=Dist(1,I)*Fac
      layout2d(2,IAT)=Dist(2,I)*Fac
      WRITE(IOUT,1028) IAT,layout2d(1,IAT),layout2d(2,IAT),
     1 IC3(IAT,1),IC3(IAT,2),IC3(IAT,3)
      if(ifs.ge.2) 
     1 Write(2,902) IAT,layout2d(1,IAT),layout2d(2,IAT),
     1 IC3(IAT,1),IC3(IAT,2),IC3(IAT,3)
      enddo
      if(ifs.ge.2) WRITE(IOUT,1032) trim(filename)//"-2D.dat"
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

C   Produce Schlegel picture
C   Atoms
 999   span=0.d0
       iflaga=0
      do i=1,number_vertices
       xabs=dabs(layout2d(1,I))
       yabs=dabs(layout2d(2,I))
       if(yabs.gt.xabs) xabs=yabs
       if(xabs.gt.span) span=xabs
      enddo
       msrs2=msrs/2
       grid=dfloat(msrs2-1)
       grids=grid/span
      do i=1,number_vertices
      ix=int(layout2d(1,I)*grids)
      ixpos=ix+msrs2
      iy=int(layout2d(2,I)*grids)
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
      go to 9999

C   Algorithm 3 (Tutte):
C   Use Tutte method and optimize graph
C   Take top ring for outer rim
 30   write (Iout,1033)
      if(is1.le.0.or.is2.le.0.or.is3.le.0) then
         Do I=1,5
            IS(I)=N5M(1,I)
         enddo
         IS(6)=0
         write (Iout,1034) (IS(I),I=1,5)
         lring=5
      else
C   Find ring
C   Search in 5-ring
      n5hit=0 ! only to remove warning
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
           go to 250
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
         go to 250
        endif
       enddo
      enddo
 250  nhit=0
      if(n5hit.eq.3) then
      nhit=1
       Do I=1,5
        IS(I)=N5M(IR,I)
       enddo
       lring=5
       IS(6)=0
       WRITE(IOUT,1038) IR,(IS(I),I=1,5),(c(i),i=1,3)
      endif
      if(n6hit.eq.3) then
      nhit=1
       Do I=1,6
        IS(I)=N6M(IR,I)
       enddo
       lring=6
       WRITE(IOUT,1039) IR,(IS(I),I=1,6),(c(i),i=1,3)
      endif
      if(nhit.eq.0) then
       WRITE(IOUT,1017)
       Return
       endif
C  End of search
      endif

      if(ISchlegel.ge.7) then
       WRITE(IOUT,1041) 
       maxl=0
      Do I=1,number_vertices
      Do J=I+1,number_vertices
       if(MDist(I,J).gt.maxl) maxl=MDist(I,J)
      enddo
      enddo
      endif
      write (Iout,1036)
      call tutte_layout_b(g,is(1),is(2),is(3),layout2d)
C     Get Barycenter of outer ring and use as origin
      write (Iout,1040)
      xc=0.d0
      yc=0.d0
      do I=1,lring
       II=IS(I)
       xc=xc+layout2d(1,II)
       yc=yc+layout2d(2,II)
      enddo
      denom=dfloat(lring)
      xc=xc/denom
      yc=yc/denom
      Do I=1,number_vertices
       layout2d(1,I)=layout2d(1,I)-xc
       layout2d(2,I)=layout2d(2,I)-yc
      enddo
      if(ISchlegel.le.4) then
C  Radially scale Tutte graph
       if(ISchlegel.eq.4) then
c        CALL ScaleTutteB(g,number_vertices,Iout,IS,lring,layout2d)
         CALL ScaleTutte(number_vertices,Iout,IS,lring,fscale,layout2d)
       endif
C     Write to unit 2
      write (Iout,1037)
       Do I=1,number_vertices
       WRITE(IOUT,1028) I,layout2d(1,I),layout2d(2,I),
     1  IC3(I,1),IC3(I,2),IC3(I,3)
      if(ifs.eq.2.or.ifs.eq.3) 
     1  Write(2,902) I,layout2d(1,I),layout2d(2,I),
     1  IC3(I,1),IC3(I,2),IC3(I,3)
       enddo
       Close(unit=2)
C  Optimize Tutte graph, spring embedding
      else
C  IOP=1: Spring embedding
C  IOP=2: Spring embedding + Coulomb repulsion from barycenter
C  IOP=3: Pisanski-Plestenjak-Graovac embedding using the distance matrix MDist
C  IOP=4: Kamada-Kawai embedding using the distance matrix MDist
       IOP=ISchlegel-4
       CALL OptGraph(IOP,Iout,IDA,IS,IC3,MDist,maxl,
     1  scalePPG,layout2d)
       write (Iout,1037)
       Do I=1,number_vertices
       WRITE(IOUT,1028) I,layout2d(1,I),layout2d(2,I),
     1  IC3(I,1),IC3(I,2),IC3(I,3)
      if(ifs.ge.2) Write(2,902) I,layout2d(1,I),layout2d(2,I),
     1  IC3(I,1),IC3(I,2),IC3(I,3)
       enddo
      endif

 9999 layout_is_spherical = 0
      call set_layout2d(g, layout2d, layout_is_spherical)
      if(ifs.eq.1.or.ifs.eq.3) then

      if(nhamcyc.ne.0) then
       hamname=trim(filename)//".ham"
       Open(unit=8,file=hamname,form='formatted')
       Read(8,*) numberham
       if(numberham.ne.number_vertices) then
        Write(Iout,1035) numberham,number_vertices
        Return
       endif
        Do I=1,nhamcyc-1
         Read(8,'(A1)') Dummy
        enddo
        Read(8,*) (mhamfield(i),I=1,numberham)
        Write(Iout,1042) nhamcyc,(mhamfield(i),I=1,numberham)
        close(unit=8)

C     Call format: draw_graph_with_path(filename, format (string), dimensions ((w,h) in cm), 
C     edge_colour (x'rrggbb'), path_colour (x'rrggbb), vertex_colour (x'rrggbb), 
C     edge_width (in mm), path_width (in mm), vertex_diameter (in mm) )
            call draw_graph_with_path(g,filename,"tex",(/15.d0,15.d0/),
     1           x'270470',x'bb2510',x'458b00', 0.3d0, 1.0d0, 2.5d0, 
     2           numberham,mhamfield)
      else
C     Call format: draw_graph(filename, format (string),ndual (0|1), dimensions ((w,h) in cm), 
C     edge_colour (x'rrggbb'), vertex_colour (x'rrggbb), 
C     edge_width (in mm), vertex_diameter (in mm) )
         call draw_graph(g,filename,"tex",ndual, nlabel_vertices,
     1     (/15.d0,15.d0/), x'274070', x'458b00', 0.3d0, 2.5d0)
      endif
c Output to POVRay raytracer. Commented out currently.
c$$$      call draw_graph(g,filename, "pov",0, (/10.d0,10.d0/), 
c$$$     1                x'bb7755', x'8899bb', 0.3d0, .8d0)
      endif
      call delete_fullerene_graph(g)

      if(ifs.ge.2) Close(unit=2)
      Return
  901 Format(I6,2F12.6)
  902 Format(I6,2(1X,F12.6),1X,3(1X,I6))
 1000 Format(/1X,'Schlegel diagram for selected projection point',
     1 /1X,' First sort z-distances and produce NC circles'
     2  ' with NA atoms on it: Tolerance= ',F12.6,/1X,
     3  '  Atom       X            Y            Z      NC   NA')
 1001 Format(1X,I4,3(1X,F12.6),2(1X,I4))
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
 1008 Format(1X,'Vertex of atom ',I4,' chosen for alignment.',
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
     3  '  Ring       X            Y            Z      NC   NA  PG',
     4  ' (P: pentagon, H: hexagon)')
 1016 Format(1X,'6-Ring center ',I4,' defined by atoms ',I4,' , ',
     1 I4,' and ',I4,' chosen for Schlegel alignment.',
     2 /1X,' Coordinates of ring center: X= ',
     3 F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1017 Format(1X,'ERROR: No ring found with given atom numbers')
 1018 Format(/1X,'Finally project the ring centers:',/1X,
     2  '  Ring         X            Y      Spiral Angle  Scale factor')
 1019 Format(1X,I4,3(1X,F12.6),2(1X,I4),3X,A1)
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
     1 '(^ for 5 and * for 6) on a (',I4,','I4,') matrix:')
 1025 Format(2X,200A1)
 1026 Format(1X,'Fullerene too large for Schlegel picture, ',
     1 'could not print ',I4,' atom symbols and ',I4,' ring symbols: ',
     1 /2X,'Either increase parameter  msrs  or use proper ',
     1 'plotting program')
 1027 Format(1X,I4,1X,A1,2(1X,F12.6))
 1028 Format(1X,I4,2(1X,F12.6),1X,3(1X,I4))
 1029 Format(/1X,'Finally project the ring centers:',/1X,
     2  '  Ring         X            Y')
 1030 Format(/1X,'Now project the vertices and print adjacent vertices:'
     2 /1X,'  Atom       X            Y       N1   N2   N3')
 1031 Format(/1X,'Reset focal point distance to nearest ring to ',F12.6)
 1032 Format(/1X,'File ',A20,' written out for input into',
     1 ' 2D graph program')
 1033 FORMAT(/1X,'Using the Tutte-embedding algorithm for fullerene ',
     1 'graph as a starting point for embedding algorithms')
 1034 Format(1X,'No input for specifying circumfencing ring chosen, ',
     1 'first pentagon taken instead including atoms ',5I4)
 1035 Format(1X,'Number of vertices in Hamiltonian file ',I5,
     1 ' not identical to number of vertices ',I5,' ===> RETURN')
 1036 Format(1X,'Calculating Tutte-embedding and shift to barycenter')
 1037 FORMAT(1X,'Fullerene graph deleted',/,' Graph coordinates:',
     1 /,'  Atom       X            Y        N1   N2   N3') 
 1038 Format(1X,'5-Ring center 'I4,' defined by atoms ',3(I4,' , '),
     1 I4,' and ',I4,' chosen for Schlegel alignment.',
     2 /1X,' Coordinates of ring center: X= ',
     3 F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1039 Format(1X,'6-Ring center ',I4,' defined by atoms ',4(I4,' , '),
     1 I4,' and ',I4,' chosen for Schlegel alignment.',
     2 /1X,' Coordinates of ring center: X= ',
     3 F12.6,', Y= ',F12.6,', Z= ',F12.6)
 1040 Format(1X,'Set Tutte graph to barycenter of outer ring')
 1041 Format(1X,'Distance matrix produced')
 1042 Format(/1X,'Hamilton cycles read in from external file:',
     1 /1X,'Cycle number: ',I5,/1X,'Vertex numbers: ',500(I3,'-'))
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
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Rot(3,3),vec(3),vec1(3)
      data eps/1.d-8/
C     Construct rotation matrix for turning vec into (0,0,1)
C     First rotate around y-axis to produce (0,y1,z1)
      beta=datan(-vec(1)/vec(3))
      betaw=beta*rad2deg
      cosb=dcos(beta)
      sinb=dsin(beta)
      xi=cosb*vec(1)+sinb*vec(3)
      yi=vec(2)
      zi=cosb*vec(3)-sinb*vec(1)
C     Then  rotate around x-axis to produce (0, 0,z2)
      alpha=datan(yi/zi)
      alphaw=alpha*rad2deg
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

      SUBROUTINE func2d(IOP,A,IS,MDist,maxd,p,fc,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Embedding algorithms for fullerene graph, energy
      Real*8 p(NMAX*2)
      Integer A(NMAX,NMAX),IS(6),MDist(NMAX,NMAX)
      Data r,f,coulomb/2.0d0,1.d-1,.3d0/
      fc=0.d0
C     simple spring embedding
      if(IOP.le.2) then
      ehook=0.d0
      Do I=1,2*number_vertices,2
        I1=(I+1)/2
      Do J=I+2,2*number_vertices,2
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
      Do I=1,2*number_vertices,2
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
      Do I=1,2*number_vertices,2
        I1=(I+1)/2
        maxu=1000000
        do K=1,iloop
         K1=IS(K)
         if(Mdist(K1,I1).lt.maxu) maxu=Mdist(K1,I1)
        enddo
      Do J=I+2,2*number_vertices,2
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
      Do I=1,2*number_vertices,2
        I1=(I+1)/2
      Do J=I+2,2*number_vertices,2
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

      SUBROUTINE dfunc2d(IOP,A,IS,MDist,maxd,p,x,RAA)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Embedding algorithms for fullerene graph, gradient
      Real*8 p(NMAX*2),x(NMAX*2)
      Integer A(NMAX,NMAX),IS(6),MDist(NMAX,NMAX)
      Data r,f,coulomb/2.0d0,1.d-1,.3d0/
C     simple spring embedding
      if(IOP.le.2) then
      Do I=1,2*number_vertices,2
        ehookx=0.d0
        ehooky=0.d0
        I1=(I+1)/2
      Do J=1,2*number_vertices,2
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
      Do I=1,2*number_vertices,2
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
      Do I=1,2*number_vertices,2
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
       Do J=1,2*number_vertices,2
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
      Do I=1,2*number_vertices,2
        ehookx=0.d0
        ehooky=0.d0
        I1=(I+1)/2
      Do J=1,2*number_vertices,2
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

      SUBROUTINE ScaleTutte(M,Iout,IS,lring,fs,Dist)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Dist(2,Nmax),IS(6)
C Routine to move inner vertices of a Tutte graph outwards
C through a linear relation
C   Determine min distance of outer ring
      Write(Iout,1000) fs
      rmin=1.d10
      do i=1,lring
       II=IS(I)
       r=dsqrt(Dist(1,II)**2+Dist(2,II)**2)
       if(r.lt.rmin) rmin=r
      enddo
C   Scale
      do i=1,M
       r=dsqrt(Dist(1,I)**2+Dist(2,I)**2)
       scale=1.d0+.5d0*fs*(rmin-r)/rmin      
       Dist(1,I)=Dist(1,I)*scale
       Dist(2,I)=Dist(2,I)*scale
      enddo
 1000 Format(1X,'Linear scaling of Tutte graph to move inner ',
     1 'vertices out (Pre-factor for scaling: ',F12.2,')')
      return
      END

c     Layout based on Tutte embedding: Keep angles
c     to barycentrum, throw away radii and distribute
c     evenly.
      SUBROUTINE ScaleTutteB(g,M,Iout,IS,lring,Dist)
      use config
      use iso_c_binding
      integer is(6), depths(M), d_max,i
      real*8 Dist(2,Nmax), c(2), radius, dr, d
      type(c_ptr)::g
      Write(Iout,1000) 

      call vertex_depth(g,IS,lring,depths,d_max)

C     Calculate barycentre 
      c = (/0,0/)
      do i=1,M
         c(1) = c(1) + Dist(1,i)*(1.d0/M)
         c(2) = c(2) + Dist(2,i)*(1.d0/M)
      end do

C   Scale
      do i=1,M
         d = sqrt((Dist(1,i)-c(1))**2 + (Dist(2,i)-c(2))**2)
         dr = .8/REAL(d_max-.5d0) ! Slightly less than radius 0.809 of circle inscribed in pentagon

         if(d.gt.0.d0) then
            if(depths(i).eq.0) then
               radius = 1.0
            else
               radius = dr*(d_max-depths(i)+.5d0)
            endif
            write (*,*) "Scaling vertex", i, "at depth",depths(i),
     1           "by",radius
            Dist(1,I)=(Dist(1,I)-c(1))*radius/d+c(1)
            Dist(2,I)=(Dist(2,I)-c(2))*radius/d+c(2)
         endif
      enddo
 1000 Format(1X,'Graph layout from Tutte-embedding angles, linearly 
     1           distributed raii.')
      return
      END
