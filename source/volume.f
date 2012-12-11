      SUBROUTINE VolumeAndArea(points,IDA,exactarea,
     1     exactvolume,hullarea,hullvolume)!,filename)
      use config
      use iso_c_binding
      IMPLICIT REAL*8 (A-H,O-Z)
c      character*50 filename
      dimension points(3,nmax),IDA(nmax,nmax)
      type(c_ptr)::g,p,hull,new_polyhedron,convex_hull,
     1             new_fullerene_graph

      g = new_fullerene_graph(NMax,number_vertices,IDA)
      p = new_polyhedron(g,points)
      exactarea = get_surface_area(p)
      exactvolume      = get_volume(p)

      hull = convex_hull(p)
      hullarea   = get_surface_area(hull)
      hullvolume = get_volume(hull)

C     POV-Ray output. Commented out for now.
c$$$      call draw_polyhedron(p,trim(filename)//"-3D ","pov",
c$$$     1 x'bb9988', x'8899bb',x'99bb88',
c$$$     1 0.4d0, 1.d0, 0.6d0)
c$$$
c$$$      call draw_polyhedron(hull,trim(filename)//"-hull ","pov",
c$$$     1 x'bb9988', x'8899bb',x'99bb88',
c$$$     1 0.2d0, .4d0, 0.6d0)
c$$$

      call delete_polyhedron(p)
      call delete_polyhedron(hull)
      call delete_fullerene_graph(g)
      return
      end 

      SUBROUTINE Volume(Iout,N5MEM,N6MEM,
     1 IDA,N5Ring,N6Ring,DIST,CRing5,CRing6,VolSphere,Asphere,
     2 exactarea,exactvolume,Rmin5,Rmin6,Rmax5,Rmax6)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Calculate the volume (and surface) of the cage molecule by
C     1) Exact polyhedron volume using the normal of faces
C     2) Convex hull 
C     3) Summing over all trigonal pyramids (areas) defined by the vectors
C     CM-CR  (center of cage to the center of ring)
C     CM-CA1 (center of cage to atom 1 in ring)
C     CM-CA2 (center of cage to atom 2 in ring)
C     There are 5 such trigonal pyramids in a 5-ring and 6 in a 6-ring
C     Note that CM is already in the origin
C     Let CR=(X1,Y1,Z1) , CA1=(X2,Y2,Z2) , and CA2(X3,Y3,Z3)
C     Then the volume V for a irregular trigonal pyramid is given by 
C       the determinant
C
C                                 | X1 Y1 Z1 |
C     V = abs(Vdet)  ,   V =  1/6 | X2 Y2 Z2 |
C                                 | X3 Y3 Z3 |
C
      DIMENSION Dist(3,Nmax)
      DIMENSION CRing5(3,Mmax),CRing6(3,Mmax)
      DIMENSION N5MEM(Mmax,5),N6MEM(Mmax,6),IDA(Nmax,Nmax)

C     Sum over all 5-rings
      Vol5=0.d0
      Area5=0.d0
      ITH=0
      DO I=1,N5Ring
      DO J=1,5
      IAT1=N5MEM(I,J)
      J1=J+1
      IF(J1.eq.6) J1=1
      IAT2=N5MEM(I,J1)
      VolTes=VolTP(CRing5(1,I),CRing5(2,I),CRing5(3,I),
     1 Dist(1,IAT1),Dist(2,IAT1),Dist(3,IAT1),
     1 Dist(1,IAT2),Dist(2,IAT2),Dist(3,IAT2))
      AreaTes=AreaTP(CRing5(1,I),CRing5(2,I),CRing5(3,I),
     1 Dist(1,IAT1),Dist(2,IAT1),Dist(3,IAT1),
     1 Dist(1,IAT2),Dist(2,IAT2),Dist(3,IAT2))
      Vol5=Vol5+VolTes
      Area5=Area5+AreaTes
      ITH=ITH+1
      enddo
      enddo
      Vol5=Vol5/6.d0
      Area5=Area5/2.d0
C     Sum over all 6-rings
      Area6=0.d0
      Vol6=0.d0
      If(N6Ring.eq.0) Go to 10
      DO I=1,N6Ring
      DO J=1,6
      IAT1=N6MEM(I,J)
      J1=J+1
      IF(J1.eq.7) J1=1
      IAT2=N6MEM(I,J1)
      VolTes=VolTP(CRing6(1,I),CRing6(2,I),CRing6(3,I),
     1 Dist(1,IAT1),Dist(2,IAT1),Dist(3,IAT1),
     1 Dist(1,IAT2),Dist(2,IAT2),Dist(3,IAT2))
      AreaTes=AreaTP(CRing6(1,I),CRing6(2,I),CRing6(3,I),
     1 Dist(1,IAT1),Dist(2,IAT1),Dist(3,IAT1),
     1 Dist(1,IAT2),Dist(2,IAT2),Dist(3,IAT2))
      Vol6=Vol6+VolTes
      Area6=Area6+AreaTes
      ITH=ITH+1
      enddo
      enddo
      Vol6=Vol6/6.d0
      Area6=Area6/2.d0
C     Now do the total volume calculation
  10  VTol=Vol5+Vol6
      ATol=Area5+Area6
      DVSphere=VTol-VolSphere
      VPC=DVSphere/VolSphere*1.d2
      DVArea=ATol-Asphere
      APC=DVArea/Asphere*1.d2
C     Convex Hull and exact polyhedral Volume and Area
      Call VolumeAndArea(Dist,IDA,exactarea,
     1 exactvolume,hullarea,hullvolume)!,filename)
C     Calculate the volume and area for C60 using R5 and R6
      If(number_vertices.eq.60) then
      R5=(Rmin5+Rmax5)*.5d0
      R5div=dabs(Rmin5-Rmax5)
      R6div1=dabs(R5-Rmin6)
      R6div2=dabs(R5-Rmax6)
      R6=Rmin6
      if(R6div2.gt.R6div1) R6=Rmax6
      sqrt5=dsqrt(5.d0)
      Volico=(5.d0/12.d0)*(3.d0+sqrt5)*(R5+R5+R6)**3
      Volpenta=.5d0*(5.d0+sqrt5)*(R5)**3
      VolIcocap=Volico-Volpenta
      DVIcocap=VTol-VolIcocap
      fac1=dsqrt(sqrt5*10.d0+25.d0)*3.d0
      fac2=5.d0*dsqrt(3.d0)
      r5s=R5**2
      AreaIcocap=fac1*r5s+fac2*(R6**2+r5s+4.d0*R5*R6)
      Write(Iout,1002) VolIco,VolIcocap,AreaIcocap,R5,R6
      If(R5div.gt.1.d-5) Write(Iout,1003) R5div,R5
      endif
      If(number_vertices.eq.20) then
      R5=(Rmin5+Rmax5)*.5d0
      R5div=dabs(Rmin5-Rmax5)
      sqrt5=dsqrt(5.d0)
      Voldode=.25d0*(15.d0+7.d0*sqrt5)*R5**3
      Areadode=dsqrt(sqrt5*10.d0+25.d0)*3.d0*R5**2
      Write(Iout,1005) Voldode,Areadode,R5
      If(R5div.gt.1.d-5) Write(Iout,1003) R5div,R5
      endif
      ratioVhullexact=hullvolume/exactvolume
      ratioAhullexact=hullarea/exactarea
      Write(Iout,1000) VTol,Vol5,Vol6,ITH,DVSphere,VPC
      If(number_vertices.eq.60) Write(Iout,1004) DVIcocap
      Write(Iout,1006) exactvolume,hullvolume,ratioVhullexact
      Write(Iout,1001) ATol,Area5,Area6,DVArea,APC
      Write(Iout,1007) exactarea,hullarea,ratioAhullexact
      Write(Iout,1008) Atol/VTol,exactarea/exactvolume,
     1 hullarea/hullvolume
C    Measure for convexity
      convex1=1.d2*dsqrt(exactarea/hullarea)
      convex2=1d2*(exactvolume/hullvolume)**(1./3.)
      Write(Iout,1009) convex1,convex2
 1000 Format(/1X,'Final result:',/1X,
     1 'Total calculated volume of fullerene through tesselation ',
     1 'from barycenter of polyhedron (valid for convex polyhedra): ',
     1 D14.8,/10X,'(Contributions from 5-ring: ',D14.8,
     1 ' , and from 6-ring: ',D14.8,')'
     1 /3X,' using ',I6,' trigonal pyramidal tesselations',/1X,
     1 ' Deviation from spherical central covering: ',D12.6,1X,
     1 ' in percent: ',D12.6)
 1001 Format(/1X,'Total surface area of fullerene taken ',
     1 ' from triangulation from barycenter of ring: ',D14.8,
     1 /10X,'(Contributions from 5-ring: ',D14.8,
     1 ' , and from 6-ring: ',D14.8,')',/1X,
     1 ' Deviation from spherical central covering: ',D12.6,1X,
     1 ' in percent: ',D12.6)
 1002 Format(/1X,' Volume of the original icosahedron: ',D14.8,
     1 ' and volume and surface area of the capped icosahedron: V= ',
     1 D14.8,1X,', A= ',D14.8,/1X,' Distances taken: R5= ',D12.6,
     1 ', R6= ',D12.6)
 1003 Format(1X,' Warning: Distances in 5-ring not the same, ',
     1 ' difference is ',D12.6,
     1 ', taking the average distance of ',D12.6)
 1004 Format(1X,' Deviation from capped Ih icosahedron: ',D12.6)
 1005 Format(/1X,' Volume and surface area of the dodecahedron: V= ',
     1 D14.8,1X,', A= ',D14.8,/1X,' Distance taken: R5= ',D12.6)
 1006 Format(1X,' Exact volume of fullerene: ',D12.6,
     1      /1X,' Convex hull volume of fullerene: ',D12.6,
     1      /1X,' Ratio V(hull)/V(exact): ',D12.6)
 1007 Format(1X,' Exact surface area of fullerene: ',D12.6,
     1      /1X,' Convex hull surface area of fullerene: ',D12.6,
     1      /1X,' Ratio A(hull)/A(exact): ',D12.6)
 1008 Format(1X,' Ratio Area/Volume: ',D12.6,' (Tesselation), ',
     1 D12.6,' (exact), ',D12.6,' (convex hull)')
 1009 Format(1X,' Percent convexity: ',F12.3,' (area) ',F12.3,
     1 ' (volume)')
      Return
      END 

      FUNCTION VolTP(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)
C     Volume of trigonal pyramid enclosed between three vectors from origin
      IMPLICIT REAL*8 (A-H,O-Z)
      AMat=X1*Y2*Z3+Y1*Z2*X3+Z1*X2*Y3-Z1*Y2*X3-Y1*X2*Z3-X1*Z2*Y3
      VolTP=dabs(AMat)
      Return
      END
 
      FUNCTION AreaTP(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)
C     Area of triangle at end of three vectors from origin
      IMPLICIT REAL*8 (A-H,O-Z)
      C1=Y1*Z2+Y2*Z3+Y3*Z1-Z1*Y2-Z2*Y3-Z3*Y1
      C2=Z1*X2+Z2*X3+Z3*X1-X1*Z2-X2*Z3-X3*Z1
      C3=X1*Y2+X2*Y3+X3*Y1-Y1*X2-Y2*X3-Y3*X1
      ATP=C1*C1+C2*C2+C3*C3
      AreaTP=dsqrt(ATP)
      Return
      END
