      SUBROUTINE func3d(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt)
      use config
c n=MATOM*3
      IMPLICIT REAL*8 (A-H,O-Z)
      integer iopt

      select case(iopt)
        case(1)
          CALL wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt)
        case(2)
          CALL wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt)
        case(3)
          CALL extwu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force)
      end select

      return
      END


c subroutine dist takes 6 reals (=2 coordinates) and yields a positive distance
      SUBROUTINE DIST(ax,ay,az,bx,by,bz,dist_ab)
      real*8 ax,ay,az,bx,by,bz,dist_ab
      dist_ab=dsqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)
      return
      END  


c subroutine angle takes 9 reals (=3 coordinates) and yields an angel between 0 and +\pi (in radians)
      SUBROUTINE ANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,
     2 r1L,r2L,r2M,r1R,r2R,angle_abc
      r2L=(ax-bx)**2 + (ay-by)**2 + (az-bz)**2
      r1L=dsqrt(r2L)
      r2M=(ax-cx)**2 + (ay-cy)**2 + (az-cz)**2
      r2R=(bx-cx)**2 + (by-cy)**2 + (bz-cz)**2
      r1R=dsqrt(r2R)
      angle_abc=dacos((r2L+r2R-r2M)/(2.0*r1L*r1R))
      return
      END  


c subroutine dist takes 12 reals (=4 coordinates) and yields an angel between -\pi/2 and +\pi/2 (in radians)
      SUBROUTINE DIHEDRAL(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2 dihedral_abcd)
      use config
      IMPLICIT REAL*8 (a-z)

c normal vectors on abc and bcd
      abc_x=-az*by+ay*bz+az*cy-bz*cy-ay*cz+by*cz
      abc_y= az*bx-ax*bz-az*cx+bz*cx+ax*cz-bx*cz
      abc_z=-ay*bx+ax*by+ay*cx-by*cx-ax*cy+bx*cy
      bcd_x=-bz*cy+by*cz+bz*dy-cz*dy-by*dz+cy*dz
      bcd_y= bz*cx-bx*cz-bz*dx+cz*dx+bx*dz-cx*dz
      bcd_z=-by*cx+bx*cy+by*dx-cy*dx-bx*dy+cx*dy
c their respective lengths
      abc_length=dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      bcd_length=dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
c normal vectors (length 1) on abc and bcd
      abc_1_x=abc_x/abc_length
      abc_1_y=abc_y/abc_length
      abc_1_z=abc_z/abc_length
      bcd_1_x=bcd_x/bcd_length
      bcd_1_y=bcd_y/bcd_length
      bcd_1_z=bcd_z/bcd_length
c      write(*,*)dsqrt(abc_1_x**2 + abc_1_y**2 + abc_1_z**2)
c      write(*,*)dsqrt(bcd_1_x**2 + bcd_1_y**2 + bcd_1_z**2)
c two auxiliary vectors
      bc_x=bx-cx
      bc_y=by-cy
      bc_z=bz-cz
      bc_length=dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      bc_1_x=bc_x/bc_length
      bc_1_y=bc_y/bc_length
      bc_1_z=bc_z/bc_length
      aux_x=abc_1_y*bc_1_z-bc_1_y*abc_1_z
      aux_y=abc_1_z*bc_1_x-bc_1_z*abc_1_x
      aux_z=abc_1_x*bc_1_y-bc_1_x*abc_1_y
c two auxiliary reals
      aux_1=abc_1_x*bcd_1_x + abc_1_y*bcd_1_y + abc_1_z*bcd_1_z
      aux_2=aux_x*bcd_1_x + aux_y*bcd_1_y + aux_z*bcd_1_z
c the result
      dihedral_abcd=atan2(aux_2, aux_1)
      write(*,*)abc_length,bcd_length,aux_1,aux_2,dihedral_abcd
      return
      END  


      SUBROUTINE wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)

C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, energy
      Real*8 p(nmax*3),force(ffmaxdim)
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6)
      IERR=0
      rp=force(1)
      rh=force(2)
      ap=force(3)
      ah=force(4)
      frp=force(5)
      frh=force(6)
      fap=force(7)
      fah=force(8)
      fco=force(9)

C     Stretching
      ehookrp=0.d0
      ehookrh=0.d0
      Do I=1,n,3
        I1=(I+2)/3
        Do J=I+3,n,3
          J1=(J+2)/3
          if(A(I1,J1).ne.0) then
            call dist(p(i),p(i+1),p(i+2),p(j),p(j+1),p(j+2),ratom)
C           Check if bond is part of 5-ring
            do IB=1,12
              ir1=0
              ir2=0
              do JB=1,5
                if(I1.eq.N5M(IB,JB)) ir1=1
                if(J1.eq.N5M(IB,JB)) ir2=1
              enddo
              if(ir1.eq.1.and.ir2.eq.1) then
C               5-ring
                ehookrp=ehookrp+(ratom-rp)**2
                go to 1
              endif
            enddo
C           6-ring
            ehookrh=ehookrh+(ratom-rh)**2
          endif
  1       continue
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
        call angle(p(JL),p(JL+1),p(JL+2),p(JM),p(JM+1),p(JM+2),
     2   p(JR),p(JR+1),p(JR+2),angle_p)
        ehookap=ehookap+(angle_p-ap)**2
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
        call angle(p(JL),p(JL+1),p(JL+2),p(JM),p(JM+1),p(JM+2),
     2   p(JR),p(JR+1),p(JR+2),angle_h)
        ehookah=ehookah+(angle_h-ah)**2
      enddo
      enddo

C     Coulomb repulsion from origin
      ecoulomb=0.d0
      if (iopt.eq.2 .and. fco.ne.0.d0)  then
       Do I=1,n,3
        rinv=1.d0/dsqrt(p(I)**2+p(I+1)**2+p(I+2)**2)
        ecoulomb=ecoulomb+rinv
       enddo
      endif

C     total energy  
  2   fc=frp*ehookrp+frh*ehookrh+fap*ehookap+fah*ehookah+fco*ecoulomb
      Return
      END



      SUBROUTINE extwu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 p(nmax*3),force(ffmaxdim),increment
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),neighbour_atoms(3),
     2 neighbour_faces_h(3),neighbour_faces_p(3),pentagoncount,
     3 hexagoncount,arbitrary_index
      IERR=0
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

C Stretching
c we distinguish between bonds between two hexagons, two pentagons and hex/pent
      ehookrhh=0.d0
      ehookrhp=0.d0
      ehookrpp=0.d0
      Do I=1,n,3 ! n = number of atoms * 3 !!!
        I1=(I+2)/3 ! I1 = 1, 2, ... (n+2)/3
        Do J=I+3,n,3
          J1=(J+2)/3 ! J1 = I1, I1+1 ... (n+2)/3
          if(A(I1,J1).ne.0) then ! if connected
c get distance
            call dist(p(i),p(i+1),p(i+2),p(j),p(j+1),p(j+2),ratom)
C Check if bond is part of 5-ring
            pentagoncount=0
            do IB=1,12! number of pentagons
              ir1=0
              ir2=0
              do JB=1,5!number of atoms per pentagons
                if(I1.eq.N5M(IB,JB)) ir1=1 !
                if(J1.eq.N5M(IB,JB)) ir2=1 ! if I1 and J2 happen to be in the same pentagon
              enddo
              if(ir1.eq.1 .and. ir2.eq.1) then
                pentagoncount=pentagoncount+1
              endif
            enddo
            if(pentagoncount.eq.0) then
C             6-ring, 6-ring
              ehookrhh=ehookrhh+(ratom-rhh)**2
            else if(pentagoncount.eq.1) then
C             5-ring, 6-ring
              ehookrhp=ehookrhp+(ratom-rhp)**2
            else
C             5-ring, 5-ring
              ehookrpp=ehookrpp+(ratom-rpp)**2
            endif
          endif ! connected
        enddo
      enddo


C Bending
c we distinguish between angles of pentagons and hexagons
C Loop over 5-rings
      ehookap=0.d0
      Do I=1,N5 ! and N5=12
        Do J=1,5
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=5 ! pseudo cyclic sequence
          if(JRX.eq.6) JRX=1
          JM=3*N5M(I,J)-2 ! middle
          JL=3*N5M(I,JLX)-2 ! left
          JR=3*N5M(I,JRX)-2 ! right
          call angle(p(JL),p(JL+1),p(JL+2),p(JM),p(JM+1),p(JM+2),
     2     p(JR),p(JR+1),p(JR+2),angle_p)
          ehookap=ehookap+(angle_p-ap)**2
        enddo
      enddo
C     Loop over 6-rings
      ehookah=0.d0
      Do I=1,N6
        Do J=1,6
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=6
          if(JRX.eq.7) JRX=1
          JM=3*N6M(I,J)  -2
          JL=3*N6M(I,JLX)-2
          JR=3*N6M(I,JRX)-2
          call angle(p(JL),p(JL+1),p(JL+2),p(JM),p(JM+1),p(JM+2),
     2     p(JR),p(JR+1),p(JR+2),angle_h)
          ehookah=ehookah+(angle_h-ah)**2
        enddo
      enddo


C dihedrals 
      ehookdppp=0.d0
      ehookdhpp=0.d0
      ehookdhhp=0.d0
      ehookdhhh=0.d0
      Do I=1,n/3 ! iterate over atoms
c count adjacent pentagons (0 to 3)
        pentagoncount=0
        do IB=1,12 ! iterate over pentagons
          do JB=1,5 ! iterate over atoms in pentagon
            if(I.eq.N5M(IB,JB)) then
              pentagoncount=pentagoncount+1 ! find if I is part of 0,...,3 pentagons
              neighbour_faces_p(pentagoncount)=IB
            endif
          enddo
        enddo
c count adjacent hexagons (0 to 3)
        hexagoncount=0
        do IB=1,n/3-12 ! because n=vertex_count * 3
          do JB=1,6
            if(i.eq.N6M(IB,JB)) then
              hexagoncount=hexagoncount+1
              neighbour_faces_h(hexagoncount)=IB
            endif
          enddo
        enddo
c find neighbouring atoms (3)
        neighbour_atom_count=1
        do j=1,n/3
          if(A(I,j).ne.0) then
            neighbour_atoms(neighbour_atom_count)=j
            neighbour_atom_count=neighbour_atom_count+1
          endif
        enddo
c sort neighbours
c we make use of the fact, that for any dihedral ijkl=-ikjl (hence: ijkl=abs(ikjl))
c therefore we only need to find the special vertex and dont care about the others which are much harder to distinguish
        if(pentagoncount.eq.1) then
          do k=1,3 ! iterate over neighbour atoms
            arbitrary_index=0
            do l=1,hexagoncount ! iterate over neighbour hexagons
              do m=1,6 ! iterate over atoms in hexagon
                if (neighbour_atoms(k).eq.N6M(l,m)) then
                  arbitrary_index=arbitrary_index+1
                endif
              enddo
            enddo
            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two hexagons
              buffer=neighbour_atoms(k)
              neighbour_atoms(k)=neighbour_atoms(1)
              neighbour_atoms(1)=buffer
            endif
          enddo
        endif
        if(pentagoncount.eq.2) then
          do k=1,3 ! iterate over neighbour atoms
            arbitrary_index=0
            do l=1,pentagoncount! iterate over neighbour pentagons
              do m=1,5 ! iterate over atoms in hexagon
                if (neighbour_atoms(k).eq.N5M(l,m)) then
                  arbitrary_index=arbitrary_index+1
                endif
              enddo
            enddo
            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two pentagons
              buffer=neighbour_atoms(k)
              neighbour_atoms(k)=neighbour_atoms(1)
              neighbour_atoms(1)=buffer
            endif
          enddo
        endif
c atoms
        J1=neighbour_atoms(1)
        J2=neighbour_atoms(2)
        J3=neighbour_atoms(3)
        J4=I
c        write(*,*)j1,j2,j3,j4,"atoms of dihedral"
c coordinates
        call dihedral(p(J1*3-2),p(J1*3-1),p(J1*3),p(J2*3-2),p(J2*3-1),
     2   p(J2*3),p(J3*3-2),p(J3*3-1),p(J3*3),p(J4*3-2),p(J4*3-1),
     3   p(J4*3),angle_abcd)
c        write(*,*)i,angle_abcd,"dihedral angle (in radians)"
        select case(pentagoncount)
          case(0)
            zero_value=dhhh
          case(1)
            zero_value=dhhp
          case(2)
            zero_value=dhpp
          case(3)
            zero_value=dppp
        end select
        angle_abcd=dabs(angle_abcd)
        increment=(angle_abcd-zero_value)**2
        write(*,*)i,angle_abcd*57.29,zero_value*57.29,"dihedral angle"
        select case(pentagoncount)
          case(0)
            ehookdhhh=ehookdhhh+increment
          case(1)
            ehookdhhp=ehookdhhp+increment
          case(2)
            ehookdhpp=ehookdhpp+increment
          case(3)
            ehookdppp=ehookdppp+increment
        end select
      enddo


C     total energy  
      fc=frpp*ehookrpp+frhp*ehookrhp+frhh*ehookrhh ! stretching
     2 +fap*ehookap+fah*ehookah ! bending
     3 +fdppp*ehookdppp+fdhpp*ehookdhpp+fdhhp*ehookdhhp+fdhhh*ehookdhhh! dihedral
c      write(*,*)fc,"energy"
      Return
      END



      SUBROUTINE dfunc3d(n,A,N5,N6,N5M,N6M,p,x,force,iopt)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      integer iopt
c      write(*,*)"entering dfunc3d"

      select case(iopt)
        case(1)
          CALL dwu(n,A,N5,N6,N5M,N6M,p,x,force,iopt)
        case(2)
          CALL dwu(n,A,N5,N6,N5M,N6M,p,x,force,iopt)
        case(3)
          CALL dextwu(n,A,N5,N6,N5M,N6M,p,x,force)
      end select

c      write(*,*)"leaving dfunc3d"
      return
      END


      SUBROUTINE dwu(n,A,N5,N6,N5M,N6M,p,x,force,iopt)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, gradient
      Real*8 p(nmax*3),x(nmax*3),force(ffmaxdim)
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6)
      rp=force(1)
      rh=force(2)
      ap=force(3)
      ah=force(4)
      frp=force(5)
      frh=force(6)
      fap=force(7)
      fah=force(8)
      fco=force(9)
C     Stretching
      Do I=1,n,3
        ehookx=0.d0
        ehooky=0.d0
        ehookz=0.d0
        I1=(I+2)/3
        Do J=1,n,3
          J1=(J+2)/3
          if(A(I1,J1).ne.0) then
            px=p(i)-p(j)
            py=p(i+1)-p(j+1)
            pz=p(i+2)-p(j+2)
            ratom=dsqrt(px**2 + py**2 + pz**2)
            ratominv=1.d0/ratom
C           Check if bond is part of 5-ring
            do IB=1,12
              ir1=0
              ir2=0
              do JB=1,5
                if(I1.eq.N5M(IB,JB)) ir1=1
                if(J1.eq.N5M(IB,JB)) ir2=1
              enddo
              if(ir1.eq.1.and.ir2.eq.1) then
C               5-ring
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
  1       continue
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
      if (iopt.eq.2 .and. fco.ne.0.d0)  then
       Do I=1,n/3
        rinv=(p(I)**2+p(I+1)**2+p(I+2)**2)**(-1.5d0)
        x(I*3-2)=x(I*3-2)-fco*rinv*p(I)
        x(I*3-1)=x(I*3-1)-fco*rinv*p(I+1)
        x(I*3)=x(I*3)-fco*rinv*p(I+2)
       enddo
      endif

      return
      END



      SUBROUTINE dextwu(n,A,N5,N6,N5M,N6M,p,x,force)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 p(nmax*3),x(nmax*3),force(ffmaxdim)
      real*8 J1x,J1y,J1z,J2x,J2y,J2z,J3x,J3y,J3z,J4x,J4y,J4z
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),pentagoncount,
     2 hexagoncount,arbitrary_index,neighbour_atoms(3),
     3 neighbour_faces_h(3),neighbour_faces_p(3)
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

C     Stretching
c we distinguish between bonds between two hexagons, two pentagons and hex/pent
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
C           Check if bond is part of 5-ring
            pentagoncount=0
            do IB=1,12! number of pentagons
              ir1=0
              ir2=0
              do JB=1,5!number of atoms per pentagons
                if(I1.eq.N5M(IB,JB)) ir1=1 !
                if(J1.eq.N5M(IB,JB)) ir2=1 ! if I1 and J2 happen to be in the same pentagon
              enddo
              if(ir1.eq.1 .and. ir2.eq.1) then
                pentagoncount=pentagoncount+1
              endif
            enddo
            if(pentagoncount.eq.0) then
C             6-ring, 6-ring
              fac=frhh*ratominv*(ratom-rhh)
            else if(pentagoncount.eq.1) then
C             5-ring, 6-ring
              fac=frhp*ratominv*(ratom-rhp)
            else
C             5-ring, 5-ring
              fac=frpp*ratominv*(ratom-rpp)
            endif
            ehookx=ehookx+fac*px! add up forces on a single atom, in x, y, z
            ehooky=ehooky+fac*py
            ehookz=ehookz+fac*pz
C           Check if bond is part of 5-ring
          endif
        enddo
        x(I)  =2.d0*ehookx
        x(I+1)=2.d0*ehooky
        x(I+2)=2.d0*ehookz
      enddo
c      write(*,*)x,"displacement after dists"

C     Bending
      Do I=1,N5 ! Loop over 5-rings (and N5 == 12)
        Do J=1,5 ! loop over atoms in pentagon
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=5
          if(JRX.eq.6) JRX=1
          JM=3*N5M(I,J)-2! position of x coordinate of middle atom in p
          JL=3*N5M(I,JLX)-2
          JR=3*N5M(I,JRX)-2
c left bond
           pxL=p(JM)  -p(JL)
           pyL=p(JM+1)-p(JL+1)
           pzL=p(JM+2)-p(JL+2)
          r2L=pxL*pxL+pyL*pyL+pzL*pzL
          r1L=dsqrt(r2L)
          r3L=r1L*r2L
c right bond
           pxR=p(JM)  -p(JR)
           pyR=p(JM+1)-p(JR+1)
           pzR=p(JM+2)-p(JR+2)
          r2R=pxR*pxR+pyR*pyR+pzR*pzR
          r1R=dsqrt(r2R)
          r3R=r1R*r2R
c no bond
           pxM=p(JL)  -p(JR)
           pyM=p(JL+1)-p(JR+1)
           pzM=p(JL+2)-p(JR+2)
          r2M=pxM*pxM+pyM*pyM+pzM*pzM
          r1M=dsqrt(r2M)
          r3M=r1M*r2M
c law of cosines
          cosarg=.5d0*(r2L+r2R-r2M)/(r1L*r1R)
          if(cosarg.gt.1.d0) cosarg=1.d0
          if(cosarg.lt.-1.d0) cosarg=-1.d0
          anglep=dacos(cosarg)
c          write(*,*)anglep,"pentagon angle"
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
C Derivative of left atom
           fac6=-fac/(r3L*r1R)
          x(JL)  =x(JL)  +fac6*(pxL*fac3-2.d0*pxM*r2L)
          x(JL+1)=x(JL+1)+fac6*(pyL*fac3-2.d0*pyM*r2L)
          x(JL+2)=x(JL+2)+fac6*(pzL*fac3-2.d0*pzM*r2L)
C Derivative of right atom
           fac7=-fac/(r3R*r1L)
          x(JR)  =x(JR)  +fac7*(pxR*fac2+2.d0*pxM*r2R)
          x(JR+1)=x(JR+1)+fac7*(pyR*fac2+2.d0*pyM*r2R)
          x(JR+2)=x(JR+2)+fac7*(pzR*fac2+2.d0*pzM*r2R)
        enddo
      enddo
c      write(*,*)x,"displacement after pentagons angles"
      
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
           angleh=dacos(cosarg)! angle in hexagon
c          write(*,*)angleh,"hexagon angle"
           anglesin=dabs(dsin(angleh)) ! sine of the angle
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
C       Derivative of left atom
           fac6=-fac/(r3L*r1R)
          x(JL)  =x(JL)  +fac6*(pxL*fac3-2.d0*pxM*r2L)
          x(JL+1)=x(JL+1)+fac6*(pyL*fac3-2.d0*pyM*r2L)
          x(JL+2)=x(JL+2)+fac6*(pzL*fac3-2.d0*pzM*r2L)
C       Derivative of right atom
           fac7=-fac/(r3R*r1L)
          x(JR)  =x(JR)  +fac7*(pxR*fac2+2.d0*pxM*r2R)
          x(JR+1)=x(JR+1)+fac7*(pyR*fac2+2.d0*pyM*r2R)
          x(JR+2)=x(JR+2)+fac7*(pzR*fac2+2.d0*pzM*r2R)
        enddo
      enddo
c      write(*,*)x,"displacement after hexagon angles"


C dihedrals 
      Do I=1,n/3 ! iterate over atoms
c classify vertex acording to adjacent faces
c find neighbouring faces (3)
        pentagoncount=0
        do IB=1,12 ! iterate over pentagons
          do JB=1,5 ! iterate over atoms in pentagon
            if(I.eq.N5M(IB,JB)) then
              pentagoncount=pentagoncount+1 ! find if I is part of 0,...,3 pentagons
              neighbour_faces_p(pentagoncount)=IB
            endif
          enddo
        enddo
        hexagoncount=0
        do IB=1,n/3-12 ! because n=vertex_count * 3 ()
          do JB=1,6
            if(i.eq.N6M(IB,JB)) then
              hexagoncount=hexagoncount+1
              neighbour_faces_h(hexagoncount)=IB
            endif
          enddo
        enddo
c find neighbouring atoms (3)
        neighbour_atom_count=1
        do j=1,n/3
          if(A(I,j).ne.0) then
            neighbour_atoms(neighbour_atom_count)=j
            neighbour_atom_count=neighbour_atom_count+1
          endif
        enddo
c        write(*,*)"p:",pentagoncount, "h:",hexagoncount
c sort neighbours
c we make use of the fact, that for any dihedral ijkl=ikjl
c therefore we only need to find the special vertex and dont care about the others which are much harder to distinguish
        if(pentagoncount.eq.1) then
          do k=1,3 ! iterate over neighbour atoms
            arbitrary_index=0
            do l=1,hexagoncount ! iterate over neighbour hexagons
              do m=1,6 ! iterate over atoms in hexagon
                if (neighbour_atoms(k).eq.N6M(l,m)) then
                  arbitrary_index=arbitrary_index+1
                endif
              enddo
            enddo
            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two hexagons
              buffer=neighbour_atoms(k)
              neighbour_atoms(k)=neighbour_atoms(1)
              neighbour_atoms(1)=buffer
            endif
          enddo
        endif
        if(pentagoncount.eq.2) then
          do k=1,3 ! iterate over neighbour atoms
            arbitrary_index=0
            do l=1,pentagoncount ! iterate over neighbour pentagons
              do m=1,5 ! iterate over atoms in hexagon
                if (neighbour_atoms(k).eq.N5M(l,m)) then
                  arbitrary_index=arbitrary_index+1
                endif
              enddo
            enddo
            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two pentagons
              buffer=neighbour_atoms(k)
              neighbour_atoms(k)=neighbour_atoms(1)
              neighbour_atoms(1)=buffer
            endif
          enddo
        endif
c atoms
        J1=neighbour_atoms(1)
        J2=neighbour_atoms(2)
        J3=neighbour_atoms(3)
        J4=I
c        write(*,*)j1,j2,j3,j4
c coordinates
        J1x=p(J1*3-2)
        J1y=p(J1*3-1)
        J1z=p(J1*3)
        J2x=p(J2*3-2)
        J2y=p(J2*3-1)
        J2z=p(J2*3)
        J3x=p(J3*3-2)
        J3y=p(J3*3-1)
        J3z=p(J3*3)
        J4x=p(J4*3-2)
        J4y=p(J4*3-1)
        J4z=p(J4*3)
c        write(*,*)J1x,J1y,J1z,J2x,J2y,J2z,J3x,J3y,J3z,J4x,J4y,J4z
c some auxiliary factors without any physical meaning
        fac_A=(-J1y*J2x+J1x*J2y+J1y*J3x-J2y*J3x-J1x*J3y+J2x*J3y)
        fac_B=(-J2y*J3x+J2x*J3y+J2y*J4x-J3y*J4x-J2x*J4y+J3x*J4y)
        fac_C=(J1z*J2x-J1x*J2z-J1z*J3x+J2z*J3x+J1x*J3z-J2x*J3z)
        fac_D=(J2z*J3x-J2x*J3z-J2z*J4x+J3z*J4x+J2x*J4z-J3x*J4z)
        fac_E=(-J1z*J2y+J1y*J2z+J1z*J3y-J2z*J3y-J1y*J3z+J2y*J3z)
        fac_F=(-J2z*J3y+J2y*J3z+J2z*J4y-J3z*J4y-J2y*J4z+J3y*J4z)
        fac_G=(fac_A*fac_B+fac_C*fac_D+fac_E*fac_F)
        fac_H=(fac_A**2+fac_C**2+fac_E**2)
        fac_I=(fac_B**2+fac_D**2+fac_F**2)
        fac_J=1/(fac_H * fac_I)
        fac_K=1/(fac_H * fac_I**2)
        fac_L=1/(fac_H**2 * fac_I)
        fac_AC=fac_G * fac_J
        write(*,*)"abcdef",fac_a, fac_b, fac_c, fac_d, fac_e, fac_f
        if(fac_AC .ge. dpi) fac_AC=fac_AC-2*dpi
        fac_ac=dabs(fac_ac)
        write(*,*)"ac",fac_ac
        select case(pentagoncount)
          case(0)
          fac_M=2*(fdhhh*(dhhh-dacos(fac_AC)))/
     2          (dsqrt(1-fac_AC**2))
          case(1)
          fac_M=2*(fdhhp*(dhhp-dacos(fac_AC)))/
     2          (dsqrt(1-fac_AC**2))
          case(2)
          fac_M=2*(fdhpp*(dhpp-dacos(fac_AC)))/
     2          (dsqrt(1-fac_AC**2))
          case(3)
          fac_M=2*(fdppp*(dppp-dacos(fac_AC)))/
     2          (dsqrt(1-fac_AC**2))
        end select
        fac_N=(J2y-J3y)
        fac_O=(J2z-J3z)
        fac_P=(J2x-J3x)
        fac_Q=(J3y-J4y)
        fac_R=(J3z-J4z)
        fac_S=(J1y-J3y)
        fac_T=(J1z-J3z)
        fac_U=(J3x-J4x)
        fac_V=(J1x-J3x)
        fac_W=(J2y-J4y)
        fac_X=(J2z-J4z)
        fac_Y=(J1y-J2y)
        fac_Z=(J1z-J2z)
        fac_AA=(J2x-J4x)
        fac_AB=(J1x-J2x)
        write(*,*)"mnopqr",fac_m, fac_n, fac_o, fac_p, fac_q, fac_r
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(J1*3-2)=x(J1*3-2)+
     2     ((fac_N*fac_B-fac_O*fac_D)*fac_J+
     3     (-2*fac_N*fac_A+2*fac_O*fac_C)*fac_G*fac_L)*fac_M
        x(J1*3-1)=x(J1*3-1)+
     2     ((-fac_P*fac_B+fac_O*fac_F)*fac_J-
     3     ((-2*fac_P*fac_A+2*fac_O*fac_E)*fac_G)*fac_L)*fac_M
        x(J1*3)=x(J1*3)+
     2     ((fac_P*fac_D-fac_N*fac_F)*fac_J-
     3     ((2*fac_P*fac_C-2*fac_N*fac_E)*fac_G)*fac_L)*fac_M

        x(J2*3-2)=x(J2*3-2)+
     2     (-(2*fac_Q*fac_B-2*fac_R*fac_D)*fac_G*fac_K+
     3     (fac_A*fac_Q-fac_S*fac_B-fac_C*fac_R+fac_T*fac_D)*fac_J-
     4     (-2*fac_S*fac_A+2*fac_T*fac_C)*fac_G*fac_L)*fac_M
        x(J2*3-1)=x(J2*3-1)+
     2     fac_M*(fac_E*fac_J*fac_R-2*fac_F*fac_G*fac_K*fac_R-fac_F*
     3     fac_J*fac_T+2*fac_E*fac_G*fac_L*fac_T-fac_A*fac_J*fac_U+
     4     2*fac_B*fac_G*fac_K*fac_U+fac_B*fac_J*fac_V-
     5     2*fac_A*fac_G*fac_L*fac_V)
        x(J2*3)=x(J2*3)+
     2     fac_M*(2*fac_F*fac_G*fac_K*fac_Q+fac_F*fac_J*fac_S-
     3     fac_E*(fac_J*fac_Q+2*fac_G*fac_L*fac_S)+fac_C*fac_J*fac_U-
     4     2*fac_D*fac_G*fac_K*fac_U-fac_D*fac_J*fac_V+
     5     2*fac_C*fac_G*fac_L*fac_V)

        x(J3*3-2)=x(J3*3-2)+
     2     fac_M*(2*fac_B*fac_G*fac_K*fac_W+fac_C*fac_J*fac_X-
     3     2*fac_D*fac_G*fac_K*fac_X+fac_B*fac_J*fac_Y-
     4     fac_A*(fac_J*fac_W+2*fac_G*fac_L*fac_Y)-fac_D*fac_J*fac_Z+
     5     2*fac_C*fac_G*fac_L*fac_Z)
        x(J3*3-1)=x(J3*3-1)+
     2     fac_M*(fac_A*fac_AA*fac_J-fac_AB*fac_B*fac_J-
     3     2*fac_AA*fac_B*fac_G*fac_K+2*fac_A*fac_AB*fac_G*fac_L-
     4     fac_E*fac_J*fac_X+2*fac_F*fac_G*fac_K*fac_X+
     5     fac_F*fac_J*fac_Z-2*fac_E*fac_G*fac_L*fac_Z)
        x(J3*3)=x(J3*3)+
     2     fac_M*(-fac_AA*fac_C*fac_J+fac_AB*fac_D*fac_J+
     3     2*fac_AA*fac_D*fac_G*fac_K-2*fac_AB*fac_C*fac_G*fac_L+
     4     fac_E*fac_J*fac_W-2*fac_F*fac_G*fac_K*fac_W-
     5     fac_F*fac_J*fac_Y+2*fac_E*fac_G*fac_L*fac_Y)

        x(J4*3-2)=x(J4*3-2)+
     2     fac_M*(fac_A*fac_J*fac_N-2*fac_B*fac_G*fac_K*fac_N-
     3     fac_C*fac_J*fac_O+2*fac_D*fac_G*fac_K*fac_O)
        x(J4*3-1)=x(J4*3-1)+
     2     fac_M*(fac_E*fac_J*fac_O-2*fac_F*fac_G*fac_K*fac_O-
     3     fac_A*fac_J*fac_P+2*fac_B*fac_G*fac_K*fac_P)
        x(J4*3)=x(J4*3)+
     2     fac_M*(-fac_E*fac_J*fac_N+2*fac_F*fac_G*fac_K*fac_N+
     3     fac_C*fac_J*fac_P-2*fac_D*fac_G*fac_K*fac_P)
      enddo
c      write(*,*)"x: ",x
c      write(*,*)"leaving dextwu"
      return
      END
