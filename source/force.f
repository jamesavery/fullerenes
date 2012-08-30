      SUBROUTINE func3d(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt)
c n=MATOM*3
      use config
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


      SUBROUTINE wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)

C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, energy
      Real*8 p(nmax*3),force(ffmaxdim)
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),iopt
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
      Do I=1,n,3! iterate over half the adjacency mtx
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
        JM=3*N5M(I,J)
        JL=3*N5M(I,JLX)
        JR=3*N5M(I,JRX)
        call angle(p(JL-2),p(JL-1),p(JL),p(JM-2),p(JM-1),p(JM),
     2   p(JR-2),p(JR-1),p(JR),angle_p)
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
        JM=3*N6M(I,J)
        JL=3*N6M(I,JLX)
        JR=3*N6M(I,JRX)
        call angle(p(JL-2),p(JL-1),p(JL),p(JM-2),p(JM-1),p(JM),
     2   p(JR-2),p(JR-1),p(JR),angle_h)
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
          JM=3*N5M(I,J) ! middle
          JL=3*N5M(I,JLX) ! left
          JR=3*N5M(I,JRX) ! right
          call angle(p(JL-2),p(JL-1),p(JL),p(JM-2),p(JM-1),p(JM),
     2     p(JR-2),p(JR-1),p(JR),angle_p)
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
          JM=3*N6M(I,J)
          JL=3*N6M(I,JLX)
          JR=3*N6M(I,JRX)
          call angle(p(JL-2),p(JL-1),p(JL),p(JM-2),p(JM-1),p(JM),
     2     p(JR-2),p(JR-1),p(JR),angle_h)
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
        J1=3*neighbour_atoms(1)
        J2=3*neighbour_atoms(2)
        J3=3*neighbour_atoms(3)
        J4=3*I
c coordinates
        call dihedral(p(J1-2),p(J1-1),p(J1),p(J2-2),p(J2-1),
     2   p(J2),p(J3-2),p(J3-1),p(J3),p(J4-2),p(J4-1),
     3   p(J4),angle_abcd)
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
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        increment=(angle_abcd-zero_value)**2
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
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),iopt
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
        Do J=I+3,n,3
          J1=(J+2)/3
          if(A(I1,J1).ne.0) then
            ax=p(i)
            ay=p(i+1)
            az=p(i+2)
            bx=p(j)
            by=p(j+1)
            bz=p(j+2)
            call DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,ratom)
C           Check if bond is part of 5-ring
            pentagoncount=0
            do IB=1,12
              ir1=0
              ir2=0
              do JB=1,5
                if(I1.eq.N5M(IB,JB)) ir1=1
                if(J1.eq.N5M(IB,JB)) ir2=1
              enddo
              if(ir1.eq.1.and.ir2.eq.1) then
C               5-ring
                pentagoncount=pentagoncount+1
                go to 1
              endif
            enddo
C           6-ring
 1          if (pentagoncount.eq.0)then
              zero_value=rh
              force_constant=frh
            else
              zero_value=rp
              force_constant=frp
            end if
            dE_over_dc=force_constant*(ratom-zero_value)
            x(i  )=x(i  )+dax*dE_over_dc
            x(i+1)=x(i+1)+day*dE_over_dc
            x(i+2)=x(i+2)+daz*dE_over_dc
            x(j  )=x(j  )+dbx*dE_over_dc
            x(j+1)=x(j+1)+dby*dE_over_dc
            x(j+2)=x(j+2)+dbz*dE_over_dc
          endif
        enddo
      enddo
        
C     Bending
C     Loop over 5-rings
      Do I=1,N5 !=12
        Do J=1,5
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=5
          if(JRX.eq.6) JRX=1
          JL=3*N5M(I,JLX)
          JM=3*N5M(I,J)
          JR=3*N5M(I,JRX)
          ax=p(JL-2)
          ay=p(JL-1)
          az=p(JL)
          bx=p(JM-2)
          by=p(JM-1)
          bz=p(JM)
          cx=p(JR-2)
          cy=p(JR-1)
          cz=p(JR)
          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3     angle_abc)
          zero_value=ap
          force_constant=fap
          dE_over_dc=force_constant*(angle_abc-zero_value)
          x(JL-2)=x(JL-2)+dax*dE_over_dc
          x(JL-1)=x(JL-1)+day*dE_over_dc
          x(JL)  =x(JL)  +daz*dE_over_dc
          x(JM-2)=x(JM-2)+dbx*dE_over_dc
          x(JM-1)=x(JM-1)+dby*dE_over_dc
          x(JM)  =x(JM)  +dbz*dE_over_dc
          x(JR-2)=x(JR-2)+dcx*dE_over_dc
          x(JR-1)=x(JR-1)+dcy*dE_over_dc
          x(JR)  =x(JR)  +dcz*dE_over_dc
        enddo
      enddo
      
C     Loop over 6-rings
      Do I=1,N6
        Do J=1,6
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=6
          if(JRX.eq.7) JRX=1
          JL=3*N6M(I,JLX)
          JM=3*N6M(I,J)
          JR=3*N6M(I,JRX)
          ax=p(JL-2)
          ay=p(JL-1)
          az=p(JL)
          bx=p(JM-2)
          by=p(JM-1)
          bz=p(JM)
          cx=p(JR-2)
          cy=p(JR-1)
          cz=p(JR)
          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3     angle_abc)
          zero_value=ah
          force_constant=fah
          dE_over_dc=force_constant*(angle_abc-zero_value)
c      write(*,*)"a,z,de",angle_abc,zero_value,dE_over_dc
c      write(*,*)"dd",dax*dE_over_dc,day*dE_over_dc,
c     2  daz*dE_over_dc,
c     2  dbx*dE_over_dc,dby*dE_over_dc,dbz*dE_over_dc,
c     2  dcx*dE_over_dc,dcy*dE_over_dc,dcz*dE_over_dc
          x(JL-2)=x(JL-2)+dax*dE_over_dc
          x(JL-1)=x(JL-1)+day*dE_over_dc
          x(JL)  =x(JL)  +daz*dE_over_dc
          x(JM-2)=x(JM-2)+dbx*dE_over_dc
          x(JM-1)=x(JM-1)+dby*dE_over_dc
          x(JM)  =x(JM)  +dbz*dE_over_dc
          x(JR-2)=x(JR-2)+dcx*dE_over_dc
          x(JR-1)=x(JR-1)+dcy*dE_over_dc
          x(JR)  =x(JR)  +dcz*dE_over_dc
        enddo
      enddo

C     Coulomb repulsion from origin
      if (iopt.eq.2 .and. fco.ne.0.d0)  then
        Do I=1,n,3
          rinv=(p(I)**2+p(I+1)**2+p(I+2)**2)**(-1.5d0)
          x(I)  =x(I)  +fco*rinv*p(I)
          x(I+1)=x(I+1)+fco*rinv*p(I+1)
          x(I+2)=x(I+2)+fco*rinv*p(I+2)
        enddo
      endif

      return
      END



      SUBROUTINE dextwu(n,A,N5,N6,N5M,N6M,p,x,force)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 p(nmax*3),x(nmax*3),force(ffmaxdim)
c      real*8 J1x,J1y,J1z,J2x,J2y,J2z,J3x,J3y,J3z,J4x,J4y,J4z
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),pentagoncount,
     2 hexagoncount,arbitrary_index,neighbour_atoms(3),
     3 neighbour_faces_h(3),neighbour_faces_p(3),buffer
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
C     Stretching
      Do I=1,n,3
        i1=(i+2)/3
        Do J=I+3,n,3
          j1=(j+2)/3
c check if bond exists
          if(A(i1,j1).ne.0) then
c get coordinates
            ax=p(I)
            ay=p(I+1)
            az=p(I+2)
            bx=p(J)
            by=p(J+1)
            bz=p(J+2)
            call DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,ratom)
C           Check if bond is part of 5-ring
            pentagoncount=0
            do IB=1,12
              ir1=0
              ir2=0
              do JB=1,5
                if(i1.eq.N5M(IB,JB)) ir1=1
                if(j1.eq.N5M(IB,JB)) ir2=1
              enddo
              if(ir1.eq.1.and.ir2.eq.1) then
                pentagoncount=pentagoncount+1
              endif
            enddo
            select case(pentagoncount)
            case(0)
              zero_value=rhh
              force_constant=frhh
            case(1)
              zero_value=rhp
              force_constant=frhp
            case(2)
              zero_value=rpp
              force_constant=frpp
            end select
            dE_over_dc=force_constant*(ratom-zero_value)
            x(I)  =x(I)  +dax*dE_over_dc
            x(I+1)=x(I+1)+day*dE_over_dc
            x(I+2)=x(I+2)+daz*dE_over_dc
            x(J)  =x(J)  +dbx*dE_over_dc
            x(J+1)=x(J+1)+dby*dE_over_dc
            x(J+2)=x(J+2)+dbz*dE_over_dc
          endif
        enddo
      enddo


C     Bending        
C     Loop over 5-rings
      Do I=1,N5 !and n5==12
        Do J=1,5
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=5
          if(JRX.eq.6) JRX=1
          JL=3*N5M(I,JLX)
          JM=3*N5M(I,J)
          JR=3*N5M(I,JRX)
          ax=p(JL-2)
          ay=p(JL-1)
          az=p(JL)
          bx=p(JM-2)
          by=p(JM-1)
          bz=p(JM)
          cx=p(JR-2)
          cy=p(JR-1)
          cz=p(JR)
          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3     angle_abc)
c          call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
          zero_value=ap
          force_constant=fap
          dE_over_dc=force_constant*(angle_abc-zero_value)
          x(JL-2)=x(JL-2)+dax*dE_over_dc
          x(JL-1)=x(JL-1)+day*dE_over_dc
          x(JL)  =x(JL)  +daz*dE_over_dc
          x(JM-2)=x(JM-2)+dbx*dE_over_dc
          x(JM-1)=x(JM-1)+dby*dE_over_dc
          x(JM)  =x(JM)  +dbz*dE_over_dc
          x(JR-2)=x(JR-2)+dcx*dE_over_dc
          x(JR-1)=x(JR-1)+dcy*dE_over_dc
          x(JR)  =x(JR)  +dcz*dE_over_dc
        enddo
      enddo
      
C     Loop over 6-rings
      Do I=1,N6
        Do J=1,6
          JLX=J-1
          JRX=J+1
          if(JLX.eq.0) JLX=6
          if(JRX.eq.7) JRX=1
          JL=3*N6M(I,JLX)
          JM=3*N6M(I,J)
          JR=3*N6M(I,JRX)
          ax=p(JL-2)
          ay=p(JL-1)
          az=p(JL)
          bx=p(JM-2)
          by=p(JM-1)
          bz=p(JM)
          cx=p(JR-2)
          cy=p(JR-1)
          cz=p(JR)
          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3     angle_abc)
          zero_value=ah
          force_constant=fah
          dE_over_dc=force_constant*(angle_abc-zero_value)
          x(JL-2)=x(JL-2)+dax*dE_over_dc
          x(JL-1)=x(JL-1)+day*dE_over_dc
          x(JL)  =x(JL)  +daz*dE_over_dc
          x(JM-2)=x(JM-2)+dbx*dE_over_dc
          x(JM-1)=x(JM-1)+dby*dE_over_dc
          x(JM)  =x(JM)  +dbz*dE_over_dc
          x(JR-2)=x(JR-2)+dcx*dE_over_dc
          x(JR-1)=x(JR-1)+dcy*dE_over_dc
          x(JR)  =x(JR)  +dcz*dE_over_dc
        enddo
      enddo

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
c we make use of the fact, that for any dihedral ijkl=-ikjl
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
        else if(pentagoncount.eq.2) then
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
        J1=3*neighbour_atoms(1)
        J2=3*neighbour_atoms(2)
        J3=3*neighbour_atoms(3)
        J4=3*I
c        write(*,*)j1,j2,j3,j4
c coordinates
        ax=p(J1-2)
        ay=p(J1-1)
        az=p(J1)
        bx=p(J2-2)
        by=p(J2-1)
        bz=p(J2)
        cx=p(J3-2)
        cy=p(J3-1)
        cz=p(J3)
        dx=p(J4-2)
        dy=p(J4-1)
        dz=p(J4)
        call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        select case(pentagoncount)
          case(0)
            zero_value=dhhh
            force_constant=fdhhh
          case(1)
            zero_value=dhhp
            force_constant=fdhhp
          case(2)
            zero_value=dhpp
            force_constant=fdhpp
          case(3)
            zero_value=dppp
            force_constant=fdppp
        end select
        dE_over_dc=force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(J1-2)=x(J1-2)+dax*dE_over_dc
        x(J1-1)=x(J1-1)+day*dE_over_dc
        x(J1)  =x(J1)  +daz*dE_over_dc
        x(J2-2)=x(J2-2)+dbx*dE_over_dc
        x(J2-1)=x(J2-1)+dby*dE_over_dc
        x(J2)  =x(J2)  +dbz*dE_over_dc
        x(J3-2)=x(J3-2)+dcx*dE_over_dc
        x(J3-1)=x(J3-1)+dcy*dE_over_dc
        x(J3)  =x(J3)  +dcz*dE_over_dc
        x(J4-2)=x(J4-2)+ddx*dE_over_dc
        x(J4-1)=x(J4-1)+ddy*dE_over_dc
        x(J4)  =x(J4)  +ddz*dE_over_dc
      enddo
c      write(*,*)"d,0: ",angle_abcd,zero_value," (should be similar)"
      return
      END
