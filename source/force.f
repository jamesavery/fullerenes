      SUBROUTINE func3d(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p,
     1  d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
c n=MATOM*3
      use config
c      IMPLICIT REAL*8 (A-H,O-Z)
      integer iopt

c      write(*,*)'entering func3d'

      select case(iopt)
        case(1)
          CALL wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt,
     1      e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1      a_h,a_p)
        case(2)
          CALL wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt,
     1      e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1      a_h,a_p)
        case(3)
          CALL extwu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,
     1      e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     2      a_h,a_p,
     1      d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      end select

c      write(*,*)'leaving func3d'
      
      return
      END SUBROUTINE


      SUBROUTINE wu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,iopt,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p)
c n=MATOM*3
      use config
      IMPLICIT REAL*8 (A-H,O-Z)

C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, energy
      Real*8 p(nmax*3),force(ffmaxdim)
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),iopt

c     edges with 0, 1, 2 pentagons
      integer e_hh(2,N/2), e_hp(2,N/2), e_pp(2,N/2)
      integer a_p(3,60), a_h(3,N-60)
c     counter for edges with 0, 1, 2 pentagons neighbours
      integer ne_hh,ne_hp,ne_pp

c      write(*,*)'entering wu'

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

      if(ne_hh .ne. 0) then
        do i=1,ne_hh
          call dist(p(3*e_hh(1,i)-2),p(3*e_hh(1,i)-1),p(3*e_hh(1,i)),
     1              p(3*e_hh(2,i)-2),p(3*e_hh(2,i)-1),p(3*e_hh(2,i)),
     1              ratom)
          ehookrh=ehookrh+(ratom-rh)**2
        enddo
      endif

c     in the wu ff there is no difference between 1 and 2 pentagons
      if(ne_hp .ne. 0) then
        do i=1,ne_hp
          call dist(p(3*e_hp(1,i)-2),p(3*e_hp(1,i)-1),p(3*e_hp(1,i)),
     1              p(3*e_hp(2,i)-2),p(3*e_hp(2,i)-1),p(3*e_hp(2,i)),
     2              ratom)
          ehookrp=ehookrp+(ratom-rp)**2
        enddo
      endif

      if(ne_pp .ne. 0) then
        do i=1,ne_pp
          call dist(p(3*e_pp(1,i)-2),p(3*e_pp(1,i)-1),p(3*e_pp(1,i)),
     1              p(3*e_pp(2,i)-2),p(3*e_pp(2,i)-1),p(3*e_pp(2,i)),
     1              ratom)
          ehookrp=ehookrp+(ratom-rp)**2
        enddo
      endif

c      Do I=1,n,3! iterate over half the adjacency mtx
c        I1=(I+2)/3
c        Do J=I+3,n,3
c          J1=(J+2)/3
c          if(A(I1,J1).ne.0) then
c            call dist(p(i),p(i+1),p(i+2),p(j),p(j+1),p(j+2),ratom)
C           Check if bond is part of 5-ring
c            do IB=1,12
c              ir1=0
c              ir2=0
c              do JB=1,5
c                if(I1.eq.N5M(IB,JB)) ir1=1
c                if(J1.eq.N5M(IB,JB)) ir2=1
c              enddo
c              if(ir1.eq.1.and.ir2.eq.1) then
c               5-ring
c                ehookrp=ehookrp+(ratom-rp)**2
c                go to 1
c              endif
c            enddo
c           6-ring
c            ehookrh=ehookrh+(ratom-rh)**2
c          endif
c  1       continue
c        enddo
c      enddo
     
c      write(*,*)'before angle - p'
c      write(*,*)'p',p
c      write(*,*)'a_p', a_p

C     Bending
      ehookap=0.d0
      ehookah=0.d0
C     Loop over 5-rings
      do i=1,60
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
        call angle(p(3*a_p(1,i)-2),p(3*a_p(1,i)-1),p(3*a_p(1,i)),
     1             p(3*a_p(2,i)-2),p(3*a_p(2,i)-1),p(3*a_p(2,i)),
     1             p(3*a_p(3,i)-2),p(3*a_p(3,i)-1),p(3*a_p(3,i)),
     1             angle_p)
        ehookap=ehookap+(angle_p-ap)**2
c        write(*,*)'diff',angle_p,ap
      enddo
c      write(*,*)'ehookap',ehookap

c      write(*,*)'before angle - h'

c      ehookap=0.d0

c      Do I=1,N5
c      Do J=1,5
c        JLX=J-1
c        JRX=J+1
c        if(JLX.eq.0) JLX=5
c        if(JRX.eq.6) JRX=1
c        JM=3*N5M(I,J)
c        JL=3*N5M(I,JLX)
c        JR=3*N5M(I,JRX)
c        call angle(p(JL-2),p(JL-1),p(JL),p(JM-2),p(JM-1),p(JM),
c     2   p(JR-2),p(JR-1),p(JR),angle_p)
c        ehookap=ehookap+(angle_p-ap)**2
c        write(*,*)'diff',angle_p,ap
c      enddo
c      enddo

c      write(*,*)'n',n

      if(n/3 .gt. 20) then
C     Loop over 6-rings
      do i=1,n-60
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
        call angle(p(3*a_h(1,i)-2),p(3*a_h(1,i)-1),p(3*a_h(1,i)),
     1             p(3*a_h(2,i)-2),p(3*a_h(2,i)-1),p(3*a_h(2,i)),
     1             p(3*a_h(3,i)-2),p(3*a_h(3,i)-1),p(3*a_h(3,i)),
     1             angle_h)
        ehookah=ehookah+(angle_h-ah)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif

c      if(N6.eq.0) go to 2
c      Do I=1,N6
c        Do J=1,6
c          JLX=J-1
c          JRX=J+1
c          if(JLX.eq.0) JLX=6
c          if(JRX.eq.7) JRX=1
c          JM=3*N6M(I,J)
c          JL=3*N6M(I,JLX)
c          JR=3*N6M(I,JRX)
c          call angle(p(JL-2),p(JL-1),p(JL),p(JM-2),p(JM-1),p(JM),
c     2      p(JR-2),p(JR-1),p(JR),angle_h)
c          ehookah=ehookah+(angle_h-ah)**2
c        enddo
c      enddo


C     Coulomb repulsion from origin
      ecoulomb=0.d0
      if (iopt.eq.2 .and. fco.ne.0.d0)  then
        Do I=1,n,3
          rinv=1.d0/dsqrt(p(I)**2+p(I+1)**2+p(I+2)**2)
          ecoulomb=ecoulomb+rinv
        enddo
      endif
 

C     total energy  
      fc=frp*ehookrp+frh*ehookrh+fap*ehookap+fah*ehookah+fco*ecoulomb
c  2   fc=frp*ehookrp+frh*ehookrh+fap*ehookap+fah*ehookah+fco*ecoulomb
c      write(*,*)'leaving wu'
      Return
      END SUBROUTINE



      SUBROUTINE extwu(n,IERR,A,N5,N6,N5M,N6M,p,fc,force,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p,
     1  d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
c n=MATOM*3
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 p(nmax*3),force(ffmaxdim),increment
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),neighbour_atoms(3),
     2 neighbour_faces_h(3),neighbour_faces_p(3),pentagoncount,
     3 hexagoncount,arbitrary_index

c     edges with 0, 1, 2 pentagons
      integer e_hh(2,N/2), e_hp(2,N/2), e_pp(2,N/2)
      integer a_p(3,60), a_h(3,N-60)
      integer d_hhh(4,n/3),d_hpp(4,n/3),d_hhp(4,n/3),d_ppp(4,n/3)
c     counter for edges with 0, 1, 2 pentagons neighbours
      integer ne_hh,ne_hp,ne_pp
c     counter for dihedrals with 0, 1, 2, 3 pentagons neighbours
      integer nd_hhh,nd_hhp,nd_hpp,nd_ppp

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

      if(ne_hh .ne. 0) then
        do i=1,ne_hh
          call dist(p(3*e_hh(1,i)-2),p(3*e_hh(1,i)-1),p(3*e_hh(1,i)),
     1              p(3*e_hh(2,i)-2),p(3*e_hh(2,i)-1),p(3*e_hh(2,i)),
     1              ratom)
          ehookrhh=ehookrhh+(ratom-rhh)**2
        enddo
      endif

      if(ne_hp .ne. 0) then
        do i=1,ne_hp
          call dist(p(3*e_hp(1,i)-2),p(3*e_hp(1,i)-1),p(3*e_hp(1,i)),
     1              p(3*e_hp(2,i)-2),p(3*e_hp(2,i)-1),p(3*e_hp(2,i)),
     2              ratom)
          ehookrhp=ehookrhp+(ratom-rhp)**2
        enddo
      endif

      if(ne_pp .ne. 0) then
        do i=1,ne_pp
          call dist(p(3*e_pp(1,i)-2),p(3*e_pp(1,i)-1),p(3*e_pp(1,i)),
     1              p(3*e_pp(2,i)-2),p(3*e_pp(2,i)-1),p(3*e_pp(2,i)),
     1              ratom)
          ehookrpp=ehookrpp+(ratom-rpp)**2
        enddo
      endif

c      Do I=1,n,3 ! n = number of atoms * 3 !!!
c        I1=(I+2)/3 ! I1 = 1, 2, ... (n+2)/3
c        Do J=I+3,n,3
c          J1=(J+2)/3 ! J1 = I1, I1+1 ... (n+2)/3
c          if(A(I1,J1).ne.0) then ! if connected
c get distance
c            call dist(p(i),p(i+1),p(i+2),p(j),p(j+1),p(j+2),ratom)
C Check if bond is part of 5-ring
c            pentagoncount=0
c            do IB=1,12! number of pentagons
c              ir1=0
c              ir2=0
c              do JB=1,5!number of atoms per pentagons
c                if(I1.eq.N5M(IB,JB)) ir1=1 !
c                if(J1.eq.N5M(IB,JB)) ir2=1 ! if I1 and J2 happen to be in the same pentagon
c              enddo
c              if(ir1.eq.1 .and. ir2.eq.1) then
c                pentagoncount=pentagoncount+1
c              endif
c            enddo
c            if(pentagoncount.eq.0) then
C             6-ring, 6-ring
c              ehookrhh=ehookrhh+(ratom-rhh)**2
c            else if(pentagoncount.eq.1) then
C             5-ring, 6-ring
c              ehookrhp=ehookrhp+(ratom-rhp)**2
c            else
C             5-ring, 5-ring
c              ehookrpp=ehookrpp+(ratom-rpp)**2
c            endif
c          endif ! connected
c        enddo
c      enddo


C     Bending
c     we distinguish between angles of pentagons and hexagons
C     Loop over 5-rings
      ehookap=0.d0
      do i=1,60
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
        call angle(p(3*a_p(1,i)-2),p(3*a_p(1,i)-1),p(3*a_p(1,i)),
     1             p(3*a_p(2,i)-2),p(3*a_p(2,i)-1),p(3*a_p(2,i)),
     1             p(3*a_p(3,i)-2),p(3*a_p(3,i)-1),p(3*a_p(3,i)),
     1             angle_abc)
        ehookap=ehookap+(angle_abc-ap)**2
c        write(*,*)'diff',angle_p,ap
      enddo

      if(n/3 .gt. 20) then
C     Loop over 6-rings
      ehookah=0.d0
      do i=1,n-60
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
        call angle(p(3*a_h(1,i)-2),p(3*a_h(1,i)-1),p(3*a_h(1,i)),
     1             p(3*a_h(2,i)-2),p(3*a_h(2,i)-1),p(3*a_h(2,i)),
     1             p(3*a_h(3,i)-2),p(3*a_h(3,i)-1),p(3*a_h(3,i)),
     1             angle_abc)
        ehookah=ehookah+(angle_abc-ah)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif


C dihedrals 
      ehookdppp=0.d0
      ehookdhpp=0.d0
      ehookdhhp=0.d0
      ehookdhhh=0.d0
c      write(*,*)nd_hhh,nd_hhp,nd_hpp,nd_ppp

      if(nd_hhh .ne. 0) then
C     3 hexagons
      do i=1,nd_hhh
      call dihedral(p(3*d_hhh(1,i)-2),p(3*d_hhh(1,i)-1),p(3*d_hhh(1,i)),
     1             p(3*d_hhh(2,i)-2),p(3*d_hhh(2,i)-1),p(3*d_hhh(2,i)),
     1             p(3*d_hhh(3,i)-2),p(3*d_hhh(3,i)-1),p(3*d_hhh(3,i)),
     1             p(3*d_hhh(4,i)-2),p(3*d_hhh(4,i)-1),p(4*d_hhh(3,i)),
     1             angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        ehookdhhh=ehookdhhh+(angle_abcd-dhhh)**2
      enddo
      endif


      if(nd_hhp .ne. 0) then
C     2 hexagons, 1 pentagon
      do i=1,nd_hhp
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
      call dihedral(p(3*d_hhp(1,i)-2),p(3*d_hhp(1,i)-1),p(3*d_hhp(1,i)),
     1             p(3*d_hhp(2,i)-2),p(3*d_hhp(2,i)-1),p(3*d_hhp(2,i)),
     1             p(3*d_hhp(3,i)-2),p(3*d_hhp(3,i)-1),p(3*d_hhp(3,i)),
     1             p(3*d_hhp(4,i)-2),p(3*d_hhp(4,i)-1),p(4*d_hhp(3,i)),
     1             angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        ehookdhhp=ehookdhhp+(angle_abcd-dhhp)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif


      if(nd_hpp .ne. 0) then
C     1 hexagon, 2 pentagons
      do i=1,nd_hpp
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
      call dihedral(p(3*d_hpp(1,i)-2),p(3*d_hpp(1,i)-1),p(3*d_hpp(1,i)),
     1             p(3*d_hpp(2,i)-2),p(3*d_hpp(2,i)-1),p(3*d_hpp(2,i)),
     1             p(3*d_hpp(3,i)-2),p(3*d_hpp(3,i)-1),p(3*d_hpp(3,i)),
     1             p(3*d_hpp(4,i)-2),p(3*d_hpp(4,i)-1),p(4*d_hpp(3,i)),
     1             angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        ehookdhpp=ehookdhpp+(angle_abcd-dhpp)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif


      if(nd_ppp .ne. 0) then
C     3 pentagons
      do i=1,nd_ppp
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
      call dihedral(p(3*d_ppp(1,i)-2),p(3*d_ppp(1,i)-1),p(3*d_ppp(1,i)),
     1             p(3*d_ppp(2,i)-2),p(3*d_ppp(2,i)-1),p(3*d_ppp(2,i)),
     1             p(3*d_ppp(3,i)-2),p(3*d_ppp(3,i)-1),p(3*d_ppp(3,i)),
     1             p(3*d_ppp(4,i)-2),p(3*d_ppp(4,i)-1),p(4*d_ppp(3,i)),
     1             angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        ehookdppp=ehookdppp+(angle_abcd-dppp)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif



c      Do I=1,n/3 ! iterate over atoms
c count adjacent pentagons (0 to 3)
c        pentagoncount=0
c        do IB=1,12 ! iterate over pentagons
c          do JB=1,5 ! iterate over atoms in pentagon
c            if(I.eq.N5M(IB,JB)) then
c              pentagoncount=pentagoncount+1 ! find if I is part of 0,...,3 pentagons
c              neighbour_faces_p(pentagoncount)=IB
c            endif
c          enddo
c        enddo
c count adjacent hexagons (0 to 3)
c        hexagoncount=0
c        do IB=1,n/3-12 ! because n=vertex_count * 3
c          do JB=1,6
c            if(i.eq.N6M(IB,JB)) then
c              hexagoncount=hexagoncount+1
c              neighbour_faces_h(hexagoncount)=IB
c            endif
c          enddo
c        enddo
c find neighbouring atoms (3)
c        neighbour_atom_count=1
c        do j=1,n/3
c          if(A(I,j).ne.0) then
c            neighbour_atoms(neighbour_atom_count)=j
c            neighbour_atom_count=neighbour_atom_count+1
c          endif
c        enddo
c sort neighbours
c we make use of the fact, that for any dihedral ijkl=-ikjl (hence: ijkl=abs(ikjl))
c therefore we only need to find the special vertex and dont care about the others which are much harder to distinguish
c        if(pentagoncount.eq.1) then
c          do k=1,3 ! iterate over neighbour atoms
c            arbitrary_index=0
c            do l=1,hexagoncount ! iterate over neighbour hexagons
c              do m=1,6 ! iterate over atoms in hexagon
c                if (neighbour_atoms(k).eq.N6M(l,m)) then
c                  arbitrary_index=arbitrary_index+1
c                endif
c              enddo
c            enddo
c            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two hexagons
c              buffer=neighbour_atoms(k)
c              neighbour_atoms(k)=neighbour_atoms(1)
c              neighbour_atoms(1)=buffer
c            endif
c          enddo
c        endif
c        if(pentagoncount.eq.2) then
c          do k=1,3 ! iterate over neighbour atoms
c            arbitrary_index=0
c            do l=1,pentagoncount! iterate over neighbour pentagons
c              do m=1,5 ! iterate over atoms in hexagon
c                if (neighbour_atoms(k).eq.N5M(l,m)) then
c                  arbitrary_index=arbitrary_index+1
c                endif
c              enddo
c            enddo
c            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two pentagons
c              buffer=neighbour_atoms(k)
c              neighbour_atoms(k)=neighbour_atoms(1)
c              neighbour_atoms(1)=buffer
c            endif
c          enddo
c        endif
c atoms
c        J1=3*neighbour_atoms(1)
c        J2=3*neighbour_atoms(2)
c        J3=3*neighbour_atoms(3)
c        J4=3*I
c coordinates
c        call dihedral(p(J1-2),p(J1-1),p(J1),p(J2-2),p(J2-1),
c     2   p(J2),p(J3-2),p(J3-1),p(J3),p(J4-2),p(J4-1),
c     3   p(J4),angle_abcd)
c        write(*,*)i,angle_abcd,"dihedral angle (in radians)"
c        select case(pentagoncount)
c          case(0)
c            zero_value=dhhh
c          case(1)
c            zero_value=dhhp
c          case(2)
c            zero_value=dhpp
c          case(3)
c            zero_value=dppp
c        end select
c        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
c        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
c        angle_abcd=dabs(angle_abcd)
c        increment=(angle_abcd-zero_value)**2
c        select case(pentagoncount)
c          case(0)
c            ehookdhhh=ehookdhhh+increment
c          case(1)
c            ehookdhhp=ehookdhhp+increment
c          case(2)
c            ehookdhpp=ehookdhpp+increment
c          case(3)
c            ehookdppp=ehookdppp+increment
c        end select
c      enddo


C     total energy  
      fc=frpp*ehookrpp+frhp*ehookrhp+frhh*ehookrhh ! stretching
     2 +fap*ehookap+fah*ehookah ! bending
     3 +fdppp*ehookdppp+fdhpp*ehookdhpp+fdhhp*ehookdhhp+fdhhh*ehookdhhh! dihedral
c      write(*,*)fc,"energy"
      Return
      END SUBROUTINE



      SUBROUTINE dfunc3d(n,A,N5,N6,N5M,N6M,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      integer iopt
c      write(*,*)"entering dfunc3d"
      
      select case(iopt)
        case(1)
          CALL dwu(n,A,N5,N6,N5M,N6M,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p)
        case(2)
          CALL dwu(n,A,N5,N6,N5M,N6M,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p)
        case(3)
          CALL dextwu(n,A,N5,N6,N5M,N6M,p,x,force,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      end select

c      write(*,*)"leaving dfunc3d"
      return
      END SUBROUTINE


      SUBROUTINE dwu(n,A,N5,N6,N5M,N6M,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p)
c n=MATOM*3
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, gradient
      Real*8 p(nmax*3),x(nmax*3),force(ffmaxdim)
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),iopt

c     edges with 0, 1, 2 pentagons
      integer e_hh(2,N/2), e_hp(2,N/2), e_pp(2,N/2)
      integer a_p(3,60), a_h(3,N-60)
c     counter for edges with 0, 1, 2 pentagons neighbours
      integer ne_hh,ne_hp,ne_pp

c      write(*,*)'entering dwu'

      rp=force(1)
      rh=force(2)
      ap=force(3)
      ah=force(4)
      frp=force(5)
      frh=force(6)
      fap=force(7)
      fah=force(8)
      fco=force(9)

      if(ne_hh .ne. 0) then
        do i=1,ne_hh
          call ddist(p(3*e_hh(1,i)-2),p(3*e_hh(1,i)-1),p(3*e_hh(1,i)),
     1      p(3*e_hh(2,i)-2),p(3*e_hh(2,i)-1),p(3*e_hh(2,i)),
     1      dax,day,daz,dbx,dby,dbz,ratom)
          zero_value=rh
          force_constant=frh
          dE_over_dc=force_constant*(ratom-zero_value)
          x(3*e_hh(1,i)-2)=x(3*e_hh(1,i)-2) + dax*dE_over_dc
          x(3*e_hh(1,i)-1)=x(3*e_hh(1,i)-1) + day*dE_over_dc
          x(3*e_hh(1,i))  =x(3*e_hh(1,i))   + daz*dE_over_dc
          x(3*e_hh(2,i)-2)=x(3*e_hh(2,i)-2) + dbx*dE_over_dc
          x(3*e_hh(2,i)-1)=x(3*e_hh(2,i)-1) + dby*dE_over_dc
          x(3*e_hh(2,i))  =x(3*e_hh(2,i))   + dbz*dE_over_dc
        enddo
      endif


c     in the wu ff there is no difference between 1 and 2 pentagons
      if(ne_hp .ne. 0) then
        do i=1,ne_hp
          call ddist(p(3*e_hp(1,i)-2),p(3*e_hp(1,i)-1),p(3*e_hp(1,i)),
     1      p(3*e_hp(2,i)-2),p(3*e_hp(2,i)-1),p(3*e_hp(2,i)),
     1      dax,day,daz,dbx,dby,dbz,ratom)
          zero_value=rp
          force_constant=frp
          dE_over_dc=force_constant*(ratom-zero_value)
          x(3*e_hp(1,i)-2)=x(3*e_hp(1,i)-2) + dax*dE_over_dc
          x(3*e_hp(1,i)-1)=x(3*e_hp(1,i)-1) + day*dE_over_dc
          x(3*e_hp(1,i))  =x(3*e_hp(1,i))   + daz*dE_over_dc
          x(3*e_hp(2,i)-2)=x(3*e_hp(2,i)-2) + dbx*dE_over_dc
          x(3*e_hp(2,i)-1)=x(3*e_hp(2,i)-1) + dby*dE_over_dc
          x(3*e_hp(2,i))  =x(3*e_hp(2,i))   + dbz*dE_over_dc
        enddo
      endif

      if(ne_pp .ne. 0) then
        do i=1,ne_pp
          call ddist(p(3*e_pp(1,i)-2),p(3*e_pp(1,i)-1),p(3*e_pp(1,i)),
     1      p(3*e_pp(2,i)-2),p(3*e_pp(2,i)-1),p(3*e_pp(2,i)),
     1      dax,day,daz,dbx,dby,dbz,ratom)
          zero_value=rp
          force_constant=frp
          dE_over_dc=force_constant*(ratom-zero_value)
          x(3*e_pp(1,i)-2)=x(3*e_pp(1,i)-2) + dax*dE_over_dc
          x(3*e_pp(1,i)-1)=x(3*e_pp(1,i)-1) + day*dE_over_dc
          x(3*e_pp(1,i))  =x(3*e_pp(1,i))   + daz*dE_over_dc
          x(3*e_pp(2,i)-2)=x(3*e_pp(2,i)-2) + dbx*dE_over_dc
          x(3*e_pp(2,i)-1)=x(3*e_pp(2,i)-1) + dby*dE_over_dc
          x(3*e_pp(2,i))  =x(3*e_pp(2,i))   + dbz*dE_over_dc
        enddo
      endif


C     Stretching
c      Do I=1,n,3
c        ehookx=0.d0
c        ehooky=0.d0
c        ehookz=0.d0
c        I1=(I+2)/3
c        Do J=I+3,n,3
c          J1=(J+2)/3
c          if(A(I1,J1).ne.0) then
c            ax=p(i)
c            ay=p(i+1)
c            az=p(i+2)
c            bx=p(j)
c            by=p(j+1)
c            bz=p(j+2)
c            call DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,ratom)
C           Check if bond is part of 5-ring
c            pentagoncount=0
c            do IB=1,12
c              ir1=0
c              ir2=0
c              do JB=1,5
c                if(I1.eq.N5M(IB,JB)) ir1=1
c                if(J1.eq.N5M(IB,JB)) ir2=1
c              enddo
c              if(ir1.eq.1.and.ir2.eq.1) then
C               5-ring
c                pentagoncount=pentagoncount+1
c                go to 1
c              endif
c            enddo
C           6-ring
c 1          if (pentagoncount.eq.0)then
c              zero_value=rh
c              force_constant=frh
c            else
c              zero_value=rp
c              force_constant=frp
c            end if
c            dE_over_dc=force_constant*(ratom-zero_value)
c            x(i  )=x(i  )+dax*dE_over_dc
c            x(i+1)=x(i+1)+day*dE_over_dc
c            x(i+2)=x(i+2)+daz*dE_over_dc
c            x(j  )=x(j  )+dbx*dE_over_dc
c            x(j+1)=x(j+1)+dby*dE_over_dc
c            x(j+2)=x(j+2)+dbz*dE_over_dc
c          endif
c        enddo
c      enddo
        
C     Bending
C     Loop over 5-rings
      do i=1,60
        ax=p(3*a_p(1,i)-2)
        ay=p(3*a_p(1,i)-1)
        az=p(3*a_p(1,i))
        bx=p(3*a_p(2,i)-2)
        by=p(3*a_p(2,i)-1)
        bz=p(3*a_p(2,i))
        cx=p(3*a_p(3,i)-2)
        cy=p(3*a_p(3,i)-1)
        cz=p(3*a_p(3,i))  
        call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3   angle_abc)
        zero_value=ap
        force_constant=fap
        dE_over_dc=force_constant*(angle_abc-zero_value)
        x(3*a_p(1,i)-2)=x(3*a_p(1,i)-2)+dax*dE_over_dc
        x(3*a_p(1,i)-1)=x(3*a_p(1,i)-1)+day*dE_over_dc
        x(3*a_p(1,i))  =x(3*a_p(1,i))  +daz*dE_over_dc
        x(3*a_p(2,i)-2)=x(3*a_p(2,i)-2)+dbx*dE_over_dc
        x(3*a_p(2,i)-1)=x(3*a_p(2,i)-1)+dby*dE_over_dc
        x(3*a_p(2,i))  =x(3*a_p(2,i))  +dbz*dE_over_dc
        x(3*a_p(3,i)-2)=x(3*a_p(3,i)-2)+dcx*dE_over_dc
        x(3*a_p(3,i)-1)=x(3*a_p(3,i)-1)+dcy*dE_over_dc
        x(3*a_p(3,i))  =x(3*a_p(3,i))  +dcz*dE_over_dc
      enddo

c      Do I=1,N5 !=12
c        Do J=1,5
c          JLX=J-1
c          JRX=J+1
c          if(JLX.eq.0) JLX=5
c          if(JRX.eq.6) JRX=1
c          JL=3*N5M(I,JLX)
c          JM=3*N5M(I,J)
c          JR=3*N5M(I,JRX)
c          ax=p(JL-2)
c          ay=p(JL-1)
c          az=p(JL)
c          bx=p(JM-2)
c          by=p(JM-1)
c          bz=p(JM)
c          cx=p(JR-2)
c          cy=p(JR-1)
c          cz=p(JR)
c          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
c     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
c     3     angle_abc)
c          zero_value=ap
c          force_constant=fap
c          dE_over_dc=force_constant*(angle_abc-zero_value)
c          x(JL-2)=x(JL-2)+dax*dE_over_dc
c          x(JL-1)=x(JL-1)+day*dE_over_dc
c          x(JL)  =x(JL)  +daz*dE_over_dc
c          x(JM-2)=x(JM-2)+dbx*dE_over_dc
c          x(JM-1)=x(JM-1)+dby*dE_over_dc
c          x(JM)  =x(JM)  +dbz*dE_over_dc
c          x(JR-2)=x(JR-2)+dcx*dE_over_dc
c          x(JR-1)=x(JR-1)+dcy*dE_over_dc
c          x(JR)  =x(JR)  +dcz*dE_over_dc
c        enddo
c      enddo
      
C     Loop over 6-rings
      if(n/3 .gt. 20) then
        do i=1,n-60
          ax=p(3*a_h(1,i)-2)
          ay=p(3*a_h(1,i)-1)
          az=p(3*a_h(1,i))
          bx=p(3*a_h(2,i)-2)
          by=p(3*a_h(2,i)-1)
          bz=p(3*a_h(2,i))
          cx=p(3*a_h(3,i)-2)
          cy=p(3*a_h(3,i)-1)
          cz=p(3*a_h(3,i))  
          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3     angle_abc)
          zero_value=ah
          force_constant=fah
          dE_over_dc=force_constant*(angle_abc-zero_value)
          x(3*a_h(1,i)-2)=x(3*a_h(1,i)-2)+dax*dE_over_dc
          x(3*a_h(1,i)-1)=x(3*a_h(1,i)-1)+day*dE_over_dc
          x(3*a_h(1,i))  =x(3*a_h(1,i))  +daz*dE_over_dc
          x(3*a_h(2,i)-2)=x(3*a_h(2,i)-2)+dbx*dE_over_dc
          x(3*a_h(2,i)-1)=x(3*a_h(2,i)-1)+dby*dE_over_dc
          x(3*a_h(2,i))  =x(3*a_h(2,i))  +dbz*dE_over_dc
          x(3*a_h(3,i)-2)=x(3*a_h(3,i)-2)+dcx*dE_over_dc
          x(3*a_h(3,i)-1)=x(3*a_h(3,i)-1)+dcy*dE_over_dc
          x(3*a_h(3,i))  =x(3*a_h(3,i))  +dcz*dE_over_dc
        enddo
      endif

c      Do I=1,N6
c        Do J=1,6
c          JLX=J-1
c          JRX=J+1
c          if(JLX.eq.0) JLX=6
c          if(JRX.eq.7) JRX=1
c          JL=3*N6M(I,JLX)
c          JM=3*N6M(I,J)
c          JR=3*N6M(I,JRX)
c          ax=p(JL-2)
c          ay=p(JL-1)
c          az=p(JL)
c          bx=p(JM-2)
c          by=p(JM-1)
c          bz=p(JM)
c          cx=p(JR-2)
c          cy=p(JR-1)
c          cz=p(JR)
c          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
c     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
c     3     angle_abc)
c          zero_value=ah
c          force_constant=fah
c          dE_over_dc=force_constant*(angle_abc-zero_value)
c      write(*,*)"a,z,de",angle_abc,zero_value,dE_over_dc
c      write(*,*)"dd",dax*dE_over_dc,day*dE_over_dc,
c     2  daz*dE_over_dc,
c     2  dbx*dE_over_dc,dby*dE_over_dc,dbz*dE_over_dc,
c     2  dcx*dE_over_dc,dcy*dE_over_dc,dcz*dE_over_dc
c          x(JL-2)=x(JL-2)+dax*dE_over_dc
c          x(JL-1)=x(JL-1)+day*dE_over_dc
c          x(JL)  =x(JL)  +daz*dE_over_dc
c          x(JM-2)=x(JM-2)+dbx*dE_over_dc
c          x(JM-1)=x(JM-1)+dby*dE_over_dc
c          x(JM)  =x(JM)  +dbz*dE_over_dc
c          x(JR-2)=x(JR-2)+dcx*dE_over_dc
c          x(JR-1)=x(JR-1)+dcy*dE_over_dc
c          x(JR)  =x(JR)  +dcz*dE_over_dc
c        enddo
c      enddo
c      write(*,*)x

C     Coulomb repulsion from origin
      if (iopt.eq.2 .and. fco.ne.0.d0)  then
        Do I=1,n,3
          rinv=(p(I)**2+p(I+1)**2+p(I+2)**2)**(-1.5d0)
          x(I)  =x(I)  +fco*rinv*p(I)
          x(I+1)=x(I+1)+fco*rinv*p(I+1)
          x(I+2)=x(I+2)+fco*rinv*p(I+2)
        enddo
      endif

c      write(*,*)'leaving dwu'
      return
      END


      SUBROUTINE dextwu(n,A,N5,N6,N5M,N6M,p,x,force,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 p(nmax*3),x(nmax*3),force(ffmaxdim)
c      real*8 J1x,J1y,J1z,J2x,J2y,J2z,J3x,J3y,J3z,J4x,J4y,J4z
      Integer A(NMAX,NMAX),N5M(MMAX,5),N6M(MMAX,6),pentagoncount,
     2 hexagoncount,arbitrary_index,neighbour_atoms(3),
     3 neighbour_faces_h(3),neighbour_faces_p(3),buffer

c     edges with 0, 1, 2 pentagons
      integer e_hh(2,N/2), e_hp(2,N/2), e_pp(2,N/2)
      integer a_p(3,60), a_h(3,N-60)
      integer d_hhh(4,n/3),d_hpp(4,n/3),d_hhp(4,n/3),d_ppp(4,n/3)
c     counter for edges with 0, 1, 2 pentagons neighbours
      integer ne_hh,ne_hp,ne_pp
c     counter for dihedrals with 0, 1, 2, 3 pentagons neighbours
      integer nd_hhh,nd_hhp,nd_hpp,nd_ppp

c      write(*,*)'entering dextwu'
c      write(*,*)'n',n

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
c     we distinguish between bonds between two hexagons, two pentagons and hex/pent
      
c      write(*,*)'before stretch'
      if(ne_hh .ne. 0) then
        do i=1,ne_hh
          call ddist(p(3*e_hh(1,i)-2),p(3*e_hh(1,i)-1),p(3*e_hh(1,i)),
     1      p(3*e_hh(2,i)-2),p(3*e_hh(2,i)-1),p(3*e_hh(2,i)),
     1      dax,day,daz,dbx,dby,dbz,ratom)
          zero_value=rhh
          force_constant=frhh
          dE_over_dc=force_constant*(ratom-zero_value)
          x(3*e_hh(1,i)-2)=x(3*e_hh(1,i)-2) + dax*dE_over_dc
          x(3*e_hh(1,i)-1)=x(3*e_hh(1,i)-1) + day*dE_over_dc
          x(3*e_hh(1,i))  =x(3*e_hh(1,i))   + daz*dE_over_dc
          x(3*e_hh(2,i)-2)=x(3*e_hh(2,i)-2) + dbx*dE_over_dc
          x(3*e_hh(2,i)-1)=x(3*e_hh(2,i)-1) + dby*dE_over_dc
          x(3*e_hh(2,i))  =x(3*e_hh(2,i))   + dbz*dE_over_dc
        enddo
      endif

      if(ne_hp .ne. 0) then
        do i=1,ne_hp
          call ddist(p(3*e_hp(1,i)-2),p(3*e_hp(1,i)-1),p(3*e_hp(1,i)),
     1      p(3*e_hp(2,i)-2),p(3*e_hp(2,i)-1),p(3*e_hp(2,i)),
     1      dax,day,daz,dbx,dby,dbz,ratom)
          zero_value=rhp
          force_constant=frhp
          dE_over_dc=force_constant*(ratom-zero_value)
          x(3*e_hp(1,i)-2)=x(3*e_hp(1,i)-2) + dax*dE_over_dc
          x(3*e_hp(1,i)-1)=x(3*e_hp(1,i)-1) + day*dE_over_dc
          x(3*e_hp(1,i))  =x(3*e_hp(1,i))   + daz*dE_over_dc
          x(3*e_hp(2,i)-2)=x(3*e_hp(2,i)-2) + dbx*dE_over_dc
          x(3*e_hp(2,i)-1)=x(3*e_hp(2,i)-1) + dby*dE_over_dc
          x(3*e_hp(2,i))  =x(3*e_hp(2,i))   + dbz*dE_over_dc
        enddo
      endif

      if(ne_pp .ne. 0) then
        do i=1,ne_pp
          call ddist(p(3*e_pp(1,i)-2),p(3*e_pp(1,i)-1),p(3*e_pp(1,i)),
     1      p(3*e_pp(2,i)-2),p(3*e_pp(2,i)-1),p(3*e_pp(2,i)),
     1      dax,day,daz,dbx,dby,dbz,ratom)
          zero_value=rpp
          force_constant=frpp
          dE_over_dc=force_constant*(ratom-zero_value)
          x(3*e_pp(1,i)-2)=x(3*e_pp(1,i)-2) + dax*dE_over_dc
          x(3*e_pp(1,i)-1)=x(3*e_pp(1,i)-1) + day*dE_over_dc
          x(3*e_pp(1,i))  =x(3*e_pp(1,i))   + daz*dE_over_dc
          x(3*e_pp(2,i)-2)=x(3*e_pp(2,i)-2) + dbx*dE_over_dc
          x(3*e_pp(2,i)-1)=x(3*e_pp(2,i)-1) + dby*dE_over_dc
          x(3*e_pp(2,i))  =x(3*e_pp(2,i))   + dbz*dE_over_dc
        enddo
      endif




c      Do I=1,n,3
c        i1=(i+2)/3
c        Do J=I+3,n,3
c          j1=(j+2)/3
c check if bond exists
c          if(A(i1,j1).ne.0) then
c get coordinates
c            ax=p(I)
c            ay=p(I+1)
c            az=p(I+2)
c            bx=p(J)
c            by=p(J+1)
c            bz=p(J+2)
c            call DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,ratom)
C           Check if bond is part of 5-ring
c            pentagoncount=0
c            do IB=1,12
c              ir1=0
c              ir2=0
c              do JB=1,5
c                if(i1.eq.N5M(IB,JB)) ir1=1
c                if(j1.eq.N5M(IB,JB)) ir2=1
c              enddo
c              if(ir1.eq.1.and.ir2.eq.1) then
c                pentagoncount=pentagoncount+1
c              endif
c            enddo
c            select case(pentagoncount)
c            case(0)
c              zero_value=rhh
c              force_constant=frhh
c            case(1)
c              zero_value=rhp
c              force_constant=frhp
c            case(2)
c              zero_value=rpp
c              force_constant=frpp
c            end select
c            dE_over_dc=force_constant*(ratom-zero_value)
c            x(I)  =x(I)  +dax*dE_over_dc
c            x(I+1)=x(I+1)+day*dE_over_dc
c            x(I+2)=x(I+2)+daz*dE_over_dc
c            x(J)  =x(J)  +dbx*dE_over_dc
c            x(J+1)=x(J+1)+dby*dE_over_dc
c            x(J+2)=x(J+2)+dbz*dE_over_dc
c          endif
c        enddo
c      enddo

c      enddo
c      write(*,*)'before bend-p'
        
C     Bending
C     Loop over 5-rings
      do i=1,60
        ax=p(3*a_p(1,i)-2)
        ay=p(3*a_p(1,i)-1)
        az=p(3*a_p(1,i))
        bx=p(3*a_p(2,i)-2)
        by=p(3*a_p(2,i)-1)
        bz=p(3*a_p(2,i))
        cx=p(3*a_p(3,i)-2)
        cy=p(3*a_p(3,i)-1)
        cz=p(3*a_p(3,i))  
        call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3   angle_abc)
        zero_value=ap
        force_constant=fap
        dE_over_dc=force_constant*(angle_abc-zero_value)
        x(3*a_p(1,i)-2)=x(3*a_p(1,i)-2)+dax*dE_over_dc
        x(3*a_p(1,i)-1)=x(3*a_p(1,i)-1)+day*dE_over_dc
        x(3*a_p(1,i))  =x(3*a_p(1,i))  +daz*dE_over_dc
        x(3*a_p(2,i)-2)=x(3*a_p(2,i)-2)+dbx*dE_over_dc
        x(3*a_p(2,i)-1)=x(3*a_p(2,i)-1)+dby*dE_over_dc
        x(3*a_p(2,i))  =x(3*a_p(2,i))  +dbz*dE_over_dc
        x(3*a_p(3,i)-2)=x(3*a_p(3,i)-2)+dcx*dE_over_dc
        x(3*a_p(3,i)-1)=x(3*a_p(3,i)-1)+dcy*dE_over_dc
        x(3*a_p(3,i))  =x(3*a_p(3,i))  +dcz*dE_over_dc
      enddo

c      write(*,*)'before bend-h'
C     Loop over 6-rings
      if(n/3 .gt. 20) then
        do i=1,n-60
          ax=p(3*a_h(1,i)-2)
          ay=p(3*a_h(1,i)-1)
          az=p(3*a_h(1,i))
          bx=p(3*a_h(2,i)-2)
          by=p(3*a_h(2,i)-1)
          bz=p(3*a_h(2,i))
          cx=p(3*a_h(3,i)-2)
          cy=p(3*a_h(3,i)-1)
          cz=p(3*a_h(3,i))  
          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3     angle_abc)
          zero_value=ah
          force_constant=fah
          dE_over_dc=force_constant*(angle_abc-zero_value)
          x(3*a_h(1,i)-2)=x(3*a_h(1,i)-2)+dax*dE_over_dc
          x(3*a_h(1,i)-1)=x(3*a_h(1,i)-1)+day*dE_over_dc
          x(3*a_h(1,i))  =x(3*a_h(1,i))  +daz*dE_over_dc
          x(3*a_h(2,i)-2)=x(3*a_h(2,i)-2)+dbx*dE_over_dc
          x(3*a_h(2,i)-1)=x(3*a_h(2,i)-1)+dby*dE_over_dc
          x(3*a_h(2,i))  =x(3*a_h(2,i))  +dbz*dE_over_dc
          x(3*a_h(3,i)-2)=x(3*a_h(3,i)-2)+dcx*dE_over_dc
          x(3*a_h(3,i)-1)=x(3*a_h(3,i)-1)+dcy*dE_over_dc
          x(3*a_h(3,i))  =x(3*a_h(3,i))  +dcz*dE_over_dc
        enddo
      endif

C     Bending        
C     Loop over 5-rings
c      Do I=1,N5 !and n5==12
c        Do J=1,5
c          JLX=J-1
c          JRX=J+1
c          if(JLX.eq.0) JLX=5
c          if(JRX.eq.6) JRX=1
c          JL=3*N5M(I,JLX)
cc          JM=3*N5M(I,J)
c          JR=3*N5M(I,JRX)
c          ax=p(JL-2)
c          ay=p(JL-1)
c          az=p(JL)
c          bx=p(JM-2)
c          by=p(JM-1)
c          bz=p(JM)
c          cx=p(JR-2)
c          cy=p(JR-1)
c          cz=p(JR)
c          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
c     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
c     3     angle_abc)
c          call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
c          zero_value=ap
c          force_constant=fap
c          dE_over_dc=force_constant*(angle_abc-zero_value)
c          x(JL-2)=x(JL-2)+dax*dE_over_dc
c          x(JL-1)=x(JL-1)+day*dE_over_dc
c          x(JL)  =x(JL)  +daz*dE_over_dc
c          x(JM-2)=x(JM-2)+dbx*dE_over_dc
c          x(JM-1)=x(JM-1)+dby*dE_over_dc
c          x(JM)  =x(JM)  +dbz*dE_over_dc
c          x(JR-2)=x(JR-2)+dcx*dE_over_dc
c          x(JR-1)=x(JR-1)+dcy*dE_over_dc
c          x(JR)  =x(JR)  +dcz*dE_over_dc
c        enddo
c      enddo
      
C     Loop over 6-rings
c      Do I=1,N6
c        Do J=1,6
c          JLX=J-1
c          JRX=J+1
c          if(JLX.eq.0) JLX=6
c          if(JRX.eq.7) JRX=1
c          JL=3*N6M(I,JLX)
c          JM=3*N6M(I,J)
c          JR=3*N6M(I,JRX)
c          ax=p(JL-2)
c          ay=p(JL-1)
c          az=p(JL)
c          bx=p(JM-2)
c          by=p(JM-1)
c          bz=p(JM)
c          cx=p(JR-2)
c          cy=p(JR-1)
c          cz=p(JR)
c          call DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
c     2     dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
c     3     angle_abc)
c          zero_value=ah
c          force_constant=fah
c          dE_over_dc=force_constant*(angle_abc-zero_value)
c          x(JL-2)=x(JL-2)+dax*dE_over_dc
c          x(JL-1)=x(JL-1)+day*dE_over_dc
c          x(JL)  =x(JL)  +daz*dE_over_dc
c          x(JM-2)=x(JM-2)+dbx*dE_over_dc
c          x(JM-1)=x(JM-1)+dby*dE_over_dc
c          x(JM)  =x(JM)  +dbz*dE_over_dc
c          x(JR-2)=x(JR-2)+dcx*dE_over_dc
c          x(JR-1)=x(JR-1)+dcy*dE_over_dc
c          x(JR)  =x(JR)  +dcz*dE_over_dc
c        enddo
c      enddo

C     dihedrals 
c     3 hexagons
      if (nd_hhh .ne. 0) then
      do i=1,nd_hhh
        ax=p(3*d_hhh(1,i)-2)
        ay=p(3*d_hhh(1,i)-1)
        az=p(3*d_hhh(1,i))
        bx=p(3*d_hhh(2,i)-2)
        by=p(3*d_hhh(2,i)-1)
        bz=p(3*d_hhh(2,i))
        cx=p(3*d_hhh(3,i)-2)
        cy=p(3*d_hhh(3,i)-1)
        cz=p(3*d_hhh(3,i))  
        dx=p(3*d_hhh(4,i)-2)
        dy=p(3*d_hhh(4,i)-1)
        dz=p(3*d_hhh(4,i))  
        call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        zero_value=dhhh
        force_constant=fdhhh
        dE_over_dc=force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(3*d_hhh(1,i)-2)=x(3*d_hhh(1,i)-2)+dax*dE_over_dc
        x(3*d_hhh(1,i)-1)=x(3*d_hhh(1,i)-1)+day*dE_over_dc
        x(3*d_hhh(1,i))  =x(3*d_hhh(1,i))  +daz*dE_over_dc
        x(3*d_hhh(2,i)-2)=x(3*d_hhh(2,i)-2)+dbx*dE_over_dc
        x(3*d_hhh(2,i)-1)=x(3*d_hhh(2,i)-1)+dby*dE_over_dc
        x(3*d_hhh(2,i))  =x(3*d_hhh(2,i))  +dbz*dE_over_dc
        x(3*d_hhh(3,i)-2)=x(3*d_hhh(3,i)-2)+dcx*dE_over_dc
        x(3*d_hhh(3,i)-1)=x(3*d_hhh(3,i)-1)+dcy*dE_over_dc
        x(3*d_hhh(3,i))  =x(3*d_hhh(3,i))  +dcz*dE_over_dc
        x(3*d_hhh(4,i)-2)=x(3*d_hhh(4,i)-2)+ddx*dE_over_dc
        x(3*d_hhh(4,i)-1)=x(3*d_hhh(4,i)-1)+ddy*dE_over_dc
        x(3*d_hhh(4,i))  =x(3*d_hhh(4,i))  +ddz*dE_over_dc
      enddo
      endif

c     2 hexagon, 1 pentagon
      if (nd_hhp .ne. 0) then
      do i=1,nd_hhp
        ax=p(3*d_hhp(1,i)-2)
        ay=p(3*d_hhp(1,i)-1)
        az=p(3*d_hhp(1,i))
        bx=p(3*d_hhp(2,i)-2)
        by=p(3*d_hhp(2,i)-1)
        bz=p(3*d_hhp(2,i))
        cx=p(3*d_hhp(3,i)-2)
        cy=p(3*d_hhp(3,i)-1)
        cz=p(3*d_hhp(3,i))  
        dx=p(3*d_hhp(4,i)-2)
        dy=p(3*d_hhp(4,i)-1)
        dz=p(3*d_hhp(4,i))  
        call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        zero_value=dhhp
        force_constant=fdhhp
        dE_over_dc=force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(3*d_hhp(1,i)-2)=x(3*d_hhp(1,i)-2)+dax*dE_over_dc
        x(3*d_hhp(1,i)-1)=x(3*d_hhp(1,i)-1)+day*dE_over_dc
        x(3*d_hhp(1,i))  =x(3*d_hhp(1,i))  +daz*dE_over_dc
        x(3*d_hhp(2,i)-2)=x(3*d_hhp(2,i)-2)+dbx*dE_over_dc
        x(3*d_hhp(2,i)-1)=x(3*d_hhp(2,i)-1)+dby*dE_over_dc
        x(3*d_hhp(2,i))  =x(3*d_hhp(2,i))  +dbz*dE_over_dc
        x(3*d_hhp(3,i)-2)=x(3*d_hhp(3,i)-2)+dcx*dE_over_dc
        x(3*d_hhp(3,i)-1)=x(3*d_hhp(3,i)-1)+dcy*dE_over_dc
        x(3*d_hhp(3,i))  =x(3*d_hhp(3,i))  +dcz*dE_over_dc
        x(3*d_hhp(4,i)-2)=x(3*d_hhp(4,i)-2)+ddx*dE_over_dc
        x(3*d_hhp(4,i)-1)=x(3*d_hhp(4,i)-1)+ddy*dE_over_dc
        x(3*d_hhp(4,i))  =x(3*d_hhp(4,i))  +ddz*dE_over_dc
      enddo
      endif

c     1 hexagon, 2 pentagons
      if (nd_hpp .ne. 0) then
      do i=1,nd_hpp
        ax=p(3*d_hpp(1,i)-2)
        ay=p(3*d_hpp(1,i)-1)
        az=p(3*d_hpp(1,i))
        bx=p(3*d_hpp(2,i)-2)
        by=p(3*d_hpp(2,i)-1)
        bz=p(3*d_hpp(2,i))
        cx=p(3*d_hpp(3,i)-2)
        cy=p(3*d_hpp(3,i)-1)
        cz=p(3*d_hpp(3,i))  
        dx=p(3*d_hpp(4,i)-2)
        dy=p(3*d_hpp(4,i)-1)
        dz=p(3*d_hpp(4,i))  
        call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        zero_value=dhpp
        force_constant=fdhpp
        dE_over_dc=force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(3*d_hpp(1,i)-2)=x(3*d_hpp(1,i)-2)+dax*dE_over_dc
        x(3*d_hpp(1,i)-1)=x(3*d_hpp(1,i)-1)+day*dE_over_dc
        x(3*d_hpp(1,i))  =x(3*d_hpp(1,i))  +daz*dE_over_dc
        x(3*d_hpp(2,i)-2)=x(3*d_hpp(2,i)-2)+dbx*dE_over_dc
        x(3*d_hpp(2,i)-1)=x(3*d_hpp(2,i)-1)+dby*dE_over_dc
        x(3*d_hpp(2,i))  =x(3*d_hpp(2,i))  +dbz*dE_over_dc
        x(3*d_hpp(3,i)-2)=x(3*d_hpp(3,i)-2)+dcx*dE_over_dc
        x(3*d_hpp(3,i)-1)=x(3*d_hpp(3,i)-1)+dcy*dE_over_dc
        x(3*d_hpp(3,i))  =x(3*d_hpp(3,i))  +dcz*dE_over_dc
        x(3*d_hpp(4,i)-2)=x(3*d_hpp(4,i)-2)+ddx*dE_over_dc
        x(3*d_hpp(4,i)-1)=x(3*d_hpp(4,i)-1)+ddy*dE_over_dc
        x(3*d_hpp(4,i))  =x(3*d_hpp(4,i))  +ddz*dE_over_dc
      enddo
      endif

c     3 pentagons
      if (nd_ppp .ne. 0) then
      do i=1,nd_ppp
        ax=p(3*d_ppp(1,i)-2)
        ay=p(3*d_ppp(1,i)-1)
        az=p(3*d_ppp(1,i))
        bx=p(3*d_ppp(2,i)-2)
        by=p(3*d_ppp(2,i)-1)
        bz=p(3*d_ppp(2,i))
        cx=p(3*d_ppp(3,i)-2)
        cy=p(3*d_ppp(3,i)-1)
        cz=p(3*d_ppp(3,i))  
        dx=p(3*d_ppp(4,i)-2)
        dy=p(3*d_ppp(4,i)-1)
        dz=p(3*d_ppp(4,i))  
        call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        zero_value=dppp
        force_constant=fdppp
        dE_over_dc=force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(3*d_ppp(1,i)-2)=x(3*d_ppp(1,i)-2)+dax*dE_over_dc
        x(3*d_ppp(1,i)-1)=x(3*d_ppp(1,i)-1)+day*dE_over_dc
        x(3*d_ppp(1,i))  =x(3*d_ppp(1,i))  +daz*dE_over_dc
        x(3*d_ppp(2,i)-2)=x(3*d_ppp(2,i)-2)+dbx*dE_over_dc
        x(3*d_ppp(2,i)-1)=x(3*d_ppp(2,i)-1)+dby*dE_over_dc
        x(3*d_ppp(2,i))  =x(3*d_ppp(2,i))  +dbz*dE_over_dc
        x(3*d_ppp(3,i)-2)=x(3*d_ppp(3,i)-2)+dcx*dE_over_dc
        x(3*d_ppp(3,i)-1)=x(3*d_ppp(3,i)-1)+dcy*dE_over_dc
        x(3*d_ppp(3,i))  =x(3*d_ppp(3,i))  +dcz*dE_over_dc
        x(3*d_ppp(4,i)-2)=x(3*d_ppp(4,i)-2)+ddx*dE_over_dc
        x(3*d_ppp(4,i)-1)=x(3*d_ppp(4,i)-1)+ddy*dE_over_dc
        x(3*d_ppp(4,i))  =x(3*d_ppp(4,i))  +ddz*dE_over_dc
      enddo
      endif



c      Do I=1,n/3 ! iterate over atoms
c classify vertex acording to adjacent faces
c find neighbouring faces (3)
c        pentagoncount=0
c        do IB=1,12 ! iterate over pentagons
c          do JB=1,5 ! iterate over atoms in pentagon
c            if(I.eq.N5M(IB,JB)) then
c              pentagoncount=pentagoncount+1 ! find if I is part of 0,...,3 pentagons
c              neighbour_faces_p(pentagoncount)=IB
c            endif
c          enddo
c        enddo
c        hexagoncount=0
c        do IB=1,n/3-12 ! because n=vertex_count * 3 ()
c          do JB=1,6
c            if(i.eq.N6M(IB,JB)) then
c              hexagoncount=hexagoncount+1
c              neighbour_faces_h(hexagoncount)=IB
c            endif
c          enddo
c        enddo
c find neighbouring atoms (3)
c        neighbour_atom_count=1
c        do j=1,n/3
c          if(A(I,j).ne.0) then
c            neighbour_atoms(neighbour_atom_count)=j
c            neighbour_atom_count=neighbour_atom_count+1
c          endif
c        enddo
c        write(*,*)"p:",pentagoncount, "h:",hexagoncount
c sort neighbours
c we make use of the fact, that for any dihedral ijkl=-ikjl
c therefore we only need to find the special vertex and dont care about the others which are much harder to distinguish
c        if(pentagoncount.eq.1) then
c          do k=1,3 ! iterate over neighbour atoms
c            arbitrary_index=0
c            do l=1,hexagoncount ! iterate over neighbour hexagons
c              do m=1,6 ! iterate over atoms in hexagon
c                if (neighbour_atoms(k).eq.N6M(l,m)) then
c                  arbitrary_index=arbitrary_index+1
c                endif
c              enddo
c            enddo
c            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two hexagons
c              buffer=neighbour_atoms(k)
c              neighbour_atoms(k)=neighbour_atoms(1)
c              neighbour_atoms(1)=buffer
c            endif
c          enddo
c        else if(pentagoncount.eq.2) then
c          do k=1,3 ! iterate over neighbour atoms
c            arbitrary_index=0
c            do l=1,pentagoncount ! iterate over neighbour pentagons
c              do m=1,5 ! iterate over atoms in hexagon
c                if (neighbour_atoms(k).eq.N5M(l,m)) then
c                  arbitrary_index=arbitrary_index+1
c                endif
c              enddo
c            enddo
c            if(arbitrary_index.eq.2) then ! we found the atom that lies between the two pentagons
c              buffer=neighbour_atoms(k)
c              neighbour_atoms(k)=neighbour_atoms(1)
c              neighbour_atoms(1)=buffer
c            endif
c          enddo
c        endif
c atoms
c        J1=3*neighbour_atoms(1)
c        J2=3*neighbour_atoms(2)
c        J3=3*neighbour_atoms(3)
c        J4=3*I
c        write(*,*)j1,j2,j3,j4
c coordinates
c        ax=p(J1-2)
c        ay=p(J1-1)
c        az=p(J1)
c        bx=p(J2-2)
c        by=p(J2-1)
c        bz=p(J2)
c        cx=p(J3-2)
c        cy=p(J3-1)
c        cz=p(J3)
c        dx=p(J4-2)
c        dy=p(J4-1)
c        dz=p(J4)
c        call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
c     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
c     3   angle_abcd)
c        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
c        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
c        angle_abcd=dabs(angle_abcd)
c        select case(pentagoncount)
c          case(0)
c            zero_value=dhhh
c            force_constant=fdhhh
c          case(1)
c            zero_value=dhhp
c            force_constant=fdhhp
c          case(2)
c            zero_value=dhpp
c            force_constant=fdhpp
c          case(3)
c            zero_value=dppp
c            force_constant=fdppp
c        end select
c        dE_over_dc=force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
c        x(J1-2)=x(J1-2)+dax*dE_over_dc
c        x(J1-1)=x(J1-1)+day*dE_over_dc
c        x(J1)  =x(J1)  +daz*dE_over_dc
c        x(J2-2)=x(J2-2)+dbx*dE_over_dc
c        x(J2-1)=x(J2-1)+dby*dE_over_dc
c        x(J2)  =x(J2)  +dbz*dE_over_dc
c        x(J3-2)=x(J3-2)+dcx*dE_over_dc
c        x(J3-1)=x(J3-1)+dcy*dE_over_dc
c        x(J3)  =x(J3)  +dcz*dE_over_dc
c        x(J4-2)=x(J4-2)+ddx*dE_over_dc
c        x(J4-1)=x(J4-1)+ddy*dE_over_dc
c        x(J4)  =x(J4)  +ddz*dE_over_dc
c      enddo
c      write(*,*)"d,0: ",angle_abcd,zero_value," (should be similar)"
c      write(*,*)'leaving dextwu'
      return
      END
