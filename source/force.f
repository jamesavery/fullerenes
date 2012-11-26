      SUBROUTINE func3d(n,IERR,p,fc,force,iopt,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p,
     1  d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
c n=MATOM*3
      use config
      integer iopt

      select case(iopt)
        case(1, 2)
          CALL wu(n,IERR,p,fc,force,iopt,
     1      e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1      a_h,a_p)
        case(3, 4)
          CALL extwu(n,IERR,p,fc,force,iopt,
     1      e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     2      a_h,a_p,
     1      d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        case default
          IERR=1
      end select

      return
      END SUBROUTINE


      SUBROUTINE wu(n,IERR,p,fc,force,iopt,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p)
c n=MATOM*3
      use config
      implicit real(8) (A-H,O-Z)

C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, energy
      real(8) p(nmax*3),force(ffmaxdim)
      real(8) fc
      integer iopt

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

C     Bending
      ehookap=0.d0
      ehookah=0.d0
C     Loop over 5-rings
      do i=1,60
        call angle(p(3*a_p(1,i)-2),p(3*a_p(1,i)-1),p(3*a_p(1,i)),
     1             p(3*a_p(2,i)-2),p(3*a_p(2,i)-1),p(3*a_p(2,i)),
     1             p(3*a_p(3,i)-2),p(3*a_p(3,i)-1),p(3*a_p(3,i)),
     1             angle_p)
        ehookap=ehookap+(angle_p-ap)**2
      enddo

      if(n/3 .gt. 20) then
C     Loop over 6-rings
      do i=1,n-60
        call angle(p(3*a_h(1,i)-2),p(3*a_h(1,i)-1),p(3*a_h(1,i)),
     1             p(3*a_h(2,i)-2),p(3*a_h(2,i)-1),p(3*a_h(2,i)),
     1             p(3*a_h(3,i)-2),p(3*a_h(3,i)-1),p(3*a_h(3,i)),
     1             angle_h)
        ehookah=ehookah+(angle_h-ah)**2
      enddo
      endif

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
      fc=0.5d0 * fc
c      write(*,*)'leaving wu'
      Return
      END SUBROUTINE



      SUBROUTINE extwu(n,IERR,p,fc,force,iopt,
     1  e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1  a_h,a_p,
     1  d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
c n=MATOM*3
      use config
      implicit real(8) (A-H,O-Z)
      real(8) p(nmax*3),force(ffmaxdim)
      real(8) fc

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
      fco=force(19)

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

C     Bending
c     we distinguish between angles of pentagons and hexagons
C     Loop over 5-rings
      ehookap=0.d0
      ehookah=0.d0
      do i=1,60
        call angle(p(3*a_p(1,i)-2),p(3*a_p(1,i)-1),p(3*a_p(1,i)),
     1             p(3*a_p(2,i)-2),p(3*a_p(2,i)-1),p(3*a_p(2,i)),
     1             p(3*a_p(3,i)-2),p(3*a_p(3,i)-1),p(3*a_p(3,i)),
     1             angle_abc)
        ehookap=ehookap+(angle_abc-ap)**2
      enddo

      if(n/3 .gt. 20) then
C     Loop over 6-rings
      do i=1,n-60
        call angle(p(3*a_h(1,i)-2),p(3*a_h(1,i)-1),p(3*a_h(1,i)),
     1             p(3*a_h(2,i)-2),p(3*a_h(2,i)-1),p(3*a_h(2,i)),
     1             p(3*a_h(3,i)-2),p(3*a_h(3,i)-1),p(3*a_h(3,i)),
     1             angle_abc)
        ehookah=ehookah+(angle_abc-ah)**2
      enddo
      endif


C dihedrals 
      ehookdppp=0.d0
      ehookdhpp=0.d0
      ehookdhhp=0.d0
      ehookdhhh=0.d0

      if(nd_hhh .ne. 0) then
C     3 hexagons
      do i=1,nd_hhh
      call dihedral(p(3*d_hhh(1,i)-2),p(3*d_hhh(1,i)-1),p(3*d_hhh(1,i)),
     1              p(3*d_hhh(2,i)-2),p(3*d_hhh(2,i)-1),p(3*d_hhh(2,i)),
     1              p(3*d_hhh(3,i)-2),p(3*d_hhh(3,i)-1),p(3*d_hhh(3,i)),
     1              p(3*d_hhh(4,i)-2),p(3*d_hhh(4,i)-1),p(3*d_hhh(4,i)),
     1              angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
c        write(*,*)angle_abcd
        ehookdhhh=ehookdhhh+(angle_abcd-dhhh)**2
      enddo
      endif

      if(nd_hhp .ne. 0) then
C     2 hexagons, 1 pentagon
      do i=1,nd_hhp
      call dihedral(p(3*d_hhp(1,i)-2),p(3*d_hhp(1,i)-1),p(3*d_hhp(1,i)),
     1              p(3*d_hhp(2,i)-2),p(3*d_hhp(2,i)-1),p(3*d_hhp(2,i)),
     1              p(3*d_hhp(3,i)-2),p(3*d_hhp(3,i)-1),p(3*d_hhp(3,i)),
     1              p(3*d_hhp(4,i)-2),p(3*d_hhp(4,i)-1),p(3*d_hhp(4,i)),
     1              angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
c        write(*,*)angle_abcd
        ehookdhhp=ehookdhhp+(angle_abcd-dhhp)**2
      enddo
      endif

      if(nd_hpp .ne. 0) then
C     1 hexagon, 2 pentagons
      do i=1,nd_hpp
      call dihedral(p(3*d_hpp(1,i)-2),p(3*d_hpp(1,i)-1),p(3*d_hpp(1,i)),
     1              p(3*d_hpp(2,i)-2),p(3*d_hpp(2,i)-1),p(3*d_hpp(2,i)),
     1              p(3*d_hpp(3,i)-2),p(3*d_hpp(3,i)-1),p(3*d_hpp(3,i)),
     1              p(3*d_hpp(4,i)-2),p(3*d_hpp(4,i)-1),p(3*d_hpp(4,i)),
     1              angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
c        write(*,*)angle_abcd
        ehookdhpp=ehookdhpp+(angle_abcd-dhpp)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif

      if(nd_ppp .ne. 0) then
C     3 pentagons
      do i=1,nd_ppp
c        write(*,*)a_p(1,i),a_p(2,i),a_p(3,i)
      call dihedral(p(3*d_ppp(1,i)-2),p(3*d_ppp(1,i)-1),p(3*d_ppp(1,i)),
     1              p(3*d_ppp(2,i)-2),p(3*d_ppp(2,i)-1),p(3*d_ppp(2,i)),
     1              p(3*d_ppp(3,i)-2),p(3*d_ppp(3,i)-1),p(3*d_ppp(3,i)),
     1              p(3*d_ppp(4,i)-2),p(3*d_ppp(4,i)-1),p(3*d_ppp(4,i)),
     1              angle_abcd)
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
c        write(*,*)angle_abcd
        ehookdppp=ehookdppp+(angle_abcd-dppp)**2
c        write(*,*)'diff',angle_p,ap
      enddo
      endif

C     Coulomb repulsion from origin
      ecoulomb=0.d0
      if (iopt.eq.4 .and. fco.ne.0.d0)  then
        Do I=1,n,3
          rinv=1.d0/dsqrt(p(I)**2+p(I+1)**2+p(I+2)**2)
          ecoulomb=ecoulomb+rinv
        enddo
      endif

C     total energy  
      fc=frpp*ehookrpp+frhp*ehookrhp+frhh*ehookrhh ! stretching
     2 +fap*ehookap+fah*ehookah ! bending
     3 +fdppp*ehookdppp+fdhpp*ehookdhpp+fdhhp*ehookdhhp+fdhhh*ehookdhhh! dihedral
     4 +fco*ecoulomb
      fc=0.5d0 * fc
c      write(*,*)fc,"energy"
      Return
      END SUBROUTINE



      SUBROUTINE dfunc3d(n,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
c      implicit real(8) (A-H,O-Z)
      integer iopt
c      write(*,*)"entering dfunc3d"
      
      select case(iopt)
        case(1, 2)
          CALL dwu(n,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p)
        case(3, 4)
          CALL dextwu(n,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
        case default
      end select

c      write(*,*)"leaving dfunc3d"
      return
      END SUBROUTINE


      SUBROUTINE dwu(n,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p)
c n=MATOM*3
      use config
      implicit real(8) (A-H,O-Z)
C     Wu force field in terms of harmonic oscillators for stretching
C     and bending, gradient
      real(8) p(nmax*3),x(nmax*3),force(ffmaxdim)
      integer n,iopt

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
     1               p(3*e_hh(2,i)-2),p(3*e_hh(2,i)-1),p(3*e_hh(2,i)),
     1               dax,day,daz,dbx,dby,dbz,ratom)
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
     1               p(3*e_hp(2,i)-2),p(3*e_hp(2,i)-1),p(3*e_hp(2,i)),
     1               dax,day,daz,dbx,dby,dbz,ratom)
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
     1               p(3*e_pp(2,i)-2),p(3*e_pp(2,i)-1),p(3*e_pp(2,i)),
     1               dax,day,daz,dbx,dby,dbz,ratom)
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


      SUBROUTINE dextwu(n,p,x,force,iopt,
     1 e_hh,e_hp,e_pp,ne_hh,ne_hp,ne_pp,
     1 a_h,a_p,
     1 d_hhh,d_hpp,d_hhp,d_ppp,nd_hhh,nd_hhp,nd_hpp,nd_ppp)
      use config
      implicit real(8) (A-H,O-Z)
      real(8) p(nmax*3),x(nmax*3),force(ffmaxdim)

c     edges with 0, 1, 2 pentagons
      integer e_hh(2,N/2), e_hp(2,N/2), e_pp(2,N/2)
      integer a_p(3,60), a_h(3,N-60)
      integer d_hhh(4,n/3),d_hpp(4,n/3),d_hhp(4,n/3),d_ppp(4,n/3)
c     counter for edges with 0, 1, 2 pentagons neighbours
      integer ne_hh,ne_hp,ne_pp
c     counter for dihedrals with 0, 1, 2, 3 pentagons neighbours
      integer nd_hhh,nd_hhp,nd_hpp,nd_ppp

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
      fco=force(19)

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
c        angle_abcd=dabs(angle_abcd)
        zero_value=sign(dhhh,angle_abcd)
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
c        angle_abcd=dabs(angle_abcd)
        zero_value=sign(dhhp,angle_abcd)
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
c        angle_abcd=dabs(angle_abcd)
        zero_value=sign(dhpp,angle_abcd)
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
c        angle_abcd=dabs(angle_abcd)
        zero_value=dsign(dppp, angle_abcd)
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

C     Coulomb repulsion from origin
      if (iopt.eq.4 .and. fco.ne.0.d0)  then
        Do I=1,n,3
          rinv=(p(I)**2+p(I+1)**2+p(I+2)**2)**(-1.5d0)
          x(I)  =x(I)  +fco*rinv*p(I)
          x(I+1)=x(I+1)+fco*rinv*p(I+1)
          x(I+2)=x(I+2)+fco*rinv*p(I+2)
        enddo
      endif

      return
      END
