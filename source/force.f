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


c subroutine dist takes 6 reals (=2 coordinates) and yields a positive distance
      SUBROUTINE DIST(ax,ay,az,bx,by,bz,dist_ab)
      implicit real*8 (a-z)
      dist_ab=dsqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)
      return
      END  


c subroutine ddist takes 6 reals (=2 coordinates) and yields all 6 first derivations of the distance
      SUBROUTINE DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,
     2  dist_ab)
      implicit real*8 (a-z)
      dist_ab=dsqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)
      dist_ab_inv=1/dist_ab
      aux_1=ax-bx
      aux_2=ay-by
      aux_3=az-bz
      dax=aux_1*dist_ab_inv
      dbx=-dax
      day=aux_2*dist_ab_inv
      dby=-day
      daz=aux_3*dist_ab_inv
      dbz=-daz
      return
      END  


c subroutine angle takes 9 reals (=3 coordinates) and yields an angel between 0 and +\pi (in radians)
c via law of cosines
c mult: 11, div: 1, root: 2, add/sub: 8, arccos: 1
      SUBROUTINE ANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      implicit real*8 (a-z)
      r2L=(ax-bx)**2 + (ay-by)**2 + (az-bz)**2
      r1L=dsqrt(r2L)
      r2M=(ax-cx)**2 + (ay-cy)**2 + (az-cz)**2
      r2R=(bx-cx)**2 + (by-cy)**2 + (bz-cz)**2
      r1R=dsqrt(r2R)
      angle_abc=dacos((r2L+r2R-r2M)/(2.0*r1L*r1R))
      return
      END


c subroutine angle takes 9 reals (=3 coordinates) and yields an angel between 0 and +\pi (in radians)
c via vector definition of the cosine
c mult: 10, div: 1, root: 2, add/sub: 12, arccos: 1
c      SUBROUTINE ANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
c      implicit real*8 (a-z)
c      aux_ax=ax-bx
c      aux_ay=ay-by
c      aux_az=az-bz
c      aux_bx=bx-cx
c      aux_by=by-cy
c      aux_bz=bz-cz
c      r2L=aux_ax**2 + aux_ay**2 + aux_az**2
c      r1L=dsqrt(r2L)
c      r2R=aux_bx**2 + aux_by**2 + aux_bz**2
c      r1R=dsqrt(r2R)
c      aux=dabs(aux_ax*aux_bx + aux_ay*aux_by + aux_az*aux_bz)
c      angle_abc=dacos(aux/(r1l*r1r))
c      return
c      END


c subroutine dangle takes 9 reals (=3 coordinates) and yields all 9 first derivations of the angle
c via law of cosines (calculating the derivative of Abs[foo] is rather troublesome)
      SUBROUTINE DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3 angle_abc)
      implicit real*8 (a-z)
c vectors from a to b and b to c and a to c
      aux_lx=ax-bx
      aux_ly=ay-by
      aux_lz=az-bz
      aux_rx=bx-cx
      aux_ry=by-cy
      aux_rz=bz-cz
      aux_mx=ax-cx
      aux_my=ay-cy
      aux_mz=az-cz
c length of a-b and b-c
      r2L=aux_lx**2 + aux_ly**2 + aux_lz**2
      r2R=aux_rx**2 + aux_ry**2 + aux_rz**2
      r2M=aux_mx**2 + aux_my**2 + aux_mz**2
      r1L=dsqrt(r2L)
      r1R=dsqrt(r2R)
      r1M=dsqrt(r2M)
      r3L=r2L*r1L
      r3R=r2R*r1R
      r3M=r2M*r1M
c some auxiliary products
      aux_11_inv=1/(2*r1L*r1R)
      aux_31_inv=1/(2*r3L*r1R)
      aux_13_inv=1/(2*r1L*r3R)
      aux_1=r2L + r2R - r2M
      arccos_arg=aux_1*aux_11_inv
c the actual angle, because it will always be required
      angle_abc=dacos(arccos_arg)
c not sure which is faster
      den_inv=-1/dsqrt(1-arccos_arg**2)
c      den_inv=-1/dabs(dsin(angle_abc))              
c more auxiliary products
      aux_2=2*aux_11_inv
      aux_3=aux_1*aux_31_inv
      aux_3x=aux_lx*aux_3
      aux_3y=aux_ly*aux_3
      aux_3z=aux_lz*aux_3
      aux_4=aux_1*aux_13_inv
      aux_4x=aux_rx*aux_4
      aux_4y=aux_ry*aux_4
      aux_4z=aux_rz*aux_4
c the derivations
      dax=((aux_lx-aux_mx)*aux_2-aux_3x)*den_inv
      day=((aux_ly-aux_my)*aux_2-aux_3y)*den_inv
      daz=((aux_lz-aux_mz)*aux_2-aux_3z)*den_inv
      dbx=((aux_rx-aux_lx)*aux_2-aux_4x+aux_3x)*den_inv
      dby=((aux_ry-aux_ly)*aux_2-aux_4y+aux_3y)*den_inv
      dbz=((aux_rz-aux_lz)*aux_2-aux_4z+aux_3z)*den_inv
      dcx=((aux_mx-aux_rx)*aux_2+aux_4x)*den_inv
      dcy=((aux_my-aux_ry)*aux_2+aux_4y)*den_inv
      dcz=((aux_mz-aux_rz)*aux_2+aux_4z)*den_inv
      return
      END


c subroutine dist takes 12 reals (=4 coordinates) and yields an angel between -\pi and +\pi (in radians)
      SUBROUTINE DIHEDRAL(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2 dihedral_abcd)
      IMPLICIT REAL*8 (a-z)
      ab_x=ax-bx
      ab_y=ay-by
      ab_z=az-bz
      ac_x=ax-cx
      ac_y=ay-cy
      ac_z=az-cz
      bc_x=bx-cx
      bc_y=by-cy
      bc_z=bz-cz
      bd_x=bx-dx
      bd_y=by-dy
      bd_z=bz-dz
      cd_x=cx-dx
      cd_y=cy-dy
      cd_z=cz-dz
c normal vectors on abc and bcd
c      abc_x=-az*by+ay*bz+az*cy-bz*cy-ay*cz+by*cz
c      abc_x=(cy-by)*az+(ay-cy)*bz+(by-ay)*cz
      abc_x=-bc_y*az + ac_y*bz - ab_y*cz
c      abc_y= az*bx-ax*bz-az*cx+bz*cx+ax*cz-bx*cz
c      abc_y=(bx-cx)*az+(cx-ax)*bz+(ax-bx)*cz
      abc_y=bc_x*az - ac_x*bz + ab_x*cz
c      abc_z=-ay*bx+ax*by+ay*cx-by*cx-ax*cy+bx*cy
c      abc_z=(cx-bx)*ay+(ax-cx)*by+(bx-ax)*cy
      abc_z=-bc_x*ay + ac_x*by - ab_x*cy
c      bcd_x=-bz*cy+by*cz+bz*dy-cz*dy-by*dz+cy*dz
c      bcd_x=(dy-cy)*bz+(by-dy)*cz+(cy-by)*dz
      bcd_x=-cd_y*bz + bd_y*cz - bc_y*dz
c      bcd_y= bz*cx-bx*cz-bz*dx+cz*dx+bx*dz-cx*dz
c      bcd_y=(cx-dx)*bz+(dx-bx)*cz+(bx-cx)*dz
      bcd_y=cd_x*bz - bd_x*cz + bc_x*dz
c      bcd_z=-by*cx+bx*cy+by*dx-cy*dx-bx*dy+cx*dy
c      bcd_z=(cy-dy)*bx+(dy-by)*cx+(by-cy)*dx
      bcd_z=cd_y*bx - bd_y*cx + bc_y*dx
c their respective lengths
      abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
c normal vectors (length 1) on abc and bcd
      abc_1_x=abc_x*abc_length_inv
      abc_1_y=abc_y*abc_length_inv
      abc_1_z=abc_z*abc_length_inv
      bcd_1_x=bcd_x*bcd_length_inv
      bcd_1_y=bcd_y*bcd_length_inv
      bcd_1_z=bcd_z*bcd_length_inv
c two auxiliary vectors
      bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      bc_1_x=bc_x*bc_length_inv
      bc_1_y=bc_y*bc_length_inv
      bc_1_z=bc_z*bc_length_inv
      aux_x=abc_1_y*bc_1_z-bc_1_y*abc_1_z
      aux_y=abc_1_z*bc_1_x-bc_1_z*abc_1_x
      aux_z=abc_1_x*bc_1_y-bc_1_x*abc_1_y
c two auxiliary reals
      aux_1=abc_1_x*bcd_1_x + abc_1_y*bcd_1_y + abc_1_z*bcd_1_z
      aux_2=aux_x*bcd_1_x + aux_y*bcd_1_y + aux_z*bcd_1_z
c the result
      dihedral_abcd=atan2(aux_2, aux_1)
      return
      END  


      SUBROUTINE DDIHEDRAL(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2  dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz)
      implicit real*8 (a-z)

      abx=ax-bx
      aby=ay-by
      abz=az-bz
      acx=ax-cx
      acy=ay-cy
      acz=az-cz
      bcx=bx-cx
      bcy=by-cy
      bcz=bz-cz
      bdx=bx-dx
      bdy=by-dy
      bdz=bz-dz
      cdx=cx-dx
      cdy=cy-dy
      cdz=cz-dz
c      write(*,*)"a-b",abx,aby,abz
c      write(*,*)"a-c",acx,acy,acz
c      write(*,*)"b-c",bcx,bcy,bcz
c      write(*,*)"b-d",bdx,bdy,bdz
c      write(*,*)"c-d",cdx,cdy,cdz

c normal vectors on abc and bcd
c      abc_x=-az*by+ay*bz+az*cy-bz*cy-ay*cz+by*cz
c      abc_x=(cy-by)*az+(ay-cy)*bz+(by-ay)*cz
      abcx=-bcy*az + acy*bz - aby*cz
c      abc_y= az*bx-ax*bz-az*cx+bz*cx+ax*cz-bx*cz
c      abc_y=(bx-cx)*az+(cx-ax)*bz+(ax-bx)*cz
      abcy= bcx*az - acx*bz + abx*cz
c      abc_z=-ay*bx+ax*by+ay*cx-by*cx-ax*cy+bx*cy
c      abc_z=(cx-bx)*ay+(ax-cx)*by+(bx-ax)*cy
      abcz=-bcx*ay + acx*by - abx*cy
c      bcd_x=-bz*cy+by*cz+bz*dy-cz*dy-by*dz+cy*dz
c      bcd_x=(dy-cy)*bz+(by-dy)*cz+(cy-by)*dz
      bcdx=-cdy*bz + bdy*cz - bcy*dz
c      bcd_y= bz*cx-bx*cz-bz*dx+cz*dx+bx*dz-cx*dz
c      bcd_y=(cx-dx)*bz+(dx-bx)*cz+(bx-cx)*dz
      bcdy= cdx*bz - bdx*cz + bcx*dz
c      bcd_z=-by*cx+bx*cy+by*dx-cy*dx-bx*dy+cx*dy
c      bcd_z=(cy-dy)*bx+(dy-by)*cx+(by-cy)*dx
      bcdz= cdy*bx - bdy*cx + bcy*dx
c      write(*,*)"abc",abcx,abcy,abcz
c      write(*,*)"bcd",bcdx,bcdy,bcdz

      BH= 2*abcz*acx - 2*abcx*acz
      BI=-2*bcdz*cdx + 2*bcdx*cdz
c      write(*,*)"b-hi",bh,bi

      abc2n=abcx**2 + abcy**2 + abcz**2
      bcd2n=bcdx**2 + bcdy**2 + bcdz**2
      bc2=bcx**2 + bcy**2 + bcz**2
c      write(*,*)"abc2n,bcd2n,bc2n",abc2n,bcd2n,bc2n

   
       dax= (abcx*bcdx*bcdy*bx*by + abcy*bcdy**2*bx*by + 
     -    abcz*bcdy*bcdz*bx*by - ax*bcdy*bcdz*bcy*bx*by + 
     -    ax*bcdy**2*bcz*bx*by - abcx*bcdx**2*by**2 - 
     -    abcy*bcdx*bcdy*by**2 - abcz*bcdx*bcdz*by**2 + 
     -    ax*bcdx*bcdz*bcy*by**2 - ax*bcdx*bcdy*bcz*by**2 + 
     -    abcx*bcdx*bcdz*bx*bz + abcy*bcdy*bcdz*bx*bz + 
     -    abcz*bcdz**2*bx*bz - ax*bcdz**2*bcy*bx*bz + 
     -    ax*bcdy*bcdz*bcz*bx*bz - abcx*bcdx**2*bz**2 - 
     -    abcy*bcdx*bcdy*bz**2 - abcz*bcdx*bcdz*bz**2 + 
     -    ax*bcdx*bcdz*bcy*bz**2 - ax*bcdx*bcdy*bcz*bz**2 - 
     -    abcx*bcdx*bcdy*by*cx - abcy*bcdy**2*by*cx - 
     -    abcz*bcdy*bcdz*by*cx + ax*bcdy*bcdz*bcy*by*cx - 
     -    ax*bcdy**2*bcz*by*cx + bcdy*bcdz*bcy*bx*by*cx - 
     -    bcdy**2*bcz*bx*by*cx - bcdx*bcdz*bcy*by**2*cx + 
     -    bcdx*bcdy*bcz*by**2*cx - abcx*bcdx*bcdz*bz*cx - 
     -    abcy*bcdy*bcdz*bz*cx - abcz*bcdz**2*bz*cx + 
     -    ax*bcdz**2*bcy*bz*cx - ax*bcdy*bcdz*bcz*bz*cx + 
     -    bcdz**2*bcy*bx*bz*cx - bcdy*bcdz*bcz*bx*bz*cx - 
     -    bcdx*bcdz*bcy*bz**2*cx + bcdx*bcdy*bcz*bz**2*cx - 
     -    bcdy*bcdz*bcy*by*cx**2 + bcdy**2*bcz*by*cx**2 - 
     -    bcdz**2*bcy*bz*cx**2 + bcdy*bcdz*bcz*bz*cx**2 - 
     -    abcx*bcdx*bcdy*bx*cy - abcy*bcdy**2*bx*cy - 
     -    abcz*bcdy*bcdz*bx*cy + ax*bcdy*bcdz*bcy*bx*cy - 
     -    ax*bcdy**2*bcz*bx*cy - bcdy*bcdz*bcy*bx**2*cy + 
     -    bcdy**2*bcz*bx**2*cy + 2*abcx*bcdx**2*by*cy + 
     -    2*abcy*bcdx*bcdy*by*cy + 2*abcz*bcdx*bcdz*by*cy - 
     -    2*ax*bcdx*bcdz*bcy*by*cy + 2*ax*bcdx*bcdy*bcz*by*cy + 
     -    bcdx*bcdz*bcy*bx*by*cy - bcdx*bcdy*bcz*bx*by*cy + 
     -    bcdz**2*bcy*by*bz*cy - bcdy*bcdz*bcz*by*bz*cy - 
     -    bcdy*bcdz*bcy*bz**2*cy + bcdy**2*bcz*bz**2*cy + 
     -    abcx*bcdx*bcdy*cx*cy + abcy*bcdy**2*cx*cy + 
     -    abcz*bcdy*bcdz*cx*cy - ax*bcdy*bcdz*bcy*cx*cy + 
     -    ax*bcdy**2*bcz*cx*cy + bcdy*bcdz*bcy*bx*cx*cy - 
     -    bcdy**2*bcz*bx*cx*cy + bcdx*bcdz*bcy*by*cx*cy - 
     -    bcdx*bcdy*bcz*by*cx*cy - abcx*bcdx**2*cy**2 - 
     -    abcy*bcdx*bcdy*cy**2 - abcz*bcdx*bcdz*cy**2 + 
     -    ax*bcdx*bcdz*bcy*cy**2 - ax*bcdx*bcdy*bcz*cy**2 - 
     -    bcdx*bcdz*bcy*bx*cy**2 + bcdx*bcdy*bcz*bx*cy**2 - 
     -    bcdz**2*bcy*bz*cy**2 + bcdy*bcdz*bcz*bz*cy**2 + 
     -    az*(bcdz*bcy - bcdy*bcz)*
     -     (bcdz*(bcx**2 + bcy**2) + 
     -       bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) - 
     -    abcx*bcdx*bcdz*bx*cz - abcy*bcdy*bcdz*bx*cz - 
     -    abcz*bcdz**2*bx*cz + ax*bcdz**2*bcy*bx*cz - 
     -    ax*bcdy*bcdz*bcz*bx*cz - bcdz**2*bcy*bx**2*cz + 
     -    bcdy*bcdz*bcz*bx**2*cz - bcdz**2*bcy*by**2*cz + 
     -    bcdy*bcdz*bcz*by**2*cz + 2*abcx*bcdx**2*bz*cz + 
     -    2*abcy*bcdx*bcdy*bz*cz + 2*abcz*bcdx*bcdz*bz*cz - 
     -    2*ax*bcdx*bcdz*bcy*bz*cz + 2*ax*bcdx*bcdy*bcz*bz*cz + 
     -    bcdx*bcdz*bcy*bx*bz*cz - bcdx*bcdy*bcz*bx*bz*cz + 
     -    bcdy*bcdz*bcy*by*bz*cz - bcdy**2*bcz*by*bz*cz + 
     -    abcx*bcdx*bcdz*cx*cz + abcy*bcdy*bcdz*cx*cz + 
     -    abcz*bcdz**2*cx*cz - ax*bcdz**2*bcy*cx*cz + 
     -    ax*bcdy*bcdz*bcz*cx*cz + bcdz**2*bcy*bx*cx*cz - 
     -    bcdy*bcdz*bcz*bx*cx*cz + bcdx*bcdz*bcy*bz*cx*cz - 
     -    bcdx*bcdy*bcz*bz*cx*cz + bcdz**2*bcy*by*cy*cz - 
     -    bcdy*bcdz*bcz*by*cy*cz + bcdy*bcdz*bcy*bz*cy*cz - 
     -    bcdy**2*bcz*bz*cy*cz - abcx*bcdx**2*cz**2 - 
     -    abcy*bcdx*bcdy*cz**2 - abcz*bcdx*bcdz*cz**2 + 
     -    ax*bcdx*bcdz*bcy*cz**2 - ax*bcdx*bcdy*bcz*cz**2 - 
     -    bcdx*bcdz*bcy*bx*cz**2 + bcdx*bcdy*bcz*bx*cz**2 - 
     -    bcdy*bcdz*bcy*by*cz**2 + bcdy**2*bcz*by*cz**2 - 
     -    ay*(-(bcdz*bcy) + bcdy*bcz)*
     -     (bcdy*(bcx**2 + bcz**2) + 
     -       bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))/
     -  (Sqrt(bc2)*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy))
     -           + ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))
     -         **2/bc2))


       day=(-(abcz*bcdy*bcdz*bx**2) - ay*bcdy*bcdz*bcx*bx**2 - 
     -    az*bcdz**2*bcx*bx**2 + ay*bcdx*bcdy*bcz*bx**2 + 
     -    az*bcdx*bcdz*bcz*bx**2 + abcz*bcdx*bcdz*bx*by + 
     -    ay*bcdx*bcdz*bcx*bx*by + ax*bcdy*bcdz*bcx*bx*by - 
     -    ay*bcdx**2*bcz*bx*by - ax*bcdx*bcdy*bcz*bx*by - 
     -    ax*bcdx*bcdz*bcx*by**2 - az*bcdz**2*bcx*by**2 + 
     -    ax*bcdx**2*bcz*by**2 + az*bcdx*bcdz*bcz*by**2 + 
     -    az*bcdx*bcdz*bcx*bx*bz + ax*bcdz**2*bcx*bx*bz - 
     -    az*bcdx**2*bcz*bx*bz - ax*bcdx*bcdz*bcz*bx*bz + 
     -    abcz*bcdz**2*by*bz + az*bcdy*bcdz*bcx*by*bz + 
     -    ay*bcdz**2*bcx*by*bz - az*bcdx*bcdy*bcz*by*bz - 
     -    ay*bcdx*bcdz*bcz*by*bz - abcz*bcdy*bcdz*bz**2 - 
     -    ax*bcdx*bcdz*bcx*bz**2 - ay*bcdy*bcdz*bcx*bz**2 + 
     -    ax*bcdx**2*bcz*bz**2 + ay*bcdx*bcdy*bcz*bz**2 + 
     -    2*abcz*bcdy*bcdz*bx*cx + 2*ay*bcdy*bcdz*bcx*bx*cx + 
     -    2*az*bcdz**2*bcx*bx*cx - 2*ay*bcdx*bcdy*bcz*bx*cx - 
     -    2*az*bcdx*bcdz*bcz*bx*cx - abcz*bcdx*bcdz*by*cx - 
     -    ay*bcdx*bcdz*bcx*by*cx - ax*bcdy*bcdz*bcx*by*cx + 
     -    ay*bcdx**2*bcz*by*cx + ax*bcdx*bcdy*bcz*by*cx - 
     -    bcdy*bcdz*bcx*bx*by*cx + bcdx*bcdy*bcz*bx*by*cx + 
     -    bcdx*bcdz*bcx*by**2*cx - bcdx**2*bcz*by**2*cx - 
     -    az*bcdx*bcdz*bcx*bz*cx - ax*bcdz**2*bcx*bz*cx + 
     -    az*bcdx**2*bcz*bz*cx + ax*bcdx*bcdz*bcz*bz*cx - 
     -    bcdz**2*bcx*bx*bz*cx + bcdx*bcdz*bcz*bx*bz*cx + 
     -    bcdx*bcdz*bcx*bz**2*cx - bcdx**2*bcz*bz**2*cx - 
     -    abcz*bcdy*bcdz*cx**2 - ay*bcdy*bcdz*bcx*cx**2 - 
     -    az*bcdz**2*bcx*cx**2 + ay*bcdx*bcdy*bcz*cx**2 + 
     -    az*bcdx*bcdz*bcz*cx**2 + bcdy*bcdz*bcx*by*cx**2 - 
     -    bcdx*bcdy*bcz*by*cx**2 + bcdz**2*bcx*bz*cx**2 - 
     -    bcdx*bcdz*bcz*bz*cx**2 - abcz*bcdx*bcdz*bx*cy - 
     -    ay*bcdx*bcdz*bcx*bx*cy - ax*bcdy*bcdz*bcx*bx*cy + 
     -    ay*bcdx**2*bcz*bx*cy + ax*bcdx*bcdy*bcz*bx*cy + 
     -    bcdy*bcdz*bcx*bx**2*cy - bcdx*bcdy*bcz*bx**2*cy + 
     -    2*ax*bcdx*bcdz*bcx*by*cy + 2*az*bcdz**2*bcx*by*cy - 
     -    2*ax*bcdx**2*bcz*by*cy - 2*az*bcdx*bcdz*bcz*by*cy - 
     -    bcdx*bcdz*bcx*bx*by*cy + bcdx**2*bcz*bx*by*cy - 
     -    abcz*bcdz**2*bz*cy - az*bcdy*bcdz*bcx*bz*cy - 
     -    ay*bcdz**2*bcx*bz*cy + az*bcdx*bcdy*bcz*bz*cy + 
     -    ay*bcdx*bcdz*bcz*bz*cy - bcdz**2*bcx*by*bz*cy + 
     -    bcdx*bcdz*bcz*by*bz*cy + bcdy*bcdz*bcx*bz**2*cy - 
     -    bcdx*bcdy*bcz*bz**2*cy + abcz*bcdx*bcdz*cx*cy + 
     -    ay*bcdx*bcdz*bcx*cx*cy + ax*bcdy*bcdz*bcx*cx*cy - 
     -    ay*bcdx**2*bcz*cx*cy - ax*bcdx*bcdy*bcz*cx*cy - 
     -    bcdy*bcdz*bcx*bx*cx*cy + bcdx*bcdy*bcz*bx*cx*cy - 
     -    bcdx*bcdz*bcx*by*cx*cy + bcdx**2*bcz*by*cx*cy - 
     -    ax*bcdx*bcdz*bcx*cy**2 - az*bcdz**2*bcx*cy**2 + 
     -    ax*bcdx**2*bcz*cy**2 + az*bcdx*bcdz*bcz*cy**2 + 
     -    bcdx*bcdz*bcx*bx*cy**2 - bcdx**2*bcz*bx*cy**2 + 
     -    bcdz**2*bcx*bz*cy**2 - bcdx*bcdz*bcz*bz*cy**2 - 
     -    az*bcdx*bcdz*bcx*bx*cz - ax*bcdz**2*bcx*bx*cz + 
     -    az*bcdx**2*bcz*bx*cz + ax*bcdx*bcdz*bcz*bx*cz + 
     -    bcdz**2*bcx*bx**2*cz - bcdx*bcdz*bcz*bx**2*cz - 
     -    abcz*bcdz**2*by*cz - az*bcdy*bcdz*bcx*by*cz - 
     -    ay*bcdz**2*bcx*by*cz + az*bcdx*bcdy*bcz*by*cz + 
     -    ay*bcdx*bcdz*bcz*by*cz + bcdz**2*bcx*by**2*cz - 
     -    bcdx*bcdz*bcz*by**2*cz + 2*abcz*bcdy*bcdz*bz*cz + 
     -    2*ax*bcdx*bcdz*bcx*bz*cz + 2*ay*bcdy*bcdz*bcx*bz*cz - 
     -    2*ax*bcdx**2*bcz*bz*cz - 2*ay*bcdx*bcdy*bcz*bz*cz - 
     -    bcdx*bcdz*bcx*bx*bz*cz + bcdx**2*bcz*bx*bz*cz - 
     -    bcdy*bcdz*bcx*by*bz*cz + bcdx*bcdy*bcz*by*bz*cz + 
     -    az*bcdx*bcdz*bcx*cx*cz + ax*bcdz**2*bcx*cx*cz - 
     -    az*bcdx**2*bcz*cx*cz - ax*bcdx*bcdz*bcz*cx*cz - 
     -    bcdz**2*bcx*bx*cx*cz + bcdx*bcdz*bcz*bx*cx*cz - 
     -    bcdx*bcdz*bcx*bz*cx*cz + bcdx**2*bcz*bz*cx*cz + 
     -    abcz*bcdz**2*cy*cz + az*bcdy*bcdz*bcx*cy*cz + 
     -    ay*bcdz**2*bcx*cy*cz - az*bcdx*bcdy*bcz*cy*cz - 
     -    ay*bcdx*bcdz*bcz*cy*cz - bcdz**2*bcx*by*cy*cz + 
     -    bcdx*bcdz*bcz*by*cy*cz - bcdy*bcdz*bcx*bz*cy*cz + 
     -    bcdx*bcdy*bcz*bz*cy*cz - abcz*bcdy*bcdz*cz**2 - 
     -    ax*bcdx*bcdz*bcx*cz**2 - ay*bcdy*bcdz*bcx*cz**2 + 
     -    ax*bcdx**2*bcz*cz**2 + ay*bcdx*bcdy*bcz*cz**2 + 
     -    bcdx*bcdz*bcx*bx*cz**2 - bcdx**2*bcz*bx*cz**2 + 
     -    bcdy*bcdz*bcx*by*cz**2 - bcdx*bcdy*bcz*by*cz**2 + 
     -    abcx*bcdx*(-(bcdy*(bcx**2 + bcz**2)) + 
     -       bcy*(bcdx*bx + bcdz*bz - bcdx*cx - bcdz*cz)) - 
     -    abcy*bcdy*(bcdy*(bcx**2 + bcz**2) + 
     -       bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))/
     -  (Sqrt(bc2)*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy))
     -           + ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))
     -         **2/bc2))


       daz= (-(abcz*bcdz**2*bx**2) + ay*bcdy**2*bcx*bx**2 + 
     -    az*bcdy*bcdz*bcx*bx**2 - ay*bcdx*bcdy*bcy*bx**2 - 
     -    az*bcdx*bcdz*bcy*bx**2 - ay*bcdx*bcdy*bcx*bx*by - 
     -    ax*bcdy**2*bcx*bx*by + ay*bcdx**2*bcy*bx*by + 
     -    ax*bcdx*bcdy*bcy*bx*by - abcz*bcdz**2*by**2 + 
     -    ax*bcdx*bcdy*bcx*by**2 + az*bcdy*bcdz*bcx*by**2 - 
     -    ax*bcdx**2*bcy*by**2 - az*bcdx*bcdz*bcy*by**2 + 
     -    abcz*bcdx*bcdz*bx*bz - az*bcdx*bcdy*bcx*bx*bz - 
     -    ax*bcdy*bcdz*bcx*bx*bz + az*bcdx**2*bcy*bx*bz + 
     -    ax*bcdx*bcdz*bcy*bx*bz + abcz*bcdy*bcdz*by*bz - 
     -    az*bcdy**2*bcx*by*bz - ay*bcdy*bcdz*bcx*by*bz + 
     -    az*bcdx*bcdy*bcy*by*bz + ay*bcdx*bcdz*bcy*by*bz + 
     -    ax*bcdx*bcdy*bcx*bz**2 + ay*bcdy**2*bcx*bz**2 - 
     -    ax*bcdx**2*bcy*bz**2 - ay*bcdx*bcdy*bcy*bz**2 + 
     -    2*abcz*bcdz**2*bx*cx - 2*ay*bcdy**2*bcx*bx*cx - 
     -    2*az*bcdy*bcdz*bcx*bx*cx + 2*ay*bcdx*bcdy*bcy*bx*cx + 
     -    2*az*bcdx*bcdz*bcy*bx*cx + ay*bcdx*bcdy*bcx*by*cx + 
     -    ax*bcdy**2*bcx*by*cx - ay*bcdx**2*bcy*by*cx - 
     -    ax*bcdx*bcdy*bcy*by*cx + bcdy**2*bcx*bx*by*cx - 
     -    bcdx*bcdy*bcy*bx*by*cx - bcdx*bcdy*bcx*by**2*cx + 
     -    bcdx**2*bcy*by**2*cx - abcz*bcdx*bcdz*bz*cx + 
     -    az*bcdx*bcdy*bcx*bz*cx + ax*bcdy*bcdz*bcx*bz*cx - 
     -    az*bcdx**2*bcy*bz*cx - ax*bcdx*bcdz*bcy*bz*cx + 
     -    bcdy*bcdz*bcx*bx*bz*cx - bcdx*bcdz*bcy*bx*bz*cx - 
     -    bcdx*bcdy*bcx*bz**2*cx + bcdx**2*bcy*bz**2*cx - 
     -    abcz*bcdz**2*cx**2 + ay*bcdy**2*bcx*cx**2 + 
     -    az*bcdy*bcdz*bcx*cx**2 - ay*bcdx*bcdy*bcy*cx**2 - 
     -    az*bcdx*bcdz*bcy*cx**2 - bcdy**2*bcx*by*cx**2 + 
     -    bcdx*bcdy*bcy*by*cx**2 - bcdy*bcdz*bcx*bz*cx**2 + 
     -    bcdx*bcdz*bcy*bz*cx**2 + ay*bcdx*bcdy*bcx*bx*cy + 
     -    ax*bcdy**2*bcx*bx*cy - ay*bcdx**2*bcy*bx*cy - 
     -    ax*bcdx*bcdy*bcy*bx*cy - bcdy**2*bcx*bx**2*cy + 
     -    bcdx*bcdy*bcy*bx**2*cy + 2*abcz*bcdz**2*by*cy - 
     -    2*ax*bcdx*bcdy*bcx*by*cy - 2*az*bcdy*bcdz*bcx*by*cy + 
     -    2*ax*bcdx**2*bcy*by*cy + 2*az*bcdx*bcdz*bcy*by*cy + 
     -    bcdx*bcdy*bcx*bx*by*cy - bcdx**2*bcy*bx*by*cy - 
     -    abcz*bcdy*bcdz*bz*cy + az*bcdy**2*bcx*bz*cy + 
     -    ay*bcdy*bcdz*bcx*bz*cy - az*bcdx*bcdy*bcy*bz*cy - 
     -    ay*bcdx*bcdz*bcy*bz*cy + bcdy*bcdz*bcx*by*bz*cy - 
     -    bcdx*bcdz*bcy*by*bz*cy - bcdy**2*bcx*bz**2*cy + 
     -    bcdx*bcdy*bcy*bz**2*cy - ay*bcdx*bcdy*bcx*cx*cy - 
     -    ax*bcdy**2*bcx*cx*cy + ay*bcdx**2*bcy*cx*cy + 
     -    ax*bcdx*bcdy*bcy*cx*cy + bcdy**2*bcx*bx*cx*cy - 
     -    bcdx*bcdy*bcy*bx*cx*cy + bcdx*bcdy*bcx*by*cx*cy - 
     -    bcdx**2*bcy*by*cx*cy - abcz*bcdz**2*cy**2 + 
     -    ax*bcdx*bcdy*bcx*cy**2 + az*bcdy*bcdz*bcx*cy**2 - 
     -    ax*bcdx**2*bcy*cy**2 - az*bcdx*bcdz*bcy*cy**2 - 
     -    bcdx*bcdy*bcx*bx*cy**2 + bcdx**2*bcy*bx*cy**2 - 
     -    bcdy*bcdz*bcx*bz*cy**2 + bcdx*bcdz*bcy*bz*cy**2 + 
     -    abcx*bcdx*(-(bcdz*(bcx**2 + bcy**2)) + 
     -       bcz*(bcdx*bx + bcdy*by - bcdx*cx - bcdy*cy)) + 
     -    abcy*bcdy*(-(bcdz*(bcx**2 + bcy**2)) + 
     -       bcz*(bcdx*bx + bcdy*by - bcdx*cx - bcdy*cy)) - 
     -    abcz*bcdx*bcdz*bx*cz + az*bcdx*bcdy*bcx*bx*cz + 
     -    ax*bcdy*bcdz*bcx*bx*cz - az*bcdx**2*bcy*bx*cz - 
     -    ax*bcdx*bcdz*bcy*bx*cz - bcdy*bcdz*bcx*bx**2*cz + 
     -    bcdx*bcdz*bcy*bx**2*cz - abcz*bcdy*bcdz*by*cz + 
     -    az*bcdy**2*bcx*by*cz + ay*bcdy*bcdz*bcx*by*cz - 
     -    az*bcdx*bcdy*bcy*by*cz - ay*bcdx*bcdz*bcy*by*cz - 
     -    bcdy*bcdz*bcx*by**2*cz + bcdx*bcdz*bcy*by**2*cz - 
     -    2*ax*bcdx*bcdy*bcx*bz*cz - 2*ay*bcdy**2*bcx*bz*cz + 
     -    2*ax*bcdx**2*bcy*bz*cz + 2*ay*bcdx*bcdy*bcy*bz*cz + 
     -    bcdx*bcdy*bcx*bx*bz*cz - bcdx**2*bcy*bx*bz*cz + 
     -    bcdy**2*bcx*by*bz*cz - bcdx*bcdy*bcy*by*bz*cz + 
     -    abcz*bcdx*bcdz*cx*cz - az*bcdx*bcdy*bcx*cx*cz - 
     -    ax*bcdy*bcdz*bcx*cx*cz + az*bcdx**2*bcy*cx*cz + 
     -    ax*bcdx*bcdz*bcy*cx*cz + bcdy*bcdz*bcx*bx*cx*cz - 
     -    bcdx*bcdz*bcy*bx*cx*cz + bcdx*bcdy*bcx*bz*cx*cz - 
     -    bcdx**2*bcy*bz*cx*cz + abcz*bcdy*bcdz*cy*cz - 
     -    az*bcdy**2*bcx*cy*cz - ay*bcdy*bcdz*bcx*cy*cz + 
     -    az*bcdx*bcdy*bcy*cy*cz + ay*bcdx*bcdz*bcy*cy*cz + 
     -    bcdy*bcdz*bcx*by*cy*cz - bcdx*bcdz*bcy*by*cy*cz + 
     -    bcdy**2*bcx*bz*cy*cz - bcdx*bcdy*bcy*bz*cy*cz + 
     -    ax*bcdx*bcdy*bcx*cz**2 + ay*bcdy**2*bcx*cz**2 - 
     -    ax*bcdx**2*bcy*cz**2 - ay*bcdx*bcdy*bcy*cz**2 - 
     -    bcdx*bcdy*bcx*bx*cz**2 + bcdx**2*bcy*bx*cz**2 - 
     -    bcdy**2*bcx*by*cz**2 + bcdx*bcdy*bcy*by*cz**2)/
     -  (Sqrt(bc2)*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy))
     -           + ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))
     -         **2/bc2))


       dbx=(-2*abcz*ay*bc2*bcdy*bcdz*bx - 2*abcz*az*bc2*bcdz**2*bx + 
     -    ay*az*bc2*bcdy**2*bx**2 - ay**2*bc2*bcdy*bcdz*bx**2 + 
     -    az**2*bc2*bcdy*bcdz*bx**2 - ay*az*bc2*bcdz**2*bx**2 + 
     -    abcz*ay*bcdy*bcdz*bcx*bx**2 + abcz*az*bcdz**2*bcx*bx**2 + 
     -    abcz*ay*bc2*bcdx*bcdz*by + abcz*ax*bc2*bcdy*bcdz*by - 
     -    ay*az*bc2*bcdx*bcdy*bx*by - ax*az*bc2*bcdy**2*bx*by + 
     -    ay**2*bc2*bcdx*bcdz*bx*by + ax*ay*bc2*bcdy*bcdz*bx*by - 
     -    abcz*ay*bcdx*bcdz*bcx*bx*by - abcz*ax*bcdy*bcdz*bcx*bx*by + 
     -    ax*az*bc2*bcdx*bcdy*by**2 - ax*ay*bc2*bcdx*bcdz*by**2 + 
     -    az**2*bc2*bcdy*bcdz*by**2 - ay*az*bc2*bcdz**2*by**2 + 
     -    abcz*ax*bcdx*bcdz*bcx*by**2 + abcz*az*bcdz**2*bcx*by**2 + 
     -    abcz*az*bc2*bcdx*bcdz*bz + abcz*ax*bc2*bcdz**2*bz - 
     -    az**2*bc2*bcdx*bcdy*bx*bz + ay*az*bc2*bcdx*bcdz*bx*bz - 
     -    ax*az*bc2*bcdy*bcdz*bx*bz + ax*ay*bc2*bcdz**2*bx*bz - 
     -    abcz*az*bcdx*bcdz*bcx*bx*bz - abcz*ax*bcdz**2*bcx*bx*bz - 
     -    az**2*bc2*bcdy**2*by*bz + ay**2*bc2*bcdz**2*by*bz - 
     -    abcz*az*bcdy*bcdz*bcx*by*bz - abcz*ay*bcdz**2*bcx*by*bz + 
     -    ax*az*bc2*bcdx*bcdy*bz**2 + ay*az*bc2*bcdy**2*bz**2 - 
     -    ax*ay*bc2*bcdx*bcdz*bz**2 - ay**2*bc2*bcdy*bcdz*bz**2 + 
     -    abcz*ax*bcdx*bcdz*bcx*bz**2 + abcz*ay*bcdy*bcdz*bcx*bz**2 + 
     -    abcz*ay*bc2*bcdy*bx**2*cdy - abcz*ay*bc2*bcdx*bx*by*cdy - 
     -    abcz*ax*bc2*bcdy*bx*by*cdy + abcz*ax*bc2*bcdx*by**2*cdy - 
     -    abcz*az*bc2*bcdx*bx*bz*cdy - abcz*az*bc2*bcdy*by*bz*cdy + 
     -    abcz*ax*bc2*bcdx*bz**2*cdy + abcz*ay*bc2*bcdy*bz**2*cdy + 
     -    abcz*ay*bc2*bcdz*bx**2*cdz - abcz*ax*bc2*bcdz*bx*by*cdz - 
     -    abcz*az*bc2*bcdz*by*bz*cdz + abcz*ay*bc2*bcdz*bz**2*cdz + 
     -    2*abcz*ay*bc2*bcdy*bcdz*cx + 2*abcz*az*bc2*bcdz**2*cx - 
     -    2*ay*az*bc2*bcdy**2*bx*cx + 2*ay**2*bc2*bcdy*bcdz*bx*cx - 
     -    2*az**2*bc2*bcdy*bcdz*bx*cx + 2*ay*az*bc2*bcdz**2*bx*cx - 
     -    2*abcz*ay*bcdy*bcdz*bcx*bx*cx - 2*abcz*az*bcdz**2*bcx*bx*cx + 
     -    ay*az*bc2*bcdx*bcdy*by*cx + ax*az*bc2*bcdy**2*by*cx - 
     -    ay**2*bc2*bcdx*bcdz*by*cx - abcz*bc2*bcdy*bcdz*by*cx - 
     -    ax*ay*bc2*bcdy*bcdz*by*cx + abcz*ay*bcdx*bcdz*bcx*by*cx + 
     -    abcz*ax*bcdy*bcdz*bcx*by*cx + az*bc2*bcdy**2*bx*by*cx - 
     -    ay*bc2*bcdy*bcdz*bx*by*cx + abcz*bcdy*bcdz*bcx*bx*by*cx - 
     -    az*bc2*bcdx*bcdy*by**2*cx + ay*bc2*bcdx*bcdz*by**2*cx - 
     -    abcz*bcdx*bcdz*bcx*by**2*cx + az**2*bc2*bcdx*bcdy*bz*cx - 
     -    ay*az*bc2*bcdx*bcdz*bz*cx + ax*az*bc2*bcdy*bcdz*bz*cx - 
     -    abcz*bc2*bcdz**2*bz*cx - ax*ay*bc2*bcdz**2*bz*cx + 
     -    abcz*az*bcdx*bcdz*bcx*bz*cx + abcz*ax*bcdz**2*bcx*bz*cx + 
     -    az*bc2*bcdy*bcdz*bx*bz*cx - ay*bc2*bcdz**2*bx*bz*cx + 
     -    abcz*bcdz**2*bcx*bx*bz*cx - az*bc2*bcdx*bcdy*bz**2*cx + 
     -    ay*bc2*bcdx*bcdz*bz**2*cx - abcz*bcdx*bcdz*bcx*bz**2*cx - 
     -    2*abcz*ay*bc2*bcdy*bx*cdy*cx + abcz*ay*bc2*bcdx*by*cdy*cx + 
     -    abcz*ax*bc2*bcdy*by*cdy*cx + abcz*bc2*bcdy*bx*by*cdy*cx - 
     -    abcz*bc2*bcdx*by**2*cdy*cx + abcz*az*bc2*bcdx*bz*cdy*cx - 
     -    abcz*bc2*bcdx*bz**2*cdy*cx - 2*abcz*ay*bc2*bcdz*bx*cdz*cx + 
     -    abcz*ax*bc2*bcdz*by*cdz*cx + abcz*bc2*bcdz*bx*by*cdz*cx + 
     -    ay*az*bc2*bcdy**2*cx**2 - ay**2*bc2*bcdy*bcdz*cx**2 + 
     -    az**2*bc2*bcdy*bcdz*cx**2 - ay*az*bc2*bcdz**2*cx**2 + 
     -    abcz*ay*bcdy*bcdz*bcx*cx**2 + abcz*az*bcdz**2*bcx*cx**2 - 
     -    az*bc2*bcdy**2*by*cx**2 + ay*bc2*bcdy*bcdz*by*cx**2 - 
     -    abcz*bcdy*bcdz*bcx*by*cx**2 - az*bc2*bcdy*bcdz*bz*cx**2 + 
     -    ay*bc2*bcdz**2*bz*cx**2 - abcz*bcdz**2*bcx*bz*cx**2 + 
     -    abcz*ay*bc2*bcdy*cdy*cx**2 - abcz*bc2*bcdy*by*cdy*cx**2 + 
     -    abcz*ay*bc2*bcdz*cdz*cx**2 - abcz*bc2*bcdz*by*cdz*cx**2 - 
     -    abcz*ay*bc2*bcdx*bcdz*cy - abcz*ax*bc2*bcdy*bcdz*cy + 
     -    ay*az*bc2*bcdx*bcdy*bx*cy + ax*az*bc2*bcdy**2*bx*cy - 
     -    ay**2*bc2*bcdx*bcdz*bx*cy + 2*abcz*bc2*bcdy*bcdz*bx*cy - 
     -    ax*ay*bc2*bcdy*bcdz*bx*cy + abcz*ay*bcdx*bcdz*bcx*bx*cy + 
     -    abcz*ax*bcdy*bcdz*bcx*bx*cy - az*bc2*bcdy**2*bx**2*cy + 
     -    2*ay*bc2*bcdy*bcdz*bx**2*cy + az*bc2*bcdz**2*bx**2*cy - 
     -    abcz*bcdy*bcdz*bcx*bx**2*cy - 2*ax*az*bc2*bcdx*bcdy*by*cy - 
     -    abcz*bc2*bcdx*bcdz*by*cy + 2*ax*ay*bc2*bcdx*bcdz*by*cy - 
     -    2*az**2*bc2*bcdy*bcdz*by*cy + 2*ay*az*bc2*bcdz**2*by*cy - 
     -    2*abcz*ax*bcdx*bcdz*bcx*by*cy - 2*abcz*az*bcdz**2*bcx*by*cy + 
     -    az*bc2*bcdx*bcdy*bx*by*cy - 2*ay*bc2*bcdx*bcdz*bx*by*cy - 
     -    ax*bc2*bcdy*bcdz*bx*by*cy + abcz*bcdx*bcdz*bcx*bx*by*cy + 
     -    ax*bc2*bcdx*bcdz*by**2*cy + az*bc2*bcdz**2*by**2*cy + 
     -    az**2*bc2*bcdy**2*bz*cy - ay**2*bc2*bcdz**2*bz*cy + 
     -    abcz*az*bcdy*bcdz*bcx*bz*cy + abcz*ay*bcdz**2*bcx*bz*cy - 
     -    az*bc2*bcdx*bcdz*bx*bz*cy - ax*bc2*bcdz**2*bx*bz*cy - 
     -    2*ay*bc2*bcdz**2*by*bz*cy + abcz*bcdz**2*bcx*by*bz*cy - 
     -    az*bc2*bcdy**2*bz**2*cy + ax*bc2*bcdx*bcdz*bz**2*cy + 
     -    2*ay*bc2*bcdy*bcdz*bz**2*cy - abcz*bcdy*bcdz*bcx*bz**2*cy + 
     -    abcz*ay*bc2*bcdx*bx*cdy*cy + abcz*ax*bc2*bcdy*bx*cdy*cy - 
     -    abcz*bc2*bcdy*bx**2*cdy*cy - 2*abcz*ax*bc2*bcdx*by*cdy*cy + 
     -    abcz*bc2*bcdx*bx*by*cdy*cy + abcz*az*bc2*bcdy*bz*cdy*cy - 
     -    abcz*bc2*bcdy*bz**2*cdy*cy + abcz*ax*bc2*bcdz*bx*cdz*cy - 
     -    abcz*bc2*bcdz*bx**2*cdz*cy + abcz*az*bc2*bcdz*bz*cdz*cy - 
     -    abcz*bc2*bcdz*bz**2*cdz*cy - ay*az*bc2*bcdx*bcdy*cx*cy - 
     -    ax*az*bc2*bcdy**2*cx*cy + ay**2*bc2*bcdx*bcdz*cx*cy - 
     -    abcz*bc2*bcdy*bcdz*cx*cy + ax*ay*bc2*bcdy*bcdz*cx*cy - 
     -    abcz*ay*bcdx*bcdz*bcx*cx*cy - abcz*ax*bcdy*bcdz*bcx*cx*cy + 
     -    az*bc2*bcdy**2*bx*cx*cy - 3*ay*bc2*bcdy*bcdz*bx*cx*cy - 
     -    2*az*bc2*bcdz**2*bx*cx*cy + abcz*bcdy*bcdz*bcx*bx*cx*cy + 
     -    az*bc2*bcdx*bcdy*by*cx*cy + ax*bc2*bcdy*bcdz*by*cx*cy + 
     -    abcz*bcdx*bcdz*bcx*by*cx*cy + bc2*bcdy*bcdz*bx*by*cx*cy - 
     -    bc2*bcdx*bcdz*by**2*cx*cy + az*bc2*bcdx*bcdz*bz*cx*cy + 
     -    ax*bc2*bcdz**2*bz*cx*cy + bc2*bcdz**2*bx*bz*cx*cy - 
     -    bc2*bcdx*bcdz*bz**2*cx*cy - abcz*ay*bc2*bcdx*cdy*cx*cy - 
     -    abcz*ax*bc2*bcdy*cdy*cx*cy + abcz*bc2*bcdy*bx*cdy*cx*cy + 
     -    abcz*bc2*bcdx*by*cdy*cx*cy - abcz*ax*bc2*bcdz*cdz*cx*cy + 
     -    abcz*bc2*bcdz*bx*cdz*cx*cy + ay*bc2*bcdy*bcdz*cx**2*cy + 
     -    az*bc2*bcdz**2*cx**2*cy - bc2*bcdy*bcdz*by*cx**2*cy - 
     -    bc2*bcdz**2*bz*cx**2*cy + ax*az*bc2*bcdx*bcdy*cy**2 + 
     -    abcz*bc2*bcdx*bcdz*cy**2 - ax*ay*bc2*bcdx*bcdz*cy**2 + 
     -    az**2*bc2*bcdy*bcdz*cy**2 - ay*az*bc2*bcdz**2*cy**2 + 
     -    abcz*ax*bcdx*bcdz*bcx*cy**2 + abcz*az*bcdz**2*bcx*cy**2 - 
     -    az*bc2*bcdx*bcdy*bx*cy**2 + 2*ay*bc2*bcdx*bcdz*bx*cy**2 + 
     -    ax*bc2*bcdy*bcdz*bx*cy**2 - abcz*bcdx*bcdz*bcx*bx*cy**2 - 
     -    bc2*bcdy*bcdz*bx**2*cy**2 - 2*ax*bc2*bcdx*bcdz*by*cy**2 - 
     -    2*az*bc2*bcdz**2*by*cy**2 + bc2*bcdx*bcdz*bx*by*cy**2 + 
     -    2*ay*bc2*bcdz**2*bz*cy**2 - abcz*bcdz**2*bcx*bz*cy**2 + 
     -    bc2*bcdz**2*by*bz*cy**2 - bc2*bcdy*bcdz*bz**2*cy**2 + 
     -    abcz*ax*bc2*bcdx*cdy*cy**2 - abcz*bc2*bcdx*bx*cdy*cy**2 - 
     -    ay*bc2*bcdx*bcdz*cx*cy**2 - ax*bc2*bcdy*bcdz*cx*cy**2 + 
     -    bc2*bcdy*bcdz*bx*cx*cy**2 + bc2*bcdx*bcdz*by*cx*cy**2 + 
     -    ax*bc2*bcdx*bcdz*cy**3 + az*bc2*bcdz**2*cy**3 - 
     -    bc2*bcdx*bcdz*bx*cy**3 - bc2*bcdz**2*bz*cy**3 - 
     -    abcz*az*bc2*bcdx*bcdz*cz - abcz*ax*bc2*bcdz**2*cz + 
     -    az**2*bc2*bcdx*bcdy*bx*cz - ay*az*bc2*bcdx*bcdz*bx*cz + 
     -    ax*az*bc2*bcdy*bcdz*bx*cz + 2*abcz*bc2*bcdz**2*bx*cz - 
     -    ax*ay*bc2*bcdz**2*bx*cz + abcz*az*bcdx*bcdz*bcx*bx*cz + 
     -    abcz*ax*bcdz**2*bcx*bx*cz - ay*bc2*bcdy**2*bx**2*cz - 
     -    2*az*bc2*bcdy*bcdz*bx**2*cz + ay*bc2*bcdz**2*bx**2*cz - 
     -    abcz*bcdz**2*bcx*bx**2*cz + az**2*bc2*bcdy**2*by*cz - 
     -    ay**2*bc2*bcdz**2*by*cz + abcz*az*bcdy*bcdz*bcx*by*cz + 
     -    abcz*ay*bcdz**2*bcx*by*cz + ay*bc2*bcdx*bcdy*bx*by*cz + 
     -    ax*bc2*bcdy**2*bx*by*cz - ax*bc2*bcdx*bcdy*by**2*cz - 
     -    2*az*bc2*bcdy*bcdz*by**2*cz + ay*bc2*bcdz**2*by**2*cz - 
     -    abcz*bcdz**2*bcx*by**2*cz - 2*ax*az*bc2*bcdx*bcdy*bz*cz - 
     -    2*ay*az*bc2*bcdy**2*bz*cz - abcz*bc2*bcdx*bcdz*bz*cz + 
     -    2*ax*ay*bc2*bcdx*bcdz*bz*cz + 2*ay**2*bc2*bcdy*bcdz*bz*cz - 
     -    2*abcz*ax*bcdx*bcdz*bcx*bz*cz - 
     -    2*abcz*ay*bcdy*bcdz*bcx*bz*cz + 2*az*bc2*bcdx*bcdy*bx*bz*cz - 
     -    ay*bc2*bcdx*bcdz*bx*bz*cz + ax*bc2*bcdy*bcdz*bx*bz*cz + 
     -    abcz*bcdx*bcdz*bcx*bx*bz*cz + 2*az*bc2*bcdy**2*by*bz*cz + 
     -    abcz*bcdy*bcdz*bcx*by*bz*cz - ax*bc2*bcdx*bcdy*bz**2*cz - 
     -    ay*bc2*bcdy**2*bz**2*cz + abcz*az*bc2*bcdx*bx*cdy*cz + 
     -    abcz*az*bc2*bcdy*by*cdy*cz - 2*abcz*ax*bc2*bcdx*bz*cdy*cz - 
     -    2*abcz*ay*bc2*bcdy*bz*cdy*cz + abcz*bc2*bcdx*bx*bz*cdy*cz + 
     -    abcz*bc2*bcdy*by*bz*cdy*cz + abcz*az*bc2*bcdz*by*cdz*cz - 
     -    2*abcz*ay*bc2*bcdz*bz*cdz*cz + abcz*bc2*bcdz*by*bz*cdz*cz - 
     -    az**2*bc2*bcdx*bcdy*cx*cz + ay*az*bc2*bcdx*bcdz*cx*cz - 
     -    ax*az*bc2*bcdy*bcdz*cx*cz - abcz*bc2*bcdz**2*cx*cz + 
     -    ax*ay*bc2*bcdz**2*cx*cz - abcz*az*bcdx*bcdz*bcx*cx*cz - 
     -    abcz*ax*bcdz**2*bcx*cx*cz + 2*ay*bc2*bcdy**2*bx*cx*cz + 
     -    3*az*bc2*bcdy*bcdz*bx*cx*cz - ay*bc2*bcdz**2*bx*cx*cz + 
     -    abcz*bcdz**2*bcx*bx*cx*cz - ay*bc2*bcdx*bcdy*by*cx*cz - 
     -    ax*bc2*bcdy**2*by*cx*cz - bc2*bcdy**2*bx*by*cx*cz + 
     -    bc2*bcdx*bcdy*by**2*cx*cz - ay*bc2*bcdx*bcdz*bz*cx*cz - 
     -    ax*bc2*bcdy*bcdz*bz*cx*cz + abcz*bcdx*bcdz*bcx*bz*cx*cz - 
     -    bc2*bcdy*bcdz*bx*bz*cx*cz + bc2*bcdx*bcdy*bz**2*cx*cz - 
     -    abcz*az*bc2*bcdx*cdy*cx*cz + abcz*bc2*bcdx*bz*cdy*cx*cz - 
     -    ay*bc2*bcdy**2*cx**2*cz - az*bc2*bcdy*bcdz*cx**2*cz + 
     -    bc2*bcdy**2*by*cx**2*cz + bc2*bcdy*bcdz*bz*cx**2*cz - 
     -    az**2*bc2*bcdy**2*cy*cz + ay**2*bc2*bcdz**2*cy*cz - 
     -    abcz*az*bcdy*bcdz*bcx*cy*cz - abcz*ay*bcdz**2*bcx*cy*cz - 
     -    ay*bc2*bcdx*bcdy*bx*cy*cz - ax*bc2*bcdy**2*bx*cy*cz + 
     -    az*bc2*bcdx*bcdz*bx*cy*cz + ax*bc2*bcdz**2*bx*cy*cz + 
     -    bc2*bcdy**2*bx**2*cy*cz - bc2*bcdz**2*bx**2*cy*cz + 
     -    2*ax*bc2*bcdx*bcdy*by*cy*cz + 4*az*bc2*bcdy*bcdz*by*cy*cz + 
     -    abcz*bcdz**2*bcx*by*cy*cz - bc2*bcdx*bcdy*bx*by*cy*cz - 
     -    bc2*bcdz**2*by**2*cy*cz - 2*ax*bc2*bcdx*bcdz*bz*cy*cz - 
     -    4*ay*bc2*bcdy*bcdz*bz*cy*cz + abcz*bcdy*bcdz*bcx*bz*cy*cz + 
     -    bc2*bcdx*bcdz*bx*bz*cy*cz + bc2*bcdy**2*bz**2*cy*cz - 
     -    abcz*az*bc2*bcdy*cdy*cy*cz + abcz*bc2*bcdy*bz*cdy*cy*cz - 
     -    abcz*az*bc2*bcdz*cdz*cy*cz + abcz*bc2*bcdz*bz*cdz*cy*cz + 
     -    ay*bc2*bcdx*bcdy*cx*cy*cz + ax*bc2*bcdy**2*cx*cy*cz - 
     -    az*bc2*bcdx*bcdz*cx*cy*cz - ax*bc2*bcdz**2*cx*cy*cz - 
     -    bc2*bcdy**2*bx*cx*cy*cz + bc2*bcdz**2*bx*cx*cy*cz - 
     -    bc2*bcdx*bcdy*by*cx*cy*cz + bc2*bcdx*bcdz*bz*cx*cy*cz - 
     -    ax*bc2*bcdx*bcdy*cy**2*cz - 2*az*bc2*bcdy*bcdz*cy**2*cz - 
     -    ay*bc2*bcdz**2*cy**2*cz + bc2*bcdx*bcdy*bx*cy**2*cz + 
     -    bc2*bcdz**2*by*cy**2*cz + 2*bc2*bcdy*bcdz*bz*cy**2*cz + 
     -    ax*az*bc2*bcdx*bcdy*cz**2 + ay*az*bc2*bcdy**2*cz**2 + 
     -    abcz*bc2*bcdx*bcdz*cz**2 - ax*ay*bc2*bcdx*bcdz*cz**2 - 
     -    ay**2*bc2*bcdy*bcdz*cz**2 + abcz*ax*bcdx*bcdz*bcx*cz**2 + 
     -    abcz*ay*bcdy*bcdz*bcx*cz**2 - 2*az*bc2*bcdx*bcdy*bx*cz**2 + 
     -    ay*bc2*bcdx*bcdz*bx*cz**2 - ax*bc2*bcdy*bcdz*bx*cz**2 - 
     -    abcz*bcdx*bcdz*bcx*bx*cz**2 + bc2*bcdy*bcdz*bx**2*cz**2 - 
     -    2*az*bc2*bcdy**2*by*cz**2 - abcz*bcdy*bcdz*bcx*by*cz**2 + 
     -    bc2*bcdy*bcdz*by**2*cz**2 + 2*ax*bc2*bcdx*bcdy*bz*cz**2 + 
     -    2*ay*bc2*bcdy**2*bz*cz**2 - bc2*bcdx*bcdy*bx*bz*cz**2 - 
     -    bc2*bcdy**2*by*bz*cz**2 + abcz*ax*bc2*bcdx*cdy*cz**2 + 
     -    abcz*ay*bc2*bcdy*cdy*cz**2 - abcz*bc2*bcdx*bx*cdy*cz**2 - 
     -    abcz*bc2*bcdy*by*cdy*cz**2 + abcz*ay*bc2*bcdz*cdz*cz**2 - 
     -    abcz*bc2*bcdz*by*cdz*cz**2 + az*bc2*bcdx*bcdy*cx*cz**2 + 
     -    ax*bc2*bcdy*bcdz*cx*cz**2 - bc2*bcdy*bcdz*bx*cx*cz**2 - 
     -    bc2*bcdx*bcdy*bz*cx*cz**2 + az*bc2*bcdy**2*cy*cz**2 + 
     -    ax*bc2*bcdx*bcdz*cy*cz**2 + 2*ay*bc2*bcdy*bcdz*cy*cz**2 - 
     -    bc2*bcdx*bcdz*bx*cy*cz**2 - 2*bc2*bcdy*bcdz*by*cy*cz**2 - 
     -    bc2*bcdy**2*bz*cy*cz**2 - ax*bc2*bcdx*bcdy*cz**3 - 
     -    ay*bc2*bcdy**2*cz**3 + bc2*bcdx*bcdy*bx*cz**3 + 
     -    bc2*bcdy**2*by*cz**3 + 
     -    abcx*bcdx*(ax*bc2*bcdy*by - ax*bcdy*bcx*bx*by + 
     -       ax*bcdx*bcx*by**2 + ax*bc2*bcdz*bz - ax*bcdz*bcx*bx*bz + 
     -       ax*bcdx*bcx*bz**2 + ax*bc2*bx*bz*cdy - ax*bc2*bx*by*cdz - 
     -       bc2*bcdy*by*cx + ax*bcdy*bcx*by*cx + bcdy*bcx*bx*by*cx - 
     -       bcdx*bcx*by**2*cx - bc2*bcdz*bz*cx + ax*bcdz*bcx*bz*cx + 
     -       bcdz*bcx*bx*bz*cx - bcdx*bcx*bz**2*cx - ax*bc2*bz*cdy*cx - 
     -       bc2*bx*bz*cdy*cx + ax*bc2*by*cdz*cx + bc2*bx*by*cdz*cx - 
     -       bcdy*bcx*by*cx**2 - bcdz*bcx*bz*cx**2 + bc2*bz*cdy*cx**2 - 
     -       bc2*by*cdz*cx**2 - ax*bc2*bcdy*cy + 2*bc2*bcdy*bx*cy + 
     -       ax*bcdy*bcx*bx*cy - bcdy*bcx*bx**2*cy - bc2*bcdx*by*cy - 
     -       2*ax*bcdx*bcx*by*cy + bcdx*bcx*bx*by*cy + 
     -       bcdz*bcx*by*bz*cy - bcdy*bcx*bz**2*cy - bc2*by*bz*cdy*cy + 
     -       ax*bc2*bx*cdz*cy - bc2*bx**2*cdz*cy - bc2*bz**2*cdz*cy - 
     -       bc2*bcdy*cx*cy - ax*bcdy*bcx*cx*cy + bcdy*bcx*bx*cx*cy + 
     -       bcdx*bcx*by*cx*cy - ax*bc2*cdz*cx*cy + bc2*bx*cdz*cx*cy + 
     -       bc2*bcdx*cy**2 + ax*bcdx*bcx*cy**2 - bcdx*bcx*bx*cy**2 - 
     -       bcdz*bcx*bz*cy**2 + bc2*bz*cdy*cy**2 - ax*bc2*bcdz*cz + 
     -       2*bc2*bcdz*bx*cz + ax*bcdz*bcx*bx*cz - bcdz*bcx*bx**2*cz - 
     -       bcdz*bcx*by**2*cz - bc2*bcdx*bz*cz - 2*ax*bcdx*bcx*bz*cz + 
     -       bcdx*bcx*bx*bz*cz + bcdy*bcx*by*bz*cz - ax*bc2*bx*cdy*cz + 
     -       bc2*bx**2*cdy*cz + bc2*by**2*cdy*cz + bc2*by*bz*cdz*cz - 
     -       bc2*bcdz*cx*cz - ax*bcdz*bcx*cx*cz + bcdz*bcx*bx*cx*cz + 
     -       bcdx*bcx*bz*cx*cz + ax*bc2*cdy*cx*cz - bc2*bx*cdy*cx*cz + 
     -       bcdz*bcx*by*cy*cz + bcdy*bcx*bz*cy*cz - bc2*by*cdy*cy*cz + 
     -       bc2*bz*cdz*cy*cz + bc2*bcdx*cz**2 + ax*bcdx*bcx*cz**2 - 
     -       bcdx*bcx*bx*cz**2 - bcdy*bcx*by*cz**2 - bc2*by*cdz*cz**2 + 
     -       az*(bcx*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(bcdx*(-bx + cx) + bcdy*(-by + cy))) - 
     -          bc2*(2*abx*bcdz + bx**2*cdy + by**2*cdy + by*bz*cdz - 
     -             2*bx*cdy*cx + cdy*cx**2 - 2*by*cdy*cy - bz*cdz*cy + 
     -             cdy*cy**2 - by*cdz*cz + cdz*cy*cz + bcdx*(-bz + cz)))
     -        + ay*(bc2*(-2*bcdy*bx + bcdx*by + by*bz*cdy + bx**2*cdz + 
     -             bz**2*cdz + 2*bcdy*cx - 2*bx*cdz*cx + cdz*cx**2 - 
     -             bcdx*cy - bz*cdy*cy - by*cdy*cz - 2*bz*cdz*cz + 
     -             cdy*cy*cz + cdz*cz**2) + 
     -          bcx*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(bcdx*(-bx + cx) + bcdz*(-bz + cz))))) - 
     -    abcy*(-(ax*bc2*bcdy**2*by) + ax*bcdy**2*bcx*bx*by - 
     -       ax*bcdx*bcdy*bcx*by**2 - ax*bc2*bcdy*bcdz*bz + 
     -       ax*bcdy*bcdz*bcx*bx*bz - ax*bcdx*bcdy*bcx*bz**2 - 
     -       ax*bc2*bcdy*bx*bz*cdy + ax*bc2*bcdx*by**2*cdz - 
     -       ax*bc2*bcdz*bx*bz*cdz + ax*bc2*bcdx*bz**2*cdz + 
     -       az*(bcdy*bcx*(-(bcdz*(bcx**2 + bcy**2)) + 
     -             (abx*bcdx + bcdy*bcy)*bcz) + 
     -          bc2*(-(bcdx*bcdy*bcz) + 
     -             (bcx**2 + bcy**2)*(bcdy*cdy + bcdz*cdz) + 
     -             abx*(2*bcdy*bcdz - bcdx*bcz*cdz))) + 
     -       bc2*bcdy**2*by*cx - ax*bcdy**2*bcx*by*cx - 
     -       bcdy**2*bcx*bx*by*cx + bcdx*bcdy*bcx*by**2*cx + 
     -       bc2*bcdy*bcdz*bz*cx - ax*bcdy*bcdz*bcx*bz*cx - 
     -       bcdy*bcdz*bcx*bx*bz*cx + bcdx*bcdy*bcx*bz**2*cx + 
     -       ax*bc2*bcdy*bz*cdy*cx + bc2*bcdy*bx*bz*cdy*cx - 
     -       bc2*bcdx*by**2*cdz*cx + ax*bc2*bcdz*bz*cdz*cx + 
     -       bc2*bcdz*bx*bz*cdz*cx - bc2*bcdx*bz**2*cdz*cx + 
     -       bcdy**2*bcx*by*cx**2 + bcdy*bcdz*bcx*bz*cx**2 - 
     -       bc2*bcdy*bz*cdy*cx**2 - bc2*bcdz*bz*cdz*cx**2 + 
     -       ax*bc2*bcdy**2*cy - 2*bc2*bcdy**2*bx*cy - 
     -       ax*bcdy**2*bcx*bx*cy + bcdy**2*bcx*bx**2*cy + 
     -       bc2*bcdx*bcdy*by*cy + 2*ax*bcdx*bcdy*bcx*by*cy - 
     -       bcdx*bcdy*bcx*bx*by*cy - bcdy*bcdz*bcx*by*bz*cy + 
     -       bcdy**2*bcx*bz**2*cy + bc2*bcdy*by*bz*cdy*cy - 
     -       2*ax*bc2*bcdx*by*cdz*cy + bc2*bcdx*bx*by*cdz*cy + 
     -       bc2*bcdz*by*bz*cdz*cy + bc2*bcdy**2*cx*cy + 
     -       ax*bcdy**2*bcx*cx*cy - bcdy**2*bcx*bx*cx*cy - 
     -       bcdx*bcdy*bcx*by*cx*cy + bc2*bcdx*by*cdz*cx*cy - 
     -       bc2*bcdx*bcdy*cy**2 - ax*bcdx*bcdy*bcx*cy**2 + 
     -       bcdx*bcdy*bcx*bx*cy**2 + bcdy*bcdz*bcx*bz*cy**2 - 
     -       bc2*bcdy*bz*cdy*cy**2 + ax*bc2*bcdx*cdz*cy**2 - 
     -       bc2*bcdx*bx*cdz*cy**2 - bc2*bcdz*bz*cdz*cy**2 + 
     -       ax*bc2*bcdy*bcdz*cz - 2*bc2*bcdy*bcdz*bx*cz - 
     -       ax*bcdy*bcdz*bcx*bx*cz + bcdy*bcdz*bcx*bx**2*cz + 
     -       bcdy*bcdz*bcx*by**2*cz + bc2*bcdx*bcdy*bz*cz + 
     -       2*ax*bcdx*bcdy*bcx*bz*cz - bcdx*bcdy*bcx*bx*bz*cz - 
     -       bcdy**2*bcx*by*bz*cz + ax*bc2*bcdy*bx*cdy*cz - 
     -       bc2*bcdy*bx**2*cdy*cz - bc2*bcdy*by**2*cdy*cz + 
     -       ax*bc2*bcdz*bx*cdz*cz - bc2*bcdz*bx**2*cdz*cz - 
     -       bc2*bcdz*by**2*cdz*cz - 2*ax*bc2*bcdx*bz*cdz*cz + 
     -       bc2*bcdx*bx*bz*cdz*cz + bc2*bcdy*bcdz*cx*cz + 
     -       ax*bcdy*bcdz*bcx*cx*cz - bcdy*bcdz*bcx*bx*cx*cz - 
     -       bcdx*bcdy*bcx*bz*cx*cz - ax*bc2*bcdy*cdy*cx*cz + 
     -       bc2*bcdy*bx*cdy*cx*cz - ax*bc2*bcdz*cdz*cx*cz + 
     -       bc2*bcdz*bx*cdz*cx*cz + bc2*bcdx*bz*cdz*cx*cz - 
     -       bcdy*bcdz*bcx*by*cy*cz - bcdy**2*bcx*bz*cy*cz + 
     -       bc2*bcdy*by*cdy*cy*cz + bc2*bcdz*by*cdz*cy*cz - 
     -       bc2*bcdx*bcdy*cz**2 - ax*bcdx*bcdy*bcx*cz**2 + 
     -       bcdx*bcdy*bcx*bx*cz**2 + bcdy**2*bcx*by*cz**2 + 
     -       ax*bc2*bcdx*cdz*cz**2 - bc2*bcdx*bx*cdz*cz**2 + 
     -       ay*(bc2*(2*abx*bcdy**2 - 
     -             bcy*(bcdx*(bcdy + bx*cdz - cdz*cx) + 
     -                (bcdy*cdy + bcdz*cdz)*(bz - cz))) - 
     -          bcdy*bcx*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(bcdx*(-bx + cx) + bcdz*(-bz + cz))))))/
     -  (bc2**1.5*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**2/
     -       bc2))


       dby= (abcz*ay*bc2*bcdx*bcdz*bx + abcz*ax*bc2*bcdy*bcdz*bx - 
     -    ay*az*bc2*bcdx*bcdy*bx**2 - az**2*bc2*bcdx*bcdz*bx**2 + 
     -    ax*ay*bc2*bcdy*bcdz*bx**2 + ax*az*bc2*bcdz**2*bx**2 + 
     -    abcz*ay*bcdy*bcdz*bcy*bx**2 + abcz*az*bcdz**2*bcy*bx**2 - 
     -    2*abcz*ax*bc2*bcdx*bcdz*by - 2*abcz*az*bc2*bcdz**2*by + 
     -    ay*az*bc2*bcdx**2*bx*by + ax*az*bc2*bcdx*bcdy*bx*by - 
     -    ax*ay*bc2*bcdx*bcdz*bx*by - ax**2*bc2*bcdy*bcdz*bx*by - 
     -    abcz*ay*bcdx*bcdz*bcy*bx*by - abcz*ax*bcdy*bcdz*bcy*bx*by - 
     -    ax*az*bc2*bcdx**2*by**2 + ax**2*bc2*bcdx*bcdz*by**2 - 
     -    az**2*bc2*bcdx*bcdz*by**2 + ax*az*bc2*bcdz**2*by**2 + 
     -    abcz*ax*bcdx*bcdz*bcy*by**2 + abcz*az*bcdz**2*bcy*by**2 + 
     -    abcz*az*bc2*bcdy*bcdz*bz + abcz*ay*bc2*bcdz**2*bz + 
     -    az**2*bc2*bcdx**2*bx*bz - ax**2*bc2*bcdz**2*bx*bz - 
     -    abcz*az*bcdx*bcdz*bcy*bx*bz - abcz*ax*bcdz**2*bcy*bx*bz + 
     -    az**2*bc2*bcdx*bcdy*by*bz + ay*az*bc2*bcdx*bcdz*by*bz - 
     -    ax*az*bc2*bcdy*bcdz*by*bz - ax*ay*bc2*bcdz**2*by*bz - 
     -    abcz*az*bcdy*bcdz*bcy*by*bz - abcz*ay*bcdz**2*bcy*by*bz - 
     -    ax*az*bc2*bcdx**2*bz**2 - ay*az*bc2*bcdx*bcdy*bz**2 + 
     -    ax**2*bc2*bcdx*bcdz*bz**2 + ax*ay*bc2*bcdy*bcdz*bz**2 + 
     -    abcz*ax*bcdx*bcdz*bcy*bz**2 + abcz*ay*bcdy*bcdz*bcy*bz**2 - 
     -    abcz*ay*bc2*bcdy*bx**2*cdx + abcz*ay*bc2*bcdx*bx*by*cdx + 
     -    abcz*ax*bc2*bcdy*bx*by*cdx - abcz*ax*bc2*bcdx*by**2*cdx + 
     -    abcz*az*bc2*bcdx*bx*bz*cdx + abcz*az*bc2*bcdy*by*bz*cdx - 
     -    abcz*ax*bc2*bcdx*bz**2*cdx - abcz*ay*bc2*bcdy*bz**2*cdx + 
     -    abcz*ay*bc2*bcdz*bx*by*cdz - abcz*ax*bc2*bcdz*by**2*cdz + 
     -    abcz*az*bc2*bcdz*bx*bz*cdz - abcz*ax*bc2*bcdz*bz**2*cdz - 
     -    abcz*ay*bc2*bcdx*bcdz*cx - abcz*ax*bc2*bcdy*bcdz*cx + 
     -    2*ay*az*bc2*bcdx*bcdy*bx*cx + 2*az**2*bc2*bcdx*bcdz*bx*cx - 
     -    abcz*bc2*bcdy*bcdz*bx*cx - 2*ax*ay*bc2*bcdy*bcdz*bx*cx - 
     -    2*ax*az*bc2*bcdz**2*bx*cx - 2*abcz*ay*bcdy*bcdz*bcy*bx*cx - 
     -    2*abcz*az*bcdz**2*bcy*bx*cx - ay*bc2*bcdy*bcdz*bx**2*cx - 
     -    az*bc2*bcdz**2*bx**2*cx - ay*az*bc2*bcdx**2*by*cx - 
     -    ax*az*bc2*bcdx*bcdy*by*cx + 2*abcz*bc2*bcdx*bcdz*by*cx + 
     -    ax*ay*bc2*bcdx*bcdz*by*cx + ax**2*bc2*bcdy*bcdz*by*cx + 
     -    abcz*ay*bcdx*bcdz*bcy*by*cx + abcz*ax*bcdy*bcdz*bcy*by*cx - 
     -    az*bc2*bcdx*bcdy*bx*by*cx + ay*bc2*bcdx*bcdz*bx*by*cx + 
     -    2*ax*bc2*bcdy*bcdz*bx*by*cx + abcz*bcdy*bcdz*bcy*bx*by*cx + 
     -    az*bc2*bcdx**2*by**2*cx - 2*ax*bc2*bcdx*bcdz*by**2*cx - 
     -    az*bc2*bcdz**2*by**2*cx - abcz*bcdx*bcdz*bcy*by**2*cx - 
     -    az**2*bc2*bcdx**2*bz*cx + ax**2*bc2*bcdz**2*bz*cx + 
     -    abcz*az*bcdx*bcdz*bcy*bz*cx + abcz*ax*bcdz**2*bcy*bz*cx + 
     -    2*ax*bc2*bcdz**2*bx*bz*cx + abcz*bcdz**2*bcy*bx*bz*cx + 
     -    az*bc2*bcdy*bcdz*by*bz*cx + ay*bc2*bcdz**2*by*bz*cx + 
     -    az*bc2*bcdx**2*bz**2*cx - 2*ax*bc2*bcdx*bcdz*bz**2*cx - 
     -    ay*bc2*bcdy*bcdz*bz**2*cx - abcz*bcdx*bcdz*bcy*bz**2*cx + 
     -    2*abcz*ay*bc2*bcdy*bx*cdx*cx - abcz*ay*bc2*bcdx*by*cdx*cx - 
     -    abcz*ax*bc2*bcdy*by*cdx*cx - abcz*bc2*bcdy*bx*by*cdx*cx + 
     -    abcz*bc2*bcdx*by**2*cdx*cx - abcz*az*bc2*bcdx*bz*cdx*cx + 
     -    abcz*bc2*bcdx*bz**2*cdx*cx - abcz*ay*bc2*bcdz*by*cdz*cx + 
     -    abcz*bc2*bcdz*by**2*cdz*cx - abcz*az*bc2*bcdz*bz*cdz*cx + 
     -    abcz*bc2*bcdz*bz**2*cdz*cx - ay*az*bc2*bcdx*bcdy*cx**2 - 
     -    az**2*bc2*bcdx*bcdz*cx**2 + abcz*bc2*bcdy*bcdz*cx**2 + 
     -    ax*ay*bc2*bcdy*bcdz*cx**2 + ax*az*bc2*bcdz**2*cx**2 + 
     -    abcz*ay*bcdy*bcdz*bcy*cx**2 + abcz*az*bcdz**2*bcy*cx**2 + 
     -    2*ay*bc2*bcdy*bcdz*bx*cx**2 + 2*az*bc2*bcdz**2*bx*cx**2 + 
     -    az*bc2*bcdx*bcdy*by*cx**2 - ay*bc2*bcdx*bcdz*by*cx**2 - 
     -    2*ax*bc2*bcdy*bcdz*by*cx**2 - abcz*bcdy*bcdz*bcy*by*cx**2 - 
     -    bc2*bcdy*bcdz*bx*by*cx**2 + bc2*bcdx*bcdz*by**2*cx**2 - 
     -    2*ax*bc2*bcdz**2*bz*cx**2 - abcz*bcdz**2*bcy*bz*cx**2 - 
     -    bc2*bcdz**2*bx*bz*cx**2 + bc2*bcdx*bcdz*bz**2*cx**2 - 
     -    abcz*ay*bc2*bcdy*cdx*cx**2 + abcz*bc2*bcdy*by*cdx*cx**2 - 
     -    ay*bc2*bcdy*bcdz*cx**3 - az*bc2*bcdz**2*cx**3 + 
     -    bc2*bcdy*bcdz*by*cx**3 + bc2*bcdz**2*bz*cx**3 + 
     -    2*abcz*ax*bc2*bcdx*bcdz*cy + 2*abcz*az*bc2*bcdz**2*cy - 
     -    ay*az*bc2*bcdx**2*bx*cy - ax*az*bc2*bcdx*bcdy*bx*cy - 
     -    abcz*bc2*bcdx*bcdz*bx*cy + ax*ay*bc2*bcdx*bcdz*bx*cy + 
     -    ax**2*bc2*bcdy*bcdz*bx*cy + abcz*ay*bcdx*bcdz*bcy*bx*cy + 
     -    abcz*ax*bcdy*bcdz*bcy*bx*cy + az*bc2*bcdx*bcdy*bx**2*cy - 
     -    ax*bc2*bcdy*bcdz*bx**2*cy - abcz*bcdy*bcdz*bcy*bx**2*cy + 
     -    2*ax*az*bc2*bcdx**2*by*cy - 2*ax**2*bc2*bcdx*bcdz*by*cy + 
     -    2*az**2*bc2*bcdx*bcdz*by*cy - 2*ax*az*bc2*bcdz**2*by*cy - 
     -    2*abcz*ax*bcdx*bcdz*bcy*by*cy - 2*abcz*az*bcdz**2*bcy*by*cy - 
     -    az*bc2*bcdx**2*bx*by*cy + ax*bc2*bcdx*bcdz*bx*by*cy + 
     -    abcz*bcdx*bcdz*bcy*bx*by*cy - az**2*bc2*bcdx*bcdy*bz*cy - 
     -    ay*az*bc2*bcdx*bcdz*bz*cy + ax*az*bc2*bcdy*bcdz*bz*cy - 
     -    abcz*bc2*bcdz**2*bz*cy + ax*ay*bc2*bcdz**2*bz*cy + 
     -    abcz*az*bcdy*bcdz*bcy*bz*cy + abcz*ay*bcdz**2*bcy*bz*cy - 
     -    az*bc2*bcdx*bcdz*by*bz*cy + ax*bc2*bcdz**2*by*bz*cy + 
     -    abcz*bcdz**2*bcy*by*bz*cy + az*bc2*bcdx*bcdy*bz**2*cy - 
     -    ax*bc2*bcdy*bcdz*bz**2*cy - abcz*bcdy*bcdz*bcy*bz**2*cy - 
     -    abcz*ay*bc2*bcdx*bx*cdx*cy - abcz*ax*bc2*bcdy*bx*cdx*cy + 
     -    abcz*bc2*bcdy*bx**2*cdx*cy + 2*abcz*ax*bc2*bcdx*by*cdx*cy - 
     -    abcz*bc2*bcdx*bx*by*cdx*cy - abcz*az*bc2*bcdy*bz*cdx*cy + 
     -    abcz*bc2*bcdy*bz**2*cdx*cy - abcz*ay*bc2*bcdz*bx*cdz*cy + 
     -    2*abcz*ax*bc2*bcdz*by*cdz*cy - abcz*bc2*bcdz*bx*by*cdz*cy + 
     -    ay*az*bc2*bcdx**2*cx*cy + ax*az*bc2*bcdx*bcdy*cx*cy - 
     -    abcz*bc2*bcdx*bcdz*cx*cy - ax*ay*bc2*bcdx*bcdz*cx*cy - 
     -    ax**2*bc2*bcdy*bcdz*cx*cy - abcz*ay*bcdx*bcdz*bcy*cx*cy - 
     -    abcz*ax*bcdy*bcdz*bcy*cx*cy - az*bc2*bcdx*bcdy*bx*cx*cy - 
     -    ay*bc2*bcdx*bcdz*bx*cx*cy + abcz*bcdy*bcdz*bcy*bx*cx*cy + 
     -    bc2*bcdy*bcdz*bx**2*cx*cy - az*bc2*bcdx**2*by*cx*cy + 
     -    3*ax*bc2*bcdx*bcdz*by*cx*cy + 2*az*bc2*bcdz**2*by*cx*cy + 
     -    abcz*bcdx*bcdz*bcy*by*cx*cy - bc2*bcdx*bcdz*bx*by*cx*cy - 
     -    az*bc2*bcdy*bcdz*bz*cx*cy - ay*bc2*bcdz**2*bz*cx*cy - 
     -    bc2*bcdz**2*by*bz*cx*cy + bc2*bcdy*bcdz*bz**2*cx*cy + 
     -    abcz*ay*bc2*bcdx*cdx*cx*cy + abcz*ax*bc2*bcdy*cdx*cx*cy - 
     -    abcz*bc2*bcdy*bx*cdx*cx*cy - abcz*bc2*bcdx*by*cdx*cx*cy + 
     -    abcz*ay*bc2*bcdz*cdz*cx*cy - abcz*bc2*bcdz*by*cdz*cx*cy + 
     -    ay*bc2*bcdx*bcdz*cx**2*cy + ax*bc2*bcdy*bcdz*cx**2*cy - 
     -    bc2*bcdy*bcdz*bx*cx**2*cy - bc2*bcdx*bcdz*by*cx**2*cy - 
     -    ax*az*bc2*bcdx**2*cy**2 + ax**2*bc2*bcdx*bcdz*cy**2 - 
     -    az**2*bc2*bcdx*bcdz*cy**2 + ax*az*bc2*bcdz**2*cy**2 + 
     -    abcz*ax*bcdx*bcdz*bcy*cy**2 + abcz*az*bcdz**2*bcy*cy**2 + 
     -    az*bc2*bcdx**2*bx*cy**2 - ax*bc2*bcdx*bcdz*bx*cy**2 - 
     -    abcz*bcdx*bcdz*bcy*bx*cy**2 + az*bc2*bcdx*bcdz*bz*cy**2 - 
     -    ax*bc2*bcdz**2*bz*cy**2 - abcz*bcdz**2*bcy*bz*cy**2 - 
     -    abcz*ax*bc2*bcdx*cdx*cy**2 + abcz*bc2*bcdx*bx*cdx*cy**2 - 
     -    abcz*ax*bc2*bcdz*cdz*cy**2 + abcz*bc2*bcdz*bx*cdz*cy**2 - 
     -    ax*bc2*bcdx*bcdz*cx*cy**2 - az*bc2*bcdz**2*cx*cy**2 + 
     -    bc2*bcdx*bcdz*bx*cx*cy**2 + bc2*bcdz**2*bz*cx*cy**2 - 
     -    abcz*az*bc2*bcdy*bcdz*cz - abcz*ay*bc2*bcdz**2*cz - 
     -    az**2*bc2*bcdx**2*bx*cz + ax**2*bc2*bcdz**2*bx*cz + 
     -    abcz*az*bcdx*bcdz*bcy*bx*cz + abcz*ax*bcdz**2*bcy*bx*cz + 
     -    ay*bc2*bcdx*bcdy*bx**2*cz + 2*az*bc2*bcdx*bcdz*bx**2*cz - 
     -    ax*bc2*bcdz**2*bx**2*cz - abcz*bcdz**2*bcy*bx**2*cz - 
     -    az**2*bc2*bcdx*bcdy*by*cz - ay*az*bc2*bcdx*bcdz*by*cz + 
     -    ax*az*bc2*bcdy*bcdz*by*cz + 2*abcz*bc2*bcdz**2*by*cz + 
     -    ax*ay*bc2*bcdz**2*by*cz + abcz*az*bcdy*bcdz*bcy*by*cz + 
     -    abcz*ay*bcdz**2*bcy*by*cz - ay*bc2*bcdx**2*bx*by*cz - 
     -    ax*bc2*bcdx*bcdy*bx*by*cz + ax*bc2*bcdx**2*by**2*cz + 
     -    2*az*bc2*bcdx*bcdz*by**2*cz - ax*bc2*bcdz**2*by**2*cz - 
     -    abcz*bcdz**2*bcy*by**2*cz + 2*ax*az*bc2*bcdx**2*bz*cz + 
     -    2*ay*az*bc2*bcdx*bcdy*bz*cz - 2*ax**2*bc2*bcdx*bcdz*bz*cz - 
     -    abcz*bc2*bcdy*bcdz*bz*cz - 2*ax*ay*bc2*bcdy*bcdz*bz*cz - 
     -    2*abcz*ax*bcdx*bcdz*bcy*bz*cz - 
     -    2*abcz*ay*bcdy*bcdz*bcy*bz*cz - 2*az*bc2*bcdx**2*bx*bz*cz + 
     -    abcz*bcdx*bcdz*bcy*bx*bz*cz - 2*az*bc2*bcdx*bcdy*by*bz*cz - 
     -    ay*bc2*bcdx*bcdz*by*bz*cz + ax*bc2*bcdy*bcdz*by*bz*cz + 
     -    abcz*bcdy*bcdz*bcy*by*bz*cz + ax*bc2*bcdx**2*bz**2*cz + 
     -    ay*bc2*bcdx*bcdy*bz**2*cz - abcz*az*bc2*bcdx*bx*cdx*cz - 
     -    abcz*az*bc2*bcdy*by*cdx*cz + 2*abcz*ax*bc2*bcdx*bz*cdx*cz + 
     -    2*abcz*ay*bc2*bcdy*bz*cdx*cz - abcz*bc2*bcdx*bx*bz*cdx*cz - 
     -    abcz*bc2*bcdy*by*bz*cdx*cz - abcz*az*bc2*bcdz*bx*cdz*cz + 
     -    2*abcz*ax*bc2*bcdz*bz*cdz*cz - abcz*bc2*bcdz*bx*bz*cdz*cz + 
     -    az**2*bc2*bcdx**2*cx*cz - ax**2*bc2*bcdz**2*cx*cz - 
     -    abcz*az*bcdx*bcdz*bcy*cx*cz - abcz*ax*bcdz**2*bcy*cx*cz - 
     -    2*ay*bc2*bcdx*bcdy*bx*cx*cz - 4*az*bc2*bcdx*bcdz*bx*cx*cz + 
     -    abcz*bcdz**2*bcy*bx*cx*cz + bc2*bcdz**2*bx**2*cx*cz + 
     -    ay*bc2*bcdx**2*by*cx*cz + ax*bc2*bcdx*bcdy*by*cx*cz - 
     -    az*bc2*bcdy*bcdz*by*cx*cz - ay*bc2*bcdz**2*by*cx*cz + 
     -    bc2*bcdx*bcdy*bx*by*cx*cz - bc2*bcdx**2*by**2*cx*cz + 
     -    bc2*bcdz**2*by**2*cx*cz + 4*ax*bc2*bcdx*bcdz*bz*cx*cz + 
     -    2*ay*bc2*bcdy*bcdz*bz*cx*cz + abcz*bcdx*bcdz*bcy*bz*cx*cz - 
     -    bc2*bcdy*bcdz*by*bz*cx*cz - bc2*bcdx**2*bz**2*cx*cz + 
     -    abcz*az*bc2*bcdx*cdx*cx*cz - abcz*bc2*bcdx*bz*cdx*cx*cz + 
     -    abcz*az*bc2*bcdz*cdz*cx*cz - abcz*bc2*bcdz*bz*cdz*cx*cz + 
     -    ay*bc2*bcdx*bcdy*cx**2*cz + 2*az*bc2*bcdx*bcdz*cx**2*cz + 
     -    ax*bc2*bcdz**2*cx**2*cz - bc2*bcdz**2*bx*cx**2*cz - 
     -    bc2*bcdx*bcdy*by*cx**2*cz - 2*bc2*bcdx*bcdz*bz*cx**2*cz + 
     -    az**2*bc2*bcdx*bcdy*cy*cz + ay*az*bc2*bcdx*bcdz*cy*cz - 
     -    ax*az*bc2*bcdy*bcdz*cy*cz - abcz*bc2*bcdz**2*cy*cz - 
     -    ax*ay*bc2*bcdz**2*cy*cz - abcz*az*bcdy*bcdz*bcy*cy*cz - 
     -    abcz*ay*bcdz**2*bcy*cy*cz + ay*bc2*bcdx**2*bx*cy*cz + 
     -    ax*bc2*bcdx*bcdy*bx*cy*cz - bc2*bcdx*bcdy*bx**2*cy*cz - 
     -    2*ax*bc2*bcdx**2*by*cy*cz - 3*az*bc2*bcdx*bcdz*by*cy*cz + 
     -    ax*bc2*bcdz**2*by*cy*cz + abcz*bcdz**2*bcy*by*cy*cz + 
     -    bc2*bcdx**2*bx*by*cy*cz + ay*bc2*bcdx*bcdz*bz*cy*cz + 
     -    ax*bc2*bcdy*bcdz*bz*cy*cz + abcz*bcdy*bcdz*bcy*bz*cy*cz + 
     -    bc2*bcdx*bcdz*by*bz*cy*cz - bc2*bcdx*bcdy*bz**2*cy*cz + 
     -    abcz*az*bc2*bcdy*cdx*cy*cz - abcz*bc2*bcdy*bz*cdx*cy*cz - 
     -    ay*bc2*bcdx**2*cx*cy*cz - ax*bc2*bcdx*bcdy*cx*cy*cz + 
     -    az*bc2*bcdy*bcdz*cx*cy*cz + ay*bc2*bcdz**2*cx*cy*cz + 
     -    bc2*bcdx*bcdy*bx*cx*cy*cz + bc2*bcdx**2*by*cx*cy*cz - 
     -    bc2*bcdz**2*by*cx*cy*cz - bc2*bcdy*bcdz*bz*cx*cy*cz + 
     -    ax*bc2*bcdx**2*cy**2*cz + az*bc2*bcdx*bcdz*cy**2*cz - 
     -    bc2*bcdx**2*bx*cy**2*cz - bc2*bcdx*bcdz*bz*cy**2*cz - 
     -    ax*az*bc2*bcdx**2*cz**2 - ay*az*bc2*bcdx*bcdy*cz**2 + 
     -    ax**2*bc2*bcdx*bcdz*cz**2 + abcz*bc2*bcdy*bcdz*cz**2 + 
     -    ax*ay*bc2*bcdy*bcdz*cz**2 + abcz*ax*bcdx*bcdz*bcy*cz**2 + 
     -    abcz*ay*bcdy*bcdz*bcy*cz**2 + 2*az*bc2*bcdx**2*bx*cz**2 - 
     -    abcz*bcdx*bcdz*bcy*bx*cz**2 - bc2*bcdx*bcdz*bx**2*cz**2 + 
     -    2*az*bc2*bcdx*bcdy*by*cz**2 + ay*bc2*bcdx*bcdz*by*cz**2 - 
     -    ax*bc2*bcdy*bcdz*by*cz**2 - abcz*bcdy*bcdz*bcy*by*cz**2 - 
     -    bc2*bcdx*bcdz*by**2*cz**2 - 2*ax*bc2*bcdx**2*bz*cz**2 - 
     -    2*ay*bc2*bcdx*bcdy*bz*cz**2 + bc2*bcdx**2*bx*bz*cz**2 + 
     -    bc2*bcdx*bcdy*by*bz*cz**2 - abcz*ax*bc2*bcdx*cdx*cz**2 - 
     -    abcz*ay*bc2*bcdy*cdx*cz**2 + abcz*bc2*bcdx*bx*cdx*cz**2 + 
     -    abcz*bc2*bcdy*by*cdx*cz**2 - abcz*ax*bc2*bcdz*cdz*cz**2 + 
     -    abcz*bc2*bcdz*bx*cdz*cz**2 - az*bc2*bcdx**2*cx*cz**2 - 
     -    2*ax*bc2*bcdx*bcdz*cx*cz**2 - ay*bc2*bcdy*bcdz*cx*cz**2 + 
     -    2*bc2*bcdx*bcdz*bx*cx*cz**2 + bc2*bcdy*bcdz*by*cx*cz**2 + 
     -    bc2*bcdx**2*bz*cx*cz**2 - az*bc2*bcdx*bcdy*cy*cz**2 - 
     -    ay*bc2*bcdx*bcdz*cy*cz**2 + bc2*bcdx*bcdz*by*cy*cz**2 + 
     -    bc2*bcdx*bcdy*bz*cy*cz**2 + ax*bc2*bcdx**2*cz**3 + 
     -    ay*bc2*bcdx*bcdy*cz**3 - bc2*bcdx**2*bx*cz**3 - 
     -    bc2*bcdx*bcdy*by*cz**3 + 
     -    abcy*bcdy*(az*bcdz*bcy*bx**2 - 2*az*bc2*bcdz*by + 
     -       az*bcdz*bcy*by**2 + az*bc2*bcdy*bz - az*bcdx*bcy*bx*bz - 
     -       az*bcdy*bcy*by*bz + az*bc2*bx**2*cdx + az*bc2*by**2*cdx + 
     -       az*bc2*bx*bz*cdz - bc2*bcdy*bx*cx - 2*az*bcdz*bcy*bx*cx + 
     -       2*bc2*bcdx*by*cx + bcdy*bcy*bx*by*cx - bcdx*bcy*by**2*cx + 
     -       az*bcdx*bcy*bz*cx + bcdz*bcy*bx*bz*cx - 
     -       bcdx*bcy*bz**2*cx - 2*az*bc2*bx*cdx*cx + 
     -       bc2*bx*bz*cdx*cx + bc2*by**2*cdz*cx - az*bc2*bz*cdz*cx + 
     -       bc2*bz**2*cdz*cx + bc2*bcdy*cx**2 + az*bcdz*bcy*cx**2 - 
     -       bcdy*bcy*by*cx**2 - bcdz*bcy*bz*cx**2 + az*bc2*cdx*cx**2 - 
     -       bc2*bz*cdx*cx**2 + 2*az*bc2*bcdz*cy - bc2*bcdx*bx*cy - 
     -       bcdy*bcy*bx**2*cy - 2*az*bcdz*bcy*by*cy + 
     -       bcdx*bcy*bx*by*cy - bc2*bcdz*bz*cy + az*bcdy*bcy*bz*cy + 
     -       bcdz*bcy*by*bz*cy - bcdy*bcy*bz**2*cy - 
     -       2*az*bc2*by*cdx*cy + bc2*by*bz*cdx*cy - bc2*bx*by*cdz*cy - 
     -       bc2*bcdx*cx*cy + bcdy*bcy*bx*cx*cy + bcdx*bcy*by*cx*cy - 
     -       bc2*by*cdz*cx*cy + az*bcdz*bcy*cy**2 - bcdx*bcy*bx*cy**2 - 
     -       bcdz*bcy*bz*cy**2 + az*bc2*cdx*cy**2 - bc2*bz*cdx*cy**2 + 
     -       bc2*bx*cdz*cy**2 - az*bc2*bcdy*cz + az*bcdx*bcy*bx*cz - 
     -       bcdz*bcy*bx**2*cz + 2*bc2*bcdz*by*cz + az*bcdy*bcy*by*cz - 
     -       bcdz*bcy*by**2*cz - bc2*bcdy*bz*cz + bcdx*bcy*bx*bz*cz + 
     -       bcdy*bcy*by*bz*cz - bc2*bx**2*cdx*cz - bc2*by**2*cdx*cz - 
     -       az*bc2*bx*cdz*cz - bc2*bx*bz*cdz*cz - az*bcdx*bcy*cx*cz + 
     -       bcdz*bcy*bx*cx*cz + bcdx*bcy*bz*cx*cz + bc2*bx*cdx*cx*cz + 
     -       az*bc2*cdz*cx*cz - bc2*bz*cdz*cx*cz - bc2*bcdz*cy*cz - 
     -       az*bcdy*bcy*cy*cz + bcdz*bcy*by*cy*cz + 
     -       bcdy*bcy*bz*cy*cz + bc2*by*cdx*cy*cz + bc2*bcdy*cz**2 - 
     -       bcdx*bcy*bx*cz**2 - bcdy*bcy*by*cz**2 + bc2*bx*cdz*cz**2 + 
     -       ax*(-(abx*bcy*(bcdy*bcy + bcdz*bcz)) + 
     -          bcdx*(bcy**3 + bcy*bcz**2 - 2*bc2*by + 2*bc2*cy) + 
     -          bc2*(bcdy*bx - bx*bz*cdx - by**2*cdz - bz**2*cdz - 
     -             bcdy*cx + bz*cdx*cx + 2*by*cdz*cy - cdz*cy**2 + 
     -             bx*cdx*cz + 2*bz*cdz*cz - cdx*cx*cz - cdz*cz**2)) + 
     -       ay*(bcy*(-(bcy*(abx*bcdx + bcdz*bcz)) + 
     -             bcdy*(bcx**2 + bcz**2)) + 
     -          bc2*(abx*bcdx + bcdz*bcz + 
     -             bcy*(-(bz*cdx) + bx*cdz - cdz*cx + cdx*cz)))) + 
     -    abcx*(az*bcdx*bcdz*bcy*bx**2 - 2*az*bc2*bcdx*bcdz*by + 
     -       az*bcdx*bcdz*bcy*by**2 + az*bc2*bcdx*bcdy*bz - 
     -       az*bcdx**2*bcy*bx*bz - az*bcdx*bcdy*bcy*by*bz + 
     -       az*bc2*bcdx*bx**2*cdx + az*bc2*bcdx*by**2*cdx + 
     -       az*bc2*bcdz*bx**2*cdz + az*bc2*bcdz*by**2*cdz - 
     -       az*bc2*bcdy*by*bz*cdz + 
     -       ax*(bcdx**2*bcy*(-2*bc2 + bcy**2 + bcz**2) - 
     -          abx*(bcdx*bcy*(bcdy*bcy + bcdz*bcz) + 
     -             bc2*(-(bcdx*bcdy) + bcdx*bcz*cdx + bcdy*bcy*cdz + 
     -                bcdz*bcz*cdz))) - bc2*bcdx*bcdy*bx*cx - 
     -       2*az*bcdx*bcdz*bcy*bx*cx + 2*bc2*bcdx**2*by*cx + 
     -       bcdx*bcdy*bcy*bx*by*cx - bcdx**2*bcy*by**2*cx + 
     -       az*bcdx**2*bcy*bz*cx + bcdx*bcdz*bcy*bx*bz*cx - 
     -       bcdx**2*bcy*bz**2*cx - 2*az*bc2*bcdx*bx*cdx*cx + 
     -       bc2*bcdx*bx*bz*cdx*cx - 2*az*bc2*bcdz*bx*cdz*cx + 
     -       bc2*bcdy*bx*by*cdz*cx + bc2*bcdz*bx*bz*cdz*cx + 
     -       bc2*bcdx*bcdy*cx**2 + az*bcdx*bcdz*bcy*cx**2 - 
     -       bcdx*bcdy*bcy*by*cx**2 - bcdx*bcdz*bcy*bz*cx**2 + 
     -       az*bc2*bcdx*cdx*cx**2 - bc2*bcdx*bz*cdx*cx**2 + 
     -       az*bc2*bcdz*cdz*cx**2 - bc2*bcdy*by*cdz*cx**2 - 
     -       bc2*bcdz*bz*cdz*cx**2 + 2*az*bc2*bcdx*bcdz*cy - 
     -       bc2*bcdx**2*bx*cy - bcdx*bcdy*bcy*bx**2*cy - 
     -       2*az*bcdx*bcdz*bcy*by*cy + bcdx**2*bcy*bx*by*cy - 
     -       bc2*bcdx*bcdz*bz*cy + az*bcdx*bcdy*bcy*bz*cy + 
     -       bcdx*bcdz*bcy*by*bz*cy - bcdx*bcdy*bcy*bz**2*cy - 
     -       2*az*bc2*bcdx*by*cdx*cy + bc2*bcdx*by*bz*cdx*cy - 
     -       bc2*bcdy*bx**2*cdz*cy - 2*az*bc2*bcdz*by*cdz*cy + 
     -       az*bc2*bcdy*bz*cdz*cy + bc2*bcdz*by*bz*cdz*cy - 
     -       bc2*bcdy*bz**2*cdz*cy - bc2*bcdx**2*cx*cy + 
     -       bcdx*bcdy*bcy*bx*cx*cy + bcdx**2*bcy*by*cx*cy + 
     -       bc2*bcdy*bx*cdz*cx*cy + az*bcdx*bcdz*bcy*cy**2 - 
     -       bcdx**2*bcy*bx*cy**2 - bcdx*bcdz*bcy*bz*cy**2 + 
     -       az*bc2*bcdx*cdx*cy**2 - bc2*bcdx*bz*cdx*cy**2 + 
     -       az*bc2*bcdz*cdz*cy**2 - bc2*bcdz*bz*cdz*cy**2 - 
     -       az*bc2*bcdx*bcdy*cz + az*bcdx**2*bcy*bx*cz - 
     -       bcdx*bcdz*bcy*bx**2*cz + 2*bc2*bcdx*bcdz*by*cz + 
     -       az*bcdx*bcdy*bcy*by*cz - bcdx*bcdz*bcy*by**2*cz - 
     -       bc2*bcdx*bcdy*bz*cz + bcdx**2*bcy*bx*bz*cz + 
     -       bcdx*bcdy*bcy*by*bz*cz - bc2*bcdx*bx**2*cdx*cz - 
     -       bc2*bcdx*by**2*cdx*cz - bc2*bcdz*bx**2*cdz*cz + 
     -       az*bc2*bcdy*by*cdz*cz - bc2*bcdz*by**2*cdz*cz + 
     -       bc2*bcdy*by*bz*cdz*cz - az*bcdx**2*bcy*cx*cz + 
     -       bcdx*bcdz*bcy*bx*cx*cz + bcdx**2*bcy*bz*cx*cz + 
     -       bc2*bcdx*bx*cdx*cx*cz + bc2*bcdz*bx*cdz*cx*cz - 
     -       bc2*bcdx*bcdz*cy*cz - az*bcdx*bcdy*bcy*cy*cz + 
     -       bcdx*bcdz*bcy*by*cy*cz + bcdx*bcdy*bcy*bz*cy*cz + 
     -       bc2*bcdx*by*cdx*cy*cz - az*bc2*bcdy*cdz*cy*cz + 
     -       bc2*bcdz*by*cdz*cy*cz + bc2*bcdy*bz*cdz*cy*cz + 
     -       bc2*bcdx*bcdy*cz**2 - bcdx**2*bcy*bx*cz**2 - 
     -       bcdx*bcdy*bcy*by*cz**2 - bc2*bcdy*by*cdz*cz**2 + 
     -       ay*(bc2*(abx*bcdx**2 + bcdx*bcz*(bcdz - bcy*cdx) + 
     -             (-(bcdz*bcy*bcz) + bcdy*(bcx**2 + bcz**2))*cdz) + 
     -          bcdx*bcy*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(bcdx*(-bx + cx) + bcdz*(-bz + cz))))))/
     -  (bc2**1.5*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**2/
     -       bc2))

       dbz= (abcz*az*bc2*bcdx*bcdz*bx + abcz*ax*bc2*bcdz**2*bx + 
     -    ay**2*bc2*bcdx*bcdy*bx**2 - ax*ay*bc2*bcdy**2*bx**2 + 
     -    ay*az*bc2*bcdx*bcdz*bx**2 - ax*az*bc2*bcdy*bcdz*bx**2 + 
     -    abcz*ay*bcdy*bcdz*bcz*bx**2 + abcz*az*bcdz**2*bcz*bx**2 + 
     -    abcz*az*bc2*bcdy*bcdz*by + abcz*ay*bc2*bcdz**2*by - 
     -    ay**2*bc2*bcdx**2*bx*by + ax**2*bc2*bcdy**2*bx*by - 
     -    abcz*ay*bcdx*bcdz*bcz*bx*by - abcz*ax*bcdy*bcdz*bcz*bx*by + 
     -    ax*ay*bc2*bcdx**2*by**2 - ax**2*bc2*bcdx*bcdy*by**2 + 
     -    ay*az*bc2*bcdx*bcdz*by**2 - ax*az*bc2*bcdy*bcdz*by**2 + 
     -    abcz*ax*bcdx*bcdz*bcz*by**2 + abcz*az*bcdz**2*bcz*by**2 - 
     -    2*abcz*ax*bc2*bcdx*bcdz*bz - 2*abcz*ay*bc2*bcdy*bcdz*bz - 
     -    ay*az*bc2*bcdx**2*bx*bz + ax*az*bc2*bcdx*bcdy*bx*bz - 
     -    ax*ay*bc2*bcdx*bcdz*bx*bz + ax**2*bc2*bcdy*bcdz*bx*bz - 
     -    abcz*az*bcdx*bcdz*bcz*bx*bz - abcz*ax*bcdz**2*bcz*bx*bz - 
     -    ay*az*bc2*bcdx*bcdy*by*bz + ax*az*bc2*bcdy**2*by*bz - 
     -    ay**2*bc2*bcdx*bcdz*by*bz + ax*ay*bc2*bcdy*bcdz*by*bz - 
     -    abcz*az*bcdy*bcdz*bcz*by*bz - abcz*ay*bcdz**2*bcz*by*bz + 
     -    ax*ay*bc2*bcdx**2*bz**2 - ax**2*bc2*bcdx*bcdy*bz**2 + 
     -    ay**2*bc2*bcdx*bcdy*bz**2 - ax*ay*bc2*bcdy**2*bz**2 + 
     -    abcz*ax*bcdx*bcdz*bcz*bz**2 + abcz*ay*bcdy*bcdz*bcz*bz**2 - 
     -    abcz*ay*bc2*bcdz*bx**2*cdx + abcz*ax*bc2*bcdz*bx*by*cdx + 
     -    abcz*az*bc2*bcdz*by*bz*cdx - abcz*ay*bc2*bcdz*bz**2*cdx - 
     -    abcz*ay*bc2*bcdz*bx*by*cdy + abcz*ax*bc2*bcdz*by**2*cdy - 
     -    abcz*az*bc2*bcdz*bx*bz*cdy + abcz*ax*bc2*bcdz*bz**2*cdy - 
     -    abcz*az*bc2*bcdx*bcdz*cx - abcz*ax*bc2*bcdz**2*cx - 
     -    2*ay**2*bc2*bcdx*bcdy*bx*cx + 2*ax*ay*bc2*bcdy**2*bx*cx - 
     -    2*ay*az*bc2*bcdx*bcdz*bx*cx + 2*ax*az*bc2*bcdy*bcdz*bx*cx - 
     -    abcz*bc2*bcdz**2*bx*cx - 2*abcz*ay*bcdy*bcdz*bcz*bx*cx - 
     -    2*abcz*az*bcdz**2*bcz*bx*cx + ay*bc2*bcdy**2*bx**2*cx + 
     -    az*bc2*bcdy*bcdz*bx**2*cx + ay**2*bc2*bcdx**2*by*cx - 
     -    ax**2*bc2*bcdy**2*by*cx + abcz*ay*bcdx*bcdz*bcz*by*cx + 
     -    abcz*ax*bcdy*bcdz*bcz*by*cx - 2*ax*bc2*bcdy**2*bx*by*cx + 
     -    abcz*bcdy*bcdz*bcz*bx*by*cx - ay*bc2*bcdx**2*by**2*cx + 
     -    2*ax*bc2*bcdx*bcdy*by**2*cx + az*bc2*bcdy*bcdz*by**2*cx - 
     -    abcz*bcdx*bcdz*bcz*by**2*cx + ay*az*bc2*bcdx**2*bz*cx - 
     -    ax*az*bc2*bcdx*bcdy*bz*cx + 2*abcz*bc2*bcdx*bcdz*bz*cx + 
     -    ax*ay*bc2*bcdx*bcdz*bz*cx - ax**2*bc2*bcdy*bcdz*bz*cx + 
     -    abcz*az*bcdx*bcdz*bcz*bz*cx + abcz*ax*bcdz**2*bcz*bz*cx - 
     -    az*bc2*bcdx*bcdy*bx*bz*cx + ay*bc2*bcdx*bcdz*bx*bz*cx - 
     -    2*ax*bc2*bcdy*bcdz*bx*bz*cx + abcz*bcdz**2*bcz*bx*bz*cx - 
     -    az*bc2*bcdy**2*by*bz*cx - ay*bc2*bcdy*bcdz*by*bz*cx - 
     -    ay*bc2*bcdx**2*bz**2*cx + 2*ax*bc2*bcdx*bcdy*bz**2*cx + 
     -    ay*bc2*bcdy**2*bz**2*cx - abcz*bcdx*bcdz*bcz*bz**2*cx + 
     -    2*abcz*ay*bc2*bcdz*bx*cdx*cx - abcz*ax*bc2*bcdz*by*cdx*cx - 
     -    abcz*bc2*bcdz*bx*by*cdx*cx + abcz*ay*bc2*bcdz*by*cdy*cx - 
     -    abcz*bc2*bcdz*by**2*cdy*cx + abcz*az*bc2*bcdz*bz*cdy*cx - 
     -    abcz*bc2*bcdz*bz**2*cdy*cx + ay**2*bc2*bcdx*bcdy*cx**2 - 
     -    ax*ay*bc2*bcdy**2*cx**2 + ay*az*bc2*bcdx*bcdz*cx**2 - 
     -    ax*az*bc2*bcdy*bcdz*cx**2 + abcz*bc2*bcdz**2*cx**2 + 
     -    abcz*ay*bcdy*bcdz*bcz*cx**2 + abcz*az*bcdz**2*bcz*cx**2 - 
     -    2*ay*bc2*bcdy**2*bx*cx**2 - 2*az*bc2*bcdy*bcdz*bx*cx**2 + 
     -    2*ax*bc2*bcdy**2*by*cx**2 - abcz*bcdy*bcdz*bcz*by*cx**2 + 
     -    bc2*bcdy**2*bx*by*cx**2 - bc2*bcdx*bcdy*by**2*cx**2 + 
     -    az*bc2*bcdx*bcdy*bz*cx**2 - ay*bc2*bcdx*bcdz*bz*cx**2 + 
     -    2*ax*bc2*bcdy*bcdz*bz*cx**2 - abcz*bcdz**2*bcz*bz*cx**2 + 
     -    bc2*bcdy*bcdz*bx*bz*cx**2 - bc2*bcdx*bcdy*bz**2*cx**2 - 
     -    abcz*ay*bc2*bcdz*cdx*cx**2 + abcz*bc2*bcdz*by*cdx*cx**2 + 
     -    ay*bc2*bcdy**2*cx**3 + az*bc2*bcdy*bcdz*cx**3 - 
     -    bc2*bcdy**2*by*cx**3 - bc2*bcdy*bcdz*bz*cx**3 - 
     -    abcz*az*bc2*bcdy*bcdz*cy - abcz*ay*bc2*bcdz**2*cy + 
     -    ay**2*bc2*bcdx**2*bx*cy - ax**2*bc2*bcdy**2*bx*cy + 
     -    abcz*ay*bcdx*bcdz*bcz*bx*cy + abcz*ax*bcdy*bcdz*bcz*bx*cy - 
     -    2*ay*bc2*bcdx*bcdy*bx**2*cy + ax*bc2*bcdy**2*bx**2*cy - 
     -    az*bc2*bcdx*bcdz*bx**2*cy - abcz*bcdy*bcdz*bcz*bx**2*cy - 
     -    2*ax*ay*bc2*bcdx**2*by*cy + 2*ax**2*bc2*bcdx*bcdy*by*cy - 
     -    2*ay*az*bc2*bcdx*bcdz*by*cy + 2*ax*az*bc2*bcdy*bcdz*by*cy - 
     -    abcz*bc2*bcdz**2*by*cy - 2*abcz*ax*bcdx*bcdz*bcz*by*cy - 
     -    2*abcz*az*bcdz**2*bcz*by*cy + 2*ay*bc2*bcdx**2*bx*by*cy + 
     -    abcz*bcdx*bcdz*bcz*bx*by*cy - ax*bc2*bcdx**2*by**2*cy - 
     -    az*bc2*bcdx*bcdz*by**2*cy + ay*az*bc2*bcdx*bcdy*bz*cy - 
     -    ax*az*bc2*bcdy**2*bz*cy + ay**2*bc2*bcdx*bcdz*bz*cy + 
     -    2*abcz*bc2*bcdy*bcdz*bz*cy - ax*ay*bc2*bcdy*bcdz*bz*cy + 
     -    abcz*az*bcdy*bcdz*bcz*bz*cy + abcz*ay*bcdz**2*bcz*bz*cy + 
     -    az*bc2*bcdx**2*bx*bz*cy + ax*bc2*bcdx*bcdz*bx*bz*cy + 
     -    az*bc2*bcdx*bcdy*by*bz*cy + 2*ay*bc2*bcdx*bcdz*by*bz*cy - 
     -    ax*bc2*bcdy*bcdz*by*bz*cy + abcz*bcdz**2*bcz*by*bz*cy - 
     -    ax*bc2*bcdx**2*bz**2*cy - 2*ay*bc2*bcdx*bcdy*bz**2*cy + 
     -    ax*bc2*bcdy**2*bz**2*cy - abcz*bcdy*bcdz*bcz*bz**2*cy - 
     -    abcz*ax*bc2*bcdz*bx*cdx*cy + abcz*bc2*bcdz*bx**2*cdx*cy - 
     -    abcz*az*bc2*bcdz*bz*cdx*cy + abcz*bc2*bcdz*bz**2*cdx*cy + 
     -    abcz*ay*bc2*bcdz*bx*cdy*cy - 2*abcz*ax*bc2*bcdz*by*cdy*cy + 
     -    abcz*bc2*bcdz*bx*by*cdy*cy - ay**2*bc2*bcdx**2*cx*cy + 
     -    ax**2*bc2*bcdy**2*cx*cy - abcz*ay*bcdx*bcdz*bcz*cx*cy - 
     -    abcz*ax*bcdy*bcdz*bcz*cx*cy + 4*ay*bc2*bcdx*bcdy*bx*cx*cy + 
     -    2*az*bc2*bcdx*bcdz*bx*cx*cy + abcz*bcdy*bcdz*bcz*bx*cx*cy - 
     -    bc2*bcdy**2*bx**2*cx*cy - 4*ax*bc2*bcdx*bcdy*by*cx*cy - 
     -    2*az*bc2*bcdy*bcdz*by*cx*cy + abcz*bcdx*bcdz*bcz*by*cx*cy + 
     -    bc2*bcdx**2*by**2*cx*cy - az*bc2*bcdx**2*bz*cx*cy + 
     -    az*bc2*bcdy**2*bz*cx*cy - ax*bc2*bcdx*bcdz*bz*cx*cy + 
     -    ay*bc2*bcdy*bcdz*bz*cx*cy - bc2*bcdx*bcdz*bx*bz*cx*cy + 
     -    bc2*bcdy*bcdz*by*bz*cx*cy + bc2*bcdx**2*bz**2*cx*cy - 
     -    bc2*bcdy**2*bz**2*cx*cy + abcz*ax*bc2*bcdz*cdx*cx*cy - 
     -    abcz*bc2*bcdz*bx*cdx*cx*cy - abcz*ay*bc2*bcdz*cdy*cx*cy + 
     -    abcz*bc2*bcdz*by*cdy*cx*cy - 2*ay*bc2*bcdx*bcdy*cx**2*cy - 
     -    ax*bc2*bcdy**2*cx**2*cy - az*bc2*bcdx*bcdz*cx**2*cy + 
     -    bc2*bcdy**2*bx*cx**2*cy + 2*bc2*bcdx*bcdy*by*cx**2*cy + 
     -    bc2*bcdx*bcdz*bz*cx**2*cy + ax*ay*bc2*bcdx**2*cy**2 - 
     -    ax**2*bc2*bcdx*bcdy*cy**2 + ay*az*bc2*bcdx*bcdz*cy**2 - 
     -    ax*az*bc2*bcdy*bcdz*cy**2 + abcz*bc2*bcdz**2*cy**2 + 
     -    abcz*ax*bcdx*bcdz*bcz*cy**2 + abcz*az*bcdz**2*bcz*cy**2 - 
     -    2*ay*bc2*bcdx**2*bx*cy**2 - abcz*bcdx*bcdz*bcz*bx*cy**2 + 
     -    bc2*bcdx*bcdy*bx**2*cy**2 + 2*ax*bc2*bcdx**2*by*cy**2 + 
     -    2*az*bc2*bcdx*bcdz*by*cy**2 - bc2*bcdx**2*bx*by*cy**2 - 
     -    az*bc2*bcdx*bcdy*bz*cy**2 - 2*ay*bc2*bcdx*bcdz*bz*cy**2 + 
     -    ax*bc2*bcdy*bcdz*bz*cy**2 - abcz*bcdz**2*bcz*bz*cy**2 - 
     -    bc2*bcdx*bcdz*by*bz*cy**2 + bc2*bcdx*bcdy*bz**2*cy**2 + 
     -    abcz*ax*bc2*bcdz*cdy*cy**2 - abcz*bc2*bcdz*bx*cdy*cy**2 + 
     -    ay*bc2*bcdx**2*cx*cy**2 + 2*ax*bc2*bcdx*bcdy*cx*cy**2 + 
     -    az*bc2*bcdy*bcdz*cx*cy**2 - 2*bc2*bcdx*bcdy*bx*cx*cy**2 - 
     -    bc2*bcdx**2*by*cx*cy**2 - bc2*bcdy*bcdz*bz*cx*cy**2 - 
     -    ax*bc2*bcdx**2*cy**3 - az*bc2*bcdx*bcdz*cy**3 + 
     -    bc2*bcdx**2*bx*cy**3 + bc2*bcdx*bcdz*bz*cy**3 + 
     -    2*abcz*ax*bc2*bcdx*bcdz*cz + 2*abcz*ay*bc2*bcdy*bcdz*cz + 
     -    ay*az*bc2*bcdx**2*bx*cz - ax*az*bc2*bcdx*bcdy*bx*cz - 
     -    abcz*bc2*bcdx*bcdz*bx*cz + ax*ay*bc2*bcdx*bcdz*bx*cz - 
     -    ax**2*bc2*bcdy*bcdz*bx*cz + abcz*az*bcdx*bcdz*bcz*bx*cz + 
     -    abcz*ax*bcdz**2*bcz*bx*cz - ay*bc2*bcdx*bcdz*bx**2*cz + 
     -    ax*bc2*bcdy*bcdz*bx**2*cz - abcz*bcdz**2*bcz*bx**2*cz + 
     -    ay*az*bc2*bcdx*bcdy*by*cz - ax*az*bc2*bcdy**2*by*cz + 
     -    ay**2*bc2*bcdx*bcdz*by*cz - abcz*bc2*bcdy*bcdz*by*cz - 
     -    ax*ay*bc2*bcdy*bcdz*by*cz + abcz*az*bcdy*bcdz*bcz*by*cz + 
     -    abcz*ay*bcdz**2*bcz*by*cz - ay*bc2*bcdx*bcdz*by**2*cz + 
     -    ax*bc2*bcdy*bcdz*by**2*cz - abcz*bcdz**2*bcz*by**2*cz - 
     -    2*ax*ay*bc2*bcdx**2*bz*cz + 2*ax**2*bc2*bcdx*bcdy*bz*cz - 
     -    2*ay**2*bc2*bcdx*bcdy*bz*cz + 2*ax*ay*bc2*bcdy**2*bz*cz - 
     -    2*abcz*ax*bcdx*bcdz*bcz*bz*cz - 
     -    2*abcz*ay*bcdy*bcdz*bcz*bz*cz + ay*bc2*bcdx**2*bx*bz*cz - 
     -    ax*bc2*bcdx*bcdy*bx*bz*cz + abcz*bcdx*bcdz*bcz*bx*bz*cz + 
     -    ay*bc2*bcdx*bcdy*by*bz*cz - ax*bc2*bcdy**2*by*bz*cz + 
     -    abcz*bcdy*bcdz*bcz*by*bz*cz - abcz*az*bc2*bcdz*by*cdx*cz + 
     -    2*abcz*ay*bc2*bcdz*bz*cdx*cz - abcz*bc2*bcdz*by*bz*cdx*cz + 
     -    abcz*az*bc2*bcdz*bx*cdy*cz - 2*abcz*ax*bc2*bcdz*bz*cdy*cz + 
     -    abcz*bc2*bcdz*bx*bz*cdy*cz - ay*az*bc2*bcdx**2*cx*cz + 
     -    ax*az*bc2*bcdx*bcdy*cx*cz - abcz*bc2*bcdx*bcdz*cx*cz - 
     -    ax*ay*bc2*bcdx*bcdz*cx*cz + ax**2*bc2*bcdy*bcdz*cx*cz - 
     -    abcz*az*bcdx*bcdz*bcz*cx*cz - abcz*ax*bcdz**2*bcz*cx*cz + 
     -    az*bc2*bcdx*bcdy*bx*cx*cz + ay*bc2*bcdx*bcdz*bx*cx*cz + 
     -    abcz*bcdz**2*bcz*bx*cx*cz - bc2*bcdy*bcdz*bx**2*cx*cz + 
     -    az*bc2*bcdy**2*by*cx*cz + ay*bc2*bcdy*bcdz*by*cx*cz - 
     -    bc2*bcdy*bcdz*by**2*cx*cz + ay*bc2*bcdx**2*bz*cx*cz - 
     -    3*ax*bc2*bcdx*bcdy*bz*cx*cz - 2*ay*bc2*bcdy**2*bz*cx*cz + 
     -    abcz*bcdx*bcdz*bcz*bz*cx*cz + bc2*bcdx*bcdy*bx*bz*cx*cz + 
     -    bc2*bcdy**2*by*bz*cx*cz - abcz*az*bc2*bcdz*cdy*cx*cz + 
     -    abcz*bc2*bcdz*bz*cdy*cx*cz - az*bc2*bcdx*bcdy*cx**2*cz - 
     -    ax*bc2*bcdy*bcdz*cx**2*cz + bc2*bcdy*bcdz*bx*cx**2*cz + 
     -    bc2*bcdx*bcdy*bz*cx**2*cz - ay*az*bc2*bcdx*bcdy*cy*cz + 
     -    ax*az*bc2*bcdy**2*cy*cz - ay**2*bc2*bcdx*bcdz*cy*cz - 
     -    abcz*bc2*bcdy*bcdz*cy*cz + ax*ay*bc2*bcdy*bcdz*cy*cz - 
     -    abcz*az*bcdy*bcdz*bcz*cy*cz - abcz*ay*bcdz**2*bcz*cy*cz - 
     -    az*bc2*bcdx**2*bx*cy*cz - ax*bc2*bcdx*bcdz*bx*cy*cz + 
     -    bc2*bcdx*bcdz*bx**2*cy*cz - az*bc2*bcdx*bcdy*by*cy*cz - 
     -    ax*bc2*bcdy*bcdz*by*cy*cz + abcz*bcdz**2*bcz*by*cy*cz + 
     -    bc2*bcdx*bcdz*by**2*cy*cz + 2*ax*bc2*bcdx**2*bz*cy*cz + 
     -    3*ay*bc2*bcdx*bcdy*bz*cy*cz - ax*bc2*bcdy**2*bz*cy*cz + 
     -    abcz*bcdy*bcdz*bcz*bz*cy*cz - bc2*bcdx**2*bx*bz*cy*cz - 
     -    bc2*bcdx*bcdy*by*bz*cy*cz + abcz*az*bc2*bcdz*cdx*cy*cz - 
     -    abcz*bc2*bcdz*bz*cdx*cy*cz + az*bc2*bcdx**2*cx*cy*cz - 
     -    az*bc2*bcdy**2*cx*cy*cz + ax*bc2*bcdx*bcdz*cx*cy*cz - 
     -    ay*bc2*bcdy*bcdz*cx*cy*cz - bc2*bcdx*bcdz*bx*cx*cy*cz + 
     -    bc2*bcdy*bcdz*by*cx*cy*cz - bc2*bcdx**2*bz*cx*cy*cz + 
     -    bc2*bcdy**2*bz*cx*cy*cz + az*bc2*bcdx*bcdy*cy**2*cz + 
     -    ay*bc2*bcdx*bcdz*cy**2*cz - bc2*bcdx*bcdz*by*cy**2*cz - 
     -    bc2*bcdx*bcdy*bz*cy**2*cz + ax*ay*bc2*bcdx**2*cz**2 - 
     -    ax**2*bc2*bcdx*bcdy*cz**2 + ay**2*bc2*bcdx*bcdy*cz**2 - 
     -    ax*ay*bc2*bcdy**2*cz**2 + abcz*ax*bcdx*bcdz*bcz*cz**2 + 
     -    abcz*ay*bcdy*bcdz*bcz*cz**2 - ay*bc2*bcdx**2*bx*cz**2 + 
     -    ax*bc2*bcdx*bcdy*bx*cz**2 - abcz*bcdx*bcdz*bcz*bx*cz**2 - 
     -    ay*bc2*bcdx*bcdy*by*cz**2 + ax*bc2*bcdy**2*by*cz**2 - 
     -    abcz*bcdy*bcdz*bcz*by*cz**2 - abcz*ay*bc2*bcdz*cdx*cz**2 + 
     -    abcz*bc2*bcdz*by*cdx*cz**2 + abcz*ax*bc2*bcdz*cdy*cz**2 - 
     -    abcz*bc2*bcdz*bx*cdy*cz**2 + ax*bc2*bcdx*bcdy*cx*cz**2 + 
     -    ay*bc2*bcdy**2*cx*cz**2 - bc2*bcdx*bcdy*bx*cx*cz**2 - 
     -    bc2*bcdy**2*by*cx*cz**2 - ax*bc2*bcdx**2*cy*cz**2 - 
     -    ay*bc2*bcdx*bcdy*cy*cz**2 + bc2*bcdx**2*bx*cy*cz**2 + 
     -    bc2*bcdx*bcdy*by*cy*cz**2 + 
     -    abcy*(ay*bcdy**2*bcz*bx**2 + ay*bc2*bcdy*bcdz*by - 
     -       ay*bcdx*bcdy*bcz*bx*by - 2*ay*bc2*bcdy**2*bz - 
     -       ay*bcdy*bcdz*bcz*by*bz + ay*bcdy**2*bcz*bz**2 - 
     -       ay*bc2*bcdx*bx*by*cdx - ay*bc2*bcdz*by*bz*cdx - 
     -       ay*bc2*bcdy*bx*by*cdy + 
     -       az*(bcdy*bcz*(bcdz*(bcx**2 + bcy**2) - 
     -             (abx*bcdx + bcdy*bcy)*bcz) + 
     -          bc2*(bcdy**2*bcy + bcdz*(bcx**2 + bcy**2)*cdx + 
     -             abx*bcdx*(bcdy - bcz*cdx) - abx*bcdy*bcz*cdy)) + 
     -       ax*(bcdy*bcz*(-(abx*(bcdy*bcy + bcdz*bcz)) + 
     -             bcdx*(bcy**2 + bcz**2)) + 
     -          bc2*(abx*bcdz*(bcdy - bcz*cdx) + 
     -             bcdx*(-2*bcdy*bcz + (bcy**2 + bcz**2)*cdx) + 
     -             bcdy*(bcy**2 + bcz**2)*cdy)) - bc2*bcdy*bcdz*bx*cx - 
     -       2*ay*bcdy**2*bcz*bx*cx + ay*bcdx*bcdy*bcz*by*cx + 
     -       bcdy**2*bcz*bx*by*cx - bcdx*bcdy*bcz*by**2*cx + 
     -       2*bc2*bcdx*bcdy*bz*cx + bcdy*bcdz*bcz*bx*bz*cx - 
     -       bcdx*bcdy*bcz*bz**2*cx + ay*bc2*bcdx*by*cdx*cx - 
     -       bc2*bcdx*by**2*cdx*cx + bc2*bcdz*bx*bz*cdx*cx - 
     -       bc2*bcdx*bz**2*cdx*cx + ay*bc2*bcdy*by*cdy*cx - 
     -       bc2*bcdy*by**2*cdy*cx - bc2*bcdy*bz**2*cdy*cx + 
     -       bc2*bcdy*bcdz*cx**2 + ay*bcdy**2*bcz*cx**2 - 
     -       bcdy**2*bcz*by*cx**2 - bcdy*bcdz*bcz*bz*cx**2 - 
     -       bc2*bcdz*bz*cdx*cx**2 - ay*bc2*bcdy*bcdz*cy + 
     -       ay*bcdx*bcdy*bcz*bx*cy - bcdy**2*bcz*bx**2*cy - 
     -       bc2*bcdy*bcdz*by*cy + bcdx*bcdy*bcz*bx*by*cy + 
     -       2*bc2*bcdy**2*bz*cy + ay*bcdy*bcdz*bcz*bz*cy + 
     -       bcdy*bcdz*bcz*by*bz*cy - bcdy**2*bcz*bz**2*cy + 
     -       ay*bc2*bcdx*bx*cdx*cy + bc2*bcdx*bx*by*cdx*cy + 
     -       ay*bc2*bcdz*bz*cdx*cy + bc2*bcdz*by*bz*cdx*cy + 
     -       ay*bc2*bcdy*bx*cdy*cy + bc2*bcdy*bx*by*cdy*cy - 
     -       ay*bcdx*bcdy*bcz*cx*cy + bcdy**2*bcz*bx*cx*cy + 
     -       bcdx*bcdy*bcz*by*cx*cy - ay*bc2*bcdx*cdx*cx*cy + 
     -       bc2*bcdx*by*cdx*cx*cy - ay*bc2*bcdy*cdy*cx*cy + 
     -       bc2*bcdy*by*cdy*cx*cy + bc2*bcdy*bcdz*cy**2 - 
     -       bcdx*bcdy*bcz*bx*cy**2 - bcdy*bcdz*bcz*bz*cy**2 - 
     -       bc2*bcdx*bx*cdx*cy**2 - bc2*bcdz*bz*cdx*cy**2 - 
     -       bc2*bcdy*bx*cdy*cy**2 + 2*ay*bc2*bcdy**2*cz - 
     -       bc2*bcdx*bcdy*bx*cz - bcdy*bcdz*bcz*bx**2*cz - 
     -       bc2*bcdy**2*by*cz + ay*bcdy*bcdz*bcz*by*cz - 
     -       bcdy*bcdz*bcz*by**2*cz - 2*ay*bcdy**2*bcz*bz*cz + 
     -       bcdx*bcdy*bcz*bx*bz*cz + bcdy**2*bcz*by*bz*cz - 
     -       bc2*bcdz*bx**2*cdx*cz + ay*bc2*bcdz*by*cdx*cz - 
     -       bc2*bcdz*by**2*cdx*cz + bc2*bcdx*bx*bz*cdx*cz + 
     -       bc2*bcdy*bx*bz*cdy*cz - bc2*bcdx*bcdy*cx*cz + 
     -       bcdy*bcdz*bcz*bx*cx*cz + bcdx*bcdy*bcz*bz*cx*cz + 
     -       bc2*bcdz*bx*cdx*cx*cz + bc2*bcdx*bz*cdx*cx*cz + 
     -       bc2*bcdy*bz*cdy*cx*cz - bc2*bcdy**2*cy*cz - 
     -       ay*bcdy*bcdz*bcz*cy*cz + bcdy*bcdz*bcz*by*cy*cz + 
     -       bcdy**2*bcz*bz*cy*cz - ay*bc2*bcdz*cdx*cy*cz + 
     -       bc2*bcdz*by*cdx*cy*cz + ay*bcdy**2*bcz*cz**2 - 
     -       bcdx*bcdy*bcz*bx*cz**2 - bcdy**2*bcz*by*cz**2 - 
     -       bc2*bcdx*bx*cdx*cz**2 - bc2*bcdy*bx*cdy*cz**2) + 
     -    abcx*(ay*bcdx*bcdy*bcz*bx**2 + ay*bc2*bcdx*bcdz*by - 
     -       ay*bcdx**2*bcz*bx*by - 2*ay*bc2*bcdx*bcdy*bz - 
     -       ay*bcdx*bcdz*bcz*by*bz + ay*bcdx*bcdy*bcz*bz**2 - 
     -       ay*bc2*bcdx*bx**2*cdx - ay*bc2*bcdx*bz**2*cdx - 
     -       ay*bc2*bcdy*bx**2*cdy + ay*bc2*bcdz*by*bz*cdy - 
     -       ay*bc2*bcdy*bz**2*cdy - 
     -       ax*(-(bcdx**2*bcz*(-2*bc2 + bcy**2 + bcz**2)) - 
     -          abx*(-(bcdx*bcz*(bcdy*bcy + bcdz*bcz)) + 
     -             bc2*(bcdx*bcdz + bcdx*bcy*cdx + bcdy*bcy*cdy + 
     -                bcdz*bcz*cdy))) - bc2*bcdx*bcdz*bx*cx - 
     -       2*ay*bcdx*bcdy*bcz*bx*cx + ay*bcdx**2*bcz*by*cx + 
     -       bcdx*bcdy*bcz*bx*by*cx - bcdx**2*bcz*by**2*cx + 
     -       2*bc2*bcdx**2*bz*cx + bcdx*bcdz*bcz*bx*bz*cx - 
     -       bcdx**2*bcz*bz**2*cx + 2*ay*bc2*bcdx*bx*cdx*cx - 
     -       bc2*bcdx*bx*by*cdx*cx + 2*ay*bc2*bcdy*bx*cdy*cx - 
     -       bc2*bcdy*bx*by*cdy*cx - bc2*bcdz*bx*bz*cdy*cx + 
     -       bc2*bcdx*bcdz*cx**2 + ay*bcdx*bcdy*bcz*cx**2 - 
     -       bcdx*bcdy*bcz*by*cx**2 - bcdx*bcdz*bcz*bz*cx**2 - 
     -       ay*bc2*bcdx*cdx*cx**2 + bc2*bcdx*by*cdx*cx**2 - 
     -       ay*bc2*bcdy*cdy*cx**2 + bc2*bcdy*by*cdy*cx**2 + 
     -       bc2*bcdz*bz*cdy*cx**2 - ay*bc2*bcdx*bcdz*cy + 
     -       ay*bcdx**2*bcz*bx*cy - bcdx*bcdy*bcz*bx**2*cy - 
     -       bc2*bcdx*bcdz*by*cy + bcdx**2*bcz*bx*by*cy + 
     -       2*bc2*bcdx*bcdy*bz*cy + ay*bcdx*bcdz*bcz*bz*cy + 
     -       bcdx*bcdz*bcz*by*bz*cy - bcdx*bcdy*bcz*bz**2*cy + 
     -       bc2*bcdx*bx**2*cdx*cy + bc2*bcdx*bz**2*cdx*cy + 
     -       bc2*bcdy*bx**2*cdy*cy - ay*bc2*bcdz*bz*cdy*cy - 
     -       bc2*bcdz*by*bz*cdy*cy + bc2*bcdy*bz**2*cdy*cy - 
     -       ay*bcdx**2*bcz*cx*cy + bcdx*bcdy*bcz*bx*cx*cy + 
     -       bcdx**2*bcz*by*cx*cy - bc2*bcdx*bx*cdx*cx*cy - 
     -       bc2*bcdy*bx*cdy*cx*cy + bc2*bcdx*bcdz*cy**2 - 
     -       bcdx**2*bcz*bx*cy**2 - bcdx*bcdz*bcz*bz*cy**2 + 
     -       bc2*bcdz*bz*cdy*cy**2 + 
     -       az*(bc2*(abx*bcdx**2 + bcdx*bcy*(bcdy + bcz*cdx) - 
     -             (bcdz*(bcx**2 + bcy**2) - bcdy*bcy*bcz)*cdy) + 
     -          bcdx*bcz*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(bcdx*(-bx + cx) + bcdy*(-by + cy)))) + 
     -       2*ay*bc2*bcdx*bcdy*cz - bc2*bcdx**2*bx*cz - 
     -       bcdx*bcdz*bcz*bx**2*cz - bc2*bcdx*bcdy*by*cz + 
     -       ay*bcdx*bcdz*bcz*by*cz - bcdx*bcdz*bcz*by**2*cz - 
     -       2*ay*bcdx*bcdy*bcz*bz*cz + bcdx**2*bcz*bx*bz*cz + 
     -       bcdx*bcdy*bcz*by*bz*cz + 2*ay*bc2*bcdx*bz*cdx*cz - 
     -       bc2*bcdx*by*bz*cdx*cz + bc2*bcdz*bx**2*cdy*cz - 
     -       ay*bc2*bcdz*by*cdy*cz + bc2*bcdz*by**2*cdy*cz + 
     -       2*ay*bc2*bcdy*bz*cdy*cz - bc2*bcdy*by*bz*cdy*cz - 
     -       bc2*bcdx**2*cx*cz + bcdx*bcdz*bcz*bx*cx*cz + 
     -       bcdx**2*bcz*bz*cx*cz - bc2*bcdz*bx*cdy*cx*cz - 
     -       bc2*bcdx*bcdy*cy*cz - ay*bcdx*bcdz*bcz*cy*cz + 
     -       bcdx*bcdz*bcz*by*cy*cz + bcdx*bcdy*bcz*bz*cy*cz - 
     -       bc2*bcdx*bz*cdx*cy*cz + ay*bc2*bcdz*cdy*cy*cz - 
     -       bc2*bcdz*by*cdy*cy*cz - bc2*bcdy*bz*cdy*cy*cz + 
     -       ay*bcdx*bcdy*bcz*cz**2 - bcdx**2*bcz*bx*cz**2 - 
     -       bcdx*bcdy*bcz*by*cz**2 - ay*bc2*bcdx*cdx*cz**2 + 
     -       bc2*bcdx*by*cdx*cz**2 - ay*bc2*bcdy*cdy*cz**2 + 
     -       bc2*bcdy*by*cdy*cz**2))/
     -  (bc2**1.5*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**2/
     -       bc2))

       dcx= (2*abcz*ay*bc2*bcdy*bcdz*bx + 2*abcz*az*bc2*bcdz**2*bx - 
     -    abz*ay*bc2*bcdy**2*bx**2 + aby*ay*bc2*bcdy*bcdz*bx**2 - 
     -    abz*az*bc2*bcdy*bcdz*bx**2 + aby*az*bc2*bcdz**2*bx**2 - 
     -    abcz*ay*bcdy*bcdz*bcx*bx**2 - abcz*az*bcdz**2*bcx*bx**2 - 
     -    abcz*ay*bc2*bcdx*bcdz*by - abcz*ax*bc2*bcdy*bcdz*by + 
     -    abz*ay*bc2*bcdx*bcdy*bx*by + abz*ax*bc2*bcdy**2*bx*by - 
     -    aby*ay*bc2*bcdx*bcdz*bx*by - abcz*bc2*bcdy*bcdz*bx*by - 
     -    aby*ax*bc2*bcdy*bcdz*bx*by + abcz*ay*bcdx*bcdz*bcx*bx*by + 
     -    abcz*ax*bcdy*bcdz*bcx*bx*by - abcz*ay*bc2*bcdy*bx**2*by - 
     -    abz*ax*bc2*bcdx*bcdy*by**2 + abcz*bc2*bcdx*bcdz*by**2 + 
     -    aby*ax*bc2*bcdx*bcdz*by**2 - abz*az*bc2*bcdy*bcdz*by**2 + 
     -    aby*az*bc2*bcdz**2*by**2 - abcz*ax*bcdx*bcdz*bcx*by**2 - 
     -    abcz*az*bcdz**2*bcx*by**2 + abcz*ay*bc2*bcdx*bx*by**2 + 
     -    abcz*ax*bc2*bcdy*bx*by**2 - abcz*ax*bc2*bcdx*by**3 - 
     -    abcz*az*bc2*bcdx*bcdz*bz - abcz*ax*bc2*bcdz**2*bz + 
     -    abz*az*bc2*bcdx*bcdy*bx*bz - aby*az*bc2*bcdx*bcdz*bx*bz + 
     -    abz*ax*bc2*bcdy*bcdz*bx*bz - abcz*bc2*bcdz**2*bx*bz - 
     -    aby*ax*bc2*bcdz**2*bx*bz + abcz*az*bcdx*bcdz*bcx*bx*bz + 
     -    abcz*ax*bcdz**2*bcx*bx*bz - abcz*ay*bc2*bcdz*bx**2*bz + 
     -    abz*az*bc2*bcdy**2*by*bz + abz*ay*bc2*bcdy*bcdz*by*bz - 
     -    aby*az*bc2*bcdy*bcdz*by*bz - aby*ay*bc2*bcdz**2*by*bz + 
     -    abcz*az*bcdy*bcdz*bcx*by*bz + abcz*ay*bcdz**2*bcx*by*bz + 
     -    abcz*az*bc2*bcdx*bx*by*bz + abcz*ax*bc2*bcdz*bx*by*bz + 
     -    abcz*az*bc2*bcdy*by**2*bz - abz*ax*bc2*bcdx*bcdy*bz**2 - 
     -    abz*ay*bc2*bcdy**2*bz**2 + abcz*bc2*bcdx*bcdz*bz**2 + 
     -    aby*ax*bc2*bcdx*bcdz*bz**2 + aby*ay*bc2*bcdy*bcdz*bz**2 - 
     -    abcz*ax*bcdx*bcdz*bcx*bz**2 - abcz*ay*bcdy*bcdz*bcx*bz**2 - 
     -    abcz*ax*bc2*bcdx*by*bz**2 - abcz*ay*bc2*bcdy*by*bz**2 + 
     -    abcz*az*bc2*bcdz*by*bz**2 - abcz*ay*bc2*bcdz*bz**3 - 
     -    2*abcz*ay*bc2*bcdy*bcdz*cx - 2*abcz*az*bc2*bcdz**2*cx + 
     -    2*abz*ay*bc2*bcdy**2*bx*cx - 2*aby*ay*bc2*bcdy*bcdz*bx*cx + 
     -    2*abz*az*bc2*bcdy*bcdz*bx*cx - 2*aby*az*bc2*bcdz**2*bx*cx + 
     -    2*abcz*ay*bcdy*bcdz*bcx*bx*cx + 2*abcz*az*bcdz**2*bcx*bx*cx - 
     -    abz*ay*bc2*bcdx*bcdy*by*cx - abz*ax*bc2*bcdy**2*by*cx + 
     -    aby*ay*bc2*bcdx*bcdz*by*cx + 2*abcz*bc2*bcdy*bcdz*by*cx + 
     -    aby*ax*bc2*bcdy*bcdz*by*cx - abcz*ay*bcdx*bcdz*bcx*by*cx - 
     -    abcz*ax*bcdy*bcdz*bcx*by*cx + 2*abcz*ay*bc2*bcdy*bx*by*cx - 
     -    abz*bc2*bcdy**2*bx*by*cx + aby*bc2*bcdy*bcdz*bx*by*cx - 
     -    abcz*bcdy*bcdz*bcx*bx*by*cx - abcz*ay*bc2*bcdx*by**2*cx - 
     -    abcz*ax*bc2*bcdy*by**2*cx + abz*bc2*bcdx*bcdy*by**2*cx - 
     -    aby*bc2*bcdx*bcdz*by**2*cx + abcz*bcdx*bcdz*bcx*by**2*cx - 
     -    abcz*bc2*bcdy*bx*by**2*cx + abcz*bc2*bcdx*by**3*cx - 
     -    abz*az*bc2*bcdx*bcdy*bz*cx + aby*az*bc2*bcdx*bcdz*bz*cx - 
     -    abz*ax*bc2*bcdy*bcdz*bz*cx + 2*abcz*bc2*bcdz**2*bz*cx + 
     -    aby*ax*bc2*bcdz**2*bz*cx - abcz*az*bcdx*bcdz*bcx*bz*cx - 
     -    abcz*ax*bcdz**2*bcx*bz*cx + 2*abcz*ay*bc2*bcdz*bx*bz*cx - 
     -    abz*bc2*bcdy*bcdz*bx*bz*cx + aby*bc2*bcdz**2*bx*bz*cx - 
     -    abcz*bcdz**2*bcx*bx*bz*cx - abcz*az*bc2*bcdx*by*bz*cx - 
     -    abcz*ax*bc2*bcdz*by*bz*cx - abcz*bc2*bcdz*bx*by*bz*cx + 
     -    abz*bc2*bcdx*bcdy*bz**2*cx - aby*bc2*bcdx*bcdz*bz**2*cx + 
     -    abcz*bcdx*bcdz*bcx*bz**2*cx + abcz*bc2*bcdx*by*bz**2*cx - 
     -    abz*ay*bc2*bcdy**2*cx**2 + aby*ay*bc2*bcdy*bcdz*cx**2 - 
     -    abz*az*bc2*bcdy*bcdz*cx**2 + aby*az*bc2*bcdz**2*cx**2 - 
     -    abcz*ay*bcdy*bcdz*bcx*cx**2 - abcz*az*bcdz**2*bcx*cx**2 - 
     -    abcz*ay*bc2*bcdy*by*cx**2 + abz*bc2*bcdy**2*by*cx**2 - 
     -    aby*bc2*bcdy*bcdz*by*cx**2 + abcz*bcdy*bcdz*bcx*by*cx**2 + 
     -    abcz*bc2*bcdy*by**2*cx**2 - abcz*ay*bc2*bcdz*bz*cx**2 + 
     -    abz*bc2*bcdy*bcdz*bz*cx**2 - aby*bc2*bcdz**2*bz*cx**2 + 
     -    abcz*bcdz**2*bcx*bz*cx**2 + abcz*bc2*bcdz*by*bz*cx**2 + 
     -    abcz*ay*bc2*bcdx*bcdz*cy + abcz*ax*bc2*bcdy*bcdz*cy - 
     -    abz*ay*bc2*bcdx*bcdy*bx*cy - abz*ax*bc2*bcdy**2*bx*cy + 
     -    aby*ay*bc2*bcdx*bcdz*bx*cy - abcz*bc2*bcdy*bcdz*bx*cy + 
     -    aby*ax*bc2*bcdy*bcdz*bx*cy - abcz*ay*bcdx*bcdz*bcx*bx*cy - 
     -    abcz*ax*bcdy*bcdz*bcx*bx*cy + abz*bc2*bcdy**2*bx**2*cy - 
     -    aby*bc2*bcdy*bcdz*bx**2*cy + abcz*bcdy*bcdz*bcx*bx**2*cy + 
     -    2*abz*ax*bc2*bcdx*bcdy*by*cy - abcz*bc2*bcdx*bcdz*by*cy - 
     -    2*aby*ax*bc2*bcdx*bcdz*by*cy + 2*abz*az*bc2*bcdy*bcdz*by*cy - 
     -    2*aby*az*bc2*bcdz**2*by*cy + 2*abcz*ax*bcdx*bcdz*bcx*by*cy + 
     -    2*abcz*az*bcdz**2*bcx*by*cy - abcz*ay*bc2*bcdx*bx*by*cy - 
     -    abcz*ax*bc2*bcdy*bx*by*cy - abz*bc2*bcdx*bcdy*bx*by*cy + 
     -    aby*bc2*bcdx*bcdz*bx*by*cy - abcz*bcdx*bcdz*bcx*bx*by*cy + 
     -    abcz*bc2*bcdy*bx**2*by*cy + 2*abcz*ax*bc2*bcdx*by**2*cy - 
     -    abcz*bc2*bcdx*bx*by**2*cy - abz*az*bc2*bcdy**2*bz*cy - 
     -    abz*ay*bc2*bcdy*bcdz*bz*cy + aby*az*bc2*bcdy*bcdz*bz*cy + 
     -    aby*ay*bc2*bcdz**2*bz*cy - abcz*az*bcdy*bcdz*bcx*bz*cy - 
     -    abcz*ay*bcdz**2*bcx*bz*cy - abcz*ax*bc2*bcdz*bx*bz*cy + 
     -    abcz*bc2*bcdz*bx**2*bz*cy - abcz*az*bc2*bcdy*by*bz*cy - 
     -    abz*bc2*bcdy*bcdz*by*bz*cy + aby*bc2*bcdz**2*by*bz*cy - 
     -    abcz*bcdz**2*bcx*by*bz*cy + abz*bc2*bcdy**2*bz**2*cy - 
     -    abcz*az*bc2*bcdz*bz**2*cy - aby*bc2*bcdy*bcdz*bz**2*cy + 
     -    abcz*bcdy*bcdz*bcx*bz**2*cy + abcz*bc2*bcdy*by*bz**2*cy + 
     -    abcz*bc2*bcdz*bz**3*cy + abz*ay*bc2*bcdx*bcdy*cx*cy + 
     -    abz*ax*bc2*bcdy**2*cx*cy - aby*ay*bc2*bcdx*bcdz*cx*cy - 
     -    aby*ax*bc2*bcdy*bcdz*cx*cy + abcz*ay*bcdx*bcdz*bcx*cx*cy + 
     -    abcz*ax*bcdy*bcdz*bcx*cx*cy - abz*bc2*bcdy**2*bx*cx*cy + 
     -    aby*bc2*bcdy*bcdz*bx*cx*cy - abcz*bcdy*bcdz*bcx*bx*cx*cy + 
     -    abcz*ay*bc2*bcdx*by*cx*cy + abcz*ax*bc2*bcdy*by*cx*cy - 
     -    abz*bc2*bcdx*bcdy*by*cx*cy + aby*bc2*bcdx*bcdz*by*cx*cy - 
     -    abcz*bcdx*bcdz*bcx*by*cx*cy - abcz*bc2*bcdy*bx*by*cx*cy - 
     -    abcz*bc2*bcdx*by**2*cx*cy + abcz*ax*bc2*bcdz*bz*cx*cy - 
     -    abcz*bc2*bcdz*bx*bz*cx*cy - abz*ax*bc2*bcdx*bcdy*cy**2 + 
     -    aby*ax*bc2*bcdx*bcdz*cy**2 - abz*az*bc2*bcdy*bcdz*cy**2 + 
     -    aby*az*bc2*bcdz**2*cy**2 - abcz*ax*bcdx*bcdz*bcx*cy**2 - 
     -    abcz*az*bcdz**2*bcx*cy**2 + abz*bc2*bcdx*bcdy*bx*cy**2 - 
     -    aby*bc2*bcdx*bcdz*bx*cy**2 + abcz*bcdx*bcdz*bcx*bx*cy**2 - 
     -    abcz*ax*bc2*bcdx*by*cy**2 + abcz*bc2*bcdx*bx*by*cy**2 + 
     -    abz*bc2*bcdy*bcdz*bz*cy**2 - aby*bc2*bcdz**2*bz*cy**2 + 
     -    abcz*bcdz**2*bcx*bz*cy**2 + abcz*az*bc2*bcdx*bcdz*cz + 
     -    abcz*ax*bc2*bcdz**2*cz - abz*az*bc2*bcdx*bcdy*bx*cz + 
     -    aby*az*bc2*bcdx*bcdz*bx*cz - abz*ax*bc2*bcdy*bcdz*bx*cz - 
     -    abcz*bc2*bcdz**2*bx*cz + aby*ax*bc2*bcdz**2*bx*cz - 
     -    abcz*az*bcdx*bcdz*bcx*bx*cz - abcz*ax*bcdz**2*bcx*bx*cz + 
     -    abz*bc2*bcdy*bcdz*bx**2*cz - aby*bc2*bcdz**2*bx**2*cz + 
     -    abcz*bcdz**2*bcx*bx**2*cz - abz*az*bc2*bcdy**2*by*cz - 
     -    abz*ay*bc2*bcdy*bcdz*by*cz + aby*az*bc2*bcdy*bcdz*by*cz + 
     -    aby*ay*bc2*bcdz**2*by*cz - abcz*az*bcdy*bcdz*bcx*by*cz - 
     -    abcz*ay*bcdz**2*bcx*by*cz - abcz*az*bc2*bcdx*bx*by*cz - 
     -    abcz*az*bc2*bcdy*by**2*cz + abz*bc2*bcdy*bcdz*by**2*cz - 
     -    aby*bc2*bcdz**2*by**2*cz + abcz*bcdz**2*bcx*by**2*cz + 
     -    2*abz*ax*bc2*bcdx*bcdy*bz*cz + 2*abz*ay*bc2*bcdy**2*bz*cz - 
     -    abcz*bc2*bcdx*bcdz*bz*cz - 2*aby*ax*bc2*bcdx*bcdz*bz*cz - 
     -    2*aby*ay*bc2*bcdy*bcdz*bz*cz + 
     -    2*abcz*ax*bcdx*bcdz*bcx*bz*cz + 
     -    2*abcz*ay*bcdy*bcdz*bcx*bz*cz - abz*bc2*bcdx*bcdy*bx*bz*cz + 
     -    aby*bc2*bcdx*bcdz*bx*bz*cz - abcz*bcdx*bcdz*bcx*bx*bz*cz + 
     -    2*abcz*ax*bc2*bcdx*by*bz*cz + 2*abcz*ay*bc2*bcdy*by*bz*cz - 
     -    abz*bc2*bcdy**2*by*bz*cz - abcz*az*bc2*bcdz*by*bz*cz + 
     -    aby*bc2*bcdy*bcdz*by*bz*cz - abcz*bcdy*bcdz*bcx*by*bz*cz - 
     -    abcz*bc2*bcdx*bx*by*bz*cz - abcz*bc2*bcdy*by**2*bz*cz + 
     -    2*abcz*ay*bc2*bcdz*bz**2*cz - abcz*bc2*bcdz*by*bz**2*cz + 
     -    abz*az*bc2*bcdx*bcdy*cx*cz - aby*az*bc2*bcdx*bcdz*cx*cz + 
     -    abz*ax*bc2*bcdy*bcdz*cx*cz - aby*ax*bc2*bcdz**2*cx*cz + 
     -    abcz*az*bcdx*bcdz*bcx*cx*cz + abcz*ax*bcdz**2*bcx*cx*cz - 
     -    abz*bc2*bcdy*bcdz*bx*cx*cz + aby*bc2*bcdz**2*bx*cx*cz - 
     -    abcz*bcdz**2*bcx*bx*cx*cz + abcz*az*bc2*bcdx*by*cx*cz - 
     -    abz*bc2*bcdx*bcdy*bz*cx*cz + aby*bc2*bcdx*bcdz*bz*cx*cz - 
     -    abcz*bcdx*bcdz*bcx*bz*cx*cz - abcz*bc2*bcdx*by*bz*cx*cz + 
     -    abz*az*bc2*bcdy**2*cy*cz + abz*ay*bc2*bcdy*bcdz*cy*cz - 
     -    aby*az*bc2*bcdy*bcdz*cy*cz - aby*ay*bc2*bcdz**2*cy*cz + 
     -    abcz*az*bcdy*bcdz*bcx*cy*cz + abcz*ay*bcdz**2*bcx*cy*cz + 
     -    abcz*az*bc2*bcdy*by*cy*cz - abz*bc2*bcdy*bcdz*by*cy*cz + 
     -    aby*bc2*bcdz**2*by*cy*cz - abcz*bcdz**2*bcx*by*cy*cz - 
     -    abz*bc2*bcdy**2*bz*cy*cz + abcz*az*bc2*bcdz*bz*cy*cz + 
     -    aby*bc2*bcdy*bcdz*bz*cy*cz - abcz*bcdy*bcdz*bcx*bz*cy*cz - 
     -    abcz*bc2*bcdy*by*bz*cy*cz - abcz*bc2*bcdz*bz**2*cy*cz - 
     -    abz*ax*bc2*bcdx*bcdy*cz**2 - abz*ay*bc2*bcdy**2*cz**2 + 
     -    aby*ax*bc2*bcdx*bcdz*cz**2 + aby*ay*bc2*bcdy*bcdz*cz**2 - 
     -    abcz*ax*bcdx*bcdz*bcx*cz**2 - abcz*ay*bcdy*bcdz*bcx*cz**2 + 
     -    abz*bc2*bcdx*bcdy*bx*cz**2 - aby*bc2*bcdx*bcdz*bx*cz**2 + 
     -    abcz*bcdx*bcdz*bcx*bx*cz**2 - abcz*ax*bc2*bcdx*by*cz**2 - 
     -    abcz*ay*bc2*bcdy*by*cz**2 + abz*bc2*bcdy**2*by*cz**2 - 
     -    aby*bc2*bcdy*bcdz*by*cz**2 + abcz*bcdy*bcdz*bcx*by*cz**2 + 
     -    abcz*bc2*bcdx*bx*by*cz**2 + abcz*bc2*bcdy*by**2*cz**2 - 
     -    abcz*ay*bc2*bcdz*bz*cz**2 + abcz*bc2*bcdz*by*bz*cz**2 + 
     -    abcz*ay*bc2*bcdy*bx**2*dy - abcz*ay*bc2*bcdx*bx*by*dy - 
     -    abcz*ax*bc2*bcdy*bx*by*dy + abcz*ax*bc2*bcdx*by**2*dy - 
     -    abcz*az*bc2*bcdx*bx*bz*dy - abcz*az*bc2*bcdy*by*bz*dy + 
     -    abcz*ax*bc2*bcdx*bz**2*dy + abcz*ay*bc2*bcdy*bz**2*dy - 
     -    2*abcz*ay*bc2*bcdy*bx*cx*dy + abcz*ay*bc2*bcdx*by*cx*dy + 
     -    abcz*ax*bc2*bcdy*by*cx*dy + abcz*bc2*bcdy*bx*by*cx*dy - 
     -    abcz*bc2*bcdx*by**2*cx*dy + abcz*az*bc2*bcdx*bz*cx*dy - 
     -    abcz*bc2*bcdx*bz**2*cx*dy + abcz*ay*bc2*bcdy*cx**2*dy - 
     -    abcz*bc2*bcdy*by*cx**2*dy + abcz*ay*bc2*bcdx*bx*cy*dy + 
     -    abcz*ax*bc2*bcdy*bx*cy*dy - abcz*bc2*bcdy*bx**2*cy*dy - 
     -    2*abcz*ax*bc2*bcdx*by*cy*dy + abcz*bc2*bcdx*bx*by*cy*dy + 
     -    abcz*az*bc2*bcdy*bz*cy*dy - abcz*bc2*bcdy*bz**2*cy*dy - 
     -    abcz*ay*bc2*bcdx*cx*cy*dy - abcz*ax*bc2*bcdy*cx*cy*dy + 
     -    abcz*bc2*bcdy*bx*cx*cy*dy + abcz*bc2*bcdx*by*cx*cy*dy + 
     -    abcz*ax*bc2*bcdx*cy**2*dy - abcz*bc2*bcdx*bx*cy**2*dy + 
     -    abcz*az*bc2*bcdx*bx*cz*dy + abcz*az*bc2*bcdy*by*cz*dy - 
     -    2*abcz*ax*bc2*bcdx*bz*cz*dy - 2*abcz*ay*bc2*bcdy*bz*cz*dy + 
     -    abcz*bc2*bcdx*bx*bz*cz*dy + abcz*bc2*bcdy*by*bz*cz*dy - 
     -    abcz*az*bc2*bcdx*cx*cz*dy + abcz*bc2*bcdx*bz*cx*cz*dy - 
     -    abcz*az*bc2*bcdy*cy*cz*dy + abcz*bc2*bcdy*bz*cy*cz*dy + 
     -    abcz*ax*bc2*bcdx*cz**2*dy + abcz*ay*bc2*bcdy*cz**2*dy - 
     -    abcz*bc2*bcdx*bx*cz**2*dy - abcz*bc2*bcdy*by*cz**2*dy + 
     -    abcz*ay*bc2*bcdz*bx**2*dz - abcz*ax*bc2*bcdz*bx*by*dz - 
     -    abcz*az*bc2*bcdz*by*bz*dz + abcz*ay*bc2*bcdz*bz**2*dz - 
     -    2*abcz*ay*bc2*bcdz*bx*cx*dz + abcz*ax*bc2*bcdz*by*cx*dz + 
     -    abcz*bc2*bcdz*bx*by*cx*dz + abcz*ay*bc2*bcdz*cx**2*dz - 
     -    abcz*bc2*bcdz*by*cx**2*dz + abcz*ax*bc2*bcdz*bx*cy*dz - 
     -    abcz*bc2*bcdz*bx**2*cy*dz + abcz*az*bc2*bcdz*bz*cy*dz - 
     -    abcz*bc2*bcdz*bz**2*cy*dz - abcz*ax*bc2*bcdz*cx*cy*dz + 
     -    abcz*bc2*bcdz*bx*cx*cy*dz + abcz*az*bc2*bcdz*by*cz*dz - 
     -    2*abcz*ay*bc2*bcdz*bz*cz*dz + abcz*bc2*bcdz*by*bz*cz*dz - 
     -    abcz*az*bc2*bcdz*cy*cz*dz + abcz*bc2*bcdz*bz*cy*cz*dz + 
     -    abcz*ay*bc2*bcdz*cz**2*dz - abcz*bc2*bcdz*by*cz**2*dz - 
     -    abcy*(ax*bc2*bcdy**2*by + bc2*bcdy**2*bx*by - 
     -       ax*bcdy**2*bcx*bx*by - bc2*bcdx*bcdy*by**2 + 
     -       ax*bcdx*bcdy*bcx*by**2 + ax*bc2*bcdy*bcdz*bz + 
     -       bc2*bcdy*bcdz*bx*bz - ax*bcdy*bcdz*bcx*bx*bz + 
     -       ax*bc2*bcdy*bx*by*bz - ax*bc2*bcdx*by**2*bz - 
     -       bc2*bcdx*bcdy*bz**2 + ax*bcdx*bcdy*bcx*bz**2 + 
     -       ax*bc2*bcdz*bx*bz**2 - ax*bc2*bcdx*bz**3 - 
     -       2*bc2*bcdy**2*by*cx + ax*bcdy**2*bcx*by*cx + 
     -       bcdy**2*bcx*bx*by*cx - bcdx*bcdy*bcx*by**2*cx - 
     -       2*bc2*bcdy*bcdz*bz*cx + ax*bcdy*bcdz*bcx*bz*cx + 
     -       bcdy*bcdz*bcx*bx*bz*cx - ax*bc2*bcdy*by*bz*cx - 
     -       bc2*bcdy*bx*by*bz*cx + bc2*bcdx*by**2*bz*cx - 
     -       ax*bc2*bcdz*bz**2*cx - bcdx*bcdy*bcx*bz**2*cx - 
     -       bc2*bcdz*bx*bz**2*cx + bc2*bcdx*bz**3*cx - 
     -       bcdy**2*bcx*by*cx**2 - bcdy*bcdz*bcx*bz*cx**2 + 
     -       bc2*bcdy*by*bz*cx**2 + bc2*bcdz*bz**2*cx**2 - 
     -       ax*bc2*bcdy**2*cy + bc2*bcdy**2*bx*cy + 
     -       ax*bcdy**2*bcx*bx*cy - bcdy**2*bcx*bx**2*cy + 
     -       bc2*bcdx*bcdy*by*cy - 2*ax*bcdx*bcdy*bcx*by*cy + 
     -       bcdx*bcdy*bcx*bx*by*cy + 2*ax*bc2*bcdx*by*bz*cy + 
     -       bcdy*bcdz*bcx*by*bz*cy - bc2*bcdx*bx*by*bz*cy - 
     -       bc2*bcdy*by**2*bz*cy - bcdy**2*bcx*bz**2*cy - 
     -       bc2*bcdz*by*bz**2*cy - ax*bcdy**2*bcx*cx*cy + 
     -       bcdy**2*bcx*bx*cx*cy + bcdx*bcdy*bcx*by*cx*cy - 
     -       bc2*bcdx*by*bz*cx*cy + ax*bcdx*bcdy*bcx*cy**2 - 
     -       bcdx*bcdy*bcx*bx*cy**2 - ax*bc2*bcdx*bz*cy**2 - 
     -       bcdy*bcdz*bcx*bz*cy**2 + bc2*bcdx*bx*bz*cy**2 + 
     -       bc2*bcdy*by*bz*cy**2 + bc2*bcdz*bz**2*cy**2 - 
     -       ax*bc2*bcdy*bcdz*cz + bc2*bcdy*bcdz*bx*cz + 
     -       ax*bcdy*bcdz*bcx*bx*cz - bcdy*bcdz*bcx*bx**2*cz - 
     -       ax*bc2*bcdy*bx*by*cz + bc2*bcdy*bx**2*by*cz - 
     -       bcdy*bcdz*bcx*by**2*cz + bc2*bcdy*by**3*cz + 
     -       bc2*bcdx*bcdy*bz*cz - 2*ax*bcdx*bcdy*bcx*bz*cz - 
     -       ax*bc2*bcdz*bx*bz*cz + bcdx*bcdy*bcx*bx*bz*cz + 
     -       bc2*bcdz*bx**2*bz*cz + bcdy**2*bcx*by*bz*cz + 
     -       bc2*bcdz*by**2*bz*cz + 2*ax*bc2*bcdx*bz**2*cz - 
     -       bc2*bcdx*bx*bz**2*cz - ax*bcdy*bcdz*bcx*cx*cz + 
     -       bcdy*bcdz*bcx*bx*cx*cz + ax*bc2*bcdy*by*cx*cz - 
     -       bc2*bcdy*bx*by*cx*cz + ax*bc2*bcdz*bz*cx*cz + 
     -       bcdx*bcdy*bcx*bz*cx*cz - bc2*bcdz*bx*bz*cx*cz - 
     -       bc2*bcdx*bz**2*cx*cz + bcdy*bcdz*bcx*by*cy*cz - 
     -       bc2*bcdy*by**2*cy*cz + bcdy**2*bcx*bz*cy*cz - 
     -       bc2*bcdz*by*bz*cy*cz + ax*bcdx*bcdy*bcx*cz**2 - 
     -       bcdx*bcdy*bcx*bx*cz**2 - bcdy**2*bcx*by*cz**2 - 
     -       ax*bc2*bcdx*bz*cz**2 + bc2*bcdx*bx*bz*cz**2 - 
     -       ax*bc2*bcdy*bx*bz*dy + ax*bc2*bcdy*bz*cx*dy + 
     -       bc2*bcdy*bx*bz*cx*dy - bc2*bcdy*bz*cx**2*dy + 
     -       bc2*bcdy*by*bz*cy*dy - bc2*bcdy*bz*cy**2*dy + 
     -       ax*bc2*bcdy*bx*cz*dy - bc2*bcdy*bx**2*cz*dy - 
     -       bc2*bcdy*by**2*cz*dy - ax*bc2*bcdy*cx*cz*dy + 
     -       bc2*bcdy*bx*cx*cz*dy + bc2*bcdy*by*cy*cz*dy + 
     -       az*(bcdy*bcx*(bcdz*(bcx**2 + bcy**2) - 
     -             (abx*bcdx + bcdy*bcy)*bcz) + 
     -          bc2*(bcdy*(-2*abx*bcdz - by**3 + bcdx*bz - by*cx**2 + 
     -                2*by**2*cy - by*cy**2 - bcdx*cz + 
     -                2*bx*cx*(by - dy) + by**2*dy + cx**2*dy - 
     -                2*by*cy*dy + cy**2*dy + bx**2*(-by + dy)) - 
     -             (bcdz*(bcx**2 + bcy**2) - abx*bcdx*bcz)*(bz - dz)))
     -        + ay*(bcdy*bcx*
     -           (-(bcy*(abx*bcdx + bcdz*bcz)) + bcdy*(bcx**2 + bcz**2))
     -            + bc2*(-2*abx*bcdy**2 + 
     -             bcdy*bcy*(bcdx + bcz*(by - dy)) + 
     -             bcy*(abx*bcdx + bcdz*bcz)*(bz - dz))) + 
     -       ax*bc2*bcdx*by**2*dz - ax*bc2*bcdz*bx*bz*dz + 
     -       ax*bc2*bcdx*bz**2*dz - bc2*bcdx*by**2*cx*dz + 
     -       ax*bc2*bcdz*bz*cx*dz + bc2*bcdz*bx*bz*cx*dz - 
     -       bc2*bcdx*bz**2*cx*dz - bc2*bcdz*bz*cx**2*dz - 
     -       2*ax*bc2*bcdx*by*cy*dz + bc2*bcdx*bx*by*cy*dz + 
     -       bc2*bcdz*by*bz*cy*dz + bc2*bcdx*by*cx*cy*dz + 
     -       ax*bc2*bcdx*cy**2*dz - bc2*bcdx*bx*cy**2*dz - 
     -       bc2*bcdz*bz*cy**2*dz + ax*bc2*bcdz*bx*cz*dz - 
     -       bc2*bcdz*bx**2*cz*dz - bc2*bcdz*by**2*cz*dz - 
     -       2*ax*bc2*bcdx*bz*cz*dz + bc2*bcdx*bx*bz*cz*dz - 
     -       ax*bc2*bcdz*cx*cz*dz + bc2*bcdz*bx*cx*cz*dz + 
     -       bc2*bcdx*bz*cx*cz*dz + bc2*bcdz*by*cy*cz*dz + 
     -       ax*bc2*bcdx*cz**2*dz - bc2*bcdx*bx*cz**2*dz) - 
     -    abcx*bcdx*(ax*bc2*bcdy*by + bc2*bcdy*bx*by - 
     -       ax*bcdy*bcx*bx*by - bc2*bcdx*by**2 + ax*bcdx*bcx*by**2 + 
     -       ax*bc2*bcdz*bz + bc2*bcdz*bx*bz - ax*bcdz*bcx*bx*bz - 
     -       bc2*bcdx*bz**2 + ax*bcdx*bcx*bz**2 - 2*bc2*bcdy*by*cx + 
     -       ax*bcdy*bcx*by*cx + bcdy*bcx*bx*by*cx - 
     -       bcdx*bcx*by**2*cx - 2*bc2*bcdz*bz*cx + ax*bcdz*bcx*bz*cx + 
     -       bcdz*bcx*bx*bz*cx - bcdx*bcx*bz**2*cx - 
     -       bcdy*bcx*by*cx**2 - bcdz*bcx*bz*cx**2 - ax*bc2*bcdy*cy + 
     -       bc2*bcdy*bx*cy + ax*bcdy*bcx*bx*cy - bcdy*bcx*bx**2*cy + 
     -       bc2*bcdx*by*cy - 2*ax*bcdx*bcx*by*cy + bcdx*bcx*bx*by*cy + 
     -       ax*bc2*bx*bz*cy - bc2*bx**2*bz*cy + bcdz*bcx*by*bz*cy - 
     -       bc2*by**2*bz*cy - bcdy*bcx*bz**2*cy - bc2*bz**3*cy - 
     -       ax*bcdy*bcx*cx*cy + bcdy*bcx*bx*cx*cy + 
     -       bcdx*bcx*by*cx*cy - ax*bc2*bz*cx*cy + bc2*bx*bz*cx*cy + 
     -       ax*bcdx*bcx*cy**2 - bcdx*bcx*bx*cy**2 - 
     -       bcdz*bcx*bz*cy**2 + bc2*by*bz*cy**2 - ax*bc2*bcdz*cz + 
     -       bc2*bcdz*bx*cz + ax*bcdz*bcx*bx*cz - bcdz*bcx*bx**2*cz - 
     -       ax*bc2*bx*by*cz + bc2*bx**2*by*cz - bcdz*bcx*by**2*cz + 
     -       bc2*by**3*cz + bc2*bcdx*bz*cz - 2*ax*bcdx*bcx*bz*cz + 
     -       bcdx*bcx*bx*bz*cz + bcdy*bcx*by*bz*cz + bc2*by*bz**2*cz - 
     -       ax*bcdz*bcx*cx*cz + bcdz*bcx*bx*cx*cz + ax*bc2*by*cx*cz - 
     -       bc2*bx*by*cx*cz + bcdx*bcx*bz*cx*cz + bcdz*bcx*by*cy*cz - 
     -       bc2*by**2*cy*cz + bcdy*bcx*bz*cy*cz + bc2*bz**2*cy*cz + 
     -       ax*bcdx*bcx*cz**2 - bcdx*bcx*bx*cz**2 - 
     -       bcdy*bcx*by*cz**2 - bc2*by*bz*cz**2 - ax*bc2*bx*bz*dy + 
     -       ax*bc2*bz*cx*dy + bc2*bx*bz*cx*dy - bc2*bz*cx**2*dy + 
     -       bc2*by*bz*cy*dy - bc2*bz*cy**2*dy + ax*bc2*bx*cz*dy - 
     -       bc2*bx**2*cz*dy - bc2*by**2*cz*dy - ax*bc2*cx*cz*dy + 
     -       bc2*bx*cx*cz*dy + bc2*by*cy*cz*dy + ax*bc2*bx*by*dz - 
     -       ax*bc2*by*cx*dz - bc2*bx*by*cx*dz + bc2*by*cx**2*dz - 
     -       ax*bc2*bx*cy*dz + bc2*bx**2*cy*dz + bc2*bz**2*cy*dz + 
     -       ax*bc2*cx*cy*dz - bc2*bx*cx*cy*dz - bc2*by*bz*cz*dz - 
     -       bc2*bz*cy*cz*dz + bc2*by*cz**2*dz + 
     -       az*(bcx*(bcdz*(bcx**2 + bcy**2) - 
     -             (abx*bcdx + bcdy*bcy)*bcz) + 
     -          bc2*(-2*abx*bcdz - by**3 + bcdx*bz - by*bz**2 - 
     -             by*cx**2 + 2*by**2*cy + bz**2*cy - by*cy**2 - 
     -             bcdx*cz + by*bz*cz - bz*cy*cz + 2*bx*cx*(by - dy) + 
     -             by**2*dy + cx**2*dy - 2*by*cy*dy + cy**2*dy + 
     -             bx**2*(-by + dy) + by*bz*dz - bz*cy*dz - by*cz*dz + 
     -             cy*cz*dz)) + 
     -       ay*(bcx*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(bcdx*(-bx + cx) + bcdz*(-bz + cz))) + 
     -          bc2*(-2*bcdy*bx + bcdx*by + bx**2*bz + by**2*bz + 
     -             bz**3 + 2*bcdy*cx - 2*bx*bz*cx + bz*cx**2 - 
     -             bcdx*cy - by*bz*cy - by**2*cz - 2*bz**2*cz + 
     -             by*cy*cz + bz*cz**2 - by*bz*dy + bz*cy*dy + 
     -             by*cz*dy - cy*cz*dy - bx**2*dz - bz**2*dz + 
     -             2*bx*cx*dz - cx**2*dz + 2*bz*cz*dz - cz**2*dz))))/
     -  (bc2**1.5*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**2/
     -       bc2))


       dcy= (-(abcz*ay*bc2*bcdx*bcdz*bx) - abcz*ax*bc2*bcdy*bcdz*bx + 
     -    abz*ay*bc2*bcdx*bcdy*bx**2 + abz*az*bc2*bcdx*bcdz*bx**2 + 
     -    abcz*bc2*bcdy*bcdz*bx**2 - abx*ay*bc2*bcdy*bcdz*bx**2 - 
     -    abx*az*bc2*bcdz**2*bx**2 - abcz*ay*bcdy*bcdz*bcy*bx**2 - 
     -    abcz*az*bcdz**2*bcy*bx**2 + abcz*ay*bc2*bcdy*bx**3 + 
     -    2*abcz*ax*bc2*bcdx*bcdz*by + 2*abcz*az*bc2*bcdz**2*by - 
     -    abz*ay*bc2*bcdx**2*bx*by - abz*ax*bc2*bcdx*bcdy*bx*by - 
     -    abcz*bc2*bcdx*bcdz*bx*by + abx*ay*bc2*bcdx*bcdz*bx*by + 
     -    abx*ax*bc2*bcdy*bcdz*bx*by + abcz*ay*bcdx*bcdz*bcy*bx*by + 
     -    abcz*ax*bcdy*bcdz*bcy*bx*by - abcz*ay*bc2*bcdx*bx**2*by - 
     -    abcz*ax*bc2*bcdy*bx**2*by + abz*ax*bc2*bcdx**2*by**2 - 
     -    abx*ax*bc2*bcdx*bcdz*by**2 + abz*az*bc2*bcdx*bcdz*by**2 - 
     -    abx*az*bc2*bcdz**2*by**2 - abcz*ax*bcdx*bcdz*bcy*by**2 - 
     -    abcz*az*bcdz**2*bcy*by**2 + abcz*ax*bc2*bcdx*bx*by**2 - 
     -    abcz*az*bc2*bcdy*bcdz*bz - abcz*ay*bc2*bcdz**2*bz - 
     -    abz*az*bc2*bcdx**2*bx*bz - abz*ax*bc2*bcdx*bcdz*bx*bz + 
     -    abx*az*bc2*bcdx*bcdz*bx*bz + abx*ax*bc2*bcdz**2*bx*bz + 
     -    abcz*az*bcdx*bcdz*bcy*bx*bz + abcz*ax*bcdz**2*bcy*bx*bz - 
     -    abcz*az*bc2*bcdx*bx**2*bz - abz*az*bc2*bcdx*bcdy*by*bz - 
     -    abz*ay*bc2*bcdx*bcdz*by*bz + abx*az*bc2*bcdy*bcdz*by*bz - 
     -    abcz*bc2*bcdz**2*by*bz + abx*ay*bc2*bcdz**2*by*bz + 
     -    abcz*az*bcdy*bcdz*bcy*by*bz + abcz*ay*bcdz**2*bcy*by*bz - 
     -    abcz*az*bc2*bcdy*bx*by*bz - abcz*ay*bc2*bcdz*bx*by*bz + 
     -    abcz*ax*bc2*bcdz*by**2*bz + abz*ax*bc2*bcdx**2*bz**2 + 
     -    abz*ay*bc2*bcdx*bcdy*bz**2 - abx*ax*bc2*bcdx*bcdz*bz**2 + 
     -    abcz*bc2*bcdy*bcdz*bz**2 - abx*ay*bc2*bcdy*bcdz*bz**2 - 
     -    abcz*ax*bcdx*bcdz*bcy*bz**2 - abcz*ay*bcdy*bcdz*bcy*bz**2 + 
     -    abcz*ax*bc2*bcdx*bx*bz**2 + abcz*ay*bc2*bcdy*bx*bz**2 - 
     -    abcz*az*bc2*bcdz*bx*bz**2 + abcz*ax*bc2*bcdz*bz**3 + 
     -    abcz*ay*bc2*bcdx*bcdz*cx + abcz*ax*bc2*bcdy*bcdz*cx - 
     -    2*abz*ay*bc2*bcdx*bcdy*bx*cx - 
     -    2*abz*az*bc2*bcdx*bcdz*bx*cx - abcz*bc2*bcdy*bcdz*bx*cx + 
     -    2*abx*ay*bc2*bcdy*bcdz*bx*cx + 2*abx*az*bc2*bcdz**2*bx*cx + 
     -    2*abcz*ay*bcdy*bcdz*bcy*bx*cx + 
     -    2*abcz*az*bcdz**2*bcy*bx*cx - 2*abcz*ay*bc2*bcdy*bx**2*cx + 
     -    abz*ay*bc2*bcdx**2*by*cx + abz*ax*bc2*bcdx*bcdy*by*cx - 
     -    abcz*bc2*bcdx*bcdz*by*cx - abx*ay*bc2*bcdx*bcdz*by*cx - 
     -    abx*ax*bc2*bcdy*bcdz*by*cx - abcz*ay*bcdx*bcdz*bcy*by*cx - 
     -    abcz*ax*bcdy*bcdz*bcy*by*cx + abcz*ay*bc2*bcdx*bx*by*cx + 
     -    abcz*ax*bc2*bcdy*bx*by*cx + abz*bc2*bcdx*bcdy*bx*by*cx - 
     -    abx*bc2*bcdy*bcdz*bx*by*cx - abcz*bcdy*bcdz*bcy*bx*by*cx + 
     -    abcz*bc2*bcdy*bx**2*by*cx - abz*bc2*bcdx**2*by**2*cx + 
     -    abx*bc2*bcdx*bcdz*by**2*cx + abcz*bcdx*bcdz*bcy*by**2*cx - 
     -    abcz*bc2*bcdx*bx*by**2*cx + abz*az*bc2*bcdx**2*bz*cx + 
     -    abz*ax*bc2*bcdx*bcdz*bz*cx - abx*az*bc2*bcdx*bcdz*bz*cx - 
     -    abx*ax*bc2*bcdz**2*bz*cx - abcz*az*bcdx*bcdz*bcy*bz*cx - 
     -    abcz*ax*bcdz**2*bcy*bz*cx + abcz*az*bc2*bcdx*bx*bz*cx + 
     -    abz*bc2*bcdx*bcdz*bx*bz*cx - abx*bc2*bcdz**2*bx*bz*cx - 
     -    abcz*bcdz**2*bcy*bx*bz*cx + abcz*ay*bc2*bcdz*by*bz*cx - 
     -    abcz*bc2*bcdz*by**2*bz*cx - abz*bc2*bcdx**2*bz**2*cx + 
     -    abcz*az*bc2*bcdz*bz**2*cx + abx*bc2*bcdx*bcdz*bz**2*cx + 
     -    abcz*bcdx*bcdz*bcy*bz**2*cx - abcz*bc2*bcdx*bx*bz**2*cx - 
     -    abcz*bc2*bcdz*bz**3*cx + abz*ay*bc2*bcdx*bcdy*cx**2 + 
     -    abz*az*bc2*bcdx*bcdz*cx**2 - abx*ay*bc2*bcdy*bcdz*cx**2 - 
     -    abx*az*bc2*bcdz**2*cx**2 - abcz*ay*bcdy*bcdz*bcy*cx**2 - 
     -    abcz*az*bcdz**2*bcy*cx**2 + abcz*ay*bc2*bcdy*bx*cx**2 - 
     -    abz*bc2*bcdx*bcdy*by*cx**2 + abx*bc2*bcdy*bcdz*by*cx**2 + 
     -    abcz*bcdy*bcdz*bcy*by*cx**2 - abcz*bc2*bcdy*bx*by*cx**2 - 
     -    abz*bc2*bcdx*bcdz*bz*cx**2 + abx*bc2*bcdz**2*bz*cx**2 + 
     -    abcz*bcdz**2*bcy*bz*cx**2 - 2*abcz*ax*bc2*bcdx*bcdz*cy - 
     -    2*abcz*az*bc2*bcdz**2*cy + abz*ay*bc2*bcdx**2*bx*cy + 
     -    abz*ax*bc2*bcdx*bcdy*bx*cy + 2*abcz*bc2*bcdx*bcdz*bx*cy - 
     -    abx*ay*bc2*bcdx*bcdz*bx*cy - abx*ax*bc2*bcdy*bcdz*bx*cy - 
     -    abcz*ay*bcdx*bcdz*bcy*bx*cy - abcz*ax*bcdy*bcdz*bcy*bx*cy + 
     -    abcz*ay*bc2*bcdx*bx**2*cy + abcz*ax*bc2*bcdy*bx**2*cy - 
     -    abz*bc2*bcdx*bcdy*bx**2*cy + abx*bc2*bcdy*bcdz*bx**2*cy + 
     -    abcz*bcdy*bcdz*bcy*bx**2*cy - abcz*bc2*bcdy*bx**3*cy - 
     -    2*abz*ax*bc2*bcdx**2*by*cy + 2*abx*ax*bc2*bcdx*bcdz*by*cy - 
     -    2*abz*az*bc2*bcdx*bcdz*by*cy + 2*abx*az*bc2*bcdz**2*by*cy + 
     -    2*abcz*ax*bcdx*bcdz*bcy*by*cy + 
     -    2*abcz*az*bcdz**2*bcy*by*cy - 2*abcz*ax*bc2*bcdx*bx*by*cy + 
     -    abz*bc2*bcdx**2*bx*by*cy - abx*bc2*bcdx*bcdz*bx*by*cy - 
     -    abcz*bcdx*bcdz*bcy*bx*by*cy + abcz*bc2*bcdx*bx**2*by*cy + 
     -    abz*az*bc2*bcdx*bcdy*bz*cy + abz*ay*bc2*bcdx*bcdz*bz*cy - 
     -    abx*az*bc2*bcdy*bcdz*bz*cy + 2*abcz*bc2*bcdz**2*bz*cy - 
     -    abx*ay*bc2*bcdz**2*bz*cy - abcz*az*bcdy*bcdz*bcy*bz*cy - 
     -    abcz*ay*bcdz**2*bcy*bz*cy + abcz*az*bc2*bcdy*bx*bz*cy + 
     -    abcz*ay*bc2*bcdz*bx*bz*cy - 2*abcz*ax*bc2*bcdz*by*bz*cy + 
     -    abz*bc2*bcdx*bcdz*by*bz*cy - abx*bc2*bcdz**2*by*bz*cy - 
     -    abcz*bcdz**2*bcy*by*bz*cy + abcz*bc2*bcdz*bx*by*bz*cy - 
     -    abz*bc2*bcdx*bcdy*bz**2*cy + abx*bc2*bcdy*bcdz*bz**2*cy + 
     -    abcz*bcdy*bcdz*bcy*bz**2*cy - abcz*bc2*bcdy*bx*bz**2*cy - 
     -    abz*ay*bc2*bcdx**2*cx*cy - abz*ax*bc2*bcdx*bcdy*cx*cy + 
     -    abx*ay*bc2*bcdx*bcdz*cx*cy + abx*ax*bc2*bcdy*bcdz*cx*cy + 
     -    abcz*ay*bcdx*bcdz*bcy*cx*cy + abcz*ax*bcdy*bcdz*bcy*cx*cy - 
     -    abcz*ay*bc2*bcdx*bx*cx*cy - abcz*ax*bc2*bcdy*bx*cx*cy + 
     -    abz*bc2*bcdx*bcdy*bx*cx*cy - abx*bc2*bcdy*bcdz*bx*cx*cy - 
     -    abcz*bcdy*bcdz*bcy*bx*cx*cy + abcz*bc2*bcdy*bx**2*cx*cy + 
     -    abz*bc2*bcdx**2*by*cx*cy - abx*bc2*bcdx*bcdz*by*cx*cy - 
     -    abcz*bcdx*bcdz*bcy*by*cx*cy + abcz*bc2*bcdx*bx*by*cx*cy - 
     -    abcz*ay*bc2*bcdz*bz*cx*cy + abcz*bc2*bcdz*by*bz*cx*cy + 
     -    abz*ax*bc2*bcdx**2*cy**2 - abx*ax*bc2*bcdx*bcdz*cy**2 + 
     -    abz*az*bc2*bcdx*bcdz*cy**2 - abx*az*bc2*bcdz**2*cy**2 - 
     -    abcz*ax*bcdx*bcdz*bcy*cy**2 - abcz*az*bcdz**2*bcy*cy**2 + 
     -    abcz*ax*bc2*bcdx*bx*cy**2 - abz*bc2*bcdx**2*bx*cy**2 + 
     -    abx*bc2*bcdx*bcdz*bx*cy**2 + abcz*bcdx*bcdz*bcy*bx*cy**2 - 
     -    abcz*bc2*bcdx*bx**2*cy**2 + abcz*ax*bc2*bcdz*bz*cy**2 - 
     -    abz*bc2*bcdx*bcdz*bz*cy**2 + abx*bc2*bcdz**2*bz*cy**2 + 
     -    abcz*bcdz**2*bcy*bz*cy**2 - abcz*bc2*bcdz*bx*bz*cy**2 + 
     -    abcz*az*bc2*bcdy*bcdz*cz + abcz*ay*bc2*bcdz**2*cz + 
     -    abz*az*bc2*bcdx**2*bx*cz + abz*ax*bc2*bcdx*bcdz*bx*cz - 
     -    abx*az*bc2*bcdx*bcdz*bx*cz - abx*ax*bc2*bcdz**2*bx*cz - 
     -    abcz*az*bcdx*bcdz*bcy*bx*cz - abcz*ax*bcdz**2*bcy*bx*cz + 
     -    abcz*az*bc2*bcdx*bx**2*cz - abz*bc2*bcdx*bcdz*bx**2*cz + 
     -    abx*bc2*bcdz**2*bx**2*cz + abcz*bcdz**2*bcy*bx**2*cz + 
     -    abz*az*bc2*bcdx*bcdy*by*cz + abz*ay*bc2*bcdx*bcdz*by*cz - 
     -    abx*az*bc2*bcdy*bcdz*by*cz - abcz*bc2*bcdz**2*by*cz - 
     -    abx*ay*bc2*bcdz**2*by*cz - abcz*az*bcdy*bcdz*bcy*by*cz - 
     -    abcz*ay*bcdz**2*bcy*by*cz + abcz*az*bc2*bcdy*bx*by*cz - 
     -    abz*bc2*bcdx*bcdz*by**2*cz + abx*bc2*bcdz**2*by**2*cz + 
     -    abcz*bcdz**2*bcy*by**2*cz - 2*abz*ax*bc2*bcdx**2*bz*cz - 
     -    2*abz*ay*bc2*bcdx*bcdy*bz*cz + 
     -    2*abx*ax*bc2*bcdx*bcdz*bz*cz - abcz*bc2*bcdy*bcdz*bz*cz + 
     -    2*abx*ay*bc2*bcdy*bcdz*bz*cz + 
     -    2*abcz*ax*bcdx*bcdz*bcy*bz*cz + 
     -    2*abcz*ay*bcdy*bcdz*bcy*bz*cz - 
     -    2*abcz*ax*bc2*bcdx*bx*bz*cz + abz*bc2*bcdx**2*bx*bz*cz - 
     -    2*abcz*ay*bc2*bcdy*bx*bz*cz + abcz*az*bc2*bcdz*bx*bz*cz - 
     -    abx*bc2*bcdx*bcdz*bx*bz*cz - abcz*bcdx*bcdz*bcy*bx*bz*cz + 
     -    abcz*bc2*bcdx*bx**2*bz*cz + abz*bc2*bcdx*bcdy*by*bz*cz - 
     -    abx*bc2*bcdy*bcdz*by*bz*cz - abcz*bcdy*bcdz*bcy*by*bz*cz + 
     -    abcz*bc2*bcdy*bx*by*bz*cz - 2*abcz*ax*bc2*bcdz*bz**2*cz + 
     -    abcz*bc2*bcdz*bx*bz**2*cz - abz*az*bc2*bcdx**2*cx*cz - 
     -    abz*ax*bc2*bcdx*bcdz*cx*cz + abx*az*bc2*bcdx*bcdz*cx*cz + 
     -    abx*ax*bc2*bcdz**2*cx*cz + abcz*az*bcdx*bcdz*bcy*cx*cz + 
     -    abcz*ax*bcdz**2*bcy*cx*cz - abcz*az*bc2*bcdx*bx*cx*cz + 
     -    abz*bc2*bcdx*bcdz*bx*cx*cz - abx*bc2*bcdz**2*bx*cx*cz - 
     -    abcz*bcdz**2*bcy*bx*cx*cz + abz*bc2*bcdx**2*bz*cx*cz - 
     -    abcz*az*bc2*bcdz*bz*cx*cz - abx*bc2*bcdx*bcdz*bz*cx*cz - 
     -    abcz*bcdx*bcdz*bcy*bz*cx*cz + abcz*bc2*bcdx*bx*bz*cx*cz + 
     -    abcz*bc2*bcdz*bz**2*cx*cz - abz*az*bc2*bcdx*bcdy*cy*cz - 
     -    abz*ay*bc2*bcdx*bcdz*cy*cz + abx*az*bc2*bcdy*bcdz*cy*cz + 
     -    abx*ay*bc2*bcdz**2*cy*cz + abcz*az*bcdy*bcdz*bcy*cy*cz + 
     -    abcz*ay*bcdz**2*bcy*cy*cz - abcz*az*bc2*bcdy*bx*cy*cz + 
     -    abz*bc2*bcdx*bcdz*by*cy*cz - abx*bc2*bcdz**2*by*cy*cz - 
     -    abcz*bcdz**2*bcy*by*cy*cz + abz*bc2*bcdx*bcdy*bz*cy*cz - 
     -    abx*bc2*bcdy*bcdz*bz*cy*cz - abcz*bcdy*bcdz*bcy*bz*cy*cz + 
     -    abcz*bc2*bcdy*bx*bz*cy*cz + abz*ax*bc2*bcdx**2*cz**2 + 
     -    abz*ay*bc2*bcdx*bcdy*cz**2 - abx*ax*bc2*bcdx*bcdz*cz**2 - 
     -    abx*ay*bc2*bcdy*bcdz*cz**2 - abcz*ax*bcdx*bcdz*bcy*cz**2 - 
     -    abcz*ay*bcdy*bcdz*bcy*cz**2 + abcz*ax*bc2*bcdx*bx*cz**2 - 
     -    abz*bc2*bcdx**2*bx*cz**2 + abcz*ay*bc2*bcdy*bx*cz**2 + 
     -    abx*bc2*bcdx*bcdz*bx*cz**2 + abcz*bcdx*bcdz*bcy*bx*cz**2 - 
     -    abcz*bc2*bcdx*bx**2*cz**2 - abz*bc2*bcdx*bcdy*by*cz**2 + 
     -    abx*bc2*bcdy*bcdz*by*cz**2 + abcz*bcdy*bcdz*bcy*by*cz**2 - 
     -    abcz*bc2*bcdy*bx*by*cz**2 + abcz*ax*bc2*bcdz*bz*cz**2 - 
     -    abcz*bc2*bcdz*bx*bz*cz**2 - abcz*ay*bc2*bcdy*bx**2*dx + 
     -    abcz*ay*bc2*bcdx*bx*by*dx + abcz*ax*bc2*bcdy*bx*by*dx - 
     -    abcz*ax*bc2*bcdx*by**2*dx + abcz*az*bc2*bcdx*bx*bz*dx + 
     -    abcz*az*bc2*bcdy*by*bz*dx - abcz*ax*bc2*bcdx*bz**2*dx - 
     -    abcz*ay*bc2*bcdy*bz**2*dx + 2*abcz*ay*bc2*bcdy*bx*cx*dx - 
     -    abcz*ay*bc2*bcdx*by*cx*dx - abcz*ax*bc2*bcdy*by*cx*dx - 
     -    abcz*bc2*bcdy*bx*by*cx*dx + abcz*bc2*bcdx*by**2*cx*dx - 
     -    abcz*az*bc2*bcdx*bz*cx*dx + abcz*bc2*bcdx*bz**2*cx*dx - 
     -    abcz*ay*bc2*bcdy*cx**2*dx + abcz*bc2*bcdy*by*cx**2*dx - 
     -    abcz*ay*bc2*bcdx*bx*cy*dx - abcz*ax*bc2*bcdy*bx*cy*dx + 
     -    abcz*bc2*bcdy*bx**2*cy*dx + 2*abcz*ax*bc2*bcdx*by*cy*dx - 
     -    abcz*bc2*bcdx*bx*by*cy*dx - abcz*az*bc2*bcdy*bz*cy*dx + 
     -    abcz*bc2*bcdy*bz**2*cy*dx + abcz*ay*bc2*bcdx*cx*cy*dx + 
     -    abcz*ax*bc2*bcdy*cx*cy*dx - abcz*bc2*bcdy*bx*cx*cy*dx - 
     -    abcz*bc2*bcdx*by*cx*cy*dx - abcz*ax*bc2*bcdx*cy**2*dx + 
     -    abcz*bc2*bcdx*bx*cy**2*dx - abcz*az*bc2*bcdx*bx*cz*dx - 
     -    abcz*az*bc2*bcdy*by*cz*dx + 2*abcz*ax*bc2*bcdx*bz*cz*dx + 
     -    2*abcz*ay*bc2*bcdy*bz*cz*dx - abcz*bc2*bcdx*bx*bz*cz*dx - 
     -    abcz*bc2*bcdy*by*bz*cz*dx + abcz*az*bc2*bcdx*cx*cz*dx - 
     -    abcz*bc2*bcdx*bz*cx*cz*dx + abcz*az*bc2*bcdy*cy*cz*dx - 
     -    abcz*bc2*bcdy*bz*cy*cz*dx - abcz*ax*bc2*bcdx*cz**2*dx - 
     -    abcz*ay*bc2*bcdy*cz**2*dx + abcz*bc2*bcdx*bx*cz**2*dx + 
     -    abcz*bc2*bcdy*by*cz**2*dx + abcz*ay*bc2*bcdz*bx*by*dz - 
     -    abcz*ax*bc2*bcdz*by**2*dz + abcz*az*bc2*bcdz*bx*bz*dz - 
     -    abcz*ax*bc2*bcdz*bz**2*dz - abcz*ay*bc2*bcdz*by*cx*dz + 
     -    abcz*bc2*bcdz*by**2*cx*dz - abcz*az*bc2*bcdz*bz*cx*dz + 
     -    abcz*bc2*bcdz*bz**2*cx*dz - abcz*ay*bc2*bcdz*bx*cy*dz + 
     -    2*abcz*ax*bc2*bcdz*by*cy*dz - abcz*bc2*bcdz*bx*by*cy*dz + 
     -    abcz*ay*bc2*bcdz*cx*cy*dz - abcz*bc2*bcdz*by*cx*cy*dz - 
     -    abcz*ax*bc2*bcdz*cy**2*dz + abcz*bc2*bcdz*bx*cy**2*dz - 
     -    abcz*az*bc2*bcdz*bx*cz*dz + 2*abcz*ax*bc2*bcdz*bz*cz*dz - 
     -    abcz*bc2*bcdz*bx*bz*cz*dz + abcz*az*bc2*bcdz*cx*cz*dz - 
     -    abcz*bc2*bcdz*bz*cx*cz*dz - abcz*ax*bc2*bcdz*cz**2*dz + 
     -    abcz*bc2*bcdz*bx*cz**2*dz - 
     -    abcx*(-(bc2*bcdx*bcdy*bx**2) + az*bcdx*bcdz*bcy*bx**2 + 
     -       az*bc2*bcdx*bx**3 - 2*az*bc2*bcdx*bcdz*by + 
     -       bc2*bcdx**2*bx*by + az*bcdx*bcdz*bcy*by**2 + 
     -       az*bc2*bcdx*bx*by**2 + az*bc2*bcdx*bcdy*bz - 
     -       az*bcdx**2*bcy*bx*bz + az*bc2*bcdz*bx**2*bz + 
     -       bc2*bcdx*bcdz*by*bz - az*bcdx*bcdy*bcy*by*bz + 
     -       az*bc2*bcdz*by**2*bz - bc2*bcdx*bcdy*bz**2 - 
     -       az*bc2*bcdy*by*bz**2 + bc2*bcdx*bcdy*bx*cx - 
     -       2*az*bcdx*bcdz*bcy*bx*cx - 2*az*bc2*bcdx*bx**2*cx + 
     -       bc2*bcdx**2*by*cx + bcdx*bcdy*bcy*bx*by*cx - 
     -       bcdx**2*bcy*by**2*cx + az*bcdx**2*bcy*bz*cx - 
     -       2*az*bc2*bcdz*bx*bz*cx + bcdx*bcdz*bcy*bx*bz*cx + 
     -       bc2*bcdx*bx**2*bz*cx + bc2*bcdy*bx*by*bz*cx - 
     -       bcdx**2*bcy*bz**2*cx + bc2*bcdz*bx*bz**2*cx + 
     -       az*bcdx*bcdz*bcy*cx**2 + az*bc2*bcdx*bx*cx**2 - 
     -       bcdx*bcdy*bcy*by*cx**2 + az*bc2*bcdz*bz*cx**2 - 
     -       bcdx*bcdz*bcy*bz*cx**2 - bc2*bcdx*bx*bz*cx**2 - 
     -       bc2*bcdy*by*bz*cx**2 - bc2*bcdz*bz**2*cx**2 + 
     -       2*az*bc2*bcdx*bcdz*cy - 2*bc2*bcdx**2*bx*cy - 
     -       bcdx*bcdy*bcy*bx**2*cy - 2*az*bcdx*bcdz*bcy*by*cy - 
     -       2*az*bc2*bcdx*bx*by*cy + bcdx**2*bcy*bx*by*cy - 
     -       2*bc2*bcdx*bcdz*bz*cy + az*bcdx*bcdy*bcy*bz*cy - 
     -       bc2*bcdy*bx**2*bz*cy - 2*az*bc2*bcdz*by*bz*cy + 
     -       bcdx*bcdz*bcy*by*bz*cy + bc2*bcdx*bx*by*bz*cy + 
     -       az*bc2*bcdy*bz**2*cy - bcdx*bcdy*bcy*bz**2*cy + 
     -       bc2*bcdz*by*bz**2*cy - bc2*bcdy*bz**3*cy + 
     -       bcdx*bcdy*bcy*bx*cx*cy + bcdx**2*bcy*by*cx*cy + 
     -       bc2*bcdy*bx*bz*cx*cy + az*bcdx*bcdz*bcy*cy**2 + 
     -       az*bc2*bcdx*bx*cy**2 - bcdx**2*bcy*bx*cy**2 + 
     -       az*bc2*bcdz*bz*cy**2 - bcdx*bcdz*bcy*bz*cy**2 - 
     -       bc2*bcdx*bx*bz*cy**2 - bc2*bcdz*bz**2*cy**2 - 
     -       az*bc2*bcdx*bcdy*cz + az*bcdx**2*bcy*bx*cz - 
     -       bcdx*bcdz*bcy*bx**2*cz - bc2*bcdx*bx**3*cz + 
     -       bc2*bcdx*bcdz*by*cz + az*bcdx*bcdy*bcy*by*cz - 
     -       bcdx*bcdz*bcy*by**2*cz - bc2*bcdx*bx*by**2*cz + 
     -       bc2*bcdx*bcdy*bz*cz + bcdx**2*bcy*bx*bz*cz - 
     -       bc2*bcdz*bx**2*bz*cz + az*bc2*bcdy*by*bz*cz + 
     -       bcdx*bcdy*bcy*by*bz*cz - bc2*bcdz*by**2*bz*cz + 
     -       bc2*bcdy*by*bz**2*cz - az*bcdx**2*bcy*cx*cz + 
     -       bcdx*bcdz*bcy*bx*cx*cz + bc2*bcdx*bx**2*cx*cz + 
     -       bcdx**2*bcy*bz*cx*cz + bc2*bcdz*bx*bz*cx*cz - 
     -       az*bcdx*bcdy*bcy*cy*cz + bcdx*bcdz*bcy*by*cy*cz + 
     -       bc2*bcdx*bx*by*cy*cz - az*bc2*bcdy*bz*cy*cz + 
     -       bcdx*bcdy*bcy*bz*cy*cz + bc2*bcdz*by*bz*cy*cz + 
     -       bc2*bcdy*bz**2*cy*cz - bcdx**2*bcy*bx*cz**2 - 
     -       bcdx*bcdy*bcy*by*cz**2 - bc2*bcdy*by*bz*cz**2 - 
     -       az*bc2*bcdx*bx**2*dx - az*bc2*bcdx*by**2*dx + 
     -       2*az*bc2*bcdx*bx*cx*dx - bc2*bcdx*bx*bz*cx*dx - 
     -       az*bc2*bcdx*cx**2*dx + bc2*bcdx*bz*cx**2*dx + 
     -       2*az*bc2*bcdx*by*cy*dx - bc2*bcdx*by*bz*cy*dx - 
     -       az*bc2*bcdx*cy**2*dx + bc2*bcdx*bz*cy**2*dx + 
     -       bc2*bcdx*bx**2*cz*dx + bc2*bcdx*by**2*cz*dx - 
     -       bc2*bcdx*bx*cx*cz*dx - bc2*bcdx*by*cy*cz*dx + 
     -       ax*(bcdx**2*bcy*(-2*bc2 + bcy**2 + bcz**2) + 
     -          abx*(-(bcdx*bcy*(bcdy*bcy + bcdz*bcz)) + 
     -             bc2*(bcdx*(bcdy - bcz*bx + bcz*dx) - 
     -                (bcdy*bcy + bcdz*bcz)*(bz - dz)))) + 
     -       ay*(bcdx*bcy*(-(bcy*(abx*bcdx + bcdz*bcz)) + 
     -             bcdy*(bcx**2 + bcz**2)) + 
     -          bc2*(abx*bcdx**2 + bcdx*bcz*(bcdz + bcy*(-bx + dx)) + 
     -             (-(bcdz*bcy*bcz) + bcdy*(bcx**2 + bcz**2))*
     -              (bz - dz))) - az*bc2*bcdz*bx**2*dz - 
     -       az*bc2*bcdz*by**2*dz + az*bc2*bcdy*by*bz*dz + 
     -       2*az*bc2*bcdz*bx*cx*dz - bc2*bcdy*bx*by*cx*dz - 
     -       bc2*bcdz*bx*bz*cx*dz - az*bc2*bcdz*cx**2*dz + 
     -       bc2*bcdy*by*cx**2*dz + bc2*bcdz*bz*cx**2*dz + 
     -       bc2*bcdy*bx**2*cy*dz + 2*az*bc2*bcdz*by*cy*dz - 
     -       az*bc2*bcdy*bz*cy*dz - bc2*bcdz*by*bz*cy*dz + 
     -       bc2*bcdy*bz**2*cy*dz - bc2*bcdy*bx*cx*cy*dz - 
     -       az*bc2*bcdz*cy**2*dz + bc2*bcdz*bz*cy**2*dz + 
     -       bc2*bcdz*bx**2*cz*dz - az*bc2*bcdy*by*cz*dz + 
     -       bc2*bcdz*by**2*cz*dz - bc2*bcdy*by*bz*cz*dz - 
     -       bc2*bcdz*bx*cx*cz*dz + az*bc2*bcdy*cy*cz*dz - 
     -       bc2*bcdz*by*cy*cz*dz - bc2*bcdy*bz*cy*cz*dz + 
     -       bc2*bcdy*by*cz**2*dz) - 
     -    abcy*bcdy*(-(bc2*bcdy*bx**2) + az*bcdz*bcy*bx**2 + 
     -       az*bc2*bx**3 - 2*az*bc2*bcdz*by + bc2*bcdx*bx*by + 
     -       az*bcdz*bcy*by**2 + az*bc2*bx*by**2 + az*bc2*bcdy*bz - 
     -       az*bcdx*bcy*bx*bz + bc2*bcdz*by*bz - az*bcdy*bcy*by*bz - 
     -       bc2*bcdy*bz**2 + az*bc2*bx*bz**2 + bc2*bcdy*bx*cx - 
     -       2*az*bcdz*bcy*bx*cx - 2*az*bc2*bx**2*cx + 
     -       bc2*bcdx*by*cx + bcdy*bcy*bx*by*cx - bcdx*bcy*by**2*cx + 
     -       az*bcdx*bcy*bz*cx + bcdz*bcy*bx*bz*cx + 
     -       bc2*bx**2*bz*cx + bc2*by**2*bz*cx - az*bc2*bz**2*cx - 
     -       bcdx*bcy*bz**2*cx + bc2*bz**3*cx + az*bcdz*bcy*cx**2 + 
     -       az*bc2*bx*cx**2 - bcdy*bcy*by*cx**2 - 
     -       bcdz*bcy*bz*cx**2 - bc2*bx*bz*cx**2 + 2*az*bc2*bcdz*cy - 
     -       2*bc2*bcdx*bx*cy - bcdy*bcy*bx**2*cy - 
     -       2*az*bcdz*bcy*by*cy - 2*az*bc2*bx*by*cy + 
     -       bcdx*bcy*bx*by*cy - 2*bc2*bcdz*bz*cy + 
     -       az*bcdy*bcy*bz*cy + bcdz*bcy*by*bz*cy - 
     -       bcdy*bcy*bz**2*cy + bcdy*bcy*bx*cx*cy + 
     -       bcdx*bcy*by*cx*cy - bc2*by*bz*cx*cy + 
     -       az*bcdz*bcy*cy**2 + az*bc2*bx*cy**2 - 
     -       bcdx*bcy*bx*cy**2 - bcdz*bcy*bz*cy**2 - az*bc2*bcdy*cz + 
     -       az*bcdx*bcy*bx*cz - bcdz*bcy*bx**2*cz - bc2*bx**3*cz + 
     -       bc2*bcdz*by*cz + az*bcdy*bcy*by*cz - bcdz*bcy*by**2*cz - 
     -       bc2*bx*by**2*cz + bc2*bcdy*bz*cz - az*bc2*bx*bz*cz + 
     -       bcdx*bcy*bx*bz*cz + bcdy*bcy*by*bz*cz - 
     -       bc2*bx*bz**2*cz - az*bcdx*bcy*cx*cz + 
     -       bcdz*bcy*bx*cx*cz + bc2*bx**2*cx*cz + az*bc2*bz*cx*cz + 
     -       bcdx*bcy*bz*cx*cz - bc2*bz**2*cx*cz - 
     -       az*bcdy*bcy*cy*cz + bcdz*bcy*by*cy*cz + 
     -       bc2*bx*by*cy*cz + bcdy*bcy*bz*cy*cz - 
     -       bcdx*bcy*bx*cz**2 - bcdy*bcy*by*cz**2 + 
     -       bc2*bx*bz*cz**2 - az*bc2*bx**2*dx - az*bc2*by**2*dx + 
     -       2*az*bc2*bx*cx*dx - bc2*bx*bz*cx*dx - az*bc2*cx**2*dx + 
     -       bc2*bz*cx**2*dx + 2*az*bc2*by*cy*dx - bc2*by*bz*cy*dx - 
     -       az*bc2*cy**2*dx + bc2*bz*cy**2*dx + bc2*bx**2*cz*dx + 
     -       bc2*by**2*cz*dx - bc2*bx*cx*cz*dx - bc2*by*cy*cz*dx - 
     -       az*bc2*bx*bz*dz - bc2*by**2*cx*dz + az*bc2*bz*cx*dz - 
     -       bc2*bz**2*cx*dz + bc2*bx*by*cy*dz + bc2*by*cx*cy*dz - 
     -       bc2*bx*cy**2*dz + az*bc2*bx*cz*dz + bc2*bx*bz*cz*dz - 
     -       az*bc2*cx*cz*dz + bc2*bz*cx*cz*dz - bc2*bx*cz**2*dz + 
     -       ax*(abx*(bc2*bcdy - bcy*(bcdy*bcy + bcdz*bcz)) + 
     -          bcdx*(bcy**3 + bcy*bcz**2 + 2*bc2*(-by + cy)) + 
     -          bc2*(-bz**3 - bz*cy**2 + 2*bz**2*cz - bz*cz**2 + 
     -             bx**2*(-bz + cz) - bz*cx*dx + cx*cz*dx + 
     -             bx*(bz - cz)*(cx + dx) + 2*by*cy*(bz - dz) + 
     -             bz**2*dz + cy**2*dz - 2*bz*cz*dz + cz**2*dz + 
     -             by**2*(-bz + dz))) + 
     -       ay*(bcy*(-(bcy*(abx*bcdx + bcdz*bcz)) + 
     -             bcdy*(bcx**2 + bcz**2)) + 
     -          bc2*(abx*bcdx + bcdz*bcz - 
     -             bcy*(bz*cdx - bx*cz + cz*dx + bx*dz - cx*dz)))))/
     -  (bc2**1.5*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**2
     -        /bc2))



       dcz= (-(abcz*az*bc2*bcdx*bcdz*bx) - abcz*ax*bc2*bcdz**2*bx - 
     -    aby*ay*bc2*bcdx*bcdy*bx**2 + abx*ay*bc2*bcdy**2*bx**2 - 
     -    aby*az*bc2*bcdx*bcdz*bx**2 + abx*az*bc2*bcdy*bcdz*bx**2 + 
     -    abcz*bc2*bcdz**2*bx**2 - abcz*ay*bcdy*bcdz*bcz*bx**2 - 
     -    abcz*az*bcdz**2*bcz*bx**2 + abcz*ay*bc2*bcdz*bx**3 - 
     -    abcz*az*bc2*bcdy*bcdz*by - abcz*ay*bc2*bcdz**2*by + 
     -    aby*ay*bc2*bcdx**2*bx*by + aby*ax*bc2*bcdx*bcdy*bx*by - 
     -    abx*ay*bc2*bcdx*bcdy*bx*by - abx*ax*bc2*bcdy**2*bx*by + 
     -    abcz*ay*bcdx*bcdz*bcz*bx*by + abcz*ax*bcdy*bcdz*bcz*bx*by - 
     -    abcz*ax*bc2*bcdz*bx**2*by - aby*ax*bc2*bcdx**2*by**2 + 
     -    abx*ax*bc2*bcdx*bcdy*by**2 - aby*az*bc2*bcdx*bcdz*by**2 + 
     -    abx*az*bc2*bcdy*bcdz*by**2 + abcz*bc2*bcdz**2*by**2 - 
     -    abcz*ax*bcdx*bcdz*bcz*by**2 - abcz*az*bcdz**2*bcz*by**2 + 
     -    abcz*ay*bc2*bcdz*bx*by**2 - abcz*ax*bc2*bcdz*by**3 + 
     -    2*abcz*ax*bc2*bcdx*bcdz*bz + 2*abcz*ay*bc2*bcdy*bcdz*bz + 
     -    aby*az*bc2*bcdx**2*bx*bz - abx*az*bc2*bcdx*bcdy*bx*bz - 
     -    abcz*bc2*bcdx*bcdz*bx*bz + aby*ax*bc2*bcdx*bcdz*bx*bz - 
     -    abx*ax*bc2*bcdy*bcdz*bx*bz + abcz*az*bcdx*bcdz*bcz*bx*bz + 
     -    abcz*ax*bcdz**2*bcz*bx*bz + aby*az*bc2*bcdx*bcdy*by*bz - 
     -    abx*az*bc2*bcdy**2*by*bz + aby*ay*bc2*bcdx*bcdz*by*bz - 
     -    abcz*bc2*bcdy*bcdz*by*bz - abx*ay*bc2*bcdy*bcdz*by*bz + 
     -    abcz*az*bcdy*bcdz*bcz*by*bz + abcz*ay*bcdz**2*bcz*by*bz - 
     -    aby*ax*bc2*bcdx**2*bz**2 + abx*ax*bc2*bcdx*bcdy*bz**2 - 
     -    aby*ay*bc2*bcdx*bcdy*bz**2 + abx*ay*bc2*bcdy**2*bz**2 - 
     -    abcz*ax*bcdx*bcdz*bcz*bz**2 - abcz*ay*bcdy*bcdz*bcz*bz**2 + 
     -    abcz*ay*bc2*bcdz*bx*bz**2 - abcz*ax*bc2*bcdz*by*bz**2 + 
     -    abcz*az*bc2*bcdx*bcdz*cx + abcz*ax*bc2*bcdz**2*cx + 
     -    2*aby*ay*bc2*bcdx*bcdy*bx*cx - 2*abx*ay*bc2*bcdy**2*bx*cx + 
     -    2*aby*az*bc2*bcdx*bcdz*bx*cx - 
     -    2*abx*az*bc2*bcdy*bcdz*bx*cx - abcz*bc2*bcdz**2*bx*cx + 
     -    2*abcz*ay*bcdy*bcdz*bcz*bx*cx + 
     -    2*abcz*az*bcdz**2*bcz*bx*cx - 2*abcz*ay*bc2*bcdz*bx**2*cx - 
     -    aby*ay*bc2*bcdx**2*by*cx - aby*ax*bc2*bcdx*bcdy*by*cx + 
     -    abx*ay*bc2*bcdx*bcdy*by*cx + abx*ax*bc2*bcdy**2*by*cx - 
     -    abcz*ay*bcdx*bcdz*bcz*by*cx - abcz*ax*bcdy*bcdz*bcz*by*cx - 
     -    aby*bc2*bcdx*bcdy*bx*by*cx + abx*bc2*bcdy**2*bx*by*cx + 
     -    abcz*ax*bc2*bcdz*bx*by*cx - abcz*bcdy*bcdz*bcz*bx*by*cx + 
     -    abcz*bc2*bcdz*bx**2*by*cx + aby*bc2*bcdx**2*by**2*cx - 
     -    abx*bc2*bcdx*bcdy*by**2*cx - abcz*ay*bc2*bcdz*by**2*cx + 
     -    abcz*bcdx*bcdz*bcz*by**2*cx + abcz*bc2*bcdz*by**3*cx - 
     -    aby*az*bc2*bcdx**2*bz*cx + abx*az*bc2*bcdx*bcdy*bz*cx - 
     -    abcz*bc2*bcdx*bcdz*bz*cx - aby*ax*bc2*bcdx*bcdz*bz*cx + 
     -    abx*ax*bc2*bcdy*bcdz*bz*cx - abcz*az*bcdx*bcdz*bcz*bz*cx - 
     -    abcz*ax*bcdz**2*bcz*bz*cx - aby*bc2*bcdx*bcdz*bx*bz*cx + 
     -    abx*bc2*bcdy*bcdz*bx*bz*cx - abcz*bcdz**2*bcz*bx*bz*cx - 
     -    abcz*az*bc2*bcdz*by*bz*cx + aby*bc2*bcdx**2*bz**2*cx - 
     -    abx*bc2*bcdx*bcdy*bz**2*cx + abcz*bcdx*bcdz*bcz*bz**2*cx + 
     -    abcz*bc2*bcdz*by*bz**2*cx - aby*ay*bc2*bcdx*bcdy*cx**2 + 
     -    abx*ay*bc2*bcdy**2*cx**2 - aby*az*bc2*bcdx*bcdz*cx**2 + 
     -    abx*az*bc2*bcdy*bcdz*cx**2 - abcz*ay*bcdy*bcdz*bcz*cx**2 - 
     -    abcz*az*bcdz**2*bcz*cx**2 + abcz*ay*bc2*bcdz*bx*cx**2 + 
     -    aby*bc2*bcdx*bcdy*by*cx**2 - abx*bc2*bcdy**2*by*cx**2 + 
     -    abcz*bcdy*bcdz*bcz*by*cx**2 - abcz*bc2*bcdz*bx*by*cx**2 + 
     -    aby*bc2*bcdx*bcdz*bz*cx**2 - abx*bc2*bcdy*bcdz*bz*cx**2 + 
     -    abcz*bcdz**2*bcz*bz*cx**2 + abcz*az*bc2*bcdy*bcdz*cy + 
     -    abcz*ay*bc2*bcdz**2*cy - aby*ay*bc2*bcdx**2*bx*cy - 
     -    aby*ax*bc2*bcdx*bcdy*bx*cy + abx*ay*bc2*bcdx*bcdy*bx*cy + 
     -    abx*ax*bc2*bcdy**2*bx*cy - abcz*ay*bcdx*bcdz*bcz*bx*cy - 
     -    abcz*ax*bcdy*bcdz*bcz*bx*cy + aby*bc2*bcdx*bcdy*bx**2*cy - 
     -    abx*bc2*bcdy**2*bx**2*cy + abcz*ax*bc2*bcdz*bx**2*cy + 
     -    abcz*bcdy*bcdz*bcz*bx**2*cy - abcz*bc2*bcdz*bx**3*cy + 
     -    2*aby*ax*bc2*bcdx**2*by*cy - 2*abx*ax*bc2*bcdx*bcdy*by*cy + 
     -    2*aby*az*bc2*bcdx*bcdz*by*cy - 
     -    2*abx*az*bc2*bcdy*bcdz*by*cy - abcz*bc2*bcdz**2*by*cy + 
     -    2*abcz*ax*bcdx*bcdz*bcz*by*cy + 
     -    2*abcz*az*bcdz**2*bcz*by*cy - aby*bc2*bcdx**2*bx*by*cy + 
     -    abx*bc2*bcdx*bcdy*bx*by*cy - abcz*ay*bc2*bcdz*bx*by*cy - 
     -    abcz*bcdx*bcdz*bcz*bx*by*cy + 2*abcz*ax*bc2*bcdz*by**2*cy - 
     -    abcz*bc2*bcdz*bx*by**2*cy - aby*az*bc2*bcdx*bcdy*bz*cy + 
     -    abx*az*bc2*bcdy**2*bz*cy - aby*ay*bc2*bcdx*bcdz*bz*cy - 
     -    abcz*bc2*bcdy*bcdz*bz*cy + abx*ay*bc2*bcdy*bcdz*bz*cy - 
     -    abcz*az*bcdy*bcdz*bcz*bz*cy - abcz*ay*bcdz**2*bcz*bz*cy + 
     -    abcz*az*bc2*bcdz*bx*bz*cy - aby*bc2*bcdx*bcdz*by*bz*cy + 
     -    abx*bc2*bcdy*bcdz*by*bz*cy - abcz*bcdz**2*bcz*by*bz*cy + 
     -    aby*bc2*bcdx*bcdy*bz**2*cy - abx*bc2*bcdy**2*bz**2*cy + 
     -    abcz*bcdy*bcdz*bcz*bz**2*cy - abcz*bc2*bcdz*bx*bz**2*cy + 
     -    aby*ay*bc2*bcdx**2*cx*cy + aby*ax*bc2*bcdx*bcdy*cx*cy - 
     -    abx*ay*bc2*bcdx*bcdy*cx*cy - abx*ax*bc2*bcdy**2*cx*cy + 
     -    abcz*ay*bcdx*bcdz*bcz*cx*cy + abcz*ax*bcdy*bcdz*bcz*cx*cy - 
     -    aby*bc2*bcdx*bcdy*bx*cx*cy + abx*bc2*bcdy**2*bx*cx*cy - 
     -    abcz*ax*bc2*bcdz*bx*cx*cy - abcz*bcdy*bcdz*bcz*bx*cx*cy + 
     -    abcz*bc2*bcdz*bx**2*cx*cy - aby*bc2*bcdx**2*by*cx*cy + 
     -    abx*bc2*bcdx*bcdy*by*cx*cy + abcz*ay*bc2*bcdz*by*cx*cy - 
     -    abcz*bcdx*bcdz*bcz*by*cx*cy - abcz*bc2*bcdz*by**2*cx*cy - 
     -    aby*ax*bc2*bcdx**2*cy**2 + abx*ax*bc2*bcdx*bcdy*cy**2 - 
     -    aby*az*bc2*bcdx*bcdz*cy**2 + abx*az*bc2*bcdy*bcdz*cy**2 - 
     -    abcz*ax*bcdx*bcdz*bcz*cy**2 - abcz*az*bcdz**2*bcz*cy**2 + 
     -    aby*bc2*bcdx**2*bx*cy**2 - abx*bc2*bcdx*bcdy*bx*cy**2 + 
     -    abcz*bcdx*bcdz*bcz*bx*cy**2 - abcz*ax*bc2*bcdz*by*cy**2 + 
     -    abcz*bc2*bcdz*bx*by*cy**2 + aby*bc2*bcdx*bcdz*bz*cy**2 - 
     -    abx*bc2*bcdy*bcdz*bz*cy**2 + abcz*bcdz**2*bcz*bz*cy**2 - 
     -    2*abcz*ax*bc2*bcdx*bcdz*cz - 2*abcz*ay*bc2*bcdy*bcdz*cz - 
     -    aby*az*bc2*bcdx**2*bx*cz + abx*az*bc2*bcdx*bcdy*bx*cz + 
     -    2*abcz*bc2*bcdx*bcdz*bx*cz - aby*ax*bc2*bcdx*bcdz*bx*cz + 
     -    abx*ax*bc2*bcdy*bcdz*bx*cz - abcz*az*bcdx*bcdz*bcz*bx*cz - 
     -    abcz*ax*bcdz**2*bcz*bx*cz + aby*bc2*bcdx*bcdz*bx**2*cz - 
     -    abx*bc2*bcdy*bcdz*bx**2*cz + abcz*bcdz**2*bcz*bx**2*cz - 
     -    aby*az*bc2*bcdx*bcdy*by*cz + abx*az*bc2*bcdy**2*by*cz - 
     -    aby*ay*bc2*bcdx*bcdz*by*cz + 2*abcz*bc2*bcdy*bcdz*by*cz + 
     -    abx*ay*bc2*bcdy*bcdz*by*cz - abcz*az*bcdy*bcdz*bcz*by*cz - 
     -    abcz*ay*bcdz**2*bcz*by*cz + aby*bc2*bcdx*bcdz*by**2*cz - 
     -    abx*bc2*bcdy*bcdz*by**2*cz + abcz*bcdz**2*bcz*by**2*cz + 
     -    2*aby*ax*bc2*bcdx**2*bz*cz - 2*abx*ax*bc2*bcdx*bcdy*bz*cz + 
     -    2*aby*ay*bc2*bcdx*bcdy*bz*cz - 2*abx*ay*bc2*bcdy**2*bz*cz + 
     -    2*abcz*ax*bcdx*bcdz*bcz*bz*cz + 
     -    2*abcz*ay*bcdy*bcdz*bcz*bz*cz - aby*bc2*bcdx**2*bx*bz*cz + 
     -    abx*bc2*bcdx*bcdy*bx*bz*cz - 2*abcz*ay*bc2*bcdz*bx*bz*cz - 
     -    abcz*bcdx*bcdz*bcz*bx*bz*cz - aby*bc2*bcdx*bcdy*by*bz*cz + 
     -    abx*bc2*bcdy**2*by*bz*cz + 2*abcz*ax*bc2*bcdz*by*bz*cz - 
     -    abcz*bcdy*bcdz*bcz*by*bz*cz + aby*az*bc2*bcdx**2*cx*cz - 
     -    abx*az*bc2*bcdx*bcdy*cx*cz + aby*ax*bc2*bcdx*bcdz*cx*cz - 
     -    abx*ax*bc2*bcdy*bcdz*cx*cz + abcz*az*bcdx*bcdz*bcz*cx*cz + 
     -    abcz*ax*bcdz**2*bcz*cx*cz - aby*bc2*bcdx*bcdz*bx*cx*cz + 
     -    abx*bc2*bcdy*bcdz*bx*cx*cz - abcz*bcdz**2*bcz*bx*cx*cz + 
     -    abcz*az*bc2*bcdz*by*cx*cz - aby*bc2*bcdx**2*bz*cx*cz + 
     -    abx*bc2*bcdx*bcdy*bz*cx*cz - abcz*bcdx*bcdz*bcz*bz*cx*cz - 
     -    abcz*bc2*bcdz*by*bz*cx*cz + aby*az*bc2*bcdx*bcdy*cy*cz - 
     -    abx*az*bc2*bcdy**2*cy*cz + aby*ay*bc2*bcdx*bcdz*cy*cz - 
     -    abx*ay*bc2*bcdy*bcdz*cy*cz + abcz*az*bcdy*bcdz*bcz*cy*cz + 
     -    abcz*ay*bcdz**2*bcz*cy*cz - abcz*az*bc2*bcdz*bx*cy*cz - 
     -    aby*bc2*bcdx*bcdz*by*cy*cz + abx*bc2*bcdy*bcdz*by*cy*cz - 
     -    abcz*bcdz**2*bcz*by*cy*cz - aby*bc2*bcdx*bcdy*bz*cy*cz + 
     -    abx*bc2*bcdy**2*bz*cy*cz - abcz*bcdy*bcdz*bcz*bz*cy*cz + 
     -    abcz*bc2*bcdz*bx*bz*cy*cz - aby*ax*bc2*bcdx**2*cz**2 + 
     -    abx*ax*bc2*bcdx*bcdy*cz**2 - aby*ay*bc2*bcdx*bcdy*cz**2 + 
     -    abx*ay*bc2*bcdy**2*cz**2 - abcz*ax*bcdx*bcdz*bcz*cz**2 - 
     -    abcz*ay*bcdy*bcdz*bcz*cz**2 + aby*bc2*bcdx**2*bx*cz**2 - 
     -    abx*bc2*bcdx*bcdy*bx*cz**2 + abcz*ay*bc2*bcdz*bx*cz**2 + 
     -    abcz*bcdx*bcdz*bcz*bx*cz**2 + aby*bc2*bcdx*bcdy*by*cz**2 - 
     -    abx*bc2*bcdy**2*by*cz**2 - abcz*ax*bc2*bcdz*by*cz**2 + 
     -    abcz*bcdy*bcdz*bcz*by*cz**2 - abcz*ay*bc2*bcdz*bx**2*dx + 
     -    abcz*ax*bc2*bcdz*bx*by*dx + abcz*az*bc2*bcdz*by*bz*dx - 
     -    abcz*ay*bc2*bcdz*bz**2*dx + 2*abcz*ay*bc2*bcdz*bx*cx*dx - 
     -    abcz*ax*bc2*bcdz*by*cx*dx - abcz*bc2*bcdz*bx*by*cx*dx - 
     -    abcz*ay*bc2*bcdz*cx**2*dx + abcz*bc2*bcdz*by*cx**2*dx - 
     -    abcz*ax*bc2*bcdz*bx*cy*dx + abcz*bc2*bcdz*bx**2*cy*dx - 
     -    abcz*az*bc2*bcdz*bz*cy*dx + abcz*bc2*bcdz*bz**2*cy*dx + 
     -    abcz*ax*bc2*bcdz*cx*cy*dx - abcz*bc2*bcdz*bx*cx*cy*dx - 
     -    abcz*az*bc2*bcdz*by*cz*dx + 2*abcz*ay*bc2*bcdz*bz*cz*dx - 
     -    abcz*bc2*bcdz*by*bz*cz*dx + abcz*az*bc2*bcdz*cy*cz*dx - 
     -    abcz*bc2*bcdz*bz*cy*cz*dx - abcz*ay*bc2*bcdz*cz**2*dx + 
     -    abcz*bc2*bcdz*by*cz**2*dx - abcz*ay*bc2*bcdz*bx*by*dy + 
     -    abcz*ax*bc2*bcdz*by**2*dy - abcz*az*bc2*bcdz*bx*bz*dy + 
     -    abcz*ax*bc2*bcdz*bz**2*dy + abcz*ay*bc2*bcdz*by*cx*dy - 
     -    abcz*bc2*bcdz*by**2*cx*dy + abcz*az*bc2*bcdz*bz*cx*dy - 
     -    abcz*bc2*bcdz*bz**2*cx*dy + abcz*ay*bc2*bcdz*bx*cy*dy - 
     -    2*abcz*ax*bc2*bcdz*by*cy*dy + abcz*bc2*bcdz*bx*by*cy*dy - 
     -    abcz*ay*bc2*bcdz*cx*cy*dy + abcz*bc2*bcdz*by*cx*cy*dy + 
     -    abcz*ax*bc2*bcdz*cy**2*dy - abcz*bc2*bcdz*bx*cy**2*dy + 
     -    abcz*az*bc2*bcdz*bx*cz*dy - 2*abcz*ax*bc2*bcdz*bz*cz*dy + 
     -    abcz*bc2*bcdz*bx*bz*cz*dy - abcz*az*bc2*bcdz*cx*cz*dy + 
     -    abcz*bc2*bcdz*bz*cx*cz*dy + abcz*ax*bc2*bcdz*cz**2*dy - 
     -    abcz*bc2*bcdz*bx*cz**2*dy - 
     -    abcx*(-(bc2*bcdx*bcdz*bx**2) + ay*bcdx*bcdy*bcz*bx**2 - 
     -       ay*bc2*bcdx*bx**3 + ay*bc2*bcdx*bcdz*by - 
     -       ay*bcdx**2*bcz*bx*by - ay*bc2*bcdy*bx**2*by - 
     -       bc2*bcdx*bcdz*by**2 - 2*ay*bc2*bcdx*bcdy*bz + 
     -       bc2*bcdx**2*bx*bz + bc2*bcdx*bcdy*by*bz - 
     -       ay*bcdx*bcdz*bcz*by*bz + ay*bc2*bcdz*by**2*bz + 
     -       ay*bcdx*bcdy*bcz*bz**2 - ay*bc2*bcdx*bx*bz**2 - 
     -       ay*bc2*bcdy*by*bz**2 + bc2*bcdx*bcdz*bx*cx - 
     -       2*ay*bcdx*bcdy*bcz*bx*cx + 2*ay*bc2*bcdx*bx**2*cx + 
     -       ay*bcdx**2*bcz*by*cx + 2*ay*bc2*bcdy*bx*by*cx + 
     -       bcdx*bcdy*bcz*bx*by*cx - bc2*bcdx*bx**2*by*cx - 
     -       bcdx**2*bcz*by**2*cx - bc2*bcdy*bx*by**2*cx + 
     -       bc2*bcdx**2*bz*cx + bcdx*bcdz*bcz*bx*bz*cx - 
     -       bc2*bcdz*bx*by*bz*cx - bcdx**2*bcz*bz**2*cx + 
     -       ay*bcdx*bcdy*bcz*cx**2 - ay*bc2*bcdx*bx*cx**2 - 
     -       ay*bc2*bcdy*by*cx**2 - bcdx*bcdy*bcz*by*cx**2 + 
     -       bc2*bcdx*bx*by*cx**2 + bc2*bcdy*by**2*cx**2 - 
     -       bcdx*bcdz*bcz*bz*cx**2 + bc2*bcdz*by*bz*cx**2 - 
     -       ay*bc2*bcdx*bcdz*cy + ay*bcdx**2*bcz*bx*cy - 
     -       bcdx*bcdy*bcz*bx**2*cy + bc2*bcdx*bx**3*cy + 
     -       bc2*bcdx*bcdz*by*cy + bcdx**2*bcz*bx*by*cy + 
     -       bc2*bcdy*bx**2*by*cy + bc2*bcdx*bcdy*bz*cy + 
     -       ay*bcdx*bcdz*bcz*bz*cy - ay*bc2*bcdz*by*bz*cy + 
     -       bcdx*bcdz*bcz*by*bz*cy - bc2*bcdz*by**2*bz*cy - 
     -       bcdx*bcdy*bcz*bz**2*cy + bc2*bcdx*bx*bz**2*cy + 
     -       bc2*bcdy*by*bz**2*cy - ay*bcdx**2*bcz*cx*cy + 
     -       bcdx*bcdy*bcz*bx*cx*cy - bc2*bcdx*bx**2*cx*cy + 
     -       bcdx**2*bcz*by*cx*cy - bc2*bcdy*bx*by*cx*cy - 
     -       bcdx**2*bcz*bx*cy**2 - bcdx*bcdz*bcz*bz*cy**2 + 
     -       bc2*bcdz*by*bz*cy**2 + 2*ay*bc2*bcdx*bcdy*cz - 
     -       2*bc2*bcdx**2*bx*cz - bcdx*bcdz*bcz*bx**2*cz - 
     -       2*bc2*bcdx*bcdy*by*cz + ay*bcdx*bcdz*bcz*by*cz + 
     -       bc2*bcdz*bx**2*by*cz - ay*bc2*bcdz*by**2*cz - 
     -       bcdx*bcdz*bcz*by**2*cz + bc2*bcdz*by**3*cz - 
     -       2*ay*bcdx*bcdy*bcz*bz*cz + 2*ay*bc2*bcdx*bx*bz*cz + 
     -       bcdx**2*bcz*bx*bz*cz + 2*ay*bc2*bcdy*by*bz*cz + 
     -       bcdx*bcdy*bcz*by*bz*cz - bc2*bcdx*bx*by*bz*cz - 
     -       bc2*bcdy*by**2*bz*cz + bcdx*bcdz*bcz*bx*cx*cz - 
     -       bc2*bcdz*bx*by*cx*cz + bcdx**2*bcz*bz*cx*cz - 
     -       ay*bcdx*bcdz*bcz*cy*cz + ay*bc2*bcdz*by*cy*cz + 
     -       bcdx*bcdz*bcz*by*cy*cz - bc2*bcdz*by**2*cy*cz + 
     -       bcdx*bcdy*bcz*bz*cy*cz - bc2*bcdx*bx*bz*cy*cz - 
     -       bc2*bcdy*by*bz*cy*cz + ay*bcdx*bcdy*bcz*cz**2 - 
     -       ay*bc2*bcdx*bx*cz**2 - bcdx**2*bcz*bx*cz**2 - 
     -       ay*bc2*bcdy*by*cz**2 - bcdx*bcdy*bcz*by*cz**2 + 
     -       bc2*bcdx*bx*by*cz**2 + bc2*bcdy*by**2*cz**2 + 
     -       ay*bc2*bcdx*bx**2*dx + ay*bc2*bcdx*bz**2*dx - 
     -       2*ay*bc2*bcdx*bx*cx*dx + bc2*bcdx*bx*by*cx*dx + 
     -       ay*bc2*bcdx*cx**2*dx - bc2*bcdx*by*cx**2*dx - 
     -       bc2*bcdx*bx**2*cy*dx - bc2*bcdx*bz**2*cy*dx + 
     -       bc2*bcdx*bx*cx*cy*dx - 2*ay*bc2*bcdx*bz*cz*dx + 
     -       bc2*bcdx*by*bz*cz*dx + bc2*bcdx*bz*cy*cz*dx + 
     -       ay*bc2*bcdx*cz**2*dx - bc2*bcdx*by*cz**2*dx + 
     -       ax*(bcdx**2*bcz*(-2*bc2 + bcy**2 + bcz**2) + 
     -          abx*(-(bcdx*bcz*(bcdy*bcy + bcdz*bcz)) + 
     -             bc2*(bcdx*(bcdz + bcy*bx - bcy*dx) + 
     -                (bcdy*bcy + bcdz*bcz)*(by - dy)))) + 
     -       az*(bcdx*bcz*(bcdz*(bcx**2 + bcy**2) - 
     -             (abx*bcdx + bcdy*bcy)*bcz) + 
     -          bc2*(abx*bcdx**2 + bcdx*bcy*(bcdy + bcz*(bx - dx)) - 
     -             (bcdz*(bcx**2 + bcy**2) - bcdy*bcy*bcz)*(by - dy)))
     -        + ay*bc2*bcdy*bx**2*dy - ay*bc2*bcdz*by*bz*dy + 
     -       ay*bc2*bcdy*bz**2*dy - 2*ay*bc2*bcdy*bx*cx*dy + 
     -       bc2*bcdy*bx*by*cx*dy + bc2*bcdz*bx*bz*cx*dy + 
     -       ay*bc2*bcdy*cx**2*dy - bc2*bcdy*by*cx**2*dy - 
     -       bc2*bcdz*bz*cx**2*dy - bc2*bcdy*bx**2*cy*dy + 
     -       ay*bc2*bcdz*bz*cy*dy + bc2*bcdz*by*bz*cy*dy - 
     -       bc2*bcdy*bz**2*cy*dy + bc2*bcdy*bx*cx*cy*dy - 
     -       bc2*bcdz*bz*cy**2*dy - bc2*bcdz*bx**2*cz*dy + 
     -       ay*bc2*bcdz*by*cz*dy - bc2*bcdz*by**2*cz*dy - 
     -       2*ay*bc2*bcdy*bz*cz*dy + bc2*bcdy*by*bz*cz*dy + 
     -       bc2*bcdz*bx*cx*cz*dy - ay*bc2*bcdz*cy*cz*dy + 
     -       bc2*bcdz*by*cy*cz*dy + bc2*bcdy*bz*cy*cz*dy + 
     -       ay*bc2*bcdy*cz**2*dy - bc2*bcdy*by*cz**2*dy) - 
     -    abcy*(-(bc2*bcdy*bcdz*bx**2) + ay*bcdy**2*bcz*bx**2 + 
     -       ay*bc2*bcdy*bcdz*by - ay*bcdx*bcdy*bcz*bx*by - 
     -       ay*bc2*bcdx*bx**2*by - bc2*bcdy*bcdz*by**2 - 
     -       ay*bc2*bcdy*bx*by**2 - 2*ay*bc2*bcdy**2*bz + 
     -       bc2*bcdx*bcdy*bx*bz + bc2*bcdy**2*by*bz - 
     -       ay*bcdy*bcdz*bcz*by*bz - ay*bc2*bcdz*bx*by*bz + 
     -       ay*bcdy**2*bcz*bz**2 + bc2*bcdy*bcdz*bx*cx - 
     -       2*ay*bcdy**2*bcz*bx*cx + ay*bcdx*bcdy*bcz*by*cx + 
     -       ay*bc2*bcdx*bx*by*cx + bcdy**2*bcz*bx*by*cx + 
     -       ay*bc2*bcdy*by**2*cx - bcdx*bcdy*bcz*by**2*cx - 
     -       bc2*bcdx*bx*by**2*cx - bc2*bcdy*by**3*cx + 
     -       bc2*bcdx*bcdy*bz*cx + bcdy*bcdz*bcz*bx*bz*cx + 
     -       bc2*bcdz*bx**2*bz*cx - bcdx*bcdy*bcz*bz**2*cx - 
     -       bc2*bcdx*bx*bz**2*cx - bc2*bcdy*by*bz**2*cx + 
     -       ay*bcdy**2*bcz*cx**2 - bcdy**2*bcz*by*cx**2 - 
     -       bcdy*bcdz*bcz*bz*cx**2 - bc2*bcdz*bx*bz*cx**2 - 
     -       ay*bc2*bcdy*bcdz*cy + ay*bcdx*bcdy*bcz*bx*cy + 
     -       ay*bc2*bcdx*bx**2*cy - bcdy**2*bcz*bx**2*cy + 
     -       bc2*bcdy*bcdz*by*cy + ay*bc2*bcdy*bx*by*cy + 
     -       bcdx*bcdy*bcz*bx*by*cy + bc2*bcdx*bx**2*by*cy + 
     -       bc2*bcdy*bx*by**2*cy + bc2*bcdy**2*bz*cy + 
     -       ay*bcdy*bcdz*bcz*bz*cy + ay*bc2*bcdz*bx*bz*cy + 
     -       bcdy*bcdz*bcz*by*bz*cy + bc2*bcdz*bx*by*bz*cy - 
     -       bcdy**2*bcz*bz**2*cy - ay*bcdx*bcdy*bcz*cx*cy - 
     -       ay*bc2*bcdx*bx*cx*cy + bcdy**2*bcz*bx*cx*cy - 
     -       ay*bc2*bcdy*by*cx*cy + bcdx*bcdy*bcz*by*cx*cy + 
     -       bc2*bcdx*bx*by*cx*cy + bc2*bcdy*by**2*cx*cy - 
     -       bcdx*bcdy*bcz*bx*cy**2 - bc2*bcdx*bx**2*cy**2 - 
     -       bc2*bcdy*bx*by*cy**2 - bcdy*bcdz*bcz*bz*cy**2 - 
     -       bc2*bcdz*bx*bz*cy**2 + 2*ay*bc2*bcdy**2*cz - 
     -       2*bc2*bcdx*bcdy*bx*cz - bcdy*bcdz*bcz*bx**2*cz - 
     -       bc2*bcdz*bx**3*cz - 2*bc2*bcdy**2*by*cz + 
     -       ay*bcdy*bcdz*bcz*by*cz + ay*bc2*bcdz*bx*by*cz - 
     -       bcdy*bcdz*bcz*by**2*cz - bc2*bcdz*bx*by**2*cz - 
     -       2*ay*bcdy**2*bcz*bz*cz + bcdx*bcdy*bcz*bx*bz*cz + 
     -       bc2*bcdx*bx**2*bz*cz + bcdy**2*bcz*by*bz*cz + 
     -       bc2*bcdy*bx*by*bz*cz + bcdy*bcdz*bcz*bx*cx*cz + 
     -       bc2*bcdz*bx**2*cx*cz + bcdx*bcdy*bcz*bz*cx*cz + 
     -       bc2*bcdx*bx*bz*cx*cz + bc2*bcdy*by*bz*cx*cz - 
     -       ay*bcdy*bcdz*bcz*cy*cz - ay*bc2*bcdz*bx*cy*cz + 
     -       bcdy*bcdz*bcz*by*cy*cz + bc2*bcdz*bx*by*cy*cz + 
     -       bcdy**2*bcz*bz*cy*cz + ay*bcdy**2*bcz*cz**2 - 
     -       bcdx*bcdy*bcz*bx*cz**2 - bc2*bcdx*bx**2*cz**2 - 
     -       bcdy**2*bcz*by*cz**2 - bc2*bcdy*bx*by*cz**2 + 
     -       ay*bc2*bcdx*bx*by*dx + ay*bc2*bcdz*by*bz*dx - 
     -       ay*bc2*bcdx*by*cx*dx + bc2*bcdx*by**2*cx*dx - 
     -       bc2*bcdz*bx*bz*cx*dx + bc2*bcdx*bz**2*cx*dx + 
     -       bc2*bcdz*bz*cx**2*dx - ay*bc2*bcdx*bx*cy*dx - 
     -       bc2*bcdx*bx*by*cy*dx - ay*bc2*bcdz*bz*cy*dx - 
     -       bc2*bcdz*by*bz*cy*dx + ay*bc2*bcdx*cx*cy*dx - 
     -       bc2*bcdx*by*cx*cy*dx + bc2*bcdx*bx*cy**2*dx + 
     -       bc2*bcdz*bz*cy**2*dx + bc2*bcdz*bx**2*cz*dx - 
     -       ay*bc2*bcdz*by*cz*dx + bc2*bcdz*by**2*cz*dx - 
     -       bc2*bcdx*bx*bz*cz*dx - bc2*bcdz*bx*cx*cz*dx - 
     -       bc2*bcdx*bz*cx*cz*dx + ay*bc2*bcdz*cy*cz*dx - 
     -       bc2*bcdz*by*cy*cz*dx + bc2*bcdx*bx*cz**2*dx + 
     -       ax*(bcdx*(bcdy*(bcy**2*bcz + bcz**3 - 2*bc2*bz + 
     -                2*bc2*cz) + bc2*(bcy**2 + bcz**2)*(bx - dx)) + 
     -          abx*(-(bcdy*bcz*(bcdy*bcy + bcdz*bcz)) + 
     -             bc2*bcdz*(bcdy - bcz*bx + bcz*dx)) + 
     -          bc2*bcdy*(by**2 + bz**2 - 2*by*cy + cy**2 - 2*bz*cz + 
     -             cz**2)*(by - dy)) + ay*bc2*bcdy*bx*by*dy - 
     -       ay*bc2*bcdy*by*cx*dy + bc2*bcdy*by**2*cx*dy + 
     -       bc2*bcdy*bz**2*cx*dy - ay*bc2*bcdy*bx*cy*dy - 
     -       bc2*bcdy*bx*by*cy*dy + ay*bc2*bcdy*cx*cy*dy - 
     -       bc2*bcdy*by*cx*cy*dy + bc2*bcdy*bx*cy**2*dy - 
     -       bc2*bcdy*bx*bz*cz*dy - bc2*bcdy*bz*cx*cz*dy + 
     -       bc2*bcdy*bx*cz**2*dy + 
     -       az*(bcdy*bcz*(bcdz*(bcx**2 + bcy**2) - 
     -             (abx*bcdx + bcdy*bcy)*bcz) + 
     -          bc2*(bcdz*(bcx**2 + bcy**2)*(bx - dx) + 
     -             abx*bcdx*(bcdy + bcz*(-bx + dx)) + 
     -             bcdy*(bcdy*bcy + abx*bcz*(-by + dy))))))/
     -  (bc2**1.5*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**2/
     -       bc2))


       ddx=(-(abcy*(-(ay*bcy*bcz*(abx*bcdx + bcdy*bcy + bcdz*bcz)) + 
     -         az*(bcdy*bcy*(bcx**2 + bcy**2) + 
     -            bcz*(bcdz*(bcx**2 + bcy**2) - abx*bcdx*bcz)) + 
     -         ax*bcdx*bcz*by**2 - ax*bcdy*bcy*bx*bz - 
     -         ax*bcdz*bcz*bx*bz + ax*bcdx*bcz*bz**2 - 
     -         bcdx*bcz*by**2*cx + ax*bcdy*bcy*bz*cx + 
     -         ax*bcdz*bcz*bz*cx + bcdy*bcy*bx*bz*cx + 
     -         bcdz*bcz*bx*bz*cx - bcdx*bcz*bz**2*cx - 
     -         bcdy*bcy*bz*cx**2 - bcdz*bcz*bz*cx**2 - 
     -         2*ax*bcdx*bcz*by*cy + bcdx*bcz*bx*by*cy + 
     -         bcdy*bcy*by*bz*cy + bcdz*bcz*by*bz*cy + 
     -         bcdx*bcz*by*cx*cy + ax*bcdx*bcz*cy**2 - 
     -         bcdx*bcz*bx*cy**2 - bcdy*bcy*bz*cy**2 - 
     -         bcdz*bcz*bz*cy**2 + ax*bcdy*bcy*bx*cz + 
     -         ax*bcdz*bcz*bx*cz - bcdy*bcy*bx**2*cz - 
     -         bcdz*bcz*bx**2*cz - bcdy*bcy*by**2*cz - 
     -         bcdz*bcz*by**2*cz - 2*ax*bcdx*bcz*bz*cz + 
     -         bcdx*bcz*bx*bz*cz - ax*bcdy*bcy*cx*cz - 
     -         ax*bcdz*bcz*cx*cz + bcdy*bcy*bx*cx*cz + 
     -         bcdz*bcz*bx*cx*cz + bcdx*bcz*bz*cx*cz + 
     -         bcdy*bcy*by*cy*cz + bcdz*bcz*by*cy*cz + 
     -         ax*bcdx*bcz*cz**2 - bcdx*bcz*bx*cz**2)) + 
     -    abcx*bcdx*(-(az*bcy*(bcx**2 + bcy**2 + bcz**2)) + 
     -       ay*bcz*(bcx**2 + bcy**2 + bcz**2) - ax*bcz*bx*by + 
     -       ax*bcy*bx*bz + ax*bcz*by*cx + bcz*bx*by*cx - 
     -       ax*bcy*bz*cx - bcy*bx*bz*cx - bcz*by*cx**2 + 
     -       bcy*bz*cx**2 + ax*bcz*bx*cy - bcz*bx**2*cy - 
     -       bcy*by*bz*cy - bcz*bz**2*cy - ax*bcz*cx*cy + 
     -       bcz*bx*cx*cy + bcy*bz*cy**2 - ax*bcy*bx*cz + 
     -       bcy*bx**2*cz + bcy*by**2*cz + bcz*by*bz*cz + 
     -       ax*bcy*cx*cz - bcy*bx*cx*cz - bcy*by*cy*cz + 
     -       bcz*bz*cy*cz - bcz*by*cz**2) + 
     -    abcz*(ay*(-(abx*bcdx*bcy**2) + 
     -          bcdy*bcy*(bcx**2 + bcz**2) + 
     -          bcdz*bcz*(bcx**2 + bcz**2)) - az*bcdx*bcy*bx*bz - 
     -       az*bcdy*bcy*by*bz - az*bcdz*bcz*by*bz + 
     -       bcdy*bcy*bx*by*cx + bcdz*bcz*bx*by*cx - 
     -       bcdx*bcy*by**2*cx + az*bcdx*bcy*bz*cx - 
     -       bcdx*bcy*bz**2*cx - bcdy*bcy*by*cx**2 - 
     -       bcdz*bcz*by*cx**2 - bcdy*bcy*bx**2*cy - 
     -       bcdz*bcz*bx**2*cy + bcdx*bcy*bx*by*cy + 
     -       az*bcdy*bcy*bz*cy + az*bcdz*bcz*bz*cy - 
     -       bcdy*bcy*bz**2*cy - bcdz*bcz*bz**2*cy + 
     -       bcdy*bcy*bx*cx*cy + bcdz*bcz*bx*cx*cy + 
     -       bcdx*bcy*by*cx*cy - bcdx*bcy*bx*cy**2 + 
     -       az*bcdx*bcy*bx*cz + az*bcdy*bcy*by*cz + 
     -       az*bcdz*bcz*by*cz + bcdx*bcy*bx*bz*cz + 
     -       bcdy*bcy*by*bz*cz + bcdz*bcz*by*bz*cz - 
     -       az*bcdx*bcy*cx*cz + bcdx*bcy*bz*cx*cz - 
     -       az*bcdy*bcy*cy*cz - az*bcdz*bcz*cy*cz + 
     -       bcdy*bcy*bz*cy*cz + bcdz*bcz*bz*cy*cz - 
     -       bcdx*bcy*bx*cz**2 - bcdy*bcy*by*cz**2 - 
     -       bcdz*bcz*by*cz**2 + 
     -       ax*(-(bcdy*bcy*(bx - cx)*(by - cy)) - 
     -          bcdz*bcz*(bx - cx)*(by - cy) + 
     -          bcdx*bcy*(by**2 + bz**2 - 2*by*cy + cy**2 - 
     -             2*bz*cz + cz**2))))/
     -  (Sqrt(bc2)*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**
     -        2/bc2))


       ddy=(-(abcz*(ay*(-(abx*bcy*(bcdx*bcx + bcdz*bcz)) + 
     -            bcdy*bcx*(bcx**2 + bcz**2)) + 
     -         ax*(-(abx*bcdy*bcx*bcy) + 
     -            (bcdx*bcx + bcdz*bcz)*(bcy**2 + bcz**2)) - 
     -         az*bcdx*bcx*bx*bz - az*bcdz*bcz*bx*bz - 
     -         az*bcdy*bcx*by*bz + bcdy*bcx*bx*by*cx - 
     -         bcdx*bcx*by**2*cx - bcdz*bcz*by**2*cx + 
     -         az*bcdx*bcx*bz*cx + az*bcdz*bcz*bz*cx - 
     -         bcdx*bcx*bz**2*cx - bcdz*bcz*bz**2*cx - 
     -         bcdy*bcx*by*cx**2 - bcdy*bcx*bx**2*cy + 
     -         bcdx*bcx*bx*by*cy + bcdz*bcz*bx*by*cy + 
     -         az*bcdy*bcx*bz*cy - bcdy*bcx*bz**2*cy + 
     -         bcdy*bcx*bx*cx*cy + bcdx*bcx*by*cx*cy + 
     -         bcdz*bcz*by*cx*cy - bcdx*bcx*bx*cy**2 - 
     -         bcdz*bcz*bx*cy**2 + az*bcdx*bcx*bx*cz + 
     -         az*bcdz*bcz*bx*cz + az*bcdy*bcx*by*cz + 
     -         bcdx*bcx*bx*bz*cz + bcdz*bcz*bx*bz*cz + 
     -         bcdy*bcx*by*bz*cz - az*bcdx*bcx*cx*cz - 
     -         az*bcdz*bcz*cx*cz + bcdx*bcx*bz*cx*cz + 
     -         bcdz*bcz*bz*cx*cz - az*bcdy*bcx*cy*cz + 
     -         bcdy*bcx*bz*cy*cz - bcdx*bcx*bx*cz**2 - 
     -         bcdz*bcz*bx*cz**2 - bcdy*bcx*by*cz**2)) + 
     -    abcx*(az*(bcdx*bcx*(bcx**2 + bcy**2) + 
     -          bcz*(bcdz*(bcx**2 + bcy**2) - bcdy*bcy*bcz)) + 
     -       ay*bcz*(-(bcy*(bcdx*bcx + bcdz*bcz)) + 
     -          bcdy*(bcx**2 + bcz**2)) - ax*bcdy*bcz*bx*by - 
     -       ax*bcdx*bcx*bx*bz - ax*bcdz*bcz*bx*bz + 
     -       ax*bcdy*bcz*by*cx + bcdy*bcz*bx*by*cx + 
     -       ax*bcdx*bcx*bz*cx + ax*bcdz*bcz*bz*cx + 
     -       bcdx*bcx*bx*bz*cx + bcdz*bcz*bx*bz*cx - 
     -       bcdy*bcz*by*cx**2 - bcdx*bcx*bz*cx**2 - 
     -       bcdz*bcz*bz*cx**2 + ax*bcdy*bcz*bx*cy - 
     -       bcdy*bcz*bx**2*cy + bcdx*bcx*by*bz*cy + 
     -       bcdz*bcz*by*bz*cy - bcdy*bcz*bz**2*cy - 
     -       ax*bcdy*bcz*cx*cy + bcdy*bcz*bx*cx*cy - 
     -       bcdx*bcx*bz*cy**2 - bcdz*bcz*bz*cy**2 + 
     -       ax*bcdx*bcx*bx*cz + ax*bcdz*bcz*bx*cz - 
     -       bcdx*bcx*bx**2*cz - bcdz*bcz*bx**2*cz - 
     -       bcdx*bcx*by**2*cz - bcdz*bcz*by**2*cz + 
     -       bcdy*bcz*by*bz*cz - ax*bcdx*bcx*cx*cz - 
     -       ax*bcdz*bcz*cx*cz + bcdx*bcx*bx*cx*cz + 
     -       bcdz*bcz*bx*cx*cz + bcdx*bcx*by*cy*cz + 
     -       bcdz*bcz*by*cy*cz + bcdy*bcz*bz*cy*cz - 
     -       bcdy*bcz*by*cz**2) - 
     -    abcy*bcdy*(-(az*(bcx**3 + bcx*bcy**2 + abx*bcz**2)) + 
     -       ax*bcz*by**2 + ax*bcx*bx*bz + ax*bcz*bz**2 - 
     -       bcz*by**2*cx - ax*bcx*bz*cx - bcx*bx*bz*cx - 
     -       bcz*bz**2*cx + bcx*bz*cx**2 - 2*ax*bcz*by*cy + 
     -       bcz*bx*by*cy - bcx*by*bz*cy + bcz*by*cx*cy + 
     -       ax*bcz*cy**2 - bcz*bx*cy**2 + bcx*bz*cy**2 - 
     -       ax*bcx*bx*cz + bcx*bx**2*cz + bcx*by**2*cz - 
     -       2*ax*bcz*bz*cz + bcz*bx*bz*cz + ax*bcx*cx*cz - 
     -       bcx*bx*cx*cz + bcz*bz*cx*cz - bcx*by*cy*cz + 
     -       ax*bcz*cz**2 - bcz*bx*cz**2 + 
     -       ay*bcy*(-(bcz*bx) + bcx*bz + bcz*cx - bcx*cz)))/
     -  (Sqrt(bc2)*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**
     -        2/bc2))


       ddz= (abcy*(-(ay*bcy*(abx*bcdx*bcx + abx*bcdy*bcy + 
     -            bcdz*bcx*bcz)) + 
     -       az*(bcdz*bcx*(bcx**2 + bcy**2) - 
     -          abx*(bcdx*bcx + bcdy*bcy)*bcz) + ax*bcdx*bcx*by**2 + 
     -       ax*bcdy*bcy*by**2 - ax*bcdz*bcx*bx*bz + 
     -       ax*bcdx*bcx*bz**2 + ax*bcdy*bcy*bz**2 - 
     -       bcdx*bcx*by**2*cx - bcdy*bcy*by**2*cx + 
     -       ax*bcdz*bcx*bz*cx + bcdz*bcx*bx*bz*cx - 
     -       bcdx*bcx*bz**2*cx - bcdy*bcy*bz**2*cx - 
     -       bcdz*bcx*bz*cx**2 - 2*ax*bcdx*bcx*by*cy - 
     -       2*ax*bcdy*bcy*by*cy + bcdx*bcx*bx*by*cy + 
     -       bcdy*bcy*bx*by*cy + bcdz*bcx*by*bz*cy + 
     -       bcdx*bcx*by*cx*cy + bcdy*bcy*by*cx*cy + 
     -       ax*bcdx*bcx*cy**2 + ax*bcdy*bcy*cy**2 - 
     -       bcdx*bcx*bx*cy**2 - bcdy*bcy*bx*cy**2 - 
     -       bcdz*bcx*bz*cy**2 + ax*bcdz*bcx*bx*cz - 
     -       bcdz*bcx*bx**2*cz - bcdz*bcx*by**2*cz - 
     -       2*ax*bcdx*bcx*bz*cz - 2*ax*bcdy*bcy*bz*cz + 
     -       bcdx*bcx*bx*bz*cz + bcdy*bcy*bx*bz*cz - 
     -       ax*bcdz*bcx*cx*cz + bcdz*bcx*bx*cx*cz + 
     -       bcdx*bcx*bz*cx*cz + bcdy*bcy*bz*cx*cz + 
     -       bcdz*bcx*by*cy*cz + ax*bcdx*bcx*cz**2 + 
     -       ax*bcdy*bcy*cz**2 - bcdx*bcx*bx*cz**2 - 
     -       bcdy*bcy*bx*cz**2) + 
     -    abcz*bcdz*(ax*bcy*(abx*bcx + bcy**2 + bcz**2) - 
     -       ay*(bcx**3 + abx*bcy**2 + bcx*bcz**2) - az*bcy*bx*bz + 
     -       az*bcx*by*bz - bcx*bx*by*cx - bcy*by**2*cx + 
     -       az*bcy*bz*cx - bcy*bz**2*cx + bcx*by*cx**2 + 
     -       bcx*bx**2*cy + bcy*bx*by*cy - az*bcx*bz*cy + 
     -       bcx*bz**2*cy - bcx*bx*cx*cy + bcy*by*cx*cy - 
     -       bcy*bx*cy**2 + az*bcy*bx*cz - az*bcx*by*cz + 
     -       bcy*bx*bz*cz - bcx*by*bz*cz - az*bcy*cx*cz + 
     -       bcy*bz*cx*cz + az*bcx*cy*cz - bcx*bz*cy*cz - 
     -       bcy*bx*cz**2 + bcx*by*cz**2) - 
     -    abcx*(az*bcy*(bcdz*(bcx**2 + bcy**2) - 
     -          (bcdx*bcx + bcdy*bcy)*bcz) + 
     -       ay*(bcdx*bcx*(bcx**2 + bcz**2) + 
     -          bcy*(-(bcdz*bcy*bcz) + bcdy*(bcx**2 + bcz**2))) - 
     -       ax*bcdx*bcx*bx*by - ax*bcdy*bcy*bx*by - 
     -       ax*bcdz*bcy*bx*bz + ax*bcdx*bcx*by*cx + 
     -       ax*bcdy*bcy*by*cx + bcdx*bcx*bx*by*cx + 
     -       bcdy*bcy*bx*by*cx + ax*bcdz*bcy*bz*cx + 
     -       bcdz*bcy*bx*bz*cx - bcdx*bcx*by*cx**2 - 
     -       bcdy*bcy*by*cx**2 - bcdz*bcy*bz*cx**2 + 
     -       ax*bcdx*bcx*bx*cy + ax*bcdy*bcy*bx*cy - 
     -       bcdx*bcx*bx**2*cy - bcdy*bcy*bx**2*cy + 
     -       bcdz*bcy*by*bz*cy - bcdx*bcx*bz**2*cy - 
     -       bcdy*bcy*bz**2*cy - ax*bcdx*bcx*cx*cy - 
     -       ax*bcdy*bcy*cx*cy + bcdx*bcx*bx*cx*cy + 
     -       bcdy*bcy*bx*cx*cy - bcdz*bcy*bz*cy**2 + 
     -       ax*bcdz*bcy*bx*cz - bcdz*bcy*bx**2*cz - 
     -       bcdz*bcy*by**2*cz + bcdx*bcx*by*bz*cz + 
     -       bcdy*bcy*by*bz*cz - ax*bcdz*bcy*cx*cz + 
     -       bcdz*bcy*bx*cx*cz + bcdz*bcy*by*cy*cz + 
     -       bcdx*bcx*bz*cy*cz + bcdy*bcy*bz*cy*cz - 
     -       bcdx*bcx*by*cz**2 - bcdy*bcy*by*cz**2))/
     -  (Sqrt(bc2)*((abcx*bcdx + abcy*bcdy + abcz*bcdz)**2 + 
     -      (-(ax*bcdy*bx*by) + ax*bcdx*by**2 - ax*bcdz*bx*bz + 
     -          ax*bcdx*bz**2 + ax*bcdy*by*cx + bcdy*bx*by*cx - 
     -          bcdx*by**2*cx + ax*bcdz*bz*cx + bcdz*bx*bz*cx - 
     -          bcdx*bz**2*cx - bcdy*by*cx**2 - bcdz*bz*cx**2 + 
     -          ax*bcdy*bx*cy - bcdy*bx**2*cy - 2*ax*bcdx*by*cy + 
     -          bcdx*bx*by*cy + bcdz*by*bz*cy - bcdy*bz**2*cy - 
     -          ax*bcdy*cx*cy + bcdy*bx*cx*cy + bcdx*by*cx*cy + 
     -          ax*bcdx*cy**2 - bcdx*bx*cy**2 - bcdz*bz*cy**2 + 
     -          az*(bcdz*(bcx**2 + bcy**2) + 
     -             bcz*(-(bcdx*bx) - bcdy*by + bcdx*cx + bcdy*cy)) + 
     -          ax*bcdz*bx*cz - bcdz*bx**2*cz - bcdz*by**2*cz - 
     -          2*ax*bcdx*bz*cz + bcdx*bx*bz*cz + bcdy*by*bz*cz - 
     -          ax*bcdz*cx*cz + bcdz*bx*cx*cz + bcdx*bz*cx*cz + 
     -          bcdz*by*cy*cz + bcdy*bz*cy*cz + ax*bcdx*cz**2 - 
     -          bcdx*bx*cz**2 - bcdy*by*cz**2 + 
     -          ay*(bcdy*(bcx**2 + bcz**2) + 
     -             bcy*(-(bcdx*bx) - bcdz*bz + bcdx*cx + bcdz*cz)))**
     -        2/bc2))



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
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz)
        call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
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
