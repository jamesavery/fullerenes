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
          write(*,*)"entering extwu"
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
      SUBROUTINE DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz)
      implicit real*8 (a-z)
      dist_ab_inv=1/dsqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)
      aux_1=ax-bx
      aux_2=ay-by
      aux_3=az-bz
      dax=aux_1*dist_ab_inv
      day=-dax
      daz=aux_2*dist_ab_inv
      dbx=-day
      dby=aux_3*dist_ab_inv
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
c      aux=aux_ax*aux_bx + aux_ay*aux_by + aux_az*aux_bz
c      angle_abc=dacos(aux/(r1l*r1r))
c      return
c      END


c subroutine dangle takes 9 reals (=3 coordinates) and yields all 9 first derivations of the angle
      SUBROUTINE DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz)
      implicit real*8 (a-z)
c vectors from a to b and b to c
      aux_lx=ax-bx
      aux_ly=ay-by
      aux_lz=az-bz
      aux_rx=bx-cx
      aux_ry=by-cy
      aux_rz=bz-cz
c length of a-b and b-c
      r2L=aux_lx**2 + aux_ly**2 + aux_lz**2
      r2R=aux_rx**2 + aux_ry**2 + aux_rz**2
      r1L=dsqrt(r2l)
      r1R=dsqrt(r2R)
      r3L=r2l*r1l
      r3R=r2r*r1r
c some auxiliary products
      l_dot_r=aux_lx*aux_rx + aux_ly*aux_ry + aux_lz*aux_rz
      aux_1_inv=1/dsqrt(r1L*r1R)
      aux_2_inv=1/dsqrt(r3L*r1R)
      aux_3_inv=1/dsqrt(r1L*r3R)
      den_inv=1/dsqrt(1-(l_dot_r**2)*aux_1_inv**2)
      aux_4=l_dot_r*aux_2_inv
      aux_5=l_dot_r*aux_3_inv
c the derivations
      dax=(-aux_rx*aux_1_inv + aux_lx*aux_4)*den_inv
      day=(-aux_ry*aux_1_inv + aux_ly*aux_4)*den_inv
      daz=(-aux_rz*aux_1_inv + aux_lz*aux_4)*den_inv
      dbx=(aux_rx*aux_5-(aux_lx-aux_rx)*aux_1_inv-aux_lx*aux_4)*den_inv
      dby=(aux_ry*aux_5-(aux_ly-aux_ry)*aux_1_inv-aux_ly*aux_4)*den_inv
      dbz=(aux_rz*aux_5-(aux_lz-aux_rz)*aux_1_inv-aux_lz*aux_4)*den_inv
      dcx=(aux_lx*aux_1_inv - aux_rx*aux_5)*den_inv
      dcy=(aux_ly*aux_1_inv - aux_ry*aux_5)*den_inv
      dcz=(aux_lz*aux_1_inv - aux_rz*aux_5)*den_inv
      return
      END


c subroutine dist takes 12 reals (=4 coordinates) and yields an angel between -\pi and +\pi (in radians)
      SUBROUTINE DIHEDRAL(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2 dihedral_abcd)
      IMPLICIT REAL*8 (a-z)
c normal vectors on abc and bcd
      abc_x=-az*by+ay*bz+az*cy-bz*cy-ay*cz+by*cz
      abc_y= az*bx-ax*bz-az*cx+bz*cx+ax*cz-bx*cz
      abc_z=-ay*bx+ax*by+ay*cx-by*cx-ax*cy+bx*cy
      bcd_x=-bz*cy+by*cz+bz*dy-cz*dy-by*dz+cy*dz
      bcd_y= bz*cx-bx*cz-bz*dx+cz*dx+bx*dz-cx*dz
      bcd_z=-by*cx+bx*cy+by*dx-cy*dx-bx*dy+cx*dy
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
      bc_x=bx-cx
      bc_y=by-cy
      bc_z=bz-cz
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

      AA=-az*by+ay*bz+az*cy-bz*cy-ay*cz+by*cz
      AB= az*bx-ax*bz-az*cx+bz*cx+ax*cz-bx*cz
      AC=-ay*bx+ax*by+ay*cx-by*cx-ax*cy+bx*cy
      AD=-bz*cy+by*cz+bz*dy-cz*dy-by*dz+cy*dz
      AE= bz*cx-bx*cz-bz*dx+cz*dx+bx*dz-cx*dz
      AF=-by*cx+bx*cy+by*dx-cy*dx-bx*dy+cx*dy 
c      write(*,*)"abcdef",aa,ab,ac,ad,ae,af
    
      AG=ax-bx
      AH=ay-by
      AI=az-bz
      AJ=bx-cx
      AK=by-cy
      AL=bz-cz
      AM=cx-dx
      AN=cy-dy
      AO=cz-dz
c      write(*,*)"g-o",AG,AH,AI,AJ,AK,AL,AM,AN,AO
      BB=ax-cx
      BC=ay-cy
      BD=az-cz
      BE=bx-dx
      BF=by-dy
      BG=bz-dz

      BH= 2*AC*BB - 2*AA*BD
      BI=-2*AF*AM + 2*AD*AO
c      write(*,*)"hi",bh,bi

      AP=AA**2 + AB**2 + AC**2
      AQ=AJ**2 + AK**2 + AL**2
      AR=AD**2 + AE**2 + AF**2

      AS=AI**2+AL**2
      AT=AJ**2+AK**2
      AU=AJ**2+AL**2
c      write(*,*)"p-u",ap,aq,ar,as,at,au

      AV=-(AJ*AK*ax) + AU*ay - az*by*bz + bx*by*cx - by*cx**2 -
     2    bx**2*cy + az*bz*cy - bz**2*cy + bx*cx*cy + az*by*cz +
     3    by*bz*cz - az*cy*cz + bz*cy*cz - by*cz**2
      AW=-(AJ*AL*ax) + AT*az - ay*by*bz + bx*bz*cx - bz*cx**2 +
     2    ay*bz*cy + by*bz*cy - bz*cy**2 - bx**2*cz + ay*by*cz -
     3    by**2*cz + bx*cx*cz - ay*cy*cz + by*cy*cz
      BA=AS*ax - AJ*AK*ay - az*bx*bz - by**2*cx + az*bz*cx -
     2    bz**2*cx + bx*by*cy + by*cx*cy - bx*cy**2 + az*bx*cz +
     3    bx*bz*cz - az*cx*cz + bz*cx*cz - bx*cz**2
c      write(*,*)"v-ba",av,aw,ba
      BM=-(AJ*AL*ax) + AJ**2*az + AK**2*az - ay*by*bz + bx*bz*cx -
     2    bz*cx**2 + ay*bz*cy + by*bz*cy - bz*cy**2 - bx**2*cz +
     3    ay*by*cz - by**2*cz + bx*cx*cz - ay*cy*cz + by*cy*cz
      BN=AI**2*ax + AL**2*ax - AJ*AK*ay - az*bx*bz - by**2*cx +
     2    az*bz*cx - bz**2*cx + bx*by*cy + by*cx*cy - bx*cy**2 +
     3    az*bx*cz + bx*bz*cz - az*cx*cz + bz*cx*cz - bx*cz**2
      BO=-(AJ*AK*ax) + AJ**2*ay + AL**2*ay - az*by*bz + bx*by*cx -
     2    by*cx**2 - bx**2*cy + az*bz*cy - bz**2*cy + bx*cx*cy +
     3    az*by*cz + by*bz*cz - az*cy*cz + bz*cy*cz - by*cz**2

      den_inv=1/(Sqrt(AQ)*((AA*AD + AB*AE + AC*AF)**2 +
     6 (AE*AV + AF*AW + AD*BA)**2/AQ))

      dax=(AA*AD**2*AS + AE**2*(-(AB*AJ*AK) + AL*AV) - AF**2*(AC*AJ*AL +
     2 AK*AW) + AD*AF*(-(AA*AJ*AL) + AC*AS - AK*BA) + 
     3 AE*(AF*(-(AJ*(AC*AK + AB*AL)) - AK*AV + AL*AW) + AD*(-(AA*AJ*AK)+
     4 AB*AS + AL*BA)))*den_inv
c      write(*,*)den_inv,dax

      day=(AC*AE*AF*bx**2 + AE*AF*AJ*ay*bx**2 - AD*AE*AL*ay*bx**2 + 
     -    AF**2*AJ*az*bx**2 - AD*AF*AL*az*bx**2 - AC*AD*AF*bx*by - 
     -    AE*AF*AJ*ax*bx*by + AD*AE*AL*ax*bx*by - AD*AF*AJ*ay*bx*by + 
     -    AD**2*AL*ay*bx*by + AD*AF*AJ*ax*by**2 - AD**2*AL*ax*by**2 + 
     -    AF**2*AJ*az*by**2 - AD*AF*AL*az*by**2 - AF**2*AJ*ax*bx*bz + 
     -    AD*AF*AL*ax*bx*bz - AD*AF*AJ*az*bx*bz + AD**2*AL*az*bx*bz - 
     -    AC*AF**2*by*bz - AF**2*AJ*ay*by*bz + AD*AF*AL*ay*by*bz - 
     -    AE*AF*AJ*az*by*bz + AD*AE*AL*az*by*bz + AC*AE*AF*bz**2 + 
     -    AD*AF*AJ*ax*bz**2 - AD**2*AL*ax*bz**2 + AE*AF*AJ*ay*bz**2 - 
     -    AD*AE*AL*ay*bz**2 - 2*AC*AE*AF*bx*cx - 2*AE*AF*AJ*ay*bx*cx + 
     -    2*AD*AE*AL*ay*bx*cx - 2*AF**2*AJ*az*bx*cx + 
     -    2*AD*AF*AL*az*bx*cx + AC*AD*AF*by*cx + AE*AF*AJ*ax*by*cx - 
     -    AD*AE*AL*ax*by*cx + AD*AF*AJ*ay*by*cx - AD**2*AL*ay*by*cx + 
     -    AE*AF*AJ*bx*by*cx - AD*AE*AL*bx*by*cx - AD*AF*AJ*by**2*cx + 
     -    AD**2*AL*by**2*cx + AF**2*AJ*ax*bz*cx - AD*AF*AL*ax*bz*cx + 
     -    AD*AF*AJ*az*bz*cx - AD**2*AL*az*bz*cx + AF**2*AJ*bx*bz*cx - 
     -    AD*AF*AL*bx*bz*cx - AD*AF*AJ*bz**2*cx + AD**2*AL*bz**2*cx + 
     -    AC*AE*AF*cx**2 + AE*AF*AJ*ay*cx**2 - AD*AE*AL*ay*cx**2 + 
     -    AF**2*AJ*az*cx**2 - AD*AF*AL*az*cx**2 - AE*AF*AJ*by*cx**2 + 
     -    AD*AE*AL*by*cx**2 - AF**2*AJ*bz*cx**2 + AD*AF*AL*bz*cx**2 + 
     -    AC*AD*AF*bx*cy + AE*AF*AJ*ax*bx*cy - AD*AE*AL*ax*bx*cy + 
     -    AD*AF*AJ*ay*bx*cy - AD**2*AL*ay*bx*cy - AE*AF*AJ*bx**2*cy + 
     -    AD*AE*AL*bx**2*cy - 2*AD*AF*AJ*ax*by*cy + 2*AD**2*AL*ax*by*cy-
     -    2*AF**2*AJ*az*by*cy + 2*AD*AF*AL*az*by*cy + AD*AF*AJ*bx*by*cy-
     -    AD**2*AL*bx*by*cy + AC*AF**2*bz*cy + AF**2*AJ*ay*bz*cy - 
     -    AD*AF*AL*ay*bz*cy + AE*AF*AJ*az*bz*cy - AD*AE*AL*az*bz*cy + 
     -    AF**2*AJ*by*bz*cy - AD*AF*AL*by*bz*cy - AE*AF*AJ*bz**2*cy + 
     -    AD*AE*AL*bz**2*cy - AC*AD*AF*cx*cy - AE*AF*AJ*ax*cx*cy + 
     -    AD*AE*AL*ax*cx*cy - AD*AF*AJ*ay*cx*cy + AD**2*AL*ay*cx*cy + 
     -    AE*AF*AJ*bx*cx*cy - AD*AE*AL*bx*cx*cy + AD*AF*AJ*by*cx*cy - 
     -    AD**2*AL*by*cx*cy + AD*AF*AJ*ax*cy**2 - AD**2*AL*ax*cy**2 + 
     -    AF**2*AJ*az*cy**2 - AD*AF*AL*az*cy**2 - AD*AF*AJ*bx*cy**2 + 
     -    AD**2*AL*bx*cy**2 - AF**2*AJ*bz*cy**2 + AD*AF*AL*bz*cy**2 + 
     -    AF**2*AJ*ax*bx*cz - AD*AF*AL*ax*bx*cz + AD*AF*AJ*az*bx*cz - 
     -    AD**2*AL*az*bx*cz - AF**2*AJ*bx**2*cz + AD*AF*AL*bx**2*cz + 
     -    AC*AF**2*by*cz + AF**2*AJ*ay*by*cz - AD*AF*AL*ay*by*cz + 
     -    AE*AF*AJ*az*by*cz - AD*AE*AL*az*by*cz - AF**2*AJ*by**2*cz + 
     -    AD*AF*AL*by**2*cz - 2*AC*AE*AF*bz*cz - 2*AD*AF*AJ*ax*bz*cz + 
     -    2*AD**2*AL*ax*bz*cz - 2*AE*AF*AJ*ay*bz*cz + 
     -    2*AD*AE*AL*ay*bz*cz + AD*AF*AJ*bx*bz*cz - AD**2*AL*bx*bz*cz + 
     -    AE*AF*AJ*by*bz*cz - AD*AE*AL*by*bz*cz - AF**2*AJ*ax*cx*cz + 
     -    AD*AF*AL*ax*cx*cz - AD*AF*AJ*az*cx*cz + AD**2*AL*az*cx*cz + 
     -    AF**2*AJ*bx*cx*cz - AD*AF*AL*bx*cx*cz + AD*AF*AJ*bz*cx*cz - 
     -    AD**2*AL*bz*cx*cz - AC*AF**2*cy*cz - AF**2*AJ*ay*cy*cz + 
     -    AD*AF*AL*ay*cy*cz - AE*AF*AJ*az*cy*cz + AD*AE*AL*az*cy*cz + 
     -    AF**2*AJ*by*cy*cz - AD*AF*AL*by*cy*cz + AE*AF*AJ*bz*cy*cz - 
     -    AD*AE*AL*bz*cy*cz + AC*AE*AF*cz**2 + AD*AF*AJ*ax*cz**2 - 
     -    AD**2*AL*ax*cz**2 + AE*AF*AJ*ay*cz**2 - AD*AE*AL*ay*cz**2 - 
     -    AD*AF*AJ*bx*cz**2 + AD**2*AL*bx*cz**2 - AE*AF*AJ*by*cz**2 + 
     -    AD*AE*AL*by*cz**2 + 
     -    AA*AD*(AE*AU + AK*(-(AD*bx) - AF*bz + AD*cx + AF*cz)) + 
     -    AB*AE*(AE*AU + AK*(-(AD*bx) - AF*bz + AD*cx + AF*cz)))/
     -  (Sqrt(AQ)*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))


      daz=(AC*AF**2*bx**2 - AE**2*AJ*ay*bx**2 + AD*AE*AK*ay*bx**2 - 
     -    AE*AF*AJ*az*bx**2 + AD*AF*AK*az*bx**2 + AE**2*AJ*ax*bx*by - 
     -    AD*AE*AK*ax*bx*by + AD*AE*AJ*ay*bx*by - AD**2*AK*ay*bx*by + 
     -    AC*AF**2*by**2 - AD*AE*AJ*ax*by**2 + AD**2*AK*ax*by**2 - 
     -    AE*AF*AJ*az*by**2 + AD*AF*AK*az*by**2 - AC*AD*AF*bx*bz + 
     -    AE*AF*AJ*ax*bx*bz - AD*AF*AK*ax*bx*bz + AD*AE*AJ*az*bx*bz - 
     -    AD**2*AK*az*bx*bz - AC*AE*AF*by*bz + AE*AF*AJ*ay*by*bz - 
     -    AD*AF*AK*ay*by*bz + AE**2*AJ*az*by*bz - AD*AE*AK*az*by*bz - 
     -    AD*AE*AJ*ax*bz**2 + AD**2*AK*ax*bz**2 - AE**2*AJ*ay*bz**2 + 
     -    AD*AE*AK*ay*bz**2 - 2*AC*AF**2*bx*cx + 2*AE**2*AJ*ay*bx*cx - 
     -    2*AD*AE*AK*ay*bx*cx +2*AE*AF*AJ*az*bx*cx -2*AD*AF*AK*az*bx*cx-
     -    AE**2*AJ*ax*by*cx + AD*AE*AK*ax*by*cx - AD*AE*AJ*ay*by*cx + 
     -    AD**2*AK*ay*by*cx - AE**2*AJ*bx*by*cx + AD*AE*AK*bx*by*cx + 
     -    AD*AE*AJ*by**2*cx - AD**2*AK*by**2*cx + AC*AD*AF*bz*cx - 
     -    AE*AF*AJ*ax*bz*cx + AD*AF*AK*ax*bz*cx - AD*AE*AJ*az*bz*cx + 
     -    AD**2*AK*az*bz*cx - AE*AF*AJ*bx*bz*cx + AD*AF*AK*bx*bz*cx + 
     -    AD*AE*AJ*bz**2*cx - AD**2*AK*bz**2*cx + AC*AF**2*cx**2 - 
     -    AE**2*AJ*ay*cx**2 + AD*AE*AK*ay*cx**2 - AE*AF*AJ*az*cx**2 + 
     -    AD*AF*AK*az*cx**2 + AE**2*AJ*by*cx**2 - AD*AE*AK*by*cx**2 + 
     -    AE*AF*AJ*bz*cx**2 - AD*AF*AK*bz*cx**2 - AE**2*AJ*ax*bx*cy + 
     -    AD*AE*AK*ax*bx*cy - AD*AE*AJ*ay*bx*cy + AD**2*AK*ay*bx*cy + 
     -    AE**2*AJ*bx**2*cy - AD*AE*AK*bx**2*cy - 2*AC*AF**2*by*cy + 
     -    2*AD*AE*AJ*ax*by*cy -2*AD**2*AK*ax*by*cy +2*AE*AF*AJ*az*by*cy-
     -    2*AD*AF*AK*az*by*cy - AD*AE*AJ*bx*by*cy + AD**2*AK*bx*by*cy + 
     -    AC*AE*AF*bz*cy - AE*AF*AJ*ay*bz*cy + AD*AF*AK*ay*bz*cy - 
     -    AE**2*AJ*az*bz*cy + AD*AE*AK*az*bz*cy - AE*AF*AJ*by*bz*cy + 
     -    AD*AF*AK*by*bz*cy + AE**2*AJ*bz**2*cy - AD*AE*AK*bz**2*cy + 
     -    AE**2*AJ*ax*cx*cy - AD*AE*AK*ax*cx*cy + AD*AE*AJ*ay*cx*cy - 
     -    AD**2*AK*ay*cx*cy - AE**2*AJ*bx*cx*cy + AD*AE*AK*bx*cx*cy - 
     -    AD*AE*AJ*by*cx*cy + AD**2*AK*by*cx*cy + AC*AF**2*cy**2 - 
     -    AD*AE*AJ*ax*cy**2 + AD**2*AK*ax*cy**2 - AE*AF*AJ*az*cy**2 + 
     -    AD*AF*AK*az*cy**2 + AD*AE*AJ*bx*cy**2 - AD**2*AK*bx*cy**2 + 
     -    AE*AF*AJ*bz*cy**2 - AD*AF*AK*bz*cy**2 + 
     -    AA*AD*(AF*AT + AL*(-(AD*bx) - AE*by + AD*cx + AE*cy)) + 
     -    AB*AE*(AF*AT + AL*(-(AD*bx) - AE*by + AD*cx + AE*cy)) + 
     -    AC*AD*AF*bx*cz - AE*AF*AJ*ax*bx*cz + AD*AF*AK*ax*bx*cz - 
     -    AD*AE*AJ*az*bx*cz + AD**2*AK*az*bx*cz + AE*AF*AJ*bx**2*cz - 
     -    AD*AF*AK*bx**2*cz + AC*AE*AF*by*cz - AE*AF*AJ*ay*by*cz + 
     -    AD*AF*AK*ay*by*cz - AE**2*AJ*az*by*cz + AD*AE*AK*az*by*cz + 
     -    AE*AF*AJ*by**2*cz - AD*AF*AK*by**2*cz + 2*AD*AE*AJ*ax*bz*cz - 
     -    2*AD**2*AK*ax*bz*cz +2*AE**2*AJ*ay*bz*cz -2*AD*AE*AK*ay*bz*cz-
     -    AD*AE*AJ*bx*bz*cz + AD**2*AK*bx*bz*cz - AE**2*AJ*by*bz*cz + 
     -    AD*AE*AK*by*bz*cz - AC*AD*AF*cx*cz + AE*AF*AJ*ax*cx*cz - 
     -    AD*AF*AK*ax*cx*cz + AD*AE*AJ*az*cx*cz - AD**2*AK*az*cx*cz - 
     -    AE*AF*AJ*bx*cx*cz + AD*AF*AK*bx*cx*cz - AD*AE*AJ*bz*cx*cz + 
     -    AD**2*AK*bz*cx*cz - AC*AE*AF*cy*cz + AE*AF*AJ*ay*cy*cz - 
     -    AD*AF*AK*ay*cy*cz + AE**2*AJ*az*cy*cz - AD*AE*AK*az*cy*cz - 
     -    AE*AF*AJ*by*cy*cz + AD*AF*AK*by*cy*cz - AE**2*AJ*bz*cy*cz + 
     -    AD*AE*AK*bz*cy*cz - AD*AE*AJ*ax*cz**2 + AD**2*AK*ax*cz**2 - 
     -    AE**2*AJ*ay*cz**2 + AD*AE*AK*ay*cz**2 + AD*AE*AJ*bx*cz**2 - 
     -    AD**2*AK*bx*cz**2 + AE**2*AJ*by*cz**2 - AD*AE*AK*by*cz**2)/
     -  (Sqrt(AQ)*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))

      dbx=(2*AC*AE*AF*AQ*ay*bx + 2*AC*AF**2*AQ*az*bx - 
     -    AC*AE*AF*AJ*ay*bx**2 - AC*AE*AN*AQ*ay*bx**2 - 
     -    AC*AF*AO*AQ*ay*bx**2 - AC*AF**2*AJ*az*bx**2 + 
     -    AE*AF*AQ*ay*BC*bx**2 + AF**2*AQ*az*BC*bx**2 - 
     -    AE**2*AQ*ay*BD*bx**2 - AE*AF*AQ*az*BD*bx**2 - 
     -    AC*AE*AF*AQ*ax*by - AC*AD*AF*AQ*ay*by + 
     -    AC*AE*AF*AJ*ax*bx*by + AC*AE*AN*AQ*ax*bx*by + 
     -    AC*AF*AO*AQ*ax*bx*by + AC*AD*AF*AJ*ay*bx*by + 
     -    AC*AD*AN*AQ*ay*bx*by - AE*AF*AQ*ax*BC*bx*by - 
     -    AD*AF*AQ*ay*BC*bx*by + AE**2*AQ*ax*BD*bx*by + 
     -    AD*AE*AQ*ay*BD*bx*by - AC*AD*AF*AJ*ax*by**2 - 
     -    AC*AD*AN*AQ*ax*by**2 - AC*AF**2*AJ*az*by**2 + 
     -    AD*AF*AQ*ax*BC*by**2 + AF**2*AQ*az*BC*by**2 - 
     -    AD*AE*AQ*ax*BD*by**2 - AE*AF*AQ*az*BD*by**2 - 
     -    AC*AF**2*AQ*ax*bz - AC*AD*AF*AQ*az*bz + 
     -    AC*AF**2*AJ*ax*bx*bz + AC*AD*AF*AJ*az*bx*bz + 
     -    AC*AD*AN*AQ*az*bx*bz - AF**2*AQ*ax*BC*bx*bz - 
     -    AD*AF*AQ*az*BC*bx*bz + AE*AF*AQ*ax*BD*bx*bz + 
     -    AD*AE*AQ*az*BD*bx*bz + AC*AF**2*AJ*ay*by*bz + 
     -    AC*AE*AF*AJ*az*by*bz + AC*AE*AN*AQ*az*by*bz + 
     -    AC*AF*AO*AQ*az*by*bz - AF**2*AQ*ay*BC*by*bz - 
     -    AE*AF*AQ*az*BC*by*bz + AE*AF*AQ*ay*BD*by*bz + 
     -    AE**2*AQ*az*BD*by*bz - AC*AD*AF*AJ*ax*bz**2 - 
     -    AC*AD*AN*AQ*ax*bz**2 - AC*AE*AF*AJ*ay*bz**2 - 
     -    AC*AE*AN*AQ*ay*bz**2 - AC*AF*AO*AQ*ay*bz**2 + 
     -    AD*AF*AQ*ax*BC*bz**2 + AE*AF*AQ*ay*BC*bz**2 - 
     -    AD*AE*AQ*ax*BD*bz**2 - AE**2*AQ*ay*BD*bz**2 - 
     -    2*AC*AE*AF*AQ*ay*cx - 2*AC*AF**2*AQ*az*cx + 
     -    2*AC*AE*AF*AJ*ay*bx*cx + 2*AC*AE*AN*AQ*ay*bx*cx + 
     -    2*AC*AF*AO*AQ*ay*bx*cx + 2*AC*AF**2*AJ*az*bx*cx - 
     -    2*AE*AF*AQ*ay*BC*bx*cx - 2*AF**2*AQ*az*BC*bx*cx + 
     -    2*AE**2*AQ*ay*BD*bx*cx + 2*AE*AF*AQ*az*BD*bx*cx + 
     -    AC*AE*AF*AQ*by*cx - AC*AE*AF*AJ*ax*by*cx - 
     -    AC*AE*AN*AQ*ax*by*cx - AC*AF*AO*AQ*ax*by*cx - 
     -    AC*AD*AF*AJ*ay*by*cx - AC*AD*AN*AQ*ay*by*cx + 
     -    AE*AF*AQ*ax*BC*by*cx + AD*AF*AQ*ay*BC*by*cx - 
     -    AE**2*AQ*ax*BD*by*cx - AD*AE*AQ*ay*BD*by*cx - 
     -    AC*AE*AF*AJ*bx*by*cx - AC*AE*AN*AQ*bx*by*cx - 
     -    AC*AF*AO*AQ*bx*by*cx + AE*AF*AQ*BC*bx*by*cx - 
     -    AE**2*AQ*BD*bx*by*cx + AC*AD*AF*AJ*by**2*cx + 
     -    AC*AD*AN*AQ*by**2*cx - AD*AF*AQ*BC*by**2*cx + 
     -    AD*AE*AQ*BD*by**2*cx + AC*AF**2*AQ*bz*cx - 
     -    AC*AF**2*AJ*ax*bz*cx - AC*AD*AF*AJ*az*bz*cx - 
     -    AC*AD*AN*AQ*az*bz*cx + AF**2*AQ*ax*BC*bz*cx + 
     -    AD*AF*AQ*az*BC*bz*cx - AE*AF*AQ*ax*BD*bz*cx - 
     -    AD*AE*AQ*az*BD*bz*cx - AC*AF**2*AJ*bx*bz*cx + 
     -    AF**2*AQ*BC*bx*bz*cx - AE*AF*AQ*BD*bx*bz*cx + 
     -    AC*AD*AF*AJ*bz**2*cx + AC*AD*AN*AQ*bz**2*cx - 
     -    AD*AF*AQ*BC*bz**2*cx + AD*AE*AQ*BD*bz**2*cx - 
     -    AC*AE*AF*AJ*ay*cx**2 - AC*AE*AN*AQ*ay*cx**2 - 
     -    AC*AF*AO*AQ*ay*cx**2 - AC*AF**2*AJ*az*cx**2 + 
     -    AE*AF*AQ*ay*BC*cx**2 + AF**2*AQ*az*BC*cx**2 - 
     -    AE**2*AQ*ay*BD*cx**2 - AE*AF*AQ*az*BD*cx**2 + 
     -    AC*AE*AF*AJ*by*cx**2 + AC*AE*AN*AQ*by*cx**2 + 
     -    AC*AF*AO*AQ*by*cx**2 - AE*AF*AQ*BC*by*cx**2 + 
     -    AE**2*AQ*BD*by*cx**2 + AC*AF**2*AJ*bz*cx**2 - 
     -    AF**2*AQ*BC*bz*cx**2 + AE*AF*AQ*BD*bz*cx**2 + 
     -    AC*AE*AF*AQ*ax*cy + AC*AD*AF*AQ*ay*cy - 2*AC*AE*AF*AQ*bx*cy - 
     -    AC*AE*AF*AJ*ax*bx*cy - AC*AE*AN*AQ*ax*bx*cy - 
     -    AC*AF*AO*AQ*ax*bx*cy - AC*AD*AF*AJ*ay*bx*cy - 
     -    AC*AD*AN*AQ*ay*bx*cy + AE*AF*AQ*ax*BC*bx*cy + 
     -    AD*AF*AQ*ay*BC*bx*cy - AE**2*AQ*ax*BD*bx*cy - 
     -    AD*AE*AQ*ay*BD*bx*cy + AC*AE*AF*AJ*bx**2*cy + 
     -    AC*AE*AN*AQ*bx**2*cy + AC*AF*AO*AQ*bx**2*cy - 
     -    AE*AF*AQ*BC*bx**2*cy + AE**2*AQ*BD*bx**2*cy + 
     -    AC*AD*AF*AQ*by*cy + 2*AC*AD*AF*AJ*ax*by*cy + 
     -    2*AC*AD*AN*AQ*ax*by*cy + 2*AC*AF**2*AJ*az*by*cy - 
     -    2*AD*AF*AQ*ax*BC*by*cy - 2*AF**2*AQ*az*BC*by*cy + 
     -    2*AD*AE*AQ*ax*BD*by*cy + 2*AE*AF*AQ*az*BD*by*cy - 
     -    AC*AD*AF*AJ*bx*by*cy - AC*AD*AN*AQ*bx*by*cy + 
     -    AD*AF*AQ*BC*bx*by*cy - AD*AE*AQ*BD*bx*by*cy - 
     -    AC*AF**2*AJ*ay*bz*cy - AC*AE*AF*AJ*az*bz*cy - 
     -    AC*AE*AN*AQ*az*bz*cy - AC*AF*AO*AQ*az*bz*cy + 
     -    AF**2*AQ*ay*BC*bz*cy + AE*AF*AQ*az*BC*bz*cy - 
     -    AE*AF*AQ*ay*BD*bz*cy - AE**2*AQ*az*BD*bz*cy - 
     -    AC*AF**2*AJ*by*bz*cy + AF**2*AQ*BC*by*bz*cy - 
     -    AE*AF*AQ*BD*by*bz*cy + AC*AE*AF*AJ*bz**2*cy + 
     -    AC*AE*AN*AQ*bz**2*cy + AC*AF*AO*AQ*bz**2*cy - 
     -    AE*AF*AQ*BC*bz**2*cy + AE**2*AQ*BD*bz**2*cy + 
     -    AC*AE*AF*AQ*cx*cy + AC*AE*AF*AJ*ax*cx*cy + 
     -    AC*AE*AN*AQ*ax*cx*cy + AC*AF*AO*AQ*ax*cx*cy + 
     -    AC*AD*AF*AJ*ay*cx*cy + AC*AD*AN*AQ*ay*cx*cy - 
     -    AE*AF*AQ*ax*BC*cx*cy - AD*AF*AQ*ay*BC*cx*cy + 
     -    AE**2*AQ*ax*BD*cx*cy + AD*AE*AQ*ay*BD*cx*cy - 
     -    AC*AE*AF*AJ*bx*cx*cy - AC*AE*AN*AQ*bx*cx*cy - 
     -    AC*AF*AO*AQ*bx*cx*cy + AE*AF*AQ*BC*bx*cx*cy - 
     -    AE**2*AQ*BD*bx*cx*cy - AC*AD*AF*AJ*by*cx*cy - 
     -    AC*AD*AN*AQ*by*cx*cy + AD*AF*AQ*BC*by*cx*cy - 
     -    AD*AE*AQ*BD*by*cx*cy - AC*AD*AF*AQ*cy**2 - 
     -    AC*AD*AF*AJ*ax*cy**2 - AC*AD*AN*AQ*ax*cy**2 - 
     -    AC*AF**2*AJ*az*cy**2 + AD*AF*AQ*ax*BC*cy**2 + 
     -    AF**2*AQ*az*BC*cy**2 - AD*AE*AQ*ax*BD*cy**2 - 
     -    AE*AF*AQ*az*BD*cy**2 + AC*AD*AF*AJ*bx*cy**2 + 
     -    AC*AD*AN*AQ*bx*cy**2 - AD*AF*AQ*BC*bx*cy**2 + 
     -    AD*AE*AQ*BD*bx*cy**2 + AC*AF**2*AJ*bz*cy**2 - 
     -    AF**2*AQ*BC*bz*cy**2 + AE*AF*AQ*BD*bz*cy**2 + 
     -    AC*AF**2*AQ*ax*cz + AC*AD*AF*AQ*az*cz - 2*AC*AF**2*AQ*bx*cz - 
     -    AC*AF**2*AJ*ax*bx*cz - AC*AD*AF*AJ*az*bx*cz - 
     -    AC*AD*AN*AQ*az*bx*cz + AF**2*AQ*ax*BC*bx*cz + 
     -    AD*AF*AQ*az*BC*bx*cz - AE*AF*AQ*ax*BD*bx*cz - 
     -    AD*AE*AQ*az*BD*bx*cz + AC*AF**2*AJ*bx**2*cz - 
     -    AF**2*AQ*BC*bx**2*cz + AE*AF*AQ*BD*bx**2*cz - 
     -    AC*AF**2*AJ*ay*by*cz - AC*AE*AF*AJ*az*by*cz - 
     -    AC*AE*AN*AQ*az*by*cz - AC*AF*AO*AQ*az*by*cz + 
     -    AF**2*AQ*ay*BC*by*cz + AE*AF*AQ*az*BC*by*cz - 
     -    AE*AF*AQ*ay*BD*by*cz - AE**2*AQ*az*BD*by*cz + 
     -    AC*AF**2*AJ*by**2*cz - AF**2*AQ*BC*by**2*cz + 
     -    AE*AF*AQ*BD*by**2*cz + AC*AD*AF*AQ*bz*cz + 
     -    2*AC*AD*AF*AJ*ax*bz*cz + 2*AC*AD*AN*AQ*ax*bz*cz + 
     -    2*AC*AE*AF*AJ*ay*bz*cz + 2*AC*AE*AN*AQ*ay*bz*cz + 
     -    2*AC*AF*AO*AQ*ay*bz*cz - 2*AD*AF*AQ*ax*BC*bz*cz - 
     -    2*AE*AF*AQ*ay*BC*bz*cz + 2*AD*AE*AQ*ax*BD*bz*cz + 
     -    2*AE**2*AQ*ay*BD*bz*cz - AC*AD*AF*AJ*bx*bz*cz - 
     -    AC*AD*AN*AQ*bx*bz*cz + AD*AF*AQ*BC*bx*bz*cz - 
     -    AD*AE*AQ*BD*bx*bz*cz - AC*AE*AF*AJ*by*bz*cz - 
     -    AC*AE*AN*AQ*by*bz*cz - AC*AF*AO*AQ*by*bz*cz + 
     -    AE*AF*AQ*BC*by*bz*cz - AE**2*AQ*BD*by*bz*cz + 
     -    AC*AF**2*AQ*cx*cz + AC*AF**2*AJ*ax*cx*cz + 
     -    AC*AD*AF*AJ*az*cx*cz + AC*AD*AN*AQ*az*cx*cz - 
     -    AF**2*AQ*ax*BC*cx*cz - AD*AF*AQ*az*BC*cx*cz + 
     -    AE*AF*AQ*ax*BD*cx*cz + AD*AE*AQ*az*BD*cx*cz - 
     -    AC*AF**2*AJ*bx*cx*cz + AF**2*AQ*BC*bx*cx*cz - 
     -    AE*AF*AQ*BD*bx*cx*cz - AC*AD*AF*AJ*bz*cx*cz - 
     -    AC*AD*AN*AQ*bz*cx*cz + AD*AF*AQ*BC*bz*cx*cz - 
     -    AD*AE*AQ*BD*bz*cx*cz + AC*AF**2*AJ*ay*cy*cz + 
     -    AC*AE*AF*AJ*az*cy*cz + AC*AE*AN*AQ*az*cy*cz + 
     -    AC*AF*AO*AQ*az*cy*cz - AF**2*AQ*ay*BC*cy*cz - 
     -    AE*AF*AQ*az*BC*cy*cz + AE*AF*AQ*ay*BD*cy*cz + 
     -    AE**2*AQ*az*BD*cy*cz - AC*AF**2*AJ*by*cy*cz + 
     -    AF**2*AQ*BC*by*cy*cz - AE*AF*AQ*BD*by*cy*cz - 
     -    AC*AE*AF*AJ*bz*cy*cz - AC*AE*AN*AQ*bz*cy*cz - 
     -    AC*AF*AO*AQ*bz*cy*cz + AE*AF*AQ*BC*bz*cy*cz - 
     -    AE**2*AQ*BD*bz*cy*cz - AC*AD*AF*AQ*cz**2 - 
     -    AC*AD*AF*AJ*ax*cz**2 - AC*AD*AN*AQ*ax*cz**2 - 
     -    AC*AE*AF*AJ*ay*cz**2 - AC*AE*AN*AQ*ay*cz**2 - 
     -    AC*AF*AO*AQ*ay*cz**2 + AD*AF*AQ*ax*BC*cz**2 + 
     -    AE*AF*AQ*ay*BC*cz**2 - AD*AE*AQ*ax*BD*cz**2 - 
     -    AE**2*AQ*ay*BD*cz**2 + AC*AD*AF*AJ*bx*cz**2 + 
     -    AC*AD*AN*AQ*bx*cz**2 - AD*AF*AQ*BC*bx*cz**2 + 
     -    AD*AE*AQ*BD*bx*cz**2 + AC*AE*AF*AJ*by*cz**2 + 
     -    AC*AE*AN*AQ*by*cz**2 + AC*AF*AO*AQ*by*cz**2 - 
     -    AE*AF*AQ*BC*by*cz**2 + AE**2*AQ*BD*by*cz**2 - 
     -    AA*AD*(AO*AQ*ay*bx**2 - AN*AQ*az*bx**2 + AD*AQ*ay*by - 
     -       AO*AQ*ax*bx*by - AD*AJ*ay*bx*by + AD*AJ*ax*by**2 - 
     -       AN*AQ*az*by**2 + AD*AQ*az*bz + AN*AQ*ax*bx*bz - 
     -       AD*AJ*az*bx*bz + AN*AQ*ay*by*bz - AO*AQ*az*by*bz + 
     -       AD*AJ*ax*bz**2 + AO*AQ*ay*bz**2 - 2*AO*AQ*ay*bx*cx + 
     -       2*AN*AQ*az*bx*cx + AO*AQ*ax*by*cx + AD*AJ*ay*by*cx + 
     -       AO*AQ*bx*by*cx - AD*AJ*by**2*cx - AN*AQ*ax*bz*cx + 
     -       AD*AJ*az*bz*cx - AN*AQ*bx*bz*cx - AD*AJ*bz**2*cx + 
     -       AO*AQ*ay*cx**2 - AN*AQ*az*cx**2 - AO*AQ*by*cx**2 + 
     -       AN*AQ*bz*cx**2 - AD*AQ*ay*cy + AO*AQ*ax*bx*cy + 
     -       AD*AJ*ay*bx*cy - AO*AQ*bx**2*cy - AD*AQ*by*cy - 
     -       2*AD*AJ*ax*by*cy + 2*AN*AQ*az*by*cy + AD*AJ*bx*by*cy - 
     -       AN*AQ*ay*bz*cy + AO*AQ*az*bz*cy - AN*AQ*by*bz*cy - 
     -       AO*AQ*bz**2*cy - AO*AQ*ax*cx*cy - AD*AJ*ay*cx*cy + 
     -       AO*AQ*bx*cx*cy + AD*AJ*by*cx*cy + AD*AQ*cy**2 + 
     -       AD*AJ*ax*cy**2 - AN*AQ*az*cy**2 - AD*AJ*bx*cy**2 + 
     -       AN*AQ*bz*cy**2 + 
     -       AE*(AJ*AV + AQ*(AK*ax - 2*AJ*ay - by*cx + 2*bx*cy - 
     -             cx*cy)) - AD*AQ*az*cz - AN*AQ*ax*bx*cz + 
     -       AD*AJ*az*bx*cz + AN*AQ*bx**2*cz - AN*AQ*ay*by*cz + 
     -       AO*AQ*az*by*cz + AN*AQ*by**2*cz - AD*AQ*bz*cz - 
     -       2*AD*AJ*ax*bz*cz - 2*AO*AQ*ay*bz*cz + AD*AJ*bx*bz*cz + 
     -       AO*AQ*by*bz*cz + AN*AQ*ax*cx*cz - AD*AJ*az*cx*cz - 
     -       AN*AQ*bx*cx*cz + AD*AJ*bz*cx*cz + AN*AQ*ay*cy*cz - 
     -       AO*AQ*az*cy*cz - AN*AQ*by*cy*cz + AO*AQ*bz*cy*cz + 
     -       AD*AQ*cz**2 + AD*AJ*ax*cz**2 + AO*AQ*ay*cz**2 - 
     -       AD*AJ*bx*cz**2 - AO*AQ*by*cz**2 + 
     -       AF*(AJ*AW + AQ*(AL*ax - 2*AJ*az - bz*cx + 2*bx*cz - cx*cz))
     -       ) + AB*(AO*AQ*(AF*AW + AD*BA) + 
     -       AE**2*(-(AJ*AV) + 
     -          AQ*(-(AK*ax) + 2*AJ*ay + by*cx - 2*bx*cy + cx*cy)) + 
     -       AE*(AN*AQ*AW + AF*
     -           (-(AJ*AW) + 
     -             AQ*(-(AL*ax) + 2*AJ*az + bz*cx - 2*bx*cz + cx*cz)) - 
     -          AD*(AJ*BA + AQ*
     -           (AK*ay + AL*az - by*cy + cy**2 - bz*cz + cz**2)))))/
     -  (AQ**1.5*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))


      dby=(-(AC*AE*AF*AQ*ax*bx) - AC*AD*AF*AQ*ay*bx - 
     -    AC*AE*AF*AK*ay*bx**2 + AC*AE*AM*AQ*ay*bx**2 - 
     -    AC*AF**2*AK*az*bx**2 - AE*AF*AQ*ay*BB*bx**2 - 
     -    AF**2*AQ*az*BB*bx**2 + AD*AE*AQ*ay*BD*bx**2 + 
     -    AD*AF*AQ*az*BD*bx**2 + 2*AC*AD*AF*AQ*ax*by + 
     -    2*AC*AF**2*AQ*az*by + AC*AE*AF*AK*ax*bx*by - 
     -    AC*AE*AM*AQ*ax*bx*by + AC*AD*AF*AK*ay*bx*by - 
     -    AC*AD*AM*AQ*ay*bx*by - AC*AF*AO*AQ*ay*bx*by + 
     -    AE*AF*AQ*ax*BB*bx*by + AD*AF*AQ*ay*BB*bx*by - 
     -    AD*AE*AQ*ax*BD*bx*by - AD**2*AQ*ay*BD*bx*by - 
     -    AC*AD*AF*AK*ax*by**2 + AC*AD*AM*AQ*ax*by**2 + 
     -    AC*AF*AO*AQ*ax*by**2 - AC*AF**2*AK*az*by**2 - 
     -    AD*AF*AQ*ax*BB*by**2 - AF**2*AQ*az*BB*by**2 + 
     -    AD**2*AQ*ax*BD*by**2 + AD*AF*AQ*az*BD*by**2 - 
     -    AC*AF**2*AQ*ay*bz - AC*AE*AF*AQ*az*bz + 
     -    AC*AF**2*AK*ax*bx*bz + AC*AD*AF*AK*az*bx*bz - 
     -    AC*AD*AM*AQ*az*bx*bz - AC*AF*AO*AQ*az*bx*bz + 
     -    AF**2*AQ*ax*BB*bx*bz + AD*AF*AQ*az*BB*bx*bz - 
     -    AD*AF*AQ*ax*BD*bx*bz - AD**2*AQ*az*BD*bx*bz + 
     -    AC*AF**2*AK*ay*by*bz + AC*AE*AF*AK*az*by*bz - 
     -    AC*AE*AM*AQ*az*by*bz + AF**2*AQ*ay*BB*by*bz + 
     -    AE*AF*AQ*az*BB*by*bz - AD*AF*AQ*ay*BD*by*bz - 
     -    AD*AE*AQ*az*BD*by*bz - AC*AD*AF*AK*ax*bz**2 + 
     -    AC*AD*AM*AQ*ax*bz**2 + AC*AF*AO*AQ*ax*bz**2 - 
     -    AC*AE*AF*AK*ay*bz**2 + AC*AE*AM*AQ*ay*bz**2 - 
     -    AD*AF*AQ*ax*BB*bz**2 - AE*AF*AQ*ay*BB*bz**2 + 
     -    AD**2*AQ*ax*BD*bz**2 + AD*AE*AQ*ay*BD*bz**2 + 
     -    AC*AE*AF*AQ*ax*cx + AC*AD*AF*AQ*ay*cx + 
     -    AC*AE*AF*AQ*bx*cx + 2*AC*AE*AF*AK*ay*bx*cx - 
     -    2*AC*AE*AM*AQ*ay*bx*cx + 2*AC*AF**2*AK*az*bx*cx + 
     -    2*AE*AF*AQ*ay*BB*bx*cx + 2*AF**2*AQ*az*BB*bx*cx - 
     -    2*AD*AE*AQ*ay*BD*bx*cx - 2*AD*AF*AQ*az*BD*bx*cx - 
     -    2*AC*AD*AF*AQ*by*cx - AC*AE*AF*AK*ax*by*cx + 
     -    AC*AE*AM*AQ*ax*by*cx - AC*AD*AF*AK*ay*by*cx + 
     -    AC*AD*AM*AQ*ay*by*cx + AC*AF*AO*AQ*ay*by*cx - 
     -    AE*AF*AQ*ax*BB*by*cx - AD*AF*AQ*ay*BB*by*cx + 
     -    AD*AE*AQ*ax*BD*by*cx + AD**2*AQ*ay*BD*by*cx - 
     -    AC*AE*AF*AK*bx*by*cx + AC*AE*AM*AQ*bx*by*cx - 
     -    AE*AF*AQ*BB*bx*by*cx + AD*AE*AQ*BD*bx*by*cx + 
     -    AC*AD*AF*AK*by**2*cx - AC*AD*AM*AQ*by**2*cx - 
     -    AC*AF*AO*AQ*by**2*cx + AD*AF*AQ*BB*by**2*cx - 
     -    AD**2*AQ*BD*by**2*cx - AC*AF**2*AK*ax*bz*cx - 
     -    AC*AD*AF*AK*az*bz*cx + AC*AD*AM*AQ*az*bz*cx + 
     -    AC*AF*AO*AQ*az*bz*cx - AF**2*AQ*ax*BB*bz*cx - 
     -    AD*AF*AQ*az*BB*bz*cx + AD*AF*AQ*ax*BD*bz*cx + 
     -    AD**2*AQ*az*BD*bz*cx - AC*AF**2*AK*bx*bz*cx - 
     -    AF**2*AQ*BB*bx*bz*cx + AD*AF*AQ*BD*bx*bz*cx + 
     -    AC*AD*AF*AK*bz**2*cx - AC*AD*AM*AQ*bz**2*cx - 
     -    AC*AF*AO*AQ*bz**2*cx + AD*AF*AQ*BB*bz**2*cx - 
     -    AD**2*AQ*BD*bz**2*cx - AC*AE*AF*AQ*cx**2 - 
     -    AC*AE*AF*AK*ay*cx**2 + AC*AE*AM*AQ*ay*cx**2 - 
     -    AC*AF**2*AK*az*cx**2 - AE*AF*AQ*ay*BB*cx**2 - 
     -    AF**2*AQ*az*BB*cx**2 + AD*AE*AQ*ay*BD*cx**2 + 
     -    AD*AF*AQ*az*BD*cx**2 + AC*AE*AF*AK*by*cx**2 - 
     -    AC*AE*AM*AQ*by*cx**2 + AE*AF*AQ*BB*by*cx**2 - 
     -    AD*AE*AQ*BD*by*cx**2 + AC*AF**2*AK*bz*cx**2 + 
     -    AF**2*AQ*BB*bz*cx**2 - AD*AF*AQ*BD*bz*cx**2 - 
     -    2*AC*AD*AF*AQ*ax*cy - 2*AC*AF**2*AQ*az*cy + 
     -    AC*AD*AF*AQ*bx*cy - AC*AE*AF*AK*ax*bx*cy + 
     -    AC*AE*AM*AQ*ax*bx*cy - AC*AD*AF*AK*ay*bx*cy + 
     -    AC*AD*AM*AQ*ay*bx*cy + AC*AF*AO*AQ*ay*bx*cy - 
     -    AE*AF*AQ*ax*BB*bx*cy - AD*AF*AQ*ay*BB*bx*cy + 
     -    AD*AE*AQ*ax*BD*bx*cy + AD**2*AQ*ay*BD*bx*cy + 
     -    AC*AE*AF*AK*bx**2*cy - AC*AE*AM*AQ*bx**2*cy + 
     -    AE*AF*AQ*BB*bx**2*cy - AD*AE*AQ*BD*bx**2*cy + 
     -    2*AC*AD*AF*AK*ax*by*cy - 2*AC*AD*AM*AQ*ax*by*cy - 
     -    2*AC*AF*AO*AQ*ax*by*cy + 2*AC*AF**2*AK*az*by*cy + 
     -    2*AD*AF*AQ*ax*BB*by*cy + 2*AF**2*AQ*az*BB*by*cy - 
     -    2*AD**2*AQ*ax*BD*by*cy - 2*AD*AF*AQ*az*BD*by*cy - 
     -    AC*AD*AF*AK*bx*by*cy + AC*AD*AM*AQ*bx*by*cy + 
     -    AC*AF*AO*AQ*bx*by*cy - AD*AF*AQ*BB*bx*by*cy + 
     -    AD**2*AQ*BD*bx*by*cy + AC*AF**2*AQ*bz*cy - 
     -    AC*AF**2*AK*ay*bz*cy - AC*AE*AF*AK*az*bz*cy + 
     -    AC*AE*AM*AQ*az*bz*cy - AF**2*AQ*ay*BB*bz*cy - 
     -    AE*AF*AQ*az*BB*bz*cy + AD*AF*AQ*ay*BD*bz*cy + 
     -    AD*AE*AQ*az*BD*bz*cy - AC*AF**2*AK*by*bz*cy - 
     -    AF**2*AQ*BB*by*bz*cy + AD*AF*AQ*BD*by*bz*cy + 
     -    AC*AE*AF*AK*bz**2*cy - AC*AE*AM*AQ*bz**2*cy + 
     -    AE*AF*AQ*BB*bz**2*cy - AD*AE*AQ*BD*bz**2*cy + 
     -    AC*AD*AF*AQ*cx*cy + AC*AE*AF*AK*ax*cx*cy - 
     -    AC*AE*AM*AQ*ax*cx*cy + AC*AD*AF*AK*ay*cx*cy - 
     -    AC*AD*AM*AQ*ay*cx*cy - AC*AF*AO*AQ*ay*cx*cy + 
     -    AE*AF*AQ*ax*BB*cx*cy + AD*AF*AQ*ay*BB*cx*cy - 
     -    AD*AE*AQ*ax*BD*cx*cy - AD**2*AQ*ay*BD*cx*cy - 
     -    AC*AE*AF*AK*bx*cx*cy + AC*AE*AM*AQ*bx*cx*cy - 
     -    AE*AF*AQ*BB*bx*cx*cy + AD*AE*AQ*BD*bx*cx*cy - 
     -    AC*AD*AF*AK*by*cx*cy + AC*AD*AM*AQ*by*cx*cy + 
     -    AC*AF*AO*AQ*by*cx*cy - AD*AF*AQ*BB*by*cx*cy + 
     -    AD**2*AQ*BD*by*cx*cy - AC*AD*AF*AK*ax*cy**2 + 
     -    AC*AD*AM*AQ*ax*cy**2 + AC*AF*AO*AQ*ax*cy**2 - 
     -    AC*AF**2*AK*az*cy**2 - AD*AF*AQ*ax*BB*cy**2 - 
     -    AF**2*AQ*az*BB*cy**2 + AD**2*AQ*ax*BD*cy**2 + 
     -    AD*AF*AQ*az*BD*cy**2 + AC*AD*AF*AK*bx*cy**2 - 
     -    AC*AD*AM*AQ*bx*cy**2 - AC*AF*AO*AQ*bx*cy**2 + 
     -    AD*AF*AQ*BB*bx*cy**2 - AD**2*AQ*BD*bx*cy**2 + 
     -    AC*AF**2*AK*bz*cy**2 + AF**2*AQ*BB*bz*cy**2 - 
     -    AD*AF*AQ*BD*bz*cy**2 + AC*AF**2*AQ*ay*cz + 
     -    AC*AE*AF*AQ*az*cz - AC*AF**2*AK*ax*bx*cz - 
     -    AC*AD*AF*AK*az*bx*cz + AC*AD*AM*AQ*az*bx*cz + 
     -    AC*AF*AO*AQ*az*bx*cz - AF**2*AQ*ax*BB*bx*cz - 
     -    AD*AF*AQ*az*BB*bx*cz + AD*AF*AQ*ax*BD*bx*cz + 
     -    AD**2*AQ*az*BD*bx*cz + AC*AF**2*AK*bx**2*cz + 
     -    AF**2*AQ*BB*bx**2*cz - AD*AF*AQ*BD*bx**2*cz - 
     -    2*AC*AF**2*AQ*by*cz - AC*AF**2*AK*ay*by*cz - 
     -    AC*AE*AF*AK*az*by*cz + AC*AE*AM*AQ*az*by*cz - 
     -    AF**2*AQ*ay*BB*by*cz - AE*AF*AQ*az*BB*by*cz + 
     -    AD*AF*AQ*ay*BD*by*cz + AD*AE*AQ*az*BD*by*cz + 
     -    AC*AF**2*AK*by**2*cz + AF**2*AQ*BB*by**2*cz - 
     -    AD*AF*AQ*BD*by**2*cz + AC*AE*AF*AQ*bz*cz + 
     -    2*AC*AD*AF*AK*ax*bz*cz - 2*AC*AD*AM*AQ*ax*bz*cz - 
     -    2*AC*AF*AO*AQ*ax*bz*cz + 2*AC*AE*AF*AK*ay*bz*cz - 
     -    2*AC*AE*AM*AQ*ay*bz*cz + 2*AD*AF*AQ*ax*BB*bz*cz + 
     -    2*AE*AF*AQ*ay*BB*bz*cz - 2*AD**2*AQ*ax*BD*bz*cz - 
     -    2*AD*AE*AQ*ay*BD*bz*cz - AC*AD*AF*AK*bx*bz*cz + 
     -    AC*AD*AM*AQ*bx*bz*cz + AC*AF*AO*AQ*bx*bz*cz - 
     -    AD*AF*AQ*BB*bx*bz*cz + AD**2*AQ*BD*bx*bz*cz - 
     -    AC*AE*AF*AK*by*bz*cz + AC*AE*AM*AQ*by*bz*cz - 
     -    AE*AF*AQ*BB*by*bz*cz + AD*AE*AQ*BD*by*bz*cz + 
     -    AC*AF**2*AK*ax*cx*cz + AC*AD*AF*AK*az*cx*cz - 
     -    AC*AD*AM*AQ*az*cx*cz - AC*AF*AO*AQ*az*cx*cz + 
     -    AF**2*AQ*ax*BB*cx*cz + AD*AF*AQ*az*BB*cx*cz - 
     -    AD*AF*AQ*ax*BD*cx*cz - AD**2*AQ*az*BD*cx*cz - 
     -    AC*AF**2*AK*bx*cx*cz - AF**2*AQ*BB*bx*cx*cz + 
     -    AD*AF*AQ*BD*bx*cx*cz - AC*AD*AF*AK*bz*cx*cz + 
     -    AC*AD*AM*AQ*bz*cx*cz + AC*AF*AO*AQ*bz*cx*cz - 
     -    AD*AF*AQ*BB*bz*cx*cz + AD**2*AQ*BD*bz*cx*cz + 
     -    AC*AF**2*AQ*cy*cz + AC*AF**2*AK*ay*cy*cz + 
     -    AC*AE*AF*AK*az*cy*cz - AC*AE*AM*AQ*az*cy*cz + 
     -    AF**2*AQ*ay*BB*cy*cz + AE*AF*AQ*az*BB*cy*cz - 
     -    AD*AF*AQ*ay*BD*cy*cz - AD*AE*AQ*az*BD*cy*cz - 
     -    AC*AF**2*AK*by*cy*cz - AF**2*AQ*BB*by*cy*cz + 
     -    AD*AF*AQ*BD*by*cy*cz - AC*AE*AF*AK*bz*cy*cz + 
     -    AC*AE*AM*AQ*bz*cy*cz - AE*AF*AQ*BB*bz*cy*cz + 
     -    AD*AE*AQ*BD*bz*cy*cz - AC*AE*AF*AQ*cz**2 - 
     -    AC*AD*AF*AK*ax*cz**2 + AC*AD*AM*AQ*ax*cz**2 + 
     -    AC*AF*AO*AQ*ax*cz**2 - AC*AE*AF*AK*ay*cz**2 + 
     -    AC*AE*AM*AQ*ay*cz**2 - AD*AF*AQ*ax*BB*cz**2 - 
     -    AE*AF*AQ*ay*BB*cz**2 + AD**2*AQ*ax*BD*cz**2 + 
     -    AD*AE*AQ*ay*BD*cz**2 + AC*AD*AF*AK*bx*cz**2 - 
     -    AC*AD*AM*AQ*bx*cz**2 - AC*AF*AO*AQ*bx*cz**2 + 
     -    AD*AF*AQ*BB*bx*cz**2 - AD**2*AQ*BD*bx*cz**2 + 
     -    AC*AE*AF*AK*by*cz**2 - AC*AE*AM*AQ*by*cz**2 + 
     -    AE*AF*AQ*BB*by*cz**2 - AD*AE*AQ*BD*by*cz**2 - 
     -    AB*AE*(AF*AK*az*bx**2 + AM*AQ*az*bx**2 - 2*AF*AQ*az*by + 
     -       AO*AQ*ay*bx*by - AO*AQ*ax*by**2 + AF*AK*az*by**2 + 
     -       AM*AQ*az*by**2 + AF*AQ*ay*bz - AF*AK*ax*bx*bz - 
     -       AM*AQ*ax*bx*bz + AO*AQ*az*bx*bz - AF*AK*ay*by*bz - 
     -       AM*AQ*ay*by*bz - AO*AQ*ax*bz**2 - 2*AF*AK*az*bx*cx - 
     -       2*AM*AQ*az*bx*cx - AO*AQ*ay*by*cx + AO*AQ*by**2*cx + 
     -       AF*AK*ax*bz*cx + AM*AQ*ax*bz*cx - AO*AQ*az*bz*cx + 
     -       AF*AK*bx*bz*cx + AM*AQ*bx*bz*cx + AO*AQ*bz**2*cx + 
     -       AF*AK*az*cx**2 + AM*AQ*az*cx**2 - AF*AK*bz*cx**2 - 
     -       AM*AQ*bz*cx**2 + 2*AF*AQ*az*cy - AO*AQ*ay*bx*cy + 
     -       2*AO*AQ*ax*by*cy - 2*AF*AK*az*by*cy - 
     -       2*AM*AQ*az*by*cy - AO*AQ*bx*by*cy - AF*AQ*bz*cy + 
     -       AF*AK*ay*bz*cy + AM*AQ*ay*bz*cy + AF*AK*by*bz*cy + 
     -       AM*AQ*by*bz*cy + AO*AQ*ay*cx*cy - AO*AQ*by*cx*cy - 
     -       AO*AQ*ax*cy**2 + AF*AK*az*cy**2 + AM*AQ*az*cy**2 + 
     -       AO*AQ*bx*cy**2 - AF*AK*bz*cy**2 - AM*AQ*bz*cy**2 + 
     -       AD*(AK*BA - AQ*
     -           (2*AK*ax - AJ*ay - 2*by*cx + bx*cy + cx*cy)) - 
     -       AF*AQ*ay*cz + AF*AK*ax*bx*cz + AM*AQ*ax*bx*cz - 
     -       AO*AQ*az*bx*cz - AF*AK*bx**2*cz - AM*AQ*bx**2*cz + 
     -       2*AF*AQ*by*cz + AF*AK*ay*by*cz + AM*AQ*ay*by*cz - 
     -       AF*AK*by**2*cz - AM*AQ*by**2*cz + 2*AO*AQ*ax*bz*cz - 
     -       AO*AQ*bx*bz*cz - AF*AK*ax*cx*cz - AM*AQ*ax*cx*cz + 
     -       AO*AQ*az*cx*cz + AF*AK*bx*cx*cz + AM*AQ*bx*cx*cz - 
     -       AO*AQ*bz*cx*cz - AF*AQ*cy*cz - AF*AK*ay*cy*cz - 
     -       AM*AQ*ay*cy*cz + AF*AK*by*cy*cz + AM*AQ*by*cy*cz - 
     -       AO*AQ*ax*cz**2 + AO*AQ*bx*cz**2 + 
     -       AE*(AK*AV + AQ*
     -           (AJ*ax + AL*az - bx*cx + cx**2 - bz*cz + cz**2))) + 
     -    AA*(-(AO*AQ*(AE*AV + AF*AW)) - 
     -       AD**2*(AK*BA - 
     -          AQ*(2*AK*ax - AJ*ay - 2*by*cx + bx*cy + cx*cy)) - 
     -       AD*(AM*AQ*AW + 
     -          AF*(AK*AW + 
     -             AQ*(AL*ay - 2*AK*az - bz*cy + 2*by*cz - cy*cz)) + 
     -          AE*(AK*AV + 
     -             AQ*(AJ*ax + AL*az - bx*cx + cx**2 - bz*cz + cz**2)
     -             ))))/
     -  (AQ**1.5*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))

      dbz= (-(AC*AF**2*AQ*ax*bx) - AC*AD*AF*AQ*az*bx - 
     -    AC*AE*AF*AL*ay*bx**2 + AC*AF*AM*AQ*ay*bx**2 - 
     -    AC*AF**2*AL*az*bx**2 + AE**2*AQ*ay*BB*bx**2 + 
     -    AE*AF*AQ*az*BB*bx**2 - AD*AE*AQ*ay*BC*bx**2 - 
     -    AD*AF*AQ*az*BC*bx**2 - AC*AF**2*AQ*ay*by - 
     -    AC*AE*AF*AQ*az*by + AC*AE*AF*AL*ax*bx*by - 
     -    AC*AF*AM*AQ*ax*bx*by + AC*AD*AF*AL*ay*bx*by + 
     -    AC*AF*AN*AQ*ay*bx*by - AE**2*AQ*ax*BB*bx*by - 
     -    AD*AE*AQ*ay*BB*bx*by + AD*AE*AQ*ax*BC*bx*by + 
     -    AD**2*AQ*ay*BC*bx*by - AC*AD*AF*AL*ax*by**2 - 
     -    AC*AF*AN*AQ*ax*by**2 - AC*AF**2*AL*az*by**2 + 
     -    AD*AE*AQ*ax*BB*by**2 + AE*AF*AQ*az*BB*by**2 - 
     -    AD**2*AQ*ax*BC*by**2 - AD*AF*AQ*az*BC*by**2 + 
     -    2*AC*AD*AF*AQ*ax*bz + 2*AC*AE*AF*AQ*ay*bz + 
     -    AC*AF**2*AL*ax*bx*bz + AC*AD*AF*AL*az*bx*bz + 
     -    AC*AF*AN*AQ*az*bx*bz - AE*AF*AQ*ax*BB*bx*bz - 
     -    AD*AE*AQ*az*BB*bx*bz + AD*AF*AQ*ax*BC*bx*bz + 
     -    AD**2*AQ*az*BC*bx*bz + AC*AF**2*AL*ay*by*bz + 
     -    AC*AE*AF*AL*az*by*bz - AC*AF*AM*AQ*az*by*bz - 
     -    AE*AF*AQ*ay*BB*by*bz - AE**2*AQ*az*BB*by*bz + 
     -    AD*AF*AQ*ay*BC*by*bz + AD*AE*AQ*az*BC*by*bz - 
     -    AC*AD*AF*AL*ax*bz**2 - AC*AF*AN*AQ*ax*bz**2 - 
     -    AC*AE*AF*AL*ay*bz**2 + AC*AF*AM*AQ*ay*bz**2 + 
     -    AD*AE*AQ*ax*BB*bz**2 + AE**2*AQ*ay*BB*bz**2 - 
     -    AD**2*AQ*ax*BC*bz**2 - AD*AE*AQ*ay*BC*bz**2 + 
     -    AC*AF**2*AQ*ax*cx + AC*AD*AF*AQ*az*cx + 
     -    AC*AF**2*AQ*bx*cx + 2*AC*AE*AF*AL*ay*bx*cx - 
     -    2*AC*AF*AM*AQ*ay*bx*cx + 2*AC*AF**2*AL*az*bx*cx - 
     -    2*AE**2*AQ*ay*BB*bx*cx - 2*AE*AF*AQ*az*BB*bx*cx + 
     -    2*AD*AE*AQ*ay*BC*bx*cx + 2*AD*AF*AQ*az*BC*bx*cx - 
     -    AC*AE*AF*AL*ax*by*cx + AC*AF*AM*AQ*ax*by*cx - 
     -    AC*AD*AF*AL*ay*by*cx - AC*AF*AN*AQ*ay*by*cx + 
     -    AE**2*AQ*ax*BB*by*cx + AD*AE*AQ*ay*BB*by*cx - 
     -    AD*AE*AQ*ax*BC*by*cx - AD**2*AQ*ay*BC*by*cx - 
     -    AC*AE*AF*AL*bx*by*cx + AC*AF*AM*AQ*bx*by*cx + 
     -    AE**2*AQ*BB*bx*by*cx - AD*AE*AQ*BC*bx*by*cx + 
     -    AC*AD*AF*AL*by**2*cx + AC*AF*AN*AQ*by**2*cx - 
     -    AD*AE*AQ*BB*by**2*cx + AD**2*AQ*BC*by**2*cx - 
     -    2*AC*AD*AF*AQ*bz*cx - AC*AF**2*AL*ax*bz*cx - 
     -    AC*AD*AF*AL*az*bz*cx - AC*AF*AN*AQ*az*bz*cx + 
     -    AE*AF*AQ*ax*BB*bz*cx + AD*AE*AQ*az*BB*bz*cx - 
     -    AD*AF*AQ*ax*BC*bz*cx - AD**2*AQ*az*BC*bz*cx - 
     -    AC*AF**2*AL*bx*bz*cx + AE*AF*AQ*BB*bx*bz*cx - 
     -    AD*AF*AQ*BC*bx*bz*cx + AC*AD*AF*AL*bz**2*cx + 
     -    AC*AF*AN*AQ*bz**2*cx - AD*AE*AQ*BB*bz**2*cx + 
     -    AD**2*AQ*BC*bz**2*cx - AC*AF**2*AQ*cx**2 - 
     -    AC*AE*AF*AL*ay*cx**2 + AC*AF*AM*AQ*ay*cx**2 - 
     -    AC*AF**2*AL*az*cx**2 + AE**2*AQ*ay*BB*cx**2 + 
     -    AE*AF*AQ*az*BB*cx**2 - AD*AE*AQ*ay*BC*cx**2 - 
     -    AD*AF*AQ*az*BC*cx**2 + AC*AE*AF*AL*by*cx**2 - 
     -    AC*AF*AM*AQ*by*cx**2 - AE**2*AQ*BB*by*cx**2 + 
     -    AD*AE*AQ*BC*by*cx**2 + AC*AF**2*AL*bz*cx**2 - 
     -    AE*AF*AQ*BB*bz*cx**2 + AD*AF*AQ*BC*bz*cx**2 + 
     -    AC*AF**2*AQ*ay*cy + AC*AE*AF*AQ*az*cy - 
     -    AC*AE*AF*AL*ax*bx*cy + AC*AF*AM*AQ*ax*bx*cy - 
     -    AC*AD*AF*AL*ay*bx*cy - AC*AF*AN*AQ*ay*bx*cy + 
     -    AE**2*AQ*ax*BB*bx*cy + AD*AE*AQ*ay*BB*bx*cy - 
     -    AD*AE*AQ*ax*BC*bx*cy - AD**2*AQ*ay*BC*bx*cy + 
     -    AC*AE*AF*AL*bx**2*cy - AC*AF*AM*AQ*bx**2*cy - 
     -    AE**2*AQ*BB*bx**2*cy + AD*AE*AQ*BC*bx**2*cy + 
     -    AC*AF**2*AQ*by*cy + 2*AC*AD*AF*AL*ax*by*cy + 
     -    2*AC*AF*AN*AQ*ax*by*cy + 2*AC*AF**2*AL*az*by*cy - 
     -    2*AD*AE*AQ*ax*BB*by*cy - 2*AE*AF*AQ*az*BB*by*cy + 
     -    2*AD**2*AQ*ax*BC*by*cy + 2*AD*AF*AQ*az*BC*by*cy - 
     -    AC*AD*AF*AL*bx*by*cy - AC*AF*AN*AQ*bx*by*cy + 
     -    AD*AE*AQ*BB*bx*by*cy - AD**2*AQ*BC*bx*by*cy - 
     -    2*AC*AE*AF*AQ*bz*cy - AC*AF**2*AL*ay*bz*cy - 
     -    AC*AE*AF*AL*az*bz*cy + AC*AF*AM*AQ*az*bz*cy + 
     -    AE*AF*AQ*ay*BB*bz*cy + AE**2*AQ*az*BB*bz*cy - 
     -    AD*AF*AQ*ay*BC*bz*cy - AD*AE*AQ*az*BC*bz*cy - 
     -    AC*AF**2*AL*by*bz*cy + AE*AF*AQ*BB*by*bz*cy - 
     -    AD*AF*AQ*BC*by*bz*cy + AC*AE*AF*AL*bz**2*cy - 
     -    AC*AF*AM*AQ*bz**2*cy - AE**2*AQ*BB*bz**2*cy + 
     -    AD*AE*AQ*BC*bz**2*cy + AC*AE*AF*AL*ax*cx*cy - 
     -    AC*AF*AM*AQ*ax*cx*cy + AC*AD*AF*AL*ay*cx*cy + 
     -    AC*AF*AN*AQ*ay*cx*cy - AE**2*AQ*ax*BB*cx*cy - 
     -    AD*AE*AQ*ay*BB*cx*cy + AD*AE*AQ*ax*BC*cx*cy + 
     -    AD**2*AQ*ay*BC*cx*cy - AC*AE*AF*AL*bx*cx*cy + 
     -    AC*AF*AM*AQ*bx*cx*cy + AE**2*AQ*BB*bx*cx*cy - 
     -    AD*AE*AQ*BC*bx*cx*cy - AC*AD*AF*AL*by*cx*cy - 
     -    AC*AF*AN*AQ*by*cx*cy + AD*AE*AQ*BB*by*cx*cy - 
     -    AD**2*AQ*BC*by*cx*cy - AC*AF**2*AQ*cy**2 - 
     -    AC*AD*AF*AL*ax*cy**2 - AC*AF*AN*AQ*ax*cy**2 - 
     -    AC*AF**2*AL*az*cy**2 + AD*AE*AQ*ax*BB*cy**2 + 
     -    AE*AF*AQ*az*BB*cy**2 - AD**2*AQ*ax*BC*cy**2 - 
     -    AD*AF*AQ*az*BC*cy**2 + AC*AD*AF*AL*bx*cy**2 + 
     -    AC*AF*AN*AQ*bx*cy**2 - AD*AE*AQ*BB*bx*cy**2 + 
     -    AD**2*AQ*BC*bx*cy**2 + AC*AF**2*AL*bz*cy**2 - 
     -    AE*AF*AQ*BB*bz*cy**2 + AD*AF*AQ*BC*bz*cy**2 - 
     -    2*AC*AD*AF*AQ*ax*cz - 2*AC*AE*AF*AQ*ay*cz + 
     -    AC*AD*AF*AQ*bx*cz - AC*AF**2*AL*ax*bx*cz - 
     -    AC*AD*AF*AL*az*bx*cz - AC*AF*AN*AQ*az*bx*cz + 
     -    AE*AF*AQ*ax*BB*bx*cz + AD*AE*AQ*az*BB*bx*cz - 
     -    AD*AF*AQ*ax*BC*bx*cz - AD**2*AQ*az*BC*bx*cz + 
     -    AC*AF**2*AL*bx**2*cz - AE*AF*AQ*BB*bx**2*cz + 
     -    AD*AF*AQ*BC*bx**2*cz + AC*AE*AF*AQ*by*cz - 
     -    AC*AF**2*AL*ay*by*cz - AC*AE*AF*AL*az*by*cz + 
     -    AC*AF*AM*AQ*az*by*cz + AE*AF*AQ*ay*BB*by*cz + 
     -    AE**2*AQ*az*BB*by*cz - AD*AF*AQ*ay*BC*by*cz - 
     -    AD*AE*AQ*az*BC*by*cz + AC*AF**2*AL*by**2*cz - 
     -    AE*AF*AQ*BB*by**2*cz + AD*AF*AQ*BC*by**2*cz + 
     -    2*AC*AD*AF*AL*ax*bz*cz + 2*AC*AF*AN*AQ*ax*bz*cz + 
     -    2*AC*AE*AF*AL*ay*bz*cz - 2*AC*AF*AM*AQ*ay*bz*cz - 
     -    2*AD*AE*AQ*ax*BB*bz*cz - 2*AE**2*AQ*ay*BB*bz*cz + 
     -    2*AD**2*AQ*ax*BC*bz*cz + 2*AD*AE*AQ*ay*BC*bz*cz - 
     -    AC*AD*AF*AL*bx*bz*cz - AC*AF*AN*AQ*bx*bz*cz + 
     -    AD*AE*AQ*BB*bx*bz*cz - AD**2*AQ*BC*bx*bz*cz - 
     -    AC*AE*AF*AL*by*bz*cz + AC*AF*AM*AQ*by*bz*cz + 
     -    AE**2*AQ*BB*by*bz*cz - AD*AE*AQ*BC*by*bz*cz + 
     -    AC*AD*AF*AQ*cx*cz + AC*AF**2*AL*ax*cx*cz + 
     -    AC*AD*AF*AL*az*cx*cz + AC*AF*AN*AQ*az*cx*cz - 
     -    AE*AF*AQ*ax*BB*cx*cz - AD*AE*AQ*az*BB*cx*cz + 
     -    AD*AF*AQ*ax*BC*cx*cz + AD**2*AQ*az*BC*cx*cz - 
     -    AC*AF**2*AL*bx*cx*cz + AE*AF*AQ*BB*bx*cx*cz - 
     -    AD*AF*AQ*BC*bx*cx*cz - AC*AD*AF*AL*bz*cx*cz - 
     -    AC*AF*AN*AQ*bz*cx*cz + AD*AE*AQ*BB*bz*cx*cz - 
     -    AD**2*AQ*BC*bz*cx*cz + AC*AE*AF*AQ*cy*cz + 
     -    AC*AF**2*AL*ay*cy*cz + AC*AE*AF*AL*az*cy*cz - 
     -    AC*AF*AM*AQ*az*cy*cz - AE*AF*AQ*ay*BB*cy*cz - 
     -    AE**2*AQ*az*BB*cy*cz + AD*AF*AQ*ay*BC*cy*cz + 
     -    AD*AE*AQ*az*BC*cy*cz - AC*AF**2*AL*by*cy*cz + 
     -    AE*AF*AQ*BB*by*cy*cz - AD*AF*AQ*BC*by*cy*cz - 
     -    AC*AE*AF*AL*bz*cy*cz + AC*AF*AM*AQ*bz*cy*cz + 
     -    AE**2*AQ*BB*bz*cy*cz - AD*AE*AQ*BC*bz*cy*cz - 
     -    AC*AD*AF*AL*ax*cz**2 - AC*AF*AN*AQ*ax*cz**2 - 
     -    AC*AE*AF*AL*ay*cz**2 + AC*AF*AM*AQ*ay*cz**2 + 
     -    AD*AE*AQ*ax*BB*cz**2 + AE**2*AQ*ay*BB*cz**2 - 
     -    AD**2*AQ*ax*BC*cz**2 - AD*AE*AQ*ay*BC*cz**2 + 
     -    AC*AD*AF*AL*bx*cz**2 + AC*AF*AN*AQ*bx*cz**2 - 
     -    AD*AE*AQ*BB*bx*cz**2 + AD**2*AQ*BC*bx*cz**2 + 
     -    AC*AE*AF*AL*by*cz**2 - AC*AF*AM*AQ*by*cz**2 - 
     -    AE**2*AQ*BB*by*cz**2 + AD*AE*AQ*BC*by*cz**2 - 
     -    AB*(AM*AQ*(AF*AW + AD*BA) + 
     -       AE**2*(AL*AV - 
     -          AQ*(2*AL*ay - AK*az - 2*bz*cy + by*cz + cy*cz))
     -        + AE*(AN*AQ*BA + 
     -          AF*(AL*AW + 
     -             AQ*(AJ*ax + AK*ay - bx*cx + cx**2 - by*cy + 
     -                cy**2)) + 
     -          AD*(AL*BA - 
     -             AQ*(2*AL*ax - AJ*az - 2*bz*cx + bx*cz + cx*cz)
     -             ))) + AA*
     -     (AN*AQ*(AE*AV + AF*AW) - 
     -       AD**2*(AL*BA - 
     -          AQ*(2*AL*ax - AJ*az - 2*bz*cx + bx*cz + cx*cz))
     -        - AD*(-(AM*AQ*AV) + 
     -          AF*(AL*AW + 
     -             AQ*(AJ*ax + AK*ay - bx*cx + cx**2 - by*cy + 
     -                cy**2)) + 
     -          AE*(AL*AV - 
     -             AQ*(2*AL*ay - AK*az - 2*bz*cy + by*cz + cy*cz)))))/
     -  (AQ**1.5*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))


      dcx=(-2*AC*AE*AF*AQ*ay*bx - 2*AC*AF**2*AQ*az*bx + 
     -    AC*AE*AF*AJ*ay*bx**2 - AE*AF*AH*AQ*ay*bx**2 + 
     -    AE**2*AI*AQ*ay*bx**2 + AC*AF**2*AJ*az*bx**2 - 
     -    AF**2*AH*AQ*az*bx**2 + AE*AF*AI*AQ*az*bx**2 + 
     -    AC*AE*AF*AQ*ax*by + AC*AD*AF*AQ*ay*by + 
     -    AC*AE*AF*AQ*bx*by - AC*AE*AF*AJ*ax*bx*by + 
     -    AE*AF*AH*AQ*ax*bx*by - AE**2*AI*AQ*ax*bx*by - 
     -    AC*AD*AF*AJ*ay*bx*by + AD*AF*AH*AQ*ay*bx*by - 
     -    AD*AE*AI*AQ*ay*bx*by + AC*AE*AQ*ay*bx**2*by - 
     -    AC*AD*AF*AQ*by**2 + AC*AD*AF*AJ*ax*by**2 - 
     -    AD*AF*AH*AQ*ax*by**2 + AD*AE*AI*AQ*ax*by**2 + 
     -    AC*AF**2*AJ*az*by**2 - AF**2*AH*AQ*az*by**2 + 
     -    AE*AF*AI*AQ*az*by**2 - AC*AE*AQ*ax*bx*by**2 - 
     -    AC*AD*AQ*ay*bx*by**2 + AC*AD*AQ*ax*by**3 + 
     -    AC*AF**2*AQ*ax*bz + AC*AD*AF*AQ*az*bz + 
     -    AC*AF**2*AQ*bx*bz - AC*AF**2*AJ*ax*bx*bz + 
     -    AF**2*AH*AQ*ax*bx*bz - AE*AF*AI*AQ*ax*bx*bz - 
     -    AC*AD*AF*AJ*az*bx*bz + AD*AF*AH*AQ*az*bx*bz - 
     -    AD*AE*AI*AQ*az*bx*bz + AC*AF*AQ*ay*bx**2*bz - 
     -    AC*AF**2*AJ*ay*by*bz + AF**2*AH*AQ*ay*by*bz - 
     -    AE*AF*AI*AQ*ay*by*bz - AC*AE*AF*AJ*az*by*bz + 
     -    AE*AF*AH*AQ*az*by*bz - AE**2*AI*AQ*az*by*bz - 
     -    AC*AF*AQ*ax*bx*by*bz - AC*AD*AQ*az*bx*by*bz - 
     -    AC*AE*AQ*az*by**2*bz - AC*AD*AF*AQ*bz**2 + 
     -    AC*AD*AF*AJ*ax*bz**2 - AD*AF*AH*AQ*ax*bz**2 + 
     -    AD*AE*AI*AQ*ax*bz**2 + AC*AE*AF*AJ*ay*bz**2 - 
     -    AE*AF*AH*AQ*ay*bz**2 + AE**2*AI*AQ*ay*bz**2 + 
     -    AC*AD*AQ*ax*by*bz**2 + AC*AE*AQ*ay*by*bz**2 - 
     -    AC*AF*AQ*az*by*bz**2 + AC*AF*AQ*ay*bz**3 + 
     -    2*AC*AE*AF*AQ*ay*cx + 2*AC*AF**2*AQ*az*cx - 
     -    2*AC*AE*AF*AJ*ay*bx*cx + 2*AE*AF*AH*AQ*ay*bx*cx - 
     -    2*AE**2*AI*AQ*ay*bx*cx - 2*AC*AF**2*AJ*az*bx*cx + 
     -    2*AF**2*AH*AQ*az*bx*cx - 2*AE*AF*AI*AQ*az*bx*cx - 
     -    2*AC*AE*AF*AQ*by*cx + AC*AE*AF*AJ*ax*by*cx - 
     -    AE*AF*AH*AQ*ax*by*cx + AE**2*AI*AQ*ax*by*cx + 
     -    AC*AD*AF*AJ*ay*by*cx - AD*AF*AH*AQ*ay*by*cx + 
     -    AD*AE*AI*AQ*ay*by*cx + AC*AE*AF*AJ*bx*by*cx - 
     -    AE*AF*AH*AQ*bx*by*cx + AE**2*AI*AQ*bx*by*cx - 
     -    2*AC*AE*AQ*ay*bx*by*cx - AC*AD*AF*AJ*by**2*cx + 
     -    AD*AF*AH*AQ*by**2*cx - AD*AE*AI*AQ*by**2*cx + 
     -    AC*AE*AQ*ax*by**2*cx + AC*AD*AQ*ay*by**2*cx + 
     -    AC*AE*AQ*bx*by**2*cx - AC*AD*AQ*by**3*cx - 
     -    2*AC*AF**2*AQ*bz*cx + AC*AF**2*AJ*ax*bz*cx - 
     -    AF**2*AH*AQ*ax*bz*cx + AE*AF*AI*AQ*ax*bz*cx + 
     -    AC*AD*AF*AJ*az*bz*cx - AD*AF*AH*AQ*az*bz*cx + 
     -    AD*AE*AI*AQ*az*bz*cx + AC*AF**2*AJ*bx*bz*cx - 
     -    AF**2*AH*AQ*bx*bz*cx + AE*AF*AI*AQ*bx*bz*cx - 
     -    2*AC*AF*AQ*ay*bx*bz*cx + AC*AF*AQ*ax*by*bz*cx + 
     -    AC*AD*AQ*az*by*bz*cx + AC*AF*AQ*bx*by*bz*cx - 
     -    AC*AD*AF*AJ*bz**2*cx + AD*AF*AH*AQ*bz**2*cx - 
     -    AD*AE*AI*AQ*bz**2*cx - AC*AD*AQ*by*bz**2*cx + 
     -    AC*AE*AF*AJ*ay*cx**2 - AE*AF*AH*AQ*ay*cx**2 + 
     -    AE**2*AI*AQ*ay*cx**2 + AC*AF**2*AJ*az*cx**2 - 
     -    AF**2*AH*AQ*az*cx**2 + AE*AF*AI*AQ*az*cx**2 - 
     -    AC*AE*AF*AJ*by*cx**2 + AE*AF*AH*AQ*by*cx**2 - 
     -    AE**2*AI*AQ*by*cx**2 + AC*AE*AQ*ay*by*cx**2 - 
     -    AC*AE*AQ*by**2*cx**2 - AC*AF**2*AJ*bz*cx**2 + 
     -    AF**2*AH*AQ*bz*cx**2 - AE*AF*AI*AQ*bz*cx**2 + 
     -    AC*AF*AQ*ay*bz*cx**2 - AC*AF*AQ*by*bz*cx**2 - 
     -    AC*AE*AF*AQ*ax*cy - AC*AD*AF*AQ*ay*cy + 
     -    AC*AE*AF*AQ*bx*cy + AC*AE*AF*AJ*ax*bx*cy - 
     -    AE*AF*AH*AQ*ax*bx*cy + AE**2*AI*AQ*ax*bx*cy + 
     -    AC*AD*AF*AJ*ay*bx*cy - AD*AF*AH*AQ*ay*bx*cy + 
     -    AD*AE*AI*AQ*ay*bx*cy - AC*AE*AF*AJ*bx**2*cy + 
     -    AE*AF*AH*AQ*bx**2*cy - AE**2*AI*AQ*bx**2*cy + 
     -    AC*AD*AF*AQ*by*cy - 2*AC*AD*AF*AJ*ax*by*cy + 
     -    2*AD*AF*AH*AQ*ax*by*cy - 2*AD*AE*AI*AQ*ax*by*cy - 
     -    2*AC*AF**2*AJ*az*by*cy + 2*AF**2*AH*AQ*az*by*cy - 
     -    2*AE*AF*AI*AQ*az*by*cy + AC*AD*AF*AJ*bx*by*cy - 
     -    AD*AF*AH*AQ*bx*by*cy + AD*AE*AI*AQ*bx*by*cy + 
     -    AC*AE*AQ*ax*bx*by*cy + AC*AD*AQ*ay*bx*by*cy - 
     -    AC*AE*AQ*bx**2*by*cy - 2*AC*AD*AQ*ax*by**2*cy + 
     -    AC*AD*AQ*bx*by**2*cy + AC*AF**2*AJ*ay*bz*cy - 
     -    AF**2*AH*AQ*ay*bz*cy + AE*AF*AI*AQ*ay*bz*cy + 
     -    AC*AE*AF*AJ*az*bz*cy - AE*AF*AH*AQ*az*bz*cy + 
     -    AE**2*AI*AQ*az*bz*cy + AC*AF*AQ*ax*bx*bz*cy - 
     -    AC*AF*AQ*bx**2*bz*cy + AC*AF**2*AJ*by*bz*cy - 
     -    AF**2*AH*AQ*by*bz*cy + AE*AF*AI*AQ*by*bz*cy + 
     -    AC*AE*AQ*az*by*bz*cy - AC*AE*AF*AJ*bz**2*cy + 
     -    AE*AF*AH*AQ*bz**2*cy - AE**2*AI*AQ*bz**2*cy + 
     -    AC*AF*AQ*az*bz**2*cy - AC*AE*AQ*by*bz**2*cy - 
     -    AC*AF*AQ*bz**3*cy - AC*AE*AF*AJ*ax*cx*cy + 
     -    AE*AF*AH*AQ*ax*cx*cy - AE**2*AI*AQ*ax*cx*cy - 
     -    AC*AD*AF*AJ*ay*cx*cy + AD*AF*AH*AQ*ay*cx*cy - 
     -    AD*AE*AI*AQ*ay*cx*cy + AC*AE*AF*AJ*bx*cx*cy - 
     -    AE*AF*AH*AQ*bx*cx*cy + AE**2*AI*AQ*bx*cx*cy + 
     -    AC*AD*AF*AJ*by*cx*cy - AD*AF*AH*AQ*by*cx*cy + 
     -    AD*AE*AI*AQ*by*cx*cy - AC*AE*AQ*ax*by*cx*cy - 
     -    AC*AD*AQ*ay*by*cx*cy + AC*AE*AQ*bx*by*cx*cy + 
     -    AC*AD*AQ*by**2*cx*cy - AC*AF*AQ*ax*bz*cx*cy + 
     -    AC*AF*AQ*bx*bz*cx*cy + AC*AD*AF*AJ*ax*cy**2 - 
     -    AD*AF*AH*AQ*ax*cy**2 + AD*AE*AI*AQ*ax*cy**2 + 
     -    AC*AF**2*AJ*az*cy**2 - AF**2*AH*AQ*az*cy**2 + 
     -    AE*AF*AI*AQ*az*cy**2 - AC*AD*AF*AJ*bx*cy**2 + 
     -    AD*AF*AH*AQ*bx*cy**2 - AD*AE*AI*AQ*bx*cy**2 + 
     -    AC*AD*AQ*ax*by*cy**2 - AC*AD*AQ*bx*by*cy**2 - 
     -    AC*AF**2*AJ*bz*cy**2 + AF**2*AH*AQ*bz*cy**2 - 
     -    AE*AF*AI*AQ*bz*cy**2 - AC*AF**2*AQ*ax*cz - 
     -    AC*AD*AF*AQ*az*cz + AC*AF**2*AQ*bx*cz + 
     -    AC*AF**2*AJ*ax*bx*cz - AF**2*AH*AQ*ax*bx*cz + 
     -    AE*AF*AI*AQ*ax*bx*cz + AC*AD*AF*AJ*az*bx*cz - 
     -    AD*AF*AH*AQ*az*bx*cz + AD*AE*AI*AQ*az*bx*cz - 
     -    AC*AF**2*AJ*bx**2*cz + AF**2*AH*AQ*bx**2*cz - 
     -    AE*AF*AI*AQ*bx**2*cz + AC*AF**2*AJ*ay*by*cz - 
     -    AF**2*AH*AQ*ay*by*cz + AE*AF*AI*AQ*ay*by*cz + 
     -    AC*AE*AF*AJ*az*by*cz - AE*AF*AH*AQ*az*by*cz + 
     -    AE**2*AI*AQ*az*by*cz + AC*AD*AQ*az*bx*by*cz - 
     -    AC*AF**2*AJ*by**2*cz + AF**2*AH*AQ*by**2*cz - 
     -    AE*AF*AI*AQ*by**2*cz + AC*AE*AQ*az*by**2*cz + 
     -    AC*AD*AF*AQ*bz*cz - 2*AC*AD*AF*AJ*ax*bz*cz + 
     -    2*AD*AF*AH*AQ*ax*bz*cz - 2*AD*AE*AI*AQ*ax*bz*cz - 
     -    2*AC*AE*AF*AJ*ay*bz*cz + 2*AE*AF*AH*AQ*ay*bz*cz - 
     -    2*AE**2*AI*AQ*ay*bz*cz + AC*AD*AF*AJ*bx*bz*cz - 
     -    AD*AF*AH*AQ*bx*bz*cz + AD*AE*AI*AQ*bx*bz*cz + 
     -    AC*AE*AF*AJ*by*bz*cz - AE*AF*AH*AQ*by*bz*cz + 
     -    AE**2*AI*AQ*by*bz*cz - 2*AC*AD*AQ*ax*by*bz*cz - 
     -    2*AC*AE*AQ*ay*by*bz*cz + AC*AF*AQ*az*by*bz*cz + 
     -    AC*AD*AQ*bx*by*bz*cz + AC*AE*AQ*by**2*bz*cz - 
     -    2*AC*AF*AQ*ay*bz**2*cz + AC*AF*AQ*by*bz**2*cz - 
     -    AC*AF**2*AJ*ax*cx*cz + AF**2*AH*AQ*ax*cx*cz - 
     -    AE*AF*AI*AQ*ax*cx*cz - AC*AD*AF*AJ*az*cx*cz + 
     -    AD*AF*AH*AQ*az*cx*cz - AD*AE*AI*AQ*az*cx*cz + 
     -    AC*AF**2*AJ*bx*cx*cz - AF**2*AH*AQ*bx*cx*cz + 
     -    AE*AF*AI*AQ*bx*cx*cz - AC*AD*AQ*az*by*cx*cz + 
     -    AC*AD*AF*AJ*bz*cx*cz - AD*AF*AH*AQ*bz*cx*cz + 
     -    AD*AE*AI*AQ*bz*cx*cz + AC*AD*AQ*by*bz*cx*cz - 
     -    AC*AF**2*AJ*ay*cy*cz + AF**2*AH*AQ*ay*cy*cz - 
     -    AE*AF*AI*AQ*ay*cy*cz - AC*AE*AF*AJ*az*cy*cz + 
     -    AE*AF*AH*AQ*az*cy*cz - AE**2*AI*AQ*az*cy*cz + 
     -    AC*AF**2*AJ*by*cy*cz - AF**2*AH*AQ*by*cy*cz + 
     -    AE*AF*AI*AQ*by*cy*cz - AC*AE*AQ*az*by*cy*cz + 
     -    AC*AE*AF*AJ*bz*cy*cz - AE*AF*AH*AQ*bz*cy*cz + 
     -    AE**2*AI*AQ*bz*cy*cz - AC*AF*AQ*az*bz*cy*cz + 
     -    AC*AE*AQ*by*bz*cy*cz + AC*AF*AQ*bz**2*cy*cz + 
     -    AC*AD*AF*AJ*ax*cz**2 - AD*AF*AH*AQ*ax*cz**2 + 
     -    AD*AE*AI*AQ*ax*cz**2 + AC*AE*AF*AJ*ay*cz**2 - 
     -    AE*AF*AH*AQ*ay*cz**2 + AE**2*AI*AQ*ay*cz**2 - 
     -    AC*AD*AF*AJ*bx*cz**2 + AD*AF*AH*AQ*bx*cz**2 - 
     -    AD*AE*AI*AQ*bx*cz**2 - AC*AE*AF*AJ*by*cz**2 + 
     -    AE*AF*AH*AQ*by*cz**2 - AE**2*AI*AQ*by*cz**2 + 
     -    AC*AD*AQ*ax*by*cz**2 + AC*AE*AQ*ay*by*cz**2 - 
     -    AC*AD*AQ*bx*by*cz**2 - AC*AE*AQ*by**2*cz**2 + 
     -    AC*AF*AQ*ay*bz*cz**2 - AC*AF*AQ*by*bz*cz**2 - 
     -    AC*AE*AQ*ay*bx**2*dy + AC*AE*AQ*ax*bx*by*dy + 
     -    AC*AD*AQ*ay*bx*by*dy - AC*AD*AQ*ax*by**2*dy + 
     -    AC*AD*AQ*az*bx*bz*dy + AC*AE*AQ*az*by*bz*dy - 
     -    AC*AD*AQ*ax*bz**2*dy - AC*AE*AQ*ay*bz**2*dy + 
     -    2*AC*AE*AQ*ay*bx*cx*dy - AC*AE*AQ*ax*by*cx*dy - 
     -    AC*AD*AQ*ay*by*cx*dy - AC*AE*AQ*bx*by*cx*dy + 
     -    AC*AD*AQ*by**2*cx*dy - AC*AD*AQ*az*bz*cx*dy + 
     -    AC*AD*AQ*bz**2*cx*dy - AC*AE*AQ*ay*cx**2*dy + 
     -    AC*AE*AQ*by*cx**2*dy - AC*AE*AQ*ax*bx*cy*dy - 
     -    AC*AD*AQ*ay*bx*cy*dy + AC*AE*AQ*bx**2*cy*dy + 
     -    2*AC*AD*AQ*ax*by*cy*dy - AC*AD*AQ*bx*by*cy*dy - 
     -    AC*AE*AQ*az*bz*cy*dy + AC*AE*AQ*bz**2*cy*dy + 
     -    AC*AE*AQ*ax*cx*cy*dy + AC*AD*AQ*ay*cx*cy*dy - 
     -    AC*AE*AQ*bx*cx*cy*dy - AC*AD*AQ*by*cx*cy*dy - 
     -    AC*AD*AQ*ax*cy**2*dy + AC*AD*AQ*bx*cy**2*dy - 
     -    AC*AD*AQ*az*bx*cz*dy - AC*AE*AQ*az*by*cz*dy + 
     -    2*AC*AD*AQ*ax*bz*cz*dy + 2*AC*AE*AQ*ay*bz*cz*dy - 
     -    AC*AD*AQ*bx*bz*cz*dy - AC*AE*AQ*by*bz*cz*dy + 
     -    AC*AD*AQ*az*cx*cz*dy - AC*AD*AQ*bz*cx*cz*dy + 
     -    AC*AE*AQ*az*cy*cz*dy - AC*AE*AQ*bz*cy*cz*dy - 
     -    AC*AD*AQ*ax*cz**2*dy - AC*AE*AQ*ay*cz**2*dy + 
     -    AC*AD*AQ*bx*cz**2*dy + AC*AE*AQ*by*cz**2*dy + 
     -    AB*(AE**2*(AJ*AV + 
     -          AQ*(-2*ay*bx + ax*by + bx*by + 2*ay*cx - 
     -             2*by*cx - ax*cy + bx*cy)) + 
     -       AE*(AD*(AJ*BA + 
     -             AQ*(AI*AL + AK*ay + by*(-by + cy))) + 
     -          AF*(AJ*AW + 
     -             AQ*(-2*az*bx + ax*bz + bx*bz + 2*az*cx - 
     -                2*bz*cx - ax*cz + bx*cz)) - 
     -          AQ*AW*(by - dy)) - AQ*(AF*AW + AD*BA)*(bz - dz))
     -      - AC*AF*AQ*ay*bx**2*dz + AC*AF*AQ*ax*bx*by*dz + 
     -    AC*AF*AQ*az*by*bz*dz - AC*AF*AQ*ay*bz**2*dz + 
     -    2*AC*AF*AQ*ay*bx*cx*dz - AC*AF*AQ*ax*by*cx*dz - 
     -    AC*AF*AQ*bx*by*cx*dz - AC*AF*AQ*ay*cx**2*dz + 
     -    AC*AF*AQ*by*cx**2*dz - AC*AF*AQ*ax*bx*cy*dz + 
     -    AC*AF*AQ*bx**2*cy*dz - AC*AF*AQ*az*bz*cy*dz + 
     -    AC*AF*AQ*bz**2*cy*dz + AC*AF*AQ*ax*cx*cy*dz - 
     -    AC*AF*AQ*bx*cx*cy*dz - AC*AF*AQ*az*by*cz*dz + 
     -    2*AC*AF*AQ*ay*bz*cz*dz - AC*AF*AQ*by*bz*cz*dz + 
     -    AC*AF*AQ*az*cy*cz*dz - AC*AF*AQ*bz*cy*cz*dz - 
     -    AC*AF*AQ*ay*cz**2*dz + AC*AF*AQ*by*cz**2*dz + 
     -    AA*AD*(AD*AQ*ay*by - AD*AJ*ay*bx*by - 
     -       AQ*az*bx**2*by - AD*AQ*by**2 + AD*AJ*ax*by**2 - 
     -       AQ*az*by**3 + AD*AQ*az*bz - AD*AJ*az*bx*bz + 
     -       AQ*ay*bx**2*bz + AQ*ay*by**2*bz - AD*AQ*bz**2 + 
     -       AD*AJ*ax*bz**2 - AQ*az*by*bz**2 + AQ*ay*bz**3 + 
     -       AD*AJ*ay*by*cx + 2*AQ*az*bx*by*cx - 
     -       AD*AJ*by**2*cx + AD*AJ*az*bz*cx - 
     -       2*AQ*ay*bx*bz*cx - AD*AJ*bz**2*cx - 
     -       AQ*az*by*cx**2 + AQ*ay*bz*cx**2 - AD*AQ*ay*cy + 
     -       AD*AJ*ay*bx*cy + AD*AQ*by*cy - 2*AD*AJ*ax*by*cy + 
     -       AD*AJ*bx*by*cy + 2*AQ*az*by**2*cy + 
     -       AQ*ax*bx*bz*cy - AQ*bx**2*bz*cy - AQ*ay*by*bz*cy - 
     -       AQ*by**2*bz*cy + AQ*az*bz**2*cy - AQ*bz**3*cy - 
     -       AD*AJ*ay*cx*cy + AD*AJ*by*cx*cy - AQ*ax*bz*cx*cy + 
     -       AQ*bx*bz*cx*cy + AD*AJ*ax*cy**2 - AD*AJ*bx*cy**2 - 
     -       AQ*az*by*cy**2 + AQ*by*bz*cy**2 + 
     -       AE*(AJ*AV + AQ*
     -           (-2*ay*bx + ax*by + bx*by + 2*ay*cx - 
     -             2*by*cx - ax*cy + bx*cy)) - AD*AQ*az*cz + 
     -       AD*AJ*az*bx*cz - AQ*ax*bx*by*cz + AQ*bx**2*by*cz - 
     -       AQ*ay*by**2*cz + AQ*by**3*cz + AD*AQ*bz*cz - 
     -       2*AD*AJ*ax*bz*cz + AD*AJ*bx*bz*cz + 
     -       AQ*az*by*bz*cz - 2*AQ*ay*bz**2*cz + 
     -       AQ*by*bz**2*cz - AD*AJ*az*cx*cz + AQ*ax*by*cx*cz - 
     -       AQ*bx*by*cx*cz + AD*AJ*bz*cx*cz + AQ*ay*by*cy*cz - 
     -       AQ*by**2*cy*cz - AQ*az*bz*cy*cz + AQ*bz**2*cy*cz + 
     -       AD*AJ*ax*cz**2 - AD*AJ*bx*cz**2 + AQ*ay*bz*cz**2 - 
     -       AQ*by*bz*cz**2 + 
     -       AF*(AJ*AW + AQ*
     -           (-2*az*bx + ax*bz + bx*bz + 2*az*cx - 
     -             2*bz*cx - ax*cz + bx*cz)) + AQ*az*bx**2*dy + 
     -       AQ*az*by**2*dy - AQ*ax*bx*bz*dy - AQ*ay*by*bz*dy - 
     -       2*AQ*az*bx*cx*dy + AQ*ax*bz*cx*dy + 
     -       AQ*bx*bz*cx*dy + AQ*az*cx**2*dy - AQ*bz*cx**2*dy - 
     -       2*AQ*az*by*cy*dy + AQ*ay*bz*cy*dy + 
     -       AQ*by*bz*cy*dy + AQ*az*cy**2*dy - AQ*bz*cy**2*dy + 
     -       AQ*ax*bx*cz*dy - AQ*bx**2*cz*dy + AQ*ay*by*cz*dy - 
     -       AQ*by**2*cz*dy - AQ*ax*cx*cz*dy + AQ*bx*cx*cz*dy - 
     -       AQ*ay*cy*cz*dy + AQ*by*cy*cz*dy - AQ*ay*bx**2*dz + 
     -       AQ*ax*bx*by*dz + AQ*az*by*bz*dz - AQ*ay*bz**2*dz + 
     -       2*AQ*ay*bx*cx*dz - AQ*ax*by*cx*dz - 
     -       AQ*bx*by*cx*dz - AQ*ay*cx**2*dz + AQ*by*cx**2*dz - 
     -       AQ*ax*bx*cy*dz + AQ*bx**2*cy*dz - AQ*az*bz*cy*dz + 
     -       AQ*bz**2*cy*dz + AQ*ax*cx*cy*dz - AQ*bx*cx*cy*dz - 
     -       AQ*az*by*cz*dz + 2*AQ*ay*bz*cz*dz - 
     -       AQ*by*bz*cz*dz + AQ*az*cy*cz*dz - AQ*bz*cy*cz*dz - 
     -       AQ*ay*cz**2*dz + AQ*by*cz**2*dz))/
     -  (AQ**1.5*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))

      dcy= (AC*AE*AF*AQ*ax*bx + AC*AD*AF*AQ*ay*bx - 
     -    AC*AE*AF*AQ*bx**2 + AC*AE*AF*AK*ay*bx**2 + 
     -    AE*AF*AG*AQ*ay*bx**2 - AD*AE*AI*AQ*ay*bx**2 + 
     -    AC*AF**2*AK*az*bx**2 + AF**2*AG*AQ*az*bx**2 - 
     -    AD*AF*AI*AQ*az*bx**2 - AC*AE*AQ*ay*bx**3 - 
     -    2*AC*AD*AF*AQ*ax*by - 2*AC*AF**2*AQ*az*by + 
     -    AC*AD*AF*AQ*bx*by - AC*AE*AF*AK*ax*bx*by - 
     -    AE*AF*AG*AQ*ax*bx*by + AD*AE*AI*AQ*ax*bx*by - 
     -    AC*AD*AF*AK*ay*bx*by - AD*AF*AG*AQ*ay*bx*by + 
     -    AD**2*AI*AQ*ay*bx*by + AC*AE*AQ*ax*bx**2*by + 
     -    AC*AD*AQ*ay*bx**2*by + AC*AD*AF*AK*ax*by**2 + 
     -    AD*AF*AG*AQ*ax*by**2 - AD**2*AI*AQ*ax*by**2 + 
     -    AC*AF**2*AK*az*by**2 + AF**2*AG*AQ*az*by**2 - 
     -    AD*AF*AI*AQ*az*by**2 - AC*AD*AQ*ax*bx*by**2 + 
     -    AC*AF**2*AQ*ay*bz + AC*AE*AF*AQ*az*bz - 
     -    AC*AF**2*AK*ax*bx*bz - AF**2*AG*AQ*ax*bx*bz + 
     -    AD*AF*AI*AQ*ax*bx*bz - AC*AD*AF*AK*az*bx*bz - 
     -    AD*AF*AG*AQ*az*bx*bz + AD**2*AI*AQ*az*bx*bz + 
     -    AC*AD*AQ*az*bx**2*bz + AC*AF**2*AQ*by*bz - 
     -    AC*AF**2*AK*ay*by*bz - AF**2*AG*AQ*ay*by*bz + 
     -    AD*AF*AI*AQ*ay*by*bz - AC*AE*AF*AK*az*by*bz - 
     -    AE*AF*AG*AQ*az*by*bz + AD*AE*AI*AQ*az*by*bz + 
     -    AC*AF*AQ*ay*bx*by*bz + AC*AE*AQ*az*bx*by*bz - 
     -    AC*AF*AQ*ax*by**2*bz - AC*AE*AF*AQ*bz**2 + 
     -    AC*AD*AF*AK*ax*bz**2 + AD*AF*AG*AQ*ax*bz**2 - 
     -    AD**2*AI*AQ*ax*bz**2 + AC*AE*AF*AK*ay*bz**2 + 
     -    AE*AF*AG*AQ*ay*bz**2 - AD*AE*AI*AQ*ay*bz**2 - 
     -    AC*AD*AQ*ax*bx*bz**2 - AC*AE*AQ*ay*bx*bz**2 + 
     -    AC*AF*AQ*az*bx*bz**2 - AC*AF*AQ*ax*bz**3 - 
     -    AC*AE*AF*AQ*ax*cx - AC*AD*AF*AQ*ay*cx + 
     -    AC*AE*AF*AQ*bx*cx - 2*AC*AE*AF*AK*ay*bx*cx - 
     -    2*AE*AF*AG*AQ*ay*bx*cx + 2*AD*AE*AI*AQ*ay*bx*cx - 
     -    2*AC*AF**2*AK*az*bx*cx - 2*AF**2*AG*AQ*az*bx*cx + 
     -    2*AD*AF*AI*AQ*az*bx*cx + 2*AC*AE*AQ*ay*bx**2*cx + 
     -    AC*AD*AF*AQ*by*cx + AC*AE*AF*AK*ax*by*cx + 
     -    AE*AF*AG*AQ*ax*by*cx - AD*AE*AI*AQ*ax*by*cx + 
     -    AC*AD*AF*AK*ay*by*cx + AD*AF*AG*AQ*ay*by*cx - 
     -    AD**2*AI*AQ*ay*by*cx + AC*AE*AF*AK*bx*by*cx + 
     -    AE*AF*AG*AQ*bx*by*cx - AD*AE*AI*AQ*bx*by*cx - 
     -    AC*AE*AQ*ax*bx*by*cx - AC*AD*AQ*ay*bx*by*cx - 
     -    AC*AE*AQ*bx**2*by*cx - AC*AD*AF*AK*by**2*cx - 
     -    AD*AF*AG*AQ*by**2*cx + AD**2*AI*AQ*by**2*cx + 
     -    AC*AD*AQ*bx*by**2*cx + AC*AF**2*AK*ax*bz*cx + 
     -    AF**2*AG*AQ*ax*bz*cx - AD*AF*AI*AQ*ax*bz*cx + 
     -    AC*AD*AF*AK*az*bz*cx + AD*AF*AG*AQ*az*bz*cx - 
     -    AD**2*AI*AQ*az*bz*cx + AC*AF**2*AK*bx*bz*cx + 
     -    AF**2*AG*AQ*bx*bz*cx - AD*AF*AI*AQ*bx*bz*cx - 
     -    AC*AD*AQ*az*bx*bz*cx - AC*AF*AQ*ay*by*bz*cx + 
     -    AC*AF*AQ*by**2*bz*cx - AC*AD*AF*AK*bz**2*cx - 
     -    AD*AF*AG*AQ*bz**2*cx + AD**2*AI*AQ*bz**2*cx - 
     -    AC*AF*AQ*az*bz**2*cx + AC*AD*AQ*bx*bz**2*cx + 
     -    AC*AF*AQ*bz**3*cx + AC*AE*AF*AK*ay*cx**2 + 
     -    AE*AF*AG*AQ*ay*cx**2 - AD*AE*AI*AQ*ay*cx**2 + 
     -    AC*AF**2*AK*az*cx**2 + AF**2*AG*AQ*az*cx**2 - 
     -    AD*AF*AI*AQ*az*cx**2 - AC*AE*AQ*ay*bx*cx**2 - 
     -    AC*AE*AF*AK*by*cx**2 - AE*AF*AG*AQ*by*cx**2 + 
     -    AD*AE*AI*AQ*by*cx**2 + AC*AE*AQ*bx*by*cx**2 - 
     -    AC*AF**2*AK*bz*cx**2 - AF**2*AG*AQ*bz*cx**2 + 
     -    AD*AF*AI*AQ*bz*cx**2 + 2*AC*AD*AF*AQ*ax*cy + 
     -    2*AC*AF**2*AQ*az*cy - 2*AC*AD*AF*AQ*bx*cy + 
     -    AC*AE*AF*AK*ax*bx*cy + AE*AF*AG*AQ*ax*bx*cy - 
     -    AD*AE*AI*AQ*ax*bx*cy + AC*AD*AF*AK*ay*bx*cy + 
     -    AD*AF*AG*AQ*ay*bx*cy - AD**2*AI*AQ*ay*bx*cy - 
     -    AC*AE*AF*AK*bx**2*cy - AE*AF*AG*AQ*bx**2*cy + 
     -    AD*AE*AI*AQ*bx**2*cy - AC*AE*AQ*ax*bx**2*cy - 
     -    AC*AD*AQ*ay*bx**2*cy + AC*AE*AQ*bx**3*cy - 
     -    2*AC*AD*AF*AK*ax*by*cy - 2*AD*AF*AG*AQ*ax*by*cy + 
     -    2*AD**2*AI*AQ*ax*by*cy - 2*AC*AF**2*AK*az*by*cy - 
     -    2*AF**2*AG*AQ*az*by*cy + 2*AD*AF*AI*AQ*az*by*cy + 
     -    AC*AD*AF*AK*bx*by*cy + AD*AF*AG*AQ*bx*by*cy - 
     -    AD**2*AI*AQ*bx*by*cy + 2*AC*AD*AQ*ax*bx*by*cy - 
     -    AC*AD*AQ*bx**2*by*cy - 2*AC*AF**2*AQ*bz*cy + 
     -    AC*AF**2*AK*ay*bz*cy + AF**2*AG*AQ*ay*bz*cy - 
     -    AD*AF*AI*AQ*ay*bz*cy + AC*AE*AF*AK*az*bz*cy + 
     -    AE*AF*AG*AQ*az*bz*cy - AD*AE*AI*AQ*az*bz*cy - 
     -    AC*AF*AQ*ay*bx*bz*cy - AC*AE*AQ*az*bx*bz*cy + 
     -    AC*AF**2*AK*by*bz*cy + AF**2*AG*AQ*by*bz*cy - 
     -    AD*AF*AI*AQ*by*bz*cy + 2*AC*AF*AQ*ax*by*bz*cy - 
     -    AC*AF*AQ*bx*by*bz*cy - AC*AE*AF*AK*bz**2*cy - 
     -    AE*AF*AG*AQ*bz**2*cy + AD*AE*AI*AQ*bz**2*cy + 
     -    AC*AE*AQ*bx*bz**2*cy - AC*AE*AF*AK*ax*cx*cy - 
     -    AE*AF*AG*AQ*ax*cx*cy + AD*AE*AI*AQ*ax*cx*cy - 
     -    AC*AD*AF*AK*ay*cx*cy - AD*AF*AG*AQ*ay*cx*cy + 
     -    AD**2*AI*AQ*ay*cx*cy + AC*AE*AF*AK*bx*cx*cy + 
     -    AE*AF*AG*AQ*bx*cx*cy - AD*AE*AI*AQ*bx*cx*cy + 
     -    AC*AE*AQ*ax*bx*cx*cy + AC*AD*AQ*ay*bx*cx*cy - 
     -    AC*AE*AQ*bx**2*cx*cy + AC*AD*AF*AK*by*cx*cy + 
     -    AD*AF*AG*AQ*by*cx*cy - AD**2*AI*AQ*by*cx*cy - 
     -    AC*AD*AQ*bx*by*cx*cy + AC*AF*AQ*ay*bz*cx*cy - 
     -    AC*AF*AQ*by*bz*cx*cy + AC*AD*AF*AK*ax*cy**2 + 
     -    AD*AF*AG*AQ*ax*cy**2 - AD**2*AI*AQ*ax*cy**2 + 
     -    AC*AF**2*AK*az*cy**2 + AF**2*AG*AQ*az*cy**2 - 
     -    AD*AF*AI*AQ*az*cy**2 - AC*AD*AF*AK*bx*cy**2 - 
     -    AD*AF*AG*AQ*bx*cy**2 + AD**2*AI*AQ*bx*cy**2 - 
     -    AC*AD*AQ*ax*bx*cy**2 + AC*AD*AQ*bx**2*cy**2 - 
     -    AC*AF**2*AK*bz*cy**2 - AF**2*AG*AQ*bz*cy**2 + 
     -    AD*AF*AI*AQ*bz*cy**2 - AC*AF*AQ*ax*bz*cy**2 + 
     -    AC*AF*AQ*bx*bz*cy**2 - AC*AF**2*AQ*ay*cz - 
     -    AC*AE*AF*AQ*az*cz + AC*AF**2*AK*ax*bx*cz + 
     -    AF**2*AG*AQ*ax*bx*cz - AD*AF*AI*AQ*ax*bx*cz + 
     -    AC*AD*AF*AK*az*bx*cz + AD*AF*AG*AQ*az*bx*cz - 
     -    AD**2*AI*AQ*az*bx*cz - AC*AF**2*AK*bx**2*cz - 
     -    AF**2*AG*AQ*bx**2*cz + AD*AF*AI*AQ*bx**2*cz - 
     -    AC*AD*AQ*az*bx**2*cz + AC*AF**2*AQ*by*cz + 
     -    AC*AF**2*AK*ay*by*cz + AF**2*AG*AQ*ay*by*cz - 
     -    AD*AF*AI*AQ*ay*by*cz + AC*AE*AF*AK*az*by*cz + 
     -    AE*AF*AG*AQ*az*by*cz - AD*AE*AI*AQ*az*by*cz - 
     -    AC*AE*AQ*az*bx*by*cz - AC*AF**2*AK*by**2*cz - 
     -    AF**2*AG*AQ*by**2*cz + AD*AF*AI*AQ*by**2*cz + 
     -    AC*AE*AF*AQ*bz*cz - 2*AC*AD*AF*AK*ax*bz*cz - 
     -    2*AD*AF*AG*AQ*ax*bz*cz + 2*AD**2*AI*AQ*ax*bz*cz - 
     -    2*AC*AE*AF*AK*ay*bz*cz - 2*AE*AF*AG*AQ*ay*bz*cz + 
     -    2*AD*AE*AI*AQ*ay*bz*cz + AC*AD*AF*AK*bx*bz*cz + 
     -    AD*AF*AG*AQ*bx*bz*cz - AD**2*AI*AQ*bx*bz*cz + 
     -    2*AC*AD*AQ*ax*bx*bz*cz + 2*AC*AE*AQ*ay*bx*bz*cz - 
     -    AC*AF*AQ*az*bx*bz*cz - AC*AD*AQ*bx**2*bz*cz + 
     -    AC*AE*AF*AK*by*bz*cz + AE*AF*AG*AQ*by*bz*cz - 
     -    AD*AE*AI*AQ*by*bz*cz - AC*AE*AQ*bx*by*bz*cz + 
     -    2*AC*AF*AQ*ax*bz**2*cz - AC*AF*AQ*bx*bz**2*cz - 
     -    AC*AF**2*AK*ax*cx*cz - AF**2*AG*AQ*ax*cx*cz + 
     -    AD*AF*AI*AQ*ax*cx*cz - AC*AD*AF*AK*az*cx*cz - 
     -    AD*AF*AG*AQ*az*cx*cz + AD**2*AI*AQ*az*cx*cz + 
     -    AC*AF**2*AK*bx*cx*cz + AF**2*AG*AQ*bx*cx*cz - 
     -    AD*AF*AI*AQ*bx*cx*cz + AC*AD*AQ*az*bx*cx*cz + 
     -    AC*AD*AF*AK*bz*cx*cz + AD*AF*AG*AQ*bz*cx*cz - 
     -    AD**2*AI*AQ*bz*cx*cz + AC*AF*AQ*az*bz*cx*cz - 
     -    AC*AD*AQ*bx*bz*cx*cz - AC*AF*AQ*bz**2*cx*cz - 
     -    AC*AF**2*AK*ay*cy*cz - AF**2*AG*AQ*ay*cy*cz + 
     -    AD*AF*AI*AQ*ay*cy*cz - AC*AE*AF*AK*az*cy*cz - 
     -    AE*AF*AG*AQ*az*cy*cz + AD*AE*AI*AQ*az*cy*cz + 
     -    AC*AE*AQ*az*bx*cy*cz + AC*AF**2*AK*by*cy*cz + 
     -    AF**2*AG*AQ*by*cy*cz - AD*AF*AI*AQ*by*cy*cz + 
     -    AC*AE*AF*AK*bz*cy*cz + AE*AF*AG*AQ*bz*cy*cz - 
     -    AD*AE*AI*AQ*bz*cy*cz - AC*AE*AQ*bx*bz*cy*cz + 
     -    AC*AD*AF*AK*ax*cz**2 + AD*AF*AG*AQ*ax*cz**2 - 
     -    AD**2*AI*AQ*ax*cz**2 + AC*AE*AF*AK*ay*cz**2 + 
     -    AE*AF*AG*AQ*ay*cz**2 - AD*AE*AI*AQ*ay*cz**2 - 
     -    AC*AD*AF*AK*bx*cz**2 - AD*AF*AG*AQ*bx*cz**2 + 
     -    AD**2*AI*AQ*bx*cz**2 - AC*AD*AQ*ax*bx*cz**2 - 
     -    AC*AE*AQ*ay*bx*cz**2 + AC*AD*AQ*bx**2*cz**2 - 
     -    AC*AE*AF*AK*by*cz**2 - AE*AF*AG*AQ*by*cz**2 + 
     -    AD*AE*AI*AQ*by*cz**2 + AC*AE*AQ*bx*by*cz**2 - 
     -    AC*AF*AQ*ax*bz*cz**2 + AC*AF*AQ*bx*bz*cz**2 + 
     -    AC*AE*AQ*ay*bx**2*dx - AC*AE*AQ*ax*bx*by*dx - 
     -    AC*AD*AQ*ay*bx*by*dx + AC*AD*AQ*ax*by**2*dx - 
     -    AC*AD*AQ*az*bx*bz*dx - AC*AE*AQ*az*by*bz*dx + 
     -    AC*AD*AQ*ax*bz**2*dx + AC*AE*AQ*ay*bz**2*dx - 
     -    2*AC*AE*AQ*ay*bx*cx*dx + AC*AE*AQ*ax*by*cx*dx + 
     -    AC*AD*AQ*ay*by*cx*dx + AC*AE*AQ*bx*by*cx*dx - 
     -    AC*AD*AQ*by**2*cx*dx + AC*AD*AQ*az*bz*cx*dx - 
     -    AC*AD*AQ*bz**2*cx*dx + AC*AE*AQ*ay*cx**2*dx - 
     -    AC*AE*AQ*by*cx**2*dx + AC*AE*AQ*ax*bx*cy*dx + 
     -    AC*AD*AQ*ay*bx*cy*dx - AC*AE*AQ*bx**2*cy*dx - 
     -    2*AC*AD*AQ*ax*by*cy*dx + AC*AD*AQ*bx*by*cy*dx + 
     -    AC*AE*AQ*az*bz*cy*dx - AC*AE*AQ*bz**2*cy*dx - 
     -    AC*AE*AQ*ax*cx*cy*dx - AC*AD*AQ*ay*cx*cy*dx + 
     -    AC*AE*AQ*bx*cx*cy*dx + AC*AD*AQ*by*cx*cy*dx + 
     -    AC*AD*AQ*ax*cy**2*dx - AC*AD*AQ*bx*cy**2*dx + 
     -    AC*AD*AQ*az*bx*cz*dx + AC*AE*AQ*az*by*cz*dx - 
     -    2*AC*AD*AQ*ax*bz*cz*dx - 2*AC*AE*AQ*ay*bz*cz*dx + 
     -    AC*AD*AQ*bx*bz*cz*dx + AC*AE*AQ*by*bz*cz*dx - 
     -    AC*AD*AQ*az*cx*cz*dx + AC*AD*AQ*bz*cx*cz*dx - 
     -    AC*AE*AQ*az*cy*cz*dx + AC*AE*AQ*bz*cy*cz*dx + 
     -    AC*AD*AQ*ax*cz**2*dx + AC*AE*AQ*ay*cz**2*dx - 
     -    AC*AD*AQ*bx*cz**2*dx - AC*AE*AQ*by*cz**2*dx + 
     -    AA*(AD**2*(AK*BA + 
     -          AQ*(AJ*ay - 2*ax*by + bx*by + by*cx + 2*ax*cy - 
     -             2*bx*cy)) + 
     -       AD*(AE*(AK*AV + 
     -             AQ*(AI*AL + AJ*ax + bx*(-bx + cx))) + 
     -          AF*(-(AJ*AK*AL*ax) + AJ**2*AK*az + AK**3*az + 
     -             AQ*(-2*az*by + ay*bz + by*bz + 2*az*cy - 
     -                2*bz*cy - ay*cz + by*cz) + 
     -             AK*(-(bz*cx**2) + by*bz*cy - bz*cy**2 - 
     -                ay*(by - cy)*(bz - cz) - bx**2*cz - 
     -                by**2*cz + by*cy*cz + bx*cx*(bz + cz))) + 
     -          AQ*AW*(bx - dx)) + AQ*(AE*AV + AF*AW)*(bz - dz))
     -     - AC*AF*AQ*ay*bx*by*dz + AC*AF*AQ*ax*by**2*dz - 
     -    AC*AF*AQ*az*bx*bz*dz + AC*AF*AQ*ax*bz**2*dz + 
     -    AC*AF*AQ*ay*by*cx*dz - AC*AF*AQ*by**2*cx*dz + 
     -    AC*AF*AQ*az*bz*cx*dz - AC*AF*AQ*bz**2*cx*dz + 
     -    AC*AF*AQ*ay*bx*cy*dz - 2*AC*AF*AQ*ax*by*cy*dz + 
     -    AC*AF*AQ*bx*by*cy*dz - AC*AF*AQ*ay*cx*cy*dz + 
     -    AC*AF*AQ*by*cx*cy*dz + AC*AF*AQ*ax*cy**2*dz - 
     -    AC*AF*AQ*bx*cy**2*dz + AC*AF*AQ*az*bx*cz*dz - 
     -    2*AC*AF*AQ*ax*bz*cz*dz + AC*AF*AQ*bx*bz*cz*dz - 
     -    AC*AF*AQ*az*cx*cz*dz + AC*AF*AQ*bz*cx*cz*dz + 
     -    AC*AF*AQ*ax*cz**2*dz - AC*AF*AQ*bx*cz**2*dz + 
     -    AB*AE*(AF*AK*az*bx**2 + AQ*az*bx**3 - 2*AF*AQ*az*by + 
     -       AF*AK*az*by**2 + AQ*az*bx*by**2 + AF*AQ*ay*bz - 
     -       AF*AK*ax*bx*bz - AQ*ax*bx**2*bz + AF*AQ*by*bz - 
     -       AF*AK*ay*by*bz - AQ*ax*by**2*bz + AQ*az*bx*bz**2 - 
     -       AQ*ax*bz**3 - 2*AF*AK*az*bx*cx - 2*AQ*az*bx**2*cx + 
     -       AF*AK*ax*bz*cx + AF*AK*bx*bz*cx + AQ*ax*bx*bz*cx + 
     -       AQ*bx**2*bz*cx - AQ*ay*by*bz*cx + AQ*by**2*bz*cx - 
     -       AQ*az*bz**2*cx + AQ*bz**3*cx + AF*AK*az*cx**2 + 
     -       AQ*az*bx*cx**2 - AF*AK*bz*cx**2 - AQ*bx*bz*cx**2 + 
     -       AE*(AK*AV + AQ*(AI*AL + AJ*ax + bx*(-bx + cx))) + 
     -       2*AF*AQ*az*cy - 2*AF*AK*az*by*cy - 
     -       2*AQ*az*bx*by*cy - 2*AF*AQ*bz*cy + AF*AK*ay*bz*cy + 
     -       AF*AK*by*bz*cy + 2*AQ*ax*by*bz*cy + 
     -       AQ*ay*bz*cx*cy - AQ*by*bz*cx*cy + AF*AK*az*cy**2 + 
     -       AQ*az*bx*cy**2 - AF*AK*bz*cy**2 - AQ*ax*bz*cy**2 + 
     -       AD*(AK*BA + AQ*
     -           (AJ*ay - 2*ax*by + bx*by + by*cx + 2*ax*cy - 
     -             2*bx*cy)) - AF*AQ*ay*cz + AF*AK*ax*bx*cz - 
     -       AF*AK*bx**2*cz + AQ*ax*bx**2*cz - AQ*bx**3*cz + 
     -       AF*AQ*by*cz + AF*AK*ay*by*cz + AQ*ay*bx*by*cz - 
     -       AF*AK*by**2*cz - AQ*bx*by**2*cz - AQ*az*bx*bz*cz + 
     -       2*AQ*ax*bz**2*cz - AQ*bx*bz**2*cz - 
     -       AF*AK*ax*cx*cz + AF*AK*bx*cx*cz - AQ*ax*bx*cx*cz + 
     -       AQ*bx**2*cx*cz + AQ*az*bz*cx*cz - AQ*bz**2*cx*cz - 
     -       AF*AK*ay*cy*cz - AQ*ay*bx*cy*cz + AF*AK*by*cy*cz + 
     -       AQ*bx*by*cy*cz - AQ*ax*bz*cz**2 + AQ*bx*bz*cz**2 - 
     -       AQ*az*bx**2*dx - AQ*az*by**2*dx + AQ*ax*bx*bz*dx + 
     -       AQ*ay*by*bz*dx + 2*AQ*az*bx*cx*dx - 
     -       AQ*ax*bz*cx*dx - AQ*bx*bz*cx*dx - AQ*az*cx**2*dx + 
     -       AQ*bz*cx**2*dx + 2*AQ*az*by*cy*dx - 
     -       AQ*ay*bz*cy*dx - AQ*by*bz*cy*dx - AQ*az*cy**2*dx + 
     -       AQ*bz*cy**2*dx - AQ*ax*bx*cz*dx + AQ*bx**2*cz*dx - 
     -       AQ*ay*by*cz*dx + AQ*by**2*cz*dx + AQ*ax*cx*cz*dx - 
     -       AQ*bx*cx*cz*dx + AQ*ay*cy*cz*dx - AQ*by*cy*cz*dx - 
     -       AQ*ay*bx*by*dz + AQ*ax*by**2*dz - AQ*az*bx*bz*dz + 
     -       AQ*ax*bz**2*dz + AQ*ay*by*cx*dz - AQ*by**2*cx*dz + 
     -       AQ*az*bz*cx*dz - AQ*bz**2*cx*dz + AQ*ay*bx*cy*dz - 
     -       2*AQ*ax*by*cy*dz + AQ*bx*by*cy*dz - 
     -       AQ*ay*cx*cy*dz + AQ*by*cx*cy*dz + AQ*ax*cy**2*dz - 
     -       AQ*bx*cy**2*dz + AQ*az*bx*cz*dz - 
     -       2*AQ*ax*bz*cz*dz + AQ*bx*bz*cz*dz - 
     -       AQ*az*cx*cz*dz + AQ*bz*cx*cz*dz + AQ*ax*cz**2*dz - 
     -       AQ*bx*cz**2*dz))/
     -  (AQ**1.5*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))

      dcz=        (AC*AF**2*AQ*ax*bx + AC*AD*AF*AQ*az*bx - 
     -    AC*AF**2*AQ*bx**2 + AC*AE*AF*AL*ay*bx**2 - 
     -    AE**2*AG*AQ*ay*bx**2 + AD*AE*AH*AQ*ay*bx**2 + 
     -    AC*AF**2*AL*az*bx**2 - AE*AF*AG*AQ*az*bx**2 + 
     -    AD*AF*AH*AQ*az*bx**2 - AC*AF*AQ*ay*BE*bx**2 + 
     -    AC*AF**2*AQ*ay*by + AC*AE*AF*AQ*az*by - 
     -    AC*AE*AF*AL*ax*bx*by + AE**2*AG*AQ*ax*bx*by - 
     -    AD*AE*AH*AQ*ax*bx*by - AC*AD*AF*AL*ay*bx*by + 
     -    AD*AE*AG*AQ*ay*bx*by - AD**2*AH*AQ*ay*bx*by + 
     -    AC*AF*AQ*ax*BE*bx*by - AC*AF*AQ*ay*BF*bx*by - 
     -    AC*AF**2*AQ*by**2 + AC*AD*AF*AL*ax*by**2 - 
     -    AD*AE*AG*AQ*ax*by**2 + AD**2*AH*AQ*ax*by**2 + 
     -    AC*AF**2*AL*az*by**2 - AE*AF*AG*AQ*az*by**2 + 
     -    AD*AF*AH*AQ*az*by**2 + AC*AF*AQ*ax*BF*by**2 - 
     -    2*AC*AD*AF*AQ*ax*bz - 2*AC*AE*AF*AQ*ay*bz + 
     -    AC*AD*AF*AQ*bx*bz - AC*AF**2*AL*ax*bx*bz + 
     -    AE*AF*AG*AQ*ax*bx*bz - AD*AF*AH*AQ*ax*bx*bz - 
     -    AC*AD*AF*AL*az*bx*bz + AD*AE*AG*AQ*az*bx*bz - 
     -    AD**2*AH*AQ*az*bx*bz - AC*AF*AQ*az*BF*bx*bz + 
     -    AC*AE*AF*AQ*by*bz - AC*AF**2*AL*ay*by*bz + 
     -    AE*AF*AG*AQ*ay*by*bz - AD*AF*AH*AQ*ay*by*bz - 
     -    AC*AE*AF*AL*az*by*bz + AE**2*AG*AQ*az*by*bz - 
     -    AD*AE*AH*AQ*az*by*bz + AC*AF*AQ*az*BE*by*bz + 
     -    AC*AD*AF*AL*ax*bz**2 - AD*AE*AG*AQ*ax*bz**2 + 
     -    AD**2*AH*AQ*ax*bz**2 + AC*AE*AF*AL*ay*bz**2 - 
     -    AE**2*AG*AQ*ay*bz**2 + AD*AE*AH*AQ*ay*bz**2 - 
     -    AC*AF*AQ*ay*BE*bz**2 + AC*AF*AQ*ax*BF*bz**2 - 
     -    AC*AF**2*AQ*ax*cx - AC*AD*AF*AQ*az*cx + 
     -    AC*AF**2*AQ*bx*cx - 2*AC*AE*AF*AL*ay*bx*cx + 
     -    2*AE**2*AG*AQ*ay*bx*cx - 2*AD*AE*AH*AQ*ay*bx*cx - 
     -    2*AC*AF**2*AL*az*bx*cx + 2*AE*AF*AG*AQ*az*bx*cx - 
     -    2*AD*AF*AH*AQ*az*bx*cx + 2*AC*AF*AQ*ay*BE*bx*cx + 
     -    AC*AE*AF*AL*ax*by*cx - AE**2*AG*AQ*ax*by*cx + 
     -    AD*AE*AH*AQ*ax*by*cx + AC*AD*AF*AL*ay*by*cx - 
     -    AD*AE*AG*AQ*ay*by*cx + AD**2*AH*AQ*ay*by*cx - 
     -    AC*AF*AQ*ax*BE*by*cx + AC*AF*AQ*ay*BF*by*cx + 
     -    AC*AE*AF*AL*bx*by*cx - AE**2*AG*AQ*bx*by*cx + 
     -    AD*AE*AH*AQ*bx*by*cx - AC*AF*AQ*BE*bx*by*cx - 
     -    AC*AD*AF*AL*by**2*cx + AD*AE*AG*AQ*by**2*cx - 
     -    AD**2*AH*AQ*by**2*cx - AC*AF*AQ*BF*by**2*cx + 
     -    AC*AD*AF*AQ*bz*cx + AC*AF**2*AL*ax*bz*cx - 
     -    AE*AF*AG*AQ*ax*bz*cx + AD*AF*AH*AQ*ax*bz*cx + 
     -    AC*AD*AF*AL*az*bz*cx - AD*AE*AG*AQ*az*bz*cx + 
     -    AD**2*AH*AQ*az*bz*cx + AC*AF*AQ*az*BF*bz*cx + 
     -    AC*AF**2*AL*bx*bz*cx - AE*AF*AG*AQ*bx*bz*cx + 
     -    AD*AF*AH*AQ*bx*bz*cx - AC*AD*AF*AL*bz**2*cx + 
     -    AD*AE*AG*AQ*bz**2*cx - AD**2*AH*AQ*bz**2*cx - 
     -    AC*AF*AQ*BF*bz**2*cx + AC*AE*AF*AL*ay*cx**2 - 
     -    AE**2*AG*AQ*ay*cx**2 + AD*AE*AH*AQ*ay*cx**2 + 
     -    AC*AF**2*AL*az*cx**2 - AE*AF*AG*AQ*az*cx**2 + 
     -    AD*AF*AH*AQ*az*cx**2 - AC*AF*AQ*ay*BE*cx**2 - 
     -    AC*AE*AF*AL*by*cx**2 + AE**2*AG*AQ*by*cx**2 - 
     -    AD*AE*AH*AQ*by*cx**2 + AC*AF*AQ*BE*by*cx**2 - 
     -    AC*AF**2*AL*bz*cx**2 + AE*AF*AG*AQ*bz*cx**2 - 
     -    AD*AF*AH*AQ*bz*cx**2 - AC*AF**2*AQ*ay*cy - 
     -    AC*AE*AF*AQ*az*cy + AC*AE*AF*AL*ax*bx*cy - 
     -    AE**2*AG*AQ*ax*bx*cy + AD*AE*AH*AQ*ax*bx*cy + 
     -    AC*AD*AF*AL*ay*bx*cy - AD*AE*AG*AQ*ay*bx*cy + 
     -    AD**2*AH*AQ*ay*bx*cy - AC*AF*AQ*ax*BE*bx*cy + 
     -    AC*AF*AQ*ay*BF*bx*cy - AC*AE*AF*AL*bx**2*cy + 
     -    AE**2*AG*AQ*bx**2*cy - AD*AE*AH*AQ*bx**2*cy + 
     -    AC*AF*AQ*BE*bx**2*cy + AC*AF**2*AQ*by*cy - 
     -    2*AC*AD*AF*AL*ax*by*cy + 2*AD*AE*AG*AQ*ax*by*cy - 
     -    2*AD**2*AH*AQ*ax*by*cy - 2*AC*AF**2*AL*az*by*cy + 
     -    2*AE*AF*AG*AQ*az*by*cy - 2*AD*AF*AH*AQ*az*by*cy - 
     -    2*AC*AF*AQ*ax*BF*by*cy + AC*AD*AF*AL*bx*by*cy - 
     -    AD*AE*AG*AQ*bx*by*cy + AD**2*AH*AQ*bx*by*cy + 
     -    AC*AF*AQ*BF*bx*by*cy + AC*AE*AF*AQ*bz*cy + 
     -    AC*AF**2*AL*ay*bz*cy - AE*AF*AG*AQ*ay*bz*cy + 
     -    AD*AF*AH*AQ*ay*bz*cy + AC*AE*AF*AL*az*bz*cy - 
     -    AE**2*AG*AQ*az*bz*cy + AD*AE*AH*AQ*az*bz*cy - 
     -    AC*AF*AQ*az*BE*bz*cy + AC*AF**2*AL*by*bz*cy - 
     -    AE*AF*AG*AQ*by*bz*cy + AD*AF*AH*AQ*by*bz*cy - 
     -    AC*AE*AF*AL*bz**2*cy + AE**2*AG*AQ*bz**2*cy - 
     -    AD*AE*AH*AQ*bz**2*cy + AC*AF*AQ*BE*bz**2*cy - 
     -    AC*AE*AF*AL*ax*cx*cy + AE**2*AG*AQ*ax*cx*cy - 
     -    AD*AE*AH*AQ*ax*cx*cy - AC*AD*AF*AL*ay*cx*cy + 
     -    AD*AE*AG*AQ*ay*cx*cy - AD**2*AH*AQ*ay*cx*cy + 
     -    AC*AF*AQ*ax*BE*cx*cy - AC*AF*AQ*ay*BF*cx*cy + 
     -    AC*AE*AF*AL*bx*cx*cy - AE**2*AG*AQ*bx*cx*cy + 
     -    AD*AE*AH*AQ*bx*cx*cy - AC*AF*AQ*BE*bx*cx*cy + 
     -    AC*AD*AF*AL*by*cx*cy - AD*AE*AG*AQ*by*cx*cy + 
     -    AD**2*AH*AQ*by*cx*cy + AC*AF*AQ*BF*by*cx*cy + 
     -    AC*AD*AF*AL*ax*cy**2 - AD*AE*AG*AQ*ax*cy**2 + 
     -    AD**2*AH*AQ*ax*cy**2 + AC*AF**2*AL*az*cy**2 - 
     -    AE*AF*AG*AQ*az*cy**2 + AD*AF*AH*AQ*az*cy**2 + 
     -    AC*AF*AQ*ax*BF*cy**2 - AC*AD*AF*AL*bx*cy**2 + 
     -    AD*AE*AG*AQ*bx*cy**2 - AD**2*AH*AQ*bx*cy**2 - 
     -    AC*AF*AQ*BF*bx*cy**2 - AC*AF**2*AL*bz*cy**2 + 
     -    AE*AF*AG*AQ*bz*cy**2 - AD*AF*AH*AQ*bz*cy**2 + 
     -    2*AC*AD*AF*AQ*ax*cz + 2*AC*AE*AF*AQ*ay*cz - 
     -    2*AC*AD*AF*AQ*bx*cz + AC*AF**2*AL*ax*bx*cz - 
     -    AE*AF*AG*AQ*ax*bx*cz + AD*AF*AH*AQ*ax*bx*cz + 
     -    AC*AD*AF*AL*az*bx*cz - AD*AE*AG*AQ*az*bx*cz + 
     -    AD**2*AH*AQ*az*bx*cz + AC*AF*AQ*az*BF*bx*cz - 
     -    AC*AF**2*AL*bx**2*cz + AE*AF*AG*AQ*bx**2*cz - 
     -    AD*AF*AH*AQ*bx**2*cz - 2*AC*AE*AF*AQ*by*cz + 
     -    AC*AF**2*AL*ay*by*cz - AE*AF*AG*AQ*ay*by*cz + 
     -    AD*AF*AH*AQ*ay*by*cz + AC*AE*AF*AL*az*by*cz - 
     -    AE**2*AG*AQ*az*by*cz + AD*AE*AH*AQ*az*by*cz - 
     -    AC*AF*AQ*az*BE*by*cz - AC*AF**2*AL*by**2*cz + 
     -    AE*AF*AG*AQ*by**2*cz - AD*AF*AH*AQ*by**2*cz - 
     -    2*AC*AD*AF*AL*ax*bz*cz + 2*AD*AE*AG*AQ*ax*bz*cz - 
     -    2*AD**2*AH*AQ*ax*bz*cz - 2*AC*AE*AF*AL*ay*bz*cz + 
     -    2*AE**2*AG*AQ*ay*bz*cz - 2*AD*AE*AH*AQ*ay*bz*cz + 
     -    2*AC*AF*AQ*ay*BE*bz*cz - 2*AC*AF*AQ*ax*BF*bz*cz + 
     -    AC*AD*AF*AL*bx*bz*cz - AD*AE*AG*AQ*bx*bz*cz + 
     -    AD**2*AH*AQ*bx*bz*cz + AC*AF*AQ*BF*bx*bz*cz + 
     -    AC*AE*AF*AL*by*bz*cz - AE**2*AG*AQ*by*bz*cz + 
     -    AD*AE*AH*AQ*by*bz*cz - AC*AF*AQ*BE*by*bz*cz - 
     -    AC*AF**2*AL*ax*cx*cz + AE*AF*AG*AQ*ax*cx*cz - 
     -    AD*AF*AH*AQ*ax*cx*cz - AC*AD*AF*AL*az*cx*cz + 
     -    AD*AE*AG*AQ*az*cx*cz - AD**2*AH*AQ*az*cx*cz - 
     -    AC*AF*AQ*az*BF*cx*cz + AC*AF**2*AL*bx*cx*cz - 
     -    AE*AF*AG*AQ*bx*cx*cz + AD*AF*AH*AQ*bx*cx*cz + 
     -    AC*AD*AF*AL*bz*cx*cz - AD*AE*AG*AQ*bz*cx*cz + 
     -    AD**2*AH*AQ*bz*cx*cz + AC*AF*AQ*BF*bz*cx*cz - 
     -    AC*AF**2*AL*ay*cy*cz + AE*AF*AG*AQ*ay*cy*cz - 
     -    AD*AF*AH*AQ*ay*cy*cz - AC*AE*AF*AL*az*cy*cz + 
     -    AE**2*AG*AQ*az*cy*cz - AD*AE*AH*AQ*az*cy*cz + 
     -    AC*AF*AQ*az*BE*cy*cz + AC*AF**2*AL*by*cy*cz - 
     -    AE*AF*AG*AQ*by*cy*cz + AD*AF*AH*AQ*by*cy*cz + 
     -    AC*AE*AF*AL*bz*cy*cz - AE**2*AG*AQ*bz*cy*cz + 
     -    AD*AE*AH*AQ*bz*cy*cz - AC*AF*AQ*BE*bz*cy*cz + 
     -    AC*AD*AF*AL*ax*cz**2 - AD*AE*AG*AQ*ax*cz**2 + 
     -    AD**2*AH*AQ*ax*cz**2 + AC*AE*AF*AL*ay*cz**2 - 
     -    AE**2*AG*AQ*ay*cz**2 + AD*AE*AH*AQ*ay*cz**2 - 
     -    AC*AF*AQ*ay*BE*cz**2 + AC*AF*AQ*ax*BF*cz**2 - 
     -    AC*AD*AF*AL*bx*cz**2 + AD*AE*AG*AQ*bx*cz**2 - 
     -    AD**2*AH*AQ*bx*cz**2 - AC*AF*AQ*BF*bx*cz**2 - 
     -    AC*AE*AF*AL*by*cz**2 + AE**2*AG*AQ*by*cz**2 - 
     -    AD*AE*AH*AQ*by*cz**2 + AC*AF*AQ*BE*by*cz**2 + 
     -    AB*(AQ*(AF*AW + AD*BA)*BE + 
     -       AE**2*(-(AJ*AK*AL*ax) + AJ**2*AL*ay + AL**3*ay + 
     -          AQ*(AK*az - 2*ay*bz + by*bz + bz*cy + 
     -             2*ay*cz - 2*by*cz) + 
     -          AL*(-(by*cx**2) - bx**2*cy - bz**2*cy + 
     -             bx*cx*(by + cy) - az*(by - cy)*(bz - cz) + 
     -             by*bz*cz + bz*cy*cz - by*cz**2)) + 
     -       AE*(AQ*BA*BF + 
     -          AF*(AL*AW + 
     -             AQ*(AH*AK + AJ*ax + bx*(-bx + cx))) + 
     -          AD*(AL*BA + 
     -             AQ*(AJ*az - 2*ax*bz + bx*bz + bz*cx + 
     -                2*ax*cz - 2*bx*cz)))) + 
     -    AA*(-(AQ*(AE*AV + AF*AW)*BF) + 
     -       AD**2*(AL*BA + 
     -          AQ*(AJ*az - 2*ax*bz + bx*bz + bz*cx + 
     -             2*ax*cz - 2*bx*cz)) + 
     -       AD*(-(AQ*AV*BE) + 
     -          AF*(AL*AW + 
     -             AQ*(AH*AK + AJ*ax + bx*(-bx + cx))) + 
     -          AE*(-(AJ*AK*AL*ax) + AJ**2*AL*ay + AL**3*ay + 
     -             AQ*(AK*az - 2*ay*bz + by*bz + bz*cy + 
     -                2*ay*cz - 2*by*cz) + 
     -             AL*(-(by*cx**2) - bx**2*cy - bz**2*cy + 
     -                bx*cx*(by + cy) - 
     -                az*(by - cy)*(bz - cz) + by*bz*cz + 
     -                bz*cy*cz - by*cz**2)))))/
     -  (AQ**1.5*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))
    

      ddx=   (AB*((AE*AK + AF*AL)*BM + AD*AL*BN) - 
     -    AC*(AD*AK*BN + (AE*AK + AF*AL)*BO) + 
     -    AA*AD*(-(AL**3*ay) + AJ**2*(-(AL*ay) + AK*az) + 
     -       AK*(AK**2*az + bx*bz*cx - bz*cx**2 + by*bz*cy - 
     -          bz*cy**2 - ay*(by - cy)*(bz - cz) - bx**2*cz - 
     -          by**2*cz + bx*cx*cz + by*cy*cz) + 
     -       AL*(by*cx**2 + bx**2*cy + bz**2*cy - 
     -          bx*cx*(by + cy) + az*(by - cy)*(bz - cz) - 
     -          by*bz*cz - bz*cy*cz + by*cz**2)))/
     -  (Sqrt(AQ)*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))

      ddy=        (-(AA*(AE*AL*AV + (AD*AJ + AF*AL)*AW)) - 
     -    AB*AE*(AJ*AW - AL*BA) + 
     -    AC*(AE*AJ*AV + (AD*AJ + AF*AL)*BA))/
     -  (Sqrt(AQ)*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))

      ddz=        (AA*(AD*AJ*AV + AK*(AE*AV + AF*AW)) + 
     -    AC*AF*(AJ*AV - AK*BA) - 
     -    AB*(AF*AJ*AW + (AD*AJ + AE*AK)*BA))/
     -  (Sqrt(AQ)*((AA*AD + AB*AE + AC*AF)**2 + 
     -      (AE*AV + AF*AW + AD*BA)**2/AQ))
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
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        increment=(angle_abcd-zero_value)**2
c        write(*,*)i,angle_abcd*57.29,zero_value*57.29,"dihedral angle"
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
          write(*,*)"entering dextwu"
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
c      real*8 J1x,J1y,J1z,J2x,J2y,J2z,J3x,J3y,J3z,J4x,J4y,J4z
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
      write(*,*)"stretches done"
      
C     Loop over 6-rings
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
      write(*,*)"angles done"

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
        ax=p(J1*3-2)
        ay=p(J1*3-1)
        az=p(J1*3)
        bx=p(J2*3-2)
        by=p(J2*3-1)
        bz=p(J2*3)
        cx=p(J3*3-2)
        cy=p(J3*3-1)
        cz=p(J3*3)
        dx=p(J4*3-2)
        dy=p(J4*3-1)
        dz=p(J4*3)
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2  dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz)
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
c        if(fac_AC .ge. dpi) fac_AC=fac_AC-2*dpi
c        fac_ac=dabs(fac_ac)
c        write(*,*)"ac",fac_ac
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
        if(angle_abcd.gt.dpi)angle_abcd=angle_abcd-2*dpi
        if(angle_abcd.lt.-dpi)angle_abcd=angle_abcd+2*dpi
        angle_abcd=dabs(angle_abcd)
        dE_over_dc=2*force_constant*(angle_abcd-zero_value)
c derivations of the energy with respect the x,y,z of each of the four atoms
        x(J1*3-2)=x(J1*3-2)+dax*dE_over_dc
        x(J1*3-1)=x(J1*3-1)+day*dE_over_dc
        x(J1*3)=x(J1*3)+daz*dE_over_dc
        x(J2*3-2)=x(J2*3-2)+dbx*dE_over_dc
        x(J2*3-1)=x(J2*3-1)+dby*dE_over_dc
        x(J2*3)=x(J2*3)+dbz*dE_over_dc
        x(J3*3-2)=x(J3*3-2)+dcx*dE_over_dc
        x(J3*3-1)=x(J3*3-1)+dcy*dE_over_dc
        x(J3*3)=x(J3*3)+dcz*dE_over_dc
        x(J4*3-2)=x(J4*3-2)+ddx*dE_over_dc
        x(J4*3-1)=x(J4*3-1)+ddy*dE_over_dc
        x(J4*3)=x(J4*3)+ddz*dE_over_dc
      enddo
      write(*,*)"d,0: ",angle_abcd,zero_value," (should be similar)"
      write(*,*)"dihedrals done"
      return
      END
