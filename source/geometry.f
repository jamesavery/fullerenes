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
c vectors ab, bc and cd
      ab_x=ax-bx
      ab_y=ay-by
      ab_z=az-bz
      bc_x=bx-cx
      bc_y=by-cy
      bc_z=bz-cz
      cd_x=cx-dx
      cd_y=cy-dy
      cd_z=cz-dz
c vector bc normed to length 1
      bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      bc_1_x=bc_x*bc_length_inv
      bc_1_y=bc_y*bc_length_inv
      bc_1_z=bc_z*bc_length_inv
c normal vectors on abc and bcd
c and the signs are this way because one of the two vectors points in the wrong direction
      abc_x=-ab_y*bc_z + ab_z*bc_y
      abc_y=-ab_z*bc_x + ab_x*bc_z
      abc_z=-ab_x*bc_y + ab_y*bc_x
      bcd_x=-bc_y*cd_z + bc_z*cd_y
      bcd_y=-bc_z*cd_x + bc_x*cd_z
      bcd_z=-bc_x*cd_y + bc_y*cd_x
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
c abc \times bcd
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
     2  dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3  dihedral_abcd)
      implicit real*8 (a-z)
C at first the dihedral (copied from above)
c vectors ab, bc and cd
      ab_x=ax-bx
      ab_y=ay-by
      ab_z=az-bz
      bc_x=bx-cx
      bc_y=by-cy
      bc_z=bz-cz
      cd_x=cx-dx
      cd_y=cy-dy
      cd_z=cz-dz
c vector bc normed to length 1
      bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      bc1_x=bc_x*bc_length_inv
      bc1_y=bc_y*bc_length_inv
      bc1_z=bc_z*bc_length_inv
c normal vectors on abc and bcd
c and the signs are this way because one of the two vectors points in the wrong direction
      abc_x=-ab_y*bc_z + ab_z*bc_y
      abc_y=-ab_z*bc_x + ab_x*bc_z
      abc_z=-ab_x*bc_y + ab_y*bc_x
      bcd_x=-bc_y*cd_z + bc_z*cd_y
      bcd_y=-bc_z*cd_x + bc_x*cd_z
      bcd_z=-bc_x*cd_y + bc_y*cd_x
c their respective lengths
      abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
c normal vectors (length 1) on abc and bcd
      abc1_x=abc_x*abc_length_inv
      abc1_y=abc_y*abc_length_inv
      abc1_z=abc_z*abc_length_inv
      bcd1_x=bcd_x*bcd_length_inv
      bcd1_y=bcd_y*bcd_length_inv
      bcd1_z=bcd_z*bcd_length_inv
c abc \times bcd
      aux_x=abc1_y*bc1_z-bc1_y*abc1_z
      aux_y=abc1_z*bc1_x-bc1_z*abc1_x
      aux_z=abc1_x*bc1_y-bc1_x*abc1_y
c two auxiliary reals
c     x=\vec abc1 \cdot \vec bcd_1
      x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
c     y=\vec aux  \cdot \vec bcd_1
      y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
c the result
      dihedral_abcd=atan2(y, x)

C THE DERIVATIVES
c to be read from bottom to top

c derivatives of single vectors
c all other combinations are zero
      dab_x__dax=1
      dab_x__dbx=-1
      dab_y__day=1
      dab_y__dby=-1
      dab_z__daz=1
      dab_z__dbz=-1

      dbc_x__dbx=1
      dbc_x__dcx=-1
      dbc_y__dby=1
      dbc_y__dcy=-1
      dbc_z__dbz=1
      dbc_z__dcz=-1

      dcd_x__dcx=1
      dcd_x__ddx=-1
      dcd_y__dcy=1
      dcd_y__ddy=-1
      dcd_z__dcz=1
      dcd_z__ddz=-1

c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      bc_length_inv_cub=bc_length_inv**3
      dbc_length_inv__dbc_x=-bc_x*bc_length_inv_cub
      dbc_length_inv__dbc_y=-bc_y*bc_length_inv_cub
      dbc_length_inv__dbc_z=-bc_z*bc_length_inv_cub

c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
c and the other terms are zero
      dbc_length_inv__dbx = dbc_length_inv__dbc_x*dbc_x__dbx
      dbc_length_inv__dby = dbc_length_inv__dbc_y*dbc_y__dby
      dbc_length_inv__dbz = dbc_length_inv__dbc_z*dbc_z__dbz
      dbc_length_inv__dcx = dbc_length_inv__dbc_x*dbc_x__dcx
      dbc_length_inv__dcy = dbc_length_inv__dbc_y*dbc_y__dcy
      dbc_length_inv__dcz = dbc_length_inv__dbc_z*dbc_z__dcz

c bc1_y=bc_y*bc_length_inv
c bc1_z=bc_z*bc_length_inv
      dbc1_x__dbc_x=bc_length_inv
      dbc1_y__dbc_y=bc_length_inv
      dbc1_z__dbc_z=bc_length_inv

      dbc1_x__dbc_length_inv=bc_x
      dbc1_y__dbc_length_inv=bc_y
      dbc1_z__dbc_length_inv=bc_z

c bc1_x=bc_x*bc_length_inv
      dbc1_x__dbx=
     2 dbc1_x__dbc_x*dbc_x__dbx +
     2 dbc1_x__dbc_length_inv*dbc_length_inv__dbx
      dbc1_x__dby=
     2 dbc1_x__dbc_length_inv*dbc_length_inv__dby
      dbc1_x__dbz=
     2 dbc1_x__dbc_length_inv*dbc_length_inv__dbz
      dbc1_x__dcx=
     2 dbc1_x__dbc_x*dbc_x__dcx +
     2 dbc1_x__dbc_length_inv*dbc_length_inv__dcx
      dbc1_x__dcy=
     2 dbc1_x__dbc_length_inv*dbc_length_inv__dcy
      dbc1_x__dcz=
     2 dbc1_x__dbc_length_inv*dbc_length_inv__dcz
c bc1_y=bc_y*bc_length_inv
      dbc1_y__dbx=
     2 dbc1_y__dbc_length_inv*dbc_length_inv__dbx
      dbc1_y__dby=
     2 dbc1_y__dbc_y*dbc_y__dby +
     2 dbc1_y__dbc_length_inv*dbc_length_inv__dby
      dbc1_y__dbz=
     2 dbc1_y__dbc_length_inv*dbc_length_inv__dbz
      dbc1_y__dcx=
     2 dbc1_y__dbc_length_inv*dbc_length_inv__dcx
      dbc1_y__dcy=
     2 dbc1_y__dbc_y*dbc_y__dcy +
     2 dbc1_y__dbc_length_inv*dbc_length_inv__dcy
      dbc1_y__dcz=
     2 dbc1_y__dbc_length_inv*dbc_length_inv__dcz
c bc1_z=bc_z*bc_length_inv
      dbc1_z__dbx=
     2 dbc1_z__dbc_length_inv*dbc_length_inv__dbx
      dbc1_z__dby=
     2 dbc1_z__dbc_length_inv*dbc_length_inv__dby
      dbc1_z__dbz=
     2 dbc1_z__dbc_z*dbc_z__dbz +
     2 dbc1_z__dbc_length_inv*dbc_length_inv__dbz
      dbc1_z__dcx=
     2 dbc1_z__dbc_length_inv*dbc_length_inv__dcx
      dbc1_z__dcy=
     2 dbc1_z__dbc_length_inv*dbc_length_inv__dcy
      dbc1_z__dcz=
     2 dbc1_z__dbc_z*dbc_z__dcz +
     2 dbc1_z__dbc_length_inv*dbc_length_inv__dcz

c abc_x=-ab_y*bc_z + ab_z*bc_y
      dabc_x__dab_y=-bc_z
      dabc_x__dbc_z=-ab_y
      dabc_x__dab_z=bc_y
      dabc_x__dbc_y=ab_z
c abc_y=-ab_z*bc_x + ab_x*bc_z
      dabc_y__dab_z=-bc_x
      dabc_y__dbc_x=-ab_z
      dabc_y__dab_x=bc_z
      dabc_y__dbc_z=ab_x
c abc_z=-ab_x*bc_y + ab_y*bc_x
      dabc_z__dab_x=-bc_y
      dabc_z__dbc_y=-ab_x
      dabc_z__dab_y=bc_x
      dabc_z__dbc_x=ab_y
c bcd_x=-bc_y*cd_z + bc_z*cd_y
      dbcd_x__dbc_y=-cd_z
      dbcd_x__dcd_z=-bc_y
      dbcd_x__dbc_z=cd_y
      dbcd_x__dcd_y=bc_z
c bcd_y=-bc_z*cd_x + bc_x*cd_z
      dbcd_y__dbc_z=-cd_x
      dbcd_y__dcd_x=-bc_z
      dbcd_y__dbc_x=cd_z
      dbcd_y__dcd_z=bc_x
c bcd_z=-bc_x*cd_y + bc_y*cd_x
      dbcd_z__dbc_x=-cd_y
      dbcd_z__dcd_y=-bc_x
      dbcd_z__dbc_y=cd_x
      dbcd_z__dcd_x=bc_y

c abc_x=-ab_y*bc_z + ab_z*bc_y
c      dabc_x__dax=0
      dabc_x__day=dabc_x__dab_y*dab_y__day
      dabc_x__daz=dabc_x__dab_z*dab_z__daz
c      dabc_x__dbx=0
      dabc_x__dby=dabc_x__dab_y*dab_y__dby + dabc_x__dbc_y*dbc_y__dby
      dabc_x__dbz=dabc_x__dbc_z*dbc_z__dbz + dabc_x__dab_z*dab_z__dbz
c      dabc_x__dcx=0
      dabc_x__dcy=dabc_x__dbc_y*dbc_y__dcy
      dabc_x__dcz=dabc_x__dbc_z*dbc_z__dcz
c abc_y=-ab_z*bc_x + ab_x*bc_z
      dabc_y__dax=dabc_y__dab_x*dab_x__dax 
c      dabc_y__day=0
      dabc_y__daz=dabc_y__dab_z*dab_z__daz 
      dabc_y__dbx=dabc_y__dbc_x*dbc_x__dbx + dabc_y__dab_x*dab_x__dbx
c      dabc_y__dby=0
      dabc_y__dbz=dabc_y__dab_z*dab_z__dbz + dabc_y__dbc_z*dbc_z__dbz
      dabc_y__dcx=dabc_y__dbc_x*dbc_x__dcx
c      dabc_y__dcy=0
      dabc_y__dcz=dabc_y__dbc_z*dbc_z__dcz
c abc_z=-ab_x*bc_y + ab_y*bc_x
      dabc_z__dax=dabc_z__dab_x*dab_x__dax
      dabc_z__day=dabc_z__dab_y*dab_y__day 
c      dabc_z__daz=0
      dabc_z__dbx=dabc_z__dbc_x*dbc_x__dbx + dabc_z__dab_x*dab_x__dbx
      dabc_z__dby=dabc_z__dbc_y*dbc_y__dby + dabc_z__dab_y*dab_y__dby 
c      dabc_z__dbz=0 
      dabc_z__dcx=dabc_z__dbc_x*dbc_x__dcx
      dabc_z__dcy=dabc_z__dbc_y*dbc_y__dcy
c      dabc_z__dcz=0
c bcd_x=-bc_y*cd_z + bc_z*cd_y
c      dbcd_x__dbx=0
      dbcd_x__dby=dbcd_x__dbc_y*dbc_y__dby 
      dbcd_x__dbz=dbcd_x__dbc_z*dbc_z__dbz 
c      dbcd_x__dcx=0
      dbcd_x__dcy=dbcd_x__dbc_y*dbc_y__dcy + dbcd_x__dcd_y*dcd_y__dcy
      dbcd_x__dcz=dbcd_x__dcd_z*dcd_z__dcz + dbcd_x__dbc_z*dbc_z__dcz 
c      dbcd_x__ddx=0
      dbcd_x__ddy=dbcd_x__dcd_y*dcd_y__ddy
      dbcd_x__ddz=dbcd_x__dcd_z*dcd_z__ddz
c bcd_y=-bc_z*cd_x + bc_x*cd_z
      dbcd_y__dbx=dbcd_y__dbc_x*dbc_x__dbx
c      dbcd_y__dby=0
      dbcd_y__dbz=dbcd_y__dbc_z*dbc_z__dbz
      dbcd_y__dcx=dbcd_y__dcd_x*dcd_x__dcx + dbcd_y__dbc_x*dbc_x__dcx 
c      dbcd_y__dcy=0
      dbcd_y__dcz=dbcd_y__dbc_z*dbc_z__dcz + dbcd_y__dcd_z*dcd_z__dcz
      dbcd_y__ddx=dbcd_y__dcd_x*dcd_x__ddx
c      dbcd_y__ddy=0
      dbcd_y__ddz=dbcd_y__dcd_z*dcd_z__ddz
c bcd_z=-bc_x*cd_y + bc_y*cd_x
      dbcd_z__dbx=dbcd_z__dbc_x*dbc_x__dbx 
      dbcd_z__dby=dbcd_z__dbc_y*dbc_y__dby
c      dbcd_z__dbz=0
      dbcd_z__dcx=dbcd_z__dbc_x*dbc_x__dcx + dbcd_z__dcd_x*dcd_x__dcx
      dbcd_z__dcy=dbcd_z__dcd_y*dcd_y__dcy + dbcd_z__dbc_y*dbc_y__dcy
c      dbcd_z__dcz=0
      dbcd_z__ddx=dbcd_z__dcd_x*dcd_x__ddx
      dbcd_z__ddy=dbcd_z__dcd_y*dcd_y__ddy
c      dbcd_z__ddz=0

c
c abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      abc_length_inv_cub=abc_length_inv**3
      dabc_length_inv__dabc_x=-abc_x*abc_length_inv_cub
      dabc_length_inv__dabc_y=-abc_y*abc_length_inv_cub
      dabc_length_inv__dabc_z=-abc_z*abc_length_inv_cub

c bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
      bcd_length_inv_cub=bcd_length_inv**3
      dbcd_length_inv__dbcd_x=-bcd_x*bcd_length_inv_cub
      dbcd_length_inv__dbcd_y=-bcd_y*bcd_length_inv_cub
      dbcd_length_inv__dbcd_z=-bcd_z*bcd_length_inv_cub

c abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      dabc_length_inv__dax=dabc_length_inv__dabc_y*dabc_y__dax
     4 + dabc_length_inv__dabc_z*dabc_z__dax
      dabc_length_inv__day=dabc_length_inv__dabc_x*dabc_x__day
     4 + dabc_length_inv__dabc_z*dabc_z__day
      dabc_length_inv__daz=dabc_length_inv__dabc_x*dabc_x__daz
     3 + dabc_length_inv__dabc_y*dabc_y__daz
      dabc_length_inv__dbx=dabc_length_inv__dabc_y*dabc_y__dbx
     4 + dabc_length_inv__dabc_z*dabc_z__dbx
      dabc_length_inv__dby=dabc_length_inv__dabc_x*dabc_x__dby
     4 + dabc_length_inv__dabc_z*dabc_z__dby
      dabc_length_inv__dbz=dabc_length_inv__dabc_x*dabc_x__dbz
     3 + dabc_length_inv__dabc_y*dabc_y__dbz
      dabc_length_inv__dcx=dabc_length_inv__dabc_y*dabc_y__dcx
     4 + dabc_length_inv__dabc_z*dabc_z__dcx
      dabc_length_inv__dcy=dabc_length_inv__dabc_x*dabc_x__dcy
     4 + dabc_length_inv__dabc_z*dabc_z__dcy
      dabc_length_inv__dcz=dabc_length_inv__dabc_x*dabc_x__dcz
     3 + dabc_length_inv__dabc_y*dabc_y__dcz

c bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
c derivatives according to dax, day, daz
      dbcd_length_inv__dbx=dbcd_length_inv__dbcd_y*dbcd_y__dbx
     4 + dbcd_length_inv__dbcd_z*dbcd_z__dbx
      dbcd_length_inv__dby=dbcd_length_inv__dbcd_x*dbcd_x__dby
     4 + dbcd_length_inv__dbcd_z*dbcd_z__dby
      dbcd_length_inv__dbz=dbcd_length_inv__dbcd_x*dbcd_x__dbz
     3 + dbcd_length_inv__dbcd_y*dbcd_y__dbz
      dbcd_length_inv__dcx=dbcd_length_inv__dbcd_y*dbcd_y__dcx
     4 + dbcd_length_inv__dbcd_z*dbcd_z__dcx
      dbcd_length_inv__dcy=dbcd_length_inv__dbcd_x*dbcd_x__dcy
     4 + dbcd_length_inv__dbcd_z*dbcd_z__dcy
      dbcd_length_inv__dcz=dbcd_length_inv__dbcd_x*dbcd_x__dcz
     3 + dbcd_length_inv__dbcd_y*dbcd_y__dcz
      dbcd_length_inv__ddx=dbcd_length_inv__dbcd_y*dbcd_y__ddx
     4 + dbcd_length_inv__dbcd_z*dbcd_z__ddx
      dbcd_length_inv__ddy=dbcd_length_inv__dbcd_x*dbcd_x__ddy
     4 + dbcd_length_inv__dbcd_z*dbcd_z__ddy
      dbcd_length_inv__ddz=dbcd_length_inv__dbcd_x*dbcd_x__ddz
     3 + dbcd_length_inv__dbcd_y*dbcd_y__ddz

c abc1_x=abc_x*abc_length_inv
c abc1_y=abc_y*abc_length_inv
c abc1_z=abc_z*abc_length_inv
      dabc1_x__dabc_x=abc_length_inv
      dabc1_y__dabc_y=abc_length_inv
      dabc1_z__dabc_z=abc_length_inv

      dabc1_x__dabc_length_inv=abc_x
      dabc1_y__dabc_length_inv=abc_y
      dabc1_z__dabc_length_inv=abc_z

c derivation of the components of the normals
c abc1_x=abc_x*abc_length_inv
c abc1_y=abc_y*abc_length_inv
c abc1_z=abc_z*abc_length_inv
      dabc1_x__dax=
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dax
      dabc1_y__dax=dabc1_y__dabc_y*dabc_y__dax +
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dax
      dabc1_z__dax=dabc1_z__dabc_z*dabc_z__dax +
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dax
      dabc1_x__day=dabc1_x__dabc_x*dabc_x__day +
     2 dabc1_x__dabc_length_inv*dabc_length_inv__day
      dabc1_y__day=
     2 dabc1_y__dabc_length_inv*dabc_length_inv__day
      dabc1_z__day=dabc1_z__dabc_z*dabc_z__day +
     2 dabc1_z__dabc_length_inv*dabc_length_inv__day
      dabc1_x__daz=dabc1_x__dabc_x*dabc_x__daz +
     2 dabc1_x__dabc_length_inv*dabc_length_inv__daz
      dabc1_y__daz=dabc1_y__dabc_y*dabc_y__daz +
     2 dabc1_y__dabc_length_inv*dabc_length_inv__daz
      dabc1_z__daz=
     2 dabc1_z__dabc_length_inv*dabc_length_inv__daz

      dabc1_x__dbx=
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dbx
      dabc1_y__dbx=dabc1_y__dabc_y*dabc_y__dbx +
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dbx
      dabc1_z__dbx=dabc1_z__dabc_z*dabc_z__dbx +
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dbx
      dabc1_x__dby=dabc1_x__dabc_x*dabc_x__dby +
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dby
      dabc1_y__dby=
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dby
      dabc1_z__dby=dabc1_z__dabc_z*dabc_z__dby +
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dby
      dabc1_x__dbz=dabc1_x__dabc_x*dabc_x__dbz +
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dbz
      dabc1_y__dbz=dabc1_y__dabc_y*dabc_y__dbz +
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dbz
      dabc1_z__dbz=
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dbz

      dabc1_x__dcx=
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dcx
      dabc1_y__dcx=dabc1_y__dabc_y*dabc_y__dcx +
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dcx
      dabc1_z__dcx=dabc1_z__dabc_z*dabc_z__dcx +
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dcx
      dabc1_x__dcy=dabc1_x__dabc_x*dabc_x__dcy +
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dcy
      dabc1_y__dcy=
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dcy
      dabc1_z__dcy=dabc1_z__dabc_z*dabc_z__dcy +
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dcy
      dabc1_x__dcz=dabc1_x__dabc_x*dabc_x__dcz +
     2 dabc1_x__dabc_length_inv*dabc_length_inv__dcz
      dabc1_y__dcz=dabc1_y__dabc_y*dabc_y__dcz +
     2 dabc1_y__dabc_length_inv*dabc_length_inv__dcz
      dabc1_z__dcz=
     2 dabc1_z__dabc_length_inv*dabc_length_inv__dcz

c      dabc1_x__ddx=0
c      dabc1_y__ddx=0
c      dabc1_z__ddx=0
c      dabc1_x__ddy=0
c      dabc1_y__ddy=0
c      dabc1_z__ddy=0
c      dabc1_x__ddz=0
c      dabc1_y__ddz=0
c      dabc1_z__ddz=0

c bcd1_x=bcd_x*bcd_length_inv
c bcd1_y=bcd_y*bcd_length_inv
c bcd1_z=bcd_z*bcd_length_inv
      dbcd1_x__dbcd_x=bcd_length_inv
      dbcd1_y__dbcd_y=bcd_length_inv
      dbcd1_z__dbcd_z=bcd_length_inv

      dbcd1_x__dbcd_length_inv=bcd_x
      dbcd1_y__dbcd_length_inv=bcd_y
      dbcd1_z__dbcd_length_inv=bcd_z

c bcd1_x=bcd_x*bcd_length_inv
c bcd1_y=bcd_y*bcd_length_inv
c bcd1_z=bcd_z*bcd_length_inv

c      dbcd1_x__dax=0
c      dbcd1_y__dax=0
c      dbcd1_z__dax=0
c      dbcd1_x__day=0
c      dbcd1_y__day=0
c      dbcd1_z__day=0
c      dbcd1_x__daz=0
c      dbcd1_y__daz=0
c      dbcd1_z__daz=0

      dbcd1_x__dbx=
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dbx
      dbcd1_y__dbx=dbcd1_y__dbcd_y*dbcd_y__dbx +
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dbx
      dbcd1_z__dbx=dbcd1_z__dbcd_z*dbcd_z__dbx +
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dbx
      dbcd1_x__dby=dbcd1_x__dbcd_x*dbcd_x__dby +
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dby
      dbcd1_y__dby=
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dby
      dbcd1_z__dby=dbcd1_z__dbcd_z*dbcd_z__dby +
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dby
      dbcd1_x__dbz=dbcd1_x__dbcd_x*dbcd_x__dbz +
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dbz
      dbcd1_y__dbz=dbcd1_y__dbcd_y*dbcd_y__dbz +
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dbz
      dbcd1_z__dbz=
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dbz

      dbcd1_x__dcx=
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dcx
      dbcd1_y__dcx=dbcd1_y__dbcd_y*dbcd_y__dcx +
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dcx
      dbcd1_z__dcx=dbcd1_z__dbcd_z*dbcd_z__dcx +
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dcx
      dbcd1_x__dcy=dbcd1_x__dbcd_x*dbcd_x__dcy +
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dcy
      dbcd1_y__dcy=
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dcy
      dbcd1_z__dcy=dbcd1_z__dbcd_z*dbcd_z__dcy +
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dcy
      dbcd1_x__dcz=dbcd1_x__dbcd_x*dbcd_x__dcz +
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dcz
      dbcd1_y__dcz=dbcd1_y__dbcd_y*dbcd_y__dcz +
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dcz
      dbcd1_z__dcz=
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dcz

      dbcd1_x__ddx=
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__ddx
      dbcd1_y__ddx=dbcd1_y__dbcd_y*dbcd_y__ddx +
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__ddx
      dbcd1_z__ddx=dbcd1_z__dbcd_z*dbcd_z__ddx +
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__ddx
      dbcd1_x__ddy=dbcd1_x__dbcd_x*dbcd_x__ddy +
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__ddy
      dbcd1_y__ddy=
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__ddy
      dbcd1_z__ddy=dbcd1_z__dbcd_z*dbcd_z__ddy +
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__ddy
      dbcd1_x__ddz=dbcd1_x__dbcd_x*dbcd_x__ddz +
     2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__ddz
      dbcd1_y__ddz=dbcd1_y__dbcd_y*dbcd_y__ddz +
     2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__ddz
      dbcd1_z__ddz=
     2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__ddz

c aux_x=abc1_y*bc1_z-bc1_y*abc1_z
c aux_y=abc1_z*bc1_x-bc1_z*abc1_x
c aux_z=abc1_x*bc1_y-bc1_x*abc1_y
      daux_x__dabc1_y=bc1_z
      daux_x__dabc1_z=-bc1_y
      daux_y__dabc1_z=bc1_x
      daux_y__dabc1_x=-bc1_z
      daux_z__dabc1_x=bc1_y
      daux_z__dabc1_y=-bc1_x

      daux_x__dbc1_z=abc1_y
      daux_x__dbc1_y=-abc1_z
      daux_y__dbc1_x=abc1_z
      daux_y__dbc1_z=-abc1_x
      daux_z__dbc1_y=abc1_x
      daux_z__dbc1_x=-abc1_y

c aux_x=abc1_y*bc1_z-bc1_y*abc1_z
      daux_x__dax=
     2 daux_x__dabc1_y*dabc1_y__dax + daux_x__dabc1_z*dabc1_z__dax
      daux_x__day=
     2 daux_x__dabc1_y*dabc1_y__day + daux_x__dabc1_z*dabc1_z__day
      daux_x__daz=
     2 daux_x__dabc1_y*dabc1_y__daz + daux_x__dabc1_z*dabc1_z__daz
      daux_x__dbx=
     2 daux_x__dabc1_y*dabc1_y__dbx + daux_x__dbc1_z*dbc1_z__dbx
     4 + daux_x__dbc1_y*dbc1_y__dbx + daux_x__dabc1_z*dabc1_z__dbx
      daux_x__dby=
     2 daux_x__dabc1_y*dabc1_y__dby + daux_x__dbc1_z*dbc1_z__dby
     4 + daux_x__dbc1_y*dbc1_y__dby + daux_x__dabc1_z*dabc1_z__dby
      daux_x__dbz=
     2 daux_x__dabc1_y*dabc1_y__dbz + daux_x__dbc1_z*dbc1_z__dbz
     4 + daux_x__dbc1_y*dbc1_y__dbz + daux_x__dabc1_z*dabc1_z__dbz
      daux_x__dcx=
     2 daux_x__dabc1_y*dabc1_y__dcx + daux_x__dbc1_z*dbc1_z__dcx
     4 + daux_x__dbc1_y*dbc1_y__dcx + daux_x__dabc1_z*dabc1_z__dcx
      daux_x__dcy=
     2 daux_x__dabc1_y*dabc1_y__dcy + daux_x__dbc1_z*dbc1_z__dcy
     4 + daux_x__dbc1_y*dbc1_y__dcy + daux_x__dabc1_z*dabc1_z__dcy
      daux_x__dcz=
     2 daux_x__dabc1_y*dabc1_y__dcz + daux_x__dbc1_z*dbc1_z__dcz
     4 + daux_x__dbc1_y*dbc1_y__dcz + daux_x__dabc1_z*dabc1_z__dcz
c      daux_x__ddx=0
c      daux_x__ddy=0
c      daux_x__ddz=0
c aux_y=abc1_z*bc1_x-bc1_z*abc1_x
      daux_y__dax=
     2 daux_y__dabc1_z*dabc1_z__dax + daux_y__dabc1_x*dabc1_x__dax
      daux_y__day=
     2 daux_y__dabc1_z*dabc1_z__day + daux_y__dabc1_x*dabc1_x__day
      daux_y__daz=
     2 daux_y__dabc1_z*dabc1_z__daz + daux_y__dabc1_x*dabc1_x__daz
      daux_y__dbx=
     2 daux_y__dabc1_z*dabc1_z__dbx + daux_y__dbc1_x*dbc1_x__dbx
     4 + daux_y__dbc1_z*dbc1_z__dbx + daux_y__dabc1_x*dabc1_x__dbx
      daux_y__dby=
     2 daux_y__dabc1_z*dabc1_z__dby + daux_y__dbc1_x*dbc1_x__dby
     4 + daux_y__dbc1_z*dbc1_z__dby + daux_y__dabc1_x*dabc1_x__dby
      daux_y__dbz=
     2 daux_y__dabc1_z*dabc1_z__dbz + daux_y__dbc1_x*dbc1_x__dbz
     4 + daux_y__dbc1_z*dbc1_z__dbz + daux_y__dabc1_x*dabc1_x__dbz
      daux_y__dcx=
     2 daux_y__dabc1_z*dabc1_z__dcx + daux_y__dbc1_x*dbc1_x__dcx
     4 + daux_y__dbc1_z*dbc1_z__dcx + daux_y__dabc1_x*dabc1_x__dcx
      daux_y__dcy=
     2 daux_y__dabc1_z*dabc1_z__dcy + daux_y__dbc1_x*dbc1_x__dcy
     4 + daux_y__dbc1_z*dbc1_z__dcy + daux_y__dabc1_x*dabc1_x__dcy
      daux_y__dcz=
     2 daux_y__dabc1_z*dabc1_z__dcz + daux_y__dbc1_x*dbc1_x__dcz
     4 + daux_y__dbc1_z*dbc1_z__dcz + daux_y__dabc1_x*dabc1_x__dcz
c      daux_y__ddx=0
c      daux_y__ddy=0
c      daux_y__ddz=0
c aux_z=abc1_x*bc1_y-bc1_x*abc1_y
      daux_z__dax=
     2 daux_z__dabc1_x*dabc1_x__dax + daux_z__dabc1_y*dabc1_y__dax
      daux_z__day=
     2 daux_z__dabc1_x*dabc1_x__day + daux_z__dabc1_y*dabc1_y__day
      daux_z__daz=
     2 daux_z__dabc1_x*dabc1_x__daz + daux_z__dabc1_y*dabc1_y__daz
      daux_z__dbx=
     2 daux_z__dabc1_x*dabc1_x__dbx + daux_z__dbc1_y*dbc1_y__dbx
     4 + daux_z__dbc1_x*dbc1_x__dbx + daux_z__dabc1_y*dabc1_y__dbx
      daux_z__dby=
     2 daux_z__dabc1_x*dabc1_x__dby + daux_z__dbc1_y*dbc1_y__dby
     4 + daux_z__dbc1_x*dbc1_x__dby + daux_z__dabc1_y*dabc1_y__dby
      daux_z__dbz=
     2 daux_z__dabc1_x*dabc1_x__dbz + daux_z__dbc1_y*dbc1_y__dbz
     4 + daux_z__dbc1_x*dbc1_x__dbz + daux_z__dabc1_y*dabc1_y__dbz
      daux_z__dcx=
     2 daux_z__dabc1_x*dabc1_x__dcx + daux_z__dbc1_y*dbc1_y__dcx
     4 + daux_z__dbc1_x*dbc1_x__dcx + daux_z__dabc1_y*dabc1_y__dcx
      daux_z__dcy=
     2 daux_z__dabc1_x*dabc1_x__dcy + daux_z__dbc1_y*dbc1_y__dcy
     4 + daux_z__dbc1_x*dbc1_x__dcy + daux_z__dabc1_y*dabc1_y__dcy
      daux_z__dcz=
     2 daux_z__dabc1_x*dabc1_x__dcz + daux_z__dbc1_y*dbc1_y__dcz
     4 + daux_z__dbc1_x*dbc1_x__dcz + daux_z__dabc1_y*dabc1_y__dcz
c      daux_z__ddx=0
c      daux_z__ddy=0
c      daux_z__ddz=0

c y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
c x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
      dy__daux_x=bcd1_x
      dy__daux_y=bcd1_y
      dy__daux_z=bcd1_z

      dy__dbcd1_x=aux_x
      dy__dbcd1_y=aux_y
      dy__dbcd1_z=aux_z

      dx__dabc1_x=bcd1_x
      dx__dabc1_y=bcd1_y
      dx__dabc1_z=bcd1_z

      dx__dbcd1_x=abc1_x
      dx__dbcd1_y=abc1_y
      dx__dbcd1_z=abc1_z

c derivation of y
c y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
      dy__dax=
     2 dy__daux_x*daux_x__dax + 
     3 dy__daux_y*daux_y__dax + 
     4 dy__daux_z*daux_z__dax  
      dy__day=
     2 dy__daux_x*daux_x__day + 
     3 dy__daux_y*daux_y__day + 
     4 dy__daux_z*daux_z__day 
      dy__daz=
     2 dy__daux_x*daux_x__daz + 
     3 dy__daux_y*daux_y__daz + 
     4 dy__daux_z*daux_z__daz 
      dy__dbx=
     2 dy__daux_x*daux_x__dbx + dy__dbcd1_x*dbcd1_x__dbx +
     3 dy__daux_y*daux_y__dbx + dy__dbcd1_y*dbcd1_y__dbx +
     4 dy__daux_z*daux_z__dbx + dy__dbcd1_z*dbcd1_z__dbx
      dy__dby=
     2 dy__daux_x*daux_x__dby + dy__dbcd1_x*dbcd1_x__dby +
     3 dy__daux_y*daux_y__dby + dy__dbcd1_y*dbcd1_y__dby +
     4 dy__daux_z*daux_z__dby + dy__dbcd1_z*dbcd1_z__dby
      dy__dbz=
     2 dy__daux_x*daux_x__dbz + dy__dbcd1_x*dbcd1_x__dbz +
     3 dy__daux_y*daux_y__dbz + dy__dbcd1_y*dbcd1_y__dbz +
     4 dy__daux_z*daux_z__dbz + dy__dbcd1_z*dbcd1_z__dbz
      dy__dcx=
     2 dy__daux_x*daux_x__dcx + dy__dbcd1_x*dbcd1_x__dcx +
     3 dy__daux_y*daux_y__dcx + dy__dbcd1_y*dbcd1_y__dcx +
     4 dy__daux_z*daux_z__dcx + dy__dbcd1_z*dbcd1_z__dcx
      dy__dcy=
     2 dy__daux_x*daux_x__dcy + dy__dbcd1_x*dbcd1_x__dcy +
     3 dy__daux_y*daux_y__dcy + dy__dbcd1_y*dbcd1_y__dcy +
     4 dy__daux_z*daux_z__dcy + dy__dbcd1_z*dbcd1_z__dcy
      dy__dcz=
     2 dy__daux_x*daux_x__dcz + dy__dbcd1_x*dbcd1_x__dcz +
     3 dy__daux_y*daux_y__dcz + dy__dbcd1_y*dbcd1_y__dcz +
     4 dy__daux_z*daux_z__dcz + dy__dbcd1_z*dbcd1_z__dcz
      dy__ddx=
     2 dy__dbcd1_x*dbcd1_x__ddx +
     3 dy__dbcd1_y*dbcd1_y__ddx +
     4 dy__dbcd1_z*dbcd1_z__ddx
      dy__ddy=
     2 dy__dbcd1_x*dbcd1_x__ddy +
     3 dy__dbcd1_y*dbcd1_y__ddy +
     4 dy__dbcd1_z*dbcd1_z__ddy
      dy__ddz=
     2 dy__dbcd1_x*dbcd1_x__ddz +
     3 dy__dbcd1_y*dbcd1_y__ddz +
     4 dy__dbcd1_z*dbcd1_z__ddz

c derivation of x
c x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
      dx__dax=
     2 dx__dabc1_x*dabc1_x__dax + 
     3 dx__dabc1_y*dabc1_y__dax + 
     4 dx__dabc1_z*dabc1_z__dax  
      dx__day=
     2 dx__dabc1_x*dabc1_x__day +
     3 dx__dabc1_y*dabc1_y__day + 
     4 dx__dabc1_z*dabc1_z__day 
      dx__daz=
     2 dx__dabc1_x*dabc1_x__daz + 
     3 dx__dabc1_y*dabc1_y__daz + 
     4 dx__dabc1_z*dabc1_z__daz 
      dx__dbx=
     2 dx__dabc1_x*dabc1_x__dbx + dx__dbcd1_x*dbcd1_x__dbx +
     3 dx__dabc1_y*dabc1_y__dbx + dx__dbcd1_y*dbcd1_y__dbx +
     4 dx__dabc1_z*dabc1_z__dbx + dx__dbcd1_z*dbcd1_z__dbx
      dx__dby=
     2 dx__dabc1_x*dabc1_x__dby + dx__dbcd1_x*dbcd1_x__dby +
     3 dx__dabc1_y*dabc1_y__dby + dx__dbcd1_y*dbcd1_y__dby +
     4 dx__dabc1_z*dabc1_z__dby + dx__dbcd1_z*dbcd1_z__dby
      dx__dbz=
     2 dx__dabc1_x*dabc1_x__dbz + dx__dbcd1_x*dbcd1_x__dbz +
     3 dx__dabc1_y*dabc1_y__dbz + dx__dbcd1_y*dbcd1_y__dbz +
     4 dx__dabc1_z*dabc1_z__dbz + dx__dbcd1_z*dbcd1_z__dbz
      dx__dcx=
     2 dx__dabc1_x*dabc1_x__dcx + dx__dbcd1_x*dbcd1_x__dcx +
     3 dx__dabc1_y*dabc1_y__dcx + dx__dbcd1_y*dbcd1_y__dcx +
     4 dx__dabc1_z*dabc1_z__dcx + dx__dbcd1_z*dbcd1_z__dcx
      dx__dcy=
     2 dx__dabc1_x*dabc1_x__dcy + dx__dbcd1_x*dbcd1_x__dcy +
     3 dx__dabc1_y*dabc1_y__dcy + dx__dbcd1_y*dbcd1_y__dcy +
     4 dx__dabc1_z*dabc1_z__dcy + dx__dbcd1_z*dbcd1_z__dcy
      dx__dcz=
     2 dx__dabc1_x*dabc1_x__dcz + dx__dbcd1_x*dbcd1_x__dcz +
     3 dx__dabc1_y*dabc1_y__dcz + dx__dbcd1_y*dbcd1_y__dcz +
     4 dx__dabc1_z*dabc1_z__dcz + dx__dbcd1_z*dbcd1_z__dcz
      dx__ddx=
     2 dx__dbcd1_x*dbcd1_x__ddx +
     3 dx__dbcd1_y*dbcd1_y__ddx +
     4 dx__dbcd1_z*dbcd1_z__ddx
      dx__ddy=
     2 dx__dbcd1_x*dbcd1_x__ddy +
     3 dx__dbcd1_y*dbcd1_y__ddy +
     4 dx__dbcd1_z*dbcd1_z__ddy
      dx__ddz=
     2 dx__dbcd1_x*dbcd1_x__ddz +
     3 dx__dbcd1_y*dbcd1_y__ddz +
     4 dx__dbcd1_z*dbcd1_z__ddz

c derivation atan2(y,x) according to x and y
      df__dx=-y/(x**2 + y**2)
      df__dy=x/(x**2 + y**2)

c derive f according to all 12 cartesion components of the four points
c f=atan2(y,x)
      df__dax=df__dx*dx__dax + df__dy*dy__dax
      df__day=df__dx*dx__day + df__dy*dy__day
      df__daz=df__dx*dx__daz + df__dy*dy__daz
      df__dbx=df__dx*dx__dbx + df__dy*dy__dbx
      df__dby=df__dx*dx__dby + df__dy*dy__dby
      df__dbz=df__dx*dx__dbz + df__dy*dy__dbz
      df__dcx=df__dx*dx__dcx + df__dy*dy__dcx
      df__dcy=df__dx*dx__dcy + df__dy*dy__dcy
      df__dcz=df__dx*dx__dcz + df__dy*dy__dcz
      df__ddx=df__dx*dx__ddx + df__dy*dy__ddx
      df__ddy=df__dx*dx__ddy + df__dy*dy__ddy
      df__ddz=df__dx*dx__ddz + df__dy*dy__ddz

c actually it would have to be df__dax etc
      dax=df__dax
      day=df__day
      daz=df__daz
      dbx=df__dbx
      dby=df__dby
      dbz=df__dbz
      dcx=df__dcx
      dcy=df__dcy
      dcz=df__dcz
      ddx=df__ddx
      ddy=df__ddy
      ddz=df__ddz

c f is the function atan2(y, x)
c fortran (and most other sources) use 'atan2(y,x)' while mathematica uses 'atan2(x,y)'

      return
      END

