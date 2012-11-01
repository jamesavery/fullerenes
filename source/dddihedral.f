      SUBROUTINE DDDIHEDRAL(
     1 ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz,
     1 df__dax, df__day, df__daz, df__dbx, df__dby, df__dbz, df__dcx,
     1 df__dcy, df__dcz, df__ddx, df__ddy, df__ddz,
     1 ddf11dax__dax, ddf11dax__day, ddf11dax__daz, ddf11dax__dbx,
     1 ddf11dax__dby, ddf11dax__dbz, ddf11dax__dcx, ddf11dax__dcy,
     1 ddf11dax__dcz, ddf11dax__ddx, ddf11dax__ddy, ddf11dax__ddz,
     1 ddf11day__day, ddf11day__daz, ddf11day__dbx, ddf11day__dby,
     1 ddf11day__dbz, ddf11day__dcx, ddf11day__dcy, ddf11day__dcz,
     1 ddf11day__ddx, ddf11day__ddy, ddf11day__ddz, ddf11daz__daz,
     1 ddf11daz__dbx, ddf11daz__dby, ddf11daz__dbz, ddf11daz__dcx,
     1 ddf11daz__dcy, ddf11daz__dcz, ddf11daz__ddx, ddf11daz__ddy,
     1 ddf11daz__ddz, ddf11dbx__dbx, ddf11dbx__dby, ddf11dbx__dbz,
     1 ddf11dbx__dcx, ddf11dbx__dcy, ddf11dbx__dcz, ddf11dbx__ddx,
     1 ddf11dbx__ddy, ddf11dbx__ddz, ddf11dby__dby, ddf11dby__dbz,
     1 ddf11dby__dcx, ddf11dby__dcy, ddf11dby__dcz, ddf11dby__ddx,
     1 ddf11dby__ddy, ddf11dby__ddz, ddf11dbz__dbz, ddf11dbz__dcx,
     1 ddf11dbz__dcy, ddf11dbz__dcz, ddf11dbz__ddx, ddf11dbz__ddy,
     1 ddf11dbz__ddz, ddf11dcx__dcx, ddf11dcx__dcy, ddf11dcx__dcz,
     1 ddf11dcx__ddx, ddf11dcx__ddy, ddf11dcx__ddz, ddf11dcy__dcy,
     1 ddf11dcy__dcz, ddf11dcy__ddx, ddf11dcy__ddy, ddf11dcy__ddz,
     1 ddf11dcz__dcz, ddf11dcz__ddx, ddf11dcz__ddy, ddf11dcz__ddz,
     1 ddf11ddx__ddx, ddf11ddx__ddy, ddf11ddx__ddz, ddf11ddy__ddy,
     1 ddf11ddy__ddz, ddf11ddz__ddz,
     1 dihedral_abcd)
      IMPLICIT REAL*8 (a-z)

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

C THE FISRT DERIVATIVES
c to be read from bottom to top

c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      bc_length_inv_cub=bc_length_inv**3
      dbc_length_inv__dbc_x=-bc_x*bc_length_inv_cub
      dbc_length_inv__dbc_y=-bc_y*bc_length_inv_cub
      dbc_length_inv__dbc_z=-bc_z*bc_length_inv_cub

c bc1_x=bc_x*bc_length_inv
      dbc1_x__dbx=
     2 bc_length_inv + bc_x*dbc_length_inv__dbc_x
      dbc1_x__dby=
     2 bc_x*dbc_length_inv__dbc_y
      dbc1_x__dbz=
     2 bc_x*dbc_length_inv__dbc_z
      dbc1_x__dcx=
     2 -bc_length_inv - bc_x*dbc_length_inv__dbc_x
      dbc1_x__dcy=
     2 -bc_x*dbc_length_inv__dbc_y
      dbc1_x__dcz=
     2 -bc_x*dbc_length_inv__dbc_z
c bc1_y=bc_y*bc_length_inv
      dbc1_y__dbx=
     2 bc_y*dbc_length_inv__dbc_x
      dbc1_y__dby=
     2 bc_length_inv + bc_y*dbc_length_inv__dbc_y
      dbc1_y__dbz=
     2 bc_y*dbc_length_inv__dbc_z
      dbc1_y__dcx=
     2 -bc_y*dbc_length_inv__dbc_x
      dbc1_y__dcy=
     2 -bc_length_inv - bc_y*dbc_length_inv__dbc_y
      dbc1_y__dcz=
     2 -bc_y*dbc_length_inv__dbc_z
c bc1_z=bc_z*bc_length_inv
      dbc1_z__dbx=
     2 bc_z*dbc_length_inv__dbc_x
      dbc1_z__dby=
     2 bc_z*dbc_length_inv__dbc_y
      dbc1_z__dbz=
     2 bc_length_inv + bc_z*dbc_length_inv__dbc_z
      dbc1_z__dcx=
     2 -bc_z*dbc_length_inv__dbc_x
      dbc1_z__dcy=
     2 -bc_z*dbc_length_inv__dbc_y
      dbc1_z__dcz=
     2 -bc_length_inv - bc_z*dbc_length_inv__dbc_z

c abc_x=-ab_y*bc_z + ab_z*bc_y
      dabc_x__dby=bc_z + ab_z
      dabc_x__dbz=-ab_y - bc_y
c abc_y=-ab_z*bc_x + ab_x*bc_z
      dabc_y__dbx=-ab_z - bc_z
      dabc_y__dbz=bc_x + ab_x
c abc_z=-ab_x*bc_y + ab_y*bc_x
      dabc_z__dbx=ab_y + bc_y
      dabc_z__dby=-ab_x - bc_x 
c bcd_x=-bc_y*cd_z + bc_z*cd_y
      dbcd_x__dcy=cd_z + bc_z
      dbcd_x__dcz=-bc_y - cd_y 
c bcd_y=-bc_z*cd_x + bc_x*cd_z
      dbcd_y__dcx=-bc_z - cd_z 
      dbcd_y__dcz=cd_x + bc_x
c bcd_z=-bc_x*cd_y + bc_y*cd_x
      dbcd_z__dcx=cd_y + bc_y
      dbcd_z__dcy=-bc_x - cd_x

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
      dabc_length_inv__dax=dabc_length_inv__dabc_y*bc_z
     4 - dabc_length_inv__dabc_z*bc_y
      dabc_length_inv__day=-dabc_length_inv__dabc_x*bc_z
     4 + dabc_length_inv__dabc_z*bc_x
      dabc_length_inv__daz=dabc_length_inv__dabc_x*bc_y
     3 - dabc_length_inv__dabc_y*bc_x
      dabc_length_inv__dbx=dabc_length_inv__dabc_y*dabc_y__dbx
     4 + dabc_length_inv__dabc_z*dabc_z__dbx
      dabc_length_inv__dby=dabc_length_inv__dabc_x*dabc_x__dby
     4 + dabc_length_inv__dabc_z*dabc_z__dby
      dabc_length_inv__dbz=dabc_length_inv__dabc_x*dabc_x__dbz
     3 + dabc_length_inv__dabc_y*dabc_y__dbz
      dabc_length_inv__dcx=dabc_length_inv__dabc_y*ab_z
     4 - dabc_length_inv__dabc_z*ab_y
      dabc_length_inv__dcy=-dabc_length_inv__dabc_x*ab_z
     4 + dabc_length_inv__dabc_z*ab_x
      dabc_length_inv__dcz=dabc_length_inv__dabc_x*ab_y
     3 - dabc_length_inv__dabc_y*ab_x

c bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
c derivatives according to dax, day, daz
      dbcd_length_inv__dbx=dbcd_length_inv__dbcd_y*cd_z
     4 - dbcd_length_inv__dbcd_z*cd_y
      dbcd_length_inv__dby=-dbcd_length_inv__dbcd_x*cd_z
     4 + dbcd_length_inv__dbcd_z*cd_x
      dbcd_length_inv__dbz=dbcd_length_inv__dbcd_x*cd_y
     3 - dbcd_length_inv__dbcd_y*cd_x
      dbcd_length_inv__dcx=dbcd_length_inv__dbcd_y*dbcd_y__dcx
     4 + dbcd_length_inv__dbcd_z*dbcd_z__dcx
      dbcd_length_inv__dcy=dbcd_length_inv__dbcd_x*dbcd_x__dcy
     4 + dbcd_length_inv__dbcd_z*dbcd_z__dcy
      dbcd_length_inv__dcz=dbcd_length_inv__dbcd_x*dbcd_x__dcz
     3 + dbcd_length_inv__dbcd_y*dbcd_y__dcz
      dbcd_length_inv__ddx=dbcd_length_inv__dbcd_y*bc_z
     4 - dbcd_length_inv__dbcd_z*bc_y
      dbcd_length_inv__ddy=-dbcd_length_inv__dbcd_x*bc_z
     4 + dbcd_length_inv__dbcd_z*bc_x
      dbcd_length_inv__ddz=dbcd_length_inv__dbcd_x*bc_y
     3 - dbcd_length_inv__dbcd_y*bc_x


c derivation of the components of the normals
c abc1_x=abc_x*abc_length_inv
c abc1_y=abc_y*abc_length_inv
c abc1_z=abc_z*abc_length_inv
      dabc1_x__dax=
     2 abc_x*dabc_length_inv__dax
      dabc1_y__dax=abc_length_inv*bc_z +
     2 abc_y*dabc_length_inv__dax
      dabc1_z__dax=-abc_length_inv*bc_y +
     2 abc_z*dabc_length_inv__dax
      dabc1_x__day=- abc_length_inv*bc_z +
     2 abc_x*dabc_length_inv__day
      dabc1_y__day=
     2 abc_y*dabc_length_inv__day
      dabc1_z__day=abc_length_inv*bc_x +
     2 abc_z*dabc_length_inv__day
      dabc1_x__daz=abc_length_inv*bc_y +
     2 abc_x*dabc_length_inv__daz
      dabc1_y__daz=-abc_length_inv*bc_x +
     2 abc_y*dabc_length_inv__daz
      dabc1_z__daz=
     2 abc_z*dabc_length_inv__daz

      dabc1_x__dbx=
     2 abc_x*dabc_length_inv__dbx
      dabc1_y__dbx=abc_length_inv*dabc_y__dbx +
     2 abc_y*dabc_length_inv__dbx
      dabc1_z__dbx=abc_length_inv*dabc_z__dbx +
     2 abc_z*dabc_length_inv__dbx
      dabc1_x__dby=abc_length_inv*dabc_x__dby +
     2 abc_x*dabc_length_inv__dby
      dabc1_y__dby=
     2 abc_y*dabc_length_inv__dby
      dabc1_z__dby=abc_length_inv*dabc_z__dby +
     2 abc_z*dabc_length_inv__dby
      dabc1_x__dbz=abc_length_inv*dabc_x__dbz +
     2 abc_x*dabc_length_inv__dbz
      dabc1_y__dbz=abc_length_inv*dabc_y__dbz +
     2 abc_y*dabc_length_inv__dbz
      dabc1_z__dbz=
     2 abc_z*dabc_length_inv__dbz

      dabc1_x__dcx=
     2 abc_x*dabc_length_inv__dcx
      dabc1_y__dcx=abc_length_inv*ab_z +
     2 abc_y*dabc_length_inv__dcx
      dabc1_z__dcx=-abc_length_inv*ab_y +
     2 abc_z*dabc_length_inv__dcx
      dabc1_x__dcy=-abc_length_inv*ab_z +
     2 abc_x*dabc_length_inv__dcy
      dabc1_y__dcy=
     2 abc_y*dabc_length_inv__dcy
      dabc1_z__dcy=abc_length_inv*ab_x +
     2 abc_z*dabc_length_inv__dcy
      dabc1_x__dcz=abc_length_inv*ab_y +
     2 abc_x*dabc_length_inv__dcz
      dabc1_y__dcz=-abc_length_inv*ab_x +
     2 abc_y*dabc_length_inv__dcz
      dabc1_z__dcz=
     2 abc_z*dabc_length_inv__dcz


c bcd1_x=bcd_x*bcd_length_inv
c bcd1_y=bcd_y*bcd_length_inv
c bcd1_z=bcd_z*bcd_length_inv

      dbcd1_x__dbx=
     2 bcd_x*dbcd_length_inv__dbx
      dbcd1_y__dbx=bcd_length_inv*cd_z +
     2 bcd_y*dbcd_length_inv__dbx
      dbcd1_z__dbx=-bcd_length_inv*cd_y +
     2 bcd_z*dbcd_length_inv__dbx
      dbcd1_x__dby=-bcd_length_inv*cd_z +
     2 bcd_x*dbcd_length_inv__dby
      dbcd1_y__dby=
     2 bcd_y*dbcd_length_inv__dby
      dbcd1_z__dby=bcd_length_inv*cd_x +
     2 bcd_z*dbcd_length_inv__dby
      dbcd1_x__dbz=bcd_length_inv*cd_y +
     2 bcd_x*dbcd_length_inv__dbz
      dbcd1_y__dbz=-bcd_length_inv*cd_x +
     2 bcd_y*dbcd_length_inv__dbz
      dbcd1_z__dbz=
     2 bcd_z*dbcd_length_inv__dbz

      dbcd1_x__dcx=
     2 bcd_x*dbcd_length_inv__dcx
      dbcd1_y__dcx=bcd_length_inv*dbcd_y__dcx +
     2 bcd_y*dbcd_length_inv__dcx
      dbcd1_z__dcx=bcd_length_inv*dbcd_z__dcx +
     2 bcd_z*dbcd_length_inv__dcx
      dbcd1_x__dcy=bcd_length_inv*dbcd_x__dcy +
     2 bcd_x*dbcd_length_inv__dcy
      dbcd1_y__dcy=
     2 bcd_y*dbcd_length_inv__dcy
      dbcd1_z__dcy=bcd_length_inv*dbcd_z__dcy +
     2 bcd_z*dbcd_length_inv__dcy
      dbcd1_x__dcz=bcd_length_inv*dbcd_x__dcz +
     2 bcd_x*dbcd_length_inv__dcz
      dbcd1_y__dcz=bcd_length_inv*dbcd_y__dcz +
     2 bcd_y*dbcd_length_inv__dcz
      dbcd1_z__dcz=
     2 bcd_z*dbcd_length_inv__dcz

      dbcd1_x__ddx=
     2 bcd_x*dbcd_length_inv__ddx
      dbcd1_y__ddx=bcd_length_inv*bc_z +
     2 bcd_y*dbcd_length_inv__ddx
      dbcd1_z__ddx=-bcd_length_inv*bc_y +
     2 bcd_z*dbcd_length_inv__ddx
      dbcd1_x__ddy=-bcd_length_inv*bc_z +
     2 bcd_x*dbcd_length_inv__ddy
      dbcd1_y__ddy=
     2 bcd_y*dbcd_length_inv__ddy
      dbcd1_z__ddy=bcd_length_inv*bc_x +
     2 bcd_z*dbcd_length_inv__ddy
      dbcd1_x__ddz=bcd_length_inv*bc_y +
     2 bcd_x*dbcd_length_inv__ddz
      dbcd1_y__ddz=-bcd_length_inv*bc_x +
     2 bcd_y*dbcd_length_inv__ddz
      dbcd1_z__ddz=
     2 bcd_z*dbcd_length_inv__ddz

      daux_x__dax=
     2 bc1_z*dabc1_y__dax - bc1_y*dabc1_z__dax
      daux_x__day=
     2 bc1_z*dabc1_y__day - bc1_y*dabc1_z__day
      daux_x__daz=
     2 bc1_z*dabc1_y__daz - bc1_y*dabc1_z__daz
      daux_x__dbx=
     2 bc1_z*dabc1_y__dbx + abc1_y*dbc1_z__dbx
     4 - abc1_z*dbc1_y__dbx - bc1_y*dabc1_z__dbx
      daux_x__dby=
     2 bc1_z*dabc1_y__dby + abc1_y*dbc1_z__dby
     4 - abc1_z*dbc1_y__dby - bc1_y*dabc1_z__dby
      daux_x__dbz=
     2 bc1_z*dabc1_y__dbz + abc1_y*dbc1_z__dbz
     4 - abc1_z*dbc1_y__dbz - bc1_y*dabc1_z__dbz
      daux_x__dcx=
     2 bc1_z*dabc1_y__dcx + abc1_y*dbc1_z__dcx
     4 - abc1_z*dbc1_y__dcx - bc1_y*dabc1_z__dcx
      daux_x__dcy=
     2 bc1_z*dabc1_y__dcy + abc1_y*dbc1_z__dcy
     4 - abc1_z*dbc1_y__dcy - bc1_y*dabc1_z__dcy
      daux_x__dcz=
     2 bc1_z*dabc1_y__dcz + abc1_y*dbc1_z__dcz
     4 - abc1_z*dbc1_y__dcz - bc1_y*dabc1_z__dcz
c aux_y=abc1_z*bc1_x-bc1_z*abc1_x
      daux_y__dax=
     2 bc1_x*dabc1_z__dax - bc1_z*dabc1_x__dax
      daux_y__day=
     2 bc1_x*dabc1_z__day - bc1_z*dabc1_x__day
      daux_y__daz=
     2 bc1_x*dabc1_z__daz - bc1_z*dabc1_x__daz
      daux_y__dbx=
     2 bc1_x*dabc1_z__dbx + abc1_z*dbc1_x__dbx
     4 - abc1_x*dbc1_z__dbx - bc1_z*dabc1_x__dbx
      daux_y__dby=
     2 bc1_x*dabc1_z__dby + abc1_z*dbc1_x__dby
     4 - abc1_x*dbc1_z__dby - bc1_z*dabc1_x__dby
      daux_y__dbz=
     2 bc1_x*dabc1_z__dbz + abc1_z*dbc1_x__dbz
     4 - abc1_x*dbc1_z__dbz - bc1_z*dabc1_x__dbz
      daux_y__dcx=
     2 bc1_x*dabc1_z__dcx + abc1_z*dbc1_x__dcx
     4 - abc1_x*dbc1_z__dcx - bc1_z*dabc1_x__dcx
      daux_y__dcy=
     2 bc1_x*dabc1_z__dcy + abc1_z*dbc1_x__dcy
     4 - abc1_x*dbc1_z__dcy - bc1_z*dabc1_x__dcy
      daux_y__dcz=
     2 bc1_x*dabc1_z__dcz + abc1_z*dbc1_x__dcz
     4 - abc1_x*dbc1_z__dcz - bc1_z*dabc1_x__dcz
c      daux_y__ddx=0
c      daux_y__ddy=0
c      daux_y__ddz=0
c aux_z=abc1_x*bc1_y-bc1_x*abc1_y
      daux_z__dax=
     2 bc1_y*dabc1_x__dax - bc1_x*dabc1_y__dax
      daux_z__day=
     2 bc1_y*dabc1_x__day - bc1_x*dabc1_y__day
      daux_z__daz=
     2 bc1_y*dabc1_x__daz - bc1_x*dabc1_y__daz
      daux_z__dbx=
     2 bc1_y*dabc1_x__dbx + abc1_x*dbc1_y__dbx
     4 - abc1_y*dbc1_x__dbx - bc1_x*dabc1_y__dbx
      daux_z__dby=
     2 bc1_y*dabc1_x__dby + abc1_x*dbc1_y__dby
     4 - abc1_y*dbc1_x__dby - bc1_x*dabc1_y__dby
      daux_z__dbz=
     2 bc1_y*dabc1_x__dbz + abc1_x*dbc1_y__dbz
     4 - abc1_y*dbc1_x__dbz - bc1_x*dabc1_y__dbz
      daux_z__dcx=
     2 bc1_y*dabc1_x__dcx + abc1_x*dbc1_y__dcx
     4 - abc1_y*dbc1_x__dcx - bc1_x*dabc1_y__dcx
      daux_z__dcy=
     2 bc1_y*dabc1_x__dcy + abc1_x*dbc1_y__dcy
     4 - abc1_y*dbc1_x__dcy - bc1_x*dabc1_y__dcy
      daux_z__dcz=
     2 bc1_y*dabc1_x__dcz + abc1_x*dbc1_y__dcz
     4 - abc1_y*dbc1_x__dcz - bc1_x*dabc1_y__dcz
c      daux_z__ddx=0
c      daux_z__ddy=0
c      daux_z__ddz=0

c y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
c x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
c      dy__daux_x=bcd1_x
c      dy__daux_y=bcd1_y
c      dy__daux_z=bcd1_z
c
c      dy__dbcd1_x=aux_x
c      dy__dbcd1_y=aux_y
c      dy__dbcd1_z=aux_z
c
c      dx__dabc1_x=bcd1_x
c      dx__dabc1_y=bcd1_y
c      dx__dabc1_z=bcd1_z
c
c      dx__dbcd1_x=abc1_x
c      dx__dbcd1_y=abc1_y
c      dx__dbcd1_z=abc1_z

c derivation of y
c y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
      dy__dax=
     2 bcd1_x*daux_x__dax + 
     3 bcd1_y*daux_y__dax + 
     4 bcd1_z*daux_z__dax  
      dy__day=
     2 bcd1_x*daux_x__day + 
     3 bcd1_y*daux_y__day + 
     4 bcd1_z*daux_z__day 
      dy__daz=
     2 bcd1_x*daux_x__daz + 
     3 bcd1_y*daux_y__daz + 
     4 bcd1_z*daux_z__daz 
      dy__dbx=
     2 bcd1_x*daux_x__dbx + aux_x*dbcd1_x__dbx +
     3 bcd1_y*daux_y__dbx + aux_y*dbcd1_y__dbx +
     4 bcd1_z*daux_z__dbx + aux_z*dbcd1_z__dbx
      dy__dby=
     2 bcd1_x*daux_x__dby + aux_x*dbcd1_x__dby +
     3 bcd1_y*daux_y__dby + aux_y*dbcd1_y__dby +
     4 bcd1_z*daux_z__dby + aux_z*dbcd1_z__dby
      dy__dbz=
     2 bcd1_x*daux_x__dbz + aux_x*dbcd1_x__dbz +
     3 bcd1_y*daux_y__dbz + aux_y*dbcd1_y__dbz +
     4 bcd1_z*daux_z__dbz + aux_z*dbcd1_z__dbz
      dy__dcx=
     2 bcd1_x*daux_x__dcx + aux_x*dbcd1_x__dcx +
     3 bcd1_y*daux_y__dcx + aux_y*dbcd1_y__dcx +
     4 bcd1_z*daux_z__dcx + aux_z*dbcd1_z__dcx
      dy__dcy=
     2 bcd1_x*daux_x__dcy + aux_x*dbcd1_x__dcy +
     3 bcd1_y*daux_y__dcy + aux_y*dbcd1_y__dcy +
     4 bcd1_z*daux_z__dcy + aux_z*dbcd1_z__dcy
      dy__dcz=
     2 bcd1_x*daux_x__dcz + aux_x*dbcd1_x__dcz +
     3 bcd1_y*daux_y__dcz + aux_y*dbcd1_y__dcz +
     4 bcd1_z*daux_z__dcz + aux_z*dbcd1_z__dcz
      dy__ddx=
     2 aux_x*dbcd1_x__ddx +
     3 aux_y*dbcd1_y__ddx +
     4 aux_z*dbcd1_z__ddx
      dy__ddy=
     2 aux_x*dbcd1_x__ddy +
     3 aux_y*dbcd1_y__ddy +
     4 aux_z*dbcd1_z__ddy
      dy__ddz=
     2 aux_x*dbcd1_x__ddz +
     3 aux_y*dbcd1_y__ddz +
     4 aux_z*dbcd1_z__ddz

c derivation of x
c x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
      dx__dax=
     2 bcd1_x*dabc1_x__dax + 
     3 bcd1_y*dabc1_y__dax + 
     4 bcd1_z*dabc1_z__dax  
      dx__day=
     2 bcd1_x*dabc1_x__day +
     3 bcd1_y*dabc1_y__day + 
     4 bcd1_z*dabc1_z__day 
      dx__daz=
     2 bcd1_x*dabc1_x__daz + 
     3 bcd1_y*dabc1_y__daz + 
     4 bcd1_z*dabc1_z__daz 
      dx__dbx=
     2 bcd1_x*dabc1_x__dbx + abc1_x*dbcd1_x__dbx +
     3 bcd1_y*dabc1_y__dbx + abc1_y*dbcd1_y__dbx +
     4 bcd1_z*dabc1_z__dbx + abc1_z*dbcd1_z__dbx
      dx__dby=
     2 bcd1_x*dabc1_x__dby + abc1_x*dbcd1_x__dby +
     3 bcd1_y*dabc1_y__dby + abc1_y*dbcd1_y__dby +
     4 bcd1_z*dabc1_z__dby + abc1_z*dbcd1_z__dby
      dx__dbz=
     2 bcd1_x*dabc1_x__dbz + abc1_x*dbcd1_x__dbz +
     3 bcd1_y*dabc1_y__dbz + abc1_y*dbcd1_y__dbz +
     4 bcd1_z*dabc1_z__dbz + abc1_z*dbcd1_z__dbz
      dx__dcx=
     2 bcd1_x*dabc1_x__dcx + abc1_x*dbcd1_x__dcx +
     3 bcd1_y*dabc1_y__dcx + abc1_y*dbcd1_y__dcx +
     4 bcd1_z*dabc1_z__dcx + abc1_z*dbcd1_z__dcx
      dx__dcy=
     2 bcd1_x*dabc1_x__dcy + abc1_x*dbcd1_x__dcy +
     3 bcd1_y*dabc1_y__dcy + abc1_y*dbcd1_y__dcy +
     4 bcd1_z*dabc1_z__dcy + abc1_z*dbcd1_z__dcy
      dx__dcz=
     2 bcd1_x*dabc1_x__dcz + abc1_x*dbcd1_x__dcz +
     3 bcd1_y*dabc1_y__dcz + abc1_y*dbcd1_y__dcz +
     4 bcd1_z*dabc1_z__dcz + abc1_z*dbcd1_z__dcz
      dx__ddx=
     2 abc1_x*dbcd1_x__ddx +
     3 abc1_y*dbcd1_y__ddx +
     4 abc1_z*dbcd1_z__ddx
      dx__ddy=
     2 abc1_x*dbcd1_x__ddy +
     3 abc1_y*dbcd1_y__ddy +
     4 abc1_z*dbcd1_z__ddy
      dx__ddz=
     2 abc1_x*dbcd1_x__ddz +
     3 abc1_y*dbcd1_y__ddz +
     4 abc1_z*dbcd1_z__ddz

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

c f is the function atan2(y, x)
c fortran (and most other sources) use 'atan2(y,x)' while mathematica uses 'atan2(x,y)'

c////////////////////////////////////////////////////////////////////////
c////////////////////////////////////////////////////////////////////////
c// SECOND DERIVATIVES
c////////////////////////////////////////////////////////////////////////
c///////////////////////////////////////////////////////////////////////




      ddf11dx__dax
      ddf11dx__day
      ddf11dx__daz
      ddf11dx__dbx
      ddf11dx__dby
      ddf11dx__dbz
      ddf11dx__dcx
      ddf11dx__dcy
      ddf11dx__dcz
      ddf11dx__ddx
      ddf11dx__ddy
      ddf11dx__ddz


      ddf11dy__dax
      ddf11dy__day
      ddf11dy__daz
      ddf11dy__dbx
      ddf11dy__dby
      ddf11dy__dbz
      ddf11dy__dcx
      ddf11dy__dcy
      ddf11dy__dcz
      ddf11dy__ddx
      ddf11dy__ddy
      ddf11dy__ddz


      ddx11dax__dax
      ddx11dax__day
      ddx11dax__daz
      ddx11dax__dbx
      ddx11dax__dby
      ddx11dax__dbz
      ddx11dax__dcx
      ddx11dax__dcy
      ddx11dax__dcz
      ddx11dax__ddx
      ddx11dax__ddy
      ddx11dax__ddz
      ddx11day__day
      ddx11day__daz
      ddx11day__dbx
      ddx11day__dby
      ddx11day__dbz
      ddx11day__dcx
      ddx11day__dcy
      ddx11day__dcz
      ddx11day__ddx
      ddx11day__ddy
      ddx11day__ddz
      ddx11daz__daz
      ddx11daz__dbx
      ddx11daz__dby
      ddx11daz__dbz
      ddx11daz__dcx
      ddx11daz__dcy
      ddx11daz__dcz
      ddx11daz__ddx
      ddx11daz__ddy
      ddx11daz__ddz
      ddx11dbx__dbx
      ddx11dbx__dby
      ddx11dbx__dbz
      ddx11dbx__dcx
      ddx11dbx__dcy
      ddx11dbx__dcz
      ddx11dbx__ddx
      ddx11dbx__ddy
      ddx11dbx__ddz
      ddx11dby__dby
      ddx11dby__dbz
      ddx11dby__dcx
      ddx11dby__dcy
      ddx11dby__dcz
      ddx11dby__ddx
      ddx11dby__ddy
      ddx11dby__ddz
      ddx11dbz__dbz
      ddx11dbz__dcx
      ddx11dbz__dcy
      ddx11dbz__dcz
      ddx11dbz__ddx
      ddx11dbz__ddy
      ddx11dbz__ddz
      ddx11dcx__dcx
      ddx11dcx__dcy
      ddx11dcx__dcz
      ddx11dcx__ddx
      ddx11dcx__ddy
      ddx11dcx__ddz
      ddx11dcy__dcy
      ddx11dcy__dcz
      ddx11dcy__ddx
      ddx11dcy__ddy
      ddx11dcy__ddz
      ddx11dcz__dcz
      ddx11dcz__ddx
      ddx11dcz__ddy
      ddx11dcz__ddz
      ddx11ddx__ddx
      ddx11ddx__ddy
      ddx11ddx__ddz
      ddx11ddy__ddy
      ddx11ddy__ddz
      ddx11ddz__ddz


c dy__dax= bcd1_x*daux_x__dax + bcd1_y*daux_y__dax + bcd1_z*daux_z__dax
      ddy11dax__dax=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dax
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dax
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dax
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dax
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dax
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dax
      ddy11dax__day=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__day
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__day
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__day
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__day
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__day
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__day
      ddy11dax__daz=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__daz
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__daz
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__daz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__daz
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__daz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__daz
      ddy11dax__dbx=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dbx
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dbx
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dbx
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dbx
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dbx
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dbx
      ddy11dax__dby=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dby
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dby
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dby
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dby
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dby
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dby
      ddy11dax__dbz=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dbz
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dbz
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dbz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dbz
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dbz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dbz
      ddy11dax__dcx=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dcx
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dcx
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dcx
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dcx
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dcx
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dcx
      ddy11dax__dcy=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dcy
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dcy
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dcy
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dcy
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dcy
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dcy
      ddy11dax__dcz=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__dcz
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dcz
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__dcz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dcz
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__dcz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dcz
      ddy11dax__ddx=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__ddx
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__ddx
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__ddx
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__ddx
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__ddx
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__ddx
      ddy11dax__ddy=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__ddy
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__ddy
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__ddy
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__ddy
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__ddy
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__ddy
      ddy11dax__ddz=
     1              ddy11dax__ddy11daux_x*ddy11daux_x__ddz
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__ddz
     1            + ddy11dax__ddy11daux_y*ddy11daux_y__ddz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__ddz
     1            + ddy11dax__ddy11daux_z*ddy11daux_z__ddz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__ddz
c dy__day= bcd1_x*daux_x__day + bcd1_y*daux_y__day + bcd1_z*daux_z__day
      ddy11day__day=
     1              ddy11day__ddy11daux_x*ddy11daux_x__day
     1            + ddy11day__ddaux_x11day*ddaux_x11day__day
     1            + ddy11day__ddy11daux_y*ddy11daux_y__day
     1            + ddy11day__ddaux_y11day*ddaux_y11day__day
     1            + ddy11day__ddy11daux_z*ddy11daux_z__day
     1            + ddy11day__ddaux_z11day*ddaux_z11day__day
      ddy11day__daz=
     1              ddy11day__ddy11daux_x*ddy11daux_x__daz
     1            + ddy11day__ddaux_x11day*ddaux_x11day__daz
     1            + ddy11day__ddy11daux_y*ddy11daux_y__daz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__daz
     1            + ddy11day__ddy11daux_z*ddy11daux_z__daz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__daz
      ddy11day__dbx=
     1              ddy11day__ddy11daux_x*ddy11daux_x__dbx
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dbx
     1            + ddy11day__ddy11daux_y*ddy11daux_y__dbx
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dbx
     1            + ddy11day__ddy11daux_z*ddy11daux_z__dbx
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dbx
      ddy11day__dby=
     1              ddy11day__ddy11daux_x*ddy11daux_x__dby
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dby
     1            + ddy11day__ddy11daux_y*ddy11daux_y__dby
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dby
     1            + ddy11day__ddy11daux_z*ddy11daux_z__dby
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dby
      ddy11day__dbz=
     1              ddy11day__ddy11daux_x*ddy11daux_x__dbz
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dbz
     1            + ddy11day__ddy11daux_y*ddy11daux_y__dbz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dbz
     1            + ddy11day__ddy11daux_z*ddy11daux_z__dbz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dbz
      ddy11day__dcx=
     1              ddy11day__ddy11daux_x*ddy11daux_x__dcx
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dcx
     1            + ddy11day__ddy11daux_y*ddy11daux_y__dcx
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dcx
     1            + ddy11day__ddy11daux_z*ddy11daux_z__dcx
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dcx
      ddy11day__dcy=
     1              ddy11day__ddy11daux_x*ddy11daux_x__dcy
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dcy
     1            + ddy11day__ddy11daux_y*ddy11daux_y__dcy
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dcy
     1            + ddy11day__ddy11daux_z*ddy11daux_z__dcy
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dcy
      ddy11day__dcz=
     1              ddy11day__ddy11daux_x*ddy11daux_x__dcz
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dcz
     1            + ddy11day__ddy11daux_y*ddy11daux_y__dcz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dcz
     1            + ddy11day__ddy11daux_z*ddy11daux_z__dcz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dcz
      ddy11day__ddx=
     1              ddy11day__ddy11daux_x*ddy11daux_x__ddx
     1            + ddy11day__ddaux_x11day*ddaux_x11day__ddx
     1            + ddy11day__ddy11daux_y*ddy11daux_y__ddx
     1            + ddy11day__ddaux_y11day*ddaux_y11day__ddx
     1            + ddy11day__ddy11daux_z*ddy11daux_z__ddx
     1            + ddy11day__ddaux_z11day*ddaux_z11day__ddx
      ddy11day__ddy=
     1              ddy11day__ddy11daux_x*ddy11daux_x__ddy
     1            + ddy11day__ddaux_x11day*ddaux_x11day__ddy
     1            + ddy11day__ddy11daux_y*ddy11daux_y__ddy
     1            + ddy11day__ddaux_y11day*ddaux_y11day__ddy
     1            + ddy11day__ddy11daux_z*ddy11daux_z__ddy
     1            + ddy11day__ddaux_z11day*ddaux_z11day__ddy
      ddy11day__ddz=
     1              ddy11day__ddy11daux_x*ddy11daux_x__ddz
     1            + ddy11day__ddaux_x11day*ddaux_x11day__ddz
     1            + ddy11day__ddy11daux_y*ddy11daux_y__ddz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__ddz
     1            + ddy11day__ddy11daux_z*ddy11daux_z__ddz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__ddz
c dy__daz= bcd1_x*daux_x__daz + bcd1_y*daux_y__daz + bcd1_z*daux_z__daz
      ddy11daz__daz=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__daz
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__daz
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__daz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__daz
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__daz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__daz
      ddy11daz__dbx=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dbx
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dbx
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dbx
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dbx
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dbx
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dbx
      ddy11daz__dby=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dby
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dby
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dby
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dby
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dby
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dby
      ddy11daz__dbz=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dbz
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dbz
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dbz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dbz
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dbz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dbz
      ddy11daz__dcx=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dcx
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dcx
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dcx
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dcx
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dcx
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dcx
      ddy11daz__dcy=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dcy
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dcy
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dcy
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dcy
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dcy
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dcy
      ddy11daz__dcz=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dcz
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dcz
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dcz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dcz
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dcz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dcz
      ddy11daz__ddx=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__ddx
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__ddx
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__ddx
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__ddx
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__ddx
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__ddx
      ddy11daz__ddy=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__ddy
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__ddy
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__ddy
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__ddy
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__ddy
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__ddy
      ddy11daz__ddz=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__ddz
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__ddz
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__ddz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__ddz
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__ddz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__ddz
c dy__dbx= bcd1_x*daux_x__dbx + aux_x*dbcd1_x__dbx + bcd1_y*daux_y__dbx + aux_y*dbcd1_y__dbx + bcd1_z*daux_z__dbx + aux_z*dbcd1_z__dbx
      ddy11dbx__dbx=
     1              ddy11daz__ddy11daux_x*ddy11daux_x__dbx
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dbx
     1            + ddy11daz__ddy11daux_y*ddy11daux_y__dbx
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dbx
     1            + ddy11daz__ddy11daux_z*ddy11daux_z__dbx
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dbx
      ddy11dbx__dby=
      ddy11dbx__dbz=
      ddy11dbx__dcx=
      ddy11dbx__dcy=
      ddy11dbx__dcz=
      ddy11dbx__ddx=
      ddy11dbx__ddy=
      ddy11dbx__ddz=
c dy__dby= bcd1_x*daux_x__dby + aux_x*dbcd1_x__dby + bcd1_y*daux_y__dby + aux_y*dbcd1_y__dby + bcd1_z*daux_z__dby + aux_z*dbcd1_z__dby
c dy__dbz= bcd1_x*daux_x__dbz + aux_x*dbcd1_x__dbz + bcd1_y*daux_y__dbz + aux_y*dbcd1_y__dbz + bcd1_z*daux_z__dbz + aux_z*dbcd1_z__dbz
c dy__dcx= bcd1_x*daux_x__dcx + aux_x*dbcd1_x__dcx + bcd1_y*daux_y__dcx + aux_y*dbcd1_y__dcx + bcd1_z*daux_z__dcx + aux_z*dbcd1_z__dcx
c dy__dcy= bcd1_x*daux_x__dcy + aux_x*dbcd1_x__dcy + bcd1_y*daux_y__dcy + aux_y*dbcd1_y__dcy + bcd1_z*daux_z__dcy + aux_z*dbcd1_z__dcy
c dy__dcz= bcd1_x*daux_x__dcz + aux_x*dbcd1_x__dcz + bcd1_y*daux_y__dcz + aux_y*dbcd1_y__dcz + bcd1_z*daux_z__dcz + aux_z*dbcd1_z__dcz
c dy__ddx= aux_x*dbcd1_x__ddx + aux_y*dbcd1_y__ddx + aux_z*dbcd1_z__ddx
c dy__ddy= aux_x*dbcd1_x__ddy + aux_y*dbcd1_y__ddy + aux_z*dbcd1_z__ddy
c dy__ddz= aux_x*dbcd1_x__ddz + aux_y*dbcd1_y__ddz + aux_z*dbcd1_z__ddz
      ddy11dby__dby=
      ddy11dby__dbz=
      ddy11dby__dcx=
      ddy11dby__dcy=
      ddy11dby__dcz=
      ddy11dby__ddx=
      ddy11dby__ddy=
      ddy11dby__ddz=
      ddy11dbz__dbz=
      ddy11dbz__dcx=
      ddy11dbz__dcy=
      ddy11dbz__dcz=
      ddy11dbz__ddx=
      ddy11dbz__ddy=
      ddy11dbz__ddz=
      ddy11dcx__dcx=
      ddy11dcx__dcy=
      ddy11dcx__dcz=
      ddy11dcx__ddx=
      ddy11dcx__ddy=
      ddy11dcx__ddz=
      ddy11dcy__dcy=
      ddy11dcy__dcz=
      ddy11dcy__ddx=
      ddy11dcy__ddy=
      ddy11dcy__ddz=
      ddy11dcz__dcz=
      ddy11dcz__ddx=
      ddy11dcz__ddy=
      ddy11dcz__ddz=
      ddy11ddx__ddx=
      ddy11ddx__ddy=
      ddy11ddx__ddz=
      ddy11ddy__ddy=
      ddy11ddy__ddz=
      ddy11ddz__ddz=


c df__dax=df__dx*dx__dax + df__dy*dy__dax
      ddf11dax__ddf11dx =dx__dax
      ddf11dax__ddf11dy =dy__dax
      ddf11dax__ddx11dax=df__dx
      ddf11dax__ddy11dax=df__dy
c df__day=df__dx*dx__day + df__dy*dy__day
      ddf11day__ddf11dx =dx__day
      ddf11day__ddf11dy =dy__day
      ddf11day__ddx11day=df__dx
      ddf11day__ddy11day=df__dy
c df__daz=df__dx*dx__daz + df__dy*dy__daz
      ddf11daz__ddf11dx =dx__daz
      ddf11daz__ddf11dy =dy__daz
      ddf11daz__ddx11daz=df__dx
      ddf11daz__ddy11daz=df__dy
c df__dbx=df__dx*dx__dbx + df__dy*dy__dbx
      ddf11dbx__ddf11dx =dx__dbx
      ddf11dbx__ddf11dy =dy__dbx
      ddf11dbx__ddx11dbx=df__dx
      ddf11dbx__ddy11dbx=df__dy
c df__dby=df__dx*dx__dby + df__dy*dy__dby
      ddf11dby__ddf11dx =dx__dby
      ddf11dby__ddf11dy =dy__dby
      ddf11dby__ddx11dby=df__dx
      ddf11dby__ddy11dby=df__dy
c df__dbz=df__dx*dx__dbz + df__dy*dy__dbz
      ddf11dbz__ddf11dx =dx__dbz
      ddf11dbz__ddf11dy =dy__dbz
      ddf11dbz__ddx11dbz=df__dx
      ddf11dbz__ddy11dbz=df__dy
c df__dcx=df__dx*dx__dcx + df__dy*dy__dcx
      ddf11dcx__ddf11dx =dx__dcx
      ddf11dcx__ddf11dy =dy__dcx
      ddf11dcx__ddx11dcx=df__dx
      ddf11dcx__ddy11dcx=df__dy
c df__dcy=df__dx*dx__dcy + df__dy*dy__dcy
      ddf11dcy__ddf11dx =dx__dcy
      ddf11dcy__ddf11dy =dy__dcy
      ddf11dcy__ddx11dcy=df__dx
      ddf11dcy__ddy11dcy=df__dy
c df__dcz=df__dx*dx__dcz + df__dy*dy__dcz
      ddf11dcz__ddf11dx =dx__dcz
      ddf11dcz__ddf11dy =dy__dcz
      ddf11dcz__ddx11dcz=df__dx
      ddf11dcz__ddy11dcz=df__dy
c df__ddx=df__dx*dx__ddx + df__dy*dy__ddx
      ddf11ddx__ddf11dx =dx__ddx
      ddf11ddx__ddf11dy =dy__ddx
      ddf11ddx__ddx11ddx=df__dx
      ddf11ddx__ddy11ddx=df__dy
c df__ddy=df__dx*dx__ddy + df__dy*dy__ddy
      ddf11ddy__ddf11dx =dx__ddy
      ddf11ddy__ddf11dy =dy__ddy
      ddf11ddy__ddx11ddy=df__dx
      ddf11ddy__ddy11ddy=df__dy
c df__ddz=df__dx*dx__ddz + df__dy*dy__ddz
      ddf11ddz__ddf11dx =dx__ddz
      ddf11ddz__ddf11dy =dy__ddz
      ddf11ddz__ddx11ddz=df__dx
      ddf11ddz__ddy11ddz=df__dy


c df__dax=df__dx*dx__dax + df__dy*dy__dax
      ddf11dax__dax=
     1   ddf11dax__ddf11dx*ddf11dx__dax
     1 + ddf11dax__ddx11dax*ddx11dax__dax
     1 + ddf11dax__ddf11dy*ddf11dy__dax
     1 + ddf11dax__ddy11dax*ddy11dax__dax
      ddf11dax__day=
     1   ddf11dax__ddf11dx*ddf11dx__day
     1 + ddf11dax__ddx11dax*ddx11dax__day
     1 + ddf11dax__ddf11dy*ddf11dy__day
     1 + ddf11dax__ddy11dax*ddy11dax__day
      ddf11dax__daz=
     1   ddf11dax__ddf11dx*ddf11dx__daz
     1 + ddf11dax__ddx11dax*ddx11dax__daz
     1 + ddf11dax__ddf11dy*ddf11dy__daz
     1 + ddf11dax__ddy11dax*ddy11dax__daz
      ddf11dax__dbx=
     1   ddf11dax__ddf11dx*ddf11dx__dbx
     1 + ddf11dax__ddx11dax*ddx11dax__dbx
     1 + ddf11dax__ddf11dy*ddf11dy__dbx
     1 + ddf11dax__ddy11dax*ddy11dax__dbx
      ddf11dax__dby=
     1   ddf11dax__ddf11dx*ddf11dx__dby
     1 + ddf11dax__ddx11dax*ddx11dax__dby
     1 + ddf11dax__ddf11dy*ddf11dy__dby
     1 + ddf11dax__ddy11dax*ddy11dax__dby
      ddf11dax__dbz=
     1   ddf11dax__ddf11dx*ddf11dx__dbz
     1 + ddf11dax__ddx11dax*ddx11dax__dbz
     1 + ddf11dax__ddf11dy*ddf11dy__dbz
     1 + ddf11dax__ddy11dax*ddy11dax__dbz
      ddf11dax__dcx=
     1   ddf11dax__ddf11dx*ddf11dx__dcx
     1 + ddf11dax__ddx11dax*ddx11dax__dcx
     1 + ddf11dax__ddf11dy*ddf11dy__dcx
     1 + ddf11dax__ddy11dax*ddy11dax__dcx
      ddf11dax__dcy=
     1   ddf11dax__ddf11dx*ddf11dx__dcy
     1 + ddf11dax__ddx11dax*ddx11dax__dcy
     1 + ddf11dax__ddf11dy*ddf11dy__dcy
     1 + ddf11dax__ddy11dax*ddy11dax__dcy
      ddf11dax__dcz=
     1   ddf11dax__ddf11dx*ddf11dx__dcz
     1 + ddf11dax__ddx11dax*ddx11dax__dcz
     1 + ddf11dax__ddf11dy*ddf11dy__dcz
     1 + ddf11dax__ddy11dax*ddy11dax__dcz
      ddf11dax__ddx=
     1   ddf11dax__ddf11dx*ddf11dx__ddx
     1 + ddf11dax__ddx11dax*ddx11dax__ddx
     1 + ddf11dax__ddf11dy*ddf11dy__ddx
     1 + ddf11dax__ddy11dax*ddy11dax__ddx
      ddf11dax__ddy=
     1   ddf11dax__ddf11dx*ddf11dx__ddy
     1 + ddf11dax__ddx11dax*ddx11dax__ddy
     1 + ddf11dax__ddf11dy*ddf11dy__ddy
     1 + ddf11dax__ddy11dax*ddy11dax__ddy
      ddf11dax__ddz=
     1   ddf11dax__ddf11dx*ddf11dx__ddz
     1 + ddf11dax__ddx11dax*ddx11dax__ddz
     1 + ddf11dax__ddf11dy*ddf11dy__ddz
     1 + ddf11dax__ddy11dax*ddy11dax__ddz
c df__day=df__dx*dx__day + df__dy*dy__day
      ddf11day__day=
     1   ddf11day__ddf11dx*ddf11dx__day
     1 + ddf11day__ddx11day*ddx11day__day
     1 + ddf11day__ddf11dy*ddf11dy__day
     1 + ddf11day__ddy11day*ddy11day__day
      ddf11day__daz=
     1   ddf11day__ddf11dx*ddf11dx__daz
     1 + ddf11day__ddx11day*ddx11day__daz
     1 + ddf11day__ddf11dy*ddf11dy__daz
     1 + ddf11day__ddy11day*ddy11day__daz
      ddf11day__dbx=
     1   ddf11day__ddf11dx*ddf11dx__dbx
     1 + ddf11day__ddx11day*ddx11day__dbx
     1 + ddf11day__ddf11dy*ddf11dy__dbx
     1 + ddf11day__ddy11day*ddy11day__dbx
      ddf11day__dby=
     1   ddf11day__ddf11dx*ddf11dx__dby
     1 + ddf11day__ddx11day*ddx11day__dby
     1 + ddf11day__ddf11dy*ddf11dy__dby
     1 + ddf11day__ddy11day*ddy11day__dby
      ddf11day__dbz=
     1   ddf11day__ddf11dx*ddf11dx__dbz
     1 + ddf11day__ddx11day*ddx11day__dbz
     1 + ddf11day__ddf11dy*ddf11dy__dbz
     1 + ddf11day__ddy11day*ddy11day__dbz
      ddf11day__dcx=
     1   ddf11day__ddf11dx*ddf11dx__dcx
     1 + ddf11day__ddx11day*ddx11day__dcx
     1 + ddf11day__ddf11dy*ddf11dy__dcx
     1 + ddf11day__ddy11day*ddy11day__dcx
      ddf11day__dcy=
     1   ddf11day__ddf11dx*ddf11dx__dcy
     1 + ddf11day__ddx11day*ddx11day__dcy
     1 + ddf11day__ddf11dy*ddf11dy__dcy
     1 + ddf11day__ddy11day*ddy11day__dcy
      ddf11day__dcz=
     1   ddf11day__ddf11dx*ddf11dx__dcz
     1 + ddf11day__ddx11day*ddx11day__dcz
     1 + ddf11day__ddf11dy*ddf11dy__dcz
     1 + ddf11day__ddy11day*ddy11day__dcz
      ddf11day__ddx=
     1   ddf11day__ddf11dx*ddf11dx__ddx
     1 + ddf11day__ddx11day*ddx11day__ddx
     1 + ddf11day__ddf11dy*ddf11dy__ddx
     1 + ddf11day__ddy11day*ddy11day__ddx
      ddf11day__ddy=
     1   ddf11day__ddf11dx*ddf11dx__ddy
     1 + ddf11day__ddx11day*ddx11day__ddy
     1 + ddf11day__ddf11dy*ddf11dy__ddy
     1 + ddf11day__ddy11day*ddy11day__ddy
      ddf11day__ddz=
     1   ddf11day__ddf11dx*ddf11dx__ddz
     1 + ddf11day__ddx11day*ddx11day__ddz
     1 + ddf11day__ddf11dy*ddf11dy__ddz
     1 + ddf11day__ddy11day*ddy11day__ddz
c df__daz=df__dx*dx__daz + df__dy*dy__daz
      ddf11daz__daz=
     1   ddf11daz__ddf11dx*ddf11dx__daz
     1 + ddf11daz__ddx11daz*ddx11daz__daz
     1 + ddf11daz__ddf11dy*ddf11dy__daz
     1 + ddf11daz__ddy11daz*ddy11daz__daz
      ddf11daz__dbx=
     1   ddf11daz__ddf11dx*ddf11dx__dbx
     1 + ddf11daz__ddx11daz*ddx11daz__dbx
     1 + ddf11daz__ddf11dy*ddf11dy__dbx
     1 + ddf11daz__ddy11daz*ddy11daz__dbx
      ddf11daz__dby=
     1   ddf11daz__ddf11dx*ddf11dx__dby
     1 + ddf11daz__ddx11daz*ddx11daz__dby
     1 + ddf11daz__ddf11dy*ddf11dy__dby
     1 + ddf11daz__ddy11daz*ddy11daz__dby
      ddf11daz__dbz=
     1   ddf11daz__ddf11dx*ddf11dx__dbz
     1 + ddf11daz__ddx11daz*ddx11daz__dbz
     1 + ddf11daz__ddf11dy*ddf11dy__dbz
     1 + ddf11daz__ddy11daz*ddy11daz__dbz
      ddf11daz__dcx=
     1   ddf11daz__ddf11dx*ddf11dx__dcx
     1 + ddf11daz__ddx11daz*ddx11daz__dcx
     1 + ddf11daz__ddf11dy*ddf11dy__dcx
     1 + ddf11daz__ddy11daz*ddy11daz__dcx
      ddf11daz__dcy=
     1   ddf11daz__ddf11dx*ddf11dx__dcy
     1 + ddf11daz__ddx11daz*ddx11daz__dcy
     1 + ddf11daz__ddf11dy*ddf11dy__dcy
     1 + ddf11daz__ddy11daz*ddy11daz__dcy
      ddf11daz__dcz=
     1   ddf11daz__ddf11dx*ddf11dx__dcz
     1 + ddf11daz__ddx11daz*ddx11daz__dcz
     1 + ddf11daz__ddf11dy*ddf11dy__dcz
     1 + ddf11daz__ddy11daz*ddy11daz__dcz
      ddf11daz__ddx=
     1   ddf11daz__ddf11dx*ddf11dx__ddx
     1 + ddf11daz__ddx11daz*ddx11daz__ddx
     1 + ddf11daz__ddf11dy*ddf11dy__ddx
     1 + ddf11daz__ddy11daz*ddy11daz__ddx
      ddf11daz__ddy=
     1   ddf11daz__ddf11dx*ddf11dx__ddy
     1 + ddf11daz__ddx11daz*ddx11daz__ddy
     1 + ddf11daz__ddf11dy*ddf11dy__ddy
     1 + ddf11daz__ddy11daz*ddy11daz__ddy
      ddf11daz__ddz=
     1   ddf11daz__ddf11dx*ddf11dx__ddz
     1 + ddf11daz__ddx11daz*ddx11daz__ddz
     1 + ddf11daz__ddf11dy*ddf11dy__ddz
     1 + ddf11daz__ddy11daz*ddy11daz__ddz
c df__dbx=df__dx*dx__dbx + df__dy*dy__dbx
      ddf11dbx__dbx=
     1   ddf11dbx__ddf11dx*ddf11dx__dbx
     1 + ddf11dbx__ddx11dbx*ddx11dbx__dbx
     1 + ddf11dbx__ddf11dy*ddf11dy__dbx
     1 + ddf11dbx__ddy11dbx*ddy11dbx__dbx
      ddf11dbx__dby=
     1   ddf11dbx__ddf11dx*ddf11dx__dby
     1 + ddf11dbx__ddx11dbx*ddx11dbx__dby
     1 + ddf11dbx__ddf11dy*ddf11dy__dby
     1 + ddf11dbx__ddy11dbx*ddy11dbx__dby
      ddf11dbx__dbz=
     1   ddf11dbx__ddf11dx*ddf11dx__dbz
     1 + ddf11dbx__ddx11dbx*ddx11dbx__dbz
     1 + ddf11dbx__ddf11dy*ddf11dy__dbz
     1 + ddf11dbx__ddy11dbx*ddy11dbx__dbz
      ddf11dbx__dcx=
     1   ddf11dbx__ddf11dx*ddf11dx__dcx
     1 + ddf11dbx__ddx11dbx*ddx11dbx__dcx
     1 + ddf11dbx__ddf11dy*ddf11dy__dcx
     1 + ddf11dbx__ddy11dbx*ddy11dbx__dcx
      ddf11dbx__dcy=
     1   ddf11dbx__ddf11dx*ddf11dx__dcy
     1 + ddf11dbx__ddx11dbx*ddx11dbx__dcy
     1 + ddf11dbx__ddf11dy*ddf11dy__dcy
     1 + ddf11dbx__ddy11dbx*ddy11dbx__dcy
      ddf11dbx__dcz=
     1   ddf11dbx__ddf11dx*ddf11dx__dcz
     1 + ddf11dbx__ddx11dbx*ddx11dbx__dcz
     1 + ddf11dbx__ddf11dy*ddf11dy__dcz
     1 + ddf11dbx__ddy11dbx*ddy11dbx__dcz
      ddf11dbx__ddx=
     1   ddf11dbx__ddf11dx*ddf11dx__ddx
     1 + ddf11dbx__ddx11dbx*ddx11dbx__ddx
     1 + ddf11dbx__ddf11dy*ddf11dy__ddx
     1 + ddf11dbx__ddy11dbx*ddy11dbx__ddx
      ddf11dbx__ddy=
     1   ddf11dbx__ddf11dx*ddf11dx__ddy
     1 + ddf11dbx__ddx11dbx*ddx11dbx__ddy
     1 + ddf11dbx__ddf11dy*ddf11dy__ddy
     1 + ddf11dbx__ddy11dbx*ddy11dbx__ddy
      ddf11dbx__ddz=
     1   ddf11dbx__ddf11dx*ddf11dx__ddz
     1 + ddf11dbx__ddx11dbx*ddx11dbx__ddz
     1 + ddf11dbx__ddf11dy*ddf11dy__ddz
     1 + ddf11dbx__ddy11dbx*ddy11dbx__ddz
c df__dby=df__dx*dx__dby + df__dy*dy__dby
      ddf11dby__dby=
     1   ddf11dby__ddf11dx*ddf11dx__dby
     1 + ddf11dby__ddx11dby*ddx11dby__dby
     1 + ddf11dby__ddf11dy*ddf11dy__dby
     1 + ddf11dby__ddy11dby*ddy11dby__dby
      ddf11dby__dbz=
     1   ddf11dby__ddf11dx*ddf11dx__dbz
     1 + ddf11dby__ddx11dby*ddx11dby__dbz
     1 + ddf11dby__ddf11dy*ddf11dy__dbz
     1 + ddf11dby__ddy11dby*ddy11dby__dbz
      ddf11dby__dcx=
     1   ddf11dby__ddf11dx*ddf11dx__dcx
     1 + ddf11dby__ddx11dby*ddx11dby__dcx
     1 + ddf11dby__ddf11dy*ddf11dy__dcx
     1 + ddf11dby__ddy11dby*ddy11dby__dcx
      ddf11dby__dcy=
     1   ddf11dby__ddf11dx*ddf11dx__dcy
     1 + ddf11dby__ddx11dby*ddx11dby__dcy
     1 + ddf11dby__ddf11dy*ddf11dy__dcy
     1 + ddf11dby__ddy11dby*ddy11dby__dcy
      ddf11dby__dcz=
     1   ddf11dby__ddf11dx*ddf11dx__dcz
     1 + ddf11dby__ddx11dby*ddx11dby__dcz
     1 + ddf11dby__ddf11dy*ddf11dy__dcz
     1 + ddf11dby__ddy11dby*ddy11dby__dcz
      ddf11dby__ddx=
     1   ddf11dby__ddf11dx*ddf11dx__ddx
     1 + ddf11dby__ddx11dby*ddx11dby__ddx
     1 + ddf11dby__ddf11dy*ddf11dy__ddx
     1 + ddf11dby__ddy11dby*ddy11dby__ddx
      ddf11dby__ddy=
     1   ddf11dby__ddf11dx*ddf11dx__ddy
     1 + ddf11dby__ddx11dby*ddx11dby__ddy
     1 + ddf11dby__ddf11dy*ddf11dy__ddy
     1 + ddf11dby__ddy11dby*ddy11dby__ddy
      ddf11dby__ddz=
     1   ddf11dby__ddf11dx*ddf11dx__ddz
     1 + ddf11dby__ddx11dby*ddx11dby__ddz
     1 + ddf11dby__ddf11dy*ddf11dy__ddz
     1 + ddf11dby__ddy11dby*ddy11dby__ddz
c df__dbz=df__dx*dx__dbz + df__dy*dy__dbz
      ddf11dbz__dbz=
     1   ddf11dbz__ddf11dx*ddf11dx__dbz
     1 + ddf11dbz__ddx11dbz*ddx11dbz__dbz
     1 + ddf11dbz__ddf11dy*ddf11dy__dbz
     1 + ddf11dbz__ddy11dbz*ddy11dbz__dbz
      ddf11dbz__dcx=
     1   ddf11dbz__ddf11dx*ddf11dx__dcx
     1 + ddf11dbz__ddx11dbz*ddx11dbz__dcx
     1 + ddf11dbz__ddf11dy*ddf11dy__dcx
     1 + ddf11dbz__ddy11dbz*ddy11dbz__dcx
      ddf11dbz__dcy=
     1   ddf11dbz__ddf11dx*ddf11dx__dcy
     1 + ddf11dbz__ddx11dbz*ddx11dbz__dcy
     1 + ddf11dbz__ddf11dy*ddf11dy__dcy
     1 + ddf11dbz__ddy11dbz*ddy11dbz__dcy
      ddf11dbz__dcz=
     1   ddf11dbz__ddf11dx*ddf11dx__dcz
     1 + ddf11dbz__ddx11dbz*ddx11dbz__dcz
     1 + ddf11dbz__ddf11dy*ddf11dy__dcz
     1 + ddf11dbz__ddy11dbz*ddy11dbz__dcz
      ddf11dbz__ddx=
     1   ddf11dbz__ddf11dx*ddf11dx__ddx
     1 + ddf11dbz__ddx11dbz*ddx11dbz__ddx
     1 + ddf11dbz__ddf11dy*ddf11dy__ddx
     1 + ddf11dbz__ddy11dbz*ddy11dbz__ddx
      ddf11dbz__ddy=
     1   ddf11dbz__ddf11dx*ddf11dx__ddy
     1 + ddf11dbz__ddx11dbz*ddx11dbz__ddy
     1 + ddf11dbz__ddf11dy*ddf11dy__ddy
     1 + ddf11dbz__ddy11dbz*ddy11dbz__ddy
      ddf11dbz__ddz=
     1   ddf11dbz__ddf11dx*ddf11dx__ddz
     1 + ddf11dbz__ddx11dbz*ddx11dbz__ddz
     1 + ddf11dbz__ddf11dy*ddf11dy__ddz
     1 + ddf11dbz__ddy11dbz*ddy11dbz__ddz
c df__dcx=df__dx*dx__dcx + df__dy*dy__dcx
      ddf11dcx__dcx=
     1   ddf11dcx__ddf11dx*ddf11dx__dcx
     1 + ddf11dcx__ddx11dcx*ddx11dcx__dcx
     1 + ddf11dcx__ddf11dy*ddf11dy__dcx
     1 + ddf11dcx__ddy11dcx*ddy11dcx__dcx
      ddf11dcx__dcy=
     1   ddf11dcx__ddf11dx*ddf11dx__dcy
     1 + ddf11dcx__ddx11dcx*ddx11dcx__dcy
     1 + ddf11dcx__ddf11dy*ddf11dy__dcy
     1 + ddf11dcx__ddy11dcx*ddy11dcx__dcy
      ddf11dcx__dcz=
     1   ddf11dcx__ddf11dx*ddf11dx__dcz
     1 + ddf11dcx__ddx11dcx*ddx11dcx__dcz
     1 + ddf11dcx__ddf11dy*ddf11dy__dcz
     1 + ddf11dcx__ddy11dcx*ddy11dcx__dcz
      ddf11dcx__ddx=
     1   ddf11dcx__ddf11dx*ddf11dx__ddx
     1 + ddf11dcx__ddx11dcx*ddx11dcx__ddx
     1 + ddf11dcx__ddf11dy*ddf11dy__ddx
     1 + ddf11dcx__ddy11dcx*ddy11dcx__ddx
      ddf11dcx__ddy=
     1   ddf11dcx__ddf11dx*ddf11dx__ddy
     1 + ddf11dcx__ddx11dcx*ddx11dcx__ddy
     1 + ddf11dcx__ddf11dy*ddf11dy__ddy
     1 + ddf11dcx__ddy11dcx*ddy11dcx__ddy
      ddf11dcx__ddz=
     1   ddf11dcx__ddf11dx*ddf11dx__ddz
     1 + ddf11dcx__ddx11dcx*ddx11dcx__ddz
     1 + ddf11dcx__ddf11dy*ddf11dy__ddz
     1 + ddf11dcx__ddy11dcx*ddy11dcx__ddz
c df__dcy=df__dx*dx__dcy + df__dy*dy__dcy
      ddf11dcy__dcy=
     1   ddf11dcy__ddf11dx*ddf11dx__dcy
     1 + ddf11dcy__ddx11dcy*ddx11dcy__dcy
     1 + ddf11dcy__ddf11dy*ddf11dy__dcy
     1 + ddf11dcy__ddy11dcy*ddy11dcy__dcy
      ddf11dcy__dcz=
     1   ddf11dcy__ddf11dx*ddf11dx__dcz
     1 + ddf11dcy__ddx11dcy*ddx11dcy__dcz
     1 + ddf11dcy__ddf11dy*ddf11dy__dcz
     1 + ddf11dcy__ddy11dcy*ddy11dcy__dcz
      ddf11dcy__ddx=
     1   ddf11dcy__ddf11dx*ddf11dx__ddx
     1 + ddf11dcy__ddx11dcy*ddx11dcy__ddx
     1 + ddf11dcy__ddf11dy*ddf11dy__ddx
     1 + ddf11dcy__ddy11dcy*ddy11dcy__ddx
      ddf11dcy__ddy=
     1   ddf11dcy__ddf11dx*ddf11dx__ddy
     1 + ddf11dcy__ddx11dcy*ddx11dcy__ddy
     1 + ddf11dcy__ddf11dy*ddf11dy__ddy
     1 + ddf11dcy__ddy11dcy*ddy11dcy__ddy
      ddf11dcy__ddz=
     1   ddf11dcy__ddf11dx*ddf11dx__ddz
     1 + ddf11dcy__ddx11dcy*ddx11dcy__ddz
     1 + ddf11dcy__ddf11dy*ddf11dy__ddz
     1 + ddf11dcy__ddy11dcy*ddy11dcy__ddz
c df__dcz=df__dx*dx__dcz + df__dy*dy__dcz
      ddf11dcz__dcz=
     1   ddf11dcz__ddf11dx*ddf11dx__dcz
     1 + ddf11dcz__ddx11dcz*ddx11dcz__dcz
     1 + ddf11dcz__ddf11dy*ddf11dy__dcz
     1 + ddf11dcz__ddy11dcz*ddy11dcz__dcz
      ddf11dcz__ddx=
     1   ddf11dcz__ddf11dx*ddf11dx__ddx
     1 + ddf11dcz__ddx11dcz*ddx11dcz__ddx
     1 + ddf11dcz__ddf11dy*ddf11dy__ddx
     1 + ddf11dcz__ddy11dcz*ddy11dcz__ddx
      ddf11dcz__ddy=
     1   ddf11dcz__ddf11dx*ddf11dx__ddy
     1 + ddf11dcz__ddx11dcz*ddx11dcz__ddy
     1 + ddf11dcz__ddf11dy*ddf11dy__ddy
     1 + ddf11dcz__ddy11dcz*ddy11dcz__ddy
      ddf11dcz__ddz=
     1   ddf11dcz__ddf11dx*ddf11dx__ddz
     1 + ddf11dcz__ddx11dcz*ddx11dcz__ddz
     1 + ddf11dcz__ddf11dy*ddf11dy__ddz
     1 + ddf11dcz__ddy11dcz*ddy11dcz__ddz
c df__ddx=df__dx*dx__ddx + df__dy*dy__ddx
      ddf11ddx__ddx=
     1   ddf11ddx__ddf11dx*ddf11dx__ddx
     1 + ddf11ddx__ddx11ddx*ddx11ddx__ddx
     1 + ddf11ddx__ddf11dy*ddf11dy__ddx
     1 + ddf11ddx__ddy11ddx*ddy11ddx__ddx
      ddf11ddx__ddy=
     1   ddf11ddx__ddf11dx*ddf11dx__ddy
     1 + ddf11ddx__ddx11ddx*ddx11ddx__ddy
     1 + ddf11ddx__ddf11dy*ddf11dy__ddy
     1 + ddf11ddx__ddy11ddx*ddy11ddx__ddy
      ddf11ddx__ddz=
     1   ddf11ddx__ddf11dx*ddf11dx__ddz
     1 + ddf11ddx__ddx11ddx*ddx11ddx__ddz
     1 + ddf11ddx__ddf11dy*ddf11dy__ddz
     1 + ddf11ddx__ddy11ddx*ddy11ddx__ddz
c df__ddy=df__dx*dx__ddy + df__dy*dy__ddy
      ddf11ddy__ddy=
     1   ddf11ddy__ddf11dx*ddf11dx__ddy
     1 + ddf11ddy__ddx11ddy*ddx11ddy__ddy
     1 + ddf11ddy__ddf11dy*ddf11dy__ddy
     1 + ddf11ddy__ddy11ddy*ddy11ddy__ddy
      ddf11ddy__ddz=
     1   ddf11ddy__ddf11dx*ddf11dx__ddz
     1 + ddf11ddy__ddx11ddy*ddx11ddy__ddz
     1 + ddf11ddy__ddf11dy*ddf11dy__ddz
     1 + ddf11ddy__ddy11ddy*ddy11ddy__ddz
c df__ddz=df__dx*dx__ddz + df__dy*dy__ddz
      ddf11ddz__ddz=
     1   ddf11ddz__ddf11dx*ddf11dx__ddz
     1 + ddf11ddz__ddx11ddz*ddx11ddz__ddz
     1 + ddf11ddz__ddf11dy*ddf11dy__ddz
     1 + ddf11ddz__ddy11ddz*ddy11ddz__ddz

      return
      END SUBROUTINE
