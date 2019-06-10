      SUBROUTINE DDDIHEDRAL(
     1 ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz,
     1 dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz, ddx, ddy, ddz,
     1 axax, axay, axaz, axbx, axby, axbz, axcx, axcy, axcz, axdx, axdy,
     1 axdz, ayay, ayaz, aybx, ayby, aybz, aycx, aycy, aycz, aydx, aydy,
     1 aydz, azaz, azbx, azby, azbz, azcx, azcy, azcz, azdx, azdy, azdz,
     1 bxbx, bxby, bxbz, bxcx, bxcy, bxcz, bxdx, bxdy, bxdz, byby, bybz,
     1 bycx, bycy, bycz, bydx, bydy, bydz, bzbz, bzcx, bzcy, bzcz, bzdx,
     1 bzdy, bzdz, cxcx, cxcy, cxcz, cxdx, cxdy, cxdz, cycy, cycz, cydx,
     1 cydy, cydz, czcz, czdx, czdy, czdz, dxdx, dxdy, dxdz, dydy, dydz,
     1 dzdz,
     1 abcd)
      IMPLICIT REAL*8 (a-z)

c In the following routine all vectors ab_x etc are defined reversed.  As a fix,
c all results are reversed as well (yes I know).  The values returned by this
c routine are tested (against mathematica) and correct although the intermediates
c are generally not.

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
c abc \times bc
      aux_x=abc1_y*bc1_z-bc1_y*abc1_z
      aux_y=abc1_z*bc1_x-bc1_z*abc1_x
      aux_z=abc1_x*bc1_y-bc1_x*abc1_y
c two auxiliary reals
c     x=\vec abc1 \cdot \vec bcd_1
      x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
c     y=\vec aux  \cdot \vec bcd_1
      y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
c the result
      dihedral_abcd=datan2(y,x)
c reverse the sign
      abcd = -dihedral_abcd


c//////////////////////////////////////////////////////////////
C//////////////////////////////////////////////////////////////
C// THE FISRT DERIVATIVES
c// to be read from bottom to top
c//////////////////////////////////////////////////////////////
c//////////////////////////////////////////////////////////////

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


c bc1_x=bc_x*bc_length_inv
      dbc1_x__dbx= bc_length_inv + bc_x*dbc_length_inv__dbc_x
      dbc1_x__dby= bc_x*dbc_length_inv__dbc_y
      dbc1_x__dbz= bc_x*dbc_length_inv__dbc_z
      dbc1_x__dcx= -bc_length_inv - bc_x*dbc_length_inv__dbc_x
      dbc1_x__dcy= -bc_x*dbc_length_inv__dbc_y
      dbc1_x__dcz= -bc_x*dbc_length_inv__dbc_z
c bc1_y=bc_y*bc_length_inv
      dbc1_y__dbx= bc_y*dbc_length_inv__dbc_x
      dbc1_y__dby= bc_length_inv + bc_y*dbc_length_inv__dbc_y
      dbc1_y__dbz= bc_y*dbc_length_inv__dbc_z
      dbc1_y__dcx= -bc_y*dbc_length_inv__dbc_x
      dbc1_y__dcy= -bc_length_inv - bc_y*dbc_length_inv__dbc_y
      dbc1_y__dcz= -bc_y*dbc_length_inv__dbc_z
c bc1_z=bc_z*bc_length_inv
      dbc1_z__dbx= bc_z*dbc_length_inv__dbc_x
      dbc1_z__dby= bc_z*dbc_length_inv__dbc_y
      dbc1_z__dbz= bc_length_inv + bc_z*dbc_length_inv__dbc_z
      dbc1_z__dcx= -bc_z*dbc_length_inv__dbc_x
      dbc1_z__dcy= -bc_z*dbc_length_inv__dbc_y
      dbc1_z__dcz= -bc_length_inv - bc_z*dbc_length_inv__dbc_z


c abc_x=-ab_y*bc_z + ab_z*bc_y
      dabc_x__dby=bc_z + ab_z
      dabc_x__dbz=-ab_y - bc_y
c abc_y=-ab_z*bc_x + ab_x*bc_z
      dabc_y__dbx=-ab_z - bc_z
      dabc_y__dbz=bc_x + ab_x
c abc_z=-ab_x*bc_y + ab_y*bc_x
      dabc_z__dbx=ab_y + bc_y
      dabc_z__dby=-ab_x - bc_x 


c abc_x=-ab_y*bc_z + ab_z*bc_y
      dabc_x__dab_y=-bc_z
      dabc_x__dab_z=bc_y
      dabc_x__dbc_y=ab_z
      dabc_x__dbc_z=-ab_y
c abc_y=-ab_z*bc_x + ab_x*bc_z
      dabc_y__dab_x=bc_z
      dabc_y__dab_z=-bc_x
      dabc_y__dbc_x=-ab_z
      dabc_y__dbc_z=ab_x
c abc_z=-ab_x*bc_y + ab_y*bc_x
      dabc_z__dab_x=-bc_y
      dabc_z__dab_y=bc_x
      dabc_z__dbc_x=ab_y
      dabc_z__dbc_y=-ab_x


c abc_x=-ab_y*bc_z + ab_z*bc_y
      dabc_x__day= dabc_x__dab_y*dab_y__day
      dabc_x__daz= dabc_x__dab_z*dab_z__daz
      dabc_x__dby= dabc_x__dab_y*dab_y__dby + dabc_x__dbc_y*dbc_y__dby
      dabc_x__dbz= dabc_x__dbc_z*dbc_z__dbz + dabc_x__dab_z*dab_z__dbz
      dabc_x__dcy= dabc_x__dbc_y*dbc_y__dcy
      dabc_x__dcz= dabc_x__dbc_z*dbc_z__dcz
c abc_y=-ab_z*bc_x + ab_x*bc_z
      dabc_y__dax= dabc_y__dab_x*dab_x__dax
      dabc_y__daz= dabc_y__dab_z*dab_z__daz
      dabc_y__dbx= dabc_y__dbc_x*dbc_x__dbx + dabc_y__dab_x*dab_x__dbx
      dabc_y__dbz= dabc_y__dab_z*dab_z__dbz + dabc_y__dbc_z*dbc_z__dbz
      dabc_y__dcx= dabc_y__dbc_x*dbc_x__dcx
      dabc_y__dcz= dabc_y__dbc_z*dbc_z__dcz
c abc_z=-ab_x*bc_y + ab_y*bc_x
      dabc_z__dax= dabc_z__dab_x*dab_x__dax
      dabc_z__day= dabc_z__dab_y*dab_y__day
      dabc_z__dbx= dabc_z__dab_x*dab_x__dbx + dabc_z__dbc_x*dbc_x__dbx
      dabc_z__dby= dabc_z__dbc_y*dbc_y__dby + dabc_z__dab_y*dab_y__dby
      dabc_z__dcx= dabc_z__dbc_x*dbc_x__dcx
      dabc_z__dcy= dabc_z__dbc_y*dbc_y__dcy


c bcd_x=-bc_y*cd_z + bc_z*cd_y
      dbcd_x__dbc_y=-cd_z
      dbcd_x__dbc_z=cd_y
      dbcd_x__dcd_y=bc_z
      dbcd_x__dcd_z=-bc_y
c bcd_y=-bc_z*cd_x + bc_x*cd_z
      dbcd_y__dbc_x=cd_z
      dbcd_y__dbc_z=-cd_x
      dbcd_y__dcd_x=-bc_z
      dbcd_y__dcd_z=bc_x
c bcd_z=-bc_x*cd_y + bc_y*cd_x
      dbcd_z__dbc_x=-cd_y
      dbcd_z__dbc_y=cd_x
      dbcd_z__dcd_x=bc_y
      dbcd_z__dcd_y=-bc_x


c bcd_x=-bc_y*cd_z + bc_z*cd_y
      dbcd_x__dby= dbcd_x__dbc_y*dbc_y__dby 
      dbcd_x__dbz= dbcd_x__dbc_z*dbc_z__dbz
      dbcd_x__dcy= dbcd_x__dbc_y*dbc_y__dcy + dbcd_x__dcd_y*dcd_y__dcy
      dbcd_x__dcz= dbcd_x__dcd_z*dcd_z__dcz + dbcd_x__dbc_z*dbc_z__dcz 
      dbcd_x__ddy= dbcd_x__dcd_y*dcd_y__ddy
      dbcd_x__ddz= dbcd_x__dcd_z*dcd_z__ddz
c bcd_y=-bc_z*cd_x + bc_x*cd_z
      dbcd_y__dbx= dbcd_y__dbc_x*dbc_x__dbx
      dbcd_y__dbz= dbcd_y__dbc_z*dbc_z__dbz
      dbcd_y__dcx= dbcd_y__dcd_x*dcd_x__dcx + dbcd_y__dbc_x*dbc_x__dcx
      dbcd_y__dcz= dbcd_y__dbc_z*dbc_z__dcz + dbcd_y__dcd_z*dcd_z__dcz
      dbcd_y__ddx= dbcd_y__dcd_x*dcd_x__ddx
      dbcd_y__ddz= dbcd_y__dcd_z*dcd_z__ddz
c bcd_z=-bc_x*cd_y + bc_y*cd_x
      dbcd_z__dbx= dbcd_z__dbc_x*dbc_x__dbx
      dbcd_z__dby= dbcd_z__dbc_y*dbc_y__dby
      dbcd_z__dcx= dbcd_z__dbc_x*dbc_x__dcx + dbcd_z__dcd_x*dcd_x__dcx
      dbcd_z__dcy= dbcd_z__dcd_y*dcd_y__dcy + dbcd_z__dbc_y*dbc_y__dcy
      dbcd_z__ddx= dbcd_z__dcd_x*dcd_x__ddx
      dbcd_z__ddy= dbcd_z__dcd_y*dcd_y__ddy


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


c abc1_x=abc_x*abc_length_inv
c abc1_y=abc_y*abc_length_inv
c abc1_z=abc_z*abc_length_inv
      dabc1_x__dabc_x=abc_length_inv
      dabc1_y__dabc_y=abc_length_inv
      dabc1_z__dabc_z=abc_length_inv

      dabc1_x__dabc_length_inv=abc_x
      dabc1_y__dabc_length_inv=abc_y
      dabc1_z__dabc_length_inv=abc_z


c derivative of the components of the normals
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


c derivative of y
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


c derivative of x
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


c derivative of atan2(y,x) according to x and y
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

c reverse the signs
      dax=-df__dax
      day=-df__day
      daz=-df__daz
      dbx=-df__dbx
      dby=-df__dby
      dbz=-df__dbz
      dcx=-df__dcx
      dcy=-df__dcy
      dcz=-df__dcz
      ddx=-df__ddx
      ddy=-df__ddy
      ddz=-df__ddz
      

c f is the function atan2(y, x)
c fortran (and most other sources) use 'atan2(y,x)' while mathematica uses 'atan2(x,y)'





c////////////////////////////////////////////////////////////////////////
c////////////////////////////////////////////////////////////////////////
c// SECOND DERIVATIVES
c// to be read from bottom to top
c////////////////////////////////////////////////////////////////////////
c////////////////////////////////////////////////////////////////////////


c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      dbc_length_inv__dbx= dbc_length_inv__dbc_x*dbc_x__dbx
      dbc_length_inv__dby= dbc_length_inv__dbc_y*dbc_y__dby
      dbc_length_inv__dbz= dbc_length_inv__dbc_z*dbc_z__dbz
      dbc_length_inv__dcx= dbc_length_inv__dbc_x*dbc_x__dcx
      dbc_length_inv__dcy= dbc_length_inv__dbc_y*dbc_y__dcy
      dbc_length_inv__dcz= dbc_length_inv__dbc_z*dbc_z__dcz


c bcd_length_inv_cub=bcd_length_inv**3
      dbcd_length_inv_cub__dbcd_length_inv=3*bcd_length_inv**2

c bc_length_inv_cub=bc_length_inv**3
      dbc_length_inv_cub__dbc_length_inv=3*bc_length_inv**2

c abc_length_inv_cub=abc_length_inv**3
      dabc_length_inv_cub__dabc_length_inv=3*abc_length_inv**2


c abc_length_inv_cub=abc_length_inv**3
      dabc_length_inv_cub__dax=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dax
      dabc_length_inv_cub__day=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__day
      dabc_length_inv_cub__daz=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__daz
      dabc_length_inv_cub__dbx=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dbx
      dabc_length_inv_cub__dby=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dby
      dabc_length_inv_cub__dbz=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dbz
      dabc_length_inv_cub__dcx=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dcx
      dabc_length_inv_cub__dcy=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dcy
      dabc_length_inv_cub__dcz=
     1   dabc_length_inv_cub__dabc_length_inv*dabc_length_inv__dcz


c dabc_length_inv__dabc_x=-abc_x*abc_length_inv_cub
      ddabc_length_inv11dabc_x__dabc_length_inv_cub=-abc_x
      ddabc_length_inv11dabc_x__dabc_x=-abc_length_inv_cub
c dabc_length_inv__dabc_y=-abc_y*abc_length_inv_cub
      ddabc_length_inv11dabc_y__dabc_length_inv_cub=-abc_y
      ddabc_length_inv11dabc_y__dabc_y=-abc_length_inv_cub
c dabc_length_inv__dabc_z=-abc_z*abc_length_inv_cub
      ddabc_length_inv11dabc_z__dabc_length_inv_cub=-abc_z
      ddabc_length_inv11dabc_z__dabc_z=-abc_length_inv_cub


c bc_length_inv_cub=bc_length_inv**3
      dbc_length_inv_cub__dbx=
     1        dbc_length_inv_cub__dbc_length_inv*dbc_length_inv__dbx
      dbc_length_inv_cub__dby=
     1        dbc_length_inv_cub__dbc_length_inv*dbc_length_inv__dby
      dbc_length_inv_cub__dbz=
     1        dbc_length_inv_cub__dbc_length_inv*dbc_length_inv__dbz
      dbc_length_inv_cub__dcx=
     1        dbc_length_inv_cub__dbc_length_inv*dbc_length_inv__dcx
      dbc_length_inv_cub__dcy=
     1        dbc_length_inv_cub__dbc_length_inv*dbc_length_inv__dcy
      dbc_length_inv_cub__dcz=
     1        dbc_length_inv_cub__dbc_length_inv*dbc_length_inv__dcz


c dbc_length_inv__dbc_x=-bc_x*bc_length_inv_cub
      ddbc_length_inv11dbc_x__dbc_length_inv_cub=-bc_x
      ddbc_length_inv11dbc_x__dbc_x=-bc_length_inv_cub
c dbc_length_inv__dbc_y=-bc_y*bc_length_inv_cub
      ddbc_length_inv11dbc_y__dbc_length_inv_cub=-bc_y
      ddbc_length_inv11dbc_y__dbc_y=-bc_length_inv_cub
c dbc_length_inv__dbc_z=-bc_z*bc_length_inv_cub
      ddbc_length_inv11dbc_z__dbc_length_inv_cub=-bc_z
      ddbc_length_inv11dbc_z__dbc_z=-bc_length_inv_cub


c dbcd_x__dcy=cd_z + bc_z
      ddbcd_x11dcy__dbc_z=1
      ddbcd_x11dcy__dcd_z=1
c dbcd_x__dcz=-bc_y - cd_y 
      ddbcd_x11dcz__dbc_y=-1
      ddbcd_x11dcz__dcd_y=-1
c dbcd_y__dcx=-bc_z - cd_z 
      ddbcd_y11dcx__dbc_z=-1
      ddbcd_y11dcx__dcd_z=-1
c dbcd_y__dcz=cd_x + bc_x
      ddbcd_y11dcz__dbc_x=1
      ddbcd_y11dcz__dcd_x=1
c dbcd_z__dcx=cd_y + bc_y
      ddbcd_z11dcx__dbc_y=1
      ddbcd_z11dcx__dcd_y=1
c dbcd_z__dcy=-bc_x - cd_x
      ddbcd_z11dcy__dbc_x=-1
      ddbcd_z11dcy__dcd_x=-1


c bcd_length_inv_cub=bcd_length_inv**3
      dbcd_length_inv_cub__dbx=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__dbx
      dbcd_length_inv_cub__dby=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__dby
      dbcd_length_inv_cub__dbz=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__dbz
      dbcd_length_inv_cub__dcx=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__dcx
      dbcd_length_inv_cub__dcy=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__dcy
      dbcd_length_inv_cub__dcz=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__dcz
      dbcd_length_inv_cub__ddx=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__ddx
      dbcd_length_inv_cub__ddy=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__ddy
      dbcd_length_inv_cub__ddz=
     1         dbcd_length_inv_cub__dbcd_length_inv*dbcd_length_inv__ddz


c dbcd_length_inv__dbcd_x=-bcd_x*bcd_length_inv_cub
      ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub=-bcd_x
      ddbcd_length_inv11dbcd_x__dbcd_x=-bcd_length_inv_cub
c dbcd_length_inv__dbcd_y=-bcd_y*bcd_length_inv_cub
      ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub=-bcd_y
      ddbcd_length_inv11dbcd_y__dbcd_y=-bcd_length_inv_cub
c dbcd_length_inv__dbcd_z=-bcd_z*bcd_length_inv_cub
      ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub=-bcd_z
      ddbcd_length_inv11dbcd_z__dbcd_z=-bcd_length_inv_cub


c dabc_length_inv__dabc_x=-abc_x*abc_length_inv_cub
      ddabc_length_inv11dabc_x__day=
     1   ddabc_length_inv11dabc_x__dabc_x*dabc_x__day
     1 + ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__day
      ddabc_length_inv11dabc_x__daz=
     1   ddabc_length_inv11dabc_x__dabc_x*dabc_x__daz
     1 + ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__daz
      ddabc_length_inv11dabc_x__dbx=
     1  ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dbx
      ddabc_length_inv11dabc_x__dby=
     1   ddabc_length_inv11dabc_x__dabc_x*dabc_x__dby
     1 + ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dby
      ddabc_length_inv11dabc_x__dbz=
     1   ddabc_length_inv11dabc_x__dabc_x*dabc_x__dbz
     1 + ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dbz
      ddabc_length_inv11dabc_x__dcx=
     1  ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcx
      ddabc_length_inv11dabc_x__dcy=
     1   ddabc_length_inv11dabc_x__dabc_x*dabc_x__dcy
     1 + ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcy
      ddabc_length_inv11dabc_x__dcz=
     1   ddabc_length_inv11dabc_x__dabc_x*dabc_x__dcz
     1 + ddabc_length_inv11dabc_x__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcz
c dabc_length_inv__dabc_y=-abc_y*abc_length_inv_cub
      ddabc_length_inv11dabc_y__dax=
     1   ddabc_length_inv11dabc_y__dabc_y*dabc_y__dax
     1 + ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dax
      ddabc_length_inv11dabc_y__day=
     1  ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__day
      ddabc_length_inv11dabc_y__daz=
     1   ddabc_length_inv11dabc_y__dabc_y*dabc_y__daz
     1 + ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__daz
      ddabc_length_inv11dabc_y__dbx=
     1   ddabc_length_inv11dabc_y__dabc_y*dabc_y__dbx
     1 + ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dbx
      ddabc_length_inv11dabc_y__dby=
     1  ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dby
      ddabc_length_inv11dabc_y__dbz=
     1   ddabc_length_inv11dabc_y__dabc_y*dabc_y__dbz
     1 + ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dbz
      ddabc_length_inv11dabc_y__dcx=
     1   ddabc_length_inv11dabc_y__dabc_y*dabc_y__dcx
     1 + ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcx
      ddabc_length_inv11dabc_y__dcy=
     1  ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcy
      ddabc_length_inv11dabc_y__dcz=
     1   ddabc_length_inv11dabc_y__dabc_y*dabc_y__dcz
     1 + ddabc_length_inv11dabc_y__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcz
c dabc_length_inv__dabc_z=-abc_z*abc_length_inv_cub
      ddabc_length_inv11dabc_z__dax=
     1   ddabc_length_inv11dabc_z__dabc_z*dabc_z__dax
     1 + ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dax
      ddabc_length_inv11dabc_z__day=
     1   ddabc_length_inv11dabc_z__dabc_z*dabc_z__day
     1 + ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__day
      ddabc_length_inv11dabc_z__daz=
     1  ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__daz
      ddabc_length_inv11dabc_z__dbx=
     1   ddabc_length_inv11dabc_z__dabc_z*dabc_z__dbx
     1 + ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dbx
      ddabc_length_inv11dabc_z__dby=
     1   ddabc_length_inv11dabc_z__dabc_z*dabc_z__dby
     1 + ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dby
      ddabc_length_inv11dabc_z__dbz=
     1  ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dbz
      ddabc_length_inv11dabc_z__dcx=
     1   ddabc_length_inv11dabc_z__dabc_z*dabc_z__dcx
     1 + ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcx
      ddabc_length_inv11dabc_z__dcy=
     1   ddabc_length_inv11dabc_z__dabc_z*dabc_z__dcy
     1 + ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcy
      ddabc_length_inv11dabc_z__dcz=
     1  ddabc_length_inv11dabc_z__dabc_length_inv_cub*
     1                               dabc_length_inv_cub__dcz


c dabc_length_inv__dax=dabc_length_inv__dabc_y*bc_z - dabc_length_inv__dabc_z*bc_y
      ddabc_length_inv11dax__dbc_y=-dabc_length_inv__dabc_z
      ddabc_length_inv11dax__dbc_z=dabc_length_inv__dabc_y
      ddabc_length_inv11dax__ddabc_length_inv11dabc_y=bc_z
      ddabc_length_inv11dax__ddabc_length_inv11dabc_z=-bc_y
c dabc_length_inv__day=-dabc_length_inv__dabc_x*bc_z + dabc_length_inv__dabc_z*bc_x
      ddabc_length_inv11day__dbc_x=dabc_length_inv__dabc_z
      ddabc_length_inv11day__dbc_z=-dabc_length_inv__dabc_x
      ddabc_length_inv11day__ddabc_length_inv11dabc_x=-bc_z
      ddabc_length_inv11day__ddabc_length_inv11dabc_z=bc_x
c dabc_length_inv__daz=dabc_length_inv__dabc_x*bc_y - dabc_length_inv__dabc_y*bc_x
      ddabc_length_inv11daz__dbc_x=-dabc_length_inv__dabc_y
      ddabc_length_inv11daz__dbc_y=dabc_length_inv__dabc_x
      ddabc_length_inv11daz__ddabc_length_inv11dabc_x=bc_y
      ddabc_length_inv11daz__ddabc_length_inv11dabc_y=-bc_x
c dabc_length_inv__dbx=dabc_length_inv__dabc_y*dabc_y__dbx + dabc_length_inv__dabc_z*dabc_z__dbx
      ddabc_length_inv11dbx__ddabc_length_inv11dabc_y=dabc_y__dbx
      ddabc_length_inv11dbx__ddabc_length_inv11dabc_z=dabc_z__dbx
      ddabc_length_inv11dbx__ddabc_y11dbx=dabc_length_inv__dabc_y
      ddabc_length_inv11dbx__ddabc_z11dbx=dabc_length_inv__dabc_z
c dabc_length_inv__dby=dabc_length_inv__dabc_x*dabc_x__dby + dabc_length_inv__dabc_z*dabc_z__dby
      ddabc_length_inv11dby__ddabc_length_inv11dabc_x=dabc_x__dby
      ddabc_length_inv11dby__ddabc_length_inv11dabc_z=dabc_z__dby
      ddabc_length_inv11dby__ddabc_x11dby=dabc_length_inv__dabc_x
      ddabc_length_inv11dby__ddabc_z11dby=dabc_length_inv__dabc_z
c dabc_length_inv__dbz=dabc_length_inv__dabc_x*dabc_x__dbz + dabc_length_inv__dabc_y*dabc_y__dbz
      ddabc_length_inv11dbz__ddabc_length_inv11dabc_x=dabc_x__dbz
      ddabc_length_inv11dbz__ddabc_length_inv11dabc_y=dabc_y__dbz
      ddabc_length_inv11dbz__ddabc_x11dbz=dabc_length_inv__dabc_x
      ddabc_length_inv11dbz__ddabc_y11dbz=dabc_length_inv__dabc_y
c dabc_length_inv__dcx=dabc_length_inv__dabc_y*ab_z - dabc_length_inv__dabc_z*ab_y
      ddabc_length_inv11dcx__dab_y=-dabc_length_inv__dabc_z
      ddabc_length_inv11dcx__dab_z=dabc_length_inv__dabc_y
      ddabc_length_inv11dcx__ddabc_length_inv11dabc_y=ab_z
      ddabc_length_inv11dcx__ddabc_length_inv11dabc_z=-ab_y
c dabc_length_inv__dcy=-dabc_length_inv__dabc_x*ab_z + dabc_length_inv__dabc_z*ab_x
      ddabc_length_inv11dcy__dab_x=dabc_length_inv__dabc_z
      ddabc_length_inv11dcy__dab_z=-dabc_length_inv__dabc_x
      ddabc_length_inv11dcy__ddabc_length_inv11dabc_x=-ab_z
      ddabc_length_inv11dcy__ddabc_length_inv11dabc_z=ab_x
c dabc_length_inv__dcz=dabc_length_inv__dabc_x*ab_y - dabc_length_inv__dabc_y*ab_x
      ddabc_length_inv11dcz__dab_x=-dabc_length_inv__dabc_y
      ddabc_length_inv11dcz__dab_y=dabc_length_inv__dabc_x
      ddabc_length_inv11dcz__ddabc_length_inv11dabc_x=ab_y
      ddabc_length_inv11dcz__ddabc_length_inv11dabc_y=-ab_x


c dabc_x__dby=bc_z + ab_z
      ddabc_x11dby__dab_z=1
      ddabc_x11dby__dbc_z=1
c dabc_x__dbz=-ab_y - bc_y
      ddabc_x11dbz__dab_y=-1
      ddabc_x11dbz__dbc_y=-1
c dabc_y__dbx=-ab_z - bc_z
      ddabc_y11dbx__dab_z=-1
      ddabc_y11dbx__dbc_z=-1
c dabc_y__dbz=bc_x + ab_x
      ddabc_y11dbz__dab_x=1
      ddabc_y11dbz__dbc_x=1
c dabc_z__dbx=ab_y + bc_y
      ddabc_z11dbx__dab_y=1
      ddabc_z11dbx__dbc_y=1
c dabc_z__dby=-ab_x - bc_x 
      ddabc_z11dby__dab_x=-1
      ddabc_z11dby__dbc_x=-1


c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
      dbc_length_inv__dbx=dbc_length_inv__dbc_x*dbc_x__dbx
      dbc_length_inv__dby=dbc_length_inv__dbc_y*dbc_y__dby
      dbc_length_inv__dbz=dbc_length_inv__dbc_z*dbc_z__dbz
      dbc_length_inv__dcx=dbc_length_inv__dbc_x*dbc_x__dcx
      dbc_length_inv__dcy=dbc_length_inv__dbc_y*dbc_y__dcy
      dbc_length_inv__dcz=dbc_length_inv__dbc_z*dbc_z__dcz


c dbc_length_inv__dbc_x=-bc_x*bc_length_inv_cub
      ddbc_length_inv11dbc_x__dbx=
     1   ddbc_length_inv11dbc_x__dbc_x*dbc_x__dbx
     1 + ddbc_length_inv11dbc_x__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dbx
      ddbc_length_inv11dbc_x__dby=
     1  ddbc_length_inv11dbc_x__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dby
      ddbc_length_inv11dbc_x__dbz=
     1  ddbc_length_inv11dbc_x__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dbz
      ddbc_length_inv11dbc_x__dcx=
     1   ddbc_length_inv11dbc_x__dbc_x*dbc_x__dcx
     1 + ddbc_length_inv11dbc_x__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcx
      ddbc_length_inv11dbc_x__dcy=
     1  ddbc_length_inv11dbc_x__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcy
      ddbc_length_inv11dbc_x__dcz=
     1  ddbc_length_inv11dbc_x__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcz
c dbc_length_inv__dbc_y=-bc_y*bc_length_inv_cub
      ddbc_length_inv11dbc_y__dby=
     1   ddbc_length_inv11dbc_y__dbc_y*dbc_y__dby
     1 + ddbc_length_inv11dbc_y__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dby
      ddbc_length_inv11dbc_y__dbz=
     1  ddbc_length_inv11dbc_y__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dbz
      ddbc_length_inv11dbc_y__dcx=
     1  ddbc_length_inv11dbc_y__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcx
      ddbc_length_inv11dbc_y__dcy=
     1   ddbc_length_inv11dbc_y__dbc_y*dbc_y__dcy
     1 + ddbc_length_inv11dbc_y__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcy
      ddbc_length_inv11dbc_y__dcz=
     1  ddbc_length_inv11dbc_y__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcz
c dbc_length_inv__dbc_z=-bc_z*bc_length_inv_cub
      ddbc_length_inv11dbc_z__dbz=
     1   ddbc_length_inv11dbc_z__dbc_z*dbc_z__dbz
     1 + ddbc_length_inv11dbc_z__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dbz
      ddbc_length_inv11dbc_z__dcx=
     1  ddbc_length_inv11dbc_z__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcx
      ddbc_length_inv11dbc_z__dcy=
     1  ddbc_length_inv11dbc_z__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcy
      ddbc_length_inv11dbc_z__dcz=
     1   ddbc_length_inv11dbc_z__dbc_z*dbc_z__dcz
     1 + ddbc_length_inv11dbc_z__dbc_length_inv_cub*
     1                                     dbc_length_inv_cub__dcz


c dbc1_x__dbx= bc_length_inv + bc_x*dbc_length_inv__dbc_x
      ddbc1_x11dbx__dbc_length_inv=1
      ddbc1_x11dbx__dbc_x=dbc_length_inv__dbc_x
      ddbc1_x11dbx__ddbc_length_inv11dbc_x=bc_x
c dbc1_x__dby= bc_x*dbc_length_inv__dbc_y
      ddbc1_x11dby__dbc_x=dbc_length_inv__dbc_y
      ddbc1_x11dby__ddbc_length_inv11dbc_y=bc_x
c dbc1_x__dbz= bc_x*dbc_length_inv__dbc_z
      ddbc1_x11dbz__dbc_x=dbc_length_inv__dbc_z
      ddbc1_x11dbz__ddbc_length_inv11dbc_z=bc_x
c dbc1_x__dcx= -bc_length_inv - bc_x*dbc_length_inv__dbc_x
      ddbc1_x11dcx__dbc_length_inv=-1
      ddbc1_x11dcx__dbc_x=-dbc_length_inv__dbc_x
      ddbc1_x11dcx__ddbc_length_inv11dbc_x=-bc_x
c dbc1_x__dcy= -bc_x*dbc_length_inv__dbc_y
      ddbc1_x11dcy__dbc_x=-dbc_length_inv__dbc_y
      ddbc1_x11dcy__ddbc_length_inv11dbc_y=-bc_x
c dbc1_x__dcz= -bc_x*dbc_length_inv__dbc_z
      ddbc1_x11dcz__dbc_x=-dbc_length_inv__dbc_z
      ddbc1_x11dcz__ddbc_length_inv11dbc_z=-bc_x
c dbc1_y__dbx= bc_y*dbc_length_inv__dbc_x
      ddbc1_y11dbx__dbc_y=dbc_length_inv__dbc_x
      ddbc1_y11dbx__ddbc_length_inv11dbc_x=bc_y
c dbc1_y__dby= bc_length_inv + bc_y*dbc_length_inv__dbc_y
      ddbc1_y11dby__dbc_length_inv=1
      ddbc1_y11dby__dbc_y=dbc_length_inv__dbc_y
      ddbc1_y11dby__ddbc_length_inv11dbc_y=bc_y
c dbc1_y__dbz= bc_y*dbc_length_inv__dbc_z
      ddbc1_y11dbz__dbc_y=dbc_length_inv__dbc_z
      ddbc1_y11dbz__ddbc_length_inv11dbc_z=bc_y
c dbc1_y__dcx= -bc_y*dbc_length_inv__dbc_x
      ddbc1_y11dcx__dbc_y=-dbc_length_inv__dbc_x
      ddbc1_y11dcx__ddbc_length_inv11dbc_x=-bc_y
c dbc1_y__dcy= -bc_length_inv - bc_y*dbc_length_inv__dbc_y
      ddbc1_y11dcy__dbc_length_inv=-1
      ddbc1_y11dcy__dbc_y=-dbc_length_inv__dbc_y
      ddbc1_y11dcy__ddbc_length_inv11dbc_y=-bc_y
c dbc1_y__dcz= -bc_y*dbc_length_inv__dbc_z
      ddbc1_y11dcz__dbc_y=-dbc_length_inv__dbc_z
      ddbc1_y11dcz__ddbc_length_inv11dbc_z=-bc_y
c dbc1_z__dbx= bc_z*dbc_length_inv__dbc_x
      ddbc1_z11dbx__dbc_z=dbc_length_inv__dbc_x
      ddbc1_z11dbx__ddbc_length_inv11dbc_x=bc_z
c dbc1_z__dby= bc_z*dbc_length_inv__dbc_y
      ddbc1_z11dby__dbc_z=dbc_length_inv__dbc_y
      ddbc1_z11dby__ddbc_length_inv11dbc_y=bc_z
c dbc1_z__dbz= bc_length_inv + bc_z*dbc_length_inv__dbc_z
      ddbc1_z11dbz__dbc_length_inv=1
      ddbc1_z11dbz__dbc_z=dbc_length_inv__dbc_z
      ddbc1_z11dbz__ddbc_length_inv11dbc_z=bc_z
c dbc1_z__dcx= -bc_z*dbc_length_inv__dbc_x
      ddbc1_z11dcx__dbc_z=-dbc_length_inv__dbc_x
      ddbc1_z11dcx__ddbc_length_inv11dbc_x=-bc_z
c dbc1_z__dcy= -bc_z*dbc_length_inv__dbc_y
      ddbc1_z11dcy__dbc_z=-dbc_length_inv__dbc_y
      ddbc1_z11dcy__ddbc_length_inv11dbc_y=-bc_z
c dbc1_z__dcz= -bc_length_inv - bc_z*dbc_length_inv__dbc_z
      ddbc1_z11dcz__dbc_length_inv=-1
      ddbc1_z11dcz__dbc_z=-dbc_length_inv__dbc_z
      ddbc1_z11dcz__ddbc_length_inv11dbc_z=-bc_z


c dbcd_x__dcy=cd_z + bc_z
      ddbcd_x11dcy__dcz=
     1   ddbcd_x11dcy__dcd_z*dcd_z__dcz
     1 + ddbcd_x11dcy__dbc_z*dbc_z__dcz
      ddbcd_x11dcy__ddz=
     1   ddbcd_x11dcy__dcd_z*dcd_z__ddz
c dbcd_x__dcz=-bc_y - cd_y 
      ddbcd_x11dcz__ddy=
     1  ddbcd_x11dcz__dcd_y*dcd_y__ddy
c dbcd_y__dcx=-bc_z - cd_z 
      ddbcd_y11dcx__dcz=
     1    ddbcd_y11dcx__dbc_z*dbc_z__dcz
     1 +  ddbcd_y11dcx__dcd_z*dcd_z__dcz
      ddbcd_y11dcx__ddz=
     1   ddbcd_y11dcx__dcd_z*dcd_z__ddz
c dbcd_y__dcz=cd_x + bc_x
      ddbcd_y11dcz__ddx=
     1    ddbcd_y11dcz__dcd_x*dcd_x__ddx
c dbcd_z__dcx=cd_y + bc_y
      ddbcd_z11dcx__dcy=
     1   ddbcd_z11dcx__dcd_y*dcd_y__dcy
     1 + ddbcd_z11dcx__dbc_y*dbc_y__dcy
      ddbcd_z11dcx__ddy=
     1   ddbcd_z11dcx__dcd_y*dcd_y__ddy
c dbcd_z__dcy=-bc_x - cd_x
      ddbcd_z11dcy__ddx=
     1  ddbcd_z11dcy__dcd_x*dcd_x__ddx


c dbcd_length_inv__dbcd_x=-bcd_x*bcd_length_inv_cub
      ddbcd_length_inv11dbcd_x__dby=
     1   ddbcd_length_inv11dbcd_x__dbcd_x*dbcd_x__dby
     1 + ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dby
      ddbcd_length_inv11dbcd_x__dbz=
     1   ddbcd_length_inv11dbcd_x__dbcd_x*dbcd_x__dbz
     1 + ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dbz
      ddbcd_length_inv11dbcd_x__dcx=
     1  ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcx
      ddbcd_length_inv11dbcd_x__dcy=
     1   ddbcd_length_inv11dbcd_x__dbcd_x*dbcd_x__dcy
     1 + ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcy
      ddbcd_length_inv11dbcd_x__dcz=
     1   ddbcd_length_inv11dbcd_x__dbcd_x*dbcd_x__dcz
     1 + ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcz
      ddbcd_length_inv11dbcd_x__ddx=
     1  ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddx
      ddbcd_length_inv11dbcd_x__ddy=
     1   ddbcd_length_inv11dbcd_x__dbcd_x*dbcd_x__ddy
     1 + ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddy
      ddbcd_length_inv11dbcd_x__ddz=
     1   ddbcd_length_inv11dbcd_x__dbcd_x*dbcd_x__ddz
     1 + ddbcd_length_inv11dbcd_x__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddz
c dbcd_length_inv__dbcd_y=-bcd_y*bcd_length_inv_cub
      ddbcd_length_inv11dbcd_y__dbx=
     1   ddbcd_length_inv11dbcd_y__dbcd_y*dbcd_y__dbx
     1 + ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dbx
      ddbcd_length_inv11dbcd_y__dby=
     1  ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dby
      ddbcd_length_inv11dbcd_y__dbz=
     1   ddbcd_length_inv11dbcd_y__dbcd_y*dbcd_y__dbz
     1 + ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dbz
      ddbcd_length_inv11dbcd_y__dcx=
     1   ddbcd_length_inv11dbcd_y__dbcd_y*dbcd_y__dcx
     1 + ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcx
      ddbcd_length_inv11dbcd_y__dcy=
     1  ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcy
      ddbcd_length_inv11dbcd_y__dcz=
     1   ddbcd_length_inv11dbcd_y__dbcd_y*dbcd_y__dcz
     1 + ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcz
      ddbcd_length_inv11dbcd_y__ddx=
     1   ddbcd_length_inv11dbcd_y__dbcd_y*dbcd_y__ddx
     1 + ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddx
      ddbcd_length_inv11dbcd_y__ddy=
     1  ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddy
      ddbcd_length_inv11dbcd_y__ddz=
     1   ddbcd_length_inv11dbcd_y__dbcd_y*dbcd_y__ddz
     1 + ddbcd_length_inv11dbcd_y__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddz
c dbcd_length_inv__dbcd_z=-bcd_z*bcd_length_inv_cub
      ddbcd_length_inv11dbcd_z__dbx=
     1   ddbcd_length_inv11dbcd_z__dbcd_z*dbcd_z__dbx
     1 + ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dbx
      ddbcd_length_inv11dbcd_z__dby=
     1   ddbcd_length_inv11dbcd_z__dbcd_z*dbcd_z__dby
     1 + ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dby
      ddbcd_length_inv11dbcd_z__dbz=
     1  ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dbz
      ddbcd_length_inv11dbcd_z__dcx=
     1   ddbcd_length_inv11dbcd_z__dbcd_z*dbcd_z__dcx
     1 + ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcx
      ddbcd_length_inv11dbcd_z__dcy=
     1   ddbcd_length_inv11dbcd_z__dbcd_z*dbcd_z__dcy
     1 + ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcy
      ddbcd_length_inv11dbcd_z__dcz=
     1  ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__dcz
      ddbcd_length_inv11dbcd_z__ddx=
     1   ddbcd_length_inv11dbcd_z__dbcd_z*dbcd_z__ddx
     1 + ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddx
      ddbcd_length_inv11dbcd_z__ddy=
     1   ddbcd_length_inv11dbcd_z__dbcd_z*dbcd_z__ddy
     1 + ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddy
      ddbcd_length_inv11dbcd_z__ddz=
     1  ddbcd_length_inv11dbcd_z__dbcd_length_inv_cub*
     1                               dbcd_length_inv_cub__ddz


c dbcd_length_inv__dbx=dbcd_length_inv__dbcd_y*cd_z - dbcd_length_inv__dbcd_z*cd_y
      ddbcd_length_inv11dbx__dcd_y=-dbcd_length_inv__dbcd_z
      ddbcd_length_inv11dbx__dcd_z=dbcd_length_inv__dbcd_y
      ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y=cd_z
      ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z=-cd_y
c dbcd_length_inv__dby=-dbcd_length_inv__dbcd_x*cd_z + dbcd_length_inv__dbcd_z*cd_x
      ddbcd_length_inv11dby__dcd_x=dbcd_length_inv__dbcd_z
      ddbcd_length_inv11dby__dcd_z=-dbcd_length_inv__dbcd_x
      ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x=-cd_z
      ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z=cd_x
c dbcd_length_inv__dbz=dbcd_length_inv__dbcd_x*cd_y - dbcd_length_inv__dbcd_y*cd_x
      ddbcd_length_inv11dbz__dcd_x=-dbcd_length_inv__dbcd_y
      ddbcd_length_inv11dbz__dcd_y=dbcd_length_inv__dbcd_x
      ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x=cd_y
      ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y=-cd_x
c dbcd_length_inv__dcx=dbcd_length_inv__dbcd_y*dbcd_y__dcx + dbcd_length_inv__dbcd_z*dbcd_z__dcx
      ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y=dbcd_y__dcx
      ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z=dbcd_z__dcx
      ddbcd_length_inv11dcx__ddbcd_y11dcx=dbcd_length_inv__dbcd_y
      ddbcd_length_inv11dcx__ddbcd_z11dcx=dbcd_length_inv__dbcd_z
c dbcd_length_inv__dcy=dbcd_length_inv__dbcd_x*dbcd_x__dcy + dbcd_length_inv__dbcd_z*dbcd_z__dcy
      ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_x=dbcd_x__dcy
      ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_z=dbcd_z__dcy
      ddbcd_length_inv11dcy__ddbcd_x11dcy=dbcd_length_inv__dbcd_x
      ddbcd_length_inv11dcy__ddbcd_z11dcy=dbcd_length_inv__dbcd_z
c dbcd_length_inv__dcz=dbcd_length_inv__dbcd_x*dbcd_x__dcz + dbcd_length_inv__dbcd_y*dbcd_y__dcz
      ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_x=dbcd_x__dcz
      ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_y=dbcd_y__dcz
      ddbcd_length_inv11dcz__ddbcd_x11dcz=dbcd_length_inv__dbcd_x
      ddbcd_length_inv11dcz__ddbcd_y11dcz=dbcd_length_inv__dbcd_y
c dbcd_length_inv__ddx=dbcd_length_inv__dbcd_y*bc_z - dbcd_length_inv__dbcd_z*bc_y
      ddbcd_length_inv11ddx__dbc_y=-dbcd_length_inv__dbcd_z
      ddbcd_length_inv11ddx__dbc_z=dbcd_length_inv__dbcd_y
      ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_y=bc_z
      ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_z=-bc_y
c dbcd_length_inv__ddy=-dbcd_length_inv__dbcd_x*bc_z + dbcd_length_inv__dbcd_z*bc_x
      ddbcd_length_inv11ddy__dbc_x=dbcd_length_inv__dbcd_z
      ddbcd_length_inv11ddy__dbc_z=-dbcd_length_inv__dbcd_x
      ddbcd_length_inv11ddy__ddbcd_length_inv11dbcd_x=-bc_z
      ddbcd_length_inv11ddy__ddbcd_length_inv11dbcd_z=bc_x
c dbcd_length_inv__ddz=dbcd_length_inv__dbcd_x*bc_y - dbcd_length_inv__dbcd_y*bc_x
      ddbcd_length_inv11ddz__dbc_x=-dbcd_length_inv__dbcd_y
      ddbcd_length_inv11ddz__dbc_y=dbcd_length_inv__dbcd_x
      ddbcd_length_inv11ddz__ddbcd_length_inv11dbcd_x=bc_y
      ddbcd_length_inv11ddz__ddbcd_length_inv11dbcd_y=-bc_x


c dabc_x__dby=bc_z + ab_z
      ddabc_x11dby__dbz=
     1   ddabc_x11dby__dbc_z*dbc_z__dbz + ddabc_x11dby__dab_z*dab_z__dbz
      ddabc_x11dby__dcz=
     1   ddabc_x11dby__dbc_z*dbc_z__dcz 
c dabc_x__dbz=-ab_y - bc_y
      ddabc_x11dbz__dcy=
     1   ddabc_x11dbz__dbc_y*dbc_y__dcy

c dabc_y__dbx=-ab_z - bc_z
      ddabc_y11dbx__dbz=
     1   ddabc_y11dbx__dab_z*dab_z__dbz + ddabc_y11dbx__dbc_z*dbc_z__dbz
      ddabc_y11dbx__dcz=
     1    ddabc_y11dbx__dbc_z*dbc_z__dcz
c dabc_y__dbz=bc_x + ab_x
      ddabc_y11dbz__dcx=
     1   ddabc_y11dbz__dbc_x*dbc_x__dcx 

c dabc_z__dbx=ab_y + bc_y
      ddabc_z11dbx__dby=
     1   ddabc_z11dbx__dab_y*dab_y__dby + ddabc_z11dbx__dbc_y*dbc_y__dby
      ddabc_z11dbx__dcy=
     1    ddabc_z11dbx__dbc_y*dbc_y__dcy
c dabc_z__dby=-ab_x - bc_x 
      ddabc_z11dby__dcx=
     1   ddabc_z11dby__dbc_x*dbc_x__dcx


c dabc_length_inv__dax=dabc_length_inv__dabc_y*bc_z - dabc_length_inv__dabc_z*bc_y
      ddabc_length_inv11dax__dax=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dax
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dax
      ddabc_length_inv11dax__day=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__day
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__day
      ddabc_length_inv11dax__daz=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__daz
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__daz
      ddabc_length_inv11dax__dbx=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbx
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbx
      ddabc_length_inv11dax__dby=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dby
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dby
     1 + ddabc_length_inv11dax__dbc_y*dbc_y__dby
      ddabc_length_inv11dax__dbz=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbz
     1 + ddabc_length_inv11dax__dbc_z*dbc_z__dbz
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbz
      ddabc_length_inv11dax__dcx=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcx
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcx
      ddabc_length_inv11dax__dcy=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcy
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcy
     1 + ddabc_length_inv11dax__dbc_y*dbc_y__dcy
      ddabc_length_inv11dax__dcz=
     1   ddabc_length_inv11dax__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcz
     1 + ddabc_length_inv11dax__dbc_z*dbc_z__dcz
     1 + ddabc_length_inv11dax__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcz
c dabc_length_inv__day=-dabc_length_inv__dabc_x*bc_z + dabc_length_inv__dabc_z*bc_x
      ddabc_length_inv11day__day=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__day
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__day
      ddabc_length_inv11day__daz=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__daz
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__daz
      ddabc_length_inv11day__dbx=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dbx
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbx
     1 + ddabc_length_inv11day__dbc_x*dbc_x__dbx
      ddabc_length_inv11day__dby=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dby
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dby
      ddabc_length_inv11day__dbz=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dbz
     1 + ddabc_length_inv11day__dbc_z*dbc_z__dbz
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbz
      ddabc_length_inv11day__dcx=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcx
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcx
     1 + ddabc_length_inv11day__dbc_x*dbc_x__dcx
      ddabc_length_inv11day__dcy=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcy
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcy
      ddabc_length_inv11day__dcz=
     1   ddabc_length_inv11day__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcz
     1 + ddabc_length_inv11day__dbc_z*dbc_z__dcz
     1 + ddabc_length_inv11day__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcz
c dabc_length_inv__daz=dabc_length_inv__dabc_x*bc_y - dabc_length_inv__dabc_y*bc_x
      ddabc_length_inv11daz__daz=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__daz
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__daz
      ddabc_length_inv11daz__dbx=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dbx
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbx
     1 + ddabc_length_inv11daz__dbc_x*dbc_x__dbx
      ddabc_length_inv11daz__dby=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dby
     1 + ddabc_length_inv11daz__dbc_y*dbc_y__dby
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dby
      ddabc_length_inv11daz__dbz=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dbz
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbz
      ddabc_length_inv11daz__dcx=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcx
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcx
     1 + ddabc_length_inv11daz__dbc_x*dbc_x__dcx
      ddabc_length_inv11daz__dcy=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcy
     1 + ddabc_length_inv11daz__dbc_y*dbc_y__dcy
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcy
      ddabc_length_inv11daz__dcz=
     1   ddabc_length_inv11daz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcz
     1 + ddabc_length_inv11daz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcz
c dabc_length_inv__dbx=dabc_length_inv__dabc_y*dabc_y__dbx + dabc_length_inv__dabc_z*dabc_z__dbx
      ddabc_length_inv11dbx__dbx=
     1   ddabc_length_inv11dbx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbx
     1 + ddabc_length_inv11dbx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbx
      ddabc_length_inv11dbx__dby=
     1   ddabc_length_inv11dbx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dby
     1 + ddabc_length_inv11dbx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dby
     1 + ddabc_length_inv11dbx__ddabc_z11dbx*ddabc_z11dbx__dby
      ddabc_length_inv11dbx__dbz=
     1   ddabc_length_inv11dbx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbz
     1 + ddabc_length_inv11dbx__ddabc_y11dbx*ddabc_y11dbx__dbz
     1 + ddabc_length_inv11dbx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbz
      ddabc_length_inv11dbx__dcx=
     1   ddabc_length_inv11dbx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcx
     1 + ddabc_length_inv11dbx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcx
      ddabc_length_inv11dbx__dcy=
     1   ddabc_length_inv11dbx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcy
     1 + ddabc_length_inv11dbx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcy
     1 + ddabc_length_inv11dbx__ddabc_z11dbx*ddabc_z11dbx__dcy
      ddabc_length_inv11dbx__dcz=
     1   ddabc_length_inv11dbx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcz
     1 + ddabc_length_inv11dbx__ddabc_y11dbx*ddabc_y11dbx__dcz
     1 + ddabc_length_inv11dbx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcz
c dabc_length_inv__dby=dabc_length_inv__dabc_x*dabc_x__dby + dabc_length_inv__dabc_z*dabc_z__dby
      ddabc_length_inv11dby__dby=
     1   ddabc_length_inv11dby__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dby
     1 + ddabc_length_inv11dby__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dby
      ddabc_length_inv11dby__dbz=
     1   ddabc_length_inv11dby__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dbz
     1 + ddabc_length_inv11dby__ddabc_x11dby*ddabc_x11dby__dbz
     1 + ddabc_length_inv11dby__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dbz
      ddabc_length_inv11dby__dcx=
     1   ddabc_length_inv11dby__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcx
     1 + ddabc_length_inv11dby__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcx
     1 + ddabc_length_inv11dby__ddabc_z11dby*ddabc_z11dby__dcx
      ddabc_length_inv11dby__dcy=
     1   ddabc_length_inv11dby__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcy
     1 + ddabc_length_inv11dby__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcy
      ddabc_length_inv11dby__dcz=
     1   ddabc_length_inv11dby__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcz
     1 + ddabc_length_inv11dby__ddabc_x11dby*ddabc_x11dby__dcz
     1 + ddabc_length_inv11dby__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcz
c dabc_length_inv__dbz=dabc_length_inv__dabc_x*dabc_x__dbz + dabc_length_inv__dabc_y*dabc_y__dbz
      ddabc_length_inv11dbz__dbz=
     1   ddabc_length_inv11dbz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dbz
     1 + ddabc_length_inv11dbz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dbz
      ddabc_length_inv11dbz__dcx=
     1   ddabc_length_inv11dbz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcx
     1 + ddabc_length_inv11dbz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcx
     1 + ddabc_length_inv11dbz__ddabc_y11dbz*ddabc_y11dbz__dcx
      ddabc_length_inv11dbz__dcy=
     1   ddabc_length_inv11dbz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcy
     1 + ddabc_length_inv11dbz__ddabc_x11dbz*ddabc_x11dbz__dcy
     1 + ddabc_length_inv11dbz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcy
      ddabc_length_inv11dbz__dcz=
     1   ddabc_length_inv11dbz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcz
     1 + ddabc_length_inv11dbz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcz
c dabc_length_inv__dcx=dabc_length_inv__dabc_y*ab_z - dabc_length_inv__dabc_z*ab_y
      ddabc_length_inv11dcx__dcx=
     1   ddabc_length_inv11dcx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcx
     1 + ddabc_length_inv11dcx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcx
      ddabc_length_inv11dcx__dcy=
     1   ddabc_length_inv11dcx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcy
     1 + ddabc_length_inv11dcx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcy
      ddabc_length_inv11dcx__dcz=
     1   ddabc_length_inv11dcx__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcz
     1 + ddabc_length_inv11dcx__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcz
c dabc_length_inv__dcy=-dabc_length_inv__dabc_x*ab_z + dabc_length_inv__dabc_z*ab_x
      ddabc_length_inv11dcy__dcy=
     1   ddabc_length_inv11dcy__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcy
     1 + ddabc_length_inv11dcy__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcy
      ddabc_length_inv11dcy__dcz=
     1   ddabc_length_inv11dcy__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcz
     1 + ddabc_length_inv11dcy__ddabc_length_inv11dabc_z*
     1                              ddabc_length_inv11dabc_z__dcz
c dabc_length_inv__dcz=dabc_length_inv__dabc_x*ab_y - dabc_length_inv__dabc_y*ab_x
      ddabc_length_inv11dcz__dcz=
     1   ddabc_length_inv11dcz__ddabc_length_inv11dabc_x*
     1                              ddabc_length_inv11dabc_x__dcz
     1 + ddabc_length_inv11dcz__ddabc_length_inv11dabc_y*
     1                              ddabc_length_inv11dabc_y__dcz


c dabc1_x__dax= abc_x*dabc_length_inv__dax
      ddabc1_x11dax__dabc_x=dabc_length_inv__dax
      ddabc1_x11dax__ddabc_length_inv11dax=abc_x
c dabc1_x__day=- abc_length_inv*bc_z + abc_x*dabc_length_inv__day
      ddabc1_x11day__dabc_length_inv=-bc_z
      ddabc1_x11day__dabc_x=dabc_length_inv__day
      ddabc1_x11day__dbc_z=-abc_length_inv
      ddabc1_x11day__ddabc_length_inv11day=abc_x
c dabc1_x__daz=abc_length_inv*bc_y + abc_x*dabc_length_inv__daz
      ddabc1_x11daz__dabc_length_inv=bc_y
      ddabc1_x11daz__dabc_x=dabc_length_inv__daz
      ddabc1_x11daz__dbc_y=abc_length_inv
      ddabc1_x11daz__ddabc_length_inv11daz=abc_x
c dabc1_x__dbx= abc_x*dabc_length_inv__dbx
      ddabc1_x11dbx__dabc_x=dabc_length_inv__dbx
      ddabc1_x11dbx__ddabc_length_inv11dbx=abc_x
c dabc1_x__dby=abc_length_inv*dabc_x__dby + abc_x*dabc_length_inv__dby
      ddabc1_x11dby__dabc_length_inv=dabc_x__dby
      ddabc1_x11dby__dabc_x=dabc_length_inv__dby
      ddabc1_x11dby__ddabc_length_inv11dby=abc_x
      ddabc1_x11dby__ddabc_x11dby=abc_length_inv
c dabc1_x__dbz=abc_length_inv*dabc_x__dbz + abc_x*dabc_length_inv__dbz
      ddabc1_x11dbz__dabc_length_inv=dabc_x__dbz
      ddabc1_x11dbz__dabc_x=dabc_length_inv__dbz
      ddabc1_x11dbz__ddabc_length_inv11dbz=abc_x
      ddabc1_x11dbz__ddabc_x11dbz=abc_length_inv
c dabc1_x__dcx= abc_x*dabc_length_inv__dcx
      ddabc1_x11dcx__dabc_x=dabc_length_inv__dcx
      ddabc1_x11dcx__ddabc_length_inv11dcx=abc_x
c dabc1_x__dcy=-abc_length_inv*ab_z + abc_x*dabc_length_inv__dcy
      ddabc1_x11dcy__dab_z=-abc_length_inv
      ddabc1_x11dcy__dabc_length_inv=-ab_z
      ddabc1_x11dcy__dabc_x=dabc_length_inv__dcy
      ddabc1_x11dcy__ddabc_length_inv11dcy=abc_x
c dabc1_x__dcz=abc_length_inv*ab_y + abc_x*dabc_length_inv__dcz
      ddabc1_x11dcz__dab_y=abc_length_inv
      ddabc1_x11dcz__dabc_length_inv=ab_y
      ddabc1_x11dcz__dabc_x=dabc_length_inv__dcz
      ddabc1_x11dcz__ddabc_length_inv11dcz=abc_x

c dabc1_y__dax=abc_length_inv*bc_z + abc_y*dabc_length_inv__dax
      ddabc1_y11dax__dabc_length_inv=bc_z
      ddabc1_y11dax__dabc_y=dabc_length_inv__dax
      ddabc1_y11dax__dbc_z=abc_length_inv
      ddabc1_y11dax__ddabc_length_inv11dax=abc_y
c dabc1_y__day= abc_y*dabc_length_inv__day
      ddabc1_y11day__dabc_y=dabc_length_inv__day
      ddabc1_y11day__ddabc_length_inv11day=abc_y
c dabc1_y__daz=-abc_length_inv*bc_x + abc_y*dabc_length_inv__daz
      ddabc1_y11daz__dabc_length_inv=-bc_x
      ddabc1_y11daz__dabc_y=dabc_length_inv__daz
      ddabc1_y11daz__dbc_x=-abc_length_inv
      ddabc1_y11daz__ddabc_length_inv11daz=abc_y
c dabc1_y__dbx=abc_length_inv*dabc_y__dbx + abc_y*dabc_length_inv__dbx
      ddabc1_y11dbx__dabc_length_inv=dabc_y__dbx
      ddabc1_y11dbx__dabc_y=dabc_length_inv__dbx
      ddabc1_y11dbx__ddabc_length_inv11dbx=abc_y
      ddabc1_y11dbx__ddabc_y11dbx=abc_length_inv
c dabc1_y__dby= abc_y*dabc_length_inv__dby
      ddabc1_y11dby__dabc_y=dabc_length_inv__dby
      ddabc1_y11dby__ddabc_length_inv11dby=abc_y
c dabc1_y__dbz=abc_length_inv*dabc_y__dbz + abc_y*dabc_length_inv__dbz
      ddabc1_y11dbz__dabc_length_inv=dabc_y__dbz
      ddabc1_y11dbz__dabc_y=dabc_length_inv__dbz
      ddabc1_y11dbz__ddabc_length_inv11dbz=abc_y
      ddabc1_y11dbz__ddabc_y11dbz=abc_length_inv
c dabc1_y__dcx=abc_length_inv*ab_z + abc_y*dabc_length_inv__dcx
      ddabc1_y11dcx__dab_z=abc_length_inv
      ddabc1_y11dcx__dabc_length_inv=ab_z
      ddabc1_y11dcx__dabc_y=dabc_length_inv__dcx
      ddabc1_y11dcx__ddabc_length_inv11dcx=abc_y
c dabc1_y__dcy= abc_y*dabc_length_inv__dcy
      ddabc1_y11dcy__dabc_y=dabc_length_inv__dcy
      ddabc1_y11dcy__ddabc_length_inv11dcy=abc_y
c dabc1_y__dcz=-abc_length_inv*ab_x + abc_y*dabc_length_inv__dcz
      ddabc1_y11dcz__dab_x=-abc_length_inv
      ddabc1_y11dcz__dabc_length_inv=-ab_x
      ddabc1_y11dcz__dabc_y=dabc_length_inv__dcz
      ddabc1_y11dcz__ddabc_length_inv11dcz=abc_y

c dabc1_z__dax=-abc_length_inv*bc_y + abc_z*dabc_length_inv__dax
      ddabc1_z11dax__dabc_length_inv=-bc_y
      ddabc1_z11dax__dabc_z=dabc_length_inv__dax
      ddabc1_z11dax__dbc_y=-abc_length_inv
      ddabc1_z11dax__ddabc_length_inv11dax=abc_z
c dabc1_z__day=abc_length_inv*bc_x + abc_z*dabc_length_inv__day
      ddabc1_z11day__dabc_length_inv=bc_x
      ddabc1_z11day__dabc_z=dabc_length_inv__day
      ddabc1_z11day__dbc_x=abc_length_inv
      ddabc1_z11day__ddabc_length_inv11day=abc_z
c dabc1_z__daz= abc_z*dabc_length_inv__daz
      ddabc1_z11daz__dabc_z=dabc_length_inv__daz
      ddabc1_z11daz__ddabc_length_inv11daz=abc_z
c dabc1_z__dbx=abc_length_inv*dabc_z__dbx + abc_z*dabc_length_inv__dbx
      ddabc1_z11dbx__dabc_length_inv=dabc_z__dbx
      ddabc1_z11dbx__dabc_z=dabc_length_inv__dbx
      ddabc1_z11dbx__ddabc_length_inv11dbx=abc_z
      ddabc1_z11dbx__ddabc_z11dbx=abc_length_inv
c dabc1_z__dby=abc_length_inv*dabc_z__dby + abc_z*dabc_length_inv__dby
      ddabc1_z11dby__dabc_length_inv=dabc_z__dby
      ddabc1_z11dby__dabc_z=dabc_length_inv__dby
      ddabc1_z11dby__ddabc_length_inv11dby=abc_z
      ddabc1_z11dby__ddabc_z11dby=abc_length_inv
c dabc1_z__dbz= abc_z*dabc_length_inv__dbz
      ddabc1_z11dbz__dabc_z=dabc_length_inv__dbz
      ddabc1_z11dbz__ddabc_length_inv11dbz=abc_z
c dabc1_z__dcx=-abc_length_inv*ab_y + abc_z*dabc_length_inv__dcx
      ddabc1_z11dcx__dab_y=-abc_length_inv
      ddabc1_z11dcx__dabc_length_inv=-ab_y
      ddabc1_z11dcx__dabc_z=dabc_length_inv__dcx
      ddabc1_z11dcx__ddabc_length_inv11dcx=abc_z
c dabc1_z__dcy=abc_length_inv*ab_x + abc_z*dabc_length_inv__dcy
      ddabc1_z11dcy__dab_x=abc_length_inv
      ddabc1_z11dcy__dabc_length_inv=ab_x
      ddabc1_z11dcy__dabc_z=dabc_length_inv__dcy
      ddabc1_z11dcy__ddabc_length_inv11dcy=abc_z
c dabc1_z__dcz= abc_z*dabc_length_inv__dcz
      ddabc1_z11dcz__dabc_z=dabc_length_inv__dcz
      ddabc1_z11dcz__ddabc_length_inv11dcz=abc_z



c dbc1_x__dbx= bc_length_inv + bc_x*dbc_length_inv__dbc_x
      ddbc1_x11dbx__dbx=
     1  ddbc1_x11dbx__dbc_length_inv*dbc_length_inv__dbx
     1 +ddbc1_x11dbx__dbc_x*dbc_x__dbx
     1 +ddbc1_x11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dbx
      ddbc1_x11dbx__dby=
     1  ddbc1_x11dbx__dbc_length_inv*dbc_length_inv__dby
     1 +ddbc1_x11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dby
      ddbc1_x11dbx__dbz=
     1  ddbc1_x11dbx__dbc_length_inv*dbc_length_inv__dbz
     1 +ddbc1_x11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dbz
      ddbc1_x11dbx__dcx=
     1  ddbc1_x11dbx__dbc_length_inv*dbc_length_inv__dcx
     1 +ddbc1_x11dbx__dbc_x*dbc_x__dcx
     1 +ddbc1_x11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcx
      ddbc1_x11dbx__dcy=
     1  ddbc1_x11dbx__dbc_length_inv*dbc_length_inv__dcy
     1 +ddbc1_x11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcy
      ddbc1_x11dbx__dcz=
     1  ddbc1_x11dbx__dbc_length_inv*dbc_length_inv__dcz
     1 +ddbc1_x11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcz
c dbc1_x__dby= bc_x*dbc_length_inv__dbc_y
      ddbc1_x11dby__dby=
     1 ddbc1_x11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dby
      ddbc1_x11dby__dbz=
     1 ddbc1_x11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dbz
      ddbc1_x11dby__dcx=
     1  ddbc1_x11dby__dbc_x*dbc_x__dcx
     1 +ddbc1_x11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcx
      ddbc1_x11dby__dcy=
     1 ddbc1_x11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcy
      ddbc1_x11dby__dcz=
     1 ddbc1_x11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcz
c dbc1_x__dbz= bc_x*dbc_length_inv__dbc_z
      ddbc1_x11dbz__dbz=
     1 ddbc1_x11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dbz
      ddbc1_x11dbz__dcx=
     1  ddbc1_x11dbz__dbc_x*dbc_x__dcx
     1 +ddbc1_x11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcx
      ddbc1_x11dbz__dcy=
     1 ddbc1_x11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcy
      ddbc1_x11dbz__dcz=
     1 ddbc1_x11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcz
c dbc1_x__dcx= -bc_length_inv - bc_x*dbc_length_inv__dbc_x
      ddbc1_x11dcx__dcx=
     1  ddbc1_x11dcx__dbc_length_inv*dbc_length_inv__dcx
     1 +ddbc1_x11dcx__dbc_x*dbc_x__dcx
     1 +ddbc1_x11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcx
      ddbc1_x11dcx__dcy=
     1  ddbc1_x11dcx__dbc_length_inv*dbc_length_inv__dcy
     1 +ddbc1_x11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcy
      ddbc1_x11dcx__dcz=
     1  ddbc1_x11dcx__dbc_length_inv*dbc_length_inv__dcz
     1 +ddbc1_x11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcz
c dbc1_x__dcy= -bc_x*dbc_length_inv__dbc_y
      ddbc1_x11dcy__dcy=
     1 ddbc1_x11dcy__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcy
      ddbc1_x11dcy__dcz=
     1 ddbc1_x11dcy__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcz
c dbc1_x__dcz= -bc_x*dbc_length_inv__dbc_z
      ddbc1_x11dcz__dcz=
     1 ddbc1_x11dcz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcz

c dbc1_y__dbx= bc_y*dbc_length_inv__dbc_x
      ddbc1_y11dbx__dbx=
     1 ddbc1_y11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dbx
      ddbc1_y11dbx__dby=
     1  ddbc1_y11dbx__dbc_y*dbc_y__dby
     1 +ddbc1_y11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dby
      ddbc1_y11dbx__dbz=
     1 ddbc1_y11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dbz
      ddbc1_y11dbx__dcx=
     1 ddbc1_y11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcx
      ddbc1_y11dbx__dcy=
     1  ddbc1_y11dbx__dbc_y*dbc_y__dcy
     1 +ddbc1_y11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcy
      ddbc1_y11dbx__dcz=
     1 ddbc1_y11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcz
c dbc1_y__dby= bc_length_inv + bc_y*dbc_length_inv__dbc_y
      ddbc1_y11dby__dby=
     1  ddbc1_y11dby__dbc_length_inv*dbc_length_inv__dby
     1 +ddbc1_y11dby__dbc_y*dbc_y__dby
     1 +ddbc1_y11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dby
      ddbc1_y11dby__dbz=
     1  ddbc1_y11dby__dbc_length_inv*dbc_length_inv__dbz
     1 +ddbc1_y11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dbz
      ddbc1_y11dby__dcx=
     1  ddbc1_y11dby__dbc_length_inv*dbc_length_inv__dcx
     1 +ddbc1_y11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcx
      ddbc1_y11dby__dcy=
     1  ddbc1_y11dby__dbc_length_inv*dbc_length_inv__dcy
     1 +ddbc1_y11dby__dbc_y*dbc_y__dcy
     1 +ddbc1_y11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcy
      ddbc1_y11dby__dcz=
     1  ddbc1_y11dby__dbc_length_inv*dbc_length_inv__dcz
     1 +ddbc1_y11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcz
c dbc1_y__dbz= bc_y*dbc_length_inv__dbc_z
      ddbc1_y11dbz__dbz=
     1 ddbc1_y11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dbz
      ddbc1_y11dbz__dcx=
     1 ddbc1_y11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcx
      ddbc1_y11dbz__dcy=
     1  ddbc1_y11dbz__dbc_y*dbc_y__dcy
     1 +ddbc1_y11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcy
      ddbc1_y11dbz__dcz=
     1 ddbc1_y11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcz
c dbc1_y__dcx= -bc_y*dbc_length_inv__dbc_x
      ddbc1_y11dcx__dcx=
     1 ddbc1_y11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcx
      ddbc1_y11dcx__dcy=
     1  ddbc1_y11dcx__dbc_y*dbc_y__dcy
     1 +ddbc1_y11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcy
      ddbc1_y11dcx__dcz=
     1 ddbc1_y11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcz
c dbc1_y__dcy= -bc_length_inv - bc_y*dbc_length_inv__dbc_y
      ddbc1_y11dcy__dcy=
     1  ddbc1_y11dcy__dbc_length_inv*dbc_length_inv__dcy
     1 +ddbc1_y11dcy__dbc_y*dbc_y__dcy
     1 +ddbc1_y11dcy__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcy
      ddbc1_y11dcy__dcz=
     1  ddbc1_y11dcy__dbc_length_inv*dbc_length_inv__dcz
     1 +ddbc1_y11dcy__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcz
c dbc1_y__dcz= -bc_y*dbc_length_inv__dbc_z
      ddbc1_y11dcz__dcz=
     1 ddbc1_y11dcz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcz

c dbc1_z__dbx= bc_z*dbc_length_inv__dbc_x
      ddbc1_z11dbx__dbx=
     1 ddbc1_z11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dbx
      ddbc1_z11dbx__dby=
     1 ddbc1_z11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dby
      ddbc1_z11dbx__dbz=
     1  ddbc1_z11dbx__dbc_z*dbc_z__dbz
     1 +ddbc1_z11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dbz
      ddbc1_z11dbx__dcx=
     1 ddbc1_z11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcx
      ddbc1_z11dbx__dcy=
     1 ddbc1_z11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcy
      ddbc1_z11dbx__dcz=
     1  ddbc1_z11dbx__dbc_z*dbc_z__dcz
     1 +ddbc1_z11dbx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcz
c dbc1_z__dby= bc_z*dbc_length_inv__dbc_y
      ddbc1_z11dby__dby=
     1 ddbc1_z11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dby
      ddbc1_z11dby__dbz=
     1  ddbc1_z11dby__dbc_z*dbc_z__dbz
     1 +ddbc1_z11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dbz
      ddbc1_z11dby__dcx=
     1 ddbc1_z11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcx
      ddbc1_z11dby__dcy=
     1 ddbc1_z11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcy
      ddbc1_z11dby__dcz=
     1  ddbc1_z11dby__dbc_z*dbc_z__dcz
     1 +ddbc1_z11dby__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcz
c dbc1_z__dbz= bc_length_inv + bc_z*dbc_length_inv__dbc_z
      ddbc1_z11dbz__dbz=
     1  ddbc1_z11dbz__dbc_length_inv*dbc_length_inv__dbz
     1 +ddbc1_z11dbz__dbc_z*dbc_z__dbz
     1 +ddbc1_z11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dbz
      ddbc1_z11dbz__dcx=
     1  ddbc1_z11dbz__dbc_length_inv*dbc_length_inv__dcx
     1 +ddbc1_z11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcx
      ddbc1_z11dbz__dcy=
     1  ddbc1_z11dbz__dbc_length_inv*dbc_length_inv__dcy
     1 +ddbc1_z11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcy
      ddbc1_z11dbz__dcz=
     1  ddbc1_z11dbz__dbc_length_inv*dbc_length_inv__dcz
     1 +ddbc1_z11dbz__dbc_z*dbc_z__dcz
     1 +ddbc1_z11dbz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcz
c dbc1_z__dcx= -bc_z*dbc_length_inv__dbc_x
      ddbc1_z11dcx__dcx=
     1 ddbc1_z11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcx
      ddbc1_z11dcx__dcy=
     1 ddbc1_z11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcy
      ddbc1_z11dcx__dcz=
     1  ddbc1_z11dcx__dbc_z*dbc_z__dcz
     1 +ddbc1_z11dcx__ddbc_length_inv11dbc_x*ddbc_length_inv11dbc_x__dcz
c dbc1_z__dcy= -bc_z*dbc_length_inv__dbc_y
      ddbc1_z11dcy__dcy=
     1 ddbc1_z11dcy__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcy
      ddbc1_z11dcy__dcz=
     1  ddbc1_z11dcy__dbc_z*dbc_z__dcz
     1 +ddbc1_z11dcy__ddbc_length_inv11dbc_y*ddbc_length_inv11dbc_y__dcz
c dbc1_z__dcz= -bc_length_inv - bc_z*dbc_length_inv__dbc_z
      ddbc1_z11dcz__dcz=
     1  ddbc1_z11dcz__dbc_length_inv*dbc_length_inv__dcz
     1 +ddbc1_z11dcz__dbc_z*dbc_z__dcz
     1 +ddbc1_z11dcz__ddbc_length_inv11dbc_z*ddbc_length_inv11dbc_z__dcz


c daux_x__dax= bc1_z*dabc1_y__dax - bc1_y*dabc1_z__dax
      ddaux_x11dax__dbc1_y=-dabc1_z__dax
      ddaux_x11dax__dbc1_z=dabc1_y__dax
      ddaux_x11dax__ddabc1_y11dax=bc1_z
      ddaux_x11dax__ddabc1_z11dax=-bc1_y
c daux_x__day= bc1_z*dabc1_y__day - bc1_y*dabc1_z__day
      ddaux_x11day__dbc1_y=-dabc1_z__day
      ddaux_x11day__dbc1_z=dabc1_y__day
      ddaux_x11day__ddabc1_y11dax=bc1_z
      ddaux_x11day__ddabc1_z11dax=-bc1_y
c daux_x__daz= bc1_z*dabc1_y__daz - bc1_y*dabc1_z__daz
      ddaux_x11daz__dbc1_y=-dabc1_z__daz
      ddaux_x11daz__dbc1_z=dabc1_y__daz
      ddaux_x11daz__ddabc1_y11dax=bc1_z
      ddaux_x11daz__ddabc1_z11dax=-bc1_y
c daux_x__dbx= bc1_z*dabc1_y__dbx + abc1_y*dbc1_z__dbx - abc1_z*dbc1_y__dbx - bc1_y*dabc1_z__dbx
      ddaux_x11dbx__dabc1_y=dbc1_z__dbx
      ddaux_x11dbx__dabc1_z=-dbc1_y__dbx
      ddaux_x11dbx__dbc1_y=-dabc1_z__dbx
      ddaux_x11dbx__dbc1_z=dabc1_y__dbx
      ddaux_x11dbx__ddabc1_y11dbx=bc1_z
      ddaux_x11dbx__ddabc1_z11dbx=-bc1_y
      ddaux_x11dbx__ddbc1_y11dbx=-abc1_z
      ddaux_x11dbx__ddbc1_z11dbx=abc1_y
c daux_x__dby= bc1_z*dabc1_y__dby + abc1_y*dbc1_z__dby - abc1_z*dbc1_y__dby - bc1_y*dabc1_z__dby
      ddaux_x11dby__dabc1_y=dbc1_z__dby
      ddaux_x11dby__dabc1_z=-dbc1_y__dby
      ddaux_x11dby__dbc1_y=-dabc1_z__dby
      ddaux_x11dby__dbc1_z=dabc1_y__dby
      ddaux_x11dby__ddabc1_y11dby=bc1_z
      ddaux_x11dby__ddabc1_z11dby=-bc1_y
      ddaux_x11dby__ddbc1_y11dby=-abc1_z
      ddaux_x11dby__ddbc1_z11dby=abc1_y
c daux_x__dbz= bc1_z*dabc1_y__dbz + abc1_y*dbc1_z__dbz - abc1_z*dbc1_y__dbz - bc1_y*dabc1_z__dbz
      ddaux_x11dbz__dabc1_y=dbc1_z__dbz
      ddaux_x11dbz__dabc1_z=-dbc1_y__dbz
      ddaux_x11dbz__dbc1_y=-dabc1_z__dbz
      ddaux_x11dbz__dbc1_z=dabc1_y__dbz
      ddaux_x11dbz__ddabc1_y11dbz=bc1_z
      ddaux_x11dbz__ddabc1_z11dbz=-bc1_y
      ddaux_x11dbz__ddbc1_y11dbz=-abc1_z
      ddaux_x11dbz__ddbc1_z11dbz=abc1_y
c daux_x__dcx= bc1_z*dabc1_y__dcx + abc1_y*dbc1_z__dcx - abc1_z*dbc1_y__dcx - bc1_y*dabc1_z__dcx
      ddaux_x11dcx__dabc1_y=dbc1_z__dcx
      ddaux_x11dcx__dabc1_z=-dbc1_y__dcx
      ddaux_x11dcx__dbc1_y=-dabc1_z__dcx
      ddaux_x11dcx__dbc1_z=dabc1_y__dcx
      ddaux_x11dcx__ddabc1_y11dcx=bc1_z
      ddaux_x11dcx__ddabc1_z11dcx=-bc1_y
      ddaux_x11dcx__ddbc1_y11dcx=-abc1_z
      ddaux_x11dcx__ddbc1_z11dcx=abc1_y
c daux_x__dcy= bc1_z*dabc1_y__dcy + abc1_y*dbc1_z__dcy - abc1_z*dbc1_y__dcy - bc1_y*dabc1_z__dcy
      ddaux_x11dcy__dabc1_y=dbc1_z__dcy
      ddaux_x11dcy__dabc1_z=-dbc1_y__dcy
      ddaux_x11dcy__dbc1_y=-dabc1_z__dcy
      ddaux_x11dcy__dbc1_z=dabc1_y__dcy
      ddaux_x11dcy__ddabc1_y11dcy=bc1_z
      ddaux_x11dcy__ddabc1_z11dcy=-bc1_y
      ddaux_x11dcy__ddbc1_y11dcy=-abc1_z
      ddaux_x11dcy__ddbc1_z11dcy=abc1_y
c daux_x__dcz= bc1_z*dabc1_y__dcz + abc1_y*dbc1_z__dcz - abc1_z*dbc1_y__dcz - bc1_y*dabc1_z__dcz
      ddaux_x11dcz__dabc1_y=dbc1_z__dcz
      ddaux_x11dcz__dabc1_z=-dbc1_y__dcz
      ddaux_x11dcz__dbc1_y=-dabc1_z__dcz
      ddaux_x11dcz__dbc1_z=dabc1_y__dcz
      ddaux_x11dcz__ddabc1_y11dcz=bc1_z
      ddaux_x11dcz__ddabc1_z11dcz=-bc1_y
      ddaux_x11dcz__ddbc1_y11dcz=-abc1_z
      ddaux_x11dcz__ddbc1_z11dcz=abc1_y


c daux_y__dax= bc1_x*dabc1_z__dax - bc1_z*dabc1_x__dax
      ddaux_y11dax__dbc1_x=dabc1_z__dax
      ddaux_y11dax__dbc1_z=-dabc1_x__dax
      ddaux_y11dax__ddabc1_x11dax=-bc1_z
      ddaux_y11dax__ddabc1_z11dax=bc1_x
c daux_y__day= bc1_x*dabc1_z__day - bc1_z*dabc1_x__day
      ddaux_y11day__dbc1_x=dabc1_z__day
      ddaux_y11day__dbc1_z=-dabc1_x__day
      ddaux_y11day__ddabc1_x11day=-bc1_z
      ddaux_y11day__ddabc1_z11day=bc1_x
c daux_y__daz= bc1_x*dabc1_z__daz - bc1_z*dabc1_x__daz
      ddaux_y11daz__dbc1_x=dabc1_z__daz
      ddaux_y11daz__dbc1_z=-dabc1_x__daz
      ddaux_y11daz__ddabc1_x11daz=-bc1_z
      ddaux_y11daz__ddabc1_z11daz=bc1_x
c daux_y__dbx= bc1_x*dabc1_z__dbx + abc1_z*dbc1_x__dbx - abc1_x*dbc1_z__dbx - bc1_z*dabc1_x__dbx
      ddaux_y11dbx__dabc1_x=-dbc1_z__dbx
      ddaux_y11dbx__dabc1_z=dbc1_x__dbx
      ddaux_y11dbx__dbc1_x=dabc1_z__dbx
      ddaux_y11dbx__dbc1_z=-dabc1_x__dbx
      ddaux_y11dbx__ddabc1_x11dbx=-bc1_z
      ddaux_y11dbx__ddabc1_z11dbx=bc1_x
      ddaux_y11dbx__ddbc1_x11dbx=abc1_z
      ddaux_y11dbx__ddbc1_z11dbx=-abc1_x
c daux_y__dby= bc1_x*dabc1_z__dby + abc1_z*dbc1_x__dby - abc1_x*dbc1_z__dby - bc1_z*dabc1_x__dby
      ddaux_y11dby__dabc1_x=-dbc1_z__dby
      ddaux_y11dby__dabc1_z=dbc1_x__dby
      ddaux_y11dby__dbc1_x=dabc1_z__dby
      ddaux_y11dby__dbc1_z=-dabc1_x__dby
      ddaux_y11dby__ddabc1_x11dby=-bc1_z
      ddaux_y11dby__ddabc1_z11dby=bc1_x
      ddaux_y11dby__ddbc1_x11dby=abc1_z
      ddaux_y11dby__ddbc1_z11dby=-abc1_x
c daux_y__dbz= bc1_x*dabc1_z__dbz + abc1_z*dbc1_x__dbz - abc1_x*dbc1_z__dbz - bc1_z*dabc1_x__dbz
      ddaux_y11dbz__dabc1_x=-dbc1_z__dbz
      ddaux_y11dbz__dabc1_z=dbc1_x__dbz
      ddaux_y11dbz__dbc1_x=dabc1_z__dbz
      ddaux_y11dbz__dbc1_z=-dabc1_x__dbz
      ddaux_y11dbz__ddabc1_x11dbz=-bc1_z
      ddaux_y11dbz__ddabc1_z11dbz=bc1_x
      ddaux_y11dbz__ddbc1_x11dbz=abc1_z
      ddaux_y11dbz__ddbc1_z11dbz=-abc1_x
c daux_y__dcx= bc1_x*dabc1_z__dcx + abc1_z*dbc1_x__dcx - abc1_x*dbc1_z__dcx - bc1_z*dabc1_x__dcx
      ddaux_y11dcx__dabc1_x=-dbc1_z__dcx
      ddaux_y11dcx__dabc1_z=dbc1_x__dcx
      ddaux_y11dcx__dbc1_x=dabc1_z__dcx
      ddaux_y11dcx__dbc1_z=-dabc1_x__dcx
      ddaux_y11dcx__ddabc1_x11dcx=-bc1_z
      ddaux_y11dcx__ddabc1_z11dcx=bc1_x
      ddaux_y11dcx__ddbc1_x11dcx=abc1_z
      ddaux_y11dcx__ddbc1_z11dcx=-abc1_x
c daux_y__dcy= bc1_x*dabc1_z__dcy + abc1_z*dbc1_x__dcy - abc1_x*dbc1_z__dcy - bc1_z*dabc1_x__dcy
      ddaux_y11dcy__dabc1_x=-dbc1_z__dcy
      ddaux_y11dcy__dabc1_z=dbc1_x__dcy
      ddaux_y11dcy__dbc1_x=dabc1_z__dcy
      ddaux_y11dcy__dbc1_z=-dabc1_x__dcy
      ddaux_y11dcy__ddabc1_x11dcy=-bc1_z
      ddaux_y11dcy__ddabc1_z11dcy=bc1_x
      ddaux_y11dcy__ddbc1_x11dcy=abc1_z
      ddaux_y11dcy__ddbc1_z11dcy=-abc1_x
c daux_y__dcz= bc1_x*dabc1_z__dcz + abc1_z*dbc1_x__dcz - abc1_x*dbc1_z__dcz - bc1_z*dabc1_x__dcz
      ddaux_y11dcz__dabc1_x=-dbc1_z__dcz
      ddaux_y11dcz__dabc1_z=dbc1_x__dcz
      ddaux_y11dcz__dbc1_x=dabc1_z__dcz
      ddaux_y11dcz__dbc1_z=-dabc1_x__dcz
      ddaux_y11dcz__ddabc1_x11dcz=-bc1_z
      ddaux_y11dcz__ddabc1_z11dcz=bc1_x
      ddaux_y11dcz__ddbc1_x11dcz=abc1_z
      ddaux_y11dcz__ddbc1_z11dcz=-abc1_x


c daux_z__dax= bc1_y*dabc1_x__dax - bc1_x*dabc1_y__dax
      ddaux_z11dax__dbc1_x=-dabc1_y__dax
      ddaux_z11dax__dbc1_y=dabc1_x__dax
      ddaux_z11dax__ddabc1_x11dax=bc1_y
      ddaux_z11dax__ddabc1_y11dax=-bc1_x
c daux_z__day= bc1_y*dabc1_x__day - bc1_x*dabc1_y__day
      ddaux_z11day__dbc1_x=-dabc1_y__day
      ddaux_z11day__dbc1_y=dabc1_x__day
      ddaux_z11day__ddabc1_x11day=bc1_y
      ddaux_z11day__ddabc1_y11day=-bc1_x
c daux_z__daz= bc1_y*dabc1_x__daz - bc1_x*dabc1_y__daz
      ddaux_z11daz__dbc1_x=-dabc1_y__daz
      ddaux_z11daz__dbc1_y=dabc1_x__daz
      ddaux_z11daz__ddabc1_x11daz=bc1_y
      ddaux_z11daz__ddabc1_y11daz=-bc1_x
c daux_z__dbx= bc1_y*dabc1_x__dbx + abc1_x*dbc1_y__dbx - abc1_y*dbc1_x__dbx - bc1_x*dabc1_y__dbx
      ddaux_z11dbx__dabc1_x=dbc1_y__dbx
      ddaux_z11dbx__dabc1_y=-dbc1_x__dbx
      ddaux_z11dbx__dbc1_x=-dabc1_y__dbx
      ddaux_z11dbx__dbc1_y=dabc1_x__dbx
      ddaux_z11dbx__ddabc1_x11dbx=bc1_y
      ddaux_z11dbx__ddabc1_y11dbx=-bc1_x
      ddaux_z11dbx__ddbc1_x11dbx=-abc1_y
      ddaux_z11dbx__ddbc1_y11dbx=abc1_x
c daux_z__dby= bc1_y*dabc1_x__dby + abc1_x*dbc1_y__dby - abc1_y*dbc1_x__dby - bc1_x*dabc1_y__dby
      ddaux_z11dby__dabc1_x=dbc1_y__dby
      ddaux_z11dby__dabc1_y=-dbc1_x__dby
      ddaux_z11dby__dbc1_x=-dabc1_y__dby
      ddaux_z11dby__dbc1_y=dabc1_x__dby
      ddaux_z11dby__ddabc1_x11dby=bc1_y
      ddaux_z11dby__ddabc1_y11dby=-bc1_x
      ddaux_z11dby__ddbc1_x11dby=-abc1_y
      ddaux_z11dby__ddbc1_y11dby=abc1_x
c daux_z__dbz= bc1_y*dabc1_x__dbz + abc1_x*dbc1_y__dbz - abc1_y*dbc1_x__dbz - bc1_x*dabc1_y__dbz
      ddaux_z11dbz__dabc1_x=dbc1_y__dbz
      ddaux_z11dbz__dabc1_y=-dbc1_x__dbz
      ddaux_z11dbz__dbc1_x=-dabc1_y__dbz
      ddaux_z11dbz__dbc1_y=dabc1_x__dbz
      ddaux_z11dbz__ddabc1_x11dbz=bc1_y
      ddaux_z11dbz__ddabc1_y11dbz=-bc1_x
      ddaux_z11dbz__ddbc1_x11dbz=-abc1_y
      ddaux_z11dbz__ddbc1_y11dbz=abc1_x
c daux_z__dcx= bc1_y*dabc1_x__dcx + abc1_x*dbc1_y__dcx - abc1_y*dbc1_x__dcx - bc1_x*dabc1_y__dcx
      ddaux_z11dcx__dabc1_x=dbc1_y__dcx
      ddaux_z11dcx__dabc1_y=-dbc1_x__dcx
      ddaux_z11dcx__dbc1_x=-dabc1_y__dcx
      ddaux_z11dcx__dbc1_y=dabc1_x__dcx
      ddaux_z11dcx__ddabc1_x11dcx=bc1_y
      ddaux_z11dcx__ddabc1_y11dcx=-bc1_x
      ddaux_z11dcx__ddbc1_x11dcx=-abc1_y
      ddaux_z11dcx__ddbc1_y11dcx=abc1_x
c daux_z__dcy= bc1_y*dabc1_x__dcy + abc1_x*dbc1_y__dcy - abc1_y*dbc1_x__dcy - bc1_x*dabc1_y__dcy
      ddaux_z11dcy__dabc1_x=dbc1_y__dcy
      ddaux_z11dcy__dabc1_y=-dbc1_x__dcy
      ddaux_z11dcy__dbc1_x=-dabc1_y__dcy
      ddaux_z11dcy__dbc1_y=dabc1_x__dcy
      ddaux_z11dcy__ddabc1_x11dcy=bc1_y
      ddaux_z11dcy__ddabc1_y11dcy=-bc1_x
      ddaux_z11dcy__ddbc1_x11dcy=-abc1_y
      ddaux_z11dcy__ddbc1_y11dcy=abc1_x
c daux_z__dcz= bc1_y*dabc1_x__dcz + abc1_x*dbc1_y__dcz - abc1_y*dbc1_x__dcz - bc1_x*dabc1_y__dcz
      ddaux_z11dcz__dabc1_x=dbc1_y__dcz
      ddaux_z11dcz__dabc1_y=-dbc1_x__dcz
      ddaux_z11dcz__dbc1_x=-dabc1_y__dcz
      ddaux_z11dcz__dbc1_y=dabc1_x__dcz
      ddaux_z11dcz__ddabc1_x11dcz=bc1_y
      ddaux_z11dcz__ddabc1_y11dcz=-bc1_x
      ddaux_z11dcz__ddbc1_x11dcz=-abc1_y
      ddaux_z11dcz__ddbc1_y11dcz=abc1_x


c dbcd_length_inv__dbx=dbcd_length_inv__dbcd_y*cd_z - dbcd_length_inv__dbcd_z*cd_y
      ddbcd_length_inv11dbx__dbx=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dbx
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dbx
      ddbcd_length_inv11dbx__dby=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dby
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dby
      ddbcd_length_inv11dbx__dbz=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dbz
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dbz
      ddbcd_length_inv11dbx__dcx=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcx
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcx
      ddbcd_length_inv11dbx__dcy=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcy
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcy
     1 + ddbcd_length_inv11dbx__dcd_y*dcd_y__dcy
      ddbcd_length_inv11dbx__dcz=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcz
     1 + ddbcd_length_inv11dbx__dcd_z*dcd_z__dcz
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcz
      ddbcd_length_inv11dbx__ddx=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddx
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddx
      ddbcd_length_inv11dbx__ddy=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddy
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddy
     1 + ddbcd_length_inv11dbx__dcd_y*dcd_y__ddy
      ddbcd_length_inv11dbx__ddz=
     1   ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddz
     1 + ddbcd_length_inv11dbx__dcd_z*dcd_z__ddz
     1 + ddbcd_length_inv11dbx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddz
c dbcd_length_inv__dby=-dbcd_length_inv__dbcd_x*cd_z + dbcd_length_inv__dbcd_z*cd_x
      ddbcd_length_inv11dby__dby=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dby
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dby
      ddbcd_length_inv11dby__dbz=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dbz
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dbz
      ddbcd_length_inv11dby__dcx=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcx
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcx
     1 + ddbcd_length_inv11dby__dcd_x*dcd_x__dcx
      ddbcd_length_inv11dby__dcy=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcy
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcy
      ddbcd_length_inv11dby__dcz=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcz
     1 + ddbcd_length_inv11dby__dcd_z*dcd_z__dcz
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcz
      ddbcd_length_inv11dby__ddx=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddx
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddx
     1 + ddbcd_length_inv11dby__dcd_x*dcd_x__ddx
      ddbcd_length_inv11dby__ddy=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddy
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddy
      ddbcd_length_inv11dby__ddz=
     1   ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddz
     1 + ddbcd_length_inv11dby__dcd_z*dcd_z__ddz
     1 + ddbcd_length_inv11dby__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddz
c dbcd_length_inv__dbz=dbcd_length_inv__dbcd_x*cd_y - dbcd_length_inv__dbcd_y*cd_x
      ddbcd_length_inv11dbz__dbz=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dbz
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dbz
      ddbcd_length_inv11dbz__dcx=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcx
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcx
     1 + ddbcd_length_inv11dbz__dcd_x*dcd_x__dcx
      ddbcd_length_inv11dbz__dcy=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcy
     1 + ddbcd_length_inv11dbz__dcd_y*dcd_y__dcy
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcy
      ddbcd_length_inv11dbz__dcz=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcz
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcz
      ddbcd_length_inv11dbz__ddx=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddx
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddx
     1 + ddbcd_length_inv11dbz__dcd_x*dcd_x__ddx
      ddbcd_length_inv11dbz__ddy=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddy
     1 + ddbcd_length_inv11dbz__dcd_y*dcd_y__ddy
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddy
      ddbcd_length_inv11dbz__ddz=
     1   ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddz
     1 + ddbcd_length_inv11dbz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddz
c dbcd_length_inv__dcx=dbcd_length_inv__dbcd_y*dbcd_y__dcx + dbcd_length_inv__dbcd_z*dbcd_z__dcx
      ddbcd_length_inv11dcx__dcx=
     1   ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcx
     1 + ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcx
      ddbcd_length_inv11dcx__dcy=
     1   ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcy
     1 + ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcy
     1 + ddbcd_length_inv11dcx__ddbcd_z11dcx*
     1                                   ddbcd_z11dcx__dcy
      ddbcd_length_inv11dcx__dcz=
     1   ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcz
     1 + ddbcd_length_inv11dcx__ddbcd_y11dcx*
     1                                   ddbcd_y11dcx__dcz
     1 + ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcz
      ddbcd_length_inv11dcx__ddx=
     1   ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddx
     1 + ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddx
      ddbcd_length_inv11dcx__ddy=
     1   ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddy
     1 + ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddy
     1 + ddbcd_length_inv11dcx__ddbcd_z11dcx*
     1                                   ddbcd_z11dcx__ddy
      ddbcd_length_inv11dcx__ddz=
     1   ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddz
     1 + ddbcd_length_inv11dcx__ddbcd_y11dcx*
     1                                   ddbcd_y11dcx__ddz
     1 + ddbcd_length_inv11dcx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddz
c dbcd_length_inv__dcy=dbcd_length_inv__dbcd_x*dbcd_x__dcy + dbcd_length_inv__dbcd_z*dbcd_z__dcy
      ddbcd_length_inv11dcy__dcy=
     1   ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcy
     1 + ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcy
      ddbcd_length_inv11dcy__dcz=
     1   ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcz
     1 + ddbcd_length_inv11dcy__ddbcd_x11dcy*
     1                                   ddbcd_x11dcy__dcz
     1 + ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__dcz
      ddbcd_length_inv11dcy__ddx=
     1   ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddx
     1 + ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddx
     1 + ddbcd_length_inv11dcy__ddbcd_z11dcy*
     1                                   ddbcd_z11dcy__ddx
c     1 ddbcd_length_inv11dcy__ddx
      ddbcd_length_inv11dcy__ddy=
     1   ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddy
     1 + ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddy
      ddbcd_length_inv11dcy__ddz=
     1   ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddz
     1 + ddbcd_length_inv11dcy__ddbcd_x11dcy*
     1                                   ddbcd_x11dcy__ddz
     1 + ddbcd_length_inv11dcy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddz
c dbcd_length_inv__dcz=dbcd_length_inv__dbcd_x*dbcd_x__dcz + dbcd_length_inv__dbcd_y*dbcd_y__dcz
      ddbcd_length_inv11dcz__dcz=
     1   ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__dcz
     1 + ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__dcz
      ddbcd_length_inv11dcz__ddx=
     1   ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddx
     1 + ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddx
     1 + ddbcd_length_inv11dcz__ddbcd_y11dcz*
     1                                   ddbcd_y11dcz__ddx
      ddbcd_length_inv11dcz__ddy=
     1   ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddy
     1 + ddbcd_length_inv11dcz__ddbcd_x11dcz*
     1                                   ddbcd_x11dcz__ddy
     1 + ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddy
      ddbcd_length_inv11dcz__ddz=
     1   ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddz
     1 + ddbcd_length_inv11dcz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddz
c dbcd_length_inv__ddx=dbcd_length_inv__dbcd_y*bc_z - dbcd_length_inv__dbcd_z*bc_y
      ddbcd_length_inv11ddx__ddx=
     1   ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddx
     1 + ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddx
      ddbcd_length_inv11ddx__ddy=
     1   ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddy
     1 + ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddy
      ddbcd_length_inv11ddx__ddz=
     1   ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddz
     1 + ddbcd_length_inv11ddx__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddz
c dbcd_length_inv__ddy=-dbcd_length_inv__dbcd_x*bc_z + dbcd_length_inv__dbcd_z*bc_x
      ddbcd_length_inv11ddy__ddy=
     1   ddbcd_length_inv11ddy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddy
     1 + ddbcd_length_inv11ddy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddy
      ddbcd_length_inv11ddy__ddz=
     1   ddbcd_length_inv11ddy__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddz
     1 + ddbcd_length_inv11ddy__ddbcd_length_inv11dbcd_z*
     1                                   ddbcd_length_inv11dbcd_z__ddz
c dbcd_length_inv__ddz=dbcd_length_inv__dbcd_x*bc_y - dbcd_length_inv__dbcd_y*bc_x
      ddbcd_length_inv11ddz__ddz=
     1   ddbcd_length_inv11ddz__ddbcd_length_inv11dbcd_x*
     1                                   ddbcd_length_inv11dbcd_x__ddz
     1 + ddbcd_length_inv11ddz__ddbcd_length_inv11dbcd_y*
     1                                   ddbcd_length_inv11dbcd_y__ddz


c dbcd1_x__dbx= bcd_x*dbcd_length_inv__dbx
      ddbcd1_x11dbx__dbcd_x=dbcd_length_inv__dbx
      ddbcd1_x11dbx__ddbcd_length_inv11dbx=bcd_x
c dbcd1_x__dby=-bcd_length_inv*cd_z + bcd_x*dbcd_length_inv__dby
      ddbcd1_x11dby__dbcd_length_inv=-cd_z
      ddbcd1_x11dby__dbcd_x=dbcd_length_inv__dby
      ddbcd1_x11dby__dcd_z=-bcd_length_inv
      ddbcd1_x11dby__ddbcd_length_inv11dby=bcd_x
c dbcd1_x__dbz=bcd_length_inv*cd_y + bcd_x*dbcd_length_inv__dbz
      ddbcd1_x11dbz__dbcd_length_inv=cd_y
      ddbcd1_x11dbz__dbcd_x=dbcd_length_inv__dbz
      ddbcd1_x11dbz__dcd_y=bcd_length_inv
      ddbcd1_x11dbz__ddbcd_length_inv11dbz=bcd_x
c dbcd1_x__dcx= bcd_x*dbcd_length_inv__dcx
      ddbcd1_x11dcx__dbcd_x=dbcd_length_inv__dcx
      ddbcd1_x11dcx__ddbcd_length_inv11dcx=bcd_x
c dbcd1_x__dcy=bcd_length_inv*dbcd_x__dcy + bcd_x*dbcd_length_inv__dcy
      ddbcd1_x11dcy__dbcd_length_inv=dbcd_x__dcy
      ddbcd1_x11dcy__dbcd_x=dbcd_length_inv__dcy
      ddbcd1_x11dcy__ddbcd_length_inv11dcy=bcd_x
      ddbcd1_x11dcy__ddbcd_x11dcy=bcd_length_inv
c dbcd1_x__dcz=bcd_length_inv*dbcd_x__dcz + bcd_x*dbcd_length_inv__dcz
      ddbcd1_x11dcz__dbcd_length_inv=dbcd_x__dcz
      ddbcd1_x11dcz__dbcd_x=dbcd_length_inv__dcz
      ddbcd1_x11dcz__ddbcd_length_inv11dcz=bcd_x
      ddbcd1_x11dcz__ddbcd_x11dcz=bcd_length_inv
c dbcd1_x__ddx= bcd_x*dbcd_length_inv__ddx
      ddbcd1_x11ddx__dbcd_x=dbcd_length_inv__ddx
      ddbcd1_x11ddx__ddbcd_length_inv11ddx=bcd_x
c dbcd1_x__ddy=-bcd_length_inv*bc_z + bcd_x*dbcd_length_inv__ddy
      ddbcd1_x11ddy__dbc_z=-bcd_length_inv
      ddbcd1_x11ddy__dbcd_length_inv=-bc_z
      ddbcd1_x11ddy__dbcd_x=dbcd_length_inv__ddy
      ddbcd1_x11ddy__ddbcd_length_inv11ddy=bcd_x
c dbcd1_x__ddz=bcd_length_inv*bc_y + bcd_x*dbcd_length_inv__ddz
      ddbcd1_x11ddz__dbc_y=bcd_length_inv
      ddbcd1_x11ddz__dbcd_length_inv=bc_y
      ddbcd1_x11ddz__dbcd_x=dbcd_length_inv__ddz
      ddbcd1_x11ddz__ddbcd_length_inv11ddz=bcd_x

c dbcd1_y__dbx=bcd_length_inv*cd_z + bcd_y*dbcd_length_inv__dbx
      ddbcd1_y11dbx__dbcd_length_inv=cd_z
      ddbcd1_y11dbx__dbcd_y=dbcd_length_inv__dbx
      ddbcd1_y11dbx__dcd_z=bcd_length_inv
      ddbcd1_y11dbx__ddbcd_length_inv11dbx=bcd_y
c dbcd1_y__dby= bcd_y*dbcd_length_inv__dby
      ddbcd1_y11dby__dbcd_y=dbcd_length_inv__dby
      ddbcd1_y11dby__ddbcd_length_inv11dby=bcd_y
c dbcd1_y__dbz=-bcd_length_inv*cd_x + bcd_y*dbcd_length_inv__dbz
      ddbcd1_y11dbz__dbcd_length_inv=-cd_x
      ddbcd1_y11dbz__dbcd_y=dbcd_length_inv__dbz
      ddbcd1_y11dbz__dcd_x=-bcd_length_inv
      ddbcd1_y11dbz__ddbcd_length_inv11dbz=bcd_y
c dbcd1_y__dcx=bcd_length_inv*dbcd_y__dcx + bcd_y*dbcd_length_inv__dcx
      ddbcd1_y11dcx__dbcd_length_inv=dbcd_y__dcx
      ddbcd1_y11dcx__dbcd_y=dbcd_length_inv__dcx
      ddbcd1_y11dcx__ddbcd_length_inv11dcx=bcd_y
      ddbcd1_y11dcx__ddbcd_y11dcx=bcd_length_inv
c dbcd1_y__dcy= bcd_y*dbcd_length_inv__dcy
      ddbcd1_y11dcy__dbcd_y=dbcd_length_inv__dcy
      ddbcd1_y11dcy__ddbcd_length_inv11dcy=bcd_y
c dbcd1_y__dcz=bcd_length_inv*dbcd_y__dcz + bcd_y*dbcd_length_inv__dcz
      ddbcd1_y11dcz__dbcd_length_inv=dbcd_y__dcz
      ddbcd1_y11dcz__dbcd_y=dbcd_length_inv__dcz
      ddbcd1_y11dcz__ddbcd_length_inv11dcz=bcd_y
      ddbcd1_y11dcz__ddbcd_y11dcz=bcd_length_inv
c dbcd1_y__ddx=bcd_length_inv*bc_z + bcd_y*dbcd_length_inv__ddx
      ddbcd1_y11ddx__dbc_z=bcd_length_inv
      ddbcd1_y11ddx__dbcd_length_inv=bc_z
      ddbcd1_y11ddx__dbcd_y=dbcd_length_inv__ddx
      ddbcd1_y11ddx__ddbcd_length_inv11ddx=bcd_y
c dbcd1_y__ddy= bcd_y*dbcd_length_inv__ddy
      ddbcd1_y11ddy__dbcd_y=dbcd_length_inv__ddy
      ddbcd1_y11ddy__ddbcd_length_inv11ddy= bcd_y
c dbcd1_y__ddz=-bcd_length_inv*bc_x + bcd_y*dbcd_length_inv__ddz
      ddbcd1_y11ddz__dbc_x=-bcd_length_inv
      ddbcd1_y11ddz__dbcd_length_inv=-bc_x
      ddbcd1_y11ddz__dbcd_y=dbcd_length_inv__ddz
      ddbcd1_y11ddz__ddbcd_length_inv11ddz=bcd_y

c dbcd1_z__dbx=-bcd_length_inv*cd_y + bcd_z*dbcd_length_inv__dbx
      ddbcd1_z11dbx__dbcd_length_inv=-cd_y
      ddbcd1_z11dbx__dbcd_z=dbcd_length_inv__dbx
      ddbcd1_z11dbx__dcd_y=-bcd_length_inv
      ddbcd1_z11dbx__ddbcd_length_inv11dbx=bcd_z
c dbcd1_z__dby=bcd_length_inv*cd_x + bcd_z*dbcd_length_inv__dby
      ddbcd1_z11dby__dbcd_length_inv=cd_x
      ddbcd1_z11dby__dbcd_z=dbcd_length_inv__dby
      ddbcd1_z11dby__dcd_x=bcd_length_inv
      ddbcd1_z11dby__ddbcd_length_inv11dby=bcd_z
c dbcd1_z__dbz= bcd_z*dbcd_length_inv__dbz
      ddbcd1_z11dbz__dbcd_z=dbcd_length_inv__dbz
      ddbcd1_z11dbz__ddbcd_length_inv11dbz=bcd_z
c dbcd1_z__dcx=bcd_length_inv*dbcd_z__dcx + bcd_z*dbcd_length_inv__dcx
      ddbcd1_z11dcx__dbcd_length_inv=dbcd_z__dcx
      ddbcd1_z11dcx__dbcd_z=dbcd_length_inv__dcx
      ddbcd1_z11dcx__ddbcd_length_inv11dcx=bcd_z
      ddbcd1_z11dcx__ddbcd_z11dcx=bcd_length_inv
c dbcd1_z__dcy=bcd_length_inv*dbcd_z__dcy + bcd_z*dbcd_length_inv__dcy
      ddbcd1_z11dcy__dbcd_length_inv=dbcd_z__dcy
      ddbcd1_z11dcy__dbcd_z=dbcd_length_inv__dcy
      ddbcd1_z11dcy__ddbcd_length_inv11dcy=bcd_z
      ddbcd1_z11dcy__ddbcd_z11dcy=bcd_length_inv
c dbcd1_z__dcz= bcd_z*dbcd_length_inv__dcz
      ddbcd1_z11dcz__dbcd_z=dbcd_length_inv__dcz
      ddbcd1_z11dcz__ddbcd_length_inv11dcz=bcd_z
c dbcd1_z__ddx=-bcd_length_inv*bc_y + bcd_z*dbcd_length_inv__ddx
      ddbcd1_z11ddx__dbc_y=-bcd_length_inv
      ddbcd1_z11ddx__dbcd_length_inv=-bc_y
      ddbcd1_z11ddx__dbcd_z=dbcd_length_inv__ddx
      ddbcd1_z11ddx__ddbcd_length_inv11ddx=bcd_z
c dbcd1_z__ddy=bcd_length_inv*bc_x + bcd_z*dbcd_length_inv__ddy
      ddbcd1_z11ddy__dbc_x=bcd_length_inv
      ddbcd1_z11ddy__dbcd_length_inv=bc_x
      ddbcd1_z11ddy__dbcd_z=dbcd_length_inv__ddy
      ddbcd1_z11ddy__ddbcd_length_inv11ddy=bcd_z
c dbcd1_z__ddz= bcd_z*dbcd_length_inv__ddz
      ddbcd1_z11ddz__dbcd_z=dbcd_length_inv__ddz
      ddbcd1_z11ddz__ddbcd_length_inv11ddz=bcd_z


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c block 3a

c df__dx=-y/(x**2 + y**2)
      ddf11dx__dx= -2*x/(x**2 + y**2) * df__dx
      ddf11dx__dy= (-2*y/(x**2 + y**2) + 1/y ) * df__dx
c df__dy=x/(x**2 + y**2)
      ddf11dy__dx= ddf11dx__dy
      ddf11dy__dy= -ddf11dx__dx


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c block 3b

c dabc1_x__dax= abc_x*dabc_length_inv__dax
      ddabc1_x11dax__dax=
     1  ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dax
      ddabc1_x11dax__day=
     1   ddabc1_x11dax__dabc_x*dabc_x__day
     1 + ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__day
      ddabc1_x11dax__daz=
     1   ddabc1_x11dax__dabc_x*dabc_x__daz
     1 + ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__daz
      ddabc1_x11dax__dbx=
     1  ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dbx
      ddabc1_x11dax__dby=
     1   ddabc1_x11dax__dabc_x*dabc_x__dby
     1 + ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dby
      ddabc1_x11dax__dbz=
     1   ddabc1_x11dax__dabc_x*dabc_x__dbz
     1 + ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dbz
      ddabc1_x11dax__dcx=
     1  ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcx
      ddabc1_x11dax__dcy=
     1   ddabc1_x11dax__dabc_x*dabc_x__dcy
     1 + ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcy
      ddabc1_x11dax__dcz=
     1   ddabc1_x11dax__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcz
c dabc1_x__day=- abc_length_inv*bc_z + abc_x*dabc_length_inv__day
      ddabc1_x11day__day=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__day
     1 + ddabc1_x11day__dabc_x*dabc_x__day
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__day
      ddabc1_x11day__daz=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__daz
     1 + ddabc1_x11day__dabc_x*dabc_x__daz
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__daz
      ddabc1_x11day__dbx=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__dbx
      ddabc1_x11day__dby=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_x11day__dabc_x*dabc_x__dby
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__dby
      ddabc1_x11day__dbz=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_x11day__dbc_z*dbc_z__dbz
     1 + ddabc1_x11day__dabc_x*dabc_x__dbz
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__dbz
      ddabc1_x11day__dcx=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__dcx
      ddabc1_x11day__dcy=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_x11day__dabc_x*dabc_x__dcy
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__dcy
      ddabc1_x11day__dcz=
     1   ddabc1_x11day__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_x11day__dbc_z*dbc_z__dcz
     1 + ddabc1_x11day__dabc_x*dabc_x__dcz
     1 + ddabc1_x11day__ddabc_length_inv11day*ddabc_length_inv11day__dcz
c dabc1_x__daz=abc_length_inv*bc_y + abc_x*dabc_length_inv__daz
      ddabc1_x11daz__daz=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__daz
     1 + ddabc1_x11daz__dabc_x*dabc_x__daz
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__daz
      ddabc1_x11daz__dbx=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dbx
      ddabc1_x11daz__dby=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_x11daz__dbc_y*dbc_y__dby
     1 + ddabc1_x11daz__dabc_x*dabc_x__dby
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dby
      ddabc1_x11daz__dbz=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_x11daz__dabc_x*dabc_x__dbz
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dbz
      ddabc1_x11daz__dcx=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcx
      ddabc1_x11daz__dcy=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_x11daz__dbc_y*dbc_y__dcy
     1 + ddabc1_x11daz__dabc_x*dabc_x__dcy
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcy
      ddabc1_x11daz__dcz=
     1   ddabc1_x11daz__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_x11daz__dabc_x*dabc_x__dcz
     1 + ddabc1_x11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcz
c dabc1_x__dbx= abc_x*dabc_length_inv__dbx
      ddabc1_x11dbx__dbx=
     1  ddabc1_x11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dbx
      ddabc1_x11dbx__dby=
     1   ddabc1_x11dbx__dabc_x*dabc_x__dby
     1 + ddabc1_x11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dby
      ddabc1_x11dbx__dbz=
     1   ddabc1_x11dbx__dabc_x*dabc_x__dbz
     1 + ddabc1_x11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dbz
      ddabc1_x11dbx__dcx=
     1  ddabc1_x11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcx
      ddabc1_x11dbx__dcy=
     1   ddabc1_x11dbx__dabc_x*dabc_x__dcy
     1 + ddabc1_x11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcy
      ddabc1_x11dbx__dcz=
     1   ddabc1_x11dbx__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcz
c dabc1_x__dby=abc_length_inv*dabc_x__dby + abc_x*dabc_length_inv__dby
      ddabc1_x11dby__dby=
     1   ddabc1_x11dby__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_x11dby__dabc_x*dabc_x__dby
     1 + ddabc1_x11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dby
      ddabc1_x11dby__dbz=
     1   ddabc1_x11dby__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_x11dby__ddabc_x11dby*ddabc_x11dby__dbz
     1 + ddabc1_x11dby__dabc_x*dabc_x__dbz
     1 + ddabc1_x11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dbz
      ddabc1_x11dby__dcx=
     1   ddabc1_x11dby__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_x11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcx
      ddabc1_x11dby__dcy=
     1   ddabc1_x11dby__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_x11dby__dabc_x*dabc_x__dcy
     1 + ddabc1_x11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcy
      ddabc1_x11dby__dcz=
     1   ddabc1_x11dby__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_x11dby__ddabc_x11dby*ddabc_x11dby__dcz
     1 + ddabc1_x11dby__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcz
c dabc1_x__dbz=abc_length_inv*dabc_x__dbz + abc_x*dabc_length_inv__dbz
      ddabc1_x11dbz__dbz=
     1   ddabc1_x11dbz__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_x11dbz__dabc_x*dabc_x__dbz
     1 + ddabc1_x11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dbz
      ddabc1_x11dbz__dcx=
     1   ddabc1_x11dbz__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_x11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcx
      ddabc1_x11dbz__dcy=
     1   ddabc1_x11dbz__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_x11dbz__ddabc_x11dbz*ddabc_x11dbz__dcy
     1 + ddabc1_x11dbz__dabc_x*dabc_x__dcy
     1 + ddabc1_x11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcy
      ddabc1_x11dbz__dcz=
     1   ddabc1_x11dbz__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_x11dbz__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcz
c dabc1_x__dcx= abc_x*dabc_length_inv__dcx
      ddabc1_x11dcx__dcx=
     1  ddabc1_x11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcx
      ddabc1_x11dcx__dcy=
     1   ddabc1_x11dcx__dabc_x*dabc_x__dcy
     1 + ddabc1_x11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcy
      ddabc1_x11dcx__dcz=
     1   ddabc1_x11dcx__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcz
c dabc1_x__dcy=-abc_length_inv*ab_z + abc_x*dabc_length_inv__dcy
      ddabc1_x11dcy__dcy=
     1   ddabc1_x11dcy__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_x11dcy__dabc_x*dabc_x__dcy
     1 + ddabc1_x11dcy__ddabc_length_inv11dcy*ddabc_length_inv11dcy__dcy
      ddabc1_x11dcy__dcz=
     1   ddabc1_x11dcy__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_x11dcy__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dcy__ddabc_length_inv11dcy*ddabc_length_inv11dcy__dcz
c dabc1_x__dcz=abc_length_inv*ab_y + abc_x*dabc_length_inv__dcz
      ddabc1_x11dcz__dcz=
     1   ddabc1_x11dcz__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_x11dcz__dabc_x*dabc_x__dcz
     1 + ddabc1_x11dcz__ddabc_length_inv11dcz*ddabc_length_inv11dcz__dcz

c dabc1_y__dax=abc_length_inv*bc_z + abc_y*dabc_length_inv__dax
      ddabc1_y11dax__dax=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dax
     1 + ddabc1_y11dax__dabc_y*dabc_y__dax
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dax
      ddabc1_y11dax__day=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__day
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__day
      ddabc1_y11dax__daz=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__daz
     1 + ddabc1_y11dax__dabc_y*dabc_y__daz
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__daz
      ddabc1_y11dax__dbx=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_y11dax__dabc_y*dabc_y__dbx
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dbx
      ddabc1_y11dax__dby=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dby
      ddabc1_y11dax__dbz=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_y11dax__dbc_z*dbc_z__dbz
     1 + ddabc1_y11dax__dabc_y*dabc_y__dbz
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dbz
      ddabc1_y11dax__dcx=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_y11dax__dabc_y*dabc_y__dcx
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcx
      ddabc1_y11dax__dcy=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcy
      ddabc1_y11dax__dcz=
     1   ddabc1_y11dax__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_y11dax__dbc_z*dbc_z__dcz
     1 + ddabc1_y11dax__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcz
c dabc1_y__day= abc_y*dabc_length_inv__day
      ddabc1_y11day__day=
     1  ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__day
      ddabc1_y11day__daz=
     1   ddabc1_y11day__dabc_y*dabc_y__daz
     1 + ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__daz
      ddabc1_y11day__dbx=
     1   ddabc1_y11day__dabc_y*dabc_y__dbx
     1 + ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__dbx
      ddabc1_y11day__dby=
     1  ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__dby
      ddabc1_y11day__dbz=
     1   ddabc1_y11day__dabc_y*dabc_y__dbz
     1 + ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__dbz
      ddabc1_y11day__dcx=
     1   ddabc1_y11day__dabc_y*dabc_y__dcx
     1 + ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__dcx
      ddabc1_y11day__dcy=
     1  ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__dcy
      ddabc1_y11day__dcz=
     1   ddabc1_y11day__dabc_y*dabc_y__dcz
     1 + ddabc1_y11day__ddabc_length_inv11day*ddabc_length_inv11day__dcz
c dabc1_y__daz=-abc_length_inv*bc_x + abc_y*dabc_length_inv__daz
      ddabc1_y11daz__daz=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__daz
     1 + ddabc1_y11daz__dabc_y*dabc_y__daz
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__daz
      ddabc1_y11daz__dbx=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_y11daz__dbc_x*dbc_x__dbx
     1 + ddabc1_y11daz__dabc_y*dabc_y__dbx
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dbx
      ddabc1_y11daz__dby=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dby
      ddabc1_y11daz__dbz=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_y11daz__dabc_y*dabc_y__dbz
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dbz
      ddabc1_y11daz__dcx=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_y11daz__dbc_x*dbc_x__dcx
     1 + ddabc1_y11daz__dabc_y*dabc_y__dcx
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcx
      ddabc1_y11daz__dcy=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcy
      ddabc1_y11daz__dcz=
     1   ddabc1_y11daz__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_y11daz__dabc_y*dabc_y__dcz
     1 + ddabc1_y11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcz
c dabc1_y__dbx=abc_length_inv*dabc_y__dbx + abc_y*dabc_length_inv__dbx
      ddabc1_y11dbx__dbx=
     1   ddabc1_y11dbx__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_y11dbx__dabc_y*dabc_y__dbx
     1 + ddabc1_y11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dbx
      ddabc1_y11dbx__dby=
     1   ddabc1_y11dbx__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_y11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dby
      ddabc1_y11dbx__dbz=
     1   ddabc1_y11dbx__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_y11dbx__ddabc_y11dbx*ddabc_y11dbx__dbz
     1 + ddabc1_y11dbx__dabc_y*dabc_y__dbz
     1 + ddabc1_y11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dbz
      ddabc1_y11dbx__dcx=
     1   ddabc1_y11dbx__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_y11dbx__dabc_y*dabc_y__dcx
     1 + ddabc1_y11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcx
      ddabc1_y11dbx__dcy=
     1   ddabc1_y11dbx__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_y11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcy
      ddabc1_y11dbx__dcz=
     1   ddabc1_y11dbx__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_y11dbx__ddabc_y11dbx*ddabc_y11dbx__dcz
     1 + ddabc1_y11dbx__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcz
c dabc1_y__dby= abc_y*dabc_length_inv__dby
      ddabc1_y11dby__dby=
     1  ddabc1_y11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dby
      ddabc1_y11dby__dbz=
     1   ddabc1_y11dby__dabc_y*dabc_y__dbz
     1 + ddabc1_y11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dbz
      ddabc1_y11dby__dcx=
     1   ddabc1_y11dby__dabc_y*dabc_y__dcx
     1 + ddabc1_y11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcx
      ddabc1_y11dby__dcy=
     1  ddabc1_y11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcy
      ddabc1_y11dby__dcz=
     1   ddabc1_y11dby__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcz
c dabc1_y__dbz=abc_length_inv*dabc_y__dbz + abc_y*dabc_length_inv__dbz
      ddabc1_y11dbz__dbz=
     1   ddabc1_y11dbz__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_y11dbz__dabc_y*dabc_y__dbz
     1 + ddabc1_y11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dbz
      ddabc1_y11dbz__dcx=
     1   ddabc1_y11dbz__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_y11dbz__ddabc_y11dbz*ddabc_y11dbz__dcx
     1 + ddabc1_y11dbz__dabc_y*dabc_y__dcx
     1 + ddabc1_y11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcx
      ddabc1_y11dbz__dcy=
     1   ddabc1_y11dbz__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_y11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcy
      ddabc1_y11dbz__dcz=
     1   ddabc1_y11dbz__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_y11dbz__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcz
c dabc1_y__dcx=abc_length_inv*ab_z + abc_y*dabc_length_inv__dcx
      ddabc1_y11dcx__dcx=
     1   ddabc1_y11dcx__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_y11dcx__dabc_y*dabc_y__dcx
     1 + ddabc1_y11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcx
      ddabc1_y11dcx__dcy=
     1   ddabc1_y11dcx__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_y11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcy
      ddabc1_y11dcx__dcz=
     1   ddabc1_y11dcx__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_y11dcx__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcz
c dabc1_y__dcy= abc_y*dabc_length_inv__dcy
      ddabc1_y11dcy__dcy=
     1  ddabc1_y11dcy__ddabc_length_inv11dcy*ddabc_length_inv11dcy__dcy
      ddabc1_y11dcy__dcz=
     1   ddabc1_y11dcy__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dcy__ddabc_length_inv11dcy*ddabc_length_inv11dcy__dcz
c dabc1_y__dcz=-abc_length_inv*ab_x + abc_y*dabc_length_inv__dcz
      ddabc1_y11dcz__dcz=
     1   ddabc1_y11dcz__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_y11dcz__dabc_y*dabc_y__dcz
     1 + ddabc1_y11dcz__ddabc_length_inv11dcz*ddabc_length_inv11dcz__dcz

c dabc1_z__dax=-abc_length_inv*bc_y + abc_z*dabc_length_inv__dax
      ddabc1_z11dax__dax=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dax
     1 + ddabc1_z11dax__dabc_z*dabc_z__dax
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dax
      ddabc1_z11dax__day=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__day
     1 + ddabc1_z11dax__dabc_z*dabc_z__day
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__day
      ddabc1_z11dax__daz=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__daz
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__daz
      ddabc1_z11dax__dbx=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_z11dax__dabc_z*dabc_z__dbx
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dbx
      ddabc1_z11dax__dby=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_z11dax__dbc_y*dbc_y__dby
     1 + ddabc1_z11dax__dabc_z*dabc_z__dby
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dby
      ddabc1_z11dax__dbz=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dbz
      ddabc1_z11dax__dcx=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_z11dax__dabc_z*dabc_z__dcx
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcx
      ddabc1_z11dax__dcy=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_z11dax__dbc_y*dbc_y__dcy
     1 + ddabc1_z11dax__dabc_z*dabc_z__dcy
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcy
      ddabc1_z11dax__dcz=
     1   ddabc1_z11dax__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_z11dax__ddabc_length_inv11dax*ddabc_length_inv11dax__dcz
c dabc1_z__day=abc_length_inv*bc_x + abc_z*dabc_length_inv__day
      ddabc1_z11day__day=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__day
     1 + ddabc1_z11day__dabc_z*dabc_z__day
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__day
      ddabc1_z11day__daz=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__daz
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__daz
      ddabc1_z11day__dbx=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_z11day__dbc_x*dbc_x__dbx
     1 + ddabc1_z11day__dabc_z*dabc_z__dbx
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__dbx
      ddabc1_z11day__dby=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_z11day__dabc_z*dabc_z__dby
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__dby
      ddabc1_z11day__dbz=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__dbz
      ddabc1_z11day__dcx=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_z11day__dbc_x*dbc_x__dcx
     1 + ddabc1_z11day__dabc_z*dabc_z__dcx
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__dcx
      ddabc1_z11day__dcy=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_z11day__dabc_z*dabc_z__dcy
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__dcy
      ddabc1_z11day__dcz=
     1   ddabc1_z11day__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_z11day__ddabc_length_inv11day*ddabc_length_inv11day__dcz
c dabc1_z__daz= abc_z*dabc_length_inv__daz
      ddabc1_z11daz__daz=
     1  ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__daz
      ddabc1_z11daz__dbx=
     1   ddabc1_z11daz__dabc_z*dabc_z__dbx
     1 + ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dbx
      ddabc1_z11daz__dby=
     1   ddabc1_z11daz__dabc_z*dabc_z__dby
     1 + ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dby
      ddabc1_z11daz__dbz=
     1  ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dbz
      ddabc1_z11daz__dcx=
     1   ddabc1_z11daz__dabc_z*dabc_z__dcx
     1 + ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcx
      ddabc1_z11daz__dcy=
     1   ddabc1_z11daz__dabc_z*dabc_z__dcy
     1 + ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcy
      ddabc1_z11daz__dcz=
     1  ddabc1_z11daz__ddabc_length_inv11daz*ddabc_length_inv11daz__dcz
c dabc1_z__dbx=abc_length_inv*dabc_z__dbx + abc_z*dabc_length_inv__dbx
      ddabc1_z11dbx__dbx=
     1   ddabc1_z11dbx__dabc_length_inv*dabc_length_inv__dbx
     1 + ddabc1_z11dbx__dabc_z*dabc_z__dbx
     1 + ddabc1_z11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dbx
      ddabc1_z11dbx__dby=
     1   ddabc1_z11dbx__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_z11dbx__ddabc_z11dbx*ddabc_z11dbx__dby
     1 + ddabc1_z11dbx__dabc_z*dabc_z__dby
     1 + ddabc1_z11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dby
      ddabc1_z11dbx__dbz=
     1   ddabc1_z11dbx__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_z11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dbz
      ddabc1_z11dbx__dcx=
     1   ddabc1_z11dbx__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_z11dbx__dabc_z*dabc_z__dcx
     1 + ddabc1_z11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcx
      ddabc1_z11dbx__dcy=
     1   ddabc1_z11dbx__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_z11dbx__ddabc_z11dbx*ddabc_z11dbx__dcy
     1 + ddabc1_z11dbx__dabc_z*dabc_z__dcy
     1 + ddabc1_z11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcy
      ddabc1_z11dbx__dcz=
     1   ddabc1_z11dbx__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_z11dbx__ddabc_length_inv11dbx*ddabc_length_inv11dbx__dcz
c dabc1_z__dby=abc_length_inv*dabc_z__dby + abc_z*dabc_length_inv__dby
      ddabc1_z11dby__dby=
     1   ddabc1_z11dby__dabc_length_inv*dabc_length_inv__dby
     1 + ddabc1_z11dby__dabc_z*dabc_z__dby
     1 + ddabc1_z11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dby
      ddabc1_z11dby__dbz=
     1   ddabc1_z11dby__dabc_length_inv*dabc_length_inv__dbz
     1 + ddabc1_z11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dbz
      ddabc1_z11dby__dcx=
     1   ddabc1_z11dby__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_z11dby__ddabc_z11dby*ddabc_z11dby__dcx
     1 + ddabc1_z11dby__dabc_z*dabc_z__dcx
     1 + ddabc1_z11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcx
      ddabc1_z11dby__dcy=
     1   ddabc1_z11dby__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_z11dby__dabc_z*dabc_z__dcy
     1 + ddabc1_z11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcy
      ddabc1_z11dby__dcz=
     1   ddabc1_z11dby__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_z11dby__ddabc_length_inv11dby*ddabc_length_inv11dby__dcz
c dabc1_z__dbz= abc_z*dabc_length_inv__dbz
      ddabc1_z11dbz__dbz=
     1  ddabc1_z11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dbz
      ddabc1_z11dbz__dcx=
     1   ddabc1_z11dbz__dabc_z*dabc_z__dcx
     1 + ddabc1_z11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcx
      ddabc1_z11dbz__dcy=
     1   ddabc1_z11dbz__dabc_z*dabc_z__dcy
     1 + ddabc1_z11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcy
      ddabc1_z11dbz__dcz=
     1  ddabc1_z11dbz__ddabc_length_inv11dbz*ddabc_length_inv11dbz__dcz
c dabc1_z__dcx=-abc_length_inv*ab_y + abc_z*dabc_length_inv__dcx
      ddabc1_z11dcx__dcx=
     1   ddabc1_z11dcx__dabc_length_inv*dabc_length_inv__dcx
     1 + ddabc1_z11dcx__dabc_z*dabc_z__dcx
     1 + ddabc1_z11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcx
      ddabc1_z11dcx__dcy=
     1   ddabc1_z11dcx__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_z11dcx__dabc_z*dabc_z__dcy
     1 + ddabc1_z11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcy
      ddabc1_z11dcx__dcz=
     1   ddabc1_z11dcx__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_z11dcx__ddabc_length_inv11dcx*ddabc_length_inv11dcx__dcz
c dabc1_z__dcy=abc_length_inv*ab_x + abc_z*dabc_length_inv__dcy
      ddabc1_z11dcy__dcy=
     1   ddabc1_z11dcy__dabc_length_inv*dabc_length_inv__dcy
     1 + ddabc1_z11dcy__dabc_z*dabc_z__dcy
     1 + ddabc1_z11dcy__ddabc_length_inv11dcy*ddabc_length_inv11dcy__dcy
      ddabc1_z11dcy__dcz=
     1   ddabc1_z11dcy__dabc_length_inv*dabc_length_inv__dcz
     1 + ddabc1_z11dcy__ddabc_length_inv11dcy*ddabc_length_inv11dcy__dcz
c dabc1_z__dcz= abc_z*dabc_length_inv__dcz
      ddabc1_z11dcz__dcz=
     1  ddabc1_z11dcz__ddabc_length_inv11dcz*ddabc_length_inv11dcz__dcz


c dx__dax= bcd1_x*dabc1_x__dax + bcd1_y*dabc1_y__dax + bcd1_z*dabc1_z__dax  
      ddx11dax__dbcd1_x=dabc1_x__dax
      ddx11dax__dbcd1_y=dabc1_y__dax
      ddx11dax__dbcd1_z=dabc1_z__dax
      ddx11dax__ddabc1_x11dax=bcd1_x
      ddx11dax__ddabc1_y11dax=bcd1_y
      ddx11dax__ddabc1_z11dax=bcd1_z

c dx__day= bcd1_x*dabc1_x__day + bcd1_y*dabc1_y__day + bcd1_z*dabc1_z__day 
      ddx11day__dbcd1_x=dabc1_x__day
      ddx11day__dbcd1_y=dabc1_y__day
      ddx11day__dbcd1_z=dabc1_z__day
      ddx11day__ddabc1_x11day=bcd1_x
      ddx11day__ddabc1_y11day=bcd1_y
      ddx11day__ddabc1_z11day=bcd1_z

c dx__daz= bcd1_x*dabc1_x__daz + bcd1_y*dabc1_y__daz + bcd1_z*dabc1_z__daz 
      ddx11daz__dbcd1_x=dabc1_x__daz
      ddx11daz__dbcd1_y=dabc1_y__daz
      ddx11daz__dbcd1_z=dabc1_z__daz
      ddx11daz__ddabc1_x11daz=bcd1_x
      ddx11daz__ddabc1_y11daz=bcd1_y
      ddx11daz__ddabc1_z11daz=bcd1_z

c dx__dbx= bcd1_x*dabc1_x__dbx + abc1_x*dbcd1_x__dbx + bcd1_y*dabc1_y__dbx + abc1_y*dbcd1_y__dbx + bcd1_z*dabc1_z__dbx + abc1_z*dbcd1_z__dbx
      ddx11dbx__dabc1_x=dbcd1_x__dbx
      ddx11dbx__dabc1_y=dbcd1_y__dbx
      ddx11dbx__dabc1_z=dbcd1_z__dbx
      ddx11dbx__dbcd1_x=dabc1_x__dbx
      ddx11dbx__dbcd1_y=dabc1_y__dbx
      ddx11dbx__dbcd1_z=dabc1_z__dbx
      ddx11dbx__ddabc1_x11dbx=bcd1_x
      ddx11dbx__ddabc1_y11dbx=bcd1_y
      ddx11dbx__ddabc1_z11dbx=bcd1_z
      ddx11dbx__ddbcd1_x11dbx=abc1_x
      ddx11dbx__ddbcd1_y11dbx=abc1_y
      ddx11dbx__ddbcd1_z11dbx=abc1_z

c dx__dby= bcd1_x*dabc1_x__dby + abc1_x*dbcd1_x__dby + bcd1_y*dabc1_y__dby + abc1_y*dbcd1_y__dby + bcd1_z*dabc1_z__dby + abc1_z*dbcd1_z__dby
      ddx11dby__dabc1_x=dbcd1_x__dby
      ddx11dby__dabc1_y=dbcd1_y__dby
      ddx11dby__dabc1_z=dbcd1_z__dby
      ddx11dby__dbcd1_x=dabc1_x__dby
      ddx11dby__dbcd1_y=dabc1_y__dby
      ddx11dby__dbcd1_z=dabc1_z__dby
      ddx11dby__ddabc1_x11dby=bcd1_x
      ddx11dby__ddabc1_y11dby=bcd1_y
      ddx11dby__ddabc1_z11dby=bcd1_z
      ddx11dby__ddbcd1_x11dby=abc1_x
      ddx11dby__ddbcd1_y11dby=abc1_y
      ddx11dby__ddbcd1_z11dby=abc1_z

c dx__dbz= bcd1_x*dabc1_x__dbz + abc1_x*dbcd1_x__dbz + bcd1_y*dabc1_y__dbz + abc1_y*dbcd1_y__dbz + bcd1_z*dabc1_z__dbz + abc1_z*dbcd1_z__dbz
      ddx11dbz__dabc1_x=dbcd1_x__dbz
      ddx11dbz__dabc1_y=dbcd1_y__dbz
      ddx11dbz__dabc1_z=dbcd1_z__dbz
      ddx11dbz__dbcd1_x=dabc1_x__dbz
      ddx11dbz__dbcd1_y=dabc1_y__dbz
      ddx11dbz__dbcd1_z=dabc1_z__dbz
      ddx11dbz__ddabc1_x11dbz=bcd1_x
      ddx11dbz__ddabc1_y11dbz=bcd1_y
      ddx11dbz__ddabc1_z11dbz=bcd1_z
      ddx11dbz__ddbcd1_x11dbz=abc1_x
      ddx11dbz__ddbcd1_y11dbz=abc1_y
      ddx11dbz__ddbcd1_z11dbz=abc1_z

c dx__dcx= bcd1_x*dabc1_x__dcx + abc1_x*dbcd1_x__dcx + bcd1_y*dabc1_y__dcx + abc1_y*dbcd1_y__dcx + bcd1_z*dabc1_z__dcx + abc1_z*dbcd1_z__dcx
      ddx11dcx__dabc1_x=dbcd1_x__dcx
      ddx11dcx__dabc1_y=dbcd1_y__dcx
      ddx11dcx__dabc1_z=dbcd1_z__dcx
      ddx11dcx__dbcd1_x=dabc1_x__dcx
      ddx11dcx__dbcd1_y=dabc1_y__dcx
      ddx11dcx__dbcd1_z=dabc1_z__dcx
      ddx11dcx__ddabc1_x11dcx=bcd1_x
      ddx11dcx__ddabc1_y11dcx=bcd1_y
      ddx11dcx__ddabc1_z11dcx=bcd1_z
      ddx11dcx__ddbcd1_x11dcx=abc1_x
      ddx11dcx__ddbcd1_y11dcx=abc1_y
      ddx11dcx__ddbcd1_z11dcx=abc1_z

c dx__dcy= bcd1_x*dabc1_x__dcy + abc1_x*dbcd1_x__dcy + bcd1_y*dabc1_y__dcy + abc1_y*dbcd1_y__dcy + bcd1_z*dabc1_z__dcy + abc1_z*dbcd1_z__dcy
      ddx11dcy__dabc1_x=dbcd1_x__dcy
      ddx11dcy__dabc1_y=dbcd1_y__dcy
      ddx11dcy__dabc1_z=dbcd1_z__dcy
      ddx11dcy__dbcd1_x=dabc1_x__dcy
      ddx11dcy__dbcd1_y=dabc1_y__dcy
      ddx11dcy__dbcd1_z=dabc1_z__dcy
      ddx11dcy__ddabc1_x11dcy=bcd1_x
      ddx11dcy__ddabc1_y11dcy=bcd1_y
      ddx11dcy__ddabc1_z11dcy=bcd1_z
      ddx11dcy__ddbcd1_x11dcy=abc1_x
      ddx11dcy__ddbcd1_y11dcy=abc1_y
      ddx11dcy__ddbcd1_z11dcy=abc1_z

c dx__dcz= bcd1_x*dabc1_x__dcz + abc1_x*dbcd1_x__dcz + bcd1_y*dabc1_y__dcz + abc1_y*dbcd1_y__dcz + bcd1_z*dabc1_z__dcz + abc1_z*dbcd1_z__dcz
      ddx11dcz__dabc1_x=dbcd1_x__dcz
      ddx11dcz__dabc1_y=dbcd1_y__dcz
      ddx11dcz__dabc1_z=dbcd1_z__dcz
      ddx11dcz__dbcd1_x=dabc1_x__dcz
      ddx11dcz__dbcd1_y=dabc1_y__dcz
      ddx11dcz__dbcd1_z=dabc1_z__dcz
      ddx11dcz__ddabc1_x11dcz=bcd1_x
      ddx11dcz__ddabc1_y11dcz=bcd1_y
      ddx11dcz__ddabc1_z11dcz=bcd1_z
      ddx11dcz__ddbcd1_x11dcz=abc1_x
      ddx11dcz__ddbcd1_y11dcz=abc1_y
      ddx11dcz__ddbcd1_z11dcz=abc1_z

c dx__ddx= abc1_x*dbcd1_x__ddx + abc1_y*dbcd1_y__ddx + abc1_z*dbcd1_z__ddx
      ddx11ddx__dabc1_x=dbcd1_x__ddx
      ddx11ddx__dabc1_y=dbcd1_y__ddx
      ddx11ddx__dabc1_z=dbcd1_z__ddx
      ddx11ddx__ddbcd1_x11ddx=abc1_x
      ddx11ddx__ddbcd1_y11ddx=abc1_y
      ddx11ddx__ddbcd1_z11ddx=abc1_z

c dx__ddy= abc1_x*dbcd1_x__ddy + abc1_y*dbcd1_y__ddy + abc1_z*dbcd1_z__ddy
      ddx11ddy__dabc1_x=dbcd1_x__ddy
      ddx11ddy__dabc1_y=dbcd1_y__ddy
      ddx11ddy__dabc1_z=dbcd1_z__ddy
      ddx11ddy__ddbcd1_x11ddy=abc1_x
      ddx11ddy__ddbcd1_y11ddy=abc1_y
      ddx11ddy__ddbcd1_z11ddy=abc1_z

c dx__ddz= abc1_x*dbcd1_x__ddz + abc1_y*dbcd1_y__ddz + abc1_z*dbcd1_z__ddz
      ddx11ddz__dabc1_x=dbcd1_x__ddz
      ddx11ddz__dabc1_y=dbcd1_y__ddz
      ddx11ddz__dabc1_z=dbcd1_z__ddz
      ddx11ddz__ddbcd1_x11ddz=abc1_x
      ddx11ddz__ddbcd1_y11ddz=abc1_y
      ddx11ddz__ddbcd1_z11ddz=abc1_z


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c block 3c

c daux_x__dax= bc1_z*dabc1_y__dax - bc1_y*dabc1_z__dax
      ddaux_x11dax__dax=
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dax
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dax
      ddaux_x11dax__day=
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__day
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__day
      ddaux_x11dax__daz=
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__daz
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__daz
      ddaux_x11dax__dbx=
     1   ddaux_x11dax__dbc1_z*dbc1_z__dbx
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dbx
     1 + ddaux_x11dax__dbc1_y*dbc1_y__dbx
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dbx
      ddaux_x11dax__dby=
     1   ddaux_x11dax__dbc1_z*dbc1_z__dby
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dby
     1 + ddaux_x11dax__dbc1_y*dbc1_y__dby
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dby
      ddaux_x11dax__dbz=
     1   ddaux_x11dax__dbc1_z*dbc1_z__dbz
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dbz
     1 + ddaux_x11dax__dbc1_y*dbc1_y__dbz
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dbz
      ddaux_x11dax__dcx=
     1   ddaux_x11dax__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dcx
     1 + ddaux_x11dax__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dcx
      ddaux_x11dax__dcy=
     1   ddaux_x11dax__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dcy
     1 + ddaux_x11dax__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dcy
      ddaux_x11dax__dcz=
     1   ddaux_x11dax__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dax__ddabc1_y11dax*ddabc1_y11dax__dcz
     1 + ddaux_x11dax__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dax__ddabc1_z11dax*ddabc1_z11dax__dcz
c daux_x__day= bc1_z*dabc1_y__day - bc1_y*dabc1_z__day
      ddaux_x11day__day=
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__day
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__day
      ddaux_x11day__daz=
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__daz
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__daz
      ddaux_x11day__dbx=
     1   ddaux_x11day__dbc1_z*dbc1_z__dbx
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__dbx
     1 + ddaux_x11day__dbc1_y*dbc1_y__dbx
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__dbx
      ddaux_x11day__dby=
     1   ddaux_x11day__dbc1_z*dbc1_z__dby
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__dby
     1 + ddaux_x11day__dbc1_y*dbc1_y__dby
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__dby
      ddaux_x11day__dbz=
     1   ddaux_x11day__dbc1_z*dbc1_z__dbz
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__dbz
     1 + ddaux_x11day__dbc1_y*dbc1_y__dbz
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__dbz
      ddaux_x11day__dcx=
     1   ddaux_x11day__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__dcx
     1 + ddaux_x11day__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__dcx
      ddaux_x11day__dcy=
     1   ddaux_x11day__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__dcy
     1 + ddaux_x11day__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__dcy
      ddaux_x11day__dcz=
     1   ddaux_x11day__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11day__ddabc1_y11dax*ddabc1_y11day__dcz
     1 + ddaux_x11day__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11day__ddabc1_z11dax*ddabc1_z11day__dcz
c daux_x__daz= bc1_z*dabc1_y__daz - bc1_y*dabc1_z__daz
      ddaux_x11daz__daz=
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__daz
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__daz
      ddaux_x11daz__dbx=
     1   ddaux_x11daz__dbc1_z*dbc1_z__dbx
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__dbx
     1 + ddaux_x11daz__dbc1_y*dbc1_y__dbx
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__dbx
      ddaux_x11daz__dby=
     1   ddaux_x11daz__dbc1_z*dbc1_z__dby
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__dby
     1 + ddaux_x11daz__dbc1_y*dbc1_y__dby
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__dby
      ddaux_x11daz__dbz=
     1   ddaux_x11daz__dbc1_z*dbc1_z__dbz
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__dbz
     1 + ddaux_x11daz__dbc1_y*dbc1_y__dbz
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__dbz
      ddaux_x11daz__dcx=
     1   ddaux_x11daz__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__dcx
     1 + ddaux_x11daz__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__dcx
      ddaux_x11daz__dcy=
     1   ddaux_x11daz__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__dcy
     1 + ddaux_x11daz__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__dcy
      ddaux_x11daz__dcz=
     1   ddaux_x11daz__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11daz__ddabc1_y11dax*ddabc1_y11daz__dcz
     1 + ddaux_x11daz__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11daz__ddabc1_z11dax*ddabc1_z11daz__dcz
c daux_x__dbx= bc1_z*dabc1_y__dbx + abc1_y*dbc1_z__dbx - abc1_z*dbc1_y__dbx - bc1_y*dabc1_z__dbx
      ddaux_x11dbx__dbx=
     1   ddaux_x11dbx__dbc1_z*dbc1_z__dbx
     1 + ddaux_x11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dbx
     1 + ddaux_x11dbx__dabc1_y*dabc1_y__dbx
     1 + ddaux_x11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dbx
     1 + ddaux_x11dbx__dabc1_z*dabc1_z__dbx
     1 + ddaux_x11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dbx
     1 + ddaux_x11dbx__dbc1_y*dbc1_y__dbx
     1 + ddaux_x11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dbx
      ddaux_x11dbx__dby=
     1   ddaux_x11dbx__dbc1_z*dbc1_z__dby
     1 + ddaux_x11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dby
     1 + ddaux_x11dbx__dabc1_y*dabc1_y__dby
     1 + ddaux_x11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dby
     1 + ddaux_x11dbx__dabc1_z*dabc1_z__dby
     1 + ddaux_x11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dby
     1 + ddaux_x11dbx__dbc1_y*dbc1_y__dby
     1 + ddaux_x11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dby
      ddaux_x11dbx__dbz=
     1   ddaux_x11dbx__dbc1_z*dbc1_z__dbz
     1 + ddaux_x11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dbz
     1 + ddaux_x11dbx__dabc1_y*dabc1_y__dbz
     1 + ddaux_x11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dbz
     1 + ddaux_x11dbx__dabc1_z*dabc1_z__dbz
     1 + ddaux_x11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dbz
     1 + ddaux_x11dbx__dbc1_y*dbc1_y__dbz
     1 + ddaux_x11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dbz
      ddaux_x11dbx__dcx=
     1   ddaux_x11dbx__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcx
     1 + ddaux_x11dbx__dabc1_y*dabc1_y__dcx
     1 + ddaux_x11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dcx
     1 + ddaux_x11dbx__dabc1_z*dabc1_z__dcx
     1 + ddaux_x11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dcx
     1 + ddaux_x11dbx__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcx
      ddaux_x11dbx__dcy=
     1   ddaux_x11dbx__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcy
     1 + ddaux_x11dbx__dabc1_y*dabc1_y__dcy
     1 + ddaux_x11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dcy
     1 + ddaux_x11dbx__dabc1_z*dabc1_z__dcy
     1 + ddaux_x11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dcy
     1 + ddaux_x11dbx__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcy
      ddaux_x11dbx__dcz=
     1   ddaux_x11dbx__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcz
     1 + ddaux_x11dbx__dabc1_y*dabc1_y__dcz
     1 + ddaux_x11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dcz
     1 + ddaux_x11dbx__dabc1_z*dabc1_z__dcz
     1 + ddaux_x11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dcz
     1 + ddaux_x11dbx__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcz
c daux_x__dby= bc1_z*dabc1_y__dby + abc1_y*dbc1_z__dby - abc1_z*dbc1_y__dby - bc1_y*dabc1_z__dby
      ddaux_x11dby__dby=
     1   ddaux_x11dby__dbc1_z*dbc1_z__dby
     1 + ddaux_x11dby__ddabc1_y11dby*ddabc1_y11dby__dby
     1 + ddaux_x11dby__dabc1_y*dabc1_y__dby
     1 + ddaux_x11dby__ddbc1_z11dby*ddbc1_z11dby__dby
     1 + ddaux_x11dby__dabc1_z*dabc1_z__dby
     1 + ddaux_x11dby__ddbc1_y11dby*ddbc1_y11dby__dby
     1 + ddaux_x11dby__dbc1_y*dbc1_y__dby
     1 + ddaux_x11dby__ddabc1_z11dby*ddabc1_z11dby__dby
      ddaux_x11dby__dbz=
     1   ddaux_x11dby__dbc1_z*dbc1_z__dbz
     1 + ddaux_x11dby__ddabc1_y11dby*ddabc1_y11dby__dbz
     1 + ddaux_x11dby__dabc1_y*dabc1_y__dbz
     1 + ddaux_x11dby__ddbc1_z11dby*ddbc1_z11dby__dbz
     1 + ddaux_x11dby__dabc1_z*dabc1_z__dbz
     1 + ddaux_x11dby__ddbc1_y11dby*ddbc1_y11dby__dbz
     1 + ddaux_x11dby__dbc1_y*dbc1_y__dbz
     1 + ddaux_x11dby__ddabc1_z11dby*ddabc1_z11dby__dbz
      ddaux_x11dby__dcx=
     1   ddaux_x11dby__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11dby__ddabc1_y11dby*ddabc1_y11dby__dcx
     1 + ddaux_x11dby__dabc1_y*dabc1_y__dcx
     1 + ddaux_x11dby__ddbc1_z11dby*ddbc1_z11dby__dcx
     1 + ddaux_x11dby__dabc1_z*dabc1_z__dcx
     1 + ddaux_x11dby__ddbc1_y11dby*ddbc1_y11dby__dcx
     1 + ddaux_x11dby__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11dby__ddabc1_z11dby*ddabc1_z11dby__dcx
      ddaux_x11dby__dcy=
     1   ddaux_x11dby__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11dby__ddabc1_y11dby*ddabc1_y11dby__dcy
     1 + ddaux_x11dby__dabc1_y*dabc1_y__dcy
     1 + ddaux_x11dby__ddbc1_z11dby*ddbc1_z11dby__dcy
     1 + ddaux_x11dby__dabc1_z*dabc1_z__dcy
     1 + ddaux_x11dby__ddbc1_y11dby*ddbc1_y11dby__dcy
     1 + ddaux_x11dby__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11dby__ddabc1_z11dby*ddabc1_z11dby__dcy
      ddaux_x11dby__dcz=
     1   ddaux_x11dby__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dby__ddabc1_y11dby*ddabc1_y11dby__dcz
     1 + ddaux_x11dby__dabc1_y*dabc1_y__dcz
     1 + ddaux_x11dby__ddbc1_z11dby*ddbc1_z11dby__dcz
     1 + ddaux_x11dby__dabc1_z*dabc1_z__dcz
     1 + ddaux_x11dby__ddbc1_y11dby*ddbc1_y11dby__dcz
     1 + ddaux_x11dby__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dby__ddabc1_z11dby*ddabc1_z11dby__dcz
c daux_x__dbz= bc1_z*dabc1_y__dbz + abc1_y*dbc1_z__dbz - abc1_z*dbc1_y__dbz - bc1_y*dabc1_z__dbz
      ddaux_x11dbz__dbz=
     1   ddaux_x11dbz__dbc1_z*dbc1_z__dbz
     1 + ddaux_x11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dbz
     1 + ddaux_x11dbz__dabc1_y*dabc1_y__dbz
     1 + ddaux_x11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dbz
     1 + ddaux_x11dbz__dabc1_z*dabc1_z__dbz
     1 + ddaux_x11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dbz
     1 + ddaux_x11dbz__dbc1_y*dbc1_y__dbz
     1 + ddaux_x11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dbz
      ddaux_x11dbz__dcx=
     1   ddaux_x11dbz__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcx
     1 + ddaux_x11dbz__dabc1_y*dabc1_y__dcx
     1 + ddaux_x11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dcx
     1 + ddaux_x11dbz__dabc1_z*dabc1_z__dcx
     1 + ddaux_x11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dcx
     1 + ddaux_x11dbz__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcx
      ddaux_x11dbz__dcy=
     1   ddaux_x11dbz__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcy
     1 + ddaux_x11dbz__dabc1_y*dabc1_y__dcy
     1 + ddaux_x11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dcy
     1 + ddaux_x11dbz__dabc1_z*dabc1_z__dcy
     1 + ddaux_x11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dcy
     1 + ddaux_x11dbz__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcy
      ddaux_x11dbz__dcz=
     1   ddaux_x11dbz__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcz
     1 + ddaux_x11dbz__dabc1_y*dabc1_y__dcz
     1 + ddaux_x11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dcz
     1 + ddaux_x11dbz__dabc1_z*dabc1_z__dcz
     1 + ddaux_x11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dcz
     1 + ddaux_x11dbz__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcz
c daux_x__dcx= bc1_z*dabc1_y__dcx + abc1_y*dbc1_z__dcx - abc1_z*dbc1_y__dcx - bc1_y*dabc1_z__dcx
      ddaux_x11dcx__dcx=
     1   ddaux_x11dcx__dbc1_z*dbc1_z__dcx
     1 + ddaux_x11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcx
     1 + ddaux_x11dcx__dabc1_y*dabc1_y__dcx
     1 + ddaux_x11dcx__ddbc1_z11dcx*ddbc1_z11dcx__dcx
     1 + ddaux_x11dcx__dabc1_z*dabc1_z__dcx
     1 + ddaux_x11dcx__ddbc1_y11dcx*ddbc1_y11dcx__dcx
     1 + ddaux_x11dcx__dbc1_y*dbc1_y__dcx
     1 + ddaux_x11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcx
      ddaux_x11dcx__dcy=
     1   ddaux_x11dcx__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcy
     1 + ddaux_x11dcx__dabc1_y*dabc1_y__dcy
     1 + ddaux_x11dcx__ddbc1_z11dcx*ddbc1_z11dcx__dcy
     1 + ddaux_x11dcx__dabc1_z*dabc1_z__dcy
     1 + ddaux_x11dcx__ddbc1_y11dcx*ddbc1_y11dcx__dcy
     1 + ddaux_x11dcx__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcy
      ddaux_x11dcx__dcz=
     1   ddaux_x11dcx__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcz
     1 + ddaux_x11dcx__dabc1_y*dabc1_y__dcz
     1 + ddaux_x11dcx__ddbc1_z11dcx*ddbc1_z11dcx__dcz
     1 + ddaux_x11dcx__dabc1_z*dabc1_z__dcz
     1 + ddaux_x11dcx__ddbc1_y11dcx*ddbc1_y11dcx__dcz
     1 + ddaux_x11dcx__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcz
c daux_x__dcy= bc1_z*dabc1_y__dcy + abc1_y*dbc1_z__dcy - abc1_z*dbc1_y__dcy - bc1_y*dabc1_z__dcy
      ddaux_x11dcy__dcy=
     1   ddaux_x11dcy__dbc1_z*dbc1_z__dcy
     1 + ddaux_x11dcy__ddabc1_y11dcy*ddabc1_y11dcy__dcy
     1 + ddaux_x11dcy__dabc1_y*dabc1_y__dcy
     1 + ddaux_x11dcy__ddbc1_z11dcy*ddbc1_z11dcy__dcy
     1 + ddaux_x11dcy__dabc1_z*dabc1_z__dcy
     1 + ddaux_x11dcy__ddbc1_y11dcy*ddbc1_y11dcy__dcy
     1 + ddaux_x11dcy__dbc1_y*dbc1_y__dcy
     1 + ddaux_x11dcy__ddabc1_z11dcy*ddabc1_z11dcy__dcy
      ddaux_x11dcy__dcz=
     1   ddaux_x11dcy__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dcy__ddabc1_y11dcy*ddabc1_y11dcy__dcz
     1 + ddaux_x11dcy__dabc1_y*dabc1_y__dcz
     1 + ddaux_x11dcy__ddbc1_z11dcy*ddbc1_z11dcy__dcz
     1 + ddaux_x11dcy__dabc1_z*dabc1_z__dcz
     1 + ddaux_x11dcy__ddbc1_y11dcy*ddbc1_y11dcy__dcz
     1 + ddaux_x11dcy__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dcy__ddabc1_z11dcy*ddabc1_z11dcy__dcz
c daux_x__dcz= bc1_z*dabc1_y__dcz + abc1_y*dbc1_z__dcz - abc1_z*dbc1_y__dcz - bc1_y*dabc1_z__dcz
      ddaux_x11dcz__dcz=
     1   ddaux_x11dcz__dbc1_z*dbc1_z__dcz
     1 + ddaux_x11dcz__ddabc1_y11dcz*ddabc1_y11dcz__dcz
     1 + ddaux_x11dcz__dabc1_y*dabc1_y__dcz
     1 + ddaux_x11dcz__ddbc1_z11dcz*ddbc1_z11dcz__dcz
     1 + ddaux_x11dcz__dabc1_z*dabc1_z__dcz
     1 + ddaux_x11dcz__ddbc1_y11dcz*ddbc1_y11dcz__dcz
     1 + ddaux_x11dcz__dbc1_y*dbc1_y__dcz
     1 + ddaux_x11dcz__ddabc1_z11dcz*ddabc1_z11dcz__dcz


c daux_y__dax= bc1_x*dabc1_z__dax - bc1_z*dabc1_x__dax
      ddaux_y11dax__dax=
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dax
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dax
      ddaux_y11dax__day=
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__day
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__day
      ddaux_y11dax__daz=
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__daz
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__daz
      ddaux_y11dax__dbx=
     1   ddaux_y11dax__dbc1_x*dbc1_x__dbx
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dbx
     1 + ddaux_y11dax__dbc1_z*dbc1_z__dbx
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dbx
      ddaux_y11dax__dby=
     1   ddaux_y11dax__dbc1_x*dbc1_x__dby
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dby
     1 + ddaux_y11dax__dbc1_z*dbc1_z__dby
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dby
      ddaux_y11dax__dbz=
     1   ddaux_y11dax__dbc1_x*dbc1_x__dbz
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dbz
     1 + ddaux_y11dax__dbc1_z*dbc1_z__dbz
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dbz
      ddaux_y11dax__dcx=
     1   ddaux_y11dax__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dcx
     1 + ddaux_y11dax__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dcx
      ddaux_y11dax__dcy=
     1   ddaux_y11dax__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dcy
     1 + ddaux_y11dax__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dcy
      ddaux_y11dax__dcz=
     1   ddaux_y11dax__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dax__ddabc1_z11dax*ddabc1_z11dax__dcz
     1 + ddaux_y11dax__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dax__ddabc1_x11dax*ddabc1_x11dax__dcz
c daux_y__day= bc1_x*dabc1_z__day - bc1_z*dabc1_x__day
      ddaux_y11day__day=
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__day
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__day
      ddaux_y11day__daz=
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__daz
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__daz
      ddaux_y11day__dbx=
     1   ddaux_y11day__dbc1_x*dbc1_x__dbx
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__dbx
     1 + ddaux_y11day__dbc1_z*dbc1_z__dbx
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__dbx
      ddaux_y11day__dby=
     1   ddaux_y11day__dbc1_x*dbc1_x__dby
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__dby
     1 + ddaux_y11day__dbc1_z*dbc1_z__dby
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__dby
      ddaux_y11day__dbz=
     1   ddaux_y11day__dbc1_x*dbc1_x__dbz
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__dbz
     1 + ddaux_y11day__dbc1_z*dbc1_z__dbz
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__dbz
      ddaux_y11day__dcx=
     1   ddaux_y11day__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__dcx
     1 + ddaux_y11day__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__dcx
      ddaux_y11day__dcy=
     1   ddaux_y11day__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__dcy
     1 + ddaux_y11day__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__dcy
      ddaux_y11day__dcz=
     1   ddaux_y11day__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11day__ddabc1_z11day*ddabc1_z11day__dcz
     1 + ddaux_y11day__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11day__ddabc1_x11day*ddabc1_x11day__dcz
c daux_y__daz= bc1_x*dabc1_z__daz - bc1_z*dabc1_x__daz
      ddaux_y11daz__daz=
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__daz
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__daz
      ddaux_y11daz__dbx=
     1   ddaux_y11daz__dbc1_x*dbc1_x__dbx
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__dbx
     1 + ddaux_y11daz__dbc1_z*dbc1_z__dbx
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__dbx
      ddaux_y11daz__dby=
     1   ddaux_y11daz__dbc1_x*dbc1_x__dby
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__dby
     1 + ddaux_y11daz__dbc1_z*dbc1_z__dby
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__dby
      ddaux_y11daz__dbz=
     1   ddaux_y11daz__dbc1_x*dbc1_x__dbz
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__dbz
     1 + ddaux_y11daz__dbc1_z*dbc1_z__dbz
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__dbz
      ddaux_y11daz__dcx=
     1   ddaux_y11daz__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__dcx
     1 + ddaux_y11daz__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__dcx
      ddaux_y11daz__dcy=
     1   ddaux_y11daz__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__dcy
     1 + ddaux_y11daz__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__dcy
      ddaux_y11daz__dcz=
     1   ddaux_y11daz__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11daz__ddabc1_z11daz*ddabc1_z11daz__dcz
     1 + ddaux_y11daz__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11daz__ddabc1_x11daz*ddabc1_x11daz__dcz
c daux_y__dbx= bc1_x*dabc1_z__dbx + abc1_z*dbc1_x__dbx - abc1_x*dbc1_z__dbx - bc1_z*dabc1_x__dbx
      ddaux_y11dbx__dbx=
     1   ddaux_y11dbx__dbc1_x*dbc1_x__dbx
     1 + ddaux_y11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dbx
     1 + ddaux_y11dbx__dabc1_z*dabc1_z__dbx
     1 + ddaux_y11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dbx
     1 + ddaux_y11dbx__dabc1_x*dabc1_x__dbx
     1 + ddaux_y11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dbx
     1 + ddaux_y11dbx__dbc1_z*dbc1_z__dbx
     1 + ddaux_y11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dbx
      ddaux_y11dbx__dby=
     1   ddaux_y11dbx__dbc1_x*dbc1_x__dby
     1 + ddaux_y11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dby
     1 + ddaux_y11dbx__dabc1_z*dabc1_z__dby
     1 + ddaux_y11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dby
     1 + ddaux_y11dbx__dabc1_x*dabc1_x__dby
     1 + ddaux_y11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dby
     1 + ddaux_y11dbx__dbc1_z*dbc1_z__dby
     1 + ddaux_y11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dby
      ddaux_y11dbx__dbz=
     1   ddaux_y11dbx__dbc1_x*dbc1_x__dbz
     1 + ddaux_y11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dbz
     1 + ddaux_y11dbx__dabc1_z*dabc1_z__dbz
     1 + ddaux_y11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dbz
     1 + ddaux_y11dbx__dabc1_x*dabc1_x__dbz
     1 + ddaux_y11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dbz
     1 + ddaux_y11dbx__dbc1_z*dbc1_z__dbz
     1 + ddaux_y11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dbz
      ddaux_y11dbx__dcx=
     1   ddaux_y11dbx__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcx
     1 + ddaux_y11dbx__dabc1_z*dabc1_z__dcx
     1 + ddaux_y11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dcx
     1 + ddaux_y11dbx__dabc1_x*dabc1_x__dcx
     1 + ddaux_y11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dcx
     1 + ddaux_y11dbx__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcx
      ddaux_y11dbx__dcy=
     1   ddaux_y11dbx__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcy
     1 + ddaux_y11dbx__dabc1_z*dabc1_z__dcy
     1 + ddaux_y11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dcy
     1 + ddaux_y11dbx__dabc1_x*dabc1_x__dcy
     1 + ddaux_y11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dcy
     1 + ddaux_y11dbx__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcy
      ddaux_y11dbx__dcz=
     1   ddaux_y11dbx__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcz
     1 + ddaux_y11dbx__dabc1_z*dabc1_z__dcz
     1 + ddaux_y11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dcz
     1 + ddaux_y11dbx__dabc1_x*dabc1_x__dcz
     1 + ddaux_y11dbx__ddbc1_z11dbx*ddbc1_z11dbx__dcz
     1 + ddaux_y11dbx__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcz
c daux_y__dby= bc1_x*dabc1_z__dby + abc1_z*dbc1_x__dby - abc1_x*dbc1_z__dby - bc1_z*dabc1_x__dby
      ddaux_y11dby__dby=
     1   ddaux_y11dby__dbc1_x*dbc1_x__dby
     1 + ddaux_y11dby__ddabc1_z11dby*ddabc1_z11dby__dby
     1 + ddaux_y11dby__dabc1_z*dabc1_z__dby
     1 + ddaux_y11dby__ddbc1_x11dby*ddbc1_x11dby__dby
     1 + ddaux_y11dby__dabc1_x*dabc1_x__dby
     1 + ddaux_y11dby__ddbc1_z11dby*ddbc1_z11dby__dby
     1 + ddaux_y11dby__dbc1_z*dbc1_z__dby
     1 + ddaux_y11dby__ddabc1_x11dby*ddabc1_x11dby__dby
      ddaux_y11dby__dbz=
     1   ddaux_y11dby__dbc1_x*dbc1_x__dbz
     1 + ddaux_y11dby__ddabc1_z11dby*ddabc1_z11dby__dbz
     1 + ddaux_y11dby__dabc1_z*dabc1_z__dbz
     1 + ddaux_y11dby__ddbc1_x11dby*ddbc1_x11dby__dbz
     1 + ddaux_y11dby__dabc1_x*dabc1_x__dbz
     1 + ddaux_y11dby__ddbc1_z11dby*ddbc1_z11dby__dbz
     1 + ddaux_y11dby__dbc1_z*dbc1_z__dbz
     1 + ddaux_y11dby__ddabc1_x11dby*ddabc1_x11dby__dbz
      ddaux_y11dby__dcx=
     1   ddaux_y11dby__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11dby__ddabc1_z11dby*ddabc1_z11dby__dcx
     1 + ddaux_y11dby__dabc1_z*dabc1_z__dcx
     1 + ddaux_y11dby__ddbc1_x11dby*ddbc1_x11dby__dcx
     1 + ddaux_y11dby__dabc1_x*dabc1_x__dcx
     1 + ddaux_y11dby__ddbc1_z11dby*ddbc1_z11dby__dcx
     1 + ddaux_y11dby__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11dby__ddabc1_x11dby*ddabc1_x11dby__dcx
      ddaux_y11dby__dcy=
     1   ddaux_y11dby__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11dby__ddabc1_z11dby*ddabc1_z11dby__dcy
     1 + ddaux_y11dby__dabc1_z*dabc1_z__dcy
     1 + ddaux_y11dby__ddbc1_x11dby*ddbc1_x11dby__dcy
     1 + ddaux_y11dby__dabc1_x*dabc1_x__dcy
     1 + ddaux_y11dby__ddbc1_z11dby*ddbc1_z11dby__dcy
     1 + ddaux_y11dby__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11dby__ddabc1_x11dby*ddabc1_x11dby__dcy
      ddaux_y11dby__dcz=
     1   ddaux_y11dby__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dby__ddabc1_z11dby*ddabc1_z11dby__dcz
     1 + ddaux_y11dby__dabc1_z*dabc1_z__dcz
     1 + ddaux_y11dby__ddbc1_x11dby*ddbc1_x11dby__dcz
     1 + ddaux_y11dby__dabc1_x*dabc1_x__dcz
     1 + ddaux_y11dby__ddbc1_z11dby*ddbc1_z11dby__dcz
     1 + ddaux_y11dby__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dby__ddabc1_x11dby*ddabc1_x11dby__dcz
c daux_y__dbz= bc1_x*dabc1_z__dbz + abc1_z*dbc1_x__dbz - abc1_x*dbc1_z__dbz - bc1_z*dabc1_x__dbz
      ddaux_y11dbz__dbz=
     1   ddaux_y11dbz__dbc1_x*dbc1_x__dbz
     1 + ddaux_y11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dbz
     1 + ddaux_y11dbz__dabc1_z*dabc1_z__dbz
     1 + ddaux_y11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dbz
     1 + ddaux_y11dbz__dabc1_x*dabc1_x__dbz
     1 + ddaux_y11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dbz
     1 + ddaux_y11dbz__dbc1_z*dbc1_z__dbz
     1 + ddaux_y11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dbz
      ddaux_y11dbz__dcx=
     1   ddaux_y11dbz__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcx
     1 + ddaux_y11dbz__dabc1_z*dabc1_z__dcx
     1 + ddaux_y11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dcx
     1 + ddaux_y11dbz__dabc1_x*dabc1_x__dcx
     1 + ddaux_y11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dcx
     1 + ddaux_y11dbz__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcx
      ddaux_y11dbz__dcy=
     1   ddaux_y11dbz__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcy
     1 + ddaux_y11dbz__dabc1_z*dabc1_z__dcy
     1 + ddaux_y11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dcy
     1 + ddaux_y11dbz__dabc1_x*dabc1_x__dcy
     1 + ddaux_y11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dcy
     1 + ddaux_y11dbz__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcy
      ddaux_y11dbz__dcz=
     1   ddaux_y11dbz__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcz
     1 + ddaux_y11dbz__dabc1_z*dabc1_z__dcz
     1 + ddaux_y11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dcz
     1 + ddaux_y11dbz__dabc1_x*dabc1_x__dcz
     1 + ddaux_y11dbz__ddbc1_z11dbz*ddbc1_z11dbz__dcz
     1 + ddaux_y11dbz__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcz
c daux_y__dcx= bc1_x*dabc1_z__dcx + abc1_z*dbc1_x__dcx - abc1_x*dbc1_z__dcx - bc1_z*dabc1_x__dcx
      ddaux_y11dcx__dcx=
     1   ddaux_y11dcx__dbc1_x*dbc1_x__dcx
     1 + ddaux_y11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcx
     1 + ddaux_y11dcx__dabc1_z*dabc1_z__dcx
     1 + ddaux_y11dcx__ddbc1_x11dcx*ddbc1_x11dcx__dcx
     1 + ddaux_y11dcx__dabc1_x*dabc1_x__dcx
     1 + ddaux_y11dcx__ddbc1_z11dcx*ddbc1_z11dcx__dcx
     1 + ddaux_y11dcx__dbc1_z*dbc1_z__dcx
     1 + ddaux_y11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcx
      ddaux_y11dcx__dcy=
     1   ddaux_y11dcx__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcy
     1 + ddaux_y11dcx__dabc1_z*dabc1_z__dcy
     1 + ddaux_y11dcx__ddbc1_x11dcx*ddbc1_x11dcx__dcy
     1 + ddaux_y11dcx__dabc1_x*dabc1_x__dcy
     1 + ddaux_y11dcx__ddbc1_z11dcx*ddbc1_z11dcx__dcy
     1 + ddaux_y11dcx__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcy
      ddaux_y11dcx__dcz=
     1   ddaux_y11dcx__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcz
     1 + ddaux_y11dcx__dabc1_z*dabc1_z__dcz
     1 + ddaux_y11dcx__ddbc1_x11dcx*ddbc1_x11dcx__dcz
     1 + ddaux_y11dcx__dabc1_x*dabc1_x__dcz
     1 + ddaux_y11dcx__ddbc1_z11dcx*ddbc1_z11dcx__dcz
     1 + ddaux_y11dcx__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcz
c daux_y__dcy= bc1_x*dabc1_z__dcy + abc1_z*dbc1_x__dcy - abc1_x*dbc1_z__dcy - bc1_z*dabc1_x__dcy
      ddaux_y11dcy__dcy=
     1   ddaux_y11dcy__dbc1_x*dbc1_x__dcy
     1 + ddaux_y11dcy__ddabc1_z11dcy*ddabc1_z11dcy__dcy
     1 + ddaux_y11dcy__dabc1_z*dabc1_z__dcy
     1 + ddaux_y11dcy__ddbc1_x11dcy*ddbc1_x11dcy__dcy
     1 + ddaux_y11dcy__dabc1_x*dabc1_x__dcy
     1 + ddaux_y11dcy__ddbc1_z11dcy*ddbc1_z11dcy__dcy
     1 + ddaux_y11dcy__dbc1_z*dbc1_z__dcy
     1 + ddaux_y11dcy__ddabc1_x11dcy*ddabc1_x11dcy__dcy
      ddaux_y11dcy__dcz=
     1   ddaux_y11dcy__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dcy__ddabc1_z11dcy*ddabc1_z11dcy__dcz
     1 + ddaux_y11dcy__dabc1_z*dabc1_z__dcz
     1 + ddaux_y11dcy__ddbc1_x11dcy*ddbc1_x11dcy__dcz
     1 + ddaux_y11dcy__dabc1_x*dabc1_x__dcz
     1 + ddaux_y11dcy__ddbc1_z11dcy*ddbc1_z11dcy__dcz
     1 + ddaux_y11dcy__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dcy__ddabc1_x11dcy*ddabc1_x11dcy__dcz
c daux_y__dcz= bc1_x*dabc1_z__dcz + abc1_z*dbc1_x__dcz - abc1_x*dbc1_z__dcz - bc1_z*dabc1_x__dcz
      ddaux_y11dcz__dcz=
     1   ddaux_y11dcz__dbc1_x*dbc1_x__dcz
     1 + ddaux_y11dcz__ddabc1_z11dcz*ddabc1_z11dcz__dcz
     1 + ddaux_y11dcz__dabc1_z*dabc1_z__dcz
     1 + ddaux_y11dcz__ddbc1_x11dcz*ddbc1_x11dcz__dcz
     1 + ddaux_y11dcz__dabc1_x*dabc1_x__dcz
     1 + ddaux_y11dcz__ddbc1_z11dcz*ddbc1_z11dcz__dcz
     1 + ddaux_y11dcz__dbc1_z*dbc1_z__dcz
     1 + ddaux_y11dcz__ddabc1_x11dcz*ddabc1_x11dcz__dcz


c daux_z__dax= bc1_y*dabc1_x__dax - bc1_x*dabc1_y__dax
      ddaux_z11dax__dax=
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dax
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dax
      ddaux_z11dax__day=
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__day
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__day
      ddaux_z11dax__daz=
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__daz
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__daz
      ddaux_z11dax__dbx=
     1   ddaux_z11dax__dbc1_y*dbc1_y__dbx
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dbx
     1 + ddaux_z11dax__dbc1_x*dbc1_x__dbx
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dbx
      ddaux_z11dax__dby=
     1   ddaux_z11dax__dbc1_y*dbc1_y__dby
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dby
     1 + ddaux_z11dax__dbc1_x*dbc1_x__dby
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dby
      ddaux_z11dax__dbz=
     1   ddaux_z11dax__dbc1_y*dbc1_y__dbz
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dbz
     1 + ddaux_z11dax__dbc1_x*dbc1_x__dbz
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dbz
      ddaux_z11dax__dcx=
     1   ddaux_z11dax__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dcx
     1 + ddaux_z11dax__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dcx
      ddaux_z11dax__dcy=
     1   ddaux_z11dax__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dcy
     1 + ddaux_z11dax__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dcy
      ddaux_z11dax__dcz=
     1   ddaux_z11dax__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dax__ddabc1_x11dax*ddabc1_x11dax__dcz
     1 + ddaux_z11dax__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dax__ddabc1_y11dax*ddabc1_y11dax__dcz
c daux_z__day= bc1_y*dabc1_x__day - bc1_x*dabc1_y__day
      ddaux_z11day__day=
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__day
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__day
      ddaux_z11day__daz=
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__daz
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__daz
      ddaux_z11day__dbx=
     1   ddaux_z11day__dbc1_y*dbc1_y__dbx
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__dbx
     1 + ddaux_z11day__dbc1_x*dbc1_x__dbx
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__dbx
      ddaux_z11day__dby=
     1   ddaux_z11day__dbc1_y*dbc1_y__dby
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__dby
     1 + ddaux_z11day__dbc1_x*dbc1_x__dby
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__dby
      ddaux_z11day__dbz=
     1   ddaux_z11day__dbc1_y*dbc1_y__dbz
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__dbz
     1 + ddaux_z11day__dbc1_x*dbc1_x__dbz
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__dbz
      ddaux_z11day__dcx=
     1   ddaux_z11day__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__dcx
     1 + ddaux_z11day__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__dcx
      ddaux_z11day__dcy=
     1   ddaux_z11day__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__dcy
     1 + ddaux_z11day__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__dcy
      ddaux_z11day__dcz=
     1   ddaux_z11day__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11day__ddabc1_x11day*ddabc1_x11day__dcz
     1 + ddaux_z11day__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11day__ddabc1_y11day*ddabc1_y11day__dcz
c daux_z__daz= bc1_y*dabc1_x__daz - bc1_x*dabc1_y__daz
      ddaux_z11daz__daz=
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__daz
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__daz
      ddaux_z11daz__dbx=
     1   ddaux_z11daz__dbc1_y*dbc1_y__dbx
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__dbx
     1 + ddaux_z11daz__dbc1_x*dbc1_x__dbx
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__dbx
      ddaux_z11daz__dby=
     1   ddaux_z11daz__dbc1_y*dbc1_y__dby
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__dby
     1 + ddaux_z11daz__dbc1_x*dbc1_x__dby
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__dby
      ddaux_z11daz__dbz=
     1   ddaux_z11daz__dbc1_y*dbc1_y__dbz
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__dbz
     1 + ddaux_z11daz__dbc1_x*dbc1_x__dbz
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__dbz
      ddaux_z11daz__dcx=
     1   ddaux_z11daz__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__dcx
     1 + ddaux_z11daz__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__dcx
      ddaux_z11daz__dcy=
     1   ddaux_z11daz__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__dcy
     1 + ddaux_z11daz__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__dcy
      ddaux_z11daz__dcz=
     1   ddaux_z11daz__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11daz__ddabc1_x11daz*ddabc1_x11daz__dcz
     1 + ddaux_z11daz__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11daz__ddabc1_y11daz*ddabc1_y11daz__dcz
c daux_z__dbx= bc1_y*dabc1_x__dbx + abc1_x*dbc1_y__dbx - abc1_y*dbc1_x__dbx - bc1_x*dabc1_y__dbx
      ddaux_z11dbx__dbx=
     1   ddaux_z11dbx__dbc1_y*dbc1_y__dbx
     1 + ddaux_z11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dbx
     1 + ddaux_z11dbx__dabc1_x*dabc1_x__dbx
     1 + ddaux_z11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dbx
     1 + ddaux_z11dbx__dabc1_y*dabc1_y__dbx
     1 + ddaux_z11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dbx
     1 + ddaux_z11dbx__dbc1_x*dbc1_x__dbx
     1 + ddaux_z11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dbx
      ddaux_z11dbx__dby=
     1   ddaux_z11dbx__dbc1_y*dbc1_y__dby
     1 + ddaux_z11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dby
     1 + ddaux_z11dbx__dabc1_x*dabc1_x__dby
     1 + ddaux_z11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dby
     1 + ddaux_z11dbx__dabc1_y*dabc1_y__dby
     1 + ddaux_z11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dby
     1 + ddaux_z11dbx__dbc1_x*dbc1_x__dby
     1 + ddaux_z11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dby
      ddaux_z11dbx__dbz=
     1   ddaux_z11dbx__dbc1_y*dbc1_y__dbz
     1 + ddaux_z11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dbz
     1 + ddaux_z11dbx__dabc1_x*dabc1_x__dbz
     1 + ddaux_z11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dbz
     1 + ddaux_z11dbx__dabc1_y*dabc1_y__dbz
     1 + ddaux_z11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dbz
     1 + ddaux_z11dbx__dbc1_x*dbc1_x__dbz
     1 + ddaux_z11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dbz
      ddaux_z11dbx__dcx=
     1   ddaux_z11dbx__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcx
     1 + ddaux_z11dbx__dabc1_x*dabc1_x__dcx
     1 + ddaux_z11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dcx
     1 + ddaux_z11dbx__dabc1_y*dabc1_y__dcx
     1 + ddaux_z11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dcx
     1 + ddaux_z11dbx__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcx
      ddaux_z11dbx__dcy=
     1   ddaux_z11dbx__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcy
     1 + ddaux_z11dbx__dabc1_x*dabc1_x__dcy
     1 + ddaux_z11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dcy
     1 + ddaux_z11dbx__dabc1_y*dabc1_y__dcy
     1 + ddaux_z11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dcy
     1 + ddaux_z11dbx__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcy
      ddaux_z11dbx__dcz=
     1   ddaux_z11dbx__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcz
     1 + ddaux_z11dbx__dabc1_x*dabc1_x__dcz
     1 + ddaux_z11dbx__ddbc1_y11dbx*ddbc1_y11dbx__dcz
     1 + ddaux_z11dbx__dabc1_y*dabc1_y__dcz
     1 + ddaux_z11dbx__ddbc1_x11dbx*ddbc1_x11dbx__dcz
     1 + ddaux_z11dbx__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcz
c daux_z__dby= bc1_y*dabc1_x__dby + abc1_x*dbc1_y__dby - abc1_y*dbc1_x__dby - bc1_x*dabc1_y__dby
      ddaux_z11dby__dby=
     1   ddaux_z11dby__dbc1_y*dbc1_y__dby
     1 + ddaux_z11dby__ddabc1_x11dby*ddabc1_x11dby__dby
     1 + ddaux_z11dby__dabc1_x*dabc1_x__dby
     1 + ddaux_z11dby__ddbc1_y11dby*ddbc1_y11dby__dby
     1 + ddaux_z11dby__dabc1_y*dabc1_y__dby
     1 + ddaux_z11dby__ddbc1_x11dby*ddbc1_x11dby__dby
     1 + ddaux_z11dby__dbc1_x*dbc1_x__dby
     1 + ddaux_z11dby__ddabc1_y11dby*ddabc1_y11dby__dby
      ddaux_z11dby__dbz=
     1   ddaux_z11dby__dbc1_y*dbc1_y__dbz
     1 + ddaux_z11dby__ddabc1_x11dby*ddabc1_x11dby__dbz
     1 + ddaux_z11dby__dabc1_x*dabc1_x__dbz
     1 + ddaux_z11dby__ddbc1_y11dby*ddbc1_y11dby__dbz
     1 + ddaux_z11dby__dabc1_y*dabc1_y__dbz
     1 + ddaux_z11dby__ddbc1_x11dby*ddbc1_x11dby__dbz
     1 + ddaux_z11dby__dbc1_x*dbc1_x__dbz
     1 + ddaux_z11dby__ddabc1_y11dby*ddabc1_y11dby__dbz
      ddaux_z11dby__dcx=
     1   ddaux_z11dby__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11dby__ddabc1_x11dby*ddabc1_x11dby__dcx
     1 + ddaux_z11dby__dabc1_x*dabc1_x__dcx
     1 + ddaux_z11dby__ddbc1_y11dby*ddbc1_y11dby__dcx
     1 + ddaux_z11dby__dabc1_y*dabc1_y__dcx
     1 + ddaux_z11dby__ddbc1_x11dby*ddbc1_x11dby__dcx
     1 + ddaux_z11dby__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11dby__ddabc1_y11dby*ddabc1_y11dby__dcx
      ddaux_z11dby__dcy=
     1   ddaux_z11dby__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11dby__ddabc1_x11dby*ddabc1_x11dby__dcy
     1 + ddaux_z11dby__dabc1_x*dabc1_x__dcy
     1 + ddaux_z11dby__ddbc1_y11dby*ddbc1_y11dby__dcy
     1 + ddaux_z11dby__dabc1_y*dabc1_y__dcy
     1 + ddaux_z11dby__ddbc1_x11dby*ddbc1_x11dby__dcy
     1 + ddaux_z11dby__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11dby__ddabc1_y11dby*ddabc1_y11dby__dcy
      ddaux_z11dby__dcz=
     1   ddaux_z11dby__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dby__ddabc1_x11dby*ddabc1_x11dby__dcz
     1 + ddaux_z11dby__dabc1_x*dabc1_x__dcz
     1 + ddaux_z11dby__ddbc1_y11dby*ddbc1_y11dby__dcz
     1 + ddaux_z11dby__dabc1_y*dabc1_y__dcz
     1 + ddaux_z11dby__ddbc1_x11dby*ddbc1_x11dby__dcz
     1 + ddaux_z11dby__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dby__ddabc1_y11dby*ddabc1_y11dby__dcz
c daux_z__dbz= bc1_y*dabc1_x__dbz + abc1_x*dbc1_y__dbz - abc1_y*dbc1_x__dbz - bc1_x*dabc1_y__dbz
      ddaux_z11dbz__dbz=
     1   ddaux_z11dbz__dbc1_y*dbc1_y__dbz
     1 + ddaux_z11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dbz
     1 + ddaux_z11dbz__dabc1_x*dabc1_x__dbz
     1 + ddaux_z11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dbz
     1 + ddaux_z11dbz__dabc1_y*dabc1_y__dbz
     1 + ddaux_z11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dbz
     1 + ddaux_z11dbz__dbc1_x*dbc1_x__dbz
     1 + ddaux_z11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dbz
      ddaux_z11dbz__dcx=
     1   ddaux_z11dbz__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcx
     1 + ddaux_z11dbz__dabc1_x*dabc1_x__dcx
     1 + ddaux_z11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dcx
     1 + ddaux_z11dbz__dabc1_y*dabc1_y__dcx
     1 + ddaux_z11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dcx
     1 + ddaux_z11dbz__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcx
      ddaux_z11dbz__dcy=
     1   ddaux_z11dbz__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcy
     1 + ddaux_z11dbz__dabc1_x*dabc1_x__dcy
     1 + ddaux_z11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dcy
     1 + ddaux_z11dbz__dabc1_y*dabc1_y__dcy
     1 + ddaux_z11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dcy
     1 + ddaux_z11dbz__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcy
      ddaux_z11dbz__dcz=
     1   ddaux_z11dbz__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcz
     1 + ddaux_z11dbz__dabc1_x*dabc1_x__dcz
     1 + ddaux_z11dbz__ddbc1_y11dbz*ddbc1_y11dbz__dcz
     1 + ddaux_z11dbz__dabc1_y*dabc1_y__dcz
     1 + ddaux_z11dbz__ddbc1_x11dbz*ddbc1_x11dbz__dcz
     1 + ddaux_z11dbz__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcz
c daux_z__dcx= bc1_y*dabc1_x__dcx + abc1_x*dbc1_y__dcx - abc1_y*dbc1_x__dcx - bc1_x*dabc1_y__dcx
      ddaux_z11dcx__dcx=
     1   ddaux_z11dcx__dbc1_y*dbc1_y__dcx
     1 + ddaux_z11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcx
     1 + ddaux_z11dcx__dabc1_x*dabc1_x__dcx
     1 + ddaux_z11dcx__ddbc1_y11dcx*ddbc1_y11dcx__dcx
     1 + ddaux_z11dcx__dabc1_y*dabc1_y__dcx
     1 + ddaux_z11dcx__ddbc1_x11dcx*ddbc1_x11dcx__dcx
     1 + ddaux_z11dcx__dbc1_x*dbc1_x__dcx
     1 + ddaux_z11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcx
      ddaux_z11dcx__dcy=
     1   ddaux_z11dcx__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcy
     1 + ddaux_z11dcx__dabc1_x*dabc1_x__dcy
     1 + ddaux_z11dcx__ddbc1_y11dcx*ddbc1_y11dcx__dcy
     1 + ddaux_z11dcx__dabc1_y*dabc1_y__dcy
     1 + ddaux_z11dcx__ddbc1_x11dcx*ddbc1_x11dcx__dcy
     1 + ddaux_z11dcx__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcy
      ddaux_z11dcx__dcz=
     1   ddaux_z11dcx__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcz
     1 + ddaux_z11dcx__dabc1_x*dabc1_x__dcz
     1 + ddaux_z11dcx__ddbc1_y11dcx*ddbc1_y11dcx__dcz
     1 + ddaux_z11dcx__dabc1_y*dabc1_y__dcz
     1 + ddaux_z11dcx__ddbc1_x11dcx*ddbc1_x11dcx__dcz
     1 + ddaux_z11dcx__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcz
c daux_z__dcy= bc1_y*dabc1_x__dcy + abc1_x*dbc1_y__dcy - abc1_y*dbc1_x__dcy - bc1_x*dabc1_y__dcy
      ddaux_z11dcy__dcy=
     1   ddaux_z11dcy__dbc1_y*dbc1_y__dcy
     1 + ddaux_z11dcy__ddabc1_x11dcy*ddabc1_x11dcy__dcy
     1 + ddaux_z11dcy__dabc1_x*dabc1_x__dcy
     1 + ddaux_z11dcy__ddbc1_y11dcy*ddbc1_y11dcy__dcy
     1 + ddaux_z11dcy__dabc1_y*dabc1_y__dcy
     1 + ddaux_z11dcy__ddbc1_x11dcy*ddbc1_x11dcy__dcy
     1 + ddaux_z11dcy__dbc1_x*dbc1_x__dcy
     1 + ddaux_z11dcy__ddabc1_y11dcy*ddabc1_y11dcy__dcy
      ddaux_z11dcy__dcz=
     1   ddaux_z11dcy__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dcy__ddabc1_x11dcy*ddabc1_x11dcy__dcz
     1 + ddaux_z11dcy__dabc1_x*dabc1_x__dcz
     1 + ddaux_z11dcy__ddbc1_y11dcy*ddbc1_y11dcy__dcz
     1 + ddaux_z11dcy__dabc1_y*dabc1_y__dcz
     1 + ddaux_z11dcy__ddbc1_x11dcy*ddbc1_x11dcy__dcz
     1 + ddaux_z11dcy__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dcy__ddabc1_y11dcy*ddabc1_y11dcy__dcz
c daux_z__dcz= bc1_y*dabc1_x__dcz + abc1_x*dbc1_y__dcz - abc1_y*dbc1_x__dcz - bc1_x*dabc1_y__dcz
      ddaux_z11dcz__dcz=
     1   ddaux_z11dcz__dbc1_y*dbc1_y__dcz
     1 + ddaux_z11dcz__ddabc1_x11dcz*ddabc1_x11dcz__dcz
     1 + ddaux_z11dcz__dabc1_x*dabc1_x__dcz
     1 + ddaux_z11dcz__ddbc1_y11dcz*ddbc1_y11dcz__dcz
     1 + ddaux_z11dcz__dabc1_y*dabc1_y__dcz
     1 + ddaux_z11dcz__ddbc1_x11dcz*ddbc1_x11dcz__dcz
     1 + ddaux_z11dcz__dbc1_x*dbc1_x__dcz
     1 + ddaux_z11dcz__ddabc1_y11dcz*ddabc1_y11dcz__dcz
 
 
c dbcd1_x__dbx= bcd_x*dbcd_length_inv__dbx
      ddbcd1_x11dbx__dbx=
     1  ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dbx
      ddbcd1_x11dbx__dby=
     1   ddbcd1_x11dbx__dbcd_x*dbcd_x__dby
     1 + ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dby
      ddbcd1_x11dbx__dbz=
     1   ddbcd1_x11dbx__dbcd_x*dbcd_x__dbz
     1 + ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dbz
      ddbcd1_x11dbx__dcx=
     1  ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcx
      ddbcd1_x11dbx__dcy=
     1   ddbcd1_x11dbx__dbcd_x*dbcd_x__dcy
     1 + ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcy
      ddbcd1_x11dbx__dcz=
     1   ddbcd1_x11dbx__dbcd_x*dbcd_x__dcz
     1 + ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcz
      ddbcd1_x11dbx__ddx=
     1  ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddx
      ddbcd1_x11dbx__ddy=
     1   ddbcd1_x11dbx__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddy
      ddbcd1_x11dbx__ddz=
     1   ddbcd1_x11dbx__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddz
c dbcd1_x__dby=-bcd_length_inv*cd_z + bcd_x*dbcd_length_inv__dby
      ddbcd1_x11dby__dby=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__dby
     1 + ddbcd1_x11dby__dbcd_x*dbcd_x__dby
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dby
      ddbcd1_x11dby__dbz=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__dbz
     1 + ddbcd1_x11dby__dbcd_x*dbcd_x__dbz
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dbz
      ddbcd1_x11dby__dcx=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcx
      ddbcd1_x11dby__dcy=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_x11dby__dbcd_x*dbcd_x__dcy
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcy
      ddbcd1_x11dby__dcz=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_x11dby__dcd_z*dcd_z__dcz
     1 + ddbcd1_x11dby__dbcd_x*dbcd_x__dcz
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcz
      ddbcd1_x11dby__ddx=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddx
      ddbcd1_x11dby__ddy=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_x11dby__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddy
      ddbcd1_x11dby__ddz=
     1   ddbcd1_x11dby__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_x11dby__dcd_z*dcd_z__ddz
     1 + ddbcd1_x11dby__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddz
c dbcd1_x__dbz=bcd_length_inv*cd_y + bcd_x*dbcd_length_inv__dbz
      ddbcd1_x11dbz__dbz=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__dbz
     1 + ddbcd1_x11dbz__dbcd_x*dbcd_x__dbz
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dbz
      ddbcd1_x11dbz__dcx=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcx
      ddbcd1_x11dbz__dcy=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_x11dbz__dcd_y*dcd_y__dcy
     1 + ddbcd1_x11dbz__dbcd_x*dbcd_x__dcy
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcy
      ddbcd1_x11dbz__dcz=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_x11dbz__dbcd_x*dbcd_x__dcz
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcz
      ddbcd1_x11dbz__ddx=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddx
      ddbcd1_x11dbz__ddy=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_x11dbz__dcd_y*dcd_y__ddy
     1 + ddbcd1_x11dbz__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddy
      ddbcd1_x11dbz__ddz=
     1   ddbcd1_x11dbz__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_x11dbz__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddz
c dbcd1_x__dcx= bcd_x*dbcd_length_inv__dcx
      ddbcd1_x11dcx__dcx=
     1  ddbcd1_x11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcx
      ddbcd1_x11dcx__dcy=
     1   ddbcd1_x11dcx__dbcd_x*dbcd_x__dcy
     1 + ddbcd1_x11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcy
      ddbcd1_x11dcx__dcz=
     1   ddbcd1_x11dcx__dbcd_x*dbcd_x__dcz
     1 + ddbcd1_x11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcz
      ddbcd1_x11dcx__ddx=
     1  ddbcd1_x11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddx
      ddbcd1_x11dcx__ddy=
     1   ddbcd1_x11dcx__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddy
      ddbcd1_x11dcx__ddz=
     1   ddbcd1_x11dcx__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddz
c dbcd1_x__dcy=bcd_length_inv*dbcd_x__dcy + bcd_x*dbcd_length_inv__dcy
      ddbcd1_x11dcy__dcy=
     1   ddbcd1_x11dcy__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_x11dcy__dbcd_x*dbcd_x__dcy
     1 + ddbcd1_x11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__dcy
      ddbcd1_x11dcy__dcz=
     1   ddbcd1_x11dcy__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_x11dcy__ddbcd_x11dcy*ddbcd_x11dcy__dcz
     1 + ddbcd1_x11dcy__dbcd_x*dbcd_x__dcz
     1 + ddbcd1_x11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__dcz
      ddbcd1_x11dcy__ddx=
     1   ddbcd1_x11dcy__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_x11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddx
      ddbcd1_x11dcy__ddy=
     1   ddbcd1_x11dcy__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_x11dcy__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddy
      ddbcd1_x11dcy__ddz=
     1   ddbcd1_x11dcy__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_x11dcy__ddbcd_x11dcy*ddbcd_x11dcy__ddz
     1 + ddbcd1_x11dcy__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddz
c dbcd1_x__dcz=bcd_length_inv*dbcd_x__dcz + bcd_x*dbcd_length_inv__dcz
      ddbcd1_x11dcz__dcz=
     1   ddbcd1_x11dcz__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_x11dcz__dbcd_x*dbcd_x__dcz
     1 + ddbcd1_x11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__dcz
      ddbcd1_x11dcz__ddx=
     1   ddbcd1_x11dcz__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_x11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddx
      ddbcd1_x11dcz__ddy=
     1   ddbcd1_x11dcz__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_x11dcz__ddbcd_x11dcz*ddbcd_x11dcz__ddy
     1 + ddbcd1_x11dcz__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddy
      ddbcd1_x11dcz__ddz=
     1   ddbcd1_x11dcz__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_x11dcz__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddz
c dbcd1_x__ddx= bcd_x*dbcd_length_inv__ddx
      ddbcd1_x11ddx__ddx=
     1  ddbcd1_x11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddx
      ddbcd1_x11ddx__ddy=
     1   ddbcd1_x11ddx__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddy
      ddbcd1_x11ddx__ddz=
     1   ddbcd1_x11ddx__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddz
c dbcd1_x__ddy=-bcd_length_inv*bc_z + bcd_x*dbcd_length_inv__ddy
      ddbcd1_x11ddy__ddy=
     1   ddbcd1_x11ddy__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_x11ddy__dbcd_x*dbcd_x__ddy
     1 + ddbcd1_x11ddy__ddbcd_length_inv11ddy*ddbcd_length_inv11ddy__ddy
      ddbcd1_x11ddy__ddz=
     1   ddbcd1_x11ddy__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_x11ddy__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11ddy__ddbcd_length_inv11ddy*ddbcd_length_inv11ddy__ddz
c dbcd1_x__ddz=bcd_length_inv*bc_y + bcd_x*dbcd_length_inv__ddz
      ddbcd1_x11ddz__ddz=
     1   ddbcd1_x11ddz__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_x11ddz__dbcd_x*dbcd_x__ddz
     1 + ddbcd1_x11ddz__ddbcd_length_inv11ddz*ddbcd_length_inv11ddz__ddz


c dbcd1_y__dbx=bcd_length_inv*cd_z + bcd_y*dbcd_length_inv__dbx
      ddbcd1_y11dbx__dbx=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__dbx
     1 + ddbcd1_y11dbx__dbcd_y*dbcd_y__dbx
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dbx
      ddbcd1_y11dbx__dby=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__dby
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dby
      ddbcd1_y11dbx__dbz=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__dbz
     1 + ddbcd1_y11dbx__dbcd_y*dbcd_y__dbz
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dbz
      ddbcd1_y11dbx__dcx=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_y11dbx__dbcd_y*dbcd_y__dcx
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcx
      ddbcd1_y11dbx__dcy=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcy
      ddbcd1_y11dbx__dcz=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_y11dbx__dcd_z*dcd_z__dcz
     1 + ddbcd1_y11dbx__dbcd_y*dbcd_y__dcz
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcz
      ddbcd1_y11dbx__ddx=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_y11dbx__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddx
      ddbcd1_y11dbx__ddy=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddy
      ddbcd1_y11dbx__ddz=
     1   ddbcd1_y11dbx__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_y11dbx__dcd_z*dcd_z__ddz
     1 + ddbcd1_y11dbx__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddz
c dbcd1_y__dby= bcd_y*dbcd_length_inv__dby
      ddbcd1_y11dby__dby=
     1  ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dby
      ddbcd1_y11dby__dbz=
     1   ddbcd1_y11dby__dbcd_y*dbcd_y__dbz
     1 + ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dbz
      ddbcd1_y11dby__dcx=
     1   ddbcd1_y11dby__dbcd_y*dbcd_y__dcx
     1 + ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcx
      ddbcd1_y11dby__dcy=
     1  ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcy
      ddbcd1_y11dby__dcz=
     1   ddbcd1_y11dby__dbcd_y*dbcd_y__dcz
     1 + ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcz
      ddbcd1_y11dby__ddx=
     1   ddbcd1_y11dby__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddx
      ddbcd1_y11dby__ddy=
     1  ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddy
      ddbcd1_y11dby__ddz=
     1   ddbcd1_y11dby__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddz
c dbcd1_y__dbz=-bcd_length_inv*cd_x + bcd_y*dbcd_length_inv__dbz
      ddbcd1_y11dbz__dbz=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__dbz
     1 + ddbcd1_y11dbz__dbcd_y*dbcd_y__dbz
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dbz
      ddbcd1_y11dbz__dcx=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_y11dbz__dcd_x*dcd_x__dcx
     1 + ddbcd1_y11dbz__dbcd_y*dbcd_y__dcx
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcx
      ddbcd1_y11dbz__dcy=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcy
      ddbcd1_y11dbz__dcz=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_y11dbz__dbcd_y*dbcd_y__dcz
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcz
      ddbcd1_y11dbz__ddx=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_y11dbz__dcd_x*dcd_x__ddx
     1 + ddbcd1_y11dbz__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddx
      ddbcd1_y11dbz__ddy=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddy
      ddbcd1_y11dbz__ddz=
     1   ddbcd1_y11dbz__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_y11dbz__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddz
c dbcd1_y__dcx=bcd_length_inv*dbcd_y__dcx + bcd_y*dbcd_length_inv__dcx
      ddbcd1_y11dcx__dcx=
     1   ddbcd1_y11dcx__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_y11dcx__dbcd_y*dbcd_y__dcx
     1 + ddbcd1_y11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcx
      ddbcd1_y11dcx__dcy=
     1   ddbcd1_y11dcx__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_y11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcy
      ddbcd1_y11dcx__dcz=
     1   ddbcd1_y11dcx__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_y11dcx__ddbcd_y11dcx*ddbcd_y11dcx__dcz
     1 + ddbcd1_y11dcx__dbcd_y*dbcd_y__dcz
     1 + ddbcd1_y11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcz
      ddbcd1_y11dcx__ddx=
     1   ddbcd1_y11dcx__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_y11dcx__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddx
      ddbcd1_y11dcx__ddy=
     1   ddbcd1_y11dcx__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_y11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddy
      ddbcd1_y11dcx__ddz=
     1   ddbcd1_y11dcx__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_y11dcx__ddbcd_y11dcx*ddbcd_y11dcx__ddz
     1 + ddbcd1_y11dcx__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddz
c dbcd1_y__dcy= bcd_y*dbcd_length_inv__dcy
      ddbcd1_y11dcy__dcy=
     1  ddbcd1_y11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__dcy
      ddbcd1_y11dcy__dcz=
     1   ddbcd1_y11dcy__dbcd_y*dbcd_y__dcz
     1 + ddbcd1_y11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__dcz
      ddbcd1_y11dcy__ddx=
     1   ddbcd1_y11dcy__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddx
      ddbcd1_y11dcy__ddy=
     1  ddbcd1_y11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddy
      ddbcd1_y11dcy__ddz=
     1   ddbcd1_y11dcy__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddz
c dbcd1_y__dcz=bcd_length_inv*dbcd_y__dcz + bcd_y*dbcd_length_inv__dcz
      ddbcd1_y11dcz__dcz=
     1   ddbcd1_y11dcz__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_y11dcz__dbcd_y*dbcd_y__dcz
     1 + ddbcd1_y11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__dcz
      ddbcd1_y11dcz__ddx=
     1   ddbcd1_y11dcz__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_y11dcz__ddbcd_y11dcz*ddbcd_y11dcz__ddx
     1 + ddbcd1_y11dcz__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddx
      ddbcd1_y11dcz__ddy=
     1   ddbcd1_y11dcz__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_y11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddy
      ddbcd1_y11dcz__ddz=
     1   ddbcd1_y11dcz__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_y11dcz__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddz
c dbcd1_y__ddx=bcd_length_inv*bc_z + bcd_y*dbcd_length_inv__ddx
      ddbcd1_y11ddx__ddx=
     1   ddbcd1_y11ddx__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_y11ddx__dbcd_y*dbcd_y__ddx
     1 + ddbcd1_y11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddx
      ddbcd1_y11ddx__ddy=
     1   ddbcd1_y11ddx__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_y11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddy
      ddbcd1_y11ddx__ddz=
     1   ddbcd1_y11ddx__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_y11ddx__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddz
c dbcd1_y__ddy= bcd_y*dbcd_length_inv__ddy
      ddbcd1_y11ddy__ddy=
     1  ddbcd1_y11ddy__ddbcd_length_inv11ddy*ddbcd_length_inv11ddy__ddy
      ddbcd1_y11ddy__ddz=
     1   ddbcd1_y11ddy__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11ddy__ddbcd_length_inv11ddy*ddbcd_length_inv11ddy__ddz
c dbcd1_y__ddz=-bcd_length_inv*bc_x + bcd_y*dbcd_length_inv__ddz
      ddbcd1_y11ddz__ddz=
     1   ddbcd1_y11ddz__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_y11ddz__dbcd_y*dbcd_y__ddz
     1 + ddbcd1_y11ddz__ddbcd_length_inv11ddz*ddbcd_length_inv11ddz__ddz


c dbcd1_z__dbx=-bcd_length_inv*cd_y + bcd_z*dbcd_length_inv__dbx
      ddbcd1_z11dbx__dbx=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__dbx
     1 + ddbcd1_z11dbx__dbcd_z*dbcd_z__dbx
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dbx
      ddbcd1_z11dbx__dby=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__dby
     1 + ddbcd1_z11dbx__dbcd_z*dbcd_z__dby
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dby
      ddbcd1_z11dbx__dbz=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__dbz
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dbz
      ddbcd1_z11dbx__dcx=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_z11dbx__dbcd_z*dbcd_z__dcx
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcx
      ddbcd1_z11dbx__dcy=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_z11dbx__dcd_y*dcd_y__dcy
     1 + ddbcd1_z11dbx__dbcd_z*dbcd_z__dcy
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcy
      ddbcd1_z11dbx__dcz=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__dcz
      ddbcd1_z11dbx__ddx=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_z11dbx__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddx
      ddbcd1_z11dbx__ddy=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_z11dbx__dcd_y*dcd_y__ddy
     1 + ddbcd1_z11dbx__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddy
      ddbcd1_z11dbx__ddz=
     1   ddbcd1_z11dbx__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_z11dbx__ddbcd_length_inv11dbx*ddbcd_length_inv11dbx__ddz
c dbcd1_z__dby=bcd_length_inv*cd_x + bcd_z*dbcd_length_inv__dby
      ddbcd1_z11dby__dby=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__dby
     1 + ddbcd1_z11dby__dbcd_z*dbcd_z__dby
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dby
      ddbcd1_z11dby__dbz=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__dbz
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dbz
      ddbcd1_z11dby__dcx=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_z11dby__dcd_x*dcd_x__dcx
     1 + ddbcd1_z11dby__dbcd_z*dbcd_z__dcx
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcx
      ddbcd1_z11dby__dcy=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_z11dby__dbcd_z*dbcd_z__dcy
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcy
      ddbcd1_z11dby__dcz=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__dcz
      ddbcd1_z11dby__ddx=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_z11dby__dcd_x*dcd_x__ddx
     1 + ddbcd1_z11dby__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddx
      ddbcd1_z11dby__ddy=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_z11dby__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddy
      ddbcd1_z11dby__ddz=
     1   ddbcd1_z11dby__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_z11dby__ddbcd_length_inv11dby*ddbcd_length_inv11dby__ddz
c dbcd1_z__dbz= bcd_z*dbcd_length_inv__dbz
      ddbcd1_z11dbz__dbz=
     1  ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dbz
      ddbcd1_z11dbz__dcx=
     1   ddbcd1_z11dbz__dbcd_z*dbcd_z__dcx
     1 + ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcx
      ddbcd1_z11dbz__dcy=
     1   ddbcd1_z11dbz__dbcd_z*dbcd_z__dcy
     1 + ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcy
      ddbcd1_z11dbz__dcz=
     1  ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__dcz
      ddbcd1_z11dbz__ddx=
     1   ddbcd1_z11dbz__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddx
      ddbcd1_z11dbz__ddy=
     1   ddbcd1_z11dbz__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddy
      ddbcd1_z11dbz__ddz=
     1  ddbcd1_z11dbz__ddbcd_length_inv11dbz*ddbcd_length_inv11dbz__ddz
c dbcd1_z__dcx=bcd_length_inv*dbcd_z__dcx + bcd_z*dbcd_length_inv__dcx
      ddbcd1_z11dcx__dcx=
     1   ddbcd1_z11dcx__dbcd_length_inv*dbcd_length_inv__dcx
     1 + ddbcd1_z11dcx__dbcd_z*dbcd_z__dcx
     1 + ddbcd1_z11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcx
      ddbcd1_z11dcx__dcy=
     1   ddbcd1_z11dcx__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_z11dcx__ddbcd_z11dcx*ddbcd_z11dcx__dcy
     1 + ddbcd1_z11dcx__dbcd_z*dbcd_z__dcy
     1 + ddbcd1_z11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcy
      ddbcd1_z11dcx__dcz=
     1   ddbcd1_z11dcx__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_z11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__dcz
      ddbcd1_z11dcx__ddx=
     1   ddbcd1_z11dcx__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_z11dcx__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddx
      ddbcd1_z11dcx__ddy=
     1   ddbcd1_z11dcx__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_z11dcx__ddbcd_z11dcx*ddbcd_z11dcx__ddy
     1 + ddbcd1_z11dcx__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddy
      ddbcd1_z11dcx__ddz=
     1   ddbcd1_z11dcx__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_z11dcx__ddbcd_length_inv11dcx*ddbcd_length_inv11dcx__ddz
c dbcd1_z__dcy=bcd_length_inv*dbcd_z__dcy + bcd_z*dbcd_length_inv__dcy
      ddbcd1_z11dcy__dcy=
     1   ddbcd1_z11dcy__dbcd_length_inv*dbcd_length_inv__dcy
     1 + ddbcd1_z11dcy__dbcd_z*dbcd_z__dcy
     1 + ddbcd1_z11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__dcy
      ddbcd1_z11dcy__dcz=
     1   ddbcd1_z11dcy__dbcd_length_inv*dbcd_length_inv__dcz
     1 + ddbcd1_z11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__dcz
      ddbcd1_z11dcy__ddx=
     1   ddbcd1_z11dcy__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_z11dcy__ddbcd_z11dcy*ddbcd_z11dcy__ddx
     1 + ddbcd1_z11dcy__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddx
      ddbcd1_z11dcy__ddy=
     1   ddbcd1_z11dcy__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_z11dcy__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddy
      ddbcd1_z11dcy__ddz=
     1   ddbcd1_z11dcy__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_z11dcy__ddbcd_length_inv11dcy*ddbcd_length_inv11dcy__ddz
c dbcd1_z__dcz= bcd_z*dbcd_length_inv__dcz
      ddbcd1_z11dcz__dcz=
     1  ddbcd1_z11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__dcz
      ddbcd1_z11dcz__ddx=
     1   ddbcd1_z11dcz__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddx
      ddbcd1_z11dcz__ddy=
     1   ddbcd1_z11dcz__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddy
      ddbcd1_z11dcz__ddz=
     1  ddbcd1_z11dcz__ddbcd_length_inv11dcz*ddbcd_length_inv11dcz__ddz
c dbcd1_z__ddx=-bcd_length_inv*bc_y + bcd_z*dbcd_length_inv__ddx
      ddbcd1_z11ddx__ddx=
     1   ddbcd1_z11ddx__dbcd_length_inv*dbcd_length_inv__ddx
     1 + ddbcd1_z11ddx__dbcd_z*dbcd_z__ddx
     1 + ddbcd1_z11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddx
      ddbcd1_z11ddx__ddy=
     1   ddbcd1_z11ddx__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_z11ddx__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddy
      ddbcd1_z11ddx__ddz=
     1   ddbcd1_z11ddx__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_z11ddx__ddbcd_length_inv11ddx*ddbcd_length_inv11ddx__ddz
c dbcd1_z__ddy=bcd_length_inv*bc_x + bcd_z*dbcd_length_inv__ddy
      ddbcd1_z11ddy__ddy=
     1   ddbcd1_z11ddy__dbcd_length_inv*dbcd_length_inv__ddy
     1 + ddbcd1_z11ddy__dbcd_z*dbcd_z__ddy
     1 + ddbcd1_z11ddy__ddbcd_length_inv11ddy*ddbcd_length_inv11ddy__ddy
      ddbcd1_z11ddy__ddz=
     1   ddbcd1_z11ddy__dbcd_length_inv*dbcd_length_inv__ddz
     1 + ddbcd1_z11ddy__ddbcd_length_inv11ddy*ddbcd_length_inv11ddy__ddz
c dbcd1_z__ddz= bcd_z*dbcd_length_inv__ddz
      ddbcd1_z11ddz__ddz=
     1  ddbcd1_z11ddz__ddbcd_length_inv11ddz*ddbcd_length_inv11ddz__ddz


c dy__dax= bcd1_x*daux_x__dax + bcd1_y*daux_y__dax + bcd1_z*daux_z__dax
      ddy11dax__dbcd1_x=daux_x__dax
      ddy11dax__dbcd1_y=daux_y__dax
      ddy11dax__dbcd1_z=daux_z__dax
      ddy11dax__ddaux_x11dax=bcd1_x
      ddy11dax__ddaux_y11dax=bcd1_y
      ddy11dax__ddaux_z11dax=bcd1_z
c dy__day= bcd1_x*daux_x__day + bcd1_y*daux_y__day + bcd1_z*daux_z__day
      ddy11day__dbcd1_x=daux_x__day
      ddy11day__dbcd1_y=daux_y__day
      ddy11day__dbcd1_z=daux_z__day
      ddy11day__ddaux_x11day=bcd1_x
      ddy11day__ddaux_y11day=bcd1_y
      ddy11day__ddaux_z11day=bcd1_z
c dy__daz= bcd1_x*daux_x__daz + bcd1_y*daux_y__daz + bcd1_z*daux_z__daz
      ddy11daz__dbcd1_x=daux_x__daz
      ddy11daz__dbcd1_y=daux_y__daz 
      ddy11daz__dbcd1_z=daux_z__daz
      ddy11daz__ddaux_x11daz=bcd1_x
      ddy11daz__ddaux_y11daz=bcd1_y
      ddy11daz__ddaux_z11daz=bcd1_z
c dy__dbx= bcd1_x*daux_x__dbx + aux_x*dbcd1_x__dbx + bcd1_y*daux_y__dbx + aux_y*dbcd1_y__dbx + bcd1_z*daux_z__dbx + aux_z*dbcd1_z__dbx
      ddy11dbx__daux_x=dbcd1_x__dbx
      ddy11dbx__daux_y=dbcd1_y__dbx
      ddy11dbx__daux_z=dbcd1_z__dbx
      ddy11dbx__dbcd1_x=daux_x__dbx
      ddy11dbx__dbcd1_y=daux_y__dbx
      ddy11dbx__dbcd1_z=daux_z__dbx
      ddy11dbx__ddaux_x11dbx=bcd1_x
      ddy11dbx__ddaux_y11dbx=bcd1_y
      ddy11dbx__ddaux_z11dbx=bcd1_z
      ddy11dbx__ddbcd1_x11dbx=aux_x
      ddy11dbx__ddbcd1_y11dbx=aux_y
      ddy11dbx__ddbcd1_z11dbx=aux_z
c dy__dby= bcd1_x*daux_x__dby + aux_x*dbcd1_x__dby + bcd1_y*daux_y__dby + aux_y*dbcd1_y__dby + bcd1_z*daux_z__dby + aux_z*dbcd1_z__dby
      ddy11dby__daux_x=dbcd1_x__dby
      ddy11dby__daux_y=dbcd1_y__dby
      ddy11dby__daux_z=dbcd1_z__dby
      ddy11dby__dbcd1_x=daux_x__dby
      ddy11dby__dbcd1_y=daux_y__dby
      ddy11dby__dbcd1_z=daux_z__dby
      ddy11dby__ddaux_x11dby=bcd1_x
      ddy11dby__ddaux_y11dby=bcd1_y
      ddy11dby__ddaux_z11dby=bcd1_z
      ddy11dby__ddbcd1_x11dby=aux_x
      ddy11dby__ddbcd1_y11dby=aux_y
      ddy11dby__ddbcd1_z11dby=aux_z
c dy__dbz= bcd1_x*daux_x__dbz + aux_x*dbcd1_x__dbz + bcd1_y*daux_y__dbz + aux_y*dbcd1_y__dbz + bcd1_z*daux_z__dbz + aux_z*dbcd1_z__dbz
      ddy11dbz__daux_x=dbcd1_x__dbz
      ddy11dbz__daux_y=dbcd1_y__dbz
      ddy11dbz__daux_z=dbcd1_z__dbz
      ddy11dbz__dbcd1_x=daux_x__dbz
      ddy11dbz__dbcd1_y=daux_y__dbz
      ddy11dbz__dbcd1_z=daux_z__dbz
      ddy11dbz__ddaux_x11dbz=bcd1_x
      ddy11dbz__ddaux_y11dbz=bcd1_y
      ddy11dbz__ddaux_z11dbz=bcd1_z
      ddy11dbz__ddbcd1_x11dbz=aux_x
      ddy11dbz__ddbcd1_y11dbz=aux_y
      ddy11dbz__ddbcd1_z11dbz=aux_z
c dy__dcx= bcd1_x*daux_x__dcx + aux_x*dbcd1_x__dcx + bcd1_y*daux_y__dcx + aux_y*dbcd1_y__dcx + bcd1_z*daux_z__dcx + aux_z*dbcd1_z__dcx
      ddy11dcx__daux_x=dbcd1_x__dcx
      ddy11dcx__daux_y=dbcd1_y__dcx
      ddy11dcx__daux_z=dbcd1_z__dcx
      ddy11dcx__dbcd1_x=daux_x__dcx
      ddy11dcx__dbcd1_y=daux_y__dcx
      ddy11dcx__dbcd1_z=daux_z__dcx
      ddy11dcx__ddaux_x11dcx=bcd1_x
      ddy11dcx__ddaux_y11dcx=bcd1_y
      ddy11dcx__ddaux_z11dcx=bcd1_z
      ddy11dcx__ddbcd1_x11dcx=aux_x
      ddy11dcx__ddbcd1_y11dcx=aux_y
      ddy11dcx__ddbcd1_z11dcx=aux_z
c dy__dcy= bcd1_x*daux_x__dcy + aux_x*dbcd1_x__dcy + bcd1_y*daux_y__dcy + aux_y*dbcd1_y__dcy + bcd1_z*daux_z__dcy + aux_z*dbcd1_z__dcy
      ddy11dcy__daux_x=dbcd1_x__dcy
      ddy11dcy__daux_y=dbcd1_y__dcy
      ddy11dcy__daux_z=dbcd1_z__dcy
      ddy11dcy__dbcd1_x=daux_x__dcy
      ddy11dcy__dbcd1_y=daux_y__dcy
      ddy11dcy__dbcd1_z=daux_z__dcy
      ddy11dcy__ddaux_x11dcy=bcd1_x
      ddy11dcy__ddaux_y11dcy=bcd1_y
      ddy11dcy__ddaux_z11dcy=bcd1_z
      ddy11dcy__ddbcd1_x11dcy=aux_x
      ddy11dcy__ddbcd1_y11dcy=aux_y
      ddy11dcy__ddbcd1_z11dcy=aux_z
c dy__dcz= bcd1_x*daux_x__dcz + aux_x*dbcd1_x__dcz + bcd1_y*daux_y__dcz + aux_y*dbcd1_y__dcz + bcd1_z*daux_z__dcz + aux_z*dbcd1_z__dcz
      ddy11dcz__daux_x=dbcd1_x__dcz
      ddy11dcz__daux_y=dbcd1_y__dcz
      ddy11dcz__daux_z=dbcd1_z__dcz
      ddy11dcz__dbcd1_x=daux_x__dcz
      ddy11dcz__dbcd1_y=daux_y__dcz
      ddy11dcz__dbcd1_z=daux_z__dcz
      ddy11dcz__ddaux_x11dcz=bcd1_x
      ddy11dcz__ddaux_y11dcz=bcd1_y
      ddy11dcz__ddaux_z11dcz=bcd1_z
      ddy11dcz__ddbcd1_x11dcz=aux_x
      ddy11dcz__ddbcd1_y11dcz=aux_y
      ddy11dcz__ddbcd1_z11dcz=aux_z
c dy__ddx= aux_x*dbcd1_x__ddx + aux_y*dbcd1_y__ddx + aux_z*dbcd1_z__ddx
      ddy11ddx__daux_x=dbcd1_x__ddx
      ddy11ddx__daux_y=dbcd1_y__ddx
      ddy11ddx__daux_z=dbcd1_z__ddx
      ddy11ddx__ddbcd1_x11ddx=aux_x
      ddy11ddx__ddbcd1_y11ddx=aux_y
      ddy11ddx__ddbcd1_z11ddx=aux_z
c dy__ddy= aux_x*dbcd1_x__ddy + aux_y*dbcd1_y__ddy + aux_z*dbcd1_z__ddy
      ddy11ddy__daux_x=dbcd1_x__ddy
      ddy11ddy__daux_y=dbcd1_y__ddy
      ddy11ddy__daux_z=dbcd1_z__ddy
      ddy11ddy__ddbcd1_x11ddy=aux_x
      ddy11ddy__ddbcd1_y11ddy=aux_y
      ddy11ddy__ddbcd1_z11ddy=aux_z
c dy__ddz= aux_x*dbcd1_x__ddz + aux_y*dbcd1_y__ddz + aux_z*dbcd1_z__ddz
      ddy11ddz__daux_x=dbcd1_x__ddz
      ddy11ddz__daux_y=dbcd1_y__ddz
      ddy11ddz__daux_z=dbcd1_z__ddz
      ddy11ddz__ddbcd1_x11ddz=aux_x
      ddy11ddz__ddbcd1_y11ddz=aux_y
      ddy11ddz__ddbcd1_z11ddz=aux_z


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c block 2

c df__dx=-y/(x**2 + y**2)
      ddf11dx__dax=ddf11dx__dx*dx__dax + ddf11dx__dy*dy__dax
      ddf11dx__day=ddf11dx__dx*dx__day + ddf11dx__dy*dy__day
      ddf11dx__daz=ddf11dx__dx*dx__daz + ddf11dx__dy*dy__daz
      ddf11dx__dbx=ddf11dx__dx*dx__dbx + ddf11dx__dy*dy__dbx
      ddf11dx__dby=ddf11dx__dx*dx__dby + ddf11dx__dy*dy__dby
      ddf11dx__dbz=ddf11dx__dx*dx__dbz + ddf11dx__dy*dy__dbz
      ddf11dx__dcx=ddf11dx__dx*dx__dcx + ddf11dx__dy*dy__dcx
      ddf11dx__dcy=ddf11dx__dx*dx__dcy + ddf11dx__dy*dy__dcy
      ddf11dx__dcz=ddf11dx__dx*dx__dcz + ddf11dx__dy*dy__dcz
      ddf11dx__ddx=ddf11dx__dx*dx__ddx + ddf11dx__dy*dy__ddx
      ddf11dx__ddy=ddf11dx__dx*dx__ddy + ddf11dx__dy*dy__ddy
      ddf11dx__ddz=ddf11dx__dx*dx__ddz + ddf11dx__dy*dy__ddz


c df__dy=x/(x**2 + y**2)
      ddf11dy__dax=ddf11dy__dx*dx__dax + ddf11dy__dy*dy__dax
      ddf11dy__day=ddf11dy__dx*dx__day + ddf11dy__dy*dy__day
      ddf11dy__daz=ddf11dy__dx*dx__daz + ddf11dy__dy*dy__daz
      ddf11dy__dbx=ddf11dy__dx*dx__dbx + ddf11dy__dy*dy__dbx
      ddf11dy__dby=ddf11dy__dx*dx__dby + ddf11dy__dy*dy__dby
      ddf11dy__dbz=ddf11dy__dx*dx__dbz + ddf11dy__dy*dy__dbz
      ddf11dy__dcx=ddf11dy__dx*dx__dcx + ddf11dy__dy*dy__dcx
      ddf11dy__dcy=ddf11dy__dx*dx__dcy + ddf11dy__dy*dy__dcy
      ddf11dy__dcz=ddf11dy__dx*dx__dcz + ddf11dy__dy*dy__dcz
      ddf11dy__ddx=ddf11dy__dx*dx__ddx + ddf11dy__dy*dy__ddx
      ddf11dy__ddy=ddf11dy__dx*dx__ddy + ddf11dy__dy*dy__ddy
      ddf11dy__ddz=ddf11dy__dx*dx__ddz + ddf11dy__dy*dy__ddz


c dx__dax= bcd1_x*dabc1_x__dax + bcd1_y*dabc1_y__dax + bcd1_z*dabc1_z__dax  
      ddx11dax__dax=
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dax
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dax
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dax
      ddx11dax__day=
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__day
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__day
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__day
      ddx11dax__daz=
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__daz
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__daz
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__daz
      ddx11dax__dbx=
     1              ddx11dax__dbcd1_x*dbcd1_x__dbx
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dbx
     1            + ddx11dax__dbcd1_y*dbcd1_y__dbx
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dbx
     1            + ddx11dax__dbcd1_z*dbcd1_z__dbx
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dbx
      ddx11dax__dby=
     1              ddx11dax__dbcd1_x*dbcd1_x__dby
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dby
     1            + ddx11dax__dbcd1_y*dbcd1_y__dby
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dby
     1            + ddx11dax__dbcd1_z*dbcd1_z__dby
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dby
      ddx11dax__dbz=
     1              ddx11dax__dbcd1_x*dbcd1_x__dbz
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dbz
     1            + ddx11dax__dbcd1_y*dbcd1_y__dbz
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dbz
     1            + ddx11dax__dbcd1_z*dbcd1_z__dbz
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dbz
      ddx11dax__dcx=
     1              ddx11dax__dbcd1_x*dbcd1_x__dcx
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dcx
     1            + ddx11dax__dbcd1_y*dbcd1_y__dcx
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dcx
     1            + ddx11dax__dbcd1_z*dbcd1_z__dcx
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dcx
      ddx11dax__dcy=
     1              ddx11dax__dbcd1_x*dbcd1_x__dcy
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dcy
     1            + ddx11dax__dbcd1_y*dbcd1_y__dcy
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dcy
     1            + ddx11dax__dbcd1_z*dbcd1_z__dcy
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dcy
      ddx11dax__dcz=
     1              ddx11dax__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dax__ddabc1_x11dax*ddabc1_x11dax__dcz
     1            + ddx11dax__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dax__ddabc1_y11dax*ddabc1_y11dax__dcz
     1            + ddx11dax__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dax__ddabc1_z11dax*ddabc1_z11dax__dcz
      ddx11dax__ddx=
     1              ddx11dax__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dax__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dax__dbcd1_z*dbcd1_z__ddx
      ddx11dax__ddy=
     1              ddx11dax__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dax__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dax__dbcd1_z*dbcd1_z__ddy
      ddx11dax__ddz=
     1              ddx11dax__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dax__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dax__dbcd1_z*dbcd1_z__ddz
c      dx__day= bcd1_x*dabc1_x__day + bcd1_y*dabc1_y__day + bcd1_z*dabc1_z__day 
      ddx11day__day=
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__day
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__day
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__day
      ddx11day__daz=
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__daz
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__daz
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__daz
      ddx11day__dbx=
     1              ddx11day__dbcd1_x*dbcd1_x__dbx
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__dbx
     1            + ddx11day__dbcd1_y*dbcd1_y__dbx
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__dbx
     1            + ddx11day__dbcd1_z*dbcd1_z__dbx
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__dbx
      ddx11day__dby=
     1              ddx11day__dbcd1_x*dbcd1_x__dby
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__dby
     1            + ddx11day__dbcd1_y*dbcd1_y__dby
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__dby
     1            + ddx11day__dbcd1_z*dbcd1_z__dby
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__dby
      ddx11day__dbz=
     1              ddx11day__dbcd1_x*dbcd1_x__dbz
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__dbz
     1            + ddx11day__dbcd1_y*dbcd1_y__dbz
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__dbz
     1            + ddx11day__dbcd1_z*dbcd1_z__dbz
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__dbz
      ddx11day__dcx=
     1              ddx11day__dbcd1_x*dbcd1_x__dcx
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__dcx
     1            + ddx11day__dbcd1_y*dbcd1_y__dcx
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__dcx
     1            + ddx11day__dbcd1_z*dbcd1_z__dcx
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__dcx
      ddx11day__dcy=
     1              ddx11day__dbcd1_x*dbcd1_x__dcy
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__dcy
     1            + ddx11day__dbcd1_y*dbcd1_y__dcy
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__dcy
     1            + ddx11day__dbcd1_z*dbcd1_z__dcy
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__dcy
      ddx11day__dcz=
     1              ddx11day__dbcd1_x*dbcd1_x__dcz
     1            + ddx11day__ddabc1_x11day*ddabc1_x11day__dcz
     1            + ddx11day__dbcd1_y*dbcd1_y__dcz
     1            + ddx11day__ddabc1_y11day*ddabc1_y11day__dcz
     1            + ddx11day__dbcd1_z*dbcd1_z__dcz
     1            + ddx11day__ddabc1_z11day*ddabc1_z11day__dcz
      ddx11day__ddx=
     1              ddx11day__dbcd1_x*dbcd1_x__ddx
     1            + ddx11day__dbcd1_y*dbcd1_y__ddx
     1            + ddx11day__dbcd1_z*dbcd1_z__ddx
      ddx11day__ddy=
     1              ddx11day__dbcd1_x*dbcd1_x__ddy
     1            + ddx11day__dbcd1_y*dbcd1_y__ddy
     1            + ddx11day__dbcd1_z*dbcd1_z__ddy
      ddx11day__ddz=
     1              ddx11day__dbcd1_x*dbcd1_x__ddz
     1            + ddx11day__dbcd1_y*dbcd1_y__ddz
     1            + ddx11day__dbcd1_z*dbcd1_z__ddz
c      dx__daz= bcd1_x*dabc1_x__daz + bcd1_y*dabc1_y__daz + bcd1_z*dabc1_z__daz 
      ddx11daz__daz=
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__daz
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__daz
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__daz
      ddx11daz__dbx=
     1              ddx11daz__dbcd1_x*dbcd1_x__dbx
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__dbx
     1            + ddx11daz__dbcd1_y*dbcd1_y__dbx
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__dbx
     1            + ddx11daz__dbcd1_z*dbcd1_z__dbx
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__dbx
      ddx11daz__dby=
     1              ddx11daz__dbcd1_x*dbcd1_x__dby
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__dby
     1            + ddx11daz__dbcd1_y*dbcd1_y__dby
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__dby
     1            + ddx11daz__dbcd1_z*dbcd1_z__dby
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__dby
      ddx11daz__dbz=
     1              ddx11daz__dbcd1_x*dbcd1_x__dbz
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__dbz
     1            + ddx11daz__dbcd1_y*dbcd1_y__dbz
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__dbz
     1            + ddx11daz__dbcd1_z*dbcd1_z__dbz
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__dbz
      ddx11daz__dcx=
     1              ddx11daz__dbcd1_x*dbcd1_x__dcx
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__dcx
     1            + ddx11daz__dbcd1_y*dbcd1_y__dcx
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__dcx
     1            + ddx11daz__dbcd1_z*dbcd1_z__dcx
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__dcx
      ddx11daz__dcy=
     1              ddx11daz__dbcd1_x*dbcd1_x__dcy
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__dcy
     1            + ddx11daz__dbcd1_y*dbcd1_y__dcy
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__dcy
     1            + ddx11daz__dbcd1_z*dbcd1_z__dcy
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__dcy
      ddx11daz__dcz=
     1              ddx11daz__dbcd1_x*dbcd1_x__dcz
     1            + ddx11daz__ddabc1_x11daz*ddabc1_x11daz__dcz
     1            + ddx11daz__dbcd1_y*dbcd1_y__dcz
     1            + ddx11daz__ddabc1_y11daz*ddabc1_y11daz__dcz
     1            + ddx11daz__dbcd1_z*dbcd1_z__dcz
     1            + ddx11daz__ddabc1_z11daz*ddabc1_z11daz__dcz
      ddx11daz__ddx=
     1              ddx11daz__dbcd1_x*dbcd1_x__ddx
     1            + ddx11daz__dbcd1_y*dbcd1_y__ddx
     1            + ddx11daz__dbcd1_z*dbcd1_z__ddx
      ddx11daz__ddy=
     1              ddx11daz__dbcd1_x*dbcd1_x__ddy
     1            + ddx11daz__dbcd1_y*dbcd1_y__ddy
     1            + ddx11daz__dbcd1_z*dbcd1_z__ddy
      ddx11daz__ddz=
     1              ddx11daz__dbcd1_x*dbcd1_x__ddz
     1            + ddx11daz__dbcd1_y*dbcd1_y__ddz
     1            + ddx11daz__dbcd1_z*dbcd1_z__ddz
c      dx__dbx= bcd1_x*dabc1_x__dbx + abc1_x*dbcd1_x__dbx + bcd1_y*dabc1_y__dbx + abc1_y*dbcd1_y__dbx + bcd1_z*dabc1_z__dbx + abc1_z*dbcd1_z__dbx
      ddx11dbx__dbx=
     1              ddx11dbx__dbcd1_x*dbcd1_x__dbx
     1            + ddx11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dbx
     1            + ddx11dbx__dabc1_x*dabc1_x__dbx
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dbx
     1            + ddx11dbx__dbcd1_y*dbcd1_y__dbx
     1            + ddx11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dbx
     1            + ddx11dbx__dabc1_y*dabc1_y__dbx
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dbx
     1            + ddx11dbx__dbcd1_z*dbcd1_z__dbx
     1            + ddx11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dbx
     1            + ddx11dbx__dabc1_z*dabc1_z__dbx
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dbx
      ddx11dbx__dby=
     1              ddx11dbx__dbcd1_x*dbcd1_x__dby
     1            + ddx11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dby
     1            + ddx11dbx__dabc1_x*dabc1_x__dby
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dby
     1            + ddx11dbx__dbcd1_y*dbcd1_y__dby
     1            + ddx11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dby
     1            + ddx11dbx__dabc1_y*dabc1_y__dby
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dby
     1            + ddx11dbx__dbcd1_z*dbcd1_z__dby
     1            + ddx11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dby
     1            + ddx11dbx__dabc1_z*dabc1_z__dby
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dby
      ddx11dbx__dbz=
     1              ddx11dbx__dbcd1_x*dbcd1_x__dbz
     1            + ddx11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dbz
     1            + ddx11dbx__dabc1_x*dabc1_x__dbz
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dbz
     1            + ddx11dbx__dbcd1_y*dbcd1_y__dbz
     1            + ddx11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dbz
     1            + ddx11dbx__dabc1_y*dabc1_y__dbz
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dbz
     1            + ddx11dbx__dbcd1_z*dbcd1_z__dbz
     1            + ddx11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dbz
     1            + ddx11dbx__dabc1_z*dabc1_z__dbz
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dbz
      ddx11dbx__dcx=
     1              ddx11dbx__dbcd1_x*dbcd1_x__dcx
     1            + ddx11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcx
     1            + ddx11dbx__dabc1_x*dabc1_x__dcx
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dcx
     1            + ddx11dbx__dbcd1_y*dbcd1_y__dcx
     1            + ddx11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcx
     1            + ddx11dbx__dabc1_y*dabc1_y__dcx
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dcx
     1            + ddx11dbx__dbcd1_z*dbcd1_z__dcx
     1            + ddx11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcx
     1            + ddx11dbx__dabc1_z*dabc1_z__dcx
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dcx
      ddx11dbx__dcy=
     1              ddx11dbx__dbcd1_x*dbcd1_x__dcy
     1            + ddx11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcy
     1            + ddx11dbx__dabc1_x*dabc1_x__dcy
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dcy
     1            + ddx11dbx__dbcd1_y*dbcd1_y__dcy
     1            + ddx11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcy
     1            + ddx11dbx__dabc1_y*dabc1_y__dcy
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dcy
     1            + ddx11dbx__dbcd1_z*dbcd1_z__dcy
     1            + ddx11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcy
     1            + ddx11dbx__dabc1_z*dabc1_z__dcy
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dcy
      ddx11dbx__dcz=
     1              ddx11dbx__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dbx__ddabc1_x11dbx*ddabc1_x11dbx__dcz
     1            + ddx11dbx__dabc1_x*dabc1_x__dcz
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dcz
     1            + ddx11dbx__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dbx__ddabc1_y11dbx*ddabc1_y11dbx__dcz
     1            + ddx11dbx__dabc1_y*dabc1_y__dcz
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dcz
     1            + ddx11dbx__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dbx__ddabc1_z11dbx*ddabc1_z11dbx__dcz
     1            + ddx11dbx__dabc1_z*dabc1_z__dcz
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dcz
      ddx11dbx__ddx=
     1              ddx11dbx__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__ddx
     1            + ddx11dbx__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__ddx
     1            + ddx11dbx__dbcd1_z*dbcd1_z__ddx
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__ddx
      ddx11dbx__ddy=
     1              ddx11dbx__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__ddy
     1            + ddx11dbx__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__ddy
     1            + ddx11dbx__dbcd1_z*dbcd1_z__ddy
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__ddy
      ddx11dbx__ddz=
     1              ddx11dbx__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__ddz
     1            + ddx11dbx__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__ddz
     1            + ddx11dbx__dbcd1_z*dbcd1_z__ddz
     1            + ddx11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__ddz
c      dx__dby= bcd1_x*dabc1_x__dby + abc1_x*dbcd1_x__dby + bcd1_y*dabc1_y__dby + abc1_y*dbcd1_y__dby + bcd1_z*dabc1_z__dby + abc1_z*dbcd1_z__dby
      ddx11dby__dby=
     1              ddx11dby__dbcd1_x*dbcd1_x__dby
     1            + ddx11dby__ddabc1_x11dby*ddabc1_x11dby__dby
     1            + ddx11dby__dabc1_x*dabc1_x__dby
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__dby
     1            + ddx11dby__dbcd1_y*dbcd1_y__dby
     1            + ddx11dby__ddabc1_y11dby*ddabc1_y11dby__dby
     1            + ddx11dby__dabc1_y*dabc1_y__dby
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__dby
     1            + ddx11dby__dbcd1_z*dbcd1_z__dby
     1            + ddx11dby__ddabc1_z11dby*ddabc1_z11dby__dby
     1            + ddx11dby__dabc1_z*dabc1_z__dby
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__dby
      ddx11dby__dbz=
     1              ddx11dby__dbcd1_x*dbcd1_x__dbz
     1            + ddx11dby__ddabc1_x11dby*ddabc1_x11dby__dbz
     1            + ddx11dby__dabc1_x*dabc1_x__dbz
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__dbz
     1            + ddx11dby__dbcd1_y*dbcd1_y__dbz
     1            + ddx11dby__ddabc1_y11dby*ddabc1_y11dby__dbz
     1            + ddx11dby__dabc1_y*dabc1_y__dbz
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__dbz
     1            + ddx11dby__dbcd1_z*dbcd1_z__dbz
     1            + ddx11dby__ddabc1_z11dby*ddabc1_z11dby__dbz
     1            + ddx11dby__dabc1_z*dabc1_z__dbz
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__dbz
      ddx11dby__dcx=
     1              ddx11dby__dbcd1_x*dbcd1_x__dcx
     1            + ddx11dby__ddabc1_x11dby*ddabc1_x11dby__dcx
     1            + ddx11dby__dabc1_x*dabc1_x__dcx
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__dcx
     1            + ddx11dby__dbcd1_y*dbcd1_y__dcx
     1            + ddx11dby__ddabc1_y11dby*ddabc1_y11dby__dcx
     1            + ddx11dby__dabc1_y*dabc1_y__dcx
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__dcx
     1            + ddx11dby__dbcd1_z*dbcd1_z__dcx
     1            + ddx11dby__ddabc1_z11dby*ddabc1_z11dby__dcx
     1            + ddx11dby__dabc1_z*dabc1_z__dcx
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__dcx
      ddx11dby__dcy=
     1              ddx11dby__dbcd1_x*dbcd1_x__dcy
     1            + ddx11dby__ddabc1_x11dby*ddabc1_x11dby__dcy
     1            + ddx11dby__dabc1_x*dabc1_x__dcy
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__dcy
     1            + ddx11dby__dbcd1_y*dbcd1_y__dcy
     1            + ddx11dby__ddabc1_y11dby*ddabc1_y11dby__dcy
     1            + ddx11dby__dabc1_y*dabc1_y__dcy
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__dcy
     1            + ddx11dby__dbcd1_z*dbcd1_z__dcy
     1            + ddx11dby__ddabc1_z11dby*ddabc1_z11dby__dcy
     1            + ddx11dby__dabc1_z*dabc1_z__dcy
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__dcy
      ddx11dby__dcz=
     1              ddx11dby__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dby__ddabc1_x11dby*ddabc1_x11dby__dcz
     1            + ddx11dby__dabc1_x*dabc1_x__dcz
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__dcz
     1            + ddx11dby__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dby__ddabc1_y11dby*ddabc1_y11dby__dcz
     1            + ddx11dby__dabc1_y*dabc1_y__dcz
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__dcz
     1            + ddx11dby__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dby__ddabc1_z11dby*ddabc1_z11dby__dcz
     1            + ddx11dby__dabc1_z*dabc1_z__dcz
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__dcz
      ddx11dby__ddx=
     1              ddx11dby__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__ddx
     1            + ddx11dby__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__ddx
     1            + ddx11dby__dbcd1_z*dbcd1_z__ddx
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__ddx
      ddx11dby__ddy=
     1              ddx11dby__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__ddy
     1            + ddx11dby__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__ddy
     1            + ddx11dby__dbcd1_z*dbcd1_z__ddy
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__ddy
      ddx11dby__ddz=
     1              ddx11dby__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dby__ddbcd1_x11dby*ddbcd1_x11dby__ddz
     1            + ddx11dby__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dby__ddbcd1_y11dby*ddbcd1_y11dby__ddz
     1            + ddx11dby__dbcd1_z*dbcd1_z__ddz
     1            + ddx11dby__ddbcd1_z11dby*ddbcd1_z11dby__ddz
c      dx__dbz= bcd1_x*dabc1_x__dbz + abc1_x*dbcd1_x__dbz + bcd1_y*dabc1_y__dbz + abc1_y*dbcd1_y__dbz + bcd1_z*dabc1_z__dbz + abc1_z*dbcd1_z__dbz
      ddx11dbz__dbz=
     1              ddx11dbz__dbcd1_x*dbcd1_x__dbz
     1            + ddx11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dbz
     1            + ddx11dbz__dabc1_x*dabc1_x__dbz
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dbz
     1            + ddx11dbz__dbcd1_y*dbcd1_y__dbz
     1            + ddx11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dbz
     1            + ddx11dbz__dabc1_y*dabc1_y__dbz
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dbz
     1            + ddx11dbz__dbcd1_z*dbcd1_z__dbz
     1            + ddx11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dbz
     1            + ddx11dbz__dabc1_z*dabc1_z__dbz
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dbz
      ddx11dbz__dcx=
     1              ddx11dbz__dbcd1_x*dbcd1_x__dcx
     1            + ddx11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcx
     1            + ddx11dbz__dabc1_x*dabc1_x__dcx
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dcx
     1            + ddx11dbz__dbcd1_y*dbcd1_y__dcx
     1            + ddx11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcx
     1            + ddx11dbz__dabc1_y*dabc1_y__dcx
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dcx
     1            + ddx11dbz__dbcd1_z*dbcd1_z__dcx
     1            + ddx11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcx
     1            + ddx11dbz__dabc1_z*dabc1_z__dcx
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dcx
      ddx11dbz__dcy=
     1              ddx11dbz__dbcd1_x*dbcd1_x__dcy
     1            + ddx11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcy
     1            + ddx11dbz__dabc1_x*dabc1_x__dcy
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dcy
     1            + ddx11dbz__dbcd1_y*dbcd1_y__dcy
     1            + ddx11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcy
     1            + ddx11dbz__dabc1_y*dabc1_y__dcy
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dcy
     1            + ddx11dbz__dbcd1_z*dbcd1_z__dcy
     1            + ddx11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcy
     1            + ddx11dbz__dabc1_z*dabc1_z__dcy
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dcy
      ddx11dbz__dcz=
     1              ddx11dbz__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dbz__ddabc1_x11dbz*ddabc1_x11dbz__dcz
     1            + ddx11dbz__dabc1_x*dabc1_x__dcz
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dcz
     1            + ddx11dbz__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dbz__ddabc1_y11dbz*ddabc1_y11dbz__dcz
     1            + ddx11dbz__dabc1_y*dabc1_y__dcz
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dcz
     1            + ddx11dbz__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dbz__ddabc1_z11dbz*ddabc1_z11dbz__dcz
     1            + ddx11dbz__dabc1_z*dabc1_z__dcz
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dcz
      ddx11dbz__ddx=
     1              ddx11dbz__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__ddx
     1            + ddx11dbz__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__ddx
     1            + ddx11dbz__dbcd1_z*dbcd1_z__ddx
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__ddx
      ddx11dbz__ddy=
     1              ddx11dbz__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__ddy
     1            + ddx11dbz__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__ddy
     1            + ddx11dbz__dbcd1_z*dbcd1_z__ddy
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__ddy
      ddx11dbz__ddz=
     1              ddx11dbz__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__ddz
     1            + ddx11dbz__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__ddz
     1            + ddx11dbz__dbcd1_z*dbcd1_z__ddz
     1            + ddx11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__ddz
c      dx__dcx= bcd1_x*dabc1_x__dcx + abc1_x*dbcd1_x__dcx + bcd1_y*dabc1_y__dcx + abc1_y*dbcd1_y__dcx + bcd1_z*dabc1_z__dcx + abc1_z*dbcd1_z__dcx
      ddx11dcx__dcx=
     1              ddx11dcx__dbcd1_x*dbcd1_x__dcx
     1            + ddx11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcx
     1            + ddx11dcx__dabc1_x*dabc1_x__dcx
     1            + ddx11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__dcx
     1            + ddx11dcx__dbcd1_y*dbcd1_y__dcx
     1            + ddx11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcx
     1            + ddx11dcx__dabc1_y*dabc1_y__dcx
     1            + ddx11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__dcx
     1            + ddx11dcx__dbcd1_z*dbcd1_z__dcx
     1            + ddx11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcx
     1            + ddx11dcx__dabc1_z*dabc1_z__dcx
     1            + ddx11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__dcx
      ddx11dcx__dcy=
     1              ddx11dcx__dbcd1_x*dbcd1_x__dcy
     1            + ddx11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcy
     1            + ddx11dcx__dabc1_x*dabc1_x__dcy
     1            + ddx11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__dcy
     1            + ddx11dcx__dbcd1_y*dbcd1_y__dcy
     1            + ddx11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcy
     1            + ddx11dcx__dabc1_y*dabc1_y__dcy
     1            + ddx11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__dcy
     1            + ddx11dcx__dbcd1_z*dbcd1_z__dcy
     1            + ddx11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcy
     1            + ddx11dcx__dabc1_z*dabc1_z__dcy
     1            + ddx11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__dcy
      ddx11dcx__dcz=
     1              ddx11dcx__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dcx__ddabc1_x11dcx*ddabc1_x11dcx__dcz
     1            + ddx11dcx__dabc1_x*dabc1_x__dcz
     1            + ddx11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__dcz
     1            + ddx11dcx__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dcx__ddabc1_y11dcx*ddabc1_y11dcx__dcz
     1            + ddx11dcx__dabc1_y*dabc1_y__dcz
     1            + ddx11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__dcz
     1            + ddx11dcx__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dcx__ddabc1_z11dcx*ddabc1_z11dcx__dcz
     1            + ddx11dcx__dabc1_z*dabc1_z__dcz
     1            + ddx11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__dcz
      ddx11dcx__ddx=
     1              ddx11dcx__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__ddx
     1            + ddx11dcx__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__ddx
     1            + ddx11dcx__dbcd1_z*dbcd1_z__ddx
     1            + ddx11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__ddx
      ddx11dcx__ddy=
     1              ddx11dcx__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__ddy
     1            + ddx11dcx__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__ddy
     1            + ddx11dcx__dbcd1_z*dbcd1_z__ddy
     1            + ddx11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__ddy
      ddx11dcx__ddz=
     1              ddx11dcx__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__ddz
     1            + ddx11dcx__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__ddz
     1            + ddx11dcx__dbcd1_z*dbcd1_z__ddz
     1            + ddx11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__ddz
c      dx__dcy= bcd1_x*dabc1_x__dcy + abc1_x*dbcd1_x__dcy + bcd1_y*dabc1_y__dcy + abc1_y*dbcd1_y__dcy + bcd1_z*dabc1_z__dcy + abc1_z*dbcd1_z__dcy
      ddx11dcy__dcy=
     1              ddx11dcy__dbcd1_x*dbcd1_x__dcy
     1            + ddx11dcy__ddabc1_x11dcy*ddabc1_x11dcy__dcy
     1            + ddx11dcy__dabc1_x*dabc1_x__dcy
     1            + ddx11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__dcy
     1            + ddx11dcy__dbcd1_y*dbcd1_y__dcy
     1            + ddx11dcy__ddabc1_y11dcy*ddabc1_y11dcy__dcy
     1            + ddx11dcy__dabc1_y*dabc1_y__dcy
     1            + ddx11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__dcy
     1            + ddx11dcy__dbcd1_z*dbcd1_z__dcy
     1            + ddx11dcy__ddabc1_z11dcy*ddabc1_z11dcy__dcy
     1            + ddx11dcy__dabc1_z*dabc1_z__dcy
     1            + ddx11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__dcy
      ddx11dcy__dcz=
     1              ddx11dcy__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dcy__ddabc1_x11dcy*ddabc1_x11dcy__dcz
     1            + ddx11dcy__dabc1_x*dabc1_x__dcz
     1            + ddx11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__dcz
     1            + ddx11dcy__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dcy__ddabc1_y11dcy*ddabc1_y11dcy__dcz
     1            + ddx11dcy__dabc1_y*dabc1_y__dcz
     1            + ddx11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__dcz
     1            + ddx11dcy__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dcy__ddabc1_z11dcy*ddabc1_z11dcy__dcz
     1            + ddx11dcy__dabc1_z*dabc1_z__dcz
     1            + ddx11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__dcz
      ddx11dcy__ddx=
     1              ddx11dcy__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__ddx
     1            + ddx11dcy__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__ddx
     1            + ddx11dcy__dbcd1_z*dbcd1_z__ddx
     1            + ddx11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__ddx
      ddx11dcy__ddy=
     1              ddx11dcy__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__ddy
     1            + ddx11dcy__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__ddy
     1            + ddx11dcy__dbcd1_z*dbcd1_z__ddy
     1            + ddx11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__ddy
      ddx11dcy__ddz=
     1              ddx11dcy__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__ddz
     1            + ddx11dcy__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__ddz
     1            + ddx11dcy__dbcd1_z*dbcd1_z__ddz
     1            + ddx11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__ddz
c      dx__dcz= bcd1_x*dabc1_x__dcz + abc1_x*dbcd1_x__dcz + bcd1_y*dabc1_y__dcz + abc1_y*dbcd1_y__dcz + bcd1_z*dabc1_z__dcz + abc1_z*dbcd1_z__dcz
      ddx11dcz__dcz=
     1              ddx11dcz__dbcd1_x*dbcd1_x__dcz
     1            + ddx11dcz__ddabc1_x11dcz*ddabc1_x11dcz__dcz
     1            + ddx11dcz__dabc1_x*dabc1_x__dcz
     1            + ddx11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__dcz
     1            + ddx11dcz__dbcd1_y*dbcd1_y__dcz
     1            + ddx11dcz__ddabc1_y11dcz*ddabc1_y11dcz__dcz
     1            + ddx11dcz__dabc1_y*dabc1_y__dcz
     1            + ddx11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__dcz
     1            + ddx11dcz__dbcd1_z*dbcd1_z__dcz
     1            + ddx11dcz__ddabc1_z11dcz*ddabc1_z11dcz__dcz
     1            + ddx11dcz__dabc1_z*dabc1_z__dcz
     1            + ddx11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__dcz
      ddx11dcz__ddx=
     1              ddx11dcz__dbcd1_x*dbcd1_x__ddx
     1            + ddx11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__ddx
     1            + ddx11dcz__dbcd1_y*dbcd1_y__ddx
     1            + ddx11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__ddx
     1            + ddx11dcz__dbcd1_z*dbcd1_z__ddx
     1            + ddx11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__ddx
      ddx11dcz__ddy=
     1              ddx11dcz__dbcd1_x*dbcd1_x__ddy
     1            + ddx11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__ddy
     1            + ddx11dcz__dbcd1_y*dbcd1_y__ddy
     1            + ddx11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__ddy
     1            + ddx11dcz__dbcd1_z*dbcd1_z__ddy
     1            + ddx11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__ddy
      ddx11dcz__ddz=
     1              ddx11dcz__dbcd1_x*dbcd1_x__ddz
     1            + ddx11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__ddz
     1            + ddx11dcz__dbcd1_y*dbcd1_y__ddz
     1            + ddx11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__ddz
     1            + ddx11dcz__dbcd1_z*dbcd1_z__ddz
     1            + ddx11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__ddz
c      dx__ddx= abc1_x*dbcd1_x__ddx + abc1_y*dbcd1_y__ddx + abc1_z*dbcd1_z__ddx
      ddx11ddx__ddx=
     1            + ddx11ddx__ddbcd1_x11ddx*ddbcd1_x11ddx__ddx
     1            + ddx11ddx__ddbcd1_y11ddx*ddbcd1_y11ddx__ddx
     1            + ddx11ddx__ddbcd1_z11ddx*ddbcd1_z11ddx__ddx
      ddx11ddx__ddy=
     1            + ddx11ddx__ddbcd1_x11ddx*ddbcd1_x11ddx__ddy
     1            + ddx11ddx__ddbcd1_y11ddx*ddbcd1_y11ddx__ddy
     1            + ddx11ddx__ddbcd1_z11ddx*ddbcd1_z11ddx__ddy
      ddx11ddx__ddz=
     1            + ddx11ddx__ddbcd1_x11ddx*ddbcd1_x11ddx__ddz
     1            + ddx11ddx__ddbcd1_y11ddx*ddbcd1_y11ddx__ddz
     1            + ddx11ddx__ddbcd1_z11ddx*ddbcd1_z11ddx__ddz
c      dx__ddy= abc1_x*dbcd1_x__ddy + abc1_y*dbcd1_y__ddy + abc1_z*dbcd1_z__ddy
      ddx11ddy__ddy=
     1            + ddx11ddy__ddbcd1_x11ddy*ddbcd1_x11ddy__ddy
     1            + ddx11ddy__ddbcd1_y11ddy*ddbcd1_y11ddy__ddy
     1            + ddx11ddy__ddbcd1_z11ddy*ddbcd1_z11ddy__ddy
      ddx11ddy__ddz=
     1            + ddx11ddy__ddbcd1_x11ddy*ddbcd1_x11ddy__ddz
     1            + ddx11ddy__ddbcd1_y11ddy*ddbcd1_y11ddy__ddz
     1            + ddx11ddy__ddbcd1_z11ddy*ddbcd1_z11ddy__ddz
c      dx__ddz= abc1_x*dbcd1_x__ddz + abc1_y*dbcd1_y__ddz + abc1_z*dbcd1_z__ddz
      ddx11ddz__ddz=
     1            + ddx11ddz__ddbcd1_x11ddz*ddbcd1_x11ddz__ddz
     1            + ddx11ddz__ddbcd1_y11ddz*ddbcd1_y11ddz__ddz
     1            + ddx11ddz__ddbcd1_z11ddz*ddbcd1_z11ddz__ddz


c dy__dax= bcd1_x*daux_x__dax + bcd1_y*daux_y__dax + bcd1_z*daux_z__dax
      ddy11dax__dax=
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dax
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dax
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dax
      ddy11dax__day=
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__day
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__day
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__day
      ddy11dax__daz=
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__daz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__daz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__daz
      ddy11dax__dbx=
     1              ddy11dax__dbcd1_x*dbcd1_x__dbx
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dbx
     1            + ddy11dax__dbcd1_y*dbcd1_y__dbx
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dbx
     1            + ddy11dax__dbcd1_z*dbcd1_z__dbx
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dbx
      ddy11dax__dby=
     1              ddy11dax__dbcd1_x*dbcd1_x__dby
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dby
     1            + ddy11dax__dbcd1_y*dbcd1_y__dby
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dby
     1            + ddy11dax__dbcd1_z*dbcd1_z__dby
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dby
      ddy11dax__dbz=
     1              ddy11dax__dbcd1_x*dbcd1_x__dbz
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dbz
     1            + ddy11dax__dbcd1_y*dbcd1_y__dbz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dbz
     1            + ddy11dax__dbcd1_z*dbcd1_z__dbz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dbz
      ddy11dax__dcx=
     1              ddy11dax__dbcd1_x*dbcd1_x__dcx
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dcx
     1            + ddy11dax__dbcd1_y*dbcd1_y__dcx
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dcx
     1            + ddy11dax__dbcd1_z*dbcd1_z__dcx
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dcx
      ddy11dax__dcy=
     1              ddy11dax__dbcd1_x*dbcd1_x__dcy
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dcy
     1            + ddy11dax__dbcd1_y*dbcd1_y__dcy
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dcy
     1            + ddy11dax__dbcd1_z*dbcd1_z__dcy
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dcy
      ddy11dax__dcz=
     1              ddy11dax__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dax__ddaux_x11dax*ddaux_x11dax__dcz
     1            + ddy11dax__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dax__ddaux_y11dax*ddaux_y11dax__dcz
     1            + ddy11dax__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dax__ddaux_z11dax*ddaux_z11dax__dcz
      ddy11dax__ddx=
     1              ddy11dax__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dax__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dax__dbcd1_z*dbcd1_z__ddx
      ddy11dax__ddy=
     1              ddy11dax__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dax__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dax__dbcd1_z*dbcd1_z__ddy
      ddy11dax__ddz=
     1              ddy11dax__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dax__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dax__dbcd1_z*dbcd1_z__ddz
c dy__day= bcd1_x*daux_x__day + bcd1_y*daux_y__day + bcd1_z*daux_z__day
      ddy11day__day=
     1            + ddy11day__ddaux_x11day*ddaux_x11day__day
     1            + ddy11day__ddaux_y11day*ddaux_y11day__day
     1            + ddy11day__ddaux_z11day*ddaux_z11day__day
      ddy11day__daz=
     1            + ddy11day__ddaux_x11day*ddaux_x11day__daz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__daz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__daz
      ddy11day__dbx=
     1              ddy11day__dbcd1_x*dbcd1_x__dbx
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dbx
     1            + ddy11day__dbcd1_y*dbcd1_y__dbx
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dbx
     1            + ddy11day__dbcd1_z*dbcd1_z__dbx
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dbx
      ddy11day__dby=
     1              ddy11day__dbcd1_x*dbcd1_x__dby
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dby
     1            + ddy11day__dbcd1_y*dbcd1_y__dby
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dby
     1            + ddy11day__dbcd1_z*dbcd1_z__dby
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dby
      ddy11day__dbz=
     1              ddy11day__dbcd1_x*dbcd1_x__dbz
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dbz
     1            + ddy11day__dbcd1_y*dbcd1_y__dbz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dbz
     1            + ddy11day__dbcd1_z*dbcd1_z__dbz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dbz
      ddy11day__dcx=
     1              ddy11day__dbcd1_x*dbcd1_x__dcx
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dcx
     1            + ddy11day__dbcd1_y*dbcd1_y__dcx
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dcx
     1            + ddy11day__dbcd1_z*dbcd1_z__dcx
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dcx
      ddy11day__dcy=
     1              ddy11day__dbcd1_x*dbcd1_x__dcy
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dcy
     1            + ddy11day__dbcd1_y*dbcd1_y__dcy
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dcy
     1            + ddy11day__dbcd1_z*dbcd1_z__dcy
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dcy
      ddy11day__dcz=
     1              ddy11day__dbcd1_x*dbcd1_x__dcz
     1            + ddy11day__ddaux_x11day*ddaux_x11day__dcz
     1            + ddy11day__dbcd1_y*dbcd1_y__dcz
     1            + ddy11day__ddaux_y11day*ddaux_y11day__dcz
     1            + ddy11day__dbcd1_z*dbcd1_z__dcz
     1            + ddy11day__ddaux_z11day*ddaux_z11day__dcz
      ddy11day__ddx=
     1              ddy11day__dbcd1_x*dbcd1_x__ddx
     1            + ddy11day__dbcd1_y*dbcd1_y__ddx
     1            + ddy11day__dbcd1_z*dbcd1_z__ddx
      ddy11day__ddy=
     1              ddy11day__dbcd1_x*dbcd1_x__ddy
     1            + ddy11day__dbcd1_y*dbcd1_y__ddy
     1            + ddy11day__dbcd1_z*dbcd1_z__ddy
      ddy11day__ddz=
     1              ddy11day__dbcd1_x*dbcd1_x__ddz
     1            + ddy11day__dbcd1_y*dbcd1_y__ddz
     1            + ddy11day__dbcd1_z*dbcd1_z__ddz
c dy__daz= bcd1_x*daux_x__daz + bcd1_y*daux_y__daz + bcd1_z*daux_z__daz
      ddy11daz__daz=
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__daz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__daz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__daz
      ddy11daz__dbx=
     1              ddy11daz__dbcd1_x*dbcd1_x__dbx
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dbx
     1            + ddy11daz__dbcd1_y*dbcd1_y__dbx
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dbx
     1            + ddy11daz__dbcd1_z*dbcd1_z__dbx
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dbx
      ddy11daz__dby=
     1              ddy11daz__dbcd1_x*dbcd1_x__dby
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dby
     1            + ddy11daz__dbcd1_y*dbcd1_y__dby
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dby
     1            + ddy11daz__dbcd1_z*dbcd1_z__dby
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dby
      ddy11daz__dbz=
     1              ddy11daz__dbcd1_x*dbcd1_x__dbz
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dbz
     1            + ddy11daz__dbcd1_y*dbcd1_y__dbz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dbz
     1            + ddy11daz__dbcd1_z*dbcd1_z__dbz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dbz
      ddy11daz__dcx=
     1              ddy11daz__dbcd1_x*dbcd1_x__dcx
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dcx
     1            + ddy11daz__dbcd1_y*dbcd1_y__dcx
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dcx
     1            + ddy11daz__dbcd1_z*dbcd1_z__dcx
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dcx
      ddy11daz__dcy=
     1              ddy11daz__dbcd1_x*dbcd1_x__dcy
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dcy
     1            + ddy11daz__dbcd1_y*dbcd1_y__dcy
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dcy
     1            + ddy11daz__dbcd1_z*dbcd1_z__dcy
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dcy
      ddy11daz__dcz=
     1              ddy11daz__dbcd1_x*dbcd1_x__dcz
     1            + ddy11daz__ddaux_x11daz*ddaux_x11daz__dcz
     1            + ddy11daz__dbcd1_y*dbcd1_y__dcz
     1            + ddy11daz__ddaux_y11daz*ddaux_y11daz__dcz
     1            + ddy11daz__dbcd1_z*dbcd1_z__dcz
     1            + ddy11daz__ddaux_z11daz*ddaux_z11daz__dcz
      ddy11daz__ddx=
     1              ddy11daz__dbcd1_x*dbcd1_x__ddx
     1            + ddy11daz__dbcd1_y*dbcd1_y__ddx
     1            + ddy11daz__dbcd1_z*dbcd1_z__ddx
      ddy11daz__ddy=
     1              ddy11daz__dbcd1_x*dbcd1_x__ddy
     1            + ddy11daz__dbcd1_y*dbcd1_y__ddy
     1            + ddy11daz__dbcd1_z*dbcd1_z__ddy
      ddy11daz__ddz=
     1              ddy11daz__dbcd1_x*dbcd1_x__ddz
     1            + ddy11daz__dbcd1_y*dbcd1_y__ddz
     1            + ddy11daz__dbcd1_z*dbcd1_z__ddz
c dy__dbx= bcd1_x*daux_x__dbx + aux_x*dbcd1_x__dbx + bcd1_y*daux_y__dbx + aux_y*dbcd1_y__dbx + bcd1_z*daux_z__dbx + aux_z*dbcd1_z__dbx
      ddy11dbx__dbx=
     1              ddy11dbx__dbcd1_x*dbcd1_x__dbx
     1            + ddy11dbx__ddaux_x11dbx*ddaux_x11dbx__dbx
     1            + ddy11dbx__daux_x*daux_x__dbx
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dbx
     1            + ddy11dbx__dbcd1_y*dbcd1_y__dbx
     1            + ddy11dbx__ddaux_y11dbx*ddaux_y11dbx__dbx
     1            + ddy11dbx__daux_y*daux_y__dbx
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dbx
     1            + ddy11dbx__dbcd1_z*dbcd1_z__dbx
     1            + ddy11dbx__ddaux_z11dbx*ddaux_z11dbx__dbx
     1            + ddy11dbx__daux_z*daux_z__dbx
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dbx
      ddy11dbx__dby=
     1              ddy11dbx__dbcd1_x*dbcd1_x__dby
     1            + ddy11dbx__ddaux_x11dbx*ddaux_x11dbx__dby
     1            + ddy11dbx__daux_x*daux_x__dby
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dby
     1            + ddy11dbx__dbcd1_y*dbcd1_y__dby
     1            + ddy11dbx__ddaux_y11dbx*ddaux_y11dbx__dby
     1            + ddy11dbx__daux_y*daux_y__dby
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dby
     1            + ddy11dbx__dbcd1_z*dbcd1_z__dby
     1            + ddy11dbx__ddaux_z11dbx*ddaux_z11dbx__dby
     1            + ddy11dbx__daux_z*daux_z__dby
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dby
      ddy11dbx__dbz=
     1              ddy11dbx__dbcd1_x*dbcd1_x__dbz
     1            + ddy11dbx__ddaux_x11dbx*ddaux_x11dbx__dbz
     1            + ddy11dbx__daux_x*daux_x__dbz
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dbz
     1            + ddy11dbx__dbcd1_y*dbcd1_y__dbz
     1            + ddy11dbx__ddaux_y11dbx*ddaux_y11dbx__dbz
     1            + ddy11dbx__daux_y*daux_y__dbz
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dbz
     1            + ddy11dbx__dbcd1_z*dbcd1_z__dbz
     1            + ddy11dbx__ddaux_z11dbx*ddaux_z11dbx__dbz
     1            + ddy11dbx__daux_z*daux_z__dbz
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dbz
      ddy11dbx__dcx=
     1              ddy11dbx__dbcd1_x*dbcd1_x__dcx
     1            + ddy11dbx__ddaux_x11dbx*ddaux_x11dbx__dcx
     1            + ddy11dbx__daux_x*daux_x__dcx
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dcx
     1            + ddy11dbx__dbcd1_y*dbcd1_y__dcx
     1            + ddy11dbx__ddaux_y11dbx*ddaux_y11dbx__dcx
     1            + ddy11dbx__daux_y*daux_y__dcx
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dcx
     1            + ddy11dbx__dbcd1_z*dbcd1_z__dcx
     1            + ddy11dbx__ddaux_z11dbx*ddaux_z11dbx__dcx
     1            + ddy11dbx__daux_z*daux_z__dcx
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dcx
      ddy11dbx__dcy=
     1              ddy11dbx__dbcd1_x*dbcd1_x__dcy
     1            + ddy11dbx__ddaux_x11dbx*ddaux_x11dbx__dcy
     1            + ddy11dbx__daux_x*daux_x__dcy
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dcy
     1            + ddy11dbx__dbcd1_y*dbcd1_y__dcy
     1            + ddy11dbx__ddaux_y11dbx*ddaux_y11dbx__dcy
     1            + ddy11dbx__daux_y*daux_y__dcy
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dcy
     1            + ddy11dbx__dbcd1_z*dbcd1_z__dcy
     1            + ddy11dbx__ddaux_z11dbx*ddaux_z11dbx__dcy
     1            + ddy11dbx__daux_z*daux_z__dcy
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dcy
      ddy11dbx__dcz=
     1              ddy11dbx__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dbx__ddaux_x11dbx*ddaux_x11dbx__dcz
     1            + ddy11dbx__daux_x*daux_x__dcz
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__dcz
     1            + ddy11dbx__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dbx__ddaux_y11dbx*ddaux_y11dbx__dcz
     1            + ddy11dbx__daux_y*daux_y__dcz
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__dcz
     1            + ddy11dbx__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dbx__ddaux_z11dbx*ddaux_z11dbx__dcz
     1            + ddy11dbx__daux_z*daux_z__dcz
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__dcz
      ddy11dbx__ddx=
     1              ddy11dbx__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__ddx
     1            + ddy11dbx__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__ddx
     1            + ddy11dbx__dbcd1_z*dbcd1_z__ddx
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__ddx
      ddy11dbx__ddy=
     1              ddy11dbx__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__ddy
     1            + ddy11dbx__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__ddy
     1            + ddy11dbx__dbcd1_z*dbcd1_z__ddy
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__ddy
      ddy11dbx__ddz=
     1              ddy11dbx__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dbx__ddbcd1_x11dbx*ddbcd1_x11dbx__ddz
     1            + ddy11dbx__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dbx__ddbcd1_y11dbx*ddbcd1_y11dbx__ddz
     1            + ddy11dbx__dbcd1_z*dbcd1_z__ddz
     1            + ddy11dbx__ddbcd1_z11dbx*ddbcd1_z11dbx__ddz
c dy__dby= bcd1_x*daux_x__dby + aux_x*dbcd1_x__dby + bcd1_y*daux_y__dby + aux_y*dbcd1_y__dby + bcd1_z*daux_z__dby + aux_z*dbcd1_z__dby
      ddy11dby__dby=
     1              ddy11dby__dbcd1_x*dbcd1_x__dby
     1            + ddy11dby__ddaux_x11dby*ddaux_x11dby__dby
     1            + ddy11dby__daux_x*daux_x__dby
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__dby
     1            + ddy11dby__dbcd1_y*dbcd1_y__dby
     1            + ddy11dby__ddaux_y11dby*ddaux_y11dby__dby
     1            + ddy11dby__daux_y*daux_y__dby
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__dby
     1            + ddy11dby__dbcd1_z*dbcd1_z__dby
     1            + ddy11dby__ddaux_z11dby*ddaux_z11dby__dby
     1            + ddy11dby__daux_z*daux_z__dby
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__dby
      ddy11dby__dbz=
     1              ddy11dby__dbcd1_x*dbcd1_x__dbz
     1            + ddy11dby__ddaux_x11dby*ddaux_x11dby__dbz
     1            + ddy11dby__daux_x*daux_x__dbz
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__dbz
     1            + ddy11dby__dbcd1_y*dbcd1_y__dbz
     1            + ddy11dby__ddaux_y11dby*ddaux_y11dby__dbz
     1            + ddy11dby__daux_y*daux_y__dbz
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__dbz
     1            + ddy11dby__dbcd1_z*dbcd1_z__dbz
     1            + ddy11dby__ddaux_z11dby*ddaux_z11dby__dbz
     1            + ddy11dby__daux_z*daux_z__dbz
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__dbz
      ddy11dby__dcx=
     1              ddy11dby__dbcd1_x*dbcd1_x__dcx
     1            + ddy11dby__ddaux_x11dby*ddaux_x11dby__dcx
     1            + ddy11dby__daux_x*daux_x__dcx
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__dcx
     1            + ddy11dby__dbcd1_y*dbcd1_y__dcx
     1            + ddy11dby__ddaux_y11dby*ddaux_y11dby__dcx
     1            + ddy11dby__daux_y*daux_y__dcx
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__dcx
     1            + ddy11dby__dbcd1_z*dbcd1_z__dcx
     1            + ddy11dby__ddaux_z11dby*ddaux_z11dby__dcx
     1            + ddy11dby__daux_z*daux_z__dcx
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__dcx
      ddy11dby__dcy=
     1              ddy11dby__dbcd1_x*dbcd1_x__dcy
     1            + ddy11dby__ddaux_x11dby*ddaux_x11dby__dcy
     1            + ddy11dby__daux_x*daux_x__dcy
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__dcy
     1            + ddy11dby__dbcd1_y*dbcd1_y__dcy
     1            + ddy11dby__ddaux_y11dby*ddaux_y11dby__dcy
     1            + ddy11dby__daux_y*daux_y__dcy
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__dcy
     1            + ddy11dby__dbcd1_z*dbcd1_z__dcy
     1            + ddy11dby__ddaux_z11dby*ddaux_z11dby__dcy
     1            + ddy11dby__daux_z*daux_z__dcy
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__dcy
      ddy11dby__dcz=
     1              ddy11dby__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dby__ddaux_x11dby*ddaux_x11dby__dcz
     1            + ddy11dby__daux_x*daux_x__dcz
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__dcz
     1            + ddy11dby__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dby__ddaux_y11dby*ddaux_y11dby__dcz
     1            + ddy11dby__daux_y*daux_y__dcz
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__dcz
     1            + ddy11dby__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dby__ddaux_z11dby*ddaux_z11dby__dcz
     1            + ddy11dby__daux_z*daux_z__dcz
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__dcz
      ddy11dby__ddx=
     1              ddy11dby__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__ddx
     1            + ddy11dby__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__ddx
     1            + ddy11dby__dbcd1_z*dbcd1_z__ddx
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__ddx
      ddy11dby__ddy=
     1              ddy11dby__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__ddy
     1            + ddy11dby__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__ddy
     1            + ddy11dby__dbcd1_z*dbcd1_z__ddy
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__ddy
      ddy11dby__ddz=
     1              ddy11dby__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dby__ddbcd1_x11dby*ddbcd1_x11dby__ddz
     1            + ddy11dby__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dby__ddbcd1_y11dby*ddbcd1_y11dby__ddz
     1            + ddy11dby__dbcd1_z*dbcd1_z__ddz
     1            + ddy11dby__ddbcd1_z11dby*ddbcd1_z11dby__ddz
c dy__dbz= bcd1_x*daux_x__dbz + aux_x*dbcd1_x__dbz + bcd1_y*daux_y__dbz + aux_y*dbcd1_y__dbz + bcd1_z*daux_z__dbz + aux_z*dbcd1_z__dbz
      ddy11dbz__dbz=
     1              ddy11dbz__dbcd1_x*dbcd1_x__dbz
     1            + ddy11dbz__ddaux_x11dbz*ddaux_x11dbz__dbz
     1            + ddy11dbz__daux_x*daux_x__dbz
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dbz
     1            + ddy11dbz__dbcd1_y*dbcd1_y__dbz
     1            + ddy11dbz__ddaux_y11dbz*ddaux_y11dbz__dbz
     1            + ddy11dbz__daux_y*daux_y__dbz
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dbz
     1            + ddy11dbz__dbcd1_z*dbcd1_z__dbz
     1            + ddy11dbz__ddaux_z11dbz*ddaux_z11dbz__dbz
     1            + ddy11dbz__daux_z*daux_z__dbz
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dbz
      ddy11dbz__dcx=
     1              ddy11dbz__dbcd1_x*dbcd1_x__dcx
     1            + ddy11dbz__ddaux_x11dbz*ddaux_x11dbz__dcx
     1            + ddy11dbz__daux_x*daux_x__dcx
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dcx
     1            + ddy11dbz__dbcd1_y*dbcd1_y__dcx
     1            + ddy11dbz__ddaux_y11dbz*ddaux_y11dbz__dcx
     1            + ddy11dbz__daux_y*daux_y__dcx
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dcx
     1            + ddy11dbz__dbcd1_z*dbcd1_z__dcx
     1            + ddy11dbz__ddaux_z11dbz*ddaux_z11dbz__dcx
     1            + ddy11dbz__daux_z*daux_z__dcx
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dcx
      ddy11dbz__dcy=
     1              ddy11dbz__dbcd1_x*dbcd1_x__dcy
     1            + ddy11dbz__ddaux_x11dbz*ddaux_x11dbz__dcy
     1            + ddy11dbz__daux_x*daux_x__dcy
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dcy
     1            + ddy11dbz__dbcd1_y*dbcd1_y__dcy
     1            + ddy11dbz__ddaux_y11dbz*ddaux_y11dbz__dcy
     1            + ddy11dbz__daux_y*daux_y__dcy
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dcy
     1            + ddy11dbz__dbcd1_z*dbcd1_z__dcy
     1            + ddy11dbz__ddaux_z11dbz*ddaux_z11dbz__dcy
     1            + ddy11dbz__daux_z*daux_z__dcy
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dcy
      ddy11dbz__dcz=
     1              ddy11dbz__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dbz__ddaux_x11dbz*ddaux_x11dbz__dcz
     1            + ddy11dbz__daux_x*daux_x__dcz
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__dcz
     1            + ddy11dbz__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dbz__ddaux_y11dbz*ddaux_y11dbz__dcz
     1            + ddy11dbz__daux_y*daux_y__dcz
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__dcz
     1            + ddy11dbz__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dbz__ddaux_z11dbz*ddaux_z11dbz__dcz
     1            + ddy11dbz__daux_z*daux_z__dcz
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__dcz
      ddy11dbz__ddx=
     1              ddy11dbz__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__ddx
     1            + ddy11dbz__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__ddx
     1            + ddy11dbz__dbcd1_z*dbcd1_z__ddx
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__ddx
      ddy11dbz__ddy=
     1              ddy11dbz__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__ddy
     1            + ddy11dbz__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__ddy
     1            + ddy11dbz__dbcd1_z*dbcd1_z__ddy
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__ddy
      ddy11dbz__ddz=
     1              ddy11dbz__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dbz__ddbcd1_x11dbz*ddbcd1_x11dbz__ddz
     1            + ddy11dbz__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dbz__ddbcd1_y11dbz*ddbcd1_y11dbz__ddz
     1            + ddy11dbz__dbcd1_z*dbcd1_z__ddz
     1            + ddy11dbz__ddbcd1_z11dbz*ddbcd1_z11dbz__ddz
c dy__dcx= bcd1_x*daux_x__dcx + aux_x*dbcd1_x__dcx + bcd1_y*daux_y__dcx + aux_y*dbcd1_y__dcx + bcd1_z*daux_z__dcx + aux_z*dbcd1_z__dcx
      ddy11dcx__dcx=
     1              ddy11dcx__dbcd1_x*dbcd1_x__dcx
     1            + ddy11dcx__ddaux_x11dcx*ddaux_x11dcx__dcx
     1            + ddy11dcx__daux_x*daux_x__dcx
     1            + ddy11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__dcx
     1            + ddy11dcx__dbcd1_y*dbcd1_y__dcx
     1            + ddy11dcx__ddaux_y11dcx*ddaux_y11dcx__dcx
     1            + ddy11dcx__daux_y*daux_y__dcx
     1            + ddy11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__dcx
     1            + ddy11dcx__dbcd1_z*dbcd1_z__dcx
     1            + ddy11dcx__ddaux_z11dcx*ddaux_z11dcx__dcx
     1            + ddy11dcx__daux_z*daux_z__dcx
     1            + ddy11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__dcx
      ddy11dcx__dcy=
     1              ddy11dcx__dbcd1_x*dbcd1_x__dcy
     1            + ddy11dcx__ddaux_x11dcx*ddaux_x11dcx__dcy
     1            + ddy11dcx__daux_x*daux_x__dcy
     1            + ddy11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__dcy
     1            + ddy11dcx__dbcd1_y*dbcd1_y__dcy
     1            + ddy11dcx__ddaux_y11dcx*ddaux_y11dcx__dcy
     1            + ddy11dcx__daux_y*daux_y__dcy
     1            + ddy11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__dcy
     1            + ddy11dcx__dbcd1_z*dbcd1_z__dcy
     1            + ddy11dcx__ddaux_z11dcx*ddaux_z11dcx__dcy
     1            + ddy11dcx__daux_z*daux_z__dcy
     1            + ddy11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__dcy
      ddy11dcx__dcz=
     1              ddy11dcx__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dcx__ddaux_x11dcx*ddaux_x11dcx__dcz
     1            + ddy11dcx__daux_x*daux_x__dcz
     1            + ddy11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__dcz
     1            + ddy11dcx__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dcx__ddaux_y11dcx*ddaux_y11dcx__dcz
     1            + ddy11dcx__daux_y*daux_y__dcz
     1            + ddy11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__dcz
     1            + ddy11dcx__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dcx__ddaux_z11dcx*ddaux_z11dcx__dcz
     1            + ddy11dcx__daux_z*daux_z__dcz
     1            + ddy11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__dcz
      ddy11dcx__ddx=
     1              ddy11dcx__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__ddx
     1            + ddy11dcx__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__ddx
     1            + ddy11dcx__dbcd1_z*dbcd1_z__ddx
     1            + ddy11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__ddx
      ddy11dcx__ddy=
     1              ddy11dcx__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__ddy
     1            + ddy11dcx__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__ddy
     1            + ddy11dcx__dbcd1_z*dbcd1_z__ddy
     1            + ddy11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__ddy
      ddy11dcx__ddz=
     1              ddy11dcx__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dcx__ddbcd1_x11dcx*ddbcd1_x11dcx__ddz
     1            + ddy11dcx__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dcx__ddbcd1_y11dcx*ddbcd1_y11dcx__ddz
     1            + ddy11dcx__dbcd1_z*dbcd1_z__ddz
     1            + ddy11dcx__ddbcd1_z11dcx*ddbcd1_z11dcx__ddz
c dy__dcy= bcd1_x*daux_x__dcy + aux_x*dbcd1_x__dcy + bcd1_y*daux_y__dcy + aux_y*dbcd1_y__dcy + bcd1_z*daux_z__dcy + aux_z*dbcd1_z__dcy
      ddy11dcy__dcy=
     1              ddy11dcy__dbcd1_x*dbcd1_x__dcy
     1            + ddy11dcy__ddaux_x11dcy*ddaux_x11dcy__dcy
     1            + ddy11dcy__daux_x*daux_x__dcy
     1            + ddy11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__dcy
     1            + ddy11dcy__dbcd1_y*dbcd1_y__dcy
     1            + ddy11dcy__ddaux_y11dcy*ddaux_y11dcy__dcy
     1            + ddy11dcy__daux_y*daux_y__dcy
     1            + ddy11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__dcy
     1            + ddy11dcy__dbcd1_z*dbcd1_z__dcy
     1            + ddy11dcy__ddaux_z11dcy*ddaux_z11dcy__dcy
     1            + ddy11dcy__daux_z*daux_z__dcy
     1            + ddy11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__dcy
      ddy11dcy__dcz=
     1              ddy11dcy__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dcy__ddaux_x11dcy*ddaux_x11dcy__dcz
     1            + ddy11dcy__daux_x*daux_x__dcz
     1            + ddy11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__dcz
     1            + ddy11dcy__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dcy__ddaux_y11dcy*ddaux_y11dcy__dcz
     1            + ddy11dcy__daux_y*daux_y__dcz
     1            + ddy11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__dcz
     1            + ddy11dcy__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dcy__ddaux_z11dcy*ddaux_z11dcy__dcz
     1            + ddy11dcy__daux_z*daux_z__dcz
     1            + ddy11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__dcz
      ddy11dcy__ddx=
     1              ddy11dcy__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__ddx
     1            + ddy11dcy__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__ddx
     1            + ddy11dcy__dbcd1_z*dbcd1_z__ddx
     1            + ddy11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__ddx
      ddy11dcy__ddy=
     1              ddy11dcy__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__ddy
     1            + ddy11dcy__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__ddy
     1            + ddy11dcy__dbcd1_z*dbcd1_z__ddy
     1            + ddy11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__ddy
      ddy11dcy__ddz=
     1              ddy11dcy__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dcy__ddbcd1_x11dcy*ddbcd1_x11dcy__ddz
     1            + ddy11dcy__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dcy__ddbcd1_y11dcy*ddbcd1_y11dcy__ddz
     1            + ddy11dcy__dbcd1_z*dbcd1_z__ddz
     1            + ddy11dcy__ddbcd1_z11dcy*ddbcd1_z11dcy__ddz
c dy__dcz= bcd1_x*daux_x__dcz + aux_x*dbcd1_x__dcz + bcd1_y*daux_y__dcz + aux_y*dbcd1_y__dcz + bcd1_z*daux_z__dcz + aux_z*dbcd1_z__dcz
      ddy11dcz__dcz=
     1              ddy11dcz__dbcd1_x*dbcd1_x__dcz
     1            + ddy11dcz__ddaux_x11dcz*ddaux_x11dcz__dcz
     1            + ddy11dcz__daux_x*daux_x__dcz
     1            + ddy11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__dcz
     1            + ddy11dcz__dbcd1_y*dbcd1_y__dcz
     1            + ddy11dcz__ddaux_y11dcz*ddaux_y11dcz__dcz
     1            + ddy11dcz__daux_y*daux_y__dcz
     1            + ddy11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__dcz
     1            + ddy11dcz__dbcd1_z*dbcd1_z__dcz
     1            + ddy11dcz__ddaux_z11dcz*ddaux_z11dcz__dcz
     1            + ddy11dcz__daux_z*daux_z__dcz
     1            + ddy11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__dcz
      ddy11dcz__ddx=
     1              ddy11dcz__dbcd1_x*dbcd1_x__ddx
     1            + ddy11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__ddx
     1            + ddy11dcz__dbcd1_y*dbcd1_y__ddx
     1            + ddy11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__ddx
     1            + ddy11dcz__dbcd1_z*dbcd1_z__ddx
     1            + ddy11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__ddx
      ddy11dcz__ddy=
     1              ddy11dcz__dbcd1_x*dbcd1_x__ddy
     1            + ddy11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__ddy
     1            + ddy11dcz__dbcd1_y*dbcd1_y__ddy
     1            + ddy11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__ddy
     1            + ddy11dcz__dbcd1_z*dbcd1_z__ddy
     1            + ddy11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__ddy
      ddy11dcz__ddz=
     1              ddy11dcz__dbcd1_x*dbcd1_x__ddz
     1            + ddy11dcz__ddbcd1_x11dcz*ddbcd1_x11dcz__ddz
     1            + ddy11dcz__dbcd1_y*dbcd1_y__ddz
     1            + ddy11dcz__ddbcd1_y11dcz*ddbcd1_y11dcz__ddz
     1            + ddy11dcz__dbcd1_z*dbcd1_z__ddz
     1            + ddy11dcz__ddbcd1_z11dcz*ddbcd1_z11dcz__ddz
c dy__ddx= aux_x*dbcd1_x__ddx + aux_y*dbcd1_y__ddx + aux_z*dbcd1_z__ddx
      ddy11ddx__ddx=
     1            + ddy11ddx__ddbcd1_x11ddx*ddbcd1_x11ddx__ddx
     1            + ddy11ddx__ddbcd1_y11ddx*ddbcd1_y11ddx__ddx
     1            + ddy11ddx__ddbcd1_z11ddx*ddbcd1_z11ddx__ddx
      ddy11ddx__ddy=
     1            + ddy11ddx__ddbcd1_x11ddx*ddbcd1_x11ddx__ddy
     1            + ddy11ddx__ddbcd1_y11ddx*ddbcd1_y11ddx__ddy
     1            + ddy11ddx__ddbcd1_z11ddx*ddbcd1_z11ddx__ddy
      ddy11ddx__ddz=
     1            + ddy11ddx__ddbcd1_x11ddx*ddbcd1_x11ddx__ddz
     1            + ddy11ddx__ddbcd1_y11ddx*ddbcd1_y11ddx__ddz
     1            + ddy11ddx__ddbcd1_z11ddx*ddbcd1_z11ddx__ddz
c dy__ddy= aux_x*dbcd1_x__ddy + aux_y*dbcd1_y__ddy + aux_z*dbcd1_z__ddy
      ddy11ddy__ddy=
     1            + ddy11ddy__ddbcd1_x11ddy*ddbcd1_x11ddy__ddy
     1            + ddy11ddy__ddbcd1_y11ddy*ddbcd1_y11ddy__ddy
     1            + ddy11ddy__ddbcd1_z11ddy*ddbcd1_z11ddy__ddy
      ddy11ddy__ddz=
     1            + ddy11ddy__ddbcd1_x11ddy*ddbcd1_x11ddy__ddz
     1            + ddy11ddy__ddbcd1_y11ddy*ddbcd1_y11ddy__ddz
     1            + ddy11ddy__ddbcd1_z11ddy*ddbcd1_z11ddy__ddz
c dy__ddz= aux_x*dbcd1_x__ddz + aux_y*dbcd1_y__ddz + aux_z*dbcd1_z__ddz
      ddy11ddz__ddz=
     1            + ddy11ddz__ddbcd1_x11ddz*ddbcd1_x11ddz__ddz
     1            + ddy11ddz__ddbcd1_y11ddz*ddbcd1_y11ddz__ddz
     1            + ddy11ddz__ddbcd1_z11ddz*ddbcd1_z11ddz__ddz


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c block 1

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
c the 9 derivatives with respect to a{xyz} and d{xyz} are 0
      ddf11dax__ddx=0
c     1   ddf11dax__ddf11dx*ddf11dx__ddx
c     1 + ddf11dax__ddx11dax*ddx11dax__ddx
c     1 + ddf11dax__ddf11dy*ddf11dy__ddx
c     1 + ddf11dax__ddy11dax*ddy11dax__ddx
      ddf11dax__ddy=0
c     1   ddf11dax__ddf11dx*ddf11dx__ddy
c     1 + ddf11dax__ddx11dax*ddx11dax__ddy
c     1 + ddf11dax__ddf11dy*ddf11dy__ddy
c     1 + ddf11dax__ddy11dax*ddy11dax__ddy
      ddf11dax__ddz=0
c     1   ddf11dax__ddf11dx*ddf11dx__ddz
c     1 + ddf11dax__ddx11dax*ddx11dax__ddz
c     1 + ddf11dax__ddf11dy*ddf11dy__ddz
c     1 + ddf11dax__ddy11dax*ddy11dax__ddz
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
      ddf11day__ddx=0
c     1   ddf11day__ddf11dx*ddf11dx__ddx
c     1 + ddf11day__ddx11day*ddx11day__ddx
c     1 + ddf11day__ddf11dy*ddf11dy__ddx
c     1 + ddf11day__ddy11day*ddy11day__ddx
      ddf11day__ddy=0
c     1   ddf11day__ddf11dx*ddf11dx__ddy
c     1 + ddf11day__ddx11day*ddx11day__ddy
c     1 + ddf11day__ddf11dy*ddf11dy__ddy
c     1 + ddf11day__ddy11day*ddy11day__ddy
      ddf11day__ddz=0
c     1   ddf11day__ddf11dx*ddf11dx__ddz
c     1 + ddf11day__ddx11day*ddx11day__ddz
c     1 + ddf11day__ddf11dy*ddf11dy__ddz
c     1 + ddf11day__ddy11day*ddy11day__ddz
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
      ddf11daz__ddx=0
c     1   ddf11daz__ddf11dx*ddf11dx__ddx
c     1 + ddf11daz__ddx11daz*ddx11daz__ddx
c     1 + ddf11daz__ddf11dy*ddf11dy__ddx
c     1 + ddf11daz__ddy11daz*ddy11daz__ddx
      ddf11daz__ddy=0
c     1   ddf11daz__ddf11dx*ddf11dx__ddy
c     1 + ddf11daz__ddx11daz*ddx11daz__ddy
c     1 + ddf11daz__ddf11dy*ddf11dy__ddy
c     1 + ddf11daz__ddy11daz*ddy11daz__ddy
      ddf11daz__ddz=0
c     1   ddf11daz__ddf11dx*ddf11dx__ddz
c     1 + ddf11daz__ddx11daz*ddx11daz__ddz
c     1 + ddf11daz__ddf11dy*ddf11dy__ddz
c     1 + ddf11daz__ddy11daz*ddy11daz__ddz
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

c reverse the signs
      axax=-ddf11dax__dax
      axay=-ddf11dax__day
      axaz=-ddf11dax__daz
      axbx=-ddf11dax__dbx
      axby=-ddf11dax__dby
      axbz=-ddf11dax__dbz
      axcx=-ddf11dax__dcx
      axcy=-ddf11dax__dcy
      axcz=-ddf11dax__dcz
      axdx=-ddf11dax__ddx
      axdy=-ddf11dax__ddy
      axdz=-ddf11dax__ddz

      ayay=-ddf11day__day
      ayaz=-ddf11day__daz
      aybx=-ddf11day__dbx
      ayby=-ddf11day__dby
      aybz=-ddf11day__dbz
      aycx=-ddf11day__dcx
      aycy=-ddf11day__dcy
      aycz=-ddf11day__dcz
      aydx=-ddf11day__ddx
      aydy=-ddf11day__ddy
      aydz=-ddf11day__ddz

      azaz=-ddf11daz__daz
      azbx=-ddf11daz__dbx
      azby=-ddf11daz__dby
      azbz=-ddf11daz__dbz
      azcx=-ddf11daz__dcx
      azcy=-ddf11daz__dcy
      azcz=-ddf11daz__dcz
      azdx=-ddf11daz__ddx
      azdy=-ddf11daz__ddy
      azdz=-ddf11daz__ddz

      bxbx=-ddf11dbx__dbx
      bxby=-ddf11dbx__dby
      bxbz=-ddf11dbx__dbz
      bxcx=-ddf11dbx__dcx
      bxcy=-ddf11dbx__dcy
      bxcz=-ddf11dbx__dcz
      bxdx=-ddf11dbx__ddx
      bxdy=-ddf11dbx__ddy
      bxdz=-ddf11dbx__ddz

      byby=-ddf11dby__dby
      bybz=-ddf11dby__dbz
      bycx=-ddf11dby__dcx
      bycy=-ddf11dby__dcy
      bycz=-ddf11dby__dcz
      bydx=-ddf11dby__ddx
      bydy=-ddf11dby__ddy
      bydz=-ddf11dby__ddz

      bzbz=-ddf11dbz__dbz
      bzcx=-ddf11dbz__dcx
      bzcy=-ddf11dbz__dcy
      bzcz=-ddf11dbz__dcz
      bzdx=-ddf11dbz__ddx
      bzdy=-ddf11dbz__ddy
      bzdz=-ddf11dbz__ddz

      cxcx=-ddf11dcx__dcx
      cxcy=-ddf11dcx__dcy
      cxcz=-ddf11dcx__dcz
      cxdx=-ddf11dcx__ddx
      cxdy=-ddf11dcx__ddy
      cxdz=-ddf11dcx__ddz

      cycy=-ddf11dcy__dcy
      cycz=-ddf11dcy__dcz
      cydx=-ddf11dcy__ddx
      cydy=-ddf11dcy__ddy
      cydz=-ddf11dcy__ddz

      czcz=-ddf11dcz__dcz
      czdx=-ddf11dcz__ddx
      czdy=-ddf11dcz__ddy
      czdz=-ddf11dcz__ddz

      dxdx=-ddf11ddx__ddx
      dxdy=-ddf11ddx__ddy
      dxdz=-ddf11ddx__ddz

      dydy=-ddf11ddy__ddy
      dydz=-ddf11ddy__ddz

      dzdz=-ddf11ddz__ddz

      return
      END SUBROUTINE DDDIHEDRAL
