c subroutine dist takes 6 reals (=2 coordinates) and yields a positive distance
      SUBROUTINE DIST(ax,ay,az,bx,by,bz,dist_ab)
      implicit real(8) (a-z)
      dist_ab=dsqrt((bx-ax)**2 + (by-ay)**2 + (bz-az)**2)
      return
      END SUBROUTINE DIST


c subroutine ddist takes 6 reals (=2 coordinates) and yields all 6 first derivatives of the distance
      SUBROUTINE DDIST(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,
     2  dist_ab)
      implicit real(8) (a-z)
      ab_x=bx-ax
      ab_y=by-ay
      ab_z=bz-az
      dist_ab=dsqrt((ab_x)**2 + (ab_y)**2 + (ab_z)**2)
      dist_ab_inv=1.0/dist_ab
      dax=-ab_x*dist_ab_inv
      dbx=-dax
      day=-ab_y*dist_ab_inv
      dby=-day
      daz=-ab_z*dist_ab_inv
      dbz=-daz
      return
      END SUBROUTINE DDIST


c subroutine dddist takes 6 reals (=2 coordinates) and yields all 6 first derivatives,
c all 21 second derivatives of the distance and the the distance itself
      SUBROUTINE DDDIST(ax,ay,az,bx,by,bz,
     2 df__dax,df__day,df__daz,df__dbx,df__dby,df__dbz,
     3 ddf11dax__dax,ddf11dax__day,ddf11dax__daz,ddf11dax__dbx,
     4 ddf11dax__dby,ddf11dax__dbz,ddf11day__day,ddf11day__daz,
     5 ddf11day__dbx,ddf11day__dby,ddf11day__dbz,ddf11daz__daz,
     6 ddf11daz__dbx,ddf11daz__dby,ddf11daz__dbz,ddf11dbx__dbx,
     7 ddf11dbx__dby,ddf11dbx__dbz,ddf11dby__dby,ddf11dby__dbz,
     5 ddf11dbz__dbz,
     6 f)
      implicit real(8) (a-z)
      ab_x=bx-ax
      ab_y=by-ay
      ab_z=bz-az
      f=dsqrt((ab_x)**2 + (ab_y)**2 + (ab_z)**2)
      f_inv=1/f
      f_inv_cub=f_inv**3

c first derivatives
      df__dax=-ab_x*f_inv
      df__dbx= ab_x*f_inv
      df__day=-ab_y*f_inv
      df__dby= ab_y*f_inv
      df__daz=-ab_z*f_inv
      df__dbz= ab_z*f_inv

c second derivatives 
c f_inv=1/dsqrt((ab_x)**2 + (ab_y)**2 + (ab_z)**2)
      df_inv__dab_x=-ab_x*f_inv_cub
      df_inv__dab_y=-ab_y*f_inv_cub
      df_inv__dab_z=-ab_z*f_inv_cub

c df__dax=df__dab_x*dab_x__dax
      ddf11dax__dax= f_inv + ab_x*df_inv__dab_x
      ddf11dax__day= ab_x*df_inv__dab_y
      ddf11dax__daz= ab_x*df_inv__dab_z
      ddf11dax__dbx=-f_inv - ab_x*df_inv__dab_x
      ddf11dax__dby=-ab_x*df_inv__dab_y
      ddf11dax__dbz=-ab_x*df_inv__dab_z

c df__day=df__dab_y*dab_y__day
      ddf11day__day= f_inv + ab_y*df_inv__dab_y
      ddf11day__daz= ab_y*df_inv__dab_z
      ddf11day__dbx=-ab_y*df_inv__dab_x
      ddf11day__dby=-f_inv - ab_y*df_inv__dab_y
      ddf11day__dbz=-ab_y*df_inv__dab_z

c df__daz=df__dab_z*dab_z__daz
      ddf11daz__daz= f_inv + ab_z*df_inv__dab_z
      ddf11daz__dbx=-ab_z*df_inv__dab_x
      ddf11daz__dby=-ab_z*df_inv__dab_y
      ddf11daz__dbz=-f_inv - ab_z*df_inv__dab_z

c df__dbx=df__dab_x*dab_x__dbx
      ddf11dbx__dbx= f_inv + ab_x*df_inv__dab_x
      ddf11dbx__dby= ab_x*df_inv__dab_y
      ddf11dbx__dbz= ab_x*df_inv__dab_z

c df__dby=df__dab_y*dab_y__dby
      ddf11dby__dby= f_inv + ab_y*df_inv__dab_y
      ddf11dby__dbz= ab_y*df_inv__dab_z

c df__dbz=df__dab_z*dab_z__dbz
      ddf11dbz__dbz= f_inv + ab_z*df_inv__dab_z

      return
      END SUBROUTINE DDDIST


c subroutine angle takes 9 reals (=3 coordinates) and yields an angel between 0 and +\pi (in radians)
c via law of cosines
c mult: 11, div: 1, root: 2, add/sub: 8, arccos: 1
      SUBROUTINE ANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      implicit real(8) (a-z)
      r2L=(ax-bx)**2 + (ay-by)**2 + (az-bz)**2
      r2M=(ax-cx)**2 + (ay-cy)**2 + (az-cz)**2
      r2R=(bx-cx)**2 + (by-cy)**2 + (bz-cz)**2
      den=2.0 * dsqrt(r2L * r2R)
      arg=(r2L+r2R-r2M)/den
c the following two exceptions may be called in case of rounding errors
      arg = min(arg, 1.d0)
      arg = max(arg, -1.d0)
c      if(arg .gt. 1.d0) arg=1.d0
c      if(arg .lt. -1.d0) arg=-1.d0
      angle_abc=dacos(arg)
      return
      END SUBROUTINE ANGLE


c subroutine angle takes 9 reals (=3 coordinates) and yields an angel between 0 and +\pi (in radians)
c via vector definition of the cosine
c mult: 10, div: 1, root: 2, add/sub: 12, arccos: 1
c      SUBROUTINE ANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
c      implicit real(8) (a-z)
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


c subroutine dangle takes 9 reals (=3 coordinates) and yields all 9 first derivatives of the angle
c via law of cosines (calculating the derivative of Abs[foo] is rather troublesome)
      SUBROUTINE DANGLE(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     3 angle_abc)
      implicit real(8) (a-z)
c the direction of vectors doesn't matter here, only the length does.
      ab_x=bx-ax
      ab_y=by-ay
      ab_z=bz-az
      bc_x=cx-bx
      bc_y=cy-by
      bc_z=cz-bz
      ca_x=ax-cx
      ca_y=ay-cy
      ca_z=az-cz
c length of a-b, b-c and a-c
      r2L=ab_x**2 + ab_y**2 + ab_z**2
      r2R=bc_x**2 + bc_y**2 + bc_z**2
      r2M=ca_x**2 + ca_y**2 + ca_z**2
      r1L=dsqrt(r2L)
      r1R=dsqrt(r2R)
      r1M=dsqrt(r2M)
      r3L=r2L*r1L
      r3R=r2R*r1R
      r3M=r2M*r1M
c some auxiliary products
      aux_11_inv=1.d0/(2.d0*r1L*r1R)
      aux_31_inv=1.d0/(2.d0*r3L*r1R)
      aux_13_inv=1.d0/(2.d0*r1L*r3R)
      aux_1=r2L + r2R - r2M
      arccos_arg=aux_1*aux_11_inv
c if the angle is defined as the smallest angle (which it should be) there is a singularity
c at 180 degrees.  Depending on the coordinate system the derivative is between 1 and -1
c (in case of sides of length 1) 
      if(arccos_arg .ge. 1.d0) then
        angle_abc=dacos(1.d0)
        dax=0
        day=0
        daz=0
        dbx=0
        dby=0
        dbz=0
        dcx=0
        dcy=0
        dcz=0
        return
      elseif(arccos_arg .le. -1.d0) then
        angle_abc=dacos(-1.d0)
        dax=0
        day=0
        daz=0
        dbx=0
        dby=0
        dbz=0
        dcx=0
        dcy=0
        dcz=0
        return
      endif
      arccos_arg=aux_1*aux_11_inv
c the actual angle, because it will always be required
      angle_abc=dacos(arccos_arg)
c not sure which is faster
      den_inv=-1.d0/dsqrt(1-arccos_arg**2)
c      den_inv=-1/dabs(dsin(angle_abc))            
c more auxiliary products
      aux_2=2*aux_11_inv
      aux_3=aux_1*aux_31_inv
      aux_3x=-ab_x*aux_3
      aux_3y=-ab_y*aux_3
      aux_3z=-ab_z*aux_3
      aux_4=aux_1*aux_13_inv
      aux_4x=-bc_x*aux_4
      aux_4y=-bc_y*aux_4
      aux_4z=-bc_z*aux_4
c the derivatives
      dax=(bc_x*aux_2-aux_3x)*den_inv
      day=(bc_y*aux_2-aux_3y)*den_inv
      daz=(bc_z*aux_2-aux_3z)*den_inv
      dbx=((-bc_x+ab_x)*aux_2-aux_4x+aux_3x)*den_inv
      dby=((-bc_y+ab_y)*aux_2-aux_4y+aux_3y)*den_inv
      dbz=((-bc_z+ab_z)*aux_2-aux_4z+aux_3z)*den_inv
      dcx=(-ab_x*aux_2+aux_4x)*den_inv
      dcy=(-ab_y*aux_2+aux_4y)*den_inv
      dcz=(-ab_z*aux_2+aux_4z)*den_inv

      return
      END SUBROUTINE DANGLE


c subroutine dangle takes 9 reals (=3 coordinates) and yields all 9 first derivatives of the angle
c via law of cosines (calculating the derivative of Abs[foo] is rather troublesome)
      SUBROUTINE DDANGLE(ax, ay, az, bx, by, bz, cx, cy, cz,
     1 df__dax, df__day, df__daz, df__dbx, df__dby, df__dbz,
     1 df__dcx, df__dcy, df__dcz,
     1 ddf11dax__dax, ddf11dax__day, ddf11dax__daz, ddf11dax__dbx,
     1 ddf11dax__dby, ddf11dax__dbz, ddf11dax__dcx, ddf11dax__dcy,
     1 ddf11dax__dcz, ddf11day__day, ddf11day__daz, ddf11day__dbx,
     1 ddf11day__dby, ddf11day__dbz, ddf11day__dcx, ddf11day__dcy,
     1 ddf11day__dcz, ddf11daz__daz, ddf11daz__dbx, ddf11daz__dby,
     1 ddf11daz__dbz, ddf11daz__dcx, ddf11daz__dcy, ddf11daz__dcz,
     1 ddf11dbx__dbx, ddf11dbx__dby, ddf11dbx__dbz, ddf11dbx__dcx,
     1 ddf11dbx__dcy, ddf11dbx__dcz, ddf11dby__dby, ddf11dby__dbz,
     1 ddf11dby__dcx, ddf11dby__dcy, ddf11dby__dcz, ddf11dbz__dbz,
     1 ddf11dbz__dcx, ddf11dbz__dcy, ddf11dbz__dcz, ddf11dcx__dcx,
     1 ddf11dcx__dcy, ddf11dcx__dcz, ddf11dcy__dcy, ddf11dcy__dcz,
     1 ddf11dcz__dcz,
     1 angle_abc)
      implicit real(8) (a-z)

c vectors from a to b and b to c and a to c
      ab_x=ax-bx
      ab_y=ay-by
      ab_z=az-bz
      bc_x=bx-cx
      bc_y=by-cy
      bc_z=bz-cz
      ac_x=ax-cx
      ac_y=ay-cy
      ac_z=az-cz
c length of a-b and b-c
      r2ab=ab_x**2 + ab_y**2 + ab_z**2
      r2bc=bc_x**2 + bc_y**2 + bc_z**2
      r2ac=ac_x**2 + ac_y**2 + ac_z**2
      r1ab=dsqrt(r2ab)
      r1bc=dsqrt(r2bc)
      r1ac=dsqrt(r2ac)
      r3ab=r2ab*r1ab
      r3bc=r2bc*r1bc
      r3ac=r2ac*r1ac
c some auxiliary products
      aux_11_inv=1/(2*r1ab*r1bc)
      aux_31_inv=1/(2*r3ab*r1bc)
      aux_13_inv=1/(2*r1ab*r3bc)
      aux_1=r2ab + r2bc - r2ac
      arccos_arg=aux_1*aux_11_inv
cc the actual angle, because it will always be required
cc also referred to as 'f'
c if the angle is defined as the smallest angle (which it should be) the derivatives
c feature a singularity at 180 degrees.  Depending on the coordinate system the
c derivative is between 1 and -1 (in case of sides of length 1).  The first derivatives
c are artificially set to zero.  The second redivatives are artificially set to zero
c as well which is far less justified.  Some second derivatives are in fact infinite ...
c we could also gracefully terminate ...
      if(arccos_arg .le. -1.d0 .or. arccos_arg .ge. 1.d0) then
        df__dax=0
        df__day=0
        df__daz=0
        df__dbx=0
        df__dby=0
        df__dbz=0
        df__dcx=0
        df__dcy=0
        df__dcz=0
        ddf11dax__dax=0
        ddf11dax__day=0
        ddf11dax__daz=0
        ddf11dax__dbx=0
        ddf11dax__dby=0
        ddf11dax__dbz=0
        ddf11dax__dcx=0
        ddf11dax__dcy=0
        ddf11dax__dcz=0
        ddf11day__day=0
        ddf11day__daz=0
        ddf11day__dbx=0
        ddf11day__dby=0
        ddf11day__dbz=0
        ddf11day__dcx=0
        ddf11day__dcy=0
        ddf11day__dcz=0
        ddf11daz__daz=0
        ddf11daz__dbx=0
        ddf11daz__dby=0
        ddf11daz__dbz=0
        ddf11daz__dcx=0
        ddf11daz__dcy=0
        ddf11daz__dcz=0
        ddf11dbx__dbx=0
        ddf11dbx__dby=0
        ddf11dbx__dbz=0
        ddf11dbx__dcx=0
        ddf11dbx__dcy=0
        ddf11dbx__dcz=0
        ddf11dby__dby=0
        ddf11dby__dbz=0
        ddf11dby__dcx=0
        ddf11dby__dcy=0
        ddf11dby__dcz=0
        ddf11dbz__dbz=0
        ddf11dbz__dcx=0
        ddf11dbz__dcy=0
        ddf11dbz__dcz=0
        ddf11dcx__dcx=0
        ddf11dcx__dcy=0
        ddf11dcx__dcz=0
        ddf11dcy__dcy=0
        ddf11dcy__dcz=0
        ddf11dcz__dcz=0
        if(arccos_arg .le. -1.d0) then
          angle_abc=dacos(-1.d0)
          return
        elseif(arccos_arg .ge. 1.d0) then
          angle_abc=dacos(1.d0)
          return
        endif
      endif
      angle_abc=dacos(arccos_arg)
      den_inv=-1/dsqrt(1-arccos_arg**2)
c more auxiliary products
      aux_2=2*aux_11_inv
      aux_3=aux_1*aux_31_inv
      aux_3x=ab_x*aux_3
      aux_3y=ab_y*aux_3
      aux_3z=ab_z*aux_3
      aux_4=aux_1*aux_13_inv
      aux_4x=bc_x*aux_4
      aux_4y=bc_y*aux_4
      aux_4z=bc_z*aux_4
c the derivatives
      df__dax=(-bc_x*aux_2 - aux_3x)*den_inv
      df__day=(-bc_y*aux_2 - aux_3y)*den_inv
      df__daz=(-bc_z*aux_2 - aux_3z)*den_inv
      df__dbx=((bc_x-ab_x)*aux_2 - aux_4x + aux_3x)*den_inv
      df__dby=((bc_y-ab_y)*aux_2 - aux_4y + aux_3y)*den_inv
      df__dbz=((bc_z-ab_z)*aux_2 - aux_4z + aux_3z)*den_inv
      df__dcx=(ab_x*aux_2 + aux_4x)*den_inv
      df__dcy=(ab_y*aux_2 + aux_4y)*den_inv
      df__dcz=(ab_z*aux_2 + aux_4z)*den_inv

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SECOND DERIVATIVES

      dr2ab__dab_x=2.d0*ab_x
      dr2ab__dab_y=2.d0*ab_y
      dr2ab__dab_z=2.d0*ab_z

      dr2bc__dbc_x=2.d0*bc_x
      dr2bc__dbc_y=2.d0*bc_y
      dr2bc__dbc_z=2.d0*bc_z

      dr2ac__dac_x=2.d0*ac_x
      dr2ac__dac_y=2.d0*ac_y
      dr2ac__dac_z=2.d0*ac_z

      dr1ab__dr2ab=1.d0/(2.d0*dsqrt(r2ab))

c r1ab=dsqrt(r2ab)
      dr1ab__dax=dr1ab__dr2ab*dr2ab__dab_x
      dr1ab__day=dr1ab__dr2ab*dr2ab__dab_y
      dr1ab__daz=dr1ab__dr2ab*dr2ab__dab_z
      dr1ab__dbx=-dr1ab__dr2ab*dr2ab__dab_x
      dr1ab__dby=-dr1ab__dr2ab*dr2ab__dab_y
      dr1ab__dbz=-dr1ab__dr2ab*dr2ab__dab_z
c      dr1ab__dcx=0
c      dr1ab__dcy=0
c      dr1ab__dcz=0

      dr1bc__dr2bc=1.d0/(2.d0*dsqrt(r2bc))

c r1bc=dsqrt(r2bc)
c      dr1bc__dax=0
c      dr1bc__day=0
c      dr1bc__daz=0
      dr1bc__dbx=dr1bc__dr2bc*dr2bc__dbc_x
      dr1bc__dby=dr1bc__dr2bc*dr2bc__dbc_y
      dr1bc__dbz=dr1bc__dr2bc*dr2bc__dbc_z
      dr1bc__dcx=-dr1bc__dr2bc*dr2bc__dbc_x
      dr1bc__dcy=-dr1bc__dr2bc*dr2bc__dbc_y
      dr1bc__dcz=-dr1bc__dr2bc*dr2bc__dbc_z

      dr1ac__dr2ac=1.d0/(2.d0*dsqrt(r2ac))

c r1ac=dsqrt(r2ac)
      dr1ac__dax=dr1ac__dr2ac*dr2ac__dac_x
      dr1ac__day=dr1ac__dr2ac*dr2ac__dac_y
      dr1ac__daz=dr1ac__dr2ac*dr2ac__dac_z
c      dr1ac__dbx=0
c      dr1ac__dby=0
c      dr1ac__dbz=0
      dr1ac__dcx=-dr1ac__dr2ac*dr2ac__dac_x
      dr1ac__dcy=-dr1ac__dr2ac*dr2ac__dac_y
      dr1ac__dcz=-dr1ac__dr2ac*dr2ac__dac_z

c r3ab=r2ab*r1ab
      dr3ab__dax=r1ab*dr2ab__dab_x + r2ab*dr1ab__dax
      dr3ab__day=r1ab*dr2ab__dab_y + r2ab*dr1ab__day
      dr3ab__daz=r1ab*dr2ab__dab_z + r2ab*dr1ab__daz
      dr3ab__dbx=-r1ab*dr2ab__dab_x + r2ab*dr1ab__dbx
      dr3ab__dby=-r1ab*dr2ab__dab_y + r2ab*dr1ab__dby
      dr3ab__dbz=-r1ab*dr2ab__dab_z + r2ab*dr1ab__dbz
c      dr3ab__dcx=0
c      dr3ab__dcy=0
c      dr3ab__dcz=0

c r3bc=r2bc*r1bc
c      dr3bc__dax=0
c      dr3bc__day=0
c      dr3bc__daz=0
      dr3bc__dbx=r1bc*dr2bc__dbc_x + r2bc*dr1bc__dbx
      dr3bc__dby=r1bc*dr2bc__dbc_y + r2bc*dr1bc__dby
      dr3bc__dbz=r1bc*dr2bc__dbc_z + r2bc*dr1bc__dbz
      dr3bc__dcx=-r1bc*dr2bc__dbc_x + r2bc*dr1bc__dcx
      dr3bc__dcy=-r1bc*dr2bc__dbc_y + r2bc*dr1bc__dcy
      dr3bc__dcz=-r1bc*dr2bc__dbc_z + r2bc*dr1bc__dcz

c r3ac=r2ac*r1ac
      dr3ac__dax=r1ac*dr2ac__dac_x + r2ac*dr1ac__dax
      dr3ac__day=r1ac*dr2ac__dac_y + r2ac*dr1ac__day
      dr3ac__daz=r1ac*dr2ac__dac_z + r2ac*dr1ac__daz
c      dr3ac__dbx=0
c      dr3ac__dby=0
c      dr3ac__dbz=0
      dr3ac__dcx=-r1ac*dr2ac__dac_x + r2ac*dr1ac__dcx
      dr3ac__dcy=-r1ac*dr2ac__dac_y + r2ac*dr1ac__dcy
      dr3ac__dcz=-r1ac*dr2ac__dac_z + r2ac*dr1ac__dcz

      daux_11_inv__dr1ab=-1.d0/(2.d0 * r1ab**2 * r1bc)
      daux_11_inv__dr1bc=-1.d0/(2.d0 * r1ab * r1bc**2)

c aux_11_inv=1/(2*r1ab*r1bc)
      daux_11_inv__dax=daux_11_inv__dr1ab*dr1ab__dax
      daux_11_inv__day=daux_11_inv__dr1ab*dr1ab__day
      daux_11_inv__daz=daux_11_inv__dr1ab*dr1ab__daz
      daux_11_inv__dbx=daux_11_inv__dr1ab*dr1ab__dbx
     1 + daux_11_inv__dr1bc*dr1bc__dbx
      daux_11_inv__dby=daux_11_inv__dr1ab*dr1ab__dby
     1 + daux_11_inv__dr1bc*dr1bc__dby
      daux_11_inv__dbz=daux_11_inv__dr1ab*dr1ab__dbz
     1 + daux_11_inv__dr1bc*dr1bc__dbz
      daux_11_inv__dcx=daux_11_inv__dr1bc*dr1bc__dcx
      daux_11_inv__dcy=daux_11_inv__dr1bc*dr1bc__dcy
      daux_11_inv__dcz=daux_11_inv__dr1bc*dr1bc__dcz

      daux_31_inv__dr3ab=-1.d0/(2.d0 * r3ab**2 * r1bc)
      daux_31_inv__dr1bc=-1.d0/(2.d0 * r3ab * r1bc**2)

c aux_31_inv=1/(2*r3ab*r1bc)
      daux_31_inv__dax=daux_31_inv__dr3ab*dr3ab__dax
      daux_31_inv__day=daux_31_inv__dr3ab*dr3ab__day
      daux_31_inv__daz=daux_31_inv__dr3ab*dr3ab__daz
      daux_31_inv__dbx=daux_31_inv__dr3ab*dr3ab__dbx
     1 + daux_31_inv__dr1bc*dr1bc__dbx
      daux_31_inv__dby=daux_31_inv__dr3ab*dr3ab__dby
     1 + daux_31_inv__dr1bc*dr1bc__dby
      daux_31_inv__dbz=daux_31_inv__dr3ab*dr3ab__dbz
     1 + daux_31_inv__dr1bc*dr1bc__dbz
      daux_31_inv__dcx=daux_31_inv__dr1bc*dr1bc__dcx
      daux_31_inv__dcy=daux_31_inv__dr1bc*dr1bc__dcy
      daux_31_inv__dcz=daux_31_inv__dr1bc*dr1bc__dcz

c aux_13_inv=1/(2*r1ab*r3bc)
      daux_13_inv__dr1ab=-1.d0/(2.d0 * r1ab**2 * r3bc)
      daux_13_inv__dr3bc=-1.d0/(2.d0 * r1ab * r3bc**2)

c aux_13_inv=1/(2*r1ab*r3bc)
      daux_13_inv__dax=daux_13_inv__dr1ab*dr1ab__dax
      daux_13_inv__day=daux_13_inv__dr1ab*dr1ab__day
      daux_13_inv__daz=daux_13_inv__dr1ab*dr1ab__daz
      daux_13_inv__dbx=daux_13_inv__dr1ab*dr1ab__dbx
     1 + daux_13_inv__dr3bc*dr3bc__dbx
      daux_13_inv__dby=daux_13_inv__dr1ab*dr1ab__dby
     1 + daux_13_inv__dr3bc*dr3bc__dby
      daux_13_inv__dbz=daux_13_inv__dr1ab*dr1ab__dbz
     1 + daux_13_inv__dr3bc*dr3bc__dbz
      daux_13_inv__dcx=daux_13_inv__dr3bc*dr3bc__dcx
      daux_13_inv__dcy=daux_13_inv__dr3bc*dr3bc__dcy
      daux_13_inv__dcz=daux_13_inv__dr3bc*dr3bc__dcz

c aux_1=r2ab + r2bc - r2ac
      daux_1__dax=dr2ab__dab_x - dr2ac__dac_x
      daux_1__day=dr2ab__dab_y - dr2ac__dac_y
      daux_1__daz=dr2ab__dab_z - dr2ac__dac_z
      daux_1__dbx=-dr2ab__dab_x + dr2bc__dbc_x 
      daux_1__dby=-dr2ab__dab_y + dr2bc__dbc_y 
      daux_1__dbz=-dr2ab__dab_z + dr2bc__dbc_z 
      daux_1__dcx=-dr2bc__dbc_x + dr2ac__dac_x
      daux_1__dcy=-dr2bc__dbc_y + dr2ac__dac_y
      daux_1__dcz=-dr2bc__dbc_z + dr2ac__dac_z

c arccos_arg=aux_1*aux_11_inv
      darccos_arg__dax=aux_11_inv*daux_1__dax + aux_1*daux_11_inv__dax
      darccos_arg__day=aux_11_inv*daux_1__day + aux_1*daux_11_inv__day
      darccos_arg__daz=aux_11_inv*daux_1__daz + aux_1*daux_11_inv__daz
      darccos_arg__dbx=aux_11_inv*daux_1__dbx + aux_1*daux_11_inv__dbx
      darccos_arg__dby=aux_11_inv*daux_1__dby + aux_1*daux_11_inv__dby
      darccos_arg__dbz=aux_11_inv*daux_1__dbz + aux_1*daux_11_inv__dbz
      darccos_arg__dcx=aux_11_inv*daux_1__dcx + aux_1*daux_11_inv__dcx
      darccos_arg__dcy=aux_11_inv*daux_1__dcy + aux_1*daux_11_inv__dcy
      darccos_arg__dcz=aux_11_inv*daux_1__dcz + aux_1*daux_11_inv__dcz

c den_inv=-1/dsqrt(1-arccos_arg**2)
      dden_inv__darccos_arg=-arccos_arg*((1-arccos_arg**2)**(-1.5d0))

c den_inv=-1/dsqrt(1-arccos_arg**2)
      dden_inv__dax=dden_inv__darccos_arg*darccos_arg__dax
      dden_inv__day=dden_inv__darccos_arg*darccos_arg__day
      dden_inv__daz=dden_inv__darccos_arg*darccos_arg__daz
      dden_inv__dbx=dden_inv__darccos_arg*darccos_arg__dbx
      dden_inv__dby=dden_inv__darccos_arg*darccos_arg__dby
      dden_inv__dbz=dden_inv__darccos_arg*darccos_arg__dbz
      dden_inv__dcx=dden_inv__darccos_arg*darccos_arg__dcx
      dden_inv__dcy=dden_inv__darccos_arg*darccos_arg__dcy
      dden_inv__dcz=dden_inv__darccos_arg*darccos_arg__dcz

c aux_2=2*aux_11_inv
      daux_2__dax=2.d0*daux_11_inv__dax
      daux_2__day=2.d0*daux_11_inv__day
      daux_2__daz=2.d0*daux_11_inv__daz
      daux_2__dbx=2.d0*daux_11_inv__dbx
      daux_2__dby=2.d0*daux_11_inv__dby
      daux_2__dbz=2.d0*daux_11_inv__dbz
      daux_2__dcx=2.d0*daux_11_inv__dcx
      daux_2__dcy=2.d0*daux_11_inv__dcy
      daux_2__dcz=2.d0*daux_11_inv__dcz
      
c aux_3=aux_1*aux_31_inv
      daux_3__dax=aux_31_inv*daux_1__dax + aux_1*daux_31_inv__dax
      daux_3__day=aux_31_inv*daux_1__day + aux_1*daux_31_inv__day
      daux_3__daz=aux_31_inv*daux_1__daz + aux_1*daux_31_inv__daz
      daux_3__dbx=aux_31_inv*daux_1__dbx + aux_1*daux_31_inv__dbx
      daux_3__dby=aux_31_inv*daux_1__dby + aux_1*daux_31_inv__dby
      daux_3__dbz=aux_31_inv*daux_1__dbz + aux_1*daux_31_inv__dbz
      daux_3__dcx=aux_31_inv*daux_1__dcx + aux_1*daux_31_inv__dcx
      daux_3__dcy=aux_31_inv*daux_1__dcy + aux_1*daux_31_inv__dcy
      daux_3__dcz=aux_31_inv*daux_1__dcz + aux_1*daux_31_inv__dcz

c aux_3x=ab_x*aux_3
      daux_3x__dax=aux_3+ab_x*daux_3__dax
      daux_3x__day=ab_x*daux_3__day
      daux_3x__daz=ab_x*daux_3__daz
      daux_3x__dbx=-aux_3+ab_x*daux_3__dbx
      daux_3x__dby=ab_x*daux_3__dby
      daux_3x__dbz=ab_x*daux_3__dbz
      daux_3x__dcx=ab_x*daux_3__dcx
      daux_3x__dcy=ab_x*daux_3__dcy
      daux_3x__dcz=ab_x*daux_3__dcz

c aux_3y=ab_y*aux_3
      daux_3y__day=aux_3+ab_y*daux_3__day
      daux_3y__daz=ab_y*daux_3__daz
      daux_3y__dbx=ab_y*daux_3__dbx
      daux_3y__dby=-aux_3+ab_y*daux_3__dby
      daux_3y__dbz=ab_y*daux_3__dbz
      daux_3y__dcx=ab_y*daux_3__dcx
      daux_3y__dcy=ab_y*daux_3__dcy
      daux_3y__dcz=ab_y*daux_3__dcz

c aux_3z=ab_z*aux_3
      daux_3z__daz=aux_3+ab_z*daux_3__daz
      daux_3z__dbx=ab_z*daux_3__dbx
      daux_3z__dby=ab_z*daux_3__dby
      daux_3z__dbz=-aux_3+ab_z*daux_3__dbz
      daux_3z__dcx=ab_z*daux_3__dcx
      daux_3z__dcy=ab_z*daux_3__dcy
      daux_3z__dcz=ab_z*daux_3__dcz

c aux_4=aux_1*aux_13_inv
      daux_4__dbx=aux_13_inv*daux_1__dbx + aux_1*daux_13_inv__dbx
      daux_4__dby=aux_13_inv*daux_1__dby + aux_1*daux_13_inv__dby
      daux_4__dbz=aux_13_inv*daux_1__dbz + aux_1*daux_13_inv__dbz
      daux_4__dcx=aux_13_inv*daux_1__dcx + aux_1*daux_13_inv__dcx
      daux_4__dcy=aux_13_inv*daux_1__dcy + aux_1*daux_13_inv__dcy
      daux_4__dcz=aux_13_inv*daux_1__dcz + aux_1*daux_13_inv__dcz

c aux_4x=bc_x*aux_4
      daux_4x__dbx=aux_4+bc_x*daux_4__dbx
      daux_4x__dby=bc_x*daux_4__dby
      daux_4x__dbz=bc_x*daux_4__dbz
      daux_4x__dcx=-aux_4+bc_x*daux_4__dcx
      daux_4x__dcy=bc_x*daux_4__dcy
      daux_4x__dcz=bc_x*daux_4__dcz

c aux_4y=bc_y*aux_4
      daux_4y__dby=aux_4+bc_y*daux_4__dby
      daux_4y__dbz=bc_y*daux_4__dbz
      daux_4y__dcx=bc_y*daux_4__dcx
      daux_4y__dcy=-aux_4+bc_y*daux_4__dcy
      daux_4y__dcz=bc_y*daux_4__dcz

c aux_4z=bc_z*aux_4
      daux_4z__dbz=aux_4+bc_z*daux_4__dbz
      daux_4z__dcx=bc_z*daux_4__dcx
      daux_4z__dcy=bc_z*daux_4__dcy
      daux_4z__dcz=-aux_4+bc_z*daux_4__dcz

CCCCC
      ddf11dax__dden_inv=(-bc_x*aux_2 - aux_3x)
      ddf11day__dden_inv=(-bc_y*aux_2 - aux_3y)
      ddf11daz__dden_inv=(-bc_z*aux_2 - aux_3z)
      ddf11dbx__dden_inv=((bc_x-ab_x)*aux_2 - aux_4x + aux_3x)
      ddf11dby__dden_inv=((bc_y-ab_y)*aux_2 - aux_4y + aux_3y)
      ddf11dbz__dden_inv=((bc_z-ab_z)*aux_2 - aux_4z + aux_3z)
      ddf11dcx__dden_inv=(ab_x*aux_2 + aux_4x)
      ddf11dcy__dden_inv=(ab_y*aux_2 + aux_4y)
      ddf11dcz__dden_inv=(ab_z*aux_2 + aux_4z)

      ddf11dax__daux_2=-bc_x*den_inv
      ddf11day__daux_2=-bc_y*den_inv
      ddf11daz__daux_2=-bc_z*den_inv
      ddf11dbx__daux_2=(bc_x-ab_x)*den_inv
      ddf11dby__daux_2=(bc_y-ab_y)*den_inv
      ddf11dbz__daux_2=(bc_z-ab_z)*den_inv
      ddf11dcx__daux_2=ab_x*den_inv
      ddf11dcy__daux_2=ab_y*den_inv
      ddf11dcz__daux_2=ab_z*den_inv

      ddf11dax__dbc_x=-aux_2*den_inv
      ddf11day__dbc_y=-aux_2*den_inv
      ddf11daz__dbc_z=-aux_2*den_inv
      ddf11dbx__dab_x=-aux_2*den_inv
      ddf11dbx__dbc_x=aux_2*den_inv
      ddf11dby__dab_y=-aux_2*den_inv
      ddf11dby__dbc_y=aux_2*den_inv
      ddf11dbz__dab_z=-aux_2*den_inv
      ddf11dbz__dbc_z=aux_2*den_inv

c df__dax=((-bc_x)*aux_2-aux_3x)*den_inv
      ddf11dax__dax=ddf11dax__daux_2*daux_2__dax
     1 - den_inv*daux_3x__dax
     1 + ddf11dax__dden_inv*dden_inv__dax
      ddf11dax__day=ddf11dax__daux_2*daux_2__day
     1 - den_inv*daux_3x__day
     1 + ddf11dax__dden_inv*dden_inv__day
      ddf11dax__daz=ddf11dax__daux_2*daux_2__daz
     1 - den_inv*daux_3x__daz
     1 + ddf11dax__dden_inv*dden_inv__daz
      ddf11dax__dbx=ddf11dax__dbc_x
     1 + ddf11dax__daux_2*daux_2__dbx
     1 - den_inv*daux_3x__dbx
     1 + ddf11dax__dden_inv*dden_inv__dbx
      ddf11dax__dby=ddf11dax__daux_2*daux_2__dby
     1 - den_inv*daux_3x__dby
     1 + ddf11dax__dden_inv*dden_inv__dby
      ddf11dax__dbz=ddf11dax__daux_2*daux_2__dbz
     1 - den_inv*daux_3x__dbz
     1 + ddf11dax__dden_inv*dden_inv__dbz
      ddf11dax__dcx=-ddf11dax__dbc_x
     1 + ddf11dax__daux_2*daux_2__dcx
     1 - den_inv*daux_3x__dcx
     1 + ddf11dax__dden_inv*dden_inv__dcx
      ddf11dax__dcy=ddf11dax__daux_2*daux_2__dcy
     1 - den_inv*daux_3x__dcy
     1 + ddf11dax__dden_inv*dden_inv__dcy
      ddf11dax__dcz=ddf11dax__daux_2*daux_2__dcz
     1 - den_inv*daux_3x__dcz
     1 + ddf11dax__dden_inv*dden_inv__dcz

c df__day=((-bc_y)*aux_2-aux_3y)*den_inv
      ddf11day__day=ddf11day__daux_2*daux_2__day
     1 - den_inv*daux_3y__day
     1 + ddf11day__dden_inv*dden_inv__day
      ddf11day__daz=ddf11day__daux_2*daux_2__daz
     1 - den_inv*daux_3y__daz
     1 + ddf11day__dden_inv*dden_inv__daz
      ddf11day__dbx=ddf11day__daux_2*daux_2__dbx
     1 - den_inv*daux_3y__dbx
     1 + ddf11day__dden_inv*dden_inv__dbx
      ddf11day__dby=ddf11day__dbc_y
     1 + ddf11day__daux_2*daux_2__dby
     1 - den_inv*daux_3y__dby
     1 + ddf11day__dden_inv*dden_inv__dby
      ddf11day__dbz=ddf11day__daux_2*daux_2__dbz
     1 - den_inv*daux_3y__dbz
     1 + ddf11day__dden_inv*dden_inv__dbz
      ddf11day__dcx=ddf11day__daux_2*daux_2__dcx
     1 - den_inv*daux_3y__dcx
     1 + ddf11day__dden_inv*dden_inv__dcx
      ddf11day__dcy=-ddf11day__dbc_y
     1 + ddf11day__daux_2*daux_2__dcy
     1 - den_inv*daux_3y__dcy
     1 + ddf11day__dden_inv*dden_inv__dcy
      ddf11day__dcz=ddf11day__daux_2*daux_2__dcz
     1 - den_inv*daux_3y__dcz
     1 + ddf11day__dden_inv*dden_inv__dcz

c df__daz=((-bc_z)*aux_2-aux_3z)*den_inv
      ddf11daz__daz=ddf11daz__daux_2*daux_2__daz
     1 - den_inv*daux_3z__daz
     1 + ddf11daz__dden_inv*dden_inv__daz
      ddf11daz__dbx=ddf11daz__daux_2*daux_2__dbx
     1 - den_inv*daux_3z__dbx
     1 + ddf11daz__dden_inv*dden_inv__dbx
      ddf11daz__dby=ddf11daz__daux_2*daux_2__dby
     1 - den_inv*daux_3z__dby
     1 + ddf11daz__dden_inv*dden_inv__dby
      ddf11daz__dbz=ddf11daz__dbc_z
     1 + ddf11daz__daux_2*daux_2__dbz
     1 - den_inv*daux_3z__dbz
     1 + ddf11daz__dden_inv*dden_inv__dbz
      ddf11daz__dcx=ddf11daz__daux_2*daux_2__dcx
     1 - den_inv*daux_3z__dcx
     1 + ddf11daz__dden_inv*dden_inv__dcx
      ddf11daz__dcy=ddf11daz__daux_2*daux_2__dcy
     1 - den_inv*daux_3z__dcy
     1 + ddf11daz__dden_inv*dden_inv__dcy
      ddf11daz__dcz=-ddf11daz__dbc_z
     1 + ddf11daz__daux_2*daux_2__dcz
     1 - den_inv*daux_3z__dcz
     1 + ddf11daz__dden_inv*dden_inv__dcz

c df__dbx=((bc_x-ab_x)*aux_2-aux_4x+aux_3x)*den_inv
      ddf11dbx__dbx=ddf11dbx__dbc_x
     1 - ddf11dbx__dab_x
     1 + ddf11dbx__daux_2*daux_2__dbx
     1 + den_inv*daux_3x__dbx
     1 - den_inv*daux_4x__dbx
     1 + ddf11dbx__dden_inv*dden_inv__dbx
      ddf11dbx__dby=ddf11dbx__daux_2*daux_2__dby
     1 + den_inv*daux_3x__dby
     1 - den_inv*daux_4x__dby
     1 + ddf11dbx__dden_inv*dden_inv__dby
      ddf11dbx__dbz=ddf11dbx__daux_2*daux_2__dbz
     1 + den_inv*daux_3x__dbz
     1 - den_inv*daux_4x__dbz
     1 + ddf11dbx__dden_inv*dden_inv__dbz
      ddf11dbx__dcx=-ddf11dbx__dbc_x
     1 + ddf11dbx__daux_2*daux_2__dcx
     1 + den_inv*daux_3x__dcx
     1 - den_inv*daux_4x__dcx
     1 + ddf11dbx__dden_inv*dden_inv__dcx
      ddf11dbx__dcy=ddf11dbx__daux_2*daux_2__dcy
     1 + den_inv*daux_3x__dcy
     1 - den_inv*daux_4x__dcy
     1 + ddf11dbx__dden_inv*dden_inv__dcy
      ddf11dbx__dcz=ddf11dbx__daux_2*daux_2__dcz
     1 + den_inv*daux_3x__dcz
     1 - den_inv*daux_4x__dcz
     1 + ddf11dbx__dden_inv*dden_inv__dcz

c df__dby=((bc_y-ab_y)*aux_2-aux_4y+aux_3y)*den_inv
      ddf11dby__dby=ddf11dby__dbc_y
     1 - ddf11dby__dab_y
     1 + ddf11dby__daux_2*daux_2__dby
     1 + den_inv*daux_3y__dby
     1 - den_inv*daux_4y__dby
     1 + ddf11dby__dden_inv*dden_inv__dby
      ddf11dby__dbz=ddf11dby__daux_2*daux_2__dbz
     1 + den_inv*daux_3y__dbz
     1 - den_inv*daux_4y__dbz
     1 + ddf11dby__dden_inv*dden_inv__dbz
      ddf11dby__dcx=ddf11dby__daux_2*daux_2__dcx
     1 + den_inv*daux_3y__dcx
     1 - den_inv*daux_4y__dcx
     1 + ddf11dby__dden_inv*dden_inv__dcx
      ddf11dby__dcy=-ddf11dby__dbc_y
     1 + ddf11dby__daux_2*daux_2__dcy
     1 + den_inv*daux_3y__dcy
     1 - den_inv*daux_4y__dcy
     1 + ddf11dby__dden_inv*dden_inv__dcy
      ddf11dby__dcz=ddf11dby__daux_2*daux_2__dcz
     1 + den_inv*daux_3y__dcz
     1 - den_inv*daux_4y__dcz
     1 + ddf11dby__dden_inv*dden_inv__dcz

c df__dbz=((bc_z-ab_z)*aux_2-aux_4z+aux_3z)*den_inv
      ddf11dbz__dbz=ddf11dbz__dbc_z
     1 - ddf11dbz__dab_z
     1 + ddf11dbz__daux_2*daux_2__dbz
     1 + den_inv*daux_3z__dbz
     1 - den_inv*daux_4z__dbz
     1 + ddf11dbz__dden_inv*dden_inv__dbz
      ddf11dbz__dcx=ddf11dbz__daux_2*daux_2__dcx
     1 + den_inv*daux_3z__dcx
     1 - den_inv*daux_4z__dcx
     1 + ddf11dbz__dden_inv*dden_inv__dcx
      ddf11dbz__dcy=ddf11dbz__daux_2*daux_2__dcy
     1 + den_inv*daux_3z__dcy
     1 - den_inv*daux_4z__dcy
     1 + ddf11dbz__dden_inv*dden_inv__dcy
      ddf11dbz__dcz=-ddf11dbz__dbc_z
     1 + ddf11dbz__daux_2*daux_2__dcz
     1 + den_inv*daux_3z__dcz
     1 - den_inv*daux_4z__dcz
     1 + ddf11dbz__dden_inv*dden_inv__dcz

c df__dcx=((ab_x)*aux_2+aux_4x)*den_inv
      ddf11dcx__dcx=ddf11dcx__daux_2*daux_2__dcx
     1 + den_inv*daux_4x__dcx
     1 + ddf11dcx__dden_inv*dden_inv__dcx
      ddf11dcx__dcy=ddf11dcx__daux_2*daux_2__dcy
     1 + den_inv*daux_4x__dcy
     1 + ddf11dcx__dden_inv*dden_inv__dcy
      ddf11dcx__dcz=ddf11dcx__daux_2*daux_2__dcz
     1 + den_inv*daux_4x__dcz
     1 + ddf11dcx__dden_inv*dden_inv__dcz

c df__dcy=((ab_y)*aux_2+aux_4y)*den_inv
      ddf11dcy__dcy=ddf11dcy__daux_2*daux_2__dcy
     1 + den_inv*daux_4y__dcy
     1 + ddf11dcy__dden_inv*dden_inv__dcy
      ddf11dcy__dcz=ddf11dcy__daux_2*daux_2__dcz
     1 + den_inv*daux_4y__dcz
     1 + ddf11dcy__dden_inv*dden_inv__dcz

c df__dcz=((ab_z)*aux_2+aux_4z)*den_inv
      ddf11dcz__dcz=ddf11dcz__daux_2*daux_2__dcz
     1 + den_inv*daux_4z__dcz
     1 + ddf11dcz__dden_inv*dden_inv__dcz

      return
      END SUBROUTINE DDANGLE




c subroutine dist takes 12 reals (=4 coordinates) and yields an angel between -\pi and +\pi (in radians)
      SUBROUTINE DIHEDRAL(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2 dihedral_abcd)
      IMPLICIT real(8) (a-z)
c vectors ab, bc and cd
      ab_x=bx-ax
      ab_y=by-ay
      ab_z=bz-az
      bc_x=cx-bx
      bc_y=cy-by
      bc_z=cz-bz
      cd_x=dx-cx
      cd_y=dy-cy
      cd_z=dz-cz
c vector bc normalized to length 1
      bc_length_inv=1.d0/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
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
      abc_length_inv=1.d0/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      bcd_length_inv=1.d0/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
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
      dihedral_abcd=datan2(aux_2, aux_1)
      return
      END SUBROUTINE DIHEDRAL


      SUBROUTINE DDIHEDRAL(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2  dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3  dihedral_abcd)
      implicit real(8) (a-z)
C at first the dihedral (copied from above)
c vectors ab, bc and cd
      ab_x=bx-ax
      ab_y=by-ay
      ab_z=bz-az
      bc_x=cx-bx
      bc_y=cy-by
      bc_z=cz-bz
      cd_x=dx-cx
      cd_y=dy-cy
      cd_z=dz-cz
c vector bc normalized to length 1
      bc_length_inv=1.d0/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
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
      abc_length_inv=1.d0/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
      bcd_length_inv=1.d0/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
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
      dihedral_abcd=datan2(y, x)

C THE DERIVATIVES
c to be read from bottom to top

c derivatives of single vectors
c all other combinations are zero
      dab_x__dax=-1
      dab_x__dbx=1
      dab_y__day=-1
      dab_y__dby=1
      dab_z__daz=-1
      dab_z__dbz=1

      dbc_x__dbx=-1
      dbc_x__dcx=1
      dbc_y__dby=-1
      dbc_y__dcy=1
      dbc_z__dbz=-1
      dbc_z__dcz=1

      dcd_x__dcx=-1
      dcd_x__ddx=1
      dcd_y__dcy=-1
      dcd_y__ddy=1
      dcd_z__dcz=-1
      dcd_z__ddz=1

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

c derivative of y
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

c derivative of x
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
      END SUBROUTINE DDIHEDRAL


c probably not helpful
c      SUBROUTINE COULOMB(ax,ay,az,c)
c      implicit real(8) (a-z)
c      c=1/dsqrt(ax**2 + ay**2 + az**2)
c      return
c      END SUBROUTINE


      SUBROUTINE DCOULOMB(ax,ay,az,dax,day,daz,c)
      implicit real(8) (a-z)
      dist_2=ax**2 + ay**2 + az**2
      c=1.d0/dsqrt(dist_2)
      dist_3_2_inv=c/dist_2
      dax=-ax*dist_3_2_inv
      day=-ay*dist_3_2_inv
      daz=-az*dist_3_2_inv
      return
      END SUBROUTINE DCOULOMB


      SUBROUTINE DDCOULOMB(ax,ay,az,dax,day,daz,
     1 daxax,daxay,daxaz,dayay,dayaz,dazaz,c)
      implicit real(8) (a-z)
      dist_2=ax**2 + ay**2 + az**2
      c=1.d0/dsqrt(dist_2)
      dist_3_2_inv=c/dist_2
      dist_5_2_inv=dist_3_2_inv/dist_2

      dax=-ax*dist_3_2_inv
      day=-ay*dist_3_2_inv
      daz=-az*dist_3_2_inv

      daxax=3*ax**2.0*dist_5_2_inv - dist_3_2_inv
      daxay=3*ax*ay**dist_5_2_inv
      daxaz=3*ax*az**dist_5_2_inv
      dayay=3*ay**2.0*dist_5_2_inv - dist_3_2_inv
      dayaz=3*ay*az**dist_5_2_inv
      dazaz=3*az**2.0*dist_5_2_inv - dist_3_2_inv
      return
      END SUBROUTINE DDCOULOMB


