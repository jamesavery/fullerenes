C to be linked against source/geometry.f and source/dddihedral.f

      PROGRAM testderivatives
      use config
      IMPLICIT REAL*8 (a-z)
      integer iout,testno,minor
      iout=6
      e=1.d-6


C######################################################################
C---DIHEDRAL---
C######################################################################
      write(iout,1000)
      write(*,*)"checking dihedral, ddihedral and dddihedral"
      write(iout,1000)

C######################################################################
      testno=11
      ax=1.1
      ay=0
      az=0.2
      bx=0.1
      by=0
      bz=0.2
      cx=0.1
      cy=1
      cz=0.2
      dx=0.1
      dy=1
      dz=1.2

      minor=1
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(dpi/2)
      end if

      minor=2
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
c      write(*,*)'first derivatives',
c     1        dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(1)).lt.e.and.
     2   dabs(dcx-(1)).lt.e.and.
     2   dabs(dcy-(0)).lt.e.and.
     2   dabs(dcz-(0)).lt.e.and.
     2   dabs(ddx-(-1)).lt.e.and.
     2   dabs(ddy-(0.)).lt.e.and.
     2   dabs(ddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

      minor=3
      if(dabs(angle_abcd - (dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(dpi/2)
      end if
      call dddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     1 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
     1 daxdcz, daxddx, daxddy, daxddz, dayday, daydaz, daydbx, daydby,
     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz, dazdaz,
     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
     1 dazddz, dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
     1 dbxddy, dbxddz, dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
     1 dbyddy, dbyddz, dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
     1 dbzddz, dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz, dcydcy,
     1 dcydcz, dcyddx, dcyddy, dcyddz, dczdcz, dczddx, dczddy, dczddz,
     1 ddxddx, ddxddy, ddxddz, ddyddy, ddyddz, ddzddz,
     1 angle_abcd)
      minor=4
c      write(*,*)'first derivatives',
c     1    dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
c       write(*,*)'2nd, batch 1',
c      1 '(should be 0, 0, -1 -- 0, 0, 1 -- 0, 0, 0 -- 0, 0, 0)',
c      1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
c      1 daxdcz, daxddx, daxddy, daxddz
c       write(*,*)'2nd, batch 2',
c      1 '(should be 0, 0 -- 0, 0, 1 -- 0, 0, -1 -- 0, 0, 0)',
c      1 dayday, daydaz, daydbx, daydby, daydbz, daydcx, daydcy,
c      1 daydcz, dayddx, dayddy, dayddz
c       write(*,*)'2nd, batch 3',
c      1 '(should be 0 -- 1, 0, 0 -- 0, 0, 0 -- 0, 0, 0)',
c      1 dazdaz, dazdbx, dazdby, dazdbz, dazdcx, dazdcy,
c      1 dazdcz, dazddx, dazddy, dazddz
c       write(*,*)'2nd, batch 4',
c      1 '(should be 0, 0, -2 -- 0, -1, 1 -- 0, 1, 0)',
c      1 dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy,
c      1 dbxdcz, dbxddx, dbxddy, dbxddz
c       write(*,*)'2nd, batch 5',
c      1 '(should be  0, -1 -- 0, 0, 0 -- 0, 0, 0)',
c      1 dbydby, dbydbz, dbydcx, dbydcy,
c      1 dbydcz, dbyddx, dbyddy, dbyddz
c       write(*,*)'2nd, batch 6',
c      1 '(should be 0 -- 1, 0, 0 -- 0, 0, 0)',
c      1 dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy, dbzddz
c       write(*,*)'2nd, batch 7',
c      1 '(should be 0, 0, -2 -- 0, -1, 1)',
c      1 dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz
c       write(*,*)'2nd, batch 8',
c      1 '(should be 0, 0 -- 0, 0, 0)',
c      1 dcydcy, dcydcz, dcyddx, dcyddy, dcyddz
c       write(*,*)'2nd, batch 9',
c      1 '(should be 0 -- 1, 0, 0)',
c      1 dczdcz, dczddx, dczddy, dczddz
c       write(*,*)'2nd, batch 10',
c      1 '(should be 0, 0, -1)',
c      1 ddxddx, ddxddy, ddxddz
c       write(*,*)'2nd, batch 11',
c      1 '(should be 0, 0)',
c      1 ddyddy, ddyddz
c       write(*,*)'2nd, batch 12',
c      1 '(should be 0)',
c      1 ddzddz
      if(dabs(daxdax-(0.)).lt.e.and.
     2   dabs(daxday-(0.)).lt.e.and.
     2   dabs(daxdaz-(1.)).lt.e.and.
     2   dabs(daxdbx-(0.)).lt.e.and.
     2   dabs(daxdby-(0.)).lt.e.and.
     2   dabs(daxdbz-(-1.)).lt.e.and.
     2   dabs(daxdcx-(0.)).lt.e.and.
     2   dabs(daxdcy-(0.)).lt.e.and.
     2   dabs(daxdcz-(0.)).lt.e.and.
     2   dabs(daxddx-(0.)).lt.e.and.
     2   dabs(daxddy-(0.)).lt.e.and.
     2   dabs(daxddz-(0.)).lt.e.and.

     2   dabs(dayday-(0.)).lt.e.and.
     2   dabs(daydaz-(0.)).lt.e.and.
     2   dabs(daydbx-(0.)).lt.e.and.
     2   dabs(daydby-(0.)).lt.e.and.
     2   dabs(daydbz-(-1.)).lt.e.and.
     2   dabs(daydcx-(0.)).lt.e.and.
     2   dabs(daydcy-(0.)).lt.e.and.
     2   dabs(daydcz-(1.)).lt.e.and.
     2   dabs(dayddx-(0.)).lt.e.and.
     2   dabs(dayddy-(0.)).lt.e.and.
     2   dabs(dayddz-(0.)).lt.e.and.

     2   dabs(dazdaz-(0.)).lt.e.and.
     2   dabs(dazdbx-(-1.)).lt.e.and.
     2   dabs(dazdby-(0.)).lt.e.and.
     2   dabs(dazdbz-(0.)).lt.e.and.
     2   dabs(dazdcx-(0.)).lt.e.and.
     2   dabs(dazdcy-(0.)).lt.e.and.
     2   dabs(dazdcz-(0.)).lt.e.and.
     2   dabs(dazddx-(0.)).lt.e.and.
     2   dabs(dazddy-(0.)).lt.e.and.
     2   dabs(dazddz-(0.)).lt.e.and.

     2   dabs(dbxdbx-(0.)).lt.e.and.
     2   dabs(dbxdby-(0.)).lt.e.and.
     2   dabs(dbxdbz-(2.)).lt.e.and.
     2   dabs(dbxdcx-(0.)).lt.e.and.
     2   dabs(dbxdcy-(1.)).lt.e.and.
     2   dabs(dbxdcz-(-1.)).lt.e.and.
     2   dabs(dbxddx-(0.)).lt.e.and.
     2   dabs(dbxddy-(-1.)).lt.e.and.
     2   dabs(dbxddz-(0.)).lt.e.and.

     2   dabs(dbydby-(0.)).lt.e.and.
     2   dabs(dbydbz-(1.)).lt.e.and.
     2   dabs(dbydcx-(0.)).lt.e.and.
     2   dabs(dbydcy-(0.)).lt.e.and.
     2   dabs(dbydcz-(-1.)).lt.e.and.
     2   dabs(dbyddx-(0.)).lt.e.and.
     2   dabs(dbyddy-(0.)).lt.e.and.
     2   dabs(dbyddz-(0.)).lt.e.and.

     2   dabs(dbzdbz-(0.)).lt.e.and.
     2   dabs(dbzdcx-(-1.)).lt.e.and.
     2   dabs(dbzdcy-(0.)).lt.e.and.
     2   dabs(dbzdcz-(0.)).lt.e.and.
     2   dabs(dbzddx-(0.)).lt.e.and.
     2   dabs(dbzddy-(0.)).lt.e.and.
     2   dabs(dbzddz-(0.)).lt.e.and.

     2   dabs(dcxdcx-(0.)).lt.e.and.
     2   dabs(dcxdcy-(-1.)).lt.e.and.
     2   dabs(dcxdcz-(2.)).lt.e.and.
     2   dabs(dcxddx-(0.)).lt.e.and.
     2   dabs(dcxddy-(1.)).lt.e.and.
     2   dabs(dcxddz-(-1.)).lt.e.and.

     2   dabs(dcydcy-(0.)).lt.e.and.
     2   dabs(dcydcz-(0.)).lt.e.and.
     2   dabs(dcyddx-(0.)).lt.e.and.
     2   dabs(dcyddy-(0.)).lt.e.and.
     2   dabs(dcyddz-(0.)).lt.e.and.

     2   dabs(dczdcz-(0.)).lt.e.and.
     2   dabs(dczddx-(-1.)).lt.e.and.
     2   dabs(dczddy-(0.)).lt.e.and.
     2   dabs(dczddz-(0.)).lt.e.and.

     2   dabs(ddxddx-(0.)).lt.e.and.
     2   dabs(ddxddy-(0.)).lt.e.and.
     2   dabs(ddxddz-(1.)).lt.e.and.

     2   dabs(ddyddy-(0.)).lt.e.and.
     2   dabs(ddyddz-(0.)).lt.e.and.

     2   dabs(ddzddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=5
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(1)).lt.e.and.
     2   dabs(dcx-(1)).lt.e.and.
     2   dabs(dcy-(0)).lt.e.and.
     2   dabs(dcz-(0)).lt.e.and.
     2   dabs(ddx-(-1)).lt.e.and.
     2   dabs(ddy-(0.)).lt.e.and.
     2   dabs(ddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abcd - (dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(dpi/2)
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=12
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=1
      ay=1
      az=0
      bx=0
      by=0
      bz=0
      cx=-1
      cy=1
      cz=0
      dx=-1
      dy=1
      dz=1
      minor=1
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(dpi/2)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcx-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcy-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcz-(0)).lt.e.and.
     2   dabs(ddx-(-1/sqrt(2.))).lt.e.and.
     2   dabs(ddy-(-1/sqrt(2.))).lt.e.and.
     2   dabs(ddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(dpi/2)
      end if
      call dddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     1 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
     1 daxdcz, daxddx, daxddy, daxddz, dayday, daydaz, daydbx, daydby,
     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz, dazdaz,
     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
     1 dazddz, dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
     1 dbxddy, dbxddz, dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
     1 dbyddy, dbyddz, dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
     1 dbzddz, dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz, dcydcy,
     1 dcydcz, dcyddx, dcyddy, dcyddz, dczdcz, dczddx, dczddy, dczddz,
     1 ddxddx, ddxddy, ddxddz, ddyddy, ddyddz, ddzddz,
     1 angle_abcd)
      minor=4
c      write(*,*)'first derivatives',
c     1    dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
c      write(*,*)'2nd, batch 1',
c     1 '(should be 0, 0, -.35 -- 0, 0, 0 -- 0, 0, 0.35 -- 0, 0, 0)',
c     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
c     1 daxdcz, daxddx, daxddy, daxddz
c      write(*,*)'2nd, batch 2',
c     1 '(should be 0, -.35 -- 0, 0, .707 -- 0, 0, -.35 -- 0, 0, 0)',
c     1 dayday, daydaz, daydbx, daydby, daydbz, daydcx, daydcy,
c     1 daydcz, dayddx, dayddy, dayddz
c      write(*,*)'2nd, batch 3',
c     1 '(should be 0 -- .35, .35, 0 -- 0, 0, 0 -- 0, 0, 0)',
c     1 dazdaz, dazdbx, dazdby, dazdbz, dazdcx, dazdcy,
c     1 dazdcz, dazddx, dazddy, dazddz
c      write(*,*)'2nd, batch 4',
c     1 '(should be 0, 0, -.35 -- .35, -.35, 0 -- -.35, .35, 0)',
c     1 dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy,
c     1 dbxdcz, dbxddx, dbxddy, dbxddz
c      write(*,*)'2nd, batch 5',
c     1 '(should be  0, -1.06 -- .35, -.35, .707 -- -.35, .35, 0)',
c     1 dbydby, dbydbz, dbydcx, dbydcy,
c     1 dbydcz, dbyddx, dbyddy, dbyddz
c      write(*,*)'2nd, batch 6',
c     1 '(should be 0 -- .35, .35, 0 -- 0, 0, 0)',
c     1 dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy, dbzddz
c      write(*,*)'2nd, batch 7',
c     1 '(should be -.707, 0, -1.06 -- .35, -.35, .707)',
c     1 dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz
c      write(*,*)'2nd, batch 8',
c     1 '(should be .707, -1.06 -- .35, -.35, .707)',
c     1 dcydcy, dcydcz, dcyddx, dcyddy, dcyddz
c      write(*,*)'2nd, batch 9',
c     1 '(should be 0 -- .707, .707, 0)',
c     1 dczdcz, dczddx, dczddy, dczddz
c      write(*,*)'2nd, batch 10',
c     1 '(should be 0, 0, -.707)',
c     1 ddxddx, ddxddy, ddxddz
c      write(*,*)'2nd, batch 11',
c     1 '(should be 0, -0.707)',
c     1 ddyddy, ddyddz
c      write(*,*)'2nd, batch 12',
c     1 '(should be 0)',
c     1 ddzddz
      sqrteight=2*sqrt(2.)
      if(dabs(daxdax-(0.)).lt.e.and.
     2   dabs(daxday-(0.)).lt.e.and.
     2   dabs(daxdaz-(1/sqrt(8.))).lt.e.and.
     2   dabs(daxdbx-(0.)).lt.e.and.
     2   dabs(daxdby-(0.)).lt.e.and.
     2   dabs(daxdbz-(0.)).lt.e.and.
     2   dabs(daxdcx-(0.)).lt.e.and.
     2   dabs(daxdcy-(0.)).lt.e.and.
     2   dabs(daxdcz-(-1/sqrt(8.))).lt.e.and.
     2   dabs(daxddx-(0.)).lt.e.and.
     2   dabs(daxddy-(0.)).lt.e.and.
     2   dabs(daxddz-(0.)).lt.e.and.

     2   dabs(dayday-(0.)).lt.e.and.
     2   dabs(daydaz-(1/sqrt(8.))).lt.e.and.
     2   dabs(daydbx-(0.)).lt.e.and.
     2   dabs(daydby-(0.)).lt.e.and.
     2   dabs(daydbz-(-1/sqrt(2.))).lt.e.and.
     2   dabs(daydcx-(0.)).lt.e.and.
     2   dabs(daydcy-(0.)).lt.e.and.
     2   dabs(daydcz-(1/sqrt(8.))).lt.e.and.
     2   dabs(dayddx-(0.)).lt.e.and.
     2   dabs(dayddy-(0.)).lt.e.and.
     2   dabs(dayddz-(0.)).lt.e.and.

     2   dabs(dazdaz-(0.)).lt.e.and.
     2   dabs(dazdbx-(-1/sqrt(8.))).lt.e.and.
     2   dabs(dazdby-(-1/sqrt(8.))).lt.e.and.
     2   dabs(dazdbz-(0.)).lt.e.and.
     2   dabs(dazdcx-(0.)).lt.e.and.
     2   dabs(dazdcy-(0.)).lt.e.and.
     2   dabs(dazdcz-(0.)).lt.e.and.
     2   dabs(dazddx-(0.)).lt.e.and.
     2   dabs(dazddy-(0.)).lt.e.and.
     2   dabs(dazddz-(0.)).lt.e.and.

     2   dabs(dbxdbx-(0.)).lt.e.and.
     2   dabs(dbxdby-(0.)).lt.e.and.
     2   dabs(dbxdbz-(1/sqrteight)).lt.e.and.
     2   dabs(dbxdcx-(-1/sqrteight)).lt.e.and.
     2   dabs(dbxdcy-(1/sqrteight)).lt.e.and.
     2   dabs(dbxdcz-(0.)).lt.e.and.
     2   dabs(dbxddx-(1/sqrteight)).lt.e.and.
     2   dabs(dbxddy-(-1/sqrteight)).lt.e.and.
     2   dabs(dbxddz-(0.)).lt.e.and.

     2   dabs(dbydby-(0.)).lt.e.and.
     2   dabs(dbydbz-(3/sqrteight)).lt.e.and.
     2   dabs(dbydcx-(-1/sqrteight)).lt.e.and.
     2   dabs(dbydcy-(1/sqrteight)).lt.e.and.
     2   dabs(dbydcz-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dbyddx-(1/sqrteight)).lt.e.and.
     2   dabs(dbyddy-(-1/sqrteight)).lt.e.and.
     2   dabs(dbyddz-(0.)).lt.e.and.

     2   dabs(dbzdbz-(0.)).lt.e.and.
     2   dabs(dbzdcx-(-1/sqrteight)).lt.e.and.
     2   dabs(dbzdcy-(-1/sqrteight)).lt.e.and.
     2   dabs(dbzdcz-(0.)).lt.e.and.
     2   dabs(dbzddx-(0.)).lt.e.and.
     2   dabs(dbzddy-(0.)).lt.e.and.
     2   dabs(dbzddz-(0.)).lt.e.and.

     2   dabs(dcxdcx-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcxdcy-(0.)).lt.e.and.
     2   dabs(dcxdcz-(3/sqrteight)).lt.e.and.
     2   dabs(dcxddx-(-1/sqrteight)).lt.e.and.
     2   dabs(dcxddy-(1/sqrteight)).lt.e.and.
     2   dabs(dcxddz-(-1/sqrt(2.))).lt.e.and.

     2   dabs(dcydcy-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dcydcz-(3/sqrteight)).lt.e.and.
     2   dabs(dcyddx-(-1/sqrteight)).lt.e.and.
     2   dabs(dcyddy-(1/sqrteight)).lt.e.and.
     2   dabs(dcyddz-(-1/sqrt(2.))).lt.e.and.

     2   dabs(dczdcz-(0.)).lt.e.and.
     2   dabs(dczddx-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dczddy-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dczddz-(0.)).lt.e.and.

     2   dabs(ddxddx-(0.)).lt.e.and.
     2   dabs(ddxddy-(0.)).lt.e.and.
     2   dabs(ddxddz-(1/sqrt(2.))).lt.e.and.

     2   dabs(ddyddy-(0.)).lt.e.and.
     2   dabs(ddyddz-(1/sqrt(2.))).lt.e.and.

     2   dabs(ddzddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=5
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcx-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcy-(1/sqrt(2.))).lt.e.and.
     2   dabs(dcz-(0)).lt.e.and.
     2   dabs(ddx-(-1/sqrt(2.))).lt.e.and.
     2   dabs(ddy-(-1/sqrt(2.))).lt.e.and.
     2   dabs(ddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abcd - (dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(dpi/2)
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=1
      ay=0
      az=0
      bx=0
      by=1
      bz=0
      cx=-1
      cy=-1
      cz=0
      dx=0
      dy=0
      dz=.2
      minor=1
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (-0.4205343353)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd, (-0.4205343353)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(0.7453559925)).lt.e.and.
     2   dabs(dbx-(-0.4472135955)).lt.e.and.
     2   dabs(dby-(0.2236067977)).lt.e.and.
     2   dabs(dbz-(0.5217491947)).lt.e.and.
     2   dabs(dcx-(-0.298142397)).lt.e.and.
     2   dabs(dcy-(0.1490711985)).lt.e.and.
     2   dabs(dcz-(0.596284794)).lt.e.and.
     2   dabs(ddx-(0.7453559925)).lt.e.and.
     2   dabs(ddy-(-0.3726779962)).lt.e.and.
     2   dabs(ddz-(-1.863389981)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (-0.4205343353)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(-0.4205343353)
      end if
      call dddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     1 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
     1 daxdcz, daxddx, daxddy, daxddz, dayday, daydaz, daydbx, daydby,
     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz, dazdaz,
     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
     1 dazddz, dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
     1 dbxddy, dbxddz, dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
     1 dbyddy, dbyddz, dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
     1 dbzddz, dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz, dcydcy,
     1 dcydcz, dcyddx, dcyddy, dcyddz, dczdcz, dczddx, dczddy, dczddz,
     1 ddxddx, ddxddy, ddxddz, ddyddy, ddyddz, ddzddz,
     1 angle_abcd)
      minor=4
c      write(*,*)'first derivatives',
c     1    dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
c      write(*,*)'2nd, batch 1',
c     1 '(should be 0, 0, 0.496903995 -- 0, 0, -0.2484519975',
c     1 ' -- 0, 0, -0.2484519975 -- 0, 0, 0)',
c     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
c     1 daxdcz, daxddx, daxddy, daxddz
c      write(*,*)'2nd, batch 2',
c     1 '(should be 0, -.35 -- 0, 0, .707 -- 0, 0, -.35 -- 0, 0, 0)',
c     1 dayday, daydaz, daydbx, daydby, daydbz, daydcx, daydcy,
c     1 daydcz, dayddx, dayddy, dayddz
c      write(*,*)'2nd, batch 3',
c     1 '(should be 0 -- .35, .35, 0 -- 0, 0, 0 -- 0, 0, 0)',
c     1 dazdaz, dazdbx, dazdby, dazdbz, dazdcx, dazdcy,
c     1 dazdcz, dazddx, dazddy, dazddz
c      write(*,*)'2nd, batch 4',
c     1 '(should be .77, 0, -.35 -- .35, -.35, 0 -- -.35, .35, 0)',
c     1 dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy,
c     1 dbxdcz, dbxddx, dbxddy, dbxddz
c      write(*,*)'2nd, batch 5',
c     1 '(should be  ,  -- ,   -- , , )',
c     1 dbydby, dbydbz, dbydcx, dbydcy,
c     1 dbydcz, dbyddx, dbyddy, dbyddz
c      write(*,*)'2nd, batch 6',
c     1 '(should be -1.1925 -- -0.4273, 0.5118, -0.6708 --',
c     1 '1.118033989, -1.4907, 1.8633)',
c     1 dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy, dbzddz
c      write(*,*)'2nd, batch 7',
c     1 '(should be 0.5764, -0.1391, -0.4670 -- -1.1428, 0.1987,0.9938)',
c     1 dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz
c      write(*,*)'2nd, batch 8',
c     1 '(should be -0.00496, -0.06459 -- 0.571439, -0.09938, -0.49690)',
c     1 dcydcy, dcydcz, dcyddx, dcyddy, dcyddz
c      write(*,*)'2nd, batch 9',
c     1 '(should be -0.57143959 -- 1.3664859, 0.24845199, 1.2422599)',
c     1 dczdcz, dczddx, dczddy, dczddz
c      write(*,*)'2nd, batch 10',
c     1 '(should be 2.484519975, -1.242259987, -2.484519975)',
c     1 ddxddx, ddxddy, ddxddz
c      write(*,*)'2nd, batch 11',
c     1 '(should be 0.6211299937, 1.242259987)',
c     1 ddyddy, ddyddz
c      write(*,*)'2nd, batch 12',
c     1 '(should be -3.105649969)',
c     1 ddzddz
      if(dabs(daxdax-(0.)).lt.e.and.
     2   dabs(daxday-(0.)).lt.e.and.
     2   dabs(daxdaz-(-0.496903995)).lt.e.and.
     2   dabs(daxdbx-(0.)).lt.e.and.
     2   dabs(daxdby-(0.)).lt.e.and.
     2   dabs(daxdbz-(0.2484519975)).lt.e.and.
     2   dabs(daxdcx-(0.)).lt.e.and.
     2   dabs(daxdcy-(0.)).lt.e.and.
     2   dabs(daxdcz-(0.2484519975)).lt.e.and.
     2   dabs(daxddx-(0.)).lt.e.and.
     2   dabs(daxddy-(0.)).lt.e.and.
     2   dabs(daxddz-(0.)).lt.e.and.

     2   dabs(dayday-(0.)).lt.e.and.
     2   dabs(daydaz-(0.2484519975)).lt.e.and.
     2   dabs(daydbx-(0.)).lt.e.and.
     2   dabs(daydby-(0.)).lt.e.and.
     2   dabs(daydbz-(-0.496903995)).lt.e.and.
     2   dabs(daydcx-(0.)).lt.e.and.
     2   dabs(daydcy-(0.)).lt.e.and.
     2   dabs(daydcz-(0.2484519975)).lt.e.and.
     2   dabs(dayddx-(0.)).lt.e.and.
     2   dabs(dayddy-(0.)).lt.e.and.
     2   dabs(dayddz-(0.)).lt.e.and.

     2   dabs(dazdaz-(0.)).lt.e.and.
     2   dabs(dazdbx-(0.397523196)).lt.e.and.
     2   dabs(dazdby-(-0.198761598)).lt.e.and.
     2   dabs(dazdbz-(0.)).lt.e.and.
     2   dabs(dazdcx-(0.099380799)).lt.e.and.
     2   dabs(dazdcy-(-0.0496903995)).lt.e.and.
     2   dabs(dazdcz-(0.)).lt.e.and.
     2   dabs(dazddx-(0.)).lt.e.and.
     2   dabs(dazddy-(0.)).lt.e.and.
     2   dabs(dazddz-(0.)).lt.e.and.

     2   dabs(dbxdbx-(-0.7751702322)).lt.e.and.
     2   dabs(dbxdby-(0.6111919138)).lt.e.and.
     2   dabs(dbxdbz-(0.4422445555)).lt.e.and.
     2   dabs(dbxdcx-(-0.5664705543)).lt.e.and.
     2   dabs(dbxdcy-(0.4323064756)).lt.e.and.
     2   dabs(dbxdcz-(0.6509442334)).lt.e.and.
     2   dabs(dbxddx-(1.341640786)).lt.e.and.
     2   dabs(dbxddy-(-1.043498389)).lt.e.and.
     2   dabs(dbxddz-(-1.490711985)).lt.e.and.

     2   dabs(dbydby-(-0.4173993558)).lt.e.and.
     2   dabs(dbydbz-(-0.4819968751)).lt.e.and.
     2   dabs(dbydcx-(0.0596284794)).lt.e.and.
     2   dabs(dbydcy-(-0.1043498389)).lt.e.and.
     2   dabs(dbydcz-(-0.06459751935)).lt.e.and.
     2   dabs(dbyddx-(-0.6708203932)).lt.e.and.
     2   dabs(dbyddy-(0.5217491947)).lt.e.and.
     2   dabs(dbyddz-(0.7453559925)).lt.e.and.

     2   dabs(dbzdbz-(1.192569588)).lt.e.and.
     2   dabs(dbzdcx-(0.4273374357)).lt.e.and.
     2   dabs(dbzdcy-(-0.5118111148)).lt.e.and.
     2   dabs(dbzdcz-(0.6708203932)).lt.e.and.
     2   dabs(dbzddx-(-1.118033989)).lt.e.and.
     2   dabs(dbzddy-(1.490711985)).lt.e.and.
     2   dabs(dbzddz-(-1.863389981)).lt.e.and.

     2   dabs(dcxdcx-(-0.5764086342)).lt.e.and.
     2   dabs(dcxdcy-(0.1391331186)).lt.e.and.
     2   dabs(dcxdcz-(0.4670897553)).lt.e.and.
     2   dabs(dcxddx-(1.142879188)).lt.e.and.
     2   dabs(dcxddy-(-0.198761598)).lt.e.and.
     2   dabs(dcxddz-(-0.99380799)).lt.e.and.

     2   dabs(dcydcy-(0.00496903995)).lt.e.and.
     2   dabs(dcydcz-(0.06459751935)).lt.e.and.
     2   dabs(dcyddx-(-0.5714395942)).lt.e.and.
     2   dabs(dcyddy-(0.099380799)).lt.e.and.
     2   dabs(dcyddz-(0.496903995)).lt.e.and.

     2   dabs(dczdcz-(0.5714395942)).lt.e.and.
     2   dabs(dczddx-(-1.366485986)).lt.e.and.
     2   dabs(dczddy-(-0.2484519975)).lt.e.and.
     2   dabs(dczddz-(-1.242259987)).lt.e.and.

     2   dabs(ddxddx-(-2.484519975)).lt.e.and.
     2   dabs(ddxddy-(1.242259987)).lt.e.and.
     2   dabs(ddxddz-(2.484519975)).lt.e.and.

     2   dabs(ddyddy-(-0.6211299937)).lt.e.and.
     2   dabs(ddyddz-(-1.242259987)).lt.e.and.

     2   dabs(ddzddz-(3.105649969)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=5
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(0.7453559925)).lt.e.and.
     2   dabs(dbx-(-0.4472135955)).lt.e.and.
     2   dabs(dby-(0.2236067977)).lt.e.and.
     2   dabs(dbz-(0.5217491947)).lt.e.and.
     2   dabs(dcx-(-0.298142397)).lt.e.and.
     2   dabs(dcy-(0.1490711985)).lt.e.and.
     2   dabs(dcz-(0.596284794)).lt.e.and.
     2   dabs(ddx-(0.7453559925)).lt.e.and.
     2   dabs(ddy-(-0.3726779962)).lt.e.and.
     2   dabs(ddz-(-1.863389981)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abcd - (-0.4205343353)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(-0.4205343353)
      end if



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=14
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=.1
      ay=-.2
      az=.3
      bx=-1.4
      by=.5
      bz=-.6
      cx=.7
      cy=-1.8
      cz=.9
      dx=-.1
      dy=.11
      dz=-1.12
      minor=1
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (1.177791214)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd, (1.177791214)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(-0.6926793076)).lt.e.and.
     2   dabs(day-(0.2444750498)).lt.e.and.
     2   dabs(daz-(1.344612774)).lt.e.and.
     2   dabs(dbx-(0.6121546021)).lt.e.and.
     2   dabs(dby-(0.3479097709)).lt.e.and.
     2   dabs(dbz-(-0.3235547943)).lt.e.and.
     2   dabs(dcx-(0.4397469771)).lt.e.and.
     2   dabs(dcy-(0.02117730752)).lt.e.and.
     2   dabs(dcz-(-0.5831738965)).lt.e.and.
     2   dabs(ddx-(-0.3592222716)).lt.e.and.
     2   dabs(ddy-(-0.6135621281)).lt.e.and.
     2   dabs(ddz-(-0.4378840829)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (1.177791214)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(1.177791214)
      end if
      call dddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     1 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
     1 daxdcz, daxddx, daxddy, daxddz, dayday, daydaz, daydbx, daydby,
     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz, dazdaz,
     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
     1 dazddz, dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
     1 dbxddy, dbxddz, dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
     1 dbyddy, dbyddz, dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
     1 dbzddz, dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz, dcydcy,
     1 dcydcz, dcyddx, dcyddy, dcyddz, dczdcz, dczddx, dczddy, dczddz,
     1 ddxddx, ddxddy, ddxddz, ddyddy, ddyddz, ddzddz,
     1 angle_abcd)
      minor=4
c      write(*,*)'ax',
c     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
c     1 daxdcz, daxddx, daxddy, daxddz
c      write(*,*)'ay',
c     1 dayday, daydaz, daydbx, daydby,
c     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz
c      write(*,*)'az', dazdaz,
c     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
c     1 dazddz
c      write(*,*)'bx',
c     1 dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
c     1 dbxddy, dbxddz
c      write(*,*)'by', dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
c     1 dbyddy, dbyddz
c      write(*,*)'bz',dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
c     1 dbzddz
c      write(*,*)'cx',dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz
c      write(*,*)'cy',dcydcy, dcydcz, dcyddx, dcyddy, dcyddz
c      write(*,*)'cz',dczdcz, dczddx, dczddy, dczddz
c      write(*,*)'dx',ddxddx, ddxddy, ddxddz
c      write(*,*)'dy',ddyddy, ddyddz
c      write(*,*)'dz',ddzddz
      if(dabs(daxdax-(1.38633836)).lt.e.and.
     2   dabs(daxday-(0.529350149)).lt.e.and.
     2   dabs(daxdaz-(-1.129203475)).lt.e.and.
     2   dabs(daxdbx-(-0.7992336876)).lt.e.and.
     2   dabs(daxdby-(-0.2157328256)).lt.e.and.
     2   dabs(daxdbz-(0.7881368301)).lt.e.and.
     2   dabs(daxdcx-(-0.5871046721)).lt.e.and.
     2   dabs(daxdcy-(-0.3136173234)).lt.e.and.
     2   dabs(daxdcz-(0.3410666451)).lt.e.and.
     2   dabs(daxddx-(-0.0)).lt.e.and.
     2   dabs(daxddy-(-0.0)).lt.e.and.
     2   dabs(daxddz-(-0.0)).lt.e.and.

     2   dabs(dayday-(-0.5463515963)).lt.e.and.
     2   dabs(daydaz-(-1.578829323)).lt.e.and.
     2   dabs(daydbx-(-0.1253759383)).lt.e.and.
     2   dabs(daydby-(0.2199498501)).lt.e.and.
     2   dabs(daydbz-(0.5127827503)).lt.e.and.
     2   dabs(daydcx-(-0.4039742107)).lt.e.and.
     2   dabs(daydcy-(0.3264017463)).lt.e.and.
     2   dabs(daydcz-(1.066046573)).lt.e.and.
     2   dabs(dayddx-(-0.)).lt.e.and.
     2   dabs(dayddy-(-0.)).lt.e.and.
     2   dabs(dayddz-(-0.)).lt.e.and.

     2   dabs(dazdaz-(-0.8399867634)).lt.e.and.
     2   dabs(dazdbx-(0.4648978522)).lt.e.and.
     2   dabs(dazdby-(0.802265759)).lt.e.and.
     2   dabs(dazdbz-(0.5792838375)).lt.e.and.
     2   dabs(dazdcx-(0.664305623)).lt.e.and.
     2   dabs(dazdcy-(0.7765635639)).lt.e.and.
     2   dabs(dazdcz-(0.2607029258)).lt.e.and.
     2   dabs(dazddx-(-0.)).lt.e.and.
     2   dabs(dazddy-(-0.)).lt.e.and.
     2   dabs(dazddz-(-0.)).lt.e.and.

     2   dabs(dbxdbx-(0.8274497251)).lt.e.and.
     2   dabs(dbxdby-(0.3235559827)).lt.e.and.
     2   dabs(dbxdbz-(-0.2542073736)).lt.e.and.
     2   dabs(dbxdcx-(0.3400620634)).lt.e.and.
     2   dabs(dbxdcy-(0.02342186139)).lt.e.and.
     2   dabs(dbxdcz-(-0.1470087165)).lt.e.and.
     2   dabs(dbxddx-(-0.3682781009)).lt.e.and.
     2   dabs(dbxddy-(-0.2216019058)).lt.e.and.
     2   dabs(dbxddz-(-0.06368176204)).lt.e.and.

     2   dabs(dbydby-(-0.1876127313)).lt.e.and.
     2   dabs(dbydbz-(-0.5087113831)).lt.e.and.
     2   dabs(dbydcx-(0.2907403967)).lt.e.and.
     2   dabs(dbydcy-(-0.04747758385)).lt.e.and.
     2   dabs(dbydcz-(-0.4657173123)).lt.e.and.
     2   dabs(dbyddx-(-0.3985635538)).lt.e.and.
     2   dabs(dbyddy-(0.01514046505)).lt.e.and.
     2   dabs(dbyddz-(0.1721629363)).lt.e.and.

     2   dabs(dbzdbz-(-0.6398369938)).lt.e.and.
     2   dabs(dbzdcx-(-0.4383880152)).lt.e.and.
     2   dabs(dbzdcy-(-0.3375294151)).lt.e.and.
     2   dabs(dbzdcz-(-0.2925844795)).lt.e.and.
     2   dabs(dbzddx-(-0.09554144126)).lt.e.and.
     2   dabs(dbzddy-(0.3334580478)).lt.e.and.
     2   dabs(dbzddz-(0.3531376358)).lt.e.and.

     2   dabs(dcxdcx-(0.2793529441)).lt.e.and.
     2   dabs(dcxdcy-(0.2733032252)).lt.e.and.
     2   dabs(dcxdcz-(-0.2651938278)).lt.e.and.
     2   dabs(dcxddx-(-0.03231033536)).lt.e.and.
     2   dabs(dcxddy-(-0.1600694112)).lt.e.and.
     2   dabs(dcxddz-(0.03927622002)).lt.e.and.

     2   dabs(dcydcy-(-0.1286343689)).lt.e.and.
     2   dabs(dcydcz-(-0.5939820859)).lt.e.and.
     2   dabs(dcyddx-(0.01689223679)).lt.e.and.
     2   dabs(dcyddy-(-0.1502897936)).lt.e.and.
     2   dabs(dcyddz-(0.1549479371)).lt.e.and.

     2   dabs(dczdcz-(-0.1507185752)).lt.e.and.
     2   dabs(dczddx-(0.07113589924)).lt.e.and.
     2   dabs(dczddy-(-0.0063471744)).lt.e.and.
     2   dabs(dczddz-(0.1826001289)).lt.e.and.

     2   dabs(ddxddx-(0.4005884362)).lt.e.and.
     2   dabs(ddxddy-(0.381671317)).lt.e.and.
     2   dabs(ddxddz-(0.02440554202)).lt.e.and.

     2   dabs(ddyddy-(0.1351493285)).lt.e.and.
     2   dabs(ddyddz-(-0.3271108734)).lt.e.and.

     2   dabs(ddzddz-(-0.5357377647)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=5
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(-0.6926793076)).lt.e.and.
     2   dabs(day-(0.2444750498)).lt.e.and.
     2   dabs(daz-(1.344612774)).lt.e.and.
     2   dabs(dbx-(0.6121546021)).lt.e.and.
     2   dabs(dby-(0.3479097709)).lt.e.and.
     2   dabs(dbz-(-0.3235547943)).lt.e.and.
     2   dabs(dcx-(0.4397469771)).lt.e.and.
     2   dabs(dcy-(0.02117730752)).lt.e.and.
     2   dabs(dcz-(-0.5831738965)).lt.e.and.
     2   dabs(ddx-(-0.3592222716)).lt.e.and.
     2   dabs(ddy-(-0.6135621281)).lt.e.and.
     2   dabs(ddz-(-0.4378840829)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abcd - (1.177791214)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(1.177791214)
      end if


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=15
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=-0.50
      ay=2.16
      az=-1.08
      bx=0.93
      by=-1.76
      bz=-1.63
      cx=2.74
      cy=0.70
      cz=1.12
      dx=2.05
      dy=2.53
      dz=-1.09
      minor=1
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (0.6462092662)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd, (0.6462092662)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(-0.1715882641)).lt.e.and.
     2   dabs(day-(-0.08969841577)).lt.e.and.
     2   dabs(daz-(0.1931755857)).lt.e.and.
     2   dabs(dbx-(0.03283161139)).lt.e.and.
     2   dabs(dby-(0.05458710137)).lt.e.and.
     2   dabs(dbz-(-0.07043981309)).lt.e.and.
     2   dabs(dcx-(-0.1705127382)).lt.e.and.
     2   dabs(dcy-(0.097224564)).lt.e.and.
     2   dabs(dcz-(0.02525659229)).lt.e.and.
     2   dabs(ddx-(0.3092693909)).lt.e.and.
     2   dabs(ddy-(-0.0621132496)).lt.e.and.
     2   dabs(ddz-(-0.1479923649)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (0.6462092662)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(0.6462092662)
      end if
      call dddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     1 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
     1 daxdcz, daxddx, daxddy, daxddz, dayday, daydaz, daydbx, daydby,
     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz, dazdaz,
     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
     1 dazddz, dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
     1 dbxddy, dbxddz, dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
     1 dbyddy, dbyddz, dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
     1 dbzddz, dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz, dcydcy,
     1 dcydcz, dcyddx, dcyddy, dcyddz, dczdcz, dczddx, dczddy, dczddz,
     1 ddxddx, ddxddy, ddxddz, ddyddy, ddyddz, ddzddz,
     1 angle_abcd)
      minor=4
c      write(*,*)'ax',
c     1 daxdax, daxday, daxdaz, daxdbx, daxdby, daxdbz, daxdcx, daxdcy,
c     1 daxdcz, daxddx, daxddy, daxddz
c      write(*,*)'ay',
c     1 dayday, daydaz, daydbx, daydby,
c     1 daydbz, daydcx, daydcy, daydcz, dayddx, dayddy, dayddz
c      write(*,*)'az', dazdaz,
c     1 dazdbx, dazdby, dazdbz, dazdcx, dazdcy, dazdcz, dazddx, dazddy,
c     1 dazddz
c      write(*,*)'bx',
c     1 dbxdbx, dbxdby, dbxdbz, dbxdcx, dbxdcy, dbxdcz, dbxddx,
c     1 dbxddy, dbxddz
c      write(*,*)'by', dbydby, dbydbz, dbydcx, dbydcy, dbydcz, dbyddx,
c     1 dbyddy, dbyddz
c      write(*,*)'bz',dbzdbz, dbzdcx, dbzdcy, dbzdcz, dbzddx, dbzddy,
c     1 dbzddz
c      write(*,*)'cx',dcxdcx, dcxdcy, dcxdcz, dcxddx, dcxddy, dcxddz
c      write(*,*)'cy',dcydcy, dcydcz, dcyddx, dcyddy, dcyddz
c      write(*,*)'cz',dczdcz, dczddx, dczddy, dczddz
c      write(*,*)'dx',ddxddx, ddxddy, ddxddz
c      write(*,*)'dy',ddyddy, ddyddz
c      write(*,*)'dz',ddzddz
      if(dabs(daxdax-(-0.06027909794)).lt.e.and.
     2   dabs(daxday-(0.01854379247)).lt.e.and.
     2   dabs(daxdaz-(0.02308634102)).lt.e.and.
     2   dabs(daxdbx-(0.01131520754)).lt.e.and.
     2   dabs(daxdby-(-0.01874995018)).lt.e.and.
     2   dabs(daxdbz-(0.009325218834)).lt.e.and.
     2   dabs(daxdcx-(0.0489638904)).lt.e.and.
     2   dabs(daxdcy-(0.0002061577102)).lt.e.and.
     2   dabs(daxdcz-(-0.03241155985)).lt.e.and.
     2   dabs(daxddx-(-0.)).lt.e.and.
     2   dabs(daxddy-(-0.)).lt.e.and.
     2   dabs(daxddz-(-0.)).lt.e.and.

     2   dabs(dayday-(0.03586023683)).lt.e.and.
     2   dabs(daydaz-(-0.0442837989)).lt.e.and.
     2   dabs(daydbx-(-0.03412886796)).lt.e.and.
     2   dabs(daydby-(-0.03073472676)).lt.e.and.
     2   dabs(daydbz-(0.04995661049)).lt.e.and.
     2   dabs(daydcx-(0.01558507549)).lt.e.and.
     2   dabs(daydcy-(-0.005125510075)).lt.e.and.
     2   dabs(daydcz-(-0.005672811584)).lt.e.and.
     2   dabs(dayddx-(-0.)).lt.e.and.
     2   dabs(dayddy-(-0.)).lt.e.and.
     2   dabs(dayddz-(-0.)).lt.e.and.

     2   dabs(dazdaz-(0.02441886111)).lt.e.and.
     2   dabs(dazdbx-(-0.03931337256)).lt.e.and.
     2   dabs(dazdby-(0.007216880684)).lt.e.and.
     2   dabs(dazdbz-(0.01941951922)).lt.e.and.
     2   dabs(dazdcx-(0.01622703154)).lt.e.and.
     2   dabs(dazdcy-(0.03706691822)).lt.e.and.
     2   dabs(dazdcz-(-0.04383838033)).lt.e.and.
     2   dabs(dazddx-(-0.)).lt.e.and.
     2   dabs(dazddy-(-0.)).lt.e.and.
     2   dabs(dazddz-(-0.)).lt.e.and.

     2   dabs(dbxdbx-(-0.02354906597)).lt.e.and.
     2   dabs(dbxdby-(-0.001195534247)).lt.e.and.
     2   dabs(dbxdbz-(0.02850779456)).lt.e.and.
     2   dabs(dbxdcx-(-0.01604432683)).lt.e.and.
     2   dabs(dbxdcy-(-0.02428221552)).lt.e.and.
     2   dabs(dbxdcz-(-0.02972300236)).lt.e.and.
     2   dabs(dbxddx-(0.02827818526)).lt.e.and.
     2   dabs(dbxddy-(0.05960661773)).lt.e.and.
     2   dabs(dbxddz-(0.04052858037)).lt.e.and.

     2   dabs(dbydby-(0.04798733298)).lt.e.and.
     2   dabs(dbydbz-(-0.02229011665)).lt.e.and.
     2   dabs(dbydcx-(0.01203912386)).lt.e.and.
     2   dabs(dbydcy-(-0.002552756634)).lt.e.and.
     2   dabs(dbydcz-(0.02971401132)).lt.e.and.
     2   dabs(dbyddx-(0.007906360568)).lt.e.and.
     2   dabs(dbyddy-(-0.01469984959)).lt.e.and.
     2   dabs(dbyddz-(-0.01464077536)).lt.e.and.

     2   dabs(dbzdbz-(-0.02443826701)).lt.e.and.
     2   dabs(dbzdcx-(-0.0121482271)).lt.e.and.
     2   dabs(dbzdcy-(-0.001584185439)).lt.e.and.
     2   dabs(dbzdcz-(0.01859708346)).lt.e.and.
     2   dabs(dbzddx-(-0.0256847863)).lt.e.and.
     2   dabs(dbzddy-(-0.0260823084)).lt.e.and.
     2   dabs(dbzddz-(-0.01357833567)).lt.e.and.

     2   dabs(dcxdcx-(0.02444358697)).lt.e.and.
     2   dabs(dcxdcy-(-0.05509716564)).lt.e.and.
     2   dabs(dcxdcz-(0.09520322666)).lt.e.and.
     2   dabs(dcxddx-(-0.05736315055)).lt.e.and.
     2   dabs(dcxddy-(0.0274729663)).lt.e.and.
     2   dabs(dcxddz-(-0.0992820311)).lt.e.and.

     2   dabs(dcydcy-(0.02678313479)).lt.e.and.
     2   dabs(dcydcz-(-0.02304916573)).lt.e.and.
     2   dabs(dcyddx-(0.07917322346)).lt.e.and.
     2   dabs(dcyddy-(-0.01910486808)).lt.e.and.
     2   dabs(dcyddz-(-0.01243356705)).lt.e.and.

     2   dabs(dczdcz-(-0.05122672176)).lt.e.and.
     2   dabs(dczddx-(-0.03306866444)).lt.e.and.
     2   dabs(dczddy-(-0.0009920340064)).lt.e.and.
     2   dabs(dczddz-(0.07646801863)).lt.e.and.

     2   dabs(ddxddx-(0.02908496529)).lt.e.and.
     2   dabs(ddxddy-(-0.08707958402)).lt.e.and.
     2   dabs(ddxddz-(0.05875345074)).lt.e.and.

     2   dabs(ddyddy-(0.03380471767)).lt.e.and.
     2   dabs(ddyddz-(0.0270743424)).lt.e.and.

     2   dabs(ddzddz-(-0.06288968296)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=5
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(-0.1715882641)).lt.e.and.
     2   dabs(day-(-0.08969841577)).lt.e.and.
     2   dabs(daz-(0.1931755857)).lt.e.and.
     2   dabs(dbx-(0.03283161139)).lt.e.and.
     2   dabs(dby-(0.05458710137)).lt.e.and.
     2   dabs(dbz-(-0.07043981309)).lt.e.and.
     2   dabs(dcx-(-0.1705127382)).lt.e.and.
     2   dabs(dcy-(0.097224564)).lt.e.and.
     2   dabs(dcz-(0.02525659229)).lt.e.and.
     2   dabs(ddx-(0.3092693909)).lt.e.and.
     2   dabs(ddy-(-0.0621132496)).lt.e.and.
     2   dabs(ddz-(-0.1479923649)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abcd - (0.6462092662)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(0.6462092662)
      end if


 1    FORMAT('passed test ',I2,'.',I1)
 2    FORMAT('FAILED test ',I2,'.',I1,' ...')
 3    FORMAT('result is ',D6.5,' but should be ',D6.5)
 1000 FORMAT(72('-'))


      stop
      END PROGRAM testderivatives

