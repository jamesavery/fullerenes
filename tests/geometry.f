      PROGRAM testderivations
      use config
      IMPLICIT REAL*8 (a-z)
      integer iout,testno,minor
      iout=6
      e=1.d-6
 
      write(iout,1000)
      write(*,*)"checking dist and ddist"
      write(iout,1000)

      testno=1
      minor=1
      ax=0
      ay=0
      az=0
      bx=1
      by=0
      bz=0
      call dist(ax,ay,az,bx,by,bz,dist_ab)
      if(dabs(dist_ab - 1).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call ddist(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,dist_ab)
      minor=2
      if(dabs(dax-(-1)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(0)).lt.e.and.
     2   dabs(dbx-(1)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(dist_ab - 1).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

      testno=2
      minor=1
      ax=0
      ay=0
      az=0
      bx=0
      by=1
      bz=0
      call dist(ax,ay,az,bx,by,bz,dist_ab)
      if(dabs(dist_ab - 1).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call ddist(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,dist_ab)
      minor=2
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(-1)).lt.e.and.
     2   dabs(daz-(0)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(1)).lt.e.and.
     2   dabs(dbz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(dist_ab - 1).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

      testno=3
      minor=1
      ax=0
      ay=0
      az=0
      bx=0
      by=0
      bz=1
      call dist(ax,ay,az,bx,by,bz,dist_ab)
      if(dabs(dist_ab - 1).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=2
      call ddist(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,dist_ab)
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e .and.
     2   dabs(dbz-(1)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(dist_ab - 1).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

      testno=4
      minor=1
      ax=0
      ay=0
      az=0
      bx=-2
      by=-2
      bz=-2
      call dist(ax,ay,az,bx,by,bz,dist_ab)
      if(dabs(dist_ab - 2*sqrt(3.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call ddist(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,dist_ab)
      minor=2
      if(dabs(dax-(1/sqrt(3.))).lt.e .and.
     2   dabs(day-(1/sqrt(3.))).lt.e .and.
     2   dabs(daz-(1/sqrt(3.))).lt.e.and.
     3   dabs(dbx-(-1/sqrt(3.))).lt.e.and.
     2   dabs(dby-(-1/sqrt(3.))).lt.e.and.
     4   dabs(dbz-(-1/sqrt(3.))).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(dist_ab - 2*sqrt(3.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

      testno=5
      minor=1
      ax=-1
      ay=1
      az=1
      bx=1
      by=1
      bz=-1
      call dist(ax,ay,az,bx,by,bz,dist_ab)
      if(dabs(dist_ab - 2*sqrt(2.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call ddist(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,dist_ab)
      minor=2
      if(dabs(dax-(-1/sqrt(2.))).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(1/sqrt(2.))).lt.e.and.
     2   dabs(dbx-(1/sqrt(2.))).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(-1/sqrt(2.))).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(dist_ab - 2*sqrt(2.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if



       
      write(iout,1000)
      write(*,*)"checking angle and dangle"
      write(iout,1000)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=6
      minor=1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=2
      ay=0
      az=0
      bx=0
      by=0
      bz=0
      cx=0
      cy=1
      cz=0
      call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      if(dabs(angle_abc - dpi/2).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call dangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,angle_abc)
      minor=2
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(-.5)).lt.e.and.
     2   dabs(daz-(0.)).lt.e.and.
     2   dabs(dbx-(1.)).lt.e.and.
     2   dabs(dby-(.5)).lt.e.and.
     2   dabs(dbz-(0.)).lt.e.and.
     2   dabs(dcx-(-1.)).lt.e.and.
     2   dabs(dcy-(0.)).lt.e.and.
     2   dabs(dcz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abc - dpi/2).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=7
      minor=1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=1
      ay=0
      az=0
      bx=0
      by=0
      bz=0
      cx=-1
      cy=1
      cz=0
      call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      if(dabs(angle_abc - 3*dpi/4).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call dangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,angle_abc)
      minor=2
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(-1)).lt.e.and.
     2   dabs(daz-(0)).lt.e.and.
     2   dabs(dbx-(.5)).lt.e.and.
     2   dabs(dby-(1.5)).lt.e.and.
     2   dabs(dbz-(0)).lt.e.and.
     2   dabs(dcx-(-.5)).lt.e.and.
     2   dabs(dcy-(-.5)).lt.e.and.
     2   dabs(dcz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abc - 3*dpi/4).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write(iout,1000)
      write(*,*)"checking dihedral and ddihedral"
      write(iout,1000)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=8
      minor=1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (-dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(-dpi/2)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(1)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(-1)).lt.e.and.
     2   dabs(dcx-(-1)).lt.e.and.
     2   dabs(dcy-(0)).lt.e.and.
     2   dabs(dcz-(0)).lt.e.and.
     2   dabs(ddx-(1)).lt.e.and.
     2   dabs(ddy-(0.)).lt.e.and.
     2   dabs(ddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (-dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(-dpi/2)
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=9
      minor=1
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
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (-dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(-dpi/2)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(1/sqrt(2.))).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dcx-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dcy-(-1/sqrt(2.))).lt.e.and.
     2   dabs(dcz-(0)).lt.e.and.
     2   dabs(ddx-(1/sqrt(2.))).lt.e.and.
     2   dabs(ddy-(1/sqrt(2.))).lt.e.and.
     2   dabs(ddz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (-dpi/2)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(-dpi/2)
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=10
      minor=1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (0.4205343353)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd, (0.4205343353)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
c      write(*,*)dax,' dax ','0'
c      write(*,*)day,' day ','0'
c      write(*,*)daz,' daz ','-0.7453559925'
c      write(*,*)dbx,' dbx ','0.4472135955'
c      write(*,*)dby,' dby ','-0.2236067977'
c      write(*,*)dbz,' dbz ','-0.5217491947'
c      write(*,*)dcx,' dcx ','0.298142397'
c      write(*,*)dcy,' dcy ','-0.1490711985'
c      write(*,*)dcz,' dcz ','-0.596284794'
c      write(*,*)ddx,' ddx ','-0.7453559925'
c      write(*,*)ddy,' ddy ','0.3726779962'
c      write(*,*)ddz,' ddz ','1.863389981'
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-0.7453559925)).lt.e.and.
     2   dabs(dbx-(0.4472135955)).lt.e.and.
     2   dabs(dby-(-0.2236067977)).lt.e.and.
     2   dabs(dbz-(-0.5217491947)).lt.e.and.
     2   dabs(dcx-(0.298142397)).lt.e.and.
     2   dabs(dcy-(-0.1490711985)).lt.e.and.
     2   dabs(dcz-(-0.596284794)).lt.e.and.
     2   dabs(ddx-(-0.7453559925)).lt.e.and.
     2   dabs(ddy-(0.3726779962)).lt.e.and.
     2   dabs(ddz-(1.863389981)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (0.4205343353)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(0.4205343353)
      end if


 1    FORMAT('passed test ',I2,'.',I1)
 2    FORMAT('FAILED test ',I2,'.',I1,' ...')
 3    FORMAT('result is ',D10.8,' but should be ',D10.8)
 1000 FORMAT(70('-'))

      stop
      END
