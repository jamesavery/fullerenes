C to be linked against source/geometry.f and source/dddihedral.f

      PROGRAM testderivatives
      use config
      IMPLICIT REAL*8 (a-z)
      integer iout,testno,minor
      iout=6
      e=1.d-6


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C---DIST---
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      write(iout,1000)
      write(*,*)"checking dist, ddist and dddist"
      write(iout,1000)

C######################################################################
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

      minor=2
      call ddist(ax,ay,az,bx,by,bz,dax,day,daz,dbx,dby,dbz,dist_ab)
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
      call dddist(ax,ay,az,bx,by,bz,
     2 dax,day,daz,dbx,dby,dbz,
     3 daxax,daxay,daxaz,daxbx,daxby,daxbz,dayay,dayaz,daybx,dayby,
     4 daybz,dazaz,dazbx,dazby,dazbz,dbxbx,dbxby,dbxbz,dbyby,dbybz,
     5 dbzbz,
     6 dist_ab)

      minor=4
      if(dabs(daxax-(0)).lt.e.and.
     2   dabs(daxay-(0)).lt.e.and.
     2   dabs(daxaz-(0)).lt.e.and.
     2   dabs(daxbx-(0)).lt.e.and.
     2   dabs(daxby-(0)).lt.e.and.
     2   dabs(daxbz-(0)).lt.e.and.

     2   dabs(dayay-(1)).lt.e.and.
     2   dabs(dayaz-(0)).lt.e.and.
     2   dabs(daybx-(0)).lt.e.and.
     2   dabs(dayby-(-1)).lt.e.and.
     2   dabs(daybz-(0)).lt.e.and.

     2   dabs(dazaz-(1)).lt.e.and.
     2   dabs(dazbx-(0)).lt.e.and.
     2   dabs(dazby-(0)).lt.e.and.
     2   dabs(dazbz-(-1)).lt.e.and.

     2   dabs(dbxbx-(0)).lt.e.and.
     2   dabs(dbxby-(0)).lt.e.and.
     2   dabs(dbxbz-(0)).lt.e.and.

     2   dabs(dbyby-(1)).lt.e.and.
     2   dabs(dbybz-(0)).lt.e.and.

     2   dabs(dbzbz-(1)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

      minor=5
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

      minor=6
      if(dabs(dist_ab - (1)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

C######################################################################
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
      call dddist(ax,ay,az,bx,by,bz,
     2 dax,day,daz,dbx,dby,dbz,
     3 daxax,daxay,daxaz,daxbx,daxby,daxbz,dayay,dayaz,daybx,dayby,
     4 daybz,dazaz,dazbx,dazby,dazbz,dbxbx,dbxby,dbxbz,dbyby,dbybz,
     5 dbzbz,
     6 dist_ab)
      minor=4
      if(dabs(daxax-(1)).lt.e.and.
     2   dabs(daxay-(0)).lt.e.and.
     2   dabs(daxaz-(0)).lt.e.and.
     2   dabs(daxbx-(-1)).lt.e.and.
     2   dabs(daxby-(0)).lt.e.and.
     2   dabs(daxbz-(0)).lt.e.and.

     2   dabs(dayay-(0)).lt.e.and.
     2   dabs(dayaz-(0)).lt.e.and.
     2   dabs(daybx-(0)).lt.e.and.
     2   dabs(dayby-(0)).lt.e.and.
     2   dabs(daybz-(0)).lt.e.and.

     2   dabs(dazaz-(1)).lt.e.and.
     2   dabs(dazbx-(0)).lt.e.and.
     2   dabs(dazby-(0)).lt.e.and.
     2   dabs(dazbz-(-1)).lt.e.and.

     2   dabs(dbxbx-(1)).lt.e.and.
     2   dabs(dbxby-(0)).lt.e.and.
     2   dabs(dbxbz-(0)).lt.e.and.

     2   dabs(dbyby-(0)).lt.e.and.
     2   dabs(dbybz-(0)).lt.e.and.

     2   dabs(dbzbz-(1)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax, daxay, daxaz, daxbx, daxby, daxbz
        write(*,*)dayay, dayaz, daybx, dayby, daybz
        write(*,*)dazaz, dazbx, dazby, dazbz
        write(*,*)dbxbx, dbxby, dbxbz
        write(*,*)dbyby, dbybz
        write(*,*)dbzbz
      end if
      minor=5
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
      minor=6
      if(dabs(dist_ab - (1)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

C######################################################################
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

C######################################################################
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

C######################################################################
      testno=5
      minor=1
      ax=-1
      ay=1
      az=1
      bx=1
      by=1
      bz=-1
      call dist(ax,ay,az,bx,by,bz,dist_ab)
      if(dabs(dist_ab - (2*sqrt(2.))).lt.e) then
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
      call dddist(ax,ay,az,bx,by,bz,
     2 dax,day,daz,dbx,dby,dbz,
     3 daxax,daxay,daxaz,daxbx,daxby,daxbz,dayay,dayaz,daybx,dayby,
     4 daybz,dazaz,dazbx,dazby,dazbz,dbxbx,dbxby,dbxbz,dbyby,dbybz,
     5 dbzbz,
     6 dist_ab)
      minor=4
      if(dabs(daxax-(1/(4*sqrt(2.)))).lt.e.and.
     2   dabs(daxay-(0)).lt.e.and.
     2   dabs(daxaz-(1/(4*sqrt(2.)))).lt.e.and.
     2   dabs(daxbx-(-1/(4*sqrt(2.)))).lt.e.and.
     2   dabs(daxby-(0)).lt.e.and.
     2   dabs(daxbz-(-1/(4*sqrt(2.)))).lt.e.and.

     2   dabs(dayay-(1/(2*sqrt(2.)))).lt.e.and.
     2   dabs(dayaz-(0)).lt.e.and.
     2   dabs(daybx-(0)).lt.e.and.
     2   dabs(dayby-(-1/(2*sqrt(2.)))).lt.e.and.
     2   dabs(daybz-(0)).lt.e.and.

     2   dabs(dazaz-(1/(4*sqrt(2.)))).lt.e.and.
     2   dabs(dazbx-(-1/(4*sqrt(2.)))).lt.e.and.
     2   dabs(dazby-(0)).lt.e.and.
     2   dabs(dazbz-(-1/(4*sqrt(2.)))).lt.e.and.

     2   dabs(dbxbx-(1/(4*sqrt(2.)))).lt.e.and.
     2   dabs(dbxby-(0)).lt.e.and.
     2   dabs(dbxbz-(1/(4*sqrt(2.)))).lt.e.and.

     2   dabs(dbyby-(1/(2*sqrt(2.)))).lt.e.and.
     2   dabs(dbybz-(0)).lt.e.and.

     2   dabs(dbzbz-(1/(4*sqrt(2.)))).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax, daxay, daxaz, daxbx, daxby, daxbz
        write(*,*)dayay, dayaz, daybx, dayby, daybz
        write(*,*)dazaz, dazbx, dazby, dazbz
        write(*,*)dbxbx, dbxby, dbxbz
        write(*,*)dbyby, dbybz
        write(*,*)dbzbz
      end if
      minor=5
      if(dabs(dax-(-1/sqrt(2.))).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(1/sqrt(2.))).lt.e.and.
     2   dabs(dbx-(1/sqrt(2.))).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(-1/sqrt(2.))).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)dax, day, daz, dbx, dby, dbz
      end if
      minor=6
      if(dabs(dist_ab - (2*sqrt(2.))).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)dist_ab
      end if


 1    FORMAT('passed test ',I2,'.',I1)
 2    FORMAT('FAILED test ',I2,'.',I1,' ...')
 3    FORMAT('result is ',D6.5,' but should be ',D6.5)
 1000 FORMAT(72('-'))


      stop
      END PROGRAM testderivatives

