      PROGRAM testderivations
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


C######################################################################
C---ANGLE---
C######################################################################
      write(iout,1000)
      write(*,*)"checking angle, dangle and ddangle"
      write(iout,1000)

c######################################################################
      testno=6
      minor=1
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
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,angle_abc)
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
      call ddangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     2 daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz,
     2 dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz,
     2 dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz,
     2 dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz,
     2 dbyby,dbybz,dbycx,dbycy,dbycz,
     2 dbzbz,dbzcx,dbzcy,dbzcz,
     2 dcxcx,dcxcy,dcxcz,
     2 dcycy,dcycz,
     2 dczcz,
     3 angle_abc)
      minor=4
      if(dabs(daxax-(0.)).lt.e.and.
     2   dabs(daxay-(0.25)).lt.e.and.
     2   dabs(daxaz-(0.)).lt.e.and.
     2   dabs(daxbx-(0.)).lt.e.and.
     2   dabs(daxby-(-0.25)).lt.e.and.
     2   dabs(daxbz-(0.)).lt.e.and.
     2   dabs(daxcx-(0.)).lt.e.and.
     2   dabs(daxcy-(0.)).lt.e.and.
     2   dabs(daxcz-(0.)).lt.e.and.

     2   dabs(dayay-(0.)).lt.e.and.
     2   dabs(dayaz-(0.)).lt.e.and.
     2   dabs(daybx-(-0.25)).lt.e.and.
     2   dabs(dayby-(0.)).lt.e.and.
     2   dabs(daybz-(0.)).lt.e.and.
     2   dabs(daycx-(0.)).lt.e.and.
     2   dabs(daycy-(0.)).lt.e.and.
     2   dabs(daycz-(0.)).lt.e.and.

     2   dabs(dazaz-(0.)).lt.e.and.
     2   dabs(dazbx-(0.)).lt.e.and.
     2   dabs(dazby-(0.)).lt.e.and.
     2   dabs(dazbz-(0.5)).lt.e.and.
     2   dabs(dazcx-(0.)).lt.e.and.
     2   dabs(dazcy-(0.)).lt.e.and.
     2   dabs(dazcz-(-0.5)).lt.e.and.

     2   dabs(dbxbx-(0.)).lt.e.and.
     2   dabs(dbxby-(1.25)).lt.e.and.
     2   dabs(dbxbz-(0.)).lt.e.and.
     2   dabs(dbxcx-(0.)).lt.e.and.
     2   dabs(dbxcy-(-1.)).lt.e.and.
     2   dabs(dbxcz-(0.)).lt.e.and.

     2   dabs(dbyby-(0.)).lt.e.and.
     2   dabs(dbybz-(0.)).lt.e.and.
     2   dabs(dbycx-(-1.)).lt.e.and.
     2   dabs(dbycy-(0.)).lt.e.and.
     2   dabs(dbycz-(0.)).lt.e.and.

     2   dabs(dbzbz-(-1.)).lt.e.and.
     2   dabs(dbzcx-(0.)).lt.e.and.
     2   dabs(dbzcy-(0.)).lt.e.and.
     2   dabs(dbzcz-(0.5)).lt.e.and.

     2   dabs(dcxcx-(0.)).lt.e.and.
     2   dabs(dcxcy-(1.)).lt.e.and.
     2   dabs(dcxcz-(0.)).lt.e.and.

     2   dabs(dcycy-(0.)).lt.e.and.
     2   dabs(dcycz-(0.)).lt.e.and.

     2   dabs(dczcz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz
        write(*,*)dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz
        write(*,*)dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz
        write(*,*)dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz
        write(*,*)dbyby,dbybz,dbycx,dbycy,dbycz
        write(*,*)dbzbz,dbzcx,dbzcy,dbzcz
        write(*,*)dcxcx,dcxcy,dcxcz
        write(*,*)dcycy,dcycz
        write(*,*)dczcz
      end if
      minor=5
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
      minor=6
      if(dabs(angle_abc - dpi/2).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

C######################################################################
      testno=7
      minor=1
      ax=0
      ay=1
      az=0
      bx=0
      by=0
      bz=0
      cx=0
      cy=0
      cz=1
      call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      if(dabs(angle_abc - dpi/2).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call dangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,angle_abc)
      minor=2
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1.)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(1)).lt.e.and.
     2   dabs(dbz-(1)).lt.e.and.
     2   dabs(dcx-(0)).lt.e.and.
     2   dabs(dcy-(-1)).lt.e.and.
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
      call ddangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     2 daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz,
     2 dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz,
     2 dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz,
     2 dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz,
     2 dbyby,dbybz,dbycx,dbycy,dbycz,
     2 dbzbz,dbzcx,dbzcy,dbzcz,
     2 dcxcx,dcxcy,dcxcz,
     2 dcycy,dcycz,
     2 dczcz,
     3 angle_abc)
      minor=4
      if(dabs(daxax-(0.)).lt.e.and.
     2   dabs(daxay-(0)).lt.e.and.
     2   dabs(daxaz-(0.)).lt.e.and.
     2   dabs(daxbx-(1.)).lt.e.and.
     2   dabs(daxby-(0)).lt.e.and.
     2   dabs(daxbz-(0.)).lt.e.and.
     2   dabs(daxcx-(-1.)).lt.e.and.
     2   dabs(daxcy-(0.)).lt.e.and.
     2   dabs(daxcz-(0.)).lt.e.and.

     2   dabs(dayay-(0.)).lt.e.and.
     2   dabs(dayaz-(1.)).lt.e.and.
     2   dabs(daybx-(0)).lt.e.and.
     2   dabs(dayby-(0.)).lt.e.and.
     2   dabs(daybz-(-1.)).lt.e.and.
     2   dabs(daycx-(0.)).lt.e.and.
     2   dabs(daycy-(0.)).lt.e.and.
     2   dabs(daycz-(0.)).lt.e.and.

     2   dabs(dazaz-(0.)).lt.e.and.
     2   dabs(dazbx-(0.)).lt.e.and.
     2   dabs(dazby-(-1.)).lt.e.and.
     2   dabs(dazbz-(0)).lt.e.and.
     2   dabs(dazcx-(0.)).lt.e.and.
     2   dabs(dazcy-(0.)).lt.e.and.
     2   dabs(dazcz-(0)).lt.e.and.

     2   dabs(dbxbx-(-2.)).lt.e.and.
     2   dabs(dbxby-(0.)).lt.e.and.
     2   dabs(dbxbz-(0.)).lt.e.and.
     2   dabs(dbxcx-(1.)).lt.e.and.
     2   dabs(dbxcy-(0.)).lt.e.and.
     2   dabs(dbxcz-(0.)).lt.e.and.

     2   dabs(dbyby-(0.)).lt.e.and.
     2   dabs(dbybz-(2.)).lt.e.and.
     2   dabs(dbycx-(0.)).lt.e.and.
     2   dabs(dbycy-(0.)).lt.e.and.
     2   dabs(dbycz-(-1.)).lt.e.and.

     2   dabs(dbzbz-(0.)).lt.e.and.
     2   dabs(dbzcx-(0.)).lt.e.and.
     2   dabs(dbzcy-(-1.)).lt.e.and.
     2   dabs(dbzcz-(0)).lt.e.and.

     2   dabs(dcxcx-(0.)).lt.e.and.
     2   dabs(dcxcy-(0.)).lt.e.and.
     2   dabs(dcxcz-(0.)).lt.e.and.

     2   dabs(dcycy-(0.)).lt.e.and.
     2   dabs(dcycz-(1.)).lt.e.and.

     2   dabs(dczcz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz
        write(*,*)dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz
        write(*,*)dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz
        write(*,*)dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz
        write(*,*)dbyby,dbybz,dbycx,dbycy,dbycz
        write(*,*)dbzbz,dbzcx,dbzcy,dbzcz
        write(*,*)dcxcx,dcxcy,dcxcz
        write(*,*)dcycy,dcycz
        write(*,*)dczcz
      end if
      minor=5
      if(dabs(dax-(0.)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(-1.)).lt.e.and.
     2   dabs(dbx-(0.)).lt.e.and.
     2   dabs(dby-(1)).lt.e.and.
     2   dabs(dbz-(1)).lt.e.and.
     2   dabs(dcx-(0.)).lt.e.and.
     2   dabs(dcy-(-1.)).lt.e.and.
     2   dabs(dcz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abc - dpi/2).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

C######################################################################
      testno=8
      minor=1
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
      call ddangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     2 daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz,
     2 dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz,
     2 dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz,
     2 dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz,
     2 dbyby,dbybz,dbycx,dbycy,dbycz,
     2 dbzbz,dbzcx,dbzcy,dbzcz,
     2 dcxcx,dcxcy,dcxcz,
     2 dcycy,dcycz,
     2 dczcz,
     3 angle_abc)
      minor=4
      if(dabs(daxax-(0.)).lt.e.and.
     2   dabs(daxay-(1.)).lt.e.and.
     2   dabs(daxaz-(0.)).lt.e.and.
     2   dabs(daxbx-(0.)).lt.e.and.
     2   dabs(daxby-(-1.)).lt.e.and.
     2   dabs(daxbz-(0.)).lt.e.and.
     2   dabs(daxcx-(0.)).lt.e.and.
     2   dabs(daxcy-(0.)).lt.e.and.
     2   dabs(daxcz-(0.)).lt.e.and.

     2   dabs(dayay-(0.)).lt.e.and.
     2   dabs(dayaz-(0.)).lt.e.and.
     2   dabs(daybx-(-1.)).lt.e.and.
     2   dabs(dayby-(0.)).lt.e.and.
     2   dabs(daybz-(0.)).lt.e.and.
     2   dabs(daycx-(0.)).lt.e.and.
     2   dabs(daycy-(0.)).lt.e.and.
     2   dabs(daycz-(0.)).lt.e.and.

     2   dabs(dazaz-(-1.)).lt.e.and.
     2   dabs(dazbx-(0.)).lt.e.and.
     2   dabs(dazby-(0.)).lt.e.and.
     2   dabs(dazbz-(2.)).lt.e.and.
     2   dabs(dazcx-(0.)).lt.e.and.
     2   dabs(dazcy-(0.)).lt.e.and.
     2   dabs(dazcz-(-1)).lt.e.and.

     2   dabs(dbxbx-
     3      (0.25-sqrt(2.)*(-5./(4.*sqrt(2.))+sqrt(2.)))).lt.e.and.
     2   dabs(dbxby-(1.)).lt.e.and.
     2   dabs(dbxbz-(0.)).lt.e.and.
     2   dabs(dbxcx-(0.5)).lt.e.and.
     2   dabs(dbxcy-(0.)).lt.e.and.
     2   dabs(dbxcz-(0.)).lt.e.and.

     2   dabs(dbyby-
     3      (2.25-sqrt(2.)*(-1./(4.*sqrt(2.))+sqrt(2.)))).lt.e.and.
     2   dabs(dbybz-(0.)).lt.e.and.
     2   dabs(dbycx-(0.)).lt.e.and.
     2   dabs(dbycy-(-.5)).lt.e.and.
     2   dabs(dbycz-(0.)).lt.e.and.

     2   dabs(dbzbz-(-3.5)).lt.e.and.
     2   dabs(dbzcx-(0.)).lt.e.and.
     2   dabs(dbzcy-(0.)).lt.e.and.
     2   dabs(dbzcz-(1.5)).lt.e.and.

     2   dabs(dcxcx-(-.5)).lt.e.and.
     2   dabs(dcxcy-(0.)).lt.e.and.
     2   dabs(dcxcz-(0.)).lt.e.and.

     2   dabs(dcycy-(0.5)).lt.e.and.
     2   dabs(dcycz-(0.)).lt.e.and.

     2   dabs(dczcz-(-.5)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz
        write(*,*)dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz
        write(*,*)dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz
        write(*,*)dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz
        write(*,*)dbyby,dbybz,dbycx,dbycy,dbycz
        write(*,*)dbzbz,dbzcx,dbzcy,dbzcz
        write(*,*)dcxcx,dcxcy,dcxcz
        write(*,*)dcycy,dcycz
        write(*,*)dczcz
      end if
      minor=5
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
      minor=6
      if(dabs(angle_abc - 3*dpi/4).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if


C######################################################################
      testno=9
      minor=1
      ax=1
      ay=0
      az=0
      bx=0
      by=0
      bz=0
      cx=1
      cy=1
      cz=0
      call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      if(dabs(angle_abc - dpi/4).lt.e) then
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
     2   dabs(dby-(.5)).lt.e.and.
     2   dabs(dbz-(0)).lt.e.and.
     2   dabs(dcx-(-.5)).lt.e.and.
     2   dabs(dcy-(.5)).lt.e.and.
     2   dabs(dcz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abc - dpi/4).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call ddangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     2 daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz,
     2 dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz,
     2 dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz,
     2 dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz,
     2 dbyby,dbybz,dbycx,dbycy,dbycz,
     2 dbzbz,dbzcx,dbzcy,dbzcz,
     2 dcxcx,dcxcy,dcxcz,
     2 dcycy,dcycz,
     2 dczcz,
     3 angle_abc)
      minor=4
      if(dabs(daxax-(0.)).lt.e.and.
     2   dabs(daxay-(1.)).lt.e.and.
     2   dabs(daxaz-(0.)).lt.e.and.
     2   dabs(daxbx-(0.)).lt.e.and.
     2   dabs(daxby-(-1.)).lt.e.and.
     2   dabs(daxbz-(0.)).lt.e.and.
     2   dabs(daxcx-(0.)).lt.e.and.
     2   dabs(daxcy-(0.)).lt.e.and.
     2   dabs(daxcz-(0.)).lt.e.and.

     2   dabs(dayay-(0.)).lt.e.and.
     2   dabs(dayaz-(0.)).lt.e.and.
     2   dabs(daybx-(-1.)).lt.e.and.
     2   dabs(dayby-(0.)).lt.e.and.
     2   dabs(daybz-(0.)).lt.e.and.
     2   dabs(daycx-(0.)).lt.e.and.
     2   dabs(daycy-(0.)).lt.e.and.
     2   dabs(daycz-(0.)).lt.e.and.

     2   dabs(dazaz-(1.)).lt.e.and.
     2   dabs(dazbx-(0.)).lt.e.and.
     2   dabs(dazby-(0.)).lt.e.and.
     2   dabs(dazbz-(0.)).lt.e.and.
     2   dabs(dazcx-(0.)).lt.e.and.
     2   dabs(dazcy-(0.)).lt.e.and.
     2   dabs(dazcz-(-1)).lt.e.and.

     2   dabs(dbxbx-(0.5)).lt.e.and.
     2   dabs(dbxby-(1.)).lt.e.and.
     2   dabs(dbxbz-(0.)).lt.e.and.
     2   dabs(dbxcx-(-.5)).lt.e.and.
     2   dabs(dbxcy-(0.)).lt.e.and.
     2   dabs(dbxcz-(0.)).lt.e.and.

     2   dabs(dbyby-(-.5)).lt.e.and.
     2   dabs(dbybz-(0.)).lt.e.and.
     2   dabs(dbycx-(0.)).lt.e.and.
     2   dabs(dbycy-(.5)).lt.e.and.
     2   dabs(dbycz-(0.)).lt.e.and.

     2   dabs(dbzbz-(-.5)).lt.e.and.
     2   dabs(dbzcx-(0.)).lt.e.and.
     2   dabs(dbzcy-(0.)).lt.e.and.
     2   dabs(dbzcz-(.5)).lt.e.and.

     2   dabs(dcxcx-(.5)).lt.e.and.
     2   dabs(dcxcy-(0.)).lt.e.and.
     2   dabs(dcxcz-(0.)).lt.e.and.

     2   dabs(dcycy-(-0.5)).lt.e.and.
     2   dabs(dcycz-(0.)).lt.e.and.

     2   dabs(dczcz-(.5)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz
        write(*,*)dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz
        write(*,*)dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz
        write(*,*)dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz
        write(*,*)dbyby,dbybz,dbycx,dbycy,dbycz
        write(*,*)dbzbz,dbzcx,dbzcy,dbzcz
        write(*,*)dcxcx,dcxcy,dcxcz
        write(*,*)dcycy,dcycz
        write(*,*)dczcz
      end if
      minor=5
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(-1)).lt.e.and.
     2   dabs(daz-(0)).lt.e.and.
     2   dabs(dbx-(.5)).lt.e.and.
     2   dabs(dby-(.5)).lt.e.and.
     2   dabs(dbz-(0)).lt.e.and.
     2   dabs(dcx-(-.5)).lt.e.and.
     2   dabs(dcy-(.5)).lt.e.and.
     2   dabs(dcz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abc - dpi/4).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if

C######################################################################
      testno=10
      minor=1
      ax=2
      ay=0
      az=0
      bx=0
      by=0
      bz=0
      cx=-1
      cy=0
      cz=0
      call angle(ax,ay,az,bx,by,bz,cx,cy,cz,angle_abc)
      if(dabs(angle_abc - dpi).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call dangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,angle_abc)
      minor=2
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(0)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(0)).lt.e.and.
     2   dabs(dcx-(0)).lt.e.and.
     2   dabs(dcy-(0)).lt.e.and.
     2   dabs(dcz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abc - dpi).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      call ddangle(ax,ay,az,bx,by,bz,cx,cy,cz,
     2 dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,
     2 daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz,
     2 dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz,
     2 dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz,
     2 dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz,
     2 dbyby,dbybz,dbycx,dbycy,dbycz,
     2 dbzbz,dbzcx,dbzcy,dbzcz,
     2 dcxcx,dcxcy,dcxcz,
     2 dcycy,dcycz,
     2 dczcz,
     3 angle_abc)
      minor=4
      if(dabs(daxax-(0.)).lt.e.and.
     2   dabs(daxay-(0.)).lt.e.and.
     2   dabs(daxaz-(0.)).lt.e.and.
     2   dabs(daxbx-(0.)).lt.e.and.
     2   dabs(daxby-(0.)).lt.e.and.
     2   dabs(daxbz-(0.)).lt.e.and.
     2   dabs(daxcx-(0.)).lt.e.and.
     2   dabs(daxcy-(0.)).lt.e.and.
     2   dabs(daxcz-(0.)).lt.e.and.

     2   dabs(dayay-(0.)).lt.e.and.
     2   dabs(dayaz-(0.)).lt.e.and.
     2   dabs(daybx-(0.)).lt.e.and.
     2   dabs(dayby-(0.)).lt.e.and.
     2   dabs(daybz-(0.)).lt.e.and.
     2   dabs(daycx-(0.)).lt.e.and.
     2   dabs(daycy-(0.)).lt.e.and.
     2   dabs(daycz-(0.)).lt.e.and.

     2   dabs(dazaz-(0.)).lt.e.and.
     2   dabs(dazbx-(0.)).lt.e.and.
     2   dabs(dazby-(0.)).lt.e.and.
     2   dabs(dazbz-(0.)).lt.e.and.
     2   dabs(dazcx-(0.)).lt.e.and.
     2   dabs(dazcy-(0.)).lt.e.and.
     2   dabs(dazcz-(0)).lt.e.and.

     2   dabs(dbxbx-(0.)).lt.e.and.
     2   dabs(dbxby-(0.)).lt.e.and.
     2   dabs(dbxbz-(0.)).lt.e.and.
     2   dabs(dbxcx-(0.)).lt.e.and.
     2   dabs(dbxcy-(0.)).lt.e.and.
     2   dabs(dbxcz-(0.)).lt.e.and.

     2   dabs(dbyby-(0.)).lt.e.and.
     2   dabs(dbybz-(0.)).lt.e.and.
     2   dabs(dbycx-(0.)).lt.e.and.
     2   dabs(dbycy-(0.)).lt.e.and.
     2   dabs(dbycz-(0.)).lt.e.and.

     2   dabs(dbzbz-(0.)).lt.e.and.
     2   dabs(dbzcx-(0.)).lt.e.and.
     2   dabs(dbzcy-(0.)).lt.e.and.
     2   dabs(dbzcz-(0.)).lt.e.and.

     2   dabs(dcxcx-(0.)).lt.e.and.
     2   dabs(dcxcy-(0.)).lt.e.and.
     2   dabs(dcxcz-(0.)).lt.e.and.

     2   dabs(dcycy-(0.)).lt.e.and.
     2   dabs(dcycz-(0.)).lt.e.and.

     2   dabs(dczcz-(0.)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)daxax,daxay,daxaz,daxbx,daxby,daxbz,daxcx,daxcy,daxcz
        write(*,*)dayay,dayaz,daybx,dayby,daybz,daycx,daycy,daycz
        write(*,*)dazaz,dazbx,dazby,dazbz,dazcx,dazcy,dazcz
        write(*,*)dbxbx,dbxby,dbxbz,dbxcx,dbxcy,dbxcz
        write(*,*)dbyby,dbybz,dbycx,dbycy,dbycz
        write(*,*)dbzbz,dbzcx,dbzcy,dbzcz
        write(*,*)dcxcx,dcxcy,dcxcz
        write(*,*)dcycy,dcycz
        write(*,*)dczcz
      end if
      minor=5
      if(dabs(dax-(0)).lt.e.and.
     2   dabs(day-(0)).lt.e.and.
     2   dabs(daz-(0)).lt.e.and.
     2   dabs(dbx-(0)).lt.e.and.
     2   dabs(dby-(0)).lt.e.and.
     2   dabs(dbz-(0)).lt.e.and.
     2   dabs(dcx-(0)).lt.e.and.
     2   dabs(dcy-(0)).lt.e.and.
     2   dabs(dcz-(0)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=6
      if(dabs(angle_abc - dpi).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(*,*)angle_abc
      end if


C######################################################################
C---DIHEDRAL---
C######################################################################
      write(iout,1000)
      write(*,*)"checking dihedral, ddihedral"
      write(iout,1000)

C######################################################################
      testno=11
      minor=1
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
      testno=12
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
      testno=13
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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testno=14
      minor=1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax=0.123
      ay=0.234
      az=0.348
      bx=0.321
      by=1.432
      bz=0.543
      cx=-1.7
      cy=-1.8
      cz=0.2
      dx=0.3
      dy=0.4
      dz=.25
      call dihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,angle_abcd)
      if(dabs(angle_abcd - (-0.1835399309)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd, (-0.1835399309)
      end if
      call ddihedral(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
     2   dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz,
     3   angle_abcd)
      minor=2
c      write(*,*) dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz,ddx,ddy,ddz
      write(*,*)dax,' dax '
      write(*,*)day,' day '
      write(*,*)daz,' daz '
      write(*,*)dbx,' dbx '
      write(*,*)dby,' dby '
      write(*,*)dbz,' dbz '
      write(*,*)dcx,' dcx '
      write(*,*)dcy,' dcy '
      write(*,*)dcz,' dcz '
      write(*,*)ddx,' ddx '
      write(*,*)ddy,' ddy '
      write(*,*)ddz,' ddz '
      if(dabs(dax-(-0.2522856842)).lt.e.and.
     2   dabs(day-(0.3751985482)).lt.e.and.
     2   dabs(daz-(-2.048898)).lt.e.and.
     2   dabs(dbx-(-0.185615)).lt.e.and.
     2   dabs(dby-(0.0941814)).lt.e.and.
     2   dabs(dbz-(0.2062218)).lt.e.and.
     2   dabs(dcx-(-0.0383647)).lt.e.and.
     2   dabs(dcy-(0.0004204)).lt.e.and.
     2   dabs(dcz-(0.2220884)).lt.e.and.
     2   dabs(ddx-(0.4762657)).lt.e.and.
     2   dabs(ddy-(-.469800)).lt.e.and.
     2   dabs(ddz-(1.620588)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
      end if
      minor=3
      if(dabs(angle_abcd - (-0.1835399309)).lt.e) then
        write(iout,1)testno,minor
      else
        write(iout,2)testno,minor
        write(iout,3)angle_abcd,(0.4205343353)
      end if


 1    FORMAT('passed test ',I2,'.',I1)
 2    FORMAT('FAILED test ',I2,'.',I1,' ...')
 3    FORMAT('result is ',D10.8,' but should be ',D10.8)
 1000 FORMAT(72('-'))


      stop
      END

