C to be linked against source/geometry.f and source/dddihedral.f

      PROGRAM testderivatives
      use config
      IMPLICIT REAL*8 (a-z)
      integer iout,testno,minor
      iout=6
      e=1.d-6


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


 1    FORMAT('passed test ',I2,'.',I1)
 2    FORMAT('FAILED test ',I2,'.',I1,' ...')
 3    FORMAT('result is ',D6.5,' but should be ',D6.5)
 1000 FORMAT(72('-'))


      stop
      END PROGRAM testderivatives

