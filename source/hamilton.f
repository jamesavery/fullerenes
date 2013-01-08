      SUBROUTINE HamiltonCyc(maxiter,Iout,nbatch,A,Nhamilton)
      use config
C     Back-track algorithm from Darko Babic to create Hamitonian cycles
C      optimized for program Isomer
      integer list(Nmax,3,3),path(0:Nmax+1),stack(3*Nmax)
      integer x(0:Nmax)
      integer i,j,k,l,m,last,next,ngb1,ngb2,jlast,jnext,jngb1,jngb2
      integer ptr,prev,oldptr
      logical occ(Nmax),pass(Nmax),end(Nmax),flag
      integer ic3(Nmax,3),A(Nmax,Nmax)
      ifirst=0 
      nhamilton=0
      maxN=30
      nbatch=0

C Prepare field list
      do i=1,number_vertices
        k=0
         do j=1,number_vertices
           if(A(I,J).eq.1) then
            k=k+1
            ic3(i,k)=j
            endif
         end do
      end do
      
      do i=1,number_vertices
         do j=1,3
            list(i,j,1)=ic3(i,j)
         end do
      end do
      
      do i=1,number_vertices
         do j=1,3
            k=list(i,j,1) 
            l=1
            do m=1,3
               if (list(k,m,1).ne.i) then
                  l=l+1
                  list(i,j,l)=m
               endif
            end do
         end do
      end do

c Start algorithm
      do i=1,number_vertices
         pass(i)=.false.
         occ(i)=.false.
         end(i)=.false.
      end do
      do i=0,Nmax+1
         path(i)=0
      end do
      x(0)=number_vertices

      stack(1)=0
      stack(2)=0
      stack(3)=1
      stack(4)=1
      stack(5)=2
      oldptr=3
      ptr=6
      occ(1)=.true.
      last=1
      next=list(1,1,1)
      jnext=1
      end(list(1,2,1))=.true.
      end(list(1,3,1))=.true.
      flag=.false.
      l=1
      path(1)=1
      goto 5

    4 jngb1=list(prev,jlast,2)
      jngb2=list(prev,jlast,3)
      ngb1=list(last,jngb1,1)
      ngb2=list(last,jngb2,1)

      if (occ(ngb1)) then
         next=ngb2
         jnext=jngb2
      else
         if (occ(ngb2)) then
            next=ngb1
            jnext=jngb1
         else
            if (pass(ngb1)) then
               if (pass(ngb2).or.end(ngb1)) go to 6
               next=ngb1
               jnext=jngb1
               pass(ngb2)=.true.
               stack(ptr)=ngb2
               ptr=ptr+1
            else
               if (pass(ngb2)) then
                  if (end(ngb2)) go to 6
                  next=ngb2
                  jnext=jngb2
                  pass(ngb1)=.true.
                  stack(ptr)=ngb1
                  ptr=ptr+1
               else
                  next=ngb1
                  jnext=jngb1
                  pass(ngb2)=.true.
                  stack(ptr)=oldptr
                  stack(ptr+1)=l
                  stack(ptr+2)=jngb2
                  oldptr=ptr
                  ptr=ptr+3
               endif
            endif
         endif
      endif

    5 path(l+1)=next
      if (l.eq.number_vertices-1) then
      nhamilton=nhamilton+1
      if(nhamilton.gt.maxiter) then
       if(maxiter.lt.1000000000) then
        nbatch=nbatch+1
        Write(Iout,1000) nbatch,nhamilton
        nhamilton=0
       else 
        Return
       endif
      endif
      endif

      if (end(next)) then
         if (flag.or.pass(next)) go to 6
         flag=.true.
      endif

      l=l+1
      occ(next)=.true.

      prev=last
      last=next
      jlast=jnext
      go to 4


    6 l1=stack(oldptr+1)
      do i=oldptr+3,ptr-1
         pass(stack(i))=.false.
      end do
      do i=l1+1,l
         occ(path(i))=.false.
         if (end(path(i))) flag=.false.
      end do
      ptr=oldptr
      oldptr=stack(ptr)
C     put this one in to avoid segmentation fault
      if(oldptr.le.0) return
      pass(path(l1+1))=.true.
      stack(ptr)=path(l1+1)
      ptr=ptr+1
      last=path(l1)
      jnext=stack(ptr+1)
      next=list(last,jnext,1)
      pass(next)=.false.

      l=l1
      go to 5
C     if (oldptr.gt.0) go to 5

1000  Format(1X,'Batch ',I3,' of Hamiltonian cycles: ',I10)
      return
      END

      SUBROUTINE Hamilton(Iout,iprint,ihamstore,maxiter,IC3,filename)
C     Back-track algorithm from Darko Babic to create Hamitonian cycles
C     and the IUPAC name of a fullerene
      use config
      integer list(Nmax,3,3),path(0:Nmax+1),stack(3*Nmax),pos(Nmax)
      integer bridge(Nmax),x(0:Nmax),y(0:Nmax),saved(Nmax)
      integer i,j,k,l,m,last,next,ngb1,ngb2,jlast,jnext,jngb1,jngb2
      integer ptr,prev,oldptr,cur,prv,nxt,ngb,diff,maxdif,relk,relbr
      logical occ(Nmax),pass(Nmax),end(Nmax),flag,better
      integer ic3(Nmax,3)
      CHARACTER*50 filename,hamname
      
      if(ihamstore.ne.0) then
       hamname=trim(filename)//".ham"
       Open(unit=8,file=hamname,form='formatted')
       Write(8,*) number_vertices
      endif

      ifirst=0 
      nhamilton=0
      maxN=30
      if(number_vertices.lt.maxN) maxN=number_vertices
      write (Iout,1009) maxiter

C Set field list
      do i=1,number_vertices
         do j=1,3
            list(i,j,1)=ic3(i,j)
         end do
      end do

      do i=1,number_vertices
         do j=1,3
            k=list(i,j,1) 
            l=1
            do m=1,3
               if (list(k,m,1).ne.i) then
                  l=l+1
                  list(i,j,l)=m
               endif
            end do
         end do
      end do

C Start algorithm
      do i=1,number_vertices
         pass(i)=.false.
         occ(i)=.false.
         end(i)=.false.
      end do
      do i=0,Nmax+1
         path(i)=0
      end do
      x(0)=number_vertices

      stack(1)=0
      stack(2)=0
      stack(3)=1
      stack(4)=1
      stack(5)=2
      oldptr=3
      ptr=6
      occ(1)=.true.
      last=1
      next=list(1,1,1)
      jnext=1
      end(list(1,2,1))=.true.
      end(list(1,3,1))=.true.
      flag=.false.
      l=1
      path(1)=1
      goto 5

    4 jngb1=list(prev,jlast,2)
      jngb2=list(prev,jlast,3)
      ngb1=list(last,jngb1,1)
      ngb2=list(last,jngb2,1)

      if (occ(ngb1)) then
         next=ngb2
         jnext=jngb2
      else
         if (occ(ngb2)) then
            next=ngb1
            jnext=jngb1
         else
            if (pass(ngb1)) then
               if (pass(ngb2).or.end(ngb1)) go to 6
               next=ngb1
               jnext=jngb1
               pass(ngb2)=.true.
               stack(ptr)=ngb2
               ptr=ptr+1
            else
               if (pass(ngb2)) then
                  if (end(ngb2)) go to 6
                  next=ngb2
                  jnext=jngb2
                  pass(ngb1)=.true.
                  stack(ptr)=ngb1
                  ptr=ptr+1
               else
                  next=ngb1
                  jnext=jngb1
                  pass(ngb2)=.true.
                  stack(ptr)=oldptr
                  stack(ptr+1)=l
                  stack(ptr+2)=jngb2
                  oldptr=ptr
                  ptr=ptr+3
               endif
            endif
         endif
      endif

    5 path(l+1)=next
      if (l.eq.number_vertices-1) then
      nhamilton=nhamilton+1
      if(nhamilton.gt.maxiter) then
      write (Iout,1010) maxiter
      Return
      endif
      if(ifirst.eq.0) then
      write (Iout,1011)
      endif
      if(iprint.ne.0) then
      if(ifirst.eq.0) write (Iout,1005)
      write (Iout,1004) nhamilton,(path(j),j=1,maxN)
      if(ihamstore.ne.0) write (8,1020) (path(j),j=1,maxN)
      if(number_vertices.gt.30) then
      do I=31,number_vertices,30
      jmax=I+29
      if(jmax.gt.number_vertices) jmax=number_vertices
      write (Iout,1001) (path(j),j=I,jmax)
      enddo
      endif
      endif
      ifirst=1
         do j=1,number_vertices
            pos(path(j))=j
         end do
         path(0)=path(number_vertices)
         path(number_vertices+1)=path(1)

         maxdif=number_vertices
         do j=1,number_vertices
            cur=path(j)
            prv=path(j-1)
            nxt=path(j+1)
            do k=1,3
               ngb=list(cur,k,1)
               if (ngb.ne.prv.and.ngb.ne.nxt) then
                  bridge(j)=pos(ngb)
                  diff=abs(bridge(j)-j)
                  if (number_vertices-diff.gt.diff)
     1              diff=number_vertices-diff
                  if (maxdif.gt.diff) maxdif=diff
                  go to 11
               endif
            end do
   11    continue
         end do

         maxdif=maxdif-1
         if (maxdif.gt.x(0)) go to 6

         do 13 j=1,number_vertices
            better=.false.
            diff=mod(number_vertices+bridge(j)-j,number_vertices)-1
            if (diff.eq.maxdif) then
               if (maxdif.lt.x(0)) then
                  x(0)=maxdif
                  y(0)=number_vertices-maxdif-2
                  better=.true.
               endif
               i=0
               k=j
               do 14 m=1,number_vertices-1
                  k=mod(k,number_vertices)+1
                  relk=mod(number_vertices+k-j,number_vertices)+1
                  relbr=
     1                mod(number_vertices+bridge(k)-j,number_vertices)+1
                  if (relbr.lt.relk) go to 14
                  i=i+1
                  if (.not.better) then
                     if (x(i)-relk) 17,15,16
   15                if (y(i)-relbr) 17,14,16
   16                better=.true.
                  endif
                  x(i)=relk
                  y(i)=relbr
   14          continue
            endif

            if (better) then
               do m=1,number_vertices
                  saved(m)=path(m)
               end do
            end if
   17       better=.false.
            diff=number_vertices-diff-2
            if (diff.eq.maxdif) then
               if (maxdif.lt.x(0)) then
                  x(0)=maxdif
                  y(0)=number_vertices-maxdif-2
                  better=.true.
               endif
               i=0
               k=j
               do 18 m=1,number_vertices-1
                  k=mod(number_vertices+k-2,number_vertices)+1
                  relk=mod(number_vertices+j-k,number_vertices)+1
                  relbr=
     1               mod(number_vertices+j-bridge(k),number_vertices)+1
                  if (relbr.lt.relk) go to 18
                  i=i+1
                  if (.not.better) then
                     if (x(i)-relk) 13,19,20
   19                if (y(i)-relbr) 13,18,20
   20                better=.true.
                  endif
                  x(i)=relk
                  y(i)=relbr
   18          continue
            endif
            if (better) then
               do m=1,number_vertices
                  saved(m)=path(m)
               end do
            end if
   13    continue
         go to 6
      endif
      if (end(next)) then
         if (flag.or.pass(next)) go to 6
         flag=.true.
      endif

      l=l+1
      occ(next)=.true.

      prev=last
      last=next
      jlast=jnext
      go to 4


    6 l1=stack(oldptr+1)
      do i=oldptr+3,ptr-1
         pass(stack(i))=.false.
      end do
      do i=l1+1,l
         occ(path(i))=.false.
         if (end(path(i))) flag=.false.
      end do
      ptr=oldptr
      oldptr=stack(ptr)
C     put this one in to avoid segmentation fault
      if (oldptr.le.0) return
      pass(path(l1+1))=.true.
      stack(ptr)=path(l1+1)
      ptr=ptr+1
      last=path(l1)
      jnext=stack(ptr+1)
      next=list(last,jnext,1)
      pass(next)=.false.

      l=l1
C     if (oldptr.gt.0) go to 5
      go to 5

      if (x(0).eq.0) then
         write (Iout,1002) 
         Return
      endif
 
      write (Iout,1003) nhamilton
      write (Iout,1000) (x(i),i=0,(number_vertices-2)/2)
      write (Iout,1000) (y(i),i=0,(number_vertices-2)/2)
      write (Iout,1006)
      write (Iout,1008) (saved(j),j=1,number_vertices)
 1000 format (1X,i3,i5,98I3)
 1001 format (9X,30(I4,'-'))
 1002 format (/1X,'There is no Hamiltonian cycle in the graph!')
 1003 format (/1X,I10,' Hamiltonian cycles found',
     1 //1X,'Half-ring sizes and IUPAC superscripts:')
 1004 format (1X,I10,': ',30(I4,'-'))
 1005 format (/1X,'Distinct Hamiltonian cycles:',/6X,'NH    vertices',
     1 /1X,20('-'))
 1006 format (/1X,'The best Hamiltonian cycle:')
 1008 format (30(I4,'-'))
 1009 format (/1X,'Calculate Hamiltonian cycles, half-ring sizes'
     1 ' and IUPAC superscripts: (D. Babic, J. Chem. Inf. Comput. Sci.'
     1 ' 35, 515-526 (1995).)',/1X,'Maximum allowed iteration: ',I10)
 1010 format (I10,' Maximum Hamiltonian cycles reached: Return')
 1011 format (1X,' Hamiltonian cycle detected')
 1020 format (500I3)
      return
      END

      SUBROUTINE PathStatistic(IOUT,iprintf,IA,A,evec,df)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C Calculate the number of all possible distinct paths from the
C adjaceny matric A(i,j) by producing the (n-1)-th power
C of A. In this case n is number_vertices. This gives all possible walks
C of length (n-1).
      DIMENSION A(Nmax,Nmax),evec(Nmax),df(Nmax)
      DIMENSION IA(Nmax,Nmax),IM(Nmax,Nmax)
      DIMENSION IMF(Nmax,Nmax),IMF1(Nmax,Nmax)
      DIMENSION IS1(10),IS2(10),APN(10)
      Integer IHamCycmax(42),IHamCycmin(42)
      Integer IHAMIPRmin(32),IHAMIPRmax(32)
      Data Ihuge,over,explimit/180,1.d-10,45.d0/
      Data IHamCycmin/30,0,34,24,18,20,40,28,42,24,68,44,120,76,152,80,
     1 262,66,440,173,618,288,1062,197,1750,320,2688,1182,4230,1596,
     1 7110,2400,10814,1980,17905,1280,29944,7930,46231,13307,72168,
     1 20754/
      Data IHamCycmax/30,0,34,24,43,32,76,66,128,96,280,150,327,260,512,
     1 410,806,642,1746,1068,3040,1802,3340,3096,6018,4818,10428,7832,
     1 15926,12226,35200,20856,39067,33427,76063,51586,117106,90221,
     1 209692,156288,417280,249148/
      Data IHAMIPRmin/ 1090,0,0,0,0,2790,3852,4794,6078,6988,9004,11226,
     1 14748,17853,22661,29277,36949,44730,60070,71950,93986,35907,
     1 149920,180243,237580,244254,383218,457235,630059,723505,1038971,
     1 1368498/
      Data IHAMIPRmax/1090,0,0,0,0,2790,3852,4794,6643,8244,10970,13614,
     1 18260,21756,28652,36852,47054,59118,78044,95694,131690,161148,
     1 207165,257746,351976,426750,571622,699908,1013844,1151918,
     1 1590875,1888558/

      dAtom=dfloat(number_vertices)
C General fullerenes
C     Correct upper and lower limit
      if(number_vertices.le.102) then
       ifield=number_vertices/2-9
       Write(Iout,1016) IHamCycmin(ifield),IHamCycmax(ifield)

      else

C     Conjectured upper and lower limit obtained from D5H and D5d nanotubes
C     Lower limt
       alow=1.d-1*dAtom+1.33d0
       exp1=dAtom/1.d1-1.d0
       if(alow.le.explimit) then
        alowerNT=5.d0*2.d0**exp1
        write (Iout,1012) dint(alowerNT+over)
       else
        write (Iout,1013) exp1
       endif

C     Upper limt
       ahigh=1.8d-1*dAtom+2.33d0
       exp1=dAtom/1.d1-1.d0
       exp2=dAtom/2.d1-1.d0
       if(ahigh.le.explimit) then
        ahigherNT=5.d0*(2.d0**exp1)*(2.d0*3.d0**exp2 + 1)
        write (Iout,1011) dint(ahigherNT+over)
       else
        write (Iout,1021) exp1,exp2
       endif

      endif

C IPR fullerenes
C     Correct upper and lower limit
      if(number_vertices.eq.60.or.number_vertices.ge.70) then
      if(number_vertices.le.122) then
       ifield=number_vertices/2-29
       Write(Iout,1010) IHamIPRmin(ifield),IHamIPRmax(ifield)

      else

C     Approximate number of IPR fullerenes
       explow=1.136977d-1*dAtom
       exphigh=1.230625d-1*dAtom
       If(exphigh.lt.6.d2) then
        fullIPRh=5.698541d-1*dexp(exphigh)*1.2d0
        fullIPRl=1.050204d0*dexp(explow)/1.2d0
       endif
       If(exphigh.lt.3.2d1) then
        write (Iout,1014) dint(fullIPRl),dint(fullIPRh)
       else
        if(exphigh.lt.6.d2) then
         write (Iout,1015) fullIPRl,fullIPRh
        else
         fullexph=dlog(5.698541d-1*1.2d0)*exphigh
         fullexpl=dlog(1.050204d0/1.2d0)*explow
         write (Iout,1025) fullexpl,fullexph
        endif
       endif
 
      endif
      endif

C     Epstein upper limit
      power=dAtom/3.d0
      if(power.lt.explimit) then
       ulepstein=2.d0**power
       write (Iout,1005) dint(ulepstein+over)
      else
       write (Iout,1000) power
      endif

C     Limit for number of atoms
      if(iprintf.eq.0) RETURN
      if(number_vertices.gt.Ihuge) then
       write (Iout,1009) Ihuge 
       RETURN
      endif
      write (Iout,1001) 
C     This is only good for C20, C24 already causes integer overflow
      if(number_vertices.eq.20) then
C     Produce symmetric A^2
      do i=1,number_vertices
      do j=1,i
      IM(i,j)=0
      do k=1,number_vertices
      IM(i,j)=IM(i,j)+IA(i,k)*IA(k,j)
      IM(j,i)=IM(i,j)
      enddo
      enddo
      enddo
C     Now loop do (A^2)^8
      do i=1,number_vertices
      do j=1,i
      IMF(i,j)=IM(i,j)
      IMF(j,i)=IM(j,i)
      IMF1(i,j)=0
      IMF1(j,i)=0
      end do
      end do

      nloop=(number_vertices-2)/2-1
      do loop=1,nloop
       do i=1,number_vertices
       do j=1,number_vertices
       do k=1,number_vertices
       IMF1(i,j)=IMF1(i,j)+IMF(i,k)*IM(k,j)
       enddo
       enddo
       enddo
       do i=1,number_vertices
       do j=1,number_vertices
       IMF(i,j)=IMF1(i,j)
       IMF1(i,j)=0
       enddo
       enddo
      enddo

C     Now last multiplication with A
      do i=1,number_vertices
      do j=1,number_vertices
      do k=1,number_vertices
      IMF1(i,j)=IMF1(i,j)+IMF(i,k)*IA(k,j)
      enddo
      enddo
      enddo
C     Now print number of paths for all vertices
      write (Iout,1002)
      do i=1,number_vertices
      write (Iout,1003) i,(IMF1(i,j),j=1,10)
      enddo
      write (Iout,1004)
      do i=1,number_vertices
      write (Iout,1003) i,(IMF1(i,j),j=11,20)
      enddo

      else
C     For larger matrices using Jordan decomposition is recommended where
C     A=LDL' where D is the diagonal matrix containing the eigenvalues of A
C     L is the matrix of eigenvector, L' its transposed
C     Then A^n = L D^n L' where D^n is the diagonal matric containing all
C     eigenvalues d^n (that is to the nth power)
C     This procedure is much faster and good for large graphs, but real 
C     matrix elements instead of integers are produced.
C     Diagonalize adjacency matrix A
      write (Iout,1007) 
      do i=1,number_vertices
      do j=1,i
      A(i,j)=dfloat(IA(i,j))
      A(j,i)=A(i,j)
      enddo
      enddo
      call tred2(A,number_vertices,Nmax,evec,df)
      call tqli(evec,df,number_vertices,Nmax,A)
C     Calculate A^(n-1) = L D^(n-1) L but only printing the adjacent vertices
C     NP values
      mpower=number_vertices-1
      do i=1,number_vertices
      evec(i)=evec(i)**mpower
      enddo

      limit=5
      ic=0
      do i=1,number_vertices
      do j=1,i-1
      if(IA(i,j).eq.1) then
      ic=ic+1
      IS1(ic)=i
      IS2(ic)=j
      amat=0.d0
      do k=1,number_vertices
      amat=amat+A(i,k)*evec(k)*A(j,k)
      enddo
      APN(ic)=amat
      endif
      if(ic.eq.limit) then
      ic=0
      Write(IOUT,1006) (IS1(l),IS2(l),APN(l),l=1,limit)
      endif
      enddo
      enddo
      if(ic.ne.0) Write(IOUT,1006) (IS1(l),IS2(l),APN(l),l=1,ic)
      endif
 1000 Format(1X,'Epstein upper limit for Hamiltonian cycles in '
     1 'cubic graphs: 2**',F16.4)
 1001 Format(/1X,'Calculate the number of paths (PN) of length (n-1) '
     1 '(n= number of vertices) between vertices i and j'
     1 /1X,'(elements of the (n-1) th power of the adjacency matrix)',
     1 /1X,'Matrix contains numbers PN of distinct paths ',
     1 '(including all Hamiltonian cycles):'
     1 /1X,'(Note this includes the possibility going through vertices'
     1 ' and edges multiple times)')
 1002 Format(/,'   i',5X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',
     1 9X,'8',9X,'9',9X,'10')
 1003 Format(1X,I3,10I10)
 1004 Format(/,'   i',5X,'11',8X,'12',8X,'13',8X,'14',8X,'15',8X,'16',
     1 8X,'17',8X,'18',8X,'19',8X,'20')
 1005 Format(1X,'Epstein upper limit for Hamiltonian cycles in '
     1 'cubic graphs: ',F20.0)
 1006 Format(1X,5('('I3,',',I3,')',D21.14,','))
 1007 Format(1X,'Only matrix elements of adjacent vertices are printed')
 1009 Format(1X,'Number of atoms exceeds ',I3,', change Ihuge value ',
     1 ' (if you dare)') 
 1010 Format(1X,'Exact limits for Hamiltonian cycles for IPR ',
     1 'fullerenes. Upper limit =',I7,', lower limit= ',I7)
 1011 Format(1X,'Estimated upper limit for Hamiltonian cycles in '
     1 'fullerene graphs: ',F20.0)
 1012 Format(1X,'Estimated lower limit for Hamiltonian cycles in '
     1 'fullerene graphs: ',F20.0)
 1013 Format(1X,'Estimated lower limit for Hamiltonian cycles in '
     1 'fullerene graphs: 5*2**',F16.4)
 1014 Format(1X,'Approximate number of Hamiltonian cycles in IPR '
     1 'fullerene graphs: between',F20.0,' and',F20.0)
 1015 Format(1X,'Approximate number of Hamiltonian cycles in IPR '
     1 'fullerene graphs: between ',D20.10,' and ',D20.10)
 1016 Format(1X,'Exact limits for Hamiltonian cycles. Upper ',
     1 'limit =',I7,', lower limit= ',I7)
 1021 Format(1X,'Estimated upper limit for Hamiltonian cycles in '
     1 'fullerene graphs: 5*2**',F16.4,' * (2*3**',F16.4,' + 1)')
 1025 Format(1X,'Approximate number of Hamiltonian cycles in IPR '
     1 'fullerene graphs: between appr. e**a and e**b with a= ',
     1  D22.14,' and b= ',D22.14)
      RETURN
      END
