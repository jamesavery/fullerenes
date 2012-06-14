      SUBROUTINE HamiltonCyc(N,maxiter,Iout,nbatch,
     1 A,Nhamilton)
      use config
C     Subroutine from Darko Babic to create Hamitonian cycles
C      optimized for program Isomer
      integer list(Nmax,3,3),path(0:Nmax+1),stack(3*Nmax),pos(Nmax)
      integer x(0:Nmax),y(0:Nmax),saved(Nmax)
      integer i,j,k,l,n,m,last,next,ngb1,ngb2,jlast,jnext,jngb1,jngb2
      integer ptr,prev,oldptr,cur,prv,nxt,ngb,diff,maxdif,relk,relbr
      logical occ(Nmax),pass(Nmax),end(Nmax),flag,better
      integer ic3(Nmax,3),A(Nmax,Nmax)
      ifirst=0 
      nhamilton=0
      maxN=30
      nbatch=0

      do i=1,n
        k=0
         do j=1,n
           if(A(I,J).eq.1) then
            k=k+1
            ic3(i,k)=j
            endif
         end do
      end do
      
      do i=1,n
         do j=1,3
            list(i,j,1)=ic3(i,j)
         end do
      end do
      
      do i=1,n
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

      do i=1,n
         pass(i)=.false.
         occ(i)=.false.
         end(i)=.false.
      end do
      do i=0,Nmax+1
         path(i)=0
      end do
      x(0)=n

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
      if (l.eq.n-1) then
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
      end

      SUBROUTINE Hamilton(n,Iout,iprint,maxiter,IC3)
C     Subroutine from Darko Babic to create Hamitonian cycles
C     and the IUPAC name of a fullerene
      use config
      integer list(Nmax,3,3),path(0:Nmax+1),stack(3*Nmax),pos(Nmax)
      integer bridge(Nmax),x(0:Nmax),y(0:Nmax),saved(Nmax)
      integer i,j,k,l,n,m,last,next,ngb1,ngb2,jlast,jnext,jngb1,jngb2
      integer ptr,prev,oldptr,cur,prv,nxt,ngb,diff,maxdif,relk,relbr
      logical occ(Nmax),pass(Nmax),end(Nmax),flag,better
      integer ic3(Nmax,3)
      
      ifirst=0 
      nhamilton=0
      maxN=30
      if(n.lt.maxN) maxN=n
      write (Iout,1009) maxiter
      do i=1,n
         do j=1,3
            list(i,j,1)=ic3(i,j)
         end do
      end do

      do i=1,n
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

      do i=1,n
         pass(i)=.false.
         occ(i)=.false.
         end(i)=.false.
      end do
      do i=0,Nmax+1
         path(i)=0
      end do
      x(0)=n

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
      if (l.eq.n-1) then
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
      if(n.gt.30) then
      do I=31,n,30
      jmax=I+29
      if(jmax.gt.n) jmax=n
      write (Iout,1001) (path(j),j=I,jmax)
      enddo
      endif
      endif
      ifirst=1
         do j=1,n           
            pos(path(j))=j
         end do
         path(0)=path(n)
         path(n+1)=path(1)

         maxdif=n
         do j=1,n
            cur=path(j)
            prv=path(j-1)
            nxt=path(j+1)
            do k=1,3
               ngb=list(cur,k,1)
               if (ngb.ne.prv.and.ngb.ne.nxt) then
                  bridge(j)=pos(ngb)
                  diff=abs(bridge(j)-j)
                  if (n-diff.gt.diff) diff=n-diff
                  if (maxdif.gt.diff) maxdif=diff
                  go to 11
               endif
            end do
   11    continue
         end do

         maxdif=maxdif-1
         if (maxdif.gt.x(0)) go to 6

         do 13 j=1,n
            better=.false.
            diff=mod(n+bridge(j)-j,n)-1
            if (diff.eq.maxdif) then
               if (maxdif.lt.x(0)) then
                  x(0)=maxdif
                  y(0)=n-maxdif-2
                  better=.true.
               endif
               i=0
               k=j
               do 14 m=1,n-1
                  k=mod(k,n)+1
                  relk=mod(n+k-j,n)+1
                  relbr=mod(n+bridge(k)-j,n)+1
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
               do m=1,n
                  saved(m)=path(m)
               end do
            end if
   17       better=.false.
            diff=n-diff-2
            if (diff.eq.maxdif) then
               if (maxdif.lt.x(0)) then
                  x(0)=maxdif
                  y(0)=n-maxdif-2
                  better=.true.
               endif
               i=0
               k=j
               do 18 m=1,n-1
                  k=mod(n+k-2,n)+1
                  relk=mod(n+j-k,n)+1
                  relbr=mod(n+j-bridge(k),n)+1
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
               do m=1,n
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
      write (Iout,1000) (x(i),i=0,(n-2)/2)
      write (Iout,1000) (y(i),i=0,(n-2)/2)
      write (Iout,1006)
      write (Iout,1008) (saved(j),j=1,n)
 1000 format (1X,i3,i5,98I3)
 1001 format (9X,30(I4,'-'))
 1002 format (/1X,'There is no Hamiltonian cycle in the graph!')
 1003 format (/1X,I10,' Hamiltonian cycles found',
     1 //1X,'Half-ring sizes and IUPAC superscripts:')
 1004 format (1X,I10,': ',30(I4,'-'))
 1005 format (/1X,'Distinct Hamiltonian cycles:',/6X,'NH    vertices',
     1 /1X,20('-'))
 1006 format (/1X,'The best Hamiltonian cycle:')
 1007 format (1X,'Next Hamiltonian cycles:')
 1008 format (30(I4,'-'))
 1009 format (/1X,'Calculate Hamiltonian cycles, half-ring sizes'
     1 ' and IUPAC superscripts: (D. Babic, J. Chem. Inf. Comput. Sci.'
     1 ' 35, 515-526 (1995).)',/1X,'Maximum allowed iteration: ',I10)
 1010 format (I10,' Maximum Hamiltonian cycles reached: Return')
 1011 format (1X,' Hamiltonian cycle detected')
      return
      end

      SUBROUTINE Paths(MAtom,IOUT,IA,A,evec,df)
      use config
      IMPLICIT REAL*8 (A-H,O-Z)
C Calculate the number of all possible distinct paths from the
C adjaceny matric A(i,j) by producing the (n-1)-th power
C of A. In this case n is Matom. This gives all possible walks
C of length (n-1).
      DIMENSION A(Nmax,Nmax),evec(Nmax),df(Nmax)
      DIMENSION IA(Nmax,Nmax),IM(Nmax,Nmax),IM2(Nmax,Nmax)
      DIMENSION IMF(Nmax,Nmax),IMF1(Nmax,Nmax)
      DIMENSION IS1(10),IS2(10),APN(10)
      DIMENSION B(Nmax,Nmax)
      Real*8 upperschwerd,lowerschwerd
      Data Ihuge/180/
C     Epstein upper limit
      dAtom=dfloat(MAtom)
      power=dAtom/3.d0
      if(power.lt.1.d3) ulepstein=2.d0**power
      If(power.lt.31.01d0) then
       ilepstein=dint(ulepstein)
       write (Iout,1005) ilepstein
      else
      if(power.lt.1.d3) write (Iout,1000) ulepstein
      if(power.gt.1.d3) write (Iout,1016) power
      endif
C     Schwerdtfeger upper and lower limit
      aupper=1.7205d-1
      bupper=1.466d0
      alower=1.188d-1
      blower=-3.720d-1
      powerupper=aupper*datom+bupper
      powerlower=alower*datom-blower
      upperschwerd=2.d0**powerupper
      lowerschwerd=2.d0**powerlower
      If(powerupper.lt.31.01d0) then
       iupperschwerd=dint(upperschwerd)
       write (Iout,1010) iupperschwerd
      else
       write (Iout,1011) upperschwerd
      endif
      If(powerlower.lt.31.01d0) then
       ilowerschwerd=dint(lowerschwerd)
       If(ilowerschwerd.lt.18) ilowerschwerd=18
        write (Iout,1012) ilowerschwerd
      else
       write (Iout,1013) lowerschwerd
      endif

C     Approximate number of IPR fullerenes
      if(MAtom.eq.60.or.MAtom.ge.70) then
      exphigh=1.230625d-1*dAtom
      explow=1.136977d-1*dAtom
      fullIPRh=5.698541d-1*dexp(exphigh)*1.2d0
      fullIPRl=1.050204d0*dexp(explow)/1.2d0
      If(fullIPRh.lt.2.d9) then
       ifullIPRh=dint(fullIPRh)
       ifullIPRl=dint(fullIPRl)
       write (Iout,1014) ifullIPRl,ifullIPRh
      else
       write (Iout,1015) fullIPRl,fullIPRh
       if(fullIPRh.gt.1.d10) then
        write (Iout,1008)
        Return
       endif
      endif
      endif
C     Limit for number of atoms
      if(Matom.gt.Ihuge) then
       write (Iout,1009) Ihuge 
       RETURN
      endif
      write (Iout,1001) 
C     This is only good for C20, C24 already causes integer overflow
      if(MAtom.eq.20) then
C     Produce symmetric A^2
      do i=1,MAtom
      do j=1,i
      IM(i,j)=0
      do k=1,MAtom
      IM(i,j)=IM(i,j)+IA(i,k)*IA(k,j)
      IM(j,i)=IM(i,j)
      enddo
      enddo
      enddo
C     Now loop do (A^2)^8
      do i=1,MAtom
      do j=1,i
      IMF(i,j)=IM(i,j)
      IMF(j,i)=IM(j,i)
      IMF1(i,j)=0
      IMF1(j,i)=0
      end do
      end do

      nloop=(MAtom-2)/2-1
      do loop=1,nloop
       do i=1,MAtom
       do j=1,MAtom
       do k=1,MAtom
       IMF1(i,j)=IMF1(i,j)+IMF(i,k)*IM(k,j)
       enddo
       enddo
       enddo
       do i=1,MAtom
       do j=1,MAtom
       IMF(i,j)=IMF1(i,j)
       IMF1(i,j)=0
       enddo
       enddo
      enddo

C     Now last multiplication with A
      do i=1,MAtom
      do j=1,MAtom
      do k=1,MAtom
      IMF1(i,j)=IMF1(i,j)+IMF(i,k)*IA(k,j)
      enddo
      enddo
      enddo
C     Now print number of paths for all vertices
      write (Iout,1002)
      do i=1,MAtom
      write (Iout,1003) i,(IMF1(i,j),j=1,10)
      enddo
      write (Iout,1004)
      do i=1,MAtom
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
      do i=1,MAtom
      do j=1,i
      A(i,j)=dfloat(IA(i,j))
      A(j,i)=A(i,j)
      enddo
      enddo
      call tred2(A,Matom,Nmax,evec,df)
      call tqli(evec,df,Matom,Nmax,A)
C     Calculate A^(n-1) = L D^(n-1) L but only printing the adjacent vertices
C     NP values
      mpower=Matom-1
      do i=1,MAtom
      evec(i)=evec(i)**mpower
      enddo

      limit=5
      ic=0
      do i=1,MAtom
      do j=1,i-1
      if(IA(i,j).eq.1) then
      ic=ic+1
      IS1(ic)=i
      IS2(ic)=j
      amat=0.d0
      do k=1,Matom
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
 1000 Format(/1X,'Epstein upper limit for Hamiltonian cycles in '
     1 'cubic graphs: ',D22.14)
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
 1005 Format(/1X,'Epstein upper limit for Hamiltonian cycles in '
     1 'cubic graphs: ',I12)
 1006 Format(1X,5('('I3,',',I3,')',D21.14,','))
 1007 Format(1X,'Only matrix elements of adjacent vertices are printed')
 1008 Format(1X,'Number of paths of length (n-1) exceeds computer'
     1 ' real number limit --> Return') 
 1009 Format(1X,'Number of atoms exceeds ',I3,', change Ihuge value ',
     1 ' (if you dare)') 
 1010 Format(1X,'Schwerdtfeger upper limit for Hamiltonian cycles in '
     1 'fullerene graphs: ',I12)
 1011 Format(1X,'Schwerdtfeger upper limit for Hamiltonian cycles in '
     1 'fullerene graphs: ',D22.14)
 1012 Format(1X,'Schwerdtfeger lower limit for Hamiltonian cycles in '
     1 'fullerene graphs: ',I12)
 1013 Format(1X,'Schwerdtfeger lower limit for Hamiltonian cycles in '
     1 'fullerene graphs: ',D22.14)
 1014 Format(1X,'Approximate number of Hamiltonian cycles in IPR '
     1 'fullerene graphs: between appr.',I12,' and ',I12)
 1015 Format(1X,'Approximate number of Hamiltonian cycles in IPR '
     1 'fullerene graphs: between appr.',D22.14,' and ',D22.14)
 1016 Format(/1X,'Epstein upper limit for Hamiltonian cycles in '
     1 'cubic graphs (only power to base 2 given): ',D22.14)
      RETURN
      END
