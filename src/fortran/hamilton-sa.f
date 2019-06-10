      SUBROUTINE HamiltonCyc(maxiter,Iout,nbatch,IC3,nhamilton)
      use config
!---------------------------------------------------------------------------!
!  This routine counts the number of Hamiltonian cycles using the           !
!  back-track algorithm of Darko Babic. It has been optimized for speed     ! 
!  and is called by the following subroutines:                              !
!      Main program Fullerene                                               !
!      Printdatabase                                                        !
!      SpiralRestart                                                        !
!                                                                           !
!  Input: maxiter,Iout,A                                                    !
!  Output: nbatch,nhamilton                                                 !
!---------------------------------------------------------------------------!
      implicit none
      integer maxiter
      integer iout
      integer nbatch ! number of batches
      integer IC3(nmax,3) ! is only of size N, but is input as nmax, which makes this size necassary at this point
      integer nhamilton ! number of found cycles
      integer list(number_vertices,3,3)
      integer path(0:number_vertices+1)
      integer stack(2*number_vertices) !FIXME: find correct minimal size. N+5?
      integer i,j,k,l,l1,m
      integer last,next,ngb1,ngb2
      integer jlast,jnext,jngb1,jngb2
      integer ptr,prev,oldptr
      logical occ_bool(number_vertices) ! true for vertices that have been added to the path, else false
      logical pass_bool(number_vertices) ! true for vertices that have been passed by the path and thus have one less valency left, else false
      logical end_bool(number_vertices) ! true for the two vertices that are neighbours of the first one (and could be the end of the path), else false
      logical flag ! second neighbour of starting point has been used, i.e. there is a defined last vertex

      nhamilton=0
      nbatch=0

c copy adjacency list from IC3(*,*) to list(*,*,1)
      do i=1,number_vertices
        do j=1,3
          list(i,j,1)=IC3(i,j)
        enddo
      enddo
      
c list(*,*,1): neighbour list
c list(*,*,2): which second neighbour (dist=2 in the graph) is the first != self
c list(*,*,3): which second neighbour (dist=2 in the graph) is the second != self
      do i=1,number_vertices
        do j=1,3
          k=list(i,j,1) 
          l=1
          do m=1,3
            if (list(k,m,1).ne.i) then
              l=l+1
              list(i,j,l)=m
            endif
          enddo
        enddo
      enddo

c init
      do i=1,number_vertices
        pass_bool(i)=.false.
        occ_bool(i)=.false.
        end_bool(i)=.false.
      enddo
      do i=0,number_vertices+1
        path(i)=0
      enddo
      stack(1)=0
      stack(2)=0
      stack(3)=1
      stack(4)=1
      stack(5)=2
      oldptr=3
      ptr=6
      occ_bool(1)=.true.
      last=1
      next=list(1,1,1)
      jnext=1
      end_bool(list(1,2,1))=.true.
      end_bool(list(1,3,1))=.true.
      flag=.false.
      l=1
      path(1)=1

c Start algorithm
      goto 5

c choose next vertex in path
    4 jngb1=list(prev,jlast,2)
      jngb2=list(prev,jlast,3)
      ngb1=list(last,jngb1,1)
      ngb2=list(last,jngb2,1)

      if (occ_bool(ngb1)) then
        next=ngb2
        jnext=jngb2
      else
        if (occ_bool(ngb2)) then
          next=ngb1
          jnext=jngb1
        else
          if (pass_bool(ngb1)) then
            if (pass_bool(ngb2).or.end_bool(ngb1)) goto 6
            next=ngb1
            jnext=jngb1
            pass_bool(ngb2)=.true.
            stack(ptr)=ngb2
            ptr=ptr+1
          else
            if (pass_bool(ngb2)) then
              if (end_bool(ngb2)) goto 6
              next=ngb2
              jnext=jngb2
              pass_bool(ngb1)=.true.
              stack(ptr)=ngb1
              ptr=ptr+1
            else
              next=ngb1
              jnext=jngb1
              pass_bool(ngb2)=.true.
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
        nhamilton=nhamilton+1 ! one circuit completed
         if(nhamilton.ge.maxiter) then
           if(maxiter.lt.1000000000) then
             nbatch=nbatch+1
             Write(Iout,1000) nbatch,nhamilton
             nhamilton=0
           else 
             Return
           endif
        endif
      endif

      if (end_bool(next)) then ! is the next vertex adjacent to the starting point?
         ! is 'next' the second of the two ends or has the next vertex been passed before? : backtrack (also backtrack if a cycle was completet)
         if (flag.or.pass_bool(next)) goto 6
         flag=.true.
      endif

      l=l+1
      occ_bool(next)=.true.

      prev=last
      last=next
      jlast=jnext
      goto 4

c backtrack
    6 l1=stack(oldptr+1)
      do i=oldptr+3,ptr-1
        pass_bool(stack(i))=.false.
      enddo
      do i=l1+1,l
        occ_bool(path(i))=.false.
        if (end_bool(path(i))) flag=.false.
      enddo
      ptr=oldptr
      oldptr=stack(ptr)
      pass_bool(path(l1+1))=.true.
      stack(ptr)=path(l1+1)
      ptr=ptr+1
      last=path(l1)
      jnext=stack(ptr+1)
      next=list(last,jnext,1)
      pass_bool(next)=.false.

      l=l1
      if(oldptr.gt.0) goto 5

1000  Format(1X,'Batch ',I3,' of Hamiltonian cycles: ',I10)

      return
      END SUBROUTINE HamiltonCyc
