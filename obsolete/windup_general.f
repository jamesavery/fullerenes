      PROGRAM TEST

      Integer MMAX, M, N, IPR, IER, S, D, Z, Y
c spiral
      DIMENSION S(192)
c dual matrix
      DIMENSION D(192,192)
c jumps
      DIMENSION Z(10,2)

C regular
c      S=(/5,5,5,5,5,6,5,6,6,5,6,5,5,5,5,5/)
C requires skip, jump, leap, whatever
c      S=(/5,5,5,6,5,6,6,5,5,6,5,5,5,5,5,5/) ! 28
C requires skip, jump, leap, whatever
      S=(/5,5,5,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,
     2    6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,
     3    6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,
     4    6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,6,6,6,6,
     5    6,6,6,6,6,6,6,6,6,6, 6,6,6,6,6,6,5,5,6,6, 6,6,6,6,6,6,5,5,6,6,
     6    6,6,6,6,6,6,5,5,6,6, 6,6,6,6,6,5,6,6,6,6, 6,6,5,6,6,6,6,6,6,5,
     7    6,6,6,6,6,6,6,6,6,6, 6,6/) ! 192
      MMAX=192
      N=380
      M=N/2+2
      IPR=0
C number of required jumps
      Y=0
C Z stores position and length of jumps

      CALL WindupGeneral(MMAX,M,IPR,IER,S,D,Z,Y)

      if(ier .eq. 0)then
         write(*,*)'success! ',y,' jumps required.' 
         WRITE(*,*)'z',Z
      else
         WRITE(*,*)'Fail after ',IER,' pentagons.'
      endif
      END



      SUBROUTINE WindupGeneral(MMAX,M,IPR,IER,S,D,Z,Y)
c it should be possible to write windup_general in a much simpler way if there was 
c an equivalent of c++-vectors in fortran:  there would be vector holding all vertices
c whith open valencies in the correct order.  j would be the first entry, j+1 the
c second etc.  l is the last entry, l-1 last but one.  whenever a vertex runs out of
c valencies it is erased from the vector.  Jumping simply moves a number of vertices
c from the beginning to the end (before k is added).  In the end of each cycle k is 
c appended to the list.

c depending of the performance it might be a good idea to have just one 'windup' and
c activate general spiral search by passing a flag.
      IMPLICIT INTEGER (A-Z)
      DIMENSION S(MMAX),D(MMAX,MMAX)
      DIMENSION R(MMAX),C(MMAX,6)
      DIMENSION Z(10,2)
      DIMENSION W(MMAX)
C       This subroutine attempts to wind up an input spiral S into
C       a fullerene dual (face) adjacency matrix D. It returns with      
C       IER = P if the spiral shorts or is discovered to be open-ended      
C       after P pentagons have been added. Otherwise IER = 0 on return.
      J=1
C sparse adjacency matrix of the dual
      C(1,1)=2
      C(2,1)=1
C used valencies per face +1
      R(1)=2
      R(2)=2
      E=1
C auxirialy variables
      X=0
      U=1 !number of skipped valencies +1
C number of used pentagons
      P=12-S(1)-S(2)
      DO 5 K=3,M-1
c         WRITE(*,*)'nn',S(J)-R(J)+1, S(J+1)-R(J+1)+1, S(J+2)-R(J+2)+1,
c     2       S(J+3)-R(J+3)+1, x
         IF(S(J)-R(J)+1.EQ.2 .and. k.lt.m-1) THEN
            W(U)=K-1
            U=U+1
            DO WHILE(S(J+X)-R(J+X)+1.EQ.2)
c               write(*,*)'ijkx',i,j,k,x
c               write(*,*)'found 2 ','s(j+x)',S(J+x),'r(j+x)',R(J+x),
c     2                   'x',X
               W(U)=J+X ! there are two open valencies, ergo twice
               U=U+1
               W(U)=J+X
               U=U+1
               X=X+1
            ENDDO
            U=U-1
c            write(*,*)'s(j)',S(J),'r(j)',R(J),'x',X,'u',u,'w',w
            IF(S(J+X)-R(J+X)+1.EQ.1 .AND. S(J+X+1)-R(J+X+1)+1.EQ.1 .AND.
     2            S(J+X+2)-R(J+X+2)+1.NE.1 .AND. S(K).EQ.5) THEN
c               write(*,*)'jump5 ',K,X+1,
c     2                 'U is ',u,' (should be ',2*x+1,')'
c               W(U)=J+X
c               U=U+1
               P=P+1
               R(K)=1
               I=J+X-1
               J=J+X
c               write(*,*)'pijkxy',p,i,j,k,x,y
               Y=Y+1
               Z(Y,1)=K
               Z(Y,2)=X+1
c               write(*,*)'used pentagons: ',P,'connection 1: ',I,
c     2            'connection 2: ',J,'x: ',x
               X=0
               GO TO 1
            ELSEIF(S(J+X)-R(J+X)+1.EQ.1 .AND. S(J+X+1)-R(J+X+1)+1.EQ.1
     2            .AND. S(J+X+2)-R(J+X+2)+1.EQ.1
     3            .AND. S(J+X+3)-R(J+X+3)+1.NE.1 .AND. S(K).EQ.6) THEN
c               write(*,*)'jump6 ',K,X+1
c               W(U)=J+X
c               U=U+1
               R(K)=1
               I=J+X-1
               J=J+X
               Y=Y+1
               Z(Y,1)=K
               Z(Y,2)=X+1
c               write(*,*)'used pentagons: ',P,'connection 1: ',I,
c     2            'connection 2: ',J,'x: ',x
               X=0
               GO TO 1
            ELSEIF(S(J+X)-R(J+X)+1.EQ.1 .AND. S(J+X+1)-R(J+X+1)+1.EQ.1
     2            .AND. S(J+X+2)-R(J+X+2)+1.EQ.1
     3            .AND. S(K).EQ.5) THEN
c               WRITE(*,*)'fail 1.0'
               GO TO 10
            ELSEIF(S(J+X)-R(J+X)+1.EQ.1 .AND. S(J+X+1)-R(J+X+1)+1.EQ.1
     2            .AND. S(J+X+2)-R(J+X+2)+1.EQ.1
     3            .AND. S(J+X+3)-R(J+X+3)+1.EQ.1) THEN
c               WRITE(*,*)'fail 1.1'
               GO TO 10
            ENDIF
            X=0
            U=1
         ENDIF
         P=P+6-S(K)
         R(K)=1
         I=K-1
C connection to K-1 (and possibly K-2 etc)
 1       IF(IPR.EQ.1 .AND. S(I).EQ.5 .AND. S(K).EQ.5) GO TO 10 
         IF(R(K).GE.S(K)) GO TO 10 
c         WRITE(*,*)'section1, ijkxuw',i,j,k,x,u,w
         C(I,R(I))=K    ! Connect face K to the last open face I
         C(K,R(K))=I    ! in the preceding spiral      
         R(I)=R(I)+1
         R(K)=R(K)+1
         IF(R(I).GT.S(I) .and. u.eq.1) THEN    ! we are walking out of a cavity
                                ! If this closes face I update I and go
            DO I=K-2,J+1,-1     ! back to connect face K to the new I      
               IF(R(I).LE.S(I)) GO TO 1
            enddo
c            WRITE(*,*)'fail 2'
            GO TO 10
         ENDIF
C connection to valencies that were skipped while jumping.  These are stored in w, w(u-1) is to be connected next.
         DO WHILE(U.gt.1 .and. C(W(U-1),R(W(U-1))-1).ne.K) ! if there are skipped valencies and if the vertices are not connected yet
            IF(IPR.EQ.1 .AND. S(W(U-1)).EQ.5 .AND. S(K).EQ.5) GO TO 10 ! the usual checking ...
            IF(R(K).GE.S(K)) GO TO 10
c            WRITE(*,*)'section2.0, ijkxuw',i,j,k,x,u,w
c            WRITE(*,*)'section2.01,',c
            C(W(U-1),R(W(U-1)))=K ! the two bond entries
            C(K,R(K))=W(U-1)
            R(W(U-1))=R(W(U-1))+1 ! increment the used valencies
            R(K)=R(K)+1
            U=U-1                 ! decrement the number of skipped valencies
                                    ! we are walking backwards
                                    ! If this closes face J update J and go
                                    ! back to connect face K to the new J 
            IF(R(K) .gt. S(K))THEN
c               WRITE(*,*)'fail 3'
               GO TO 10
            ENDIF
c            WRITE(*,*)'section2.01,',c
c            WRITE(*,*)'section2.1, ijkxuw',i,j,k,x,u,w
         ENDDO
C connection to J (and possibly J+1 etc)
 3       IF(IPR.EQ.1 .AND. S(J).EQ.5 .AND. S(K).EQ.5) GO TO 10
         IF(R(K).GE.S(K)) GO TO 10
c         WRITE(*,*)'section3, ijkxu',i,j,k,x,u
         C(J,R(J))=K    ! Connect face K to the first open face J
         C(K,R(K))=J    ! in the preceding spiral 
         R(J)=R(J)+1
         R(K)=R(K)+1
         IF (R(J).GT.S(J)) THEN   ! we are walking into a cavity
                                ! If this closes face J update J and go
            DO J=J+1,K-2,+1     ! back to connect face K to the new J 
               IF(R(J).LE.S(J)) GO TO 3
            enddo
c            WRITE(*,*)'fail 4'
            GO TO 10
         ENDIF
         H=K-P ! used hexagons, used faces, used pentagons
c         write(*,*)'hkp',h,k,p
         E=E+R(K)-1     ! Use Euler's theorem to streamline the ! used edges
c         write(*,*)'er(k)',e,r(k)
         V=3+2*P+3*H-E  ! search. F is a lower bound on the # of 
c         write(*,*)'vphe',v,p,h,e
         F=(V+1)/2+1    ! additional faces required for closure 
c         write(*,*)'fvmk',f,v,m,k,m-k
         IF(F.GT.M-K)THEN
c            WRITE(*,*)'fail detected via euler'
            GO TO 10
         ENDIF
 5    CONTINUE
      P=12
      R(M)=1
      DO 6 K=J,M-1
         IF (R(K).LT.S(K)) GO TO 10
         IF (R(K).GT.S(K)) GO TO 6
         IF (R(M).GT.S(M)) GO TO 10
         IF (IPR.EQ.1.AND.S(K).EQ.5.AND.S(M).EQ.5) GO TO 10
         C(K,R(K))=M    ! Connect face M to all remaining
         C(M,R(M))=K    ! open faces (including face M-1)
         R(K)=R(K)+1
         R(M)=R(M)+1
 6    CONTINUE
      IF(R(M).LE.S(M)) GO TO 10
      P=0               ! Successful spiral
      DO J=1,M        ! Form dual adjacency matrix in D
         DO I=1,M
            D(I,J)=0
         enddo
         DO K=1,S(J)
            I=C(J,K)
            D(I,J)=1
         enddo
      enddo
 10   IER=P
c      WRITE(*,*)'d',c
      RETURN
      END

