      SUBROUTINE spwindup(IER,number_faces,MP,D1,S,JP,FreeRing)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION D1(6,MMAX),S(MMAX),JP(12),FreeRing(6)
C This routine tries to find the spiral from a preset of
C three connected rings stored in S(1), S(2) and S(3) using
C information on ring fusions from the reduced dual matrix D1
C If pentagon is found MP increases by 1.
C IER=1 implies that no spiral found
C This algorithm is much faster than the one in version 4.2

      IER=0
C   Big loop from ring 4 onwards to the end of the spiral
C---- Start loop
      nsmall=1
      faces: Do i=4,number_faces

C       First collect next possible faces connected to the previous 
C       5- or 6-ring (nloop=5 or 6) numbered IP=S(i-1) 
C       and make sure it is not one of the previous rings
C       in the spiral
        IP=S(i-1)
        nloop=6
        nring=0
        if(D1(6,IP).eq.0) nloop=5
        loops: do j=1,nloop
C         Get all adjacent rings to previous one
          nr=D1(j,IP)
C         Make sure it is not one in the existing spiral
          do j1=I-2,1,-1
            if(nr.eq.S(j1)) cycle loops
          enddo
C       Collect them
        nring=nring+1
        FreeRing(nring)=nr
        enddo loops
          
C       Now it needs to be connected to a previous ring
C       Last 2 are not needed
        do k=nsmall,i-3
          KP=S(k)
          do k1=1,6
            do j2=1,nring
              nr=FreeRing(j2)
              if(D1(k1,KP).eq.nr) then
                nsmall=k
                S(i)=nr
                if(D1(6,nr).eq.0) then
                  MP=MP+1
                  JP(MP)=i
                endif
                cycle faces
              endif
            enddo
          enddo
        enddo

        if(S(i).eq.0) then
          IER=1   ! Spiral has dead end
          Return
        endif
      enddo faces

C  Finally success, spiral found
      Return
      END


      SUBROUTINE Windup(M,IPR,IER,S,D)
      use config
      IMPLICIT INTEGER (A-Z)
      DIMENSION S(MMAX),D(MMAX,MMAX)
      DIMENSION R(MMAX),C(MMAX,6)
C     This subroutine attempts to wind up an input spiral S into
C     a fullerene dual (face) adjacency matrix D. It returns with      
C     IER = P if the spiral shorts or is discovered to be open-ended      
C     after P pentagons have been added. Otherwise IER = 0 on return.
      J=1
      C(1,1)=2
      C(2,1)=1
      R(1)=2
      R(2)=2
      E=1
      P=12-S(1)-S(2)
      DO 5 K=3,M-1
         P=P+6-S(K)
         R(K)=1
         I=K-1
 1       IF(IPR.EQ.1.AND.S(I).EQ.5.AND.S(K).EQ.5) GO TO 10      
         IF(R(K).GE.S(K)) GO TO 10      
         C(I,R(I))=K            ! Connect face K to the last open face I
         C(K,R(K))=I            ! in the preceding spiral      
         R(I)=R(I)+1
         R(K)=R(K)+1
         IF(R(I).GT.S(I)) THEN
            L=I-1               ! If this closes face I update I and go
            DO I=L,J+1,-1       ! back to connect face K to the new I      
               IF(R(I).LE.S(I)) GO TO 1
            enddo
            GO TO 10
         ENDIF
 3       IF(IPR.EQ.1.AND.S(J).EQ.5.AND.S(K).EQ.5) GO TO 10
         IF(R(K).GE.S(K)) GO TO 10
         C(J,R(J))=K            ! Connect face K to the first open face J
         C(K,R(K))=J            ! in the preceding spiral 
         R(J)=R(J)+1
         R(K)=R(K)+1
         IF (R(J).GT.S(J)) THEN
            L=J+1               ! If this closes face J update J and go
            DO J=L,I-1,+1       ! back to connect face K to the new J 
               IF(R(J).LE.S(J)) GO TO 3
            enddo
            GO TO 10
         ENDIF
         H=K-P
         E=E+R(K)-1             ! Use Euler's theorem to streamline the 
         V=3+2*P+3*H-E          ! search. F is a lower bound on the # of 
         F=(V+1)/2+1            ! additional faces required for closure 
         IF(F.GT.M-K) GO TO 10
 5    CONTINUE
      P=12
      R(M)=1
      DO 6 K=J,M-1
         IF (R(K).LT.S(K)) GO TO 10
         IF (R(K).GT.S(K)) GO TO 6
         IF (R(M).GT.S(M)) GO TO 10
         IF (IPR.EQ.1.AND.S(K).EQ.5.AND.S(M).EQ.5) GO TO 10
         C(K,R(K))=M            ! Connect face M to all remaining
         C(M,R(M))=K            ! open faces (including face M-1)
         R(K)=R(K)+1
         R(M)=R(M)+1
 6    CONTINUE
      IF(R(M).LE.S(M)) GO TO 10
      P=0                       ! Successful spiral
      DO J=1,M                  ! Form dual adjacency matrix in D
         DO I=1,M
            D(I,J)=0
         enddo
         DO K=1,S(J)
            I=C(J,K)
            D(I,J)=1
         enddo
      enddo
 10   IER=P
      RETURN
      END
