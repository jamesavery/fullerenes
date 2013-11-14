c tridiagonal matrix from a real symmetric matrix
c reduced by all features that are required for eigenvectors
c borrowed and modified from numerical recipes
      SUBROUTINE tred2l(a,n,np,d,e)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 a(np,np),d(np),e(np)
C     a: real symmetric matrix of dimension n
C     np: dimension of array
      do 18 i=n,2,-1
        l=i-1
        h=0.d0
        scale=0.
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+dabs(a(i,k))
11        continue
          if(scale.eq.0.d0)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(dsqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.
            do 15 j=1,l
              g=0.
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
      e(1)=0.d0
      do 24 i=1,n
        d(i)=a(i,i)
24    continue
      return
      END

c eigenvalues from tridiagonal matrix
c reduced by all features that are required for eigenvectors
c borrowed and modified from numerical recipes
      SUBROUTINE tqlil(d,e,n,np)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 d(np),e(np)
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=dabs(d(m))+dabs(d(m+1))
          if (dabs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.300) then
          Print*,' ERROR: Too many iterations in tqlil'
          return
          endif
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=dsqrt(g*g+1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=dsqrt(g*g+f*f)
            e(i+1)=r
            if(r.eq.0.d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.d0
          goto 1
        endif
15    continue
      return
      END

c tridiagonal matrix from real symmetric matrix
c borrowed from numerical recipes
      SUBROUTINE tred2(a,n,np,d,e)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 a(np,np),d(np),e(np)
C     a: real symmetric matrix of dimension n
C     np: dimension of array
      do 18 i=n,2,-1
        l=i-1
        h=0.d0
        scale=0.d0
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+dabs(a(i,k))
11        continue
          if(scale.eq.0.d0)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(dsqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.
            do 15 j=1,l
              a(j,i)=a(i,j)/h
              g=0.d0
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
      d(1)=0.d0
      e(1)=0.d0
      do 24 i=1,n
        l=i-1
        if(d(i).ne.0.)then
          do 22 j=1,l
            g=0.d0
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
        d(i)=a(i,i)
        a(i,i)=1.
        do 23 j=1,l
          a(i,j)=0.d0
          a(j,i)=0.d0
23      continue
24    continue
      return
      END

c eigenvectors and eigenvalues from tridiagonal matrix
c borrowed from numerical recipes
      SUBROUTINE tqli(d,e,n,np,z)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 d(np),e(np),z(np,np)
c	-------------------------------------------
c	QL Algorithm with implicit shifts, to determine the 
c	eigenvalues and eigenvectors of a real, symmetric, tridiagonal matrix
c	D: The N diagonal elements of the matrix.
c	E: Subdiagonal elements of the matrix, with e(1) arbitrary.
c          (On output is destroyed).
c	For eigenvectors: Z must be the identity matrix. 
c	   On output the kth column returns the normalized eigenvector
c	   corresponding to d(k)
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.d0
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=dabs(d(m))+dabs(d(m+1))
          if (dabs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.300) then
          Print*,' ERROR: Too many iterations in tqli'
          return
          endif
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=dsqrt(g*g+1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.d0
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=dsqrt(g*g+f*f)
            e(i+1)=r
            if(r.eq.0.)then
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.d0
          goto 1
        endif
15    continue
      return
      END
