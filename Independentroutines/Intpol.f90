!///////////////////////////////////////////////////////////////////////////////////////////////
!////////                                                                               ////////
!////////     Intpol.f90                                                                ////////
!////////                                                                               ////////
!////////     contains:                                                                 ////////
!////////                                                                               ////////
!////////           SUBROUTINE polin3()                                                 ////////
!////////           SUBROUTINE polint()                                                 ////////
!////////                                                                               ////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!      
      SUBROUTINE polin3(x1a,x2a,x3a,ya,m,n,o,x1,x2,x3,y,dy)
!
!     3-d Polynomial interpolation
!     date: 9/9/02
!     routines called: polint
!     uses: none
!     applicability: all
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      INTEGER :: m,n,o,NMAX,MMAX,OMAX
      DOUBLE PRECISION    ::  dy,x1,x2,x3,y,x1a(m),x2a(n),x3a(o),ya(m,n,o)
      PARAMETER (NMAX=20,MMAX=20,OMAX=20)
! CU    USES polint
      INTEGER :: j,k
      DOUBLE PRECISION    :: ymtmp(MMAX),yntmp(NMAX),yotmp(OMAX)
      do 12 j=1,m
        do 11 k=1,n
	    do 14 l=1,o
            yotmp(l)=ya(j,k,l)
14        continue
          call polint(x3a,yotmp,o,x3,yntmp(k),dy)
11      continue
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
12    continue
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      END
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE polint(xa,ya,n,x,y,dy)
    
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      INTEGER :: n,NMAX
      DOUBLE PRECISION  ::  dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=20)
      INTEGER i,m,ns
      DOUBLE PRECISION :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)

          den=ho-hp
          if(den.eq.0.) pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
