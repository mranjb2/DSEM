!> \file dtcalc.f90
!! Calculation of timestep sizes
!///////////////////////////////////////////////////////////////////////
!
!> @brief applicability: all nonlinear Euler eqns stac3m version.
!!
!!     note: computed this way, the metric terms are off by an error
!!     of order delta_x, but this shouldn't affect things too much.
      SUBROUTINE dtcalc(dt,d) 
!
!......................................................................
!     date: 11/24/98
!     routines called: none                                             
!     applicability: all nonlinear Euler eqns stac3m version
!     note: computed this way, the metric terms are off by an error
!     of order delta_x, but this shouldn't affect things too much.
!......................................................................
!
      USE domain_definition
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d 
!
!     ------------------
!     inviscid time step
!     ------------------
!
      DO l = 1,d%ncg(3)
         dcz = d%cx(l+1,3) - d%cx(l,3) 
         DO m = 1,d%ncg(2) 
            dcy = d%cx(m+1,2) - d%cx(m,2) 
            DO n = 1,d%ncg(1)
               dcx = d%cx(n+1,1) - d%cx(n,1)

               u0 = abs(d%Q(n,m,l,2)/d%Q(n,m,l,1))
               v0 = abs(d%Q(n,m,l,3)/d%Q(n,m,l,1))
               w0 = abs(d%Q(n,m,l,4)/d%Q(n,m,l,1))

               p = (gamma-1.d0)*(d%Q(n,m,l,5) - 0.5d0*d%Q(n,m,l,1)*(u0**2 + v0**2 + w0**2))                       
               a = sqrt(gamma*p/d%Q(n,m,l,1)) 

               uu = (u0 + a)*d%gmet(1,1,n,m,l) + (v0 + a)*d%gmet(1,2,n,m,l) + (w0 + a)*d%gmet(1,3,n,m,l)
               vv = (u0 + a)*d%gmet(2,1,n,m,l) + (v0 + a)*d%gmet(2,2,n,m,l) + (w0 + a)*d%gmet(2,3,n,m,l)
               ww = (u0 + a)*d%gmet(3,1,n,m,l) + (v0 + a)*d%gmet(3,2,n,m,l) + (w0 + a)*d%gmet(3,3,n,m,l)

               xlam = abs(uu*d%jacob(n,m,l)) 
               ylam = abs(vv*d%jacob(n,m,l))
               zlam = abs(ww*d%jacob(n,m,l))

               dti = xlam/dcx + ylam/dcy + zlam/dcz + 1.d-14 
               dt =  min(dt,1.d0/dti)
            END DO
         END DO
      END DO
!
      RETURN 
      END SUBROUTINE dtcalc
