!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE wallbc_isoth(mrtr,flux) 
!
!......................................................................
!     date: 03/14/01
!
!     Compute the wall bc using characteristic splitting,                 
!     by specifying an equal but opposite flow as the external condition.
!     Compute Q(5) using twall
!......................................................................
!
      USE  mortar_definition
      USE  physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property1 = 0.0d0, property2 = 0.0d0
!
!
!     -----------------------------------------
!     Generate the external flow along the face
!     -----------------------------------------
!
      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)

            qnorm  = (mrtr%n_hat(1,i,j)*mrtr%Q(i,j,2,1) + &   ! normal momentum
                      mrtr%n_hat(2,i,j)*mrtr%Q(i,j,3,1) + &
                      mrtr%n_hat(3,i,j)*mrtr%Q(i,j,4,1))
            q_tanx = mrtr%Q(i,j,2,1) - qnorm*mrtr%n_hat(1,i,j)
            q_tany = mrtr%Q(i,j,3,1) - qnorm*mrtr%n_hat(2,i,j)
            q_tanz = mrtr%Q(i,j,4,1) - qnorm*mrtr%n_hat(3,i,j)

            mrtr%Q(i,j,1,2) = mrtr%Q(i,j,1,1)
            mrtr%Q(i,j,2,2) = q_tanx - qnorm*mrtr%n_hat(1,i,j)
            mrtr%Q(i,j,3,2) = q_tany - qnorm*mrtr%n_hat(2,i,j)
            mrtr%Q(i,j,4,2) = q_tanz - qnorm*mrtr%n_hat(3,i,j)

            rho = mrtr%Q(i,j,1,2)
            u   = mrtr%Q(i,j,2,2)/rho
            v   = mrtr%Q(i,j,3,2)/rho
            w   = mrtr%Q(i,j,4,2)/rho
            mrtr%Q(i,j,5,2) = rho*twall/(mach*mach*gamma*(gamma-1.d0)) + 0.5d0*rho*(u**2 + v**2 + w**2)

         END DO
      END DO
!
!     -------------------------
!     Solve the riemann problem
!     -------------------------
!
      CALL interface_flux(mrtr,flux,property1,property2,0)
!
      RETURN 
      END SUBROUTINE wallbc_isoth
