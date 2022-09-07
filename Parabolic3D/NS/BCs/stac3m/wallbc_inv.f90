!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the wall bc using characteristic splitting,
!> by specifying an equal but opposite flow as the external condition.
!!
!! This bc is equivalent to an interface mortar condition
!! where the neighbor element has the exterior flow.
!
      SUBROUTINE wallbc_Inv(mrtr,flux) 
!
!......................................................................
!     date: 12/1/98
!
!     Compute the wall bc using characteristic splitting,                 
!     by specifying an equal but opposite flow as the external condition.
!     This bc is equivalent to an interface mortar condition
!     where the neighbor element has the exterior flow.
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property1 = 0.0d0, property2 = 0.0d0
!!
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
            mrtr%Q(i,j,5,2) = mrtr%Q(i,j,5,1)

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
      END SUBROUTINE wallbc_Inv
