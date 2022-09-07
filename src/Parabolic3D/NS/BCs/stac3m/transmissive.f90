!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!
!> @brief Apply the supersonic transmissive boundary condition by fixing
!! the pressure and updating all the values of the outside of the the domain
!! with the values inside the domain. This routine is used to calculate fluxes.
      SUBROUTINE transmissive_boundary_flux(mrtr,flux) 
!
!......................................................................
!
      USE mortar_definition
      USE physics
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property1 = 0.0d0, property2 = 0.0d0
!-----
!
      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
           
           DO nv = 1,neq
              mrtr%Q(i,j,nv,2) = mrtr%Q(i,j,nv,1)
           END DO
           
           p = pout
           rho = mrtr%Q(i,j,1,1)
           u   = mrtr%Q(i,j,2,1)/rho
           v   = mrtr%Q(i,j,3,1)/rho
           w   = mrtr%Q(i,j,4,1)/rho
           mrtr%Q(i,j,5,2) = p/(gamma-1.0d0) + 0.5d0*rho*(u**2 + v**2 + w**2)
           
         END DO
      END DO
      
      CALL interface_flux(mrtr,flux,property1,property2,0.0d0)
!

      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Apply the supersonic transmissive boundary condition by fixing
!! the pressure and updating all the values of the outside of the the domain
!! with the values inside the domain. This routine is used to calculate
!! value of the variables at the boundary.
      SUBROUTINE transmissive_boundary_values(mrtr,Qbound)
!
!......................................................................
!
      USE mortar_definition
      USE physics
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Ql,Qr,Q

!-----
!
!
!     ------------------------------------------------------
!     Rotate and project the solutions onto the mortar space
!     ------------------------------------------------------
!
      IF ( mrtr%orient /= 0 )     CALL rotate_to_mortar(mrtr)
      DO iside = 1,2
         IF ( .NOT.mrtr%conforming(1,iside)   .OR. &
              .NOT.mrtr%conforming(2,iside) )      &
            CALL map_to_mortar_space(iside,mrtr)
      END DO
!
!     ---------------------------------------
!     Solve the Riemann problem on the mortar
!     ---------------------------------------
!
       SELECT CASE ( mrtr%iface(1) )

          CASE (4)
             DO j = 1,mrtr%lenmortar(2)
                DO i = 1,mrtr%lenmortar(1)

                   Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                   Qbound(i,j,2)  = mrtr%Q(i,j,2,1)
                   Qbound(i,j,3)  = mrtr%Q(i,j,3,1)
                   Qbound(i,j,4)  = mrtr%Q(i,j,4,1)
                   
                   p = pout
                   rho = mrtr%Q(i,j,1,1)
                   u   = mrtr%Q(i,j,2,1)/rho
                   v   = mrtr%Q(i,j,3,1)/rho
                   w   = mrtr%Q(i,j,4,1)/rho
                   Qbound(i,j,5) = p/(gamma-1.0d0) + 0.5d0*rho*(u**2 + v**2 + w**2)

                END DO
             END DO
       END SELECT
		
      RETURN
      END
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Apply the supersonic transmissive boundary condition by fixing
!! the pressure and updating all the values of the outside of the the domain
!! with the values inside the domain. This routine is used to calculate
!! viscous fluxes.
      SUBROUTINE transmissive_boundary_vis_flux(mrtr,Qbound)
!
!......................................................................
!     date: 03/14/01
!
!     Compute the inflow/outflow viscous flux from the exterior
!     and interior solutions
!......................................................................
!
      USE mortar_definition
      USE physics
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Q


      IF ( mrtr%orient /= 0 )     CALL rotate_fv_to_mortar(mrtr)

      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
            DO nv = 1,neq
               
               IF (nv==5) THEN

                 p = pout
                 rho = mrtr%Q(i,j,1,1)
                 u   = mrtr%Q(i,j,2,1)/rho
                 v   = mrtr%Q(i,j,3,1)/rho
                 w   = mrtr%Q(i,j,4,1)/rho
                 mrtr%fv(i,j,nv,2) = p/(gamma-1.0d0) + 0.5d0*rho*(u**2 + v**2 + w**2)
               ELSE
                 mrtr%fv(i,j,nv,2) = mrtr%fv(i,j,nv,1)
               END IF
               
               IF (mrtr%nsign == -1) THEN
                  Qs = -mrtr%fv(i,j,nv,2)
               ELSE IF (mrtr%nsign == 1) THEN
                 Qs = mrtr%fv(i,j,nv,2)
               ENDIF
               Qm = mrtr%fv(i,j,nv,1)
               Q(nv)  = 0.5d0*(Qm + Qs)*mrtr%norm(i,j)
            END DO

            DO nv = 1,neq
               Qbound(i,j,nv) = Q(nv)
            END DO

        END DO
      END DO

!
      RETURN
      END
!