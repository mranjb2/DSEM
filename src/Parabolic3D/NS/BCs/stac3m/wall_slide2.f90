!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE wallbc_slide(mrtr,flux) 
!
!......................................................................
!     date: 03/14/01
!
!     Compute the wall bc using characteristic splitting,                 
!     by specifying an equal but opposite flow as the external condition.
!     Compute Q(5) using twall
!     Compute the wall flux: For a *moving* isothermal wall
!     velocity of moving wall is assumed to be in User_Data
!     as the variable "upper_velocity"
!......................................................................
!
      USE  mortar_definition
      USE  physics
      USE  User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property1 = 0.0d0, property2 = 0.0d0
      DOUBLE PRECISION, DIMENSION(nmax,nmax)     :: uu
      DOUBLE PRECISION, DIMENSION(neq)           :: Ql, Qr, f
!
!
!     -----------------------------------------
!     Generate the external flow along the face
!     -----------------------------------------
!
      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
           uu(i,j) = mrtr%Q(i,j,2,2)
         ENDDO
      ENDDO

      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)

            qnorm  = (mrtr%n_hat(1,i,j)*mrtr%Q(i,j,2,1) + &   ! normal momentum
                      mrtr%n_hat(2,i,j)*mrtr%Q(i,j,3,1) + &
                      mrtr%n_hat(3,i,j)*mrtr%Q(i,j,4,1))
            !q_tanx = mrtr%Q(i,j,2,1) - qnorm*mrtr%n_hat(1,i,j)
            q_tanx =  mrtr%Q(i,j,1,1)*uu(i,j)
            !q_tanx = mrtr%Q(i,j,1,1)*0.0d0
            !q_tany = mrtr%Q(i,j,1,1)*upper_velocity
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
!     -----------------------------------------------
!     Solve the riemann problem using an osher solver
!     -----------------------------------------------
!
      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                 DO nv = 1,neq
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                  END DO

                  CALL osher(-mrtr%n_hat(:,i,j),-mrtr%norm(i,j),Ql,Qr,f)

                  DO nv = 1,neq
                     flux(i,j,nv) = f(nv)
                  END DO
	       END DO
            END DO 

         CASE (2,4,5)

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                 DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                  END DO

                  CALL osher(mrtr%n_hat(:,i,j),mrtr%norm(i,j),Ql,Qr,f)

                  DO nv = 1,neq
                     flux(i,j,nv) = f(nv)
                  END DO
	       END DO
            END DO 

      END SELECT
!
      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
           mrtr%Q(i,j,2,2) =  uu(i,j)
         ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE wallbc_slide
