!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine uses the mortar master and slave solution to determine the 
!> solution at the interface with the linearized riemann solver.
!!
!! Project the solution at the face values of the Gauss/Lobatto points
      SUBROUTINE boundary_values(mrtr,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     Use the mortar master and slave solution to determine the 
!     solution at the interface with the linearized riemann solver.
!     Project the solution at the face values of the Gauss/Lobatto points
!
!     Subroutines called
!                    riemann_solution
!......................................................................
!
      USE  mortar_definition
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
 
       SELECT CASE(mrtr%iface(1))
         CASE(1,3,6) 
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,2)
                     Qr(nv) = mrtr%Q(i,j,nv,1)
                  END DO

                  CALL riemann_solution(mrtr%n_hat(:,i,j),Ql,Qr,Q)

                  Qbound(i,j,:) = Q(:)
  

               END DO
            END DO


        CASE(2,5)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                  END DO

                  CALL riemann_solution(mrtr%n_hat(:,i,j),Ql,Qr,Q)

                  Qbound(i,j,:) = Q(:)
  

               END DO
            END DO
            
        CASE(4)
        
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  DO nv = 1,neq
!                     mrtr%Q(i,j,nv,2) = mrtr%Q(i,j,nv,1)
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                  END DO

                  CALL riemann_solution(mrtr%n_hat(:,i,j),Ql,Qr,Q)

                  Qbound(i,j,:) = Q(:)
  
               END DO
            END DO
         
      END SELECT


!
      RETURN
      END
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine applies dirichlet wall conditions for  an adiabatic wall 
!> conditions.
!!
!! u = v = z = 0 ==> rhou = rhov = rhoz = 0
      SUBROUTINE set_adiabatic_wall(mrtr,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     apply dirichlet wall conditions for  an adiabatic wall
!     conditions:
!
!     u = v = z = 0 ==> rhou = rhov = rhoz = 0
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound

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

      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  
                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  Qbound(i,j,5)  = mrtr%Q(i,j,5,1)

                  DO nv = 2,neq-1
                    Qbound(i,j,nv) = 0.0d0
                  END DO

               END DO
            END DO

         CASE (2,4,5)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  Qbound(i,j,5)  = mrtr%Q(i,j,5,1)

                  DO nv = 2,neq-1
                    Qbound(i,j,nv) = 0.0d0
                  END DO

               END DO
            END DO
      END SELECT
      RETURN
      END

!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the inflow/outflow viscous flux from the exterior
!> and interior solutions.
!!
!
      SUBROUTINE boundary_vis_flux(mrtr,Qbound)
!
!......................................................................
!     date: 03/14/01
!
!     Compute the inflow/outflow viscous flux from the exterior
!     and interior solutions
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Q


      IF ( mrtr%orient /= 0 )     CALL rotate_fv_to_mortar(mrtr)

      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
            DO nv = 1,neq
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
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the mortar viscous flux from the exterior
!> and interior solutions.
!!
!
      SUBROUTINE mortar_vis_flux(mrtr,Qbound)
!
!......................................................................
!     date: 11/31/01
!
!     Compute the mortar viscous flux from the exterior
!     and interior solutions
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Q


      IF ( mrtr%orient /= 0 )     CALL rotate_fv_to_mortar(mrtr)

      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
            DO nv = 1,neq
               IF (mrtr%nsign == -1) THEN
                  Qs = -mrtr%fv(i,j,nv,2)
               ELSEIF (mrtr%nsign == 1) THEN
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
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sets the normal energy flux to zero, rotate and 
!> project the solutions onto the mortar space.
!!
!
      SUBROUTINE wall_vis_flux(mrtr,Qbound,bcond)
!
!......................................................................
!     date: 03/12/01
!
!     set the normal energy flux to zero
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      CHARACTER(LEN = 9)                         :: bcond



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

      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                      Qbound(i,j,1) = mrtr%fv(i,j,1,1)*mrtr%norm(i,j)
                      Qbound(i,j,2) = mrtr%fv(i,j,2,1)*mrtr%norm(i,j)
                      Qbound(i,j,3) = mrtr%fv(i,j,3,1)*mrtr%norm(i,j)
                      Qbound(i,j,4) = mrtr%fv(i,j,4,1)*mrtr%norm(i,j)
                   
                      IF (bcond == 'walladiab') THEN
                        Qbound(i,j,5) = 0.0d0
                      ELSE
                        Qbound(i,j,5) = mrtr%fv(i,j,5,1)*mrtr%norm(i,j)
                      END IF
               END DO
            END DO

         CASE (2,4,5)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                      Qbound(i,j,1) = mrtr%fv(i,j,1,1)*mrtr%norm(i,j)
                      Qbound(i,j,2) = mrtr%fv(i,j,2,1)*mrtr%norm(i,j)
                      Qbound(i,j,3) = mrtr%fv(i,j,3,1)*mrtr%norm(i,j)
                      Qbound(i,j,4) = mrtr%fv(i,j,4,1)*mrtr%norm(i,j)

                      IF (bcond == 'walladiab') THEN
                        Qbound(i,j,5) = 0.0d0
                      ELSE
                        Qbound(i,j,5) = mrtr%fv(i,j,5,1)*mrtr%norm(i,j)
                      END IF

               END DO
            END DO
      END SELECT
!
      RETURN
      END
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine applies dirichlet wall conditions for  an adiabatic 
!> wall conditions.
!!
!! u = v = z = 0 ==> rhou = rhov = rhoz = 0, T = twall
      SUBROUTINE set_isothermal_wall(mrtr,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     apply dirichlet wall conditions for  an adiabatic wall
!     conditions:
!
!     u = v = z = 0 ==> rhou = rhov = rhoz = 0
!     T = twall
!......................................................................
!
      USE  mortar_definition
      USE  physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound

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

      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  DO nv = 2,neq-1
                    Qbound(i,j,nv) = 0.0d0
                  END DO
                  Qbound(i,j,5) = mrtr%Q(i,j,1,1)*twall/(mach*mach*gamma*(gamma-1.d0)) 

               END DO
            END DO

         CASE (2,4,5)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  DO nv = 2,neq-1
                    Qbound(i,j,nv) = 0.0d0
                  END DO
                  Qbound(i,j,5) = mrtr%Q(i,j,1,1)*twall/(mach*mach*gamma*(gamma-1.d0)) 

               END DO
            END DO
      END SELECT
      RETURN
      END

!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine applies dirichlet wall conditions for  an isothermal
!>  sliding wall.
!!
!! u = upper_velocity  ==> rhou = rho*upper_velocity
!! v = z = 0 ==> rhov = rhoz = 0
!! T = twall
!
      SUBROUTINE set_sliding_wall(mrtr,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     apply dirichlet wall conditions for  an isothermal
!     sliding wall
!
!     u = upper_velocity  ==> rhou = rho*upper_velocity
!     v = z = 0 ==> rhov = rhoz = 0
!     T = twall
!......................................................................
!
      USE  mortar_definition
      USE  physics
      USE  User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound

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

      SELECT CASE ( mrtr%iface(1) )

         CASE (1,3,6)

            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  DO nv = 2,neq-1
                    Qbound(i,j,nv) = 0.0d0
                  END DO
                  Qbound(i,j,3) = Qbound(i,j,1)*upper_velocity
                  Qbound(i,j,5)  = mrtr%Q(i,j,1,1)*twall/(mach*mach*gamma*(gamma-1.d0)) 

               END DO
            END DO

         CASE (2,4,5)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  DO nv = 2,neq-1
                    Qbound(i,j,nv) = 0.0d0
                  END DO
                  Qbound(i,j,2) = Qbound(i,j,1)*mrtr%Q(i,j,2,2)
                  !Qbound(i,j,2) = Qbound(i,j,1)*0.0d0
                  Qbound(i,j,5) = Qbound(i,j,1)*twall/(mach*mach*gamma*(gamma-1.d0)) + 0.5d0*(Qbound(i,j,2)**2)/Qbound(i,j,1)

               END DO
            END DO
      END SELECT
      RETURN
      END

!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sets the periodic boundary condition.
!!
!
      SUBROUTINE set_periodm(mrtr,Qbound)
!
!......................................................................
!     date: 12/18/01
!     
!......................................................................
!     
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Ql,Qr,Q
!                    
      !CALL V_Avg_Soln(mrtr,Qbound)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,2)
                     Qr(nv) = mrtr%Q(i,j,nv,1)
                  END DO

                  CALL riemann_solution(mrtr%n_hat(:,i,j),Ql,Qr,Q)

                  Qbound(i,j,:) = Q(:)
 
      
               END DO
            END DO
      
      RETURN 
      END
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sets the periodic boundary condition.
!!
!           
!              
      SUBROUTINE set_periods(mrtr,Qbound)
!                    
!......................................................................
!     date: 12/18/01
!
!......................................................................
!
      USE  mortar_definition
      USE  physics
      USE  input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      DOUBLE PRECISION, DIMENSION(neq)           :: Ql,Qr,Q
!

      !CALL V_Avg_Soln(mrtr,Qbound)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)
                  DO nv = 1,neq
                     Ql(nv) = mrtr%Q(i,j,nv,1)
                     Qr(nv) = mrtr%Q(i,j,nv,2)
                  END DO

                  CALL riemann_solution(mrtr%n_hat(:,i,j),Ql,Qr,Q)

                  Qbound(i,j,:) = Q(:)


               END DO
            END DO


      DO j = 1,mrtr%lenmortar(2)
          DO i = 1,mrtr%lenmortar(1)
             rho  = Qbound(i,j,1)
             rhou = Qbound(i,j,2)
             rhov = Qbound(i,j,3)
             rhow = Qbound(i,j,4)
             rhoe = Qbound(i,j,5)

             u = rhou/rho
             v = rhov/rho
             w = rhow/rho
             pold = (gamma-1.d0)*(rhoe - 0.5d0*rho*(u**2 + v**2 + w**2))
             dpdx    = -3.0d0/re
             dpdx    = dpdxturb
             pnew    = pold + dpdx*char_length
             rhonew  = pnew*rho/pold

             Qbound(i,j,1)  = rhonew
             Qbound(i,j,2)  = rhonew*u
             Qbound(i,j,3)  = rhonew*v
             Qbound(i,j,4)  = rhonew*w
             Qbound(i,j,5)  = pnew/(gamma - 1.d0) +  0.5d0*rhonew*(u**2 + v**2 + w**2)
!            write(*,*) 'periods',Qbound(i,j,:)
          END DO
      END DO
!
      RETURN
      END
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine applies dirichlet wall conditions for  an adiabatic wal
!> conditions.
!!
!! u = v = z = 0 ==> rhou = rhov = rhoz = 0
      SUBROUTINE Symm_Soln(mrtr,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     apply dirichlet wall conditions for  an adiabatic wall
!     conditions:
!
!     u = v = z = 0 ==> rhou = rhov = rhoz = 0
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound

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

!       SELECT CASE ( mrtr%iface(1) )


!          CASE (2)
            DO j = 1,mrtr%lenmortar(2)
               DO i = 1,mrtr%lenmortar(1)

                  Qbound(i,j,1)  = mrtr%Q(i,j,1,1)
                  Qbound(i,j,2)  = mrtr%Q(i,j,2,1)
                  Qbound(i,j,3)  = -mrtr%Q(i,j,3,1)
                  Qbound(i,j,4)  = mrtr%Q(i,j,4,1)
                  Qbound(i,j,5)  = mrtr%Q(i,j,5,1)


               END DO
            END DO
!       END SELECT
      RETURN
      END

!
!///////////////////////////////////////////////////////////////////////
