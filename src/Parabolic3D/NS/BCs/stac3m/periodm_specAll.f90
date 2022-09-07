!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes periodic bc with pressure jump, solve the 
!> riemann problem with the linearized riemann problem
!!
!
      SUBROUTINE periodm_specAll(mrtr,flux) 
!
!......................................................................
!     date: 12/18/01
!
!     compute periodic bc with pressure jump
!......................................................................
!
      USE  mortar_definition
      USE  physics
      USE  input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(neq)           :: Ql,Qr,f
 
!-----
!
!
      DO j = 1,mrtr%lenmortar(2)
          DO i = 1,mrtr%lenmortar(1)
             rho  = mrtr%Q(i,j,1,2)
             rhou = mrtr%Q(i,j,2,2)
             rhov = mrtr%Q(i,j,3,2)
             rhow = mrtr%Q(i,j,4,2)
             rhoe = mrtr%Q(i,j,5,2)

             u = rhou/rho
             v = rhov/rho
             w = rhow/rho
             pold = (gamma-1.d0)*(rhoe - 0.5d0*rho*(u**2 + v**2 + w**2))
             dpdx = -3.0d0/re
             dpdx = dpdxturb
             pnew = pold - dpdx*char_length
             rhonew = rho*pnew/pold
            
             mrtr%Q(i,j,1,2)  = rhonew 
             mrtr%Q(i,j,2,2)  = rhonew*u
             mrtr%Q(i,j,3,2)  = rhonew*v
             mrtr%Q(i,j,4,2)  = rhonew*w
             mrtr%Q(i,j,5,2)  = pnew/(gamma - 1.d0) +  0.5d0*rhonew*(u**2 + v**2 + w**2)
          END DO
      END DO
 
!   Solve the riemann problem with the linearized riemann problem

      DO j = 1,mrtr%lenmortar(2)
         DO i = 1,mrtr%lenmortar(1)
            DO nv = 1,neq
               Qr(nv) = mrtr%Q(i,j,nv,1)
               Ql(nv) = mrtr%Q(i,j,nv,2)
            END DO

            CALL riemann_alt(mrtr%n_hat(:,i,j),mrtr%norm(i,j),Ql,Qr,f)

            DO nv = 1,neq
               flux(i,j,nv) = f(nv)
            END DO
         END DO
      END DO

!
      RETURN 
      END                                           
