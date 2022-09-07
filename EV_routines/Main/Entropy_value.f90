!>\file Entropy_value.f90
!! This file contains all subroutines to calculate the entropy value.

!> @brief Calculate the Entropy 'Sh' on Guass-Gauss collocation points
!! and also calculate the entropy norm in the entire domain.
!!
!! The entropy is calculated using the follwing formula:
!!
!! \f$ e = \frac{\rho}{\gamma^2 - \gamma} \times \log (\frac{p}{\rho^{\gamma}}) \f$
      SUBROUTINE Entropy_Value(d,ngrid,ent_nm_d)    
!      
      USE domain_definition
      USE physics
      USE mpi_par
      USE mpi    
!
      IMPLICIT NONE
!
      TYPE (domain)                          :: d(ngp)
!
      INTEGER                                :: ngrid,id,n,m,l,ntot
      DOUBLE PRECISION                       :: rho,u,v,w,p  
      DOUBLE PRECISION                       :: sum_entropy,ent_nm_d,ent_nm_d_loc,savg,savga
      DOUBLE PRECISION, DIMENSION(ngp)       :: ent_nm
!
!   -------------------------------------------------------------------
!    Calculate the entropy 
!   -------------------------------------------------------------------
      sum_entropy = 0.0d0
      DO id = 1,ngrid
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
               DO n = 1,d(id)%ncg(1)
!
                 rho             =  d(id)%Q(n,m,l,1)*d(id)%jacob(n,m,l)
                 u               =  d(id)%Q(n,m,l,2)/d(id)%Q(n,m,l,1)
                 v               =  d(id)%Q(n,m,l,3)/d(id)%Q(n,m,l,1)
                 w               =  d(id)%Q(n,m,l,4)/d(id)%Q(n,m,l,1)
                 p               = (d(id)%Q(n,m,l,5)*d(id)%jacob(n,m,l)    -  &
                                    0.5d0*rho*(u**2 + v**2 + w**2))*(gamma-1.0d0)
!                                   
                 d(id)%Sh(n,m,l) = (rho/(gamma**2-gamma))*log(p/rho**gamma)
                 sum_entropy     = d(id)%Sh(n,m,l) + sum_entropy 
!                
               END DO
            END DO
         END DO
      END DO 
!   -------------------------------------------------------------------
!    Calculate space-averaged entropy
!   -------------------------------------------------------------------      
      CALL MPI_ALLREDUCE(sum_entropy,savga,1,MPI_DOUBLE_PRECISION, MPI_SUM, &
                          comm1d, ierr)
      CALL MPI_ALLREDUCE(ngrid,ntot,1,MPI_DOUBLE_PRECISION, MPI_SUM, &
                          comm1d, ierr)
    
      IF (myid==0) savg = savga/(ntot*d(1)%ncg(1)*d(1)%ncg(2)*d(1)%ncg(3)) 
      
      CALL MPI_BCAST(savg,1,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
!   -------------------------------------------------------------------
!     Calculate entropy norm on the entire domain
!   -------------------------------------------------------------------
      ent_nm_d_loc = 0.0d0
      DO id = 1,ngrid
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
               DO n = 1,d(id)%ncg(1)
                   ent_nm_d_loc = max(ent_nm_d_loc,abs(d(id)%Sh(n,m,l)-savg))   
               ENDDO
            ENDDO
         ENDDO
      END DO
      CALL MPI_ALLREDUCE(ent_nm_d_loc,ent_nm_d,1,MPI_DOUBLE_PRECISION, MPI_MAX, &
                         comm1d, ierr)
      RETURN
      END SUBROUTINE Entropy_Value

    
!
!> @brief Stores entropy and density values for the past two time steps.
      SUBROUTINE EV_store(time,dt,d,ngrid) 
!                                                               
      USE domain_definition
      USE mpi_par
      USE User_Data
      USE rk_coefs
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)           :: d(ngp)  
      DOUBLE PRECISION, DIMENSION(ngp) :: dt

!   -------------------------------------------------------------------
!     Non-self-starting time schemes, so set S_old to current Sh for the 
!     first two time steps, to ensure that dS/dt and d(rho)/dt are zero. 
!   -------------------------------------------------------------------
!
        IF (time < 2.0*dt(1)) THEN
          DO id = 1,ngrid
             d(id)%Sh_old(:,:,:)   =  d(id)%Sh(:,:,:)
             d(id)%Sh_old2(:,:,:)  =  d(id)%Sh(:,:,:)
             d(id)%rho_old(:,:,:)  =  d(id)%Q(:,:,:,1)*d(id)%jacob(:,:,:)
             d(id)%rho_old2(:,:,:) =  d(id)%Q(:,:,:,1)*d(id)%jacob(:,:,:)
          ENDDO
!
         ELSE
!
          DO id = 1,ngrid
!
             d(id)%Sh_old2  (:,:,:)  = d(id)%Sh_old(:,:,:)
             d(id)%Sh_old   (:,:,:)  = d(id)%Sh    (:,:,:)                
             d(id)%rho_old2 (:,:,:)  = d(id)%rho_old (:,:,:)
             d(id)%rho_old  (:,:,:)  = d(id)%Q(:,:,:,1)*d(id)%jacob(:,:,:)
!
          ENDDO
!
         ENDIF
      RETURN
      END SUBROUTINE EV_store
!
!///////////////////////////////////////////////////////////////////////////////