!>\file Driver.f90
!! This file is the driver for for the entropy viscosity method. It is being called in the
!! tstep subroutine and calculates and stores the entropy viscosity.

!/////////////////////////////////////////////////////////////////////////////////////////
!                                                                                 ////////
! SUBROUTINE   Entropy(time,dt,d,ngrid,nmort,mrtr)                                ////////
!                                                                                 ////////
!    Main routine for Computation of Entropy viscosity according                  ////////
!    to Guermond et al. 2011                                                      ////////
!                                                                                 ////////
!      Contains:                                                                  ////////
!          SUBROUTINE Entropy_Value     (time,d,ngrid,dt)                         ////////
!          SUBROUTINE Entropy_Fluxes    (time,d,nmort,mrtr,ngrid)                 ////////
!          SUBROUTINE Entropy_Viscosity (d,ngrid,dt)                              ////////
!          SUBROUTINE EV_Prolong        (d,ngrid)                                 ////////
!                                                                                 ////////
!//////////////////////////////////////////////////////////////////////////////////////////
!                                                                       
!> @brief Main routine for Computation of Entropy viscosity according to Guermond et al. (2011).
      SUBROUTINE Entropy(time,dt,d,ngrid,nmort,mrtr)

      USE domain_definition
      USE mortar_definition
      USE mpi_par
      USE mpi
      USE input_data
!
      IMPLICIT NONE 

      TYPE (domain)                   :: d(ngp)         !< domain
      TYPE (mortar)                   :: mrtr(nmp)      !< mortar
!
      INTEGER                         :: ngrid    !< number of elements
      INTEGER                         :: nmort    !< number of mortars
      DOUBLE PRECISION                :: time     !< time
      DOUBLE PRECISION                :: ent_nm_d !< entropy norm
      DOUBLE PRECISION, DIMENSION(ngp):: dt       !< time step size

!     ---------------------------------------
!     Calculate the EV residual
!     ---------------------------------------
      CALL RHS_ent_trans     (d,ngrid)  
!     ---------------------------------------
!     Construct EV on Gauss points
!     ---------------------------------------
      CALL Strain            (d,ngrid)
      CALL Entropy_Viscosity (d,ngrid)
!     ---------------------------------------
!     Perform smoothing on EV
!     ---------------------------------------      
      IF (Smoothing) CALL entropy_viscosity_smooting(d,ngrid)
!     ---------------------------------------
!     Prolong EV from Gauus to Labatto points
!     ---------------------------------------     
      CALL EV_Prolong        (d,ngrid)  
!            
      RETURN
!
      END SUBROUTINE Entropy
!
!/////////////////////////////////////////////////////////////////////////////////////////