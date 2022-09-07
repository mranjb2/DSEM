!> \file NS_data.f90
!! Navier-Stokes physics modules
!//////////////////////////////////////////////////////////////////////////////////////
!> Contain Navier-Stokes physics values
      MODULE physics
!
         DOUBLE PRECISION  :: gamma                   !< the gas gamma (nonlinear)
         DOUBLE PRECISION  :: sqrtgam                 !< sqrt(gamma)
         DOUBLE PRECISION  :: mach                    !< non-dimensional mach number
         DOUBLE PRECISION  :: machfluid               !< free-stream mach number
         DOUBLE PRECISION  :: pr                      !< prandtl number
         DOUBLE PRECISION  :: re                      !< reynolds number
         double precision  :: Sc_t                    !< turbulent Schmidt number
         double precision  :: Dam                     !< Damkohler number
         double precision  :: Ce                      !< Rate of reaction term
         double precision  :: Ze                      !< Zeldovich number
         double precision  :: c_omega                 !< FMDF sgs constant
         double precision  :: constOmegaLeft          !< Left hand side of Omega_m
         double precision  :: constOmegaRight         !< Right hand side of Omega_m
         double precision  :: Sc                      !< molecular Schmidt number
         DOUBLE PRECISION  :: twall                   !< wall temperature
         DOUBLE PRECISION  :: dpdxturb                !< turbulent pressure gradient
         DOUBLE PRECISION  :: C_s                     !< Smagorinsky Turbulent Coefficient
         DOUBLE PRECISION  :: numax                   !< Maximum Turbulent [kinematic] Viscosity
         DOUBLE PRECISION  :: Cmu                     !< EV method coefficient
         DOUBLE PRECISION  :: mumax                   !< Maximum Entropy [dynamic] Viscosity
         CHARACTER (LEN = 14) :: source_type = 'none' !< source terms in the equations (none,fixed,time_dependent)
!
!
!
!    Particle Physics
!
         DOUBLE PRECISION  :: RHOP                    !< particle density
         DOUBLE PRECISION  :: TBOIL                   !< boiling temperature
         DOUBLE PRECISION  :: PBOIL                   !< boiling pressure
         DOUBLE PRECISION  :: TAUP0                   !< particle relaxation time
         DOUBLE PRECISION  :: A1                      !< normalized latent heat of evaporation
         DOUBLE PRECISION  :: A2                      !< specific heat of liquid/specific heat of particle

      END MODULE physics
!
!///////////////////////////////////////////////////////////////////////
!
!
!>    @brief Store equation dependent keywords needed to read the input file.
!!     Any number of physics keywords can be added to be accessible from
!!     the data module "physics"
!!
!!     Communication is through the data defined in the module "physics"
!!     The values in physics are set below in the routine set_physics_values()
   MODULE physics_keywords
      SAVE
!
      INTEGER, PARAMETER                   :: mxpk = 21 !< maximum physics keywords
      INTEGER                              :: num_physics_keywords = mxpk
      CHARACTER (LEN=10), DIMENSION(mxpk)  :: physics_keyword_list = (/ &
                                                'mach      ',&
                                                'maf       ',&
                                                'reynolds  ',&
                                                'wall      ',&
                                                'prandtl   ',&
                                                'RHOP      ',&
                                                'TBOIL     ',&
                                                'PBOIL     ',&
                                                'TAUP0     ',&
                                                'A1        ',&
                                                'A2        ',&
                                                'Sc        ',&
                                                'sgsSh     ',&
                                                'COmega    ',&
                                                'Da        ',&
                                                'Ze        ',&
                                                'Ce        ',&
                                                'Cs        ',&
                                                'numax     ',&
                                                'Cmu       ',&
                                                'mumax     ' &
                                           /)
      LOGICAL, DIMENSION(mxpk) :: physics_keyword_set = .false.
   END MODULE physics_keywords
!
!///////////////////////////////////////////////////////////////////////
!
!> Reads in the input data-file and saves the values in physics
   SUBROUTINE set_physics_values(the_word,input_line,kword)
      USE physics
      USE physics_keywords
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION      :: get_dp_value   !< extracts a double prec number
      INTEGER               :: get_int_value  !< extracts an integer
      CHARACTER (LEN = 132) :: input_line
      CHARACTER (LEN = 10)  :: the_word
!
      SELECT CASE(the_word)

         CASE('mach')
            mach = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('maf')
            machfluid = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('wall')
            twall = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('reynolds')
            re = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('prandtl')
            pr = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('RHOP')
            RHOP = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('TBOIL')
            TBOIL = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('PBOIL')
            PBOIL = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('TAUP0')
            TAUP0 = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('A1')
            A1 = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         CASE('A2')
            A2 = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         case('sgsSh')
            Sc_t = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

         case('Sc')
            Sc = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.
        
         case('COmega')
            c_omega = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.
            
        case('Da')
            Dam = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.
        
        case('Ce')
            Ce = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.
            
        case('Ze')
            Ze = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

        case('Cs')
            C_s = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

        case('numax')
            numax = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

        case('Cmu')
            Cmu = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

        case('mumax')
            mumax = get_dp_value(input_line)
            physics_keyword_set(kword) = .true.

      END SELECT
!
      RETURN
   END SUBROUTINE set_physics_values
!
!//////////////////////////////////////////////////////////////////////////////////////
!>     this routine will set up equation dependent physical parameters
!!     in this case, set the gass gamma and its square root
   SUBROUTINE set_up_physics
!
      USE physics

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      gamma     = 1.4d0
      sqrtgam   = sqrt(gamma)
      dpdxturb  = -3.0d0/re

      RETURN
   END SUBROUTINE set_up_physics
