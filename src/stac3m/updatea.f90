!< \file updatea.f90
!! Contains explicit timestep solution update routines
!///////////////////////////////////////////////////////////////////////////////
!////////								////////
!////////	update.f90						////////
!////////								////////
!////////	contains:						////////
!////////								////////
!////////	   SUBROUTINE update(time,dt,d,kk,resid)		////////
!////////	   SUBROUTINE update_interior(time,dt,d,kk,resid)	////////
!////////								////////
!///////////////////////////////////////////////////////////////////////////////
!                    
!> @brief Call the appropriate routine to update the solution on a grid based on a case switch
      SUBROUTINE update(time,dt,d,kk,resid) 
!                    
!.......................................................................
!     date: 11/5/98 
!     routines called: update_interior
!                      update_w_temporal_damping                                            
!                    
!     Call the appropriate routine to update the solution on a grid                                         
!.......................................................................
!                    
      USE domain_definition
      USE File_Units
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d
      integer       :: id
!
!     added by ksengu1
!      d%domain_type='damping'

      SELECT CASE (d%domain_type)

         CASE ('interior','buffer')
            CALL update_interior(time,dt,d,kk,resid)

         CASE ('damping')
            CALL update_w_temporal_damping(time,dt,d,kk,resid)

         CASE DEFAULT
            WRITE(err_unit,*) 'Unknown domain type: ',d%domain_type
            STOP 'Unknown domain type'
      END SELECT

      RETURN
   END SUBROUTINE update
!                    
!///////////////////////////////////////////////////////////////////////
!                    
!> @brief Update the solution on a grid.
!!
!! Updates the interior point solutions explicitly. Also includes switches for different types of source terms that are
!! specified in the bump file. 
   SUBROUTINE update_interior(time,dt,d,kk,resid) 
!                    
!.......................................................................
!     date: 11/5/98 
!     routines called: dx3dg                                            
!                      dy3dg                                            
!                      dz3dg                                            
!                      compute_source                                            
!
!     applicability: STAC3M hyperbolic equations                             
!                    
!     update the solution on a grid                                         
!     .........................................................................
!                    
      USE domain_definition
      USE physics
      USE rk_coefs
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d
      integer       :: id
!
!     local arrays
!
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qt,source
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: flux_der_X,flux_der_Y,flux_der_Z
      DOUBLE PRECISION, DIMENSION(neq)          :: Q,source_loc
      
      double precision, dimension(4)            :: MW     ! dimensional molecular weights for source term
      double precision, dimension(4)            :: cpT, hT! temporary specific heat and enthalpy variables
      double precision                          :: gammas ! variable ratio of specific heats locally calculated
      double precision                          :: tempT  ! temporary temperature storage variable
      
      LOGICAL                          :: Is_NaN,Is_INF

	  ! Zero out source term to prevent issues
	  source(:,:,:,:) = 0.0d0
!--
!
!     ----------------------------------------------------------------
!     Compute derivatives at each spatial point and compute tendencies
!     ----------------------------------------------------------------
!
      DO nv = 1,neq
         CALL dx3dga(d%f(:,:,:,nv),d%nc,flux_der_X,d%ncg,-0.5d0,d%dmx)
         CALL dy3dga(d%g(:,:,:,nv),d%nc,flux_der_Y,d%ncg,-0.5d0,d%dmy)
         CALL dz3dga(d%h(:,:,:,nv),d%nc,flux_der_Z,d%ncg,-0.5d0,d%dmz)
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2) 
               DO n = 1,d%ncg(1)
                  Qt(n,m,l,nv)    = -(flux_der_X(n,m,l) + &
                                      flux_der_Y(n,m,l) + &
                                      flux_der_Z(n,m,l))
               END DO
            END DO
         END DO
      END DO
!
!     --------------------------
!     Add in source if necessary
!     --------------------------
!
      IF ( source_type == 'uniform' )     THEN

         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               DO n = 1,d%ncg(1)
                  DO nv = 1,neq
                     Q(nv) = d%Q(n,m,l,nv)
                  END DO
                  IF (l==1 .and. m==1 .and. n==1) THEN ! source is not space dependent
                    deltatt = dt*brk(kk)
                    CALL compute_source(d%xg(1,n,m,l),d%xg(2,n,m,l) &
                                     ,d%xg(3,n,m,l),time,source_loc,deltatt)
                  ENDIF
                  DO nv = 1,neq
                     source(n,m,l,nv) = source_loc(nv)/d%jacob(n,m,l)
                  END DO
               END DO
            END DO
         END DO
     end if
     
     if (source_type == 'time_dependent') then
         do l = 1,d%ncg(3)
             do m = 1,d%ncg(2)
                 do n = 1,d%ncg(1)
                     do nv=1,neq
                         Q(nv) = d%Q(n,m,l,nv)
                     end do
                     
                     deltatt = dt*brk(kk)
                     
                     call compute_source( d%xg(1,n,m,l), d%xg(2,n,m,l), &
                                          d%xg(3,n,m,l), time, source_loc, deltatt)
                    do nv=1,neq
                        source(n,m,l,nv) = source_loc(nv)/d%jacob(n,m,l)
                    end do
                end do
            end do
        end do
     end if
    
    if( source_type == 'fmdf' ) then
        ! set some constants
        MW(1) = 2.02d0  / 2.02d0  ! molecular weight of hydrogen
        MW(2) = 32.0d0  / 2.02d0  ! molecular weight of oxygen
        MW(3) = 18.02d0 / 2.02d0  ! molecular weight of nitrogen
        MW(4) = 28.01d0 / 2.02d0  ! molecular weight of products
        
        do l=1,d%ncg(3)
            do m=1,d%ncg(2)
                do n=1,d%ncg(1)
                    do nv=1,neq
                         Q(nv) = d%Q(n,m,l,nv)
                    end do
                    
                    deltatt = dt*brk(kk)
                    
                    ! get the temperature at the gauss point
                    tempT  = d%s(n,m,l,1)
                    
                    ! calculate the ratio of specific heats at the gauss point location
                    call thermodynamics(tempT,d%s(n,m,l,2:5),MW,cpT,hT,gammas)
                    
                    ! bandaiding an error
                    if (kk < 3) gammas = 1.40d0
                    
                    !call fmdfTESource( d, id, n, m, l, dt, source(n,m,l,5))
                    !source(n,m,l,5) = source(n,m,l,5) / d%jacob(n,m,l)
                    !source(n,m,l,5) = d%react(n,m,l)*d%Q(n,m,l,1) / (d%jacob(1,1,1)**2*mach*mach*gamma*(gamma-1.0d0))
                    source(n,m,l,5) = d%react(n,m,l)*d%Q(n,m,l,1) / (gammas*(gammas-1)*mach*mach)
                    !source(n,m,l,5) = d%react(n,m,l)*d%Q(n,m,l,1) / (gamma*(gamma-1)*mach*mach*d%jacob(1,1,1))
                    !write(*,*) source(n,m,l,5)
                end do
            end do
        end do
    end if
         
    if( source_type /= 'none' ) then
         DO nv = 1,neq
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Qt(n,m,l,nv) = Qt(n,m,l,nv) + source(n,m,l,nv)
                  END DO
               END DO
            END DO
         END DO

      END IF
!
!
!     Compute new solution over all interior (gauss) points
!     -----------------------------------------------------
!                    
      resid = 0.0d0
      DO nv = 1,neq
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2) 
               DO n = 1,d%ncg(1)
!                    
                  resid = max(abs(Qt(n,m,l,nv)),resid) 
!                    
                  d%g_Q(n,m,l,nv) = ark(kk)*d%g_Q(n,m,l,nv) + Qt(n,m,l,nv) 
                  d%Q(n,m,l,nv)   = d%Q(n,m,l,nv) + crk(kk)*dt*d%g_Q(n,m,l,nv)
!
               END DO
            END DO
         END DO
      END DO
!                    
      RETURN 
   END SUBROUTINE update_interior
!                    
!///////////////////////////////////////////////////////////////////////
!> @brief Updates the interior solution using temporal damping terms.
!!
!! This subroutines does not include any sourcing terms.             
   SUBROUTINE update_w_temporal_damping(time,dt,d,kk,resid)
!                    
!.......................................................................
!     date: 11/5/98 
!     routines called: dx3dg                                            
!                      dy3dg                                            
!                      dz3dg                                            
!                      compute_source                                            
!
!     applicability: STAC3M hyperbolic equations                             
!                    
!     update the solution on a grid                                         
!     .........................................................................
!                    
      USE domain_definition
      USE User_Data
      USE physics
      USE rk_coefs
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d
!
!     local arrays
!
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qt,source
      DOUBLE PRECISION, DIMENSION(nx,ny,nz)     :: flux_der_X,flux_der_Y,flux_der_Z
      DOUBLE PRECISION, DIMENSION(neq)          :: Q,source_loc
      DOUBLE PRECISION, DIMENSION(nx)          ::  dampx
      DOUBLE PRECISION, DIMENSION(ny)          ::  dampy
      DOUBLE PRECISION, DIMENSION(nz)          ::  dampz
!--
!
!     ----------------------------------------------------------------
!     Compute derivatives at each spatial point and compute tendencies
!     ----------------------------------------------------------------
!
      DO nv = 1,neq
         CALL dx3dga(d%f(:,:,:,nv),d%nc,flux_der_X,d%ncg,-0.5d0,d%dmx)
         CALL dy3dga(d%g(:,:,:,nv),d%nc,flux_der_Y,d%ncg,-0.5d0,d%dmy)
         CALL dz3dga(d%h(:,:,:,nv),d%nc,flux_der_Z,d%ncg,-0.5d0,d%dmz)
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2) 
               DO n = 1,d%ncg(1)
                  Qt(n,m,l,nv)    = -(flux_der_X(n,m,l) + &
                                      flux_der_Y(n,m,l) + &
                                      flux_der_Z(n,m,l))
               END DO
            END DO
         END DO
      END DO
!
!     --------------------------
!     Add in source if necessary
!     --------------------------
!
      IF ( source_type /= 'none' )     THEN

         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               DO n = 1,d%ncg(1)
                  DO nv = 1,neq
                     Q(nv) = d%Q(n,m,l,nv)
                  END DO
                  IF (l==1 .and. m==1 .and. n==1) THEN ! source is not space dependent
                    deltatt = dt*brk(kk)
                    CALL compute_source(d%xg(1,n,m,l),d%xg(2,n,m,l) &
                                     ,d%xg(3,n,m,l),time,source_loc,deltatt)
                  ENDIF
                  DO nv = 1,neq
                     source(n,m,l,nv) = source_loc(nv)/d%jacob(n,m,l)
                  END DO
               END DO
            END DO
         END DO

         DO nv = 1,neq
            DO l = 1,d%ncg(3)
               DO m = 1,d%ncg(2)
                  DO n = 1,d%ncg(1)
                     Qt(n,m,l,nv) = Qt(n,m,l,nv) + source(n,m,l,nv)
                  END DO
               END DO
            END DO
         END DO

      END IF
!
      IF ( d%domain_type == 'damping' )     THEN
!
!        --------------------------------
!        compute temporal damping factors
!        --------------------------------
!
         dampx = 0.0d0
         dampy = 0.0d0
         dampz = 0.0d0
         damping_factor = 10.0d0
         power          = 2.0d0
         IF(d%bcond(4) == 'outflow')     THEN
            DO n = 1,d%ncg(1)
               dampx(n) = damping_factor*d%cxg(n,1)**power
            END DO
         END IF

         nv = 5
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               DO n = 1,d%ncg(1)
                 press= (d%Q(n,m,l,nv) - (d%Q(n,m,l,2)**2+d%Q(n,m,l,3)**2+d%Q(n,m,l,4)**2) &
                     /(d%Q(n,m,l,1)*2.0d0))*(gamma-1.0d0)*d%jacob(n,m,l)
                 Qt(n,m,l,nv) = Qt(n,m,l,nv) - (dampx(n) + dampy(m))*(press-pout)/(gamma-1.0d0)/d%jacob(n,m,l)
               END DO
            END DO
         END DO
!
      END IF

!
!
!     -----------------------------------------------------
!     Compute new solution over all interior (gauss) points
!     -----------------------------------------------------
!                    
      resid = 0.0d0
      DO nv = 1,neq
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2) 
               DO n = 1,d%ncg(1)
!                    
                  resid = max(abs(Qt(n,m,l,nv)),resid) 
!                    
                  d%g_Q(n,m,l,nv) = ark(kk)*d%g_Q(n,m,l,nv) + Qt(n,m,l,nv) 
                  d%Q(n,m,l,nv)   = d%Q(n,m,l,nv) + crk(kk)*dt*d%g_Q(n,m,l,nv)
!
               END DO
            END DO
         END DO
      END DO
!                    
      RETURN 
   END SUBROUTINE update_w_temporal_damping
