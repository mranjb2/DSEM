!> \file part_integrate_adbash.f90
!! Particle integration routines using Euler and Adams-Bashforth
!//////////////////////////////////////////////////////////////////////////////
!////////                                                               ////////
!////////       part_int.f90                                            ////////
!////////                                                               ////////
!////////       contains:                                               ////////
!////////                                                               ////////
!////////          SUBROUTINE part_integrate(k,dt,d,drop,////////
!////////				     ngrid,nprt)                ////////
!////////          SUBROUTINE integr_onepart_euler(dt,drop)             ////////
!////////                                                               ////////
!///////////////////////////////////////////////////////////////////////////////
!> Particle integration routine
   SUBROUTINE part_integrate(k,dt,d,drop,ngrid,nprt,time)
!
!     date: 01/14/00
!     
!     contains subroutines;
!                         lin_Pinterp
!                         integr_onepart
!
!     applicability: mortar versions, euler/navier-stokes
!
!     integrate the particles for one time step
!
!
!
      USE size
      USE domain_definition
      USE particle_definition
      USE part_par
      USE input_data
      USE mpi_par

!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain),DIMENSION(ngp)     :: d
      TYPE(particle),DIMENSION(npart)  :: drop
      DOUBLE PRECISION,DIMENSION(ngp)  :: dt

      LOGICAL          :: dombound
      INTEGER          :: nsurr,msurr,ngrid,pgrid
!
!     --------------------------------------------------------------------------
!     interpolate the solution properties to properties at the particle position
!     and integrate the particle in time
!     the following boundary condition is applied:
!     if the particle is outside the domain then stop integrating
!     --------------------------------------------------------------------------
!
         DO ip=1,nprt
  	    IF (drop(ip)%onoff == 1) THEN
!
 	       pgrid = drop(ip)%ngrid
               !IF (d(pgrid)%ncg(1) < 7 ) THEN !!!!NOTE THIS SHOULD BE 7
                 CALL spec_Pinterp(d(pgrid),drop(ip),myid,ip)
               !ELSE
               !  CALL Lagr_Pinterp(d(pgrid),drop(ip),myid,ip)
               !END IF
!
!              drop(ip)%Vp = drop(ip)%Vfp ! for particle tracers
               !IF (k == 1) THEN
               !   CALL integr_onepart_euler(dt,drop(ip))
               !ELSEIF (k > 1) THEN
                !  nrms = 0
                !  CALL integr_onepart_adbash(dt,drop(ip),time,k,nrms)
               !END IF
               call integr_onepart_sde(dt,drop(ip))
!
	       CALL part_map(d(pgrid),drop(ip))
!
               IF (drop(ip)%Xpmap(1) > 1.0d0 .or. drop(ip)%Xpmap(1)<0.0d0 .or.      & 
                   drop(ip)%Xpmap(2) > 1.0d0 .or. drop(ip)%Xpmap(2)<0.0d0 .or.     &
                   drop(ip)%Xpmap(3) > 1.0d0 .or. drop(ip)%Xpmap(3)<0.0d0 ) THEN 
                   CALL part_bc(d,drop,ip,ngrid)
               END IF
            END IF
         ENDDO
        
      RETURN
   END SUBROUTINE part_integrate
!
!///////////////////////////////////////////////////////////////////////////////
!
!> Integration using Euler scheme
       SUBROUTINE integr_onepart_euler(dt,dr)
!
!      this subroutine integrates one particle forward in time using a 
!      first order Euler scheme
!
  
      USE size
      USE particle_definition
      USE constants
      USE part_par
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      TYPE(particle)  :: dr

      DOUBLE PRECISION, DIMENSION(ngp) :: dt
     
!
!     define right-hand-side of particle equations
!

      rhs1 = dr%Vp(1)
      rhs2 = dr%Vp(2)
      rhs3 = dr%Vp(3)
      rhs4 = dr%Vfp(1)
      rhs5 = dr%Vfp(2)
      rhs6 = dr%Vfp(3)
      rhs7 = dr%Tfp
      rhs8 = -f4*yfps*sqrt(taup)

!
!     integration of right-hand side
!

      dr%Xp(1) 	= dr%Xp(1) + dt(1)*rhs1 
      dr%Xp(2) 	= dr%Xp(2) + dt(1)*rhs2
      dr%Xp(3) 	= dr%Xp(3) + dt(1)*rhs3 
      dr%Vp(1) 	= dr%Vfp(1) 
      dr%Vp(2) 	= dr%Vfp(2) 
      dr%Vp(3) 	= dr%Vfp(3) 
      
      dr%Tp     = dr%Tfp


      dr%Xpnm(1)  = rhs1
      dr%Xpnm(2)  = rhs2
      dr%Xpnm(3)  = rhs3
      dr%Vp1nm(1) = rhs4
      dr%Vp1nm(2) = rhs5
      dr%Vp1nm(3) = rhs6
      dr%Tpnm     = rhs7
      dr%Mpnm     = rhs8

      RETURN
   END SUBROUTINE integr_onepart_euler 

   !> @brief Integration of particles using Euler scheme
   !! 
   !! Integrates the FMDF stocastic differential equations in addition of particle motion
    SUBROUTINE integr_onepart_sde(dt,p)
!
!      this subroutine integrates one particle forward in time using a 
!      first order Euler scheme
!
        use random
        USE size
        USE particle_definition
        USE constants
        USE part_par
        USE physics
!
        implicit none

        type(particle)                      :: p
        double precision, dimension(ngp)    :: dt
        double precision, dimension(3)      :: velo,drift,diff
        double precision                    :: diffpart,ran1,ran2,ran3, constdr
        
        !p%Vp=p%Vfp
        
        !   Calculate the convective portion
        velo(1) = p%Xp(1) + dt(1)*p%Vp(1)
        velo(2) = p%Xp(2) + dt(1)*p%Vp(2)
        velo(3) = p%Xp(3) + dt(1)*p%Vp(3)
        
        !   Calculate the drift portion
        constdr  = (dt(1)/(re*p%rhofp*Sc))
        !constdr  = dt(1)/p%rhofp
        drift(1) = constdr*p%nu_d(1)
        drift(2) = constdr*p%nu_d(2)
        drift(3) = constdr*p%nu_d(3)
        
        !   Calculate the diffusion portion
        ran1        = random_normal()/100.0
        ran2        = random_normal()/100.0
        ran3        = random_normal()/100.0
!        diffpart    = dsqrt(2.0d0*dt(1)/re)*dsqrt( (1.0d0/(re*Sc)) + p%nu_t/Sc_t )
!        write(*,*) 'diffpart',diffpart
!        write(*,*)'dt',dt(1),'re',re,'Sc',Sc,'Sc_t',Sc_t,'nu_t',p%nu_t
        diffpart = dsqrt(2.0d0*dt(1)*((1.0d0/(re*Sc) + (p%nu_t)/Sc_t)))
!        diffpart = dsqrt(2.0d0*dt(1)*((p%nu_t)/(re*Sc_t)))/20.0d0
        diff(1)     = diffpart*ran1
        diff(2)     = diffpart*ran2
        diff(3)     = diffpart*ran3
        
        !   Calculate new particle position
        p%Xp(1) = velo(1) !+ drift(1) !+ diff(1)
        p%Xp(2) = velo(2) !+ drift(2) !+ diff(2)
        p%Xp(3) = velo(3) !+ drift(3)! + diff(3)
        
        !   store stuff
        p%Vp = p%Vfp
  END SUBROUTINE integr_onepart_sde
