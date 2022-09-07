!> \file particleOperations.f90
!! Contains particle routines pertaining to FMDF Monte Carlo particles

!> @brief This subroutine counts the number of particles that belong to an ensemble element
!!
!! It then assigns the particle to a list of pointers which belong to the ensemble.
!! This prevents the need to search multiple times.
subroutine countParticles(d,p,eP,ngridl,nprt)
    use domain_definition
    use particle_definition
    use size
    use mpi
    implicit none
    type(domain), dimension(ngp)            :: d
    type(particle), dimension(npart),target :: p
    type(ePointers), dimension(ngp,npart) :: eP
    integer, dimension(ngp)                 :: count
    integer, intent(in)                     :: ngridl, nprt
    integer                                 :: np, id
    
    ! Clear out the counter array
    count(:) = 0
    
    write(*,*) 'in counting nprt',nprt,'npart',npart
    
    ! Loop over all the particles
    do np=1,nprt
        if( p(np)%onoff == 1 ) then ! Check if particle is live
            id = p(np)%ngrid
            
            ! create a pointer to the particle
            count(id) = count(id) + 1
            eP(id,count(id))%p => p(np) 
        end if
    end do
    
    ! Tell each ensemble how many particles it has
    do id=1,ngridl
        d(id)%np = count(id)
    end do
    
end subroutine countParticles

!> @brief Initially calculates the weighting of the particles seeded into the domain
!!
!! \f[
!!      w = \frac{\rho^{\alpha} V}{N_p}
!! \f]
!! Where \f$ \rho^{\alpha} \f$ is the density assigned to a particle, \f$ V \f$ 
!! is the volume of the element (in our case the Jacobian), and \f$ N_p \f$ is the 
!! quantity of particles within the element
subroutine calcWeight(d,eP,id)
    use size
    use particle_definition
    use domain_definition
    implicit none
    type(domain)                                :: d
    type(ePointers), dimension(ngp,npart)       :: eP
    double precision                            :: w, vol
    integer                                     :: np, id
    
    ! Set the volume of the element to the jacobian
    vol = 1.0d0/d%jacob(1,1,1)
    !write(*,*) 'vol',vol
    
    !write(*,*) 'number of ensemble particles',e%np
    
    ! calculate the weight and assign it to the particle
    do np=1,d%np
        !write(*,*) 'rhofp',eP(id,np)%p%Rhofp
        w = eP(id,np)%p%Rhofp * vol / d%np
        !write(*,*) 'd%np',d%np
        eP(id,np)%p%w = w
        !write(*,*) 'w',w
    end do
    
end subroutine calcWeight

!!> This subroutine initializes the scalars based on their intiial condition and properties. Must be called prior to initiallizing FMDF
! subroutine initScalars(d,eP,id)
!     use size
!     use particle_definition
!     use domain_definition
!     implicit none
!     type(domain)                                :: d
!     type(ePointers), dimension(ngp,npart)       :: eP
!     double precision                            :: x,y,z
!     integer                                     :: np, id
!     double precision, parameter                 :: pie        = 3.14159265358979323846264338327d0
!     double precision                            :: delta
!     double precision                            :: T,F,O,N,P,spike
!
!
!     do np=1,d%np
!         x = eP(id,np)%p%Xp(1)
!         y = eP(id,np)%p%Xp(2)
!         z = eP(id,np)%p%Xp(3)
!
!         !For CV Reactor
!         !eP(id,np)%p%scalar(1) = eP(id,np)%p%T
!         !eP(id,np)%p%scalar(2) = 0.01152d0
!         !eP(id,np)%p%scalar(3) = 0.23042d0
!         !eP(id,np)%p%scalar(4) = 0.0d0
!         !ep(id,np)%p%scalar(5) = 0.75806d0
!
! !         !for CV reactor
! !         F     = 0.01152d0
! !         O     = 0.23042d0
! !         N     = 0.75806d0
! !         P   = 0.0d0
!         !for ramp cavity with injector
!         F     = 0.000001d0
!         O     = 0.24194d0
!         N     = 0.75806d0
!         P     = 0.0d0
!
!         !for Sod Problem
!         !if (x > 0.50d0) then
!         !    F     = 0.0115d0
!         !    O     = 0.2304d0
!         !    N     = 0.7581d0
!         !    P   = 0.0d0
!         !else
!         !    F = 0.0d0
!         !    O = 0.23310d0
!         !    N = 0.76690d0
!         !    P = 0.0d0
!         !end if
!
!         eP(id,np)%p%scalar(1) = eP(id,np)%p%T
!         eP(id,np)%p%scalar(2) = F
!         eP(id,np)%p%scalar(3) = O
!         eP(id,np)%p%scalar(4) = P
!         eP(id,np)%p%scalar(5) = N
!         eP(id,np)%p%react = 0.0d0
!
!         ! for ramp cavity
!         !T = eP(id,np)%p%Tfp
!         !
!         !if(x.le.3) then
!         !    !T = 1.0d0*(spike-(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
!         !    F = 0.0d0
!         !    O = 0.2331d0
!         !    N = 0.7669d0
!         !    P = 1.0d0-F-O-N
!         !else if(x.ge.6) then
!         !    F = 0.0d0
!         !    O = 0.2331d0
!         !    N = 0.7669d0
!         !    P = 1.0d0-F-O-N
!         !!else if(y.le.dsqrt(1.50d0**2.0-(x-4.50d0)**2.0)+0.650d0) then
!         !else if(y.le.0.650d0) then
!         !    !T = 1.0d0*(spike+(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
!         !    F = 0.0115d0
!         !    O = 0.2304d0
!         !    N = 0.7581d0
!         !    P = 1.0d0-F-O-N
!         !else
!         !    F = 0.0d0
!         !    O = 0.2331d0
!         !    N = 0.7669d0
!         !    P = 1.0d0-F-O-N
!         !end if
!         !
!         !eP(id,np)%p%scalar(1) = eP(id,np)%p%T
!         !eP(id,np)%p%scalar(2) = F
!         !eP(id,np)%p%scalar(3) = O
!         !eP(id,np)%p%scalar(4) = P
!         !eP(id,np)%p%scalar(5) = N
!         ! end of for cavity
!
!         !eP(id,np)%p%scalar(1) = 3.0
!         !eP(id,np)%p%scalar(2) = 4.0
!
!         !if (y .ge. 0.0d0) then
!         !    ep(id,np)%p%scalar(1) = erf(3.14159d0**0.5d0 * y / delta)
!         !else
!         !    ep(id,np)%p%scalar(1) = 0.0d0
!         !end if
!
!         !if (y .ge. 0.0d0) then
!         !    eP(id,np)%p%scalar(1) = 0.0d0
!         !    ep(id,np)%p%scalar(2) = 1.0!erf(sqrt(pie)*y/delta)
!         !end if
!         !if (y .lt. 0.0d0) then
!         !    eP(id,np)%p%scalar(2) = 0.0d0
!         !    ep(id,np)%p%scalar(1) = 1.0!- erf(sqrt(pie)*y/delta)
!         !end if
!         !
!         !eP(id,np)%p%scalar(3) = eP(id,np)%p%Rhofp
!
!     end do
!
! end subroutine initScalars

!> This subroutine initializes the particles with their scalar quantities at
!! the initial condition. 
subroutine initScalarsShear(d,eP,id)
    use size
    use particle_definition
    use domain_definition
    implicit none
    type(domain)                                :: d
    type(ePointers), dimension(ngp,npart)       :: eP
    double precision                            :: x,y,z
    integer                                     :: np, id
    double precision, parameter                 :: pie		= 3.14159265358979323846264338327d0
    double precision                            :: delta
    double precision                            :: T,F,O,N,P,spike
    
    delta = 1.0d0
    spike = 2.088d0
    
    do np=1,d%np
        x = eP(id,np)%p%Xp(1)
        y = eP(id,np)%p%Xp(2)
        z = eP(id,np)%p%Xp(3)
        
!        eP(id,np)%p%scalar(1) = eP(id,np)%p%Tfp
!        eP(id,np)%p%scalar(2) = 0.01152d0
!        eP(id,np)%p%scalar(3) = 0.23042d0
!        eP(id,np)%p%scalar(4) = 0.0d0
!        ep(id,np)%p%scalar(5) = 0.75806d0
        
        
        T = eP(id,np)%p%Tfp
        
		! for shear layer from here
        !if(y>0.0d0) then
        !    !T = 1.0d0*(spike-(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
        !    F = 0.0d0
        !    O = 0.2331d0*erf(sqrt(pie)*y/delta)
        !    N = 0.7669d0*erf(sqrt(pie)*y/delta)
        !    P = 1.0d0-F-O-N
        !else
        !    !T = 1.0d0*(spike+(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
        !    F = -0.0115d0*erf(sqrt(pie)*y/delta)
        !    O = -0.2304d0*erf(sqrt(pie)*y/delta)
        !    N = -0.7581d0*erf(sqrt(pie)*y/delta)
        !    P = 1.0d0-F-O-N
        !end if
        
		!for CV reactor
		! F     = 0.01152d0
!         O     = 0.23042d0
!         N     = 0.75806d0
!         P   = 0.0d0

        !for ramp cavity with injector
        F     = 0.000001d0
        O     = 0.24194d0
        N     = 0.75806d0
        P     = 0.0d0

		!for Sod Problem
		!if (x > 0.50d0) then
		!	F 	= 0.0115d0
		!	O 	= 0.2304d0
		!	N 	= 0.7581d0
		!	P   = 0.0d0
		!else
		!	F = 0.0d0
		!	O = 0.23310d0
		!	N = 0.76690d0
		!	P = 0.0d0
		!end if

        eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        eP(id,np)%p%scalar(2) = F
        eP(id,np)%p%scalar(3) = O
        eP(id,np)%p%scalar(4) = P
        eP(id,np)%p%scalar(5) = N
        
        !eP(id,np)%p%scalar(1) = 3.0
        !eP(id,np)%p%scalar(2) = 4.0
        
        !if (y .ge. 0.0d0) then
        !    ep(id,np)%p%scalar(1) = erf(3.14159d0**0.5d0 * y / delta)
        !else
        !    ep(id,np)%p%scalar(1) = 0.0d0
        !end if
        
        !if (y .ge. 0.0d0) then
        !    eP(id,np)%p%scalar(1) = 0.0d0
        !    ep(id,np)%p%scalar(2) = 1.0!erf(sqrt(pie)*y/delta)
        !end if
        !if (y .lt. 0.0d0) then
        !    eP(id,np)%p%scalar(2) = 0.0d0
        !    ep(id,np)%p%scalar(1) = 1.0!- erf(sqrt(pie)*y/delta)
        !end if
        !
        !eP(id,np)%p%scalar(3) = eP(id,np)%p%Rhofp
        
    end do
    
end subroutine initScalarsShear


!> @brief Calculation of constant portion of subgrid mixing frequency
!!
!! Calculates the constant portion of \f$ \Omega_m \f$ in FMDF. Done in preprocessing to prevent
!! additional floating point calculations on each particle
subroutine calcOmegaConstants
    use physics
    implicit none
    
    constOmegaLeft  = (c_omega * (1.0d0/Sc))/(re * Sc)
    constOmegaRight = (c_omega) / (re * Sc_t)
    !constOmegaRight = 0.0d0
    write(*,*)'Values in calcOmegaConstants:'
    write(*,*)'constOmegaLeft',constOmegaLeft
    write(*,*)'constOmegaRight',constOmegaRight
    !write(*,*)'Sc',Sc
    !write(*,*)'Sc_t',Sc_t
    
end subroutine calcOmegaConstants

!> @brief Calculates the variable portion of subgrid mixing frequency
!!
!! \f$ \Omega_m \f$ using the constant portion defined earlier.
subroutine calcOmega(p,d)
    use particle_definition
    use domain_definition
    use physics
    implicit none
    type(particle)                  :: p    !< particles
    type(domain)                    :: d    !< domain
    double precision                :: oL, oR, delta, nu
    double precision                :: delta_g
    
    delta = delta_g(d,1,1,1)
    nu = p%nu_t
    !write(*,*)'delta:',delta
    !write(*,*)'nu:',nu
    
    oL = constOmegaLeft / (p%eRho * delta* delta)
    !write(*,*)'prho:',p%eRho
    !write(*,*)'oL:',oL
    oR = constOmegaRight * nu / (delta * delta)
    !write(*,*)'oR:',oR
    p%omega_m = oL + oR
    !write(*,*)'nu: ',nu
    !write(*,*)'omega', p%omega_m

end subroutine calcOmega