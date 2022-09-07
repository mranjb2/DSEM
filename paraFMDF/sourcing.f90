!> \file 
!! Chemical reaction source term routines
!
!> @brief Source term added to the FMDF scalar transport equations
!! This source term is modeled as:
!! \f[
!!      S_\alpha = \gamma Ce Da \rho Y_f Y_o \exp{- \frac{Ze}{T}}
!! \f]
subroutine singleStepSource(p, dt)
    use particle_definition
    use physics
    implicit none
    type(particle)                  :: p
    double precision, intent(in)    :: dt
    double precision                :: s
    double precision                :: Yf, Yo
    
    Yf = p%scalar(3)
    Yo = p%scalar(4)
    
    s = Dam * p%scalar(2) * Yf * Yo * exp( - Ze / p%scalar(1) )
!    p%react = s
    p%scalar(1) = (Ce*gamma*p%react*dt) + p%scalar(1)
    p%scalar(3) = -s*dt + p%scalar(3)
    p%scalar(4) = -s*dt + p%scalar(4)
    p%react = s
end subroutine singleStepSource

!> @brief Source term added to the FMDF scalar transport equations
!! This source term is modeled as:
!! \f[
!!      S_\alpha = \gamma Ce Da \rho Y_f Y_o \exp{- \frac{Ze}{T}}
!! \f]
subroutine singleStepSourceCp(p, dt)
    use particle_definition
    use physics
    implicit none
    type(particle)                  :: p
    double precision, intent(in)    :: dt
    double precision                :: s
    double precision                :: Y(4), rho,a,b,T
    double precision                :: omega_dot, rco(4), MW_dim(4), MW(4), k
    double precision                :: sum1,sum2,sum3,cp(4),h(4),gamma_var
    integer                         :: i
    
    !---- Dimensional quantities
    MW_dim(1) = 2.02d0
    MW_dim(2) = 32.0d0
    MW_dim(3) = 18.02d0
    MW_dim(4) = 28.01d0
    
    rco(1) = -1.0d0
    rco(2) = -.5d0
    rco(3) = 1.0d0
    rco(4) = 0.0d0
    
    !---- Nondimensionalization
    do i=1,4
        MW(i) = MW_dim(i) / MW_dim(1)
    end do
    
    a = -rco(1)
    b = -rco(2)
    
    Y(1) = p%scalar(2)
    Y(2) = p%scalar(3)
    Y(3) = p%scalar(4)
    Y(4) = p%scalar(5)
    !rho = p%scalar(2)
    !rho = p%rho_a
    T = p%scalar(1)
    
    !---- Calc variable properties
    call thermodynamics(T,Y,MW_dim,cp,h,gamma_var)
    !---- Reaction calc
    rho = p%P * gamma_var / T
    k = Dam*(rho**(a+b-1.0d0))*(Y(1)**a)*(Y(2)**b)*exp(-Ze/T) !rate of change
	k = max(k,0.0d0)
    
    sum1=0.0
    sum2=0.0
    sum3=0.0
    do i=1,4 !summation terms for source term
        sum1 = sum1 + T*rco(i)
        sum2 = sum2 + rco(i)*h(i)
        sum3 = sum3 + Y(i)*cp(i)/MW(i)
    end do
    
    p%react = gamma_var*k*((sum1-sum2)/sum3) ! reaction source term
    !if(p%react .ge. 10.0) write(*,*)'dt,',dt,'preact',p%react,'tb',p%scalar(1)
    p%scalar(1) = p%scalar(1) + p%react*dt       ! integrated reaction temperature
    !if(p%react .ge. 10.0) write(*,*)'tf',p%scalar(1)
    p%rho_a     = rho
    p%scalar(2) = p%scalar(2) + MW(1)*rco(1)*k*dt! integrated fuel
    p%scalar(3) = p%scalar(3) + MW(2)*rco(2)*k*dt
    p%scalar(4) = p%scalar(4) + MW(3)*rco(3)*k*dt
    if(p%scalar(2) < 0.0) p%scalar(2) = 0.0d0
    if(p%scalar(3) < 0.0) p%scalar(3) = 0.0d0
    if(p%scalar(4) < 0.0) p%scalar(4) = 0.0d0
    !write(*,*) p%scalar(5)
!    p%scalar(6) = p%scalar(6)! + MW(4)*rco(4)*k*dt    
    
!    s = Da * p%scalar(2) * Yf * Yo * exp( - Ze / p%scalar(1) )
!    p%react = s
!    p%scalar(1) = (Ce*gamma*p%react*dt) + p%scalar(1)
!    p%scalar(3) = -s*dt + p%scalar(3)
!    p%scalar(4) = -s*dt + p%scalar(4)
!    p%react = s
end subroutine singleStepSourceCp

!> @Brief Source term to be added to the carrier phase energy equation for a single step reaction
!! This source term is modelled as:
!! \f[
!!      S = \frac{Ce Da}{Ma_f^2 (\gamma -1)} \rho^2 Y_f Y_o \exp{- \frac{Ze}{T} }
!! \f]
subroutine fmdfTESource(d, id, n, m, l, dt, source)
    use domain_definition
    use physics
    implicit none
    type(domain)    , intent(in)    :: d
    double precision, intent(in)    :: dt
    double precision, intent(out)   :: source
    double precision                :: T, rhof
    integer         , intent(in)    :: n,m,l,id
    
    rhof = d%Q(n,m,l,1) * d%jacob(n,m,l)
    source = (rhof / (mach*mach*(gamma-1.0d0))) *  d%react(n,m,l)!(rhof / (mach*mach*gamma*(gamma-1.0d0))) *  d%react(n,m,l)
!    write(*,*)'rhof',rhof,'react',d%react(n,m,l),'source',source
end subroutine fmdfTESource

!///////////////////////////////////////////////////////////////////////
!

!-----------------------------------------------------------------------
!>   Thermodynamic Properties of Species
!!   This subroutine take non-dimentional temperature, mass fractions,
!!   dimensional molecular masses, and gives non-dimensional Cp's,
!!   non-dimensional enthalpy, and real mixture gamma
!-----------------------------------------------------------------------
!
  SUBROUTINE thermodynamics(T_nondim,Y,M,cp,h,gammas)
  USE size
  USE physics
!
  IMPLICIT NONE
  INTEGER :: i
  DOUBLE PRECISION  :: T,T_nondim,gammas,cp_mix,R_mix,M_mix,M_mix_1
  double precision  :: tref
  DOUBLE PRECISION, DIMENSION(4)      :: cp,h,M,Y
  DOUBLE PRECISION, DIMENSION(4,2,7)  :: a
  double precision, parameter::ru = 8.314d0
!
  tref = 800.0d0
  T = T_nondim * tref
!
! -------------- Curvfit Coefficients (1:lower, 2:higher)
!
  a(1,1,:) = (/ +3.29812400d+00, +8.24944200d-04, -8.14301500d-07, -9.47543400d-11, &
                +4.13487200d-13, -1.01252100d+03, -3.29409400d+00 /)  ! H2 lower
  a(1,2,:) = (/ +2.99142300d+00, +7.00064400d-04, -5.63382900d-08, -9.23157800d-12, &
                +1.58275200d-15, -8.35034000d+02, -1.35511000d+00 /)  ! H2 higher
  a(2,1,:) = (/ +3.21293600d+00, +1.12748600d-03, -5.75615000d-07, +1.31387700d-09, &
                -8.76855400d-13, -1.00524900d+03, +6.03473800d+00 /)  ! O2 lower
  a(2,2,:) = (/ +3.69757800d+00, +6.13519700d-04, -1.25884200d-07, +1.77528100d-11, &
                -1.13643500d-15, -1.23393000d+03, +3.18916600d+00 /)  ! O2 higher
  a(3,1,:) = (/ +3.38684200d+00, +3.47498200d-03, -6.35469600d-06, +6.96858100d-09, &
                -2.50658800d-12, -3.02081100d+04, +2.59023300d+00 /)  ! H2O lower
  a(3,2,:) = (/ +2.67214600d+00, +3.05629300d-03, -8.73026000d-07, +1.20099600d-10, &
                -6.39161800d-15, -2.98992100d+04, +6.86281700d+00 /)  ! H2O higher
  a(4,1,:) = (/ +3.29867700d+00, +1.40824000d-03, -3.96322200d-06, +5.64151500d-09, &
                -2.44485500d-12, -1.02090000d+03, +3.95037200d+00 /)  ! N2 lower
  a(4,2,:) = (/ +2.92664000d+00, +1.48797700d-03, -5.68476100d-07, +1.00970400d-10, &
                -6.75335100d-15, -9.22797700d+02, +5.98052800d+00 /)  ! N2 higher
!
! ------------- Specific Heat Capacity (Cp) and Sensible Enthalpy (hs)
!               Turn's Combustion Book (Appendix table A.13)
!
  DO i = 1,4
    IF (T < 1000d0) THEN
      cp(i) = a(i,1,1) + a(i,1,2)*T + a(i,1,3)*T**2 + a(i,1,4)*T**3 + a(i,1,5)*T**4
      h(i)  = a(i,1,1) + a(i,1,2)*T/2.0d0 + a(i,1,3)*T**2/3.0d0 &
            + a(i,1,4)*T**3/4.0d0 + a(i,1,5)*T**4/5.0d0 + a(i,1,6)/T
    ELSE
      cp(i) = a(i,2,1) + a(i,2,2)*T + a(i,2,3)*T**2 + a(i,2,4)*T**3 + a(i,2,5)*T**4
      h(i)  = a(i,2,1) + a(i,2,2)*T/2.0d0 + a(i,2,3)*T**2/3.0d0 &
            + a(i,2,4)*T**3/4.0d0 + a(i,2,5)*T**4/5.0d0 + a(i,2,6)/T
    END IF
    h(i) = h(i) * T_nondim
  END DO
!
! ------------- gamma
!
  cp_mix = 0.0d0
  DO i = 1,4
    cp_mix = cp_mix + Y(i)*cp(i)*ru/M(i) ! mixture heat capacity
  END DO
  M_mix_1 = 0.0d0
  DO i = 1,4
    M_mix_1 = M_mix_1 + Y(i)/M(i)
  END DO
  M_mix = 1.0d0/M_mix_1               ! mixture molecular weight
  R_mix = ru/M_mix                    ! mixture specific gas constant
  gammas = cp_mix/(cp_mix-R_mix)      ! mixture gamma
!
  END SUBROUTINE thermodynamics
!
!///////////////////////////////////////////////////////////////////////
