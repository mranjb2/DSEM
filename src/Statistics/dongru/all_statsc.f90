!///////////////////////////////////////////////////////////////////////////////
!////////                                                               ////////
!////////       all_stats.f90                                           ////////
!////////                                                               ////////
!////////       contains:                                               ////////
!////////                                                               ////////
!////////             SUBROUTINE Stats_All                              ////////
!////////             SUBROUTINE Stats_Average                          ////////
!////////             SUBROUTINE Stats_Fluct                            ////////
!////////                                                               ////////
!///////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Stats_All(d,ngrid,stats,statsfluct,mrtr,nmort)
!
!.......................................................................
!     compute the average fluid velocity
!
!     DATE: 06/19/01
!
!.......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE stats_definition
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (domain)     :: d(ngp)
      TYPE (mortar)     :: mrtr(nmp)
      TYPE (statfl)     :: stats(ngp)
      TYPE (statfluct)  :: statsfluct(ngp)
!
      INTEGER           :: ngrid,nmort
!
!   Compute fluid averaged statistics  
!
      IF (average) THEN
        nsample = nsample + 1
        DO id = 1,ngrid
          CALL Stats_Average(d(id),stats(id))      
        END DO
      ENDIF
!   Compute fluctuating average
      IF (rms) THEN
         nsample = nsample + 1
         CALL Stats_Fluct_Reynolds(d,ngrid,stats,statsfluct,mrtr,nmort)       
         CALL Stats_Fluct_Favre(d,ngrid,stats,statsfluct,mrtr,nmort)      
      ENDIF
!
      RETURN
      END SUBROUTINE Stats_All
!
!///////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Stats_Average(d,stats)
!
!.......................................................................
!     compute the Reynolds and Favre average fluid properties
!     assuming periodicity in y-direction and assuming
!     a square domain
!
!     Q_av(1)  = rho       ! Reynolds averaged density
!     Q_av(2-4)= rho*u_i   ! Favre averaged velocities
!     Q_av(5)  = p         ! Reynolds averaged pressure
!     Q_av(6-8)= u_i       ! Reynolds averaged velocities
!     Q_av(9)  = T         ! Reynolds averaged temperature
!
!     DATE: 06/19/01
!
!.......................................................................
!
      USE stats_definition
      USE domain_definition
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (statfl)   :: stats
      TYPE (domain) :: d
!
      DOUBLE PRECISION :: velocity(ny,3),pressure(ny),temperature(ny)
!
     pressure = 0.0d0
     velocity = 0.0d0
     temperature = 0.0d0
!
     ncg2 = d%ncg(2)
!
     DO n = 1,d%ncg(1)
       DO m = 1,d%ncg(3)
         DO j= 1,d%ncg(2)
!
! Compute the instanteneous velocity, pressure and temperature for purpose
! of determining the Reynolds ensemble average of these quantities
!
           DO i=1,3
              velocity(j,i) = d%Q(n,j,m,i+1)/d%Q(n,j,m,1)
           ENDDO
           pressure(j) = (d%Q(n,j,m,5) - (d%Q(n,j,m,2)**2+d%Q(n,j,m,3)**2+d%Q(n,j,m,4)**2) &
                   /(d%Q(n,j,m,1)*2.0d0))*(gamma-1.0d0) ! the jacobian is NOT scaled out
           
           temperature(j) = pressure(j)*gamma/d%Q(n,j,m,1)
!          ENDDO
! 
! Determine the ensemble averages
!
         DO nv = 1,neq-1
            stats%Q_av(n,j,m,nv) = stats%Q_av(n,j,m,nv) + d%Q(n,j,m,nv)*d%jacob(n,j,m)
         ENDDO
         stats%Q_av(n,j,m,5)     = stats%Q_av(n,j,m,5)  + pressure(j)*d%jacob(n,j,m)
         DO i=1,3
           stats%Q_av(n,j,m,5+i) = stats%Q_av(n,j,m,5+i)+ velocity(j,i)
         ENDDO
         stats%Q_av(n,j,m,9)     = stats%Q_av(n,j,m,9)  + temperature(j)
         
        END DO
       END DO
     END DO
!      
      !write(*,*) stats%Q_av(1,1,1,1)
      RETURN
      END SUBROUTINE Stats_Average

!
!///////////////////////////////////////////////////////////////////////////////
!
!
      SUBROUTINE Stats_Fluct_Reynolds(d,ngrid,stats,sfl,mrtr,nmort)
!
!.......................................................................
!
!     Compute the Reynolds fluctuating averages (turbulence statistics)
!
!     Q_fluct(3) = d( u_i'*p')/dx_i      ! transport through pressure fluct.
!     Q_fluct(4) = p'd(u_i')/dx_i        ! pressure-dillation 
!     Q_fluct(5) = tau_ik'*d(u_i')/dx_k  ! rest term turbulence dissipation
!                  -Q_fluct(6) - Q_fluct(7)
!     Q_fluct(6) = omega_i'^2            ! solenoidal dissipation
!     Q_fluct(7) = divergence'^2         ! dillational dissipation
!
!     DATE: 05/05/03
!
!.......................................................................
!
      USE stats_definition
      USE mortar_definition
      USE domain_definition
      USE physics
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (domain)      :: d(ngp)
      TYPE (mortar)      :: mrtr(nmp)
      TYPE (statfl)      :: stats(ngp)
      TYPE (statfluct)   :: sfl(ngp)
!
      ncg1 = d(1)%ncg(1)
      ncg2 = d(1)%ncg(2)
      ncg3 = d(1)%ncg(3)
     
!
!  Determine Reynolds fluctuations
!
     r=0.0d0
     s=0.0d0
     t=0.0d0
     DO id=1,ngrid
      DO n=1,ncg1
        DO j=1,ncg2
          DO m=1,ncg3
            pressure = (d(id)%Q(n,j,m,5) - (d(id)%Q(n,j,m,2)**2+d(id)%Q(n,j,m,3)**2+d(id)%Q(n,j,m,4)**2) &
                   /(d(id)%Q(n,j,m,1)*2.0d0))*(gamma-1.0d0)*d(id)%jacob(n,j,m)
            u        = d(id)%Q(n,j,m,2)/d(id)%Q(n,j,m,1)
            v        = d(id)%Q(n,j,m,3)/d(id)%Q(n,j,m,1)
            w        = d(id)%Q(n,j,m,4)/d(id)%Q(n,j,m,1)
            r(id,n,j,m,1) = stats(id)%Q_av(n,j,m,5) - pressure      ! pressure fluctuation
            r(id,n,j,m,2) = stats(id)%Q_av(n,j,m,6) - u             ! x-vel fluctuation Reynolds
            r(id,n,j,m,3) = stats(id)%Q_av(n,j,m,7) - v             ! y-vel fluctuation Reynolds
            r(id,n,j,m,4) = stats(id)%Q_av(n,j,m,8) - w             ! z-vel fluctuation Reynolds
          ENDDO
        ENDDO
      ENDDO
     ENDDO
!
! Determine gradients of Reynolds fluctuations
!
     CALL compute_derivative(d,ngrid,mrtr,nmort)

!
! Determine terms with Reynolds fluctuations
!
     DO id=1,ngrid
       DO n=1,ncg1
         DO m=1,ncg3

           DO j=1,ncg2
            pressure = (d(id)%Q(n,j,m,5) - (d(id)%Q(n,j,m,2)**2+d(id)%Q(n,j,m,3)**2+d(id)%Q(n,j,m,4)**2) &
                   /(d(id)%Q(n,j,m,1)*2.0d0))*(gamma-1.0d0)*d(id)%jacob(n,j,m)
            u        = d(id)%Q(n,j,m,2)/d(id)%Q(n,j,m,1)
            v        = d(id)%Q(n,j,m,3)/d(id)%Q(n,j,m,1)
            w        = d(id)%Q(n,j,m,4)/d(id)%Q(n,j,m,1)
            Temp     = pressure*gamma/d(id)%Q(n,j,m,1)
            pprime   = stats(id)%Q_av(n,j,m,5) - pressure    
            uprime   = stats(id)%Q_av(n,j,m,6) - u     
            vprime   = stats(id)%Q_av(n,j,m,7) - v    
            wprime   = stats(id)%Q_av(n,j,m,8) - w   
            Tprime   = stats(id)%Q_av(n,j,m,9) - Temp  

            om1prime = t(id,n,j,m,3) - s(id,n,j,m,4)
            om2prime = t(id,n,j,m,2) - r(id,n,j,m,4)
            om3prime = - s(id,n,j,m,2) + r(id,n,j,m,3)
            divprime = r(id,n,j,m,2) + s(id,n,j,m,3) + t(id,n,j,m,4)
            restprime= -r(id,n,j,m,2)*(s(id,n,j,m,3) + t(id,n,j,m,4)) - &
                       s(id,n,j,m,3)*t(id,n,j,m,4) +  &      
                       s(id,n,j,m,2)*r(id,n,j,m,3) +  &      
                       t(id,n,j,m,2)*r(id,n,j,m,4) +  &      
                       t(id,n,j,m,3)*s(id,n,j,m,4)
! eps_1
            eps_1      =  r(id,n,j,m,2)* (r(id,n,j,m,2)+r(id,n,j,m,2)) + &
                          r(id,n,j,m,3)* (r(id,n,j,m,3)+s(id,n,j,m,2)) + &
                          r(id,n,j,m,4)* (r(id,n,j,m,4)+t(id,n,j,m,2)) - &
                          (2.0/3.0)*r(id,n,j,m,2)*divprime
            eps_1      =  eps_1 +                                        &
                          s(id,n,j,m,2)* (s(id,n,j,m,2)+r(id,n,j,m,3)) + &
                          s(id,n,j,m,3)* (s(id,n,j,m,3)+s(id,n,j,m,3)) + &
                          s(id,n,j,m,4)* (s(id,n,j,m,4)+t(id,n,j,m,3)) - &
                          (2.0/3.0)*s(id,n,j,m,3)*divprime
            eps_1      =  eps_1 +                                        &
                          t(id,n,j,m,2)* (t(id,n,j,m,2)+r(id,n,j,m,4)) + &
                          t(id,n,j,m,3)* (t(id,n,j,m,3)+s(id,n,j,m,4)) + &
                          t(id,n,j,m,4)* (t(id,n,j,m,4)+t(id,n,j,m,4)) - &
                          (2.0/3.0)*t(id,n,j,m,4)*divprime
! tau_prime
            tau11prime = 2.0d0*r(id,n,j,m,2) - 2.0d0*(divprime)/3.0d0
            tau22prime = 2.0d0*s(id,n,j,m,3) - 2.0d0*(divprime)/3.0d0
            tau33prime = 2.0d0*t(id,n,j,m,4) - 2.0d0*(divprime)/3.0d0
            tau13prime = r(id,n,j,m,4) + t(id,n,j,m,2)
            tau12prime = r(id,n,j,m,3) + s(id,n,j,m,2)
            tau23prime = s(id,n,j,m,4) + t(id,n,j,m,3)
!
!             sum1     = sum1 + pprime*divprime
!             sum2     = sum2 + uprime*r(id,n,j,m,1) +vprime*s(id,n,j,m,1) +wprime*t(id,n,j,m,1) 
!             sum3     = sum3 + 4.0d0*restprime
!             sum4     = sum4 + (om1prime**2+om2prime**2+om3prime**2)
!             sum5     = sum5 + (4.0/3.0)*(divprime**2)
! ! Reynolds averaged Reynolds stresses
!             sum6     = sum6 + uprime*uprime
!             sum7     = sum7 + vprime*vprime
!             sum8     = sum8 + wprime*wprime
!             sum9     = sum9 + uprime*vprime
!             sum10    = sum10 + uprime*wprime
!             sum11    = sum11 + vprime*wprime
! ! Reynolds average velocity temp average
!             sum12    = sum12 + uprime*Tprime
!             sum13    = sum13 + vprime*Tprime
!             sum14    = sum14 + wprime*Tprime
! ! Reynolds average velocity pressure average
!             sum15    = sum15 + uprime*pprime
!             sum16    = sum16 + vprime*pprime
!             sum17    = sum17 + wprime*pprime
! ! stress velocity average
!             sum18    = sum18 + uprime*tau11prime + vprime*tau12prime + wprime*tau13prime
!             sum19    = sum19 + uprime*tau12prime + vprime*tau22prime + wprime*tau23prime
!             sum20    = sum20 + uprime*tau13prime + vprime*tau23prime + wprime*tau33prime
! ! turbulence kinetic energy
!             sum21    = sum21 + (uprime*uprime + vprime*vprime + wprime*wprime)*0.5
! ! eps_1
!             sum22    = sum22 + eps_1

!            ENDDO
!
           sfl(id)%Q_fluct(n,j,m,3) =  sfl(id)%Q_fluct(n,j,m,3) + pprime*divprime
           sfl(id)%Q_fluct(n,j,m,4) =  sfl(id)%Q_fluct(n,j,m,4) + &
           uprime*r(id,n,j,m,1) +vprime*s(id,n,j,m,1) +wprime*t(id,n,j,m,1) 
           sfl(id)%Q_fluct(n,j,m,5) =  sfl(id)%Q_fluct(n,j,m,5) + (4.0d0*restprime)/re
           sfl(id)%Q_fluct(n,j,m,6) =  sfl(id)%Q_fluct(n,j,m,6) + &
           (om1prime**2+om2prime**2+om3prime**2)/re
           sfl(id)%Q_fluct(n,j,m,7) =  sfl(id)%Q_fluct(n,j,m,7) + ((4.0/3.0)*(divprime**2))/re
!
           sfl(id)%Q_fluct(n,j,m,19) =  sfl(id)%Q_fluct(n,j,m,19) + uprime*uprime
           sfl(id)%Q_fluct(n,j,m,20) =  sfl(id)%Q_fluct(n,j,m,20) + vprime*vprime
           sfl(id)%Q_fluct(n,j,m,21) =  sfl(id)%Q_fluct(n,j,m,21) + wprime*wprime
           sfl(id)%Q_fluct(n,j,m,22) =  sfl(id)%Q_fluct(n,j,m,22) + uprime*vprime
           sfl(id)%Q_fluct(n,j,m,23) =  sfl(id)%Q_fluct(n,j,m,23) + uprime*wprime
           sfl(id)%Q_fluct(n,j,m,24) =  sfl(id)%Q_fluct(n,j,m,24) + vprime*wprime
!
           sfl(id)%Q_fluct(n,j,m,25) =  sfl(id)%Q_fluct(n,j,m,25) + uprime*Tprime
           sfl(id)%Q_fluct(n,j,m,26) =  sfl(id)%Q_fluct(n,j,m,26) + vprime*Tprime
           sfl(id)%Q_fluct(n,j,m,27) =  sfl(id)%Q_fluct(n,j,m,27) + wprime*Tprime
!
           sfl(id)%Q_fluct(n,j,m,28) =  sfl(id)%Q_fluct(n,j,m,28) + uprime*pprime
           sfl(id)%Q_fluct(n,j,m,29) =  sfl(id)%Q_fluct(n,j,m,29) + vprime*pprime
           sfl(id)%Q_fluct(n,j,m,30) =  sfl(id)%Q_fluct(n,j,m,30) + wprime*pprime
!
           sfl(id)%Q_fluct(n,j,m,31) =  sfl(id)%Q_fluct(n,j,m,31) + &
           (uprime*tau11prime + vprime*tau12prime + wprime*tau13prime)/re
           
           sfl(id)%Q_fluct(n,j,m,32) =  sfl(id)%Q_fluct(n,j,m,32) + &
           (uprime*tau12prime + vprime*tau22prime + wprime*tau23prime)/re
           
           sfl(id)%Q_fluct(n,j,m,33) =  sfl(id)%Q_fluct(n,j,m,33) + &
           (uprime*tau13prime + vprime*tau23prime + wprime*tau33prime)/re
!
           sfl(id)%Q_fluct(n,j,m,34) =  sfl(id)%Q_fluct(n,j,m,34) + &
           (uprime*uprime + vprime*vprime + wprime*wprime)*0.5
!
           sfl(id)%Q_fluct(n,j,m,42) =  sfl(id)%Q_fluct(n,j,m,42) + eps_1/re
           ENDDO
         ENDDO
       ENDDO
     ENDDO
!
     r=0.0d0
     s=0.0d0
     t=0.0d0

!
        RETURN
      END SUBROUTINE Stats_Fluct_Reynolds
!  
!///////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Stats_Fluct_Favre(d,ngrid,stats,sfl,mrtr,nmort)
!
!.......................................................................
!
!     Compute the Favre fluctuating averages (turbulence statistics)
!
!     Q_fluct(1)    = rho*u''*T''            ! turbulent heat flux in x-direction
!     Q_fluct(2)    = rho*v''*T''            ! turbulent heat flux in y-direction
!     Q_fluct(8-13) = rho*u_i''*uj''         ! Reynolds stress tensor
!     Q_fluct(14)   = rho*w''*T''            ! turbulent heat flux in z-direction
!
!     DATE: 05/05/03
!
!.......................................................................
!
      USE stats_definition
      USE mortar_definition
      USE domain_definition
      USE physics
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (domain)      :: d(ngp)
      TYPE (mortar)      :: mrtr(nmp)
      TYPE (statfl)      :: stats(ngp)
      TYPE (statfluct)   :: sfl(ngp)
!
      ncg1 = d(1)%ncg(1)
      ncg2 = d(1)%ncg(2)
      ncg3 = d(1)%ncg(3)
     
!
! Determine terms with Favre fluctuations
!
     DO id=1,ngrid
       DO n=1,ncg1
         DO m=1,ncg3
 
!
!
           DO j=1,ncg2
            rho      = d(id)%Q(n,j,m,1)*d(id)%jacob(n,j,m)
            u        = d(id)%Q(n,j,m,2)/d(id)%Q(n,j,m,1)
            v        = d(id)%Q(n,j,m,3)/d(id)%Q(n,j,m,1)
            w        = d(id)%Q(n,j,m,4)/d(id)%Q(n,j,m,1)
            uav      = stats(id)%Q_av(n,j,m,2)/stats(id)%Q_av(n,j,m,1)
            vav      = stats(id)%Q_av(n,j,m,3)/stats(id)%Q_av(n,j,m,1)
            wav      = stats(id)%Q_av(n,j,m,4)/stats(id)%Q_av(n,j,m,1)
            uprime   = u - stats(id)%Q_av(n,j,m,2)/stats(id)%Q_av(n,j,m,1)  
            vprime   = v - stats(id)%Q_av(n,j,m,3)/stats(id)%Q_av(n,j,m,1)  
            wprime   = w - stats(id)%Q_av(n,j,m,4)/stats(id)%Q_av(n,j,m,1)  
            temp_av  = stats(id)%Q_av(n,j,m,5)*1.4d0/stats(id)%Q_av(n,j,m,1)
!            press    = d(id)%Q(n,j,m,5) - 0.5d0*rho*(u**2+v**2+w**2)
            press    = d(id)%Q(n,j,m,5) - 0.5d0*d(id)%Q(n,j,m,1)*(u**2+v**2+w**2)
            press    = press*d(id)%jacob(n,j,m)*0.4d0
            Tprime   = press*1.4d0/rho - temp_av 
            bKprime   =  uav*uprime + vav*vprime + wav*wprime
            turbk    = (uprime*uprime + vprime*vprime + wprime*wprime)*0.5d0
! moving average determination
            smallkprime = turbk - (sfl(id)%Q_fluct(n,j,m,35)/stats(id)%Q_av(n,j,m,1))/(nsample) 
!
! Favre Reynolds stresses
!             sum1     = sum1 + rho*uprime*uprime 
!             sum2     = sum2 +  rho*vprime*vprime
!             sum3     = sum3 +  rho*wprime*wprime
!             sum4     = sum4 +  rho*uprime*vprime
!             sum5     = sum5 +  rho*uprime*wprime
!             sum6     = sum6 +  rho*vprime*wprime
! ! Favre velocity temp average (turbulent heat flux)
!             sum7     = sum7 +  rho*uprime*Tprime
!             sum8     = sum8 +  rho*vprime*Tprime
!             sum9     = sum9 +  rho*wprime*Tprime
! ! Favre fluctuating average
!             sum10    = sum10 + uprime
!             sum11    = sum11 + vprime
!             sum12    = sum12 + wprime
!             sum13    = sum13 + Tprime
! ! turbulence kinetic energy
!             sum14    = sum14 + rho*turbk
! ! u''K''
!             sum15    = sum15 + rho*uprime*bKprime
!             sum16    = sum16 + rho*vprime*bKprime
!             sum17    = sum17 + rho*wprime*bKprime
! ! u''k''
!             sum18    = sum18 + rho*uprime*smallkprime
!             sum19    = sum19 + rho*vprime*smallkprime
!             sum20    = sum20 + rho*wprime*smallkprime
!

!            ENDDO
!
           sfl(id)%Q_fluct(n,j,m,8)  = sfl(id)%Q_fluct(n,j,m,8)  + rho*uprime*uprime 
           sfl(id)%Q_fluct(n,j,m,9)  = sfl(id)%Q_fluct(n,j,m,9)  + rho*vprime*vprime
           sfl(id)%Q_fluct(n,j,m,10) = sfl(id)%Q_fluct(n,j,m,10) + rho*wprime*wprime
           sfl(id)%Q_fluct(n,j,m,11) = sfl(id)%Q_fluct(n,j,m,11) + rho*uprime*vprime
           sfl(id)%Q_fluct(n,j,m,12) = sfl(id)%Q_fluct(n,j,m,12) + rho*uprime*wprime
           sfl(id)%Q_fluct(n,j,m,13) = sfl(id)%Q_fluct(n,j,m,13) + rho*vprime*wprime
           sfl(id)%Q_fluct(n,j,m,1)  = sfl(id)%Q_fluct(n,j,m,1)  + rho*uprime*Tprime
           sfl(id)%Q_fluct(n,j,m,2)  = sfl(id)%Q_fluct(n,j,m,2)  + rho*vprime*Tprime
           sfl(id)%Q_fluct(n,j,m,14) = sfl(id)%Q_fluct(n,j,m,14) + rho*wprime*Tprime
           sfl(id)%Q_fluct(n,j,m,15) = sfl(id)%Q_fluct(n,j,m,15) + uprime
           sfl(id)%Q_fluct(n,j,m,16) = sfl(id)%Q_fluct(n,j,m,16) + vprime
           sfl(id)%Q_fluct(n,j,m,17) = sfl(id)%Q_fluct(n,j,m,17) + wprime
           sfl(id)%Q_fluct(n,j,m,18) = sfl(id)%Q_fluct(n,j,m,18) + Tprime
!
           sfl(id)%Q_fluct(n,j,m,35) = sfl(id)%Q_fluct(n,j,m,35) + rho*turbk
!
           sfl(id)%Q_fluct(n,j,m,36) = sfl(id)%Q_fluct(n,j,m,36) + rho*uprime*bKprime
           sfl(id)%Q_fluct(n,j,m,37) = sfl(id)%Q_fluct(n,j,m,37) + rho*vprime*bKprime
           sfl(id)%Q_fluct(n,j,m,38) = sfl(id)%Q_fluct(n,j,m,38) + rho*wprime*bKprime
!           
           sfl(id)%Q_fluct(n,j,m,39) = sfl(id)%Q_fluct(n,j,m,39) + rho*uprime*smallkprime
           sfl(id)%Q_fluct(n,j,m,40) = sfl(id)%Q_fluct(n,j,m,40) + rho*vprime*smallkprime
           sfl(id)%Q_fluct(n,j,m,41) = sfl(id)%Q_fluct(n,j,m,41) + rho*wprime*smallkprime
           ENDDO
         ENDDO
       ENDDO
     ENDDO
!
      RETURN
      
    END SUBROUTINE Stats_Fluct_Favre
