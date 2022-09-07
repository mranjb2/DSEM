!> \file baselib.f90
!! Contains polynomial functions and Runge-Kutta coefficients
!///////////////////////////////////////////////////////////////////////////////////////////////
!////////										////////
!////////	baselib.f90								////////
!////////										////////
!////////	contains:								////////
!////////										////////
!////////	   DOUBLE PRECISION FUNCTION polyn(k,x,n,z)				////////
!////////	   DOUBLE PRECISION FUNCTION polynder(k,x,n,z)				////////
!////////	   SUBROUTINE rkcoefs(iord)						////////
!////////	   DOUBLE PRECISION FUNCTION elapsed_time(start_time,end_time,elapsed)	////////
!////////										////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!                                                                       
!> @brief Compute at \f$ x \f$ the lagrange polynomial \f$ k \f$ of degree 
!! \f$ n-1 \f$ whose zeros are given by the \f$ z(i) \f$.
      DOUBLE PRECISION FUNCTION polyn(k,x,n,z) 
!                                                                       
!     date: 2/13/95                                                     
!     routines called: none                                             
!     includes: none                                                    
!     applicability:all                                                 
!                                                                       
!                                                                       
!     compute at x the lagrange polynomial k of degree n-1              
!     whose zeros are given by the z(i)                                 
!                                                                       
      implicit double precision (a-h,o-z) 
      dimension z(n) 
!                                                                       
      if(k.eq.1)     then 
         polyn = (x - z(2))/(z(k) - z(2)) 
         do 100 j = 3,n 
            polyn = polyn*(x - z(j))/(z(k) - z(j)) 
  100    continue 
      else 
         polyn = (x - z(1))/(z(k) - z(1)) 
         do 200 j = 2,k-1 
            polyn = polyn*(x - z(j))/(z(k) - z(j)) 
  200    continue 
         do 300 j = k+1,n 
            polyn = polyn*(x - z(j))/(z(k) - z(j)) 
  300    continue 
      endif 
!                                                                       
      return 
      END                                           
!                                                                       
!///////////////////////////////////////////////////////////////////////
!                             
!> @brief Compute the derivative of the polynomial, \f$ h' \f$.                                          
      DOUBLE PRECISION FUNCTION polynder(k,x,n,z) 
!                                                                       
!     date: 2/13/95                                                     
!     routines called: none                                             
!     includes: none                                                    
!     applicability: all                                                
!                                                                       
!     compute the derivative of the polynomial, h'.                     
!                                                                       
      implicit double precision (a-h,o-z) 
      dimension z(n) 
      double precision hp,poly 
!                                                                       
      hp = 0.0d0 
      do 200 l = 1,n 
         if(l.eq.k)     go to 200 
         poly = 1.0d0 
         do 100 m = 1,n 
            if(m.eq.l)     go to 100 
            if(m.eq.k)     go to 100 
            poly = poly*(x - z(m))/(z(k) - z(m)) 
  100    continue 
         hp = hp + poly/(z(k) - z(l)) 
  200 continue 
      polynder = hp 
!                                                                       
      return 
      END                                           
!                                                                       
!///////////////////////////////////////////////////////////////////////
!                                                                       
!> @brief Set the Runge-Kutta integration coefficients for low storage Runge-Kutta scheme
!! @param iord is the order of the RK scheme, i.e. RK4 is iord=4
      SUBROUTINE rkcoefs(iord) 
!                                                                       
! SET THE RUNGE-KUTTA INTEGRATION COEFS FOR LOW STORAGE rUNGE-kUTTAS    
!                                                                       
!     date: 6/10/96                                                     
!     routines called: none                                             
!     includes: none                                                    
!     uses: rkcoefs                                                     
!     applicability: all                                                
!                                                                       
      USE rk_coefs
!                                                                       
      ark = 0.0d0 
      brk = 0.0d0 
      crk = 0.0d0 
!                                                                       
      select case(iord) 
!                                                                       
                                        ! forward Euler                 
         case(1) 
            kord = 1 
            ark(1) = 0.0d0 
            brk(1) = 0.0d0 
            crk(1) = 1.0d0
            dtconst = 1.d0
                                        ! second order modified euler method
         case(2) 
            kord = 2 
            ark(1) = 0.0d0 
            ark(2) = -1.0d0 
            brk(1) = 0.0d0 
            brk(2) = 1.0d0 
            crk(1) = 1.0d0 
            crk(2) = 0.5d0 
            dtconst = 1.d0
                                        ! Williamson 3rd order method   
         case(3) 
            kord = 3 
            ark(1) = 0.0d0 
            ark(2) = -5.d0/9.d0 
            ark(3) = -153.d0/128.d0 
            brk(1) = 0.0d0 
            brk(2) = 1.0d0/3.d0 
            brk(3) = 3.0d0/4.0d0
            crk(1) = 1.0d0/3.0d0 
            crk(2) = 15.0d0/16.0d0 
            crk(3) = 8.0d0/15.0d0
            dtconst = 1.7d0
                                        ! Carpenter/Kennedy 4th Order Method
         case(4) 
            kord = 5 
            ark(1) = 0.0d0 
            ark(2) = -0.4178904745d0 
            ark(3) = -1.192151694643d0 
            ark(4) = -1.697784692471d0 
            ark(5) = -1.514183444257d0 
            brk(1) = 0.0d0 
            brk(2) = 0.1496590219993d0 
            brk(3) = 0.3704009573644d0 
            brk(4) = 0.6222557631345d0 
            brk(5) = 0.9582821306748d0 
            crk(1) = 0.1496590219993d0 
            crk(2) = 0.3792103129999d0 
            crk(3) = 0.8229550293869d0 
            crk(4) = 0.6994504559488d0 
            crk(5) = 0.1530572479681d0
            dtconst = 2.8d0
                      ! method not implemented                          
         case default 
            write(6,*) 'maximum time step order is 4, choose another' 
            stop 
      end select 
!                                                                       
      return 
      END                                           
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Compute the elapsed time in minutes for a compuation.
!!
!! Will give the wrong answer for leap year near march 1 or any computation run
!! over new year's eve.
   DOUBLE PRECISION FUNCTION elapsed_time(start_time,end_time,elapsed)
!
!     compute the elapsed time in minutes for a compuation
!     will give the wrong answer for leap year near march 1 or any computation run
!     over new year's eve
!
      INTEGER,DIMENSION(8)  :: start_time,end_time,elapsed
      INTEGER,DIMENSION(12) :: no_of_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
!
!@
      elapsed = end_time - start_time
!
      elapsed_time = 0.0d0
      DO month = start_time(2),end_time(2)-1
         elapsed_time = no_of_days(month) + elapsed_time
      END DO
      elapsed_time = elapsed_time*24*60

      elapsed_time = elapsed_time + (elapsed(3)*24.d0 + elapsed(5))*60.d0 &
                                  + elapsed(6) &
                                  + (elapsed(7) + elapsed(8)/1000.d0)/60.d0
!
      RETURN
   END FUNCTION elapsed_time
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Returns .TRUE. if the argument is an NaN.
      LOGICAL FUNCTION Is_NaN(x)
!
!.......................................................................
!     Returns .TRUE. if the argument is an NaN
!.......................................................................
!
      DOUBLE PRECISION, INTENT(IN) :: x
      Is_NaN = .FALSE.
      IF (x /= x)     Is_NaN = .TRUE.
      RETURN
      END FUNCTION Is_NaN
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Returns .TRUE. if the argument is an INF.
      LOGICAL FUNCTION Is_INF(x)
!
!.......................................................................
!     Returns .TRUE. if the argument is an INF
!.......................................................................
!
      DOUBLE PRECISION, INTENT(IN) :: x
      Is_INF = .FALSE.
      IF (ABS(x) > HUGE(1.0d0))     Is_INF = .TRUE.
      RETURN
      END FUNCTION Is_INF
