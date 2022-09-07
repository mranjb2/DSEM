!///////////////////////////////////////////////////////////////////////////////
!////////								////////
!////////	tstep.f90						////////
!////////								////////
!////////	contains:						////////
!////////								////////
!////////	   SUBROUTINE tstep(time,dt,d,ngrid,mrtr,nmort,resid)	////////
!////////								////////
!///////////////////////////////////////////////////////////////////////////////
!
!> @brief Take one time step over all of the grids.
      SUBROUTINE tstep(time,dt,d,ngrid,mrtr,nmort,resid) 
!
!     ..........................................................................
!     date: 11/17/98                                                     
!     routines called: fluxes                                           
!                      update                                           
!
!     applicability: All 3d programs
!
!     take one time step over all of the grids.                              
!     ..........................................................................
!                                                               
      USE domain_definition
      USE mortar_definition
      USE rk_coefs
      USE mpi_par
      USE User_Data
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain)                    :: d(ngp) 
      TYPE (mortar)                    :: mrtr(nmp) 
      DOUBLE PRECISION, DIMENSION(ngp) :: dt
!
	  IF (shock) THEN
          CALL Entropy(time,dt,d,ngrid,nmort,mrtr)
      END IF
!
      IF (smagorinsky) THEN
          CALL simpleSmag(d,ngridl)
          CALL nut_prolong(d,ngridl)
      END IF
!
!     -------------------------------------------------------------------
!     start of step.  zero "g" variables
!     -----------------------------------
!
      DO id = 1,ngrid
         DO nv = 1,neq
            DO l = 1,d(id)%nc(3)
               DO m = 1,d(id)%nc(2) 
                  DO n = 1,d(id)%nc(1)
                     d(id)%g_Q(n,m,l,nv) = 0.0d0
                  END DO
               END DO
            END DO
         END DO
      END DO
!
!
!     ---------------------------
!     runge-kutta to step in time
!     ---------------------------
!
      DO kk = 1,kord 
        krkiteration  = kk

         tk = time + dt(1)*brk(kk) 
!
!        ------------------
!        compute the fluxes                                   
!        ------------------
!
         CALL Fluxes(d,ngrid,mrtr,nmort,tk,dt)
!
!        -----------------------------------
!        compute time derivatives and update
!        -----------------------------------
!
         resid = 0.0d0
         DO id = 1,ngrid 
           IF (ngridedge(id) == 0) THEN

            CALL update(tk,dt(id),d(id),kk,resid,id)

           ENDIF
         END DO
!
! wait + add the viscous fluxes
!
      CALL wait_send_q_to_mortar()                                             
!
      ndummort = nmort-nmortd
      DO j = 1,nmort-nmortd
!
        IF (nmortedge(j,1) == 1) THEN                                          
          ndummort = ndummort+1
          IF (myid == nmortedge(j,3)) THEN
            jl   = nmortmpi(ndummort,1)
            itag = j
!
            CALL add_flux_mpi(mrtr(jl),d,itag)
          ENDIF
        ENDIF                                                                  
      ENDDO


      DO id = 1,ngrid 
        IF (ngridedge(id) == 1) THEN

          CALL update(tk,dt(id),d(id),kk,resid)

         ENDIF
       END DO
     ENDDO


      RETURN 
      END SUBROUTINE tstep
