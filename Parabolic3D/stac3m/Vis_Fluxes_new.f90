!///////////////////////////////////////////////////////////////////////////////
!////////								////////
!////////	Vis_Fluxes.f90						////////
!////////								////////
!////////	contains:						//////// 
!////////	      SUBROUTINE Vis_Fluxes(d,ngrid,mrtr,nmort,time)	////////
!////////								////////
!///////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the viscous fluxes and add them to the inviscid fluxes 
!!
!! @code
!! USE domain_definition
!! @endcode
!
      SUBROUTINE Vis_Fluxes(d,ngrid,mrtr,nmort,time,dt) 
!
!.......................................................................
!     compute the viscous fluxes and add them to the inviscid fluxes    
!
!     DATE: 3/12/01                                                    
!
!     Modifications:                                                    
!
!     Routines called: DirichletBC
!                      V_Avg_Soln
!                      Gauss_Deriv
!                      Clear_Mortar
!
!     Applicability: parabolic equations
!.......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE mpi_par
      USE User_Data
      USE mpi
      USE input_data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d(ngp) 
      TYPE (mortar) :: mrtr(nmp) 
!
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq):: Qbound
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxlgg, Qylgg,Qzlgg    !< lobatto/gauss/gauss point derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxglg, Qyglg,Qzglg    !< gauss/lobatto/gauss point derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qxggl, Qyggl,Qzggl    !< gauss/gauss/lobatto point derivatives
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq) :: Qx, Qy, Qz            !< gauss/gauss/gauss   point derivatives
      !DOUBLE PRECISION, DIMENSION(ngp,nx,ny,nz,neq) :: r,s,t        !< viscous fluxes
      !DOUBLE PRECISION, DIMENSION(1,nx,ny,nz,neq) :: r,s,t        !< viscous fluxes
      DOUBLE PRECISION, DIMENSION(neq)          :: Q,Qx_loc,Qy_loc,Qz_loc !< pointwise quantities
      DOUBLE PRECISION, DIMENSION(neq)          :: f,g,h
      DOUBLE PRECISION, DIMENSION(ngp) :: dt

!
      DO 110 intproc = 1,2
!
!     ------------------------------------
!>     Compute the viscous fluxes by domain
!     ------------------------------------
!
      DO id = 1,ngrid

        IF ( ngridedge(id) == (-intproc+2)) THEN
!
!
!        ------------------------------------------------
!        Compute solution gradients at gauss/gauss points
!        ------------------------------------------------
!
!!!!! NOTE THAT GMETG IS USED INSTEAD OF GMET
         DO nv = 1,neq
            CALL Gauss_Deriv(d(id)%gmetg,d(id)%Qlgg(:,:,:,nv),d(id)%Qglg(:,:,:,nv), &
                             d(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv), &
			     d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)
         END DO
!
!        -----------------------------------------
!        Prolong the derivatives to the cell faces
!        -----------------------------------------
!

      DO nv = 1,neq
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
                CALL interp(nx,d(id)%bx,Qx(:,m,l,nv),d(id)%ncg(1),Qxlgg(:,m,l,nv),d(id)%nc(1))
                CALL interp(nx,d(id)%bx,Qy(:,m,l,nv),d(id)%ncg(1),Qylgg(:,m,l,nv),d(id)%nc(1))
                CALL interp(nx,d(id)%bx,Qz(:,m,l,nv),d(id)%ncg(1),Qzlgg(:,m,l,nv),d(id)%nc(1))
            END DO
         END DO

         DO l = 1,d(id)%ncg(3)
            DO n = 1,d(id)%ncg(1)
                CALL interp(ny,d(id)%by,Qx(n,:,l,nv),d(id)%ncg(2),Qxglg(n,:,l,nv),d(id)%nc(2))
               CALL interp(ny,d(id)%by,Qy(n,:,l,nv),d(id)%ncg(2),Qyglg(n,:,l,nv),d(id)%nc(2))
              CALL interp(ny,d(id)%by,Qz(n,:,l,nv),d(id)%ncg(2),Qzglg(n,:,l,nv),d(id)%nc(2))
            END DO
         END DO

         DO m = 1,d(id)%ncg(2)
            DO n = 1,d(id)%ncg(1)
               CALL interp(nz,d(id)%bz,Qx(n,m,:,nv),d(id)%ncg(3),Qxggl(n,m,:,nv),d(id)%nc(3))
               CALL interp(nz,d(id)%bz,Qy(n,m,:,nv),d(id)%ncg(3),Qyggl(n,m,:,nv),d(id)%nc(3))
               CALL interp(nz,d(id)%bz,Qz(n,m,:,nv),d(id)%ncg(3),Qzggl(n,m,:,nv),d(id)%nc(3))
            END DO
         END DO

      END DO
!
!        ------------------------------------------
!        Compute the viscous fluxes in X-direction
!        ------------------------------------------
!
         DO m = 1,d(id)%ncg(2)
           DO n = 1,d(id)%ncg(1)+1
             DO l = 1,d(id)%ncg(3)
               Q(1:neq)      = d(id)%Qlgg(n,m,l,1:neq)
               Qx_loc(1:neq) = Qxlgg(n,m,l,1:neq)
               Qy_loc(1:neq) = Qylgg(n,m,l,1:neq)
               Qz_loc(1:neq) = Qzlgg(n,m,l,1:neq)
               
               amu = 1.0d0
               IF (smagorinsky) amu = amu + d(id)%nu_t_lgg(n,m,l)*Q(1)
               IF (shock)       amu = amu + d(id)%muhlgg(n,m,l)

               CALL x_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,f, amu)
               CALL y_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,g, amu)
               CALL z_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,h, amu)

               r(id,n,m,l,1:neq)=d(id)%gmet(1,1,n,m,l)*f(1:neq) + d(id)%gmet(1,2,n,m,l)*g(1:neq) + &
                            d(id)%gmet(1,3,n,m,l)*h(1:neq) 
              END DO
            END DO
         END DO
!
!        ------------------------------------------
!        Compute the viscous fluxes in Y-direction
!        ------------------------------------------
!
         DO m = 1,d(id)%ncg(2) +1
           DO n = 1,d(id)%ncg(1)
             DO l = 1,d(id)%ncg(3)
               Q(1:neq)      = d(id)%Qglg(n,m,l,1:neq)
               Qx_loc(1:neq) = Qxglg(n,m,l,1:neq)
               Qy_loc(1:neq) = Qyglg(n,m,l,1:neq)
               Qz_loc(1:neq) = Qzglg(n,m,l,1:neq)
!
               amu = 1.0d0
               IF (smagorinsky) amu = amu + d(id)%nu_t_glg(n,m,l)*Q(1)
               IF (shock)       amu = amu + d(id)%muhglg(n,m,l)
               
               CALL x_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,f, amu)
               CALL y_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,g, amu)
               CALL z_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,h, amu)
!
               s(id,n,m,l,1:neq)=d(id)%gmet(2,1,n,m,l)*f(1:neq) + d(id)%gmet(2,2,n,m,l)*g(1:neq) + &
                            d(id)%gmet(2,3,n,m,l)*h(1:neq) 
               
              END DO
            END DO
         END DO

!
!        ------------------------------------------
!        Compute the viscous fluxes in Z-direction
!        ------------------------------------------
!
         DO m = 1,d(id)%ncg(2)
           DO n = 1,d(id)%ncg(1)
             DO l = 1,d(id)%ncg(3) +1
               Q(1:neq)      = d(id)%Qggl(n,m,l,1:neq)
               Qx_loc(1:neq) = Qxggl(n,m,l,1:neq)
               Qy_loc(1:neq) = Qyggl(n,m,l,1:neq)
               Qz_loc(1:neq) = Qzggl(n,m,l,1:neq)
!
               amu = 1.0d0
               IF (smagorinsky) amu = amu + d(id)%nu_t_ggl(n,m,l)*Q(1)
               IF (shock)       amu = amu + d(id)%muhggl(n,m,l)
               
               CALL x_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,f, amu)
               CALL y_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,g, amu)
               CALL z_vis_flux(Q,Qx_loc,Qy_loc,Qz_loc,h, amu)
!
               t(id,n,m,l,1:neq)=d(id)%gmet(3,1,n,m,l)*f(1:neq) + d(id)%gmet(3,2,n,m,l)*g(1:neq) + &
                           d(id)%gmet(3,3,n,m,l)*h(1:neq) 
              END DO
            END DO
         END DO
!
!        -----------------------------------------------------
!        Compute the full interior fluxes: Advective + Viscous
!        -----------------------------------------------------
!
        DO nv = 1,neq 

!         X-direction


          DO n = 2,d(id)%ncg(1)
            DO m = 1,d(id)%ncg(2)
              DO l = 1,d(id)%ncg(3)
                d(id)%f(n,m,l,nv) = d(id)%f(n,m,l,nv) - r(id,n,m,l,nv)/re
              END DO
            END DO
          END DO

!         Y-direction

          DO n = 1,d(id)%ncg(1)
            DO m = 2,d(id)%ncg(2)
              DO l = 1,d(id)%ncg(3)
                d(id)%g(n,m,l,nv) = d(id)%g(n,m,l,nv) - s(id,n,m,l,nv)/re
              END DO
            END DO
          END DO

!         Z-direction

          DO n = 1,d(id)%ncg(1)
            DO m = 1,d(id)%ncg(2)
              DO l = 2,d(id)%ncg(3)
                d(id)%h(n,m,l,nv) = d(id)%h(n,m,l,nv) - t(id,n,m,l,nv)/re
              END DO
            END DO
          END DO

        END DO

       ENDIF
      END DO ! done with all domains


      IF (intproc == 1) THEN
!
!     -------------------------
!     send fluxes to mortar with
!     domains on different processors
!     -------------------------
!
      nreq = 0
      ndummort = nmort-nmortd
      DO j = 1,nmort-nmortd
        IF (nmortedge(j,1) == 1) THEN
          ndummort = ndummort+1
          IF (myid == nmortedge(j,3)) THEN
            nreq = nreq + 1
            jl   = nmortmpi(ndummort,1)
            itag = j
            !CALL send_flux_to_mortars_mpi(mrtr(jl),r,s,t,itag,d)
            CALL send_flux_to_mortars_mpi(mrtr(jl),itag,d)
          ELSEIF (myid == nmortedge(j,2)) THEN
            nreq = nreq + 1
            jl   = nmortmpi(j,1)
            itag = j
            !CALL send_flux_to_mortars_mpi(mrtr(jl),r,s,t,itag,d)
            CALL send_flux_to_mortars_mpi(mrtr(jl),itag,d)
          ENDIF
        ENDIF
      ENDDO
!
      ENDIF
       
      IF (intproc == 2) THEN
!

!
!     -------------------------
!     send fluxes to mortar with
!     domains on equal processors
!     -------------------------
!
      DO j = 1,nmort-nmortd
        IF (nmortedge(j,1) == 0) THEN
          IF (myid == nmortedge(j,2)) THEN
            jl   = nmortmpi(j,1)
            itag = j
            !CALL send_flux_to_mortars(mrtr(jl),r,s,t,itag,d)
            CALL send_flux_to_mortars(mrtr(jl),itag,d)
          ENDIF
        ENDIF
      END DO                                                                  
!
      ENDIF

110   CONTINUE
      CALL wait_send_q_to_mortar()

!  
!
!       ---------------------------------------------------------------------- 
!       Implement Neumann boundary conditions and average the interface fluxes
!       ----------------------------------------------------------------------
!
      ndummort = nmort-nmortd
      nreq = 0
      DO j = 1,nmort-nmortd
!
        IF (nmortedge(j,1) == 1) THEN
          ndummort = ndummort+1
          IF (myid == nmortedge(j,2)) THEN
            nreq = nreq + 1
!
            jl   = nmortmpi(j,1)
            itag = j
!
!     ---------------------------------
!     Average the solution at the mortar
!     and send the flux to the copied
!     mortar
!     ---------------------------------
!
            CALL mortar_vis_flux(mrtr(jl),Qbound)
            CALL send_flux_to_faces_add_mpi(mrtr(jl),d,Qbound,itag)

          ELSEIF (myid == nmortedge(j,3)) THEN

            nreq = nreq + 1
            jl   = nmortmpi(ndummort,1)
            itag = j
!
            CALL send_flux_to_faces_add_mpi(mrtr(jl),d,Qbound,itag)
!
          ENDIF
        ENDIF
!
      END DO

      DO j = 1,nmort-nmortd
!
        IF (nmortedge(j,1) == 0) THEN
          IF (myid == nmortedge(j,2)) THEN
!
            jl   = nmortmpi(j,1)
            itag = j
!
            IF (mrtr(jl)%id(2) == 0) THEN ! this is a boundary face
               iface = mrtr(jl)%iface(1)
               idm   = nmortmpi(j,2)
               idml  = ngridmpi(idm,1)

               CALL Vis_Neumann(time,mrtr(jl),d(idml)%bcond(iface),Qbound)

            ELSE
!
!     ---------------------------------
!     Average the solution at the mortar
!     ---------------------------------
!
               CALL mortar_vis_flux(mrtr(jl),Qbound)
               !IF (myid==0 .and. j==60) THEN
               !  DO k=1,neq
               !  write(25,*) 'this is k'
               !  DO i=1,mrtr(jl)%len(2,1)
               !  DO j2=1,mrtr(jl)%len(2,1)
               !  write(25,*) Qbound(i,j2,k),mrtr(jl)%fv(i,j2,k,1),mrtr(jl)%fv(i,j2,k,2)
               !  ENDDO
               !  ENDDO
               !  ENDDO
               !ENDIF

            END IF
!
!     ---------------------------------------------------------------------
!     Send mortar fluxes to the domain faces and add to the inviscid fluxes
!     ---------------------------------------------------------------------
!
           
           CALL send_flux_to_faces_add(mrtr(jl),d,Qbound,itag)
           
     
      
!
          ENDIF
        ENDIF
!
      END DO



      RETURN
      END
