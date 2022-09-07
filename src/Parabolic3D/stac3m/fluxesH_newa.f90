!> \file fluxesH_newa.f90
!! Computation of interface fluxes
!///////////////////////////////////////////////////////////////////////
!////////							////////
!////////	fluxesH.f90					////////
!////////							////////
!////////	contains:					////////
!////////							////////
!////////	      SUBROUTINE Fluxes(d,ngrid,mrtr,nmort,time)////////
!////////	      SUBROUTINE Int_fluxes(d)    		////////
!////////							////////
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the inviscid fluxes at all points
!!
!   
      SUBROUTINE Fluxes(d,ngrid,mrtr,nmort,time) 
!
!.......................................................................
!     compute the inviscid fluxes at all points
!
!     date: 11/9/98                                                     
!     routines called: Int_fluxes                                             
!                      Soln_to_Mrtr                                         
!                      interfaces                                         
!
!     applicability: hyperbolic equations
!.......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d(ngp) 
      TYPE (mortar) :: mrtr(nmp) 
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq):: Qbound 
!
!     --------------------------------------------------
!     initialize stuff   
!     --------------------------------------------------
!
      Qbound = 0.0d0

!
!     --------------------------------------------------
!     compute interior point fluxes and save edge values
!     on the mortars
!     --------------------------------------------------
!
      DO id = 1,ngrid 
       IF (ngridedge(id) == 1) THEN
         CALL Int_fluxes(d(id))
       ENDIF
      END DO
!
!     -------------------------
!     send fluxes to mortar with 
!     domains on different processors
!     -------------------------
!
      nreq = 0
      DO j = 1,nmort-nmortd
        IF (nmortedge(j,1) == 1) THEN
          IF (myid == nmortedge(j,3)) THEN
            nreq = nreq + 1
            jl   = nmortmpi(nmortdouble(j),1)
            itag = j
            CALL send_q_to_mortars_mpi(mrtr(jl),d,itag) 
            nreq = nreq + 1
            jl   = nmortmpi(nmortdouble(j),1)
            itag = nmortdouble(j)
            CALL send_q_to_mortars_mpi(mrtr(jl),d,itag) 
          ELSEIF (myid == nmortedge(j,2) ) THEN
            nreq = nreq + 1
            jl   = nmortmpi(j,1)
            itag = nmortdouble(j)
            CALL send_q_to_mortars_mpi(mrtr(jl),d,itag) 
            nreq = nreq + 1
            jl   = nmortmpi(j,1)
            itag = j
            CALL send_q_to_mortars_mpi(mrtr(jl),d,itag) 
          ENDIF
        ENDIF
      END DO
!
!     --------------------------------------------------
!     compute interior point fluxes and save edge values
!     on the mortars
!     --------------------------------------------------
!
      DO id = 1,ngrid 
       IF (ngridedge(id) == 0) THEN
         CALL Int_fluxes(d(id))
       ENDIF
      END DO
!
!     -------------------------
!     send fluxes to mortar with
!     domains on equal processors
!     -------------------------
!
      DO j = 1,nmort
        IF (nmortedge(j,1) == 0) THEN
          IF (myid == nmortedge(j,2)) THEN
            jl   = nmortmpi(j,1)
            itag = j
            CALL send_q_to_mortars(mrtr(jl),d,itag) 
          ENDIF
        ENDIF
      END DO

!
!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
      DO j = 1,nmort
          IF (myid == nmortedge(j,2) .and. nmortedge(j,1) /= 1) THEN
            jl   = nmortmpi(j,1)
            itag = j
            CALL mortar_fluxes(mrtr(jl),d,time,itag) 
            IF (mrtr(jl)%id(2) == 0) THEN ! this is a boundary face
               iface = mrtr(jl)%iface(1)
               idm   = nmortmpi(j,2)
               idml  = ngridmpi(idm,1)
 
               CALL DirichletBC(time,mrtr(jl),d(idml)%bcond(iface),Qbound)
      
            ELSE
!
!     ---------------------------------
!     Average the solution at the mortar
!     ---------------------------------
!
               CALL V_Avg_Soln(mrtr(jl),Qbound)
            END IF
!
!     ----------------------------------
!     Send mortar Q to the domain faces
!     ----------------------------------
 
           CALL send_q_to_faces(mrtr(jl),d,Qbound,itag)
 
          ENDIF
      END DO
!
!     ------------------------------------
!     Wait until Q is send to the mortar
!     for mortars with domains on different
!     processors
!     ------------------------------------
!
      CALL wait_send_q_to_mortar() 
!
!
!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
      DO j = 1,nmort
          IF (myid == nmortedge(j,2) .and. nmortedge(j,1) == 1) THEN
            jl   = nmortmpi(j,1)
            itag = j
            CALL mortar_fluxes(mrtr(jl),d,time,itag) 
            IF (mrtr(jl)%id(2) == 0) THEN ! this is a boundary face
               iface = mrtr(jl)%iface(1)
               idm   = nmortmpi(j,2)
               idml  = ngridmpi(idm,1)
 
               CALL DirichletBC(time,mrtr(jl),d(idml)%bcond(iface),Qbound)
 
            ELSE
!
!     ---------------------------------
!     Average the solution at the mortar
!     ---------------------------------
!
               CALL V_Avg_Soln(mrtr(jl),Qbound)
            END IF
!
!     ----------------------------------
!     Send mortar Q to the domain faces
!     ----------------------------------
 
           CALL send_q_to_faces(mrtr(jl),d,Qbound,itag)
 
          ENDIF
      END DO
!
!     ----------------------
!     compute viscous fluxes
!     ----------------------
!
      CALL Vis_Fluxes(d,ngrid,mrtr,nmort,time,dt)

!
      RETURN 
      END                                           
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the computational space fluxes at the interior lobatto points.
!!
!!     Compute the "X" direction fluxes on the lobatto-gauss-gauss points
!!     by first interpolating the solutions to the lobatto-gauss-gauss points
!!     and then computing the fluxes. Save the boundary and interface
!!     solutions along the faces to be later sent to the mortars. This code
!!     block is for between faces 6 and 4.
!!     Compute the "Y" direction fluxes on the lobatto/gauss points
!!     by first interpolating the solutions to the lobatto/gauss points
!!     and then computing the fluxes. Save the boundary and interface
!!     solutions along the faces to be later sent to the mortars. This block
!!     is for between faces 1 and 2.
!!     Compute the "Z" direction fluxes on the lobatto/gauss points
!!     by first interpolating the solutions to the lobatto/gauss points
!!     and then computing the fluxes. Save the boundary and interface
!!     solutions along the faces to be later sent to the mortars. This code
!!     block is for between faces 3 and 5.
!
      SUBROUTINE Int_fluxes(d) 
!
!.......................................................................
!     date: 11/9/98
!     routines called: fflux
!                      gflux
!                      hflux
!
!     applicability: hyperbolic equations, stac3m versions of the code            
!
!     compute the computational space fluxes at the interior lobatto points
!.......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE Material_Properties
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d
      DOUBLE PRECISION, DIMENSION(nx,ny,nz,neq)   :: Q
      DOUBLE PRECISION, DIMENSION(neq + mxprops)  :: Qv
      DOUBLE PRECISION, DIMENSION(neq)            :: f,g,h
!
!     Initialize
!
      Q = 0.0d0
!
!
!     ---------------------------------------------------------------------
!     Compute the "X" direction fluxes on the lobatto-gauss-gauss points
!     by first interpolating the solutions to the lobatto-gauss-gauss points
!     and then computing the fluxes. Save the boundary and interface
!     solutions along the faces to be later sent to the mortars. This code
!     block is for between faces 6 and 4.
!     ---------------------------------------------------------------------
!
      DO np = 1,num_prop
         Qv(neq+np) = material_property(np,d%material_id)
      END DO

      DO nv = 1,neq

         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2) 
               DO n = 1,d%ncg(1)
                  Q(n,m,l,nv)  = d%Q(n,m,l,nv)*d%jacob(n,m,l)
               END DO
            END DO
         END DO

      END DO
!
      CALL interpxa(nx,d%bx,Q,d%ncg(:),d%Qlgg,d%nc(:))
!
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2) 
            DO n = 2,d%nc(1)-1
               DO nv = 1,neq
                  Qv(nv) = d%Qlgg(n,m,l,nv)
               END DO
               CALL fflux(Qv,f)
               CALL gflux(Qv,g)
               CALL hflux(Qv,h)

               DO nv = 1,neq
                  d%f(n,m,l,nv) = d%gmet(1,1,n,m,l)*f(nv) + &
                                  d%gmet(1,2,n,m,l)*g(nv) + &
                                  d%gmet(1,3,n,m,l)*h(nv)
               END DO
            END DO
         END DO
      END DO
!
!     ---------------------------------------------------------------------
!     Compute the "Y" direction fluxes on the lobatto/gauss points
!     by first interpolating the solutions to the lobatto/gauss points
!     and then computing the fluxes. Save the boundary and interface
!     solutions along the faces to be later sent to the mortars. This block
!     is for between faces 1 and 2.
!     ---------------------------------------------------------------------
!
!
      CALL interpya(ny,d%by,Q,d%ncg(:),d%Qglg,d%nc(:))
!
      DO l = 1,d%ncg(3)
         DO m = 2,d%nc(2)-1 
            DO n = 1,d%ncg(1)
               DO nv = 1,neq
                  Qv(nv) = d%Qglg(n,m,l,nv)
               END DO
               CALL fflux(Qv,f)
               CALL gflux(Qv,g)
               CALL hflux(Qv,h)

               DO nv = 1,neq
                  d%g(n,m,l,nv) = d%gmet(2,1,n,m,l)*f(nv) + &
                                  d%gmet(2,2,n,m,l)*g(nv) + &
                                  d%gmet(2,3,n,m,l)*h(nv)
               END DO
            END DO
         END DO
      END DO
!
!     ---------------------------------------------------------------------
!     Compute the "Z" direction fluxes on the lobatto/gauss points
!     by first interpolating the solutions to the lobatto/gauss points
!     and then computing the fluxes. Save the boundary and interface
!     solutions along the faces to be later sent to the mortars. This code
!     block is for between faces 3 and 5.
!     ---------------------------------------------------------------------
!
      CALL interpza(nz,d%bz,Q,d%ncg(:),d%Qggl,d%nc(:))
!
      DO l = 2,d%nc(3)-1
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               DO nv = 1,neq
                  Qv(nv) = d%Qggl(n,m,l,nv)
               END DO
               CALL fflux(Qv,f)
               CALL gflux(Qv,g)
               CALL hflux(Qv,h)

               DO nv = 1,neq
                  d%h(n,m,l,nv) = d%gmet(3,1,n,m,l)*f(nv) + &
                                  d%gmet(3,2,n,m,l)*g(nv) + &
                                  d%gmet(3,3,n,m,l)*h(nv)
               END DO
            END DO
         END DO
      END DO
!
      RETURN
      END
