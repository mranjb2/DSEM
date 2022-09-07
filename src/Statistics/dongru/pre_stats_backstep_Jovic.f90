!
! 
      SUBROUTINE compute_average_stats(stats,ncg,ngrid)
!

!.......................................................................
!     compute the average along a plane
!
!     DATE: 06/05/03
!
!.......................................................................
!
      USE stats_definition
      USE input_data
      USE User_Data
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (statfl)          :: stats(ngp)
!
      INTEGER                :: ngrid

        DO id = 1,ngrid
           stats(id)%Q_av(:,:,:,:) = stats(id)%Q_av(:,:,:,:)/(nsample)
        ENDDO

        RETURN
      END SUBROUTINE compute_average_stats
! 
!////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE compute_average_statsfluct(statsfluct,stats,ncg,ngrid)
!

!.......................................................................
!     compute the average along a plane
!
!     DATE: 06/05/03
!
!.......................................................................
!
      USE stats_definition
      USE input_data
      USE User_Data
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (statfluct)          :: statsfluct(ngp)
      TYPE (statfl)             :: stats(ngp)
!
      INTEGER                :: ngrid
!
!
         DO id = 1,ngrid
             statsfluct(id)%Q_fluct(:,:,:,:) = statsfluct(id)%Q_fluct(:,:,:,:)/(nsample)
         ENDDO
!
!
! Once the ensemble averaging is complete
! divide the necessary statsfluct by rho_av
! to obtain favre averaged quantities

      DO id=1,ngrid
        DO n=1,nx-1
          DO j=1,ny-1
            DO m=1,nz-1
             statsfluct(id)%Q_fluct(n,j,m,8)  = statsfluct(id)%Q_fluct(n,j,m,8)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,9)  = statsfluct(id)%Q_fluct(n,j,m,9)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,10) = statsfluct(id)%Q_fluct(n,j,m,10)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,11) = statsfluct(id)%Q_fluct(n,j,m,11)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,12) = statsfluct(id)%Q_fluct(n,j,m,12)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,13) = statsfluct(id)%Q_fluct(n,j,m,13)/stats(id)%Q_av(n,j,m,1)
!
             statsfluct(id)%Q_fluct(n,j,m,1)  = statsfluct(id)%Q_fluct(n,j,m,1)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,2)  = statsfluct(id)%Q_fluct(n,j,m,2)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,14) = statsfluct(id)%Q_fluct(n,j,m,14)/stats(id)%Q_av(n,j,m,1)
!
             statsfluct(id)%Q_fluct(n,j,m,35) = statsfluct(id)%Q_fluct(n,j,m,35)/stats(id)%Q_av(n,j,m,1)
!
             statsfluct(id)%Q_fluct(n,j,m,36) = statsfluct(id)%Q_fluct(n,j,m,36)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,37) = statsfluct(id)%Q_fluct(n,j,m,37)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,38) = statsfluct(id)%Q_fluct(n,j,m,38)/stats(id)%Q_av(n,j,m,1)
!
             statsfluct(id)%Q_fluct(n,j,m,39) = statsfluct(id)%Q_fluct(n,j,m,39)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,40) = statsfluct(id)%Q_fluct(n,j,m,40)/stats(id)%Q_av(n,j,m,1)
             statsfluct(id)%Q_fluct(n,j,m,41) = statsfluct(id)%Q_fluct(n,j,m,41)/stats(id)%Q_av(n,j,m,1)
           END DO                                                              
         END DO
        END DO
       END DO

        RETURN
      END SUBROUTINE compute_average_statsfluct
!
!////////////////////////////////////////////////////////////////////////////
!
!
     SUBROUTINE compute_velocity_gradients(stats,d,ngrid,mrtr,nmort)
!
      USE domain_definition
      USE mortar_definition
      USE stats_definition
      USE input_data
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (domain)          :: d(ngp)
      TYPE (mortar)          :: mrtr(ngp)
      TYPE (statfl)          :: stats(ngp)
!
      INTEGER                :: ngrid
!
       r=0.0d0
       s=0.0d0
       t=0.0d0
!
       DO id=1,ngrid
         DO  j=1,d(id)%ncg(2)
           DO i=1,d(id)%ncg(1)
             DO k=1,d(id)%ncg(3)
               DO nv=2,4
                 r(id,i,j,k,nv) = stats(id)%Q_av(i,j,k,nv)/stats(id)%Q_av(i,j,k,1)
               ENDDO
             ENDDO
            ENDDO
         ENDDO
       ENDDO
!
       CALL compute_derivative(d,ngrid,mrtr,nmort)
!
       DO nv=1,4
       IF (nv<3) THEN
         IF (nv==1) nv1 = 2
         IF (nv==2) nv1 = 4
         DO id=1,ngrid
           DO i=1,d(id)%ncg(1)
             DO  j=1,d(id)%ncg(2)
               DO k=1,d(id)%ncg(3)
                 stats(id)%Q_xyz(i,j,k,nv) = r(id,i,j,k,nv1)  ! store du/dx and dw/dx
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ELSE
         IF (nv==3) nv1 = 2
         IF (nv==4) nv1 = 4
         DO id=1,ngrid
           DO i=1,d(id)%ncg(2)
             DO  j=1,d(id)%ncg(2)
               DO k=1,d(id)%ncg(2)
                 stats(id)%Q_xyz(i,j,k,nv) = t(id,i,j,k,nv1) ! store du/dz and dw/dz
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDIF
       ENDDO
!
       r=0.0d0
       s=0.0d0
       t=0.0d0
!
        
       RETURN
     END SUBROUTINE compute_velocity_gradients
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE compute_derivative(d,ngrid,mrtr,nmort) 
!
!.......................................................................
!     compute the derivative of a specified variable
!
!     date:  05/07/02                                       
!.......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE User_Data
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      TYPE (domain) :: d(ngp) 
      TYPE (mortar) :: mrtr(ngp) ! nmp
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq):: Qbound 
      LOGICAL :: Is_NaN,Is_INF

!
!      INCLUDE 'mpif.h'
!
!     
!     -------------------------
!     Initializations
!     -------------------------
!
      Qbound = 0.0d0
!
!     -------------------------
!     Interpolate variable to  
!     Lobatto grid
!     -------------------------
!
      DO id=1,ngrid
        CALL interpxa(nx,d(id)%bx,r(id,:,:,:,:),d(id)%ncg(:),d(id)%Qlgg(:,:,:,:),d(id)%nc(:))
        CALL interpya(ny,d(id)%by,r(id,:,:,:,:),d(id)%ncg(:),d(id)%Qglg(:,:,:,:),d(id)%nc(:))
        CALL interpza(nz,d(id)%bz,r(id,:,:,:,:),d(id)%ncg(:),d(id)%Qggl(:,:,:,:),d(id)%nc(:))
      ENDDO

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

      CALL wait_send_q_to_mortar()

!
!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
      DO j = 1,nmort
          IF (myid == nmortedge(j,2) .and. nmortedge(j,1) /= 1) THEN
            jl   = nmortmpi(j,1)
            itag = j
            !Qbound = 0.0d0
            IF (mrtr(jl)%id(2) == 0) THEN ! this is a boundary face
               iface = mrtr(jl)%iface(1)
               idm   = nmortmpi(j,2)
               idml  = ngridmpi(idm,1)
 
!              CALL StatsBC(time,mrtr(jl),d(idml)%bcond(iface),Qbound)
      
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

!     -------------------------
!     compute the mortar fluxes
!     -------------------------
!
      DO j = 1,nmort
          IF (myid == nmortedge(j,2) .and. nmortedge(j,1) == 1) THEN
            jl   = nmortmpi(j,1)
            itag = j
            !Qbound = 0.0d0
            IF (mrtr(jl)%id(2) == 0) THEN ! this is a boundary face
               iface = mrtr(jl)%iface(1)
               idm   = nmortmpi(j,2)
               idml  = ngridmpi(idm,1)
 
!              CALL StatsBC(time,mrtr(jl),d(idml)%bcond(iface),Qbound)
 
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
!     compute the derivatives
!     ----------------------
!
      DO id=1,ngrid
         DO nv = 1,neq
            CALL Gauss_Deriv(d(id)%gmetg,d(id)%Qlgg(:,:,:,nv),d(id)%Qglg(:,:,:,nv), &
                             d(id)%Qggl(:,:,:,nv),r(id,:,:,:,nv),s(id,:,:,:,nv),t(id,:,:,:,nv), &
			     d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz)
         END DO
      ENDDO

!
      RETURN 
      END SUBROUTINE compute_derivative
