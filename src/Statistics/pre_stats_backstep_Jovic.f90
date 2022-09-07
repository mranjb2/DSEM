!> \file
!! Statistics routines for Eulerian solver
!
!> @brief
!> This subroutine computes fluid averaged statistics along a plane.
!!
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
      DOUBLE PRECISION       :: dum(neq+nav),Q(nx,ny,nz,neq+nav),Q1(nz,neq+nav)
      DOUBLE PRECISION       :: QQ(nx,ny,nz,neq+nav),QQ1(nz,neq+nav)
!
!     INCLUDE "mpif.h"
!

      ntot    = neq+nav
          
      ndomy       = 4   
      ndomz       = 4 
      DO 1000 no_block = 1,3
       IF (no_block ==1) THEN
        ndom_offset = 0
        ndomx       = 5 
       ELSEIF (no_block ==2) THEN
        ndom_offset = 80 
        ndomx       = 18 
       ELSEIF (no_block ==3) THEN
        ndom_offset = 368 
        ndomx       = 18
       ENDIF
      
!
!  Sum Q_av in y-direction, and average with nsample and ncg*ndomy
!
      Q  = 0.0d0
      QQ = 0.0d0
      DO idx=1,ndomx
        DO idz=1,ndomz
         id1 = (idx-1)*ndomy + (idz-1)*ndomx*ndomy+1 + ndom_offset
         Q  = 0.0d0
         QQ = 0.0d0
         DO idy=1,ndomy
           id  = id1+idy-1
           idl = ngridmpi(id,1)
           nno = ngridmpi(id,2)
           IF (myid == nno) THEN
             Q(:,:,:,:) = Q(:,:,:,:) + stats(idl)%Q_av(:,:,:,:)
           ENDIF
         ENDDO
!
         CALL MPI_ALLREDUCE(Q,QQ,nx*nz*ntot,MPI_DOUBLE_PRECISION, MPI_SUM, &
                        comm1d, ierr)
!
         DO idy=1,ndomy
           id  = id1+idy-1
           idl = ngridmpi(id,1)
           nno = ngridmpi(id,2)
           IF (myid == nno) THEN
             stats(idl)%Q_av(:,:,:,:) = QQ(:,:,:,:)/(nsample*ncg*ndomy)
           ENDIF
         ENDDO
!
        ENDDO
      ENDDO
!
1000  CONTINUE

        RETURN
      END SUBROUTINE compute_average_stats
! 
!////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes Reynolds and Favre fluctuating averages
!> (turbulence statistics) along a plane.
!!
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
      DOUBLE PRECISION       :: dum(nfluct),Q(nx,ny,nz,nfluct),Q1(nz,nfluct)
      DOUBLE PRECISION       :: QQ(nx,ny,nz,nfluct),QQ1(nz,nfluct)
!
!      INCLUDE "mpif.h"
!

      ntot    = nfluct
          
      ndomy       = 4 
      ndomz       = 4
      DO 1000 no_block = 1,3
       IF (no_block ==1) THEN
        ndom_offset = 0
        ndomx       = 5 
       ELSEIF (no_block ==2) THEN
        ndom_offset = 80 
        ndomx       = 18
       ELSEIF (no_block ==3) THEN
        ndom_offset = 368
        ndomx       = 18 
       ENDIF
!
!  Sum Q_fluct in y-direction, and average with nsample and ncg*ndomy
!
      Q  = 0.0d0
      QQ = 0.0d0
      DO idx=1,ndomx
        DO idz=1,ndomz
         id1 = (idx-1)*ndomy + (idz-1)*ndomx*ndomy+1 + ndom_offset
         Q  = 0.0d0
         QQ = 0.0d0
         DO idy=1,ndomy
           id  = id1+idy-1
           idl = ngridmpi(id,1)
           nno = ngridmpi(id,2)
           IF (myid == nno) THEN
             Q(:,:,:,:) = Q(:,:,:,:) + statsfluct(idl)%Q_fluct(:,:,:,:)
           ENDIF
         ENDDO
!
         CALL MPI_ALLREDUCE(Q,QQ,nx*nz*ntot,MPI_DOUBLE_PRECISION, MPI_SUM, &
                        comm1d, ierr)
!
         DO idy=1,ndomy
           id  = id1+idy-1
           idl = ngridmpi(id,1)
           nno = ngridmpi(id,2)
           IF (myid == nno) THEN
             statsfluct(idl)%Q_fluct(:,:,:,:) = QQ(:,:,:,:)/(nsample*ncg*ndomy)
           ENDIF
         ENDDO
!
        ENDDO
      ENDDO
!
1000  CONTINUE
!
! Once the ensemble averaging is complete
! divide the necessary statsfluct by rho_av
! to obtain favre averaged quantities

      DO id=1,ngrid
        DO n=1,nx-1
         DO m=1,ny-1
         do l=1,nz-1
           statsfluct(id)%Q_fluct(n,m,l,8) = statsfluct(id)%Q_fluct(n,m,l,8)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,9) = statsfluct(id)%Q_fluct(n,m,l,9)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,10) = statsfluct(id)%Q_fluct(n,m,l,10)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,11) = statsfluct(id)%Q_fluct(n,m,l,11)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,12) = statsfluct(id)%Q_fluct(n,m,l,12)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,13) = statsfluct(id)%Q_fluct(n,m,l,13)/stats(id)%Q_av(n,m,l,1)
!
           statsfluct(id)%Q_fluct(n,m,l,1) = statsfluct(id)%Q_fluct(n,m,l,1)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,2) = statsfluct(id)%Q_fluct(n,m,l,2)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,14) = statsfluct(id)%Q_fluct(n,m,l,14)/stats(id)%Q_av(n,m,l,1)
!
           statsfluct(id)%Q_fluct(n,m,l,35) = statsfluct(id)%Q_fluct(n,m,l,35)/stats(id)%Q_av(n,m,l,1)
!
           statsfluct(id)%Q_fluct(n,m,l,36) = statsfluct(id)%Q_fluct(n,m,l,36)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,37) = statsfluct(id)%Q_fluct(n,m,l,37)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,38) = statsfluct(id)%Q_fluct(n,m,l,38)/stats(id)%Q_av(n,m,l,1)
!
           statsfluct(id)%Q_fluct(n,m,l,39) = statsfluct(id)%Q_fluct(n,m,l,39)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,40) = statsfluct(id)%Q_fluct(n,m,l,40)/stats(id)%Q_av(n,m,l,1)
           statsfluct(id)%Q_fluct(n,m,l,41) = statsfluct(id)%Q_fluct(n,m,l,41)/stats(id)%Q_av(n,m,l,1)
          end do                                                           
         END DO
        END DO
       END DO

        RETURN
      END SUBROUTINE compute_average_statsfluct
!
!////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes velocity gradient over domain.
!!
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
             DO k=1,d(id)%ncg(3)
               stats(id)%Q_xyz(i,j,k,nv) = r(id,i,1,k,nv1)  ! store du/dx and dw/dx
             ENDDO
           ENDDO
         ENDDO
       ELSE
         IF (nv==3) nv1 = 2
         IF (nv==4) nv1 = 4
         DO id=1,ngrid
           DO i=1,d(id)%ncg(2)
             DO k=1,d(id)%ncg(2)
               stats(id)%Q_xyz(i,j,k,nv) = t(id,i,1,k,nv1) ! store du/dz and dw/dz
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
!> @brief
!> This subroutine computes the derivative of a specified variable.
!!
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
      TYPE (mortar) :: mrtr(nmp) 
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
