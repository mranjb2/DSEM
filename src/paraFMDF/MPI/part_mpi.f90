!> \file part_mpi.f90
!! Particle communication routines for MPI
!//////////////////////////////////////////////////////////////////////////////////////
!
!     part_mpi.f90 contains
!         SUBROUTINE gather_part_matrix()
!         SUBROUTINE exchange_particles(drop,nprt)
!         SUBROUTINE part_integrate_mpi(k,dt,d,drop,ngrid,nprt,time)
!
!
!     This file takes care of the parallel programming issues of the dispersed phase
!
!     Data created : 08/04/03
!     Programmer   : Guus Jacobs
!
!//////////////////////////////////////////////////////////////////////////////////////
!> @brief Creates matrix of particles on cores
       SUBROUTINE gather_part_matrix()
!
!        Gather the local particle matrix into the global particle matrix and sort
!
         USE mpi_par_part
         USE mpi_par
         USE mpi
!  
         IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
         INTEGER           :: dummat(3*numprocs*npart_mpi)
!
!
!        Gather the local matrix into the global matrix
!
         !write(*,*) '** in gather_part_matrix'
         CALL MPI_ALLGATHER(local_part_mpi(1),npart_mpi*3,MPI_INTEGER,   &
                            dummat,npart_mpi*3,MPI_INTEGER, &
                            comm1d,ierr)
!         write(*,*) 'local_part_mpi',local_part_mpi(:)
!         write(*,*) 'local_part_mpi',dummat(:)

!
!        Sort the global matrix
!
         npart_send_tot = 0
         DO i=0,numprocs-2
           DO j=i+1,numprocs-1
             DO k=1,numprocs*npart_mpi
                npos = (k-1)*3+1
                IF ((i == dummat(npos+1) .AND. j == dummat(npos+2))  .OR. &
                   (i == dummat(npos+2) .AND. j == dummat(npos+1)) ) THEN 
                    npart_send_tot = npart_send_tot + 1
                    global_part_mpi(npart_send_tot,:) = dummat(npos:npos+2) 
                ENDIF               
             ENDDO
           ENDDO
         ENDDO
!
         RETURN
       END SUBROUTINE gather_part_matrix
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Passes particles between cores
       SUBROUTINE exchange_particles(drop,nprt)
!
!     This subroutine exhange particle between nodes
!
         USE mpi_par
         USE mpi_par_part
         USE particle_definition
         USE mpi
!
         IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
         TYPE (particle), DIMENSION(npart)   :: drop          
!
!
         !write(*,*) '** in exchange_particles'
         nreq = 0
!
         DO 1000 i = 1,npart_send_tot
            IF (i==1) THEN
            nmaster_node = global_part_mpi(1,2)
            nslave_node  = global_part_mpi(1,3)
            npart_send  = 1
            npart_recv  = 1
            nsend1       = 0
            nsend2       = 0
            ENDIF
! problem with the cycling, doesn't stop at the right point!!!!!
!
!           write(*,*) 'master node',nmaster_node
            IF ((nmaster_node == global_part_mpi(i,2) .and. nslave_node == global_part_mpi(i,3)) .or.  &
                (nmaster_node == global_part_mpi(i,3) .and. nslave_node == global_part_mpi(i,2)) ) THEN
            IF (nmaster_node == global_part_mpi(i,2) .and. nslave_node == global_part_mpi(i,3)) THEN
                IF (myid == nmaster_node) THEN ! copy particle object into send array
                    ip     = global_part_mpi(i,1)
                    nsend1 = nsend1 + 1
                    npos   = (nsend1-1)*send_p_length
                    
                    send_drop(npos+1:npos+3)    = drop(ip)%Xp
                    send_drop(npos+4:npos+6)    = drop(ip)%Vp
                    send_drop(npos+7)           = drop(ip)%Tp
                    send_drop(npos+8)           = drop(ip)%Mp
                    send_drop(npos+9:npos+11)   = drop(ip)%Xpnm
                    send_drop(npos+12:npos+14)  = drop(ip)%Vp1nm
                    send_drop(npos+15)          = drop(ip)%Tpnm
                    send_drop(npos+16)          = drop(ip)%Mpnm
                    send_drop(npos+17)          = drop(ip)%w
                    send_drop(npos+18)          = drop(ip)%rho
                    send_drop(npos+19)          = drop(ip)%nu_t
                    send_drop(npos+20)          = drop(ip)%omega_m
                    send_drop(npos+21)          = drop(ip)%P
                    send_drop(npos+22)          = drop(ip)%rho_a
                    send_drop(npos+23)          = drop(ip)%react
                    send_drop(npos+24:npos+33)  = drop(ip)%scalar(1:10)
                    drop(ip)%onoff = 0 ! turn droplet off on node that it is leaving !@todo: KILL PARTICLES ONLY HERE
                ELSEIF (myid == nslave_node) THEN
                    nsend2 = nsend2 + 1
                ENDIF
            ENDIF
!
            IF (nmaster_node == global_part_mpi(i,3) .and. nslave_node == global_part_mpi(i,2)) THEN
                IF (myid == nslave_node) THEN ! copy particle object into send array
                    ip     = global_part_mpi(i,1)
                    nsend1 = nsend1 + 1
                    npos   = (nsend1-1)*send_p_length
                    
                    send_drop(npos+1:npos+3)    = drop(ip)%Xp
                    send_drop(npos+4:npos+6)    = drop(ip)%Vp
                    send_drop(npos+7)           = drop(ip)%Tp
                    send_drop(npos+8)           = drop(ip)%Mp
                    send_drop(npos+9:npos+11)   = drop(ip)%Xpnm
                    send_drop(npos+12:npos+14)  = drop(ip)%Vp1nm
                    send_drop(npos+15)          = drop(ip)%Tpnm
                    send_drop(npos+16)          = drop(ip)%Mpnm
                    send_drop(npos+17)          = drop(ip)%w
                    send_drop(npos+18)          = drop(ip)%rho
                    send_drop(npos+19)          = drop(ip)%nu_t
                    send_drop(npos+20)          = drop(ip)%omega_m
                    send_drop(npos+21)          = drop(ip)%P
                    send_drop(npos+22)          = drop(ip)%rho_a
                    send_drop(npos+23)          = drop(ip)%react
                    send_drop(npos+24:npos+33)  = drop(ip)%scalar(1:10)
                    drop(ip)%onoff = 0 ! turn droplet off on node that it is leaving
                ELSEIF (myid == nmaster_node) THEN
                    nsend2 = nsend2 + 1
                ENDIF
            ENDIF
             IF ( (nmaster_node == global_part_mpi(i+1,3) .and. nslave_node == global_part_mpi(i+1,2) ) .or. &
                 (nmaster_node == global_part_mpi(i+1,2) .and. nslave_node == global_part_mpi(i+1,3)) )THEN 
               IF( i< npart_send_tot) CYCLE
             ENDIF
            ENDIF

            IF (myid == nmaster_node) THEN
                    
                    IF (nsend1>0) THEN
                      nreq = nreq + 1
                      npos=(npart_send-1)*send_p_length+1
                      CALL MPI_ISEND (                                          & 
                        send_drop(npos),nsend1*send_p_length, MPI_DOUBLE_PRECISION,   	& 	
                        nslave_node,2*i,comm1d,req(nreq),ierr)
                    ENDIF
                    IF (nsend2>0) THEN
                      nreq = nreq + 1
                      npos=(npart_recv-1)*send_p_length+1
                     CALL MPI_IRECV (                                          & 
                        recv_drop(npos),nsend2*send_p_length, MPI_DOUBLE_PRECISION,   	& 	
                        nslave_node,2*i+1,comm1d,req(nreq),ierr)
                    ENDIF
                    npart_send = npart_send + nsend1
                    npart_recv = npart_recv + nsend2
                 ELSEIF (myid == nslave_node) THEN
                    IF (nsend2>0) THEN
                      nreq = nreq + 1
                      npos=(npart_recv-1)*send_p_length+1
                      CALL MPI_IRECV (                                          & 
                        recv_drop(npos),nsend2*send_p_length, MPI_DOUBLE_PRECISION,   	& 	
                        nmaster_node,2*i,comm1d,req(nreq),ierr)
                    ENDIF
                    IF (nsend1>0) THEN
                      nreq = nreq + 1
                      npos=(npart_send-1)*send_p_length+1
                      CALL MPI_ISEND (                                          & 
                        send_drop(npos),nsend1*send_p_length, MPI_DOUBLE_PRECISION,   	& 	
                        nmaster_node,2*i+1,comm1d,req(nreq),ierr)
                    ENDIF
                    npart_send = npart_send + nsend1
                    npart_recv = npart_recv + nsend2
              
                 ENDIF

                 nsend1 = 0
                 nsend2 = 0
                 nmaster_node = global_part_mpi(i+1,2) 
                 nslave_node  = global_part_mpi(i+1,3) 

1000     CONTINUE
!
         local_part_mpi =0
         global_part_mpi =0
         !write(*,*) '** end exchange_particles'
         RETURN
      END SUBROUTINE exchange_particles
!
!///////////////////////////////////////////////////////////////////////////////
!> Particle integration when passing particles with MPI
   SUBROUTINE part_integrate_mpi(k,dt,d,drop,ngrid,nprt,time)
!
!
      USE size
      USE domain_definition
      USE particle_definition
      USE part_par
      USE input_data
      USE mpi_par_part
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
         npart_send = 0
         nrecv1     = 0 
         DO ip=1,npart
            !IF (ip == nprt) nprt = nprt + 1 ! check this later, if nprt to small to fit all particles, increase it
  	    IF (drop(ip)%onoff == 0) THEN
               nrecv1= nrecv1 + 1
               IF (npart_recv-1 < nrecv1) THEN
                  npart_recv = 0
                  RETURN   ! return to main program if received particles have filled up the array
               ELSEIF (ip > nprt) THEN
                   nprt = nprt + 1
               ELSEIF (ip == npart) THEN
                   write(*,*) 'out of bounds in part_mpi.f90' 
               ENDIF
!
! unpack the receive particle
!
               npos = (nrecv1-1)*send_p_length
               drop(ip)%Xp          = recv_drop(npos+1:npos+3)  
               drop(ip)%Vp          = recv_drop(npos+4:npos+6)   
               drop(ip)%Tp          = recv_drop(npos+7)     
               drop(ip)%Mp          = recv_drop(npos+8)     
               drop(ip)%Xpnm        = recv_drop(npos+9:npos+11)  
               drop(ip)%Vp1nm       = recv_drop(npos+12:npos+14) 
               drop(ip)%Tpnm        = recv_drop(npos+15)    
               drop(ip)%Mpnm        = recv_drop(npos+16)    
               drop(ip)%w           = recv_drop(npos+17)
               drop(ip)%rho         = recv_drop(npos+18)
               drop(ip)%nu_t        = recv_drop(npos+19)
               drop(ip)%omega_m     = recv_drop(npos+20)
               drop(ip)%P           = recv_drop(npos+21)
               drop(ip)%rho_a       = recv_drop(npos+22)
               drop(ip)%react       = recv_drop(npos+23)
               drop(ip)%scalar(1:10)= recv_drop(npos+24:npos+33)
               drop(ip)%onoff = 1
               drop(ip)%Rhofp = 1.0d0
               drop(ip)%Yffp  = 1.0d0
!
! find subdomain and map particle
!
              CALL Part_Interf(drop(ip),ngrid,d)
                
              IF (drop(ip)%onoff ==1) THEN
 	       pgrid = drop(ip)%ngrid
               IF (d(pgrid)%ncg(1) <7) THEN !!!!NOTE THIS SHOULD BE 7
                 CALL spec_Pinterp(d(pgrid),drop(ip),myid,ip)
               ELSE
                 CALL Lagr_Pinterp(d(pgrid),drop(ip),myid,ip)
               END IF
!
               IF (k == 1) THEN
                  CALL integr_onepart_euler(dt,drop(ip))
               ELSEIF (k > 1) THEN
                  nrms = 0
                  IF (nprt>nprtmax) nrms =1
                  CALL integr_onepart_euler(dt,drop(ip))
               END IF
	       CALL part_map(d(pgrid),drop(ip))
               IF (drop(ip)%Xpmap(1) > 1.0d0 .or. drop(ip)%Xpmap(1)<0.0d0 .or.      & 
                   drop(ip)%Xpmap(2) > 1.0d0 .or. drop(ip)%Xpmap(2)<0.0d0 .or.     &
                   drop(ip)%Xpmap(3) > 1.0d0 .or. drop(ip)%Xpmap(3)<0.0d0 ) THEN 
                   CALL part_bc(d,drop,ip,ngrid)
               END IF
              ENDIF
            END IF
         ENDDO
         npart_recv = 0
        
      RETURN
   END SUBROUTINE part_integrate_mpi

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Unpacks particle list from MPI send/recieve
    subroutine partUnpackMPI(k,dt,d,drop,ngrid,nprt,time)
        USE size
        USE domain_definition
        USE particle_definition
        USE part_par
        USE input_data
        USE mpi_par_part
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
          npart_send = 0
          nrecv1     = 0 
          DO ip=1,npart
             !IF (ip == nprt) nprt = nprt + 1 ! check this later, if nprt to small to fit all particles, increase it
    	    IF (drop(ip)%onoff == 0) THEN
                nrecv1= nrecv1 + 1
                IF (npart_recv-1 < nrecv1) THEN
                   npart_recv = 0
                   RETURN   ! return to main program if received particles have filled up the array
                ELSEIF (ip > nprt) THEN
                    nprt = nprt + 1
                ELSEIF (ip == npart) THEN
                    write(*,*) 'out of bounds in part_mpi.f90' 
                ENDIF
    !
    ! unpack the receive particle
    !
                npos = (nrecv1-1)*send_p_length
                drop(ip)%Xp          = recv_drop(npos+1:npos+3)  
                drop(ip)%Vp          = recv_drop(npos+4:npos+6)   
                drop(ip)%Tp          = recv_drop(npos+7)     
                drop(ip)%Mp          = recv_drop(npos+8)     
                drop(ip)%Xpnm        = recv_drop(npos+9:npos+11)  
                drop(ip)%Vp1nm       = recv_drop(npos+12:npos+14) 
                drop(ip)%Tpnm        = recv_drop(npos+15)    
                drop(ip)%Mpnm        = recv_drop(npos+16)    
                drop(ip)%w           = recv_drop(npos+17)
                drop(ip)%rho         = recv_drop(npos+18)
                drop(ip)%nu_t        = recv_drop(npos+19)
                drop(ip)%omega_m     = recv_drop(npos+20)
                drop(ip)%P           = recv_drop(npos+21)
                drop(ip)%rho_a       = recv_drop(npos+22)
                drop(ip)%react       = recv_drop(npos+23)
                drop(ip)%scalar(1:10)      = recv_drop(npos+24:npos+33)
                drop(ip)%onoff = 1
!                drop(ip)%Rhofp = 1.0d0
!                drop(ip)%Yffp  = 1.0d0
    !
    ! find subdomain and map particle
    !
               CALL Part_Interf(drop(ip),ngrid,d)

               IF (drop(ip)%onoff ==1) THEN
    	            pgrid = drop(ip)%ngrid
                    !IF (d(pgrid)%ncg(1) <7) THEN !!!!NOTE THIS SHOULD BE 7
                        CALL spec_Pinterp(d(pgrid),drop(ip),myid,ip)
                    !ELSE
                    !    CALL Lagr_Pinterp(d(pgrid),drop(ip),myid,ip)
                    !END IF
    !
                    !IF (k == 1) THEN
                    !    CALL integr_onepart_euler(dt,drop(ip))
                    !ELSEIF (k > 1) THEN
                    !    nrms = 0
                    !    IF (nprt>nprtmax) nrms =1
                    !    CALL integr_onepart_euler(dt,drop(ip))
                    !END IF
    	            CALL part_map(d(pgrid),drop(ip))
                    IF (drop(ip)%Xpmap(1) > 1.0d0 .or. drop(ip)%Xpmap(1)<0.0d0 .or.      & 
                        drop(ip)%Xpmap(2) > 1.0d0 .or. drop(ip)%Xpmap(2)<0.0d0 .or.     &
                        drop(ip)%Xpmap(3) > 1.0d0 .or. drop(ip)%Xpmap(3)<0.0d0 ) THEN 
                        CALL part_bc(d,drop,ip,ngrid)
                    END IF
                END IF
            END IF
        END DO
        npart_recv = 0

    end subroutine partUnpackMPI


!
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Sends particles between cores for plotting
      SUBROUTINE send_drop_plot(drop,itag,inode)
!
!      sends the droplet array to the inode processor
!
       USE mpi_par
       use mpi_par_part
       USE particle_definition
       USE mpi
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       TYPE( particle ) :: drop(npart)
!
       double precision, allocatable :: send_dropa(:)
!
!
	   if (.not. allocated(send_dropa)) allocate(send_dropa(npart*send_p_length))

       DO ip=1,npart
           npos   = (ip-1)*send_p_length
                 
           send_dropa(npos+1:npos+3)     = drop(ip)%Xp
           send_dropa(npos+4:npos+6)     = drop(ip)%Vp
           send_dropa(npos+7)            = drop(ip)%Tp
           send_dropa(npos+8)            = drop(ip)%Mp
           send_dropa(npos+9:npos+11)    = drop(ip)%Xpnm
           send_dropa(npos+12:npos+14)   = drop(ip)%Vp1nm
           send_dropa(npos+15)           = drop(ip)%Tpnm
           send_dropa(npos+16)           = drop(ip)%Mpnm
           send_dropa(npos+17)           = drop(ip)%onoff
           send_dropa(npos+18)           = drop(ip)%ngrid
           send_dropa(npos+19)           = drop(ip)%nu_t
           send_dropa(npos+20)           = drop(ip)%omega_m
           send_dropa(npos+21)           = drop(ip)%P
           send_dropa(npos+22)           = drop(ip)%rho_a
           send_dropa(npos+23)           = drop(ip)%react
           send_dropa(npos+24:npos+33)   = drop(ip)%scalar(1:10)
        ENDDO
 
        CALL MPI_SEND(send_dropa,npart*send_p_length,MPI_DOUBLE_PRECISION,inode,itag,comm1d,ierr)
	    
	    if (allocated(send_dropa)) deallocate(send_dropa)
!
       RETURN
       END SUBROUTINE send_drop_plot
!
!//////////////////////////////////////////////////////////////////////////////////////
!> Recieves particles between cores for plotting
      SUBROUTINE recv_drop_plot(dropa,itag,inode)
!
!      receives the droplet array to the inode  processor 
!
       USE mpi_par
       use mpi_par_part
       USE particle_definition
       USE mpi
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       TYPE( particle ) :: dropa(npart)
!
       double precision, allocatable :: send_dropa(:)
!
!
	   if(.not. allocated(send_dropa)) allocate(send_dropa(npart*send_p_length))
	
       CALL MPI_RECV(send_dropa,npart*send_p_length,MPI_DOUBLE_PRECISION,inode,itag,comm1d,stat,ierr) 
!
       DO ip=1,npart
           npos = (ip-1)*send_p_length
                  
           dropa(ip)%Xp         = send_dropa(npos+1:npos+3)  
           dropa(ip)%Vp         = send_dropa(npos+4:npos+6)   
           dropa(ip)%Tp         = send_dropa(npos+7)     
           dropa(ip)%Mp         = send_dropa(npos+8)     
           dropa(ip)%Xpnm       = send_dropa(npos+9:npos+11)  
           dropa(ip)%Vp1nm      = send_dropa(npos+12:npos+14) 
           dropa(ip)%Tpnm       = send_dropa(npos+15)    
           dropa(ip)%Mpnm       = send_dropa(npos+16)    
           dropa(ip)%onoff      = send_dropa(npos+17)    
           dropa(ip)%ngrid      = send_dropa(npos+18)
           dropa(ip)%nu_t       = send_dropa(npos+19)
           dropa(ip)%omega_m    = send_dropa(npos+20)
           dropa(ip)%P          = send_dropa(npos+21)
           dropa(ip)%rho_a      = send_dropa(npos+22)
           dropa(ip)%react      = send_dropa(npos+23)
           dropa(ip)%scalar(1:10)     = send_dropa(npos+24:npos+33)    
        ENDDO
      if(allocated(send_dropa)) deallocate(send_dropa)
!
      RETURN
      END SUBROUTINE recv_drop_plot
!




 
