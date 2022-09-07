!////////////////////////////////////////////////////////////////////////////////////////////////
!////////										////////
!////////	restart.f90								////////
!////////										////////
!////////	contains:								////////
!////////										////////
!////////	      SUBROUTINE write_restart(rst_unit,ngrid,d,nmort,mrtr)		////////
!////////	      SUBROUTINE read_restart(rst_unit,ngrid,d,nmort,mrtr)		////////
!////////	      SUBROUTINE restart_mortar_matrices(mrtr)				////////
!////////	      LOGICAL FUNCTION copy_available(idir,idirt,mface,which_copy)	////////
!////////										////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart(rst_unit,ngrid,d,nmort,mrtr,stats,statsfluct,drop,time)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     write out binary file for restarting the computation
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE stats_definition
      USE particle_definition
      USE Order_Matrices
      USE Material_Properties
      USE mpi_par
      USE input_data
      USE User_Data
      USE physics
      USE mpi
!
      IMPLICIT none
!
!      INCLUDE "mpif.h"
!
      TYPE (mortar)    :: mrtr(nmp) 
      TYPE (mortar)    :: mrtra
      TYPE (domain)    :: d(ngp)
      TYPE (domain)    :: da
      TYPE (statfl)    :: stats(ngp)
      TYPE (statfl)    :: statsa
      TYPE (statfluct) :: statsfluct(ngp)
      TYPE (statfluct) :: statsflucta
      TYPE (particle)  :: drop(npart)
!
      INTEGER          :: rst_unit,ngrid,nmort,id,idm,j,mid,nprta,nprtb
      INTEGER          :: idl,nsnode,nrnode,ns,jl,itag,handshake,ip,i
      DOUBLE PRECISION :: time
!
      IF (myid == 0) THEN
        WRITE(rst_unit) ngrid,nmort-nmortd
        IF (defined_resolution) THEN
           CALL write_restart_matrices(rst_unit)
        ELSE
!
!     ---------------------------
!     write out computed matrices
!     ---------------------------
!
        WRITE(rst_unit) num_ordersx,num_ordersy,num_ordersz
        WRITE(rst_unit) order_mapx,order_mapy,order_mapz
        WRITE(rst_unit) bx
        WRITE(rst_unit) by
        WRITE(rst_unit) bz
        WRITE(rst_unit) dmx
        WRITE(rst_unit) dmy
        WRITE(rst_unit) dmz

        WRITE(rst_unit) num_mrtr_matrices
        WRITE(rst_unit) mrtr_mat_map
        WRITE(rst_unit) prolong
        WRITE(rst_unit) restrict
      ENDIF
      ENDIF
!
!     ----------------------------
!     write out domain information
!     ----------------------------
!
      DO id = 1,ngrid
         idl  = ngridmpi(id,1)
         IF (myid == 0 .and. ngridmpi(id,2) == 0) THEN
            da = d(idl)
         ELSEIF (numprocs >1) THEN
            nsnode = ngridmpi(id,2)                      ! node from which one is sending
            nrnode = 0                                   ! receiving node (master in this case)
            itag = id
            DO ns=1,3
              IF (myid /= 0 .and. myid==nsnode) THEN
                 CALL send_domain(d(idl),ns,nrnode,itag)
              ELSEIF (myid == 0 .and. nsnode /= 0) THEN
                 CALL recv_domain(da,ns,nsnode,itag)
              ENDIF
            ENDDO
         ENDIF
!
         IF (myid ==0) THEN
           IF (defined_resolution) THEN 
              CALL write_restart_domain(rst_unit,da)
           ELSE
           WRITE(rst_unit) da%type            ! (8 = hex8,... number of nodes in hex)
           WRITE(rst_unit) da%node
           WRITE(rst_unit) da%orientation 
           WRITE(rst_unit) da%mortar          ! mortar # and side for each face
           WRITE(rst_unit) da%which_surface   ! surface a face is attached to
           WRITE(rst_unit) da%nc,da%ncg       ! grid size
           WRITE(rst_unit) da%Q               ! dependent variables
           WRITE(rst_unit) da%Qlgg            ! dependent variables
           WRITE(rst_unit) da%Qglg            ! dependent variables
           WRITE(rst_unit) da%Qggl            ! dependent variables
           WRITE(rst_unit) da%f               ! fluxes
           WRITE(rst_unit) da%g               ! fluxes
           WRITE(rst_unit) da%h               ! fluxes
           WRITE(rst_unit) da%g_Q             ! rk arrays
!
           WRITE(rst_unit) da%material_id
!
           WRITE(rst_unit) da%cx,da%cxg
           WRITE(rst_unit) da%xg
           WRITE(rst_unit) da%jacob
           WRITE(rst_unit) da%gmet
           WRITE(rst_unit) da%gmetg
!
           WRITE(rst_unit) da%corner
!
           WRITE(rst_unit) da%ibtype
           WRITE(rst_unit) da%bcond
           WRITE(rst_unit) da%domain_type
         ENDIF
         ENDIF

         IF (myid==0) handshake = 1
         CALL MPI_BCAST(handshake,1,MPI_INTEGER,0,comm1d,ierr)
      END DO
!
!     -----------------
!     write out mortars
!     -----------------
!
      DO idm = 1,nmort-nmortd
         jl    = nmortmpi(idm,1)
         id    = nmortmpi(idm,2)
         itag  = idm
         IF (myid == 0 .and. ngridmpi(id,2)==0) THEN
            mrtra = mrtr(jl)
         ELSEIF (numprocs >1) THEN
           nsnode = ngridmpi(id,2)                      ! node from which one is sending
           nrnode = 0                                   ! receiving node (master in this case)
           IF (myid /= 0 .and. myid==nsnode) THEN
               CALL send_mrtr(mrtr(jl),nrnode,itag)
           ELSEIF (myid == 0 .and. nsnode /= 0) THEN
               CALL recv_mrtr(mrtra,nsnode,itag)
           ENDIF
         ENDIF

         IF (myid ==0) THEN
         IF (defined_resolution) THEN
           CALL write_restart_mortar(rst_unit,mrtra)
         ELSE
           WRITE(rst_unit) mrtra%lenmortar    ! dimension of the mortar
           WRITE(rst_unit) mrtra%len          ! dimension of element face (2 directions x 2 sides)
           WRITE(rst_unit) mrtra%id           ! id of the contributing elements
           WRITE(rst_unit) mrtra%iface        ! face # of contributing elements
           WRITE(rst_unit) mrtra%orient       ! orientation of #2 element face (orient*pi/2)
           WRITE(rst_unit) mrtra%sign         ! orientation of normals (+1/-1)
           WRITE(rst_unit) mrtra%nsign         ! orientation of normals (+1/-1)
!
           WRITE(rst_unit) mrtra%n_hat
           WRITE(rst_unit) mrtra%norm         ! geometry arrays along mortar
           WRITE(rst_unit) mrtra%Q            ! solutions along mortar
           WRITE(rst_unit) mrtra%fv            ! solutions along mortar
           WRITE(rst_unit) mrtra%xg           ! positions along mortar
           WRITE(rst_unit) mrtra%conforming   ! mortar is conforming or not (2 directions x 2 sides)
         ENDIF
         ENDIF

         IF (myid==0) handshake = 1
         CALL MPI_BCAST(handshake,1,MPI_INTEGER,0,comm1d,ierr)
      END DO
!
!     -----------------------------
!     write out material properties
!     -----------------------------
!
      IF (myid == 0) THEN
        WRITE(rst_unit) num_prop
        DO mid = 1,max_materials
           WRITE(rst_unit) (material_property(j,mid),j=1,num_prop)
        END DO
      ENDIF
!
!     -----------------------------
!     write case specific stuff
!     -----------------------------
!
      CALL write_restart_case(rst_unit,time)
!
!     ----------------------------
!     write out stats information
!     ----------------------------
!
      IF  (statistics) THEN
!           write(*,*) 'write stats'
!            IF (average) THEN
!               CALL compute_average_stats(stats,d(1)%ncg(2),ngridl)
!            ELSEIF (rms) THEN
!               CALL compute_average_statsfluct(statsfluct,stats,d(1)%ncg(2),ngridl)
!            ENDIF
! for printing the result only, one has to turn this off to start of continue the rms
!         CALL compute_average_statsfluct(statsfluct,stats,d(1)%ncg(2),ngridl)
!----------------------------------------------------------------------------        
!         write(*,*) 'first stats'
       IF (myid==0) write(rst_unit) nsample
      DO id = 1,ngrid
         idl  = ngridmpi(id,1)
         IF (myid == 0 .and. ngridmpi(id,2) == 0) THEN
            statsa = stats(idl)
            IF (rms) statsflucta = statsfluct(idl)
         ELSEIF (numprocs >1) THEN
            nsnode = ngridmpi(id,2)                      ! node from which one is sending
            nrnode = 0                                   ! receiving node (master in this case)
            itag = id
            DO ns   = 2,3
            IF (myid /= 0 .and. myid==nsnode) THEN
               IF (ns==3) CALL send_stats_plot(stats(idl),ns,nrnode,itag)
               IF (rms .and. ns==2)   CALL send_statsfluct_plot(statsfluct(idl),ns,nrnode,itag)
            ELSEIF (myid == 0 .and. nsnode /= 0) THEN
               IF (ns==3) CALL recv_stats_plot(statsa,ns,nsnode,itag)
               IF (rms .and. ns==2)   CALL recv_statsfluct_plot(statsflucta,ns,nsnode,itag)
            ENDIF
            ENDDO
         ENDIF
     
         IF (myid==0) THEN
           write(rst_unit) statsa%Q_av
           IF (rms) THEN 
             write(rst_unit) statsflucta%Q_fluct
           ENDIF
           
           CALL zero_stats(statsa,statsflucta)
         ENDIF
         IF (myid==0) handshake = 1
         CALL MPI_BCAST(handshake,1,MPI_INTEGER,0,comm1d,ierr)

      ENDDO
      ENDIF
!
!     ----------------------------
!     write out droplet information
!     ----------------------------
!
     IF (drops) THEN
!
!  Determine the total number of particles per processor, reduce later
!

      nprta = 0
      DO i=0,numprocs-1
        IF (myid==i) THEN
          DO ip=1,npart
            IF (drop(ip)%onoff==1) THEN
              nprta=nprta+1
             ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!   Send the droplets from other processors to the root processer and write
!   to file
!
      IF (numprocs>1) THEN
!        convert local to global number
         DO ip=1,nprta
          DO id=1,ngrid
            IF (myid==ngridmpi(id,2) .and. drop(ip)%ngrid==ngridmpi(id,1)) THEN
              drop(ip)%ngrid = id
            END IF 
           END DO
         END DO 
         CALL MPI_REDUCE(nprta,nprtb,1,MPI_INTEGER, MPI_SUM,0, &
              comm1d, ierr)

        DO i=0,numprocs-1
           IF (i==0) THEN
             IF (myid==0) THEN
              WRITE(rst_unit) nprtb
              CALL write_restart_drop(drop,rst_unit)
            ENDIF
           ELSE
              IF (myid==i) THEN
                CALL send_drop_plot(drop,i,0)
              ELSEIF (myid==0) THEN
               CALL recv_drop_plot(drop,i,i)
               CALL write_restart_drop(drop,rst_unit)
              ENDIF
           ENDIF
        ENDDO
       ELSE   ! one processor
         WRITE(rst_unit) nprta
         CALL write_restart_drop(drop,rst_unit)
       ENDIF  ! if loop on numprocessors
 
      ENDIF
!
      RETURN
      END SUBROUTINE write_restart
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE read_restart(rst_unit,ngrid,d,nmort,mrtr,stats,statsfluct,drop,nprt,time)
!
!......................................................................
!     date: 11/24/98
!
!     stac3m version
!
!     read in binary file for restarting the computation
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE stats_definition
      USE particle_definition
      USE Order_Matrices
      USE Material_Properties
      USE User_Data
      USE input_data
      USE mpi_par
      USE constants
      USE physics
      USE mpi
!
      IMPLICIT none

!      INCLUDE "mpif.h"
!
      TYPE (mortar)     :: mrtr(nmp) 
      TYPE (mortar)     :: mrtra
      TYPE (domain)     :: d(ngp)
      TYPE (domain)     :: da
      TYPE (statfl)     :: stats(ngp)
      TYPE (statfl)     :: statsa
      TYPE (statfluct)  :: statsfluct(ngp)
      TYPE (statfluct)  :: statsflucta
      TYPE (particle)   :: drop(npart)
      TYPE (particle)   :: dropa

      INTEGER                               :: rst_unit,ngrid,nmort,id,idm,k,j,mid,epsa,inode,i
      INTEGER                               :: nsnode,nrnode,ip,ns,idl,nprt,nprta
      INTEGER                               :: ndumdouble,ndummort,id1,id2,iface1,iface2
      INTEGER, DIMENSION(numprocs)          :: idp
      INTEGER                               :: itag,nno
      INTEGER                               :: ngriddum
      DOUBLE PRECISION                      :: xmid,ymid,zmid,time
!
!     -----------------
!     read in run sizes
!     -----------------
!
      IF (myid == 0 ) THEN
        READ(rst_unit) ngrid,nmort
      ENDIF
!
      CALL MPI_BCAST(ngrid,1,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(nmort,1,MPI_INTEGER,0,comm1d,ierr)
 
!
!   Set stuff for domain decomposition for MPI (equidistant
!   parts in x-direction)
!
      idp(:) = 0

!
!     -------------------------
!     read in computed matrices
!     -------------------------
!
      IF (myid == 0) THEN
        READ(rst_unit) num_ordersx,num_ordersy,num_ordersz
        READ(rst_unit) order_mapx,order_mapy,order_mapz
        READ(rst_unit) bx
        READ(rst_unit) by
        READ(rst_unit) bz
        READ(rst_unit) dmx
        READ(rst_unit) dmy
        READ(rst_unit) dmz

        READ(rst_unit) num_mrtr_matrices
        READ(rst_unit) mrtr_mat_map
        READ(rst_unit) prolong
        READ(rst_unit) restrict
      ENDIF

      CALL MPI_BCAST(num_ordersx,1,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(num_ordersy,1,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(num_ordersz,1,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(order_mapx,max_orders,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(order_mapy,max_orders,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(order_mapz,max_orders,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(bx ,nx*nx*max_orders,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(by ,ny*ny*max_orders,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(bz ,nz*nz*max_orders,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(dmx,nx*nx*max_orders,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(dmy,ny*ny*max_orders,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(dmz,nz*nz*max_orders,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(num_mrtr_matrices,1,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(mrtr_mat_map,2*max_mrtr_mat,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(prolong ,nmax*nmax*max_mrtr_mat,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      CALL MPI_BCAST(restrict ,nmax*nmax*max_mrtr_mat,MPI_DOUBLE_PRECISION,0,comm1d,ierr)

!
!     --------------------------
!     read in domain information
!     --------------------------
!
      IF (myid == 0) THEN
        DO id = 1,ngrid
          CALL zero_domain(da)
          READ(rst_unit) da%type            ! (8 = hex8,... number of nodes in hex)
          READ(rst_unit) da%node
          READ(rst_unit) da%orientation 
          READ(rst_unit) da%mortar          ! mortar # and side for each face
          READ(rst_unit) da%which_surface   ! surface a face is attached to

          READ(rst_unit) da%nc,da%ncg    ! grid size

          READ(rst_unit) da%Q               ! dependent variables
          READ(rst_unit) da%Qlgg            ! dependent variables
          READ(rst_unit) da%Qglg            ! dependent variables
          READ(rst_unit) da%Qggl            ! dependent variables
          READ(rst_unit) da%f               ! fluxes
          READ(rst_unit) da%g               ! fluxes
          READ(rst_unit) da%h               ! fluxes
          READ(rst_unit) da%g_Q             ! rk arrays
!
          READ(rst_unit) da%material_id
!
          READ(rst_unit) da%cx,da%cxg
          READ(rst_unit) da%xg
          READ(rst_unit) da%jacob
          READ(rst_unit) da%gmet
          READ(rst_unit) da%gmetg
!
          READ(rst_unit) da%corner
!
          READ(rst_unit) da%ibtype
          READ(rst_unit) da%bcond
          READ(rst_unit) da%domain_type
!
!   Compute approximate middle point of element
!   for use of grid partitioning and the node no. (nno)
!
          xmid = SUM(da%corner(1,:))/8.0d0
          ymid = SUM(da%corner(2,:))/8.0d0
          zmid = SUM(da%corner(3,:))/8.0d0
!
!   Compute approximate middle point of element
!   for use of grid partitioning and the node no. (nno)
!       
          CALL find_processor(xmid,ymid,zmid,id, nno, ngrid)

          IF (nno==0) THEN        ! this is the master node
              inode                = 0
              idp(inode+1)         = idp(inode+1) + 1
              ngridmpi(id,1)       = idp(inode+1)
              ngridmpi(id,2)       = inode
              d(idp(inode+1))      = da
           ELSE                   ! this a slave node
              inode                = nno
              idp(inode+1)         = idp(inode+1) + 1
              ngridmpi(id,1)       = idp(inode+1)
              ngridmpi(id,2)       = inode
           ENDIF
       ENDDO
     ENDIF
!
       CALL MPI_BCAST(ngridmpi,ngp*numprocs*2,MPI_INTEGER,0,comm1d,ierr)
       CALL MPI_BCAST(idp,numprocs,MPI_INTEGER,0,comm1d,ierr)

      ngridl = idp(myid+1)
      IF (myid ==0) THEN
        OPEN(unit=22,file="grid_distribution.dat")
        DO id=1,ngrid
          write(22,*) 'ngridl ',ngridmpi(id,:)
        ENDDO
        CLOSE(22)
      ENDIF

      IF (myid==0) THEN
      DO idm = 1,ngrid
         DO i=1,24
           BACKSPACE(rst_unit)
         ENDDO
      ENDDO
      ENDIF
      
!
        DO id = 1,ngrid
          IF (myid==0) THEN
          CALL zero_domain(da)
          READ(rst_unit) da%type            ! (8 = hex8,... number of nodes in hex)
          READ(rst_unit) da%node
          READ(rst_unit) da%orientation 
          READ(rst_unit) da%mortar          ! mortar # and side for each face
          READ(rst_unit) da%which_surface   ! surface a face is attached to

          READ(rst_unit) da%nc,da%ncg    ! grid size

          READ(rst_unit) da%Q               ! dependent variables
          READ(rst_unit) da%Qlgg            ! dependent variables
          READ(rst_unit) da%Qglg            ! dependent variables
          READ(rst_unit) da%Qggl            ! dependent variables
          READ(rst_unit) da%f               ! fluxes
          READ(rst_unit) da%g               ! fluxes
          READ(rst_unit) da%h               ! fluxes
          READ(rst_unit) da%g_Q             ! rk arrays
!
          READ(rst_unit) da%material_id
!
          READ(rst_unit) da%cx,da%cxg
          READ(rst_unit) da%xg
          READ(rst_unit) da%jacob
          READ(rst_unit) da%gmet
          READ(rst_unit) da%gmetg
!
          READ(rst_unit) da%corner
!
          READ(rst_unit) da%ibtype
          READ(rst_unit) da%bcond
          READ(rst_unit) da%domain_type
       ENDIF
       IF (numprocs >1) THEN
         DO ip = 1,numprocs-1
           nsnode = 0
           nrnode = ip
           itag   = id
           DO ns=1,3
             IF (myid == 0 .and. ngridmpi(id,2) == nrnode) THEN
               CALL send_domain(da,ns,nrnode,itag)
             ELSEIF (myid == nrnode .and. ngridmpi(id,2) == nrnode) THEN
               idl = ngridmpi(id,1)
               CALL recv_domain(d(idl),ns,nsnode,itag)
             ENDIF
           ENDDO
         ENDDO
       ENDIF
!
      END DO ! end grid loop
!
      ngridl = idp(myid+1)
!
!     ---------------
!     read in mortars
!     ---------------
!
      IF (myid ==0) THEN
      idp(:) = 0
      DO idm = 1,nmort
        j = idm
          CALL zero_mortar(mrtra)
          READ(rst_unit) mrtra%lenmortar ! dimension of the mortar
          READ(rst_unit) mrtra%len     ! dimension of element face (2 directions x 2 sides)
          READ(rst_unit) mrtra%id        ! id of the contributing elements
          READ(rst_unit) mrtra%iface     ! face # of contributing elements
          READ(rst_unit) mrtra%orient       ! orientation of #2 element face (orient*pi/2)
          READ(rst_unit) mrtra%sign         ! orientation of normals (+1/-1)
          READ(rst_unit) mrtra%nsign         ! orientation of normals (+1/-1)
!
          READ(rst_unit) mrtra%n_hat
          READ(rst_unit) mrtra%norm         ! geometry arrays along mortar
          READ(rst_unit) mrtra%Q            ! solutions along mortar
          READ(rst_unit) mrtra%fv            ! solutions along mortar
          READ(rst_unit) mrtra%xg           ! positions along mortar
          READ(rst_unit) mrtra%conforming   ! mortar is conforming or not (2 directions x 2 sides)
          id = mrtra%id(1)
          IF (ngridmpi(id,2) == 0) THEN
            idp(1)         = idp(1) + 1
            nmortmpi(j,1)  = idp(1)
            nmortmpi(j,2)  = id
            nmortmpi(j,3)  = mrtra%id(2)
            nslavempi(j,1) = mrtra%iface(2)
            nslavempi(j,2) = mrtra%len(1,2)
            nslavempi(j,3) = mrtra%len(2,2)
            mrtr(idp(1))   = mrtra
          ELSE
            inode          = ngridmpi(id,2)
            idp(inode+1)   = idp(inode+1)+1
            nmortmpi(j,1)  = idp(inode+1)
            nmortmpi(j,2)  = id
            nmortmpi(j,3)  = mrtra%id(2)
            nslavempi(j,1) = mrtra%iface(2)
            nslavempi(j,2) = mrtra%len(1,2)
            nslavempi(j,3) = mrtra%len(2,2)
          ENDIF
        ENDDO
      ENDIF
!
!  Rewind file
!
      IF (myid==0) THEN
      DO idm = 1,nmort
         DO i=1,13
           BACKSPACE(rst_unit)
         ENDDO
      ENDDO
      ENDIF
! 
      CALL MPI_BCAST(nmortmpi,nmp*numprocs*3,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(idp,numprocs,MPI_INTEGER,0,comm1d,ierr)
!
      nmortl = idp(myid +1)
!
      CALL preproc_mpi(nmort) 
!
      ndumdouble = 0
      DO idm = 1,nmort-nmortd
        j = idm
        IF (myid == 0) THEN
          CALL zero_mortar(mrtra)
          READ(rst_unit) mrtra%lenmortar 
          READ(rst_unit) mrtra%len     
          READ(rst_unit) mrtra%id     
          READ(rst_unit) mrtra%iface 
          READ(rst_unit) mrtra%orient 
          READ(rst_unit) mrtra%sign   
          READ(rst_unit) mrtra%nsign  
!
          READ(rst_unit) mrtra%n_hat
          READ(rst_unit) mrtra%norm  
          READ(rst_unit) mrtra%Q            
          READ(rst_unit) mrtra%fv           
          READ(rst_unit) mrtra%xg          
          READ(rst_unit) mrtra%conforming  
      
        ENDIF

        IF ( nmortedge(j,1) == 1) THEN 
          ndumdouble = ndumdouble + 1
          ndummort   = nmort-nmortd+ ndumdouble
        ENDIF

        IF (numprocs >1) THEN
          DO ip =1,numprocs-1
            id     = nmortmpi(idm,2)
            nsnode = 0
            nrnode = ip
            itag   = j
            IF (myid == 0 .and.  ngridmpi(id,2)== nrnode) THEN
              CALL send_mrtr(mrtra,nrnode,itag)
            ELSEIF (myid == nrnode .and. ngridmpi(id,2) == nrnode) THEN
              idl  = nmortmpi(idm,1)
              CALL recv_mrtr(mrtr(idl),nsnode,itag)
              id1  = mrtr(idl)%id(1)
              id2  = ngridmpi(id1,1)
              IF  (nmortedge(j,1) == 1) ngridedge(id2) = 1
            ENDIF
          ENDDO
        ENDIF

        IF ( nmortedge(j,1) == 1) THEN 
         IF (myid == 0 .and. nmortedge(j,2) ==0 ) THEN
           id1                  = mrtra%id(1)
           idl                  = ngridmpi(id1,1)
           ngridedge(idl)   = 1
         ENDIF
         IF (myid == 0 .and.  nmortedge(j,3) == 0) THEN           ! mortar is copied on the master node
              !idp(nrnode+1) = idp(nrnode+1) + 1
              idp(1) = idp(1)+1
              nmortl               = nmortl+1
              nmortmpi(ndummort,1) = nmortl
              nslavempi(ndummort,1) = mrtra%iface(1)
              nslavempi(ndummort,2) = mrtra%len(2,2)
              nslavempi(ndummort,3) = mrtra%len(1,2)
              mrtr(nmortl)         = mrtra
              id1                  = mrtr(nmortl)%id(1)
              id2                  = mrtr(nmortl)%id(2)
              iface1               = mrtr(nmortl)%iface(1)
              iface2               = mrtr(nmortl)%iface(2)
              mrtr(nmortl)%id(1)   = id2
              mrtr(nmortl)%id(2)   = id1
              mrtr(nmortl)%iface(1)= iface2
              mrtr(nmortl)%iface(2)= iface1
              idl                  = ngridmpi(id2,1)
              ngridedge(idl)       = 1
         ELSE                    ! mortar is copied on a slave node
          DO ip =1,numprocs-1
            id     = nmortmpi(ndummort,2)
            nsnode = 0
            nrnode = ip
            itag   = ndummort
            IF (myid == 0 .and.  ngridmpi(id,2)== nrnode) THEN
              idp(nrnode+1) = idp(nrnode+1) + 1
              nmortmpi(ndummort,1) = idp(nrnode+1)
              CALL send_mrtr(mrtra,nrnode,itag)
              nslavempi(ndummort,1) = mrtra%iface(1)
              nslavempi(ndummort,2) = mrtra%len(2,2)
              nslavempi(ndummort,3) = mrtra%len(1,2)
            ELSEIF (myid == nrnode .and. ngridmpi(id,2) == nrnode) THEN
              nmortl               = nmortl+1
              CALL recv_mrtr(mrtr(nmortl),nsnode,itag)
              id1                  = mrtr(nmortl)%id(1)
              id2                  = mrtr(nmortl)%id(2)
              iface1               = mrtr(nmortl)%iface(1)
              iface2               = mrtr(nmortl)%iface(2)
              mrtr(nmortl)%id(1)   = id2
              mrtr(nmortl)%id(2)   = id1
              mrtr(nmortl)%iface(1)= iface2
              mrtr(nmortl)%iface(2)= iface1
              idl                  = ngridmpi(id2,1)
              ngridedge(idl)       = 1
            ENDIF
          ENDDO
         ENDIF
        ENDIF
!
      END DO ! end mortar loop
!
      ngriddum = 0
      DO id=1,ngridl
        IF (ngridedge(id) ==1) ngriddum = ngriddum +1
      ENDDO

      CALL MPI_BCAST(nmortmpi,nmp*numprocs*3,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(nslavempi,nmp*8*3,MPI_INTEGER,0,comm1d,ierr)
!
      IF (myid ==0) THEN
        OPEN(unit=21,file="mortar_distribution.dat")
        write(21,*) '#grid,node,ngridedge '
        DO id=1,ngrid
          write(21,*) id,ngridmpi(id,2),ngridedge(id)
        ENDDO
        write(21,*) '#mort,nmortmpi(1),(2),(2) '
        DO i=1,nmort
          write(21,*) i,nmortmpi(i,1),nmortmpi(i,2),nmortmpi(i,3)
        ENDDO 
        CLOSE(21)
      ENDIF
!     -----------------------------
!     read in material properties
!     -----------------------------

      IF (myid ==0 ) THEN
        READ(rst_unit) num_prop
        DO mid = 1,max_materials
           READ(rst_unit) (material_property(j,mid),j=1,num_prop)
        END DO
      ENDIF
!
!     -----------------------------
!     read case specific stuff
!     -----------------------------
!
      CALL read_restart_case(rst_unit,time)
!
!-------------------------------------------------
!  Read statistics information
!-------------------------------------------------
!! read the average data from previous simulation
!------------------------------------------------------
      IF  (statistics .and. average .and. rmscont) THEN
        IF (myid==0) THEN
           READ(rst_unit) nsample
!            IF (.not.rmscont) nsample = 0
        ENDIF
! send the number of samples to every node, modified by Dongru
        CALL MPI_BCAST(nsample,1,MPI_INTEGER,0,comm1d,ierr)
        DO id = 1,ngrid
          IF (myid==0) THEN
           CALL zero_stats(statsa,statsflucta)
           READ(rst_unit) statsa%Q_av    
!            IF (rmscont) THEN 
!              READ(rst_unit) statsflucta%Q_fluct
!            ENDIF
          ENDIF
         IF (numprocs >1) THEN
          IF (ngridmpi(id,2) == 0 .and. myid == 0) THEN
            idl = ngridmpi(id,1)
            stats(idl)%Q_av =statsa%Q_av
!             IF (rmscont) statsfluct(idl) =statsflucta
          ELSE
           DO ip = 1,numprocs-1
             nsnode = 0
             nrnode = ip
             itag   = id
             DO ns=1,3
               IF (myid == 0 .and. ngridmpi(id,2) == nrnode) THEN
                 CALL send_stats_plot(statsa,ns,nrnode,itag)
!                  IF (ns==2 .and. rmscont) CALL send_statsfluct_plot(statsflucta,ns,nrnode,itag)
               ELSEIF (myid == nrnode .and. ngridmpi(id,2) == nrnode) THEN
                 idl = ngridmpi(id,1)
                 CALL recv_stats_plot(stats(idl),ns,nsnode,itag)
!                  IF (ns==2 .and. rmscont) CALL recv_statsfluct_plot(statsfluct(idl),ns,nsnode,itag)
               ENDIF
             ENDDO
           ENDDO
          ENDIF
         ELSE
            idl = ngridmpi(id,1)
            stats(idl)%Q_av =statsa%Q_av
!             IF (rmscont) statsfluct(idl) =statsflucta
         ENDIF
!
      END DO ! end grid loop
! 
!
     ENDIF
     
!------------------------------------------------------------------------------------
!!! read the rms data from the previous simulation
!------------------------------------------------------------------------------------
!
      IF  (statistics .and. rms) THEN
        IF (myid==0) THEN
           READ(rst_unit) nsample
           IF (.not.rmscont) nsample = 0
        ENDIF
! send the number of samples to every node, modified by Dongru
        CALL MPI_BCAST(nsample,1,MPI_INTEGER,0,comm1d,ierr)
        DO id = 1,ngrid
          IF (myid==0) THEN
           CALL zero_stats(statsa,statsflucta)
           READ(rst_unit) statsa%Q_av    
           IF (rmscont) THEN 
             READ(rst_unit) statsflucta%Q_fluct
           ENDIF
          ENDIF
         IF (numprocs >1) THEN
          IF (ngridmpi(id,2) == 0 .and. myid == 0) THEN
            idl = ngridmpi(id,1)
            stats(idl)%Q_av =statsa%Q_av
            IF (rmscont) statsfluct(idl) =statsflucta
          ELSE
           DO ip = 1,numprocs-1
             nsnode = 0
             nrnode = ip
             itag   = id
             DO ns=1,3
               IF (myid == 0 .and. ngridmpi(id,2) == nrnode) THEN
                 CALL send_stats_plot(statsa,ns,nrnode,itag)
                 IF (ns==2 .and. rmscont) CALL send_statsfluct_plot(statsflucta,ns,nrnode,itag)
               ELSEIF (myid == nrnode .and. ngridmpi(id,2) == nrnode) THEN
                 idl = ngridmpi(id,1)
                 CALL recv_stats_plot(stats(idl),ns,nsnode,itag)
                 IF (ns==2 .and. rmscont) CALL recv_statsfluct_plot(statsfluct(idl),ns,nsnode,itag)
               ENDIF
             ENDDO
           ENDDO
          ENDIF
         ELSE
            idl = ngridmpi(id,1)
            stats(idl)%Q_av =statsa%Q_av
            IF (rmscont) statsfluct(idl) =statsflucta
         ENDIF
!
      END DO ! end grid loop
! 
!
     ENDIF
!
!     -----------------------------
!     read in droplets
!     -----------------------------

      IF (drops .and. (.not. restfldrop)) THEN
        DO 1000 ip=0,numprocs-1
         id    = numprocs-1-ip
         IF (myid==0) THEN
          nprt = 0
          DO i=1,npart
            CALL zero_drops(drop(i))
          ENDDO
          READ(rst_unit) nprta
          DO j=1,nprta
            CALL zero_drops(dropa)
            READ(rst_unit) dropa%Xp
            READ(rst_unit) dropa%Vp
            READ(rst_unit) dropa%Tp
            READ(rst_unit) dropa%Mp
            READ(rst_unit) dropa%Xpnm
            READ(rst_unit) dropa%Vp1nm
            READ(rst_unit) dropa%Tpnm
            READ(rst_unit) dropa%Mpnm
            READ(rst_unit) dropa%ngrid
            READ(rst_unit) dropa%onoff

!
            IF (dropa%onoff == 1) THEN
              CALL find_processor(dropa%Xp(1),dropa%Xp(2),dropa%Xp(3),dropa%ngrid, nno, ngrid)
!              IF (ngridmpi(dropa%ngrid,2)==id) THEN
              IF (nno==id) THEN
               nprt             = nprt + 1
               drop(nprt)       = dropa
               drop(nprt)%onoff = 1
               drop(nprt)%Rhofp = 1.0d0
               drop(nprt)%Yffp  = 1.0d0
               drop(nprt)%ngrid  = ngridmpi(drop(nprt)%ngrid,1)
              ENDIF
            ENDIF
           ENDDO
!
!
           IF (id/=0) CALL send_drop_plot(drop,id,id)

           DO j=1,nprta
             DO i=1,10
               BACKSPACE(rst_unit)
             ENDDO
           ENDDO
           BACKSPACE(rst_unit)
         ELSEIF (myid==id) THEN
!
           IF (id/=0) CALL recv_drop_plot(drop,id,0)
!
           nprt = 0
           DO j=1,npart
             IF (drop(j)%onoff ==1) THEN
               nprt = nprt + 1
               drop(nprt)%Rhofp = 1.0d0
               drop(nprt)%Yffp  = 1.0d0
               drop(nprt)%ngrid  = ngridmpi(drop(nprt)%ngrid,1)
             ENDIF
           ENDDO
         ENDIF
1000    CONTINUE
      ENDIF

      
      IF (.NOT. defined_resolution ) THEN
!
!
!     --------------
!     reset pointers
!     --------------
!
      DO id = 1,ngridl
         DO k = 1,num_ordersx
            IF ( order_mapx(k) == d(id)%nc(1) )     THEN
               d(id)%dmx => dmx(:,:,k)
               d(id)%bx  => bx (:,:,k)
               EXIT
            END IF
         END DO
      END DO
      DO id = 1,ngridl
         DO k = 1,num_ordersy
            IF ( order_mapy(k) == d(id)%nc(2) )     THEN
               d(id)%dmy => dmy(:,:,k)
               d(id)%by  => by (:,:,k)
               EXIT
            END IF
         END DO
      END DO
      DO id = 1,ngridl
         DO k = 1,num_ordersz
            IF ( order_mapz(k) == d(id)%nc(3) )     THEN
               d(id)%dmz => dmz(:,:,k)
               d(id)%bz  => bz (:,:,k)
               EXIT
            END IF
         END DO
      END DO

      DO j = 1,nmortl
         CALL restart_mortar_matrices(mrtr(j))
      END DO

      ENDIF
!
  
      IF  (statistics .and. rms) THEN
        CALL compute_velocity_gradients(stats,d,ngridl,mrtr,nmort)
      ENDIF
!
      RETURN
      END SUBROUTINE read_restart
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE restart_mortar_matrices(mrtr) 
!
!......................................................................
!     date: 12/11/98
!     routines called: none                                           
!
!     sets up arrays for mortar projections. to save on storage,
!     only unique prolongation and restriction arrays are computed
!     and stored in an array.
!......................................................................
!
      USE mortar_definition
      USE Order_Matrices
      USE File_Units
      USE Edge_Mappings

      IMPLICIT none

      TYPE (mortar)             :: mrtr
      INTEGER                   :: k,listloc,n,m,l,jj
      INTEGER                   :: mface,idir,idirs,idirm
!
!     --------------------
!     side 1 (master) face
!     --------------------
!
      mface = 1
      idir = 1
      idirm = 1
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirm,mface,listloc) )     THEN
            mrtr%prolong1_dir1 => prolong(:,:,listloc)
            mrtr%restrict1_dir1 => restrict(:,:,listloc)
         ELSE
            WRITE(err_unit,*) 'mortar matrix not available'
            STOP 'mortar problem on restart'
         END IF
      END IF

      idir = 2
      idirm = 2
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirm,mface,listloc) )     THEN
            mrtr%prolong1_dir2 => prolong(:,:,listloc)
            mrtr%restrict1_dir2 => restrict(:,:,listloc)
         ELSE
            WRITE(err_unit,*) 'mortar matrix not available'
            STOP 'mortar problem on restart'
         END IF
      END IF
!
!     ------------------------------------------------
!     side 2 (slave) face (may be rotated wrt mortar)
!     ------------------------------------------------
!
      mface = 2
      idir = 1
      SELECT CASE (mrtr%orient)
         CASE (DFLT,B1F2,B1B2,F1B2)
            idirs = 1
         CASE (F2F1,B2F1,B2B1,F2B1)
            idirs = 2
      END SELECT
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirs,mface,listloc) )     THEN
            mrtr%prolong2_dir1 => prolong(:,:,listloc)
            mrtr%restrict2_dir1 => restrict(:,:,listloc)
         ELSE
            WRITE(err_unit,*) 'mortar matrix not available'
            STOP 'mortar problem on restart'
         END IF
      END IF

      idir = 2
      SELECT CASE (mrtr%orient)
         CASE (DFLT,B1F2,B1B2,F1B2)
            idirs = 2
         CASE (F2F1,B2F1,B2B1,F2B1)
            idirs = 1
      END SELECT
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirs,mface,listloc) )     THEN
            mrtr%prolong2_dir2 => prolong(:,:,listloc)
            mrtr%restrict2_dir2 => restrict(:,:,listloc)
         ELSE
            WRITE(err_unit,*) 'mortar matrix not available'
            STOP 'mortar problem on restart'
         END IF
      END IF
!
      RETURN
!     --------
      CONTAINS
!     --------

         LOGICAL FUNCTION copy_available(idir,idirt,mface,which_copy)
!
!        Go through the list of arrays and see if we've already computed
!        this one.
!
         INTEGER :: which_copy,k,idir,mface,idirt
!
         copy_available = .FALSE.
         DO k = 1,num_mrtr_matrices
            IF ( mrtr_mat_map(1,k) == mrtr%len(idirt,mface) .AND. &
                 mrtr_mat_map(2,k) == mrtr%lenmortar(idir) )  THEN
               copy_available = .TRUE.
               which_copy = k
               EXIT
            END IF
         END DO
         RETURN
         END FUNCTION copy_available
!
      END SUBROUTINE restart_mortar_matrices
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart_matrices(rst_unit)
!
!......................................................................
!     date: 01/22/02
!......................................................................
!
      USE Order_Matrices

      IMPLICIT none
!
      INTEGER       :: rst_unit
      DOUBLE PRECISION :: bxdum(nx+nref,nx+nref,max_orders)
      DOUBLE PRECISION :: prolongdum(nmax+nref,nmax+nref,max_mrtr_mat)


      WRITE(rst_unit) num_ordersx,num_ordersy,num_ordersz
      WRITE(rst_unit) order_mapx,order_mapy,order_mapz


      bxdum = 0.0d0
      bxdum(1:nx,1:nx,:) = bx(:,:,:)
      WRITE(rst_unit) bxdum
      bxdum = 0.0d0
      bxdum(1:nx,1:nx,:) = by(:,:,:)
      WRITE(rst_unit) bxdum
      bxdum = 0.0d0
      bxdum(1:nx,1:nx,:) = bz(:,:,:)
      WRITE(rst_unit) bxdum
      bxdum = 0.0d0
      bxdum(1:nx,1:nx,:) = dmx(:,:,:)
      WRITE(rst_unit) bxdum
      bxdum = 0.0d0
      bxdum(1:nx,1:nx,:) = dmy(:,:,:)
      WRITE(rst_unit) bxdum
      bxdum = 0.0d0
      bxdum(1:nx,1:nx,:) = dmz(:,:,:)
      WRITE(rst_unit) bxdum

      WRITE(rst_unit) num_mrtr_matrices
      WRITE(rst_unit) mrtr_mat_map


      prolongdum = 0.0d0
      prolongdum(1:nmax,1:nmax,:) = prolong(:,:,:)
      WRITE(rst_unit) prolongdum
      prolongdum = 0.0d0
      prolongdum(1:nmax,1:nmax,:) = restrict(:,:,:)
      WRITE(rst_unit) prolongdum

      RETURN
      END SUBROUTINE write_restart_matrices

!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart_domain(rst_unit,da)
!
!......................................................................
!     date: 01/22/02
!......................................................................
!
      USE domain_definition
!
      IMPLICIT none
!
      TYPE (domain) :: da
      INTEGER       :: rst_unit
      DOUBLE PRECISION :: Qdum(nx+nref,nx+nref,nx+nref,neq)
      DOUBLE PRECISION :: cxdum(nx+nref,3)
      DOUBLE PRECISION :: cxduma(nx+nref,3)
      DOUBLE PRECISION :: xgdum(3,nx+nref,nx+nref,nx+nref)
      DOUBLE PRECISION :: jacobdum(nx+nref,nx+nref,nx+nref)
      DOUBLE PRECISION :: gmetdum(3,3,nx+nref,nx+nref,nx+nref)
      

      WRITE(rst_unit) da%type            ! (8 = hex8,... number of nodes in hex)
      WRITE(rst_unit) da%node
      WRITE(rst_unit) da%orientation 
      WRITE(rst_unit) da%mortar          ! mortar # and side for each face
      WRITE(rst_unit) da%which_surface   ! surface a face is attached to

      WRITE(rst_unit) da%nc,da%ncg    ! grid size

      Qdum = 0.0d0
      Qdum(1:nx,1:nx,1:nx,:) = da%Q(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%Qlgg(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%Qglg(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%Qggl(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%f(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%g(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%h(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
      Qdum(1:nx,1:nx,1:nx,:) = da%g_Q(:,:,:,:)
      WRITE(rst_unit) Qdum               ! dependent variables
!
      WRITE(rst_unit) da%material_id
!
      cxdum(:,:) = 0.0d0
      cxdum(1:nx,:) = da%cx(:,:)
      cxduma(:,:) = 0.0d0
      cxduma(1:nx,:) = da%cxg(:,:)
      WRITE(rst_unit) cxdum,cxduma
      xgdum(:,:,:,:) = 0.0d0
      xgdum(:,1:nx,1:nx,1:nx) = da%xg(:,:,:,:)
      WRITE(rst_unit) xgdum
      jacobdum(:,:,:) = 0.0d0
      jacobdum(1:nx,1:nx,1:nx) = da%jacob(:,:,:)
      WRITE(rst_unit) jacobdum
      gmetdum(:,:,:,:,:) = 0.0d0
      gmetdum(:,:,1:nx,1:nx,1:nx) = da%gmet(:,:,:,:,:) 
      WRITE(rst_unit) gmetdum
      gmetdum(:,:,1:nx,1:nx,1:nx) = da%gmetg(:,:,:,:,:) 
      WRITE(rst_unit) gmetdum
!
      WRITE(rst_unit) da%corner
!
      WRITE(rst_unit) da%ibtype
      WRITE(rst_unit) da%bcond
      WRITE(rst_unit) da%domain_type

      RETURN
      END SUBROUTINE write_restart_domain

!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart_mortar(rst_unit,mrtra)
!
!......................................................................
!     date: 01/22/02
!......................................................................
!
      USE mortar_definition
!
      IMPLICIT none
!
      TYPE (mortar) :: mrtra
      INTEGER       :: rst_unit
      DOUBLE PRECISION :: n_hatdum(3,nx+nref,nx+nref)
      DOUBLE PRECISION :: normdum(nx+nref,nx+nref)
      DOUBLE PRECISION :: Qdum(nx+nref,nx+nref,neq,2)

      WRITE(rst_unit) mrtra%lenmortar ! dimension of the mortar
      WRITE(rst_unit) mrtra%len     ! dimension of element face (2 directions x 2 sides)
      WRITE(rst_unit) mrtra%id        ! id of the contributing elements
      WRITE(rst_unit) mrtra%iface     ! face # of contributing elements
      WRITE(rst_unit) mrtra%orient       ! orientation of #2 element face (orient*pi/2)
      WRITE(rst_unit) mrtra%sign         ! orientation of normals (+1/-1)
      WRITE(rst_unit) mrtra%nsign         ! orientation of normals (+1/-1)
!
      n_hatdum(:,:,:) = 0.0d0
      n_hatdum(:,1:nx,1:nx) = mrtra%n_hat
      WRITE(rst_unit) n_hatdum
      normdum(:,:) = 0.0d0
      normdum(1:nx,1:nx) = mrtra%norm(:,:)
      WRITE(rst_unit) normdum         ! geometry arrays along mortar
      Qdum(:,:,:,:) = 0.0d0
      Qdum(1:nx,1:nx,:,:) = mrtra%Q(:,:,:,:)
      WRITE(rst_unit) Qdum            ! solutions along mortar
      Qdum(1:nx,1:nx,:,:) = mrtra%fv(:,:,:,:)
      WRITE(rst_unit) Qdum            ! solutions along mortar
      n_hatdum(:,1:nx,1:nx) = mrtra%xg(:,:,:)
      WRITE(rst_unit) n_hatdum           ! positions along mortar
      WRITE(rst_unit) mrtra%conforming   ! mortar is conforming or not (2 directions x 2 sides)

      RETURN
      END SUBROUTINE write_restart_mortar
!
!///////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart_drop(drop,rst_unit)
!
      USE particle_definition
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
      TYPE(particle) :: drop(npart)
!
      INTEGER :: rst_unit
!

      DO ip=1,npart
        IF (drop(ip)%onoff ==1) THEN
          write(rst_unit) drop(ip)%Xp
          write(rst_unit) drop(ip)%Vp
          write(rst_unit) drop(ip)%Tp
          write(rst_unit) drop(ip)%Mp
          write(rst_unit) drop(ip)%Xpnm
          write(rst_unit) drop(ip)%Vp1nm
          write(rst_unit) drop(ip)%Tpnm
          write(rst_unit) drop(ip)%Mpnm
          write(rst_unit) drop(ip)%ngrid
          write(rst_unit) drop(ip)%onoff
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE
