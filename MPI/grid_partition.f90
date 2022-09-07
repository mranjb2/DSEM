!> \file grid_partition.f90
!! For partitioning meshes using METIS
!///////////////////////////////////////////////////////////////////////
!
!> @brief Partition the grid and send the various partitions to
!! different nodes
     SUBROUTINE grid_partitioning(d,mrtr,ngrid,nmort)
!
!......................................................................
!     date: 01/18/02
!
!     stac3m version
!
!     partition the grid and send the various partitions to
!     different nodes
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE Order_Matrices
      USE Material_Properties
      USE mpi_par
      USE input_data
      USE mpi
!
      IMPLICIT none

!
      TYPE (mortar) :: mrtr(nmp)
      TYPE (mortar) :: mrtra
      TYPE (domain) :: d(ngp)
      TYPE (domain) :: da

      INTEGER       :: ngrid,nmort,id,idm,k,j,mid,epsa,inode,i
      INTEGER       :: nsnode,nrnode,ip,ns,idl

      INTEGER, DIMENSION(10) :: idp
      INTEGER                :: itag
      DOUBLE PRECISION, DIMENSION(numprocs) :: L
      DOUBLE PRECISION                      :: L1,domainlength

!
      CALL MPI_BCAST(ngrid,1,MPI_INTEGER,0,comm1d,ierr)
      CALL MPI_BCAST(nmort,1,MPI_INTEGER,0,comm1d,ierr)
!
!   Set stuff for domain decomposition for MPI (equidistant
!   parts in x-direction)
!
      domainlength = char_length
      idp(:) = 0
      L1   =  domainlength/REAL(numprocs)
      epsa =  L1/1000.0d0
      write(*,*) L1,epsa
      DO i =1,numprocs
        L(i) = REAL(i)*L1
      ENDDO
!
      DO id = 1,ngrid
        CALL zero_domain(da)
        da = d(id)
        IF (myid == 0) THEN
          IF (da%corner(1,1) < L(1)-epsa) THEN ! this is the master node
              inode                = 0
              idp(inode+1)         = idp(inode+1) + 1
              ngridmpi(id,1)       = idp(inode+1)
              ngridmpi(id,2)       = inode
              d(idp(inode+1))      = da
              if(idp(inode+1)==6) write(*,*) 'hello'
          ENDIF
          DO ip = 1,numprocs-1
            IF ((da%corner(1,1)>= L(ip)-epsa)                         &
                .AND. (da%corner(1,1)<L(ip+1)-epsa)) THEN  ! this a slave node
              inode                = ip
              idp(inode+1)         = idp(inode+1) + 1
              ngridmpi(id,1)       = idp(inode+1)
              ngridmpi(id,2)       = inode
            ENDIF
         ENDDO
       ENDIF
!
       CALL MPI_BCAST(ngridmpi,ngp*8*2,MPI_INTEGER,0,comm1d,ierr)
       CALL MPI_BCAST(idp(:),10,MPI_INTEGER,0,comm1d,ierr)
!
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
!
      ngridl = idp(myid+1)

      write(*,*) 'myid,ngridl',myid,ngridl
      !IF (myid ==0) THEN
      !  write(*,*) 'ngridl ',ngridmpi(:,1)
      !  write(*,*) 'nprocs ',ngridmpi(:,2)
      !ENDIF

!
!     ---------------
!     read in mortars
!     ---------------
!
      idp(:) = 0
      DO idm = 1,nmort
        j = idm
        IF (myid ==0) THEN
          CALL zero_mortar(mrtra)
          mrtra=mrtr(idm)
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
        ENDIF
!
        CALL MPI_BCAST(nmortmpi(:,:),nmp*8*3,MPI_INTEGER,0,comm1d,ierr)
        CALL MPI_BCAST(nslavempi(:,:),nmp*8*3,MPI_INTEGER,0,comm1d,ierr)
        CALL MPI_BCAST(idp(:),10,MPI_INTEGER,0,comm1d,ierr)
!
        IF (numprocs >1) THEN
          DO ip =1,numprocs-1
            id     = nmortmpi(idm,2)
            nsnode = 0
            nrnode = ip
            IF (myid == 0 .and.  ngridmpi(id,2)== nrnode) THEN
              CALL send_mrtr(mrtra,nrnode,j)
            ELSEIF (myid == nrnode .and. ngridmpi(id,2) == nrnode) THEN
              idl  = nmortmpi(idm,1)
              CALL recv_mrtr(mrtr(idl),nsnode,j)
            ENDIF
          ENDDO
        ENDIF
!
      END DO ! end mortar loop
!
      nmortl = idp(myid +1)

      RETURN
     END SUBROUTINE grid_partitioning
