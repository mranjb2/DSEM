!> \file preprocessing.f90
!! Preprocessing performed for MPI communication
!///////////////////////////////////////////////////////////////////////
!
!> @brief Find edge mortars needing preprocessing for programming convenience.
      SUBROUTINE preproc_mpi(nmort) 
!
!......................................................................
!     date: 05/23/02
!
!     stac3m version
!
!     find edge mortars needing preprocessing 
!     for programming convenience
!......................................................................
!
!
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)

!

      nedge = 0
      DO j = 1,nmort
         jl     = nmortmpi(j  ,1)
         idm    = nmortmpi(j  ,2)
         idml   = ngridmpi(idm,1)
         nprocm = ngridmpi(idm,2)
         ids    = nmortmpi(j  ,3)
         IF (ids /= 0) THEN
           nprocs = ngridmpi(ids,2)
           idsl   = ngridmpi(ids,1)
         ENDIF
         IF (ids /= 0 .AND. nprocm /= nprocs) THEN ! this is and edge mortar
           nedge = nedge + 1
           nmortedge(j,1) = 1 
           nmortedge(j,2) = nprocm
           nmortedge(j,3) = nprocs
         ELSE
           nmortedge(j,1) = 0
           nmortedge(j,2) = nprocm
           nmortedge(j,3) = nprocs
         ENDIF
      END DO
!
! retrieve information needed for copying mortars that
! are between processors
!
      nmortd = 0
      DO j=1,nmort
         IF ( nmortedge(j,1) == 1) THEN
           nmortd             = nmortd + 1
           nmort              = nmort+1
           nmortedge(nmort,1) = 1 
           nmortedge(nmort,2) = nmortedge(j,3) 
           nmortedge(nmort,3) = nmortedge(j,2) 
           nmortmpi(nmort,2)  = nmortmpi(j,3)
           nmortmpi(nmort,3)  = nmortmpi(j,2)
           nmortdouble(j)     = nmort
           nmortdouble(nmort)= j
         ENDIF
      ENDDO
!
      END SUBROUTINE preproc_mpi
!
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Find edge mortars needing preprocessing for programming convenience
      SUBROUTINE preproc_mpi_driver(nmort) 
!
!......................................................................
!     date: 05/23/02
!
!     stac3m version
!
!     find edge mortars needing preprocessing 
!     for programming convenience
!......................................................................
!
!
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION(a-h,o-z)

!

      nedge = 0
      DO j = 1,nmort
         jl     = nmortmpi(j  ,1)
         idm    = nmortmpi(j  ,2)
         idml   = ngridmpi(idm,1)
         nprocm = ngridmpi(idm,2)
         ids    = nmortmpi(j  ,3)
         IF (ids /= 0) THEN
           nprocs = ngridmpi(ids,2)
           idsl   = ngridmpi(ids,1)
         ENDIF
         IF (ids /= 0 .AND. nprocm /= nprocs) THEN ! this is and edge mortar
           nedge = nedge + 1
           nmortedge(j,1) = 1 
           nmortedge(j,2) = nprocm
           nmortedge(j,3) = nprocs
         ELSE
           nmortedge(j,1) = 0
           nmortedge(j,2) = nprocm
           nmortedge(j,3) = nprocs
         ENDIF
      END DO
!
      END SUBROUTINE preproc_mpi_driver
!
