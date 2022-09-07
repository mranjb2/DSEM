!> \file gatherbcast_dt.f90
!! Communication of timestep
!//////////////////////////////////////////////////////////////
!
!> @brief Gathers and finds the minimum of time step size from all cores.
      SUBROUTINE gathbcast_dt(dt)
!
       USE mpi_par
       USE mpi
!
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
!
       DOUBLE PRECISION :: dtmin

       CALL MPI_ALLREDUCE(dt,dtmin,1,MPI_DOUBLE_PRECISION, MPI_MIN, &
                          comm1d, ierr)
       dt = dtmin

       RETURN
      END SUBROUTINE
