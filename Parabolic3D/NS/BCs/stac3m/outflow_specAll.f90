!> @file
!! Boundary conditions
!///////////////////////////////////////////////////////////////////////
!> @brief
!> Compute ouflow bc using characteristic splitting, given all variables specified 
!> in the external field. This bc is equivalent to an interface mortar condition
!> where the neighbor element is the exterior.
!!
!
      SUBROUTINE outflow_specAll(mrtr,flux) 
!
!......................................................................
!     date: 11/18/98
!
!     compute ouflow bc using characteristic splitting,                 
!     given all variables specified in the external field
!     this bc is equivalent to an interface mortar condition
!     where the neighbor element is the exterior.
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property1 = 0.0d0, property2 = 0.0d0
!-----
!
      !p = (mrtr%Q(2,2,5,2) - 0.5d0*(mrtr%Q(2,2,2,2)**2+mrtr%Q(2,2,3,2)**2+mrtr%Q(2,2,4,2)**2)/ &
      !    mrtr%Q(2,2,1,2))*0.4d0
      !write(*,*) 'chick p',p
      CALL interface_flux(mrtr,flux,property1,property2,0.0d0)
!

      RETURN 
      END                                           
