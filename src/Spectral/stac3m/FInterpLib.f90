!> \file FInterpLib.f90
!! Interpolation library
!///////////////////////////////////////////////////////////////////////
!> @brief
!> Perform interpolation from Gauss to Lobatto grid.
!!
!
      SUBROUTINE interp(ncol,b,u,ncg,ugl,nc) 
!
!     date: 2/13/95                                                     
!     routines called: none                                             
!     includes: none                                                    
!     applicability:all                                                 
!
!
!     PERFORM INTERPOLATION FROM GAUSS TO LOBATTO GRID                  
!
      implicit double precision (a-h,o-z) 
      dimension b(ncol,*),u(*),ugl(*) 
!
      call mxv(ncol,b,nc,u,ncg,ugl) 
!
      return 
      END                                           
