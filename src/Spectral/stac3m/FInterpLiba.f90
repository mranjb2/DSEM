!>\file FinterpLiba.f90
!! Spectral interpolation library
!
!> @brief
!> This subroutine finds interpolation for computing x derivatives 
!> on gauss grid from gauss-lobatto data.
!!
!
!
      SUBROUTINE interpx(ncol,b,u,ncg,ugl,nc) 
!
!     .................................................................
!     date: 11/5/98                                                     
!     routines called: none                                             
!     applicability: 3d stac3m codes                                                 
!
!     Compute x derivatives on gauss grid from gauss-lobatto data       
!     Note: u is assumed to be on a semi-staggered (nc x mcg x lcg) grid      
!           ud is returned on an (ncg x mcg x lcg) grid                       
!     .................................................................
!
      USE size
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(nx,nx),u(nx,ny,nz,neq),ugl(nx,ny,nz,neq) 
      DOUBLE PRECISION      :: uugl(nc(1),ncg(2)*ncg(3)*neq),uu(ncg(1),ncg(2)*ncg(3)*neq)
!
!     DO i = 1,ncg(1)
!        DO j = 1,ncg(2)
!          DO k = 1,ncg(3)
!              kk = j + (k-1)*ncg(2)
!              uugl(i,kk) = 0.0d0
!              uu(i,kk)  = 0.0d0
!              ugl(i,j,k) = 0.0d0
!          ENDDO
!        END DO
!     END DO
!
      nca = nc(1)
      ncga = ncg(1)
      DO nv=1,neq
      DO i=1,ncga
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga +(nv-1)*ncga*ncga
            uu(i,kk) = u(i,j,k,nv)
          ENDDO
        ENDDO
      ENDDO
      ENDDO

      CALL mxm(b(1:nc(1),1:ncg(1)),nc(1),uu,ncg(1),uugl,ncg(2)*ncg(3)*neq)

      DO nv=1,neq
      DO i=1,nca
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga +(nv-1)*ncga*ncga
            ugl(i,j,k,nv) = uugl(i,kk)
          ENDDO
        ENDDO
      ENDDO
      ENDDO
!
      RETURN 
      END SUBROUTINE interpx
!
!///////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine finds interpolation for computing y derivatives 
!> on gauss grid from gauss-lobatto data.  
!!
!
   SUBROUTINE interpy(ncol,b,u,ncg,ud,nc)  
!
!     .................................................................
!     date: 11/5/98                                                     
!     routines called: none                                             
!     applicability: 3d stac3m codes                                                 
!
!
!     compute y derivatives on gauss grid from gauss-lobatto data       
!     .................................................................
!
      USE size
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(ny,ny),u(nx,ny,nz,neq),ud(nx,ny,nz,neq),dmatt(ny,ny) 
!

      CALL transpose( dmatt(1:ncg(2),1:nc(2)) ,ncg(2), b(1:nc(2),1:ncg(2)) ,nc(2))    
      DO nv=1,neq
      DO k=1,ncg(3)                 
         CALL mxm(u(1:ncg(1),1:ncg(2),k,nv),ncg(1),dmatt(1:ncg(2),1:nc(2)),ncg(2),ud(1:ncg(1),1:nc(2),k,nv),nc(2))
      ENDDO

      ENDDO

      RETURN 
   END SUBROUTINE interpy
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine finds interpolation for computing z derivatives 
!> on gauss grid from gauss-lobatto data.  
!!
!
   SUBROUTINE interpz(ncol,b,u,ncg,ud,nc) 
!
!     .................................................................
!     date: 11/5/98                                                     
!     routines called: none                                             
!     applicability: 3d stac3m codes                                                 
!
!
!     compute Z derivatives on gauss grid from gauss-lobatto data       
!     .................................................................
!
      USE size
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(nz,nz),u(nx,ny,nz,neq),ud(nx,ny,nz,neq),dmatt(nz,nz)
      DOUBLE PRECISION      :: uu(ncg(1)*ncg(2)*neq,ncg(3)),uud(ncg(1)*ncg(2)*neq,nc(3))
!

      nca = nc(3)
      ncga = ncg(1)
      DO nv=1,neq
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga + (nv-1)*ncga*ncga
      DO i=1,nca
            uu(kk,i) = u(j,k,i,nv)
      ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL transpose(dmatt(1:ncg(3),1:nc(3)),ncg(3),b(1:nc(3),1:ncg(3)),nc(3))    
      CALL mxm(uu,ncg(1)*ncg(2)*neq,dmatt(1:ncg(3),1:nc(3)),ncg(3),uud,nc(3))

      DO nv=1,neq
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga + (nv-1)*ncga*ncga
      DO i=1,nca
            ud(j,k,i,nv) = uud(kk,i)
      ENDDO
          ENDDO
        ENDDO
      ENDDO

!
      RETURN 
   END SUBROUTINE interpz

!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine finds interpolation for computing x derivatives 
!> on gauss grid from gauss-lobatto data.  
!!
!
      SUBROUTINE interpxa(ncol,b,u,ncg,ugl,nc) 
!
!     .................................................................
!     date: 11/5/98                                                     
!     routines called: none                                             
!     applicability: 3d stac3m codes                                                 
!
!     Compute x derivatives on gauss grid from gauss-lobatto data       
!     Note: u is assumed to be on a semi-staggered (nc x mcg x lcg) grid      
!           ud is returned on an (ncg x mcg x lcg) grid                       
!     .................................................................
!
      USE size
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(nx,nx),u(nx,ny,nz,neq),ugl(nx,ny,nz,neq) 
      DO nv=1,neq
        CALL mxm(b,nx,u(1,1,1,nv),nx,ugl(1,1,1,nv),ny*nz)
      ENDDO
!
      RETURN 
      END SUBROUTINE interpxa
!
!///////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine finds interpolation for computing y derivatives 
!> on gauss grid from gauss-lobatto data.  
!!
!
   SUBROUTINE interpya(ncol,b,u,ncg,ud,nc)  
!
!     .................................................................
!     date: 11/5/98                                                     
!     routines called: none                                             
!     applicability: 3d stac3m codes                                                 
!
!
!     compute y derivatives on gauss grid from gauss-lobatto data       
!     .................................................................
!
      USE size
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(ny,ny),u(nx,ny,nz,neq),ud(nx,ny,nz,neq),dmatt(ny,ny) 
!

      CALL transpose( dmatt(1:ny,1:ny) ,ny, b(1:ny,1:ny) ,ny)    
      DO nv=1,neq
      DO k=1,nz                 
         CALL mxm(u(1,1,k,nv),ny,dmatt,ny,ud(1,1,k,nv),nx)
      ENDDO
      ENDDO

      RETURN 
   END SUBROUTINE interpya
!
!///////////////////////////////////////////////////////////////////////

!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine finds interpolation for computing z derivatives 
!> on gauss grid from gauss-lobatto data.  
!!
!
   SUBROUTINE interpza(ncol,b,u,ncg,ud,nc) 
!
!     .................................................................
!     date: 11/5/98                                                     
!     routines called: none                                             
!     applicability: 3d stac3m codes                                                 
!
!
!     compute Z derivatives on gauss grid from gauss-lobatto data       
!     .................................................................
!
      USE size
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(nz,nz),u(nx,ny,nz,neq),ud(nx,ny,nz,neq),dmatt(nz,nz)
!


      CALL transpose(dmatt(1:nz,1:nz),nz,b(1:nz,1:nz),nz)    
      DO nv=1,neq
       CALL mxm(u(1,1,1,nv),nx*ny,dmatt,nz,ud(1,1,1,nv),nz)
      ENDDO


!
      RETURN 
   END SUBROUTINE interpza

