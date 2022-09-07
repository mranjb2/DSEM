!> \file FDerivLibb.f90
!! Spectral derivative library
!///////////////////////////////////////////////////////////////////////
!////////							////////
!////////	DerivLib.f90					////////
!////////							////////
!////////	contains:					////////
!////////							////////
!////////	   SUBROUTINE dx3dg(u,nc,ud,ncg,scale,dmat)	////////
!////////	   SUBROUTINE dy3dg(u,nc,ud,ncg,scale,dmat)	////////
!////////	   SUBROUTINE dz3dg(u,nc,ud,ncg,scale,dmat)	////////
!////////							////////
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes x derivatives on gauss grid from gauss-lobatto data.
!!
!! Note: u is assumed to be on a semi-staggered (nc x mcg x lcg) grid
!! ud is returned on an (ncg x mcg x lcg) grid  
!
      SUBROUTINE dx3dg(u,nc,ud,ncg,scale,dmat) 
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
      DOUBLE PRECISION      :: dmat(nx,nx),u(nx,ny,nz),ud(nx,ny,nz) 
      DOUBLE PRECISION      :: uu(nc(1),ncg(1)*ncg(1)),uud(ncg(1),ncg(1)*ncg(1))
!
      nca = nc(1)
      ncga = ncg(1)
      DO i=1,nca
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga
            uu(i,kk) = u(i,j,k)
          ENDDO
        ENDDO
      ENDDO

      CALL mxm(dmat(1:ncg(1),1:nc(1)),ncg(1),uu,nc(1),uud,ncg(2)*ncg(3))

      DO i=1,ncga
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga
            ud(i,j,k) = uud(i,kk)
          ENDDO
        ENDDO
      ENDDO
!
      RETURN 
      END SUBROUTINE dx3dg
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes y derivatives on gauss grid from gauss-lobatto data.
!!
!
   SUBROUTINE dy3dg(u,nc,ud,ncg,scale,dmat) 
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
      DOUBLE PRECISION      :: dmat(ny,ny),u(nx,ny,nz),ud(nx,ny,nz),dmatt(ny,ny) 
      DOUBLE PRECISION      :: uu(nc(1),ncg(1)*ncg(1)),uud(ncg(1),ncg(1)*ncg(1))
!
      CALL transpose(dmatt(1:nc(2),1:ncg(2)),nc(2),dmat(1:ncg(2),1:nc(2)),ncg(2))    
      DO k=1,ncg(3)                 
         CALL mxm(u(1:ncg(1),1:nc(2),k),ncg(1),dmatt(1:nc(2),1:ncg(2)),nc(2),ud(1:ncg(1),1:ncg(2),k),ncg(2))
      ENDDO

      RETURN 
   END SUBROUTINE dy3dg
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes z derivatives on gauss grid from gauss-lobatto data.
!!
!
   SUBROUTINE dz3dg(u,nc,ud,ncg,scale,dmat) 
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
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: dmat(nz,nz),u(nx,ny,nz),ud(nx,ny,nz),dmatt(nz,nz)
      DOUBLE PRECISION      :: uu(ncg(1)*ncg(2),nc(3)),uud(ncg(1)*ncg(2),ncg(3)+1)
!
      nca = nc(3)
      ncga = ncg(1)
      DO i=1,nca
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga
            uu(kk,i) = u(j,k,i)
          ENDDO
        ENDDO
      ENDDO

      CALL transpose(dmatt(1:nc(3),1:ncg(3)),nc(3),dmat(1:ncg(3),1:nc(3)),ncg(3))    
      CALL mxm(uu,ncg(1)*ncg(2),dmatt(1:nc(3),1:ncg(3)),nc(3),uud,ncg(3))

      DO i=1,nca
        DO j=1,ncga
          DO k=1,ncga
            kk = j + (k-1)*ncga
            ud(j,k,i) = uud(kk,i)
          ENDDO
        ENDDO
      ENDDO
!
      RETURN 
   END SUBROUTINE dz3dg
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes x derivatives on gauss grid from gauss-lobatto data.
!!
!! Note: u is assumed to be on a semi-staggered (nc x mcg x lcg) grid
!! ud is returned on an (ncg x mcg x lcg) grid.
!
      SUBROUTINE dx3dga(u,nc,ud,ncg,scale,dmat) 
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
      DOUBLE PRECISION      :: dmat(nx,nx),u(nx,ny,nz),ud(nx,ny,nz) 
!

      CALL mxm(dmat,nx,u(1,1,1),nx,ud(1,1,1),ny*nz)
 
!
      RETURN 
      END SUBROUTINE dx3dga
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes y derivatives on gauss grid from gauss-lobatto data.
!!
!
   SUBROUTINE dy3dga(u,nc,ud,ncg,scale,dmat) 
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
      DOUBLE PRECISION      :: dmat(ny,ny),u(nx,ny,nz),ud(nx,ny,nz),dmatt(ny,ny) 

      CALL transpose( dmatt(1:ny,1:ny) ,ny, dmat(1:ny,1:ny) ,ny)    
      DO k=1,nz                 
         CALL mxm(u(1,1,k),ny,dmatt,ny,ud(1,1,k),nx)
      ENDDO



      RETURN 
   END SUBROUTINE dy3dga
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes z derivatives on gauss grid from gauss-lobatto data.
!!
!
   SUBROUTINE dz3dga(u,nc,ud,ncg,scale,dmat) 
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
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: dmat(nz,nz),u(nx,ny,nz),ud(nx,ny,nz),dmatt(nz,nz)
!
      CALL transpose(dmatt(1:nz,1:nz),nz,dmat(1:nz,1:nz),nz)    
      CALL mxm(u(1,1,1),nx*ny,dmatt,nz,ud(1,1,1),nz)


!
      RETURN 
   END SUBROUTINE dz3dga
!-----------------------------------------------------------------------
      subroutine transpose(a,m,b,n)
      implicit double precision(a-h,o-z)
      double precision a(m,n),b(n,m)
!
      do i=1,m
      do j=1,n
         a(i,j) = b(j,i)
      enddo
      enddo
      return
      end
!-----------------------------------------------------------------------

