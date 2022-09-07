!#define 

!>\file staticSmag.f90
!! This file must be compiled when using any FMDF routines as it calculates turbulent
!! viscosity,\f$ \nu_t \f$, for use in the calculation of the SGS mixing frequency,
!! \f$ \Omega_m \f$, used in the scalar transport equation. 

!> Module for globally storing turbulence quantities
module turbulence
    use domain_definition
    implicit none
    double precision, parameter                 :: C_s = 0.17   !< Smagorinsky coefficient \f$C_s\f$ taken from Pope P.588
    double precision, dimension(ngp,nx,ny,nz)   :: nu_t         !< \f$\nu_t\f$, turbulent viscosity array
    double precision, dimension(ngp,nx,ny,nz)   :: nu_x         !< x-derivative of viscosity
    double precision, dimension(ngp,nx,ny,nz)   :: nu_y         !< y-derivative of viscosity
    double precision, dimension(ngp,nx,ny,nz)   :: nu_z         !< z-derivative of viscosity
    double precision, dimension(ngp,nx,ny,nz)   :: gamma_t      !< \f$\gamma_t\f$
end module turbulence

!> \brief
!! Subroutine to calculate the the modelled eddy viscosity for the simple Smagorinsky model
!!
!! This function makes use of filteredStrain() to calculate the sub-grid scale eddy 
!! viscosity. The turbulent viscosity is modelled as:
!! \f{eqnarray*}{
!!      \nu_t &=& \mathcal{l}^2_S \bar S \\
!!            &=& \left ( C_s \Delta \right)^2 \bar S 
!! \f}
!! where \f$ \mathcal{l}_S \f$ is the Smagorinsky lengthscale, \f$ \Delta \f$ is the filter width
!! and \f$ C_s \f$ is the Smagorinksy coefficient (defined in \ref turbulence ).
!! \callergraph \callgraph
subroutine simpleSmag(d,ngridl)
    use domain_definition
    use physics
    USE input_data
    implicit none
    
    type(domain),dimension(ngridl)              :: d            !< Domain
    integer,intent(in)                          :: ngridl       !< Number of grid points
    integer                                     :: id, eqs, x, y, z, nv, n, m, l
    double precision, dimension(nx,ny,nz,neq)   :: Qx, Qy, Qz
    double precision, dimension(nx,ny,nz)       :: nulgg, nuglg, nuggl
    double precision                            :: delta_g
    double precision                            :: filteredStrain
    double precision                            :: xloc,center,width,far,minimum,maximum
    
    !We loop over the elements of the grid
    do id=1,ngridl
        
        ! Now calculate the solution derivatives for each of the equations (to be used)
        do nv=1,neq
            call Gauss_Deriv(d(id)%gmetg,d(id)%Qlgg(:,:,:,nv),d(id)%Qglg(:,:,:,nv),     &
                         d(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv),   &
                         d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz                        )
!#ifdef verbose write(*,*)'qx',Qx
        end do
        
        ! finally calculate the sgs viscosity
        do z=1,d(id)%ncg(3)
            do y=1,d(id)%ncg(2)
                do x=1,d(id)%ncg(1)

                    d(id)%nu_t(x,y,z) = re * ( C_s * delta_g(d(id),x,y,z) ) *  &
                                   &         ( C_s * delta_g(d(id),x,y,z) ) *  &
                                   &         filteredStrain(d(id),Qx,Qy,Qz,x,y,z)
                    
                    ! modify smagorinsky viscosity:
                    if (shockSensor) d(id)%nu_t(x,y,z) = d(id)%nu_t(x,y,z) * (1.d0 - d(id)%ducros(x,y,z))
                    if (rhoSensor)   d(id)%nu_t(x,y,z) = d(id)%nu_t(x,y,z) * d(id)%rhoSensor(x,y,z)
                    
                    xloc = d(id)%xg(1,x,y,z)
                    
                    d(id)%nu_t(x,y,z) = min(max(d(id)%nu_t(x,y,z),0.d0),numax)
                    
                end do
            end do
        end do
        ! interpolate sgs viscosity to lobatto grid
        DO l = 1,d(id)%ncg(3)
           DO m = 1,d(id)%ncg(2)
               CALL interp(nx,d(id)%bx,d(id)%nu_t(:,m,l),d(id)%ncg(1),nulgg(:,m,l),d(id)%nc(1))
           END DO
        END DO
                                                                                                            
        DO l = 1,d(id)%ncg(3)
           DO n = 1,d(id)%ncg(1)
               CALL interp(ny,d(id)%by,d(id)%nu_t(n,:,l),d(id)%ncg(2),nuglg(n,:,l),d(id)%nc(2))
           END DO
        END DO
                                                                                                            
        DO m = 1,d(id)%ncg(2)
           DO n = 1,d(id)%ncg(1)
               CALL interp(nz,d(id)%bz,d(id)%nu_t(n,m,:),d(id)%ncg(3),nuggl(n,m,:),d(id)%nc(3))
           END DO
        END DO
        
        ! calculate the derivatives of the sgs viscosity on the gauss grid
        call Gauss_Deriv(d(id)%gmetg,nulgg(:,:,:),nuglg(:,:,:),                 &
                     nuggl(:,:,:),d(id)%nu_x(:,:,:),d(id)%nu_y(:,:,:),d(id)%nu_z(:,:,:), &
                     d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz                    )
    end do
    
end subroutine simpleSmag

!> \brief 
!! This function calculates the filtered rate of strain by using a static Smagorinsky model
!! as described on page 578 of Pope's "Turbulent Flows". 
!!
!! We begin by using the velocity gradients ( \f$ \partial \bar U_i / \partial x_j \f$ )
!! to calculate the filtered rate of strain tensor:
!! \f[
!!      \bar S_{ij} \equiv \frac{1}{2} \left ( \frac{\partial \bar U_i}{\partial x_j} +
!!                         \frac{\partial \bar U_j}{\partial x_i} \right )
!! \f]
!! Once we have the tensor terms, we then calculate the charactersic filtered rate of strain,
!! which is defined as:
!! \f[
!!      \bar S \equiv \left( 2 \bar S_{ij} \bar S_{ij} \right)^{1/2}
!! \f]
!!
!! @param S11 first tensor term
!! @param Snn consecutive tensor terms
!! @param S characterstic filtered rate of strain
!! \callergraph
!! \callgraph
double precision function filteredStrain(d,Qx,Qy,Qz,i,j,k) result(S)
    use domain_definition
    implicit none
    type(domain)            :: d                                    !< Domain
    double precision, dimension(nx,ny,nz,neq)       :: Qx           !< Spectral derivative of d%Q with respect to x
    double precision, dimension(nx,ny,nz,neq)       :: Qy           !< Spectral derivative of d%Q with respect to y
    double precision, dimension(nx,ny,nz,neq)       :: Qz           !< Spectral derivative of d%Q with respect to z
    double precision        :: rho,u,v,w                            !< Instantaneous solutions
    double precision        :: rho_x, rho_y, rho_z                  !< 
    double precision        :: u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z  !< Velocity Derivatives
    double precision        :: S11,S12,S13,S21,S22,S23,S31,S32,S33  !< Tensor Terms
    integer                 :: i,j,k
    
    ! Register the instantaneous velocities
    rho = d%Q(i,j,k,1) * d%jacob(1,1,1)
    u   = d%Q(i,j,k,2) / rho
    v   = d%Q(i,j,k,3) / rho
    w   = d%Q(i,j,k,4) / rho
!#ifdef verbose write(*,*)'in simpleSmag: rho=',rho
!#ifdef verbose write(*,*)'in simpleSmag: u=  ',u
!#ifdef verbose write(*,*)'in simpleSmag: v=  ',v
!#ifdef verbose write(*,*)'in simpleSmag: w=  ',w
    
    ! Register the instantaneous derivatives
    rho_x   = Qx(i,j,k,1) !* d%jacob(1,1,1)
    rho_y   = Qy(i,j,k,1) !* d%jacob(1,1,1)
    rho_z   = Qz(i,j,k,1) !* d%jacob(1,1,1)
    u_x     = ( Qx(i,j,k,2) - rho_x*u ) / rho
    u_y     = ( Qy(i,j,k,2) - rho_y*u ) / rho
    u_z     = ( Qz(i,j,k,2) - rho_z*u ) / rho
    v_x     = ( Qx(i,j,k,3) - rho_x*v ) / rho
    v_y     = ( Qy(i,j,k,3) - rho_y*v ) / rho
    v_z     = ( Qz(i,j,k,3) - rho_z*v ) / rho
    w_x     = ( Qx(i,j,k,4) - rho_x*w ) / rho
    w_y     = ( Qy(i,j,k,4) - rho_y*w ) / rho
    w_z     = ( Qz(i,j,k,4) - rho_z*w ) / rho
    
    ! Calculate the tensor components - from Pope "turbulent flows" PP.578 eq. 13.73
    S11 = u_x
    S12 = 0.50d0 * (u_y + v_x)
    S13 = 0.50d0 * (u_z + w_x)
    S21 = 0.50d0 * (v_x + u_y)
    S22 = v_y
    S23 = 0.50d0 * (v_z + w_y)
    S31 = 0.50d0 * (w_x + u_z)
    S32 = 0.50d0 * (w_y + v_z)
    S33 = w_z
    
    ! Calculate the characteristic filtered rate of strain - Pope 'turbulent flows' eq 13.74
    S   = sqrt( 2*S11*S11 + 2*S12*S12 + 2*S13*S13 + &
        &       2*S21*S21 + 2*S22*S22 + 2*S23*S23 + &
        &       2*S31*S31 + 2*S32*S32 + 2*S33*S33   )
    
end function filteredStrain

!> @brief Calculates \f$ \Delta_G \f$ for use in simpleSmag()
!!
!! This subroutine uses the Jacobian of the element to calculate the filter width,
!! \f$ \Delta_G \f$, used for the Smagorinsky-Lilly model.
!! \f[
!!      \Delta_G = \frac{ | J |^{1/3} }{ (P_x P_y P_z)^{1/3} }
!! \f]
!! Where \f$ P_i \f$ are the polynomial orders in each direction within the element.
!! \date 3-8-12
!! \callergraph \callgraph
double precision function delta_g(d,x,y,z) result(delta)
    use domain_definition
    implicit none
    type(domain)        :: d    !< Domain
    integer,intent(in)  :: x,y,z
    
    !delta_g is defined as the cube root of the jacobian
    delta = ( abs( 1.0d0/d%jacob(x,y,z) )**(1.0d0/3.0d0)) / &
          & ( ( d%ncg(1) ) * ( d%ncg(2) ) * ( d%ncg(3) ) )**(1.0d0/3.0d0)
!#ifdef verbose write(*,*)'delta',delta
end function delta_g
!
!=====================================================================================
!   Prolong nu_t values from ggg points to ggl points   
!  
      SUBROUTINE nut_prolong(d,ngrid)
!
!      USE size                                                         
      USE domain_definition 
!      USE mpi                                           
!
      IMPLICIT NONE
!
      TYPE (domain)   :: d(ngp)
!
      INTEGER         :: ngrid,id,m,n,l                
! 
      
      DO id= 1,ngrid
!
         DO l = 1,d(id)%ncg(3)
            DO m = 1,d(id)%ncg(2)
               CALL interp(nx,d(id)%bx,                           &
                              d(id)%nu_t    (:,m,l)   ,d(id)%ncg(1), &
                              d(id)%nu_t_lgg(:,m,l)   ,d(id)%nc(1)) 
!                              nu_t    (id,:,m,l)   ,d(id)%ncg(1), &
!                              nu_t_lgg(id,:,m,l)   ,d(id)%nc(1)) 
               DO n = 1,d(id)%nc(1)
                  IF (d(id)%nu_t_lgg(n,m,l) .lt. 0.d0) d(id)%nu_t_lgg(n,m,l)=0.d0
!                  IF (nu_t_lgg(id,n,m,l) .lt. 0.d0) nu_t_lgg(id,n,m,l)=0.d0
               END DO
            END DO
         END DO
!
         DO l = 1,d(id)%ncg(3)
            DO n = 1,d(id)%ncg(1)
               CALL interp(ny,d(id)%by,                            &
                              d(id)%nu_t    (n,:,l)   ,d(id)%ncg(2), &
                              d(id)%nu_t_glg(n,:,l)   ,d(id)%nc(2))
			   DO m = 1,d(id)%nc(2)
                  IF (d(id)%nu_t_glg(n,m,l) .lt. 0.d0) d(id)%nu_t_glg(n,m,l)=0.d0
               END DO
            END DO
         END DO
!
         DO m = 1,d(id)%ncg(2)
            DO n = 1,d(id)%ncg(1)
               CALL interp(nz,d(id)%bz,                            &
                              d(id)%nu_t    (n,m,:)   ,d(id)%ncg(3), &
                              d(id)%nu_t_ggl(n,m,:)   ,d(id)%nc(3))
               DO l = 1,d(id)%nc(3)
                  IF (d(id)%nu_t_ggl(n,m,l) .lt. 0.d0) d(id)%nu_t_ggl(n,m,l)=0.d0
               END DO
            END DO
         END DO         
!
      END DO 
!
      RETURN
!
      END SUBROUTINE nut_prolong
!
