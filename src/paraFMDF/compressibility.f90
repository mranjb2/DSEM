!> \file compressibility.f90
!! Contains compressibility terms for FMDF

!> Initializes the compressibility term for FMDF in preprocessing section of code
!! Calculates the material derivative of pressure and stores it in the domain
subroutine fmdfCompressInit(d,ngridl)
	use domain_definition
	use size
	use physics
	implicit none
	type(domain),dimension(ngp)	:: d !< element list
	double precision			:: u,v,w,rho,rhoe
	integer, intent(in)			:: ngridl !< number of local elements on processor
	integer						:: i,j,k,id
	
	
	do id=1,ngridl
		do k=1,d(id)%ncg(3)
			do j=1,d(id)%ncg(2)
				do i=1,d(id)%ncg(1)
					d(id)%DP(i,j,k) = 0.0d0
					rho				= d(id)%Q(i,j,k,1)*d(id)%jacob(i,j,k)
					u 				= d(id)%Q(i,j,k,2)/d(id)%Q(i,j,k,1)
					v 				= d(id)%Q(i,j,k,3)/d(id)%Q(i,j,k,1)
					w				= d(id)%Q(i,j,k,4)/d(id)%Q(i,j,k,1)
					rhoe			= d(id)%Q(i,j,k,5)*d(id)%jacob(i,j,k)
					d(id)%Po(i,j,k) = (rhoe-rho*0.50d0*(u*u+v*v+w*w))*(gamma-1.0d0)
					d(id)%Poo(i,j,k)= (rhoe-rho*0.50d0*(u*u+v*v+w*w))*(gamma-1.0d0)
				end do
			end do
		end do
	end do
end subroutine fmdfCompressInit

!> Updates the temporal compressibility term for FMDF
!! Caclulates the material derivative pressure term over all elements in the domain and stores values at the Gauss
!! quadrature points.
subroutine fmdfCompressUpdate(d,ngridl,dt)
	use domain_definition
	use size
	use physics
	implicit none
	type(domain), dimension(ngp)	:: d !< element list
	double precision, intent(in)	:: dt !< timestep size
	double precision				:: u,v,w,rho,p,rhoe
	integer, intent(in)				:: ngridl !< number of local elements on processor
	integer							:: i,j,k,id
	
	do id=1,ngridl
		do k=1,d(id)%ncg(3)
			do j=1,d(id)%ncg(2)
				do i=1,d(id)%ncg(1)
					rho				= d(id)%Q(i,j,k,1)*d(id)%jacob(i,j,k)
					u 				= d(id)%Q(i,j,k,2)/d(id)%Q(i,j,k,1)
					v 				= d(id)%Q(i,j,k,3)/d(id)%Q(i,j,k,1)
					w				= d(id)%Q(i,j,k,4)/d(id)%Q(i,j,k,1)
					rhoe			= d(id)%Q(i,j,k,5)*d(id)%jacob(i,j,k)
					p				= (rhoe-rho*0.50d0*(u*u+v*v+w*w))*(gamma-1.0d0)
					
					d(id)%DP(i,j,k) = (d(id)%Po(i,j,k) - d(id)%Poo(i,j,k))/dt
					d(id)%Poo(i,j,k)= d(id)%Po(i,j,k)
					d(id)%Po(i,j,k) = p
				end do
			end do
		end do
	end do
	
end subroutine fmdfCompressUpdate

!> @brief Updates the compressibility term for FMDF including spatial portion of the material derivative
!!
!! \f[ 
!! DP = \frac{\partial P}{\partial t} + u_i \frac{\partial P}{\partial x}
!! \f]
subroutine fmdfCompressSpace(d,ngridl,dt)
	use domain_definition
	use size
	use physics
	implicit none
	type(domain), dimension(ngp)	:: d !< element list
	double precision, intent(in)	:: dt !< timestep size
	double precision, dimension(nx,ny,nz)	:: rho, u, v, w, rhoe, p
	double precision, dimension(nx,ny,nz)	:: p_lgg, p_glg, p_ggl
	double precision, dimension(nx,ny,nz)	:: p_x, p_y, p_z
	integer, intent(in)				:: ngridl !< number of local elements on processor
	integer							:: i,j,k,id,l,m,n
	
	do id=1,ngridl
		rho(:,:,:) 	= d(id)%Q(:,:,:,1)*d(id)%jacob(:,:,:)
		u(:,:,:)   	= d(id)%Q(:,:,:,2)/d(id)%Q(:,:,:,1)
		v(:,:,:)	= d(id)%Q(:,:,:,3)/d(id)%Q(:,:,:,1)
		w(:,:,:)	= d(id)%Q(:,:,:,4)/d(id)%Q(:,:,:,1)
		rhoe(:,:,:)	= d(id)%Q(:,:,:,5)*d(id)%jacob(:,:,:)
		p(:,:,:)	= (rhoe(:,:,:) - rho(:,:,:)*0.50d0*(u(:,:,:)**2+v(:,:,:)**2+w(:,:,:)**2))*(gamma-1.0d0)
		
		!--- Calculate the spatial derivatives of pressure from this point
		! First interpolate to lobatto grid
		do l = 1,d(id)%ncg(3)
           do m = 1,d(id)%ncg(2)
               call interp(nx,d(id)%bx,p(:,m,l),d(id)%ncg(1),p_lgg(:,m,l),d(id)%nc(1))
           end do
        end do                                              
        do l = 1,d(id)%ncg(3)
           do n = 1,d(id)%ncg(1)
               call interp(ny,d(id)%by,p(n,:,l),d(id)%ncg(2),p_glg(n,:,l),d(id)%nc(2))
           end do
        end do                                                                                                  
        do m = 1,d(id)%ncg(2)
           do n = 1,d(id)%ncg(1)
               call interp(nz,d(id)%bz,p(n,m,:),d(id)%ncg(3),p_ggl(n,m,:),d(id)%nc(3))
           end do
        end do
		!now calculate the derivatives at gauss points
		call Gauss_Deriv(d(id)%gmetg,p_lgg(:,:,:),p_glg(:,:,:),     &
                     p_ggl(:,:,:),p_x(:,:,:),p_y(:,:,:),p_z(:,:,:), &
                     d(id)%ncg,d(id)%dmx,d(id)%dmy,d(id)%dmz		)
		
		!--- construct the material derivative of pressure term
		d(id)%DP(:,:,:) = ((p(:,:,:) - d(id)%Po(:,:,:))/dt) + u(:,:,:)*p_x(:,:,:) + v(:,:,:)*p_y(:,:,:) + w(:,:,:)*p_z(:,:,:)
		d(id)%Po(:,:,:) = p(:,:,:)
	end do
	
end subroutine fmdfCompressSpace
