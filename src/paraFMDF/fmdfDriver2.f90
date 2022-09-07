!> \file fmdfDriver2.f90
!! Main driver file for FMDF logic, required for FMDF

!> @brief Sets up the FMDF solution in preprocessing
!!
!! This must be called before any FMDF routines are called. It counts the MC particles and places initial condition values on them.
subroutine fmdfPreprocessing(d,p,eP,ngridl,nprt,myid)
    use domain_definition
    use size 
    use physics
    use particle_definition
    use turbulence
    implicit none
    
    ! Domain Types
    type(domain), dimension(ngp)                    :: d            !< Carrier Domain list
    ! Particle Types
    type(particle), dimension(npart)                :: p            !< Particle
    type(ePointers), dimension(ngridl,nprt)         :: eP           !< Pointer to a local particle
    
    double precision                                :: x,y,z
    double precision, dimension(nx)                 :: hx, hy, hz
    integer, intent(in)                             :: ngridl       !< Local elements
    integer, intent(inout)                          :: nprt			!< number of local particles on the core
    integer, intent(in)                             :: myid			!< The processor
    integer                                         :: id, i, j, k, np  !< counters

    
!--- Begin logic

    ! Calculate sgs viscosity
    call simpleSmag(d,ngridl)
    
    ! Calculate subgrid constants
    call calcOmegaConstants
    
    ! Count the particles in each element and assign pointers to them
    call countParticles(d,p,eP,ngridl,nprt)
    
    ! Interpolate carrier phase quantities onto the particles
    !do np=1,nprt
    !    id = p(np)%ngrid
    !    !call Lagr_Pinterp(d(id),p(np),myid,np)
    !end do
    
    
    ! Calculate the particle weights and set the scalar IC's
    do id=1,ngridl
        !write(*,*) 'begin weight calc'
        call calcWeight(d(id),eP,id)
        !write(*,*) 'end weight calc'
        
        !write(*,*) 'begin scalar init'
        call initScalars(d(id),eP,id)
        !write(*,*) 'end scalar init'
    end do
    
end subroutine fmdfPreprocessing


!> @brief Main FMDF driver routine which calls all FMDF subroutines and processes
subroutine fmdfDriver(d,mrtr,p,eP,ngridl,nprt,tk,dt,myid)
    use domain_definition
    use mortar_definition
    use size
    use physics
    use particle_definition
    use constants
    use mpi
    use turbulence
    use random

    implicit none
    
    type(domain),       dimension(ngp)              :: d !< Carrier domain list
    type(mortar),       dimension(nmp)              :: mrtr !< Mortar list
    type(particle),     dimension(npart),   target  :: p !< particle list
    type(ePointers),    dimension(ngridl,nprt)      :: eP !< particle pointers that are subset of element
    
    double precision, intent(in)                    :: dt !< timestep size
    integer, intent(in)                             :: ngridl !< number of local elements to the core
    integer, intent(inout)                          :: nprt !< number of local particles to the core
    integer, intent(in)                             :: myid, tk
    integer                                         :: id, i, j, k, l,s,np
    
    double precision                                :: w, ex, ey , ez, scalar, short, b(3,2)
    double precision, dimension(3)                  :: Xmap
    integer                                         :: il, ir, jl, jr, kl, kr
    integer, dimension(3)                           :: ncg
    double precision, dimension(nx+2,ny+2,nz+2)     :: sumW,eRho,sumR
    double precision, dimension(nx+2,ny+2,nz+2,20)  :: sumT,asum
    double precision, dimension(ngridl,6,nx+2,ny+2,20)  :: mortSumT
    double precision, dimension(ngridl,6,nx+2,ny+2)     :: mortSumW
    double precision                                :: mysum
    integer, dimension(nx+2,ny+2,nz+2)              :: count
    integer, dimension(ngridl+2)                    :: mycount
    double precision, dimension(3)                  :: Xp
    double precision, parameter                     :: pie = 3.14159265358979323846264338327d0
    integer                                         :: loc,loc2
    double precision                                :: delta_g
    
    if(tk==1) then
        do id=1,ngridl
            do k=1,d(id)%ncg(3)
                do j=1,d(id)%ncg(2)
                    do i=1,d(id)%ncg(1)
                        do s=1,20
                            d(id)%s(i,j,k,s) = 0.0d0
                        end do
                    end do
                end do
            end do
        end do
    end if
    
    !call countParticles(d,p,eP,ngridl,nprt)
    if(tk==1) then
        do id=1,ngridl
            call initScalars(d(id),eP,id)
        end do
    end if
    
    if(tk==2) then
        do id=1,ngridl
            call initScalars(d(id),eP,id)
        end do
    end if
    
    ! Calculate eddy viscosity for sgs terms
    call simpleSmag(d,ngridl)

    ! Loop over all the particles
    ! Creates particle pointer array
    mycount(:) = 0
    do np=1,nprt
       !write(*,*) 'np',np,'on',p(np)%onoff,'ngrid',p(np)%ngrid,'z',p(np)%Xp(3)
        if( p(np)%onoff == 1 ) then ! Check if particle is live
            id = p(np)%ngrid
            Xmap = p(np)%Xpmap
            if ((Xmap(1) .ge. 0.0d0) .and. (Xmap(1) .le. 1.0d0)) then
                if((Xmap(2) .ge. 0.0d0) .and. (Xmap(2) .le. 1.0d0)) then
                    if((Xmap(3) .ge. 0.0d0) .and. (Xmap(3) .le. 1.0d0)) then
                        ! create a pointer to the particle
                        mycount(id) = mycount(id) + 1
                        eP(id,mycount(id))%p => p(np) 
                    end if
                end if
            end if
        end if
    end do
    
    ! Tell each ensemble how many particles it has
    do id=1,ngridl
        d(id)%np = mycount(id)
    end do

    ! Loop over all elements
    do id=1,ngridl
        eRho(:,:,:)   = 0.0d0
        sumW(:,:,:)   = 0.0d0
        sumT(:,:,:,:) = 0.0d0
        sumR(:,:,:)   = 0.0d0
        count(:,:,:)  = 0
        
        !call calcFMDFSums(e(id),eP(id,:))
        ncg = d(id)%ncg
        
        ! Check for an empty ensemble
        if( d(id)%np == 0 ) write(*,*) 'ENCOUNTERED AN EMPTY ENSEMBLE ELEMENT! ID=',id
        
        ! Loop over all particles belonging to the element
        do np=1,d(id)%np
            
            ! Get local quantities to the particle
            Xmap = eP(id,np)%p%Xpmap
            w    = eP(id,np)%p%w
            !write(*,*) 'w',w

            ! Calculate which gauss points the particle belongs to
            ex = (dble(ncg(1)) * acos(1.0d0 - 2.0d0 * Xmap(1)))/pie - 0.50d0
            ey = (dble(ncg(2)) * acos(1.0d0 - 2.0d0 * Xmap(2)))/pie - 0.50d0
            ez = (dble(ncg(3)) * acos(1.0d0 - 2.0d0 * Xmap(3)))/pie - 0.50d0

            ! Get the gauss point to the left and right of ex,ey,ez
            ! This is +1 because we do not want to yield gauss point 0
            ! This is corrected when the ends are chopped off in storage
!             il      = floor(ex) + 2
!             ir      = il + 1
!             b(1,1)  = ep(id,np)%p%h(1,:)
!             b(1,2)  = 1.0d0 !- b(1,1)
!             jl      = floor(ey) + 2
!             jr      = jl + 1
!             b(2,1)  = 1.0d0 !- (ey - dble(floor(ey)))
!             b(2,2)  = 1.0d0 !- b(2,1)
!             kl      = floor(ez) + 2
!             kr      = kl + 1
!             b(3,1)  = 1.0d0 !- (ez - dble(floor(ez)))
!             b(3,2)  = 1.0d0 !- b(3,1)
            
            ! Storage of basis values
!             eP(id,np)%p%b(1,1) = b(1,1)
!             eP(id,np)%p%b(1,2) = b(1,2)
!             eP(id,np)%p%b(2,1) = b(2,1)
!             eP(id,np)%p%b(2,2) = b(2,2)
!             eP(id,np)%p%b(3,1) = b(3,1)
!             eP(id,np)%p%b(3,2) = b(3,2)
            
            ! Tell the particle which gauss points it belongs to
!             eP(id,np)%p%il = il
!             eP(id,np)%p%jl = jl
!             eP(id,np)%p%kl = kl
            !write(*,*) 'il',il,'ir',ir
            !write(*,*) 'jl',jl,'jr',jr
            !write(*,*) 'kl',kl,'kr',kr
            !write(*,*) 'il',il,'jl',jl,'kl',kl
            
            ! Construct denominator sums on all gauss point permutations of these values
!             sumW(il,jl,kl) = sumW(il,jl,kl) + w*b(1,1)*b(2,1)*b(3,1)
!             sumW(ir,jl,kl) = sumW(ir,jl,kl) + w*b(1,2)*b(2,1)*b(3,1)
!             sumW(il,jr,kl) = sumW(il,jr,kl) + w*b(1,1)*b(2,2)*b(3,1)
!             sumW(ir,jr,kl) = sumW(ir,jr,kl) + w*b(1,2)*b(2,2)*b(3,1)
!             sumW(il,jl,kr) = sumW(il,jl,kr) + w*b(1,1)*b(2,1)*b(3,2)
!             sumW(ir,jl,kr) = sumW(ir,jl,kr) + w*b(1,2)*b(2,1)*b(3,2)
!             sumW(il,jr,kr) = sumW(il,jr,kr) + w*b(1,1)*b(2,2)*b(3,2)
!             sumW(ir,jr,kr) = sumW(ir,jr,kr) + w*b(1,2)*b(2,2)*b(3,2)
            do i=1,ncg(1)
               do j=1,ncg(2)
                  do k=1,ncg(3)
                     sumW(i,j,k) = sumW(i,j,k) + w*ep(id,np)%p%h(1,i)*ep(id,np)%p%h(2,j)*ep(id,np)%p%h(3,k)
                  end do
               end do
            end do 
            !write(*,*) sumW

            ! Count that the particle contributed to the gauss point
!             count(il,jl,kl) = count(il,jl,kl) + 1
!             count(ir,jl,kl) = count(ir,jl,kl) + 1
!             count(il,jr,kl) = count(il,jr,kl) + 1
!             count(ir,jr,kl) = count(ir,jr,kl) + 1
!             count(il,jl,kr) = count(il,jl,kr) + 1
!             count(ir,jl,kr) = count(ir,jl,kr) + 1
!             count(il,jr,kr) = count(il,jr,kr) + 1
!             count(ir,jr,kr) = count(ir,jr,kr) + 1

            ! Construct numerator sums for all scalars
            do s=1,20
                scalar  = eP(id,np)%p%scalar(s)
                do i=1,ncg(1)
                   do j=1,ncg(2)
                      do k=1,ncg(3)
                         short = w*scalar*ep(id,np)%p%h(1,i)*ep(id,np)%p%h(2,j)*ep(id,np)%p%h(3,k)
                         sumT(i,j,k,s) = sumT(i,j,k,s) + short
                      end do
                   end do
                end do
                
                ! sumT(il,jl,kl,s) = sumT(il,jl,kl,s) + short*b(1,1)*b(2,1)*b(3,1)
!                 sumT(ir,jl,kl,s) = sumT(ir,jl,kl,s) + short*b(1,2)*b(2,1)*b(3,1)
!                 sumT(il,jr,kl,s) = sumT(il,jr,kl,s) + short*b(1,1)*b(2,2)*b(3,1)
!                 sumT(ir,jr,kl,s) = sumT(ir,jr,kl,s) + short*b(1,2)*b(2,2)*b(3,1)
!                 sumT(il,jl,kr,s) = sumT(il,jl,kr,s) + short*b(1,1)*b(2,1)*b(3,2)
!                 sumT(ir,jl,kr,s) = sumT(ir,jl,kr,s) + short*b(1,2)*b(2,1)*b(3,2)
!                 sumT(il,jr,kr,s) = sumT(il,jr,kr,s) + short*b(1,1)*b(2,2)*b(3,2)
!                 sumT(ir,jr,kr,s) = sumT(ir,jr,kr,s) + short*b(1,2)*b(2,2)*b(3,2)
            end do
            
            ! Calculation of ensemble averaged source term on gauss points
            short          = eP(id,np)%p%react*w
            sumR(il,jl,kl) = sumR(il,jl,kl) + short*b(1,1)*b(2,1)*b(3,1)
            sumR(ir,jl,kl) = sumR(ir,jl,kl) + short*b(1,2)*b(2,1)*b(3,1)
            sumR(il,jr,kl) = sumR(il,jr,kl) + short*b(1,1)*b(2,2)*b(3,1)
            sumR(ir,jr,kl) = sumR(ir,jr,kl) + short*b(1,2)*b(2,2)*b(3,1)
            sumR(il,jl,kr) = sumR(il,jl,kr) + short*b(1,1)*b(2,1)*b(3,2)
            sumR(ir,jl,kr) = sumR(ir,jl,kr) + short*b(1,2)*b(2,1)*b(3,2)
            sumR(il,jr,kr) = sumR(il,jr,kr) + short*b(1,1)*b(2,2)*b(3,2)
            sumR(ir,jr,kr) = sumR(ir,jr,kr) + short*b(1,2)*b(2,2)*b(3,2)
            
        end do ! Ends loop over particles
!        write(*,*) 'hi at 162'
        
        !loop over counts to determine if we need to add or remove particles
        do k=2,ncg(3)+1
            do j=2,ncg(2)+1
                do i=2,ncg(1)+1
                    
                    ! temp particle adding
                    ! if (count(i,j,k) .lt. 1) then
!                         write(*,*) 'Element ',id,' had only ',count(i,j,k), 'particles'
!                         do l=1,(20-count(i,j,k))
!                             !find an empty particle
!                             loc2 = 0
!                             do np=1,npart
!                                 if (p(np)%onoff==0) then
!                                     loc2 = np
!                                     if (np>nprt) nprt = np
!                                     goto 124
!                                 end if
!                             end do
! 124                         continue
!                             if(loc2==0) write(*,*) '** NO EMPTY PARTICLES TO FILL!!!!!!'
!
!                             !initialize the particle
!                             p(loc2)%Xp(1) = d(id)%xg(1,i-1,j-1,k-1) + random_normal()*delta_g(d(id),1,1,1)/500.0d0
!                             p(loc2)%Xp(2) = d(id)%xg(2,i-1,j-1,k-1) + random_normal()*delta_g(d(id),1,1,1)/500.0d0
!                             p(loc2)%Xp(3) = 0.050d0
!                             p(loc2)%onoff = 1
!                             p(loc2)%ngrid = id
!
!                             call part_map(d(id),p(loc2))
!                             CALL spec_Pinterp(d(id),p(loc2),myid,loc2)
!
!                             ! set the weighting
!                             p(loc2)%w = p(loc2)%Rhofp / (dble(d(id)%np)*d(id)%jacob(1,1,1))
!
!                             ! set the scalar properties
!                             p(loc2)%scalar(1) = p(loc2)%T
!                             p(loc2)%scalar(2) = d(id)%s(i-1,j-1,k-1,2)
!                             p(loc2)%scalar(3) = d(id)%s(i-1,j-1,k-1,3)
!                             p(loc2)%scalar(4) = d(id)%s(i-1,j-1,k-1,4)
!                             p(loc2)%scalar(5) = d(id)%s(i-1,j-1,k-1,5)
!                             p(loc2)%eRho      = p(loc2)%P*1.40d0/p(loc2)%scalar(1)
!                             p(loc2)%react = 0.0d0
!
!                         end do
!                     end if
                    
                    !if we need to add particles
                    ! if (count(i,j,k) .lt. 10) then
!                         write(*,*) 'Element ',id,' had only ',count(i,j,k), 'particles'
!                         do l=1,(10-count(i,j,k))
!                             ! find the highest weight particle
!                             do np=2,d(id)%np
!                                 loc = 1
!                                 if(max(eP(id,np)%p%w,eP(id,np-1)%p%w) == eP(id,np)%p%w) then
!                                     loc = np
!                                 end if
!                             end do
!
!                             ! find an empty particle
!                             loc2=0
!                             do np=1,npart !changed to npart from nprt
!                                 if (p(np)%onoff==0) then
!                                     loc2 = np
!                                     if (np>nprt) nprt = np
!                                     goto 123
!                                 end if
!                             end do
! 123                         continue
!                             if(loc2==0) write(*,*) '** NO EMPTY PARTICLES TO FILL!!!!!!'
!
!                             ! Copy over the particle so they are a mirror of each other
!                             p(loc2) = eP(id,loc)%p
!
!                             ! Find the new position of the particle
!                             Xp(:)           = eP(id,loc)%p%Xp(:)
!                             !write(*,*)'xp before',Xp
!                             !Xp(1)           = Xp(1) + random_normal()*d(id)%jacob(1,1,1)**(1.0d0/3.0d0)
!                             !Xp(2)           = Xp(2) + random_normal()*d(id)%jacob(1,1,1)**(1.0d0/3.0d0)
!                             !Xp(3)           = Xp(3) + random_normal()*d(id)%jacob(1,1,1)**(1.0d0/3.0d0)
!                             Xp(1)           = Xp(1) + random_normal()*delta_g(d(id),1,1,1)/40.0d0
!                             Xp(2)           = Xp(2) + random_normal()*delta_g(d(id),1,1,1)/40.0d0
!                             Xp(3)           = Xp(3) + random_normal()*delta_g(d(id),1,1,1)/40.0d0
!                             !write(*,*)'xp after',Xp
!                             p(loc2)%Xp(:)   = Xp(:)
!
!                             ! split the weights
!                             p(loc2)%w       = p(loc2)%w/2.0d0
!                             eP(id,loc)%p%w  = p(loc2)%w
!
!                             ! turn on the new particle
!                             p(loc2)%onoff   = 1
!                         end do
!                     end if
                    
                    !if we need to remove particles
!                    if (count(i,j,k) .gt. 100) then
!                        write(*,*) 'Element ',id,' had too many particles, removing ',count(i,j,k)-100, 'particles'
!                        do l=1,(count(i,j,k)-100)
!                            ! find the lowest weight particle
!                            loc = 1
!                            do np=2,d(id)%np
!                                if (eP(id,np)%p%onoff == 1) then
!                                    if(min(eP(id,np)%p%w,eP(id,np-1)%p%w) == eP(id,np)%p%w) then
!                                        loc = np
!                                    end if
!                                end if
!                            end do
!                            loc2 = 2
!                            do np=3,d(id)%np
!                                if (eP(id,np)%p%onoff == 1) then
!                                    if(np .ne. loc) then
!                                        if(min(eP(id,np)%p%w,eP(id,np-1)%p%w) == eP(id,np)%p%w) then
!                                            loc2 = np
!                                        end if
!                                    end if
!                                end if
!                            end do
!
!
!                            ! Find the new position of the particle
!                            Xp(1)   = eP(id,loc )%p%Xp(1)*eP(id,loc )%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w) &
!                               &    + eP(id,loc2)%p%Xp(1)*eP(id,loc2)%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w)
!                            Xp(2)   = eP(id,loc )%p%Xp(2)*eP(id,loc )%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w) &
!                               &    + eP(id,loc2)%p%Xp(2)*eP(id,loc2)%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w)
!                            Xp(3)   = eP(id,loc )%p%Xp(3)*eP(id,loc )%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w) &
!                               &    + eP(id,loc2)%p%Xp(3)*eP(id,loc2)%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w)
!                            !write(*,*)'xp after',Xp
!                            eP(id,loc)%p%Xp(:) = Xp(:)
!
!                            ! add the scalars
!                            do s=1,5
!                                eP(id,loc)%p%scalar(s) = eP(id,loc )%p%scalar(s)*eP(id,loc )%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w) &
!                                                       + eP(id,loc2)%p%scalar(s)*eP(id,loc2)%p%w/(eP(id,loc)%p%w+eP(id,loc2)%p%w)
!                            end do
!
!                            ! add weights
!                            eP(id,loc)%p%w  = eP(id,loc)%p%w + eP(id,loc2)%p%w
!
!                            ! turn on the new particle, off the old particle
!                            eP(id,loc)%p%onoff  = 1
!                            eP(id,loc2)%p%onoff = 0
!                        end do
!                    end if
                    
                    
                end do
            end do
        end do!ends loop over counts
        
        ! Store the mortar values for passing
!        mortSumT(id,1,:,:,:)    =   sumT(    :,    1,    :, :)                  !y-dir bottom face
!        mortSumT(id,2,:,:,:)    =   sumT(    :, ny+2,    :, :)                  !y-dir top face
!        mortSumT(id,3,:,:,:)    =   sumT(    :,    :,    1, :)                  !z-dir in face
!        mortSumT(id,4,:,:,:)    =   sumT( nx+2,    :,    :, :)                  !x-dir right face
!        mortSumT(id,5,:,:,:)    =   sumT(    :,    :, nz+2, :)                  !z-dir out face
!        mortSumT(id,6,:,:,:)    =   sumT(    1,    :,    :, :)                  !x-dir left face
        
!        mortSumW(id,1,:,:)      =   sumW(    :,    1,    :)                     !y-dir bottom face
!        mortSumW(id,2,:,:)      =   sumW(    :, ny+2,    :)                     !y-dir top face
!        mortSumW(id,3,:,:)      =   sumW(    :,    :,    1)                     !z-dir in face
!        mortSumW(id,4,:,:)      =   sumW( nx+2,    :,    :)                     !x-dir right face
!        mortSumW(id,5,:,:)      =   sumW(    :,    :, nz+2)                     !z-dir out face
!        mortSumW(id,6,:,:)      =   sumW(    1,    :     :)                     !x-dir left face
        
        d(id)%sumT = sumT
        d(id)%sumW = sumW
        d(id)%sumR = sumR
        
    end do ! End first loop over all elements

    
!    il = 0
!    ir = 0
!    jl = 0
!    jr = 0
!    do id=1,nmp ! Second loop over all mortars
!        ! Here we take the partial sums from the element mortars and add them to interior points
!        
!        ! get the two elements involved
!        il  =   mrtr(id)%id(1)
!        ir  =   mrtr(id)%id(2)
!        
!        ! determine their faces involved
!        jl  =   mrtr(id)%iface(1)
!        jr  =   mrtr(id)%iface(2)
!        
!        
!        select case(jl)
!             case(1)
!                 d(ir)%sumT(:,ncg(2)+1,:,:)  = d(ir)%sumT(:,ncg(2)+1,:,:)    + d(il)%sumT(:,1,:,:)
!                 d(ir)%sumW(:,ncg(2)+1,:)    = d(ir)%sumW(:,ncg(2)+1,:)      + d(il)%sumW(:,1,:)
!                 d(ir)%sumR(:,ncg(2)+1,:)    = d(ir)%sumR(:,ncg(2)+1,:)      + d(il)%sumR(:,1,:)
!                 d(il)%sumT(:,       1,:,:)  = d(ir)%sumT(:,ncg(2)+1,:,:)
!                 d(il)%sumW(:,       1,:)    = d(ir)%sumW(:,ncg(2)+1,:)
!                 d(il)%sumR(:,       1,:)    = d(ir)%sumR(:,ncg(2)+1,:)
!             case(2)
!                 d(ir)%sumT(:,       2,:,:)  = d(ir)%sumT(:,2,:,:)           + d(il)%sumT(:,ncg(2)+2,:,:)
!                 d(ir)%sumW(:,       2,:)    = d(ir)%sumW(:,2,:)             + d(il)%sumW(:,ncg(2)+2,:)
!                 d(ir)%sumR(:,       2,:)    = d(ir)%sumR(:,2,:)             + d(il)%sumR(:,ncg(2)+2,:)
!                 d(il)%sumT(:,ncg(2)+2,:,:)  = d(ir)%sumT(:,2,:,:)
!                 d(il)%sumW(:,ncg(2)+2,:)    = d(ir)%sumW(:,2,:)
!                 d(il)%sumR(:,ncg(2)+2,:)    = d(ir)%sumR(:,2,:)
!             case(3)
!                 d(ir)%sumT(:,:,ncg(3)+1,:)  = d(ir)%sumT(:,:,ncg(3)+1,:)    + d(il)%sumT(:,:,1,:)
!                 d(ir)%sumW(:,:,ncg(3)+1)    = d(ir)%sumW(:,:,ncg(3)+1)      + d(il)%sumW(:,:,1)
!                 d(ir)%sumR(:,:,ncg(3)+1)    = d(ir)%sumR(:,:,ncg(3)+1)      + d(il)%sumR(:,:,1)
!                 d(il)%sumT(:,:,       1,:)  = d(ir)%sumT(:,:,ncg(3)+1,:)
!                 d(il)%sumW(:,:,       1)    = d(ir)%sumW(:,:,ncg(3)+1)
!                 d(il)%sumR(:,:,       1)    = d(ir)%sumR(:,:,ncg(3)+1)
!             case(4)
!                 d(ir)%sumT(       2,:,:,:)  = d(ir)%sumT(2,:,:,:)           + d(il)%sumT(ncg(1)+2,:,:,:)
!                 d(ir)%sumW(       2,:,:)    = d(ir)%sumW(2,:,:)             + d(il)%sumW(ncg(1)+2,:,:)
!                 d(ir)%sumR(       2,:,:)    = d(ir)%sumR(2,:,:)             + d(il)%sumR(ncg(1)+2,:,:)
!                 d(il)%sumT(ncg(1)+2,:,:,:)  = d(ir)%sumT(2,:,:,:)
!                 d(il)%sumW(ncg(1)+2,:,:)    = d(ir)%sumW(2,:,:)
!                 d(il)%sumR(ncg(1)+2,:,:)    = d(ir)%sumR(2,:,:)
!             case(5)
!                 d(ir)%sumT(:,:,       2,:)  = d(ir)%sumT(:,:,2,:)           + d(il)%sumT(:,:,ncg(3)+2,:)
!                 d(ir)%sumW(:,:,       2)    = d(ir)%sumW(:,:,2)             + d(il)%sumW(:,:,ncg(3)+2)
!                 d(ir)%sumR(:,:,       2)    = d(ir)%sumR(:,:,2)             + d(il)%sumR(:,:,ncg(3)+2)
!                 d(il)%sumT(:,:,ncg(3)+2,:)  = d(ir)%sumT(:,:,2,:)
!                 d(il)%sumW(:,:,ncg(3)+2)    = d(ir)%sumW(:,:,2)
!                 d(il)%sumR(:,:,ncg(3)+2)    = d(ir)%sumR(:,:,2)
!             case(6)
!                 d(ir)%sumT(ncg(1)+1,:,:,:)  = d(ir)%sumT(ncg(1)+1,:,:,:)    + d(il)%sumT(1,:,:,:)
!                 d(ir)%sumW(ncg(1)+1,:,:)    = d(ir)%sumW(ncg(1)+1,:,:)      + d(il)%sumW(1,:,:)
!                 d(ir)%sumR(ncg(1)+1,:,:)    = d(ir)%sumR(ncg(1)+1,:,:)      + d(il)%sumR(1,:,:)
!                 d(il)%sumT(       1,:,:,:)  = d(ir)%sumT(ncg(1)+1,:,:,:)
!                 d(il)%sumW(       1,:,:)    = d(ir)%sumW(ncg(1)+1,:,:)
!                 d(il)%sumR(       1,:,:)    = d(ir)%sumR(ncg(1)+1,:,:)
!         end select
!
!         select case(jr)
!             case(1)
!                 d(il)%sumT(:,ncg(2)+1,:,:)  = d(il)%sumT(:,ncg(2)+1,:,:)    + d(ir)%sumT(:,1,:,:)
!                 d(il)%sumW(:,ncg(2)+1,:)    = d(il)%sumW(:,ncg(2)+1,:)      + d(ir)%sumW(:,1,:)
!                 d(il)%sumR(:,ncg(2)+1,:)    = d(il)%sumR(:,ncg(2)+1,:)      + d(ir)%sumR(:,1,:)
!                 d(ir)%sumT(:,       1,:,:)  = d(il)%sumT(:,ncg(2)+1,:,:)
!                 d(ir)%sumW(:,       1,:)    = d(il)%sumW(:,ncg(2)+1,:)
!                 d(ir)%sumR(:,       1,:)    = d(il)%sumR(:,ncg(2)+1,:)
!             case(2)
!                 d(il)%sumT(:,       2,:,:)  = d(il)%sumT(:,2,:,:)           + d(ir)%sumT(:,ncg(2)+2,:,:)
!                 d(il)%sumW(:,       2,:)    = d(il)%sumW(:,2,:)             + d(ir)%sumW(:,ncg(2)+2,:)
!                 d(il)%sumR(:,       2,:)    = d(il)%sumR(:,2,:)             + d(ir)%sumR(:,ncg(2)+2,:)
!                 d(ir)%sumT(:,ncg(2)+2,:,:)  = d(il)%sumT(:,2,:,:)
!                 d(ir)%sumW(:,ncg(2)+2,:)    = d(il)%sumW(:,2,:)
!                 d(ir)%sumR(:,ncg(2)+2,:)    = d(il)%sumR(:,2,:)
!             case(3)
!                 d(il)%sumT(:,:,ncg(3)+1,:)  = d(il)%sumT(:,:,ncg(3)+1,:)    + d(ir)%sumT(:,:,1,:)
!                 d(il)%sumW(:,:,ncg(3)+1)    = d(il)%sumW(:,:,ncg(3)+1)      + d(ir)%sumW(:,:,1)
!                 d(il)%sumR(:,:,ncg(3)+1)    = d(il)%sumR(:,:,ncg(3)+1)      + d(ir)%sumR(:,:,1)
!                 d(ir)%sumT(:,:,       1,:)  = d(il)%sumT(:,:,ncg(3)+1,:)
!                 d(ir)%sumW(:,:,       1)    = d(il)%sumW(:,:,ncg(3)+1)
!                 d(ir)%sumR(:,:,       1)    = d(il)%sumR(:,:,ncg(3)+1)
!             case(4)
!                 d(il)%sumT(       2,:,:,:)  = d(il)%sumT(2,:,:,:)           + d(ir)%sumT(ncg(1)+2,:,:,:)
!                 d(il)%sumW(       2,:,:)    = d(il)%sumW(2,:,:)             + d(ir)%sumW(ncg(1)+2,:,:)
!                 d(il)%sumR(       2,:,:)    = d(il)%sumR(2,:,:)             + d(ir)%sumR(ncg(1)+2,:,:)
!                 d(ir)%sumT(ncg(1)+2,:,:,:)  = d(il)%sumT(2,:,:,:)
!                 d(ir)%sumW(ncg(1)+2,:,:)    = d(il)%sumW(2,:,:)
!                 d(ir)%sumR(ncg(1)+2,:,:)    = d(il)%sumR(2,:,:)
!             case(5)
!                 d(il)%sumT(:,:,       2,:)  = d(il)%sumT(:,:,2,:)           + d(ir)%sumT(:,:,ncg(3)+2,:)
!                 d(il)%sumW(:,:,       2)    = d(il)%sumW(:,:,2)             + d(ir)%sumW(:,:,ncg(3)+2)
!                 d(il)%sumR(:,:,       2)    = d(il)%sumR(:,:,2)             + d(ir)%sumR(:,:,ncg(3)+2)
!                 d(ir)%sumT(:,:,ncg(3)+2,:)  = d(il)%sumT(:,:,2,:)
!                 d(ir)%sumW(:,:,ncg(3)+2)    = d(il)%sumW(:,:,2)
!                 d(ir)%sumR(:,:,ncg(3)+2)    = d(il)%sumR(:,:,2)
!             case(6)
!                 d(il)%sumT(ncg(1)+1,:,:,:)  = d(il)%sumT(ncg(1)+1,:,:,:)    + d(ir)%sumT(1,:,:,:)
!                 d(il)%sumW(ncg(1)+1,:,:)    = d(il)%sumW(ncg(1)+1,:,:)      + d(ir)%sumW(1,:,:)
!                 d(il)%sumR(ncg(1)+1,:,:)    = d(il)%sumR(ncg(1)+1,:,:)      + d(ir)%sumR(1,:,:)
!                 d(ir)%sumT(       1,:,:,:)  = d(il)%sumT(ncg(1)+1,:,:,:)
!                 d(ir)%sumW(       1,:,:)    = d(il)%sumW(ncg(1)+1,:,:)
!                 d(ir)%sumR(       1,:,:)    = d(il)%sumR(ncg(1)+1,:,:)
!         end select
!
!     end do ! end the mortar loop
        
    do id=1,ngridl
        
        ncg = d(id)%ncg
        
        ! Store these as element values
        do k=1,ncg(3)
            do j=1,ncg(2)
                do i=1,ncg(1)
                    
                    do s=1,20 !CHANGE THIS TO NUMBER OF SCALARS
                        d(id)%s(i,j,k,s) = d(id)%sumT(i,j,k,s) / d(id)%sumW(i,j,k)
                    end do
                    d(id)%pCount(i,j,k) = count(i,j,k)
                    d(id)%react(i,j,k) = d(id)%sumR(i,j,k) / d(id)%sumW(i,j,k)
                    !d(id)%s(i,j,k,3) = dble(d(id)%pCount(i,j,k))
                    
                    ! Ensure non-zero
                    if (d(id)%s(i,j,k,2) .lt. 0.0d0) then
                        write(*,*) '**Encountered Ensemble with negative fuel!'
                        d(id)%s(i,j,k,2) = 0.0d0
                    end if
                    
                end do
            end do
        end do ! ends storage to element loop
        
        ! construct </rho>_l
        do k=1,ncg(3)+2
            do j=1,ncg(2)+2
                do i=1,ncg(1)+2
                    eRho(i,j,k) = d(id)%sumW(i,j,k) * d(id)%jacob(1,1,1)
                    do s=1,20
                        asum(i,j,k,s) = d(id)%sumT(i,j,k,s) / d(id)%sumW(i,j,k)
                    end do
                end do
            end do
        end do
        
        
        !loop over particles in the element once again
        do np=1,d(id)%np
!             b(1,1) = eP(id,np)%p%b(1,1)
!             b(1,2) = eP(id,np)%p%b(1,2)
!             b(2,1) = eP(id,np)%p%b(2,1)
!             b(2,2) = eP(id,np)%p%b(2,2)
!             b(3,1) = eP(id,np)%p%b(3,1)
!             b(3,2) = eP(id,np)%p%b(3,2)
            
            ! get bounding gauss point of particle
!             il = eP(id,np)%p%il
!             jl = eP(id,np)%p%jl
!             kl = eP(id,np)%p%kl
            
            ! Calculate <\rho>_l on the particle
            !eP(id,np)%p%eRho =  (   eRho(il  ,jl  ,kl  )    *b(1,1) *b(2,1) *b(3,1)     & 
            !                 &  +   eRho(il+1,jl  ,kl  )    *b(1,2) *b(2,1) *b(3,1)     &
            !                 &  +   eRho(il+1,jl+1,kl  )    *b(1,2) *b(2,2) *b(3,1)     &
            !                 &  +   eRho(il+1,jl+1,kl+1)    *b(1,2) *b(2,2) *b(3,2)     & 
            !                 &  +   eRho(il  ,jl+1,kl  )    *b(1,1) *b(2,2) *b(3,1)     &
            !                 &  +   eRho(il  ,jl+1,kl+1)    *b(1,1) *b(2,2) *b(3,2)     &
            !                 &  +   eRho(il  ,jl  ,kl+1)    *b(1,1) *b(2,1) *b(3,2)     &
            !                 &  +   eRho(il+1,jl  ,kl+1)    *b(1,2) *b(2,1) *b(3,2) )/  &
            !                 &  (   b(1,1) *b(2,1) *b(3,1)  +   b(1,2) *b(2,1) *b(3,1)  &
            !                 &  +   b(1,2) *b(2,2) *b(3,1)  +   b(1,2) *b(2,2) *b(3,2)  &
            !                 &  +   b(1,1) *b(2,2) *b(3,1)  +   b(1,1) *b(2,2) *b(3,2)  &
            !                 &  +   b(1,1) *b(2,1) *b(3,2)  +   b(1,2) *b(2,1) *b(3,2)  )
            
            eP(id,np)%p%eRho = eP(id,np)%p%P*1.40d0/eP(id,np)%p%scalar(1)
            
            ! calculate omega_m on the particle
            call calcOmega(eP(id,np)%p,d(id))
            
            ! Calculate scalars on particle
            do s=1,20
                ! mysum = (   asum(il  ,jl  ,kl  ,s) *b(1,1) *b(2,1) *b(3,1)      &
!                       & +   asum(il+1,jl  ,kl  ,s) *b(1,2) *b(2,1) *b(3,1)      &
!                       & +   asum(il+1,jl+1,kl  ,s) *b(1,2) *b(2,2) *b(3,1)      &
!                       & +   asum(il+1,jl+1,kl+1,s) *b(1,2) *b(2,2) *b(3,2)      &
!                       & +   asum(il  ,jl+1,kl  ,s) *b(1,1) *b(2,2) *b(3,1)      &
!                       & +   asum(il  ,jl+1,kl+1,s) *b(1,1) *b(2,2) *b(3,2)      &
!                       & +   asum(il  ,jl  ,kl+1,s) *b(1,1) *b(2,1) *b(3,2)      &
!                       & +   asum(il+1,jl  ,kl+1,s) *b(1,2) *b(2,1) *b(3,2)  )/  &
!                       & (   b(1,1) *b(2,1) *b(3,1)  +   b(1,2) *b(2,1) *b(3,1)  &
!                       & +   b(1,2) *b(2,2) *b(3,1)  +   b(1,2) *b(2,2) *b(3,2)  &
!                       & +   b(1,1) *b(2,2) *b(3,1)  +   b(1,1) *b(2,2) *b(3,2)  &
!                       & +   b(1,1) *b(2,1) *b(3,2)  +   b(1,2) *b(2,1) *b(3,2)  )
               mysum = 0.0d0
               do i=1,ncg(1)
                  do j=1,ncg(2)
                     do k=1,ncg(3)
                        short = eP(id,np)%p%scalar(s)
                        mysum = mysum + short*eP(id,np)%p%h(1,i)*eP(id,np)%p%h(2,j)*eP(id,np)%p%h(3,k)
                     end do
                  end do
               end do
                
                ! Numerical integration
                eP(id,np)%p%scalar(s) = ( - eP(id,np)%p%omega_m * ( eP(id,np)%p%scalar(s) - mysum )*dt + &
                                      &     eP(id,np)%p%scalar(s) )
                ! Analytical solution
                !eP(id,np)%p%scalar(s) = ( eP(id,np)%p%scalar(s) - mysum )*exp(-eP(id,np)%p%omega_m*dt) + mysum 
            end do
            
            !adds the compressibility term to FMDF temperature
            eP(id,np)%p%scalar(1) = eP(id,np)%p%scalar(1) + ((1.40d0-1.0d0)/eP(id,np)%p%eRho)*eP(id,np)%p%dp*dt
            
            if (eP(id,np)%p%scalar(2) .lt. 0.0d0) then
                write(*,*) '**Encountered Negative Fuel MF on particle ',np
                eP(id,np)%p%scalar(2) = 0.0d0
            end if
			! Add the reaction source term to the particle
            if(source_type=='fmdf') then
                !call singleStepSourceCp(eP(id,np)%p, dt)
                !call mechanismInterface(eP(id,np)%p, dt) !uncomment and fix this one
                !write(*,*)'* FMDF Source term applied to particles'
            end if
			
			! Add the compressibility term
			!eP(id,np)%p%scalar(1) = eP(id,np)%p%scalar(1) + ((1.40d0-1.0d0)/eP(id,np)%p%eRho)*eP(id,np)%p%dp*dt
			
        end do
        
    end do
    
end subroutine fmdfDriver