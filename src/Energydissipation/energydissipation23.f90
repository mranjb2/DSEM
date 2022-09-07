
subroutine energydissipation(timea, da, ngridla)
   use mpi
   use mpi_par
   use size
   use constants
   use physics
   use domain_definition
   use input_data
   use File_Units
   implicit DOUBLE PRECISION (a-h,o-z)

   type(domain), dimension(ngridla)            :: da
   integer                                     :: ngridla
   double precision                            :: timea, sum, sum2, dissipatedenergy, Jacob
   double precision, dimension(nx,ny,nz,neq)   :: Qx, Qy, Qz
   double precision, dimension(neq)            :: Q, Qx2, Qy2, Qz2
   double precision, dimension(nx,ny,nz)       :: curlu1, curlu2, curlu3, vorticity

   sum = 0.d0
   sum2 =0.d0
   do id=1,ngridla
        
      do nv=1,neq
         call Gauss_Deriv(da(id)%gmetg,da(id)%Qlgg(:,:,:,nv),da(id)%Qglg(:,:,:,nv),     &
                          da(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv),   &
                          da(id)%ncg,da(id)%dmx,da(id)%dmy,da(id)%dmz                        )
      end do
        
      do l = 1,da(id)%ncg(1)
         do m = 1,da(id)%ncg(2)
            do n =1,da(id)%ncg(3)
   !------------------------------Jon-----------------------------------------------
   !            Q(1:neq)   = da(id)%Q(l,m,n,1:neq)*da(id)%jacob(l,m,n)
   !            Qx2(1:neq) = Qx(l,m,n,1:neq)
   !            Qy2(1:neq) = Qy(l,m,n,1:neq)
   !            Qz2(1:neq) = Qz(l,m,n,1:neq)
   !            rho  = Q(1)
   !            rhou = Q(2)
   !            rhov = Q(3)
   !            rhow = Q(4)
   !            rhoe = Q(5)
   !            u = rhou/rho
   !            v = rhov/rho
   !            w = rhow/rho
   !            ux = (Qx2(2) - u*Qx2(1))/rho
   !            uy = (Qy2(2) - u*Qy2(1))/rho
   !            uz = (Qz2(2) - u*Qz2(1))/rho
   !            vx = (Qx2(3) - v*Qx2(1))/rho
   !            vy = (Qy2(3) - v*Qy2(1))/rho
   !            vz = (Qz2(3) - v*Qz2(1))/rho
   !            wx = (Qx2(4) - w*Qx2(1))/rho
   !            wy = (Qy2(4) - w*Qy2(1))/rho
   !            wz = (Qz2(4) - w*Qz2(1))/rho
   !            cur = (wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2 
   !            vorticity(l,m,n) = cur/da(id)%jacob(l,m,n)
   !            vorticity(l,m,n) = cur*(lx/(nex*(nx-1)))*(ly/(ney*(ny-1)))*(lz/(nez*(nz-1)))
   !------------------------me-----------------------------------------------------------------
               udery            = ((Qy(l,m,n,2)-(da(id)%Q(l,m,n,2)/da(id)%Q(l,m,n,1))*Qy(l,m,n,1)))/da(id)%Q(l,m,n,1)!*da(id)%jacob(l,m,n))
               uderz            = ((Qz(l,m,n,2)-(da(id)%Q(l,m,n,2)/da(id)%Q(l,m,n,1))*Qz(l,m,n,1)))/da(id)%Q(l,m,n,1)!*da(id)%jacob(l,m,n))
               vderx            = ((Qx(l,m,n,3)-(da(id)%Q(l,m,n,3)/da(id)%Q(l,m,n,1))*Qx(l,m,n,1)))/da(id)%Q(l,m,n,1)!*da(id)%jacob(l,m,n))
               vderz            = ((Qz(l,m,n,3)-(da(id)%Q(l,m,n,3)/da(id)%Q(l,m,n,1))*Qz(l,m,n,1)))/da(id)%Q(l,m,n,1)!*da(id)%jacob(l,m,n))
               wderx            = ((Qx(l,m,n,4)-(da(id)%Q(l,m,n,4)/da(id)%Q(l,m,n,1))*Qx(l,m,n,1)))/da(id)%Q(l,m,n,1)!*da(id)%jacob(l,m,n))
               wdery            = ((Qy(l,m,n,4)-(da(id)%Q(l,m,n,4)/da(id)%Q(l,m,n,1))*Qy(l,m,n,1)))/da(id)%Q(l,m,n,1)!*da(id)%jacob(l,m,n))
               curlu1(l,m,n)    = (wdery-vderz)/da(id)%jacob(l,m,n)
               curlu2(l,m,n)    = (uderz-wderx)/da(id)%jacob(l,m,n)
               curlu3(l,m,n)    = (vderx-udery)/da(id)%jacob(l,m,n)
               vorticity(l,m,n) = curlu1(l,m,n)**2+curlu2(l,m,n)**2+curlu3(l,m,n)**2
               vorticity(l,m,n) = vorticity(l,m,n)/da(id)%jacob(l,m,n)
               ! vorticity(l,m,n) = vorticity(l,m,n)*(lx/(nex*(nx-1)))*(ly/(ney*(ny-1)))*(lz/(nez*(nz-1)))
   !------------------------------------------------------------------------------------------------------
              sum2             =sum2 +1.d0/da(id)%jacob(l,m,n)
    !         sum2             =sum2 +da(id)%jacob(l,m,n)
              sum              =sum + vorticity(l,m,n)
            end do
         end do
      end do

   end do

   call MPI_reduce(sum, dissipatedenergy , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm1d, ierr)
   call MPI_reduce(sum2, Jacob , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm1d, ierr)
   if (myid == 0) then
      dissipatedenergy = (dissipatedenergy/Jacob) * (1.d0/re)
!     open (unit =dissipated_unit, file = 'dissipated.plt')
!     call writedissipatedenergy(dissipated_unit, timea, dissipatedenergy)
      write(20,*) timea,dissipatedenergy
   end if

end subroutine energydissipation

subroutine writedissipatedenergy(unit, t, dissen)

   double precision :: t, dissen
   integer          :: unit

   write (unit,900)  t,dissen
   900 format (1x,f8.2 ,1x, f8.2)


end subroutine writedissipatedenergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!///////////////////////////////////////////////////////////////////////////////////////////////
!%This subroutines calculates the turbulent viscousity required to match DNS and LES.
!%Mathematical relations can be found in Dongro's Thesis page 171 to 173.
subroutine vcalculation(timea, da, ngridla)
   use mpi
   use mpi_par
   use size
   use domain_definition
   use File_Units

   implicit DOUBLE PRECISION (a-h,o-z)
   type(domain), dimension(ngridla)            :: da
   integer                                     :: ngridla,id,ngridaa,p,q,r
   double precision                            :: timea
   double precision, dimension(ngridla,9)      :: Ts
   double precision, dimension(9)              :: Tss,Tss0
   double precision, dimension(nx,ny,nz,3)     :: uf,ur
   double precision, dimension(nx,ny,nz)       :: trace
   double precision, dimension(nx,ny,nz)       :: count
   double precision, dimension(nx,ny,nz,4)     :: pav
   double precision, dimension(nx,ny,nz,3,3)   :: puuav,uuf
   double precision, dimension(nx,ny,nz,3)     :: pufav,purav,uff,urf
   double precision, dimension(nx,ny,nz,3,3)   :: puurav,uurf
   double precision, dimension(nx,ny,nz,3,3)   :: puruav,uruf
   double precision, dimension(nx,ny,nz,3,3)   :: pururav,ururf
   double precision, dimension(nx,ny,nz,3,3)   :: LL,MM,RR
   double precision, dimension(nx,ny,nz,9)     :: T,sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!//////////////////////////////////////////////////////////////////////////////////////////////
Ts  = 0.0d0
Tss = 0.0d0
do id=1,ngridla

  pav  = 0.0d0
!%%% Equation 6.45 of Dongru's Thesis %%%!

  do l=1,da(id)%ncg(3)
    do m=1,da(id)%ncg(2)
      do n=1,da(id)%ncg(1)
        do i=1,4
          pav(l,m,n,i) = da(id)%Q(l,m,n,i)*da(id)%jacob(l,m,n)
        end do

        x0 = da(id)%xg(1,l,m,n)
        y0 = da(id)%xg(2,l,m,n)
        z0 = da(id)%xg(3,l,m,n)

        do p=1,da(id)%ncg(3)
           do q=1,da(id)%ncg(2)
              do r=1,da(id)%ncg(1)
                x = da(id)%xg(1,p,q,r)
                y = da(id)%xg(2,p,q,r)
                z = da(id)%xg(3,p,q,r)
                distance = sqrt( (x-x0)**2+(y-y0)**2+(z-z0)**2 )
                if (0 <distance .and. distance<= 0.4d0 ) then
                    do i=1,4
                       pav(l,m,n,i) = pav(l,m,n,i)+da(id)%Q(p,q,r,i)*da(id)%jacob(p,q,r)
                    end do
                    count(l,m,n) = count(l,m,n)+1
                end if 
              end do
          end do
        end do

        do i=1,3
          uf(l,m,n,i) = pav(l,m,n,i+1)/pav(l,m,n,1)
        end do

      end do
    end do
  end do

! if(myid ==0) then
!   write(24,*) 'puav values '
!   write(24,*) pav(1,1,1,1), pav(1,1,1,2), pav(1,1,1,3), pav(1,1,1,4), count(1,1,1)
!   write(24,*) pav(2,2,2,1), pav(2,2,2,2), pav(2,2,2,3), pav(2,2,2,4), count(2,2,2)
! end if 

! if(myid ==0) then
!   write(24,*) 'uf values '
!   write(24,*) uf(1,1,1,1), uf(1,1,1,2), uf(1,1,1,3), uf(2,2,2,1), uf(2,2,2,2), uf(2,2,2,3)
! end if 

!%%%%% u =uf+ur !%%%%%
do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
        
        do i=1,3
          ur(l,m,n,i)=( da(id)%Q(l,m,n,i+1)/da(id)%Q(l,m,n,1)-uf(l,m,n,i) )
                      
        end do

     end do
   end do
end do

! if(myid ==0) then
!   write(24,*) 'u values'
!   write(24,*) da(id)%Q(1,1,1,2)/da(id)%Q(1,1,1,1), da(id)%Q(1,1,1,3)/da(id)%Q(1,1,1,1), da(id)%Q(1,1,1,4)/da(id)%Q(1,1,1,1)
! end if 

! if(myid ==0) then
!   write(24,*) 'ur values'
!   write(24,*) ur(1,1,1,1), ur(1,1,1,2), ur(1,1,1,3), ur(2,2,2,1), ur(2,2,2,2), ur(2,2,2,3)
! end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!% computation of Leonard Stress tensor, Equation 6-51 %!!!!
puuav = 0.0d0
pufav = 0.0d0
purav = 0.0d0


 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)

      rho = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)

      do i=1,3
        do j=1,3
          puuav(l,m,n,i,j) = rho*uf(l,m,n,i)*uf(l,m,n,j)
        end do
      end do

      do i=1,3
        pufav(l,m,n,i) = rho*uf(l,m,n,i)
        purav(l,m,n,i) = rho*ur(l,m,n,i)
      end do

      x0 = da(id)%xg(1,l,m,n)
      y0 = da(id)%xg(2,l,m,n)
      z0 = da(id)%xg(3,l,m,n)

      do p=1,da(id)%ncg(3)
        do q=1,da(id)%ncg(2)
          do r=1,da(id)%ncg(1)
            x = da(id)%xg(1,p,q,r)
            y = da(id)%xg(2,p,q,r)
            z = da(id)%xg(3,p,q,r)

            distance =sqrt( (x-x0)**2+(y-y0)**2+(z-z0)**2 )

            if (0 < distance .and. distance<= 0.4d0) then

              rhor = da(id)%Q(p,q,r,1)*da(id)%jacob(p,q,r)

              do i=1,3
                do j=1,3
                  puuav(l,m,n,i,j) = puuav(l,m,n,i,j)+rhor*uf(p,q,r,i)*uf(p,q,r,j)
                end do
              end do

              do i=1,3
                pufav(l,m,n,i) = pufav(l,m,n,i)+rhor*uf(p,q,r,i)
                purav(l,m,n,i) = purav(l,m,n,i)+rhor*ur(p,q,r,i)
              end do

            end if 

          end do
        end do
      end do

      uuf(l,m,n,:,:) = puuav(l,m,n,:,:)/pav(l,m,n,1)
      uff(l,m,n,:)   = pufav(l,m,n,:)/pav(l,m,n,1)
      urf(l,m,n,:)   = purav(l,m,n,:)/pav(l,m,n,1)

      end do
   end do
end do


! if(myid ==0 ) then
!   write(24,*) 'puuav             ','puvav            ','puwav           ','pvuav          '
!   write(24,*) puuav, puvav, puwav, pvuav, pvvav, pvwav, pwuav, pwvav, pwwav
!  end if 


! if(myid ==0) then
!   write(24,*) 'uuf values'
!   write(24,*) uuf(1,1,1,1,1), uuf(1,1,1,1,2), uuf(1,1,1,1,3), uuf(2,2,2,1,1), uuf(2,2,2,1,2), uuf(2,2,2,2,1)
! end if 

! if(myid ==0) then
!   write(24,*) 'uff values'
!   write(24,*) uff(1,1,1,1), uff(1,1,1,2), uff(1,1,1,3)
! end if 

! if(myid ==0) then
!   write(24,*) 'urf values'
!    do l=1,da(id)%ncg(3)
!     do m=1,da(id)%ncg(2)
!      do n=1,da(id)%ncg(1)
!        write(24,*) da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n,:),ur(l,m,n,:)
!      end do
!     end do
!   end do
!   write(24,*) da(id)%Q(1,1,1,1), da(id)%jacob(1,1,1),ur(1,1,1,1), purav(1,1,1,1), pav(1,1,1,1), urf(1,1,1,1), urf(1,1,1,2),&
!               urf(1,1,1,3), urf(3,3,3,1), urf(3,3,3,2), urf(3,3,3,3)
! end if 

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)

      do i=1,3
        do j=1,3
          LL(l,m,n,i,j)= ( pav(l,m,n,1)/count(l,m,n) )*(uuf(l,m,n,i,j)-uff(l,m,n,i)*uff(l,m,n,j))
        end do
      end do

     end do
   end do
end do

if(myid ==0) then
  write(24,*) 'LL values'
  write(24,*) LL(2,2,2,1,1), LL(2,2,2,1,2), LL(2,2,2,1,3), LL(2,2,2,2,1), LL(2,2,2,2,2), LL(2,2,2,2,3),&
              LL(2,2,2,3,1), LL(2,2,2,3,2), LL(2,2,2,3,3)
end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!/////computes the cross tensor!  equation 6-52////
puurav = 0.0d0
puruav = 0.0d0


 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)

      rho = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)

      do i=1,3
        do j=1,3
          puurav(l,m,n,i,j)= rho*uf(l,m,n,i)*ur(l,m,n,j)
          puruav(l,m,n,i,j)= rho*ur(l,m,n,i)*uf(l,m,n,j)
        end do
      end do

      x0 = da(id)%xg(1,l,m,n)
      y0 = da(id)%xg(2,l,m,n)
      z0 = da(id)%xg(3,l,m,n)

      do p=1,da(id)%ncg(3)
        do q=1,da(id)%ncg(2)
          do r=1,da(id)%ncg(1)
            x = da(id)%xg(1,p,q,r)
            y = da(id)%xg(2,p,q,r)
            z = da(id)%xg(3,p,q,r)
            distance = sqrt( (x-x0)**2+(y-y0)**2+(z-z0)**2 )

            if (0 < distance .and. distance<= 0.4d0 ) then

              rhor = da(id)%Q(p,q,r,1)*da(id)%jacob(p,q,r)

              puurav(l,m,n,i,j) = puurav(l,m,n,i,j)+rhor*uf(p,q,r,i)*ur(p,q,r,j)
              puruav(l,m,n,i,j) = puruav(l,m,n,i,j)+rhor*ur(p,q,r,i)*uf(p,q,r,j)

            end if 
          end do
        end do
      end do
      uurf(l,m,n,:,:)  = puurav(l,m,n,:,:)/pav(l,m,n,1)
      uruf(l,m,n,:,:)  = puruav(l,m,n,:,:)/pav(l,m,n,1)
     end do
   end do
end do


! if (myid==0) then
!   write(24,*) 'uur values'
!   write(24,*) '**********************'
!   write(24,*) uurf(2,2,2,1,1), uurf(2,2,2,1,2), uurf(2,2,2,1,3), uurf(2,2,2,2,1), uurf(2,2,2,2,2), uurf(2,2,2,2,3),&
!               uurf(2,2,2,3,1), uurf(2,2,2,3,2), uurf(2,2,2,3,3)
!   write(24,*) uurf(3,3,3,1,1), uurf(3,3,3,1,2), uurf(3,3,3,1,3), uurf(3,3,3,2,1), uurf(3,3,3,2,2),&
!               uurf(3,3,3,2,3), uurf(3,3,3,3,1), uurf(3,3,3,3,2), uurf(3,3,3,3,3)
! end if 

! if (myid==0) then
!   write(24,*) 'uru values'
!   write(24,*) uruf(2,2,2,1,1), uruf(2,2,2,1,2), uruf(2,2,2,1,3), uruf(2,2,2,2,1), uruf(2,2,2,2,2), uruf(2,2,2,2,3),&
!               uruf(2,2,2,3,1), uruf(2,2,2,3,2), uruf(2,2,2,3,3)
!   write(24,*) uruf(3,3,3,1,1), uruf(3,3,3,1,2), uruf(3,3,3,1,3), uruf(3,3,3,2,1), uruf(3,3,3,2,2),&
!               uruf(3,3,3,2,3), uruf(3,3,3,3,1), uruf(3,3,3,3,2), uruf(3,3,3,3,3)
! end if 


 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)

        do i=1,3
          do j=1,3
            MM(l,m,n,i,j)= ( pav(l,m,n,1)/count(l,m,n) )*(uurf(l,m,n,i,j)+uruf(l,m,n,i,j)-uff(l,m,n,i)*urf(l,m,n,j)&
                             &-urf(l,m,n,i)*uff(l,m,n,j))
          end do
        end do

     end do
   end do
end do

if (myid==0) then
  write(24,*) 'MM values'
  write(24,*) MM(2,2,2,1,1), MM(2,2,2,1,2), MM(2,2,2,1,3), MM(2,2,2,2,1), MM(2,2,2,2,2), MM(2,2,2,2,3),&
              & MM(2,2,2,3,1), MM(2,2,2,3,2), MM(2,2,2,3,3)
  ! write(24,*) MM(3,3,3,1,1), MM(3,3,3,1,2), MM(3,3,3,1,3), MM(3,3,3,2,1), MM(3,3,3,2,2), MM(3,3,3,2,3),&
  !             & MM(3,3,3,3,1), MM(3,3,3,3,2), MM(3,3,3,3,3)
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!/ computation of Reynolds strees Equation 6-53/!
pururav = 0.0d0

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      rho = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)

      do i=1,3
        do j=1,3
          pururav(l,m,n,i,j)= rho*ur(l,m,n,i)*ur(l,m,n,j)
        end do
      end do 

      x0 = da(id)%xg(1,l,m,n)
      y0 = da(id)%xg(2,l,m,n)
      z0 = da(id)%xg(3,l,m,n)

      do p=1,da(id)%ncg(3)
        do q=1,da(id)%ncg(2)
          do r=1,da(id)%ncg(1)
            x = da(id)%xg(1,p,q,r)
            y = da(id)%xg(2,p,q,r)
            z = da(id)%xg(3,p,q,r)
            distance = sqrt( (x-x0)**2+(y-y0)**2+(z-z0)**2 )

            if ( 0 < distance .and. distance<= 0.4d0 ) then
              rhor = da(id)%Q(p,q,r,1)*da(id)%jacob(p,q,r)
              do i=1,3
                do j=1,3
                  pururav(l,m,n,i,j) = pururav(l,m,n,i,j)+rhor*ur(p,q,r,i)*ur(p,q,r,j)
                end do
              end do 
  
            end if 

          end do
        end do
      end do

    ururf(l,m,n,:,:) = pururav(l,m,n,:,:)/pav(l,m,n,1)

    end do
  end do
end do


! if (myid==0) then
!   write(24,*) 'urur values'
!   write(24,*) ururf(2,2,2), urvrf(2,2,2), urwrf(2,2,2), vrurf(2,2,2), vrvrf(2,2,2), vrwrf(2,2,2), wrurf(2,2,2), wrvrf(2,2,2),&
!   &wrwrf(2,2,2)
!   write(24,*) ururf(3,3,3), urvrf(3,3,3), urwrf(3,3,3), vrurf(3,3,3), vrvrf(3,3,3), vrwrf(3,3,3), wrurf(3,3,3), wrvrf(3,3,3),&
!   &wrwrf(3,3,3)
! end if 

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)

       do i=1,3
          do j=1,3
            RR(l,m,n,i,j)= ( pav(l,m,n,1)/count(l,m,n) )*( (ururf(l,m,n,i,j)-( urf(l,m,n,i)*urf(l,m,n,j) ) ) )
          end do
        end do

     end do
   end do
 end do

 if (myid==0) then
  write(24,*) 'RR values'
  write(24,*) RR(2,2,2,1,1), RR(2,2,2,1,2), RR(2,2,2,1,3), RR(2,2,2,2,1), RR(2,2,2,2,2), RR(2,2,2,2,3),&
              & RR(2,2,2,3,1), RR(2,2,2,3,2), RR(2,2,2,3,3)
  ! write(24,*) RR(3,3,3,1,1), RR(3,3,3,1,2), RR(3,3,3,1,3), RR(3,3,3,2,1), RR(3,3,3,2,2), RR(3,3,3,2,3),&
  !             & RR(3,3,3,3,1), RR(3,3,3,3,2), RR(3,3,3,3,3)
 end if
 T = 0.0d0
 trace = 0.0d0
 do i=1,3
    do j=1,3
      T(:,:,:,j+(i-1)*3) = LL(:,:,:,i,j)+MM(:,:,:,i,j)+RR(:,:,:,i,j)
    end do
  end do

 trace = ( T(:,:,:,1)+T(:,:,:,5)+T(:,:,:,9) )
 T(:,:,:,1) = T(:,:,:,1)-1/3*trace
 T(:,:,:,5) = T(:,:,:,5)-1/3*trace
 T(:,:,:,9) = T(:,:,:,9)-1/3*trace
 da(id)%Tensor = T
  do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      do i=1,9
       Ts(id,i) = Ts(id,i)+T(l,m,n,i)
      end do 
     end do
   end do
 end do
  Ts = Ts/( (nx-1)*(ny-1)*(nz-1) )!%averaging over number of solution points
  da(id)%tt = Ts(id,:)
  do l=1,9
     Tss(l)= Tss(l)+Ts(id,l)
  end do
 end do

 call MPI_reduce(Tss, Tss0 , 9, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm1d, ierr)
 call MPI_Allreduce(ngridla, ngridaa , 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm1d, ierr)

 ! Tss0 = Tss0/ngridaa   !%averaging over number of elements.

if (myid==0) then
  Tss0 = Tss0/ngridaa
  write(*,*) 'done'
  write(24,*) ngridaa, (nx-1)*(ny-1)*(nz-1),numprocs
  write(24,*) '*********************************************'
  write(24,10) (Tss0(i), i=1,9),timea 
end if
 10 format (1x,9f12.8 ,1x, f12.8)
Return
end subroutine vcalculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smagvcalculation(timea, da, ngridla)
   use mpi
   use mpi_par
   use size
   use domain_definition
   use File_Units

   implicit DOUBLE PRECISION (a-h,o-z)
   type(domain), dimension(ngridla)            :: da
   integer                                     :: ngridla,id,ngridaa
   double precision                            :: timea
   double precision                            :: delta_g,norms
   double precision, dimension(nx,ny,nz,neq)   :: Qx , Qy, Qz
   double precision, dimension(nx,ny,nz,9)     :: s
   double precision, dimension(9)              :: sum2
   double precision, dimension(ngridla,9)      :: sum
   double precision, dimension(nx,ny,nz)       :: trace

   sum2 = 0.0d0
   sum = 0.0d0
   ! if (myid==0) then
   !  write(*,*) 'smagvcalculation has been called'
   ! end if
   !  if (myid==0) then
   !    write(*,*) '22smagvcalculation has been called'
   !  end if

   do id=1,ngridla

      do nv=1,neq
        call Gauss_Deriv(da(id)%gmetg,da(id)%Qlgg(:,:,:,nv),da(id)%Qglg(:,:,:,nv),&
                         da(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv),&
                         da(id)%ncg,da(id)%dmx,da(id)%dmy,da(id)%dmz)
      end do
    ! if (myid==0) then
    !   write(*,*) '32smagvcalculation has been called'
    ! end if
    pav = 0
   do l=1,da(id)%ncg(1)
      do m=1,da(id)%ncg(2)
        do n=1,da(id)%ncg(3)
          pav = pav + ( da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n) )
        end do
      end do
   end do
    
  ! if (myid==0) then
  !   write(25,*) pav
  ! end if 


    do l=1,da(id)%ncg(1)
      do m=1,da(id)%ncg(2)
        do n=1,da(id)%ncg(3)
          call filteredStraintensor(da(id),Qx,Qy,Qz,l,m,n,s,norms)
          
          s(l,m,n,:) = s(l,m,n,:) *norms*pav*delta_g( da(id),l,m,n)& 
                          &*delta_g(da(id),l,m,n) 
          trace      = 1.0d0/3.0d0*(s(l,m,n,1)+s(l,m,n,5)+s(l,m,n,9))
          s(l,m,n,1) = s(l,m,n,1)-tarce
          s(l,m,n,5) = s(l,m,n,5)-tarce
          s(l,m,n,9) = s(l,m,n,9)-tarce
          da(id)%Tensor = s

          do i=1,9
            sum(id,i) = sum(id,i)+s(l,m,n,i)
          end do
          end do
      end do
    end do
      if (myid==0) then
        write(25,*) 'Svalues'
        write(25,*) (s(2,2,2,i), i=1,9)
      end if
    sum  = sum /((nx-1)*(ny-1)*(nz-1))
    da(id)%tt = sum(id,:)
    sum2 = sum2 + sum(id,:)
   end do
    ! if (myid==0) then
    !   write(25,*) 'sum1'
    !   write(25,*) (sum(2,i), i=1,9)
    !   write(25,*) 'sum2'
    !   write(25,*) (sum2(i), i=1,9)
    ! end if
   call MPI_Allreduce(ngridla, ngridaa , 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm1d, ierr)
   sum2 = sum2/ngridaa
   ! if (myid==0) then
   !  write(25,*)  'averaged sum2'
   !  write(25,11) (sum2(i), i=1,9)
   ! end if
   11 format (1x,9f20.8 ,1x, f20.8)
 
end subroutine smagvcalculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filteredStraintensor(d,Qx,Qy,Qz,i,j,k,s,norms)
    use domain_definition
    implicit none
    type(domain)            :: d                                    !< Domain
    double precision, dimension(nx,ny,nz,9)         :: s 
    double precision, dimension(nx,ny,nz,neq)       :: Qx           !< Spectral derivative of d%Q with respect to x
    double precision, dimension(nx,ny,nz,neq)       :: Qy           !< Spectral derivative of d%Q with respect to y
    double precision, dimension(nx,ny,nz,neq)       :: Qz           !< Spectral derivative of d%Q with respect to z
    double precision        :: rho,u,v,w                            !< Instantaneous solutions
    double precision        :: norms                      
    double precision        :: rho_x, rho_y, rho_z                  !< 
    double precision        :: u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z  !< Velocity Derivatives
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
    s(i,j,k,1) = u_x
    s(i,j,k,2) = 0.50d0 * (u_y + v_x)
    s(i,j,k,3) = 0.50d0 * (u_z + w_x)
    s(i,j,k,4) = 0.50d0 * (v_x + u_y)
    s(i,j,k,5) = v_y
    s(i,j,k,6) = 0.50d0 * (v_z + w_y)
    s(i,j,k,7) = 0.50d0 * (w_x + u_z)
    s(i,j,k,8) = 0.50d0 * (w_y + v_z)
    s(i,j,k,9) = w_z
    norms   = sqrt( 2*s(i,j,k,1)*s(i,j,k,1) + 2*s(i,j,k,2)*s(i,j,k,2) + 2*s(i,j,k,3)*s(i,j,k,3) + &
                   &2*s(i,j,k,4)*s(i,j,k,4) + 2*s(i,j,k,5)*s(i,j,k,5) + 2*s(i,j,k,6)*s(i,j,k,6) + &
                   &2*s(i,j,k,7)*s(i,j,k,7) + 2*s(i,j,k,8)*s(i,j,k,8) + 2*s(i,j,k,9)*s(i,j,k,9)   )
    
end subroutine filteredStraintensor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!