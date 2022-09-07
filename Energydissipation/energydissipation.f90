
subroutine energydissipation(timea, da, ngridla)
   use mpi
   use mpi_par
   use size
   use constants
   use physics
   use domain_definition
   use input_data
   use File_Units
   implicit none

   type(domain), dimension(ngridla)            :: da
   integer                                     :: id, ngridla,l,m,n,nv
   double precision                            :: timea, sum, sum2, dissipatedenergy, Jacob
   double precision                            ::udery,uderz,vderx,vderz,wderx,wdery
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
               ! vorticity(l,m,n) = vorticity(l,m,n)*(2*pi/(11*(nx-1)))*(2*pi/(11*(ny-1)))*(2*pi/(11*(nz-1)))
   !------------------------------------------------------------------------------------------------------
              sum2             =sum2 +1.d0/da(id)%jacob(l,m,n)
              sum              =sum + vorticity(l,m,n)
            end do
         end do
      end do

   end do

   call MPI_reduce(sum, dissipatedenergy , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm1d, ierr)
   call MPI_reduce(sum2, Jacob , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm1d, ierr)
   if (myid == 0) then
      ! dissipatedenergy = (dissipatedenergy/((2*pi)**3)) * (1.d0/re)
      dissipatedenergy = 100.0d0 * (dissipatedenergy/Jacob) * (1.d0/re)
      write(20,*) timea,dissipatedenergy
   end if

end subroutine energydissipation

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

   implicit double precision(a-h,o-z)
   type(domain), dimension(ngridla)            :: da
   integer                                     :: ngridla,id,ngridaa,l,m,n
   double precision                            :: timea
   double precision, dimension(ngridla,9)      :: Ts
   double precision, dimension(9)              :: Tss,Tss0
   double precision, dimension(nx,ny,nz) :: uf,vf,wf,ur, vr,wr, uuf,uvf,uwf,vuf,vvf,vwf,wuf,wvf,wwf,uff,vff,wff,urf,vrf,wrf
   double precision, dimension(nx,ny,nz) :: uurf,uvrf,uwrf,vurf,vvrf,vwrf,wurf,wvrf,wwrf,uruf,urvf,urwf,vruf,vrvf,vrwf
   double precision, dimension(nx,ny,nz) :: wruf, wrvf, wrwf, ururf, urvrf,urwrf, vrurf,vrvrf,vrwrf,wrurf,wrvrf,wrwrf
   double precision, dimension(nx,ny,nz) :: trace
   double precision, dimension(nx,ny,nz,9) :: LL,MM,RR,T,sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!//////////////////////////////////////////////////////////////////////////////////////////////
Ts  = 0.0d0
Tss = 0.0d0
do id=1,ngridla
  puav = 0.0d0
  pvav = 0.0d0
  pwav = 0.0d0
  pav  = 0.0d0
!%%% Equation 6.45 of Dongru's Thesis %%%!
  do l=1,da(id)%ncg(3)
    do m=1,da(id)%ncg(2)
      do n=1,da(id)%ncg(1)
        puav= puav+da(id)%Q(l,m,n,2)*da(id)%jacob(l,m,n)
        pvav= pvav+da(id)%Q(l,m,n,3)*da(id)%jacob(l,m,n)
        pwav= pwav+da(id)%Q(l,m,n,4)*da(id)%jacob(l,m,n)
        pav = pav+da(id)%Q(l,m,n,1) *da(id)%jacob(l,m,n)
     end do
   end do
end do
! if(myid ==0) then
!   write(24,*) 'puav values '
!   write(24,*) puav, pvav, pwav, pav
! end if 

 uf(: ,:, :) = puav/pav
 vf(: ,:, :) = pvav/pav
 wf(: ,:, :) = pwav/pav

! if(myid ==0) then
!   write(24,*) 'uf values '
!   write(24,*) uf(1,1,1), vf(1,1,1), wf(1,1,1), uf(2,2,2), vf(2,2,2), wf(2,2,2)
! end if 

!%%%%% u =uf+ur !%%%%%
do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
        ur(l,m,n)=(da(id)%Q(l,m,n,2)/da(id)%Q(l,m,n,1))-uf(l,m,n)
        vr(l,m,n)=(da(id)%Q(l,m,n,3)/da(id)%Q(l,m,n,1))-vf(l,m,n)
        wr(l,m,n)=(da(id)%Q(l,m,n,4)/da(id)%Q(l,m,n,1))-wf(l,m,n)
     end do
   end do
end do

! if(myid ==0) then
!   write(24,*) 'u values'
!   write(24,*) da(id)%Q(1,1,1,2)/da(id)%Q(1,1,1,1), da(id)%Q(1,1,1,3)/da(id)%Q(1,1,1,1), da(id)%Q(1,1,1,4)/da(id)%Q(1,1,1,1)
! end if 

! if(myid ==0) then
!   write(24,*) 'ur values'
!   write(24,*) ur(1,1,1), vr(1,1,1), wr(1,1,1), ur(2,2,2), vr(2,2,2), wr(2,2,2)
! end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!% computation of Leonard Stress tensor, Equation 6-51 %!!!!
puuav = 0.0d0
puvav = 0.0d0
puwav = 0.0d0
pvuav = 0.0d0
pvvav = 0.0d0
pvwav = 0.0d0
pwuav = 0.0d0 
pwvav = 0.0d0
pwwav = 0.0d0
!!!!
pufav = 0.0d0
pvfav = 0.0d0
pwfav = 0.0d0
purav = 0.0d0
pvrav = 0.0d0
pwrav = 0.0d0

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      puuav= puuav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)*uf(l,m,n)
      puvav= puvav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)*vf(l,m,n)
      puwav= puwav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)*wf(l,m,n)
      pvuav= pvuav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)*uf(l,m,n)
      pvvav= pvvav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)*vf(l,m,n)
      pvwav= pvwav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)*wf(l,m,n)
      pwuav= pwuav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)*uf(l,m,n)
      pwvav= pwvav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)*vf(l,m,n)
      pwwav= pwwav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)*wf(l,m,n)

!   
      pufav=pufav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)
      pvfav=pvfav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)
      pwfav=pwfav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)
      purav=purav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)
      pvrav=pvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)
      pwrav=pwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)
      end do
   end do
end do


! if(myid ==0 ) then
!   write(24,*) 'puuav             ','puvav            ','puwav           ','pvuav          '
!   write(24,*) puuav, puvav, puwav, pvuav, pvvav, pvwav, pwuav, pwvav, pwwav
!  end if 


uuf(: , : , :) = puuav/pav
uvf(: , : , :) = puvav/pav
uwf(: , : , :) = puwav/pav
vuf(: , : , :) = pvuav/pav
vvf(: , : , :) = pvvav/pav
vwf(: , : , :) = pvwav/pav
wuf(: , : , :) = pwuav/pav
wvf(: , : , :) = pwvav/pav
wwf(: , : , :) = pwwav/pav
!
uff(: , : , :) = pufav/pav
vff(: , : , :) = pvfav/pav
wff(: , : , :) = pwfav/pav
urf(: , : , :) = purav/pav
vrf(: , : , :) = pvrav/pav
wrf(: , : , :) = pwrav/pav
! if(myid ==0) then
!   write(24,*) 'uuf values'
!   write(24,*) uuf(1,1,1), uvf(1,1,1), uwf(1,1,1), uuf(2,2,2), vuf(2,2,2), uwf(2,2,2)
! end if 

! if(myid ==0) then
!   write(24,*) 'uff values'
!   write(24,*) uff(1,1,1), vff(1,1,1), wff(1,1,1)
! end if 

! if(myid ==0) then
!   write(24,*) 'urf values'
!    do l=1,da(id)%ncg(3)
!     do m=1,da(id)%ncg(2)
!      do n=1,da(id)%ncg(1)
!        write(24,*) da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n),ur(l,m,n)
!      end do
!     end do
!   end do
!   write(24,*) da(id)%Q(1,1,1,1), da(id)%jacob(1,1,1),ur(1,1,1), purav, pav, urf(1,1,1), vrf(1,1,1), wrf(1,1,1),&
!               & urf(3,3,3), vrf(3,3,3), wrf(3,3,3)
! end if 

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
       LL(l,m,n,1)= pav*(uuf(l,m,n)-uff(l,m,n)*uff(l,m,n))
       LL(l,m,n,2)= pav*(uvf(l,m,n)-uff(l,m,n)*vff(l,m,n))
       LL(l,m,n,3)= pav*(uwf(l,m,n)-uff(l,m,n)*wff(l,m,n))
       LL(l,m,n,4)= pav*(vuf(l,m,n)-vff(l,m,n)*uff(l,m,n))
       LL(l,m,n,5)= pav*(vvf(l,m,n)-vff(l,m,n)*vff(l,m,n))
       LL(l,m,n,6)= pav*(vwf(l,m,n)-vff(l,m,n)*wff(l,m,n))
       LL(l,m,n,7)= pav*(wuf(l,m,n)-wff(l,m,n)*uff(l,m,n))
       LL(l,m,n,8)= pav*(wvf(l,m,n)-wff(l,m,n)*vff(l,m,n))
       LL(l,m,n,9)= pav*(wwf(l,m,n)-wff(l,m,n)*wff(l,m,n))
     end do
   end do
end do
if(myid ==0) then
  write(24,*) 'LL values'
  write(24,*) LL(2,2,2,1), LL(2,2,2,2), LL(2,2,2,3), LL(2,2,2,4), LL(2,2,2,5), LL(2,2,2,6), LL(2,2,2,7), LL(2,2,2,8), LL(2,2,2,9)
end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!/////computes the cross tensor!  equation 6-52////
puurav = 0.0d0
puvrav = 0.0d0
puwrav = 0.0d0
pvurav = 0.0d0
pvvrav = 0.0d0
pvwrav = 0.0d0
pwurav = 0.0d0 
pwvrav = 0.0d0
pwwrav = 0.0d0
!
puruav = 0.0d0
purvav = 0.0d0
purwav = 0.0d0
pvruav = 0.0d0
pvrvav = 0.0d0
pvrwav = 0.0d0
pwruav = 0.0d0 
pwrvav = 0.0d0
pwrwav = 0.0d0
 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      puurav= puurav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)*ur(l,m,n)
      puvrav= puvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)*vr(l,m,n)
      puwrav= puwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*uf(l,m,n)*wr(l,m,n)
      pvurav= pvurav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)*ur(l,m,n)
      pvvrav= pvvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)*vr(l,m,n)
      pvwrav= pvwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vf(l,m,n)*wr(l,m,n)
      pwurav= pwurav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)*ur(l,m,n)
      pwvrav= pwvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)*vr(l,m,n)
      pwwrav= pwwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wf(l,m,n)*wr(l,m,n)
!   
      puruav= puruav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)*uf(l,m,n)
      purvav= purvav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)*vf(l,m,n)
      purwav= purwav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)*wf(l,m,n)
      pvruav= pvruav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)*uf(l,m,n)
      pvrvav= pvrvav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)*vf(l,m,n)
      pvrwav= pvrwav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)*wf(l,m,n)
      pwruav= pwruav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)*uf(l,m,n)
      pwrvav= pwrvav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)*vf(l,m,n)
      pwrwav= pwrwav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)*wf(l,m,n)

     end do
   end do
end do
uurfarv = puurav/pav
uvrfarv = puvrav/pav
uwrfarv = puwrav/pav
vurfarv = pvurav/pav
vvrfarv = pvvrav/pav
vwrfarv = pvwrav/pav
wurfarv = pwurav/pav
wvrfarv = pwvrav/pav
wwrfarv = pwwrav/pav
!
urufarv = puruav/pav
urvfarv = purvav/pav
urwfarv = purwav/pav
vrufarv = pvruav/pav
vrvfarv = pvrvav/pav
vrwfarv = pvrwav/pav
wrufarv = pwruav/pav
wrvfarv = pwrvav/pav
wrwfarv = pwrwav/pav



uurf(: , : , :) = uurfarv
uvrf(: , : , :) = uvrfarv
uwrf(: , : , :) = uwrfarv 
vurf(: , : , :) = vurfarv
vvrf(: , : , :) = vvrfarv
vwrf(: , : , :) = vwrfarv
wurf(: , : , :) = wurfarv
wvrf(: , : , :) = wvrfarv
wwrf(: , : , :) = wwrfarv
!
uruf(: , : , :) = urufarv
urvf(: , : , :) = urvfarv
urwf(: , : , :) = urwfarv 
vruf(: , : , :) = vrufarv
vrvf(: , : , :) = vrvfarv
vrwf(: , : , :) = vrwfarv
wruf(: , : , :) = wrufarv
wrvf(: , : , :) = wrvfarv
wrwf(: , : , :) = wrwfarv

! if (myid==0) then
!   write(24,*) 'uur values'
!   write(24,*) '**********************'
!   write(24,*) uf(1,1,1),   ur(1,1,1)
!   write(24,*) uurf(2,2,2), uvrf(2,2,2), uwrf(2,2,2), vurf(2,2,2), vvrf(2,2,2), vwrf(2,2,2), wurf(2,2,2), wvrf(2,2,2), wwrf(2,2,2)
!   write(24,*) uurf(3,3,3), uvrf(3,3,3), uwrf(3,3,3), vurf(3,3,3), vvrf(3,3,3), vwrf(3,3,3), wurf(3,3,3), wvrf(3,3,3), wwrf(3,3,3)
! end if 

! if (myid==0) then
!   write(24,*) 'uru values'
!   write(24,*) uruf(2,2,2), urvf(2,2,2), urwf(2,2,2), vruf(2,2,2), vrvf(2,2,2), vrwf(2,2,2), wruf(2,2,2), wrvf(2,2,2), wrwf(2,2,2)
!   write(24,*) uruf(3,3,3), urvf(3,3,3), urwf(3,3,3), vruf(3,3,3), vrvf(3,3,3), vrwf(3,3,3), wruf(3,3,3), wrvf(3,3,3), wrwf(3,3,3)
! end if 


 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
       MM(l,m,n,1)= pav*(uurf(l,m,n)+uruf(l,m,n)-uff(l,m,n)*urf(l,m,n)-urf(l,m,n)*uff(l,m,n))
       MM(l,m,n,2)= pav*(uvrf(l,m,n)+urvf(l,m,n)-uff(l,m,n)*vrf(l,m,n)-urf(l,m,n)*vff(l,m,n))
       MM(l,m,n,3)= pav*(uwrf(l,m,n)+wurf(l,m,n)-uff(l,m,n)*wrf(l,m,n)-urf(l,m,n)*wff(l,m,n))
       MM(l,m,n,4)= pav*(vurf(l,m,n)+uvrf(l,m,n)-vff(l,m,n)*urf(l,m,n)-vrf(l,m,n)*uff(l,m,n))
       MM(l,m,n,5)= pav*(vvrf(l,m,n)+vrvf(l,m,n)-vff(l,m,n)*vrf(l,m,n)-vrf(l,m,n)*vff(l,m,n))
       MM(l,m,n,6)= pav*(vwrf(l,m,n)+wvrf(l,m,n)-vff(l,m,n)*wrf(l,m,n)-vrf(l,m,n)*wff(l,m,n))
       MM(l,m,n,7)= pav*(wurf(l,m,n)+uwrf(l,m,n)-urf(l,m,n)*wff(l,m,n)-wrf(l,m,n)*uff(l,m,n))
       MM(l,m,n,8)= pav*(wvrf(l,m,n)+vwrf(l,m,n)-wff(l,m,n)*vrf(l,m,n)-wrf(l,m,n)*vff(l,m,n))
       MM(l,m,n,9)= pav*(wwrf(l,m,n)+wrwf(l,m,n)-wff(l,m,n)*wrf(l,m,n)-wrf(l,m,n)*wff(l,m,n))
     end do
   end do
end do

if (myid==0) then
  write(24,*) 'MM values'
  write(24,*) MM(2,2,2,1), MM(2,2,2,2), MM(2,2,2,3), MM(2,2,2,4), MM(2,2,2,5), MM(2,2,2,6), MM(2,2,2,7), MM(2,2,2,8), MM(2,2,2,9)
  write(24,*) MM(3,3,3,1), MM(3,3,3,2), MM(3,3,3,3), MM(3,3,3,4), MM(3,3,3,5), MM(3,3,3,6), MM(3,3,3,7), MM(3,3,3,8), MM(3,3,3,9)
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!/ computation of Reynolds strees Equation 6-53/!
pururav = 0.0d0
purvrav = 0.0d0
purwrav = 0.0d0
pvrurav = 0.0d0
pvrvrav = 0.0d0
pvrwrav = 0.0d0
pwrurav = 0.0d0 
pwrvrav = 0.0d0
pwrwrav = 0.0d0
 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      pururav= pururav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)*ur(l,m,n)
      purvrav= purvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)*vr(l,m,n)
      purwrav= purwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n)*wr(l,m,n)
      pvrurav= pvrurav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)*ur(l,m,n)
      pvrvrav= pvrvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)*vr(l,m,n)
      pvrwrav= pvrwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*vr(l,m,n)*wr(l,m,n)
      pwrurav= pwrurav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)*ur(l,m,n)
      pwrvrav= pwrvrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)*vr(l,m,n)
      pwrwrav= pwrwrav+da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*wr(l,m,n)*wr(l,m,n)
     end do
   end do
 end do
ururf(: , : , :) = pururav/pav
urvrf(: , : , :) = purvrav/pav
urwrf(: , : , :) = purwrav/pav
vrurf(: , : , :) = pvrurav/pav
vrvrf(: , : , :) = pvrvrav/pav
vrwrf(: , : , :) = pvrwrav/pav
wrurf(: , : , :) = pwrurav/pav
wrvrf(: , : , :) = pwrvrav/pav
wrwrf(: , : , :) = pwrwrav/pav

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
       RR(l,m,n,1)= pav*( (ururf(l,m,n)-( urf(l,m,n)*urf(l,m,n) ) ) )
       RR(l,m,n,2)= pav*( (urvrf(l,m,n)-( urf(l,m,n)*vrf(l,m,n) ) ) )
       RR(l,m,n,3)= pav*( (urwrf(l,m,n)-( urf(l,m,n)*wrf(l,m,n) ) ) )
       RR(l,m,n,4)= pav*( (vrurf(l,m,n)-( vrf(l,m,n)*urf(l,m,n) ) ) )
       RR(l,m,n,5)= pav*( (vrvrf(l,m,n)-( vrf(l,m,n)*vrf(l,m,n) ) ) )
       RR(l,m,n,6)= pav*( (vrwrf(l,m,n)-( vrf(l,m,n)*wrf(l,m,n) ) ) )
       RR(l,m,n,7)= pav*( (wrurf(l,m,n)-( wrf(l,m,n)*urf(l,m,n) ) ) )
       RR(l,m,n,8)= pav*( (wrvrf(l,m,n)-( wrf(l,m,n)*vrf(l,m,n) ) ) )
       RR(l,m,n,9)= pav*( (wrwrf(l,m,n)-( wrf(l,m,n)*wrf(l,m,n) ) ) )
     end do
   end do
 end do
 if (myid==0) then
  write(24,*) 'RR values'
  write(24,*) RR(2,2,2,1), RR(2,2,2,2), RR(2,2,2,3), RR(2,2,2,4), RR(2,2,2,5), RR(2,2,2,6), RR(2,2,2,7), RR(2,2,2,8), RR(2,2,2,9)
  write(24,*) RR(3,3,3,1), RR(3,3,3,2), RR(3,3,3,3), RR(3,3,3,4), RR(3,3,3,5), RR(3,3,3,6), RR(3,3,3,7), RR(3,3,3,8), RR(3,3,3,9)
end if 
 T   = 0.0d0
 trace = 0.0d0
 T = LL+MM+RR
 trace = ( T(:,:,:,1)+T(:,:,:,5)+T(:,:,:,9) )
 T(:,:,:,1) = T(:,:,:,1)-1/3*trace
 T(:,:,:,5) = T(:,:,:,5)-1/3*trace
 T(:,:,:,9) = T(:,:,:,9)-1/3*trace
  ! da(id)%Tensor = T
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
  ! da(id)%tt = Ts(id,:)
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

   implicit double precision(a-h,o-z)
   type(domain), dimension(ngridla)            :: da
   integer                                     :: ngridla,id,nv,i, l,m,n,ngridaa
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
          ! da(id)%Tensor = s

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
    ! da(id)%tt = sum(id,:)
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
subroutine nodalfilter(da,ngridla,timea)

  use size
  use domain_definition
  use mpi
  use mpi_par
  implicit none
  integer                            :: id, ngridla
  double precision                   ::timea
  double precision                   ::alpha
  type(domain) , dimension(ngridla)  ::da
  type (domain)                      ::d
  double precision, dimension(nx,nx) ::cx,cy,cz         !%interpolating to lower p matrices
  double precision, dimension(nx,nx) ::ccx,ccy,ccz      !%interpolating back matrices
  double precision, dimension(nx,nx,nx,neq)::Qs

   call zero_domain(d)
   d%nc  = nx-2
   d%ncg = nx-3
   call set_points(d)
   CALL intrpmat(nx,da(1)%ncg(1),da(1)%cxg(:,1),d%ncg(1),d%cxg(:,1),cx(:,:)) 
   CALL intrpmat(ny,da(1)%ncg(2),da(1)%cxg(:,2),d%ncg(2),d%cxg(:,2),cy(:,:)) 
   CALL intrpmat(nz,da(1)%ncg(3),da(1)%cxg(:,3),d%ncg(3),d%cxg(:,3),cz(:,:)) 

   CALL intrpmat(nx,d%ncg(1),d%cxg(:,1),da(1)%ncg(1),da(1)%cxg(:,1),ccx(:,:))
   CALL intrpmat(ny,d%ncg(2),d%cxg(:,2),da(1)%ncg(2),da(1)%cxg(:,2),ccy(:,:)) 
   CALL intrpmat(nz,d%ncg(3),d%cxg(:,3),da(1)%ncg(3),da(1)%cxg(:,3),ccz(:,:))  


    do id = 1,ngridla
       Qs = da(id)%Q(:,:,:,:)
      call interpx( nx, cx, da(id)%Q(:,:,:,:), da(1)%ncg, d%Q(:,:,:,:), d%ncg )
      call interpy( ny, cy, da(id)%Q(:,:,:,:), da(1)%ncg, d%Q(:,:,:,:), d%ncg )
      call interpzz( nz, cz, da(id)%Q(:,:,:,:), da(1)%ncg, d%Q(:,:,:,:), d%ncg )

      call interpx( nx, ccx, d%Q(:,:,:,:), d%ncg, da(id)%Q(:,:,:,:), da(1)%ncg )
      call interpy( ny, ccy, d%Q(:,:,:,:), d%ncg, da(id)%Q(:,:,:,:), da(1)%ncg )
      call interpzz( nz, ccz, d%Q(:,:,:,:), d%ncg, da(id)%Q(:,:,:,:), da(1)%ncg)

      da(id)%Q = (1-0.01)*Qs + 0.01*da(id)%Q
      ! da(id)%Q = da(id)%Q
   end do
   if (myid==0 ) write(*,*) 'filtering is done'
   Return
end subroutine nodalfilter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE interpzz(ncol,b,u,ncg,ud,nc)  
      USE size
      implicit none
      integer               ::nv, k,ncol 
      INTEGER, DIMENSION(3) :: nc,ncg
      DOUBLE PRECISION      :: b(ny,ny),u(nx,ny,nz,neq),ud(nx,ny,nz,neq),dmatt(ny,ny) 

      
      CALL transpose( dmatt(1:ncg(3),1:nc(3)) ,ncg(3), b(1:nc(3),1:ncg(3)) ,nc(3))    
      DO nv=1,neq
      DO k=1,ncg(2)                 
         CALL mxm(u(1:ncg(1),k,1:ncg(3),nv),ncg(1),dmatt(1:ncg(3),1:nc(3)),ncg(3),ud(1:ncg(1),k,1:nc(3),nv),nc(3))
      ENDDO

      ENDDO

      RETURN 
   END SUBROUTINE interpzz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine modalfilter(da,ngridla,timea,step)
    use size
    use constants
    use domain_definition
    USE, INTRINSIC :: iso_c_binding
    use mpi
    use mpi_par
    use File_Units
    implicit none
    type (domain), dimension(ngridla)   :: da
    integer                             :: ngridla
    double precision                    :: timea,maximum
    double precision,dimension(9,3)       :: former
    INCLUDE 'fftw3.f03'

    INTEGER          :: id, nv
    integer          :: i, j, k, N, m, ii, jj, kk,step
!   double precision :: x1, x2, x3, pi=3.14159265358979323846264338327950
    type (C_PTR)     :: plan1, plan2
    real (C_DOUBLE), dimension(:,:,:), allocatable :: in, mid, out,fout
!
    N = da(1)%ncg(1)
    allocate (in(N,N,N),mid(N,N,N),out(N,N,N),fout(N,N,N))
      
    !------------------------ ----Prepare FFTW plans-----------------------------------------

    plan1 = fftw_plan_r2r_3d(N,N,N,in,mid,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE)
    plan2 = fftw_plan_r2r_3d(N,N,N,mid,out,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE)
      
    do id = 1, ngridla
      do nv = 2, 4
          
      !------------------------Transform from nodal to modal-----------------------------------
        in(:,:,:) = da(id)%Q(1:N,1:N,1:N,nv)*da(id)%jacob(1,1,1)
        ! fout(:,:,:) = in(:,:,:)  
 
          
        call fftw_execute_r2r(plan1, in, mid)
   
      do i = 1, N
         do j = 1, N
          mid(1,i,j) = mid(1,i,j) / (2.*N)
          mid(i,1,j) = mid(i,1,j) / (2.*N)
          mid(i,j,1) = mid(i,j,1) / (2.*N)
         end do
      end do
      do i = 1, N
         do j = 1, N
          do k = 2, N
           mid(i,j,k) = mid(i,j,k) / N
           mid(j,k,i) = mid(j,k,i) / N
           mid(k,i,j) = mid(k,i,j) / N
          end do
         end do
      end do
      
      !---------------------------------- Make changes in modes-------------------------------
      
      ! do m = 8,9
      !     mid(m,:,:) = 0.0d0
      !     mid(:,m,:) = 0.0d0
      !     mid(:,:,m) = 0.0d0
      !    if (myid==29 .and. id==12 .and. nv==4 ) write(*,*) 'third is done'
      ! end do

      !------------------------------ Transform back from modal to nodal----------------------
   
      do i = 1, N
         do j = 1, N
          mid(1,i,j) = mid(1,i,j)
          mid(i,1,j) = mid(i,1,j)
          mid(i,j,1) = mid(i,j,1)
         end do
      end do
      do i = 1, N
         do j = 1, N
          do k = 2, N
           mid(i,j,k) = mid(i,j,k) / 2
           mid(j,k,i) = mid(j,k,i) / 2
           mid(k,i,j) = mid(k,i,j) / 2
          end do
         end do
      end do

       if(myid==50 .and. id==108) then
       ! if (step == 50) then

        ! if(nv==2) then
        !   write(30,*) timea,(mid(i,i,i),i=1,N)
        !   ! do i=1,N
        !   !   former(i,nv-1) = mid(i,i,i)
        !   ! end do
        ! end if

        ! if(nv==3) then
        !   write(31,*) timea,(mid(i,i,i),i=1,N)
        !   ! do i=1,N
        !   !   former(i,nv-1) = mid(i,i,i)
        !   ! end do
        ! end if

        ! if(nv==4) then
        !   write(32,*) timea,(mid(i,i,i),i=1,N)
        !   ! do i=1,N
        !   !   former(i,nv-1) = mid(i,i,i)
        !   ! end do
        ! end if

       ! else

        if(nv==2) then
          write(30,*) timea,(mid(i,i,i),i=1,N)
        end if

        if(nv==3) then
          write(31,*) timea,(mid(i,i,i),i=1,N)
        end if

        if(nv==4) then
          write(32,*) timea,(mid(i,i,i),i=1,N)
        end if

       

       maximum = mid(1,1,1)
       ii = 1
       jj = 1
       kk = 1
       do i=1,N
        do j=1,N
         do k=1,N
           if (abs(mid(i,j,k)) > abs(maximum)) then
            maximum = mid(i,j,k) 
            ii = i
            jj = j
            kk = k
           end if
        end do
       end do
      end do

        if(nv==2) write(33,*) timea,ii,jj,kk,da(id)%xg(1,ii,jj,kk),da(id)%xg(2,ii,jj,kk),da(id)%xg(3,ii,jj,kk),mid(ii,jj,kk)
        if(nv==3) write(34,*) timea,ii,jj,kk,da(id)%xg(1,ii,jj,kk),da(id)%xg(2,ii,jj,kk),da(id)%xg(3,ii,jj,kk),mid(ii,jj,kk)
        if(nv==4) write(35,*) timea,ii,jj,kk,da(id)%xg(1,ii,jj,kk),da(id)%xg(2,ii,jj,kk),da(id)%xg(3,ii,jj,kk),mid(ii,jj,kk)
      if(nv==3 .and. step==50) then
        write(43,*) da(id)%xg(1,1,1,1)
        write(43,*) da(id)%xg(2,1,1,1)
        write(43,*) da(id)%xg(3,1,1,1)
      end if

      end if
      
      ! end if

      !  if(myid==23 .and. id==56) then
      !  ! if (step == 50) then

      !   if(nv==2) then
      !     write(36,*) timea,(mid(i,i,i),i=1,N)
      !   end if

      !   if(nv==3) then
      !     write(37,*) timea,(mid(i,i,i),i=1,N)
      !   end if

      !   if(nv==4) then
      !     write(38,*) timea,(mid(i,i,i),i=1,N)
      !   end if

      !  ! else

      !  !  if(nv==2) then
      !  !    write(30,*) timea,(mid(i,i,i)-former(i,nv-1),i=1,N)
      !  !    do i=1,N
      !  !      former(i,nv-1) = mid(i,i,i)
      !  !    end do
      !  !  end if

      !  !  if(nv==3) then
      !  !    write(31,*) timea,(mid(i,i,i)-former(i,nv-1),i=1,N)
      !  !    do i=1,N
      !  !      former(i,nv-1) = mid(i,i,i)
      !  !    end do
      !  !  end if

      !  !  if(nv==4) then
      !  !    write(32,*) timea,(mid(i,i,i)-former(i,nv-1),i=1,N)
      !  !    do i=1,N
      !  !      former(i,nv-1) = mid(i,i,i)
      !  !    end do
      !  !  end if

      !  ! end if

      !  maximum = mid(1,1,1)
      !  ii = 1
      !  jj = 1
      !  kk = 1
      !  do i=1,N
      !   do j=1,N
      !    do k=1,N
      !      if (abs(mid(i,j,k)) > abs(maximum)) then
      !       maximum = mid(i,j,k) 
      !       ii = i
      !       jj = j
      !       kk = k
      !      end if
      !   end do
      !  end do
      ! end do

      !   if(nv==2) write(39,*) timea,ii,jj,kk,da(id)%xg(1,ii,jj,kk),da(id)%xg(2,ii,jj,kk),da(id)%xg(3,ii,jj,kk),mid(ii,jj,kk)
      !   if(nv==3) write(40,*) timea,ii,jj,kk,da(id)%xg(1,ii,jj,kk),da(id)%xg(2,ii,jj,kk),da(id)%xg(3,ii,jj,kk),mid(ii,jj,kk)
      !   if(nv==4) write(41,*) timea,ii,jj,kk,da(id)%xg(1,ii,jj,kk),da(id)%xg(2,ii,jj,kk),da(id)%xg(3,ii,jj,kk),mid(ii,jj,kk)
      ! if(nv==3 .and. step==50) then
      !   write(42,*) da(id)%xg(1,1,1,1)
      !   write(42,*) da(id)%xg(2,1,1,1)
      !   write(42,*) da(id)%xg(3,1,1,1)
      ! end if
      
      ! end if

   
      call fftw_execute_r2r(plan2, mid, out)
      
      ! da(id)%Q(1:N,1:N,1:N,nv) = (0.01*out(:,:,:)+0.99*fout(:,:,:))/da(id)%jacob(1,1,1)
      da(id)%Q(1:N,1:N,1:N,nv) = out(:,:,:)/da(id)%jacob(1,1,1)
      
      end do
    end do
      
    
      
      !--------------------------------------Destroy FFTW plans----------------------
   
    call fftw_destroy_plan(plan1)
    call fftw_destroy_plan(plan2)
    if (myid==29) write(*,*) 'filtering is done'

END SUBROUTINE modalfilter
