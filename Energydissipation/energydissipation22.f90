
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
   double precision, dimension(nx,ny,nz) :: uf,vf,wf,ur, vr,wr, uuf,uvf,uwf,vuf,vvf,vwf,wuf,wvf,wwf,uff,vff,wff,urf,vrf,wrf
   double precision, dimension(nx,ny,nz) :: uurf,uvrf,uwrf,vurf,vvrf,vwrf,wurf,wvrf,wwrf,uruf,urvf,urwf,vruf,vrvf,vrwf
   double precision, dimension(nx,ny,nz) :: wruf, wrvf, wrwf, ururf, urvrf,urwrf, vrurf,vrvrf,vrwrf,wrurf,wrvrf,wrwrf
   double precision, dimension(nx,ny,nz) :: trace
   double precision, dimension(nx,ny,nz) :: count
   double precision, dimension(nx,ny,nz) :: pav,puav,pvav,pwav
   double precision, dimension(nx,ny,nz) :: puuav,puvav,puwav,pvuav,pvvav,pvwav,pwuav,pwvav,pwwav
   double precision, dimension(nx,ny,nz) :: pufav,pvfav,pwfav,purav,pvrav,pwrav
   double precision, dimension(nx,ny,nz) :: puurav,puvrav,puwrav,pvurav,pvvrav,pvwrav,pwurav,pwvrav,pwwrav
   double precision, dimension(nx,ny,nz) :: puruav,purvav,purwav,pvruav,pvrvav,pvrwav,pwruav,pwrvav,pwrwav
   double precision, dimension(nx,ny,nz) :: pururav,purvrav,purwrav,pvrurav,pvrvrav,pvrwrav,pwrurav,pwrvrav,pwrwrav
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
        puav(l,m,n) = da(id)%Q(l,m,n,2)*da(id)%jacob(l,m,n)
        pvav(l,m,n) = da(id)%Q(l,m,n,3)*da(id)%jacob(l,m,n)
        pwav(l,m,n) = da(id)%Q(l,m,n,4)*da(id)%jacob(l,m,n)
        pav(l,m,n)  = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)
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
                    puav(l,m,n) = puav(l,m,n)+da(id)%Q(p,q,r,2)*da(id)%jacob(p,q,r)
                    pvav(l,m,n) = pvav(l,m,n)+da(id)%Q(p,q,r,3)*da(id)%jacob(p,q,r)
                    pwav(l,m,n) = pwav(l,m,n)+da(id)%Q(p,q,r,4)*da(id)%jacob(p,q,r)
                    pav(l,m,n)  = pav(l,m,n)+da(id)%Q(p,q,r,1)*da(id)%jacob(p,q,r)
                    count(l,m,n) = count(l,m,n)+1
                end if 
            end do
          end do
        end do
        uf(l,m,n) = puav(l,m,n)/pav(l,m,n)
        vf(l,m,n) = pvav(l,m,n)/pav(l,m,n)
        wf(l,m,n) = pwav(l,m,n)/pav(l,m,n)
      end do
    end do
  end do
if(myid ==0) then
  write(24,*) 'puav values '
  write(24,*) puav(1,1,1), pvav(1,1,1), pwav(1,1,1), pav(1,1,1), count(1,1,1)
  write(24,*) puav(2,2,2), pvav(2,2,2), pwav(2,2,2), pav(2,2,2), count(2,2,2)
end if 

if(myid ==0) then
  write(24,*) 'uf values '
  write(24,*) uf(1,1,1), vf(1,1,1), wf(1,1,1), uf(2,2,2), vf(2,2,2), wf(2,2,2)
end if 

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

if(myid ==0) then
  write(24,*) 'u values'
  write(24,*) da(id)%Q(1,1,1,2)/da(id)%Q(1,1,1,1), da(id)%Q(1,1,1,3)/da(id)%Q(1,1,1,1), da(id)%Q(1,1,1,4)/da(id)%Q(1,1,1,1)
end if 

if(myid ==0) then
  write(24,*) 'ur values'
  write(24,*) ur(1,1,1), vr(1,1,1), wr(1,1,1), ur(2,2,2), vr(2,2,2), wr(2,2,2)
end if 
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
      rho = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)
      puuav(l,m,n) = rho*uf(l,m,n)*uf(l,m,n)
      puvav(l,m,n) = rho*uf(l,m,n)*vf(l,m,n)
      puwav(l,m,n) = rho*uf(l,m,n)*wf(l,m,n)
      pvuav(l,m,n) = rho*vf(l,m,n)*uf(l,m,n)
      pvvav(l,m,n) = rho*vf(l,m,n)*vf(l,m,n)
      pvwav(l,m,n) = rho*vf(l,m,n)*wf(l,m,n)
      pwuav(l,m,n) = rho*wf(l,m,n)*uf(l,m,n)
      pwvav(l,m,n) = rho*wf(l,m,n)*vf(l,m,n)
      pwwav(l,m,n) = rho*wf(l,m,n)*wf(l,m,n)

!   
      pufav(l,m,n) = rho*uf(l,m,n)
      pvfav(l,m,n) = rho*vf(l,m,n)
      pwfav(l,m,n) = rho*wf(l,m,n)
      purav(l,m,n) = rho*ur(l,m,n)
      pvrav(l,m,n) = rho*vr(l,m,n)
      pwrav(l,m,n) = rho*wr(l,m,n)
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
              puuav(l,m,n) = puuav(l,m,n)+rhor*uf(p,q,r)*uf(p,q,r)
              puvav(l,m,n) = puvav(l,m,n)+rhor*uf(p,q,r)*vf(p,q,r)
              puwav(l,m,n) = puwav(l,m,n)+rhor*uf(p,q,r)*wf(p,q,r)
              pvuav(l,m,n) = pvuav(l,m,n)+rhor*vf(p,q,r)*uf(p,q,r)
              pvvav(l,m,n) = pvvav(l,m,n)+rhor*vf(p,q,r)*vf(p,q,r)
              pvwav(l,m,n) = pvwav(l,m,n)+rhor*vf(p,q,r)*wf(p,q,r)
              pwuav(l,m,n) = pwuav(l,m,n)+rhor*wf(p,q,r)*uf(p,q,r)
              pwvav(l,m,n) = pwvav(l,m,n)+rhor*wf(p,q,r)*vf(p,q,r)
              pwwav(l,m,n) = pwwav(l,m,n)+rhor*wf(p,q,r)*wf(p,q,r)

!   
              pufav(l,m,n) = pufav(l,m,n)+rhor*uf(p,q,r)
              pvfav(l,m,n) = pvfav(l,m,n)+rhor*vf(p,q,r)
              pwfav(l,m,n) = pwfav(l,m,n)+rhor*wf(p,q,r)
              purav(l,m,n) = purav(l,m,n)+rhor*ur(p,q,r)
              pvrav(l,m,n) = pvrav(l,m,n)+rhor*vr(p,q,r)
              pwrav(l,m,n) = pwrav(l,m,n)+rhor*wr(p,q,r)
            end if 
          end do
        end do
      end do
      uuf(l,m,n) = puuav(l,m,n)/pav(l,m,n)
      uvf(l,m,n) = puvav(l,m,n)/pav(l,m,n)
      uwf(l,m,n) = puwav(l,m,n)/pav(l,m,n)
      vuf(l,m,n) = pvuav(l,m,n)/pav(l,m,n)
      vvf(l,m,n) = pvvav(l,m,n)/pav(l,m,n)
      vwf(l,m,n) = pvwav(l,m,n)/pav(l,m,n)
      wuf(l,m,n) = pwuav(l,m,n)/pav(l,m,n)
      wvf(l,m,n) = pwvav(l,m,n)/pav(l,m,n)
      wwf(l,m,n) = pwwav(l,m,n)/pav(l,m,n)
!
      uff(l,m,n) = pufav(l,m,n)/pav(l,m,n)
      vff(l,m,n) = pvfav(l,m,n)/pav(l,m,n)
      wff(l,m,n) = pwfav(l,m,n)/pav(l,m,n)
      urf(l,m,n) = purav(l,m,n)/pav(l,m,n)
      vrf(l,m,n) = pvrav(l,m,n)/pav(l,m,n)
      wrf(l,m,n) = pwrav(l,m,n)/pav(l,m,n)
      end do
   end do
end do


! if(myid ==0 ) then
!   write(24,*) 'puuav             ','puvav            ','puwav           ','pvuav          '
!   write(24,*) puuav, puvav, puwav, pvuav, pvvav, pvwav, pwuav, pwvav, pwwav
!  end if 


if(myid ==0) then
  write(24,*) 'uuf values'
  write(24,*) uuf(1,1,1), uvf(1,1,1), uwf(1,1,1), uuf(2,2,2), vuf(2,2,2), uwf(2,2,2)
end if 

if(myid ==0) then
  write(24,*) 'uff values'
  write(24,*) uff(1,1,1), vff(1,1,1), wff(1,1,1)
end if 

if(myid ==0) then
  write(24,*) 'urf values'
   do l=1,da(id)%ncg(3)
    do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
       write(24,*) da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)*ur(l,m,n),ur(l,m,n)
     end do
    end do
  end do
  write(24,*) da(id)%Q(1,1,1,1), da(id)%jacob(1,1,1),ur(1,1,1), purav(1,1,1), pav(1,1,1), urf(1,1,1), vrf(1,1,1), wrf(1,1,1),&
              & urf(3,3,3), vrf(3,3,3), wrf(3,3,3)
end if 

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
       LL(l,m,n,1)= ( pav(l,m,n)/count(l,m,n) )*(uuf(l,m,n)-uff(l,m,n)*uff(l,m,n))
       LL(l,m,n,2)= ( pav(l,m,n)/count(l,m,n) )*(uvf(l,m,n)-uff(l,m,n)*vff(l,m,n))
       LL(l,m,n,3)= ( pav(l,m,n)/count(l,m,n) )*(uwf(l,m,n)-uff(l,m,n)*wff(l,m,n))
       LL(l,m,n,4)= ( pav(l,m,n)/count(l,m,n) )*(vuf(l,m,n)-vff(l,m,n)*uff(l,m,n))
       LL(l,m,n,5)= ( pav(l,m,n)/count(l,m,n) )*(vvf(l,m,n)-vff(l,m,n)*vff(l,m,n))
       LL(l,m,n,6)= ( pav(l,m,n)/count(l,m,n) )*(vwf(l,m,n)-vff(l,m,n)*wff(l,m,n))
       LL(l,m,n,7)= ( pav(l,m,n)/count(l,m,n) )*(wuf(l,m,n)-wff(l,m,n)*uff(l,m,n))
       LL(l,m,n,8)= ( pav(l,m,n)/count(l,m,n) )*(wvf(l,m,n)-wff(l,m,n)*vff(l,m,n))
       LL(l,m,n,9)= ( pav(l,m,n)/count(l,m,n) )*(wwf(l,m,n)-wff(l,m,n)*wff(l,m,n))
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

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      rho = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)
      puurav(l,m,n)= rho*uf(l,m,n)*ur(l,m,n)
      puvrav(l,m,n)= rho*uf(l,m,n)*vr(l,m,n)
      puwrav(l,m,n)= rho*uf(l,m,n)*wr(l,m,n)
      pvurav(l,m,n)= rho*vf(l,m,n)*ur(l,m,n)
      pvvrav(l,m,n)= rho*vf(l,m,n)*vr(l,m,n)
      pvwrav(l,m,n)= rho*vf(l,m,n)*wr(l,m,n)
      pwurav(l,m,n)= rho*wf(l,m,n)*ur(l,m,n)
      pwvrav(l,m,n)= rho*wf(l,m,n)*vr(l,m,n)
      pwwrav(l,m,n)= rho*wf(l,m,n)*wr(l,m,n)
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
              puurav(l,m,n) = puurav(l,m,n)+rhor*uf(p,q,r)*ur(p,q,r)
              puvrav(l,m,n) = puvrav(l,m,n)+rhor*uf(p,q,r)*vr(p,q,r)
              puwrav(l,m,n) = puwrav(l,m,n)+rhor*uf(p,q,r)*wr(p,q,r)
              pvurav(l,m,n) = pvurav(l,m,n)+rhor*vf(p,q,r)*ur(p,q,r)
              pvvrav(l,m,n) = pvvrav(l,m,n)+rhor*vf(p,q,r)*vr(p,q,r)
              pvwrav(l,m,n) = pvwrav(l,m,n)+rhor*vf(p,q,r)*wr(p,q,r)
              pwurav(l,m,n) = pwurav(l,m,n)+rhor*wf(p,q,r)*ur(p,q,r)
              pwvrav(l,m,n) = pwvrav(l,m,n)+rhor*wf(p,q,r)*vr(p,q,r)
              pwwrav(l,m,n) = pwwrav(l,m,n)+rhor*wf(p,q,r)*wr(p,q,r)

            end if 
          end do
        end do
      end do
    uurf(l,m,n) = puurav(l,m,n)/pav(l,m,n)
    uvrf(l,m,n) = puvrav(l,m,n)/pav(l,m,n)
    uwrf(l,m,n) = puwrav(l,m,n)/pav(l,m,n)
    vurf(l,m,n) = pvurav(l,m,n)/pav(l,m,n)
    vvrf(l,m,n) = pvvrav(l,m,n)/pav(l,m,n)
    vwrf(l,m,n) = pvwrav(l,m,n)/pav(l,m,n)
    wurf(l,m,n) = pwurav(l,m,n)/pav(l,m,n)
    wvrf(l,m,n) = pwvrav(l,m,n)/pav(l,m,n)
    wwrf(l,m,n) = pwwrav(l,m,n)/pav(l,m,n)
     end do
   end do
end do

uruf(: , : , :) = uurf(:,:,:)
urvf(: , : , :) = vurf(:,:,:)
urwf(: , : , :) = uwrf(:,:,:)
vruf(: , : , :) = uvrf(:,:,:)
vrvf(: , : , :) = vvrf(:,:,:)
vrwf(: , : , :) = wvrf(:,:,:)
wruf(: , : , :) = uwrf(:,:,:)
wrvf(: , : , :) = vwrf(:,:,:)
wrwf(: , : , :) = wwrf(:,:,:)

if (myid==0) then
  write(24,*) 'uur values'
  write(24,*) '**********************'
  write(24,*) uf(1,1,1),   ur(1,1,1)
  write(24,*) uurf(2,2,2), uvrf(2,2,2), uwrf(2,2,2), vurf(2,2,2), vvrf(2,2,2), vwrf(2,2,2), wurf(2,2,2), wvrf(2,2,2), wwrf(2,2,2)
  write(24,*) uurf(3,3,3), uvrf(3,3,3), uwrf(3,3,3), vurf(3,3,3), vvrf(3,3,3), vwrf(3,3,3), wurf(3,3,3), wvrf(3,3,3), wwrf(3,3,3)
end if 

if (myid==0) then
  write(24,*) 'uru values'
  write(24,*) uruf(2,2,2), urvf(2,2,2), urwf(2,2,2), vruf(2,2,2), vrvf(2,2,2), vrwf(2,2,2), wruf(2,2,2), wrvf(2,2,2), wrwf(2,2,2)
  write(24,*) uruf(3,3,3), urvf(3,3,3), urwf(3,3,3), vruf(3,3,3), vrvf(3,3,3), vrwf(3,3,3), wruf(3,3,3), wrvf(3,3,3), wrwf(3,3,3)
end if 


 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
       MM(l,m,n,1)= ( pav(l,m,n)/count(l,m,n) )*(uurf(l,m,n)+uruf(l,m,n)-uff(l,m,n)*urf(l,m,n)-urf(l,m,n)*uff(l,m,n))
       MM(l,m,n,2)= ( pav(l,m,n)/count(l,m,n) )*(uvrf(l,m,n)+urvf(l,m,n)-uff(l,m,n)*vrf(l,m,n)-urf(l,m,n)*vff(l,m,n))
       MM(l,m,n,3)= ( pav(l,m,n)/count(l,m,n) )*(uwrf(l,m,n)+wurf(l,m,n)-uff(l,m,n)*wrf(l,m,n)-urf(l,m,n)*wff(l,m,n))
       MM(l,m,n,4)= ( pav(l,m,n)/count(l,m,n) )*(vurf(l,m,n)+uvrf(l,m,n)-vff(l,m,n)*urf(l,m,n)-vrf(l,m,n)*uff(l,m,n))
       MM(l,m,n,5)= ( pav(l,m,n)/count(l,m,n) )*(vvrf(l,m,n)+vrvf(l,m,n)-vff(l,m,n)*vrf(l,m,n)-vrf(l,m,n)*vff(l,m,n))
       MM(l,m,n,6)= ( pav(l,m,n)/count(l,m,n) )*(vwrf(l,m,n)+wvrf(l,m,n)-vff(l,m,n)*wrf(l,m,n)-vrf(l,m,n)*wff(l,m,n))
       MM(l,m,n,7)= ( pav(l,m,n)/count(l,m,n) )*(wurf(l,m,n)+uwrf(l,m,n)-urf(l,m,n)*wff(l,m,n)-wrf(l,m,n)*uff(l,m,n))
       MM(l,m,n,8)= ( pav(l,m,n)/count(l,m,n) )*(wvrf(l,m,n)+vwrf(l,m,n)-wff(l,m,n)*vrf(l,m,n)-wrf(l,m,n)*vff(l,m,n))
       MM(l,m,n,9)= ( pav(l,m,n)/count(l,m,n) )*(wwrf(l,m,n)+wrwf(l,m,n)-wff(l,m,n)*wrf(l,m,n)-wrf(l,m,n)*wff(l,m,n))
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
pururav(:,:,:) = 0.0d0
purvrav(:,:,:) = 0.0d0
purwrav(:,:,:) = 0.0d0
pvrurav(:,:,:) = 0.0d0
pvrvrav(:,:,:) = 0.0d0
pvrwrav(:,:,:) = 0.0d0
pwrurav(:,:,:) = 0.0d0 
pwrvrav(:,:,:) = 0.0d0
pwrwrav(:,:,:) = 0.0d0
 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
      rho = da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n)
      pururav(l,m,n)= rho*ur(l,m,n)*ur(l,m,n)
      purvrav(l,m,n)= rho*ur(l,m,n)*vr(l,m,n)
      purwrav(l,m,n)= rho*ur(l,m,n)*wr(l,m,n)
      pvrurav(l,m,n)= rho*vr(l,m,n)*ur(l,m,n)
      pvrvrav(l,m,n)= rho*vr(l,m,n)*vr(l,m,n)
      pvrwrav(l,m,n)= rho*vr(l,m,n)*wr(l,m,n)
      pwrurav(l,m,n)= rho*wr(l,m,n)*ur(l,m,n)
      pwrvrav(l,m,n)= rho*wr(l,m,n)*vr(l,m,n)
      pwrwrav(l,m,n)= rho*wr(l,m,n)*wr(l,m,n)
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
              pururav(l,m,n) = pururav(l,m,n)+rhor*ur(p,q,r)*ur(p,q,r)
              purvrav(l,m,n) = purvrav(l,m,n)+rhor*ur(p,q,r)*vr(p,q,r)
              purwrav(l,m,n) = purwrav(l,m,n)+rhor*ur(p,q,r)*wr(p,q,r)
              pvrurav(l,m,n) = pvrurav(l,m,n)+rhor*vr(p,q,r)*ur(p,q,r)
              pvrvrav(l,m,n) = pvrvrav(l,m,n)+rhor*vr(p,q,r)*vr(p,q,r)
              pvrwrav(l,m,n) = pvrwrav(l,m,n)+rhor*vr(p,q,r)*wr(p,q,r)
              pwrurav(l,m,n) = pwrurav(l,m,n)+rhor*wr(p,q,r)*ur(p,q,r)
              pwrvrav(l,m,n) = pwrvrav(l,m,n)+rhor*wr(p,q,r)*vr(p,q,r)
              pwrwrav(l,m,n) = pwrwrav(l,m,n)+rhor*wr(p,q,r)*wr(p,q,r)

            end if 
          end do
        end do
      end do
    ururf(l,m,n) = pururav(l,m,n)/pav(l,m,n)
    urvrf(l,m,n) = purvrav(l,m,n)/pav(l,m,n)
    urwrf(l,m,n) = purwrav(l,m,n)/pav(l,m,n)
    vrurf(l,m,n) = pvrurav(l,m,n)/pav(l,m,n)
    vrvrf(l,m,n) = pvrvrav(l,m,n)/pav(l,m,n)
    vrwrf(l,m,n) = pvrwrav(l,m,n)/pav(l,m,n)
    wrurf(l,m,n) = pwrurav(l,m,n)/pav(l,m,n)
    wrvrf(l,m,n) = pwrvrav(l,m,n)/pav(l,m,n)
    wrwrf(l,m,n) = pwrwrav(l,m,n)/pav(l,m,n)
    end do
  end do
end do


if (myid==0) then
  write(24,*) 'urur values'
  write(24,*) ururf(2,2,2), urvrf(2,2,2), urwrf(2,2,2), vrurf(2,2,2), vrvrf(2,2,2), vrwrf(2,2,2), wrurf(2,2,2), wrvrf(2,2,2),&
  &wrwrf(2,2,2)
  write(24,*) ururf(3,3,3), urvrf(3,3,3), urwrf(3,3,3), vrurf(3,3,3), vrvrf(3,3,3), vrwrf(3,3,3), wrurf(3,3,3), wrvrf(3,3,3),&
  &wrwrf(3,3,3)
end if 

 do l=1,da(id)%ncg(3)
   do m=1,da(id)%ncg(2)
     do n=1,da(id)%ncg(1)
       RR(l,m,n,1)= ( pav(l,m,n)/count(l,m,n) )*( (ururf(l,m,n)-( urf(l,m,n)*urf(l,m,n) ) ) )
       RR(l,m,n,2)= ( pav(l,m,n)/count(l,m,n) )*( (urvrf(l,m,n)-( urf(l,m,n)*vrf(l,m,n) ) ) )
       RR(l,m,n,3)= ( pav(l,m,n)/count(l,m,n) )*( (urwrf(l,m,n)-( urf(l,m,n)*wrf(l,m,n) ) ) )
       RR(l,m,n,4)= ( pav(l,m,n)/count(l,m,n) )*( (vrurf(l,m,n)-( vrf(l,m,n)*urf(l,m,n) ) ) )
       RR(l,m,n,5)= ( pav(l,m,n)/count(l,m,n) )*( (vrvrf(l,m,n)-( vrf(l,m,n)*vrf(l,m,n) ) ) )
       RR(l,m,n,6)= ( pav(l,m,n)/count(l,m,n) )*( (vrwrf(l,m,n)-( vrf(l,m,n)*wrf(l,m,n) ) ) )
       RR(l,m,n,7)= ( pav(l,m,n)/count(l,m,n) )*( (wrurf(l,m,n)-( wrf(l,m,n)*urf(l,m,n) ) ) )
       RR(l,m,n,8)= ( pav(l,m,n)/count(l,m,n) )*( (wrvrf(l,m,n)-( wrf(l,m,n)*vrf(l,m,n) ) ) )
       RR(l,m,n,9)= ( pav(l,m,n)/count(l,m,n) )*( (wrwrf(l,m,n)-( wrf(l,m,n)*wrf(l,m,n) ) ) )
     end do
   end do
 end do
 if (myid==0) then
  write(24,*) 'RR values'
  write(24,*) RR(2,2,2,1), RR(2,2,2,2), RR(2,2,2,3), RR(2,2,2,4), RR(2,2,2,5), RR(2,2,2,6), RR(2,2,2,7), RR(2,2,2,8), RR(2,2,2,9)
  write(24,*) RR(3,3,3,1), RR(3,3,3,2), RR(3,3,3,3), RR(3,3,3,4), RR(3,3,3,5), RR(3,3,3,6), RR(3,3,3,7), RR(3,3,3,8), RR(3,3,3,9)
end if 
 T = LL+MM+RR
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
   if (myid==0) then
    write(*,*) 'smagvcalculation has been called'
   end if
   do id=1,ngridla
    if (myid==0) then
      write(*,*) '22smagvcalculation has been called'
    end if
      do nv=1,neq
        call Gauss_Deriv(da(id)%gmetg,da(id)%Qlgg(:,:,:,nv),da(id)%Qglg(:,:,:,nv),&
                         da(id)%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv),&
                         da(id)%ncg,da(id)%dmx,da(id)%dmy,da(id)%dmz)
      end do
    if (myid==0) then
      write(*,*) '32smagvcalculation has been called'
    end if
      pav = 0
   do l=1,da(id)%ncg(1)
      do m=1,da(id)%ncg(2)
        do n=1,da(id)%ncg(3)
          pav = pav + ( da(id)%Q(l,m,n,1)*da(id)%jacob(l,m,n) )
        end do
      end do
   end do
    
  if (myid==0) then
    write(25,*) pav
  end if 


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
    if (myid==0) then
      write(25,*) 'sum1'
      write(25,*) (sum(2,i), i=1,9)
      write(25,*) 'sum2'
      write(25,*) (sum2(i), i=1,9)
    end if
   call MPI_Allreduce(ngridla, ngridaa , 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm1d, ierr)
   sum2 = sum2/ngridaa
   if (myid==0) then
    write(25,*)  'averaged sum2'
    write(25,11) (sum2(i), i=1,9)
   end if
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