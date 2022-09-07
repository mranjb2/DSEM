!
!///////////////////////////////////////////////////////////////////////
!
      MODULE user_keywords
!
!......................................................................
!     store user modifiable keywords needed to read the input file.
!     any number of user keywords can be added to be accessible from
!     the user modifiable routines:
!
!        user_setup
!        do_user_stuff
!        cleanup
!        initc
!        + other user supplied routines called by them
!
!     communication is through the data defined in the module "User_Data"
!     The values in User_Data are set below in the routine "set_user_values"
!......................................................................
!
      SAVE
!
      INTEGER                           :: num_user_keywords = 3
      CHARACTER (LEN=10), DIMENSION(3)  :: user_keyword_list = (/ &
                                                'theta     '    , &
                                                'phi       '    ,  &
                                                'outlet    '      &
                                           /)
      LOGICAL, DIMENSION(3) :: user_keyword_set = .false.
      END MODULE user_keywords
!
!
!///////////////////////////////////////////////////////////////////////
!
MODULE User_Data
   USE size
   SAVE

   DOUBLE PRECISION          :: pout         ! outlet pressure
   DOUBLE PRECISION, DIMENSION(ngp,nx,ny,nz,neq) :: r,s,t        ! viscous fluxes
   DOUBLE PRECISION :: xxin(38),uuin(38),xx1in(94),uu1in(94)
   DOUBLE PRECISION :: xx2in(61),uu2in(61),xx3in(61),uu3in(61)
   INTEGER          :: itimeinflowbc
   DOUBLE PRECISION :: intensity = 0.025d0,prestep 
   DOUBLE PRECISION :: ufluct(nx-1,nx-1,2),vfluct(nx-1,nx-1,2),wfluct(nx-1,nx-1,2)
   DOUBLE PRECISION :: urms(nx-1,nx-1),vrms(nx-1,nx-1),wrms(nx-1,nx-1)
   DOUBLE PRECISION :: umean(nx-1,nx-1),vmean(nx-1,nx-1),wmean(nx-1,nx-1)
   INTEGER, PARAMETER :: Nfilt=2
   INTEGER, PARAMETER :: nc_new=5
   DOUBLE PRECISION :: bz_filt(Nfilt,5), bz_proj(5,Nfilt)
   DOUBLE PRECISION :: bxltog(nx-1,nx),bxgtol(nx,nx-1)
   DOUBLE PRECISION,ALLOCATABLE :: diss(:,:,:,:,:)
   DOUBLE PRECISION,ALLOCATABLE :: subdiss(:,:,:,:) 
   DOUBLE PRECISION,ALLOCATABLE :: ima(:,:,:)   !minimum of six domains!!!!!!!!!!!!!!!
   DOUBLE PRECISION,ALLOCATABLE :: imb(:,:,:)  
   DOUBLE PRECISION,ALLOCATABLE :: imc(:,:,:)                          !!!!!!!!!!!!!!!
   DOUBLE PRECISION,ALLOCATABLE :: btohigh(:,:)!(N_high,nx)
   DOUBLE PRECISION,ALLOCATABLE :: btolow(:,:)
   DOUBLE PRECISION,ALLOCATABLE :: bgtol_high(:,:)                     !!!!!!!!!!!!!!!


END MODULE User_Data
!
!///////////////////////////////////////////////////////////////////////
!
      MODULE Skin_Data
         USE size
         SAVE
         DOUBLE PRECISION, ALLOCATABLE          :: bmatu1b(:,:,:)      ! integration matrix
         DOUBLE PRECISION, ALLOCATABLE          :: bmataa(:,:)      ! integration matrix

      END MODULE Skin_Data
!
!///////////////////////////////////////////////////////////////////////
module particleInlets
    use size
    integer                          :: inletCount
    integer, allocatable             :: inletList(:)
    double precision, allocatable    :: cumPart(:)
    
end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE set_user_values(the_word,input_line,kword)
!
!......................................................................
!     read in the input data-file and save the values in User_Data
!......................................................................
!
      USE physics
      USE User_Data
      USE user_keywords
      USE constants
!
      DOUBLE PRECISION      :: get_dp_value   ! extracts a double prec number
      INTEGER               :: get_int_value  ! extracts an integer
      CHARACTER (LEN = 132) :: input_line
      CHARACTER (LEN = 10)  :: the_word
!
!     ---------------------------------
!     two user keywords for this problem
!     ---------------------------------
!
      SELECT CASE(the_word)
         CASE('theta')
            theta = get_dp_value(input_line)
            theta = theta*180.d0/pi
            user_keyword_set(kword) = .true.
         CASE('phi')
            phi = get_dp_value(input_line)
            phi = phi*180.d0/pi
            user_keyword_set(kword) = .true.
         CASE('outlet')
            pout = get_dp_value(input_line)
            user_keyword_set(kword) = .true.


      END SELECT
!
      RETURN
      END SUBROUTINE SET_USER_VALUES
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE user_setup(d,ngrid,mrtr,nmort,restart)
!
!......................................................................
!     ALLOWS PROBLEM SPECIFIC COMPUTATIONS TO BE DONE
!......................................................................
      USE size      !!!!!!!!!!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE Skin_Data
      USE mpi_par   !!!!!!!!!!
      USE mpi       !!!!!!!!!!
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
      TYPE (mortar) :: mrtr(nmp)

!
      DOUBLE PRECISION                  :: xwgt(nmax),cxg(nmax)
      INTEGER:: o, fh=126252
      INTEGER,DIMENSION(ngp)            :: npoin
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
      DOUBLE PRECISION                  :: polyn
      INTEGER, PARAMETER                :: legendre = 1
      INTEGER,DIMENSION(nmort)          :: mortmat
      LOGICAL                           :: nodouble,restart
      integer                           :: yperiodic, xperiodic, zperiodic

      ALLOCATE(bmatu1b(nmax,nmax,nmax))   !!!!!!!!!
      ALLOCATE(bmataa(nmax-1,nmax-1))
      ALLOCATE(btohigh(nxo-1,nx-1))
      ALLOCATE(btolow(nx-1,nxo-1))
      ALLOCATE(bgtol_high(Nh+1,nx))       !!!!!!!!!


IF (.not. restart) THEN
      mortmat = 0
      chick = 1
      IF (chick ==1) THEN
!
!  set periodic boundary conditions
!

!  find the interfaces  at which periodic boundary conditions
!  apply, adjust domain and mortar parameters
      zperiodic = 1
      if(zperiodic == 1) then
    !write(*,*)'hi jon line 142'
            kkk = 0
            do idm = 1,ngrid
                  if(d(idm)%bcond(3) == 'outflow') then
                        d(idm)%bcond(3)  = 'interface'
                        d(idm)%ibtype(3) = 1
                        mrtr_no          = d(idm)%mortar(1,3)
                        d(idm)%which_surface(3) = 0
                        do ids=1,ngrid
                              if (d(ids)%bcond(5)    =='outflow'            .and.  &
                                  d(idm)%corner(1,1) == d(ids)%corner(1,5)  .and.  &
                                  d(idm)%corner(2,1) == d(ids)%corner(2,5))  then
                        !write(*,*)'hello jon line 154'
                        kkk          = kkk + 1
                        mrtr_no_sl   = d(ids)%mortar(1,5)
                        mortmat(kkk) = mrtr_no_sl 
                        
                        mrtr(mrtr_no)%len(:,2)  = mrtr(mrtr_no_sl)%len(:,1)
                        mrtr(mrtr_no)%id(2)     = ids
                        mrtr(mrtr_no)%iface(2)  = 5
                        d(ids)%bcond(5)         = 'interface'
                        d(ids)%ibtype(5)        = 1
                        d(ids)%which_surface(5) = 0
                        d(ids)%mortar(1,5)      = mrtr_no
                        d(ids)%mortar(2,5)      = 2
                    end if
                end do
            end if
        end do
    end if


      yperiodic = 1
      if(yperiodic == 1) then
      !This is periodic in y
      !kkk = 0 !a counter
      DO idm = 1,ngrid
         IF (d(idm)%bcond(1) == 'outflow') THEN  !if inflow on face 1
            d(idm)%bcond(1)  = 'interface'      !Change to interface
            d(idm)%ibtype(1) = 1
            mrtr_no          = d(idm)%mortar(1,1) !Finds the corresponding mortar
            d(idm)%which_surface(1) = 0 
            DO ids = 1,ngrid                      !Over all the elements
               IF (d(ids)%bcond(2) == 'outflow'              .and.  & !is the face on the other side of the element inflow on face 2?
                   d(idm)%corner(1,1) == d(ids)%corner(1,4) .and.  & !do corner points match?
                   d(idm)%corner(3,1) == d(ids)%corner(3,4)) THEN    !matches the x & z -corner 1 is x and corner 3 is z
                                                                     !if 1 and 4 match then same element

                   kkk = kkk + 1
                   mrtr_no_sl   = d(ids)%mortar(1,2)
                   mortmat(kkk) = mrtr_no_sl
                  
                   mrtr(mrtr_no)%len(:,2) = mrtr(mrtr_no_sl)%len(:,1)
                   mrtr(mrtr_no)%id(2)    = ids 
                   mrtr(mrtr_no)%iface(2) = 2
                   d(ids)%bcond(2)  = 'interface'
                   d(ids)%ibtype(2) = 1
                   d(ids)%which_surface(2) = 0
                   d(ids)%mortar(1,2) = mrtr_no
                   d(ids)%mortar(2,2) = 2
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      end if
!
!  find the interfaces  at which periodic boundary conditions
!  apply, adjust domain and mortar parameters
!
    chicka=1 !periodic in x
    IF (chicka==1) THEN
      DO idm = 1,ngrid
         IF (d(idm)%bcond(6) == 'outflow') THEN
             d(idm)%bcond(6)  = 'interface'
             d(idm)%ibtype(6) = 1
             mrtr_no          = d(idm)%mortar(1,6)
             d(idm)%which_surface(6) = 0

            DO ids = 1,ngrid
               IF (d(ids)%bcond(4) == 'outflow'              .and. &
                   d(idm)%corner(3,1) == d(ids)%corner(3,2) .and.  &
                   d(idm)%corner(2,1) == d(ids)%corner(2,2)) THEN

                   kkk = kkk + 1
                   mrtr_no_sl   = d(ids)%mortar(1,4)
                   mortmat(kkk) = mrtr_no

                   mrtr(mrtr_no)%len(:,2) = mrtr(mrtr_no_sl)%len(:,1)
                   mrtr(mrtr_no)%id(2)    = ids
                   mrtr(mrtr_no)%iface(2) = 4
                   mrtr(mrtr_no_sl)%len(:,2) = mrtr(mrtr_no)%len(:,1)
                   mrtr(mrtr_no_sl)%id(2)    = idm
                   mrtr(mrtr_no_sl)%iface(2) = 6

                   d(ids)%bcond(4)  = 'interface'
                   d(ids)%ibtype(4) = 1
                   d(ids)%which_surface(4) = 0
                   d(ids)%mortar(1,4) = mrtr_no
                   d(ids)%mortar(2,4) = 1
              ENDIF
            ENDDO
         ENDIF
      ENDDO
    ENDIF


      nmortdouble = kkk
      ENDIF
!
!  Clean up mortars
!
      nmorta = 0
      DO im = 1,nmort
        nodouble = .true.
        DO imd = 1,kkk !nmortdouble
           IF (im==mortmat(imd)) THEN
             nodouble = .false.
           ENDIF
        ENDDO
        IF (nodouble) THEN 
          nmorta = nmorta +1
          mrtr(nmorta) = mrtr(im)
          idm          = mrtr(im)%id(1) 
          ids          = mrtr(im)%id(2)
          ifacem       = mrtr(im)%iface(1) 
          ifaces       = mrtr(im)%iface(2) 
          d(idm)%mortar(1,ifacem) = nmorta 
          IF (ids /= 0) THEN
            d(ids)%mortar(1,ifaces) = nmorta 
          ENDIF
        ENDIF
 
      ENDDO
      write(*,*) 'chick nmort old and new due to periodic boundary conditions',nmort,nmorta
      nmort = nmorta
ENDIF
!      ama      = 0.4
!      twall    = 1/ama**2
      itimeinflowbc=0
!
!  Initialize matrices
!
      bmatu1b = 0.0d0
      bmataa  = 0.0d0
      r       = 0.0d0
      s       = 0.0d0
      t       = 0.0d0
!
!  Compute integration matrices
!
      cxg  = 0.0d0
      ncg  = d(1)%ncg(1)
      CALL gaussq(legendre,ncg,0.0d0,0.0d0,0,endpts,work,cxg,xwgt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DO i = 1,ncg
      !    DO j = 1,ncg
      !      DO l = 1,ncg
      !       sum = 0.0d0
      !       DO k = 1,ncg
      !          hl = polyn(l,0.5d0*(cxg(k) + 1.d0),ncg,d(1)%cxg(1:ncg,1))
      !          hj = polyn(j,0.5d0*(cxg(k) + 1.d0),ncg,d(1)%cxg(1:ncg,1))
      !          hi = polyn(i,0.5d0*(cxg(k) + 1.d0),ncg,d(1)%cxg(1:ncg,1))
      !          sum = sum + 0.5d0*xwgt(k)*hi*hj*hl
      !       END DO
      !       bmatu1b(i,j,l) = sum
      !      END DO
      !    END DO
      !END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO i = 1,ncg
         DO j = 1,ncg
           DO l = 1,ncg
            sum = 0.0d00
            DO k = 1,ncg
              DO m = 1,ncg
                DO n= 1,ncg
                   hl = polyn(l,0.5d0*(cxg(n) + 1.0d0),ncg,d(1)%cxg(1:ncg,1))
                   hj = polyn(j,0.5d0*(cxg(m) + 1.0d0),ncg,d(1)%cxg(1:ncg,1))
                   hi = polyn(i,0.5d0*(cxg(k) + 1.0d0),ncg,d(1)%cxg(1:ncg,1))
                   sum = sum + 0.1250d00*xwgt(k)*xwgt(m)*xwgt(n)*hi*hj*hl
                END DO
              ENDDO
            ENDDO
            bmatu1b(i,j,l) = sum
           END DO
         END DO
      END DO

      CALL compute_bmata(bmataa,d(1)%cxg(1:ncg,1),ncg,ncg)
      bgtol_high=0.0d00000
      ALLOCATE(ima((nfour/4+1)**3,ngrid,nx))
      ALLOCATE(imb((nfour/4+1)**3,ngrid,nx))
      ALLOCATE(imc((nfour/4+1)**3,ngrid,nx))
!

      RETURN
      END SUBROUTINE user_setup
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE do_user_stuff(time,dt,ngrid,d,mrtr,nmort,stats,stat_unit,kk)
!
!......................................................................
    ! allow the user to do periodic computations/io,
    ! such as a probe measurement. Frequency of calling is set by the
    ! variable "nout" i.e. by the keyword "frequency"
!......................................................................
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       USE domain_definition
!       USE mortar_definition
!       USE physics                                                                                         
!       USE constants                                                               
!       USE User_Data
!       USE Skin_Data
!       USE mpi_par
!       use mpi
!       USE size
!       USE particle_definition
!       USE part_par
      
! !
!       IMPLICIT DOUBLE PRECISION (a-h,o-z)
!       TYPE (domain)   :: d(ngp)
!       TYPE (mortar)   :: mrtr(nmp)
!       type (particle) :: drop(1)
! !
!       INTEGER :: stat_unit
! !
!     Write stuff to file with extension .sta

!      IF (myid ==0) THEN
!        u = d(1)%Q(3,3,3,2)
!        v = d(1)%Q(3,3,3,3)
!        w = d(1)%Q(3,3,3,4)
!        rho = d(1)%Q(3,3,3,1)
!        p = d(1)%Q(3,3,3,5)*(gamma-1)-0.5*rho*(u*u+v*v+w*w)
!        write(stat_unit,*)  time, p, itimeinflowbc
!      END IF
!      RETURN

!   do i=1,ngrid
!       ncg = d(i)%ncg(1)
!       write(*,*)'n:',i
!       write(*,*)'max:',maxval(d(i)%Q(1:ncg,1:ncg,1:ncg,2)),'min:',minval(d(i)%Q(1:ncg,1:ncg,1:ncg,2))
!   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      USE size
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE Skin_Data
      USE mpi_par
      USE mpi
      ! USE stats_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngrid)
      TYPE (mortar) :: mrtr(nmort)
      ! TYPE (statfl) :: stats(ngp)

      DOUBLE PRECISION                  :: xwgt(nmax),cxg(nmax)
      DOUBLE PRECISION                  :: bmat(nmax,nmax,nmax)
      INTEGER, PARAMETER                :: legendre = 1
      INTEGER                           :: o,fh=31231
      INTEGER,DIMENSION(ngp)            :: npoin
     
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
      DOUBLE PRECISION, DIMENSION(10)   :: statall,statalla
      DOUBLE PRECISION                  :: polyn,polyn_ami
      DOUBLE PRECISION                  :: som,som1,som2,dt
      CHARACTER(LEN=32)                 :: fname,vname
      COMPLEX, DIMENSION(nfour,nfour,nfour) :: ufour,vfour,wfour,u_test,v_test,w_test,udum
!
      INTEGER :: stat_unit,unt,kk
!
!     Write stuff to file with extension .sta
!
      ! Allocate dynamic arrays
      douser = 1

    if(kk==1) then
        prestep = 0.0d0
    end if

      IF( douser ==1 .and. MOD( kk, 1 ) == 0) THEN


      ALLOCATE(diss(ngrid,nx,ny,nz,3))
      diss=0.0d00
      ! subdiss=0.0d00
      ! ima=0.0d0
      ! imb=0.0d0
      ! imc=0.0d0
      statall=0.0d000
      statalla=0.0d000

      
!     ------------------------------
!     set up gauss grids and weights
!     ------------------------------
!     
     ! CALL linetrack(d,id,time)
     
      bmat = 0.0d00
      cxg  = 0.0d00
      ncg  = d(1)%ncg(1)
      CALL gaussq(legendre,ncg,0.0d00,0.0d00,0,endpts,work,cxg,xwgt)

      bmat=bmatu1b

      DO id=1,ngridl
         CALL compute_diss(d(id),id)
      ENDDO
!
!    ----------------------
!    compute averages
!    ----------------------
!
      uav = 0.0d000
      vav = 0.0d000
      wav = 0.0d00
      eav = 0.0d00
      rav = 0.0d00 
      tav = 0.0d00
      epav= 0.0d00
      ediss1= 0.0d0
      ediss2= 0.0d0
      sdiss = 0.0d0 
   !
      DO id =1, ngridl
         DO i = 1,ncg
           DO j = 1,ncg
             DO l = 1,ncg
               u    =  d(id)%Q(i,j,l,2)/d(id)%Q(i,j,l,1)
               v    =  d(id)%Q(i,j,l,3)/d(id)%Q(i,j,l,1)
               w    =  d(id)%Q(i,j,l,4)/d(id)%Q(i,j,l,1)
               rho  =  d(id)%Q(i,j,l,1)*d(id)%jacob(i,j,l)
               temp = d(id)%Q(i,j,l,5)/d(id)%Q(i,j,l,1)
               temp = temp - 0.50*(u**2+v**2+w**2)
               temp = temp*gamma*(gamma-1.0d00)
               press = temp*rho/gamma

               uav = uav + u*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
               vav = vav + v*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
               wav = wav + w*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
               eav = eav + d(id)%Q(i,j,l,5)*bmat(i,j,l)*(d(id)%jacob(i,j,l))
               rav = rav + rho*bmat(i,j,l)
               tav = tav + temp*bmat(i,j,l)
               epav= epav + press/(gamma-1.0d00)*bmat(i,j,l)
               ediss1 = ediss1 + diss(id,i,j,l,1)*bmat(i,j,l)
               ediss2 = ediss2 + diss(id,i,j,l,2)*bmat(i,j,l)
               ! if (id==1 .and. myid==0 .and. abs(diss(id,i,j,l,1))>1000000.0d0) write(*,*)'diss=',diss(id,i,j,l,1),i,j,l
               ! if (id==1 .and. myid==0 .and. abs(diss(id,i,j,l,1))>1000000.0d0) write(*,*)'Qlgg=',d(id)%Qlgg(i,j,l,2),i,j,l
               sdiss = sdiss + diss(id,i,j,l,3)*bmat(i,j,l)
             ENDDO
           ENDDO
         ENDDO
      ENDDO

   !
      uav = uav/d(1)%jacob(1,1,1)
      vav = vav/d(1)%jacob(1,1,1)
      wav = wav/d(1)%jacob(1,1,1)
      eav = eav/d(1)%jacob(1,1,1)
      rav = rav/d(1)%jacob(1,1,1)
      tav = tav/d(1)%jacob(1,1,1)
      epav= epav/d(1)%jacob(1,1,1)
      ediss1= ediss1/d(1)%jacob(1,1,1)
      ediss2= ediss2/d(1)%jacob(1,1,1)
      sdiss= sdiss/d(1)%jacob(1,1,1)

      ! write(*,*)'uave=',uav,vav,wav,eav,rav,tav,epav,ediss1,ediss2,sdiss
      
   !
      statall(1) = uav
      statall(2) = vav
      statall(3) = wav
      statall(4) = eav
      statall(5) = rav
      statall(6) = tav
      statall(7) = epav
      statall(8) = ediss1
      statall(9) = ediss2
      statall(10) = sdiss
     
     ! write(*,*)statall,myid

      ! write(*,*)myid,'=',statall
      ! CALL MPI_BARRIER(comm1d, ierr)
      ! write(*,*)'DOUBLE PRECISION=',MPI_DOUBLE_PRECISION,MPI_SUM
   !
      CALL MPI_ALLREDUCE(statall,statalla,10,MPI_DOUBLE_PRECISION, MPI_SUM,comm1d,ierr)

        ! write(*,*)'521.',myid


      IF (myid==0) THEN
        uav = statalla(1)
        vav = statalla(2)
        wav = statalla(3)
        eav = statalla(4)
        rav = statalla(5)
        tav = statalla(6)
        epav= statalla(7)
        ediss1= statalla(8)
        ediss2= statalla(9)
        sdiss= statalla(10)

        uav = uav/(8.0d0*pi**3)
        vav = vav/(8.0d0*pi**3)
        wav = wav/(8.0d0*pi**3)
        eav = eav/(8.0d0*pi**3)
        rav = rav/(8.0d0*pi**3)
        tav = tav/(8.0d0*pi**3)
        epav= epav/(8.0d0*pi**3)
        ediss1= ediss1/(8.0d0*pi**3)
        ediss2= ediss2/(8.0d0*pi**3)
        sdiss= sdiss/(8.0d0*pi**3)

        statall(1) = uav
        statall(2) = vav
        statall(3) = wav
        statall(4) = eav
        statall(5) = rav
        statall(6) = tav
        statall(7) = epav
        statall(8) = ediss1
        statall(9) = ediss2
        statall(10) = sdiss
     ENDIF
     ! write(*,*)'542.',myid
!
     CALL MPI_BCAST(statall,10,MPI_DOUBLE_PRECISION,0,comm1d,ierr)

     !    -------------------------
!    compute root mean squares
!    -------------------------
!
      rav = statall(5)
      tav = statall(6)
      epav= statall(7)
   !
      urmss = 0.0d0
      vrmss = 0.0d0
      wrmss = 0.0d0
      rrms = 0.0d0
      trms = 0.0d0
      amarms= 0.0d0
      eprms= 0.0d0
   !
      DO id =1, ngrid
         DO i = 1,ncg
           DO j = 1,ncg
             DO l = 1,ncg
               u    =  d(id)%Q(i,j,l,2)/d(id)%Q(i,j,l,1)
               v    =  d(id)%Q(i,j,l,3)/d(id)%Q(i,j,l,1)
               w    =  d(id)%Q(i,j,l,4)/d(id)%Q(i,j,l,1)
               rho  =  d(id)%Q(i,j,l,1)*d(id)%jacob(i,j,l)
               temp = d(id)%Q(i,j,l,5)/d(id)%Q(i,j,l,1)
               temp = temp - 0.5d0*(u**2+v**2+w**2)
               temp = temp*gamma*(gamma-1.0d0)
               press= temp*rho/gamma

               rrms = rrms + ((rho-rav)**2)*bmat(i,j,l)
               trms = trms + ((temp-tav)**2)*bmat(i,j,l)
               urmss = urmss + (u**2)*bmat(i,j,l)
               vrmss = vrmss + (v**2)*bmat(i,j,l)
               wrmss = wrmss + (w**2)*bmat(i,j,l)
               eprms= eprms + ((press/(gamma-1.0d0)-epav)**2)*bmat(i,j,l)
               amarms= amarms + (u**2+v**2+w**2)*bmat(i,j,l)/tav
             ENDDO
           ENDDO
         ENDDO
      ENDDO
   !
      rrms = rrms/d(1)%jacob(1,1,1)
      trms = trms/d(1)%jacob(1,1,1)
      urmss = urmss/d(1)%jacob(1,1,1)
      vrmss = vrmss/d(1)%jacob(1,1,1)
      wrmss = wrmss/d(1)%jacob(1,1,1)
      eprms = eprms/d(1)%jacob(1,1,1)
      amarms = amarms/d(1)%jacob(1,1,1)
   !
      statall(1) = urmss
      statall(2) = vrmss
      statall(3) = wrmss
      statall(4) = rrms
      statall(5) = trms
      statall(6) = eprms
      statall(7) = amarms
   !
      CALL MPI_ALLREDUCE(statall,statalla,10,MPI_DOUBLE_PRECISION, MPI_SUM, &
                             comm1d, ierr)
   !
      IF (myid==0) THEN
        urmss = statalla(1)
        vrmss = statalla(2)
        wrmss = statalla(3)
        rrms = statalla(4)
        trms = statalla(5)
        eprms = statalla(6)
        amarms = statalla(7)

        turbk= 0.5d0*(urmss+vrmss+wrmss)/(8.0d0*pi**3)
        urmss = sqrt(urmss/(8.0d0*pi**3))
        vrmss = sqrt(vrmss/(8.0d0*pi**3))
        wrmss = sqrt(wrmss/(8.0d0*pi**3))
        rrms  = sqrt(rrms/(8.0d0*pi**3))
        trms  = sqrt(trms/(8.0d0*pi**3))
        eprms = sqrt(eprms/(8.0d0*pi**3))
        amarms = sqrt(amarms/(8.0d0*pi**3))
        ! write(*,*) prestep
        ! write(*,*) turbk
        dissp = (prestep-turbk)/(dt)

        prestep = turbk
        ! write(*,*) prestep
        ! no=INT(time*100.0)
        ! WRITE(fname, fmt='(a4,i3,a3)') 'average',no,'.23'
        ! OPEN(unit=21,file=fname)
        ! WRITE(vname, fmt='(a4,i3,a3)') 'rms',no,'.23'
        ! OPEN(unit=22,file=vname)
        ! WRITE(21,*) time,uav,vav,wav,tav,rav,eav,epav,ediss1,ediss2,sdiss
        write(22,*) time,turbk
        write(23,*) time,dissp
        ! WRITE(22,*) time,urmss,vrmss,wrmss,turbk,trms,rrms,eprms,amarms
        ! CLOSE(21)
        ! CLOSE(22)
      ENDIF
      DEALLOCATE(diss)
    ENDIF

       nspec=0
IF (nspec==1 .and. MOD(kk,45350)==0) THEN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(diss(ngrid,nx,ny,nz,3))

      diss=0.0d00
      statall= 0.0d000
      statalla=0.0d000

!     ------------------------------
!     set up gauss grids and weights
!     ------------------------------
!     
      bmat = 0.0d00
      cxg  = 0.0d00
      ncg  = d(1)%ncg(1)
      CALL gaussq(legendre,ncg,0.0d00,0.0d00,0,endpts,work,cxg,xwgt)

      bmat=bmatu1b
      DO id=1,ngridl
         CALL compute_diss(d(id),id)
      ENDDO


      

      ! write(*,*)'435.',myid
!

!    ----------------------
!    compute averages
!    ----------------------
!
      uav = 0.0d000
      vav = 0.0d000
      wav = 0.0d00
      eav = 0.0d00
      rav = 0.0d00 
      tav = 0.0d00
      epav= 0.0d00
      ediss1= 0.0d0
      ediss2= 0.0d0
      sdiss = 0.0d0 
   !
      DO id =1, ngridl
         DO i = 1,ncg
           DO j = 1,ncg
             DO l = 1,ncg
               u    =  d(id)%Q(i,j,l,2)/d(id)%Q(i,j,l,1)
               v    =  d(id)%Q(i,j,l,3)/d(id)%Q(i,j,l,1)
               w    =  d(id)%Q(i,j,l,4)/d(id)%Q(i,j,l,1)
               rho  =  d(id)%Q(i,j,l,1)*d(id)%jacob(i,j,l)
               temp = d(id)%Q(i,j,l,5)/d(id)%Q(i,j,l,1)
               temp = temp - 0.50*(u**2+v**2+w**2)
               temp = temp*gamma*(gamma-1.0d00)
               press = temp*rho/gamma

               uav = uav + u*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
               vav = vav + v*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
               wav = wav + w*bmat(i,j,l)/(d(id)%Q(i,j,l,1))
               eav = eav + d(id)%Q(i,j,l,5)*bmat(i,j,l)*(d(id)%jacob(i,j,l))
               rav = rav + rho*bmat(i,j,l)
               tav = tav + temp*bmat(i,j,l)
               epav= epav + press/(gamma-1.0d00)*bmat(i,j,l)
               ediss1 = ediss1 + diss(id,i,j,l,1)*bmat(i,j,l)
               ediss2 = ediss2 + diss(id,i,j,l,2)*bmat(i,j,l)
               ! if (id==1 .and. myid==0 .and. abs(diss(id,i,j,l,1))>1000000.0d0) write(*,*)'diss=',diss(id,i,j,l,1),i,j,l
               ! if (id==1 .and. myid==0 .and. abs(diss(id,i,j,l,1))>1000000.0d0) write(*,*)'Qlgg=',d(id)%Qlgg(i,j,l,2),i,j,l
               sdiss = sdiss + diss(id,i,j,l,3)*bmat(i,j,l)
             ENDDO
           ENDDO
         ENDDO
      ENDDO

   !
      uav = uav/d(1)%jacob(1,1,1)
      vav = vav/d(1)%jacob(1,1,1)
      wav = wav/d(1)%jacob(1,1,1)
      eav = eav/d(1)%jacob(1,1,1)
      rav = rav/d(1)%jacob(1,1,1)
      tav = tav/d(1)%jacob(1,1,1)
      epav= epav/d(1)%jacob(1,1,1)
      ediss1= ediss1/d(1)%jacob(1,1,1)
      ediss2= ediss2/d(1)%jacob(1,1,1)
      sdiss= sdiss/d(1)%jacob(1,1,1)

      ! write(*,*)'uave=',uav,vav,wav,eav,rav,tav,epav,ediss1,ediss2,sdiss
      
   !
      statall(1) = uav
      statall(2) = vav
      statall(3) = wav
      statall(4) = eav
      statall(5) = rav
      statall(6) = tav
      statall(7) = epav
      statall(8) = ediss1
      statall(9) = ediss2
      statall(10) = sdiss
   !
      CALL MPI_ALLREDUCE(statall,statalla,10,MPI_DOUBLE_PRECISION, MPI_SUM,comm1d,ierr)

        ! write(*,*)'521.',myid


      IF (myid==0) THEN
        uav = statalla(1)
        vav = statalla(2)
        wav = statalla(3)
        eav = statalla(4)
        rav = statalla(5)
        tav = statalla(6)
        epav= statalla(7)
        ediss1= statalla(8)
        ediss2= statalla(9)
        sdiss= statalla(10)

        uav = uav/(8.0d0*pi**3)
        vav = vav/(8.0d0*pi**3)
        wav = wav/(8.0d0*pi**3)
        eav = eav/(8.0d0*pi**3)
        rav = rav/(8.0d0*pi**3)
        tav = tav/(8.0d0*pi**3)
        epav= epav/(8.0d0*pi**3)
        ediss1= ediss1/(8.0d0*pi**3)
        ediss2= ediss2/(8.0d0*pi**3)
        sdiss= sdiss/(8.0d0*pi**3)

        statall(1) = uav
        statall(2) = vav
        statall(3) = wav
        statall(4) = eav
        statall(5) = rav
        statall(6) = tav
        statall(7) = epav
        statall(8) = ediss1
        statall(9) = ediss2
        statall(10) = sdiss
     ENDIF
     ! write(*,*)'542.',myid
!
     CALL MPI_BCAST(statall,10,MPI_DOUBLE_PRECISION,0,comm1d,ierr)

     !    -------------------------
!    compute root mean squares
!    -------------------------
!
      rav = statall(5)
      tav = statall(6)
      epav= statall(7)
   !
      urmss = 0.0d0
      vrmss = 0.0d0
      wrmss = 0.0d0
      rrms = 0.0d0
      trms = 0.0d0
      amarms= 0.0d0
      eprms= 0.0d0
   !
      DO id =1, ngrid
         DO i = 1,ncg
           DO j = 1,ncg
             DO l = 1,ncg
               u    =  d(id)%Q(i,j,l,2)/d(id)%Q(i,j,l,1)
               v    =  d(id)%Q(i,j,l,3)/d(id)%Q(i,j,l,1)
               w    =  d(id)%Q(i,j,l,4)/d(id)%Q(i,j,l,1)
               rho  =  d(id)%Q(i,j,l,1)*d(id)%jacob(i,j,l)
               temp = d(id)%Q(i,j,l,5)/d(id)%Q(i,j,l,1)
               temp = temp - 0.5d0*(u**2+v**2+w**2)
               temp = temp*gamma*(gamma-1.0d0)
               press= temp*rho/gamma
               rrms = rrms + ((rho-rav)**2)*bmat(i,j,l)
               trms = trms + ((temp-tav)**2)*bmat(i,j,l)
               urmss = urmss + (u**2)*bmat(i,j,l)
               vrmss = vrmss + (v**2)*bmat(i,j,l)
               wrmss = wrmss + (w**2)*bmat(i,j,l)
               eprms= eprms + ((press/(gamma-1.0d0)-epav)**2)*bmat(i,j,l)
               amarms= amarms + (u**2+v**2+w**2)*bmat(i,j,l)/tav
             ENDDO
           ENDDO
         ENDDO
      ENDDO
   !
      rrms = rrms/d(1)%jacob(1,1,1)
      trms = trms/d(1)%jacob(1,1,1)
      urmss = urmss/d(1)%jacob(1,1,1)
      vrmss = vrmss/d(1)%jacob(1,1,1)
      wrmss = wrmss/d(1)%jacob(1,1,1)
      eprms = eprms/d(1)%jacob(1,1,1)
      amarms = amarms/d(1)%jacob(1,1,1)
   !
      statall(1) = urmss
      statall(2) = vrmss
      statall(3) = wrmss
      statall(4) = rrms
      statall(5) = trms
      statall(6) = eprms
      statall(7) = amarms
   !
      CALL MPI_ALLREDUCE(statall,statalla,10,MPI_DOUBLE_PRECISION, MPI_SUM, &
                             comm1d, ierr)
   !
      IF (myid==0) THEN
        urmss = statalla(1)
        vrmss = statalla(2)
        wrmss = statalla(3)
        rrms = statalla(4)
        trms = statalla(5)
        eprms = statalla(6)
        amarms = statalla(7)

        turbk= 0.5d0*(urmss+vrmss+wrmss)/(8.0d0*pi**3)
        urmss = sqrt(urmss/(8.0d0*pi**3))
        vrmss = sqrt(vrmss/(8.0d0*pi**3))
        wrmss = sqrt(wrmss/(8.0d0*pi**3))
        rrms = sqrt(rrms/(8.0d0*pi**3))
        trms = sqrt(trms/(8.0d0*pi**3))
        eprms = sqrt(eprms/(8.0d0*pi**3))
        amarms = sqrt(amarms/(8.0d0*pi**3))
      

        ! WRITE(20,*) time,uav,vav,wav,tav,rav,eav,epav,ediss1,ediss2,sdiss
        ! WRITE(22,*) time,urms,vrms,wrms,turbk,trms,rrms,eprms,amarms
        
      ENDIF
      DEALLOCATE(diss)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ufour=(0.,0.)
      vfour=(0.,0.)
      wfour=(0.,0.)
      npoin=0
      DO i=1,nfour
        DO j=1,nfour
          DO k=1,nfour
           xfour =DBLE(i-1)*2.0d0*pi/DBLE(nfour-1)
           yfour =DBLE(j-1)*2.0d0*pi/DBLE(nfour-1)
           zfour =DBLE(k-1)*2.0d0*pi/DBLE(nfour-1)
           ! write(100+myid,*)xfour,yfour,zfour,i,j,k
           DO id=1,ngrid
             IF (xfour>=d(id)%corner(1,1) .and. xfour<=d(id)%corner(1,2) .and. &
                 yfour>=d(id)%corner(2,1) .and. yfour<=d(id)%corner(2,3) .and. &
                 zfour>=d(id)%corner(3,1) .and. zfour<=d(id)%corner(3,5) ) THEN
             npoin(id) = npoin(id) + 1
             npoin1  = npoin(id)
             som =0.0d0
             som1=0.0d0
             som2=0.0d0
             xm=(xfour-d(id)%corner(1,1))/(d(id)%corner(1,2)-d(id)%corner(1,1))
             ym=(yfour-d(id)%corner(2,1))/(d(id)%corner(2,3)-d(id)%corner(2,1))
             zm=(zfour-d(id)%corner(3,1))/(d(id)%corner(3,5)-d(id)%corner(3,1))
             DO m=1,ncg 
               DO n=1,ncg
                  DO o=1,ncg
                    u    = d(id)%Q(m,n,o,2)/d(id)%Q(m,n,o,1)
                    v    = d(id)%Q(m,n,o,3)/d(id)%Q(m,n,o,1)
                    w    = d(id)%Q(m,n,o,4)/d(id)%Q(m,n,o,1)
                    ww   = ima(npoin1,id,m)* &
                           imb(npoin1,id,n)* &
                           imc(npoin1,id,o)
                    som  = som  + u*ww
                    som1 = som1 + v*ww
                    som2 = som2 + w*ww
                  ENDDO
               ENDDO
             ENDDO
             ufour(i,j,k)=CMPLX(som,0)
             vfour(i,j,k)=CMPLX(som1,0)
             wfour(i,j,k)=CMPLX(som2,0)
             ENDIF
           ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL MPI_ALLREDUCE(ufour,udum,nfour*nfour*nfour,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm1d, ierr)
      IF (myid==0) ufour=udum
      IF (myid==0) u_test=REAL(ufour)
      CALL MPI_ALLREDUCE(vfour,udum,nfour*nfour*nfour,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm1d, ierr)
      IF (myid==0) vfour=udum
      IF (myid==0) v_test=REAL(vfour)
      CALL MPI_ALLREDUCE(wfour,udum,nfour*nfour*nfour,MPI_DOUBLE_PRECISION, MPI_SUM, &
                         comm1d, ierr)
      IF (myid==0) wfour=udum
      IF (myid==0) w_test=REAL(wfour)

       IF (myid==0) THEN
      !   no=INT(time*100.0) 
      !   WRITE(fname, fmt='(a4,i3,a3)') 'fort',no,'.23'
      !   OPEN(unit=23,file=fname)
           CALL compute_espec(ufour,vfour,wfour,nfour,nfour,nfour,nfour,time)
      !   CLOSE(23)
       ENDIF
ENDIF
      ! DEALLOCATE(subdiss)
      ! DEALLOCATE(ima)    !!!!!!!!!!!!!
      ! DEALLOCATE(imb)    !!!!!!!!!!!!!
      ! DEALLOCATE(imc)    !!!!!!!!!!!!

     RETURN
      END SUBROUTINE do_user_stuff
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !//////////////////////////////////////////////////////////////////////////

 SUBROUTINE compute_diss(d,id)
 !
       USE domain_definition
       USE User_Data
       USE constants
       USE physics
       USE mpi_par
 !
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
 !
       TYPE (domain) :: d
 !
       DOUBLE PRECISION, DIMENSION(nx,ny,nz,5) :: Qx, Qy, Qz     
       INTEGER                 :: id
 !     storage for sub-grid dissipation calculation
       DOUBLE PRECISION :: u,v,w,rho,p,div
       DOUBLE PRECISION :: ux,uy,uz,vx,vy,vz,wx,wy,wz
       DOUBLE PRECISION :: S11,S22,S33,S12,S13,S23,SRT

       INTEGER          :: i,j,k,n,m,l

 !
 ! compute derivatives
      DO nv = 1,4
          CALL Gauss_Deriv(d%gmetg,d%Qlgg(:,:,:,nv),d%Qglg(:,:,:,nv), &
                          d%Qggl(:,:,:,nv),Qx(:,:,:,nv),Qy(:,:,:,nv),Qz(:,:,:,nv), &
                          d%ncg,d%dmx,d%dmy,d%dmz)
      END DO
                                                                       
       DO k = 1,d%ncg(3)
         DO j = 1,d%ncg(2)
          DO i = 1,d%ncg(1)
           rho = d%Q(i,j,k,1)*d%jacob(i,j,k)
           u   = d%Q(i,j,k,2)/d%Q(i,j,k,1)
           v   = d%Q(i,j,k,3)/d%Q(i,j,k,1)
           w   = d%Q(i,j,k,4)/d%Q(i,j,k,1)  
           p   = (d%Q(i,j,k,5)*d%jacob(i,j,k)-0.50d0*rho*(u**2+v**2+w**2))*(gamma-1.0d000)-twall/gamma                                                        
           ux = (Qx(i,j,k,2) - u*Qx(i,j,k,1))/rho
           uy = (Qy(i,j,k,2) - u*Qy(i,j,k,1))/rho
           uz = (Qz(i,j,k,2) - u*Qz(i,j,k,1))/rho
           vx = (Qx(i,j,k,3) - v*Qx(i,j,k,1))/rho
           vy = (Qy(i,j,k,3) - v*Qy(i,j,k,1))/rho
           vz = (Qz(i,j,k,3) - v*Qz(i,j,k,1))/rho
           wx = (Qx(i,j,k,4) - w*Qx(i,j,k,1))/rho
           wy = (Qy(i,j,k,4) - w*Qy(i,j,k,1))/rho
           wz = (Qz(i,j,k,4) - w*Qz(i,j,k,1))/rho
           
           om1 = wy - vz
           om2 = uz - wx
           om3 = vx - uy
           div = ux+vy+wz
           om1 = om1**2
           om2 = om2**2
           om3 = om3**2
                                                                        
           S11  = ux
           S22  = vy
           S33  = wz
           S12  = 0.5d0*(uy+vx)
           S13  = 0.5d0*(uz+wx)
           S23  = 0.5d0*(vz+wy)

           SRT    =  2*(S11*S11+S22*S22+S33*S33 + &
                       2*S12*S12+2*S13*S13+2*S23*S23)
                   
           diss(id,i,j,k,1)= (om1+om2+om3)/re
           diss(id,i,j,k,2)= -p*div!4.0d0*(ux+vy+wz)**2/(3.0d0*re)
           diss(id,i,j,k,3)= SRT/re
                                                                        
          END DO
        END DO
       END DO
                                                                       
       RETURN
       END SUBROUTINE compute_diss
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !//////////////////////////////////////////////////////////////////////////////////

       SUBROUTINE compute_espec(u1,u2,u3,nnx,nny,nnz,nfour,tt)
 !
       USE constants
       USE physics

 !
       parameter(nknum=256)
       integer   :: nnx,nyy,nnz,nfour,i,j,k,ii,jj,kk,isign,ndat,KMAX
       complex   :: u1(1:nfour,1:nfour,1:nfour)
       complex   :: u2(1:nfour,1:nfour,1:nfour)
       complex   :: u3(1:nfour,1:nfour,1:nfour)
       complex   :: s7(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! dill. Fourier coef
       complex   :: s6(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! vmag  Fourier coef
       complex   :: s3(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! u-vel  Fourier coef
       complex   :: s2(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! v-vel  Fourier coef
       complex   :: s1(-nfour/2+1:nfour/2,-nfour/2+1:nfour/2,-nfour/2+1:nfour/2)  ! w-vel  Fourier coef
       real      :: espec(nknum),espec1(nknum),espec2(nknum),espec3(nknum)
       real      :: especd(nknum),dspec(nknum)
       real      :: ur(3),ui(3),an1(3)
       double precision :: tt
 ! 
 !     ----------------------------------------
 !     tranfourrm the velocities to Fourier space
 !     ----------------------------------------
 !
       s6=(0.,0.)
       s3=(0.,0.)
       s2=(0.,0.)
       s1=(0.,0.)

       isign=-1
       ndat=nnx*nny*nnz
       ndim=3
       call fft(s1,u1,ndim,ndat,nnx,nny,nnz,isign)
       call fft(s2,u2,ndim,ndat,nnx,nny,nnz,isign)
       call fft(s3,u3,ndim,ndat,nnx,nny,nnz,isign)
 !

       do k=-nnz/2+1,nnz/2
         do j=-nny/2+1,nny/2
           do i=-nnx/2+1,nnx/2
 !          write(28,*) u3(i+nnx/2,j+nnx/2,k+nnx/2) 
       an1(1)   = real(i)
       an1(2)   = real(j)
       an1(3)   = real(k)
       ur(1)    = REAL(s1(i,j,k))    
       ur(2)    = REAL(s2(i,j,k))    
       ur(3)    = REAL(s3(i,j,k))      
       ui(1)    = AIMAG(s1(i,j,k))    
       ui(2)    = AIMAG(s2(i,j,k))    
       ui(3)    = AIMAG(s3(i,j,k))    
       dum      = 0.
       wmag     = an1(1)**2 + an1(2)**2 + an1(3)**2
       do ii =1,3
         do jj=1,3
           dum  = dum+an1(ii)*an1(jj)*(ur(ii)*ur(jj)+ui(ii)*ui(jj))
         enddo
       enddo
       dum      = dum/wmag
       s7(i,j,k)=cmplx(dum,0.)
       s1(i,j,k)=cmplx(real(s1(i,j,k))**2 &
                 +aimag(s1(i,j,k))**2,0.) 
       s2(i,j,k)=cmplx(real(s2(i,j,k))**2 &
                 +aimag(s2(i,j,k))**2,0.)
       s3(i,j,k)=cmplx(real(s3(i,j,k))**2 &
                 +aimag(s3(i,j,k))**2,0.)
       s6(i,j,k)=s1(i,j,k)+s2(i,j,k)+s3(i,j,k)
           enddo
         enddo
       enddo
 !
 !     -------------------------------------
 !     compute the spectra in Fourrier space
 !     -------------------------------------
 !
       espec=0.
       espec1=0.
       espec2=0.
       espec3=0.
       especd=0.
       const=4.*pi/3.
       KMAX =nfour
       KMAX=INT(SQRT(2.)*REAL(nfour)/3.)                                 
       do 40 knum=1,KMAX
       r=float(knum)
       num=0
       do 30 k=-nnz/2+1,nnz/2
       do 30 j=-nny/2+1,nny/2
       do 30 i=-nnx/2+1,nnx/2
       w12=real(i)*real(i)
       w22=real(j)*real(j)
       w32=real(k)*real(k)
       wmag=sqrt(w12+w22+w32)
       if((wmag.le.r-.5).or.(wmag.gt.r+.5)) go to 30
 20  continue
       num=num+1
       espec1(knum)=espec1(knum)+real(s1(i,j,k))
       espec2(knum)=espec2(knum)+real(s2(i,j,k))
       espec3(knum)=espec3(knum)+real(s3(i,j,k))
       espec(knum)=espec(knum)+real(s6(i,j,k))
       especd(knum)=especd(knum)+real(s7(i,j,k))
 30  continue
       if(num.eq.0) go to 40
       corec=float(num)/(const*((r+0.5)**3-(r-0.5)**3))
       if(knum.eq.1) corec=float(num)/(const*(1.5)**3)
       espec(knum) =espec(knum)/corec
       espec1(knum)=espec1(knum)/corec
       espec2(knum)=espec2(knum)/corec
       espec3(knum)=espec3(knum)/corec
       especd(knum)=especd(knum)/corec
 40  continue

       do 50 k=1,KMAX
       dspec(k)=2.0*espec(k)*REAL(k)**2 /re
       !write(23,*) k,espec(k)
        write(23,*) k,espec(k),espec1(k),espec2(k),espec3(k),especd(k),dspec(k)
 50  continue
 !
 !     -------------------------------------
 !     postprocess the spectra to obtain 
 !     turbulence statistics
 !     -------------------------------------
 !
       sum1=0.0
       sum2=0.0
       sum3=0.0
       sum4=0.0
       do 70 k=1,KMAX
       sum1=sum1+espec(k)
       sum2=sum2+dspec(k)
       sum3=sum3+espec(k)/float(k)
       sum4=sum4+especd(k)
 70    continue
       q=sqrt(2.0*sum1)
       qk=sum1
       qd=sum4
       uvar=q**2 /3.
       epsi=sum2
       xlength=3.0*pi*sum3/(2.0*q**2)
       xlam=sqrt(15.0*uvar/(re*epsi))
       eta=(1./(epsi*re**3))**(1./4.)
       taue=xlength/q*sqrt(3.)
       tauen=qk/epsi
       relamv=sqrt(uvar)*xlam*re
       ! write(90,*) tt,q,qd,qd/qk,epsi,xlength,xlam,eta,taue,tauen

       RETURN
       END SUBROUTINE compute_espec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!/////////////////////////////////////////////////////////////////////////


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE cleanup(time,ngrid,d)
!
!......................................................................
!     perform problem dependent computations
!......................................................................
!
      USE domain_definition
      USE physics
      USE constants
      USE User_Data
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain) :: d(ngp)
!
!     ----------------------------
!     do nothing
!     ----------------------------
!
!
      RETURN
      END SUBROUTINE cleanup
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE initc(id,d,mrtr)
!
!......................................................................
!     Set up the initial conditions as a function of x and y
!     these conditions are for a point source solution
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d
!   type(bounds)             :: dbnd
      TYPE (mortar),DIMENSION(nmp) :: mrtr
      DOUBLE PRECISION             :: fv(5)
      INTEGER                      :: neighbor
      double precision              :: u,v,w,gOffset
    
!
!     -------------------------
!     Inputs 
!     -------------------------
!



!
!     ----------------------------
!     compute interior point values
!     ----------------------------
!
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               x = d%xg(1,n,m,l)
               y = d%xg(2,n,m,l)
               z = d%xg(3,n,m,l)
               CALL exactsol(x,y,z,d%Q(n,m,l,1),d%Q(n,m,l,2), &
                             d%Q(n,m,l,3),d%Q(n,m,l,4),d%Q(n,m,l,5),fv(:))
           !d%Q(n,m,l,1) = rho
           !d%Q(n,m,l,2) = u*rho
           !d%Q(n,m,l,3) = v*rho
           !d%Q(n,m,l,4) = w*rho
           !d%Q(n,m,l,5) = press/(gamma-1.0d0)+0.5d0*d%Q(n,m,l,1)*(u**2+v**2+w**2 )
           !fv(:) = 0.0d0
            END DO 
         END DO
      END DO
!
!     --------------------------
!     compute boundary solutions
!     --------------------------
!
!     ----------
!     face 1 & 2
!     ----------
!
      DO iface = 1,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         mrtr_no = nmortmpi(mrtr_no,1)
         neighbor  = mrtr(mrtr_no)%id(2)
         IF ( neighbor /= 0 )     CYCLE
         DO l = 1,d%ncg(3)
            DO n = 1,d%ncg(1)
               x = mrtr(mrtr_no)%xg(1,n,l)
               y = mrtr(mrtr_no)%xg(2,n,l)
               z = mrtr(mrtr_no)%xg(3,n,l)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(n,l,1,2),mrtr(mrtr_no)%Q(n,l,2,2), &
                                   mrtr(mrtr_no)%Q(n,l,3,2),mrtr(mrtr_no)%Q(n,l,4,2), &
                                   mrtr(mrtr_no)%Q(n,l,5,2),mrtr(mrtr_no)%fv(n,l,:,2))
               !mrtr(mrtr_no)%Q(n,l,1,2) = rho
               !mrtr(mrtr_no)%Q(n,l,2,2) = u*rho
               !mrtr(mrtr_no)%Q(n,l,3,2) = v*rho
               !mrtr(mrtr_no)%Q(n,l,4,2) = w*rho
               !mrtr(mrtr_no)%Q(n,l,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(n,l,:,2) = 0
            END DO
         END DO
      END DO
!
!     ----------
!     face 3 & 5
!     ----------
!
      DO iface = 3,5,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         mrtr_no = nmortmpi(mrtr_no,1)
         neighbor  = mrtr(mrtr_no)%id(2)
         IF ( neighbor /= 0 )     CYCLE
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               x = mrtr(mrtr_no)%xg(1,n,m)
               y = mrtr(mrtr_no)%xg(2,n,m)
               z = mrtr(mrtr_no)%xg(3,n,m)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(n,m,1,2),mrtr(mrtr_no)%Q(n,m,2,2), &
                                   mrtr(mrtr_no)%Q(n,m,3,2),mrtr(mrtr_no)%Q(n,m,4,2), &
                                   mrtr(mrtr_no)%Q(n,m,5,2),mrtr(mrtr_no)%fv(n,m,:,2))

               !mrtr(mrtr_no)%Q(n,m,1,2) = rho
               !mrtr(mrtr_no)%Q(n,m,2,2) = u*rho
               !mrtr(mrtr_no)%Q(n,m,3,2) = v*rho
               !mrtr(mrtr_no)%Q(n,m,4,2) = w*rho
               !mrtr(mrtr_no)%Q(n,m,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(n,m,:,2) = 0.0d0
            END DO
         END DO
      END DO
!
!     ----------
!     face 4 & 6
!     ----------
!
      DO iface = 4,6,2
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side == 2 )     CYCLE
         mrtr_no = nmortmpi(mrtr_no,1)
         neighbor  = mrtr(mrtr_no)%id(2)
         IF ( neighbor /= 0 )     CYCLE
         IF (iface==6) THEN
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               x = mrtr(mrtr_no)%xg(1,m,l)
               y = mrtr(mrtr_no)%xg(2,m,l)
               z = mrtr(mrtr_no)%xg(3,m,l)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(m,l,1,2),mrtr(mrtr_no)%Q(m,l,2,2), &
                                   mrtr(mrtr_no)%Q(m,l,3,2),mrtr(mrtr_no)%Q(m,l,4,2), &
                                   mrtr(mrtr_no)%Q(m,l,5,2),mrtr(mrtr_no)%fv(m,l,:,2))

               !mrtr(mrtr_no)%Q(m,l,1,2) = rho
               !mrtr(mrtr_no)%Q(m,l,2,2) = u*rho
               !mrtr(mrtr_no)%Q(m,l,3,2) = v*rho
               !mrtr(mrtr_no)%Q(m,l,4,2) = w*rho
               !mrtr(mrtr_no)%Q(m,l,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(m,l,:,2) = 0.0d0
            END DO
         END DO
         ELSEIF (iface==4) THEN ! outflow boundary
         DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               x = mrtr(mrtr_no)%xg(1,m,l)
               y = mrtr(mrtr_no)%xg(2,m,l)
               z = mrtr(mrtr_no)%xg(3,m,l)
               CALL exactsol(x,y,z,mrtr(mrtr_no)%Q(m,l,1,2),mrtr(mrtr_no)%Q(m,l,2,2), &
                                   mrtr(mrtr_no)%Q(m,l,3,2),mrtr(mrtr_no)%Q(m,l,4,2), &
                                   mrtr(mrtr_no)%Q(m,l,5,2),mrtr(mrtr_no)%fv(m,l,:,2))

               !mrtr(mrtr_no)%Q(m,l,1,2) = rho
               !mrtr(mrtr_no)%Q(m,l,2,2) = u*rho
               !mrtr(mrtr_no)%Q(m,l,3,2) = v*rho
               !mrtr(mrtr_no)%Q(m,l,4,2) = w*rho
               !mrtr(mrtr_no)%Q(m,l,5,2) = press/(gamma-1.0d0)+0.5d0*rho*(u**2+v**2)
               !mrtr(mrtr_no)%fv(m,l,:,2) = 0.0d0
            END DO
         END DO
         ENDIF
      END DO
      
!
      RETURN
      END SUBROUTINE initc
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Xinitc(id,d,mrtr)
!
!......................................................................
!     Interpolate solution values to new grid
!......................................................................
!
      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d
      TYPE (mortar),DIMENSION(nmp) :: mrtr
      DOUBLE PRECISION             :: cxg(nx),Q(nx,ny,nz,neq),polyn
      INTEGER                      :: p
!
      Q = 0.0d0
      cxg = 0.0d0
!
!     Define old Chebyshev grid parameters
!
      ncg1 = 5
      ncg2 = 5
      ncg3 = 5
      cxg = 0.0d0
      DO n = 1,ncg1
         xx = cos((2.d0*DBLE(n)-1.d0)*pi/(2.d0*(DBLE(ncg1) - 1.d0) + 2.d0))
         cxg(n) = 0.5d0*(1.d0 - xx)
      END DO


!
!     Interpolate Q to refined grid within each domain
!
      DO i=1,d%ncg(1)
        DO j=1,d%ncg(2)
          DO k=1,d%ncg(3)
            X =d%cxg(i,1)
            Y =d%cxg(j,2)
            Z =d%cxg(k,3)
!
            DO nv=1,neq
              som  = 0.0d0
              DO m=1,ncg1
                DO n=1,ncg2
                  DO p=1,ncg3
                    hx=polyn(m,X,ncg1,cxg(:))
                    hy=polyn(n,Y,ncg2,cxg(:))
                    hz=polyn(p,Z,ncg3,cxg(:))
                    som  = som  + d%Q(m,n,p,nv)*hx*hy*hz
                  ENDDO
                ENDDO
              ENDDO
!
              Q(i,j,k,nv)=som
!
            ENDDO
!
          ENDDO
        ENDDO
      ENDDO
!
!   Multiply by jacobian: this one is scaled out again in the the main program
!
      DO l = 1,d%ncg(3)
         DO m = 1,d%ncg(2)
            DO n = 1,d%ncg(1)
               DO nv = 1,neq
                   d%Q(n,m,l,nv) = Q(n,m,l,nv)*d%jacob(n,m,l)
               END DO
            END DO
        END DO
      END DO
!
!
!
      RETURN
      END SUBROUTINE Xinitc
!
!///////////////////////////////////////////////////////////////////////
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE exactsol(x,y,z,rho,rhou,rhov,rhow,rhoe,fv)
!
!......................................................................
!     compute the initial for the channel flow with periodic bc
!......................................................................

      USE constants
      USE User_Data
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION :: fv(5)
      DOUBLE PRECISION :: hlow, hupp, uu, vv, ww, p, rho1, som, initvel, Tref
      DOUBLE PRECISION :: xx(16),uuu(16)
      DOUBLE PRECISION :: xx1(38),uuu1(38)
      DOUBLE PRECISION :: polyn, hz
      double precision, parameter :: pie = 3.141592653589793238462643383279502d0
      double precision  :: radius, smooth
      double precision  :: rr,Machy,Temp
      
      V0=1.0d0!6.2900
      uu   = V0*sin(x)*cos(y)*cos(z)
      vv   = -V0*cos(x)*sin(y)*cos(z)
      ww   = 0.0d0
      umax = 1.5d0
      
      twall = twall
      rho1 = 1.0d0
      p    = twall*rho1/gamma+(cos(2.0*x)+cos(2.0*y))*(cos(2.0*z)+2.0)/16.0d0 !twall*rho1/gamma!
      
!
      rho  = rho1
      rhou = rho*uu
      rhov = rho*vv
      rhow = rho*ww
      rhoe = p/(gamma-1.0d0) + 0.5d0*rho*(uu**2 + vv**2 + ww**2)
      fv(:) = 0.0d0
      END SUBROUTINE exactsol
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE compute_source(x,y,z,time,source,delta_t)
!
!......................................................................
!     physical source term
!......................................................................
!
      USE size
      USE physics
      USE User_Data
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION, DIMENSION(neq) :: source
      double precision                  :: coef, xshift, yshift, scale
!
!         coef    = 10.0d0
!         scale   = 0.020d0
!         xshift  = 1.0d0
!         yshift  = 1.0d0
!
!         source(1) = 0.0d0
!         source(2) = coef * exp( -( (x-xshift)**2 + (y-yshift)**2 ) / scale )
!         source(3) = 0.0d0
!         source(4) = 0.0d0
!         source(5) = 10.0d0*coef * exp( -( (x-xshift)**2 + (y-yshift)**2 ) / scale )
      RETURN
      END SUBROUTINE compute_source
!
!///////////////////////////////////////////////////////////////////////
!
!
   SUBROUTINE boundary_flux(time,mrtr,bcond,property,num_props,flux)
!
!......................................................................
!     date: 11/17/98
!
!     routines called:
!                      wallbc
!                      inflow
!                      outflow
!
!     applicability: all
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: flux
      DOUBLE PRECISION, DIMENSION(mxprops)       :: property
      CHARACTER(LEN = 9)                         :: bcond
!
!     ------------------------------------------------------
!     call appropriate boundary conditions for this problem
!     -----------------------------------------------------
!
      SELECT CASE(bcond)

         CASE('walladiab')
            CALL wallbc_Inv(mrtr,flux)
         CASE('wallisoth')
            CALL wallbc_isoth(mrtr,flux)
         CASE('inflow')
            CALL inflow_specAll(mrtr,flux)
         CASE('outflow')
            CALL outflow_specAll(mrtr,flux)
         CASE('periodm')
            CALL periodm_specAll(mrtr,flux)
         CASE('periods')
            CALL periods_specAll(mrtr,flux)
      END SELECT
!
      RETURN
      END SUBROUTINE boundary_flux

!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE DirichletBC(time,mrtr,bcond,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     routines called:
!                      boundary_values
!                      set_adiabatic wall
!
!     applicability: 3D Navier-Stokes
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      CHARACTER(LEN = 9)                         :: bcond
!
!     ------------------------------------------------------
!     call appropriate boundary conditions for this problem
!     -----------------------------------------------------
!
      SELECT CASE(bcond)

         CASE('wallisoth')
            CALL set_isothermal_wall(mrtr,Qbound)
         CASE('walladiab')
            CALL set_adiabatic_wall(mrtr,Qbound)
         CASE('inflow')
            CALL boundary_values(mrtr,Qbound)
         CASE('outflow')
            CALL boundary_values(mrtr,Qbound)
         CASE('periodm')
            CALL set_periodm(mrtr,Qbound)
         CASE('periods')
            CALL set_periods(mrtr,Qbound)
      END SELECT
!
      RETURN
      END SUBROUTINE DirichletBC
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Vis_Neumann(time,mrtr,bcond,Qbound)
!
!......................................................................
!     date: 03/12/01
!
!     routines called:
!                      boundary_values
!                      set_adiabatic wall
!
!     applicability: 3D Navier-Stokes
!......................................................................
!
      USE  mortar_definition
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (mortar)                              :: mrtr
      DOUBLE PRECISION, DIMENSION(nmax,nmax,neq) :: Qbound
      CHARACTER(LEN = 9)                         :: bcond
!
!     ------------------------------------------------------
!     call appropriate boundary conditions for this problem
!     -----------------------------------------------------
!
      SELECT CASE(bcond)

         CASE('wallisoth')
            CALL wall_vis_flux(mrtr,Qbound,bcond)
         CASE('walladiab')
            CALL wall_vis_flux(mrtr,Qbound,bcond)
         CASE('inflow')
            CALL boundary_vis_flux(mrtr,Qbound)
         CASE('outflow')
            CALL boundary_vis_flux(mrtr,Qbound)
      END SELECT
!
      RETURN
      END SUBROUTINE Vis_Neumann
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE set_part_par()

!
!    initialization of particle parameters
!
      USE size
      USE constants
      USE part_par
      USE physics
!
!
!     parameters not defined in the euler code
!
      IEVAP  = 0

      SCF    = 0.7d0
      SCO    = 0.7d0
!      DA     = 6.0d0
!      ZE     = 4.0d0
!      CE     = 40.0d0
      RCOEF  = 1.0d0
      YFINIT = 0.0d0
!
!     compute particle parameters (see formulation)
!
      PEF    = SCF*re
      PD0    = dsqrt(18.0d0*TAUP0/(re*RHOP))
      PM0    = pi*RHOP*PD0**3 /6.d0
      PSI    = dble(nprt)*PM0/(2.d0*pi)**3
      CRE    = re*(6.d0/(pi*RHOP))**(1.0d0/3.0d0)
      CTAUP  = (RHOP**(1.0d0/3.0d0) *re/18.0d0)*(6.0d0/pi)**(2.0d0/3.0d0)
      CYFPS  = gamma*A1/((gamma-1.)*TBOIL)
      CF2    = 1.0d0/(3.0d0*pr*A2)
      CF3    = A1/(3.0d0*SCF*A2)
      CF4    = pi*dsqrt(18.0d0/RHOP)/(re**(3.0d0/2.0d0) *SCF)

!   
      nprtmax = 200000

      RETURN
   END SUBROUTINE set_part_par
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE part_initc(drop,nprt)

!
!    initialization of particle parameters
!
      USE size
!
      USE particle_definition
!
      USE constants
      USE part_par
      USE physics
      USE User_Data
!
      TYPE(particle)  :: drop
!
!    initialize particle field
!
      DO j=1,3
         drop%Vp(j)     = 0.0d0
         drop%Vfp(j)    = 0.0d0
      END DO
         drop%Tp       = 14.0d0
         drop%Mp       = PM0
         drop%Tfp      = 14.0d0
         drop%Rhofp    = 1.0d0
         drop%Yffp     = 1.0d0
!
         drop%onoff     = 1
!

      RETURN
   END SUBROUTINE part_initc
!
!///////////////////////////////////////////////////////////////////////
!
   SUBROUTINE find_processor(xmid,ymid,zmid,id, nno,ngrid)
!
       USE constants
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
! structured uniform grid in x-direction, and y-direction
! cosine distributiono in z-directions
           CALL read_metis_file(id, nno,ngrid)
      ! ELSE
      !     domlx = 4.0d0
      !     domly = 2.0d0
      !     domlz = 2.0d0
      !     nprocx = 0
      !     nprocy = 1
      !     nprocz = 1
      !     zmid   = 2.0d0*(ACOS(1.0d0-2.0d0*(zmid/2.0d0))/pi)
!
!           n1   = INT(xmid*REAL(nprocx)/domlx)+1
!           n2   = INT(ymid*REAL(nprocy)/domly)+1
!           n3   = INT(zmid*REAL(nprocz)/domlz)+1
!           nno  = (n1-1)*nprocy+n2-1+(n3-1)*nprocx*nprocy  
!
!           IF (ymid < 3.0) THEN
!             nno = 0
!           ELSE
!             nno = 1
!           ENDIF
!     ENDIF

!
      RETURN
   END SUBROUTINE find_processor
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE read_restart_case(rst_unit,time)
!
! read case specific stuff from the restart file
!
       USE User_Data
       USE mpi_par
!
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       INTEGER :: rst_unit
!
       INCLUDE "mpif.h"
!
      IF (myid == 0) THEN
         READ(rst_unit) time
         write(*,*) 'time',time
      ENDIF
!
      CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
!

       RETURN
      END SUBROUTINE read_restart_case
 
!
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE write_restart_case(rst_unit,time)
!
! Write case specific stuff to the restart file
!
       USE User_Data
       USE mpi_par
!       
       IMPLICIT DOUBLE PRECISION(a-h,o-z)
!
       INTEGER :: rst_unit
!
       INCLUDE "mpif.h"
    
      IF (myid == 0) THEN
         WRITE(rst_unit) time
         write(*,*) 'final time',time
      ENDIF
!
       RETURN
      END SUBROUTINE write_restart_case


!/////////////////////////////////////////////////////////////
      !! Particle inflow subroutines
subroutine partInflow(d,p,dt,ngridl,nprt,myid)
    use domain_definition
    use size 
    use physics
    use particle_definition
    use particleInlets
    implicit none
    ! Domain Types
    type(domain), dimension(ngp)                 :: d            !< Carrier Domain
    ! Particle Types
    type(particle), dimension(npart)             :: p            !< Particle
    
    double precision                             :: uniformRandom
    
    integer, intent(in)                          :: ngridl       !< Local elements
    integer, intent(inout)                       :: nprt
    integer, intent(in)                          :: myid
    double precision, intent(in)                 :: dt
    integer                                      :: i,j,k,id,iface,padd,np
    double precision                             :: u, paddD, yposm, ypos, y1, y2, dy, x1, x2, dx, xpos
    double precision                             :: Fui, Oxi, Nii, Pdi, Fuc, Oxc, Nic, Pdc
    double precision, parameter                  :: Lx1 = 0.38448375335127d0 !< length of the first element
    double precision, parameter                  :: Lx2 = 0.04250d0
    integer, parameter                           :: ideal = 80
    

    
end subroutine partInflow

subroutine partOutflow(d,p,dt,ngridl,nprt,myid)
    use domain_definition
    use size 
    use physics
    use particle_definition
    use particleInlets
    implicit none
    ! Domain Types
    type(domain), dimension(ngp)                 :: d            !< Carrier Domain
    ! Particle Types
    type(particle), dimension(npart)             :: p            !< Particle
    
    integer, intent(in)                          :: ngridl       !< Local elements
    integer, intent(inout)                       :: nprt
    integer, intent(in)                          :: myid
    double precision, intent(in)                 :: dt
    
    integer                                      :: np
    
    
    
end subroutine partOutflow

subroutine initScalars(d,eP,id)
    use size
    use particle_definition
    use domain_definition
    implicit none
    type(domain)                                :: d
    type(ePointers), dimension(ngp,npart)       :: eP
    double precision                            :: x,y,z
    integer                                     :: np, id
    double precision, parameter                 :: pie		= 3.14159265358979323846264338327d0
    double precision                            :: delta
    double precision                            :: T,F,O,N,P,spike
    
    
    do np=1,d%np
        x = eP(id,np)%p%Xp(1)
        y = eP(id,np)%p%Xp(2)
        z = eP(id,np)%p%Xp(3)
        
		!For CV Reactor
        !eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        !eP(id,np)%p%scalar(2) = 0.01152d0
        !eP(id,np)%p%scalar(3) = 0.23042d0
        !eP(id,np)%p%scalar(4) = 0.0d0
        !ep(id,np)%p%scalar(5) = 0.75806d0
        
!         !for CV reactor
!         F     = 0.01152d0
!         O     = 0.23042d0
!         N     = 0.75806d0
!         P   = 0.0d0
        !for ramp cavity with injector
        ! F     = 0.000001d0
!         O     = 0.24194d0
!         N     = 0.75806d0
!         P     = 0.0d0

        !for Sod Problem
        if (x > 0.50d0) then
            F     = 0.0115d0
            O     = 0.2304d0
            N     = 0.7581d0
            P   = 0.0d0
        else
            F = 0.0d0
            O = 0.23310d0
            N = 0.76690d0
            P = 0.0d0
        end if

        eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        eP(id,np)%p%scalar(2) = F
        eP(id,np)%p%scalar(3) = O
        eP(id,np)%p%scalar(4) = P
        eP(id,np)%p%scalar(5) = N
        eP(id,np)%p%react = 0.0d0

        ! for ramp cavity
        !T = eP(id,np)%p%Tfp
        !
        !if(x.le.3) then
        !    !T = 1.0d0*(spike-(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
        !    F = 0.0d0
        !    O = 0.2331d0
        !    N = 0.7669d0
        !    P = 1.0d0-F-O-N
        !else if(x.ge.6) then
        !    F = 0.0d0
        !    O = 0.2331d0
        !    N = 0.7669d0
        !    P = 1.0d0-F-O-N
        !!else if(y.le.dsqrt(1.50d0**2.0-(x-4.50d0)**2.0)+0.650d0) then
        !else if(y.le.0.650d0) then
        !    !T = 1.0d0*(spike+(spike-1.0d0)*erf(pie**0.5d0*(y)/delta))
        !    F = 0.0115d0
        !    O = 0.2304d0
        !    N = 0.7581d0
        !    P = 1.0d0-F-O-N
        !else
        !    F = 0.0d0
        !    O = 0.2331d0
        !    N = 0.7669d0
        !    P = 1.0d0-F-O-N 
        !end if
        !
        !eP(id,np)%p%scalar(1) = eP(id,np)%p%T
        !eP(id,np)%p%scalar(2) = F
        !eP(id,np)%p%scalar(3) = O
        !eP(id,np)%p%scalar(4) = P
        !eP(id,np)%p%scalar(5) = N
        ! end of for cavity

        !eP(id,np)%p%scalar(1) = 3.0
        !eP(id,np)%p%scalar(2) = 4.0
        
        !if (y .ge. 0.0d0) then
        !    ep(id,np)%p%scalar(1) = erf(3.14159d0**0.5d0 * y / delta)
        !else
        !    ep(id,np)%p%scalar(1) = 0.0d0
        !end if
        
        !if (y .ge. 0.0d0) then
        !    eP(id,np)%p%scalar(1) = 0.0d0
        !    ep(id,np)%p%scalar(2) = 1.0!erf(sqrt(pie)*y/delta)
        !end if
        !if (y .lt. 0.0d0) then
        !    eP(id,np)%p%scalar(2) = 0.0d0
        !    ep(id,np)%p%scalar(1) = 1.0!- erf(sqrt(pie)*y/delta)
        !end if
        !
        !eP(id,np)%p%scalar(3) = eP(id,np)%p%Rhofp
        
    end do
    
end subroutine initScalars

!------------------------------------
!-------------------------------------
!--- Stochastic inlet starts here-----
!-------------------------------------
!-------------------------------------
!-------------------------------------

!///////////////////////////////////////////////////////////////////

      SUBROUTINE inflow_time_bc(id,d,mrtr,k,dt)

      USE domain_definition
      USE mortar_definition
      USE physics
      USE constants
      USE User_Data
      USE mpi_par
      USE mpi
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      TYPE (domain)                :: d
      TYPE (mortar),DIMENSION(nmp) :: mrtr
      DOUBLE PRECISION             :: fv(5)
      INTEGER                      :: neighbor
      
      hlow     = 1.0d0 ! height of profile
      neighbor = 1
      utau     = 1.0d0!0.1596d0 / 2.0d0 ! 0.1596d0 ! intensity

      IF(k==1 .and. id==1) THEN
        IF(myid==0) THEN 
        
!        read mean profile
         OPEN(unit=22,file='flame-meanv.txt')
         DO i=1,25
           READ(22,*) xxin(i),uuin(i)
         ENDDO
         CLOSE(22)
!        read v rms profile
         OPEN(unit=224,file='flame-vrms.txt')
         DO i=1,25
           READ(224,*) xx2in(i),uu2in(i)
           !xx2in(i) = 1.2d0*xx2in(i) + hlow ! scale from z/delta99 to z/h
           uu2in(i) = utau*uu2in(i)
           write(*,*) 'vrms: ',xx2in(i),uu2in(i)
         ENDDO
         CLOSE(224)      
         
        END IF
        uu1in = 0.0d0
        uu3in = 0.0d0
! Broadcast now
        CALL MPI_BCAST(xxin,38,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(uuin,38,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(xx1in,94,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(uu1in,94,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(xx2in,61,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(uu2in,61,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(xx3in,61,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
        CALL MPI_BCAST(uu3in,61,MPI_DOUBLE_PRECISION,0,comm1d,ierr)
      END IF

     DO iface = 1,6
       IF (d%bcond(iface)  == 'inflow') THEN

!      DO iface = 4,6,2
!         iface = 6
         mrtr_no   = d%mortar(1,iface)
         mrtr_side = d%mortar(2,iface)
         IF ( mrtr_side /= 2 )  THEN
           mrtr_no = nmortmpi(mrtr_no,1)
           neighbor  = mrtr(mrtr_no)%id(2)
         END IF
         IF ( neighbor == 0 )  THEN
          DO l = 1,d%ncg(3)
            DO m = 1,d%ncg(2)
               x = mrtr(mrtr_no)%xg(1,m,l)
               y = mrtr(mrtr_no)%xg(2,m,l)
               z = mrtr(mrtr_no)%xg(3,m,l)
               IF (x<2.0d0) THEN ! specify the inflow area, avoid messing up with the injector.
                    CALL stoch_init(x,y,z,m,l,mrtr(mrtr_no)%Q(m,l,1,2), &
                           mrtr(mrtr_no)%Q(m,l,2,2), mrtr(mrtr_no)%Q(m,l,3,2),&
                           mrtr(mrtr_no)%Q(m,l,4,2))
                    CALL stoch(ufluct(m,l,:),vfluct(m,l,:),wfluct(m,l,:),urms(m,l), &
                          vrms(m,l),wrms(m,l),dt)
                    CALL inflow_cond(x,y,z,m,l,mrtr(mrtr_no)%Q(m,l,1,2),mrtr(mrtr_no)%Q(m,l,2,2), &
                                   mrtr(mrtr_no)%Q(m,l,3,2),mrtr(mrtr_no)%Q(m,l,4,2), &
                                   mrtr(mrtr_no)%Q(m,l,5,2),mrtr(mrtr_no)%fv(m,l,:,2),ufluct(m,l,2), &
                                   vfluct(m,l,2), wfluct(m,l,2),k)
                END IF
            END DO
          END DO
         END IF
       END IF
     END DO

      return
      END SUBROUTINE inflow_time_bc
!/////////////////////////////////////////////////////////////

      SUBROUTINE inflow_cond(x,y,z,m,l,rho,rhou,rhov,rhow,rhoe,fv,uprime,vprime,wprime,k)

      USE constants
      USE User_Data
      USE physics
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION :: fv(5)
      DOUBLE PRECISION :: hlow, hupp, uu, vv, ww, p, rho1, som
      DOUBLE PRECISION :: polyn, hz

      rho1     = 1.0d0
      hlow     = 1.0d0
      hupp     = 6.0d0
      p        = twall*rho1/gamma

!      CALL RANDOM_NUMBER(dum)
!      dum  = (dum-0.5d0)/2.0d0

!IF(MOD(k,10)==0 .OR. MOD(k-1,10)==0 .OR. k<3) THEN
!write(*,*)'X',x,y,z
!IF(z>0.9d0 .AND. z<1.1d0) write(*,*) 'probe',k,x,y,z,umean(m,l),uprime
!END IF
      uu = umean(m,l) + uprime
      vv = vmean(m,l) + vprime 
      ww = wmean(m,l) + wprime
                                                                                                               
!
      rho  = rho1
      rhou = rho*uu
      rhov = rho*vv
      rhow = rho*ww
      rhoe = p/(gamma-1.0d0) + 0.5d0*rho*(uu**2 + vv**2 + ww**2)
      fv(:) = 0.0d0

      return       
      END SUBROUTINE inflow_cond

!/////////////////////////////////////////////////////

      SUBROUTINE stoch_init(x,y,z,m,l,rho,rhou,rhov,rhow)

      USE User_Data
      
      IMPLICIT NONE
      INTEGER          :: m,l 
      DOUBLE PRECISION :: x,y,z
      INTEGER          :: i 
      DOUBLE PRECISION :: rho,rhou,rhov,rhow
      DOUBLE PRECISION :: u,v,w, uvel



      ! IF (x<2.0d0) THEN
! !       Spalart boundary layer urms profile
!            DO i=1,92
!              IF (y .ge. xx1in(i) .and. y.le.xx1in(i+1)) THEN
!                urms(m,l)=((y-xx1in(i))/(xx1in(i+1)-xx1in(i)))*(uu1in(i+1)-uu1in(i))+uu1in(i)
!              END IF
!            ENDDO
!       ELSE
!             urms(m,l)=uu1in(93)
!       END IF

      IF (x<2.0d0) THEN
!       Spalart boundary layer vrms profile 
           DO i=1,60
             IF (y .ge. xx2in(i) .and. y.le.xx2in(i+1)) THEN
               vrms(m,l)=((y-xx2in(i))/(xx2in(i+1)-xx2in(i)))*(uu2in(i+1)-uu2in(i))+uu2in(i)
             END IF
           ENDDO
      ELSE
           vrms(m,l)=uu2in(61) 
      END IF

!       IF (x<2.0d0) THEN
! !       Spalart boundary layer wrms profile
!            DO i=1,60
!              IF (y .ge. xx3in(i) .and. y.le.xx3in(i+1)) THEN
!                wrms(m,l)=((y-xx3in(i))/(xx3in(i+1)-xx3in(i)))*(uu3in(i+1)-uu3in(i))+uu3in(i)
!              END IF
!            ENDDO
!       ELSE
!            wrms(m,l)=uu3in(61)
!       END IF

!     Mean profile
      IF (x<2.0d0) THEN
!       Spalart boundary layer mean profile + random disturbance
           DO i=1,35
             IF (y .ge. xxin(i) .and. y.le.xxin(i+1)) THEN
               vmean(m,l)=((y-xxin(i))/(xxin(i+1)-xxin(i)))*(uuin(i+1)-uuin(i))+uuin(i)
             END IF
           ENDDO
      ELSE
           umean(m,l)=1.0d0
      END IF
      umean(m,l) = 0.0d0
      wmean(m,l) = 0.0d0 
      
      ! zia recommended this to fix problem
      ! uvel = 0.4d0
      !       urms(m,l) = urms(m,l) * uvel
      !       vrms(m,l) = vrms(m,l) * uvel
      !       wrms(m,l) = wrms(m,l) * uvel
      !       umean     = umean(m,l) * uvel
      
!     initialize uflucts
      u = rhou/rho
      v = rhov/rho
      w = rhow/rho
      ufluct(m,l,1) = u - umean(m,l) 
      vfluct(m,l,1) = v - vmean(m,l) 
      wfluct(m,l,1) = w - wmean(m,l) 

!      urms(m,l)=0.15d0
!      vrms(m,l)=0.15d0
!      wrms(m,l)=0.15d0

!       write(*,*) z, urms(m,l),vrms(m,l),wrms(m,l)

      return

      END SUBROUTINE stoch_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE stoch(uprime,vprime,wprime,urms,vrms,wrms,tstep)
       

        IMPLICIT NONE
          
        DOUBLE PRECISION  :: uprime(2),vprime(2),wprime(2),tstep
        DOUBLE PRECISION  :: urms,vrms,wrms
        DOUBLE PRECISION  :: beta(3)        
        DOUBLE PRECISION  :: l(3,3),d(3,3),z(3),d_t(3)
        DOUBLE PRECISION  :: sum1,sum2,uu,vv,ww,ke
        DOUBLE PRECISION  :: Ruu,Rvv,Rww,dum,T1,dum1,dum2,dum3
        INTEGER           :: i,j,k
        DOUBLE PRECISION  :: intensity,Tuu,Tvv,Tww

!        intensity = 0.025d0     
!       integral time scales
        Tuu = 1.08d0    
        Tvv = 0.366d0
        Tww = 0.394d0

        uu = urms**2
        vv = vrms**2
        ww = wrms**2
        
!        ke = (intensity**2)/1.07d0
!        CALL RANDOM_NUMBER(dum)
!        IF(dum > 0.0d0) ke = dum*ke

!       set the random vector z
            
        CALL RANDOM_NUMBER(dum1)
        z(1)=dum1-0.5d0
        CALL RANDOM_NUMBER(dum2)
        z(2)=dum2-0.5d0
        CALL RANDOM_NUMBER(dum3)
        z(3)=dum3-0.5d0
         
!        urms(2) = sqrt(1.07d0*ke)
!        vrms(2) = sqrt(0.37d0*ke)
!        wrms(2) = sqrt(0.56d0*ke)

!        Ruu = exp(-tstep/T1)   
!        Rvv = exp(-tstep/T1)   
!        Rww = exp(-tstep/T1)
   
        Ruu = exp(-tstep/Tuu)   
        Rvv = exp(-tstep/Tvv)   
        Rww = exp(-tstep/Tww)   
      
        beta(1) = Ruu
        beta(2) = Rvv
        beta(3) = Rww

        d(1,1) = urms**2-beta(1)*Ruu*urms**2
        d(2,2) = vrms**2-beta(2)*Rvv*vrms**2
        d(3,3) = wrms**2-beta(3)*Rww*wrms**2
        d(2,1) = 0.0d0
        d(3,1) = 0.0d0
        d(3,2) = 0.0d0

! Cholesky factorization for a 3x3 matrix 

        DO i=1,3
           sum1=0.0d0
           DO k=1,i-1
            sum1 = sum1 + l(i,k)**2
           END DO
           l(i,i) = sqrt(d(i,i)-sum1)
          DO j=i+1,3
            sum2=0.0d0
            DO k=1,i-1
             sum2 = sum2 + l(j,k)*l(i,k)
            END DO
            l(j,i) = (d(j,i) - sum2)/l(i,i)
          END DO
        END DO

        d_t(1) = l(1,1)*z(1)
        d_t(2) = l(2,1)*z(1) + l(2,2)*z(2)
        d_t(3) = l(3,1)*z(1) + l(3,2)*z(2) + l(3,3)*z(3)

! fluctuations at the current time step

        uprime(2) = beta(1)*uprime(1) + d_t(1)
        vprime(2) = beta(2)*vprime(1) + d_t(2)
        wprime(2) = beta(3)*wprime(1) + d_t(3)
 

!       update 
!        urms(1) = urms(2) 
!        vrms(1) = vrms(2) 
!        wrms(1) = wrms(2) 
!        uprime(1) = uprime(2)      
!        vprime(1) = vprime(2)      
!        wprime(1) = wprime(2)      

        return
 
      END SUBROUTINE stoch 

    
