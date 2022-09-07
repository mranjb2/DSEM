!> \file mortarLib.f90
!! All routines pertaining to finite element mortar method
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!////////											////////
!////////	mortarLib.f90									////////
!////////											////////
!////////	contains:									////////
!////////											////////
!////////	   SUBROUTINE set_mortar_matrices(mrtr)						////////
!////////	   LOGICAL FUNCTION copy_available(idir,idirt,mface,which_copy)			////////
!////////	   SUBROUTINE compute_mortar_matrices(mrtr,mface,idir,idirt,prolong,restrict)	////////
!////////	   SUBROUTINE compute_pmat(ltype,pmat,cxin,ni,cxo,no,oset,scale,ndim)		////////
!////////	   SUBROUTINE compute_bmat(bmat,cx,n,ndim)					////////
!////////	   SUBROUTINE compute_resmat(amat,cxin,ni,cxo,no,oset,scale,ndim)		////////
!////////	   SUBROUTINE compute_promat(amat,cxin,ni,cxo,no,oset,scale,ndim)		////////
!////////	   SUBROUTINE mrtrmat(ncol,no,xo,js,je,xnew,b,scale)				////////
!////////	   SUBROUTINE decomp(ndim,n,a,cond,ipvt,work)					////////
!////////	   SUBROUTINE solve(ndim, n, a, b, ipvt)					////////
!////////											////////
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine sets up arrays for mortar projections.
!> To save on storage, only unique prolongation and restriction 
!> arrays are computed and stored in an array.
!!
!
      SUBROUTINE set_mortar_matrices(mrtr) 
!
!......................................................................
!     date: 12/9/98
!     routines called: compute_pmat                                           
!
!     sets up arrays for mortar projections. to save on storage,
!     only unique prolongation and restriction arrays are computed
!     and stored in an array.
!......................................................................
!
      USE mortar_definition
      USE Order_Matrices
      USE File_Units
      USE Edge_Mappings

      IMPLICIT none

      TYPE (mortar)             :: mrtr
      INTEGER                   :: k,listloc,n,m,l,jj
      INTEGER                   :: mface,idir,idirs,idirm
!
!     --------------------
!     side 1 (master) face
!     --------------------
!
      mface = 1
      idir = 1
      idirm = 1
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirm,mface,listloc) )     THEN
            mrtr%prolong1_dir1  => prolong(:,:,listloc)
            mrtr%restrict1_dir1 => restrict(:,:,listloc)
         ELSE
            num_mrtr_matrices = num_mrtr_matrices + 1
            IF ( num_mrtr_matrices > max_mrtr_mat )     THEN
               WRITE(err_unit,*) 'Too many mortar matrices. increase max_mrtr_mat in size'
               STOP 'Too many mortar matrices'
            END IF
            mrtr_mat_map(1,num_mrtr_matrices) = mrtr%len(idir,mface)
            mrtr_mat_map(2,num_mrtr_matrices) = mrtr%lenmortar(idir)
            CALL compute_mortar_matrices(mrtr,mface,idir,idirm, &
                 prolong(:,:,num_mrtr_matrices),restrict(:,:,num_mrtr_matrices)) 
            mrtr%prolong1_dir1  => prolong(:,:,num_mrtr_matrices)
            mrtr%restrict1_dir1 => restrict(:,:,num_mrtr_matrices)
         END IF
      END IF

      idir = 2
      idirm = 2
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirm,mface,listloc) )     THEN
            mrtr%prolong1_dir2  => prolong(:,:,listloc)
            mrtr%restrict1_dir2 => restrict(:,:,listloc)
         ELSE
            num_mrtr_matrices = num_mrtr_matrices + 1
            IF ( num_mrtr_matrices > max_mrtr_mat )     THEN
               WRITE(err_unit,*) 'Too many mortar matrices. increase max_mrtr_mat in size'
               STOP 'Too many mortar matrices'
            END IF
            mrtr_mat_map(1,num_mrtr_matrices) = mrtr%len(idir,mface)
            mrtr_mat_map(2,num_mrtr_matrices) = mrtr%lenmortar(idir)
            CALL compute_mortar_matrices(mrtr,mface,idir,idirm, &
                 prolong(:,:,num_mrtr_matrices),restrict(:,:,num_mrtr_matrices)) 
            mrtr%prolong1_dir2  => prolong(:,:,num_mrtr_matrices)
            mrtr%restrict1_dir2 => restrict(:,:,num_mrtr_matrices)
         END IF
      END IF
!
!     ------------------------------------------------
!     side 2 (slave) face (may be rotated wrt mortar)
!     ------------------------------------------------
!
      mface = 2
      idir = 1
      SELECT CASE (mrtr%orient)
         CASE (DFLT,B1F2,B1B2,F1B2)
            idirs = 1
         CASE (F2F1,B2F1,B2B1,F2B1)
            idirs = 2
      END SELECT
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirs,mface,listloc) )     THEN
            mrtr%prolong2_dir1  => prolong(:,:,listloc)
            mrtr%restrict2_dir1 => restrict(:,:,listloc)
         ELSE
            num_mrtr_matrices = num_mrtr_matrices + 1
            IF ( num_mrtr_matrices > max_mrtr_mat )     THEN
               WRITE(err_unit,*) 'Too many mortar matrices. increase max_mrtr_mat in size'
               STOP 'Too many mortar matrices'
            END IF
            mrtr_mat_map(1,num_mrtr_matrices) = mrtr%len(idirs,mface)
            mrtr_mat_map(2,num_mrtr_matrices) = mrtr%lenmortar(idir)
            CALL compute_mortar_matrices(mrtr,mface,idir,idirs, &
                 prolong(:,:,num_mrtr_matrices),restrict(:,:,num_mrtr_matrices)) 
            mrtr%prolong2_dir1 => prolong(:,:,num_mrtr_matrices)
            mrtr%restrict2_dir1 => restrict(:,:,num_mrtr_matrices)
         END IF
      END IF

      idir = 2
      SELECT CASE (mrtr%orient)
         CASE (DFLT,B1F2,B1B2,F1B2)
            idirs = 2
         CASE (F2F1,B2F1,B2B1,F2B1)
            idirs = 1
      END SELECT
      IF ( .NOT.mrtr%conforming(idir,mface) )     THEN
         IF ( copy_available(idir,idirs,mface,listloc) )     THEN
            mrtr%prolong2_dir2  => prolong(:,:,listloc)
            mrtr%restrict2_dir2 => restrict(:,:,listloc)
         ELSE
            num_mrtr_matrices = num_mrtr_matrices + 1
            IF ( num_mrtr_matrices > max_mrtr_mat )     THEN
               WRITE(err_unit,*) 'Too many mortar matrices. increase max_mrtr_mat in size'
               STOP 'Too many mortar matrices'
            END IF
            mrtr_mat_map(1,num_mrtr_matrices) = mrtr%len(idirs,mface)
            mrtr_mat_map(2,num_mrtr_matrices) = mrtr%lenmortar(idir)
            CALL compute_mortar_matrices(mrtr,mface,idir,idirs, &
                 prolong(:,:,num_mrtr_matrices),restrict(:,:,num_mrtr_matrices)) 
            mrtr%prolong2_dir2 => prolong(:,:,num_mrtr_matrices)
            mrtr%restrict2_dir2 => restrict(:,:,num_mrtr_matrices)
         END IF
      END IF
!
      RETURN
!     --------
      CONTAINS
!     --------

         LOGICAL FUNCTION copy_available(idir,idirt,mface,which_copy)
!
!        Go through the list of arrays and see if we've already computed
!        this one.
!
         INTEGER :: which_copy,k,idir,mface,idirt
!
         copy_available = .false.
         DO k = 1,num_mrtr_matrices
            IF ( mrtr_mat_map(1,k) == mrtr%len(idirt,mface) .AND. &
                 mrtr_mat_map(2,k) == mrtr%lenmortar(idir) )  THEN
               copy_available = .true.
               which_copy = k
               EXIT
            END IF
         END DO
         RETURN
         END FUNCTION copy_available
!
      END SUBROUTINE set_mortar_matrices
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the projection matrices for the mortars,
!> the prolongation matrices, which map the faces onto the 
!> mortars, are just the interpolation matrices.
!!
!! The restriction matrices, which project from the mortars onto the faces
!! are the L2 projection matrices.
      SUBROUTINE compute_mortar_matrices(mrtr,mface,idir,idirt,prolong,restrict) 
!
!......................................................................
!     date: 12/9/98                                                     
!
!     compute the projection matrices for the mortars
!     the prolongation matrices, which map the faces onto the
!     mortars are just the interpolation matrices. The restriction
!     matrices, which project from the mortars onto the faces
!     are the L2 projection matrices
!......................................................................
!
      USE mortar_definition
      USE physics
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      TYPE (mortar)                          :: mrtr
      INTEGER                                :: legendre = 1
      DOUBLE PRECISION, DIMENSION(nmax)      :: work
      DOUBLE PRECISION, DIMENSION(2)         :: endpts
      DOUBLE PRECISION, DIMENSION(nmax,2)    :: cx,z,weight
      DOUBLE PRECISION, DIMENSION(nmax,nmax) :: prolong,restrict
      DOUBLE PRECISION                       :: offset = 0.0d0,scalem = 1.0d0
!
!     ---------------------------------------------
!     compute gauss points and weights along mortar
!     and scale them for the interval [0,1]
!     ---------------------------------------------
!
      CALL mortar_gauss_points(mrtr%lenmortar,z,weight)
      CALL mortar_gauss_points(mrtr%len(:,mface),cx,weight)
!
!     -------------------------------------------------------------
!     Compute the prolongation matrix. Note that we assume that the
!     meshes in 3D are geometrically conforming, so offset and
!     scale are always 0.0d0 and 1.0d0
!     -------------------------------------------------------------
!
      CALL compute_pmat(1,prolong,cx(:,idirt),mrtr%len(idirt,mface),&
                                   z(:,idir),mrtr%lenmortar(idir),&
                                   offset,scalem,nmax)
      CALL compute_pmat(2,restrict,z(:,idir),mrtr%lenmortar(idir),&
                                   cx(:,idirt),mrtr%len(idirt,mface),&
                                   offset,scalem,nmax)
!
      RETURN 
   END SUBROUTINE compute_mortar_matrices
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the projection matrix for order refinement. 
!!
!
   SUBROUTINE compute_pmat(ltype,pmat,cxin,ni,cxo,no,oset,scale,ndim) 
!
!     date: 2/22/95                                                     
!     routines called: compute_bmat                                     
!                      compute_amat                                     
!                      decomp                                           
!                      solve                                            
!     applicability:mortar routines                                     
!
!     Compute the projection matrix for order refinement                
!
      USE size
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
      DIMENSION bmat(nmax,nmax),cxo(*),amat(nmax,nmax),cxin(*) 
      DIMENSION pmat(ndim,ndim) 
      DIMENSION ipvt(nmax),work(nmax) 
 
!
!     --------------------------------------
!     compute "b" matrix for output variable                            
!     --------------------------------------
!
      CALL compute_bmat(bmat,cxo,no,nmax)
!
!     --------------------------------------
!     compute "a" matrix for input variables                            
!     ltype = 1 <=> prolongation                                        
!     ltype = 2 <=> restriction                                         
!     --------------------------------------
!
      IF(ltype.eq.2)     THEN 
         CALL compute_resmat(amat,cxin,ni,cxo,no,oset,scale,nmax) 
      ELSE 
         CALL compute_promat(amat,cxin,ni,cxo,no,oset,scale,nmax) 
      ENDIF 
!
!     -------------------------
!     compute projection matrix
!     -------------------------
!
      CALL decomp(nmax,no,bmat,cond,ipvt,work) 
      DO j = 1,ni 
         CALL solve(nmax,no,bmat,amat(1,j),ipvt) 
      END DO
      eps = 100*epsilon(eps)
      DO j = 1,ni 
          DO i = 1,no 
             pmat(i,j) = amat(i,j)
             IF( ABS(pmat(i,j)) <= eps )     pmat(i,j) = 0.0d0
          END DO 
      END DO 
!
      RETURN 
   END SUBROUTINE compute_pmat
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the b matrix for order refinement.
!!
!
   SUBROUTINE compute_bmat(bmat,cx,n,ndim) 
!
!     date: 2/21/95                                                     
!     routines called: polyn                                            
!     applicability:mortar routines                                     
!
!     compute the b matrix for order refinement                         
!
      USE size
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      DIMENSION w(nmax),cxg(nmax),bmat(ndim,ndim),cx(ndim)
      INTEGER                           :: legendre = 1
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
!
      CALL gaussq(legendre,n,0.0d0,0.0d0,0,endpts,work,cxg,w)
!
      DO i = 1,n 
         DO j = 1,n 
            sum = 0.0d0 
            DO k = 1,n 
               hj = polyn(j,0.5d0*(cxg(k) + 1.d0),n,cx) 
               hi = polyn(i,0.5d0*(cxg(k) + 1.d0),n,cx) 
               sum = sum + 0.5d0*w(k)*hi*hj 
            END DO 
            bmat(i,j) = sum 
         END DO 
      END DO 
!
      RETURN 
   END SUBROUTINE compute_bmat
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the restriction rhs matrix.
!!
!
   SUBROUTINE compute_resmat(amat,cxin,ni,cxo,no,oset,scale,ndim) 
!
!     date: 2/21/95                                                     
!     routines called: polyn                                            
!     applicability:mortar routines                                     
!
!     compute the restriction rhs matrix                                
!
      USE size
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
      DIMENSION w(nmax),cxg(nmax),amat(ndim,ndim),cxo(ndim),cxin(ndim) 
      INTEGER                           :: legendre = 1
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
!
      ncc = max(ni,no)
      CALL gaussq(legendre,ncc,0.0d0,0.0d0,0,endpts,work,cxg,w)
!
      do i = 1,no 
         do j = 1,ni 
            sum = 0.0d0 
            do k = 1,ncc 
               xi = oset + scale*0.5d0*(cxg(k) + 1.d0)
               hin = polyn(j,0.5d0*(cxg(k) + 1.d0),ni,cxin) 
               hout = polyn(i,xi,no,cxo) 
               sum = sum + 0.5d0*w(k)*hin*hout 
            end do 
            amat(i,j) = sum 
         end do 
      end do 
!
      return 
   END SUBROUTINE compute_resmat
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the prolongation rhs matrix.
!!
!
   SUBROUTINE compute_promat(amat,cxin,ni,cxo,no,oset,scale,ndim) 
!
!     date: 2/21/95                                                     
!     routines called: polyn                                            
!     includes: size.inc                                                
!     applicability:mortar routines                                     
!
!     compute the prolongation rhs matrix                               
!
      USE size
      USE constants
!
      implicit double precision (a-h,o-z) 
      dimension w(nmax),cxg(nmax),amat(ndim,ndim),cxo(ndim),cxin(ndim) 
      INTEGER :: legendre = 1
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
!
      ncc = max(ni,no)
      call gaussq(legendre,ncc,0.0d0,0.0d0,0,endpts,work,cxg,w)
!
      do i = 1,no 
         do j = 1,ni 
            sum = 0.0d0 
            do k = 1,ncc 
               xi = oset + scale*0.5d0*(cxg(k) + 1.d0) 
               hin = polyn(j,xi,ni,cxin) 
               hout = polyn(i,0.5d0*(cxg(k) + 1.d0),no,cxo) 
               sum = sum + 0.5d0*w(k)*hin*hout 
            end do 
            amat(i,j) = sum 
         end do 
      end do 
!
      return 
   END SUBROUTINE compute_promat
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the interpolation matrix.
!!
!
   SUBROUTINE mrtrmat(ncol,no,xo,js,je,xnew,b,scale) 
!
!     date: 2/22/95                                                     
!     routines called: polyn                                            
!     includes: none                                                    
!     applicability:mortar versions                                     
!
!
!     compute the interpolation matrix.                                 
!
      implicit double precision (a-h,o-z) 
      dimension xo(*),xnew(*) 
      double precision b(ncol,*),polyn 
!
      s = 1.d0/scale 
      do 100 k = js,je 
         do 100 j = 1,no 
            b(k,j) = s*polyn(j,xnew(k),no,xo) 
  100 continue 
!
      return 
   END SUBROUTINE mrtrmat  
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine decomposes a double precision matrix by gaussian elimination
!> and estimates the condition of the matrix.
!!
!
   SUBROUTINE decomp(ndim,n,a,cond,ipvt,work) 
!
      integer ndim,n 
      double precision a(ndim,n),cond,work(n) 
      integer ipvt(n) 
!
!     decomposes a double precision matrix by gaussian elimination      
!     and estimates the condition of the matrix.                        
!
!     use solve to compute solutions to linear systems.                 
!
!     input..                                                           
!
!        ndim = declared row dimension of the array containing  a.      
!
!        n = order of the matrix.                                       
!
!        a = matrix to be triangularized%                               
!
!     output..                                                          
!
!        a  contains an upper triangular matrix  u  and a permuted      
!          version of a lower triangular matrix  i-l  so that           
!          (permutation matrix)*a = l*u .                               
!
!        cond = an estimate of the condition of  a .                    
!           for the linear system  a*x = b, changes in  a  and  b       
!           may cause changes  cond  times as large in  x .             
!           if  cond+1.0 .eq. cond , a is singular to working           
!           precision.  cond  is set to  1.0d+32  if exact              
!           singularity is detected%                                    
!
!        ipvt = the pivot vector.                                       
!           ipvt(k) = the index of the k-th pivot row                   
!           ipvt(n) = (-1)**(number of interchanges)                    
!
!     work space..  the vector  work  must be declared and included     
!                   in the call.  its input contents are ignored%       
!                   its output contents are usually unimportant.        
!
!     the determinant of a can be obtained on output by                 
!        det(a) = ipvt(n) * a(1,1) * a(2,2) * ... * a(n,n)%             
!
      double precision ek, t, anorm, ynorm, znorm 
      integer nm1, i, j, k, kp1, kb, km1, m 
      double precision dabs, dsign 
!
      ipvt(n) = 1 
      if (n .eq. 1) go to 80 
      nm1 = n - 1 
!
!     compute 1-norm of a                                               
!
      anorm = 0.0d0 
      do 10 j = 1, n 
         t = 0.0d0 
         do 5 i = 1, n 
            t = t + dabs(a(i,j)) 
    5    continue 
         if (t .gt. anorm) anorm = t 
   10 continue 
!
!     gaussian elimination with partial pivoting                        
!
      do 35 k = 1,nm1 
         kp1= k+1 
!
!        find pivot                                                     
!
         m = k 
         do 15 i = kp1,n 
            if (dabs(a(i,k)) .gt. dabs(a(m,k))) m = i 
   15    continue 
         ipvt(k) = m 
         if (m .ne. k) ipvt(n) = -ipvt(n) 
         t = a(m,k) 
         a(m,k) = a(k,k) 
         a(k,k) = t 
!
!        skip step if pivot is zero                                     
!
         if (t .eq. 0.0d0) go to 35 
!
!        compute multipliers                                            
!
         do 20 i = kp1,n 
             a(i,k) = -a(i,k)/t 
   20    continue 
!
!        interchange and eliminate by columns                           
!
         do 30 j = kp1,n 
             t = a(m,j) 
             a(m,j) = a(k,j) 
             a(k,j) = t 
             if (t .eq. 0.0d0) go to 30 
             do 25 i = kp1,n 
                a(i,j) = a(i,j) + a(i,k)*t 
   25        continue 
   30    continue 
   35 continue 
!
!     cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)         
!     estimate obtained by one step of inverse iteration for the        
!     small singular vector.  this involves solving two systems         
!     of equations, (a-transpose)*y = e  and  a*z = y  where  e         
!     is a vector of +1 or -1 chosen to cause growth in y.              
!     estimate = (1-norm of z)/(1-norm of y)                            
!
!     solve (a-transpose)*y = e                                         
!
      do 50 k = 1, n 
         t = 0.0d0 
         if (k .eq. 1) go to 45 
         km1 = k-1 
         do 40 i = 1, km1 
            t = t + a(i,k)*work(i) 
   40    continue 
   45    ek = 1.0d0 
         if (t .lt. 0.0d0) ek = -1.0d0 
         if (a(k,k) .eq. 0.0d0) go to 90 
         work(k) = -(ek + t)/a(k,k) 
   50 continue 
      do 60 kb = 1, nm1 
         k = n - kb 
         t = 0.0d0 
         kp1 = k+1 
         do 55 i = kp1, n 
            t = t + a(i,k)*work(i) 
   55    continue 
         work(k) = t + work(k) 
         m = ipvt(k) 
         if (m .eq. k) go to 60 
         t = work(m) 
         work(m) = work(k) 
         work(k) = t 
   60 continue 
!
      ynorm = 0.0d0 
      do 65 i = 1, n 
         ynorm = ynorm + dabs(work(i)) 
   65 continue 
!
!     solve a*z = y                                                     
!
      call solve(ndim, n, a, work, ipvt) 
!
      znorm = 0.0d0 
      do 70 i = 1, n 
         znorm = znorm + dabs(work(i)) 
   70 continue 
!
!     estimate condition                                                
!
      cond = anorm*znorm/ynorm 
      if (cond .lt. 1.0d0) cond = 1.0d0 
      return 
!
!     1-by-1                                                            
!
   80 cond = 1.0d0 
      if (a(1,1) .ne. 0.0d0) return 
!
!     exact singularity                                                 
!
   90 cond = 1.0d+32 
      return 
   END SUBROUTINE decomp
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine solves for a solution of linear system, a*x = b.
!!
!! Do not use if decomp has detected singularity. 
!
   SUBROUTINE solve(ndim, n, a, b, ipvt) 
!
      integer ndim, n, ipvt(n) 
      double precision a(ndim,n),b(n) 
!
!   solution of linear system, a*x = b .                                
!   do not use if decomp has detected singularity.                      
!
!   input..                                                             
!
!     ndim = declared row dimension of array containing a .             
!
!     n = order of matrix.                                              
!
!     a = triangularized matrix obtained from decomp .                  
!
!     b = right hand side vector.                                       
!
!     ipvt = pivot vector obtained from decomp .                        
!
!   output..                                                            
!
!     b = solution vector, x .                                          
!
      integer kb, km1, nm1, kp1, i, k, m 
      double precision t 
!
!     forward elimination                                               
!
      if (n .eq. 1) go to 50 
      nm1 = n-1 
      do 20 k = 1, nm1 
         kp1 = k+1 
         m = ipvt(k) 
         t = b(m) 
         b(m) = b(k) 
         b(k) = t 
         do 10 i = kp1, n 
             b(i) = b(i) + a(i,k)*t 
   10    continue 
   20 continue 
!
!     back substitution                                                 
!
      do 40 kb = 1,nm1 
         km1 = n-kb 
         k = km1+1 
         b(k) = b(k)/a(k,k) 
         t = -b(k) 
         do 30 i = 1, km1 
             b(i) = b(i) + a(i,k)*t 
   30    continue 
   40 continue 
   50 b(1) = b(1)/a(1,1) 
      return 
   END SUBROUTINE solve
!
!///////////////////////////////////////////////////////////////////////
!> @brief
!> This subroutine computes the b matrix for order refinement
!!
!
   SUBROUTINE compute_bmata(bmat,cx,n,ndim) 
!
!     date: 2/21/95                                                     
!     routines called: polyn                                            
!     applicability:mortar routines                                     
!
!     compute the b matrix for order refinement                         
!
      USE size
      USE constants
!
      IMPLICIT DOUBLE PRECISION (a-h,o-z) 
!
      DIMENSION w(nmax),cxg(nmax),bmat(ndim,ndim),cx(ndim)
      INTEGER                           :: legendre = 1
      DOUBLE PRECISION, DIMENSION(nmax) :: work
      DOUBLE PRECISION, DIMENSION(2)    :: endpts
!
      CALL gaussq(legendre,n,0.0d0,0.0d0,0,endpts,work,cxg,w)
!
      DO i = 1,n 
         DO j = 1,n 
            sum = 0.0d0 
            DO k = 1,n 
               hj = polyn(j,0.5d0*(cxg(k) + 1.d0),n,cx) 
               hi = polyn(i,0.5d0*(cxg(k) + 1.d0),n,cx) 
               sum = sum + 0.5d0*w(k)*hi*hj 
            END DO 
            bmat(i,j) = sum 
         END DO 
      END DO 
!
      RETURN 
   END SUBROUTINE compute_bmata
!
!///////////////////////////////////////////////////////////////////////
!
