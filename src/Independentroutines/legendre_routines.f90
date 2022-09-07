!> \file legendre_routines.f90
!! Contains routines for calculation of Gauss Legendre quadrature rules
!///////////////////////////////////////////////////////////////////////////////////////////////
!////////										////////
!////////	legendre_routines.f90							////////
!////////										////////
!////////	contains:								////////
!////////										////////
!////////	      SUBROUTINE gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)	////////
!////////										////////
!///////////////////////////////////////////////////////////////////////////////////////////////
!
!> @brief This set of routines computes the nodes \f$ t(j) \f$ and weights    
!!       \f$ w(j) \f$ for gaussian-type quadrature rules with pre-assigned      
!!        nodes. these are used when one wishes to approximate
!!
!!                 \f$ \int_a^b \! f(x) w(x) \, dx \f$
!!
!! by \f$ \sum\limits_{j=1}^n w_j f(t_j) \f$
!!
!!    (note \f$ w(x) \f$ and \f$ w(j) \f$ have no connection with each other.)       
!!    here w(x) is one of six possible non-negative weight           
!!    functions (listed below), and \f$ f(x) \f$ is the                      
!!    function to be integrated.  gaussian quadrature is particularly
!!    useful on infinite intervals (with appropriate weight          
!!    functions), since then other techniques often fail.            
!!
!!       associated with each weight function \f$ w(x) \f$ is a set of       
!!    orthogonal polynomials.  the nodes \f$ t(j) \f$ are just the zeroes    
!!    of the proper n-th degree polynomial.                          
!!
!! Input parameters (all real numbers are in double precision):
!!
!!    kind     an integer between 1 and 6 giving the type of         
!!             quadrature rule:                                      
!!
!!    kind = 1:  legendre quadrature, \f$ w(x) = 1 \f$ on \f$ (-1, 1) \f$            
!!
!!    kind = 2:  chebyshev quadrature of the first kind              
!!               \f$ w(x) = \frac{1}{\sqrt{1-x^2}} \f$ on \f$ (-1, +1) \f$                  
!!
!!    kind = 3:  chebyshev quadrature of the second kind             
!!               \f$ w(x) = \sqrt{1-x^2} \f$ on \f$ (-1, 1) \f$                     
!!
!!    kind = 4:  hermite quadrature, \f$ w(x) = \exp(-x^2) \f$ on             
!!               (- \infty, + \infty)                              
!!
!!    n:        the number of points used for the quadrature rule     
!!
!!    alpha:    real parameter used only for gauss-jacobi and gauss-  
!!              laguerre quadrature (otherwise use 0.d0).             
!!
!!    beta:     real parameter used only for gauss-jacobi quadrature--
!!              (otherwise use 0.d0)                                  
!!
!!    kpts:     (integer) normally 0, unless the left or right end-   
!!              point (or both) of the interval is required to be a   
!!              node (this is called gauss-radau or gauss-lobatto     
!!              quadrature).  then kpts is the number of fixed        
!!              endpoints (1 or 2).                                   
!!
!!    endpts:   real array of length 2.  contains the values of       
!!              any fixed endpoints, if kpts = 1 or 2.                
!!
!!    b:        real scratch array of length n                        
!!
!! Output parameters (both double precision arrays of length n):
!!
!!    t:        will contain the desired nodes.                       
!!    w:        will contain the desired weights \f$ w(j) \f$.                
!!
!! Underflow may sometimes occur, but is harmless.
!!
!! References:
!!                                                        
!!    1.  golub, g. h., and welsch, j. h., "calculation of gaussian  
!!        quadrature rules," mathematics of computation 23 (april,   
!!        1969), pp. 221-230.                                        
!!
!!    2.  golub, g. h., "some modified matrix eigenvalue problems,"  
!!        siam review 15 (april, 1973), pp. 318-334 (section 7).     
!!
!!    3.  stroud and secrest, gaussian quadrature formulas, prentice-
!!        hall, englewood cliffs, n.j., 1966.                        
!!
!!    original version 20 jan 1975 from stanford.                     
!!    modified 21 dec 1983 by eric grosse.
!!                            
!!      imtql2 => gausq2                                             
!!
!!      hex constant => d1mach (from core library)                   
!!
!!      compute pi using datan                                       
!!
!!      removed accuracy claims, description of method               
!!
!!      added single precision version
      SUBROUTINE gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w) 
!                               
!
      double precision b(n), t(n), w(n), endpts(2), muzero, t1,         &
     & gam, solveg, dsqrt, alpha, beta                                   
!
      call class (kind, n, alpha, beta, b, t, muzero) 
!
!           the matrix of coefficients is assumed to be symmetric.      
!           the array t contains the diagonal elements, the array       
!           b the off-diagonal elements.                                
!           make appropriate changes in the lower right 2 by 2          
!           submatrix.                                                  
!
      if (kpts.eq.0)  go to 100 
      if (kpts.eq.2)  go to  50 
!
!           if kpts=1, only t(n) must be changed                        
!
      t(n) = solveg(endpts(1), n, t, b)*b(n-1)**2 + endpts(1) 
      go to 100 
!
!           if kpts=2, t(n) and b(n-1) must be recomputed               
!
   50 gam = solveg(endpts(1), n, t, b) 
      t1 = ((endpts(1) - endpts(2))/(solveg(endpts(2), n, t, b) - gam)) 
      b(n-1) = dsqrt(t1) 
      t(n) = endpts(1) + gam*t1 
!
!           note that the indices of the elements of b run from 1 to n-1
!           and thus the value of b(n) is arbitrary.                    
!           now compute the eigenvalues of the symmetric tridiagonal    
!           matrix, which has been modified as necessary.               
!           the method used is a ql-type method with origin shifting    
!
  100 w(1) = 1.0d0 
      do 105 i = 2, n 
  105    w(i) = 0.0d0 
!
      call gausq2 (n, t, b, w, ierr) 
      do 110 i = 1, n 
  110    w(i) = muzero * w(i) * w(i) 
!
      return 
      END                                           
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief This procedure performs elimination to solve for the            
!! n-th component of the solution delta to the equation            
!!
!!             (jn - shift*identity) * delta  = en,                      
!!
!!       where en is the vector of all zeroes except for 1 in            
!!       the n-th position.                                              
!!
!!       the matrix jn is symmetric tridiagonal, with diagonal           
!!       elements \f$ a(i) \f$, off-diagonal elements \f$ b(i) \f$.  this equation       
!!       must be solved to obtain the appropriate changes in the lower   
!!       2 by 2 submatrix of coefficients for orthogonal polynomials.
      double precision function solveg(shift, n, a, b)     
!
!
      double precision shift, a(n), b(n), alpha 
!
      alpha = a(1) - shift 
      nm1 = n - 1 
      do 10 i = 2, nm1 
   10    alpha = a(i) - shift - b(i-1)**2/alpha 
      solveg = 1.0d0/alpha 
      return 
      END                                           
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief This procedure supplies the coefficients \f$ a(j) \f$ , \f$ b(j) \f$ of the  
!!        recurrence relation                                            
!!
!!        \f$ b_j p_j(x) = (x-a_j) p_{j-1} (x) - b_{j-1} p_{j-2} (x) \f$
!!
!! for the various classical (normalized) orthogonal polynomials, 
!! and the zero-th moment                                         
!!
!!             \f$ muzero = \int \! w(x) dx \f$
!!
!! of the given polynomial's weight function \f$ w(x) \f$ .  since the     
!! polynomials are orthonormalized, the tridiagonal matrix is     
!! guaranteed to be symmetric.                                    
!!
!!    the input parameter alpha is used only for laguerre and     
!! jacobi polynomials, and the parameter beta is used only for    
!! jacobi polynomials.  the laguerre and jacobi polynomials       
!! require the gamma function. (DISABLED)
      subroutine class(kind, n, alpha, beta, b, a, muzero) 
!                         
!
      double precision a(n), b(n), muzero, alpha, beta 
      double precision abi, a2b2, pi, dsqrt, ab 
!
      pi = 4.0d0 * datan(1.0d0) 
      nm1 = n - 1 
      go to (10, 20, 30, 40), kind 
!
!              kind = 1:  legendre polynomials p(x)                     
!              on (-1, +1), w(x) = 1.                                   
!
   10 muzero = 2.0d0 
      do 11 i = 1, nm1 
         a(i) = 0.0d0 
         abi = i 
   11    b(i) = abi/dsqrt(4*abi*abi - 1.0d0) 
      a(n) = 0.0d0 
      return 
!
!              kind = 2:  chebyshev polynomials of the first kind t(x)  
!              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)                    
!
   20 muzero = pi 
      do 21 i = 1, nm1 
         a(i) = 0.0d0 
   21    b(i) = 0.5d0 
      b(1) = dsqrt(0.5d0) 
      a(n) = 0.0d0 
      return 
!
!              kind = 3:  chebyshev polynomials of the second kind u(x) 
!              on (-1, +1), w(x) = sqrt(1 - x*x)                        
!
   30 muzero = pi/2.0d0 
      do 31 i = 1, nm1 
         a(i) = 0.0d0 
   31    b(i) = 0.5d0 
      a(n) = 0.0d0 
      return 
!
!              kind = 4:  hermite polynomials h(x) on (-infinity,       
!              +infinity), w(x) = exp(-x**2)                            
!
   40 muzero = dsqrt(pi) 
      do 41 i = 1, nm1 
         a(i) = 0.0d0 
   41    b(i) = dsqrt(i/2.0d0) 
      a(n) = 0.0d0 
      return 
!
      return 
      END                                           
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief This subroutine is a translation of an algol procedure,           
!!     num. math. 12, 377-383(1968) by martin and wilkinson,             
!!     as modified in num. math. 15, 450(1970) by dubrulle.              
!!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).   
!!     this is a modified version of the 'eispack' routine imtql2.       
!!
!!     this subroutine finds the eigenvalues and first components of the 
!!     eigenvectors of a symmetric tridiagonal matrix by the implicit ql 
!!     method.                                                           
!!
!!     on input:                                                         
!!
!!        n is the order of the matrix;                                  
!!
!!        d contains the diagonal elements of the input matrix;          
!!
!!        e contains the subdiagonal elements of the input matrix        
!!          in its first n-1 positions.  e(n) is arbitrary;              
!!
!!        z contains the first row of the identity matrix.               
!!
!!      on output:                                                       
!!
!!        d contains the eigenvalues in ascending order.  if an          
!!          error exit is made, the eigenvalues are correct but          
!!          unordered for indices 1, 2, ..., ierr-1;                     
!!
!!        e has been destroyed;                                          
!!
!!        z contains the first components of the orthonormal eigenvectors
!!          of the symmetric tridiagonal matrix.  if an error exit is    
!!          made, z contains the eigenvectors associated with the stored 
!!          eigenvalues;                                                 
!!
!!        ierr is set to                                                 
!!          zero       for normal return,                                
!!          j          if the j-th eigenvalue has not been               
!!                     determined after 30 iterations. 
      subroutine gausq2(n, d, e, z, ierr) 
!
!     ------------------------------------------------------------------
!
      integer i, j, k, l, m, n, ii, mml, ierr 
      double precision d(n), e(n), z(n), b, c, f, g, p, r, s, machep 
      double precision dsqrt, dabs, dsign 
!
      machep = epsilon(1.d0) 
!
      ierr = 0 
      if (n .eq. 1) go to 1001 
!
      e(n) = 0.0d0 
      do 240 l = 1, n 
         j = 0 
!     :::::::::: look for small sub-diagonal element ::::::::::         
  105    do 110 m = l, n 
            if (m .eq. n) go to 120 
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1))))   &
     &         go to 120                                                
  110    continue 
!
  120    p = d(l) 
         if (m .eq. l) go to 240 
         if (j .eq. 30) go to 1000 
         j = j + 1 
!     :::::::::: form shift ::::::::::                                  
         g = (d(l+1) - p) / (2.0d0 * e(l)) 
         r = dsqrt(g*g+1.0d0) 
         g = d(m) - p + e(l) / (g + dsign(r, g)) 
         s = 1.0d0 
         c = 1.0d0 
         p = 0.0d0 
         mml = m - l 
!
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::             
         do 200 ii = 1, mml 
            i = m - ii 
            f = s * e(i) 
            b = c * e(i) 
            if (dabs(f) .lt. dabs(g)) go to 150 
            c = g / f 
            r = dsqrt(c*c+1.0d0) 
            e(i+1) = f * r 
            s = 1.0d0 / r 
            c = c * s 
            go to 160 
  150       s = f / g 
            r = dsqrt(s*s+1.0d0) 
            e(i+1) = g * r 
            c = 1.0d0 / r 
            s = s * c 
  160       g = d(i+1) - p 
            r = (d(i) - g) * s + 2.0d0 * c * b 
            p = s * r 
            d(i+1) = g + p 
            g = c * r - b 
!     :::::::::: form first component of vector ::::::::::              
            f = z(i+1) 
            z(i+1) = s * z(i) + c * f 
  200       z(i) = c * z(i) - s * f 
!
         d(l) = d(l) - p 
         e(l) = g 
         e(m) = 0.0d0 
         go to 105 
  240 continue 
!
!     :::::::::: order eigenvalues and eigenvectors ::::::::::          
      do 300 ii = 2, n 
         i = ii - 1 
         k = i 
         p = d(i) 
!
         do 260 j = ii, n 
            if (d(j) .ge. p) go to 260 
            k = j 
            p = d(j) 
  260    continue 
!
         if (k .eq. i) go to 300 
         d(k) = d(i) 
         d(i) = p 
         p = z(i) 
         z(i) = z(k) 
         z(k) = p 
  300 continue 
!
      go to 1001 
!     :::::::::: set error -- no convergence to an                      
!                eigenvalue after 30 iterations ::::::::::              
 1000 ierr = l 
 1001 return 
!     :::::::::: last card of gausq2 ::::::::::                         
      END                                           
