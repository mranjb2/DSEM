module filterStorage
    use size
    use physics
    use domain_definition
    implicit none
    
    
    integer, dimension(3)                   :: hncg !< number of gauss points on high p-order grid
    integer, dimension(3)                   :: lncg !< number of gauss points on low p-order grid
    double precision, dimension(nmax,3)     :: hcxg !< mapped gauss point location on high p-order grid
    double precision, dimension(nmax,3)     :: lcxg !< mapped gauss point location on low p-order grid
    double precision, dimension(nx,ny,nz)   :: fhx  !< filter basis in x
    double precision, dimension(nx,ny,nz)   :: fhy  !< filter basis in y
    double precision, dimension(nx,ny,nz)   :: fhz  !< filter basis in z
    
    
contains
    
    subroutine filterPreCompute(d,phigh,plow,alpha)
        
        implicit none
        type(domain), intent(in)    :: d !< Element storage
        integer                     :: i,j,k,l,m,n
            
        
        ! Compute the low p grids constants
        hncg = d%ncg
        lncg = hncg - (phigh-plow)
        
        ! Set up the low-p gauss points
        lcxg = 0.0d0
        DO k = 1,3
           DO n = 1,lncg(k)
              xx = cos((2.d0*n-1.d0)*pi/(2.d0*(lncg(k) - 1.d0) + 2.d0)) 
              lcxg(n,k) = 0.5d0*(1.d0 - xx)
           END DO
        END DO
        
        ! Compute the basis functions between grids
        hx = 0.0d0
        hy = 0.0d0
        hz = 0.0d0
        
        do n=1,hncg(1)
            hx(n) = polyn(n,)
        end do
        DO m=1,d%ncg(2)
           hy(m)=polyn(m,dr%Xpmap(2),d%ncg(2),d%cxg(1:d%ncg(2),2))
        ENDDO
        DO l=1,d%ncg(3)
           hz(l)=polyn(l,dr%Xpmap(3),d%ncg(3),d%cxg(1:d%ncg(3),3))
        ENDDO
        
    
    end subroutine 
    
    subroutine filter3d(Qold,Qnew,plow,phigh,alpha)
        
        implicit none
        ! variable definitions
        double precision, intent(in)    :: Qold     !< Q variable to be filtered
        double precision, intent(inout) :: Qnew     !< Q variable post filtering
        double precision, intent(in)    :: alpha    !< Filter blending quantity
        integer, intent(in)             :: phigh    !< high polynomial order
        integer, intent(in)             :: plow     !< low polynomial order
        integer                         :: sub
    
        ! begin code
        sub = phigh - plow
    
    end subroutine filter3d

end module filterStorage




