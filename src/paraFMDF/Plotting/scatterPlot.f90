!> \file
!! Consistency plotting routines

!> Creates a scatter plot comparing the carrier phase value to the FMDF calculated scalar quantity
subroutine scatterPlot(d,tk)
    use size
    use domain_definition
    use physics
    implicit none
    type(domain),dimension(ngp)                             :: d        !< domain
    double precision                                        :: tempT, rho, p, u, v, w, temp
    integer,    intent(in)                                  :: tk
    integer                                                 :: n,m,l,id
    character(len=32)                                       :: filename
    
    !generate a filename
    write(filename, fmt='(a7,i8.8,a4)') 'scatter',tk,'.csv'
    open(unit=81,file=filename,status="new")
    write(81,*)'Carrier Temp, FMDF Temp'
    do id=1,ngp
        !find average temperature in element
        do l=1,d(id)%ncg(3)
            do m=1,d(id)%ncg(2)
                do n=1,d(id)%ncg(1)
                    rho  = d(id)%Q(n,m,l,1)
                    u    = d(id)%Q(n,m,l,2)/rho
                    v    = d(id)%Q(n,m,l,3)/rho
                    w    = d(id)%Q(n,m,l,4)/rho
                    p    = (gamma-1.d0)*(d(id)%Q(n,m,l,5) - 0.5d0*rho*(u*u + v*v + w*w))
                    temp = p*gamma*mach*mach/rho
                    write(81,*)temp,',',d(id)%s(n,m,l,1)
                end do
            end do
        end do
        
    end do
    close(81)
    
end subroutine scatterPlot