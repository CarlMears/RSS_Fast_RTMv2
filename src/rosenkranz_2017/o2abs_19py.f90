
subroutine abs_o2_rosen_2017(n,t,p,pv,f,o2abs_out_rosen)

    use o2abs_19, only: o2abs

    implicit none

    integer(4) :: n
    real(4), dimension(n), intent(in) :: t,p,pv,f
    real(4), dimension(n) :: o2abs_out_rosen

    integer(4) :: i
    real(4) :: vapor_density
    real(4), parameter :: rd=287.05,epsilon=0.622,rv=rd/epsilon

    do i=1,n   
        ! convert from pv to vapor density in g/m**3
        ! which obabs expects
        vapor_density = (1.0e5/rv)*pv(i)/t(i)

        o2abs_out_rosen(i) = o2abs(t(i),p(i),vapor_density,f(i))
        !print *, 'i,t,p,pv,dens,f,o2abs_out',i, t(i),p(i),pv(i),vapor_density,f(i),o2abs_out_rosen(i)
    enddo

end subroutine abs_o2_rosen_2017



