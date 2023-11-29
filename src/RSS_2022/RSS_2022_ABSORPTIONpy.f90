
subroutine abs_o2_rss_2022(n,t,p,pv,f,o2abs_rss_out)

    use rss_2022_absorption, only: absoxy_rss_2022_func

    implicit none

    integer(4), intent(in) :: n
    real(4), dimension(n), intent(in) :: t,p,pv,f
    real(4), dimension(n), intent(inout) :: o2abs_rss_out

    integer(4) :: i
    print *, 'n = ',n
    do i=1,n
        o2abs_rss_out(i) = absoxy_rss_2022_func(p(i),t(i),pv(i),f(i))
        !print *, 'i,t,p,pv,f,o2abs_rss_out',i, t(i),p(i),pv(i),f(i),o2abs_rss_out(i)
    enddo

end subroutine abs_o2_rss_2022

subroutine abs_h2o_rss_2022(n,t,p,pv,f,h2oabs_rss_out)

    use rss_2022_absorption, only: absh2o_rss_2022

    implicit none

    integer(4), intent(in) :: n
    real(4), dimension(n), intent(in) :: t,p,pv,f
    real(4), dimension(n), intent(inout) :: h2oabs_rss_out

    integer(4) :: i
    real(4) :: temp
    print *, 'n = ',n
    do i=1,n
        call absh2o_rss_2022(p(i),t(i),pv(i),f(i),  temp)
        h2oabs_rss_out(i) = temp
    enddo

end subroutine abs_h2o_rss_2022

subroutine goff_gratch_rss_2022(n,t,p,relhum,pv)

    use goff_gratch, only: goff_gratch_vap 

    implicit none

    integer(4), intent(in) :: n
    real(4), dimension(n), intent(in) :: t,p,relhum
    real(4), dimension(n), intent(inout) :: pv

    integer(4) :: i
    real(4) :: P_V,RHO_V

    print *, 'goff gratch rss 2022, n = ',n 

    do i=1,n
        call goff_gratch_vap(t(i),p(i),relhum(i),P_V,RHO_V)
        pv(i) = P_V
    enddo
end subroutine goff_gratch_rss_2022





