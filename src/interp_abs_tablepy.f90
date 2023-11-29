

! real(4) function interp_abs_table(abs,q_values,p_values,t_values,q_interp,p_interp,t_interp)

!     real(kind=4), intent(in) :: abs(:,:,:)
!     real(kind=4), intent(in) :: q_values(:)
!     real(kind=4), intent(in) :: p_values(:)
!     real(kind=4), intent(in) :: t_values(:)
!     real(kind=4), intent(in) :: q_interp
!     real(kind=4), intent(in) :: p_interp
!     real(kind=4), intent(in) :: t_interp

!     interp_abs_table = -9999.0
!     if (q_interp < q_values(1) .or. q_interp > q_values(size(q_values)) .or. &
!         p_interp < p_values(1) .or. p_interp > p_values(size(p_values)) .or. &
!         t_interp < t_values(1) .or. t_interp > t_values(size(t_values))) then
!         return
!     end if

!     call trilinear_interpolation(abs, q_values, p_values, t_values, &
!                                  (/q_interp, p_interp, t_interp/), &
!                                  interp_abs_table)

! end function interp_abs_table




    !abs_table = np.zeros((num_q+1,num_p+1,num_T+1),dtype=np.float32,order='F')
subroutine trilinear_interpolation(array, xvalues, yvalues, zvalues, interp_point, result)
    real(kind=4), intent(in) :: array(:,:,:)
    real(kind=4), intent(in) :: xvalues(:)
    real(kind=4), intent(in) :: yvalues(:)
    real(kind=4), intent(in) :: zvalues(:)
    real(kind=4), intent(in) :: interp_point(3)
    real(kind=4), intent(out) :: result

    integer :: i, j, k
    real(kind=8) :: xd, yd, zd
    real(kind=8) :: c00, c01, c10, c11, c0, c1

    ! Find the indices of the nearest lower grid points
    i = 1
    do while (i < size(xvalues) - 1 .and. interp_point(1) > xvalues(i + 1))
        i = i + 1
    end do

    j = 1
    do while (j < size(yvalues) - 1 .and. interp_point(2) > yvalues(j + 1))
        j = j + 1
    end do

    k = 1
    do while (k < size(zvalues) - 1 .and. interp_point(3) > zvalues(k + 1))
        k = k + 1
    end do

    ! Compute fractional distances
    xd = (interp_point(1) - xvalues(i)) / (xvalues(i + 1) - xvalues(i))
    yd = (interp_point(2) - yvalues(j)) / (yvalues(j + 1) - yvalues(j))
    zd = (interp_point(3) - zvalues(k)) / (zvalues(k + 1) - zvalues(k))

    ! Trilinear interpolation with scaling and offset
    c00 = array(i,j,k) * (1.0d0 - xd) + array(i+1,j,k) * xd
    c01 = array(i,j+1,k) * (1.0d0 - xd) + array(i+1,j+1,k) * xd
    c10 = array(i,j,k+1) * (1.0d0 - xd) + array(i+1,j,k+1) * xd
    c11 = array(i,j+1,k+1) * (1.0d0 - xd) + array(i+1,j+1,k+1) * xd

    c0 = c00 * (1.0d0 - yd) + c01 * yd
    c1 = c10 * (1.0d0 - yd) + c11 * yd

    result = c0 * (1.0d0 - zd) + c1 * zd
end subroutine trilinear_interpolation
