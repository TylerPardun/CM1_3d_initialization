
module cm1_interp3d
  implicit none
  private
  public :: bilinear_resample2d, vert_linear_to_zh
contains

  pure subroutine find_bracket_and_weight(coord, val, i0, w)
    ! coord(:) must be monotonic increasing
    real, intent(in)  :: coord(:), val
    integer, intent(out) :: i0     ! left index so that coord(i0) <= val <= coord(i0+1)
    real,    intent(out) :: w      ! weight in [0,1] toward i0+1
    integer :: n, lo, hi, mid

    n = size(coord)
    if (val <= coord(1)) then
      i0 = 1; w = 0.0
      return
    else if (val >= coord(n)) then
      i0 = n-1; if (i0 < 1) i0 = 1
      w = 1.0
      return
    end if

    ! binary search
    lo = 1; hi = n
    do
      if (hi - lo <= 1) exit
      mid = (lo + hi) / 2
      if (coord(mid) <= val) then
        lo = mid
      else
        hi = mid
      end if
    end do
    i0 = lo
    if (coord(i0+1) > coord(i0)) then
      w = (val - coord(i0)) / (coord(i0+1) - coord(i0))
    else
      w = 0.0
    end if
  end subroutine find_bracket_and_weight

  pure logical function is_nan(x)
    real, intent(in) :: x
    is_nan = (x /= x)
  end function is_nan

  subroutine bilinear_resample2d(xsrc, ysrc, gsrc, xtgt, ytgt, gdst)
    ! gsrc: (nx_s, ny_s) → gdst: (nx_t, ny_t)
    real, intent(in)  :: xsrc(:), ysrc(:), gsrc(:,:)
    real, intent(in)  :: xtgt(:), ytgt(:)
    real, intent(out) :: gdst(:,:)
    integer :: i,j, ix, iy
    real :: wx, wy
    integer :: nx_s, ny_s

    nx_s = size(xsrc); ny_s = size(ysrc)
    if (nx_s < 2 .or. ny_s < 2) then
      gdst = 0.0; return
    end if

    do j=1,size(ytgt)
      call find_bracket_and_weight(ysrc, ytgt(j), iy, wy)
      if (iy < 1) iy = 1
      if (iy > ny_s-1) iy = ny_s-1
      do i=1,size(xtgt)
        call find_bracket_and_weight(xsrc, xtgt(i), ix, wx)
        if (ix < 1) ix = 1
        if (ix > nx_s-1) ix = nx_s-1

        ! corners
        gdst(i,j) = (1-wx)*(1-wy)*gsrc(ix,  iy  ) + wx*(1-wy)*gsrc(ix+1,iy  ) + &
                    (1-wx)*wy    *gsrc(ix,  iy+1) + wx*wy    *gsrc(ix+1,iy+1)
        if (gdst(i,j) /= gdst(i,j)) gdst(i,j) = 0.0  ! NaN guard
      end do
    end do
  end subroutine bilinear_resample2d
  
subroutine vert_linear_to_zh(zsrc, prof_src, ztgt_col, prof_tgt)
  real, intent(in)  :: zsrc(:), prof_src(:)          ! length nzsnd
  real, intent(in)  :: ztgt_col(:)                   ! length nk (k varies with (i,j))
  real, intent(out) :: prof_tgt(:)                   ! length nk
  integer :: nzs, nk, k, kk, kdn, kup
  real :: w

  nzs = size(zsrc); nk = size(ztgt_col)
  do k=1, nk
    kk=2
    do while (kk<=nzs .and. zsrc(kk) < ztgt_col(k)); kk=kk+1; end do
    kdn = max(1, min(nzs-1, kk-1)); kup = kdn+1
    w = (ztgt_col(k) - zsrc(kdn)) / max(1.e-9, (zsrc(kup)-zsrc(kdn)))
    prof_tgt(k) = prof_src(kdn) + w*(prof_src(kup)-prof_src(kdn))
  end do
end subroutine

end module cm1_interp3d
