
module cm1_init3d_nc
  implicit none
  private
  public :: read_init3d_nc
contains

subroutine disp_err(where, status)
  character(len=*), intent(in) :: where
  integer,          intent(in) :: status
end subroutine disp_err
            
subroutine read_init3d_nc( fname, &
     thsnd3, qvsnd3, usnd3, vsnd3, &
     p_sfc2d, th_sfc2d, qv_sfc2d, &
     zsnd, xsnd, ysnd, ierr, nxsnd, nysnd, nzsnd )

  use netcdf
  implicit none
  !----- interface -----
  character(len=*), intent(in)   :: fname
  real, allocatable, intent(out) :: thsnd3(:,:,:), qvsnd3(:,:,:)
  real, allocatable, intent(out) :: usnd3(:,:,:),  vsnd3(:,:,:)
  real, allocatable, intent(out) :: p_sfc2d(:,:), th_sfc2d(:,:), qv_sfc2d(:,:)
  real, allocatable, intent(out) :: zsnd(:), xsnd(:), ysnd(:)
  real, allocatable :: theta_f(:,:,:), qv_f(:,:,:), u_f(:,:,:), v_f(:,:,:)
  integer,          intent(out)  :: ierr, nxsnd, nysnd, nzsnd

  !----- locals -----
  integer :: ncid, status
  integer :: dimid_xh, dimid_yh, dimid_zh
  integer :: varid
  integer :: nx, ny, nz
  character(len=NF90_MAX_NAME) :: dimname
  integer :: dlen

  ierr = 0

  !---------------------------------------------
  ! open
  status = nf90_open(trim(fname), NF90_NOWRITE, ncid)
  if (status /= nf90_noerr) then
    call disp_err('open('//trim(fname)//')', status); ierr=1; return
  end if

  !---------------------------------------------
  ! dims
  status = nf90_inq_dimid(ncid, "xh", dimid_xh)
  if (status /= nf90_noerr) then; call disp_err('inq_dimid(xh)', status); ierr=2; goto 900; end if
  status = nf90_inq_dimid(ncid, "yh", dimid_yh)
  if (status /= nf90_noerr) then; call disp_err('inq_dimid(yh)', status); ierr=3; goto 900; end if
  status = nf90_inq_dimid(ncid, "zh", dimid_zh)
  if (status /= nf90_noerr) then; call disp_err('inq_dimid(zh)', status); ierr=4; goto 900; end if

  status = nf90_inquire_dimension(ncid, dimid_xh, name=dimname, len=nx)
  if (status /= nf90_noerr) then; call disp_err('inq_dim(xh)', status); ierr=5; goto 900; end if
  status = nf90_inquire_dimension(ncid, dimid_yh, name=dimname, len=ny)
  if (status /= nf90_noerr) then; call disp_err('inq_dim(yh)', status); ierr=6; goto 900; end if
  status = nf90_inquire_dimension(ncid, dimid_zh, name=dimname, len=nz)
  if (status /= nf90_noerr) then; call disp_err('inq_dim(zh)', status); ierr=7; goto 900; end if

  nxsnd = nx; nysnd = ny; nzsnd = nz

  !---------------------------------------------
  ! coords: xh,yh,zh are double in your file
  allocate(xsnd(nx), ysnd(ny), zsnd(nz))

  status = nf90_inq_varid(ncid, "xh", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(xh)', status); ierr=8; goto 900; end if
  status = nf90_get_var(ncid, varid, xsnd)
  if (status /= nf90_noerr) then; call disp_err('get_var(xh)', status); ierr=9; goto 900; end if

  status = nf90_inq_varid(ncid, "yh", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(yh)', status); ierr=10; goto 900; end if
  status = nf90_get_var(ncid, varid, ysnd)
  if (status /= nf90_noerr) then; call disp_err('get_var(yh)', status); ierr=11; goto 900; end if

  status = nf90_inq_varid(ncid, "zh", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(zh)', status); ierr=12; goto 900; end if
  status = nf90_get_var(ncid, varid, zsnd)
  if (status /= nf90_noerr) then; call disp_err('get_var(zh)', status); ierr=13; goto 900; end if

  !---------------------------------------------
  ! 2-D surfaces (yh,xh)
  allocate(p_sfc2d(nx,ny), th_sfc2d(nx,ny), qv_sfc2d(nx,ny))  ! match your Fortran storage (i=x, j=y)

  status = nf90_inq_varid(ncid, "psfc", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(psfc)', status); ierr=14; goto 900; end if
  status = nf90_get_var(ncid, varid, p_sfc2d)
  if (status /= nf90_noerr) then; call disp_err('get_var(psfc)', status); ierr=15; goto 900; end if

  status = nf90_inq_varid(ncid, "th0", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(th0)', status); ierr=16; goto 900; end if
  status = nf90_get_var(ncid, varid, th_sfc2d)
  if (status /= nf90_noerr) then; call disp_err('get_var(th0)', status); ierr=17; goto 900; end if
  
  status = nf90_inq_varid(ncid, "qv0", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(qv0)', status); ierr=18; goto 900; end if
  status = nf90_get_var(ncid, varid, qv_sfc2d)
  if (status /= nf90_noerr) then; call disp_err('get_var(qv0)', status); ierr=19; goto 900; end if
  
  !---------------------------------------------
  ! 3-D fields are stored as (zh,yh,xh) in file; we keep them that way for now, then reorder
  ! File read buffers in file order:
  
  allocate(thsnd3(nx,ny,nz), qvsnd3(nx,ny,nz), usnd3(nx,ny,nz), vsnd3(nx,ny,nz))

  status = nf90_inq_varid(ncid, "theta", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(theta)', status); ierr=20; goto 900; end if
  status = nf90_get_var(ncid, varid, thsnd3)
  if (status /= nf90_noerr) then; call disp_err('get_var(theta)', status); ierr=21; goto 900; end if
  
  !write(*,'(A,1X,ES12.5,1X,ES12.5)') 'theta_f min/max =', minval(thsnd3), maxval(thsnd3)
  
  status = nf90_inq_varid(ncid, "qv", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(qv)', status); ierr=22; goto 900; end if
  status = nf90_get_var(ncid, varid, qvsnd3)
  if (status /= nf90_noerr) then; call disp_err('get_var(qv)', status); ierr=23; goto 900; end if

  !write(*,'(A,1X,ES12.5,1X,ES12.5)') 'qv_f min/max =', minval(qv_f(:,:,1)), maxval(qv_f(:,:,1))

  status = nf90_inq_varid(ncid, "u", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(u)', status); ierr=24; goto 900; end if
  status = nf90_get_var(ncid, varid, usnd3)
  if (status /= nf90_noerr) then; call disp_err('get_var(u)', status); ierr=25; goto 900; end if
  
  status = nf90_inq_varid(ncid, "v", varid)
  if (status /= nf90_noerr) then; call disp_err('inq_varid(v)', status); ierr=26; goto 900; end if
  status = nf90_get_var(ncid, varid, vsnd3)
  if (status /= nf90_noerr) then; call disp_err('get_var(v)', status); ierr=27; goto 900; end if

  ! close
900 continue
  status = nf90_close(ncid)
  if (status /= nf90_noerr) then
    call disp_err('close('//trim(fname)//')', status)
    if (ierr == 0) ierr = -99
  end if


end subroutine read_init3d_nc

end module cm1_init3d_nc
