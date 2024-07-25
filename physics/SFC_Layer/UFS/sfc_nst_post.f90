!>  \file sfc_nst_post.f90
!!  This file contains code to be executed after the GFS NSST model.

module sfc_nst_post

  use machine               , only : kind_phys, kp => kind_phys
  use module_nst_water_prop , only : get_dtzm_2d

  implicit none

contains

  ! \defgroup GFS_NSST_POST GFS Near-Surface Sea Temperature Post

  !> \section arg_table_sfc_nst_post_run Argument Table
  !! \htmlinclude sfc_nst_post_run.html
  !!
  ! \section NSST_general_post_algorithm General Algorithm
  !
  ! \section NSST_detailed_post_algorithm Detailed Algorithm
  ! @{
  subroutine sfc_nst_post_run                                    &
       ( im, kdt, rlapse, tgice, wet, oceanfrac, use_lake_model, &
       icy, oro, oro_uf, nstf_name1,                             &
       nstf_name4, nstf_name5, xt, xz, dt_cool, z_c, tref, xlon, &
       xlat, tsurf_wat, tsfc_wat, nthreads, dtzm, errmsg, errflg &
       )
    !  ---  inputs:
    integer, intent(in) :: im, kdt, nthreads
    logical, dimension(:), intent(in) :: wet, icy
    integer, dimension(:), intent(in) :: use_lake_model
    real (kind=kind_phys), intent(in) :: rlapse, tgice
    real (kind=kind_phys), dimension(:), intent(in) :: oro, oro_uf
    integer, intent(in) :: nstf_name1, nstf_name4, nstf_name5
    real (kind=kind_phys), dimension(:), intent(in) :: xlat, xlon, oceanfrac
    real (kind=kind_phys), dimension(:), intent(in), optional :: xt, xz, dt_cool, z_c, tref

    !  ---  input/outputs:
    real (kind=kind_phys), dimension(:), intent(inout) :: tsurf_wat, tsfc_wat

    !  ---  outputs:
    real (kind=kind_phys), dimension(:), intent(out) :: dtzm

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !  ---  locals
    integer :: i
    real(kind=kind_phys) :: zsea1, zsea2

    real(kind=kind_phys) :: frz=273.15, small=.05, testlon, testlat,	&
                            alon, alat
    real,parameter :: rad2deg = 57.2957795
    logical doprint
    doprint(alon,alat)=abs(testlon-alon).lt.small .and.		&
                       abs(testlat-alat).lt.small

    call get_testpt(testlon,testlat)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    !     if (lprnt) print *,' tseaz2=',tseal(ipr),' tref=',tref(ipr),
    !    &     ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*xt(ipr)/xz(ipr),
    !    &     ' kdt=',kdt

    !      do i = 1, im
    !        if (wet(i) .and. .not. icy(i)) then
    !          tsurf_wat(i) = tsurf_wat(i) - (oro(i)-oro_uf(i)) * rlapse
    !        endif
    !      enddo

    !  --- ...  run nsst model  ... ---

    if (nstf_name1 > 1) then
       zsea1 = 0.001_kp*real(nstf_name4)
       zsea2 = 0.001_kp*real(nstf_name5)
       call get_dtzm_2d (xt, xz, dt_cool, z_c, wet, zsea1, zsea2, im, 1, nthreads, dtzm)
       do i = 1, im
          !         if (wet(i) .and. .not.icy(i)) then
          !         if (wet(i) .and. (frac_grid .or. .not. icy(i))) then
          if (wet(i) .and. use_lake_model(i) /=1 .and. oceanfrac(i)==0.)&
                                                                then
             tsfc_wat(i) = max(tgice, tref(i) + dtzm(i))
             !           tsfc_wat(i) = max(271.2, tref(i) + dtzm(i)) -  &
             !                           (oro(i)-oro_uf(i))*rlapse
          endif

          alon=xlon(i)*rad2deg
          alat=xlat(i)*rad2deg
          if (doprint(alon,alat)) then
           print 98,'exiting sfc_nst_post_run, lola=',                   &
             alon,alat,                                                  &
           'tsfc_w',tsfc_wat(i)-frz
          end if
 99       format (/a,2f7.2/(5(a8,"=",f7.2)))
 98       format (/a,2f7.2/(4(a8,"=",es11.4)))

       enddo
    endif

    !     if (lprnt) print *,' tseaz2=',tsea(ipr),' tref=',tref(ipr),   &
    !    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

    return
  end subroutine sfc_nst_post_run

end module sfc_nst_post
