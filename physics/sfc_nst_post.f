!>  \file sfc_nst_post.f
!!  This file contains code to be executed after the GFS NSST model.

      module sfc_nst_post

      contains

! \defgroup GFS_NSST_POST GFS Near-Surface Sea Temperature Post

!> \section arg_table_sfc_nst_post_run Argument Table
!! \htmlinclude sfc_nst_post_run.html
!!
! \section NSST_general_post_algorithm General Algorithm
!
! \section NSST_detailed_post_algorithm Detailed Algorithm
! @{
      subroutine sfc_nst_post_run                                       &
     &     ( im, kdt, rlapse, tgice, wet, use_flake, icy, oro, oro_uf,  &
     &       nstf_name1,                                                &
     &       nstf_name4, nstf_name5, xt, xz, dt_cool, z_c, tref, xlon,  &
     &       xlat, tsfc_wat, nthreads, dtzm, errmsg, errflg             &
     &       )

      use machine , only : kind_phys
      use module_nst_water_prop, only: get_dtzm_2d
      use module_nst_parameters, only: rad2deg

      implicit none

      integer, parameter :: kp = kind_phys

!  ---  inputs:
      integer, intent(in) :: im, kdt, nthreads
      logical, dimension(:), intent(in) :: wet, icy, use_flake
      real (kind=kind_phys), intent(in) :: rlapse, tgice
      real (kind=kind_phys), dimension(:), intent(in) :: oro, oro_uf
      integer, intent(in) :: nstf_name1, nstf_name4, nstf_name5
      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &      dt_cool, z_c, tref, xlon, xlat

!  ---  input/outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: xt, xz

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(out) :: dtzm, tsfc_wat

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals
      integer :: i
      real (kind=kind_phys), dimension(im) :: dt_warm
      real(kind=kind_phys) :: zsea1, zsea2, alon, alat
      real(kind=kind_phys) :: frz=273.15, small=.05, testlon, testlat
      common /testpt/ testlon,testlat           ! (values defined in dcyc2t3.f)
      logical doprint
      doprint(alon,alat)=abs(testlon-alon).lt.small .and.
     &                   abs(testlat-alat).lt.small


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
        call get_dtzm_2d (xt, xz, dt_cool, z_c, wet, zsea1, zsea2,      &
     &                    im, 1, nthreads, dtzm)
        do i = 1, im
!         if (wet(i) .and. .not.icy(i)) then
!         if (wet(i) .and. (frac_grid .or. .not. icy(i))) then
          if (wet(i) .and. .not. use_flake(i)) then
            tsfc_wat(i) = max(tgice, tref(i) + dtzm(i))
!           tsfc_wat(i) = max(271.2, tref(i) + dtzm(i)) -  &
!                           (oro(i)-oro_uf(i))*rlapse
            dt_warm(i) = 2.*xt(i)/xz(i)
          endif
        enddo
      endif

      do i = 1, im
        alon=xlon(i)*rad2deg
        alat=xlat(i)*rad2deg
        if (doprint(alon,alat)) then
         print 98,'exiting sfc_nst_post   lon,lat=',alon,alat,
     &    'xt',xt(i),
     &    'xz',xz(i),
     &    'dt_warm',dt_warm(i),
     &    'dt_cool',dt_cool(i),
     &    'z_c',z_c(i),
     &    'dtzm',dtzm(i),
     &    'tref',tref(i)-frz,
     &    'tsfc_wa',tsfc_wat(i)
 99      format (/a,2f7.2/(5(a8,"=",f7.2)))
 98      format (/a,2f7.2/(4(a8,"=",es11.4)))
         print '(a,3i2)','nstf_name1/name4/name5=',
     &          nstf_name1,nstf_name4,nstf_name5
        end if
      enddo

!     if (lprnt) print *,' tseaz2=',tsea(ipr),' tref=',tref(ipr),   &
!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

      return
      end subroutine sfc_nst_post_run

      end module sfc_nst_post
