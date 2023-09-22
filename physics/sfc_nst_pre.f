!>  \file sfc_nst_pre.f
!!  This file contains preparation for the GFS NSST model.

      module sfc_nst_pre

      contains

! \defgroup GFS_NSST_PRE GFS Near-Surface Sea Temperature Pre
!!
!! The NSST scheme is one of the three schemes used to represent the
!! surface in the GFS physics suite. The other two are the Noah land
!! surface model and the sice simplified ice model.
!!
      
!! \section arg_table_sfc_nst_pre_run Argument Table
!! \htmlinclude sfc_nst_pre_run.html
!!
!> \section NSST_general_pre_algorithm General Algorithm
!! @{
      subroutine sfc_nst_pre_run
     &    (im, wet, tgice, tsfco, tsurf_wat,
     &     tseal, xt, xz, dt_cool, z_c, tref, cplflx,
     &     xlon, xlat, oceanfrac, nthreads, errmsg, errflg)

      use machine , only : kind_phys
      use module_nst_water_prop, only: get_dtzm_2d
      use module_nst_parameters, only : rad2deg

      implicit none

      integer, parameter :: kp = kind_phys

!  ---  inputs:
      integer, intent(in) :: im, nthreads
      logical, dimension(:), intent(in) :: wet
      real (kind=kind_phys), intent(in) :: tgice
      real (kind=kind_phys), dimension(:), intent(in) ::
     &     xlon, xlat, tsfco, xt, xz, dt_cool, z_c, oceanfrac
      logical, intent(in) :: cplflx

!  ---  input/outputs:
      real (kind=kind_phys), dimension(:), intent(out) ::
     &    tsurf_wat, tseal, tref

!  ---  outputs:
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals
      integer :: i
      real(kind=kind_phys), parameter :: zero = 0.0_kp,
     &                                   one  = 1.0_kp,
     &                                   half = 0.5_kp,
     &                                   omz1 = 2.0_kp
      real(kind=kind_phys) :: tem1, tem2, dnsst, alon, alat
      real(kind=kind_phys), dimension(im) :: dtzm, z_c_0

      logical :: doprint
      real(kind=kind_phys) :: frz=273.15, small=.05, testlon, testlat
      common /testpt/ testlon,testlat           ! (values defined in dcyc2t3.f)
      doprint(alon,alat)=abs(testlon-alon).lt.small .and.
     &                   abs(testlat-alat).lt.small

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (wet(i) .and. oceanfrac(i) > 0.0) then
!          tem         = (oro(i)-oro_uf(i)) * rlapse
          ! DH* 20190927 simplyfing this code because tem is zero
          !tem          = zero
          !tseal(i)     = tsfco(i)  + tem
          tseal(i)      = tsfco(i)
          !tsurf_wat(i) = tsurf_wat(i) + tem
          ! *DH
        endif
      enddo
!
!   update tsfc & tref with T1 from OGCM & NSST Profile if coupled
!
      if (cplflx) then
        z_c_0 = zero
        call get_dtzm_2d (xt,  xz, dt_cool,                             &
     &                    z_c_0, wet, zero, omz1, im, 1, nthreads, dtzm)
        do i=1,im
         if (wet(i) .and. oceanfrac(i) > zero ) then
!           dnsst   = tsfc_wat(i) - tref(i)                 !  retrive/get difference of Ts and Tf
            tref(i) = max(tgice, tsfco(i) - dtzm(i))        !  update Tf with T1 and NSST T-Profile
!           tsfc_wat(i) = max(271.2,tref(i) + dnsst)        !  get Ts updated due to Tf update
!           tseal(i)    = tsfc_wat(i)
            if (abs(xz(i)) > zero) then
              tem2 = one / xz(i)
            else
              tem2 = zero
            endif
            tseal(i)     = tref(i) + (xt(i)+xt(i)) * tem2 - dt_cool(i)
            tsurf_wat(i) = tseal(i)
          endif
        enddo
      endif

      do i=1,im
        alon=xlon(i)*rad2deg
        alat=xlat(i)*rad2deg
        if (doprint(alon,alat))
     &    print 98,'exiting sfc_nst_run_pre   lon,lat=',alon,alat,
     &     'tsfco',tsfco(i)-frz,
     &     'tseal',tseal(i)-frz,
     &     'tref',tref(i)-frz,
     &     'tsurf_w',tsurf_wat(i)-frz
 99     format (/a,2f7.2/(5(a8,"=",f7.2)))
 98     format (/a,2f7.2/(4(a8,"=",es11.4)))
      end do

      return
      end subroutine sfc_nst_pre_run
!! @}
      end module sfc_nst_pre
