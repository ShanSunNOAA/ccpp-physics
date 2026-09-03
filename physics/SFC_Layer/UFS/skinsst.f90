!>\file skinsst.f90
!! This file contains Rainer's skin temperature scheme.

module state_eqn
  implicit none
  ! --- coefficients for sigma-0 (based on Brydon & Sun fit, JGR 1999)
  real, parameter, dimension(7) :: coef = (/ &
       -1.36471E-01,  4.68181E-02,  8.07004E-01, -7.45353E-03, -2.94418E-03, &
        3.43570E-05,  3.48658E-05 /)
contains

   real function sig(t,s)
! --- sea water density (sigma = density - 1000) at p=0
   real, intent(in) :: t, s
!  sig = coef(1)+s*coef(3)+				&
!     t*(coef(2)+s*coef(5)+				&
!     t*(coef(4)+s*coef(7)+t*coef(6)))
   sig = t*(t*(s*coef(7)+t*coef(6)+coef(4))		& ! alt.grouping
              +s*coef(5)          +coef(2))		&
              +s*coef(3)          +coef(1)
   return
   end function sig

   real function dsigdt(t,s)
! --- thermal expansion coefficient
   real, intent(in) :: t, s
!  dsigdt = coef(2)+s*coef(5)+2.*t*(coef(4)+s*coef(7)+1.5*t*coef(6))
   dsigdt = 2.*t*(1.5*t*coef(6)+s*coef(7)+coef(4))+s*coef(5)+coef(2) ! alt.grouping
   return
   end function dsigdt

   real function dsigds(t,s)
! --- saline contraction coefficient
   real, intent(in) :: t, s
!  dsigds = coef(3)+t*(coef(5)+t*coef(7))
   dsigds = t*(t*coef(7)+coef(5))+coef(3)		 ! alt.grouping
   return
   end function dsigds
end module state_eqn


module skinsst
  use machine,                 only : kind_phys
  use state_eqn
  use funcphys,                only : fpvs
  use module_nst_parameters,   only : kw => tc_w, visw, cp_w, z_c_max, z_c_ini, ustar_a_min

  implicit none
  private
  public :: skinsst_init, skinsst_run, skinsst_finalize

contains

  subroutine skinsst_init()
  end subroutine skinsst_init

  subroutine skinsst_finalize()
  end subroutine skinsst_finalize

!>\defgroup gfs_ocean_main GFS Simple Ocean Scheme Module
!! This subroutine calculates thermodynamical properties over
!! open water.
!! \section arg_table_skinsst_run Argument Table
!! \htmlinclude skinsst_run.html
!!

  subroutine skinsst_run(                                         &
    im,              & ! horiz. loop extent                 in
    iter,            & ! ccpp loop counter                  in
    wet,             & ! .true. at ocean & lake points      inout
    oceanfrac,       & ! cell portion covered by ocean      in
    timestep,        & ! model timestep                     in
    xlon, xlat,      & ! longitude, latitude                in
    sfcemis,         & ! sea surface emissivity             in
    ulwflx,          & ! upwelling LW flux                  inout
    dlwflx,          & ! absorbed downwelling LW flux       in
    sfcnsw,          & ! net SW flux, pos.down              in
    tsfco,           & ! ocean/lake top layer temperature   inout
    psfc,            & ! surface pressure                   in
    wind,            & ! atm. mid-layer 1 wind              in
    stress,          & ! wind stress (N/m^2)                in
    plyr1,           & ! atm. mid-layer 1 pressure          in
    tlyr1,           & ! atm. mid-layer 1 temperature       in
    qlyr1,           & ! atm. mid-layer 1 humidity          in
    ulyr1,           & ! atm. mid-layer 1 zonal wind        in
    vlyr1,           & ! atm. mid-layer 1 meridional wind   in
    compres,         & ! adiabat.compression factor         in
    cm,              & ! drag coeff for momentum            in
    ch,              & ! drag coeff for heat and moisture   in
    hvap,            & ! latent heat of evaporation         in
    cp,              & ! specif.heat of air                 in
    rd,              & ! gas constant of dry air            in
    eps,             & ! ratio of gas constants, rd/rv      in
    sbc,             & ! stefan-boltzmann constant          in
    lakedepth,       & ! lakedepth                          in
    tskin,           & ! skin temp                          inout
    skinold,         & ! previous tskin                     inout
    temwat,          & ! lake mixed layer temperature       inout
    xtinct,          & ! extinction coefficient             inout
    hice,            & ! sea/lake ice thickness             inout
    fice,            & ! sea/lake ice fraction              inout
    tisfc,           & ! sea/lake ice temperature           inout
    islmsk,          & ! ocean/land/ice = 0/1/2             inout
    islmsk_cice,     & ! ocean/land/ice = 0/1/2             inout
    flag_lakefreeze, & ! flag for lake freeze               inout
    icy,             & ! flag_nonzero_sea_ice_sfc_fraction  inout 
    flxold,          & ! previous lake sfc heat flux        inout
    dt_cool,         & ! skin layer cooling amount          inout
    dt_warm,         & ! warm-layer surface warming amount  out
    z_c,             & ! sub-layer cooling thickness        out
    qsat,            & ! saturation specif. humidity        out
    evap,            & ! kinematic latent heat flux, pos.up out
    hflx,            & ! kinematic sensible heat flux       out
    ep,              & ! potential latent heat flux         out
    cmm,             & ! momentum exchange coeff            out
    chh,             & ! thermal exchange coeff             out
    lseaspray,       & ! sea spray flag                     in
    fm,              & ! Monin-Obukhov function at surface  in
    fm10,            & ! Monin-Obukhov function at 10m      in
    errmsg, errflg)

    implicit none

    ! --- input:
    integer, intent(in)                                 :: im, iter
    logical, intent(in)                                 :: lseaspray
    real(kind=kind_phys), dimension(:), intent(in)      :: xlon, xlat, &
         sfcemis, dlwflx, sfcnsw, wind, psfc, plyr1, tlyr1, qlyr1,     &
         ulyr1, vlyr1, cm, ch, compres, stress, fm, fm10, oceanfrac, lakedepth
    real(kind=kind_phys), intent(in)                    :: timestep, hvap, cp, rd, eps, sbc

    ! --- inout:
    real(kind=kind_phys), dimension(:), intent(inout)   :: ulwflx, tsfco, tskin, dt_cool, hice, fice, tisfc
    real(kind=kind_phys), dimension(:), intent(inout)   :: &
         skinold,   & ! previous skin temperature
         xtinct,    & ! SW extinction coefficient
         temwat,    & ! lake mixed layer temperature
         flxold       ! previous lake surface heatflux
    integer, dimension(:), intent(inout)                :: islmsk, islmsk_cice
    logical, dimension(:), intent(inout)                :: wet, icy, flag_lakefreeze

    ! --- output:
    real(kind=kind_phys), dimension(:), intent(out)     :: evap, hflx, ep, qsat, cmm, chh, dt_warm, z_c
    character(len=*), intent(out)                       :: errmsg
    integer, intent(out)                                :: errflg

    ! --- locals:
    integer :: i, n, loop
    real :: alon, alat, virt, rho_air, rho_wat, pvap, tsq, piston, vel, &
               vertdf, & ! vertical temperature difference
               dfloss, & ! heat loss by downward diffusion
               nonsol    ! sum of nonsolar air-sea fluxes (pos.up)
               
    real, parameter :: spcifh = 3990. ! seawater specific heat
    real, parameter :: grav   = 9.806 ! gravity
    real, parameter :: sss    = 34.7  ! sea surface salinity

    integer, parameter :: itmax = 5   ! regula falsi iterations
    real :: rnl_ts, hs_ts, rf_ts, alpha, beta, rch, ustar,              &
            hist(0:itmax) = 0., x1, x2, x3, y1, y2, dif1, dif2, dif3

    ! variables for sea spray effect
    real(kind=kind_phys) :: f10m, u10m, v10m, ws10, ru10, qss1, bb1, hflxs, evaps, ptem, tem
    real(kind=kind_phys), parameter :: alps=0.75, bets=0.75, gams=0.15, &
                                       ws10cr=30., conlf=7.2e-9, consf=6.4e-8

    logical :: doprint, details, frstrip
    real(kind=kind_phys) :: kd_par, frz=273.15, small=.05, totflx,      &
                            oldflx, testlon, testlat, hice_old, hice_min=0.1
    external kd_par

    real, parameter :: rad2deg = 57.2957795
    real, parameter :: dz      = 2.0        ! nominal z increment in diffusion eqn
    real, parameter :: dffus   = 1.43e-7    ! thermal diffusivity (m^2/sec)
    real, parameter :: wipout  = 900.       ! relax.time (sec) for warm-lyr wipeout
    real, parameter :: homog   = 8.         ! wind speed needed for homogenization

    real(kind=kind_phys), parameter :: alon1=107.17, alat1=52.79
    real(kind=kind_phys), parameter :: alon2=273.58, alat2=47.57

    ! --- piston velocity: molecular diffusion when vel = 0.
    ! --- piston vel. set to wipe out warm layer when vel > homog
    !piston(vel) = dffus / dz + min(1., vel / homog)      * dz / wipout  ! linear
     piston(vel) = dffus / dz + min(1., vel / homog)**1.5 * dz / wipout  ! non-linear
    !piston(vel) = dffus / dz + min(1., vel / homog)**2   * dz / wipout  ! quadratic

    doprint(alon, alat) = abs(testlon - alon) < small .and. abs(testlat - alat) < small

    if (iter > 1) return

    errflg = 0
    errmsg = ""
    call get_testpt(testlon, testlat)

    do i = 1, im
      if (wet(i)) then
        alon = xlon(i) * rad2deg
        alat = xlat(i) * rad2deg

      if (doprint(alon,alat)) then
        print 99,'entering skinsst_run   lon,lat=',alon,alat,             &
       'temwat',temwat(i)-frz,           & ! lake water temperature
       'xtinct',xtinct(i),               & ! extinction coefficient
       'hice',hice(i),                   & ! ice thickness
       'ocnfrac',oceanfrac(i),           & ! ocean fraction
       'stress',stress(i),               & ! wind stress (N/m^2)
       'sfcemis',sfcemis(i),             & ! sfc emissivity
       'wind',wind(i),                   & ! surface wind
       'pstonE3',piston(wind(i))*1.e3,   & ! piston velocity
       'sfcnsw',sfcnsw(i),               & ! total sky net SW flx into ocean
       'dlwflx',dlwflx(i),               & ! absorbed downwelling LW flux
       'ulwflx',ulwflx(i),               & ! upwelling LW flux
       'psfc',psfc(i)*.01,               & ! surface pressure (mb)
       'plyr1',plyr1(i)*.01,             & ! atm.layer 1 presure
       'tlyr1',tlyr1(i)-frz,             & ! atm.layer 1 air temp
       'qlyr1',qlyr1(i)*1.e3,            & ! atm.layer 1 humidity (g/kg)
       'sigma_t',sig(tlyr1(i)-frz,sss),  & ! sea water density - 1000
       'compres',compres(i),             & ! midlyr-to-sfc adiab.compression
       'skinold',skinold(i)-frz,         & ! previous tskin
       'dcoolE2',dt_cool(i)*100.,        & ! previous dtcool
       'tsfco',tsfco(i)-frz                ! ocean top layer temperature
        print '(5(a13,"=",l2))','lseaspray',lseaspray
        if (oceanfrac(i)==0.) print '(2f7.2,a)',alon,alat,' is -lake- point'
      end if
 99  format (/a,2f7.2/(15(a8,"=",f7.2)))
 98  format (/a,2f7.2/(4(a8,"=",es11.4)))
 97  format (/a,2f7.2/(4(a8,"=",f11.6)))

        virt    = tlyr1(i) * (1. + (eps - 1.) * qlyr1(i))
        rho_air = plyr1(i) / (rd * virt)
        rch     = rho_air * cp * ch(i) * wind(i)  ! W/m^2/deg
        cmm(i)  = cm(i) * wind(i)
        chh(i)  = rho_air * ch(i) * wind(i)
        ep(i)   = 0.
        rho_wat = 1000. + sig(tsfco(i) - frz, sss)
        alpha   = -dsigdt(tsfco(i) - frz, sss) / rho_wat
        beta    =  dsigds(tsfco(i) - frz, sss) / rho_wat
        ustar   = sqrt(stress(i) / rho_air)       ! air friction velocity

        if (skinold(i) == 0.) then              ! use skinold=0 as indicator for t=0
          frstrip   = .true.
          dt_cool(i)= 0.
          tskin(i)  = tsfco(i)
          temwat(i) = tsfco(i)                    ! lake temp
          xtinct(i) = kd_par(alon,alat)           ! from Son & Wang (2015)
          flxold(i) = 0.                          ! old heat flux over lake
          vertdf  = 0.
          dfloss  = 0.

          if (oceanfrac(i) == 0. ) then           ! lake points
            hice(i) = fice(i)*hice(i)
            if (hice(i) < 0.015) then             ! ignore ice if volume < 0.1m at fice=15%
              hice(i)=0.
              fice(i)=0.
            else
              hice(i) = max(hice(i),hice_min)
              fice(i) = 0.99
            end if
          end if

        else
          frstrip  = .false.
          tskin(i) = skinold(i)                   ! previous tskin
        end if

        details = .false.

        ! --- bypass warm layer calculation if tskin is below zero or less than tsfco
        if (tskin(i) < frz .or. tskin(i) <= tsfco(i) - dt_cool(i)) then
          tskin(i) = tsfco(i) - dt_cool(i)
          dfloss   = 0.
        else
          ! --- surface cooling by downward heat diffusion
          vertdf = tskin(i) + dt_cool(i) - tsfco(i)
          dfloss = vertdf * timestep * piston(wind(i)) / dz

          if (sfcnsw(i) <= 0.) dfloss = max(0.01, dfloss)
          dfloss   = min(dfloss, vertdf)
          tskin(i) = tskin(i) - dfloss
        end if

        dt_warm(i) = 0.
        if (sfcnsw(i) > 0.) then               ! daytime
          ! --- evaluate warm-layer increment
          dt_warm(i) = sfcnsw(i) * timestep * xtinct(i) / (rho_wat * spcifh)
          tskin(i)   = tskin(i) + dt_warm(i)      ! dt_warm is cumulative
        end if

        ! --- start cool-skin iteration, using REGULA FALSI (aiming for x_n = y_n)
        ! --- x1, x2, x3, y1, y2 are consecutive dt_cool approximations.
        x1 = -.5
        x2 = +.5

        call surflx(nonsol, tskin(i) + x1, tlyr1(i) * compres(i), qlyr1(i),   &
             psfc(i), hflx(i), qsat(i), evap(i), hvap / cp, eps, rch, sbc,    &
             sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, details)

        call coolskin(ustar, nonsol, sfcnsw(i), evap(i), sss, alpha,          &
                      beta, rho_wat, rho_air, tskin(i) + x1, grav, hvap,      &
                      y1, z_c(i), alon, alat, details)

        dif1 = y1 - x1
        if (y1 .ne. 0.) then
          do loop = 1, itmax
            call surflx(nonsol, tskin(i) + x2, tlyr1(i) * compres(i), qlyr1(i), &
                 psfc(i), hflx(i), qsat(i), evap(i), hvap / cp, eps, rch, sbc,  &
                 sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, details)

            call coolskin(ustar, nonsol, sfcnsw(i), evap(i), sss, alpha,        &
                          beta, rho_wat, rho_air, tskin(i) + x2, grav, hvap,    &
                          y2, z_c(i), alon, alat, details)
            
            dif2 = y2 - x2

            if (details) print '(a,3es11.3,i7)','(skinsst)  x1,y1,y1-x1 =', x1, y1, dif1, loop
            if (details) print '(a,3es11.3,i7)','(skinsst)  x2,y2,y2-x2 =', x2, y2, dif2, loop

            x3 = (x1 * dif2 - x2 * dif1) / (dif2 - dif1)  ! regula falsi

            if (abs(dif2) > 1.e-4) then                ! test for convergence
              if (abs(dif1) > abs(dif2)) then
                x1   = x2
                y1   = y2
                dif1 = dif2
              end if
              x2 = x3
            else
              ! --- dt_cool is not cumulative => subtract previous dt_cool from tskin
              tskin(i)   = tskin(i) + dt_cool(i) - y2     ! new minus old
              dt_cool(i) = y2                             ! save new dt_cool for next step
              exit                                        ! all done
            end if
            hist(loop) = y2
          end do

          if (abs(dif2) > 1.e-4) then
            print '(a,3f8.2/(11f7.2))', 'dt_cool not converging at lon,lat', &
                  alon, alat, hist(loop), (hist(n), n=1,loop)
          end if
        end if ! y1 nonzero

        if (oceanfrac(i) == 0. .and. .not. frstrip) then
          call surflx(nonsol, tskin(i), tlyr1(i) * compres(i), qlyr1(i),     &
             psfc(i), hflx(i), qsat(i), evap(i), hvap / cp, eps, rch, sbc,   &
             sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, details)

          ! --- use rudimentary energy loan lake model 'enloan' 
          ! --- it forms new ice or melt ice, let ice3lay calculate ice thickness
          totflx = sfcnsw(i) - nonsol

          ! --- average totflx over 2 time steps to suppress comput.mode in enloan
          oldflx    = flxold(i)
          flxold(i) = totflx
          totflx    = .5 * (totflx + oldflx)

          if (hice(i) > 0.) hice(i) = max(hice(i),hice_min) !follow sfc_sice
          hice_old = hice(i)
          call enloan(timestep, totflx, hice(i), tskin(i), skinold(i),       &
                      temwat(i), lakedepth(i), alon, alat, doprint(alon,alat))
          tsfco(i) = temwat(i)

       !  if ((hice(i) >  0. .and. hice_old == 0.) .or.              &
       !      (hice(i) == 0. .and. hice_old >  0.))                  &
       !    print 90, 'lakeice change lon,lat=', alon, alat,         &
       !         'hice_old', hice_old,                               &
       !         'hice', hice(i),                                    &
       !         'temwat', temwat(i) - frz,                          &
       !         'tskin', tskin(i) - frz,                            &
       !         'tisfc', tisfc(i) - frz,                            &
       !         'islmsk', islmsk(i),                                &
       !         'wet', wet(i),                                      &
       !         'icy', icy(i)
90        format (/a,2f7.2,5(a8,"=",f9.4),a8,"=",i2,2(a8,"=",l2))

       !  if ((abs(alon-alon1) < small .and. abs(alat-alat1) < small) .or. &
       !      (abs(alon-alon2) < small .and. abs(alat-alat2) < small))     &
       !    print 91,'qq1 skin  lon,lat=', alon, alat,                     &
       !    'hice',hice(i),                                                &
       !    'hice_o',hice_old,                                             &
       !    'fice',fice(i),                                                &
       !    'islmsk',islmsk(i)
91        format (/a,2f7.2,3(a8,"=",f9.4),a8,"=",i2)

          flag_lakefreeze(i) = .false.
          if (hice_old == 0. .and. hice(i) > 0.) then ! new ice
            hice(i) = max(hice(i), hice_min)
            fice(i) = .99
            tisfc(i) = tskin(i)
            flag_lakefreeze(i) = .true.
          end if

          if (hice_old > 0. ) then
                  if (hice(i)  >= hice_min) then ! no change in hice if icy both bef & aft
              hice(i) = hice_old
            else ! melting
              hice(i) = 0.
              fice(i) = 0.
              flag_lakefreeze(i) = .true.
            end if
          end if

          if (hice(i) > 0. ) then
            evap(i)        = 0.
            fice(i)        = 0.99
            islmsk(i)      = 2
            islmsk_cice(i) = 2
            icy(i)         = .true.
          else
            fice(i)        = 0.
            islmsk(i)      = 0
            islmsk_cice(i) = 0
            icy(i)         = .false.
          end if
         
         !  if ((abs(alon-alon1) < small .and. abs(alat-alat1) < small) .or. &
         !    (abs(alon-alon2) < small .and. abs(alat-alat2) < small))       &
         !  print 91,'qq2 skin  lon,lat=', alon, alat,                       &
         !  'hice',hice(i),                                                  &
         !  'hice_o',hice_old,                                               &
         !  'fice',fice(i),                                                  &
         !  'islmsk',islmsk(i)
        end if ! enloan

        skinold(i) = tskin(i)

 ! --- according to GFS_surface_composites_inter.F90,
 ! --- dlwflx is the absorbed portion of downwelling LW flux.
 ! --- hence, the total downwelling flux is dlwflx/sfcemis
 ! --- and the reflected part is (1-sfcemis)*dlwflx/sfcemis

        ulwflx(i) = ulwflx(i) + dlwflx(i) * (1. - sfcemis(i)) / sfcemis(i)

        if (lseaspray) then
          f10m = fm10(i) / fm(i)
          u10m = f10m * ulyr1(i)
          v10m = f10m * vlyr1(i)
          ws10 = max(min(sqrt(u10m*u10m + v10m*v10m), ws10cr), 1.)
          
          tem  = .015 * ws10 * ws10
          ru10 = 1. - .087 * log(10. / tem)
          qss1 = fpvs(tlyr1(i))
          qss1 = eps * qss1 / (plyr1(i) + (eps - 1.) * qss1)
          
          tem  = rd * cp * tlyr1(i) * tlyr1(i)
          tem  = 1. + eps * hvap * hvap * qss1 / tem
          bb1  = 1. / tem
          
          evaps = conlf * (ws10**5.4) * ru10 * bb1
          evaps = evaps * rho_air * hvap * (qss1 - qlyr1(i))
          evap(i) = evap(i) + alps * evaps
          
          hflxs = consf * (ws10**3.4) * ru10
          hflxs = hflxs * rho_air * cp * (tskin(i) - tlyr1(i))
          ptem  = alps - gams
          hflx(i) = hflx(i) + bets * hflxs - ptem * evaps
        endif

        if (doprint(alon,alat))                                           &
        print 99,'exiting skinsst_run   lon,lat=',alon,alat,              &
        'virt',virt-frz,                  & ! virtual air temp
        'rho_air',rho_air,                & ! air density
        'pvap',pvap,                      & ! satur. vapor pressure (mb)
        'qsat',qsat(i),                   & ! satur. specif.humidity
        'hflx',hflx(i),                   & ! sensible heat flux
        'evap',evap(i),                   & ! latent heat flux
        'nonsol',nonsol,                  & ! net non-solar surface flux
        'sfcnsw',sfcnsw(i),               & ! net solar surface flux
        'ulwflx',ulwflx(i),               & ! upwelling LW flux
        'dwarmE2',dt_warm(i)*100.,        & ! temperature increment due to SW
        'dcoolE2',dt_cool(i)*100.,        & ! cool-skin temperature correction
        'tskin',tskin(i)-frz,             & ! skin temperature
        'vertdE2',vertdf*100.,            & ! difference tskin - tsfco
        'dflossE2',dfloss*100.,           & ! heat loss by dnwd diffusion
        'tsfco',tsfco(i)-frz                ! ocean top layer temperature

        ! --- save (tskin - top layer T) for diagnostic purposes
        dt_warm(i) = tskin(i) - tsfco(i)

        ! --- convert fluxes from W/m^2 to "kinematic" (velocity x fluxed variable)
        hflx(i) = hflx(i) / (rho_air * cp)   ! deg m/sec
        evap(i) = evap(i) / (rho_air * hvap) ! m/sec
      end if ! wet
    end do   ! im loop

    do i = 1, im
      if (fice(i) == 1. .and. oceanfrac(i) == 0.) then ! 100% lake ice
        fice(i) = 0.99    ! Allow for a tiny amount of "wet" by setting fice(i) = 0.99
        wet(i) = .true.
        tskin(i) = frz
        flag_lakefreeze(i) = .true.
      end if ! Checking for 100% lake ice
    end do   ! im loop

    return
  end subroutine skinsst_run


  subroutine coolskin(ustar_a, f_nsol, f_sol_0, evap, sss, alpha, beta, &
                      rho_w, rho_a, ts, grav, latnt, deltat_c, z_c, alon, alat, doprint)
    ! upper ocean cool-skin parameterization, Fairall et al, 1996.
    implicit none

    logical, intent(in) :: doprint
    real(kind=kind_phys), intent(in)  :: ustar_a, f_nsol, f_sol_0, evap, &
                                         sss, alpha, beta, rho_w, rho_a, ts, grav, latnt, alon, alat
    real(kind=kind_phys), intent(out) :: deltat_c, z_c

    real(kind=kind_phys), parameter :: frz=273.15
    real(kind=kind_phys) :: xi, hb, ustar1_a, bigc, deltaf, fxp

    if (doprint) print 98, 'entering coolskin   lon,lat=', alon, alat, &
         'ustar_a', ustar_a, 'f_nsol', f_nsol, 'f_sol_0', f_sol_0,     &
         'evap', evap, 'rho_w', rho_w, 'rho_a', rho_a, 'ts', ts - frz, 'zc_ini', z_c_ini * 1.e3

98  format (/a,2f7.2/(4(a8,"=",es11.4)))

    z_c      = z_c_ini ! initial guess
    ustar1_a = max(ustar_a, ustar_a_min)

    call sw_rad_skin(z_c, fxp)
    deltaf = f_sol_0 * fxp
    hb     = alpha * (f_nsol - deltaf) + beta * sss * cp_w * evap / latnt
    bigc   = 16. * grav * cp_w * (rho_w * visw)**3 / (rho_a * kw)**2

    if (hb > 0) then
      xi = 6. / (1. + (bigc * hb / ustar1_a**4)**0.75)**0.3333333
    else
      xi = 6.0
    endif
    
    z_c = min(z_c_max, xi * visw / (sqrt(rho_a / rho_w) * ustar1_a))
    call sw_rad_skin(z_c, fxp)

    deltaf = f_sol_0 * fxp
    deltaf = f_nsol - deltaf
    
    if (deltaf > 0) then
      deltat_c = deltaf * z_c / kw
    else
      deltat_c = 0.
      z_c      = 0.
    endif

    if (doprint) print 98, 'exiting coolskin   lon,lat=', alon, alat, &
         'fxp', fxp, 'deltaf', deltaf, 'hb', hb, 'bigc', bigc, 'xi', xi, &
         'delt_c', deltat_c, 'z_c', z_c

    return
  end subroutine coolskin


  subroutine sw_rad_skin(z, fxp)
  ! original name: elemental subroutine sw_ohlmann_v1(z,fxp)
  ! fraction of the solar radiation absorbed by the ocean at the depth z
    implicit none
    real(kind=kind_phys), intent(in)  :: z
    real(kind=kind_phys), intent(out) :: fxp

    if (z > 0) then
       fxp = .065 + 11. * z - 6.6e-5 / z * (1. - exp(-z / 8.0e-4))
    else
       fxp = 0.
    endif
  end subroutine sw_rad_skin


  subroutine surflx(    &
   nonsol,              & ! sum of nonsolar heat fluxes, pos.up
   tskin,               & ! skin temperature
   tlyr1,               & ! temperature in lowest atmo layer
   qlyr1,               & ! sfc.humidity in lowest atmo layer
   psfc,                & ! surface pressure
   hflx,                & ! sensible heat flux, pos.up                  (out)
   qsat,                & ! satur.specf. humidity                       (out)
   evap,                & ! latent heat flux, pos.up                    (out)
   elocp,               & ! heat of evaporation over specif.heat, hvap/cp
   eps,                 & ! ratio of air/vapor gas constants
   rch,                 & ! rho * cp * ch * wind  [W/deg]
   sbc,                 & ! stefan-boltzmann constant
   sfcemis,             & ! sea surface emissivity
   dlwflx,              & ! absorbed downwelling LW flux, pos.down
   ulwflx,              & ! surface-emitted LW flux, pos.up             (out)
   alon,alat,doprint)

! --- compute sum of nonsolar air-sea fluxes
! --- watch out for nonstandard sign convention:
! --- dlwflx is pos.down, all other fluxes pos.up.

    implicit none
    
    real, intent(in)    :: tskin, tlyr1, qlyr1, psfc, elocp, eps, rch, sbc, sfcemis, dlwflx, alon, alat
    logical, intent(in) :: doprint
    real, intent(out)   :: nonsol, qsat, evap, hflx, ulwflx
    real                :: pvap, frz=273.15

    if (doprint) print 99, 'entering surflx   lon,lat=', alon, alat, &
         'nonsol', nonsol, 'tskin', tskin - frz, 'tlyr1', tlyr1 - frz

    pvap   = fpvs(tskin)
    qsat   = eps * pvap / (psfc + (eps - 1.) * pvap)
    evap   = elocp * rch * (qsat - qlyr1)
    hflx   = rch * (tskin - tlyr1)
    ulwflx = sfcemis * sbc * tskin**4
    nonsol = hflx + evap + ulwflx - dlwflx

    if (doprint) print 99,'exiting surflx   lon,lat=',alon,alat,          &
     'tskin',tskin-frz, & ! skin temperature
     'psfc',psfc*.01,   & ! surface pressure (mb)
     'pvap',pvap,       & ! saturation vapor pressure
     'qsat',qsat*1.e3,  & !saturation specif. humidity (g/kg)
     'evap',evap,       & ! latent heat flux
     'hflx',hflx,       & ! sensible heat flux
     'ulwflx',ulwflx,   & ! upwelling long-wave flux
     'dlwflx',dlwflx,   & ! downwelling long-wave flux
     'nonsol',nonsol      ! sum of nonsolar heat fluxes (pos.up)
 99  format (/a,2f7.2/(5(a8,"=",f7.2)))
 98  format (/a,2f7.2/(4(a8,"=",es11.4)))

    return
  end subroutine surflx


  subroutine enloan(delt, surflx, hice, temice, skinold, temwat, lakedepth, alon, alat, doprint)
! --- single-column version of 'energy loan' ice model.
! --- ice amount represents energy 'loaned' to water column to prevent
! --- wintertime cooling below freezing level. 'loan' is paid back in summer.

    implicit none
    
    real, parameter     :: frz=273.15
    logical, intent(in) :: doprint
    real, intent(in)    :: delt, alon, alat, lakedepth
    real, intent(inout) ::                              &
    surflx,      & ! net total heat flux between atm and ice (W/m^2)
    hice,        & ! grid-box averaged ice thickness (m)
    temwat,      & ! mixed layer temperaure
    temice,      & ! ice surface temperature
    skinold        ! previous ice surface temperature

    real ::              &
    tmelt=frz-.2,        & ! melting point (deg K)
    thin=.01,            & ! min.ice thickness
    rhoice=917.,         & ! ice density (kg/m^3)
    rhowat=1000.,        & ! water density (kg/m^3)
    kice=2.04,           & ! heat conductivity in ice (W/m/deg)
    fusion=334.e3,       & ! latent heat of fusion (J/kg)
    rate=.2/3600.,       & ! max. ice melting rate (m/sec)
!   fluctn=3./3600,      & ! limit on temice fluctuation (deg/sec)
    fluctn=2./3600,      & ! limit on temice fluctuation (deg/sec)
    spcifh=4190.,        & ! specific heat of water (J/kg/deg)
    dpmin=40.              ! nominal mixed layer depth (m)
    real :: tnew, borrow, paybak, avail, dpth

! --- energy loan: add extra energy to the ocean to keep SST from dropping
! --- below tmelt in winter. return this borrowed energy to the 'energy bank'
! --- in summer as quickly as surflx > 0 allows.

   if (doprint) print 97,'entering enloan     lon,lat=',alon,alat,      &
    'surflx',surflx,                                                    &
    'temwat',temwat-frz,                                                &
    'hice',hice,                                                        &
    'temice',temice-frz,                                                &
    'skinold',skinold - frz
99  format (/a,2f7.2/(5(a8,"=",f7.2)))
98  format (/a,2f7.2/(5(a8,"=",es11.4)))
97  format (/a,2f7.2/(5(a8,"=",f11.6)))

    borrow = 0.
    paybak = 0.
    dpth   = min(dpmin, max(20., lakedepth))
    tnew   = temwat + surflx * delt / (rhowat * spcifh * dpth)

    if (surflx < 0.) then        ! cooling
      if (tnew > tmelt) then 
        if (doprint) print 97, 'enloan action 1     lon,lat=', alon, alat, 'tnew', tnew - frz, 'temwat', temwat - frz
        temwat = tnew
      else
        ! --- borrow energy to keep temwat from dropping below tmelt
        borrow = (tmelt - tnew) * rhowat * spcifh * dpth / delt
        temwat = tmelt
        hice   = hice + borrow * delt / (rhoice * fusion)
        if (doprint) print 97,'enloan action 2     lon,lat=',alon,alat,    &
          'tnew',tnew-frz,'borrow',borrow,'surflx',surflx,'hice',hice,     &
          'temwat',temwat-frz
      end if
    else                            ! warming
      if (hice > 0.) then
        ! --- return the borrowed amount whenever tnew > tmelt
        avail  = (tnew - tmelt) * rhowat * spcifh * dpth / delt                    ! W/m2
        paybak = min(hice * rhoice * fusion / delt, avail, rate * rhoice * fusion) ! W/m2
        hice   = hice - paybak * delt / (rhoice * fusion)                          ! m
        temwat = tnew - paybak * delt / (rhowat * spcifh * dpth)                   ! K
        if (doprint) print 97,'enloan action 3     lon,lat=',alon,alat,    &
          'tnew',tnew-frz,'temwat',temwat-frz,'temice',temice-frz,         &
          'paybak',paybak,'surflx',surflx,'hice',hice
      else
        temwat = tnew
      end if
    end if         ! surflx

    ! --- compute ice surface temperature
    if (hice > thin) then
! --- assume zero flux divergence at ice surface, so
! --- surflx = (temice-temwat)*kice/hice
      temice = min(tmelt, temwat + hice * surflx / kice)
! --- put limits on temice tendency
      temice = max(skinold - fluctn * delt, min(skinold + fluctn * delt, temice))

     if (doprint) print 97,'enloan action 4     lon,lat=',alon,alat,     &
      'tnew',tnew-frz,'flxice',surflx,'hice',hice,'temice',temice-frz
    else
      temice = temwat
    end if
    skinold = temice

    if (doprint) print 97,'exiting enloan     lon,lat=',alon,alat,       &
    'temwat',temwat-frz,                                                 &
    'hice',hice,                                                         &
    'temice',temice-frz,                                                 &
    'tnew',tnew-frz

    return
  end subroutine enloan

end module skinsst


subroutine get_testpt(testlon, testlat)
  ! --- define test point for detailed diagnostic
  implicit none
  real, intent(out) :: testlon, testlat

  ! Initialize to dummy values since original test points were removed
  testlon = -999.0
  testlat = -999.0

!  testlon = 279.12 ; testlat = 42.00  ! high-lat lake
!  testlon = 273.58 ; testlat = 47.57  ! high-lat lake
!  testlon = 274.14 ; testlat = 46.89  ! lake superior
!  testlon = 278.88 ; testlat = 42.24  ! lake erie
!  testlon = 48.48  ; testlat = 44.02  ! caspian sea
!  testlon = 53.44  ; testlat = 37.48  ! caspian sea
!  testlon = 272.60 ; testlat = 48.70  ! lake superior
!  testlon = 273.39 ; testlat = 47.79  ! lake superior
!  testlon = 271.92 ; testlat = 47.15  ! lake superior
!  testlon = 292.34 ; testlat = 66.86
!  testlon = 245.31 ; testlat = 61.72  ! great slave lake
!  testlon =  60.28 ; testlat = 67.99
!  testlon = 107.04 ; testlat = 53.05  ! lake baikal
!  testlon =  73.20 ; testlat = 69.03  ! ob river
!  testlon =  73.53 ; testlat = 68.07  ! ob river
!  testlon = 269.95 ; testlat =-44.96
!  testlon =   0.84 ; testlat =-15.63
!  testlon = 233.41 ; testlat = 43.47
!  testlon =  80.65 ; testlat =  0.13  ! c384
!  testlon =  80.78 ; testlat =  0.26  ! c192
!  testlon =  52.25 ; testlat = 38.27  ! lake c192

!  print '(a,2f8.2)','(get_testpt) set test point location',testlon,testlat

   return
   end subroutine get_testpt


   real function kd_par(deglon,deglat)
!  Son and Wang (2015): Diffuse attenuation coefficient of the photosynthetically available
!  radiation Kd(PAR) for global open ocean and coastal waters
! --- ipar = aa + bb *log(kd_par) at 1x1 deg resolution; kd_par: 0.03-0.3, ipar:1-9
! --- kd_par = exp((ipar-aa)/bb)

   implicit none
   real,intent(IN) :: deglon,deglat
   real, parameter:: aa=13.1830, bb=3.4744
   integer :: ipar(360,180),i,j,k
   character(len=6480) :: char(0:9)

   char(0)=	&
"999999999999999999999999999999999999999999999555566655544444433333344445555544444555555555555567999999999999999999999999999999555599999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555566655544444443333444445555544444555566555555557999999999999999999999999999999955559999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555566555544444444444444455555554555555666555555557999999999999999999999999999999955656999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555566655544444444444444455555555555556666555555555999999999999999999999999999999955555599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555566655544444444444445555555555555556666655555998999999999999999999999999999999955555599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555566655544444444444455555555555555556666665555689999999999999999999999999999999955555599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555566555544444444445555555666555555566666666655699999999999999999999999999999999955555599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555566555554444444455555566666655555666776666656799999999999999999999999999999999955555559999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995955666555555444455555556666666655556666776668666799999999999999999999999999999999955995559999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555566555555554555555566666776665556667777799987999999999999999999999999999999999955995859998999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666555555555555555666677776666666667887999999999999999999999999999999999999989955455599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555666555555555555566666777786666666667999999999999999999999999999999999999999655555455999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666555555555555666777799999976666679999999999999999999999999999999999999999445555457998999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666655555555566677889999999999779999999999999999999999999999999999999999999444495459996799999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666655555555666899999999999999999999999999999999999999999999999999999999999444495459656999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666655555557789999999999999999999999999999999999999999999999999999999999994444555599559999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666665555667999999999999999999999999999999999999999999999999999999999999944444449996589999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666666567899999999999999999999999999999999999999999999999999999999999999944444445595599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665666679999999999999999999999999999999999999999999999999999999999999999934444444555999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665566699999999999999999999999999999999999999999999999999999999999999999444444444498999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665566689999999999999999999999999999999999999999999999999999999999999999997444445999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665556679999999999999999999999999999999999999999999999999999999999999999999444459999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995566665555679999999999999999999999999999999999999999999999999999999999999999999444499999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666555679999999999999999999999999999999999999999999999999999999999999999994445478599999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995566666655579999999999999999999999999999999999999999999999999999999999999999994444455569999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666665579999999999999999999999999999999999999999999999999999999999999999994444455569999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666665559999999999999999999999999999999999999999999999999999999999999999964444457569999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666655557999999999999999999999999999999999999999999999999999999999999999944444499979999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666555555999999999999999999999999999999999999999999999999999999999999999944444499987799999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665555555599999999999999999999999997669999999999999999999999999999999999944444999997678999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665555555559999999999999999999999879999999999999999999999999999999999999994444699997667799999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666555555555699999999999999999999999999999999999999999999999999999999999994444499997666799999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666655555555559999999999999999999999999999999999999999999999999999999996974444999999766799999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666555555555555559999999999999999999999999999999999999999999999999999999944444999999767999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555555555544459999999999986679999999999999999999999999999999999994499994444599999777999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556665555555555544445696699999999999999999999999999999999999999999999444999989944699999777999999999999999999999999999999999999999999999"
   char(1)=	&
"999999999999999999999999999999999999999999999556665555555554444444444459999999999999999999999999999999999999999954459999999999999997777999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555555554444444444444999999999999999999999999999999999999965544999999999999999997777999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995556666555555544444444444444599999999999999999999999999999999999765559999999999999989997779999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666555555444444444444444499999999785589999999999999999999977855999999999999999999997779999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555554444444444444444449999894444445999999999999999999987569999999999999999999997799999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556665555554444444444444444444444444444455599999999999999997766999999999999999999999997999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556665555554444444544444444444444444444555558999999999999976999999999999999999999998999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556665555555444444555599445554444444444555555599999999998699999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555665555555444444569999999994444444444555555569999999986699999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555665555544444444599999999995444444444455555559999999997699999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555665555544444444469999999999444445444455555555999999997699999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555665555544444444445999999999644445444455555555699999996569999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999665665555544444444444334999999996445544455555555569999996569999999999999989999999999999789999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999655665555544444444444333334799999445555555555555555699997667999999999999877999999997898778999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999665665555544444444444433333337344455555555555555555566999666999999999989779999999977797779999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999655665555544444444444333333333444445555555555555555556677666999999999997799999999977777799999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999656655555544444444444333333333444444555545555555555555666667799999999777899999999977777799999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999656655555544444444444333333333444444555555555555555555566967799999999777999999999988999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999656555555544444444444333333333444444555555555555555555566666799999999978999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999956555555544443344444333333333444444555555455555555555556666699999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556555555544443334444333333333444444555555445555555555556666667999999967999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556555555544443333443333933333444444555554444555555555555666666799999669999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433333343333333344444444555554444455555555555666666979999699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995565555555544433333333333333354444445544444444445555555555669699967996699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995665555555544433333333333333334444454444444444444455555555666969996677699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999655555555544433333333333333333344454444444444444444555555566969996966699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555544333333333333333333344444444444444444444445555566999996666999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996955555555544333333333333333333344444444444444444444445555566999999996699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555444333333333333333333344444444444444444444444555556999999669699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333344444444444444444444444455555999969696699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333344444444444444444444444455559599959666799999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333344444444444444444444444455559999999967999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544333333333333333333344444444444444444444444445559999959699999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333444444444444444444444444445555999996999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333444444444444444444444444455555999967999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333444444445444444445554445555555999969999999999999999999999999999999999999999999999999999999999999999999999"
   char(2)=	&
"999999999999999999999999999999999999999999999555555555444333333333333333333444444445544444445555555555555999799999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995955555555444333333333333333333444444444444444645666555555599999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333444444444444444445555555556679999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555444333333333333333333444444444444444445555555668999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544333333332333333333444444444444444444555567999999999999999999999999999999999799999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544333333322333333333444444444444444444555569999999999999999999999999999999999599999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544333332232333333333444444444444444445555566999999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555544333322222333333333444444444444444444555666899999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544333322222333333333444444444444444444556999955565799999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544333322223333333333444444444444444445555995555555599999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433322223333333333444444444444444445555555555555579999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433333322333333333444444444444444444555555555555556999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433332222333333333444444444444444444555555555555555699999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433333322333333333444444444444444444555555555555555579999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433333333333333333444444444444444444455555555555555569999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433333323333333333444444444444444444445555555555555567999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554433333333333333333444444444444444444445555554444555556999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554433333333333333333444444444444444444445555544444455555999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544433333333333333333444444444444444444445555544444455556999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333333444444444444444444445555554444555556999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333334444444444444444444445555556555555559999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333334444444444444444444445566555555555699999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555544443333333333333334444444444444444444445555555555797999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555544443333333333333334444444444444444445559555555555999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333334444444444444444455595555555556999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555554443333333333333334444444444444445566996555555556999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333334444444444444459569966669998999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333334444444444444555999866799969999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554443333333333333344444444444445569999799986555699999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544443333333333333344444444444455999998999955555799999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555544444333333333333344444444444559999999999555556999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444333333333333444444444445599999999975556689999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999955555555544444444443343334444444444455999999775555559799999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444444444444444444444454456999986655545567999999999999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444444444444444444444444569677865555545556999999999669999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444444444444444444444444599666655655444555699999996667999999999999999999999999999999999999999999999999999999999999999999999"
   char(3)=	&
"999999999999999999999999999999999999999999999955555555544444444444444444444444444599666667765544445569999965976999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444444444444444444444444596566699995544444455665545999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444444444444444444444445596569999996544444455444444699999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555544444444444444444444444455696569999999644444444444454556999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554444444555444444444444455695569999999655444444444444455999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555554444455989654444444444455765569999999965444444444444455799999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555445579999995444444444455655679999999995544444444444455689999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555568669999999999844444444455955569999999999554444444444455589999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555599999999999999754444444455955567999999999954444444444455569999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555599999999999999754444444455965556669999999957444444444455569999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555599999999999999964444444455955655555595599555944444444555567999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555559999999999999965444444559755969995555555545596444444455557799999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555559999999999999965444444555955566994555456545555579999555499779999999998999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555559999999999999996565444555655579994955445755555559999854446966999999997998999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555559999999999999999965444455955676555854444565989795444434444466679999977997999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555559999999999999999987545465955655575954444559576694333333444456667799777777899999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555999999999999999998555549655555555654445999559933333333444445667799777778999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555555999999999999999999865545655555555555544999859333333333344445566797797899999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555699999999999999999966555955559555555443394333333333333344444556676999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555699999999999999999966555565555599665444333333333333333334445455666699999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555699999999999999999966655555557556999544333333333333333334444455556679997699999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555555699999999999999999976755555555555555444433333333333333333444455556667766669999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555555599999999999999999999955555559567554444433333333333333333344445599966666666999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555555599999999999999999999955595555695544444433333333333333333344445569996666666699999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555555599999999999999999999965555656899544444333333333333333333344444455996666666699999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555699999999999999999999865565569995444444333333333333333333334444445998666666699999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555556999999999999999999999665579969955444444494333333333333333334444445679666666699999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555599999999999999999999999665667689555444444433333333333333333334444455698666666669999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555555555569999999999999999998679766877699955444444433333333333333333334444455599966666666999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555697999999999999999996666666998999985444444434333333333333333334444555569966666666699999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555556999999999999999999975555666999999964444444433333233333333333334444455569966666666679999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556555556999999999999999999996556666799999954444444333332233333333333334444455559997666666667999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566655599999999999999999999976666667999999954444444333332223333333333334444445558999996766666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655699999999999999999999999997668999999554444444333332223333333333334444445555668996796666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655699999999999999999999999999978999999554444443333332223333333333334444445555566666999779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655999999999999999999999999997666999999554444443333332222333333333333444445555566667997799999999999999999999999999999999999999999999"
   char(4)=	&
"999999999999999999999999999999999999999999999666666699999999999999999999999944455699998554444443333322232333333333333444445555556666797779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996669965999999999999999999999986544445699965554444433333322222333333333333344444555566666787779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999699966999999999999999999996644444445999555544444443333322222223333333333344444555556666677779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996669966699999999999999999995333444459969555544444433333322222223333333333344444555556666677779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666676699999999999999999953433444459555755444444433334322222223333333333344444555555666667779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666699999999999999997553343444456555754554444433333322222222333333333344444555555666666679999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666565699999999999995733433444499655954744444433333322222222333333333334444555555666666679999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555556999999999666433433444454445955444444333353322222222333333333334444555555666666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555555799999996556333333444545444554444444333333322222222333333333334444555555666666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555555555698554443333333444944444444444444334333322222222333333333334444555555666666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666665555555544444444443333333444444444544444444333333322222222323333333334444555555666666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666655555555554444444443333333444444459444444444333333322222222222333333334444555555566666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655555555554444444433333333344444455444444444433333322222222223233333334444555555566666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555555555554444444433333333444444564444444444333333322222222232333333334444555555566666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555555555554444444333443333444445574444444444433333322222222222233333334444555555566666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555555555554444444334333333444445644444444444333333322222222222223333334444555555566666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999665555555555554444444333433334444446544444444444433333322222222222233333334444555555566666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555555555554444444333444434444499444444444444433333322222222222233333334444555555566666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999665555555555544444444333344444444444444444444444433333323222222222233333334444555555566666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555555555544444444333454444444444444444444444443333322222222222333333334444555555566666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555555555544444444334944444444444444444444444433333322222222222233333333444455555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999665555555555544444444439444444444444444444444444443333333222222222333333333444455555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996665555555555544444444464644444443444443444444444443333323222222222333333333444455555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999998665555555555544444444454444494343444443444444444443333333222222222333333333444555555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999965555555555544444444444444543333344433444444444443343332222222222333333333444555555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999975555555555544444444444454433333333334444444544443333332222222222333333333444555555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999996555555555544444444444444433333333334444444554444333332222222222333333333444455555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996799755555555544444444444444333333333334444445554444333332222222223333333333444455555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996699975555565544444444444334333333333334444445554444333332222222222333333333444455555556666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996668966555955444444444333333333333333334444455554444333332222222223333333334444455555555666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666669679755444444444333333333333333334444455554444333332222222223333333333444455555556666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666899997555444444433333443333333333334444455544444433332222222233333333333444455555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666699965555444444333333443333333333334444455544444433333222222233333333333444455555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666955554444444333333444843333333334444445544444433333222222233333333333444455555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666555665554444443333334444663333333344444455554444433333222222333333333333444455555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655555555444443333334444493333333344444455554444433333222222333333333333444445555555666669999999999999999999999999999999999999999999"
   char(5)=	&
"999999999999999999999999999999999999999999999666655555554444443333334444453333333344444455554444433333222233333333333333444445555555666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666555555544444443333334445433333333344444455554444433333222223333333333333444445555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666555554444444443333344444333333333344444455554444443333222233333333333333444445555555566666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666555544444444443333344444333333333344444455554444443333222333333333333333444455555555666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999665555544444444433333344444333333333344444455554444443333322333333333333333444445555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999565555444444444433333344444333333333344444455554444443333323333333333333333444445555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555444444444333333344443333333333344444455555444443333323333333333333333444445555555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996555555444444444333333344333333393333444444455555444443333333333333333333333444445555555565669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555554444444444333333333333333333333444444455555444443333333333333333333333344444555555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555554444444433333333333333333333333444444455555444443333333333333333353333444444555555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555554444433333333333333333333333333444444455555544443333333333333333333333344444555555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555554444333333333333333322233333333444444455555544444333333333333333333333344444555555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555554444333333333333333222222333333444444455555544443333333333333333333333344444555555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544444333333333333332222222233333444444455555544444333333333333333333333344444555555555666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333333333332222222233334444444455555544444333333333333333333333344444555555555666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333333333322222222333334444444455555544444333333333333334333333344444555555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333333333222222222233334444444455555544444333333333333333333333344444555555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333333332222222222333334444444455555544444333333333333333333333344444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333333222222222222333334444444455555544444333333333333333333333344444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333332222222222222333334444444455555544444333333333333333333333344444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555544443333333332222222222222333334444444455555554444433333333333333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555444443333333322222222222223333334444444455555554444433333333333333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555444443333333322222222222223333334444444455555554444433333333333333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555444443333333322222222222223333344444444455555554444433333333333333333333334444455555555666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555444443333333222222222222233333344444444455555544444433333333393333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555444443333333222222222222233333344444444455555544444443333333333333333333333444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444443333333222222222222233333444444444455555554444443333333333333333333334444445555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995554444443333332222222222223333333444444444455555554444443333333333333333333333444445555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555444443333332222222222223353333444444444455555554444443333333333333333333333444445555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995554444443333332222222222233333334444444444455555555444444333333333333333333333444445555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444443333332222222322233633334444444444455555554444444333333333333333333333444445555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444443333332222222222233333334444444444455555555444444333333333333333333333444445555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995554444443333322222222222233333344444444444455555555444444433333333333333333333444445555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444443333322222222222233334444444444444455555555444444433333333333333333333444445555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444443333322222222222223555344444444444455555555444444443333333333333333333444444555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444443333322223222222223343344444444444455555555444444443333333333333333333444444555555569999999999999999999999999999999999999999999"
   char(6)=	&
"999999999999999999999999999999999999999999999444444443333322222222222323333444444544444455555555444444443333333333333333333344444555555659999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999544444443333222222222222223333444445544444455555555444444444333333333333333333444444555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222223333444445544444455555555444444444333333333333333333444444555555566999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222233333444455544444455555555444444444333333333333333333444444555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222223333444455444444455555555444444444333333333333333333444444455555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222223334444454444444455555555544444444433333333333333333344444455555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222223334444444444444455555555544444444433333333333333333444444455555556999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222224334444444444444455555555544444444433333333333333333444444455555556999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222223233334444444444444455555555544444444433333333333333333344444455555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222233334444444444444455555555544444444433333333333333333344444455555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444443333222222222222233344444444444444455555555544444444443333333333333333344444455555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444443333222222222222233344444444444444455555555544444444443333333333333333344444455555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222233344444444444444455555555544444444443333333333333333344444555555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222233344444444444444455555555544444444443333333333333333444445555555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444443333222222222222233344444444444444455555555544444444444333333333333333444455555555555999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444443333222222222222233344444444444444455555555544444444444333333333333333444555555555556999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444443333222222222222233344444444444444455555555544444444444333333333333334444555556665556999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222222222222233344444444444444555555555544444444444333333333333344445555666666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222222222222333344444444444444555555555444444444444433333333333344455566666666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222222222222333344444444444444555555555444444444444433333333333444555666687877789999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444444333222222222222333444444444444444555555555444444444444433333333333445556677999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222222111222333444444444444444555555555444444444444433333333333445566799999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222221111222333444444444444444555555555444444444444433333333334455667999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222221111222233444444444444444555555555444444444444443333333444455789999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222221111222233444444444444444555555555444444444444443333334444456799999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222221111222333444444444444444555555555444444444444443333344444456699999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222221111222333444444444444444555555555444444444444443333444444556999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222221111222333444444444444444555555555444444444444444344444555599999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433222222111222333444444444444445555555555444444444444444444445557999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433222221111222333444444444444445555555555444444444444444444445969799999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433222222222222333444444444444455555555555444444444444444444457998799999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433222222222222333444444444444455555555555444444444444444444569978999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444333222222222222333444444444444555555555554444455444444444455999699999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433322222222222333444444444444555555555554444555444444444556666999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433322222222222334444444444445555555555554444555544444444559667999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433322222222222334444444444445555555555554445555544444445555699999999999999999999999999999999999999999999999999999999999999999"
   char(7)=	&
"999999999999999999999999999999999999999999999444444444433322222222222334444444444455555555555554445555544444555555999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444433322222222222334444444444455555555555554445555544444555659999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444443322222222222334444444444455555555555554455555554445568799999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444443322222222222334444444444455555555555554455555554445599999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444443332222222222334444444444555555555555554455555554445999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444443332222222222334444444444555555555555554455555554455999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444444443332222222223334444444444555555555555554455555555556999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444444443332222222223334444444444555555555555554455555555559999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444443332222222223334444444445555555555555555555555555569999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444333222222223334444444455555555555555555555555555599999999999999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444333222222233344444454555555555555555555555555555599999956699999999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444333222222233344444555555555555555555555555555565699995444456699999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444333322222333344444555555555555555555555555555666699954444445579999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444433322222333344445555555555555555555555555655566699655444444569999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444433333233333344445555555555555556655555555655555699855444444569999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444433333333333344445555555555555666655555556665556999955544444569999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444444444443333333333444455555555555555667655555556666559999965544444569999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444444444443333333333444455555555555555667655555556666569999999644444579999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444443333333333444455555555555555566555555556666599999999754444569999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444443333333334444555555555555555566555555556666699999999954444558999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444443333333334444555555555555555666555555555666699945599954444556999999999999669999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444444333333334444555555555555555666555555555566999944444544444556999999999999776799999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999444444444444444433333344444555555555555666666555555555599999944444444444557999999999999999699999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999994444444444454444444334444444555555556666666666555555555799999944444464444557999999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999554444445555444444444444445555555566666666666555555559995999754444494445569999999999999999699999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555454555555444444444444445555556666666777666655555559944445554444694657999999999999999986699999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555544444444444455555666666677777766655555569544444354444494689999999999999999997699999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555544444444444455555666666779999998665566699544444454444494499996669999999999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555555555555554444444444455556666667999999999965566667944444444333795444444456999999999979999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999655555555555555555444444555566667778999999999999556668644444444434997784444445799999999979999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666666655555555555555555566678999999999999999997699444444444537977547344445699999999979999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666666666555555555555555667899999999999999999999999544444444339477443344444589999999989999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999776666666666655555555555555667999999999999999999999999995444444339346533334444459799999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999998877776677766666655555555556679999999999999999999999999996544444534343333333444445667999999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999878998998776666666655555556689999999999999999999999999996544444933333333333444444566799999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999999999999977666666666555566699999999999999999999999999999544444993433333333344444456679999999999999999999999999999999999999999999999999"
   char(8)=	&
"999999999999999999999999999999999999999999999999989999999899776666666666999999999999999999999999999999944444993933333333344444456679999999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999999999999999999999999976778999999999999999999999999999997644444993333333333344444456679899999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999999999999999999999999999999997999999999999999999999999999554444934333333333344444456677779999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999995555444933333333333344444456677779999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999995555444433333333333344444455668779999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997999999999999999999999999999999999999999999999999999999995555444633333333333344444455567778999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997999999999999999999999999999999999999999999999999999999955555444533333333333344444455566799999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997778799999999999999999999999999999999999999999999999999976655444533333333333344444455556679999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997778899999999999999999999999999999999999999999999999999997655445433333333333334444455556677999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997777899999999999999999999999999999999999999999999999999997555444333333333333334444455556667999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997777778999999999999999999999999999999999999999999999999997655554333333322333334444455556667999999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999997787777999999999999999999999999999999999999999999999999966554444333333333333334444455556667799999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996678777999999999999999999999999999999999999999999999999555554444333333333333334444455555667779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996667777999999999999999999999999999999999999999999999995555554443333333333333334444455555666779999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666877799999999999999999999999999999999999999999999965555554443333333333333334444455555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666677789999999999999999999999999999999999999999999755555554443333333333333334444445555666667999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666667778999999999999999999999999999999999999999999765555544443333333333333334444455555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666777899999999999999999999999999999999999999999765555544433333333333333334444445555666668999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666666667799999999999999999999999999999999999999999765555444333333333333333334444445555666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555679999999999999999999999999999999999999998765554443333333333333333334444445555666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555557899999999999999999999999999999999999998655554443333333333333333334444445555666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555555679999999999999999999999999999999999995555544433333333333333333334444445556666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666666555555556799999999999999999999999999999999854555544433333333333333333334444455556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555555445579799999999999999999999999999996444555544433333333333333333334444555556666666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666666555555444455569999999999999999999999999864444555444433333333333333333334444555556666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555555444444556999999999999999999999998544444555444433333333333333333334444555555566669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555554444444456999999999999999999999995544444555444443333333333333333334444555555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666666555554444444444699999999999999999999995444444554444443333333333333333334444555555556669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555554444444444699999999999999999999954444444544444433333333333333333334444555555556669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666666555554444444444599999999999999999999744444444444444433333333333333333334444555555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555554444444444569999999999999999999644444444444444433333333333333333334444555555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555544444434444456999999999999999999644444444444444443333333333333333334444555555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666665555544444433444444466999999999999999544444444454444443333333333333333334444455555555559999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666655555544444433334344445543349999999995444444444444444443333333333333333334444555555555659999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655555544444433333334444433334999999974444444444444444443333333333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666655555544444433333333433333333399999944444444444444444444333333333333333334444455555555569999999999999999999999999999999999999999999"
   char(9)=	&
"999999999999999999999999999999999999999999999666666555544444433333333333333333337999944444444444444444444333333333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555444444333333333333333333333343344444444444444444444433333333333333333444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555444444333333333333333333333333444444444444444444444444433333333333334444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555444443333333333333332233333333444444444444444444444444444433333333334444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666665555444443333333333333222223333333444445444444444444444444444443333333334444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655555444443333333333332222222333333444445544444444444444444444444333333334444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555444433333333333332222222333333444445544444444444444444444444433333334444455555555566999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655555444433333333333322222222333334444445544444444444444444444444433333344444456555555566999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666655555444433333333333222222222333334444445544444444444444444444444443333344444455555555569999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666665555444433333333332222222222333334444445544444544444454444444444443333344444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666655555544433333333332222222222333334444445544444554445445555444444444333344444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999996666665555544433333333322222222223333334444455554444554455555555554444444333344444455555555666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555544433333333322222222223333344444455544444555555555555555544444433344444445555555666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999666665555444333333333322222222223333344444455544444555556655555555554444443444444445555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665555444333333333222222222233333344444555554455555557755556665555444444444444445555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665555444333333333222222222333333344445555554455555566566666666655444444444444445555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566665555444333333332222222233333333344445555554455555555666666677665444444444444445555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999566666555444333333333222222333333333444445555544455555566667667777765544444444444445555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555444333333333222223333333333444455555554455555567778877779876554444444444445555555666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555444333333333222333333333333444455555554455556679999998899998655454445444445555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555444333333332223333333333333444455555554455556789999999999999855544444444445555556666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555444333333333223333333333333444455555554455566799999999999999996554444444445555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555444333333333323333333333333444455555554555577999999999999999999755444444445555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666555444333333333333333333333333444455555554555799999999999999999999965544444455555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999556666655444333333333333333333333334444455555554555699999999999999999999995554444555555555669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666655444333333333333333344333334444455555554556999999999999999999999996654444555555656669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555666655444333333333333333444443334444455555554559999999999999999999999999996444557666666669999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555666655444333333333333333444444344444455555554569999999999999999999999999999644599999976566999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999995555666655444333333333333334444444444444455555554599999999999999999999999999999944599999996666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555565655444333333333333344444444444444555555554579999999999999999999999999999995799999995666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666655544433333333333344444444444444555555544569999999999999999999999999999999999999995666999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666655544433333333333344444444444445555555554569999999999999999999999999999995699999965667999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666655544433333333333444444444444445555555555569999999999999999999999999999996699999986667999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555666655544443333333334444444444444455555555555569999999999999999999999999999995699999966679999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999555566655544443333333334444445444444455555555555569999999999999999999999999999995599999997899999999999999999999999999999999999999999999&
999999999999999999999999999999999999999999999955566655544444333333334444555444444555555555555569999999999999999999999999999999559999999999999999999999999999999999999999999999999999"

! (1,1) => (0.5E,89.5S); (360,180) => (359.E,89.5N)
   do k = 0,9
     read(char(k),'(36(180i1))') ((ipar(i,:)),i=(k*36+1),(k+1)*36)
   end do 

   i=max(1.,min(360.,deglon+1.))
   j=max(-89.5,min(89.5,deglat))+91.

   kd_par = exp((float(ipar(i,j))-aa)/bb) * 7.

   return
   end function kd_par
