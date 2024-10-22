!>\file skinsst.f90
!! This file contains Rainer's skin temperature scheme.

module state_eqn
   implicit none
! --- coefficients for sigma-0 (based on Brydon & Sun fit, JGR 1999)
   real, parameter, dimension(7) :: coef = (/				&
   -1.36471E-01, 4.68181E-02, 8.07004E-01,-7.45353E-03,-2.94418E-03,	&
     3.43570E-05, 3.48658E-05 /)
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
   end

   real function dsigds(t,s)
! --- saline contraction coefficient
   real, intent(in) :: t, s
!  dsigds = coef(3)+t*(coef(5)+t*coef(7))
   dsigds = t*(t*coef(7)+coef(5))+coef(3)		 ! alt.grouping
   return
   end
end module state_eqn


module skinsst
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

   subroutine skinsst_run(					&
!   c_0, c_d, w_0, w_d, d_conv, ifd,				&
    im,			& ! horiz. loop extent				in
    iter,		& ! ccpp loop counter				in
    wet,		& ! .true. at ocean & lake points		in
    oceanfrac,		& ! cell portion covered by ocean		in
    timestep,		& ! model timestep				in
    xlon, xlat,		& ! longitude, latitude				in
    sfcemis,		& ! sea surface emissivity			in
    ulwflx,             & ! upwelling LW flux				inout
    dlwflx,		& ! absorbed downwelling LW flux		in
    sfcnsw,		& ! net SW flux, pos.down			in
    tsfco,		& ! ocean/lake top layer temperature		inout
    tref,		& ! foundation temperature			in
    psfc,		& ! surface pressure				in
    wind,	 	& ! atm. mid-layer 1 wind			in
    stress,		& ! wind stress (N/m^2)				in
    plyr1, 		& ! atm. mid-layer 1 pressure			in
    tlyr1, 		& ! atm. mid-layer 1 temperature		in
    qlyr1,		& ! atm. mid-layer 1 humidity			in
    ulyr1,		& ! atm. mid-layer 1 zonal wind			in
    vlyr1,		& ! atm. mid-layer 1 meridional wind		in
    compres,		& ! adiabat.compression factor, layer 1 => sfc	in
    cm,			& ! drag coeff for momentum			in
    ch,			& ! drag coeff for heat and moisture		in
    hvap,		& ! latent heat of evaporation			in
    cp,			& ! specif.heat of air				in
    rd,			& ! gas constant of dry air			in
    eps,		& ! ratio of gas constants, rd/rv		in
    sbc,		& ! stefan-boltzmann constant			in
    tskin,		& ! skin temp					inout
    xzts,		& ! holding place for tskin			inout
    xt,			& ! lake mixed layer temperature		inout
    xv,			& ! seawifs extinction coefficient		inout
    xz,			& ! lake ice thickness				inout
    xs,			& ! not used, held in reserve			inout
    xu,			& ! previous lake ice surf. temp		inout
    zm,			& ! previous lake sfc heat flux			inout
    qsat,		& ! saturation specif. humidity			out
    z_c,		& ! sub-layer cooling thickness			out
    dt_warm,		& ! warm-layer surface warming amount		out
    dt_cool,		& ! skin layer cooling amount			inout
    evap,		& ! kinematic latent heat flux, pos.up		out
    hflx,		& ! kinematic sensible heat flux, pos.up	out
    ep,			& ! potential latent heat flux, pos.up		out
    cmm,		& ! momentum exchange coeff			out
    chh,		& ! thermal exchange coeff			out
    lseaspray,		& ! sea spray flag				in
    fm,			& ! Monin-Obukhov function at surface		in
    fm10,		& ! Monin-Obukhov function at 10m		in
    errmsg, errflg)

   use machine , only : kind_phys
   use state_eqn
   use funcphys, only : fpvs		! vapor pressure

   implicit none

! --- input:
   integer, intent(in) :: im,iter
   logical, dimension(:), intent(in) :: wet
   logical, intent(in) :: lseaspray
   real (kind=kind_phys), dimension(:), intent(in) :: xlon,xlat,	&
      sfcemis, dlwflx, sfcnsw, wind, psfc, plyr1, tlyr1, qlyr1,		&
      ulyr1, vlyr1, cm, ch, compres, stress, fm, fm10,oceanfrac,tref
   real (kind=kind_phys), intent(in) :: timestep, hvap, cp, rd, eps,	&
      sbc

! --- inout:
   real (kind=kind_phys), dimension(:), intent(inout) :: ulwflx,	&
      tsfco, tskin, dt_cool, xzts, xs, xt, xu, xv, xz, zm
!     c_0, c_d, w_0, w_d, d_conv, ifd
!     temice, thkice

! --- output:
   real (kind=kind_phys), dimension(:), intent(out) :: evap, hflx,	&
      ep, qsat, cmm, chh, dt_warm, z_c
   character(len=*),intent(out) :: errmsg
   integer,         intent(out) :: errflg

! --- locals:
   integer :: i, n, loop
   real :: alon, alat, virt, rho_air, rho_wat, pvap, tsq, piston, vel,	&
     nonsol,			& ! sum of nonsolar air-sea fluxes (pos.up)
     spcifh = 3990.,		& ! seawater specific heat
     grav  = 9.806,		& ! gravity
     sss = 34.7			  ! sea surface salinity
   integer,parameter :: itmax = 5		! regula falsi iterations
   real :: rnl_ts, hs_ts, rf_ts, alpha, beta, rch, ustar, 		&
      hist(0:itmax) = 0., x1, x2, x3, y1, y2, y3, dif1, dif2, dif3

!  variables for sea spray effect
   real (kind=kind_phys) :: f10m, u10m, v10m, ws10, ru10, qss1,		&
                            bb1, hflxs, evaps, ptem, tem
   real (kind=kind_phys), parameter :: alps=0.75,bets=0.75,gams=0.15,	&
                            ws10cr=30., conlf=7.2e-9, consf=6.4e-8

   logical :: doprint, details, frstrip
   real(kind=kind_phys) :: xtinct, frz=273.15, small=.05, totflx,	&
      oldflx, testlon, testlat
   external xtinct
   real,parameter :: rad2deg = 57.2957795
   doprint(alon,alat)=abs(testlon-alon).lt.small .and.			&
                      abs(testlat-alat).lt.small

! --- piston velocity: at 0 m/s, molecular diffusion only.
! ---                  at 8 m/s, destroy warm layer during 1800 sec time step
! ---                  based on 1 m thickness scale.
!  piston(vel)= 1.4e-7 + 5.5555e-4*(vel/8.)**2	! quadratic interpolation
   piston(vel)= 1.4e-7 + 5.5555e-4*(vel/8.)	! linear interpolation

   if (iter.gt.1) return

   call get_testpt(testlon,testlat)
! --- temporary:
!  print '(a,2f8.2,l5)','entering skinsst_run, testpt =',testlon,testlaa

   do i = 1,im

    details=.false.

! --- check which arrays inherited from NSST are touched outside skinsst.
!    if (xzts(i).ne.-.03125) details=.true.
     if (xs(i).ne.-.03125) details=.true.
!    if (xt(i).ne.-.03125) details=.true.
!    if (xu(i).ne.-.03125) details=.true.
!    if (xv(i).ne.-.03125) details=.true.
!    if (xz(i).ne.-.03125) details=.true.
!    if (zm(i).ne.-.03125) details=.true.
!   if (c_0(i).ne.-.03125) details=.true.
!   if (c_d(i).ne.-.03125) details=.true.
!   if (w_0(i).ne.-.03125) details=.true.
!   if (w_d(i).ne.-.03125) details=.true.
!if (d_conv(i).ne.-.03125) details=.true.
!   if (ifd(i).ne.-.03125) details=.true.
!   if (details) then
!    alon=xlon(i)*rad2deg
!    alat=xlat(i)*rad2deg
!    if (doprint(alon,alat))						&
!    print '(a,2f8.2/13f6.3)','problem at lon,lat=',alon,alat,		& 
!     xzts(i), xs(i), xt(i), xu(i), xv(i), xz(i), zm(i), c_0(i),	&
!     c_d(i), w_0(i), w_d(i), d_conv(i), ifd(i)
!   end if

    if (wet(i)) then

     alon=xlon(i)*rad2deg
     alat=xlat(i)*rad2deg

! --- temporary:
!    if (xzts(i).eq.0.)							&
!     print '(a,3f8.3,l5)','now at lon,lat',alon,alat,			&
!     oceanfrac(i),doprint(alon,alat)

     if (doprint(alon,alat)) then
      print 97,'entering skinsst_run   lon,lat=',alon,alat,		&
!     'xs',xs(i),'xz',xz(i),'c_0',c_0(i),'c_d',c_d(i),			&
!     'w_0',w_0(i),'w_d',w_d(i),'d_conv',d_conv(i),'ifd',ifd(i),	&
!     'temwat',xt(i)-frz,		& ! lake water temperature
      'xtinct',xv(i),			& ! extinction coefficient
!     'thkice',xz(i),			& ! lake ice thickness
      'ocnfrac',oceanfrac(i),		& ! ocean fraction
      'stress',stress(i),		& ! wind stress (N/m^2)
!     'sfcemis',sfcemis(i),		& ! sfc emissivity
      'wind',wind(i),			& ! surface wind
      'pstonE3',piston(wind(i))*1.e3,	& ! piston velocity
      'sfcnsw',sfcnsw(i),		& ! total sky net SW flx into ocean
      'dlwflx',dlwflx(i),		& ! absorbed downwelling LW flux
      'ulwflx',ulwflx(i),		& ! upwelling LW flux
      'psfc',psfc(i)*.01,		& ! surface pressure (mb)
      'plyr1',plyr1(i)*.01,		& ! atm.layer 1 presure
      'tlyr1',tlyr1(i)-frz,		& ! atm.layer 1 air temp
      'qlyr1',qlyr1(i)*1.e3,		& ! atm.layer 1 humidity (g/kg)
!     'sigma_t',sig(tlyr1(i)-frz,sss),	& ! sea water density - 1000
!     'compres',compres(i),		& ! midlyr-to-sfc adiab.compression
      'xzts',xzts(i)-frz,		& ! previous tskin
      'dcoolE2',dt_cool(i)*100.,	& ! previous dtcool
!     'tref',tref(i)-frz,		& ! foundation temp
      'tsfco',tsfco(i)-frz		  ! ocean top layer temperature
      print '(5(a13,"=",l2))','lseaspray',lseaspray
!     if (oceanfrac(i).eq.0.) print '(2f7.2,a)',alon,alat,' is lake point'
     end if
 99  format (/a,2f7.2/(5(a8,"=",f7.2)))
 98  format (/a,2f7.2/(4(a8,"=",es11.4)))
 97  format (/a,2f7.2/(4(a8,"=",f11.6)))

     if (xzts(i).ne.0.) then
      if (xzts (i)-frz.lt.-40. .or. xzts (i)-frz.gt.40.)		&
       print '(a,4f8.2)','questionable xzts at',alon,alat,		&
       xzts (i)-frz,oceanfrac(i)
     end if

     virt = tlyr1(i) * (1. + (eps-1.)*qlyr1(i))
     rho_air = plyr1(i) / (rd*virt)
     rch = rho_air * cp * ch(i) * wind(i)		! W/m^2/deg
     cmm(i) = cm(i) * wind(i)
     chh(i) = rho_air * ch(i) * wind(i)
     ep(i) = 0.
     rho_wat = 1000. + sig(tsfco(i)-frz,sss)
     alpha = -dsigdt(tsfco(i)-frz,sss)/rho_wat
     beta  =  dsigds(tsfco(i)-frz,sss)/rho_wat
     ustar = sqrt(stress(i)/rho_air)			! air friction velocity

     if (xzts(i).eq.0.) then		! use xzts=0 as indicator for t=0
      frstrip = .true.
      dt_cool(i) = 0.
      tskin(i) = tsfco(i)
      xt(i) = tsfco(i)			! lake ftemp
      xu(i) = tsfco(i)			! previous lake ftemp
      xv(i) = xtinct(alon,alat)		! seawifs extinction coeffient
      xz(i) = 0.			! lake ice thickness
      zm(i) = 0.			! old heat flux
     else
      frstrip = .false.
      tskin(i) = xzts(i)			! tskin from previous time step
     end if

!     details = doprint(alon,alat)
      details = .false.


     if (oceanfrac(i).gt.0.) then

! --- apply warm layer correction
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- tskin from last time step has been saved in xzts
! --- lake variables (water temp, thickness) are saved in xt,xz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (tskin(i)-frz.lt.-25.)						&
       print '(a,2f8.2,es13.3)', 'excessively cold tskin at lon,lat',	&
       alon,alat,tskin(i)-frz

! --- bypass warm layer calculation if SW flux is up or tskin is below freezing
      if (sfcnsw(i).lt.0. .or. tskin(i).lt.frz) then
       tskin(i)=tsfco(i)
      else 					! SW flux > 0

! --- evaluate warm-layer increment during current time step

       dt_warm(i) = sfcnsw(i) * timestep * xv(i)/(rho_wat * spcifh)
! --- note: dt_warm is cumulative.
       tskin(i) = tskin(i) + dt_warm(i)					&
! --- subtract previous   dt_cool (since dt_cool is not cumulative)
         + dt_cool(i)				! sign convention: dt_cool > 0
! --- cooling by heat diffusion
       tskin(i) = tskin(i) + (tsfco(i)-tskin(i))			&
          * min(1.,timestep*piston(wind(i)))
      end if					! SW flux > 0

! --- save (tskin - top layer T) for diagnostic purposes
      dt_warm(i) = tskin(i) - tsfco(i)

! --- start cool-skin iteration, using REGULA FALSI (aiming for x_n = y_n)
! --- x1,x2,x3,y1,y2,y3 are consecutive dt_cool approximations.

      x1 = -.5
      x2 = +.5
      call surflx(nonsol, tskin(i)+x1, tlyr1(i)*compres(i), qlyr1(i),	&
           psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
           sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, details)

      call coolskin(ustar, nonsol, sfcnsw(i), evap(i), sss, alpha,	&
                     beta, rho_wat, rho_air, tskin(i)+x1, grav, hvap, 	&
                     y1, z_c(i), alon, alat, details)

      dif1 = y1 - x1
      if (y1.ne.0.) then
      
       do loop = 1,itmax

        call surflx(nonsol, tskin(i)+x2, tlyr1(i)*compres(i), qlyr1(i),	&
            psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
            sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, details)

        call coolskin(ustar, nonsol, sfcnsw(i), evap(i), sss, alpha,	&
                      beta, rho_wat, rho_air, tskin(i)+x2, grav, hvap, 	&
                      y2, z_c(i), alon, alat, details)
        dif2 = y2 - x2

        if (details) print '(a,3es11.3,i7)','(skinsst)  x1,y1,y1-x1 =',	&
          x1,y1,dif1,loop
        if (details) print '(a,3es11.3,i7)','(skinsst)  x2,y2,y2-x2 =',	&
          x2,y2,dif2,loop

        x3 = (x1*dif2-x2*dif1)/(dif2-dif1)		! regula falsi

        if (abs(dif2).gt.1.e-5) then

         if (abs(dif1).gt.abs(dif2)) then
          x1 = x2
          y1 = y2
          dif1 = dif2
         end if
         x2 = x3

        else				! we have convergence
         dt_cool(i) = x3 
         exit				! all done
        end if				! convergence
        hist(loop) = y2

       end do		! iteration loop

       if (loop.eq.itmax) then 
        if (abs(hist(loop)).gt..5) then
         print '(a,3f8.2/(11f7.2))','tskin not converging at lon,lat',	&
            alon,alat,hist(loop),(hist(n),n=1,loop)
           end if
       end if
       tskin(i) = tskin(i)-dt_cool(i)	! apply cold-skin correction
      end if				! y1 nonzero

     else				! oceanfrac = 0 => call sea ice model

      call surflx(nonsol, tskin(i), tlyr1(i)*compres(i), qlyr1(i),	&
          psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
          sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, details)

      if (.not.frstrip) then		! skip on 1st time step

! --- use rudimentary energy loan lake model 'enloan'.
! --- arguments: time step, air-sea heat flux, ice thickness ('xz'),
! --- water temp ('xt'),ice sfc.temp ('tskin')
       totflx = sfcnsw(i) - nonsol	! pos.down

! --- average totflx over 2 time steps to suppress comput.mode in enloan
       oldflx = zm(i)
       zm(i) = totflx
       totflx = .5*(totflx +oldflx)

       call enloan(timestep, totflx, xz(i), tskin(i), xu(i), xt(i),	 	&
         alon, alat, doprint(alon,alat))

       tsfco(i) = xt(i)
       if (xz(i).gt.0.) evap(i) = 0.
        
        
      end if

     end if				! oceanfrac zero or nonzero

! --- save tskin for next call to skinsst

     xzts(i) = tskin(i)			! for use at next time step

! --- according to  GFS_surface_composites_inter.F90,
! --- dlwflx is the absorbed portion of downwelling LW flux.
! --- hence, the total downwelling flux is dlwflx/sfcemis
! --- and the reflected part is (1-sfcemis)*dlwflx/sfcemis

     ulwflx(i) = ulwflx(i) + dlwflx(i)*(1.-sfcemis(i))/sfcemis(i)

     if (lseaspray) then
       f10m = fm10(i) / fm(i)
       u10m = f10m * ulyr1(i)
       v10m = f10m * vlyr1(i)
       ws10 = sqrt(u10m*u10m + v10m*v10m)
       ws10 = max(ws10,1.)
       ws10 = min(ws10,ws10cr)
       tem = .015 * ws10 * ws10
       ru10 = 1. - .087 * log(10./tem)
       qss1 = fpvs(tlyr1(i))
       qss1 = eps * qss1 / (plyr1(i) + (eps-1.) * qss1)
       tem = rd * cp * tlyr1(i) * tlyr1(i)
       tem = 1. + eps * hvap * hvap * qss1 / tem
       bb1 = 1. / tem
       evaps = conlf * (ws10**5.4) * ru10 * bb1
       evaps = evaps * rho_air * hvap * (qss1 - qlyr1(i))
       evap(i) = evap(i) + alps * evaps
       hflxs = consf * (ws10**3.4) * ru10
       hflxs = hflxs * rho_air * cp * (tskin(i) - tlyr1(i))
       ptem = alps - gams
       hflx(i) = hflx(i) + bets * hflxs - ptem * evaps
     endif

     if (doprint(alon,alat)) then
      print 97,'exiting skinsst_run   lon,lat=',alon,alat,		&
      'virt',virt-frz,			& ! virtual temp
      'rho_air',rho_air,		& ! air density
      'pvap',pvap,			& ! satur. vapor pressure (mb)
      'qsat',qsat(i),		 	& ! satur. specif.humidity
      'hflx',hflx(i),			& ! sensible heat flux
      'evap',evap(i),			& ! latent heat flux
      'nonsol',nonsol,			& ! net non-solar surface flux
      'sfcnsw',sfcnsw(i),		& ! net solar surface flux
      'ulwflx',ulwflx(i),		& ! upwelling LW flux
      'dwarmE2',dt_warm(i)*100.,	& ! temperature increment due to SW
      'dcoolE2',dt_cool(i)*100.,	& ! cool-skin temperature correction
      'tskin',tskin(i)-frz,		& ! skin temperature
      'tsfco',tsfco(i)-frz 		  ! ocean top layer temperature
      if (oceanfrac(i).eq.0.) print '(2f7.2,a)',alon,alat,' is lake point'
     end if

! --- convert fluxes from W/m^2 to "kinematic" i.e. velocity x fluxed variable
     hflx(i) = hflx(i)/(rho_air * cp)				! deg m/sec
     evap(i) = evap(i)/(rho_air * hvap)				! m/sec

    end if				! wet

! --- check which arrays inherited from NSST are touched outside skinsst.
!    xzts(i) = -.03125
     xs(i) = -.03125
!    xt(i) = -.03125
!    xu(i) = -.03125
!    xv(i) = -.03125
!    xz(i) = -.03125
!    zm(i) = -.03125
!   c_0(i) = -.03125
!   c_d(i) = -.03125
!   w_0(i) = -.03125
!   w_d(i) = -.03125
!d_conv(i) = -.03125
!   ifd(i) = -.03125

   end do		! im loop

   return
   end subroutine skinsst_run


  subroutine coolskin(ustar_a,f_nsol,f_sol_0,evap,sss,alpha,beta,	&
                      rho_w,rho_a,ts,grav,latnt,deltat_c,z_c,alon,alat,doprint)

! upper ocean cool-skin parameterizaion, Fairall et al, 1996.
! code extracted from NSST package

! input:
! ustar_a : atmosphreic friction velocity at the air-sea interface (m/s)
! f_nsol  : the "nonsolar" part of the surface heat flux (w/m^s)
! f_sol_0 : solar radiation at the ocean surface (w/m^2)
! evap    : latent heat flux (w/m^2)
! sss     : ocean upper mixed layer salinity (ppu)
! alpha   : thermal expansion coefficient
! beta    : saline contraction coefficient
! rho_w   : oceanic density
! rho_a   : atmospheric density
! ts      : oceanic surface temperature
! grav    : gravity 
! latnt   : latent heat of evaporation 

! output:
! deltat_c: cool-skin temperature correction (>0 means cooling)
! z_c     : molecular sublayer (cool-skin) thickness (m)

  use machine , only : kind_phys
  use module_nst_parameters, only: kw => tc_w,visw,cp_w, z_c_max,	&
    z_c_ini,ustar_a_min
  implicit none
  logical,intent(in) :: doprint
  real(kind=kind_phys), intent(in) :: ustar_a,f_nsol,f_sol_0,evap,	&
     sss,alpha,beta,rho_w,rho_a,ts,grav,latnt,alon,alat
  real(kind=kind_phys), intent(out):: deltat_c,z_c
! local variables:
  real(kind=kind_phys), parameter :: frz=273.15
  real(kind=kind_phys) :: xi,hb,ustar1_a,bigc,deltaf,fxp

  if (doprint) print 98,'entering coolskin   lon,lat=',alon,alat,	&
     'ustar_a',ustar_a,		&	! atmospheric friction velocity (m/s)
     'f_nsol',f_nsol,		&	! "nonsolar" part of the surface heat flux
     'f_sol_0',f_sol_0,		&	! solar radiation at the ocean surface
     'evap',evap,		&	! latent heat flux (w/m^2)
     'rho_w',rho_w,		&	! sea water density
     'rho_a',rho_a,		&	! air density
     'ts',ts-frz,		&	! ocean surface temperature
     'zc_ini',z_c_ini*1.e3		! sublayer thickness (mm)

 99 format (/a,2f7.2/(5(a8,"=",f7.2)))
 98 format (/a,2f7.2/(4(a8,"=",es11.4)))

  z_c = z_c_ini                 ! initial guess

  ustar1_a = max(ustar_a,ustar_a_min)

  call sw_rad_skin(z_c,fxp)
  deltaf = f_sol_0 * fxp

  hb   = alpha * (f_nsol-deltaf) + beta * sss * cp_w * evap/latnt
  bigc = 16. * grav * cp_w * (rho_w * visw)**3 / (rho_a * kw)**2

  if ( hb > 0 ) then
    xi = 6./(1+(bigc * hb/ustar1_a**4)**0.75)**0.3333333
  else
    xi = 6.0
  endif
  z_c = min(z_c_max,xi * visw/(sqrt(rho_a/rho_w) * ustar1_a ))

  call sw_rad_skin(z_c,fxp)

  deltaf = f_sol_0 * fxp
  deltaf = f_nsol - deltaf
  if ( deltaf > 0 ) then
    deltat_c =  deltaf * z_c / kw
  else
    deltat_c = 0.
    z_c      = 0.
  endif

  if (doprint) print 98,'exiting coolskin   lon,lat=',alon,alat,	&
     'fxp',fxp,			&
     'deltaf',deltaf,		&
     'hb',hb,			&
     'bigc',bigc,		&
     'xi',xi,			&
     'delt_c',deltat_c,		&	! skin layer temperature correction
     'z_c',z_c				! skin layer thickness

  return
  end subroutine coolskin


  subroutine sw_rad_skin(z,fxp)
  ! original name: elemental subroutine sw_ohlmann_v1(z,fxp)
  ! fraction of the solar radiation absorbed by the ocean at the depth z
  
  ! input:
  ! z: depth (m)
 
  ! output:
  ! fxp: fraction of solar radiation absorbed by the ocean at depth z (w/m^2)

  use machine , only : kind_phys
  implicit none
  real(kind=kind_phys),intent(in):: z
  real(kind=kind_phys),intent(out):: fxp

  if(z>0) then
     fxp=.065+11.*z-6.6e-5/z*(1.-exp(-z/8.0e-4))
  else
     fxp=0.
  endif

! end subroutine sw_ohlmann_v1
  end subroutine sw_rad_skin


   subroutine surflx(	&
   nonsol,		& ! sum of nonsolar heat fluxes, pos.up
   tsfc,		& ! skin temperature
   tlyr1,		& ! temperature in lowest atmo layer
   qlyr1,		& ! sfc.humidity in lowest atmo layer
   psfc,		& ! surface pressure
   hflx,		& ! sensible heat flux, pos.up			(out)
   qsat,		& ! satur.specf. humidity			(out)
   evap,		& ! latent heat flux, pos.up			(out)
   elocp,		& ! heat of evaporation over specif.heat, hvap/cp
   eps,			& ! ratio of air/vapor gas constants
   rch,			& ! rho * cp * ch * wind  [W/deg]
   sbc,			& ! stefan-boltzmann constant
   sfcemis,		& ! sea surface emissivity 
   dlwflx,		& ! absorbed downwelling LW flux, pos.down
   ulwflx, 		& ! surface-emitted LW flux, pos.up		(out)
   alon,alat,doprint)

! --- compute sum of nonsolar air-sea fluxes
! --- watch out for nonstandard sign convention:
! --- dlwflx is pos.down, all other fluxes pos.up.

   use funcphys, only : fpvs		! vapor pressure
   implicit none
   real, intent(in) :: tsfc, tlyr1, qlyr1, psfc, elocp, eps, rch, sbc,	&
                       sfcemis, dlwflx, alon, alat
   logical, intent(in) :: doprint
   real, intent(out) :: nonsol, qsat, evap, hflx, ulwflx
   real :: pvap, frz=273.15

   if (doprint) print 99,'entering surflx   lon,lat=',alon,alat,          &
    'nonsol',nonsol,	&
    'tsfc',tsfc-frz,	&
    'tlyr1',tlyr1-frz

   pvap = fpvs(tsfc)             ! saturation vapor pressure (pa)
   qsat = eps*pvap / (psfc + (eps-1.)*pvap)
   evap = elocp * rch * (qsat - qlyr1)
   hflx = rch * (tsfc - tlyr1)
   ulwflx = sfcemis * sbc * tsfc**4
   nonsol = hflx + evap + ulwflx - dlwflx

   if (doprint) print 99,'exiting surflx   lon,lat=',alon,alat,          &
     'tsfc',tsfc-frz,	& ! skin temperature
     'psfc',psfc*.01,	& ! surface pressure (mb)
     'pvap',pvap,	& ! saturation vapor pressure
     'qsat',qsat*1.e3,	& !saturation specif. humidity (g/kg)
     'evap',evap,	& ! latent heat flux
     'hflx',hflx,	& ! sensible heat flux
     'ulwflx',ulwflx,	& ! upwelling long-wave flux
     'dlwflx',dlwflx,	& ! downwelling long-wave flux
     'nonsol',nonsol      ! sum of nonsolar heat fluxes (pos.up)
 99 format (/a,2f7.2/(5(a8,"=",f7.2)))
 98 format (/a,2f7.2/(4(a8,"=",es11.4)))

   return
   end subroutine surflx


   subroutine enloan(delt,surflx,thkice,temice,ticold,temwat,alon,alat,doprint)

! --- single-column version of 'energy loan' ice model.
! --- ice amount represents energy 'loaned' to water column to prevent
! --- wintertime cooling below freezing level. 'loan' is paid back in summer.

   implicit none

   real,parameter :: frz=273.15
   logical,intent(IN) :: doprint
   real,   intent(IN) :: delt,alon,alat		! time step
   real,intent(INOUT) ::		&
   surflx,	& ! net total heat flux between atm and ice (W/m^2)
   thkice,	& ! grid-box averaged ice thickness (m)
   temwat,	& ! mixed layer temperaure
   temice,	& ! ice surface temperature
   ticold  	  ! previous ice surface temperature

   real ::		&
   tmelt=frz-.2,	& ! melting point (deg)
   thin=.01,		& ! min.ice thickness
   rhoice=917.,		& ! ice density (kg/m^3)
   rhowat=1000.,	& ! water density (kg/m^3)
   kice=2.04,		& ! heat conductivity in ice (W/m/deg)
   fusion=334.e3,	& ! latent heat of fusion (J/kg)
   rate=.2/3600.,	& ! max. ice melting rate (m/sec)
   fluctn=3./3600,	& ! limit on temice fluctuation (deg/sec)
   spcifh=4190.,	& ! specific heat of water (J/kg/deg)
   dpth=40.		  ! nominal mixed layer depth (m)
   real :: tnew,borrow,paybak,avail

! --- energy loan: add extra energy to the ocean to keep SST from dropping
! --- below tmelt in winter. return this borrowed energy to the 'energy bank'
! --- in summer as quickly as surflx > 0 allows.

   if (doprint) print 97,'entering enloan     lon,lat=',alon,alat,	&
    'surflx',surflx,			&
    'temwat',temwat-frz,		&
    'temice',temice-frz,		&
    'ticold',ticold-frz,		&
    'thkice',thkice
99  format (/a,2f7.2/(5(a8,"=",f7.2)))
98  format (/a,2f7.2/(4(a8,"=",es11.4)))
97  format (/a,2f7.2/(4(a8,"=",f11.6)))

   borrow=0.
   paybak=0.
   tnew=temwat+surflx*delt/(rhowat*spcifh*dpth)				! deg

   if (surflx.lt.0.) then		! cooling
    if (tnew.gt.tmelt) then 		! no action
     temwat=tnew

     if (doprint) print 97,'enloan action 1     lon,lat=',alon,alat,	&
      'tnew',tnew-frz,'temwat',temwat-frz
    else				! tnew < tmelt

! --- borrow energy to keep temwat from dropping below tmelt
     borrow=(tmelt-tnew)*rhowat*spcifh*dpth/delt			! W/m^2
     temwat=tmelt
     thkice=thkice+borrow*delt/(rhoice*fusion)				! m

     if (doprint) print 97,'enloan action 2     lon,lat=',alon,alat,	&
      'tnew',tnew-frz,'borrow',borrow,'surflx',surflx,'thkice',thkice,	&
      'temwat',temwat-frz
    end if

   else 				! warming
    if (thkice.gt.0.) then

! --- return the borrowed amount whenever tnew > tmelt

     avail=(tnew-tmelt)*rhowat*spcifh*dpth/delt				!  W/m^2
     paybak=min(thkice*rhoice*fusion/delt,avail,			&
                rate*rhoice*fusion)                          		!  W/m^2
     thkice=thkice-paybak*delt/(rhoice*fusion)				! m

     temwat=tnew-paybak*delt/(rhowat*spcifh*dpth)			! deg

     if (doprint) print 97,'enloan action 3     lon,lat=',alon,alat,	&
      'tnew',tnew-frz,'temwat',temwat-frz,'temice',temice-frz,		&
      'paybak',paybak,'surflx',surflx,'thkice',thkice

    else				! thkice = 0
     temwat=tnew
    end if
   end if				! surflx

! --- compute ice surface temperature
   if (thkice.gt.thin) then

! --- assume zero flux divergence at ice surface, so
! --- surflx = (temice-temwat)*kice/thkice

    temice=min(tmelt,temwat+thkice*surflx/kice)
! --- put limits on temice tendency
    temice=max(ticold-fluctn*delt,min(ticold+fluctn*delt,temice))

    if (doprint) print 97,'enloan action 4     lon,lat=',alon,alat,	&
     'tnew',tnew-frz,'flxice',surflx,'thkice',thkice,'temice',temice-frz

   else
    temice=temwat
   end if
   ticold=temice

   if (doprint) print 97,'exiting enloan     lon,lat=',alon,alat,	&
    'temwat',temwat-frz,	&
    'thkice',thkice,		&
    'temice',temice-frz,	&
    'ticold',ticold-frz,	&
    'tnew',tnew-frz

   return
   end subroutine enloan
end module skinsst


   subroutine get_testpt(testlon,testlat)

! --- define test point for detailed diagnostic

   implicit none
   real,intent(OUT) :: testlon,testlat

!  testlon=237.58  ; testlat= 65.58
!  testlon=273.03  ; testlat= 47.03
!  testlon= 74.06  ; testlat= 46.02
!  testlon= 30.71  ; testlat= 61.03	! high-lat lake
!  testlon= 35.00  ; testlat= 61.71	! high-lat lake
!  testlon= 32.20  ; testlat= 60.47	! high-lat lake
!  testlon=245.98  ; testlat= 67.59	! high-lat lake
!  testlon=247.21  ; testlat= 61.91	! high-lat lake
!  testlon=239.76  ! testlat= 65.90	! high-lat lake
!  testlon=237.58  ! testlat= 65.58	! high-lat lake

!  testlon=274.14  ; testlat=  46.89	! lake superior
!  testlon=278.88  ; testlat=  42.24	! lake erie
!  testlon= 48.48  ; testlat=  44.02	! caspian sea
!  testlon= 53.44  ; testlat=  37.48	! caspian sea
!  testlon=272.88  ; testlat=  43.09	! lake
!  testlon= 29.41  ; testlat=  60.37
!  testlon=302.23  ; testlat= -34.65
!  testlon=228.23  ; testlat=  55.55 
!  testlon=145.56  ; testlat=  -4.02
!  testlon=201.82  ; testlat=  21.54
!  testlon=169.22  ; testlat= -45.19
!  testlon=354.39  ; testlat=  83.24
!  testlon=155.29  ; testlat= -68.52

!  testlon= 180.84  ; testlat= -0.51
!  testlon= 180.84  ; testlat=  0.51
!  testlon= 179.82  ; testlat=  0.51
!  testlon= 179.82  ; testlat= -0.51
!  testlon= 215.00  ; testlat=-76.24

!  testlon= 272.60  ; testlat= 48.70	! lake superior
!  testlon= 273.39  ; testlat= 47.79	! lake superior
!  testlon= 271.92  ; testlat= 47.15	! lake superior
!  testlon= 292.34  ; testlat= 66.86
!  testlon= 245.31  ; testlat= 61.72	! great slave lake
!  testlon=  60.28  ; testlat= 67.99

!  testlon= 107.04  ; testlat= 53.05 	! lake baikal
!  testlon=  73.20  ; testlat=  69.0	! ob river
!  testlon=  73.53  ; testlat= 68.07	! ob river

!  print '(a,2f8.2)','(get_testpt) set test point location',testlon,testlat

   return
   end subroutine get_testpt


   real function xtinct(deglon,deglat)

! --- extinction coefficient (x 10) from SeaWIFS at 5x5 deg resolution

   implicit none
   real,intent(IN) :: deglon,deglat
   integer i,j,ib,jb
   real x,y,alon,alat
   real :: excoef(72,36)

 data excoef(37:72,1:36)		&
 / 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
  16,16,16,16,18,20,22,22,22,22, 9, 9, 9,19,20,20,20, 9, 9, 9, 9, 9,22,22,22,19,19,19,19,19,19,19, 9, 9, 9, 9, &
  15,16,15,15,16,18,20,22,22,22,18,18,18,19,19,20,20,19,17,17,19,19,22,22,18,19,18,19,18,18,19,19,19,19,20,21, &
  19,19,18,18,18,18,19,19,18,17,18,18,17,17,16,18,19,17,15,16,17,19,19,21,19,20,24,23,20,19,17,17,18,19,19,20, &
  18,19,18,18,18,18,18,19,19,19,19,19,20,20,20,20,20,20,19,19,18,17,17,20,20,21,20,19,19,19,18,18,18,19,18,18, &
  19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,19,19,18,18,18,18,18,18,19,19,19,19,19,19, &
  20,20,20,19,19,19,19,19,20,19,20,20,20,20,20,20,20,20,21,21,20,20,18,20,19,19,17,17,17,18,17,18,18,18,19,18, &
  20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,18,11,10,13,16,19,19,17,16,17,17,17,17,18,19,19, &
  19,18,19,19,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,21,14, 9, 9,11,16,16,16,17,17,18,18,18,18,19,18,19, &
  17,18,19,20,21,21,21,21,22,22,22,22,21,22,22,22,22,22,21,21,17, 9, 9, 9,13,15,15,16,17,17,17,18,18,18,19,19, &
  19,20,21,22,22,22,22,23,23,23,23,23,22,22,22,22,22,22,22,21,17,10,14, 9, 9,14,17,18,18,19,19,19,19,19,19,19, &
  21,22,22,23,23,24,24,24,24,24,25,25,25,24,24,24,23,23,22,21,18,14,16, 9, 9, 9,19,20,21,21,21,22,22,22,22,22, &
  22,23,24,24,25,25,25,25,25,26,26,26,26,26,26,25,25,24,23,22,19,16,17, 9, 9, 9,17,22,22,22,23,23,23,23,23,23, &
  23,24,25,26,26,26,26,26,26,26,26,26,26,26,26,25,24,24,23,21,20,17,17, 9, 9,17,12,18,23,23,24,24,24,24,23,22, &
  24,25,26,27,26,26,25,25,25,25,25,24,24,24,23,23,22,22,21,20,18,15, 9, 9, 9,12,18, 9,22,24,25,24,24,23,22,21, &
  25,26,26,26,25,25,24,23,23,23,23,23,23,22,22,22,21,21,20,18,12,18, 9, 9, 9, 9, 9,22,22,24,24,24,24,23,22,21, &
  24,24,24,23,23,23,22,22,22,22,22,22,22,21,21,21,21,21,19,15, 9,12, 9, 9, 9, 9,18,22,16,23,23,22,21,21,20,19, &
  22,22,21,21,21,21,21,21,21,21,21,21,20,20,20,20,20,19,19,15, 9,15, 9, 9, 9, 9, 9,18,21,21,20,19,18,17,16,16, &
  22,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,20,20,18,15,15, 9,12, 9, 9,15,21,21,21,20,19,19,18,17,17, &
  24,24,23,23,23,23,22,22,22,22,22,22,22,22,22,21,20,19,19,17,14, 9, 9, 9,12,12,18,20,20,21,20,19,17,13, 9, 9, &
  25,25,25,25,24,24,23,23,22,22,22,22,22,22,21,20,18,17,14,19,20,17,17,13,17,19,21,22,21,21,20,16, 9, 9,13,10, &
  25,25,25,24,24,24,24,23,23,23,23,23,22,22,20,19,16,14,19,20,21,20,20,20,21,22,23,23,23,22,20,15, 9, 9, 9, 9, &
  24,24,24,24,23,24,24,24,23,23,23,22,21,19,16, 9,18,19,20,19,16,22,23,23,23,23,23,23,23,22,21,18, 9, 9,15, 9, &
  23,23,23,23,23,23,23,23,23,22,22,21,19,12,11, 9,14,16,17,14,20,22,22,23,23,23,23,23,23,22,22,21,19,15, 9,16, &
  20,21,21,21,21,21,22,22,22,21,20,17,14,10,12,14,16, 9, 9,10,18,21,21,21,21,21,21,21,21,21,21,21,20,19,16,17, &
  19,19,19,19,19,19,19,20,20,19,17,11, 9, 9,10, 9, 9,11,11,18, 9,16,18,19,19,19,19,19,19,19,19,19,19,18,15,15, &
  17,17,17,17,17,17,17,18,18,18,15, 9, 9, 9, 9, 9, 9,11,11, 9, 9, 9,11,13,15,15,15,17,18,18,18,18,17,17,14,14, &
  16,16,17,17,16,17,17,17,17,17,14, 9, 9, 9, 9, 9, 9, 9,10,10, 9, 9, 9, 9,11,12,14,15,16,16,16,16,17,16,15,11, &
  14,14,14,14,15,16,17,17,16,11, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,19, 9, 9,10,14,14,15,15,15,16,16,16,15,10, 9, &
  14,12,11, 9, 9,11,14,14,12, 9,11, 9, 9, 9, 9, 9, 9,10,12,11, 9, 9,11,13,14,15,15,16,16,16,16,16,16,16,13,11, &
  11,12, 9, 9, 9, 9,10, 9, 9,12, 9, 9, 9, 9, 9, 9, 9,13,14,13,13,13,14,15,14,13,13,16,16,16,16,14,15,15,15,15, &
  11,12,11, 9, 9,15, 9, 9, 9, 9, 9, 9,11,10, 9,11,10, 9,13,13,11,10,15,15,15,13,19,19,15,15,15,14,14,15,16,16, &
  14,15,15,14,14,14,15,15,16,16,16,15,15,15,16,16,15,15,14,14,14,15,15,15,16,13,19,19,19,15,11,13,15,15,15,16, &
  16,17,17,17,17,17,17,18,18,19,21,19,17,16,18,17,15,14,12,13,15,15,15,17,17,19,19,18,18,17,16, 9,16,16,16,16, &
  24,17,17,17,16,17,22,23,22,22,21,21,19,23,23,23,19,16,13,14,11,15,13,16,18,19,18,18,17,16,15,16,18,18,17,17, &
  24,24,17,17,17,22,23,23,23,22,21,21, 9,23,23,23,23,19,16,14,15,15,16,18,19,19,19,18,18,17,16,18,18,18,18,18  /
 data excoef( 1:36,1:36)		&
 / 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
   9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,17,17,17,16,15, &
  22,22,22,21,21,21,21, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,21,17,16,15,15, &
  21,22,21,21,20,21,21,20,20,20,20,20,20,19,19,20,20,20,20,21,21,21,20,20,21,21,21,20,20,21,21,21,21,19,18,18, &
  18,18,18,19,20,19,19,20,20,20,20,20,19,19,17,18,20,20,18,20,21,19,20,19,20,21,20,20,19,19,21,20,19,18,19,19, &
  18,18,18,19,19,20,20,20,20,20,20,20,20,20,20,19,18,18,19,19,19,19,19,20,19,20,20,19,19,19,18,18,18,19,19,19, &
  19,19,18,19,19,19,19,19,20,20,20,20,20,20,20,19,19,19,19,19,19,20,20,20,20,20,20,20,20,19,19,20,20,20,20,20, &
  18,19,19,18,19,19,19,19,19,19,19,20,20,20,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20, &
  19,19,19,19,19,19,19,19,19,19,19,19,19,19,18,18,19,19,20,20,20,20,20,20,20,20,20,20,20,19,19,19,18,17,19,19, &
  19,19,19,18,18,18,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,19,18,15,17,18,18,17,13,15, &
  19,19,18,17,17,18,18,19,19,19,19,19,19,20,20,20,20,20,20,21,20,20,20,19,19,19,19,18,16,15,18,18,19,19,17,17, &
  21,20,19,13,10,16,19,20,20,21,21,21,22,22,22,22,22,22,22,21,21,21,20,15,17,18,17,11, 9,19,18,20,20,20,20,20, &
  22,20,15, 9,16,19,18,20,20,20,22,22,22,23,23,24,24,24,23,23,22,21,19,20,18,18,18,17,13,20,19,21,22,21,22,22, &
  21,18,10,15, 9,18,10,20,20,21,23,23,23,23,24,24,24,24,23,23,22,22,20,12,20,16,13,13,20,13,20,23,22,22,23,23, &
  19,14, 9,12, 9, 9, 9,18,19,17,23,23,23,23,23,23,23,23,23,23,22,22,21,20,16, 9,15,13,12,20,23,23,23,23,24,23, &
  20,17,12,12, 9, 9,12, 9,21,20,22,22,22,23,22,22,22,22,22,22,22,21,20,20,19,15,13,15,16,21,22,23,23,24,25,25, &
  16,12, 9,12, 9, 9, 9,17,21,21,21,21,22,22,22,22,22,22,22,22,20,17,16,17,18,19,17,12,11,18,21,22,23,24,24,24, &
  14,10, 9, 9, 9, 9, 9,17,19,19,20,21,22,22,22,22,22,22,22,21,15,13,10,15,19,19,17,17,20,21,22,23,23,23,22,22, &
  17,13,13, 9, 9, 9, 9, 9,18,17,18,20,21,21,20,20,20,21,21,17,11,18,14,16,21,21,21,22,23,23,23,23,23,22,22,22, &
  10, 9,13, 9, 9, 9, 9, 9,18,14,16,18,20,20,19,15,18,19,19,16,16,18,20,18,19,23,24,24,24,24,25,25,25,25,24,24, &
   9,10, 9, 9, 9, 9, 9, 9,10,13,13,16,17,18,17,10,18,19,19,16,14,16,20,21,18,25,25,26,26,26,26,26,26,26,26,25, &
   9, 9, 9, 9, 9, 9,18,13,10,13,12,12,14,15,13, 9,14,17,16, 9, 9,14,20,20,22,24,25,25,25,26,26,26,26,26,25,25, &
   9, 9, 9, 9, 9,17, 9,18,11,12, 9, 9, 9,10, 9, 9, 9, 9, 9,16, 9, 9,14,16,21,23,24,24,24,25,25,25,25,25,25,25, &
   9,16,20,20,20,20,17,18,18, 9, 9, 9, 9, 9,10,12,12,12,14,14, 9,14, 9, 9,13,19,21,22,22,23,23,23,23,23,23,23, &
  17,18,16,20,20,20,19,17, 9, 9,10,10, 9,16,16, 9,12,12,12, 9,12,11, 9, 9, 9,13,15,18,19,20,20,19,20,20,20,20, &
  17,17,18,18,17,18,18,17,12, 9, 9,10, 9,16,16, 9,11, 9, 9,12,11,11, 9, 9, 9, 9,13,13,15,16,17,17,18,18,18,18, &
  16,17,16,15,12, 9,11,11,10, 9, 9, 9, 9, 9, 9,14,11, 9,12,12,12,11, 9, 9, 9,10,12,13,12,13,15,16,16,16,17,17, &
   9, 9,10, 9,15, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,12, 9,11,11,13,14,15,16,16,16,16, &
   9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,11, 9, 9,14,14,11,14,15,15,15, &
  13, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,12,12,12,10,11,13,14,14, &
  14,12, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,12,12,12, 9,10, 9, 9,12,11, &
  16,15,13,13, 9, 9,11,10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
  16,16,16,15,15,14,14,14,14,14,15,14,12, 9, 9, 9, 9, 9,14, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,12,12,13, &
  16,16,16,15,15,16,16,16,16,16,16,16,16,15,13,13,13,13,14,13,14,14,14,14,13,13,12,12,11,12,13,13,14,15,15,15, &
  18,18,17,17,17,17,18,18,18,18,18,18,18,18,16,17,18,18,17,18,18,16,16,15,15,14,14,11,10, 9,14,17,17,18,17,16, &
  18,18,18,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,16,16,15,15,14,14,11,14,17,17,18,18,18,24  /

   alon=(deglon+2.5)/5.
   alat=(max(-87.5,min(87.5,deglat))+92.5)/5.
   i=alon
   j=alat
   ib=mod(i,72)+1
   jb=min(j+1,36)
   x=alon-i
   y=alat-j
   if (x.lt.0. .or. x.gt.1. .or. y.lt.0. .or. y.gt.1.) then
    print '(a,4f8.2,4i5)','(xtinct) error: x or y out of range',	&
     alon,alat,x,y,i,j,ib,jb
    stop
   end if
   xtinct=((excoef(i ,j )*(1.-x) + excoef(ib,j )*x)*(1.-y)		&
          +(excoef(i ,jb)*(1.-x) + excoef(ib,jb)*x)*y)

   if (xtinct.lt.9.) 							&
   print '(a,5f8.2/8x,2i4,4f8.2)','xtinct error:',deglon,deglat,xtinct,	&
     excoef(ib,j),excoef(ib,jb),i,j,x,y,excoef(i,j),excoef(i,jb)

   xtinct = xtinct * 0.1
   return
   end function xtinct

