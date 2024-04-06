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

   subroutine skinsst_run(				&
    im,			& ! horiz. loop extent				in
    iter,		& ! ccpp loop counter				in
    wet,		& ! .true. at ocean points			in
    timestep,		& ! model timestep				in
    xlon, xlat,		& ! longitude, latitude				in
    sfcemis,		& ! sea surface emissivity			in
    dlwflx,		& ! absorbed downwelling LW flux		in
    sfcnsw,		& ! net SW flux					in
    tsfco,		& ! ocean top layer temperature			in
    psfc,		& ! surface pressure				in
    stress,		& ! wind stress (N/m^2)				in
    wind,	 	& ! atm. mid-layer 1 wind			in
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
    ulwflx,             & ! upwelling LW flux				inout
    tskin,		& ! skin temp					inout
    qsat,		& ! saturation specif. humidity			out
    z_c,		& ! sub-layer cooling thickness			out
    dt_warm,		& ! warm-layer surface warming amount		out
    dt_cool,		& ! sub-layer cooling amount			inout
    evap,		& ! kinematic latent heat flux, pos.up		out
    hflx,		& ! kinematic sensible heat flux, pos.up	out
    ep,			& ! potential latent heat flux, pos.up		out
    gflux,		& ! heat flux over ground			out
    cmm,		& ! momentum exchange coeff			out
    chh,		& ! thermal exchange coeff			out
    lseaspray,		& ! sea spray flag				in
    fm,			& ! Monin-Obukhov function at surface		in
    fm10,		& ! Monin-Obukhov function at 10m		in
    errmsg, errflg)

   use machine , only : kind_phys
!  use module_nst_parameters, only : rad2deg
   use state_eqn
   use funcphys, only : fpvs		! vapor pressure

   implicit none

! --- input:
   integer, intent(in) :: im,iter
   logical, dimension(:), intent(in) :: wet
   logical, intent(in) :: lseaspray
   real (kind=kind_phys), dimension(:), intent(in) :: xlon,xlat,	&
      sfcemis, dlwflx, sfcnsw, tsfco, wind, psfc, plyr1, tlyr1, qlyr1,	&
      ulyr1, vlyr1, cm, ch, compres, stress, fm, fm10
   real (kind=kind_phys), intent(in) :: timestep, hvap, cp, rd, eps,	&
      sbc

! --- inout:
   real (kind=kind_phys), dimension(:), intent(inout) :: ulwflx,	&
      tskin, dt_cool

! --- output:
   real (kind=kind_phys), dimension(:), intent(out) :: evap, hflx,	&
      ep, qsat, gflux, cmm, chh, dt_warm, z_c
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! --- locals:
   integer :: i, n, loop
   real :: alon, alat, virt, rho_air, rho_wat, pvap, tsq, piston, vel,	&
     nonsol,			& ! sum of nonsolar air-sea fluxes
     spcifh = 3990.,		& ! seawater specific heat
     grav  = 9.806,		& ! gravity
     sss = 34.7,		& ! sea surface salinity
     penetr = 9.,		& ! shortwave penetration depth (m)
     grnblu
   integer,parameter :: itmax = 5		! regula falsi iterations
   real :: rnl_ts, hs_ts, rf_ts, alpha, beta, rch, ustar,		&
      hist(0:itmax) = 0., x1, x2, x3, y1, y2, y3, dif1, dif2, dif3

!  variables for sea spray effect
   real (kind=kind_phys) :: f10m, u10m, v10m, ws10, ru10, qss1,		&
                            bb1, hflxs, evaps, ptem, tem
   real (kind=kind_phys), parameter :: alps=0.75,bets=0.75,gams=0.15,	&
                            ws10cr=30., conlf=7.2e-9, consf=6.4e-8

   logical :: doprint, details
   real, parameter :: rad2deg = 57.2957795
   real(kind=kind_phys) :: frz=273.15, small=.05, testlon, testlat
   common /testpt/ testlon,testlat		! (values defined in dcyc2t3.f)
   doprint(alon,alat)=abs(testlon-alon).lt.small .and.			&
                      abs(testlat-alat).lt.small

! --- latitude-dependent green/blue water transition factor affecting pendep
!  grnblu(alat)= 1. - (alat/90.)**2 * (1. - .5*(alat/90.)**2)
!  grnblu(alat)= 1. - .25*(alat/90.)**2 * (3. - (alat/90.)**2 * (alat/90.)**2)
   grnblu(alat)= 1. - (alat/90.)**2 * (alat/90.)**2 * (1.5 - (alat/90.)**2)

! --- piston velocity: at 0 m/s, molecular diffusion only.
! ---                  at 8 m/s, homogenize top 3 m during 600 sec time step
   piston(vel)=1.e-4 + 1.6e-3*(vel/8.)**2

   if (iter.gt.1) return

   do i = 1,im
    if (wet(i)) then

     alon=xlon(i)*rad2deg
     alat=xlat(i)*rad2deg
     if (doprint(alon,alat)) then
      print 98,'entering skinsst_run   lon,lat=',alon,alat,		&
      'stress',stress(i),		& ! wind stress (N/m^2)
!     'sfcemis',sfcemis(i),		& ! sfc emissivity
      'wind',wind(i),			& ! surface wind
      'sfcnsw',sfcnsw(i),		& ! total sky net SW flx into ocean
      'dlwflx',dlwflx(i),		& ! absorbed downwelling LW flux
      'ulwflx',ulwflx(i),		& ! upwelling LW flux
      'psfc',psfc(i),			& ! surface pressure (mb)
      'plyr1',plyr1(i),			& ! atm.layer 1 presure
      'tlyr1',tlyr1(i)-frz,		& ! atm.layer 1 air temp
      'qlyr1',qlyr1(i),			& ! atm.layer 1 humidity
!     'sigma_t',sig(tlyr1(i)-frz,sss),	& ! sea water density - 1000
      'compres',compres(i),		& ! midlyr-to-sfc adiab.compression
      'tskin',tskin(i)-frz,		& ! surface skin temperature
      'tsfco',tsfco(i)-frz		  ! ocean top layer temperature
      print '(5(a13,"=",l2))','lseaspray',lseaspray
     end if
 99  format (/a,2f7.2/(5(a8,"=",f7.2)))
 98  format (/a,2f7.2/(4(a8,"=",es11.4)))

     details = .false.
     if (doprint(alon,alat)) details = .true.

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

! --- apply warm layer correction
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- dt_warm currently archived under the name "xtts".
! --- other unused arrays: c0,w0,wd,xs,xt,xu,xv
! - - - - - - - - - - - - - - - - - - - - - - - - - - -

     dt_warm(i) = sfcnsw(i) * timestep 					&
         / (2. * penetr * grnblu(alat) * rho_wat * spcifh)
! --- note:  dt_warm is cumulative, dt_cool is not. remember: dt_cool > 0
     tskin(i) = tskin(i) + dt_warm(i)					&
! --- subtract old dt_cool
        + dt_cool(i)
! --- allow cooling by vertical diffusion
     tskin(i) = tskin(i)						&
        + (tsfco(i)-tskin(i)) * min(1.,piston(wind(i)) * timestep)

! --- start cool-skin iteration, using regula falsi (aiming for x_n = y_n)
! --- x1,x2,x3,y1,y2,y3 are consecutive dt_cool approximations.

     x1 = -.5
     x2 = +.5
     call surflx(nonsol, tskin(i)+x1, tlyr1(i)*compres(i), qlyr1(i),	&
         psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
         sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, doprint(alon,alat))

     call coolskin(ustar,nonsol,sfcnsw(i),evap(i),sss,alpha,beta,	&
                   rho_wat,rho_air,tskin(i)+x1,grav,hvap,		&
                   y1,z_c(i),alon,alat,doprint(alon,alat))

     call surflx(nonsol, tskin(i)+x2, tlyr1(i)*compres(i), qlyr1(i),	&
         psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
         sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, doprint(alon,alat))

     call coolskin(ustar,nonsol,sfcnsw(i),evap(i),sss,alpha,beta,	&
                   rho_wat,rho_air,tskin(i)+x2,grav,hvap,		&
                   y2,z_c(i),alon,alat,doprint(alon,alat))

     if (details) print '(2f7.2,a,3es11.3)',alon,alat,'  x1,y1 =',x1,y1,y1-x1
     if (details) print '(2f7.2,a,3es11.3)',alon,alat,'  x2,y2 =',x2,y2,y2-x2

     do loop = 1,itmax

      dif1 = y1 - x1
      dif2 = y2 - x2
      if (abs(dif2-dif1).gt..001*(abs(dif1)+abs(dif2))) then
       if (dif1*dif2.gt.0.) then
         x3 = x1-dif1*(x2-x1)/(dif2-dif1)
       else
         x3 = (dif1*x2-dif2*x1)/(dif1-dif2)
       end if

       call surflx(nonsol, tskin(i)+x3, tlyr1(i)*compres(i), qlyr1(i),	&
         psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
         sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, doprint(alon,alat))

       call coolskin(ustar,nonsol,sfcnsw(i),evap(i),sss,alpha,beta,	&
                     rho_wat,rho_air,tskin(i)+x3,grav,hvap,		&
                     y3,z_c(i),alon,alat,doprint(alon,alat))

      else			! we have convergence
        x3 = x2
        y3 = y2
      end if
      if (details)			&
      print '(2f7.2,a,3es11.3,i4)',alon,alat,'  x3,y3 =',x3,y3,y3-x3,iter
      dt_cool(i) = y3
      hist(loop) = y3
      dif3 = y3 - x3
      if (abs(dif3).lt..001*(abs(dif1)+abs(dif2))) exit

      if (abs(dif1).lt.abs(dif2)) then
        x2 = x3
        y2 = y3
      else
        x1 = x3
        y1 = y3
      end if

      if (loop.eq.itmax) then 
       if (abs(hist(loop)).gt..5) then
        print '(a,3f8.2/(11f7.2))','tskin not converging at lon,lat',	&
        alon,alat,hist(loop),(hist(n),n=1,loop)
       end if
      end if

     end do		! iteration loop

     tskin(i) = tskin(i)-dt_cool(i)		! dt_cool > 0 => cooling

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

     if (doprint(alon,alat))						&
      print 98,'exiting skinsst_run   lon,lat=',alon,alat,		&
      'loop',float(loop),		&
      'virt',virt-frz,			& ! virtual temp
      'rho_air',rho_air,		& ! air density
      'pvap',pvap,			& ! satur. vapor pressure (mb)
      'qsat',qsat(i),		 	& ! satur. specif.humidity
      'hflx',hflx(i),			& ! sensible heat flux
      'evap',evap(i),			& ! latent heat flux
      'nonsol',nonsol,			& ! net non-solar surface flux
      'sfcnsw',sfcnsw(i),		& ! net solar surface flux
      'ulwflx',ulwflx(i),		& ! upwelling LW flux
      'dt_cool',dt_cool(i),		& ! cool-surface correction (set > 0)
      'dt_warm',dt_warm(i),		& ! temperature increment due to SW
      'tskin',tskin(i)-frz,		& ! skin temperature
      'cmm',cmm(i),			&
      'chh',chh(i),			&
      'spray',bets*hflxs-ptem*evaps	  ! spray contrib. to heat flux

! --- convert fluxes from W/m^2 to "kinematic" i.e. velocity x fluxed variable
     hflx(i) = hflx(i)/(rho_air * cp)				! deg m/sec
     evap(i) = evap(i)/(rho_air * hvap)				! m/sec

    end if		! wet
   end do		! im loop

   return
   end subroutine skinsst_run


  subroutine coolskin(ustar_a,f_nsol,f_sol_0,evap,sss,alpha,beta,	&
                      rho_w,rho_a,ts,grav,le,deltat_c,z_c,alon,alat,doprint)

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
! le      : latent heat of evaporation 

! output:
! deltat_c: cool-skin temperature correction (>0 means cooling)
! z_c     : molecular sublayer (cool-skin) thickness (m)

  use machine , only : kind_phys
  use module_nst_parameters, only: kw => tc_w,visw,cp_w, z_c_max,	&
    z_c_ini,ustar_a_min,rad2deg
  implicit none
  logical,intent(in) :: doprint
  real(kind=kind_phys), intent(in) :: ustar_a,f_nsol,f_sol_0,evap,	&
     sss,alpha,beta,rho_w,rho_a,ts,grav,le,alon,alat
  real(kind=kind_phys), intent(out):: deltat_c,z_c
! local variables:
  real(kind=kind_phys), parameter :: frz=273.15
  real(kind=kind_phys) :: xi,hb,ustar1_a,bigc,deltaf,fxp

  if (doprint) print 98,'entering coolskin   lon,lat=',alon,alat,	&
     'ustar_a',ustar_a, &	! atmospheric friction velocity (m/s)
     'f_nsol',f_nsol,   &	! "nonsolar" part of the surface heat flux
     'f_sol_0',f_sol_0, &	! solar radiation at the ocean surface
     'evap',evap,	&	! latent heat flux (w/m^2)
     'rho_w',rho_w,	&	! sea water density
     'rho_a',rho_a,	&	! air density
     'ts',ts-frz		! ocean surface temperature
 99 format (/a,2f7.2/(5(a8,"=",f7.2)))
 98 format (/a,2f7.2/(4(a8,"=",es11.4)))

  z_c = z_c_ini                 ! initial guess

  ustar1_a = max(ustar_a,ustar_a_min)

  call sw_rad_skin(z_c,fxp)
  deltaf = f_sol_0 * fxp

  hb   = alpha * (f_nsol-deltaf) + beta * sss * cp_w * evap/le
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
!    'fxp',fxp,                 &
!    'deltaf',deltaf,           &
!    'hb',hb,                   &
!    'bigc',bigc,               &
!    'xi',xi,                   &
     'deltat',deltat_c,         &	! skin layer temperature correction
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
   tsfc,		& ! sea surface temperature
   tlyr1,		& ! temperature in lowest atmo layer
   qlyr1,		& ! spfc.humidity in lowest atmo layer
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
! --- watch out for nonstandard sign convention: upwd fluxes are treated as > 0

   use funcphys, only : fpvs		! vapor pressure
   implicit none
   real, intent(in) :: tsfc, tlyr1, qlyr1, psfc, elocp, eps, rch, sbc,	&
                       sfcemis, dlwflx, alon, alat
   logical, intent(in) :: doprint
   real, intent(out) :: nonsol, qsat, evap, hflx, ulwflx
   real :: pvap, frz=273.15

   pvap = fpvs(tsfc)             ! saturation vapor pressure (pa)
   qsat = eps*pvap / (psfc + (eps-1.)*pvap)
   evap = elocp * rch * (qsat - qlyr1)
   hflx = rch * (tsfc - tlyr1)
   ulwflx = sfcemis * sbc * tsfc**4
   nonsol = hflx + evap + ulwflx - dlwflx

  if (doprint) print 98,'exiting surflx   lon,lat=',alon,alat,          &
     'tsfc',tsfc-frz,	& ! surface temperature
     'psfc',psfc,	& ! surface pressure
     'pvap',pvap,	& ! saturation vapor pressure
     'qsat',qsat,	& !saturation specif. humidity
     'evap',evap,	& ! latent heat flux
     'hflx',hflx,	& ! sensible heat flux
     'ulwflx',ulwflx,	& ! upwelling long-wave flux
     'dlwflx',dlwflx,	& ! downwelling long-wave flux
     'nonsol',nonsol      ! sum of nonsolar heat fluxes
 99 format (/a,2f7.2/(5(a8,"=",f7.2)))
 98 format (/a,2f7.2/(4(a8,"=",es11.4)))

   return
   end subroutine surflx

end module skinsst
