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
    wet,		& ! .true. at ocean & lake points		in
    oceanfrac,		& ! cell portion covered by ocean		in
    timestep,		& ! model timestep				in
    xlon, xlat,		& ! longitude, latitude				in
    sfcemis,		& ! sea surface emissivity			in
    dlwflx,		& ! absorbed downwelling LW flux		in
    sfcnsw,		& ! net SW flux					in
    tsfco,		& ! ocean/lake top layer temperature		in
    tref,		& ! foundation temperature			in
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
    xzts,		& ! holding place for tskin			inout
    xs,			& ! holding place for dt_cool			inout
    qsat,		& ! saturation specif. humidity			out
    z_c,		& ! sub-layer cooling thickness			out
    dt_warm,		& ! warm-layer surface warming amount		out
    dt_cool,		& ! skin layer cooling amount			inout
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
   use state_eqn
   use funcphys, only : fpvs		! vapor pressure

   implicit none

! --- input:
   integer, intent(in) :: im,iter
   logical, dimension(:), intent(in) :: wet
   logical, intent(in) :: lseaspray
   real (kind=kind_phys), dimension(:), intent(in) :: xlon,xlat,	&
      sfcemis, dlwflx, sfcnsw, tsfco, wind, psfc, plyr1, tlyr1, qlyr1,	&
      ulyr1, vlyr1, cm, ch, compres, stress, fm, fm10,oceanfrac,tref
   real (kind=kind_phys), intent(in) :: timestep, hvap, cp, rd, eps,	&
      sbc

! --- inout:
   real (kind=kind_phys), dimension(:), intent(inout) :: ulwflx,	&
      tskin, dt_cool, xzts, xs

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
     penetr = 6.,		& ! shortwave penetration depth (m)
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
   real(kind=kind_phys) :: frz=273.15, small=.25, testlon, testlat
   real,parameter :: rad2deg = 57.2957795
   doprint(alon,alat)=abs(testlon-alon).lt.small .and.			&
                      abs(testlat-alat).lt.small

! --- latitude-dependent green/blue water transition factor affecting pendep
!  grnblu(alat)= 1. - (alat/90.)**2 * (1. - .5*(alat/90.)**2)
!  grnblu(alat)= 1. - .25*(alat/90.)**2 * (3. - (alat/90.)**2 * (alat/90.)**2)
   grnblu(alat)= 1. - (alat/90.)**2 * (alat/90.)**2 * (1.5 - (alat/90.)**2)

! --- piston velocity: at 0 m/s, molecular diffusion only.
! ---                  at 8 m/s, destroy warm layer during 600 sec time step
! ---                  based on 1 m thickness scale.
!  piston(vel)= 1.4e-7 + 1.66667e-3*(vel/8.)**2		! quadratic interpolation
   piston(vel)= 1.4e-7 + 1.66667e-3*(vel/8.)		! linear interpolation

   if (iter.gt.1) return

    call get_testpt(testlon,testlat)
! --- temporary
!   print '(a,2f8.2)','entering skinsst_run, testpt =',testlon,testlat

   do i = 1,im
    if (wet(i)) then

     alon=xlon(i)*rad2deg
     alat=xlat(i)*rad2deg
     if (doprint(alon,alat)) then
      print 97,'entering skinsst_run   lon,lat=',alon,alat,		&
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
      'xzts',xzts(i)-frz,		& ! holding place for tskin
      'xsE2',xs(i)*100.,		& ! holding place for dt_cool x 100
      'tref',tref(i)-frz,		& ! foundation temperature
      'ocnfrc',oceanfrac(i),		& ! ocean fraction
      'tsfco',tsfco(i)-frz		  ! ocean top layer temperature
     end if
 99  format (/a,2f7.2/(5(a8,"=",f7.2)))
 98  format (/a,2f7.2/(4(a8,"=",es11.4)))
 97  format (/a,2f7.2/(4(a8,"=",f11.6)))

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
! --- tskin from last time step has been saved in xzts
! --- dt_cool from last time step has been saved in xs
! --- other unused arrays: c0,w0,wd,xt,xu,xv
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if (xzts(i).eq.0.) then			! initial xzts assumed to be 0
      tskin(i) = tsfco(i)
     else
      tskin(i) = xzts(i)			! tskin from previous time step
      dt_cool(i) = xs(i)
     end if
     dt_warm(i)=0.

     if (tskin(i)-frz.lt.-5.) print '(a,2f8.2,es13.3)',			&
         'excessively cold tskin at lon,lat',alon,alat,tskin(i)-frz

     if (sfcnsw(i).lt.0.) then
      tskin(i)=tsfco(i)
     else 					! SW flux > 0

! --- evaluate warmr-layer increment during current time step

      dt_warm(i) = sfcnsw(i) * timestep 				&
         * 2./(penetr * grnblu(alat) * rho_wat * spcifh)
! --- note: dt_warm is cumulative.
      tskin(i) = tskin(i) + dt_warm(i)					&
! --- subtract old dt_cool (since dt_cool is not cumulative)
        + dt_cool(i)				! sign convention: dt_cool > 0
! --- cooling by heat diffusion
      tskin(i) = tskin(i) + (tsfco(i)-tskin(i))				&
         * min(1.,timestep*piston(wind(i)))
     end if					! SW flux > 0

! --- save (tskin - top layer T) for diagnostic purposes
     dt_warm(i) = tskin(i) - tsfco(i)

! --- start cool-skin iteration, using regula falsi (aiming for x_n = y_n)
! --- x1,x2,x3,y1,y2,y3 are consecutive dt_cool approximations.

     details = doprint(alon,alat)

     x1 = -.5
     x2 = +.5
     call surflx(nonsol, tskin(i)+x1, tlyr1(i)*compres(i), qlyr1(i),	&
         psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
         sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, doprint(alon,alat))

     call coolskin(ustar,nonsol,sfcnsw(i),evap(i),sss,alpha,beta,	&
                   rho_wat,rho_air,tskin(i)+x1,grav,hvap,		&
                   y1,z_c(i),alon,alat,doprint(alon,alat))
     dif1 = y1 - x1

     if (y1.ne.0.) then
     
      do loop = 1,itmax

       call surflx(nonsol, tskin(i)+x2, tlyr1(i)*compres(i), qlyr1(i),	&
           psfc(i), hflx(i), qsat(i), evap(i), hvap/cp, eps, rch, sbc,	&
           sfcemis(i), dlwflx(i), ulwflx(i), alon, alat, doprint(alon,alat))

       call coolskin(ustar,nonsol,sfcnsw(i),evap(i),sss,alpha,beta,	&
                     rho_wat,rho_air,tskin(i)+x2,grav,hvap,		&
                     y2,z_c(i),alon,alat,doprint(alon,alat))
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

      end do		! regula falsi iteration loop

      if (loop.eq.itmax) then 
       if (abs(hist(loop)).gt..5) then
        print '(a,3f8.2/(11f7.2))','tskin not converging at lon,lat',	&
           alon,alat,hist(loop),(hist(n),n=1,loop)
          end if
      end if

! --- apply cool-skin correction
      tskin(i) = tskin(i)-dt_cool(i)
     end if				! y1 nonzero

! --- if no lake model, do this:
     if (oceanfrac(i).eq.0.) then
      tskin(i) = .333333*(tlyr1(i)+2.*tref(i))
      dt_cool(i) = 0.
     end if
 
! --- save tskin and dt_cool for next call to skinsst

     xzts(i) = tskin(i)			! for use at next time step
     xs(i) = dt_cool(i)			! for use at next time step

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
      'cmm',cmm(i),			&
      'chh',chh(i),			&
      'spray',bets*hflxs-ptem*evaps	  ! spray contrib. to heat flux
      if (oceanfrac(i).eq.0.) print '(2f7.2,a)',alon,alat,' is lake point'
     end if

! --- convert fluxes from W/m^2 to "kinematic" i.e. velocity x fluxed variable
     hflx(i) = hflx(i)/(rho_air * cp)				! deg m/sec
     evap(i) = evap(i)/(rho_air * hvap)				! m/sec

    end if		! wet
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

  if (doprint) print 99,'entering coolskin   lon,lat=',alon,alat,	&
     'ustar_a',ustar_a, &	! atmospheric friction velocity (m/s)
     'f_nsol',f_nsol,   &	! "nonsolar" part of the surface heat flux
     'f_sol_0',f_sol_0, &	! solar radiation at the ocean surface
     'evap',evap,	&	! latent heat flux (w/m^2)
     'rho_w',rho_w,	&	! sea water density
     'rho_a',rho_a,	&	! air density
     'ts',ts-frz,	&	! ocean surface temperature
     'zc_ini',z_c_ini*1.e3	! sublayer thickness (mm)

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

  if (doprint) print 99,'exiting coolskin   lon,lat=',alon,alat,	&
     'fxp',fxp,                 &
     'deltaf',deltaf,           &
     'hb',hb,                   &
     'bigc',bigc,               &
     'xi',xi,                   &
     'delt_c',deltat_c,         &	! skin layer temperature correction
     'z_c',z_c*1.e3			! skin layer thickness (mm)

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

  if (doprint) print 99,'exiting surflx   lon,lat=',alon,alat,          &
     'tsfc',tsfc-frz,	& ! skin temperature
     'psfc',psfc*.01,	& ! surface pressure (mb)
     'pvap',pvap,	& ! saturation vapor pressure
     'qsat',qsat*1.e3,	& !saturation specif. humidity (g/kg)
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


   subroutine get_testpt(testlon,testlat)
   real,intent(OUT) :: testlon,testlat

!  testlon=302.23  ; testlat= -34.65
!  testlon= 48.48  ; testlat=  44.02	! lake
!  testlon= 53.44  ; testlat=  37.48	! lake
!  testlon=272.88  ; testlat=  43.09	! lake
!  testlon=228.23  ; testlat=  55.55 
!  testlon=145.56  ; testlat=  -4.02
!  testlon=201.82  ; testlat=  21.54
!  testlat=169.22  ; testlon= -45.19

!  testlon= 180.84  ; testlat= -0.51
!  testlon= 180.84  ; testlat=  0.51
!  testlon= 179.82  ; testlat=  0.51
!  testlon= 179.82  ; testlat= -0.51
!  testlon= 215.00  ; testlat=-76.24

!  testlon= 272.60  ; testlat= 48.70
!  testlon= 273.39  ; testlat= 47.79	! lake
!  testlon= 271.92  ; testlat= 47.15	! lake
!  testlon= 292.34  ; testlat= 66.86
!  testlon= 245.31  ; testlat= 61.72

!  print '(a,2f8.2)','(get_testpt) set test point location',testlon,testlat
   return
   end subroutine get_testpt
