!>\file gsd_chem_driver.F90
!! This file is GSD Chemistry driver with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 01/2020

 module gsd_chem_driver

   use physcons,        only : grav => con_g
   use machine ,        only : kind_phys
   use gsd_chem_prep,   only : gsd_chem_prep_run
   use gsd_chem_config
   use dep_dry_mod
   use dep_wet_ls_mod
   use gocart_settling_mod
   use gocart_aerosols_mod
   use gocart_dmsemis_mod
   use gocart_chem_mod
   use opt_mod
   use seas_mod,        only : gocart_seasalt_driver
   use dust_gocart_mod, only : gocart_dust_driver
   use dust_afwa_mod,   only : gocart_dust_afwa_driver
   use dust_fengsha_mod,only : gocart_dust_fengsha_driver
   use dust_data_mod
   use plume_rise_mod
   use vash_settling_mod

   implicit none

   private

   public :: gsd_chem_driver_init, gsd_chem_driver_run, gsd_chem_driver_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_driver_init()
      end subroutine gsd_chem_driver_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_driver_finalize Argument Table
!!
      subroutine gsd_chem_driver_finalize()
      end subroutine gsd_chem_driver_finalize

!> \defgroup gsd_chem_group GSD Chem driver Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_driver GSD Chem driver Module  
!> \ingroup gsd_chem_group
!! This is the GSD Chem driver Module
!! \section arg_table_gsd_chem_driver_run Argument Table
!! \htmlinclude gsd_chem_driver_run.html
!!
!>\section gsd_chem_driver GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_driver_run(im, kte, kme, ktau, dt, garea, land, &
                   u10m, v10m, ustar, rlat, rlon, tskin,julian,xcosz,   &
                   rain_cpl, rainc_cpl, hf2d, pb2d,                     &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, exch, dqdt,                                       &
                   nsoil, smc, vegtype, soiltyp, sigmaf,jdate,          & 
                   dswsfc, zorl,snow_cpl,                               &
                   dust_in,emi_in,emi2_in,fire_GBBEPx,fire_MODIS,       &
                   nseasalt,ntrac,                                      &
                   ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                     &
                   ntbc1,ntbc2,ntoc1,ntoc2,                             &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,      &
                   gq0,ebu,tile_num,                                    &
                   lmk,nbdsw,faersw1,faersw2,faersw3,                   & ! radiation feedback
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im, kte, kme, ktau, nsoil, jdate(8), tile_num
    integer,        intent(in) :: nseasalt,ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind_phys),intent(in) :: dt,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype, soiltyp        
    real(kind_phys), dimension(im,nsoil), intent(in) :: smc
    real(kind_phys), dimension(im,    5), intent(in) :: dust_in
    real(kind_phys), dimension(im,   10), intent(in) :: emi_in
    real(kind_phys), dimension(im,64, 3), intent(in) :: emi2_in
    real(kind_phys), dimension(im,    5), intent(in) :: fire_GBBEPx
    real(kind_phys), dimension(im,   13), intent(in) :: fire_MODIS
    real(kind_phys), dimension(im), intent(in) :: u10m, v10m, ustar,              &
                garea, rlat,rlon, tskin, rain_cpl, rainc_cpl,                     &
                hf2d, pb2d, sigmaf, dswsfc, zorl, snow_cpl, xcosz
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, exch, dqdt
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ebu), intent(inout) :: ebu
    integer,        intent(in) :: lmk,nbdsw
    real(kind_phys), dimension(im, lmk, nbdsw),intent(inout) :: faersw1, faersw2, faersw3
!    real(kind_phys), dimension(im,nseasalt), intent(inout) :: ssem
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!   integer, parameter :: num_moist=3, num_chem=20, num_emis_seas=5, num_emis_dust=5
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, zmid,        &
                     exch_h, dqdti

    real(kind_phys), dimension(ims:im, jms:jme) :: u10, v10, ust, tsk,            &
                     xland, xlat, xlong, dxy, rcav, rnav, hfx, pbl

!>- sea salt & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  ::                  &
                     var_rmv, dry_fall, tr_fall, sedim
    real(kind_phys), dimension(ims:im, 1, jms:jme, 1:num_emis_seas  ) :: emis_seas
    real(kind_phys), dimension(ims:im, jms:jme) :: seashelp

    integer :: ide, ime, ite, kde, julday

!   integer, parameter :: SEAS_OPT_DEFAULT = 1
    integer, parameter :: chem_in_opt = 0  ! 0 for coldstart, 1 for restart
    logical, parameter :: readrestart = .false.
    integer, parameter :: nvl_gocart  = 64  ! number of input levels from gocart file
   
!>- dust & chemistry variables
    real(kind_phys), dimension(ims:im, jms:jme, 3) ::    erod ! read from input?
    real(kind_phys), dimension(ims:im, jms:jme) :: ssm, rdrag, uthr, snowh  ! fengsha dust
    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, gsw, znt, clayf, sandf, dms_0
    real(kind_phys), dimension(ims:im, nsoil, jms:jme) :: smois
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:num_emis_dust) :: emis_dust
    real(kind_phys), dimension(ims:im, 1:1, jms:jme, 1:5)             :: srce_dust
    real(kind_phys), dimension(ims:im, jms:jme) :: dusthelp
    integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp

    integer :: current_month

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: pm10, pm2_5_dry, pm2_5_dry_ec
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ac3, ahno3, anh3, asulf, cor3, h2oai, h2oaj, nu3
    real(kind_phys), dimension(ims:im, jms:jme) :: dep_vel_o3, e_co, ash_fall

!>- chemical background variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: backg_oh,backg_h2o2,backg_no3

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: oh_t, h2o2_t, no3_t
    real(kind_phys), dimension(ims:im, jms:jme) :: ttday, tcosz

!>- optical variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: relhum
    real(kind_phys), dimension(ims:im,         jms:jme) :: aod
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: extt
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: ssca
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:nbands) :: asympar
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:4) ::                  &
        tauaersw, gaersw, waersw, bscoefsw,                                        &
        l2aer,  l3aer, l4aer, l5aer, l6aer, l7aer           
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:16) :: tauaerlw
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ext_coef) :: ext_coeff
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_bscat_coef) :: bscat_coeff
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_asym_par)   :: asym_par
!>- optical variables
    real(kind_phys), dimension(im) :: aod2d
    real(kind_phys), dimension(im, kte, 1:nbands) :: ext_cof, sscal, asymp

!>- plume variables
    ! -- buffers
   !real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_ebu) :: ebu
    real(kind_phys), dimension(ims:im, jms:jme, num_ebu_in) :: ebu_in
    real(kind_phys), dimension(ims:im, jms:jme) ::                             &
         mean_fct_agef, mean_fct_aggr, mean_fct_agsv, mean_fct_agtf,            &
         firesize_agef, firesize_aggr, firesize_agsv, firesize_agtf
    real(kind_phys), dimension(ims:im, jms:jme, num_frp_plume ) :: plume_frp
    real(kind_phys), dimension(ims:im, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(kind_phys) :: dtstep, gmt
    integer,parameter :: plumerise_flag = 2  ! 1=MODIS, 2=GBBEPx
    logical :: call_plume, scale_fire_emiss
    logical, save :: firstfire = .true.
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    ! -- output tracers
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem) :: wet_dep
!*  real(kind_phys), dimension(ims:im, jms:jme, 1:nvl, 1:nbands) :: ext_cof, sscal, asymp
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: p10, pm25, ebu_oc 
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: oh_bg, h2o2_bg, no3_bg


!>-- local variables
    real(kind_phys) :: curr_secs
    real(kind_phys) :: factor, factor2, factor3
    logical :: call_gocart, call_radiation
    logical :: store_arrays
    integer :: nbegin, nv, nvv
    integer :: i, j, jp, k, kp, n
  

    errmsg = ''
    errflg = 0

!    print*,'hli test1 ktau',ktau
    h2oai = 0.
    h2oaj = 0.
    nu3   = 0.
    ac3   = 0.
    cor3  = 0.
    asulf = 0.
    ahno3 = 0.
    anh3  = 0.
    e_co  = 0.
    dep_vel_o3 = 0.
    ash_fall   = 0.
    extt =0.
    ssca   =0.
    asympar=0.


    gmt = real(jdate(5))
    julday = real(julian)                                       

    current_month=jdate(2)
    curr_secs = ktau * dt

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- initialize large-sacle wet depostion
    if (ktau==1) then
     call dep_wet_ls_init()
    endif

    ! -- set control flags
!    call_plume = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq > 0)
!    if (call_plume) &
!       call_plume    = (mod(ktau, max(1, int(60*plumerisefire_frq/dt))) == 0) &
!                        .or. (ktau == 1) .or. firstfire
!    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
!    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or. (ktau == 1)
    call_plume       = (biomass_burn_opt == BURN_OPT_ENABLE) .and. (plumerisefire_frq > 0)
    if (call_plume) &
       call_plume    = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) &
                        .or. (ktau == 1)
    call_gocart      = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)
    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq)) == 0) .or. (ktau == 1)
    scale_fire_emiss = .false.

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
     ! -- initialize buffers
!      ebu   = 0._kind_phys
    end if

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cpl(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
    enddo

!!!

!>- get ready for chemistry run
!   print*,'GSD vegtype,soiltyp,sigmaf',vegtype,soiltyp,sigmaf
    call gsd_chem_prep_run(                                             &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                      &
        u10m,v10m,ustar,land,garea,rlat,rlon,tskin,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                 &
        exch,dqdt,                                                      &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                   &
        snow_cpl,dust_in,emi_in,emi2_in,                                &
        fire_GBBEPx,fire_MODIS,hf2d,pb2d,                               &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                   &
        t8w,exch_h,dqdti,                                               &
        z_at_w,vvel,zmid,                                               &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,                                        &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,gq0,                                                      &
        num_chem, num_moist,num_ebu_in,                                 &
        call_gocart,nvl_gocart,                                         &
        ttday,tcosz,gmt,julday,                                         &
        backg_oh,backg_h2o2,backg_no3,                                  &
        plumerise_flag,num_plume_data,num_emis_ant,                     &
        emis_ant,ppm2ugkg,                                              &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
        moist,chem,plume_frp,ebu_in,                                    &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,               &
        relhum,snowh,clayf,rdrag,sandf,ssm,uthr,dms_0,erod,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
!    print*,'hli test2 ktau',ktau

!>- compute sea-salt
    ! -- compute sea salt
    if (seas_opt >= SEAS_OPT_DEFAULT) then
    call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist,                 &
        u_phy,v_phy,chem,rho_phy,dz8w,u10,v10,ust,p8w,tsk,              &
        xland,xlat,xlong,dxy,grav,emis_seas,                            &
        seashelp,num_emis_seas,num_moist,num_chem,seas_opt,             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
    endif

    !-- compute dust
    !store_arrays = .false.
    !print*,'hli dust_opt',dust_opt
    !print*,'hli DUST_OPT_AFWA',DUST_OPT_AFWA
    select case (dust_opt)
      case (DUST_OPT_AFWA)
        dust_alpha = afwa_alpha
        dust_gamma = afwa_gamma
        call gocart_dust_afwa_driver(ktau,dt,rri,t_phy,moist,u_phy,     &
          v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp, &
          vegfrac,xland,xlat,xlong,gsw,dxy,grav,emis_dust,srce_dust,    &
          dusthelp,ust,znt,clayf,sandf,                                 &
          num_emis_dust,num_moist,num_chem,nsoil,                       &
          ids,ide, jds,jde, kds,kde,                                    &
          ims,ime, jms,jme, kms,kme,                                    &
          its,ite, jts,jte, kts,kte)
       !store_arrays = .true.
      case (DUST_OPT_FENGSHA)
       dust_alpha    = 2.0  !fengsha_alpha
       dust_gamma    = 1.8  !fengsha_gamma
       call gocart_dust_fengsha_driver(dt,chem,rho_phy,smois,p8w,ssm,   &
            isltyp,vegfrac,snowh,xland,dxy,grav,emis_dust,ust,znt,      &
            clayf,sandf,rdrag,uthr,                                     &
            num_emis_dust,num_moist,num_chem,nsoil,                     &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte)
       !store_arrays = .true.
      case (DUST_OPT_GOCART)
        dust_alpha = gocart_alpha
        dust_gamma = gocart_gamma
        call gocart_dust_driver(chem_opt,ktau,dt,rri,t_phy,moist,u_phy, &
          v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp, &
          vegfrac,xland,xlat,xlong,gsw,dxy,grav,emis_dust,srce_dust,    &
          dusthelp,num_emis_dust,num_moist,num_chem,nsoil,              &
          current_month,                                                &
          ids,ide, jds,jde, kds,kde,                                    &
          ims,ime, jms,jme, kms,kme,                                    &
          its,ite, jts,jte, kts,kte)
       !store_arrays = .true.
    end select

    ! compute wild-fire plumes
    if (call_plume) then
      call plumerise_driver (ktau,dtstep,num_chem,num_ebu,num_ebu_in,   &
        ebu,ebu_in,                                                     &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
        'GOCART','BIOMASSB', t_phy,moist(:,:,:,p_qv),                   &
        rho_phy,vvel,u_phy,v_phy,p_phy,                                 &
        z_at_w,scale_fire_emiss,plume_frp,plumerise_flag,               &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte                                )
    end if

    if (dmsemis_opt == DMSE_OPT_ENABLE) then
      call gocart_dmsemis(dt,rri,t_phy,u_phy,v_phy,                     &
         chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,                       &
         ivgtyp,isltyp,xland,dxy,grav,mwdry,                            &
         num_chem,p_dms,                                                &
         ids,ide, jds,jde, kds,kde,                                     &
         ims,ime, jms,jme, kms,kme,                                     &
         its,ite, jts,jte, kts,kte)
    endif

    if ((dust_opt /= DUST_OPT_NONE) .or.                                &
        (seas_opt /= SEAS_OPT_NONE)) then
      call gocart_settling_driver(dt,t_phy,moist,                       &
        chem,rho_phy,dz8w,p8w,p_phy,sedim,                              &
        dusthelp,seashelp,dxy,grav,                                     &
        num_moist,num_chem,                                             &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)
    end if

      ! -- 4 volcanic size bins
      call vashshort_settling_driver(dt,t_phy,moist, &
           chem,rho_phy,dz8w,p8w,p_phy,dxy,          &
           ash_fall,grav,num_moist,num_chem,         &
           ids,ide, jds,jde, kds,kde,                &
           ims,ime, jms,jme, kms,kme,                &
           its,ite, jts,jte, kts,kte)


    ! -- add biomass burning emissions at every timestep
    if (biomass_burn_opt == BURN_OPT_ENABLE) then
      jp = jte
      factor3 = 0._kind_phys
      select case (plumerise_flag)
        case (FIRE_OPT_MODIS)
          factor3 = 4.828e-04_kind_phys/60.
          kp = kte    ! full column
        case (FIRE_OPT_GBBEPx)
          factor3 = 1.e-03_kind_phys * mwdry / mw_so2_aer
          if (plumerisefire_frq > 0) then
            kp = kte  ! full column
          else
            kp = kts  ! surface only
          end if
        case default
          ! -- no further options available, skip this step
          jp = jts - 1
      end select

      if (kp == kts) then
        ! -- only include surface emissions
        k = kts
        do j = jts, jp
          do i = its, ite
            ! -- factor for pm emissions, factor2 for burn emissions
            factor  = dt*rri(i,k,j)/dz8w(i,k,j)
            factor2 = factor * factor3
            chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu_in(i,j,p_ebu_in_oc  )
            chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu_in(i,j,p_ebu_in_bc  )
            chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu_in(i,j,p_ebu_in_pm25)
            chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu_in(i,j,p_ebu_in_pm10)
            chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu_in(i,j,p_ebu_in_so2 )
          end do
        end do

      else
        ! -- use full-column emissions
        do j = jts, jp
          do k = kts, kp
            do i = its, ite
              ! -- factor for pm emissions, factor2 for burn emissions
              factor  = dt*rri(i,k,j)/dz8w(i,k,j)
              factor2 = factor * factor3
              chem(i,k,j,p_oc1) = chem(i,k,j,p_oc1) + factor  * ebu(i,k,j,p_ebu_oc  )
              chem(i,k,j,p_bc1) = chem(i,k,j,p_bc1) + factor  * ebu(i,k,j,p_ebu_bc  )
              chem(i,k,j,p_p25) = chem(i,k,j,p_p25) + factor  * ebu(i,k,j,p_ebu_pm25)
              chem(i,k,j,p_p10) = chem(i,k,j,p_p10) + factor  * ebu(i,k,j,p_ebu_pm10)
              chem(i,k,j,p_so2) = chem(i,k,j,p_so2) + factor2 * ebu(i,k,j,p_ebu_so2 )
            end do
          end do
        end do
      end if

    end if

    ! -- subgrid convective transport
!    if (chem_conv_tr == CTRA_OPT_GRELL) then
!      call grelldrvct(dt,ktau,                  &
!        rho_phy,rcav,chem,tr_fall,              &
!        u_phy,v_phy,t_phy,moist,dz8w,p_phy,p8w, &
!        pbl,xlv,cp,grvity,rv,z_at_w,cu_co_ten,  &
!        numgas,chem_opt,                        &
!        num_chem,num_moist,tile,                &
!        ids,ide, jds,jde, kds,kde,              &
!        ims,ime, jms,jme, kms,kme,              &
!        its,ite, jts,jte, kts,kte)
!     endif


     !>-- compute dry deposition
     call dry_dep_driver(ktau,dt,julday,current_month,t_phy,p_phy,&
       moist,p8w,rmol,rri,gmt,t8w,rcav,                           &
       chem,rho_phy,dz8w,exch_h,hfx,                              &
       ivgtyp,tsk,gsw,vegfrac,pbl,ust,znt,zmid,z_at_w,            &
       xland,xlat,xlong,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,     &
       anh3,dry_fall,dep_vel_o3,grav,                             &
       e_co,kemit,snowh,numgas,                                   &
       num_chem,num_moist,                                        &
       ids,ide, jds,jde, kds,kde,                                 &
       ims,ime, jms,jme, kms,kme,                                 &
       its,ite, jts,jte, kts,kte)

     ! -- ls wet deposition
     select case (wetdep_ls_opt)
       case (WDLS_OPT_GSD)
         call wetdep_ls(dt,chem,rnav,moist,rho_phy,var_rmv,             &
                        num_moist,num_chem,p_qc,p_qi,dz8w,vvel,         &
                        ids,ide, jds,jde, kds,kde,                      &
                        ims,ime, jms,jme, kms,kme,                      &
                        its,ite, jts,jte, kts,kte)
       case (WDLS_OPT_NGAC)
         call WetRemovalGOCART(its,ite, jts,jte, kts,kte, 1,1, dt,      &
                               num_chem,var_rmv,chem,p_phy,t_phy,       &
                               rho_phy,dqdti,rcav,rnav,                 &
                               ims,ime, jms,jme, kms,kme)
         !if (chem_rc_check(localrc, msg="Failure in NGAC wet removal scheme", &
         !  file=__FILE__, line=__LINE__, rc=rc)) return
       case default
         ! -- no further option implemented
    end select

    if (call_gocart) then
      call gocart_chem_driver(ktau,dt,dtstep,gmt,julday,xcosz,          &
           t_phy,moist,chem,rho_phy,dz8w,p8w,backg_oh,oh_t,             &
           backg_h2o2,h2o2_t,backg_no3,no3_t,                           &
           dxy,grav,xlat,xlong,ttday,tcosz,                             &
           chem_opt,num_chem,num_moist,                                 &
           ids,ide, jds,jde, kds,kde,                                   &
           ims,ime, jms,jme, kms,kme,                                   &
           its,ite, jts,jte, kts,kte                        )
      call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,              &
           chem,rho_phy,dz8w,p8w,dxy,grav,                              &
           chem_opt,num_chem,num_moist,                                 &
           ids,ide, jds,jde, kds,kde,                                   &
           ims,ime, jms,jme, kms,kme,                                   &
           its,ite, jts,jte, kts,kte                        )
    endif

    if (call_radiation) then
      store_arrays = .false.
      select case (aer_ra_feedback)
        case (1)
          call optical_driver(curr_secs,dtstep,             &
               chem,dz8w,rri,relhum,                        &
               h2oai,h2oaj,                                 &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,    &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,         &
               num_chem,chem_opt,ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme,                   &
               its,ite, jts,jte, kts,kte)
          call aer_opt_out(aod,dz8w,                        &
               ext_coeff,bscat_coeff,asym_par,              &
               tauaersw,gaersw,waersw,tauaerlw,             &
               num_ext_coef,num_bscat_coef,num_asym_par,    &
               ids,ide, jds,jde, kds,kde,                   &
               ims,ime, jms,jme, kms,kme,                   &
               its,ite, jts,jte, kts,kte)
          call aer_ra(dz8w                                  &
               ,extt,ssca,asympar,nbands                    &
               ,tauaersw,gaersw,waersw,tauaerlw             &
               ,ids,ide, jds,jde, kds,kde                   &
               ,ims,ime, jms,jme, kms,kme                   &
               ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case (2)
          call aero_opt('sw',dz8w,chem                                  &
                   ,rri,relhum,aod                                      &
                   ,extt,ssca,asympar                                   &
!                  ,extt,ssca,asympar,num_chem                          &
                   ,ids,ide, jds,jde, kds,kde                           &
                   ,ims,ime, jms,jme, kms,kme                           &
                   ,its,ite, jts,jte, kts,kte)
          store_arrays = .true.
        case default
          ! -- no feedback
      end select
      if (store_arrays) then
        do nv = 1, nbands
          do k = kts, kte
            do i = its, ite
              ext_cof(i,k,nv) = extt   (i,k,1,nv)
              sscal  (i,k,nv) = ssca   (i,k,1,nv)
              asymp  (i,k,nv) = asympar(i,k,1,nv)
            end do
          end do
        end do
        aod2d(its:ite) = aod(its:ite,1)
      end if
    endif

    call sum_pm_gocart (                                                &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                       &
         num_chem,chem_opt,                                             &
         ids,ide, jds,jde, kds,kde,                                     &
         ims,ime, jms,jme, kms,kme,                                     &
         its,ite, jts,jte, kts,kte)

    ! -- pm25 and pm10 for output , not for tracer options
    do j = jts, jte
      do k = kts, kte
        do i = its, ite
          pm25  (i,j,k) = pm2_5_dry(i,k,j)
          p10   (i,j,k) = pm10     (i,k,j)
          ebu_oc(i,j,k) = ebu      (i,k,j,p_ebu_oc)
        end do
      end do
    end do

    if (call_gocart) then
      do j = jts, jte
        do k = kts, kte
          do i = its, ite
            oh_bg  (i,j,k) = max(0., oh_t  (i,k,j))
            h2o2_bg(i,j,k) = max(0., h2o2_t(i,k,j))
            no3_bg (i,j,k) = max(0., no3_t (i,k,j))
          end do
        end do
      end do
    end if


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntso2  )=ppm2ugkg(p_so2   ) * max(epsilc,chem(i,k,1,p_so2))
       gq0(i,k,ntsulf )=ppm2ugkg(p_sulf  ) * max(epsilc,chem(i,k,1,p_sulf))
       gq0(i,k,ntdms  )=ppm2ugkg(p_dms   ) * max(epsilc,chem(i,k,1,p_dms)) 
       gq0(i,k,ntmsa  )=ppm2ugkg(p_msa   ) * max(epsilc,chem(i,k,1,p_msa))
       gq0(i,k,ntpp25 )=ppm2ugkg(p_p25   ) * max(epsilc,chem(i,k,1,p_p25))
       gq0(i,k,ntbc1  )=ppm2ugkg(p_bc1   ) * max(epsilc,chem(i,k,1,p_bc1))
       gq0(i,k,ntbc2  )=ppm2ugkg(p_bc2   ) * max(epsilc,chem(i,k,1,p_bc2))
       gq0(i,k,ntoc1  )=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
       gq0(i,k,ntoc2  )=ppm2ugkg(p_oc2   ) * max(epsilc,chem(i,k,1,p_oc2))
       gq0(i,k,ntdust1)=ppm2ugkg(p_dust_1) * max(epsilc,chem(i,k,1,p_dust_1))
       gq0(i,k,ntdust2)=ppm2ugkg(p_dust_2) * max(epsilc,chem(i,k,1,p_dust_2))
       gq0(i,k,ntdust3)=ppm2ugkg(p_dust_3) * max(epsilc,chem(i,k,1,p_dust_3))
       gq0(i,k,ntdust4)=ppm2ugkg(p_dust_4) * max(epsilc,chem(i,k,1,p_dust_4))
       gq0(i,k,ntdust5)=ppm2ugkg(p_dust_5) * max(epsilc,chem(i,k,1,p_dust_5))
       gq0(i,k,ntss1  )=ppm2ugkg(p_seas_1) * max(epsilc,chem(i,k,1,p_seas_1))
       gq0(i,k,ntss2  )=ppm2ugkg(p_seas_2) * max(epsilc,chem(i,k,1,p_seas_2))
       gq0(i,k,ntss3  )=ppm2ugkg(p_seas_3) * max(epsilc,chem(i,k,1,p_seas_3))
       gq0(i,k,ntss4  )=ppm2ugkg(p_seas_4) * max(epsilc,chem(i,k,1,p_seas_4))
       gq0(i,k,ntss5  )=ppm2ugkg(p_seas_5) * max(epsilc,chem(i,k,1,p_seas_5))
       gq0(i,k,ntpp10 )=ppm2ugkg(p_p10   ) * max(epsilc,chem(i,k,1,p_p10))
     enddo
    enddo

!>---- feedback to radiation
    do nv = 1, nbands
     do k = kts, kte
      do i = its, ite
       faersw1(i,k,nv) =  ext_cof(i,k,nv)
       faersw2(i,k,nv) =  sscal  (i,k,nv)
       faersw3(i,k,nv) =  asymp  (i,k,nv)
      end do
     end do
    end do
!   print*,'hli rad',maxval(ext_cof),maxval(sscal),maxval(asymp)
!>---- feedback to radiation
!    print*,'hli test3 ktau',ktau
!   call gsd_chem_post() ! postprocessing for diagnostics

!
   end subroutine gsd_chem_driver_run
!> @}
  end module gsd_chem_driver
