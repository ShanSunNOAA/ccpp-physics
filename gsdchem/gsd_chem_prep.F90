 module gsd_chem_prep

  use machine ,       only : kind_phys
  use physcons, g => con_g, pi => con_pi
  use gsd_chem_config
  use plume_data_mod

  implicit none


contains

  subroutine gsd_chem_prep_run(                                        &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                     &
        u10m,v10m,ustar,land,garea,rlat,rlon,ts2d,                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,                &
        exch,dqdt,                                                     &
        nsoil,smc,vegtype,soiltyp,sigmaf,dswsfc,zorl,                  &
        snow_cpl,dust_in,emi_in,emi2_in,                               &
        fire_GBBEPx,fire_MODIS,hf2d,pb2d,                              &
        u10,v10,ust,tsk,xland,xlat,xlong,dxy,                          &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,exch_h,dqdti,                                              &
        z_at_w,vvel,zmid,                                              &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,                                       &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,gq0,                                                     &
        num_chem, num_moist,num_ebu_in,                                &
        call_gocart,nvl_gocart,                                        &
        ttday,tcosz,gmt,julday,                                        &
        backg_oh,backg_h2o2,backg_no3,                                 &
        plumerise_flag,num_plume_data,num_emis_ant,                    &
        emis_ant,ppm2ugkg,                                             &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,       &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,       &
        moist,chem,plumedist,ebu_in,                                   &
        smois,ivgtyp,isltyp,vegfrac,rmol,gsw,znt,hfx,pbl,              &
        relhum,snowh,clayf,rdrag,sandf,ssm,uthr,dms_0,erod,            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    logical, intent(in) :: readrestart
    integer, intent(in) :: chem_in_opt, ktau, julday
    real(kind=kind_phys), intent(in) :: dtstep, gmt

    !FV3 input variables
    integer, intent(in) :: nsoil
    integer, dimension(ims:ime), intent(in) :: land, vegtype, soiltyp
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         u10m, v10m, ustar, garea, rlat, rlon, ts2d, sigmaf, dswsfc,       &
         zorl, snow_cpl, hf2d, pb2d,xcosz
    real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc 
    real(kind=kind_phys), dimension(ims:ime,     5),   intent(in) :: dust_in
    real(kind=kind_phys), dimension(ims:ime,    10),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, 64, 3),   intent(in) :: emi2_in
    real(kind=kind_phys), dimension(ims:ime,     5),   intent(in) :: fire_GBBEPx
    real(kind=kind_phys), dimension(ims:ime,    13),   intent(in) :: fire_MODIS
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,exch,dqdt
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist, num_ebu_in,                &
                           plumerise_flag, num_plume_data, num_emis_ant,   &
                           nvl_gocart
    logical,intent(in) ::  call_gocart
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme),    intent(out) ::          &
                           backg_oh,backg_h2o2,backg_no3

    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in),intent(out) :: ebu_in
    real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), intent(inout) :: emis_ant
    
    integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
    real(kind_phys), dimension(ims:ime, jms:jme, 3), intent(inout) :: erod
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, dz8w, p8w, t8w, vvel, zmid,         &
         exch_h,dqdti
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         u10, v10, ust, tsk, xland, xlat, xlong, dxy, vegfrac, rmol, gsw, znt, hfx,    &
         pbl, snowh, clayf, rdrag, sandf, ssm, uthr, dms_0, ttday, tcosz
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w, relhum
    real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois
    real(kind_phys), dimension(ims:ime, jms:jme, num_frp_plume), intent(out) :: plumedist
    real(kind_phys), dimension(ims:ime, jms:jme   ), intent(out) ::                    &
                   mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,            &
                   firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr       
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_ab
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_abu
    real(kind_phys), dimension(ims:ime, jms:jme, num_plume_data) :: plume
    real(kind_phys), parameter :: frac_so2_ant = 0.5_kind_phys  ! antropogenic so2 fraction
    real(kind_phys), parameter :: frp2plume = 1.e+06_kind_phys  ! FRP-to-plume conversion factor
    real(kind_phys), parameter :: frpc  = 1.e+09_kind_phys      ! FRP conversion factor
    real(kind_phys), dimension(nvl_gocart) :: p_gocart

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys), dimension(ims:ime, jms:jme, nvl_gocart) :: oh_backgd,h2o2_backgd,no3_backgd
    real(kind_phys) ::  factor,factor2,pu,pl,aln,pwant
    real(kind_phys) ::  xhour,xmin,gmtp,xlonn,xtime,real_time
    real(kind_phys), DIMENSION (1,1) :: sza,cosszax
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour

    p_gocart = (/ 1000., 992.5, 985., 977.5, 970., 955., 940., 925., 910.,               &
          895., 880., 865., 850., 825., 800., 775., 750., 712.5,  675., 637.5, 600.,     &
          562.5, 525., 487.5, 450., 412.5, 375., 337.5, 288.08, 244.88, 208.15, 176.93,  &
          150.39, 127.84, 108.66, 92.37, 78.51, 66.6, 56.39, 47.64, 40.18, 33.81, 28.37, &
          23.73, 19.79,  16.46, 13.64, 11.28, 9.29, 7.62, 6.22, 5.05, 4.08, 3.28, 2.62,  &
          2.08, 1.65, 1.3, 1.02, 0.8, 0.62, 0.48, 0.37, 0.28 /)

    ! -- initialize output arrays
    backg_oh       = 0._kind_phys
    backg_h2o2     = 0._kind_phys 
    backg_no3      = 0._kind_phys
    ebu_in         = 0._kind_phys
    isltyp         = 0._kind_phys
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    zmid           = 0._kind_phys
    exch_h         = 0._kind_phys
    dqdti          = 0._kind_phys
    u10            = 0._kind_phys
    v10            = 0._kind_phys
    ust            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    dxy            = 0._kind_phys
    vegfrac        = 0._kind_phys
    rmol           = 0._kind_phys
    gsw            = 0._kind_phys
    znt            = 0._kind_phys
    hfx            = 0._kind_phys
    pbl            = 0._kind_phys
    snowh          = 0._kind_phys
    clayf          = 0._kind_phys
    rdrag          = 0._kind_phys
    sandf          = 0._kind_phys 
    ssm            = 0._kind_phys
    uthr           = 0._kind_phys
    dms_0          = 0._kind_phys
    ttday          = 0._kind_phys
    tcosz          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    relhum         = 0._kind_phys

    ! -- initialize fire emissions
    plume          = 0._kind_phys
    plumedist      = 0._kind_phys
    mean_fct_agtf  = 0._kind_phys
    mean_fct_agef  = 0._kind_phys
    mean_fct_agsv  = 0._kind_phys
    mean_fct_aggr  = 0._kind_phys
    firesize_agtf  = 0._kind_phys
    firesize_agef  = 0._kind_phys
    firesize_agsv  = 0._kind_phys
    firesize_aggr  = 0._kind_phys


    do i=its,ite
     u10  (i,1)=u10m (i)
     v10  (i,1)=v10m (i)
     tsk  (i,1)=ts2d (i)
     ust  (i,1)=ustar(i)
     dxy  (i,1)=garea(i)
     xland(i,1)=real(land(i))
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     gsw  (i,1)=dswsfc(i)
     znt  (i,1)=zorl(i)*0.01
     hfx  (i,1)=hf2d(i)
     pbl  (i,1)=pb2d(i)
     snowh(i,1)=snow_cpl(i)*0.001
     clayf(i,1)=dust_in(i,1)
     rdrag(i,1)=dust_in(i,2)
     sandf(i,1)=dust_in(i,3)
     ssm  (i,1)=dust_in(i,4)
     uthr (i,1)=dust_in(i,5)
     ivgtyp (i,1)=vegtype(i)
     isltyp (i,1)=soiltyp(i)
     vegfrac(i,1)=sigmaf (i)
     dms_0(i,1  )=emi_in(i,7) ! --dm0
     erod (i,1,1)=emi_in(i,8) ! --ero1
     erod (i,1,2)=emi_in(i,9) ! --ero2
     erod (i,1,3)=emi_in(i,10)! --ero3
    enddo
   
    rmol=0.

    do k=1,nsoil
     do j=jts,jte
      do i=its,ite
       smois(i,k,j)=smc(i,k)
      enddo
     enddo
    enddo

    if (ktau <= 1) then
      emis_ant = 0.
     !emis_vol = 0.
    end if

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          u_phy(i,k,j)=us3d(ip,kkp)
          dqdti(i,k,j)=dqdt(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          relhum(i,k,j) = .95
          relhum(i,k,j) = MIN( .95, moist(i,k,j,1) / &
            (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
            (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
          relhum(i,k,j)=max(0.1,relhum(i,k,j))
          !--
          zmid(i,k,j)=phl3d(ip,kkp)/g
        enddo
      enddo
    enddo

    ! -- the imported atmospheric heat diffusivity is only available up to kte-1
    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte-1
        kkp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          exch_h(i,k,j)=exch(ip,kkp)
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

    ! -- fire
    emiss_ab  = 0.   ! background
    emiss_abu = 0.   ! fire emission
    do j=jts,jte
     do i=its,ite
      emiss_ab(i,j,p_e_bc)   =emi_in(i,1)
      emiss_ab(i,j,p_e_oc)   =emi_in(i,2)
      emiss_ab(i,j,p_e_sulf) =emi_in(i,3)
      emiss_ab(i,j,p_e_pm_25)=emi_in(i,4)
      emiss_ab(i,j,p_e_so2)  =emi_in(i,5)
      emiss_ab(i,j,p_e_pm_10)=emi_in(i,6)
     enddo
    enddo

    !print*,'hli ',plumerise_flag,FIRE_OPT_MODIS,FIRE_OPT_GBBEPx
    select case (plumerise_flag)
      case (FIRE_OPT_MODIS)
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_MODIS(i,1)
          emiss_abu(i,j,p_e_oc)   =fire_MODIS(i,2)
          emiss_abu(i,j,p_e_pm_25)=fire_MODIS(i,3)
          emiss_abu(i,j,p_e_so2)  =fire_MODIS(i,4)
          emiss_abu(i,j,p_e_pm_10)=fire_MODIS(i,5)
          plume(i,j,1)            =fire_MODIS(i,6)
          plume(i,j,2)            =fire_MODIS(i,7)
          plume(i,j,3)            =fire_MODIS(i,8)
          plume(i,j,4)            =fire_MODIS(i,9)
          plume(i,j,5)            =fire_MODIS(i,10)
          plume(i,j,6)            =fire_MODIS(i,11)
          plume(i,j,7)            =fire_MODIS(i,12)
          plume(i,j,8)            =fire_MODIS(i,13)
         enddo
        enddo
      case (FIRE_OPT_GBBEPx)
        do j=jts,jte
         do i=its,ite
          emiss_abu(i,j,p_e_bc)   =fire_GBBEPx(i,1)
          emiss_abu(i,j,p_e_oc)   =fire_GBBEPx(i,2)
          emiss_abu(i,j,p_e_pm_25)=fire_GBBEPx(i,3)
          emiss_abu(i,j,p_e_so2)  =fire_GBBEPx(i,4)
          plume(i,j,1)            =fire_GBBEPx(i,5)
         enddo
        enddo
!        print*,'hli GBBEPx plume',maxval(plume(:,:,1))
      case default
          ! -- no further option available
    end select


    factor=0.
    jmax=0
    jmaxi=0
    k=kts
    if (p_bc2 > 1) then
      do j=jts,jte
        do i=its,ite
          emis_ant(i,k,j,p_e_bc)=emiss_ab(i,j,p_e_bc)
          emis_ant(i,k,j,p_e_oc)=emiss_ab(i,j,p_e_oc) + emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_sulf)=emiss_ab(i,j,p_e_sulf)
          emis_ant(i,k,j,p_e_so2)=frac_so2_ant * emiss_ab(i,j,p_e_so2)
          emis_ant(i,k,j,p_e_dms)= 0. !emiss_ab(j,p_e_dms)
          emis_ant(i,k,j,p_e_pm_25)=emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_pm_10)=emiss_ab(i,j,p_e_pm_10)


          ebu_in(i,j,p_ebu_in_pm10)=emiss_abu(i,j,p_e_pm_10)
          ebu_in(i,j,p_ebu_in_dms)= 0._kind_phys

          select case (plumerise_flag)
            case (FIRE_OPT_MODIS)
              ebu_in(i,j,p_ebu_in_oc)   = emiss_abu(i,j,p_e_oc)
              ebu_in(i,j,p_ebu_in_bc)   = emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = emiss_abu(i,j,p_e_pm_25)
              ebu_in(i,j,p_ebu_in_so2)  = emiss_abu(i,j,p_e_so2)
              mean_fct_agtf(i,j)=plume(i,j,1)
              mean_fct_agef(i,j)=plume(i,j,2)
              mean_fct_agsv(i,j)=plume(i,j,3)
              mean_fct_aggr(i,j)=plume(i,j,4)
              firesize_agtf(i,j)=plume(i,j,5)
              firesize_agef(i,j)=plume(i,j,6)
              firesize_agsv(i,j)=plume(i,j,7)
              firesize_aggr(i,j)=plume(i,j,8)
            case (FIRE_OPT_GBBEPx)
              ebu_in(i,j,p_ebu_in_oc)   = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc))
              ebu_in(i,j,p_ebu_in_bc)   = frpc * emiss_abu(i,j,p_e_bc)
              ebu_in(i,j,p_ebu_in_pm25) = frpc * (emiss_abu(i,j,p_e_pm_25) - emiss_abu(i,j,p_e_bc) - emiss_abu(i,j,p_e_oc))
              ebu_in(i,j,p_ebu_in_so2)  = frpc * emiss_abu(i,j,p_e_so2)
              plumedist(i,j,p_frp_flam_frac) = flaming(catb(ivgtyp(i,j)))
              plumedist(i,j,p_frp_mean     ) = frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_std      ) = 0.3_kind_phys   * frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_mean_size) = msize(ivgtyp(i,j)) * frp2plume * plume(i,j,1)
              plumedist(i,j,p_frp_std_size ) = 0.5_kind_phys * plumedist(i,j,p_frp_mean_size)
            case default
              ! -- no further option available
          end select
        enddo
      enddo
    endif
!    print*,'hli plumedist(:,:,p_frp_mean)',maxval(plumedist(:,:,p_frp_mean))

 
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_so2   )=max(epsilc,gq0(i,k,ntso2  )/ppm2ugkg(p_so2))
       chem(i,k,jts,p_sulf  )=max(epsilc,gq0(i,k,ntsulf )/ppm2ugkg(p_sulf))
       chem(i,k,jts,p_dms   )=max(epsilc,gq0(i,k,ntdms  )/ppm2ugkg(p_dms))
       chem(i,k,jts,p_msa   )=max(epsilc,gq0(i,k,ntmsa  )/ppm2ugkg(p_msa))
       chem(i,k,jts,p_p25   )=max(epsilc,gq0(i,k,ntpp25 )/ppm2ugkg(p_p25))
       chem(i,k,jts,p_bc1   )=max(epsilc,gq0(i,k,ntbc1  )/ppm2ugkg(p_bc1))
       chem(i,k,jts,p_bc2   )=max(epsilc,gq0(i,k,ntbc2  )/ppm2ugkg(p_bc2))
       chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntoc1  )/ppm2ugkg(p_oc1))
       chem(i,k,jts,p_oc2   )=max(epsilc,gq0(i,k,ntoc2  )/ppm2ugkg(p_oc2))
       chem(i,k,jts,p_dust_1)=max(epsilc,gq0(i,k,ntdust1)/ppm2ugkg(p_dust_1))
       chem(i,k,jts,p_dust_2)=max(epsilc,gq0(i,k,ntdust2)/ppm2ugkg(p_dust_2))
       chem(i,k,jts,p_dust_3)=max(epsilc,gq0(i,k,ntdust3)/ppm2ugkg(p_dust_3))
       chem(i,k,jts,p_dust_4)=max(epsilc,gq0(i,k,ntdust4)/ppm2ugkg(p_dust_4))
       chem(i,k,jts,p_dust_5)=max(epsilc,gq0(i,k,ntdust5)/ppm2ugkg(p_dust_5))
       chem(i,k,jts,p_seas_1)=max(epsilc,gq0(i,k,ntss1  )/ppm2ugkg(p_seas_1))
       chem(i,k,jts,p_seas_2)=max(epsilc,gq0(i,k,ntss2  )/ppm2ugkg(p_seas_2))
       chem(i,k,jts,p_seas_3)=max(epsilc,gq0(i,k,ntss3  )/ppm2ugkg(p_seas_3))
       chem(i,k,jts,p_seas_4)=max(epsilc,gq0(i,k,ntss4  )/ppm2ugkg(p_seas_4))
       chem(i,k,jts,p_seas_5)=max(epsilc,gq0(i,k,ntss5  )/ppm2ugkg(p_seas_5))
       chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
     enddo
    enddo

    if (.NOT. readrestart) then
      if (chem_in_opt == 0 ) then
        if(ktau.le.1)then
!           if(chem_opt > 0 ) then
          do j=jts,jte
            jp = j - jts + 1
            do k=kts,kte
              do i=its,ite
                ip = i - its + 1
                if (chem_opt == CHEM_OPT_GOCART) then
                  do n=1,num_chem
                    chem(i,k,j,n)=1.e-12
                  enddo
                endif  ! chem_opt==300
                chem(i,k,j,p_so2)=5.e-6
                chem(i,k,j,p_sulf)=3.e-6
                if ((chem_opt >= CHEM_OPT_GOCART) .and. (chem_opt < CHEM_OPT_MAX)) then
                  chem(i,k,j,p_msa)=0.1e-6
                  chem(i,k,j,p_dms)=0.1e-6
                  chem(i,k,j,p_bc1)=0.1e-3
                  chem(i,k,j,p_bc2)=0.1e-3
                  chem(i,k,j,p_oc1)=0.1e-3
                  chem(i,k,j,p_oc2)=0.1e-3
                  chem(i,k,j,p_p25)=0.1e-3 !lzhang
                  chem(i,k,j,p_p10)=0.1e-3 !lzhang
                endif !chem_opt >= 300 .and. chem_opt <  500

!                if ((chem_opt == CHEM_OPT_GOCART_RACM) .or. (chem_opt == CHEM_OPT_RACM_SOA_VBS)) then  !added o3 background !lzhang
!                  kk=min(k,kte)
!                  kkp = kk - kts + 1
!                  ! -- add initial constant into O3,CH4 and CO ect.
!                  chem(i,k,j,p_o3)=epsilc
!                  ! -- this section needs to be revisited before enabling the
!                  ! corresponding chem_opt options
!                  ! maxth=min(400.,th_pvsrf(i,j))
!                  ! if (tr3d(ip,jp,kkp,p_atm_ptem) > maxth) then
!                  !   chem(i,k,j,p_o3)=(airmw/48.)*tr3d(ip,jp,kkp,p_atm_o3mr)*1e6
!                  !   !convert kg/kg to ppm
!                  ! else
!                  !   chem(i,k,j,p_o3)=0.03 !ppm
!                  ! endif
!                  chem(i,k,j,p_ch4)=1.85 !ppm
!                  chem(i,k,j,p_co)=0.06 !ppm
!                  chem(i,k,j,p_co2)=380.
!                  chem(i,k,j,p_ete)=epsilc
!                  chem(i,k,j,p_udd)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_hket)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_api)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_lim)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_dien)=chem(i,k,j,p_ete)
!                  chem(i,k,j,p_macr)=chem(i,k,j,p_ete)
!                endif !( (chem_opt == 301.or.chem_opt==108))
              enddo
            enddo
          enddo
        endif !(ktau<=1)

      else !(chem_in_opt == 0 )

        if ((ktau<=1).and.((chem_opt == CHEM_OPT_GOCART_RACM).or.(chem_opt == CHEM_OPT_RACM_SOA_VBS))) then  !added GFS o3 background above 380K!lzhang
          do j=jts,jte
            jp = j - jts + 1
            do k=kts,kte+1
              kk=min(k,kte)
              kkp = kk - kts + 1
              do i=its,ite
                ip = i - its + 1
                ! -- this section needs to be revisited before enabling the
                ! corresponding chem_opt options
                ! maxth=min(400.,th_pvsrf(i,j))
                ! if (tr3d(ip,jp,kkp,p_atm_ptem) >= maxth) then
                !   chem(i,k,j,p_o3)=(airmw/48.)*tr3d(ip,jp,kkp,p_atm_o3mr)*1e6 !convert kg/kg to ppm
                ! endif !380K
              enddo
            enddo
          enddo
        endif ! chem_opt == 301.or.chem_opt==108

      endif !(chem_in_opt == 1 )
     endif ! readrestart

     !-- assgin read in 3D background chemical species
     do i=its,ite
       do k=1,nvl_gocart
          h2o2_backgd(i,1,k)=emi2_in(i,k,1)
          no3_backgd (i,1,k)=emi2_in(i,k,2)
          oh_backgd  (i,1,k)=emi2_in(i,k,3)
       enddo
     enddo

    !
    ! -- gocart background fields only if gocart is called
    !
    !if (.NOT. readrestart) then
    if (call_gocart .and. (chem_opt == CHEM_OPT_GOCART))then
      do j=jts,jte
        do i=its,ite
          do k=kts,kte
            do ll=2,nvl_gocart
              l=ll
              if (p_gocart(l) < .01*p_phy(i,k,j)) exit
            enddo
            pu=alog(p_gocart(l))
            pl=alog(p_gocart(l-1))
            pwant=alog(.01*p_phy(i,k,j))
            if (pwant > pl)then
              backg_oh(i,k,j)=oh_backgd(i,j,l)
              backg_h2o2(i,k,j)=h2o2_backgd(i,j,l)
              backg_no3(i,k,j)=no3_backgd(i,j,l)
            else
              aln=(oh_backgd(i,j,l)*(pwant-pl)+            &
                oh_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_oh(i,k,j)=aln
              aln=(h2o2_backgd(i,j,l)*(pwant-pl)+            &
                h2o2_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_h2o2(i,k,j)=aln
              aln=(no3_backgd(i,j,l)*(pwant-pl)+            &
                no3_backgd(i,j,l-1)*(pu-pwant))/(pu-pl)
              backg_no3(i,k,j)=aln
            endif
          enddo
        enddo
      enddo
    endif   ! end gocart stuff
    !endif !restart



!   emis_ant=0.
    nv=1
    k=kts
    factor2=0.
    factor=0.
    if (p_bc2 > 1)then
      if (chem_opt == CHEM_OPT_GOCART) then
        do j=jts,jte
          do i=its,ite
            factor=dtstep*rri(i,k,j)/dz8w(i,k,j)
            factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+emis_ant(i,k,j,p_e_bc)*factor
            chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+emis_ant(i,k,j,p_e_oc)*factor
            chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+emis_ant(i,k,j,p_e_pm_25)*factor
            chem(i,k,j,p_p10)=chem(i,k,j,p_p10)+emis_ant(i,k,j,p_e_pm_10)*factor
            chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+emis_ant(i,k,j,p_e_sulf)*factor
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+emis_ant(i,k,j,p_e_so2)*factor2
          enddo
        enddo
      endif
    else if (p_tr2 > 1)then    !co2 here
      do j=jts,jte
        do i=its,ite
!           factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
          !factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
          !chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr1)*factor2
          !chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
        enddo
      enddo
    else if ((p_tr2 > 1) .and. (p_bc2 > 1))then
      !call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
      !  file=__FILE__, line=__LINE__, rc=rc)
      return
    endif

    ! -- real-time application, keeping eruption constant
!
!    if (ktau <= 2) then
!      ! -- volcanic emissions
!      if (num_emis_vol > 0) then
!         !------------------
!      end if
!     end if

      if ((chem_opt == CHEM_OPT_RACM_SOA_VBS) .or. (chem_opt >= CHEM_OPT_GOCART .and. chem_opt < CHEM_OPT_MAX)) then
      !ndystep=86400/ifix(dtstepc)
      ndystep=86400/ifix(dtstep)
      do j=jts,jte
        do i=its,ite
          tcosz(i,j)=0.
          ttday(i,j)=0.
!         rlat=xlat(i,j)*3.1415926535590/180.
          xlonn=xlong(i,j)
          do n=1,ndystep
            xtime=n*dtstep/60.
            ixhour=ifix(gmt+.01)+ifix(xtime/60.)
            xhour=float(ixhour)
            xmin=60.*gmt+(xtime-xhour*60.)
            gmtp=mod(xhour,24.)
            gmtp=gmtp+xmin/60.
            CALL szangle(1, 1, julday, gmtp, sza, cosszax,xlonn,rlat(i))
            TCOSZ(i,j)=TCOSZ(I,J)+cosszax(1,1)
            if (cosszax(1,1) > 0.) ttday(i,j)=ttday(i,j)+dtstep
            !--use physics inst cosine zenith -- hli 03/06/2020
!            TCOSZ(i,j)=TCOSZ(I,J)+xcosz(i)
!            if (xcosz(i) > 0.) ttday(i,j)=ttday(i,j)+dtstep
          enddo
        enddo
      enddo
    endif !chem_opt >= 300 .and. chem_opt <  500


  end subroutine gsd_chem_prep_run


  SUBROUTINE szangle(imx, jmx, doy, xhour, sza, cossza,xlon,rlat)

!
! ****************************************************************************
! **                                                                        **
! **  This subroutine computes solar zenith angle (SZA):                    **
! **                                                                        **
! **      cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)         **
! **                                                                        **
! **  where LAT is the latitude angle, DEC is the solar declination angle,  **
! **  and AHR is the hour angle, all in radius.                             **
! **                                                                        **
! **  DOY = day-of-year, XHOUR = UT time (hrs).                             **
! **  XLON = longitude in degree, RLAT = latitude in radian.                **
! ****************************************************************************
!

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: imx, jmx
  INTEGER, INTENT(IN)    :: doy
  REAL(kind_phys),    INTENT(IN)    :: xhour
  REAL(kind_phys),    INTENT(OUT)   :: sza(imx,jmx), cossza(imx,jmx)

  REAL(kind_phys)    :: a0, a1, a2, a3, b1, b2, b3, r, dec, timloc, ahr,xlon,rlat
  real(kind_phys), parameter :: pi=3.14
  INTEGER :: i, j

  ! executable statements

  ! ***************************************************************************
  ! *  Solar declination angle:                                               *
  ! ***************************************************************************
  a0 = 0.006918
  a1 = 0.399912
  a2 = 0.006758
  a3 = 0.002697
  b1 = 0.070257
  b2 = 0.000907
  b3 = 0.000148
  r  = 2.0* pi * REAL(doy-1)/365.0
  !
  dec = a0 - a1*COS(  r)   + b1*SIN(  r)   &
           - a2*COS(2.0*r) + b2*SIN(2.0*r) &
           - a3*COS(3.0*r) + b3*SIN(3.0*r)
  !
  DO i = 1,imx
     ! ************************************************************************
     ! *  Hour angle (AHR) is a function of longitude.  AHR is zero at        *
     ! *  solar noon, and increases by 15 deg for every hour before or        *
     ! *  after solar noon.                                                   *
     ! ************************************************************************
     ! -- Local time in hours
     timloc  = xhour + xlon/15.0
     !      IF (timloc < 0.0) timloc = 24.0 + timloc
     IF (timloc > 24.0) timloc = timloc - 24.0
     !
     ! -- Hour angle
     ahr = ABS(timloc - 12.0) * 15.0 * pi/180.0

     DO j = 1,jmx
        ! -- Solar zenith angle      
        cossza(i,j) = SIN(rlat) * SIN(dec) + &
                      COS(rlat) * COS(dec) * COS(ahr)
        sza(i,j)    = ACOS(cossza(i,j)) * 180.0/pi
        IF (cossza(i,j) < 0.0)   cossza(i,j) = 0.0
        !
     END do
  END DO

  END SUBROUTINE szangle

 end module gsd_chem_prep
