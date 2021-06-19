!>\file gsd_chem_anthropogenic_wrapper.F90
!! This file is GSDChem anthropogenic emission wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 07/2020

 module gsd_chem_anthropogenic_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use gsd_chem_config

   implicit none

   private

   public :: gsd_chem_anthropogenic_wrapper_init, gsd_chem_anthropogenic_wrapper_run, gsd_chem_anthropogenic_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_anthropogenic_wrapper_init()
      end subroutine gsd_chem_anthropogenic_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_anthropogenic_wrapper_finalize Argument Table
!!
      subroutine gsd_chem_anthropogenic_wrapper_finalize()
      end subroutine gsd_chem_anthropogenic_wrapper_finalize

!> \defgroup gsd_chem_anthropogenic_group GSD Chem seas wrapper Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_anthropogenic_wrapper GSD Chem seas wrapper Module  
!> \ingroup gsd_chem_anthropogenic_group
!! This is the GSD Chem seas wrapper Module
!! \section arg_table_gsd_chem_anthropogenic_wrapper_run Argument Table
!! \htmlinclude gsd_chem_anthropogenic_wrapper_run.html
!!
!>\section gsd_chem_anthropogenic_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_anthropogenic_wrapper_run(im, kte, kme, ktau, dt,               &
                   pr3d, ph3d,phl3d, prl3d, tk3d, spechum,emi_in,                       &
                   ntrac,ntso2,ntsulf,ntpp25,ntbc1,ntoc1,ntpp10,                        &
                   gq0,qgrs,abem,chem_opt_in,kemit_in,pert_scale_anthro,                &
                   emis_amp_anthro,do_sppt_emis,sppt_wts,errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau
    integer,        intent(in) :: ntrac
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf
    real(kind_phys),intent(in) :: dt, emis_amp_anthro, pert_scale_anthro
    real(kind_phys), optional, intent(in) :: sppt_wts(:,:)
    logical,        intent(in) :: do_sppt_emis

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im, 10), intent(in) :: emi_in
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d, spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,7        ), intent(inout) :: abem
    integer,           intent(in) :: chem_opt_in, kemit_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,       &
                     p_phy, z_at_w, dz8w, p8w, rho_phy

!>- chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    integer :: ide, ime, ite, kde
    real(kind_phys), dimension(ims:im, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg
    real(kind_phys), parameter :: ugkg = 1.e-09_kind_phys !lzhang

    integer :: i, j, jp, k, kp, n
    real(kind_phys) :: random_factor(ims:im,jms:jme)
  

    errmsg = ''
    errflg = 0

    chem_opt          = chem_opt_in
    kemit             = kemit_in

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    if(do_sppt_emis) then
      random_factor(:,jms) = pert_scale_anthro*max(min(1+(sppt_wts(:,kme/2)-1)*emis_amp_anthro,2.0),0.0)
    else
      random_factor = 1.0
    endif

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

!>- get ready for chemistry run
    call gsd_chem_prep_anthropogenic(                                   &
        ktau,dtstep,                                                    &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,emi_in,                      &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,z_at_w,                        & 
        ntso2,ntsulf,ntpp25,ntbc1,ntoc1,ntpp10,ntrac,gq0,               &
        num_chem, num_ebu_in,num_emis_ant,                              &
        emis_ant,ppm2ugkg,chem,random_factor,                           &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntso2  )=ppm2ugkg(p_so2   ) * max(epsilc,chem(i,k,1,p_so2))
       gq0(i,k,ntsulf )=ppm2ugkg(p_sulf  ) * max(epsilc,chem(i,k,1,p_sulf))
       gq0(i,k,ntpp25 )=ppm2ugkg(p_p25   ) * max(epsilc,chem(i,k,1,p_p25))
       gq0(i,k,ntbc1  )=ppm2ugkg(p_bc1   ) * max(epsilc,chem(i,k,1,p_bc1))
       gq0(i,k,ntoc1  )=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
       gq0(i,k,ntpp10 )=ppm2ugkg(p_p10   ) * max(epsilc,chem(i,k,1,p_p10))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntso2  )=gq0(i,k,ntso2  )
       qgrs(i,k,ntsulf )=gq0(i,k,ntsulf )
       qgrs(i,k,ntpp25 )=gq0(i,k,ntpp25 )
       qgrs(i,k,ntbc1  )=gq0(i,k,ntbc1  )
       qgrs(i,k,ntoc1  )=gq0(i,k,ntoc1  )
       qgrs(i,k,ntpp10 )=gq0(i,k,ntpp10 )
     enddo
    enddo

    abem(:,1)=ugkg*emis_ant(:,kts,1,p_e_bc )
    abem(:,2)=ugkg*emis_ant(:,kts,1,p_e_oc )
    abem(:,3)=ugkg*emis_ant(:,kts,1,p_e_so2)

!
   end subroutine gsd_chem_anthropogenic_wrapper_run
!> @}

   subroutine gsd_chem_prep_anthropogenic(                               &
        ktau,dtstep,                                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,emi_in,                       &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,z_at_w,                         &
        ntso2,ntsulf,ntpp25,ntbc1,ntoc1,ntpp10,ntrac,gq0,                &
        num_chem, num_ebu_in,num_emis_ant,                               &
        emis_ant,ppm2ugkg,chem,random_factor,                            &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: ktau
    real(kind=kind_phys), intent(in) :: dtstep

    !Stochastic physics variables
    real(kind_phys), intent(in) :: random_factor(ims:ime,jms:jme)

    !FV3 input variables
    integer, intent(in) :: ntrac
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf
    real(kind=kind_phys), dimension(ims:ime,    10),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: phl3d,tk3d,prl3d,spechum
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_ebu_in,num_emis_ant
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

    real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, num_emis_ant), intent(inout) :: emis_ant
    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, p_phy, rho_phy, dz8w, p8w
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, jms:jme, num_ebu_in) :: emiss_ab
    real(kind_phys), parameter :: frac_so2_ant = 0.5_kind_phys  ! antropogenic so2 fraction

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2
    integer i,ip,j,jp,k,kp,kk,kkp,nv,l,n

    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys

    ! -- initialize fire emissions


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
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          !--
        enddo
      enddo
    enddo

    ! -- anthropagenic
    emiss_ab  = 0.   ! background
    do j=jts,jte
     do i=its,ite
      emiss_ab(i,j,p_e_bc)   =emi_in(i,1)*random_factor(i,j)
      emiss_ab(i,j,p_e_oc)   =emi_in(i,2)*random_factor(i,j)
      emiss_ab(i,j,p_e_sulf) =emi_in(i,3)*random_factor(i,j)
      emiss_ab(i,j,p_e_pm_25)=emi_in(i,4)*random_factor(i,j)
      emiss_ab(i,j,p_e_so2)  =emi_in(i,5)*random_factor(i,j)
      emiss_ab(i,j,p_e_pm_10)=emi_in(i,6)*random_factor(i,j)
     enddo
    enddo


    factor=0.
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
        enddo
      enddo
    endif

 
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_so2   )=max(epsilc,gq0(i,k,ntso2  )/ppm2ugkg(p_so2))
       chem(i,k,jts,p_sulf  )=max(epsilc,gq0(i,k,ntsulf )/ppm2ugkg(p_sulf))
       chem(i,k,jts,p_p25   )=max(epsilc,gq0(i,k,ntpp25 )/ppm2ugkg(p_p25))
       chem(i,k,jts,p_bc1   )=max(epsilc,gq0(i,k,ntbc1  )/ppm2ugkg(p_bc1))
       chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntoc1  )/ppm2ugkg(p_oc1))
       chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
     enddo
    enddo

    !
    ! -- gocart background fields only if gocart is called
    !

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


  end subroutine gsd_chem_prep_anthropogenic

!> @}
  end module gsd_chem_anthropogenic_wrapper
