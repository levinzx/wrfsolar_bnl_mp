  module module_ra_farms

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                 !
    ! Purpose: Couples the FARMS model to WRF         !
    !                                                 !
    ! Author: Yu Xie coded the FARMS model            !
    !         Pedro A. Jimenez coupled FARMS to WRF   !
    !                                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use module_model_constants, only : G

    implicit none

    private

    public :: Farms_driver

    real, parameter :: THREE_OVER_TWO = 3.0 / 2.0
    integer, parameter :: TAU_ICE_METHOD = 1
    logical, parameter :: USE_REST2 = .false.
    logical, parameter :: DEBUG = .false.

    real, parameter :: PTAU_MIN = 0.000001
    real, parameter :: DELTA_MIN = 0.000001
    real, parameter :: DE_ICE_MIN = 5.0, DE_ICE_MAX = 140.0
    real, parameter :: DE_CLOUD_MIN = 5.0, DE_CLOUD_MAX = 120.0
    real, parameter :: TAU_MIN = 0.0001, TAU_MAX = 300.0
    real, parameter :: AOD550_VAL = 0.12, ANGEXP_VAL = 1.3, AERSSA_VAL = 0.85, AERASY_VAL = 0.9
    real, parameter :: RE_CLOUD_CLIM = 8.E-6, RE_ICE_CLIM = 24.E-6, RE_SNOW_CLIM = 24.E-6

  contains

    subroutine Farms_driver (ims, ime, jms, jme, its, ite, jts, jte, kms, kme, kts, kte,  &
           p8w, rho, dz8w, albedo, aer_opt, aerssa2d, aerasy2d, aod5502d, angexp2d,  &
           coszen_loc, qv, qi, qs, qc, re_cloud, re_ice, re_snow,    &
           julian, swdown2, swddir2, swddni2, swddif2, swdownc2, swddnic2, &
           has_reqc, has_reqi, has_reqs, cldfra, swddir2_new, tauc_sw2)   ! X. Zhou: tauc_sw2 added

      implicit None

      integer, intent(in) :: ims, ime, jms, jme, its, ite, jts, jte, kms, kme, &
          kts, kte

      integer, intent(in) :: aer_opt, has_reqc, has_reqi, has_reqs
      real,    intent(in) :: julian

      real, dimension(ims:ime, jms:jme), intent(in) :: albedo, coszen_loc

      real, dimension(ims:ime, kms:kme, jms:jme ), intent(in) :: qv, qi, qs, qc, &
          p8w, rho, dz8w, cldfra
      real, dimension(ims:ime, kms:kme, jms:jme ), intent(in) :: re_cloud, re_ice, re_snow

      real, dimension(ims:ime, jms:jme), intent(inout) :: aerssa2d, aerasy2d, aod5502d, angexp2d
      real, dimension(ims:ime,jms:jme), intent(inout) :: swddir2, swdown2, &
                                                         swddif2, swddni2, &
                                                         swdownc2, swddnic2, &
                                                         swddir2_new, &
														 tauc_sw2 ! X. Zhou

        ! Local
      integer :: i, j, k
      real    :: tau_qv, tau_qi, tau_qs, pmw, swp, iwp, lwp, beta
      real    :: re_cloud_path, re_ice_path, re_snow_path, q_aux, cldfra_h
      real, dimension(kms:kme) :: rhodz, re_cloud_k, re_ice_k, re_snow_k


      j_loop: do j = jts, jte
         i_loop: do i = its, ite
           daytime_if: if (coszen_loc(i, j) <= 0.0 ) then
             swdown2(i, j) = 0.0
             swddni2(i, j) = 0.0
             swddir2(i, j) = 0.0
             swddif2(i, j) = 0.0
             swdownc2(i, j) = 0.0
             swddnic2(i, j) = 0.0
             swddir2_new(i, j) = 0.0
			 tauc_sw2(i, j) = 0.0 ! X. Zhou
           else
             rhodz(:) = rho(i, :, j) * dz8w(i, :, j) / (1. + qv(i, :, j))
             re_cloud_k(:) = re_cloud(i, :, j)
             re_ice_k(:) = re_ice(i, :, j)
             re_snow_k(:) = re_snow(i, :, j)

             if (has_reqc == 1) then
               do k = kts, kte
                 if (CLDFRA (i, k, j) > 0.0 .and. re_cloud_k(k) < 2.5E-6) re_cloud_k(k) = RE_CLOUD_CLIM
               end do
             else
               re_cloud_k(:) = RE_CLOUD_CLIM
             end if

             if (has_reqi == 1) then
               do k = kts, kte
                 if (cldfra(i, k, j) > 0.0 .and. re_ice_k(k) < 5.0E-6) re_ice_k(k) = RE_ICE_CLIM
               end do
             else
               re_ice_k(:) = RE_ICE_CLIM
             end if

             if (has_reqs == 1) then
               do k = kts, kte
                 if (cldfra(i, k, j) > 0.0 .and. re_snow_k(k) < 10.0E-6) re_snow_k(k) = RE_SNOW_CLIM
               end do
             else
               re_snow_k(:) = RE_SNOW_CLIM
             end if

               ! PMW
             pmw = integrate_1var (rhodz, qv(i, :, j), kms, kme, kts, kte)

               ! Calc effective radius water
             q_aux = integrate_1var (rhodz, qc(i, :, j), kms, kme, kts, kte)
             lwp = q_aux

             if (q_aux > 0.0) then
               re_cloud_path = integrate_2var (rhodz, qc(i, :, j), &
                   re_cloud_k, kms, kme, kts, kte)
               re_cloud_path = re_cloud_path / q_aux
             else
               re_cloud_path = 0.0
             end if

               ! Calc effective radius ice
             q_aux = integrate_1var (rhodz, qi(i, :, j), kms, kme, kts, kte)
             iwp = q_aux

             if (q_aux > 0.0) then
               re_ice_path = integrate_2var (rhodz, qi(i, :, j), &
                   re_ice_k, kms, kme, kts, kte)
               re_ice_path = re_ice_path / q_aux
             else
               re_ice_path = 0.0
             end if

               ! Calc effective radius snow
             q_aux = integrate_1var (rhodz, qs(i, :, j), kms, kme, kts, kte)
             swp = q_aux

             if (q_aux > 0.0) then
               re_snow_path = integrate_2var (rhodz, qs(i, :, j), &
                   re_snow_k, kms, kme, kts, kte)
               re_snow_path = re_snow_path / q_aux
             else
               re_snow_path = 0.0
             end if

               ! Calc horizontal cloud fraction 
             q_aux = integrate_1var (rhodz, qc(i, :, j) + qi(i, :, j) + qs(i, :, j), kms, kme, kts, kte)
             if (q_aux > 0.0) then
               cldfra_h = integrate_2var (rhodz, qc(i, :, j) + qi(i, :, j) + qs(i, :, j), &
                   cldfra(i, :, j), kms, kme, kts, kte)
               cldfra_h = cldfra_h / q_aux
             else
               cldfra_h = 0.0
             end if

               ! optical thickness  water
             if (re_cloud_path > 0.0) then
               tau_qv = THREE_OVER_TWO * lwp / re_cloud_path / 1000.0
             else
               tau_qv = 0.0
             end if

               ! Optical thickness ice 
             if (re_ice_path > 0.0) then
               if (TAU_ICE_METHOD == 1) then
                   ! Eq 10 in Matrosov et al. (2002)
                 tau_qi = iwp * 1000.0 * (0.02 + 4.2 / (2.0 * re_ice_path * 1.0e+6))
               else
                 tau_qi = iwp * 1000.0 * (-0.006656 + 3.686 / (2.0 * re_ice_path * 1.0e+6))
               end if
             else
               tau_qi = 0.0
             end if

               ! Optical thickness snow 
             if (re_snow_path > 0.0) then
               if (TAU_ICE_METHOD == 1) then
                   ! Eq 10 in Matrosov et al. (2002)
                 tau_qs = swp * 1000.0 * (0.02 + 4.2 / (2.0 * re_snow_path * 1.0e+6))
               else
                 tau_qs = swp * 1000.0 * (-0.006656 + 3.686 / (2.0 * re_snow_path * 1.0e+6))
               end if
             else
               tau_qs = 0.0
             end if

               ! Aerosols
             if (aer_opt == 1) then
                angexp2d(i, j) = ANGEXP_VAL
                aerssa2d(i, j) = AERSSA_VAL
                aerasy2d(i, j) = AERASY_VAL
             else if (aer_opt == 0) then
                aod5502d(i, j) = 0.0
                angexp2d(i, j) = 0.0
                aerssa2d(i, j) = 0.0
                aerasy2d(i, j) = 0.0
             end if

             beta = aod5502d(i, j) * (1000.0/ 550.0) ** (- angexp2d(i, j))

             Call Farms (p8w(i, 1, j), albedo(i, j), aerssa2d(i, j),         &
                 aerasy2d(i, j), coszen_loc(i, j), beta,                     &
                 angexp2d(i, j), pmw, tau_qv, tau_qi, tau_qs, cldfra_h,      &
                 re_cloud_path, re_ice_path, re_snow_path, int(julian),      &
                 swdown2(i, j), swddni2(i, j), swddif2(i, j), swddir2(i, j), &
                 swdownc2(i, j), swddnic2(i, j), swddir2_new(i, j),          &
				 tauc_sw2(i, j)) ! X. Zhou

          end if daytime_if
        end do i_loop
      end do j_loop

    end subroutine Farms_driver


    function Integrate_1var(rhodz, var1_1d, kms, kme, kts, kte) &
        result (return_value)
      
      implicit none

      integer, intent(in) :: kts, kte, kms, kme
      real, dimension(kms:kme), intent(in) :: var1_1d, rhodz

        ! Local
      real :: return_value
      integer :: k

      return_value = 0.0
      do k = kts, kte - 1
        return_value = return_value + var1_1d(k) * rhodz(k)
      end do

    end function Integrate_1var


    function Integrate_2var(rhodz, var1_1d, var2_1d, kms, kme, kts, kte) &
        result (return_value)

      implicit none

      integer, intent(in) :: kts, kte, kms, kme
      real, dimension(kms:kme), intent(in) :: var1_1d, var2_1d, rhodz

        ! Local
      real :: return_value
      integer :: k

      return_value = 0.0
      do k = kts, kte - 1
        return_value = return_value + var1_1d(k) * var2_1d(k) * rhodz(k)
      end do

    end function Integrate_2var


    subroutine FARMS (p_pa, albdo, ssa, g, solarangle, beta, alpha, w_mm, &
        tau_qv, tau_qi, tau_qs, cldfra_h, re_cloud_path_m, re_ice_path_m, re_snow_path_m, &
        juday, ghi, dni, dif, dir, ghi_clear, dni_clear, dir_new, tau_tot) ! X. Zhou: tau_tot added

        !!!!!! This Fast All-sky Radiation Model for Solar applications (FARMS) was developed by
        !!!!!! Yu Xie (Yu.Xie@nrel.gov). Please contact him for more
        !!!!!! information.

        ! Adapted by PAJ to couple it with WRF

        ! Input values
        ! p_pa: surface air pressure (Pa)
        ! albdo: surface albedo
        ! SSA: single-scattering albedo for aerosol
        ! g: asymmetric factor of aerosol
        ! Z: solar zenith angle 
        ! beta: Angstrom turbidity coeff., i.e. AOD at 1000 nm.
        ! alpha: Angstrom wavelength exponent
        !        For a cloudy sky, the suggested values are beta=0.1 and alpha=1.2,
        !        unless you have more accurate measurements.
        ! w: PWV (cm)
        ! tau_qv: cloud optical thickness (liquid water)
        ! tau_qi: cloud optical thickness (ice)
        ! tau_qs: cloud optical thickness (snow)
        ! re_cloud_path_m: effective radious liquid
        ! re_ice_path_m: effective radious ice
        ! re_snow_path_m: effective radious snow
        ! juday: day of year (1-366)
        ! ozone: ozone amount (1000DU, i.e. 0.4=400DU)
        ! phase: cloud thermodynamic phase (1=water, 2=ice)
        ! De: Cloud effective particle size (micron) 
        !     (For a water cloud, De=2*effective radius)
        !     (For an ice cloud, the defination of De follows Fu. 1996.
        !     Another useful defination is given by "Determination of ice cloud models using MODIS and MISR data" ) 


      implicit none
  
      real, intent(in) :: p_pa, albdo, ssa, g, solarangle, beta, alpha, w_mm, &
          tau_qv, tau_qi, tau_qs, re_cloud_path_m, re_ice_path_m, re_snow_path_m, cldfra_h
      integer, intent(in) :: juday
      real, intent(out) :: ghi, dni, dir, dif, ghi_clear, dni_clear, dir_new, tau_tot !X. Zhou: tau_tot added

        ! Local vars
      real :: de_cloud, de_cloud2, de_ice, de_ice2, de_snow, de_snow2, f0, f1, ftotal, ozone, &
          p, radius, Ruucld, pi, Ruucld_water, Ruucld_ice, Ruucld_snow, Ruuclr, tau_qv2, tau_qi2, &
          tau_qs2, tau, Tddcld0, Tddclr, Tducld, Tducld_ice, tducld_snow, tducld_water, Tduclr, Tuuclr, &
          ! w, Z, tau_tot, &
		  w, Z, & ! X. Zhou: tau_tot moved to out
          taudni, Tddcld1, Tddcld1_water, Tddcld1_ice, Tddcld1_snow, Fd, Fd0, Fd1, Fd2, Fd_all


      PI = acos(-1.0)

      p = p_pa / 100.0
      ozone = 0.265
      w = w_mm / 10.0

      de_cloud = 2.0 * re_cloud_path_m / 1.0e-6
      de_ice = 2.0 * re_ice_path_m / 1.0e-6
      de_snow = 2.0 * re_snow_path_m / 1.0e-6

      de_cloud = 2.0 * re_cloud_path_m / 1.0e-6
      de_cloud2 = Max (de_cloud, DE_CLOUD_MIN)
      de_cloud2 = Min (de_cloud2, DE_CLOUD_MAX)

      de_ice = 2.0 * re_ice_path_m / 1.0e-6
      de_ice2 = Max (de_ice, DE_ICE_MIN)
      de_ice2 = Min (de_ice2, DE_ICE_MAX)

      de_snow = 2.0 * re_snow_path_m / 1.0e-6
      de_snow2 = Max (de_snow, DE_ICE_MIN)
      de_snow2 = Min (de_snow2, DE_ICE_MAX)

      tau_qv2 = Max (tau_qv, TAU_MIN)
      tau_qv2 = Min (tau_qv2, TAU_MAX)
      tau_qi2 = Max (tau_qi, TAU_MIN)
      tau_qi2 = Min (tau_qi2, TAU_MAX)
      tau_qs2 = Max (tau_qs, TAU_MIN)
      tau_qs2 = Min (tau_qs2, TAU_MAX)

      tau_tot = 0.0
      if (tau_qv > TAU_MIN) tau_tot = tau_tot + tau_qv2
      if (tau_qi > TAU_MIN) tau_tot = tau_tot + tau_qi2
      if (tau_qs > TAU_MIN) tau_tot = tau_tot + tau_qs2
      tau_tot = Min (tau_tot, TAU_MAX)

      Z = acos(solarangle) * 180.0 / PI

      call SUNEARTH( juday, Radius )
      F0 = 1361.2/(Radius*Radius)

      call CLEARSKYALL(p, albdo, SSA, g, Z, Radius, beta,&
          alpha, ozone, w, Tddclr, Tduclr, Ruuclr, Tuuclr)

        ! Liquid hydrometeors
      if (tau_qv > TAU_MIN) then
        Call Watermodel (tau_tot, de_cloud2, solarangle, Tducld_water, Ruucld_water)
      else
        Tducld_water = 1.0
        Ruucld_water = 0.0
        tau_qv2 = 0.0
      end if

        ! Ice hydrometeors
      if (tau_qi > TAU_MIN) then
        call Icemodel (tau_tot, de_ice2, solarangle, Tducld_ice, Ruucld_ice)
      else
        Tducld_ice = 1.0
        Ruucld_ice = 0.0
        tau_qi2 = 0.0
      end if

        ! Snow hydrometeors
      if (tau_qs > TAU_MIN) then
        call Icemodel (tau_tot, de_snow2, solarangle, Tducld_snow, Ruucld_snow)
      else
        Tducld_snow = 1.0
        Ruucld_snow = 0.0
        tau_qs2 = 0.0
      end if

      if (tau_tot > 0.0) then
        Tducld =  (tau_qv2 * Tducld_water + tau_qi2 * Tducld_ice + tau_qs2 * Tducld_snow) / tau_tot
        Tducld = min(Tducld, 1.0)
        Tducld = max(Tducld, 0.0)

        Ruucld = (tau_qv2 * Ruucld_water + tau_qi2 * Ruucld_ice + tau_qs2 * Ruucld_snow)  / tau_tot
        Ruucld = min(Ruucld, 1.0) 
        Ruucld = max(Ruucld, 0.0) 
      end if

!      tau = tau_qv2 + tau_qi2 + tau_qs2 
!      tau = Min (tau, TAU_MAX)
!      tau = Max (tau, TAU_MIN)

      Tddcld0 = exp(-tau_tot/solarangle)

      dni_clear = F0 * Tddclr
      dni = dni_clear * Tddcld0
      dir = dni * solarangle 
      Fd0 = solarangle*F0*Tddcld0*Tddclr

      F1 = solarangle*F0*( Tddcld0*(Tddclr+Tduclr) + Tducld*Tuuclr )

      Ftotal = F1/(1.0-albdo*(Ruuclr + Ruucld*Tuuclr*Tuuclr))
      ghi = Ftotal
      ghi_clear =  solarangle * F0 * ((Tddclr + Tduclr) / (1.0 - albdo * Ruuclr))

      ghi = cldfra_h * ghi + (1.0 - cldfra_h) * ghi_clear
      dni = cldfra_h * dni + (1.0 - cldfra_h) * dni_clear
      
	  if (DEBUG) print *, 'F0 =', F0
      if (DEBUG) print *, 'F1 =', F1
      if (DEBUG) print *, 'Albedo =', albdo
      if (DEBUG) print *, 'Ruuclr =', Ruuclr
      if (DEBUG) print *, 'Tuuclr =', Tuuclr
      if (DEBUG) print *, 'Reflection factor = ', (1.0-albdo*(Ruuclr + Ruucld*Tuuclr*Tuuclr))
	  
      dif = ghi - dir
	  
	  if (DEBUG) print *, 'dni = ', dni
      if (DEBUG) print *, 'dir = ', dir
      if (DEBUG) print *, 'dif = ', dif
      if (DEBUG) print *, 'ghi = ', ghi

      !!!! New FARMS dir
        ! Liquid hydrometeors
      if (tau_qv > TAU_MIN) then
         call TAUDNIWATER( tau_tot, taudni)
         call Pwater( taudni, Z, de_cloud2, Tddcld1_water )
      else
         Tddcld1_water = 1.0
      end if

        ! Ice hydrometeors
      if (tau_qi > TAU_MIN) then
         call TAUDNIICE( tau_tot, taudni )
         call Pice( taudni, Z, de_ice2, Tddcld1_ice )
      else
         Tddcld1_ice = 1.0
      end if

        ! Snow hydrometeors
      if (tau_qs > TAU_MIN) then
         call TAUDNIICE( tau_tot, taudni )
         call Pice( taudni, Z, de_snow2, Tddcld1_snow )
      else
        Tddcld1_snow = 1.0
      end if

     if (tau_tot > 0.0) then
        Tddcld1 =  (tau_qv2 * Tddcld1_water + tau_qi2 * Tddcld1_ice + tau_qs2 * Tddcld1_snow) / tau_tot
        Tddcld1 = min(Tddcld1, 1.0)
        Tddcld1 = max(Tddcld1, 0.0)
     end if


      Fd1 = solarangle*F0*Tddclr*Tddcld1

      call TDD2(Z, Ftotal, F1, Fd2)
      Fd = Fd0 + Fd1 + Fd2

      if (tau_tot .eq. 0.0 ) then
         Fd = 0.0
      end if

      Fd_all = cldfra_h * Fd + (1.0 - cldfra_h) * dni_clear*solarangle

      dir_new = Fd_all

    end subroutine farms

    subroutine TAUDNIWATER( tau, taudni )

       implicit none

       real, intent (in) :: tau
       real, intent (out) :: taudni

       if (tau.lt.8.0) then 
          taudni = (0.254825*tau - 0.00232717*tau*tau + (5.19320e-06)*tau**3.0)*(1.0+(8.0-tau)*0.07 )
       else 
          taudni = 0.2*(tau-8.0)**1.5 + 2.10871
       end if
    end

    subroutine TAUDNIICE(tau, taudni)

       implicit none

       real, intent (in) :: tau
       real, intent (out) :: taudni

       if (tau.lt.8.0) then
          taudni = (0.345353*tau - 0.00244671*tau*tau + (4.74263e-06)*tau**3.0)
       else
          taudni = 0.2*(tau-8.0)**1.5 + 2.91345
       end if
    end

    subroutine TDD2( Z, Ftotal, F1, Fd2)
     
       implicit none

       real, intent (in) :: Z, Ftotal, F1
       real, intent (out) :: Fd2
       real*16  pi, a, b, c, muomega

      PI = acos(-1.0)
      a = 5.94991536e-03
      b = 5.42116600e-01
      c = 331280.9859904468
      muomega = exp(-((Z-b)**3.0)/c )
      muomega = a*muomega
      Fd2 = (Ftotal - F1)*muomega/PI

    end subroutine TDD2

    subroutine Pwater(tau, Z, De, Tddcld2)
       implicit none

       real, intent (in) :: tau, Z, De
       real, intent (out) :: Tddcld2 

       real :: pi, umu0
       real :: taup, Tddp, h, a, b
 
       PI = acos(-1.0)
       umu0 = cos(Z*PI/180.0)

       if (De .lt. 10.0) then
       if (umu0 .lt. 0.1391) then
       taup = 0.1
       endif
       if (umu0 .ge. 0.1391 .and. umu0 .lt. 0.2419) then 
       taup = 0.2
       endif
       if (umu0 .ge. 0.2419 .and. umu0 .lt. 0.3090) then
       taup = 0.3
       endif
       if (umu0 .ge. 0.3090 .and. umu0 .lt. 0.4067) then
       taup = 0.4
       endif
       if (umu0 .ge. 0.4067 .and. umu0 .lt. 0.6156) then 
       taup = 0.5
       endif
       if (umu0 .ge. 0.6156) then 
       taup = 1.0
       endif
       endif

       if (De .ge. 10) then 
       if (umu0 .lt. 0.1391) then 
       taup = 0.1
       endif
       if (umu0 .ge. 0.1391 .and. umu0 .lt. 0.2079) then 
       taup = 0.2
       endif
       if (umu0 .ge. 0.2079 .and. umu0 .lt. 0.3090) then 
       taup = 0.3
       endif
       if (umu0 .ge. 0.3090 .and. umu0 .lt. 0.3746) then 
       taup = 0.4
       endif
       if (umu0 .ge. 0.3746 .and. umu0 .lt. 0.6156) then 
       taup = 0.5
       endif
       if (umu0 .ge. 0.6156) then 
       taup = 1.0
       endif
       endif


       !!!Tddp
       h = 0.005553*log(De) + 0.002503

       if ( umu0 .ge. 0.0 .and. umu0 .lt. 0.342) then 
       Tddp = h*(-0.1787*umu0*umu0 + 0.2207*umu0 + 0.977)
       endif
       if ( umu0 .ge. 0.342 .and. umu0 .lt. 0.4694) then 
       Tddp = h
       endif
       if ( umu0 .ge. 0.4694 .and. umu0 .lt. 0.7193) then 
       Tddp = h*( 2.6399*umu0*umu0 - 3.2111*umu0 + 1.9434 )
       endif
       if ( umu0 .ge. 0.7193 .and. umu0 .lt. 0.8829) then 
       Tddp = h*( -0.224*umu0*umu0 + 0.0835*umu0 + 1.056 )
       endif
       if ( umu0 .ge. 0.8829 .and. umu0 .lt. 0.9396) then 
       Tddp = h*( -94.381*umu0*umu0 + 170.32*umu0 - 75.843 )
       endif
       if ( umu0 .ge. 0.9396 .and. umu0 .lt. 0.9945) then
       Tddp = h*( -12.794*umu0*umu0 + 22.686*umu0 - 8.9392 )
       endif
       if ( umu0 .ge. 0.9945 .and. umu0 .lt. 0.9999) then 
       Tddp = h*( 11248.61*umu0*umu0 - 22441.07*umu0 + 11193.59 )
       endif
       if ( umu0 .ge. 0.9999 ) then 
       Tddp = 0.76*h
       endif

       !!!Tddcld2
       a = 2.0339*umu0**(-0.927)
       b = 6.6421*umu0**(2.0672)

       if ( tau .le. 0.9*taup ) then
       Tddcld2 = Tddp*tanh(a*tau)
       endif
       if ( tau .gt. 0.9*taup .and. tau .lt. taup ) then
       Tddcld2 = Tddp*tanh(0.9*a*taup) + &
             Tddp*(tanh(b/taup**2.0)-tanh(0.9*a*taup))*(tau-0.9*taup)/(0.1*taup)
       endif
       if ( tau .ge. taup ) then
       Tddcld2 = Tddp*tanh(b/tau**2.0)
       endif

    end

    subroutine Pice(tau, Z, De, Tddcld2)
       implicit none

       real, intent (in) :: tau, Z, De
       real, intent (out) :: Tddcld2 

       real :: pi, umu0
       real :: taup, Tddp, a, b, ade, bde

      PI = acos(-1.0)
      umu0 = cos(Z*PI/180.0)

      !!!!taup
      if ( De .ge. 5.0 .and. De .lt. 14.0 ) then

      if (umu0 .lt. 0.1391) then
      taup = 0.1
      endif
      if (umu0 .ge. 0.1391 .and. umu0 .lt. 0.2079) then
      taup = 0.2
      endif
      if (umu0 .ge. 0.2079 .and. umu0 .lt. 0.3090) then
      taup = 0.3
      endif
      if (umu0 .ge. 0.3090 .and. umu0 .lt. 0.3746) then
      taup = 0.4
      endif
      if (umu0 .ge. 0.3746 .and. umu0 .lt. 0.6156) then
      taup = 0.5
      endif
      if (umu0 .ge. 0.6156 .and. umu0 .lt. 0.9994) then
      taup = 1.0
      endif
      if (umu0 .ge. 0.9994) then
      taup = 1.5
      endif

      endif

      if ( De .ge. 14.0 .and. De .lt. 50.0 ) then

      if (umu0 .lt. 0.139173) then
      taup = 0.1
      endif
      if (umu0 .ge. 0.139173 .and. umu0 .lt. -0.0011*De+0.2307) then
      taup = 0.2
      endif
      if (umu0 .ge. -0.0011*De+0.2307 .and. umu0 .lt. -0.0022*De+0.3340) then
      taup = 0.3
      endif
      if (umu0 .ge. -0.0022*De+0.3340 .and. umu0 .lt. -0.0020*De+0.4096) then
      taup = 0.4
      endif
      if (umu0 .ge. -0.0020*De+0.4096 .and. umu0 .lt. -0.0033*De+0.6461) then
      taup = 0.5
      endif
      if (umu0 .ge. -0.0033*De+0.6461 .and. umu0 .lt. -0.0049*De+1.0713) then
      taup = 1.0
      endif
      if (umu0 .ge. -0.0049*De+1.0713) then
      taup = 1.5
      endif

      endif


      if ( De .ge. 50.0 ) then 

      if (umu0 .lt. -0.0006*De+0.2109) then
      taup = 0.2
      endif
      if (umu0 .ge. -0.0006*De+0.2109 .and. umu0 .lt. -0.0005*De+0.2581) then
      taup = 0.3
      endif
      if (umu0 .ge. -0.0005*De+0.2581 .and. umu0 .lt. -0.0010*De+0.3907) then
      taup = 0.4
      endif
      if (umu0 .ge. -0.0010*De+0.3907 .and. umu0 .lt. -0.0008*De+0.4900) then
      taup = 0.5
      endif
      if (umu0 .ge. -0.0008*De+0.4900 .and. umu0 .lt. -0.0017*De+0.8708) then
      taup = 1.0
      endif
      if (umu0 .ge. -0.0017*De+0.8708 .and. umu0 .lt. -0.0006*De+1.0367) then
      taup = 1.5
      endif
      if (umu0 .ge. -0.0006*De+1.0367) then
      taup = 2.0
      endif

      endif


      !!!Tddp
      if ( umu0 .ge. 0.9994 ) then
  
      if ( De .le. 10.0 ) then
      Tddp = 0.12269
      endif
      if ( De .gt. 10.0 .and. De .le. 16.0 ) then
      Tddp = 0.0015*De + 0.1078
      endif
      if ( De .gt. 16.0 ) then
      Tddp = 0.1621*exp(-0.016*De)
      endif

      endif

      if ( umu0 .lt. 0.9994 ) then 

      if ( De .le. 10.0 ) then
      if ( umu0 .lt. 0.9396 ) then 
      Tddp = 0.14991
      endif
      if ( umu0 .ge. 0.9396 .and. umu0 .lt. 0.9945) then
      Tddp = -4.5171*umu0**2.0 + 8.3056*umu0 - 3.6476
      endif
      if ( umu0 .ge. 0.9945 ) then 
      Tddp = 298.45*umu0**2.0 - 601.33*umu0 + 303.04
      endif
      endif



      if ( De .gt. 10.0 .and. De .le. 30.0) then
      ade = -0.000232338*De**2.0 + 0.012748726*De + 0.046745083
      if ( umu0 .le. 0.2419 ) then 
      Tddp = (-8.454*umu0**2.0 + 2.4095*umu0 +0.8425)*ade
      endif
      if ( umu0 .gt. 0.2419 .and. umu0 .le. 0.3746) then 
      Tddp = (-13.528*umu0**2.0 + 7.8403*umu0 -0.1221)*ade
      endif
      if ( umu0 .gt. 0.3746 .and. umu0 .le. 0.4694) then 
      Tddp = (19.524*umu0**2.0 - 16.5*umu0 +4.4612)*ade
      endif
      if ( umu0 .gt. 0.4694 .and. umu0 .le. 0.5877) then
      Tddp = (16.737*umu0**2.0 - 17.419*umu0 +5.4881)*ade
      endif
      if ( umu0 .gt. 0.5877 .and. umu0 .le. 0.6691) then
      Tddp = (-39.493*umu0**2.0 + 48.963*umu0 -14.175)*ade
      endif
      if ( umu0 .gt. 0.6691 .and. umu0 .le. 0.7660) then 
      Tddp = (0.4017*umu0**2.0 - 0.243*umu0 +0.9609)*ade
      endif
      if ( umu0 .gt. 0.7660 .and. umu0 .le. 0.8480) then 
      Tddp = (-11.183*umu0**2.0 + 18.126*umu0 -6.3417)*ade
      endif
      if ( umu0 .gt. 0.8480 .and. umu0 .le. 0.8987) then 
      Tddp = (-163.36*umu0**2.0 + 283.35*umu0 -121.91)*ade
      endif
      if ( umu0 .gt. 0.8987 .and. umu0 .le. 0.9396) then 
      Tddp = (-202.72*umu0**2.0 + 368.75*umu0 -166.75)*ade
      endif
      if ( umu0 .gt. 0.9396 .and. umu0 .le. 0.9702) then
      Tddp = (-181.72*umu0**2.0 + 343.59*umu0 -161.3)*ade
      endif
      if ( umu0 .gt. 0.9702 .and. umu0 .le. 0.9945) then 
      Tddp = (127.66*umu0**2.0 - 255.73*umu0 +129.03)*ade
      endif
      if ( umu0 .gt. 0.9945 ) then 
      Tddp = (908.66*umu0**2.0 - 1869.3*umu0 +961.63)*ade
      endif
      endif


      if ( De .gt. 30.0 ) then
      bde = 0.0000166112*De**2.0 - 0.00410998*De + 0.352026619
      if ( umu0 .le. 0.2419 ) then 
      Tddp = (-4.362*umu0**2.0 - 0.0878*umu0 +1.1218)*bde
      endif
      if ( umu0 .gt. 0.2419 .and. umu0 .le. 0.3746) then
      Tddp = (-49.566*umu0**2.0 + 28.767*umu0 -3.1299)*bde
      endif
      if ( umu0 .gt. 0.3746 .and. umu0 .le. 0.4694) then
      Tddp = (58.572*umu0**2.0 - 49.5*umu0 +11.363)*bde
      endif
      if ( umu0 .gt. 0.4694 .and. umu0 .le. 0.5877) then
      Tddp = (62.118*umu0**2.0 - 63.037*umu0 +16.875)*bde
      endif
      if ( umu0 .gt. 0.5877 .and. umu0 .le. 0.6691) then 
      Tddp = (-237.68*umu0**2.0 + 293.21*umu0 -89.328)*bde
      endif
      if ( umu0 .gt. 0.6691 .and. umu0 .le. 0.7660) then
      Tddp = (1.2051*umu0**2.0 - 0.7291*umu0 +0.8826)*bde
      endif
      if ( umu0 .gt. 0.7660 .and. umu0 .le. 0.8480) then 
      Tddp = (-55.6*umu0**2.0 + 90.698*umu0 -35.905)*bde
      endif
      if ( umu0 .gt. 0.8480 .and. umu0 .le. 0.8987) then
      Tddp = (-422.36*umu0**2.0 + 733.97*umu0 -317.89)*bde
      endif
      if ( umu0 .gt. 0.8987 .and. umu0 .le. 0.9396) then
      Tddp = (-457.09*umu0**2.0 + 831.11*umu0 -376.85)*bde
      endif
      if ( umu0 .gt. 0.9396 .and. umu0 .le. 0.9702) then
      Tddp = (-344.91*umu0**2.0 + 655.67*umu0 -310.5)*bde
      endif
      if ( umu0 .gt. 0.9702 .and. umu0 .le. 0.9945) then 
      Tddp = (622.85*umu0**2.0 - 1227.6*umu0 +605.97)*bde
      endif
      if ( umu0 .gt. 0.9945 ) then 
      Tddp = (6309.63*umu0**2.0 - 12654.78*umu0 +6346.15)*bde
      endif
      endif

      endif


      !!!Tddcld2
      a = 1.7686*umu0**(-0.95)
      b = 7.117*umu0**(1.9658)

      if ( tau .le. 0.9*taup ) then
      Tddcld2 = Tddp*tanh(a*tau)
      endif
      if ( tau .gt. 0.9*taup .and. tau .lt. taup ) then
      Tddcld2 = Tddp*tanh(0.9*a*taup) + &
      Tddp*(tanh(b/taup**2.0)-tanh(0.9*a*taup))*(tau-0.9*taup)/(0.1*taup)
      endif
      if ( tau .ge. taup ) then
      Tddcld2 = Tddp*tanh(b/tau**2.0)
      endif


    end

    subroutine SUNEARTH( juday, R )

      implicit none

      integer, intent(in) :: juday
      real, intent(out) :: R

      real :: pi, b, R1


      PI = acos(-1.0)

      b = 2.0*PI*juday/365.0
      R1 = 1.00011 + 0.034221*cos(b) + 0.001280*sin(b) + &
           0.000719*cos(2.0*b) +0.000077*sin(2.0*b)

      R = R1**(-0.5)

    end subroutine SUNEARTH


    subroutine CLEARSKYALL( p, albdo, SSA, g, Z, Radius, beta,&
        alpha, ozone, w, Tddclr, Tduclr, Ruuclr, Tuuclr )

        implicit none

        real, intent(in) :: p, albdo, SSA, g, Z, Radius, beta, alpha, ozone, w
        real, intent(out) :: Tddclr, Tduclr, Ruuclr, Tuuclr

          ! Local vars
        integer, parameter :: nangle = 10
        integer :: i
        real :: mu(nangle), angle(nangle), aa(nangle)
        real :: PI


        PI = acos(-1.0)

        do i=1, nangle
          mu(i) = (i-1.0)*0.1+0.1
          angle(i) = acos(mu(i))*180.0/PI
        end do

        do i=1, nangle
            ! PAJ: can use both REST2 and BIRD
          if (USE_REST2) then
            call CLEARSKY (p, albdo, SSA, g, angle(i), Radius, beta,&
            alpha, ozone, w, Tddclr, Tduclr, Ruuclr)
          else
            call BIRD(p, albdo, SSA, g, angle(i), Radius, beta,&
            alpha, ozone, w, Tddclr, Tduclr, Ruuclr)
          end if
          aa(i) = Tddclr
        end do

        Tuuclr = 0.0
        do i=1, nangle
          Tuuclr = Tuuclr + 2.0*mu(i)*aa(i)*0.1
        end do 
 
          ! PAJ: can use both REST2 and BIRD
        if (USE_REST2) then
          call CLEARSKY(p, albdo, SSA, g, Z, Radius, beta,&
          alpha, ozone, w, Tddclr, Tduclr, Ruuclr)
        else
          call BIRD(p, albdo, SSA, g, Z, Radius, beta,&
          alpha, ozone, w, Tddclr, Tduclr, Ruuclr)
        end if
      
    end subroutine CLEARSKYALL 


    subroutine CLEARSKY(p, albdo, SSA, g, Z, Radius, beta,&
        alpha, ozone, w, Tddclr, Tduclr, Ruuclr)

      implicit none

      real, intent(in) :: p, albdo, SSA, g, Z, Radius, beta, alpha, ozone, w
      real, intent(out) :: Tddclr, Tduclr, Ruuclr
     

      return

    end subroutine CLEARSKY


    subroutine BIRD(p, albdo, SSA, g, Z, Radius, beta,&
        alpha, ozone, w, Tddclr, Tduclr, Ruuclr)

        ! This clear-sky model follows the equations given by Bird (1981)
        ! This subroutine for the all-sky fast model is given by
        ! Yu Xie (Yu.Xie@nrel.gov)

      implicit none

      real, intent (in) :: p, albdo, SSA, g, Z, Radius, beta, alpha, ozone, w
      real, intent (out) :: Tddclr, Tduclr, Ruuclr

        ! Local vars
      real :: degrad, airmass, airmassp, T_rayleigh, x0, T_o, T_gases, xw, T_water, &
          tau038, tau050, taua, T_aerosol, T_AA, T_AS, F0, Fdif, Ftotal, Fddclr


      degrad=.017453293d+00
      airmass = 1/(cos(Z*degrad)+0.15*(93.885-Z)**(-1.25))
      airmassp = p*airmass/1013.0

      T_rayleigh = exp(-0.0903*(airmassp**0.84)*(1.0 + airmassp -&
          airmassp**1.01))

      x0 = ozone*airmass
      T_o = 1.0 - 0.1611*x0*(1.0+139.48*x0)**(-0.3035) -&
          0.002715*x0/(1.0+0.044*x0 + 0.0003*x0**2.0)
      T_gases = exp( -0.0127*(airmassp**0.26) )
      xw = w*airmass
      T_water = 1.0 - 2.4959*xw/( (1.0+79.034*xw)**0.6828 + 6.385*xw )

      tau038 = beta*(0.38**(-alpha))
      tau050 = beta*(0.50**(-alpha))
      taua = 0.2758*tau038 + 0.35*tau050

      T_aerosol = exp( -(taua**0.873)*(1.0+taua-taua**0.7088)&
          *(airmass**0.9108) )
      T_AA = 1.0 - 0.1*(1.0 - airmass + airmass**1.06)*&
          (1.0 - T_aerosol)
      T_AS = T_aerosol/T_AA

      Ruuclr = 0.0685 + (1.0-g)*(1.0-T_AS)

      F0 = 1.0
      Fddclr = F0*(cos(Z*degrad))*0.9662*T_rayleigh*T_o*T_gases*T_water*T_aerosol

      Fdif = F0*( cos(Z*degrad) )*0.79*T_o*T_gases*T_water*T_AA*&
          (0.5*(1.0-T_rayleigh) + g*(1.0-T_AS) )&
         /(1.0 - airmass + airmass**1.02 )

      Ftotal = ( Fddclr + Fdif )/(1.0 - albdo*Ruuclr)

      Tddclr = Fddclr/( F0*cos(Z*degrad) )  
      Tduclr = (Ftotal - Fddclr)/(F0*cos(Z*degrad))

    end subroutine BIRD


    subroutine WATERMODEL( tau, De, solarangle, Tducld, Ruucld )
        
      implicit none

      real, intent (in) :: tau, De, solarangle
      real, intent (out) :: Tducld, Ruucld

        ! Local vars
      real :: solarconst, Ptau, PDHI, delta, y, PPDHI


      solarconst = 1385.72180
      Ptau = (2.8850+0.002*(De-60.0))*solarangle-0.007347
      Ptau = max(Ptau, PTAU_MIN)
      PDHI = (1.0+(De-60.0)*0.0002)*1087.24*solarangle**1.1605

      delta = -0.644531*solarangle+1.20117+0.129807/solarangle &
          -0.00121096/(solarangle*solarangle) + &
          1.52587e-07/(solarangle*solarangle*solarangle)
      delta = Max (delta, DELTA_MIN)

      y = 0.012*(tau-Ptau)*solarangle
      PPDHI = (1.0+SINH(y))*PDHI*&
          exp(-( (log10(tau)-log10(Ptau))**2.0 )/delta)

      Tducld = PPDHI/(solarconst*solarangle)

      if (tau .lt. 1.0) then 
        Ruucld = 0.107359*tau
      endif

      if (tau .ge. 1.0) then 
        Ruucld = 1.03 - exp(-(0.5+log10(tau))*(0.5+log10(tau))/3.105 )
      endif

    end subroutine WATERMODEL


    subroutine ICEMODEL( tau, De, solarangle, Tducld, Ruucld )

      implicit none

      real, intent (in) :: tau, De, solarangle
      real, intent (out) :: Tducld, Ruucld

        ! Local vars

      real :: solarconst, Ptau, PDHI, delta, y, PPDHI

      solarconst = 1385.72180

      if ( De .le. 26.0 ) then 
        Ptau = 2.8487*solarangle- 0.0029
      endif
      if ( De .gt. 26.0 ) then 
        Ptau =  (2.8355 + (100.0-De)*0.006)*solarangle - 0.00612
        Ptau = max(Ptau, PTAU_MIN)
      endif

      PDHI = 1047.6367*solarangle**1.0883

      delta = -0.0549531*solarangle+0.617632+0.17876/(solarangle) &
          -0.002174/(solarangle*solarangle) 
      delta = Max (delta, DELTA_MIN)


      y = 0.01*(tau-Ptau)*solarangle
      PPDHI = (1.0+SINH(y))*PDHI*&
          exp(-( (log10(tau)-log10(Ptau))**2.0 )/delta)

      Tducld = PPDHI/(solarconst*solarangle)

      if (tau .lt. 1.0) then 
        Ruucld = 0.094039*tau
      endif

      if (tau .ge. 1.0) then 
        Ruucld = 1.02 - exp(-(0.5+log10(tau))*(0.5+log10(tau))/3.25 )
      endif

    end subroutine ICEMODEL

  end module module_ra_farms



