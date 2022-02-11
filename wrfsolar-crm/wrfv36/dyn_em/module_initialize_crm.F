!IDEAL:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!-----------------------------------------------------------------------
!
! * module_initialize_crm.F is based on module_initialize_les.F.
!   Some additional settings are available, including: 
!
!       - lat/lon for Coriolis factors and radiation 
!       - control initial random perturbation
!       - initialize albedo, tsk (and Nccn)
!       - non-zero model surface height
!       - some other land surface parameters for radiation
!       - constant dz, streched dz (only for flat surface), specified z
!     
!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain
   USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_timing
   USE module_configure
   USE module_init_utilities
   USE module_soil_pre ! added
#ifdef DM_PARALLEL
   USE module_dm
#endif


CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the 
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid 
   !  Local data.
   INTEGER :: idum1, idum2

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include <actual_new_args.inc>
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include <dummy_new_args.inc>
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include <dummy_new_decl.inc>

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k

   ! Local data

   INTEGER, PARAMETER :: nl_max = 30000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in
   INTEGER :: nl_in


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: z_scale, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm
   REAL    :: pi

   character (len=256) :: mminlu2 ! landuse data

!  stuff from original initialization that has been dropped from the Registry 
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zw
   REAL    :: qvf1, qvf2, pd_surf
   INTEGER :: it
   real :: thtmp, ptmp, temp(3)

   LOGICAL :: moisture_init
   !LOGICAL :: stretch_grid, dry_sounding
   LOGICAL :: dry_sounding
   INTEGER :: stretch_grid ! changed

  !INTEGER :: xs , xe , ys , ye ! for mountain. no need
  !REAL :: mtn_ht

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor

!  For CRM
   integer :: seedsize,n
   integer,allocatable :: seed(:)

   real,allocatable :: randxy(:,:)
   real,allocatable :: randz(:)

   integer :: nl_grid_in
   real :: z_level  ! for stretched dz grid
   real, dimension(nl_max) :: zw_grid_in 

   

#ifdef DM_PARALLEL
#    include <data_calls.inc>
#endif



   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT



   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi

   ! * for position of bubble. no need
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

  
   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

   ! here we check to see if the boundary conditions are set properly
   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif


! * nl_set: namelist set. see Registry-generated inc/nl_config.inc
!    ex  SUBROUTINE nl_set_cen_lat ( id_id , cen_lat )
!        USE module_configure, ONLY : model_config_rec 
!        model_config_rec%cen_lat(id_id) = cen_lat
 

    ! * for surface model
    !mminlu2 = '    '           ! not just empty. 4digits. work in phy_init in wrf.exe
    !mminlu2(1:4) = 'USGS'
    !CALL nl_set_mminlu(1, mminlu2) ! land use dataset

    CALL nl_set_mminlu(1, '    ')
    CALL nl_set_iswater(1,16)      ! land use index of water
    CALL nl_set_isice(1,3)         ! land use index of ice

    !CALL nl_set_iswater(1,0)

    ! * for map proj. 
    CALL nl_set_cen_lat(1,30.)
    CALL nl_set_cen_lon(1,-140.)

    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.) ! center latitude of the most coarse grid"      "degrees, negative is south"
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
    ! * MAP: "0=none (Cylindrical), 1=Lambert, 2=polar, 3=Mercator"

    CALL nl_get_iswater(1,grid%iswater)

!  here we initialize data that is not initialized 
!  in the input data

    DO j = jts, jte
      DO i = its, ite

         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.

         ! * map scale factor on mass, u and v grid *******************
         !grid%msft(i,j)     = 1.
         !grid%msfu(i,j)     = 1.
         !grid%msfv(i,j)     = 1.
 
         ! * latlon for radiation and Coriolis force *****************
         grid%xlat(i,j)     = config_flags%crm_lat
         grid%xlong(i,j)    = config_flags%crm_lon


         !grid%e(i,j)        = 2.0*EOMEG*cos(config_flags%crm_lat*DEGRAD) 
         !grid%f(i,j)        = 2.0*EOMEG*sin(config_flags%crm_lat*DEGRAD)
                                                    ! DEGRAD = pi/180
         grid%e(i,j)        = 0. 
         !grid%f(i,j)        = 8.5e-5
         grid%f(i,j)        = 0.
                                                   
         ! * for surface model *********************    
         grid%xland(i,j)     = 1.  ! land mask (1 for land, 2 for water), real
         grid%lu_index(i,j)  = 18  ! land use category

         grid%albedo(i,j)   = config_flags%crm_albedo
         grid%tsk(i,j)      = config_flags%crm_tsk

	 grid%tmn(i,j) = grid%tsk(i,j) - 5  ! soil temperature at lower boundary

      END DO
    END DO

    write(6,*) ' Coriolis e is ', grid%e(2,2), '  f is ', grid%f(2,2)
    write(6,*) ' xlat is ', grid%xlat(2,2), '  xlong is ', grid%xlong(2,2)
    write(6,*) ' Initial albedo is ', grid%albedo(2,2)

 
   ! * for surface model 

    other_masked_fields : SELECT CASE ( model_config_rec%sf_surface_physics(grid%id) )

    CASE (SLABSCHEME) ! thermal diffusion scheme: 1

       ! do nothing

    CASE (LSMSCHEME)  ! Unified Noah land-surface model: 2

        DO j = jts , MIN(jde-1,jte)
           DO i = its , MIN(ide-1,ite)
              IF (grid%xland(i,j) .lt. 1.5) THEN
                 grid%vegfra(i,j) = 0.5
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 18
                 grid%isltyp(i,j) = 8
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ELSE
                 grid%vegfra(i,j) = 0.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 16
                 grid%isltyp(i,j) = 14
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ENDIF
           END DO
        END DO

     CASE (RUCLSMSCHEME) ! RUC land-surface model

        ! do nothing

   END SELECT other_masked_fields


   DO j = jts, jte
   DO k = kts, kte
   DO i = its, ite
   
      grid%ww(i,k,j)     = 0. ! "ww" "mu-coupled eta-dot" "Pa s-1"

   END DO
   END DO
   END DO

   grid%step_number = 0 ! * current number of timestep

! Process the soil; note that there are some things hard-wired into share/module_soil_pre.F
   CALL process_soil_ideal(grid%xland,grid%xice,grid%vegfra,grid%snow,grid%canwat, &
                                     grid%ivgtyp,grid%isltyp,grid%tslb,grid%smois, &
              grid%tsk,grid%tmn,grid%zs,grid%dzs,model_config_rec%num_soil_layers, &
                model_config_rec%sf_surface_physics(grid%id),&
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )


!---------------------------------------------------------------------
! Set up the grid.
!
! modified to read the dry sounding first and set the grid by using z
!
! config_flags%crm_stretch =
! 1 : constant grid in eta (em_les)
! 2 : stretched grid in eta (em_les)
! 3 : constant grid in z
! 4 : read z profile from external file 
!

   stretch_grid = config_flags%crm_stretch 

! inverse horizontal grids

   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy 

!  get the sounding from the ascii sounding file, first get dry sounding and 
!  calculate base state

  dry_sounding = .true.
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '
  ENDIF

  ! use config_flags%crm_zsfc for calculation of the lowest pm

  CALL get_sounding(config_flags%crm_zsfc, zk, p_in, pd_in, theta, rho, &
                    u, v, qv, dry_sounding, nl_max, nl_in )
!  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 ) 

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

  ! Surface height 
  DO j=jts,jte
  DO i=its,ite
    grid%ht(i,j) = config_flags%crm_zsfc ! "HGT" "Terrain Height"   "m"
  ENDDO
  ENDDO
  
  write(6,*) 'surface height is ',grid%ht(1,1), ' m'

  ! "phb" "base-state geopotential"  "m2 s-2"
  ! "ph0" "initial geopotential"     "m2 s-2"
  ! "mub" "base state dry air mass in column" "Pa"
  ! "pb"  "BASE STATE PRESSURE " "Pa"

  DO j=jts,jte
  DO i=its,ite
    grid%phb(i,1,j) = g * grid%ht(i,j)   ! ht to phb only at surface
    grid%ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO



  ! ***** set vertical grid *****

  if (stretch_grid .eq. 0 .or. stretch_grid .eq. 1) then

    !  find ptop for the desired ztop (ztop is input from the namelist),
    !  and set surface level
    grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

    ! ***** Option 1: set znw first. then, calc phb *****   

    IF (stretch_grid .eq. 0) THEN ! uniform eta grid spacing

      DO k=1, kde
        grid%znw(k) = 1. - float(k-1)/float(kde-1)
      ENDDO

      write(6,*) ' uniform eta grids',stretch_grid

    ELSE  ! exponential stretch for eta (nearly constant dz) 

      z_scale = .40          ! * for stretch.

      DO k=1, kde
        grid%znw(k) = (exp(-float(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
      ENDDO

      write(6,*) ' stretched eta grids',stretch_grid

    ENDIF

   ! "znw"  "eta values on full (w) levels"
   ! "znu"  "eta values on half (mass) levels"
   ! "dnw"  "d(eta) values between full (w) levels"
   ! "dn "  "d(eta) values between half (mass) levels"
   ! "fnp"  "lower weight for vertical stretching"
   ! "fnm"  "upper weight for vertical stretching" 
   ! "cf1"  "2nd order extrapolation constant"
   ! "cfn"  "extrapolation constant"

   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k) 
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO

   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)

  ! calc phb 

   DO J = jts, jte
   DO I = its, ite

     p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
     grid%mub(i,j) = p_surf-grid%p_top

!  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state

     DO K = 1, kte-1
       p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
       grid%pb(i,k,j) = p_level
       grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0
       grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
     ENDDO

!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.

     DO k  = 2,kte
       grid%phb(i,k,j) = grid%phb(i,k-1,j) - grid%dnw(k-1)*grid%mub(i,j)*grid%alb(i,k-1,j)
     ENDDO

    ENDDO
    ENDDO



  elseif (stretch_grid .eq. 2 .or. stretch_grid .eq. 3) then

  ! ***** Option 2: set phb first. then, calc znw etc *****

    IF (stretch_grid .eq. 2) THEN ! uniform z grid spacing (base state geopotential height)

      !  find ptop for the desired ztop (ztop is input from the namelist),
      !  and set surface level
      grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

      DO k=2, kde ! phb(i,1,j) has been set

        z_level = grid%ht(its,jts) + (config_flags%ztop - grid%ht(its,jts)) * float(k-1)/float(kde-1)

        DO J = jts, jte
        DO I = its, ite
          grid%phb(i,k,j) = z_level * g
        ENDDO
        ENDDO

      ENDDO

      write(6,*) ' uniform phb grids',stretch_grid

    ELSE  ! read external ascii z file

      call get_zw_grid(zw_grid_in, nl_grid_in, nl_max)

      call wrf_dm_bcast_real( zw_grid_in , nl_max )
      call wrf_dm_bcast_integer ( nl_grid_in , 1 ) 

      write(6,*) ' read zw(phb) grids',stretch_grid

      if (zw_grid_in(1) .lt. config_flags%crm_zsfc) then 
        write(6,*) ' ERROR: zw_grid_in(1), crm_zsfc ', zw_grid_in(1), config_flags%crm_zsfc
        call wrf_error_fatal ( ' first layer < crm_zsfc in namelist' )
      endif
      if (nl_grid_in .ne. kte) then
        write(6,*) ' ERROR: nl_grid_in, kte ', nl_grid_in, kte
        call wrf_error_fatal ( ' Num of lines in grid file != Num of vertical levels in namelist' )
      endif
	  
      if (config_flags%ztop .ne. zw_grid_in(kde)) then 
        write(6,*) 'model top in namelist (ztop):', config_flags%ztop
        write(6,*) 'zw_grid_in(kde  )', zw_grid_in(kde)
        write(6,*) 'zw_grid_in(kde-1)', zw_grid_in(kde-1)
        write(6,*) 'zw_grid_in(kde) will be used for model top'
      endif

      !  find ptop from zw_grid_in(kde) instead of ztop in namelist 
      grid%p_top = interp_0( p_in, zk, zw_grid_in(kde), nl_in )


  
      DO k=2, kde ! phb(i,1,j) has been set

        DO J = jts, jte
        DO I = its, ite
          grid%phb(i,k,j) = zw_grid_in(k) * g
        ENDDO
        ENDDO

      ENDDO

    ENDIF ! stretch_grid 2 or 3

	do k=1,kde
	    write(6,*) 'k phb/g',k, grid%phb(its,k,jts)/g
	end do 

  ! ***** calc znw etc *****

   DO J = jts, jte
   DO I = its, ite

     p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
     grid%mub(i,j) = p_surf - grid%p_top

     grid%pb(i,1,j) = p_surf 
     grid%t_init(i,1,j) = interp_0( theta, p_in, p_surf, nl_in ) - t0
     grid%alb(i,1,j) = (r_d/p1000mb)*(grid%t_init(i,1,j)+t0)*(grid%pb(i,1,j)/p1000mb)**cvpm

     ! based on dry hydrostatic sounding, for given phb 
     ! interp theta (from interp) and compute 1/rho from eqn. of state

     DO K = 2, kte-1

       grid%pb(i,k,j) = grid%pb(i,k-1,j) - (grid%phb(i,k,j) - grid%phb(i,k-1,j)) / grid%alb(i,k-1,j)
       grid%t_init(i,k,j) = interp_0( theta, p_in, grid%pb(i,k,j), nl_in ) - t0
       grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm

     ENDDO

    ENDDO
    ENDDO

     DO K = 1, kte-1 

       grid%znw(k) = (grid%pb(its,k,jts) - grid%p_top) / grid%mub(its,jts)

     ENDDO


   ! "znw"  "eta values on full (w) levels"
   ! "znu"  "eta values on half (mass) levels"
   ! "dnw"  "d(eta) values between full (w) levels"
   ! "dn "  "d(eta) values between half (mass) levels"
   ! "fnp"  "lower weight for vertical stretching"
   ! "fnm"  "upper weight for vertical stretching" 
   ! "cf1"  "2nd order extrapolation constant"
   ! "cfn"  "extrapolation constant"

    DO k=1, kde-1
      grid%dnw(k) = grid%znw(k+1) - grid%znw(k) 
      grid%rdnw(k) = 1./grid%dnw(k)
      grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
    ENDDO

    DO k=2, kde-1
      grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
      grid%rdn(k) = 1./grid%dn(k)
      grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
      grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
    ENDDO

    cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
    cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
    grid%cf1  = grid%fnp(2) + cof1
    grid%cf2  = grid%fnm(2) - cof1 - cof2
    grid%cf3  = cof2       

    grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
    grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)

  else

    write(6,*) ' ERROR: crm_stretch = ', stretch_grid
    call wrf_error_fatal ( ' crm_stretch has to be in 0-3 ' )

  endif ! stretch_grid 


  IF ( wrf_dm_on_monitor() ) THEN
    write(6,*) ' ptop is ',grid%p_top
    write(6,*) ' base state grid%mub(1,1), p_surf is ',grid%mub(1,1),grid%mub(1,1)+grid%p_top
  ENDIF


! ***** full state is being used as perturbation *****

!  calculate full state for each column - this includes moisture.

  write(6,*) ' getting moist sounding for full state '
  dry_sounding = .false.
  CALL get_sounding(config_flags%crm_zsfc, zk, p_in, pd_in, theta, rho, &
                       u, v, qv, dry_sounding, nl_max, nl_in )

  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 ) 

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)

!  At this point grid%p_top is already set. find the DRY mass in the column 
!  by interpolating the DRY pressure.  

   pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )

!  compute the perturbation mass and the full mass

    grid%mu_1(i,j) = pd_surf-grid%p_top - grid%mub(i,j)
    grid%mu_2(i,j) = grid%mu_1(i,j)
    grid%mu0(i,j) = grid%mu_1(i,j) + grid%mub(i,j)

! given the dry pressure and coordinate system, interp the potential
! temperature and qv

    do k=1,kde-1

      p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top

      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure

    k = kte-1  ! top level

    qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2

!    grid%p(i,k,j) = - 0.5*grid%mu_1(i,j)/grid%rdnw(k)
    grid%p(i,k,j) = - 0.5*(grid%mu_1(i,j)+qvf1*grid%mub(i,j))/grid%rdnw(k)/qvf2
    qvf = 1. + rvovrd*moist(i,k,j,P_QV)
    grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
    grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)

!  down the column

    do k=kte-2,1,-1
      qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%p(i,k,j) = grid%p(i,k+1,j) - (grid%mu_1(i,j) + qvf1*grid%mub(i,j))/qvf2/grid%rdn(k+1)
      qvf = 1. + rvovrd*moist(i,k,j,P_QV)
      grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                  (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
      grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
    enddo

!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.


    grid%ph_1(i,1,j) = 0.
    DO k  = 2,kte
      grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (1./grid%rdnw(k-1))*(       &
                   (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
                    grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
      grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
      grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
    ENDDO

    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     write(6,*) ' grid%ph_1 calc ',grid%ph_1(2,1,2),grid%ph_1(2,2,2),&
                              grid%mu_1(2,2)+grid%mub(2,2),grid%mu_1(2,2), &
                              grid%alb(2,1,2),grid%al(1,2,1),grid%rdnw(1)
    endif
    ENDIF

  ENDDO
  ENDDO



!  Thermal perturbation to kick off turbulence. 
!  For LES, change the initial random perturbations
!  put in preturbation theta and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

  ! Update random seeds. Generate random numbers in X-Y space and 
  ! then use them as the seeds to generate the perturbations in 
  ! each column. 

  call random_seed(size=seedsize)
  allocate(seed(seedsize))

  seed = config_flags%crm_pert_random_seed
  call random_seed(put=seed)

  allocate(randxy(ids:ide-1,jds:jde-1))
  allocate(randz(1:config_flags%crm_num_pert_layers))

  call random_number(randxy)

  DO j = jts, min(jde-1,jte)

    DO i = its, min(ide-1,ite)

      seed = int(randxy(i,j)*1E9)
      call random_seed(put=seed)
      call random_number(randz)
      randz = (randz - 0.5) * 2  
      ! random perturbation, -1 < randz < 1 
      
      DO k = 1, config_flags%crm_num_pert_layers

           grid%t_1(i,k,j) = grid%t_1(i,k,j) + config_flags%crm_pert_amp * randz(k)
           grid%t_2(i,k,j) = grid%t_1(i,k,j)
           qvf = 1. + rvovrd * moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb) * (grid%t_1(i,k,j)+t0) * qvf * &
                        (((grid%p(i,k,j) + grid%pb(i,k,j)) / p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
           ! * alb: inverse base density m3
           
      ENDDO

      ! rebalance hydrostatically
      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (1./grid%rdnw(k-1))*(       &
                     (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
                      grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
        grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO

  write(6,*) ' perturbation amplitude pert_amp is ', config_flags%crm_pert_amp
  write(6,*) ' perturbation range 1 to ', config_flags%crm_num_pert_layers




   IF ( wrf_dm_on_monitor() ) THEN
   write(6,*) ' grid%mu_1 from comp ', grid%mu_1(1,1)
   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i5,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1)+grid%phb(1,k,1), &
                                      grid%p(1,k,1)+grid%pb(1,k,1), grid%alt(1,k,1), &
                                      grid%t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, alp, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i5,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1), &
                                      grid%p(1,k,1), grid%al(1,k,1), &
                                      grid%t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF

! interp v

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%phb(i,1,j)+grid%phb(i,1,j-1))/g
    END IF
    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

! interp u

  DO J = jts, min(jde-1,jte)
  DO I = its, ite

    IF (i == ids) THEN
      z_at_u = grid%phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%phb(i,1,j)+grid%phb(i-1,1,j))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_u, nl_in )

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO


! ***** for wdm, initial ccn ****** 

  micro_select: SELECT CASE(config_flags%mp_physics)

   CASE (WDM5SCHEME, WDM6SCHEME)

      DO J = jts, min(jde-1,jte)
      DO K = kts, kte
      DO I = its, min(ide-1,ite)

        scalar(i,k,j,P_QNN) = config_flags%crm_init_ccn * 1000000 / grid%alb(i,k,j)  ! kg-1

      ENDDO
      ENDDO
      ENDDO

   END SELECT micro_select



!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)

    grid%h_diabatic(i,k,j) = 0.
    ! * PREVIOUS TIMESTEP CONDENSATIONAL HEATING K

  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
    DO k=1,kte-1
      grid%t_base(k) = grid%t_1(1,k,1)
      grid%qv_base(k) = moist(1,k,1,P_QV)
      grid%u_base(k) = grid%u_1(1,k,1)
      grid%v_base(k) = grid%v_1(1,k,1)
      grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g
    ENDDO
  ENDIF

  CALL wrf_dm_bcast_real( grid%t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base, kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)

     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!     grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
!     grid%tmn(I,J)=grid%tsk(I,J)-0.5

  ENDDO
  ENDDO

 END SUBROUTINE init_domain_rk

!---------------------------------------------------------------------

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( zsfc, zk, p, p_dry, theta, rho, &
                               u, v, qv, dry, nl_max, nl_in )
      implicit none

      integer nl_max, nl_in

      real zsfc ! added
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n
      parameter(n=30000)
      logical debug
      parameter( debug = .true.)

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      real r
      parameter (r = r_d)
      integer k, it, nl
      real qvf, qvf1, dz

!  first, read the sounding

      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1) 
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)

      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1) 
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1) - zsfc

          if(dz < 0.0) then

             write(6,*) "Your first input sounding level is below the surface level"
             stop "module_initialize_crm: get_sounding"
             
          endif

          do it=1,10
            pm_input(1) = p_surf &
                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
 
            do it=1,10
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
            enddo
          enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i5,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
      logical end_of_file
      logical debug

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        k = k+1
        if(debug) write(6,'(1x,i5,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

!     subroutine get_zw_grid( zw_grid_in, nl_grid_in, nl_max )
      subroutine get_zw_grid( zw, nl, nl_max )

!---------------------------------------------------------------------------
!
! This subroutine reads vertical grid levels in geopotential height
!  from external ascii file.  07/27/2011 sendo
! 
!---------------------------------------------------------------------------

      implicit none

      integer nl_max 
      !parameter (nl_max = 1000)
      integer nl

      real zw(nl_max)

      logical end_of_file
      logical debug

      integer k
      parameter( debug = .true.)

      open(unit=10,file='input_zw_grid',form='formatted',status='old')
      rewind(10)

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) zw(k+1)
        k = k+1
!        if(debug) write(6,'(1x,i3,1x,f8.2)') k, zw(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine get_zw_grid


END MODULE module_initialize_ideal


