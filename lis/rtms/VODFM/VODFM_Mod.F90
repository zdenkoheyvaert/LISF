!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module VODFM_Mod
    !BOP
    !
    ! !MODULE: VODFM_Mod
    !
    ! !DESCRIPTION:
    !    This modules provides implements a linear forward model and a
    !    SVR-based forward model for VOD.
    !    Although it's not a real RTM it programmatically works similarly and is
    !    therefore implemented analogously to a RTM.
    !
    !    The module provides the following options for lis.config:
    !
    !    VODFM model type: 
    !        "linear" or "SVR"
    !    VODFM parameter file:
    !        Pathname to the parameter file
    !
    !    The linear model calculates VOD as linear function of LAI and SM_i,
    !    i=1...4. To use it, you need to provide a netCDF file with 6
    !    variables, "intercept", "LAI_coef", and "SM{i}_coef" (i = 1...4).
    !    Each variable must have the dimensions (ntimes, ngrid), where ntimes
    !    is either 1 or 12, depending on whether a single model for the full
    !    year should be used or monthly varying models.
    !    ngrid is the dimension of the model tile space (i.e. a flat 1D array
    !    with only land points such that east_west is changing fastest)
    !
    !    The SVR based model calculates VOD via:
    !
    !        VOD[lat,lon] = \sum_{j}^{n_SV_[lat,lon]} dual_coef_[lat,lon,j]
    !                       * K(x[lat,lon], support_vectors_[lat,lon,j])
    !
    !    with x = (LAI, SM1, .., SM4)^T
    !    and K(x, s) = exp(-\sum_k gamma_[k] (x[k] - s[k])**2)
    !
    !    To use the SVR model, you need to provide a netCDF file with 5 variables,
    !    intercept_, n_SV_, dual_coef_, support_vectors_, gamma_.
    !    The support vectors are 5-d vectors with dimensions (LAI, SM1, SM2, SM3,
    !    SM4).
    ! 
    !    The variables on file have to have the following specifications:
    !    - intercept_: float, (ngrid,)
    !    - n_SV_: float, (ngrid,)
    !    - dual_coef_: float, (ngrid, n_SV) (n_SV == maximum number of support vectors)
    !    - gamma_: float, (ngrid, 5)
    !    - support_vectors_: float, (ngrid, n_SV, 5)
    !
    ! !HISTORY:
    ! 02 Mar 2022: Samuel Scherrer; initial contribution based on WCMRTM
    !
    ! !USES:

#if (defined RTMS)

    use ESMF
    use LIS_coreMod
    use LIS_RTMMod
    use LIS_logMod

    implicit none

    PRIVATE

    !-----------------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    !-----------------------------------------------------------------------------
    public :: VODFM_initialize
    public :: VODFM_f2t
    public :: VODFM_run
    public :: VODFM_output
    public :: VODFM_geometry
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: vodfm_struc
    !EOP

    ! The actual implementation of the model equations is done in the
    ! subclasses below, this is just to provide a common interface.
    ! Using instances of this type will raise an error.
    type, public ::  vodfm_type_dec
        character*256 :: parameter_fname
        !-------output------------!
        real, allocatable :: VOD(:)
        contains
        procedure, pass(self) :: initialize => VODFM_initialize_default
        procedure, pass(self) :: run => VODFM_run_default
    end type vodfm_type_dec

    type, extends(vodfm_type_dec) :: vodfm_linear_model
        integer           :: ntimes
        integer           :: timeidx
        real, allocatable :: intercept(:,:)
        real, allocatable :: laicoef(:,:)
        real, allocatable :: sm1coef(:,:)
        real, allocatable :: sm2coef(:,:)
        real, allocatable :: sm3coef(:,:)
        real, allocatable :: sm4coef(:,:)
        contains
            procedure, pass(self) :: initialize => VODFM_initialize_linear_model
            procedure, pass(self) :: run => VODFM_run_linear_model
    end type vodfm_linear_model

    type, extends(vodfm_type_dec) :: vodfm_svr_model
        integer           :: n_SV
        real, allocatable :: intercept(:)
        real, allocatable :: actual_n_SV(:)
        ! dual coef has n_SV as additional dim
        real, allocatable :: dual_coef(:,:)
        ! for each pixel, this is a matrix n_SV x 5
        real, allocatable :: support_vectors(:,:,:)
        ! for each pixel, this is a vector of length 5
        real, allocatable :: gam(:,:)
        contains
            procedure, pass(self) :: initialize => VODFM_initialize_svr_model
            procedure, pass(self) :: run => VODFM_run_svr_model
    end type vodfm_svr_model

    type, extends(vodfm_type_dec) :: vodfm_X_model
        real, allocatable :: intercept(:)
        real, allocatable :: cwccoef(:)
        real, allocatable :: laicoef(:)
        real, allocatable :: laipsicoef(:)
        real, allocatable :: laisqtvegcoef(:)
        contains
            procedure, pass(self) :: initialize => VODFM_initialize_X_model
            procedure, pass(self) :: run => VODFM_run_X_model
    end type vodfm_X_model


    class(vodfm_type_dec), allocatable :: vodfm_struc(:)

    SAVE

contains
    !BOP
    !
    ! !ROUTINE: VODFM_initialize
    ! \label{VODFM_initialize}
    !
    ! !INTERFACE:
    subroutine VODFM_initialize()
        ! !DESCRIPTION:
        !
        !  This routine creates the datatypes and allocates memory for noahMP3.6-specific
        !  variables. It also invokes the routine to read the runtime specific options
        !  for noahMP3.6 from the configuration file.
        !
        !  The routines invoked are:
        !  \begin{description}
        !   \item[readVODFMcrd](\ref{readVODFMcrd}) \newline
        !
        !EOP
        implicit none

        integer :: rc, ios
        integer :: n, nid, ngrid, ngridId
        character*100 :: modeltype(LIS_rc%nnest)
        character*256 :: parameter_fname(LIS_rc%nnest)

        write(LIS_logunit,*) "[INFO] Starting VODFM setup"


        ! read config from file
        call ESMF_ConfigFindLabel(LIS_config, "VODFM model type:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, modeltype(n), rc=rc)
            call LIS_verify(rc, "VODFM model type: not defined")
            if (modeltype(n) .ne. "linear"&
                 .and. modeltype(n) .ne. "SVR"&
                 .and. modeltype(n) .ne. "X") then
                write(LIS_logunit, *)&
                     "[ERR] VODFM model type must be 'linear', 'SVR' or 'X'"
                call LIS_endrun
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config, "VODFM parameter file:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, parameter_fname(n), rc=rc)
            call LIS_verify(rc, "VODFM parameter file: not defined")
        enddo

        if (modeltype(1) == "linear") then
            allocate(vodfm_linear_model :: vodfm_struc(LIS_rc%nnest))
        elseif (modeltype(1) == "SVR") then
            allocate(vodfm_svr_model :: vodfm_struc(LIS_rc%nnest))
        elseif (modeltype(1) == "X") then
            allocate(vodfm_X_model :: vodfm_struc(LIS_rc%nnest))
        else
            write(LIS_logunit, *)&
                 "[ERR] VODFM model type must be 'linear', 'SVR' or 'X'"
            call LIS_endrun
        endif

        do n=1,LIS_rc%nnest
            vodfm_struc(n)%parameter_fname = parameter_fname(n)
            allocate(vodfm_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 1")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 2")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 3")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 4")
            call add_sfc_fields(n,LIS_sfcState(n), "Canopy Water Content")
            call add_sfc_fields(n,LIS_sfcState(n), "Vegetation Transpiration")
            call add_sfc_fields(n,LIS_sfcState(n), "Root Zone Soil Water Potential")

            call add_sfc_fields(n,LIS_forwardState(n),"VODFM_VOD")
        enddo


        ! read parameters/initialize model
        do n=1, LIS_rc%nnest
            call vodfm_struc(n)%initialize(n)
        enddo

        write(LIS_logunit,*) '[INFO] Finished VODFM setup'
    end subroutine VODFM_initialize

    subroutine VODFM_initialize_default(self, n)
        class(vodfm_type_dec), intent(inout) :: self
        integer, intent(in) :: n
        write(LIS_logunit,*) "[ERR] VODFM should use linear or svr model"
        call LIS_endrun
    end subroutine VODFM_initialize_default

    subroutine VODFM_initialize_linear_model(self, n)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vodfm_linear_model), intent(inout) :: self
        integer, intent(in) :: n

        integer :: ios, nid, npatch
        integer :: ngridId, ntimesId, ngrid, ntimes

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        ! try opening the parameter file
        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(self%parameter_fname)
        ios = nf90_open(path=trim(self%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(vodfm_struc(n)%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(self%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! read ntimes
        ios = nf90_inq_dimid(nid, "ntimes", ntimesId)
        call LIS_verify(ios, "Error nf90_inq_varid: ntimes")
        ios = nf90_inquire_dimension(nid, ntimesId, len=ntimes)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ntimes")
        self%ntimes = ntimes

        ! now that we have ntimes, we can allocate the coefficient arrays
        ! for the nest
        npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
        allocate(self%intercept(npatch, ntimes))
        allocate(self%laicoef(npatch, ntimes))
        allocate(self%sm1coef(npatch, ntimes))
        allocate(self%sm2coef(npatch, ntimes))
        allocate(self%sm3coef(npatch, ntimes))
        allocate(self%sm4coef(npatch, ntimes))

        call read_2d_coef_from_file(n, nid, ngrid, ntimes, "intercept", &
            self%intercept)
        call read_2d_coef_from_file(n, nid, ngrid, ntimes, "LAI_coef", &
            self%laicoef)
        call read_2d_coef_from_file(n, nid, ngrid, ntimes, "SM1_coef", &
            self%sm1coef)
        call read_2d_coef_from_file(n, nid, ngrid, ntimes, "SM2_coef", &
            self%sm2coef)
        call read_2d_coef_from_file(n, nid, ngrid, ntimes, "SM3_coef", &
            self%sm3coef)
        call read_2d_coef_from_file(n, nid, ngrid, ntimes, "SM4_coef", &
            self%sm4coef)
#endif
    end subroutine VODFM_initialize_linear_model

    subroutine VODFM_initialize_svr_model(self, n)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vodfm_svr_model), intent(inout) :: self
        integer, intent(in) :: n

        integer :: ios, nid, npatch
        integer :: ngridId, ngrid, n_SV, n_SVId

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        ! try opening the parameter file
        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(self%parameter_fname)
        ios = nf90_open(path=trim(self%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(vodfm_struc(n)%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(self%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! read n_SV, because we don't know it's length a priori
        ios = nf90_inq_dimid(nid, "n_SV", n_SVId)
        call LIS_verify(ios, "Error nf90_inq_varid: n_SV")
        ios = nf90_inquire_dimension(nid, n_SVId, len=n_SV)
        call LIS_verify(ios, "Error nf90_inquire_dimension: n_SV")
        self%n_SV = n_SV

        ! now that we have n_SV, we can allocate the coefficient arrays
        ! for the nest
        npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
        allocate(self%intercept(npatch))
        allocate(self%actual_n_SV(npatch))
        allocate(self%dual_coef(n_SV, npatch))
        allocate(self%support_vectors(5, n_SV, npatch))
        allocate(self%gam(5, npatch))

        ! read the arrays from the file
        call read_1d_coef_from_file(n, nid, ngrid, "intercept_", self%intercept)
        call read_1d_coef_from_file(n, nid, ngrid, "n_SV_", self%actual_n_SV)
        call read_2d_coef_from_file(n, nid, n_SV, ngrid, "dual_coef_", &
            self%dual_coef, ngrid_first=.false.)
        call read_2d_coef_from_file(n, nid, 5, ngrid, "gamma_", self%gam, &
            ngrid_first=.false.)
        call read_3d_coef_from_file(n, nid, 5, n_SV, ngrid, "support_vectors_",&
             self%support_vectors, ngrid_first=.false.)
#endif
    end subroutine VODFM_initialize_svr_model

    subroutine VODFM_initialize_X_model(self, n)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vodfm_X_model), intent(inout) :: self
        integer, intent(in) :: n

        integer :: ios, nid, npatch
        integer :: ngridId, ngrid

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        ! try opening the parameter file
        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(self%parameter_fname)
        ios = nf90_open(path=trim(self%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(vodfm_struc(n)%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(self%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! now that we have ntimes, we can allocate the coefficient arrays
        ! for the nest
        npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
        allocate(self%intercept(npatch))
        allocate(self%cwccoef(npatch))
        allocate(self%laicoef(npatch))
        allocate(self%laipsicoef(npatch))
        allocate(self%laisqtvegcoef(npatch))

        call read_1d_coef_from_file(n, nid, ngrid, "intercept", self%intercept)
        call read_1d_coef_from_file(n, nid, ngrid, "CWC_coef", self%cwccoef)
        call read_1d_coef_from_file(n, nid, ngrid, "LAI_coef", self%laicoef)
        call read_1d_coef_from_file(n, nid, ngrid, "LAIPSI_coef", self%laipsicoef)
        call read_1d_coef_from_file(n, nid, ngrid, "LAI2TVEG_coef", self%laisqtvegcoef)
#endif
    end subroutine VODFM_initialize_X_model

    subroutine read_1d_coef_from_file(n, nid, ngrid, varname, coef)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, ngrid
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:)

        real, allocatable :: coef_file(:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else

        allocate(coef_file(ngrid))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        call LIS_convertVarToLocalSpace(n, coef_file, coef_grid)
        call gridvar_to_patchvar(&
             n, LIS_rc%lsm_index, coef_grid, coef)
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_1d_coef_from_file

    subroutine read_2d_coef_from_file(n, nid, n1, n2, varname, coef, ngrid_first)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, n1, n2
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:)
        logical, intent(in), optional :: ngrid_first

        logical :: ngrid_is_first
        real, allocatable :: coef_file(:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        if (present(ngrid_first)) then
            ngrid_is_first = ngrid_first
        else
            ngrid_is_first = .true.
        endif

        allocate(coef_file(n1, n2))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))

        allocate(coef_grid(LIS_rc%ngrid(n)))
        if (ngrid_is_first) then
            do j=1,n2
                call LIS_convertVarToLocalSpace(n, coef_file(:,j), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(:, j))
            enddo
        else ! ngrid is the second dimension, we have to loop over n1
            do j=1,n1
                call LIS_convertVarToLocalSpace(n, coef_file(j, :), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(j, :))
            enddo
        endif

        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_2d_coef_from_file

    subroutine read_3d_coef_from_file(n, nid, n1, n2, n3, varname, coef, &
            ngrid_first)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, n1, n2, n3
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:,:)
        logical, intent(in), optional :: ngrid_first

        logical :: ngrid_is_first
        real, allocatable :: coef_file(:,:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j2, j3

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        if (present(ngrid_first)) then
            ngrid_is_first = ngrid_first
        else
            ngrid_is_first = .true.
        endif

        allocate(coef_file(n1, n2, n3))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        if (ngrid_is_first) then
            do j2=1,n2
                do j3=1,n3
                    call LIS_convertVarToLocalSpace(n, coef_file(:,j2,j3), coef_grid)
                    call gridvar_to_patchvar(&
                         n, LIS_rc%lsm_index, coef_grid, coef(:, j2,j3))
                 enddo
            enddo
        else  ! ngrid is last
            do j2=1,n1
                do j3=1,n2
                    call LIS_convertVarToLocalSpace(n, coef_file(j2,j3, :), coef_grid)
                    call gridvar_to_patchvar(&
                         n, LIS_rc%lsm_index, coef_grid, coef(j2,j3, :))
                 enddo
            enddo
        endif
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_3d_coef_from_file

    subroutine gridvar_to_patchvar(n,m,gvar,tvar)
        ! Converts a variable in local gridspace (length LIS_rc%ngrid(n))
        ! to local patch space (length LIS_rc%npatch(n,m))
        ! patch space = ensembles * ngrid

        implicit none

        integer, intent(in) :: n 
        integer, intent(in) :: m
        real, intent(in)    :: gvar(LIS_rc%ngrid(n))
        real, intent(inout) :: tvar(LIS_rc%npatch(n,m))
        integer             :: t,r,c

        do t=1,LIS_rc%npatch(n,m)
            r = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            c = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            if (LIS_domain(n)%gindex(c,r).ge.0) then
                tvar(t) = gvar(LIS_domain(n)%gindex(c,r))
            else
                tvar(t) = LIS_rc%udef
            endif
        enddo

    end subroutine gridvar_to_patchvar
    !!--------------------------------------------------------------------------------



    subroutine add_sfc_fields(n, sfcState,varname)

        implicit none

        integer            :: n
        type(ESMF_State)   :: sfcState
        character(len=*)   :: varname

        type(ESMF_Field)     :: varField
        type(ESMF_ArraySpec) :: arrspec
        integer              :: status
        real :: sum
        call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        varField = ESMF_FieldCreate(arrayspec=arrSpec, &
             grid=LIS_vecTile(n), name=trim(varname), &
             rc=status)
        call LIS_verify(status, 'Error in field_create of '//trim(varname))

        call ESMF_StateAdd(sfcState, (/varField/), rc=status)
        call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

    end subroutine add_sfc_fields


    subroutine VODFM_f2t(n)

        implicit none

        integer, intent(in)    :: n

    end subroutine VODFM_f2t


    subroutine VODFM_geometry(n)
        implicit none
        integer, intent(in)    :: n

    end subroutine VODFM_geometry

    subroutine VODFM_run(n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        integer, intent(in) :: n

        call vodfm_struc(n)%run(n)

    end subroutine VODFM_run

    subroutine VODFM_run_default(self, n)
        class(vodfm_type_dec), intent(inout) :: self
        integer, intent(in) :: n
        write(LIS_logunit,*) "[ERR] VODFM should use linear or svr model"
        call LIS_endrun
    end subroutine VODFM_run_default

    subroutine VODFM_run_linear_model(self, n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodfm_linear_model), intent(inout) :: self
        integer, intent(in) :: n

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
        real, pointer       :: vodval(:)
        real                :: intercept
        real                :: laicoef, sm1coef, sm2coef, sm3coef, sm4coef
        logical             :: coefs_valid, values_valid

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 1", sm1)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 2", sm2)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 3", sm3)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 4", sm4)

        if (self%ntimes == 1) then
            self%timeidx = 1
        elseif (self%ntimes == 12) then
            self%timeidx = LIS_rc%mo
        else
            write(LIS_logunit, *) "[ERR] ntimes in "//trim(self%parameter_fname)&
                 //" must be 1 or 12, but is ", self%ntimes
            call LIS_endrun
        endif

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

            intercept = self%intercept(t, self%timeidx)
            laicoef = self%laicoef(t, self%timeidx)
            sm1coef = self%sm1coef(t, self%timeidx)
            sm2coef = self%sm2coef(t, self%timeidx)
            sm3coef = self%sm3coef(t, self%timeidx)
            sm4coef = self%sm4coef(t, self%timeidx)

            ! For some pixels no VOD was available and therefore no forward
            ! model was fitted. No assimilation will take place over these
            ! pixels anyways, so it's no problem to not predict anything here
            coefs_valid = (.not.isnan(intercept).and.intercept.ne.LIS_rc%udef&
                 .and..not.isnan(sm1coef).and.sm1coef.ne.LIS_rc%udef&
                 .and..not.isnan(sm2coef).and.sm2coef.ne.LIS_rc%udef&
                 .and..not.isnan(sm3coef).and.sm3coef.ne.LIS_rc%udef&
                 .and..not.isnan(sm4coef).and.sm4coef.ne.LIS_rc%udef)

            ! normally the modelled values should not be invalid, but just to
            ! be on the safe side
            values_valid = (.not.isnan(lai(t)).and.lai(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm1(t)).and.sm1(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm2(t)).and.sm2(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm3(t)).and.sm3(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm4(t)).and.sm4(t).ne.LIS_rc%udef)

            if (coefs_valid.and.values_valid) then
                self%VOD(t) = intercept&
                     + laicoef * lai(t)&
                     + sm1coef * sm1(t)&
                     + sm2coef * sm2(t)&
                     + sm3coef * sm3(t)&
                     + sm4coef * sm4(t)
            else
                self%VOD(t)=LIS_rc%udef
            endif

            if (self%VOD(t).ne.LIS_rc%udef.and.self%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=self%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODFM_VOD", vodval)
        vodval = self%VOD

    end subroutine VODFM_run_linear_model

    subroutine VODFM_run_svr_model(self, n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodfm_svr_model), intent(inout) :: self
        integer, intent(in) :: n

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
        real, pointer       :: vodval(:)
        integer             :: j
        real                :: intercept, actual_n_SV, dual_coef
        real                :: gam(5), support_vectors(5)
        real                :: kernelval
        logical             :: coefs_valid, values_valid


        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 1", sm1)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 2", sm2)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 3", sm3)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 4", sm4)

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

            intercept = self%intercept(t)
            actual_n_SV = self%actual_n_SV(t)
            gam = self%gam(:, t)

            ! For some pixels no VOD was available and therefore no forward
            ! model was fitted. No assimilation will take place over these
            ! pixels anyways, so it's no problem to not predict anything here
            coefs_valid = (.not.isnan(intercept).and.intercept.ne.LIS_rc%udef)

            ! normally the modelled values should not be invalid, but just to
            ! be on the safe side
            values_valid = (.not.isnan(lai(t)).and.lai(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm1(t)).and.sm1(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm2(t)).and.sm2(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm3(t)).and.sm3(t).ne.LIS_rc%udef&
                 .and..not.isnan(sm4(t)).and.sm4(t).ne.LIS_rc%udef)

            if (coefs_valid.and.values_valid) then
                ! VOD = \sum_j dual_coef[j] * K(x, s[j])
                self%VOD(t) = intercept
                do j=1, actual_n_SV
                    dual_coef = self%dual_coef(j, t)
                    support_vectors = self%support_vectors(:, j, t)
                    ! calculate kernel value
                    ! K = exp(-\sum_k gamma[k] * (x[k] - s[k])**2)
                    kernelval = exp(-gam(1) * (lai(t) - support_vectors(1))**2&
                         - gam(2) * (sm1(t) - support_vectors(2))**2&
                         - gam(3) * (sm2(t) - support_vectors(3))**2&
                         - gam(4) * (sm3(t) - support_vectors(4))**2&
                         - gam(5) * (sm4(t) - support_vectors(5))**2)
                    self%VOD(t) = self%VOD(t) + dual_coef * kernelval
                end do
            else
                self%VOD(t) = LIS_rc%udef
            endif

            if (self%VOD(t).ne.LIS_rc%udef.and.self%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=self%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODFM_VOD", vodval)
        vodval = self%VOD

    end subroutine VODFM_run_svr_model


    subroutine VODFM_run_X_model(self, n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodfm_X_model), intent(inout) :: self
        integer, intent(in) :: n

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: cwc(:), lai(:), psi(:), tveg(:)
        real, pointer       :: vodval(:)
        real                :: intercept, cwccoef, laicoef, laipsicoef, laisqtvegcoef
        logical             :: coefs_valid, values_valid


        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Canopy Water Content", cwc)
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Root Zone Soil Water Potential", psi)
        call getsfcvar(LIS_sfcState(n), "Vegetation Transpiration", tveg)

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

            intercept = self%intercept(t)
            cwccoef = self%cwccoef(t)
            laicoef = self%laicoef(t)
            laipsicoef = self%laipsicoef(t)
            laisqtvegcoef = self%laisqtvegcoef(t)

            ! For some pixels no VOD was available and therefore no forward
            ! model was fitted. No assimilation will take place over these
            ! pixels anyways, so it's no problem to not predict anything here
            coefs_valid = (.not.isnan(intercept).and.intercept.ne.LIS_rc%udef)

            ! normally the modelled values should not be invalid, but just to
            ! be on the safe side
            values_valid = (.not.isnan(lai(t)).and.lai(t).ne.LIS_rc%udef&
                 .and..not.isnan(cwc(t)).and.cwc(t).ne.LIS_rc%udef&
                 .and..not.isnan(psi(t)).and.psi(t).ne.LIS_rc%udef&
                 .and..not.isnan(tveg(t)).and.tveg(t).ne.LIS_rc%udef)

            if (coefs_valid.and.values_valid) then
                ! VOD = \sum_j dual_coef[j] * K(x, s[j])
                self%VOD(t) = intercept&
                    + cwccoef * cwc(t)&
                    + laicoef * lai(t)&
                    + laipsicoef * lai(t) * psi(t)&
                    + laisqtvegcoef * lai(t) * lai(t) * psi(t)
            else
                self%VOD(t) = LIS_rc%udef
            endif

            if (self%VOD(t).ne.LIS_rc%udef.and.self%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=self%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODFM_VOD", vodval)
        vodval = self%VOD

    end subroutine VODFM_run_X_model


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine getsfcvar(sfcState, varname, var)
        ! !USES:

        implicit none

        type(ESMF_State)      :: sfcState
        type(ESMF_Field)      :: varField
        character(len=*)      :: varname
        real, pointer         :: var(:)
        integer               :: status

        call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
        call LIS_verify(status, "Error in StateGet: VODFM_getsfcvar "//trim(varname))
        call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
        call LIS_verify(status, "Error in FieldGet: VODFM_getsfcvar "//trim(varname))

    end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: VODFM_output
!!!! \label{VODFM_output}
!!!!
!!!! !INTERFACE:
    subroutine VODFM_output(n)
        integer, intent(in) :: n
    end subroutine VODFM_output
#endif
end module VODFM_Mod
