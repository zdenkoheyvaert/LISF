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
    !    VODFM predictors:
    !        "mech" or "stat". "mech": (CWC, LAI, PSI*LAI, VPD*LAI**2),
    !        "stat": (CWC, LAI, SM1, SM2, SM3, SM4, VPD)
    !    VODFM parameter file:
    !        Pathname to the parameter file
    !
    !    The linear model calculates VOD as linear function of the predictors.
    !    To use it, you need to provide a netCDF file with variables,
    !    "intercept_" (ngrid,) and "coef_" (npred,ngrid (C-notation)).
    !
    !    The SVR based model calculates VOD via:
    !
    !        VOD[lat,lon] = \sum_{j}^{n_SV_[lat,lon]} dual_coef_[lat,lon,j]
    !                       * K(x[lat,lon], support_vectors_[lat,lon,j])
    !
    !    with x = (LAI, SM1, .., SM4)^T
    !    and K(x, s) = exp(-\sum_k gamma_[k] (x[k] - s[k])**2)
    !
    !    To use the SVR model, you need to provide a netCDF file with
    !    npredictors variables: intercept_, n_SV_, dual_coef_,
    !    support_vectors_, gamma_.
    !    The support vectors are npredictors-d vectors.
    ! 
    !    The variables on file have to have the following specifications:
    !    - intercept_: float, (ngrid,)
    !    - n_SV_: float, (ngrid,)
    !    - dual_coef_: float, (ngrid, n_SV) (n_SV == maximum number of support vectors)
    !    - gamma_: float, (ngrid, npredictors)
    !    - support_vectors_: float, (ngrid, n_SV, npredictors)
    !
    ! !HISTORY:
    ! 02 Mar 2022: Samuel Scherrer; initial contribution based on WCMRTM
    ! 25 Aug 2022: Samuel Scherrer; support for multiple models
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
        character*10   :: predictors
        integer       :: npred
        character*256 :: parameter_fname
        !-------output------------!
        real, allocatable :: VOD(:)
        contains
        procedure, pass(self) :: initialize => VODFM_initialize_default
        procedure, pass(self) :: run => VODFM_run_default
    end type vodfm_type_dec

    type, extends(vodfm_type_dec) :: vodfm_linear_model
        real, allocatable :: intercept(:)
        real, allocatable :: coefs(:, :)
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
        character*10 :: predictors(LIS_rc%nnest)
        character*256 :: parameter_fname(LIS_rc%nnest)

        write(LIS_logunit,*) "[INFO] Starting VODFM setup"


        ! read config from file
        call ESMF_ConfigFindLabel(LIS_config, "VODFM model type:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, modeltype(n), rc=rc)
            call LIS_verify(rc, "VODFM model type: not defined")
            if (modeltype(n) .ne. "linear"&
                 .and. modeltype(n) .ne. "SVR") then
                write(LIS_logunit, *)&
                     "[ERR] VODFM model type must be 'linear' or 'SVR'"
                call LIS_endrun
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config, "VODFM predictors:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, predictors(n), rc=rc)
            call LIS_verify(rc, "VODFM predictors: not defined")
            if (predictors(n) .ne. "mech"&
                 .and. predictors(n) .ne. "stat_sm"&
                 .and. predictors(n) .ne. "stat_rzsm") then
                write(LIS_logunit, *)&
                     "[ERR] VODFM predictors must be 'mech', 'stat_sm' or 'stat_rzsm'"
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
        else
            write(LIS_logunit, *)&
                 "[ERR] VODFM model type must be 'linear' or 'SVR'"
            call LIS_endrun
        endif

        do n=1,LIS_rc%nnest
            vodfm_struc(n)%predictors = predictors(n)
            if (predictors(n) == "mech") then
                vodfm_struc(n)%npred = 5
            else if (predictors(n) == "stat_sm") then
                vodfm_struc(n)%npred = 7
            else if (predictors(n) == "stat_rzsm") then
                vodfm_struc(n)%npred = 4
            endif
            vodfm_struc(n)%parameter_fname = parameter_fname(n)
            allocate(vodfm_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Canopy Water Content")
            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 1")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 2")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 3")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 4")
            call add_sfc_fields(n,LIS_sfcState(n), "Root Zone Soil Water Potential")
            call add_sfc_fields(n,LIS_sfcState(n), "Root Zone Soil Moisture")
            call add_sfc_fields(n,LIS_sfcState(n), "Canopy Vapor Pressure Deficit")
            call add_sfc_fields(n,LIS_sfcState(n), "Canopy Temperature")
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

        integer :: ios, nid, npatch, npred
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
        allocate(self%coefs(self%npred, npatch))

        call read_1d_coef_from_file(n, nid, ngrid, "intercept_", self%intercept)
        call read_2d_coef_from_file(n, nid, self%npred, ngrid, "coef_", &
            self%coefs, ngrid_first=.false.)
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
        allocate(self%support_vectors(self%npred, n_SV, npatch))
        allocate(self%gam(self%npred, npatch))

        ! read the arrays from the file
        call read_1d_coef_from_file(n, nid, ngrid, "intercept_", self%intercept)
        call read_1d_coef_from_file(n, nid, ngrid, "n_SV_", self%actual_n_SV)
        call read_2d_coef_from_file(n, nid, n_SV, ngrid, "dual_coef_", &
            self%dual_coef, ngrid_first=.false.)
        call read_2d_coef_from_file(n, nid, self%npred, ngrid, "gamma_", self%gam, &
            ngrid_first=.false.)
        call read_3d_coef_from_file(n, nid, self%npred, n_SV, ngrid, "support_vectors_",&
             self%support_vectors, ngrid_first=.false.)
#endif
    end subroutine VODFM_initialize_svr_model

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

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:), psi(:), cwc(:), tc(:), vpd(:)
        real, pointer       :: vodval(:), rzsm(:)
        real                :: intercept
        real                :: eps_water, deltaT
        real                :: pred(vodfm_struc(n)%npred)

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Canopy Water Content", cwc)
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 1", sm1)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 2", sm2)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 3", sm3)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 4", sm4)
        call getsfcvar(LIS_sfcState(n), "Root Zone Soil Water Potential", psi)
        call getsfcvar(LIS_sfcState(n), "Root Zone Soil Moisture", rzsm)
        call getsfcvar(LIS_sfcState(n), "Canopy Temperature", tc)
        call getsfcvar(LIS_sfcState(n), "Canopy Vapor Pressure Deficit", vpd)

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

            if (vodfm_struc(n)%predictors .eq. "mech") then
                pred(1) = cwc(t)
                pred(2) = lai(t)
                pred(3) = psi(t) * lai(t)
                pred(4) = vpd(t) * lai(t)**2
                pred(5) = tc(t)
            else if (vodfm_struc(n)%predictors .eq. "stat_sm") then
                pred(1) = cwc(t)
                pred(2) = lai(t)
                pred(3) = sm1(t)
                pred(4) = sm2(t)
                pred(5) = sm3(t)
                pred(6) = sm4(t)
                pred(7) = vpd(t)
            else if (vodfm_struc(n)%predictors .eq. "stat_rzsm") then
                pred(1) = cwc(t)
                pred(2) = lai(t)
                pred(3) = rzsm(t)
                pred(4) = vpd(t)
            endif

            call vodfm_struc(n)%run(n, t, pred)

            if (vodfm_struc(n)%VOD(t).ne.LIS_rc%udef.and.vodfm_struc(n)%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
                vodfm_struc(n)%VOD(t) = LIS_rc%udef
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vodfm_struc(n)%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODFM_VOD", vodval)
        vodval = vodfm_struc(n)%VOD


    end subroutine VODFM_run

    subroutine VODFM_run_default(self, n, t, pred)
        class(vodfm_type_dec), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%npred)
        write(LIS_logunit,*) "[ERR] VODFM should use linear or svr model"
        call LIS_endrun
    end subroutine VODFM_run_default

    subroutine VODFM_run_linear_model(self, n, t, pred)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodfm_linear_model), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%npred)

        real                :: intercept
        real                :: coefs(self%npred)
        logical             :: coefs_valid

        intercept = self%intercept(t)
        coefs = self%coefs(:, t)

        ! For some pixels no VOD was available and therefore no forward
        ! model was fitted. No assimilation will take place over these
        ! pixels anyways, so it's no problem to not predict anything here
        coefs_valid = (.not.isnan(intercept).and.intercept.ne.LIS_rc%udef)
        if (coefs_valid) then
            self%VOD(t) = intercept + sum(coefs * pred)
        else
            self%VOD(t)=LIS_rc%udef
        endif
    end subroutine VODFM_run_linear_model

    subroutine VODFM_run_svr_model(self, n, t, pred)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodfm_svr_model), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%npred)

        integer             :: j
        real                :: intercept, actual_n_SV, dual_coef, vod
        real                :: gam(self%npred), support_vectors(self%npred)
        real                :: kernelval
        logical             :: coefs_valid, values_valid


        intercept = self%intercept(t)
        actual_n_SV = self%actual_n_SV(t)
        gam = self%gam(:, t)

        ! For some pixels no VOD was available and therefore no forward
        ! model was fitted. No assimilation will take place over these
        ! pixels anyways, so it's no problem to not predict anything here
        coefs_valid = (.not.isnan(intercept).and.intercept.ne.LIS_rc%udef)

        if (coefs_valid) then
            ! VOD = \sum_j dual_coef[j] * K(x, s[j])
            vod = intercept
            do j=1, actual_n_SV
                dual_coef = self%dual_coef(j, t)
                support_vectors = self%support_vectors(:, j, t)
                ! calculate kernel value
                ! K = exp(-\sum_k gamma[k] * (x[k] - s[k])**2)
                kernelval = exp(-sum(gam * (pred - support_vectors)**2))
                vod = vod + dual_coef * kernelval
            end do
            ! if (vod.lt.0) then
            !     ! write(LIS_logunit, *) "[WARN] VOD lower than 0"
            !     vod = LIS_rc%udef
            ! endif
        else
            vod = LIS_rc%udef
        endif
        self%VOD(t) = vod
    end subroutine VODFM_run_svr_model


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
