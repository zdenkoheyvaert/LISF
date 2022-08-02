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

    ! The actual implementation is done via derived types from a generic model
    ! type. The derived types contain the parameters and provide a initialize
    ! and a run method.
    type, abstract :: vod_forward_model
        character*256 :: parameter_fname
    contains
        procedure(initialize_forward_model), deferred, pass(self) :: initialize
        procedure(run_forward_model), deferred, pass(self) :: run
    end type vod_forward_model

    abstract interface
        subroutine initialize_forward_model(self, n, fname)
            import vod_forward_model
            class(vod_forward_model), intent(inout) :: self
            integer, intent(in) :: n
            character(len=*), intent(in) :: fname
        end subroutine initialize_forward_model
    end interface

    abstract interface
        subroutine run_forward_model(self, n, vod, lai, sm1, sm2, sm3,&
            sm4)
            import vod_forward_model
            class(vod_forward_model), intent(in) :: self
            integer, intent(in) :: n
            real, intent(inout) :: vod(:)
            real, intent(in)    :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
        end subroutine run_forward_model
    end interface

    type, extends(vod_forward_model) :: vod_linear_model
        integer           :: ntimes
        real, allocatable :: intercept(:,:)
        real, allocatable :: laicoef(:,:)
        real, allocatable :: sm1coef(:,:)
        real, allocatable :: sm2coef(:,:)
        real, allocatable :: sm3coef(:,:)
        real, allocatable :: sm4coef(:,:)
        contains
            procedure, pass(self) :: initialize => VODFM_initialize_linear_model
            procedure, pass(self) :: run => VODFM_run_linear_model
    end type vod_linear_model

    type, extends(vod_forward_model) :: vod_svr_model
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
    end type vod_svr_model

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
    type, public ::  vodfm_type_dec
        character*100 :: modeltype
        character*256 :: parameter_fname
        class(vod_forward_model), pointer :: forward_model

        !-------output------------!
        real, allocatable :: VOD(:)
    end type vodfm_type_dec

    type(vodfm_type_dec), allocatable :: vodfm_struc(:)

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
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        implicit none

        integer :: rc, ios
        integer :: n, nid, ngrid, ngridId

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        write(LIS_logunit,*) "[INFO] Starting VODFM setup"
        !allocate memory for nest
        allocate(vodfm_struc(LIS_rc%nnest))

        do n=1,LIS_rc%nnest
            allocate(vodfm_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 1")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 2")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 3")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 4")

            call add_sfc_fields(n,LIS_forwardState(n),"VODFM_VOD")
        enddo


        ! read config from file
        call ESMF_ConfigFindLabel(LIS_config, "VODFM model type:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,vodfm_struc(n)%modeltype, rc=rc)
            call LIS_verify(rc, "VODFM model type: not defined")
            if (vodfm_struc(n)%modeltype .ne. "linear" .and. vodfm_struc(n)%modeltype .ne. "SVR") then
                write(LIS_logunit,*) "[ERR] VODFM model type must be 'linear' or 'svr'"
                call LIS_endrun
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config, "VODFM parameter file:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,vodfm_struc(n)%parameter_fname, rc=rc)
            call LIS_verify(rc, "VODFM parameter file: not defined")
        enddo


        ! read parameters/initialize model
        do n=1, LIS_rc%nnest

            ! try opening the parameter file
            write(LIS_logunit,*) '[INFO] Reading ',&
                trim(vodfm_struc(n)%parameter_fname)
            ios = nf90_open(path=trim(vodfm_struc(n)%parameter_fname),&
                mode=NF90_NOWRITE,ncid=nid)
            call LIS_verify(ios,'Error opening file '&
                //trim(vodfm_struc(n)%parameter_fname))

            ! check if ngrid is as expected
            ios = nf90_inq_dimid(nid, "ngrid", ngridId)
            call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
            ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
            call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
            if (ngrid /= LIS_rc%glbngrid_red(n)) then
                write(LIS_logunit, *) "[ERR] ngrid in "//trim(vodfm_struc(n)%parameter_fname)&
                     //" not consistent with expected ngrid: ", ngrid,&
                     " instead of ",LIS_rc%glbngrid_red(n)
                call LIS_endrun
            endif

            ! create the model instance
            call create_model(vodfm_struc(n), vodfm_struc(n)%modeltype, n,&
                vodfm_struc(n)%parameter_fname)
        enddo

        write(LIS_logunit,*) '[INFO] Finished VODFM setup'
#endif
    end subroutine VODFM_initialize

    subroutine VODFM_initialize_linear_model(self, n, fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vod_linear_model), intent(inout) :: self
        integer, intent(in) :: n
        character(len=*), intent(in) :: fname

        integer :: ios, nid, npatch
        integer :: ngridId, ntimesId, ngrid, ntimes

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        self%parameter_fname = fname

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

    subroutine VODFM_initialize_svr_model(self, n, fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vod_svr_model), intent(inout) :: self
        integer, intent(in) :: n
        character(len=*), intent(in) :: fname

        integer :: ios, nid, npatch
        integer :: ngridId, ngrid, n_SV, n_SVId

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        self%parameter_fname = fname

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
        allocate(self%dual_coef(npatch, n_SV))
        allocate(self%support_vectors(npatch, n_SV, 5))
        allocate(self%gam(npatch, 5))

        ! read the arrays from the file
        call read_1d_coef_from_file(n, nid, ngrid, "intercept_", self%intercept)
        call read_1d_coef_from_file(n, nid, ngrid, "n_SV_", self%actual_n_SV)
        call read_2d_coef_from_file(n, nid, ngrid, n_SV, "dual_coef_", self%dual_coef)
        call read_2d_coef_from_file(n, nid, ngrid, 5, "gamma_", self%gam)
        call read_3d_coef_from_file(n, nid, ngrid, n_SV, 5, "support_vectors_",&
             self%support_vectors)
#endif
    end subroutine VODFM_initialize_svr_model

    subroutine create_model(vodstruc, modeltype, n, fname)
        type(vodfm_type_dec), intent(inout) :: vodstruc
        character(len=*), intent(in) :: modeltype
        integer, intent(in) :: n
        character(len=*), intent(in) :: fname

        if (modeltype == "linear") then
            allocate(vod_linear_model :: vodstruc%forward_model)
        elseif (modeltype == "SVR") then
            allocate(vod_svr_model :: vodstruc%forward_model)
        else
            write(LIS_logunit, *) "[ERR] VODFM model type must be 'linear' or 'SVR'"
            call LIS_endrun
        endif
        call vodstruc%forward_model%initialize(n, fname)
    end subroutine create_model

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

    subroutine read_2d_coef_from_file(n, nid, ngrid, n2, varname, coef)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, ngrid, n2
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:)

        real, allocatable :: coef_file(:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        allocate(coef_file(ngrid, n2))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        do j=1,n2
            call LIS_convertVarToLocalSpace(n, coef_file(:,j), coef_grid)
            call gridvar_to_patchvar(&
                 n, LIS_rc%lsm_index, coef_grid, coef(:, j))
        enddo
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_2d_coef_from_file

    subroutine read_3d_coef_from_file(n, nid, ngrid, n2, n3, varname, coef)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, ngrid, n2, n3
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:,:)

        real, allocatable :: coef_file(:,:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j2, j3

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM requires NETCDF"
        call LIS_endrun
#else
        allocate(coef_file(ngrid, n2, n3))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        do j2=1,n2
            do j3=1,n3
                call LIS_convertVarToLocalSpace(n, coef_file(:,j2,j3), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(:, j2,j3))
             enddo
        enddo
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

    subroutine VODFM_run_linear_model(self, n, vod, lai, sm1, sm2, sm3, sm4)
        use LIS_histDataMod
        class(vod_linear_model), intent(in) :: self
        integer, intent(in) :: n
        real, intent(inout) :: vod(:)
        real, intent(in)    :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)

        integer :: t, timeidx, j
        real                :: intercept
        real                :: laicoef, sm1coef, sm2coef, sm3coef, sm4coef
        logical             :: coefs_valid, values_valid

        if (self%ntimes == 1) then
            timeidx = 1
        elseif (self%ntimes == 12) then
            timeidx = LIS_rc%mo
        else
            write(LIS_logunit, *) "[ERR] ntimes in "//trim(self%parameter_fname)&
                 //" must be 1 or 12, but is ", self%ntimes
            call LIS_endrun
        endif

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            intercept = self%intercept(t, timeidx)
            laicoef = self%laicoef(t, timeidx)
            sm1coef = self%sm1coef(t, timeidx)
            sm2coef = self%sm2coef(t, timeidx)
            sm3coef = self%sm3coef(t, timeidx)
            sm4coef = self%sm4coef(t, timeidx)

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
                vod(t) = intercept&
                     + laicoef * lai(t)&
                     + sm1coef * sm1(t)&
                     + sm2coef * sm2(t)&
                     + sm3coef * sm3(t)&
                     + sm4coef * sm4(t)
                ! if (self%VOD(t) .lt. 0.0) then
                !     self%VOD(t) = LIS_rc%udef
                ! endif
            else
                vod(t)=LIS_rc%udef
            endif

            if (vod(t).ne.LIS_rc%udef.and.vod(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vod(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo
    end subroutine VODFM_run_linear_model

    subroutine VODFM_run_svr_model(self, n, vod, lai, sm1, sm2, sm3, sm4)
        ! calculates VOD based on modelled input values
        use LIS_histDataMod
        implicit none

        class(vod_svr_model), intent(in) :: self
        integer, intent(in) :: n
        real, intent(inout) :: vod(:)
        real, intent(in)    :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)

        integer :: t, timeidx, j
        real                :: intercept
        real                :: actual_n_SV
        real, pointer       :: dual_coef(:), gam(:), support_vectors(:,:)
        real                :: kernelval
        logical             :: coefs_valid, values_valid


        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            intercept = self%intercept(t)
            actual_n_SV = self%actual_n_SV(t)
            dual_coef = self%dual_coef(t, :)
            gam = self%gam(t, :)
            support_vectors = self%support_vectors(t, :, :)

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
                vod(t) = 0.0
                do j=1, actual_n_SV
                    ! calculate kernel value
                    ! K = exp(-\sum_k gamma[k] * (x[k] - s[k])**2)
                    kernelval = exp(-gam(1) * (lai(t) - support_vectors(j,1))**2&
                         - gam(2) * (sm1(t) - support_vectors(j,2))**2&
                         - gam(3) * (sm2(t) - support_vectors(j,3))**2&
                         - gam(4) * (sm3(t) - support_vectors(j,4))**2&
                         - gam(5) * (sm4(t) - support_vectors(j,5))**2)
                    vod(t) = vod(t) + dual_coef(j) * kernelval
                end do
            else
                vod(t) = LIS_rc%udef
            endif

            if (vod(t).ne.LIS_rc%udef.and.vod(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vod(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo
    end subroutine VODFM_run_svr_model


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

        integer             :: t, timeidx
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
        real                :: intercept, laicoef, sm1coef, sm2coef, sm3coef, sm4coef
        real, pointer       :: vodval(:)
        logical             :: coefs_valid, values_valid

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 1", sm1)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 2", sm2)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 3", sm3)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 4", sm4)


        call vodfm_struc(n)%forward_model%run(n, vodfm_struc(n)%VOD, lai, sm1, sm2, sm3, sm4)

        call getsfcvar(LIS_forwardState(n),"VODFM_VOD", vodval)
        vodval = vodfm_struc(n)%VOD

    end subroutine VODFM_run


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
