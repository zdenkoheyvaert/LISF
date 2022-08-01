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
    type vod_forward_model
        character*256 :: parameter_fname
    contains
        procedure :: initialize => VODFM_initialize_vod_forward_model
        procedure :: run => VODFM_run_vod_forward_model
    end type vod_forward_model

    type, extends(vod_forward_model) :: vod_linear_model
        integer           :: ntimes
        real, allocatable :: intercept(:,:)
        real, allocatable :: laicoef(:,:)
        real, allocatable :: sm1coef(:,:)
        real, allocatable :: sm2coef(:,:)
        real, allocatable :: sm3coef(:,:)
        real, allocatable :: sm4coef(:,:)
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
        implicit none

        integer :: rc
        integer :: n

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
            if (vodfm_struc(n)%modeltype .ne. "linear" .and. vodfm_struc(n)%modeltype .ne. "SVR")
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
            call vodfm_struc(n)%forward_model%initialize(vodfm_struc(n)%parameter_fname)
        enddo

        write(LIS_logunit,*) '[INFO] Finished VODFM setup'
    end subroutine VODFM_initialize


    ! This is the type bound procedure to set up the specific model
    ! implementations which is called in the public initialize functions
    ! that the VODFM module exposes.
    subroutine VODFM_initialize_vod_forward_model(forward_model, fname)
        ! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        implicit none

        class(vod_forward_model) :: forward_model
        character(len=*) :: fname

        integer :: ios, nid, m
        integer :: ngridId, ntimesId, ngrid, ntimes, n_SV

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODFM_initialize requires NETCDF"
        call LIS_endrun
#else

        forward_model%parameter_fname = fname

        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(forward_model%parameter_fname)
        ios = nf90_open(path=trim(forward_model%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(forward_model%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(forward_model%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! model specific setup
        select type(forward_model)
        type is (vod_linear_model)
            ! read ntimes
            ios = nf90_inq_dimid(nid, "ntimes", ntimesId)
            call LIS_verify(ios, "Error nf90_inq_varid: ntimes")
            ios = nf90_inquire_dimension(nid, ntimesId, len=ntimes)
            call LIS_verify(ios, "Error nf90_inquire_dimension: ntimes")
            forward_model%ntimes = ntimes

            ! now that we have ntimes, we can allocate the coefficient arrays
            ! for the nest
            m = LIS_rc%lsm_index
            allocate(forward_model%intercept(LIS_rc%npatch(n, m), ntimes))
            allocate(forward_model%laicoef(LIS_rc%npatch(n, m), ntimes))
            allocate(forward_model%sm1coef(LIS_rc%npatch(n, m), ntimes))
            allocate(forward_model%sm2coef(LIS_rc%npatch(n, m), ntimes))
            allocate(forward_model%sm3coef(LIS_rc%npatch(n, m), ntimes))
            allocate(forward_model%sm4coef(LIS_rc%npatch(n, m), ntimes))

            call read_2d_coef_from_file(nid, ngrid, ntimes, "intercept", &
                forward_model%intercept)
            call read_2d_coef_from_file(nid, ngrid, ntimes, "LAI_coef", &
                forward_model%laicoef)
            call read_2d_coef_from_file(nid, ngrid, ntimes, "SM1_coef", &
                forward_model%sm1coef)
            call read_2d_coef_from_file(nid, ngrid, ntimes, "SM2_coef", &
                forward_model%sm2coef)
            call read_2d_coef_from_file(nid, ngrid, ntimes, "SM3_coef", &
                forward_model%sm3coef)
            call read_2d_coef_from_file(nid, ngrid, ntimes, "SM4_coef", &
                forward_model%sm4coef)
        type is (vod_svr_model)
            ! read n_SV, because we don't know it's length a priori
            ios = nf90_inq_dimid(nid, "n_SV", n_SVId)
            call LIS_verify(ios, "Error nf90_inq_varid: n_SV")
            ios = nf90_inquire_dimension(nid, n_SVId, len=ntimes)
            call LIS_verify(ios, "Error nf90_inquire_dimension: n_SV")
            forward_model%n_SV = n_SV

            ! now that we have n_SV, we can allocate the coefficient arrays
            ! for the nest
            m = LIS_rc%lsm_index
            allocate(forward_model%intercept(LIS_rc%npatch(n, m)))
            allocate(forward_model%actual_n_SV(LIS_rc%npatch(n, m)))
            allocate(forward_model%dual_coef(LIS_rc%npatch(n, m), n_SV))
            allocate(forward_model%support_vectors(LIS_rc%npatch(n, m), n_SV, 5))
            allocate(forward_model%gam(LIS_rc%npatch(n, m), 5))

            ! read the arrays from the file
            call read_1d_coef_from_file(nid, ngrid, "intercept_", forward_model%intercept)
            call read_1d_coef_from_file(nid, ngrid, "n_SV_", forward_model%actual_n_SV)
            call read_2d_coef_from_file(nid, ngrid, n_SV, "dual_coef_", forward_model%dual_coef)
            call read_2d_coef_from_file(nid, ngrid, 5, "gamma_", forward_model%gam)
            call read_3d_coef_from_file(nid, ngrid, n_SV, 5, "support_vectors_",&
                 forward_model%support_vectors)
        type default
            write(LIS_logunit, *) "[ERR] invalid type in VOD forward model run routine"
            call LIS_endrun
        end select

    contains

        subroutine read_1d_coef_from_file(nid, ngrid, varname, coef)
            use LIS_historyMod, only: LIS_convertVarToLocalSpace
            integer, intent(in) :: nid, ngrid
            character(len=*), intent(in) :: varname
            real, intent(inout) :: coef(:)

            real, allocatable :: coef_file(:)
            real, allocatable :: coef_grid(:)
            integer :: ios, varid, j

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
        end subroutine read_1d_coef_from_file

        subroutine read_2d_coef_from_file(nid, ngrid, n2, varname, coef)
            use LIS_historyMod, only: LIS_convertVarToLocalSpace
            integer, intent(in) :: nid, ngrid, n2
            character(len=*), intent(in) :: varname
            real, intent(inout) :: coef(:,:)

            real, allocatable :: coef_file(:,:)
            real, allocatable :: coef_grid(:)
            integer :: ios, varid, j

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
        end subroutine read_2d_coef_from_file

        subroutine read_3d_coef_from_file(nid, ngrid, n2, n3, varname, coef)
            use LIS_historyMod, only: LIS_convertVarToLocalSpace
            integer, intent(in) :: nid, ngrid, n2, n3
            character(len=*), intent(in) :: varname
            real, intent(inout) :: coef(:,:,:)

            real, allocatable :: coef_file(:,:,:)
            real, allocatable :: coef_grid(:,:)
            integer :: ios, varid, j

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
#endif
    end subroutine VODFM_initialize_vod_forward_model
    !!--------------------------------------------------------------------------------

    subroutine VODFM_run_vod_forward_model(forward_model, n, lai, sm1, sm2, sm3, sm4)
        ! calculates VOD based on modelled input values
        implicit none

        class(vod_forward_model) :: forward_model
        integer, intent(in) :: n
        real, intent(in)    :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)

        integer :: t, timeidx
        real                :: intercept
        real                :: 
        real                :: intercept
        real                :: laicoef, sm1coef, sm2coef, sm3coef, sm4coef
        real                :: actual_n_SV, dual_coef(:), gam(:), support_vectors(:,:)
        real                :: kernelval
        real, pointer       :: vodval(:)
        logical             :: coefs_valid, values_valid


        select type
        type is (vod_linear_model)
            if (forward_model%ntimes == 1) then
                timeidx = 1
            elseif (forward_model%ntimes == 12) then
                timeidx = LIS_rc%mo
            else
                write(LIS_logunit, *) "[ERR] ntimes in "//trim(forward_model%parameter_fname)&
                     //" must be 1 or 12, but is ", forward_model%ntimes
                call LIS_endrun
        end select

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        endif
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            select type (forward_model)
            type is (linear_fm_params)
                intercept = forward_model%intercept(t, timeidx)
                laicoef = forward_model%laicoef(t, timeidx)
                sm1coef = forward_model%sm1coef(t, timeidx)
                sm2coef = forward_model%sm2coef(t, timeidx)
                sm3coef = forward_model%sm3coef(t, timeidx)
                sm4coef = forward_model%sm4coef(t, timeidx)

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
                    forward_model%VOD(t) = intercept&
                         + laicoef * lai(t)&
                         + sm1coef * sm1(t)&
                         + sm2coef * sm2(t)&
                         + sm3coef * sm3(t)&
                         + sm4coef * sm4(t)
                    ! if (forward_model%VOD(t) .lt. 0.0) then
                    !     forward_model%VOD(t) = LIS_rc%udef
                    ! endif
                else
                    forward_model%VOD(t)=LIS_rc%udef
                endif
            type is (svr_fm_params)
                intercept = forward_model%intercept(t)
                actual_n_SV = forward_model%actual_n_SV(t)
                dual_coef = forward_model%dual_coef(t, :)
                gam = forward_model%gam(t, :)
                support_vectors = forward_model%support_vectors(t, :, :)

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
                    vod = 0.0
                    do j=1, actual_n_SV
                        ! calculate kernel value
                        ! K = exp(-\sum_k gamma[k] * (x[k] - s[k])**2)
                        kernelval = exp(-gam(1) * (lai(t) - support_vectors(j,1))**2&
                             - gam(2) * (sm1(t) - support_vectors(j,2))**2&
                             - gam(3) * (sm2(t) - support_vectors(j,3))**2&
                             - gam(4) * (sm3(t) - support_vectors(j,4))**2&
                             - gam(5) * (sm4(t) - support_vectors(j,5))**2)
                        vod = vod + dual_coef(j) * kernelval
                    end do
                    forward_model%VOD(t) = vod
                else
                    forward_model%VOD(t)=LIS_rc%udef
                endif
            type default
                write(LIS_logunit, *) "[ERR] invalid type in VOD forward model run routine"
                call LIS_endrun
            end select

            if (forward_model%VOD(t).ne.LIS_rc%udef.and.forward_model%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=forward_model%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

    end subroutine VODFM_run_vod_forward_model

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


        call vodfm_struc(n)%forward_model%run(n, lai, sm1, sm2, sm3, sm4)

        call getsfcvar(LIS_forwardState(n),"VODFM_VOD", vodval)
        vodval = vodfm_struc(n)%forward_model%VOD

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
