!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module VODLinFM_Mod
    !BOP
    !
    ! !MODULE: VODLinFM_Mod
    !
    ! !DESCRIPTION:
    !    This modules provides implements a linear forward model for VOD.
    !    Although it's not a real RTM it programmatically works similarly and is
    !    therefore implemented analogously to a RTM.
    !
    ! To use this module, you need to provide a netCDF file with 6 variables,
    ! intercept, LAI_coef, and SM1_coef to SM4_coef. Each of these variables
    ! must have the dimensions (ntimes, ngrid), where ntimes is either 1 or 12,
    ! depending on whether a single model for the full year should be used, or
    ! monthly varying models, and ngrid is the dimension of the model tile
    ! space (i.e. a flat 1D array with only land points such that east_west is
    ! changing fastest).
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
    public :: VODLinFM_initialize
    public :: VODLinFM_f2t
    public :: VODLinFM_run
    public :: VODLinFM_output
    public :: VODLinFM_geometry
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: vodlinfm_struc
    !EOP
    type, public ::  vodlinfm_type_dec

        character*256     :: parameter_fname

        integer           :: ntimes
        real, allocatable :: intercept(:,:)
        real, allocatable :: laicoef(:,:)
        real, allocatable :: sm1coef(:,:)
        real, allocatable :: sm2coef(:,:)
        real, allocatable :: sm3coef(:,:)
        real, allocatable :: sm4coef(:,:)

        !-------output------------!
        real, allocatable :: VOD(:)
    end type vodlinfm_type_dec

    type(vodlinfm_type_dec), allocatable :: vodlinfm_struc(:)

    SAVE

contains
    !BOP
    !
    ! !ROUTINE: VODLinFM_initialize
    ! \label{VODLinFM_initialize}
    !
    ! !INTERFACE:
    subroutine VODLinFM_initialize()
        ! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif

        ! !DESCRIPTION:
        !
        !  This routine creates the datatypes and allocates memory for noahMP3.6-specific
        !  variables. It also invokes the routine to read the runtime specific options
        !  for noahMP3.6 from the configuration file.
        !
        !  The routines invoked are:
        !  \begin{description}
        !   \item[readVODLinFMcrd](\ref{readVODLinFMcrd}) \newline
        !
        !EOP
        implicit none

        integer :: rc
        integer :: n,t, m
        integer :: ierr, ios, nid
        integer :: ngridId, ntimesId, ngrid, ntimes
        integer , parameter :: OPEN_OK = 0
        character*128 :: message

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODLinFM_initialize requires NETCDF"
        call LIS_endrun
#else

        !allocate memory for nest
        allocate(vodlinfm_struc(LIS_rc%nnest))

        call ESMF_ConfigFindLabel(LIS_config, "VODLinFM parameter file:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,vodlinfm_struc(n)%parameter_fname, rc=rc)
            call LIS_verify(rc, "VODLinFM parameter file: not defined")
        enddo

        do n=1,LIS_rc%nnest
            allocate(vodlinfm_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 1")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 2")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 3")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 4")

        enddo



        do n=1, LIS_rc%nnest
            write(LIS_logunit,*) '[INFO] Reading ',&
                trim(vodlinfm_struc(n)%parameter_fname)

            ios = nf90_open(path=trim(vodlinfm_struc(n)%parameter_fname),&
                mode=NF90_NOWRITE,ncid=nid)
            call LIS_verify(ios,'Error opening file '&
                //trim(vodlinfm_struc(n)%parameter_fname))

            ! check if ngrid is as expected
            ios = nf90_inq_dimid(nid, "ngrid", ngridId)
            call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
            ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
            call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
            if (ngrid /= LIS_rc%glbngrid_red(n)) then
                write(LIS_logunit, *) "[ERR] ngrid in "//trim(vodlinfm_struc(n)%parameter_fname)&
                     //" not consistent with expected ngrid: ", ngrid,&
                     " instead of ",LIS_rc%glbngrid_red(n)
                call LIS_endrun
            endif

            ! read ntimes
            ios = nf90_inq_dimid(nid, "ntimes", ntimesId)
            call LIS_verify(ios, "Error nf90_inq_varid: ntimes")
            ios = nf90_inquire_dimension(nid, ntimesId, len=ntimes)
            call LIS_verify(ios, "Error nf90_inquire_dimension: ntimes")
            vodlinfm_struc(n)%ntimes = ntimes

            ! now that we have ntimes, we can allocate the coefficient arrays
            ! for the nest
            m = LIS_rc%lsm_index
            allocate(vodlinfm_struc(n)%intercept(LIS_rc%npatch(n, m), ntimes))
            allocate(vodlinfm_struc(n)%laicoef(LIS_rc%npatch(n, m), ntimes))
            allocate(vodlinfm_struc(n)%sm1coef(LIS_rc%npatch(n, m), ntimes))
            allocate(vodlinfm_struc(n)%sm2coef(LIS_rc%npatch(n, m), ntimes))
            allocate(vodlinfm_struc(n)%sm3coef(LIS_rc%npatch(n, m), ntimes))
            allocate(vodlinfm_struc(n)%sm4coef(LIS_rc%npatch(n, m), ntimes))

            call read_coef_from_file(nid, ngrid, ntimes, "intercept", &
                vodlinfm_struc(n)%intercept)
            call read_coef_from_file(nid, ngrid, ntimes, "LAI_coef", &
                vodlinfm_struc(n)%laicoef)
            call read_coef_from_file(nid, ngrid, ntimes, "SM1_coef", &
                vodlinfm_struc(n)%sm1coef)
            call read_coef_from_file(nid, ngrid, ntimes, "SM2_coef", &
                vodlinfm_struc(n)%sm1coef)
            call read_coef_from_file(nid, ngrid, ntimes, "SM3_coef", &
                vodlinfm_struc(n)%sm1coef)
            call read_coef_from_file(nid, ngrid, ntimes, "SM4_coef", &
                vodlinfm_struc(n)%sm1coef)
        enddo

        do n=1,LIS_rc%nnest
            call add_fields_toState(n,LIS_forwardState(n),"VOD")
        enddo

    contains

        subroutine read_coef_from_file(nid, ngrid, ntimes, varname, coef)
            use LIS_historyMod, only: LIS_convertVarToLocalSpace
            integer, intent(in) :: nid, ngrid, ntimes
            character(len=*), intent(in) :: varname
            real, intent(inout) :: coef(:,:)

            real, allocatable :: coef_file(:,:)
            real, allocatable :: coef_grid(:)
            integer :: ios, varid, j

            allocate(coef_file(ngrid, ntimes))
            ios = nf90_inq_varid(nid, trim(varname), varid)
            call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
            ios = nf90_get_var(nid, varid, coef_file)
            call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
            allocate(coef_grid(LIS_rc%ngrid(n)))
            do j=1,ntimes
                call LIS_convertVarToLocalSpace(n, coef_file(:,j), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(:, j))
            enddo
            deallocate(coef_grid)
            deallocate(coef_file)

        end subroutine read_coef_from_file

        subroutine gridvar_to_patchvar(n,m,gvar,tvar)
            ! Converts a variable in local gridspace (length LIS_rc%ngrid(n))
            ! to local patch space (length LIS_rc%npatch(n,m))

            implicit none

            integer, intent(in) :: n 
            integer, intent(in) :: m
            real, intent(in)    :: gvar(LIS_rc%ngrid(n))
            real, intent(inout) :: tvar(LIS_rc%npatch(n,m))
            integer             :: t,r,c

            do t=1,LIS_rc%npatch(n,m)
                r = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
                c = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
                tvar(t) = gvar(LIS_domain(n)%gindex(c,r))
            enddo

        end subroutine gridvar_to_patchvar

        subroutine add_fields_toState(n, inState,varname)

            use LIS_logMod,   only : LIS_verify
            use LIS_coreMod,  only : LIS_vecTile

            implicit none

            integer            :: n
            type(ESMF_State)   :: inState
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
            call LIS_verify(status, "Error in field_create of "//trim(varname))

            call ESMF_StateAdd(inState, (/varField/), rc=status)
            call LIS_verify(status, "Error in StateAdd of "//trim(varname))

        end subroutine add_fields_toState
#endif
    end subroutine VODLinFM_initialize
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


    subroutine VODLinFM_f2t(n)

        implicit none

        integer, intent(in)    :: n

    end subroutine VODLinFM_f2t


    subroutine VODLinFM_geometry(n)
        implicit none
        integer, intent(in)    :: n

    end subroutine VODLinFM_geometry
    !Do nothing for now
    subroutine VODLinFM_run(n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        integer, intent(in) :: n

        integer             :: t, timeidx
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
        real                :: intercept, laicoef, sm1coef, sm2coef, sm3coef, sm4coef

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 1", sm1)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 2", sm2)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 3", sm3)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 4", sm4)


        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        if (vodlinfm_struc(n)%ntimes == 1) then
            timeidx = 1
        elseif (vodlinfm_struc(n)%ntimes == 12) then
            timeidx = LIS_rc%mo
        else
            write(LIS_logunit, *) "[ERR] ntimes in "//trim(vodlinfm_struc(n)%parameter_fname)&
                 //" must be 1 or 12, but is ", vodlinfm_struc(n)%ntimes
            call LIS_endrun
        endif
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            intercept = vodlinfm_struc(n)%intercept(t, t)
            laicoef = vodlinfm_struc(n)%laicoef(t, t)
            sm1coef = vodlinfm_struc(n)%sm1coef(t, t)
            sm2coef = vodlinfm_struc(n)%sm2coef(t, t)
            sm3coef = vodlinfm_struc(n)%sm3coef(t, t)
            sm4coef = vodlinfm_struc(n)%sm4coef(t, t)

            if(.not.isnan(lai(t)).and.lai(t).ne.LIS_rc%udef&
                 .and. .not.isnan(sm1(t)).and.sm1(t).ne.LIS_rc%udef&
                 .and. .not.isnan(sm2(t)).and.sm2(t).ne.LIS_rc%udef&
                 .and. .not.isnan(sm3(t)).and.sm3(t).ne.LIS_rc%udef&
                 .and. .not.isnan(sm4(t)).and.sm4(t).ne.LIS_rc%udef) then
                vodlinfm_struc(n)%VOD(t) = intercept + laicoef * lai(t)&
                     + sm1coef * sm1(t)&
                     + sm2coef * sm2(t)&
                     + sm3coef * sm3(t)&
                     + sm4coef * sm4(t)
            else
                vodlinfm_struc(n)%VOD(t)=LIS_rc%udef
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vodlinfm_struc(n)%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="+")
        enddo

    end subroutine VODLinFM_run


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
        call LIS_verify(status, "Error in StateGet: VODLinFM_getsfcvar "//trim(varname))
        call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
        call LIS_verify(status, "Error in FieldGet: VODLinFM_getsfcvar "//trim(varname))

    end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: VODLinFM_output
!!!! \label{VODLinFM_output}
!!!!
!!!! !INTERFACE:
    subroutine VODLinFM_output(n)
        integer, intent(in) :: n
    end subroutine VODLinFM_output
#endif
end module VODLinFM_Mod
