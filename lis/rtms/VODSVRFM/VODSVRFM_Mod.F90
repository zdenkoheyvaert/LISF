!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module VODSVRFM_Mod
    !BOP
    !
    ! !MODULE: VODSVRFM_Mod
    !
    ! !DESCRIPTION:
    !    This modules provides implements a SVR-based forward model for VOD.
    !    Although it's not a real RTM it programmatically works similarly and is
    !    therefore implemented analogously to a RTM.
    !
    ! To use this module, you need to provide a netCDF file with 5 variables,
    ! intercept_, n_SV_, dual_coef_, support_vectors_, gamma_.
    ! The support vectors are 5-d vectors with dimensions (LAI, SM1, SM2, SM3,
    ! SM4).
    !
    ! The predicted VOD is calculated via
    !
    ! VOD[lat,lon] = \sum_{j}^{n_SV_[lat,lon]} dual_coef_[lat,lon,j]
    !                   * K(x[lat,lon], support_vectors_[lat,lon,j])
    !
    ! with x = (LAI, SM1, .., SM4)^T
    ! and K(x, s) = exp(-\sum_k gamma_[k] (x[k] - s[k])**2)
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
    public :: VODSVRFM_initialize
    public :: VODSVRFM_f2t
    public :: VODSVRFM_run
    public :: VODSVRFM_output
    public :: VODSVRFM_geometry
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: vodsvrfm_struc
    !EOP
    type, public ::  vodsvrfm_type_dec

        character*256     :: parameter_fname

        integer           :: n_SV
        real, allocatable :: intercept(:)
        real, allocatable :: actual_n_SV(:)
        ! dual coef has n_SV as additional dim
        real, allocatable :: dual_coef(:,:)
        ! for each pixel, this is a matrix n_SV x 5
        real, allocatable :: support_vectors(:,:,:)
        ! for each pixel, this is a vector of length 5
        real, allocatable :: gam(:,:)

        !-------output------------!
        real, allocatable :: VOD(:)
    end type vodsvrfm_type_dec

    type(vodsvrfm_type_dec), allocatable :: vodsvrfm_struc(:)

    SAVE

contains
    !BOP
    !
    ! !ROUTINE: VODSVRFM_initialize
    ! \label{VODSVRFM_initialize}
    !
    ! !INTERFACE:
    subroutine VODSVRFM_initialize()
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
        !   \item[readVODSVRFMcrd](\ref{readVODSVRFMcrd}) \newline
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
        write(LIS_logunit,*) "[ERR] VODSVRFM_initialize requires NETCDF"
        call LIS_endrun
#else

        !allocate memory for nest
        allocate(vodsvrfm_struc(LIS_rc%nnest))

        call ESMF_ConfigFindLabel(LIS_config, "VODSVRFM parameter file:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,vodsvrfm_struc(n)%parameter_fname, rc=rc)
            call LIS_verify(rc, "VODSVRFM parameter file: not defined")
        enddo

        do n=1,LIS_rc%nnest
            allocate(vodsvrfm_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 1")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 2")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 3")
            call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 4")

            call add_sfc_fields(n,LIS_forwardState(n),"VODSVRFM_VOD")
        enddo



        ! read in parameter files
        do n=1, LIS_rc%nnest
            write(LIS_logunit,*) '[INFO] Reading ',&
                trim(vodsvrfm_struc(n)%parameter_fname)

            ios = nf90_open(path=trim(vodsvrfm_struc(n)%parameter_fname),&
                mode=NF90_NOWRITE,ncid=nid)
            call LIS_verify(ios,'Error opening file '&
                //trim(vodsvrfm_struc(n)%parameter_fname))

            ! check if ngrid is as expected
            ios = nf90_inq_dimid(nid, "ngrid", ngridId)
            call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
            ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
            call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
            if (ngrid /= LIS_rc%glbngrid_red(n)) then
                write(LIS_logunit, *) "[ERR] ngrid in "//trim(vodsvrfm_struc(n)%parameter_fname)&
                     //" not consistent with expected ngrid: ", ngrid,&
                     " instead of ",LIS_rc%glbngrid_red(n)
                call LIS_endrun
            endif

            ! read n_SV, because we don't know it's length a priori
            ios = nf90_inq_dimid(nid, "n_SV", n_SVId)
            call LIS_verify(ios, "Error nf90_inq_varid: n_SV")
            ios = nf90_inquire_dimension(nid, n_SVId, len=ntimes)
            call LIS_verify(ios, "Error nf90_inquire_dimension: n_SV")
            vodsvrfm_struc(n)%n_SV = n_SV

            ! now that we have n_SV, we can allocate the coefficient arrays
            ! for the nest
            m = LIS_rc%lsm_index
            allocate(vodsvrfm_struc(n)%intercept(LIS_rc%npatch(n, m)))
            allocate(vodsvrfm_struc(n)%actual_n_SV(LIS_rc%npatch(n, m)))
            allocate(vodsvrfm_struc(n)%dual_coef(LIS_rc%npatch(n, m), n_SV))
            allocate(vodsvrfm_struc(n)%support_vectors(LIS_rc%npatch(n, m), n_SV, 5))
            allocate(vodsvrfm_struc(n)%gam(LIS_rc%npatch(n, m), 5))

            ! read the arrays from the file
            call read_1d_coef_from_file(nid, ngrid, "intercept_", vodsvrfm_struc(n)%intercept)
            call read_1d_coef_from_file(nid, ngrid, "n_SV_", vodsvrfm_struc(n)%actual_n_SV)
            call read_2d_coef_from_file(nid, ngrid, n_SV, "dual_coef_", vodsvrfm_struc(n)%dual_coef)
            call read_2d_coef_from_file(nid, ngrid, 5, "gamma_", vodsvrfm_struc(n)%gam)
            call read_3d_coef_from_file(nid, ngrid, n_SV, 5, "support_vectors_",&
                 vodsvrfm_struc(n)%support_vectors)
        enddo


        write(LIS_logunit,*) '[INFO] Finished VODSVRFM setup'

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
    end subroutine VODSVRFM_initialize
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


    subroutine VODSVRFM_f2t(n)

        implicit none

        integer, intent(in)    :: n

    end subroutine VODSVRFM_f2t


    subroutine VODSVRFM_geometry(n)
        implicit none
        integer, intent(in)    :: n

    end subroutine VODSVRFM_geometry
    !Do nothing for now
    subroutine VODSVRFM_run(n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        integer, intent(in) :: n

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
        real                :: intercept, actual_n_SV
        real                :: dual_coef(:), gam(:)
        real,               :: support_vectors(:,:)
        real, pointer       :: vodval(:)
        logical             :: coefs_valid, values_valid
        real                :: kernelval

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 1", sm1)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 2", sm2)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 3", sm3)
        call getsfcvar(LIS_sfcState(n), "Soil Moisture Layer 4", sm4)


        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
            write(LIS_logunit, *) "[ERR] ntimes in "//trim(vodsvrfm_struc(n)%parameter_fname)&
                 //" must be 1 or 12, but is ", vodsvrfm_struc(n)%ntimes
            call LIS_endrun
        endif
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            intercept = vodsvrfm_struc(n)%intercept(t)
            actual_n_SV = vodsvrfm_struc(n)%actual_n_SV(t)
            dual_coef = vodsvrfm_struc(n)%dual_coef(t, :)
            gam = vodsvrfm_struc(n)%gam(t, :)
            support_vectors = vodsvrfm_struc(n)%support_vectors(t, :, :)

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
                vodsvrfm_struc(n)%VOD(t) = vod
            else
                vodsvrfm_struc(n)%VOD(t)=LIS_rc%udef
            endif

            if (vodsvrfm_struc(n)%VOD(t).ne.LIS_rc%udef.and.vodsvrfm_struc(n)%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vodsvrfm_struc(n)%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODSVRFM_VOD", vodval)
        vodval = vodsvrfm_struc(n)%VOD

    end subroutine VODSVRFM_run


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
        call LIS_verify(status, "Error in StateGet: VODSVRFM_getsfcvar "//trim(varname))
        call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
        call LIS_verify(status, "Error in FieldGet: VODSVRFM_getsfcvar "//trim(varname))

    end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: VODSVRFM_output
!!!! \label{VODSVRFM_output}
!!!!
!!!! !INTERFACE:
    subroutine VODSVRFM_output(n)
        integer, intent(in) :: n
    end subroutine VODSVRFM_output
#endif
end module VODSVRFM_Mod
