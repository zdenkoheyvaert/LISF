!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: CustomLAI_Mod
!
! !DESCRIPTION:
!   This module contains a custom reader for global LAI netCDF products on a
!   regular latitude longitude grid stored in netCDF files, with all necessary
!   flags already applied beforehand during preprocessing.
!
!   The custom reader provides the following options for lis.config:
!
!   Custom LAI data directory:
!       Path to directory that contains image netCDF files.
!   Custom LAI netCDF variable name:
!       Name of the variable in the netCDF files.
!   Custom LAI netCDF name prefix:
!       Prefix for the netCDF image file name. The image file names must follow
!       the pattern <prefix>_YYYY_MM_DD.nc
!   Custom LAI spatial resolution
!       Spatial resolution of the product. It is assumed that the grid starts
!       at -180 + res/2 longitude and -90 + res/2 latitude.
!   Data assimilation scaling strategy:
!       Options are "none", "CDF matching", "seasonal", "seasonal multiplicative"
!   Custom LAI model scaling file:
!       Path to the model CDF/mean file (optional, only required if scaling is applied)
!   Custom LAI observation scaling file:
!       Path to the observation CDF/mean file (optional, only required if scaling is applied)
!   Custom LAI varname in model scaling file:
!       Prefix of the names of the variables in the model scaling files,
!       e.g. <prefix>_mu, <prefix>_sigma, ... (optional, only required if
!       scaling is applied)
!   Custom LAI varname in observation scaling file:
!       Prefix of the names of the variables in the observation scaling files,
!       e.g. <prefix>_mu, <prefix>_sigma, ... (optional, only required if
!       scaling is applied)
!   Custom LAI number of bins in the CDF: (optional, only required if CDF scaling is applied)
!       Number of bins in the CDFs.
!
! !REVISION HISTORY:
!  02 Mar 2022    Samuel Scherrer; initial reader based on MODIS LAI reader
!
module CustomLAI_Mod
    ! !USES:
    use CustomNcReader_Mod, only: CustomNcReader_dec

    implicit none

    public :: CustomLAI_setup
    public :: CustomLAI_struc

    type(CustomNcReader_dec),allocatable :: CustomLAI_struc(:)

contains

    subroutine CustomLAI_setup(k, OBS_State, OBS_Pert_State)

        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomLAI_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomLAI_struc(n)%varname = "LAI"
            CustomLAI_struc(n)%min_value = 0.0001
            CustomLAI_struc(n)%max_value = 10.0
            CustomLAI_struc(n)%qcmin_value = 0.0
            CustomLAI_struc(n)%qcmax_value = 100.0
        enddo

        call CustomNcReader_setup(k, OBS_State, OBS_Pert_State, CustomLAI_struc)

    end subroutine CustomLAI_setup

    subroutine read_CustomLAI(n, k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: read_CustomNetCDF
        integer, intent(in)       :: n, k
        type(ESMF_State)          :: OBS_State
        type(ESMF_State)          :: OBS_Pert_State
        call read_CustomNetCDF(n, k, OBS_State, OBS_Pert_State, CustomLAI_struc)
    end subroutine read_CustomLAI

end module CustomLAI_Mod
