!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: GenericLAI_Mod
! 
! !DESCRIPTION: 
!   This is a reader for generic global LAI products on a regular latitude
!   longitude grid stored in netCDF files, with all necessary flags already
!   applied beforehand during preprocessing.
!
!   The reader provides the following options for lis.config:
!
!   Generic LAI data directory:
!       Path to directory that contains image netCDF files.
!   Generic LAI netCDF variable name:
!       Name of the LAI variable in the netCDF files.
!   Generic LAI netCDF name prefix:
!       Prefix for the netCDF image file name. The image file names must follow
!       the pattern <prefix>_YYYY_MM_DD.nc
!   Generic LAI spatial resolution
!       Spatial resolution of the product. It is assumed that the grid starts
!       at -180 + res/2 longitude and -90 + res/2 latitude.
!
!   Data assimilation scaling strategy:
!       Options are "none", "CDF matching", "seasonal", "seasonal multiplicative"
!   Generic LAI model scaling file:
!       Path to the model LAI CDF/mean file (optional, only required if scaling is applied)
!   Generic LAI observation scaling file:
!       Path to the observation LAI CDF/mean file (optional, only required if scaling is applied)
!   Generic LAI number of bins in the CDF: (optional, only required if CDF scaling is applied)
!       Number of bins in the CDFs.
!     
! 
! !REVISION HISTORY: 
!  02 Mar 2022    Samuel Scherrer; initial reader based on MODIS LAI reader
! 
module GenericLAI_Mod
    ! !USES: 
    use ESMF
    use map_utils
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    PRIVATE

    !-----------------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    !-----------------------------------------------------------------------------
    public :: GenericLAI_setup, GenericLAI_rescale_with_seasonal_scaling,&
         GenericLAI_updateSsdev
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: GenericLAI_struc
    !EOP
    type, public:: GenericLAI_dec

        character*100          :: version
        logical                :: startMode
        integer                :: nc
        integer                :: nr
        integer                :: mi
        real                   :: gridDesci(50)    
        real*8                 :: time1, time2
        integer                :: fnd
        integer                :: useSsdevScal
        logical                :: mult_scaling
        character*20           :: nc_varname
        character*100          :: nc_prefix
        real*8                 :: spatialres
        real*8                 :: dlat, dlon
        real*8                 :: lat_lower_left, lon_lower_left
        real*8                 :: lat_upper_right, lon_upper_right
        real,    allocatable :: rlat(:)
        real,    allocatable :: rlon(:)
        integer, allocatable :: n11(:)
        integer, allocatable :: n12(:)
        integer, allocatable :: n21(:)
        integer, allocatable :: n22(:)
        real,    allocatable :: w11(:)
        real,    allocatable :: w12(:)
        real,    allocatable :: w21(:)
        real,    allocatable :: w22(:)

        real                       :: ssdev_inp
        real,    allocatable       :: model_xrange(:,:,:)
        real,    allocatable       :: obs_xrange(:,:,:)
        real,    allocatable       :: model_cdf(:,:,:)
        real,    allocatable       :: obs_cdf(:,:,:)
        real,    allocatable       :: model_mu(:,:)
        real,    allocatable       :: obs_mu(:,:)
        real,    allocatable       :: model_sigma(:,:)
        real,    allocatable       :: obs_sigma(:,:)

        integer                :: nbins
        integer                :: ntimes

    end type GenericLAI_dec

    type(GenericLAI_dec),allocatable :: GenericLAI_struc(:)

contains

    !BOP
    ! 
    ! !ROUTINE: GenericLAI_setup
    ! \label{GenericLAI_setup}
    ! 
    ! !INTERFACE: 
    subroutine GenericLAI_setup(k, OBS_State, OBS_Pert_State)
        ! !USES: 
        use ESMF
        use LIS_coreMod
        use LIS_timeMgrMod
        use LIS_historyMod
        use LIS_dataAssimMod
        use LIS_perturbMod
        use LIS_DAobservationsMod
        use LIS_logmod

        implicit none 

        ! !ARGUMENTS: 
        integer                ::  k
        type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
        type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
        ! 
        ! !DESCRIPTION: 
        !   
        !   This routine completes the runtime initializations and 
        !   creation of data structures required for handling Generic LAI data.
        !  
        !   The arguments are: 
        !   \begin{description}
        !    \item[k] number of observation state 
        !    \item[OBS\_State]   observation state 
        !    \item[OBS\_Pert\_State] observation perturbations state
        !   \end{description}
        !EOP
        integer                ::  n,i
        integer                ::  ftn
        integer                ::  status
        type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
        type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
        type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
        type(ESMF_ArraySpec)   ::  pertArrSpec
        character*100          ::  laiobsdir
        character*100          ::  temp
        real, parameter        ::  minssdev =0.001
        real,  allocatable     ::  ssdev(:)
        character*1            ::  vid(2)
        type(pert_dec_type)    ::  obs_pert
        real, pointer          ::  obs_temp(:,:)
        character*40, allocatable  ::  vname(:)
        real        , allocatable  ::  varmin(:)
        real        , allocatable  ::  varmax(:)
        character*100          :: modelscalingfile(LIS_rc%nnest)
        character*100          :: obsscalingfile(LIS_rc%nnest)
        integer                :: c,r
        integer                :: ngrid
        integer                :: timeidx

        allocate(GenericLAI_struc(LIS_rc%nnest))

        call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        !------------------------------------------------------------
        ! File options
        !------------------------------------------------------------
        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI data directory:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,laiobsdir,&
                 rc=status)
            call LIS_verify(status, 'Generic LAI data directory: is missing')

            call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
                 laiobsdir, rc=status)
            call LIS_verify(status)
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI netCDF variable name:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,GenericLAI_struc(n)%nc_varname,&
                 rc=status)
            call LIS_verify(status, 'Generic LAI netCDF variable name: is missing')
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI netCDF name prefix:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,GenericLAI_struc(n)%nc_prefix,&
                 rc=status)
            call LIS_verify(status, 'Generic LAI netCDF name prefix: is missing')
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI spatial resolution:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,Genericlai_struc(n)%spatialres,&
                 rc=status)
            call LIS_verify(status, 'Generic LAI spatial resolution: is missing')
        enddo

        !------------------------------------------------------------
        ! Options for scaling
        !------------------------------------------------------------
        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI use scaled standard deviation model:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 
                call ESMF_ConfigGetAttribute(LIS_config,GenericLAI_struc(n)%useSsdevScal, &
                     rc=status)
                call LIS_verify(status, "Generic LAI use scaled standard deviation model: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI model scaling file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 
                call ESMF_ConfigGetAttribute(LIS_config,modelscalingfile(n),rc=status)
                call LIS_verify(status, 'Generic LAI model scaling file: not defined')
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI observation scaling file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 
                call ESMF_ConfigGetAttribute(LIS_config,obsscalingfile(n),rc=status)
                call LIS_verify(status, 'Generic LAI observation scaling file: not defined')
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config, "Generic LAI number of bins in the CDF:", rc=status)
        do n=1, LIS_rc%nnest
            if(LIS_rc%dascaloption(k).eq."CDF matching") then 
                call ESMF_ConfigGetAttribute(LIS_config,Genericlai_struc(n)%nbins, rc=status)
                call LIS_verify(status, "Generic LAI number of bins in the CDF: not defined")
            endif
        enddo


        do n=1,LIS_rc%nnest
            call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
                 .false., rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
                 -99.0, rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
                 .false., rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
                 LIS_rc%obs_ngrid(k),rc=status)
            call LIS_verify(status)

        enddo

        write(LIS_logunit,*)&
             '[INFO] read Generic LAI data specifications'       

        !------------------------------------------------------------
        ! Read perturbation settings
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest

            write(unit=temp,fmt='(i2.2)') 1
            read(unit=temp,fmt='(2a1)') vid

            obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
                 grid=LIS_obsVecGrid(n,k),&
                 name="Observation"//vid(1)//vid(2), rc=status)
            call LIS_verify(status)

            !Perturbations State
            write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
                 trim(LIS_rc%obsattribfile(k))
            ftn = LIS_getNextUnitNumber()
            open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
            read(ftn,*)
            read(ftn,*) LIS_rc%nobtypes(k)
            read(ftn,*)

            allocate(vname(LIS_rc%nobtypes(k)))
            allocate(varmax(LIS_rc%nobtypes(k)))
            allocate(varmin(LIS_rc%nobtypes(k)))

            do i=1,LIS_rc%nobtypes(k)
                read(ftn,fmt='(a40)') vname(i)
                read(ftn,*) varmin(i),varmax(i)
                write(LIS_logunit,*) '[INFO] ',vname(i),varmin(i),varmax(i)
            enddo
            call LIS_releaseUnitNumber(ftn)   

            allocate(ssdev(LIS_rc%obs_ngrid(k)))

            if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
                allocate(obs_pert%vname(1))
                allocate(obs_pert%perttype(1))
                allocate(obs_pert%ssdev(1))
                allocate(obs_pert%stdmax(1))
                allocate(obs_pert%zeromean(1))
                allocate(obs_pert%tcorr(1))
                allocate(obs_pert%xcorr(1))
                allocate(obs_pert%ycorr(1))
                allocate(obs_pert%ccorr(1,1))

                call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
                     obs_pert)

                ! Set obs err to be uniform (will be rescaled later for each grid point). 
                ssdev = obs_pert%ssdev(1)

                pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                     grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
                     rc=status)
                call LIS_verify(status)

                ! initializing the perturbations to be zero 
                call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
                call LIS_verify(status)
                obs_temp(:,:) = 0 

                call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
                     obs_pert%perttype(1), rc=status)
                call LIS_verify(status)

                if(LIS_rc%obs_ngrid(k).gt.0) then 
                    call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                         ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                    call LIS_verify(status)
                endif

                call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
                     obs_pert%stdmax(1), rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
                     obs_pert%zeromean(1),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
                     obs_pert%tcorr(1),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
                     obs_pert%xcorr(1),rc=status)

                call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
                     obs_pert%ycorr(1),rc=status)

                call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
                     obs_pert%ccorr(1,:),itemCount=1,rc=status)

                call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
                call LIS_verify(status)         

            endif

            deallocate(vname)
            deallocate(varmax)
            deallocate(varmin)
            deallocate(ssdev)   

        enddo

        !------------------------------------------------------------
        ! Initialize scaling
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest

            if(LIS_rc%dascaloption(k).ne."none") then

                allocate(ssdev(LIS_rc%obs_ngrid(k)))
                ssdev = obs_pert%ssdev(1)

                allocate(GenericLAI_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
                     GenericLAI_struc(n)%ntimes))
                allocate(GenericLAI_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
                     GenericLAI_struc(n)%ntimes))
                allocate(GenericLAI_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
                     GenericLAI_struc(n)%ntimes))
                allocate(GenericLAI_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
                     GenericLAI_struc(n)%ntimes))

                if(LIS_rc%dascaloption(k).eq."CDF matching") then 
                    !------------------------------------------------------------
                    ! CDF matching
                    !------------------------------------------------------------

                    call LIS_getCDFattributes(k,modelscalingfile(n),&
                         GenericLAI_struc(n)%ntimes, ngrid)

                    if(GenericLAI_struc(n)%ntimes.eq.1) then 
                        timeidx = 1
                    else
                        timeidx = LIS_rc%mo
                    endif

                    allocate(GenericLAI_struc(n)%model_xrange(&
                         LIS_rc%obs_ngrid(k), GenericLAI_struc(n)%ntimes, &
                         GenericLAI_struc(n)%nbins))
                    allocate(GenericLAI_struc(n)%obs_xrange(&
                         LIS_rc%obs_ngrid(k), GenericLAI_struc(n)%ntimes, &
                         GenericLAI_struc(n)%nbins))
                    allocate(GenericLAI_struc(n)%model_cdf(&
                         LIS_rc%obs_ngrid(k), GenericLAI_struc(n)%ntimes, &
                         GenericLAI_struc(n)%nbins))
                    allocate(GenericLAI_struc(n)%obs_cdf(&
                         LIS_rc%obs_ngrid(k), GenericLAI_struc(n)%ntimes, & 
                         GenericLAI_struc(n)%nbins))

                    !----------------------------------------------------------------------------
                    ! Read the model and observation CDF data
                    !----------------------------------------------------------------------------

                    ! model mean sigma
                    call LIS_readMeanSigmaData(n,k,&
                         GenericLAI_struc(n)%ntimes, & 
                         LIS_rc%obs_ngrid(k), &
                         modelscalingfile(n), &
                         "LAI",&
                         GenericLAI_struc(n)%model_mu,&
                         GenericLAI_struc(n)%model_sigma)

                    ! observation mean sigma
                    call LIS_readMeanSigmaData(n,k,&
                         GenericLAI_struc(n)%ntimes, & 
                         LIS_rc%obs_ngrid(k), &
                         obsscalingfile(n), &
                         "LAI",&
                         GenericLAI_struc(n)%obs_mu,&
                         GenericLAI_struc(n)%obs_sigma)

                    ! model CDF
                    call LIS_readCDFdata(n,k,&
                         GenericLAI_struc(n)%nbins,&
                         GenericLAI_struc(n)%ntimes, & 
                         LIS_rc%obs_ngrid(k), &
                         modelscalingfile(n), &
                         "LAI",&
                         GenericLAI_struc(n)%model_xrange,&
                         GenericLAI_struc(n)%model_cdf)

                    ! observation CDF
                    call LIS_readCDFdata(n,k,&
                         GenericLAI_struc(n)%nbins,&
                         GenericLAI_struc(n)%ntimes, & 
                         LIS_rc%obs_ngrid(k), &
                         obsscalingfile(n), &
                         "LAI",&
                         GenericLAI_struc(n)%obs_xrange,&
                         GenericLAI_struc(n)%obs_cdf)

                    !------------------------------------------------------------
                    ! end CDF matching
                    !------------------------------------------------------------

                elseif (LIS_rc%dascaloption(k).eq."seasonal"&
                     .or.LIS_rc%dascaloption(k).eq."seasonal multiplicative") then

                    !------------------------------------------------------------
                    ! seasonal scaling
                    !------------------------------------------------------------

                    GenericLAI_struc(n)%mult_scaling = &
                         LIS_rc%dascaloption(k).eq."seasonal multiplicative"

                    GenericLAI_struc(n)%ntimes = 366
                    timeidx = LIS_rc%da

                    call GenericLAI_readSeasonalScalingData(n,k,&
                         GenericLAI_struc(n)%ntimes, & 
                         LIS_rc%obs_ngrid(k), &
                         modelscalingfile(n), &
                         "LAI",&
                         GenericLAI_struc(n)%model_mu,&
                         GenericLAI_struc(n)%model_sigma)

                    if (GenericLAI_struc(n)%mult_scaling) then
                        ! When doing the multiplicative scaling, the standard
                        ! deviation in the file is not used, and the
                        ! perturbation standard deviations should be corrected
                        ! by the ratio of means instead of the ratio of sigmas.
                        GenericLAI_struc(n)%model_sigma = GenericLAI_struc(n)%model_mu
                        GenericLAI_struc(n)%obs_sigma = GenericLAI_struc(n)%obs_mu
                    endif

                endif

                call GenericLAI_updateSsdev(k,&
                     GenericLAI_struc(n)%obs_sigma(:, timeidx),&
                     GenericLAI_struc(n)%model_sigma(:, timeidx),&
                     ssdev)

                if(LIS_rc%obs_ngrid(k).gt.0) then 
                    call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                         ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                    call LIS_verify(status, &
                         "updating perturbation standard deviation failed")
                endif

                deallocate(ssdev)

            endif
        enddo

        !------------------------------------------------------------
        ! Read grid information and resample
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest

            if(LIS_rc%lis_obs_map_proj(k).eq."latlon") then
                GenericLAI_struc(n)%lat_lower_left = 90 - 0.5 * GenericLAI_struc(n)%spatialres
                GenericLAI_struc(n)%lat_upper_right = -90 + 0.5 * GenericLAI_struc(n)%spatialres
                GenericLAI_struc(n)%lon_lower_left = -180 + 0.5 * GenericLAI_struc(n)%spatialres
                GenericLAI_struc(n)%lon_upper_right = 180 - 0.5 * GenericLAI_struc(n)%spatialres
                ! dlat is positive since LIS will figure out that latitude is decreasing
                GenericLAI_struc(n)%dlat = GenericLAI_struc(n)%spatialres
                GenericLAI_struc(n)%dlon = GenericLAI_struc(n)%spatialres
                GenericLAI_struc(n)%nr = nint(180.0 / GenericLAI_struc(n)%spatialres)
                GenericLAI_struc(n)%nc = 2 * GenericLAI_struc(n)%nr
            elseif(LIS_rc%lis_obs_map_proj(k).eq."lambert") then
                write(unit=error_unit, fmt=*) &
                     'The Generic_LAI module only works with latlon projection'
                stop 1
            endif

            GenericLAI_struc(n)%gridDesci(1) = 0  ! regular lat-lon grid
            GenericLAI_struc(n)%gridDesci(2) = GenericLAI_struc(n)%nc
            GenericLAI_struc(n)%gridDesci(3) = GenericLAI_struc(n)%nr 
            GenericLAI_struc(n)%gridDesci(4) = GenericLAI_struc(n)%lat_lower_left
            GenericLAI_struc(n)%gridDesci(5) = GenericLAI_struc(n)%lon_lower_left
            GenericLAI_struc(n)%gridDesci(6) = 128
            GenericLAI_struc(n)%gridDesci(7) = GenericLAI_struc(n)%lat_upper_right
            GenericLAI_struc(n)%gridDesci(8) = GenericLAI_struc(n)%lon_upper_right
            GenericLAI_struc(n)%gridDesci(9) = GenericLAI_struc(n)%dlat
            GenericLAI_struc(n)%gridDesci(10) = GenericLAI_struc(n)%dlon
            GenericLAI_struc(n)%gridDesci(20) = 64

            GenericLAI_struc(n)%mi = GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr

            !-----------------------------------------------------------------------------
            !   Use interpolation if LIS is running finer than native resolution. 
            !-----------------------------------------------------------------------------
            if(LIS_rc%obs_gridDesc(k,10).lt.GenericLAI_struc(n)%dlon) then

                allocate(GenericLAI_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(GenericLAI_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

                write(LIS_logunit,*)&
                     '[INFO] create interpolation input for Generic LAI'       

                call bilinear_interp_input_withgrid(GenericLAI_struc(n)%gridDesci(:), &
                     LIS_rc%obs_gridDesc(k,:),&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                     GenericLAI_struc(n)%rlat, GenericLAI_struc(n)%rlon,&
                     GenericLAI_struc(n)%n11, GenericLAI_struc(n)%n12, &
                     GenericLAI_struc(n)%n21, GenericLAI_struc(n)%n22, &
                     GenericLAI_struc(n)%w11, GenericLAI_struc(n)%w12, &
                     GenericLAI_struc(n)%w21, GenericLAI_struc(n)%w22)
            else

                allocate(GenericLAI_struc(n)%n11(&
                     GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr))

                write(LIS_logunit,*)&
                     '[INFO] create upscaling input for Generic LAI'       

                call upscaleByAveraging_input(GenericLAI_struc(n)%gridDesci(:),&
                     LIS_rc%obs_gridDesc(k,:),&
                     GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr, &
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), GenericLAI_struc(n)%n11)

                write(LIS_logunit,*)&
                     '[INFO] finished creating upscaling input for Generic LAI'       
            endif

            call LIS_registerAlarm("Generic LAI read alarm",&
                 86400.0, 86400.0)

            GenericLAI_struc(n)%startMode = .true. 

            call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
            call LIS_verify(status, "Adding observation field failed")

        enddo
    end subroutine GenericLAI_setup

    !BOP
    ! 
    ! !ROUTINE: GenericLAI_rescale_with_seasonal_scaling
    ! \label{GenericLAI_rescale_with_seasonal_scaling}
    !
    ! !INTERFACE:
    subroutine GenericLAI_rescale_with_seasonal_scaling(&
         n,             & 
         k,             & 
         timeidx,        & 
         ntimes,         & 
         max_obs_value, &
         min_obs_value, &
         multiplicative, &
         model_mu,    &       
         model_sigma,    &       
         obs_mu,       &
         obs_sigma,       &
         obs_value) 

        use LIS_coreMod, only: LIS_rc
        use LIS_DAobservationsMod, only: LIS_obs_domain

        implicit none
        ! 
        ! !ARGUMENTS: 
        integer, intent(in)      :: n 
        integer, intent(in)      :: k
        integer, intent(in)      :: ntimes
        integer, intent(in)      :: timeidx
        real, intent(in)         :: max_obs_value
        real, intent(in)         :: min_obs_value
        logical, intent(in)      :: multiplicative
        real, intent(in)         :: model_mu(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(in)         :: model_sigma(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(in)         :: obs_mu(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(in)         :: obs_sigma(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(inout)      :: obs_value(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
        !
        ! !DESCRIPTION: 
        ! 
        !   This routine rescales the input observation data by subtracting the
        !   observation mean for the given time index and adding the model mean
        !   for the same time index

        !  The arguments are: 
        !  \begin{description}
        !  \item[n]               index of the nest
        !  \item[k]               index of observation state
        !  \item[timeidx]         time index for mean value to use, .e.g. month
        !                         or day
        !  \item[ntimes]          number of times in the model and obs mean files
        !  \item[max\_obs\_value] maximum allowable value of observation
        !  \item[min\_obs\_value] minimum allowable value of observation
        !  \item[model\_mu]       mean values of the model OL run
        !  \item[model\_sigma]    std.dev. values of the model OL run
        !  \item[obs\_mu]         mean values of the observations
        !  \item[obs\_sigma]      std.dev. values of the observations
        !  \item[obs\_value]      observation value to be rescaled. 
        ! \end{description}
        !EOP

        integer             :: grididx,i
        integer             :: binval
        integer             :: col,row
        real                :: cdf_obsval
        real                :: obs_tmp

        do grididx=1,LIS_rc%obs_ngrid(k)

            col = LIS_obs_domain(n,k)%col(grididx)
            row = LIS_obs_domain(n,k)%row(grididx)

            if(obs_value(col,row).ne.-9999.0 &
                 .and. obs_mu(grididx, timeidx).ne.-9999.0 &
                 .and. obs_sigma(grididx, timeidx).gt.epsilon(0.0)) then 

                if (multiplicative.and.obs_mu(grididx, timeidx).ne.epsilon(0.0)) then
                    obs_tmp = obs_value(col, row) / obs_mu(grididx, timeidx)&
                         * model_mu(grididx, timeidx)
                else if ((.not.multiplicative)&
                        .and. obs_sigma(grididx, timeidx).gt.(epsilon(0.0))) then
                    obs_tmp = (obs_value(col,row) - obs_mu(grididx,timeidx)) &
                         / obs_sigma(grididx,timeidx)
                    obs_tmp = obs_tmp * model_sigma(grididx,timeidx)&
                         + model_mu(grididx, timeidx)
                else
                    obs_tmp = LIS_rc%udef
                endif

                if (obs_tmp < min_obs_value .or. obs_tmp > max_obs_value) then
                    obs_tmp = LIS_rc%udef
                endif

                obs_value(col, row) = obs_tmp
            else
                obs_value(col,row) = LIS_rc%udef
            endif
        enddo
    end subroutine GenericLAI_rescale_with_seasonal_scaling

    !BOP
    ! !ROUTINE: GenericLAI_readSeasonalScalingData
    ! \label{GenericLAI_readSeasonalScalingData}
    !
    ! !INTERFACE: 
    subroutine GenericLAI_readSeasonalScalingData(&
         n, k, ntimes, ngrid, filename, varname, mu, sigma)

        use netcdf
        use LIS_coreMod, only: LIS_rc
        use LIS_logMod, only: LIS_logunit, LIS_verify
        use LIS_DAobservationsMod, only: LIS_convertObsVarToLocalSpace

        implicit none
        ! !ARGUMENTS:      
        integer,   intent(in)    :: n
        integer,   intent(in)    :: k
        integer,   intent(in)    :: ntimes
        integer,   intent(in)    :: ngrid
        character(len=*)         :: filename
        character(len=*)         :: varname
        real                     :: mu(ngrid, ntimes)
        real                     :: sigma(ngrid, ntimes)
        ! 
        ! !DESCRIPTION: 
        !  This routine reads the input seasonal mean file.
        ! 
        !  The arguments are: 
        !  \begin{description}
        !  \item[n]             index of the nest
        !  \item[k]             index of observation state
        !  \item[filename]      name of the CDF file
        !  \item[varname]       name of the variable being extracted.
        !  \item[mu]            mean values
        ! \end{description}
        !EOP
        integer                  :: j
        integer                  :: nlevsId, gId
        integer                  :: ngrid_file, nlevs_file
        integer                  :: muid, sigmaid
        real, allocatable        :: mu_file(:,:,:), sigma_file(:,:,:)
        integer                  :: nid

        write(LIS_logunit,*) '[INFO] Reading mean from seasonal seasonal file ',trim(filename)
        call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
             ncid=nid),'failed to open file '//trim(filename))

        call LIS_verify(nf90_inq_dimid(nid,trim(varname)//"_levels",nlevsId), &
             'nf90_inq_dimid failed for '//trim(varname)//"_levels")

        call LIS_verify(nf90_inquire_dimension(nid,nlevsId, len=nlevs_file),&
             'nf90_inquire_dimension failed for nlevsId')

        call LIS_verify(nf90_inq_dimid(nid, 'ngrid',gId), &
             'Error nf90_inq_dimid: ngrid')

        call LIS_verify(nf90_inquire_dimension(nid, gId, len=ngrid_file), &
             'Error nf90_inquire_dimension:ngrid') 

        allocate(mu_file(ngrid_file,ntimes,nlevs_file))
        allocate(sigma_file(ngrid_file,ntimes,nlevs_file))

        call LIS_verify(nf90_inq_varid(nid,trim(varname)//'_mu',muid),&
             'nf90_inq_varid failed for for '//trim(varname)//'_mu')

        call LIS_verify(nf90_get_var(nid,muid,mu_file),&
             'nf90_get_var failed for '//trim(varname)//'_mu')

        call LIS_verify(nf90_inq_varid(nid,trim(varname)//'_sigma',sigmaid),&
             'nf90_inq_varid failed for for '//trim(varname)//'_sigma')

        call LIS_verify(nf90_get_var(nid,sigmaid,sigma_file),&
             'nf90_get_var failed for '//trim(varname)//'_sigma')

        if(LIS_rc%obs_ngrid(k).gt.0) then 
            do j=1,ntimes
                call LIS_convertObsVarToLocalSpace(n,k,mu_file(:,j,1), mu(:,j))
                call LIS_convertObsVarToLocalSpace(n,k,sigma_file(:,j,1), sigma(:,j))
            enddo
        endif


        deallocate(mu_file)
        deallocate(sigma_file)

        call LIS_verify(nf90_close(nid),&
             'failed to close file '//trim(filename))
        write(LIS_logunit,*)&
             '[INFO] Successfully read mean and seasonal scaling file ',&
             trim(filename)
    end subroutine GenericLAI_readSeasonalScalingData

    subroutine GenericLAI_updateSsdev(k, obs_sigma, model_sigma, ssdev)
        use LIS_coreMod
        implicit none

        integer, intent(in)     :: k
        real, intent(in)        :: obs_sigma(:), model_sigma(:)
        real, intent(inout)     :: ssdev(:)

        integer                 :: grididx
        real, parameter         ::  minssdev =0.001


        do grididx=1,LIS_rc%obs_ngrid(k)
            if(obs_sigma(grididx).ne.LIS_rc%udef) then 
                if(obs_sigma(grididx).ne.0) then 
                    ssdev(grididx) = ssdev(grididx)&
                         * model_sigma(grididx)&
                         / obs_sigma(grididx)
                endif

                if(ssdev(grididx).lt.minssdev) then 
                    ssdev(grididx) = minssdev
                endif
            endif
        enddo

    end subroutine GenericLAI_updateSsdev


end module GenericLAI_Mod
