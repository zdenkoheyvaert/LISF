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
!   Generic LAI model CDF file: (optional, only required if scaling is applied)
!       Path to the model LAI CDF file
!   Generic LAI observation CDF file: (optional, only required if scaling is applied)
!       Path to the observation LAI CDF file
!   Generic LAI number of bins in the CDF: (optional, only required if scaling is applied)
!       Number of bins in the CDFs.
!   Generic LAI spatial resolution:
!       Spatial resolution of the product. It is assumed that the grid starts
!       at -180 + res/2 longitude and -90 + res/2 latitude.
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
    public :: GenericLAI_setup
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
        integer                ::  n,i,t,kk,jj
        integer                ::  ftn
        integer                ::  status
        type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
        type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
        type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
        type(ESMF_ArraySpec)   ::  pertArrSpec
        character*100          ::  laiobsdir
        character*100          ::  temp
        real,  allocatable         ::  ssdev(:)
        character*1            ::  vid(2)
        type(pert_dec_type)    ::  obs_pert
        real, pointer          ::  obs_temp(:,:)
        character*40, allocatable  ::  vname(:)
        real        , allocatable  ::  varmin(:)
        real        , allocatable  ::  varmax(:)
        integer                :: c,r

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
            if (Genericlai_struc(n)%isresampled) then
                call ESMF_ConfigGetAttribute(LIS_config,Genericlai_struc(n)%spatialres,&
                     rc=status)
                call LIS_verify(status, 'Generic LAI spatial resolution: is missing')
            endif
        enddo

        !------------------------------------------------------------
        ! Options for CDF matching
        !------------------------------------------------------------
        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI model CDF file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 
                call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
                call LIS_verify(status, 'Generic LAI model CDF file: not defined')
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Generic LAI observation CDF file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 
                call ESMF_ConfigGetAttribute(LIS_config,obscdffile(n),rc=status)
                call LIS_verify(status, 'Generic LAI observation CDF file: not defined')
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config, "Generic LAI number of bins in the CDF:", rc=status)
        do n=1, LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 
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
        ! Initialize CDF scaling
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then 

                call LIS_getCDFattributes(k,modelcdffile(n),&
                     GenericLAI_struc(n)%ntimes, ngrid)

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
                call LIS_readMeanSigmaData(n,k,&
                     GenericLAI_struc(n)%ntimes, & 
                     LIS_rc%obs_ngrid(k), &
                     modelcdffile(n), &
                     "LAI",&
                     GenericLAI_struc(n)%model_mu,&
                     GenericLAI_struc(n)%model_sigma)

                call LIS_readMeanSigmaData(n,k,&
                     GenericLAI_struc(n)%ntimes, & 
                     LIS_rc%obs_ngrid(k), &
                     obscdffile(n), &
                     "LAI",&
                     GenericLAI_struc(n)%obs_mu,&
                     GenericLAI_struc(n)%obs_sigma)

                call LIS_readCDFdata(n,k,&
                     GenericLAI_struc(n)%nbins,&
                     GenericLAI_struc(n)%ntimes, & 
                     LIS_rc%obs_ngrid(k), &
                     modelcdffile(n), &
                     "LAI",&
                     GenericLAI_struc(n)%model_xrange,&
                     GenericLAI_struc(n)%model_cdf)

                call LIS_readCDFdata(n,k,&
                     GenericLAI_struc(n)%nbins,&
                     GenericLAI_struc(n)%ntimes, & 
                     LIS_rc%obs_ngrid(k), &
                     obscdffile(n), &
                     "LAI",&
                     GenericLAI_struc(n)%obs_xrange,&
                     GenericLAI_struc(n)%obs_cdf)

                if(GenericLAI_struc(n)%useSsdevScal.eq.1) then 
                    if(GenericLAI_struc(n)%ntimes.eq.1) then 
                        jj = 1
                    else
                        jj = LIS_rc%mo
                    endif
                    do t=1,LIS_rc%obs_ngrid(k)
                        if(GenericLAI_struc(n)%obs_sigma(t,jj).ne.LIS_rc%udef) then 
                            print*, ssdev(t),GenericLAI_struc(n)%model_sigma(t,jj),&
                                 GenericLAI_struc(n)%obs_sigma(t,jj)
                            if(GenericLAI_struc(n)%obs_sigma(t,jj).ne.0) then 
                                ssdev(t) = ssdev(t)*GenericLAI_struc(n)%model_sigma(t,jj)/&
                                     GenericLAI_struc(n)%obs_sigma(t,jj)
                            endif

                            if(ssdev(t).lt.minssdev) then 
                                ssdev(t) = minssdev
                            endif
                        endif
                    enddo
                endif

                if(LIS_rc%obs_ngrid(k).gt.0) then 
                    call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                         ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                    call LIS_verify(status)
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
            call LIS_verify(status)

        enddo
    end subroutine GenericLAI_setup
end module GenericLAI_Mod
