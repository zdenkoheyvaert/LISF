!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: VODCAlai_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   read VODCA X data that has been rescaled to a reference LAI.
! 
! !REVISION HISTORY: 
!  18 Nov 2021    Samuel Scherrer; initial reader based on VODCA LAI reader
! 
module VODCAlai_Mod
    ! !USES: 
    use ESMF
    use map_utils
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    PRIVATE

    !-----------------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    !-----------------------------------------------------------------------------
    public :: VODCAlai_setup
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: VODCAlai_struc
    !EOP
    type, public:: VODCAlai_dec

        character*100          :: version
        logical                :: startMode
        integer                :: nc
        integer                :: nr
        integer                :: mi
        real,     allocatable  :: laiobs1(:)
        real,     allocatable  :: laiobs2(:)
        real                   :: gridDesci(50)    
        real*8                 :: time1, time2
        integer                :: fnd
        character*20           :: nc_varname
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

    end type VODCAlai_dec

    type(VODCAlai_dec),allocatable :: VODCAlai_struc(:)

contains

    !BOP
    ! 
    ! !ROUTINE: VODCAlai_setup
    ! \label{VODCAlai_setup}
    ! 
    ! !INTERFACE: 
    subroutine VODCAlai_setup(k, OBS_State, OBS_Pert_State)
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
        !   creation of data structures required for handling VODCA LAI data.
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

        allocate(VODCAlai_struc(LIS_rc%nnest))

        call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ConfigFindLabel(LIS_config,"VODCA LAI data directory:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,laiobsdir,&
                 rc=status)
            call LIS_verify(status, 'VODCA LAI data directory: is missing')

            call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
                 laiobsdir, rc=status)
            call LIS_verify(status)
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"VODCA LAI netCDF variable name:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,VODCAlai_struc(n)%nc_varname,&
                 rc=status)
            call LIS_verify(status, 'VODCA LAI netCDF variable name: is missing')
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
             '[INFO] read VODCA LAI data specifications'       

        do n=1,LIS_rc%nnest

            allocate(VODCAlai_struc(n)%laiobs1(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
            allocate(VODCAlai_struc(n)%laiobs2(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

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

        do n=1,LIS_rc%nnest

            if(LIS_rc%lis_obs_map_proj(k).eq."latlon") then
                ! original spatial resolution of the downloaded data product
                VODCAlai_struc(n)%lat_lower_left = 89.875
                VODCAlai_struc(n)%lon_lower_left = -179.875
                VODCAlai_struc(n)%lat_upper_right = -89.875
                VODCAlai_struc(n)%lon_upper_right = 179.875
                ! dlat is positive since LIS will figure out that latitude is
                ! decreasing
                VODCAlai_struc(n)%dlat = 0.25
                VODCAlai_struc(n)%dlon = 0.25
                VODCAlai_struc(n)%nc = 720
                VODCAlai_struc(n)%nr = 1440
            elseif(LIS_rc%lis_obs_map_proj(k).eq."lambert") then
                write(unit=error_unit, fmt=*) &
                     'The VODCA_LAI module only works with latlon projection'
                stop 1
            endif

            ! from netCDF files

            VODCAlai_struc(n)%gridDesci(1) = 0  ! regular lat-lon grid
            VODCAlai_struc(n)%gridDesci(2) = VODCAlai_struc(n)%nc
            VODCAlai_struc(n)%gridDesci(3) = VODCAlai_struc(n)%nr 
            VODCAlai_struc(n)%gridDesci(4) = VODCAlai_struc(n)%lat_lower_left
            VODCAlai_struc(n)%gridDesci(5) = VODCAlai_struc(n)%lon_lower_left
            VODCAlai_struc(n)%gridDesci(6) = 128
            VODCAlai_struc(n)%gridDesci(7) = VODCAlai_struc(n)%lat_upper_right
            VODCAlai_struc(n)%gridDesci(8) = VODCAlai_struc(n)%lon_upper_right
            VODCAlai_struc(n)%gridDesci(9) = VODCAlai_struc(n)%dlat
            VODCAlai_struc(n)%gridDesci(10) = VODCAlai_struc(n)%dlon
            VODCAlai_struc(n)%gridDesci(20) = 64

            VODCAlai_struc(n)%mi = VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr

            !-----------------------------------------------------------------------------
            !   Use interpolation if LIS is running finer than native resolution. 
            !-----------------------------------------------------------------------------
            if(LIS_rc%obs_gridDesc(k,10).le.VODCAlai_struc(n)%dlon) then 

                allocate(VODCAlai_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(VODCAlai_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

                write(LIS_logunit,*)&
                     '[INFO] create interpolation input for VODCA LAI'       

                call bilinear_interp_input_withgrid(VODCAlai_struc(n)%gridDesci(:), &
                     LIS_rc%obs_gridDesc(k,:),&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                     VODCAlai_struc(n)%rlat, VODCAlai_struc(n)%rlon,&
                     VODCAlai_struc(n)%n11, VODCAlai_struc(n)%n12, &
                     VODCAlai_struc(n)%n21, VODCAlai_struc(n)%n22, &
                     VODCAlai_struc(n)%w11, VODCAlai_struc(n)%w12, &
                     VODCAlai_struc(n)%w21, VODCAlai_struc(n)%w22)
            else

                allocate(VODCAlai_struc(n)%n11(&
                     VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr))

                write(LIS_logunit,*)&
                     '[INFO] create upscaling input for VODCA LAI'       

                call upscaleByAveraging_input(VODCAlai_struc(n)%gridDesci(:),&
                     LIS_rc%obs_gridDesc(k,:),&
                     VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr, &
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), VODCAlai_struc(n)%n11)

                write(LIS_logunit,*)&
                     '[INFO] finished creating upscaling input for VODCA LAI'       
            endif

            call LIS_registerAlarm("VODCA LAI read alarm",&
                 86400.0, 86400.0)

            VODCAlai_struc(n)%startMode = .true. 

            call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
            call LIS_verify(status)

        enddo
    end subroutine VODCAlai_setup
end module VODCAlai_Mod
