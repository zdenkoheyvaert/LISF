!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_GenericLAI
! \label{read_GenericLAI}
!
! !REVISION HISTORY:
!  02 Mar 2022    Samuel Scherrer; initial reader based on MODIS LAI reader
!
! !INTERFACE: 
subroutine read_GenericLAI(n, k, OBS_State, OBS_Pert_State)
    ! !USES: 
    use ESMF
    use LIS_mpiMod
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod
    use LIS_dataAssimMod
    use LIS_DAobservationsMod
    use map_utils
    use LIS_pluginIndices
    use GenericLAI_Mod, only : GenericLAI_struc,&
                               GenericLAI_rescale_with_seasonal_scaling

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: k
    type(ESMF_State)    :: OBS_State
    type(ESMF_State)    :: OBS_Pert_State
    !
    ! !DESCRIPTION:
    !  
    !  reads the Generic LAI observations from NETCDF files.

    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[n] index of the nest
    !  \item[k] number of observation state
    !  \item[OBS\_State] observations state
    !  \item[OBS\_Pert\_State] observation perturbations state
    !  \end{description}
    !
    !EOP
    real,  parameter       :: MAX_LAI_VALUE=10.0, MIN_LAI_VALUE=0.0001
    integer                :: status
    character*100          :: laiobsdir
    character*300          :: fname
    integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
    real                   :: wt1, wt2,ts
    integer                :: count
    real                   :: cgmt
    real*8                 :: time
    logical                :: alarmCheck, file_exists,dataCheck
    integer                :: c,r,i,j,p,t
    real,          pointer :: obsl(:)
    type(ESMF_Field)       :: laifield
    integer                :: gid(LIS_rc%obs_ngrid(k))
    integer                :: assimflag(LIS_rc%obs_ngrid(k))
    logical                :: data_update
    logical                :: data_upd_flag(LIS_npes)
    logical                :: data_upd_flag_local
    logical                :: data_upd
    real                   :: laiobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                :: fnd
    real                   :: timenow
    real                   :: ssdev(LIS_rc%obs_ngrid(k))
    type(ESMF_Field)       :: pertfield
    integer                :: timeidx


    call ESMF_AttributeGet(OBS_State,"Data Directory",&
         laiobsdir, rc=status)
    call LIS_verify(status)
    call ESMF_AttributeGet(OBS_State,"Data Update Status",&
         data_update, rc=status)
    call LIS_verify(status)

    data_upd = .false. 

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Generic LAI read alarm")

    if(alarmCheck.or.GenericLAI_struc(n)%startMode) then 
        GenericLAI_struc(n)%startMode = .false.

        call create_GenericLAI_filename(GenericLAI_struc(n)%nc_prefix,&
             laiobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)

        inquire(file=fname,exist=file_exists)          
        if(file_exists) then 
            write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
            call read_Generic_LAI_data(n,k, fname,laiobs)
            fnd = 1
        else
            fnd = 0 
            write(LIS_logunit,*) '[WARN] Missing LAI file: ',trim(fname)
        endif
    else
        fnd = 0 
        laiobs = LIS_rc%udef
    endif

    dataCheck = .false.
    if(alarmCheck) then 
        if(fnd.eq.1) then 
            dataCheck = .true. 
        endif
    else
        fnd = 0 
        dataCheck = .false.
    endif

    if(dataCheck) then 

        call ESMF_StateGet(OBS_State,"Observation01",laifield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        call ESMF_FieldGet(laifield,localDE=0,farrayPtr=obsl,rc=status)
        call LIS_verify(status, 'Error: FieldGet')

        obsl = LIS_rc%udef 
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                         laiobs(c+(r-1)*LIS_rc%obs_lnc(k))
                endif
            enddo
        enddo

        !-------------------------------------------------------------------------
        !  Transform data to the LSM climatology using a CDF-scaling approach
        !-------------------------------------------------------------------------     
        
        if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then

            call LIS_rescale_with_CDF_matching(     &
                 n,k,                               & 
                 GenericLAI_struc(n)%nbins,         & 
                 GenericLAI_struc(n)%ntimes,        & 
                 MAX_LAI_VALUE,                      & 
                 MIN_LAI_VALUE,                      & 
                 GenericLAI_struc(n)%model_xrange,  &
                 GenericLAI_struc(n)%obs_xrange,    &
                 GenericLAI_struc(n)%model_cdf,     &
                 GenericLAI_struc(n)%obs_cdf,       &
                 laiobs)
        elseif (LIS_rc%dascaloption(k).eq."seasonal".and.fnd.ne.0) then

            call GenericLAI_rescale_with_seasonal_scaling(&
                 n,k,&
                 LIS_rc%da,&
                 GenericLAI_struc(n)%ntimes,        & 
                 MAX_LAI_VALUE,                      & 
                 MIN_LAI_VALUE,                      & 
                 GenericLAI_struc(n)%mult_scaling, &
                 GenericLAI_struc(n)%model_mu,  &
                 GenericLAI_struc(n)%model_sigma,  &
                 GenericLAI_struc(n)%obs_mu,  &
                 GenericLAI_struc(n)%obs_sigma,  &
                 laiobs)

        endif

        obsl = LIS_rc%udef 
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                         laiobs(c+(r-1)*LIS_rc%obs_lnc(k))
                endif
            enddo
        enddo

        if(fnd.eq.0) then 
            data_upd_flag_local = .false. 
        else
            data_upd_flag_local = .true. 
        endif

#if (defined SPMD)
        call MPI_ALLGATHER(data_upd_flag_local,1, &
             MPI_LOGICAL, data_upd_flag(:),&
             1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
        data_upd = .false.
        do p=1,LIS_npes
            data_upd = data_upd.or.data_upd_flag(p)
        enddo

        if(data_upd) then 

            do t=1,LIS_rc%obs_ngrid(k)
                gid(t) = t
                if(obsl(t).ne.-9999.0) then 
                    assimflag(t) = 1
                else
                    assimflag(t) = 0
                endif
            enddo

            call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                 .true. , rc=status)
            call LIS_verify(status)

            if(LIS_rc%obs_ngrid(k).gt.0) then 
                call ESMF_AttributeSet(laifield,"Grid Number",&
                     gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(laifield,"Assimilation Flag",&
                     assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

            endif

            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                     rc=status)
                call LIS_verify(status, 'Error: StateGet Observation01')

                ssdev = GenericLAI_struc(n)%ssdev_inp 

                if (LIS_rc%dascaloption(k).eq."CDF matching") then
                    if(GenericLAI_struc(n)%ntimes.eq.1) then 
                        timeidx = 1
                    else
                        timeidx = LIS_rc%mo
                    endif
                elseif(LIS_rc%dascaloption(k).eq."seasonal"&
                     .or.LIS_rc%dascaloption(k).eq."seasonal multiplicative") then
                    timeidx = LIS_rc%da
                endif

                call GenericLAI_updateSsdev(k,&
                     GenericLAI_struc(n)%obs_sigma(:, timeidx),&
                     GenericLAI_struc(n)%model_sigma(:, timeidx),&
                     ssdev)

                if(LIS_rc%obs_ngrid(k).gt.0) then 
                    call ESMF_AttributeSet(pertfield,"Standard Deviation",&
                         ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                    call LIS_verify(status)
                endif
            endif


        else
            call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                 .false., rc=status)
            call LIS_verify(status)     
        endif
    else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)     
    endif

end subroutine read_GenericLAI

!BOP
! 
! !ROUTINE: read_Generic_LAI_data
! \label{read_Generic_LAI_data}
!
! !INTERFACE:
subroutine read_Generic_LAI_data(n, k, fname, laiobs_ip)
    ! 
    ! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod
    use LIS_timeMgrMod
    use GenericLAI_Mod, only : GenericLAI_struc

    implicit none
    !
    ! !INPUT PARAMETERS: 
    ! 
    integer                       :: n 
    integer                       :: k
    character (len=*)             :: fname
    real                          :: laiobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real*8                        :: cornerlat(2), cornerlon(2)

    ! !OUTPUT PARAMETERS:
    !
    !
    ! !DESCRIPTION: 
    !  This subroutine reads the Generic LAI file and applies the data
    !  quality flags to filter the data. 
    !
    !  The arguments are: 
    !  \begin{description}
    !  \item[n]            index of the nest
    !  \item[k]            number of observation state
    !  \item[k]            number of observation state
    !  \item[fname]        name of the Generic LAI file
    !  \item[laiobs\_ip]   Generic LAI data processed to the LIS domain
    !  \end{description}
    !
    !
    !EOP

    integer                 :: lat_offset, lon_offset
    real                    :: lai(GenericLAI_struc(n)%nc,GenericLAI_struc(n)%nr)
    real                    :: lai_in(GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr)
    logical*1               :: lai_data_b(GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr)
    logical*1               :: laiobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                 :: c,r,t
    integer                 :: nid
    integer                 :: laiid, flagid
    integer                 :: ios
    character(len=20)       :: lai_name
    integer                 :: numvalidobs

    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                                :: numLons, numLats

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    !values

    cornerlat(1)=GenericLAI_struc(n)%gridDesci(4)
    cornerlon(1)=GenericLAI_struc(n)%gridDesci(5)
    cornerlat(2)=GenericLAI_struc(n)%gridDesci(7)
    cornerlon(2)=GenericLAI_struc(n)%gridDesci(8)

    lai_data_b = .false.

    lat_offset = 1  ! no offset
    lon_offset = 1


    ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
    call LIS_verify(ios,'Error opening file '//trim(fname))

    ios = nf90_inq_varid(nid, trim(GenericLAI_struc(n)%nc_varname), laiid)
    call LIS_verify(ios, 'Error nf90_inq_varid: '//GenericLAI_struc(n)%nc_varname)

    ios = nf90_get_var(nid, laiid, lai, &
         start=(/lon_offset,lat_offset/), &
         count=(/GenericLAI_struc(n)%nc,GenericLAI_struc(n)%nr/)) 

    call LIS_verify(ios, 'Error nf90_get_var: '//GenericLAI_struc(n)%nc_varname)

    ios = nf90_close(ncid=nid)
    call LIS_verify(ios,'Error closing file '//trim(fname))

    ! the data is already read into 'lai', but we have to replace
    ! NaNs/invalid values with LIS_rc%udef
    do r=1, GenericLAI_struc(n)%nr
        do c=1, GenericLAI_struc(n)%nc
            if (isnan(lai(c, r))) then
                lai(c, r) = LIS_rc%udef
            else if (.not. (0.0 < lai(c, r) .and. lai(c, r) < 20.0)) then
                lai(c, r) = LIS_rc%udef
            endif
        end do
    end do


    ! fill lai_in and lai_data_b, which are required further on
    do r=1, GenericLAI_struc(n)%nr
        do c=1, GenericLAI_struc(n)%nc
            lai_in(c+(r-1)*GenericLAI_struc(n)%nc) = lai(c,r)
            if(lai(c,r).ne.LIS_rc%udef) then
                lai_data_b(c+(r-1)*GenericLAI_struc(n)%nc) = .true.
            else
                lai_data_b(c+(r-1)*GenericLAI_struc(n)%nc) = .false.
            endif
        enddo
    enddo

    if(LIS_rc%obs_gridDesc(k,10).lt.GenericLAI_struc(n)%dlon) then 
        write(LIS_logunit,*) '[INFO] interpolating Generic LAI',trim(fname)
        !--------------------------------------------------------------------------
        ! Interpolate to the LIS running domain if model has finer resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
             GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             GenericLAI_struc(n)%rlat,GenericLAI_struc(n)%rlon,&
             GenericLAI_struc(n)%w11,GenericLAI_struc(n)%w12,&
             GenericLAI_struc(n)%w21,GenericLAI_struc(n)%w22,&
             GenericLAI_struc(n)%n11,GenericLAI_struc(n)%n12,&
             GenericLAI_struc(n)%n21,GenericLAI_struc(n)%n22,LIS_rc%udef,ios)
     else
        write(LIS_logunit,*) '[INFO] upscaling Generic LAI',trim(fname)
        !--------------------------------------------------------------------------
        ! Upscale to the LIS running domain if model has coarser resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call upscaleByAveraging(GenericLAI_struc(n)%nc*GenericLAI_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             LIS_rc%udef, GenericLAI_struc(n)%n11,&
             lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
    endif

#endif
end subroutine read_Generic_LAI_data


!BOP
! !ROUTINE: create_GenericLAI_filename
! \label{create_GenericLAI_filename}
! 
! !INTERFACE: 
subroutine create_GenericLAI_filename(prefix, ndir, year, month, day, filename)
    ! !USES:   

    implicit none
    ! !ARGUMENTS: 
    character (len=*), intent(in)     :: prefix
    character (len=*)                 :: ndir
    integer, value                    :: year, month, day
    character(len=*)                  :: filename
    ! 
    ! !DESCRIPTION: 
    !  This subroutine creates the Generic LAI filename
    !  based on the time and date 
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[varname] variable name of the netCDF variable
    !  \item[ndir] name of the Generic LAI data directory
    !  \item[version] version of the Generic LAI data
    !  \item[year]  current year
    !  \item[month]  current month
    !  \item[day]  current day
    !  \item[filename] Generated Generic LAI filename
    !  \end{description}
    !
    !EOP

    !  The naming scheme is
    !    <varname>_<year>_<month>_<day>.nc


    character(len=4) :: yearstr
    character(len=2) :: monthstr
    character(len=2) :: daystr

    write(unit=yearstr, fmt='(i4.4)') year
    write(unit=monthstr, fmt='(i2.2)') month
    write(unit=daystr, fmt='(i2.2)') day


    filename = trim(ndir)//'/'//&
         trim(prefix)//'_'//yearstr//'_'//monthstr//'_'//daystr//'.nc'

end subroutine create_GenericLAI_filename
