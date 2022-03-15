!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_CustomNetCDF
! \label{read_CustomNetCDF}
!
! !REVISION HISTORY:
!  02 Mar 2022    Samuel Scherrer; initial reader based on MODIS LAI reader
!
! !INTERFACE: 
subroutine read_CustomNetCDF(n, k, OBS_State, OBS_Pert_State, reader_struc)
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
    use CustomNcReader_Mod, only:&
         CustomNcReader_rescale_with_seasonal_scaling, &
         CustomNcReader_updateSsdev

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: k
    type(ESMF_State)    :: OBS_State
    type(ESMF_State)    :: OBS_Pert_State
    !
    ! !DESCRIPTION:
    !  
    !  reads the Custom observations from NETCDF files.

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
    integer                :: status
    character*100          :: obsdir
    character*300          :: fname
    integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
    real                   :: wt1, wt2,ts
    integer                :: count
    real                   :: cgmt
    real*8                 :: time
    logical                :: alarmCheck, file_exists,dataCheck
    integer                :: c,r,i,j,p,t
    real,          pointer :: obsl(:)
    type(ESMF_Field)       :: varfield
    integer                :: gid(LIS_rc%obs_ngrid(k))
    integer                :: assimflag(LIS_rc%obs_ngrid(k))
    logical                :: data_update
    logical                :: data_upd_flag(LIS_npes)
    logical                :: data_upd_flag_local
    logical                :: data_upd
    real                   :: observations(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                :: fnd
    real                   :: timenow
    real                   :: ssdev(LIS_rc%obs_ngrid(k))
    type(ESMF_Field)       :: pertfield
    integer                :: timeidx


    call ESMF_AttributeGet(OBS_State,"Data Directory",&
         obsdir, rc=status)
    call LIS_verify(status)
    call ESMF_AttributeGet(OBS_State,"Data Update Status",&
         data_update, rc=status)
    call LIS_verify(status)

    data_upd = .false. 

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Custom "&
         //trim(reader_struc(n)%varname//" read alarm")

    if(alarmCheck) then 
        call create_CustomNetCDF_filename(reader_struc(n)%nc_prefix,&
             obsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)

        inquire(file=fname,exist=file_exists)          
        if(file_exists) then 
            write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
            call read_CustomNetCDF_data(n,k, fname,observations)
            fnd = 1
        else
            fnd = 0 
            write(LIS_logunit,*) '[WARN] Missing observation file: ',trim(fname)
        endif
    else
        fnd = 0 
        observations = LIS_rc%udef
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

        call ESMF_StateGet(OBS_State,"Observation01",varfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        call ESMF_FieldGet(varfield,localDE=0,farrayPtr=obsl,rc=status)
        call LIS_verify(status, 'Error: FieldGet')

        obsl = LIS_rc%udef 
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                         observations(c+(r-1)*LIS_rc%obs_lnc(k))
                endif
            enddo
        enddo

        !-------------------------------------------------------------------------
        !  Transform data to the LSM climatology using a CDF-scaling approach
        !-------------------------------------------------------------------------     
        
        if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then

            call LIS_rescale_with_CDF_matching(     &
                 n,k,                               & 
                 reader_struc(n)%nbins,         & 
                 reader_struc(n)%ntimes,        & 
                 reader_struc(n)%max_value, &
                 reader_struc(n)%min_value, &
                 reader_struc(n)%model_xrange,  &
                 reader_struc(n)%obs_xrange,    &
                 reader_struc(n)%model_cdf,     &
                 reader_struc(n)%obs_cdf,       &
                 observations)
        elseif ((LIS_rc%dascaloption(k).eq."seasonal"&
                 .or.LIS_rc%dascaloption(k).eq."seasonal multiplicative")&
             .and.fnd.ne.0) then

            call CustomNcReader_rescale_with_seasonal_scaling(&
                 n,k,&
                 nint(LIS_get_curr_calday(LIS_rc, 0)), &
                 reader_struc(n)%ntimes,        & 
                 reader_struc(n)%max_value, &
                 reader_struc(n)%min_value, &
                 reader_struc(n)%mult_scaling, &
                 reader_struc(n)%model_mu,  &
                 reader_struc(n)%model_sigma,  &
                 reader_struc(n)%obs_mu,  &
                 reader_struc(n)%obs_sigma,  &
                 observations)

        endif

        obsl = LIS_rc%udef 
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                         observations(c+(r-1)*LIS_rc%obs_lnc(k))
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
                call ESMF_AttributeSet(varfield,"Grid Number",&
                     gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(varfield,"Assimilation Flag",&
                     assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

            endif

            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                     rc=status)
                call LIS_verify(status, 'Error: StateGet Observation01')

                ssdev = reader_struc(n)%ssdev_inp 

                if (LIS_rc%dascaloption(k).eq."CDF matching") then
                    if(reader_struc(n)%ntimes.eq.1) then 
                        timeidx = 1
                    else
                        timeidx = LIS_rc%mo
                    endif
                elseif(LIS_rc%dascaloption(k).eq."seasonal"&
                     .or.LIS_rc%dascaloption(k).eq."seasonal multiplicative") then
                    timeidx = nint(LIS_get_curr_calday(LIS_rc, 0))
                endif

                call CustomNcReader_updateSsdev(k,&
                     reader_struc(n)%obs_sigma(:, timeidx),&
                     reader_struc(n)%model_sigma(:, timeidx),&
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

end subroutine read_CustomNetCDF

!BOP
! 
! !ROUTINE: read_CustomNetCDF_data
! \label{read_CustomNetCDF_data}
!
! !INTERFACE:
subroutine read_CustomNetCDF_data(n, k, fname, observations_ip, reader_struc)
    ! 
    ! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod
    use LIS_timeMgrMod

    implicit none
    !
    ! !INPUT PARAMETERS: 
    ! 
    integer                       :: n 
    integer                       :: k
    character (len=*)             :: fname
    real                          :: observations_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real*8                        :: cornerlat(2), cornerlon(2)

    ! !OUTPUT PARAMETERS:
    !
    !
    ! !DESCRIPTION: 
    !  This subroutine reads the Custom netCDF file and applies the data
    !  quality flags to filter the data. 
    !
    !  The arguments are: 
    !  \begin{description}
    !  \item[n]            index of the nest
    !  \item[k]            number of observation state
    !  \item[k]            number of observation state
    !  \item[fname]        name of the netCDF file
    !  \item[observations\_ip]   Observation data processed to the LIS domain
    !  \end{description}
    !
    !
    !EOP

    integer                 :: lat_offset, lon_offset
    real                    :: observation(reader_struc(n)%nc,reader_struc(n)%nr)
    real                    :: obs_in(reader_struc(n)%nc*reader_struc(n)%nr)
    logical*1               :: obs_data_b(reader_struc(n)%nc*reader_struc(n)%nr)
    logical*1               :: observations_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                 :: c,r,t
    integer                 :: nid
    integer                 :: obsid, flagid
    integer                 :: ios
    integer                 :: numvalidobs

    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                                :: numLons, numLats

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    !values

    cornerlat(1)=reader_struc(n)%gridDesci(4)
    cornerlon(1)=reader_struc(n)%gridDesci(5)
    cornerlat(2)=reader_struc(n)%gridDesci(7)
    cornerlon(2)=reader_struc(n)%gridDesci(8)

    obs_data_b = .false.

    lat_offset = 1  ! no offset
    lon_offset = 1


    ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
    call LIS_verify(ios,'Error opening file '//trim(fname))

    ios = nf90_inq_varid(nid, trim(reader_struc(n)%nc_varname), obsid)
    call LIS_verify(ios, 'Error nf90_inq_varid: '//reader_struc(n)%nc_varname)

    ios = nf90_get_var(nid, obsid, observation, &
         start=(/lon_offset,lat_offset/), &
         count=(/reader_struc(n)%nc,reader_struc(n)%nr/)) 

    call LIS_verify(ios, 'Error nf90_get_var: '//reader_struc(n)%nc_varname)

    ios = nf90_close(ncid=nid)
    call LIS_verify(ios,'Error closing file '//trim(fname))

    ! the data is already read into 'observation', but we have to replace
    ! NaNs/invalid values with LIS_rc%udef
    do r=1, reader_struc(n)%nr
        do c=1, reader_struc(n)%nc
            if (isnan(observation(c, r))) then
                observation(c, r) = LIS_rc%udef
            else if (.not. (reader_struc(n)%min_value < observation(c, r) &
                 .and. observation(c, r) < reader_struc(n)%max_value)) then
                observation(c, r) = LIS_rc%udef
            endif
        end do
    end do


    ! fill obs_in and obs_data_b, which are required further on
    do r=1, reader_struc(n)%nr
        do c=1, reader_struc(n)%nc
            obs_in(c+(r-1)*reader_struc(n)%nc) = observation(c,r)
            if(observation(c,r).ne.LIS_rc%udef) then
                obs_data_b(c+(r-1)*reader_struc(n)%nc) = .true.
            else
                obs_data_b(c+(r-1)*reader_struc(n)%nc) = .false.
            endif
        enddo
    enddo

    if(LIS_rc%obs_gridDesc(k,10).lt.reader_struc(n)%dlon) then 
        write(LIS_logunit,*) '[INFO] interpolating Custom',&
             trim(reader_struc(n)%varname),&
             trim(fname)
        !--------------------------------------------------------------------------
        ! Interpolate to the LIS running domain if model has finer resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             obs_data_b, obs_in, observations_b_ip, observations_ip, &
             reader_struc(n)%nc*reader_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             reader_struc(n)%rlat,reader_struc(n)%rlon,&
             reader_struc(n)%w11,reader_struc(n)%w12,&
             reader_struc(n)%w21,reader_struc(n)%w22,&
             reader_struc(n)%n11,reader_struc(n)%n12,&
             reader_struc(n)%n21,reader_struc(n)%n22,LIS_rc%udef,ios)
     else
        write(LIS_logunit,*) '[INFO] upscaling Custom',&
             trim(reader_struc(n)%varname),&
             trim(fname)
        !--------------------------------------------------------------------------
        ! Upscale to the LIS running domain if model has coarser resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call upscaleByAveraging(reader_struc(n)%nc*reader_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             LIS_rc%udef, reader_struc(n)%n11,&
             obs_data_b,obs_in, observations_b_ip, observations_ip)
    endif

#endif
end subroutine read_CustomNetCDF_data


!BOP
! !ROUTINE: create_CustomNetCDF_filename
! \label{create_CustomNetCDF_filename}
! 
! !INTERFACE: 
subroutine create_CustomNetCDF_filename(prefix, ndir, year, month, day, filename)
    ! !USES:   

    implicit none
    ! !ARGUMENTS: 
    character (len=*), intent(in)     :: prefix
    character (len=*), intent(in)     :: ndir
    integer, intent(in)               :: year, month, day
    character(len=*), intent(out)     :: filename
    ! 
    ! !DESCRIPTION: 
    !  This subroutine creates the netCDF filename
    !  based on the time and date 
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[varname] variable name of the netCDF variable
    !  \item[ndir] name of the data directory
    !  \item[year]  current year
    !  \item[month]  current month
    !  \item[day]  current day
    !  \item[filename] Generated filename
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

end subroutine create_CustomNetCDF_filename
