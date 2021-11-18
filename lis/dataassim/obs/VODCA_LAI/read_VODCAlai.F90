!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_VODCAlai
! \label{read_VODCAlai}
!
! !REVISION HISTORY:
!  18 Nov 2021    Samuel Scherrer; initial reader based on CGLS LAI reader
!
! !INTERFACE: 
subroutine read_VODCAlai(n, k, OBS_State, OBS_Pert_State)
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
    use VODCAlai_Mod, only : VODCAlai_struc

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: k
    type(ESMF_State)    :: OBS_State
    type(ESMF_State)    :: OBS_Pert_State
    !
    ! !DESCRIPTION:
    !  
    !  reads the VODCA LAI observations from NETCDF files.

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
    integer                :: grid_index
    character*100          :: laiobsdir
    character*300          :: fname
    integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
    real                   :: wt1, wt2,ts
    integer                :: count
    real                   :: cgmt
    real*8                 :: time
    logical                :: alarmCheck, file_exists,dataCheck
    integer                :: t,c,r,i,j,p,jj
    real,          pointer :: obsl(:)
    type(ESMF_Field)       :: laifield, pertField
    integer                :: gid(LIS_rc%obs_ngrid(k))
    integer                :: assimflag(LIS_rc%obs_ngrid(k))
    logical                :: data_update
    logical                :: data_upd_flag(LIS_npes)
    logical                :: data_upd_flag_local
    logical                :: data_upd
    real                   :: laiobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                :: fnd
    real                   :: timenow


    call ESMF_AttributeGet(OBS_State,"Data Directory",&
         laiobsdir, rc=status)
    call LIS_verify(status)
    call ESMF_AttributeGet(OBS_State,"Data Update Status",&
         data_update, rc=status)
    call LIS_verify(status)

    data_upd = .false. 

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "VODCA LAI read alarm")

    if(alarmCheck.or.VODCAlai_struc(n)%startMode) then 
        VODCAlai_struc(n)%startMode = .false.

        call create_VODCAlai_filename(VODCAlai_struc(n)%band,&
             laiobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)

        inquire(file=fname,exist=file_exists)          
        if(file_exists) then 
            write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
            call read_VODCA_LAI_data(n,k, fname,laiobs)
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
end subroutine read_VODCAlai

!BOP
! 
! !ROUTINE: read_VODCA_LAI_data
! \label{read_VODCA_LAI_data}
!
! !INTERFACE:
subroutine read_VODCA_LAI_data(n, k, fname, laiobs_ip)
    ! 
    ! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod
    use LIS_timeMgrMod
    use VODCAlai_Mod, only : VODCAlai_struc

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
    !  This subroutine reads the VODCA LAI file and applies the data
    !  quality flags to filter the data. 
    !
    !  The arguments are: 
    !  \begin{description}
    !  \item[n]            index of the nest
    !  \item[k]            number of observation state
    !  \item[k]            number of observation state
    !  \item[fname]        name of the VODCA LAI file
    !  \item[laiobs\_ip]   VODCA LAI data processed to the LIS domain
    !  \end{description}
    !
    !
    !EOP

    integer                 :: lat_offset, lon_offset
    real                    :: lai(VODCAlai_struc(n)%nc,VODCAlai_struc(n)%nr)
    real                    :: lai_in(VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr)
    logical*1               :: lai_data_b(VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr)
    logical*1               :: laiobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                 :: c,r,t
    integer                 :: nid
    integer                 :: laiid, flagid
    integer                 :: ios

    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                                :: numLons, numLats

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    !values

    cornerlat(1)=VODCAlai_struc(n)%gridDesci(4)
    cornerlon(1)=VODCAlai_struc(n)%gridDesci(5)
    cornerlat(2)=VODCAlai_struc(n)%gridDesci(7)
    cornerlon(2)=VODCAlai_struc(n)%gridDesci(8)

    lai_data_b = .false.

    lat_offset = 1  ! no offset
    lon_offset = 1


    ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
    call LIS_verify(ios,'Error opening file '//trim(fname))

    ios = nf90_inq_varid(nid, 'VODCA_LAI',laiid)
    call LIS_verify(ios, 'Error nf90_inq_varid: VODCA_LAI')

    ios = nf90_get_var(nid, laiid, lai, &
         start=(/lon_offset,lat_offset/), &
         count=(/VODCAlai_struc(n)%nc,VODCAlai_struc(n)%nr/)) 

    call LIS_verify(ios, 'Error nf90_get_var: VODCA_LAI')

    ios = nf90_close(ncid=nid)
    call LIS_verify(ios,'Error closing file '//trim(fname))

    ! the data is already read into 'lai', but we have to replace
    ! NaNs/invalid values with LIS_rc%udef
    do r=1, VODCAlai_struc(n)%nr
        do c=1, VODCAlai_struc(n)%nc
            if (isnan(lai(c, r))) then
                lai(c, r) = LIS_rc%udef
            else if (.not. (0.0 < lai(c, r) .and. lai(c, r) < 20.0)) then
                lai(c, r) = LIS_rc%udef
            endif
        end do
    end do


    ! fill lai_in and lai_data_b, which are required further on
    do r=1, VODCAlai_struc(n)%nr
        do c=1, VODCAlai_struc(n)%nc
            lai_in(c+(r-1)*VODCAlai_struc(n)%nc) = lai_flagged(c,r)
            if(lai_flagged(c,r).ne.LIS_rc%udef) then
                lai_data_b(c+(r-1)*VODCAlai_struc(n)%nc) = .true.
            else
                lai_data_b(c+(r-1)*VODCAlai_struc(n)%nc) = .false.
            endif
        enddo
    enddo

    if(LIS_rc%obs_gridDesc(k,10).lt.VODCAlai_struc(n)%dlon) then 
        write(LIS_logunit,*) '[INFO] interpolating VODCA LAI',trim(fname)
        !--------------------------------------------------------------------------
        ! Interpolate to the LIS running domain if model has finer resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
             VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             VODCAlai_struc(n)%rlat,VODCAlai_struc(n)%rlon,&
             VODCAlai_struc(n)%w11,VODCAlai_struc(n)%w12,&
             VODCAlai_struc(n)%w21,VODCAlai_struc(n)%w22,&
             VODCAlai_struc(n)%n11,VODCAlai_struc(n)%n12,&
             VODCAlai_struc(n)%n21,VODCAlai_struc(n)%n22,LIS_rc%udef,ios)
     else if(LIS_rc%obs_gridDesc(k,10).gt.VODCAlai_struc(n)%dlon) then
        write(LIS_logunit,*) '[INFO] upscaling VODCA LAI',trim(fname)
        !--------------------------------------------------------------------------
        ! Upscale to the LIS running domain if model has coarser resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call upscaleByAveraging(VODCAlai_struc(n)%nc*VODCAlai_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             LIS_rc%udef, VODCAlai_struc(n)%n11,&
             lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
    else
        write(LIS_logunit,*) '[INFO] VODCA LAI already at correct resolution',trim(fname)
    endif

#endif
end subroutine read_VODCA_LAI_data


!BOP
! !ROUTINE: create_VODCAlai_filename
! \label{create_VODCAlai_filename}
! 
! !INTERFACE: 
subroutine create_VODCAlai_filename(band, ndir, year, month, day, filename)
    ! !USES:   

    implicit none
    ! !ARGUMENTS: 
    character, value     :: band
    character (len=*)    :: ndir
    integer, value       :: year, month, day
    character(len=*)     :: filename
    ! 
    ! !DESCRIPTION: 
    !  This subroutine creates the VODCA LAI filename
    !  based on the time and date 
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[band] frequency band
    !  \item[ndir] name of the VODCA LAI data directory
    !  \item[version] version of the VODCA LAI data
    !  \item[year]  current year
    !  \item[month]  current month
    !  \item[day]  current day
    !  \item[filename] Generated VODCA LAI filename
    !  \end{description}
    !
    !EOP

    !  The naming scheme is
    !    VODCA_<band>_LAI_<year>_<month>_<day>.nc


    character(len=4) :: yearstr
    character(len=2) :: monthstr
    character(len=2) :: daystr

    write(unit=yearstr, fmt='(i4.4)') year
    write(unit=monthstr, fmt='(i2.2)') month
    write(unit=daystr, fmt='(i2.2)') day


    filename = trim(ndir)//'/'//&
         'VODCA_'//band//'_LAI_'//yearstr//'_'//monthstr//'_'//daystr//'.nc'

end subroutine create_VODCAlai_filename
