!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_qc_LAIobs
! \label{noahmp401_qc_LAIobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 1 Aug 2016: Mahdi Navari; Modified for Noahmp401 
! 14 Apr 2022: Samuel Scherrer; added check for large updates
!
! !INTERFACE:
subroutine noahmp401_qc_LAIobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE


  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil moisture observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen or (3) ground is fully or partially covered 
!  with snow MN:(4) ground is covered with vegatation (more than 50%). 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)         :: obs_field
  real, pointer            :: obs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status
  real                     :: ninnov, mu, mu_old, std, count, val
  real                     :: forecast(LIS_rc%ngrid(n))
  real                     :: spread(LIS_rc%ngrid(n))
  real                     :: forecast_obsspace(LIS_rc%obs_ngrid(k))
  real                     :: spread_obsspace(LIS_rc%obs_ngrid(k))



  call ESMF_StateGet(OBS_State,"Observation01",obs_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in NoahMP401_qc_LAIobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=obs,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in NoahMP401_qc_LAIobs")

  !-------------------------------------------------------------------
  ! FORECAST AND SPREAD CALCULATION
  !-------------------------------------------------------------------

  do i=1, LIS_rc%npatch(n,LIS_rc%lsm_index), LIS_rc%nensem
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 

     ! calculate mean and standard deviation using Welford's algorithm
     mu = 0
     count = 0
     std = 0
     do m=1,LIS_rc%nensem(n)
         t = i+m-1
         val = noahmp401_struc(n)%noahmp401(t)%lai
         if (val .ne. LIS_rc%udef) then
             mu_old = mu
             mu = mu + (val - mu) / count
             std = std + (val - mu_old) * (val - mu)
             count = count + 1
         endif
     enddo

     if (count .ne. 0) then
         forecast(gid) = mu
         spread(gid) = std / (count - 1)
     else
         forecast(gid) = LIS_rc%udef
         spread(gid) = LIS_rc%udef
     endif

  enddo

  gridSpaceToObsSpace(n, k, forecast, forecast_obsspace)
  gridSpaceToObsSpace(n, k, spread, spread_obsspace)

  !-------------------------------------------------------------------
  ! REJECT UPDATES WITH HIGH NINNOV VALUE
  !-------------------------------------------------------------------

  do t = 1,LIS_rc%obs_ngrid(k)

      if (obs(t).ne.LIS_rc%udef.and.forecast_obsspace(t).ne.LIS_rc%udef &
           .and.spread_obsspace(t).ne.LIS_rc%udef) then
          ninnov = (obs(t) - forecast_obsspace(t)) / spread_obsspace(t)

          if (ninnov > 10) then
              obs(t) = LIS_rc%udef
          endif
      endif
  enddo


contains

    ! adapted from LIS_convertPatchSpaceToObsSpace
    subroutine gridSpaceToObsSpace(
       n,&
       k,&
       gvar, &
       ovar)
!
! !ARGUMENTS: 
    integer,          intent(in) :: n 
    integer,          intent(in) :: k
    real                         :: gvar(LIS_rc%ngrid(n))
    real                         :: ovar(LIS_rc%obs_ngrid(k))

    integer                      :: g
    logical*1                    :: li(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real                         :: obs_gvar(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                      :: iret

    li = .false. 
    do g=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
       if(gvar(g).ne.LIS_rc%udef) then 
          li(g) = .true. 
       endif
    enddo

    if(LIS_isatAfinerResolution(n,LIS_obs_domain(n,k)%datares)) then     
       call upscaleByWeightedAveraging(&
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
            LIS_rc%udef, &
            LIS_obs_domain(n,k)%nbr_index, &
            LIS_obs_domain(n,k)%weight, &
            li, &
            gvar, &
            lo, &
            obs_gvar)
    else
       call neighbor_interp(LIS_rc%obs_gridDesc(k,:), &
            li, lis_gvar, lo, obs_gvar,&
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
            LIS_obs_domain(n,k)%rlat, &
            LIS_obs_domain(n,k)%rlon, &
            LIS_obs_domain(n,k)%nbr_index, &
            LIS_rc%udef,iret)
    endif
    ovar = LIS_rc%udef
    do r=1,LIS_rc%obs_lnr(k)
       do c=1,LIS_rc%obs_lnc(k)          
          if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
             ovar(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                  obs_gvar(c+(r-1)*LIS_rc%obs_lnc(k))
          endif
       enddo
    enddo

    end subroutine gridSpaceToObsSpace


end subroutine noahmp401_qc_LAIobs

