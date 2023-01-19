!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp36_setsnowsimplevars
! \label{noahmp36_setsnowsimplevars}
!
! !REVISION HISTORY:
! 15 Aug 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_setsnowsimplevars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noahmp36_lsmMod
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. The state vector consists of total SWE and snow depth. 
!  This routine also updates other model prognostics (snice, snliq,
!  snow thickness, snow temperature) based on the update. 
! 
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField

  integer                :: t
  integer                :: status
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh

  real                   :: swe_new, snow_dens
  
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

!  write(LIS_logunit,*)'[INFO] upd in:SWE',NOAHMP36_struc(1)%noahmp36(10)%sneqv,'SD',NOAHMP36_struc(1)%noahmp36(10)%snowh
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     dsneqv = swe(t)*1000.0 - noahmp36_struc(n)%noahmp36(t)%sneqv !in mm
     dsnowh = snod(t) - noahmp36_struc(n)%noahmp36(t)%snowh  !in m

!NOTE: if snow becomes unphysical after applying the deltas, then 
! the update is rejected. This may cause inconsistent ensemble 
! updates. TBD
!
#if 0 
!only update dshowh, update SWE based on density
     dsnowh = snod(t) - noahmp36_struc(n)%noahmp36(t)%snowh  !in m
     snow_dens = noahmp36_struc(n)%noahmp36(t)%sneqv/&
          noahmp36_struc(n)%noahmp36(t)%snowh
     swe_new = snow_dens*snod(t)
     
     dsneqv = swe_new - noahmp36_struc(n)%noahmp36(t)%sneqv !in mm
#endif
     
     call noahmp36_snowsimple_update(n, t, dsneqv, dsnowh)
 
     if(noahmp36_struc(n)%noahmp36(t)%sneqv.lt.0.or.&
          noahmp36_struc(n)%noahmp36(t)%snowh.lt.0) then 
        print*, dsneqv, dsnowh
        print*, swe(t), snod(t)
        stop
     endif
  enddo
!  write(LIS_logunit,*)'[INFO] upd out:SWE',NOAHMP36_struc(1)%noahmp36(10)%sneqv,'SD',NOAHMP36_struc(1)%noahmp36(10)%snowh
end subroutine noahmp36_setsnowsimplevars


