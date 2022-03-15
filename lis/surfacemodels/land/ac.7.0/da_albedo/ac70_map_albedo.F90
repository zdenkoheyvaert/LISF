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
! !ROUTINE: ac70_map_albedo
! \label{ac70_map_albedo}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac70_map_albedo(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for albedo data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: albdIncrField
  type(ESMF_Field)         :: albiIncrField
  type(ESMF_Field)         :: obs_ALBD_field
  type(ESMF_Field)         :: obs_ALBI_field
  real, pointer            :: albdincr(:)
  real, pointer            :: albiincr(:)
  real, pointer            :: ALBDobs(:)
  real, pointer            :: ALBIobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)
  real, allocatable            :: ac70_albd(:)
  real, allocatable            :: ac70_albi(:)

  allocate(ac70_albd(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(ac70_albi(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(LSM_Incr_State,"ALBD",albdIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(albdIncrField,localDE=0,farrayPtr=albdincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"ALBI",albiIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(albiIncrField,localDE=0,farrayPtr=albiincr,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  
  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_ALBD_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_ALBD_field,localDE=0,farrayPtr=ALBDobs,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_ALBI_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_ALBI_field,localDE=0,farrayPtr=ALBIobs,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ac70_albd(t)  = ac70_struc(n)%ac70(t)%albd(1)
     ac70_albi(t)  = ac70_struc(n)%ac70(t)%albi(1)
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(ALBDobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index).ge.0) then 
        albdincr(t) = &
             ALBDobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)&
             - ac70_albd(t)
        albiincr(t) = &
             ALBIobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)&
             - ac70_albi(t)
     else
        albdincr(t) = 0 
        albiincr(t) = 0 
     endif
  enddo
!  stop
  deallocate(obs_state_objs)
  deallocate(ac70_albd)
  deallocate(ac70_albi)

end subroutine ac70_map_albedo
   
