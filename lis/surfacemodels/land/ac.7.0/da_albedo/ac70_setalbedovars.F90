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
! !ROUTINE: ac70_setalbedovars
! \label{ac70_setalbedovars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac70_setalbedovars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use ac70_lsmMod
  use AC_VEG_PARAMETERS_70
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the progognostic variables to noah's
!  model space. The state vector consists of veg
! 
!EOP

  type(ESMF_Field)       :: albdField,albiField
  integer                :: t
  integer                :: status
  real, pointer          :: albd(:), albi(:)
 
  call ESMF_StateGet(LSM_State,"ALBD",albdField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(albdField,localDE=0,farrayPtr=albd,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"ALBI",albiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(albiField,localDE=0,farrayPtr=albi,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     AC70_struc(n)%ac70(t)%albd(:) = albd(t) 
     AC70_struc(n)%ac70(t)%albi(:) = albi(t) 
!     AC70_struc(n)%ac70(t)%albd(1) = albd(t) 
!     AC70_struc(n)%ac70(t)%albd(2) = albi(t) 
!     AC70_struc(n)%ac70(t)%albi(1) = albd(t) 
!     AC70_struc(n)%ac70(t)%albi(2) = albi(t) 
     if(albd(t).ne.-9999.0.and.albi(t).ne.-9999.0) then 
        AC70_struc(n)%ac70(t)%alb_upd_flag = .true.
     else
        AC70_struc(n)%ac70(t)%alb_upd_flag = .false. 
     endif
  enddo

    
end subroutine ac70_setalbedovars


