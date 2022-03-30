!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: Noahmp401_sfc2vod
!  \label{Noahmp401_sfc2vod}
!
! !REVISION HISTORY:
!  4 Sep 2020: Sara Modanesi; Initial Specification
! 16 Mar 2022: Samuel Scherrer; adaption for VOD observation operator
! !INTERFACE:
subroutine noahmp401_sfc2vod(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod
  use LIS_logMod,    only : LIS_verify
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use Noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState

! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the noahmp401 specific surface variables
! to the VOD observation operator. 
!
!EOP
  type(ESMF_Field)    :: laiField, sm1field, sm2field, sm3field, sm4field
  real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
  integer             :: t,status

  call ESMF_StateGet(sfcState,"Leaf Area Index",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Moisture Layer 1",sm1field,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm1field,localDE=0,farrayPtr=sm1, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Moisture Layer 2",sm2field,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2field,localDE=0,farrayPtr=sm2, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Moisture Layer 3",sm3field,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3field,localDE=0,farrayPtr=sm3, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Moisture Layer 4",sm4field,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm4field,localDE=0,farrayPtr=sm4, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      lai(t) = noahmp401_struc(n)%noahmp401(t)%lai
      sm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
      sm2(t) = noahmp401_struc(n)%noahmp401(t)%smc(2)
      sm3(t) = noahmp401_struc(n)%noahmp401(t)%smc(3)
      sm4(t) = noahmp401_struc(n)%noahmp401(t)%smc(4)
  enddo

end subroutine noahmp401_sfc2vod
