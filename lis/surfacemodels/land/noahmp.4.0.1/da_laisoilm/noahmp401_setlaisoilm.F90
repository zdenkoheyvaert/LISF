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
! !ROUTINE: NoahMP401_setlaisoilm
!  \label{NoahMP401_setlaisoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari: Modified for NoahMP401 
! 19 Jan 2022: Samuel Scherrer: Modified from da_soilm for da_laisoilm
! 
! Apply the update if it met the update conditions
! Update conditions: 
!                  1- Prior SM(sh2o) + increment > MIN_THRESHOLD 
!                  2- Prior SM(sh2o) + increment < sm_threshold
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states    


! !INTERFACE:
subroutine NoahMP401_setlaisoilm(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP401_lsmMod
  !use module_sf_noahlsm_36  !MN
  !use module_sf_noahmpdrv_401, only: parameters
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE


  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the LAI and soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP
  type(ESMF_Field)       :: laiField
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  real, pointer          :: lai(:)
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  integer                :: t, status
  real                   :: lfmass
  real                   :: delta1, delta2, delta3, delta4

  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: LAI failed in NoahMP401_setlaisoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in NoahMP401_setlaisoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in NoahMP401_setlaisoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in NoahMP401_setlaisoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in NoahMP401_setlaisoilm")

  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: LAI failed in NoahMP401_setlaisoilm")
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in NoahMP401_setlaisoilm")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in NoahMP401_setlaisoilm")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in NoahMP401_setlaisoilm")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in NoahMP401_setlaisoilm")


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      ! sm updates
      delta1 = soilm1(t) - noahmp_struc(n)%noahmp(t)%smc(1)
      delta2 = soilm2(t) - noahmp_struc(n)%noahmp(t)%smc(2)
      delta3 = soilm3(t) - noahmp_struc(n)%noahmp(t)%smc(3)
      delta4 = soilm4(t) - noahmp_struc(n)%noahmp(t)%smc(4)
      noahmp_struc(n)%noahmp(t)%smc(1) = soilm1(t)
      noahmp_struc(n)%noahmp(t)%smc(2) = soilm2(t)
      noahmp_struc(n)%noahmp(t)%smc(3) = soilm3(t)
      noahmp_struc(n)%noahmp(t)%smc(4) = soilm4(t)
      noahmp_struc(n)%noahmp(t)%sh2o(1) = noahmp_struc(n)%noahmp(t)%sh2o(1) + delta1
      noahmp_struc(n)%noahmp(t)%sh2o(2) = noahmp_struc(n)%noahmp(t)%sh2o(2) + delta2
      noahmp_struc(n)%noahmp(t)%sh2o(3) = noahmp_struc(n)%noahmp(t)%sh2o(3) + delta3
      noahmp_struc(n)%noahmp(t)%sh2o(4) = noahmp_struc(n)%noahmp(t)%sh2o(4) + delta4

      ! lai updates
      noahmp401_struc(n)%noahmp401(t)%lai = lai(t)
      lfmass = lai(t)*1000.0/(NOAHMP401_struc(n)%noahmp401(t)%param%sla)
      noahmp401_struc(n)%noahmp401(t)%lfmass = lfmass
  enddo

end subroutine NoahMP401_setlaisoilm

