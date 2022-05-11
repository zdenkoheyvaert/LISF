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
! !ROUTINE: ac70_getsoilmLAI
! \label{ac70_getsoilmLAI}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for ac70 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
! !INTERFACE:
subroutine ac70_getsoilmLAI(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture and LAI related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: laiField,lfmassField
  type(ESMF_Field)       :: ac70sm1Field
  type(ESMF_Field)       :: ac70sm2Field
  type(ESMF_Field)       :: ac70sm3Field
  type(ESMF_Field)       :: ac70sm4Field
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: lai(:)
  real, pointer          :: ac70soilm1(:)
  real, pointer          :: ac70soilm2(:)
  real, pointer          :: ac70soilm3(:)
  real, pointer          :: ac70soilm4(:)
  character*100          :: lsm_state_objs(5)


  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm4 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for LAI in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"AC70 Soil Moisture Layer 1",ac70sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"AC70 Soil Moisture Layer 2",ac70sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"AC70 Soil Moisture Layer 3",ac70sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"AC70 Soil Moisture Layer 4",ac70sm4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm4 in ac70_getsoilmLAI')

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in ac70_getsoilmLAI')
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm2 in ac70_getsoilmLAI')
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm3 in ac70_getsoilmLAI')
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm4 in ac70_getsoilmLAI')
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for LAI in ac70_getsoilmLAI')
  call ESMF_FieldGet(ac70sm1Field,localDE=0,farrayPtr=ac70soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for ac70sm1 in ac70_getsoilmLAI')
  call ESMF_FieldGet(ac70sm2Field,localDE=0,farrayPtr=ac70soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for ac70sm2 in ac70_getsoilmLAI')
  call ESMF_FieldGet(ac70sm3Field,localDE=0,farrayPtr=ac70soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for ac70sm3 in ac70_getsoilmLAI')
  call ESMF_FieldGet(ac70sm4Field,localDE=0,farrayPtr=ac70soilm4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for ac70sm4 in ac70_getsoilmLAI')


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = AC70_struc(n)%ac70(t)%smc(1)
     soilm2(t) = AC70_struc(n)%ac70(t)%smc(2)
     soilm3(t) = AC70_struc(n)%ac70(t)%smc(3)
     soilm4(t) = AC70_struc(n)%ac70(t)%smc(4)
     lai(t) = AC70_struc(n)%ac70(t)%lai
     ac70soilm1(t) = AC70_struc(n)%ac70(t)%ac70smc(1)
     ac70soilm2(t) = AC70_struc(n)%ac70(t)%ac70smc(2)
     ac70soilm3(t) = AC70_struc(n)%ac70(t)%ac70smc(3)
     ac70soilm4(t) = AC70_struc(n)%ac70(t)%ac70smc(4)
  enddo

end subroutine ac70_getsoilmLAI

