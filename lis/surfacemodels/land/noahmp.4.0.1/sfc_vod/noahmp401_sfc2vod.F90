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
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE, PSISAT_TABLE, BEXP_TABLE

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
  type(ESMF_Field)    :: cwcfield, psifield, cvpdfield, tcfield, rzsmfield
  real, pointer       :: lai(:), sm1(:), sm2(:), sm3(:), sm4(:)
  real, pointer       :: cwc(:), psi(:), cvpd(:), tc(:), rzsm(:)
  integer             :: t,status, soiltyp, nroot, iz
  real, parameter     :: PSIWLT = -150.
  real                :: smcmax, smcwlt, psisat, sh2o, tmp_psi, dz, ztotal, bexp, eah, esat

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

  call ESMF_StateGet(sfcState,"Canopy Water Content",cwcfield,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(cwcfield,localDE=0,farrayPtr=cwc, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Canopy Temperature",tcfield,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(tcfield,localDE=0,farrayPtr=tc, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Canopy Vapor Pressure Deficit",cvpdfield,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(cvpdfield,localDE=0,farrayPtr=cvpd, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Root Zone Soil Water Potential",psifield,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(psifield,localDE=0,farrayPtr=psi, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Root Zone Soil Moisture",rzsmfield,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzsmfield,localDE=0,farrayPtr=rzsm, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      lai(t) = noahmp401_struc(n)%noahmp401(t)%lai
      sm1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
      sm2(t) = noahmp401_struc(n)%noahmp401(t)%smc(2)
      sm3(t) = noahmp401_struc(n)%noahmp401(t)%smc(3)
      sm4(t) = noahmp401_struc(n)%noahmp401(t)%smc(4)
      cwc(t) = noahmp401_struc(n)%noahmp401(t)%canliq
      tc(t) = noahmp401_struc(n)%noahmp401(t)%tah
      eah = noahmp401_struc(n)%noahmp401(t)%eah
      ! magnus formula
      esat = 610.94 * exp(17.625 * (tc(t) - 273.15) / ((tc(t) - 273.15) + 243.04))
      cvpd(t) = esat - eah

      ! calculate root-zone averaged PSI
      soiltyp = noahmp401_struc(n)%noahmp401(t)%soiltype
      smcmax = SMCMAX_TABLE(soiltyp)
      smcwlt = SMCWLT_TABLE(soiltyp)
      psisat = PSISAT_TABLE(soiltyp)
      bexp = BEXP_TABLE(soiltyp)
      nroot = noahmp401_struc(n)%noahmp401(t)%param%nroot
      ztotal = 0.0
      psi(t) = 0.0
      rzsm(t) = 0.0
      do iz = 1, nroot
          dz = noahmp401_struc(n)%sldpth(iz)
          ztotal = ztotal + dz
          sh2o = noahmp401_struc(n)%noahmp401(t)%sh2o(iz)
          tmp_psi = max(PSIWLT, -psisat * (max(0.01, sh2o)/smcmax)**(-bexp))
          psi(t) = psi(t) + dz * tmp_psi
          rzsm(t) = rzsm(t) + dz * sh2o
      enddo
      if (ztotal.gt.0) then
          psi(t) = psi(t) / ztotal
          rzsm(t) = rzsm(t) / ztotal
      endif
  enddo

end subroutine noahmp401_sfc2vod
