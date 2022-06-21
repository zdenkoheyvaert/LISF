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
! !ROUTINE: ac70_updatesoilmLAI
!  \label{ac70_updatesoilmLAI}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari; Modified for ac70 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
!
! !INTERFACE:
subroutine ac70_updatesoilmLAI(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: CCiActualField
  type(ESMF_Field)       :: sm1IncrField
  type(ESMF_Field)       :: CCiActualIncrField

  real, pointer          :: soilm1(:)
  real, pointer          :: CCiActual(:)
  real, pointer          :: soilmIncr1(:)
  real, pointer          :: CCiActualincr(:)
  integer                :: t,i,m,gid
  integer                :: status

  real                   :: CCiActualtmp,CCiActualmax,CCiActualmin

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: CCiActualmean(LIS_rc%ngrid(n))
  integer                :: nCCiActualmean(LIS_rc%ngrid(n))

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")

  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")
  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")

  call ESMF_StateGet(LSM_State,"CCiActual",CCiActualField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: LSM_State, failed in ac70_updatesoilmLAI")
  call ESMF_FieldGet(CCiActualField,localDE=0,farrayPtr=CCiActual,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: CCiActualField failed in ac70_updatesoilmLAI")
 
  call ESMF_StateGet(LSM_Incr_State,"CCiActual",CCiActualIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: LSM_Incr_State CCiActual failed in ac70_updatesoilmLAI")
  call ESMF_FieldGet(CCiActualIncrField,localDE=0,farrayPtr=CCiActualincr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: CCiActualIncrField failed in ac70_updatesoilmLAI")

  call ESMF_AttributeGet(CCiActualField,"Max Value",CCiActualmax,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: CCiActualField Max Value failed in ac70_updatesoilmLAI")
  call ESMF_AttributeGet(CCiActualField,"Min Value",CCiActualmin,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: CCiActualField Min Value failed in ac70_updatesoilmLAI")


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) + soilmIncr1(t)
  enddo


  update_flag    = .true.
  perc_violation = 0.0
  CCiActualmean       = 0.0
  nCCiActualmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     CCiActualtmp =  CCiActual(t) + CCiActualincr(t)


     if(CCiActualtmp.lt.CCiActualmin.or.CCiActualtmp.gt.CCiActualmax) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid)/LIS_rc%nensem(n)
  enddo

! For ensembles that are unphysical, compute the
! ensemble average after excluding them. This
! is done only if the majority of the ensemble
! members are good (>60%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
           if((CCiActual(t)+CCiActualincr(t).gt.CCiActualmin).and.&
                (CCiActual(t)+CCiActualincr(t).lt.CCiActualmax)) then 
              CCiActualmean(gid) = CCiActualmean(gid) + &
                   CCiActual(t) + CCiActualincr(t)
              nCCiActualmean(gid) = nCCiActualmean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nCCiActualmean(gid).gt.0) then
        CCiActualmean(gid) = CCiActualmean(gid)/nCCiActualmean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     CCiActualtmp =  CCiActual(t) + CCiActualincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        CCiActual(t) = CCiActualtmp
     elseif(perc_violation(gid).lt.0.8) then
        if(CCiActualtmp.lt.CCiActualmin.or.CCiActualtmp.gt.CCiActualmax) then
           CCiActual(t) = CCiActualmean(gid)
        else
           CCiActual(t) = CCiActual(t) + CCiActualincr(t)
        endif
     endif
  enddo

end subroutine ac70_updatesoilmLAI

