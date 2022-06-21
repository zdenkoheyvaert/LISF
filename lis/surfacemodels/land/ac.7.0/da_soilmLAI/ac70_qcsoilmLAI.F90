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
! !ROUTINE: ac70_qcsoilmLAIveg
! \label{ac70_qcsoilmLAIveg}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for ac70 
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
!
! !INTERFACE:
subroutine ac70_qcsoilmLAI(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
!  type(ESMF_Field)       :: sm2Field
!  type(ESMF_Field)       :: sm3Field
!  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: CCiActualField
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
!  real, pointer          :: soilm2(:)
!  real, pointer          :: soilm3(:)
!  real, pointer          :: soilm4(:)
  real, pointer          :: CCiActual(:)
  real                   :: smmax1!,smmax2,smmax3,smmax4
  real                   :: smmin1!,smmin2,smmin3,smmin4
  real                   :: CCiActualmax
  real                   :: CCiActualmin
  integer                :: gid
  real                   :: CCiActualtmp

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))
  real                   :: CCiActualmean(LIS_rc%ngrid(n))
  integer                :: nCCiActualmean(LIS_rc%ngrid(n))
  integer                :: N_ens
  real                   :: state_tmp(LIS_rc%nensem(n)),state_mean

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in ac70_qcsoilmLAI")
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 1 failed in ac70_qcsoilmLAI")

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in ac70_qcsoilmLAI")
  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in ac70_qcsoilmLAI")

  call ESMF_StateGet(LSM_State,"CCiActual",CCiActualField,rc=status)
  call LIS_verify(status,&
           "ESMF_StateGet for CCiActual failed in ac70_qcsoilmLAI")
  call ESMF_FieldGet(CCiActualField,localDE=0,farrayPtr=CCiActual,rc=status)
  call LIS_verify(status,&
           "ESMF_FieldGet for CCiActual failed in ac70_qcsoilmLAI")

  call ESMF_AttributeGet(CCiActualField,"Max Value",CCiActualmax,rc=status)
  call LIS_verify(status,&
           "ESMF_AttributeGet for CCiActual Max Value failed in ac70_qcsoilmLAI")
  call ESMF_AttributeGet(CCiActualField,"Min Value",CCiActualmin,rc=status)
  call LIS_verify(status,&
           "ESMF_AttributeGet for CCiActual Min Value failed in ac70_qcsoilmLAI")



  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(soilm1(t).gt.smmax1) soilm1(t) = smmax1
     if(soilm1(t).lt.smmin1) soilm1(t) = smmin1
  enddo

  update_flag    = .true.
  perc_violation = 0.0
  CCiActualmean       = 0.0
  nCCiActualmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     CCiActualtmp =  CCiActual(t)

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
           if((CCiActual(t).gt.CCiActualmin).and.&
                (CCiActual(t).lt.CCiActualmax)) then 
              CCiActualmean(gid) = CCiActualmean(gid) + &
                   CCiActual(t) 
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

     CCiActualtmp =  CCiActual(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        CCiActual(t) = CCiActualtmp
     elseif(perc_violation(gid).lt.0.8) then
        if(CCiActualtmp.lt.CCiActualmin.or.CCiActualtmp.gt.CCiActualmax) then
           CCiActual(t) = CCiActualmean(gid)
        else
           CCiActual(t) = CCiActual(t) 
        endif
     endif
  enddo

#if 0 
  N_ens = LIS_rc%nensem(n)
  do t=1,N_ens
     state_tmp(t) = CCiActual(t)
  enddo
  state_mean =sum(state_tmp)/N_ens

  write(113,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,21F8.3)') &
       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
       LIS_rc%mn, LIS_rc%ss, &
       state_mean, &
       state_tmp
#endif
end subroutine ac70_qcsoilmLAI

