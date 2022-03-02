!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module VODObservationOperator_Mod
!BOP
!
! !MODULE: VODObservationOperator_Mod
!
! !DESCRIPTION:
!    This modules provides routines for an observation operator for VOD.
!    Although it's not a real RTM it programmatically works similarly and is
!    therefore implemented analogously to a RTM.
!
! !HISTORY:
! 02 Mar 2022: Samuel Scherrer; based on VODObservationOperator
!
! !USES:        

#if (defined RTMS)

  use ESMF
  use LIS_coreMod
  use LIS_RTMMod
  use LIS_logMod

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: VODObservationOperator_initialize
  public :: VODObservationOperator_f2t
  public :: VODObservationOperator_run
  public :: VODObservationOperator_output
  public :: VODObservationOperator_geometry 
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: vod_obsop_struc
!EOP
  type, public ::  vod_obsop_type_dec 
   
     character*256      :: parameter_fname

     real, allocatable :: laicoef(:)
     real, allocatable :: sm1coef(:)
     real, allocatable :: sm2coef(:)
     real, allocatable :: sm3coef(:)
     real, allocatable :: sm4coef(:)

     real, allocatable :: lone(:)
     real, allocatable :: late(:)
     !-------output------------!   
     real, allocatable :: VOD(:)
  end type vod_obsop_type_dec

  type(vod_obsop_type_dec), allocatable :: vod_obsop_struc(:) 

  SAVE

contains
!BOP
! 
! !ROUTINE: VODObservationOperator_initialize
! \label{VODObservationOperator_initialize}
! 
! !INTERFACE:
   subroutine VODObservationOperator_initialize()
! !USES:

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noahMP3.6-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noahMP3.6 from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readVODObservationOperatorcrd](\ref{readVODObservationOperatorcrd}) \newline
!
!EOP
   implicit none
    
   integer :: rc
   integer :: n,t
   integer :: ierr, ios
   integer , parameter :: OPEN_OK = 0
   character*128 :: message
   ! typically the parameters are available on the same grid as the observations
   real          :: laicoef(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

   !allocate memory for nest
   allocate(vod_obsop_struc(LIS_rc%nnest))

   do n=1,LIS_rc%nnest
       !allocate memory for all tile in current nest

      allocate(vod_obsop_struc(n)%laicoef(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))
      allocate(vod_obsop_struc(n)%sm1coef(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))
      allocate(vod_obsop_struc(n)%sm2coef(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))
      allocate(vod_obsop_struc(n)%sm3coef(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))
      allocate(vod_obsop_struc(n)%sm4coef(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))

      allocate(vod_obsop_struc(n)%lone(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))
      allocate(vod_obsop_struc(n)%late(LIS_rc%glbnpatch(n,LIS_rc%lsm_index)))

      allocate(vod_obsop_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

      call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
      call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 1")
      call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 2")
      call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 3")
      call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Layer 4")

   enddo 


   call ESMF_ConfigFindLabel(LIS_config, "VOD Observation Operator parameter file:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,vod_obsop_struc(n)%parameter_fname, rc=rc)
      call LIS_verify(rc, "VOD Observation Operator parameter file: not defined")
   enddo


   do n=1,LIS_rc%nnest

       !------------------------------------------------------------
       ! Read grid information and initialize resampling
       !------------------------------------------------------------
       lat_lower_left = 90 - 0.5 * vod_obsop_struc(n)%spatialres
       lat_upper_right = -90 + 0.5 * vod_obsop_struc(n)%spatialres
       lon_lower_left = -180 + 0.5 * vod_obsop_struc(n)%spatialres
       lon_upper_right = 180 - 0.5 * vod_obsop_struc(n)%spatialres
       ! dlat is positive since LIS will figure out that latitude is decreasing
       dlat = vod_obsop_struc(n)%spatialres
       dlon = vod_obsop_struc(n)%spatialres
       nr = nint(180.0 / vod_obsop_struc(n)%spatialres)
       nc = 2 * nr

       gridDesci(1) = 0  ! regular lat-lon grid
       gridDesci(2) = vod_obsop_struc(n)%nc
       gridDesci(3) = vod_obsop_struc(n)%nr 
       gridDesci(4) = vod_obsop_struc(n)%lat_lower_left
       gridDesci(5) = vod_obsop_struc(n)%lon_lower_left
       gridDesci(6) = 128
       gridDesci(7) = vod_obsop_struc(n)%lat_upper_right
       gridDesci(8) = vod_obsop_struc(n)%lon_upper_right
       gridDesci(9) = vod_obsop_struc(n)%dlat
       gridDesci(10) = vod_obsop_struc(n)%dlon
       gridDesci(20) = 64
       mi = nc*nr

       !-----------------------------------------------------------------------------
       !   Use interpolation if LIS is running finer than native resolution. 
       !-----------------------------------------------------------------------------
       if(LIS_rc%obs_gridDesc(k,10).lt.dlon) then

           allocate(rlat(nc*nr))
           allocate(rlon(nc*nr))
           allocate(n11(nc*nr))
           allocate(n12(nc*nr))
           allocate(n21(nc*nr))
           allocate(n22(nc*nr))
           allocate(w11(nc*nr))
           allocate(w12(nc*nr))
           allocate(w21(nc*nr))
           allocate(w22(nc*nr))

           write(LIS_logunit,*)&
                '[INFO] create interpolation input for Generic LAI'       

           call bilinear_interp_input_withgrid(vod_obsop_struc(n)%gridDesci(:), &
                LIS_rc%obs_gridDesc(k,:),&
                LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                vod_obsop_struc(n)%rlat, vod_obsop_struc(n)%rlon,&
                vod_obsop_struc(n)%n11, vod_obsop_struc(n)%n12, &
                vod_obsop_struc(n)%n21, vod_obsop_struc(n)%n22, &
                vod_obsop_struc(n)%w11, vod_obsop_struc(n)%w12, &
                vod_obsop_struc(n)%w21, vod_obsop_struc(n)%w22)
       else

           allocate(vod_obsop_struc(n)%n11(&
                vod_obsop_struc(n)%nc*vod_obsop_struc(n)%nr))

           write(LIS_logunit,*)&
                '[INFO] create upscaling input for Generic LAI'       

           call upscaleByAveraging_input(vod_obsop_struc(n)%gridDesci(:),&
                LIS_rc%obs_gridDesc(k,:),&
                vod_obsop_struc(n)%nc*vod_obsop_struc(n)%nr, &
                LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), vod_obsop_struc(n)%n11)

           write(LIS_logunit,*)&
                '[INFO] finished creating upscaling input for Generic LAI'       
       endif
   enddo

   

    do n=1, LIS_rc%nnest
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        ios = nf90_open(path=trim(vod_obsop_struc(n)%parameter_fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(vod_obsop_struc(n)%parameter_fname))

        ios = nf90_inq_varid(nid, trim(vod_obsop_struc(n)%laicoef_varname), laiid)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//vod_obsop_struc(n)%laicoef_varname)

        ios = nf90_get_var(nid, laiid, laicoef, &
             count=(/GenericLAI_struc(n)%nc,GenericLAI_struc(n)%nr/)) 


        ios = nf90_inq_varid(nid, trim(vod_obsop_struc(n)%sm1coef_varname), sm1id)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//vod_obsop_struc(n)%sm1coef_varname)
        ios = nf90_inq_varid(nid, trim(vod_obsop_struc(n)%sm2coef_varname), sm2id)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//vod_obsop_struc(n)%sm2coef_varname)
        ios = nf90_inq_varid(nid, trim(vod_obsop_struc(n)%sm3coef_varname), sm3id)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//vod_obsop_struc(n)%sm3coef_varname)
        ios = nf90_inq_varid(nid, trim(vod_obsop_struc(n)%sm4coef_varname), sm4id)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//vod_obsop_struc(n)%sm4coef_varname)
    enddo




!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%AA_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
           READ (19,*)vod_obsop_struc(n)%AA_VV(t),vod_obsop_struc(n)%lone(t),&
          vod_obsop_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%BB_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
           READ (19,*)vod_obsop_struc(n)%BB_VV(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%CC_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)  
           READ (19,*) vod_obsop_struc(n)%CC_VV(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%DD_VV_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
           READ (19,*) vod_obsop_struc(n)%DD_VV(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%AA_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
           READ (19,*)vod_obsop_struc(n)%AA_VH(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%BB_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
           READ (19,*)vod_obsop_struc(n)%BB_VH(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%CC_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)  
           READ (19,*) vod_obsop_struc(n)%CC_VH(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(vod_obsop_struc(n)%DD_VH_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
           READ (19,*) vod_obsop_struc(n)%DD_VH(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   end subroutine VODObservationOperator_initialize 
!!--------------------------------------------------------------------------------

   subroutine add_sfc_fields(n, sfcState,varname)

   implicit none 

   integer            :: n 
   type(ESMF_State)   :: sfcState
   character(len=*)   :: varname

   type(ESMF_Field)     :: varField
   type(ESMF_ArraySpec) :: arrspec
   integer              :: status
   real :: sum
   call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
   call LIS_verify(status)

   varField = ESMF_FieldCreate(arrayspec=arrSpec, & 
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
   call LIS_verify(status, 'Error in field_create of '//trim(varname))
    
   call ESMF_StateAdd(sfcState, (/varField/), rc=status)
   call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

   end subroutine add_sfc_fields


   subroutine VODObservationOperator_f2t(n)

   implicit none

   integer, intent(in)    :: n 

   end subroutine VODObservationOperator_f2t


   subroutine VODObservationOperator_geometry(n)
   implicit none
   integer, intent(in)    :: n

   end subroutine VODObservationOperator_geometry 
  !Do nothing for now
   subroutine VODObservationOperator_run(n)
   use LIS_histDataMod
! !USES: 
   implicit none

   integer, intent(in) :: n

   integer             :: t,p
   integer             :: status
   integer             :: col,row
   real                :: A_VV_cal,B_VV_cal,C_VV_cal,D_VV_cal,A_VH_cal,&
                          B_VH_cal, C_VH_cal, D_VH_cal,lon,lat,lon1,lat1
   real, pointer       :: sm(:), lai(:)
   real                :: sigmabare_VV,sigmabare_VH,s0VV_s_db, s0VH_s_dB, &
                        sigmacan_VV, sigmacan_VH,sigmasoil_VV,sigmasoil_VH,&
                        tt_VV, tt_VH

   real                :: theta, ctheta


   theta = 0.6458 !incidence angle in radians (37 deg)
   ctheta = cos(theta)


!   map surface properties to SFC    
   call getsfcvar(LIS_sfcState(n), "Soil Moisture Content",&
         sm)
   call getsfcvar(LIS_sfcState(n), "Leaf Area Index", &
         lai)


!---------------------------------------------
! Tile loop 
!--------------------------------------------
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbnpatch(n,LIS_rc%lsm_index)
          lon1= vod_obsop_struc(n)%lone(p)
          lat1= vod_obsop_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_cal=vod_obsop_struc(n)%AA_VV(p)
             B_VV_cal=vod_obsop_struc(n)%BB_VV(p)
             C_VV_cal=vod_obsop_struc(n)%CC_VV(p)
             D_VV_cal=vod_obsop_struc(n)%DD_VV(p)        
             A_VH_cal=vod_obsop_struc(n)%AA_VH(p)
             B_VH_cal=vod_obsop_struc(n)%BB_VH(p)
             C_VH_cal=vod_obsop_struc(n)%CC_VH(p)
             D_VH_cal=vod_obsop_struc(n)%DD_VH(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_s_db=C_VV_cal+D_VV_cal*sm(t)
          s0VH_s_db=C_VH_cal+D_VH_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV=10.**(s0VV_s_db/10.)
          sigmabare_VH=10.**(s0VH_s_db/10.)
         !attenuation
          tt_VV=exp(-2.*B_VV_cal*lai(t)/ctheta)
          tt_VH=exp(-2.*B_VH_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV=tt_VV*sigmabare_VV
          sigmasoil_VH=tt_VH*sigmabare_VH
         !vegetation backscatter in linear units
          sigmacan_VV=(1.-tt_VV)*ctheta*(A_VV_cal*lai(t))
          sigmacan_VH=(1.-tt_VH)*ctheta*(A_VH_cal*lai(t))
         !total backscatter
          vod_obsop_struc(n)%Sig0VV(t)=10.*log10(sigmacan_VV+sigmasoil_VV)
          vod_obsop_struc(n)%Sig0VH(t)=10.*log10(sigmacan_VH+sigmasoil_VH)
       else
          vod_obsop_struc(n)%Sig0VV(t)=LIS_rc%udef
          vod_obsop_struc(n)%Sig0VH(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          vod_obsop_struc(n)%Sig0VV(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          vod_obsop_struc(n)%Sig0VH(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   end subroutine VODObservationOperator_run


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine getsfcvar(sfcState, varname, var)
! !USES: 
    
   implicit none
    
   type(ESMF_State)      :: sfcState
   type(ESMF_Field)      :: varField
   character(len=*)      :: varname
   real, pointer         :: var(:)
   integer               :: status

   call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
   call LIS_verify(status, 'Error in StateGet: CMEM3_handlerMod '//trim(varname))
   call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
   call LIS_verify(status, 'Error in FieldGet: CMEM3_handlerMod '//trim(varname))

   end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: VODObservationOperator_output
!!!! \label{VODObservationOperator_output}
!!!!
!!!! !INTERFACE: 
   subroutine VODObservationOperator_output(n)
   integer, intent(in) :: n 
   end subroutine VODObservationOperator_output
#endif
end module VODObservationOperator_Mod


