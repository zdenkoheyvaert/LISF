!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac70_coldstart
! \label{Ac70_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and Ac70
!
! !INTERFACE:
subroutine Ac70_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use Ac70_lsmMod
   
!
! !DESCRIPTION:
!
!  This routine initializes the Ac70 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    
    ! added by Shugong Wang
    integer :: isnow
    real, allocatable, dimension(:) :: zsnso 
    real, allocatable, dimension(:) :: tsno
    real, allocatable, dimension(:) :: snice
    real, allocatable, dimension(:) :: snliq 
    real, allocatable, dimension(:) :: zsoil 
    ! end add

    !EMK...Temporary arrays.
    real :: tmp_swe(1,1),tmp_tgxy(1,1),tmp_snodep(1,1)
    integer :: tmp_isnowxy(1,1)

    do n=1, LIS_rc%nnest
       isnow = -AC70_struc(n)%nsnow
        ! added by shugong
        allocate(zsnso(-AC70_struc(n)%nsnow+1:AC70_struc(n)%nsoil))
        allocate(tsno(-AC70_struc(n)%nsnow+1:0))
        allocate(snice(-AC70_struc(n)%nsnow+1:0))
        allocate(snliq(-AC70_struc(n)%nsnow+1:0))
        allocate(zsoil(AC70_struc(n)%nsoil))
        zsoil(1) = -AC70_struc(n)%sldpth(1)
        do l=2, AC70_struc(n)%nsoil
          zsoil(l) = zsoil(l-1) - AC70_struc(n)%sldpth(l) 
        enddo
        ! end add 

        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: Ac70_coldstart -- cold-starting Ac70"
            do t=1, LIS_rc%npatch(n,mtype)
                AC70_struc(n)%ac70(t)%albold = AC70_struc(n)%init_albold
                AC70_struc(n)%ac70(t)%sneqvo = AC70_struc(n)%init_sneqvo
                ! only soil temperature is intialized, snow temperature is calculated by snow_init 
                do l=1, AC70_struc(n)%nsoil
                    AC70_struc(n)%ac70(t)%sstc(AC70_struc(n)%nsnow+l) = AC70_struc(n)%init_stc(l)
                enddo
                do l=1, AC70_struc(n)%nsoil
                    AC70_struc(n)%ac70(t)%sh2o(l) = AC70_struc(n)%init_sh2o(l)
                enddo
                do l=1, AC70_struc(n)%nsoil
                    AC70_struc(n)%ac70(t)%smc(l) = AC70_struc(n)%init_smc(l)
                enddo
                AC70_struc(n)%ac70(t)%tah = AC70_struc(n)%init_tah
                AC70_struc(n)%ac70(t)%eah = AC70_struc(n)%init_eah
                AC70_struc(n)%ac70(t)%fwet = AC70_struc(n)%init_fwet
                AC70_struc(n)%ac70(t)%canliq = AC70_struc(n)%init_canliq
                AC70_struc(n)%ac70(t)%canice = AC70_struc(n)%init_canice
                AC70_struc(n)%ac70(t)%tv = AC70_struc(n)%init_tv
                AC70_struc(n)%ac70(t)%tg = AC70_struc(n)%init_tg
                AC70_struc(n)%ac70(t)%qsnow = AC70_struc(n)%init_qsnow
                !do l=1, AC70_struc(n)%nsoil + AC70_struc(n)%nsnow
                !    AC70_struc(n)%ac70(t)%zss(l) = AC70_struc(n)%init_zss(l)
                !enddo
                AC70_struc(n)%ac70(t)%snowh = AC70_struc(n)%init_snowh
                AC70_struc(n)%ac70(t)%sneqv = AC70_struc(n)%init_sneqv
                !do l=1, AC70_struc(n)%nsnow
                !    AC70_struc(n)%ac70(t)%snowice(l) = AC70_struc(n)%init_snowice(l)
                !enddo
                !do l=1, AC70_struc(n)%nsnow
                !    AC70_struc(n)%ac70(t)%snowliq(l) = AC70_struc(n)%init_snowliq(l)
                !enddo
                AC70_struc(n)%ac70(t)%zwt = AC70_struc(n)%init_zwt
                AC70_struc(n)%ac70(t)%wa = AC70_struc(n)%init_wa
                AC70_struc(n)%ac70(t)%wt = AC70_struc(n)%init_wt
                AC70_struc(n)%ac70(t)%wslake = AC70_struc(n)%init_wslake
                AC70_struc(n)%ac70(t)%lfmass = AC70_struc(n)%init_lfmass
                AC70_struc(n)%ac70(t)%rtmass = AC70_struc(n)%init_rtmass
                AC70_struc(n)%ac70(t)%stmass = AC70_struc(n)%init_stmass
                AC70_struc(n)%ac70(t)%wood = AC70_struc(n)%init_wood
                AC70_struc(n)%ac70(t)%stblcp = AC70_struc(n)%init_stblcp
                AC70_struc(n)%ac70(t)%fastcp = AC70_struc(n)%init_fastcp
                AC70_struc(n)%ac70(t)%lai = AC70_struc(n)%init_lai
                AC70_struc(n)%ac70(t)%sai = AC70_struc(n)%init_sai
                AC70_struc(n)%ac70(t)%cm = AC70_struc(n)%init_cm
                AC70_struc(n)%ac70(t)%ch = AC70_struc(n)%init_ch
                AC70_struc(n)%ac70(t)%tauss = AC70_struc(n)%init_tauss
                AC70_struc(n)%ac70(t)%smcwtd = AC70_struc(n)%init_smcwtd
                AC70_struc(n)%ac70(t)%deeprech = AC70_struc(n)%init_deeprech
                AC70_struc(n)%ac70(t)%rech = AC70_struc(n)%init_rech
                AC70_struc(n)%ac70(t)%zlvl = AC70_struc(n)%init_zlvl 

                ! added by shugong 
                zsnso = 0.0 
                !EMK...snow_init_70 is expecting several arrays which
                !are being passed as scalars.  Although no memory corruption
                !occurs here because of the declared array dimensions (all 1),
                !this is still technically a syntax error.  So, we will
                !copy the required fields to temporary arrays with the
                !correct declarations and pass those instead.
                !call snow_init_70(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                !                  AC70_struc(n)%nsnow,          & !input 
                !                  AC70_struc(n)%nsoil,          & !input 
                !                  zsoil,                            & !input
                !                  AC70_struc(n)%init_sneqv,     & !input
                !                  AC70_struc(n)%init_tg,        & !input
                !                  AC70_struc(n)%init_snowh,     & !input
                !                  zsnso, tsno, snice, snliq, isnow) ! output 
                tmp_swe(1,1) = AC70_struc(n)%init_sneqv
                tmp_tgxy(1,1) = AC70_struc(n)%init_tg
                tmp_snodep(1,1) = AC70_struc(n)%init_snowh
                call snow_init_70(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                                  AC70_struc(n)%nsnow,          & !input 
                                  AC70_struc(n)%nsoil,          & !input 
                                  zsoil,                            & !input
                                  tmp_swe,         & !input
                                  tmp_tgxy,        & !input
                                  tmp_snodep,      & !input
                                  zsnso, tsno, snice, snliq, &
                                  tmp_isnowxy) ! output
                isnow = tmp_isnowxy(1,1)
                AC70_struc(n)%ac70(t)%snowice(1:AC70_struc(n)%nsnow) = snice(-AC70_struc(n)%nsnow+1:0)
                AC70_struc(n)%ac70(t)%snowliq(1:AC70_struc(n)%nsnow) = snliq(-AC70_struc(n)%nsnow+1:0)
                AC70_struc(n)%ac70(t)%zss(1:AC70_struc(n)%nsnow+AC70_struc(n)%nsoil) = zsnso(-AC70_struc(n)%nsnow+1:AC70_struc(n)%nsoil) 
                AC70_struc(n)%ac70(t)%sstc(AC70_struc(n)%nsnow+isnow+1:AC70_struc(n)%nsnow) = tsno(isnow+1:0) 
                AC70_struc(n)%ac70(t)%isnow = isnow                

                
                ! end add
            enddo
        endif
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "MSG: Ac70_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
        deallocate(zsnso)
        deallocate(tsno)
        deallocate(snice)
        deallocate(snliq)
        deallocate(zsoil) 
    enddo
end subroutine Ac70_coldstart
