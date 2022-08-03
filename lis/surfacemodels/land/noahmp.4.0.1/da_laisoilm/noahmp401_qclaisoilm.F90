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
! !ROUTINE: NoahMP401_qclaisoilm
! \label{NoahMP401_qclaisoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 19 Jan 2022: Samuel Scherrer: Modified from da_soilm for da_laisoilm
!
! !INTERFACE:
subroutine NoahMP401_qclaisoilm(n, LSM_State)

    ! !USES:
    use ESMF
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
    use LIS_logMod,  only  : LIS_verify
    use NoahMP401_lsmMod
    use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in)    :: n
    type(ESMF_State)       :: LSM_State
    !
    ! !DESCRIPTION:
    !
    !  Returns the LAI and soilmoisture related state prognostic variables for
    !  data assimilation
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[n] index of the nest \newline
    !  \item[LSM\_State] ESMF State container for LSM state variables \newline
    !  \end{description}
    !EOP
    integer                :: t, gid, status
    integer                :: status
    type(ESMF_Field)       :: laiField
    type(ESMF_Field)       :: sm1Field
    type(ESMF_Field)       :: sm2Field
    type(ESMF_Field)       :: sm3Field
    type(ESMF_Field)       :: sm4Field
    real, pointer          :: lai(:)
    real, pointer          :: soilm1(:), soilm2(:), soilm3(:), soilm4(:)
    real                   :: delta(4), bias(4)
    real                   :: laimax
    real                   :: laimin
    integer                :: SOILTYP           ! soil type index [-]
    real                   :: SMCMAX , SMCWLT
    real, parameter        :: MIN_THRESHOLD = 0.02 
    real                   :: MAX_threshold
    real                   :: sm_threshold
    real                   :: tmp

    logical                :: update_flag(LIS_rc%ngrid(n))
    logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
    logical                :: update_flag_ens(LIS_rc%ngrid(n))

    logical                :: bound_violation
    integer                :: niter
    logical                :: ens_flag(LIS_rc%nensem-1)


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

    call ESMF_AttributeGet(laiField,"Max Value",laimax,rc=status)
    call LIS_verify(status,&
         "ESMF_AttributeGet: LAI Max Value failed in NoahMP401_qclaisoilm")
    call ESMF_AttributeGet(laiField,"Min Value",laimin,rc=status)
    call LIS_verify(status,&
         "ESMF_AttributeGet: LAI Min Value failed in NoahMP401_qclaisoilm")

    update_flag = .true. 
    update_flag_tile= .true. 

    ! in a first loop we figure out which of the new updated values are within
    ! the allowed bounds
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

        gid = LIS_domain(n)%gindex(&
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
             LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 

        SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype        
        MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)   ! MAXSMC (SOILTYP)
        sm_threshold = SMCMAX_TABLE(SOILTYP) - 0.02  ! MAXSMC (SOILTYP) - 0.02

        ! the updatelaisoilm routine already added the increments to soilmX, the
        ! full volumetric water content. Since we also check the liquid water
        ! content we need to have the increments.
        delta(1) = soilm1(t)-noahmp401_struc(n)%noahmp401(t)%smc(1)
        delta(2) = soilm2(t)-noahmp401_struc(n)%noahmp401(t)%smc(2)
        delta(3) = soilm3(t)-noahmp401_struc(n)%noahmp401(t)%smc(3)
        delta(4) = soilm4(t)-noahmp401_struc(n)%noahmp401(t)%smc(4)

        do j=1,4
            ! updated liquid water content
            tmp = noahmp401_struc(n)%noahmp401(t)%sh2o(j) + delta(j)
            ! if the updated liquid water content is between 0.02 and smcmax-0.02,
            ! we accept the update
            if (tmp.gt.MIN_THRESHOLD .and. tmp.lt.sm_threshold) then
                update_flag(gid) = update_flag(gid).and.(.true.)
            else
                update_flag(gid) = update_flag(gid).and.(.false.)
                update_flag_tile(t) = update_flag_tile(t).and.(.false.)
            endif
        enddo
        if(lai(t).ge.laimin .and. lai(t).le.laimax) then
            update_flag(gid) = update_flag(gid).and.(.true.)
            update_flag_tile(t) = update_flag_tile(t).and.(.true.)
        else
            update_flag(gid) = update_flag(gid).and.(.false.)
            update_flag_tile(t) = update_flag_tile(t).and.(.false.)
        endif
    enddo

    ! in case there are bounds violations for a grid cell, we can now count the
    ! number of violations and if there are more than 50% valid updates we can
    ! still do an update and replace the invalid tiles with the mean of the valid
    ! tiles
    update_flag_ens = .True.
    do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
        gid = LIS_domain(n)%gindex(&
             LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
             LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
        ! update_flag_ens = (num(valid tiles) > 50%) 
        update_flag_ens(gid) = (COUNT(update_flag_tile(i:i+LIS_rc%nensem(n)-1)).ge.LIS_rc%nensem(n)*0.5)
    enddo


    ! now we can adapt the values if necessary
    do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
        gid = LIS_domain(n)%gindex(&
             LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
             LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 

        if (update_flag(gid)) then
            ! lai and sm updates are all okay, so we don't have to apply any
            ! changes and can jump to the next grid cell
            cycle
        elseif (update_flag_ens(gid)) then
            ! there are some tiles that have invalid values, but we replace them
            ! with the mean of the valid cells
            meansoilm1 = 0
            meansoilm2 = 0
            meansoilm3 = 0
            meansoilm4 = 0
            meanlai = 0.
            nmean = 0
            do m=1,LIS_rc%nensem(n)
                t = i+m-1
                if (update_flag_tile(t)) then
                    nmean = nmean + 1
                    meanlai = meanlai  + (lai(t) - meanlai) / n
                    meansoilm1 = meansoilm1 + (soilm1(t) - meansoilm1) / n
                    meansoilm2 = meansoilm2 + (soilm2(t) - meansoilm2) / n
                    meansoilm3 = meansoilm3 + (soilm3(t) - meansoilm3) / n
                    meansoilm4 = meansoilm4 + (soilm4(t) - meansoilm4) / n
                endif
            enddo

            do m=1,LIS_rc%nensem(n)
                t = i+m-1
                if (.not.update_flag_tile(t)) then
                    lai(t) = meanlai
                    soilm1(t) = meansoilm1
                    soilm2(t) = meansoilm2
                    soilm3(t) = meansoilm3
                    soilm4(t) = meansoilm4
                endif
            enddo
        else
            ! If there are too many invalid tiles, we don't want to update, but
            ! we might need to adapt the states in case of perturbatoion bias
            ! correction.
            if(LIS_rc%pert_bias_corr.eq.1) then           
                ! if perturbation bias correction is turned on we readjust the ensemble

                bounds_violation = .true. 
                nIter = 0
                ens_flag = .true. 

                do while(bounds_violation) 
                    niter = niter + 1
                    t_unpert = i+LIS_rc%nensem(n)-1

                    ! calculate bias between perturbed and unperturbed
                    do j=1,4
                        bias(j) = 0.0
                        nmean = 0
                        do m=1,LIS_rc%nensem(n)-1
                            t = i+m-1
                            tmp = noahmp401_struc(n)%noahmp401(t)%sh2o(j)&
                                 - noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j)
                            ! uses Welford online mean, which should not be a
                            ! problem as long as nensem is small (<1000)
                            nmean = nmean + 1
                            bias(j) = bias(j) + (tmp - bias(j)) / nmean
                        enddo
                    enddo

                    do j=1,4
                        do m=1,LIS_rc%nensem(n)-1
                            t = i+m-1
                            SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype        
                            sm_threshold = SMCMAX_TABLE(SOILTYP) - 0.02  ! MAXSMC (SOILTYP) - 0.02

                            ! bias corrected value
                            tmp = noahmp401_struc(n)%noahmp401(t)%sh2o(j) - bias(j)

                            ! if the bias corrected value violates the threshold, we
                            ! set the current state to be within the bounds
                            if(tmp.le.MIN_THRESHOLD) then 
                                noahmp401_struc(n)%noahmp401(t)%sh2o(j) = &
                                     max(noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j),&
                                     MIN_THRESHOLD)
                                noahmp401_struc(n)%noahmp401(t)%smc(j) = &
                                     max(noahmp401_struc(n)%noahmp401(t_unpert)%smc(j),&
                                     MIN_THRESHOLD)
                                ens_flag(m) = .false. 
                            elseif(tmp.ge.sm_threshold) then
                                noahmp401_struc(n)%noahmp401(t)%sh2o(j) = &
                                     min(noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j),&
                                     sm_threshold)
                                noahmp401_struc(n)%noahmp401(t)%smc(j) = &
                                     min(noahmp401_struc(n)%noahmp401(t_unpert)%smc(j),&
                                     sm_threshold)
                                ens_flag(m) = .false. 
                            endif
                        enddo
                    enddo

                    !--------------------------------------------------------------------------
                    ! Recalculate the deltas and adjust the ensemble
                    !--------------------------------------------------------------------------
                    do j=1,4
                        bias(j) = 0.0
                        nmean = 0
                        do m=1,LIS_rc%nensem(n)-1
                            t = i+m-1
                            tmp = noahmp401_struc(n)%noahmp401(t)%sh2o(j)&
                                 - noahmp401_struc(n)%noahmp401(t_unpert)%sh2o(j)
                            ! uses Welford online mean, which should not be a
                            ! problem as long as nensem is small (<1000)
                            nmean = nmean + 1
                            bias(j) = bias(j) + (tmp - bias(j)) / nmean
                        enddo
                    enddo

                    do j=1,4
                        do m=1,LIS_rc%nensem(n)-1
                            t = i+m-1

                            SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype  
                            MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)           

                            if(ens_flag(m)) then 
                                tmp = noahmp401_struc(n)%noahmp401(t)%sh2o(j) - bias(j)

                                if(.not.(tmp.le.0.0 .or.&
                                     tmp.gt.(MAX_THRESHOLD))) then 

                                    noahmp401_struc(n)%noahmp401(t)%smc(j) = &
                                         noahmp401_struc(n)%noahmp401(t)%smc(j) - bias(j)
                                    noahmp401_struc(n)%noahmp401(t)%sh2o(j) = &
                                         noahmp401_struc(n)%noahmp401(t)%sh2o(j) - bias(j)
                                    bounds_violation = .false.
                                endif
                            endif

                            tmp = noahmp401_struc(n)%noahmp401(t)%sh2o(j)
                            if(tmp.le.0.0 .or.&
                                 tmp.gt.(MAX_THRESHOLD)) then 
                                bounds_violation = .true. 
                            else
                                bounds_violation = .false.
                            endif
                        enddo
                    enddo

                    if(nIter.gt.10.and.bounds_violation) then 
                        !--------------------------------------------------------------------------
                        ! All else fails, set to the bounds
                        !--------------------------------------------------------------------------

                        write(LIS_logunit,*) '[ERR] Ensemble structure violates physical bounds '
                        write(LIS_logunit,*) '[ERR] Please adjust the perturbation settings ..'

                        do j=1,4
                            do m=1,LIS_rc%nensem(n)
                                t = i+m-1
                                !t = (i-1)*LIS_rc%nensem(n)+m

                                SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype  
                                MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)           

                                if(noahmp401_struc(n)%noahmp401(t)%sh2o(j).gt.MAX_THRESHOLD.or.&
                                     noahmp401_struc(n)%noahmp401(t)%smc(j).gt.MAX_THRESHOLD) then                        
                                    noahmp401_struc(n)%noahmp401(t)%sh2o(j) = MAX_THRESHOLD
                                    noahmp401_struc(n)%noahmp401(t)%smc(j) = MAX_THRESHOLD
                                endif

                                if(noahmp401_struc(n)%noahmp401(t)%sh2o(j).lt.MIN_THRESHOLD.or.&
                                     noahmp401_struc(n)%noahmp401(t)%smc(j).lt.MIN_THRESHOLD) then                        
                                    noahmp401_struc(n)%noahmp401(t)%sh2o(j) = MIN_THRESHOLD
                                    noahmp401_struc(n)%noahmp401(t)%smc(j) = MIN_THRESHOLD
                                endif
                            enddo
                        enddo ! j
                    endif ! niter > 10
                enddo ! while bounds_violation
            endif !pert bias correction

            ! in the case of no updates we just set the updated values to the
            ! ones currently set
            do m=1,LIS_rc%nensem(n)
                t = i+m-1
                if (.not.update_flag_tile(t)) then
                    lai(t) = noahmp_struc(n)%noahmp401(t)%lai
                    soilm1(t) = noahmp_struc(n)%noahmp401(t)%smc(1)
                    soilm2(t) = noahmp_struc(n)%noahmp401(t)%smc(2)
                    soilm3(t) = noahmp_struc(n)%noahmp401(t)%smc(3)
                    soilm4(t) = noahmp_struc(n)%noahmp401(t)%smc(4)
                endif
            enddo
        endif !update_flags
    enddo ! patch loop

end subroutine NoahMP401_qclaisoilm

