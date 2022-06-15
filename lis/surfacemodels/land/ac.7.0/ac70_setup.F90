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
! !ROUTINE: Ac70_setup
! \label{Ac70_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and Ac70
!   2/7/18: Soni Yatheendradas; code added for OPTUE to work 
!
! !INTERFACE:
subroutine Ac70_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    !use AC_VEG_PARAMETERS_70, only: read_mp_veg_parameters
    !use MODULE_SF_ACLSM_70, only: read_mp_veg_parameters
    use MODULE_SF_ACLSM_70, only: read_mp_veg_parameters, &
           SLCATS, LUCATS, CSOIL_DATA, BB, SATDK, SATDW, &
           SATPSI, QTZ, MAXSMC, REFSMC, WLTSMC, &
           ! MB: AC70
           OC, WP, SAT, FC, INFRATE, SD, CL, SI, & 
           CZIL_DATA, FRZK_DATA, REFDK_DATA, REFKDT_DATA, SLOPE_DATA, &
           TOPT_DATA, RGLTBL, RSMAX_DATA, RSTBL, HSTBL, NROTBL, &
           CH2OP, DLEAF, Z0MVT, HVT, HVB, RC, RHOL, RHOS, TAUL, TAUS, &
           XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &         
           LTOVRC,  DILEFC,  DILEFW,  RMF25,  SLA,  FRAGR,  TMIN, &
           VCMX25,  TDLEF,  BP, MP, QE25, RMS25, RMR25, ARM, &
           FOLNMX, WDPOOL, WRRAT, MRP, DRYSMC ! SY: adding 
           ! calibratable parameters + DVEG to this list for OPTUE to work by 
           ! updating further below the corresponding values in Ac70_module 
           ! after call to subroutine SOIL_VEG_GEN_PARM_70   
           ! SY: Not used by Ac70 from REDPRM:  F11, DRYSMC            
           ! SY: Not used by Ac70 from read_mp_veg_parameters : AQE and 
           !     SLAREA (since subroutine BVOCFLUX not used), DEN
           ! SY: Also added DRYSMC used as constraint by OPTUE
 

    use Ac70_lsmMod

    !!! MB:
    use ac_global, only: typeproject_typeprm, &
                         typeproject_typepro, &
                         GetSimulParam_ThicknessTopSWC, &
                         GetRootZoneWC_Actual,&
                         GetRootZoneWC_FC,&
                         GetRootZoneWC_WP,&
                         GetRootZoneWC_SAT,&
                         GetRootZoneWC_Leaf,&
                         GetRootZoneWC_Thresh,&
                         GetRootZoneWC_Sen,&
                         GetRootZoneWC_ZtopAct,&
                         GetRootZoneWC_ZtopFC,&
                         GetRootZoneWC_ZtopWP,&
                         GetRootZoneWC_ZtopThresh,&
                         SetSimulParam_ThicknessTopSWC, &
                         SetRootZoneWC_Actual,&
                         SetRootZoneWC_FC,&
                         SetRootZoneWC_WP,&
                         SetRootZoneWC_SAT,&
                         SetRootZoneWC_Leaf,&
                         SetRootZoneWC_Thresh,&
                         SetRootZoneWC_Sen,&
                         SetRootZoneWC_ZtopAct,&
                         SetRootZoneWC_ZtopFC,&
                         SetRootZoneWC_ZtopWP,&
                         SetRootZoneWC_ZtopThresh,&
                         GetCompartment,&
                         SetCompartment,&
                         GetSoilLayer,&
                         SetSoilLayer,&
                         GetTotalSaltContent,&
                         GetTotalWaterContent,&
                         Geteffectiverain,&
                         GetSumWaBal,&
                         GetRootZoneSalt,&
                         GetSimulation,&
            GetSimulation_SumGDD,&
            GetSimulation_SumGDDfromDay1,&
                         SetTotalSaltContent,&
                         SetTotalWaterContent,&
                         Seteffectiverain,&
                         SetSumWaBal,&
                         SetRootZoneSalt,&
                         SetSimulation,&
                         GetIrrigation,&
                         SetIrrigation,&
                         GetCompartment_theta,&
                         SetCompartment_theta,&
                         GetIrriECw,&
                         GetManagement,&
                         GetPerennialPeriod,&
                         GetSimulParam,&
                         GetManagement_Cuttings,&
                         GetOnset,&
                         GetEndSeason,&
                         GetCrop,&
                         GetSoil,&
                         GetTemperatureRecord,&
                         GetClimRecord,&
                         GetRainRecord,&
                         GetEToRecord,&
                         SetIrriECw,&
                         SetManagement,&
                         SetPerennialPeriod,&
                         SetSimulParam,&
                         SetManagement_Cuttings,&
                         SetOnset,&
                         SetEndSeason,&
                         SetCrop,&
                         SetSoil,&
                         SetTemperatureRecord,&
                         SetClimRecord,&
                         SetRainRecord,&
                         SetEToRecord,&
                         GetSimulParam_GDDMethod,&

                         GetGenerateTimeMode,&
                        GetGenerateDepthMode,&
                        GetIrriMode,&
                        GetIrriMethod,&
                        GetDaySubmerged,&
                        GetMaxPlotNew,&
                        GetNrCompartments,&
                        GetIrriFirstDayNr,&
                        GetZiAqua,&
                        GetIniPercTAW,&
                        GetMaxPlotTr,&
                        GetOutputAggregate,&
                        GetEvapoEntireSoilSurface,&
                        GetPreDay,&
                        GetOutDaily,&
                        GetOut1Wabal,&
                        GetOut2Crop,&
                        GetOut3Prof,&
                        GetOut4Salt,&
                        GetOut5CompWC,&
                        GetOut6CompEC,&
                        GetOut7Clim,&
                        GetPart1Mult,&
                        GetPart2Eval,&
                        SetGenerateTimeMode,&
                        SetGenerateDepthMode,&
                        SetIrriMode,&
                        SetIrriMethod,&
                        SetDaySubmerged,&
                        SetMaxPlotNew,&
                        SetNrCompartments,&
                        SetIrriFirstDayNr,&
                        SetZiAqua,&
                        SetIniPercTAW,&
                        SetMaxPlotTr,&
                        SetOutputAggregate,&
                        SetEvapoEntireSoilSurface,&
                        SetPreDay,&
                        SetOutDaily,&
                        SetOut1Wabal,&
                        SetOut2Crop,&
                        SetOut3Prof,&
                        SetOut4Salt,&
                        SetOut5CompWC,&
                        SetOut6CompEC,&
                        SetOut7Clim,&
                        SetPart1Mult,&
                        SetPart2Eval,&


                        GetCCiActual,&
                        GetCCiprev,&
                        GetCCiTopEarlySen,&
                        GetCRsalt,&  
                        GetCRwater,& 
                        GetECdrain,& 
                        GetECiAqua,& 
                        GetECstorage,& 
                        GetEact,& 
                        GetEpot,& 
                        GetETo,&
                        GetDrain,&  
                        GetInfiltrated,&
                        GetRain,& 
                        GetRootingDepth,&
                        GetRunoff,& 
                        GetSaltInfiltr,&
                        GetSurf0,& 
                        GetSurfaceStorage,&
                        GetTact,&
                        GetTpot,&
                        GetTactWeedInfested,&
                        GetTmax,& 
                        GetTmin,& 
                        SetCCiActual,&
                        SetCCiprev,&
                        SetCCiTopEarlySen,&
                        SetCRsalt,&  
                        SetCRwater,& 
                        SetECdrain,& 
                        SetECiAqua,& 
                        SetECstorage,& 
                        SetEact,& 
                        SetEpot,& 
                        SetETo,&
                        SetDrain,&  
                        SetInfiltrated,&
                        SetRain,& 
                        SetRootingDepth,&
                        SetRunoff,& 
                        SetSaltInfiltr,&
                        SetSurf0,& 
                        SetSurfaceStorage,&
                        SetTact,&
                        SetTpot,&
                        SetTactWeedInfested,&
                        SetTmax,& 
                        SetTmin,&
                        GetIrriBeforeSeason,&
                        SetIrriBeforeSeason,&
                        GetIrriAfterSeason,&
                        SetIrriAfterSeason,&
                        GetCrop_Day1,&
                        DegreesDay,&
                         GetSimulation_ToDayNr, &
                         SetSimulation_ToDayNr, &
                        SetSimulation_SumGDDfromDay1,&
                        SetSimulation_SumGDD, &
                        SetSimulation_NrRuns, &
                        SetFullFileNameProgramParameters, &
                        SetSimulation_MultipleRun, &
                        GetSimulation_MultipleRunWithKeepSWC, &
                        GetSimulation_MultipleRunConstZrx, &
                        SetSimulation_MultipleRunConstZrx, &
                        CheckForKeepSWC, &
                        GetMultipleProjectFileFull, &
                        GetMultipleProjectFile, &
                        SetSimulation_MultipleRunWithKeepSWC, &
                        GetFullFileNameProgramParameters, &
                        SetMultipleProjectDescription, &
                        GetNumberSimulationRuns, &
                        SetMultipleProjectFile, &
                        SetMultipleProjectFileFull, &
                        SetPathNameOutp, &
                        SetPathNameSimul, &
                        SetPathNameList, &
                        SetPathNameParam, &
                        SetPathNameProg, &
                        GetPathNameList
           !!! MB_AC70

    use ac_run, only:    SetDayNri,&
                         GetIrriInterval,&
                         GetIrriInfoRecord1,&
                         GetIrriInfoRecord2,&
                         SetIrriInterval,&
                         SetIrriInfoRecord1,&
                         SetIrriInfoRecord2,&
                         GetTheProjectFile,&
                        GetGwTable,&
                        GetPlotVarCrop,&
                        GetStressTot,&
                        GetCutInfoRecord1,&
                        GetCutInfoRecord2,&
                        GetTransfer,&
                        GetPreviousSum,&
                        GetTadj,&
                        GetGDDTadj,&
                        GetDayLastCut,&
                        GetNrCut,&
                        GetSumInterval,&
                        GetPreviousStressLevel,&
                        GetStressSFadjNEW,&
                        GetBin,&
                        GetBout,&
                        GetGDDayi,&
                        GetCO2i,&
                        GetFracBiomassPotSF,&
                        GetSumETo,&
                        GetSumGDD,&
                        GetZiprev,&
                        GetSumGDDPrev,&
                        GetCCxWitheredTpot,&
                        GetCCxWitheredTpotNoS,&
                        GetCoeffb0,&
                        GetCoeffb1,&
                        GetCoeffb2,&
                        GetCoeffb0Salt,&
                        GetCoeffb1Salt,&
                        GetCoeffb2Salt,&
                        GetStressLeaf,&
                        GetStressSenescence ,&
                        GetDayFraction,&
                        GetGDDayFraction,&
                        GetCGCref,&
                        GetGDDCGCref ,&
                        GetTimeSenescence ,&
                        GetSumKcTop,&
                        GetSumKcTopStress,&
                        GetSumKci,&
                        GetCCoTotal,&
                        GetCCxTotal,&
                        GetCDCTotal,&
                        GetGDDCDCTotal,&
                        GetCCxCropWeedsNoSFstress,&
                        GetWeedRCi,&
                        GetCCiActualWeedInfested,&
                        GetfWeedNoS,&
                        GetZeval,&
                        GetBprevSum,&
                        GetYprevSum,&
                        GetSumGDDcuts,&
                        GetHItimesBEF,&
                        GetScorAT1,&
                        GetScorAT2,&
                        GetHItimesAT1,&
                        GetHItimesAT2,&
                        GetHItimesAT,&
                        GetalfaHI,&
                        GetalfaHIAdj,&
                        GetNextSimFromDayNr ,&
                        GetDayNr1Eval,&
                        GetDayNrEval,&
                        GetLineNrEval,&
                        GetPreviousSumETo,&
                        GetPreviousSumGDD,&
                        GetPreviousBmob,&
                        GetPreviousBsto,&
                        GetStageCode,&
                        GetPreviousDayNr,&
                        GetNoYear,&
                        GetWaterTableInProfile,&
                        GetStartMode,&
                        GetNoMoreCrop,&
                        GetCGCadjustmentAfterCutting,&

                        SetGwTable,&
                        SetPlotVarCrop,&
                        SetStressTot,&
                        SetCutInfoRecord1,&
                        SetCutInfoRecord2,&
                        SetTransfer,&
                        SetPreviousSum,&
                        SetTadj,&
                        SetGDDTadj,&
                        SetDayLastCut,&
                        SetNrCut,&
                        SetSumInterval,&
                        SetPreviousStressLevel,&
                        SetStressSFadjNEW,&
                        SetBin,&
                        SetBout,&
                        SetGDDayi,&
                        SetCO2i,&
                        SetFracBiomassPotSF,&
                        SetSumETo,&
                        SetSumGDD,&
                        SetZiprev,&
                        SetSumGDDPrev,&
                        SetCCxWitheredTpot,&
                        SetCCxWitheredTpotNoS,&
                        SetCoeffb0,&
                        SetCoeffb1,&
                        SetCoeffb2,&
                        SetCoeffb0Salt,&
                        SetCoeffb1Salt,&
                        SetCoeffb2Salt,&
                        SetStressLeaf,&
                        SetStressSenescence ,&
                        SetDayFraction,&
                        SetGDDayFraction,&
                        SetCGCref,&
                        SetGDDCGCref ,&
                        SetTimeSenescence ,&
                        SetSumKcTop,&
                        SetSumKcTopStress,&
                        SetSumKci,&
                        SetCCoTotal,&
                        SetCCxTotal,&
                        SetCDCTotal,&
                        SetGDDCDCTotal,&
                        SetCCxCropWeedsNoSFstress,&
                        SetWeedRCi,&
                        SetCCiActualWeedInfested,&
                        SetfWeedNoS,&
                        SetZeval,&
                        SetBprevSum,&
                        SetYprevSum,&
                        SetSumGDDcuts,&
                        SetHItimesBEF,&
                        SetScorAT1,&
                        SetScorAT2,&
                        SetHItimesAT1,&
                        SetHItimesAT2,&
                        SetHItimesAT,&
                        SetalfaHI,&
                        SetalfaHIAdj,&
                        SetNextSimFromDayNr ,&
                        SetDayNr1Eval,&
                        SetDayNrEval,&
                        SetLineNrEval,&
                        SetPreviousSumETo,&
                        SetPreviousSumGDD,&
                        SetPreviousBmob,&
                        SetPreviousBsto,&
                        SetStageCode,&
                        SetPreviousDayNr,&
                        SetNoYear,&
                        SetWaterTableInProfile,&
                        SetStartMode,&
                        SetNoMoreCrop,&
                        SetCGCadjustmentAfterCutting,&
                        AdvanceOneTimeStep, &
                        ReadClimateNextDay, &
                        SetGDDVariablesNextDay, &
                         FinalizeRun1, &
                         FinalizeRun2, &
                         GetDayNri,&
                         GetCrop_Tbase, &
                         GetCrop_Tupper, &
                         FinalizeSimulation, &
                         InitializeSimulation ,&
                         InitializeRunPart1 ,&
                         InitializeRunPart2 ,&
                         InitializeClimate,&
                         GetSimulation_ToDayNr, &
                         SetRain,& 
                         SetTmin,&
                         SetTmax,& 
                         SetETo,&
                         InitializeClimate, &
                         SetTheProjectFile

    use ac_kinds, only: intEnum, &
                        int32, &
                        int8, &
                        dp,&
                        sp
                         
    use ac_startunit, only:  GetListProjectsFile, &
                         GetNumberOfProjects, &
                         GetProjectFileName, &
                         GetProjectType, &
                         GetSimulation_NrRuns, &
                         InitializeTheProgram, &
                         InitializeProject, &
                         WriteProjectsInfo, &
                         GetTimeAggregationResults, &
                         GetRequestDailyResults, &
                         GetRequestParticularResults, &
                         !PrepareReport, &
                        LoadProgramParametersProjectPlugIn, &
                        ComposeFileForProgramParameters

    use ac_initialsettings, only: InitializeSettings

     !!! MB:
!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for Ac70.  These include: 
!    vegetype     - land cover type index [-]
!    soiltype     - soil type index [-]
!    slopetype    - slope type for Noah baseflow [-]
!    tbot         - deep-layer soil temperature [K]
!    pblh         - planetary boundary layer height [m]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data from NetCDF file
!  \item[AC70\_read\_MULTILEVEL\_param](\ref{AC70_read_MULTILEVEL_param}) \newline
!    retrieves MULTILEVEL spatial parameter from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t, k, n, l
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    ! MB
    real              :: Z_surf, cl_tmp, si_tmp, sd_tmp, InfRate_tmp
    integer           :: REW, descr
    
    !!! MB_AC70
    integer :: daynr, todaynr, iproject, nprojects, NrRuns
    integer(intEnum) :: TheProjectType
    logical :: ListProjectFileExist
    character(len=:), allocatable :: ListProjectsFile, TheProjectFile

    logical ::  ProgramParametersAvailable 
    integer(int32) :: TotalSimRuns
    character(len=:), allocatable :: FullFileNameProgramParametersLocal       
    logical :: MultipleRunWithKeepSWC_temp    
    real(dp) :: MultipleRunConstZrx_temp    
    !!! MB_AC70

    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder for #n nest
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        
        !------------------------------------!
        ! reading spatial spatial parameters !
        !------------------------------------!
        ! vegetype takes value from the LIS built-in parameter vegt
        !TODO: convert vegetation data source into vegetation types
        if(LIS_rc%uselcmap(n) .ne. 'none') then
            write(LIS_logunit,*) "Ac70: retrieve parameter VEGETYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                AC70_struc(n)%ac70(t)%vegetype= LIS_surface(n, mtype)%tile(t)%vegt
            enddo
        else 
            ! read: vegetype
            write(LIS_logunit,*) "Ac70: reading parameter VEGETYPE from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(AC70_struc(n)%LDT_ncvar_vegetype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AC70_struc(n)%ac70(t)%vegetype = placeholder(col, row)
            enddo 
        endif
        ! read: soiltype
        !write(LIS_logunit,*) "Ac70: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
        !call LIS_read_param(n, trim(AC70_struc(n)%LDT_ncvar_soiltype), placeholder)
        !do t = 1, LIS_rc%npatch(n, mtype)
        !    col = LIS_surface(n, mtype)%tile(t)%col
        !    row = LIS_surface(n, mtype)%tile(t)%row
        !    AC70_struc(n)%ac70(t)%soiltype = placeholder(col, row)
        !enddo 

        ! soiltype takes value from the LIS built-in parameter soilt
        !TODO: convert soil texture into soil types according to scheme
        if(LIS_rc%usetexturemap(n) .ne. 'none') then
            write(LIS_logunit,*) "Ac70: retrieve parameter SOILTYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                AC70_struc(n)%ac70(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
            enddo
        else 
            ! read: soiltype
            write(LIS_logunit,*) "Ac70: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(AC70_struc(n)%LDT_ncvar_soiltype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AC70_struc(n)%ac70(t)%soiltype = placeholder(col, row)
            enddo 
        endif

        ! read: slopetype
        write(LIS_logunit,*) "Ac70: reading parameter SLOPETYPE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AC70_struc(n)%LDT_ncvar_slopetype), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AC70_struc(n)%ac70(t)%slopetype = placeholder(col, row)
        enddo 

        ! read: tbot
        write(LIS_logunit,*) "Ac70: reading parameter TBOT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AC70_struc(n)%LDT_ncvar_tbot), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AC70_struc(n)%ac70(t)%tbot = placeholder(col, row)
        enddo 

        ! read: pblh
        write(LIS_logunit,*) "Ac70: reading parameter PBLH from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(AC70_struc(n)%LDT_ncvar_pblh), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            AC70_struc(n)%ac70(t)%pblh = placeholder(col, row)
        enddo 

        !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters !
        !----------------------------------------------!
        ! read: shdfac_monthly
        write(LIS_logunit,*) "Ac70: reading parameter SHDFAC_MONTHLY from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call AC70_read_MULTILEVEL_param(n, AC70_struc(n)%LDT_ncvar_shdfac_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                AC70_struc(n)%ac70(t)%shdfac_monthly(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: smceq for (opt_run=5)  Miguez-Macho & Fan groundwater with equilibrium water table
        if(AC70_struc(n)%run_opt .eq. 5) then
            write(LIS_logunit,*) "Ac70: reading parameter SMCEQ from ", trim(LIS_rc%paramfile(n))
            do k = 1, AC70_struc(n)%nsoil
                call AC70_read_MULTILEVEL_param(n, AC70_struc(n)%LDT_ncvar_smceq, k, placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    AC70_struc(n)%ac70(t)%smceq(k) = placeholder(col, row)
                enddo 
            enddo 
        endif

        deallocate(placeholder)
        call SOIL_VEG_GEN_PARM_70(AC70_struc(n)%landuse_tbl_name,   & 
                                  AC70_struc(n)%soil_tbl_name,      &
                                  AC70_struc(n)%gen_tbl_name,       &
                                  AC70_struc(n)%landuse_scheme_name,& 
                                  AC70_struc(n)%soil_scheme_name)
        ! SY: Begin for enabling OPTUE 
        do t = 1, LIS_rc%npatch(n, mtype)
            
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row

            ! SY: Begin lines following those in REDPRM

            ! SY: Begin SOIL PARAMETERS
            IF (AC70_struc(n)%ac70(t)%soiltype .gt. SLCATS) THEN
               write(LIS_logunit, *) 'SOILTYP must be less than SLCATS:'
               write(LIS_logunit, '("t = ", I6, "; SOILTYP = ", I6, ";    SLCATS = ", I6)') &
                         t, AC70_struc(n)%ac70(t)%soiltype, SLCATS
               write(LIS_logunit, *) 'Ac70_setup: Error: too many input soil types'
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            END IF
            AC70_struc(n)%ac70(t)%csoil = CSOIL_DATA
            AC70_struc(n)%ac70(t)%bexp = BB(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%dksat = SATDK(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%dwsat = SATDW(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%psisat = SATPSI(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%quartz = QTZ(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%smcmax = MAXSMC(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%smcref = REFSMC(AC70_struc(n)%ac70(t)%soiltype)
            AC70_struc(n)%ac70(t)%smcwlt = WLTSMC(AC70_struc(n)%ac70(t)%soiltype)
            ! SY: End SOIL PARAMETERS

            ! MB: AC70


            !!! Start the Program AC70
            ! Everything below is the equivalent of "StartTheProgram()"
            ! call InitializeTheProgram()

            call SetPathNameOutp(trim(AC70_struc(n)%PathNameOutp))
            call SetPathNameSimul(trim(AC70_struc(n)%PathNameSimul))
            call SetPathNameList(trim(AC70_struc(n)%PathNameList))
            call SetPathNameParam(trim(AC70_struc(n)%PathNameParam))
            call SetPathNameProg('')

            call GetTimeAggregationResults()
            call GetRequestDailyResults()
            call GetRequestParticularResults()
            !call PrepareReport()

            ListProjectsFile = GetListProjectsFile()
            inquire(file=trim(ListProjectsFile), exist=ListProjectFileExist)

            call WriteProjectsInfo('')
            call WriteProjectsInfo('Projects handled:')
            
            iproject = 1
            TheProjectFile = GetProjectFileName(iproject)
            call GetProjectType(TheProjectFile, TheProjectType)

        
            ! set AC70 soil parameters based on soiltype
            !AC70_struc(n)%ac70(t)%SoilLayer(1)%oc = OC(AC70_struc(n)%ac70(t)%soiltype) * 100
            AC70_struc(n)%ac70(t)%SoilLayer(1)%wp = WP(AC70_struc(n)%ac70(t)%soiltype) * 100
            AC70_struc(n)%ac70(t)%SoilLayer(1)%sat = SAT(AC70_struc(n)%ac70(t)%soiltype) * 100
            AC70_struc(n)%ac70(t)%SoilLayer(1)%fc = FC(AC70_struc(n)%ac70(t)%soiltype) * 100
            AC70_struc(n)%ac70(t)%SoilLayer(1)%InfRate = INFRATE(AC70_struc(n)%ac70(t)%soiltype) * 86400000
            sd_tmp = sd(AC70_struc(n)%ac70(t)%soiltype)
            cl_tmp = cl(AC70_struc(n)%ac70(t)%soiltype)
            si_tmp = si(AC70_struc(n)%ac70(t)%soiltype)
                
            ! define default CN
            if (AC70_struc(n)%ac70(t)%SoilLayer(1)%InfRate>864) then
                AC70_struc(n)%ac70(t)%Soil%CNvalue = 46
            elseif (AC70_struc(n)%ac70(t)%SoilLayer(1)%InfRate>=347) then
                AC70_struc(n)%ac70(t)%Soil%CNvalue = 61
            elseif (AC70_struc(n)%ac70(t)%SoilLayer(1)%InfRate>=36) then
                AC70_struc(n)%ac70(t)%Soil%CNvalue = 72
            else
                AC70_struc(n)%ac70(t)%Soil%CNvalue = 77
            endif
            
            ! define default REW (only top layer will be used)
            Z_surf = 0.04_dp
            REW=nint(10.0_dp*(AC70_struc(n)%ac70(t)%SoilLayer(1)%fc-(AC70_struc(n)%ac70(t)%SoilLayer(1)%wp/2.0))*Z_surf)
            if (REW < 0) REW = 0
            if (REW > 15) REW = 15
            AC70_struc(n)%ac70(t)%Soil%REW = REW

            !  associate soil class with USDA soil type for soil description
            if ((1.5 * cl_tmp + si_tmp) < 15.0) then
                descr = 0
            elseif (((1.5 * cl_tmp + si_tmp) >= 15) .and. ((2 * cl_tmp + si_tmp) <= 30)) then
                descr = 1
            elseif ((cl_tmp >= 7 .and. cl_tmp < 20 .and. sd_tmp >= 52) .and. (2 * cl_tmp + si_tmp >= 30)) then
                descr = 2
            elseif ((cl_tmp < 7) .and. (si_tmp < 50) .and. (2 * cl_tmp + si_tmp >= 30)) then
                descr = 2
            elseif ((cl_tmp >= 7) .and. (cl_tmp < 27) .and. (si_tmp >= 28) .and. (si_tmp < 50) .and. (sd_tmp < 52)) then
                descr = 3
            elseif ((si_tmp >= 50) .and. (cl_tmp >= 12) .and. (cl_tmp < 27)) then
                descr = 4
            elseif ((si_tmp >= 50) .and. (si_tmp < 80) .and. (cl_tmp < 12)) then
                descr = 4
            elseif ((si_tmp >= 80) .and. (cl_tmp < 12)) then
                descr = 5
            elseif ((cl_tmp >= 20) .and. (cl_tmp < 35) .and. (si_tmp < 28) .and. (sd_tmp > 45)) then
                descr = 6
            elseif ((cl_tmp >= 27) .and. (cl_tmp < 40) .and. (sd_tmp >= 20) .and. (sd_tmp < 45)) then
                descr = 7
            elseif ((cl_tmp >= 27) .and. (cl_tmp < 40) .and. (sd_tmp < 20)) then
                descr = 8
            elseif ((cl_tmp >= 35) .and. (cl_tmp < 55) .and. (sd_tmp >= 45) .and. (sd_tmp < 65)) then
                descr = 9
            elseif ((cl_tmp >= 40) .and. (si_tmp >= 40)) then
                descr = 10
            elseif ((cl_tmp >= 40) .and. (sd_tmp < 45) .and. (si_tmp < 40)) then
                descr = 11
            else
               write(LIS_logunit, *) 'no soil texture found'
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            end if
            ! soil_type = ['sand', 'loamy sand', 'sandy loam', 'loam', 'silt loam', 'silt', 'sandy clay loam',
            !              'clay loam',
            !              'silty clay loam', 'sandy clay', 'silty clay', 'clay']
            InfRate_tmp = AC70_struc(n)%ac70(t)%SoilLayer(1)%InfRate
            if ((descr == 0) .or. (descr == 1) .or. (descr == 2)) then
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRa = -0.3112-10**(-5)*InfRate_tmp
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRb = -1.4936+0.2416*log(InfRate_tmp)
            elseif ((descr == 3) .or. (descr == 4) .or. (descr == 5)) then
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRa = -0.4986-9*10**(-5)*InfRate_tmp
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRb = -2.1320+0.4778*log(InfRate_tmp)
            elseif ((descr == 6) .or. (descr == 7) .or. (descr == 9)) then
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRa = -0.5677-4*10**(-5)*InfRate_tmp
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRb = -3.7189+0.5922*log(InfRate_tmp)
            elseif ((descr == 8) .or. (descr == 10) .or. (descr == 11)) then
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRa = -0.6366+8*10**(-4)*InfRate_tmp
                AC70_struc(n)%ac70(t)%SoilLayer(1)%CRb = -1.9165+0.7063*log(InfRate_tmp)
            endif
            ! Set GravelMass and Penetrability
            AC70_struc(n)%ac70(t)%SoilLayer(1)%Penetrability = 100.0
            AC70_struc(n)%ac70(t)%SoilLayer(1)%GravelMass = 10.0
            AC70_struc(n)%ac70(t)%SoilLayer(1)%Description = 'soil type from LIS'
            AC70_struc(n)%ac70(t)%SoilLayer(2) = AC70_struc(n)%ac70(t)%SoilLayer(1)
            AC70_struc(n)%ac70(t)%SoilLayer(2)%GravelMass = 0.0

            AC70_struc(n)%ac70(t)%ProfDescription = 'soil profile from LIS'
            
            AC70_struc(n)%ac70(t)%SoilLayer(1)%Thickness = AC70_struc(n)%Thickness(1)
            AC70_struc(n)%ac70(t)%SoilLayer(2)%Thickness = AC70_struc(n)%Thickness(2)
            
            AC70_struc(n)%ac70(t)%Soil%NrSoilLayers = AC70_struc(n)%NrSoilLayers

            ! Set soil global variables 
            call SetSoilLayer(AC70_struc(n)%ac70(t)%SoilLayer)
            call SetSoil(AC70_struc(n)%ac70(t)%Soil)
            !!! Set all constant parameters from LIS_config
             

            !SetNumberSimulationRuns(AC70_struc(n)%NumberSimulationRuns)
            !call InitializeProject(iproject, TheProjectFile, TheProjectType)
            call InitializeSettings()
            
            ! paths need to be reset after Initialize Settings
            call SetPathNameOutp(trim(AC70_struc(n)%PathNameOutp))
            call SetPathNameSimul(trim(AC70_struc(n)%PathNameSimul))
            call SetPathNameList(trim(AC70_struc(n)%PathNameList))
            call SetPathNameParam(trim(AC70_struc(n)%PathNameParam))
            call SetPathNameProg('')
            !

            call SetMultipleProjectFile(TheProjectFile)
            call SetMultipleProjectFileFull(GetPathNameList() // &
                                    GetMultipleProjectFile())
            call GetNumberSimulationRuns(GetMultipleProjectFileFull(), &
                                                        TotalSimRuns)
            call SetMultipleProjectDescription('undefined')
            FullFileNameProgramParametersLocal = GetFullFileNameProgramParameters()
            call ComposeFileForProgramParameters(GetMultipleProjectFile(), &
                                      FullFileNameProgramParametersLocal)
            call SetFullFileNameProgramParameters(FullFileNameProgramParametersLocal)
            call LoadProgramParametersProjectPlugIn(&
                  GetFullFileNameProgramParameters(), &
                         ProgramParametersAvailable)
            call SetSimulation_MultipleRun(.true.)
            call SetSimulation_NrRuns(TotalSimRuns) 
            MultipleRunWithKeepSWC_temp = GetSimulation_MultipleRunWithKeepSWC() 
            MultipleRunConstZrx_temp = GetSimulation_MultipleRunConstZrx()
            call CheckForKeepSWC(GetMultipleProjectFileFull(), &  
                GetSimulation_NrRuns(), &       
                MultipleRunWithKeepSWC_temp, & 
                MultipleRunConstZrx_temp)  
            call SetSimulation_MultipleRunWithKeepSWC(MultipleRunWithKeepSWC_temp) 
            call SetSimulation_MultipleRunConstZrx(MultipleRunConstZrx_temp)        

            ! The remainder of this loop is the equivalent of "RunSimulation()"
            ! call InitializeSimulation(TheProjectFile, TheProjectType)
            call SetTheProjectFile(trim(TheProjectFile))

            ! Overwrite all AC70_struc after Initialization
            AC70_struc(n)%ac70(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
            AC70_struc(n)%ac70(t)%RootZoneWC_FC = GetRootZoneWC_FC()
            AC70_struc(n)%ac70(t)%RootZoneWC_WP = GetRootZoneWC_WP()
            AC70_struc(n)%ac70(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
            AC70_struc(n)%ac70(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
            AC70_struc(n)%ac70(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
            AC70_struc(n)%ac70(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
            AC70_struc(n)%ac70(t)%Compartment = GetCompartment()
            AC70_struc(n)%ac70(t)%TotalSaltContent = GetTotalSaltContent()
            AC70_struc(n)%ac70(t)%TotalWaterContent = GetTotalWaterContent()
            AC70_struc(n)%ac70(t)%effectiverain = Geteffectiverain()
            AC70_struc(n)%ac70(t)%SumWaBal = GetSumWaBal()
            AC70_struc(n)%ac70(t)%RootZoneSalt = GetRootZoneSalt()
            AC70_struc(n)%ac70(t)%Simulation = GetSimulation()
            AC70_struc(n)%ac70(t)%IrriInterval = GetIrriInterval()
            AC70_struc(n)%ac70(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
            AC70_struc(n)%ac70(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            AC70_struc(n)%ac70(t)%Irrigation = GetIrrigation()
            AC70_struc(n)%ac70(t)%IrriBeforeSeason = GetIrriBeforeSeason()
            AC70_struc(n)%ac70(t)%IrriAfterSeason = GetIrriAfterSeason()
            AC70_struc(n)%ac70(t)%SoilLayer = GetSoilLayer()
            AC70_struc(n)%ac70(t)%daynri = GetDayNri()
            do l=1, AC70_struc(n)%ac70(t)%NrCompartments
                 AC70_struc(n)%ac70(t)%smc(l) = GetCompartment_theta(l)
            enddo
            !write(*,'(e23.15e3)') AC70_struc(n)%ac70(t)%ac70smc(1)
            AC70_struc(n)%ac70(t)%IrriECw = GetIrriECw()
            AC70_struc(n)%ac70(t)%Management = GetManagement()
            AC70_struc(n)%ac70(t)%PerennialPeriod = GetPerennialPeriod()
            AC70_struc(n)%ac70(t)%simulparam = GetSimulParam()
            AC70_struc(n)%ac70(t)%Cuttings = GetManagement_Cuttings()
            AC70_struc(n)%ac70(t)%onset = GetOnset()
            AC70_struc(n)%ac70(t)%endseason = GetEndSeason()
            AC70_struc(n)%ac70(t)%crop = GetCrop()
            AC70_struc(n)%ac70(t)%Soil = GetSoil()
            AC70_struc(n)%ac70(t)%TemperatureRecord = GetTemperatureRecord()
            AC70_struc(n)%ac70(t)%ClimRecord = GetClimRecord()
            AC70_struc(n)%ac70(t)%RainRecord = GetRainRecord()
            AC70_struc(n)%ac70(t)%EToRecord = GetEToRecord()

            AC70_struc(n)%ac70(t)%GenerateTimeMode = GetGenerateTimeMode()
            AC70_struc(n)%ac70(t)%GenerateDepthMode = GetGenerateDepthMode()
            AC70_struc(n)%ac70(t)%IrriMode = GetIrriMode()
            AC70_struc(n)%ac70(t)%IrriMethod = GetIrriMethod()
            AC70_struc(n)%ac70(t)%DaySubmerged = GetDaySubmerged()
            AC70_struc(n)%ac70(t)%MaxPlotNew = GetMaxPlotNew()
            AC70_struc(n)%ac70(t)%NrCompartments = GetNrCompartments()
            AC70_struc(n)%ac70(t)%IrriFirstDayNr = GetIrriFirstDayNr()
            AC70_struc(n)%ac70(t)%ZiAqua = GetZiAqua()
            AC70_struc(n)%ac70(t)%IniPercTAW = GetIniPercTAW()
            AC70_struc(n)%ac70(t)%MaxPlotTr = GetMaxPlotTr()
            AC70_struc(n)%ac70(t)%OutputAggregate = GetOutputAggregate()

            AC70_struc(n)%ac70(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
            AC70_struc(n)%ac70(t)%PreDay = GetPreDay()
            AC70_struc(n)%ac70(t)%OutDaily = GetOutDaily()
            AC70_struc(n)%ac70(t)%Out1Wabal = GetOut1Wabal()
            AC70_struc(n)%ac70(t)%Out2Crop = GetOut2Crop()
            AC70_struc(n)%ac70(t)%Out3Prof = GetOut3Prof()
            AC70_struc(n)%ac70(t)%Out4Salt = GetOut4Salt()
            AC70_struc(n)%ac70(t)%Out5CompWC = GetOut5CompWC()
            AC70_struc(n)%ac70(t)%Out6CompEC = GetOut6CompEC()
            AC70_struc(n)%ac70(t)%Out7Clim = GetOut7Clim()
            AC70_struc(n)%ac70(t)%Part1Mult = GetPart1Mult()
            AC70_struc(n)%ac70(t)%Part2Eval = GetPart2Eval()

            !
            AC70_struc(n)%ac70(t)%CCiActual = GetCCiActual()
            AC70_struc(n)%ac70(t)%CCiprev = GetCCiprev()
            AC70_struc(n)%ac70(t)%CCiTopEarlySen = GetCCiTopEarlySen()
            AC70_struc(n)%ac70(t)%CRsalt = GetCRsalt () ! gram/m2
            AC70_struc(n)%ac70(t)%CRwater = GetCRwater() ! mm/day
            AC70_struc(n)%ac70(t)%ECdrain = GetECdrain() ! EC drain water dS/m
            AC70_struc(n)%ac70(t)%ECiAqua = GetECiAqua() ! EC of the groundwater table in dS/m
            AC70_struc(n)%ac70(t)%ECstorage = GetECstorage() !EC surface storage dS/m
            AC70_struc(n)%ac70(t)%Eact = GetEact() ! mm/day
            AC70_struc(n)%ac70(t)%Epot = GetEpot() ! mm/day
            AC70_struc(n)%ac70(t)%ETo_ac = GetETo() ! mm/day
            AC70_struc(n)%ac70(t)%Drain = GetDrain()  ! mm/day
            AC70_struc(n)%ac70(t)%Infiltrated = GetInfiltrated() ! mm/day
            AC70_struc(n)%ac70(t)%PREC_ac = GetRain()  ! mm/day
            AC70_struc(n)%ac70(t)%RootingDepth = GetRootingDepth()
            AC70_struc(n)%ac70(t)%Runoff = GetRunoff()  ! mm/day
            AC70_struc(n)%ac70(t)%SaltInfiltr = GetSaltInfiltr() ! salt infiltrated in soil profile Mg/ha
            AC70_struc(n)%ac70(t)%Surf0 = GetSurf0()  ! surface water [mm] begin day
            AC70_struc(n)%ac70(t)%SurfaceStorage = GetSurfaceStorage() !mm/day
            AC70_struc(n)%ac70(t)%Tact = GetTact() ! mm/day
            AC70_struc(n)%ac70(t)%Tpot = GetTpot() ! mm/day
            AC70_struc(n)%ac70(t)%TactWeedInfested = GetTactWeedInfested() !mm/day
            AC70_struc(n)%ac70(t)%Tmax_ac = GetTmax() ! degC
            AC70_struc(n)%ac70(t)%Tmin_ac =GetTmin() ! degC


            AC70_struc(n)%ac70(t)%GwTable = GetGwTable()
            AC70_struc(n)%ac70(t)%PlotVarCrop = GetPlotVarCrop()
            AC70_struc(n)%ac70(t)%StressTot = GetStressTot()
            AC70_struc(n)%ac70(t)%CutInfoRecord1 = GetCutInfoRecord1()
            AC70_struc(n)%ac70(t)%CutInfoRecord2 = GetCutInfoRecord2()
            AC70_struc(n)%ac70(t)%Transfer = GetTransfer()
            AC70_struc(n)%ac70(t)%PreviousSum = GetPreviousSum()
            AC70_struc(n)%ac70(t)%Tadj = GetTadj()
            AC70_struc(n)%ac70(t)%GDDTadj = GetGDDTadj()
            AC70_struc(n)%ac70(t)%DayLastCut = GetDayLastCut()
            AC70_struc(n)%ac70(t)%NrCut = GetNrCut()
            AC70_struc(n)%ac70(t)%SumInterval = GetSumInterval()
            AC70_struc(n)%ac70(t)%PreviousStressLevel = GetPreviousStressLevel()
            AC70_struc(n)%ac70(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC70_struc(n)%ac70(t)%Bin = GetBin()
            AC70_struc(n)%ac70(t)%Bout = GetBout()
            AC70_struc(n)%ac70(t)%GDDayi = GetGDDayi()
            AC70_struc(n)%ac70(t)%CO2i = GetCO2i()
            AC70_struc(n)%ac70(t)%FracBiomassPotSF = GetFracBiomassPotSF()
            AC70_struc(n)%ac70(t)%SumETo = GetSumETo()
            AC70_struc(n)%ac70(t)%SumGDD = GetSumGDD()
            AC70_struc(n)%ac70(t)%Ziprev = GetZiprev()
            AC70_struc(n)%ac70(t)%SumGDDPrev = GetSumGDDPrev()
            AC70_struc(n)%ac70(t)%CCxWitheredTpot = GetCCxWitheredTpot()
            AC70_struc(n)%ac70(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
            AC70_struc(n)%ac70(t)%Coeffb0 = GetCoeffb0()
            AC70_struc(n)%ac70(t)%Coeffb1 = GetCoeffb1()
            AC70_struc(n)%ac70(t)%Coeffb2 = GetCoeffb2()
            AC70_struc(n)%ac70(t)%Coeffb0Salt = GetCoeffb0Salt()
            AC70_struc(n)%ac70(t)%Coeffb1Salt = GetCoeffb1Salt()
            AC70_struc(n)%ac70(t)%Coeffb2Salt = GetCoeffb2Salt()
            AC70_struc(n)%ac70(t)%StressLeaf = GetStressLeaf()
            AC70_struc(n)%ac70(t)%StressSenescence = GetStressSenescence ()
            AC70_struc(n)%ac70(t)%DayFraction = GetDayFraction()
            AC70_struc(n)%ac70(t)%GDDayFraction = GetGDDayFraction()
            AC70_struc(n)%ac70(t)%CGCref = GetCGCref()
            AC70_struc(n)%ac70(t)%GDDCGCref = GetGDDCGCref ()
            AC70_struc(n)%ac70(t)%TimeSenescence = GetTimeSenescence ()
            AC70_struc(n)%ac70(t)%SumKcTop = GetSumKcTop()
            AC70_struc(n)%ac70(t)%SumKcTopStress = GetSumKcTopStress()
            AC70_struc(n)%ac70(t)%SumKci = GetSumKci()
            AC70_struc(n)%ac70(t)%CCoTotal = GetCCoTotal()
            AC70_struc(n)%ac70(t)%CCxTotal = GetCCxTotal()
            AC70_struc(n)%ac70(t)%CDCTotal = GetCDCTotal()
            AC70_struc(n)%ac70(t)%GDDCDCTotal = GetGDDCDCTotal()
            AC70_struc(n)%ac70(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
            AC70_struc(n)%ac70(t)%WeedRCi = GetWeedRCi()
            AC70_struc(n)%ac70(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
            AC70_struc(n)%ac70(t)%fWeedNoS = GetfWeedNoS()
            AC70_struc(n)%ac70(t)%Zeval = GetZeval()
            AC70_struc(n)%ac70(t)%BprevSum = GetBprevSum()
            AC70_struc(n)%ac70(t)%YprevSum = GetYprevSum()
            AC70_struc(n)%ac70(t)%SumGDDcuts = GetSumGDDcuts()
            AC70_struc(n)%ac70(t)%HItimesBEF = GetHItimesBEF()
            AC70_struc(n)%ac70(t)%ScorAT1 = GetScorAT1()
            AC70_struc(n)%ac70(t)%ScorAT2 = GetScorAT2()
            AC70_struc(n)%ac70(t)%HItimesAT1 = GetHItimesAT1()
            AC70_struc(n)%ac70(t)%HItimesAT2 = GetHItimesAT2()
            AC70_struc(n)%ac70(t)%HItimesAT = GetHItimesAT()
            AC70_struc(n)%ac70(t)%alfaHI = GetalfaHI()
            AC70_struc(n)%ac70(t)%alfaHIAdj = GetalfaHIAdj()
            AC70_struc(n)%ac70(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
            AC70_struc(n)%ac70(t)%DayNr1Eval = GetDayNr1Eval()
            AC70_struc(n)%ac70(t)%DayNrEval = GetDayNrEval()
            AC70_struc(n)%ac70(t)%LineNrEval = GetLineNrEval()
            AC70_struc(n)%ac70(t)%PreviousSumETo = GetPreviousSumETo()
            AC70_struc(n)%ac70(t)%PreviousSumGDD = GetPreviousSumGDD()
            AC70_struc(n)%ac70(t)%PreviousBmob = GetPreviousBmob()
            AC70_struc(n)%ac70(t)%PreviousBsto = GetPreviousBsto()
            AC70_struc(n)%ac70(t)%StageCode = GetStageCode()
            AC70_struc(n)%ac70(t)%PreviousDayNr = GetPreviousDayNr()
            AC70_struc(n)%ac70(t)%NoYear = GetNoYear()
            AC70_struc(n)%ac70(t)%WaterTableInProfile = GetWaterTableInProfile()
            AC70_struc(n)%ac70(t)%StartMode = GetStartMode()
            AC70_struc(n)%ac70(t)%NoMoreCrop = GetNoMoreCrop()
            AC70_struc(n)%ac70(t)%CGCadjustmentAfterCutting = GetCGCadjustmentAfterCutting()
                

            ! 
            !!! MB_AC70
            !call GetTimeAggregationResults()
            !call GetRequestDailyResults()
            !call GetRequestParticularResults()
            !call PrepareReport()

            !ListProjectsFile = GetListProjectsFile()
            !inquire(file=trim(ListProjectsFile), exist=ListProjectFileExist)
            !nprojects = GetNumberOfProjects()

            !if (nprojects > 0) then
            !    call WriteProjectsInfo('')
            !    call WriteProjectsInfo('Projects handled:')
            !end if
            
            !do iproject = 1, nprojects
                !TheProjectFile = GetProjectFileName(iproject)
                !TheProjectFile = GetProjectFileName(1)
                !call GetProjectType(TheProjectFile, TheProjectType)
 !               TheProjectType = 1 ! for LIS always set to multiple years
 !               call InitializeProject(iproject, TheProjectFile, TheProjectType)
            
                ! The remainder of this loop is the equivalent of "RunSimulation()"
 !               call InitializeSimulation(TheProjectFile, TheProjectType)
                ! 
            !end do
            !!! MB_AC70


            ! MB:

            ! SY: Begin SOIL PARAMETER CONSTRAINT
            AC70_struc(n)%ac70(t)%smcdry = DRYSMC(AC70_struc(n)%ac70(t)%soiltype)
            ! SY: End SOIL PARAMETER CONSTRAINT

            ! SY: Begin UNIVERSAL PARAMETERS
            AC70_struc(n)%ac70(t)%czil = CZIL_DATA
            AC70_struc(n)%ac70(t)%frzk = FRZK_DATA
            AC70_struc(n)%ac70(t)%refdk = REFDK_DATA
            AC70_struc(n)%ac70(t)%refkdt = REFKDT_DATA
            AC70_struc(n)%ac70(t)%slope = SLOPE_DATA(AC70_struc(n)%ac70(t)%slopetype)
            ! SY: End UNIVERSAL PARAMETERS

            ! SY: Begin VEGETATION PARAMETERS
            IF (AC70_struc(n)%ac70(t)%vegetype .gt. LUCATS) THEN
               write(LIS_logunit, *) 'VEGTYP must be less than LUCATS:'
               write(LIS_logunit, '("t = ", I6, "; VEGTYP = ", I6, ";    LUCATS = ", I6)') &
                         t, AC70_struc(n)%ac70(t)%vegetype, LUCATS
               write(LIS_logunit, *) 'Ac70_setup: Error: too many input landuse types'
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            END IF
            AC70_struc(n)%ac70(t)%topt = TOPT_DATA
            AC70_struc(n)%ac70(t)%rgl = RGLTBL(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%rsmax = RSMAX_DATA
            AC70_struc(n)%ac70(t)%rsmin = RSTBL(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%hs = HSTBL(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%nroot = & 
                 NROTBL(AC70_struc(n)%ac70(t)%vegetype)

            IF (NROTBL(AC70_struc(n)%ac70(t)%vegetype) .gt. &
                AC70_struc(n)%nsoil) THEN
               WRITE (LIS_logunit,*) 'Warning: too many root layers'
               write (LIS_logunit,*) 'NROOT = ', NROTBL(AC70_struc(n)%ac70(t)%vegetype)
               write (LIS_logunit,*) 'NSOIL = ', AC70_struc(n)%nsoil
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            END IF
            ! SY: End VEGETATION PARAMETERS

            ! SY: End lines following those in REDPRM

        enddo ! do t = 1, LIS_rc%npatch(n, mtype)
        ! SY: End for enabling OPTUE 

        call read_mp_veg_parameters(trim(AC70_struc(n)%ac_tbl_name), &
                                    trim(AC70_struc(n)%landuse_scheme_name))

        ! SY: Begin for enabling OPTUE 
        do t = 1, LIS_rc%npatch(n, mtype)
            
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row

            ! SY: Begin lines following those in read_mp_veg_parameters
            AC70_struc(n)%ac70(t)%CH2OP = CH2OP(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%DLEAF = DLEAF(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%Z0MVT = Z0MVT(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%HVT = HVT(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%HVB = HVB(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%RC = RC(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%RHOL1 = RHOL(AC70_struc(n)%ac70(t)%vegetype,1)
            AC70_struc(n)%ac70(t)%RHOL2 = RHOL(AC70_struc(n)%ac70(t)%vegetype,2)
            AC70_struc(n)%ac70(t)%RHOS1 = RHOS(AC70_struc(n)%ac70(t)%vegetype,1)
            AC70_struc(n)%ac70(t)%RHOS2 = RHOS(AC70_struc(n)%ac70(t)%vegetype,2)
            AC70_struc(n)%ac70(t)%TAUL1 = TAUL(AC70_struc(n)%ac70(t)%vegetype,1)
            AC70_struc(n)%ac70(t)%TAUL2 = TAUL(AC70_struc(n)%ac70(t)%vegetype,2)
            AC70_struc(n)%ac70(t)%TAUS1 = TAUS(AC70_struc(n)%ac70(t)%vegetype,1)
            AC70_struc(n)%ac70(t)%TAUS2 = TAUS(AC70_struc(n)%ac70(t)%vegetype,2)
            AC70_struc(n)%ac70(t)%XL = XL(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%CWPVT = CWPVT(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%C3PSN = C3PSN(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%KC25 = KC25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%AKC = AKC(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%KO25 = KO25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%AKO = AKO(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%AVCMX = AVCMX(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%AQE = AQE(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%LTOVRC = LTOVRC(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%DILEFC = DILEFC(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%DILEFW = DILEFW(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%RMF25 = RMF25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%SLA = SLA(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%FRAGR = FRAGR(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%TMIN = TMIN(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%VCMX25 = VCMX25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%TDLEF = TDLEF(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%BP = BP(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%MP = MP(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%QE25 = QE25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%RMS25 = RMS25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%RMR25 = RMR25(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%ARM = ARM(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%FOLNMX = FOLNMX(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%WDPOOL = WDPOOL(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%WRRAT = WRRAT(AC70_struc(n)%ac70(t)%vegetype)
            AC70_struc(n)%ac70(t)%MRP = MRP(AC70_struc(n)%ac70(t)%vegetype)
            ! SY: End lines following those in read_mp_veg_parameters
            
            ! MB:
            AC70_struc(n)%ac70(t)%NrRuns = NrRuns
            AC70_struc(n)%ac70(t)%TheProjectType = TheProjectType

        enddo ! do t = 1, LIS_rc%npatch(n, mtype)
        ! SY: End for enabling OPTUE 
    enddo

end subroutine Ac70_setup

!BOP
!
! !ROUTINE: AC70_read_MULTILEVEL_param
!  \label{AC70_read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
!
! !INTERFACE:
subroutine AC70_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
! !USES:
    use netcdf
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet,   &   
                            LIS_ews_halo_ind, LIS_ewe_halo_ind, &
                            LIS_nss_halo_ind, LIS_nse_halo_ind   
    use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    implicit none
! !ARGUMENTS: 
    integer, intent(in)          :: n
    integer, intent(in)          :: level
    character(len=*), intent(in) :: ncvar_name 
    real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine reads MULTILEVEL parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[level]
!    level index (month, quarter, soil layer, snow layer) of the data to be read
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

    integer       :: ios1
    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, t, nlevel, k
    real, pointer :: level_data(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in AC70_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AC70_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in AC70_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC70_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in AC70_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in AC70_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in AC70_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
 end subroutine AC70_read_MULTILEVEL_param
                                          

