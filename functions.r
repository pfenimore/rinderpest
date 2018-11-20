#Rinderpest public functions.r

#  Read in GIS data 
initGeog <- function(wbb, HistCatRed, type) {  # type is not used here, but is a place holder for animal type.
  print("In initGeog")
  if (is.null(wbb))
    return(NULL)
  if (type == "cattle") {
    rr <- crop(cattle, wbb)/HistCatRed        # cattle is raster data - see global.r
    rr <- rr*area(rr)  # Change from density to number  area is the area of each cell in km^2
    # Note about area function: area is computed as the height (latitudial span) of a cell (which is constant among all cells) times the width 
    # (longitudinal span) in the (latitudinal) middle of a cell. The width is smaller at the poleward side than at 
    # the equator-ward side of a cell. This variation is greatest near the poles and the values are thus not very 
    # precise for very high latitudes.
  }  

return(rr)
}

# Initialize the population matrix, N 
initN <- function(Ndstates,Ndcomb,Nages,r,rVacc){
  Ns = values(r)
  nbox=length(r)
  N = array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox))

  TotalPop <- sum(Ns)
  N[1,1,,] <- unlist(Ns)
  N <- TotalPop*N/sum(N)
  
  # Assumes 97% effectiveness of vaccination
  N[9,1,1,] <- N[1,1,1,]*getValues(rVacc)/100*0.97
  N[1,1,1,] <- N[1,1,1,]*(100-getValues(rVacc)*0.97)/100
  return(N)
}

kValues <- function(dd.weight_totI, dd.weight_totHplus, betaRatio, N) {
  Ib <- (dd.weight_totI %*% N[3,1,1,] + dd.weight_totHplus %*% N[4,1,1,])  # dd.weights contain beta!!  
  totalN <- colSums(N[1:5,1,1,])+colSums(N[8:9,1,1,])
  k <-  Ib/totalN #beta*I/N
  k[which(totalN <= 1)] <- Ib[which(totalN <= 1)]                     # Discuss with Carrie?
  return(k)
}

#================================================================ NextN ========================================================
NextN <- function(N, ApproxPixelArea, MM, tstep, klist, Ks, CellsToVacc, CellsToCull, CullFrom) {
  Ndstates <- dim(N)[1]
  nages=dim(N)[3]
  Nboxes <- dim(N)[4]
  
  totalStart <- colSums(N, dims = 3)
  # print("Box 89 at start")
  # print(N[,1,1,89])
  
  dd.weight_totI <- MM$dd.weight_totI
  dd.weight_totHplus <- MM$dd.weight_totHplus
  
  ################ Read in parameters from klist and Ks #############
  kEIvalue   <- klist$kEI   #   kEI   <- kD1[2]
  kIH   <- Ks$kIH           #   kIH   <- kD1[3]
  kIR   <- Ks$kIR           #   kIR   <- (kD1[3]+Rr1[3]*kD1[3])*Rr1[2]
  kHHt  <- Ks$kHHt          #   kHHt  <- Rr1[4]*kD1[4]
  kHD   <- klist$kHD        #   kHD   <- kD1[4]
  kHR   <- klist$kHR        #   kHD   <- kD1[4]
  kHtR  <- Ks$kHtR          #   kHtR  <- Rr1[6]*kD1[6]
  kHtD  <- Ks$kHtD          #   kHtD  <- kD1[6]
  ke    <- klist$ke

  K1    <- Ks$K1  #   K1    <- kIH + kIIt + kIR   Can't physically be 0.
  K2    <- Ks$K2  #   K2    <- kHHt + kHD         Can't physically be 0.
  K4    <- Ks$K4  #   K4    <- kHtD + kHtR        Could be 0.
  
  Zs   <- rep(0,Nboxes) # zeros
  Os   <- rep(1,Nboxes) # ones
  
  VaccPerStep <- Ks$VaccPerStep
  Crate <- Ks$Crate
  betaRatio <- Ks$betaRatio
  
  InitVals <- InitialValues(N)  # In FunctionsSame.r
  
  t <- tstep
  CutOff <- animalCut/26*ApproxPixelArea
  
  kEI <- kEIValues(kEIvalue, InitVals, CutOff, Nboxes) # In FunctionsSame.r
  
  k <- kValues(dd.weight_totI, dd.weight_totHplus, betaRatio, N)
  
  CullRateFunc <- function(N, CellsToCull, Crate, CullFrom, tstep, CutOff) {
    kCS <- N[1,1,1,]*0  
    kCE <- N[1,1,1,]*0
    kCI <- N[1,1,1,]*0
    kCH <- N[1,1,1,]*0
    if (!is.null(CellsToCull)) {
      CullPerStep <- Crate*tstep  # Crate is number per day
      if (CullFrom == "Honly") totalToCull <- sum(N[4,1,1,CellsToCull])
      if (CullFrom == "IandH") totalToCull <- sum(N[3:4,1,1,CellsToCull])
      if (CullFrom == "All")  totalToCull <- sum(N[1:4,1,1,CellsToCull])
      Ccheck_rate <- N[1,1,1,]*NaN
      Ccheck_rate[CellsToCull] <- 1-CullPerStep/totalToCull
      #print(CullPerStep)
      #print(totalToCull)
      kCH[which(Ccheck_rate > 0)] <- -log(Ccheck_rate[which(Ccheck_rate>0)])/tstep  
    
      if (CullFrom == "IandH")  kCI <- kCH
      if (CullFrom == "All") {
        kCI <- kCH
        kCS <- kCH
        kCE <- kCH
      }
      
      if (CullFrom == "Honly") {
        kCH[which(Ccheck_rate <= 0)] <- -log(CutOff/N[4,1,1,which(Ccheck_rate<=0)])/tstep
        kCH[which(N[4,1,1,] < CutOff)] <- 0
      }
      
      if (CullFrom == "IandH") {
        if (length(which(Ccheck_rate<=0)) == 1) {  # Use sum not colSums
          kCH[which(Ccheck_rate <= 0)] <- -log(CutOff/sum(N[3:4,1,1,which(Ccheck_rate<=0)]))/tstep
        } else {
          kCH[which(Ccheck_rate <= 0)] <- -log(CutOff/colSums(N[3:4,1,1,which(Ccheck_rate<=0)]))/tstep
        }  
        kCH[which(N[4,1,1,] < CutOff)] <- 0
        kCI[which(Ccheck_rate <= 0)] <- kCH[which(Ccheck_rate <= 0)]
      }
      if (CullFrom == "All")   {
        cull_all_indices <- which(Ccheck_rate <= 0)
        #print(cull_all_indices)
        #print(colSums(N[1:4,1,1,cull_all_indices]))
        if (length(cull_all_indices) == 1) {
          #kCH[cull_all_indices] <- -log(CutOff/sum(N[1:4,1,1,cull_all_indices]))/tstep
          kCH[cull_all_indices] <- -log(0.0001/sum(N[1:4,1,1,cull_all_indices]))/tstep
          kCH[which(kCH<0)] <- 0
        }
        else {
          #kCH[cull_all_indices] <- -log(CutOff/colSums(N[1:4,1,1,cull_all_indices]))/tstep
          kCH[cull_all_indices] <- -log(0.0001/colSums(N[1:4,1,1,cull_all_indices]))/tstep
          kCH[which(kCH<0)] <- 0
          #print(kCH[cull_all_indices])
        }
        kCI[which(Ccheck_rate <= 0)] <- kCH[which(Ccheck_rate <= 0)]
        kCE[which(Ccheck_rate <= 0)] <- kCH[which(Ccheck_rate <= 0)]
        kCS[which(Ccheck_rate <= 0)] <- kCH[which(Ccheck_rate <= 0)]
      }  
    }
    CullRates <- list(kCS=kCS,kCE=kCE,kCI=kCI,kCH=kCH)
    return(CullRates)
  }  
  CullRates <- CullRateFunc(N, CellsToCull, Crate, CullFrom, tstep, CutOff/10)
  kCS <- CullRates$kCS
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  VaccRates <- function(N, CellsToVacc, VaccPerStep, tstep, CutOff) {
    kV <- N[1,1,1,]*0           # Initializing
    check_rate <- N[1,1,1,]*NaN
    if (!is.null(CellsToVacc)) {
      #Vacc_for_each_cell <- VaccPerStep/sum(Ng[1,CellsToVacc])*Ng[1,CellsToVacc]  #Derivation
      #check_rate[CellsToVacc] <- 1-Vacc_for_each_cell/Ng[1,CellsToVacc]
      if(sum(N[1,1,1,CellsToVacc]) > 0)  { 
        check_rate[CellsToVacc] <- 1-VaccPerStep/sum(N[1,1,1,CellsToVacc])
        kV[which(check_rate > 0)] <- -log(check_rate[which(check_rate>0)])/t  
        kV[which(check_rate <= 0)] <- -log(CutOff/N[1,1,1,which(check_rate<=0)])/t  
        kV[which(N[1,1,1,] < CutOff)] <- 0
      }
    }
    return(kV)
  }
  kV <- VaccRates(N, CellsToVacc, VaccPerStep, tstep, CutOff)
  
  # Choose boxes: Vacc or No Vacc  
    v   <- which(kV>0)
    nv  <- which(kV==0)
    
  CullRates <- AvoidRedundantEigenValues(CullRates, K1, K2, K4, kEIvalue, ke, kV)  # In FunctionsSame.r
    
  kIP <- kIndepParams(kEI, kIH, kHHt, CullRates, K1, K2, K4, tstep) # In FunctionsSame.r
  kDP <- kDepParams(kEI, kIH, kHHt, k, ke, kV, CullRates, K1, K2, K4) # In FunctionsSame.r

  E1N <- EigenVector1NoVacc(Zs, Os, k, kEI, kIH, kIR, kHD, kHR, kHHt, kHtR, kHtD, K1, K2, K4, CullRates)
  if (length(v) > 0) {
    E1 <- EigenVector1wVacc(Zs, Os, k, kDP$fkV, ke, kDP$Ka, kDP$fK1, kDP$fK2, kDP$fK3, kDP$fK5, CullRates, kHD, kHtD, kIR, kHR, kHtR)
    E2 <- EigenVector2wVacc(Zs, Os, k, kEI, kIH, kIR, kHHt, kHD, kHR, kHtD, kHtR, CullRates, ke, kDP$ks, K1, K2, K4)
  }
  # Eigenvector 3 - a place holder for when treatment in It is in place
  E4 <- EigenVector4(Zs, kEI, kIP$fkEI1, kIP$fkEI2, kIP$fkEI4, CullRates, kHD, kHtD, kIR, kHR, kHtR)
  E5 <- EigenVector5(Os, kIP$fK1a, kIP$fK1e, K1, CullRates, kHD, kHtD, kIR, kHR, kHtR)
  E6 <- EigenVector6(Os, kIP$fK4d, K2, CullRates, kHD, kHtD, kHR, kHtR)  
  E7 <- EigenVector7(Os, K4, kHtD, kHtR, Nboxes)
  
  Ntemp <- N
  if (length(v) > 0) {
    Nv <- N_New_wVacc(v, kDP$Ka, kDP$ks, N, Ndstates, Nboxes, E1, E2, E3, E4, E5, E6, E7, InitVals, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4, tstep)
    Ntemp[,,,v] <- Nv[,,,v]
  }
  Nnv <- N_NewNoVacc(nv, N, Ndstates, Nboxes, E1N, E4, E5, E6, E7, k, CullRates, tstep, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4, InitVals, ke)
  Ntemp[,,,nv] <- Nnv[,,,nv]

  k2 <- kValues(dd.weight_totI, dd.weight_totHplus, betaRatio, Ntemp)
  k <- (k+k2)/2
  kDP <- kDepParams(kEI, kIH, kHHt, k, ke, kV, CullRates, K1, K2, K4)

  E1N <- EigenVector1NoVacc(Zs, Os, k, kEI, kIH, kIR, kHD, kHR, kHHt, kHtR, kHtD, K1, K2, K4, CullRates)
  if (length(v) > 0) {
    E1 <- EigenVector1wVacc(Zs, Os, k, kDP$fkV, ke, kDP$Ka, kDP$fK1, kDP$fK2, kDP$fK3, kDP$fK5, CullRates, kHD, kHtD, kIR, kHR, kHtR)
    E2 <- EigenVector2wVacc(Zs, Os, k, kEI, kIH, kIR, kHHt, kHD, kHR, kHtD, kHtR, CullRates, ke, kDP$ks, K1, K2, K4)
    Nv <- N_New_wVacc(v, kDP$Ka, kDP$ks, N, Ndstates, Nboxes, E1, E2, E3, E4, E5, E6, E7, InitVals, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4, tstep)
    N[,,,v] <- Nv[,,,v]
  }
  Nnv <- N_NewNoVacc(nv, N, Ndstates, Nboxes, E1N, E4, E5, E6, E7, k, CullRates, tstep, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4, InitVals, ke)
  N[,,,nv] <- Nnv[,,,nv]
  
  num_acc_ind <- which(N < 0, arr.ind=TRUE)
  if (length(num_acc_ind) > 0)  {
    print("1. values of boxes with problems")
    print(N[num_acc_ind])
    print(num_acc_ind[,4])
    print(kCS[num_acc_ind[,4]])
  }
  N[num_acc_ind] = 0
  
  totalEnd <- colSums(N, dims = 3)
  
  if (any(abs(totalStart-totalEnd) > 0.001)) {
    print("population not conserved")
    badboxes <- which(abs(totalStart-totalEnd) > 0.001)
    print(paste(intersect(120,nv), totalStart[120]-totalEnd[120], N[1,1,1,120], N[8,1,1,120], N[9,1,1,120], (N[8,1,1,120] + N[9,1,1,120]), k[120], kV[120], intersect(120,v)))
    print(badboxes)
    #print(totalStart[badboxes] - totalEnd[badboxes])
    #bs <- which(kV>0)
    #print(mean(kCS[bs]))
    #print(paste(colSums(N[1:7,1,1,badboxes]),N[9,1,1,badboxes]))
  }
  
  S0 <- InitVals$S
  Vg0 <- InitVals$Vg
  incidentbybox <- (S0+Vg0)*(1-exp(-k*tstep)) # NOT QUITE CORRECT WHEN VACCINATING (Ignores Vg -> Vm, probably OK because slow compared to time step)!!!!
  #incidentbybox2 <- (N[2,1,1,] - E0) + (N[3,1,1,] - I0) + (N[4,1,1,] - H0) + (N[5,1,1,] - Ht0) + (N[6,1,1,] - D0) + (N[7,1,1,] - R0) # INCORRECT when culling!!
  results <- list(N=N,incidentbybox=incidentbybox)
  return(results)
}

################################################################# episim #############################################################
episim <- function(N,r,dl,dd,TimeParams,beta,betaRatio,brmod,Rr0,kD0,DetectSick,srmcd,srmc,VacParams, CullParams, icell, LDMoveParams, UsedParams, PixParams, ModelType, Example) {
  
  ################################ Set parameters ##############  
  nbox=dim(N)[4]
  Nages=dim(N)[3]
  Ndcomb=dim(N)[2]
  Ndstates=dim(N)[1]

  weeks_to_sim <- TimeParams$weeks_to_sim
  nyear <- TimeParams$nyear
  len.wk <- TimeParams$len.wk 
  tstep <- TimeParams$tstep
  ActualDetectTime <- NULL
  ReallySickorDead <- NULL
  
  Vrate <- VacParams$Vrate
  if (is.na(Vrate)) Vrate <- 0
  Vdelay <- VacParams$Vdelay
  Veff <- VacParams$Veff
  vacrad <- VacParams$VacRadius
  
  Crate <- CullParams$Crate
  if (is.na(Crate)) Crate <- 0
  CullFrom <- CullParams$CullFrom
  Culldelay <- CullParams$Culldelay
  Cullrad <- CullParams$Cullrad
  CdelayPInfect <- CullParams$CdelayPInfect
  
  LDType <- LDMoveParams$LDType
  LDProb <- LDMoveParams$LDProb
  cowTypeMov <- LDMoveParams$cowTypeMov
  lrmc <- LDMoveParams$lrmc
  lrmcd <- LDMoveParams$lrmcd
  
  if (LDType=="roads") {
    RL <- LDMoveParams$RL 
    index_KeptRoads <- LDMoveParams$index_KeptRoads 
    road_list <- LDMoveParams$road_list
    RoadSpread <- LDMoveParams$RoadSpread
    IntersectInfo <- LDMoveParams$IntersectInfo # Intersections(RL, index_KeptRoads)
    extras = list(LDType=LDType)
  }
  
  # Data will be reported on a weekly basis.
  wk.prev=array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,weeks_to_sim))
  wk.inc=array(data=0,dim=c(Nages,nbox,weeks_to_sim))
  StartInc <- colSums(N[2:6,1,,])  # N[states,diseases,ages,box]
  
  # Calculate pixel parameters
  ApproxPixelArea <- PixParams$pixel_width*PixParams$pixel_height
  pixel_size <- (PixParams$pixel_width + PixParams$pixel_height)/2
  
  usedvacc <- 0
  srmc_stime <- .Machine$double.xmax # start time for transmission control 
  lrmc_stime <- .Machine$double.xmax # start time for reduction of long range movement (i.e. reduce LDProb)
  Vdelay_corr <- .Machine$double.xmax
  VdelayH_corr <- .Machine$double.xmax
  Cdelay_corr <- .Machine$double.xmax 
  
  ################################ Spreading matrix set-up ##############
  middle_location <- GetMiddleCoords(r)
  MiddleCell <- cellFromXY(r,middle_location)
  ddWeight_sr <- IsoSpreadMatrix(dd, pixel_size, dl, MiddleCell)  # Isotropic nearby spread
  local_weight <- diag(1,nbox,nbox) # No isotropic spread only get animals in same cell sick; for the really sick who are unlikely to move w/o help
  
  LD_weights <- array(data=0,dim=(c(nbox,nbox))) # 
  LD_weights_con <- array(data=0,dim=(c(nbox,nbox))) 
  LD_weights_nocon <- array(data=0,dim=(c(nbox,nbox)))
  
  if (LDType=="P2P") {
    LDcells <- LDMoveParams$Lcellnums
    LD_weights_nocon[LDcells,] <- LDProb/100
    LD_weights_con[LDcells,] <- LDProb*lrmc/100
  }
  
  init_choosencells <- seq(1:nbox)
  if (LDType == "roads")  { 
      LD_weights <- road_transport(init_choosencells, RL, index_KeptRoads, road_list, RoadSpread, IntersectInfo, r)
      LD_weights_nocon <- LD_weights*(LDProb/100)*nbox/sum(LD_weights)    # nocon - no control in this case of directed "long distance" motion
      LD_weights_con <- lrmc*LD_weights*(LDProb/100)*nbox/sum(LD_weights)
  }
  
  # AbletoMoveKernels <- function(beta, srmc, LDProb, lrmc, ddWeight_sr, LD_weights_nocon) {
  #   ddWeight_nocon <- beta*((1-LDProb/100)*ddWeight_sr + LD_weights_nocon)        # No control of movement or transmission
  #   ddWeight_LDcon <- beta*(1-LDProb*lrmc/100)*ddWeight_sr + LD_weights_con       # Reduction of long range distance movement (but becomes short range spread)
  #   ddWeight_all_con <- srmc*ddWeight_LDcon                                       # Reduction of transmission (lower beta) and long range movement
  #   ddWeight_Tcon <- srmc*ddWeight_nocon                                          # Reduction of transmission (lower beta)
  #   AMK <- list(ddWeight_nocon=ddWeight_nocon, ddWeight_LDcon=ddWeight_LDcon, ddWeight_all_con=ddWeight_all_con, ddWeight_Tcon=ddWeight_Tcon)  
  #   return(AMK)
  # }
  
  # UnabletoMoveKernels <- function(beta, srmc, LDProb, lrmc, local_weight, LD_weights_nocon) {  # Animals or People too sick to spread disease themselves.
  #   ddWeight_local_nocon <- beta*((1-LDProb/100)*local_weight + LD_weights_nocon)
  #   ddWeight_local_LDcon <- beta*((1-LDProb*lrmc/100)*local_weight + LD_weights_nocon)
  #   ddWeight_local_all_con <- beta*srmc*ddWeight_local_LDcon
  #   ddWeight_local_Tcon <- beta*srmc*ddWeight_local_nocon
  #   UMK <- list(ddWeight_local_nocon=ddWeight_local_nocon, ddWeight_local_LDcon=ddWeight_local_LDcon, ddWeight_local_all_con=ddWeight_local_all_con, ddWeight_local_Tcon=ddWeight_local_Tcon)  
  #   return(UMK)
  # }
 
  # LDRestrictedKernels <- function(beta, srmc, ddWeight_sr, local_weight) {  # Long distance movement is not allowed for some groups
  #   weights <- beta*ddWeight_sr
  #   weights_con <- srmc*beta*ddWeight_sr
  #   weights_local <- beta*local_weight
  #   weights_local_con <- srmc*beta*local_weight
  #   RMK <- list(weights=weights, weights_con=weights_con, weights_local=weights_local, weights_local_con=weights_local_con)
  # }
  
  # Precompute Matrices related to disease spread to save time 
  AMK <- AbletoMoveKernels(beta, srmc, LDProb, lrmc, ddWeight_sr, LD_weights_nocon, LD_weights_con, ModelType, Example)
  UMK <- UnabletoMoveKernels(beta, srmc, LDProb, lrmc, local_weight, LD_weights_nocon, LD_weights_con, ModelType, Example)
  RMK <- LDRestrictedKernels(beta, srmc, ddWeight_sr, local_weight, ModelType, Example)
  
  ############## Rates for Diff Eqs  ##############  
  klist <- SetupRates(kD0,Rr0,1)
  Ks <- Define_Ks(klist, Vrate, Crate, 0, tstep)
  vaccinating <- FALSE
  culling <- FALSE
  
  CellsChoosen7daysAgo <- NULL
  CellsChoosen6daysAgo <- NULL
  CellsChoosen5daysAgo <- NULL
  CellsChoosen4daysAgo <- NULL
  CellsChoosen3daysAgo <- NULL
  CellsChoosen2daysAgo <- NULL
  CellsChoosenYesterday <- NULL
  CellsChoosenToday <- NULL
  
  ############## year loop  ##############  
  for (year in 1:nyear) {
    yeartime <- (year-1)*weeks_to_sim*len.wk
    withProgress(message = 'computing', value = 0, {
      for (week in 1:weeks_to_sim) { 
        incProgress(1/weeks_to_sim, detail = paste("week", week))
        #temporary arrays for saving data during a week
        record        <- array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,len.wk))
        recordH        <- array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,len.wk))
        incidence     <- array(data=0,dim=c(Nages,nbox,len.wk)) # incidence in a single day 
        incidenceH     <- array(data=0,dim=c(Nages,nbox,len.wk)) # incidence in a single day 
        step_incidence <- array(data=0,dim=c(Nages,nbox))  #incidence in a single time step
      
        print(paste("starting week", week))
        weektime <- (week-1)*len.wk

        ############## day loop  ##############
        for (day in 1:len.wk) {  #  tstep must be 1 or less!!!
          DetectTime <- day + weektime + yeartime
          if (sum(N[4:7,1,,])>DetectSick) {   
            if (srmc_stime == .Machine$double.xmax) {  
              DetectTime <- day + weektime + yeartime
              srmc_stime <- srmcd + DetectTime  #srmcd is an input: Delay of starting short range movement control after detection of disease
              Vdelay_corr <- Vdelay + DetectTime
              Cdelay_corr <- Culldelay + DetectTime
              lrmc_stime <- lrmcd + DetectTime
              ActualDetectTime <- DetectTime
              ReallySickorDead <- sum(N[4:6,1,,]) + sum(N[7,1,,])*klist$kIH/(klist$kIH+klist$kIR)
            }
          } 
        
          ############## Cull Vacc set-up  ##############
          current_time <- day+weektime+yeartime
          if ((Vdelay_corr < current_time) && (Vrate>0)) vaccinating <- TRUE

          if ((Cdelay_corr < current_time) && (Crate>0)) culling <- TRUE
        
          CellsToCull <- NULL
          if (culling) {
            CellsChoosen7daysAgo <- CellsChoosen6daysAgo
            CellsChoosen6daysAgo <- CellsChoosen5daysAgo
            CellsChoosen5daysAgo <- CellsChoosen4daysAgo
            CellsChoosen4daysAgo <- CellsChoosen3daysAgo
            CellsChoosen3daysAgo <- CellsChoosen2daysAgo
            CellsChoosen2daysAgo <- CellsChoosenYesterday
            CellsChoosenYesterday <- CellsChoosenToday
            
            sums <- colSums(N[3:5,1,1,])
            CellsChoosenNow <- which(sums > animalCut/26*ApproxPixelArea)  # Some locations with high populations and high beta still have enough infectious population to progress
            # Possibly use something like beta&S*(I+H+Ht) rather than just (I+H+Ht)
            CullCells <- unique(unlist(lapply(CellsChoosenNow, function(x) which(dd[x,] <= Cullrad))))   # Find all cells within vacrad of a cell with sick cows.
            #CullCellsswP <- CullCells[which(N[1,1,1,CullCells] > animalCut/26*ApproxPixelArea/10)]      # Cells with susceptible population (wP - with population)
            #CellsChoosenToday <- CullCellsswP
            CellsChoosenToday <- CullCells
            maxcell <- which(N[3,1,1,] == max(N[3,1,1,]) )
            #CellsToCull <- which(sums > animalCut/26*ApproxPixelArea)
            
            if (CdelayPInfect == 0) {
              CellsToCull <- CellsChoosenToday
              #CellsToCull <- c(CellsChoosenToday,CellsChoosenYesterday, CellsChoosen2daysAgo)
            }
            if (CdelayPInfect == 1) {
              CellsToCull <- CellsChoosenYesterday
            }
            if (CdelayPInfect == 2) {
              CellsToCull <- CellsChoosen2daysAgo
            }
            if (CdelayPInfect == 3) {
              CellsToCull <- CellsChoosen3daysAgo
            }
            if (CdelayPInfect == 4) {
              CellsToCull <- CellsChoosen4daysAgo
            }
            if (CdelayPInfect == 5) {
              CellsToCull <- CellsChoosen5daysAgo
            }
            if (CdelayPInfect == 6) {
              CellsToCull <- CellsChoosen6daysAgo
            }
            if (CdelayPInfect == 7) {
              CellsToCull <- CellsChoosen7daysAgo
            }
          }  
        
          CellsToVacc <- NULL
          if (vaccinating) {
            choosencells <- which(colSums(N[3:5,1,1,]) > animalCut/26*ApproxPixelArea)
            vaccells <- unique(unlist(lapply(choosencells, function(x) which(dd[x,] <= vacrad))))   # Find all cells within vacrad of a cell with sick cows.
            vaccellswP <- vaccells[which(N[1,1,1,vaccells] > animalCut/26*ApproxPixelArea/10)]      # Cells with susceptible population (wP - with population)
            CellsToVacc <- vaccellswP
          }
        
          DelayOfSecondMit <- 0
          MM <- SpreadingMatricesForDay(LDProb, current_time, lrmc_stime, srmc_stime, AMK, UMK, RMK, cowTypeMov, DelayOfSecondMit)
            
          ############## tstep loop  ##############
          for (n in 1:(1/tstep)) {  #  1/tstep must be an integer
          
            SusceptPVaccStart <- N[1,1,1,] + N[8,1,1,] 
            Nold <- N
            results <- NextN(N, ApproxPixelArea, MM, tstep, klist, Ks, CellsToVacc, CellsToCull, CullFrom) 
            N <- results$N
            incidence[1,,day] <- incidence[1,,day] + results$incidentbybox
            
          }  #End of tstep loop
          record[,,,,day] <- N  # record <- array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,len.wk))
          #print(paste("total from record", sum(record[,,,,day]), log10(sum(record[,,,,day]))))
        }  #end day loop 
        
        # bin (average) incidence and (sum) prevalence daily data into wk.prev, and wk.inc 
        # wk.prev=array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,len.year))
        # wk.inc= array(data=0,dim=c(Nages,nbox,weeks_to_sim))
        # incidence     <- array(data=0,dim=c(Nages,nbox,len.wk))

        # wk.prev[,,,,week]=rowMeans(record,dims=4)  #Mean of all days
        # wk.prev[6,,,,week] <- record[6,,,,7]  # Want dead to be total at end of week
        # wk.prev[7,,,,week] <- record[7,,,,7]  # Want recovered to be total at end of week
        # wk.prev[9,,,,week] <- record[9,,,,7]  # Want immune due vacination to be total at end of week
        
        wk.prev[,,,,week]=record[,,,,7]  # Value at end of week
        
        # print(paste("total from wk.prev", log10(sum(wk.prev[,,,,week]))))
        wk.inc[,,week]=rowSums(incidence,dims=2)                          # Summing over days in the week 
        if (week==1) {
          wk.inc[,,1]=rowSums(incidence,dims=2) + StartInc
        } else {
          wk.inc[,,week]=rowSums(incidence,dims=2)                          # Summing over days in the week 
        }  
      }  #end week loop
    })  # end withProgress
  }# end year 
  out=list(wk.prev=wk.prev,wk.inc=wk.inc, usedvacc=usedvacc, ActualDetectTime=ActualDetectTime, ReallySickorDead=ReallySickorDead, Cdelay_corr=Cdelay_corr, UsedParams=UsedParams)  
  return(out)
}

########################################################################### runSim ###################################################
runSim=function(N,r,dl,dd,icell,TimeParams,beta,betaRatio,Rr0,kD0,DetectSick,srmcDelay,srmcBetaChange,VacParams,CullParams, LDMoveParams, UsedParams, PixParams) {
  cat("In runSim", "(functions)", "\n")
  # Ndstates=dim(N)[1]
  # Ndcomb=dim(N)[2]
  # Nages=dim(N)[3]
  nbox=dim(N)[4]
  
  N[2,1,,icell] <- 5 #1  # index case(s) defined.
  Disease <- UsedParams$DiseaseModeled
  if (Disease == "Rinderpest lineage 1") {
    N[2,1,1,icell] <- 6
    N[3,1,1,icell] <- 2
    N[4,1,1,icell] <- 1
  }
  if (Disease == "Rinderpest lineage 2") {
    N[2,1,1,icell] <- 3
    N[3,1,1,icell] <- 2
    N[4,1,1,icell] <- 2
  }
  if (Disease == "Rinderpest Pak94") {
    N[2,1,1,icell] <- 3#4
    N[3,1,1,icell] <- 2#3
    N[4,1,1,icell] <- 1#2
  }
  
  startingIcellDisProg <- c(N[2,1,1,icell], N[3,1,1,icell], N[4,1,1,icell], N[5,1,1,icell], N[6,1,1,icell], N[7,1,1,icell])
  print(paste("Number in icell in E I H, Ht, D, R"))
  print(startingIcellDisProg)
  
  N[1,1,,icell] <- N[1,1,,icell] - sum(startingIcellDisProg)
  #print(N[1,1,,icell])
  
  sim<-episim(N,r,dl,dd,TimeParams,beta,betaRatio,brmod,Rr0,kD0,DetectSick,srmcDelay,srmcBetaChange,VacParams,CullParams, icell, 
              LDMoveParams, UsedParams, PixParams, "NotUsed", "NotUsed")  
}

diseasekD <- function(Disease) {
  kD2dead  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  kD2recov  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  if (Disease == "Rinderpest") {
    kD1dead   <- c(0, kEI_Rind, kIH_Rind, kHD_Rind, kItHt_Rind, kHtD_Rind, 0, 0, 0)  
    kD1recov  <- c(0, 0,        kIR_Rind, kHR_Rind, kItR_Rind,  kHtR_Rind, 1/vac_TimeToImmunity_Rind, 0, 0)  
  }
  if (Disease == "Rinderpest lineage 1") {
    kD1dead   <- c(0, 1/5.6, 0.33,  0.12, kItHt_Rind, kHtD_Rind, 0, 0, 0)  
    kD1recov  <- c(0,     0,   0, 0.2133, kItR_Rind,  kHtR_Rind, 1/vac_TimeToImmunity_Rind, 0, 0)  
  }
  if (Disease == "Rinderpest lineage 2") {
    kD1dead   <- c(0, 1/5.6, 0.14,   0.02, kItHt_Rind, kHtD_Rind, 0, 0, 0)  
    kD1recov  <- c(0,     0,    0, 0.3133, kItR_Rind,  kHtR_Rind, 1/vac_TimeToImmunity_Rind, 0, 0)  
  }
  if (Disease == "Rinderpest Pak94") {
    kD1dead   <- c(0, kEI_Pak94, kIH_Pak94, kHD_Pak94, kItHt_Rind, kHtD_Rind, 0, 0, 0)  
    kD1recov  <- c(0,         0, kIR_Pak94, kHR_Pak94, kItR_Rind,  kHtR_Rind, 1/vac_TimeToImmunity_Rind, 0, 0)  
  }
  return(cbind(kD1dead,kD1recov,kD2dead,kD2recov))
}

diseaseRr <- function(Disease) {
  # Rate Ratio           kHHt/(kHD+kHR)              
  Rr1treat <- c(0, 0, 0, Rr_HHt_Rind,       0, 0, 0, 0, 0)  # Changed Ht-R from 10 to 100 July 2016
  Rr2treat <- NULL
  return(cbind(Rr1treat,Rr2treat))
}


