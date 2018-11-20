################# Change the resolution of raster data  ############
Reduce_res <- function(RasterData,ResolutionReduction){
  if (is.null(RasterData))
    return(NULL)
  
  if (ResolutionReduction > 1) {#For integers get same result as aggregate command (which is simple aggregation)
    xPixNum <- round(ncol(RasterData)/ResolutionReduction)
    yPixNum <- round(nrow(RasterData)/ResolutionReduction)
    newRaster <- raster(ncol=xPixNum, nrow=yPixNum, xmn=extent(RasterData)@xmin, xmx=extent(RasterData)@xmax, ymn=extent(RasterData)@ymin, ymx=extent(RasterData)@ymax)
    newRaster  <- resample(RasterData, newRaster, method='bilinear')  #method='bilinear'
  }
  if (ResolutionReduction < 1) {
    IncRes <- round(1/ResolutionReduction)
    newRaster <- disaggregate(RasterData, fact=IncRes)
  }
  totalPop <- sum(getValues(RasterData))
  sumnewRaster <- sum(getValues(newRaster))
  newRaster <- newRaster*totalPop/sumnewRaster 
  return(newRaster)
}

LatLongLengths <- function(SizeOfModel) {
  if (SizeOfModel == "Small Rect: 1 x 1.5 degrees") {
    lat <- lat_lengths_SmallRect
    long <- long_lengths_SmallRect
  }
  if (SizeOfModel == "Med Rect: 1.3 x 2 degrees") {
    lat <- lat_lengths_MedRect
    long <- long_lengths_MedRect
  }
  if (SizeOfModel == "Large Rect: 3 x 4 degrees") {
    lat <- lat_lengths_LargeRect
    long <- long_lengths_LargeRect
  }
  if (SizeOfModel == "Large: 3 x 3 degrees") {
    lat <- region_lengths_Basic
    long <- region_lengths_Basic
  }
  if (SizeOfModel == "Medium: 2 x 2 degrees") {
    lat <- region_lengths_Medium
    long <- region_lengths_Medium
  }
  if (SizeOfModel == "Small: 0.8 x 0.8 degrees") {
    lat <- region_lengths_HighRes
    long <- region_lengths_HighRes
  }
  lengths <- list(lat=lat,long=long)
  return(lengths)
}

GetMiddleCoords <- function(r) {
  if (is.null(r))
    return(NULL)
  
  top <- extent(r)@ymax
  bottom <- extent(r)@ymin
  left <- extent(r)@xmin
  right <- extent(r)@xmax
  xval <- (left + right)/2.0
  yval <- (top + bottom)/2.0
  middlecoords <- c(xval,yval)
  return(middlecoords)
}

IsoSpreadMatrix <- function(dd,pixel_size,dl,icell) {
  nbox <- dim(dd)[1]
  dd_weight <- array(data=0,dim=(c(nbox,nbox)))
  dd_weight <- (dl/pixel_size)*exp(-(dd/dl))*2*sinh(pixel_size/(2*dl)) # JRM weighting over length of box
  dd_weight[which(dd==0, arr.ind = TRUE)] <- 2*(dl/pixel_size)*(1-exp(-1*pixel_size/(2*dl)))
  SRnormVal <- sum(dd_weight[,icell])
  dd_weight <- dd_weight/SRnormVal
  SRnormVal <- sum(dd_weight[,icell])
  dd_weight_sr0 <- dd_weight/SRnormVal  #Normalizing each column - for use when there is no Long Distance movement
  #The idea here is to have the same total transmissibility regardless of dl.
  return(dd_weight_sr0)
} 

SpreadXY <- function(CentLat, CentLong, r, pixel_size, dl) {
  lat <- 0.7
  long <- 0.7
  longmin <- CentLong-long/2
  longmax <- CentLong+long/2
  latmin <- CentLat-lat/2
  latmax <- CentLat+lat/2
  wb   <- extent(longmin,longmax,latmin,latmax)
  rsmall <- crop(r, wb)
  
  dd <- distm(coordinates(rsmall),fun=distMeeus)/1000
  IndexCell <- cellFromXY(rsmall,c(CentLong,CentLat))
  
  middle_location <- GetMiddleCoords(rsmall)
  MiddleCell <- cellFromXY(rsmall,middle_location)
  
  dd.weight_sr0 <- IsoSpreadMatrix(dd, pixel_size, dl, MiddleCell)
  icellProbs <- dd.weight_sr0[,MiddleCell]
  DistsToBox <- dd[,MiddleCell]
  index_order <- order(dd[,MiddleCell])  
  SpreadXYvals <- list(x=DistsToBox[index_order], y=log(icellProbs[index_order]))
  return(SpreadXYvals)
}

SetupRates <- function(kD0,Rr0,OneOrTwo)  {
  kD0dead <- kD0$kD0dead
  kD0recov <- kD0$kD0recov
  
  kEI   <- kD0dead[2,OneOrTwo,1] 
  
  kIH   <- kD0dead[3,OneOrTwo,1] #kD1[3]
  kIR   <- kD0recov[3,OneOrTwo,1]
  
  kHD   <- kD0dead[4,OneOrTwo,1] #kD1[4]
  kHR   <- kD0recov[4,OneOrTwo,1] #kHD*0
  
  kItD  <- kD0dead[5,OneOrTwo,1] #kD1[6]          # Not used yet
  kItR  <- kD0recov[5,OneOrTwo,1]                 # Not used yet
  
  kHtD  <- kD0dead[6,OneOrTwo,1] #kD1[6]
  kHtR  <- kD0recov[6,OneOrTwo,1]
  
  kIIt  <- (kIH+kIR)*Rr0[3,OneOrTwo,1]           # Not used yet
  kHHt  <- (kHD+kHR)*Rr0[4,OneOrTwo,1]          
  
  ke    <- kD0recov[7,OneOrTwo,1]
  
  klist <- list(kEI=kEI,kIH=kIH,kIR=kIR,kHD=kHD,kHR=kHR,kHtD=kHtD,kHtR=kHtR,kHHt=kHHt,ke=ke)  
  return(klist)
}

Define_Ks <- function(klist, Vrate, Crate, betaRatio, tstep) {
  kEI   <- klist$kEI    #   kEI   <- kD1[2]
  kIH   <- klist$kIH    #   kIH   <- kD1[3]
  kIR   <- klist$kIR    #   kIR   <- (kD1[3]+Rr1[3]*kD1[3])*Rr1[2]
  kHHt  <- klist$kHHt   #   kHHt  <- 
  kHD   <- klist$kHD    #   kHD   <- kD1[4]
  kHR   <- klist$kHR    #   kHD   <- kD1[4]
  kHtR  <- klist$kHtR   #   kHtR  <- 
  kHtD  <- klist$kHtD   #   kHtD  <- kD1[6]
  ke    <- klist$ke
  
  K1    <- kIH + kIR 
  K2    <- kHHt + kHD + kHR
  K4    <- kHtD + kHtR
  
  # Avoid redundant eigenvalues. 
  # Since values are not known to an accuracy of even a 10th of a day, seems OK.
  if (K1 == K2) {  # increase K1
    kIH   <- kIH + 0.0001  
    K1    <- kIH + kIR 
  }
  if (K4 == K2)  {  # increase K4
    kHtR  <- kHtR + 0.0001
    K4    <- kHtD + kHtR 
  }
  if (K1==K4) {  # increase K4 and recheck K4==K2
    kHtD  <- kHtD + 0.0001
    K4    <- kHtD + kHtR
    if (K4 == K2)  {  # increase K4
      kHtR  <- kHtR + 0.0001
      K4    <- kHtD + kHtR 
    }
  }
  
  if (kEI == K1) {  # increase KEI
    kEI <- kEI + 0.0001
  }
  if (kEI == K2) { # increase KEI and recheck K1==kEI
    kEI <- kEI + 0.0001
    if (kEI == K1) {  # increase KEI
      kEI <- kEI + 0.0001
    } 
  }  
  if (kEI == K4) {
    kEI <- kEI + 0.0001
    if (kEI == K1) {  # Check kEI==K1 (know K1 and K2 are different)
      kEI <- kEI + 0.0001
    }
    if (kEI == K2) {  # Check kEI==K2 (know K1 and K2 are different)
      kEI <- kEI + 0.0001
    }
  }
  
  VaccPerStep <- Vrate*tstep
  
  Ks <- list(K1=K1, K2=K2, K4=K4, kIH=kIH, kHtR=kHtR, kIR=kIR, kHHt=kHHt, kHtD=kHtD, VaccPerStep=VaccPerStep, Crate=Crate, betaRatio=betaRatio)
  
  return(Ks)
}  

Reduce_res <- function(RasterData,ResolutionReduction){
  if (is.null(RasterData))
    return(NULL)
  
  if (ResolutionReduction > 1) {#For integers get same result as aggregate command (which is simple aggregation)
    xPixNum <- round(ncol(RasterData)/ResolutionReduction)
    yPixNum <- round(nrow(RasterData)/ResolutionReduction)
    newRaster <- raster(ncol=xPixNum, nrow=yPixNum, xmn=extent(RasterData)@xmin, xmx=extent(RasterData)@xmax, ymn=extent(RasterData)@ymin, ymx=extent(RasterData)@ymax)
    newRaster  <- resample(RasterData, newRaster, method='bilinear')  #method='bilinear'
  }
  if (ResolutionReduction < 1) {
    IncRes <- round(1/ResolutionReduction)
    newRaster <- disaggregate(RasterData, fact=IncRes)
  }
  totalPop <- sum(getValues(RasterData))
  sumnewRaster <- sum(getValues(newRaster))
  newRaster <- newRaster*totalPop/sumnewRaster 
  return(newRaster)
}

InitialValues <- function(N) {
  S0  <- N[1,1,1,]  #disease 1 and age group 1,  Watch out for nages!!!
  E0  <- N[2,1,1,] 
  I0  <- N[3,1,1,]
  H0  <- N[4,1,1,]
  Ht0 <- N[5,1,1,]
  D0  <- N[6,1,1,]
  R0  <- N[7,1,1,]
  Vg0 <- N[8,1,1,]
  Vim0<- N[9,1,1,]
  InitVals <- list(S=S0,E=E0,I=I0,H=H0,Ht=Ht0,D=D0,R=R0,Vg=Vg0,Vim=Vim0)
  return(InitVals)
}

kEIValues <- function(kEIvalue, InitVals, CutOff, Nboxes) {
  E0 <- InitVals$E
  D0 <- InitVals$D
  R0 <- InitVals$R
  ieEnough <- which(E0 > CutOff) # index for enough E 
  ieStarted <-  which((D0+R0) > 0)    # if already started propagating then finish 
  ieE <- union(ieEnough, ieStarted)
  kEI <- rep(0,Nboxes)
  kEI[ieE] <- kEIvalue
  return(kEI)
}

AvoidRedundantEigenValues <- function(CullRates, K1, K2, K4, kEIvalue, ke, kV) {
  kCS <- CullRates$kCS
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  CrateDiffIH <- kCI-kCH
  CrateDiffEI <- kCE-kCI
  CrateDiffEH <- kCE-kCH
  kVkCSDiff <- kV-kCS
  
  ### There is no check for k = kEI, K1, K2, K4.  
  if (any(kVkCSDiff==ke)) {           # Assuming here that ke is never 0.  There is no vaccine that causes immediate immunity. !!!
    index <- which(kVkCSDiff == ke)
    kCS[index] <- kCS[index] + 0.001
  }
  
  if (any(CrateDiffIH == K2-K1))  {   # Would have redundant eigenvalue problem.  The case K1=K2 was taken care of in Define_Ks.
    index <- which(CrateDiffIH == K2-K1)
    kCH[index] <- kCH[index] + 0.001  
  }
  if (any(kCH == K4-K2))  {           # Would have redundant eigenvalue problem.  The case K2=K4 was taken care of in Define_Ks.
    index <- which(kCH == K4-K2)
    kCH[index] <- kCH[index] + 0.001
  }
  if (any(kCI == K4-K1))  {           # Would have redundant eigenvalue problem.  The case K1=K4 was taken care of in Define_Ks.
    index <- which(kCH == K4-K1)
    kCH[index] <- kCI[index] + 0.0005
  }
  
  if (any(CrateDiffEI == K1-kEIvalue)) {  # Would have redundant eigenvalue problem.  The case K1=kEIvalue was taken care of in Define_Ks.
    index <- which(CrateDiffEI == K1-kEIvalue)
    kCE[index] <- kCE[index] + 0.0002
  }
  if (any(CrateDiffEH == K2-kEIvalue)) {  # Would have redundant eigenvalue problem.  The case K2=kEIvalue was taken care of in Define_Ks.
    index <- which(CrateDiffEH == K2-kEIvalue)
    kCE[index] <- kCE[index] + 0.0002
  }
  if (any(kCE == K4-kEIvalue)) {          # Would have redundant eigenvalue problem.  The case K1=kEIvalue was taken care of in Define_Ks.
    index <- which(kCE == K4-kEIvalue)
    kCE[index] <- kCE[index] + 0.0002
  }
  CullRates <- list(kCS=kCS,kCE=kCE,kCI=kCI,kCH=kCH)
  return(CullRates)
}

kIndepParams <- function(kEI, kIH, kHHt, CullRates, K1, K2, K4, tstep) {
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  expkEI  <- exp(-(kEI+kCE)*tstep)
  expK1   <- exp(-(K1+kCI)*tstep) 
  expK2   <- exp(-(K2+kCH)*tstep)
  expK4   <- exp(-K4*tstep)
  
  if (all(kCI==kCE)) {
    fkEI1 <- kEI/(K1-kEI)
  } else {
    fkEI1 <- kEI/(K1+kCI-(kEI+kCE))     # Can't be infinite because K1 > 0, and K1+kCI=kEI+kCE was taken care of in AvoidRedundantEigenvalues
  }
  if (all(kCH==kCE)) {
    fkEI2 <- kIH/(K2-kEI)
  } else {
    fkEI2 <- kIH/(K2+kCH-(kEI+kCE))     # Can't be infinite because K2 > 0, and K2+kCH=kEI+kCE was taken care of in AvoidRedundantEigenvalues    
  }  
  fkEI4 <- kHHt/(K4-(kEI+kCE))        # K4=kEIValue+kCE was taken care of in AvoidRedundantEigenvalues. Used in Eigenvector 4 with check for kEI>0.           
  
  fK1a <- kIH/(K2+kCH-(K1+kCI))       # Can't be infinite because K1 > 0, and K2+kCH=(K1+kCI) was taken care of in AvoidRedundantEigenvalues
  fK1e <- kHHt/(K4-(K1+kCI))          # Can't be infinite because K1 > 0, and K4-K1=kCI was taken care of in AvoidRedundantEigenvalues
  fK4d <- kHHt/(K4-(K2+kCH))          # Can't be infinite because K2 > 0, and K4-K2=kCH was taken care of in AvoidRedundantEigenvalues
  
  kIP <- list(expkEI=expkEI,expK1=expK1,expK2=expK2,expK4=expK4,fkEI1=fkEI1,fkEI2=fkEI2,fkEI4=fkEI4,fK1a=fK1a,fK1e=fK1e,fK4d=fK4d) 
  return(kIP)
} 

kDepParams <- function(kEI, kIH, kHHt, k, ke, kV, CullRates, K1, K2, K4) {
  kCS <- CullRates$kCS
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  ks  <- k+ke  
  Ka  <- k+kV+kCS
  fK1 <- k/((kEI+kCE)-Ka)
  fK2 <- kEI/(K1+kCI-Ka)
  fK3 <- kIH/(K2+kCH-Ka)
  fK5 <- kHHt/(K4-Ka)
  fkV <- kV/(ke-(kV+kCS)) 
  kDP <- list(ks=ks,Ka=Ka,fK1=fK1,fK2=fK2,fK3=fK3,fK5=fK5,fkV=fkV)
  return(kDP)
}  

SpreadingMatricesForDay <- function(LDProb, current_time, lrmc_stime, srmc_stime, AMK, UMK, RMK, cowTypeMov, DelayOfSecondMit) {
  # Choose spreading matrix based on what movement controls are in place at this time
  if (LDProb > 0) {
    if (lrmc_stime <= current_time){
      if (srmc_stime <= current_time) {  # motion and transmission controls in place
        weights_w_LDmov <- AMK$ddWeight_all_con
        weights_wo_LDmov <- RMK$weights_con
        weights_no_iso_w_LDmov <- UMK$ddWeight_local_all_con
        weights_no_iso_wo_LDmov <- RMK$weights_local_con
      } else {
        weights_w_LDmov <- AMK$ddWeight_LDcon
        weights_wo_LDmov <- RMK$weights
        weights_no_iso_w_LDmov <- UMK$ddWeight_local_LDcon
        weights_no_iso_wo_LDmov <- RMK$weights_local
      } 
    }
  }
  if ((LDProb == 0) || (lrmc_stime > current_time)) {
    if (srmc_stime > current_time) {
      weights_w_LDmov <- AMK$ddWeight_nocon
      weights_wo_LDmov <- RMK$weights
      weights_no_iso_w_LDmov <- UMK$ddWeight_local_nocon  # dd.weight_local_nocon <- beta*((1-LDProb/100)*local_weight + LD_weights_nocon),    LD_weights_nocon[LDcell,] <- LDProb/100
      weights_no_iso_wo_LDmov <- RMK$weights_local        # weights_local <- beta*local_weight
    }
    # This is for some of the case study simulations. Since Tcon1 and Tcon2 are the same for everything else it has no effect.
    
    if ((srmc_stime <= current_time) && (srmc_stime+DelayOfSecondMit > current_time)) {  # Only happens for special use cases
      weights_w_LDmov <- AMK$ddWeight_Tcon1
      weights_wo_LDmov <- RMK$weights_con
      weights_no_iso_w_LDmov <- UMK$ddWeight_local_Tcon1
      weights_no_iso_wo_LDmov <- RMK$weights_local_con1
    }
    if (srmc_stime+DelayOfSecondMit <= current_time) {
      weights_w_LDmov <- AMK$ddWeight_Tcon
      weights_wo_LDmov <- RMK$weights_con
      weights_no_iso_w_LDmov <- UMK$ddWeight_local_Tcon
      weights_no_iso_wo_LDmov <- RMK$weights_local_con
    }
  }
  
  # Choose spreading matrix for exposed, infectious, and really sick animals, based on movement restrictions.
  if (cowTypeMov == "all") {
    dd.weight_totE <- weights_w_LDmov
    dd.weight_totI <- weights_w_LDmov
    dd.weight_totHplus <- weights_no_iso_w_LDmov # kept this option - simulates sick flying or being driven
  }
  if (cowTypeMov == "NoH"){
    dd.weight_totE <- weights_w_LDmov
    dd.weight_totI <- weights_w_LDmov
    dd.weight_totHplus <- weights_no_iso_wo_LDmov  
  }
  if (cowTypeMov == "NoHI"){
    dd.weight_totE <- weights_w_LDmov
    dd.weight_totI <- weights_wo_LDmov
    dd.weight_totHplus <- weights_no_iso_wo_LDmov
  }
  # dd.weight_totE is not being used June 2017
  MM <- list(dd.weight_totE=dd.weight_totE, dd.weight_totI=dd.weight_totI, dd.weight_totHplus=dd.weight_totHplus)
  return(MM)
}

EigenVector1NoVacc <- function(Zs, Os, k, kEI, kIH, kIR, kHD, kHR, kHHt, kHtR, kHtD, K1, K2, K4, CullRates) {
  kCS <- CullRates$kCS
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  #print(paste("In EigenVector1NoVacc ", kCS[89], kCE[89], kCI[89], kCH[89]))
  
  if (all(kCE==kCS)) {
    fk1 <- k/(kEI-k)
  } else {
    fk1 <- k/(kEI+kCE-(k+kCS))
  } 
  if (all(kCI == kCS)) {
    fk2 <- kEI/(K1-k)
  } else {
    fk2 <- kEI/(K1+kCI-(k+kCS))
  }
  if (all(kCH == kCS)){
    fk3 <- kIH/(K2-k)
  } else{
    fk3 <- kIH/(K2+kCH-(k+kCS))
  }
  fk5 <- kHHt/(K4-(k+kCS))
  
  #print(paste("In EigenVector1NoVacc ", fk1[89], fk2[89], fk3[89], fk5[89]))
  
  S   <- Zs           # alpha
  E   <- Zs 
  I   <- Zs 
  H   <- Zs 
  Ht  <- Zs 
  D   <- Zs
  R   <- Zs
  
  iE1 <- which(k!=0)
  
  S[iE1]   <- Os[iE1]   # alpha
  E[iE1]   <- fk1[iE1]  
  I[iE1]   <- fk1[iE1]*fk2[iE1]
  H[iE1]   <- fk1[iE1]*fk2[iE1]*fk3[iE1]  
  Ht[iE1]  <- fk1[iE1]*fk2[iE1]*fk3[iE1]*fk5[iE1]
  D[iE1]   <- -1*(kCS[iE1]*S[iE1] + kCE[iE1]*E[iE1] + kCI[iE1]*I[iE1] + (kHD+kCH[iE1])*H[iE1] + kHtD*Ht[iE1])/(k[iE1]+kCS[iE1])   
  R[iE1]   <- -1*(kIR*I[iE1] + kHR*H[iE1] + kHtR*Ht[iE1])/(k[iE1]+kCS[iE1])
  #print(paste("In EigenVector1NoVacc ", S[89], E[89], I[89], H[89], Ht[89], D[89], R[89]))  
  E1N <- list(S=S,E=E,I=I,H=H,Ht=Ht,D=D,R=R)
  return(E1N)
}  

EigenVector1wVacc <- function(Zs, Os, k, fkV, ke, Ka, fK1, fK2, fK3, fK5, CullRates, kHD, kHtD, kIR, kHR, kHtR){
  kCS <- CullRates$kCS
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  S   <- Os         
  E   <- Zs
  I   <- Zs
  H   <- Zs
  Ht  <- Zs
  D   <- Zs
  R   <- Zs
  Vg  <- fkV
  Vim <- -ke/Ka*fkV
  
  #kgt0 <- which(k!=0)  # if k=0  E, I, H, Ht , R and R are zero  
  
  kgt0 <- which(Ka!=0)  # if k=0  E, I, H, Ht , R and R are zero     Ka = k+kV+kCS
  
  E[kgt0]   <- fK1[kgt0]*(1+fkV[kgt0])  # function of k, length is Nboxes
  I[kgt0]   <- fK1[kgt0]*fK2[kgt0]*(1+fkV[kgt0])
  H[kgt0]   <- fK1[kgt0]*fK2[kgt0]*fK3[kgt0]*(1+fkV[kgt0])  
  Ht[kgt0]  <- fK1[kgt0]*fK2[kgt0]*fK3[kgt0]*fK5[kgt0]*(1+fkV[kgt0])
  D[kgt0]   <- -(kCS[kgt0]*S[kgt0] + kCE[kgt0]*E[kgt0] + kCI[kgt0]*I[kgt0] + (kHD+kCH[kgt0])*H[kgt0] + kHtD*Ht[kgt0])/(Ka[kgt0])  
  R[kgt0]   <- -(kIR*I[kgt0] + kHR*H[kgt0] + kHtR*Ht[kgt0])/(Ka[kgt0])
  E1 <- list(S=S,E=E,I=I,H=H,Ht=Ht,D=D,R=R,Vg=Vg,Vim=Vim)
  return(E1)
}

EigenVector2wVacc <- function(Zs, Os, k, kEI, kIH, kIR, kHHt, kHD, kHR, kHtD, kHtR, CullRates, ke, ks, K1, K2, K4) {
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  fks1 <- k/(kEI+kCE-ks)
  fks2 <- kEI/(K1+kCI-ks) 
  fks3 <- kIH/(K2+kCH-ks)
  fks5 <- kHHt/(K4-ks)
  
  S   <- Zs           # chi
  Vg  <- Os
  Vim <- -ke/ks
  # if k=0 all terms E thru R are 0
  #kgt0 <- which(k!=0)  # if k=0  E, I, H, Ht , R and R are zero  
  
  kgt0 <- which(ks!=0)  # if k=0  E, I, H, Ht , R and R are zero  ks = k + ke
  
  E   <- Zs
  I   <- Zs
  H   <- Zs
  Ht  <- Zs
  D   <- Zs
  R   <- Zs
  E[kgt0]   <- fks1[kgt0]
  I[kgt0]   <- fks1[kgt0]*fks2[kgt0]
  H[kgt0]   <- fks1[kgt0]*fks2[kgt0]*fks3[kgt0]
  Ht[kgt0]  <- fks1[kgt0]*fks2[kgt0]*fks3[kgt0]*fks5[kgt0]
  D[kgt0]   <- -(kCE[kgt0]*E[kgt0] + kCI[kgt0]*I[kgt0] + (kHD+kCH[kgt0])*H[kgt0] + kHtD*Ht[kgt0])/ks[kgt0]
  R[kgt0]   <- -(kIR*I[kgt0] + kHR*H[kgt0] + kHtR*Ht[kgt0])/ks[kgt0]
  E2 <- list(S=S,E=E,I=I,H=H,Ht=Ht,D=D,R=R,Vg=Vg,Vim=Vim)
  return(E2)
}

EigenVector4 <- function(Zs, kEI, fkEI1, fkEI2, fkEI4, CullRates, kHD, kHtD, kIR, kHR, kHtR)  {
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  E   <- Zs
  I   <- Zs
  H   <- Zs
  Ht  <- Zs
  D   <- Zs
  R   <- Zs
  
  kEIgt0 <- which( (kEI+kCE) > 0)
  E[kEIgt0]   <- 1  # ones           # B
  I[kEIgt0]   <- fkEI1[kEIgt0]*E[kEIgt0]           # C
  H[kEIgt0]   <- I[kEIgt0]*fkEI2[kEIgt0]           # D
  Ht[kEIgt0]  <- H[kEIgt0]*fkEI4[kEIgt0]           # E
  D[kEIgt0]   <- -1*(kCE[kEIgt0]*E[kEIgt0] + kCI[kEIgt0]*I[kEIgt0] + (kHD+kCH[kEIgt0])*H[kEIgt0] + kHtD*Ht[kEIgt0])/(kEI[kEIgt0] + kCE[kEIgt0])
  R[kEIgt0]   <- -1*(kIR*I[kEIgt0] + kHR*H[kEIgt0] + kHtR*Ht[kEIgt0])/(kEI[kEIgt0] + kCE[kEIgt0]) 
  #print(paste("In EigenVector4 ", E[89], I[89], H[89], Ht[89], D[89], R[89]))  
  E4 <- list(E=E,I=I,H=H,Ht=Ht,D=D,R=R)
  return(E4)
}

EigenVector5 <- function(Os, fK1a, fK1e, K1, CullRates, kHD, kHtD, kIR, kHR, kHtR) {
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
  
  I   <- Os  # ones
  H   <- fK1a  # When only vaccinating H, this has dimensions of nbox
  Ht  <- fK1e*fK1a
  D   <- -1*(kCI*I + (kHD+kCH)*H + kHtD*Ht)/(K1+kCI)
  R   <- -1*(kIR*I + kHR*H + kHtR*Ht)/(K1+kCI)
  E5 <- list(I=I,H=H,Ht=Ht,D=D,R=R)
  return(E5)
}

EigenVector6 <- function(Os,fK4d, K2, CullRates, kHD, kHtD, kHR, kHtR) {
  kCH <- CullRates$kCH
  H   <- Os  
  Ht  <- fK4d
  D   <- -1*((kHD+kCH)*H + kHtD*Ht)/(K2+kCH)
  R   <- -1*(kHR*H + kHtR*Ht)/(K2+kCH)
  E6 <- list(H=H,Ht=Ht,D=D,R=R)
  return(E6)
} 

EigenVector7 <- function(Os, K4, kHtD, kHtR, Nboxes) {  
  Ht   <- Os  #zeta, K4
  D   <- rep(-kHtD/K4, Nboxes)
  R   <- rep(-kHtR/K4, Nboxes)
  E7 <- list(Ht=Ht,D=D,R=R)
  return(E7)
}

N_NewNoVacc <- function(nv, N, Ndstates, Nboxes, E1N, E4, E5, E6, E7, k, CullRates, tstep, expkEI, expK1, expK2, expK4, InitVals, ke) {
  IV <- InitVals
  kCS <- CullRates$kCS
  E_coefs <- array(data = 0, c(Ndstates-2, Nboxes))
  E_coefs[1,nv] <- IV$S[nv]  # alpha
  E_coefs[2,nv] <- IV$E[nv] - E_coefs[1,nv]*E1N$E[nv] # beta
  E_coefs[3,nv] <- IV$I[nv] - E_coefs[1,nv]*E1N$I[nv] - E_coefs[2,nv]*E4$I[nv] # gamma
  E_coefs[4,nv] <- IV$H[nv] - E_coefs[1,nv]*E1N$H[nv] - E_coefs[2,nv]*E4$H[nv] - E_coefs[3,nv]*E5$H[nv] # delta
  E_coefs[5,nv] <- IV$Ht[nv] - E_coefs[1,nv]*E1N$Ht[nv] - E_coefs[2,nv]*E4$Ht[nv] - E_coefs[3,nv]*E5$Ht[nv] - E_coefs[4,nv]*E6$Ht[nv] # zeta
  
  expk   <- exp(-(k+kCS)*tstep)
  
  N[1,1,1,nv]  <- E_coefs[1,nv]*expk[nv] 
  
  N[2,1,1,nv] <- E_coefs[1,nv]*E1N$E[nv]*expk[nv] + E_coefs[2,nv]*expkEI[nv]  # Exposed
  
  N[3,1,1,nv] <- E_coefs[1,nv]*E1N$I[nv]*expk[nv] + E_coefs[2,nv]*E4$I[nv]*expkEI[nv] + E_coefs[3,nv]*expK1[nv] # Infectious
  
  N[4,1,1,nv] <- E_coefs[1,nv]*E1N$H[nv]*expk[nv] + E_coefs[2,nv]*E4$H[nv]*expkEI[nv] + E_coefs[3,nv]*E5$H[nv]*expK1[nv] + E_coefs[4,nv]*expK2[nv] # H
  
  N[5,1,1,nv]  <- E_coefs[1,nv]*E1N$Ht[nv]*expk[nv] + E_coefs[2,nv]*E4$Ht[nv]*expkEI[nv] + E_coefs[3,nv]*E5$Ht[nv]*expK1[nv]
  N[5,1,1,nv]  <- N[5,1,1,nv]+ E_coefs[4,nv]*E6$Ht[nv]*expK2[nv] + E_coefs[5,nv]*expK4
  
  N[6,1,1,nv]  <- E_coefs[1,nv]*E1N$D[nv]*(expk[nv]-1) + E_coefs[2,nv]*E4$D[nv]*(expkEI[nv]-1) + E_coefs[3,nv]*E5$D[nv]*(expK1[nv]-1) 
  N[6,1,1,nv] <- N[6,1,1,nv] + E_coefs[4,nv]*E6$D[nv]*(expK2[nv]-1) + E_coefs[5,nv]*E7$D[nv]*(expK4-1) + IV$D[nv]
  
  if (any(!is.finite(N[6,1,1,nv]))) print("Problem with D")
  
  N[7,1,1,nv]  <- E_coefs[1,nv]*E1N$R[nv]*(expk[nv]-1) + E_coefs[2,nv]*E4$R[nv]*(expkEI[nv]-1) + E_coefs[3,nv]*E5$R[nv]*(expK1[nv]-1) 
  N[7,1,1,nv]  <- N[7,1,1,nv] + E_coefs[4,nv]*E6$R[nv]*(expK2[nv]-1) + E_coefs[5,nv]*E7$R[nv]*(expK4-1) + IV$R[nv]
  
  if (any(!is.finite(N[7,1,1,nv]))) print("Problem with R")
  
  Ng8pos_nv <- intersect(nv, which(N[8,1,1,]>0))
  N8temp <- N[8,1,1,]
  N[8,1,1,Ng8pos_nv]  <- N8temp[Ng8pos_nv]*exp(-ke*tstep)
  N[9,1,1,Ng8pos_nv]  <- N[9,1,1,Ng8pos_nv] + N8temp[Ng8pos_nv]*(1-exp(-ke*tstep))
  if (any(!is.finite(N[9,1,1,Ng8pos_nv]))) print("N[9,1,1,Ng8pos_nv] not finite")
  # For the case Ng[8,nv] ==0, Ng[9,nv] also does not change.
  
  
  #print(paste(E_coefs[1,89]*expk[89], E_coefs[1,89]*E1N$D[89]*(expk[89]-1), E_coefs[2,89]*E4$D[89]*(expkEI[89]-1), E_coefs[3,89]*E5$D[89]*(expK1[89]-1), E_coefs[4,89]*E6$D[89]*(expK2[89]-1), E_coefs[5,89]*E7$D[89]*(expK4-1)))
  #print(paste(E_coefs[1,89]*E1N$D[89]*(expk[89]-1), E_coefs[1,89], E1N$D[89], expk[89]-1))
  
  return(N)
}

N_New_wVacc <- function(v, Ka, ks, N, Ndstates, Nboxes, E1, E2, E3, E4, E5, E6, E7, InitVals, expkEI, expK1, expK2, expK4, tstep) {
  
  E_coefs <- array(data = 0, c(Ndstates, Nboxes)) 
  IV <- InitVals
  
  E_coefs[1,v] <- IV$S[v]  # alpha
  E_coefs[2,v] <- IV$Vg[v]  - E_coefs[1,v]*E1$Vg[v]   # chi
  E_coefs[3,v] <- IV$Vim[v] - E_coefs[1,v]*E1$Vim[v] - E_coefs[2,v]*E2$Vim[v]  # psi
  E_coefs[4,v] <- IV$E[v]   - E_coefs[1,v]*E1$E[v]   - E_coefs[2,v]*E2$E[v] # beta
  E_coefs[5,v] <- IV$I[v]   - E_coefs[1,v]*E1$I[v]   - E_coefs[2,v]*E2$I[v]  - E_coefs[4,v]*E4$I[v] # gamma
  E_coefs[6,v] <- IV$H[v]   - E_coefs[1,v]*E1$H[v]   - E_coefs[2,v]*E2$H[v]  - E_coefs[4,v]*E4$H[v]  - E_coefs[5,v]*E5$H[v] # delta
  E_coefs[7,v] <- IV$Ht[v]  - E_coefs[1,v]*E1$Ht[v]  - E_coefs[2,v]*E2$Ht[v] - E_coefs[4,v]*E4$Ht[v] - E_coefs[5,v]*E5$Ht[v] - E_coefs[6,v]*E6$Ht[v] # zeta
  E_coefs[8,v] <- IV$D[v]   - E_coefs[1,v]*E1$D[v]   - E_coefs[2,v]*E2$D[v]  - E_coefs[4,v]*E4$D[v]  - E_coefs[5,v]*E5$D[v]  - E_coefs[6,v]*E6$D[v] - E_coefs[7,v]*E7$D[v] # eta
  E_coefs[9,v] <- IV$R[v]   - E_coefs[1,v]*E1$R[v]   - E_coefs[2,v]*E2$R[v]  - E_coefs[4,v]*E4$R[v]  - E_coefs[5,v]*E5$R[v]  - E_coefs[6,v]*E6$R[v] - E_coefs[7,v]*E7$R[v] # theta
  
  #expressions depending on k 
  expKa   <- exp(-Ka*tstep)
  expks   <- exp(-ks*tstep)
  
  N[1,1,1,v]  <- E_coefs[1,v]*expKa[v] 
  
  N[8,1,1,v]  <- E_coefs[1,v]*E1$Vg[v]*expKa[v] + E_coefs[2,v]*expks[v]
  
  N[9,1,1,v]  <- E_coefs[1,v]*E1$Vim[v]*expKa[v] + E_coefs[2,v]*E2$Vim[v]*expks[v] + E_coefs[3,v]
  if (any(!is.finite(E1$Vim[v]))) print("VimE1 not finite")
  if (any(!is.finite(E2$Vim[v]))) print("VimE2 not finite")
  if (any(!is.finite(E_coefs[3,v]))) print("E_coefs[3,v] not finite")
  
  N[2,1,1,v] <- E_coefs[1,v]*E1$E[v]*expKa[v] + E_coefs[2,v]*E2$E[v]*expks[v] + E_coefs[4,v]*expkEI[v]  # Exposed
  
  N[3,1,1,v] <- E_coefs[1,v]*E1$I[v]*expKa[v] + E_coefs[2,v]*E2$I[v]*expks[v] + E_coefs[4,v]*E4$I[v]*expkEI[v] + E_coefs[5,v]*expK1[v] # Infectious
  
  N[4,1,1,v] <- E_coefs[1,v]*E1$H[v]*expKa[v] + E_coefs[2,v]*E2$H[v]*expks[v] + E_coefs[4,v]*E4$H[v]*expkEI[v] + E_coefs[5,v]*E5$H[v]*expK1[v] + E_coefs[6,v]*expK2[v] # H
  
  N[5,1,1,v] <- E_coefs[1,v]*E1$Ht[v]*expKa[v] + E_coefs[2,v]*E2$Ht[v]*expks[v] + E_coefs[4,v]*E4$Ht[v]*expkEI[v] + E_coefs[5,v]*E5$Ht[v]*expK1[v] + E_coefs[6,v]*E6$Ht[v]*expK2[v] 
  N[5,1,1,v] <- N[5,1,1,v] + E_coefs[7,v]*expK4
  
  N[6,1,1,v]  <- E_coefs[1,v]*E1$D[v]*expKa[v] + E_coefs[2,v]*E2$D[v]*expks[v] + E_coefs[4,v]*E4$D[v]*expkEI[v] + E_coefs[5,v]*E5$D[v]*expK1[v] 
  N[6,1,1,v]  <- N[6,1,1,v] + E_coefs[6,v]*E6$D[v]*expK2[v] + E_coefs[7,v]*E7$D[v]*expK4 + E_coefs[8,v]
  
  N[7,1,1,v]  <- E_coefs[1,v]*E1$R[v]*expKa[v] + E_coefs[2,v]*E2$R[v]*expks[v] + E_coefs[4,v]*E4$R[v]*expkEI[v] + E_coefs[5,v]*E5$R[v]*expK1[v] 
  N[7,1,1,v]  <- N[7,1,1,v] + E_coefs[6,v]*E6$R[v]*expK2[v] + E_coefs[7,v]*E7$R[v]*expK4 + E_coefs[9,v]
  return(N)
}


AbletoMoveKernels <- function(beta, srmc, LDProb, lrmc, ddWeight_sr, LD_weights_nocon, LD_weights_con, ModelType, Example) {
  ddWeight_nocon <- beta*((1-LDProb/100)*ddWeight_sr + LD_weights_nocon)        # Both short and long distance movement, No controls
  ddWeight_LDcon <- beta*(1-LDProb*lrmc/100)*ddWeight_sr + LD_weights_con       # Both short and long distance movement, Reduction of long range distance movement (but becomes short range spread)
  ddWeight_all_con <- srmc*ddWeight_LDcon                                       # Both short and long distance movement, Reduction of transmission (lower beta) and long range movement
  ddWeight_Tcon <- srmc*ddWeight_nocon    # Reduction of transmission (lower beta)   # Both short and long distance movement,  Reduction of transmission for group 1
  ddWeight_Tcon1 <- ddWeight_Tcon                                                    # Both short and long distance movement,  Reduction of transmission for group 2
  if (ModelType == "Use Case") {
    if ((Example == "Measles: Netherlands 2013") || (Example == "H7N9: China 2013")) {
      ddWeight_Tcon1 <- ddWeight_nocon %*% diag(srmc1) #t(t(ddWeight_nocon)*srmc1)  # or Matrix %*% diag(Vector)    # Reduction of transmission (lower beta)
    } 
  }                                             
  AMK <- list(ddWeight_nocon=ddWeight_nocon, ddWeight_LDcon=ddWeight_LDcon, ddWeight_all_con=ddWeight_all_con, ddWeight_Tcon1=ddWeight_Tcon1, ddWeight_Tcon=ddWeight_Tcon) 
  return(AMK)
}

UnabletoMoveKernels <- function(beta, srmc, LDProb, lrmc, local_weight, LD_weights_nocon, LD_weights_con, ModelType, Example) {  
  ddWeight_local_nocon <- beta*((1-LDProb/100)*local_weight + LD_weights_nocon)         # No short range movement, no control 
  ddWeight_local_LDcon <- beta*((1-LDProb*lrmc/100)*local_weight + LD_weights_nocon)    # No short range movement, control long distance
  ddWeight_local_all_con <- srmc*ddWeight_local_LDcon                                   # No short range movement, control long distance and reduce transmission                 
  ddWeight_local_Tcon <- srmc*ddWeight_local_nocon                                      # No short range movement, No long dist control, transmission reduction for group 2
  ddWeight_local_Tcon1 <- ddWeight_local_Tcon                                           # No short range movement, No long dist control, transmission reduction for group 1
  if (ModelType == "Use Case") {
    if ((Example == "Measles: Netherlands 2013") || (Example == "H7N9: China 2013")) {
      ddWeight_local_Tcon1 <- ddWeight_local_nocon %*% diag(srmc1)   
    }
  }
  UMK <- list(ddWeight_local_nocon=ddWeight_local_nocon, ddWeight_local_LDcon=ddWeight_local_LDcon, ddWeight_local_all_con=ddWeight_local_all_con, ddWeight_local_Tcon1=ddWeight_local_Tcon1, ddWeight_local_Tcon=ddWeight_local_Tcon)  
  return(UMK)
}


LDRestrictedKernels <- function(beta, srmc, ddWeight_sr, local_weight, ModelType, Example) {  
  weights <- beta*ddWeight_sr                                                           # short range isotropic spread, no control                 
  weights_con <- srmc*beta*ddWeight_sr                                                  # short range isotropic spread, transmission reduction
  weights_local <- beta*local_weight                                                    # local spread only, no control
  weights_local_con <- srmc*beta*local_weight                                           # local spread only, transmission control for group 2
  weights_local_con1 <- weights_local_con                                               # local spread only, transmission control for group 1
  if (ModelType == "Use Case") {
    if ((Example == "Measles: Netherlands 2013") || (Example == "H7N9: China 2013")) {
      weights_local_con1 <- diag(beta*srmc1,nbox,nbox)  #(beta*local_weight) %*% diag(srmc1)  #local_weight <- diag(1,nbox,nbox)
    }
  }
  RMK <- list(weights=weights, weights_con=weights_con, weights_local=weights_local, weights_local_con1=weights_local_con1, weights_local_con=weights_local_con)
}

PlotTitleAndZlimit <- function(kind, Rout, week) {
  if (kind == "orig") {
    mainTitle <- "Original population" # at calculation resolution"
    zlimit <- c(1,NA)
  }  
  if (kind == "deaths") {
    Ptitle <- "Cumulative Dead, wk "
    maxZval <- max(values(Rout), na.rm=TRUE)
    zlimit <- c(0,maxZval)
  }  
  if (kind == "incidence") {
    Ptitle = "Incidence, wk "
    maxZval <- ceiling(max(values(Rout), na.rm = TRUE))
    zlimit <- c(0,maxZval)
  }
  if (kind == "prevalence") {
    Ptitle = "Prevalence, wk "
    if (max(values(Rout), na.rm=TRUE) >= 0.1) zlimit <- c(0,1)
    if (max(values(Rout), na.rm=TRUE) < 0.1) zlimit <- c(0,0.1)
  }
  if (kind == "vaccinated") {
    Ptitle = "Fraction of potential susceptibles \nvaccinated on wk "
    if (min(values(Rout), na.rm=TRUE) >= 0.5) zlimit <- c(0.5,max(values(Rout), na.rm=TRUE))
    if (min(values(Rout), na.rm=TRUE) < 0.5) zlimit <- c(0,max(values(Rout), na.rm=TRUE))
  }
  if (kind == "vacimmune") {
    Ptitle = "Fraction of potentially susceptibles \nimmune by vaccination on wk "
    if (min(values(Rout), na.rm=TRUE) >= 0.5) zlimit <- c(0.5,max(values(Rout), na.rm=TRUE))
    if (min(values(Rout), na.rm=TRUE) < 0.5) zlimit <- c(0,max(values(Rout), na.rm=TRUE))
  }
  if (kind == "Cumulative cases") {
    Ptitle = "Cumulative cases, wk "
    zlimit <- c(0.3,NA)
  }
  if (kind != "orig") mainTitle <- paste(Ptitle, week, sep="")
  TitleAndZlimits <- list(mainTitle=mainTitle, zlimit=zlimit)
  return(TitleAndZlimits)  
}

RasterToPlotFunc <- function(kind, Rout, Disease, HumanOrAnimal, NByWeek, NByWeekH, IncByWeek, IncByWeekH, week, Nages, nbox) {
  if((Disease=="H7N9 Avian Influenza") && (HumanOrAnimal == "human")) {
    if (kind == "deaths") ImageData <- array(data=NByWeekH[6,1,,,week],dim=c(Nages,nbox))
    if (kind == "incidence") ImageData <- array(data=IncByWeekH[,,week] ,dim=c(Nages,nbox)) 
    if (kind == "Cumulative cases") ImageData <- array(data=rowSums(IncByWeekH[,,1:week]) ,dim=c(Nages,nbox))
    if (kind == "prevalence") {
      SickPop <- array(data=colSums(NByWeekH[2:5,1,,,week]),dim=c(Nages,nbox))
      TotalPop <- array(data=colSums(NByWeekH[1:5,1,,,week])+colSums(NByWeekH[7:9,1,,,week]),dim=c(Nages,nbox))
      ImageData <- SickPop/TotalPop
    }  
    if (kind == "vaccinated") {
      Vacc <- array(data=(NByWeekH[9,1,,,week] + NByWeekH[8,1,,,week]),dim=c(Nages,nbox))
      Sdata <- array(data=(NByWeekH[1,1,,,week] + NByWeekH[9,1,,,week] + NByWeekH[8,1,,,week]),dim=c(Nages,nbox))
      ImageData <- Vacc/Sdata
    }
    if (kind == "vacimmune") {
      Immune <- array(data=NByWeekH[9,1,,,week],dim=c(Nages,nbox))
      Sdata <- array(data=(NByWeekH[1,1,,,week] + NByWeekH[9,1,,,week] + NByWeekH[8,1,,,week]),dim=c(Nages,nbox))
      ImageData <- Immune/Sdata
    }
  } else {
    if (kind == "deaths")  ImageData <- array(data=NByWeek[6,1,,,week],dim=c(Nages,nbox))
    if (kind == "incidence") ImageData <- array(data=IncByWeek[,,week] ,dim=c(Nages,nbox))
    if (kind == "Cumulative cases") {
      if (week>1) ImageData <- array(data=rowSums(IncByWeek[,,1:week]) ,dim=c(Nages,nbox))
      if (week==1) ImageData <- array(data=IncByWeek[,,1] ,dim=c(Nages,nbox))
    }
    if (kind == "prevalence") {
      SickPop <- array(data=colSums(NByWeek[2:5,1,,,week]),dim=c(Nages,nbox))
      TotalPop <- array(data=colSums(NByWeek[1:5,1,,,week])+colSums(NByWeek[7:9,1,,,week]),dim=c(Nages,nbox))
      ImageData <- SickPop/TotalPop
    }  
    if (kind == "vaccinated") {
      Vacc <- array(data=(NByWeek[9,1,,,week] + NByWeek[8,1,,,week]),dim=c(Nages,nbox))
      Sdata <- array(data=(NByWeek[1,1,,,week] + NByWeek[9,1,,,week] + NByWeek[8,1,,,week]),dim=c(Nages,nbox))
      ImageData <- Vacc/Sdata
    }
    if (kind == "vacimmune") {
      Immune <- array(data=NByWeek[9,1,,,week],dim=c(Nages,nbox))
      Sdata <- array(data=(NByWeek[1,1,,,week] + NByWeek[9,1,,,week] + NByWeek[8,1,,,week]),dim=c(Nages,nbox))
      ImageData <- Immune/Sdata
    }
  }
  values(Rout) <- colSums(ImageData)
  return(Rout)
}
