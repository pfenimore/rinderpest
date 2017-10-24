# Rinderpest public functions.r
# GPL v3.0

###################  Read in GIS data ################
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

################### Initialize the population matrix, N ###################
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

################### Roads for Pakistan use case ###################
initRoads <- function(wbb) {
  load("Pakistan/mainPakRoads.RData")
  roads <- crop(mainPakRoads,wbb)
  if (is.null(roads)) return(NULL)
  roads <- gLineMerge(roads)
  return(roads)
}

LongRoads <- function(roads,MinRoadlength) { # Only keep roads that are of significant length, Units are km
  KeptRoads <- NULL
  nLines <- length(roads@lines[[1]]@Lines)
  for (i in 1:nLines){  
    if (LineLength(coordinates(roads@lines[[1]]@Lines[[i]]), longlat=TRUE) > MinRoadlength) {
      KeptRoads <- c(KeptRoads,i)
    }
  }
  return(KeptRoads)
}

FixLines <- function(RL, index_KeptRoads) {
  # 74.29 35.92  end pt of line 2  really 74.28727 35.92151
  # End of 25 is "74.3283002000" "35.9140555000"  Need to connect these - so there is a line across Gilgit!!!!  Extending line 2.
  x <- c(RL[[2]]@coords[,1], 74.30, 74.31, 74.32, 74.3283002000)  #25 will meet 2.
  y <- c(RL[[2]]@coords[,2], 35.92, 35.92, 35.92, 35.9140555000)
  RL[[2]] <- Line(cbind(x,y))
  
  # 74.33 35.92  start pt of line 27 (74.3323203000 35.9197664000) Combine 2 and 27.  
  x <- c(RL[[2]]@coords[,1], RL[[27]]@coords[,1])
  y <- c(RL[[2]]@coords[,2], RL[[27]]@coords[,2])
  RL[[2]] <- Line(cbind(x,y))
  index_KeptRoads <- index_KeptRoads[which(index_KeptRoads!=27)] 
  
  # 37 and 31 should be joined, 37 end is same as 31 start.
  x <- c(RL[[37]]@coords[,1], RL[[31]]@coords[,1])
  y <- c(RL[[37]]@coords[,2], RL[[31]]@coords[,2])
  RL[[37]] <- Line(cbind(x,y))
  index_KeptRoads <- index_KeptRoads[which(index_KeptRoads!=31)] 
  
  #Now 37 and 43 should be joined, end of 43 is start of 37
  x <- c(RL[[43]]@coords[,1], RL[[37]]@coords[,1])
  y <- c(RL[[43]]@coords[,2], RL[[37]]@coords[,2])
  RL[[43]] <- Line(cbind(x,y))
  index_KeptRoads <- index_KeptRoads[which(index_KeptRoads!=37)] 
  
  # 43 doesn't intersect.  Connect with 4
  x <- c(RL[[43]]@coords[,1], RL[[4]]@coords[,1])
  y <- c(RL[[43]]@coords[,2], RL[[4]]@coords[,2])
  RL[[43]] <- Line(cbind(x,y))
  index_KeptRoads <- index_KeptRoads[which(index_KeptRoads!=4)] 
  
  # 47 doesn't connect.  sprintf("%2.10f", RL[[43]]@coords[740,]) "75.7389273000" "35.2995453000"
  x <- c(RL[[47]]@coords[,1], 75.7389273000)
  y <- c(RL[[47]]@coords[,2], 35.2995453000)
  RL[[47]] <- Line(cbind(x,y))
  
  Fixes <- list(RL = RL, index_KeptRoads = index_KeptRoads)
  return(Fixes) 
}  

roadInfo <- function(RL, cattle, KeptRoads) {  #returns a list that for each road has distances from road start and corresponding cell(pixel) number
  print("In roadInfo")
  number_of_roads <- length(KeptRoads)
  road_list <- NULL
  
  withProgress(message = 'Initializing road data', value = 0, {
  for(i in 1:number_of_roads) {     # Determine which cells are on road and determine travel distances from start point to each other road point 
    incProgress(1/number_of_roads)
    index <- KeptRoads[i]
    CellsOnRoad <- cellFromXY(cattle, RL[[index]]@coords)  #Cell numbers corresponding to each point along road
    ToRemove <- which(is.na(CellsOnRoad))
    ConsecDists <- NULL
    npoints <- length(RL[[index]]@coords[,1])
    for(j in 1:(npoints-1)) {
      ConsecDists <- c(ConsecDists,distm(RL[[index]]@coords[j,], RL[[index]]@coords[j+1,]))  # Great sphere distance, results in meters
    }
    distances <- c(0,cumsum(ConsecDists))/1000  # Convert from meters to km
    if (length(ToRemove>0)) {
      CellsOnRoad <- CellsOnRoad[-ToRemove]
      distances <- distances[-ToRemove]
    }  
    road_df <- data.frame(CellsOnRoad,distances)
    road_list[[i]] <- road_df
  }
  })  
  print("Finished roadInfo")
  return(road_list)  #for each road has distance to points and the pixel(cell) the point corresponds to
}

Intersections <- function(RL, index_KeptRoads) {  #*********  Determine All Intersections ***********
  Ipoints = NULL
  ILines = NULL
  nroads <- length(index_KeptRoads)
  Ilist <- vector("list", nroads)
  for (i in 1:(nroads-1)) {
    index = index_KeptRoads[i]  
    for (j in (i+1):length(index_KeptRoads)) {  #All combinations of 2 roads.
      jndex = index_KeptRoads[j]
      test1 <- SpatialLines(list(Lines(RL[[index]], ID=index)))
      test2 <- SpatialLines(list(Lines(RL[[jndex]], ID=jndex)))
      if (gIntersects(test1,test2)) {
        test <- gIntersection(test1,test2)
        loc1 <- which(RL[[index]]@coords[,1]==extent(test)@xmin & RL[[index]]@coords[,2]==extent(test)@ymin)
        loc2 <- which(RL[[jndex]]@coords[,1]==extent(test)@xmin & RL[[jndex]]@coords[,2]==extent(test)@ymin)
        Ilist[[i]] <-  c(Ilist[[i]], loc1)
        Ilist[[j]] <-  c(Ilist[[j]], loc2)
        point <- c(extent(test)@xmin, extent(test)@ymin)
        Ipoints <- rbind(Ipoints, point)
        ILines <- rbind(ILines, c(index, loc1, jndex, loc2))
      }
    }
  }  
  IntersectInfo <- list(Ilist = Ilist, ILines = ILines)
  return(IntersectInfo)
}  

#************* Get intersection road and point on the intersection road: Inputs are a road and it's point on intersection ******
#************* Called by MoveOnRoadSegment *****
IRoadAndPoint <- function(RoadIndex, point_of_Intersection, IntersectInfo, index_KeptRoads) {
  whichInt <- which(IntersectInfo$ILines[,1]==index_KeptRoads[RoadIndex])  #First which intersection, current road could be in colum 1 or 3 (below)
  if (length(whichInt) > 0) {  #i.e. make sure its not integer(0)
    temp <- which(IntersectInfo$ILines[whichInt,2]==point_of_Intersection) 
    if (length(temp) > 0) {
      MetRoad <- IntersectInfo$ILines[whichInt[temp],3]
      R2IntersectPt <- IntersectInfo$ILines[whichInt[temp],4]
    }  
  }  
  whichInt <- which(IntersectInfo$ILines[,3]==index_KeptRoads[RoadIndex])  #
  if (length(whichInt) > 0) { #i.e. make sure its not integer(0)
    temp <- which(IntersectInfo$ILines[whichInt,4]==point_of_Intersection) 
    if (length(temp) > 0) {
      MetRoad <- IntersectInfo$ILines[whichInt[temp],1]
      R2IntersectPt <- IntersectInfo$ILines[whichInt[temp],2]
    }  
  }  
  RIndex <- which(index_KeptRoads==MetRoad)
  return(c(RIndex,R2IntersectPt))
}

#************** Called by road_transport - Determine which cells contagion will move to along a particular road segment and find intersections ****
MoveOnRoadSegment <- function(RoadSegToDo, max_dist, min_dist, IntersectInfo, index_KeptRoads, road_list) {
  RoadIndex <- RoadSegToDo[1]
  #If Road Index is 43, special case it because cattle don't move from the north to the south along this road (per Ben per Rossiter)
  StartPt <- RoadSegToDo[2]
  if ((index_KeptRoads[RoadIndex] == 43) && (StartPt > 2806)) {
    cellsC = NULL
    RoadSegsToDo = NULL
    return(list(RoadSegs=RoadSegsToDo, cellsC=cellsC))
  } 
  DirectionType <- RoadSegToDo[3]
  MovedDist <- RoadSegToDo[4]  #Should be positive
  weight <- RoadSegToDo[5]
  max_dist <- max_dist - MovedDist
  min_dist <- min_dist - MovedDist
  if (min_dist < 0) min_dist <- 0
  
  StartDist <- road_list[[RoadIndex]]$distances[StartPt]
  nIntersections <- length(IntersectInfo$Ilist[[RoadIndex]])  #How
  IntersectNeg <- FALSE
  IntersectPos <- FALSE
  RoadSegsToDo <- NULL
  cellsC <- NULL
  
  if (DirectionType > -1) { # Check for Intersection Points in positive direction
    PosIntPts <- which(IntersectInfo$Ilist[[RoadIndex]] > StartPt) # Positive direction Intersection Points
    if (length(PosIntPts) > 0 ) {
      Pos_Intersect_Pt <- min(IntersectInfo$Ilist[[RoadIndex]][PosIntPts])  # which one is closest
      Pos_intersect_dist <- road_list[[RoadIndex]]$distances[Pos_Intersect_Pt]
      if (Pos_intersect_dist < StartDist+max_dist)  #Is it within the movement range?
        IntersectPos <- TRUE
    }
  }
  if (DirectionType < 1) { # Check for Intersection Points in negative direction
    NegIntPts <- which(IntersectInfo$Ilist[[RoadIndex]] < StartPt) # Negative direction Intersection Points
    if (length(NegIntPts) > 0 ) {
      Neg_Intersect_Pt <- max(IntersectInfo$Ilist[[RoadIndex]][NegIntPts])
      Neg_intersect_dist <- road_list[[RoadIndex]]$distances[Neg_Intersect_Pt]
      if (Neg_intersect_dist > StartDist-max_dist) #Is it within the movement range?
        IntersectNeg <- TRUE
    }  
  }  
  
  if (IntersectNeg)  {
    #print("There is an intersection in the negative direction ")
    DistanceMoved <- StartDist-Neg_intersect_dist
    if (Neg_Intersect_Pt == 1) { # at the start of the road

    } else { #
      #print("in the middle of the road.")
      newDirectionType <- -1 # Go in negative direction only.
      RSTD <- c(RoadIndex, Neg_Intersect_Pt, newDirectionType, DistanceMoved, weight/2)
      RoadSegsToDo <- rbind(RoadSegsToDo, RSTD)  # Will need to continue on current road.
    }
    newDirectionType <- 0     # i.e. both directions
    NewRoadAndPoint <- IRoadAndPoint(RoadIndex, Neg_Intersect_Pt, IntersectInfo, index_KeptRoads)
    RSTD <- c(NewRoadAndPoint, newDirectionType, DistanceMoved, weight/2)  #For Pakistan this works, generally won't.
    RoadSegsToDo <- rbind(RoadSegsToDo, RSTD)   # road met at Intersection
    
    #Get cell locations of boxes between StartDist and intersect_dist - use road_list
    points <- which( (road_list[[RoadIndex]]$distances > Neg_intersect_dist) & (road_list[[RoadIndex]]$distances < (StartDist-min_dist)) )
    ContagionCells <- unique(road_list[[RoadIndex]]$CellsOnRoad[points])
    cellsC <- c(cellsC, ContagionCells)
  }
  
  if (IntersectPos){
    #print("There is an intersection in POSITIVE direction")
    DistanceMoved <- Pos_intersect_dist-StartDist
    if (Pos_Intersect_Pt == length(road_list[[RoadIndex]]$distances) ) {
      newDirectionType <- 1 #Go in positive direction only.
      RSTD <- c(RoadIndex, Pos_Intersect_Pt, newDirectionType, DistanceMoved, weight/2)
      RoadSegsToDo <- rbind(RoadSegsToDo, RSTD) # Will need to continue on current road.
    }
    newDirectionType <- 0 #i.e. both directions
    NewRoadAndPoint <- IRoadAndPoint(RoadIndex, Pos_Intersect_Pt, IntersectInfo, index_KeptRoads)
    RSTD <- c(NewRoadAndPoint, newDirectionType, DistanceMoved, weight/2)
    RoadSegsToDo <- rbind(RoadSegsToDo, RSTD) #road met at Intersection
    
    #Get cell locations of boxes between StartPt and intersect_dist - use road_list
    points <- which((road_list[[RoadIndex]]$distances < Pos_intersect_dist) & (road_list[[RoadIndex]]$distances > (StartDist+min_dist)))
    ContagionCells <- unique(road_list[[RoadIndex]]$CellsOnRoad[points])
    cellsC <- c(cellsC, ContagionCells)
  }
  
  if ( (DirectionType < 1) && (!IntersectNeg))  {
    points <- which((road_list[[RoadIndex]]$distances > StartDist-max_dist) & (road_list[[RoadIndex]]$distances < StartDist-min_dist))
    ContagionCells <- unique(road_list[[RoadIndex]]$CellsOnRoad[points])
    cellsC <- c(cellsC, ContagionCells)
  }
  
  if ( (DirectionType > -1) && (!IntersectPos) ) { 
    points <- which( (road_list[[RoadIndex]]$distances > (StartDist+min_dist)) & (road_list[[RoadIndex]]$distances < (StartDist+max_dist)) )
    ContagionCells <- unique(road_list[[RoadIndex]]$CellsOnRoad[points])
    cellsC <- c(cellsC, ContagionCells)
  }  
  return(list(RoadSegs=RoadSegsToDo, cellsC=cellsC))
}

road_transport <- function(choosencells, RL, KeptRoads, road_list, RoadSpread, IntersectInfo, regioncattle) {
  print("Starting road transport")
  class(IntersectInfo)
  nbox <- dim(regioncattle)[1]*dim(regioncattle)[2]
  road_weights <- array(data=0,dim=(c(nbox,nbox)))
  rMove <- regioncattle
  values(rMove) <- 0
  number_of_roads <- length(KeptRoads)  # KeptRoads are the long roads
  
  for(i in 1:(length(choosencells)))  {
    choosencell <- choosencells[i]
    cell_location <- coordinates(regioncattle)[choosencell,]
    CellOnRoad <- FALSE
    for(j in 1:number_of_roads) {
      Roadcells <- unique(road_list[[j]]$CellsOnRoad)
      if(choosencell%in%Roadcells) {
        line <- KeptRoads[j]
        RoadIndex <- j
        CellOnRoad = TRUE
      }  
    }  # Currently ignoring possibility of being on 2 roads.
    if (CellOnRoad) {
      newmat <- sweep(RL[[line]]@coords,MARGIN=2,cell_location,FUN="-")
      lc_distances <- sqrt((newmat[,1])^2 + (newmat[,2])^2)   # distances between linepoints and cell center. Units don't matter, because
      closest_Pt <- which.min(lc_distances)                   # just looking for minimum
      
      min_dist <- RoadSpread - 7
      max_dist <- RoadSpread + 7
      
      StartPt <- closest_Pt
      DirectionType <- 0        # Both directions
      MovedDist <- 0
      weight <- 1               # Reduced by 2 everytime there's a two-pronged fork in the road
      RoadSegsToDo <- c(RoadIndex, StartPt, DirectionType, MovedDist, weight) # RoadIndex, StartPt, DirectionType, MovedDist, previous road
      notDone <- TRUE
      while(notDone) {
        if (length(RoadSegsToDo) > 5) {
          Segment <- RoadSegsToDo[1,]
        }
        else {
          Segment <- RoadSegsToDo
        }
        MoveResults <- MoveOnRoadSegment(Segment, max_dist, min_dist, IntersectInfo, KeptRoads, road_list)
        weight <- Segment[5]
        NCells <- length(MoveResults$cellsC)
        if (NCells > 0) {
          road_weights[choosencell,MoveResults$cellsC] <- road_weights[choosencell,MoveResults$cellsC] + weight/NCells  
          road_weights[MoveResults$cellsC,choosencell] <- road_weights[MoveResults$cellsC,choosencell] + weight/NCells  
          rMove[MoveResults$cellsC] <- rMove[MoveResults$cellsC] + weight/NCells  
        }
        RoadSegsToDo <- rbind(RoadSegsToDo, MoveResults$RoadSegs)  # Need to check that not done.
        if (length(RoadSegsToDo[,1]) < 2 ) {
          notDone <- FALSE
        }  
        RoadSegsToDo <- RoadSegsToDo[-1,]
      }
    } 
  }  
  return(road_weights)
}    

################# Point to point long distance transport ###########
LD_transport <- function(icell, LDcell, cattle) {
  pstart<-xyFromCell(cattle,icell)
  pend<-xyFromCell(cattle,LDcell)
  xvals <- seq(0,100)/(100/(pend[1]-pstart[1]))+pstart[1]
  yvals <- seq(0,100)/(100/(pend[2]-pstart[2]))+pstart[2]
  linepoints <- cbind(xvals,yvals)
  CellsOnLDtransport <- unique(cellFromXY(cattle,linepoints))
  to_locations <- xyFromCell(cattle,CellsOnLDtransport)
  #need to return to_locations, and CellOnLDtransport
  out=list(to_locations=to_locations,CellsOnLDtransport=CellsOnLDtransport)
  return(out)
}

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

################# Coordinates of middle pixel ################
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

################# Spreading Matrix ################
changeDD=function(r) {
  DD <- distm(coordinates(r),fun=distMeeus)/1000 #Default units are m so change to km
  return(DD)
} #One line

testfunc <- function(val, rand) {
  if (rand < val) {
    return(1)
  }
  else {
    return(0)
  }
}

MakeSpreadStochastic <- function(k, S0, tstep) {
  whichks <- which(S0*(1-exp(-k*tstep)) < tstep) # Roughly equivalent to one person getting sick in a day.
  whichS0gtTstep <- which(S0>tstep)
  indices <- intersect(whichks,whichS0gtTstep)
  ksToCheck <- k[indices]
  NtoCheck <- length(ksToCheck)
  if (NtoCheck > 0) {
    S0ofChecks <- S0[indices]
    rands <- replicate(NtoCheck, runif(1))
    # (randN*dt < S0*(1-exp(-k*dt)))
    newIsChecks <- S0ofChecks*(1-exp(-ksToCheck*tstep))
    randValChecks <- rands*tstep
    PartlyChecked <- pmin(newIsChecks,randValChecks)
    diff1 <- PartlyChecked-randValChecks       # Is 0 if rand < k
    # PartlyChecked[which(diff1==0)] <- log(1-1/S0ofChecks[which(diff1==0)])/tstep
    PartlyChecked[which(diff1==0)] <- -log(1-tstep/S0ofChecks[which(diff1==0)])/tstep
    #largeks <- setdiff(indks,whichS0gtTstep)  
    #PartlyChecked[largeks] <- 1000     ##########  Really want to move everything from S to E, Not quite right!!!!
    diff2 <- PartlyChecked-newIsChecks
    PartlyChecked[which(diff2==0)] <- 0
    k[indices] <- PartlyChecked
  }
  return(k)
}


IsoSpreadMatrix <- function(dd,pixel_size,dl,icell) {
  nbox <- dim(dd)[1]
  dd_weight <- array(data=0,dim=(c(nbox,nbox)))
  # #Put a 15 km limit on how far cattle travel in day  (Should really be 9 or 10, but pixels are too big.)
  # dd.weight[which(dd<15, arr.ind = TRUE)] <- (dl/pixel_size)*exp(-(dd[which(dd<15, arr.ind = TRUE)]/dl))*2*sinh(pixel_size/(2*dl)) #JRM weighting over length of box
  dd_weight <- (dl/pixel_size)*exp(-(dd/dl))*2*sinh(pixel_size/(2*dl)) # JRM weighting over length of box
  dd_weight[which(dd==0, arr.ind = TRUE)] <- 2*(dl/pixel_size)*(1-exp(-1*pixel_size/(2*dl)))
  SRnormVal <- sum(dd_weight[,icell])
  dd_weight <- dd_weight/SRnormVal
  SRnormVal <- sum(dd_weight[,icell])
  dd_weight_sr0 <- dd_weight/SRnormVal  #Normalizing each column - for use when there is no Long Distance movement
  #The idea here is to have the same total transmissibility regardless of dl.
  return(dd_weight_sr0)
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
  kHHt  <- (kHD+kHR)*Rr0[4,OneOrTwo,1]  #Rr1[4]*kD1[4]
  
  ke    <- kD0recov[7,OneOrTwo,1]
  
  klist <- list(kEI=kEI,kIH=kIH,kIR=kIR,kHD=kHD,kHR=kHR,kHtD=kHtD,kHtR=kHtR,kHHt=kHHt,ke=ke)  
  return(klist)
}

Define_Ks <- function(klist, Vrate, Crate, betaRatio, tstep) {
  kEI   <- klist$kEI    #   kEI   <- kD1[2]
  kIH   <- klist$kIH    #   kIH   <- kD1[3]
  kIR   <- klist$kIR    #   kIR   <- (kD1[3]+Rr1[3]*kD1[3])*Rr1[2]
  kHHt  <- klist$kHHt   #   kHHt  <- Rr1[4]*kD1[4]
  kHD   <- klist$kHD    #   kHD   <- kD1[4]
  kHR   <- klist$kHR    #   kHD   <- kD1[4]
  kHtR  <- klist$kHtR   #   kHtR  <- Rr1[6]*kD1[6]
  kHtD  <- klist$kHtD   #   kHtD  <- kD1[6]
  ke    <- klist$ke
  
  K1    <- kIH + kIR 
  K2    <- kHHt + kHD + kHR
  K4    <- kHtD + kHtR
  
  # Avoid redundant eigenvalues. 
  # Since values are not known to an accuracy of even a 10th of a day, seems OK.
  if (K1 == K2)  {
    kIH   <- kIH + 0.001  
    K1    <- kIH + kIR 
  }
  if (K4 == K2)  {
    kHtR  <- kHtR + 0.0005
    K4    <- kHtD + kHtR 
  }
  
  if (kEI == K1) {
    kIR <- kIR + 0.001
    K1    <- kIH + kIR 
  }
  if (kEI == K2) {
    kHHt  <- kHHt + 0.0007
    K2    <- kHHt + kHD + kHR 
  }  
  if (kEI == K4) {
    kHtD  <- kHtD + 0.001
    K4    <- kHtD + kHtR
  }
  
  VaccPerStep <- Vrate*tstep
  
  Ks <- list(K1=K1, K2=K2, K4=K4, kIH=kIH, kHtR=kHtR, kIR=kIR, kHHt=kHHt, kHtD=kHtD, VaccPerStep=VaccPerStep, Crate=Crate, betaRatio=betaRatio)

  return(Ks)
}  

PixelParams <- function(r) {
  middlecoords <- GetMiddleCoords(r) # Units are lat long
  p2 <- c(middlecoords[1] + res(r)[1], middlecoords[2])
  pixel_width <- distMeeus(middlecoords, p2, a=6378137, f=1/298.257223563)/1000 #Units are km 
  p3 <- c(middlecoords[1], middlecoords[2]+res(r)[2])
  pixel_height <- distMeeus(middlecoords, p3, a=6378137, f=1/298.257223563)/1000 #Units are km  
  PixParams <- list(pixel_width=pixel_width, pixel_height=pixel_height)
  return(PixParams)
}

#================================================================ NextN ========================================================
NextN <- function(N, ApproxPixelArea, MM, tstep, klist, Ks, CellsToVacc, CellsToCull, CullFrom, DorS) {
  Ndstates <- dim(N)[1]
  nages=dim(N)[3]
  Nboxes <- dim(N)[4]
  
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
  
  # #-----------------------------------  
  # #StartPop <- (sum(N))
  # LiveN <- sum(N[1:5,1,1,])+sum(N[8:9,1,1,])
  # #----------------------------------- 
  
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
  InitVals <- InitialValues(N)
  
  disc <- NULL
  if (DorS == "Stochastic") {
    disc <- which(InitVals$I < 50)
  }
  
  t <- tstep
  
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
  
  CutOff <- animalCut/26*ApproxPixelArea
  if (DorS == "Deterministic") {
    kEI <- kEIValues(kEIvalue, InitVals, CutOff, Nboxes)
  } else {
    kEI <- rep(kEIvalue, Nboxes)
  }
    
  
  kValues <- function(dd.weight_totI, dd.weight_totHplus, betaRatio, N) {
    Ib <- (dd.weight_totI %*% N[3,1,1,] + dd.weight_totHplus %*% N[4,1,1,])  # dd.weights contain beta!!  
    totalN <- colSums(N[1:5,1,1,])+colSums(N[8:9,1,1,])
    k <-  Ib/totalN #beta*I/N
    k[which(totalN <= 1)] <- Ib[which(totalN <= 1)]                     # Discuss with Carrie?
    return(k)
  }
  k <- kValues(dd.weight_totI, dd.weight_totHplus, betaRatio, N)
  
  # ----------------------------------------------------------
  # # This is where the process of developing a 2nd order Taylor series for "k" begins:
  # dIdt <- dd.weight_totI %*% (kEI*(N[2,1,1,])-(kIR+kIH)*N[3,1,1,]) #(dd.weight_totI %*% N[3,1,1,])
  # dHdt <- dd.weight_totHplus %*% (kIH*(N[3,1,1,])-(kHD+kHR)*N[4,1,1,]) #(dd.weight_totHplus %*% N[4,1,1,]) 
  # # no kHHt???
  # dkdt <- (dIdt+dHdt)*beta_Rind/LiveN  # This ignores the time dependence of N (i.e. death).  Also, ignores any mitigation involving beta.
  # Firstd <- 0.5*dkdt*tstep/k
  # #dEdt <- ((N[3,1,1,])+(N[5,1,1,]))*((N[1,1,1,]) + (N[8,1,1,]))*beta_Rind/sum(N)-kEI_Rind*(N[2,1,1,])
  # dEdt <- k*((N[1,1,1,]) + (N[8,1,1,]))-kEI*(N[2,1,1,])
  # #
  # d2kdt2 <- kEI*dEdt-kIR*dIdt-(kHD+kHR)*dHdt
  # d2kdt2 <- d2kdt2*beta_Rind/LiveN  # This ignores the time dependence of N (i.e. death).  Also, ignores any mitigation involving beta.
  # #print(paste(max(Firstd), min(Firstd), which(Firstd==max(Firstd)), k[136], dkdt[136], N[3,1,1,136], dIdt[136], N[4,1,1,136], dHdt[136])) # N[2,1,1,136], dEdt[136]
  # #print(c("beta",beta_Rind, median(dkdt), median(d2kdt2)))
  # ----------------------------------------------------------
  
  ChooseFunc <- function(Cells, Acts, NperCell){  # When there are not enough resources to vaccinate all locations that were chosen, Pick.
    CellsToComplete <- NULL
    ExtraCell <- NULL
    ExtraActs <- NULL
    indexDecreaseOrder <- order(NperCell[Cells], decreasing = TRUE)  # Do the ones with the highest population first
    CumSumofS <- cumsum(NperCell[Cells[indexDecreaseOrder]])
    if (is.na(CumSumofS[1]) || is.na(Acts)) {
      print(paste("Cells:", Cells))
      print(paste("Acts:", Acts))
      print(paste("CumSumofS", CumSumofS))
      return(list(CellsToComplete=CellsToComplete, ExtraCell=ExtraCell, ExtraActs=ExtraActs))
    }
    if (CumSumofS[1] < Acts) {
      NtoComplete <- max(which(CumSumofS < Acts))
    } else {
      NtoComplete <- 0
    }  
    if (NtoComplete > 0) {
      CellsToComplete <- Cells[indexDecreaseOrder[1:NtoComplete]] # t6events will equal S for these cells, i.e. for stochastic doing all animals
      if (NtoComplete < length(Cells)) {
        ExtraCell <- Cells[indexDecreaseOrder[NtoComplete+1]]
        ExtraActs <- Acts - CumSumofS[NtoComplete]
      }  
    } else {
      ExtraActs <- Acts  # Number of doses/culls for location where incomplete vaccination/culling is occuring.
      ExtraCell <- Cells[indexDecreaseOrder[1]]
    }
    ToDo <- list(CellsToComplete=CellsToComplete, ExtraCell=ExtraCell, ExtraActs=ExtraActs)
    return(ToDo)
  }  
  
  Calc_kVorC <- function(Cells, Acts, NperCell) {  # for boxes calculated deterministically 
    k <- NperCell*0
    Animals <- sum(NperCell[Cells])
    if (Animals > 0) {
      expNegKvcTs <- 1-Acts/Animals # check_rate[CellsToVaccDeterm] <- 1-VaccPerStepDeterm/sum(N[1,1,1,CellsToVaccDeterm])
      if (expNegKvcTs > 0) k[Cells] <- -log(expNegKvcTs)/t   # kV[which(check_rate > 0)] <- -log(check_rate[which(check_rate>0)])/t  
      if (expNegKvcTs <= 0)  k[Cells] <- 5/t
    } 
    return(k)
  }
  
  SetUpforCullorVacc <- function(Cells, disc, VorCPerStep, NperCell){  # NperCell Number of animals of pertanent type as a function of cell
    CellsDiscrete <- NULL
    AnimalsDiscrete <- 0
    if (DorS=="Stochastic") {
      CellsDiscrete   <- intersect(Cells,disc)
      AnimalsDiscrete <- sum(NperCell[CellsDiscrete])
    }
    AmimalsDeterm   <- 0
    CellsDeterm     <- setdiff(Cells,disc)
    AnimalsDeterm   <- sum(NperCell[CellsDeterm])
    totalAnimals <- AnimalsDeterm + AnimalsDiscrete 
    params <- list(CellsDiscrete=CellsDiscrete, AnimalsDiscrete=AnimalsDiscrete, CellsDeterm=CellsDeterm, AnimalsDeterm=AnimalsDeterm, totalAnimals=totalAnimals)
    return(params)
  }
  
  VaccRates <- function(N, CellsToVacc, VaccPerStep, tstep, CutOff, Nboxes, disc) {
    kV <- N[1,1,1,]*0
    CellsToCompletelyVacc <- NULL
    ExtraCell <- NULL
    ExtraDoses <- NULL
    if (!is.null(CellsToVacc)) {
      G1params <- SetUpforCullorVacc(CellsToVacc$G1, disc, VaccPerStep, N[1,1,1,])
      if (!is.null(CellsToVacc$G2)) {
        G2params <- SetUpforCullorVacc(CellsToVacc$G2, disc, VaccPerStep, N[1,1,1,])
      }
      if (!is.null(CellsToVacc$G3)) {
        G3params <- SetUpforCullorVacc(CellsToVacc$G3, disc, VaccPerStep, N[1,1,1,])      }
      if (!is.null(CellsToVacc$G4)) {
        G4params <- SetUpforCullorVacc(CellsToVacc$G4, disc, VaccPerStep, N[1,1,1,])
      }
      
      # Start Calculations #
      if (VaccPerStep < G1params$totalAnimals) {          # *** Do part of G1 ***
        #print(paste("Not enough vacc for G1", VaccPerStep, G1params$totalAnimals))
        if (DorS == "Stochastic") {            
          VaccDiscrete <- VaccPerStep*G1params$AnimalsDiscrete/G1params$totalAnimals
          DiscVacc <- ChooseFunc(G1params$CellsDiscrete, VaccDiscrete, N[1,1,1,])
          CellsToCompletelyVacc <- DiscVacc$CellsToComplete
          ExtraCell <- DiscVacc$ExtraCell
          ExtraDoses <- DiscVacc$ExtraActs
        }  
        if (length(G1params$CellsDeterm) > 0) {
          VaccDeterm <- VaccPerStep*G1params$AnimalsDeterm/G1params$totalAnimals
          DetermVacc  <- ChooseFunc(G1params$CellsDeterm, VaccDeterm, N[1,1,1,])
          VaccCells <- c(DetermVacc$CellsToCompletelyVacc, DetermVacc$ExtraCell)
          kV <- Calc_kVorC(VaccCells, VaccDeterm, N[1,1,1,])
        } else {
          kV <- N[1,1,1,]*0
        }  
      } # End of not enough to vaccinate all of Group 1
      
      if (VaccPerStep >= G1params$totalAnimals) {         # *** Do all of G1 and part of G2 if it exists ***
        if (is.null(CellsToVacc$G2)) { # Only Group 1
          if (DorS == "Stochastic") CellsToCompletelyVacc <- G1params$CellsDiscrete
          if (length(G1params$CellsDeterm) > 0) {
            kV <- Calc_kVorC(G1params$CellsDeterm, G1params$AnimalsDeterm, N[1,1,1,]) # There has to be enough
          } else {
            kV <- N[1,1,1,]*0
          }  
        } else { # Group 2 exists
          if (VaccPerStep <= sum(G1params$totalAnimals, G2params$totalAnimals, na.rm=TRUE)) { # Vaccinate all of G1 and part of G2
            if (DorS == "Stochastic") {
              VaccDiscrete <- sum(VaccPerStep,-G1params$totalAnimals,na.rm=TRUE)*G2params$AnimalsDiscrete/G2params$totalAnimals
              DiscVacc <- ChooseFunc(G2params$CellsDiscrete, VaccDiscrete, N[1,1,1,])
              CellsToCompletelyVacc <- c(G1params$CellsDiscrete, DiscVacc$CellsToComplete)
              ExtraCell <- DiscVacc$ExtraCell
              ExtraDoses <- DiscVacc$ExtraActs
            }
            if (length(G2params$CellsDeterm) > 0) {
              VaccDeterm <- sum(VaccPerStep,-G1params$totalAnimals,na.rm=TRUE)*G2params$AnimalsDeterm/G2params$totalAnimals
              DetermVacc  <- ChooseFunc(G2params$CellsDeterm, VaccDeterm, N[1,1,1,])
              CellsDeterm <- c(G1params$CellsDeterm, DetermVacc$CellsToComplete, DetermVacc$ExtraCell)
              kV <- Calc_kVorC(CellsDeterm, VaccDeterm, N[1,1,1,]) 
            } else {
              if (length(G1params$CellsDeterm) > 0) {
                TotalVaccAvail <- G1params$AnimalsDeterm
                kV <- Calc_kVorC(G1params$CellsDeterm, TotalVaccAvail, N[1,1,1,])
              } else {                        
                kV <- N[1,1,1,]*0
              }
            }  
          }
        }
      }  # End of enough to vaccinate G1, but not all of G2
      
      if (!is.null(CellsToVacc$G2)) {  # Check that G2 exists
        if (VaccPerStep > sum(G1params$totalAnimals, G2params$totalAnimals, na.rm=TRUE)) { # Vaccinate all of G2
          if (!is.null(CellsToVacc$G3)) { # Group 2 and 3 exist
            if (VaccPerStep <= sum(G1params$totalAnimals, G2params$totalAnimals, G3params$totalAnimals, na.rm=TRUE)) {  # Vaccinate all of G2 and part of G3
              if (DorS == "Stochastic") {
                VaccDiscrete <- sum(VaccPerStep,-G1params$totalAnimals,-G2params$totalAnimals,na.rm=TRUE)*G3params$AnimalsDiscrete/G3params$totalAnimals
                DiscVacc    <- ChooseFunc(G3params$CellsDiscrete, VaccDiscrete, N[1,1,1,])
                CellsToCompletelyVacc <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, DiscVacc$CellsToComplete)
                ExtraCell <- DiscVacc$ExtraCell
                ExtraDoses <- DiscVacc$ExtraActs
              }
              if (length(G3params$CellsDeterm) > 0) {
                VaccDeterm <- (VaccPerStep-G1params$totalAnimals-G2params$totalAnimals)*G3params$AnimalsDeterm/G3params$totalAnimals
                DetermVacc  <- ChooseFunc(G3params$CellsDeterm, VaccDeterm, N[1,1,1,])
                CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm, DetermVacc$CellsToComplete, DetermVacc$ExtraCell)
                kV <- Calc_kVorC(CellsDeterm, VaccDeterm, N[1,1,1,]) 
              } else {
                if ( (length(G1params$CellsDeterm) > 0) || (length(G2params$CellsDeterm) > 0)) {
                  CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm)
                  TotalVaccAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, na.rm=TRUE)
                  kV <- Calc_kVorC(CellsDeterm, TotalVaccAvail, N[1,1,1,]) 
                } else {
                  kV <- N[1,1,1,]*0
                }
              }  
            }  
          } else { # No Group 3, vaccinate G1 and G2
            if (DorS == "Stochastic") CellsToCompletelyVacc   <- c(G1params$CellsDiscrete, G2params$CellsDiscrete)
            CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm)
            if (length(CellsDeterm) > 0) {
              TotalVaccAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, na.rm=TRUE)
              kV <- Calc_kVorC(CellsDeterm, TotalVaccAvail, N[1,1,1,])
            } else {  
              kV <- N[1,1,1,]*0
            }
          }  
        }  # End of enough to vaccinate all of G2, but not all of G3
        if (!is.null(CellsToVacc$G3)) { # G3 exists     
          if (VaccPerStep >= sum(G1params$totalAnimals, G2params$totalAnimals, G3params$totalAnimals, na.rm=TRUE)) { # Vaccinate all of G3  
            if (!is.null(CellsToVacc$G4)) { # Group 2, 3 and 4 exist, Vaccinate all of 3 and part of 4
              if (VaccPerStep < sum(G1params$totalAnimals, G2params$totalAnimal, G3params$totalAnimals, G4params$totalAnimals, na.rm=TRUE)) { # Vaccinate all of G3 and part of G4
                if (DorS == "Stochastic") {
                  VaccDiscrete <- sum(VaccPerStep,-G1params$totalAnimals,-G2params$totalAnimals,-G3params$totalAnimals,na.rm=TRUE)*G4params$AnimalsDiscrete/G4params$totalAnimals
                  DiscVacc <- ChooseFunc(G4params$CellsDiscrete, VaccDiscrete, N[1,1,1,])
                  CellsToCompletelyVacc <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, G3params$CellsDiscrete, DiscVacc$CellsToComplete)
                  ExtraCell <- DiscVacc$ExtraCell
                  ExtraDoses <- DiscVacc$ExtraActs
                }
                if (length(G4params$CellsDeterm) > 0) {  
                  VaccDeterm <- sum(VaccPerStep,-G1params$totalAnimals,-G2params$totalAnimals,-G3params$totalAnimals,na.rm=TRUE)*G4params$AnimalsDeterm/G4params$totalAnimals
                  DetermVacc  <- ChooseFunc(G4params$CellsDeterm, VaccDeterm, N[1,1,1,])
                  CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm, DetermVacc$CellsToComplete, DetermVacc$ExtraCell)
                  kV <- Calc_kVorC(CellsDeterm, VaccDeterm, N[1,1,1,]) 
                } else {
                  CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm)
                  if (length(CellsDeterm) > 0) {
                    TotalVaccAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, G3params$AnimalsDeterm, na.rm=TRUE)
                    kV <- Calc_kVorC(CellsDeterm, TotalVaccAvail, N[1,1,1,])
                  } else {  
                    kV <- N[1,1,1,]*0 
                  }  
                }
              }  
              if (VaccPerStep >= sum(G1params$totalAnimals,G2params$totalAnimals,G3params$totalAnimals,G4params$totalAnimals,na.rm=TRUE)) { # Vaccinate all of G4
                if (DorS == "Stochastic") CellsToCompletelyVacc <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, G3params$CellsDiscrete, G4params$CellsDiscrete)
                CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm, G4params$CellsDeterm)
                if (length(CellsDeterm) > 0) {
                  TotalVaccAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, G3params$AnimalsDeterm, G4params$AnimalsDeterm, na.rm=TRUE)
                  if (TotalVaccAvail > 0) {
                    kV <- Calc_kVorC(CellsDeterm, TotalVaccAvail, N[1,1,1,])
                  } else { # No animals to vaccinate
                    kV <- N[1,1,1,]*0 
                  }  
                } else {
                  kV <- N[1,1,1,]*0
                }  
              }
            } else { # Vaccinate all of G3
              if (DorS == "Stochastic") CellsToCompletelyVacc   <- c(CellsToVaccDiscreteG1, CellsToVaccDiscreteG2, CellsToVaccDiscreteG3)
              CellsDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm)
              if (length(CellsDeterm) > 0) {
                TotalVaccAvail <- sum(G1params$AnimalsDeterm,G2params$AnimalsDeterm,G3params$AnimalsDeterm,na.rm=TRUE)
                kV <- Calc_kVorC(CellsDeterm, TotalVaccAvail, N[1,1,1,])
              } else {
                kV <- N[1,1,1,]*0
              }  
            }  
          }
        }  # end if G3 exists
      } # end if G2 exists  
    } 
    VaccInfo <- list(kV=kV, CellsToCompletelyVacc=CellsToCompletelyVacc, ExtraCell=ExtraCell, ExtraDoses=ExtraDoses)
    return(VaccInfo)
  }
  
  kV <- N[1,1,1,]*0   # Initializing
  VaccInfo <- VaccRates(N, CellsToVacc, VaccPerStep, tstep, CutOff, Nboxes, disc) # returns list of kV, CellsToCompletelyVacc, ExtraCell, ExtraDoses
  
  CullRates <- function(N, CellsToCull, Crate, CullFrom, tstep) {
    CellsToCompletelyCull <- NULL
    ExtraCell <- NULL
    ExtraActs <- NULL
    kCS <- N[1,1,1,]*0  
    kCE <- N[1,1,1,]*0
    kCI <- N[1,1,1,]*0
    kCH <- N[1,1,1,]*0
    if (!is.null(CellsToCull)) { 
     if (!is.null(CellsToCull$G1)) {
      CullsPerStep <- Crate*tstep  # Crate is number per day
      
      if (CullFrom == "Honly")  {
        NC <- N[4,1,1,]
      }
      if (CullFrom == "IandH")  {
        NC <- N[3,1,1,] + N[4,1,1,]
      }  
      if (CullFrom == "All")    {
        NC <- N[1,1,1,] + N[2,1,1,] + N[3,1,1,] + N[4,1,1,]
      }  
      
      G1params <- SetUpforCullorVacc(CellsToCull$G1, disc, CullsPerStep, NC)
      if (!is.null(CellsToCull$G2)) {          
        G2params <- SetUpforCullorVacc(CellsToCull$G2, disc, CullsPerStep, NC)
      }  
      if (!is.null(CellsToCull$G3)) {
        G3params <- SetUpforCullorVacc(CellsToCull$G3, disc, CullsPerStep, NC)
      }
      if (!is.null(CellsToCull$G4)) {
        G4params <- SetUpforCullorVacc(CellsToCull$G4, disc, CullsPerStep, NC)
      }
      
      if (CullsPerStep < G1params$totalAnimals) {  # Cull part of G1
        if (DorS == "Stochastic") {    
          ActsDiscrete <- CullsPerStep*G1params$AnimalsDiscrete/G1params$totalAnimals
          DiscCull <- ChooseFunc(G1params$CellsDiscrete, ActsDiscrete, NC)
          CellsToCompletelyCull <- DiscCull$CellsToComplete
          ExtraCell   <- DiscCull$ExtraCell
          ExtraActs  <- DiscCull$ExtraActs
        }  
        if (length(G1params$CellsDeterm) > 0) {
          ActsDeterm <- CullsPerStep*G1params$AnimalsDeterm/G1params$totalAnimals 
          DetermCulls   <- ChooseFunc(G1params$CellsDeterm, ActsDeterm, NC)
          CullCells   <- c(DetermCulls$CellsToComplete, DetermCulls$ExtraCell)
          kC <- Calc_kVorC(CullCells, ActsDeterm, NC)
        } else {
          kC <- N[1,1,1,]*0
        }  
      } # End of not enough culling capacity to cull Group 1
      
      
      if (CullsPerStep >= G1params$totalAnimals) {         # *** Do all of G1 and part of G2 if it exists ***
        if (is.null(CellsToCull$G2)) { # Only Group 1
          if (DorS == "Stochastic") CellsToCompletelyCull <- G1params$CellsDiscrete
          if (length(G1params$CellsDeterm) > 0) {  # Check that some cells in G1 are propating deterministically
            kC <- Calc_kVorC(G1params$CellsDeterm, G1params$AnimalsDeterm, NC) # There has to be enough
          } else {
            kC <- N[1,1,1,]*0
          }  
        } else { # Group 2 exists
          if (CullsPerStep < sum(G1params$totalAnimals, G2params$totalAnimals, na.rm=TRUE)) {     # Do all of G1 and part of G2
            if (DorS == "Stochastic") {
              ActsDiscrete <- sum(CullsPerStep,-G1params$totalAnimals,na.rm=TRUE)*G2params$AnimalsDiscrete/G2params$totalAnimals
              DiscCull <- ChooseFunc(G2params$CellsDiscrete, ActsDiscrete, NC)
              CellsToCompletelyCull <- c(G1params$CellsDiscrete, DiscCull$CellsToComplete)
              ExtraCell <- DiscCull$ExtraCell
              ExtraActs <- DiscCull$ExtraActs
            }
            if (length(G2params$CellsDeterm) > 0) {  # Checking that all cells in G2 are not propagating discretely
              ActsDeterm <- sum(CullsPerStep,-G1params$totalAnimals,na.rm=TRUE)*G2params$AnimalsDeterm/G2params$totalAnimals
              DetermCulls <- ChooseFunc(G2params$CellsDeterm, ActsDeterm, NC)
              CullCells <- c(G1params$CellsDeterm, DetermCulls$CellsToComplete, DetermCulls$ExtraCell)
              kC <- Calc_kVorC(CullCells, ActsDeterm, NC) 
            } else {
              if (length(G1params$CellsDeterm) > 0) { # All cells in G2 are propagating discretely, what about G1?
                ActsDeterm <- G1params$AnimalsDeterm  # By definition there are enough to do this.
                kC <- Calc_kVorC(G1params$CellsDeterm, ActsDeterm, NC) 
              } else {                        
                kC <- N[1,1,1,]*0
              }
            }  
          }
        } 
      }  # End of enough to cull G1, but not all of G2

      if (!is.null(CellsToCull$G2)) {  # G2 exists
        if (CullsPerStep > sum(G1params$totalAnimals, G2params$totalAnimals, na.rm=TRUE)) { # Cull all of G2
          if (!is.null(CellsToCull$G3)) { # Group3 exist
            if (CullsPerStep <= sum(G1params$totalAnimals, G2params$totalAnimals, G3params$totalAnimals, na.rm=TRUE)) {  # Cull all of G2 and part of G3
              if (DorS == "Stochastic") {
                ActsDiscrete <- sum(CullsPerStep,-G1params$totalAnimals,-G2params$totalAnimals,na.rm=TRUE)*G3params$AnimalsDiscrete/G3params$totalAnimals
                DiscCull <- ChooseFunc(G3params$CellsDiscrete, ActsDiscrete, NC)
                CellsToCompletelyCull <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, DiscCull$CellsToComplete)
                ExtraCell <- DiscCull$ExtraCell
                ExtraActs <- DiscCull$ExtraActs
              }
              if (length(G3params$CellsDeterm) > 0) { # Some cells in G3 are propatating deterministically
                ActsDeterm <- sum(CullsPerStep,-G1params$totalAnimals,-G2params$totalAnimals,na.rm=TRUE)*G3params$AnimalsDeterm/G3params$totalAnimals
                DetermCulls <- ChooseFunc(G3params$CellsDeterm, ActsDeterm, NC)
                CullCells <- c(G1params$CellsDeterm, G2params$CellsDeterm, DetermCulls$CellsToComplete, DetermCulls$ExtraCell)
                kC <- Calc_kVorC(CullCells, ActsDeterm, NC) 
              } else {
                if ( (length(G1params$CellsDeterm) > 0) || (length(G2params$CellsDeterm) > 0)) {
                  CellsToCullDeterm <- c(G1params$CellsDeterm,G2params$CellsDeterm)
                  CullsAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, na.rm=TRUE)
                  kC <- Calc_kVorC(CellsToCullDeterm, CullsAvail, NC) 
                } else {
                  kC <- N[1,1,1,]*0
                }
              }  
            }  
          } else { # No Group 3, cull G1 and G2
            if (DorS == "Stochastic") CellsToCompletelyCull   <- c(G1params$CellsDiscrete, G2params$CellsDiscrete)
            CellsToCullDeterm <- c(G1params$CellsDeterm, G2params$CellsDeterm)
            if (length(CellsToCullDeterm) > 0) {
              CullsAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, na.rm=TRUE)
              kC <- Calc_kVorC(CellsToCullDeterm, CullsAvail, NC)
            } else {  
              kC <- N[1,1,1,]*0
            }
          }  
        }  # End of enough to cull all of G2, but not all of G3
        if (!is.null(CellsToCull$G3)) { # G3 exists     
          if (CullsPerStep >= sum(G1params$totalAnimals, G2params$totalAnimals, G3params$totalAnimals, na.rm=TRUE)) { # Cull all of G3  
            if (!is.null(CellsToCull$G4)) { # Groups 2, 3 and 4 exist
              if (CullsPerStep < sum(G1params$totalAnimals, G2params$totalAnimals, G3params$totalAnimals, G4params$totalAnimals, na.rm=TRUE)) { # Cull all of G3 and part of G4
                if (DorS == "Stochastic") {
                  ActsDiscrete <- sum(CullsPerStep,-G1params$totalAnimals,-G2params$totalAnimals,-G3params$totalAnimals,na.rm=TRUE)*G4params$AnimalsDiscrete/G4params$totalAnimals
                  DiscCull    <- ChooseFunc(G4params$CellsDiscrete, ActsDiscrete, NC)
                  CellsToCompletelyCull <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, G3params$CellsDiscrete, DiscCull$CellsToComplete)
                  ExtraCell <- DiscCull$ExtraCell
                  ExtraDoses <- DiscCull$ExtraActs
                }
                if (length(G4params$CellsDeterm) > 0) {  
                  ActsDeterm <- sum(CullsPerStep,-G1params$totalAnimals,-G2params$totalAnimals,-G3params$totalAnimals,na.rm=TRUE)*G4params$AnimalsDeterm/G4params$totalAnimals
                  DetermCulls <- ChooseFunc(G4params$CellsDeterm, ActsDeterm, NC)
                  CullCells <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm, DetermCulls$CellsToComplete, DetermCulls$ExtraCell)
                  kC <- Calc_kVorC(CullCells, ActsDeterm, NC) 
                } else {
                  CullCells <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm)
                  if (length(CullCells) > 0) {
                    CullsAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, G3params$AnimalsDeterm, na.rm=TRUE)
                    kC <- Calc_kVorC(CullCells, CullsAvail, NC)
                  } else {  
                    kC <- N[1,1,1,]*0 
                  }  
                }
              }  
              if (CullsPerStep >= sum(G1params$totalAnimals,G2params$totalAnimals,G3params$totalAnimals,G4params$totalAnimals,na.rm=TRUE)) { # Cull all of G4
                if (DorS == "Stochastic") CellsToCompletelyCull <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, G3params$CellsDiscrete, G4params$CellsDiscrete)
                CullCells <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm, G4params$CellsDeterm)
                if (length(CullCells) > 0) {
                  CullsAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, G3params$AnimalsDeterm, G4params$AnimalsDeterm, na.rm=TRUE)
                  if (CullsAvail > 0) {
                    kC <- Calc_kVorC(CullCells, CullsAvail, NC)
                  } else { # No animals to cull
                    kC <- N[1,1,1,]*0  
                  }  
                } else {
                  kC <- N[1,1,1,]*0
                }  
              }
            } else { # Cull all of G3
              if (DorS == "Stochastic") CellsToCompletelyCull <- c(G1params$CellsDiscrete, G2params$CellsDiscrete, G3params$CellsDiscrete)
              CullCells <- c(G1params$CellsDeterm, G2params$CellsDeterm, G3params$CellsDeterm)
              if (length(CullCells) > 0) {
                CullsAvail <- sum(G1params$AnimalsDeterm, G2params$AnimalsDeterm, G3params$AnimalsDeterm, na.rm=TRUE)
                kC <- Calc_kVorC(CullCells, CullsAvail, NC)
              } else {
                kC <- N[1,1,1,]*0
              }  
            }  
          }
        }  # end if G3 exists
      } # end if G2 exists  
      } else { # G1 does not exist  
        kC <- N[1,1,1,]*0
      }    
    } else {  # CellsToCull is NULL
      kC <- N[1,1,1,]*0
    }
    kCH <- kC
    if (CullFrom == "IandH")  kCI <- kC
    if (CullFrom == "All") {
      kCI <- kCH
      kCS <- kCH
      kCE <- kCH
    }
    CullInfo <- list(kCH=kCH, kCI=kCI, kCE=kCE, kCS=kCS, CellsToCompletelyCull=CellsToCompletelyCull, ExtraCell=ExtraCell, ExtraActs=ExtraActs)
    return(CullInfo)
  }  
  CullInfo <- CullRates(N, CellsToCull, Crate, CullFrom, tstep)
  
  
  ################################ Choose boxes: Vacc or No Vacc  ##############
  # Want everything with kV > 0
  # If kV == 0, but was done and Ng[8,] > 0, take those with k > 0
  # Do not want anything with both k=0 and kV=0.  (dividing by (k+kV))
    v <- which(VaccInfo$kV>0)
      
    nv <- which(VaccInfo$kV==0)
    if (DorS == "Stochastic") {
      nv <- setdiff(nv,disc)
      v <- setdiff(v,disc)
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
      kCS(index) <- kCS(index) + 0.001
    }
    
    if (any(CrateDiffIH == K2-K1))  {   # Would have redundant eigenvalue problem.  The case K1=K2 was taken care of in Define_Ks.
      index <- which(CrateDiffIH == K2-K1)
      kCH(index) <- kCH(index) + 0.001  
    }
    if (any(kCH == K4-K2))  {           # Would have redundant eigenvalue problem.  The case K2=K4 was taken care of in Define_Ks.
      index <- which(kCH == K4-K2)
      kCH(index) <- kCH(index) + 0.001
    }
    if (any(kCI == K4-K1))  {           # Would have redundant eigenvalue problem.  The case K1=K4 was taken care of in Define_Ks.
      index <- which(kCH == K4-K1)
      kCH(index) <- kCI(index) + 0.0005
    }
    
    if (any(CrateDiffEI == K1-kEIvalue)) {  # Would have redundant eigenvalue problem.  The case K1=kEIvalue was taken care of in Define_Ks.
      index <- which(CrateDiffEI == K1-kEIvalue)
      kCE(index) <- kCE(index) + 0.0002
    }
    if (any(CrateDiffEH == K2-kEIvalue)) {  # Would have redundant eigenvalue problem.  The case K2=kEIvalue was taken care of in Define_Ks.
      index <- which(CrateDiffEH == K2-kEIvalue)
      kCE(index) <- kCE(index) + 0.0002
    }
    if (any(kCE == K4-kEIvalue)) {          # Would have redundant eigenvalue problem.  The case K1=kEIvalue was taken care of in Define_Ks.
      index <- which(kCE == K4-kEIvalue)
      kCE(index) <- kCE(index) + 0.0002
    }
    CullRates$kCS <- kCS 
    CullRates$kCE <- kCE
    CullRates$kCI <- kCI
    CullRates$kCH <- kCH
    return(CullRates)
  }
  CullRates <- AvoidRedundantEigenValues(CullInfo, K1, K2, K4, kEIvalue, ke, VaccInfo$kV)
  kCS <- CullRates$kCS
  kCE <- CullRates$kCE
  kCI <- CullRates$kCI
  kCH <- CullRates$kCH
    
  kIndepParams <- function(kEI, Cullrates, K1, K2, K4) {
    kCS <- CullRates$kCS
    kCE <- CullRates$kCE
    kCI <- CullRates$kCI
    kCH <- CullRates$kCH
    
    expkEI  <- exp(-(kEI+kCE)*tstep)
    expK1   <- exp(-(K1+kCI)*tstep) 
    expK2   <- exp(-(K2+kCH)*tstep)
    expK4   <- exp(-K4*tstep)
    
    fkEI1 <- kEI/(K1+kCI-(kEI+kCE))     # Can't be infinite because K1 > 0, and K1+kCI=kEI+kCE was taken care of in AvoidRedundantEigenvalues
    fkEI2 <- kIH/(K2+kCH-(kEI+kCE))     # Can't be infinite because K2 > 0, and K2+kCH=kEI+kCE was taken care of in AvoidRedundantEigenvalues            
    fkEI4 <- kHHt/(K4-(kEI+kCE))        # K4=kEIValue+kCE was taken care of in AvoidRedundantEigenvalues. Used in Eigenvector 4 with check for kEI>0.           
    
    fK1a <- kIH/(K2+kCH-(K1+kCI))       # Can't be infinite because K1 > 0, and K2+kCH=(K1+kCI) was taken care of in AvoidRedundantEigenvalues
    fK1e <- kHHt/(K4-(K1+kCI))          # Can't be infinite because K1 > 0, and K4-K1=kCI was taken care of in AvoidRedundantEigenvalues
    fK4d <- kHHt/(K4-(K2+kCH))          # Can't be infinite because K2 > 0, and K4-K2=kCH was taken care of in AvoidRedundantEigenvalues
    
    kIP <- list(expkEI=expkEI,expK1=expK1,expK2=expK2,expK4=expK4,fkEI1=fkEI1,fkEI2=fkEI2,fkEI4=fkEI4,fK1a=fK1a,fK1e=fK1e,fK4d=fK4d) 
    return(kIP)
  }
  kIP <- kIndepParams(kEI, Cullrates, K1, K2, K4)
  
  kDepParams <- function(k,ke,kV,CullRates,K1,K2,K4) {
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
  kDP <- kDepParams(k,ke,VaccInfo$kV,CullRates,K1,K2,K4)
  
  EigenVector1wVacc <- function(Zs, Os, fkV, ke, Ka, fK1, fK2, fK3, fK5, CullRates, kHD, kHtD, kIR, kHR, kHtR){
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
    if (any(!is.finite(Vim))) {
      print("VimE1 not finite")
      print(Ka[which(!is.finite(Vim))])
      print(fkV[which(!is.finite(Vim))])
    }
    
    kgt0 <- which(k!=0)  # if k=0  E, I, H, Ht , R and R are zero  
    E[kgt0]   <- fK1[kgt0]*(1+fkV[kgt0])  # function of k, length is Nboxes
    I[kgt0]   <- fK1[kgt0]*fK2[kgt0]*(1+fkV[kgt0])
    H[kgt0]   <- fK1[kgt0]*fK2[kgt0]*fK3[kgt0]*(1+fkV[kgt0])  
    Ht[kgt0]  <- fK1[kgt0]*fK2[kgt0]*fK3[kgt0]*fK5[kgt0]*(1+fkV[kgt0])
    D[kgt0]   <- -(kCS[kgt0]*S[kgt0] + kCE[kgt0]*E[kgt0] + kCI[kgt0]*I[kgt0] + (kHD+kCH[kgt0])*H[kgt0] + kHtD*Ht[kgt0])/(Ka[kgt0])  
    R[kgt0]   <- -(kIR*I[kgt0] + kHR*H[kgt0] + kHtR*Ht[kgt0])/(Ka[kgt0])
    E1 <- list(S=S,E=E,I=I,H=H,Ht=Ht,D=D,R=R,Vg=Vg,Vim=Vim)
    return(E1)
  }  
  if (length(v) > 0) {
    E1 <- EigenVector1wVacc(Zs, Os, kDP$fkV, ke, kDP$Ka, kDP$fK1, kDP$fK2, kDP$fK3, kDP$fK5, CullRates, kHD, kHtD, kIR, kHR, kHtR)
  }
  
  EigenVector2wVacc <- function(Zs, Os, k, kEI, kIH, kIR, kHHt, kHD, kHR, kHtD, kHtR, CullRates, ke, ks, K1, K2, K4) {
    kCS <- CullRates$kCS
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
    kgt0 <- which(k!=0)  # if k=0  E, I, H, Ht , R and R are zero  
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
  if (length(v) > 0) {
  E2 <- EigenVector2wVacc(Zs, Os, k, kEI, kIH, kIR, kHHt, kHD, kHR, kHtD, kHtR, CullRates, ke, kDP$ks, K1, K2, K4)
  }
  # Eigenvector 3 - a place holder for when treatment in I is in place
    
  EigenVector4 <- function(Zs, kEI, fkEI1, fkEI2, fkEI4, CullRates, kHD, kHtD, kIR, kHR, kHtR)  {
    kCS <- CullRates$kCS
    kCE <- CullRates$kCE
    kCI <- CullRates$kCI
    kCH <- CullRates$kCH

    E   <- Zs
    I   <- Zs
    H   <- Zs
    Ht  <- Zs
    D   <- Zs
    R   <- Zs
    
    kEIgt0 <- which(kEI > 0)
    E[kEIgt0]   <- 1  # ones           # B
    I[kEIgt0]   <- fkEI1[kEIgt0]*E[kEIgt0]           # C
    H[kEIgt0]   <- I[kEIgt0]*fkEI2[kEIgt0]           # D
    Ht[kEIgt0]  <- H[kEIgt0]*fkEI4[kEIgt0]           # E
    D[kEIgt0]   <- -1*(kCE[kEIgt0]*E[kEIgt0] + kCI[kEIgt0]*I[kEIgt0] + (kHD+kCH[kEIgt0])*H[kEIgt0] + kHtD*Ht[kEIgt0])/(kEI[kEIgt0] + kCE[kEIgt0])
    R[kEIgt0]   <- -1*(kIR*I[kEIgt0] + kHR*H[kEIgt0] + kHtR*Ht[kEIgt0])/(kEI[kEIgt0] + kCE[kEIgt0]) 
    E4 <- list(E=E,I=I,H=H,Ht=Ht,D=D,R=R)
    return(E4)
  }
  E4 <- EigenVector4(Zs, kEI, kIP$fkEI1, kIP$fkEI2, kIP$fkEI4, CullRates, kHD, kHtD, kIR, kHR, kHtR)

  EigenVector5 <- function(Os, fK1a, fK1e, K1, CullRates, kHD, kHtD, kIR, kHR, kHtR) {
    kCS <- CullRates$kCS
    kCE <- CullRates$kCE
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
  E5 <- EigenVector5(Os, kIP$fK1a, kIP$fK1e, K1, CullRates, kHD, kHtD, kIR, kHR, kHtR)

  EigenVector6 <- function(Os,fK4d, K2, CullRates, kHD, kHtD, kHR, kHtR) {
    H   <- Os  
    Ht  <- fK4d
    D   <- -1*((kHD+kCH)*H + kHtD*Ht)/(K2+kCH)
    R   <- -1*(kHR*H + kHtR*Ht)/(K2+kCH)
    E6 <- list(H=H,Ht=Ht,D=D,R=R)
    return(E6)
  }  
  E6 <- EigenVector6(Os, kIP$fK4d, K2, CullRates, kHD, kHtD, kHR, kHtR)  

  EigenVector7 <- function(Os, K4, kHtD, Nboxes) {  
    Ht   <- Os  #zeta, K4
    D   <- rep(-kHtD/K4, Nboxes)
    R   <- rep(-kHtR/K4, Nboxes)
    E7 <- list(Ht=Ht,D=D,R=R)
    return(E7)
  }
  E7 <- EigenVector7(Os, K4, kHtD, Nboxes)
  
  N_New_wVacc <- function(v, Ka, ks, N, E1, E2, E3, E4, E5, E6, E7, InitVals, expkEI, expK1, expK2, expK4) {

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
  
  if (length(v) > 0) {
    Nv <- N_New_wVacc(v, kDP$Ka, kDP$ks, N, E1, E2, E3, E4, E5, E6, E7, InitVals, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4)
  }
  
  Ntemp <- N  
  if (length(v) > 0) {
    Ntemp[,,,v] <- Nv[,,,v]
  }
  
  EigenVector1NoVacc <- function(k, kEI, kIH, kHHt, K1, K2, K4, CullRates) {
    kCS <- CullRates$kCS
    kCE <- CullRates$kCE
    kCI <- CullRates$kCI
    kCH <- CullRates$kCH
    
    fk1 <- k/(kEI+kCE-(k+kCS))
    fk2 <- kEI/(K1+kCI-(k+kCS))
    fk3 <- kIH/(K2+kCH-(k+kCS))
    fk5 <- kHHt/(K4-(k+kCS))

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
    E1N <- list(S=S,E=E,I=I,H=H,Ht=Ht,D=D,R=R)
    return(E1N)
  }  
  E1N <- EigenVector1NoVacc(k, kEI, kIH, kHHt, K1, K2, K4, CullRates)
    
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
    if (any(!is.finite(N[1,1,1,nv]))) print("N_NewNoVacc: Problem with S")

    N[2,1,1,nv] <- E_coefs[1,nv]*E1N$E[nv]*expk[nv] + E_coefs[2,nv]*expkEI[nv]  # Exposed
    if (any(!is.finite(N[2,1,1,nv]))) print("Problem with E")
    
    N[3,1,1,nv] <- E_coefs[1,nv]*E1N$I[nv]*expk[nv] + E_coefs[2,nv]*E4$I[nv]*expkEI[nv] + E_coefs[3,nv]*expK1[nv] # Infectious
    if (any(!is.finite(N[3,1,1,nv]))) print("Problem with I")

    N[4,1,1,nv] <- E_coefs[1,nv]*E1N$H[nv]*expk[nv] + E_coefs[2,nv]*E4$H[nv]*expkEI[nv] + E_coefs[3,nv]*E5$H[nv]*expK1[nv] + E_coefs[4,nv]*expK2[nv] # H
    if (any(!is.finite(N[4,1,1,nv]))) print("Problem with H")

    N[5,1,1,nv]  <- E_coefs[1,nv]*E1N$Ht[nv]*expk[nv] + E_coefs[2,nv]*E4$Ht[nv]*expkEI[nv] + E_coefs[3,nv]*E5$Ht[nv]*expK1[nv]
    N[5,1,1,nv]  <- N[5,1,1,nv]+ E_coefs[4,nv]*E6$Ht[nv]*expK2[nv] + E_coefs[5,nv]*expK4
    if (any(!is.finite(N[5,1,1,nv]))) print("Problem with Ht")
    
    N[6,1,1,nv]  <- E_coefs[1,nv]*E1N$D[nv]*(expk[nv]-1) + E_coefs[2,nv]*E4$D[nv]*(expkEI[nv]-1) + E_coefs[3,nv]*E5$D[nv]*(expK1[nv]-1) 
    N[6,1,1,nv] <- N[6,1,1,nv] + E_coefs[4,nv]*E6$D[nv]*(expK2[nv]-1) + E_coefs[5,nv]*E7$D[nv]*(expK4-1) + IV$D[nv]
    if (any(!is.finite(N[6,1,1,nv]))) print("Problem with D")
    
    N[7,1,1,nv]  <- E_coefs[1,nv]*E1N$R[nv]*(expk[nv]-1) + E_coefs[2,nv]*E4$R[nv]*(expkEI[nv]-1) + E_coefs[3,nv]*E5$R[nv]*(expK1[nv]-1) 
    N[7,1,1,nv]  <- N[7,1,1,nv] + E_coefs[4,nv]*E6$R[nv]*(expK2[nv]-1) + E_coefs[5,nv]*E7$R[nv]*(expK4-1) + IV$R[nv]
    
    if (any(!is.finite(N[7,1,1,nv]))) print("Problem with R")

    Ng8pos_nv <- intersect(nv, which(N[8,1,1,]>0))
    N[8,1,1,Ng8pos_nv]  <- N[8,1,1,Ng8pos_nv]*exp(-ke*tstep)
    N[9,1,1,Ng8pos_nv]  <- N[9,1,1,Ng8pos_nv] + N[8,1,1,Ng8pos_nv]*(1-exp(-ke*tstep))
    if (any(!is.finite(N[9,1,1,Ng8pos_nv]))) print("N[9,1,1,Ng8pos_nv] not finite")
    return(N)
  }
  Nnv <- N_NewNoVacc(nv, N, Ndstates, Nboxes, E1N, E4, E5, E6, E7, k, CullRates, tstep, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4, InitVals, ke)
  Ntemp[,,,nv] <- Nnv[,,,nv]

  k2 <- kValues(dd.weight_totI, dd.weight_totHplus, betaRatio, Ntemp)
  k <- (k+k2)/2

  kDP <- kDepParams(k,ke,VaccInfo$kV,CullRates,K1,K2,K4)
  
  if (length(v) > 0) {
  E1 <- EigenVector1wVacc(Zs, Os, kDP$fkV, ke, kDP$Ka, kDP$fK1, kDP$fK2, kDP$fK3, kDP$fK5, CullRates, kHD, kHtD, kIR, kHR, kHtR)
  E2 <- EigenVector2wVacc(Zs, Os, k, kEI, kIH, kIR, kHHt, kHD, kHR, kHtD, kHtR, CullRates, ke, kDP$ks, K1, K2, K4)
  Nv <- N_New_wVacc(v, kDP$Ka, kDP$ks, N, E1, E2, E3, E4, E5, E6, E7, InitVals, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4)
  N[,,,v] <- Nv[,,,v]
  }
  
  E1N <- EigenVector1NoVacc(k, kEI, kIH, kHHt, K1, K2, K4, CullRates)
  
  Nnv <- N_NewNoVacc(nv, N, Ndstates, Nboxes, E1N, E4, E5, E6, E7, k, CullRates, tstep, kIP$expkEI, kIP$expK1, kIP$expK2, kIP$expK4, InitVals, ke)
  N[,,,nv] <- Nnv[,,,nv]
  
  tauLeepCalc <- function(InitVals,VaccInfo,k,kEI,kIH,kHR,kHD,dt,disc) {
    S0s <- InitVals$S
    E0s <- InitVals$E
    I0s <- InitVals$I
    H0s <- InitVals$H
    R0s <- InitVals$R
    D0s <- InitVals$D
    Vg0s <- InitVals$Vg
    Vim0s <- InitVals$Vim
    
    # Deal with fractional animals at end, for sampling Poisson use a value of 1 even if less
    S0s[which(S0s<1)] <- ceiling(S0s[which(S0s<1)])
    E0s[which(E0s<1)] <- ceiling(E0s[which(E0s<1)])
    I0s[which(I0s<1)] <- ceiling(I0s[which(I0s<1)])
    H0s[which(H0s<1)] <- ceiling(H0s[which(H0s<1)])
    Vg0s[which(Vg0s<1)] <- ceiling(Vg0s[which(Vg0s<1)])
    
    t1events <- rpois(length(disc),k[disc]*S0s[disc]*dt)    # S -> E; returns integers, t is for transition, 
    t2events <- rpois(length(disc),kEI[disc]*E0s[disc]*dt)  # E -> I
    t3events <- rpois(length(disc),kIH*I0s[disc]*dt)        # I -> H
    t4events <- rpois(length(disc),kHR*H0s[disc]*dt)        # H -> R
    t5events <- rpois(length(disc),kHD*H0s[disc]*dt)        # H -> D
    t7events <- rpois(length(disc),ke*Vg0s[disc]*dt)        # Vg -> Vim
    t9events <- rpois(length(disc),k[disc]*Vg0s[disc]*dt)   # Vg -> E
    
    # Assure that number of Susceptibles is non-negative.
    discSs        <- S0s[disc]
    notOKindices  <- which(discSs < t1events)
    t1events[notOKindices] <- discSs[notOKindices]
    # Assure that Vg is non-negative
    discVgs       <- Vg0s[disc]
    notOKindices  <- which(discVgs < (t7events+t9events))
    discks        <- k[disc]
    t7events[notOKindices] <- ke/(ke+discks[notOKindices])*(discVgs[notOKindices])
    t9events[notOKindices] <- discks[notOKindices]/(ke+discks[notOKindices])*(discVgs[notOKindices])
    # Assure that E is non-negative
    discEs        <- E0s[disc]
    notOKindices  <- which(discEs + t1events + t9events < t2events)
    t2events[notOKindices] <- discEs[notOKindices] + t1events[notOKindices] + t9events[notOKindices]
    # Assure I is non-negative
    discIs        <- I0s[disc]
    notOKindices  <- which((discIs+t2events) < t3events)
    t3events[notOKindices] <- discIs[notOKindices] + t2events[notOKindices]
    # Assure H is non-negative
    discHs        <- H0s[disc]
    notOKindices  <- which((discHs+t3events) < (t4events+t5events))
    t4events[notOKindices] <- kHR/(kHR+kHD)*(discHs[notOKindices]+t3events[notOKindices])
    t5events[notOKindices] <- kHD/(kHR+kHD)*(discHs[notOKindices]+t3events[notOKindices])
    
    Snew    <- discSs - t1events
    Vgnew   <- discVgs - t7events - t9events
    Vimnew  <- Vim0s[disc] + t7events
    Enew    <- discEs + t1events + t9events - t2events 
    Inew    <- discIs + t2events - t3events
    Hnew    <- discHs + t3events - t4events - t5events
    Rnew    <- R0s[disc] + t4events
    Dnew    <- D0s[disc] + t5events
    
    incidenceCheck <- -1*(Snew + Vgnew + Vimnew)
    if (any((incidenceCheck>0)&(incidenceCheck<1))) {
      print("Problem with stochastic incidence")
      index <- which((incidenceCheck>0)&(incidenceCheck<1))
      print(paste(incidenceCheck[index],discSs[index],discVgs[index]))
    }  

    tauResults <- list(Snew=Snew, Enew=Enew, Inew=Inew, Hnew=Hnew, Rnew=Rnew, Dnew=Dnew, Vgnew=Vgnew, Vimnew=Vimnew)
    return(tauResults)
  }
  
  if (DorS == "Stochastic") {
    tauResults <- tauLeepCalc(InitVals,VaccInfo,k,kEI,kIH,kHR,kHD,tstep,disc)
    N[1,1,1,disc] <- tauResults$Snew
    N[2,1,1,disc] <- tauResults$Enew
    N[3,1,1,disc] <- tauResults$Inew
    N[4,1,1,disc] <- tauResults$Hnew
    N[6,1,1,disc] <- tauResults$Dnew
    N[7,1,1,disc] <- tauResults$Rnew
    N[8,1,1,disc] <- tauResults$Vgnew
    N[9,1,1,disc] <- tauResults$Vimnew

    # **** Do Vaccination ****
    if (!is.null(VaccInfo$CellsToCompletelyVacc)) {
      N[8,1,1,VaccInfo$CellsToCompletelyVacc] <- N[8,1,1,VaccInfo$CellsToCompletelyVacc] + N[1,1,1,VaccInfo$CellsToCompletelyVacc]
      N[1,1,1,VaccInfo$CellsToCompletelyVacc] <- N[1,1,1,VaccInfo$CellsToCompletelyVacc] - N[1,1,1,VaccInfo$CellsToCompletelyVacc]
    }  
    if (!is.null(VaccInfo$ExtraCell)) {
      N[8,1,1,VaccInfo$ExtraCell] <- N[8,1,1,VaccInfo$ExtraCell] + VaccInfo$ExtraDoses
      N[1,1,1,VaccInfo$ExtraCell] <- N[1,1,1,VaccInfo$ExtraCell] - VaccInfo$ExtraDoses
    }  
    
    # **** Do culling ****
    BC <- N[1,1,1,disc]
    if (!is.null(CullInfo$CellsToCompletelyCull)) {
      if (CullFrom == "all") {
        N[6,1,1,CullInfo$CellsToCompletelyCull] <- N[6,1,1,CullInfo$CellsToCompletelyCull] + colSums(N[1:4,1,1,CullInfo$CellsToCompletelyCull])
        N[1,1,1,CullInfo$CellsToCompletelyCull] <- N[1,1,1,CullInfo$CellsToCompletelyCull] - N[1,1,1,CullInfo$CellsToCompletelyCull]
        N[2,1,1,CullInfo$CellsToCompletelyCull] <- N[2,1,1,CullInfo$CellsToCompletelyCull] - N[2,1,1,CullInfo$CellsToCompletelyCull]
        N[3,1,1,CullInfo$CellsToCompletelyCull] <- N[3,1,1,CullInfo$CellsToCompletelyCull] - N[3,1,1,CullInfo$CellsToCompletelyCull]
        N[4,1,1,CullInfo$CellsToCompletelyCull] <- N[4,1,1,CullInfo$CellsToCompletelyCull] - N[4,1,1,CullInfo$CellsToCompletelyCull]
      }
      if (CullFrom == "IandH") {
        N[6,1,1,CullInfo$CellsToCompletelyCull] <- N[6,1,1,CullInfo$CellsToCompletelyCull] + colSums(N[3:4,1,1,CullInfo$CellsToCompletelyCull])
        N[3,1,1,CullInfo$CellsToCompletelyCull] <- N[3,1,1,CullInfo$CellsToCompletelyCull] - N[3,1,1,CullInfo$CellsToCompletelyCull]
        N[4,1,1,CullInfo$CellsToCompletelyCull] <- N[4,1,1,CullInfo$CellsToCompletelyCull] - N[4,1,1,CullInfo$CellsToCompletelyCull]
      }
      if (CullFrom == "Honly") {
        N[6,1,1,CullInfo$CellsToCompletelyCull] <- N[6,1,1,CullInfo$CellsToCompletelyCull] + N[4,1,1,CullInfo$CellsToCompletelyCull]
        N[4,1,1,CullInfo$CellsToCompletelyCull] <- N[4,1,1,CullInfo$CellsToCompletelyCull] - N[4,1,1,CullInfo$CellsToCompletelyCull]
      }
    }  
    if (!is.null(CullInfo$ExtraCell)) {
      if (CullFrom == "Honly") {
          N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + CullInfo$ExtraActs
          N[4,1,1,CullInfo$ExtraCell] <- N[4,1,1,CullInfo$ExtraCell] - CullInfo$ExtraActs  
      } 
      if (CullFrom == "IandH") {
        if (N[4,1,1,CullInfo$ExtraCell] < CullInfo$ExtraActs) {
          N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + N[4,1,1,CullInfo$ExtraCell]
          N[4,1,1,CullInfo$ExtraCell] <- N[4,1,1,CullInfo$ExtraCell] - N[4,1,1,CullInfo$ExtraCell]
          N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + (CullInfo$ExtraActs - N[4,1,1,CullInfo$ExtraCell])
          N[3,1,1,CullInfo$ExtraCell] <- N[3,1,1,CullInfo$ExtraCell] - (CullInfo$ExtraActs - N[4,1,1,CullInfo$ExtraCell])
        } else {
          N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + CullInfo$ExtraActs
          N[4,1,1,CullInfo$ExtraCell] <- N[4,1,1,CullInfo$ExtraCell] - CullInfo$ExtraActs
        }
      }
      if (CullFrom == "All") {
        if (sum(N[2:4,1,1,CullInfo$ExtraCell]) < CullInfo$ExtraActs) {
          N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + sum(N[2:4,1,1,CullInfo$ExtraCell])
          leftover <- CullInfo$ExtraActs - sum(N[2:4,1,1,CullInfo$ExtraCell])
          N[4,1,1,CullInfo$ExtraCell] <- 0
          N[3,1,1,CullInfo$ExtraCell] <- 0
          N[2,1,1,CullInfo$ExtraCell] <- 0
          N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + leftover
          N[1,1,1,CullInfo$ExtraCell] <- N[1,1,1,CullInfo$ExtraCell] - leftover
        } else {
          if (sum(N[3:4,1,1,CullInfo$ExtraCell]) < CullInfo$ExtraActs) {
              N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + sum(N[3:4,1,1,CullInfo$ExtraCell])
              leftover <- CullInfo$ExtraActs - sum(N[3:4,1,1,CullInfo$ExtraCell])
              N[4,1,1,CullInfo$ExtraCell] <- 0
              N[3,1,1,CullInfo$ExtraCell] <- 0
              N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + leftover
              N[2,1,1,CullInfo$ExtraCell] <- N[2,1,1,CullInfo$ExtraCell] - leftover
          } else {
            if (sum(N[4,1,1,CullInfo$ExtraCell]) < CullInfo$ExtraActs)  {
              N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + N[4,1,1,CullInfo$ExtraCell]
              leftover <- (CullInfo$ExtraActs - N[4,1,1,CullInfo$ExtraCell])
              N[4,1,1,CullInfo$ExtraCell] <- 0
              N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + leftover
              N[3,1,1,CullInfo$ExtraCell] <- N[3,1,1,CullInfo$ExtraCell] - leftover
            } else{
              N[6,1,1,CullInfo$ExtraCell] <- N[6,1,1,CullInfo$ExtraCell] + CullInfo$ExtraActs
              N[4,1,1,CullInfo$ExtraCell] <- N[4,1,1,CullInfo$ExtraCell] - CullInfo$ExtraActs
            }
          }
        }
      }
    }  
  }
  
  num_acc_ind <- which(N < 0, arr.ind=TRUE)
  N[num_acc_ind] = 0

  S0 <- InitVals$S
  Vg0 <- InitVals$Vg
  # OK for continuous model
  incidentbybox <- (S0+Vg0)*(1-exp(-k*tstep)) # NOT QUITE CORRECT WHEN VACCINATING (Ignores Vg -> Vm, probably OK because slow compared to time step)!!!!
  # For discrete model 
  if (DorS == "Stochastic") {
    incidentbybox[disc] <- (InitVals$S[disc] - BC) + ((InitVals$Vg[disc]+InitVals$Vim[disc]) - (N[8,1,1,disc]+N[9,1,1,disc]))
    incidentbybox[which(incidentbybox<1e-12)] <- 0
  }  

  results <- list(N=N,incidentbybox=incidentbybox)
  return(results)
}


################################################################# episim #############################################################
episim <- function(N,r,dl,dd,TimeParams,beta,betaRatio,Rr0,kD0,DetectSick,srmcd,srmc,VacParams, CullParams, icell, LDMoveParams, DorS, UsedParams) {
  
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

  # Calculate pixel parameters
  PixParams <- PixelParams(r)
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
  MiddleCell=cellFromXY(r,middle_location)
  ddWeight_sr <- IsoSpreadMatrix(dd, pixel_size, dl, MiddleCell)  # Isotropic nearby spread
  local_weight <- diag(1,nbox,nbox) # No isotropic spread only get animals in same cell sick; for the really sick who are unlikely to move w/o help
  
  LD_weights <- array(data=0,dim=(c(nbox,nbox))) # 
  LD_weights_con <- array(data=0,dim=(c(nbox,nbox))) 
  LD_weights_nocon <- array(data=0,dim=(c(nbox,nbox)))
  
  if (LDType=="P2P") {
    LDcell <- LDMoveParams$Lcellnum
    LD_weights_nocon[LDcell,] <- LDProb/100
    LD_weights_con[LDcell,] <- LDProb*lrmc/100
  }
  
  init_choosencells <- seq(1:nbox)
  if (LDType == "roads")  { 
      LD_weights <- road_transport(init_choosencells, RL, index_KeptRoads, road_list, RoadSpread, IntersectInfo, r)
      LD_weights_nocon <- LD_weights*(LDProb/100)*nbox/sum(LD_weights)    # nocon - no control in this case of directed "long distance" motion
      LD_weights_con <- lrmc*LD_weights*(LDProb/100)*nbox/sum(LD_weights)
  }
  
  AbletoMoveKernels <- function(beta, srmc, LDProb, lrmc, ddWeight_sr, LD_weights_nocon) {
    ddWeight_nocon <- beta*((1-LDProb/100)*ddWeight_sr + LD_weights_nocon)        # No control of movement or transmission
    ddWeight_LDcon <- beta*(1-LDProb*lrmc/100)*ddWeight_sr + LD_weights_con       # Reduction of long range distance movement (but becomes short range spread)
    ddWeight_all_con <- srmc*ddWeight_LDcon                                       # Reduction of transmission (lower beta) and long range movement
    ddWeight_Tcon <- srmc*ddWeight_nocon                                          # Reduction of transmission (lower beta)
    AMK <- list(ddWeight_nocon=ddWeight_nocon, ddWeight_LDcon=ddWeight_LDcon, ddWeight_all_con=ddWeight_all_con, ddWeight_Tcon=ddWeight_Tcon)  
    return(AMK)
  }
  
  UnabletoMoveKernels <- function(beta, srmc, LDProb, lrmc, local_weight, LD_weights_nocon) {  # Animals or People too sick to spread disease themselves.
    ddWeight_local_nocon <- beta*((1-LDProb/100)*local_weight + LD_weights_nocon)
    ddWeight_local_LDcon <- beta*((1-LDProb*lrmc/100)*local_weight + LD_weights_nocon)
    ddWeight_local_all_con <- beta*srmc*ddWeight_local_LDcon
    ddWeight_local_Tcon <- beta*srmc*ddWeight_local_nocon
    UMK <- list(ddWeight_local_nocon=ddWeight_local_nocon, ddWeight_local_LDcon=ddWeight_local_LDcon, ddWeight_local_all_con=ddWeight_local_all_con, ddWeight_local_Tcon=ddWeight_local_Tcon)  
    return(UMK)
  }
  
  LDRestrictedKernels <- function(beta, srmc, ddWeight_sr, local_weight) {  # Long distance movement is not allowed for some groups
    weights <- beta*ddWeight_sr
    weights_con <- srmc*beta*ddWeight_sr
    weights_local <- beta*local_weight
    weights_local_con <- srmc*beta*local_weight
    RMK <- list(weights=weights, weights_con=weights_con, weights_local=weights_local, weights_local_con=weights_local_con)
  }
  
  # Precompute Matrices related to disease spread to save time 
  AMK <- AbletoMoveKernels(beta, srmc, LDProb, lrmc, ddWeight_sr, LD_weights_nocon)
  UMK <- UnabletoMoveKernels(beta, srmc, LDProb, lrmc, local_weight, LD_weights_nocon)
  RMK <- LDRestrictedKernels(beta, srmc, ddWeight_sr, local_weight)
  
  ############## Rates for Diff Eqs  ##############  
  klist <- SetupRates(kD0,Rr0,1)
  Ks <- Define_Ks(klist, Vrate, Crate, 0, tstep)
  vaccinating <- FALSE
  vaccinatingH <- FALSE
  culling <- FALSE
  
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
      
        # print(paste("week", week))
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
          } # Check if number of sick and dead is greater than input threshold for detection
        
          ############## Cull Vacc set-up  ##############
          current_time <- day+weektime+yeartime
          if ((Vdelay_corr < current_time) && (Vrate>0)) vaccinating <- TRUE
          if ((VdelayH_corr < current_time) && (VrateH>0)) {
            vaccinatingH <- TRUE
          }
          if ((Cdelay_corr < current_time) && (Crate>0)) culling <- TRUE
          
          
          ChooseCells <- function(radius, N, Cutoff, PixParams)  {
            choosencellsG1 <- NULL
            choosencellsG2 <- NULL
            choosencellsG3 <- NULL
            choosencellsG4 <- NULL
            choosencellsI <- which(N[3,1,1,] >= 1) #2*Cutoff)
            choosencellsH <- which(N[4,1,1,] >= 1) #2*Cutoff)
            choosencellsG1 <- union(choosencellsI,choosencellsH)        # Group 1 are boxes that contain a sick animal 
            if (!(length(choosencellsG1)>0)) {  # if nothing was chosen so far
              if (max(N[3,1,1,]) > 1e-4) {
                choosencellsG1 <- which(N[3,1,1,]==max(N[3,1,1,]))
              }
            }
            if (radius > PixParams$pixel_height) {                      # Group 2 are boxes adjacent to those conatining a sick animal 
              choosencells <- unique(unlist(lapply(choosencellsG1, function(x) which(dd[x,] <= (PixParams$pixel_height+0.1)))))
              choosencellsG2 <- setdiff(choosencells,choosencellsG1)
            } else {
              if (radius > PixParams$pixel_width) {
                choosencells <- unique(unlist(lapply(choosencellsG1, function(x) which(dd[x,] <= (PixParams$pixel_width+0.1)))))
                choosencellsG2 <- setdiff(choosencells,choosencellsG1)
              } 
            }  
            KittyCorneredDist <- sqrt(PixParams$pixel_height^2+PixParams$pixel_width^2)
            if (radius > KittyCorneredDist) {
              choosencells <- unique(unlist(lapply(choosencellsG1, function(x) which(dd[x,] <= (KittyCorneredDist+0.1)))))
              choosencellsG3 <- setdiff(choosencells,union(choosencellsG1,choosencellsG2))  # Group 3 are boxes KittyCornered to those conatining a sick animal 
              choosencells <- unique(unlist(lapply(choosencellsG1, function(x) which(dd[x,] <= radius))))   # Find all cells within radius of a cell with sick cows.
              choosencellsG4 <- setdiff(choosencells,union(union(choosencellsG1,choosencellsG2),choosencellsG3))  # Group 4 is the rest
            } 
            CellsChoosen <- list(G1=choosencellsG1,G2=choosencellsG2,G3=choosencellsG3,G4=choosencellsG4)
            return(CellsChoosen)
          }
          
          CellsToCull <- NULL
          if (culling) {
            CellsToCull <- ChooseCells(Cullrad, N, animalCut/26*ApproxPixelArea, PixParams)
            if (!(length(CellsToCull$G1)>0)) CellsToCull <- NULL
          }  
        
          CellsToVacc <- NULL
          if (vaccinating) {
              CellsToVacc <- ChooseCells(vacrad, N, animalCut/26*ApproxPixelArea, PixParams)
              if (!(length(CellsToVacc$G1)>0)) CellsToVacc <- NULL
          }
        
          ############## spread matrix  ##############
          SpreadingMatricesForDay <- function(current_time, lrmc_stime, srmc_stime, AMK, UMK, RMK, cowTypeMov) {
            # Choose spreading matrix based on what movement controls are in place at this time
            if ((lrmc_stime <= current_time) && (srmc_stime <= current_time)) {  # motion and transmission controls in place
              weights_w_LDmov <- AMK$ddWeight_all_con
              weights_wo_LDmov <- RMK$weights_con
              weights_no_iso_w_LDmov <- UMK$ddWeight_local_all_con
              weights_no_iso_wo_LDmov <- RMK$weights_local_con
            } else {
              if (lrmc_stime <= current_time) {
                weights_w_LDmov <- AMK$ddWeight_LDcon
                weights_wo_LDmov <- RMK$weights
                weights_no_iso_w_LDmov <- UMK$ddWeight_local_LDcon
                weights_no_iso_wo_LDmov <- RMK$weights_local
              }
              if (srmc_stime <= current_time) {
                weights_w_LDmov <- AMK$ddWeight_Tcon
                weights_wo_LDmov <- RMK$weights_con
                weights_no_iso_w_LDmov <- UMK$ddWeight_local_Tcon
                weights_no_iso_wo_LDmov <- RMK$weights_local_con
              }
              if ((lrmc_stime > current_time) && (srmc_stime > current_time)) {
                weights_w_LDmov <- AMK$ddWeight_nocon 
                weights_wo_LDmov <- RMK$weights
                weights_no_iso_w_LDmov <- UMK$ddWeight_local_nocon  # dd.weight_local_nocon <- beta*((1-LDProb/100)*local_weight + LD_weights_nocon),    LD_weights_nocon[LDcell,] <- LDProb/100
                weights_no_iso_wo_LDmov <- RMK$weights_local        # weights_local <- beta*local_weight
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
            MM <- list(dd.weight_totI=dd.weight_totI, dd.weight_totHplus=dd.weight_totHplus)
            return(MM)
          }
          
          MM <- SpreadingMatricesForDay(current_time, lrmc_stime, srmc_stime, AMK, UMK, RMK, cowTypeMov)
            
          ############## tstep loop  ##############
          for (n in 1:(1/tstep)) {  #  1/tstep must be an integer
          
            Nold <- N
            results <- NextN(N, ApproxPixelArea, MM, tstep, klist, Ks, CellsToVacc, CellsToCull, CullFrom, DorS) 
            N <- results$N
            incidence[1,,day] <- incidence[1,,day] + results$incidentbybox
            
          }  #End of tstep loop
        
          record[,,,,day] <- N
        
        }  #end day loop 
        wk.prev[,,,,week]=rowMeans(record,dims=4)  #Mean of all days
        wk.inc[,,week]=rowSums(incidence,dims=2)                          # Summing over days in the week 
      }  #end week loop
    })  # end withProgress
  }# end year 
  out=list(wk.prev=wk.prev,wk.inc=wk.inc, usedvacc=usedvacc, ActualDetectTime=ActualDetectTime, ReallySickorDead=ReallySickorDead, Cdelay_corr=Cdelay_corr, DorS=DorS, UsedParams=UsedParams)  
  return(out)
}

########################################################################### runSim ###################################################
runSim=function(N,r,dl,dd,icell,TimeParams,beta,betaRatio,Rr0,kD0,DetectSick,srmcDelay,srmcBetaChange,VacParams,CullParams, LDMoveParams, DorS, Nave, UsedParams) {
  cat("In runSim", "(functions)", "\n")
  # Ndstates=dim(N)[1]
  # Ndcomb=dim(N)[2]
  # Nages=dim(N)[3]
  nbox=dim(N)[4]
  
  if (LDMoveParams$LDType == "roads") {  # Pakistan simulation example
    N[2,1,,icell] <- 1
  } else {
    # Start with a few cases very roughly in the correct porportion of E, I and H
    gamma <- kD0$kD0dead[2]  # kEI
    kappa <- kD0$kD0recov[3] + kD0$kD0dead[3]  # kIH + kIR
    #cases = H + beta/gamma*H + 1/2*beta/gamma*beta/gamma*H + 1/2*beta/kappa*H
    H <- round(10/(1 + beta/gamma + 0.5*(beta/gamma)^2 + 0.5*beta/kappa))
    if (H < 1) H <- 1
    if (H > 2) H <- 2
    
    N[2,1,,icell] <- round(0.5*((beta/gamma)^2 + beta/kappa)*H)  # index case(s) defined.  1/2*beta/gamma*I + 1/2*beta/kappa*H; gamma = kD0$kD0dead[2]  kappa = kD0$kD0recov[3] + kD0$kD0dead[3]
    N[3,1,,icell] <- round(beta/gamma*H)  # index case(s) defined.  beta/gamma*H
    N[4,1,,icell] <- H  # index case(s) defined.  
    if (sum(N[2:4,1,1,icell]) < 7) N[2:4,1,1,icell] <- N[2:4,1,1,icell]*ceiling(6/sum(N[2:4,1,1,icell]))
    if (sum(N[2:4,1,1,icell]) > 10)  {
      aa <- sum(N[2:4,1,1,icell])/7
      N[2,1,1,icell] <- ceiling(N[2,1,1,icell]/aa)
      N[3,1,1,icell] <- ceiling(N[3,1,1,icell]/aa)
      N[4,1,1,icell] <- ceiling(N[4,1,1,icell]/aa)
    }
  }
  
  if (DorS == "Stochastic") {
    RunsToAve <- Nave
  } else {
    RunsToAve <- 1
  }
  
  sim<-episim(N,r,dl,dd,TimeParams,beta,betaRatio,Rr0,kD0,DetectSick,srmcDelay,srmcBetaChange,VacParams,CullParams, icell, 
              LDMoveParams, DorS, UsedParams)  
  wkPrevOut <- sim$wk.prev
  wkIncOut <- sim$wk.inc
  
  if (RunsToAve > 1) {
    for (i in 1:(RunsToAve-1)) {
      sim <-episim(N,r,dl,dd,TimeParams,beta,betaRatio,Rr0,kD0,DetectSick,srmcDelay,srmcBetaChange,VacParams,CullParams, icell,
                LDMoveParams, DorS, UsedParams)
      wkPrevOut <- wkPrevOut + sim$wk.prev
      wkIncOut <- wkIncOut + sim$wk.inc
    }
    wkPrevOut <- wkPrevOut/RunsToAve
    wkIncOut <- wkIncOut/RunsToAve
  }
  
  out=list(wk.prev=wkPrevOut,wk.inc=wkIncOut, usedvacc=sim$usedvacc, ActualDetectTime=sim$ActualDetectTime, ReallySickorDead=sim$ReallySickorDead, Cdelay_corr=sim$Cdelay_corr, DorS=sim$DorS, UsedParams=sim$UsedParams)  
  return(out)
}

diseasekD <- function(Disease) {
  kD2dead  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  kD2recov  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  if (Disease == "Rinderpest") {
    kD1dead   <- c(0, kEI_Rind, kIH_Rind, kHD_Rind, kItHt_Rind, kHtD_Rind, 0, 0, 0)  
    kD1recov  <- c(0, 0,        kIR_Rind, kHR_Rind, kItR_Rind,  kHtR_Rind, 1/vac_TimeToImmunity_Rind, 0, 0)  
  }
  return(cbind(kD1dead,kD1recov,kD2dead,kD2recov))
}

diseaseRr <- function(Disease) {
  # Rate Ratio           kHHt/(kHD+kHR)              
  Rr1treat <- c(0, 0, 0, Rr_HHt_Rind,       0, 0, 0, 0, 0)  # Changed Ht-R from 10 to 100 July 2016
  Rr2treat <- NULL
  return(cbind(Rr1treat,Rr2treat))
}
 