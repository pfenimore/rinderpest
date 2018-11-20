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
      #distances <- diag(new_mat[2:dim(testarray)[1],1:(dim(testarray)[1]-1)])  #just off diagonal, so dist between consequtive points - CHANGE to distm
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
    if (Neg_Intersect_Pt == 1) {
      #print("at the start of the road.")
    } else {
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
    #cat("ContagionCells NEGATIVE before intersection \t", ContagionCells, "\n")
    cellsC <- c(cellsC, ContagionCells)
  }
  
  if (IntersectPos){
    #print("There is an intersection in POSITIVE direction")
    DistanceMoved <- Pos_intersect_dist-StartDist
    if (Pos_Intersect_Pt == length(road_list[[RoadIndex]]$distances) ) {
      #print("at the start of the road.")
    } else {
      #print("in the middle of the road.")
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
    #cat("ContagionCells POSITIVE \t", ContagionCells, "\n")
    cellsC <- c(cellsC, ContagionCells)
  }
  
  if ( (DirectionType < 1) && (!IntersectNeg))  {
    if (StartDist-max_dist < 0) {
      #print("Going off start of road")
    } 
    #print("Checking negative direction: no intersection, no end.")
    points <- which((road_list[[RoadIndex]]$distances > StartDist-max_dist) & (road_list[[RoadIndex]]$distances < StartDist-min_dist))
    ContagionCells <- unique(road_list[[RoadIndex]]$CellsOnRoad[points])
    #cat("ContagionCells NEGATIVE \t", ContagionCells, "\n")
    cellsC <- c(cellsC, ContagionCells)
  }
  
  if ( (DirectionType > -1) && (!IntersectPos) ) { 
    if (StartDist+max_dist > max(road_list[[RoadIndex]]$distances)) {
      #print("Going off end of road")  
    }  
    #print("Checking positive direction: no intersection, no end.")
    points <- which( (road_list[[RoadIndex]]$distances > (StartDist+min_dist)) & (road_list[[RoadIndex]]$distances < (StartDist+max_dist)) )
    ContagionCells <- unique(road_list[[RoadIndex]]$CellsOnRoad[points])
    #cat("ContagionCells POSITIVE \t", ContagionCells, "\n")
    cellsC <- c(cellsC, ContagionCells)
  }  
  return(list(RoadSegs=RoadSegsToDo, cellsC=cellsC))
}

road_transport <- function(choosencells, RL, KeptRoads, road_list, RoadSpread, IntersectInfo, regioncattle) {
  print("Starting road transport")
  #print(IntersectInfo$ILines)
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
