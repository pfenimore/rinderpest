#EpiGrid server.r

rm(list = ls(all.names = TRUE))  #Clears workspace

 source(paste(normalizePath(getwd()),"functions.r", sep="/"))
 source(paste(normalizePath(getwd()),"PakistanRoadFunctions.r", sep="/"))
source(paste(normalizePath(getwd()),"FunctionsSame.r", sep="/"))
# source(paste(normalizePath(getwd()),"functions Deriggi works.r", sep="/"))

if (!(dir.exists("www"))) {
  dir.create("www")
}
wwwDir <- "www"

Ndstates<- 9        # Healthy and Disease states S,E,I,H,Ht,D,R
Ndcomb  <- 1        # Number of diseases   
Nages   <- 1        # One age category.

# Times for epidemic simulation
tstep   <- 0.2     # units are days  1/tstep must be integer!
nyear   <- 1        # year
len.wk  <- 7        # days in a week

# Times for disease progression calculation
nmonth  <- 2        # For disease progression calculations
len_mo  <- 30       # Number of days in a month
nstep   <- nmonth*len_mo/tstep

# Starting values (the reactive expressions only keep things up to date, they don't do initial calculations)
wb  <- extent(longmin,longmax,latmin,latmax)  # See global.r
rr  <- initGeog(wb, 1, "cattle")              # returns a raster data of population
Ns  <- values(rr)  
Ns[is.na(Ns)] <- 0
values(rr) <- Ns
if (AggFact_BasicStart != 1) {
  r   <- Reduce_res(rr, AggFact_BasicStart)
} else {
  r <- rr
}
rV  <- setValues(r, Vbase_BasicStart)
N   <- initN(Ndstates,Ndcomb,Nages,r,rV)

#SimType = "regular"  

IndexCell <- cellFromXY(r,c(longval,latval))
Lcellnum  <- 1            # Long distance cell

MinRoadlength <- 13 # km  relavent to Pakistan example only.  As of March 1, 2016, making this bigger gets rid of some lines we want.
LDType <- "P2P"     # Type of long distance movement - point-to-point.

#                 E-I       I->        H->      It->       Ht ->     Vg->Vm
kD1dead  <- c(0, kEI_Rind, kIH_Rind, kHD_Rind, kItHt_Rind, kHtD_Rind, 0,                         0,  0)  
kD1recov <- c(0, 0,        kIR_Rind, kHR_Rind, kItR_Rind,  kHtR_Rind, 1/vac_TimeToImmunity_Rind, 0,  0)  

# Rate Ratio       I->It   H->Ht                      # to treatment
Rr1treat <- c(0, 0, 0,  Rr_HHt_Rind,  0, 0, 0, 0, 0)
Rr2treat <- c(0, 0, 0,  Rr_HHt_Rind, 0, 0, 0, 0, 0)

kD0dead <- array(dim=c(Ndstates, 1, Nages))   
kD0dead[,1,] <- kD1dead

kD0recov <- array(dim=c(Ndstates, 1, Nages))   
kD0recov[,1,] <- kD1recov

Rr0 <- array(dim=c(Ndstates, 1, Nages))  # rate ratio for treatment
Rr0[,1,] <- Rr1treat

###############  Begin ShinyServer #################
shinyServer(function(input, output, session) {   
  
  ##### Get ModelType and disease parameters ######
  observe({  # code is run when input$ModelType or input$Disease changes.  ModelType options are "No" (default disease parameters), "Use Case" (examples) and "load parameter"
      if (input$ModelType == "No") {  # Default Disease Parameters
        updateRadioButtons(session, "epiOutMapButton", selected = "orig")
        # Turn mitigations off
        updateNumericInput(session, "Vrate", value=0)
        updateNumericInput(session, "Crate", value=0)
        updateSliderInput(session, "srmc", value=1) 
        
        if (input$Disease == "Rinderpest") {
          updateSliderInput(session, "beta", value=beta_Rind)
          updateSliderInput(session, "IncTime", value = 1/kEI_Rind)
          updateNumericInput(session, "kIH", value = kIH_Rind)
          updateNumericInput(session, "kIR", value = kIR_Rind)
          updateNumericInput(session, "kHD", value = kHD_Rind)
          updateNumericInput(session, "kHR", value = kHR_Rind)
          #updateSliderInput(session, "mortRate", value = kHD_Rind/(kHD_Rind+kHR_Rind))
          updateSliderInput(session, "treated", value = 0)
          updateSliderInput(session, "Vbase", value=0)
          updateSliderInput(session, "vacprog", value=vac_TimeToImmunity_Rind)
        }
        if (input$Disease == "Rinderpest lineage 1") {
          updateSliderInput(session, "beta", value=1.2)
          updateSliderInput(session, "IncTime", value = 5.6)
          updateNumericInput(session, "kIH", value = 0.33)
          updateNumericInput(session, "kIR", value = 0)
          updateNumericInput(session, "kHD", value = 0.12)
          updateNumericInput(session, "kHR", value = 0.2133)
          #updateSliderInput(session, "mortRate", value = kHD_Rind/(kHD_Rind+kHR_Rind))
          updateSliderInput(session, "treated", value = 0)
          updateSliderInput(session, "Vbase", value=0)
          updateSliderInput(session, "vacprog", value=vac_TimeToImmunity_Rind)
        }
        if (input$Disease == "Rinderpest lineage 2") {
          updateSliderInput(session, "beta", value=0.18)
          updateSliderInput(session, "IncTime", value = 6.8)
          updateNumericInput(session, "kIH", value = 0.14)
          updateNumericInput(session, "kIR", value = 0)
          updateNumericInput(session, "kHD", value = 0.02)
          updateNumericInput(session, "kHR", value = 0.3133)
          #updateSliderInput(session, "mortRate", value = kHD_Rind/(kHD_Rind+kHR_Rind))
          updateSliderInput(session, "treated", value = 0)
          updateSliderInput(session, "Vbase", value=0)
          updateSliderInput(session, "vacprog", value=vac_TimeToImmunity_Rind)
        }
        if (input$Disease == "Rinderpest Pak94") {
          updateSliderInput(session, "beta", value=betaPak94)
          updateSliderInput(session, "IncTime", value = 1/kEI_Pak94)
          updateNumericInput(session, "kIH", value = kIH_Pak94)
          updateNumericInput(session, "kIR", value = kIR_Pak94)
          updateNumericInput(session, "kHD", value = kHD_Pak94)
          updateNumericInput(session, "kHR", value = kHR_Pak94)
          #updateSliderInput(session, "mortRate", value = kHD_Rind/(kHD_Rind+kHR_Rind))
          updateSliderInput(session, "treated", value = 0)
          updateSliderInput(session, "Vbase", value=0)
          updateSliderInput(session, "vacprog", value=vac_TimeToImmunity_Rind)
        }
      }
      if (input$ModelType == "Use Case") {
        updateRadioButtons(session, "epiOutMapButton", selected = "orig")
        if (input$Example  == "Rinderpest: Pakistan 1994") {
          updateSelectInput(session, "Disease", selected = "Rinderpest Pak94", choices = "Rinderpest Pak94") 
          updateSelectInput(session, "ModelSize", selected = "Large Rect: 3 x 4 degrees", choices = "Large Rect: 3 x 4 degrees")
          updateNumericInput(session, "loniC", value=PakLoniC) #Pari wikiglobe 73.96  Google maps shows it on boths sides of highway
          updateNumericInput(session, "latC", value=PakLatC) #35.81  #34.99 wikiglobe Pari
          updateNumericInput(session, "weeks_to_sim", value=41)  #42
          updateSliderInput(session, "DetectSick", value=500)  # In truth it was not detected until August when thousands of cattle had died.
          updateSelectInput(session, "AggFact", selected=1.5)
        
          updateSliderInput(session, "dl", value=1.5)  
          updateSliderInput(session, "srmc", value=0.8) 
          updateSliderInput(session, "srmcd", value=40) # 30 
          updateSliderInput(session, "LDProb", value=0.5) 
          updateSliderInput(session, "RoadDist", value=35) 
          updateSliderInput(session, "lrmcd", value=40)    # 30
          updateSliderInput(session, "lrmc", value=0.2)     
          updateRadioButtons(session, "cowTypeMov", selected = "NoH")
          
          updateNumericInput(session, "Crate", value=0)
          updateSliderInput(session, "Culldelay", value=0)
          
          updateNumericInput(session, "Vrate", value=200)
          updateSliderInput(session, "Vdelay", value=70)  # 60 Vaccination started in August  So this plus disease detection time should be ~150.
          updateSliderInput(session, "vacrad", value=30)  
          updateSliderInput(session, "Vbase", value=0)
          updateSliderInput(session, "vacprog", value=vac_TimeToImmunity_Rind)
          #updateSliderInput(session, "Veff", value=0.6)   
          
          updateSliderInput(session, "beta", value=betaPak94) # 0.35 # 0.3 - An estimate based on comments by Rossiter.
          updateSliderInput(session, "IncTime", value=1/kEI_Pak94) 
          updateNumericInput(session, "kIH", value = kIH_Pak94)
          updateNumericInput(session, "kIR", value = kIR_Pak94)
          updateNumericInput(session, "kHD", value = kHD_Pak94) #Mortality rates were usually over 80% and near 100% for Yakmo. 
          updateNumericInput(session, "kHR", value = kHR_Pak94)
          # updateSliderInput(session, "mortRate", value=0.9) #Mortality rates were usually over 80% and near 100% for Yakmo. 
          updateSliderInput(session, "treated", value = 0)
        }   #Pakistan rinderpest epidemic 1994
        if (input$Example  == "Rinderpest: Fremantle") {
          updateSelectInput(session, "Disease", selected = "Rinderpest", choices = "Rinderpest") 
          updateSelectInput(session, "ModelSize", selected = "Small: 0.8 x 0.8 degrees", choices = "Small: 0.8 x 0.8 degrees")
          updateNumericInput(session, "loniC", value=115.61) #
          updateNumericInput(session, "latC", value=-31.19) #
          updateNumericInput(session, "weeks_to_sim", value=10)  #42
          updateSliderInput(session, "DetectSick", value=5)  # 
          updateSelectInput(session, "AggFact", selected=1.0)
          
          updateSliderInput(session, "dl", value=1.3)  
          updateSliderInput(session, "srmc", value=0.5) 
          updateSliderInput(session, "srmcd", value=2) # 30 
          updateSliderInput(session, "LDProb", value=0) 
          
          updateNumericInput(session, "Crate", value=1500)
          updateSliderInput(session, "Culldelay", value=2)
          updateSliderInput(session, "CullRad", value=7)  
          updateSliderInput(session, "CdelayPInfect", value=1)
          
          updateNumericInput(session, "Vrate", value=0)

          updateSliderInput(session, "beta", value=1.0) # 
          updateSliderInput(session, "IncTime", value=1/kEI_Rind) 
          updateNumericInput(session, "kIH", value = kIH_Rind)
          updateNumericInput(session, "kIR", value = kIR_Rind)
          updateNumericInput(session, "kHD", value = 0.8/3) # 
          updateNumericInput(session, "kHR", value = 0.2/3)
          # updateSliderInput(session, "mortRate", value=0.9) #
          updateSliderInput(session, "treated", value = 0)
        }   #Fremantle
      }  
    }) #Run when input$ModelType or input$Disease changes.  Sets disease default parameters or example parameters.
  
  observe({  
    if(input$ModelType == "No") updateSelectInput(session, "Disease", choices = DiseaseList)
    if(input$ModelType == "No") updateSelectInput(session, "ModelSize", choices = SizeList)
  })  # Changes disease and modeltype sliders
  
  observeEvent(input$csvFile, {
      inFile <- input$inputs
      if (is.null(inFile))
        return(NULL)
      
      updateRadioButtons(session, "epiOutMapButton", selected = "orig")
      
      InputVals <-read.csv(inFile$datapath, stringsAsFactors=F)
      updateSelectInput(session, "Disease", choices = DiseaseList, selected = InputVals[1,2])
      updateSelectInput(session, "ModelSize", selected = InputVals[2,2], choices = SizeList)
      updateNumericInput(session, "latC", value=as.numeric(InputVals[3,2]))
      updateNumericInput(session, "loniC", value=as.numeric(InputVals[4,2])) 
      updateNumericInput(session, "weeks_to_sim", value=as.numeric(InputVals[5,2]))
      updateSliderInput(session, "DetectSick", value=as.numeric(InputVals[6,2]))
      updateSelectInput(session, "AggFact", choices = c(0.5, 1, 1.5, 2), selected=InputVals[7,2])
      
      # read in Disease parameters
      updateSliderInput(session, "beta", value=as.numeric(InputVals[8,2]))
      updateSliderInput(session, "IncTime", value=as.numeric(InputVals[9,2])) 
      updateNumericInput(session, "kIH", value = as.numeric(InputVals[10,2]))
      updateNumericInput(session, "kIR", value = as.numeric(InputVals[11,2]))
      updateNumericInput(session, "kHD", value = as.numeric(InputVals[12,2]))
      updateNumericInput(session, "kHR", value = as.numeric(InputVals[13,2]))
      kHHt <- as.numeric(InputVals[14,2])
      kHD <- as.numeric(InputVals[12,2])
      kHR <- as.numeric(InputVals[13,2])
      rateRatio <- kHHt/(kHD+kHR)
      treated = rateRatio/(rateRatio+1)
      updateSliderInput(session, "treated", value=treated)
      
      updateNumericInput(session, "kHtD", value = as.numeric(InputVals[15,2]))
      updateNumericInput(session, "kHtR", value = as.numeric(InputVals[16,2]))

      
      # read in "short" range movement and control parameters
      updateSliderInput(session, "dl", value=as.numeric(InputVals[17,2]))  
      updateSliderInput(session, "srmc", value=as.numeric(InputVals[18,2]))  
      updateSliderInput(session, "srmcd", value=as.numeric(InputVals[19,2])) 
      
      # read in "long" range movement and control parameters
      updateSliderInput(session, "LDProb", value=as.numeric(InputVals[20,2]))
      corners <- Rcorners()
      if (InputVals[20,2]>0) {
        updateSliderInput(session, "Lloni", min = corners$longmin+0.1, max = corners$longmax-0.1, value = as.numeric(InputVals[21,2]))
        updateSliderInput(session, "Llati", min = corners$latmin+0.1,  max = corners$latmax-0.1, value = as.numeric(InputVals[22,2]))
        updateRadioButtons(session, "cowTypeMov", selected = InputVals[23,2])
      }
      updateSliderInput(session, "lrmc", value=as.numeric(InputVals[24,2]))
      updateSliderInput(session, "lrmcd", value=as.numeric(InputVals[25,2]))
      
      # read in Culling parameters
      updateNumericInput(session, "Crate", value=as.numeric(InputVals[26,2]))
      if (InputVals[26,2]>0) {
        updateSliderInput(session, "Culldelay", value=as.numeric(InputVals[27,2]))
        updateRadioButtons(session, "CullFrom", selected = InputVals[28,2])
      }
      
      # read in Vaccination parameters
      updateSliderInput(session, "Vbase", value=as.numeric(InputVals[29,2]))
      updateNumericInput(session, "Vrate", value=as.numeric(InputVals[30,2]))
      if (InputVals[30,2] > 0) {
        updateSliderInput(session, "Vdelay", value=as.numeric(InputVals[31,2]))
        updateSliderInput(session, "vacrad", value=as.numeric(InputVals[32,2]))
        updateSliderInput(session, "vacprog", value=as.numeric(InputVals[33,2]))
      }
      
  }) # Loads parameters from a csv file.
  
  reactLDType = reactive({  
    if ((input$ModelType == "Use Case") && (input$Example  == "Rinderpest: Pakistan 1994")) {
      LDType <-"roads"
    } else{
      LDType <- "P2P"
    }
  }) # Sets type of long distance movement, usually point-to-point
  
  ##### Update Leaflet map ######  
  observe({  #Changes Population map when ModelSize or Disease changes, also changes long distance sliders
    corners <- Rcorners()
    if (is.null(corners)) return(NULL)
    SquareColor = "red"
    leafletProxy("mymap") %>%
      setView("mymap", lng=input$loniC, lat=input$latC, zoom=3) %>%
      removeShape(layerId = "R") %>%
      addRectangles(
        lng1=corners$longmin, lat1=corners$latmin,
        lng2=corners$longmax, lat2=corners$latmax,
        fill=FALSE, weight = 1.5, color = SquareColor, layerId = "R") %>%
      removeShape(layerId = "I") %>%
      addMarkers(lng=input$loniC, lat=input$latC, layerId = "I")
  }) # Alters Population map and long distance slider limits when either ModelType, index location, or Disease type changes.
  
  observeEvent(input$mymap_click, {
    click <- input$mymap_click
    clat <- click$lat
    clng <- click$lng
    
    updateNumericInput(session, "loniC", value=clng) #Pari wikiglobe 73.96  Google maps shows it on boths sides of highway
    updateNumericInput(session, "latC", value=clat) #35.81  #34.99 wikiglobe Pari
  }) # Move map marker and set new center coordinates
  
  observe({  #Runs when input$weeks_to_sim changes.  Updates week slider for incidence etc.
    val <- input$weeks_to_sim
    updateSliderInput(session, "P.wk", min=1.0, max=val, step=1.0)
  }) # Runs when input$weeks_to_sim changes.  Updates week slider for output maps.
  
  ##### Geography #####
  Rcorners = reactive({
    shiny::validate(
       need(input$latC < 85, "Latitude must be < 85"),
       need(input$latC > -60, "Latitude must be > -60"),
       need(input$loniC < 180, "Longitude must be < 180"),
       need(input$loniC > -180, "Longitude must be > -180")
    )
    if ((input$ModelType == "Use Case") &&  (input$Example  == "Rinderpest: Pakistan 1994")) {
      shiny::validate(
        need(input$latC <= 35.87, "Latitude must be less than or equal to 35.87"),
        need(input$latC >= 35.6, "Latitude must be greater than or equal to 35.60"),
        need(input$loniC <= 74.57, "Longitude must be less than or equal to 74.57"),
        need(input$loniC >= 74.40, "Longitude must be greater than or equal to 74.40")
      )
    }
    
    RectLengths <- LatLongLengths(input$ModelSize) # In FunctionsSame.r
    lat <- RectLengths$lat
    long <- RectLengths$long
    
    longmin <- input$loniC-long/2
    longmax <- input$loniC+long/2
    latmin <- input$latC-lat/2
    latmax <- input$latC+lat/2
    corners <- list(longmin=longmin, longmax=longmax, latmin=latmin, latmax=latmax)
    return(corners)
  }) # (reactive: ModelSize, ModelType, Disease, lat, long)
  
  Pop2 = reactive({  #Cuts out the part of the populaton map that is needed, returns rr
    corners <- Rcorners()
    if (is.null(corners)) return(NULL)
    wb   <- extent(corners$longmin,corners$longmax,corners$latmin,corners$latmax)
    HistCatRed <- 1
    if (input$ModelType == "Use Case") {  # Don't want && because Example might not be defined
      if (input$Example  == "Rinderpest: Pakistan 1994") {
        HistCatRed <- PakHistCatRed             
      }
    }
    rr <- initGeog(wb, HistCatRed, "cattle")
    return(rr)
  })  # Cuts out the part of the populaton map that is needed, returns rr
  
  rr2r=reactive({ #Does resolution reduction
    rr <- Pop2()
    if (is.null(rr)) return(NULL)
    Ns <- values(rr)  
    Ns[is.na(Ns)] <- 0
    values(rr)    <- Ns
    
    if (as.double(input$AggFact) != 1) {r=Reduce_res(rr, as.double(input$AggFact))} else {r=rr}  
    return(r)
  }) #Does resolution reduction
  
  PixelParams = reactive ({
    r <- rr2r()
    if (is.null(r)) return(NULL)
    middlecoords <- GetMiddleCoords(r) # Units are lat long
    p2 <- c(middlecoords[1] + res(r)[1], middlecoords[2])
    pixel_width <- distMeeus(middlecoords, p2, a=6378137, f=1/298.257223563)/1000 #Units are km
    p3 <- c(middlecoords[1], middlecoords[2]+res(r)[2])
    pixel_height <- distMeeus(middlecoords, p3, a=6378137, f=1/298.257223563)/1000 #Units are km 
    PixParams <- list(pixel_width=pixel_width, pixel_height=pixel_height)
    return(PixParams)
  })
  
  changeDD = reactive ({
    r <- rr2r()
    # Tried other distance functions,  distGeo takes 1.8 times longer than distMeeus
    # distVincentyEllipsoid is really slow. 
    # distCosine, distHaversine, distVincentySphere same as distMeeus
    withProgress(message='Calculating', {
      DD <- distm(coordinates(r),fun=distMeeus)/1000 #Default units are m so change to km
      incProgress(1)
    })
    return(DD)
  }) 
  
  getRoads = reactive({
    corners <- Rcorners()
    if (is.null(corners)) return(NULL)
    wb   <- extent(corners$longmin,corners$longmax,corners$latmin,corners$latmax)
    roads <- initRoads(wb)   
  }) # Pakistan roads: Loads in road data, crops, and merges (connects) lines.
  
  roadIndices = reactive({
    roads <- getRoads() 
    RoadIndices <- LongRoads(roads, MinRoadlength)
    return(RoadIndices)
  })
  
  FixRoads = reactive({
    roads <- getRoads()                                   # The Pakistan roads in given region
    index_KeptRoads <- roadIndices()                      # index_KeptRoads <- LongRoads(roads, MinRoadlength)
    RLstart <- roads@lines[[1]]@Lines                     # All of the Road Lines
    Fixes <- FixLines(RLstart, index_KeptRoads) 
    return(Fixes)
  })
  
  AboutRoads = reactive({
    r <- rr2r()
    Fixes <- FixRoads()
    RL <- Fixes$RL
    index_KeptRoads <- Fixes$index_KeptRoads
    RoadI <- roadInfo(RL, r, index_KeptRoads)
    return(RoadI)
  })
  
  Intersects = reactive({
    Fixes <- FixRoads()
    RL <- Fixes$RL
    index_KeptRoads <- Fixes$index_KeptRoads
    IntersectInfo <- Intersections(RL, index_KeptRoads) 
  })  
  
  ##### Base vaccination levels #####
  rVacc = reactive({
    r <- rr2r()
    rV <- r
    Vbase <- input$Vbase
    rV <- setValues(rV, Vbase)
    return(rV)
  })
  
  ##### Initialize N, population #####
  Agg=reactive({ #Returns N, the matrix of populations
    r <- rr2r()
    rV <- rVacc()
    N   <- initN(Ndstates,Ndcomb,Nages,r,rV)
    return(N)
  }) # Returns N, the matrix of populations
  
  ##### Index Cell  #####  
  icell=reactive({
    r=rr2r()
    IndexCell=cellFromXY(r,c(input$loniC,input$latC))
  })  # Gets index cell number
  
  ##### Long distance point #####
  ldcell=reactive({     #Sets long distance cell
    r=rr2r()
    corners <- Rcorners()
    if (is.null(corners)) return(NULL)
    # if (input$ModelType=="Load parameters") {
    #   LDLoc <- LDFileLocation()  # This will wait until file is input.
    #   if (is.null(LDLoc))
    #     return(NULL)
    #   Lcellnum=cellFromXY(r,c(LDLoc$LDlong,LDLoc$LDLat))
    # } else {
      if (is.null(input$Lloni)) { 
        Lcellnum=cellFromXY(r,c(corners$longmax-0.2,corners$latmax-0.2))
        Lcellnum2 <- cellFromXY(r,c(corners$longmin+0.2,corners$latmax-0.2))
        Lcellnum3 <- cellFromXY(r,c(corners$longmin+0.2,corners$latmin+0.2))
        Lcellnum4 <- cellFromXY(r,c(corners$longmax-0.2,corners$latmin+0.2))
      }else{ 
        Lcellnum=cellFromXY(r,c(input$Lloni,input$Llati))
        symmetricLat <- input$latC - (input$Llati-input$latC)
        symmetricLong <- input$loniC - (input$Lloni-input$loniC)
        Lcellnum2 <- cellFromXY(r,c(input$Lloni,symmetricLat))
        Lcellnum3 <- cellFromXY(r,c(symmetricLong,symmetricLat))
        Lcellnum4 <- cellFromXY(r,c(symmetricLong,input$Llati))
        #a = b - (c-b) = 2b-c
        #c = b + (b-a)  =2b-a      
      } 
    #}
    Lcellnums <- c(Lcellnum, Lcellnum2, Lcellnum3, Lcellnum4)
    return(Lcellnums)
  }) # Gets long distance cell number
  
  LDFileLocation=eventReactive(input$csvFile, {
    inFile <- input$inputs
    if (is.null(inFile))
      return(NULL)
    InputVals <-read.csv(inFile$datapath, stringsAsFactors=F)
    LDlong <- as.numeric(InputVals[12,2])
    LDlat <- as.numeric(InputVals[13,2])
    LDloc <- list(LDlong=LDlong,LDLat=LDlat)
    return(LDloc)
  }) # Gets long distance point coordinates when csv file is being loaded
  
  LDLocation=reactive({  
    r=rr2r()
    corners <- Rcorners()
    if (is.null(corners)) return(NULL)
    if (input$ModelType=="Load parameters") {
      LDLoc <- LDFileLocation()  # This will wait until file is input.
      if (is.null(LDLoc))  return(NULL)
      LDLongVal <- LDLoc$LDlong
      LDLatVal <- LDLoc$LDLat
    } else {
      LDLongVal <- corners$longmax-0.2
      LDLatVal <- corners$latmax-0.2
    }
    LDPos <- list(LDLongVal=LDLongVal,LDLatVal=LDLatVal)
    return(LDPos)
  }) # Long distance point coordinates for updating sliders below
  
  observe({
    corners <- Rcorners()
    LDPos <- LDLocation()
    updateSliderInput(session, "Lloni", min = round(corners$longmin,1)+0.1, max = round(corners$longmax,1)-0.1, value = LDPos$LDLongVal)
    updateSliderInput(session, "Llati", min = round(corners$latmin,1)+0.1,  max = round(corners$latmax,1)-0.1, value = LDPos$LDLatVal)
  })  # Updates long distance sliders

  ##### Disease parameters ####
  kD_input=reactive({
    kD <- diseasekD(input$Disease)            # Get base parameters for disease; diseasekD is in functions.r
    kD0dead[,1,]  <- kD[,1]
    kD0recov[,1,] <- kD[,2]                   # Disease 1
    
    # Change rates according to input incubation period
    kD0dead[2,1,] <- 1/input$IncTime
  
    kD0recov[3,1,] <- input$kIR 
    kD0dead[3,1,]  <- input$kIH
    
    # Change rates according to mortality rate input
    kD0recov[4,1,] <- input$kHR 
    kD0dead[4,1,]  <- input$kHD
    
    # kHR1 <- kD0recov[4,1,]
    # kHD1 <- kD0dead[4,1,]
    # kA_H1 <- kHR1+kHD1          
    # newkHD1 <- input$mortRate*kA_H1      
    # newkHR1 <- kA_H1 - newkHD1
    # kD0recov[4,1,] <- newkHR1
    # kD0dead[4,1,]  <- newkHD1
    
    # Change time for vaccine to be effective
    kD0recov[7,1,] <- 1/input$vacprog
    
    kD0 <- list(kD0dead=kD0dead, kD0recov=kD0recov)
    return(kD0)
  })  # updates kD when disease inputs change
 
  Rr_input=reactive({
    Rr <- diseaseRr(input$Disease)
    Rr0[,1,] <- Rr[,1]
    Rr_HHt1 <- input$treated/(1-input$treated)
    Rr0[4,1,] <- Rr_HHt1
    return(Rr0)
  })  # Rate ratios for treatment
  
  MakeRasterToPlot=reactive({   #<- function(kind, Rout, Disease, HumanOrAnimal, NByWeek, NByWeekH, IncByWeek, IncByWeekH) {
    sim <- sim1()
    PlotKind <- input$epiOutMapButton  # kind of plot
    Disease <- input$Disease  #sim$DiseaseModeled
    NByWeek <- sim$wk.prev
    IncByWeek <- sim$wk.inc
    week <- input$P.wk # week to generate results for
    Rout  <- rr2r()
    nbox <- length(Rout)
    if (input$epiOutMapButton != "orig")  {
      Rout <- RasterToPlotFunc(PlotKind, Rout, Disease, "Animal", NByWeek, NULL, IncByWeek, NULL, week, Nages, nbox)
    }
    return(Rout)
  })
  
  ##### Text outputs; pixel size, day detected, etc. #####
  output$pixelsize=renderText({  
    rr <-  Pop2()           # function generates geographical population
    if (is.null(rr)) return(NULL)
    r <- rr2r()             # change resolution geographical population,  
    PixParams <- PixelParams()
    paste("Average pixel width is", round(PixParams$pixel_width,1), "km")
  })  
  
  output$DetectionDay=renderText({  
    sim=sim1()
    if (sim$UsedParams$DiseaseModeled != input$Disease) return(NULL)
    if (sim$UsedParams$ModelType != input$ModelType) return(NULL)
    Day <- sim$ActualDetectTime
    if (!is.null(Day)) {
      text <- paste("Disease detected on day ", Day)
    } else {
      text <- "The disease was not detected"
    }
  })

  output$DetectionSick=renderText({  
    sim=sim1()
    if (sim$UsedParams$DiseaseModeled != input$Disease) return(NULL)
    if (sim$UsedParams$ModelType != input$ModelType) return(NULL)
    ReallySickorDead <- sim$ReallySickorDead
    Day <- sim$ActualDetectTime
    if (!is.null(ReallySickorDead)) {
      text <- paste("Very sick or dead at", Day, "days:", round(ReallySickorDead, digits=0))
    } else {
      text <- " "
    }
  })
  
  output$CullDay=renderText({  
    sim=sim1()
    if (sim$UsedParams$DiseaseModeled != input$Disease) return(NULL)
    if (sim$UsedParams$ModelType != input$ModelType) return(NULL)
    CDay <- sim$Cdelay_corr
    if (input$Crate > 0) {
      text2 <- paste("Culling began on day ", CDay)
    } else {
      text <- ""
    }
  })
  
  output$BoundaryReached=renderText({  
    sim=sim1()
    if (sim$UsedParams$DiseaseModeled != input$Disease) return(NULL)
    if (sim$UsedParams$ModelType != input$ModelType) return(NULL)
    if (sim$UsedParams$ModelSize != input$ModelSize) return(NULL)
    if (sim$UsedParams$AggFact != input$AggFact) return(NULL)
    r <- rr2r()
    nbox <- length(r)
    nrows <- dim(r)[1]
    ncols <- dim(r)[2]
    SimulationLength <- dim(sim$wk.inc)[3]  # weeks
    Incidence <- array(data=sim$wk.inc[,,SimulationLength], dim=c(Nages,nbox))
    FirstRow <- Incidence[1,1:ncols]
    LEnd <- ncols*nrows
    LStart <- ncols*nrows - ncols
    LastRow <- Incidence[1,LStart:LEnd]
    RightIndices <- (1:nrows)*ncols
    LeftIndices <- RightIndices - (ncols-1)
    FirstCol <- Incidence[1,LeftIndices]
    LastCol <- Incidence[1,RightIndices]
    PixParams <- PixelParams()
    ApproxPixelArea <- PixParams$pixel_width*PixParams$pixel_height
    cutoff <- animalCut/26*ApproxPixelArea
    if ( any(FirstRow > cutoff) || any(LastRow > cutoff) || any(FirstCol > cutoff) || any(LastCol > cutoff) ) {
      text <- paste("<font color=\"#FF0000\">", "Epidemic reached boundary of modeled region", "</font>")
    } else {
      text <- ""
    }
  })
  
  #===================================================== Leaflet Population Map =======================  
  output$mymap <- renderLeaflet({
    #leaflet() %>% setView(lng = input$loniC, lat = input$latC, zoom = 3) %>%
    leaflet() %>% setView(lng = longval, lat = latval, zoom = 3) %>%
      addWMSTiles("http://data.fao.org/maps/wms",
                  layers = "glbctd1f0503m_12713",
                  options = WMSTileOptions(
                    format = "image/png",
                    transparent = TRUE,
                    styles = "geonetwork_glbctd1f0503m_12713_style")) %>%
    addPolylines(data = admNatEarth, fill = FALSE, color = "black", weight = 1) 
      # Cmapadm
  })  #Make the leaflet population map
  
  #===================================================== Table of disease params =======================
  output$rates <- renderTable({
      Rr_IR <- (1-input$ReallySickFrac)/input$ReallySickFrac
      kD0 <- kD_input()   # This does not have to run to change value of kD0, assigns most recent value of kD_input()
      kD0dead <- kD0$kD0dead
      kD0recov <- kD0$kD0recov
      Rr0 <- Rr_input()
      values <- c(1/kD0dead[2,1,1], 1/(kD0dead[3,1,1]+kD0recov[3,1,1]), 1/(kD0dead[4,1,1]+kD0recov[4,1,1]), input$mortRate)  #
      ttest <- data.frame(values) 
      rownames(ttest) <- c("time in E (days)","time in I (days)","time in H (days)","fatality fraction from H")
      xalign(ttest, pad=TRUE)
      ttest
    }, rownames=TRUE)  
  
  #======================================================= Calculate disease progression ===============
  dis1 = reactive({
    nbox <- 9
    N = array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox))
    kD0=kD_input()    # Puts inputs into kD0
    Rr0=Rr_input()    # Puts inputs into Rr0
   
    klist <- SetupRates(kD0,Rr0,1)
    
    N[,,,]        <- 0    # N = array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox))          
    N[2,,1,9]    <- 100   # start with 100 exposed cattle for disease 
    record        <- array(data=0,dim=c(Ndstates,Ndcomb,Nages,9,nstep+1))
    record[,,,,1] <- N
    
    Vrate <- 0 # No vaccination
    Crate <- 0 # No culling
    betaRatio <- 0
    Ks <- Define_Ks(klist, Vrate, Crate, betaRatio, tstep)

    CellsToVacc <- NULL  # Runs NextN w/o vaccination
    CellsToCull <- NULL  # Runs NextN w/o culling
    
    dd.weight_totI <- array(data=0,dim=(c(nbox,nbox)))  # Just a dummy array for calling NextN
    dd.weight_totHplus <- array(data=0,dim=(c(nbox,nbox)))  # Just a dummy array for calling NextN
    MM <- list(dd.weight_totI=dd.weight_totI, dd.weight_totHplus=dd.weight_totHplus)
    
    for (month in 1:nmonth) {
      for (day in 1:len_mo) {   
        for (n in 1:(1/tstep)) { # 25 is approximate pixel area - but is only a dummy value not used in this calculation
          results <- NextN(N, 25, MM, tstep, klist, Ks, CellsToVacc, CellsToCull)
          N <- results$N
          idx <- ((day-1) + len_mo*(month-1))*(1/tstep) + n + 1 # remember we are not programming in C, 1 index is 0 timestep
          record[,,,,idx] <- N
        }
      }  #end day loop 
    } # end month loop 
    
    return(record)
  })
  
  #======================================================= Plot disease progression ====================
  output$diseaseProg=renderPlot({  
    record <- dis1()
    time <- seq(0:nstep)*tstep

    MakeDisProgPlot <- function() {
      plot(time,record[2,1,1,9,],main="Disease progression",xlim=c(0,30),ylim=c(0,100),col="green",ylab="percent",xlab="day") #E
      points(time,record[3,1,1,9,],col="blue")    # I  untreated
      points(time,record[4,1,1,9,],col="magenta") # H  untreated
      points(time,record[6,1,1,9,],col="red")     # D dead
      points(time,record[7,1,1,9,],col="black")   # R recovered
      lines(time,record[2,1,1,9,],col="green") 
      lines(time,record[3,1,1,9,],col="blue")      
      lines(time,record[4,1,1,9,],col="magenta")       
      lines(time,record[6,1,1,9,],col="red")
      lines(time,record[7,1,1,9,],col="black")  
      lines(time,colSums(record[2:7,1,1,9,]), col = "gray")  # sum of all - should be straight line
      legend("top",legend=c("Exposed","Infectious","Very Sick", "Dead", "recovered", "all"),text.col=c("green","blue","magenta", "red", "black", "gray"))
    }
    
      MakeDisProgPlot() 
      ##### Save as jpg for possible download ###
      dataplot <- paste(wwwDir, "diseaseProg.jpg", sep="/")
      try(jpeg(dataplot)) #Paul's suggestion
      MakeDisProgPlot()        
      dev.off()      
  })  
  
  #======================================================= Plot vaccine progression ====================
  output$vaccProg=renderPlot({  
    nstep <- 30/tstep
    time <- seq(0,nstep)*tstep
    
    ke <- 1/input$vacprog
    
    JV  <- 100*exp(-ke*time)
    Im  <- 100*(1-exp(-ke*time))

    MakeVaccProgPlot <- function() {
      plot(time,JV,main="Vaccine effectiveness progression",xlim=c(0,30),ylim=c(0,100),col="red",ylab="percent of vaccinated cattle",xlab="day") #Vaccinated still susceptible
      points(time,Im,col="green") # Vaccinated built up immunity
      lines(time,JV,col="red")
      lines(time,Im,col="green")
      legend("right",legend=c("Building immunity","Immune"),text.col=c("red","green"))
    }  
    
    MakeVaccProgPlot()    
    #### Save as jpg for possible download ###
    dataplot <- paste(wwwDir, "vaccProg.jpg", sep="/")
    try(jpeg(dataplot))
    MakeVaccProgPlot()  
    dev.off()
  })  
  
  #======================================================= Plot spread function ========================
  output$SpreadPlot=renderPlot({  
    rr <-  Pop2()           # function generates geographical population
    if (is.null(rr)) return(NULL)
    r <- rr2r()             # reduced res, geographical population, 
    
    shiny::validate(
      need(input$dl >= 0.1, "Characteristic spread length must be 0.1 km or larger")
    )
    
    PixParams <- PixelParams()
    pixel_size <- (PixParams$pixel_width + PixParams$pixel_height)/2
    dl <- input$dl
    
    # Assume ModelSize == 0.5 x 0.5 degrees  # 
    SpreadXYvals <- SpreadXY(input$latC, input$loniC, r, pixel_size, dl)
    
    MakeSpreadPlot <- function() {
      plot(SpreadXYvals$x,SpreadXYvals$y,main="Spread",col="gray",xlim=c(0,30),ylim=c(-8,0), ylab="Log10(Fraction of spread)",xlab="distance (km)") 
      lines(SpreadXYvals$x,SpreadXYvals$y,col="red") 
    }  
    
    MakeSpreadPlot() 
    #### Save as jpg for possible download ###
    dataplot <- paste(wwwDir, "SpreadPlot.jpg", sep="/")
    try(jpeg(dataplot))
    MakeSpreadPlot()  
    dev.off()
  }) 
  
  #======================================================= Epi simulation! ========================
  sim1=eventReactive(input$calculateEpidemics, {
    print("SIM1")
    shiny::validate(
      need(input$weeks_to_sim > 1, "Simulation must be 2 weeks or longer")
    )
    PixParams <- PixelParams()
    
    # input$ModelSize
    # if (input$ModelSize == "Small: 0.8 x 0.8 degrees") {
    #   gridCols <- 0.8/0.05
    # }
    # if (input$ModelSize == "Small Rect: 1 x 1.5 degrees") {
    #   gridCols <- 1.5/0.05
    # }
    # if (input$ModelSize == "Large Rect: 3 x 4 degrees") {
    #   gridCols <- 4/0.05
    # }
    # if (input$ModelSize == "Large: 3 x 3 degrees") {
    #   gridCols <- 3/0.05
    # }
    # if (input$ModelSize == "Medium: 2 x 2 degrees") {
    #   gridCols <- 2/0.05
    # }
    # shiny::validate(
    #   need(PixParams$pixel_width*gridCols > input$dl*4.6, "The disease spread is outside of the modeled geographic region, use a larger geographic area or reduce the characteristic distance of disease spread.")
    # )
    
    if ( (input$ModelType == "Use Case") && (input$Example  == "Rinderpest: Pakistan 1994") ) {
      shiny::validate(
        need(input$weeks_to_sim < 43, "Simulation of 1994 Pakistan epidemic is for only 42 weeks.")
      )
    }  
    Disease <- input$Disease
    N <- Agg()              # initializes population matrix of disease states
    r <- rr2r()             # geographical population,  
    dd <- changeDD()       # Function, Matrix of distances between geographical units(squares)
    IndexCell <- icell()    # number of cell where middle coordinates are, usually epidemic starts there

    kD0 <- kD_input()       # Reactive
    Rr0 <- Rr_input()       # Reactive
    beta <- input$beta 
    betaRatio <- 0
    #SimType = "regular"
    
    SpreadConstant <- input$dl
    SickCasesWhenDetect <- input$DetectSick
    
    # Parameters not passed to sim1, but needed for csv file
    UsedParams <- list(ModelSize=input$ModelSize, latC=input$latC, loniC=input$loniC, AggFact=input$AggFact, Lloni=input$Lloni, Llati=input$Llati, Vbase=input$Vbase, 
                       mortRate=input$mortRate, treated=input$treated, mortRateH=input$mortRateH, treatedH=input$treatedH, ModelType=input$ModelType, DiseaseModeled = input$Disease)      
    
    TimeParams  <- list(weeks_to_sim=input$weeks_to_sim,nyear=nyear,len.wk=len.wk,tstep=tstep)
    VacParams   <- list(Vrate=input$Vrate,Vdelay=input$Vdelay,Veff=1, VacRadius=input$vacrad, VrateH=input$VrateH, VdelayH=input$VdelayH, VacRadiusH=input$vacradH)
    Crate <- input$Crate
    CullParams  <- list(Crate=Crate, CullFrom=input$CullFrom, Culldelay=input$Culldelay, Cullrad=input$CullRad, CdelayPInfect=input$CdelayPInfect)
    
    srmcDelay       <- input$srmcd  # Delay of short range movement control after disease detection
    srmcBetaChange  <- input$srmc   # Multiplicative factor to reduce beta due to (short range) movement control

    LDType <- reactLDType()
    if (input$ModelType == "Use Case") {  # Don't want && because Example might not be defined
      if (input$Example  == "Rinderpest: Pakistan 1994") {
        roads <- getRoads()                                   # The Pakistan roads in given region
        Fixes <- FixRoads()
        RL <- Fixes$RL
        index_KeptRoads <- Fixes$index_KeptRoads
        road_list <- AboutRoads()                             # For each road a list of cells on that road and distances to the cells.
        IntersectInfo <- Intersects()
        
        RoadSpread <- input$RoadDist
        cowTypeMov <- input$cowTypeMov
        LDProb <- input$LDProb
        lrmc <- input$lrmc
        lrmcd <- input$lrmcd
        
        LDMoveParams <- list(LDType=LDType, RL=RL, index_KeptRoads=index_KeptRoads, road_list=road_list, IntersectInfo=IntersectInfo, LDProb=LDProb, RoadSpread=RoadSpread, cowTypeMov=cowTypeMov, lrmc=lrmc, lrmcd=lrmcd)
      } else {
        Lcellnums <- ldcell()
        cowTypeMov <- input$cowTypeMov
        LDProb <- input$LDProb
        lrmc <- input$lrmc
        lrmcd <- input$lrmcd
        LDMoveParams <- list(LDType=LDType, Lcellnums=Lcellnums, LDProb=LDProb, cowTypeMov=cowTypeMov, lrmc=lrmc, lrmcd=lrmcd)
      }  
    } else {
        Lcellnums <- ldcell()
        cowTypeMov <- input$cowTypeMov
        LDProb <- input$LDProb
        lrmc <- input$lrmc
        lrmcd <- input$lrmcd
        LDMoveParams <- list(LDType=LDType, Lcellnums=Lcellnums, LDProb=LDProb, cowTypeMov=cowTypeMov, lrmc=lrmc, lrmcd=lrmcd)
    }
    
    AddRow <- function(FileName, VariableName, VariableValue) {
      next_row <- paste(VariableName, VariableValue, sep = ",")
      FileNamePlus <- rbind(FileName,next_row)
      return(FileNamePlus)
    }
    
    WriteCSVFile <- function(SpreadConstant,TimeParams,beta,betaRatio,Rr0,kD0,SickCasesWhenDetect,srmcDelay,srmcBetaChange,VacParams,CullParams,LDMoveParams,UsedParams,DiseaseModeled) {
      restart_file <- "VariableName,value"
      restart_file <- AddRow(restart_file,"Disease",UsedParams$DiseaseModeled)
      restart_file <- AddRow(restart_file,"ModelSize",UsedParams$ModelSize)  
      restart_file <- AddRow(restart_file,"center latitude", UsedParams$latC) 
      restart_file <- AddRow(restart_file,"center longitude", UsedParams$loniC) 
      restart_file <- AddRow(restart_file,"weeks to simulate", TimeParams$weeks_to_sim)
      restart_file <- AddRow(restart_file,"DetectSick (cases)", SickCasesWhenDetect)
      restart_file <- AddRow(restart_file,"Aggregation", as.double(UsedParams$AggFact))   #### ?????
      
      # Disease parameters: 8 thru 13
      #kD1dead   <- c(0, kEI_Rind, kIH_Rind, kHD_Rind, kItHt_Rind, kHtD_Rind, 0, 0, 0)
      #kD0 <- list(kD0dead=kD0dead, kD0recov=kD0recov)
      restart_file <- AddRow(restart_file,"beta", beta)
      restart_file <- AddRow(restart_file,"Incubation period of disease (days)", 1/kD0$kD0dead[2,1,1])
      restart_file <- AddRow(restart_file,"kIH (1/days)", kD0$kD0dead[3,1,1])
      restart_file <- AddRow(restart_file,"kIR (1/days)", kD0$kD0recov[3,1,1])
      restart_file <- AddRow(restart_file,"kHD (1/days)", kD0$kD0dead[4,1,1])
      restart_file <- AddRow(restart_file,"kHR (1/days)", kD0$kD0recov[4,1,1])
      kHHt <- (kD0$kD0dead[4,1,1] + kD0$kD0recov[4,1,1])*Rr0[4,1,1]
      restart_file <- AddRow(restart_file,"kHHt (1/days)", kHHt)
      restart_file <- AddRow(restart_file,"kHtD (1/days)", kD0$kD0dead[6,1,1])
      restart_file <- AddRow(restart_file,"kHtR (1/days)", kD0$kD0recov[6,1,1])

      # "Short" range isotropic movement parameters: 14 thru 16
      restart_file <- AddRow(restart_file,"length for short range spread decay (km)", SpreadConstant)
      restart_file <- AddRow(restart_file,"Fractional reduction in beta due to short range movement controls (srmc's)", srmcBetaChange)
      restart_file <- AddRow(restart_file,"Delay after disease detection for implementation of srmc's (days)", srmcDelay)

      # Long Distance movement parameters: 
      restart_file <- AddRow(restart_file,"Percent of movement that is long distance", LDMoveParams$LDProb)
      restart_file <- AddRow(restart_file,"Longitude of long distance point", UsedParams$Lloni)    
      restart_file <- AddRow(restart_file,"Latitude of long distance point", UsedParams$Llati)    
      restart_file <- AddRow(restart_file,"Characteristic of animal or human moved long distance", LDMoveParams$cowTypeMov)
      restart_file <- AddRow(restart_file,"long distance movement reduction", LDMoveParams$lrmc)
      restart_file <- AddRow(restart_file,"Delay after disease detection for implementation of long distance movement control (days)", LDMoveParams$lrmcd)

      # Culling parameters: 
      restart_file <- AddRow(restart_file,"# animals culled per day", CullParams$Crate)
      if (CullParams$Crate != 0) {
        restart_file <- AddRow(restart_file,"Delay after disease detection before culling (days)", CullParams$Culldelay)
        restart_file <- AddRow(restart_file,"States culled from", CullParams$CullFrom)
      } else {
        restart_file <- AddRow(restart_file,"Delay after disease detection before culling (days)", "NA")
        restart_file <- AddRow(restart_file,"States culled from", "NA")
      } 
      
      # Vaccination parameters: 
      restart_file <- AddRow(restart_file,"Base vaccination level (%)", UsedParams$Vbase)
      restart_file <- AddRow(restart_file,"Vaccination rate (/day)", VacParams$Vrate)
      if (VacParams$Vrate == 0) {
        restart_file <- AddRow(restart_file,"Delay after disease detection before vaccination starts (days)", "NA")
        restart_file <- AddRow(restart_file,"Vaccination radius around each case (km)", "NA") 
        restart_file <- AddRow(restart_file,"Time to build immunity (days)", "NA") 
      } else {
        restart_file <- AddRow(restart_file,"Delay after disease detection before vaccination starts (days)", VacParams$Vdelay)
        restart_file <- AddRow(restart_file,"Radius around each case that is vaccinated (km)", VacParams$VacRadius) 
        restart_file <- AddRow(restart_file,"Time to build immunity (days)", 1/kD0recov[7,1,1]) 
      }
      
    return(restart_file)
    }  
    
    PixParams <- PixelParams()
    
    CSV_File <- WriteCSVFile(SpreadConstant,TimeParams,beta,betaRatio,Rr0,kD0,SickCasesWhenDetect,srmcDelay,srmcBetaChange,VacParams,CullParams,LDMoveParams,UsedParams,DiseaseModeled)
    fname <- paste(wwwDir, "inputs.csv", sep="/")
    write(CSV_File, file = fname, ncolumns = 1)
    
    # ClearGraph <- function(GraphName) {
    #   if(file.exists(GraphName)) {
    #     file.remove(GraphName)
    #   }
    # }
    # ClearGraph("www/OrigPop.jpg")      
    # ClearGraph("www/CumDeaths.jpg")
    # ClearGraph("www/Incidence.jpg")
    # ClearGraph("www/Prevalance.jpg")
    # ClearGraph("www/VaccImmune.jpg")
    # ClearGraph("www/Cumcases.jpg")
    
    sim=runSim(N,r,SpreadConstant,dd,IndexCell,TimeParams,beta,betaRatio,Rr0,kD0,SickCasesWhenDetect,srmcDelay,srmcBetaChange,
               VacParams,CullParams,LDMoveParams, UsedParams, PixParams)      

  }, ignoreNULL=FALSE)
  
  #####======================================================= Epidemic Line Plot ========================####
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  observeEvent(input$ggplotCh_dblclick, {
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    brush <- input$ggplotCh_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      sim=sim1()
      len.sim=dim(sim$wk.prev)[5]
      nbox=dim(sim$wk.prev)[4]
      
      Incident <- sim$wk.inc
      ExposedAnimals <- array(data=sim$wk.prev[2,1,,,],dim=c(Nages,nbox,len.sim))     #E
      InfectiousAnimals <- array(data=sim$wk.prev[3,1,,,],dim=c(Nages,nbox,len.sim))  #I
      ReallysickAnimals <- array(data=sim$wk.prev[4,1,,,],dim=c(Nages,nbox,len.sim))  #H
      DeadAnimals <- array(data=sim$wk.prev[6,1,,,],dim=c(Nages,nbox,len.sim))  #D
      All <- array(data=colSums(sim$wk.prev[,1,,,]),dim=c(Nages,nbox,len.sim))
      RecoveredAnimals <- array(data=sim$wk.prev[7,1,,,],dim=c(Nages,nbox,len.sim))
      VaccinatedAnimals <- array(data=sim$wk.prev[9,1,,,],dim=c(Nages,nbox,len.sim))

      maxy <- max(colSums(ExposedAnimals,dims=2), colSums(InfectiousAnimals,dims=2), colSums(ReallysickAnimals,dims=2),
                  colSums(DeadAnimals,dims=2),colSums(RecoveredAnimals,dims=2),colSums(VaccinatedAnimals,dims=2), colSums(All,dims=2), 10)  
      
      ranges$y <- c(0, log10(maxy))
      #ranges$y <- NULL
    }
  })
  
  output$ggplotCh <- renderPlot({
    sim=sim1()
    if (sim$UsedParams$DiseaseModeled != input$Disease) return(NULL)
    if (sim$UsedParams$ModelType != input$ModelType) return(NULL)
    len.sim=dim(sim$wk.prev)[5]
    nbox=dim(sim$wk.prev)[4]
    time=seq(1,len.sim)
    
    rr <-  Pop2()           # function generates geographical population
    r <- rr2r()             # reduced res, geographical population,  
    PixParams <- PixelParams()
    ApproxPixelArea <- PixParams$pixel_width*PixParams$pixel_height
    
    #wk.inc=array(data=0,dim=c(Nages,nbox,weeks_to_sim))
    #wk.prev=array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,weeks_to_sim))
    #Incident <- array(data=sim$wk.inc[,,], dim=c(Nages,nbox,len.sim))
    
    Incident <- sim$wk.inc
    ExposedAnimals <- array(data=sim$wk.prev[2,1,,,],dim=c(Nages,nbox,len.sim))     #E
    InfectiousAnimals <- array(data=sim$wk.prev[3,1,,,],dim=c(Nages,nbox,len.sim))  #I
    ReallysickAnimals <- array(data=sim$wk.prev[4,1,,,],dim=c(Nages,nbox,len.sim))  #H
    DeadAnimals <- array(data=sim$wk.prev[6,1,,,],dim=c(Nages,nbox,len.sim))  #D
    All <- array(data=colSums(sim$wk.prev[,1,,,]),dim=c(Nages,nbox,len.sim))
    RecoveredAnimals <- array(data=sim$wk.prev[7,1,,,],dim=c(Nages,nbox,len.sim))
    VaccinatedAnimals <- array(data=sim$wk.prev[9,1,,,],dim=c(Nages,nbox,len.sim))
    
    sickmax <- max(colSums(ReallysickAnimals,dims=2))  
    expmax <- max(colSums(ExposedAnimals,dims=2))  
    
    ExposedEnd <- array(data=sim$wk.prev[2,1,1,,len.sim],dim=nbox)
    DeathsEnd <- array(data=sim$wk.prev[6,1,1,,len.sim],dim=nbox)
    RecoveriesEnd <- array(data=sim$wk.prev[7,1,1,,len.sim],dim=nbox)
    notEnoughEnd <- which(ExposedEnd < animalCut/26*ApproxPixelArea)
    notStartedEnd <-  which((DeathsEnd+RecoveriesEnd) == 0)
    NP <- intersect(notEnoughEnd, notStartedEnd)
    
    print("Number exposed but not progressing at end of simulation")
    print(sum(ExposedAnimals[1,NP,len.sim]))  # Non - progressing at end
    
    y = log10(colSums(All,dims=2))
    y2 = log10(colSums(ExposedAnimals,dims=2))      # E
    y3 = log10(colSums(InfectiousAnimals,dims=2))   # I
    y4 = log10(colSums(ReallysickAnimals,dims=2))   # H
    y6 = log10(colSums(DeadAnimals,dims=2))         # D
    y7 = log10(colSums(RecoveredAnimals,dims=2))    # R
    y9 = log10(colSums(VaccinatedAnimals,dims=2))   # Vm 
    y10 = log10(colSums(Incident,dims=2))
    df <- data.frame(time,y,y2,y3,y4,y6,y7,y9,y10)
    
    test <- ranges$y
    if (is.null(test)) {
        maxy <- max(colSums(ExposedAnimals,dims=2), colSums(InfectiousAnimals,dims=2), colSums(ReallysickAnimals,dims=2),
                    colSums(DeadAnimals,dims=2),colSums(RecoveredAnimals,dims=2),colSums(VaccinatedAnimals,dims=2), colSums(All,dims=2), 10)  
        limOfy <- c(0, log10(maxy))
    }
    else {
      limOfy <- ranges$y
    }
    
    AnimalType <- "Cattle"
    labely <- paste("log10( Number of", AnimalType, ")" )
    
    epidemic_plot <- ggplot(df) + geom_line(aes(x=time,y=y, colour = "Total")) +
        geom_line(aes(x=time,y=y2, colour = "Exposed")) +
        geom_line(aes(x=time,y=y3, colour = "Infectious")) +
        geom_line(aes(x=time,y=y4, colour = "ReallySick")) + 
        geom_line(aes(x=time,y=y6, colour = "Dead")) +
        geom_line(aes(x=time,y=y7, colour = "Recovered")) +
        geom_line(aes(x=time,y=y9, colour = "Vaccinated")) + 
        geom_line(aes(x=time,y=y10, colour = "Incident")) + 
        ylab(labely) + xlab("Time (weeks)") + 
        labs(title="Epidemic Progression") + 
        scale_colour_manual("Legend", 
                          values = c(Total = "gray", Exposed = "green", Infectious = "blue", ReallySick = "purple", Dead = "red", Recovered = "brown", Vaccinated = "black", Incident = "orange")) +
        coord_cartesian(xlim = ranges$x, ylim = limOfy) +
        theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=18),
              axis.title.y = element_text(size=20), axis.text.y  = element_text(size=18),
              legend.title = element_text(size=18),
              legend.text = element_text(size = 18),
              plot.title = element_text(size=20))
    
    pdf("www/EpidemicProgPlot.pdf")
    print(epidemic_plot)
    dev.off()
    
    epidemic_plot
  })  #Cases and deaths versus time
  
  #======================================================= Output Map of Effected Population ========================
  output$epiOutMap=renderPlot({
    corners <- Rcorners()
    sim=sim1()

    if ((sim$UsedParams$DiseaseModeled != input$Disease) && (input$epiOutMapButton != "orig")) {
      return(NULL)
    }
    if ((sim$UsedParams$ModelType != input$ModelType) && (input$epiOutMapButton != "orig")) {
      return(NULL)
    }
    
    Simulation_length <- length(sim$wk.inc[1,1,])
    if ( (input$P.wk > Simulation_length) &&  (input$epiOutMapButton != "orig") ) {
      return(NULL)
    }
    nbox=dim(sim$wk.prev)[4]
    
    #Rout  <- rr2r()  
    Rout <- MakeRasterToPlot()
    
    Nbox_simmed <- length(sim$wk.inc[1,,1])
    Nbox_now <- length(Rout)
    if ((Nbox_simmed != Nbox_now) && (input$epiOutMapButton != "orig")) {
      return (NULL)
    }
    
    # Information for the map title
    disease <- input$Disease
    animal <- "Cattle"
    
    IndexCell <- icell()
    #LDCell <- ldcell()
    PixParams <- PixelParams()
    ApproxPixelArea <- PixParams$pixel_width*PixParams$pixel_height
    AspectRatio <- PixParams$pixel_height/PixParams$pixel_width
    
    FontSizesDisplay <- c(1.3,1.1,1.2)
    MarginDisplay <- c(5,5,4,2) #Units of lines #c(0.5,0.5,0.5,0.2) # Units of inches# 
    FontSizesPresent <- c(1.7,1.5,1.6) #c(2.0,1.5,2.0)
    MarginPresent <- c(5,5,4,2) #Units of lines  c(0.5,0.5,0.5,0.2) # Units of inches# 
    ColorScaleMarginDisplay <- 5
    ColorScaleMarginPresent <- 7
    
    MakePlot <- function(Rout, mainTitle, zlimit, FontSizes, MarginSizes, ColorScaleMargin, ModelType, LDProb, Lloni, Llati) {  # 
      par(cex.main = FontSizes[1])  # Make main title a larger font
      par(cex.lab = FontSizes[3])   # Make axis labels a larger font
      par(cex.axis = FontSizes[3])  # Make axis numbers a larger font
      #par(pty="m")
      par(mar = MarginSizes)
      if ((input$epiOutMapButton == "orig") || (input$epiOutMapButton == "vacimmune") || (input$epiOutMapButton == "prevalence")) {
        cols <- rev(terrain.colors(255))
      } else {
        cols <- c("white",rev(rainbow(100)[1:75]))
      }  
      #if ((input$epiOutMapButton == "vacimmune") || (input$epiOutMapButton == "incidence") || (input$epiOutMapButton == "orig") || (input$epiOutMapButton == "Cumulative cases")) { 
        plot(Rout, main=mainTitle, xlab="longitude", ylab="latitude", zlim = zlimit, col=cols, legend.mar = ColorScaleMargin)
      #} else {
      #  plot(Rout, main=mainTitle, xlab="longitude", ylab="latitude")
      #}  
      if (input$ModelType == "Use Case") {
        if (input$Example  == "Rinderpest: Pakistan 1994") {
          LDType <- "roads"
          roads <- getRoads()  #The Pakistan roads in given region
          points(xFromCell(Rout, IndexCell), yFromCell(Rout, IndexCell), pch=8, col="red") #start location for epidemic
          if (!is.null(roads)) {  
            plot(roads, add=TRUE, col=gray(0.4))  # chocolate4
            legend("topleft", bty="n", legend=c("start location"), pch=8, col=c("red"), cex=FontSizes[2]) 
            text(74.2, 35.8, "Gilgit", cex=FontSizes[2])
            text(75.2, 37.2, "Roads are Gray", cex=FontSizes[2])
            text(76.19,35.05, "Khaplu", cex=FontSizes[2])
            text(74.65,36.4, "Hunza Valley", cex=FontSizes[2])
            points(76.33,35.17, pch=20, col="blue", cex = FontSizes[2])  # Khaplu 
            points(74.30,35.93, pch=20, col="blue", cex = FontSizes[2])  # Gilgit
          }
        }  
      } else {
        points(xFromCell(Rout, IndexCell), yFromCell(Rout, IndexCell), pch=8, col="red", cex=FontSizes[2]) #start location for epidemic
        if (input$LDProb > 0) {
          points(input$Lloni,input$Llati, pch=8, col="blue") #Long distance point
          legend("topleft", bty="n", legend=c("index case","Long distance"), pch=8, col=c("red","blue"))  
        } 
        else{
          legend("topleft", bty="n", legend="index case", pch=8, col="red", cex=FontSizes[2]) 
        }
      }
    }
    
    TitleAndZlimits <- PlotTitleAndZlimit(input$epiOutMapButton, Rout, input$P.wk)
    zlimit <- TitleAndZlimits$zlimit
    mainTitle <- TitleAndZlimits$mainTitle
    MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay, sim$ModelType, input$LDProb, input$Lloni,input$Llati)
    
    # if (input$epiOutMapButton == "orig") { 
    #   mainTitle <- "Original population \n at calculation resolution"
    #   zlimit = c(1,NA)
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay)                                               
    #   
    #   # Save as jpg for possible download #
    #   dataplot <- paste(wwwDir, "OrigPop.jpg", sep="/")
    #   totalPlotHeight <- 6 # units are inches
    #   PlotHeight <- totalPlotHeight-(MarginPresent[1]+MarginPresent[3]) # inches
    #   PlotWidth <- PlotHeight/AspectRatio
    #   totalPlotWidth <- PlotWidth + (MarginPresent[2]+MarginPresent[4])
    #   #WidthHeightRatio <- (dim(Rout)[2]*PixParams$pixel_width)/(dim(Rout)[1]*PixParams$pixel_height)
    #   #print(WidthHeightRatio)
    #   #try(jpeg(dataplot, width = totalPlotWidth, height = totalPlotHeight, units="in", res=200))  # 200 ppi
    #   try(jpeg(dataplot, width=550, height=480))
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesPresent, MarginPresent, ColorScaleMarginPresent)
    #   dev.off()
    # } 
    
    # if (input$epiOutMapButton == "deaths") {  #  wk.prev=array(data=0,dim=c(Ndstates,Ndcomb,Nages,nbox,weeks_to_sim))
    #   Ptitle = "Cumulative Dead, wk "
    #   Deaths <- array(data=sim$wk.prev[6,1,,,input$P.wk],dim=c(Nages,nbox))
    #   values(Rout) <- colSums(Deaths)
    #   mainTitle <- paste(Ptitle,input$P.wk, sep="")
    #   zlimit = c(NA,NA)
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay)
    #   
    #   dataplot <- paste(wwwDir, "CumDeaths.jpg", sep="/")                                 # Save as jpg for possible download
    #   try(jpeg(dataplot, width=550, height=480))
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesPresent, MarginPresent, ColorScaleMarginPresent)
    #   dev.off()
    # } 

    # if (input$epiOutMapButton == "incidence") {  #wk.inc=array(data=0,dim=c(Ndcomb,Nages,nbox,weeks_to_sim))
    #   AnimalIncidence <- array(data=sim$wk.inc[,,input$P.wk] ,dim=c(Nages,nbox))
    #   #AnimalIncidence <- log10(AnimalIncidence)
    #   #AnimalIncidence[which(is.infinite(AnimalIncidence), arr.ind=TRUE)] <- NA   #  THIS causes RENDERPLOT to give an ERROR, but is OK.
    #   values(Rout)=colSums(AnimalIncidence) # Summing over ages
    #   #Ptitle = "log10(Incidence), wk "
    #   Ptitle = "Incidence, wk "
    #   mainTitle <- paste(Ptitle,input$P.wk, sep="")
    #   maxZval <- ceiling(max(values(Rout), na.rm = TRUE))
    #   PixParams <- PixelParams()
    #   ApproxPixelArea <- PixParams$pixel_width*PixParams$pixel_height
    #   zlimit <- c(0, maxZval) 
    #   # zlimit <- c(-1,maxZval)
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay)
    #   
    #   ##### Save as jpg for possible download ###
    #   dataplot <- paste(wwwDir, "Incidence.jpg", sep="/")
    #   try(jpeg(dataplot, width=550, height=480))
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesPresent, MarginPresent, ColorScaleMarginPresent)
    #   dev.off()
    # }
    
    # if (input$epiOutMapButton == "prevalence") {  #wk.inc=array(data=0,dim=c(Ndcomb,Nages,nbox,weeks_to_sim))
    #   SickPop <- array(data=colSums(sim$wk.prev[2:5,1,,,input$P.wk]),dim=c(Nages,nbox))
    #   TotalPop <- array(data=colSums(sim$wk.prev[1:5,1,,,input$P.wk])+colSums(sim$wk.prev[7:9,1,,,input$P.wk]),dim=c(Nages,nbox))
    #   values(Rout) <- colSums(SickPop)/colSums(TotalPop) #Summing over ages
    #   Ptitle = "Prevalence, wk "
    #   mainTitle <- paste(Ptitle,input$P.wk, sep="")
    #   zlimit <- c(0,1)
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay)
    #   ##### Save as jpg for possible download ###
    #   dataplot <- paste(wwwDir, "Prevalence.jpg", sep="/")
    #   try(jpeg(dataplot, width=550, height=480))
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesPresent, MarginPresent, ColorScaleMarginPresent)
    #   dev.off()
    # }
    
  #   if (input$epiOutMapButton == "vacimmune") {
  #     Immune <- array(data=sim$wk.prev[9,1,,,input$P.wk],dim=c(Nages,nbox))
  #     Sdata <- sim$wk.prev[1,1,,,input$P.wk] + sim$wk.prev[9,1,,,input$P.wk] + sim$wk.prev[8,1,,,input$P.wk]
  #     Susceptible <- array(data=Sdata,dim=c(Nages,nbox))
  #     values(Rout) <- Immune/colSums(Susceptible)
  #     Ptitle = "Fraction of potentially susceptibles \nimmune by vaccination on wk "
  #     mainTitle <- paste(Ptitle,input$P.wk, sep="")
  #     zlimit <- c(0,1)
  #     MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay)
  #     
  #     ##### Save as jpg for possible download ###
  #     dataplot <- paste(wwwDir, "VaccImmune.jpg", sep="/")
  #     try(jpeg(dataplot, width=550, height=480))
  #     MakePlot(Rout, mainTitle, zlimit, FontSizesPresent, MarginPresent, ColorScaleMarginPresent)
  #     dev.off()
  # }
  #   
    # if (input$epiOutMapButton == "Cumulative cases") {
    #   AnimalIncidence <- array(data=rowSums(sim$wk.inc[1,,1:input$P.wk]) ,dim=c(Nages,nbox))
    #   values(Rout)=colSums(AnimalIncidence) #Summing over ages
    #   Ptitle = "Cumulative cases, wk "
    #   mainTitle <- paste(Ptitle,input$P.wk, sep="")
    #   zlimit <- c(0.3,NA)
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesDisplay, MarginDisplay, ColorScaleMarginDisplay)
    #   
    #   ##### Save as jpg for possible download ###
    #   dataplot <- paste(wwwDir, "Cumcases.jpg", sep="/")
    #   try(jpeg(dataplot, width=550, height=480))
    #   MakePlot(Rout, mainTitle, zlimit, FontSizesPresent, MarginPresent, ColorScaleMarginPresent)
    #   dev.off()
    # }
    
  })  # Geographically resolved plot
  
  
  ####### Output Single maps #####
  output$ImageMap <- downloadHandler(
    filename = function() {
      # print("In Image Map download")
      # print(input$epiOutMapButton)
      if (input$epiOutMapButton == "orig") startname <- "OriginalPopulation"
      if (input$epiOutMapButton == "deaths") startname <- "CumDeaths_wk"
      if (input$epiOutMapButton == "incidence") startname <- "Incidence_wk"
      if (input$epiOutMapButton == "Cumulative cases") startname <- "Cumcases_wk"
      if (input$epiOutMapButton == "prevalence") startname <- "Prevalence_wk"
      if (input$epiOutMapButton == "vaccinated") startname <- "Vaccinated_wk"
      if (input$epiOutMapButton == "vacimmune") startname <- "ImmuneByVaccine_wk"
      
      fname <- paste(startname, input$P.wk, ".jpg", sep="")
      return(fname)
    },  
    content = function(file) {
      sim <- sim1()
      IndexCell <- icell()
      Rout <- MakeRasterToPlot()
      
      MakePlot <- function(Rout, mainTitle, zlimit, ModelType, Example, LDProb, Lloni, Llati) {
        par(cex.lab = 1.8) # axes labels 2 times larger than default
        par(cex.axis = 1.6)
        par(cex.main = 1.5)
        par(mar=c(5.1,5.1,4.1,2.1)) # makes left margin a little bigger
        if ((input$epiOutMapButton == "orig") || (input$epiOutMapButton == "vacimmune")) {
          cols <- rev(terrain.colors(255))
        } else {
          cols <- c("white",rev(rainbow(100)[1:75]))
        } 
        plot(Rout, main=mainTitle, xlab="longitude", ylab="latitude", col=cols, zlim = zlimit)
        if (ModelType == "Use Case") {
          if (Example  == "Rinderpest: Pakistan 1994") {
            LDType <- "roads"
            roads <- getRoads()  #The Pakistan roads in given region
            points(xFromCell(Rout, IndexCell), yFromCell(Rout, IndexCell), pch=8, col="red") #start location for epidemic
            if (!is.null(roads)) {  
              plot(roads, add=TRUE, col="gray")  #chocolate4
              legend("topleft", bty="n", legend=c("start location"), pch=8, col=c("red"), cex=1.3) 
              text(74.2, 35.8, "Gilgit", cex = 1.5)
              text(75.2, 37.2, "Roads are Gray", cex = 1.5)
              text(76.19,35.05, "Khaplu", cex = 1.5)
              text(74.65,36.4, "Hunza Valley", cex = 1.5)
              points(76.33,35.17, pch=20, col="blue")  # Khaplu 
              points(74.30,35.93, pch=20, col="blue")  # Gilgit
            }
          }  
        } else {
          points(xFromCell(Rout, IndexCell), yFromCell(Rout, IndexCell), pch=8, col="red") #start location for epidemic
          if (LDProb > 0) {
            points(Lloni,Llati, pch=8, col="blue") #Long distance point
            legend("topleft", bty="n", legend=c("starting location","Long distance"), pch=8, col=c("red","blue"), cex=1.5)  
          } 
          else{
            legend("topleft", bty="n", legend="starting location", pch=8, col="red", cex=1.5) 
          }
        }
      }
      
      TitleAndZlimits <- PlotTitleAndZlimit(input$epiOutMapButton, Rout, input$P.wk)
      zlimit <- TitleAndZlimits$zlimit
      mainTitle <- TitleAndZlimits$mainTitle
      
      #MakePlot(Rout, mainTitle, zlimit, sim$ModelType, input$Example, input$LDProb, input$Lloni,input$Llati)  # !!! Potential problem using inputs !!!!
      
      #dataplot <- paste(wwwDir, file, sep="/")
      #print(file)
      #print(sim$UsedParams$ModelType)
      try(jpeg(file))
      MakePlot(Rout, mainTitle, zlimit, sim$UsedParams$ModelType, input$Example, input$LDProb, input$Lloni,input$Llati)  # !!! Potential problem using inputs !!!!
      dev.off()
    }
  )
  
  ###### Consequence table #####
  output$consequence <- renderTable({
    sim=sim1()
    Simulation_length <- length(sim$wk.inc[1,1,])
    if (sim$UsedParams$DiseaseModeled != input$Disease) return(NULL)
    if (sim$UsedParams$ModelType != input$ModelType) return(NULL)
    if (input$P.wk > Simulation_length)  return(NULL)    
    if ((sim$UsedParams$latC != input$latC) || (sim$UsedParams$loniC != input$loniC)) {
      return(NULL)
    }
    
    rr <- Pop2()
    if (is.null(rr)) return(NULL)
    r <- rr2r()
    
    initPop <- round(sum(getValues(r), na.rm=TRUE))
    PixParams <- PixelParams()
    ApproxPixelArea <- PixParams$pixel_width*PixParams$pixel_height
    time <- input$P.wk

    cumInc  <- round(sum(sim$wk.inc[1,,1:time])) # Nages,nbox,weeks_to_sim
    totD <- round(sum(sim$wk.prev[6,1,,,time])) # Ndstates,Ndcomb,Nages,nbox,weeks_to_sim
    totV <- round(sum(sim$wk.prev[9,1,,,time]))
    totR <- round(sum(sim$wk.prev[7,1,,,time]))
    infectedBoxesandTimes <- which(sim$wk.prev[2,1,1,,1:time] >= animalCut/26*ApproxPixelArea, arr.ind=TRUE)
    
    infectedBoxes <- unique(infectedBoxesandTimes)  #First column is the boxes
    Area <- round(length(infectedBoxes)*ApproxPixelArea)
    
    results <- c(initPop, cumInc, totR, totD, totV, Area)
    ttest <- data.frame(results) 
    rownames(ttest) <- c("Initial population modeled", paste("Cumulative cases at week", time), paste("Total recovered at wk", time), paste("Total dead at wk", time), 
                         paste("Total immune by vaccination at wk", time), paste("area affected (km^2) at wk", time))
    
    totDead <- colSums(sim$wk.prev[6,1,1,,], dims=1)
    totInc <- colSums(sim$wk.inc[1,,], dims=1)
    CumCases <- cumsum(totInc)
    time_results <- data.frame(1:input$weeks_to_sim,totInc[1:input$weeks_to_sim],CumCases[1:input$weeks_to_sim],totDead[1:input$weeks_to_sim])
    
    fname <- paste(wwwDir, "results.csv", sep="/")
    write.csv(ttest, file = fname, quote=FALSE)
    
    fname <- paste(wwwDir, "time_results.csv", sep="/")
    write.csv(time_results, file = fname, quote=FALSE)
    
    ttest
  }, rownames=TRUE, digits=0)

  #===================================================== File Download =======================
  
  output$downloadData <- downloadHandler(
    filename = function() {
      start_fname <- paste("epiGridResults", sep="/")
      fname <- paste(start_fname, Sys.time(), ".zip", sep="")
      return(fname)
    },  
    content = function(file) {
      f1 <- paste(wwwDir, "EpidemicProgPlot.pdf", sep="/")
      
      #f2 <- paste(wwwDir, "OrigPop.jpg", sep="/")
      #f3 <- paste(wwwDir, "CumDeaths.jpg", sep="/")
      #f4 <- paste(wwwDir, "Incidence.jpg", sep="/")
      #f5 <- paste(wwwDir, "Prevalence.jpg", sep="/")
      #f6 <- paste(wwwDir, "VaccImmune.jpg", sep="/")
      #f7 <- paste(wwwDir, "Cumcases.jpg", sep="/")
      
      f8 <- paste(wwwDir, "inputs.csv", sep="/")
      f9 <- paste(wwwDir, "results.csv", sep="/")
      f10 <- paste(wwwDir, "time_results.csv", sep="/")
      
      #files <- c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)
      files <- c(f1,f8,f9,f10)
      print(which(file.mtime(files) >= file.mtime(f8)))
      
      fs <- files[which(file.mtime(files) >= file.mtime(f8))]
      try(zip(zipfile=file, files=fs))
    },
    contentType = "application/zip"
  )
  
  output$Documentation <- downloadHandler(   
    filename = function() {
      return("Rinderpest.pdf")
    },
    content = function(file) {
      myfile <- "Rinderpest.pdf"
      file.copy(myfile, file)
    },
    contentType = "application/pdf"
  )
  
  # output$Documentation <- downloadHandler(   
  #   filename = function() {
  #     return("EpiGridDoc.pdf")
  #   },
  #   content = function(file) {
  #     myfile <- "EpiGridDoc.pdf"
  #     file.copy(myfile, file)
  #   } 
  # )
 
})