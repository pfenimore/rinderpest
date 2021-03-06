#Rinderpest public ui.R

shinyUI(fluidPage(    
  h3("Simulate a Rinderpest epidemic", align = "center"),                                                                                                            
      fluidRow(
        ####################  Left Hand Column ################# 
          column(2,  
          actionButton("calculateEpidemics", label = "Simulate Epidemic", style="border-color: #ff0000"),   
          
          radioButtons("ModelType", label = "Parameter set:",
                      choices = c("Default disease parameters" = "No", "Examples" = "Use Case", "Load parameters" = "Load parameters"), inline = FALSE),
          
          conditionalPanel(condition = "input.ModelType == 'Use Case'",
                           selectInput("Example", "Example Use Case", choices = c("Rinderpest: Pakistan 1994", "Rinderpest: Fremantle"))),
                            
          conditionalPanel(condition = "input.ModelType == 'Load parameters'",
                           fileInput('inputs', 'Choose file', accept = '.csv'),
                           actionButton("csvFile", label = "Load csv file parameters")),
          
          selectInput("Disease", "Disease:", choices = DiseaseList),
          selectInput("ModelSize", "Size of modeled region", choices = SizeList),
          
          fluidRow(
            column(6,  
                   numericInput("latC", "latitude", (latmin+latmax)/2,
                                min = -60, max = 85, step = 0.01) #step=5009"
            ),
            column(6,
                   numericInput("loniC", "longitude", (longmin+longmax)/2,
                                min = -180, max = 180, step = 0.01) #step=500
            )
          ),  
          
          numericInput("weeks_to_sim", "Duration of simulation (weeks)", weeks_to_sim_BasicStart,
                       min = 2, max = 90, step=1
          ),
          numericInput("DetectSick", "Number sick when disease is detected", value = 5,#25,                                                                                                                       
                      min = 0, max = 500, step = 1
          ),
          selectInput("AggFact", "Multiplicative change to geo box length", choices = c(0.5, 1, 1.5, 2), selected = AggFact_BasicStart),
          textOutput("pixelsize"),
          textOutput("DetectionDay"),
          textOutput("DetectionSick"),
          br(),  # creates space
          textOutput("CullDay"),
          htmlOutput("BoundaryReached"),
          br(),  # creates space
          tags$a(href="https://www.frontiersin.org/articles/10.3389/fvets.2018.00182/full", "Download the Paper for this software"),
          # downloadButton('Documentation', 'Paper'),
          br(),  # creates space
          downloadButton('downloadData', 'Current inputs and results'),
          br()  # creates space
          #downloadButton('Documentation', 'Documentation')
         ),  #end of far left column
        
        ####################  Middle Columns ################# 
        column(7,
            column(6,
                h5("Cattle Population"),
                leafletOutput("mymap"),
                h5("Brush and double-click to zoom"),
                plotOutput("ggplotCh", height = 300,
                       dblclick = "ggplotCh_dblclick",
                       brush = brushOpts(
                         id = "ggplotCh_brush",
                         resetOnNew = TRUE
                       )
                )
            ),  #end of middle left column
            column(6,
                tableOutput("consequence"),
                plotOutput("epiOutMap", height = "400px"),
                fluidRow(
                  column(6,  
                      radioButtons("epiOutMapButton", "Choose quantity to plot",
                             c("Original Population" = "orig",
                               "Cumulative deaths" = "deaths",
                               "Weekly incidence" = "incidence",
                               "Prevalence" = "prevalence",
                               #"current susceptibles" = "current susceptibles",
                               #"Weekly fractional incidence" = "fincidence",
                               "Cumulative cases" = "Cumulative cases",
                               "immune by vaccination" = "vacimmune"),
                               selected="orig"
                    )
                  ), 
                  column(6,
                      sliderInput("P.wk", "Week of simulated epidemic to show in results table and image-map", min=1, max=20, step = 1, ticks = FALSE, value=4, animate=TRUE),
                      downloadButton('ImageMap', 'Current Image Map')
                  ) 
                )
            ) #end of middle right column
        ),  
        ####################  Right Columns ################# 
        column(3,
          tabsetPanel(
              ####################  Spread ################# 
              tabPanel("Spread", 
                    h4("Short Range"),
                    column(4,
                           numericInput("dl", "Characteristic distance (km)", value = dl_BasicStart, min = 0.5, max = 10, step = 0.1),
                           numericInput("srmcd", "Delay of control after epidemic detection (days)", value = 1, min = 0, max = 60, step = 1)
                    ),
                    column(8,
                           plotOutput("SpreadPlot", height = "300px")
                    ),
                    sliderInput("srmc", "Fractional reduction of transmission",                                                                                                                        
                                 min = 0,                                                                                                                                  
                                 max = 1,                                                                                                                                       
                                 value = srmc_BasicStart                                                                                                                                     
                    ),

                    h4("Spread to Distant Location"),
                    sliderInput("LDProb", "Percent of spread that is long distance",                                                                                                            
                                 min = 0,                                                                                                                                  
                                 max = 3.0,                                                                                                                                       
                                 value = 0,  #0.3, #15,#3,#1.5,
                                 step = 0.02
                    ),
                    conditionalPanel(condition = "(input.Example != 'Rinderpest: Pakistan 1994'  || input.ModelType != 'Use Case')  &&  input.LDProb > 0", 
                                sliderInput("Lloni", "Longitude: long range point", step = 0.01,
                                     min = round(corners$longmin,1)+0.1,#73.0,#0,                                
                                     max = round(corners$longmax,1)-0.1,#76.0,#20, 
                                     value = corners$longmax-0.4#74.3),#10), 
                                ),
                                sliderInput("Llati", "Latitude: long range point", step = 0.01,
                                   min = round(corners$latmin,1)+0.1, #round((latmin+latmax)/2,1) - 0.2, #corners$latmin+0.1,#73.0,#0,                                
                                   max = round(corners$latmax,1)-0.1, #round((latmin+latmax)/2,1) + 0.2, #corners$latmax-0.1,#76.0,#20, 
                                   value = corners$latmax-0.4#74.3),#10), 
                                )
                    ),
                    conditionalPanel(condition = "input.Example == 'Rinderpest: Pakistan 1994' && input.ModelType == 'Use Case'", 
                                sliderInput("RoadDist",                                                                                                                       
                                   "Distance cows move along roads (km)",                                                                                                                        
                                   min = 5,                                                                                                                                  
                                   max = 100,                                                                                                                                       
                                   value = 15, #15,#3,#1.5,
                                   step = 2.5
                            )
                    ),
                    conditionalPanel(condition = "input.LDProb > 0",
                          radioButtons("cowTypeMov", "Restrictions on long distance movement",
                                    c("No restrictions" = "all",
                                      "No seriously ill" = "NoH",
                                      "No infected" = "NoHI"),
                                    selected="all"
                          ),         
                          sliderInput("lrmcd", "Long-range movement control delay after epidemic detection (days)",                                                                                                                        
                                 min = 0,                                                                                                                                
                                 max = 60,                                                                                                                                       
                                 value = 0,
                                 step = 1 #5
                          ),
                          sliderInput("lrmc", "Fraction of long range movement occurring after controls in place",                                                                                                                        
                                 min = 0,                                                                                                                                  
                                 max = 1,                                                                                                                                       
                                 value = lrmc_BasicStart                                                                                                                                     
                          )
                      )
              ),
              ####################  Animal Culling ################# 
              tabPanel("A. Cull", 
                    h5("Culling infectious animals"),
                    numericInput("Crate", "Rate (/day)", 0, #1000, 
                                  min = 0, max = 10000, step=500
                    ),
                    conditionalPanel(
                        condition = "input.Crate > 0",
                        sliderInput("Culldelay", "Delay of culling after epidemic detection (days)",
                                   min = 0,
                                   max = 60,
                                   value = 14
                        ),
                        radioButtons("CullFrom", label = "Cull from which states:",
                                     choices = c("Cull really sick animals only (state H)" = "Honly", 
                                                 "Cull all infectious animals (states I and H)" = "IandH", 
                                                 "Cull all untreated and unvaccinated" = "All"), selected="All", inline = FALSE),
                        sliderInput("CullRad",
                                    "Radius of culling ring (km)",
                                    min = 1,
                                    max = 30,
                                    value = 5
                        ),
                        sliderInput("CdelayPInfect", "Delay between when animals have first detectable symptoms in a geographic block and when culling begins in and around that block (days)",
                                    min = 0,
                                    max = 7,
                                    value = 1
                        )
                        
                )  
              ),
              ####################  Vaccination ################# 
              tabPanel("Vac",  
                    h5("Vaccination"),
                    sliderInput("Vbase", "Uniform Baseline Vaccination Level (%)",
                                min = 0,  max = 100,
                                value = Vbase_BasicStart,
                                step = 1
                    ),
                    numericInput("Vrate", "Rate of vaccination (/day)", 0, #5000,
                                   min = 0, max = 20000, step=500
                    ),
                      conditionalPanel(
                        condition = "input.Vrate > 0",
                        sliderInput("Vdelay", "Start delay for vaccination after epidemic detection (days)",
                                    min = 0,  max = 280,
                                    value = Vdelay_BasicStart,
                                    step = 1
                        ),
                        sliderInput("vacrad",
                                    "Radius of vaccination ring for vaccination around cases (km)",
                                    min = 3,
                                    max = 30,
                                    value = vacrad_BasicStart #4#2    #pixel size is 3 nautical miles (unaggregated)
                        ),
                        sliderInput("vacprog", "Time to build immunity after vaccination (days)",
                                    min = 1,  
                                    max = 14,
                                    value = vac_TimeToImmunity_Rind
                        ),
                        plotOutput("vaccProg", height = "300px", width = "300px")
                    )
              ),
              ####################  Disease Progression ################# 
              tabPanel("Dis. Prog.",
                          sliderInput("beta", #"Transmissibiliy", 
                                   "Beta (as in an SIR model) (/day)",
                                   min = 0.1,#1,                                                                                                                             
                                   max = beta_max, #100,                                                                                                                                       
                                   value = beta_Rind, #0.48 #0.9 #0.26#0.35 #32
                                   step = 0.01
                          ), 
                          sliderInput("IncTime", "Non-contagious Incubation period (days)",                                                                                                                        
                                      min = 0.25,                                                                                                                                  
                                      max = 10,                                                                                                                                       
                                      value = 1/kEI_Rind
                          ),
                          # numericInput("kIH", "Rate from I to H (1/days)", val = kIH_Rind,
                          #           min = 0, max = 20),
                          # numericInput("kIR", "Rate from I to R (1/days)", val = kIR_Rind,
                          #           min = 0, max = 20),
                       
                          fluidRow(
                            column(6,     
                                numericInput("kIH", "k(I->H) (1/day)", val = kIH_Rind, min = 0, max = 20),
                                numericInput("kHD", "k(H->D) (1/day)", val = kHD_Rind, min = 0, max = 20)
                                # numericInput("kHtD", "k(Ht->D) (1/day)", val = kHtD_Rind, min = 0, max = 20)
                                # numericInput("kIH", "k(I->H) (1/day)", val = kIH_H7N9, min = 0, max = 20),
                                # numericInput("kHD", "k(H->D) (1/day)", val = kHD_H7N9, min = 0, max = 20),
                                # numericInput("kHtD", "k(Ht->D) (1/day)", val = kHtD_H7N9, min = 0, max = 20)
                                
                            ),
                            column(6, 
                                numericInput("kIR", "k(I->R) (1/day)", val = kIR_Rind, min = 0, max = 20),
                                numericInput("kHR", "k(H->R) (1/day)", val = kHR_Rind, min = 0, max = 20)
                                # numericInput("kHtR", "k(Ht->R) (1/day)", val = kHtR_Rind, min = 0, max = 20)
                                # numericInput("kIR", "k(I->R) (1/day)", val = kIR_H7N9, min = 0, max = 20),
                                # numericInput("kHR", "k(H->R) (1/day)", val = kHR_H7N9, min = 0, max = 20),
                                # numericInput("kHtR", "k(Ht->R) (1/day)", val = kHtR_H7N9, min = 0, max = 20)
                            )
                          ),  
                       
                          # sliderInput("mortRate", "Fraction of seriously ill that die (w/o treatment)",                                                                                                                        
                          #             min = 0,                                                                                                                                  
                          #             max = 1,
                          #             value = kHD_Rind/(kHD_Rind+kHR_Rind)
                          # ),
                          
                          sliderInput("treated", "Fraction of seriously ill that are treated",                                                                                                                        
                                   min = 0,                                                                                                                                  
                                   max = 1,
                                   value = Rr_HHt_Rind/(Rr_HHt_Rind+1)
                          ),
                       
                          #tableOutput("rates"),
                          plotOutput("diseaseProg", height = "300px", width = "300px")
              )
        ) # End of tabset       
      ) #end of far right column
    ) #end of fluid row
  ) #end of fluid page                                                                                                                                                                                                                                                                            
) #End of ShinyUI                                                                                                                                                                                                             
