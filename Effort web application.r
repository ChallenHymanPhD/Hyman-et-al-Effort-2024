#------------------------------------------------------------------------------#
##################################### Libraries ################################
#------------------------------------------------------------------------------#
## Syntax packages
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))

## Visualization packages
suppressMessages(library(ggplot2))
suppressMessages(library(kableExtra))
suppressMessages(library(knitr))
suppressMessages(library(jpeg))
suppressMessages(library(magrittr))
suppressMessages(library(xtable))
options(width = 80)

### Set ggplot theme for plotting
theme_set(theme_bw(base_family = 'serif')) 

suppressMessages(library(grid))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))

## Modeling
suppressMessages(library(rstan))
suppressMessages(library(brms))
suppressMessages(library(forecast))

## Application
suppressMessages(library(shiny))

options(scipen = 999)

## Set working directory to load models
setwd("C:\\Users\\ichal\\Documents\\Hyman gag grouper models\\Reef Fish Effort Shiny Application")

## Load Data (aggregated from Hyman et al (in prep))
Effort <- read.csv("Effort_data.csv")

## Load Models
Gag_effort <- read.csv("Gag_Effort.csv") ## Gag
RS_effort  <- read.csv("RS_Effort.csv")  ## Red Grouper
RG_effort  <- read.csv("RG_Effort.csv")  ## Red Snapper
GS_effort  <- read.csv("GS_Effort.csv")  ## Gray Snapper

#------------------------------------------------------------------------------#
##################################### Functions ################################
#------------------------------------------------------------------------------#
## Calculate seasonal lengths
Season_Month_Fraction <- function(
  Year = 2023,
  Season = c("January 01", "June 15"),
  Additional = NULL
){
  if(is.null(Season)){stop("Please specify a season start and end date")}
  if(is.null(Year)){stop("Please specify a year")}
  
  ## This creates a Data Frame of all days in each month of the specified year
  All_Dates_Year <- as.Date(paste(c("January 01", "December 31"), Year), format = "%B %d %Y")
  All_Dates_Year <- data.frame(Date = seq.Date(All_Dates_Year[1], All_Dates_Year[2], by = 'day'))
  All_Dates_Year$Month <- month(All_Dates_Year$Date)
  Month_days <- table(All_Dates_Year$Month)%>%as.data.frame()
  
  if(is.null(Additional)){
    ## Create vector of season dates
    Dates <- as.Date(paste(Season, Year), format = "%B %d %Y")
    Date_vector <- data.frame(Date = seq.Date(Dates[1], Dates[2], by = 'day'))
    Date_vector$Month <- month(Date_vector$Date)
    Month_days_season <- table(Date_vector$Month)%>%as.data.frame()
    Month_days$Season <- Month_days_season$Freq[match(Month_days$Var1,Month_days_season$Var1)]
    Month_days$Season[which(is.na(Month_days$Season))] <- 0
    Month_days$Fraction <- Month_days$Season/Month_days$Freq
  } else {
    Dates_1 <- as.Date(paste(Season, Year), format = "%B %d %Y")
    Date_vector_1 <- data.frame(Date = seq.Date(Dates_1[1], Dates_1[2], by = 'day'))
    Dates_2 <- as.Date(paste(Additional, Year), format = "%B %d %Y")
    Date_vector_2 <- data.frame(Date = seq.Date(Dates_2[1], Dates_2[2], by = 'day'))
    Date_vector <- rbind(Date_vector_1, Date_vector_2)
    Date_vector$Month <- month(Date_vector$Date)
    Month_days_season <- table(Date_vector$Month)%>%as.data.frame()
    Month_days$Season <- Month_days_season$Freq[match(Month_days$Var1,Month_days_season$Var1)]
    Month_days$Season[which(is.na(Month_days$Season))] <- 0
    Month_days$Fraction <- Month_days$Season/Month_days$Freq
  }
  Fraction_of_months_open <- Month_days
  Season_length <- nrow(Date_vector)
  return(list(Fraction_of_months_open, Season_length))
}

## Wrapper function for conditional expected trips
BRM_predictions <- function(Model = NULL, Data = NULL, CI = NULL, nsamples = 100){
  probs <- c(0+(1-CI)/2,0.5, 1-(1-CI)/2) ## CI based on user-defined value
  post <- Model ## extract posterior samples for each model
  Preds <- matrix(NA, ncol = 12, nrow = nsamples) ## create data frame to store monthly predictions for all scenarios
  Pred_names <- colnames(post)%>%gsub("b_", "", .)%>%gsub("hu_", "", .) ## All predictor names to match to data
  nonzero_names <- colnames(post)[grep("b_", colnames(post))][-grep("b_hu_", colnames(post))] ## Nonzero predictor names
  hurdle_names <- colnames(post)[grep("b_hu_", colnames(post))] ## Hurdle predictor names
  
  ## Run model to obtain predictions for each scan
  for (i in 1:nsamples){
    betas <- post[i, which(colnames(post) %in% nonzero_names)]%>%as.vector()
    gammas <- post[i, which(colnames(post) %in% hurdle_names)]%>%as.vector()
    Model_data <- cbind(1,1,Data[,na.omit(match(Pred_names, colnames(Data)))])
    nonzero_preds <- as.matrix(cbind(1,Data[match(gsub("b_", "", nonzero_names[-1]), colnames(Data))]))
    hurdle_preds <- as.matrix(cbind(1,Data[match(gsub("b_hu_", "", hurdle_names[-1]), colnames(Data))]))
    nonzero <- nonzero_preds%*%t(betas)
    theta <- plogis(hurdle_preds%*%t(gammas))
    ## Expected values with lognormal correction
    Preds[i,] <- (1-theta)*exp(nonzero + (post$sigma[i]^2)/2)
  }
  Timeseries_output <- apply(Preds, 2, function(x){quantile(x, probs)})%>%t()
  Plot_output <- apply(Preds, 2, function(x){quantile(x, probs)})%>%t()
  colnames(Plot_output) <- c("Low", "Med", "Up")
  return(list(Plot_output))
}

## Function to develop plots of predicted effort based on counterfactual scenarios
Effort_test <- function(
  Gag_effort = Gag_effort,
  RS_effort = RS_effort,
  RG_effort = RG_effort,
  GS_effort = GS_effort,
  Gag_Season = c("June 01", "December 31"),
  RG_Season = c("June 01", "September 29"),
  RS_Season = c("June 01", "August 15"),
  Gag_Additional = NULL,
  RG_Additional = NULL,
  RS_Additional = NULL,
  Gag_Abundance = "Current",
  RG_Abundance = "Current",
  RS_Abundance = "Current",
  GS_Abundance = "Current",
  Conditions = "Normal",
  Dataset = Effort,
  Year = 2023,
  Facet = F
){
  ## Set Fishing Conditions
  if(Conditions == "Normal"){Fishable <- tapply(Dataset$Fishable, Dataset$Month, mean)}
  if(Conditions == "Good"){Fishable <- tapply(Dataset$Fishable, Dataset$Month, function(x){quantile(x,0.9)})}
  if(Conditions == "Poor"){Fishable <- tapply(Dataset$Fishable, Dataset$Month, function(x){quantile(x,0.1)})}
  
  ## Set population
  ### Gag
  if(Gag_Abundance == "Current"){A_Gag <- Effort$A_Gag[which(Effort$Year == 2022)]}
  if(Gag_Abundance == "High"){A_Gag <- quantile(Effort$A_Gag, 0.9)}
  if(Gag_Abundance == "Low"){A_Gag <- quantile(Effort$A_Gag, 0.1)}
  
  ### Red Grouper
  if(RG_Abundance == "Current"){A_RG <- Effort$A_RG[which(Effort$Year == 2022)]}
  if(RG_Abundance == "High"){A_RG <- quantile(Effort$A_RG, 0.9)}
  if(RG_Abundance == "Low"){A_RG <- quantile(Effort$A_RG, 0.1)}
  
  ### Red Snapper
  if(RS_Abundance == "Current"){A_RS <- Effort$A_RS[which(Effort$Year == 2022)]}
  if(RS_Abundance == "High"){A_RS <- quantile(Effort$A_RS, 0.9)}
  if(RS_Abundance == "Low"){A_RS <- quantile(Effort$A_RS, 0.1)}
  
  ### Gray Snapper
  if(GS_Abundance == "Current"){A_GS <- Effort$A_GS[which(Effort$Year == 2022)]}
  if(GS_Abundance == "High"){A_GS <- quantile(Effort$A_GS, 0.9)}
  if(GS_Abundance == "Low"){A_GS <- quantile(Effort$A_GS, 0.1)}
  
  
  ## Set management conditions
  ### Gag
  M_Gag <- Season_Month_Fraction(Year = Year,
                                 Season = Gag_Season,
                                 Additional = Gag_Additional)[[1]]$Fraction
  S_Gag <- Season_Month_Fraction(Year = Year,
                                 Season = Gag_Season,
                                 Additional = Gag_Additional)[[2]]/100
  
  ### Red Snapper
  M_RS <- Season_Month_Fraction(Year = Year,
                                Season = RS_Season,
                                Additional = RS_Additional)[[1]]$Fraction
  S_RS <- Season_Month_Fraction(Year = Year,
                                Season = RS_Season,
                                Additional = RS_Additional)[[2]]/100
  
  ### Red Grouper
  M_RG <- Season_Month_Fraction(Year = Year,
                                Season = RG_Season,
                                Additional = RG_Additional)[[1]]$Fraction
  S_RG <- Season_Month_Fraction(Year = Year,
                                Season = RG_Season,
                                Additional = RG_Additional)[[2]]/100
  
  #### NOTE: Gray Snapper does not have a season
  Model_Data_frame <- data.frame(Time = rep(1:12),
                                 Fishable = Fishable,
                                 A_Gag = A_Gag,
                                 M_Gag = M_Gag,
                                 S_Gag = S_Gag,
                                 M_RS = M_RS,
                                 M_RG = M_RG,
                                 A_RG = A_RG,
                                 A_GS = A_GS,
                                 S_RG = S_RG,
                                 S_RS = S_RS,
                                 A_RS = A_RS)
  
  ## Harmonic regression terms
  per <- 12 ## period 1 is annual (12 months in a year)
  per2 <- 6 ## period 2 is semi-annual (6 months in a year)
  
  ## First harmonic
  Model_Data_frame$sin1 <- sin(2*pi/per*Model_Data_frame$Time)
  Model_Data_frame$cos1 <- cos(2*pi/per*Model_Data_frame$Time)
  
  ## Second harmonic
  Model_Data_frame$sin2 <- sin(2*pi/per2*Model_Data_frame$Time)
  Model_Data_frame$cos2 <- cos(2*pi/per2*Model_Data_frame$Time)
  
  ## Reformat management terms
  Model_Data_frame$Open <- ceiling(Model_Data_frame$M_RS)
  Model_Data_frame$logS_Gag <- log(Model_Data_frame$S_Gag)
  Model_Data_frame$logS_RG <- log(Model_Data_frame$S_RG)
  Model_Data_frame$logS_RS <- log(Model_Data_frame$S_RS)
  
  ## Obtain effort predictions with 80% confidence interval
  Gag_preds <- BRM_predictions(Model = Gag_effort, Data = Model_Data_frame, CI = 0.8)
  Model_Data_frame_Model_Gag <- Gag_preds%>%as.data.frame()
  Model_Data_frame$Estimate <- Model_Data_frame_Model_Gag$Med
  Model_Data_frame$Low <- Model_Data_frame_Model_Gag$Low
  Model_Data_frame$Up <- Model_Data_frame_Model_Gag$Up
  
  RG_preds <- BRM_predictions(Model = RG_effort, Data = Model_Data_frame, CI = 0.8)
  Model_Data_frame_RG <- Model_Data_frame
  Model_Data_frame_Model_RG <- RG_preds%>%as.data.frame()
  Model_Data_frame_RG$Estimate <- Model_Data_frame_Model_RG$Med
  Model_Data_frame_RG$Low <- Model_Data_frame_Model_RG$Low
  Model_Data_frame_RG$Up <- Model_Data_frame_Model_RG$Up
  
  RS_preds <- BRM_predictions(Model = RS_effort, Data = Model_Data_frame, CI = 0.8)
  Model_Data_frame_RS <- Model_Data_frame
  Model_Data_frame_Model_RS <- RS_preds%>%as.data.frame()
  Model_Data_frame_RS$Estimate <- Model_Data_frame_Model_RS$Med
  Model_Data_frame_RS$Low <- Model_Data_frame_Model_RS$Low
  Model_Data_frame_RS$Up <- Model_Data_frame_Model_RS$Up
  
  GS_preds <- BRM_predictions(Model = GS_effort, Data = Model_Data_frame, CI = 0.8)
  Model_Data_frame_GS <- Model_Data_frame
  Model_Data_frame_Model_GS <- GS_preds%>%as.data.frame()
  Model_Data_frame_GS$Estimate <- Model_Data_frame_Model_GS$Med
  Model_Data_frame_GS$Low <- Model_Data_frame_Model_GS$Low
  Model_Data_frame_GS$Up <- Model_Data_frame_Model_GS$Up
  
  ## Format datasets for plotting
  Model_Data_frame$Species <- "Gag grouper"
  Model_Data_frame_RG$Species <- "Red grouper"
  Model_Data_frame_RS$Species <- "Red snapper"
  Model_Data_frame_GS$Species <- "Gray snapper"
  
  Model_Data_frame$Month <- "Jan"
  Model_Data_frame$Month[which(Model_Data_frame$Time==2)] <- "Feb"
  Model_Data_frame$Month[which(Model_Data_frame$Time==3)] <- "Mar"
  Model_Data_frame$Month[which(Model_Data_frame$Time==4)] <- "Apr"
  Model_Data_frame$Month[which(Model_Data_frame$Time==5)] <- "May"
  Model_Data_frame$Month[which(Model_Data_frame$Time==6)] <- "Jun"
  Model_Data_frame$Month[which(Model_Data_frame$Time==7)] <- "Jul"
  Model_Data_frame$Month[which(Model_Data_frame$Time==8)] <- "Aug"
  Model_Data_frame$Month[which(Model_Data_frame$Time==9)] <- "Sep"
  Model_Data_frame$Month[which(Model_Data_frame$Time==10)] <- "Oct"
  Model_Data_frame$Month[which(Model_Data_frame$Time==11)] <- "Nov"
  Model_Data_frame$Month[which(Model_Data_frame$Time==12)] <- "Dec"
  Model_Data_frame_RG$Month <- Model_Data_frame_GS$Month <- Model_Data_frame_RS$Month <- Model_Data_frame$Month
  
  
  ## Generate plot outputs for each species
  Gag_plot <- ggplot(Model_Data_frame[which(Model_Data_frame$Species == "Gag grouper"),])+
      geom_line(aes(x = Time, y = Estimate*0.001), lwd = 1, show.legend = F, col = "#9E0142")+
      geom_ribbon(aes(x = Time, ymin = Low*0.001, ymax = Up*0.001), alpha = 0.5, show.legend = F, fill = "#9E0142")+
      labs(y = "Effort (millions of trips)", x = "Month")+
      scale_x_continuous(breaks = 1:12, labels = unique(Model_Data_frame$Month))+
      theme(axis.text = element_text(size = 14, color = 1, angle = 0),
            axis.text.x = element_text(size = 14, color = 1, angle = 0),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))+coord_cartesian(ylim = c(0,400))+facet_wrap(~Species)
    RG_plot <- ggplot(Model_Data_frame_RG)+
      geom_line(aes(x = Time, y = Estimate*0.001), lwd = 1, show.legend = F, col = "#798A16")+
      geom_ribbon(aes(x = Time, ymin = Low*0.001, ymax = Up*0.001), alpha = 0.5, show.legend = F, fill = "#798A16")+
      labs(y = "Effort (millions of trips)", x = "Month")+
      scale_x_continuous(breaks = 1:12, labels = unique(Model_Data_frame$Month))+
      theme(axis.text = element_text(size = 14, color = 1, angle = 0),
            axis.text.x = element_text(size = 14, color = 1, angle = 0),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))+coord_cartesian(ylim = c(0,300))+facet_wrap(~Species)
    RS_plot <- ggplot(Model_Data_frame_RS)+
      geom_line(aes(x = Time, y = Estimate*0.001), lwd = 1, show.legend = F, col = "#4575B4")+
      geom_ribbon(aes(x = Time, ymin = Low*0.001, ymax = Up*0.001), alpha = 0.5, show.legend = F, fill = "#4575B4")+
      labs(y = "Effort (millions of trips)", x = "Month")+
      scale_x_continuous(breaks = 1:12, labels = unique(Model_Data_frame$Month))+
      theme(axis.text = element_text(size = 14, color = 1, angle = 0),
            axis.text.x = element_text(size = 14, color = 1, angle = 0),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))+coord_cartesian(ylim = c(0,550))+facet_wrap(~Species)
    GS_plot <- ggplot(Model_Data_frame_GS)+
      geom_line(aes(x = Time, y = Estimate*0.001), lwd = 1, show.legend = F, col = "forestgreen")+
      geom_ribbon(aes(x = Time, ymin = Low*0.001, ymax = Up*0.001), alpha = 0.5, show.legend = F, fill = "forestgreen")+
      labs(y = "Effort (millions of trips)", x = "Month")+
      scale_x_continuous(breaks = 1:12, labels = unique(Model_Data_frame$Month))+
      theme(axis.text = element_text(size = 14, color = 1, angle = 0),
            axis.text.x = element_text(size = 14, color = 1, angle = 0),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5, size = 20))+coord_cartesian(ylim = c(50,300))+facet_wrap(~Species)
    return(list(Gag_plot, RG_plot, RS_plot, GS_plot))
}



# Define UI ----
ui <- fluidPage(
  titlePanel(h1("Reef Fish Effort Model", align = "center")),
  fluidRow(
    column(2,
           selectInput("Gray_snapper_population", 
                       label = "Select your gray snapper population size",
                       choices = c("High", "Current", "Low"),
                       selected = "Current"),
           
           selectInput("Fishing_conditions", 
                       label = "Select fishing conditions",
                       choices = c("Good", "Normal", "Poor"),
                       selected = "Normal"),
           mainPanel(
             p("Gray snapper has never has a federal recreational season. As a result, gray snapper recreational trips per month are only affected by the management of other species, 
               as well as fishing conditions and gray snapper abundance"))
    ),
    column(7, 
           plotOutput("GS")),
    column(3,style = "height:350px;background-color:#4575B4;color: #fff",
           mainPanel(width = 12,
             h2("About the data"),
             p("We estimated monthly trips (our proxy for angler effort) using survey data from the National Marine Fisheries Service Marine Recreational Information Program (MRIP). 
               A main component of MRIP is the Access Point Angler Intercept Survey (APAIS), which is used to collect catch-per-trip data from anglers fishing from shore, private boats, 
               and for-hire vessels based on in-person dockside interviews of anglers recently returning from a fishing trip. APAIS also collects information on location, distance 
               from shore, primary and secondary target species, and other pertinent information which may be used to subset databases for relevant trips of interest and/or as predictors 
               modeling angler behavior"))
    )
  ),
  
  fluidRow(
    column(2,
           selectInput("Gag_grouper_population", 
                       label = "Select your gag grouper population size",
                       choices = c("High", "Current", "Low"),
                       selected = "Current"),
           dateRangeInput('Gag_Season',
                          label = 'Set the gag grouper season',
                          start = as.Date("2023-06-01"), end = as.Date("2023-12-31")
           ),
           selectInput("Gag_reopen", 
                       label = "Do you want the Gag grouper season to reopen later in the year?",
                       choices = c("yes", "no"),
                       selected = "no"),
           dateRangeInput('Gag_Additional',
                          label = 'Set the additional gag grouper dates (if reopened)',
                          start = as.Date("2023-09-01"), end = as.Date("2023-12-31")
           ),
    ),
    column(7, 
           plotOutput("Gag")),
    column(3,
           mainPanel(width = 12,
             h2("About the model"),
             p("Hurdle-lognormal (HLN) generalized linear models based on APAIS data were used to generate these graphs. 
               The HLN models first predict whether any anglers will embark on any trip for a given species at all. This 
               is based on whether the season for that species is open or closed. If at least one angler trip is expected, 
               the model then calculated the most likely number of trips to be taken given the user-defined conditions. 
               When you set the conditions, the model uses information from past observations to make guesses about what 
               might occur in the future."))
           )
    ),
  fluidRow(
    column(2,
           selectInput("Red_grouper_population", 
                       label = "Select your red grouper population size",
                       choices = c("High", "Current", "Low"),
                       selected = "Current"),
           dateRangeInput('RG_Season',
                          label = 'Set the red grouper season',
                          start = as.Date("2023-06-01"), end = as.Date("2023-09-29")
           ),
           selectInput("RG_reopen", 
                       label = "Do you want the red grouper season to reopen later in the year?",
                       choices = c("yes", "no"),
                       selected = "no"),
           dateRangeInput('RG_Additional',
                          label = 'Set the additional red grouper dates (if reopened)',
                          start = as.Date("2023-10-01"), end = as.Date("2023-12-31")
           ),
    ),
    column(7, 
           plotOutput("RG")),
    column(3,style = "height:350px;background-color:#4575B4;color: #fff",
           mainPanel(width = 12,
             h2("About the plots"),
             p("These graphs depict predictions of the number monthly directed trips for each species given user-defined 'what if' scenarios. 
             For each plot, the solid line is the model's best guess for what total monthly effort will be in a given month. However, there is 
             inherent uncertainty in the exact number of trips which may be expected. The colored banded above and below the solid line convey 
             that uncertainty. Roughly 80% of simulated model results fell within this band, and so we can assume with 80% confidence that the \
             true expected value, given the user-defined conditions, falls within this range."))
    )
  ),
  fluidRow(
    column(2,
           selectInput("Red_snapper_population", 
                       label = "Select your red snapper population size",
                       choices = c("High", "Current", "Low"),
                       selected = "Current"),
           dateRangeInput('RS_Season',
                          label = 'Set the red snapper season',
                          start = as.Date("2023-06-01"), end = as.Date("2023-08-31")
           ),
           selectInput("RS_reopen", 
                       label = "Do you want the red snapper season to reopen later in the year?",
                       choices = c("yes", "no"),
                       selected = "no"),
           dateRangeInput('RS_Additional',
                          label = 'SSet the additional red snapper dates (if reopened)',
                          start = as.Date("2023-09-01"), end = as.Date("2023-10-31")
           ),
    ),
    column(7, 
           plotOutput("RS"))
  )
  
  
)

# Server logic ----
server <- function(input, output) {
  Plotlist <- reactive({
    Gag_Abundance <- input$Gag_grouper_population
    RG_Abundance <- input$Red_grouper_population
    RS_Abundance <- input$Red_snapper_population
    GS_Abundance <- input$Gray_snapper_population
    Conditions <- input$Fishing_conditions
    Gag_Season <- as.Date(input$Gag_Season)
    RG_Season <- as.Date(input$RG_Season)
    RS_Season <- as.Date(input$RS_Season)
    
    Gag_start_1 <- paste(month(Gag_Season[1], label = T, abbr = F), ifelse(nchar(paste(day(Gag_Season[1])))<2, paste0("0", day(Gag_Season[1])), paste0(day(Gag_Season[1]))))
    Gag_end_1<- paste(month(Gag_Season[2], label = T, abbr = F), ifelse(nchar(paste(day(Gag_Season[2])))<2, paste0("0", day(Gag_Season[2])), paste0(day(Gag_Season[2]))))
    Gag_Season <- c(Gag_start_1, Gag_end_1)
    #if(Gag_end_1>=Gag_start_1){stop("The gag grouper start date cannot be later than the end date!")}
    
    RG_start_1 <- paste(month(RG_Season[1], label = T, abbr = F), ifelse(nchar(paste(day(RG_Season[1])))<2, paste0("0", day(RG_Season[1])), paste0(day(RG_Season[1]))))
    RG_end_1<- paste(month(RG_Season[2], label = T, abbr = F), ifelse(nchar(paste(day(RG_Season[2])))<2, paste0("0", day(RG_Season[2])), paste0(day(RG_Season[2]))))
    RG_Season <- c(RG_start_1, RG_end_1)
    #if(RG_end_1>=RG_start_1){stop("The red grouper start date cannot be later than the end date!")}                  
    
    RS_start_1 <- paste(month(RS_Season[1], label = T, abbr = F), ifelse(nchar(paste(day(RS_Season[1])))<2, paste0("0", day(RS_Season[1])), paste0(day(RS_Season[1]))))
    RS_end_1<- paste(month(RS_Season[2], label = T, abbr = F), ifelse(nchar(paste(day(RS_Season[2])))<2, paste0("0", day(RS_Season[2])), paste0(day(RS_Season[2]))))
    RS_Season <- c(RS_start_1, RS_end_1)
    #if(RS_end_1>=RS_start_1){stop("The red snapper  start date cannot be later than the end date!")}
    
    Gag_Additional <- NULL
    RG_Additional <- NULL
    RS_Additional <- NULL
    
    
    if(input$Gag_reopen == "yes"){
      Gag_Additional <- as.Date(input$Gag_Additional)
      Gag_start_2 <- paste(month(Gag_Additional[1], label = T, abbr = F), ifelse(nchar(paste(day(Gag_Additional[1])))<2, paste0("0", day(Gag_Additional[1])), paste0(day(Gag_Additional[1]))))
      Gag_end_2<- paste(month(Gag_Additional[2], label = T, abbr = F), ifelse(nchar(paste(day(Gag_Additional[2])))<2, paste0("0", day(Gag_Additional[2])), paste0(day(Gag_Additional[2]))))
      Gag_Additional <- c(Gag_start_2, Gag_end_2)
      #if(Gag_end_2>=Gag_start_2){stop("The additional gag grouper start date cannot be later than the end date!")}
    }
    
    if(input$RG_reopen == "yes"){
      RG_Additional <- as.Date(input$RG_Additional)
      RG_start_2 <- paste(month(RG_Additional[1], label = T, abbr = F), ifelse(nchar(paste(day(RG_Additional[1])))<2, paste0("0", day(RG_Additional[1])), paste0(day(RG_Additional[1]))))
      RG_end_2<- paste(month(RG_Additional[2], label = T, abbr = F), ifelse(nchar(paste(day(RG_Additional[2])))<2, paste0("0", day(RG_Additional[2])), paste0(day(RG_Additional[2]))))
      RG_Additional <- c(RG_start_2, RG_end_2)
      #if(RG_end_2>=RG_start_2){stop("The additional red grouper start date cannot be later than the end date!")}
    }
    
    if(input$RS_reopen == "yes"){
      RS_Additional <- as.Date(input$RS_Additional)
      RS_start_2 <- paste(month(RS_Additional[1], label = T, abbr = F), ifelse(nchar(paste(day(RS_Additional[1])))<2, paste0("0", day(RS_Additional[1])), paste0(day(RS_Additional[1]))))
      RS_end_2<- paste(month(RS_Additional[2], label = T, abbr = F), ifelse(nchar(paste(day(RS_Additional[2])))<2, paste0("0", day(RS_Additional[2])), paste0(day(RS_Additional[2]))))
      RS_Additional <- c(RS_start_2, RS_end_2)
      #if(RS_end_2>=RS_start_2){stop("The additional red snapper start date cannot be later than the end date!")}
    }
    
    Effort_test(
      Gag_effort = Gag_effort,
      RS_effort = RS_effort,
      RG_effort = RG_effort,
      GS_effort = GS_effort,
      Gag_Season = Gag_Season,
      RG_Season = RG_Season,
      RS_Season = RS_Season,
      Gag_Additional = Gag_Additional,
      RG_Additional = RG_Additional,
      RS_Additional = RS_Additional,
      Gag_Abundance = Gag_Abundance,
      RG_Abundance = RG_Abundance,
      RS_Abundance = RS_Abundance,
      GS_Abundance = GS_Abundance,
      Conditions = Conditions,
      Dataset = Effort,
      Year = 2023)
  })
  
  
  output$Gag <- renderPlot({Plotlist()[[1]]})
  
  output$RG <- renderPlot({Plotlist()[[2]]})
  
  output$RS <- renderPlot({Plotlist()[[3]]})
  
  output$GS <- renderPlot({Plotlist()[[4]]})
  
}

# Run app ----
shinyApp(ui, server)
