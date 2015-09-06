# R Shiny app to simulate sensorgram based on user input assay and kinetic parameters
# Copyright (C) 2015 Jeff Dantzler 

library(shiny)
library(ggplot2)

ui <- fluidPage(
    
  sidebarLayout(
    sidebarPanel(
      sliderInput("dil", "dilution factor:", min = 2, max = 10, value = 2),
      sliderInput("ass", "association time:", min = 10, max = 4800, value = 240),
      sliderInput("dis", "dissociation time:", min = 10, max = 4800, value = 480),
      numericInput("kd", "dissociation constant:", 0.001, min = 0.00001, max = 0.1)
    ),
    
    mainPanel(plotOutput("sensorgramPlot"))
  )
)

server <- function(input, output){

  output$sensorgramPlot <- renderPlot({
    
    # define kinetic and assay parameters
    
    # kinetics parameters
    ka <- 100000  # association rate constant
    kd <- input$kd  		# dissociation rate constant
    Kd <- kd/ka			# affinity
    Rmax <- 100			# Rmax value
    
    # assay parameters
    dil <- input$dil		# dilution factor
    C4 <- Kd 			      # set middle (of 7) concentration to Kd
    C5 <- C4 * dil 		  # generate concentration series 
    C6 <- C4 * dil^2
    C7 <- C4 * dil^3
    C3 <- C4 / dil
    C2 <- C4 / dil^2
    C1 <- C4 / dil^3
    # create vector of assay concentrations
    concSeries <- c(C1,C2,C3,C4,C5,C6,C7)
    
    assoc <- input$ass		  # association time (seconds)
    dissoc <- input$dis		  # dissociation time (seconds)
    
    # create vectors for assoc & dissoc phases and whole assay
    assocTime <- 1:assoc
    dissocTime <- (assoc + 1):(assoc + 1 + dissoc)
    Time <- c(assocTime, dissocTime)
    
    # define funtions for calculating response
    
    # steady state response at a given concentration
    Req <- function(Conc, Rmax, Kd){
      (Conc*Rmax)/(Kd+Conc)
    }
    
    # k observed
    kobs <- function(ka, Conc, kd){
      (ka*Conc)+kd
    }
    
    # response at a given concentration and time during association
    assocResponse <- function(time){
      Req(Conc, Rmax, Kd) - Req(Conc, Rmax, Kd) * exp(-kobs(ka, Conc, kd) * time)
    }
    
    # single exponential decay of response during dissociation
    dissocResponse <- function(time){
      Ro*exp(-kd*(time - assoc))
    }
    
    # create a matrix to hold the data
    data <- matrix(nrow=length(Time), ncol=length(concSeries)+1)
    
    # push whole assay time to first column
    data[ ,1] <- Time
    
    # loop through the concentration series and populate the rest of the data matrix with responses at each concentration / timepoint
    for (value in 1:length(concSeries)) {
      # set which concentration we're working with
      Conc <- concSeries[value]
      # define dissociation starting response as response at end of association for this concentration
      Ro <- Req(Conc, Rmax, Kd) - Req(Conc, Rmax, Kd) * exp(-kobs(ka, Conc, kd)*assoc)
      # calculate response as a function of time for this concentration
      assocR <- assocResponse(assocTime)
      dissocR <- dissocResponse(dissocTime)
      Resp <- c(assocR, dissocR)
      # push these values into the data matrix
      data[ ,(value + 1)] <- Resp
    }
    
    # convert matrix to data frame to feed to ggplot2
    data <- as.data.frame(data)

    # Plot sensorgrams
    g <- ggplot(data, aes(V1))
    for (value in 1:length(concSeries))
      g <- g + geom_line(aes_string(y=colnames(data)[value+1]), colour="red")
    g <- g + xlab("Time (seconds)")
    g <- g + ylab("Response (RU)")
    g <- g + ggtitle("Sensorgram based on these input parameters")
    g
  })
}

shinyApp(ui = ui, server = server)