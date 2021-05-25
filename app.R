
# setwd("C:\\Users\\insig\\Documents\\Projects\\Internal\\P0002_DeployRShiny\\sandbox-play")

# https://gist.github.com/yihui/6091942

library(shiny)
library(openxlsx)
library(heatmap3)
library(scatterplot3d) 
library(lattice)
library(tools)


source("functions.R")

ui <- shinyUI(
    fluidPage(
	   markdown("
		Example Usable RShiny App - TMTPrepPro
		"),
		
		sidebarPanel(
		
		fileInput('datafile', 'Choose Data Excel File', accept=".xlsx"),	
		fileInput('designfile', 'Choose Design Excel File', accept=".xlsx"),
		textInput('PvalCutoff', 'P-value cutoff for significance', "0.05"),
		textInput('FCCutoff', 'Fold change cutoff for significance', "1.2"),
		textInput('DuplicateFilter', 'Remove competitor proteins from groups', 'TRUE'),
		textInput('KeepREF', 'Keep fields marked as REF in design file', "TRUE"),
		textInput('SampleLoadNorm', 'Normalisation method', "total"),
		actionButton("runscript", "Run Script")
		),
		
		mainPanel(
		 tabsetPanel(
		  tabPanel("InputData", tableOutput("inputData")),
		  tabPanel("InputDesign", tableOutput("inputDesign")),
		  tabPanel("Inputcomparisons", tableOutput("Classification")),
		  tabPanel("HeatmapDE", plotOutput("heatmapDE")),
		  tabPanel("Params", verbatimTextOutput("currentParams"))
		 )
		)
		
  )
)

server <- function(input, output) {

options(shiny.maxRequestSize=10*1024^2)


observeEvent(  input$runscript, {            
       output$currentParams <- renderText({
                  paste0("Data file = ", input$datafile$name, "\n",
				        "Design file = ", input$designfile$name, "\n",
						 "FC Cutoff - ", input$FCCutoff, "\n",
						 "Pval Cutoff - ", input$PvalCutoff, "\n",
						 "Normalisation - ", input$SampleLoadNorm, "\n",
						 "Remove duplicates from protein groups - ", input$DuplicateFilter, "\n",
						 "Keep Reference REF - ", input$KeepREF, "\n"
						 )
                })
	    }
)

				
designRead <- eventReactive(input$runscript, {

   inFile <- input$designfile
  
  if(is.null(inFile)){
    return(NULL)
    } else

    design = read.xlsx(inFile$datapath, 1)
    comparisons = read.xlsx(inFile$datapath, 2)

list(design=design, comparisons=comparisons)
  })
  
  

dataRead <- eventReactive(input$runscript, {

   inFile <- input$datafile
  
  if(is.null(inFile)){
    return(NULL)
    } else

	dat =  readWorkbook(inFile$datapath,1, startRow = 1)
	
	dat

  })
  
dataProcess <- eventReactive(input$runscript, {
 
  if( is.null(input$datafile) | is.null(input$designfile)  ){
    return(NULL)
    } else {

	design.res = designRead()
	design = design.res$design
	comparisons = design.res$comparisons
	dat = dataRead()
	DuplicateFilter = ( input$DuplicateFilter == "TRUE" )
	# BatchNorm = (input$BatchNorm == "TRUE")
	KeepREF = (input$KeepREF == "TRUE")
	PvalCutoff = as.numeric(input$PvalCutoff)
	FCCutoff = as.numeric(input$FCCutoff)
	SampleLoadNorm = input$SampleLoadNorm
	
	
	
	params = c( "datafile", "designfile", "PvalCutoff", "FCCutoff", "DuplicateFilter",
				"KeepREF", "SampleLoadNorm")
	paramValues = c(input$datafile$name, input$designfile$name, PvalCutoff, FCCutoff, DuplicateFilter,
				KeepREF, SampleLoadNorm)

	param.df = data.frame(params,paramValues)

	# res = processDataMatrix(dat, design, DuplicateFilter=TRUE)

	res = TMTPrePro(mg.all=dat, Design=design, Comparisons=comparisons, 
	 	PvalCutoff=PvalCutoff, FCCutoff=FCCutoff, 
			DuplicateFilter=DuplicateFilter, KeepREF=KeepREF, 
			SampleLoadNorm=SampleLoadNorm)	

	
	return(param.df)

  }
  })
   
  
   output$inputDesign <- renderTable({
 
	inFile <- input$designfile
	if(is.null(inFile)){
     return(NULL)
     } else
	
    designRead()$design
	
  })
  
 
   output$inputComparisons <- renderTable({
 
	inFile <- input$designfile
	if(is.null(inFile)){
     return(NULL)
     } else
	
    designRead()$comparisons
	
  })
  
  
   output$inputData <- renderTable({
 
	inFile <- input$datafile
  
	if(is.null(inFile)){
     return(NULL)
     } else
	
    head(dataRead())
	
  })
  
  
 output$heatmapDE <- renderPlot({
 
  if(is.null(input$datafile) | is.null(input$designfile) ){
      	return(plot(1, type="n", main="No data yet", axes=FALSE))
   } else 

  # dat = dataProcess()
  return(plot(iris[,2], iris[,3]))
  
	
  })
  
  
   output$Classification <- renderTable({
 
	if(is.null(input$datafile) | is.null(input$designfile) ){
      	return(NULL)
   } else 
	dataProcess()
	
  })
  

}

shinyApp(ui, server)

