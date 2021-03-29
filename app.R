
# setwd("C:\\Users\\insig\\Documents\\Projects\\Internal\\P0002_DeployRShiny\\sandbox-play")


library(shiny)
library(UpSetR)
# library(knitr)

# rmdfiles <- c("RMarkdownFile.rmd")
# sapply(rmdfiles, knit, quiet = T)

ui <- shinyUI(
    fluidPage(
	   markdown("
		![](logo.JPG)
		
		Test Shiny App example

		[Insight Stats Home](https://insightstats.wixsite.com/home) 
		
		"),
		
		sidebarPanel(
		  numericInput('n', 'Number of samples', 3),
		  fileInput('file1', 'Choose CSV File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))	  
		),
		
		mainPanel(
		 tabsetPanel(
		  tabPanel("Plot", plotOutput("distPlot")),
		  tabPanel("Table", tableOutput("userTable")),
		  tabPanel("UpSet", plotOutput("upsetPlot"))
		 )
		)
		
        # includeMarkdown("RMarkdownFile.md")
		#, theme = "InsightStats.css"
  )
)

server <- function(input, output) {

 output$distPlot <- renderPlot({

  inFile <- input$file1
  
  if (is.null(inFile))
      	return(hist(rnorm(input$n)))

  dat = read.csv(inFile$datapath)
  boxplot(as.numeric(dat[,2]) ~ dat[,1])
  
	
  })
  
  
  output$upsetPlot <- renderPlot({

  inFile <- input$file1
  
  if (is.null(inFile))
      	return(hist(rnorm(input$n)))

  dat.View = read.csv(inFile$datapath)
  
  upset(dat.View, sets = colnames(dat.View)[-1], order.by = "freq", sets.bar.color = "#56B4E9", empty.intersections = "off")

	
  })
  
   output$userTable <- renderTable({
 
	inFile <- input$file1

    if (is.null(inFile))
      	return(NULL)
    
    head(read.csv(inFile$datapath))
	
  })
  
}

shinyApp(ui, server)

