

options(shiny.maxRequestSize=50*1024^2)

# Loads the code into environment
source("R/LICAR_functions.R")

ui <- fluidPage(
  titlePanel("LICAR"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      # File upload
      fileInput("file1", "1) Choose CSV File", multiple = TRUE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      
      # Display all or only the head of the table
      radioButtons("disp", "Data preview", choices = c(Head = "head", All = "all"), selected = "head"),
      
      tags$hr(),
      
      # Choose the moiety from which the product ion in the MRM transition originates 
      selectInput("productIon", label = "2) Choose product ion type", choices = list("Head Group", "FA", "LCB")),
      
      # Choose the lipid class and corresponding MRM pattern
      uiOutput("lipidClassSelector"),
      
      #Notes
      HTML("<b>Notes</b>:<br>
           <br>
           - R packages needed: 'enviPat'<br>
           <br>
           - When correcting the intensities of very low abundance peaks, measurement variations may result in negative corrected values. LICAR replaces these negative values by 0 and does not proceed with further correction.<br>
           <br>
           <b>"),
      
      tags$a(href="https://www.github.com/SLINGhub/LICAR/manual/LICAR_Manual.pdf", "Download the user manual.")

    ),
    
    
    mainPanel(
      
      img(src='SLING.png', width="40%"),
      
      tags$hr(),
      

      HTML("<br/><br/>
           <b>Please verify the data below were uploaded correctly (Data in the same lipid class).
              If yes, click the 'Isotopic correction' button below to perform isotopic correction.
              If no, please check and re-upload your raw data file.</b>
           <br/><br/>"),
      
      #Chech whether data belong to one class
      textOutput('classError'),
      
      # Output: Raw data file
      tableOutput("contents"),
      
      # Click to perform isotopic correction
      inputPanel(actionButton("goButton", "Isotopic correction")),
      
      # Output: Isotopic corrected data file
      tableOutput("results"),
      
      HTML("<br/><br/>
           <b>Please click the 'Download Results' button below to download the results of isotopic correction.</b>
           <br/><br/>"),
      
      # Click to download the isotopic corrected data file
      inputPanel(downloadButton("downloadResults","Download Results")),
      
      HTML("<br/><br/>
           <b>Please click the 'Download Ratios' button below to calculate and download the ratio of peak areas between after and before isotopic correction.</b>
           <br/><br/>"),
      
      # Click to calculate and download the ratio of peak areas between after and before isotopic correction
      inputPanel(downloadButton("downloadRatios","Download Ratios"))
      
      )
    )
  )

server <- function(input, output) {
  
  # Fetch the appropriate data object
  rawData <- reactive({
    req(input$file1)
    df <- read.csv(input$file1$datapath, header=TRUE, row.names = 1, check.names = FALSE)
    colnames(df)[1] <- "Precursor"
    colnames(df)[2] <- "Product"
    df <- df[order(df$Precursor, df$Product), ] #Sort the df
  })
  
  
  # Check the lipid class is unique or not
  output$classError <- renderText({
    class_name <- sub( " .*", "", rownames(rawData()) )
    if( length(table(class_name)) > 1 ) {
      stop(paste("Lipid class is not unique, please check the data!"))
    }})
    
 
  # Display the selected raw data table
  output$contents <- renderTable({
    if(input$disp == "head") {return(head(rawData()))}
    else {return(rawData())}
  }, rownames = 1)
  
  
  # Display the lipid class and MRM pattern depending on the selection of product ion origin
  output$lipidClassSelector <- renderUI({
    
    #Choose the choices1 based on input$productIon
    if(input$productIon %in% "Head Group") {
      choices1 = list( "LPC (Pos) Pro=184.1" = "LPC",
                       "LPC-O (Pos) Pro=104.1" = "LPCO",
                       "LPC-O (Pos, qualifier) Pro=184.1" = "LPCOql",
                       "PC (Pos) Pro=184.1" = "PC",
                       "SM (Pos) Pro=184.1" = "SM",
                       "LPE (Pos) Pre-Pro=141" = "LPE",
                       "PE (Pos) Pre-Pro=141" = "PE",
                       "PG (Pos) Pre-Pro=189" = "PG",
                       "PI (Pos) Pre-Pro=277" = "PI",
                       "PS (Pos) Pre-Pro=185" = "PS",
                       "LPE (Neg) Pro=196" = "LPENHG",
                       "PE (Neg) Pro=196" = "PENHG",
                       "PI (Neg) Pro=241" = "PINHG",
                       "PG (Neg) Pro=153" = "PGNHG",
                       "PS (Neg) Pre-Pro=87" = "PSNHG")
    } else if(input$productIon %in% "FA") {
      choices1 = list( "PE-P (Pos) FA" = "PEP",
                       "PC (Neg) FA" = "PCNFA",
                       "PE (Neg) FA" = "PENFA",
                       "PG (Neg) FA" = "PGNFA",
                       "PI (Neg) FA" = "PINFA",
                       "PS (Neg) FA" = "PSNFA")
    } else if(input$productIon %in% "LCB") {
      choices1 = list(  "dhCer (Pos) SphB" = "dhCer",
                        "Cer (Pos) SphB-H2O" = "Cer",
                        "Hex1Cer (Pos) SphB-H2O" = "Hex1Cer",
                        "Hex2Cer (Pos) SphB-H2O" = "Hex2Cer")
    }
    
    #Choose the choices2 based on the name of the rawData()
    class_names <- sub( " .*", "", rownames(rawData()) )
    class_name <- names(table(class_names))
    numChar <- nchar(class_name)
    
    choices2 <- list()
    if( numChar == 2) {
      k=1
      name_temp <- substr(names(table(class_name)), 1, 2)
      for( i in 1:length(choices1) )
      {
        if( name_temp == substr(choices1[[i]], 1, 2) ) {
          choices2[[k]] <- choices1[[i]]
          names(choices2)[k] <- names(choices1)[i]
          k = k + 1
        }
      }
      
      if( length(choices2) == 0 ) {
        stop(paste("Lipid class is wrong, please choose the class again!"))
      }
    } else {
      k=1
      name_temp <- substr(names(table(class_name)), 1, 3)
      for( i in 1:length(choices1) )
      {
        if( name_temp == substr(choices1[[i]], 1, 3) ) {
          choices2[[k]] <- choices1[[i]]
          names(choices2)[k] <- names(choices1)[i]
          k = k + 1
        }
      }
      
      if( length(choices2) == 0 ) {
        stop(paste("Lipid class is wrong, please choose the class again!"))
      }
    }
    
    selectInput("lipidClass",label = "3) Choose lipid class and MRM pattern:",
                choices = choices2)
    
    })
  
  
  # Perform isotopic correction with "isoCorrection" from "methods.R"
  correctedResults <- reactive(
    {
      isoCorrection <- isoCorrect(rawData(), lipidClass = input$lipidClass, lipidGroup = input$productIon)
      }
    )
  
  # Display the isotopic corrected data table
  observeEvent(input$goButton,{
    output$results <- renderTable({
      if(input$disp == "head") {return(head(correctedResults()))}
      else {return(correctedResults())}
      }, rownames = 1)
    })
  
  # Download the isotopic corrected data file
  output$downloadResults <- downloadHandler(
    filename = function() {
      fpath <- input$file1$datapath
      paste("Results of isotopic correction for ", input$file1$name, sep="")
      },
    content = function(file) {write.csv(correctedResults(), file)}
    )
  
  # Calculate and download the ratio of peak areas between before and after isotopic correction
  output$downloadRatios <- downloadHandler(
    filename = function() {
      fpath <- input$file1$datapath
      paste("Ratios after_before isotopic correction for ", input$file1$name, sep="")
      },
    content = function(file) {
      ratio <- correctedResults() / rawData()
      ratio[, c("Precursor", "Product")] <- rawData()[, c("Precursor", "Product")]
      write.csv(ratio, file)
      }
    )
  
  }

shinyApp(ui = ui, server = server)
