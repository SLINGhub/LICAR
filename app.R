

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
      radioButtons("disp", "Display the data", choices = c(Head = "head", All = "all"), selected = "head"),
      
      tags$hr(),
      
      # Choose the moiety from which the product ion in the MRM transition originates 
      selectInput("productIon", label = "2) Choose origin of product ion", choices = list("Head Group", "FA", "LCB", "Neutral", "RPLC")),
      
      # Choose the lipid class and corresponding MRM pattern
      uiOutput("lipidClassSelector"),
      
      #Notes
      HTML("<b>Notes</b>:<br>
           <br>
           - R packages needed: 'enviPat'<br>
           <br>
           - When correcting the intensities of very low abundance peaks or when the interfering tranistions is very intense, measurement variations may result in negative corrected values. LICAR replaces these negative values by 0 and does not proceed with further correction.<br>
           <br>
           <b>"),
      
      tags$a(href="https://github.com/SLINGhub/LICAR/blob/main/manual/LICAR_Manual.pdf", "Download the user manual.")
      
    ),
    
    
    mainPanel(
      
      img(src='SLING.png', width="40%"),
      
      tags$hr(),
      
      
      HTML("<br/><br/>
           <b>Please verify the data below were uploaded correctly (Data in the same lipid class).
              If yes, click the 'Isotopic correction' button below to perform isotopic correction.
              If no, please check and re-upload your raw data file.</b>
           <br/><br/>"),
      
      #Check whether data belong to one class
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
           <b>Please click the 'Download Ratio' button below to calculate and download the ratio of peak areas between after and before isotopic correction.</b>
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
  
  
  # Display the lipid class and MRM pattern depending on the selection of product ion origin
  output$lipidClassSelector <- renderUI({
    
    #Choose the choices1 based on input$productIon
    if(input$productIon %in% "Head Group") {
      choices1 = list( "AcylCarnitine (Pos) Pro=85" = "AcylCarnitine",
                       "LPC (Pos) Pro=184.1" = "LPC",
                       "LPC (Pos) Pro=104.1" = "LPCql",
                       "LPC-O (Pos) Pro=104.1" = "LPCO",
                       "LPC-O (Pos, qualifier) Pro=184.1" = "LPCOql",
                       "LPE (Pos) Pre-Pro=141" = "LPE",
                       "LPE (Neg) Pro=196.1" = "LPENHG",
                       "PC (Pos) Pro=184.1" = "PC",
                       "PCO (Pos) Pro=184.1" = "PCO",
                       "PCP (Pos) Pro=184.1" = "PCP",
                       "PE (Pos) Pre-Pro=141" = "PE",
                       "PE (Neg) Pro=196.1" = "PENHG",
                       "PG (Pos) Pre-Pro=189" = "PG",
                       "PG (Neg) Pro=153" = "PGNHG",
                       "PI (Pos) Pre-Pro=277" = "PI",
                       "PI (Neg) Pro=241" = "PINHG",
                       "PS (Pos) Pre-Pro=185" = "PS",
                       "PS (Neg) Pre-Pro=87" = "PSNHG",
                       "S1P (Pos) Pro=60.1" = "S1P",
                       "S1Pql (Pos, QL) Pro=113" = "S1Pql",
                       "SM (Pos) Pro=184.1" = "SM")
    } else if(input$productIon %in% "FA") {
      choices1 = list( "CL (Neg) FA" = "CLNFA",
                       #"CL (Neg, doubly charged) FA" = "CLNFA_2",
                       "LPC (Neg, Formate adduct) FA" = "LPCNFA",
                       "LPC (Neg, acetate adduct) FA" = "LPCNFA_2",
                       "LPE (Neg) FA" = "LPENFA",
                       "LPI (Neg) FA" = "LPINFA",
                       "LPG (Neg) FA" = "LPGNFA",
                       "PC (Neg, Formate adduct) FA" = "PCNFA",
                       "PCO (Neg, Formate adduct) FA" = "PCONFA",
                       "PCP (Neg, Formate adduct) FA" = "PCPNFA",
                       "PC (Neg, acetate adduct) FA" = "PCNFA_2",
                       "PCO (Neg, acetate adduct) FA" = "PCONFA_2",
                       "PCP (Neg, acetate adduct) FA" = "PCPNFA_2",
                       "PE (Neg) FA" = "PENFA",
                       "PEO (Neg) FA" = "PEONFA",
                       "PE-P (Pos) FA" = "PEP",
                       "PG (Neg) FA" = "PGNFA",
                       "PI (Neg) FA" = "PINFA",
                       "PS (Neg) FA" = "PSNFA")
    } else if(input$productIon %in% "LCB") {
      choices1 = list(  "Hex1Sph (Pos) SphB-H2O" = "Hex1Sph",
                        "Cer (Pos) SphB-2H2O" = "Cer",
                        "deoxyCer (Pos) SphB-H2O" = "deoxyCer",
                        "dhCer (Pos) SphB-H2O" = "dhCer",
                        "dhCer (Pos) SphB-2H2O" = "dhCer_2",
                        "GM3 (Pos) SphB-2H2O" = "GM3",
                        "Hex1Cer (Pos) SphB-2H2O" = "Hex1Cer",
                        "Hex2Cer (Pos) SphB-2H2O" = "Hex2Cer",
                        "Hex3Cer (Pos) SphB-2H2O" = "Hex3Cer")
    } else if(input$productIon %in% "Neutral") {
      choices1 = list(  "CE (Pos) FANL" = "CE",
                        "DG (Pos) FANL" = "DG",
                        "TG (Pos) FANL" = "TG")
    } else if(input$productIon %in% "RPLC") {
      choices1 = list(  "SM -> PC Pro=184.1" = "PC",
                        "PC-P -> PC-O Pro=184.1" = "PCO")
    }
    
    
    
    # Check the lipid class is unique or not
    output$classError <- renderText({
      class_name <- sub( " .*", "", rownames(rawData()) )
      if( (length(table(class_name)) > 1) && (input$productIon %in% c("Head Group", "FA", "LCB", "Neutral"))) {
        stop(paste("Lipid class is not unique, please check the data! Continue if it belongs to group RPLC."))
      }})
    
    
    # Display the selected raw data table
    output$contents <- renderTable({
      if(input$disp == "head") {return(head(rawData()))}
      else {return(rawData())}
    }, rownames = 1)
    
    
    selectInput("lipidClass",label = "3) Choose lipid class and MRM pattern:",
                choices = choices1)
    
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
