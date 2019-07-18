## Interactive app to calculate cumulative PMI for a synthetic sequence
# December 2016 - June 2018
# Jacob Albrecht, Bristol-Myers Squibb

library(shiny)
library(DiagrammeR) #for graphVis, version 0.9.0
library(reshape2)
library(ggplot2)
library(mgcv) # for rmvn
library(readxl)
library(markdown)
options(encoding = 'UTF-8')

presets_file <- 'Presets.xlsx'

mv.samples <-
  function(xmean,
           xsd,
           ymean,
           ysd,
           correlation = -0.53,
           N.iter = 5000) {
    # simple function to calculate bivariate normal distributions given means, sds and correlation
    rmvn(N.iter, c(xmean, ymean), matrix(
      c(xsd ^ 2, correlation * xsd * ysd, correlation * xsd * ysd, ysd ^ 2),
      nrow = 2
    ))
  }

pmi.calc <- function(conv.df,
                     N.iter = 5000,
                     z = 5.15,
                     correlation = -0.53) {
  
  # set seed for reproducible results
  set.seed(8675309)
  
  # create system of equations,
  all.labs <- unique(c(conv.df$inlabels, conv.df$outlabels))
  n <- length(all.labs)
  A <- matrix(rep(0, n ^ 2), nrow = n)
  diag(A) <- rep(1, n)
  #A.low <- A
  #A.high <- A
  nr <- nrow(A)
  nc <- ncol(A)  # nr==nc , because the matrix should be square
  A.mc <-
    array(rep(A, N.iter), dim = c(nr, nc, N.iter))  # giant matrix of all conditions to be solved
  # identify product from gap in labels
  prod <- all.labs[which(!all.labs %in% conv.df$inlabels)]
  

  
  w <- which('Waste' == all.labs)  # for indexing in the A.mc matrix
  
  yield.mc <- matrix(nrow = length(all.labs), ncol = N.iter)
  step.pmi <- matrix(nrow = length(all.labs), ncol = N.iter)
  
  for (out.lab in unique(conv.df$outlabel)) {

    r <- which(out.lab == all.labs) # row for A matrix
    
    # collect yield and pmi metrics for bivariate sampling
    ry = which(conv.df$outlabels == out.lab &
                 conv.df$inlabels != "Waste")[1]
    mean.yield = mean(c(conv.df$yield.high[ry], conv.df$yield.low[ry]))
    sd.yield = diff(c(conv.df$yield.high[ry], conv.df$yield.low[ry])) /z
    
    rp = which(conv.df$outlabels == out.lab &
                 conv.df$inlabels == "Waste")[1]
    # for debuggin
    mean.pmi = mean(c(conv.df$yield.high[rp], conv.df$yield.low[rp]))
    sd.pmi = diff(c(conv.df$yield.high[rp], conv.df$yield.low[rp])) / z
    
    
    # pull the MC samples: 
    yield.pmi.mc <-
      mv.samples(mean.yield,
                 sd.yield,
                 mean.pmi,
                 sd.pmi,
                 correlation,
                 N.iter = N.iter)
    
    mc.violate.constraints <- which(apply(yield.pmi.mc,1,FUN=function(x){
      if ((x[1]>0 ) & (x[1]<=1) & (x[2]>=1)){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }))
    
    # replace samples that violate constraint with samples of passing conditions:
    if (length(mc.violate.constraints)>0){
      resample.rows <- sample(setdiff(seq(N.iter),mc.violate.constraints),length(mc.violate.constraints),replace = TRUE)
      yield.pmi.mc[mc.violate.constraints,] <- yield.pmi.mc[resample.rows,]
    }
    
    yield.mc[r,] <- yield.pmi.mc[,1]
    S.mc <- yield.pmi.mc[,2]
    
    step.pmi[r,] <- S.mc
    
    for (in.lab in conv.df$inlabel[conv.df$outlabel == out.lab]) {
      #loop through all input materials to calculate RM and S:
      if (in.lab != "Waste") {
        i <-
          which(conv.df$outlabels == out.lab &
                  conv.df$inlabels == in.lab) # i should be unique! - the index in the input dataframe
        c <-which(conv.df$inlabels[i] == all.labs) # the index in the A matrix
        stoich.mc <- runif(N.iter,max=max(conv.df$stoich.high[i],conv.df$stoich.low[i]),
                           min=min(conv.df$stoich.high[i],conv.df$stoich.low[i]))# uniform sampling for stoich
        rm.mc <-
          conv.df$mw.in[i] / (conv.df$mw.out[i] * yield.mc[r, ]) * stoich.mc
        A.mc[c, r, ] <- -1 * rm.mc # put RM data into A matrix
        S.mc <- S.mc - rm.mc # as the RM values are calculated, subtract thesee from PMI to equal S when the loop is done
      }
      
      A.mc[w, r, ] <- -S.mc # put S values into A matrix
    }
  }
  
  b <- matrix(rep(0, n), nrow = n)
  b[prod == all.labs] <- 1 # product is set to 1 kg basis
  
  # TODO: error/loop checking on graph, no duplicate conversions, one and only one product, nodes labeled "waste" etc.
  
  # solve linear algebra
  soln.mc <- apply(
    A.mc,
    3,
    FUN = function(X) {
      solve(X, b)
    }
  )

  row.names(soln.mc) <- all.labs

  row.names(yield.mc) <- all.labs
  row.names(step.pmi) <- all.labs
  #print('solved pmi')
  #print(summary(t(soln.mc)))
  
  soln.low <- apply(soln.mc, 1, function(x) {
    quantile(x, probs = 0.05,na.rm = TRUE)
  })
  soln.high <- apply(soln.mc, 1, function(x) {
    quantile(x, probs = 0.95,na.rm =TRUE)
  })
  
  # dump data here for debuggin:
  #save.image(file="debug.RData")
  
  return(
    list(
      ranges = data.frame(
        kg.per.kg.product.low = soln.low,
        kg.per.kg.product.high = soln.high,
        row.names = all.labs
      ),
      mc = soln.mc,
      yield.mc = yield.mc,
      step.pmi = step.pmi,
      conv.df = conv.df
    )
  )
}

ui <- shinyUI(navbarPage(
  "PMI Calculator",
  tabPanel("Information",
    fluidPage( withMathJax(includeMarkdown("README.md")))
           ),
  tabPanel("1. Define Process",
           
           mainPanel(
             fluidRow(column(
               12,
               h3('Process'),
               p(
                 'Define each step in the process here, specifying the input and output of each step as well as the stoichiometry of each input. Starting with the longest linear sequence is recommended to organize the input fields for the steps.'
               ),
               uiOutput('uiOutpt')
             )), # END fluidRow
             fluidRow(
               column(4, div()),
               column(4, actionButton("add", "Add Row")),
               column(4, actionButton("remove", "Remove Row"))

             ), # END fluidRow
             br(),
             br(),

             br(),
             fluidRow(column(6,downloadButton('downloadSpec', 'Download Process Info as RDS file'))),
             fluidRow(column(6,
                      fileInput('uploadSpec', 'Upload Process Info RDS file',accept = c('.rds')))
             ) # END fluidrow
           )),
  # end add features panel
  tabPanel(
    "2. PMI Values",
    id = "PMI_Page",
    sidebarPanel(h3('Step PMI'),
                 p('Two significant figures are used for step yield ranges'),
                 uiOutput('uiOutpt_2'),
                 width = 9),
    mainPanel(grVizOutput('diagram'))
  ),
  tabPanel(
    "3. Results",
    sidebarPanel(
      p(" Press Calculate Button to Update Totals "),
      checkboxInput("advanced","Show Advanced Options"),
      conditionalPanel(condition="input.advanced==true",
        numericInput('N.iter', 'Number of Iterations', value =
                       5000),
        selectInput('sigma_scale','Min and Max definitions',
                    list(`Normal 99% Interval`=5.15,
                         `Normal 95% Interval`=3.92,
                         `Normal 50% Interval`=1.35),
                    selected = 5.15),
        numericInput('correlation', 'PMI/Yield Correlation', value =
                       -0.53)
        ),
      actionButton("goButton", "Calculate!"),
      br(),

      br(),
      downloadButton('downloadData', 'Download Monte Carlo Results')
      
    ),
    
    mainPanel(tabsetPanel(
      tabPanel("Overall PMI", plotOutput('OverallPMI')),
      tabPanel(
        "Step Metrics",
        p("The plots below show the distribution of the intermediate compounds needed to produce
          one unit mass of final product.  The panel labeled waste is the sum of all processing 
          masses that are not chemical intermediates"),
        br(),
        textOutput("text2"),
        plotOutput('ggdensity')
        #tableOutput('tbl')
      ),
      tabPanel("Step Yield vs Step PMI", plotOutput('yield.v.pmi'))
    ))
  ),
  id = 'inNavbar'
) # END Navbarpage
)

server <- shinyServer(function(input, output, session) {
  # read in presets
  presets.df <- data.frame(read_excel(presets_file))
  
  # for brevity, define function to collect step values:
  stepinfo <- function(i) {
    svl <- paste0('stoich.low_', i)
    svh <- paste0('stoich.high_', i)
    vn <- paste0('InLabel', i)
    vo <- paste0('OutLabel', i)
    mwn <- paste0('mw_', input[[vn]])
    mwo <- paste0('mw_', input[[vo]])
    fv <- paste0('conv_', input[[vo]])
    
    # tab 1 info
    si.1 <- data.frame(
      inlabels = input[[vn]],
      outlabels = input[[vo]],
      stoich.low = input[[svl]],
      stoich.high = input[[svh]],stringsAsFactors = FALSE
    )
    # tab 2 info
    if (all(is.finite(c(input[[mwn]],input[[mwo]],input[[fv]])))){
      si.2 <- data.frame( 
        mw.in = input[[mwn]],
        mw.out = input[[mwo]],
        yield.low = input[[fv]][1],
        yield.high = input[[fv]][2]
      )
      
      if (nrow(si.1)==nrow(si.2)){
        return(cbind(si.1,si.2))
      }else{
        return(si.1)
    }
    }else{
      return(si.1)
    }
  }
  
  # initial value for features:
  features <- reactiveValues(
    renderd = c(1, 2, 3),
    conv.range = data.frame(
      row.names = c("B", 'D'),
      low = c(0.5, 0.5),
      high = c(0.95, 0.95)
    ),
    inlabels = c('A', 'B', 'C'),
    outlabels = c('B', 'D', 'D'),
    stoich.low = c(1, 1, 1),
    stoich.high = c(1, 1, 1),
    pmi.range = data.frame(
      row.names = c("B", 'D'),
      low = c(10, 30),
      high = c(30, 100)
    )
  )
  
  # initial value for reaction values:
  rxn.info <- reactiveValues()
  
  df <- eventReactive(input$goButton, {
    all.labs <- unique(c(features$inlabels, features$outlabels))
    conv.out <- lapply(features$renderd, stepinfo)
    
    pmi.out <- lapply(seq(all.labs), function(i) {
      # get ids of PMI inputs
      pv.min <- paste0('pmi_min_', all.labs[i])
      pv.max <- paste0('pmi_max_', all.labs[i])
      
      #verify that PMI values are in input list, and create entry for "WASTE"
      if ((pv.min %in% names(input)) & (pv.max %in% names(input))) {
        # "Waste" is considered as input to get the matrix created correctly
        data.frame(
          outlabels = as.character(all.labs[i]),
          inlabels = 'Waste',
          mw.in = NA,
          mw.out = NA,
          stoich.low = NA,
          stoich.high = NA,
          yield.low = input[[pv.min]],
          yield.high = input[[pv.max]]
        )
        #in.per.out.low=input[[pv]][1],
        #in.per.out.high=input[[pv]][2])
        
      }
    })
    
    stopifnot(nrow(pmi.out)==nrow(conv.out))
    
    temp.df <- do.call(rbind, c(conv.out, pmi.out))
    temp.df$inlabels <- as.character(temp.df$inlabels)
    temp.df$outlabels <- as.character(temp.df$outlabels)
    #cbind(temp.df,data.frame(out.per.in.low=1/temp.df$in.per.out.high))
    processinfo <-
      temp.df#cbind(temp.df,data.frame(out.per.in.high=1/temp.df$in.per.out.low))
    # do the MC calculation and return the resulting list of values into df():
    print('Process Info')
    print(processinfo) # spit out information for process
    pmi.calc(processinfo,
             N.iter = input$N.iter,
             z = as.numeric(input[['sigma_scale']]),
             correlation = input[['correlation']])
  })
  
  output$text2 <- renderText({
    all.labs <- unique(c(features$inlabels, features$outlabels))
    prod <- all.labs[which(!all.labs %in% features$inlabels)]
    paste(sprintf(
      "You have selected intermediates: %s \n product: %s",
      c(paste(all.labs, collapse = ", ")),
      prod
    ))
    
  })
  
  output$tbl <- renderTable({
    df()[1]$ranges
  }, rownames = TRUE)
  
  
  output$ggdensity <- renderPlot({
    mcsamps <- df()[2]$mc
    m <- data.frame(t(mcsamps))
    conv.df <- df()[5]$conv.df
    
    all.labs <- unique(c(conv.df$inlabels, conv.df$outlabels))
    prod <- all.labs[which(!all.labs %in% conv.df$inlabels)]
  
    m$pmi <-
      apply(m, 1, sum) - apply(m[names(m) %in% unique(conv.df$outlabels)], 1, sum)
    
    m[prod] <- NULL # remove column for product
    
    mm <- melt(m,variable.name="Component Mass")

    p <-
      ggplot(mm, aes(x = value, fill = `Component Mass`)) + geom_histogram(bins = 30) + facet_wrap( ~
                                                                                     `Component Mass`, scales = 'free')
    p
    
  })
  
  # buttons to download data
  output$downloadData <- downloadHandler(
    filename = 'MC_results.csv',
    content = function(con) {
      write.csv(df()[2]$mc, con)
    }
  )
  
  output$downloadSpec <- downloadHandler(
    filename = 'PMI_specification.rds',
    content = function(con) {
      print('Saving Info')
      saveRDS(list(features=isolate(reactiveValuesToList(features)),
                   rxn.info = isolate(reactiveValuesToList(rxn.info))), file=con, ascii= TRUE)
    }
  )
  
  # to upload past specification:
  observeEvent(input$uploadSpec,{
    
      d <- readRDS(input$uploadSpec$datapath)
      for (n in names(d$features)){
        features[[n]] <- d$features[[n]]
      }
      

      for (n in names(d$rxn.info)){
        rxn.info[[n]] <- d$rxn.info[[n]]
      }
      
      # remove features not present in the uploaded file:
      for (n in setdiff(names(features),names(d$features))){
        features[[n]] <- NULL
      }
      for (n in setdiff(names(rxn.info),names(d$rxn.info))){
        rxn.info[[n]] <- NULL
      }
      print('Upload complete')
      
      #output$uiOutpt_2 <- d$uiOutpt_2
  })
  
  output$OverallPMI <- renderPlot({
    mcsamps <- df()[2]$mc
    conv.df <- df()[5]$conv.df

    m <- data.frame(t(mcsamps))
    m$pmi <-
      apply(m, 1, sum) - apply(m[names(m) %in% unique(conv.df$outlabels)], 1, sum)
    Cumulative_PMI <- m$pmi
    Cpmi_mean <- round(mean(Cumulative_PMI), 1)
    Cpmi_conf_left  <-
      round(quantile(Cumulative_PMI, c(0.025, 0.975),na.rm=TRUE)[1], 1)
    Cpmi_conf_right  <-
      round(quantile(Cumulative_PMI, c(0.025, 0.975),na.rm=TRUE)[2], 1)

    hist(
      Cumulative_PMI,
      breaks = 40,
      col = "blue",
      xlab = "Cumulative_PMI Prediction",
      main = paste(
        "Estimated Cumulative_PMI\n",
        "Mean=",
        Cpmi_mean,
        ";",
        " 95% Conf Int= [",
        Cpmi_conf_left,
        ",",
        Cpmi_conf_right,
        "]"
      )
    )
    abline(v = quantile(Cumulative_PMI, c(0.025, 0.975)), col = "red")
    
  })
  
  output$yield.v.pmi <- renderPlot({
    y <- melt((df()[3]$yield.mc), value.name = 'StepYield')
    
    p <- melt((df()[4]$step.pmi), value.name = 'StepPMI')

    dat <- merge(y, p)
    dat <- dat[complete.cases(dat),]

    ggplot(dat, aes(x = StepYield, y = StepPMI, color = Var1)) + geom_point(alpha =
                                                                              0.8) +
      facet_wrap(~Var1)+ labs(color='Reaction Product')
  })
  

  
  # Increment reactive values array used to store how may rows we have rendered
  observeEvent(input$add, {
    out <- lapply(features$renderd, stepinfo)
    
    df <- do.call(dplyr::bind_rows, out)
  
    features$inlabels <- c(as.character(df$inlabels), '')
    features$outlabels <- c(as.character(df$outlabels), '')
 
    features$renderd <-
      c(features$renderd, length(features$renderd) + 1)
    features$stoich.low <- c(df$stoich.low, 1)
    features$stoich.high <- c(df$stoich.high, 1)
    
  })
  
  observeEvent(input$remove, {
    features$renderd <- features$renderd[-length(features$renderd)]
    out <- lapply(features$renderd, stepinfo)
    df <- do.call(dplyr::bind_rows, out)
  
    if (nrow(df) > 0) {
      features$outlabels <- as.character(df$outlabels)
      features$inlabels <- as.character(df$inlabels)
      features$stoich.low <- df$stoich.low
      features$stoich.high <- df$stoich.high

    }
  })
  
  observeEvent(input$inNavbar, {
    out <- lapply(features$renderd, stepinfo)

    df <- do.call(dplyr::bind_rows, out)
 
    if (nrow(df) > 0) {
      features$outlabels <- as.character(df$outlabels)
      features$inlabels <- as.character(df$inlabels)
      features$stoich.low <- df$stoich.low
      features$stoich.high <- df$stoich.high

    }
  })
  
  
  
  # this version of output$diagram for older DiagrammeR
  
  if (packageVersion('DiagrammeR')<'0.9.0'){
  
  output$diagram<-renderGrViz({
    g1 <- create_graph(nodes_df = create_nodes(nodes=features$inlabels),
                       edges_df= create_edges(from=features$inlabels,to=features$outlabels),
                       generate_dot = TRUE)
    grViz(g1$dot_code)
  })
  }else{
  
  # this version of output$diagram for DiagrammeR v0.9.0
 output$diagram <- renderGrViz({
   nodes <- unique(c(features$inlabels, features$outlabels))
   g1 <-
     create_graph(
       nodes_df = create_node_df(length(nodes), label = nodes),
       edges_df = create_edge_df(
         from = match(features$inlabels, nodes),
         to = match(features$outlabels, nodes)
       )
     )

   g1$global_attrs$value[1] <- 'dot'  # a hack to avoid 'neato' layout
   grViz(generate_dot(g1), engine = "dot")
 })
  
  }
  
  
  observeEvent({lapply(grep('^(?!mw|inNavbar|uploadSpec)',names(input),value=TRUE,perl=TRUE),
                       function(i){(input[[i]])}
                       )}, {
  
    all.labs <- unique(c(isolate(features$inlabels), isolate(features$outlabels)))
    
    # TODO: refactor this nested loop:
    for (i in seq(all.labs)){
      #print(paste0('mw_', all.labs[i]))
      for (key in c('mw_','pmi_min_','pmi_max_','conv_','Preset_')){
        if (!is.null(input[[paste0(key,all.labs[i])]])){
          rxn.info[[paste0(key, all.labs[i])]] <- input[[paste0(key, all.labs[i])]]
        }
      }
    }

  }, priority = -1)
  
  # If reactive vector updated we render the UI again
  observe({
    # create UI panels, this piece is for Tab #1
    #
    
    output$uiOutpt <- renderUI({
      # Create rows
      rows <- lapply(features$renderd, function(i) {
        fluidRow(
          # duplicate choices make selectize unhappy, use unique():
          
          column(
            2,
            numericInput(
              paste0('stoich.low_', i),
              label = 'Stoichiometry Range: Low',
              value = features$stoich.low[i]
            )
          ),
          column(
            2,
            numericInput(
              paste0('stoich.high_', i),
              label = 'Stoichiometry Range: High',
              value = features$stoich.high[i]
            )
          ),
          column(
            3,
            selectizeInput(
              paste0('InLabel', i),
              label = 'Input Name',
              selected = features$inlabels[i],
              choices = unique(c(
                features$inlabels,
                features$outlabels,
                toupper(letters)
              )),
              #unique(c(features$inlabels[i],features$outlabels[!features$outlabels %in% features$inlabels])),
              options = list(create = TRUE)
            )
          ),
          column(
            3,
            selectizeInput(
              paste0('OutLabel', i),
              label = "Output Name",
              selected = features$outlabels[i],
              choices = unique(c(
                features$inlabels,
                features$outlabels,
                toupper(letters)
              )),
              #unique(c(features$inlabels,features$outlabels)),
              options = list(create = TRUE)
            )
          )
        )
      })
      do.call(shiny::tagList, rows)
    })
    
    # RENDER UI for tab 2:
    tab2<- renderUI({
      # Create rows
      
      all.labs <- unique(c(features$inlabels, features$outlabels))
      # check to see if new prodcuts need PMI sliders:
      new.pmi <-
        (setdiff(features$outlabels, row.names(features$pmi.range)))
      for (i in new.pmi) {
        features$pmi.range <-
          rbind(features$pmi.range,
                data.frame(
                  row.names = i,
                  low = 50,
                  high = 350
                ))
        features$conv.range <-
          rbind(features$conv.range,
                data.frame(
                  row.names = i,
                  low = 0.5,
                  high = 0.8
                ))
      }
      
      rows <- lapply(seq(all.labs), function(i) {
        if (all.labs[i] %in% unique(features$outlabels)) {
          pmi.ix <- all.labs[i] == row.names(features$pmi.range)
          
          # a nested if statement to check that the Preset input exists AND it is one of the types,
          # this could be streamlined...
          # TODO: fix the overwriting of custom fields when other presets are used for other steps
          # TODO: user experience of having all rows reset with a change to one row must be addressed, lots of isolate statements are used, do they work?
          if (!is.null(isolate(rxn.info[[paste0('Preset_', all.labs[i])]]))) {  # the row has a preset
            if (isolate(rxn.info[[paste0('Preset_', all.labs[i])]]) %in% presets.df$Type) {
              ix <-
                match(isolate(rxn.info[[paste0('Preset_', all.labs[i])]]), presets.df$Type)
              pmi_min_value <- presets.df$PMI_low[ix]
              pmi_max_value <- presets.df$PMI_high[ix]
              yield_min_max_value <-
                c(presets.df$Yield_low[ix], presets.df$Yield_high[ix])
              
            }else{  #custom option
              pmi_min_value <- isolate(rxn.info[[paste0('pmi_min_', all.labs[i])]])
              pmi_max_value <-
                isolate(rxn.info[[paste0('pmi_max_', all.labs[i])]])
              yield_min_max_value <-
                isolate(rxn.info[[paste0('conv_', all.labs[i])]])
            }
          } else{
            pmi_min_value <- isolate(rxn.info[[paste0('pmi_min_', all.labs[i])]])
            pmi_max_value <-
              isolate(rxn.info[[paste0('pmi_max_', all.labs[i])]])
            yield_min_max_value <-
              c(features$conv.range[pmi.ix, 'low'],
                features$conv.range[pmi.ix, 'high'])
          }
          
          fluidRow(
            column(
              2,
              numericInput(
                paste0('mw_', all.labs[i]),
                label = paste0("MW of ",
                               all.labs[i], " : "),
                value = (rxn.info[[paste0('mw_', all.labs[i])]])
              )
            ),
            column(
              3,
              selectizeInput(
                paste0('Preset_', all.labs[i]),
                label = 'Load from Preset:',
                selected = (rxn.info[[paste0('Preset_', all.labs[i])]]),
                choices = c(list(`custom`='custom'),split(presets.df$Type,presets.df$Category)),
                options = list(create = FALSE)
              ),
              tags$style(type='text/css', 
                         ".selectize-dropdown-content {
                         max-height: 400px; font-size: 12px}"
             )
            ),
            column(
              2,
              numericInput(
                paste0('pmi_min_', all.labs[i]),
                label = paste0("Step PMI, kg/kg Product",
                               all.labs[i], " (Min) : "),
                min = 0,
                max = 20000,
                value = pmi_min_value
              )
            ),
            column(
              2,
              numericInput(
                paste0('pmi_max_', all.labs[i]),
                label = paste0("Step PMI, kg/kg Product",
                               all.labs[i], " (Max) : "),
                min = 0,
                max = 20000,
                value = pmi_max_value
              )
            ),
            column(
              3,
              sliderInput(
                paste0('conv_', all.labs[i]),
                label = "Molar Yield",
                min = 0.0,
                max = 1,
                step = 0.01,
                value = yield_min_max_value
              )
            )
          )
        } else{
          fluidRow(column(
            2,
            numericInput(
              paste0('mw_', all.labs[i]),
              label = paste0("MW of ",
                             all.labs[i], " : "),
              min = 0,
              max = 1000,
              value = (rxn.info[[paste0('mw_', all.labs[i])]])
            )
          ),
          column(6, div()))
        }
      })
      do.call(shiny::tagList, rows)
    })# render tab 2
    output$uiOutpt_2 <- tab2
  })
})

shinyApp(ui, server)  