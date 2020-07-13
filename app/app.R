library(ggplot2)
library(jpeg)
library(Matrix)
library(grid)
library(magrittr)
library(magick)
library(zeallot)
library(yaml)
library(extrafont)

library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)

source("utils.R")
options(shiny.usecairo = FALSE)

CONFIG <- read_yaml("config.yaml")
PATIENTS <- names(CONFIG)


#### pre-run setup ####

theme_empty <- make_empty_theme()

var.features.list <- readRDS("data/var.features.list")
imdims.list <- readRDS("data/imdims.list")
dims.list <- readRDS("data/dims.list")
data.list <- readRDS("data/data.list")

c(t1.list,
  t2.list,
  t3.list) %<-% lapply(1:3,
                       function(x) {readRDS(file.path("data",
                                                      paste0("t",
                                                             x,
                                                             ".list")))})

c(t1.cells,
  t2.cells,
  t3.cells) %<-% lapply(list(t1.list,
                             t2.list,
                             t3.list),
                        function(x){Reduce(intersect,
                                          lapply(x,colnames))}
                        )

#### UI ####

ui <- dashboardPage(

  header = dashboardHeader(title = "HER2+ BC"),

  sidebar = dashboardSidebar(
    width = 350,
    column(width = 12,
           do.call(sidebarMenu,
                   c(lapply(PATIENTS,
                     function(p){
                       do.call(dropdownButton,
                               c(lapply(1:CONFIG[[p]]$n_sections, function(i) {
                                 menuItem(paste0("Section ", i),
                                 tabName = paste0(p,"_section_", i),
                                 icon = icon("angle-double-right"))
                               }),
          status = "default",
          circle = FALSE,
          size = "sm",
          label = paste0("Patient ",p)))}),
          id = "tabs"
          ))),

    column(width = 12,
           uiOutput("var_features"),
           selectInput(
             inputId = "t1",
             label = "Cell type | Tier 1",
             choices = t1.cells,
             selected = "Epithelial"),
           selectInput(
             inputId = "t2",
             label = "Cell type | Tier 2",
             choices = t2.cells,
             selected = "Epithelial_cancer"),
           selectInput(
             inputId = "t3",
             label = "Cell type | Tier 3",
             choices = t3.cells,
             selected = "Epithelial_cancer_Her2_SC"),
           sliderInput(
             inputId = "alpha",
             label = "Opacity", value = 1,
             min = 0, max = 1, step = 0.01
           ),
           selectInput(
             inputId = "cscale",
             label = "colors",
             choices = names(COLORS),
             selected = "Green"),
           radioButtons(
             inputId = "edgecolor",
             label = "edgecolor",
             choices = names(EDGESTROKES),
             selected = "On")

    )

  ),

  body = dashboardBody(
    tags$head(
           tags$link(rel = "stylesheet",
                     type = "text/css",
                     href = "css/custom.css"
                     )
         ),

    uiOutput("STplot")

  )
)

### SERVER ####
server <- function(input, output, session) {


  rv <- reactiveValues(lastBtn = character())
  observeEvent(input$var, {

    if (input$var > 0 ) {
      rv$lastBtn = "gene"
    }
  })
  observeEvent(input$t1, {
    if (input$t1 > 0 ) {
      rv$lastBtn = "t1"
    }
  })
  observeEvent(input$t2, {
    if (input$t2 > 0 ) {
      rv$lastBtn = "t2"
    }
  })
  observeEvent(input$t3, {
    if (input$t3 > 0 ) {
      rv$lastBtn = "t3"
    }
  })

  output$STplot <- renderUI({

    ls.A <- lapply(PATIENTS, function(p) {
            lapply(1:CONFIG[[p]]$n_sections,
                   function(i) {
                     tabItem(tabName = paste0(p, "_section_", i),
                             class = ifelse((p == "A") & (i == 1),
                                            yes = "active",
                                            no = ""),
                             fluidRow(column(width = 12,
                                             tags$head(tags$style(paste0("#",
                                                                         p,
                                                                         "plot",
                                                                         i,
                                                                         "{width:85vh !important; height:85vh !important;}"))),
                                             div(id = "container",
                                                 height = imdims.list[[p]][[i]][1],
                                                 width = imdims.list[[p]][[i]][2],
                                                 style = "position:relative;",
                                                 div(tags$img(src = file.path(CONFIG[[p]]$image_directory,
                                                                              CONFIG[[p]]$image_names[i]),
                                                              style = paste0("width: 85vh; height: 85vh;")),
                                                     style = "position:absolute; top:0; left:0;"),
                                                 div(plotOutput(paste0(p, "plot", i)),
                                                     style = "position:absolute; top:0; left:0;")
                                                 )
                                             ))
                             )}
                   )})

      ls.A <- unlist(ls.A,recursive =  F)

    do.call(tabItems, ls.A)

  })

  output$var_features <- renderUI({

    index <- substr(input$tabs,1,1)
    selectInput(
      inputId = "var",
      label = "Gene",
      choices = var.features.list[[index]],
      selected = var.features.list[[index]][1])

  })

  get_data <- reactive({
    index <- substr(input$tabs, 1, 1)
    c(gc_data, data) %<-% data.list[[index]]


    if (rv$lastBtn == "t1") {
      sc.data <- t1.list[[index]]
      data[, input$t1] <- sc.data[, input$t1]
    } else if (rv$lastBtn == "t2") {
      sc.data <- t2.list[[index]]
      data[, input$t2] <- sc.data[, input$t2]
    } else if (rv$lastBtn == "t3") {
      sc.data <- t3.list[[index]]
      data[, input$t3] <- sc.data[, input$t3]
    } else if (rv$lastBtn == "gene") {
      data[, input$var] <- gc_data[input$var, ]
    }

    dims <- lapply(dims.list[[index]],
                   function(x) x[2:3] %>% as.numeric())


    dfs <- lapply(dims, diff)

    return(list(data,
                dims,
                dfs,
                variable = ifelse(rv$lastBtn %in% names(input),
                                  input[[rv$lastBtn]],
                                  input$var
                                  )
                  ))

  })

  # re-render the plot with the new data -------------------------

  lapply(PATIENTS,function(p){
    lapply(1:CONFIG[[p]]$n_sections, function(i) {
      output[[paste0(p,"plot", i)]] <- renderPlot({
        c(dt, dims, dfs, variable) %<-% get_data()

        ggplot() +
          geom_point(data = subset(dt, sample == paste0(i)),
                     mapping = aes_string(x = "warped_x",
                                          y = paste0("dims[[", i, "]][2] - warped_y"),
                                          fill = paste0("`", variable, "`")
                                          ),
                     stroke = EDGESTROKES[[input$edgecolor]],
                     size = session$clientData$output_Aplot1_width/150,
                     alpha = input$alpha,
                     shape = 21) +

          theme_empty +

          ggtitle(ifelse(rv$lastBtn %in% c("t1", "t2", "t3"),
                        paste0("Cell type: ", variable),
                        paste0("Gene: ", variable))) +

          labs(fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"),
                            "Cell type \nproportion",
                            "Expression\n(Normalized)")) +

          scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
          scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
          scale_fill_gradientn(colours = COLORS[[input$cscale]])
              },
              bg = "transparent")
    })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
