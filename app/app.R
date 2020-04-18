library(ggplot2)
library(jpeg)
library(grid)
library(shiny)
library(magrittr)
library(shinycssloaders)
library(shinydashboard)
library(magick)
library(zeallot)
library(shinyWidgets)

#### pre-run setup ####

# initiate a ggplot theme for use in plotting
# (just getting rid of everything so we only see the image itself)
theme_empty <- theme_bw()
theme_empty$line <- element_blank()
theme_empty$rect <- element_blank()
theme_empty$strip.text <- element_blank()
theme_empty$axis.text <- element_blank()
theme_empty$plot.title <- element_text(colour = "white", size = 20, hjust = 0, margin = margin(t = 10, b = -20)) #element_blank()
theme_empty$axis.title <- element_blank()
theme_empty$legend.position <- c(0.97, 0.5)
theme_empty$legend.background <- element_rect(fill = "darkgrey", colour = "white")
theme_empty$legend.text <- element_text(colour = "white")
theme_empty$legend.title <- element_text(colour = "white")

colorscales <- list("Blues" = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
                    "GrRd" = c("lightgray", "mistyrose", "red", "dark red"),
                    "RdBuYl" = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),
                    "Spectral" = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")),
                    "viridis" = viridis::viridis(n = 9),
                    "magma" = viridis::magma(n = 9),
                    "Blues" = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
                    "Reds" = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
                    "Greens" = RColorBrewer::brewer.pal(n = 9, name = "Greens"))

edgecolors <- list("On" = "black",
                   "Off" = "none")

numsections <- list("A"=1:6,
                    "B" = 1:4)

LS <- LETTERS[1:8]

# set the image input file
image.files <- list.files(pattern = "_HE_image_dark", path = "www", full.names = F)
image.files <- list("A" = image.files[1:6],
                    "B" = image.files[7:12],
                    "C" = image.files[13:18],
                    "D" = image.files[19:24],
                    "E" = image.files[25:27],
                    "F" = image.files[28:30],
                    "G" = image.files[31:33],
                    "H" = image.files[34:36],
                    "J" = image.files[37:39])
dims.list <- readRDS("data/dims.list")
imdims.list <- readRDS("data/imdims.list")
var.features.list <- readRDS("data/var.features.list")
data.list <- readRDS("data/data.list")
sc.list <- readRDS("data/cs.list")
scs <- Reduce(intersect, lapply(sc.list, colnames))

t1.list <- readRDS("data/t1.list")
t1.cells <- Reduce(intersect, lapply(t1.list, colnames))

t2.list <- readRDS("data/t2.list")
t2.cells <- Reduce(intersect, lapply(t2.list, colnames))

t3.list <- readRDS("data/t3.list")
t3.cells <- Reduce(intersect, lapply(t3.list, colnames))

# Start with sample A
dims <- lapply(dims.list[[1]], function(x) x[2:3] %>% as.numeric())
dfs <- lapply(dims, diff)
c(gc_data, d) %<-% data.list[[1]]
imdims <- imdims.list[[1]]
#var.features <- Reduce(union, var.features.list)
imgs <- image.files[[1]]


#### UI ####
ui <- dashboardPage(

  header = dashboardHeader(title = "HER2+ BC"),

  ## sidebar = dashboardSidebar(
    ## width = 350,

    ## column(width = 12,
    ##        sidebarMenu(selectInput("patient",
    ##                    "Patient: ",
    ##                    c("Patient A" = "A",
    ##                      "Patient B" = "B")))),
    ## column(width = 12,
    ##        selectInput("tabs1",
    ##                    "Section: ",
    ##                                c("Section 1" = "A_section_1",
    ##                                  "Section 2" = "A_section_2")
    ##                                )),

  sidebar = dashboardSidebar(
    width = 350,

    column(width = 12,
      sidebarMenu(
        id = "tabs1",
        do.call(dropdownButton, c(lapply(1:6, function(i) {
          menuItem(paste0("Section ", i),
                   tabName = paste0("A_section_", i),
                   icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient A"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs2",
        do.call(dropdownButton, c(lapply(1:6, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("B_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient B"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs3",
        do.call(dropdownButton, c(lapply(1:6, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("C_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient C"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs4",
        do.call(dropdownButton, c(lapply(1:6, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("D_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient D"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs5",
        do.call(dropdownButton, c(lapply(1:3, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("E_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient E"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs6",
        do.call(dropdownButton, c(lapply(1:3, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("F_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient F"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs7",
        do.call(dropdownButton, c(lapply(1:3, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("G_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient G"))
      )
    ),
    column(width = 12,
      sidebarMenu(
        id = "tabs8",
        do.call(dropdownButton, c(lapply(1:3, function(i) {
          menuItem(paste0("Section ", i), tabName = paste0("H_section_", i), icon = icon("angle-double-right"))
        }), status = "default", circle = FALSE, size = "sm", label = "Patient H"))
      )
    ),

    column(width = 12,
           uiOutput("var_features"),
           selectInput(
             inputId = "t1",
             label = "tier1 cell type",
             choices = t1.cells,
             selected = "Epithelial"),
           selectInput(
             inputId = "t2",
             label = "tier2 cell type",
             choices = t2.cells,
             selected = "Epithelial_cancer"),
           selectInput(
             inputId = "t3",
             label = "tier3 cell type",
             choices = t3.cells,
             selected = "Epithelial_cancer_Her2_SC"),
           sliderInput(
             inputId = "alpha",
             label = "alpha", value = 1,
             min = 0, max = 1, step = 0.01
           ),
           selectInput(
             inputId = "cscale",
             label = "colors",
             choices = names(colorscales),
             selected = "Spectral"),
           radioButtons(
             inputId = "edgecolor",
             label = "edgecolor",
             choices = names(edgecolors),
             selected = "On")

    )

  ),

  body = dashboardBody(

    tags$head(tags$style(HTML('
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #000000;
                                }

                                /* logo when hovered */
                                .skin-blue .main-header .logo:hover {
                                background-color: #000000;
                                }

                                /* navbar (rest of the header) */
                                .skin-blue .main-header .navbar {
                                background-color: #000000;
                                }

                                /* main sidebar */
                                .skin-blue .main-sidebar {
                                background-color: #101010;
                                }

                                /* active selected tab in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                background-color: #808080;
                                }

                                /* other links in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: #FFFFFF;
                                color: #000000;
                                }

                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #0000;
                                }

                                '))),

    uiOutput("STplot")

  )
)



### SERVER ####
server <- function(input, output, session) {


  ## observe({
  ##   s_options <- sapply(numsections[[input$patient]], function(x){paste0(input$patient,"_section_",x)})
  ##   names(s_options) <- as.character(sapply(numsections[[input$patient]],
  ##                                           function(x){paste0("Section ",x)}))
  ## print(s_options)

  ## updateSelectInput(session,
  ##                   "tabs1",
  ##                   choices = s_options,
  ##                   ## selected = names(s_options)[1],
  ##                   )
  ##   })

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

      print(str(tags))

      ls.A <- append(list(
        # Create active tab
        tabItem(tabName = paste0("A_section_", 1), class = "active",
                fluidRow(column(width = 12,
                  tags$head(tags$style(paste0("#Aplot", 1, "{width:85vh !important; height:85vh !important;}"))),
                  div(id = "container",
                      height = dim(img)[1],
                      width = dim(img)[2],
                      style = "position:relative;",
                      div(tags$img(src = imgs[1],
                                   style = paste0("width: 85vh; height: 85vh;")),
                          style = "position:absolute; top:0; left:0;"),
                      div(plotOutput(paste0("Aplot", 1),
                                     height = imdims[[1]][1],
                                     width = imdims[[1]][2]),
                          style = "position:absolute; top:0; left:0; ")
                  )
                ))
        )), lapply(2:length(imdims), function(i) {
          # Add additional tabs
          tabItem(tabName = paste0("A_section_", i),
                  fluidRow(column(width = 12,
                    tags$head(tags$style(paste0("#Aplot", i, "{width:85vh !important; height:85vh !important;}"))),
                    div(id = "container",
                        height = imdims[[i]][1],
                        width = imdims[[i]][2],
                        style = "position:relative;",
                        div(tags$img(src = imgs[i],
                                     style = paste0("width: 85vh; height: 85vh;")),
                            style = "position:absolute; top:0; left:0;"),
                        div(plotOutput(paste0("Aplot", i)),
                            style = "position:absolute; top:0; left:0; ")
                    )
                  ))
          )}
        )
      )


    ls.B <- Reduce(c, lapply(2:9, function(n) {

      imdims <- imdims.list[[n]]
      imgs <- image.files[[n]]

      append(list(
        # Create active tab
        tabItem(tabName = paste0(LS[n], "_section_", 1),
                fluidRow(column(width = 12,
                  tags$head(tags$style(paste0("#", LS[n], "plot", 1, "{width:85vh !important; height:85vh !important;}"))),
                  div(id = "container",
                      height = dim(imgs[1])[1],
                      width = dim(imgs[1])[2],
                      style = "position:relative;",
                      div(tags$img(src = imgs[1],
                                   style = paste0("width: 85vh; height: 85vh;")),
                          style = "position:absolute; top:0; left:0;"),
                      div(plotOutput(paste0(LS[n], "plot", 1)),
                          style = "position:absolute; top:0; left:0;")
                  )
                ))
        )), lapply(2:length(imdims), function(i) {
          # Add additional tabs
          tabItem(tabName = paste0(LS[n], "_section_", i),
                  fluidRow(column(width = 12,
                    tags$head(tags$style(paste0("#", LS[n], "plot", i, "{width:85vh !important; height:85vh !important;}"))),
                    div(id = "container",
                        height = imdims[[i]][1],
                        width = imdims[[i]][2],
                        style = "position:relative;",
                        div(tags$img(src = imgs[i],
                                     style = paste0("width: 85vh; height: 85vh;")),
                            style = "position:absolute; top:0; left:0;"),
                        div(plotOutput(paste0(LS[n], "plot", i)),
                            style = "position:absolute; top:0; left:0;")
                    )
                  ))
              )}
          )
        )
    }))

    do.call(tabItems, c(ls.A, ls.B))

  })

  output$var_features <- renderUI({

    index <- match(substr(input$tabs1, 1, 1), LS)
    var.features <- var.features.list[[index]]
    selectInput(
      inputId = "var",
      label = "gene",
      choices = var.features,
      selected = "ERBB2")

  })

  get_data <- reactive({
    index <- match(substr(input$tabs1, 1, 1), LS)
    print(index)
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

    dims <- lapply(dims.list[[index]], function(x) x[2:3] %>% as.numeric())
    dfs <- lapply(dims, diff)
    return(list(data, dims, dfs, variable = ifelse(rv$lastBtn == "t1", input$t1, ifelse(rv$lastBtn == "t2", input$t2, ifelse(rv$lastBtn == "t3", input$t3,  input$var)))))

  })


  # re-render the plot with the new data -------------------------
  lapply(1:6, function(i) {
    output[[paste0("Aplot", i)]] <- renderPlot({
      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x",
                                        y = paste0("dims[[", i, "]][2] - warped_y"),
                                        fill = paste0("`", variable, "`")
                                        ),
                   color = edgecolors[[input$edgecolor]],
                   ## color = "red",
                   size = session$clientData$output_Aplot1_width/150,
                   alpha = input$alpha,
                   shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:6, function(i) {
    output[[paste0("Bplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:6, function(i) {
    output[[paste0("Cplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:6, function(i) {
    output[[paste0("Dplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:3, function(i) {
    output[[paste0("Eplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:3, function(i) {
    output[[paste0("Fplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:3, function(i) {
    output[[paste0("Gplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })

  # re-render the plot with the new data -------------------------
  lapply(1:3, function(i) {
    output[[paste0("Hplot", i)]] <- renderPlot({

      c(dt, dims, dfs, variable) %<-% get_data()

      sb <- subset(dt, sample == paste0(i))

      ggplot() +
        geom_point(data = subset(dt, sample == paste0(i)),
                   mapping = aes_string(x = "warped_x", y = paste0("dims[[", i, "]][2] - warped_y"), fill = paste0("`", variable, "`")),
                   size = session$clientData$output_Aplot1_width/150, alpha = input$alpha, shape = 21) +
        theme_empty +
        labs(title = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), paste0("Cell type: ", variable), paste0("Gene: ", variable)), fill = ifelse(rv$lastBtn %in% c("t1", "t2", "t3"), "cell type \nproportion", "norm. gene \nexpression")) +
        scale_x_continuous(limits = c(0, dims[[i]][1]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dims[[i]][2]), expand = c(0, 0)) +
        scale_fill_gradientn(colours = colorscales[[input$cscale]])

    },
    bg = "transparent")
  })


}

# Run the application
shinyApp(ui = ui, server = server)
