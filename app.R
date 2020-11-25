## app to navigate taxonomic and spatial data
# original navigation principle taken from:
# https://stackoverflow.com/questions/49268902/get-axis-ranges-after-plotly-resize-in-shiny/49526460

library(ALA4R)
library(shiny)
library(plotly)
library(shinyWidgets)
library(sf)
library(viridisLite)

# data
aus <- readRDS("data/small_aus_map.rds")
aus_grid <- readRDS("data/aus_grid.rds")
map_counts <- readRDS("data/hex_counts.rds")
taxon_hierarchy <- c("all", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# messy way to build map
aus_points <- as.data.frame(st_coordinates(st_centroid(aus_grid)))
aus_radii <- data.frame(
  x = aus_points$X,
  y = aus_points$Y,
  radius = scales::rescale(map_counts, to = c(0.3, 2.0)))
map_circles <- make_circles(aus_radii[order(aus_radii$radius, decreasing = TRUE), ])
map_circles$count <- rep(map_counts, each = 100)
map_circles$label <- paste0("n = ", formatC(map_circles$count, big.mark = ","))
map_circles$color <- rep(
  viridis(n = nrow(aus_points), direction = -1),
  each = 100
)
# note: div should make text invisible; but doesn't


# interface
ui <- fluidPage(
  fluidRow(
    br(), br(), br() # could add ALA logo here somewhere
  ),
  fluidRow(
    column(width = 2),
    column(width = 7,
      textInput(
        inputId = "search_bar",
        label = NULL,
        width = "100%"
      )
    ),
    column(width = 1,
      actionButton(
        inputId = "search_go",
        label = "",
        icon = icon("search")
      )
    ),
    column(width = 2)
  ),
  fluidRow(
    column(width = 1),
    column(width = 5,
      br(),
      plotlyOutput("aus_map"),
      materialSwitch(
        inputId = "log_scale",
        label = "log counts",
        value = TRUE,
        right = TRUE,
        status = "success"
      )
    ),
    column(width = 5,
      br(),
      plotlyOutput("taxon_plot"),
      textOutput("taxon_current"),
      uiOutput("taxon_navigation")
    ),
    column(width = 1)
  ),
  fluidRow(
    br(),
    htmlOutput("xlims"),
    br(),
    h4("Verbatim plotly `relayout` data"),
    verbatimTextOutput("relayout"),
    verbatimTextOutput("map_x"),
    h4("Verbatim taxonomy data"),
    verbatimTextOutput("display_previous"),
    verbatimTextOutput("display_current"),
    verbatimTextOutput("display_next"),
    verbatimTextOutput("sequence"),
    verbatimTextOutput("name_previous"),
    verbatimTextOutput("name_current"),
    verbatimTextOutput("count_df")
  )
)

server <- function(input, output, session) {

  ## STORED VALUES ##
  map_extent <- reactiveValues(
    x0 = 110,
    x1 = 155,
    y0 = (-45),
    y1 = (-10),
    bbox = st_polygon(list(rbind(
      c(110, -45), # c(0,0),
      c(155, -45), # c(1,0),
      c(155, -10), # c(1,1),
      c(110, -10), # c(0,1),
      c(110, -45)
    ))), # c(0, 0)
    counts = readRDS("data/hex_counts.rds")
  )

  taxonomy <- reactiveValues(
    display_double_previous = NULL,
    display_previous = NULL,
    display_current = "kingdom",
    display_next = "phylum",
    # display_order = 1, # option to replace explicit naming of current/previous levels
    sequence = c(
      kingdom = NA,
      phylum = NA,
      class = NA,
      order = NA,
      family = NA,
      genus = NA,
      species = NA))

  selections <- reactiveValues(
    name_current = NULL,
    name_previous = NULL,
    ok = FALSE,
    run_ala_call = 0,
    redraw_taxonomy = 0)

  taxon_data <- reactiveValues(
    ala_call = readRDS("data/kingdom_counts.rds"),
    counts = readRDS("data/kingdom_plot.rds"),
    ala_call_previous = NULL,
    counts_previous = NULL)


  ## PLOTS ##
  # draw map
  output$aus_map <- renderPlotly({

    plot_ly(
      type = "scatter",
      mode = "lines",
      source = "map"
    ) %>%
      add_sf(
        color = I("grey80"),
        hoverinfo = "none",
        data = aus
      ) %>%
      add_polygons(
        x = ~x,
        y = ~y,
        split = ~label,
        # color = ~count,
        color = I(map_circles$color),
        alpha = 1,
        hoverinfo = "text",
        hoverlabel = list(
          bgcolor = grey(0.9),
          bordercolor = grey(0.9),
          namelength = 200,
          font = list(color = "black")
        ),
        data = map_circles
      ) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = FALSE)

  })

  # create the taxonomy plot
  output$taxon_plot <- renderPlotly({
    plot_ly(taxon_data$counts,
      type = "scatter",
      mode = "lines",
      fill = "toself", # needed?
      source = "taxon",
      x = ~x,
      y = ~y,
      split = ~plotly_text,
      color = I(taxon_data$counts$color),
      hoverinfo = "text",
      hoverlabel = list(
        bgcolor = grey(0.9),
        bordercolor = grey(0.9),
        namelength = 200,
        font = list(color = "black")
      )
    ) %>%
    add_polygons(
      hoveron = "fills",
      alpha = 1) %>%
    layout(
      showlegend = FALSE,
      xaxis = list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        scaleanchor = "x" # set aspect ratio =  1
      )
    ) %>%
    config(
      displayModeBar = FALSE
    )
  })


  ## PLOT REACTIVITY ##
  # capture selections of taxonomic data
  # Note: observeEvent critical here to avoid premature updating caused by observe()
  observeEvent(event_data(event = "plotly_click", source = "taxon"),
    {
      # note:
        # for some reason curves are now added twice, hence doubling of curveNumber
        # they also index alphabetically, hence sort(levels()) code given here
      click_curve <- event_data(event = "plotly_click", source = "taxon")$curveNumber
      click_curve <- click_curve - max(taxon_data$counts$group) + 1
      available_taxa <- sort(levels(taxon_data$counts$label))
      selected_taxon <- available_taxa[click_curve]
      selections$name_previous <- selections$name_current
      selections$name_current <- selected_taxon
      selections$run_ala_call <- selections$run_ala_call + 1
    }
  )

  # handle map extent updates
  observeEvent(event_data(event = "plotly_relayout", source = "map"),
    {
      map_event <- event_data(event = "plotly_relayout", source = "map")
      if(any(names(map_event) == 'xaxis.range[0]')){
        map_extent$x0 <- map_event$'xaxis.range[0]'
        map_extent$x1 <- map_event$'xaxis.range[1]'
        map_extent$y0 <- map_event$'yaxis.range[0]'
        map_extent$y1 <- map_event$'yaxis.range[1]'
      }else{
        map_extent$x0 <- 110
        map_extent$x1 <- 155
        map_extent$y0 <- (-45)
        map_extent$y1 <- (-10)
      }
      # create a bounding box
      map_extent$bbox <- st_polygon(list(rbind(
        c(map_extent$x0, map_extent$y0), # c(0,0),
        c(map_extent$x1, map_extent$y0), # c(1,0),
        c(map_extent$x1, map_extent$y1), # c(1,1),
        c(map_extent$x0, map_extent$y1), # c(0,1),
        c(map_extent$x0, map_extent$y0) # c(0, 0)
      )))
      selections$run_ala_call <- selections$run_ala_call + 1
    }
  )

  # if either event is triggered, run ala call
  observeEvent(selections$run_ala_call, {
    if(selections$run_ala_call > 0){

      if(is.null(selections$name_current)){ # at first no name is given
        taxon_data$ala_call <- ala_counts(
          breakdown = taxonomy$display_current,
          area = map_extent$bbox,
          limit = 50)
      }else{ # but if there is a name, provide counts for that taxon only
        taxon_df <- ala_taxa(selections$name_current)
        taxon_data$ala_call <- ala_counts(
          taxon_id = taxon_df$taxon_concept_id,
          breakdown = taxonomy$display_next,
          area = map_extent$bbox,
          limit = 50)
      }

      # create or update taxonomy plot data once a name is selected
      if(nrow(taxon_data$ala_call) > 0){
        if(!is.null(taxon_data$counts) & !is.null(selections$name_current)){
          # update taxon hierarchy
          temporary_sequence <- isolate(taxonomy$sequence)
          temporary_sequence[taxonomy$display_current] <- selections$name_current
          taxonomy$sequence <- temporary_sequence
        }
        selections$redraw_taxonomy <- selections$redraw_taxonomy + 1 # ok to update taxon labelling
      }else{ # if no data, say so
        showModal(
          modalDialog(
            HTML(paste0(
              taxonomy$display_current,
              " ",
              selections$name_current,
              " has no documented children. Click anywhere to exit." # OR no data for this area
            )),
            title = "No data found",
            footer = NULL,
            easyClose = TRUE
          )
        )
      }
    }
  })

  # if ok to proceed, set up df to draw new taxonomy plot
  observeEvent({
    input$log_scale
    selections$redraw_taxonomy
  }, {
    if(selections$redraw_taxonomy > 0){
      taxon_data$counts <- taxon_plot_data(
        labels = taxon_data$ala_call$name,
        counts = taxon_data$ala_call$count,
        gap = 0.1,
        log_scale = input$log_scale)
    }
  })


  ## TAXONOMY NAVIGATION ##
  # update downward navigation of taxonomic hierarchy
  observeEvent(selections$redraw_taxonomy, {
    if(#selections$ok &
      !is.null(selections$name_current)
    ){
      taxonomy$display_double_previous <- taxonomy$display_previous
      taxonomy$display_previous <- taxonomy$display_current
      taxonomy$display_current <- taxonomy$display_next
      taxonomy$display_next <- taxon_hierarchy[which(taxon_hierarchy == taxonomy$display_next) + 1]
    }
  })

  # add label to current plot
  output$taxon_current <- renderText({
    if(taxonomy$display_current == "kingdom"){
      "Showing all Kingdoms"
    }else{
      paste(
        taxonomy$display_current,
        "in",
        tools::toTitleCase(taxonomy$display_previous),
        selections$name_current,
        collapse = " ")
    }
  })

  # add 'back' button
  output$taxon_navigation <- renderUI({
    if(selections$redraw_taxonomy > 0 & !is.null(taxonomy$display_previous)){
      if(taxonomy$display_previous == "kingdom"){
        button_text <- "all Kingdoms"
      }else{
        button_text <- paste(
          tools::toTitleCase(taxonomy$display_double_previous),
          selections$name_previous,
          sep = " ")
      }
      actionButton(
        inputId = "taxon_up",
        icon = icon("chevron-circle-left"),
        label = button_text)
    }
  })

  # make back button functional


  ## SEARCHING ##


  ## TESTING ##
  # print the verbatim event_data for plotly_relayout
  output$relayout <- renderPrint({event_data(event = "plotly_click", source = "taxon")})
  output$map_x <- renderPrint({map_extent$x0})

  # show taxonomic data
  output$display_previous <- renderPrint({taxonomy$display_previous})
  output$display_current <- renderPrint({taxonomy$display_current})
  output$display_next <- renderPrint({taxonomy$display_next})
  output$sequence <- renderPrint({taxonomy$sequence})
  # output$selected_ok <- renderPrint({selections$increment})
  output$name_previous <- renderPrint({selections$name_previous})
  output$name_current <- renderPrint({selections$name_current})
  output$count_df <- renderPrint({str(taxon_data$counts)})

}

shinyApp(ui, server)