# Load necessary libraries
library(shiny)
library(leaflet)
library(sf)
library(terra)
library(sp)
library(rgdal)
library(caret)     # For confusionMatrix
library(dplyr)
library(tidyverse)
library(flextable) # For displaying the confusion matrix nicely
library(shinyjs)   # For enabling/disabling buttons
library(mgcv)
library(recipes)
library(survival)

# Set the maximum upload size to 100 MB
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

# Define UI for the application
ui <- fluidPage(
    # Initialize shinyjs
    useShinyjs(),
    
    # Fullscreen leaflet map output with a map ID
    leafletOutput("map", width = "100%", height = "100%"),
    
    # Add a panel for user controls (top right)
    absolutePanel(
        top = 10, right = 10,
        width = 350,
        style = "z-index: 1000; background-color: white; padding: 10px; border-radius: 10px;",
        # UI elements
        fileInput("rasterFile", "Choose Raster File",
                  accept = c('.tif', '.grd', '.asc', '.nc')),
        fileInput("shapefile", "Choose Shapefile (select all files)",
                  multiple = TRUE,
                  accept = c('.shp','.dbf','.sbn','.sbx','.shx','.prj','.cpg','.xml')),
        actionButton("loadData", "Load Data"),
        uiOutput("classColumnUI"),
        # Initialize buttons as disabled
        actionButton("computeMatrix", "Compute Validation Matrix"),
        # Remove the tableOutput and downloadButton from the panel
        # They will be included in the modal dialog instead
    ),
    
    # Remove extra margin/padding and ensure full-screen layout
    tags$style(
        HTML("
              html, body {
                height: 100%;
                width: 100%;
                margin: 0;
                padding: 0;
                overflow: hidden;
              }
              #map {
                height: 100vh !important;
                width: 100vw !important;
                position: absolute;
                top: 0;
                left: 0;
              }
            ")
    )
)

# Define server logic
server <- function(input, output, session) {
    # Reactive values to store raster and shapefile data
    values <- reactiveValues(
        raster = NULL, 
        rasterExtent = NULL, 
        shape = NULL, 
        classColumn = NULL, 
        confusionMatrix = NULL,
        studySiteCenter = NULL  # New reactive value to store study site center coordinates
    )
    
    # Disable buttons initially
    disable("computeMatrix")
    
    # Render the initial Leaflet map with increased maxZoom
    output$map <- renderLeaflet({
        leaflet(options = leafletOptions(maxZoom = 28)) %>%  # Set maxZoom to a higher value
            addProviderTiles("Esri.WorldImagery", options = providerTileOptions(maxZoom = 28))
    })
    
    # Observe the loadData button
    observeEvent(input$loadData, {
        # Clear previous confusion matrix
        values$confusionMatrix <- NULL
        
        # Load raster file
        if (!is.null(input$rasterFile)) {
            rasterPath <- input$rasterFile$datapath
            values$raster <- tryCatch({
                # Use terra package to read raster
                r <- terra::rast(rasterPath)
                
                # Check and transform CRS if necessary
                raster_crs <- terra::crs(r)
                wgs84_crs <- "+proj=longlat +datum=WGS84 +no_defs"
                if (!is.na(raster_crs) && raster_crs != wgs84_crs) {
                    r <- terra::project(r, wgs84_crs)
                }
                r
            }, error = function(e) {
                showNotification(paste("Error reading raster:", e$message), type = "error")
                NULL
            })
        } else {
            values$raster <- NULL
        }
        
        # Load shapefile
        if (!is.null(input$shapefile)) {
            # Create a unique temporary directory for the session
            tempDir <- tempfile()
            dir.create(tempDir)
            
            shapefilePaths <- input$shapefile$datapath
            shapefileNames <- input$shapefile$name
            
            # Copy uploaded files to the temporary directory
            for (i in 1:length(shapefilePaths)) {
                file.copy(shapefilePaths[i], file.path(tempDir, shapefileNames[i]), overwrite = TRUE)
            }
            
            # Identify the .shp file
            shpFiles <- list.files(tempDir, pattern = "\\.shp$", full.names = TRUE)
            if (length(shpFiles) == 0) {
                showNotification("No .shp file found in uploaded files.", type = "error")
                values$shape <- NULL
            } else {
                shpPath <- shpFiles[1]
                
                # Check for presence of .shx and .dbf files
                shpBaseName <- tools::file_path_sans_ext(basename(shpPath))
                requiredExtensions <- c(".shp", ".shx", ".dbf")
                missingFiles <- c()
                for (ext in requiredExtensions) {
                    filePath <- file.path(tempDir, paste0(shpBaseName, ext))
                    if (!file.exists(filePath)) {
                        missingFiles <- c(missingFiles, paste0(shpBaseName, ext))
                    }
                }
                
                if (length(missingFiles) > 0) {
                    showNotification(paste("Missing required shapefile components:", paste(missingFiles, collapse = ", ")), type = "error")
                    values$shape <- NULL
                } else {
                    # Read shapefile using sf package and handle errors
                    values$shape <- tryCatch({
                        shp <- st_read(shpPath)
                        
                        # Ensure geometries are valid
                        shp <- st_make_valid(shp)
                        
                        # Check and transform CRS if necessary
                        if (is.na(st_crs(shp))) {
                            showNotification("Shapefile CRS is missing. Cannot proceed.", type = "error")
                            return(NULL)
                        } else if (st_crs(shp) != st_crs(4326)) {
                            shp <- st_transform(shp, 4326)
                        }
                        shp
                    }, error = function(e) {
                        showNotification(paste("Error reading shapefile:", e$message), type = "error")
                        NULL
                    })
                }
            }
        } else {
            values$shape <- NULL
        }
        
        # Update the map using leafletProxy
        leafletProxy("map") %>%
            clearShapes() %>%
            clearMarkers()  # Clear any existing markers
        
        # Create and add raster extent polygon if raster is available
        if (!is.null(values$raster)) {
            # Get the raster extent
            ext <- terra::ext(values$raster)
            coords <- matrix(c(
                ext[1], ext[3],  # xmin, ymin
                ext[2], ext[3],  # xmax, ymin
                ext[2], ext[4],  # xmax, ymax
                ext[1], ext[4],  # xmin, ymax
                ext[1], ext[3]   # xmin, ymin (close polygon)
            ), ncol = 2, byrow = TRUE)
            
            # Create a polygon from the extent
            raster_extent_polygon <- st_polygon(list(coords)) %>%
                st_sfc(crs = 4326) %>%
                st_sf()
            
            values$rasterExtent <- raster_extent_polygon
            
            # Compute center coordinates of the raster extent
            center_lon <- as.numeric(((ext[1] + ext[2]) / 2))  # (xmin + xmax) / 2
            center_lat <- as.numeric(((ext[3] + ext[4]) / 2))  # (ymin + ymax) / 2
            
            # Store the center coordinates
            values$studySiteCenter <- c(center_lon, center_lat)
            
            # Set map view to the center of the raster extent
            leafletProxy("map") %>%
                setView(
                    lng = center_lon,
                    lat = center_lat,
                    zoom = 15
                ) %>% 
                addRasterImage(values$raster, colors = terrain.colors(256), opacity = 0.8)
                # addPolygons(
                #     data = values$rasterExtent,
                #     color = "black",
                #     weight = 2,
                #     opacity = 1.0,
                #     fillOpacity = 0.5,
                #     fillColor = "red",
                #     group = "RasterExtent"
                # )
        }
        
        # Add shapefile layer if available and valid
        if (!is.null(values$shape) && nrow(values$shape) > 0) {
            leafletProxy("map") %>%
                addPolygons(
                    data = values$shape,
                    color = "blue",
                    weight = 2,
                    opacity = 1.0,
                    fillOpacity = 0.5,
                    group = "Shapefile"
                )
            # Optional: Adjust the map view to include both raster and shapefile
            if (is.null(values$raster)) {
                bounds <- st_bbox(values$shape)
                center_lon <- (bounds$xmin + bounds$xmax) / 2
                center_lat <- (bounds$ymin + bounds$ymax) / 2
                
                # Store the center coordinates
                values$studySiteCenter <- c(center_lon, center_lat)
                
                leafletProxy("map") %>%
                    setView(
                        lng = as.numeric(center_lon),
                        lat = as.numeric(center_lat),
                        zoom = 15
                    )
            }
        }
        
        # Add a marker at the study site center with a marker ID
        if (!is.null(values$studySiteCenter)) {
            leafletProxy("map") %>%
                addMarkers(
                    lng = values$studySiteCenter[1],
                    lat = values$studySiteCenter[2],
                    popup = "Study Site",
                    options = markerOptions(riseOnHover = TRUE),
                    layerId = "studySiteMarker"  # Assign an ID to the marker
                )
        }
        
        # If shapefile is loaded, update the UI to select class column
        if (!is.null(values$shape)) {
            cols <- colnames(values$shape)
            output$classColumnUI <- renderUI({
                selectInput("classColumn", "Select Class Column:", choices = cols)
            })
        } else {
            output$classColumnUI <- renderUI({})
        }
        
        # Enable or disable the computeMatrix button based on data availability
        if (!is.null(values$raster) && !is.null(values$shape)) {
            enable("computeMatrix")
        } else {
            disable("computeMatrix")
        }
    })
    
    # Observe the map's zoom level and show/hide the marker accordingly
    observe({
        req(values$studySiteCenter)
        zoomLevel <- input$map_zoom  # Get the current zoom level
        # Define the zoom threshold above which the marker will disappear
        zoomThreshold <- 15  # Adjust this value as needed
        if (!is.null(zoomLevel)) {
            if (zoomLevel >= zoomThreshold) {
                # Remove the marker when zoomed in beyond the threshold
                leafletProxy("map") %>% removeMarker(layerId = "studySiteMarker")
            } else {
                # Add the marker when zoomed out
                leafletProxy("map") %>% clearMarkers() %>%
                    addMarkers(
                        lng = values$studySiteCenter[1],
                        lat = values$studySiteCenter[2],
                        popup = "Study Site",
                        options = markerOptions(riseOnHover = TRUE),
                        layerId = "studySiteMarker"
                    )
            }
        }
    })
    
    # Observe when the computeMatrix button is clicked
    observeEvent(input$computeMatrix, {
        req(values$raster)
        req(values$shape)
        req(input$classColumn)
        
        # Extract the raster values at the locations of the shapefile features
        tryCatch({
            # Ensure shapefile has the class column
            if (!(input$classColumn %in% colnames(values$shape))) {
                showNotification("Selected class column does not exist in shapefile.", type = "error")
                return(NULL)
            }
            wgs84_crs <- "+proj=longlat +datum=WGS84 +no_defs"
            
            # Convert sf object to terra SpatVector and project to WGS84
            shape_vect <- terra::vect(values$shape) %>% 
                terra::project(wgs84_crs)
            
            # Rasterize the shapefile using the selected class column
            x <- terra::rasterize(shape_vect, values$raster, field = input$classColumn)
            
            # Stack the raster and the rasterized shapefile
            stk <- terra::rast(list(values$raster, x))
            names(stk) <- c("pred", "true")
            
            # Convert to data frame and filter out NA values
            stk_df <- as.data.frame(stk, xy = FALSE, na.rm = TRUE) %>% 
                dplyr::filter(!is.na(true)) %>% 
                mutate(
                    pred = as.integer(pred),
                    true = as.integer(true)
                )
            
            # Create confusion matrix
            cm <- caret::confusionMatrix(
                factor(stk_df$pred),
                factor(stk_df$true)
            )
            
            values$confusionMatrix <- cm
            
            # Extract overall accuracy
            overall_accuracy <- cm$overall['Accuracy']
            
            # Display the confusion matrix and overall accuracy in a modal dialog
            showModal(modalDialog(
                title = "Confusion Matrix",
                size = "l",  # Large modal dialog
                easyClose = TRUE,
                footer = tagList(
                    downloadButton("downloadMatrix", "Download Confusion Matrix"),
                    modalButton("Close")
                ),
                tagList(
                    h4(sprintf("Overall Accuracy: %.2f%%", overall_accuracy * 100)),
                    br(),
                    renderTable({
                        # Convert the confusion matrix to a table with proper labels
                        cm_table <- cm$table
                        # Convert to data frame for display
                        cm_df <- as.data.frame.matrix(cm_table)
                        # Add 'Prediction' as a column
                        cm_df <- cbind(Prediction = rownames(cm_df), cm_df)
                        # Reset row names
                        rownames(cm_df) <- NULL
                        cm_df
                    }, rownames = FALSE)
                )
            ))
            
            # Enable the downloadMatrix button
            # (No need to enable/disable here since it's inside the modal dialog)
            
        }, error = function(e) {
            showNotification(paste("Error computing validation matrix", e$message), type = "error")
        })
    })
    
    # Download handler for the confusion matrix
    output$downloadMatrix <- downloadHandler(
        filename = function() {
            paste("confusion_matrix_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            req(values$confusionMatrix)
            # Write the confusion matrix to a CSV file in the desired format
            cm_table <- values$confusionMatrix$table
            cm_df <- as.data.frame.matrix(cm_table)
            cm_df <- cbind(Prediction = rownames(cm_df), cm_df)
            rownames(cm_df) <- NULL
            write.csv(cm_df, file, row.names = FALSE)
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
