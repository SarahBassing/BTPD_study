  #'  ------------------------------
  #'  Colony grid and randomization
  #'  BTPD study
  #'  Montana State University
  #'  May 2026
  #'  ------------------------------
  
  #'  Clear workspace
  rm(list = ls())  

  #'  Load libraries
  library(sf)
  library(terra)
  library(mapview)
  library(ggplot2)
  library(tidyverse)

  #'  Load spatial data
  #'  Colony polygons: Projected CRS: WGS 84 / UTM zone 14N
  pd_2021 <- st_read(dsn = "./Spatial/BRR_pdogs/2021_pdog_colonies_Shapefiles", layer = "2021_colonies") %>%
    st_zm(., drop = TRUE, what = "ZM")
  pd_2023 <- st_read(dsn = "./Spatial/BRR_pdogs/2023_pdogs/2023_complete_colonys_on_BRR", 
                     layer = "2023_pdog_brr_complete") %>%
    st_zm(., drop = TRUE, what = "ZM")
  pd_2024 <- st_read(dsn = "./Spatial/BRR_pdogs/2024", layer = "2024_incomplete1") %>%
    st_zm(., drop = TRUE, what = "ZM")
  
  # pd_colony_gdb <- "./Spatial/BRR_pdogs/Prairie_dog_colonies_share_files.gdb"
  # gdb_info <- terra::rast(pd_colony_gdb) 
  # gdb_info <- st_layers(pd_colony_gdb) 
  # print(gdb_info)
  # # rds <- sf::st_read(dsn = gdb_info, layer = "Road_OpenStreetMap")
  
  #'  Define coordinate system and reproject
  new_proj <- function(sf_obj, new_crs) {
    # print(crs(sf_obj))
    sf_obj_new_proj <- st_transform(sf_obj, crs = new_crs) 
    return(sf_obj_new_proj)
  }
  pd_colony_list <- list(pd_2023, pd_2024) #pd_2021, 
  pd_colonies <- lapply(pd_colony_list, new_proj, new_crs = "EPSG:32615") # EPSG:32614 is WGS 84 / UTM Zone 14N
  
  #'  Visualize
  mapview::mapview(pd_colonies[[1]], zcol = "ident")
  mapview::mapview(pd_colonies[[2]], zcol = "ident")
  # mapview::mapview(pd_colonies[[3]], zcol = "ident")

  #'  List colony IDs
  print(pd_colonies[[2]]$ident)
  
  #' #'  View clusters of colonies
  #' #'  Conservation colonies
  #' pd_2024_acra <- pd_2024[pd_2024$ident == "ACRA1" | pd_2024$ident == "Acra2" | pd_2024$ident == "ACRA3",]
  #' mapview(pd_2024_acra, zcol = "ident")  
  #' #'  Cooper colonies
  #' pd_2024_cooper <- pd_2024[pd_2024$ident == "Cooper #1 Main" | pd_2024$ident == "Cooper #2 East" | pd_2024$ident == "Cooper House",]
  #' mapview(pd_2024_cooper, zcol = "ident")  
  #' #'  Problem colonies in NE corner
  #' pd_2024_NE <- pd_2024[pd_2024$ident == "T950" | pd_2024$ident == "T1268",]
  #' mapview(pd_2024_NE, zcol = "ident")
  
  #'  Define spatial extent of colonies
  st_bbox(pd_2024)  # WGS 84
  st_bbox(pd_colonies[[2]]) # WGS 84 / UTM 14N
  
  #'  Function to generate and lable grid cells based on colony and colony groupings
  colony_grid_id <- function(sf_poly, buff, grid_size) {
    #'  Clean up colony polygon data
    poly <- sf_poly %>%
      select(ident) %>%
      mutate(colony_id = 1:nrow(.)) %>%
      relocate(colony_id, .before = ident) %>%
      rename("colony_name" = "ident") 
    
    #'  Buffer each polygon and create unions where buffered colonies intersect
    poly_buff <- st_buffer(poly, dist = buff) # in meters (buff * 2 = total gap between colonies)
    buff_union <- st_union(poly_buff) %>%
      st_cast(., "POLYGON") # split multipolygons into individual clusters
    
    #'  Convert clusters to sf object with unique group ID
    clusters <- st_sf(
      group_id = 1:length(buff_union),
      geometry = buff_union)
    
    #'  Assign group ID to original colony polygon with spatial join
    poly_groups <- st_join(poly, clusters, join = st_intersects)
    #'  Keep only the first instance of when clusters intersect each other (rare)
    poly_groups <- poly_groups[!duplicated(poly_groups$colony_id),]
    
    #'  Rasterize poly and poly_groups as separate layers
    poly_groups_vect <- vect(poly_groups)
    r_template <- rast(poly_groups_vect, resolution = grid_size) # resolution in m
    
    #'  Sequential cell id  layer
    r_cells <- r_template
    values(r_cells) <- 1:ncell(r_cells)
    names(r_cells) <- "cell_id"
    print(r_cells)
    
    #'  Create layer where each cell gets the ID of the colony polygon it overlaps
    r_colony_id <- rasterize(poly_groups_vect, r_template, field = "colony_id")
    names(r_colony_id) <- "colony_id"
    print(r_colony_id)
    #'  Create layer where each cell gets the ID of the grouped polygons it overlaps
    r_group_id <- rasterize(poly_groups_vect, r_template, field = "group_id")
    names(r_group_id) <- "group_id"
    print(r_group_id)
    
    #'  Stack all layers and mask to polygons
    r_stack <- c(r_cells, r_colony_id, r_group_id)
    r_mask <- mask(r_stack, poly_groups_vect)
    
    #'  Inspect
    print(r_mask)
    plot(r_mask)
    #'  Genearte data frame with grid cell centroid coordinates and ID numbers
    grid_cells <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
    head(grid_cells)
    
    out <- list(r_mask, grid_cells)
    
    return(out)
  }
  #' #'  Generate grids for each year of data, remember that buff and grid_size are in meters
  #' pd_colony_grid_2023 <- colony_grid_id(pd_colonies[[1]], buff = 25, grid_size = 1)
  #' pd_colony_grid_2024 <- colony_grid_id(pd_colonies[[2]], buff = 25, grid_size = 1)
  #' 
  #' writeRaster(pd_colony_grid_2023[[1]], "./Outputs/colony_2023_grids_1m.tif", overwrite = TRUE)
  #' save(pd_colony_grid_2023, file = "./Outputs/colony_2023_grid_1m_values.RData")
  #' writeRaster(pd_colony_grid_2024[[1]], "./Outputs/colony_2024_grids_1m.tif", overwrite = TRUE)
  #' save(pd_colony_grid_2024, file = "./Outputs/colony_2024_grid_1m_values.RData")
  #' 
  #' #'  Visualize the raster stack
  #' plot(pd_colony_grid_2024[[1]]$cell_id)    # unique grid cells
  #' plot(pd_colony_grid_2024[[1]]$colony_id)  # unique colony polygons
  #' plot(pd_colony_grid_2024[[1]]$group_id)   # unique colony groups (colonies within 50 m of each other)
  #' 
  #' #'  Review output data frame
  #' head(pd_colony_grid_2024[[2]])
  
  #'  Load raster and grid data so you don't have to remake these
  pd_colony_grid_2023 <- rast("./Outputs/colony_2023_grids_1m.tif")
  load("./Outputs/colony_2023_grid_1m_values.RData")
  pd_colony_grid_2024 <- rast("./Outputs/colony_2024_grids_1m.tif")
  load("./Outputs/colony_2024_grid_1m_values.RData")
  
  #'  Randomly sample each colony grid
  sample_grid <- function(dat, n) {
    site_select <- dat %>%
      group_by(group_id) %>%
      slice_sample(n = n) %>%
      ungroup() %>%
      st_as_sf(coords = c("x", "y"), crs = "EPSG:32615") %>%
      st_transform(crs = crs(pd_2024)) %>%
      mutate(Long = st_coordinates(.)[,1],
             Lat = st_coordinates(.)[,2]) %>%
      st_drop_geometry()
    return(site_select)
  }
  set.seed(5262026)
  random_cells_2023 <- sample_grid(pd_colony_grid_2023[[2]], n = 30) # over sampling to allow for flexibility in the field
  random_cells_2024 <- sample_grid(pd_colony_grid_2024[[2]], n = 30)
  
  #'  Save
  write_csv(random_cells_2023, file = "./Outputs/random_cells_2023_colonies.csv") 
  write_csv(random_cells_2024, file = "./Outputs/random_cells_2024_colonies.csv") 
  
  #'  Visualize random sites
  random_cells_2024_sf <- st_as_sf(random_cells_2024, coords = c("Long", "Lat"), crs = crs(pd_2024))
  plot(random_cells_2024_sf[1])
  mapview(random_cells_2024_sf, zcol = "cell_id")
  
  #'  Convert spatial data to .GPX files for OnX and Garmin GPS units
  sample_pts <- random_cells_2024_sf %>%
    #' #'  Relabel cell ID column
    #' rename("ID" = "cell_id") %>%
    #' Transform projection to WGS 84 (lat/long)
    st_transform(., crs = "EPSG:4326")
  #'  Export as .GPX file format
  write_sf(sample_pts, file = paste0("random_points_", Sys.Date(), ".gpx"), dataset_options = "GPX_USE_EXTENSIONS=YES")
  
  #'  Convert colony polygons to .KML flie for OnX and Garmin GPS units
  colony_polygons <- pd_2024 %>%
    st_transform(., crs = "EPSG:4326")
  write_sf(colony_polygons, "colony_polygons_2024.kml", dataset_options = "GPX_USE_EXTENSIONS=YES")
  
  
  