  #'  ------------------------------
  #'  Colony grid and randomization
  #'  BTPD study
  #'  Montana State University
  #'  May 2026
  #'  ------------------------------
  #'  Script to guide randomization process for trapping BTPD and conducting raptor
  #'  perch experiment. Includes formatting historic and current colony spatial data, 
  #'  generating and randomly from a grid each colony polygon to guide randomization
  #'  of trap placement, and randomly select treatment and control colonies for
  #'  experiment while allowing for "legacy" colonies to be used as treatment sites.
  #'  ------------------------------
  
  #'  Clear work space
  rm(list = ls())  

  #'  Load libraries
  library(sf)
  library(terra)
  library(nngeo)
  library(units)
  library(spsurvey)
  library(mapview)
  library(ggplot2)
  library(tidyverse)

  #'  Load historic spatial data
  pd_2021 <- st_read(dsn = "./Spatial/BRR_pdogs/2021_pdog_colonies_Shapefiles", 
                     layer = "2021_colonies")
  pd_2023 <- st_read(dsn = "./Spatial/BRR_pdogs/2023_pdogs/2023_complete_colonys_on_BRR", 
                     layer = "2023_pdog_brr_complete")
  pd_2024 <- st_read(dsn = "./Spatial/BRR_pdogs/2024", layer = "2024_incomplete1") 
  
  #'  Review data layers in .gpx file
  st_layers("./Spatial/BRR_pdogs/2026/Tracks_num.gpx")  
  st_layers("./Spatial/BRR_pdogs/2026/Colonies.gpx")
  # st_layers("./Spatial/BRR_pdogs/2026/Selected Data from Colonies.gpx")
  #'  Load tracks layer for mapping 2026 colonies
  pd_2026 <- st_read("./Spatial/BRR_pdogs/2026/Tracks_num.gpx", layer = "tracks") %>%
    dplyr::select(c(name, geometry))
  pd_2026_final <- st_read("./Spatial/BRR_pdogs/2026/Colonies.gpx", layer = "tracks") %>%
    dplyr::select(c(name, geometry))
  pd_2026_col19 <- st_read("./Spatial/BRR_pdogs/2026/Colony 19.gpx", layer = "tracks") %>%
    dplyr::select(c(name, geometry))
  # pd_2026_col35and36 <- st_read("./Spatial/BRR_pdogs/2026/Selected Data from Colonies.gpx", layer = "tracks") %>%
  #   dplyr::select(c(name, geometry))
  
  #'  Grab geographic coordinate system of pd_2026 (WGS84)
  wgs84 <- crs(pd_2026_final)
  
  #'  --------------------------
  ####  Format colony polygons  ####
  #'  --------------------------
  #'  Clean up historic colony polygon data
  colony_name_clean <- function(sf_poly) {
    poly <- sf_poly %>%
      #'  Remove extra dimension of spatial data
      st_zm(., drop = TRUE, what = "ZM") %>%
      #'  Make sure projection is consistent: WGS 84 / UTM zone 14N
      st_transform(., crs = "EPSG:32615") %>%
      #'  Grab identification info
      select(ident) %>%
      #'  Create new colony ID name
      mutate(colony_id = 1:nrow(.)) %>%
      #'  Reorganize columns
      relocate(colony_id, .before = ident) %>%
      #'  Rename column to something a little more meaningful
      rename("colony_name" = "ident") 
    return(poly)
  }
  pd_2023 <- colony_name_clean(pd_2023)
  pd_2024 <- colony_name_clean(pd_2024)
  
  #'  Reformat mapped 2026 tracks and convert to sf polygons
  pd_2026 <- pd_2026 %>%
    #'  Assign same colony ID to colonies with multiple polygons
    mutate(colony_id = name,
           colony_id = ifelse(name == "11-1", 11, colony_id),
           colony_id = ifelse(name == "16-1" | name == "16-2" | name == "16-4", 16, colony_id),
           colony_id = ifelse(name == "17-1", 17, colony_id),
           colony_id = ifelse(name == "23-1", 23, colony_id),
           colony_id = ifelse(name == "8-2", 8, colony_id),
           colony_id = ifelse(name == "7-active", 7, colony_id),
           colony_id = ifelse(name == "7", "7-inactive", colony_id)) %>%
    #'  Re-project to be consistent with historic data - EPSG:32614 is WGS 84 / UTM Zone 14N
    st_transform(., crs = "EPSG:32615") 
  pd_2026_final <- pd_2026_final %>%
    #'  Append colony 19 tracks to larger 2026 dataset
    bind_rows(pd_2026_col19) %>%
    #'  Assign same colony ID to colonies with multiple polygons
    mutate(colony = name,
           colony = ifelse(grepl("2-", name), "Colony 2", colony),
           colony = ifelse(grepl("8-", name), "Colony 8", colony),
           colony = ifelse(grepl("10-", name), "Colony 10", colony),
           colony = ifelse(grepl("11-", name), "Colony 11", colony),
           colony = ifelse(grepl("12-", name), "Colony 12", colony),
           colony = ifelse(grepl("17-", name), "Colony 17", colony),
           colony = ifelse(grepl("21-", name), "Colony 21", colony),
           colony = ifelse(grepl("22-", name), "Colony 22", colony),
           colony = ifelse(grepl("23-", name), "Colony 23", colony),
           colony = ifelse(grepl("24-", name), "Colony 24", colony),
           colony = ifelse(grepl("25-", name), "Colony 25", colony),
           colony = ifelse(grepl("26-", name), "Colony 26", colony), # NEED TO DOUBLE CHECK THIS ONE
           colony = ifelse(grepl("27-", name), "Colony 27", colony),
           colony = ifelse(grepl("29-", name), "Colony 29", colony),
           colony = ifelse(grepl("35-", name), "Colony 35", colony),
           # colony = ifelse(name == "Coony 26-3", "Colony 26???", colony),
           colony_id = substring(colony, 8),
           colony_id = as.numeric(colony_id)) %>%
    #'  Re-project to be consistent with historic data - EPSG:32614 is WGS 84 / UTM Zone 14N
    st_transform(., crs = "EPSG:32615") 
  #' pd_2026_35and36 <- pd_2026_col35and36 %>%
  #'   #'  Assign same colony ID to colonies with multiple polygons
  #'   mutate(colony = name,
  #'          colony = ifelse(grepl("35-", name), "Colony 35", colony),
  #'          colony_id = substring(colony, 8),
  #'          colony_id = as.numeric(colony_id)) %>%
  #'   #'  Re-project to be consistent with historic data - EPSG:32614 is WGS 84 / UTM Zone 14N
  #'   st_transform(., crs = "EPSG:32615") 
  #'  Fix topology issues that arise from creating geometries with GPS tracking data
  pd_2026_fix <- pd_2026 %>%
    st_make_valid(.) 
  pd_2026_final_fix <- pd_2026_final %>%
    st_make_valid(.) 
  # pd_2026_35and36_fix <- pd_2026_35and36 %>%
  #   st_make_valid(.) 
  #'  Review what types of geometries are included
  table(st_geometry_type(pd_2026_fix))          # should be all MULTILINESTRING
  table(st_geometry_type(pd_2026_final_fix))    # should be all MULTILINESTRING
  #'  Remove any tiny or degenerated bits 
  pd_2026_fix <- pd_2026_fix[!st_is_empty(pd_2026_fix), ]
  pd_2026_final_fix <- pd_2026_final_fix[!st_is_empty(pd_2026_final_fix), ]
  # pd_2026_35and36_fix <- pd_2026_35and36_fix[!st_is_empty(pd_2026_35and36_fix), ]
  #'  Convert to polygons
  pd_2026_fix <- pd_2026_fix %>%
    st_cast(., "POLYGON") 
  pd_2026_final_fix <- pd_2026_final_fix %>%
    st_cast(., "POLYGON")
  # pd_2026_35and36_fix <- pd_2026_35and36_fix %>%
  #   st_cast(., "POLYGON")
  #'  Visualize colony polygons
  mapview(pd_2026_fix, zcol = "colony_id")
  mapview(pd_2026_final_fix, zcol = "colony_id")
  # mapview(pd_2026_35and36_fix, zcol = "colony_id")
  
  #'  Join slightly non-overlapping polygons from same colony to create single 
  #'  polygon per colony (typically arises when a fence line cuts through a colony)
  #'  Define gap-bridging distance
  gap_distance <- 10 # meters
  pd_2026_union <- pd_2026_fix %>%
    #'  Make sure spatial data are valid
    st_make_valid() %>%
    #'  Buffer polygons by gap_distance
    st_buffer(dist = gap_distance) %>%
    #'  Group by colony ID so only polygons from a single colony get merged
    group_by(colony_id) %>%
    #'  Union the geometries for grouped polygons
    summarise(geometry = st_union(geometry)) %>%
    #'  Remove the buffered distance
    st_buffer(dist = -gap_distance) %>%
    #'  Calculate area in acres 
    mutate(acres = as.numeric(set_units(st_area(.), "acres")),
           acres = round(acres, 4)) %>%
    relocate(acres, .before = geometry)
  pd_2026_final_union <- pd_2026_final_fix %>%
    #'  Make sure spatial data are valid
    st_make_valid() %>%
    #'  Buffer polygons by gap_distance
    st_buffer(dist = gap_distance) %>%
    #'  Group by colony ID so only polygons from a single colony get merged
    group_by(colony_id) %>%
    #'  Union the geometries for grouped polygons
    summarise(geometry = st_union(geometry)) %>%
    #'  Remove the buffered distance
    st_buffer(dist = -gap_distance) %>%
    #'  Calculate area in acres 
    mutate(acres = as.numeric(set_units(st_area(.), "acres")),
           acres = round(acres, 4)) %>%
    relocate(acres, .before = geometry)
  
  #'  Review colony polygons before and after union
  mapview(pd_2026_fix)
  mapview(pd_2026_union, zcol = "colony_id")
  mapview(pd_2026_final_fix)
  mapview(pd_2026_final_union, zcol = "colony_id")
  
  #'  Visualize historic and mapped 2026 colonies
  mapview::mapview(pd_2023, zcol = "colony_name")  
  mapview::mapview(pd_2024, zcol = "colony_name")  
  mapview::mapview(pd_2026_final_union, zcol = "colony_id")    
  #'  Compare current colonies with previously mapped colonies by overlaying them
  mapview::mapview(pd_2024) + mapview::mapview(pd_2026_final_union, zcol = "colony_id")
  
  #'  Save cleaned up 2026 colony polygons as csv
  pd_2026_clean <- as.data.frame(pd_2026_final_union) %>%
    dplyr::select(-geometry)
  write_csv(pd_2026_clean, "./Spatial/BRR_pdogs/2026/pd_2026_colony_attributes.csv")
  
  #'  Review spatial extent of colonies 
  st_bbox(pd_2026_final) # WGS 84 / UTM 14N
  
  #'  -------------------------------------------
  ####  Superimpose grid across active colonies  ####
  #'  -------------------------------------------
  #'  Function to generate and label grid cells based on colony and colony groupings
  colony_grid_id <- function(sf_poly, buff, grid_size) {
    #'  Buffer each polygon and create unions where buffered colonies intersect
    poly_buff <- st_buffer(sf_poly, dist = buff) # in meters (buff * 2 = total gap between colonies)
    buff_union <- st_union(poly_buff) %>%
      st_cast(., "POLYGON") # split multipolygons into individual clusters
    
    #'  Convert clusters to sf object with unique group ID
    clusters <- st_sf(
      group_id = 1:length(buff_union),
      geometry = buff_union)
    
    #'  Assign group ID to original colony polygon with spatial join
    poly_groups <- st_join(sf_poly, clusters, join = st_intersects)
    #'  Keep only the first instance of when clusters intersect each other (rare)
    poly_groups <- poly_groups[!duplicated(poly_groups$colony_id),]
    
    #'  Rasterize poly and poly_groups as separate layers
    poly_groups_vect <- vect(poly_groups)
    r_template <- rast(poly_groups_vect, resolution = grid_size) # resolution in m
    
    #'  Sequential cell id layer
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
    #'  Generate data frame with grid cell centroid coordinates and ID numbers
    grid_cells <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
    head(grid_cells)
    
    out <- list(r_mask, grid_cells)
    
    return(out)
  }

  #'  Generate grids for each year of data, remember that buff and grid_size are in meters
  #'  NOTE: This can max out R's memory - must have a lot of RAM for this to work
  #'  Only do once! Load saved grids below
  pd_colony_grid_2023 <- colony_grid_id(pd_2024, buff = 25, grid_size = 1)
  pd_colony_grid_2024 <- colony_grid_id(pd_2024, buff = 25, grid_size = 1)
  pd_colony_grid_2026 <- colony_grid_id(pd_2026_final_union, buff = 25, grid_size = 1)
  
  writeRaster(pd_colony_grid_2023[[1]], "./Spatial/Exported_data/colony_2023_grids_1m.tif", overwrite = TRUE)
  save(pd_colony_grid_2023, file = "./Spatial/Exported_data/colony_2023_grid_1m_values.RData")
  writeRaster(pd_colony_grid_2024[[1]], "./Spatial/Exported_data/colony_2024_grids_1m.tif", overwrite = TRUE)
  save(pd_colony_grid_2024, file = "./Spatial/Exported_data/colony_2024_grid_1m_values.RData")
  writeRaster(pd_colony_grid_2026[[1]], "./Spatial/Exported_data/colony_2026_grids_1m.tif", overwrite = TRUE)
  save(pd_colony_grid_2026, file = "./Spatial/Exported_data/colony_2026_grid_1m_values.RData")
  
  #'  Visualize the raster stack
  plot(pd_colony_grid_2026[[1]]$cell_id)    # unique grid cells
  plot(pd_colony_grid_2026[[1]]$colony_id)  # unique colony polygons
  plot(pd_colony_grid_2026[[1]]$group_id)   # unique colony groups (colonies within 50 m of each other)
  
  #'  Review output data frame
  head(pd_colony_grid_2026[[2]])

  #'  ---------------------------------------------------
  ####  Randomly sample colonies for trapping locations  ####
  #'  ---------------------------------------------------
  #'  Load raster and grid data so you don't have to remake these
  pd_colony_grid_2023 <- rast("./Spatial/Exported_data/colony_2023_grids_1m.tif")
  load("./Spatial/Exported_data/colony_2023_grid_1m_values.RData")
  pd_colony_grid_2024 <- rast("./Spatial/Exported_data/colony_2024_grids_1m.tif")
  load("./Spatial/Exported_data/colony_2024_grid_1m_values.RData")
  pd_colony_grid_2026 <- rast("./Spatial/Exported_data/colony_2026_grids_1m.tif")
  load("./Spatial/Exported_data/colony_2026_grid_1m_values.RData")

  #'  Randomly sample each colony grid
  sample_grid <- function(dat, n) {
    site_select <- dat %>%
      #'  Group data by group_id and draw n samples following a simple random sampling approach
      group_by(group_id) %>%
      slice_sample(n = n) %>%
      ungroup() %>%
      #'  Make data spatial
      st_as_sf(coords = c("x", "y"), crs = "EPSG:32615") %>%
      #'  Transform to WGS85 (want lat/long)
      st_transform(crs = crs(wgs84)) %>%
      #'  Create columns for Lat and Long, formatted so R knows they are coordinates
      mutate(Long = st_coordinates(.)[,1],
             Lat = st_coordinates(.)[,2]) %>%
      #'  Remove sf geometry
      st_drop_geometry()
    return(site_select)
  }
  set.seed(5262026) # For reproducibility
  random_cells_2023 <- sample_grid(pd_colony_grid_2023[[2]], n = 30) # over sampling to allow for flexibility in the field
  random_cells_2024 <- sample_grid(pd_colony_grid_2024[[2]], n = 30)
  random_cells_2026 <- sample_grid(pd_colony_grid_2026[[2]], n = 30)
  
  #'  Grab random locations from specific colonies
  random_cells_colony_newpts <- filter(random_cells_2026, colony_id == 15 | colony_id == 22 | 
                                         colony_id == 23 | colony_id == 24 | colony_id == 25 | 
                                         colony_id == 26 | colony_id == 27 | colony_id == 28 | colony_id == 34)
  random_cells_colony_35and36_newpts <- filter(random_cells_2026, colony_id == 35 | colony_id == 36)
  
  #'  Save
  write_csv(random_cells_2023, file = paste0("./Spatial/Exported_data/random_cells_2023_colonies_", Sys.Date(), ".csv")) 
  write_csv(random_cells_2024, file = paste0("./Spatial/Exported_data/random_cells_2024_colonies_", Sys.Date(), ".csv")) 
  write_csv(random_cells_2026, file = paste0("./Spatial/Exported_data/random_cells_2026_colonies_", Sys.Date(), ".csv")) 
  
  #'  Convert to sf object 
  random_cells_2024_sf <- st_as_sf(random_cells_2024, coords = c("Long", "Lat"), crs = crs(wgs84))
  random_cells_2026_sf <- st_as_sf(random_cells_2026, coords = c("Long", "Lat"), crs = crs(wgs84))
  random_cells_2026_newpts_sf <- st_as_sf(random_cells_colony_newpts, coords = c("Long", "Lat"), crs = crs(wgs84))
  random_cells_2026_35and36_newpts_sf <- st_as_sf(random_cells_colony_35and36_newpts, coords = c("Long", "Lat"), crs = crs(wgs84))
  
  #'  Visualize
  plot(random_cells_2024_sf[1])
  mapview(random_cells_2024_sf, zcol = "cell_id")
  
  plot(random_cells_2026_sf[1])
  mapview(random_cells_2026_sf, zcol = "cell_id")
  
  plot(random_cells_2026_35and36_newpts_sf[1])
  mapview(pd_2026_final_union, zcol = "colony_id") + mapview(random_cells_2026_newpts_sf)
  mapview(pd_2026_final_union, zcol = "colony_id") + mapview(random_cells_2026_35and36_newpts_sf)
  
  #'  -----------------------
  ####  Export spatial data  ####
  #'  -----------------------
  #'  Convert spatial data to .GPX files for OnX and Garmin GPS units
  convert_to_gpx <- function(cells_sf) {
    sample_pts <- cells_sf %>%
      #' Make sure projection is WGS 84 (lat/long)
      st_transform(., crs = "EPSG:4326")
    return(sample_pts)
  }
  sample_pts_2024 <- convert_to_gpx(random_cells_2024_sf)
  sample_pts_2026 <- convert_to_gpx(random_cells_2026_sf)
  sample_pts_2026_newpt <- convert_to_gpx(random_cells_2026_newpts_sf)
  sample_pts_2026_35and36_newpt <- convert_to_gpx(random_cells_2026_35and36_newpts_sf)
  
  #'  Export as .GPX file format
  write_sf(sample_pts_2024, dsn = paste0("./Spatial/Exported_data/random_points_2024_", Sys.Date(), ".gpx"), dataset_options = "GPX_USE_EXTENSIONS=YES")
  write_sf(sample_pts_2026, dsn = paste0("./Spatial/Exported_data/random_points_2026_colonies_", Sys.Date(), ".gpx"), dataset_options = "GPX_USE_EXTENSIONS=YES")
  write_sf(sample_pts_2026_newpt, dsn = paste0("./Spatial/Exported_data/random_points_2026_select_colonies_", Sys.Date(), ".gpx"), dataset_options = "GPX_USE_EXTENSIONS=YES")
  write_sf(random_cells_2026_35and36_newpts_sf, dsn = paste0("./Spatial/Exported_data/random_points_2026_colonies_35and36_", Sys.Date(), ".gpx"), dataset_options = "GPX_USE_EXTENSIONS=YES")
    
  #'  Convert colony polygons to .KML file for OnX and Garmin GPS units
  colony_polygons_2024 <- pd_2024 %>%
    st_transform(., crs = "EPSG:4326")
  colony_polygons_2026 <- pd_2026_final_union %>%
    st_transform(., crs = "EPSG:4326") %>%
    dplyr::select(-acres)
  colony_polygons_2026_12 <- pd_2026_final_union %>%
    filter(colony_id == 12) %>%
    st_remove_holes(.) %>%
    st_transform(., crs = "EPSG:4326") %>%
    dplyr::select(-acres) %>%
    st_make_valid(.)
  write_sf(colony_polygons_2024, "./Spatial/Exported_data/colony_polygons_2024.kml", dataset_options = "GPX_USE_EXTENSIONS=YES")
  write_sf(colony_polygons_2026, "./Spatial/Exported_data/colony_polygons_2026.kml", dataset_options = "GPX_USE_EXTENSIONS=YES")
  write_sf(colony_polygons_2026_12, "./Spatial/Exported_data/colony_polygons_2026_12.kml", dataset_options = "GPX_USE_EXTENSIONS=YES")
 
  #'  Save polygons as shapefiles
  st_write(pd_2026_final_union, dsn = paste0("./Spatial/BRR_pdogs/2026/BTPD_colonies_2026_", Sys.Date(), ".shp"))
  
  #'  -------------------------------
  ####  Treatment and Control sites  ####
  #'  -------------------------------
  #'  Load fence line data
  fences <- read_csv("./Data/BRR_colony_fencelines_2026.csv") %>%
    rename("colony_id" = "Colony Number") %>%
    rename("fence_present" = "Fencline Presence") %>%
    rename("colony_description" = "Colony Discription") %>%
    dplyr::select(-`Colony Size`) %>%
    mutate(fence_present = ifelse(fence_present == "?", NA, fence_present)) %>%
    mutate(fence_present = ifelse(colony_id == 40, "N", fence_present))
  
  #'  "Legacy" colonies to have raptor perches (per ranch request): 2, especially 4, and most importantly 6, 13 
  colony_dat <- pd_2026_final_union %>%
    mutate(legacy_colony = ifelse(colony_id == 2 | colony_id == 4 | 
                                    colony_id == 6 | colony_id == 13, "legacy_site", NA)) %>% 
    #' #'  Remove portion of colony 7 that is inactive (as of July 2026)
    #' filter(colony_id != "7-inactive") %>%
    #'  Ensure spatial data for polygons are valid for next step
    st_make_valid() %>%
    #'  Append colony attributes to spatial data
    full_join(fences, by = "colony_id") %>%
    #'  Remove rows with fence data but no spatial data 
    filter(!is.na(acres)) %>%
    #'  Remove rows with spatial data but no fence data
    filter(!is.na(fence_present))
    
  #'  Classify colony size into small, medium, and large based on acreage
  summary(colony_dat$acres)
  hist(colony_dat$acres, breaks = 50)
  colony_dat <- colony_dat %>%
    #'  Small = 0 - 20 acres, Medium = 20 - 60 acres, Large = 60 - 100 acres
    #'  In reality, Medium = ~20-40 acres and Large = ~70-95 acres
    #'  All size classes capture about a 20-acre range in colony sizes
    mutate(colony_size = ifelse(acres < 20, "Small", "Large"),
           colony_size = ifelse(acres >= 20 & acres < 60, "Medium", colony_size)) %>%
    arrange(colony_size, acres) %>%
    relocate(colony_size, .after = acres) %>%
    relocate(geometry, .after = colony_description) 
  table(colony_dat$colony_size)
  
  #'  Turns our all medium-sized colonies have fences - need to fake a few without for
  #'  GRTS randomization approach to work below so faking some colonies lacking fences
  #'  Set seed and randomly draw which colonies do and do not have fences
  set.seed(963)
  new_fence <- colony_dat %>%
    dplyr::select(c(colony_id, colony_size)) %>%
    group_by(colony_id) %>%
    mutate(fence_present = sample(c("Y", "N"), size = 1, prob = c(0.5, 0.5))) %>%
    ungroup
  #'  Append fake fence data to colony data for only medium-sized colonies
  colony_dat <- colony_dat %>%
    mutate(fence_present_fake = ifelse(colony_size == "Medium", new_fence$fence_present[new_fence$colony_size == "Medium"], fence_present)) %>%
    relocate(fence_present_fake, .after = fence_present) 
  
  #'  Grab centroid of each polygon (needed for GRTS step below)
  colony_centroids_sf <- st_centroid(colony_dat)
  
  #'  -----------------------------
  #####  Randomly select colonies  #####
  #'  -----------------------------
  #'  Target sample: 2 sites per stratum based on full block factorial design for
  #'  colony size [S, M, L] and presence of fence line [Y/N] = 3 * 2 * 2
  #'  --> 6 treatment and 6 control sites
  #'  -----------------------------
  #'  Set seed for reproducibility
  set.seed(97146) #9746 #45975
  #'  Generalized random tessellation stratified (GRTS) sample with legacy sites
  sample_set <- spsurvey::grts(
    sframe = colony_centroids_sf,                  # finite set of spatial points to choose from
    n_base = c(Small = 4, Medium = 4, Large = 4),  # sample size for each strata
    stratum_var = "colony_size",                   # column name of strata
    caty_var = "fence_present_fake",               # column name of second strata with unequal sampling probability (even though in this case it is equal)
    caty_n = c(Y = 2, N = 2),                      # vector of sample sizes per second strata (equal sampling probability in this case)
    legacy_var = "legacy_colony",                  # column name of legacy indicator (non-legacy colonies are NA)
    n_over = 2,                                    # over sample in each size strata - Must have sufficient colonies in each n_base strata for this to work 
    DesignID = "colony_id"                         # naming structure for each site's identifier selected in the sample
  )
  #'  Review named lists within sample_set
  names(sample_set)
  
  #'  Selected colonies
  legacy_centroids <- sample_set$sites_legacy
  random_centroids <- sample_set$sites_base
  oversample_centroids <- sample_set$sites_over
  experiment_centroids <- bind_rows(legacy_centroids, random_centroids, oversample_centroids) %>%
    dplyr::select(c(colony_id, acres, colony_size, fence_present, fence_present_fake, legacy_colony, 
                    wgt, siteuse, geometry)) %>%
    arrange(colony_size, fence_present_fake)
  View(experiment_centroids)
  ####  NOTE: Double check each size class has 2 observations per fence classification!
  
  #  Redraw specifically for new small colonies (many from previous draw were
  #  recently poisoned or plagued out so not very good colonies for experiment)
  set.seed(5696)   #45975 
  #'  Generalized random tessellation stratified (GRTS) sample with legacy sites
  sample_set_v2 <- spsurvey::grts(
    sframe = colony_centroids_sf,                  # finite set of spatial points to choose from
    n_base = c(Small = 4, Medium = 4, Large = 4),  # sample size for each strata
    stratum_var = "colony_size",                   # column name of strata
    caty_var = "fence_present_fake",               # column name of second strata with unequal sampling probability (even though in this case it is equal)
    caty_n = c(Y = 2, N = 2),                      # vector of sample sizes per second strata (equal sampling probability in this case)
    legacy_var = "legacy_colony",                  # column name of legacy indicator (non-legacy colonies are NA)
    n_over = 2,                                    # over sample in each size strata - Must have sufficient colonies in each n_base strata for this to work 
    DesignID = "colony_id"                         # naming structure for each site's identifier selected in the sample
  )
  #'  Review named lists within sample_set
  names(sample_set_v2)
  
  #'  Selected colonies
  legacy_centroids_v2 <- sample_set_v2$sites_legacy
  random_centroids_v2 <- sample_set_v2$sites_base
  oversample_centroids_v2 <- sample_set_v2$sites_over
  experiment_centroids_v2 <- bind_rows(legacy_centroids_v2, random_centroids_v2, oversample_centroids_v2) %>%
    dplyr::select(c(colony_id, acres, colony_size, fence_present, fence_present_fake, legacy_colony, 
                    wgt, siteuse, geometry)) %>%
    arrange(colony_size, fence_present_fake)
  View(experiment_centroids_v2)
  #'  Grab just the small colonies
  experiment_centroids_v2_Small <- experiment_centroids_v2 %>%
    filter(colony_size == "Small")
  
  #'  Combine original draw (Medium and Large colonies) with new draw (Small colonies)
  experiment_centroids_final <- experiment_centroids %>%
    filter(colony_size == "Medium" | colony_size == "Large") %>%
    bind_rows(experiment_centroids_v2_Small)
  
  #'  Set new seed for reproducibility
  set.seed(45761) #72027
  
  #'  Randomly select which colonies will be treatment sites, including all legacy sites 
  #'  (non-selected colonies will default as control sites)
  treatment_set <- spsurvey::grts(
    sframe = experiment_centroids_final,           # finite set of spatial points to choose from
    n_base = c(Small = 2, Medium = 2, Large = 2),  # sample size for each strata
    stratum_var = "colony_size",                   # column name of strata
    caty_var = "fence_present_fake",               # column name of second strata with unequal sampling probability (even though in this case it is equal)
    caty_n = c(Y = 1, N = 1),                      # vector of sample sizes per second strata (equal sampling probability in this case)
    legacy_var = "legacy_colony",                  # column name of legacy indicator (non-legacy colonies are NA)
    n_over = 1,
    DesignID = "colony_id"                         # naming structure for each site's identifier selected in the sample
  )
  treatment_legacy <- treatment_set$sites_legacy
  treatment_random <- treatment_set$sites_base
  treatment_oversample <- treatment_set$sites_over
  treatment_centroids <- bind_rows(treatment_legacy, treatment_random, treatment_oversample) %>%
    dplyr::select(c(colony_id, acres, colony_size, fence_present, fence_present_fake, legacy_colony, 
                    wgt, siteuse, geometry)) %>%
    arrange(colony_size, fence_present)
  View(treatment_centroids)
  ####  NOTE: Double check each size class has at least 1 observation per fence classification!
  
  #'  Save polygons of randomly selected colonies
  experiment_colonies <- colony_dat[colony_dat$colony_id %in% experiment_centroids_final$colony_id,] %>%
    #' Indicate which colonies are treatment vs control 
    mutate(treatment = ifelse(colony_id %in% treatment_centroids$colony_id, TRUE, FALSE)) %>%
    #'  Arrange by colony size and fence to more easily double check selection
    arrange(colony_size, acres, fence_present, treatment) %>%
    relocate(fence_present, .after = colony_size) %>%
    relocate(treatment, .after = fence_present) %>%
    dplyr::select(-fence_present_fake)
  View(experiment_colonies)
  
  #'  Visualize 
  #'  Distribution of treatment and control colonies (make sure now obvious spatial pattern)
  mapview(pd_2026_final_union) + mapview(experiment_colonies, zcol = "treatment")
  
  #'  Save selected colonies for experiment 
  #'  Shapefile
  st_write(experiment_colonies, dsn = paste0("./Spatial/Exported_data/experimental_colonies_", Sys.Date(), ".shp"))
  #'  CSV file
  experiment_colony_set <- experiment_colonies %>%
    as.data.frame(.) %>%
    dplyr::select(-geometry)
  write_csv(experiment_colony_set, file = paste0("./Data/experimental_colonies_", Sys.Date(), ".csv"))
  
  