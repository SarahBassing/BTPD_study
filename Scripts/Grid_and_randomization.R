  #'  ------------------------------
  #'  Colony grid and randomization
  #'  BTPD study
  #'  Montana State University
  #'  May 2026
  #'  ------------------------------

  ### NEW POINTS FOR 1 & 5  


  #'  Clear workspace
  rm(list = ls())  

  #'  Load libraries
  library(sf)
  library(terra)
  library(units)
  library(spsurvey)
  library(mapview)
  library(ggplot2)
  library(tidyverse)

  #'  Load spatial data
  pd_2021 <- st_read(dsn = "./Spatial/BRR_pdogs/2021_pdog_colonies_Shapefiles", 
                     layer = "2021_colonies")
  pd_2023 <- st_read(dsn = "./Spatial/BRR_pdogs/2023_pdogs/2023_complete_colonys_on_BRR", 
                     layer = "2023_pdog_brr_complete")
  pd_2024 <- st_read(dsn = "./Spatial/BRR_pdogs/2024", layer = "2024_incomplete1") 
  
  #'  Clean up colony polygon data
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
  
  #'  2026 mapped colonies
  st_layers("./Spatial/BRR_pdogs/Tracks_num.gpx")
  pd_2026 <- st_read("./Spatial/BRR_pdogs/Tracks_num.gpx", layer = "tracks") %>%
    dplyr::select(c(name, geometry)) %>%
    #'  Assign same colony ID to colonies with multiple polygons
    mutate(colony_id = name,
           colony_id = ifelse(name == "11-1", 11, colony_id),
           colony_id = ifelse(name == "16-1" | name == "16-2" | name == "16-4", 16, colony_id),
           colony_id = ifelse(name == "17-1", 17, colony_id),
           colony_id = ifelse(name == "23-1", 23, colony_id),
           colony_id = ifelse(name == "8-2", 8, colony_id),
           colony_id = ifelse(name == "7-active", 7, colony_id),
           colony_id = ifelse(name == "7", "7-inactive", colony_id)) %>%
    # rename("colony_id" = "name") %>%
    #'  Ensure projection is consistent - EPSG:32614 is WGS 84 / UTM Zone 14N
    st_transform(., crs = "EPSG:32615") 
  #'  Fix topology issues that arise from creating geometries with GPS tracking data
  pd_2026_fix <- pd_2026 %>%
    st_make_valid(.) 
  #'  Review what types of geometries are included
  table(st_geometry_type(pd_2026_fix))
  #'  Remove and tiny or degenerated bits
  pd_2026_fix <- pd_2026_fix[!st_is_empty(pd_2026_fix), ]
  #'  Convert to polygons
  pd_2026_fix <- pd_2026_fix %>%
    st_cast(., "POLYGON") %>%
    #'  Calculate area
    mutate(acres = as.numeric(set_units(st_area(.), "acres")),
           area_km2 = as.numeric(set_units(st_area(.), "km^2"))) %>%
    relocate(geometry, .after = area_km2)
  mapview(pd_2026_fix, zcol = "colony_id")
  
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
    st_buffer(dist = -gap_distance)
  #'  Review colony polygons before and after union
  mapview(pd_2026_fix)
  mapview(pd_2026_union, zcol = "colony_id")
  
  #'  List colony sf objects
  pd_colonies <- list(pd_2023, pd_2024, pd_2026_union) 
  
  #'  List 2026 colony IDs
  print(pd_colonies[[3]]$colony_id)
  
  #'  Visualize
  mapview::mapview(pd_colonies[[1]], zcol = "colony_name")
  mapview::mapview(pd_colonies[[2]], zcol = "colony_name")
  mapview::mapview(pd_colonies[[3]], zcol = "colony_id")
  #'  Compare current colonies with previously mapped colonies
  mapview::mapview(pd_colonies[[2]]) + mapview::mapview(pd_colonies[[3]], zcol = "colony_id")
  
  #'  Define spatial extent of colonies (currently using 2024 data b/c more colonies mapped than 2026)
  st_bbox(pd_2024)  # WGS 84
  st_bbox(pd_colonies[[2]]) # WGS 84 / UTM 14N
  
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
    #'  Generate data frame with grid cell centroid coordinates and ID numbers
    grid_cells <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
    head(grid_cells)
    
    out <- list(r_mask, grid_cells)
    
    return(out)
  }
  #'  Generate grids for each year of data, remember that buff and grid_size are in meters
  pd_colony_grid_2023 <- colony_grid_id(pd_colonies[[1]], buff = 25, grid_size = 1)
  pd_colony_grid_2024 <- colony_grid_id(pd_colonies[[2]], buff = 25, grid_size = 1)
  pd_colony_grid_2026 <- colony_grid_id(pd_colonies[[3]], buff = 25, grid_size = 1)
  
  writeRaster(pd_colony_grid_2023[[1]], "./Outputs/colony_2023_grids_1m.tif", overwrite = TRUE)
  save(pd_colony_grid_2023, file = "./Outputs/colony_2023_grid_1m_values.RData")
  writeRaster(pd_colony_grid_2024[[1]], "./Outputs/colony_2024_grids_1m.tif", overwrite = TRUE)
  save(pd_colony_grid_2024, file = "./Outputs/colony_2024_grid_1m_values.RData")
  writeRaster(pd_colony_grid_2026[[1]], "./Outputs/colony_2026_grids_1m.tif", overwrite = TRUE)
  save(pd_colony_grid_2026, file = "./Outputs/colony_2026_grid_1m_values.RData")
  
  #'  Visualize the raster stack
  plot(pd_colony_grid_2024[[1]]$cell_id)    # unique grid cells
  plot(pd_colony_grid_2024[[1]]$colony_id)  # unique colony polygons
  plot(pd_colony_grid_2024[[1]]$group_id)   # unique colony groups (colonies within 50 m of each other)
  
  #'  Review output data frame
  head(pd_colony_grid_2024[[2]])
  
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
  
  #'  -------------------------------
  ####  Treatment and Control sites  ####
  #'  -------------------------------
  #'  Append colony attributes to spatial data
  set.seed(963)
  #'  "Legacy" colonies to have raptor perches (per ranch request): 10, 15, 26, 29?
  colony_dat <- pd_colonies[[3]] %>%
    mutate(legacy_colony = ifelse(colony_id == 10 | colony_id == 15 | 
                                    colony_id == 26, "legacy_site", NA)) %>% 
    #'  Remove portion of colony 7 that is inactive (as of July 2026)
    filter(colony_id != "7-inactive") %>%
    #'  Ensure spatial data for polygons are valid for next step
    st_make_valid() %>%
    #'  Calculate total acreage for colonies with multiple polygons
    group_by(colony_id) %>%
    mutate(acres_total = sum(acres),
           #'  Currently making up which colonies have a fence line until real data is added
           fence_present = sample(c("Yes", "No"), size = 1, prob = c(0.6, 0.4))) %>% #### REMOVE THIS ONCE REAL FENCE LINE DATA ARE ADDED
    ungroup() %>%
    #'  Move new column to right after the polygon-specific acreage
    relocate(acres_total, .after = acres) 
    
  #'  Classify colony size into small, medium, and large based on acreage
  summary(colony_dat$acres_total)
  hist(colony_dat$acres_total, breaks = 20)
  colony_dat <- colony_dat %>%
    mutate(colony_size = ifelse(acres_total < 15, "Small", "Large"),
           colony_size = ifelse(acres_total >= 15 & acres_total < 60, "Medium", colony_size)) %>%
    relocate(geometry, .after = colony_size) %>%
    arrange(colony_size)
  table(colony_dat$colony_size)
  
  #' #'  Make up an extra large colony with no fence line so there's at least one 
  #' #'  colony for each block combo in the data set for now
  #' extra_colony <- colony_dat[colony_dat$acres > 90,]
  #' extra_colony$colony_id <- "100"
  #' extra_colony$fence_present <- "Yes"
  #' colony_dat <- bind_rows(colony_dat, extra_colony)
  #' ####  GET RID OF THIS ONCE A 4TH LARGE SITE IS MAPPED  ####
  
  #'  Grab centroid of each polygon (needed for GRTS step below)
  colony_centroids_sf <- st_centroid(colony_dat)
  #'  Reduce data set to centroid of largest polygon for each colony_id (gets rid 
  #'  of multiple centroids for colonies with multiple polygons)
  colony_centroids_reduced <- colony_centroids_sf %>%
    group_by(colony_id) %>%
    filter(acres == max(acres)) %>%
    ungroup()
  
  #'  -----------------------------
  #####  Randomly select colonies  #####
  #'  -----------------------------
  #'  Target sample: 2 sites per stratum based on full block factorial design for
  #'  colony size [S, M, L] and presence of fence line [Y/N] = 3 * 2 * 2
  #'  --> 6 treatment and 6 control sites
  #'  -----------------------------
  #'  Set seed for reproducibility
  # set.seed(12498) 
  set.seed(69428)
  #'  Generalized random tessellation stratified (GRTS) sample with legacy sites
  sample_set <- spsurvey::grts(
    sframe = colony_centroids_reduced,             # finite set of spatial points to choose from
    n_base = c(Small = 4, Medium = 4, Large = 4),  # sample size for each strata
    stratum_var = "colony_size",                   # column name of strata
    caty_var = "fence_present",                    # column name of second strata with unequal sampling probability (even though in this case it is equal)
    caty_n = c(Yes = 2, No = 2),                   # vector of sample sizes per second strata (equal sampling probability in this case)
    legacy_var = "legacy_colony",                  # column name of legacy indicator (non-legacy colonies are NA)
    # n_over = 4, # Must have sufficient colonies in each block for this to work -- Large colonies are an issue here
    DesignID = "Colony_id"                         # naming structure for each site's identifier selected in the sample
  )
  
  #'  Selected colonies
  legacy_centroids <- sample_set$sites_legacy
  random_centroids <- sample_set$sites_base
  experiment_centroids <- bind_rows(legacy_centroids, random_centroids) %>%
    dplyr::select(c(colony_id, name, acres, area_km2, colony_size, fence_present, 
                    legacy_colony, wgt, geometry)) %>%
    arrange(colony_size, fence_present)
  View(experiment_centroids)
  ####  NOTE: Double check each size class has 2 observations per fence classification!
  
  #'  Set new seed for reproducibility
  set.seed(97364)
  
  #'  Randomly select which colonies will be treatment sites, including all legacy sites 
  #'  (non-selected colonies will default as control sites)
  treatment_set <- spsurvey::grts(
    sframe = experiment_centroids,                 # finite set of spatial points to choose from
    n_base = c(Small = 2, Medium = 2, Large = 2),  # sample size for each strata
    stratum_var = "colony_size",                   # column name of strata
    caty_var = "fence_present",                    # column name of second strata with unequal sampling probability (even though in this case it is equal)
    caty_n = c(Yes = 1, No = 1),                   # vector of sample sizes per second strata (equal sampling probability in this case)
    legacy_var = "legacy_colony",                  # column name of legacy indicator (non-legacy colonies are NA)
    DesignID = "Colony_id"                         # naming structure for each site's identifier selected in the sample
  )
  treatment_legacy <- treatment_set$sites_legacy
  treatment_random <- treatment_set$sites_base
  treatment_centroids <- bind_rows(treatment_legacy, treatment_random) %>%
    dplyr::select(c(colony_id, name, acres, area_km2, colony_size, fence_present, 
                    legacy_colony, wgt, geometry)) %>%
    arrange(colony_size, fence_present)
  View(treatment_centroids)
  ####  NOTE: Double check each size class has 1 observation per fence classification!
  
  #'  Save polygons of randomly selected colonies
  #'  Note: Expect more rows in experiment_colonies than in experiment_centroids
  #'  because there are multiple polygons for some colonies with the same colony_id
  #'  in colony_dat. 
  experiment_colonies <- colony_dat[colony_dat$colony_id %in% experiment_centroids$colony_id,] %>%
    #'  Create placeholder column to indicate which is a treatment vs control colony
    mutate(treatment = FALSE) %>%
    arrange(colony_size, fence_present)
  View(experiment_colonies)
  mapview(experiment_colonies, zcol = "colony_id")
  
  #' Indicate which colonies are treatment vs control 
  experiment_colonies <- experiment_colonies %>%
    mutate(treatment = ifelse(colony_id %in% treatment_centroids$colony_id, TRUE, treatment)) %>%
    relocate(treatment, .before = geometry)
  
  #'  Double check each size-fence block has a treatment and control colony pair
  tst <- experiment_colonies %>%
    group_by(colony_id) %>%
    slice(1L) %>%
    ungroup() %>%
    arrange(colony_size, fence_present, treatment)
  View(tst)
  
  #'  Visualize distribution of treatment and control colonies
  mapview(experiment_colonies, zcol = "treatment")
  
  
  
  
  
  