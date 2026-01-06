# Script Name: 03_st_30m_daily_mosaic.R
# Description: Mosaics spatially tiled soil temperature GeoTIFFs exported from Google
#              Earth Engine into unified quarterly datasets.
# Author: Haoze Zhang
# Requirements: R (>= 4.0.0), terra package.
# License: MIT License


# Import necessary libraries
library(terra)  # For spatial data handling and raster operations

# Define directories
# CRITICAL NOTE: Please modify this path to match your local environment.
# Input directory containing spatially tiled GeoTIFFs
str_dir_input <- "./ST_tiles_central_TP_2001_2024"
# Output directory for the merged quarterly raster mosaics
str_dir_output <- "./ST_central_TP_2001_2024"

# Ensure the output directory exists
if (!dir.exists(str_dir_output)) {
  dir.create(str_dir_output, recursive = TRUE)
}

# Define the temporal range for processing
num_years <- 2001:2024
str_quarters <- c("Q1", "Q2", "Q3", "Q4")

# Batch processing loop
# Outer loop: Iterate through years
for (num_year in num_years) {
  
  # Inner loop: Iterate through quarters
  for (str_quarter in str_quarters) {
    
    # Construct a regular expression pattern to identify all constituent tiles for
    # the specific year and quarter (e.g., matches 'img_st_2001_Q1_tile_0_0.tif')
    str_pattern_file <- sprintf("^img_st_%d_%s_tile_.*\\.tif$", num_year, str_quarter)
    
    # Retrieve the list of matching file paths
    str_tile_files <- sort(list.files(path = str_dir_input, pattern = str_pattern_file, full.names = TRUE))
    
    # Validation: Verify data completeness by ensuring all four constituent tiles are present
    # This prevents the creation of spatially incomplete mosaics due to missing upstream files
    if (length(str_tile_files) != 4) {
      warning(
        sprintf(
          "Incomplete tile set for Year: %d, Quarter: %s. Expected 4, found %d. Skipping...",
          num_year, str_quarter, length(str_tile_files)
        )
      )
      next
    }
    
    # Define the output filename
    str_filename_output <- file.path(str_dir_output, sprintf("ST_%d_%s.tif", num_year, str_quarter))
    
    # Initialize a Virtual Raster (VRT) to spatially mosaic the tiles
    # VRTs allow processing of large datasets without loading them entirely into memory
    vrt_lmfc_tiles <- vrt(str_tile_files)
    
    # Restore band names from the first source tile to ensure consistency
    # (VRT creation may occasionally reset band names to generic defaults)
    names(vrt_lmfc_tiles) <- names(rast(str_tile_files[[1]]))
    
    # Write the mosaicked raster to disk with optimized compression settings
    writeRaster(
      vrt_lmfc_tiles, 
      filename = str_filename_output, 
      datatype = "INT2S", # Signed 16-bit Integer
      NAflag = -32768,    # Explicitly define NoData value for downstream compatibility
      gdal = c(
        "COMPRESS=LZW",         # Lossless compression
        "PREDICTOR=2",          # Horizontal differencing predictor (highly efficient for Int16)
        "TILED=YES",            # Enable internal tiling for efficient random access
        "BIGTIFF=YES",          # Support files larger than 4GB
        "NUM_THREADS=ALL_CPUS"  # Utilize multi-threading for compression
      ),
      overwrite = TRUE
    )
    
    # Message for progress tracking
    message(sprintf("Successfully processed: %s", basename(str_filename_output)))
    
    # Garbage collection to free up memory
    rm(vrt_lmfc_tiles) 
    gc()
  }
}
