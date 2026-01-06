# A 30-m daily soil temperature dataset for alpine grasslands in central Tibetan Plateau during 2001–2024

## Overview

This repository provides reproducible Google Earth Engine (GEE) and R workflows for generating a **30-m daily near-surface soil temperature (ST) dataset** for **alpine grasslands in the Nyainqentanglha Mountains, central Tibetan Plateau**, during **2001–2024**.

The code implements a satellite-based fusion and ground-based calibration method to produce a spatiotemporally continuous, 30-m daily mean soil temperature dataset at -5 cm depth. The method integrates MODIS and Landsat land surface temperature (LST), and in situ observations from the Nyainqentanglha Mountains.

The repository contains:
- A **rapid validation demo** that exports a small subset directly to **Earth Engine Assets**.
- A **full production workflow** that exports spatially tiled quarterly GeoTIFFs to **Google Drive** and mosaics them locally in **R**.

---

## Repository Structure

The workflow consists of three main scripts:

-   **`01_st_30m_daily_demo.js`**: A demonstration script for GEE. It runs the core algorithm on a sample year (2020) and exports the results directly to Earth Engine Assets for rapid validation and code verification.
-   **`02_st_30m_daily_tile_export.js`**: The full production script for GEE. It generates quarterly packed ST datasets for the complete time series (2001–2024), splits the study area (150 × 150 km) into four tiles (75 × 75 km), and exports tiled GeoTIFFs to Google Drive to bypass file size limits.
-   **`03_st_30m_daily_mosaic.R`**: An R script for post-processing. It iterates through years and quarters to mosaic the spatially tiled GeoTIFFs downloaded from Google Drive into unified quarterly datasets.

---

## Requirements

### GEE

-   A valid GEE account is required to run the `.js` scripts.
-   **Input Data**: The scripts automatically access the following public datasets from the GEE Data Catalog:
    -   **MODIS MOD11A1.061**: 1-km daily LST product
    -   **Landsat 4–9 Level-2 Collection 2**: 30-m LST product with 16-day revisit cycle
-   Sufficient GEE Asset quota or Google Drive storage for export outputs.

### R Environment

-   **R Version**: >= 4.0.0
-   **Required Packages**: `terra`
-   **Hardware**: A computer with sufficient RAM (\>= 16GB recommended) is advised for processing large GeoTIFF files.

---

## Study Area and Projection

-   **Study Area**: The default configuration targets a **150 × 150 km** area covering the Nyainqentanglha Mountains. The production script (`02_st_30m_daily_tile_export.js`) automatically subdivides this area into a **2 × 2** tile grid (four tiles, no gaps) to manage export sizes. Users can adapt the workflow to other regions by modifying the centroid coordinates (`gmt_point_center`) and grid size parameters in the "Execution Block" section.

-   **Projection**: An **Albers Equal Area Conic** CRS is used by default to ensure accurate spatial tiling and export. If your project requires a different projection, you can replace the CRS WKT string (`str_clt_crs_albers`) in the "Global Parameter Definitions" section.

---

## Methodology Overview

This dataset is generated using a **satellite-based fusion and ground-based calibration method**. The core processing steps implemented in the scripts include:

1.  **Satellite Data Pre-processing**:
    -   **MODIS LST**: Gap-filling of daily 1-km LST (MOD11A1) using pixel-wise temporal linear interpolation.
    -   **Landsat LST**: Cross-sensor harmonization (Landsat 4–9) and pixel quality control.
2.  **Satellite LST Fusion**:
    -   **Core assumption**: The ratio of daily LST to annual maximum LST would not differ between 1-km and 30-m scales at a location.
    -   Derives a normalized temporal shape factor from MODIS LST time series.
    -   Estimates the theoretical Landsat annual maximum LST. Outliers are removed using the Interquartile Range (IQR) method.
3.  **Soil Temperature Calibration**:
    -   Converts the fused LST to daily mean soil temperature at -5 cm.
    -   Applies a month-specific linear calibration model derived from in situ observations at 10 weather stations (4300–5500 m) in the Nyainqentanglha Mountains during 2006–2010.

**Note**: For detailed methodologies, please refer to the **Methods** section of the manuscript.

---

## Usage Instructions

### Quick Start (Rapid Validation Demo)

Use this method to quickly verify the algorithm and inspect the results in GEE without downloading large files.

1.  Copy the content of `01_st_30m_daily_demo.js` into the GEE Code Editor.
2.  Modify the export destination:
    -   Manually create a folder in Earth Engine Assets to store the generated soil temperature dataset (e.g., `ST_central_TP_2001_2024`).
    -   Update the `str_clt_asset_folder` variable in the "Global Parameter Definitions" section to match this path.
3.  Modify the temporal range:
    -   The script defaults to generating soil temperature data for the year **2020** for demonstration.
    -   To generate soil temperature for other time periods, modify the `list_clt_num_year` and `list_clt_dict_quarter` variables in the "Global Parameter Definitions" section.
4.  Click Run.
    -   The script will generate the 30-m daily soil temperature dataset for the selected year.
    -   To minimize the number of export tasks, images are packaged by quarter. This results in **4 tasks** appearing in the "Tasks" tab.
5.  Click Run in the "Tasks" tab to start the export to your Earth Engine Assets.
6.  Visualize results:
    -   Once the export is complete, import the generated soil temperature images from Earth Engine Assets and display them on the Map.
    -   Use the "Inspector" tab to interactively query pixel values and verify the data.

### Full Production Workflow

Use this method to reproduce the complete dataset (2001–2024).

#### Step 1: Export Tiled Data from GEE

1.  Copy the content of `02_st_30m_daily_tile_export.js` into the GEE Code Editor.
2.  Create a folder in your Google Drive to store the outputs (e.g., `ST_tiles_central_TP_2001_2024`).
3.  Click Run.
    -   The script will generate 30-m daily soil temperature dataset for the complete time series (2001–2024).
    -   Images are packaged by quarter to reduce the number of tasks.
    -   To prevent the exported GeoTIFF files from exceeding the 4GB limit (which may trigger automatic and irregular file splitting in Earth Engine exports), the quarterly images are further spatially subdivided into 4 tiles.
    -   This results in a total of **384 tasks** (24 years × 4 quarters × 4 tiles) in the "Tasks" tab.
4.  Click Run in the "Tasks" tab to export the GeoTIFF files to your Google Drive.

#### Step 2: Mosaic Tiles using R

1.  Download all tiled GeoTIFF files from Google Drive to a local folder (e.g., `./ST_tiles_central_TP_2001_2024/`).
2.  Open `03_st_30m_daily_mosaic.R` in RStudio or your preferred R environment.
3.  Ensure the `terra` package is installed (`install.packages("terra")`).
4.  Update the directory paths in the script (`str_dir_input` and `str_dir_output`) to match your local environment.
5.  Run the script. It will automatically verify that all four tiles exist for each year–quarter, mosaic the tiles, and save the final quarterly datasets with lossless compression (LZW).

---

## Output Dataset Specifications

The full production workflow generates a total of **96 multiband GeoTIFF files** (24 years × 4 quarters), with a total dataset volume of approximately **119 GB** (averaging 1.2 GB per file).

### File Structure

-   **File Naming**: `ST_<Year>_<Quarter>.tif` (e.g., ST_2020_Q1.tif).
-   **Band Naming**: `ST_YYYYMMDD` (e.g., ST_20200101 represents the soil temperature on January 1, 2020).
-   **Band Structure**: Each band represents the daily mean soil temperature for one day in that quarter.

### Data Properties

-   **Data Type**: Signed 16-bit integers (`Int16`)
-   **NoData Value**: -32768
-   **Unit**: 0.01 Degrees Celsius (°C)

### Value Conversion

The pixel values represent the soil temperature multiplied by 100. To obtain the actual soil temperature in degrees Celsius, use the following conversion:

$$\text{Soil Temperature (°C)} = \text{Pixel Value} \times 0.01$$

*(e.g., a pixel value of 1500 corresponds to 15.00°C)*.

---

## Data and Code Availability

To ensure long-term accessibility and reproducibility, both the source code and the generated dataset have been archived in recognized repositories:

-   **Dataset**: The generated 30-m daily soil temperature dataset (GeoTIFFs) is available on **ScienceDB** (DOI: `<ADD_DOI_HERE>`).
-   **Source Code**: The specific version of the code used in the manuscript is archived on **Zenodo** (DOI: `<ADD_DOI_HERE>`).

---

## Citation

If you use this code or dataset, please cite the following paper:

> Zhang, H., Ma, P., et al. (2024). A 30-m daily soil temperature dataset for alpine grasslands in central Tibetan Plateau during 2001–2024. *Scientific Data* (Submitted).

---

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

---

## Contact

If you have any questions regarding the code or the dataset, please feel free to contact:

**Haoze Zhang** <haoze_z@itpcas.ac.cn>
**Pengfei Ma** (Corresponding Author) <mapengfei@itpcas.ac.cn>
