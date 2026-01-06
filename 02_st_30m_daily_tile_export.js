/**
 * Script: 02_st_30m_daily_tile_export.js
 * Description: Generates the full daily 30-m soil temperature dataset (2001–2024)
 * for alpine grasslands in the Nyainqentanglha Mountains, central Tibetan Plateau.
 * The output is spatially tiled and exported to Google Drive to support the 
 * production of the complete dataset.
 * Author: Haoze Zhang
 * License: MIT License
 */


/**
 * =============================================================================
 * Section 1: Global Parameter Definitions
 * =============================================================================
 */

// Define the Albers Equal Area Conic Coordinate Reference System (CRS)
var str_clt_crs_albers =
  'PROJCS["Albers Equal Area",GEOGCS["WGS 84",DATUM["WGS_1984",' +
  'SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],' +
  'UNIT["degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],' +
  'PARAMETER["standard_parallel_1",25],' +
  'PARAMETER["standard_parallel_2",47],' +
  'PARAMETER["latitude_of_origin",0],' +
  'PARAMETER["central_meridian",105],' +
  'PARAMETER["false_easting",0],' +
  'PARAMETER["false_northing",0],UNIT["metre",1]]';

// Define quarterly processing periods (Q1–Q4)
var list_clt_dict_quarter = [
  {'index': 'Q1', 'year_offset': 0, 'begin': '01-01', 'end': '04-01'},
  {'index': 'Q2', 'year_offset': 0, 'begin': '04-01', 'end': '07-01'},
  {'index': 'Q3', 'year_offset': 0, 'begin': '07-01', 'end': '10-01'},
  {'index': 'Q4', 'year_offset': 1, 'begin': '10-01', 'end': '01-01'}
];

// Define the temporal range for batch processing (2001–2024)
var list_clt_num_year = Array.from({length: (2024 - 2001 + 1)}, function(v, k) {return k + 2001});


/**
 * =============================================================================
 * Section 2: Core Processing Functions
 * =============================================================================
 */

/**
 * Generates a projected grid and subdivides it into tiles.
 * * Context: To comply with the file size limits of GEE GeoTIFF exports, the 150 x 150 km
 * study area covering the Nyainqentanglha Mountains is divided into four 75 km x 75 km tiles.
 * These tiles are processed and exported separately to be merged later.
 *
 * @param {ee.Geometry.Point} gmt_point_center - The center point of the study area.
 * @param {String} str_clt_crs - The WKT string of the target projection (e.g., Albers).
 * @param {Number} num_clt_length_grid_side - The side length of the total study area in meters.
 * @param {Number} num_clt_count_tiles_side - The number of tiles along one side (e.g., 2 for a 2 x 2 grid).
 * @return {ee.FeatureCollection} A flattened FeatureCollection containing the defined tiles with unique IDs.
 */
function create_grid_projected(gmt_point_center, str_clt_crs, num_clt_length_grid_side, num_clt_count_tiles_side) {
      gmt_point_center = ee.Geometry(gmt_point_center);
  // Generate the study area grid based on the center point and projection
  var num_distance_buffer = ee.Number(num_clt_length_grid_side).divide(2);
  var gmt_point_center_projected = gmt_point_center.transform(str_clt_crs, 1);
  var gmt_grid_projected = gmt_point_center_projected
    .buffer(num_distance_buffer, ee.ErrorMargin(1, 'projected'), str_clt_crs)
    .bounds(ee.ErrorMargin(1, 'projected'), str_clt_crs);
  // Extract grid vertex coordinates in the target projection
  var list_coords_grid_projected = ee.List(gmt_grid_projected.coordinates().get(0));
  var num_grid_east_min = ee.Number(ee.List(list_coords_grid_projected.get(0)).get(0));
  var num_grid_east_max = ee.Number(ee.List(list_coords_grid_projected.get(2)).get(0));
  var num_grid_north_min = ee.Number(ee.List(list_coords_grid_projected.get(0)).get(1));
  var num_grid_north_max = ee.Number(ee.List(list_coords_grid_projected.get(2)).get(1));
  // Calculate the width and height of each individual tile
  var num_tile_width = num_grid_east_max.subtract(num_grid_east_min).divide(num_clt_count_tiles_side);
  var num_tile_height = num_grid_north_max.subtract(num_grid_north_min).divide(num_clt_count_tiles_side);
  // Generate the tile grid
  var list_num_i = ee.List.sequence(0, ee.Number(num_clt_count_tiles_side).subtract(1)); // Row index (Northing step)
  var list_num_j = ee.List.sequence(0, ee.Number(num_clt_count_tiles_side).subtract(1)); // Column index (Easting step)
  var list_fc_tiles = list_num_i.map(function(num_i) {
        num_i = ee.Number(num_i);
    var list_ftr_tile_i_j = list_num_j.map(function(num_j) {
          num_j = ee.Number(num_j);
      // Create a unique ID for each tile (e.g., tile_0_0, tile_0_1)
      var str_tile_id = ee.String('tile_').cat(ee.String(num_i.toInt())).cat('_').cat(ee.String(num_j.toInt()));
      // Calculate tile boundaries
      var num_tile_north_min = num_grid_north_min.add(num_tile_height.multiply(num_i));
      var num_tile_north_max = num_tile_north_min.add(num_tile_height);
      var num_tile_east_min = num_grid_east_min.add(num_tile_width.multiply(num_j));
      var num_tile_east_max = num_tile_east_min.add(num_tile_width);
      // Construct the tile geometry and feature
      var gmt_tile_i_j = ee.Geometry.Rectangle([num_tile_east_min, num_tile_north_min, num_tile_east_max, num_tile_north_max], str_clt_crs, false);
      var ftr_tile_i_j = ee.Feature(gmt_tile_i_j, {'tile_id': str_tile_id});
      return ftr_tile_i_j;
    });
    var fc_tiles = ee.FeatureCollection(list_ftr_tile_i_j);
    return fc_tiles;
  });
  var fc_tiles = ee.FeatureCollection(list_fc_tiles).flatten();
  return fc_tiles;
}

/**
 * Preprocesses the MODIS LST product.
 * * Operations include Digital Number (DN) conversion, pixel quality control, and
 * pixel-wise temporal linear interpolation to fill data gaps caused mainly by cloud cover.
 *
 * @param {Number} num_clt_year - The target year for processing.
 * @param {Number} num_clt_window_interp - The temporal window size (in days) for linear interpolation.
 * @param {ee.Image} img_mask_pixel - A binary mask image defining the valid tile area pixels.
 * @return {ee.ImageCollection} A collection of daily, gap-filled MODIS LST images.
 */
function tidy_LST_MODIS(num_clt_year, num_clt_window_interp, img_mask_pixel) {
      img_mask_pixel = ee.Image(img_mask_pixel);
  // Setup parameter: MODIS LST preprocessing
  var num_year = ee.Number(num_clt_year).toInt();
  var num_window_interp = ee.Number(num_clt_window_interp);
  var str_date_begin = ee.Date(ee.String(num_year).cat('-01-01')).format('YYYY-MM-dd');  // inclusive
  var date_begin = ee.Date(str_date_begin);
  var date_begin_extend = date_begin.advance(-2, 'month');  // inclusive
  var str_date_end = ee.Date(str_date_begin).advance(1, 'years').format('YYYY-MM-dd');  // exclusive
  var date_end = ee.Date(str_date_end);
  var date_end_extend = date_end.advance(2, 'month');  // exclusive
  // Preprocess MODIS LST: QC filtering, temporal extension (+/- 2 months) for interpolation, and timestamp band addition
  var ic_lst_time_m_extend = ee.ImageCollection('MODIS/061/MOD11A1').filterDate(date_begin_extend, date_end_extend)
    .map(function(img) {
          img = ee.Image(img);
      var img_qc_day = img.select('QC_Day');
      var img_mask_qc_day = img_qc_day.bitwiseAnd(3).lte(1)  // Bits 0-1: Exclude 2 (Cloud) and 3 (Other errors)
        .and(img_qc_day.rightShift(4).bitwiseAnd(3).lte(1))  // Bits 4-5: Exclude 2 (Emissivity error <= 0.04) and 3 (Emissivity error > 0.04)
        .and(img_qc_day.rightShift(6).bitwiseAnd(3).lte(1));  // Bits 6-7: Exclude 2 (LST error <= 3K) and 3 (LST error > 3K)
      var img_lst = img.select('LST_Day_1km').updateMask(img_mask_qc_day).updateMask(img_mask_pixel)
        .multiply(0.02).toFloat().rename('M_LST');
      var img_time = ee.Image.constant(img.date().millis()).toFloat().updateMask(img_lst.mask()).rename('M_time');
      var img_return = img_lst.addBands(img_time).copyProperties(img, ['system:time_start']);
      return img_return;
    }).sort('system:time_start');
  // Setup parameter: Pixel-wise temporal linear interpolation
  var num_count_days = date_end.difference(date_begin, 'days').toInt();
  var list_num_index_days_year = ee.List.sequence(0, num_count_days.subtract(1));
  var list_img_empty = list_num_index_days_year.map(function(num_index) {
        num_index = ee.Number(num_index).toInt();
    var num_millis_index = date_begin.advance(num_index, 'days').millis();
    var str_date_format = ee.Date(num_millis_index).format('YYYY-MM-dd');
    var img_return = ee.Image().set({'system:time_start': num_millis_index, 'DATE_ACQUIRED': str_date_format});
    return img_return;
  });
  var ic_empty = ee.ImageCollection.fromImages(list_img_empty).sort('system:time_start');
  var num_millis_window_interp = num_window_interp.multiply(ee.Number(24)).multiply(ee.Number(60))
    .multiply(ee.Number(60)).multiply(ee.Number(1000));
  var filter_lst_m_frame = ee.Filter.maxDifference({
    'difference': num_millis_window_interp, 'leftField': 'system:time_start', 'rightField': 'system:time_start'
  });
  var filter_lst_m_before = ee.Filter.greaterThanOrEquals({'leftField': 'system:time_start', 'rightField': 'system:time_start'});
  var filter_lst_m_after = ee.Filter.lessThanOrEquals({'leftField': 'system:time_start', 'rightField': 'system:time_start'});
  var filter_lst_m_frame_before = ee.Filter.and(filter_lst_m_frame, filter_lst_m_before);
  var filter_lst_m_frame_after = ee.Filter.and(filter_lst_m_frame, filter_lst_m_after);
  var join_lst_m_before = ee.Join.saveAll({
    'matchesKey': 'ic_lst_m_before', 'ordering': 'system:time_start', 'ascending': true
  });
  var join_lst_m_after = ee.Join.saveAll({
    'matchesKey': 'ic_lst_m_after', 'ordering': 'system:time_start', 'ascending': false
  });
  var ic_empty_before = join_lst_m_before.apply({
    'primary': ic_empty, 'secondary': ic_lst_time_m_extend, 'condition': filter_lst_m_frame_before
  });
  var ic_empty_before_after = join_lst_m_after.apply({
    'primary': ic_empty_before, 'secondary': ic_lst_time_m_extend, 'condition': filter_lst_m_frame_after
  });
  // Execute pixel-wise temporal linear interpolation
  var ic_lst_m_interpolated = ic_empty_before_after.map(function(img_empty) {
        img_empty = ee.Image(img_empty);
    // Setup parameter
    var num_millis_index = img_empty.date().millis();
    var img_lst_m_mosaic_before = ee.ImageCollection.fromImages(ee.List(img_empty.get('ic_lst_m_before'))).mosaic();
    var img_lst_m_mosaic_after = ee.ImageCollection.fromImages(ee.List(img_empty.get('ic_lst_m_after'))).mosaic();
    var img_time_x1 = img_lst_m_mosaic_before.select('M_time');
    var img_lst_y1 = img_lst_m_mosaic_before.select('M_LST');
    var img_time_x2 = img_lst_m_mosaic_after.select('M_time');
    var img_lst_y2 = img_lst_m_mosaic_after.select('M_LST');
    var img_time_x0 = ee.Image.constant(num_millis_index).toFloat();
    // Apply linear interpolation
    var img_ratio = img_time_x0.subtract(img_time_x1).divide(img_time_x2.subtract(img_time_x1))
      .where(img_time_x2.eq(img_time_x1), 0);
    var img_lst_y0 = img_lst_y2.subtract(img_lst_y1).multiply(img_ratio).add(img_lst_y1).toFloat()
      .rename('M_LST_interp').copyProperties(img_empty, ['system:time_start', 'DATE_ACQUIRED']);
    return img_lst_y0;
  });
  return ee.ImageCollection(ic_lst_m_interpolated);
}

/**
 * Derives a normalized temporal shape factor from the interpolated MODIS LST time series.
 * * Calculated by normalizing the daily LST against the annual maximum LST pixel-wise.
 *
 * @param {ee.ImageCollection} ic_lst_m - The interpolated MODIS LST ImageCollection.
 * @return {ee.ImageCollection} A collection of images representing the daily MODIS LST ratio (0-1).
 */
function dynamics_LST_MODIS(ic_lst_m) {
      ic_lst_m = ee.ImageCollection(ic_lst_m);
  // Calculate the annual maximum LST
  var img_lst_m_max = ic_lst_m.max().rename('M_LST_interp_max');
  // Normalize daily LST by the annual maximum
  var ic_lst_m_ratio = ic_lst_m.map(function(img_lst_m) {
        img_lst_m = ee.Image(img_lst_m);
    var img_lst_m_ratio = img_lst_m.divide(img_lst_m_max).rename('M_LST_ratio')
      .copyProperties(img_lst_m, ['system:time_start', 'DATE_ACQUIRED']);
    return img_lst_m_ratio;
  });
  return ic_lst_m_ratio.sort('system:time_start');
}

/**
 * Preprocesses Landsat 4/5/7 thermal bands.
 * * Operations include DN conversion and pixel quality control.
 *
 * @param {ee.Image} img - A raw Landsat 4, 5, or 7 image.
 * @return {ee.Image} A preprocessed LST image in Kelvin.
 */
function tidy_band_ST_Landsat457(img) {
      img = ee.Image(img);
  // Convert Landsat surface temperature band from DN to Kelvin
  var img_st = img.select('ST_B.+').multiply(0.00341802).add(149.0).toFloat();
  // Define Quality Assessment (QA) mask
  var img_qa = img.select('QA_PIXEL');
  var img_mask_qa = img_qa.bitwiseAnd(1).eq(0)  // Bit 0: Fill
    .and(img_qa.bitwiseAnd(2).eq(0))  // Bit 1: Dilated Cloud
    .and(img_qa.bitwiseAnd(8).eq(0))  // Bit 3: Cloud
    .and(img_qa.bitwiseAnd(16).eq(0));  // Bit 4: Cloud Shadow
  // Define saturation mask
  var img_mask_saturation = img.select('QA_RADSAT').eq(0);
  // Apply QA and saturation masks to filter out low-quality pixels
  var img_st_tidy = img_st.updateMask(img_mask_qa).updateMask(img_mask_saturation).rename('L_LST')
    .copyProperties(img, ['system:time_start', 'DATE_ACQUIRED']);
  return img_st_tidy;
}

/**
 * Preprocesses Landsat 8/9 thermal bands with inter-sensor harmonization.
 * * Operations include DN conversion, inter-sensor harmonization, and pixel quality control.
 *
 * @param {ee.Image} img - A raw Landsat 8 or 9 image.
 * @return {ee.Image} A preprocessed and harmonized LST image in Kelvin.
 */
function tidy_band_ST_Landsat89(img) {
      img = ee.Image(img);
  // Convert Landsat surface temperature band from DN to Kelvin and apply inter-sensor harmonization
  var img_st_corrected = img.select('ST_B.+').multiply(0.00341802).add(149.0)
    .multiply(0.936148532204192).add(18.2329664975714).toFloat();
  // Define Quality Assessment (QA) mask
  var img_qa = img.select('QA_PIXEL');
  var img_mask_qa = img_qa.bitwiseAnd(1).eq(0)  // Bit 0: Fill
    .and(img_qa.bitwiseAnd(2).eq(0))  // Bit 1: Dilated Cloud
    .and(img_qa.bitwiseAnd(4).eq(0))  // Bit 2: Cirrus
    .and(img_qa.bitwiseAnd(8).eq(0))  // Bit 3: Cloud
    .and(img_qa.bitwiseAnd(16).eq(0));  // Bit 4: Cloud Shadow
  // Define saturation mask
  var img_mask_saturation = img.select('QA_RADSAT').eq(0);
  // Apply QA and saturation masks to filter out low-quality pixels
  var img_st_tidy = img_st_corrected.updateMask(img_mask_qa).updateMask(img_mask_saturation).rename('L_LST')
    .copyProperties(img, ['system:time_start', 'DATE_ACQUIRED']);
  return img_st_tidy;
}

/**
 * Integrates Landsat 4-9 collections and mosaics same-day images.
 * * Filters collections by bounds and date, applies sensor-specific preprocessing,
 * merges all Landsat imagery, and mosaics images acquired on the same date.
 *
 * @param {ee.Geometry} gmt_select - The geometry to filter the collections.
 * @param {Number|String} num_clt_year - The target year.
 * @param {ee.Image} img_mask_pixel - The valid area mask.
 * @param {String} str_clt_crs - The WKT string of the target projection (e.g., Albers).
 * @return {ee.ImageCollection} A collection of preprocessed, daily mosaicked Landsat LST images.
 */
function tidy_LST_Landsat_unique(gmt_select, num_clt_year, img_mask_pixel, str_clt_crs) {
      gmt_select = ee.Geometry(gmt_select);
      img_mask_pixel = ee.Image(img_mask_pixel);
  // Setup parameter
  var num_year = ee.Number(num_clt_year).toInt();
  var str_date_begin = ee.Date(ee.String(num_year).cat('-01-01')).format('YYYY-MM-dd');  // inclusive
  var str_date_end = ee.Date(str_date_begin).advance(1, 'years').format('YYYY-MM-dd');  // exclusive
  // Load Landsat Collections
  var ic_l4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2');
  var ic_l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
  // Filter L7 collection: Exclude data acquired after the orbit change (lowered orbit) on May 4, 2022
  var ic_l7_filtered = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterDate(ee.String('1999-05-27'), ee.String('2022-05-04'));
  var ic_l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
  var ic_l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');
  // Integrate Landsat 4/5/7 collections
  var ic_l457 = ic_l4.merge(ic_l5).merge(ic_l7_filtered)
    .filterBounds(gmt_select).filterDate(str_date_begin, str_date_end)
    .map(function(img) {return ee.Image(img).updateMask(img_mask_pixel);})
    .map(tidy_band_ST_Landsat457);
  // Integrate Landsat 8/9 collections
  var ic_l89 = ic_l8.merge(ic_l9)
    .filterBounds(gmt_select).filterDate(str_date_begin, str_date_end)
    .map(function(img) {return ee.Image(img).updateMask(img_mask_pixel);})
    .map(tidy_band_ST_Landsat89);
  // Integrate all Landsat collections
  var ic_lst_l_tidy = ic_l457.merge(ic_l89).sort('system:time_start');
  // Mosaic images acquired on the same date
  var ic_empty = ic_lst_l_tidy.map(function(img_lst_l_tidy) {
        img_lst_l_tidy = ee.Image(img_lst_l_tidy);
    var img_empty = ee.Image().copyProperties(img_lst_l_tidy, ['DATE_ACQUIRED']);
    return img_empty;
  }).distinct('DATE_ACQUIRED');
  var filter_same_day = ee.Filter.equals({'leftField': 'DATE_ACQUIRED', 'rightField': 'DATE_ACQUIRED'});
  var join_lst_l_same_day = ee.Join.saveAll({
    'matchesKey': 'ic_lst_l_same_day', 'ordering': 'system:time_start', 'ascending': true
  });
  var ic_lst_l_joined = join_lst_l_same_day.apply({
    'primary': ic_empty, 'secondary': ic_lst_l_tidy, 'condition': filter_same_day
  });
  var ic_lst_l_distinct = ic_lst_l_joined.map(function(img_joined) {
        img_joined = ee.Image(img_joined);
    var str_date = img_joined.get('DATE_ACQUIRED');
    var num_millis_date = ee.Date(str_date).millis();
    var img_lst_l_mosaic = ee.ImageCollection.fromImages(ee.List(img_joined.get('ic_lst_l_same_day'))).mosaic()
      .setDefaultProjection({'crs': str_clt_crs, 'scale': 30})
      .set({'system:time_start': num_millis_date, 'DATE_ACQUIRED': str_date});
    return img_lst_l_mosaic;
  });
  return ee.ImageCollection(ic_lst_l_distinct).sort('system:time_start');
}

/**
 * Estimates the theoretical Landsat annual maximum LST.
 * * Uses the sparse Landsat observations and the corresponding MODIS temporal shape
 * factor to calculate candidate the Landsat annual maximum values. Outliers are
 * removed using the Interquartile Range (IQR) method.
 *
 * @param {ee.ImageCollection} ic_lst_l - The tidy Landsat LST collection.
 * @param {ee.ImageCollection} ic_lst_m_ratio - The MODIS-derived temporal shape factor collection.
 * @param {String} str_clt_crs - The WKT string of the target projection (e.g., Albers).
 * @return {ee.Image} An image representing the theoretical Landsat annual maximum LST.
 */
function max_LST_Landsat_theoretical(ic_lst_l, ic_lst_m_ratio, str_clt_crs) {
      ic_lst_l = ee.ImageCollection(ic_lst_l);
      ic_lst_m_ratio = ee.ImageCollection(ic_lst_m_ratio);
  // Join Landsat and MODIS temporal shape factor collections based on the 'DATE_ACQUIRED' property
  var join_date_matched = ee.Join.inner({
    'primaryKey': 'img_lst_l_matched', 'secondaryKey': 'img_lst_m_ratio_matched'
  });
  var fc_lst_date_matched = join_date_matched.apply({
    'primary': ic_lst_l, 'secondary': ic_lst_m_ratio,
    'condition': ee.Filter.equals({'leftField': 'DATE_ACQUIRED', 'rightField': 'DATE_ACQUIRED'})
  });
  // Calculate candidate annual maximum values
  var ic_lst_l_maxs_theo = ee.ImageCollection(fc_lst_date_matched.map(function(ftr_lst_date_matched) {
        ftr_lst_date_matched = ee.Feature(ftr_lst_date_matched);
    var img_lst_l_matched = ee.Image(ftr_lst_date_matched.get('img_lst_l_matched'));
    var img_lst_m_ratio_matched = ee.Image(ftr_lst_date_matched.get('img_lst_m_ratio_matched')).resample('bicubic');
    var img_lst_l_maxs_theo = img_lst_l_matched.divide(img_lst_m_ratio_matched).toFloat().rename('L_LST_maxs_theo')
      .copyProperties(img_lst_l_matched, ['system:time_start', 'DATE_ACQUIRED']);
    return img_lst_l_maxs_theo;
  }));
  // Calculate pixel-wise Interquartile Range (IQR) of candidate values of the Landsat annual maximum
  var img_lst_l_maxs_theo_p25 = ic_lst_l_maxs_theo.reduce(ee.Reducer.percentile([25]));
  var img_lst_l_maxs_theo_p75 = ic_lst_l_maxs_theo.reduce(ee.Reducer.percentile([75]));
  var img_lst_l_maxs_theo_iqr = img_lst_l_maxs_theo_p75.subtract(img_lst_l_maxs_theo_p25);
  // Generate masks to exclude outliers based on the IQR
  var img_lst_l_maxs_theo_floor = img_lst_l_maxs_theo_p25.subtract(img_lst_l_maxs_theo_iqr.multiply(ee.Number(1.5)));
  var img_lst_l_maxs_theo_ceiling = img_lst_l_maxs_theo_p75.add(img_lst_l_maxs_theo_iqr.multiply(ee.Number(1.5)));
  // Apply outlier masks
  var ic_lst_l_maxs_theo_valid = ic_lst_l_maxs_theo.map(function(img_lst_l_maxs_theo) {
        img_lst_l_maxs_theo = ee.Image(img_lst_l_maxs_theo);
    var img_mask_floor = img_lst_l_maxs_theo.gt(img_lst_l_maxs_theo_floor);
    var img_mask_ceil = img_lst_l_maxs_theo.lt(img_lst_l_maxs_theo_ceiling);
    var img_lst_l_maxs_theo_valid = img_lst_l_maxs_theo.updateMask(img_mask_floor).updateMask(img_mask_ceil).toFloat()
      .rename('L_LST_maxs_theo_valid').copyProperties(img_lst_l_maxs_theo, ['system:time_start', 'DATE_ACQUIRED']);
    return img_lst_l_maxs_theo_valid;
  });
  // Compute the mean of valid candidate values to obtain the final theoretical annual maximum LST
  var img_lst_l_max_theo = ic_lst_l_maxs_theo_valid.mean().toFloat()
    .setDefaultProjection({'crs': str_clt_crs, 'scale': 30}).rename('L_LST_max_theo');
  return img_lst_l_max_theo;
}

/**
 * Reconstructs the daily 30-m fused LST.
 * * Combines the high-temporal-resolution MODIS shape factor with the high-spatial-resolution
 * Landsat theoretical maximum.
 *
 * @param {ee.ImageCollection} ic_lst_m_ratio - The MODIS-derived temporal shape factor collection.
 * @param {ee.Image} img_lst_l_max_theo - Landsat theoretical annual maximum image.
 * @return {ee.ImageCollection} A collection of daily fused 30-m LST images.
 */
function reconstruct_LST_fused(ic_lst_m_ratio, img_lst_l_max_theo) {
      ic_lst_m_ratio = ee.ImageCollection(ic_lst_m_ratio);
      img_lst_l_max_theo = ee.Image(img_lst_l_max_theo);
  var ic_lst_fused = ic_lst_m_ratio.map(function(img_lst_m_ratio) {
        img_lst_m_ratio = ee.Image(img_lst_m_ratio);
    // Reconstruct LST: Multiply the MODIS temporal shape factor by the theoretical Landsat annual maximum
    var img_lst_fused = img_lst_m_ratio.resample('bicubic').multiply(img_lst_l_max_theo).rename('LST_fused')
      .copyProperties(img_lst_m_ratio, ['system:time_start', 'DATE_ACQUIRED']);
    return img_lst_fused;
  });
  return ic_lst_fused;
}

/**
 * Calibrates fused LST to daily mean soil temperature at -5 cm
 * * Applies a monthly linear regression model derived from in situ measurements.
 *
 * @param {ee.ImageCollection} ic_lst_fused - The fused LST collection.
 * @return {ee.ImageCollection} A collection of calibrated soil temperature images.
 */
function calibrate_LST_fused(ic_lst_fused) {
      ic_lst_fused = ee.ImageCollection(ic_lst_fused);
  // Define and construct monthly calibration coefficient images (scale and offset)
  // Monthly coefficients were derived from the linear regression model using in situ soil temperature
  // measurements from the 10 weather stations along elevations of 4300–5500 m in the Nyainqentanglha
  // Mountains during 2006–2010 and the corresponding 30-m daily fused LST data
  var list_clt_dict_coef_calib = [
    {'month': 1, 'scale': 0.3372, 'offset': 173.3712},
    {'month': 2, 'scale': 0.3293, 'offset': 175.5622},
    {'month': 3, 'scale': 0.2828, 'offset': 190.0353},
    {'month': 4, 'scale': 0.2981, 'offset': 187.6972},
    {'month': 5, 'scale': 0.2391, 'offset': 208.3056},
    {'month': 6, 'scale': 0.1865, 'offset': 227.121},
    {'month': 7, 'scale': 0.1638, 'offset': 235.3277},
    {'month': 8, 'scale': 0.1829, 'offset': 229.0344},
    {'month': 9, 'scale': 0.239, 'offset': 210.4638},
    {'month': 10, 'scale': 0.2893, 'offset': 192.7716},
    {'month': 11, 'scale': 0.2747, 'offset': 194.3303},
    {'month': 12, 'scale': 0.3982, 'offset': 156.405}
  ];
  var list_img_coef_calib = list_clt_dict_coef_calib.map(function(dict_clt_coef_calib) {
    var img_coef_calib = ee.Image.constant([dict_clt_coef_calib.scale, dict_clt_coef_calib.offset])
      .rename(['scale', 'offset']).set({'month': dict_clt_coef_calib.month});
    return img_coef_calib;
  });
  var ic_coef_calib = ee.ImageCollection.fromImages(list_img_coef_calib);
  // Apply monthly calibration to convert fused LST to Soil Temperature (ST)
  var ic_st = ic_lst_fused.map(function(img_lst_fused) {
        img_lst_fused = ee.Image(img_lst_fused);
    // Retrieve calibration coefficients matching the image month
    var num_month = img_lst_fused.date().get('month');
    var img_coef_calib = ee.Image(ic_coef_calib.filter(ee.Filter.eq('month', num_month)).first());
    var img_scale = img_coef_calib.select('scale');
    var img_offset = img_coef_calib.select('offset');
    // Apply calibration: ST = scale * LST_fused + offset
    var img_st = img_lst_fused.multiply(img_scale).add(img_offset).rename('ST')
        .copyProperties(img_lst_fused, ['system:time_start', 'DATE_ACQUIRED']);
    return img_st;
  });
  return ic_st.sort('system:time_start');
}

/**
 * Main Wrapper Function: Generate and export Soil Temperature (ST) dataset
 * * Orchestrates all preprocessing, fusion, and calibration steps, then compresses
 * and exports the results as GeoTIFFs to Google Drive.
 *
 * @param {Number|String} num_clt_year - The target year for processing.
 * @param {ee.Geometry} gmt_grid - The specific tile geometry to export.
 * @param {Object} dict_clt_quarter - The quarter configuration object (index, start/end dates).
 * @param {String} str_clt_tile_id - The unique identifier for the tile.
 * @param {String} str_clt_crs - The WKT string of the target projection (e.g., Albers).
 */
function export_ST(num_clt_year, gmt_grid, dict_clt_quarter, str_clt_tile_id, str_clt_crs) {
      gmt_grid = ee.Geometry(gmt_grid);
  // Setup parameter
  var date_quarter_begin = ee.Date(ee.String(ee.Number(num_clt_year)).cat('-').cat(ee.String(dict_clt_quarter.begin)));  // inclusive
  var date_quarter_end = ee.Date(ee.String(ee.Number(num_clt_year).add(dict_clt_quarter.year_offset)).cat('-').cat(ee.String(dict_clt_quarter.end)));  // exclusive
  var img_mask_pixel = ee.Image.constant(1).clip(gmt_grid.buffer({'distance': 10000})).selfMask();
  
  // Preprocess MODIS LST product (pixel quality control, temporal linear interpolation)
  var ic_lst_m_interp = tidy_LST_MODIS(
    /* num_clt_year = */ num_clt_year,
    /* num_clt_window_interp = */ 30,
    /* img_mask_pixel = */ img_mask_pixel
  );
  
  // Derive normalized temporal shape factor from MODIS LST
  var ic_lst_m_ratio = dynamics_LST_MODIS(/* ic_lst_m = */ ic_lst_m_interp);
  
  // Preprocess Landsat LST product (pixel quality control, inter-sensor harmonization, mosaic same-day images)
  var ic_lst_l_tidy = tidy_LST_Landsat_unique(
    /* gmt_select = */ gmt_grid,
    /* num_clt_year = */ num_clt_year,
    /* img_mask_pixel = */ img_mask_pixel,
    /* str_clt_crs = */ str_clt_crs
  );
  
  // Estimate the theoretical Landsat annual maximum
  var img_lst_l_max_theo = max_LST_Landsat_theoretical(
    /* ic_lst_l = */ ic_lst_l_tidy,
    /* ic_lst_m_ratio = */ ic_lst_m_ratio,
    /* str_clt_crs = */ str_clt_crs
  );
  
  // Reconstruct daily 30-m fused LST
  var ic_lst_fused = reconstruct_LST_fused(
    /* ic_lst_m_ratio = */ ic_lst_m_ratio,
    /* img_lst_l_max_theo = */ img_lst_l_max_theo
  );
  
  // Calibrate the fused LST to daily mean soil temperature at -5 cm
  var ic_st = calibrate_LST_fused(/* ic_lst_fused = */ ic_lst_fused);
  
  // Compress and package ST images by quarter
  var ic_st_compressed = ic_st.filterDate(date_quarter_begin, date_quarter_end).map(function(img) {
        img = ee.Image(img);
    var str_date_format = ee.String('ST_').cat(img.date().format('YYYYMMdd'));
    // Convert Kelvin to Celsius, scale by 100, cast to Int16 to reduce storage size, and use -32768 for NoData
    var img_return = img.subtract(273.15).multiply(100).round().toInt16().unmask(-32768).rename(str_date_format)
      .copyProperties(img, ['system:time_start', 'DATE_ACQUIRED']).set({'band_name': str_date_format});
    return img_return;
  });
  var list_str_band_name = ic_st_compressed.aggregate_array('band_name');
  var img_st_compressed = ic_st_compressed.toBands().rename(list_str_band_name)
    .set({'year': num_clt_year, 'quarter': dict_clt_quarter.index, 'tile_id': str_clt_tile_id, 'count_bands': ic_st_compressed.size()});
  
  // Export ST dataset to Google Drive
  Export.image.toDrive({
    'image': img_st_compressed,
    'description': 'img_st_' + num_clt_year + '_' + dict_clt_quarter.index + '_' + str_clt_tile_id,
    'folder': 'ST_tiles_central_TP_2001_2024',
    'region': gmt_grid,
    'scale': 30,
    'crs': str_clt_crs,
    'maxPixels': 1e9,  // increase if exporting larger regions
    'fileFormat': 'GeoTIFF',
    'formatOptions': {'cloudOptimized': true}
  });
}


/**
 * =============================================================================
 * Section 3: Execution Block
 * =============================================================================
 */

// Define the study area: A 150 × 150 km grid covering the Nyainqentanglha Mountains
// Note: Coordinates represent the centroid of the study grid. Study area is divided into 4 tiles (75 x 75 km each) for export.
var gmt_point_center = ee.Geometry.Point(90.62, 30.32);
// Create the tile structure
var fc_tiles_albers = create_grid_projected(
  /* gmt_point_center = */ gmt_point_center,
  /* str_clt_crs = */ str_clt_crs_albers,
  /* num_clt_length_grid_side = */ 150000,
  /* num_clt_count_tiles_side = */ 2
);
// Transfer the tile list to the client-side object for batch task generation
// Note: .getInfo() is required here to iterate through tiles and create separate export tasks for each.
var list_clt_ftr_tile_i_j = fc_tiles_albers.toList(fc_tiles_albers.size()).getInfo();

// Batch Execution: Iterate through years, quarters, and tiles to generate export tasks
list_clt_num_year.forEach(function(num_clt_year) {
  list_clt_dict_quarter.forEach(function(dict_clt_quarter) {
    list_clt_ftr_tile_i_j.forEach(function(ftr_clt_tile_i_j) {
      var ftr_tile_i_j = ee.Feature(ftr_clt_tile_i_j);
      var gmt_tile_i_j = ftr_tile_i_j.geometry();
      var str_clt_tile_id = ftr_clt_tile_i_j.properties.tile_id;
      print('ST ' + num_clt_year + ' ' + dict_clt_quarter.index + ' ' + str_clt_tile_id);
      export_ST(
        /* num_clt_year = */ num_clt_year,
        /* gmt_grid = */ gmt_tile_i_j,
        /* dict_clt_quarter = */ dict_clt_quarter,
        /* str_clt_tile_id = */ str_clt_tile_id,
        /* str_clt_crs = */ str_clt_crs_albers
      );
    });
  });
});
