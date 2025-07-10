# Basic-codes-in-geospatial-analysis

This Python script is designed to process Landsat 8 (or 9) Level-1 data to derive radiometrically corrected Top-of-Atmosphere (TOA) Radiance and a set of geospatial indices and parameters including NDVI, NDBI, NDWI, NBR, Brightness Temperature (BT), Land Surface Emissivity, and Land Surface Temperature (LST) in degrees Celsius. The user is required to provide the path to a folder containing the Landsat scene, which must include the original band files (e.g., B4, B5, B10) and the associated metadata text file (*_MTL.txt).

The script begins by automatically identifying and parsing the MTL file to extract essential constants such as the radiance multiplicative and additive rescaling factors for each band, and the thermal constants (K1 and K2) for Band 10. Using these parameters, TOA Radiance is calculated for all available bands using the equation: Radiance = (ML × DN) + AL, where ML and AL are the multiplicative and additive constants, and DN is the original pixel value.

Once the TOA Radiance correction is completed, the script computes several spectral indices: NDVI (using Band 5 and Band 4), NDBI (Band 6 and Band 5), NBR (Band 5 and Band 7), and NDWI (Band 3 and Band 5, if Band 3 is available). These indices are useful for analyzing vegetation health, built-up areas, burn severity, and water bodies, respectively.

The script further calculates Brightness Temperature (BT) from the corrected Band 10 radiance using the Planck function. Based on NDVI, land surface emissivity is estimated using a quadratic function which assigns values depending on NDVI thresholds. Using both BT and emissivity, Land Surface Temperature (LST) is calculated in Celsius, accounting for the wavelength of Band 10 and Planck’s radiation constants.

All outputs are saved as GeoTIFF files in the same folder, with names like TOA_Radiance_Bx.TIF, NDVI.TIF, BT.TIF, and LST.TIF. The script handles NoData values by masking pixels where the original value is zero and assigns np.nan during computation. It also prints out the minimum and maximum values of each derived product to help with quality control and validation.

The script requires Python 3.8 or higher and depends on the rasterio and numpy libraries, which must be installed in your environment. It is especially useful for researchers and analysts working in thermal remote sensing, environmental monitoring, or urban and vegetation studies, and provides a consistent and automated approach to Landsat-based LST estimation.
