import os
import re
import rasterio
import numpy as np

# -------- USER INPUT --------
landsat_folder = r"C:\Users\User\Desktop\Goa CVI\LC09_L1TP_147049_20250105_20250105_02_T1"

# -------- FUNCTIONS --------

def find_mtl_file(folder):
    for file in os.listdir(folder):
        if file.endswith("_MTL.txt"):
            return os.path.join(folder, file)
    raise FileNotFoundError("MTL metadata file not found in the folder.")

def parse_mtl_constants(mtl_path):
    radiance_mult = {}
    radiance_add = {}
    k1_constant = None
    k2_constant = None

    with open(mtl_path, 'r') as f:
        for line in f:
            line_upper = line.upper()

            mult_match = re.match(r'\s*RADIANCE_MULT_BAND_(\d+)\s=\s([0-9.EE+-]+)', line_upper)
            add_match  = re.match(r'\s*RADIANCE_ADD_BAND_(\d+)\s=\s([0-9.EE+-]+)', line_upper)
            k1_match   = re.match(r'\s*K1_CONSTANT_BAND_10\s=\s([0-9.EE+-]+)', line_upper)
            k2_match   = re.match(r'\s*K2_CONSTANT_BAND_10\s=\s([0-9.EE+-]+)', line_upper)

            if mult_match:
                band = int(mult_match.group(1))
                radiance_mult[band] = float(mult_match.group(2))

            if add_match:
                band = int(add_match.group(1))
                radiance_add[band] = float(add_match.group(2))

            if k1_match:
                k1_constant = float(k1_match.group(1))

            if k2_match:
                k2_constant = float(k2_match.group(1))

    return radiance_mult, radiance_add, k1_constant, k2_constant

def apply_radiance_correction(band_path, mult_factor, add_factor, output_path):
    with rasterio.open(band_path) as src:
        band = src.read(1).astype(np.float32)
        profile = src.profile

        nodata_mask = band == 0
        radiance = mult_factor * band + add_factor
        radiance[nodata_mask] = np.nan

        profile.update(dtype=rasterio.float32, nodata=np.nan)

    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(radiance, 1)

def read_band(path):
    with rasterio.open(path) as src:
        band = src.read(1).astype(np.float32)
        profile = src.profile
    return band, profile

def write_output(path, array, profile):
    profile.update(dtype=rasterio.float32, nodata=np.nan)
    
    # Remove .ovr if it exists to prevent permission error
    ovr_path = path + ".ovr"
    if os.path.exists(ovr_path):
        try:
            os.remove(ovr_path)
        except Exception as e:
            print(f"Warning: Could not delete {ovr_path} — {e}")

    with rasterio.open(path, 'w', **profile) as dst:
        dst.write(array, 1)


def compute_index(band_a, band_b):
    with np.errstate(divide='ignore', invalid='ignore'):
        index = (band_a - band_b) / (band_a + band_b + 1e-10)
        index[(band_a == 0) | (band_b == 0)] = np.nan
    return index

def compute_bt(radiance_band, k1, k2):
    bt = k2 / (np.log((k1 / (radiance_band + 1e-10)) + 1))
    bt[radiance_band == 0] = np.nan
    return bt

def compute_emissivity(ndvi):
    pv = ((ndvi - 0.2) / (0.5 - 0.2)) ** 2
    emissivity = np.where(ndvi < 0.2, 0.977,
                  np.where(ndvi > 0.5, 0.990,
                  0.977 + 0.003 * pv))
    emissivity[np.isnan(ndvi)] = np.nan
    return emissivity

def compute_lst(bt, emissivity, wavelength=10.8e-6):
    rho = 1.438e-2  # m·K
    lst = bt / (1 + ((wavelength * bt) / rho) * np.log(emissivity + 1e-10))
    lst[np.isnan(bt) | np.isnan(emissivity)] = np.nan
    return lst - 273.15  # convert K to °C

def print_min_max(name, array):
    print(f"{name}: min = {np.nanmin(array):.2f}, max = {np.nanmax(array):.2f}")

# -------- EXECUTION --------

print("Starting Landsat 8 TOA Radiance and LST Processing...\n")

# Step 1: Parse MTL
mtl_path = find_mtl_file(landsat_folder)
radiance_mult, radiance_add, k1, k2 = parse_mtl_constants(mtl_path)

print("Radiance Multiplicative Factors:")
for b in sorted(radiance_mult.keys()):
    print(f"Band {b}: {radiance_mult[b]}")

print("\nRadiance Additive Factors:")
for b in sorted(radiance_add.keys()):
    print(f"Band {b}: {radiance_add[b]}")

print(f"\nThermal Constants:\nK1 = {k1}, K2 = {k2}\n")

# Step 2: Apply TOA Radiance correction
corrected_bands = {}
for band_num in radiance_mult:
    band_in = os.path.join(landsat_folder, f"LC09_L1TP_147049_20250105_20250105_02_T1_B{band_num}.TIF")
    band_out = os.path.join(landsat_folder, f"TOA_Radiance_B{band_num}.TIF")

    if os.path.exists(band_in):
        apply_radiance_correction(band_in, radiance_mult[band_num], radiance_add[band_num], band_out)
        print(f"Corrected Band {band_num} saved.")
        corrected_bands[f"B{band_num}"] = band_out
    else:
        print(f"Band {band_num} not found. Skipping...")

# Step 3: Load necessary corrected bands
needed_bands = ['B3', 'B4', 'B5', 'B6', 'B7', 'B10']
bands = {}
for b in needed_bands:
    path = os.path.join(landsat_folder, f"TOA_Radiance_{b}.TIF")
    if os.path.exists(path):
        bands[b], profile = read_band(path)
    else:
        print(f"{b} not found.")

# Step 4: Compute indices and parameters
ndvi = compute_index(bands['B5'], bands['B4'])  # NIR - RED
ndbi = compute_index(bands['B6'], bands['B5'])  # SWIR1 - NIR
nbr  = compute_index(bands['B5'], bands['B7'])  # NIR - SWIR2
if 'B3' in bands:
    ndwi = compute_index(bands['B3'], bands['B5'])  # GREEN - NIR
else:
    ndwi = None

# Step 5: LST
bt  = compute_bt(bands['B10'], k1, k2)          # Brightness Temperature
eps = compute_emissivity(ndvi)                  # Emissivity
lst = compute_lst(bt, eps)                      # LST (°C)

# Step 6: Save and print results
write_output(os.path.join(landsat_folder, "NDVI.TIF"), ndvi, profile)
write_output(os.path.join(landsat_folder, "NDBI.TIF"), ndbi, profile)
write_output(os.path.join(landsat_folder, "NBR.TIF"), nbr, profile)
write_output(os.path.join(landsat_folder, "BT.TIF"), bt, profile)
write_output(os.path.join(landsat_folder, "EMISSIVITY.TIF"), eps, profile)
write_output(os.path.join(landsat_folder, "LST.TIF"), lst, profile)

if ndwi is not None:
    write_output(os.path.join(landsat_folder, "NDWI.TIF"), ndwi, profile)

# Step 7: Print min/max values
print()
print_min_max("NDVI", ndvi)
print_min_max("NDBI", ndbi)
print_min_max("NBR", nbr)
if ndwi is not None:
    print_min_max("NDWI", ndwi)
print_min_max("Brightness Temp (BT)", bt)
print_min_max("Emissivity", eps)
print_min_max("Land Surface Temperature (°C)", lst)

print("\nAll outputs saved. Processing complete.")
