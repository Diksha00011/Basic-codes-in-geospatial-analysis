import requests

url = "https://raw.githubusercontent.com/codeforamerica/click_that_hood/master/public/data/california-counties.geojson"

response = requests.get(url)

with open("california_counties.geojson", "wb") as f:
    f.write(response.content)

print("Download complete.")


import geopandas as gpd
import matplotlib.pyplot as plt

file_path = "california_counties.geojson"

counties = gpd.read_file(file_path)

print("Data loaded successfully!")
print(counties.head())

san_diego = counties[counties["name"] == "San Diego"]

fig, ax = plt.subplots(figsize=(8,8))

counties.plot(ax=ax, color="lightgray", edgecolor="white")
san_diego.plot(ax=ax, color="red")

plt.title("San Diego County - California")
plt.show()
san_diego["area_sqkm"] = san_diego.geometry.area / 10**6
print("Area in sq km:", san_diego["area_sqkm"])
