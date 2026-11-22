# Hamersley Basin map data (public)

Files copied unchanged from map2loop (Loop3D), `map2loop/_datasets/geodata_files/hamersley/`,
commit 078f8a6dd39794796ac608622469f299fe963bc5 (31 January 2026), MIT licence
(`LICENSE_map2loop.txt`, The Loop Project).

- `geology.geojson`: geological map polygons (field `unitname`), EPSG:28350 (GDA94 / MGA zone 50).
- `faults.geojson`: mapped fault traces.
- `structures.geojson`: bedding measurements (fields `dip`, `azimuth2` = dip direction).
- `dtm_rp.tif`: elevation model, float32 GeoTIFF.

The map data originate from the Geological Survey of Western Australia (1:500 000 State
interpreted bedrock geology). Reference: Jessell, M. et al. (2021), Automated geological map
deconstruction for 3D model construction using map2loop 1.0 and map2model 1.0, Geoscientific
Model Development 14, 5063-5092.
