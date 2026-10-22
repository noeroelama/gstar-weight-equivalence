# Data

All files cover the same 40 regions in the same order: the 38 provinces of Indonesia and two
boundary regions, Papua New Guinea (`PNG_ALL`) and the non-Indonesian part of Borneo
(Sabah, Sarawak and Brunei, `BORNEO_NONIDN`). The boundary regions enter the spatial lags;
forecast errors are computed for the 38 Indonesian provinces only. There are no missing values.

## rainfall_monthly.csv

Monthly rainfall (mm/month), January 1985 to April 2025 (484 months). The first column is the
month (`YYYY-MM-DD`, first day of the month); the other columns are regions, keyed by
`region_id`.

Each value is the mean of all CHIRPS v2.0 pixels (0.05 degree) that touch the region polygon
(rasterised with `all_touched=True`), taken from the global monthly CHIRPS v2.0 NetCDF file.

- CHIRPS v2.0: Climate Hazards Center, UC Santa Barbara,
  https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_monthly/netcdf/
  Funk, C., Peterson, P., Landsfeld, M., et al. (2015). The climate hazards infrared
  precipitation with stations: a new environmental record for monitoring extremes.
  *Scientific Data* 2, 150066. https://doi.org/10.1038/sdata.2015.66
- Province polygons: "Provinsi" layer (38 provinces, 2023) of
  https://github.com/Alf-Anas/batas-administrasi-indonesia. According to that repository, the
  boundaries come from the GIS portal of the Directorate-General of Population and Civil
  Registration (Dukcapil), Ministry of Home Affairs of Indonesia.
- Boundary regions: Natural Earth 1:110m admin-0 countries (public domain),
  https://www.naturalearthdata.com

The polygons are not included in this repository.

## climate_indices.csv

Monthly Niño 3.4 anomaly (`nino34`, HadISST) and Dipole Mode Index (`dmi`, HadISST) for the
same months, from the NOAA Physical Sciences Laboratory, https://psl.noaa.gov/data/timeseries/month/

## regions.csv

| column | meaning |
|---|---|
| `region_id` | ISO 3166-2:ID code for Indonesian provinces; `PNG_ALL`, `BORNEO_NONIDN` for the boundary regions |
| `country` | `IDN`, `PNG` or `OTH` |
| `name` | English region name |
| `kode_prov` | two-digit province code (Kepmendagri No. 300.2.2-2138 of 2025); empty for the boundary regions |
| `lat`, `lon` | polygon centroid in degrees (WGS 84) |

Centroids were computed in the Equal Earth equal-area projection (EPSG:8857) and converted back
to latitude and longitude. The inverse-distance matrix of the paper uses great-circle
(haversine) distances between these centroids.

## knn_weights.csv

The geographic nearest-neighbour matrix of the paper: for each region, the five nearest regions
by haversine distance between centroids, with weights proportional to inverse distance and each
row summing to one (zero diagonal). Rows and columns are in the order of `regions.csv`. The two
data-driven matrices (training correlation and graphical lasso) are computed by
`code/weights.py` from the first 240 months.

## Terms of use

CHIRPS data are in the public domain. The rainfall series in this folder are released under the
Creative Commons Attribution 4.0 International license (CC BY 4.0); please cite Funk et al.
(2015) as well. The centroid coordinates and the nearest-neighbour matrix are derived from the
polygon layer named above, which does not state a license; they are provided only so that the
results can be reproduced, and users should credit that source.
