# LAC_Coastal_Population
Coastal Population estimates for Latin America and Caribbean region
## Introduction

We produced first detailed population estimates of the number and percentage of the population living near the coast in Latin America and the Caribbean. These estimates are stratified into 1, 5, and 10km zones as these distances are helpful to understand at regional level the concentration of population living by the coast. 

To get consistent and comparable results running the analysis across the region, we followed the same geographic and population data input and standard workflow throughout. 

Visualising information spatially tells a more powerful and compelling story when compared to the same data being hidden in statistical tables. 

## Methodology

The methodology has been scripted in R so that the analysis can be easily reproduced and improved when a new version of the modelled population estimates become available.

The path for the analysis was to generate 1, 5 and 10 km buffers around the coastal lines and then, by overlapping the modelled gridded population from WorldPop with these Buffer Coastal Areas, estimate the population living within these coastal areas.

Geographic projections have been a challenge to deal with when carrying out the geospatial analysis across the entire LAC region due to the necessity of running the GIS related operations in a wide latitude range. Using standard Latitude and longitude projection EPSG 4326 WGS 84, does not represent a big issue if the analysis was made on latitudes close to the Equator, where area and distance distortions are not very important. However, the use of this projection in countries such as Chile or Argentina negatively affects the results as area and distance measurement are compromised by the projection’s distortion. In the sections below we explain the projection related processing methods.


More detailed information is available in the technical report.
