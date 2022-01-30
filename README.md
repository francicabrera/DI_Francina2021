# DI_Francina2021
## "Extracting urban informal settlements using machine learning"
### A pixel-based approach to classify informal Settlements in the Dominican Republic's capital using PlanetScope data.

Author: Francina Cabrera Fermin

CASA UCL 2021


This document provides instructions to extract informal settlements using spectral and texture features from 3 metres resolution satellite scenes from PlanetScope. This project creates several Random Forest Classification model based on the spectral and gray level co-occurrence matrix (variance) values of satellite image pixels.

#### Description of the data:
1. Two PlanetScope Analytic Ortho scenes with 3 metres pixel resolution, along with their Unusable Data Masks (UDM). Taken on the 23/02/2021 23:20:17. Downloaded from: https://www.Planet.com
2. Visually collected labelled training and assessment sample using VHR Google maps imagery.
3. Distrito Nacional’s boundary (.shp) collected in the 2010 Population Census by National Statistics Office of the Dominican Republic. Downloaded from: https://www.one.gob.do/informaciones-cartograficas/shapefiles
4. Distrito Nacional’s circumscription 3 (.shp) collected in the 2010 Population Census by National Statistics Office of the Dominican Republic. Downloaded from: https://www.one.gob.do/informaciones-cartograficas/shapefiles
5. Dominican Republic’s provinces (.shp) collected in the 2010 Population Census by National Statistics Office of the Dominican Republic. Downloaded from: https://www.one.gob.do/informaciones-cartograficas/shapefiles

The data is stored in an online github repository: https://github.com/francicabrera/DI_Francina2021.git

Find the process in this path: Slum_detection/Slum_detection.html
Or visit the rpubs publication: 

This code follows the process from:
https://pages.cms.hu-berlin.de/EOL/gcg_eo/index.html
