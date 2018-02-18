# molecular-moviemaker
Scripts and functions for processing LCLS data using MATLAB to make Molecular Movies

## The schema of the data analysis is as follows:

### LargeRun
IS CALLED BY: none
CALLS: GetExpFina, fetch_timepoints, rdCSPADTS, goodshotsonly, getXRayWavelength, generatePedestal, image2rad, AnalysisFunc, universalmaskmaker
+ the big boss function that organizes everything
+ "photonimg_xxxx_unbinned_tot" summed from "photonimg_xxxx_unbinned", taken from AnalysisFunc (output as "photon_img")
+ then creates "photonimg_xxxx_unbinned_avg" by: "photonimg_xxxx_unbinned_tot" / "N_xxxxshots_counter" * "goodPixels" then subtracting a weighted vacuum

### image2rad
IS CALLED BY: LargeRun
CALLS: radialAverage
+ provides additional experimental parameters to radialAverage
+ takes "photonimg_xxxx_unbinned_avg" (from LargeRun) as input, provides "I_xxxx_unbinned" (output as I_q from radialAverage) as output

### radialAverage
IS CALLED BY: image2rad
CALLS: physicalQ, constructIq
+ uses "images" (photon counts at each pixel, input as "photon_avg") to create a radial average I_q (output from constructIq)

### physicalQ
IS CALLED BY: radialAverage
CALLS: none
+ uses detector geometry to convert physical space to q space 

### constructIq 
IS CALLED BY: radialAverage
CALLS: none
+ uses "image", aka "iFlat", which is an array of photon counts at each pixel
+ then applies corrections to provide I_q to radialAverage

### AnalysisFunc
IS CALLED BY: LargeRun
CALLS: rdCSPADdataXPP, CSPADtoMS, photon_map_hybrid
+ only input is from "in_struc"
+ applies common mode and gain corrections (to masks)

### photon_map_hybrid
IS CALLED BY: AnalysisFunc
CALLS: none
+ uses hybrid photon counting method (a minimum ADU threshold)
+ then time bins data to produce "binned_img"
