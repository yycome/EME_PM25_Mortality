# Content description
R scripts for manuscript "The impact of exposure measurement error on the estimated concentration-response relationship between long-term exposure to PM2.5 and mortality".

# R scripts (in directory "RLib")
 - Exposure uncertainty is estimated with script "__Step1_EstimateUncertainty.R__".
 - Spatial correlation information within the range of 0-40km is obtained with script "__Step2_EstimateSpatialCorrelation.R__".
 - Spatial correlation information within the range of 0-100km is obtained with script "__Step2_EstimateSpatialCorrelation_Spatial100km.R__".
 - Simulation datasets are generated with script "__Step3_GenerateSimulationData.R__".
 - Simulation datasets for sensitivity analysis with respect to the range of spatial correlation between 0-100km are generated with script "__Step3_GenerateSimulationData_Spatial100km.R__".
 - Simulation under linear concentration-response model is performed with script "__Step4_Simulation_Linear.R__".
 - Simulation under linear concentration-response model at low exposure levels is performed with script "__Step4_Simulation_Linear_LowLevel.R__".
 - Simulation under linear concentration-response model for sensitivity analysis with respect to the range of spatial correlation between 0-100km is performed with script "__Step4_Simulation_Linear_Spatial100km.R__".
 - Simulation under quadratic concentration-response model is performed with scripts "Step4_Simulation_Quadratic_Quadratic.R" and "__Step4_Simulation_Quadratic_Penalized.R__".
 - Simulation under soft-threshold concentration-response model is performed with script "__Step4_Simulation_Softplus.R__".
 - Results are compiled with script "__Step5_CompileResults.R__"
