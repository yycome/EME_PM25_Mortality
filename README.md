# Content description
R scripts for manuscript "The impact of exposure measurement error on the estimated concentration-response relationship between long-term exposure to PM2.5 and mortality".

# R scripts (in directory "RLib")
 - Exposure uncertainty is estimated with script "Step1_EstimateUncertainty.R".
 - Spatial correlation information is obtained with script "Step2_EstimateSpatialCorrelation.R".
 - Simulation datasets are generated with script "Step3_GenerateSimulationData.R".
 - Simulation under linear concentration-response model is performed with script "Step4_Simulation_Linear.R".
 - Simulation under linear concentration-response model at low exposure levels is performed with script "Step4_Simulation_Linear_LowLevel.R".
 - Simulation under quadratic concentration-response model is performed with scripts "Step4_Simulation_Quadratic_Quadratic.R" and "Step4_Simulation_Quadratic_Penalized.R".
 - Simulation under soft-threshold concentration-response model is performed with script "Step4_Simulation_Softplus.R".
 - Results are compiled with script "Step5_CompileResults.R"
