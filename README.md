# Content description
R scripts and full simulation results for manuscript "The impact of exposure measurement error on the estimated concentration-response relationship between long-term exposure to PM2.5 and mortality".

# R scripts (in directory "RLib")
The grid-level PM2.5 data are publicly available at https://sedac.ciesin.columbia.edu/data/set/aqdh-pm2-5-concentrations-contiguous-us-1-km-2000-2016. The ZIP Code-level PM2.5 and uncertainty data are available from the corresponding author on reasonable request (weiyg@hsph.harvard.edu). Covariate data are publicly available with sources described in the manuscript. Data availability that relates to predictors in the PM2.5 prediction model can be referred to _Di Q, Amini H, Shi L, Kloog I, Silvern R, Kelly J, et al. An ensemble-based model of PM2.5 concentration across the contiguous United States with high spatiotemporal resolution. Environ Int. 2019;130:104909_.

 - Exposure uncertainty is estimated with script "Step1_EstimateUncertainty.R".
 - Spatial correlation information is obtained with script "Step2_EstimateSpatialCorrelation.R".
 - Simulation datasets are generated with script "Step3_GenerateSimulationData.R".
 - Simulation under linear concentration-response model is performed with script "Step4_Simulation_Linear.R".
 - Simulation under linear concentration-response model at low exposure levels is performed with script "Step4_Simulation_Linear_LowLevel.R".
 - Simulation under quadratic concentration-response model is performed with scripts "Step4_Simulation_Quadratic_Quadratic.R" and "Step4_Simulation_Quadratic_Penalized.R".
 - Simulation under soft-threshold concentration-response model is performed with script "Step4_Simulation_Softplus.R".
 - Results are compiled with script "Step5_CompileResults.R"

# Simulation results (in directory "Results")



