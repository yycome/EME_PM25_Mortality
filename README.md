# Content description
R scripts and full simulation results for manuscript "The impact of exposure measurement error on the estimated concentration-response relationship between long-term exposure to PM2.5 and mortality".

# R scripts (in directory 'Rcodes')
The grid-level PM2.5 data are publicly available at https://sedac.ciesin.columbia.edu/data/set/aqdh-pm2-5-concentrations-contiguous-us-1-km-2000-2016. The ZIP Code-level PM2.5 are available from the corresponding author on reasonable request (weiyg@hsph.harvard.edu). Covariate data are publicly available with sources described in the manuscript. Data sources for PM2.5 modeling can be referred to _Di Q, Amini H, Shi L, Kloog I, Silvern R, Kelly J, et al. An ensemble-based model of PM2.5 concentration across the contiguous United States with high spatiotemporal resolution. Environ Int. 2019;130:104909_.

 - The exposure uncertainty is estimated with script "Step1_EstimateUncertainty.R".
 - The spatial correlation information is obtained with script "Step2_EstimateSpatialCorrelation.R".
 - The simulation datasets are generated with script "Step3_GenerateSimulationData.R".
 - Results are compiled with script "Step5_CompileResults.R"



