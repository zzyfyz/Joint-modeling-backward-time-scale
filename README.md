# Joint-modeling-backward-time-scale
This is the code of simulation for the manuscript titled "Joint modeling of quality of life and survival using a Bayesian approach in a retrospective time scale"
Authors: Yizhou Fei*, Elizabeth Juarez-Colunga, Areej El-Jawahri, Jean S. Kutner, Kathryn Colborn.

Data and scripts:
The Dataset folder contains two simulated datasets generated under different baseline hazard assumptions:
1. Weibuill baseline hazard
2. Piecewise-exponential baseline hazard
The file Terminal_run.R is the main script used to fit the proposed joint model.

Required inputs
Before running the script, you will need to prepare the following:
1. Analysis dataset
A subject-level dataset containing the variables required for model fitting.
2. Missingness mask (mask)
A separate dataset indicating which longitudinal outcome values are observed versus missing. This is required because Stan does not allow NA values in the input data.

User-specified settings
The following model and MCMC settings must be specified by the user before fitting the model:
1. Number of internal knots for the spline function
2. Spline degree
3. Total number of MCMC iterations, including warmup
4. Number of MCMC chains

Notes
These choices can substantially affect both computation time and model performance. In particular, increasing the spline complexity or the number of MCMC iterations may lead to considerably longer run times. Users should select these settings carefully based on their data, inferential goals, and available computing resources.


