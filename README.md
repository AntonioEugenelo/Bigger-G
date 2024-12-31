# Bigger-G
Repository with provisional MATLAB and Dynare code for endogenous fiscal policy in a Big G style DSGE model 

Starting with the original model, as shown in https://www.nber.org/papers/w27034, I assume government expenditure is:
1. At, or close to, market prices (and therefore included, at least partly, in the firm maximization problem)
2. Endogenous in its sectoral distribution under an exogenous total government expenditure

Since the IRFs from the model rely on pivotal parameters, namely the intrasectoral consumer substitutability of goods (theta) and the CES aggregator for government expenditure (rho), the code is organized as follows:

1. Simulation_loglin.mod holds the Dynare model, with parameters a, theta, and rho expressed as global variables with standard default values.
2. Optimize_loglin.m holds the Matlab code iterating the Dynare model for all suggested values of theta and rho, for each combination of which it will determine an optimal a (policy parameter), it then saves the IRFs for the optimal values
3. The IRFs figure files contain the saved results for each optimal policy 

The figures are sorted into folders by the parameter mu, representing the ratio of government spending internalized into the Firm's price-setting problem (0<mu<1).

Further branches indicate possible alterations to the Objective function:
1. Basic quadratic objective in output gap and inflation
2. Log-linearized consumer utility (default)
3. Aggregate utility internalizing a positive effect of government expenditure, as shown in https://www.nber.org/papers/w32914 
