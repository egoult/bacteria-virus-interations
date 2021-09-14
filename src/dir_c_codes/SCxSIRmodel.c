// Implementing a bacteria-virus model(SEIR x SC)


//start_globs
//end_globs

//start_globs
//end_globs

//start_toest
//end_toest

//start_fromest
//end_fromest

//start_dmeas
//end_dmeas

//start_rmeas
//end_rmeas

//start_rinit

// susceptible to bacteria
X_SS = 1 - C0 - I0;
X_SI = I0;
X_SR = 0;

// colonised by bacteria
X_CS = C0;
X_CI = 0;
X_CR = 0;

// bacteria cases
Y_bac = 0;

// virus cases
Y_bac = 0;

// virus and bacteria coinfection
Y_bac_vir = 0;
//end_rinit

//start_skel
double beta_bac, beta_vir; // bacteria transmission, virus transmission
double prev_bac, prev_vir; // Prevalence of bacteria colonisation, virus infection
double foi_bac, foi_vir; // force of infection for bacteria and virus 

// calculate the prevalence
prev_bac = X_CS + X_CI + X_CR ;
prev_vir = X_SI + theta_vir_lambda * X_CI ;

// transmission

beta_bac = delta/(1-cstar) ; 
beta_vir = R0_vir * gamma; 

// force of infection
foi_bac = prev_bac * beta_bac; 
foi_vir = prev_vir * beta_vir; 


// Derivatives
// bacteria susceptible
DX_SS = - X_SS * foi_vir - X_SS * foi_bac + delta * X_CS;
DX_SI =   X_SS * foi_vir - X_SI * foi_bac + (( 1-rho)*delta + rho * zeta) * X_CI - gamma * X_SI ;
DX_SR =  gamma * X_SI    - X_SR * foi_bac + delta * X_CR;

// bacteria colonised
DX_CS = - theta_vir_beta * X_CS * foi_vir - delta * X_CS + foi_bac * X_SS;
DX_CI =   theta_vir_beta * X_CS * foi_vir - delta * X_CI + foi_bac * X_SI - gamma * X_CI;
DX_CR =   gamma * X_CI - delta * X_CR + foi_bac * X_SR;

// Accumulator variables
// bacteria cases
DY_bac = eta_bac * gamma *(X_CS + theta_bac_eta * X_CI + X_CR);

// virus cases
DY_vir = eta_vir * gamma *(X_SI + theta_vir_eta * X_CI);

// virus and bacteria coinfection
DY_bac_vir = eta_bac * eta_vir * theta_bac_eta * theta_vir_eta * gamma * X_CI;

//end_skel

//start_rsim
//end_rsim
