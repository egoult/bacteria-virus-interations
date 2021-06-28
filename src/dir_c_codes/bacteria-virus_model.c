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
X_SS = 1 - C0 - E0;
X_SE = E0;
X_SI = 0;
X_SR = 0;

// colonised by bacteria
X_CS = C0;
X_CE = 0;
X_CI = 0;
X_CR = 0;

// bacteria cases
Y_C = 0;

// virus cases
Y_I = 0;

// virus and bacteria coinfection
Y_CI = 0;
//end_rinit

//start_skel

double prev_bac, prev_vir; // Prevalence of bacteria colonisation, virus infection
double lambda[4], beta[2]; // bacteria transmission, virus transmission
double foi_bac[4], foi_vir[2]; // force of infection for bacteria and virus 

// calculate the prevalence
prev_bac = v1 * X_CS + v2 * X_CE + v3 * X_CI + v4 * X_CR ;
prev_vir = w1 * X_SI + w2 * X_CI ;

// might be worth expressing beta, lambda in terms of R0bac R0 vir etc
lambda[0] = R0_bac/delta1; // bacteria transmission
lambda[1] = R0_bac/delta2; // bacteria transmission
lambda[2] = R0_bac/delta3; // bacteria transmission
lambda[3] = R0_bac/delta4; // bacteria transmission


beta[0] = R0_vir/gamma1;   // virus transmission
beta[1] = R0_vir/gamma2;   // virus transmission

foi_bac[0] = prev_bac * lambda[0]; //(lambda is a 4 element vector)
foi_bac[1] = prev_bac * lambda[1]; //(lambda is a 4 element vector)
foi_bac[2] = prev_bac * lambda[2]; //(lambda is a 4 element vector)
foi_bac[3] = prev_bac * lambda[3]; //(lambda is a 4 element vector)

foi_vir[0] = prev_vir * beta[0];  //(beta is a 2 element vector)
foi_vir[1] = prev_vir * beta[1];  //(beta is a 2 element vector)


// Derivatives
// bacteria susceptible
DX_SS = - X_SS * foi_bac[0] - X_SS * foi_vir[0] + delta1 * X_CS;
DX_SE =   X_SS * foi_vir[0] - X_SE * foi_bac[1] + delta2 * X_CE - gamma1 * X_SE;
DX_SI =   gamma1 * X_SE     - X_SI * foi_bac[2] + delta3 * X_CI - sigma1 * X_SI ;
DX_SR =   sigma1 * X_SI     - X_SR * foi_bac[3] + delta4 * X_CR;

// bacteria colonised
DX_CS = - X_CS * foi_vir[1] - delta1 * X_CS + foi_bac[0] * X_SS;
DX_CE =   X_CS * foi_vir[1] - delta2 * X_CE + foi_bac[1] * X_SE - gamma2 * X_CE;
DX_CI =   gamma2 * X_CE     - delta3 * X_CI + foi_bac[2] * X_SI - sigma2 * X_CI;
DX_CR =   sigma2 * X_CI     - delta4 * X_CR + foi_bac[3] * X_SR;

// Accumulator variables
// bacteria cases
DY_C = alpha1 * X_CS + alpha2 * X_CE + alpha3 * X_CI + alpha4 * X_CR;

// virus cases
DY_I = eta1 * gamma1 * X_SE + eta2 * gamma2 * X_CE;

// virus and bacteria coinfection
DY_CI = alpha3 * eta2 * gamma2 * theta * X_CE;

//end_skel

//start_rsim
//end_rsim
