##
## Create the pomp object for the virus x bacteria model
## Author: Elizabeth Goult
## Notes: All units are in days
##        Based off code by Matthieu Domenech de Cellès
##

CreateSCxSIRmod<-function(){

  
  # Extract C code from file
  mod_code <- readLines("src/dir_c_codes/SCxSIRmodel.c")
  components_nm <- c("skel", "rinit")#"globs"
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for(nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0("start_", nm)):str_which(mod_code, paste0("end_", nm))] %>% 
      str_flatten(collapse = "\n")
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  po <- pomp(data = data.frame(time = seq(from = 0, to = 400, by = 0.01), X = NA),
             times = "time",
             t0 = 0,
             obsnames = "X",
             statenames = c("X_SS", "X_SI", "X_SR", 
                            "X_CS", "X_CI", "X_CR", 
                            "Y_bac", "Y_vir", "Y_bac_vir"),
             accumvars = c("Y_bac", "Y_vir", "Y_bac_vir"),
             paramnames = c("cstar",         "R0_vir",
                            "delta",          "gamma",
                            "rho",            "zeta",
                            "eta_bac",        "eta_vir", 
                            "theta_vir_beta", "theta_vir_lambda",
                            "theta_bac_eta",  "theta_vir_eta",                                                      
                            "I0",             "C0" ),
             params = c(cstar = 0.2,         R0_vir = 2.5,
                        delta = 1/50,        gamma = 1/6,
                        rho = 0,             zeta = 1/3,
                        eta_bac = 15/100000, eta_vir = 1, 
                        theta_vir_beta = 1,  theta_vir_lambda = 1,
                        theta_bac_eta = 1,   theta_vir_eta = 1,                                                      
                        I0 = 1e-8,           C0 = 0.2 ),
            
             #globals = components_l[["globs"]],
             skeleton = vectorfield(components_l[["skel"]]),
             rinit = components_l[["rinit"]],
             verbose = F
             #cdir = "_C", cfile = "codes" 
  )
  return(po)
}