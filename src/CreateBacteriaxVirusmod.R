##
## Create the pomp object for the virus x bacteria model
## Author: Elizabeth Goult
## Notes: All units are in days
##        Based off code by Matthieu Domenech de Cellès
##

CreateBacteriaxVirusmod<-function(){

  
  # Extract C code from file
  mod_code <- readLines("src/dir_c_codes/bacteria-virus_model.c")
  components_nm <- c("skel", "rinit")#"globs"
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for(nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0("start_", nm)):str_which(mod_code, paste0("end_", nm))] %>% 
      str_flatten(collapse = "\n")
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  po <- pomp(data = data.frame(time = seq(from = 0., to = 400, by = 0.05), X = NA),
             times = "time",
             t0 = 0,
             obsnames = "X",
             statenames = c("X_SS", "X_SE", "X_SI", "X_SR", 
                            "X_CS", "X_CE", "X_CI", "X_CR", 
                            "Y_C", "Y_I", "Y_CI"),
             accumvars = c("Y_C", "Y_I", "Y_CI"),
             paramnames = c("R0_bac", "R0_vir",
                            "w1", "w2", 
                            "v1", "v2", "v3", "v4", 
                            "delta1", "delta2", "delta3", "delta4", 
                            "gamma1", "gamma2",                             
                            "sigma1", "sigma2",
                            "alpha1", "alpha2", "alpha3", "alpha4",
                            "eta1", "eta2",                           
                            "E0", "C0", "theta"), 
             params = c(R0_bac = 0.6, R0_vir = 3,
                        w1=1, w2=1, 
                        v1=1, v2=1, v3=1, v4=1, 
                        delta1 = 1/50 , delta2 = 1/50, delta3 = 1/50, delta4 = 1/50, 
                        gamma1 = 1/4, gamma2 = 1/4,
                        sigma1 = 1/5, sigma2 = 1/5,
                        # alpha1 = 15/100000, alpha2 = 15/100000, alpha3 = 15/100000, alpha4 = 15/100000,
                        alpha1 = 1, alpha2 = 1, alpha3 = 1, alpha4 = 1,
                        # eta1 = 0.7, eta2 = 0.7,
                        eta1 = 1, eta2 = 1,
                        E0 = 0.25, C0 = 0.25, theta = 1), 
             #globals = components_l[["globs"]],
             skeleton = vectorfield(components_l[["skel"]]),
             rinit = components_l[["rinit"]],
             verbose = F
             #cdir = "_C", cfile = "codes" 
  )
  return(po)
}