##
## check model
## Author: Elizabeth Goult
## Comments:

# clear workspace 
rm(list=ls())
graphics.off()

#libraries
require(pomp)
require(tidyverse)
require(wesanderson)
require(patchwork)

source("src/CreateSCxSIRmod.R")

# save data
if(!dir.exists("results/")){
    dir.create("results/")
}

if(!dir.exists("plots/")){
    dir.create("plots/")
}


# functions-------------------------------------
SimulateAtPrevalence<-function(po, cstar_val){
    # For a given value of bacteria prevalence 
    # return peak, overall and peak timing of virus
    # po - pomp object
    # cstar_value - bacterial prevalence

    stopifnot("cstar_val must be less than 1" = (cstar_val<1))
    
    coef(po, "cstar")<-cstar_val

    # simulate
    tj0<-trajectory(po, format= "data.frame")

    # calculate peak incidence/ timing/ overall incidence
    peak<- max(tj0$Y_vir)
    peak_timing<-tj0[tj0$Y_vir==peak,"time"]
    overall<- sum(tj0$Y_vir)
    
    return(c(peak = peak, time = peak_timing, overall = overall))}

PrevalenceGradient<-function(po, step = 0.01){
    #  Computes peak, time of peak and overall prevalence of vir
    # for a gradient of  bacterial prevalences
    # po - pomp obeject
    # step - step sizer of gradient

    prev_grad<- seq(0, 1-step, by = step)

    grad_res<-sapply(prev_grad, SimulateAtPrevalence, po = po) %>% t() 
    grad_res<-data.frame(grad_res, prevalence=prev_grad)

    return(grad_res)}


# create pomp object ---------------------------
po<-CreateSCxSIRmod()

sim_check<- SimulateAtPrevalence(po, cstar_val=0.1)

# gradient bacterial prevalences
base_po<-CreateSCxSIRmod()
grad_res<-PrevalenceGradient(po = base_po)

# plot the  gradient
plt<-ggplot(grad_res, aes(x = prevalence, y = round(peak, digits = 8))) +
    geom_line()+
    theme_classic()+
    labs(x = "Bacteria prevalence", y = "Virus peak prevalence")


# scenarios ---------------------------------
theta_vir_beta<-c(0.8, 1, 1.2)# virus aquisition
theta_vir_lambda<-c(0.8, 1, 1.2) # virus - transmission
theta_vir_eta<-c(0.8, 1, 1.2) # virus reporting rate

scenarios<-data.frame(acquisition = theta_vir_beta, transmission = theta_vir_lambda, reporting = theta_vir_eta)
rm(list = c("theta_vir_beta", "theta_vir_eta", "theta_vir_lambda"))

# calculate for each scenario
scenario_list<-list()
counter<-0
for(acq in scenarios$acquisition ){
    for(trans in scenarios$transmission){
        for(rep in scenarios$transmission){ 
            # storing the results
            counter<-counter + 1           
            
            # changing the parameter values
            pars<-c(theta_vir_beta = acq, theta_vir_lambda = trans, theta_vir_eta = rep)
            po_alt<-base_po
            coef(po_alt, names(pars))<-pars
            
            # simulate for the gradient
            scenario_list[[counter]]<-PrevalenceGradient(po = po_alt, step = 0.01) %>%
                mutate(theta_vir_beta = acq, theta_vir_lambda = trans, theta_vir_eta = rep, sim = counter)
            

        }
    }
}

scenarios<-bind_rows(scenario_list)
write.csv(scenarios, "results/bacteria-virus_scenarios.csv", row.names = FALSE)

# plot all together
plt.scenarios<-ggplot(scenarios, aes(x = prevalence, y = peak, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Viral prevalence")

# print(plt.scenarios)

plt.scen.grid<-plt.scenarios +
    facet_wrap(~sim)

# plot with parameters as color ------------------------------------
# peak
plt.beta.peak<-ggplot(scenarios, aes(x = prevalence, y = peak, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Peak viral prevalence", color = "Acquisition")+
    facet_grid(theta_vir_lambda ~theta_vir_eta, labeller = label_both)


plt.eta.peak<-ggplot(scenarios, aes(x = prevalence, y = peak, group = sim, color = as.factor(theta_vir_eta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Peak viral prevalence", color = "Severity")+
    facet_grid(theta_vir_lambda ~theta_vir_beta, labeller = label_both)

plt.lambda.peak<-ggplot(scenarios, aes(x = prevalence, y = peak, group = sim, color = as.factor(theta_vir_lambda)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Peak viral prevalence", color = "Transmission")+
    facet_grid(theta_vir_beta ~theta_vir_eta, labeller = label_both)

#time
plt.beta.time<-ggplot(scenarios, aes(x = prevalence, y = time, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Timing of peak viral prevalence", color = "Aquisition")+
    facet_grid(theta_vir_lambda ~theta_vir_eta, labeller = label_both)

plt.eta.time<-ggplot(scenarios, aes(x = prevalence, y = time, group = sim, color = as.factor(theta_vir_eta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Timing of peak viral prevalence", color = "Severity")+
    facet_grid(theta_vir_lambda ~theta_vir_beta, labeller = label_both)

plt.lambda.time<-ggplot(scenarios, aes(x = prevalence, y = time, group = sim, color = as.factor(theta_vir_lambda)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Timing of peak viral prevalence", color = "Transmission")+
    facet_grid(theta_vir_beta ~theta_vir_eta, labeller = label_both)

# overall
plt.beta.overall<-ggplot(scenarios, aes(x = prevalence, y = overall, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Overall viral prevalence", color = "Acquisition")+
    facet_grid(theta_vir_lambda ~theta_vir_eta, labeller = label_both)

plt.eta.overall<-ggplot(scenarios, aes(x = prevalence, y = overall, group = sim, color = as.factor(theta_vir_eta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Overall viral prevalence", color = "Severity")+
    facet_grid(theta_vir_lambda ~theta_vir_beta, labeller = label_both)

plt.lambda.overall<-ggplot(scenarios, aes(x = prevalence, y = overall, group = sim, color = as.factor(theta_vir_lambda)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Overall viral prevalence", color = "Transmission")+
    facet_grid(theta_vir_beta ~theta_vir_eta, labeller = label_both)


# difference between null and results -----------------
scen_diff<-scenarios %>%
    group_by(sim) %>%
    mutate(peak_diff = peak - grad_res$peak, time_diff = time - grad_res$time, overall_diff = overall - grad_res$overall)

plt.diff.peak<-ggplot(scen_diff, aes(x = prevalence, y = peak_diff, group = sim, color = as.factor(sim)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in peak viral prevalence", color = "Simulation")

plt.diff.time<-ggplot(scen_diff, aes(x = prevalence, y = time_diff, group = sim, color = as.factor(sim)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in timing of peak viral prevalence", color = "Simulation")

plt.diff.overall<-ggplot(scen_diff, aes(x = prevalence, y = overall_diff, group = sim, color = as.factor(sim)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in overall viral prevalence", color = "Simulation")

# facet plots
#peak
plt.diff.peak.beta<-ggplot(scen_diff, aes(x = prevalence, y = peak_diff, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in peak viral prevalence", color = "Acquisition")+
    facet_grid(theta_vir_eta~theta_vir_lambda, as.table = FALSE, labeller = label_both)

plt.diff.peak.lambda<-ggplot(scen_diff, aes(x = prevalence, y = peak_diff, group = sim, color = as.factor(theta_vir_lambda)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in peak viral prevalence", color = "Transmission")+
    facet_grid(theta_vir_eta~theta_vir_beta, as.table = FALSE, labeller = label_both)

plt.diff.peak.eta<-ggplot(scen_diff, aes(x = prevalence, y = peak_diff, group = sim, color = as.factor(theta_vir_eta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in peak viral prevalence", color = "Severity")+
    facet_grid(theta_vir_beta~theta_vir_lambda, as.table = FALSE, labeller = label_both)


#time
plt.diff.time.beta<-ggplot(scen_diff, aes(x = prevalence, y = time_diff, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in timing of peak viral prevalence", color = "Acquistion")+
    facet_grid(theta_vir_eta~theta_vir_lambda, as.table = FALSE, labeller = label_both)

plt.diff.time.lambda<-ggplot(scen_diff, aes(x = prevalence, y = time_diff, group = sim, color = as.factor(theta_vir_lambda)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in timing of peak viral prevalence", color = "Transmission")+
    facet_grid(theta_vir_eta~theta_vir_beta, as.table = FALSE, labeller = label_both)

plt.diff.time.eta<-ggplot(scen_diff, aes(x = prevalence, y = time_diff, group = sim, color = as.factor(theta_vir_eta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in timing of peak viral prevalence", color = "Severity")+
    facet_grid(theta_vir_beta~theta_vir_lambda, as.table = FALSE, labeller = label_both)


# overall
plt.diff.overall.beta<-ggplot(scen_diff, aes(x = prevalence, y = overall_diff, group = sim, color = as.factor(theta_vir_beta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in overall viral prevalence", color = "Acquisition")+
    facet_grid(theta_vir_eta~theta_vir_lambda, as.table = FALSE, labeller = label_both)

plt.diff.overall.lambda<-ggplot(scen_diff, aes(x = prevalence, y = overall_diff, group = sim, color = as.factor(theta_vir_lambda)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in overall viral prevalence", color = "Transmission")+
    facet_grid(theta_vir_eta~theta_vir_beta, as.table = FALSE, labeller = label_both)

plt.diff.overall.eta<-ggplot(scen_diff, aes(x = prevalence, y = overall_diff, group = sim, color = as.factor(theta_vir_eta)))+
    geom_line()+
    theme_classic()+
    labs(x = "Bacterial prevalence", y = "Difference in overall peak viral prevalence", color = "Severity")+
    facet_grid(theta_vir_beta~theta_vir_lambda, as.table = FALSE, labeller = label_both)

# pdf("plots/prevalence_gradient.pdf", height = 7.5, width = 22.5)
#     print(plt.beta.peak + plt.lambda.peak + plt.eta.peak)
#     print(plt.beta.time + plt.lambda.time + plt.eta.time)
#     print(plt.beta.overall + plt.lambda.overall + plt.eta.overall)
# dev.off()

# pdf("plots/prevalence_gradient_diff.pdf", height = 7.5, width = 22.5)
#     print(plt.diff.peak.beta + plt.diff.peak.lambda + plt.diff.peak.eta)
#     print(plt.diff.time.beta + plt.diff.time.lambda + plt.diff.time.eta)
#     print(plt.diff.overall.beta + plt.diff.overall.lambda + plt.diff.overall.eta)
# dev.off()

#plot each interaction alone and combined
scen_plot<- scenarios %>%
    filter( (theta_vir_beta == 1 & theta_vir_eta == 1 ) |
            (theta_vir_beta == 1 & theta_vir_lambda == 1 ) |
            (theta_vir_lambda == 1 & theta_vir_eta == 1 ) |
            (theta_vir_beta == 0.8 & theta_vir_eta == 0.8 & theta_vir_lambda == 0.8) |
            (theta_vir_beta == 1.2 & theta_vir_eta == 1.2 & theta_vir_lambda == 1.2)) %>%
            rowwise() %>%
            mutate(parameter = (theta_vir_beta+theta_vir_eta+theta_vir_lambda)/3)%>%
            mutate(parameter = case_when(parameter== max(c(theta_vir_beta, theta_vir_eta, theta_vir_lambda))~parameter,
                                    parameter> 1 ~ max(c(theta_vir_beta, theta_vir_eta, theta_vir_lambda)),
                                    parameter< 1 ~ min(c(theta_vir_beta, theta_vir_eta, theta_vir_lambda)))) %>%
            mutate(parameter_name = case_when((theta_vir_beta == theta_vir_eta & theta_vir_eta== theta_vir_lambda) ~ "All parameters" , 
                                               theta_vir_beta == parameter ~ "Acquisition",
                                               theta_vir_eta == parameter ~ "Severity", 
                                               theta_vir_lambda == parameter ~ "Transmission")) %>%
            filter(!(theta_vir_beta == 1 & theta_vir_eta == 1 & theta_vir_lambda == 1 ))
                    

# plot to check
plt<-ggplot(scen_plot %>% filter(between(prevalence,0.05,0.6)), 
            aes(x = prevalence, y = overall, group = sim, color = sim))+
            geom_line()

plt.peak<-ggplot(scen_plot %>% filter(between(prevalence,0.05,0.6)), 
                aes(x = prevalence, y = peak, group = sim, color = parameter_name))+
                geom_line(size = 1)+
                geom_line( data = scenarios  %>% 
                                 filter(between(prevalence,0.05,0.6)) %>% 
                                 filter(theta_vir_beta == 1 & theta_vir_eta == 1 & theta_vir_lambda == 1 ), 
                                 color = "grey", linetype = "dashed", alpha = 0.7) + 
                facet_grid(~parameter)+
                theme_classic()+
                labs(x = "Bacterial prevalence", y = "Peak viral incidence", color = "Parameter varied")


plt.time<-ggplot(scen_plot %>% filter(between(prevalence,0.05,0.6)), 
                aes(x = prevalence, y = time, group = sim, color = parameter_name))+
                geom_line(size = 1)+
                geom_line( data = scenarios  %>% 
                                 filter(between(prevalence,0.05,0.6)) %>% 
                                 filter(theta_vir_beta == 1 & theta_vir_eta == 1 & theta_vir_lambda == 1 ), 
                                 color = "grey", linetype = "dashed", alpha = 0.7) + 
                facet_grid(~parameter)+
                theme_classic()+
                labs(x = "Bacterial prevalence", y = "Time of peak viral incidence", color = "Parameter varied")


plt.overall<-ggplot(scen_plot %>% filter(between(prevalence,0.05,0.6)), 
                aes(x = prevalence, y = overall, group = sim, color = parameter_name))+
                geom_line(size = 1)+
                geom_line( data = scenarios  %>% 
                                 filter(between(prevalence,0.05,0.6)) %>% 
                                 filter(theta_vir_beta == 1 & theta_vir_eta == 1 & theta_vir_lambda == 1 ), 
                                 color = "grey", linetype = "dashed", alpha = 0.7) + 
                facet_grid(~parameter)+
                theme_classic()+
                labs(x = "Bacterial prevalence", y = "Overall viral incidence", color = "Parameter varied")


# pdf("plots/prevalence_gradient_updated.pdf")
#     print(plt.peak / plt.time / plt.overall)
# dev.off()



scen_plot<- scenarios %>%
    filter(theta_vir_eta == 1) %>%
    filter( (theta_vir_beta == 1 ) |
            (theta_vir_lambda == 1 ) |            
            (theta_vir_beta == 0.8 & theta_vir_lambda == 0.8) |
            (theta_vir_beta == 1.2 & theta_vir_lambda == 1.2)) %>%
            rowwise() %>%
            mutate(parameter = (theta_vir_beta+theta_vir_lambda)/2)%>%
            mutate(parameter = case_when(parameter== max(c(theta_vir_beta, theta_vir_lambda))~parameter,
                                    parameter> 1 ~ max(c(theta_vir_beta, theta_vir_lambda)),
                                    parameter< 1 ~ min(c(theta_vir_beta, theta_vir_lambda)))) %>%
            mutate(parameter_name = case_when((theta_vir_beta == theta_vir_lambda) ~ "Acquisition and Transmission" , 
                                               theta_vir_beta == parameter ~ "Acquisition",
                                            #    theta_vir_eta == parameter ~ "Severity", 
                                               theta_vir_lambda == parameter ~ "Transmission")) %>%
            filter(!(theta_vir_beta == 1 & theta_vir_lambda == 1 ))
            

independent_scen<- scenarios  %>% 
                    filter(theta_vir_beta == 1 & theta_vir_eta == 1 & theta_vir_lambda == 1 ) %>%
                    select(prevalence, peak) %>% 
                    mutate(indep_peak = peak) %>%
                    select(!peak)

scen_plot_rescaled <- scen_plot %>%
                        full_join(independent_scen) %>%
                        mutate(rescaled_peak = peak/indep_peak) %>%
                        mutate(parameter = case_when( parameter ==  0.8 ~ "Antagonistic", 
                                                      parameter == 1.2  ~ "Synergistic"))

plt.peak<-ggplot(scen_plot_rescaled %>% filter(between(prevalence,0.05,0.6)), 
                aes(x = prevalence, y = rescaled_peak, group = sim, color = parameter_name))+
                geom_line(size = 0.9, alpha = 0.5)+
                geom_point(data = scen_plot_rescaled %>% 
                                  filter(between(prevalence,0.05,0.6), prevalence == round(prevalence, 1)),                                  
                                  shape = 19, color = "white", size = 7) +
                geom_point(data = scen_plot_rescaled %>% 
                                  filter(between(prevalence,0.05,0.6), prevalence == round(prevalence, 1)),                                  
                                  aes(shape = as.factor(parameter)), alpha = 10, size = 7) +
                geom_line(aes(y = indep_peak/ indep_peak), color = "grey", linetype = "solid", alpha = 0.5) + 
                # scale_alpha_manual(values = c(0.5, 1))+
                # scale_linetype_manual(values = c("solid", "solid"))+
                scale_color_manual(values = wes_palette(name = "Darjeeling1", n = 5)[c(1,3,5)] )+
                scale_shape_manual(values = c( "_", "+"))+
                # facet_grid(~parameter)+
                theme_classic()+
                labs(x = "Bacterial prevalence", y = "Relative Peak viral incidence", color = "Interaction", shape = "Direction")

pdf("plots/prevalence_gradient_updated2.pdf", height = 10, width = 7.5)
    print(plt.peak)
dev.off()
