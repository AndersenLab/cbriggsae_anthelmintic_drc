#' Calculate EC10, EC50, and EC90 for dose-response data
#' 
#' @param data Cleaned dose-response data
#' @param drug_name Name of the drug/compound
#' @return List containing EC10, EC50, and EC90 data frames
#' @export
library(purrr)
library(tidyverse)
library(dplyr)
#library(knitr)
library(drc)

#'helper function to get length_delta csv (dependent on correct file paths and names)
#'
#'@param drug abbreviated drug name
#'@param file prefix of delta length file (either 20260617 or mean_median)
#'@return a data frame with strain, value (length_delta), and concentration
get_data <- function(drug, file){
  directory <- "projects/Olivia/cbriggsae_drc/data/"
  drc3 <- "20250602_Cbriggsae_DRC3_"
  drc4 <- "20250619_Cbriggsae_DRC4"
  drugs_drc3 <- c("ABZ", "CLO", "DEC", "EMO", "FBZ", "IVM", "LEV", "PYR")
  drugs_drc4 <- c("MOX","MILB")
  file_name <- paste0("/", file,"_length_delta_nemasize.csv")
  
  if (drug %in% drugs_drc3){
    pathway_drc <- paste0(directory,drc3,drug,file_name)
    }else{
    pathway_drc <- paste0(directory,drc4,drug, file_name)
  }
  
  
  mean_delta <- read.csv(pathway_drc)
  if (file == "20260617"){
    mean_delta <- mean_delta %>% 
      mutate(strain = ifelse(strain == "PB420", "AF16", strain))
    
    mean_delta <- as.data.frame(mean_delta) %>%
      dplyr::select(strain, value = median_wormlength_um_delta, 
                    concentration_um) %>%
      arrange(strain)
    return(mean_delta)
  }else{
    mean_delta <- mean_delta %>% 
       mutate(strain = ifelse(strain == "PB420", "AF16", strain))
    
     mean_delta <- as.data.frame(mean_delta) %>%
       dplyr::select(strain, value = mean_median_wormlength_um_delta, 
                  concentration_um) %>% 
       arrange(strain)
     
     
     
    return(mean_delta)
  }
 
}

#'helper function to get the AIC values for each drugs data
#'
#'@param data data frame with strain, value (length_delta), and concentration
#'@return a table with the AIC values for the following models:
#'LL2.4(), LL.5(), LL.4(), LL.3(), W1.3(), W1.4(), W2.4()
get_AIC <- function(data){
  
  drug_data <- data
  
  strain.model.fit <- drc::drm(data = drug_data, formula = value ~ concentration_um, curveid = strain, 
                               pmodels = list(~strain-1,  ~1, ~1, ~strain-1), 
                               fct = LL.4(fixed=c(NA, -800, NA, NA)))
  
  comparison <- drc::mselect(object = strain.model.fit, 
                             fctList = list(LL2.4(fixed = c(NA, -800, NA, NA)), 
                                            LL.5(fixed = c(NA, -800, NA, NA, NA)), 
                                            LL.3(fixed = c(-800, NA, NA)),
                                            W1.3(fixed = c(NA, -800, NA)), 
                                            W1.4(fixed = c(NA, -800, NA, NA)),
                                            W2.4(fixed = c(NA, -800, NA, NA))), 
                             sorted = "IC")
  return(comparison)
  
}

#'helper function for where to write data
#'@param drug abbreviated drug name
#'@param file_name name for the csv file you want to write
#'@param data_table the variable name for the data you want to save
write_data <- function(drug, file_name, data_table){
  drc3 <- "20250602_Cbriggsae_DRC3_"
  drc4 <- "20250619_Cbriggsae_DRC4"
  drugs_drc3 <- c("ABZ", "CLO", "DEC", "EMO", "FBZ", "IVM", "LEV", "PYR")
  drugs_drc4 <- c("MOX","MILB")
  
  if (drug %in% drugs_drc3){
    write_file <- paste0("projects/Olivia/cbriggsae_drc/data/",drc3,drug, "/", file_name, ".csv")
  }else{
    write_file <- paste0("projects/Olivia/cbriggsae_drc/data/",drc4,drug, "/", file_name, ".csv")
  }
  
  write.csv(data_table,file = write_file,row.names = TRUE)
  
}

###################################################
#                                                 #
# Identify which log-logistic model fits the best #
#                                                 #
###################################################

#using the length delta data to find the best fit model  
ABZ_data <- get_data("ABZ","20260617")
ABZ_aic <- get_AIC(drug_data <- ABZ_data)
#LL.5() with lowest IC
FBZ_data <- get_data("FBZ","20260617")
FBZ_aic <- get_AIC(drug_data <- FBZ_data)
#LL.5() with lowest IC
IVM_data <- get_data("IVM","20260617")
IVM_aic <- get_AIC(drug_data <- IVM_data)
#LL.3() with lowest IC
MILB_data <- get_data("MILB","20260617")
MILB_aic <- get_AIC(drug_data <- MILB_data)
#LL.5() with lowest IC
MOX_data <- get_data("MOX","20260617")
MOX_aic <- get_AIC(drug_data <- MOX_data)
#LL.5() with the lowest IC
LEV_data <- get_data("LEV","20260617")
LEV_aic <- get_AIC(drug_data <- LEV_data)
#LL.5() with the lowest IC
EMO_data <- get_data("EMO","20260617")
EMO_aic <- get_AIC(drug_data <- EMO_data)
#LL.5 with the lowest IC
CLO_data <- get_data("CLO","20260617")
CLO_aic <- get_AIC(drug_data <- CLO_data)
#LL.5 with the lowest IC
PYR_data <- get_data("PYR","20260617")
PYR_aic <- get_AIC(drug_data <- PYR_data)
#W1.4 with the lowest IC

# going ahead with LL.4 model since LL.5 is such a worse fit for PYR that we get a bunch of NANs#

################################
#                              #
# Calculating the EC estimates #
#                              #
################################

#' new function with new LL.5 model (did not end up using)
#' 
#' @param data data table with strain, value, and concentration columns
#' @returns output from drc::drm (if not possible returns "Unable to optimize model")
#safe_LL5_strain <- purrr::safely(.f = ~ drc::drm(data = ., 
#                                                 formula = value ~ concentration_um, strain,
#                                                pmodels=list(~strain-1,  ~1, ~1, ~strain-1),
#                                                 fct = LL.5(fixed=c(NA, -800, NA, NA, NA))),
#                                 otherwise = "Unable to optimize model")

#' Updated function with LL.4 model only thing changed is making lower limit NA instead of -800
#' (having issues with EC50 being identical between strains when lower limit was set)
#' 
#' @params data tabe with strain, value, and concentration columns
#' @returns output from drc::drm (if not possible returns "Unable to optomize model")
safe_LL4_strain <- purrr::safely(.f = ~ drc::drm(data = ., 
                                                 formula = value ~ concentration_um, 
                                                 curveid = strain,
                                                 pmodels=list(~strain-1,  ~1, ~1, ~strain-1),
                                                 fct = LL.4(fixed=c(NA, NA, NA, NA))),
                                 otherwise = "Unable to optimize model")

####################skip################
#didn't use
safe_slope_strain <- function(fit){
  if(is.character(fit$result)){print("Unable to extract slope")}
  else {
    #slopes<- data.frame(fit$result$coefficients[1:8]) %>% 8 strains for Sam? 
    slopes<- data.frame(fit$result$coefficients[1:6]) %>%
      # dplyr::rename(slope = `fit.result.coefficients.1.8.`) %>%8 strains for Sam? 
      dplyr::rename(slope = `fit.result.coefficients.1.6.`) %>%
      dplyr::mutate(strain = rownames(.),
                    strain = gsub(strain, pattern = "b:strain", replacement = ""))
    rownames(slopes) <- NULL
    return(slopes)
  }
}
#didn't use
safe_lower_asym <- function(fit){
  if(is.character(fit$result)){print("Unable to extract slope")}
  else {
    
    asyms <- data.frame(fit$result$coefficients[10:17]) %>% #data[10:17] = q90_worm_length, max_worm_length, cv_worm_length, n, metadata, metaexperiment, magnification, assay type
      dplyr::rename(lower.asymp = `fit.result.coefficients.10.17.`) %>%
      dplyr::mutate(strain = rownames(.),
                    strain = gsub(strain, pattern = "d:strain", replacement = ""))
    rownames(asyms) <- NULL
    return(asyms)
  }
}
#didn't use
safe_predict <- function(data, fit){
  if(is.character(fit$result)){print("Unable to optimize model")}
  else{
    expand.grid(conc = exp(seq(log(max(data$concentration_um)),
                               log(min(data$concentration_um[data$concentration_um!=min(data$concentration_um)])),
                               length = 100))) %>%
      dplyr::bind_cols(., tibble::as.tibble(predict(fit$result, newdata=., interval="confidence")))
  }
  
}
####################skip###############


#'Run the drc::ED function that estimates the effective doses
#'
#'@param fit output from safe_LL4_strain (model fitting from drc::drm function utilizing LL.4() model)
#'@return EC estimates calculated from drc::ED function (if possible if not returns "unable to optimize model")
safe_ED <- function(fit){
  if(is.character(fit$result)){print("Unable to optimize model")}
  else {
    drc::ED(object = fit$result, c(10,50,90), interval = "delta", pool = "FALSE", display = FALSE) #extract strain specific EC10, 50, 90 values
  }
}

#'Calculate and save the EC values for each drug
#'
#'@param drug abbreviated drug name (string)
#'@param drug_name full drug name for labeling (string)
#'@returns data table of EC values (calculated by drc::ED function)
calculate_one_EC <- function(drug, drug_name){

drug_data <- get_data(drug,"20260617")

strain.model.fit <- safe_LL4_strain(drug_data)

ec_values <- safe_ED(strain.model.fit)

ec_values <- as.data.frame(ec_values)

#add drug column
ec_values <- ec_values %>%
  dplyr::mutate(drug = drug_name)

#save ec_values table
#write_data(drug, "EC_estimates_LL.4", ec_values)

return(ec_values)
}

#run the calculate_one_EC for each drug for final table
ABZ_EC <- calculate_one_EC(drug = "ABZ", drug_name = "Albendazole")
FBZ_EC <- calculate_one_EC(drug = "FBZ", drug_name = "Fenbendazole")
CLO_EC <- calculate_one_EC(drug = "CLO", drug_name = "Closantel")
EMO_EC <- calculate_one_EC(drug = "EMO", drug_name = "Emodepside")
IVM_EC <- calculate_one_EC(drug = "IVM", drug_name = "Ivermectin")
LEV_EC <- calculate_one_EC(drug = "LEV", drug_name = "Levamisole")
PYR_EC <- calculate_one_EC(drug = "PYR", drug_name = "Pyrantel")
MOX_EC <- calculate_one_EC(drug = "MOX", drug_name = "Moxidectin")
MILB_EC <- calculate_one_EC(drug = "MILB", drug_name = "Milbemycin")

################reformat the EC estimations to get final csv##############
merged <- rbind(ABZ_EC, FBZ_EC, CLO_EC, EMO_EC, IVM_EC, LEV_EC, PYR_EC, MOX_EC, MILB_EC)
merged <- rownames_to_column(merged, var = "rowid")
# merged <- format(merged, scientific = FALSE)
merged <- merged %>% 
  group_by(rowid) %>% 
  summarise(Estimate = paste0(Estimate,' ± ',`Std. Error`), drug)

ECA2666 <- merged %>% 
  filter(str_detect(rowid,"ECA2666"))
ED3102 <- merged %>% 
  filter(str_detect(rowid,"ED3102"))
JU3237 <- merged %>% 
  filter(str_detect(rowid,"JU3237"))
NIC1667 <- merged %>% 
  filter(str_detect(rowid,"NIC1667"))
AF16 <- merged %>% 
  filter(str_detect(rowid,"AF16"))
VX34 <- merged %>% 
  filter(str_detect(rowid,"VX34"))

#'reformat the EC_vales data tables 
#'
#'@param strain_EC the EC data table from calculate one EC 
#'@param strain_name string of the strain name
#'@return organized data table without the merged row names
add_labels <- function(strain_EC,strain_name){
  strain_EC <- strain_EC %>% 
    mutate(strain = strain_name)
  EC_10 <- strain_EC %>% 
    filter(str_detect(rowid, paste0(strain_name, ":10"))) %>%
    mutate(EC = "EC10")
  EC_50 <- strain_EC %>% 
    filter(str_detect(rowid, paste0(strain_name, ":50"))) %>% 
    mutate (EC = "EC50")
  EC_90 <- strain_EC %>% 
    filter(str_detect(rowid, paste0(strain_name, ":90"))) %>% 
    mutate (EC = "EC90")
  strain_EC <- rbind(EC_10,EC_50, EC_90)
  strain_EC <- strain_EC %>% 
    dplyr::select(drug,strain,Estimate,EC) %>% 
    arrange(drug)
  return(strain_EC)
}

ECA2666 <- add_labels(ECA2666, "ECA2666")
ED3102 <- add_labels(ED3102, "ED3102")
JU3237 <- add_labels(JU3237, "JU3237")
NIC1667 <- add_labels(NIC1667, "NIC1667")
AF16 <- add_labels(AF16, "AF16")
VX34 <- add_labels(VX34, "VX34")

merged <- rbind(ECA2666, ED3102, JU3237, NIC1667, AF16, VX34)
merged <- pivot_wider(merged, names_from = strain, values_from = Estimate)
final_EC_table <- paste0("projects/Olivia/cbriggsae_drc/data/combined_data/total_EC_estimations_LL.4_final.csv")
write.csv(merged, file = final_EC_table, row.names = FALSE)
#################################################################################


#######################
#                     #
# Calculate the slope #
#                     #
#######################

#'function that takes the coefficients to use as the slope
#'
#'@param drug abbreviated drug name (string)
#'@param drug_name full drug name for labeling in final csv(string)
#'@return slope of the drug's DRC from the model fit summary coefficients
calculate_slope <- function(drug, drug_name){
  
  data <- get_data(drug, "20260617")
  strain.model.fit <- safe_LL4_strain(data)
  summarized.model.fit <- summary(strain.model.fit$result)
  strain.coefs <- data.frame(summarized.model.fit$coefficients)
  strain.coefs$params <- rownames(strain.coefs)
  rownames(strain.coefs) <- NULL
  
  slope.df <- strain.coefs %>%
    tidyr::separate(params, c("parameter","strain"), sep = ":") %>%
    dplyr::mutate(strain = gsub(strain, pattern = "strain", replacement = "")) %>%
    dplyr::filter(parameter == "b") %>%
    dplyr::rename(metric = parameter) %>%
    dplyr::mutate(drug = drug_name) %>%
    dplyr::select(Estimate, Std..Error, metric, strain, drug, p.value)
  
  return(slope.df)
  
}

#calculate the slope per strain
ABZ_slope <- calculate_slope(drug = "ABZ", drug_name = "Albendazole")
FBZ_slope <- calculate_slope(drug = "FBZ", drug_name = "Fenbendazole")
CLO_slope <- calculate_slope(drug = "CLO", drug_name = "Closantel")
EMO_slope <- calculate_slope(drug = "EMO", drug_name = "Emodepside")
IVM_slope <- calculate_slope(drug = "IVM", drug_name = "Ivermectin")
LEV_slope <- calculate_slope(drug = "LEV", drug_name = "Levamisole")
PYR_slope <- calculate_slope(drug = "PYR", drug_name = "Pyrantel")
MOX_slope <- calculate_slope(drug = "MOX", drug_name = "Moxidectin")
MILB_slope <- calculate_slope(drug = "MILB", drug_name = "Milbemycin")

#reformat the slopes together to get final csv 
total_slopes <- rbind(ABZ_slope, FBZ_slope, IVM_slope, MOX_slope, MILB_slope,
                      LEV_slope, PYR_slope, EMO_slope, CLO_slope)
total_slopes <- total_slopes %>% 
  mutate(Estimate = paste0(Estimate,' ± ',Std..Error)) %>%
  dplyr::select(drug, strain, Estimate)
total_slopes <- pivot_wider(total_slopes, names_from = strain, values_from = Estimate)
final_slope_table <- paste0("projects/Olivia/cbriggsae_drc/data/combined_data/total_slope_calculations_LL.4_final.csv")
write.csv(total_slopes, file = final_slope_table, row.names = FALSE)

 ######################
 #                    #
 # Stats Calculations #
 #                    #
 ######################

# F-test to find if variances for each strain-specific 
# dose-response model’s EC estimates were significantly different

#'function that used EDcomp to compare the effective concentrations 
#'between strains for the drug parameter (p_value < 0.05 bonferroni corrected)
#'
#'@param drug the abbreviation of the drug you want to look at 
#'@param drug_name the drug name you want to have in your outputted table
get_EC_comp <- function(drug, drug_name){
  data <- get_data(drug, "20260617")
  fit <- safe_LL4_strain(data)
  EDcomps.df <- suppressMessages(data.frame(drc::EDcomp(fit$result, c(50,50))))
  EDcomps.df$comps <- rownames(EDcomps.df)
  rownames(EDcomps.df) <- NULL
  ed.comps <- EDcomps.df %>%
    dplyr::mutate(drug = drug_name)
  ed.comps <- ed.comps %>% 
    mutate(p.value = p.adjust(p.value, method = "bonferroni"))
  print(ed.comps[1:5], row.names = FALSE)
  ed.comps <- ed.comps %>%
    mutate(sig_p_val = ifelse(p.value < 0.05, "*", ""))
  return (ed.comps)
}


#'function that used compParm to compare the slope 
#'between strains for the drug parameter (p_value < 0.05 bonferroni corrected)
#'
#'@param drug the abbreviation of the drug you want to look at 
#'@param drug_name the drug name you want to have in your outputted table
get_slope_comp <- function(drug, drug_name){
  data <- get_data(drug, "20260617")
  fit <- safe_LL4_strain(data)
  Slope.comps.df <- suppressMessages(data.frame(drc::compParm(fit$result, "b")))
  Slope.comps.df$comps <- rownames(Slope.comps.df)
  rownames(Slope.comps.df) <- NULL
  slope.comps <- Slope.comps.df %>%
    dplyr::mutate(drug = drug_name)
  slope.comps <- slope.comps %>% 
    mutate(p.value = p.adjust(p.value, method = "bonferroni"))
  print(slope.comps[1:5], row.names = FALSE)
  slope.comps <- slope.comps %>%
    mutate(sig_p_val = ifelse(p.value < 0.05, "*", ""))
  return (slope.comps)
}

ABZ_ec.comp <- get_EC_comp("ABZ", "Albendazole")
FBZ_ec.comp <- get_EC_comp("FBZ", "Fenbendazole")
IVM_ec.comp <- get_EC_comp("IVM", "Ivermectin")
MILB_ec.comp <- get_EC_comp("MILB", "Milbemycin")
MOX_ec.comp <- get_EC_comp("MOX", "Moxidectin")
LEV_ec.comp <- get_EC_comp("LEV", "Levamisole")
PYR_ec.comp <- get_EC_comp("PYR", "Pyrantel")
CLO_ec.comp <- get_EC_comp("CLO", "Closantel")
EMO_ec.comp <- get_EC_comp("EMO", "Emodepside")

ABZ_slope.comp <- get_slope_comp("ABZ", "Albendazole")
FBZ_slope.comp <- get_slope_comp("FBZ", "Fenbendazole")
IVM_slope.comp <- get_slope_comp("IVM", "Ivermectin")
MILB_slope.comp <- get_slope_comp("MILB", "Milbemycin")
MOX_slope.comp <- get_slope_comp("MOX", "Moxidectin")
LEV_slope.comp <- get_slope_comp("LEV", "Levamisole")
PYR_slope.comp <- get_slope_comp("PYR", "Pyrantel")
CLO_slope.comp <- get_slope_comp("CLO", "Closantel")
EMO_slope.comp <- get_slope_comp("EMO", "Emodepside")

total_slope_comps <- rbind(ABZ_slope.comp, FBZ_slope.comp, IVM_slope.comp, MOX_slope.comp, 
                           MILB_slope.comp, LEV_slope.comp, PYR_slope.comp, CLO_slope.comp, 
                           EMO_slope.comp)
total_slope_comps <- total_slope_comps %>% 
  dplyr::select(drug, comps, Estimate, Std..Error, t.value, p.value, sig_p_val)
final_slopecomp_table <- paste0("projects/Olivia/cbriggsae_drc/data/combined_data/total_slope_comparisons_LL.4_final.csv")
write.csv(total_slope_comps, file = final_slopecomp_table, row.names = FALSE)

total_ec_comps <- rbind(ABZ_ec.comp, FBZ_ec.comp, IVM_ec.comp, MOX_ec.comp, 
                        MILB_ec.comp, LEV_ec.comp, PYR_ec.comp, CLO_ec.comp, 
                        EMO_ec.comp)

total_ec_comps <- total_ec_comps %>% 
  dplyr::select(drug, comps, Estimate, Std..Error, t.value, p.value, sig_p_val)
final_eccomp_table <- paste0("projects/Olivia/cbriggsae_drc/data/combined_data/total_ec_comparisons_LL.4_final.csv")
write.csv(total_ec_comps, file = final_eccomp_table, row.names = FALSE)


#############################
#                           #
# Heritability Calculations #
#                           #
#############################

library(dplyr)
library(tidyr)
library(readr)
library(lme4)
library(sommer)
library(purrr)

#File path
geno_file <- "projects/Olivia/cbriggsae_drc/data/matrix/genotype_matrix.tsv"

geno_raw <- read.table(
  geno_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

geno_raw <- geno_raw %>% 
  dplyr::select(CHROM, POS, REF, ALT, ECA2666, AF16, ED3102, JU3237, NIC1667, VX34)

#' Calculate broad-sense heritability
#' 
#' @param d Data frame with columns: value (phenotype) and strain
#' @return Numeric value of H²
#' @export
calculate_H2 <- function(d) {
  
  # Fit linear mixed model with strain as random effect
  strain_fit <- lme4::lmer(data = d, formula = value ~ 1 + (1|strain))
  
  # Extract variance components
  variances <- lme4::VarCorr(x = strain_fit)
  variance_df <- as.data.frame(variances)
  
  # Genetic variance (among strains)
  Vg <- variance_df$vcov[1]
  
  # Environmental/residual variance
  Ve <- variance_df$vcov[2]
  
  # H² = Vg / (Vg + Ve)
  H2 <- Vg / (Vg + Ve)
  
  return(H2)
}

#' Calculate narrow-sense heritability
#' 
#' @param d Data frame with columns: value, strain
#' @param geno_matrix Genotype matrix for strains
#' @return Data frame with h², upper and lower bounds
#' @export
calculate_h2 <- function(d, geno_matrix) {
  
  pheno_strains <- unique(d$strain)
  # Calculate additive relationship matrix
  A <- sommer::A.mat(t(geno_matrix[, colnames(geno_matrix) %in% pheno_strains]))
  
  # Prepare phenotype data
  df_y <- d %>%
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value) %>%
    dplyr::mutate(strain = as.character(strain)) %>% 
    dplyr::mutate(strain = factor(strain, levels = unique(.$strain))) 
  
  random <- sommer::vsm((sommer::dsm(d$strain)), Gu = A)
  
  # Fit mixed model with genetic relationship matrix
  h2_res <- sommer::mmes(
    fixed = value ~ 1,
    random = ~ sommer::vsm((sommer::ism(strain)), Gu = A), 
    data = df_y
  )
  
  # Extract h² estimate and standard error
  h2 <- vpredict(h2_res, h2 ~ (V1) / (V1 + V2))[[1]][1]
  h2_SE <- vpredict(h2_res, h2 ~ (V1) / (V1 + V2))[[2]][1]
  
  h2_df <- data.frame(h2 = h2) %>%
    dplyr::mutate(
      h2.upper = h2 + h2_SE,
      h2.lower = h2 - h2_SE
    )
  
  return(h2_df)
}

#' Calculate heritability with bootstrapped confidence intervals
#' 
#' @param d Data frame with value and strain columns
#' @param geno_matrix Genotype matrix (for h² calculation)
#' @param nreps Number of bootstrap replicates (default 100)
#' @param boot Logical, whether to perform bootstrapping
#' @return List with H² CI, bootstrap values, and h² point estimate
#' @export
calculate_heritability_bootstrap <- function(d, geno_matrix, nreps = 100, boot = TRUE) {
  
  if(boot == TRUE) {
    
    # Point estimates
    H2_point <- calculate_H2(d = d)
    h2_point <- calculate_h2(d = d, geno_matrix = geno_matrix)
    
    H2_boots <- list()
    
    for(i in 1:nreps) {
      if(i %% 10 == 0) {
        message(paste0((i/nreps) * 100, "%"))
      }
      
      # Bootstrap within strain
      nested <- d %>%
        dplyr::group_by(strain) %>%
        tidyr::nest()
      
      boot_strain <- list()
      for(j in 1:length(nested$strain)) {
        boot_strain[[j]] <- nested$data[[j]][sample(seq(1:nrow(nested$data[[j]])), replace = TRUE), ] %>%
          dplyr::mutate(strain = nested$strain[[j]])
      }
      
      boot <- boot_strain %>%
        Reduce(rbind, .)
      
      # Check if at least 2 strains sampled
      check <- boot %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise(n = n())
      
      if(1 %in% check$n) {
        message("Only 1 strain sampled in bootstrap - skipping")
        next
      }
      
      # Calculate H² for bootstrap sample
      H2_est <- calculate_H2(d = boot)
      H2_boots[i] <- H2_est
    }
    
    H2_boots_vec <- unlist(H2_boots)
    H2_quantiles <- quantile(H2_boots_vec, probs = seq(0, 1, 0.05))
    
    H2_CI <- data.frame(
      H2.Point.Estimate = H2_point, 
      H2.5.perc = as.numeric(H2_quantiles[2]), 
      H2.95.perc = as.numeric(H2_quantiles[21])
    )
    
    return(list(H2_CI, H2_boots_vec, h2_point))
    
  } else {
    
    H2_point <- calculate_H2(d = d)
    H2_CI <- data.frame(
      H2.Point.Estimate = H2_point, 
      H2.5.perc = NA, 
      H2.95.perc = NA
    )
    
    return(H2_CI)
  }
}

#' Calculate heritability for each drug concentration
#' 
#' @param rep_data Data frame with replicates at one concentration
#' @param concentration_um Concentration value
#' @param geno_matrix Genotype matrix
#' @return List with H² and h² estimates
#' @export
calculate_dose_heritability <- function(rep_data, concentration_um, geno_matrix) {
  
  # Check for sufficient strain representation
  no_reps <- rep_data %>%
    dplyr::group_by(strain) %>%
    dplyr::summarise(n = n())
  
  if(1 %in% no_reps$n) {
    H2_est <- data.frame(
      H2 = NA, 
      lower = NA, 
      upper = NA
    ) %>%
      dplyr::mutate(
        drug = unique(rep_data$drug),
        concentration_um = concentration_um
      )
    
    h2_est <- data.frame(
      h2 = NA, 
      upper = NA, 
      lower = NA
    ) %>%
      dplyr::mutate(
        drug = unique(rep_data$drug),
        concentration_um = concentration_um
      )
    
    return(list(H2_est, h2_est))
    
  } else {
    
    herits <- suppressMessages(
      calculate_heritability_bootstrap(d = rep_data, geno_matrix = geno_matrix, boot = TRUE)
    )
    
    H2_est <- as.data.frame(herits[[1]]) %>%
      dplyr::mutate(
        drug = unique(rep_data$drug),
        concentration_um = concentration_um
      ) %>%
      dplyr::rename(
        H2 = H2.Point.Estimate,
        lower = H2.5.perc,
        upper = H2.95.perc
      )
    
    H2_replicates <- as.data.frame(herits[[2]]) %>%
      dplyr::mutate(
        drug = unique(rep_data$drug),
        concentration_um = concentration_um
      ) %>%
      dplyr::rename(H2.rep = `herits[[2]]`)
    
    h2_est <- as.data.frame(herits[[3]]) %>%
      dplyr::mutate(
        drug = unique(rep_data$drug),
        concentration_um = concentration_um
      ) %>%
      dplyr::rename(
        h2.upper = h2.upper,
        h2.lower = h2.lower
      )
    
    return(list(H2_est, H2_replicates, h2_est))
  }
}

#' Helper function to reformat calculate_one_EC output for compare_EC_to_heritable_dose and 
#' find_heritability_at_EC
#' 
#' @param ec_data output from calculate_one_EC
#' @return a reformatted ec table for the functions
format_ec <- function(ec_data){
  temp_merged <- rownames_to_column(ec_data, var = "rowid")
  ECA2666 <- temp_merged %>% 
    filter(str_detect(rowid,"ECA2666"))
  ED3102 <- temp_merged %>% 
    filter(str_detect(rowid,"ED3102"))
  JU3237 <- temp_merged %>% 
    filter(str_detect(rowid,"JU3237"))
  NIC1667 <- temp_merged %>% 
    filter(str_detect(rowid,"NIC1667"))
  AF16 <- temp_merged %>% 
    filter(str_detect(rowid,"AF16"))
  VX34 <- temp_merged %>% 
    filter(str_detect(rowid,"VX34"))
  
  ECA2666 <- add_labels(ECA2666, "ECA2666")
  ED3102 <- add_labels(ED3102, "ED3102")
  JU3237 <- add_labels(JU3237, "JU3237")
  NIC1667 <- add_labels(NIC1667, "NIC1667")
  AF16 <- add_labels(AF16, "AF16")
  VX34 <- add_labels(VX34, "VX34")
  
  temp_merged <- rbind(ECA2666, ED3102, JU3237, NIC1667, AF16, VX34)
  return(temp_merged)
}

#' Find the most heritable dose and compare to EC values
#' 
#' @param heritability_data Data frame with H² values at each concentration
#' @param ec_data Data frame with EC10, EC50, EC90 estimates
#' @param drug_name Name of the drug
#' @return Data frame comparing EC values to most heritable dose
#' @export
compare_EC_to_heritable_dose <- function(drug, raw_ec_data, drug_name) {
  drug_data <- get_data(drug,"20260617")
  
  ec_data <- format_ec(raw_ec_data)
  
  #nest the data by concentration values
  dose_nested <-  drug_data %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest()
  
  #calculate heritabilities
  heritabilities <- purrr::map2(
    dose_nested$data,
    dose_nested$concentration_um,
    ~calculate_dose_heritability(.x, .y, geno_raw))
  
  # Extract H² values
  heritability_data <- purrr::map(heritabilities, ~.x[[1]]) %>%
    Reduce(rbind, .)
  
  h2_data <- purrr::map(heritabilities, ~.x[[3]]) %>%
    Reduce(rbind, .)
  
  # Find concentration with maximum H²
  max_herit <- heritability_data %>%
    dplyr::filter(!is.na(H2)) %>%
    dplyr::arrange(desc(H2)) %>%
    dplyr::slice(1)
  
  most_heritable_conc <- max_herit$concentration_um
  max_H2 <- max_herit$H2
  
  # Calculate mean EC10, EC50, EC90
  mean_ec_values <- ec_data %>%
    dplyr::group_by(EC) %>%
    dplyr::summarise(mean_EC = mean(Estimate, na.rm = TRUE))
  
  # Find which EC is closest to most heritable dose
  ec_comparison <- mean_ec_values %>%
    dplyr::mutate(
      diff_from_heritable = abs(mean_EC - most_heritable_conc),
      most_heritable_conc = most_heritable_conc,
      max_H2 = max_H2,
      drug = drug_name
    ) %>%
    dplyr::arrange(diff_from_heritable)
  
  return(ec_comparison)
}

#' Find heritability at dose closest to each EC value
#' 
#' @param heritability_data Data frame with H² at each concentration
#' @param ec_data Data frame with EC estimates  
#' @param drug_name Name of drug
#' @return Data frame with combined H² and h² at doses closest to EC10, EC50, EC90
#' @export
find_heritability_at_EC <- function(drug, raw_ec_data, drug_name) {
  drug_data <- get_data(drug,"20260617")
  
  ec_data <- format_ec(raw_ec_data)
  
  #nest the data by concentration values
  dose_nested <-  drug_data %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest()
  
  #calculate heritabilities
  heritabilities <- purrr::map2(
    dose_nested$data,
    dose_nested$concentration_um,
    ~calculate_dose_heritability(.x, .y, geno_raw))
  
  # Extract H² values
  H2_data <- purrr::map(heritabilities, ~.x[[1]]) %>%
    Reduce(rbind, .)
  
  h2_data <- purrr::map(heritabilities, ~.x[[3]]) %>%
    Reduce(rbind, .)
  
  # Calculate mean EC values across strains
  mean_ecs <- ec_data %>%
    dplyr::group_by(EC) %>%
    dplyr::summarise(mean_EC = mean(Estimate, na.rm = TRUE))
  
  results_H2 <- list()
  results_h2 <- list()
  
  for(i in 1:nrow(mean_ecs)) {
    ec_level <- mean_ecs$EC[i]
    ec_value <- mean_ecs$mean_EC[i]
    
    # Find closest tested concentration to this EC
    herit_at_ec <- H2_data %>%
      dplyr::mutate(
        diff_from_EC = abs(concentration_um - ec_value)
      ) %>%
      dplyr::arrange(diff_from_EC) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        EC_level = ec_level,
        EC_value = ec_value,
        drug = drug_name
      )
    
    results_H2[[i]] <- herit_at_ec
  }
  
  for(i in 1:nrow(mean_ecs)) {
    ec_level <- mean_ecs$EC[i]
    ec_value <- mean_ecs$mean_EC[i]
    
    # Find closest tested concentration to this EC
    herit_at_ec <- h2_data %>%
      dplyr::mutate(
        diff_from_EC = abs(concentration_um - ec_value)
      ) %>%
      dplyr::arrange(diff_from_EC) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        EC_level = ec_level,
        EC_value = ec_value,
        drug = drug_name
      )
    
    results_h2[[i]] <- herit_at_ec
  }
  h2 <- as.data.frame(Reduce(rbind, results_h2))
  h2 <- h2 %>% 
    mutate(type = "h2") %>% 
    dplyr::select(drug, type, value = h2, EC_level, EC_value, upper = h2.upper, lower = h2.lower, 
           concentration_um, diff_from_EC)
  H2 <- as.data.frame(Reduce(rbind, results_H2))
  H2 <- H2 %>% 
    mutate(type = "H2") %>% 
    dplyr::select(drug, type, value = H2, EC_level, EC_value, upper = upper, lower = lower, 
                  concentration_um, diff_from_EC)
  
  return(rbind(h2,H2))
  
  print(Reduce(rbind, results_h2))
  print(Reduce(rbind, results_H2))
}

#' Complete comparison of EC values to heritability profile
#' 
#' @param cleaned_data Cleaned dose-response data
#' @param drug_name Name of drug
#' @param geno_matrix Genotype matrix
#' @return List with EC values, heritabilities, and comparison table
#' @export
complete_EC_heritability_analysis <- function(cleaned_data, drug_name, geno_matrix) {
  
  # 1. Calculate EC values
  ec_results <- calculate_all_ECs(cleaned_data, drug_name)
  
  # 2. Calculate heritability at each dose
  dose_nested <- cleaned_data %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest()
  
  heritabilities <- purrr::map2(
    dose_nested$data,
    dose_nested$concentration_um,
    ~calculate_dose_heritability(.x, .y, geno_matrix)
  )
  
  # Extract H² values
  H2_data <- purrr::map(heritabilities, ~.x[[1]]) %>%
    Reduce(rbind, .)
  
  h2_data <- purrr::map(heritabilities, ~.x[[3]]) %>%
    Reduce(rbind, .)
  
  # 3. Find most heritable dose
  max_herit_dose <- H2_data %>%
    dplyr::filter(!is.na(H2)) %>%
    dplyr::arrange(desc(H2)) %>%
    dplyr::slice(1)
  
  # 4. Compare each EC to most heritable dose
  all_ec_data <- bind_rows(
    ec_results$EC10 %>% dplyr::mutate(EC_type = "EC10"),
    ec_results$EC50 %>% dplyr::mutate(EC_type = "EC50"),
    ec_results$EC90 %>% dplyr::mutate(EC_type = "EC90")
  )
  
  ec_herit_comparison <- compare_EC_to_heritable_dose(
    H2_data, 
    all_ec_data, 
    drug_name
  )
  
  # 5. Find heritability at each EC
  herit_at_ecs <- find_heritability_at_EC(
    H2_data,
    all_ec_data,
    drug_name
  )
  
  return(list(
    EC_values = ec_results,
    heritabilities = H2_data,
    most_heritable_dose = max_herit_dose,
    EC_comparison = ec_herit_comparison,
    heritability_at_ECs = herit_at_ecs
  ))
}

#' Calculate the narrow-sense and broad-sense heritabilities for each concentration of the drug 
#' 
#' @param drug string abbreviation of the drug 
#' @param drug_name drug name to be used for final csv
#' @return a table with heritability calculations with lower and upper values sorted by conc 
get_hert_values <- function(drug, drug_name){
  drug_data <- get_data(drug,"20260617")
  
 #nest the data by concentration values
  dose_nested <-  drug_data %>%
    dplyr::group_by(concentration_um) %>%
    tidyr::nest()
  
  #calculate heritabilities
  heritabilities <- purrr::map2(
    dose_nested$data,
    dose_nested$concentration_um,
    ~calculate_dose_heritability(.x, .y, geno_raw))
  
  # Extract H² values
  H2_data <- purrr::map(heritabilities, ~.x[[1]]) %>%
    Reduce(rbind, .)
  
  h2_data <- purrr::map(heritabilities, ~.x[[3]]) %>%
    Reduce(rbind, .)
  
  max_herit_dose <- H2_data %>%
    dplyr::filter(!is.na(H2)) %>%
    dplyr::arrange(desc(H2)) %>%
    dplyr::slice(1)
  
  View(max_herit_dose)
  
  H2_data <- H2_data %>% 
    mutate(drug = drug_name) %>% 
    mutate(type = "broad-sense (H2)") %>% 
    dplyr::select(drug, concentration_um, type, value = H2, lower, upper)
  
  h2_data <- h2_data %>% 
    mutate(drug = drug_name) %>% 
    mutate(type = "narrow-sense (h2)") %>% 
    dplyr::select(drug, concentration_um, type, value = h2, lower = h2.lower, upper = h2.upper)
  
  total_hert <- rbind(H2_data, h2_data) %>% 
    arrange(concentration_um)
  
  return(total_hert)

}

#EC50 closest to most heritable dose (then EC10)
compare_EC_to_heritable_dose("ABZ", ABZ_EC, "Albendazole")
#EC90 closest by far (then EC 50)
compare_EC_to_heritable_dose("FBZ", FBZ_EC, "Fenbendazole")
#EC90 closest (then EC50) very small difference between them all
compare_EC_to_heritable_dose("IVM", IVM_EC, "Ivermectin")
#EC90 closest (then EC50) very small 
compare_EC_to_heritable_dose("MOX", MOX_EC, "Moxidectin")
#EC90 closest (then EC50)
compare_EC_to_heritable_dose("MILB", MILB_EC, "Milbemycin")
#EC50 closest (then EC10)
compare_EC_to_heritable_dose("LEV", LEV_EC, "Levamisole")
#EC50 closest (then EC10)
compare_EC_to_heritable_dose("PYR", PYR_EC, "Pyrantel")
#EC90 closest (then EC50)
compare_EC_to_heritable_dose("CLO", CLO_EC, "Closantel")
#EC50 closest (then EC10)
compare_EC_to_heritable_dose("EMO", EMO_EC, "Emodepside")


ABZ_ec.hert <- find_heritability_at_EC("ABZ", ABZ_EC, "Albendazole")
FBZ_ec.hert <- find_heritability_at_EC("FBZ", FBZ_EC, "Fenbendazole")
IVM_ec.hert <- find_heritability_at_EC("IVM", IVM_EC, "Ivermectin")
MOX_ec.hert <- find_heritability_at_EC("MOX", MOX_EC, "Moxidectin")
MILB_ec.hert <- find_heritability_at_EC("MILB", MILB_EC, "Milbemycin")
LEV_ec.hert <- find_heritability_at_EC("LEV", LEV_EC, "Levamisole")
PYR_ec.hert <- find_heritability_at_EC("PYR", PYR_EC, "Pyrantel")
CLO_ec.hert <- find_heritability_at_EC("CLO", CLO_EC, "Closantel")
EMO_ec.hert <- find_heritability_at_EC("EMO", EMO_EC, "Emodepside")

total_ec.hert <- rbind(ABZ_ec.hert, FBZ_ec.hert, IVM_ec.hert, MILB_ec.hert, 
                       MOX_ec.hert, LEV_ec.hert, PYR_ec.hert, CLO_ec.hert,
                       EMO_ec.hert)
final_ec.hert_table <- "projects/Olivia/cbriggsae_drc/data/combined_data/total_ec.herit_calculations.csv"
write.csv(total_ec.hert, file = final_ec.hert_table, row.names = FALSE)

ABZ_hert <- get_hert_values("ABZ", "Albendazole")
FBZ_hert <- get_hert_values("FBZ", "Fenbendazole")
IVM_hert <- get_hert_values("IVM", "Ivermectin")
MOX_hert <- get_hert_values("MOX", "Moxidectin")
MILB_hert <- get_hert_values("MILB", "Milbemycin")
LEV_hert <- get_hert_values("LEV", "Levamisole")
PYR_hert <- get_hert_values("PYR", "Pyrantel")
EMO_hert <- get_hert_values("EMO", "Emodepside")
CLO_hert <- get_hert_values("CLO", "Closantel")

EMO.h <- EMO_hert %>%
  filter(!EMO_hert$concentration_um == 0)
mean(EMO.h$value)

total_hert <- rbind(ABZ_hert, FBZ_hert, IVM_hert, MOX_hert, MILB_hert,
                    LEV_hert, PYR_hert, EMO_hert, CLO_hert)
final_hert_table <- "projects/Olivia/cbriggsae_drc/data/combined_data/total_herit_calculations.csv"
write.csv(total_hert, file = final_hert_table, row.names = FALSE)

