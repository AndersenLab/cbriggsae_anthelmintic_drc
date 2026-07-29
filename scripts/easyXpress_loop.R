library(tidyverse)
library(dplyr)
library(knitr)

# To go through a drug step by step to see the variables
###########################################
#                                         #
#                 NEMASIZE                #
#                                         #
###########################################

# (3 bleach) drugs in dose response 2 (ABZ, CLO, DEC, EMO, FBZ, IVM, LEV, MOX, PYR)
# (1 bleach) drugs in dose response 3 (ABZ, CLO, DEC, EMO, FBZ, IVM, LEV, MOX, PYR)
# drugs in dose response 4 (moxidectin and milbemycin mix)

# a: "20250424_Cbriggsae_DRC2","20250602_Cbriggsae_DRC3", "20250619_Cbriggsae_DRC4"
# b: "ABZ", "CLO", "DEC", "EMO", "FBZ", "IVM", "LEV", "MOX", "PYR", "MILB"

a <- "20250619_Cbriggsae_DRC4"
b <- "MILB"
combine_drc <- FALSE
data_pathway <- "projects/Amanda/projects/nemasize/projects/"

if(combine_drc){
  design_file_path_DRC2 <- paste0(data_pathway,"20250424_Cbriggsae_DRC2_",b,"/design/20240424_designfile_",b, ".csv")
  data_file_path_DRC2 <-  paste0(data_pathway,"20250424_Cbriggsae_DRC2_",b, "/NemaSize_output/skeleton/worm_lengths.csv")
  design_file_path_DRC3 <- paste0(data_pathway,"20250602_Cbriggsae_DRC3_",b,"/design/20250602_CbDRC3_",b, "_design.csv")
  data_file_path_DRC3 <- paste0(data_pathway,"20250602_Cbriggsae_DRC3_",b, "/NemaSize_output/skeleton/worm_lengths.csv")
}else if(a == "20250424_Cbriggsae_DRC2"){
  design_file_path <- paste0(data_pathway,a,"_",b,"/design/20240424_designfile_",b, ".csv")
  message(glue::glue("Here is the design file path: {design_file_path}"))
  data_file_path <-  paste0(data_pathway,a,"_",b, "/NemaSize_output/skeleton/worm_lengths.csv")
  message(glue::glue("Here is the data file path: {data_file_path}"))
} else if(a == "20250602_Cbriggsae_DRC3"){
  design_file_path <- paste0(data_pathway,a,"_",b,"/design/20250602_CbDRC3_",b, "_design.csv")
  message(glue::glue("Here is the design file path: {design_file_path}"))
  data_file_path <- paste0(data_pathway,a,"_",b, "/NemaSize_output/skeleton/worm_lengths.csv")
  message(glue::glue("Here is the data file path: {data_file_path}"))
} else{
  design_file_path <- paste0(data_pathway,a,"/design/20250619_designfile_DRC4.csv")
  message(glue::glue("Here is the design file path: {design_file_path}"))
  data_file_path <- paste0(data_pathway,a,"/NemaSize_output/skeleton/worm_lengths.csv")
  message(glue::glue("Here is the data file path: {data_file_path}"))
}

if (!combine_drc){
  message(glue::glue("***Running easyXpress on {a} looking at drug {b}***"))
  message(glue::glue("--Retrieving the design and data file--")) 
  # Define experimental directory and file name
  design_Cb_ns <- read.csv(design_file_path)
  datafile_Cb_ns <- read.csv(data_file_path)
  
  #check to see if there are the right number of bleaches (problem with having more than 3 groups)
  num_bleaches <- length(unique(design_Cb_ns$bleach))
  if (num_bleaches > 4){
    design_Cb_ns <- design_Cb_ns %>% mutate(bleach = str_sub(design_Cb_ns$bleach, 1,1))
  }
} else{
  message(glue::glue("***Running easyXpress on combined DRC2 and DRC3 looking at drug {b}***"))
  message(glue::glue("--Retrieving the design and data files--")) 
  
  # Define experimental directory and file name
  datafile_Cb_ns_DRC2 <- read.csv(data_file_path_DRC2)
  design_Cb_ns_DRC2 <- read.csv(design_file_path_DRC2)
  design_Cb_ns_DRC3 <- read.csv(design_file_path_DRC3)
  datafile_Cb_ns_DRC3 <- read.csv(data_file_path_DRC3)
  
  num_bleaches_combine <- length(unique(design_Cb_ns_DRC2$bleach))
  
  #change DRC3 bleach to 4 
  design_Cb_ns_DRC3 <- design_Cb_ns_DRC3 %>% mutate(bleach = 4)
  
  #combine the two design and data file to one of each 
  design_Cb_ns <- bind_rows(design_Cb_ns_DRC2, design_Cb_ns_DRC3)
  datafile_Cb_ns <- bind_rows(datafile_Cb_ns_DRC2, datafile_Cb_ns_DRC3)
  
  #check to see if there are the right number of bleaches (problem with having more than 3 groups)
  num_bleaches <- length(unique(design_Cb_ns$bleach))
  if (num_bleaches > 4){
    design_Cb_ns <- design_Cb_ns %>% mutate(bleach = str_sub(design_Cb_ns$bleach, 1,1))
  }
}


#drug names need to be have matching characters 
design_Cb_ns <- design_Cb_ns %>% mutate(drug = toupper(drug))

#separate design.csv based on drug for DRC4
if (b == "MOX" && a == "20250619_Cbriggsae_DRC4"){
  design_Cb_ns <- design_Cb_ns[design_Cb_ns$drug == "MOXIDECTIN",]
}

if (b == "MILB" && a == "20250619_Cbriggsae_DRC4"){
  design_Cb_ns <- design_Cb_ns[design_Cb_ns$drug == "MILBEMYCIN",]
}

if (b == "FBZ"){
  design_Cb_ns <- design_Cb_ns[!(design_Cb_ns$concentration_um == 30.000),]
}

design_Cb_ns$bleach <- as.character(design_Cb_ns$bleach)

# Join data + design
merged_Cb_ns <- datafile_Cb_ns %>%
  left_join(design_Cb_ns, by = c("Metadata_Plate", "Metadata_Well")) 

#get rid of NA drug columns (in DRC4 case)
merged_Cb_ns <- merged_Cb_ns[!is.na(merged_Cb_ns$drug),]
merged_Cb_ns <- merged_Cb_ns[!is.na(merged_Cb_ns$concentration_um),]

merged_Cb_ns <- as.data.frame(merged_Cb_ns)

message("--merged the design and merged data frame")
# first summarize to well stats - number of worms per well
well_summary_ns <- merged_Cb_ns %>%
  group_by(Metadata_Plate, Metadata_Well, bleach, strain, concentration_um, drug, diluent) %>%
  summarise(
    n = n(),  #counts worms per well
    .groups = "drop"
  )

  # ugh - need to fix titerWF its not going to affect stuff downstream
  # titerWF - new function for nemasize
  titerWF_custom <- function(data, ..., thresh = 0.68, plot = TRUE, doseR = FALSE) {
    
    # Expecting message
    if(doseR == TRUE){
      message("You set doseR = TRUE. Expecting controls to be coded as for a dose response.")
    }
    if(doseR == FALSE){
      message("You set doseR = FALSE. Not expecting controls to be coded for a dose response.")
    }
    
    # Check on control coding
    controls <- unique(data$diluent)
    drugs <- unique(data$drug)
    
    # Check on expected controls
    if(doseR == FALSE & FALSE %in% (controls %in% drugs)) {
      message(glue::glue("WARNING: the controls are not configured as expected in the design file doseR = FALSE. Do you want doseR = TRUE? If not, control conditions should have the same value for drug and diluent and a 0 for concentration_um."))
    }
    
    # Warn about doseR=T
    if(doseR == TRUE & TRUE %in% (controls %in% drugs)){
      message(glue::glue("WARNING: the controls are not configured as expected in the design file for doseR = TRUE. Do you want doseR = FALSE? If not, control conditions should have 0 for concentration_um and be named for the drug not the diluent."))
    }
    
    grouping_vars <- enquos(...)
    
    if(doseR == FALSE) {
      # Set the flags for all independent bleaches
      d <- data %>%
        dplyr::select(Metadata_Plate, Metadata_Well, !!!grouping_vars, n, drug, concentration_um, diluent) %>%
        dplyr::filter(drug == diluent & concentration_um == 0) %>%
        dplyr::group_by(!!!grouping_vars) %>%
        dplyr::mutate(num.wells = dplyr::n(),
                      cv.n = sd(n)/mean(n)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(quant.95.cv.n = quantile(cv.n, .95, na.rm = TRUE),
                      quant.97.5.cv.n = quantile(cv.n, .975, na.rm = TRUE),
                      titer_WellFlag = ifelse(cv.n > thresh, "titer", NA_character_))
    }
    
    if(doseR == TRUE) {
      # Set the flags for all independent bleaches
      d <- data %>%
        dplyr::select(Metadata_Plate, Metadata_Well, !!!grouping_vars, n, drug, concentration_um, diluent) %>%
        dplyr::filter(concentration_um == 0) %>%
        dplyr::group_by(!!!grouping_vars) %>%
        dplyr::mutate(num.wells = dplyr::n(),
                      cv.n = sd(n)/mean(n)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(quant.95.cv.n = quantile(cv.n, .95, na.rm = TRUE),
                      quant.97.5.cv.n = quantile(cv.n, .975, na.rm = TRUE),
                      titer_WellFlag = ifelse(cv.n > thresh, "titer", NA_character_))
    }
    
    # Get the wells to flag with titer
    flags <- d %>%
      dplyr::filter(titer_WellFlag == "titer") %>%
      tidyr::unite(flag_id, !!!grouping_vars, sep = "_", remove = FALSE, na.rm = TRUE) %>%
      dplyr::distinct(flag_id) %>%
      dplyr::pull(flag_id)
    
    # Set the flag by flag.id
    t.d <- data %>%
      tidyr::unite(flag_id, !!!grouping_vars, sep = "_", remove = FALSE, na.rm = TRUE) %>%
      dplyr::mutate(titer_WellFlag = dplyr::case_when(flag_id %in% flags ~ "titer",
                                                      TRUE ~ NA_character_))
    
    # Get distinct bleaches for plotting and reporting
    d.p <- d %>%
      dplyr::distinct(!!!grouping_vars, .keep_all = TRUE)
    
    # Message
    message(glue::glue("{nrow(d.p)} independent bleaches detected. The titer_WellFlag is set in the output data."))
    
    if(plot == TRUE) {
      # Plot it with thresh
      p <- ggplot2::ggplot(d.p) +
        ggplot2::aes(x = cv.n, fill = strain) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::geom_vline(xintercept = thresh, linetype = 2, color = "red") +
        ggplot2::theme_bw() +
        ggplot2::labs(title = glue::glue("Titer filter cv.n > {thresh}"),
                      y = "bleach count")
      
      # Return
      message(glue::glue("A diagnostic plot for checking cv.n threshold is returned."))
      out <- list(d = t.d, p = p)
      return(out)
    } else {
      # Return data only
      return(t.d)
    }
  }
  
  message("--using the custom titerWF function -> titer_Cb_ns --")
  
  # Use titerWF
  titer_Cb_ns <- titerWF_custom(data = well_summary_ns,
                                bleach, strain, concentration_um,
                                thresh = 0.68,
                                plot = TRUE,
                                doseR = TRUE)
  
  message("titerWF data --> titer_Cb_data_ns | titerWF plot --> titer_Cb_plot_ns")
  
  titer_Cb_data_ns <- titer_Cb_ns$d
  titer_Cb_plot_ns <- titer_Cb_ns$p
  
  # nWF
  nWF_custom <- function(data, ..., max = 75, min = 5, plot = T) {
    # add flag
    d <- data %>%
      dplyr::mutate(n_WellFlag = dplyr::case_when(n > max ~ paste0("n>", max),
                                                  n < min ~ paste0("n<", min),
                                                  n >= min & n <= max ~ NA_character_,
                                                  TRUE ~ "ERROR"))
    # message
    message("The n_WellFlag is set in the output data.")
    
    if(plot == T) {
      # plot it with thresh
      p <- ggplot2::ggplot(d) +
        ggplot2::aes(x = n) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::geom_vline(xintercept = max, linetype = 2, color = "red") +
        ggplot2::geom_vline(xintercept = min, linetype = 2, color = "red") +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(ggplot2::vars(...)) +
        ggplot2::labs(title = glue::glue("n_WellFlag: {min} <= n <= {max}"),
                      y = "Well object count (n)")
      
      # return
      message(glue::glue("A diagnostic plot for checking the object number thresholds (max, min) is returned. See <out>$.p"))
      out <- list(d = d, p = p)
      return(out)
    } else {
      # return data only
      return(d)
    }
  }
  
  message("--Running the nWF function -> n_Cb_ns --")
  n_Cb_ns <- nWF_custom(data = titer_Cb_data_ns, drug, concentration_um,
                        max = 75,
                        min =5,
                        plot = TRUE)
  
  n_Cb_plot_ns <- n_Cb_ns$p
  
  n_Cb_data_ns <- n_Cb_ns$d
  
# summarizeWells + add width
summarizeWells_custom <- function(data, length_col = "Length_um", width_col = "Width_um", OF = "filter", drop = TRUE) {
  
  if(!(OF %in% c("ignore", "filter"))) {
    stop('Invalid OF argument. Please set the OF argument to either "ignore" or "filter"')
  }
  
  if(OF == "filter") {
    # Filter the flagged objects (if you have any object-level flags)
    # Check if size_flag column exists
    if("size_flag" %in% names(data)) {
      data.of <- data %>%
        dplyr::filter(size_flag == "pass" | is.na(size_flag))
      message("All flagged objects are filtered prior to summarizing wells.")
    } else {
      data.of <- data
      message("No object flags detected. Proceeding with all objects.")
    }
    
    if(drop == TRUE) {
      # Drop object-level data (keep only metadata and measurements)
      to.summarize <- data.of %>%
        dplyr::select(-contains("flag"),
                      -contains("_flag"),
                      -contains("ObjectFlag"))
      message("Object-level flag variables are dropped from the summarized data.")
    } else {
      to.summarize <- data.of
    }
  }
  
  if(OF == "ignore") {
    if(drop == TRUE) {
      # Drop object data
      to.summarize <- data %>%
        dplyr::select(-contains("flag"),
                      -contains("_flag"),
                      -contains("ObjectFlag"))
      message("Object-level flag variables are dropped from the summarized data.")
    } else {
      to.summarize <- data
    }
  }
  
  # Do the summary for both length and width
  out <- to.summarize %>%
    dplyr::group_by(Metadata_Plate, Metadata_Well) %>%
    dplyr::mutate(
      # Length summaries
      mean_wormlength_um = mean(.data[[length_col]], na.rm = TRUE),
      min_wormlength_um = as.numeric(quantile(.data[[length_col]], na.rm = TRUE)[1]),
      q10_wormlength_um = as.numeric(quantile(.data[[length_col]], probs = 0.1, na.rm = TRUE)[1]),
      q25_wormlength_um = as.numeric(quantile(.data[[length_col]], probs = 0.25, na.rm = TRUE)[1]),
      median_wormlength_um = median(.data[[length_col]], na.rm = TRUE),
      sd_wormlength_um = sd(.data[[length_col]], na.rm = TRUE),
      q75_wormlength_um = as.numeric(quantile(.data[[length_col]], probs = 0.75, na.rm = TRUE)[1]),
      q90_wormlength_um = as.numeric(quantile(.data[[length_col]], probs = 0.90, na.rm = TRUE)[1]),
      max_wormlength_um = as.numeric(quantile(.data[[length_col]], na.rm = TRUE)[5]),
      cv_wormlength_um = (sd_wormlength_um / mean_wormlength_um),
      
      # Width summaries
      mean_wormwidth_um = mean(.data[[width_col]], na.rm = TRUE),
      min_wormwidth_um = as.numeric(quantile(.data[[width_col]], na.rm = TRUE)[1]),
      q10_wormwidth_um = as.numeric(quantile(.data[[width_col]], probs = 0.1, na.rm = TRUE)[1]),
      q25_wormwidth_um = as.numeric(quantile(.data[[width_col]], probs = 0.25, na.rm = TRUE)[1]),
      median_wormwidth_um = median(.data[[width_col]], na.rm = TRUE),
      sd_wormwidth_um = sd(.data[[width_col]], na.rm = TRUE),
      q75_wormwidth_um = as.numeric(quantile(.data[[width_col]], probs = 0.75, na.rm = TRUE)[1]),
      q90_wormwidth_um = as.numeric(quantile(.data[[width_col]], probs = 0.90, na.rm = TRUE)[1]),
      max_wormwidth_um = as.numeric(quantile(.data[[width_col]], na.rm = TRUE)[5]),
      cv_wormwidth_um = (sd_wormwidth_um / mean_wormwidth_um),
      
      # Count
      n = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Metadata_Plate, Metadata_Well, .keep_all = TRUE)
  
  # Return it
  return(out)
}

message ("-- ran summarizeWells -> raw_wells_Cb_ns --")

raw_wells_Cb_ns <- summarizeWells_custom(
  data = merged_Cb_ns,
  length_col = "Length_um",
  # width_col = "Width_um",
  OF = "filter",
  drop = TRUE
)

message("-- merge the nWF results with the summarized wells results -> merged_summarize_Cb_ns")

# need to merge n_Cb_ns with the original df
merged_summarize_Cb_ns <- raw_wells_Cb_ns %>%
  dplyr::left_join(n_Cb_ns$d, by = c("Metadata_Plate",
                                     "Metadata_Well",
                                     "diluent",
                                     "concentration_um",
                                     "drug",
                                     "strain",
                                     "bleach",
                                     "n"))


# outlierWF
iqrOutlier <- function(x, thresh = 1.5) { # Helper function: IQR outlier detection
  q25 <- quantile(x, 0.25, na.rm = TRUE)
  q75 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q75 - q25
  lower <- q25 - thresh * iqr
  upper <- q75 + thresh * iqr
  
  outliers <- x < lower | x > upper
  return(outliers)
}

sdOutlier <- function(x, thresh = 3) { # Helper function: SD outlier detection
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  lower <- mean_x - thresh * sd_x
  upper <- mean_x + thresh * sd_x
  
  outliers <- x < lower | x > upper
  return(outliers)
}

# outlierWF function - modified for multiple phenotypes
outlierWF_custom <- function(data, ..., phenotypes = c("median_wormlength_um", "median_wormwidth_um"),
                             iqr = TRUE, thresh = NULL, filterWF = TRUE) {
  grouping_vars <- enquos(...)
  
  # Set the threshold
  if(iqr == TRUE & is.null(thresh)) {
    thresh <- 1.5
  }
  if(iqr == FALSE & is.null(thresh)) {
    thresh <- 3
  }
  
  # Decide how to handle flagged wells and warn
  if(filterWF == TRUE) {
    message("Previously flagged wells will not be used when calculating outliers within the group. This is the recommended approach.")
    # Filter out previously flagged wells
    data.in <- data %>%
      dplyr::filter(is.na(titer_WellFlag) & is.na(n_WellFlag))
  } else {
    message("WARNING: Previously flagged wells will be used when calculating outliers within the group. This is NOT the recommended approach. Please consider whether filterWF = TRUE is better.")
    data.in <- data
  }
  
  # Flag outliers for each phenotype
  out <- data.in
  
  for(phenotype in phenotypes) {
    flag_col_name <- paste0("outlier_", gsub("median_worm|_um", "", phenotype), "_WellFlag")
    
    if(iqr == TRUE) {
      # Use IQR method
      message(glue::glue("Flagging outlier wells in group if {phenotype} is outside the range: median +/- ({thresh}*IQR)"))
      out <- out %>%
        dplyr::group_by(!!!grouping_vars) %>%
        dplyr::mutate(
          !!flag_col_name := iqrOutlier(.data[[phenotype]], thresh = thresh),
          !!flag_col_name := ifelse(.data[[flag_col_name]] == TRUE, "outlier", NA_character_)
        ) %>%
        dplyr::ungroup()
    }
    
    if(iqr == FALSE) {
      # Use SD method
      message(glue::glue("Flagging outlier wells in group with {phenotype} {thresh} SDs away from mean"))
      out <- out %>%
        dplyr::group_by(!!!grouping_vars) %>%
        dplyr::mutate(
          !!flag_col_name := sdOutlier(.data[[phenotype]], thresh = thresh),
          !!flag_col_name := ifelse(.data[[flag_col_name]] == TRUE, "outlier", NA_character_)
        ) %>%
        dplyr::ungroup()
    }
  }
  
  # Create combined outlier flag (if ANY phenotype is an outlier)
  outlier_cols <- paste0("outlier_", gsub("median_worm|_um", "", phenotypes), "_WellFlag")
  
  out <- out %>%
    dplyr::mutate(
      outlier_WellFlag = dplyr::case_when(
        dplyr::if_any(dplyr::all_of(outlier_cols), ~ . == "outlier") ~ "outlier",
        TRUE ~ NA_character_
      )
    )
  
  # Add back any filtered data
  suppressMessages(out <- data %>%
                     dplyr::full_join(out, by = names(data)[names(data) %in% names(out)]))
  
  # Return it
  return(out)
}

message("--running outlierWF custom -> ow_Cb_ns --")

# Flag outliers based on both length and width
ow_Cb_ns <- outlierWF_custom(
  data = merged_summarize_Cb_ns,
  bleach, drug, concentration_um, strain,
  phenotypes = c("median_wormlength_um", "median_wormwidth_um"),
  iqr = TRUE,
  thresh = 1.5,
  filterWF = TRUE
) %>%
  dplyr::mutate(assay_bleach = paste(bleach, sep = "_"))

# View the flags
ow_Cb_ns %>%
  dplyr::select(Metadata_Plate, Metadata_Well, strain, drug, concentration_um,
                median_wormlength_um, median_wormwidth_um,
                outlier_length_WellFlag, outlier_width_WellFlag, outlier_WellFlag) %>%
  dplyr::filter(!is.na(outlier_WellFlag))


# checkWF modified
checkWF_custom <- function(data, ..., plot = TRUE) {
  # Find the WellFlags in data in the order they appear
  uf <- names(data %>% dplyr::select(contains("_WellFlag")))
  uf.short <- unlist(lapply(data[uf], function(x) unique(na.omit(x))))
  names(uf.short) <- NULL
  
  # Remove duplicates while preserving order
  uf.short <- unique(uf.short)
  
  # send an error if needed
  if(length(uf) == 0 ) {
    stop("No WellFlags detected, did you specify the correct data? See checkWF_custom() help for details.")
  }
  
  # tell us about what was found
  message(glue::glue("{length(uf)} WellFlags detected in data. They were applied in the following order:"))
  for(i in 1:length(uf)) {
    message(glue::glue("{uf[i]}"))
  }
  
  # Get a single WellFlag vector based on the order in which the flags were run
  # This creates a combined flag showing the FIRST flag that was triggered
  wellFlag <- data %>%
    dplyr::select(contains("_WellFlag")) %>%
    tidyr::unite(wellFlag, sep = "_", na.rm = TRUE) %>%
    tidyr::separate(wellFlag, into = "wellFlag", sep = "_", extra = "drop", remove = TRUE) %>%
    dplyr::pull(wellFlag)
  
  # summarized by
  grouping_vars <- enquos(...)
  sum_by <- paste(names(data %>% dplyr::select(!!!grouping_vars)), collapse = ", ")
  message(glue::glue("The data are summarized by: {sum_by}"))
  
  # Add it to data and summarize by grouping variables if provided.
  summary <- data %>%
    dplyr::bind_cols(wellFlag = wellFlag) %>%
    dplyr::mutate(grand_n = dplyr::n(),
                  grouping = sum_by) %>%
    dplyr::group_by(!!!grouping_vars) %>%
    dplyr::mutate(group_n = dplyr::n()) %>%
    dplyr::group_by(!!!grouping_vars, wellFlag) %>%
    dplyr::mutate(wellFlag_n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(wellFlag = ifelse(wellFlag == "" | is.na(wellFlag), "noFlag", wellFlag),
                  wellFlag_group_perc = wellFlag_n / group_n) %>%
    # set levels from user flag order (now with unique values)
    dplyr::mutate(wellFlag = factor(wellFlag, levels = c(uf.short, "noFlag"))) %>%
    dplyr::arrange(!!!grouping_vars, wellFlag) %>%
    dplyr::distinct(!!!grouping_vars, grouping, wellFlag, wellFlag_group_perc, grand_n, group_n, wellFlag_n) %>%
    dplyr::select(!!!grouping_vars, grouping, wellFlag, wellFlag_group_perc, grand_n, group_n, wellFlag_n)
  
  
  if(plot == TRUE) {
    # make a plot
    p <- ggplot2::ggplot(summary %>% dplyr::group_by(!!!grouping_vars)) +
      ggplot2::aes(x = "group", y = wellFlag_group_perc, fill = wellFlag, label = wellFlag_n) +
      ggplot2::geom_col() +
      ggplot2::geom_text(size = 3, position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::facet_wrap(vars(!!!grouping_vars)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "", y = "fraction") +
      ggplot2::theme(strip.background = ggplot2::element_rect(
        color="black", fill="white", size=0.5, linetype="solid"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank())
    # return data and plot
    message("Returning list with elements d (the summary data frame) and p (the summary plot)")
    out <- list(d = summary, p = p)
    return(out)
  } else {
    message("No summary plots made. Set plot = TRUE to make plots")
    # return data only
    return(summary)
  }
}

# Use the checkWF_custom function
message("Running checkWF() --> cw_Cb_ns")
cw_Cb_ns <- checkWF_custom(data = ow_Cb_ns, drug, strain, concentration_um, plot = TRUE)

# Look at the plot and ask if user wants to save
cw_Cb_ns$p
cw_Cb_plot_ns <- cw_Cb_ns$p

# Look at the summary data
cw_Cb_ns$d

# View detailed breakdown
cw_Cb_ns$d %>%
  dplyr::arrange(strain, drug, concentration_um, wellFlag)

# Want to see outliers by length or width or both
summarize_outlier_types <- function(data, ...) {
  grouping_vars <- enquos(...)
  
  data %>%
    dplyr::group_by(!!!grouping_vars) %>%
    dplyr::summarise(
      total_wells = n(),
      
      # Individual flags
      titer_flagged = sum(!is.na(titer_WellFlag)),
      n_flagged = sum(!is.na(n_WellFlag)),
      
      # Outlier breakdowns
      length_only = sum(!is.na(outlier_length_WellFlag) & is.na(outlier_width_WellFlag)),
      width_only = sum(is.na(outlier_length_WellFlag) & !is.na(outlier_width_WellFlag)),
      both_length_width = sum(!is.na(outlier_length_WellFlag) & !is.na(outlier_width_WellFlag)),
      any_outlier = sum(!is.na(outlier_WellFlag)),
      
      # No flags
      no_flags = sum(is.na(titer_WellFlag) & is.na(n_WellFlag) & is.na(outlier_WellFlag)),
      
      # Percentages
      pct_length_only = round(100 * length_only / total_wells, 1),
      pct_width_only = round(100 * width_only / total_wells, 1),
      pct_both = round(100 * both_length_width / total_wells, 1),
      pct_any_outlier = round(100 * any_outlier / total_wells, 1),
      pct_no_flags = round(100 * no_flags / total_wells, 1),
      
      .groups = "drop"
    )
}

message("--running function summarize_outlier_types -> outlier_summary--")
# Create the summary
outlier_summary <- summarize_outlier_types(ow_Cb_ns, drug, strain, concentration_um)

outlier_summary %>% # view key metrics
  dplyr::select(drug, strain, concentration_um, total_wells,
                length_only, width_only, both_length_width, any_outlier, no_flags)

outlier_summary %>%
  dplyr::select(drug, strain, concentration_um, total_wells,
                length_only, width_only, both_length_width,
                pct_length_only, pct_width_only, pct_both) %>%
  kable(caption = "Outlier Wells by Type")


outlier_plot_data <- outlier_summary %>% # Reshape for plotting
  dplyr::select(drug, strain, concentration_um, length_only, width_only, both_length_width, no_flags) %>%
  tidyr::pivot_longer(cols = c(length_only, width_only, both_length_width, no_flags),
                      names_to = "flag_type",
                      values_to = "count") %>%
  dplyr::mutate(flag_type = factor(flag_type,
                                   levels = c("no_flags", "length_only", "width_only", "both_length_width"),
                                   labels = c("No Flags", "Length Only", "Width Only", "Both")))

# Plot to look at outliers
outlier_graph <- ggplot2::ggplot(outlier_plot_data) +
  ggplot2::aes(x = factor(concentration_um), y = count, fill = flag_type) +
  ggplot2::geom_col(position = "stack", width = 0.8) +
  ggplot2::facet_wrap(~drug + strain, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = c("No Flags" = "gray80",
                                        "Length Only" = "#E69F00",
                                        "Width Only" = "#56B4E9",
                                        "Both" = "#D55E00")) +
  ggplot2::labs(title = "Wells Flagged as Outliers by Type",
                x = "Concentration (µM)",
                y = "Number of Wells",
                fill = "Outlier Type") +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

outlier_graph

# ok for now I only want to remove outliers for length not width so the function is modified here:
filterWF_custom_length <- function(data, rmVars = TRUE) {
  # Find the WellFlags in data in the order they appear
  uf <- names(data %>% dplyr::select(contains("_WellFlag")))
  
  # Check if any flags exist
  if(length(uf) == 0) {
    stop("No WellFlags detected in data. Cannot filter.")
  }
  
  message(glue::glue("{length(uf)} WellFlags detected in data. They were applied in the following order:"))
  for(i in 1:length(uf)) {
    message(glue::glue("{uf[i]}"))
  }
  
  # Check if length outlier flag exists
  if(!"outlier_length_WellFlag" %in% names(data)) {
    stop("outlier_length_WellFlag not found in data. Cannot filter by length outliers.")
  }
  
  # Filter wells that have length outlier flags
  # This includes wells with ONLY length outliers OR both length and width outliers
  data.filterWF <- data %>%
    dplyr::filter(is.na(outlier_length_WellFlag))
  
  # Report filtering stats
  n_original <- nrow(data)
  n_filtered <- nrow(data.filterWF)
  n_removed <- n_original - n_filtered
  pct_removed <- round(100 * n_removed / n_original, 1)
  
  # Report what types were filtered
  length_only <- data %>%
    filter(!is.na(outlier_length_WellFlag) & is.na(outlier_width_WellFlag)) %>%
    nrow()
  
  both_outliers <- data %>%
    filter(!is.na(outlier_length_WellFlag) & !is.na(outlier_width_WellFlag)) %>%
    nrow()
  
  message(glue::glue("Filtered {n_removed} of {n_original} wells ({pct_removed}%)."))
  message(glue::glue("  - {length_only} wells with length outliers only"))
  message(glue::glue("  - {both_outliers} wells with both length and width outliers"))
  message(glue::glue("{n_filtered} wells remain."))
  
  # Catch full filter
  if(nrow(data.filterWF) == 0){
    stop("All data filtered due to flags. Please review any WF function settings to retain more data.")
  }
  
  if(rmVars == FALSE) {
    # return with flags
    return(data.filterWF)
  }
  if(rmVars == TRUE) {
    data.filterWF.rm <- data.filterWF %>%
      dplyr::select(-dplyr::contains("_WellFlag"),
                    -dplyr::contains("wellFlag"))
    
    # return without flags
    return(data.filterWF.rm)
  }
}

# Use the modified function

#changed rmVars = FALSE because we only filtered Length and both but not just width
#and there's a second filter that needs to run? 

message("ran function filteredWF_custom_length -> fw_Cb_ns_length")

fw_Cb_ns_length <- filterWF_custom_length(data = ow_Cb_ns, rmVars = FALSE)

filterWF_custom <- function(data, rmVars = TRUE) {
  # Find the WellFlags in data in the order they appear
  uf <- names(data %>% dplyr::select(contains("_WellFlag")))
  
  # Check if any flags exist
  if(length(uf) == 0) {
    stop("No WellFlags detected in data. Cannot filter.")
  }
  
  message(glue::glue("{length(uf)} WellFlags detected in data. They were applied in the following order:"))
  for(i in 1:length(uf)) {
    message(glue::glue("{uf[i]}"))
  }
  
  # Get a single WellFlag vector based on the order in which the flags were run
  wellFlag <- data %>%
    dplyr::select(contains("_WellFlag")) %>%
    tidyr::unite(wellFlag, sep = "_", na.rm = TRUE) %>%
    tidyr::separate(wellFlag, into = "wellFlag", sep = "_", extra = "drop", remove = TRUE) %>%
    dplyr::pull(wellFlag)
  
  if("wellFlag" %in% names(data)) {
    message("WARNING: The wellFlag variable already exists in the data. It will be overwritten in the output.")
    # Add it to data and summarize by grouping variables if provided.
    data.WF <- data %>%
      dplyr::select(-wellFlag) %>%
      dplyr::bind_cols(wellFlag = wellFlag)
    
  } else {
    # Add it to data.
    data.WF <- data %>%
      dplyr::bind_cols(wellFlag = wellFlag)
  }
  
  # Filter the flagged wells
  data.filterWF <- data.WF %>%
    dplyr::filter(wellFlag == "" | is.na(wellFlag))
  
  # Report filtering stats
  n_original <- nrow(data)
  n_filtered <- nrow(data.filterWF)
  n_removed <- n_original - n_filtered
  pct_removed <- round(100 * n_removed / n_original, 1)
  
  message(glue::glue("Filtered {n_removed} of {n_original} wells ({pct_removed}%). {n_filtered} wells remain."))
  
  # Catch full filter
  if(nrow(data.filterWF) == 0){
    stop("All data filtered due to flags. Please review any WF function settings to retain more data.")
  }
  
  if(rmVars == FALSE) {
    # return with flags
    return(data.filterWF)
  }
  if(rmVars == TRUE) {
    data.filterWF.rm <- data.filterWF %>%
      dplyr::select(-dplyr::contains("_WellFlag"),
                    -dplyr::contains("wellFlag"))
    
    # return without flags
    return(data.filterWF.rm)
  }
}

# Use the filterWF_custom function and drop the flagging variables afterward

#code get's snagged here works after I changed prev filters rmVars = FALSE****
message("ran filteredWF -> fw_Cb_ns")
fw_Cb_ns <- filterWF_custom(data = fw_Cb_ns_length, rmVars = TRUE)

# 5 WellFlags detected in data. They were applied in the following order:
# titer_WellFlag
# n_WellFlag
# outlier_length_WellFlag
# outlier_width_WellFlag
# outlier_WellFlag
# Filtered 240 of 1079 wells (22.2%). 839 wells remain.

# Check results
nrow(ow_Cb_ns)  # before filtering
nrow(fw_Cb_ns)  # after filtering (length and width)
nrow(fw_Cb_ns_length)  #after filtering (just length)

# Check which strains/conditions remain
fw_Cb_ns_length %>%
  count(drug, strain, concentration_um)

# check balance custom
checkBalance_custom <- function(data, ..., design, x, size = 3) {
  grouping_vars <- enquos(...)
  x_var <- enquo(x)
  
  # Get the vars
  vars <- names(data %>% dplyr::select(!!!grouping_vars, !!x_var))
  
  # Warn if design has NAs for strain or drug
  if(anyNA(design$strain)) {
    stop(glue::glue("ERROR: There are {sum(is.na(design$strain))} rows with NAs for strain in the design file. Please review the design file."))
  }
  if(anyNA(design$drug)) {
    stop(glue::glue("ERROR: There are {sum(is.na(design$drug))} rows with NAs for drug in the design file. Please review the design file."))
  }
  
  # Check if required columns exist
  if(!all(c("Metadata_Plate", "Metadata_Well") %in% names(data))) {
    stop("ERROR: Metadata_Plate and Metadata_Well columns must be in the data.")
  }
  
  if(FALSE %in% (vars %in% names(design))) {
    stop(paste("ERROR: One or more of the specified variables <", paste(vars, collapse = ", "), "> are not found in the design. Please add them to use this function."))
  }
  
  # Create well.id in data if it doesn't exist
  if(!("well.id" %in% names(data))) {
    data <- data %>%
      dplyr::mutate(well.id = paste(Metadata_Plate, Metadata_Well, sep = "_"))
  }
  
  # Get the design, select what's needed and add a well id
  d.wells <- design %>%
    dplyr::mutate(well.id = paste(Metadata_Plate, Metadata_Well, sep = "_")) %>%
    dplyr::select(well.id, !!x_var, !!!grouping_vars)
  
  # Get lost wells
  d.wells$lost <- !(d.wells$well.id %in% data$well.id)
  
  # Get counts using the design
  summary <- d.wells %>%
    dplyr::left_join(data, by = c("well.id", vars)) %>%
    dplyr::group_by(!!!grouping_vars, !!x_var) %>%
    dplyr::mutate(!!quo_name(x_var) := as.character(!!x_var)) %>%
    dplyr::mutate(g.total = dplyr::n(),
                  n.lost = sum(lost == TRUE),
                  n.kept = g.total - n.lost,
                  perc.lost = n.lost/g.total,
                  perc.kept = n.kept/g.total) %>%
    dplyr::distinct(!!!grouping_vars, !!x_var, g.total, n.lost, n.kept, perc.lost, perc.kept) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = n.lost:n.kept, names_to = "n.type", values_to = "n") %>%
    tidyr::pivot_longer(cols = perc.lost:perc.kept, names_to = "perc.type", values_to = "perc")  %>%
    dplyr::mutate(type = dplyr::case_when(n.type == "n.lost" & perc.type == "perc.lost" ~ "lost",
                                          n.type == "n.kept" & perc.type == "perc.kept" ~ "kept",
                                          TRUE ~ NA_character_)) %>%
    dplyr::filter(!is.na(type)) %>%
    dplyr::select(-n.type, -perc.type) %>%
    dplyr::mutate(type = factor(type, levels = c("lost", "kept")))
  
  p <- ggplot2::ggplot(summary) +
    ggplot2::aes(x = !!x_var, y = perc, fill = type, label = n) +
    ggplot2::geom_col() +
    ggplot2::geom_text(size = size, position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::geom_text(ggplot2::aes(x = !!x_var, y = 1.05, label = g.total), size = size) +
    ggplot2::ylim(0, 1.05) +
    ggplot2::facet_wrap(vars(!!!grouping_vars)) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "fraction") +
    ggplot2::theme(strip.background = ggplot2::element_rect(
      color="black", fill="white", size=0.5, linetype="solid"))
  
  # Return data and plot
  message("Returning list with elements d (the summary data frame) and p (the summary plot)")
  out <- list(d = summary, p = p)
  return(out)
}

# Use the checkBalance_custom function
# Add assay_bleach variable to design if needed
design_with_bleach <- design_Cb_ns %>%
  dplyr::mutate(assay_bleach = paste(bleach, sep = "_"))

message("Running checkBalance_custom -> cb_CB_ns")

cb_Cb_ns <- checkBalance_custom(
  data = fw_Cb_ns_length,
  drug, concentration_um,
  design = design_with_bleach,
  x = strain,  # or assay_bleach
  size = 3
)

# Look at the plot and add a nicer x-axis
cb_Cb_plot_ns <- cb_Cb_ns$p +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggplot2::geom_hline(yintercept = 0.75, linetype = 2, color = "red")

cb_Cb_plot_ns
# Look at the summary data
cb_Cb_ns$d

# Check overall balance
cb_Cb_ns$d %>%
  group_by(drug, concentration_um) %>%
  summarise(
    total_wells = mean(g.total),
    total_kept = sum(n[type == "kept"]),
    total_lost = sum(n[type == "lost"]),
    avg_perc_kept = mean(perc[type == "kept"]),
    .groups = "drop"
  )

checkEff_custom <- function(data, ..., x, y, fill = NULL, size = 1.5, scales = "fixed") {
  
  grouping_vars <- enquos(...)
  x_var <- enquo(x)
  y_var <- enquo(y)
  fill_var <- enquo(fill)
  
  # Check for non-continuous x, group if necessary
  x_values <- data %>% dplyr::pull(!!x_var)
  
  if(is.character(x_values) | is.factor(x_values) | is.logical(x_values)) {
    p <- ggplot2::ggplot(data) +
      ggplot2::aes(x = !!x_var, y = !!y_var) +
      ggplot2::facet_wrap(vars(!!!grouping_vars), scales = scales) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(
        color="black", fill="white", size=0.5, linetype="solid"),
        panel.grid = ggplot2::element_blank())
  } else {
    p <- ggplot2::ggplot(data) +
      ggplot2::aes(x = !!x_var, y = !!y_var, group = !!x_var) +
      ggplot2::facet_wrap(vars(!!!grouping_vars), scales = scales) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(
        color="black", fill="white", size=0.5, linetype="solid"),
        panel.grid = ggplot2::element_blank())
  }
  
  # Check on optional fill parameter
  if(rlang::quo_is_null(fill_var)){
    p.out <- p +
      ggplot2::geom_jitter(width = 0.25,
                           size = size,
                           fill = "grey",
                           shape = 21) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.35, fill = "white")
  } else {
    p.out <- p +
      ggplot2::geom_jitter(ggplot2::aes(fill = as.character(!!fill_var)),
                           position = ggplot2::position_jitter(width = 0.25),
                           size = size,
                           shape = 21) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.35, fill = "white") +
      ggplot2::labs(fill = rlang::as_label(fill_var))
  }
  
  # Return it
  return(p.out)
}

# Use the checkEff_custom function
message("running checkEff() -> checkeff_Cb_ns")
checkeff_Cb_ns <- checkEff_custom(
  data = fw_Cb_ns_length,
  drug, strain,
  x = concentration_um,
  y = median_wormlength_um,
  fill = bleach,
  size = 1.5,
  scales = "free_x"
)

# Look at the plot
checkeff_Cb_ns

axis_breaks <- function(drug_name){
  if (drug_name == "ABZ"){
    return(c(0.1170, 0.2340, 0.4680, 0.9375, 1.8750, 3.7500, 7.5000, 15.0000, 30.0000, 60.0000, 120.0000))
  }
  if (drug_name == "CLO"){
    return(c(0.000, 0.097, 0.195, 0.390, 0.780, 1.560, 3.125, 6.250, 12.500, 25.000, 50.000, 100))
  }
  if (drug_name == "DEC"){
    return(c(0.000, 0.645, 1.290, 2.580, 5.160, 10.300, 20.625, 41.250, 82.500, 165.000, 330.000, 660.000))
  }
  if (drug_name == "EMO"){
    return(c(0.00000, 0.00097, 0.00195, 0.00390, 0.00780, 0.01560, 0.03120, 0.06250, 0.12500, 0.25000, 0.50000, 1.00000))
  }
  if (drug_name == "FBZ"){
    return(c(0.000, 0.117, 0.234, 0.468, 0.937, 1.875, 3.750, 7.500, 15.000, 30.000, 60.000, 120.000))
  }
  if (drug_name == "IVM"){
    return(c(0, 0.000029, 0.000058, 0.00011, 0.00023, 0.00046, 0.00093, 0.0018, 0.0037, 0.0075, 0.015, 0.030))
  }
  if (drug_name == "LEV"){
    return(c(0.0000, 0.1360, 0.2730, 0.5460, 1.0900, 2.1875, 4.3750, 8.7500, 17.5000, 35.0000, 70.0000, 140.0000))
  }
  if (drug_name == "MOX"){
    return(c(0.00000, 0.00097, 0.00190, 0.00390, 0.00780, 0.01560, 0.03125,
             0.06250, 0.12500, 0.25000, 0.50000, 1.00000))
    }
  if (drug_name == "PYR"){
    return(c(0.000, 0.292, 0.585, 1.170, 2.340, 4.680, 9.375 , 18.750, 37.500, 75.000, 150.000, 300.000))
  }
  if (drug_name == "MILB"){
    return(c(0.00000, 0.00097, 0.00190, 0.00390, 0.00780, 0.01560, 0.03125,
             0.06250, 0.12500, 0.25000, 0.50000, 1.00000))
  }
}

axis_labels <- function(drug_name){
  if (drug_name == "ABZ"){
    return(c("0.117", "0.234", "0.468", "0.9375", "1.875", "3.75", "7.5", "15", "30", "60", "120"))
  }
  if (drug_name == "CLO"){
    return(c("0", "0.097", "0.195", "0.39", "0.78", "1.56", "3.125", "6.25"
             , "12.5", "25", "50", "100"))
  }
  if(drug_name == "DEC"){
    return(c("0", "0.645", "1.29", "2.58", "5.16", "10.3", 
             "20.625", "41.25", "82.5", "165", "330", "660"))
  }
  if (drug_name == "EMO"){
    return(c("0", "0.00097", "0.00195", "0.0039", "0.0078", "0.0156", "0.0312", "0.0625", "0.125", "0.25", "0.5", "1.0"))
  }
  if (drug_name == "FBZ"){
    return(c("0", "0.117", "0.234", "0.468", "0.937", "1.875", "3.75", "7.5", "15", "30", "60", "120"))
  }
  if (drug_name == "IVM"){
    return(c("0", "0.000029", "0.000058", "0.00011", "0.00023", "0.00046", "0.00093", "0.0018", "0.0037", "0.0075", "0.015", "0.03"))
  }
  if (drug_name == "LEV"){
    return(c( "0", "0.136", "0.273", "0.546", "1.09", "2.1875", "4.375", "8.75", "17.5", "35" , "70" , "140"))
  }
  if (drug_name == "MOX"){
    return(c("0", "0.00097", "0.0019", "0.0039", "0.0078", "0.0156", "0.03125",
             "0.0625", "0.125", "0.25", "0.5", "1"))
    }
  if (drug_name == "PYR"){
    return(c("0", "0.292", "0.585", "1.170", "2.34", "4.68", "9.375" , "18.75", "37.5" , "75" , "150" , "300"))
  }
  if (drug_name == "MILB"){
    return(c("0", "0.00097", "0.0019", "0.0039", "0.0078", "0.0156", "0.03125",
               "0.0625", "0.125", "0.25", "0.5", "1"))
  }
}

fw_Cb_ns_length %>%
  distinct(concentration_um) %>%
  arrange(concentration_um)

checkeff_Cb_ns <- checkEff_custom( # need to force breaks and change scale because of the small drug concentrations
  data = fw_Cb_ns_length %>% filter(concentration_um > 0),  # Remove 0 for log scale
  drug, strain,
  x = concentration_um,
  y = median_wormlength_um,
  fill = bleach,
  size = 1.5,
  scales = "free_x"
) +
  scale_x_log10(
    breaks = axis_breaks(b), #note this scale will need to be different for each drug - so look at drug range
    labels = axis_labels(b)
  ) +
  annotation_logticks(sides = "b")

checkeff_Cb_ns

regEff_custom <- function(data, ..., d.var, c.var) {
  grouping_vars <- enquos(...)
  d_var <- enquo(d.var)
  c_var <- enquo(c.var)
  
  # Get variable names
  d_var_name <- as_label(d_var)
  c_var_name <- as_label(c_var)
  
  # Get data classes
  d.c <- class(data %>% dplyr::pull(!!d_var))
  c.c <- class(data %>% dplyr::pull(!!c_var))
  
  # Report data classes for dependent and confounding variable
  message(glue::glue('The dependent variable `{d_var_name}` is class {d.c}. Please ensure this is correct.'))
  message(glue::glue('The confounding variable `{c_var_name}` is class {c.c}. Please ensure this is correct.'))
  
  # Regression groups
  group_by <- paste(names(data %>% dplyr::select(!!!grouping_vars)), collapse = ", ")
  message(glue::glue("The data are grouped by: `{group_by}`"))
  
  # Filter if NAs present in c and d
  f <- data %>%
    dplyr::filter(!is.na(!!d_var) & !is.na(!!c_var))
  
  # Report on filters if present
  if(nrow(f) != nrow(data)) {
    message(glue::glue("There were {nrow(data %>% dplyr::filter(is.na(!!d_var)))} rows filtered with NAs for {d_var_name}"))
    message(glue::glue("There were {nrow(data %>% dplyr::filter(is.na(!!c_var)))} rows filtered with NAs for {c_var_name}"))
  }
  
  # Apply the linear model to nested data
  m <- f %>%
    dplyr::group_by(!!!grouping_vars) %>%
    tidyr::nest() %>%
    dplyr::mutate(model = purrr::map(data, .f = ~ stats::lm(
      formula = reformulate(c_var_name, response = d_var_name, intercept = FALSE),
      data = .x
    ))) %>%
    dplyr::ungroup()
  
  # Get the model summaries and residuals to add to the data
  m.sum <- list()
  
  for(i in 1:nrow(m)) {
    # Get the model outputs
    coeffs <- tibble::as_tibble(data.table::data.table(stats::coef(summary(m[[i, "model"]][[1]])), keep.rownames = 'term')) %>%
      dplyr::mutate(term = stringr::str_replace(term, pattern = glue::glue("{c_var_name}"), replacement = "")) %>%
      dplyr::rename(
        !!c_var_name := term,
        !!paste0(c_var_name, "_coeff") := Estimate,
        !!paste0(c_var_name, "_se") := `Std. Error`,
        !!paste0(c_var_name, "_tvalue") := `t value`,
        !!paste0(c_var_name, "_pvalue") := `Pr(>|t|)`
      ) %>%
      dplyr::select(!!c_var_name, dplyr::everything())
    
    # Get the grouping vars
    group.vars.c <- m %>%
      dplyr::select(-data, -model) %>%
      dplyr::slice(rep(i, each = nrow(coeffs)))
    
    # Add these
    coeffs.out <- cbind(group.vars.c, coeffs)
    
    # Add to the list
    m.sum[[i]] <- coeffs.out
  }
  
  # Bind these data frames by row
  coeffs.d <- dplyr::bind_rows(m.sum)
  
  # Join model summary and recalculate the residuals
  suppressMessages(
    d <- f %>%
      dplyr::group_by(!!!grouping_vars) %>%
      # Recalculate residuals
      dplyr::mutate(
        !!paste0(d_var_name, "_reg") := stats::residuals(
          stats::lm(reformulate(c_var_name, response = d_var_name, intercept = FALSE), data = cur_data())
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(., coeffs.d)
  )
  
  # Plot the effect in context
  p1 <- checkEff_custom(data = f, !!!grouping_vars, x = !!c_var, y = !!d_var) +
    ggplot2::labs(subtitle = glue::glue("{c_var_name} effect on {d_var_name}")) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(color="black", fill="white", size=0.5, linetype="solid"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
    )
  
  # Plot the coefficients and their standard error from the models
  coeffs.p <- coeffs.d %>%
    dplyr::mutate(
      coeff = .data[[paste0(c_var_name, "_coeff")]],
      se = .data[[paste0(c_var_name, "_se")]],
      min = coeff - se,
      max = coeff + se
    )
  
  # Plot regression coefficients
  p2 <- ggplot2::ggplot(coeffs.p) +
    ggplot2::aes(x = !!c_var, y = coeff, ymin = min, ymax = max) +
    ggplot2::geom_pointrange(size = 0.5, shape = 21, stroke = 0.25) +
    ggplot2::facet_wrap(vars(!!!grouping_vars)) +
    ggplot2::labs(
      title = "Regression coefficients +/- s.e.",
      subtitle = glue::glue("lm({d_var_name} ~ {c_var_name} - 1)"),
      y = glue::glue("{c_var_name} coefficient")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(color="black", fill="white", size=0.5, linetype="solid"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)
    )
  
  # Return
  return(list(d = d, p1 = p1, p2 = p2, model = m))
}

delta_custom <- function(data, ..., WF = "filter", vars = c("median_wormlength_um", "cv_wormlength_um"), doseR = FALSE) {
  
  grouping_vars <- enquos(...)
  
  # Expecting message
  if(doseR == TRUE){
    message("You set doseR = TRUE. Expecting controls to be coded as for a dose response.")
  }
  if(doseR == FALSE){
    message("You set doseR = FALSE. Not expecting controls to be coded for a dose response.")
  }
  
  # Check on control coding
  controls <- unique(data$diluent)
  drugs <- unique(data$drug)
  
  # Check on expected controls
  if(doseR == FALSE & FALSE %in% (controls %in% drugs)) {
    message(glue::glue("ERROR: the controls are not configured as expected in the design file doseR = FALSE."))
    example <- tibble::tibble(drug = c("DMSO", "water", "death juice", "seizure sauce"),
                              concentration_um = c(0, 0, 10, 100),
                              diluent = c("DMSO", "water", "DMSO", "water"))
    stop(message(paste0(capture.output(knitr::kable(example)), collapse = "\n")))
  }
  
  # Warn about doseR=T
  if(doseR == TRUE & TRUE %in% (controls %in% drugs)){
    message(glue::glue("ERROR: the controls are not configured as expected in the design file for doseR = TRUE."))
    example <- tibble::tibble(drug = c("death juice", "death juice", "seizure sauce", "seizure sauce"),
                              concentration_um = c(0, 10, 0, 100),
                              diluent = c("DMSO", "DMSO", "water", "water"))
    stop(message(paste0(capture.output(knitr::kable(example)), collapse = "\n")))
  }
  
  # Get any user flags from data
  uf1 <- names(data %>% dplyr::select(contains("_WellFlag")))
  
  # Filter wells if needed
  if(WF == "filter") {
    if(length(uf1) == 0) {
      message(glue::glue("No flagged wells detected."))
      d <- data
    } else {
      d <- filterWF_custom(data = data, rmVars = TRUE)  # Use custom function instead
      message(glue::glue("Flagged wells are filtered from the data."))
    }
  } else {
    d <- data
    message(glue::glue('Flagged wells are NOT being filtered from the data. This is NOT the recommended approach. Please consider WF = "filter".'))
  }
  
  # Grouped by message
  group_by <- paste(names(d %>% dplyr::select(!!!grouping_vars)), collapse = ", ")
  message(glue::glue("The data are grouped by, {group_by}."))
  
  # Check that all vars exist in data
  missing_vars <- vars[!vars %in% names(d)]
  if(length(missing_vars) > 0) {
    stop(glue::glue("The following variables are not in the data: {paste(missing_vars, collapse = ', ')}"))
  }
  
  # Setup means
  if(doseR == FALSE) {
    # Get control values for each group
    control_values <- d %>%
      dplyr::filter(drug %in% controls) %>% # filter to control wells
      dplyr::group_by(!!!grouping_vars) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(vars), ~mean(., na.rm = TRUE), .names = "control_{.col}"), .groups = "drop")
  }
  
  if(doseR == TRUE) {
    # Get control values for each group
    control_values <- d %>%
      dplyr::filter(concentration_um == 0) %>% # filter to control wells
      dplyr::group_by(!!!grouping_vars) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(vars), ~mean(., na.rm = TRUE), .names = "control_{.col}"), .groups = "drop")
  }
  
  # Calculate the delta for each variable
  suppressMessages(
    delta <- d %>%
      dplyr::left_join(., control_values, by = names(d %>% dplyr::select(!!!grouping_vars)))
  )
  
  # Loop through each variable to calculate delta
  for(var in vars) {
    control_col <- paste0("control_", var)
    delta_col <- paste0(var, "_delta")
    
    delta <- delta %>%
      dplyr::mutate(!!delta_col := .data[[var]] - .data[[control_col]])
  }
  
  # Message
  message("The mean control value within groups has been subtracted from the well summary statistics:")
  message(glue::glue("{paste(vars, collapse = '\n')}"))
  
  # Return the data
  return(delta)
}

#DRC2 is the only one with multiple bleaches
if(a == "20250424_Cbriggsae_DRC2" || combine_drc || a == "20250619_Cbriggsae_DRC4"){
  
  # Regress the effect of independent bleaches for each drug
  message("Running regEff_custom function (fw_Cb_ns_length data) -> reg_CB_ns_by_conc")
  
  reg_Cb_ns_by_conc <- regEff_custom(
    data = fw_Cb_ns_length,
    strain, drug,
    d.var = median_wormlength_um,
    c.var = bleach
  )
  
  View(reg_Cb_ns_by_conc$d)
  
  #Look at the regression coefficients in the diagnostic plot p2
  reg_Cb_ns_by_conc$p2
  
  #Regress the effect of independent bleaches for each drug on data filtered for width outlier
  message("Running regEff_custom function (fw_Cb_ns data) -> reg_CB_ns_by_conc_width")
  reg_Cb_ns_by_conc_width <- regEff_custom(
    data = fw_Cb_ns,
    strain, drug,
    d.var = median_wormwidth_um,
    c.var = bleach
  )
  
  #View(reg_Cb_ns_by_conc_width$d)
  
  #Look at the regression coefficients in the diagnostic plot p2
  reg_Cb_ns_by_conc_width$p2
  
  # Calculate delta for length
  message("running delta for length -> delta_length_Cb_ns")
  del_length_Cb_ns <- delta_custom(
    data = reg_Cb_ns_by_conc$d,
    drug, strain,
    WF = "ignore",
    doseR = TRUE,
    vars = c("median_wormlength_um_reg"))
  
  # Plot length delta
  delta_length <- checkEff_custom(
    data = del_length_Cb_ns %>% filter(concentration_um > 0),
    drug, strain,
    x = concentration_um,
    y = median_wormlength_um_reg_delta,
    size = 1.5,
    fill = bleach, 
    scales = "free_y"  # Changed from "free_x"
  ) +
    scale_x_log10( breaks = axis_breaks(b),
                   labels = axis_labels(b)) +
    annotation_logticks(sides = "b")
  
  #Calculate delta for width 
  message("running delta for width -> delta_width_Cb_ns")
  del_width_Cb_ns <- delta_custom(
    data = reg_Cb_ns_by_conc_width$d,
    drug, strain,
    WF = "ignore",
    doseR = TRUE,
    vars = c("median_wormwidth_um_reg"))
  
  # Plot width delta
  delta_width <- checkEff_custom(
    data = del_width_Cb_ns %>% filter(concentration_um > 0),
    drug, strain,
    x = concentration_um,
    y = median_wormwidth_um_reg_delta,
    size = 1.5,
    fill = bleach,
    scales = "free_x"
  ) + ggplot2::labs(title = "Delta (Width)") +
    scale_x_log10( breaks = axis_breaks(b),
                   labels = axis_labels(b))
  
}else{
  
  # Calculate delta for length
  message("running delta for length -> delta_length_Cb_ns")
  del_length_Cb_ns <- delta_custom(
    data = fw_Cb_ns_length,
    drug, strain,
    WF = "ignore",
    doseR = TRUE,
    vars = c("median_wormlength_um"))
  
  #Calculate delta for width 
  message("running delta for width -> delta_width_Cb_ns")
  del_width_Cb_ns <- delta_custom(
    data = fw_Cb_ns,
    drug, strain,
    WF = "ignore",
    doseR = TRUE,
    vars = c("median_wormwidth_um"))
  
  # Plot length delta
  delta_length <- checkEff_custom(
    data = del_length_Cb_ns %>% filter(concentration_um > 0),
    drug, strain,
    x = concentration_um,
    y = median_wormlength_um_delta,
    size = 1.5,
    scales = "free_y"  # Changed from "free_x"
  ) +
    scale_x_log10( breaks = axis_breaks(b),
                   labels = axis_labels(b)) +
    annotation_logticks(sides = "b")
  
  # Plot width delta
  delta_width <- checkEff_custom(
    data = del_width_Cb_ns %>% filter(concentration_um > 0),
    drug, strain,
    x = concentration_um,
    y = median_wormwidth_um_delta,
    size = 1.5,
    scales = "free_x"
  ) + ggplot2::labs(title = "Delta (Width)") +
    scale_x_log10( breaks = axis_breaks(b),
                   labels = axis_labels(b))
  
}

drug_map <- function(drug_name){
  if (drug_name == "ABZ"){
    return("ALBENDAZOLE")
  }
  if (drug_name == "CLO"){
    return("CLOSANTEL")
  }
  if(drug_name == "DEC"){
    return("DEC")
  }
  if (drug_name == "EMO"){
    return("EMODEPSIDE")
  }
  if (drug_name == "FBZ"){
    return("FENBENDAZOLE")
  }
  if (drug_name == "IVM"){
    return("IVERMECTIN")
  }
  if (drug_name == "LEV"){
    return("LEVAMISOLE")
  }
  if (drug_name == "MOX"){
    return("MOXIDECTIN")
  }
  if (drug_name == "PYR"){
    return("PYRANTEL")
    }
  if (drug_name == "MILB"){
    return ("MILBEMYCIN")
  }
}

lower_limit <- function(drug_name){
  if (drug_name == "ABZ"){
    return(c(-200,50))
  }
  if (drug_name == "CLO"){
    return(c(-600,50))
  }
  if(drug_name == "DEC"){
    return(c(-300, 0))
  }
  if (drug_name == "EMO"){
    return(c(-700,50))
  }
  if (drug_name == "FBZ"){
    return(c(-300,50))
  }
  if (drug_name == "IVM"){
    return(c(-700, 100))
  }
  if (drug_name == "LEV"){
    return(c(-600,100))
  }
  if (drug_name == "MOX"){
    return(c(-800, 100))
  }
  if (drug_name == "PYR"){
    return(c(-300, 100))
  }
  if(drug_name == "MILB"){
    return(c(-700, 200))
  }
  if (drug_name == ""){
    return(c())
  }
}

# Prepare NemaSize data
nemasize_data <- del_length_Cb_ns %>%
  dplyr::filter(drug == drug_map(b)) %>%
  dplyr::filter(concentration_um == 0) %>%
  dplyr::select(strain, median_wormlength_um) %>%
  dplyr::mutate(method = "NemaSize")

#############################
#                           #
#        DRC Plots          #
#                           #
#############################
if (a == "20250424_Cbriggsae_DRC2" || a == "20250619_Cbriggsae_DRC4"){
  
  # Plot the mean regressed worm length across concentrations with a fit line for each strain
  all_Cb_response_NS <- del_length_Cb_ns %>%
    mutate(concentration_um = as.numeric(concentration_um)) %>%
    group_by(strain, concentration_um) %>%
    summarise(mean_median_wormlength_um = mean(median_wormlength_um_reg_delta)) %>% #median_wormlength_um
    ggplot(aes(x = concentration_um, y = mean_median_wormlength_um, color = strain)) +
    geom_point()
  
  all_Cb_response_NS
  
  dose_average_NS <- del_length_Cb_ns %>%
    mutate(concentration_um = as.numeric(concentration_um)) %>%
    group_by(strain, concentration_um) %>%
    summarise(mean_median_wormlength_um_delta = mean(median_wormlength_um_reg_delta),
              se_median_wormlength_um_delta = sd(median_wormlength_um_reg_delta)) %>%
    dplyr::mutate(log.concentration = round(log(concentration_um), 2))
  
  #View(dose_average_NS )
  
}else{
  # Plot the mean regressed worm length across concentrations with a fit line for each strain
  all_Cb_response_NS <- del_length_Cb_ns %>%
    mutate(concentration_um = as.numeric(concentration_um)) %>%
    group_by(strain, concentration_um) %>%
    summarise(mean_median_wormlength_um = mean(median_wormlength_um_delta)) %>% #median_wormlength_um
    ggplot(aes(x = concentration_um, y = mean_median_wormlength_um, color = strain)) +
    geom_point()
  
  all_Cb_response_NS
  
  dose_average_NS <- del_length_Cb_ns %>%
    mutate(concentration_um = as.numeric(concentration_um)) %>%
    group_by(strain, concentration_um) %>%
    summarise(mean_median_wormlength_um_delta = mean(median_wormlength_um_delta),
              se_median_wormlength_um_delta = sd(median_wormlength_um_delta)) %>%
    dplyr::mutate(log.concentration = round(log(concentration_um), 2))
  
  #View(dose_average_NS )
}

#############################################
#                                           #
#   Make a clean version of the DRC plot    #
#                                           #
#############################################
strain_colors <- c(
  "PB420" = "#53886C",
  "ECA2666" = "#466EB4",
  "ED3102" = "#00A0E1",
  "JU3237" = "#F7D42A",
  "NIC1667" = "#7C1D6F",
  "VX34" = "#CC799D"
)

SE.DRC_NS <- dose_average_NS[which(dose_average_NS$log.concentration != -Inf),] %>% ggplot(., aes(x = concentration_um, y = mean_median_wormlength_um_delta, color = strain)) +
theme_bw(base_size = 12) +
geom_line(aes(group = strain), position = position_dodge2(width = 0.03), size = 0.75) +
  geom_pointrange(aes(ymin = mean_median_wormlength_um_delta - se_median_wormlength_um_delta,
                      ymax = mean_median_wormlength_um_delta + se_median_wormlength_um_delta),
                  position = position_dodge2(width = 0.03), size = 0.4) +
  scale_color_manual(values = strain_colors, name = "Strain") +
  scale_y_continuous("Normalized Animal Length (µm)", limits = lower_limit(b), breaks = seq(-800, 0, 100)) +
  ggtitle(drug_map(b))+ 
  scale_x_log10( breaks = axis_breaks(b),
                labels = axis_labels(b))+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "plain"),
        plot.title = element_text(face = "plain", family = "Arial", size = 20, hjust = 0),  # hjust = 0 for left align
        axis.title.x = element_text(face = "plain", family = "Arial", size = 20),
        axis.title.y = element_text(face = "plain", family = "Arial", size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20),
        axis.text = element_text(color = "black", family = "Arial", size = 20)) +
  labs(x = "Concentration (μM)", y = "Normalized Animal Length (µm)", color = expression(bold("Strain")))

print(SE.DRC_NS)

