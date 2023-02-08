# Dylan Baker 01/31/22
#source("NOAA_omics_functions.R") to use
#Contains useful importing functions for plate reader data, plot settings, and extracting standard errors from growthrates calculations

#List of packages used generally
library(reshape2)
library(tidyverse)
library(patchwork)
library(ggtext)
library(broom)
library(gsubfn)
library(nlme)
library(data.table)
library(gtable)
library(cowplot)
# You can install aomisc from GitHub
# install.packages("devtools")
# devtools::install_github("onofriAndreaPG/aomisc")
library(aomisc)

read_data <- function(file_list, skip = NULL) {
  skip <- ifelse(is.null(skip), 0, skip)
  read_csv(file_list, 
           skip = skip,
           col_types = cols(chlA_150 = col_double(),
                            PC_100  = col_double(),
                            PC_150 = col_double(),
                            N_ug_L = col_factor())) |>
    dplyr::rename(read_time = `Reading Date/Time`) |>
    #Extract plate information and species from each file to aid in isolate mapping.
    mutate(filename = file_list,
           read_day = as.numeric(str_extract(str_extract(file_list, 
                                                         pattern = "(\\d\\d|\\d).csv"), 
                                             pattern = "(\\d\\d|\\d)")),
           Plate = as.numeric(str_extract(str_extract(file_list, 
                                                      pattern = "Plate_\\d\\d|Plate_\\d"),
                                          pattern = "\\d\\d|\\d")),
           read_time = mdy_hm(read_time)) |>
    dplyr::select(-c(filename))
}

extract_dfs <- function(x, full_results) {
  num <- x
  name <- names(full_results[num])
  df <- as.data.frame(full_results[[num]]$coefficients) |>
    mutate(sample = name) |>
    dplyr::select(sample, `Std. Error`, Estimate) |>
    dplyr::rename(growthrate = Estimate,
           std_error = `Std. Error`) |>
    slice(2)
  return(df)
}

extract_dfs2 <- function(x, full_results) {
  num <- x
  name <- names(full_results[num])
  df <- as.data.frame(full_results[[num]]$par) |>
    mutate(sample = name) |>
    dplyr::select(sample, `Std. Error`, Estimate) |>
    dplyr::rename(growthrate = Estimate,
           std_error = `Std. Error`) |>
    slice(2)
  return(df)
}

#This function takes in a full dataframe of all sample data, a string for the pigment being tested, and either the letter N or P, indicating the nutrient being tested
#It splits the full dataframe into individual samples, runs regressions on each sample, and then tests whether the slopes of the triplicate samples are significantly different than the blank control slope.
#Outputs a list of dataframes with slope t-test results and slopes of each sample
slope_ttests <- function(full_data, pigment, N_or_P) {
  #Determine which nutrient we will be focusing on
  N_or_P <- ifelse(N_or_P == "N", "N_ug_L", "P_ug_L")
  
  #Turn full_data into a splitted dataframe with identified information for each sample
  splitdf <- full_data |>
    mutate(ID = paste(Well, Plate, N_or_P, Culture, State, sep="_")) %>% 
    split(., .$ID)
  
  #Create an empty list that will hold regression data
  lm_dataList <- as.list(1:length(splitdf))
  
  #Run linear regression model on sample data to produce slopes, store in empty dataframe above
  for(i in 1:length(splitdf)) {
    df <- splitdf[[i]]
    model <- lm(log(df[[pigment]]) ~ df$read_timeDays, data = df)
    
    lm_dataList[[i]] <- mutate(df, slope=coef(model)[[2]],
                               lm_std_error = coef(summary(model))[4],
                               p_val_anova = anova(model)[[5]][1])
  }
  #Create a new data frame from the splitted data to create a list of triplicates that can be tested against a blank control for significance.
  sample_lm <- bind_rows(lm_dataList) |>
    filter(Culture != "Media_CTRL", !duplicated(ID)) |> 
    dplyr::select(ID, Well, Plate, N_or_P, Culture, State, Treatment_name, slope, pigment, lm_std_error) %>%
    split(., .$Treatment_name) |>
    purrr::keep(~nrow(.) >= 3)
  
  #Create dataframe of blank samples with slopes
  blank_lm <- bind_rows(lm_dataList) |>
    filter(Culture == "Media_CTRL", !duplicated(ID)) |>
    dplyr::select(ID, slope, pigment, lm_std_error)
  
  #Empty list for t-test results to populate
  t_tests <- as.list(1:length(sample_lm))
  #Test whether the slope of sample triplicates are significantly different from the blank measurements.
  for(i in 1:length(sample_lm)) {
    df <- sample_lm[[i]]
    Culture <- unique(df$Treatment_name)
    t_tests[[i]] <- tidy(t.test(df$slope,
                                blank_lm$slope,
                                alternative = "greater")) |>
      mutate(Treatment_name = Culture,
             slope = mean(df$slope))
  }
  t_tests <- bind_rows(t_tests)
  
  return (list(t_tests,sample_lm))
}

#https://stackoverflow.com/questions/1826519/how-to-assign-from-a-function-which-returns-more-than-one-value
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}

get_nls_N <- function(df) {
  
    y = NULL;
    
    culture_list <- c(unique(df$Culture))
 
    for (i in 1:length(culture_list)) {
      
      culture_data <- df |>
        filter(Culture == culture_list[i])
      
      model_full <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks), 
                   data = culture_data, 
                   params = list(ks ~ State, mumax ~ State),
                          start = list(mumax = .5, ks = c(50,150,300)))
      
      model_full2 <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks), 
                          data = culture_data, 
                          params = list(ks ~ State - 1, mumax ~ State - 1),
                          start = list(mumax = .5, ks = c(50,150,300)))
      
      model_ks <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks),
                               start = list(mumax = c(.5), ks = c(50,150)),
                               params = list(ks ~ State, mumax ~ 1),
                               data = culture_data)
      model_mu <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks),
                                 start = list(mumax = c(.5), ks = c(50,150)),
                                 params = list(ks ~ 1, mumax ~ State),
                                 data = culture_data)
      model_simple <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks),
                           start = list(mumax = c(.5), ks = c(150)),
                           params = list(ks ~ 1, mumax ~ 1),
                           data = culture_data)
      print(anova(model_ks, model_full))
      print(anova(model_mu, model_full))
      print(anova(model_simple, model_full))
      
      tmp <- data.table(coef(summary(model_full)), keep.rownames = "term")|>
        dplyr::rename(coef = Value,
               SE = Std.Error,
               t_val = `t-value`,
               p_val = `p-value`) |>
        mutate(Culture = culture_list[i],
               psig = case_when(p_val <= 1e-3 ~ "***",
                                p_val <= 1e-2 ~ "**",
                                p_val <= 0.05 ~ "*",
                                T ~ "NS")) 
      y <- rbind(y, tmp)
      
      tmp <- data.table(coef(summary(model_full2)), keep.rownames = "term")|>
        dplyr::rename(coef = Value,
               SE = Std.Error,
               t_val = `t-value`,
               p_val = `p-value`) |>
        mutate(Culture = culture_list[i],
               p_val = NA,
               psig = NA) 
      
       y <- rbind(y, tmp)
    }
    
    model_data <- y |>
      mutate(p.adj = p.adjust(p_val, method = "fdr"),
        psig.adj = case_when(p.adj <= 1e-3 ~ "***",
                              p.adj <= 1e-2 ~ "**",
                              p.adj <= 0.05 ~ "*",
                              p.adj > 0.05 ~ "NS",
                              T ~ NA_character_))
      
    return(model_data)
}

pccx <- function(df) {
  
  #Get contrasts across xenic strains for maximum growth rate and half saturation constant (ks) values under limiting phosphorus conditions.
  df$Culture <- factor(df$Culture, levels = c("MUT_7806", "WT_7806", "PCC_9701", "NIES_843"))
  
  simple_model <- nls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks),
                       start = list(mumax = .5, ks = 2),
                       data = df)
  
  culture_full_model <- nls(growthrate_corrected~mumax[Culture]*P_ug_L/(P_ug_L + ks[Culture]), 
                            subset = State == "Xenic",
                             start = list(mumax = c(.5,.5,.6,.7), ks = c(1,2,3,4)),
                             data = df)
  
  print(anova(simple_model, culture_full_model))
  
  print(summary(culture_full_model))
  print(confint(culture_full_model))
  comparisons <- list(~ks1 - ks2, ~ks1 -ks3, ~ks1 - ks4, ~ks2 - ks3, ~ks2 - ks4, ~ks3 - ks4,
                      ~mumax1 - mumax2, ~mumax1 -mumax3, ~mumax1 - mumax4, ~mumax2 - mumax3, ~mumax2 - mumax4, ~mumax3 - mumax4
                      )
  
  y <- gnlht(culture_full_model, comparisons) 
  y <- y |>
    mutate(p.adj = p.adjust(`p-value`, method = "fdr"),
           psig.adj = case_when(p.adj <= 1e-3 ~ "***",
                                p.adj <= 1e-2 ~ "**",
                                p.adj <= 0.05 ~ "*",
                                p.adj > 0.05 ~ "NS",
                                T ~ NA_character_)) |>
    separate_wider_delim(Form,"-", names = c("Left", "Right")) |>
    mutate(Left = case_when(str_detect(Left, "1") ~ str_replace(Left, "1", "MUT_7806"),
                            str_detect(Left, '2') ~ str_replace(Left, '2', 'WT_7806'),
                            str_detect(Left, '3') ~ str_replace(Left, '3', 'PCC_9701'),
                            str_detect(Left, '4') ~ str_replace(Left, '4', 'NIES_843') ),
           Right = case_when(str_detect(Right, "1") ~ str_replace(Right, "1", "MUT_7806"),
                             str_detect(Right, '2') ~ str_replace(Right, '2', 'WT_7806'),
                             str_detect(Right, '3') ~ str_replace(Right, '3', 'PCC_9701'),
                             str_detect(Right, '4') ~ str_replace(Right, '4', 'NIES_843')),
           State = "Xenic",
           Element = "Phosphorus") |>
    unite(Left:Right, col = Form, sep = ' - ')
  return(y)
  
  
}

nccx <- function(df) {
  #Get contrasts across xenic strains for maximum growth rate and half saturation constant (ks) values under limiting nitrogen conditions.
  
  df$Culture <- factor(df$Culture, levels = c("MUT_7806", "WT_7806", "PCC_9701", "NIES_843"))
  
  # simple_model <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks),
  #                     start = list(mumax = .5, ks = 200),
  #                     data = df)
  
  culture_full_model <- nls(growthrate_corrected~mumax[Culture]*N_ug_L/(N_ug_L + ks[Culture]), 
                            subset = State == "Xenic",
                            start = list(mumax = c(.5,.5,.6,.7), ks = c(50,150,250,400)),
                            data = df)
  
  # print(anova(simple_model, culture_full_model))
  # 
  # print(summary(culture_full_model))
  # print(confint(culture_full_model))
  comparisons <- list(~ks1 - ks2, ~ks1 -ks3, ~ks1 - ks4, ~ks2 - ks3, ~ks2 - ks4, ~ks3 - ks4,
                      ~mumax1 - mumax2, ~mumax1 -mumax3, ~mumax1 - mumax4, ~mumax2 - mumax3, ~mumax2 - mumax4, ~mumax3 - mumax4
  )
  y <- gnlht(culture_full_model, comparisons) 
  y <- y |>
    mutate(p.adj = p.adjust(`p-value`, method = "fdr"),
           psig.adj = case_when(p.adj <= 1e-3 ~ "***",
                                p.adj <= 1e-2 ~ "**",
                                p.adj <= 0.05 ~ "*",
                                p.adj > 0.05 ~ "NS",
                                T ~ NA_character_)) |>
    separate_wider_delim(Form,"-", names = c("Left", "Right")) |>
    mutate(Left = case_when(str_detect(Left, "1") ~ str_replace(Left, "1", "MUT_7806"),
                            str_detect(Left, '2') ~ str_replace(Left, '2', 'WT_7806'),
                            str_detect(Left, '3') ~ str_replace(Left, '3', 'PCC_9701'),
                            str_detect(Left, '4') ~ str_replace(Left, '4', 'NIES_843') ),
           Right = case_when(str_detect(Right, "1") ~ str_replace(Right, "1", "MUT_7806"),
                             str_detect(Right, '2') ~ str_replace(Right, '2', 'WT_7806'),
                             str_detect(Right, '3') ~ str_replace(Right, '3', 'PCC_9701'),
                             str_detect(Right, '4') ~ str_replace(Right, '4', 'NIES_843')),
           State = "Xenic",
           Element = "Nitrogen") |>
    unite(Left:Right, col = Form, sep = ' - ') 
  return(y)
}
pccax <- function(df) {
  #Get contrasts across AAxenic strains for maximum growth rate and half saturation constant (ks) values under limiting phosphorus conditions.
  
  df$Culture <- factor(df$Culture, levels = c("MUT_7806", "WT_7806", "PCC_9701", "NIES_843"))
  
  simple_model <- nls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks),
                      start = list(mumax = .5, ks = 2),
                      data = df)
  
  culture_full_model <- nls(growthrate_corrected~mumax[Culture]*P_ug_L/(P_ug_L + ks[Culture]), 
                            subset = State == "Axenic",
                            start = list(mumax = c(.5,.5,.6,.7), ks = c(1,2,3,4)),
                            data = df)
  
  print(anova(simple_model, culture_full_model))
  
  print(summary(culture_full_model))
  print(confint(culture_full_model))
  comparisons <- list(~ks1 - ks2, ~ks1 -ks3, ~ks1 - ks4, ~ks2 - ks3, ~ks2 - ks4, ~ks3 - ks4,
                      ~mumax1 - mumax2, ~mumax1 -mumax3, ~mumax1 - mumax4, ~mumax2 - mumax3, ~mumax2 - mumax4, ~mumax3 - mumax4
  )
  y <- gnlht(culture_full_model, comparisons) 
  y <- y |>
    mutate(p.adj = p.adjust(`p-value`, method = "fdr"),
           psig.adj = case_when(p.adj <= 1e-3 ~ "***",
                                p.adj <= 1e-2 ~ "**",
                                p.adj <= 0.05 ~ "*",
                                p.adj > 0.05 ~ "NS",
                                T ~ NA_character_)) |>
    separate_wider_delim(Form,"-", names = c("Left", "Right")) |>
    mutate(Left = case_when(str_detect(Left, "1") ~ str_replace(Left, "1", "MUT_7806"),
                            str_detect(Left, '2') ~ str_replace(Left, '2', 'WT_7806'),
                            str_detect(Left, '3') ~ str_replace(Left, '3', 'PCC_9701'),
                            str_detect(Left, '4') ~ str_replace(Left, '4', 'NIES_843') ),
           Right = case_when(str_detect(Right, "1") ~ str_replace(Right, "1", "MUT_7806"),
                             str_detect(Right, '2') ~ str_replace(Right, '2', 'WT_7806'),
                             str_detect(Right, '3') ~ str_replace(Right, '3', 'PCC_9701'),
                             str_detect(Right, '4') ~ str_replace(Right, '4', 'NIES_843')),
           State = "Axenic",
           Element = "Phosphorus") |>
    unite(Left:Right, col = Form, sep = ' - ')
  return(y)
  
  
}

nccax <- function(df) {
  #Get contrasts across axenic strains for maximum growth rate and half saturation constant (ks) values under limiting nitrogen conditions.
  
  df$Culture <- factor(df$Culture, levels = c("MUT_7806", "WT_7806", "PCC_9701", "NIES_843"))
  
  # simple_model <- gnls(growthrate_corrected~mumax*N_ug_L/(N_ug_L + ks),
  #                     start = list(mumax = .5, ks = 200),
  #                     data = df)
  
  culture_full_model <- nls(growthrate_corrected~mumax[Culture]*N_ug_L/(N_ug_L + ks[Culture]), 
                            subset = State == "Axenic",
                            start = list(mumax = c(.5,.5,.6,.7), ks = c(50,150,250,400)),
                            data = df)
  
  # print(anova(simple_model, culture_full_model))
  # 
  # print(summary(culture_full_model))
  # print(confint(culture_full_model))
  comparisons <- list(~ks1 - ks2, ~ks1 -ks3, ~ks1 - ks4, ~ks2 - ks3, ~ks2 - ks4, ~ks3 - ks4,
                      ~mumax1 - mumax2, ~mumax1 -mumax3, ~mumax1 - mumax4, ~mumax2 - mumax3, ~mumax2 - mumax4, ~mumax3 - mumax4
  )
  y <- gnlht(culture_full_model, comparisons) 
  y <- y |>
    mutate(p.adj = p.adjust(`p-value`, method = "fdr"),
           psig.adj = case_when(p.adj <= 1e-3 ~ "***",
                                p.adj <= 1e-2 ~ "**",
                                p.adj <= 0.05 ~ "*",
                                p.adj > 0.05 ~ "NS",
                                T ~ NA_character_)) |>
    separate_wider_delim(Form,"-", names = c("Left", "Right")) |>
    mutate(Left = case_when(str_detect(Left, "1") ~ str_replace(Left, "1", "MUT_7806"),
                            str_detect(Left, '2') ~ str_replace(Left, '2', 'WT_7806'),
                            str_detect(Left, '3') ~ str_replace(Left, '3', 'PCC_9701'),
                            str_detect(Left, '4') ~ str_replace(Left, '4', 'NIES_843') ),
           Right = case_when(str_detect(Right, "1") ~ str_replace(Right, "1", "MUT_7806"),
                             str_detect(Right, '2') ~ str_replace(Right, '2', 'WT_7806'),
                             str_detect(Right, '3') ~ str_replace(Right, '3', 'PCC_9701'),
                             str_detect(Right, '4') ~ str_replace(Right, '4', 'NIES_843') ),
           State = "Axenic",
           Element = "Nitrogen") |>
    unite(Left:Right, col = Form, sep = ' - ')
  return(y)
  
  
}

get_nls_P <- function(df) {
   y = NULL;
   
  culture_list <- c(unique(df$Culture))
  
  for (i in 1:length(culture_list)) {
    
    culture_data <- df |>
      filter(Culture == culture_list[i])
    
    model_full <- gnls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks), 
                       data = culture_data, 
                       params = list(mumax ~ State, ks ~ State),
                       start = list(mumax = .5, ks = c(0.5,3,7)))
    model_full2 <- gnls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks), 
                       data = culture_data, 
                       params = list(mumax ~ State - 1, ks ~ State -1),
                       start = list(mumax = .5, ks = c(0.5,3,7)))
    model_mu <- gnls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks), 
                        data = culture_data, 
                        params = list(mumax ~ State, ks ~ 1),
                        start = list(mumax = .5, ks = c(0.5,3)))
    model_ks <- gnls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks), 
                     data = culture_data, 
                     params = list(mumax ~ 1, ks ~ State),
                     start = list(mumax = .5, ks = c(0.5,3)))
    model_simple <- gnls(growthrate_corrected~mumax*P_ug_L/(P_ug_L + ks), 
                         data = culture_data, 
                         params = list(mumax ~ 1, ks ~ 1),
                         start = list(mumax = .5, ks = c(0.5)))
    
    print(culture_list[i])
    print(anova(model_mu, model_full))
    print(anova(model_ks, model_full))
    print(anova(model_simple, model_full))
    
    tmp <- data.table(coef(summary(model_full)), keep.rownames = "term")|>
      dplyr::rename(coef = Value,
             SE = Std.Error,
             t_val = `t-value`,
             p_val = `p-value`) |>
      mutate(Culture = culture_list[i],
             psig = case_when(p_val <= 1e-3 ~ "***",
                              p_val <= 1e-2 ~ "**",
                              p_val <= 0.05 ~ "*",
                              T ~ "NS")) 
    y <- rbind(y, tmp)
    
    tmp <- data.table(coef(summary(model_full2)), keep.rownames = "term")|>
      dplyr::rename(coef = Value,
             SE = Std.Error,
             t_val = `t-value`,
             p_val = `p-value`) |>
      mutate(Culture = culture_list[i],
             p_val = NA,
             psig = NA) 
    
    y <- rbind(y, tmp)
  }
  
  model_data <- y |>
    mutate(p.adj = p.adjust(p_val, method = "fdr"),
           psig.adj = case_when(p.adj <= 1e-3 ~ "***",
                                p.adj <= 1e-2 ~ "**",
                                p.adj <= 0.05 ~ "*",
                                p.adj > 0.05 ~ "NS",
                                T ~ NA_character_))
  
  return(model_data)

}

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}


