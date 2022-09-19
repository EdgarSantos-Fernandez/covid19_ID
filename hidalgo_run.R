# This is a rough draft of the code. 
# TO install intRinsic
# devtools::install_github('sarawade/mcclust.ext')
# #install from source
# install.packages("Packages/intRinsic_0.2.0.tar.gz", repos = NULL, type="source")

# Turn off warning-error-conversion regarding package versions
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

# Use Libraries
library(gridExtra)
library(Hmisc)
library(robustbase)
library(Rcpp)
library(RcppArmadillo)
library(intRinsic)
library(imputeTS)
library(WVPlots)
library(styler)
library(Rtsne)
library(devtools)
library(scales)
library(gginnards)
library(ggrepel)
library(reshape2)
library(tidyverse)
library(caret)
library(ggalluvial)
library(transformr)
library(gganimate)
library(glue)
library("fANCOVA")
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library(future.apply)
library(coda)
library(ggpubr)
library(lubridate)
library(imputeTS)



# Code used to prepare data
prep_data = function(stage = 0, fname = NULL){
  #Create folder for output
  folder_dir <- paste("Outputs/", fname, "/", sep = "")
  dir.create(folder_dir, showWarnings = FALSE)
  
  #Create filename for automatic saving
  input_fname <-  paste(fname,"_input_data_stage_",stage, sep = "")
  
  # Apply fitted y's as a new column and then standardise
  normalit<-function(m){
    (m - mean(m))/(sd(m))
  }

  primary_inf_data <- read.csv("./Data/owid-covid-data-june.csv")
 
  # Convert from character to date.
  primary_inf_data <- rename(primary_inf_data, Code = iso_code)
  #Filter dates
  
  #Set first and last date of period to clean
  first_date <- as_date("2020-03-01",tz = "UTC", "%Y-%m-%d")
  last_date <- as_date("2021-05-29",tz = "UTC", "%Y-%m-%d")
  primary_inf_data$date <- as_date(primary_inf_data$date,tz = "UTC", "%Y-%m-%d")#"%d/%m/%Y"
  primary_inf_data <- filter(primary_inf_data, date >= first_date & date <= last_date)
  
  #Reference dates vector
  # dates_vec <- seq(as.Date(first_date), as_date(last_date), by="days") %>% as.data.frame()
  # idx_dates <- 1:nrow(dates_vec) %>% as.data.frame()
  # date_checker <-  cbind(idx_dates, dates_vec)
  # colnames(date_checker) <- c("idx_ref", "date_ref")
  
  #         ## Procedure to select which columns to drop
  time_series_to_keep <- c('Code','date','new_cases_smoothed_per_million', 'new_deaths_smoothed_per_million','stringency_index','population') #,'retail_and_recreation','grocery_and_pharmacy','parks','transit_stations','workplaces','residential')
  ts_covid_data_prelim <- primary_inf_data %>% select(time_series_to_keep) %>% 
    filter(!grepl('OWID|TWN', Code)) %>% #Rm OWID and Taiwan because not recognised everywhere
    group_by(Code) %>%
    filter(population > 1e6) %>% #RM small countries
    select(-population) %>% 
    complete(Code, date = seq.Date(first_date, last_date, by = "day")) #Fill missing dates
  
  n_dur <- as.numeric(last_date-first_date, units="days") 
  num_stages <- 4
  stage_starts <- round(seq(1, n_dur, by = n_dur/num_stages))[-1] #Increments for stage filters later below
  
  dates <- c(first_date, first_date+stage_starts[1], first_date+stage_starts[2], first_date+stage_starts[3], last_date)
  
  # Number of NaNs in each column for each country
  Nans_per_col_per_country <- ts_covid_data_prelim %>% 
    group_by(Code) %>% 
    summarise_all(funs(sum(is.na(.))))
  
  Num_rows_per_country <- ts_covid_data_prelim %>% 
    group_by(Code) %>% 
    summarise(n = n())
  
  # #Remove "chr" column to enable division
  rows_to_divide <- Nans_per_col_per_country %>% 
    select( -Code)
  #
  # # The following array depicts the percentage of unavailable data
  data_unavailability <- Num_rows_per_country$Code %>% 
    bind_cols(rows_to_divide / Num_rows_per_country$n) %>% 
    rename(Code = ...1)
  
  # #find Countries to definitely drop (greater than 50% of variables have data unavailable)
  countries_to_drop <- (data_unavailability %>% filter(.$new_cases_smoothed_per_million > 0.05 | 
                                                         .$new_deaths_smoothed_per_million > 0.05 | 
                                                         .$stringency_index > 0.05))[,1]

  ts_covid_data <- ts_covid_data_prelim %>% 
    filter(!Code %in% countries_to_drop) #RM countries with less than 5% data available in any column

  model_str <- function(df) {
    df %>% 
      select(stringency_index) %>% 
      ts() %>% 
      na_interpolation(option = "linear")
  }
  model_cases <- function(df) {
    df %>% 
      select(new_cases_smoothed_per_million) %>% 
      ts() %>% 
      na_interpolation(option = "linear")
  }
  model_deaths <- function(df) {
    df %>% 
      select(new_deaths_smoothed_per_million) %>% 
      ts() %>% 
      na_interpolation(option = "linear")
  }

  models <- ts_covid_data %>%
    tidyr::nest(data = c(date, new_cases_smoothed_per_million, new_deaths_smoothed_per_million, stringency_index))  %>%
    dplyr::mutate(
      # perform imputation on dataset
      str_smooth = purrr::map(data, model_str)
    ) %>% 
    dplyr::mutate(
      # perform imputation on dataset
      new_cases_clean = purrr::map(data, model_cases)
    ) %>%
    dplyr::mutate(
      # perform imputation on dataset
      new_deaths_clean = purrr::map(data, model_deaths)
    )

  results <- models %>%
    tidyr::unnest() %>%
    select(-c(stringency_index,
              new_cases_smoothed_per_million, 
              new_deaths_smoothed_per_million)) %>%
    group_by(Code) %>%
    mutate_each(funs(normalit), new_cases_clean, new_deaths_clean, str_smooth)

  #Separate into different stages.
  first_stage <- results %>% filter(date >= first(date) , date <= first(date)+stage_starts[1])
  second_stage <- results %>% filter(date >= first(date)+stage_starts[1] , date <= first(date)+stage_starts[2])
  third_stage <- results %>% filter(date >= first(date)+stage_starts[2] , date <= first(date)+stage_starts[3])
  fourth_stage <- results %>% filter(date >= first(date)+stage_starts[3] , date <= last(date))

    #Now transform each variable
  if(stage == 1){
    dataset_choice <- first_stage
    stage_dates <- c(dates[1], dates[2])
  }else if(stage == 2){
    dataset_choice <- second_stage
    stage_dates <- c(dates[2],dates[3])
  }else if(stage == 3){
    dataset_choice <- third_stage
    stage_dates <- c(dates[3],dates[4])
  }else if(stage == 4){
    dataset_choice <- fourth_stage
    stage_dates <- c(dates[4], dates[5])
  }else if(stage == 0){
    dataset_choice <- results
    stage_dates <- c(dates[1],dates[5])
  }
  
  input_data <- dataset_choice %>% 
    pivot_wider(names_from = date, values_from = -c(Code,date))
  
  output_list <- mget(c("input_data", "folder_dir", "stage_dates"))
  saveRDS(output_list, file = paste(folder_dir,fname,"_input_data",".rds",sep = ""))
  
  return(output_list)
}

#Function to run Hidalgo
run_hidalgo = function(num_sim = 30000, num_burnin = 10, covid_data = NULL, fname = NULL, folder_dir = NULL, thinning = 10, date_range = NULL){
  
  #prepare data
  final_data_clean <- covid_data
  
  ##Compute Intrinsic Dimension
  ID_DATA <- final_data_clean %>% ungroup() %>% select(-Code) #remove country class
  
  hid_fit_L2000 <- intRinsic::Hidalgo(X = as.matrix(ID_DATA),
                                      L = 6,
                                      alpha_Dirichlet = .05,
                                      nsim = num_sim,
                                      burn_in = num_burnin,
                                      Prior_type = "Truncated",
                                      thinning = 10)
  
  
  Hidalgo_raw_chains_gg <- function(output){
    cmm <- (apply(output$intrinsic.dimension, 2, function(x) cumsum(x)/seq_along(x)  ))
    D  <- reshape2::melt(output$intrinsic.dimension)
    D1 <- reshape2::melt(cmm)
    raw_chain_plot <- ggplot()+
      geom_line(data=D,aes(x=.data$Var1,
                           y=.data$value,
                           group=.data$Var2),col="gray",alpha=.5)+
      theme_bw()+
      ylab("Raw MCMC")+
      xlab("Iteration")+
      geom_line(data=D1,aes(x=.data$Var1,
                            y=.data$value,
                            group=.data$Var2,
                            col=factor(.data$Var2)),
                alpha=1,lwd=1)+
      scale_color_viridis_d(option="inferno",end = .7)+
      theme(legend.position = "none")
    return(raw_chain_plot)
  }
  
  raw_chains_gg <- Hidalgo_raw_chains_gg(output = hid_fit_L2000)
  # ggsave(filename=paste(folder_dir,fname,"_","raw_chains",".tiff", sep = ""), plot= raw_chains_gg, device= "tiff", height= NA, width= NA, units= "cm", dpi= 300, compression = "lzw")
  
  
  post_processed_output <- intRinsic::Hidalgo_postpr_ID_chains(output = hid_fit_L2000, all.chain=F)
  full_post_processed_output <- intRinsic::Hidalgo_postpr_ID_chains(output = hid_fit_L2000, all.chain=T) %>%t()
  
  h2 <- intRinsic::Coclustering_matrix(output = hid_fit_L2000, VI = T)
  
  classed_ID_DATA <- cbind.data.frame(Code = final_data_clean %>% select(Code), class = factor(h2[["optimalCL"]]), ID_DATA)
  classed_ID_DATA$class <- as.character(classed_ID_DATA$class)
  
  nd <- length(unique(classed_ID_DATA[,2])) #Number of colors to draw
  mycolors <- hue_pal()(nd) #Draw those colors (from default palette)
  classes_stage_gg <- intRinsic::Hidalgo_IDperCLASS(output = hid_fit_L2000, class = classed_ID_DATA[,2], type = "density") + 
    scale_fill_manual("Class",values = mycolors)
  classes_stage_gg$layers[[1]]$aes_params$alpha <- 0.7
  classes_stage_gg <- classes_stage_gg + 
    labs(y = "Kernel density")+
    theme(
      legend.position = "none"
    )
  # ggsave(filename=paste(folder_dir,fname,"_","class_densities",".tiff", sep = ""), plot= classes_stage_gg, device= "tiff", height= NA, width= NA, units= "cm", dpi= 300, compression = "lzw")
  
  timeseries_plot <- melt(final_data_clean) # Format to tidy - long
  input_data <- classed_ID_DATA %>% select(c(Code, class))
  time_series_data <- left_join(timeseries_plot, input_data)
  time_vec_axis_ticks <- seq(date_range[1], date_range[2], by="day") %>% 
    rep(3)
  time_series_data <- separate(time_series_data, 
           col = 2, 
           into = c("sub_dataset","Date"), 
           sep = -10, 
           remove = FALSE, 
           convert = TRUE, 
           extra = "warn", 
           fill = "warn") %>% 
      mutate(Date = as.Date(Date, format= "%Y-%m-%d"))

  str_time_series_gg <-  time_series_data %>%
    filter(sub_dataset == "str_smooth_") %>%
    ggplot(aes(x= Date, y = value,
               group = class, fill = class, color = class)) +
    stat_summary(fun.y = 'mean', geom = 'line', size=2) +
    stat_summary(fun.data = 'mean_sdl', geom = 'ribbon',
                 alpha = 0.1) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = "1 month", labels = date_format("%b %d %Y")) +
    labs(x = "", y = "CSI* ")#--> Removed because pLoS doesn't like titles: , title = "Mean and Std. Dev per class", subtitle = "ordered stringency index, new cases, new deaths")
  # 
  # # ggsave(filename=paste(folder_dir,fname,"_","timeseries",".tiff", sep = ""), plot=time_series_gg, device="tiff", height= NA, width= NA, units="cm", dpi=300, compression = "lzw")
  # 
  case_time_series_gg <-  time_series_data %>%
    filter(sub_dataset == "new_cases_clean_") %>%
    ggplot(aes(x= Date, y = value,
               group = class, fill = class, color = class)) +
    stat_summary(fun.y = 'mean', geom = 'line', size=2) +
    stat_summary(fun.data = 'mean_sdl', geom = 'ribbon',
                 alpha = 0.1) +
    theme_classic() +
    theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = "1 month", labels = date_format("%b %d %Y")) +
    labs(x = "Date", y = "New cases (pmp)*  ")#--> Removed because pLoS doesn't like titles: , title = "Mean and Std. Dev per class", subtitle = "ordered stringency index, new cases, new deaths")
  # 
  # 
  death_time_series_gg <-  time_series_data %>%
    filter(sub_dataset == "new_deaths_clean_") %>%
    ggplot(aes(x= Date, y = value,
               group = class, fill = class, color = class)) +
    stat_summary(fun.y = 'mean', geom = 'line', size=2) +
    stat_summary(fun.data = 'mean_sdl', geom = 'ribbon',
                 alpha = 0.1) +
    theme_classic() +
    theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(vjust=1)) +
    scale_x_date(date_breaks = "1 month", labels = date_format("%b %d %Y")) +
    labs(x = "", y = "New deaths (pmp)*  ")#--> Removed because pLoS doesn't like titles: , title = "Mean and Std. Dev per class", subtitle = "ordered stringency index, new cases, new deaths")
  
  # 
  #Plot avg intrinsic dimension of each country coloured by class
  country_codes <- final_data_clean %>% select(Code)
  med_id_country <-  cbind.data.frame(Code = country_codes, medID = post_processed_output[["ID_summary"]][["MED"]], class = classed_ID_DATA[,2], post_processed_output[["ID_summary"]])
  med_id_country_gg <- ggplot(med_id_country, aes(x=medID,y = Code, color=class)) + geom_point(alpha=0) +
    geom_label_repel(label=med_id_country$Code, max.overlaps = 10, box.padding = 0.2, label.padding = 0.2) +
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none")+
    labs(x = "Median posterior ID estimate", y = "Country Code")
  # ggsave(filename=paste(folder_dir,fname,"_","med_id",".tiff", sep = ""), plot=med_id_country_gg, device="tiff", height= NA, width= NA, units="cm", dpi=300, compression = "lzw")

# TSNE STUFF THAT I HAVE COMMENTED OUT NOW
  #Attain cluster_prob
  nsim <- hid_fit_L2000$Recap$nsim
  confidence <- 1/(post_processed_output$ID_summary$SD2 - post_processed_output$ID_summary$SD1)

  top_cluster_prob <- confidence %>%
    as.data.frame %>%
    cbind.data.frame(Code = country_codes)

  colnames(top_cluster_prob)[1] <- "clust_prob"
  top_cluster_prob <- as.data.frame(top_cluster_prob)
  #Create t-sne plots on the input data and mcmc chains

  #PCA to reduce it to within 50 vars
  pca_data <- ID_DATA
  res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  pc.use <- 50 # explains 93% of variance
  trunc <- res$x[,1:pc.use] %*% t(res$rotation[,1:pc.use])

  #and add the center (and re-scale) back to data
  if(res$scale != FALSE){
    trunc <- scale(trunc, center = FALSE , scale=1/res$scale)
  }
  if(res$center != FALSE){
    trunc <- scale(trunc, center = -1 * res$center, scale=FALSE)
  }
  dim(trunc); dim(pca_data) #double check for same dimensions

  #Run t-SNE on ORIGINAL INPUT data
  clust_tSNE <- Rtsne(trunc ,dims=2, perplexity=20, max_iter = 10000)
  plot_data1 <-  cbind.data.frame(Code = country_codes, medID = post_processed_output[["ID_summary"]][["MED"]], class = classed_ID_DATA[,2], ID_DATA)

  unique_plot_data <- plot_data1
  coords <- cbind.data.frame(Code = unique_plot_data$Code, class = unique_plot_data$class, medID =  unique_plot_data$medID, x_coords = clust_tSNE$Y[,1], y_coords = clust_tSNE$Y[,2])
  tsne_plot_gg <- ggplot(coords, aes(x= x_coords,y = y_coords, color=class)) + geom_point(shape=1, aes(size = medID)) +
    geom_text_repel(label=unique_plot_data$Code,fontface = "bold", max.overlaps = 50) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_size_continuous(range = c(3, 10)) +
    labs(title="t-SNE clustering of Input Data", color = "VI Cluster Group", size = "mean Intrinsic Dim")
  # ggsave(filename=paste(folder_dir,fname,"_","tsne_on_input",".tiff", sep = ""), plot=tsne_plot_gg, device="tiff", height= NA, width= NA, units="cm", dpi=300, compression = "lzw")


  #TSNE ON POSTERIOR CHAINS
  duplicated_countries <- full_post_processed_output %>% as_tibble() %>% duplicated()
  posteriors_countries <- cbind(Code = country_codes,class = classed_ID_DATA[,2], full_post_processed_output)
  unique_data <- posteriors_countries[!duplicated_countries,]

  #PCA
  pca_data <- unique_data[,-c(1,2)] %>% t()
  res <- prcomp(pca_data, center = TRUE,scale. = TRUE)

  pc.use <- 50 # explains 93% of variance
  trunc <- res$x[,1:pc.use] %*% t(res$rotation[,1:pc.use])

  #and add the center (and re-scale) back to data
  if(res$scale != FALSE){
    trunc <- scale(trunc, center = FALSE , scale=1/res$scale)
  }
  if(res$center != FALSE){
    trunc <- scale(trunc, center = -1 * res$center, scale=FALSE)
  }
  dim(trunc); dim(pca_data) #double check for same dimensions

  #Run t-SNE on PCA data of POSTERIOR CHAINS
  clust_tSNE <- Rtsne(pca_data %>% t(),dims=2, perplexity=20, max_iter = 10000)

  unique_plot_data <- med_id_country[!duplicated_countries,]
  tsne_plot_post_gg <- ggplot(coords, aes(x= x_coords,y = y_coords, color=class)) + geom_point(shape=1, aes(size = medID)) +
    geom_text_repel(label=unique_plot_data$Code,fontface = "bold", max.overlaps = 50) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
    scale_size_continuous(range = c(3, 10))+
    labs(title="t-SNE clustering of ID chains",
         color = "VI Cluster group", size = "mean Intrinsic Dim")
#   ggsave(filename=paste(folder_dir,fname,"_","tsne_on_chains",".tiff", sep = ""), plot=tsne_plot_post_gg, device="tiff", height= NA, width= NA, units="cm", dpi=300, compression = "lzw")

  ## Working with the world indicators
  indicators <- read.csv("./Data/World Development Indicators/indicators.csv", na.string = "..")
  indicators$value <- rowMeans(indicators[ , c(5:9)], na.rm=TRUE) #find 5 year average
  indicators <- indicators %>%  as_tibble() %>% select(-c(5:9)) #convert to df
  
  #find which countries match with the clean dataset from intrinsic dimension processing
  countries_to_keep <- semi_join(indicators,coords, by = c("Country.Code" = "Code"))
  
  # spread the variables out to countries and join with the 'coords' final dataset - MAIN DATA
    final_data <- countries_to_keep %>%
    select(Country.Code,Series.Name,value) %>%
    pivot_wider(names_from = Series.Name, values_from = value) %>%
    select_if(~ !any(is.na(.))) # remove any columns that had a nan or more
  
  # Add the region country is located
  regions <- read.csv("./Data/owid-covid-data.csv") %>%
    select(continent, iso_code) %>%
    unique()
  
  gdp_per_capita <- read.csv("./Data/average-real-gdp-per-capita-across-countries-and-regions.csv") %>% select(c(2,3,4))
  names(gdp_per_capita)[3] <- "GDP per capita"
  gdp_per_capita <- gdp_per_capita[!(is.na(gdp_per_capita$Code) | gdp_per_capita$Code==""), ]
  gdp_per_capita <- gdp_per_capita %>% pivot_wider(names_from = Year, values_from = "GDP per capita")
  gdp_per_capita <- gdp_per_capita[,  colSums(is.na(gdp_per_capita)) == 0]
  gdp_per_capita <- cbind(Code = gdp_per_capita[,1],"GDP per capita" = gdp_per_capita[,ncol(gdp_per_capita)])
  names(gdp_per_capita)[2] <- "GDP per capita"
  
  HDI <- read.csv("./Data/human-development-index.csv") %>% select(c(2,3,4))
  names(HDI)[3] <- "Human Development Index"
  HDI <- HDI[!(is.na(HDI$Code) | HDI$Code==""), ]
  HDI <- HDI %>% pivot_wider(names_from = Year, values_from = "Human Development Index")
  HDI <- HDI[ , colSums(is.na(HDI)) == 0]
  HDI <- cbind(Code = HDI[,1],"Human Development Index" = HDI[,ncol(HDI)])
  names(HDI)[2] <- "Human Development Index"
  
  
  main_data <- coords %>%
    right_join(final_data, by = c("Code" = "Country.Code")) %>%
    inner_join(regions, by = c("Code" = "iso_code"), copy = F) %>%
    left_join(gdp_per_capita, by = "Code") %>%
    left_join(HDI, by = "Code") %>% 
    left_join(top_cluster_prob, by = "Code") %>% 
    relocate(`clust_prob`, .after = y_coords) %>% 
    relocate(continent, .after = y_coords) %>%
    relocate(`GDP per capita`, .after = y_coords) %>%
    relocate(`Human Development Index`, .after = y_coords)
  
  world <- ne_countries(scale = "medium", returnclass = "sf") %>%
    left_join(main_data,by = c("adm0_a3" = "Code"))
  
  
  names(world)[65] <- "continent"
  
  world_class_gg <- ggplot(data = world) + theme_bw() +
    geom_sf(aes(fill = world$class), colour = NA) + 
    labs(fill = "ID manifold\n")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
    # scale_fill_discrete(labels=c("1" = expression(paste(d[1], ' = ', 12)), "2" = expression(paste(d[2], ' = ', 9))))
  world_class_gg
  # ggsave(filename=paste(folder_dir,fname,"_","world",".tiff", sep = ""), plot=world_class_gg, device="tiff", height= 25, width= 45, units="cm", dpi=300, compression = "lzw")
  
  #Combine and output
  ts <- mget(c("time_series_data"))
  id <- mget(c("med_id_country"))
  class <- mget(c("classed_ID_DATA", "top_cluster_prob"))
  shiny <- mget(c("main_data","world"))
  # gg <- mget(c("med_id_country_gg","time_series_gg", "tsne_plot_post_gg","tsne_plot_gg","world_class_gg","classes_stage_gg"))

  # 
  # id_large <- mget(c("hid_fit_L2000","h2"))
  output_list <- mget(c("ts", "id", "class","shiny"))
  
  
  #Plot publication ready graphics
  three_in_a_line <- ggarrange(str_time_series_gg,
                               case_time_series_gg,
                               death_time_series_gg,
                               labels = c("A", "B", "C"), nrow = 1)
  last_two <- ggarrange(med_id_country_gg , classes_stage_gg, labels = c("D", "E"), heights = c(1.5,1))
  combine <-  ggarrange(three_in_a_line, last_two, nrow = 2)
  figure <- ggarrange(combine, world_class_gg, labels = c("", "F"), nrow = 2, common.legend = TRUE, legend = "bottom", heights = c(1.5,1))
  folder_dir <- "Outputs/"
  ggexport(figure, width = 2000*1.5, height = 2650*1.5, filename = paste(folder_dir,"pub_ready_", fname,".tiff", sep = ""), res = 300, pointsize = 8)
  ggexport(figure, width = 2000*1.5, height = 2650*1.5, filename = paste(folder_dir,"pub_ready_", fname,".png", sep = ""), res = 300, pointsize = 8)
  # #Save RDS - NOT SAVING AS NOT NEEDED ATM AND TOO MUCH SPACE
  # saveRDS(gg, file = paste(folder_dir,"plots_",fname,".rds", sep = ""))
  # saveRDS(output_list, file = paste(folder_dir,"stdout_", fname,".rds", sep = ""))
  # saveRDS(id_large, file = paste(folder_dir,"lrgout_", fname,".rds", sep = ""))
  return(output_list)
}

## The following code is for running hidalgo on the entire dataset and analysing it.
output_data <- prep_data(stage = 0, fname = NULL)




hidalgo_output <- run_hidalgo(covid_data = output_data$input_data, date_range = output_data$stage_dates)

## Data Analysis (all stage)
joined_data_world <- left_join(hidalgo_output$shiny$world, hidalgo_output$id$med_id_country,by = c("adm0_a3" = "Code"))


id_density_manifold <- ggplot(data = joined_data_world, aes(x = joined_data_world$`medID`, fill = joined_data_world$class.y))+
  geom_density(alpha = 0.5)
id_density_manifold

# This graph below demonstrates that elderly biased Population distribution is more likely to be in lower id. Why? Let's think.
pop_over_65_density <- ggplot(data = joined_data_world, aes(x = joined_data_world$`Population ages 65 and above (% of total population)`, fill = joined_data_world$class.y))+
  geom_density(alpha = 0.5) + 
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  # scale_x_discrete(labels=c("1" = "d = 12",
  #                           "2" = "d = 9")) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = "Kernel density", x = "Population ages 65 and above (% of total population)", fill = "id manifold\n") + 
  scale_fill_discrete(name = "ID manifold", labels=c("1" = expression(paste(d[1], ' = ', 12)), 
                                                     "2" = expression(paste(d[2], ' = ', 9)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.1))
pop_over_65_density

elderly_age_prop_by_class <- joined_data_world %>%
  st_drop_geometry() %>% 
  select(c("class.y","Population ages 65 and above (% of total population)")) %>%
  na.omit() %>% 
  group_split(class.y, .keep = F)


# kolgomorov smirnov test for ensuring different distributions.
ks.test(elderly_age_prop_by_class[[1]]$`Population ages 65 and above (% of total population)`,
        elderly_age_prop_by_class[[2]]$`Population ages 65 and above (% of total population)`)
# Two-sample Kolmogorov-Smirnov test
# 
# data:  elderly_age_prop_by_class[[1]]$`Population ages 65 and above (% of total population)` and elderly_age_prop_by_class[[2]]$`Population ages 65 and above (% of total population)`
# D = 0.4761, p-value = 2.076e-06
# alternative hypothesis: two-sided NULL hypothesis can be rejected at the p = 0.001 level.

age_distribution <- joined_data_world %>% 
  select(contains("(% of total population)") & starts_with("Population ages")|
           "class.x") %>% 
  st_drop_geometry() %>% 
  na.omit() %>% 
  gather(Population_Group, Perc_tot_pop, `Population ages 0-14 (% of total population)`:`Population ages 65 and above (% of total population)`) %>% 
  group_by(class.x, Population_Group) %>% 
  summarise_at(vars(Perc_tot_pop), funs(mean(., na.rm=TRUE)))

# Plot average difference in population distribution by class
gg_age_distrib <- ggplot(age_distribution) + 
  geom_col(aes(x = class.x, y = Perc_tot_pop, fill = Population_Group), alpha = 0.7) +
  theme_bw() + 
  labs(y = "% of total population", x = "ID manifold") + 
  theme(legend.position="bottom") + 
  scale_fill_discrete(name = "Age group", labels = c("0-14", "15-64", ">65")) + 
  ylim(0,101) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_x_discrete(labels=c("1" = expression(paste(d[1], ' = ', 12)), 
                            "2" = expression(paste(d[2], ' = ', 9)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg_age_distrib
gg_pop_dist_id <- ggarrange(pop_over_65_density, gg_age_distrib)
gg_pop_dist_id
ggexport(gg_pop_dist_id, width = 4000, height = 2000, filename = paste("Outputs/","all_stages_population_dist",".tiff", sep = ""), res = 300, pointsize = 8)
ggexport(gg_pop_dist_id, width = 4000, height = 2000, filename = paste("Outputs/","all_stages_population_dist",".png", sep = ""), res = 300, pointsize = 8)


#Marginal difference in the log sd of the two classes. Class 2 ( id = 9) interestingly has higher standard deviation in new_cases_smoothed.

# Moran's I test
library(spdep)
library(tmap)

fixed_morans_test <- st_make_valid(joined_data_world) %>% 
  drop_na("medID")
nb <- dnearneigh(cbind(fixed_morans_test$x_coords, fixed_morans_test$y_coords), 0, 2.1) # Approx 1300km.
lw <- nb2listw(nb, style="W", zero.policy=TRUE)

moran.test(fixed_morans_test$medID, lw, alternative = "two.sided")
?moran.test
?nb2listw
# Moran's I test rejects the NULL hypothesis at the p < 0.001 significance level.

# Plot income level - High income distribution
ggplot(data = joined_data_world)+
  geom_sf(aes(fill = joined_data_world$medID))

## Subset manifold class by income level
codes <- list("1. High income: OECD" = "1. High income", 
               "2. High income: nonOECD" = "1. High income",
              "3. Upper middle income"  = "2. Upper middle income",
              "4. Lower middle income" =  "3. Lower middle income",
              "5. Low income" = "4. Low income")

class_by_income_level <- joined_data_world %>%
  st_drop_geometry() %>% 
  select(c("class.y","income_grp", "medID")) %>%
  na.omit() %>% 
  mutate(income_grp = recode(income_grp, !!!codes))

plot(class_by_income_level$class.y, class_by_income_level$medID)

gg_income_grp <- ggplot(class_by_income_level, aes(x = class.y, fill = income_grp)) +
  geom_bar(position = "fill", alpha = 0.9) + 
  scale_fill_brewer(type = "seq",palette = "Greens", direction = -1) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + 
  labs(y = "Proportion", x = "ID manifold", fill = "Income Group") + 
  theme(legend.position="bottom") +
  scale_x_discrete(labels=c("1" = expression(paste(d[1], ' = ', 12)), 
                            "2" = expression(paste(d[2], ' = ', 9)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gg_income_grp
ggexport(gg_income_grp, width = 3000, height = 2000, filename = paste("Outputs/","all_stages_income_dist",".tiff", sep = ""), res = 300, pointsize = 8)
ggexport(gg_income_grp, width = 3000, height = 2000, filename = paste("Outputs/","all_stages_income_dist",".png", sep = ""), res = 300, pointsize = 8)

## Stage 1
output_data_1 <- prep_data(stage = 1, fname = NULL)
hidalgo_output <- run_hidalgo(covid_data = output_data_1$input_data, date_range = output_data_1$stage_dates, fname = "stage_1")

## Stage 2
output_data_2 <- prep_data(stage = 2, fname = NULL)
hidalgo_output <- run_hidalgo(covid_data = output_data_2$input_data, date_range = output_data_2$stage_dates, fname = "stage_2")


## Stage 3
output_data_3 <- prep_data(stage = 3, fname = NULL)
hidalgo_output <- run_hidalgo(covid_data = output_data_3$input_data, date_range = output_data_3$stage_dates, fname = "stage_3")


## Stage 4
output_data_4 <- prep_data(stage = 4, fname = NULL)
hidalgo_output <- run_hidalgo(covid_data = output_data_4$input_data, date_range = output_data_4$stage_dates,  fname = "stage_4")