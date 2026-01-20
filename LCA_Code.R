#### Preparation ####
rm(list=ls()) 
set.seed(123)

library(dplyr)
library(poLCA)
library(cluster)
library(ggplot2)
library(tidyr)
library(plotly)
library(FactoMineR)
library(RColorBrewer)
library(colorspace)
library(mclust)
library(VarSelLCM)
library(StatMatch)
library(GGally)
library(fmsb)
library(tidyr)
library(tibble)
library(purrr)
library(showtext)
library(scales)

font_add_google("Lato", "lato")
showtext_auto()

setwd("/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning")

#final_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Old Dataset/transformed_df.csv')

#### Daily weather data - dont run - old ####

#https://www.visualcrossing.com/weather-query-builder/Madrid,%20Comunidad%20de%20Madrid,%20Espa%C3%B1a/metric/last15days/#
# weather_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/weather_df.csv')
# 
# final_df$Reference.Day <- ifelse(nchar(as.character(final_df$Reference.Day))==1, paste0("0",final_df$Reference.Day,sep=""), final_df$Reference.Day)
# final_df$Reference.Month <- ifelse(nchar(as.character(final_df$Reference.Month))==1, paste0("0",final_df$Reference.Month,sep=""), final_df$Reference.Month)
# final_df$Reference.Year <- ifelse(nchar(as.character(final_df$Reference.Year))==1, paste0("0",final_df$Reference.Year,sep=""), final_df$Reference.Year)
# 
# final_df$date <- paste(paste(final_df$Reference.Year,final_df$Reference.Month,sep="-"),final_df$Reference.Day,sep="-")
# 
# final_df <- merge(final_df, weather_df, by.x = "date", by.y = "datetime", all = FALSE)
# 
# str(final_df)

#write.csv(final_df, file = "final_df_with_weather.csv", row.names = FALSE)

#### LCA test - dont run - old ####

# cluster_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/cluster_df_final.csv')
# 
# f <- cbind(sex,spanish)~1
# model <- poLCA(f,cluster_df,nclass=2)
# 
# # Create a frequency table
# reason_counts <- table(final_df$primary_reason)
# 
# # Create a barplot
# barplot(reason_counts,
#         col = rainbow(length(reason_counts)),   # different colors
#         main = "Primary Reason Distribution",
#         xlab = "Primary Reason",
#         ylab = "Count",
#         las = 2,                                # rotate x labels
#         legend.text = TRUE)                     # add legend

#### Recoding of primary reason variables - dont run ####
# # entertainment = Shopping,Leisure,Sports
# final_df$primary_reason <- ifelse(final_df$primary_reason %in% c(5, 8, 9), 5, final_df$primary_reason)
# # work = work, work-related management
# final_df$primary_reason <- ifelse(final_df$primary_reason %in% c(2, 3), 2, final_df$primary_reason)
# # other = other, home, another residence, medical
# final_df$primary_reason <- ifelse(final_df$primary_reason %in% c(1,6, 11,12), 12, final_df$primary_reason)
# 
# final_df$primary_reason <- as.numeric(factor(final_df$primary_reason, levels = sort(unique(final_df$primary_reason))))
# 
# # Create a frequency table
# reason_counts <- table(final_df$primary_reason)
# 
# # Create a barplot
# barplot(reason_counts,
#         col = rainbow(length(reason_counts)),   # different colors
#         main = "Primary Reason Distribution",
#         xlab = "Primary Reason",
#         ylab = "Count",
#         las = 2,                                # rotate x labels
#         legend.text = TRUE)                     # add legend
# 

#write.csv(final_df, file = "transformed_df_with_weather.csv", row.names = FALSE)

#### Hourly weather data - dont run - old####
# weather1 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly1.csv')
# weather2 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly2.csv')
# weather3 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly3.csv')
# weather4 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly4.csv')
# 
# weather1 <- weather1 %>% dplyr::select(datetime,temp,feelslike,precip)
# weather2 <- weather2 %>% dplyr::select(datetime,temp,feelslike,precip)
# 
# weather <- rbind(weather1, weather2, weather3, weather4)
# weather$hour <- substr(weather$datetime,12,13)
# weather$date <- substr(weather$datetime,1,10)
# 
# final_df$Reference.Day <- ifelse(nchar(as.character(final_df$Reference.Day))==1, paste0("0",final_df$Reference.Day,sep=""), final_df$Reference.Day)
# final_df$Reference.Month <- ifelse(nchar(as.character(final_df$Reference.Month))==1, paste0("0",final_df$Reference.Month,sep=""), final_df$Reference.Month)
# final_df$Reference.Year <- ifelse(nchar(as.character(final_df$Reference.Year))==1, paste0("0",final_df$Reference.Year,sep=""), final_df$Reference.Year)
# 
# final_df$date <- paste(paste(final_df$Reference.Year,final_df$Reference.Month,sep="-"),final_df$Reference.Day,sep="-")
# 
# final_df$startHour <- ifelse(nchar(as.character(final_df$start_time)) %in% c(1,2), "00", ifelse(nchar(as.character(final_df$start_time))==3, paste0("0",substr(final_df$start_time,1,1),sep=""), substr(final_df$start_time,1,2)))
# final_df$arrivalHour <- ifelse(nchar(as.character(final_df$arrival_time)) %in% c(1,2), "00", ifelse(nchar(as.character(final_df$arrival_time))==3, paste0("0",substr(final_df$arrival_time,1,1),sep=""), substr(final_df$arrival_time,1,2)))
# 
# final_df %>% dplyr::select(start_time,startHour,arrival_time,arrivalHour) #test
# 
# final_df <- merge(final_df, weather, by.x = "date", by.y = "date", all = FALSE)
# final_df <- final_df[final_df$hour >= final_df$startHour & final_df$hour <= final_df$arrivalHour, ]
# 
# final_df <- final_df %>% dplyr::select(-datetime,-hour)
# 
# final_df <- final_df %>%
#   group_by(hh_id, ind_id, trip_id) %>%
#   summarise(
#     temp = mean(temp, na.rm = TRUE),
#     feelslike = mean(feelslike, na.rm = TRUE),
#     precip = mean(precip, na.rm = TRUE),
#     across(
#       .cols = -c(temp, feelslike, precip), 
#       .fns = ~ first(.), 
#       .names = "{.col}"
#     ),
#     .groups = "drop"
#   )

#write.csv(final_df, file = "transformed_df_with_weather_hourly.csv", row.names = FALSE)

#### Hourly weather data - dont run ####
final_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/full_df.csv')

weather1 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly1.csv')
weather2 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly2.csv')
weather3 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly3.csv')
weather4 <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/Weather Data/hourly4.csv')

weather1 <- weather1 %>% dplyr::select(datetime,temp,feelslike,precip)
weather2 <- weather2 %>% dplyr::select(datetime,temp,feelslike,precip)

weather <- rbind(weather1, weather2, weather3, weather4)
weather$hour <- substr(weather$datetime,12,13)
weather$date <- substr(weather$datetime,1,10)

final_df$Reference.Day <- ifelse(nchar(as.character(final_df$Reference.Day))==1, paste0("0",final_df$Reference.Day,sep=""), final_df$Reference.Day)
final_df$Reference.Month <- ifelse(nchar(as.character(final_df$Reference.Month))==1, paste0("0",final_df$Reference.Month,sep=""), final_df$Reference.Month)
final_df$Reference.Year <- ifelse(nchar(as.character(final_df$Reference.Year))==1, paste0("0",final_df$Reference.Year,sep=""), final_df$Reference.Year)

final_df$date <- paste(paste(final_df$Reference.Year,final_df$Reference.Month,sep="-"),final_df$Reference.Day,sep="-")

final_df$startHour <- ifelse(nchar(as.character(final_df$start_time)) %in% c(1,2), "00", ifelse(nchar(as.character(final_df$start_time))==3, paste0("0",substr(final_df$start_time,1,1),sep=""), substr(final_df$start_time,1,2)))
final_df$arrivalHour <- ifelse(nchar(as.character(final_df$arrival_time)) %in% c(1,2), "00", ifelse(nchar(as.character(final_df$arrival_time))==3, paste0("0",substr(final_df$arrival_time,1,1),sep=""), substr(final_df$arrival_time,1,2)))

final_df %>% dplyr::select(start_time,startHour,arrival_time,arrivalHour) #test

final_df <- merge(final_df, weather, by.x = "date", by.y = "date", all = FALSE)
final_df <- final_df[final_df$hour >= final_df$startHour & final_df$hour <= final_df$arrivalHour, ]

final_df <- final_df %>%
  mutate(
    # Pad with zeros to ensure 4 characters
    start_time_str = sprintf("%04d", start_time),
    arrival_time_str = sprintf("%04d", arrival_time),
    
    # Convert to datetime with lubridate
    start_dt = hm(paste0(substr(start_time_str, 1, 2), ":", substr(start_time_str, 3, 4))),
    arrival_dt = hm(paste0(substr(arrival_time_str, 1, 2), ":", substr(arrival_time_str, 3, 4))),
    
    # Adjust for overnight (arrival before start)
    arrival_dt = if_else(arrival_dt < start_dt, arrival_dt + hours(24), arrival_dt),
    
    # Convert hour to numeric
    hour_num = as.numeric(hour),
    
    # Calculate start and end of the hour window
    hour_start = hm(paste0(sprintf("%02d", hour_num), ":00")),
    hour_end = hour_start + minutes(60),
    
    # Calculate overlap
    overlap_start = pmax(start_dt, hour_start),
    overlap_end = pmin(arrival_dt, hour_end),
    
    # Duration in minutes
    minutes_in_hour = as.numeric(overlap_end - overlap_start, units = "mins")
  )


final_df <- final_df %>% dplyr::select(-datetime,-hour,-temp.x,-feelslike.x,-precip.x,-hour,-start_time_str,-arrival_time_str,-start_dt,-arrival_dt,-hour_num,-hour_start,-hour_end,-overlap_start,-overlap_end)

final_df <- final_df %>%
  group_by(hh_id, ind_id, trip_id) %>%
  summarise(
    temp = weighted.mean(temp.y, minutes_in_hour, na.rm = TRUE),
    feelslike = weighted.mean(feelslike.y, minutes_in_hour, na.rm = TRUE),
    precip = weighted.mean(precip.y, minutes_in_hour, na.rm = TRUE),
    across(
      .cols = -c(temp.y, feelslike.y, precip.y,minutes_in_hour),
      .fns = ~ first(.),
      .names = "{.col}"
    ),
    .groups = "drop"
  )


#write.csv(final_df, file = "full_df.csv", row.names = FALSE)


#### Transformation of full_df ####
# full_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/full_df.csv')
# 
# full_df$pv_for_trip <- ifelse(full_df$pv_for_trip.==2,0,full_df$pv_for_trip.)
# full_df$spanish <- ifelse(full_df$spanish==2,0,full_df$spanish)
# full_df$male <- ifelse(full_df$sex==1,1,0)
# full_df$loaded <- full_df$loaded.
# full_df$mobility_issue <- ifelse(full_df$mobility_issue==2,0,full_df$mobility_issue)
# full_df$electric <- full_df$electric.
# full_df$members_hh <- full_df$ppl_hh
# full_df$private_parking <- full_df$private_parking.
# 
# full_df <- dplyr::select(full_df, -pv_for_trip.,-sex,-loaded.,-electric.,-ppl_hh,-private_parking.)
# 
# write.csv(full_df, file = "full_df_recoded.csv", row.names = FALSE)

#### Performing LCA ####

full_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/full_df.csv')

full_df$driving_license_binary<-full_df$driving_license_binary+1
full_df$transport_card_binary<-full_df$transport_card_binary+1
full_df$electric.<-ifelse(full_df$electric.==0,2,1)
full_df$private_parking.<-ifelse(full_df$private_parking.==0,2,1)
full_df$vehicles_hh_plus1<-full_df$vehicles_hh+1

cluster_df <- full_df %>% dplyr::select(hh_id,ind_id,trip_id,Reference.Day,Reference.Month,Reference.Year,sex,age_group,spanish,driving_license_binary,educ_group,main_activity_group,prof_status_group,transport_card_binary,mobility_issue,municip_code,ppl_hh,pp_4age_hh,vehicles_hh_plus1,private_parking.,electric.,trips_hh)

cols_to_factor <- c(
  "sex", "age_group", "spanish", "driving_license_binary", "educ_group",
  "main_activity_group", "prof_status_group", "transport_card_binary",
  "mobility_issue", "ppl_hh", "pp_4age_hh", "vehicles_hh_plus1",
  "private_parking.", "electric.", "trips_hh"
)
cluster_df[cols_to_factor] <- lapply(cluster_df[cols_to_factor], function(x) as.integer(x))
str(cluster_df)

#cluster_df <- cluster_df[sample(1:nrow(cluster_df),0.05*nrow(cluster_df)),]

f <- cbind(sex,age_group,spanish,driving_license_binary,educ_group,main_activity_group,prof_status_group,transport_card_binary,mobility_issue,ppl_hh,pp_4age_hh,vehicles_hh_plus1,private_parking.,electric.,trips_hh)~1

#define number of clusters to consider below
total_cluster_nr <- 12

inforcrit_results <- matrix(nrow=total_cluster_nr-1,ncol=3)
colnames(inforcrit_results)<-c("aic","bic","caic")
rownames(inforcrit_results)<-2:total_cluster_nr

for(i in 2:total_cluster_nr){
  # #Run the poLCA
  # model <- poLCA(f,cluster_df,nclass=i,nrep=10)
  # #assign(paste0("model_cluster_", i), model)
  # saveRDS(model, file = paste0("/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/LCA_model_cluster_",i,".rds"))
  # # inforcrit_results[i-1,1]<-model$aic
  # # inforcrit_results[i-1,2]<-model$bic
  # # inforcrit_results[i-1,3]<-(-2 * model$llik + model$npar * (log(nrow(cluster_df)) + 1))

  # #Load the already created model
  assign(paste0("model_cluster_", i), readRDS(paste0("/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/LCA_model_cluster_",i,".rds")))
  inforcrit_results[i-1,1]<-get(paste0("model_cluster_", i))$aic
  inforcrit_results[i-1,2]<-get(paste0("model_cluster_", i))$bic
  inforcrit_results[i-1,3]<-(-2 * get(paste0("model_cluster_", i))$llik + get(paste0("model_cluster_", i))$npar * (log(nrow(cluster_df)) + 1))
}

# save_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/full_df.csv')
# save_df$pred_class_lca5 <- model_cluster_5$predclass
# write.csv(save_df, file = "full_df_with_cluster.csv", row.names = FALSE)

#### LCA Information Criteria ####

inforcrit_results <- as.data.frame(inforcrit_results)
inforcrit_results$cluster <- 2:total_cluster_nr
inforcrit_long <- pivot_longer(inforcrit_results, 
                               cols = c(aic, bic, caic), 
                               names_to = "criterion", 
                               values_to = "value")


ggplot(inforcrit_long, aes(x = cluster, y = value, color = criterion)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Information Criterion Value", 
       title = "",
       color = "Criterion") +
  scale_color_manual(
    values = c("aic" = "#1b9e77", "bic" = "#d95f02", "caic" = "#7570b3"),
    labels = c("aic" = "Akaike", "bic" = "Bayesian", "caic" = "CAIC")
  ) +
  scale_x_continuous(breaks = c(2:12), labels = c(2:12)) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

inforcrit_results_percentage <- inforcrit_results %>% 
  arrange(cluster) %>% 
  mutate_at(c("aic","bic","caic"), function(x) (100*(lag(x)-x)/x)) %>%
  na.omit() 

inforcrit_percentage_long <- pivot_longer(inforcrit_results_percentage, 
                               cols = c(aic, bic, caic), 
                               names_to = "criterion", 
                               values_to = "value")

ggplot(inforcrit_percentage_long, aes(x = cluster, y = value, color = criterion)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_text(aes(label = round(value, 3)), 
            vjust = -0.7,
            hjust = -0.3,
            size = 3,
            family = "lato") +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Information Criterion Value % decrease", 
       title = "",
       color = "Criterion") +
  scale_color_manual(
    values = c("aic" = "#1b9e77", "bic" = "#d95f02", "caic" = "#7570b3"),
    labels = c("aic" = "Akaike", "bic" = "Bayesian", "caic" = "CAIC")
  ) +
  scale_x_continuous(breaks = 2:12, labels = 2:12) +
  theme_minimal(base_size = 12, base_family = "lato") +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

inforcrit_results_percentage_aic <- inforcrit_results %>%
  arrange(cluster) %>%
  mutate(across(c("aic", "bic", "caic"), ~ .x / max(.x, na.rm = TRUE)))


inforcrit_results_percentage_aic <- pivot_longer(inforcrit_results_percentage_aic, 
                                          cols = c(aic, bic, caic), 
                                          names_to = "criterion", 
                                          values_to = "value")

inforcrit_percentage_long_aic <- inforcrit_results_percentage_aic[inforcrit_results_percentage_aic$criterion == "aic", ]

ggplot(inforcrit_percentage_long_aic, aes(x = cluster, y = value)) +
  geom_line(linewidth = 0.8, color = '#1b9e77') +
  geom_point(size = 2, color = '#1b9e77') +
  geom_text(aes(label = paste0(100*round(value, 3),"%")), 
            vjust = -0.7,
            hjust = -0.3,
            size = 3,
            family = "lato") +
  theme_minimal(base_size = 12, base_family = "lato") +
  labs(
    x = "Number of Clusters",
    y = "Scaled Akaike Information Criterion",
    title = ""
  ) +
  scale_x_continuous(breaks = 2:12, labels = 2:12,limits = c(2, 12.5)) +
  scale_y_continuous(labels=function(x) (paste0(100*x,'%'))) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

#### LCA 3d ####

#clusplot(cluster_df,model_cluster_5$predclass,color=TRUE,shade=TRUE,labels=8,lines=0,cex=0.2,main="LCA plot",col.p = c("red", "blue", "green", "orange", "purple"))

cluster_var_df <- cluster_df %>% dplyr::select(sex,age_group,spanish,driving_license_binary,educ_group,main_activity_group,prof_status_group,transport_card_binary,mobility_issue,ppl_hh,pp_4age_hh,vehicles_hh_plus1,private_parking.,electric.,trips_hh)
cluster_var_df[c("sex","spanish","driving_license_binary","main_activity_group","prof_status_group","transport_card_binary","mobility_issue","private_parking.","electric.")] <- lapply(cluster_var_df[c("sex","spanish","driving_license_binary","main_activity_group","prof_status_group","transport_card_binary","mobility_issue","private_parking.","electric.")], function(x) factor(x))
cluster_var_df$vehicles_hh <- cluster_var_df$vehicles_hh_plus1-1
cluster_var_df <- cluster_var_df %>% dplyr::select(-vehicles_hh_plus1)


famd_result <- FAMD(cluster_var_df, graph = FALSE)
#pca_res <- prcomp(cluster_var_df, center = TRUE, scale. = TRUE) # old PCA method for continuous only
pca_3d <- as.data.frame(famd_result$ind$coord[,1:3]) 
pca_3d$cluster <- as.factor(model_cluster_5$predclass)

plot_ly(pca_3d, 
        x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, 
        color = ~cluster, 
        colors = c("#9932cc",  "#d95f02", "#377eb8",  "#e7298a",  "#66a61e"), 
        type = 'scatter3d',
        mode = 'markers',
        marker = list(size = 2)) %>%
  layout(#title = list(text=paste0("LCA Plot with ",round(famd_result$eig[3,3],2),"% variability explained"),font = list(family = "sans-serif", size = 20)),    
        showlegend = FALSE,
        scene = list(
            xaxis = list(title = paste0("PCA 1: ",round(famd_result$eig[1,2],2),"%"),font = list(size = 8),tickfont = list(size = 10),titlefont = list(family = "sans-serif", size = 14),tickfont = list(family = "sans-serif", size = 12)),
            yaxis = list(title = paste0("PCA 2: ",round(famd_result$eig[2,2],2),"%"),font = list(size = 8),tickfont = list(size = 10),titlefont = list(family = "sans-serif", size = 14),tickfont = list(family = "sans-serif", size = 12)),
            zaxis = list(title = paste0("PCA 3: ",round(famd_result$eig[3,2],2),"%"),font = list(size = 8),tickfont = list(size = 10),titlefont = list(family = "sans-serif", size = 14),tickfont = list(family = "sans-serif", size = 12))
)
,legend = list(
  title = list(text = "Cluster", font = list(family = "sans-serif", size = 14)),
  font = list(family = "sans-serif", size = 12))
)

# set.seed(123)
# sum <- 0
# nr_runs <- 1000
# for(i in 1:nr_runs){
#   sum <- sum + adjustedRandIndex(model_cluster_5$predclass,sample(1:5, size = 222744, replace = TRUE))
# }
# print(sum/nr_runs)

#-mean(rowSums(model_cluster_4$posterior * log(model_cluster_4$posterior + 1e-10)))


#test
# generate_cov_ellipsoid <- function(center, cov_matrix, n = 50, scale = 1.5) {
#   eig <- eigen(cov_matrix)
#   radii <- sqrt(eig$values) * scale
#   axes <- eig$vectors
#   
#   # Unit sphere
#   u <- seq(0, 2*pi, length.out = n)
#   v <- seq(0, pi, length.out = n)
#   x <- outer(cos(u), sin(v))
#   y <- outer(sin(u), sin(v))
#   z <- outer(rep(1, length(u)), cos(v))
#   
#   points <- rbind(as.vector(x), as.vector(y), as.vector(z))
#   ellipsoid_points <- axes %*% diag(radii) %*% points
#   
#   x_ell <- matrix(ellipsoid_points[1, ], n, n) + center[1]
#   y_ell <- matrix(ellipsoid_points[2, ], n, n) + center[2]
#   z_ell <- matrix(ellipsoid_points[3, ], n, n) + center[3]
#   
#   list(x = x_ell, y = y_ell, z = z_ell)
# }
# 
# # Original cluster colors
# cluster_colors <- c("1" = "red", "2" = "blue", "3" = "green", "4" = "orange")
# darker_cluster_colors <- sapply(cluster_colors, function(col) darken(col, amount = 0.4))
# 
# # Base scatter plot
# p <- plot_ly(pca_3d,
#              x = ~Dim.1, y = ~Dim.2, z = ~Dim.3,
#              color = ~as.factor(cluster),
#              colors = cluster_colors,
#              type = 'scatter3d',
#              mode = 'markers',
#              marker = list(size = 2, opacity = 0.2))
# 
# # Add ellipsoids per cluster
# for (clust in unique(pca_3d$cluster)) {
#   
#   cluster_data <- pca_3d %>% filter(cluster == clust)
#   center <- colMeans(cluster_data[, c("Dim.1", "Dim.2", "Dim.3")])
#   cov_matrix <- cov(cluster_data[, c("Dim.1", "Dim.2", "Dim.3")])
#   
#   ellipsoid <- generate_cov_ellipsoid(center, cov_matrix, n = 50, scale = 1.5)
#   
#   ellipsoid_color <- darker_cluster_colors[as.character(clust)]
#   
#   p <- p %>%
#     add_surface(x = ellipsoid$x,
#                 y = ellipsoid$y,
#                 z = ellipsoid$z,
#                 opacity = 0.3,
#                 showscale = FALSE,
#                 surfacecolor = matrix(1, nrow = 50, ncol = 50),
#                 colorscale = list(c(0,1), rep(ellipsoid_color, 2)))
# }
# 
# # Layout
# p <- p %>% layout(
#   title = paste0("LCA Plot with ", round(famd_result$eig[3,3], 2), "% variability explained"),
#   scene = list(
#     xaxis = list(title = paste0("PCA 1: ", round(famd_result$eig[1,2],2), "%")),
#     yaxis = list(title = paste0("PCA 2: ", round(famd_result$eig[2,2],2), "%")),
#     zaxis = list(title = paste0("PCA 3: ", round(famd_result$eig[3,2],2), "%"))
#   )
# )
# 
# p

#### LCA summary line ####

cluster_df_for_summary <- cluster_df

cluster_df_for_summary$Male <- ifelse(cluster_df_for_summary$sex==1,1,0)
cluster_df_for_summary$spanish <- ifelse(cluster_df_for_summary$spanish==1,1,0)
cluster_df_for_summary$driving_license_binary <- cluster_df_for_summary$driving_license_binary-1
cluster_df_for_summary$transport_card_binary<-cluster_df_for_summary$transport_card_binary-1
cluster_df_for_summary$vehicles_hh<-cluster_df_for_summary$vehicles_hh_plus1-1

cluster_df_for_summary$electric. <- ifelse(cluster_df_for_summary$electric.==2,0,1)
cluster_df_for_summary$private_parking. <- ifelse(cluster_df_for_summary$private_parking.==2,0,1)
cluster_df_for_summary$mobility_issue <- ifelse(cluster_df_for_summary$mobility_issue==2,0,1)

cluster_df_for_summary$cluster_5 <- model_cluster_5$predclass

summarystats_cluster5 <- cluster_df_for_summary %>% 
  dplyr::select(Male,age_group,spanish,driving_license_binary,educ_group,main_activity_group,prof_status_group,transport_card_binary,mobility_issue,ppl_hh,pp_4age_hh,vehicles_hh,private_parking.,electric.,trips_hh,cluster_5) %>%
  group_by(cluster_5) %>%
  summarise_if(is.numeric,mean,na.rm=TRUE)

colnames(summarystats_cluster5) <- c("Cluster","Male","Age","Spanish","License","Education","Activity","Status","Card","Issue","Size","Size above 4","Vehicles","Parking","Electric","Trips")
summarystats_cluster5[1] <- c("Elderly","Children","Small Family","Large Family","Unemployed")


cluster5_parcoord <- summarystats_cluster5[,c("Cluster","Male","Age","Spanish","License","Card","Issue","Size","Size above 4","Vehicles","Parking","Electric","Trips")]

ggparcoord(data = cluster5_parcoord,
           columns = 2:ncol(cluster5_parcoord),  # assuming 1st column is 'Cluster'
           groupColumn = "Cluster",
           scale = "std",
           showPoints = TRUE,
           alphaLines = 0.6,
           mapping = aes(color = as.factor(Cluster))) + 
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#1c91d0")) +  # Colorblind-friendly and professional
  labs(color = "Cluster",x="Variables",y="Standardized Value") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = "grey80")
  )

#split into two

ggparcoord(data = summarystats_cluster5,
           columns = 2:round(ncol(summarystats_cluster5)/2),  # assuming 1st column is 'Cluster'
           groupColumn = "Cluster",
           scale = "std",
           showPoints = TRUE,
           alphaLines = 0.6,
           mapping = aes(color = as.factor(Cluster))) + 
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#1c91d0")) +  # Colorblind-friendly and professional
  labs(color = "Cluster",x="Variables",y="Standardized Value") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = "grey80")
  )

ggparcoord(data = summarystats_cluster5,
           columns = (round(ncol(summarystats_cluster5)/2)+1):ncol(summarystats_cluster5),  # assuming 1st column is 'Cluster'
           groupColumn = "Cluster",
           scale = "std",
           showPoints = TRUE,
           alphaLines = 0.6,
           mapping = aes(color = as.factor(Cluster))) + 
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#1c91d0")) +  # Colorblind-friendly and professional
  labs(color = "Cluster",x="Variables",y="Standardized Value") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = "grey80")
  )

#### LCA summary spider plot ####
cluster5_spider <- summarystats_cluster5[,c("Cluster","Male","Age","Spanish","License","Card","Issue","Size","Size above 4","Vehicles","Parking","Electric","Trips")]

colnames(cluster5_spider) <- c("Cluster","Male","Age","Spanish","License","Transport Card","Mobility Issue","Household Size","Household Size above 4 y/o","Vehicles","Private Parking","Electric Car","Trips")

circle_coords <- function(r, n_axis = ncol(cluster5_spider) - 1){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r)
}

central_distance <- 0.2

step_1 <- map_df(seq(0, 1, 0.25) + central_distance, circle_coords) %>%
  ggplot(aes(x, y)) +
  geom_polygon(data = circle_coords(1 + central_distance), 
               alpha = 1, fill = "gray97") +
  geom_path(aes(group = r), lty = 2, alpha = 0.5)  +
  theme_void()

axis_coords <- function(n_axis){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2
  x1 <- central_distance*cos(fi)
  y1 <- central_distance*sin(fi)
  x2 <- (1 + central_distance)*cos(fi)
  y2 <- (1 + central_distance)*sin(fi)
  
  tibble(x = c(x1, x2), y = c(y1, y2), id = rep(1:n_axis, 2))
}

step_2 <- step_1 + geom_line(data = axis_coords(ncol(cluster5_spider) - 1), 
                             aes(x, y, group = id), alpha = 0.3)

text_data <- cluster5_spider %>%
  dplyr::select(-Cluster) %>%
  map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 0.25)) %>%
  mutate(r = seq(0, 1, 0.25)) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "value")

text_coords <- function(r, n_axis = ncol(cluster5_spider) - 1){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2 + 0.005*2*pi/r
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r = r - central_distance)
}

labels_data <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
  bind_cols(text_data %>% dplyr::select(-r))

labels_data$value <- round(labels_data$value,2)

rescaled_coords <- function(r, n_axis){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  tibble(r, fi) %>% mutate(x = r*cos(fi), y = r*sin(fi)) %>% dplyr::select(-fi)
}

rescaled_data <- cluster5_spider %>% 
  dplyr::mutate(across(-Cluster, scales::rescale)) %>%
  dplyr::mutate(copy = pull(., 2)) %>% 
  pivot_longer(-Cluster, names_to = "parameter", values_to = "value") %>%
  group_by(Cluster) %>%
  dplyr::mutate(coords = rescaled_coords(value + central_distance, ncol(cluster5_spider) - 1)) %>%
  unnest(cols=c(coords))

step_3 <- step_2 +
  geom_path(data = rescaled_data, 
            aes(x, y, group = Cluster, col = Cluster), 
            size = 0.7) +  
  geom_polygon(data = rescaled_data, 
                                         aes(x, y, group = Cluster, 
                                             col = Cluster, fill = Cluster), 
                                         size = 1, alpha = 0.05, show.legend = FALSE) +  
  scale_fill_manual(values = c(
  "Elderly"       = "#9932cc",  # green
  "Children"      = "#d95f02",  # orange
  "Small Family"  = "#377eb8",  # purple
  "Large Family"  = "#e7298a",  # pink/magenta
  "Unemployed"    = "#66a61e"   # olive green
  )) +
  scale_color_manual(values = c(
    "Elderly" = "#9932cc",
    "Children" = "#d95f02",
    "Small Family" = "#377eb8",
    "Large Family" = "#e7298a",
    "Unemployed" = "#66a61e"
  ))

step_4 <- step_3 + 
  geom_text(data = labels_data, aes(x, y, label = value), alpha = 1,size=5,color="black",fontface = "bold") +
  geom_text(data = text_coords(1 + central_distance + 0.2), 
            aes(x, y), label = labels_data$parameter[1:(ncol(cluster5_spider)-1)],size=5)

step_5 <- step_4 + 
  labs(col = "") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 12))

#split into two
n <- ncol(summarystats_cluster5)
half_cols <- floor((n - 1) / 2)
summarystats_cluster5_sub1 <- summarystats_cluster5[, c(1, 2:(1 + half_cols))]

circle_coords <- function(r, n_axis = ncol(summarystats_cluster5_sub1) - 1){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r)
}

central_distance <- 0.2

step_1 <- map_df(seq(0, 1, 0.25) + central_distance, circle_coords) %>%
  ggplot(aes(x, y)) +
  geom_polygon(data = circle_coords(1 + central_distance), 
               alpha = 1, fill = "gray97") +
  geom_path(aes(group = r), lty = 2, alpha = 0.5) +
  geom_polygon(data = rescaled_data, 
               aes(x, y, group = Cluster, 
                   col = Cluster, fill = Cluster), 
               size = 1, alpha = 0.05, show.legend = FALSE) +
  theme_void()

axis_coords <- function(n_axis){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2
  x1 <- central_distance*cos(fi)
  y1 <- central_distance*sin(fi)
  x2 <- (1 + central_distance)*cos(fi)
  y2 <- (1 + central_distance)*sin(fi)
  
  tibble(x = c(x1, x2), y = c(y1, y2), id = rep(1:n_axis, 2))
}

step_2 <- step_1 + geom_line(data = axis_coords(ncol(summarystats_cluster5_sub1) - 1), 
                             aes(x, y, group = id), alpha = 0.3)

text_data <- summarystats_cluster5_sub1 %>%
  select(-Cluster) %>%
  map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 0.25)) %>%
  mutate(r = seq(0, 1, 0.25)) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "value")

text_coords <- function(r, n_axis = ncol(summarystats_cluster5_sub1) - 1){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2 + 0.01*2*pi/r
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r = r - central_distance)
}

labels_data <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
  bind_cols(text_data %>% select(-r))

labels_data$value <- round(labels_data$value,2)

rescaled_coords <- function(r, n_axis){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  tibble(r, fi) %>% mutate(x = r*cos(fi), y = r*sin(fi)) %>% select(-fi)
}

rescaled_data <- summarystats_cluster5_sub1 %>% 
  mutate(across(-Cluster, scales::rescale)) %>%
  mutate(copy = pull(., 2)) %>% 
  pivot_longer(-Cluster, names_to = "parameter", values_to = "value") %>%
  group_by(Cluster) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(summarystats_cluster5_sub1) - 1)) %>%
  unnest(cols=c(coords))

step_3 <- step_2 +
  geom_path(data = rescaled_data, 
            aes(x, y, group = Cluster, col = Cluster), 
            size = 0.7)

step_4 <- step_3 + 
  geom_text(data = labels_data, aes(x, y, label = value), alpha = 0.65,size=4,color="black") +
  geom_text(data = text_coords(1 + central_distance + 0.2), 
            aes(x, y), label = labels_data$parameter[1:(ncol(summarystats_cluster5_sub1)-1)])

step_5 <- step_4 + 
  labs(col = "") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))



n <- ncol(summarystats_cluster5)
half_cols <- floor((n - 1) / 2)
summarystats_cluster5_sub2 <- summarystats_cluster5[, c(1, (half_cols+1):ncol(summarystats_cluster5))]

circle_coords <- function(r, n_axis = ncol(summarystats_cluster5_sub2) - 1){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r)
}

central_distance <- 0.2

step_1 <- map_df(seq(0, 1, 0.25) + central_distance, circle_coords) %>%
  ggplot(aes(x, y)) +
  geom_polygon(data = circle_coords(1 + central_distance), 
               alpha = 1, fill = "gray97") +
  geom_path(aes(group = r), lty = 2, alpha = 0.5) +
  geom_polygon(data = rescaled_data, 
               aes(x, y, group = Cluster, 
                   col = Cluster, fill = Cluster), 
               size = 1, alpha = 0.05, show.legend = FALSE) +
  theme_void()

axis_coords <- function(n_axis){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2
  x1 <- central_distance*cos(fi)
  y1 <- central_distance*sin(fi)
  x2 <- (1 + central_distance)*cos(fi)
  y2 <- (1 + central_distance)*sin(fi)
  
  tibble(x = c(x1, x2), y = c(y1, y2), id = rep(1:n_axis, 2))
}

step_2 <- step_1 + geom_line(data = axis_coords(ncol(summarystats_cluster5_sub2) - 1), 
                             aes(x, y, group = id), alpha = 0.3)

text_data <- summarystats_cluster5_sub2 %>%
  select(-Cluster) %>%
  map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 0.25)) %>%
  mutate(r = seq(0, 1, 0.25)) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "value")

text_coords <- function(r, n_axis = ncol(summarystats_cluster5_sub2) - 1){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2 + 0.01*2*pi/r
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r = r - central_distance)
}

labels_data <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
  bind_cols(text_data %>% select(-r))

labels_data$value <- round(labels_data$value,2)

rescaled_coords <- function(r, n_axis){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  tibble(r, fi) %>% mutate(x = r*cos(fi), y = r*sin(fi)) %>% select(-fi)
}

rescaled_data <- summarystats_cluster5_sub2 %>% 
  mutate(across(-Cluster, scales::rescale)) %>%
  mutate(copy = pull(., 2)) %>% 
  pivot_longer(-Cluster, names_to = "parameter", values_to = "value") %>%
  group_by(Cluster) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(summarystats_cluster5_sub2) - 1)) %>%
  unnest(cols=c(coords))

step_3 <- step_2 +
  geom_path(data = rescaled_data, 
            aes(x, y, group = Cluster, col = Cluster), 
            size = 0.7)

step_4 <- step_3 + 
  geom_text(data = labels_data, aes(x, y, label = value), alpha = 0.65,size=4,color="black") +
  geom_text(data = text_coords(1 + central_distance + 0.2), 
            aes(x, y), label = labels_data$parameter[1:(ncol(summarystats_cluster5_sub2)-1)])

step_5 <- step_4 + 
  labs(col = "") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

#### LCA summary baloon plot ####
educ_map <- c(
  "1" = "Primary Education or less",
  "2" = "Secondary Education",
  "3" = "Bachelors and higher",
  "4" = "Other"
)

main_activity_map <- c(
  "1" = "Employed",
  "2" = "Retired",
  "3" = "Unemployed",
  "4" = "Student"
)

prof_status_map <- c(
  "1" = "Public sector employee",
  "2" = "Private sector employee",
  "3" = "Own/family business",
  "4" = "Other",
  "5" = "Unemployed"
)

cluster_map <- c(
  "1" = "Elderly",
  "2" = "Children",
  "3" = "Small Family",
  "4" = "Large Family",
  "5" = "Unemployed"
)

# First, recode values as desired
cluster5_baloon <- cluster_df_for_summary %>% 
  mutate(
    educ_group = recode(as.character(educ_group), !!!educ_map),
    main_activity_group = recode(as.character(main_activity_group), !!!main_activity_map),
    prof_status_group = recode(as.character(prof_status_group), !!!prof_status_map),
    cluster_5 = recode(as.character(cluster_5), !!!cluster_map)
  ) %>%
  dplyr::select(educ_group, main_activity_group, prof_status_group, cluster_5) %>%
  pivot_longer(cols = -cluster_5, names_to = "variable", values_to = "value") %>%
  group_by(cluster_5, variable, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(variable) %>%
  complete(
    cluster_5,
    value
    ) %>%
  arrange(variable, cluster_5, value)


cluster_colors <- c(
  "Elderly"       = "#9932cc",
  "Children"      = "#d95f02",
  "Small Family"  = "#377eb8",
  "Large Family"  = "#e7298a",
  "Unemployed"    = "#66a61e"
)


cluster5_baloon$cluster_5 <- factor(cluster5_baloon$cluster_5, levels = names(cluster_colors))
cluster5_baloon$variable <- factor(
  cluster5_baloon$variable,
  levels = c("educ_group", "main_activity_group", "prof_status_group"),
  labels = c("Education", "Main Activity", "Professional Status")
)

ggplot(cluster5_baloon, aes(x = value, y = cluster_5, size = count, fill = cluster_5)) +
  geom_point(aes(color = cluster_5), alpha = 0.7, shape = 21) +
  geom_text(
    aes(label = ifelse(is.na(count), "", comma(count))),
    vjust = 4,
    size = 2.2,
    fontface = "bold",
    family = "Courier",
  )+
  facet_wrap(~ variable, scales = "free_x") +
  scale_size(range = c(1, 10)) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  #labs(x = NULL, y = "Cluster") +
  labs(x = NULL, y = NULL) +
  guides(fill = "none", size = "none", color = "none") +  
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1,size=7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 9, face = "bold")
  )

#### LCA entropy plot ####
entropies_summary <- matrix(nrow=total_cluster_nr-1,ncol=1)
colnames(entropies_summary) <- c("Normalized Entropy")
rownames(entropies_summary) <- 2:total_cluster_nr

for (i in 2:total_cluster_nr){
  current_model <- get(paste0('model_cluster_',i))
  #entropy_values <- apply(current_model$posterior, 1, compute_entropy)
  current_entropy <- -rowSums(current_model$posterior*log(current_model$posterior+1e-10))
  
  # entropy_plot <- ggplot(data.frame(Entropy = current_entropy), aes(x = Entropy)) +
  #   geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", alpha = 0.8) +
  #   labs(
  #     title = paste0("Distribution of Entropies for ",i," clusters"),
  #     x = "Entropy",
  #     y = "Frequency"
  #   ) +
  #   theme_minimal()
  # 
  # print(entropy_plot)
  
  entropies_summary[i-1,1] <- 1 - sum(current_entropy)/(nrow(cluster_df)*log(i))
}

entropies_summary <- round(entropies_summary, 3)
print(entropies_summary)

entropies_summary_df <- data.frame(
  k = 2:total_cluster_nr,
  Normalized.Entropy= entropies_summary
)

ggplot(entropies_summary_df, aes(x = k, y = Normalized.Entropy)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(
    title = "",
    x = "Number of Clusters",
    y = "Normalized Entropy"
  ) +
  scale_x_continuous(breaks = c(2:12), labels = c(2:12)) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))


#### LCA asw ####
asw <- matrix(nrow=total_cluster_nr-1,ncol=1)
colnames(asw) <- c("ASW")
rownames(asw) <- 2:total_cluster_nr

set.seed(123)
sample_idx <- sample(1:nrow(cluster_var_df), size = 20000)
sample_data <- cluster_var_df[sample_idx, ]
#dist_matrix <- dist(sample_data, method = "manhattan")
sample_data[sapply(sample_data, is.integer)] <- lapply(sample_data[sapply(sample_data, is.integer)], factor,ordered=TRUE)
sample_data[sapply(sample_data, is.numeric)] <- lapply(sample_data[sapply(sample_data, is.numeric)], factor,ordered=TRUE)
sample_data$educ_group <- factor(sample_data$educ_group,ordered = FALSE)
#dist_matrix <- gower.dist(sample_data)
dist_matrix <- as.matrix(daisy(sample_data, metric = "gower"))

for (i in 2:total_cluster_nr){
  current_model <- get(paste0('model_cluster_',i))
  
  sil <- silhouette(current_model$predclass[sample_idx], dist_matrix)
  
  asw[i-1,1] <- mean(sil[,"sil_width"])
}

asw_df <- data.frame(
  k = 2:total_cluster_nr,
  ASW = asw
)

ggplot(asw_df, aes(x = k, y = ASW)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_text(aes(label = round(ASW, 3)), 
            vjust = -0.7,
            hjust = -0.3,
            size = 3,
            family = "lato") +
  theme_minimal(base_size = 12, base_family = "lato") +
  labs(
    title = "",
    x = "Number of Clusters",
    y = "Average Silhouette Width"
  ) +
  scale_x_continuous(breaks = 2:12, labels = 2:12,limits = c(2, 12.5)) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggplot(inforcrit_percentage_long_aic, aes(x = cluster, y = value)) +
  geom_line(linewidth = 0.8, color = '#1b9e77') +
  geom_point(size = 2, color = '#1b9e77') +
  geom_text(aes(label = paste0(100*round(value, 3),"%")), 
            vjust = -0.7,
            hjust = -0.3,
            size = 3,
            family = "lato") +
  theme_minimal(base_size = 12, base_family = "lato") +
  labs(
    x = "Number of Clusters",
    y = "Scaled Akaike Information Criterion",
    title = ""
  ) +
  scale_x_continuous(breaks = 2:12, labels = 2:12,limits = c(2, 12.5)) +
  scale_y_continuous(labels=function(x) (paste0(100*x,'%'))) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# ## lca for 11 and 12
# set.seed(123)
# for(i in 11:12){
#   # #Run the poLCA
#   model <- poLCA(f,cluster_df,nclass=i,nrep=10)
#   saveRDS(model, file = paste0("/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/LCA_model_cluster_",i,".rds"))
# }

#### LCA posterior prob ####
model_posteriors_5 <- as.data.frame(model_cluster_5$posterior)
colnames(model_posteriors_5) <- 1:5
model_posteriors_5$cluster <- model_cluster_5$predclass

avg_post_prob_5 <- model_posteriors_5 %>% 
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean))

print(round(avg_post_prob_5,2))
  
#### K Prot graphs ####
k_prot_df <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/eskin_clusters.csv')
k_prot_cost_eskin <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/eskin_costs.csv')
k_prot_cost <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/kprototypes_costs.csv')

colnames(k_prot_cost_eskin) <- c("cluster","cost")
k_prot_cost_eskin$Method <- rep("Eskin",nrow(k_prot_cost_eskin))

colnames(k_prot_cost) <- c("cluster","cost")
k_prot_cost$Method <- rep("Normal",nrow(k_prot_cost_eskin))

k_prot_cost_agg <- rbind(k_prot_cost,k_prot_cost_eskin)

ggplot(k_prot_cost_eskin, aes(x = cluster, y = cost)) +
  geom_line(linewidth = 0.8, color="#1b9e77") +
  geom_point(size = 2,color="#1b9e77") +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Cost", 
       title = "") +
  scale_x_continuous(breaks = c(2:9), labels = c(2:9)) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

k_prot_cost_eskin_percentage <- k_prot_cost_eskin %>% 
  arrange(cluster) %>% 
  mutate_at(c("cost"), function(x) (100*(lag(x)-x)/lag(x))) %>%
  na.omit() 

ggplot(k_prot_cost_eskin_percentage, aes(x = cluster, y = cost)) +
  geom_line(linewidth = 0.8, color="#1b9e77") +
  geom_point(size = 2, color="#1b9e77") +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Cost", 
       title = "") +
  scale_x_continuous(breaks = c(2:9), labels = c(2:9)) +
  scale_y_continuous(labels=function(x) (paste0(sprintf("%.2f", x),'%'))) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

k_prot_cost_eskin_percentage_comp <- k_prot_cost_eskin %>% 
  arrange(cluster) %>% 
  mutate_at(c("cost"), ~ .x / max(.x, na.rm = TRUE)) %>%
  na.omit() 

ggplot(k_prot_cost_eskin_percentage_comp, aes(x = cluster, y = cost)) +
  geom_line(linewidth = 0.8, color = '#1b9e77') +
  geom_point(size = 2, color = '#1b9e77') +
  geom_text(aes(label = paste0(100*round(cost, 3),"%")), 
            vjust = -0.7,
            hjust = -0.3,
            size = 3,
            family = "lato") +
  theme_minimal(base_size = 12, base_family = "lato") +
  labs(
    x = "Number of Clusters",
    y = "Information Criterion Value",
    title = ""
  ) +
  scale_x_continuous(breaks = 2:12, labels = 2:12) +
  scale_y_continuous(labels=function(x) (paste0(100*x,'%'))) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )


# k prot - eksin and normal
ggplot(k_prot_cost_agg, aes(x = cluster, y = cost, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Information Criterion Value", 
       title = "",
       color = "Criterion") +
  scale_color_manual(
    values = c("Eskin" = "#1b9e77", "Normal" = "#d95f02"),
    labels = c("Eskin" = "Eskin", "Normal" = "Normal")
  ) +
  scale_x_continuous(breaks = c(2:12), labels = c(2:12)) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

k_prot_cost_agg_percentage <- k_prot_cost_agg %>% 
  arrange(Method, cluster) %>%  
  group_by(Method) %>%          
  mutate(
    cost_change_pct = 100 * (lag(cost) - cost) / cost
  ) %>%
  ungroup()

ggplot(k_prot_cost_agg_percentage, aes(x = cluster, y = cost_change_pct, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_text(aes(label = paste0(round(cost_change_pct, 1),"%")), 
            vjust = -0.7,
            hjust = -0.3,
            size = 3,
            family = "lato") +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Information Criterion Value % decrease", 
       title = "",
       color = "Criterion") +
  scale_color_manual(
    values = c("Eskin" = "#1b9e77", "Normal" = "#d95f02"),
    labels = c("Eskin" = "Eskin", "Normal" = "Normal")
  ) +
  scale_x_continuous(breaks = 3:9, labels = 3:9) +
  theme_minimal(base_size = 12, base_family = "lato") +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

max_cost <- max(k_prot_cost_agg_percentage$cost)

k_prot_cost_agg_percentage <- k_prot_cost_agg_percentage %>%
  group_by(Method) %>%
  #mutate(cost_scaled = cost / max(cost)) %>%
  mutate(cost_scaled = cost / max_cost) %>%
  ungroup()

ggplot(k_prot_cost_agg_percentage, aes(x = cluster, y = cost_scaled, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_text(aes(label = paste0(100*round(cost_scaled, 4),"%")), 
            vjust = -0.3,
            hjust = -0.2,
            size = 3,
            family = "lato",
            show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Scaled Cost", 
       title = "",
       color = "Criterion") +
  scale_color_manual(
    values = c("Normal" = "#d95f02","Eskin" = "#1b9e77"),
    labels = c("Normal" = "Standard","Eskin" = "Eskin")
  ) +
  scale_x_continuous(breaks = 2:9, labels = 2:9, limits = c(2,9.7)) +
  scale_y_continuous(labels=function(x) (paste0(100*x,'%'))) +
  theme_minimal(base_size = 12, base_family = "lato") +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

#average silhouette width kprot
asw_kprot <- matrix(nrow=max(k_prot_cost_eskin$cluster)-1,ncol=1)
colnames(asw_kprot) <- c("ASW")
rownames(asw_kprot) <- 2:max(k_prot_cost_eskin$cluster)

# taken same as before:
# sample_idx <- sample(1:nrow(cluster_var_df), size = 2000)
# sample_data <- cluster_var_df[sample_idx, ]
# #dist_matrix <- dist(sample_data, method = "manhattan")
# sample_data[sapply(sample_data, is.integer)] <- lapply(sample_data[sapply(sample_data, is.integer)], as.numeric)
# dist_matrix <- gower.dist(sample_data)

for (i in 2:max(k_prot_cost_eskin$cluster)){
  modelColName <- paste0("cluster_labels_k_", i)
  sil <- silhouette(k_prot_df[sample_idx, ][[modelColName]], dist_matrix)
  
  asw_kprot[i-1,1] <- mean(sil[,"sil_width"])
}

asw_kprot_df <- data.frame(
  k = 2:max(k_prot_cost_eskin$cluster),
  ASW = asw_kprot
)

ggplot(asw_kprot_df, aes(x = k, y = ASW)) +
  geom_line(color = "#1b9e77", size = 1) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(
    title = "",
    x = "Number of Clusters",
    y = "Average Silhouette Width"
  ) +
  scale_x_continuous(breaks = c(2:12), labels = c(2:12)) +
  theme_minimal(base_size = 12,base_family = "lato") +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
        panel.grid.minor = element_blank(),    
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

contingency_matrix_5 <- matrix(table(model_cluster_5$predclass,k_prot_df$cluster_labels_k_5),nrow = 5)
contingency_matrix_5

adjustedRandIndex(model_cluster_5$predclass,k_prot_df$cluster_labels_k_5)

#### Compare ASW ####
asw_kprot_normal <- matrix(nrow=max(k_prot_cost$cluster)-1,ncol=1)
colnames(asw_kprot_normal) <- c("ASW")
rownames(asw_kprot_normal) <- 2:max(k_prot_cost$cluster)

k_prot_df_normal <- read.csv('/Users/mateuszrozalski/Documents/Github/Seminar-in-Machine-Learning/kprot_clusters.csv')

for (i in 2:max(k_prot_cost$cluster)){
  modelColName <- paste0("cluster_labels_k_", i)
  sil <- silhouette(k_prot_df_normal[sample_idx, ][[modelColName]], dist_matrix)
  
  asw_kprot_normal[i-1,1] <- mean(sil[,"sil_width"])
}

asw_kprot_normal_df <- data.frame(
  k = 2:max(k_prot_cost$cluster),
  ASW = asw_kprot_normal
)

asw_df$Method <- rep("LCA",nrow(asw_df))
asw_kprot_df$Method <- rep("Eskin",nrow(asw_kprot_df))
asw_kprot_normal_df$Method <- rep("Normal",nrow(asw_kprot_normal_df))

asw_agg <- rbind(asw_df,asw_kprot_df,asw_kprot_normal_df)

asw_agg <- asw_agg[asw_agg$Method!="Normal",]

asw_agg$hjust_val <- ifelse(asw_agg$ASW< (-0.001), 0.1, -0.2)
asw_agg$vjust_val <- ifelse(asw_agg$ASW< (-0.001), -1, -0.3)

ggplot(asw_agg, aes(x = k, y = ASW, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_text(aes(label = round(ASW, 2),
            hjust = hjust_val,
            vjust = vjust_val),
            size = 3,
            family = "lato",
            show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Average Silhouette Width", 
       title = "",
       color = "Method") +
  scale_color_manual(
    values = c("Eskin" = "#1b9e77", "Normal" = "#d95f02", "LCA" = "steelblue"),
    labels = c("Eskin" = "Eskin K-Prototypes", "Normal" = "Normal", "LCA" = "Latent Class Analysis")
    ) +
  scale_x_continuous(breaks = 2:12, labels = 2:12, limits=c(2,12.5)) +
  theme_minimal(base_size = 12, base_family = "lato") +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank(),    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

## 3d plot for kprot
pca_3d$cluster <- as.factor(k_prot_df$cluster_labels_k_5)

plot_ly(pca_3d, 
        x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, 
        color = ~cluster, 
        colors = c("#4D4D4D",  "#17BECF", "#BCBD22","#E24A33","#8C564B"), 
        type = 'scatter3d',
        mode = 'markers',
        marker = list(size = 2)) %>%
  layout(#title = list(text=paste0("LCA Plot with ",round(famd_result$eig[3,3],2),"% variability explained"),font = list(family = "sans-serif", size = 20)),    
    showlegend = FALSE,
    scene = list(
      xaxis = list(title = paste0("PCA 1: ",round(famd_result$eig[1,2],2),"%"),font = list(size = 8),tickfont = list(size = 10),titlefont = list(family = "sans-serif", size = 14),tickfont = list(family = "sans-serif", size = 12)),
      yaxis = list(title = paste0("PCA 2: ",round(famd_result$eig[2,2],2),"%"),font = list(size = 8),tickfont = list(size = 10),titlefont = list(family = "sans-serif", size = 14),tickfont = list(family = "sans-serif", size = 12)),
      zaxis = list(title = paste0("PCA 3: ",round(famd_result$eig[3,2],2),"%"),font = list(size = 8),tickfont = list(size = 10),titlefont = list(family = "sans-serif", size = 14),tickfont = list(family = "sans-serif", size = 12))
    )
    ,legend = list(
      title = list(text = "Cluster", font = list(family = "sans-serif", size = 14)),
      font = list(family = "sans-serif", size = 12))
  )
