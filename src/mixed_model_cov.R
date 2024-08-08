library(ggplot2)
library(dplyr)
library(tidyr)  

#Set working directory to source
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define the mapping of columns to time points
clone_fraction_mapping <- c(
  "T0" = "cloneFraction_df1_d1",
  "T1" = "cloneFraction_df2_d1",
  "T2" = "cloneFraction_df1_d2",
  "T3" = "cloneFraction_df2_d2",
  "T4" = "cloneFraction_df1",
  "T5" = "cloneFraction_df2"
)


dir= '../results/covsig/'
files<-list.files(dir)

# Initialize an empty list to hold the combined data
combined_data <- list()

for (file in files){
  if (file == "BV182_CD8.txt"){
    next}
  print(file)
  patient_df <- read.table(paste0(dir,file),header=TRUE)

  # Create a unique identifier by combining patient ID and TCR
  patient_df <- patient_df %>%
    mutate(id = paste0(file, "_", TCR), patient = substr(file,1,5), tcrchain=substr(file,7,9))
  
  # Normalize to first expansion
  patient_df_norm <- patient_df %>%
    mutate(cloneFraction_df1_d1 = cloneFraction_df1_d1 / cloneFraction_df2_d1,
           cloneFraction_df1_d2 = cloneFraction_df1_d2 / cloneFraction_df2_d1,
           cloneFraction_df2_d2 = cloneFraction_df2_d2 / cloneFraction_df2_d1,
           cloneFraction_df1 = cloneFraction_df1 / cloneFraction_df2_d1,
           cloneFraction_df2 = cloneFraction_df2 / cloneFraction_df2_d1,
           cloneFraction_df2_d1 = cloneFraction_df2_d1 / cloneFraction_df2_d1)
  
  # Extract the cloneFraction columns based on the mapping
  clone_fraction_data <- patient_df_norm %>%
    select(all_of(clone_fraction_mapping), id, patient,tcrchain)
  
  # Convert the wide data to long format
  clone_fraction_long <- clone_fraction_data %>%
    pivot_longer(
      cols = all_of(names(clone_fraction_mapping)),  # Convert only the mapped columns
      names_to = "time_point",  # Store time points
      values_to = "cloneFraction"  # Name of the new clone fraction column
    ) 
  
  # Append the long data frame to the combined list
  combined_data[[length(combined_data) + 1]] <- clone_fraction_long
  
}

# Combine all patient data frames into one long data frame
long_data <- bind_rows(combined_data)
ggplot(long_data,aes(x=time_point,y=cloneFraction,group=id,color=tcrchain))+geom_line()+theme_bw()


require(lme4)

# Adding binary columns based on conditions
long_data <- long_data %>%
  mutate(
    # Define the conditions for B
    B = if_else(time_point %in% c("T0", "T1", "T2", "T3", "T4", "T5"), 1, 0),
    
    # Define the conditions for B1
    B1 = if_else(time_point %in% c("T2", "T3", "T4", "T5"), 1, 0),
    
    # Define the conditions for B2
    B2 = if_else(time_point %in% c("T4", "T5"), 1, 0),
    
    # Define the conditions for E
    E = if_else(time_point %in% c("T1", "T3", "T5"), 1, 0),
    
    # Define the conditions for I
    I = if_else((time_point %in% c("T1", "T2", "T3", "T4", "T5") & 
                   (patient %in% c("BV092","BV123","BV137"))) | 
                time_point %in% c("T4", "T5") # by T4, T5 they all had an infection
                , 1, 0)
  )

long_data_nonzero_cd4 <- long_data[(long_data$cloneFraction > 0 | long_data$time_point %in% c("T0")) & long_data$tcrchain %in% c("CD4"),]

main_f <-cloneFraction ~ B1 + E + I + (1 | id)

alt_f <- c(cloneFraction ~ E + I + (1 | id),
           cloneFraction ~ B1 + I + (1 | id),
           cloneFraction ~ B1 + E + (1 | id))

# Define the mixed model
model_cd4 <- lmer(
  main_f,
  data = long_data_nonzero_cd4, REML = F
)

# Summary of the model
summary(model_cd4)
AIC(model_cd4)

for (f in alt_f){
# Define the mixed model
model_cd4_alt <- lmer(
  f,
  data = long_data_nonzero_cd4, REML = F
)

# Summary of the model
#summary(model_cd4_alt)
#AIC(model_cd4_noI)
print(f)
print((AIC(model_cd4_alt) - AIC(model_cd4))*log2(exp(1)))
}

long_data_nonzero_cd8 <- long_data[(long_data$cloneFraction > 0 | long_data$time_point %in% c("T0")) & long_data$tcrchain %in% c("CD8"),]

# Define the mixed model
model_cd8 <- lmer(
  main_f,
  data = long_data_nonzero_cd8, REML = F
)

# Summary of the model
summary(model_cd8)

for (f in alt_f){
  # Define the mixed model
  model_cd8_alt <- lmer(
    f,
    data = long_data_nonzero_cd8, REML = F
  )
  
  # Summary of the model
  #summary(model_cd4_alt)
  #AIC(model_cd4_noI)
  print(f)
  print((AIC(model_cd8_alt) - AIC(model_cd8))*log2(exp(1)))
}

