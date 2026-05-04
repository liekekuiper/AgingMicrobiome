rm(list = ls())
library("data.table")
library("dplyr")
library("ggrepel")
library("tidyr")
library("ComplexUpset")
library('stringr')

source('../MicrobiomeAging_functions.R')

#### Data loading and preprocessing ####
# Load all datasets into a named list
files <- list(
  FHS16s   = "FHS/Results_FHS_16s.csv",
  FHSshot  = "FHS/Results_FHS_shot.csv",
  MrOS16s  = "MrOS/Results_MrOS_16s.csv",
  MrOSshot = "MrOS/Results_MrOS_shot.csv",
  RS16s    = "RS/Results_RS_16s.csv",
  DCS16s   = "DCS/Results_DCS_16s.csv",
  LLshot   = "Lifelines_Frailty_results/intermediatefiles/Results_age_5_DAG3_2026-02-13.csv",
  SOL16s   = "SOL/Results_SOL_16s.csv",
  SOLshot  = "SOL/Results_SOL_shot.csv"
)

dataframes <- load_and_prepare_data(
  files           = files,
  expand_fn       = expand_MrOS,
  expand_datasets = c("MrOS16s", "MrOSshot")
)

#### First impression harmonization ####

# Apply mutation to all relevant dataframes
dataframes <- lapply(dataframes, mutate_summ_data)

# Calculate max N for each dataframe and combine into a single dataframe
max_N_df <- bind_rows(lapply(names(dataframes), function(dataset_name) {
  df <- dataframes[[dataset_name]]  # Access the dataframe by name
  df %>%
    group_by(Mnumber, Outcome, Datasplit, seq_type) %>%
    summarise(max_N = max(N, na.rm = TRUE), .groups = "drop") %>%
    mutate(dataset = dataset_name)  # Add the dataset name for reference
}))

# Calculate the sum of max_N for each subgroup in max_N_df
sum_max_N_df <- max_N_df %>%
  group_by(Mnumber, Outcome, Datasplit, seq_type) %>%
  summarise(
    sum_max_N = sum(max_N, na.rm = TRUE), 
    nstudies = n(), 
    .groups = "drop"
  )
sum_max_N_df_all = subset(sum_max_N_df, seq_type == 'all' & Outcome %in% c('age', 'fi'))
sum_max_N_df_all_m3 = subset(sum_max_N_df_all, Mnumber == 'm3')


# Apply function to each dataset and store results
sig_risk_list <- lapply(names(dataframes), function(name) {
  df <- dataframes[[name]]
  df <- df %>% filter(Datasplit == "all")  # Apply filter before processing
  get_risk(df, name)
})

# Combine all datasets into one binary matrix
all_risk_data <- Reduce(function(x, y) full_join(x, y, by = "Variable"), sig_risk_list)

# Replace NA values with 0 (absence of variable in that dataset)
all_risk_data[is.na(all_risk_data)] <- 0

# Remove Variable column as it is not needed for the UpSet plot
rownames(all_risk_data) = all_risk_data$Variable
all_risk_data <- all_risk_data %>% select(-Variable)

# Ensure all columns are numeric (0 or 1)
all_risk_data[] <- lapply(all_risk_data, function(x) as.integer(x > 0))
all_risk_data = dplyr::select(all_risk_data, !contains("all"))
colnames(all_risk_data) = gsub('shot', ' MG', colnames(all_risk_data))
colnames(all_risk_data) = gsub('16s', ' 16S', colnames(all_risk_data))
colnames(all_risk_data) = gsub('SOL', 'HCHS/SOL', colnames(all_risk_data))

# Define the risk level columns
risk_vars <- colnames(all_risk_data)
risk_vars_ordered <- rev(c('DCS 16S', 
                           'FHS 16S',
                           'FHS MG', 
                           'HCHS/SOL 16S',
                           'HCHS/SOL MG',
                           'LL MG',
                           'MrOS 16S',
                           'MrOS MG',
                           'RS 16S'))


UPSET_BASE <- 15

# Generate the UpSet plot
upsetplotgenus <- ComplexUpset::upset(
  all_risk_data,
  intersect = risk_vars_ordered,
  width_ratio = 0.15,
  name = "Combination",
  group_by = "degree",
  sort_intersections = "descending",
  sort_intersections_by = c("degree", "cardinality"),
  sort_sets = FALSE,
  
  matrix = (
    intersection_matrix(geom = geom_point(shape = "circle filled", size = 3)) +
      scale_color_manual(
        values = c(
          "FHS 16S"      = "#882255",
          "FHS MG"       = "#AA4499",
          "MrOS 16S"     = "#332288",
          "MrOS MG"      = "#88ccee",
          "RS 16S"       = "#117733",
          "LL MG"        = "#cc6677",
          "DCS 16S"      = "#44AA99",
          "HCHS/SOL 16S" = "#999933",
          "HCHS/SOL MG"  = "#ddcc77"
        ),
        guide = guide_legend(override.aes = list(shape = "circle"))
      )
  ),
  
  queries = list(
    upset_query(set = "FHS 16S",      fill = "#882255"),
    upset_query(set = "FHS MG",       fill = "#AA4499"),
    upset_query(set = "MrOS 16S",     fill = "#332288"),
    upset_query(set = "MrOS MG",      fill = "#88ccee"),
    upset_query(set = "RS 16S",       fill = "#117733"),
    upset_query(set = "LL MG",        fill = "#cc6677"),
    upset_query(set = "DCS 16S",      fill = "#44AA99"),
    upset_query(set = "HCHS/SOL 16S", fill = "#999933"),
    upset_query(set = "HCHS/SOL MG",  fill = "#ddcc77")
  ),
  
  base_annotations = list(
    "Number of core genera\noverlapping for combination" =
      intersection_size(text = list(size = UPSET_BASE * 0.25), counts = TRUE) +
      theme(
        text         = element_text(size = UPSET_BASE, family = "sans"),
        axis.title.y = element_text(size = UPSET_BASE, family = "sans")
      )
  ),
  
  themes = upset_default_themes(
    text = element_text(size = UPSET_BASE, family = "sans")
  ),
  
  set_sizes =
    upset_set_size(position = "left") +
    geom_text(
      aes(label = after_stat(count)),
      stat  = "count",
      hjust = -.1,
      size  = UPSET_BASE * 0.25,
      color = "white"
    ) +
    theme(
      text           = element_text(size = UPSET_BASE, family = "sans"),
      axis.title.x   = element_text(size = UPSET_BASE, family = "sans")
    ) +
    ylab("Number of core\ngenera per dataset") +
    coord_cartesian(clip = "off")
)  


upsetplotgenus <- upsetplotgenus +
  theme(plot.title = element_blank())

save(upsetplotgenus, file = 'RDS_objects/UpsetPlotGenus.rds')

ggsave(plot = upsetplotgenus, filename = "UpsetPlotGenus.png", width = 27, height = 15, units='in', bg = "white")
ggsave(
  plot = upsetplotgenus,
  filename = "UpsetPlotGenus.pdf",
  width = 27,
  height = 15,
  units = "in",
  bg = "white"
)

#
exclusive_DCS16s <- all_risk_data %>%
  filter(`DCS 16S` == 1 & rowSums(select(., -`DCS 16S`)) == 0)
exclusive_LLshot <- all_risk_data %>%
  filter(`LL MG` == 1 & rowSums(select(., -`LL MG`)) == 0)

shot_columns <- grep("MG$", colnames(all_risk_data), value = TRUE)
present_in_all_shot <- all_risk_data %>%
  filter(rowSums(select(., all_of(shot_columns))) == length(shot_columns))

s16_columns <- grep("16S$", colnames(all_risk_data), value = TRUE)
present_in_all_16s <- all_risk_data %>%
  filter(rowSums(select(., all_of(s16_columns))) == length(s16_columns))
