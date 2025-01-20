# Simulate data for upload option in shiny
library(data.table)
# Set random seed for reproducibility
set.seed(123)

# Specify the number of genes and samples
n_genes <- 100
n_samples <- 10

# Generate a matrix of normally distributed random numbers
# with mean = 10 and standard deviation = 5
data_matrix <- matrix(
  rnorm(n_genes * n_samples, mean = 10, sd = 5),
  nrow = n_genes,
  ncol = n_samples
)
rownames(data_matrix) <- paste0("Gene", 1:n_genes)
colnames(data_matrix) <- paste0("Sample", 1:n_samples)
data_matrix = data.frame(data_matrix)
data_matrix = tibble::rownames_to_column(data_matrix, "ID")

fwrite(data_matrix,"F:/Projects/Brooke/lung_repo/Shiny/Shiny_Github/Shiny_uploadData/www/example_data/omics.csv")

set.seed(123)

# Specify the number of samples
n_samples <- 10

# Simulate categorical phenotype data
phenotype_data <- data.frame(
  SampleID = paste0("Sample", 1:n_samples),
  Sex = sample(c("Male", "Female"), n_samples, replace = TRUE),
  Experimental_Group = sample(c("Control", "Treatment"), n_samples, replace = TRUE),
  Genotype = sample(c("WT", "Mutant"), n_samples, replace = TRUE),
  Histology = sample(c("Normal", "Abnormal", "Severe"), n_samples, replace = TRUE)
)
fwrite(phenotype_data,"F:/Projects/Brooke/lung_repo/Shiny/Shiny_Github/Shiny_uploadData/www/example_data/phenotype_data.csv")

# process the test data
dt = fread("F:/Projects/Brooke/lung_repo/Shiny/Shiny_Github/Shiny_uploadData/www/example_data/omics_test.csv")
dt$ID = sapply(1:length(dt$ID ),function(x){
  str_split(dt$ID[x], "_")[[1]][1]
})
dt[1:5,1:5]
df_unique <- dt[!duplicated(dt$ID), ]
fwrite(df_unique,"F:/Projects/Brooke/lung_repo/Shiny/Shiny_Github/Shiny_uploadData/www/example_data/omics_test.csv")

dt = fread("F:/Projects/Brooke/lung_repo/Shiny/Shiny_Github/Shiny_uploadData/www/example_data/omics_test2.csv")
dt$ID = sapply(1:length(dt$ID ),function(x){
  str_split(dt$ID[x], "_")[[1]][1]
})
dt[1:5,1:5]
df_unique <- dt[!duplicated(dt$ID), ]
fwrite(df_unique,"F:/Projects/Brooke/lung_repo/Shiny/Shiny_Github/Shiny_uploadData/www/example_data/omics_test2.csv")
