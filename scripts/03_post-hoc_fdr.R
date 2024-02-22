#####################################
## Post-Hoc Benjamini-Hochberg FDR
## Linh Nguyen                      
## nguyenllpsych@gmail.com         
## December 2023  
#####################################

# meta
library(rio)
library(here)

# function to extract p-values
p_extract <- function(path_name, sheet_name, p_names) {

  # read in excel file
  table_extracted <- rio::import(file = path_name, sheet = sheet_name, skip = 2)
  
  # extract and combine all p-value columns
  p_extracted <- c()
  for(p_name in p_names) {
    p_val <- table_extracted[[p_name]]
    p_val <- p_val[!is.na(p_val)]
    p_extracted <- append(p_extracted, p_val)
  }
  
  # ensure numeric
  p_extracted <- gsub(
    pattern = "<[0]?.001", replacement = "0",
    x = gsub(pattern = " ", replacement = "", x = p_extracted))
  p_extracted <- as.numeric(p_extracted)
  p_extracted <- p_extracted[!is.na(p_extracted)]
  
  # return result as vector of p-values
  return(p_extracted)
}

# path to samples one and two files
path_one <- paste0(here(), "/data/Tables-S1.xlsx")
path_two <- paste0(here(), "/data/Tables-S2.xlsx")

# hypothesis one and three - table 3 and 4
s1_simcor <- p_extract(path_name = path_one, sheet_name = "s1_simcor", 
                       p_names = c("...6", "...15"))
s2_simcor <- p_extract(path_name = path_two, sheet_name = "s2_simcor",
                       p_names = c("...6", "...15"))

# hypothesis two - table s1 and s2
s1_cordiff <- round(p_extract(path_name = path_one, sheet_name = "s1_cor-diff", p_names = "p"),
                    3)
s2_cordiff <- round(p_extract(path_name = path_two, sheet_name = "s2_cor-diff", p_names = "p"),
                    3)

# hypothesis four - table s3 and s4
s1_benefits_baseline <- p_extract(path_name = path_one, 
                                  sheet_name = "s1_benefits_baseline", 
                                  p_names = c("...5", "...8", "...12", "...15"))
s2_benefits_baseline <- p_extract(path_name = path_two, 
                                  sheet_name = "s2_benefits_baseline", 
                                  p_names = c("...5", "...8", "...12", "...15",
                                              "...19", "...22", "...26", "...29"))

# hypothesis five - table s5 and s6
s1_benefits_baseline_long <- p_extract(path_name = path_one, 
                                       sheet_name = "s1_benefits_baseline-long",
                                       p_names = c("...5", "...8", "...12", "...15"))
s2_benefits_baseline_long <- p_extract(path_name = path_two, 
                                       sheet_name = "s2_benefits_baseline-long",
                                       p_names = c("...5", "...8", "...12", "...15",
                                                   "...19", "...22", "...26", "...29"))

# hypothesis six - table s7 and s8
s1_benefits_long_long <- p_extract(path_name = path_one, sheet_name = "s1_benefits_long-long",
                                   p_names = c("...5", "...8", "...12", "...15"))
s2_benefits_long_long <- p_extract(path_name = path_two, sheet_name = "s2_benefits_long-long",
                                   p_names = c("...5", "...8", "...12", "...15",
                                               "...19", "...22", "...26", "...29"))

# combine all p_val with corresponding name
p_vals <- matrix(ncol = 2)
colnames(p_vals) <- c("H", "p")

for(result in paste0(c("s1_", "s2_"), 
                     rep(c("simcor", "cordiff", "benefits_baseline", 
                           "benefits_baseline_long", "benefits_long_long"), 
                         each = 2))) {
  
  # name of current hypothesis
  H <- rep(result, times = length(get(result)))
  
  # p values of current hypothesis
  p <- get(result)
  
  # append to p_vals matrix
  p_vals <- rbind(p_vals, cbind(H, p))
}

# ensure correct type
p_vals <- as.data.frame(p_vals[-1,])
p_vals$p <- as.numeric(p_vals$p)

# convert to q-val using Benjamini-Hochberg
p_vals$q <- p.adjust(p = p_vals$p, method = "BH")

# number of significant tests using p-value with FPR = 0.05
sum(p_vals$p <= 0.05)

# number of significant tests using q-value with FDR = 0.10
sum(p_vals$q <= 0.10)

# highest significant p-value under FDR
max(p_vals[which(p_vals$q <= 0.10),]$p)

# view p_value data frame
View(p_vals)