library(ggplot2)
library(reshape2)
library(zoo)
library(plyr)

options(bitmapType = "cairo")

args = commandArgs(trailingOnly=TRUE)

wd = "/home/lihaotao/saofiles/"
setwd(wd)


# Initialize variables
Pred_Method = vector()
npred = 0
a = 0
saofile = ""
filtered_table_final = ""
w_table = ""
plot_1 = ""
plot_2 = ""
w_script_pymol = ""
w_script_pp = ""
gene_name = ""
Transcript = ""

# Set default values if no arguments are provided
if (length(args) < 1) {
  args = c(
    "--read_saofile", "", 
    "--read_table", "",
    "--chi_ratio", "",
    "--gene_name", "SMAD6",
    "--transcript", "NM_005585.5",
    "--write_table", "",
    "--write_plot_1", "SMAD6_pathogenic_variance",
    "--write_plot_2", "SMAD6_case_control",
    "--write_script_pymol", "SMAD6_script_pymol",
    "--write_script_pp", "SMAD6_script_protein_paint"
  )
}

# Parse arguments
while (TRUE) {
  if (a * 2 + 1 >= length(args)) break
  
  arg = args[a * 2 + 1]
  
  if (arg == "--read_saofile") {
    saofile = args[a * 2 + 2]
  } else if (arg == "--read_table") {
    filtered_table_final = args[a * 2 + 2]
  } else if (arg == "--chi_ratio"){
    chi_ratio = args[a * 2 + 2]
  } else if (arg == "--write_table") {
    w_table = args[a * 2 + 2]
  } else if (arg == "--write_plot_1") {
    plot_1 = args[a * 2 + 2]
  } else if (arg == "--write_plot_2") {
    plot_2 = args[a * 2 + 2]
  } else if (arg == "--write_script_pymol") {
    w_script_pymol = args[a * 2 + 2]
  } else if (arg == "--write_script_pp") {
    w_script_pp = args[a * 2 + 2]
  } else if (arg == "--pred_method") {
    npred = npred + 1
    Pred_Method[npred] = args[a * 2 + 2]
  } else if (arg == "--gene_name") {
    gene_name = args[a * 2 + 2]
  } else if (arg == "--transcript") {
    Transcript = args[a * 2 + 2]
  } 
  a = a + 1
}

# Input file type check
if (saofile != "" && filtered_table_final != ""){
  print("You should only give one input file. Please check.")
  stop()
} else if (saofile == "" && filtered_table_final == ""){
  print("You haven't given any input file. Please check.")
  stop()
}

# Define other functions...
# Print prediction methods if any
if (npred > 0) {
  for (a in 1:npred) {
    print(Pred_Method[a])
  }
}

print(saofile)

if (saofile != "" & w_table != "") {
  # load the file and dataframe it
  df_whole = data.frame(read.table(saofile, header = TRUE, stringsAsFactors = FALSE, fill = TRUE))
  cols_with_na = which(colSums(is.na(df_whole)) > 0)
  df_whole = na.omit(df_whole)
  df_whole[,3:cols_with_na[1]-1] = lapply(df_whole[,3:cols_with_na[1]-1], as.numeric)
  head(df_whole)
  
  # create a new table to store information needed for plotting
  # store info into the table with selected prediction methods
  default_value = c("Locus","contAB","caseAB","transcript","aa_number","aa_change")
  col_names = c(default_value, Pred_Method)
  filtered_table_final = data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(filtered_table_final) = col_names
  head(filtered_table_final)
  
  trans_link = 40
  trans_aanum = 8
  trans_aaexchange = 9
  default_trans = 7
  default_num_aachange = 3
  comments_split = strsplit(df_whole$comment, "\\|")
  
  rr = 0 
  for (r in 1:length(comments_split)) {
    sublist = comments_split[[r]]
    for (c in 1:(length(sublist) - max(trans_aanum, trans_aaexchange))) {
      if (sublist[c] == Transcript) {
        aa_num = sublist[c + trans_aanum]
        aa_exchange = sublist[c + trans_aaexchange]
        if (!is.na(aa_num) && !is.na(aa_exchange) && aa_num != "" && aa_exchange != "" && nchar(aa_exchange) == default_num_aachange) {
          rr = rr + 1
          filtered_table_final[rr, "aa_number"] = as.numeric(aa_num)
          filtered_table_final[rr, "aa_change"] = aa_exchange
          filtered_table_final[rr, "Locus"] = sublist[1]
          filtered_table_final[rr, "contAB"] = df_whole[r, "contAB"]
          filtered_table_final[rr, "caseAB"] = df_whole[r, "caseAB"]
          filtered_table_final[rr, "transcript"] = Transcript
          filtered_table_final[rr, "aa_change"] = aa_exchange
          for (i in 1:length(Pred_Method)) {
            new_col = Pred_Method[i]
            filtered_table_final[rr, new_col] = df_whole[r, new_col]
          }
        }
      }
    }
  }
  if ("AM_score" %in% Pred_Method){
    for (i in 1: nrow(filtered_table_final)){
      if (filtered_table_final[i, "AM_score"] > 0){
        filtered_table_final[i, "AM_score"] = filtered_table_final[i, "AM_score"]/10
      }	
    }
  }
  print("Producing filtered_table based on the saofile")
}


if (saofile != "" & w_table != "") {
  table_generator(saofile)
  filtered_table_final[filtered_table_final == ""] = NA
  filtered_table_final = na.omit(filtered_table_final)
  
  summary(filtered_table_final)
  head(filtered_table_final)
  
  # save the table
  table_name = sprintf("/home/lihaotao/saofiles/%s.txt",w_table)
  write.table(filtered_table_final, table_name, sep = "\t", row.names = FALSE, quote = FALSE)
}

