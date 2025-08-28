# Copyright (C) Lihao Tao and David Curtis
library(ggplot2)
library(reshape2)
library(zoo)
library(plyr)

options(bitmapType = "cairo")

args = commandArgs(trailingOnly=TRUE)

path = ""
pred_method = vector()
npred = 0
a = 0
file_data = ""
table_data = ""
w_table = ""
plot_1 = ""
plot_2 = ""
w_script_pymol = ""
w_script_pp = ""
gene_name = ""
transcript = ""
chi_ratio = 0 # must be set if table is used instead of sao file - number of cases over number of controls

# Set default values if no arguments are provided
if (length(args) < 1) {
  args = c(
    "--read_file", "/home/lihaotao/saofiles/UKBB.HT.forAnnot.20231112.SMAD6.txt", 
    "--set_wd", "/home/lihaotao/saofiles/",
    "--gene_name", "SMAD6",
    "--transcript", "NM_005585.5",
    "--write_table", "UKBB.HT.forAnnot.20231112.SMAD6.tab",
    "--write_plot_1", "SMAD6_pathogenic_variants,png",
    "--write_plot_2", "SMAD6_case_control.png",
    "--write_script_pymol", "SMAD6_script_pymol.txt",
    "--write_script_pp", "SMAD6_script_protein_paint.txt",
    "--pred_method", "AM_score"
  )
}

while (TRUE) {
  if (a * 2 + 1 >= length(args)) break
  
  arg = args[a * 2 + 1]
  
  if (arg == "--read_file") {
    file_data = args[a * 2 + 2]
  } else if (arg == "--read_table") {
    table_data = args[a * 2 + 2]
  } else if (arg == "--set_wd") {
    path = args[a * 2 + 2]
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
    pred_method[npred] = args[a * 2 + 2]
  } else if (arg == "--gene_name") {
    gene_name = args[a * 2 + 2]
  } else if (arg == "--transcript") {
    transcript = args[a * 2 + 2]
  } 
  a = a + 1
}

setwd(path)

# Input file type check
if (file_data != "" && table_data != ""){
  print("You should only give one input file. Please check.")
  stop()
} else if (file_data == "" && table_data == ""){
  print("You haven't given any input file. Please check.")
  stop()
}

# Define other functions...
# Print prediction methods if any
if (npred > 0) {
  for (a in 1:npred) {
    print(pred_method[a])
  }
}

table_generator = function(df_file_data, pred_method){
  
  # create a new table to store information needed for plotting
  # store info into the table with selected prediction methods
  default_value = c("Locus","contAB","caseAB","transcript","aa_number","aa_change")
  col_names = c(default_value, pred_method)
  filtered_table = data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(filtered_table) = col_names
  head(filtered_table)
  
  trans_link = 40
  trans_aanum = 8
  trans_aaexchange = 9
  default_trans = 7
  default_num_aachange = 3
  comments_split = strsplit(df_file_data$comment, "\\|")
  
  rr = 0
  for (r in 1:length(comments_split)) {
    sublist = comments_split[[r]]
    for (c in 1:(length(sublist) - max(trans_aanum, trans_aaexchange))) {
      if (sublist[c] == transcript) {
        aa_num = sublist[c + trans_aanum]
        aa_exchange = sublist[c + trans_aaexchange]
        if (!is.na(aa_num) && !is.na(aa_exchange) && aa_num != "" && aa_exchange != "" && nchar(aa_exchange) == default_num_aachange) {
          rr = rr + 1
          filtered_table[rr, "aa_number"] = as.numeric(aa_num)
          filtered_table[rr, "aa_change"] = aa_exchange
          filtered_table[rr, "Locus"] = sublist[1]
          filtered_table[rr, "contAB"] = df_file_data[r, "contAB"]
          filtered_table[rr, "caseAB"] = df_file_data[r, "caseAB"]
          filtered_table[rr, "transcript"] = transcript
          filtered_table[rr, "aa_change"] = aa_exchange
          for (i in 1:length(pred_method)) {
            new_col = pred_method[i]
            filtered_table[rr, new_col] = df_file_data[r, new_col]
          }
        }
      }
    }
  }
  
  if ("AM_score" %in% pred_method){
    for (i in 1: nrow(filtered_table)){
      if (filtered_table[i, "AM_score"] > 0){
        filtered_table[i, "AM_score"] = filtered_table[i, "AM_score"]/10
      }	
    }
  }
  print("Producing filtered_table based on the file data")
  return(filtered_table)
}

#define a function to generate the plots
plot_generator_1 = function(plotting_data, method){
  results = data.frame(plotting_data)
  
  png_name = sprintf("%s.png", plot_1)
  
  long_str = c(method, "aa_number")
  print(long_str)
  if(!all(long_str %in% colnames(results))) {
    stop("Not all columns in long_str exist in the results dataframe")
  }
  
  df1 = results[, long_str]
  
  df1_long = reshape2::melt(df1, id.vars = "aa_number")
  head(df1_long)
  df1_long = ddply(df1_long, .(variable), transform, ma = rollmean(value, k = 5, fill = NA))
  
  # Create the ggplot
  p1 = ggplot(df1_long, aes(aa_number, value, colour = variable)) +
    geom_point() +
    geom_smooth(span = 0.1, color = "black") + 
    facet_grid(rows = vars(variable)) +
    scale_y_continuous(limits = c(1, 10)) +
    theme_minimal()
  png(png_name, width = 6 * 300, height = 6 * 300, res = 300)
  print(p1)
  dev.off()
}

plot_generator_2  = function(plotting_data, method){
  results = data.frame(plotting_data)
  
  png_name2 = sprintf("%s.png", plot_2)
  png(png_name2, width = 6 * 300, height = 6 * 300, res = 300)
  
  long_str = c(method, "aa_number")
  
  if(!all(long_str %in% colnames(results))) {
    stop("Not all columns in long_str exist in the results dataframe")
  }
  
  r_case = 0
  df_case = data.frame(matrix(ncol = length(long_str), nrow = sum(results$caseAB > 0)))
  colnames(df_case) = c(long_str)
  for (i in 1:nrow(results)){
    if (results[i, "caseAB"] > 0){
      for (j in 1:results[i, "caseAB"]){
        r_case = r_case + 1
        df_case[r_case,] = results[i, long_str]
      }
    }
  }
  
  r_cont = 0
  df_cont = data.frame(matrix(ncol = length(long_str), nrow = sum(results$contAB > 0)))
  colnames(df_cont) = long_str
  
  for (i in 1:nrow(results)) {
    if (results[i, "contAB"] > 0) {
      for (j in 1:results[i, "contAB"]) {
        r_cont = r_cont + 1
        df_cont[r_cont, ] = results[i, long_str]
      }
    }
  }
  
  df_case_long = reshape2::melt(df_case, id.vars = "aa_number")
  df_cont_long = reshape2::melt(df_cont, id.vars = "aa_number")
  
  df_case_long$Type = "case"
  df_cont_long$Type = "cont"
  df_cc_long = rbind(df_case_long, df_cont_long)
  
  # Create the ggplot
  p2 = ggplot(df_cc_long, aes(aa_number, value, shape = Type)) +
    geom_point(aes(colour = Type), size = 2, alpha = 0.7) +
    geom_smooth(span = 0.1, aes(group = Type, colour = Type)) + 
    facet_grid(rows = vars(variable)) +
    scale_y_continuous(limits = c(1, 10)) +
    theme_minimal()
  
  print(p2)
  dev.off()
}

script_generator_pymol = function(ratio, given_table, script_name){
  # produce pymol script
  # chi_square test
  aa_num_str = vector()
  chi_value = vector()
  c_str = vector()
  
  off_set = 0
  
  df_table = given_table
  ratio = as.numeric(ratio)
  
  for (i in 1:nrow(df_table)){
    o_cont = df_table[i,"contAB"]
    o_case = df_table[i,"caseAB"]
    n_total = o_cont + o_case
    e_cont = (ratio * n_total)/(ratio + 1)
    e_case = n_total/(ratio + 1)
    chi_cont = (o_cont - e_cont)^2/e_cont
    chi_case = (o_case - e_case)^2/e_case
    chi_val = chi_cont + chi_case
    df_table[i,"chi_value"] = chi_val
    df_table[i,"cont_vs_case"] = o_cont - o_case
  }
  
  min_chi = min(df_table["chi_value"])
  max_chi = max(df_table["chi_value"])
  
  # define functions to get color
  get_color_cont = function(chi){
    paleness = 0.5 - sqrt(chi)*0.5
    color_num = paste("1,",paleness, ",",paleness)
    return(color_num)
  }
  
  get_color_case = function(chi){sqrt(chi)
    paleness = 0.5 - sqrt(chi)*0.5
    color_num = paste(paleness,",",paleness,",1")
    return(color_num)
  }
  
  # produce pymol string
  for (i in 1:nrow(df_table)){
    num = df_table[i, "aa_number"]
    chi = df_table[i,"chi_value"]/max_chi
    get_color = function(c){
      str = sprintf("set_color col%s, [%s]",i,c)
      return(str)
    }
    if (df_table[i,"cont_vs_case"] > 0){
      c = get_color_cont(chi)
      str = get_color(c)
      aa_num_str = paste(aa_num_str, sprintf("%s \n\ncolor col%s, resi %s\n",str,i, num))
    } else if (df_table[i,"cont_vs_case"] < 0){
      c = get_color_case(chi)
      str = get_color(c)
      aa_num_str = paste(aa_num_str, sprintf("%s \n\ncolor col%s, resi %s\n",str,i, num))
    }
  }
  aa_num_str_final = sprintf("color white, all\n %s", aa_num_str)
  
  script_file_name = sprintf("%s%s.txt", path, script_name)
  write.table(aa_num_str_final, script_file_name, row.names=FALSE, quote=FALSE, col.names = FALSE)
}

script_generator_pp = function(ratio, given_table, script_name){
  pred_str = vector()
  df_table = given_table
  ratio = as.numeric(ratio)
  
  # produce protein paint str
  for (i in 1:nrow(df_table)){
    if (df_table[i,"caseAB"] > 0){
      for (j in 1:df_table[i,"caseAB"]){
        pred_str = c(pred_str, sprintf("%s, %s, M",df_table[i,"aa_change"], df_table[i,"aa_number"]))
      }
    }
    if (df_table[i,"contAB"] > 0){
      for (j in 1:df_table[i,"contAB"]){
        pred_str = c(pred_str, sprintf("%s, %s, I",df_table[i,"aa_change"], df_table[i,"aa_number"]))
      }
    }
  }
  
  pred_str_clean = gsub(",", "", pred_str)
  
  script_file_name = sprintf("%s%s.txt", path, script_name)
  write.table(pred_str_clean, script_file_name, row.names=FALSE, quote=FALSE, sep = "\t")
}

if (table_data != ""){
  filtered_table = read.table(table_data, header = TRUE, stringsAsFactors = FALSE)
}

if (file_data != ""){
  df_file_data = data.frame(read.table(file_data, header = TRUE, stringsAsFactors = FALSE, fill = TRUE))
  df_whole = df_file_data
  cols_with_na = which(colSums(is.na(df_whole)) > 0)
  df_whole = na.omit(df_whole)
  df_whole[,3:cols_with_na[1]-1] = lapply(df_whole[,3:cols_with_na[1]-1], as.numeric)
  head(df_whole)
  
  filtered_table = table_generator(df_whole, pred_method)
  filtered_table[filtered_table == ""] = NA
  filtered_table = na.omit(filtered_table)
  
  summary(filtered_table)
  head(filtered_table)
}

if (w_table != "") {
  table_name = sprintf("%s%s.txt",path, w_table)
  write.table(filtered_table, table_name, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (plot_1 != "") {
  plot_table = filtered_table
  plot_generator_1(plot_table, pred_method)
}

if (plot_2 != "") {
  plot_table = filtered_table
  plot_generator_2(plot_table, pred_method)
}

if (w_script_pymol != ""){
  if (file_data != ""){
    sum_cont = sum(df_whole["contAB"])
    sum_case = sum(df_whole["caseAB"])
    ratio = sum_case/sum_cont
  }else if(table_data != ""){
    ratio = chi_ratio
  }
  df_table = filtered_table
  script_generator_pymol(ratio, df_table, w_script_pymol)
}

if (w_script_pp != ""){
  df_table = filtered_table
  script_generator_pp(ratio, df_table, w_script_pp)
}
