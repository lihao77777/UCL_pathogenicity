# Copyright (C) Lihao Tao and David Curtis
library(ggplot2)
library(reshape2)
library(zoo)
library(plyr)

options(bitmapType = "cairo")

args = commandArgs(trailingOnly=TRUE)

path = ""
Pred_Method = vector()
npred = 0
a = 0
saofile = ""
filtered_table_file = ""
w_table = ""
plot_1 = ""
plot_2 = ""
w_script_pymol = ""
w_script_pp = ""
gene_name = ""
Transcript = ""
chi_ratio = 0 # must be set if table is used instead of sao file - number of cases over number of controls

# Set default values if no arguments are provided
if (length(args) < 1) {
  args = c(
    "--read_saofile", "/home/lihaotao/saofiles/UKBB.HT.forAnnot.20231112.SMAD6.txt", 
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
  
  if (arg == "--read_saofile") {
    saofile = args[a * 2 + 2]
  } else if (arg == "--read_table") {
    filtered_table_file = args[a * 2 + 2]
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
    Pred_Method[npred] = args[a * 2 + 2]
  } else if (arg == "--gene_name") {
    gene_name = args[a * 2 + 2]
  } else if (arg == "--transcript") {
    Transcript = args[a * 2 + 2]
  } 
  a = a + 1
}

setwd(path)

# Input file type check
if (saofile != "" && filtered_table_file != ""){
  print("You should only give one input file. Please check.")
  stop()
} else if (saofile == "" && filtered_table_file == ""){
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

table_generator = function(df_file_data, Pred_Method){
  
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
  return(filtered_table_final)
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
  
  df_filtered = given_table
  ratio = as.numeric(ratio)
  for (i in 1:nrow(df_filtered)){
    o_cont = df_filtered[i,"contAB"]
    o_case = df_filtered[i,"caseAB"]
    n_total = o_cont + o_case
    e_cont = (ratio * n_total)/(ratio + 1)
    e_case = n_total/(ratio + 1)
    chi_cont = (o_cont - e_cont)^2/e_cont
    chi_case = (o_case - e_case)^2/e_case
    chi_val = chi_cont + chi_case
    df_filtered[i,"chi_value"] = chi_val
    df_filtered[i,"cont_vs_case"] = o_cont - o_case
  }
  
  min_chi = min(df_filtered["chi_value"])
  max_chi = max(df_filtered["chi_value"])
  
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
  for (i in 1:nrow(df_filtered)){
    num = df_filtered[i, "aa_number"]
    chi = df_filtered[i,"chi_value"]/max_chi
    get_color = function(c){
      str = sprintf("set_color col%s, [%s]",i,c)
      return(str)
    }
    if (df_filtered[i,"cont_vs_case"] > 0){
      c = get_color_cont(chi)
      str = get_color(c)
      aa_num_str = paste(aa_num_str, sprintf("%s \n\ncolor col%s, resi %s\n",str,i, num))
    } else if (df_filtered[i,"cont_vs_case"] < 0){
      c = get_color_case(chi)
      str = get_color(c)
      aa_num_str = paste(aa_num_str, sprintf("%s \n\ncolor col%s, resi %s\n",str,i, num))
    }
  }
  aa_num_str_final = sprintf("color white, all\n %s", aa_num_str)
  
  aa_num_str_name = sprintf("%snew_aa_num_str_%s.txt", path, gene_name)
  script_file_name = sprintf("%s%s.txt", path, script_name)
  write.table(aa_num_str_final, script_file_name, row.names=FALSE, quote=FALSE, col.names = FALSE)
}

script_generator_pp = function(ratio, given_table, script_name){
  pred_str = vector()
  df_filtered = given_table
  ratio = as.numeric(ratio)
  
  # produce protein paint str
  for (i in 1:nrow(df_filtered)){
    if (df_filtered[i,"caseAB"] > 0){
      for (j in 1:df_filtered[i,"caseAB"]){
        pred_str = c(pred_str, sprintf("%s, %s, M",df_filtered[i,"aa_change"], df_filtered[i,"aa_number"]))
      }
    }
    if (df_filtered[i,"contAB"] > 0){
      for (j in 1:df_filtered[i,"contAB"]){
        pred_str = c(pred_str, sprintf("%s, %s, I",df_filtered[i,"aa_change"], df_filtered[i,"aa_number"]))
      }
    }
  }
  
  pred_str_clean = gsub(",", "", pred_str)
  
  script_file_name = sprintf("%s%s.txt", path, script_name)
  write.table(pred_str_clean, script_file_name, row.names=FALSE, quote=FALSE, sep = "\t")
}

if (filtered_table_file != ""){
  filtered_table_final = read.table(filtered_table_file, header = TRUE, stringsAsFactors = FALSE)
}

if (saofile != ""){
  df_file_data = data.frame(read.table(saofile, header = TRUE, stringsAsFactors = FALSE, fill = TRUE))
  df_whole = df_file_data
  cols_with_na = which(colSums(is.na(df_whole)) > 0)
  df_whole = na.omit(df_whole)
  df_whole[,3:cols_with_na[1]-1] = lapply(df_whole[,3:cols_with_na[1]-1], as.numeric)
  head(df_whole)
  
  filtered_table_final = table_generator(df_whole, Pred_Method)
  filtered_table_final[filtered_table_final == ""] = NA
  filtered_table_final = na.omit(filtered_table_final)
  
  summary(filtered_table_final)
  head(filtered_table_final)
}

if (w_table != "") {
  table_name = sprintf("%s.txt",w_table)
  write.table(filtered_table_final, table_name, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (plot_1 != "") {
  plot_table = filtered_table_final
  plot_generator_1(plot_table, Pred_Method)
}

if (plot_2 != "") {
  plot_table = filtered_table_final
  plot_generator_2(plot_table, Pred_Method)
}

if (w_script_pymol != ""){
  df_filtered = filtered_table_final
  if (saofile != ""){
    sum_cont = sum(df_whole["contAB"])
    sum_case = sum(df_whole["caseAB"])
    ratio = sum_case/sum_cont
  }else if(filtered_table_file != ""){
    ratio = chi_ratio
  }
  df_filtered = filtered_table_final
  script_generator_pymol(ratio, df_filtered, w_script_pymol)
}

if (w_script_pp != ""){
  df_filtered = filtered_table_final
  script_generator_pp(ratio, df_filtered, w_script_pp)
}
