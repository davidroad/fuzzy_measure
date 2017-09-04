####User manual for prioritizing candidates genes from NGS pipeline filtration with symptoms represented genes####
#please user R >=3.2.2
#please install R library ##ontologyIndex##hash##TopGO##
#please set or add your R library to ' export lib='/USER/R/x86_64-pc-linux-gnu-library/3.2' '
#####configure need to be set in sever_all.R#####
#the work environment path was set to './data'
#There are 4 needed arguments to be set, and 1 optional argument#
#seed gene entrez id#seed_gene= read.table("../entrez_id_final_25.txt",header = TRUE) or absolute path
#cand gene entrez id#cand_gene= read.table("../candidate_entrez.txt") or absolute path
#output path#output="../result"
#output filename# is a string like "RR"##
#cutoff number# if not set,default is 0.05##

##after the configuration###
##command example###
#R CMD BATCH "--args "/scratch/cqs/udn/code/fuzzy_measure/entrez_id_final_05.txt" "/scratch/cqs/udn/code/fuzzy_measure/candidate_entrez.txt" "/scratch/cqs/udn/code/fuzzy_measure/result" "RR" "0.05"" main.R
##find the file named like 'RRfinal_table.txt' to view the output folder /output path like./result ###

