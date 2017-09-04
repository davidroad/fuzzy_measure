args=(commandArgs(TRUE))
seed_gene_file<-args[[1]]
cand_gene_file<-args[[2]]
output_location<-args[[3]]
file_name<-args[[4]]
if (length(args)==5 && is.numeric(args[[5]])) {
cutoff<-args[[5]]
}else{
cutoff = 0.05}
setwd('./data')
seed_gene= read.table(seed_gene_file,header = TRUE)
cand_gene= read.table(cand_gene_file)
seed_gene = as.character(unlist(seed_gene))
cand_gene = as.character(unlist(cand_gene))
output=output_location

library(ontologyIndex)
library(hash)

load("gene2GO.rda")
load("gene2BP.rda")
load("gene2MF.rda")
load("gene2CC.rda")
load("gene2pubmed.rda")
load("gene2pathway.rda")
load("BPALL2gene.rda")
load("MFALL2gene.rda")
load("CCALL2gene.rda")
load("HPOALL2gene.rda")
load("MPOALL2gene.rda")	
load("GO_Ontology.rda")
load("go_cc_matrix.rda")
load("go_bp_matrix.rda")
load("go_mf_matrix.rda")
load("hpo_matrix.rda")
load("mpo_matrix.rda")
load("pubmedmatrix.rda")
load("pathwaymatrix.rda")
load('pathway2gene.rda')
load('pubmed2gene.rda')
load('gene2HPO.rda')
load('gene2MPO.rda')
load('BP_ancestor_hash.rda')
load('CC_ancestor_hash.rda')
load('MF_ancestor_hash.rda')
load('MPO_ancestor_hash.rda')
load('HPO_ancestor_hash.rda')
load('HPO2gene.rda')
load('entrezid.rda')
data(go)
data(mpo)
data(hpo)

###Ontology enrichment###
enrichGOFun <- function(seed_gene,gene2GO,go,GOALL2gene) {
#terms annotated by seed gene set, each presented multiple times indicative of number of seed genes being annotated
	seed_set_profile<-{}
	GO_seed_gene_list = unlist(gene2GO[seed_gene])
	GO_seed_gene_list <- unique(GO_seed_gene_list)
	GO_seed_gene_list<-get_ancestors(go,GO_seed_gene_list)
    	seed_length_fun <- function(GO_seed_gene_list,seed_gene,GOALL2gene){
	seed_length<- intersect(seed_gene,unique(GOALL2gene[[GO_seed_gene_list]]))
	}
	seed_length <- sapply(GO_seed_gene_list,seed_length_fun,seed_gene = seed_gene,GOALL2gene = GOALL2gene)
	seed_select <- seed_length[seed_length!='NULL']
	k = length(unique(unlist(seed_select)))
	total_length = length(unique(unlist(GOALL2gene)))
	gene_length = sapply(GOALL2gene[names(seed_select)],length)
	each_length = sapply(seed_select[names(seed_select)],length)
	if (length(each_length)>0) {
		pvalue_list=phyper(q=each_length-1,m=gene_length,n=total_length-gene_length,k,lower.tail=F) 	
		pvalue_list_adjusted<- p.adjust(pvalue_list,method = "BH")
		pvalue_list_adjusted<-pvalue_list_adjusted[order(pvalue_list_adjusted)]
		seed_set_profile<- pvalue_list_adjusted[which(pvalue_list_adjusted<cutoff)]
	}
	return(seed_set_profile)
}

###pubmed enrichment###
enrichFun <- function(seed_gene,gene2pubmed,pubmed2gene) {
#terms annotated by seed gene set, each presented multiple times indicative of number of seed genes being annotated
	seed_set_profile<-{}
	pubmed_seed_gene_list = unlist(gene2pubmed[seed_gene])
	seed_length<- length(pubmed_seed_gene_list)
	seed_gene_list_set<- unique(pubmed_seed_gene_list)
	each_length =table(pubmed_seed_gene_list)
	gene_length = sapply(pubmed2gene[names(each_length)],length)
	total_length = length(unlist(pubmed2gene))
	if (length(each_length)>0) {	
		pvalue_list=phyper(q=each_length-1,m=gene_length,n=total_length-gene_length,k=seed_length,lower.tail=F) 
#		names(pvalue_list)<- names(each_length)
		pvalue_list_adjusted<- p.adjust(pvalue_list,method = "BH")
		pvalue_list_adjusted<-pvalue_list_adjusted[order(pvalue_list_adjusted)]
		seed_set_profile<- pvalue_list_adjusted[which(pvalue_list_adjusted<cutoff)]
	}
	return(seed_set_profile)
}
####fuzzy measure lambda####
f<-function(lambda,vec){prod(1+lambda*vec)-lambda-1}
###call ancestor###
fastinter<-function(x,y)
{
y[match(x,y)]
}

common_ancestor_test<-function(cand_set,seed_set_profile,profile){
common_ancestor_each<-fastinter(profile[[seed_set_profile]],profile[[cand_set]])
common_ancestor_each<-common_ancestor_each[length(common_ancestor_each)]
return(common_ancestor_each)
}

ancestor_all_fun<-function(seed_set_profile,cand_set,profile){
common_ancestor<-sapply(cand_set,common_ancestor_test,seed_set_profile=seed_set_profile,profile =profile)
return( unique(common_ancestor))
}

####new GO mainfunction###
####GO mainfunction###

GO_mainfunction<- function(data_set,seed_set_profile,gene2GO,GOmatrix,profile) {
	s_seed<- {}
	s_cand<- {}
	score_seed_set = {}
	score_cand_set = {}
#######
	cand_set = {}
	inter = {}
	if (is.element(data_set,names(gene2GO))) {
		cand_set <- gene2GO[[data_set]]
#######
#######

		common_ancestor_all <-{}
		common_ancestor_all<-sapply(seed_set_profile,ancestor_all_fun,cand_set=cand_set,profile =profile)
		common_ancestor_group<-unique(as.vector(unlist(common_ancestor_all)))
		seed_set_1<-unique(c(common_ancestor_group,seed_set_profile))
		cand_set_1<-unique(c(common_ancestor_group,cand_set))
		seed_set<-intersect(seed_set_1,names(GOmatrix))
		cand_set<-intersect(cand_set_1,names(GOmatrix))
		inter<-intersect(seed_set,cand_set)
		if (length(inter) >0) {
			seed_set <- c(inter,setdiff(seed_set,inter))
			cand_set <- c(inter,setdiff(cand_set,inter))
		
			score_seed_set<-GOmatrix[intersect(seed_set,names(GOmatrix))]
			score_cand_set<-GOmatrix[intersect(cand_set,names(GOmatrix))]
			result_seed<-{}
			if (length(score_seed_set)>1) {
				if (sum(score_seed_set)<1){
					result_seed<-uniroot(f,interval=c(0.00001,10000),vec=score_seed_set)
					lambda_seed<-result_seed$root
					}
				if (sum(score_seed_set)>1) {
					result_seed<-uniroot(f,interval=c(-1,-0.00001),vec=score_seed_set)
					lambda_seed<-result_seed$root
					}
			} else
			{
				score<- score_seed_set[1]
				return (score)
			}
			result_cand<-{}
			if (length(score_cand_set)>1) {
				if (sum(score_cand_set)<1) {
					result_cand<-uniroot(f,interval=c(0.00001,10000),vec=score_cand_set)
					lambda_cand<-result_cand$root
					}
				if (sum(score_cand_set)>1) {
					result_cand<-uniroot(f,interval=c(-1,-0.00001),vec=score_cand_set)
					lambda_cand<-result_cand$root
					}
			}else
			{
				score<- score_cand_set[1]
				return (score)
			}
			if (length(inter)>1) {
				score_inter_cand<-score_cand_set[1:length(inter)]
				s_cand =1
				for (t in 1:(length(score_inter_cand)-1)){
					s_cand = score_inter_cand[t]+score_inter_cand[t+1]+lambda_cand*score_inter_cand[t]*score_inter_cand[t+1]
					score_inter_cand[t+1] = s_cand
				}
				score_inter_seed<-score_seed_set[1:length(inter)]
				s_seed =1
				for (s in 1:(length(score_inter_seed)-1)){
					s_seed= score_inter_seed[s]+score_inter_seed[s+1]+lambda_seed*score_inter_seed[s]*score_inter_seed[s+1]
					score_inter_seed[s+1] = s_seed
				}
				score<- (s_seed+s_cand)/2 
			}else 
			{
				score<- score_seed_set[1]
			}
		}
		else {
			score<- 0
		}
	}else {
		score<- NA
	}
#	print(data_set)
	return (score)
}

###pubmed main function####
mainfunction<- function(cand_set,seed_set_profile,gene2pubmed,pubmedmatrix) {
	score_seed_set = {}
	score_cand_set = {}
#	cand_set = {}
	inter = {}
	if (length(cand_set)>1) {
#		cand_set <- gene2pubmed[data_set]
#		print (data_set)
		seed_set<-seed_set_profile
		inter<-intersect(seed_set,cand_set)
		if (length(inter) >0) {
			seed_set <- c(inter,setdiff(seed_set,inter))
			cand_set <- c(inter,setdiff(cand_set,inter))
			score_seed_set<-pubmedmatrix[seed_set]
			score_cand_set<-pubmedmatrix[cand_set]
			result_seed<-{}
			if (length(score_seed_set)>1) {
				if (sum(score_seed_set)<1){
					result_seed<-uniroot(f,interval=c(0.00001,10000),vec=score_seed_set)
					lambda_seed<-result_seed$root
					}
				if (sum(score_seed_set)>1) {
					result_seed<-uniroot(f,interval=c(-1,-0.00001),vec=score_seed_set)
					lambda_seed<-result_seed$root
					}
			} else
			{
				score<- score_seed_set[1]
				return (score)
			}
			result_cand<-{}
			if (length(score_cand_set)>1) {
				if (sum(score_cand_set)<1) {
					result_cand<-uniroot(f,interval=c(0.00001,10000),vec=score_cand_set)
					lambda_cand<-result_cand$root
					}
				if (sum(score_cand_set)>1) {
					result_cand<-uniroot(f,interval=c(-1,-0.00001),vec=score_cand_set)
					lambda_cand<-result_cand$root
					}
			}else
			{
				score<- score_cand_set[1]
				return (score)
			}
			if (length(inter)>1) {
				score_inter_cand<-score_cand_set[1:length(inter)]
				s_cand =1
				for (t in 1:(length(score_inter_cand)-1)){
					s_cand = score_inter_cand[t]+score_inter_cand[t+1]+lambda_cand*score_inter_cand[t]*score_inter_cand[t+1]
					score_inter_cand[t+1] = s_cand
				}
				score_inter_seed<-score_seed_set[1:length(inter)]
				s_seed =1
				for (s in 1:(length(score_inter_seed)-1)){
					s_seed= score_inter_seed[s]+score_inter_seed[s+1]+lambda_seed*score_inter_seed[s]*score_inter_seed[s+1]
					score_inter_seed[s+1] = s_seed
				}
				score<- (s_seed+s_cand)/2 
			}else 
			{
				score<- score_seed_set[1]
			}
		}
		else {
			score<- 0
		}
	}else {
		score<- NA
	}
#	print(data_set)
	return (score)
}

randfunction<- function(seed_gene,repeattime){
	rand_gene<- setdiff(entrezid,seed_gene)
	rand_cand_gene <- sample(rand_gene,repeattime,replace=FALSE) # UPDATE: replace set to F instead of T.
	return(rand_cand_gene) # NOTE: a vector of egIDs.
}
###change input for pubmed###
data_prepare_fun<- function(x,gene2pubmed)gene2pubmed[x]

pval_from_rand <- function(x,rand_scores) (sum(rand_scores[is.finite(rand_scores)]>=x)+1)/length(rand_scores)
overall_table<- {}
overall_table<- data.frame(cand_gene)
repeattime<- 500
rand_gene <- randfunction(seed_gene,repeattime)
######GO CC####


fuzzy_input_fun<-function(seed_gene,gene2GO,go,GOALL2gene,GOmatrix,profile,cand_gene,rand_gene,p_overall,overall_table){
seed_set_profile<- enrichGOFun(seed_gene=seed_gene,gene2GO=gene2GO,go = go, GOALL2gene= GOALL2gene)
seed_set_profile<- names(seed_set_profile)
cand_score<- sapply(cand_gene,GO_mainfunction,seed_set_profile = seed_set_profile,gene2GO=gene2GO,GOmatrix = GOmatrix,profile = profile)
rand_scores <- sapply(rand_gene,GO_mainfunction,seed_set_profile = seed_set_profile,gene2GO=gene2GO,GOmatrix = GOmatrix,profile = profile)
cand_pval <- sapply(as.vector(cand_score),pval_from_rand,as.vector(rand_scores))
p_overall<<- cbind(p_overall,as.vector(cand_pval))
overall_table<<- data.frame(overall_table,as.vector(cand_score),as.vector(cand_pval))
return(cbind(as.vector(cand_score),as.vector(cand_pval)))
}

pubmed_input_fun<-function(cand_gene,seed_gene,rand_gene,gene2pubmed,pubmed2gene,pubmedmatrix,p_overall,overall_table){
cand_set = data_prepare_fun(cand_gene,gene2pubmed=gene2pubmed)
rand_set = data_prepare_fun(rand_gene,gene2pubmed=gene2pubmed)
seed_set_profile<- enrichFun(seed_gene=seed_gene,gene2pubmed = gene2pubmed,pubmed2gene=pubmed2gene)
seed_set_profile<- names(seed_set_profile)
cand_score<- sapply(cand_set,mainfunction,seed_set_profile = seed_set_profile,gene2pubmed=gene2pubmed,pubmedmatrix = pubmedmatrix)
rand_scores <- sapply(rand_set,mainfunction,seed_set_profile = seed_set_profile,gene2pubmed=gene2pubmed,pubmedmatrix = pubmedmatrix)
cand_pval <- sapply(as.vector(cand_score),pval_from_rand,as.vector(rand_scores))
p_overall<<- cbind(p_overall,as.vector(cand_pval))
overall_table<<- data.frame(overall_table,as.vector(cand_score),as.vector(cand_pval))
return(cbind(as.vector(cand_score),as.vector(cand_pval)))
}
###command###
p_overall <- {}
overall_table<- {}
overall_table<- data.frame(cand_gene)
cand_pval_CC<-fuzzy_input_fun(seed_gene=seed_gene,rand_gene=rand_gene,cand_gene=cand_gene,gene2GO=gene2CC,go=go,GOALL2gene=CCALL2gene,GOmatrix=go_cc_matrix,profile=CC_ancestor_hash,p_overall=p_overall,overall_table=overall_table)
cand_pval_MF<-fuzzy_input_fun(seed_gene=seed_gene,rand_gene=rand_gene,cand_gene=cand_gene,gene2GO=gene2MF,go=go,GOALL2gene=MFALL2gene,GOmatrix=go_mf_matrix,profile=MF_ancestor_hash,p_overall=p_overall,overall_table=overall_table)
cand_pval_BP<-fuzzy_input_fun(seed_gene=seed_gene,rand_gene=rand_gene,cand_gene=cand_gene,gene2GO=gene2BP,go=go,GOALL2gene=BPALL2gene,GOmatrix=go_bp_matrix,profile=BP_ancestor_hash,p_overall=p_overall,overall_table=overall_table)
cand_pval_HPO<-fuzzy_input_fun(seed_gene=seed_gene,rand_gene=rand_gene,cand_gene=cand_gene,gene2GO=gene2HPO,go=hpo,GOALL2gene=HPOALL2gene,GOmatrix=hpo_matrix,profile=HPO_ancestor_hash,p_overall=p_overall,overall_table=overall_table)
cand_pval_MPO<-fuzzy_input_fun(seed_gene=seed_gene,rand_gene=rand_gene,cand_gene=cand_gene,gene2GO=gene2MPO,go=mpo,GOALL2gene=MPOALL2gene,GOmatrix=mpo_matrix,profile=MPO_ancestor_hash,p_overall=p_overall,overall_table=overall_table)
cand_pval_pubmed<-pubmed_input_fun(seed_gene=seed_gene,cand_gene=cand_gene,rand_gene=rand_gene,gene2pubmed=gene2pubmed,pubmedmatrix = pubmedmatrix,pubmed2gene=pubmed2gene,p_overall=p_overall,overall_table=overall_table)
cand_pval_pathway<-pubmed_input_fun(seed_gene=seed_gene,cand_gene=cand_gene,rand_gene=rand_gene,gene2pubmed=gene2pathway,pubmedmatrix = pathwaymatrix,pubmed2gene=pathway2gene,p_overall=p_overall,overall_table=overall_table)

###calculate overall_p_value###
pchi<-{}
overall_p_value<-{}
for (l in 1:length(p_overall[,1])){
sum = 1
p <- p_overall[l,]
p<-na.omit(as.vector(unlist(p_overall[l,])))

if (length(p)>0) {
for (i in 1:length(p)) {
sum = p[i]*sum
}
sum = as.numeric(-2*log(sum))
pchi[l]<- (1-pchisq(sum,df = 2*length(p)))
}else{
pchi[l]<-1}

}
p_overall<-cbind(p_overall,pchi)

overall_table<-cbind(overall_table,pchi)
overall_table_final<-overall_table[order(pchi),]
names(overall_table_final)<-c('Cand_gene','GO_CC','GO_CC_pvalue','GO_MF','GO_MF_pvalue','GO_BP','GO_BP_pvalue','HPO','HPO_pvalue','MPO','MPO_pvalue','Pubmed','Pubmed_pvalue','Pathway','Pathway_pvalue','Overall_pvalue')
setwd(output)
write.table(overall_table_final,file = paste(file_name,"final_table.txt",sep=''),sep = '\t',quote=FALSE,row.names=F)
