library(ontologyIndex)
library(hash)
library(topGO)
load('gene2BP.rda')
load('gene2MF.rda')
load('gene2CC.rda')
load('gene2MPO.rda')
load('gene2HPO.rda')
setwd('D:/work_in_vandy/precision medicine/data/')
###get latest ontology file for hpo,mpo,go###
hpo_obo<-get_OBO('https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo', propagate_relationships = "is_a", extract_tags = "minimal")
data(hpo_obo)
mpo_obo<-get_OBO('http://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology', propagate_relationships = "is_a", extract_tags = "minimal")
data(mpo_obo)
go_obo<-get_OBO('http://geneontology.org/ontology/go.obo', propagate_relationships = "is_a", extract_tags = "minimal")
data(go_obo)

save(hpo_obo,file='hpo_obo.rda')
save(mpo_obo,file='mpo_obo.rda')
save(go_obo,file='go_obo.rda')

##save hash ancestors###
go_ancestor<- go_obo[5]$ancestors
hpo_ancestor<- hpo_obo[5]$ancestors
mpo_ancestor<- mpo_obo[5]$ancestors
BP_ancestor_hash<-go_ancestor[names(inverseList(gene2BP))]
MF_ancestor_hash<-go_ancestor[names(inverseList(gene2MF))]
CC_ancestor_hash<-go_ancestor[names(inverseList(gene2CC))]
MPO_ancestor_hash<-mpo_ancestor[names(inverseList(gene2MPO))]
HPO_ancestor_hash<-hpo_ancestor[names(inverseList(gene2HPO))]

BP_ancestor_new<-BP_ancestor_hash[which(names(BP_ancestor_hash)!= 'NA')]
BP_ancestor_hash<-hash()
.set(BP_ancestor_hash,names(BP_ancestor_new),BP_ancestor_new)

CC_ancestor_new<-CC_ancestor_hash[which(names(CC_ancestor_hash)!= 'NA')]
CC_ancestor_hash<-hash()
.set(CC_ancestor_hash,names(CC_ancestor_new),CC_ancestor_new)

MF_ancestor_new<-MF_ancestor_hash[which(names(MF_ancestor_hash)!= 'NA')]
MF_ancestor_hash<-hash()
.set(MF_ancestor_hash,names(MF_ancestor_new),MF_ancestor_new)

HPO_ancestor_new<-HPO_ancestor_hash[which(names(HPO_ancestor_hash)!= 'NA')]
HPO_ancestor_hash<-hash()
.set(HPO_ancestor_hash,names(HPO_ancestor_new),HPO_ancestor_new)

MPO_ancestor_new<-MPO_ancestor_hash[which(names(MPO_ancestor_hash)!= 'NA')]
MPO_ancestor_hash<-hash()
.set(MPO_ancestor_hash,names(MPO_ancestor_new),MPO_ancestor_new)

save(BP_ancestor_hash,file ='BP_ancestor_hash.rda')
save(CC_ancestor_hash,file ='CC_ancestor_hash.rda')
save(MF_ancestor_hash,file ='MF_ancestor_hash.rda')
save(HPO_ancestor_hash,file='HPO_ancestor_hash.rda')
save(MPO_ancestor_hash,file='MPO_ancestor_hash.rda')

###calculater information content&&save###
info_cont <- function(go,root,gene2MF,go_mf_matrix) {
data(go)
MF2gene = inverseList(gene2MF) 
all_GO<-get_descendants(go, roots= root )
go_all<-length(unique(unlist(MF2gene)))
#give every annotated nodes the genes# 
gene_each<-{}
for (i in 1:length(all_GO)) {
  gene_each[i]<-MF2gene[all_GO[i]]
}
names(gene_each) <- all_GO
#give every node in ontology their and their children's annotation#
gene_each_go<-{}
for (i in 1:length(all_GO)){
  gene_each_go<- get_ancestors(go,all_GO[i])
  for (l in 1:length(gene_each_go)){
    gene_each[gene_each_go[l]]<- list(unique(c(unlist(gene_each[gene_each_go[l]][[1]]),unlist(gene_each[all_GO[i]][[1]]))))
  }
}
   IC_fun <- function(all_GO) {
     Inform_content=length(gene_each[[all_GO]])/go_all
    return(Inform_content)
   }
   pt<-sapply(all_GO,IC_fun)
   log_pt<- -logb(pt)
 log_pt_test<-log_pt[log_pt!=Inf]
   maxlog_mf<-max(log_pt_test)
   file_name = go_mf_matrix
   go_mf_matrix<- log_pt_test/maxlog_mf
   save(go_mf_matrix,file = paste0(file_name,".rda"))
}
info_cont(go_obo,root = "GO:0003674",gene2MF=gene2MF,go_mf_matrix='go_mf_matrix')
info_cont(go_obo,root = "GO:0005575",gene2MF=gene2CC,go_mf_matrix='go_cc_matrix')
info_cont(go_obo,root = "GO:0008150",gene2MF=gene2BP,go_mf_matrix='go_bp_matrix')
info_cont(hpo_obo,root = "HP:0000001",gene2MF=gene2HPO,go_mf_matrix='hpo_matrix')
info_cont(mpo_obo,root = "MP:0000001",gene2MF=gene2MPO,go_mf_matrix='mpo_matrix')

