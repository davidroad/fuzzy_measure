library(ontologyIndex)
####final script####
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
   go_mf_matrix<-{}
   go_mf_matrix<- log_pt_test/maxlog_mf
   save(go_mf_matrix,file = paste0(file_name,".rda"))
}
info_cont(go_obo,root = "GO:0003674",gene2MF=gene2MF,go_mf_matrix='go_mf_matrix')
info_cont(go_obo,root = "GO:0005575",gene2MF=gene2CC,go_mf_matrix='go_cc_matrix')
info_cont(go_obo,root = "GO:0008150",gene2MF=gene2BP,go_mf_matrix='go_bp_matrix')
info_cont(go_obo,root = "HP:0000001",gene2MF=gene2HPO,go_mf_matrix='hpo_matrix')
info_cont(go_obo,root = "GO:0003674",gene2MF=gene2MPO,go_mf_matrix='mpo_matrix')

