#' @param data is list of data
#'
#' @param group data group
#' @param scale scale of data
#' @return list data
#' @export
#'
#' @examples
#'
##' \dontrun{
##' 
##' x<-readxl::read_xlsx("Col-top50.xlsx")
##' y<-readxl::read_xlsx("metadata.xlsx")
##' x%>%.[,-1]%>%t%>%as.data.frame()->x2
##' colnames(x2)<-x$Gene
##' y%>%.[,-1]%>%na.omit()%>%t%>%as.data.frame()->y2
##' row.names(y2)->rownames(x2)
##' colnames(y2)<-y%>%na.omit()%>%.$Sample%>%str_replace_all(" ","")
##' Y <-c(rep("A",4),rep("C",5))
##' X <- list(mRNA =x2,  protein = y2)
##' X %>%purrr::map(t)->X
##' .get_data_analyst_plot(X,Y)
##' }
#' @importFrom mixOmics pca splsda block.splsda spls
#' @import purrr 
#' @import dplyr

.get_data_analyst_plot <- function(data="is a list of data",group=NULL,scale = T){
.get_data_analyst <- function(data=data,group=group,scale = scale) {
  if(!(data%>%is.list()&(length(data)==2))){
    if(group%>%is.null()){cat("using pca analyst");
      return(mixOmics::pca(data%>%as.data.frame(), ncomp = 2,scale = scale))}
    if(!(group%>%is.null())){ cat("using splsda analyst")
    data<-as.data.frame(data)
    return(mixOmics::splsda(X=data,Y, ncomp=2,scale = scale))}
  }else if(length(data)==2){if(group%>%is.null()){cat("useing spls analyst")
    return(mixOmics::spls(X=data[[1]],Y=data[[2]], ncomp = 2,scale = scale))}
  }
  if(length(data)>1&is.list(data)){
  if(!(group%>%is.null())){
  {cat("using block.splsda analyst")
  return(mixOmics::block.splsda(data,group,scale = scale))}}}
}
#data%>%purrr::map(~t(.x)%>%as.data.frame)->data
.get_data_analyst(data,c(group),scale)->data_analyst_data
.get_mean<-function(data,group,scale) {
  get_med <- function(data,group,scale) {
    if(group%>%is.null())group<-rep("data",nrow(data))
    if(scale){scale(data)->data
      }
    return(data%>%data.frame(group,.)%>%dplyr::group_by(group)%>%summarise(dplyr::across(where(is.numeric), ~median(.x))))
  }
  if(!(data%>%is.list()&(length(data)==2))){get_med(data,group,scale=scale)}else{
    return(purrr::map(data,~get_med(.x,
     group =group,scale=scale)))
  }
}
.get_mean(data,group,scale)->men
return(list(analyst_data=data_analyst_data,
     men=men))
}



#explem

#' pathways_analy to pathways and analy
#'
#' @param data a list of data
#' @param group data group
#' @param org looking for https://www.genome.jp/kegg/catalog/org_list.html,using fread to keggAPI
#' @param fisher.alternative alternative  of fisher.test
#' @param scale scale of data ,It well to .get_data_analyst_plot
#' @param p.adjust.methods  looking for  p.adjust
#' @param import_model import_model=c("betweenness","degree")
#' @return metaProfiler data
#' @export
#'
#' @examples
#'
##' \dontrun{
##' 
##' x<-readxl::read_xlsx(system.file("data", "Col-top50.xlsx", package = "pathways"))
##' y<-readxl::read_xlsx(system.file("data", "metadata.xlsx", package = "pathways"))
##' CTSgetR::CTSgetR(y$Sample,"PubChem CID", "KEGG")%>%na.omit()->getkeggnames
##' full_join(getkeggnames,y,by=c("id"="Sample"),keep=T)->getkegg
##' getkegg$Sample[!is.na(getkegg$KEGG)]<-getkegg$KEGG[!is.na(getkegg$KEGG)]
##' getkegg$Sample[!is.na(getkegg$KEGG)]<-getkegg$KEGG[!is.na(getkegg$KEGG)]
##' getkegg%>%.[,-c(1:2)]%>%na.omit()->data
##' data[rev(order(data$Sample)),]->data2
##' rownames(data2)<-data2$Sample
##' data2[,-1]->data2
##' library(biomaRt)
##' mart <- useMart("ensembl","mmusculus_gene_ensembl")
##' gene_name<-getBM(attributes=c("ensembl_gene_id","entrezgene_id"),
##'                  filters = "ensembl_gene_id",values = x$Gene, mart = mart)
##' full_join(gene_name,x,by=c("ensembl_gene_id"="Gene"))->getkegg
##' getkegg$ensembl_gene_id[!is.na(getkegg$entrezgene_id)]<-getkegg$entrezgene_id[!is.na(getkegg$entrezgene_id)]
##' getkegg[,-2]->genes
##' rownames(genes)<-genes$ensembl_gene_id
##' genes[,-1]->genes
##' group<-c(rep("A",4),rep("B",5))
##' colnames(data2)<-colnames(genes)
##' list(a=data2%>%t%>%as.data.frame(),b=genes%>%t%>%as.data.frame())->data
##' pathways_analy(data = data,group = group,org = "mmu",scale = T)->n
##' 
##' # or
##' load(system.file("data", "Tes.Rdata",package = "pathways"))
##' }