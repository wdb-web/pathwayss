#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @importMethodsFrom clusterProfiler enrichResult compareClusterResult
#' @import purrr
#' @import dplyr
#' @importFrom stringr str_split str_locate_all str_sub str_remove str_match
#' @importFrom methods setClass news 
setClass("metaProfiler",
         slots = c(
           kegg_analyst = "list",
           data_mixOmics_analyst = "list",
           data_mean="list",
           org="data.frame",
           org_organism="character",
           data="list",
           group="data.frame"
         ),
         prototype = list(
           kegg_analyst = list(),
           data_mixOmics_analyst = list(),
           data_mean=list(),
           org=data.frame(),
           org_organism=character(),
           data=list(),
           group=data.frame()
         )
)
#data%>%purrr::map(colnames)%>%.[[2]]->data
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
pathways_analy <- function(data,group,org="hsa",p_model=c("phyper","fisher"),
                           scale=T,
                 p.adjust.methods="holm",import_model=c("betweenness","degree")) {
  if(is.list(data)){data%>%purrr::map(colnames)%>%
    kegg_pathway1(data = . ,org = org,p.adjust.methods=p.adjust.methods,p_model = p_model,
                  import_model=import_model)->y}else{
    kegg_pathway1(data =colnames(data) ,org = org,p.adjust.methods=p.adjust.methods,p_model = p_model,
                       import_model=import_model)->y}
  .get_data_analyst_plot( data=data,group = group,scale = scale)->h
  slot(y, "data_mixOmics_analyst") <- list(h$analyst_data)
  slot(y, "data_mean") <-h$men
  slot(y, "data") <-data
  slot(y, "group") <-group%>%as.data.frame()
  
  return(y)
}

#' plot_chor is plot metaProfiler data 
#'
#' @param f using dplyr to filter data
#' @param metaProfiler  metaProfiler data
#' @param pathways_scale_color  pathways bar color
#' @param directional https://jokergoo.github.io/circlize_book/book/
#' @param annotationTrack https://jokergoo.github.io/circlize_book/book/
#' @param line_scale_col  line color
#' @param ... https://jokergoo.github.io/circlize_book/book/
#'
#' @return a list of length ,you can using packLegend and draw to add length
#' @export
#'
#' @examples
##' \dontrun{
##' load(system.file("data", "Tes.Rdata",package = "pathways"))
##'  n@kegg_analyst$compareClusterResult%>%clusterProfiler::filter(
##' Description%in%c((n@kegg_analyst$compareClusterResult%>%group_by(Description)%>%
                    ##'                     summarise(n=n()>1))%>%dplyr::filter(n==T)%>%.$Description)|qvalue     <0.05
##' )%>%as.data.frame()->f
##' circos.clear()
##' circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE,start.degree=90)
##' plot_chor(f,metaProfiler = n)->v
##' #looking for ComplexHeatmap[https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html]
##' pd = packLegend(list =  v)
##' draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
##' }




plot_chor<- function(f=f,metaProfiler,   scale_height=0.6, begin=0,pathways_scale_color=structure(c(-log(0.05),-log(0.01)),
          names =c("red","blue")),directional = 0,
                     annotationTrack = NULL,line_scale_col=NULL,...){
  f[,c("ID","geneID","Cluster","Count","pvalue","GeneRatio")]%>%
    dplyr::group_by(ID,Cluster,Count,pvalue,GeneRatio)%>%
    dplyr::summarise(geneID=geneID%>%
                       str_split("/")%>%unlist())->f
  
  f$Cluster%>%paste(f$geneID,sep = "_")->f$line
  #chord(f[,c(7,1)])
  f[,c(7,1)]%>%unlist()%>%c->b
  group2<- structure(c( f[,c(2,7)]%>%unique%>%.$Cluster ,rep("pathways",(!b%in%(f[,c(2,7)]%>%unique%>%.$line))%>%sum)),
                     names = c(f[,c(2,7)]%>%unique%>%.$line,b[!b%in%(f[,c(2,7)]%>%unique%>%.$line)] ))
  
  chord(data=f[,c(7,1)],directional =directional,group=group2,target.prop.height = mm_h(1),diffHeight = mm_h(1),
        preAllocateTracks =1,annotationTrack=annotationTrack,col = line_scale_col,small.gap=0,
        annotationTrackHeight = c(0.01, 0.01),...)->length_chord
  #f$Cluster<-factor(f$Cluster)
  #scale_height<-0.6
  f->k
 
#  plot_add(f = k,begin=parent.frame(2)$begin,scale_height=parent.frame(2)$scale_height,
#          metaProfiler=parent.frame(2)$metaProfiler,group2=group2)->length_col
  
  
  
  
  
  
  {
    metaProfiler->n
    scale_height<-scale_height
    grt <- function(data,begin,scale_height=0.4,scale_color=structure(c(-log(0.05),-log(0.01)),names =c("red","blue"))) {
      circos.lines(c(0,length( get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])])+0.5),
                   y = c( (begin-1)*scale_height, (begin-1)*scale_height),
                   sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.lines(c(0,length( get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])])+0.5),
                   y = c( (begin-1)*scale_height, (begin-1)*scale_height)+0.5*scale_height*0.95,
                   sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.lines(c(0,length( get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])])+0.5),
                   y = c( (begin-1)*scale_height, (begin-1)*scale_height)+1*scale_height*0.95,
                   sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.text(-0.5,track.index = 1,y=(begin-1)*scale_height+scale_height/2,labels =data$Cluster[1] ,
                  sector.index = get.all.sector.index()[get.all.sector.index()%in%names(group2[group2=="pathways"])][1]
      )
      circos.text(x=rep(1.5,3),track.index = 1,y=c((begin-1)*scale_height+scale_height/2,(begin-1)*scale_height+scale_height*0.8,
                                                   (begin-1)*scale_height+0.05 ),labels = c(max(data$Count)/2,max(data$Count),0),
                  sector.index = get.all.sector.index()[get.all.sector.index()%in% names(group2[group2=="pathways"])][sum(get.all.sector.index()%in%names(group2[group2=="pathways"]))]
      )
      add.do(x= data$ID ,h = data$Count,scale.col = scale_color,
             col= -log(data$pvalue), begin = (begin-1)*scale_height ,scale_height =scale_height*0.95 )%>%return()
      
    }
    
    f%>%dplyr::group_by(Cluster) %>% dplyr::group_split()%>%setNames(unique(f$Cluster))%>%
      purrr::map2(seq_along(unique(f$Cluster)),~grt(data =  .x,begin = .y))%>%.[[1]] -> length_pathways
    get.vip<-function(.x,.y) {
      tryCatch(
        .x$value%>%as.data.frame()->.x,
        error = function(e){
          return(NULL)
        },
        finally = {
        }
        
      )
      if(.y!="comp"){
        rownames(.x)<-paste(.y,"_",rownames(.x),sep = "");return(.x) }
      return(NULL)
    }
    mixOmics::selectVar(n@data_mixOmics_analyst[[1]])%>%
      purrr::map2((names(mixOmics::selectVar(n@data_mixOmics_analyst[[1]]))),~get.vip(.x,.y))->c
    c%>%purrr::reduce(rbind)->data
    
    get_add.names <- function(.x,.y) {
      colnames(.x)%>%paste(.y,.,sep = "_")->colnames(.x)
      .x[colnames(.x)%in%get.all.sector.index()]->.x
      return(.x)
    }
    #n@data%>%purrr::map2(n@data%>%names,get_add.names)%>%purrr::reduce(cbind)->data
    #scale_height=0.6
    .get_theme <- function(x,f=f,scale_height=0.6,begin=0,data=data) {
      #x="a"
      scale_height<-scale_height
      
      f[f$Cluster==x,]->k
      circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%(k$line)]%>%length()+0.5),
                   y = c((begin)*scale_height+scale_height, (begin)*scale_height+scale_height),
                   sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%(k$line)]%>%length()+0.5),
                   y = c( 0, 0),
                   sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%(k$line)]%>%length()+0.5),
                   y = c(((begin)*scale_height+scale_height)*abs(min(data[rownames(data)%in%
                                                                            get.all.sector.index(),]))/
                           (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                              max(max(data[rownames(data)%in%get.all.sector.index(),]))),
                         ((begin)*scale_height+scale_height)*abs(min(data[rownames(data)%in%get.all.sector.index(),]))/
                           (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                              max(max(data[rownames(data)%in%get.all.sector.index(),])))),
                   sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1],
                   col = "grey",track.index = 1,lty=2
      )
      
      circos.text(-1,track.index = 1,y=c((begin)*scale_height+scale_height)*0.9,labels =
                    max(data[rownames(data)%in%get.all.sector.index(),])%>%round(.,2) ,
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
      )
      
      circos.text(-1,track.index = 1,y=c((begin)*scale_height+scale_height)*0,labels =
                    min(data[rownames(data)%in%get.all.sector.index(),])%>%round(.,2) ,
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
      )
      
      circos.text(-1,track.index = 1,y = c(((begin)*scale_height+scale_height)*max(data[rownames(data)%in%get.all.sector.index(),])/
                                             (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                                                max(max(data[rownames(data)%in%get.all.sector.index(),]))),
                                           ((begin)*scale_height+scale_height)*max(data[rownames(data)%in%get.all.sector.index(),])/
                                             (abs(min(data[rownames(data)%in%get.all.sector.index(),]))+
                                                max(max(data[rownames(data)%in%get.all.sector.index(),])))),labels =0 ,
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
      )
      
      circos.text(-0.5,track.index = 1,y=-0.1,labels =
                    k$Cluster[1] ,
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%(k$line)][1]
      )
    }
    purrr::map(f$Cluster%>%unique,~.get_theme(x=.x,f=f,scale_height=scale_height,begin=0,data=data))
    #f$Cluster%>%unique%>%.[1]->.x
    #scale_height=0.6
    add.do(x=rownames(data)[rownames(data)%in%get.all.sector.index()],h=data[rownames(data)%in%get.all.sector.index(),],
           begin = scale_height/2,
           scale_height =scale_height,#scale.col = scale.col,
           col =data[rownames(data)%in%get.all.sector.index(),]
    )->length_hub
    
    
    n@data_mean%>%names%>%purrr::map(function(.x){n@data_mean%>%.[[.x]]->x;
      colnames(x)%>%str_remove("X")->colnames(x)
      colnames(x)[-1]<-paste(.x,colnames(x)[-1],sep = "_");return(x%>%as.data.frame())})%>%
      purrr::reduce(cbind)->mean
    mean[,c("group",get.all.sector.index()[get.all.sector.index()%in%f$line])]->mean
    mean->mean2
    (mean[,-1]-min(mean[,-1]))/max(mean[,-1]-min(mean[,-1]))->mean[,-1]
    
    #chord(f[,c(7,1)],directional = 0,group=group2,annotationTrack = NULL,
    #  annotationTrackHeight = c(0.01, 0.01))
    
    for (i in seq_along(table(f$Cluster))) {
      f[f$Cluster%in%(names(table(f$Cluster))[i]),]->l
      
      circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%l$line]%>%length()+0.5),
                   y = c( scale_height,scale_height)+0.05,
                   sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1],
                   col = "grey",track.index = 1,lty=2
      )
      
      circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%l$line]%>%length()+0.5),
                   y = c( scale_height+max(mean[,-1])*scale_height/2,scale_height+max(mean[,-1])*scale_height/2)+0.05,
                   sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.lines(c(0,get.all.sector.index()[get.all.sector.index()%in%l$line]%>%length()+0.5),
                   y = c( scale_height+max(mean[,-1])*scale_height/2/2,scale_height+max(mean[,-1])*scale_height/2/2)+0.05,
                   sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1],
                   col = "grey",track.index = 1,lty=2
      )
      circos.text(-0.5,track.index = 1,y=scale_height+max(mean[,-1])*scale_height/2,labels =round(max(mean2[,-1]),2) ,
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1]
      )
      circos.text(-0.5,track.index = 1,y=scale_height+0.05,labels =round(min(mean2[,-1]),2),
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1]
      )
      circos.text(-0.5,track.index = 1,y=scale_height+max(mean[,-1])*scale_height/2/2+0.05,labels =
                    round(max(mean2[,-1])/2+min(mean2[,-1])/2,2),
                  sector.index =get.all.sector.index()[get.all.sector.index()%in%l$line][1]
      )
      col=c(ggsci::pal_aaas(palette = "default",alpha = 0.9)(10))
      pch<-16:25
      for (t in seq_len(mean%>%nrow())){
        {
          circos.lines(x=seq_along(get.all.sector.index()[get.all.sector.index()%in%l$line])-0.5,
                       y=c(mean[t,get.all.sector.index()[get.all.sector.index()%in%l$line]]/2/2+0.6)%>%as.numeric()+0.05,
                       col=col[t],track.index = 1,pt.col = col[t],
                       lwd=3,cex =1.5,pch = pch[t],
                       type="o",
                       sector.index=get.all.sector.index()[get.all.sector.index()%in%l$line][1])
        }
      }
      
      
    }
    length_mean =Legend(labels = c( table(f$Cluster)%>%names()), title = "length_mean", type = "points", pch = pch[seq_along(table(f$Cluster)%>%names())],
                        legend_gp = gpar(col =(c(ggsci::pal_aaas(palette = "default",
                                                                 alpha = 0.9)(10))[seq_along(table(f$Cluster)%>%names())] )))
    for (i in get.all.sector.index()) {
      if(i%in%f$line){add.name(begin =scale_height+max(mean[,-1])*scale_height/2/2+0.05,
                               name =f$geneID[f$line==i]%>%unique(),x = i)
        
      }else{add.name(begin =scale_height+max(mean[,-1])*scale_height/2/2+0.05,
                     name =i,x = i)}
      
      
    }
    #return(list(length_mean,length_pathways,length_hub))
    list(length_mean,length_pathways,length_hub) ->length_col
  }
  
  
  
  
  
  
  
  
  
  
  return(list(chord=length_chord,mean=length_col[[1]],pathways=length_col[[2]],hub=length_col[[3]]))
  #col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  #colorRampPalette(c("white", col_vector[8]))(20)[c((seq_along(1:20)+1)%%2)%>%as.logical %>%which%>%.[-c(1)]]->g
}

