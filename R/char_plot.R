#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @import circlize
#' @import dplyr
#' @importFrom ComplexHeatmap Legend
#' @import purrr
#' @importFrom randomcoloR distinctColorPalette 



plot_add <- function(f,begin=begin,scale_height=scale_height,scale_color=scale_color,
                     metaProfiler=metaProfiler,group2=group2) {
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
  return(list(length_mean,length_pathways,length_hub))
  
}