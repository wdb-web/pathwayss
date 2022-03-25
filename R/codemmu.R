#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @importMethodsFrom clusterProfiler enrichResult compareClusterResult
#' @import purrr
#' @import dplyr
#' @import data.table
#' @importFrom stringr str_split str_locate_all str_sub str_remove str_match
#' @importFrom methods setClass news 
#' @importFrom data.table  fread
#' @importFrom ggplot2 `+` ggplot aes geom_point theme_bw
 
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
  
{#import_model=c("betweenness","degree")
  .get_import <- function(enrichKEGG,import_model="degree",kegg_org) {
    
    kegg_org%>%as.data.frame%>%.[,colnames(kegg_org)%in% c("path","Compound_ID")]->base
    base  %>%unique( )->base
    base%>%dplyr::filter( path %in%enrichKEGG@result$Description )->base
    as.data.table(base)->base
    base[ , if (.N > 1) 
      transpose(combn(  unique(path), 2L, simplify = FALSE)), by = Compound_ID
    ][ , .(Sum = .N), keyby = .(Group.1 = V1, Group.2 = V2)]->l
    colnames(l)<-c("from","to","weight")
    net_pc<-igraph::graph_from_data_frame(
      d=l, directed=F)
    
    g <- igraph::set_edge_attr(net_pc, "weight", value= l$weight)
    
    if(import_model[1]=="betweenness"){
      igraph::betweenness(net_pc)->import_min}
    if(import_model[1]=='degree'){
      igraph::degree(net_pc, mode='total', loops=FALSE)->import_min}
    import_min%>%as.data.frame()%>%mutate(name = rownames(import_min))->import_min
    import_min%>%as.data.frame()%>%mutate(name = rownames(import_min))->import_min
    
    enrichKEGG@result  ->d2
    d2%>%dplyr::group_by(Description )%>%dplyr::summarise(geneID %>% str_split("/")%>%unlist())%>%setNames(c("path", "Compound_ID" ))%>%
      ungroup()%>%as.data.table()->da
    da[ , if (.N > 1) 
      transpose(combn(path, 2L, simplify = FALSE)), by = Compound_ID
    ][ , .(Sum = .N), keyby = .(Group.1 = V1, Group.2 = V2)]->l 
    colnames(l)<-c("from","to","weight")
    
    net_pc<-igraph::graph_from_data_frame(
      d=l,
      directed=F)
    g <- igraph::set_edge_attr(net_pc, "weight", value= l$weight)
    if(import_model[1]=="betweenness"){
      igraph::betweenness(net_pc)%>%as.data.frame()->import}
    if(import_model[1]=='degree'){
      igraph::degree(net_pc, mode='total', loops=FALSE)->import}
    import%>%as.data.frame()%>%mutate(name = rownames(import))->import
    import%>%as.data.frame()%>%mutate(name = rownames(import))->import
    
    left_join(import,import_min,by="name")%>%
      {.[,1]/.[,3]}%>%setNames(left_join(import,import_min,by="name")%>%.$name)%>%as.data.frame()->import
    colnames(import)<-"import"
    import$Description <-import%>%rownames()
    enrichKEGG@keytype<-"kegg"
    enrichKEGG@result[,1:9]%>%full_join(import,by = "Description")->enrichKEGG@result
    enrichKEGG@result$import[is.na(enrichKEGG@result$import)]<-0
    return(enrichKEGG)
  }
  get_all_p <- function(g){
    g%>%dplyr::group_by(Cluster)%>%dplyr::summarise(w=BgRatio%>%str_split("/")%>%unlist()%>%.[2])%>%ungroup()->n
    g%>%dplyr::select(Cluster,ID ,pvalue)%>%dplyr::group_split(Cluster)%>%purrr::map(~.x%>%.[,-1])%>%purrr::reduce(full_join,by="ID")->k
    k[is.na(k)]<-1
    .get_p_combine(k[1,-1],method = 'fisher',w =n[,2])$p_value
    map_dfr(k[,-1]%>%t%>%as.data.frame,~.get_p_combine(.x,method = 'fisher',w =n[,2])$p_value)%>%t%>%as.data.frame()->p
    rownames(p)<-k$ID
    p.adjust(p$V1,method = "BH")->p.adjust
    p.adjust(p$V1,method = "fdr")->qvalue
    data.frame( pvalue=p,p.adjust=p.adjust,qvalue=qvalue)->p.data
    p.data$ID<-k$ID
    g%>%dplyr::group_by(ID,Description)%>%dplyr::summarise(geneID=geneID%>%
                                                             paste(sep = "/",collapse = "/"),Count=sum(Count))->data

    colnames(g)


    full_join(data,p.data)->data
    colnames(data)[5]<-"pvalue"
    return(data)
  }
  .get_p_combine <- function(p, method = c("fisher", "SL", "MG", "tippett"), w = NULL) {
    method <- match.arg(method, choices = c("fisher", "SL", "MG", "tippett"))
    p <- p[!is.na(p)]
    n <- length(p)
    if (max(p) - .Machine$double.ep > 1) {
      stop("invalid input, p > 1")
    }

    if (n == 0) {
      warning("vector of p-values is empty")
      if (method == "fisher") {
        return(list(statistic = NA, p_value = NA, method = "Fisher (1932)",
                    statistic_name = "Xsq"))
      } else if (method == "SL") {
        return(list(
          statistic = NA, p_value = NA,
          method = "Stouffer (1949), Liptak (1958)",
          statistic_name = "Z"
        ))
      } else if (method == "MG") {
        return(list(
          statistic = NA, p_value = NA,
          method = "Mudholkar and George (1979)",
          statistic_name = "L"
        ))
      } else if (method == "tippett") {
        return(list(
          statistic = NA, p_value = NA,
          method = "Tippett (1931)",
          statistic_name = "p_min"
        ))
      }
    } else {
      if (method == "fisher") {
        Xsq <- -2 * sum(log(p))
        p_val <- stats::pchisq(Xsq, df = 2 * n, lower.tail = FALSE)
        return(list(
          statistic = Xsq, p_value = p_val, method = "Fisher (1932)",
          statistic_name = "Xsq"
        ))
      } else if (method == "SL") {
        if (is.null(w)) {
          w <- rep(1, n) / n
        } else {
          if (length(w) != n) {
            stop("length of p and w must be equal")
          }
        }
        Zi <- stats::qnorm(1 - p)
        Z <- sum(w * Zi) / sqrt(sum(w^2))
        p_value <- 1 - stats::pnorm(Z)
        return(list(
          statistic = Z, p_value = p_value,
          method = "Stouffer (1949), Liptak (1958)",
          statistic_name = "Z"
        ))
      } else if (method == "MG") {
        L <- sum((-1) * log(p / (1 - p)))
        p_value <- stats::pt(L * sqrt((15 * n + 12) / (pi^2 * n * (5 * n + 2))),
                             df = 5 * n + 4, lower.tail = FALSE
        )
        return(list(
          statistic = L, p_value = p_value,
          method = "Mudholkar and George (1979)",
          statistic_name = "L"
        ))
      } else if (method == "tippett") {
        return(list(
          statistic = min(p), p_value = 1 - (1 - min(p))^n,
          method = "Tippett (1931)",
          statistic_name = "p_min"
        ))
      }
    }
  }
  get_org<-function(org=org) {
    ls<-ls(envir=.GlobalEnv)
    get_class <- function(.x) {
      tryCatch(
        return(eval(parse(text = .x))%>%class%>%{"metaProfiler"%in%.})  ,
        error = function(e){return(F)
        },finally = {
        })
    }
    ls[purrr::map(ls,~get_class(.x))%>%unlist()]->l

    get_org_1 <- function(l) {
      l%>%parse(text = .)%>%eval%>%.@org_organism%>%return()
    }
    if(length(l)>1){purrr::map(l,get_org_1)%>%unlist->k
      if(org%in%k){
        return(parse(text = l[k==org]%>%.[1])%>%eval%>%.@org) }}
    if(length(l)==1){
      l%>%parse(text = .)%>%eval%>%.@org_organism->k
      if(org%in%k){
        return(parse(text = l[k==org])%>%eval%>%.@org)
      }
    }
    return(.get_kegg_org(org))

  }

  .get_kegg_org<- function(org="hsa") {
    fread(paste("http://rest.kegg.jp/list/pathway/",org,sep = ""),header = F)->k
    name_rmends <- function(.x) {
      .x %>%str_locate_all(" - ")%>%.[[1]]->jk
      jk[nrow(jk),]%>%.[1]%>%{.-1}%>%str_sub(.x, end =.)%>%return
    }# getname
    k$V2%>%purrr::map(~name_rmends(.x))%>%unlist()%>%cbind(k[,1])->pathways# getname

    fread("http://rest.kegg.jp/list/pathway/map",header = F)->enknowpathways
    colnames(pathways)[1]<-"path"
    colnames(enknowpathways)[2]<-"path"
    dplyr::left_join(pathways,enknowpathways,by=colnames(enknowpathways)[2])->k2
    colnames(k2)[2:3]<-c("ko","pathway")
    k2$pathway%>%str_remove("path:")%>%.[1]
    fread(paste("http://rest.kegg.jp/link/cpd/",
                "map00010",sep = ""),header = F,verbose=F,showProgress=F)
    pb <- dplyr::progress_estimated(length(k2$pathway))
    get_kegg_comp <- function(.x) {
      pb$tick()$print()

      tryCatch({
        return(fread(paste("http://rest.kegg.jp/link/cpd/",
                           .x,sep = ""),header = F,verbose=F,showProgress=F))
      },
      error = function(e){return(NULL)
      },finally = {
      })
    }

    k2$pathway%>%str_remove("path:")%>%purrr::map(~get_kegg_comp(.x)
    )->k
    k%>%purrr::reduce(bind_rows)->kegg
    colnames(kegg)<-c("pathway","Compound_ID")
    kegg%>%as.data.table()%>%group_by(pathway)%>%
      dplyr::summarise(Compound_ID=Compound_ID%>%str_remove("cpd:")%>%paste(sep="/",collapse = "/"),
                       n=n())%>%right_join(k2,.)%>%.[,-c(2:3)]%>%na.omit()->pathway

    inner_join(k2,kegg)->b

    b %>% dplyr::mutate(Compound_ID=Compound_ID%>%str_remove("cpd:"),
                        ko=ko%>%str_remove("path:"),
                        pathway=pathway%>%str_remove("path:")
    )%>%return
  }

  get_anong<-function(org=org) {
    ls<-ls(envir=.GlobalEnv)
    get_class <- function(.x) {
      tryCatch(
        return(eval(parse(text = .x))%>%class%>%{"metaProfiler"%in%.})  ,
        error = function(e){return(F)
        },finally = {
        })
    }
    ls[purrr::map(ls,~get_class(.x))%>%unlist()]->l

    get_org_1 <- function(l) {
      l%>%parse(text = .)%>%eval%>%.@org_organism%>%return()
    }
    if(length(l)>1){purrr::map(l,get_org_1)%>%unlist->k
      if(org%in%k){
        return(parse(text = l[k==org]%>%.[1])%>%eval%>%.@anong) }}
    if(length(l)==1){
      l%>%parse(text = .)%>%eval%>%.@org_organism->k
      if(org%in%k){
        return(parse(text = l[k==org])%>%eval%>%.@anong)
      }
    }
    return(NULL)
  }

  .get_kegg_p<- function(value=data,kegg_pathways=kegg_pathways,p_model=p_model,
                         p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods) {
    value%>%unique()->data
    kegg_pathways%>%as.data.frame()->match.df
    kegg_pathways%>%group_by(Compound_ID)%>%summarise(n()-1)->import
    kegg_pathways%>%group_by(path)%>%
      summarise(n=sum( data%in%Compound_ID),N=n(),
                m=sum( data%in%kegg_pathways$Compound_ID),M=kegg_pathways$Compound_ID%>%unique%>%length,
                Compound_ID=paste(unique(Compound_ID[Compound_ID%in%value]),
                                  collapse = "/"))%>%ungroup()%>%
      dplyr::filter(n>0)->n
    fisher<-function(data,model=c("phyper","fisher")){
      gene.not.interest=data.frame(c(d[3]-d[1]),c(d[4]-d[2]-d[3]+d[1])) %>% t;
      gene.in.interest=data.frame(d[1],d[3]) %>% t
      d <- data.frame(gene.not.interest,gene.in.interest)
      if(model[1]=="fisher"){return(fisher.test(d, alternative ="greater" )[["p.value"]])}
      # k -》 富集到通路，n 总的差异基因，N 所有基因，M 这条通路基因
      data -> d
      if(model[1]=="phyper"){return(phyper(d[1]-1,d[3], d[4]-d[3], d[2], lower.tail=F))}
     cat("using model is ","phyper or","fisher")
    }
    n[2,2:5]->d
    apply(n[,2:5],1, fisher,model=p_model)->n$pvalue
    p.adjust(n$pvalue,method =p.adjust.methods)->n$p.adjust
    p.adjust(n$pvalue,method =qvalue.methods)->n$qvalue
    kegg_pathways$path[kegg_pathways$Compound_ID%in% import$Compound_ID[import$Compound_ID%in%value]]%>%
      {kegg_pathways$Compound_ID[kegg_pathways$path%in%.]}->every_diff_pathways_com
    kegg_pathways%>%group_by(path)%>%
      summarise(import=sum(import$`n() - 1`[import$Compound_ID%in%Compound_ID[Compound_ID%in%value]])/
                  sum(import$`n() - 1`[import$Compound_ID%in%Compound_ID]))%>%
      ungroup()->n2
    full_join(n,n2,by = "path")%>%na.omit()->data
    right_join(kegg_pathways[,-4]%>%unique,data)->data
    data$GeneRatio<-paste(data$n,data$m,sep = "/");data$BgRatio<-paste(data$N,data$M,sep = "/")
    data%>%dplyr::select(ID=ko,Description=path,GeneRatio,BgRatio,
                         pvalue,p.adjust,qvalue,geneID=Compound_ID,Count=n,import)%>%as.data.frame()->data
    rownames(data)<-data$ID
    data%>%return

  }


  #name<-data
  .get_kegg_analyst<-function(name="list Compound_ID or gean",
                              org="hsa",p_model=c("phyper","fisher"),
                              kegg_pathways=kegg_pathways,p.adjust.methods=p.adjust.methods,
                              qvalue.methods="fdr",enrichKEGG=kk,import_model=import_model,everyorg=everyorg
  ){
    if(name[1]%>%str_match("C[0-9][0-9][0-9][0-9][0-9]")%>%{!is.na(.)}){
      .get_kegg_p(value=name%>%unique(),kegg_pathways=kegg_pathways,
                  p_model=p_model,
                  p.adjust.methods=p.adjust.methods,qvalue.methods=qvalue.methods)->y
      rownames(y)<-y$ID
      enrichKEGG@result<-y
      enrichKEGG@gene<-name
      enrichKEGG@universe<-kegg_pathways$Compound_ID%>%unique()
      kegg_pathways[kegg_pathways$ko%in%kegg_pathways$ko[kegg_pathways$Compound_ID%in%name],]->kegg_pathways
      kegg_pathways[,c(2,4)]%>%dplyr::group_split(ko)->j

      purrr::map(j,.f =  function(.x) {
        return(.x$Compound_ID%>%c)
      })->set
      names(set)<-kegg_pathways$ko%>%unique()
      enrichKEGG@geneSets<-set
      enrichKEGG@organism<-org
      enrichKEGG@pvalueCutoff<-1
      enrichKEGG@qvalueCutoff<-1
      enrichKEGG@pAdjustMethod<-p.adjust.methods
      enrichKEGG@universe<-""
      .get_import(enrichKEGG,import_model=import_model,kegg_org  = everyorg)->enrichKEGG
      enrichKEGG@result[order(enrichKEGG@result$pvalue),]->enrichKEGG@result
      return(enrichKEGG)
    }
    if(name[1]%>%str_match("C[0-9][0-9][0-9][0-9][0-9]*")%>%{is.na(.)}){
      clusterProfiler::enrichKEGG(gene         = name,
                                  pAdjustMethod =p.adjust.methods,
                                  organism     = org,pvalueCutoff = 1,qvalueCutoff = 1
      )%>%return
    }

  }

  kegg_pathway1 <- function(data=c("list"),org="mmu",p_model=c("phyper","fisher"),
                            p.adjust.methods="holm",import_model=c("degree","betweenness")) {


    get_org(org = org ) ->kegg_pathways
    kegg_pathways[kegg_pathways$path !="Metabolic pathways",] -> kegg_pathways
    s <- new("metaProfiler",org_organism =org)
    s@org<-kegg_pathways%>%as.data.frame()
    # get_anong(org)%>%.[[1]]->bb
    # bb->>bb
    #s@universe<-kegg_pathways$Compound_ID%>%unique()

    kk <- new("enrichResult",organism=org,ontology="KEGG")
    if(is.list(data)&(length(data)==1)){
      .get_kegg_analyst(name = data[[1]],org = org,kegg_pathways=kegg_pathways,import_model=import_model,
                        p_model=p_model,p.adjust.methods=p.adjust.methods,everyorg = kegg_pathways,
                        enrichKEGG=kk)%>%list->s@kegg_analyst
      names(s@kegg_analyst)<-"enrichKEGG"
      return(s)
    }
    if(!is.list(data)){.get_kegg_analyst(name = data,org = org,kegg_pathways=kegg_pathways,
                                         p_model=p_model,everyorg = kegg_pathways,
                                         p.adjust.methods=p.adjust.methods,import_model=import_model,
                                         enrichKEGG=kk)%>%list->s@kegg_analyst
      names(s@kegg_analyst)<-"enrichKEGG"
      return(s)
    }
    #data[[1]]->name
    data%>%
      purrr::map(~.get_kegg_analyst(name = .x,org = org,kegg_pathways=kegg_pathways,import_model=import_model,enrichKEGG=kk,everyorg = kegg_pathways,
                                    p_model=p_model,p.adjust.methods=p.adjust.methods))->y


    x <- new("compareClusterResult")


    Result <- function(.x,.y) {
      if(ncol(.x@result)==9){.x@result$import<-""}
      data.frame(Cluster=.y,.x@result)
    }
    purrr::map2(y, names(y),~Result(.x,.y))%>%purrr::reduce(rbind)->g
    x@compareClusterResult<-g

    gene <- function(.x) {
      .x@gene
    }
    purrr::map(y,~gene(.x))->x@geneClusters
    factor(x@compareClusterResult$Cluster)->b

    x@compareClusterResult$Cluster<-factor(x@compareClusterResult$Cluster )

    s@kegg_analyst<-list(enrichKEGG=y,compareClusterResult=x ,every=get_all_p(g=x@compareClusterResult)   )
    return(s)
  }

  comp_dotplot<- function(data) {
    data@result%>%ggplot(aes(y=-log10(qvalue),
                             x=import,col=import,size=-log10(qvalue),label=Description))+
      geom_point()+ ggsci::scale_color_material("orange")+
      ggrepel::geom_text_repel()+theme_bw()+ggplot2::scale_size(range = c(4,8) )+
      ggplot2::ylim(0,3)
  }
}
