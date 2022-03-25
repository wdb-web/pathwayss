#' packagename: Package containing tools support project work.
#'
#' This package provides a central repository for methods to facilitate 
#' project work
#' @import circlize
#' @import dplyr
#' @importFrom ComplexHeatmap Legend
#' @import purrr
#' @importFrom randomcoloR distinctColorPalette 
getcol<-function (col = col, scale.col = scale.col)
{
  if (is.numeric(col)) {
    if (is.na(scale.col)) {
      scale.col <- structure(c(max(col), median(col),
                               min(col)), names = c("red", "white", "blue"))
      return(a = data.frame(group = c(scale.col), col = c(scale.col %>%
                                                            names)))
    }
    else if (!(scale.col %>% names %>% is.numeric())) {
      return(data.frame(group = c(scale.col), col = c(scale.col %>%
                                                        names)))
    }
    else {
      scale.col <- structure(c(max(col), median(col),
                               min(col)), names = c("red", "white", "blue"))
      return(a = data.frame(group = c(scale.col), col = c(scale.col %>%
                                                            names)))
    }
  }
  if (!is.numeric(col)) {
    if (is.na(scale.col[1])) {
      col_chr = data.frame(group = c(col %>% unique()),
                           col = c(distinctColorPalette(length(col %>%
                                                                 unique()))))
    }else {
      if(!length(col[!col %in% names(scale.col)] %>% unique())){
        col_chr = data.frame(group = c(col[!col %in% names(scale.col)] %>%
                                         unique(), names(scale.col)), col = c( scale.col ))
      }else{
        col_chr = data.frame(group = c(col[!col %in% names(scale.col)] %>%
                                         unique(), names(scale.col)), col = c(distinctColorPalette(length(col[!col %in%
                                                                                                                names(scale.col)] %>% unique())), scale.col ))}
      return(col_chr)
    }
  }
}


do <- function(x,y,col=col,newcol=newcol,begin=begin,plotdata=plotdata) {
  c(0,plotdata$..h)->hi
  if(col%>%is.numeric){
    colorRamp2(c(newcol$group), c(newcol$col))->plotcol
    for (i in seq_along(plotdata$..h)) {
      circos.rect(0,ybottom =begin+hi[i] ,ytop =begin+hi[i+1],
                  xright = get.cell.meta.data("xlim",sector.index = y, track.index = 1)[2],
                  sector.index = y, track.index = 1, border =NA,
                  col=ifelse(is.na(x[i]),"white",plotcol(x[i])))

    }
  }else{
    for (i in seq_along(plotdata$..h)) {
      circos.rect(0,ybottom =begin+hi[i] ,ytop =begin+hi[i+1],
                  xright = get.cell.meta.data("xlim",sector.index = y, track.index = 1)[2],
                  sector.index = y, track.index = 1, border =NA,
                  col=ifelse(is.na(x[i]),"white",newcol$col[newcol$group==x[i]]))

    }
  }

}

do_bar <- function(x=x,col=col, newcol=newcol, begin=begin,plotdata=plotdata){
  #colnames(newdata)[2]
  plotdata<-as.data.frame(plotdata)
  c(0, plotdata[,colnames(plotdata)==x]%>%as.numeric)->hi
  for (i in seq_along(plotdata$col)) {
    circos.rect(0,ybottom =begin+hi[i] ,ytop =begin+hi[i+1],
                xright = get.cell.meta.data("xlim",sector.index = x, track.index = 1)[2],
                sector.index = x, track.index = 1, border =NA,
                col=newcol$col[newcol$group==plotdata$col[i]])

  }
}
plotposition <- function(position="fill",newdata=newdata,scale_height=0.5) {
  if(position=="fill"){apply(newdata[,-1],2,function(.x) {
    cumsum(.x)/sum(.x)})}->newdata[,-1]
  if(position=="stack"){
    max(apply(newdata[,-1],2,sum))->max
    apply(newdata[,-1],2,function(.x){
      cumsum(.x)})->newdata[,-1]
    newdata[,-1]/max->newdata[,-1]
  }
  newdata[,-1]*scale_height->newdata[,-1]

  return(newdata)
}

add.name <- function(x=x,name=name,begin=0.8 ,col = "black", sector.index = NA,
                     font = 1, cex = 1.2,facing =  "clockwise",
                     adj=c(0, degree(0)), niceFacing = TRUE,...) {
  x%>%as.data.frame()%>%unlist()%>% c->x
  name%>%as.data.frame()%>%unlist()%>% c->name
  for (i in seq_along(x)) {
    circos.text(x=get.cell.meta.data("xlim",sector.index = x[i], track.index = 1)[2]/2,
                y=begin, name[i],  col =col, sector.index = x[i],
                font = font, cex = cex,facing = facing,
                adj =adj, niceFacing = niceFacing,...)
  }
}
add.pch <- function(x=x,h=h,pch=pch,begin=0.8,cex=2,
                    scale_height = 0.8,scale.pch=NA,.order=NA,...
) {
  data.frame(x=x,h,pch=pch)->data
  if (((!x%in%get.all.sector.index())%>%sum)) {    #是不是原来全部数据都有
    cat("you begin data is", get.all.sector.index())
    cat("but,now you data is" ,x)
    cat("diff data is ",x[!x%in%get.all.sector.index()])
    stop()
  }
  colnames(data)<-c("x","..h","pch")
  data%>%tidyr::pivot_wider(names_from = x,values_from=pch)->newdata
  if(!is.na(.order)){newdata$..h<-factor(newdata$..h,levels = .order)
  newdata[order(newdata$..h),]->newdata}#排序
  newdata->plotdata
  (1:nrow(plotdata))*(scale_height/nrow(plotdata))->plotdata$..h
  #e=3
  if(is.na(scale.pch)){
    structure(c(data$pch%>%unique), names = c(seq_along(data$pch%>%unique))+15)->scale.pch
  }
  c(0,plotdata$..h)->hi
  i=1
  for(e in seq_along(newdata)[-1]){
    for(i in 1:nrow(newdata)){
      circos.points(x=mean(c(  get.cell.meta.data("xlim",sector.index = colnames(plotdata)[e]))),
                    y=mean(c(hi[c(i,i+1)]))+begin,sector.index=colnames(plotdata)[e],
                    pch =names(scale.pch)[which(c(scale.pch%>%as.character())==c(plotdata[i,e]))]%>%as.numeric(),cex =cex ,...

      )

    }
  }

}

#中心图型




add.do <- function(x=x,h=h,col=NA,  .order=NA,begin=0,position="fill",
                   scale_height=1,scale.col=NA) {
  
  data.frame(x=x,h,col=col)->data
  data[x%in%get.all.sector.index(),]->data

  getcol(col=col,scale.col=scale.col)->newcol # num col is.fac
  colnames(data)<-c("x","..h","col")
  if(!is.numeric(data$..h)){#heatmap
    data%>%tidyr::pivot_wider(names_from = x,values_from=col)->newdata
    if(!is.na(.order)){newdata$..h<-factor(newdata$..h,levels = .order)
    newdata[order(newdata$..h),]->newdata}
    newdata->plotdata
    (1:nrow(plotdata))*(scale_height/nrow(plotdata))->plotdata$..h
    #e=3
    for(e in seq_along(newdata)[-1]){
      do(x=newdata[,e]%>%as.matrix()%>%c(),y=colnames(newdata)[e],col=col,newcol=newcol,begin=begin,plotdata=plotdata)
    }
    if(!is.numeric(newcol$group))  return( Legend(at = c(newcol$group),
                                                  legend_gp = gpar(col =newcol$col,fill =newcol$col),
                                                  title_position = "topleft",
                                                  title = "legend"))

    return( Legend(col_fun = colorRamp2(c(newcol$group), c(newcol$col)), title_position = "topleft",
                   title = "legend"))  #
  } #
  if(is.numeric(data$..h)){ #
    if(is.numeric(data$col)){
      #(data$..h)/max(data$..h-min(data$..h))->data$..h
      if(min(data$..h)<0)      (data$..h)/max(data$..h-min(data$..h))->data$..h
      if(min(data$..h)>0)    (data$..h)/max(abs(data$..h))->data$..h
      colorRamp2(c(newcol$group)%>%as.numeric ,c(newcol$col))->plotcol
      for(i in 1:nrow(data)){
        circos.rect(0,ybottom =begin ,ytop =begin+data$..h[i]*scale_height,
                    xright = get.cell.meta.data("xlim",sector.index = data$x[i], track.index = 1)[2],
                    sector.index = data$x[i], track.index = 1, border =NA,
                    col=plotcol(data$col[i]))
      }
      return(Legend(col_fun =plotcol, title_position = "topleft",
                    title = "legend"))
    }


    data%>%tidyr::pivot_wider(names_from = x,values_from = ..h)->newdata
    if(!is.na(.order)){newdata$col<-factor(newdata$col,levels = .order)
    newdata[order(newdata$col),]->newdata}#鎺掑???
    newdata[is.na(newdata)]<-0
    newdata->plotdata
    #plotposition(position,plotdata,scale_height)->plotdata
    #e=2
    for(e in seq_along(newdata)[-1]){
      do_bar(x=colnames(newdata)[e],
             col=col,
             newcol=newcol,
             begin=begin,
             plotdata=plotdata)
    }
    return( Legend(at = c(newcol$group),
                   legend_gp = gpar(col =newcol$col,fill =newcol$col), title_position = "topleft",
                   title = "legend"))  #杩斿洖鍥句緥
  }
}


chord<-function (data = data, col = NA, directional = F, annotationTrack = NULL,
                 big.gap = 15, small.gap = 0.1, scale = T, preAllocateTracks = list(track.height = 0.3),
                 ...)
{
  data <- as.data.frame(data)
  df <- data
  df
  data <- data[(data[, 1] %in% names(col) | data[, 2] %in%
                  names(col)) %>% {. == F }, ]
  palette <- distinctColorPalette(length(table(data[, 1])))
  grid.col = structure(c(palette, rep("grey", length(table(data[,2])))),
                       names = c(as.character(as.data.frame(table(data[,  1]))[, 1]),
                                 as.character(as.data.frame(table(data[, 2]))[, 1]) ))
  chordDiagram(df, grid.col = grid.col, directional = directional,
               annotationTrack = annotationTrack, big.gap = big.gap,
               small.gap = small.gap,
               scale = scale,
               preAllocateTracks = preAllocateTracks,...)->data2
  names((grid.col[names(grid.col)%in%get.all.sector.index()]))->gecol
  add.do(x= (grid.col[names(grid.col)%in%get.all.sector.index()])%>%names(),
         h = rep("a",length(get.all.sector.index())),
         col =(grid.col[names(grid.col)%in%get.all.sector.index()])%>%names(),
         scale.col=grid.col,
         scale_height = 0.1,begin = -0.1)


  return(Legend(at = c(df[, 1] %>% unique), type = "lines",
                legend_gp = gpar(col = grid.col[names(grid.col) %in%
                                                  (df[, 1] %>% unique)], fill = grid.col[names(grid.col) %in%
                                                                                           (df[, 1] %>% unique)], width = 3, lwd = 10), title_position = "topleft",
                title = "legend"))
}

easy.clusterProfiler<-function (keggline, bar.col = NA, line.col = NA, pathway.name = T, 
                                gean.name = T, bar.height.scale = 1, name.begin = 6, ...) 
{
  as.num <- function(.x) {
    eval(parse(text = .x))
  }
  keggline$GeneRatio <- purrr::map_dfc(keggline$GeneRatio, 
                                       as.num) %>% as.numeric()
  library(stringr)
  if (is.na(bar.col)) {
    bar.col <- structure(c(0, 5), names = c("yellow ", 
                                            "red"))
  }
  keggline.ID <- data.frame(str_split_fixed(keggline$geneID, 
                                            pattern = "/", n = Inf), keggline$ID)
  newdata <- keggline.ID %>% tidyr::pivot_longer(cols = -keggline.ID) %>% 
    .[, -2]
  if (is.na(line.col)) 
    a <- chord(newdata[newdata$value != "", c(2, 1)], 
               ...)
  if (!is.na(line.col)) 
    a <- chord(newdata[newdata$value != "", c(2, 1)], 
               col = line.col, ...)
  d <- add.do(x = keggline$ID, h = keggline$GeneRatio, col = -log10(keggline$pvalue), 
              scale.col = bar.col, begin = 0, scale_height = bar.height.scale)
  if (pathway.name) 
    add.name(x = newdata[newdata$value != "", c(1)], 
             name = newdata[newdata$value != "", c(1)], 
             begin = name.begin)
  pas <- function(.x, .y) {
    paste(rep(" ", time = .y), .x, sep = "")
  }
  namw <- map2_chr(newdata[newdata$value != "", c(2)] %>% 
                     as.data.frame() %>% c %>% .[[1]], {
                       str_length(newdata[newdata$value != "", c(2)] %>% 
                                    as.data.frame() %>% c %>% .[[1]]) %>% max - str_length(newdata[newdata$value != 
                                                                                                     "", c(2)] %>% as.data.frame() %>% c %>% .[[1]])
                     }, pas)
  if (gean.name) 
    add.name(x = newdata[newdata$value != "", c(2)], 
             name = namw, begin = name.begin)
  return(list(a, d))
}






# char_one_clusterProfiler ------------------------------------------------
#' Title
#'
#' @param metaProfiler a metaProfiler data
#' @param data data is who
#' @param show show number
#' @param pvalue  pvalue =0.05 
#' @param heatmap_scale_color  heatmap color
#' @param bar.height.scale  pathways bar height
#' @param name.begin begin names higth 
#' @param pathway.name if T will add pathway names
#' @param gean.name if T will add gean names
#' @param bar.col   pathways bar color 
#' @param scale  scale of width
#' @param ... 
#'
#' @return list heatmap
#' @export
#'
#' @examples
#' char_one_clusterProfiler(n)%>%packLegend(list = .)->lgd_list_vertical 
#' draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4,"mm"), just = c("left", "bottom"))






char_one_clusterProfiler <- function(metaProfiler,data=NA,show=125,pvalue=0.05,heatmap_scale_color=NA,
                                     bar.height.scale=0.6,name.begin =0.8,pathway.name = T,
                                     gean.name = F,bar.col =NA,scale =T,...) {
  metaProfiler->n
  if(is.na(data)){data<-"[[1]]"}
  data->a
  eval(parse(text =paste("n@kegg_analyst$enrichKEGG",a,"@result",sep = "")))->k
  k[(k$pvalue < pvalue),]->keggline
  if(nrow(keggline)>show){keggline[1:show,]->keggline}
  easy.clusterProfiler(keggline = keggline
                       ,name.begin =0.8,pathway.name = pathway.name,gean.name = gean.name ,...
  )->length
  eval(parse(text =paste("n@data",a,sep = "")))%>%t%>%as.data.frame()->heatmap
  heatmap%>%rownames()->heatmap$ENTREZID
  heatmap[heatmap$ENTREZID%in%get.all.sector.index(),]->heatmap
  heatmap%>%
    tidyr::pivot_longer(cols =  -ENTREZID)%>%dplyr::filter(ENTREZID%in%get.all.sector.index())->heatmap_go
  add.do(x=heatmap_go$ENTREZID,h=heatmap_go$name,scale_height = 0.8,begin = 0,
         col = heatmap_go$value)->heatmap2
  
  return(list(length[[1]],length[[2]],heatmap2))
}


