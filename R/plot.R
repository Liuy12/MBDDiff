## whole genome GC distribution by quantile of methylation levels. 
MethyEnrich <- function(bin_count, fa, interactive = F){
  GCcon <- CountFreqency(fa, CG =T, ATCG = F)$CG
  temp <- quantile(bin_count)
  labels <- switch(length(unique(temp)),
                   a = c('0-100%'),
                   b = c('0-75%', '75-100%'),
                   c = c('0-50%', '50-75%', '75-100%'),
                   d = c('0-25%', '25-50%', '50-75%', '75-100%'))
  cuts <- cut(bin_count, breaks = unique(temp), labels = labels)
  cols <- c("black", "blue", "purple", "gray")[1:length(unique(temp))]
  temp1 <- bind_cols(temp, as.data.frame(cuts))
  colnames(temp1) <- c('RPKM', 'CG', 'cuts')
  ggplot() + 
    geom_density(aes(x = CG, color = cuts), adjust = 2, data = temp1)+
    scale_colour_manual(values = cols)
  if(interactive){
    dataMat <- filter(temp1, cuts == '0-25%')
    denStat <-
      density(dataMat$CG, from = min(temp1$CG), to = max(temp1$CG), adjust = 2)
    denStat <- data.frame(x = denStat$x,
                          y = denStat$y)
    denStat1 <- sapply(1:length(levels(temp1$cuts)), function(i) {
      dataMat <- filter(temp1, cuts == levels(temp1$cuts)[i])
      preset <-
        density(dataMat$CG, from = min(temp1$CG), to = max(temp1$CG), adjust = 2)
      preset$y
    })
    colnames(denStat1) <- levels(temp1$cuts)
    denStat1 <- stack(as.data.frame(denStat1))
    denStat <-
      cbind(rep(denStat[,1], times = length(levels(temp1$cuts))), denStat1)
    colnames(denStat) <- c('GC', 'Density', 'ind')
    np <-
      nPlot(Density ~ GC, group = 'ind', data = denStat, type = 'lineChart')
    np$chart(useInteractiveGuideline = TRUE)
    np$xAxis(axisLabel = 'GC content')
    np$yAxis(axisLabel = 'Density')
    np
  }
}


## 3d pca plot 
pcaplot<-function (x, subset = NULL, cv.Th = 0.1, var.Th = 0, mean.Th =0, standardize = TRUE,
                   method = c("cluster", "mds","pca"), dimension = c(1,2,3), color = 'black', princurve=F,lwd=1,starts=NULL,col.curve='red', 
                   text = T, main = NULL, psi = 4, type = 'p',interactive = F, ...)
{
  
  if (is.matrix(x)) {
    dataMatrix <- x
  }
  else {
    stop("The class of \"x\" should be matrix!")
  }
  if (is.null(subset)) {
    if (!is.null(cv.Th)){
      cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
      dataMatrix <- dataMatrix[cv.gene>cv.Th,]
      subset <- 1:nrow(dataMatrix)
    }
    else{
      if (!is.null(var.Th)) {
        var.gene<-apply(dataMatrix,1,function(x) var(x))
        dataMatrix<-dataMatrix[var.gene>var.Th,]
        subset<-1:nrow(dataMatrix)
      }
      if (!is.null(mean.Th)){
        mean.gene<-apply(dataMatrix,1,function(x) mean(x))
        dataMatrix<-dataMatrix[mean.gene>mean.Th,]
        subset<-1:nrow(dataMatrix)
      }    
    }
    if (is.null(main))
      main <- paste("Sample relations based on", length(subset),
                    "genes")
  }
  else {
    if (length(subset) == 1 && is.numeric(subset)) {
      subset <- sample(1:nrow(dataMatrix), min(subset,
                                               nrow(dataMatrix)))
    }6.3877
    if (is.null(main))
      main <- paste("Sample relations based on", length(subset),
                    "selected genes")
  }
  if (standardize)
    dataMatrix <- scale(dataMatrix)
  method <- match.arg(method)
  if (method == "cluster") {
    dd <- dist(t(dataMatrix[subset, ]))
    hc = hclust(dd, "ave")
    plot(hc, xlab = "Sample", main = main, ...)
    attr(hc, "geneNum") <- length(subset)
    return(invisible(hc))
  }
  if (method=="mds") {
    dd <- dist(t(dataMatrix[subset, ]))
    mds.result <- cmdscale(dd, k = max(dimension), eig = TRUE)
    ppoints <- mds.result$points
    eig <- mds.result$eig
    percent <- round(eig/sum(eig) * 100, 1)
    if (is.null(color)) {
      color <- 1
    }
    else {
      if (!is.numeric(color)) {
        allColor <- colors()
        if (!all(is.element(color, allColor))) {
          color <- as.numeric(factor(color, levels = unique(color)))
        }
      }
    }
    plot3d(ppoints[, dimension[1]], ppoints[, dimension[2]], ppoints[,dimension[3]],
           xlab = paste("Principal Component ",
                        dimension[1], " (", percent[dimension[1]], "%)",
                        sep = ""), ylab = paste("Principal Component ",
                                                dimension[2], " (", percent[dimension[2]], "%)",
                                                sep = ""),zlab = paste("Principal Component ",
                                                                       dimension[3], " (", percent[dimension[3]], "%)",
                                                                       sep = ""), main = main,size=psi,col=color, type = type)
    if(princurve){
      start<-aggregate(ppoints[,1:3],by=list(rank(!starts)),FUN=mean)
      start <- as.matrix(start[, -1])
      fit<-principal.curve(ppoints[,1:3],start=start,plot.true=F)
      plot3d(fit$s[fit$tag,],type='l',add=T,col=col.curve,lwd=lwd)
    }
    if(text)
      text3d(ppoints[, dimension[1]], ppoints[, dimension[2]], ppoints[, dimension[3]],
             col = color, texts = colnames(dataMatrix), cex = 1)
    attr(ppoints, "geneNum") <- length(subset)
    pp<-ppoints
    if(interactive)
      scatterplot3js(ppoints[,1], ppoints[,2], ppoints[,3], 
                     labels = colnames(x), 
                     axisLabels = c(paste('PC1 (', percent[1], '%)', sep = ''),
                                    paste('PC2 (', percent[2], '%)', sep = ''),
                                    paste('PC3 (', percent[3], '%)', sep = '')),
                     color = color,
                     renderer = 'canvas', bg = 'white')
    return(pp)
  }
  if (method=="pca") {
    pca.result <- prcomp(t(dataMatrix))
    ppoints <- pca.result$x[,1:max(dimension)]
    percent<-round((pca.result$sdev^2/sum(pca.result$sdev^2))*100,1)
    if (is.null(color)) {
      color <- 1
    }
    else {
      #       if (!is.numeric(color)) {
      #         allColor <- colors()
      #         if (!all(is.element(color, allColor))) {
      #           color <- as.numeric(factor(color, levels = unique(color)))
      #         }
      #       }
    }
    plot3d(ppoints[, dimension[1]], ppoints[, dimension[2]], ppoints[,dimension[3]],
           xlab = paste("Principal Component ",
                        dimension[1], " (", percent[dimension[1]], "%)",
                        sep = ""), ylab = paste("Principal Component ",
                                                dimension[2], " (", percent[dimension[2]], "%)",
                                                sep = ""),zlab = paste("Principal Component ",
                                                                       dimension[3], " (", percent[dimension[3]], "%)",
                                                                       sep = ""), main = main,size=psi,col=color, type = type, pch =19)
    if(princurve){
      start<-aggregate(ppoints[,1:3],by=list(rank(!starts)),FUN=mean)
      start <- as.matrix(start[, -1])
      fit<-principal.curve(ppoints[,1:3],start=start,plot.true=F)
      plot3d(fit$s[fit$tag,],type='l',add=T,col=col.curve,lwd=lwd)
    }
    if(text){
      text3d(ppoints[, dimension[1]], ppoints[, dimension[2]], ppoints[, dimension[3]],
             col = color, texts = colnames(dataMatrix), cex = 1)
    }
    pp<-ppoints
    if(interactive)
      scatterplot3js(ppoints[,1], ppoints[,2], ppoints[,3], 
                     labels = colnames(x), 
                     axisLabels = c(paste('PC1 (', percent[1], '%)', sep = ''),
                                    paste('PC2 (', percent[2], '%)', sep = ''),
                                    paste('PC3 (', percent[3], '%)', sep = '')),
                     color = color,
                     renderer = 'canvas', bg = 'white')
    return(pp)
  }
  else {
    stop("the method has to be one of the three,\'cluster\',\'mds\' or \'pca\'")
  }
}


## heatmap of promoter methylation intensity
heatmap.3 <- function(Exprs, sel=F, thres_mean, thres_var, numbreaks=100, col = c("blue","white","red"), 
                       breakratio = c(2,1,2), colsidebar, Colv=F, Rowv=T, scale= 'row', labRow=F, 
                       labCol=F, dendrogram = 'row', interactive = F){
  suppressPackageStartupMessages(invisible(require('gplots', quietly=TRUE)))
  if(labRow)
    labRow <- rownames(Exprs)
  if(labCol)
    labCol <- colnames(Exprs)
  if(sel){
    gene_mean <- apply(Exprs,1,mean)
    gene_var <- apply(Exprs,1,var)
    Exprs <- Exprs[gene_mean>thres_mean & gene_var>thres_var,]
  }
  if(scale == 'row')
    Exprs_scale <- t(scale(t(Exprs)))
  else
    Exprs_scale <- Exprs
  # lmat is a matrix describing how the screen is to be broken up. By default, heatmap.2 divides the screen into a four element grid, so lmat is a 2x2 matrix. 
  # The number in each element of the matrix describes what order to plot the next four plots in. Heatmap.2 plots its elements in the following order:
  # 1 Heatmap,
  # 2 Row dendrogram,
  # 3 Column dendrogram,
  # 4 Key
  if(missing(colsidebar)){
    lmat <- rbind(c(0,4), c(0,3), c(2,1))
    lwid <- c(1,4)
    lhei <- c(1,0.1,4)
    if(class(Colv) == 'dendrogram'){
      lhei <- c(1,1,4)
      dendrogram <- 'both'
    }
  }
  else{
    if(class(Colv) == 'dendrogram'){
      # 4 is column dendrogram, 5 is key, 1 is colcolorkey
      lmat <- rbind(c(0,5),c(0,4), c(3,2),c(0,1))
      lwid <- c(1,4)
      lhei <- c(1,1, 4,0.25)
      dendrogram <- 'both'
    }
    else{
      if(Colv){
        lmat <- rbind(c(0,5),c(0,4), c(3,2),c(0,1))
        lwid <- c(1,4)
        lhei <- c(1,1, 4,0.25)
        dendrogram <- 'both'
      }
      lmat <- rbind(c(0,5),c(0, 1), c(3,2),c(0,4))
      lwid <- c(1,4)
      lhei <- c(1,0.25,4,0.1)
    }
    
  }
  rg <- quantile(Exprs_scale,na.rm=T)
  rg_diff <- rg[4]-rg[2]
  rg_max <- max(abs(rg))
  Exprs_sd <- sd(Exprs_scale)
  Exprs_mean <- mean(Exprs_scale)
  if(rg_max > max(abs(c(Exprs_mean + 3*Exprs_sd, Exprs_mean - 3*Exprs_sd)))){
    rg_iqr <- max(abs(c(rg[2], rg[4])))
    bp <- c((breakratio[1]/sum(breakratio))*rg_diff - rg_iqr, rg_iqr - (breakratio[3]/sum(breakratio))*rg_diff)
    bk <- c(seq(-rg_max, -rg_iqr, length= numbreaks), seq(-rg_iqr,bp[1],length = numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg_iqr,length=numbreaks), 
            seq(rg_iqr, rg_max, length = numbreaks))
    hmcols<- colorRampPalette(col)(length(bk)-1)
  }
  else{
    rg <- range(Exprs_scale, na.rm=T)
    bp <- c((breakratio[1]/sum(breakratio))*diff(rg) - rg_max, rg_max - (breakratio[3]/sum(breakratio))*diff(rg))
    bk <- c(seq(-rg_max,bp[1],length=numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg_max,length=numbreaks))
    hmcols<- colorRampPalette(col)(length(bk)-1)
  }
  heatmap.2(Exprs, Colv=Colv,Rowv=Rowv, dendrogram = dendrogram,trace='none',scale=scale ,density.info='none',
            lmat=lmat,lwid=lwid,lhei=lhei,labRow=labRow,labCol=labCol,col=hmcols,breaks=bk, 
            ColSideColors=colsidebar) 
  if(interactive){
    ol <- hclust(dist(Exprs_scale))$order
    vals <- unique(c(Exprs_scale))
    o <- order(vals, decreasing = FALSE)
    cols <- colorRampPalette(c("blue","white","red"))(vals)
    colz <- setNames(data.frame(vals[o], cols[o]), NULL)
    plot_ly(z = Exprs_scale[ol,],colorscale = colz, x = colnames(Exprs_scale),  y = rownames(Exprs_scale), type = "heatmap", colorbar = list(title = 'colorkey')) %>%
      layout(xaxis = list(title = ''), yaxis = list(title = ''))
  }
}

## enrichment score calculated by yidong's program