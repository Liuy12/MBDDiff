### generate dynamic visualization for breast cancer dataset 
####### GC content enrichment 
temp1 <- bind_cols(as.data.frame(Whole_genome_100bp[,4]), as.data.frame(Whole_genome_100bp_bed[, 4]))
temp2 <- c(0, 1, 2, 5, 72363)
cuts <- cut(temp1[,1], breaks = unique(temp2), labels = c('0-25%','25-50%','50-75%','75-100%'))
temp1 <- bind_cols(temp1, as.data.frame(cuts))
colnames(temp1) <- c('count', 'CG', 'cuts')

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
#np$addParams(dom = "Density")
np$chart(useInteractiveGuideline = TRUE)
np$xAxis(axisLabel = 'GC content')
np$yAxis(axisLabel = 'Density')
np$save('GCenrich.html', standalone = TRUE)
return(np)  

######## histogram promoter count vs background count
library(plotly)
dataMat <- read.delim('../Desktop/temp.txt', sep = '\t')

Background <- dataMat[dataMat$group == 'Background',2]
Promoter <- dataMat[dataMat$group == 'Promoter', 2]

ggplot() + 
  geom_histogram(aes(x = count, y = ..density.., fill=group, alpha = 0.6),  data = dataMat, position = 'identity')+  scale_fill_manual(values = c('blue', 'red')) +
 # xlim(c(0,150)) +
  scale_fill_manual(values = c('blue', 'red')) + guides(alpha=FALSE)


plot_ly(x = Background, opacity = 0.6, type = "histogram", histnorm = 'probability density', xbins = list(start = 0, end = 150, size = 5), autobinx = F, name = 'Background') %>%
  add_trace(x = Promoter, filename="overlaid-histogram", name = 'Promoter') %>%
  layout(barmode="overlay")

########## pca3d 
temp <- pcaplot(as.matrix(log2(dataMat)), method = 'mds', text =F, cv.Th = 1.3, psi=8, color = as.factor(c(rep('T', 3), rep('N', 3))))

scatter3d <- scatterplot3js(temp[[1]][,1], temp[[1]][,2], temp[[1]][,3], 
                            labels = rownames(temp[[1]]), 
                            axisLabels = c(paste('PC1 (', temp[[2]][1], '%)', sep = ''),
                                           paste('PC2 (', temp[[2]][2], '%)', sep = ''),
                                           paste('PC3 (', temp[[2]][3], '%)', sep = '')),
                            color = rep(c('#00FFFF', '#FFE4C4'), each = 3),
                            renderer = 'canvas', bg = 'white')



########### heatmap
dataMat_sel_scaled <- t(scale(t(dataMat_sel)))
ol <- hclust(dist(dataMat_sel_scaled))$order
vals <- unique(c(dataMat_sel_scaled))
o <- order(vals, decreasing = FALSE)
cols <- colorRampPalette(c("blue","white","red"))(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)
plot_ly(z = dataMat_sel_scaled[ol,],colorscale = colz, x = colnames(dataMat_sel),  y = rownames(dataMat_sel), type = "heatmap", colorbar = list(title = 'colorkey')) %>%
  layout(xaxis = list(title = ''), yaxis = list(title = ''))
#d3heatmap(dataMat_sel_scaled, scale = "row", colors = colorRampPalette(c("blue","white","red"))(1000), Colv = FALSE)

###### data.table
datatable(dataMat, options = list(pageLength = 5))

###### DE table
TestStat <- read.delim('../Desktop/Teststas.txt')
DE_index <- which(abs(TestStat$log2FoldChange) > 1 & TestStat$padj < 0.01 & TestStat$baseMean > 30)
temp <- cbind(dataMat[DE_index,], TestStat[DE_index,c(6,8)])
datatable(temp, options = list(pageLength = 5))

### MAplot
temp <- TestStat[,c(1, 2,6)]
temp$DE_Status <- 'Not DE'
temp$DE_Status[DE_index] <- 'DE'
temp$baseMean <- log2(temp$baseMean + 0.25)

mp <-
  mjs_plot(temp, baseMean, log2FoldChange, decimals = 6) %>%
  mjs_point(
    color_accessor = DE_Status, color_range = c('red', 'grey32'), color_type = "category", x_rug =
      TRUE, y_rug = TRUE
  ) %>%
  mjs_add_baseline(y_value = 0, label = 'baseline') %>%
  mjs_labs(x_label = 'Log2 methylation levels', y_label = "Log2 fold change") %>%
  mjs_add_mouseover("function(d) {
                    $('{{ID}} svg .mg-active-datapoint')
                    .text('Gene Name: ' +  d.point.id + ',' + ' Log2 intensity: ' + d.point.baseMean + ',' + ' Log2 fold change: ' + d.point.log2FoldChange);
                    }")


