library(ggplot2);library(reshape2);library(scales);library(VennDiagram)

#########################
#Define inputs
sampleName <- 'Carcinoma1'
summaryTable <- read.table('Example/Output.neoantigens.summarytable.txt', sep='\t', header=T, stringsAsFactors = F)
numRegions <- 4
epTable <- read.table('Example/Output.neoantigens.txt', sep='\t', stringsAsFactors = F)
names(epTable) <- c('Sample', paste0('Region',1:numRegions), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
recopoTable <- read.table('Example/PredictedRecognitionPotentials.txt', sep='\t', stringsAsFactors = F, header=T)
########################

colors = c('#ffd92f','#1b9e77', '#e31a1c','#1f78b4','#b07935','#a91e7d', '#525252', '#008080', '#ff6214')


#Summary by region and binding strength
region.df <- data.frame(Region = 1:numRegions,
                         WB = as.numeric(summaryTable[,paste0('Total_WB_Region_',0:(numRegions-1))]),
                         SB = as.numeric(summaryTable[,paste0('Total_SB_Region_',0:(numRegions-1))]))
region.mdf <- melt(region.df,id='Region')
p.reg <- ggplot(region.mdf, aes(x=Region, y=value, fill=variable)) + geom_bar(stat='identity') +
  scale_fill_manual(values=c('#cb6339','#a5160f')) +
  labs(y='Number or neoepitopes', fill='Binding strength') + theme_bw() + theme(text=element_text(size=14))

#Summary by binding strength and clonality
strength.df <- data.frame(Binder=c('Total', 'SB', 'WB'),
                          Clonal = as.numeric(summaryTable[,paste0('Clonal',c('','_SB', '_WB'))]),
                          Shared = as.numeric(summaryTable[,paste0('Shared',c('','_SB', '_WB'))]),
                          Subclonal = as.numeric(summaryTable[,paste0('Subclonal',c('','_SB', '_WB'))]))
strength.mdf <- melt(strength.df, id='Binder')
p.str <- ggplot(strength.mdf, aes(x=Binder, y=value, fill=variable)) + geom_bar(stat='identity', position='fill') +
  scale_y_continuous(labels=percent_format()) + scale_x_discrete(limits=c('WB', 'SB', 'Total')) +
  coord_flip() + scale_fill_manual(values=c('#ce5c55','#dba002','#4e9ed0')) +
  labs(y='Proportion of neoepitopes', fill='Clonality') + theme_bw() + theme(text=element_text(size=14))
  
#Neoepitopes in regions

numRegionsVenn <- 5

region_eps = list()
for (i in 1:numRegionsVenn){
  region_eps[[i]] <- row.names(epTable[epTable[,paste0('Region',i)]==1,])
}

venn.diagram(x = region_eps,
             category.names = paste0('Region',1:numRegionsVenn),
             filename=paste0(sampleName,'.epitopes_regions.tiff'),
             imagetype='tiff',
             units='in',
             resolution=3000,
             height=4,
             width=5.5,
             fill=colors[1:numRegionsVenn],
             cex=1.1, cat.cex=1.1, margin=0.1)

#Recognition potential in regions
recopoTable.filtered <- subset(recopoTable, NeoantigenRecognitionPotential>1e-4)
epTable.filtered <- subset(epTable, peptide %in% recopoTable.filtered$MutantPeptide)

epTable.filtered$RecognitionPotential <- recopoTable.filtered$NeoantigenRecognitionPotential[match(epTable.filtered$peptide,
                                                                                                   recopoTable.filtered$MutantPeptide)]
recopo.df <- data.frame(Region1 = epTable.filtered$Region1*epTable.filtered$RecognitionPotential)
for (i in 2:numRegions){
  recopo.df[,paste0('Region',i)] <- epTable.filtered[,paste0('Region',i)]*epTable.filtered$RecognitionPotential
}
recopo.mdf <- melt(recopo.df)
recopo.mdf <- subset(recopo.mdf, value > 1e-4)

p.recopo <- ggplot(recopo.mdf, aes(x=variable, y=value,fill=variable)) + geom_boxplot(alpha=0.7) + geom_jitter(width=0.1, height=0) +
  scale_fill_manual(values=colors) +
  labs(x='', y='Recognition potential') + guides(fill=F) + theme_bw() + theme(text=element_text(size=14))

#Save to files
pdf(paste0(sampleName,'.regions.pdf'), height=3, width=5);print(p.reg); dev.off()
pdf(paste0(sampleName,'.strength.pdf'), height=3, width=5);print(p.str); dev.off()
pdf(paste0(sampleName,'.recopo_regions.pdf'), height=3.5, width=5);print(p.recopo + scale_y_log10()); dev.off()
