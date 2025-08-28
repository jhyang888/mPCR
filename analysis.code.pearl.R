library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggpubr)
library(scales)
#############################################
#######pearl river
setwd("../pearl_river/")
zj.primer <- read.csv("success.primer.csv") ; head( zj.primer )
zjp <- zj.primer[ rowSums(is.na( zj.primer ) ) < 3, ] ; head( zjp ) ; dim( zjp )
df.zjp <- melt( zjp ) ; head( df.zjp )

z1 <- ggplot(df.zjp, aes(x = variable, y = species, size = value, color = variable)) +
  geom_point(alpha = 0.99) +  # 绘制气泡，并设置透明度
  scale_size_continuous(range = c(2, 5), labels = label_number(accuracy = 1), breaks = c(1, 2)) +  # 设置气泡大小范围
  scale_color_manual(values = c("reference" = "#f8984e", "new.design" = "#4d97cd"), guide = F) +
  theme_light() +  # 使用简洁的主题
  labs(x = NULL, y = NULL, size = "Number of \n successful \n primers")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
z1

###primer source
zjpsource <- read.csv("zhujiang_primer_source.csv") ; head( zjpsource )
zjpsource$gene <- factor( zjpsource$gene , levels = c("COI" , "Dloop" , "Cytb" , "16S" , "12S"))

mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9",
               "#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C",
               "#E0367A","#D8D155","#64495D","#7CC767")

z2 <- ggplot(data = zjpsource, aes(axis1 = gene, axis2 = type, y = value)) +
  geom_flow(width = 1/15 , aes(fill = gene), aes.bind = "flows")+
  geom_stratum(width = 1/20, alpha = .9, color = "black", lwd = 0.3) +
  scale_fill_manual(values = mycol) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum) ), size = 3 ) +
  scale_x_discrete(limits = c("gene", "type"), expand = c(0.15, 0.05)) +
  theme_void()+
  theme(legend.position = "none" )
z2

library(dbplyr)
##read pier river
mpcr <- read.csv("../pearl_river/multi_pcr.csv" , row.names = 1) ; head( mpcr ) ; dim( mpcr )
##sequence composition pie chart
dfseqcom <- data.frame( species = names(sort(rowSums(mpcr) , decreasing = T)), 
                        reads = sort(rowSums(mpcr) , decreasing = T))

dfseqtop8 <- rbind( dfseqcom[1:8,] , data.frame(species = "Other" , reads = sum(dfseqcom[9:nrow( dfseqcom),2])) )
dfseqtop8$species <- factor( dfseqtop8$species , levels = dfseqtop8$species)
z3 <- ggplot(dfseqtop8, aes(x="", y=reads, fill=species)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  #scale_fill_manual(values = mycol) +
  scale_fill_brewer(palette="Set1")+
  theme_void() +# remove background, grid, numeric labels
  theme(strip.text = element_text(size = 7 ),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "lines"), legend.position = "left")
z3

spdf <- data.frame(species = rownames(mpcr), frequency = as.numeric(100*rowSums( mpcr > 10 )/ncol( mpcr ))) ; head( spdf )
spdf <- spdf[ order( spdf$frequency , decreasing = F), ]
spdf$species <- factor( spdf$species , levels = spdf$species)

sitedf <- data.frame(site = colnames(mpcr) , species = as.numeric(colSums( mpcr > 10 )) ) ; head( sitedf )
sitedf <- sitedf[ order( as.numeric( gsub("L" ,"" , sitedf$site)) , decreasing = T), ]
sitedf$site <- factor( sitedf$site , levels = sitedf$site)

z4 <- ggplot( data = spdf , aes(x = frequency , y = species))+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label= frequency), hjust=1, color="white", size= 2)+
  labs(x = "Frequency", y = NULL) +
  theme_light()
z4

z5 <- ggplot( data = sitedf , aes(x = species , y = site))+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label= species), hjust=1, color="white", size= 2)+
  labs(x = "Species Number", y = "Site") +
  theme_light()
z5

pblank <- ggplot()+theme_void()
zz1 <- ggarrange( z4, z5, ncol = 2 , labels = c("d" , "e") , widths = c(2,1.5)) ; zz1
zz2 <- ggarrange( z3, ggarrange(pblank,z2, widths = c(0.3,2)) , widths = c(3,2) , labels = c("b" , "c")) ; zz2
#zz2 <- ggarrange( ggarrange(pblank,z2, widths = c(1,2)) , pblank, widths = c(2,1)) ; zz2
zz3 <- ggarrange(zz2 , zz1 , ncol = 1 , labels = c("b" , "d") , heights = c(1,3))
zz4 <- ggarrange( z1 , zz3 , ncol = 2, widths = c(1,1.5), labels = c("a" , "b")) ; zz4

pdf("figure2.pdf" , width = 10, height = 8)
zz4
dev.off()
###
dfsitseq1 <- mpcr[ rownames(mpcr) %in% dfseqtop8$species, ]
dfsitseq <- rbind( dfsitseq1 , colSums(mpcr[ !rownames(mpcr) %in% dfseqtop8$species, ]) )
rownames(dfsitseq)[9] <- "Other"
df.sitseq <- melt( t(dfsitseq ) )
colnames(df.sitseq) <- c("site" , "species" , "reads")
df.sitseq$site <- factor( df.sitseq$site , levels = colnames(dfsitseq )[order( as.numeric( gsub("L" , "" , colnames( dfsitseq))))])

library(RColorBrewer)
z6 <- ggplot(df.sitseq, aes(x=site, y=reads, fill=species))+
  #geom_bar(stat = 'identity', width = 0.85, linewidth=0.2, colour='gray', position="stack")+
  geom_bar(stat = 'identity', width = 0.85, linewidth=0.2, colour='gray', position="fill")+
  theme_light()+
  scale_fill_manual(values = mycol[3:12])+
  #scale_fill_manual(values = brewer.pal(9,"Set1"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+
  labs(x = "Site" , y = "Realative abundance")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"))
z6

pdf("figureS2.pdf" , width = 10, height = 4)
z6
dev.off()

##read 12s results
metapcr <- read.csv("12s_asvs.table.tax.csv", row.names = 1, fileEncoding = "GBK") ; head( metapcr)
##rename row
rownames( metapcr ) <- metapcr$species
##delete some sample
metapcr <- metapcr[, colnames( metapcr ) %in% colnames(mpcr)] ; head( metapcr ) ; dim( metapcr )
##order the smaple
metapcr <- metapcr[ , order( colnames( metapcr ) ) ]; head( metapcr ) ; dim( metapcr )
##melt data

library(reshape2)
rownames( mpcr ) %in% rownames( metapcr ) %>% sum()

df.mpcr <- melt( t(mpcr) ) ; head( df.mpcr )
colnames( df.mpcr) <- c("site" , "species" , "reads"); head( df.mpcr )
#df.mpcr$method <- "Amplicon_seq" ; head( df.mpcr )

df.metapcr <- melt( t(metapcr) ) ; head( df.metapcr )
colnames( df.metapcr) <- c("site" , "species" , "reads"); head( df.metapcr )
#df.metapcr$method <- "Metabarcoding" ; head( df.metapcr )
##merge metabarcoding and amplicon seq
dfshare <- merge( df.metapcr , df.mpcr , by = c("species" , "site")) ; head( dfshare )
colnames(dfshare ) <- c("species" , "site" , "metabarcoding" , "amplicon_seq"); head( dfshare )
dfshare$site <- factor(dfshare$site , levels = colnames(mpcr)[order( as.numeric( gsub("L" , "" , colnames(mpcr))))])
###
dfsummary <- data.frame(c(sum(metapcr) , sum(mpcr)),
           c( sum( rowSums( metapcr ) > 10 ), sum( rowSums( mpcr) > 10 )))
colnames( dfsummary ) <- c("reads" , "species")
dfsummary$type <- c("metabarcoding" , "mPCR")

k0 <- ggplot( data = dfsummary , aes(x = type , y = reads ))+
  geom_bar(stat = "identity" , fill = brewer.pal(11,"Spectral")[10] )+
  theme_light()+
  labs( x = NULL , y = "Reads number")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"))
k0

k1 <- ggplot( data = dfsummary , aes(x = type , y = species ))+
  geom_bar(stat = "identity" , fill = brewer.pal(11,"Spectral")[9] )+
  theme_light()+
  labs( x = NULL , y = "Species number")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"))
k1

d1 <- data.frame( rowSums( metapcr )/sum(metapcr ) )
colnames( d1 ) <- "prop" ; d1$species  = rownames( d1) ; d1$type  = "metabarcoding"

d2 <- data.frame( rowSums( mpcr[rownames( mpcr ) %in% rownames( d1),] )/sum(mpcr[rownames( mpcr ) %in% rownames( d1),] ) )
colnames( d2 ) <- "prop" ; d2$species  = rownames( d2 ) ; d2$type  = "mPCR"

dfprop <- rbind( d1, d2 )

dfprop$species <- factor(dfprop$species , levels = levels(dfshare$species) )
library(ggplot2)
library(ggalluvial)
k2 <- ggplot( data = dfprop , aes(x = type , y = prop , fill = species))+
  geom_bar( stat = "identity" , position = "fill")+
  geom_flow(aes(alluvium = species), alpha = 0.5) + # 添加柱状图后的条带
  scale_fill_manual(values = brewer.pal(11,"Spectral"))+
  theme_light()+
  labs( x = NULL , y = "Proportion")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"))
k2

z7 <- ggplot(dfshare, aes(x=site, y=metabarcoding, fill=species))+
  geom_bar(stat = 'identity', width = 0.85, linewidth=0.2, colour='gray', position="fill")+
  theme_light()+
  #scale_fill_manual(values = mycol)+
  scale_fill_manual(values = brewer.pal(11,"Spectral"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+
  labs(x = NULL , y = "Reads number")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        legend.position = "none")
z7

z7.2 <- ggplot(dfshare, aes(x=site, y=metabarcoding))+
  geom_bar(stat = 'identity', width = 0.85, linewidth=0, position="stack")+
  theme_light()+
  #scale_fill_manual(values = mycol)+
  scale_fill_manual(values = brewer.pal(11,"Spectral"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)), labels = scales::scientific_format(digits = 0))+
  labs(x = NULL , y = "Reads")+
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        legend.position = "none")
z7.2
z73 <- ggarrange(z7.2, z7, ncol = 1, heights = c(1,2.5), labels = c("d" , "e") ) ; z73

z8 <- ggplot(dfshare, aes(x=site, y=amplicon_seq, fill=species))+
  #geom_bar(stat = 'identity', width = 0.85, linewidth=0.2, colour='gray', position="stack")+
  geom_bar(stat = 'identity', width = 0.85, linewidth=0.2, colour='gray', position="fill")+
  theme_light()+
  #scale_fill_manual(values = mycol)+
  scale_fill_manual(values = brewer.pal(11,"Spectral"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+
  labs(x = NULL , y = "Relative abundance")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        legend.position = "none")
z8

z8.2 <- ggplot(dfshare, aes(x=site, y=amplicon_seq))+
  geom_bar(stat = 'identity', width = 0.85, linewidth=0, position="stack")+
  theme_light()+
  #scale_fill_manual(values = mycol)+
  scale_fill_manual(values = brewer.pal(11,"Spectral"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.05),add=c(0,0)))+
  labs(x = NULL , y = "Reads")+
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"))
z8.2
z83 <- ggarrange(z8.2, z8, ncol = 1, heights = c(1,2.5), labels = c("f" , "g") ) ; z83

kk1 <- ggarrange( k0, k1, k2, pblank, widths = c(1,1,1.7,1), ncol = 4 , labels = c("a" , "b" , "c" , "") ); kk1
zz9 <- ggarrange(kk1, z73, z83, ncol = 1, heights = c(1.2,2,2), labels = c("a" , "d" , "f") ) ; zz9

pdf("figure4.pdf", width = 9, height = 10)
zz9
dev.off()

##correlation
library(ggpubr)
z9 <- ggplot(data = dfshare[dfshare$species %in% c("Channa striata","Clarias gariepinus",
                                                   "Coptodon zillii","Gambusia affinis",
                                                   "Oreochromis aureus",
                                                   "Oreochromis niloticus","Pterygoplichthys pardalis"),] , 
             aes( x = metabarcoding , y = amplicon_seq))+
  geom_point(size = 2, pch = 21, color = "gray")+
  labs(x = "Metabarcoding reads number", y = "mPCR reads number") +
  theme_light() +
  facet_wrap(vars(species) , scales = "free" , ncol = 4)+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 8),
        axis.text.x = element_text(color = "black",size = 8),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "red", method = "pearson") +
  stat_regline_equation(size = 3, color = "black",label.y.npc = 1)
z9

z11 <- ggplot(data = dfshare[!dfshare$species %in% c("Channa striata","Clarias gariepinus",
                                                   "Coptodon zillii","Gambusia affinis",
                                                   "Oreochromis aureus",
                                                   "Oreochromis niloticus","Pterygoplichthys pardalis"),] , 
             aes( x = metabarcoding , y = amplicon_seq))+
  geom_point(size = 2, pch = 21, color = "gray")+
  labs(x = "Metabarcoding reads number", y = "mPCR reads number") +
  theme_light() +
  facet_wrap(vars(species) , scales = "free" , ncol = 2)+
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 8),
        axis.text.x = element_text(color = "black",size = 8),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "red", method = "pearson") +
  stat_regline_equation(size = 3, color = "black",label.y.npc = 1)
z11

pdf("fugureS3.pdf" , width = 6 , height = 5)
z11
dev.off()
##########
df.div <- data.frame(amplicon_seq = colSums(mpcr>10), metabarcoding = colSums( metapcr>110))
z10 <- ggplot(data = df.div , aes( x = metabarcoding , y = amplicon_seq))+
  geom_point(size = 2, pch = 21, color = "gray")+
  labs(x = "Non-native species richness \n detected by metabarcoding", y = "Non-native species \n detected by mPCR") +
  theme_light() +
  theme(strip.text = element_text(size = 10),
        axis.line = element_line(color = "black",size = 0.4),
        axis.text.y = element_text(color = "black",size = 8),
        axis.text.x = element_text(color = "black",size = 8),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2 , color = "#e5e5e5"),
        legend.position = "none")+
  geom_smooth(method = "lm", size = 1, se = T, color = "black", linetype = "dashed")+
  stat_cor(size = 3, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), group = 1), color = "red3", method = "pearson") +
  stat_regline_equation(size = 3, color = "black", label.y.npc = 1)
z10

zz10 <- ggarrange( z9, ggarrange(z10, pblank, widths = c(3,3)), heights = c(2, 1.5), ncol = 1, labels = c("a" , "b")) ; zz10
zz10 <- ggarrange( z9, ggarrange(z11, ggarrange(pblank, z10, ncol = 1, heights = c(1,2), labels = c("" , "c")), widths = c(1,1), labels = c("b" , "")), heights = c(2, 2), ncol = 1, labels = c("a" , "b")) ; zz10

pdf("fugure5.pdf" , width = 10 , height = 12)
zz10
dev.off()
############

