library(reshape2)
library(ggplot2)
library(scales)
setwd("lab/")
suc.primer <- read.csv("success.amplify.primer.csv") ; head( suc.primer )
df.suc.primer <- melt( suc.primer ) ; head( df.suc.primer )
pinfo <- read.csv("primer.info.csv") ; head( pinfo )
pinfo2 <- read.csv("primer.info2.csv") ; head( pinfo2 )
pinfo2$gene <- factor( pinfo2$gene , levels = c("COI" , "Dloop" , "Cytb" , "16S" , "12S"))

mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9",
           "#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C",
           "#E0367A","#D8D155","#64495D","#7CC767")

p1 <- ggplot(df.suc.primer, aes(x = variable, y = species, size = value, color = variable)) +
  geom_point(alpha = 0.9) +  # 绘制气泡，并设置透明度
  scale_size_continuous(range = c(2, 5), labels = label_number(accuracy = 1), breaks = c(1, 2) ) +  # 设置气泡大小范围
  scale_color_manual(values = c("reference" = "#f8984e", "new.design" = "#4d97cd"), guide = F) +
  theme_light() +  # 使用简洁的主题
  labs(x = NULL, y = NULL, size = "Number of \n successful \n primers")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
p1

p2 <- ggplot(data.frame(x = pinfo$Tm), aes(x = x)) +
  geom_density(fill = "gray", alpha = 0.5) +
  labs( x = "PCR Tm", y = "Density") +
  xlim(55,70)+
  theme_light()
p2

p3 <- ggplot(data.frame(x = pinfo$length), aes(x = x)) +
  geom_density(fill = "gray", alpha = 0.5 ) +
  labs( x = "Amplicon length (bp)", y = "Density") +
  xlim(-100,500)+
  theme_light()
p3

# install.packages("ggalluvial")
library(ggalluvial)
library(ggpubr)

p4 <- ggplot(data = pinfo2, aes(axis1 = gene, axis2 = type, y = value , label = gene)) +
  geom_flow(width = 1/15 , aes(fill = gene), aes.bind = "flows")+
  geom_stratum(width = 1/20, alpha = .9, lwd = 0.3) +
  scale_fill_manual(values = mycol) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum) ), size = 3 ) +
  scale_x_discrete(limits = c("gene", "type"), expand = c(0.15, 0.05)) +
  theme_void()+
  theme(legend.position = "none" )
p4

pp1 <- ggarrange( p4 , p2, p3 , ncol = 1 , heights = c(2,1.5, 1.5) , labels = c("b" , "c" , "d"))

pp2 <- ggarrange( p1, pp1 , ncol = 2 , widths = c(2,2), labels = c("a" , "b")) ; pp2

pdf("figure1.pdf")
pp2
dev.off()

