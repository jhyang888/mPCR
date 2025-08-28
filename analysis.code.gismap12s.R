# 包导入
library(showtext)
library(ggspatial)
library(ggplot2)
library(sf)
library(pacman)
library(scatterpie) # 用于绘制饼图
p_load(sf, dplyr, ggplot2)
library(ggnewscale)
#font_add("黑体", "C:\\WINDOWS\\FONTS\\SIMHEI.TTF")
showtext_auto(enable = TRUE)

# 读取shp文件，导入底图
river_shp <- st_read("pearl_river/pearl_map/Export_Output.shp")
river_shp <- st_transform(river_shp, 4326)

# 读取点位信息和对应数据
points_df_12s <- read.csv(file="pearl_river/site_12s_asvs.table.tax.csv",fileEncoding = "GBK")
points_df_12s <- points_df_12s[!points_df_12s$site %in% c("L7" , "L8" ,"L9" ,"L12", "L19" , "L20" , "L21","L22"),]
# 计算各列的占比
points_df_12s <- points_df_12s %>%
  mutate(
    # 计算总和
    total = rowSums(across(where(is.numeric) & everything() & !1:4)),
    # 计算每一列的占比
    across(where(is.numeric) & everything() & !1:4, ~ .x / total, .names = "{col}_perc")
  )
points_df_12s$lon_copy <- points_df_12s$longitude
points_df_12s$lat_copy <- points_df_12s$latitude

# 将数据框转换为 sf 对象，并指定坐标参考系统
points_sf_12s <- st_as_sf(points_df_12s, coords = c("longitude", "latitude"), crs = 4326)

# 给河流数据添加一个唯一标识符
river_shp$river_id <- seq_along(river_shp$geometry)
# 饼图半径
points_df_12s$r <- (points_df_12s$total / max(points_df_12s$total) * (10 - 6) + 3)*0.005

# 绘制地图并根据物种总和控制饼图半径
p4 <- ggplot() +
  # 绘制河流
  geom_sf(data = river_shp, aes(fill = 类型), color = NA, show.legend = "fill") +
  # 添加图例：点的大小与物种总和相关
  guides(size = guide_legend(title = "Reads number")) +
  geom_point(data = points_df_12s, aes(x = lon_copy, y = lat_copy, size = total), color = "black",alpha = 1) + # 添加透明点用于生成图例
  scale_size_continuous(range = c(4, 8)) + # 设置图例中点的大小范围
  scale_fill_manual(name = "type", values = c("水系" = "#97dbf2", "海洋" = "#97dbf2")) +
  guides(fill = "none")+
  # 开始新的填充比例尺
  new_scale_fill()+
  # 绘制饼图，并根据总物种数控制饼图半径
  geom_scatterpie(data = points_df_12s,
                  aes(x = lon_copy, y = lat_copy, r = r,fill = species), # 控制饼图半径
                  color = "gray",
                  size = 0.02,
                  cols = c("Oreochromis_niloticus_perc",
                           "Coptodon_zillii_perc",
                           "Oreochromis_aureus_perc",
                           "Pterygoplichthys_pardalis_perc",
                           "Clarias_gariepinus_perc", 
                           "Gambusia_affinis_perc", 
                           "Cirrhinus_mrigala_perc", 
                           "Channa_striata_perc",
                           "Prochilodus_lineatus_perc",
                           "Micropterus_salmoides_perc",
                           "Anguilla_rostrata_perc")
                  )+
  scale_fill_manual(name = "species",
                    values = c("Oreochromis_niloticus_perc" = "#F27970",
                               "Coptodon_zillii_perc" = "#DAA520",
                               "Oreochromis_aureus_perc" = "#54B345",
                               "Pterygoplichthys_pardalis_perc" = "#32B897",
                               "Clarias_gariepinus_perc" = "#9DC3E7", 
                               "Gambusia_affinis_perc" = "#8983BF", 
                               "Cirrhinus_mrigala_perc" = "#C76DA2", 
                               "Channa_striata_perc" = "#FCE6CF",
                               "Prochilodus_lineatus_perc" = "#A349A4",
                               "Micropterus_salmoides_perc" = "#5B84B1",
                               "Anguilla_rostrata_perc" = "#EEE8AA"
                               ),
                    labels = c("Oreochromis_niloticus_perc" = "Oreochromis niloticus",
                               "Coptodon_zillii_perc" = "Coptodon zillii",
                               "Oreochromis_aureus_perc" = "Oreochromis aureus",
                               "Pterygoplichthys_pardalis_perc" = "Pterygoplichthys pardalis",
                               "Clarias_gariepinus_perc" = "Clarias gariepinus", 
                               "Gambusia_affinis_perc" = "Gambusia affinis", 
                               "Cirrhinus_mrigala_perc" = "Cirrhinus mrigala", 
                               "Channa_striata_perc" = "Channa striata",
                               "Prochilodus_lineatus_perc" = "Prochilodus lineatus",
                               "Micropterus_salmoides_perc" = "Micropterus salmoides",
                               "Anguilla_rostrata_perc" = "Anguilla rostrata")) +
  
  theme_light() +
  labs(title = "") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) + # 添加指北针
  annotation_scale(location = "bl", width_hint = 0.1) + # 添加比例尺
  coord_sf(xlim = c(113.17, 113.83), ylim = c(22.57, 23.31)) + # 设置经纬度范围
  geom_text(
    data = points_df_12s,
    aes(x = lon_copy, y = lat_copy, label = site),
    #color = "white",    # 设置文字颜色
    size = 2.5,           # 设置字号
    fontface = "bold",  # 设置粗体
    # check_overlap = TRUE # 防止标签重叠
  )+
  theme(legend.position = "right")# 控制图例位置
p4

########################multi amplicon seq
# 读取点位信息和对应数据
points_df_mpcr <- read.csv(file="pearl_river/site_mpcr_share_species.csv",fileEncoding = "GBK")
# 计算各列的占比
points_df_mpcr <- points_df_mpcr %>%
  mutate(
    # 计算总和
    total = rowSums(across(where(is.numeric) & everything() & !1:4)),
    # 计算每一列的占比
    across(where(is.numeric) & everything() & !1:4, ~ .x / total, .names = "{col}_perc")
  )
points_df_mpcr$lon_copy <- points_df_mpcr$longitude
points_df_mpcr$lat_copy <- points_df_mpcr$latitude

# 将数据框转换为 sf 对象，并指定坐标参考系统
points_sf_mpcr <- st_as_sf(points_df_mpcr, coords = c("longitude", "latitude"), crs = 4326)

# 给河流数据添加一个唯一标识符
river_shp$river_id <- seq_along(river_shp$geometry)
# 饼图半径
points_df_mpcr$r <- (points_df_mpcr$total / max(points_df_mpcr$total) * (10 - 6) + 3)*0.005

# 绘制地图并根据物种总和控制饼图半径
p5 <- ggplot() +
  # 绘制河流
  geom_sf(data = river_shp, aes(fill = 类型), color = NA, show.legend = "fill") +
  # 添加图例：点的大小与物种总和相关
  guides(size = guide_legend(title = "Reads number")) +
  geom_point(data = points_df_mpcr, aes(x = lon_copy, y = lat_copy, size = total), color = "black", alpha = 1) + # 添加透明点用于生成图例
  scale_size_continuous(range = c(4, 8)) + # 设置图例中点的大小范围
  scale_fill_manual(name = "type", values = c("水系" = "#97dbf2", "海洋" = "#97dbf2")) +
  guides(fill = "none")+
  # 开始新的填充比例尺
  new_scale_fill() +
  # 绘制饼图，并根据总物种数控制饼图半径
  geom_scatterpie(data = points_df_mpcr,
                  aes(x = lon_copy, y = lat_copy, r = r,fill = species), # 控制饼图半径
                  color = "gray",
                  size = 0.02,
                  cols = c("Oreochromis_niloticus_perc",
                           "Coptodon_zillii_perc",
                           "Oreochromis_aureus_perc",
                           "Pterygoplichthys_pardalis_perc",
                           "Clarias_gariepinus_perc", 
                           "Gambusia_affinis_perc", 
                           "Cirrhinus_mrigala_perc", 
                           "Channa_striata_perc",
                           "Prochilodus_lineatus_perc",
                           "Micropterus_salmoides_perc",
                           "Anguilla_rostrata_perc")
                  )+

  scale_fill_manual(name = "species",
                    values = c("Oreochromis_niloticus_perc" = "#F27970",
                               "Coptodon_zillii_perc" = "#DAA520",
                               "Oreochromis_aureus_perc" = "#54B345",
                               "Pterygoplichthys_pardalis_perc" = "#32B897",
                               "Clarias_gariepinus_perc" = "#9DC3E7", 
                               "Gambusia_affinis_perc" = "#8983BF", 
                               "Cirrhinus_mrigala_perc" = "#C76DA2", 
                               "Channa_striata_perc" = "#FCE6CF",
                               "Prochilodus_lineatus_perc" = "#A349A4",
                               "Micropterus_salmoides_perc" = "#5B84B1",
                               "Anguilla_rostrata_perc" = "#EEE8AA"
                    ),
                    labels = c("Oreochromis_niloticus_perc" = "Oreochromis niloticus",
                               "Coptodon_zillii_perc" = "Coptodon zillii",
                               "Oreochromis_aureus_perc" = "Oreochromis aureus",
                               "Pterygoplichthys_pardalis_perc" = "Pterygoplichthys pardalis",
                               "Clarias_gariepinus_perc" = "Clarias gariepinus", 
                               "Gambusia_affinis_perc" = "Gambusia affinis", 
                               "Cirrhinus_mrigala_perc" = "Cirrhinus mrigala", 
                               "Channa_striata_perc" = "Channa striata",
                               "Prochilodus_lineatus_perc" = "Prochilodus lineatus",
                               "Micropterus_salmoides_perc" = "Micropterus salmoides",
                               "Anguilla_rostrata_perc" = "Anguilla rostrata")) +
  
  theme_light() +
  #guides(fill = "none")+
  labs(title = "") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) + # 添加指北针
  annotation_scale(location = "bl", width_hint = 0.1) + # 添加比例尺
  coord_sf(xlim = c(113.17, 113.83), ylim = c(22.57, 23.31)) + # 设置经纬度范围
  geom_text(
    data = points_df_mpcr,
    aes(x = lon_copy, y = lat_copy, label = site),
    #color = "white",    # 设置文字颜色
    size = 2.5,           # 设置字号
    fontface = "bold",  # 设置粗体
    # check_overlap = TRUE # 防止标签重叠
  )+
  theme(legend.position = "right")# 控制图例位置
p5

library(ggpubr)

ggsave("pearl_river/figureS3.pdf",plot = ggarrange(p4, p5), width = 16,height = 9,units = "in",dpi=300)
