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
points_df <- read.csv(file="pearl_river/site_multi_pcr.csv",fileEncoding = "GBK")
# 计算各列的占比
points_df <- points_df %>%
  mutate(
    # 计算总和
    total = rowSums(across(where(is.numeric) & everything() & !1:4)),
    # 计算每一列的占比
    across(where(is.numeric) & everything() & !1:4, ~ .x / total, .names = "{col}_perc")
  )
points_df$lon_copy <- points_df$longitude
points_df$lat_copy <- points_df$latitude
# 将数据框转换为 sf 对象，并指定坐标参考系统
points_sf <- st_as_sf(points_df, coords = c("longitude", "latitude"), crs = 4326)

# 给河流数据添加一个唯一标识符
river_shp$river_id <- seq_along(river_shp$geometry)

points_df$r <- (points_df$total / max(points_df$total) * (10 - 6) + 3)*0.005

# 绘制地图并根据物种总和控制饼图半径
p4 <- ggplot() +
  # 绘制河流
  geom_sf(data = river_shp, aes(fill = 类型), color = NA, show.legend = "fill") +
  # 添加图例：点的大小与物种总和相关
  guides(size = guide_legend(title = "Reads number")) +
  geom_point(data = points_df, aes(x = lon_copy, y = lat_copy, size = total), color = "black",alpha = 1) + # 添加透明点用于生成图例
  scale_size_continuous(range = c(4, 8)) + # 设置图例中点的大小范围
  scale_fill_manual(name = "type", values = c("水系" = "#97dbf2", "海洋" = "#97dbf2")) +
  guides(fill = "none")+
  # 开始新的填充比例尺
  new_scale_fill() +
  # 绘制饼图，并根据总物种数控制饼图半径
  geom_scatterpie(data = points_df,
                  aes(x = lon_copy, y = lat_copy, r = r, fill = species), # 控制饼图半径
                  color = "gray",
                  size = 0,
                  #lwd = 0.02,
                  cols = c("Amphibalanus_amphitrite_perc",
                           "Tinca_tinca_perc",
                           "Coptodon_zillii_perc",
                           "Oreochromis_niloticus_perc",
                           "Limnoperna_fortunei_perc",
                           "Polyodon_spathula_perc",
                           "Oncorhynchus_kisutch_perc",
                           "Sciaenops_ocellatus_perc",
                           "Cirrhinus_mrigala_perc",
                           "Crassostrea_gigas_perc",
                           "other_perc")
                  )+
  scale_fill_manual(name = "species",
                    values = c("Amphibalanus_amphitrite_perc" = "#F27970",
                               "Tinca_tinca_perc" = "#DAA520",
                               "Coptodon_zillii_perc" = "#54B345",
                               "Oreochromis_niloticus_perc" = "#32B897",
                               "Limnoperna_fortunei_perc" = "#9DC3E7",
                               "Polyodon_spathula_perc" = "#8983BF",
                               "Oncorhynchus_kisutch_perc" = "#C76DA2",
                               "Sciaenops_ocellatus_perc" = "#FCE6CF",
                               "Cirrhinus_mrigala_perc" = "#A349A4",
                               "Crassostrea_gigas_perc" = "#5B84B1",
                               "other_perc" = "#EEE8AA"
                               ),
                    labels = c("Amphibalanus_amphitrite_perc" = "Amphibalanus amphitrite",
                               "Tinca_tinca_perc" = "Tinca tinca",
                               "Coptodon_zillii_perc" = "Coptodon zillii",
                               "Oreochromis_niloticus_perc" = "Oreochromis niloticus",
                               "Limnoperna_fortunei_perc" = "Limnoperna fortunei",
                               "Polyodon_spathula_perc" = "Polyodon spathula",
                               "Oncorhynchus_kisutch_perc" = "Oncorhynchus kisutch",
                               "Sciaenops_ocellatus_perc" = "Sciaenops ocellatus",
                               "Cirrhinus_mrigala_perc" = "Cirrhinus mrigala",
                               "Crassostrea_gigas_perc" = "Crassostrea gigas",
                               "other_perc" = "Other")) +
  theme_light() +
  labs(title = "") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) + # 添加指北针
  annotation_scale(location = "bl", width_hint = 0.1) + # 添加比例尺
  coord_sf(xlim = c(113.17, 113.83), ylim = c(22.57, 23.31)) + # 设置经纬度范围
  geom_text(
    data = points_df,
    aes(x = lon_copy, y = lat_copy, label = site),
    #color = "white",    # 设置文字颜色
    size = 2.5,           # 设置字号
    fontface = "bold",  # 设置粗体
    # check_overlap = TRUE # 防止标签重叠
  )+
  theme(legend.position = "right")# 控制图例位置
p4

ggsave("pearl_river/figure3.pdf", plot = p4, width = 9, height = 9, units = "in", dpi=300)
