# Usage: Rscript plot.r qk-11-mu.gene_num.count qk-11-mu.gene_num.count.pdf qk-11-mu
# 获取输入和输出文件路径
input<- commandArgs(trailingOnly = TRUE)[1]
output<- commandArgs(trailingOnly = TRUE)[2]
sample<- commandArgs(trailingOnly = TRUE)[3]

# Load Package
library(ggplot2)
library(dplyr)
library(tidyr)

# 生成示例数据
data <- read.table(input, header = T)

# 将数据转换为长格式
data_long <- data %>%
  gather(key = "Type", value = "Value", -CHR)
head(data_long)

data_filtered <- data_long %>% filter(.[[2]]!= "AlleNumber" & .[[2]]!= "Orther")
head(data_filtered)

# 计算百分比
data_filtered <- data_filtered %>%
  group_by(CHR) %>%
  mutate(Percentage = Value / sum(Value) * 100)
head(data_filtered)

# 绘制百分比堆积柱状图
plot_title <- paste("Genotypic statistics in ", sample, sep = "")
ggplot(data_filtered, aes(x=CHR, y=Percentage, fill=Type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Hexasomic"="#289584", "Tetrasomic"="#ccb280", "Disomic"="#3d3a6b", "NoAlle"="#65b961", "Others"="#faebd7")) +
  labs(y="Proportion(%)", x="", title = plot_title, fill="Type") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill="grey90", colour="black", size=1),
    strip.text = element_text(size=12, face="bold"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(output, width = 5, height = 3)
