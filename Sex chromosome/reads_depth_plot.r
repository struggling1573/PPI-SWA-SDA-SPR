library(qqman)
library(tidyverse)
library(ggplot2)
df <- read.table('df_MeanDepth_trans.txt', header = T)
# 1)计算chr长度
chr_len <- df %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP))
chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)
Snp_pos <- chr_pos %>%
  left_join(df, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)
X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
ggplot(Snp_pos, aes(x=BPcum, y=log2(P))) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  scale_y_continuous(expand = c(0, 0),limits=c(-0.1, 0.5)) +
  labs(title = "", x = "", y = "log2 female/male read depth ratio")+
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(filename ="reads_depth.png",width = 9,height = 3,dpi=300)
