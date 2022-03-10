library(magrittr)

x <- read.table("GSE143259_raw_counts.txt", header = TRUE, row.names = 1)

x <- x[apply(x, 1, function(x) all(x>10)),]

coe <- colSums(x) / min(colSums(x))

y <- apply(x, 1, function(x) x*coe) %>% t() %>% as.data.frame()

cv_wt <- apply(y[1:4], 1, function(x) sd(x)/mean(x)*100) %>% as.data.frame()

cv_ex <- apply(y[5:8], 1, function(x) sd(x)/mean(x)*100) %>% as.data.frame()

z <- y[(cv_wt < 40 & cv_ex < 40), ]

t <- apply(z, 1, function(x) t.test(x[1:4], x[5:8]))

z_mean <- data.frame(WT = apply(z[1:4], 1, mean), EX = apply(z[5:8], 1, mean))

z_mean$logFC <- log2(z_mean$EX / z_mean$WT)

z_mean$pvalue <- t

z_mean$pvalue %<>% as.character()

z_mean$pvalue <-
  stringr::str_replace(z_mean$pvalue, ".+, p\\.value = ([0-9\\.e\\-]+), conf\\.int.+", "\\1") %>% as.numeric()

top <- z_mean[abs(z_mean$logFC) > 1 & z_mean$pvalue < 0.05, ] %>% .[3:4]

library(ggplot2)
z_plot <- z_mean[3:4]
z_plot$threshold <-
  as.factor(ifelse(abs(z_plot$logFC) >= 1,
    ifelse(z_plot$logFC > 1 , 'Up', 'Down'),
    'Not'
  ))

p <- ggplot(data = z_plot) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_point(aes(x = logFC, y = -log10(pvalue), color = threshold), size = 1) +
  xlim(c(-3, 3)) +
  geom_vline(
    xintercept = c(-1, 1),
    lty = 4,
    col = "grey",
    lwd = 0.6
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    lty = 4,
    col = "grey",
    lwd = 0.6
  ) +
  labs(x = "log2FoldChange",
       y = "-log10(p-value)",
       title = "Volcano plot of DEG")

p

plotly::ggplotly(p)
