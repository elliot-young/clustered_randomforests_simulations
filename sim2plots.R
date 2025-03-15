# Plots for Simulation 2
# Plot 1 - QQ plots for asymptotic normality verification
# Plot 2 - Coverage plots for different alpha-level tests.

# Plot 1 - QQ plots
D=1
UNW = read.csv(paste0("~Downloads/output_",D,".csv"))
CRF = read.csv(paste0("~Downloads/output_",D,".csv"))

library(ggplot2)
ggplot(data.frame(UNW$estimators_unw), aes(sample = data)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  ggtitle("Q-Q Plot")


library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

coul <- brewer.pal(10, "Set1")
SIZE <- 1.0
QQplot <- function(CLASS,D) {
  if (CLASS=="UNW") {
    UNW = read.csv(paste0("~/Downloads/output",D,".csv"))
    data <- UNW$estimators_unw
    myplot <- ggplot(data.frame(data), aes(sample = data)) +
      stat_qq(size=SIZE, color=coul[3]) +
      stat_qq_line(color = "black") +
      xlab("Theoretical Quantiles") +
      ylab("Sample Quantiles") +
      ggtitle(paste0("d = ",D)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(c(-0.35,0.35))
  }
  if (CLASS=="CRF") {
    CRF = read.csv(paste0("~/Downloads/output",D,".csv"))
    data <- CRF$estimators_crf
    myplot <- ggplot(data.frame(data), aes(sample = data)) +
      stat_qq(size=SIZE, color=coul[5]) +
      stat_qq_line(color = "black") +
      xlab("Theoretical Quantiles") +
      ylab("Sample Quantiles") +
      theme_bw() +
      ylim(c(-0.35,0.35))
  }
  return(myplot)
}
plot_UNW_D1 <- QQplot("UNW",1); plot_UNW_D3 <- QQplot("UNW",3)
plot_UNW_D5 <- QQplot("UNW",5); plot_UNW_D10 <- QQplot("UNW",10)
plot_UNW_D50 <- QQplot("UNW",50); plot_CRF_D1 <- QQplot("CRF",1)
plot_CRF_D3 <- QQplot("CRF",3); plot_CRF_D5 <- QQplot("CRF",5)
plot_CRF_D10 <- QQplot("CRF",10); plot_CRF_D50 <- QQplot("CRF",50)

grid.arrange(plot_UNW_D1, plot_UNW_D3, plot_UNW_D5, plot_UNW_D10, plot_UNW_D50,
                        plot_CRF_D1, plot_CRF_D3, plot_CRF_D5, plot_CRF_D10, plot_CRF_D50,
                        ncol=5)

legend_metadata = data.frame(Forest=c(rep("UNW",4),rep("CRF",4)), xval=1:8, yval=1:8)
legend_plot = ggplot(legend_metadata, aes(x=xval, y=yval, group=Forest, color=Forest)) +
  geom_point(position=position_dodge(width = 0.001), size=3.5) +
  theme_bw() +
  scale_color_manual(
    values=c(coul[3],coul[5]),
    labels=c("Classical random forests (RF)","Clustered random forests (CRF)"),
    name="",
    guide = guide_legend(override.aes = list(
      shape = c(16, 16),
      linetype = c("solid","solid")
    ))
  ) +
  theme(legend.direction = "horizontal",
        legend.text = element_text(size=15))
legend = get_legend(legend_plot)

cairo_pdf("~/Downloads/crf-sim2-qqplots.pdf", 15, 6)
ggarrange(plot_UNW_D1, plot_UNW_D3, plot_UNW_D5, plot_UNW_D10, plot_UNW_D50,
          plot_CRF_D1, plot_CRF_D3, plot_CRF_D5, plot_CRF_D10, plot_CRF_D50,
          ncol=5, nrow=2, heights=c(10,9),
          common.legend=TRUE, legend="bottom", legend.grob = legend)
dev.off()





# Plot 2 - Coverage plots
SEQUENCE = seq(0.01,0.99,by=0.01)
get_coverage_vector <- function(CLASS,D) {
  if (CLASS=="UNW") {
    UNW = read.csv(paste0("~/Downloads/output",D,".csv"))
    mu <- UNW$estimators_unw
    sd <- UNW$ses_unw
  }
  if (CLASS=="CRF") {
    CRF = read.csv(paste0("~/Downloads/output",D,".csv"))
    mu <- CRF$estimators_crf
    sd <- CRF$ses_crf
  }
  coverage_vec = numeric()
  for (COV in SEQUENCE) {
    coverage_vec = append(coverage_vec,
                          mean((mu-qnorm(1-(1-COV)/2)*sd)*(mu+qnorm(1-(1-COV)/2)*sd)<0)
    )
  }
  return(coverage_vec)
}
get_color <- function(CLASS) {
  if (CLASS=="UNW") color_here <- coul[3]
  if (CLASS=="CRF") color_here <- coul[5]
  return(color_here)
}
plot_coverage <- function(CLASS, D) {
  COV_VEC = get_coverage_vector(CLASS,D)
  COLOR <- get_color(CLASS)
  data <- data.frame(
    group = SEQUENCE,
    mean_value = COV_VEC,
    se = sqrt(COV_VEC*(1-COV_VEC)/1000)
  )
  cov_plot <- ggplot(data, aes(x = group, y = mean_value)) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "solid", color = "black", size = 0.2) +
    geom_line(color=COLOR) +
    geom_ribbon(aes(ymin=mean_value-qnorm(0.975)*se, ymax=mean_value+qnorm(0.975)*se), fill = COLOR, alpha = 0.4) +
    scale_color_manual(values = COLOR) +
    theme(axis.text=element_text(size=4)) +
    theme_minimal(base_size=6) +
    xlab("Theoretical coverage of confidence interval") +
    ylab("Empirical coverage of confidence interval")

  return(cov_plot)
}
COV_UNW_1 <- plot_coverage("UNW",1)
COV_UNW_3 <- plot_coverage("UNW",3)
COV_UNW_5 <- plot_coverage("UNW",5)
COV_UNW_10 <- plot_coverage("UNW",10)
COV_UNW_50 <- plot_coverage("UNW",50)
COV_CRF_1 <- plot_coverage("CRF",1)
COV_CRF_3 <- plot_coverage("CRF",3)
COV_CRF_5 <- plot_coverage("CRF",5)
COV_CRF_10 <- plot_coverage("CRF",10)
COV_CRF_50 <- plot_coverage("CRF",50)
legend_metadata = data.frame(Forest=c(rep("UNW",4),rep("CRF",4)), xval=1:8, yval=1:8)
legend_plot = ggplot(legend_metadata, aes(x=xval, y=yval, group=Forest, color=Forest)) +
  geom_line(size=3.5) +
  theme_bw() +
  scale_color_manual(
    values=c(coul[3],coul[5]),
    labels=c("Classical random forests (RF)","Clustered random forests (CRF)"),
    name="",
    guide = guide_legend(override.aes = list(
      shape = c(16, 16),
      linetype = c("solid","solid")
    ))
  ) +
  theme(legend.direction = "horizontal",
        legend.text = element_text(size=15))
legend = get_legend(legend_plot)

cairo_pdf("~/Downloads/crf-sim2-covplots.pdf", 15, 6)
ggarrange(COV_UNW_1, COV_UNW_3, COV_UNW_5, COV_UNW_10, COV_UNW_50,
          COV_CRF_1, COV_CRF_3, COV_CRF_5, COV_CRF_10, COV_CRF_50,
          ncol=5, nrow=2, heights=c(10,9),
          common.legend=TRUE, legend="bottom", legend.grob = legend)
dev.off()
