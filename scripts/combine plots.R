## utility script to combine various figures

library(ggpubr)


## Fig. 2 ------------------------------------------------------------------

png("./figs/fig_2.png", width=5, height=5.5, units='in', res=300)
ggarrange(ggarrange(fig.2a, fig.2b, ncol=2, widths = c(0.65, 0.35), labels = c("a", "b")),
          fig.2c, nrow=2, heights = c(0.4, 0.6), labels = c("a", "c", "c"))
dev.off()


## combine cod-pollock R plot, pollock seine vs. FAR and pollock model recruitment vs. FAR plots ------------

R.R <- image_read("./figs/predicted_effect_cod_poll_R.png")
FAR <- image_read("./figs/pollock_FAR_recruit_plot.png")

img <- c(R.R, FAR)

stack1 <- image_append(image_scale(img, "100%"))

R.FAR <- image_read("./figs/predicted_effect_pollock_R_FAR.png")
stack2 <- image_append(image_scale(R.FAR, "100%"))
img <- c(stack1, stack2)
stack <- image_append(image_scale(img, "100%"), stack = T)
image_write(stack, path = "figs/pollock_R_FAR_stack.png", format = "png")


## combine modeled CMIP FAR projections and cod-pollock projected R plot ------------------------
plot.nil <- ggplot() + theme_void()

png("./figs/Fig4-projected_FARandR.png", width=5, height=5, units='in', res=300)
ggpubr::ggarrange(CMIP.FAR, 
                  ggpubr::ggarrange(cod.poll.proj.R, plot.nil, widths=c(0.7, 0.3)),
                            ncol=1, heights=c(1,0.9))
dev.off()
