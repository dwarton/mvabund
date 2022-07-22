# Code to create mvabund hex sticker

windFarms <- ecostats::windFarms
ft_wind=manyglm(mvabund(windFarms$abund)~Year*Zone, family="poisson", data=windFarms$X)

o <- par()
par(bg=NA)
plot(ft_wind)

dev.copy(png,'inst/figures/manyglmplot.png')
dev.off()

#I exported the .png and used photoshop to make it white and removed text, couldn't do this in R efficiently

imgurl <- file.path("inst/figures/manyglmplot_white.png")
sticker(imgurl, package="mvabund",
        p_size=18, p_y = 1.55,
        s_x=1, s_y=.88, s_width=.6, s_height = 0.6,
        h_fill = "#273046",
        h_color = "#CCC591",
        filename="inst/figures/mvabund_hex.png")

str(wesanderson::wes_palette("BottleRocket2"))

