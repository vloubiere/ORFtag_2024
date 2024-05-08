setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# Import data ----
dat <- readRDS("/groups/stark/nemcko/work/FACS/20240123_Dox2/raw.data.RDS")
neg <- dat[name=="input"]
negds <- density(neg$GFP.A) # Negative control
negds$y <- negds$y/max(negds$y, na.rm = T)
dat <- dat[grepl("day", name)]
dat <- dat[!grepl("day6", name)]
dat[, day:= tstrsplit(name, "_", keep= 1), name]
dat[, DOX:= tstrsplit(name, "_", keep= 3)=="DOX", name]
dat[, rep:= tstrsplit(name, "_", keep= 4), name]
dat[day=="day0" & rep=="rep1", DOX:= T]
dat <- dat[!(day!="day0" & rep=="rep2")]

# Plot ----
pdf("pdf/review_FACS_Dox_induced.pdf", 4, 4)
vl_par(mfrow= c(6, 1),
       mai= c(0.05,0,0,0),
       omi= c(.9,1,.9,1),
       cex.axis= 0.7,
       mgp= c(1.5, 0.2, 0),
       cex.lab= 1,
       tck= -0.05)
dat[, {
  # xlim
  xl <- c(6, 18)
  # Density
  nd <- density(GFP.A[(!DOX)], from= xl[1], to= xl[2])
  d <- density(GFP.A[(DOX)], from= xl[1], to= xl[2])
  nd$y <- nd$y/max(nd$y, na.rm = T)
  d$y <- d$y/max(d$y, na.rm = T)
  # Compute ylim
  yl <- c(0, max(c(nd$y, d$y)))
  # Colors <- 
  Cc <- adjustcolor(c("lightgrey", "cornflowerblue", "tomato"), .3)
  # Plot
  plot(NA,
       xlim= xl,
       ylim= yl,
       xaxt= "n",
       yaxt= "n",
       type= "n")
  legend(par("usr")[1]-strwidth("M")*1.25,
         par("usr")[4]+strheight("M"),
         day,
         bty= "n")
  at <- axisTicks(yl, nint = 3, log= F)
  at <- at[c(1, length(at))]
  axis(2,
       at = at,
       lwd= .5,
       line = -.2)
  if(.GRP==1)
    legend(par("usr")[2],
           par("usr")[4],
           fill= Cc,
           legend = c("Control", "noDOX", "DOX"),
           bty= "n",
           xpd= NA)
  if(.GRP==.NGRP)
  {
    par(mgp= c(1.5, 0, 0))
    axis(1,
         lwd= .5,
         line= .2)
    title(ylab= "Scaled density(x)",
          outer= T)
    title(xlab= "GFP intensity",
          outer= T,
          line= .75)
  }
  polygon(c(xl[1], negds$x, xl[2]),
          c(0, negds$y, 0),
          col= Cc[1],
          lwd= 0.5)
  polygon(c(xl[1], nd$x, xl[2]),
          c(0, nd$y, 0),
          col= Cc[2],
          lwd= 0.5)
  if(any(DOX))
    polygon(c(xl[1], d$x, xl[2]),
            c(0, d$y, 0),
            col= Cc[3],
            lwd= 0.5)
}, day]
dev.off()