#############################################################################
#                                                                           #
# DREWNOWSKI'S INDEX TO MEASURE LIFESPAN VARIATION: REVISITING THE GINI     #
#   COEFFICIENT OF THE LIFE TABLE                                           #
#                                                                           #
# Jose Manuel ABURTO, Ugofilippo BASELLINI, Annette BAUDISCH,               #
#   Francisco VILLAVICENCIO - fvillavicencio@ced.uab.es                     #
#                                                                           #
# CODE TO REPRODUCE RESULTS AND FIGURES FROM THE PAPER                      #
# Requires: Functions.R                                                     #
#                                                                           #
#                                                                 June 2022 #
#############################################################################


############################################################
### PACKAGES AND FUNCTIONS                               ###  
############################################################

# Clear workspace
rm(list = ls())

# Prevent scientific notation
options(scipen = 999)

# Packages
library(viridisLite)

# Functions
source('Functions.R')


############################################################
### PLOTTING PARAMETERS                                  ###  
############################################################

# Number of scenarios
n <- 9

# Color palette
my.col <- c(viridis(n = ceiling(n/2), begin = .2, end = .8, option = 'magma'),
            viridis(n = floor(n/2), begin = .4, option = 'viridis'))


############################################################
### FIGURE 1. CHANGING BASELINE MORTALITY                ###  
############################################################

# Age interval
dx <- .1

# Age range 
x <- seq(0, 80, dx)

# Gompertz parameters: Changing starting level of mortality (a)
a <- exp(seq(-2, -8, length.out = n))
b <- 0.1

# Drewnowski measures for Gompertz model
out <- DrewGomp(a, b, x, dx)

# Drewnowski's index
out$D0

# Threshold ages
out$AD

# Plot
PDF <- F
if (PDF) pdf("Fig1.pdf", width = 8, height = 5)
plotGompertz(x, MX = out$MX, DX = out$DX, D0 = out$D0, WX = out$WX, AD = out$AD,
             par1 = a, parName = 'alpha', legSplit = T)
if (PDF) dev.off()

# Remove objects
rm(a, b, dx, out, PDF, x)


############################################################
### FIGURE 2. CHANGING RATE OF AGING                     ###  
############################################################

# Age interval
dx <- .1

# Age range 
x <- seq(0, 80, dx)

# Gompertz parameters: Changing rate of ageing (b)
a <- exp(-6)
b <- seq(0.25, 0.075, length.out = n)

# Drewnowski measures for Gompertz model
out <- DrewGomp(a, b, x, dx)

# Drewnowski's index
out$D0

# Threshold ages
out$AD

# Plot
PDF <- F
if (PDF) pdf("Fig2.pdf", width = 8, height = 5)
plotGompertz(x, MX = out$MX, DX = out$DX, D0 = out$D0, WX = out$WX, AD = out$AD,
             par1 = b, parName = 'beta')
if (PDF) dev.off()

# Remove objects
rm(a, b, dx, out, PDF, x)


############################################################
### FIGURE 3. CONSTANT MORTALITY                         ###  
############################################################

# Age interval
dx <- 1

# Age range 
x <- seq(0, 1500, dx)

# Gompertz parameters: Constant mortality (b = 0), changing value of a
a <- seq(5e-3, 1e-3, length.out = n)
b <- 0

# Drewnowski measures for Gompertz model
out <- DrewGomp(a, b, x, dx)

# Drewnowski's index
out$D0

# Threshold ages
out$AD

# Plot
PDF <- F
if (PDF) pdf("Fig3.pdf", width = 8, height = 5)
plotGompertz(x, MX = out$MX, DX = out$DX, D0 = out$D0, WX = out$WX, AD = out$AD,
             n1 = 1, n2 = 0, par1 = a, legend1 = NULL)
if (PDF) dev.off()

# Remove objects
rm(a, b, dx, out, PDF, x)


############################################################
### FIGURE 4. DECREASING MORTALITY WITH AGE              ###  
############################################################

# Age interval
dx <- .01

# Age range 
x <- seq(0, 80, dx)
xPlot <- 1:length(seq(0, 10, dx))

# Gompertz parameters: Decreasing mortality (b < 0), changing value of a
a <- exp(seq(-0.05, -1, length.out = n))
b <- -0.05

# Drewnowski measures for Gompertz model
out <- DrewGomp(a, b, x, dx)

# Drewnowski's index
out$D0

# Threshold ages
out$AD

# Plot
PDF <- F
if (PDF) pdf("Fig4.pdf", width = 8, height = 5)
plotGompertz(x[xPlot], MX = out$MX[xPlot, ], DX = out$DX[xPlot, ], 
             D0 = out$D0, WX = out$WX[xPlot, ], AD = out$AD, n2 = 2,
             par1 = a, parName = 'alpha', legend1 = 'topright', legSplit = T)
if (PDF) dev.off()

# Remove objects
rm(a, b, dx, out, PDF, x, xPlot)


############################################################
### FIGURE 5. DECREASING MORTALITY WITH AGE              ###  
############################################################

# Age interval
dx <- .01

# Age range 
x <- seq(0, 80, dx)
xPlot <- 1:length(seq(0, 10, dx))

# Decreasing mortality: changing value of b
a <- exp(-0.5)
b <- -seq(0.07, 0.1125, length.out = n)

# Drewnowski measures for Gompertz model
out <- DrewGomp(a, b, x, dx)

# Drewnowski's Index
out$D0

# Threshold ages
out$AD

# Plot
PDF <- F
if (PDF) pdf("Fig5.pdf", width = 8, height = 5)
plotGompertz(x[xPlot], MX = out$MX[xPlot, ], DX = out$DX[xPlot, ], 
             D0 = out$D0, WX = out$WX[xPlot, ], AD = out$AD, n2 = 2,
             par1 = b, parName = 'beta', legend1 = 'topright')
if (PDF) dev.off()

# Remove objects
rm(a, b, dx, out, PDF, x, xPlot)


############################################################
### FIGURE 6. RATE OF MORTALITY IMPROVEMENT IN SWEDEN    ###  
############################################################

# DATA: LIFE TABLES Sweden 1990-2021, 5-year age groups, from Human Mortality Database 
#  https://www.mortality.org/ (downloaded on 11 June 2022)
dat <- read.table(file = 'SwedenFemLT_5x1.txt', header = T, skip = 2)

# Age intervals
dx <- c(1, 4, rep(5, length(unique(dat$Age)) - 2))

# Ages (beginning of age interval)
x <- c(0, cumsum(dx))
x <- x[-length(x)]

# Age groups
age <- unique(dat$Age)
age <- age[-length(age)]

# Select two focus years
times <- c(1980, 2000)

# Select two time intervals for improvement
dt <- c(5, 20)

# Years of interest
Years <- sort(c(rep(times, length(dt)) - rep(dt, each = length(times)), 
                times[length(times)]))
dat <- dat[dat$Year %in% Years, ]

# Objects to store estimates
e0 <- AD <- D0 <- rep(NA, length(Years))
names(e0) <- names(AD) <- names(D0) <- c(Years)

# Drewnowski measures
for (i in 1:length(Years)) {
  # Life expectancy at birth
  e0[i] <- dat$ex[dat$Year == Years[i]][1]
  # Survival
  lx <- dat$lx[dat$Year == Years[i]]
  lx <- lx / lx[1]
  # Drewnowski's index (Equation 5)
  Dx <- rev(cumsum(rev((lx^2)*dx))) / (rev(cumsum(rev(lx*dx))) * lx)
  D0[i] <- Dx[1]
  # Weights W(x)
  Wx <- (2*lx*Dx - D0[i]) / e0[i]
  # Threshold age
  AD[i] <- CalcThresh(x = x, y = Wx)
}
e0; round(AD, 2); round(D0, 3)


#------#
# PLOT #
#------#

# Line width
lineswd <- 1.5

# Palette
colPalette <- my.col[c(5, 7)]

# Output file
PDF <- F
if (PDF) pdf("Fig6.pdf", width = 7.5, height = 6.5)

# Layout
layout(mat = rbind(c(2, 4, 3), c(1, 1, 1)),
       widths = c(1, .2, 1), heights = c(1, .06))

# X-AXIS
par(mar = rep(0, 4))
plot(c(0, 1), c(0, 1), xlab = "", ylab = "", col = NA, axes = F)
text(0.5, 0.4, expression(paste('Rate of mortality improvement', ~rho(x))), 
     cex = 1.6)

# X-axis limit
xLim <- c(-.05, .10)

# PANELS A AND B
for (i in dt) {

  # Starting year in the plot
  iniYear <- times - i
  
  # Empty plot
  if (i == dt[1]) {
    par(las = 1, mar = c(2, 1, 2.5, 0))
  } else par(las = 1, mar = c(2, 0, 2.5, 1))
  plot(NA, xlim = xLim*1.08, ylim = c(0, length(age) - 2), axes = F)
  
  # Grid
  abline(v = seq(xLim[1], xLim[2], .05), col = grey(.85), lwd = .9)
  abline(h = 1:length(age) - 1, col = grey(.85), lwd = .9)
  
  # Re-scale threshold ages for plotting
  thage <- 1 + (AD[paste(iniYear)] - 5) / 5
  
  # Threshold ages
  abline(h = thage, col = colPalette, lwd = lineswd*.8, lty = 2)
  text(x = xLim[1] + 0.002, y = thage + c(-.4, .55),
       c(as.expression(bquote(.(iniYear[1])*':'~a^D == .(format(round(AD[paste(iniYear[1])], 2), 
                                                                nsmall = 2)))),
         as.expression(bquote(.(iniYear[2])*':'~a^D == .(format(round(AD[paste(iniYear[2])], 2), 
                                                              nsmall = 2))))),
       col = colPalette, adj = c(0, .5))
  
  # Add to existing plot
  par(new = TRUE)
  plot(NA, xlim = xLim*1.08, ylim = c(0, length(age) - 2),
       xlab = '', ylab = '', xaxt = 'n', yaxt = 'n',  bty = 'U')
  
  # Horizontal line
  abline(v = 0)
  
  # Axis
  axis(1, at = seq(xLim[1], xLim[2], .05), cex.axis = 1.08)
  if (i == dt[1]) {
    axis(4, at = 1:length(age) - 1, labels = F)
  } else axis(2, at = 1:length(age) - 1, labels = F)
  
  # Add rates of improvement
  for (j in 1:length(times)) {
    rho <- CalcRho(qx1 = dat$qx[dat$Year == times[j] - i],
                   qx2 = dat$qx[dat$Year == times[j]], dt = i)
    lines(rho, 1:length(age) - 1.5, col = colPalette[j], lwd = lineswd)
    points(rho, 1:length(age) - 1.5, col = colPalette[j], lwd = lineswd, pch = 19)
  }
  
  # Title
  if (i == dt[1]) {
    title(main = substitute(paste(italic('A) 5-year time intervals'))),
          cex.main = 1.4, line = .95, adj = 0.05)
  } else title(main = substitute(paste(italic('B) 20-year time intervals'))),
               cex.main = 1.4, line = .95, adj = 0.05)
  
  # Legend
  legend(x = .075, y = length(age) - 3.5, xjust = .5, yjust = .5,
         legend = c(paste(iniYear[1], '-', times[1]), paste(iniYear[2], '-', times[2])), 
         col = colPalette, lwd = .8*lineswd, pch = 19, 
         y.intersp = 1.4, seg.len = 2.5, bg = 'white')
    
}

# Label age groups
par(mar = c(2, 0, 2.5, 0))
plot(c(0, 1), c(0, length(age) - 2), xlab = "", ylab = "", col = NA, axes = F)
for (i in 1:length(age)) text(0.5, i - 1.5, age[i], adj = .5, cex = 1.08)

# Close device
if (PDF) dev.off()

# Remove objects
rm(AD, age, colPalette, dat, D0, dt, dx, Dx, e0, i, iniYear, j, lineswd, lx, 
   PDF, rho, thage, times, Wx, x, xLim, Years)

