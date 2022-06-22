#############################################################################
#                                                                           #
# DREWNOWSKI'S INDEX TO MEASURE LIFESPAN VARIATION: REVISITING THE GINI     #
#   COEFFICIENT OF THE LIFE TABLE                                           #
#                                                                           #
# Jose Manuel ABURTO, Ugofilippo BASELLINI, Annette BAUDISCH,               #
#   Francisco VILLAVICENCIO - fvillavicencio@ced.uab.es                     #
#                                                                           #
# SET OF FUNCTIONS TO REPRODUCE RESULTS AND FIGURES FROM THE PAPER          #
# Sourced by: Figures.R                                                     #
#                                                                           #
#                                                                 June 2022 #
#############################################################################


############################################################
### MAIN FUNCTIONS                                       ###  
############################################################

# FUNCTION TO CALCULATE DREWNOWSKI'S MEASURES FOR GOMPERTZ MODEL
DrewGomp <- function(a, b, x, dx = 1) {
  
  ## a           Baseline mortality of Gompertz model
  ## b           Mortality rate Gompertz model
  ## x           Age vector
  ## dx          Length of age intervals
  
  # Error message
  if (length(a) > 1 & length(b) > 1) stop('Only one mortality parameter can have length > 1.')
  
  # Number of simulations
  n <- max(length(a), length(b))
  
  # Objects to store estimates
  MX <- DX <- wx <- WX <- matrix(NA, length(x), n)
  D0 <- AD <- rep(NA, n)
  
  for (i in 1:n) {
    
    if (length(a) > 1) {
      # Death rates (mortality hazard)
      MX[, i] <- mxGomp(x, a[i], b)
      # Survival probabilities
      lx <- lxGomp(x, a[i], b)
    } else {
      # Death rates (mortality hazard)
      MX[, i] <- mxGomp(x, a, b[i])
      # Survival probabilities
      lx <- lxGomp(x, a, b[i])
    }
    
    # Distribution of deaths: dx = mx*lx
    DX[, i] <- MX[, i]*lx
    
    # Age-specific remaining life expectancy
    ex <- rev(cumsum(rev(lx*dx))) / lx
    
    # Drewnowski's index conditional on survival to age x 
    if (length(b) > 1) {
      if (b[i] == 0) {
        Dx <- 0.5
      } else {
        # Equation 5
        Dx <- rev(cumsum(rev((lx^2)*dx))) / (rev(cumsum(rev(lx*dx))) * lx)
        Dx[is.nan(Dx)] <- 1
      }
    } else if (b == 0) {
      Dx <- 0.5
    } else {
      # Equation 5
      Dx <- rev(cumsum(rev((lx^2)*dx))) / (rev(cumsum(rev(lx*dx))) * lx)
      Dx[is.nan(Dx)] <- 1
    } 
    D0[i] <- Dx[1]
    
    # Weights w(x)
    wx[, i] <- MX[, i]*lx*ex
    wx[which(is.nan(wx[, i])), i] <- 0
    
    # Weights W(x)
    WX[, i] <- (2*lx*Dx - Dx[1]) / ex[1]
    
    # Threshold age
    AD[i] <- CalcThresh(x = x, y = WX[, i])
    
    # Remove objects
    rm(Dx, ex, lx)
    
  }
  
  # Output
  return(list(AD = AD, D0 = D0, DX = DX, MX = MX, wx = wx, WX = WX))
  
}

# FUNCTION TO CALCULATE THRESHOLD AGE
CalcThresh <- function(x, y, n.grid = 10000) {
  
  ## x         Ages
  ## y         Weights
  ## n.grid    Number of grid points for interpolation
  
  fit <- spline(x = x, y = y, n = n.grid) 
  xi <- fit$x 
  yi <- fit$y 
  xn <- xi[which.min(abs(yi))] 
  
  # Output
  return(xn) 
  
}


# FUNCTION TO CALCULATE RATE OF MORTALITY IMPROVEMENT (Vaupel 1986, Pop Studies)
CalcRho <- function(qx1, qx2, dt){
  
  ## qx1      Death probabilities at time 1
  ## qx2      Death probabilities at time 2
  ## dt       Time interval between time 1 and time 2 
  
  # Estimate rate of improvement
  rho <- (log(-log(1-qx1)) - log(-log(1-qx2))) / dt
  rho[is.infinite(rho)] <- NA
  
  # Remove last age group
  rho <- rho[-length(rho)]
  
  # Output
  return(rho)
  
}


############################################################
### GOMPERTZ FUNCTIONS                                   ###  
############################################################

# MORTALITY HAZARD GOMPERTZ MODEL
mxGomp <- function(x, a, b){
  mx <- a*exp(b*x)
  return(mx)
}

# SURVIVAL FUNCTION GOMPERTZ MODEL
lxGomp <- function(x, a, b) {
  if (b != 0) {
    lx <- exp((a/b)*(1-exp(b*x)))  
  } else lx <- exp(-a*x)
  return(lx)
}


############################################################
### FUNCTION FOR PLOTTING                                ###  
############################################################

# FUNCTION TO PLOT FIGURES 1 TO 5
plotGompertz <- function(x, MX, DX, D0, WX, AD, n1 = 3, n2 = 1, par1, parName, 
                         legend1 = 'topleft', legSplit = F) {
  
  ## x          X-Axis values
  ## MX         Death rates (mortality hazard) (matrix)
  ## DX         Density (matrix)
  ## D0         Drewnowski's index (at age 0)
  ## WX         Weights W(x) (matrix)
  ## AD         Threshold age
  ## n1         Decimal points for Drewnowski's index
  ## n2         Decimal points for threshold age 
  ## par1       Values of parameter that varies
  ## parName    Name of the parameter that varies
  ## legend1    Position of the legend in Panel A)
  ## legSplit   Split legend of Panel A) (TRUE/FALSE) 
  
  
  # Font size of the titles
  cexMain <- 1.4
  cexAxis <- 1.1
  cexLegd <- 1.15
  
  # Line width
  lineswd <- 2
  
  # Time interval x-axis
  if (max(x) < 15) {
    dtPlot <- 2
  } else if (max(x) < 31) {
    dtPlot <- 5
  } else if (max(x) < 61) {
    dtPlot <- 10
  } else if (max(x) < 201) {
    dtPlot <- 20
  } else if (max(x) < 601) {
    dtPlot <- 50
  } else if (max(x) < 1001) {
    dtPlot <- 100
  } else if (max(x) < 1500) {
    dtPlot <- 200
  } else dtPlot <- 500
  
  # LAYOUT
  layout(mat = rbind(c(1, 2, 3), c(0, 4, 0)),
         widths = rep(1, 3), heights = c(1, .08))
  
  
  #----------#
  # PANEL A) #
  #----------#
  
  # Empty plot
  par(las = 1, mar = c(2, 3, 2, 1))
  matplot(x, log(MX), t = 'n', xaxt = 'n', xlab = '', ylab = '', 
          bty = 'L', cex.axis = cexAxis)
  
  # Title
  title(main = substitute(paste(italic('A) Mortality hazard (log)'))), 
        cex.main = cexMain, line = .9, adj = 0.1)
  
  # Lines
  matlines(x, log(MX), lwd = lineswd, col = my.col, lty = 1)
  
  # X-axis
  axis(1, at = seq(min(x), max(x), dtPlot), cex.axis = cexAxis)
  
  # Legend
  if (!is.null(legend1)) {
    legNames <- c()
    if (parName == 'alpha') {
      for (i in 1:n) {
        legNames <- c(legNames, 
                      as.expression(bquote(alpha == .(format(round(par1[i], 4), 
                                                             nsmall = 4)))))
      } 
    } else {
      for (i in 1:n) {
        legNames <- c(legNames, 
                      as.expression(bquote(beta == .(format(round(par1[i], 3), 
                                                            nsmall = 3)))))
      }
    }
    if (legSplit) {
      legend(legend1, legNames[1:5],
             cex = cexLegd, text.col = my.col[1:5], bty = 'n')
      legend2 <- ifelse(legend1 == 'topleft', 'bottomright', 'bottomleft')
      legend(legend2, legNames[6:n],
             cex = cexLegd, text.col = my.col[6:n], bty = 'n')
    } else legend(legend1, legNames,
                  cex = cexLegd, text.col = my.col, bty = 'n')  
  } else {
    for (i in 1:n){
      text(x = max(x)*.83, y = log(MX[1, i])*.997, 
           bquote(alpha == .(format(round(par1[i], 4), nsmall = 4))),
           cex = cexLegd, col = my.col[i], adj = c(0.5, 0))
    }
    
  }
  
  
  #----------#
  # PANEL B) #
  #----------#
  
  # Empty plot
  par(las = 1, mar = c(2, 3.25, 2, 0.75))
  matplot(x, DX, t = 'n', xaxt = 'n', xlab = '', ylab = '', 
          bty = 'L', cex.axis = cexAxis)
  
  # Title
  title(main = substitute(paste(italic("B) Distribution of ages at death"))), 
        cex.main = cexMain, line = .9, adj = 0.5)
  
  # Lines
  matlines(x, DX, lwd = lineswd, col = my.col, lty = 1)
  
  # X-axis
  axis(1, at = seq(min(x), max(x), dtPlot), cex.axis = cexAxis)
  
  # Legend
  legend('topright', paste("D =", format(round(D0, n1), nsmall = n1)),
         cex = cexLegd, text.col = my.col, bty = 'n')
  
  
  #----------#
  # PANEL C) #
  #----------#
  
  # Empty plot
  par(las = 1, mar = c(2, 3.5, 2, .5), xpd = F)
  matplot(x, WX, t = 'n', xaxt = 'n', ylab = '', xlab = '', 
          bty = 'L', cex.axis = cexAxis)
  
  # Title
  title(main = substitute(paste(italic("C) Weights W(x)"))),
        cex.main = cexMain, line = .9, adj = 0.03)
  
  # Horizontal line
  abline(h = 0)
  
  # Vertical lines
  abline(v = AD, lty = 2, col = my.col, lwd = 0.65)
  
  # Lines
  matlines(x, WX, lwd = lineswd, col = my.col, lty = 1)
  
  # X-axis
  axis(1, at = seq(min(x), max(x), dtPlot), cex.axis = cexAxis)
  
  # Legend
  legNames <- c()
  for (i in 1:n) legNames <- c(legNames, 
                               as.expression(bquote(a^D == .(format(round(AD[i], n2), 
                                                                    nsmall = n2)))))
  legend('topright', legend = legNames, cex = cexLegd,
         text.col = my.col, bty = 'n')

  
  #--------#
  # X-AXIS #
  #--------#
  
  # X-AXIS
  par(mar = rep(0, 4))
  plot(c(0, 1), c(0, 1), xlab = "", ylab = "", col = NA, axes = F)
  text(0.5, 0.5, 'Age (years)', cex = 1.8)
  
}

