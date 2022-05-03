# PDE 2D ABC functions adjusted (time-dependent parameters)
# Author: Yunchen Xiao

# This .R file contains all the necessary functions that carry out the 
# time-dependent simulations for the SCC patterns. 

# Read in the cell densities at different sections of the domain
# at the specified time points. 

# Load the necessary packages.
library(readr)

# Read in the reference SCC densities.
scc.ref.den <- read_rds("Reference densities SCC.rds")
t3.ref.den <- scc.ref.den$t3.ref.den
t6.ref.den <- scc.ref.den$t6.ref.den
t9.ref.den <- scc.ref.den$t9.ref.den
t12.ref.den <- scc.ref.den$t12.ref.den
t14.ref.den <- scc.ref.den$t14.ref.den

# Read in the regression model coefficients
# Regression coefficients
dn.regression.model <- read_rds("dn quadratic regression model.rds")
dn.lin.coef <- unname(dn.regression.model$coefficients[2])
dn.quad.coef <- unname(dn.regression.model$coefficients[3])

gamma.regression.model <- read_rds("gamma quadratic regression model.rds")
gamma.lin.coef <- unname(gamma.regression.model$coefficients[2])
gamma.quad.coef <- unname(gamma.regression.model$coefficients[3])

rn.regression.model <- read_rds("rn quadratic regression model.rds")
rn.lin.coef <- unname(rn.regression.model$coefficients[2])
rn.quad.coef <- unname(rn.regression.model$coefficients[3])

eta.regression.model <- read_rds("eta quadratic regression model.rds")
eta.lin.coef <- unname(eta.regression.model$coefficients[2])
eta.quad.coef <- unname(eta.regression.model$coefficients[3])

alpha.regression.model <- read_rds("alpha quadratic regression model.rds")
alpha.lin.coef <- unname(alpha.regression.model$coefficients[2])
alpha.quad.coef <- unname(alpha.regression.model$coefficients[3])

prob.prof.regression.model <- read_rds("prob.prof quadratic regression model.rds")
prob.prof.lin.coef <- unname(prob.prof.regression.model$coefficients[2])
prob.prof.quad.coef <- unname(prob.prof.regression.model$coefficients[3])

# Functions that generates the time-varying parameters. 
dn.varying.paras <- function(dn.init.par, dn.lin.temp, dn.quad.temp) {
  # Purpose: after the initial parameters values are sampled from the prior, 
  # we take them as the parameter values in the first period (day 1 - 3). This 
  # function takes the initial parameter estimates, regression model 
  # coefficients as the arguments and computes the parameter values in the 4
  # later periods. It then checks if all the time-varying parameter values
  # are supported by the corresponding prior distribution. 
  
  # Based on the idea of quadratic formula, assume the value for the first 
  # period, the linear and quadratic coefficients are known. Compute the values
  # for the 4 later periods. 
  dn.p2 <- dn.init.par + 2 * dn.quad.temp * 1 + dn.quad.temp * 1 + dn.lin.temp
  dn.p3 <- dn.init.par + 2 * dn.quad.temp * 2 + dn.quad.temp * (2^2) + 
    dn.lin.temp * 2
  dn.p4 <- dn.init.par + 2 * dn.quad.temp * 3 + dn.quad.temp * (3^2) + 
    dn.lin.temp * 3
  dn.p5 <- dn.init.par + 2 * dn.quad.temp * 4 + dn.quad.temp * (4^2) + 
    dn.lin.temp * 4
  
  # Combine the parameter values together. 
  dn.full.mat <- cbind(dn.init.par, dn.p2, dn.p3, dn.p4, dn.p5, dn.quad.temp,
                       dn.lin.temp)
  
  # Check if there are any lies beyond the support of the corresponding prior. 
  dn.full.mat[dn.full.mat[,1:5] < 0.000069] <- NA
  dn.full.mat[dn.full.mat[,1:5] > 0.02] <- NA
  
  # Remove the parameter sets that are not supported by the corresponding prior.
  dn.full.mat.invalid <- unique(which(is.na(dn.full.mat), arr.ind = TRUE)[,1])
  if (length(dn.full.mat.invalid) > 0) {
    dn.full.mat.valid <- dn.full.mat[-dn.full.mat.invalid,]
  } else {
    dn.full.mat.valid <- dn.full.mat
  }
  
  return(dn.full.mat.valid)
}

gamma.varying.paras <- function(gamma.init.par, gamma.lin.temp, gamma.quad.temp) {
  gamma.p2 <- gamma.init.par + 2 * gamma.quad.temp * 1 + gamma.quad.temp * 1 + 
    gamma.lin.temp
  gamma.p3 <- gamma.init.par + 2 * gamma.quad.temp * 2 + 
    gamma.quad.temp * (2^2) + gamma.lin.temp * 2
  gamma.p4 <- gamma.init.par + 2 * gamma.quad.temp * 3 + 
    gamma.quad.temp * (3^2) + gamma.lin.temp * 3
  gamma.p5 <- gamma.init.par + 2 * gamma.quad.temp * 4 + 
    gamma.quad.temp * (4^2) + gamma.lin.temp * 4
  
  gamma.full.mat <- cbind(gamma.init.par, gamma.p2, gamma.p3, gamma.p4, 
                          gamma.p5, gamma.quad.temp, gamma.lin.temp)
  
  gamma.full.mat[gamma.full.mat[,1:5] < 0.005] <- NA
  gamma.full.mat[gamma.full.mat[,1:5] > 0.26] <- NA
  
  gamma.full.mat.invalid <- unique(which(is.na(gamma.full.mat), 
                                         arr.ind = TRUE)[,1])
  
  if (length(gamma.full.mat.invalid) > 0) {
    gamma.full.mat.valid <- gamma.full.mat[-gamma.full.mat.invalid,]
  } else {
    gamma.full.mat.valid <- gamma.full.mat
  }
  
  return(gamma.full.mat.valid)
}

rn.varying.paras <- function(rn.init.par, rn.lin.temp, rn.quad.temp) {
  rn.p2 <- rn.init.par + 2 * rn.quad.temp * 1 + rn.quad.temp * 1 + rn.lin.temp
  rn.p3 <- rn.init.par + 2 * rn.quad.temp * 2 + rn.quad.temp * (2^2) + 
    rn.lin.temp * 2
  rn.p4 <- rn.init.par + 2 * rn.quad.temp * 3 + rn.quad.temp * (3^2) + 
    rn.lin.temp * 3
  rn.p5 <- rn.init.par + 2 * rn.quad.temp * 4 + rn.quad.temp * (4^2) + 
    rn.lin.temp * 4
  rn.full.mat <- cbind(rn.init.par, rn.p2, rn.p3, rn.p4, rn.p5, rn.quad.temp, 
                   rn.lin.temp)
  
  rn.full.mat[rn.full.mat[,1:5] < 0.0008] <- NA
  rn.full.mat[rn.full.mat[,1:5] > 0.08] <- NA
  
  rn.full.mat.invalid <- unique(which(is.na(rn.full.mat), 
                                         arr.ind = TRUE)[,1])
  
  if (length(rn.full.mat.invalid) > 0) {
    rn.full.mat.valid <- rn.full.mat[-rn.full.mat.invalid,]
  } else {
    rn.full.mat.valid <- rn.full.mat
  }

  return(rn.full.mat.valid)
}

eta.varying.paras <- function(eta.init.par, eta.lin.temp, eta.quad.temp) {
  eta.p2 <- eta.init.par + 2 * eta.quad.temp * 1 + eta.quad.temp * 1 + 
    eta.lin.temp
  eta.p3 <- eta.init.par + 2 * eta.quad.temp * 2 + eta.quad.temp * (2^2) + 
    eta.lin.temp * 2
  eta.p4 <- eta.init.par + 2 * eta.quad.temp * 3 + eta.quad.temp * (3^2) + 
    eta.lin.temp * 3
  eta.p5 <- eta.init.par + 2 * eta.quad.temp * 4 + eta.quad.temp * (4^2) + 
    eta.lin.temp * 4
  eta.full.mat <- cbind(eta.init.par, eta.p2, eta.p3, eta.p4, eta.p5, eta.quad.temp,
                    eta.lin.temp)
  
  eta.full.mat[eta.full.mat[,1:5] < 7] <- NA
  eta.full.mat[eta.full.mat[,1:5] > 18] <- NA
  
  eta.full.mat.invalid <- unique(which(is.na(eta.full.mat), 
                                      arr.ind = TRUE)[,1])
  
  if (length(eta.full.mat.invalid) > 0) {
    eta.full.mat.valid <- eta.full.mat[-eta.full.mat.invalid,]
  } else {
    eta.full.mat.valid <- eta.full.mat
  }
  
  return(eta.full.mat.valid)
}

alpha.varying.paras <- function(alpha.init.par, alpha.lin.temp, alpha.quad.temp) {
  alpha.p2 <- alpha.init.par + 2 * alpha.quad.temp * 1 + alpha.quad.temp * 1 + 
    alpha.lin.temp
  alpha.p3 <- alpha.init.par + 2 * alpha.quad.temp * 2 + 
    alpha.quad.temp * (2^2) + alpha.lin.temp * 2
  alpha.p4 <- alpha.init.par + 2 * alpha.quad.temp * 3 + 
    alpha.quad.temp * (3^2) + alpha.lin.temp * 3
  alpha.p5 <- alpha.init.par + 2 * alpha.quad.temp * 4 + 
    alpha.quad.temp * (4^2) + alpha.lin.temp * 4
  alpha.full.mat <- cbind(alpha.init.par, alpha.p2, alpha.p3, alpha.p4, alpha.p5,
                      alpha.quad.temp, alpha.lin.temp)
  
  alpha.full.mat[alpha.full.mat[,1:5] < 0.07] <- NA
  alpha.full.mat[alpha.full.mat[,1:5] > 0.18] <- NA
  
  alpha.full.mat.invalid <- unique(which(is.na(alpha.full.mat), 
                                       arr.ind = TRUE)[,1])
  
  if (length(alpha.full.mat.invalid) > 0) {
    alpha.full.mat.valid <- alpha.full.mat[-alpha.full.mat.invalid,]
  } else {
    alpha.full.mat.valid <- alpha.full.mat
  }
  
  return(alpha.full.mat.valid)
}

prob.prof.varying.paras <- function(prob.prof.init.par, prob.prof.lin.temp, 
                                    prob.prof.quad.temp) {
  prob.prof.p2 <- prob.prof.init.par + 2 * prob.prof.quad.temp * 1 + 
    prob.prof.quad.temp * 1 + prob.prof.lin.temp
  prob.prof.p3 <- prob.prof.init.par + 2 * prob.prof.quad.temp * 2 + 
    prob.prof.quad.temp * (2^2) + prob.prof.lin.temp * 2
  prob.prof.p4 <- prob.prof.init.par + 2 * prob.prof.quad.temp * 3 + 
    prob.prof.quad.temp * (3^2) + prob.prof.lin.temp * 3
  prob.prof.p5 <- prob.prof.init.par + 2 * prob.prof.quad.temp * 4 + 
    prob.prof.quad.temp * (4^2) + prob.prof.lin.temp * 4
  prob.prof.full.mat <- cbind(prob.prof.init.par, prob.prof.p2, prob.prof.p3, 
                          prob.prof.p4, prob.prof.p5, prob.prof.quad.temp,
                          prob.prof.lin.temp)
  
  prob.prof.full.mat[prob.prof.full.mat[,1:5] < 0.2] <- NA
  prob.prof.full.mat[prob.prof.full.mat[,1:5] > 1] <- NA
  
  prob.prof.full.mat.invalid <- unique(which(is.na(prob.prof.full.mat), 
                                         arr.ind = TRUE)[,1])
  
  if (length(prob.prof.full.mat.invalid) > 0) {
    prob.prof.full.mat.valid <- prob.prof.full.mat[-prob.prof.full.mat.invalid,]
  } else {
    prob.prof.full.mat.valid <- prob.prof.full.mat
  }
  
  return(prob.prof.full.mat.valid)
}


################################################################################
generate.pattern <- function(par, init.cells.cols, prob.death, prob.prof, slopes) {
  # Purpose: Generate SCC invasion pattern with the given
  # parameters.
  
  # Arguments:
  #    par: parameters of the PDE model
  #    init.cells.cols: initial columns of cells being set at the left boundary 
  #                     of the domain.
  #    prob.death: proportion of cells that will be undergo extinction at the 
  #    end of every day.
  #    prob.prof: proportion of cells that will undergo mitosis at the end 
  #    of every day.
  #    slopes: regression coefficients that are used to generate time-varying
  #    parameters. 
  
  # Space discretization 
  h <- 1/47
  space.length.y <- (1/h) + 1
  space.length.x <- round(space.length.y * (280/480))
  x <- seq(0, h * (space.length.x - 1), length.out = space.length.x)
  y <- seq(0, 1, length.out = space.length.y)

  # Time discretization
  T <- 21.1
  dt <- 0.004
  timesteps <- round(T/dt)
  day.timesteps <- 375
  
  # Parameters of PDE
  dn <- par[1]
  gamma <- par[2]
  r <- par[3]
  eta <- par[4]
  dm <- par[5]
  alpha <- par[6]
  beta <- 0
  
  # First check if the variation of parameter values go beyond the bounds
  dn.5.paras <- dn.varying.paras(dn.init.par = dn, 
                                 dn.quad.temp = slopes[1], 
                                 dn.lin.temp = slopes[2])
  
  gamma.5.paras <- gamma.varying.paras(gamma.init.par = gamma, 
                                       gamma.quad.temp = slopes[3],
                                       gamma.lin.temp = slopes[4])
  
  r.5.paras <- rn.varying.paras(rn.init.par = r, 
                                 rn.quad.temp = slopes[5],
                                 rn.lin.temp = slopes[6])
  
  eta.5.paras <- eta.varying.paras(eta.init.par = eta, 
                                   eta.quad.temp = slopes[7],
                                   eta.lin.temp = slopes[8])
  
  alpha.5.paras <- alpha.varying.paras(alpha.init.par = alpha, 
                                       alpha.quad.temp = slopes[9],
                                       alpha.lin.temp = slopes[10])
  
  prob.prof.5.paras <- prob.prof.varying.paras(prob.prof.init.par = prob.prof, 
                                         prob.prof.quad.temp = slopes[11],
                                         prob.prof.lin.temp = slopes[12])
  
  paras.varying.mat <- rbind(dn.5.paras[1:5], gamma.5.paras[1:5], r.5.paras[1:5],
                             eta.5.paras[1:5], alpha.5.paras[1:5], prob.prof.5.paras[1:5])
  
  # If any of the values is non-valid, break it and return NaNs. 
  if (any(is.na(paras.varying.mat))) {
      return(list(diff.tot = NaN, err.mess = "Invalid time-varying parameter"))
      break
  }
  
  # Initial condition
  n0 <- matrix(0, nrow = space.length.y, ncol = space.length.x)
  
  for (i in 1:space.length.x) {
    if (x[i] <= 0.1) {
      n0[1, i] <- cos(pi * x[i] * 5)
    } else {
      n0[1, i] <- 0
    }
  }
  
  for (i in 2:space.length.y) {
    n0[i, ] <- n0[1, ]
  }
  
  f0 <-  1 - 0.5 * n0
  
  m0 <-  0.5 * n0
  
  n <- n0
  f <- f0
  m <- m0
  
  # Sort the initial cells
  n0.sort <- n0
  
  n.cells <- round(space.length.y * init.cells.cols)
  x.coord <- vector()
  y.coord <- vector()
  
  # Initial cancer cells at the locations with the highest densities in
  # the domain.
  while (length(x.coord) < n.cells) {
    ind <- which(n0.sort == max(n0.sort), arr.ind = TRUE)
    if (length(ind[,1]) > 1) {
      ind.unique <- sample(seq(1, length(ind[,1]), by = 1), 1 ,
                           replace = FALSE)
      ind <- ind[ind.unique,]
      x.coord <- c(x.coord, ind[1])
      y.coord <- c(y.coord, ind[2])
      n0.sort[ind[1], ind[2]] <- -Inf
    } else {
      x.coord <- c(x.coord, ind[1])
      y.coord <- c(y.coord, ind[2])
      n0.sort[ind[1], ind[2]] <- -Inf
    }
  }
  
  ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
  
  for (z in 1:length(x.coord)) {
    ind.position[x.coord[z], y.coord[z]] <- 1
  }
  
  # Set the cut points for the domain (to be used later for density matching &
  # discrepancy calculation.)
  mat.size <- space.length.y/12
  y.cut <- seq(1, space.length.y, by = mat.size)
  x.cut <- seq(1, space.length.x, by = mat.size)
  
  # Numerical scheme which solves the PDE system
  for (i in 1:timesteps) {
    
    # At the end of everyday, some of the current cells
    # in the domain will undergo extinction or mitosis.
    
    if (i %% day.timesteps == 0) {
      # Extract the density values at the positions of current cells
      cell.den <- rep(0, length(x.coord))
      for (u in 1:length(cell.den)) {
        cell.den[u] <- n[x.coord[u], y.coord[u]]      
      }
      
      # Cell extinction: some of the current cells at the locations with the 
      # lowest cell densities will undergo extinction.
      
      dead.cells.num <- round(length(x.coord) * prob.death)
      dead.cells <- rep(0, dead.cells.num)
      
      for (uu in 1:length(dead.cells)) {
        # Find the cells at locations with the lowest densities.
        ind.dead <- which(cell.den == min(cell.den))
        if(length(ind.dead) > 1) {
          ind.dead <- sample(ind.dead, 1, replace = FALSE)
        }
        dead.cells[uu] <- ind.dead
        cell.den[ind.dead] <- Inf
      }
      
      # Update the cell coordinates and the density vector.
      x.coord <- x.coord[-dead.cells]
      y.coord <- y.coord[-dead.cells]
      cell.den <- cell.den[-dead.cells]
      
      if (length(cell.den) == 0) {
        # If all cells are dead, terminate the algorithm.
        return(list(diff.tot = NaN, err.mess = "All cells in the grid are dead."))
        break
      }
      
      # Cell mitosis: some of the current cells at the locations 
      # with the highest densities will undergo mitosis.
      
      # Choose the right parameter value for prob.prof based on the progress
      # of the invasion, since prob.prof has now become a time-varying
      # parameter. 
      prof.cells.num <- round(length(x.coord) * 
                                prob.prof.5.paras[ceiling(i/day.timesteps/3)])
      prof.cells <- rep(0, prof.cells.num)
      
      for(uu in 1:length(prof.cells)) {
        # Find the cells at the locations with the highest densities.
        ind.prof <- which(cell.den == max(cell.den))
        if (length(ind.prof) > 1) {
          ind.prof <- sample(ind.prof, 1, replace = FALSE)
        }
        prof.cells[uu] <- ind.prof
        cell.den[ind.prof] <- -Inf
      }
      
      # Record the current positions of the remaining cells.
      ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
      
      for (z in 1:length(x.coord)) {
        ind.position[x.coord[z], y.coord[z]] <- 1
      }
      
      # Proliferation mechanism
      for (q in 1:length(prof.cells)) {
        cell.position <- c(x.coord[prof.cells[q]], y.coord[prof.cells[q]])
        
        
        # Possible locations for daughter cells  (8 surrounding points.)
        right.position <- c(cell.position[1], (cell.position[2] + 1))
        right.down.position <- c((cell.position[1] + 1), (cell.position[2] + 1))
        down.position <- c((cell.position[1] + 1),  cell.position[2])
        up.right.position <- c((cell.position[1] - 1), (cell.position[2] + 1))
        up.position <- c((cell.position[1] - 1),  cell.position[2])
        up.left.position <- c((cell.position[1] - 1), (cell.position[2] - 1))  
        left.position <- c(cell.position[1], (cell.position[2] - 1))
        left.down.position <- c((cell.position[1] + 1), (cell.position[2] - 1))
        
        if (cell.position[1] == 1) {
          # Special case: top left corner
          if (cell.position[2] == 1) {
            # Possible directions to move, check if these points are occupied.
            right <- ind.position[right.position[1], right.position[2]]
            right.down <- ind.position[right.down.position[1], right.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            # Neighbouring status
            neighbouring.temp <- c(right, right.down, down)
            neighbouring.xcoord <- c(right.position[1], right.down.position[1], down.position[1]) # Neighbouring x coords
            neighbouring.ycoord <- c(right.position[2], right.down.position[2], down.position[2]) # Neighbouring y coords
            
            # If the cell has more than two neighbouring positions which are 
            # not occupied, it will proliferate. The original cell will vanish 
            # and split into two daughter cells, which will be randomly 
            # distributed into two unoccupied neighbouring locations.
            if(length(which(neighbouring.temp == 0)) >= 2) {
              # Indices of neighbouring status.
              empty.space <- which(neighbouring.temp == 0) 
              # Sample from the indices of unoccupied neighbouring coordinates.
              daughter.cells <- sample(empty.space, 2, replace = FALSE) 
              # Original cell vanishes.
              ind.position[cell.position[1], cell.position[2]] <- 0
              # Daughter cells being allocated.
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1 # Unoccupied neighbouring position now becomes occupied with the daughter cell.
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1 # Unoccupied neighbouring position now becomes occupied with the daughter cell. 
            }
            # Special case: top right corner (same reasoning is followed...) 
          } else if (cell.position[2] == space.length.x) {
            left <- ind.position[left.position[1], left.position[2]]
            left.down <- ind.position[left.down.position[1], left.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            neighbouring.temp <- c(left, left.down, down)
            neighbouring.xcoord <- c(left.position[1], left.down.position[1], down.position[1])
            neighbouring.ycoord <- c(left.position[2], left.down.position[2], down.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
          } else {
            # Other cells at the top boundary (same reasoning is followed...) 
            left <- ind.position[left.position[1], left.position[2]]
            right <- ind.position[right.position[1], right.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            left.down <- ind.position[left.down.position[1], left.down.position[2]]
            right.down <- ind.position[right.down.position[1], right.down.position[2]]
            
            neighbouring.temp <- c(left, right, down, left.down, right.down)
            neighbouring.xcoord <- c(left.position[1], right.position[1], down.position[1], left.down.position[1], right.down.position[1])
            neighbouring.ycoord <- c(left.position[2], right.position[2], down.position[2], left.down.position[2], right.down.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
          }
          # Lower boundary
        } else if (cell.position[1] == space.length.y) {
          # Special case: bottom left corner (same reasoning is followed...)
          if (cell.position[2] == 1) {
            up <- ind.position[up.position[1], up.position[2]]
            up.right <- ind.position[up.right.position[1], up.right.position[2]]
            right <- ind.position[right.position[1], right.position[2]]
            
            neighbouring.temp <- c(up, right, up.right)
            neighbouring.xcoord <- c(up.position[1], right.position[1], up.right.position[1])
            neighbouring.ycoord <- c(up.position[2], right.position[2], up.right.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
            # Special case: bottom right corner (same reasoning is followed...)
          } else if (cell.position[2] == space.length.x) {
            left <- ind.position[left.position[1], left.position[2]]
            up <- ind.position[up.position[1], up.position[2]]
            up.left <- ind.position[up.left.position[1], up.left.position[2]]
            
            neighbouring.temp <- c(left, up, up.left)
            neighbouring.xcoord <- c(left.position[1], up.position[1], up.left.position[1])
            neighbouring.ycoord <- c(left.position[2], up.position[2], up.left.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
            # Other cells on the lower boundary (same reasoning is followed...)
          } else {
            left <- ind.position[left.position[1], left.position[2]]
            right <- ind.position[right.position[1], right.position[2]]
            up <- ind.position[up.position[1], up.position[2]]
            up.left <- ind.position[up.left.position[1], up.left.position[2]]
            up.right <- ind.position[up.right.position[1], up.right.position[2]]
            
            neighbouring.temp <- c(left, up.left, up, up.right, right)
            neighbouring.xcoord <- c(left.position[1], up.left.position[1], 
                                     up.position[1], up.right.position[1], right.position[1])
            neighbouring.ycoord <- c(left.position[2], up.left.position[2],
                                     up.position[2], up.right.position[2], right.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
          }
          # Left boundary  
        } else if (cell.position[2] == 1) {
          # Special case: top left corner again 
          if (cell.position[1] == 1) {
            right <- ind.position[right.position[1], right.position[2]]
            right.down <- ind.position[right.down.position[1], right.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            # Do nothing actually, since the case of top left corner has been
            # evaluated already... 
            
            # Special case: bottom left corner again
          } else if (cell.position[1] == space.length.y) {
            up <- ind.position[up.position[1], up.position[2]]
            up.right <- ind.position[up.right.position[1], up.right.position[2]]
            right <- ind.position[right.position[1], right.position[2]]
            
            # Do nothing actually, since the case of bottom left corner has 
            # been evaluated already...
            
            # Other cells in the left boundary (same reasoning is followed...)
          } else {
            up <- ind.position[up.position[1], up.position[2]]
            up.right <- ind.position[up.right.position[1], up.right.position[2]]
            right <- ind.position[right.position[1], right.position[2]]
            right.down <- ind.position[right.down.position[1], right.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            neighbouring.temp <- c(up, up.right, right, right.down, down)
            neighbouring.xcoord <- c(up.position[1], up.right.position[1], right.position[1],
                                     right.down.position[1], down.position[1])
            neighbouring.ycoord <- c(up.position[2], up.right.position[2], right.position[2],
                                     right.down.position[2], down.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
          }
          # Right boundary
        } else if (cell.position[2] == space.length.x) {
          # Special case: top right corner again 
          if (cell.position[1] == 1) {
            left <- ind.position[left.position[1], left.position[2]]
            left.down <- ind.position[left.down.position[1], left.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            # Do nothing actually, since the case of top right corner has been
            # evaluated already.
            
            # Special case: bottom right corner again
          } else if (cell.position[1] == space.length.y) {
            left <- ind.position[left.position[1], left.position[2]]
            up.left <- ind.position[up.left.position[1], up.left.position[2]]
            up <- ind.position[up.position[1], up.position[2]]
            
            # Do nothing actually, since the case of bottom right corner has 
            # been evaluated already. 
            
            # Other cells at the right boundary (same reasoning is followed...)
          } else {
            up <- ind.position[up.position[1], up.position[2]]
            up.left <- ind.position[up.left.position[1], up.left.position[2]]
            left <- ind.position[left.position[1], left.position[2]]
            left.down <- ind.position[left.down.position[1], left.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            neighbouring.temp <- c(up, up.left, left, left.down, down)
            neighbouring.xcoord <- c(up.position[1], up.left.position[1], left.position[1],
                                     left.down.position[1], down.position[1])
            neighbouring.ycoord <- c(up.position[2], up.left.position[2], left.position[2],
                                     left.down.position[2], down.position[2])
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0)
              daughter.cells <- sample(empty.space, 2, replace = FALSE)
              ind.position[cell.position[1], cell.position[2]] <- 0
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
            }
          }
          # Other cells that are not at the boundary
          # (same reasoning is followed...)
        } else {
          left <- ind.position[left.position[1], left.position[2]]
          right <- ind.position[right.position[1], right.position[2]]
          up <- ind.position[up.position[1], up.position[2]]
          down <- ind.position[down.position[1], down.position[2]]
          up.left <- ind.position[up.left.position[1], up.left.position[2]]
          up.right <- ind.position[up.right.position[1], up.right.position[2]]
          left.down <- ind.position[left.down.position[1], left.down.position[2]]
          right.down <- ind.position[right.down.position[1], right.down.position[2]]
          
          neighbouring.temp <- c(left, right, up, down, up.left, up.right, left.down, right.down)
          neighbouring.xcoord <- c(left.position[1], right.position[1], up.position[1],
                                   down.position[1], up.left.position[1], up.right.position[1],
                                   left.down.position[1], right.down.position[1])
          neighbouring.ycoord <- c(left.position[2], right.position[2], up.position[2],
                                   down.position[2], up.left.position[2], up.right.position[2],
                                   left.down.position[2], right.down.position[2])
          
          if(length(which(neighbouring.temp == 0)) >= 2) {
            empty.space <- which(neighbouring.temp == 0)
            daughter.cells <- sample(empty.space, 2, replace = FALSE)
            ind.position[cell.position[1], cell.position[2]] <- 0
            ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1
            ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1
          }
        }
      }
      
      # Update the cell coordinates.
      cell.position.postprof <- which(ind.position == 1, arr.ind = TRUE)
      x.coord <- cell.position.postprof[,1]
      y.coord <- cell.position.postprof[,2]
      
      # Update the cell positions.
      ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
      
      for (z in 1:length(x.coord)) {
        ind.position[x.coord[z], y.coord[z]] <- 1
      }
    }
    
    
    # Solving the density PDE model
    if (i > (round(day.timesteps * (95/96)))) {
      # Diffusion starts having an impact after a certain amount of time in 
      # day 1.
      f[2:(length(y) - 1), 2:(length(x) - 1)] <- f[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - 
                                                                                              dt * eta.5.paras[ceiling(i/day.timesteps/3)] * m[2:(length(y) - 1), 2:(length(x) - 1)])
      
      n[2:(length(y) - 1), 2:(length(x) - 1)] <- n[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - 
                                                                                              (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - 
                                                                                              (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f[2:(length(y) - 1), 3:length(x)] + 
                                                                                                                                                             f[2:(length(y) - 1), 1:(length(x) - 2)] - 
                                                                                                                                                             4 * f[2:(length(y)-1), 2:(length(x)-1)] + 
                                                                                                                                                             f[1:(length(y) - 2),2:(length(x)-1)] + 
                                                                                                                                                             f[3:length(y),2:(length(x) - 1)])) + 
                                                                                              r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[2:(length(y) - 1), 2:(length(x) - 1)] - 
                                                                                                     f[2:(length(y) - 1), 2:(length(x) - 1)]) * dt) + 
        n[2:(length(y) - 1), 3:length(x)] * (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) - 
                                               (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f[2:(length(y) - 1), 3:length(x)] - 
                                                                                                                    f[2:(length(y) - 1), 1:(length(x) - 2)]))) + 
        n[2:(length(y) - 1), 1:(length(x) - 2)] * (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) +
                                                     (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f[2:(length(y) - 1), 3:length(x)] - 
                                                                                                                          f[2:(length(y) - 1), 1:(length(x) - 2)]))) + 
        n[1:(length(y) - 2), 2:(length(x) - 1)] * (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) - 
                                                     (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f[1:(length(y) - 2), 2:(length(x) - 1)] - 
                                                                                                                          f[3:length(y), 2:(length(x) - 1)]))) + 
        n[3:length(y), 2:(length(x) - 1)] * (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) + 
                                               (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f[1:(length(y) - 2), 2:(length(x) - 1)] - 
                                                                                                                    f[3:length(y), 2:(length(x)-1)])))
      
      m[2:(length(y) - 1), 2:(length(x) - 1)] <- m[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - (4 * dt * dm / (h ^ 2))) + 
        dt * alpha.5.paras[ceiling(i/day.timesteps/3)] * n[2:(length(y) - 1), 2:(length(x) - 1)] + 
        dt * dm / (h ^ 2) * (m[2:(length(y) - 1), 3:length(x)] + 
                              m[2:(length(y) - 1), 1:(length(x) - 2)] + 
                              m[1:(length(y) - 2), 2:(length(x) - 1)] + 
                              m[3:length(y), 2:(length(x)-1)])
      
    } else {
      # Otherwise the cell stay stationary. 
      f[2:(length(y) - 1), 2:(length(x) - 1)] <- f[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - 
                                                                                              dt * eta * m[2:(length(y) - 1), 2:(length(x) - 1)])
      
      n[2:(length(y) - 1), 2:(length(x) - 1)] <- n[2:(length(y) - 1), 2:(length(x) - 1)] * (1 + 
                                                                                              r * (1 - n[2:(length(y) - 1), 2:(length(x) - 1)] - 
                                                                                                     f[2:(length(y) - 1), 2:(length(x) - 1)]) * dt)
      
      m[2:(length(y) - 1), 2:(length(x) - 1)] <- m[2:(length(y) - 1), 2:(length(x) - 1)] + alpha * dt * n[2:(length(y) - 1), 2:(length(x) - 1)]
    }
    
    # Boundary condition
    n[1, ] <- n[2, ]
    n[, 1] <- n[, 2]
    n[length(y), ] <- n[(length(y) - 1), ]
    n[, length(x)] <- n[, (length(x) - 1)]
    
    f[1, ] <- f[2, ]
    f[, 1] <- f[, 2]
    f[length(y), ] <- f[(length(y) - 1), ]
    f[, length(x)] <- f[, (length(x) - 1)]
    
    m[1, ] <- m[2, ]
    m[, 1] <- m[, 2]
    m[length(y), ] <- m[(length(y) - 1), ]
    m[, length(x)] <- m[, (length(x) - 1)]
    
    # If singularity occurs while solving the PDE model, 
    # terminate the simulation. 
    den.data <- rbind(n,f,m)
    if (any(is.na(den.data)) || any(den.data < 0)) {
      return(list(diff.tot = NaN, err.mess = "Numerical singularity"))
      break
    }
    
    # Movement 
    
    if (i > (round(day.timesteps * (95/96)))) {
      # Movements only occur after diffusion starts having an impact on cells.
    for (j in 1:length(x.coord)) {
      if (y.coord[j] == 1) {
        # Special case: top left corner
        if (x.coord[j] == 1) {
          # Restricted from moving left and up
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1+f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Special case: bottom left corner
        } else if (x.coord[j] == space.length.y) {
          # Restricted from moving left and down
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- 0
          p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Other cells on the left boundary
        } else {
          # Restricted from moving left
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <-  (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <-  (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <-  (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        }
      } else if (y.coord[j] == space.length.x) {
        # Special case: top right corner
        if (x.coord[j] == 1) {
          # Restricted from moving right and up
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Special case: bottom right corner
        } else if (x.coord[j] == space.length.y) {
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- 0
          p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Other cells at the right boundary
        } else {
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        } 
      } else if (x.coord[j] == 1) {
        # Special case: top left corner again
        if (y.coord[j] == 1) {
          # Restricted from moving left and up
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Special case: top right corner again
        } else if (y.coord[j] == space.length.x) {
          # Restricted from moving right and up
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Other cells at the top boundary
        } else {
          # Restricted from moving up
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
        }
        # Bottom boundary
      } else if (x.coord[j] == space.length.y) {
        # Bottom left corner again
        if (y.coord[j] == 1) {
          # Restricted from moving left and down
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- 0
          p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Special case: bottom right corner again
        } else if (y.coord[j] == space.length.x) {
          # Restricted from moving right or down
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- 0
          p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Other cells at the bottom boundary
        } else {
          # Restricted from moving down
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- 0
          p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        }
        # Cells not at the boundary
      } else {
        # These cells can move in all directions
        f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
        f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
        f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
        f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
        
        p0 <- 1 - (4 * dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2) * (f.ip1j + 
                                                                       f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                       f.ijp1 + f.ijm1)) + 
          r.5.paras[ceiling(i/day.timesteps/3)] * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
        p1 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
        p2 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
        p3 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) - (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        p4 <- (dt * dn.5.paras[ceiling(i/day.timesteps/3)] / (h ^ 2)) + (dt * gamma.5.paras[ceiling(i/day.timesteps/3)] / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
      }
      
      # Roll the dice based on the calculated probabilities, see which 
      # direction the cell will choose to move to. 
      p <- c(p0, p1, p2, p3, p4)
      p[p < 0] <- 0
      p[p > 1] <- 1

      if (all(p == 0)) {
        mvment <- 0
      } else {
        mvment <- sample.int(5, 1, prob = p)
      }
      
      # Once the direction of movement is confirmed, if the designated 
      # position is not occupied, the cell will move to that location.
      if (mvment == 1) {
        x.coord[j] <- x.coord[j]
        y.coord[j] <- y.coord[j]
        ind.position[x.coord[j], y.coord[j]] <- 1
      } else if (mvment == 2) {
        if (ind.position[x.coord[j], (y.coord[j] - 1)] == 0) {
          ind.position[x.coord[j], y.coord[j]] <- 0
          x.coord[j] <- x.coord[j]
          y.coord[j] <- y.coord[j] - 1
          ind.position[x.coord[j], y.coord[j]] <- 1
        } else {
          x.coord[j] <- x.coord[j]
          y.coord[j] <- y.coord[j]
          ind.position[x.coord[j], y.coord[j]] <- 1
        }
      } else if (mvment == 3) {
        if (ind.position[x.coord[j], (y.coord[j] + 1)] == 0) {
          ind.position[x.coord[j], y.coord[j]] <- 0
          x.coord[j] <- x.coord[j]
          y.coord[j] <- y.coord[j] + 1
          ind.position[x.coord[j], y.coord[j]] <- 1
        } else {
          x.coord[j] <- x.coord[j]
          y.coord[j] <- y.coord[j]
          ind.position[x.coord[j], y.coord[j]] <- 1
        }
      } else if (mvment == 4) {
        if (ind.position[(x.coord[j] + 1), y.coord[j]] == 0) {
          ind.position[x.coord[j], y.coord[j]] <- 0
          x.coord[j] <- x.coord[j] + 1
          y.coord[j] <- y.coord[j]
          ind.position[x.coord[j], y.coord[j]] <- 1
        } else {
          x.coord[j] <- x.coord[j]
          y.coord[j] <- y.coord[j]
          ind.position[x.coord[j], y.coord[j]] <- 1
        }
      } else if (mvment == 5) {
        if (ind.position[(x.coord[j] - 1), y.coord[j]] == 0) {
          ind.position[x.coord[j], y.coord[j]] <- 0
          x.coord[j] <- x.coord[j] - 1
          y.coord[j] <- y.coord[j]
          ind.position[x.coord[j], y.coord[j]] <- 1
        } else {
          x.coord[j] <- x.coord[j]
          y.coord[j] <- y.coord[j]
          ind.position[x.coord[j], y.coord[j]] <- 1
        }
      }
      
    }
    }
 
    # Compute the density matrix at the end of day 3. 
    if (i == (day.timesteps * 3)) {
      density.mat.d3 <- matrix(0, nrow = length(y.cut), ncol = length(x.cut))
      ind.position.d3 <- ind.position
      for (ii in 1:length(y.cut)) {
        for (jj in 1:length(x.cut)) {
          ind.position.sub <- ind.position[y.cut[ii]:(y.cut[ii] + (mat.size - 1)), 
                                           x.cut[jj]:(x.cut[jj] + (mat.size - 1))]
          
          density.mat.d3[ii, jj] <- length(which(ind.position.sub == 1)) / (mat.size ^ 2)
        }
      }
      n.day3 <- n
      f.day3 <- f
      m.day3 <- m
    }
    
    # Compute the density matrix at the end of day 6. 
    if (i == (day.timesteps * 6)) {
      density.mat.d6 <- matrix(0, nrow = length(y.cut), ncol = length(x.cut))
      ind.position.d6 <- ind.position
      for (ii in 1:length(y.cut)) {
        for (jj in 1:length(x.cut)) {
          ind.position.sub <- ind.position[y.cut[ii]:(y.cut[ii] + (mat.size - 1)), 
                                           x.cut[jj]:(x.cut[jj] + (mat.size - 1))]
          
          density.mat.d6[ii, jj] <- length(which(ind.position.sub == 1)) / (mat.size ^ 2)
        }
      }
      n.day6 <- n
      f.day6 <- f
      m.day6 <- m
    }
    
    # Compute the density matrix at the end of day 9. 
    if (i == (day.timesteps * 9)) {
      density.mat.d9 <- matrix(0, nrow = length(y.cut), ncol = length(x.cut))
      ind.position.d9 <- ind.position
      for (ii in 1:length(y.cut)) {
        for (jj in 1:length(x.cut)) {
          ind.position.sub <- ind.position[y.cut[ii]:(y.cut[ii] + (mat.size - 1)), 
                                           x.cut[jj]:(x.cut[jj] + (mat.size - 1))]
          
          density.mat.d9[ii, jj] <- length(which(ind.position.sub == 1)) / (mat.size ^ 2)
        }
      }
      n.day9 <- n
      f.day9 <- f
      m.day9 <- m
    }
    
    # Compute the density matrix at the end of day 12. 
    if (i == (day.timesteps * 12)) {
      density.mat.d12 <- matrix(0, nrow = length(y.cut), ncol = length(x.cut))
      ind.position.d12 <- ind.position
      for (ii in 1:length(y.cut)) {
        for (jj in 1:length(x.cut)) {
          ind.position.sub <- ind.position[y.cut[ii]:(y.cut[ii] + (mat.size - 1)), 
                                           x.cut[jj]:(x.cut[jj] + (mat.size - 1))]
          
          density.mat.d12[ii, jj] <- length(which(ind.position.sub == 1)) / (mat.size ^ 2)
        }
      }
      n.day12 <- n
      f.day12 <- f
      m.day12 <- m
    }
    
    # Compute the density matrix at the end of day 14. 
    if (i == (day.timesteps * 14)) {
      density.mat.d14 <- matrix(0, nrow = length(y.cut), ncol = length(x.cut))
      ind.position.d14 <- ind.position
      for (ii in 1:length(y.cut)) {
        for (jj in 1:length(x.cut)) {
          ind.position.sub <- ind.position[y.cut[ii]:(y.cut[ii] + (mat.size - 1)), 
                                           x.cut[jj]:(x.cut[jj] + (mat.size - 1))]
          
          density.mat.d14[ii, jj] <- length(which(ind.position.sub == 1)) / (mat.size ^ 2)
        }
      }
      n.day14 <- n
      f.day14 <- f
      m.day14 <- m
    }
    
    #}
    # Optional, can be used to track the progress of the simulation
    # print(i)
  }
  
  # Sum of square measurements 
  sse.t3 <- sum((density.mat.d3 - t3.ref.den)^2)
  sse.t6 <- sum((density.mat.d6 - t6.ref.den)^2)
  sse.t9 <- sum((density.mat.d9 - t9.ref.den)^2)
  sse.t12 <- sum((density.mat.d12 - t12.ref.den)^2)
  sse.t14 <- sum((density.mat.d14 - t14.ref.den)^2)
  diff.tot <- sum(sse.t3, sse.t6, sse.t9, sse.t12, sse.t14)
  
  # Store the necessary simulation results as a list, which is also the output
  # of the function. 
  return(list(diff.tot = diff.tot, den.mat.d3 = density.mat.d3, 
              ind.position.d3 = ind.position.d3, sse.t3 = sse.t3,
              n.den.d3 = n.day3, f.den.d3 = f.day3, m.den.d3 = m.day3,
              
              den.mat.d6 = density.mat.d6, ind.position.d6 = ind.position.d6, 
              sse.t6 = sse.t6, n.den.d6 = n.day6, f.den.d6 = f.day6, 
              m.den.d6 = m.day6,
              
              den.mat.d9 = density.mat.d9, ind.position.d9 = ind.position.d9, 
              sse.t9 = sse.t9, n.den.d9 = n.day9, f.den.d9 = f.day9, 
              m.den.d9 = m.day9,
              
              den.mat.d12 = density.mat.d12, 
              ind.position.d12 = ind.position.d12, sse.t12 = sse.t12,
              n.den.d12 = n.day12, f.den.d12 = f.day12, m.den.d12 = m.day12,
              
              den.mat.d14 = density.mat.d14, 
              ind.position.d14 = ind.position.d14, sse.t14 = sse.t14,
              n.den.d14 = n.day14, f.den.d14 = f.day14, m.den.d14 = m.day14
  ))
}

# Calculate the correct bw that would yield the desirable ESS.
calculate.bw <- function(ss.mat, lb.bw, ub.bw, ess.target, step.size) {
  # Purpose: calculate the correct bandwidth factor which can yield the
  # desirable effective sample size (ESS).
  
  # Arguments: 
  # ss.mat: matrix of least square differences. 
  # lb.bw: lower bound of the decimal number sequence to be searched for
  #        the correct bandwidth factor. 
  # ub.bw: upper bound of the decimal number sequence to be searched for 
  #        the correct bandwidth factor. 
  # ess.target: desirable ESS. 
  # step.size: step of increment for the search of correct bandwidth factor
  #            in the sequence. 
  
  # Take away the singular outputs in the least square matrix. 
  ss.mat <- as.matrix(ss.mat)
  ind.nan <- which(is.na(ss.mat[,2]))
  if (length(ind.nan) == 0) {
    ss.mat.valid <- ss.mat
  } else {
    ss.mat.valid <- ss.mat[-ind.nan, ]
  }
  
  # Construct the sequence of bandwidth factor to be searched. 
  power <- seq(lb.bw, ub.bw, by = step.size)
  
  # An empty vector to store the ESS for each bandwidth factor. 
  ess.vec <- rep(0, length = length(power))
  
  # For every bandwidth factor in the sequence:
  for (i in 1:length(power)) {
    # Calculate the weights. 
    wt.temp <- ss.mat.valid[,2]^(-power[i])
    # Rescale the weights to resampling probabilities. 
    resamp.prob.temp <- rep(0, length = length(ss.mat.valid[,1]))
    
    for (j in 1:length(wt.temp)) {
      if (wt.temp[j] == min(wt.temp)) {
        resamp.prob.temp[j] <- 0
      } else if (wt.temp[j] == max(wt.temp)) {
        resamp.prob.temp[j] <- 1
      } else {
        resamp.prob.temp[j] <- (wt.temp[j]-min(wt.temp))/(max(wt.temp)-min(wt.temp))
      }
    }
    # ESS calculation based on the rescaled weights obtained using the current
    # bandwidth factor. 
    ess.temp <- ((sum(resamp.prob.temp))^2)/sum(resamp.prob.temp^2)
    ess.vec[i] <- ess.temp
    # Optional, can be used to track the progress of the search. 
    print(i)
  }
  
  # Matrix of bandwidth factors and their corresponding ESSs'.
  ess.mat <- cbind(power, ess.vec)
  # Take away the singular result. 
  ind.nan.ess <- which(is.na(ess.mat[,2]))
  
  if (length(ind.nan.ess) == 0) {
    ess.mat.valid <- ess.mat
  } else {
    ess.mat.valid <- ess.mat[-ind.nan.ess,]
  }
  
  # Locate the bandwidth factor which gives an ESS closest to the 
  # desirable one. 
  ess.obj.ind <- which(abs(ess.mat.valid[,2] - ess.target) == 
                         min(abs(ess.mat.valid[,2] - ess.target)))
  ess.obj <- ess.mat.valid[ess.obj.ind,2]
  bw.obj <- ess.mat.valid[ess.obj.ind,1]
  
  # Calculate the corresponding weights and resampling probabilities using 
  # the bandwidth factor located previously. 
  wt.obj <- ss.mat.valid[,2]^(-bw.obj)
  resamp.prob.obj <- rep(0, length = length(ss.mat.valid[,1]))
  
  for (k in 1:length(wt.obj)) {
    if (wt.obj[k] == min(wt.obj)) {
      resamp.prob.obj[k] <- 0
    } else if (wt.obj[k] == max(wt.obj)) {
      resamp.prob.obj[k] <- 1
    } else {
      resamp.prob.obj[k] <- (wt.obj[k]-min(wt.obj))/(max(wt.obj)-min(wt.obj))
    }
  }
  
  # Output. 
  info.list <- list(info.mat = cbind(ss.mat.valid, wt.obj, resamp.prob.obj),
                    ess.mat = ess.mat,
                    ess.obj = ess.obj,
                    bw.obj = bw.obj)
  
  return(info.list)
}

# The ABC scheme
abc_bcd <- function(info.mat, paras) {
  # Purpose: Given the current set of parameters and the corresponding results 
  # of summary statistics, generate a new set of parameters to be evaluated in 
  # the next round. 
  
  # Arguments: 
  # info.mat: information regarding the current set of least square differences.
  # paras: parameters to be resampled and perturbed. 
  
  # Lower and upper bounds for the parameters. 
  paras_lb <- c(0.000069, 0.005, 0.0008, 7, 0.0001, 0.07, 1, 0.01, 0.2, 
                dn.quad.coef/2, dn.lin.coef*2, gamma.quad.coef*2, gamma.lin.coef/2,
                rn.quad.coef*2, rn.lin.coef/2, eta.quad.coef/2, eta.lin.coef*2,
                alpha.quad.coef/2, alpha.lin.coef*2, 
                prob.prof.quad.coef*2, prob.prof.lin.coef/2)
  paras_ub <- c(0.02, 0.26, 0.08, 18, 0.033, 0.18, 5, 0.1, 1, 
                dn.quad.coef*2, dn.lin.coef/2, gamma.quad.coef/2, gamma.lin.coef*2,
                rn.quad.coef/2, rn.lin.coef*2, eta.quad.coef*2, eta.lin.coef/2,
                alpha.quad.coef*2, alpha.lin.coef/2, 
                prob.prof.quad.coef/2, prob.prof.lin.coef*2)
  
  # Resample the indices based on the resampling probabilities. 
  resamp_ind <- sample(info.mat[,1],size = length(paras[,1]), 
                       replace = TRUE, prob = info.mat[,length(info.mat[1,])]) 
  
  # Resampled parameter vectors, without perturbation.
  paras_nr_unperturbed <- paras[resamp_ind,]
  
  # An empty matrix used to store the perturbed parameter values.
  paras_nr_perturbed <- matrix(0,nrow = nrow(paras),ncol = ncol(paras))
  
  # Perturbation step (parameter values first) 
  for (i in 1:length(paras[1,])) {
    for (j in 1:length(paras[,1])) {
      h = sqrt(1-0.05^2)
      paras_nr_perturbed[j,i] <- rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                       0.05*sd(paras_nr_unperturbed[,i]))
      while ((paras_nr_perturbed[j,i] > paras_ub[i]) || (paras_lb[i] > paras_nr_perturbed[j,i])) {
        paras_nr_perturbed[j,i] <- rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                         0.05*sd(paras_nr_unperturbed[,i]))
      }
      #print(j)
    }
    print(i)
  }  
  
  # Perturbation step (slopes)
  # for (i in 1:length(paras[,1])) {
  #  for (j in 10:11) {
  #    h = sqrt(1-0.05^2)
  #    paras_nr_perturbed[i,j] <- rnorm(1, h*paras_nr_unperturbed[i,j]+(1-h)*mean(paras_nr_unperturbed[,j]),
  #                                     0.05*sd(paras_nr_unperturbed[,j]))
  #  }
  #  dn.init.par.temp <- paras_nr_perturbed[i,1]
  #  while((paras_nr_perturbed[i,j] > 
  #         paras_ub[j]) || (paras_lb[j] > 
  #                          paras_nr_perturbed[i,j]) || (length(dn.varying.paras(dn.init.par = 
  #                                                                        dn.init.par.temp, 
  #                                                                        dn.quad.temp = 
  #                                                                        paras_nr_perturbed[i, 10], 
  #                                                                        dn.lin.temp = 
  #                                                                        paras_nr_perturbed[i, 11])) < 7)) {
  #    paras_nr_perturbed[i,j] <- rnorm(1, h*paras_nr_unperturbed[i,j]+(1-h)*mean(paras_nr_unperturbed[,j]),
  #                                     0.05*sd(paras_nr_unperturbed[,j]))
  #  }
  #  print(i)
  #}
  
  # Parameter values for the next round. 
  return(paras_nr_perturbed)
}

reduction.disp <- function(prev, curr) {
  # Purpose: check if the decrease in averaged least square differences has 
  # dropped below 5%. 
  
  # Arguments:
  # prev: averaged least square difference in the previous round. 
  # curr: averaged least square difference in the current round. 
  
  # Calculate and return the percentage of decrease. 
  per <- (prev-curr)/prev*100
  return(per)
}