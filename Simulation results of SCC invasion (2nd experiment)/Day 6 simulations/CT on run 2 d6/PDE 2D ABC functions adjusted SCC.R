# PDE 2D ABC functions adjusted
# Author: Yunchen Xiao
# This .R file contains all the necessary functions for the simulation study
# of the post-day 3 SCC invasion patterns. 

# Read in the cell densities at day 6. 
library(readr)
t6.dat.p1 <- read.csv("Summary pic day 6 p1.csv")
t6.dat.p2 <- read.csv("Summary pic day 6 p2.csv")
t6.dat.p3 <- read.csv("Summary pic day 6 p3.csv")
t6.dat.p4 <- read.csv("Summary pic day 6 p4.csv")
t6.dat.p5 <- read.csv("Summary pic day 6 p5.csv")
t6.dat.p6 <- read.csv("Summary pic day 6 p6.csv")
t6.dat.p7 <- read.csv("Summary pic day 6 p7.csv")

t6.ref.den <- cbind(rev(t6.dat.p1[, 5]), rev(t6.dat.p2[, 5]), 
                    rev(t6.dat.p3[, 5]), rev(t6.dat.p4[, 5]), 
                    rev(t6.dat.p5[, 5]), rev(t6.dat.p6[, 5]),
                    rev(t6.dat.p7[, 5]))/100

# Initial condition (final simulation output of post-day 3 pattern)
res.d3 <- read_rds("Initial condition day 6 simulation full.rds")
n.den.d3 <- res.d3$n.den.d3
f.den.d3 <- res.d3$f.den.d3
m.den.d3 <- res.d3$m.den.d3
ind.position.d3 <- res.d3$ind.position.d3

################################################################################

generate.pattern <- function(par, prob.death, prob.prof) {
  # Purpose: Generate a SCC invasion pattern with the given
  # parameters.
  
  # Arguments:
  #    par: parameters of the PDE model
  #    prob.death: proportion of cells that will be undergo extinction at the 
  #    end of every day.
  #    prob.prof: proportion of cells that will undergo mitosis at the end 
  #    of every day.
  
  
  # Space discretization
  # Create a 60*35 domain. 
  h <- 1/59
  space.length.y <- (1/h) + 1
  space.length.x <- round(space.length.y * (280/480))
  x <- seq(0, h * (space.length.x - 1), length.out = space.length.x)
  y <- seq(0, 1, length.out = space.length.y)
  
  # Time discretization
  T <- 4.52
  # dt chosen as 0.0025 so the stability condition (dn < ((h^2)/4dt)) 
  # can be maintained.
  dt <- 0.0025
  # Total dimensionless timesteps and the dimensionless timesteps for one 
  # single day. 
  timesteps <- T/dt
  int.timesteps <- 1/dt
  day.timesteps <- 600
  
  # Parameters of PDE
  dn <- par[1]
  gamma <- par[2]
  r <- par[3]
  eta <- par[4]
  dm <- par[5]
  alpha <- par[6]
  beta <- 0
  
  # Initial condition
  n0 <- n.den.d3
  
  f0 <- f.den.d3
  
  m0 <- m.den.d3
  
  n <- n0
  f <- f0
  m <- m0
  
  # Sort the initial cells
  ind.cell <- which(ind.position.d3 == 1, arr.ind = TRUE)
  x.coord <- ind.cell[, 1]
  y.coord <- ind.cell[, 2]
  ind.position <- ind.position.d3
  
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
        # If all cells are dead, terminate the algorithm
        return(NaN)
        break
      }
      
      # Cell mitosis: some of the current cells at the locations 
      # with the highest densities will undergo mitosis.
      
      prof.cells.num <- round(length(x.coord) * prob.prof)
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
        
        
        # Possible locations for daughter cells (8 surrounding points.)
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
            neighbouring.xcoord <- c(right.position[1], right.down.position[1], down.position[1]) 
            neighbouring.ycoord <- c(right.position[2], right.down.position[2], down.position[2]) 
            
            # If the cell has more than two neighbouring positions which are 
            # not occupied, it will proliferate. The original cell will vanish 
            # and split into two daughter cells, which will be randomly 
            # distributed into two unoccupied neighbouring locations.
            if(length(which(neighbouring.temp == 0)) >= 2) {
              # Indices of neighbouring status
              empty.space <- which(neighbouring.temp == 0) 
              # Sample from the indices of unoccupied neighbouring coordinates.
              daughter.cells <- sample(empty.space, 2, replace = FALSE) 
              # Original cell vanishes.
              ind.position[cell.position[1], cell.position[2]] <- 0 
              # Daughter cells being allocated.
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1 
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
    
    # Solving the PDE model numerically
    # if (i > (round(day.timesteps * (95/96)))) {
      
      f[2:(length(y) - 1), 2:(length(x) - 1)] <- f[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - 
                                                                                              dt * eta * m[2:(length(y) - 1), 2:(length(x) - 1)])
      
      n[2:(length(y) - 1), 2:(length(x) - 1)] <- n[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - 
                                                                                              (4 * dt * dn / (h ^ 2)) - 
                                                                                              (dt * gamma / (h ^ 2) * (f[2:(length(y) - 1), 3:length(x)] + 
                                                                                                                         f[2:(length(y) - 1), 1:(length(x) - 2)] - 
                                                                                                                         4 * f[2:(length(y)-1), 2:(length(x)-1)] + 
                                                                                                                         f[1:(length(y) - 2),2:(length(x)-1)] + 
                                                                                                                         f[3:length(y),2:(length(x) - 1)])) + 
                                                                                              r * (1 - n[2:(length(y) - 1), 2:(length(x) - 1)] - 
                                                                                                     f[2:(length(y) - 1), 2:(length(x) - 1)]) * dt) + 
        n[2:(length(y) - 1), 3:length(x)] * (dt * dn / (h ^ 2) - 
                                               (dt * gamma / (4 * (h ^ 2)) * (f[2:(length(y) - 1), 3:length(x)] - 
                                                                                f[2:(length(y) - 1), 1:(length(x) - 2)]))) + 
        n[2:(length(y) - 1), 1:(length(x) - 2)] * (dt * dn / (h ^ 2) +
                                                     (dt * gamma / (4 * (h ^ 2)) * (f[2:(length(y) - 1), 3:length(x)] - 
                                                                                      f[2:(length(y) - 1), 1:(length(x) - 2)]))) + 
        n[1:(length(y) - 2), 2:(length(x) - 1)] * (dt * dn / (h ^ 2) - 
                                                     (dt * gamma / (4 * (h ^ 2)) * (f[1:(length(y) - 2), 2:(length(x) - 1)] - 
                                                                                      f[3:length(y), 2:(length(x) - 1)]))) + 
        n[3:length(y), 2:(length(x) - 1)] * (dt * dn / (h ^ 2) + 
                                               (dt * gamma / (4 * (h ^ 2)) * (f[1:(length(y) - 2), 2:(length(x) - 1)] - 
                                                                                f[3:length(y), 2:(length(x)-1)])))
      
      m[2:(length(y) - 1), 2:(length(x) - 1)] <- m[2:(length(y) - 1), 2:(length(x) - 1)] * (1 - (4 * dt * dm / (h ^ 2))) + 
        dt * alpha * n[2:(length(y) - 1), 2:(length(x) - 1)] + 
        dt * dm / (h ^ 2) * (m[2:(length(y) - 1), 3:length(x)] + 
                               m[2:(length(y) - 1), 1:(length(x) - 2)] + 
                               m[1:(length(y) - 2), 2:(length(x) - 1)] + 
                               m[3:length(y),2:(length(x)-1)])

    # }
    
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
      return(NaN)
      break
    }
    
    # Movement 
    #if (i > (round(day.timesteps * (95/96)))) {
    for (j in 1:length(x.coord)) {
      if (y.coord[j] == 1) {
        # Special case: top left corner
        if (x.coord[j] == 1) {
          # Restricted from moving left and up
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1+f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Special case: bottom left corner
        } else if (x.coord[j] == space.length.y) {
          # Restricted from moving left and down
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- 0
          p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Other cells on the left boundary.
        } else {
          # Restricted from moving left
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <-  (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <-  (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <-  (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        }
      } else if (y.coord[j] == space.length.x) {
        # Special case: top right corner
        if (x.coord[j] == 1) {
          # Restricted from moving right and up
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Special case: bottom right corner
        } else if (x.coord[j] == space.length.y) {
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- 0
          p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Other cells at the right boundary
        } else {
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        } 
      } else if (x.coord[j] == 1) {
        # Special case: top left corner again
        if (y.coord[j] == 1) {
          # Restricted from moving left and up
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1+f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Special case: top right corner again
        } else if (y.coord[j] == space.length.x) {
          # Restricted from moving right and up
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
          # Other cells at the top boundary
        } else {
          # Restricted from moving up
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- 0
          f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          p4 <- 0
        }
        # Bottom boundary
      } else if (x.coord[j] == space.length.y) {
        # Special case: bottom left corner again
        if (y.coord[j] == 1) {
          # Restricted from moving left and down
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- 0
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j], y.coord[j]]) * dt
          p1 <- 0
          p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- 0
          p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Special case: bottom right corner again
        } else if (y.coord[j] == space.length.x) {
          # Restricted from moving right or down
          f.ip1j <- 0
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- 0
          p3 <- 0
          p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
          # Other cells at the bottom boundary
        } else {
          # Restricted from moving down
          f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
          f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
          f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
          f.ijm1 <- 0
          
          p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                         f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                         f.ijp1 + f.ijm1)) + 
            r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
          p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
          p3 <- 0
          p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        }
        # Cells not at the boundary
      } else {
        # These cells can move in all directions
        f.ip1j <- f[x.coord[j], (y.coord[j] + 1)]
        f.im1j <- f[x.coord[j], (y.coord[j] - 1)]
        f.ijp1 <- f[(x.coord[j] - 1), y.coord[j]]
        f.ijm1 <- f[(x.coord[j] + 1), y.coord[j]]
        
        p0 <- 1 - (4 * dt * dn / (h ^ 2)) - (dt * gamma / (h ^ 2) * (f.ip1j + 
                                                                       f.im1j - 4 * f[x.coord[j], y.coord[j]] + 
                                                                       f.ijp1 + f.ijm1)) + 
          r * (1 - n[x.coord[j], y.coord[j]] - f[x.coord[j],y.coord[j]]) * dt
        p1 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
        p2 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ip1j - f.im1j))
        p3 <- (dt * dn / (h ^ 2)) - (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
        p4 <- (dt * dn / (h ^ 2)) + (dt * gamma / (4 * (h ^ 2)) * (f.ijp1 - f.ijm1))
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
    # }
    
    # Compute the density matrix at the end of day 6. 
    if (i == (day.timesteps * 3)) {
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
    #}
    # Optional, can be used to track the progress of the simulation
    # print(i)
  }
  
  return(list(den.mat.d6 = density.mat.d6, ind.position.d6 = ind.position.d6,
              n.den.d6 = n.day6,
              f.den.d6 = f.day6, 
              m.den.d6 = m.day6
  ))
}

calculate.sse <- function(pars, prob.death, prob.prof) {
  # Purpose: calculates the least square difference between the simulated
  # invasion pattern and the reference invasion pattern.
  # Arguments: 
  #    par: parameters of the PDE model
  #    prob.death: proportion of cells that will be undergo extinction at the 
  #    end of every day.
  #    prob.prof: proportion of cells that will undergo mitosis at the end 
  #    of every day.
  
  # Call the "generate.pattern" function.
  temp.den.table <- generate.pattern(par = pars,
                                     prob.death = prob.death, prob.prof = prob.prof)
  
  if (is.na(temp.den.table)) {
    return(list(diff = NaN))
    # If the result from "generate.pattern" is NaN, return NaN.
  } else {
    # Calculate the least square differences.
    sse.t6 <- sum((temp.den.table$den.mat.d6 - t6.ref.den)^2)
    # Output. 
    return(list(diff = sse.t6,
                den.mat.d6 = temp.den.table$den.mat.d6, 
                ind.position.d6 = temp.den.table$ind.position.d6,
                n.den.d6 = temp.den.table$n.den.d6,
                f.den.d6 = temp.den.table$f.den.d6, 
                m.den.d6 = temp.den.table$m.den.d6
                ))
  }
}

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
abc_bcd <- function(info.mat,paras) {
  # Purpose: Given the current set of parameters and the corresponding results 
  # of summary statistics, generate a new set of parameters to be evaluated in 
  # the next round. 
  
  # Arguments: 
  # info.mat: information regarding the current set of least square differences.
  # paras: parameters to be resampled and perturbed. 
  
  # Lower and upper bounds for the parameters. 
  paras_lb <- c(0.000069, 0.005, 0.0008, 7, 0.0001, 0.07, 0.01, 0.2)
  paras_ub <- c(0.02, 0.26, 0.08, 18, 0.033, 0.18, 0.1, 1)
  
  # Resample the indices based on the resampling probabilities. 
  resamp_ind <- sample(info.mat[,1],size = length(paras[,1]), 
                       replace = TRUE, prob = info.mat[,length(info.mat[1,])])
  
  # Resampled parameter vectors, without perturbation.
  paras_nr_unperturbed <- paras[resamp_ind,] 
  
  # An empty matrix used to store the perturbed parameter values.
  paras_nr_perturbed <- matrix(0,nrow = nrow(paras),ncol = ncol(paras)) 
  
  # Perturbation step. 
  for (i in 1:length(paras[1,])) {
    for (j in 1:length(paras[,1])){
      h = sqrt(1-0.05^2)
      paras_nr_perturbed[j,i] <- rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                       0.05*sd(paras_nr_unperturbed[,i]))
      while ((paras_nr_perturbed[j,i] > paras_ub[i]) || (paras_lb[i] > paras_nr_perturbed[j,i])) {
        paras_nr_perturbed[j,i] <- rnorm(1,h*paras_nr_unperturbed[j,i]+(1-h)*mean(paras_nr_unperturbed[,i]),
                                         0.05*sd(paras_nr_unperturbed[,i]))
      }
      print(j)
    }
  }  
  
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