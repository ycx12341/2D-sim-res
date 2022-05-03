# PDE 2D ABC functions adjusted

# Author: Yunchen Xiao

# Read in the cell densities at different sections of the domain
# at the specified time points. 

library(readr)
# t12 reference density 
t12.dat.p1 <- read.csv("Summary pic day 12 p1.csv")
t12.dat.p2 <- read.csv("Summary pic day 12 p2.csv")
t12.dat.p3 <- read.csv("Summary pic day 12 p3.csv")
t12.dat.p4 <- read.csv("Summary pic day 12 p4.csv")
t12.dat.p5 <- read.csv("Summary pic day 12 p5.csv")
t12.dat.p6 <- read.csv("Summary pic day 12 p6.csv")
t12.dat.p7 <- read.csv("Summary pic day 12 p7.csv")

t12.ref.den <- cbind(rev(t12.dat.p1[, 5]), rev(t12.dat.p2[, 5]), rev(t12.dat.p3[, 5]),
                    rev(t12.dat.p4[, 5]), rev(t12.dat.p5[, 5]), rev(t12.dat.p6[, 5]),
                    rev(t12.dat.p7[, 5]))/100

# Initial condition
res.d9 <- read_rds("Initial condition day 12 simulation full.rds")
n.den.d9 <- res.d9$n.den.d9
f.den.d9 <- res.d9$f.den.d9
m.den.d9 <- res.d9$m.den.d9
ind.position.d9 <- res.d9$ind.position.d9

################################################################################

generate.pattern <- function(par, prob.death, prob.prof) {
  # Purpose: Generate a 2D cancer cells invasion pattern with the given
  # parameters.
  
  # Arguments:
  #    par: parameters of the PDE model
  #    init.cells.cols: number of cell columns being set at the left of the domain
  #                     initially (we assume cells invade from left to right)
  #    prob.death: proportion of cells that will be dead at the end of every day.
  #    prob.prof: proportion of cells that will proliferate at the end of every day.
  
  # set.seed(123)
  # RNGkind(sample.kind = "Rejection")
  # Space discretization
  # h <- 1/47
  # space.length.y <- (1/h) + 1
  # space.length.x <- round(space.length.y * (280/480))
  # x <- seq(0, h * (space.length.x - 1), length.out = space.length.x)
  # y <- seq(0, 1, length.out = space.length.y)

  # Time discretization
  # T <- 4.5
  # dt <- 0.004
  # timesteps <- T/dt
  # int.timesteps <- 1/dt
  # day.timesteps <- 375
  
  # Space discretization
  h <- 1/59
  space.length.y <- (1/h) + 1
  space.length.x <- round(space.length.y * (280/480))
  x <- seq(0, h * (space.length.x - 1), length.out = space.length.x)
  y <- seq(0, 1, length.out = space.length.y)
  
  # Time discretization
  T <- 4.5
  dt <- 0.0025
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
  n0 <- n.den.d9

  f0 <- f.den.d9
  
  m0 <- m.den.d9
  
  n <- n0
  f <- f0
  m <- m0
  
  # Sort the initial cells
  ind.cell <- which(ind.position.d9 == 1, arr.ind = TRUE)
  x.coord <- ind.cell[, 1]
  y.coord <- ind.cell[, 2]
  ind.position <- ind.position.d9
  
  mat.size <- space.length.y/12
  y.cut <- seq(1, space.length.y, by = mat.size)
  x.cut <- seq(1, space.length.x, by = mat.size)
  
  # Numerical scheme which solves the PDE system
  for (i in 1:timesteps) {
    
    # At the end of everyday, some of the current cells
    # in the domain will die or proliferate
    
    if (i %% day.timesteps == 0) {
      cell.den <- rep(0, length(x.coord))
      
      for (u in 1:length(cell.den)) {
        cell.den[u] <- n[x.coord[u], y.coord[u]]      
      }
      
      # Cell death: some of the current cells at the locations with the 
      # lowest cell densities will extinct.
      
      dead.cells.num <- round(length(x.coord) * prob.death)
      dead.cells <- rep(0, dead.cells.num)
      
      for (uu in 1:length(dead.cells)) {
        ind.dead <- which(cell.den == min(cell.den))
        if(length(ind.dead) > 1) {
          ind.dead <- sample(ind.dead, 1, replace = FALSE)
        }
        dead.cells[uu] <- ind.dead
        cell.den[ind.dead] <- Inf
      }
      
      x.coord <- x.coord[-dead.cells]
      y.coord <- y.coord[-dead.cells]
      cell.den <- cell.den[-dead.cells]
      
      if (length(cell.den) == 0) {
        # If all cells are dead, end the algorithm
        return(NaN)
        break
      }
      
      # Cell proliferation: some of the current cells at the locations 
      # with the highest densities will proliferate.
      
      prof.cells.num <- round(length(x.coord) * prob.prof)
      prof.cells <- rep(0, prof.cells.num)
      
      for(uu in 1:length(prof.cells)) {
        ind.prof <- which(cell.den == max(cell.den))
        if (length(ind.prof) > 1) {
          ind.prof <- sample(ind.prof, 1, replace = FALSE)
        }
        prof.cells[uu] <- ind.prof
        cell.den[ind.prof] <- -Inf
      }
      
      # Current positions of the remaining cells
      ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
      
      for (z in 1:length(x.coord)) {
        ind.position[x.coord[z], y.coord[z]] <- 1
      }
      
      # Proliferation mechanism
      for (q in 1:length(prof.cells)) {
        cell.position <- c(x.coord[prof.cells[q]], y.coord[prof.cells[q]])
        
        
        # Possible locations for daughter cells
        right.position <- c(cell.position[1], (cell.position[2] + 1))
        right.down.position <- c((cell.position[1] + 1), (cell.position[2] + 1))
        down.position <- c((cell.position[1] + 1),  cell.position[2])
        up.right.position <- c((cell.position[1] - 1), (cell.position[2] + 1))
        up.position <- c((cell.position[1] - 1),  cell.position[2])
        up.left.position <- c((cell.position[1] - 1), (cell.position[2] - 1))  
        left.position <- c(cell.position[1], (cell.position[2] - 1))
        left.down.position <- c((cell.position[1] + 1), (cell.position[2] - 1))
        
        if (cell.position[1] == 1) {
          # Top left corner
          if (cell.position[2] == 1) {
            # Possible directions to move, check if they are occupied or not
            right <- ind.position[right.position[1], right.position[2]]
            right.down <- ind.position[right.down.position[1], right.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            neighbouring.temp <- c(right, right.down, down) # Neighbouring status
            neighbouring.xcoord <- c(right.position[1], right.down.position[1], down.position[1]) # Neighbouring x coords
            neighbouring.ycoord <- c(right.position[2], right.down.position[2], down.position[2]) # Neighbouring y coords
            
            # If the cell has more than two neighbouring positions which are not occupied, 
            # it will proliferate. The original cell will vanish and split into two daughter
            # cells, which are randomly distributed into two unoccupied neighbouring locations.
            
            if(length(which(neighbouring.temp == 0)) >= 2) {
              empty.space <- which(neighbouring.temp == 0) # Indices of neighbouring status
              daughter.cells <- sample(empty.space, 2, replace = FALSE) # Sample from the indices of unoccupied neighbouring coordinates.
              ind.position[cell.position[1], cell.position[2]] <- 0 # Original cell vanish
              ind.position[neighbouring.xcoord[daughter.cells[1]], neighbouring.ycoord[daughter.cells[1]]] <- 1 # Unoccupied neighbouring position now becomes occupied with the daughter cell.
              ind.position[neighbouring.xcoord[daughter.cells[2]], neighbouring.ycoord[daughter.cells[2]]] <- 1 # Unoccupied neighbouring position now becomes occupied with the daughter cell. 
            }
            # Top right corner  
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
            # Other cells at the top boundary
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
          # Bottom left corner
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
            # Bottom right corner
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
            # Other cells on the lower boundary  
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
          # Top left corner again 
          if (cell.position[1] == 1) {
            right <- ind.position[right.position[1], right.position[2]]
            right.down <- ind.position[right.down.position[1], right.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            # Do nothing actually... 
            
            # Bottom left corner again
          } else if (cell.position[1] == space.length.y) {
            up <- ind.position[up.position[1], up.position[2]]
            up.right <- ind.position[up.right.position[1], up.right.position[2]]
            right <- ind.position[right.position[1], right.position[2]]
            
            # Do nothing actually... 
            
            # Other cells in the left boundary
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
          # Top right corner again
          if (cell.position[1] == 1) {
            left <- ind.position[left.position[1], left.position[2]]
            left.down <- ind.position[left.down.position[1], left.down.position[2]]
            down <- ind.position[down.position[1], down.position[2]]
            
            # Do nothing actually... 
            
            # Bottom right corner again
          } else if (cell.position[1] == space.length.y) {
            left <- ind.position[left.position[1], left.position[2]]
            up.left <- ind.position[up.left.position[1], up.left.position[2]]
            up <- ind.position[up.position[1], up.position[2]]
            
            # Do nothing actually... 
            
            # Other cells at the right boundary
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
      
      cell.position.postprof <- which(ind.position == 1, arr.ind = TRUE)
      x.coord <- cell.position.postprof[,1]
      y.coord <- cell.position.postprof[,2]
      
      ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
      
      for (z in 1:length(x.coord)) {
        ind.position[x.coord[z], y.coord[z]] <- 1
      }
    }
    
    # Note that it is alright to count the cells at the corners twice, since they have
    # at most three possible locations to distribute its daughter cells if they 
    # were chosen to proliferate, if they have already proliferated once, then it 
    # is impossible for them to proliferate when they are investigated again. 
    
    # Current positions of the remaining cells
    #ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
    
    #for (z in 1:length(x.coord)) {
    #  ind.position[x.coord[z], y.coord[z]] <- 1
    #}
    
    # Solving the density PDE model
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
    
    #n.nan <- length(which(n == "NaN"))
    #f.nan <- length(which(f == "NaN"))
    #m.nan <- length(which(m == "NaN"))
    
    #if (sum(n.nan, f.nan, f.nan) > 0) {
    #  return(NaN)
    # If any singularity occurs in the algorithm, return
    # NaN and stop the algorithm.
    #  break
    #}
    
    den.data <- rbind(n,f,m)
    if (any(is.na(den.data)) || any(den.data < 0)) {
      return(NaN)
      break
    }
    # Movement 
    
    #if (i > (round(day.timesteps * (95/96)))) {
    for (j in 1:length(x.coord)) {
      if (y.coord[j] == 1) {
        # Top left corner
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
          # Bottom left corner
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
          # Other cells on the left boundary
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
        # Top right corner
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
          # Bottom right corner
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
        # Top left corner again
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
          # Top right corner again
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
        # Bottom left corner again
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
          # Bottom right corner again
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
      
      # Roll the dice based on the calculated probabilities, see which direction the
      # cell will choose to move to. 
      p <- c(p0, p1, p2, p3, p4)
      p[p < 0] <- 0
      p[p > 1] <- 1

      
      if (all(p == 0)) {
        mvment <- 0
      } else {
        mvment <- sample.int(5, 1, prob = p)
      }
      
      # If the designated position is not occupied, the cell will move to that location.
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
    
    # Update the coordinates 
    # ind.position <- matrix(0, nrow = space.length.y, ncol = space.length.x)
    
    # for (z in 1:length(x.coord)) {
    #  ind.position[x.coord[z], y.coord[z]] <- 1
    # }
    
    # Plots & summarising the densities
    
    #if (i %% (int.timesteps * 0.864) == 0) {
    #plot(x[y.coord[1]],y[space.length.y + 1 - x.coord[1]], 
    #     xlim = c(x[1], x[length(x)]), ylim = c(y[1], y[length(y)]), 
    #     pch = 21, bg = "blue", cex = 2, xlab = "x", 
    #     ylab = "y", main = paste("Invasion pattern after ", 
    #                              i/(int.timesteps * 0.864), " day", sep = ""))
    
    #for (k in 2:length(x.coord)) {
    #  points(x[y.coord[k]], y[space.length.y + 1 - x.coord[k]],
    #         pch = 21, cex = 2, bg = "blue")
    #}
    
    if (i == (day.timesteps * 3)) {
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
    #}
    # Optional, can be used to track the progress of the simulation
    # print(i)
  }
  
  return(list(den.mat.d12 = density.mat.d12, ind.position.d12 = ind.position.d12,
              n.den.d12 = n.day12,
              f.den.d12 = f.day12, 
              m.den.d12 = m.day12
  ))
}

calculate.sse <- function(pars, prob.death, prob.prof) {
  # Purpose: calculates the least square difference between the simulated
  # invasion pattern and the reference invasion pattern.
  # Arguments: 
  #    pars: parameters of the PDE model
  #    init.cells.cols: number of cells columns being placed at the left
  #    of the domain initially.
  #    prob.death: proportion of cells that will extinct at the end of the
  #                day.
  #    prob.prof: proportion of cells that will proliferate at the end of the
  #               day.
  temp.den.table <- generate.pattern(par = pars,
                                     prob.death = prob.death, prob.prof = prob.prof)
  
  if (is.na(temp.den.table)) {
    return(list(diff = NaN))
    # If the result from "generate.pattern" is NaN, return NaN.
  } else {
    # I was trying the same trick I used in the 1D case - minimize the discrepancy
    # at t = 3 first, then bring t = 6 into account and eventually all the 5 datasets,
    # but this is going to be extremely time consuming as we can expect...
    
    sse.t12 <- sum((temp.den.table$den.mat.d12 - t12.ref.den)^2)
    
    return(list(diff = sse.t12,
                den.mat.d12 = temp.den.table$den.mat.d12, 
                ind.position.d12 = temp.den.table$ind.position.d12,
                n.den.d12 = temp.den.table$n.den.d12,
                f.den.d12 = temp.den.table$f.den.d12, 
                m.den.d12 = temp.den.table$m.den.d12
                ))
  }
}

# The ABC scheme
calculate.bw <- function(ss.mat, lb.bw, ub.bw, ess.target, step.size) {
  #ss.mat <- read.table("Round 13 Least Square.txt", sep = "",
  #                     header = TRUE)
  #  lb.bw <- 0
  #  ub.bw <- 5
  #  ess.target <- 2000
  
  ss.mat <- as.matrix(ss.mat)
  ind.nan <- which(is.na(ss.mat[,2]))
  if (length(ind.nan) == 0) {
    ss.mat.valid <- ss.mat
  } else {
    ss.mat.valid <- ss.mat[-ind.nan, ]
  }
  
  power <- seq(lb.bw, ub.bw, by = step.size)
  ess.vec <- rep(0, length = length(power))
  
  for (i in 1:length(power)) {
    wt.temp <- ss.mat.valid[,2]^(-power[i])
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
    
    ess.temp <- ((sum(resamp.prob.temp))^2)/sum(resamp.prob.temp^2)
    ess.vec[i] <- ess.temp
    print(i)
  }
  
  ess.mat <- cbind(power, ess.vec)
  ind.nan.ess <- which(is.na(ess.mat[,2]))
  
  if (length(ind.nan.ess) == 0) {
    ess.mat.valid <- ess.mat
  } else {
    ess.mat.valid <- ess.mat[-ind.nan.ess,]
  }
  
  ess.obj.ind <- which(abs(ess.mat.valid[,2] - ess.target) == 
                         min(abs(ess.mat.valid[,2] - ess.target)))
  ess.obj <- ess.mat.valid[ess.obj.ind,2]
  bw.obj <- ess.mat.valid[ess.obj.ind,1]
  
  
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
  
  info.list <- list(info.mat = cbind(ss.mat.valid, wt.obj, resamp.prob.obj),
                    ess.mat = ess.mat,
                    ess.obj = ess.obj,
                    bw.obj = bw.obj)
  
  return(info.list)
}

# The ABC scheme
abc_bcd <- function(info.mat,paras) {
  # Purpose: Given the current set of parameters and the corresponding results of 
  # summary statistics, generate a new set of parameters. 
  
  paras_lb <- c(0.000069, 0.005, 0.0008, 7, 0.0001, 0.07, 0.01, 0.2)
  paras_ub <- c(0.02, 0.26, 0.08, 18, 0.033, 0.18, 0.1, 1)
  
  resamp_ind <- sample(info.mat[,1],size = length(paras[,1]), 
                       replace = TRUE, prob = info.mat[,length(info.mat[1,])]) # Resample the indices. 
  
  paras_nr_unperturbed <- paras[resamp_ind,] # Resampled parameter vectors, without perturbation.
  
  paras_nr_perturbed <- matrix(0,nrow = nrow(paras),ncol = ncol(paras)) # An empty matrix used to store the perturbed parameter values.
  
  #set.seed(123)
  
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
  }  # Perturbation, see the manuscript for more details.
  
  # paras_nr_perturbed <- as.data.frame(paras_nr_perturbed) # Change the format back to a data frame. 
  
  return(paras_nr_perturbed) # Parameter values for the next round. 
}

reduction.disp <- function(prev, curr) {
  per <- (prev-curr)/prev*100
  return(per)
}