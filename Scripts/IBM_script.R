####################################################
#### Predator, prey and strandings: An IBM #########
#### approach to modelling the 'Blue Fleet' ########
####################################################

### QBIO7004 Final Project - J.M. Sassen 46339212 ##

#### 1. Set up the Environment. ####

# The environment is a 100x100 grid. Each grid cell has certain
# properties, such as wind and current direction. Our first
# model iteration will keep this VERY simple. But we will 
# add stochasticity later.

# Spatial extent of grid
n_rows <- 1000
n_cols <- 1000

# We are interested in the behavioral ecology of this system,
# as well as the large-scale strandings. We thus need to add
# 'land' and 'water' attributes to our grid cells.

cell_type <- matrix(c(rep('Land',n_rows*n_cols*(1/n_cols)), rep('Water',n_rows*n_cols)),
                    nrow = n_rows, byrow = F) 

# Assign wind and current speed and direction to each grid cell.
# This is simple for now.

# Let's srength as constants for now.
current_strength <- matrix(rep(0.5,n_rows*n_cols), 
                           nrow=n_rows)
wind_strength <- matrix(rep(0.5,n_rows*n_cols), 
                           nrow=n_rows)

# Only to the left, or 'west'. So 1.5*pi gives us that direction in radians.
wind_direction <- matrix(rep(1.5*pi, n_rows*n_cols), 
                         nrow=n_rows, ncol=n_cols)
# For current, we will set them to the 'south' (primitively emulating the EAC)
current_direction <- matrix(rep(1*pi, n_rows*n_cols), 
                            nrow=n_rows, ncol=n_cols)


#### Step 2. Set up the Agents (Bluebottles and Dragons) ####

# Glaucus atlanticus individuals
nGlaucus <- 15
glaucus <- list()

# Our individuals also have attributes. Glaucus are our predators.
# They will seek out Physalia, and are capable of (very limited)
# powered movement. Status 0 = alive, 1 = beached, 2 = died natural
for (i in 1:nGlaucus) {
  glaucus[[i]] <- list(
    x = round(runif(1, 200, n_cols), digits = 0),
    y = round(runif(1, 1, n_rows),digits = 0),
    status = 0
  )
}


# Each entity needs movement rules, we set those up here. 
# could maybe go into a function, 'glaucusmove'
for (i in 1:nGlaucus) {
  
  if(glaucus[[i]]$x > 200){
    glaucus[[i]]$col <- 'blue'
    
    # Positional update rules: this IS the movement of a Glaucus
    glaucus[[i]]$x <- max(glaucus[[i]]$x + wind_strength[round(glaucus[[i]]$x,digits = 2), 
                                                         round(glaucus[[i]]$y, digits = 2)] * sin(wind_direction[round(glaucus[[i]]$x,digits = 2), 
                                                                                                                 round(glaucus[[i]]$y, digits = 2)]) +
                            current_strength[round(glaucus[[i]]$x,digits = 2), 
                                             round(glaucus[[i]]$y, digits = 2)] * sin(current_direction[round(glaucus[[i]]$x,digits = 2), 
                                                                                                        round(glaucus[[i]]$y, digits = 2)]),0)
    
    # Y movement (north - south) uses cosine function
    glaucus[[i]]$y <- max(glaucus[[i]]$y + wind_strength[round(glaucus[[i]]$x,digits = 2), 
                                                         round(glaucus[[i]]$y, digits = 2)] * cos(wind_direction[round(glaucus[[i]]$x,digits = 2), 
                                                                                                                 round(glaucus[[i]]$y, digits = 2)]) +
                            current_strength[round(glaucus[[i]]$x,digits = 2), 
                                             round(glaucus[[i]]$y, digits = 2)] * cos(current_direction[round(glaucus[[i]]$x,digits = 2), 
                                                                                                        round(glaucus[[i]]$y, digits = 2)]),0)
    
  } else{
    glaucus[[i]]$col <- 'red'
    glaucus[[i]]$status <- 1 # status 1 is beached
  }
}





  # Apply periodic boundary conditions
  if (glau[[i]]$x < 1) glau[[i]]$x <- n_cols
  if (glau[[i]]$x > n_cols) glau[[i]]$x <- 1
  if (glau[[i]]$y < 1) glau[[i]]$y <- n_rows
  if (glau[[i]]$y > n_rows) glau[[i]]$y <- 1







# Set up bluebottle movement and functions

n_physalia <- 10
physalia <- list()
for (i in 1:n_physalia) {
  physalia[[i]] <- list(
    x = runif(1, 1, n_cols),
    y = runif(1, 1, n_rows),
    buoyancy = rnorm(1, -1, 0.2),
    velocity = 0,
    direction = runif(1, 0, 2*pi)
  )
}
for (i in 1:n_physalia) {
  # Calculate velocity based on current speed and direction and buoyancy
  physalia[[i]]$velocity <- current_speed[round(physalia[[i]]$y), round(physalia[[i]]$x)]
  physalia[[i]]$direction <- current_direction[round(physalia[[i]]$y), round(physalia[[i]]$x)]
  physalia[[i]]$velocity <- physalia[[i]]$velocity + physalia[[i]]$buoyancy
  
  # Update position based on velocity
  physalia[[i]]$x <- physalia[[i]]$x + physalia[[i]]$velocity*cos(physalia[[i]]$direction)
  physalia[[i]]$y <- physalia[[i]]$y + physalia[[i]]$velocity*sin(physalia[[i]]$direction)
  
  # Apply periodic boundary conditions
  if (physalia[[i]]$x < 1) physalia[[i]]$x <- n_cols
  if (physalia[[i]]$x > n_cols) physalia[[i]]$x <- 1
  if (physalia[[i]]$y < 1) physalia[[i]]$y <- n_rows
  if (physalia[[i]]$y > n_rows) physalia[[i]]$y <- 1
}
