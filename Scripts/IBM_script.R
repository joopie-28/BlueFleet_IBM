####################################################
#### Predator, prey and strandings: An IBM #########
#### approach to modelling the 'Blue Fleet' ########
####################################################

### QBIO7004 Final Project - J.M. Sassen 46339212 ##

# Install packages
library(tidyverse)
library(sf)



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


#### Step 2. Set up the Agents (Bluebottles and Dragons) and movement rules ####

# We need to convert our lists to dataframes to make use of vectorization

# Glaucus atlanticus individuals
nGlaucus <- 2
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

# Set up bluebottle movement and functions

nPhysalia <- 30
physalia <- list()
for (i in 1:nPhysalia) {
  physalia[[i]] <- list(
    x = runif(1, 1, n_cols),
    y = runif(1, 1, n_rows),
    orientation = sample(c('left', 'right'),1)
  )
}

# Each entity needs movement rules, we set those up here. 
# could maybe go into a function, 'glaucusmove'
glaucusMovement <- function(glaucus){

  if(glaucus$x > 200){
    glaucus$col <- 'steelblue'
    
    # Positional update rules: this IS the movement of a Glaucus
    glaucus$x <- max(glaucus$x + wind_strength[round(glaucus$x,digits = 2), 
                                                         round(glaucus$y, digits = 2)] * sin(wind_direction[round(glaucus$x,digits = 2), 
                                                                                                                 round(glaucus$y, digits = 2)]) +
                            current_strength[round(glaucus$x,digits = 2), 
                                             round(glaucus$y, digits = 2)] * sin(current_direction[round(glaucus$x,digits = 2), 
                                                                                                        round(glaucus$y, digits = 2)]),0)
    
    # Y movement (north - south) uses cosine function
    glaucus$y <- max(glaucus$y + wind_strength[round(glaucus$x,digits = 2), 
                                                         round(glaucus$y, digits = 2)] * cos(wind_direction[round(glaucus$x,digits = 2), 
                                                                                                                 round(glaucus$y, digits = 2)]) +
                            current_strength[round(glaucus$x,digits = 2), 
                                             round(glaucus$y, digits = 2)] * cos(current_direction[round(glaucus$x,digits = 2), 
                                                                                                        round(glaucus$y, digits = 2)]),0)
    
  } else{
    glaucus$col <- 'red'
    glaucus$status <- 1 # status 1 is beached
  }
  return(glaucus)
}

# Physalia movement
physaliaMovement <- function(physalia){
  
  # Boundary conditions
  if(physalia$x >= 1000 | physalia$y >= 1000 | physalia$y <= 0){ 
    physalia$col <- 'green' 
    physalia$status <- 2
    return(physalia)}
  
  # Live condition
    if(physalia$x > 200){
    physalia$col <- 'purple'
    
    # This is a very important part. We have right and left=handed bluebottles
    # They drift in opposite directions - presumably to sustain populations.
    # We need to account for this properly. In addition, we want to add some
    # stochasticity to the movement. This is due to inherent variability
    # in the shape and size of bluebottles, but also due to waves etc.
    
    if(physalia$orientation == 'right') {direction_offset <- rnorm(1,1,0.5)*pi/3} # right-handed drift at 50 degrees from wind direction.
    if(physalia$orientation == 'left') {direction_offset <- rnorm(1,-1,0.5)*pi/3} # left handed drift at

    # Positional update rules: this IS the movement of physalia
    # Wind larger impact on physalia due to sail
    # physalia also have an offset - their sails change the way they interact
    # with wind - this is one of the cool parts in the model.
    physalia$x <- max(physalia$x + 2*wind_strength[round(physalia$x,digits = 2), 
                                                   round(physalia$y, digits = 2)] * sin(wind_direction[round(physalia$x,digits = 2), 
                                                                                                       round(physalia$y, digits = 2)]+ direction_offset) +
                        current_strength[round(physalia$x,digits = 2), 
                                         round(physalia$y, digits = 2)] * sin(current_direction[round(physalia$x,digits = 2), 
                                                                                                round(physalia$y, digits = 2)]),0)
    
    # Y movement (north - south) uses cosine function
    physalia$y <- max(physalia$y + 2*wind_strength[round(physalia$x,digits = 2), 
                                                   round(physalia$y, digits = 2)] * cos(wind_direction[round(physalia$x,digits = 2), 
                                                                                                       round(physalia$y, digits = 2)] + direction_offset) +
                        current_strength[round(physalia$x,digits = 2), 
                                         round(physalia$y, digits = 2)] * cos(current_direction[round(physalia$x,digits = 2), 
                                                                                                round(physalia$y, digits = 2)]),0)
    
  } else{
    physalia$col <- 'red'
    physalia$status <- 1 # status 1 is beached
  }
  return(physalia)
}

plot(1:1000, 1:1000, type = 'n',
     xlab = '<- West - East ->',
     ylab = '<- North - South ->')
abline(v=200, lty=2, col='red')
for(i in 1:1000){
  for(j in 1:nGlaucus){
    glaucus[[j]] <- glaucusMovement(glaucus[[j]])
    points(x=glaucus[[j]]$x,
           y=glaucus[[j]]$y,
           col = glaucus[[j]]$col,
           pch = 19,
           cex = 0.1)
  }
  for(k in 1:nPhysalia){
    physalia[[k]] <- physaliaMovement(physalia[[k]])
    points(x=physalia[[k]]$x,
           y=physalia[[k]]$y,
           col = physalia[[k]]$col,
           pch = 19,
           cex = 0.1)
  }
  
}

#### Step 3. Set up interactions between Glaucus and Physalia ####

# scan area for each glaucus
# Convert to spatial geometries to allow geometric operations.
spat.point <- st_point(c(glaucus$x, glaucus$y))

# We can decide what a reasonable buffer is. This might be a parameter tied to the 
# individual Glaucus!
detection.zone <- st_buffer(spat.point, 500)

# Find the Physalia that are within the detection zone.
# needs to be in df format!!!1 vectorise all at once.
for (individual in physalia){
  spat.point.physalia <- st_point(c(individual$x, individual$y))
  if(st_intersects(spat.point.physalia, detection.zone, sparse = F)){
    # might need to add the CLOSEST here.
    glaucus.target <- as.data.frame(st_coordinates(spat.point.physalia))
  }
  else{
    next
  }
}

# Now make it move towards the target. We want to use a random walk
# like movement with bias towards the physalia.
target_angle <- atan2(glaucus.target$Y - glaucus$y, glaucus.target$X - glaucus$x)


plot(spat.point, add=T, col ='blue', pch =19)
plot(spat.point.physalia,add=T, col ='purple', pch =19)
plot(detection.zone, add =T)

for(i in 1:1000){
  target_angle <- atan2(glaucus.target$Y - glaucus$y, glaucus.target$X - glaucus$x)
  glaucus$x <- glaucus$x+ rnorm(1, 1, 0.5) * cos(target_angle)
  glaucus$y <- glaucus$y +  rnorm(1, 1, 0.5) * sin(target_angle)
  points(x=glaucus$x,
         y=glaucus$y,
         col='blue',
         pch=19,
         cex =0.3)
}

