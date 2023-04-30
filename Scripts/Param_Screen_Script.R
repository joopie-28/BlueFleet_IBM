# Install packages
library(units, lib.loc=Sys.getenv("R_LIBS_USER"))
library(tidyverse)
library(sf)
library(data.table)

# Movement of Glaucus with extra interaction
glaucusMovement.ext <- function(glaucus, physalia){
  # Boundary conditions for movement, i.e. out of frame
  if(glaucus$y >= 100 | glaucus$y <= 1){ 
    return(glaucus)}
  
  # First check if glaucus are 'latched on' to physalia which will 
  # affect their movement.
  if(!is.na(glaucus$target_ID) & physalia[glaucus$target_ID,'status'][1] != 'EATEN'){
    glaucus$x <- physalia[glaucus$target_ID,'x']
    glaucus$y <- physalia[glaucus$target_ID,'y']
    glaucus$latch_time <- glaucus$latch_time + 1
    
    # Probability of glaucus detaching follows a beta prime distribution.
    # Shape of distribtution set to match expected behaviour.
    probs <- c(1-pbetapr(glaucus$latch_time, 24, 5), pbetapr(glaucus$latch_time, 24, 5))
    glaucus$target_ID <- sample(size=1, x = c(glaucus$target_ID, NA), prob = probs) # NA corresponds to letting go.
    
  }else{
    
    # Check the glaucus is not beached
    if(glaucus$x > 1){
      glaucus$col <- 'steelblue'
      
      
      ### Predator module
      
      # Glaucus are predators. They can detect prey from distance using chemical cues. We want to simulate these capabilities
      # by allowing glaucus limited movement to a target, if the target is within reasonable range.
      
      # scan area for each glaucus
      # Convert to spatial geometries to allow geometric operations.
      spat.point <- st_point(c(glaucus$x, glaucus$y))
      
      # We can decide what a reasonable buffer is. This might be a parameter tied to the 
      # individual Glaucus!
      detection.zone <- st_buffer(spat.point, glaucus$chemodetection)
      
      # Find the Physalia that are within the detection zone.
      
      physalia.df.spat <- st_as_sf(physalia , coords = c('x', 'y'))
      
      # Any prey in the detection zone?
      if(any(st_intersects(physalia.df.spat, detection.zone, sparse = F))){
        # Find the nearest Physalia
        print('Attack!')
        glaucus.target <- physalia.df.spat[[st_nearest_feature(spat.point, physalia.df.spat),'geometry']]
        # Now make it move towards the target. We want to use a random walk-esque
        # movement with bias towards the physalia.
        
        target_angle <- atan2(glaucus.target[2] - glaucus$y, glaucus.target[1] - glaucus$x) # atan2 calculates the angle to get from y to x
        glaucus$x <- glaucus$x + glaucus$speed * sin(target_angle)
        glaucus$y <- glaucus$y +  glaucus$speed * cos(target_angle)
        
        # Latch-on module
        eating.zone <- st_buffer(spat.point, 0.0005) # 0.5 meters is our spatial resolution
        if(any(st_intersects(physalia.df.spat, eating.zone, sparse = F))){
          # Identify the prey to latch on to
          glaucus$target_ID = which(st_intersects(physalia.df.spat, eating.zone, sparse = F) == T)
          glaucus$latch_time = 1
          
        }
      }
      
      
      # Now account for effect of current and wind on glaucus.  
      
      # Positional update rules: this IS the movement of a Glaucus
      glaucus$x <- max(glaucus$x + 0.001*wind_strength[round(glaucus$y,digits = 2), 
                                                       round(glaucus$x, digits = 2)] * sin(wind_direction[round(glaucus$y,digits = 2), 
                                                                                                          round(glaucus$x, digits = 2)]) +
                         0.001 *current_strength[round(glaucus$y,digits = 2), 
                                                 round(glaucus$x, digits = 2)] * sin(current_direction[round(glaucus$y,digits = 2), 
                                                                                                       round(glaucus$x, digits = 2)]),0)
      
      # Y movement (north - south) uses cosine function
      glaucus$y <- max(glaucus$y + 0.001*wind_strength[round(glaucus$y,digits = 2), 
                                                       round(glaucus$x, digits = 2)] * cos(wind_direction[round(glaucus$y,digits = 2), 
                                                                                                          round(glaucus$x, digits = 2)]) +
                         0.001 * current_strength[round(glaucus$y,digits = 2), 
                                                  round(glaucus$x, digits = 2)] * cos(current_direction[round(glaucus$y,digits = 2), 
                                                                                                        round(glaucus$x, digits = 2)]),0)
      
    } else{
      glaucus$status <- 'BEACHED' # status 1 is beached
    }
  }
  return(glaucus)
}

# Movement of Physalia
physaliaMovement <- function(physalia, glaucus){
  
  # Boundary conditions for movement, i.e. out of frame
  if(physalia$y >= 100 | physalia$y <= 1){ 
    return(physalia)}
  
  # Live condition
  if(physalia$x > 1 & physalia$status != 'EATEN'){
    physalia$col <- 'purple'
    
    # We have right and left=handed bluebottles
    # They drift in opposite directions - presumably to sustain populations.
    # We need to account for this properly. In addition, we want to add some
    # stochasticity to the movement. This is due to inherent variability
    # in the shape and size of bluebottles, but also due to waves etc.
    
    if(physalia$orientation == 'right') {direction_offset <- rnorm(1,1,0.15)*pi/3} # right-handed drift at 50 degrees from wind direction.
    if(physalia$orientation == 'left') {direction_offset <- rnorm(1,-1,0.15)*pi/3} # left handed drift at
    
    # Positional update rules: this IS the movement of physalia
    # Wind larger impact on physalia due to sail
    # physalia also have an offset - their sails change the way they interact
    # with wind - this is one of the cool parts in the model.
    # 0.0266 - see Lee, Schaeffer, Groeskamp (2021)
    physalia$x <- max(physalia$x + 0.0266*(wind_strength[round(physalia$y,digits = 2), 
                                                         round(physalia$x, digits = 2)]) * sin(wind_direction[round(physalia$y,digits = 2), 
                                                                                                              round(physalia$x, digits = 2)]+ direction_offset) +
                        0.001* current_strength[round(physalia$y,digits = 2), 
                                                round(physalia$x, digits = 2)] * sin(current_direction[round(physalia$y,digits = 2), 
                                                                                                       round(physalia$x, digits = 2)]),0)
    
    # Y movement (north - south) uses cosine function
    physalia$y <- max(physalia$y + 0.0266*(wind_strength[round(physalia$y,digits = 2), 
                                                         round(physalia$x, digits = 2)]) * cos(wind_direction[round(physalia$y,digits = 2), 
                                                                                                              round(physalia$x, digits = 2)] + direction_offset) +
                        0.001* current_strength[round(physalia$y,digits = 2), 
                                                round(physalia$x, digits = 2)] * cos(current_direction[round(physalia$y,digits = 2), 
                                                                                                       round(physalia$x, digits = 2)]),0)
    
    ### Predator module
    
    # We want to simulate damage to the physalia, i.e. being killed by a predator.
    # We use a simple rule: if the glaucus is within feeding range of the physalai
    # for 3 timesteps, the physalia has been eaten completely. The glaucus will then also
    # move on as we delete the physalia. 
    
    # Convert to spatial geometries to allow geometric operations.
    spat.point <- st_point(c(physalia$x, physalia$y))
    
    # We can decide what a reasonable buffer is. This is the same for
    # all physalia. 50 cm seems reasonable
    under.attack.zone <- st_buffer(spat.point, 0.0005)
    
    # Check if there is a predator nearby
    glaucus.df.spat <- st_as_sf(glaucus , coords = c('x', 'y'))
    
    # Any predators in the detection zone?
    if(any(st_intersects(glaucus.df.spat, under.attack.zone, sparse = F))){
      # If yes, the physalia is 'under attack'
      physalia$underattack <- physalia$underattack + 1
    }
    
    if(physalia$underattack >= 24){
      physalia$status <- 'EATEN'
    }
    
  } else{
    physalia$col <- 'red'
    if(physalia$x <= 1){
      physalia$status <- 'BEACHED'
    }
  }
  return(physalia)
}

# Main function that runs the simulation
simBlueFleet <- function(nTimes,n_rows, n_cols, nPhysalia, nGlaucus,
                         strength_current, strength_wind,
                         dir_wind, dir_current,
                         glaucus_Chemodetection, glaucus_Speed, iniSpace){
  
  # Assign wind and current speed and direction to each grid cell.
  # Based on param input
  current_strength <<- matrix(rep(strength_current,n_rows*n_cols), 
                             nrow=n_rows)
  wind_strength <<- matrix(rep(strength_wind,n_rows*n_cols), 
                          nrow=n_rows)
  
  # Wind and current directions
  
  wind_direction <<- matrix(rep(dir_wind, n_rows*n_cols), 
                           nrow=n_rows, ncol=n_cols)
  current_direction <<- matrix(rep(dir_current, n_rows*n_cols), 
                              nrow=n_rows, ncol=n_cols)
  
  
  ### Instantiate model agents
  
  # Glaucus atlanticus individuals
  print('Generating animals')
  # Our individuals also have attributes. Glaucus are our predators.
  # They will seek out Physalia, and are capable of (very limited)
  # powered movement. 
  
  glaucus <- data.frame(matrix(nrow=nGlaucus, ncol=8))
  colnames(glaucus) <- c("ID", "x", "y", "chemodetection", "speed", "latch_time", "target_ID", "status")
  for (i in 1:nGlaucus) {
    glaucus[i,] <- list(
      ID = i,
      x = runif(1, iniSpace$xmin, iniSpace$xmax),
      y = runif(1,iniSpace$ymin, iniSpace$ymax),
      chemodetection = glaucus_Chemodetection, 
      speed = rnorm(1,glaucus_Speed, sd = glaucus_Speed/4),
      latch_time = 0,
      target_ID = NA,
      status = 'ALIVE'
    )
  }
  
  # Set up bluebottle movement and functions. Bluebottles can be right
  # or left handed.
  physalia <- data.frame(matrix(nrow=nPhysalia, ncol=6))
  colnames(physalia) <- c("ID", "x", "y", "underattack", "status", "orientation")
  
  for (i in 1:nPhysalia) {
    physalia[i,] <- list(
      ID = i,
      x = runif(1, iniSpace$xmin, iniSpace$xmax),
      y = runif(1,iniSpace$ymin, iniSpace$ymax),
      underattack = 0,
      status = 'ALIVE',
      orientation = sample(c('left', 'right'),1)
    )
  }
  
  print('Simulating Dynamics')
  GlaucusSim <- list()
  PhysaliaSim <- list()
  
  # Run model over nTimes steps.
  for (i in 1:nTimes){
    
    # Update Physalia movement
    for(k in 1:nPhysalia){
      
      # This function updates the position of all simulated
      # Physalia. Also tracks their 'status', which may be one
      # of 'ALIVE', 'BEACHED', or 'EATEN'.
      physalia[k, ] <- physaliaMovement(physalia[k,], glaucus)
    }
    
    PhysaliaSim[[i]] <- physalia
    
    
    # Update Glaucus movement. Status includes 'BEACHED' or 'ALIVE'.
    for(j in 1:nGlaucus){
      glaucus[j,] <- glaucusMovement.ext(glaucus[j,], physalia)
      
    }
    GlaucusSim[[i]] <- glaucus
  }
  simResults <- list('GlaucusSim' = GlaucusSim,
                     'PhysaliaSim' = PhysaliaSim)
  return(simResults)
}


#### Model Analyses ####

# Initialise the starting area for the fleet. I
# did not think it made sense to just have them randomly
# throughout the ocean, they should be relatively close
# together to start. We can then see how they diverge.

iniSpace <- list(xmin = 8,
                 xmax = 9.5,
                 ymin = 90,
                 ymax = 95)

# Grid spatial extent, each unit represents a kilometer
n_rows <- 100
n_cols <- 10

# Number of simulation iterations
nTimes = 1000

# Number of Glaucus and Pysalia individuals
nPhysalia = 1000
nGlaucus = 500

# Direction parameters fixed
dir_wind = (45*pi)/180
dir_current = pi

# Dragon Chemodetection. The detection range for a predator.
# These are REALLY small animals, and they are limited in 
# movement. It really is an unknown.
# Parameter in meters. Sensible would be anything up to 1m (km units)
glaucus_Chemodetection <- 0.001

# Dragon speed. How fast can they swim? Again, not much is known
# Reasonable parameter values range between 1-25 cm an hour.
# We add some stochasticity to this in the model.
glaucus_Speed <- 0.00015

# We have a working model of the bahevorial ecology of two important animals in the blue fleet.
# What do we want to know?
# What output files do we need?
# Max beaching for different wind strengths
# Max beaching for different wind and direction (densityplot)

# How likely are interactions between glaucus and physalia for different densities? (What to do with wind)
# Is glaucus just as likely to strand, given wind should affect it less?


# Max strandings over parameter space
wind_strength_range <- seq(0,50,0.1)
current_strength_range <- seq(0,25,0.1)

max_strandings <- data.frame('strength_wind' = vector(mode='numeric', length(current_strength_range) * length(wind_strength_range)),
                             'direction_wind' = vector(mode='numeric',length(current_strength_range) * length(wind_strength_range)),
                             'current_strength' = vector(mode='numeric', length(current_strength_range) * length(wind_strength_range)),
                             'current_direction' = vector(mode='numeric', length(current_strength_range) * length(wind_strength_range)),
                             'beached_glaucus_max' = vector(mode='numeric', length(current_strength_range) * length(wind_strength_range)),
                             'beached_right_phys_max' = vector(mode='numeric', length(current_strength_range) * length(wind_strength_range)),
                             'beached_left_phys_max' = vector(mode='numeric',length(current_strength_range) * length(wind_strength_range)))
df_row = 0
for(strength_wind in wind_strength_range){
  for(strength_current in current_strength_range){
    print(strength_wind)
    df_row = df_row + 1
    BFS<-simBlueFleet(nTimes=nTimes, n_rows = n_rows, n_cols = n_cols, nPhysalia = nPhysalia, nGlaucus = nGlaucus,
                      strength_current = strength_current, strength_wind = strength_wind,
                      dir_wind = dir_wind, dir_current = dir_current, glaucus_Chemodetection = glaucus_Chemodetection,
                      glaucus_Speed = glaucus_Speed, iniSpace = iniSpace)
    
    beached.glaucus = vector(mode='numeric', length(BFS$GlaucusSim))
    beached.phys.left = vector(mode='numeric', length(BFS$GlaucusSim))
    beached.phys.right = vector(mode='numeric', length(BFS$GlaucusSim))
    
    # ACCOUNT FOR DIFFERENT STATUS LATER
    for(current.iter in 1:length(BFS$GlaucusSim)){
      beached.glaucus[current.iter] = sum(BFS$GlaucusSim[[current.iter]]$status == 'BEACHED') / nrow(BFS$GlaucusSim[[current.iter]])
      beached.phys.left[current.iter] = sum(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'left')$status == 'BEACHED') / nrow(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'left'))
      beached.phys.right[current.iter] = sum(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'right')$status == 'BEACHED')/ nrow(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'right'))
    }
    max_strandings[df_row,] <- c(strength_wind,
                                 dir_wind,
                                 strength_current,
                                 dir_current,
                                 max(beached.glaucus),
                                 max(beached.phys.right),
                                 max(beached.phys.left))
  }
}

# Return output
write.csv(max_strandings, paste0("BF_simulation_MaxStrandings.csv"))

