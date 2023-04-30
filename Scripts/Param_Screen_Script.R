# Install packages
library(tidyverse)
library(sf)
library(data.table)

# Movement of Physalia
physaliaMovement <- function(physalia, glaucus){
  
  # Boundary conditions for movement
  if(physalia$x >= 100 | physalia$y >= 100 | physalia$y <= 0){ 
    physalia$col <- 'green' 
    physalia$status <- 2
    return(physalia)}
  
  # Live condition
  if(physalia$x > 1 & physalia$status != 'EATEN'){
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
    # 0.0266 - see Lee, Schaeffer, Groeskamp (2021)
    physalia$x <- max(physalia$x + 0.0266*(wind_strength[round(physalia$x,digits = 2), 
                                                         round(physalia$y, digits = 2)]) * sin(wind_direction[round(physalia$x,digits = 2), 
                                                                                                              round(physalia$y, digits = 2)]+ direction_offset) +
                        current_strength[round(physalia$x,digits = 2), 
                                         round(physalia$y, digits = 2)] * sin(current_direction[round(physalia$x,digits = 2), 
                                                                                                round(physalia$y, digits = 2)]),0)
    
    # Y movement (north - south) uses cosine function
    physalia$y <- max(physalia$y + 0.0266*(wind_strength[round(physalia$x,digits = 2), 
                                                         round(physalia$y, digits = 2)]) * cos(wind_direction[round(physalia$x,digits = 2), 
                                                                                                              round(physalia$y, digits = 2)] + direction_offset) +
                        current_strength[round(physalia$x,digits = 2), 
                                         round(physalia$y, digits = 2)] * cos(current_direction[round(physalia$x,digits = 2), 
                                                                                                round(physalia$y, digits = 2)]),0)
    
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
    glaucus.df <- rbindlist(glaucus, fill = T)
    glaucus.df.spat <- st_as_sf(glaucus.df , coords = c('x', 'y'))
    
    # Any predators in the detection zone?
    if(any(st_intersects(glaucus.df.spat, under.attack.zone, sparse = F))){
      # If yes, the physalia is 'under attack'
      physalia$underattack <- physalia$underattack + 1
    }
    
    if(physalia$underattack >= 3){
      physalia$status <- 'EATEN'
    }
    
  } else{
    physalia$col <- 'red'
    if(physalia$x < 1){
      physalia$status <- 'BEACHED' # status 1 is beached
    }
  }
  return(physalia)
}

# Movement of Glaucus
glaucusMovement <- function(glaucus, physalia){
  
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
    # needs to be in df format!!!1 vectorise all at once.
    
    physalia.df <- rbindlist(physalia, fill = T)
    
    physalia.df.spat <- st_as_sf(physalia.df , coords = c('x', 'y'))
    
    # Any prey in the detection zone?
    if(any(st_intersects(physalia.df.spat, detection.zone, sparse = F))){
      # Find the nearest Physalia
      print('Attack!')
      glaucus.target <- physalia.df.spat[[st_nearest_feature(spat.point, physalia.df.spat),'geometry']]
      # Now make it move towards the target. We want to use a random walk-esque
      # movement with bias towards the physalia.
      
      target_angle <- atan2(glaucus.target[2] - glaucus$y, glaucus.target[1] - glaucus$x) # atan2 calculates the angle to get from y to x
      glaucus$x <- glaucus$x + glaucus$speed * cos(target_angle)
      glaucus$y <- glaucus$y +  glaucus$speed * sin(target_angle)
      
    }
    
    
    # Now account for effect of current and wind on glaucus.  
    
    # Positional update rules: this IS the movement of a Glaucus
    glaucus$x <- max(glaucus$x + 0.005*wind_strength[round(glaucus$x,digits = 2), 
                                                     round(glaucus$y, digits = 2)] * sin(wind_direction[round(glaucus$x,digits = 2), 
                                                                                                        round(glaucus$y, digits = 2)]) +
                       current_strength[round(glaucus$x,digits = 2), 
                                        round(glaucus$y, digits = 2)] * sin(current_direction[round(glaucus$x,digits = 2), 
                                                                                              round(glaucus$y, digits = 2)]),0)
    
    # Y movement (north - south) uses cosine function
    glaucus$y <- max(glaucus$y + 0.005*wind_strength[round(glaucus$x,digits = 2), 
                                                     round(glaucus$y, digits = 2)] * cos(wind_direction[round(glaucus$x,digits = 2), 
                                                                                                        round(glaucus$y, digits = 2)]) +
                       current_strength[round(glaucus$x,digits = 2), 
                                        round(glaucus$y, digits = 2)] * cos(current_direction[round(glaucus$x,digits = 2), 
                                                                                              round(glaucus$y, digits = 2)]),0)
    
  } else{
    glaucus$col <- 'red'
    glaucus$status <- 'BEACHED' # status 1 is beached
  }
  return(glaucus)
}

#### Model Analyses ####

# Initialise the starting area for the fleet. I
# did not think it made sense to just have them randomly
# throughout the ocean, they should be relatively close
# together to start. We can then see how they diverge.

iniSpace <- list(xmin = 8,
                 xmax = 10,
                 ymin = 90,
                 ymax = 100)

# Grid spatial extent, each unit represents a kilometer
n_rows <- 100
n_cols <- 100

# Number of simulation iterations
nTimes = 1000

# Number of Glaucus and Pysalia individuals
nPhysalia = 1000
nGlaucus = 25

# Current parameters
dir_wind = 1.5*pi
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
wind_direction_range <- seq(0,2,0.05)

max_strandings <- data.frame('strength_wind' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)),
                             'direction_wind' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)),
                             'current_strength' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)),
                             'current_direction' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)),
                             'beached_glaucus_max' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)),
                             'beached_right_phys_max' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)),
                             'beached_left_phys_max' = vector(mode='numeric', length(wind_direction_range) * length(wind_strength_range)))
df_row = 0
for(strength.wind in wind_strength_range){
  for(dir.wind in wind_direction_range){
    df_row = df_row + 1
    BFS<-simBlueFleet(nTimes=nTimes, n_rows = n_rows, n_cols = n_cols, nPhysalia = nPhysalia, nGlaucus = nGlaucus,
                      strength_current = strength_current, strength_wind = strength_wind,
                      dir_wind = dir_wind, dir_current = dir_current, glaucus_Chemodetection = glaucus_Chemodetection,
                      glaucus_Speed = glaucus_Speed, iniSpace = iniSpace)
    
    beached.glaucus = vector(mode='numeric', length(BFS$GlaucusStats))
    beached.phys.left = vector(mode='numeric', length(BFS$GlaucusStats))
    beached.phys.right = vector(mode='numeric', length(BFS$GlaucusStats))
    
    # ACCOUNT FOR DIFFERENT STATUS LATER
    for(current.iter in 1:length(BFS$GlaucusStats)){
      beached.glaucus[current.iter] = sum(BFS$GlaucusStats[[current.iter]]$status == 'BEACHED') / nrow(BFS$GlaucusStats[[current.iter]])
      beached.phys.left[current.iter] = sum(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'left')$status == 'BEACHED') / nrow(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'left'))
      beached.phys.right[current.iter] = sum(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'right')$status == 'BEACHED')/ nrow(subset(BFS$PhysaliaSim[[current.iter]], orientation == 'right'))
    }
    max_strandings[df_row,] <- c(strength.wind,
                                 dir.wind,
                                 strength_current,
                                 dir_current,
                                 max(beached.glaucus),
                                 max(beached.phys.right),
                                 max(beached.phys.left))
  }
}

# Return output
write.csv(max_strandings, paste0("BF_simulation_MaxStrandings.csv"))

