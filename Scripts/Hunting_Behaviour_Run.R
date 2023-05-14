#### Understanding hunting behaviour ###

# Install packages
library(units, lib.loc=Sys.getenv("R_LIBS_USER"))
library(rgeos, lib.loc=Sys.getenv("R_LIBS_USER"))
library(tidyverse)
library(sf)
library(data.table)
library(extraDistr, lib.loc=Sys.getenv("R_LIBS_USER"))

# step 1: obtain batch number from command line argument:
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) { # no arguments provided
  stop("Parameter for batch number needs to be provided.")
} else {
  batch <- as.integer(args[1])
}

# Use batch number to set different densities.
if(batch == 1){
  density<-1
}
if(batch==2){
  density <- 3
}
if(batch==3){
  density <- 5
}


#### 1. Define Functions ####

# Movement of Physalia
physaliaMovement <- function(physalia, glaucus){
  
  # Boundary conditions for movement, i.e. out of frame
  if(physalia$y >= 100 | physalia$y <= 1 |physalia$x >= 100){
    
    return(physalia)}
  
  
  # Live condition
  if(physalia$x > 1 & physalia$status != 'EATEN'){
    
    
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
    physalia$x <- physalia$x + 0.0266*(wind_strength[ceiling(physalia$y), 
                                                     ceiling(physalia$x)]) * sin(wind_direction[ceiling(physalia$y), 
                                                                                                          ceiling(physalia$x)]+ direction_offset) +
      0.001* current_strength[ceiling(physalia$y), 
                              ceiling(physalia$x)] * sin(current_direction[ceiling(physalia$y), 
                                                                                     ceiling(physalia$x)])
    
    # Y movement (north - south) uses cosine function
    physalia$y <- physalia$y + 0.0266*(wind_strength[ceiling(physalia$y), 
                                                     ceiling(physalia$x)]) * cos(wind_direction[ceiling(physalia$y), 
                                                                                                          ceiling(physalia$x)] + direction_offset) +
      0.001* current_strength[ceiling(physalia$y), 
                              ceiling(physalia$x)] * cos(current_direction[ceiling(physalia$y), 
                                                                                     ceiling(physalia$x)])
    
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
  
    if(physalia$x <= 1){
      physalia$status <- 'BEACHED'
    }
  }
  return(physalia)
}

# Movement of Glaucus with extra interaction
glaucusMovement <- function(glaucus, physalia){
  # Boundary conditions for movement, i.e. out of frame
  if(glaucus$y >= 100 | glaucus$y <= 1 | glaucus$x >= 10){ 
    return(glaucus)}
  
  # First check if glaucus are 'latched on' to physalia which will 
  # affect their movement.
  if(!is.na(glaucus$target_ID) & physalia[glaucus$target_ID,'status'][1] != 'EATEN'){
    glaucus$x <- physalia[glaucus$target_ID,'x']
    glaucus$y <- physalia[glaucus$target_ID,'y']
    glaucus$latch_time <- glaucus$latch_time + 1
    glaucus$total_latch_time <- glaucus$total_latch_time + 1
    # Probability of glaucus detaching follows a beta prime distribution.
    # Shape of distribtution set to match expected behaviour.
    probs <- c(1-pbetapr(glaucus$latch_time, 24, 5), pbetapr(glaucus$latch_time, 24, 5))
    glaucus$target_ID <- sample(size=1, x = c(glaucus$target_ID, NA), prob = probs) # NA corresponds to letting go.
    
  }else{
    
    
    # Check the glaucus is not beached
    if(glaucus$x > 1){
      glaucus$col <- 'steelblue'
      
      
      
      # Now account for effect of current and wind on glaucus.  
      
      # Positional update rules: this IS the movement of a Glaucus
      glaucus$x <- glaucus$x + 0.001*wind_strength[ceiling(glaucus$y), 
                                                   ceiling(glaucus$x)] * sin(wind_direction[ceiling(glaucus$y), 
                                                                                            ceiling(glaucus$x)])+
        0.001 * current_strength[ceiling(glaucus$y), 
                                 ceiling(glaucus$x)] * sin(current_direction[ceiling(glaucus$y), 
                                                                             ceiling(glaucus$x)])
      
      # Y movement (north - south) uses cosine function
      glaucus$y <- glaucus$y + 0.001*wind_strength[ceiling(glaucus$y), 
                                                   ceiling(glaucus$x)] * cos(wind_direction[ceiling(glaucus$y), 
                                                                                            ceiling(glaucus$x)]) +
        0.001 * current_strength[ceiling(glaucus$y), 
                                 ceiling(glaucus$x)] * cos(current_direction[ceiling(glaucus$y), 
                                                                             ceiling(glaucus$x)])
      
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
        glaucus.target <- physalia.df.spat[[which.min(st_distance(spat.point, physalia.df.spat)),'geometry']]
        
        # Now make it move towards the target. We want to use a random walk-esque
        # movement with bias towards the physalia.
        
        target_angle <- atan2(glaucus.target[2] - glaucus$y, glaucus.target[1] - glaucus$x) # atan2 calculates the angle to get from y to x
        glaucus$x <- glaucus$x + glaucus$speed * cos(target_angle) # switch????
        glaucus$y <- glaucus$y +  glaucus$speed * sin(target_angle)
        
        # Latch-on module
        eating.zone <- st_buffer(spat.point, 0.0005) # 0.5 meters is our spatial resolution
        if(any(st_intersects(physalia.df.spat, eating.zone, sparse = F))){
          # Identify the prey to latch on to
          glaucus$target_ID = which(st_intersects(physalia.df.spat, eating.zone, sparse = F) == T)[1] # problem, which min, gives double..
          glaucus$latch_time = 1
          
        }
      }
      
      
      
      
      
      
    } else{
      glaucus$status <- 'BEACHED' # status 1 is beached
    }
  }
  return(glaucus)
}

# Stochastic parameter wind strength
wind_strength_calc <-function(hour) {
  # Parameters for the sinusoidal function
  peak_hour <- 20        # Hour of the day with maximum wind
  amplitude <- 5       # Amplitude of the sinusoidal function (half of the variation)
  
  
  # Calculate daylight hours using a sinusoidal function
  strength = amplitude * sin((hour - peak_hour) * 2 * pi / 24) + 10 
  strength = (strength + rnorm(length(hour), 0, sd = 3))
  strength[which(strength<0)] <- 0
  return(strength)
}

# Simulation version that allows wind and current directions to vary stochastically at each timestep.
simBlueFleet.stochPhys <- function(nTimes,n_rows, n_cols, nPhysalia, nGlaucus,
                                   strength_current,
                                   dir_wind, dir_current,
                                   glaucus_Chemodetection, glaucus_Speed, iniSpace){
  
  # Assign wind and current speed and direction to each grid cell.
  # Based on param input
  current_strength <<- matrix(rep(strength_current,n_rows*n_cols), 
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
  
  glaucus <- data.frame(matrix(nrow=nGlaucus, ncol=9))
  colnames(glaucus) <- c("ID", "x", "y", "chemodetection", "speed", "latch_time","total_latch_time", "target_ID", "status")
  for (i in 1:nGlaucus) {
    glaucus[i,] <- list(
      ID = i,
      x = runif(1, iniSpace$xmin, iniSpace$xmax),
      y = runif(1,iniSpace$ymin, iniSpace$ymax),
      chemodetection = glaucus_Chemodetection, 
      speed = glaucus_Speed,
      latch_time = 0,
      total_latch_time = 0,
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

    wind_strength <<- matrix(rep(wind_strength_calc(i),n_rows*n_cols), 
                             nrow=n_rows)
    # Wind and current directions vary each hour, but within boundaries
    dir_wind = ifelse(rnorm(1, dir_wind, 0.1) > 1.75*pi, 1.75*pi, 
                      ifelse(rnorm(1, dir_wind, 0.1) < 1.25*pi, 1.25*pi,
                             rnorm(1, dir_wind, 0.1)))
    
    dir_current = ifelse(rnorm(1, dir_current, 0.1) > 3.5, 3.5, 
                         ifelse(rnorm(1, dir_current, 0.1) < 2.9, 2.9,
                                rnorm(1, dir_current, 0.1)))
    
    wind_direction <<- matrix(rep(dir_wind, n_rows*n_cols), 
                              nrow=n_rows, ncol=n_cols)
    current_direction <<- matrix(rep(dir_current, n_rows*n_cols), 
                                 nrow=n_rows, ncol=n_cols)
    
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
      glaucus[j,] <- glaucusMovement(glaucus[j,], physalia)
      
    }
    GlaucusSim[[i]] <- glaucus
  }
  simResults <- list('GlaucusSim' = GlaucusSim,
                     'PhysaliaSim' = PhysaliaSim)
  
  return(simResults)
}


### REMOVE PREY DENSITY AND DO THREE ITERATIONS INSTEAD 1, 3, 5

## Read in the pre-made LHS 
LHS<-read.csv('LHS_Hunting_Behaviour.csv')

# Number of predators
nGlaucus = 250

LHS$nGlaucus = nGlaucus

LHS$Total_Attacks <- 0
LHS$Total_LatchTime <- 0
LHS$Total_eaten <- 0
LHS$nPhysalia <- density*nGlaucus
LHS$density <- density

# Initialise the starting area for the fleet. I
# did not think it made sense to just have them randomly
# throughout the ocean, they should be relatively close
# together to start. We can then see how they diverge.

iniSpace <- list(xmin = 8,
                 xmax = 8.2,
                 ymin = 45,
                 ymax = 50)

# Grid spatial extent, each unit represents a kilometer
n_rows <- 100
n_cols <- 100

# Number of simulation iterations
nTimes = 1000

# Direction and strength parameters intialisation
dir_wind = (45*pi)/180
dir_current = pi
strength_current = 2

# Run 100 simulations based on the LHS

for(i in 1:nrow(LHS)){
  print(paste0(i, ' of 100'))
  # Define parameters
  nPhysalia = LHS$nPhysalia[i]
  glaucus_Chemodetection = LHS$glaucus_Chemodetection[i]
  glaucus_Speed = LHS$glaucus_Speed[i]
  
  BFS <- simBlueFleet.stochPhys(nTimes=nTimes, n_rows = n_rows, n_cols = n_cols, nPhysalia = nPhysalia, 
                                nGlaucus = nGlaucus, strength_current = strength_current,
                                dir_wind = dir_wind, dir_current = dir_current, 
                                glaucus_Chemodetection = glaucus_Chemodetection,
                                glaucus_Speed = glaucus_Speed, iniSpace = iniSpace)
  
  
  # how many attacks
  
  LHS$Total_Attacks[i] <- sum(BFS$PhysaliaSim[[nTimes]]$underattack)
  
  # how many eaten
  LHS$Total_eaten[i] <- sum(BFS$PhysaliaSim[[nTimes]]$status == 'EATEN')
  
  # latch_time
  LHS$Total_LatchTime[i] <- sum(BFS$GlaucusSim[[nTimes]]$total_latch_time)
  
}
  
  
# Return output
write.csv(LHS, paste0("Hunting_Simulations", batch, ".csv"))

  
  

