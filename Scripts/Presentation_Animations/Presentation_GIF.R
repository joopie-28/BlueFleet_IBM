#########################################################
#### Generating GIF's to communicate model equations ####
#### in Presentation ####################################
library(animation)


# Movement of Physalia
physaliaMovement <- function(physalia, glaucus){
  
  # Boundary conditions for movement, i.e. out of frame
  if(physalia$y >= 100 | physalia$y <= 1 |physalia$x >= 100){
    
    return(physalia)}
  
  
  # Live condition
  if(physalia$x > 1 & physalia$status != 'EATEN'){
    
    
    ### Predator module
    
    # We want to simulate damage to the physalia, i.e. being killed by a predator.
    # We use a simple rule: if the glaucus is within feeding range of the physalai
    # for 3 timesteps, the physalia has been eaten completely. The glaucus will then also
    # move on as we delete the physalia. 
    
    # Convert to spatial geometries to allow geometric operations.
    spat.point <- st_point(c(physalia$x, physalia$y))
    
    # We can decide what a reasonable buffer is. This is the same for
    # all physalia. 2m
    under.attack.zone <- st_buffer(spat.point, 0.002)
    
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
    
    
    
    
    
    
  } else{
    
    if(physalia$x <= 1){
      physalia$status <- 'BEACHED'
    }
  }
  return(physalia)
}

# Movement of Glaucus with extra interaction. -extended boundary
glaucusMovement <- function(glaucus, physalia){
  # Boundary conditions for movement, i.e. out of frame
  if(glaucus$y >= 100 | glaucus$y <= 1 | glaucus$x >= 100){ 
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
      
      
      
      # Now account for effect of current and wind on glaucus.  
      
      # Positional update rules: this IS the movement of a Glaucus
      glaucus$x <- glaucus$x + 0.001*wind_strength[ceiling(glaucus$y), 
                                                   ceiling(glaucus$x)] * sin(wind_direction[ceiling(glaucus$y), 
                                                                                            ceiling(glaucus$x)])+
        0.01 * current_strength[ceiling(glaucus$y), 
                                 ceiling(glaucus$x)] * sin(current_direction[ceiling(glaucus$y), 
                                                                             ceiling(glaucus$x)])
      
      # Y movement (north - south) uses cosine function
      glaucus$y <- glaucus$y + 0.001*wind_strength[ceiling(glaucus$y), 
                                                   ceiling(glaucus$x)] * cos(wind_direction[ceiling(glaucus$y), 
                                                                                            ceiling(glaucus$x)]) +
        0.01 * current_strength[ceiling(glaucus$y), 
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
        glaucus.target.id <- physalia.df.spat[[which.min(st_distance(spat.point, physalia.df.spat)),'ID']]
        
        # Now make it move towards the target. We want to use a random walk-esque
        # movement with bias towards the physalia.
        
        target_angle <- atan2(glaucus.target[2] - glaucus$y, glaucus.target[1] - glaucus$x) # atan2 calculates the angle to get from y to x
        glaucus$x <- glaucus$x + glaucus$speed * cos(target_angle) # switch????
        glaucus$y <- glaucus$y +  glaucus$speed * sin(target_angle)
        
        # track interaction - defined as entities being close enough to change behaviour
        glaucus$interaction <- as.numeric(glaucus$interaction) + 1
        
        # Latch-on module
        eating.zone <- st_buffer(spat.point, 0.002) # 2 meters is our spatial resolution
        if(st_intersects(physalia.df.spat[glaucus.target.id,], eating.zone, sparse = F)){
          print('LATCHED')
          # Identify the prey to latch on to
          glaucus$target_ID = glaucus.target.id 
          
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

# Custom functions for drawing pictograms
draw.Glaucus <- function(center, scaling){
  coords_body = as.matrix(data.frame('x'=(c(5,7,8.2,8.5,7,6,8,10,10,10,8,6,6,5,4,4,2,0,0,0,2,4,3,1.5,1.8,3,5)-5)*scaling + center[1],
                                     "y"=(c(2, 8,7,8,8,9,10.5,10,11,12,11,12,13.5,14,13.5,12,11,12,11,10,10.5,9,8,8,7,8,2)-14)*scaling + center[2]))
  coords_stripe = as.matrix(data.frame('x'=(c(5,4.5,4.5,4.2,5.8,5.5,5.5,5)-5)*scaling + center[1],
                                       "y"=(c(4.5,8,12.5,13,13,12.5,8,4.5)-14)*scaling+ center[2]))
  Glaucus_Body <- st_polygon(list(coords_body))
  Glaucus_Stripe <- st_polygon(list(coords_stripe))
  plot(Glaucus_Body, col = 'skyblue', add =T)
  plot(Glaucus_Stripe, col = 'blue', add = T)
}

draw.Physalia<- function(center, orientation,scaling=1){
  if(orientation == 'left'){
    coords_body = as.matrix(data.frame('x'=(c(4,4.5, 5, 5.5, 6.3, 6.6, 6.3, 6, 4.8, 4)-5.8)*scaling ,
                                       "y"=(c(5,6, 6.2, 6.2, 6, 5.9, 5.4,4.9, 4.8,5)-6)*scaling))
  }else{
    coords_body = as.matrix(data.frame('x'=((5-(c(4,4.5, 5, 5.5, 6.3, 6.6, 6.3, 6, 4.8,4)-5)-5.8))*scaling,
                                       "y"=(c(5,6, 6.2, 6.2, 6, 5.9, 5.4,4.9, 4.8,5)-6)*scaling))
  }
  
  phys_Body <- st_polygon(list(coords_body))
  
  
  for (i in seq(0,1.8,0.2)){
    phys_Body <-  as.matrix(data.frame('x'=(((c(4,4,4.2,4))+i-5.8)*scaling),
                                       
                                       'y'=(c(5,4,3.6,5)+0.2-6)*scaling)) |>
      list()|>
      st_polygon() |>
      st_union(phys_Body)
  }
  
  phys_Body[[1]][,1] <- phys_Body[[1]][,1] + center[1]
  phys_Body[[1]][,2] <- phys_Body[[1]][,2] + center[2]
  
  plot(phys_Body, col = 'blue', add =T)
}


#### GIF 1. Movement due to physical forces ####

# Wind and current speeds (in km/H).
# Param recommended options:
# Grid spatial extent, each unit represents a kilometer
n_cols <- 100
n_rows <- 100

strength_current = 25
strength_wind = 15

current_strength <- matrix(rep(strength_current,n_rows*n_cols), 
                           nrow=n_rows)
wind_strength <- matrix(rep(strength_wind,n_rows*n_cols), 
                        nrow=n_rows)

# Wind and current directions

wind_direction <- matrix(rep(dir_wind, n_rows*n_cols), 
                         nrow=n_rows, ncol=n_cols)
current_direction <- matrix(rep(dir_current, n_rows*n_cols), 
                            nrow=n_rows, ncol=n_cols)

# Wind and current directions
# We use radians to denote directions. Feasible parameter values are therefore:
# any number between 1 and 2, where:
# 'west' = 1.5 pi
# 'east' = .5*pi
# 'north' = 0
# 'south' = pi

dir_wind = 1.5*pi
dir_current = pi
glaucus_Chemodetection <- 0.5
glaucus_Speed <- 0.00015


nGlaucus <- 1
glaucus <- data.frame(matrix(nrow=nGlaucus, ncol=10))
colnames(glaucus) <- c("ID", "x", "y", "chemodetection", "speed", "latch_time","total_latch_time", "target_ID", "interaction", "status")
for (i in 1:nGlaucus) {
  glaucus[i,] <- list(
    ID = i,
    x = 8.96,
    y = 5.88,
    chemodetection = glaucus_Chemodetection, 
    speed = glaucus_Speed,
    latch_time = 0,
    total_latch_time = 0,
    target_ID = NA,
    interaction = 0,
    status = 'ALIVE'
  )
}

# Set up bluebottle movement and functions. Bluebottles can be right
# or left handed.
nPhysalia <- 5
physalia <- data.frame(matrix(nrow=nPhysalia, ncol=6))
colnames(physalia) <- c("ID", "x", "y", "underattack", "status", "orientation")

for (i in 1:nPhysalia) {
  physalia[i,] <- list(
    ID = i,
    x = round(runif(1, 8.8, 9),2),
    y = round(runif(1,5.5, 6),2),
    underattack = 0,
    status = 'ALIVE',
    orientation = sample(c('left', 'right'),1)
  )
}

# Movie
n_frames <- 100
interval = 0.1
fps <- 5

# create GIF
saveGIF({
  for(i in 1:100){
    plot(1:105 ,1:105, type = 'n',
         xlab = 'West - East ',
         ylab = 'South - North', xlim=c(0,100))
    text(x=20, y= 3,paste0('Time (Hours): ', i))
    
    for(t in seq(10,100,10)){
      arrows(x0=seq(100,5, by=-10), x1=seq(95,0, by=-10), y0 = 110-t, y1=110-t, col = rgb(0,0.5,1,0.8),length=0.1)
    }
    for(t in seq(10,100,10)){
      arrows(y0=seq(100,5, by=-10), y1=seq(95,0, by=-10), x0 = 112.5-t, x1=112.5-t, col = rgb(0.5,0.3,0.8,0.5),length=0.1)
    }
    for(k in 1:nPhysalia){
      physalia[k,] <- physaliaMovement(physalia[k,], glaucus)
     # draw.Physalia(c(physalia$x[k], physalia$y[k]), orientation = physalia$orientation[k])
    }
    for(j in 1:nGlaucus){
      glaucus[j,] <- glaucusMovement(glaucus[j,], physalia)
      #draw.Glaucus(c(glaucus$x[j], glaucus$y[j]), scaling=0.5)
    }
  }
}, movie.name = "./Presentation/Movement.gif",interval = interval, ani.width = 500, ani.height = 500, fps = fps)


#### GIF 2. Hunting Behaviour ####

# Wind and current speeds (in km/H).
# Param recommended options:
# Grid spatial extent, each unit represents a kilometer
n_cols <- 100
n_rows <- 100
 


strength_current = 0
strength_wind = 4 

current_strength <- matrix(rep(strength_current,n_rows*n_cols), 
                           nrow=n_rows)
wind_strength <- matrix(rep(strength_wind,n_rows*n_cols), 
                        nrow=n_rows)
#Wind and current directions
# We use radians to denote directions. Feasible parameter values are therefore:
# any number between 1 and 2, where:
# 'west' = 1.5 pi
# 'east' = .5*pi
# 'north' = 0
# 'south' = pi

dir_wind = 1.5*pi
dir_current = pi
# Wind and current directions

wind_direction <- matrix(rep(dir_wind, n_rows*n_cols), 
                         nrow=n_rows, ncol=n_cols)
current_direction <- matrix(rep(dir_current, n_rows*n_cols), 
                            nrow=n_rows, ncol=n_cols)


nGlaucus <- 1

glaucus <- data.frame(matrix(nrow=nGlaucus, ncol=10))
colnames(glaucus) <- c("ID", "x", "y", "chemodetection", "speed", "latch_time","total_latch_time", "target_ID", "interaction", "alive")
for (i in 1:nGlaucus) {
  glaucus[i,] <- list(
    ID = i,
    x = 52,
    y = 63,
    chemodetection = 5, 
    speed = .1,
    latch_time = 0,
    total_latch_time = 0,
    target_ID = NA,
    interaction = 0,
    status = 'ALIVE'
  )
}
# Set up bluebottle movement and functions. Bluebottles can be right
# or left handed.
nPhysalia <- 1

physalia <- data.frame(matrix(nrow=nPhysalia, ncol=6))
colnames(physalia) <- c("ID", "x", "y", "underattack", "status", "orientation")

for (i in 1:nPhysalia) {
  physalia[i,] <- list(
    ID = i,
    x = 58,
    y = 56.5,
    underattack = 0,
    status = 'ALIVE',
    orientation = 'right'
  )
}

saveGIF({
  for(i in 1:100){
    plot(1:100 ,1:100, type = 'n',
         xlab = 'West - East ',
         ylab = 'South - North',xlim=c(45,60), ylim=c(50,65))
    text(x=0.0010, y= 0.0005,paste0('Time (Hours): ', i))
    
 
    for(k in 1:nPhysalia){
      physalia[k,] <- physaliaMovement(physalia[k,], glaucus)
      draw.Physalia(c(physalia$x[k], physalia$y[k]), orientation = physalia$orientation[k], scaling = .2)
    }
    for(j in 1:nGlaucus){
      glaucus[j,] <- glaucusMovement(glaucus[j,], physalia)
      draw.Glaucus(c(glaucus$x[j], glaucus$y[j]), scaling = 0.1)
      detection.zone <- st_buffer(st_point(c(glaucus$x[j], glaucus$y[j])), glaucus$chemodetection[j])
      colvec <- ifelse(st_within(st_point(c(physalia$x[k], physalia$y[k])), detection.zone, sparse = F), rgb(0,0,1,0.3), 
                       rgb(0,0,1,0))
      plot(detection.zone, col = colvec, lty=2, add = T)
      
    }
  }
}, movie.name = "./Presentation/Single_hunt.gif",interval = interval, ani.width = 500, ani.height = 500, fps = fps)


#### GIF 3. 'Stochastic directions ####

# Wind and current speeds (in km/H).
# Param recommended options:
# Grid spatial extent, each unit represents a kilometer
n_cols <- 100
n_rows <- 100

strength_current = 15
strength_wind = 10

# Wind and current directions
# We use radians to denote directions. Feasible parameter values are therefore:
# any number between 1 and 2, where:
# 'west' = 1.5 pi
# 'east' = .5*pi
# 'north' = 0
# 'south' = pi

current_strength <- matrix(rep(strength_current,n_rows*n_cols), 
                           nrow=n_rows)
wind_strength <- matrix(rep(strength_wind,n_rows*n_cols), 
                        nrow=n_rows)
#Wind and current directions
# We use radians to denote directions. Feasible parameter values are therefore:
# any number between 1 and 2, where:
# 'west' = 1.5 pi
# 'east' = .5*pi
# 'north' = 0
# 'south' = pi

dir_wind = 1.5*pi
dir_current = pi

nGlaucus <- 10
glaucus <- data.frame(matrix(nrow=nGlaucus, ncol=8))
colnames(glaucus) <- c("ID", "x", "y", "chemodetection", "speed", "latch_time", "target_ID", "status")
for (i in 1:nGlaucus) {
  glaucus[i,] <- list(
    ID = i,
    x = runif(1,80, 90),
    y = runif(1,50, 60),
    chemodetection = glaucus_Chemodetection, 
    speed = rnorm(1,glaucus_Speed, sd = glaucus_Speed/4),
    latch_time = 0,
    target_ID = NA,
    status = 'ALIVE'
  )
}

# Set up bluebottle movement and functions. Bluebottles can be right
# or left handed.
nPhysalia <- 50
physalia <- data.frame(matrix(nrow=nPhysalia, ncol=6))
colnames(physalia) <- c("ID", "x", "y", "underattack", "status", "orientation")

for (i in 1:nPhysalia) {
  physalia[i,] <- list(
    ID = i,
    x = runif(1, 80, 90),
    y = runif(1,50, 60),
    underattack = 0,
    status = 'ALIVE',
    orientation = sample(c('left', 'right'),1)
  )
}

rad2deg <- function(rad) {(rad * 180) / (pi)}
saveGIF({
  for(i in 1:100){
    print(i)
    plot(1:105 ,1:105, type = 'n',
         xlab = 'West - East ',
         ylab = 'South - North', xlim=c(0,100), ylim=c(0,100))
    text(x=20, y= 3,paste0('Time (Hours): ', i))
    #Wind and current directions vary each hour, but within boundaries
    # Wind and current directions vary each hour, but within boundaries
    dir_wind = ifelse(rnorm(1, dir_wind, 0.1) > 3*pi, 3*pi, 
                      ifelse(rnorm(1, dir_wind, 0.1) < 0*pi, 0*pi,
                             rnorm(1, dir_wind, 0.1)))
    
    dir_current = ifelse(rnorm(1, dir_current, 0.1) > 3.5, 3.5, 
                         ifelse(rnorm(1, dir_current, 0.1) < 2.9, 2.9,
                                rnorm(1, dir_current, 0.1)))
    
    wind_direction <- matrix(rep(dir_wind, n_rows*n_cols), 
                             nrow=n_rows, ncol=n_cols)
    current_direction <-  matrix(rep(dir_current, n_rows*n_cols), 
                                 nrow=n_rows, ncol=n_cols)
           
    for(t in seq(0,100,10)) {                       
      arrows(x0= seq(0,100,5), x1 = seq(0,100,5) + 3* sin(dir_wind), y0=rep(100-t,20), y1=rep(100-t+3*cos(dir_wind), 20), col = rgb(0,0.5,1,0.8),length=.1)
    }
    for(t in seq(0,100,10)) { 
      arrows(x0= seq(5,95,5), x1 =seq(5,95,5) + 3* sin(dir_current), y0=rep(95-t,19), y1=rep(95-t+3*cos(dir_current),19), col = rgb(0.5,0.3,0.8,0.5),length=.1)
    }
    
    for(j in 1:nGlaucus){
      glaucus[j,] <- glaucusMovement(glaucus[j,], physalia)
      # draw.Glaucus(c(glaucus$x[j], glaucus$y[j]), scaling = .7)
    }
    for(k in 1:nPhysalia){
      physalia[k,] <- physaliaMovement(physalia[k,], glaucus)
     # draw.Physalia(c(physalia$x[k], physalia$y[k]), orientation = physalia$orientation[k])
    }
  
  }
}, movie.name = "Stochastic_Directions.gif",interval = interval, ani.width = 500, ani.height = 500, fps = fps)



#### Wind strength plot####

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


pdf('./Presentation/wind_trajectory.pdf',
    width=10, height=8)
plot(wind_strength_calc(c(1:72)), ylab = 'Wind Speed (km/H)', xlab = 'Time (H)', type='l')
dev.off()












