#########################################################
#### Generating GIF's to communicate model equations ####
#### in Presentation ####################################

#### GIF 1. Movement due to physical forces ####

# Wind and current speeds (in km/H).
# Param recommended options:
# Grid spatial extent, each unit represents a kilometer
n_cols <- 100
n_rows <- 100

strength_current = 2
strength_wind = 20

# Wind and current directions
# We use radians to denote directions. Feasible parameter values are therefore:
# any number between 1 and 2, where:
# 'west' = 1.5 pi
# 'east' = .5*pi
# 'north' = 0
# 'south' = pi

dir_wind = 1.5*pi
dir_current = pi

nGlaucus <- 1
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
nPhysalia <- 2
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


saveGIF({
  for(i in 1:100){
    plot(1:105 ,1:105, type = 'n',
         xlab = 'West - East ',
         ylab = 'South - North', xlim=c(5,100))
    text(x=20, y= 3,paste0('Time (Hours): ', i))
    
    for(t in seq(10,100,10)){
      arrows(x0=seq(100,5, by=-10), x1=seq(95,0, by=-10), y0 = 110-t, y1=110-t, col = rgb(0,0.5,1,0.8),length=0.1)
    }
    for(t in seq(10,100,10)){
      arrows(y0=seq(100,5, by=-10), y1=seq(95,0, by=-10), x0 = 112.5-t, x1=112.5-t, col = rgb(0.5,0.3,0.8,0.5),length=0.1)
    }
    for(k in 1:nPhysalia){
      physalia[k,] <- physaliaMovement(physalia[k,], glaucus)
      draw.Physalia(c(physalia$x[k], physalia$y[k]), orientation = physalia$orientation[k])
    }
    for(j in 1:nGlaucus){
      glaucus[j,] <- glaucusMovement.ext(glaucus[j,], physalia)
      draw.Glaucus(c(glaucus$x[j], glaucus$y[j]), scaling=0.5)
    }
  }
}, movie.name = "Hunting.gif",interval = interval, ani.width = 500, ani.height = 500, fps = fps)


#### GIF 2. Hunting Behaviour ####

# Wind and current speeds (in km/H).
# Param recommended options:
# Grid spatial extent, each unit represents a kilometer
n_cols <- 100
n_rows <- 100
 


strength_current = 3
strength_wind = 3

# Wind and current directions
# We use radians to denote directions. Feasible parameter values are therefore:
# any number between 1 and 2, where:
# 'west' = 1.5 pi
# 'east' = .5*pi
# 'north' = 0
# 'south' = pi

dir_wind = 1.5*pi
dir_current = pi

nGlaucus <- 1
glaucus <- data.frame(matrix(nrow=nGlaucus, ncol=8))
colnames(glaucus) <- c("ID", "x", "y", "chemodetection", "speed", "latch_time", "target_ID", "status")
for (i in 1:nGlaucus) {
  glaucus[i,] <- list(
    ID = i,
    x = 52,
    y = 63,
    chemodetection = glaucus_Chemodetection*5000, 
    speed = 0.1,
    latch_time = 0,
    target_ID = NA,
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
      glaucus[j,] <- glaucusMovement.ext(glaucus[j,], physalia)
      draw.Glaucus(c(glaucus$x[j], glaucus$y[j]), scaling = 0.1)
      detection.zone <- st_buffer(st_point(c(glaucus$x[j], glaucus$y[j])), glaucus$chemodetection[j])
      colvec <- ifelse(st_within(st_point(c(physalia$x[k], physalia$y[k])), detection.zone, sparse = F), rgb(0,0,1,0.3), 
                       rgb(0,0,1,0))
      plot(detection.zone, col = colvec, lty=2, add = T)
      
    }
  }
}, movie.name = "Hunting.gif",interval = interval, ani.width = 500, ani.height = 500, fps = fps)






#### GIF 3. 'Latchin on' behaviour ####

