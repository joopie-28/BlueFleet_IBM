###Plotting script ###

# Plot 1 Beaching rates

# Import output files
simOutput_list <- list()

for(i in 1:32){
  simOutput_list[[i]] <- read_csv(paste0("./Outputs/HPC_Outputs/PhysicalParamSims/BF_simulation_MaxStrandings",i,".csv"))[, c(2,3,4)]
}

# Add legend to top right, outside plot region
pdf('./Plots/glaucusbeaching.pdf')
plot(simOutput_list[[1]]$beached_glaucus, type = 'n',
     xlab = 'Time (Hours)', ylab = 'Proportion of Beached Glaucus', xlim=c(0,300),ylim=c(0,1))
for(i in 1:32){
 lines(simOutput_list[[i]]$beached_glaucus, col=alpha('steelblue',0.5), lwd = 2)
}
dev.off()

par(mar=c(5.1, 4.1, 4.1, 8.3), xpd=TRUE)
pdf('./Plots/physbeaching.pdf')
plot(simOutput_list[[1]]$beached.phys.left, type = 'n',ylim=c(0,1),
     xlab = 'Time (Hours)', ylab = 'Proportion of Beached Physalia', xlim=c(0,500))
for(i in 1:32){
  lines(simOutput_list[[i]]$beached.phys.left, col=alpha('purple',0.5), lwd = 2)
  lines(simOutput_list[[i]]$beached.phys.right, col=alpha('blue',0.5), lwd = 2)
}

legend("topright", legend=c("Left-Sailed","Right-Sailed"), lty=c(1,1), col = c('purple', 'blue'), title="Orientation", bg= 'white')
dev.off()

## Plot 2 hunting behaviour

huntOutput_list <- list()

for(i in 1:3){
  huntOutput_list[[i]] <- read_csv(paste0("./Outputs/HPC_Outputs/HuntingParamSims/Hunting_Simulations",i,".csv"))[,-1]
  huntOutput_list[[i]]$interactionPerInd <- huntOutput_list[[i]]$interaction/huntOutput_list[[i]]$nGlaucus
  huntOutput_list[[i]]$chemoM <- huntOutput_list[[i]]$glaucus_Chemodetection*1000
  huntOutput_list[[i]]$speedM <- huntOutput_list[[i]]$glaucus_Speed*1000
}

pdf('./Plots/Chemodetection.pdf',
    width=10, height=7)
plot(interactionPerInd~chemoM,data = huntOutput_list[[3]], pch=19, type = 'n',
     ylab= "Interactions per Individual", xlab = 'Chemodetection Range (m)')
cols <- c('blue', 'red', 'green')
for(i in 1:3){
 points(interactionPerInd~chemoM,data = huntOutput_list[[i]], bg = cols[i], pch=21, col='black')

}
legend('topleft',title = 'Prey Density' ,legend = c('1x','3x', '5x'), pch = 19,col=c('blue', 'green', 'red'))
dev.off()



pdf('./Plots/Swimming.pdf',
    width=10, height=7)
plot(interactionPerInd~speedM,data = huntOutput_list[[3]], pch=19, type = 'n',
     ylab= "Interactions per Individual", xlab = 'Swimming Speed (m)')
cols <- c('blue', 'red', 'green')
for(i in 1:3){
  points(interactionPerInd~speedM,data = huntOutput_list[[i]], bg = cols[i], pch=21, col='black')
  
}
legend('topright',title = 'Prey Density' ,legend = c('1x','3x', '5x'), pch = 19,col=c('blue', 'green', 'red'))

dev.off()





