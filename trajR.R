library("trajr")
#https://cran.r-project.org/web/packages/trajr/trajr.pdf

#trajectory= set of 2-dimensional spatial coordinates with a third temporal dimension.
#OR series of steps, each comprised of a length (Li), a turning angle (Δi), and a time
#correlated random walks: turning angle of each step is the direction of the previous step +-some error
#directed walks: or compass-based navigation, in which the angular errors at each step are added to the “ideal” or compass direction

#from the list trial_ls, take the x element (x trial) and thenselect those colums
track_coord <- trial_ls[[8]] %>% 
  select(x, y, time)
trj <- TrajFromCoords(track_coord, fps = 30, spatialUnits = "cm")
plot(trj)

#TrajFromCoords assumes that the first column in coords contains x values, 
#the second contains y values, and that points were sampled at a constant sampling frequency of 50 points/second
#The help page for TrajFromCoords describes how to alter these defaults
?TrajFromCoords


###############
#Analyse speed#
###############

# Create a smoothed trajectory, filter order 3, length 31
smoothed <- TrajSmoothSG(trj, p=3, n=31)

#plot both
plot(trj, lwd = 1, lty = 1)
lines(smoothed, col = "#FF0000A0", lwd = 2)
legend("topright", c("Original", "Smoothed"), lwd = c(1, 2), lty = c(1, 1), col = c("black", "red"), inset = 0.01)

# Calculate speed and acceleration
derivs <- TrajDerivatives(smoothed)
derivs
# Plot change-in-speed and speed
plot(derivs$acceleration ~ derivs$accelerationTimes, type = 'l', col = 'red', 
     yaxt = 'n',
     xlab = 'Time (s)',
     ylab = expression(paste('Change in speed (', m/s^2, ')')))
axis(side = 2, col = "red")
lines(derivs$speed ~ derivs$speedTimes, col = 'blue')
axis(side = 4, col = "blue")
mtext('Speed (m/s)', side = 4, line = 3)
abline(h = 0, col = 'lightGrey')
#in blue è la speed ma manca nel pot, dovrebbe essere scritto alla destra

#Once the trajectory speed has been obtained, it is simple to calculate values such as mean, maximum, minimum or standard deviation of speed: 
mean(derivs$speed)
max(derivs$speed)
min(derivs$speed)
sd(derivs$speed)

# Calculate hovering intervals
#slowerThan = here is 10 m/s or cm/s?
intervals <- TrajSpeedIntervals(trj, slowerThan = 10)
print(intervals)
# Plot speed over time with hovering intervals highlighted
plot(intervals)

######################
#Analyse straightness#
######################
#straightness index: simplest is D/L, where D is the distance from the start to the end 
#of the trajectory, and L is the length of the trajectory
#function TrajStraightness is a number ranging from 0 to 1, where 1 indicates a straight line
#straightness index is approximation of r, which is the length of the mean vector of turning angles after rediscretizing to a constant step length. 
#r can be calculated by calling Mod(TrajMeanVectorOfTurningAngles(trj)), assuming trj is a Trajectory with constant step length
Mod(TrajMeanVectorOfTurningAngles(trj))

#TrajDistance Calculates the distance between the start and end of a trajectory (or a portion of a trajectory)
#would like to check the distance from frame 1 to "food"... all the "trip_to"
#add the colum food journey to x, that has the coordinates of the shrew
track_sf_2$food_journey -> x$food_journey 
track1 <- x %>% 
  select(x, y, time, food_journey) %>% 
  #subselect only the trip to
  filter(food_journey == 'trip_to')

trj1 <- TrajFromCoords(track1, fps = 30, spatialUnits = "cm")
TrajDistance(trj1, startIndex = 1, endIndex = nrow(trj1))
TrajStraightness(trj1)


#check orientation of traj_to, direction the vector is pointing at any time
#forloop, take current fram, previous frame and food position
#angle from food distance from food, 