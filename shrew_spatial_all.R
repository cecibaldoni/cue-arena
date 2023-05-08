#script for spatial analysis of shrew movement
#Elham Nourani, PhD. May 19.2022. Konstanz, Germany.
#update to analyse all trials in all seasons: July 4. 2022
#update Oct. 12. 2022
#--------------------------------------------------------------------------

library(tidyverse)
library(sf) #for spatial manipulation
library(mapview) #for interactive visualization of the spatial data
library(parallel)
library(ggplot2)
library(trajr)
library(brms)
library(gifski)
library(png)
library(transformr)
library(gganimate)

# STEP 1: open data for all trials ------------------------------------------
setwd("~/data/cue_learning")
doors <- read.csv("~/data/cue_learning/trial_door.csv") %>%
  mutate(Trial = paste0("T",trial_n))

coords <- read.csv("~/data/cue_learning/food_door_coordinates.csv", sep = ",",header = TRUE) %>%
  mutate(Trial = paste0("T",TRIAL)) %>% 
  mutate(unique_trial_ID = paste(SEASON, Trial, ID, sep = "_"))

# #read all csv in the folder, put them together, and create a unique_trial_id
# tracking <- lapply(list.files("/home/ceci/data/cue_learning/tracking", full.names = T), read.csv) %>%
#   reduce(rbind) %>%
#   drop_na(frame) %>%
#   mutate(unique_trial_ID = paste(Season, Trial, ID, sep = "_")) #create unique trial IDs. there are 10 trials per season per shrew
# #save this 'tracking' dataframe as a csv. from now on I can just load it
# write.csv(tracking, "/home/ceci/data/cue_learning/data.csv", row.names=FALSE)

tracking <- read.csv("/home/ceci/data/cue_learning/data.csv")
tracking[c('ANGLE')][sapply(tracking[c('ANGLE')], is.infinite)] <- NA #transform inf values in NA, then drop them
tracking <- tracking %>% 
  drop_na(ANGLE)

#empty dataframe for saving info on visit to other doors on trip back
other_door_visits <- data.frame(ID = NULL, door = NULL, Trial = NULL, Season = NULL) #we can also add the length of time spent at the door if needed
no_visits <- data.frame(unique_trial_ID = NULL, other_door_visits = NULL, ID = NULL, trial = NULL, season = NULL)
trial_ls <- split(tracking, tracking$unique_trial_ID)

#prepare cluster for parallel computation
#mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)
mycl <- makeCluster(4)

clusterExport(mycl, c("coords", "trial_ls", "doors", "other_door_visits")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
})

b <- Sys.time()
other_door_visits_ls <- lapply(trial_ls, function(x){
  #to call all blocks one at a time
  #x = trial_ls[[1]]
  
  #extract food coordinates for this trial AND convert to a sf object
  food_coords <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("FOOD_x", "FOOD_y")) %>%
    st_as_sf(coords = c("FOOD_x", "FOOD_y"))
  
  food_buffer <- food_coords %>%
    st_buffer(dist = 4) #half of the length of the largest possible shrew
  
  #extract door coordinates for this trial AND convert to a sf object
  trial_door_ID <- doors %>%
    filter(Trial == unique(x$Trial)) %>%
    pull(door)
  
  doors_x <-  coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    select(4:11) %>%
    pivot_longer(cols = contains("x"), names_to = "door", values_to = "x") %>%
    select(x)
  
  doors_coords <-  coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    select(4:11) %>%
    pivot_longer(cols = contains("y"), names_to = "door", values_to = "y") %>%
    mutate(door_ID = substr(door,1,1)) %>%
    select(c("door_ID", "y")) %>%
    bind_cols(doors_x)
  
  all_doors_buffer <- doors_coords %>%
    st_as_sf(coords = c("x","y")) %>%
    st_buffer(dist = 4)
  
  trial_door_buffer <- all_doors_buffer %>%
    filter(door_ID == trial_door_ID)
  
  #convert the track into an sf object
  track_sf <- x %>%
     st_as_sf(coords = c("x", "y"))
  #separated_coord <- track_sf %>%
  #  extract(geometry, c('x', 'y'), '\\((.*), (.*)\\)', convert = TRUE) # %>% 
    
  #find the points of overlap between track and food
  #WARNING MESSAGE
  at_food <- track_sf %>%
    st_intersection(food_buffer) %>%
    as.data.frame() %>%
    arrange(frame) %>% #arrange by time/frame
    mutate(timediff = frame - lag(frame)) %>%
    mutate(new_timediff = ifelse(is.na(timediff) | timediff != 1, 1,0)) %>%
    mutate(visit_seq = cumsum(new_timediff))

  #add food journey info to x IF food was reached
  #if (nrow(at_food) > 0){
  track_sf_2 <- track_sf %>%
    full_join(at_food[c("frame", "visit_seq")]) %>%
    arrange(frame) %>%
    mutate(old_food_journey = ifelse(frame == head(at_food$frame,1), "arrival", #the first point in the at_food data is the arrival
                                 ifelse(frame == tail(at_food[at_food$visit_seq == 1, "frame"],1), "departure", #the last point in the first visit to food
                                        ifelse(frame < head(at_food$frame,1), "trip_to", #points before point of arrival are trip to food
                                               ifelse(between(frame, head(at_food[at_food$visit_seq == 1, "frame"],1), tail(at_food[at_food$visit_seq == 1,"frame"],1)), "at_food", #time spent at food, between arrival and departure from food
                                                      ifelse(frame %in% at_food[at_food$visit_seq != 1,"frame"],
                                                             paste("trip_back_revisit", visit_seq, sep = "_"),
                                                              "trip_back")))))) #all other points are trip back from food
  #THIS DOES NOT WORK  
  #I need to make "being at the exit" sequences, and use this information to label "exploration"
  #if "being at the exit" only happens once, keep it as trip_back. from the frame after the first "being at the exit sequence" it's exploration.
   #find overlap exit and track
  at_exit <- track_sf_2 %>% 
    filter(old_food_journey == "trip_back") %>% 
    st_intersection(all_doors_buffer %>% filter(door_ID == trial_door_ID)) %>% # Find frames with intersection with the buffer of the door
    as.data.frame() %>% 
    mutate(timediff = frame - lag(frame)) %>%
    mutate(new_timediff = ifelse(is.na(timediff) | timediff != 1, 1,0)) %>%
    mutate(exit_seq = cumsum(new_timediff))
  
  # trip_back_indices <- which(track_sf_2$old_food_journey == "trip_back") #which values of food_journey are labeled "trip_back"
  # intersection_frames <- st_intersection(track_sf_2[trip_back_indices,], 
  #                                        all_doors_buffer %>% filter(door_ID == trial_door_ID)) # Find frames with intersection with the buffer of the door
  # first_intersection <- min(trip_back_indices[trip_back_indices %in% intersection_frames$frame])# Get the index of the first intersection
  # last_intersection <- max(trip_back_indices[trip_back_indices %in% intersection_frames$frame]) # Get the index of the last intersection
  

  if (nrow(at_exit) > 0) {
    track_sf_2 <- track_sf_2 %>% 
      full_join(at_exit[c("frame", "exit_seq")]) %>%
      arrange(frame) %>% 
      mutate(food_journey = ifelse(frame > tail(at_exit[at_exit$exit_seq == 1, "frame"], 1), "exploration", as.character(old_food_journey))) %>%
      mutate(food_journey = ifelse(food_journey == "exploration" & !is.na(visit_seq),paste0("exploration_revisit_", visit_seq), food_journey))

  } else {
    track_sf_2 <- track_sf_2 %>%
      mutate(food_journey = old_food_journey)
  }
    #I need the x,y in the .csv sheet. so I drop the geometry before saving the .csv file. but track_sf2 has to keep the geometry, so I create a new variable to save the file
  track_save <- track_sf_2 %>% 
    extract(geometry, c('x', 'y'), '\\((.*), (.*)\\)', convert = TRUE) %>% 
    relocate(x, .after = frame) %>% 
    relocate(y, .after = x) %>% 
    relocate(unique_trial_ID, .before = ID) #%>% 
    #select(-old_food_journey)
    
  write.csv(track_save, file = paste0("/home/ceci/data/cue_learning/results/", unique(x$unique_trial_ID),".csv"))
  
  if (sum(track_sf_2$food_journey == "exploration") > 0) {
    track_sf_2$food_journey -> x$food_journey 
    trip_to <- x %>% 
      select(x, y, time, food_journey) %>%
      filter(food_journey == "trip_to")
    trip_back <- x %>% 
      select(x, y, time, food_journey) %>%
      filter(food_journey == "trip_back")
    exploration <- x %>% 
      select(x, y, time, food_journey) %>%
      filter(food_journey == "exploration")
    trj1 <- TrajFromCoords(trip_to, fps = 30, spatialUnits = "cm")
    #if to use smoothed_to: change in TrajDistance(smoothed_to, startIndex = 1, endIndex = nrow(trj1))
    dist_doorfood <- TrajDistance(trj1, startIndex = 1, endIndex = nrow(trj1))
    #distance_to <- TrajDistance(trj1, startIndex = 1, endIndex = nrow(trj1))
    trj2 <- TrajFromCoords(trip_back, fps = 30, spatialUnits = "cm")
    #smoothed_back <- TrajSmoothSG(trj2, p=3, n=19)
    trj3 <- TrajFromCoords(exploration, fps = 30, spatialUnits = "cm")
    walk_to <- TrajLength(trj1, startIndex = 1, endIndex = nrow(trj1))
    walk_back <- TrajLength(trj2, startIndex = 1, endIndex = nrow(trj2))
    straight_index_to <- TrajStraightness(trj1)
    straight_index_back <- TrajStraightness(trj2)
    straight_exploration <- TrajStraightness(trj3)
    ############
    #Added visit to previous door, check if it works
    df = data.frame(x[1,1],x[1,3],x[1,2] ,  dist_doorfood, walk_to, walk_back, straight_index_to, straight_index_back, straight_exploration)
    colnames(df)[1] = "ID"
    colnames(df)[2] = "season"
    colnames(df)[3] = "trial"
    colnames(df)[4] = "food_door"
    colnames(df)[5] = "walked_to"
    colnames(df)[6] = "walked_back"
    colnames(df)[7] = "straightness_to_food"
    colnames(df)[8] = "straightness_back"
    colnames(df)[9] = "straightness_exploration"
    
    write.csv(df, file = paste0("/home/ceci/data/cue_learning/distance/", unique(x$unique_trial_ID),".csv"))

  } else {
    track_sf_2$food_journey -> x$food_journey 
    trip_to <- x %>% 
      select(x, y, time, food_journey) %>%
      filter(food_journey == "trip_to")
    trip_back <- x %>% 
      select(x, y, time, food_journey) %>%
      filter(food_journey == "trip_back")
    trj1 <- TrajFromCoords(trip_to, fps = 30, spatialUnits = "cm")
    #if to use smoothed_to: change in TrajDistance(smoothed_to, startIndex = 1, endIndex = nrow(trj1))
    dist_doorfood <- TrajDistance(trj1, startIndex = 1, endIndex = nrow(trj1))
    #distance_to <- TrajDistance(trj1, startIndex = 1, endIndex = nrow(trj1))
    trj2 <- TrajFromCoords(trip_back, fps = 30, spatialUnits = "cm")
    #smoothed_back <- TrajSmoothSG(trj2, p=3, n=19)
    walk_to <- TrajLength(trj1, startIndex = 1, endIndex = nrow(trj1))
    walk_back <- TrajLength(trj2, startIndex = 1, endIndex = nrow(trj2))
    straight_index_to <- TrajStraightness(trj1)
    straight_index_back <- TrajStraightness(trj2)
    ############
    #Added visit to previous door, check if it works
    df = data.frame(x[1,1],x[1,3],x[1,2] ,  dist_doorfood, walk_to, walk_back, straight_index_to, straight_index_back)
    colnames(df)[1] = "ID"
    colnames(df)[2] = "season"
    colnames(df)[3] = "trial"
    colnames(df)[4] = "food_door"
    colnames(df)[5] = "walked_to"
    colnames(df)[6] = "walked_back"
    colnames(df)[7] = "straightness_to_food"
    colnames(df)[8] = "straightness_back"
    
    write.csv(df, file = paste0("/home/ceci/data/cue_learning/distance/", unique(x$unique_trial_ID),".csv"))
  }
  
    # check for overlap between the return trip and other doors
  #}
  #just put this whole block out of the previous loop
  other_doors <- track_sf_2 %>%
    filter(food_journey %in% c("trip_back", "exploration")) %>%
    st_intersection(all_doors_buffer %>% filter(door_ID != trial_door_ID))
  
  #if there was a visit to another door, save info to other_door_visits dataframe
  if(nrow(other_doors) > 0){
    new_visits <- other_doors %>%
      group_by(door_ID) %>%  slice(1) %>%  #this will give you one row per other door visited
      dplyr::select(c("ID", "Season", "Trial", "door_ID", "food_journey")) %>%
      st_drop_geometry() #convert back to non-spatial object
    
    #append to other_door_visits
    other_door_visits <<- rbind(other_door_visits,new_visits) #double arrow assignment operator allows to modify the dataframe in the global environment  
    #return other door visits
    return(new_visits)
      
  } else {
    
    empty_data <- data.frame(unique_trial_ID = unique(x$unique_trial_ID), other_door_visits = 0) %>% 
      separate(unique_trial_ID, into = c("season", "trial", "ID"), sep = "_")
    no_visits <<- rbind(no_visits, empty_data)
   
  }
  
  print(paste0("trial ", unique(x$unique_trial_ID), " completed."))
  
})
write.csv(no_visits, "/home/ceci/data/cue_learning/results/no_visits.csv", row.names = FALSE)
write.csv(other_door_visits, "/home/ceci/data/cue_learning/other_door_visit.csv", row.names = FALSE)

staSys.time() - b
stopCluster(mycl)

# other_door_visits <- other_door_visits_ls %>%
#   reduce(rbind)

#saveRDS(sp_prep, file = "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_sp.rds")
saveRDS(other_door_visits_ls, file = "/home/ceci/Documents/data/cue_learning/trials_sp.rds")
#saveRDS(other_door_visits, "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_other_doors.rds")
saveRDS(other_door_visits, "/home/ceci/Documents/data/trials_other_doors.rds")

#################
##TURNING ANGLE##
#################

# Loop through each dataframe in the list
for (i in 1:length(trial_ls)) {
  # Extract x and y columns from the dataframe
  x <- trial_ls[[i]]$x
  y <- trial_ls[[i]]$y
  # Calculate direction and turning angle
  delta_x <- diff(x)
  delta_y <- diff(y)
  direction <- atan2(delta_y, delta_x)
  turning_angle <- c(0, diff(direction))
  turning_angle <- c(turning_angle, NA) # Add missing value to match number of rows
  turning_angle_degrees <- turning_angle * 180 / pi
  # Add turning angle as a new column to the dataframe
  trial_ls[[i]]$turning_angle <- turning_angle_degrees
}

# i should define threshold value for turning angle
angle_threshold <- 125 #???

# Loop through each dataframe in the list
for (i in 1:length(trial_ls)) {
  # Extract turning angle column from the dataframe
  turning_angle <- trial_ls[[i]]$turning_angle
  # Find indices where turning angle is above the threshold
  abrupt_turns <- which(abs(turning_angle) > angle_threshold)
  # Add abrupt turns as a new column to the dataframe
  trial_ls[[i]]$abrupt_turns <- FALSE
  trial_ls[[i]]$abrupt_turns[abrupt_turns] <- TRUE
}


for (i in 1:length(trial_ls)) {
  # Extract x, y, and abrupt_turns columns from the dataframe
  #i <- trial_ls[[i]]
  x <- trial_ls[[i]]$x
  y <- trial_ls[[i]]$y
  abrupt_turns <- trial_ls[[i]]$abrupt_turns
  
  # Create a data frame for the plot
  plot_data <- data.frame(x = x, y = y, abrupt_turns = abrupt_turns)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = abrupt_turns)) + 
    geom_point() +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = paste("Path with Abrupt Turns (Dataframe", i, ")"))
  
  print(p)
}


#########
##SPEED##
#########

for (i in 1:length(trial_ls)) {
  x <- trial_ls[[i]]$x
  y <- trial_ls[[i]]$y
  time <- trial_ls[[i]]$time
  
  #calculate speed based on x, y, and time
  x_diff <- c(0, diff(x))
  y_diff <- c(0, diff(y))
  time_diff <- c(0, diff(time))
  distance <- sqrt(x_diff^2 + y_diff^2)
  speed <- distance / time_diff
  
  #add speed as a new column
  trial_ls[[i]]$speed_new <- speed
  
  # wxpected speed based on frame difference
  expected_speed <- distance / time_diff
  #threshold as the 90th percentile of the actual speed
  #it returns the value below which 90% of the data fall
  speed_threshold <- quantile(speed, 0.80, na.rm = TRUE)
  
  # Find indices where actual speed is above threshold multiple of expected speed
  high_speed <- which(speed > speed_threshold * expected_speed[-1])
  
  # Add high speed as a new column to the dataframe
  trial_ls[[i]]$high_speed <- FALSE
  trial_ls[[i]]$high_speed[high_speed] <- TRUE
}

for (i in 1:length(trial_ls)) {
  # Extract x, y, highspeed columns from the dataframe
  x <- trial_ls[[i]]$x
  y <- trial_ls[[i]]$y
  high_speed <- trial_ls[[i]]$high_speed
  
  # Create a data frame for the plot
  plot_data <- data.frame(x = x, y = y, high_speed = high_speed)
  
  # Create the plot
  s <- ggplot(plot_data, aes(x = x, y = y, color = high_speed)) + 
    geom_point() +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = paste("Path with High Speed (Dataframe", i, ")"))
  
  print(s)
}

################
##Previous door#
################

other_door_visits <- read.csv("/home/ceci/data/cue_learning/other_door_visit.csv", header = TRUE)
previous_door <- other_door_visits %>%
  mutate(unique_trial_ID = paste(Season, Trial, ID, sep = "_")) %>% 
  mutate(unique_trial_ID = as.factor(unique_trial_ID)) %>% 
  mutate(trial_n = str_remove(Trial, "^T")) %>% 
  mutate(trial_n = as.integer(trial_n))

merged_data <- inner_join(doors, previous_door, by = "trial_n", multiple = "all") %>% 
  select(-Trial.x, -Trial.y, -ID, -Season)

result <- merged_data %>%
  group_by(trial_n) %>%
  filter(trial_n != 1) %>% 
  mutate(previous_door_ID = (trial_n == 2 & door_ID == "A") |
           (trial_n == 3 & door_ID == "C") | 
           (trial_n == 4 & door_ID == "A") |
           (trial_n == 5 & door_ID == "B") |
           (trial_n == 6 & door_ID == "D") |
           (trial_n == 7 & door_ID == "C") |
           (trial_n == 8 & door_ID == "B") |
           (trial_n == 9 & door_ID == "D") |
           (trial_n == 10 & door_ID == "B"))

result <- result %>% 
  group_by(unique_trial_ID) %>%
  #filter rows in which uniquetrialID appears more than once
  filter(n() == 1 | n() > 1 & !duplicated(unique_trial_ID)) %>%
  ungroup() %>% 
  #separate the column uniquetrial id in 3 columns again
  separate(unique_trial_ID, into = c("season", "trial", "ID"), sep = "_") %>% 
  #get out columns I don't want
  select(-trial, -door, -door_ID) 

result <- result[order(result$trial_n), ]
result$trial <- result$trial_n
library(scales)
#relevel the factors, t10 comes always after t1
#result$trial <- factor(result$trial, levels = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"), ordered = TRUE)
result$season <- factor(result$season, levels = c("summer", "winter", "spring"))
ggplot(result, aes(x = trial, fill = factor(previous_door_ID))) +
  geom_bar() +
  scale_x_continuous(breaks= pretty_breaks()) +
  facet_grid(food_journey ~ season) +
  scale_fill_manual(values = c("red", "darkgreen")) +
  xlab("Trial") +
  ylab("Count") +
  ggtitle("Visit to previous door, by Season and Trial")

result$season <- factor(result$season, levels = c("summer", "winter", "spring"))
ggplot(result, aes(x = trial, fill = factor(previous_door_ID))) +
  geom_bar() +
  facet_grid(ID ~ season) +
  scale_fill_manual(values = c("red", "darkgreen")) +
  xlab("Trial") +
  ylab("Count") +
  ggtitle("Visit to previous door, by Season and Trial")

###########
###PLOTS###
###########

library(viridis)
library(ggsci)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

track_all <- read.csv("/home/ceci/data/cue_learning/results/all_track.csv", header = TRUE) %>% 
  select(-visit_seq, -exit_seq, -X)

#there are some trials with more than 20 revisits, I want to filter them
df_revisits <- dplyr::filter(track_all, grepl('revisit_', food_journey))
df_revisits['unique_trial_ID']

#List based on unique_trial_ID
track_all <- split(track_all, track_all$unique_trial_ID)
#select.list(track_all, dplyr::starts_with("spring_T10_20201106-1"))

b <-track_all$`spring_T3_20201031-1`
b <- track_all$`spring_T1_20201031-1`

#b$food_journey <- factor(b$food_journey, levels = c("trip_to", "arrival", "at_food", "departure", "trip_back", "exploration"), ordered = TRUE)
plotfood <- coords %>%
  as.data.frame(plotfood, row.names = NULL) %>% 
  filter(unique_trial_ID == unique(b$unique_trial_ID)) %>%
  dplyr::select(c("FOOD_x", "FOOD_y", "unique_trial_ID"))
plotdoor <- coords %>%
  filter(unique_trial_ID == unique(b$unique_trial_ID)) %>%
  select(4:11, 16)
plot <- b %>% ggplot(aes(x, y, colour = food_journey)) +
  ggtitle(b$unique_trial_ID) +
  geom_point(x = plotdoor$A_x, y = plotdoor$A_y,  size = 3, colour = "black") +
  geom_point(x = plotdoor$B_x, y = plotdoor$B_y,  size = 3, colour = "black") +
  geom_point(x = plotdoor$C_x, y = plotdoor$C_y,  size = 3, colour = "black") +
  geom_point(x = plotdoor$D_x, y = plotdoor$D_y,  size = 3, colour = "black") +
  geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 8, colour = "darkgreen", alpha = 1/20) +
  geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 3, colour = "green") +
  geom_point(size = 5) + transition_time(frame) +
  labs(title = "Frame: {frame_time}") +
  shadow_mark(alpha = 0.8, size = 2)
anim_save("trial.gif")
plot


#List based on ID
# track_all <- split(track_all, track_all$ID)
track_all$food_journey <- as.factor(track_all$food_journey) 
#OR 
##track_all[, 'food_journey'] <- as.factor(track_all[, 'food_journey'])

sapply(track_all, levels)

testa <- head(track_all)
lapply(testa, function(i){
  #i = track_all[[c("spring_T1_20201131-1")]]
  #i = testa[[1]]
  plotfood <- coords %>%
    as.data.frame(plotfood, row.names = NULL) %>% 
    filter(unique_trial_ID == unique(i$unique_trial_ID)) %>%
    dplyr::select(c("FOOD_x", "FOOD_y", "unique_trial_ID"))
  plotdoor <- coords %>%
    filter(unique_trial_ID == unique(i$unique_trial_ID)) %>%
    select(4:11, 16)
  
  plots <- i %>%
    #filter(food_journey == c("trip_to", "at_food", "trip_back")) %>% 
    ggplot(aes(x, y, colour = food_journey)) +
    ggtitle(i$unique_trial_ID) +
    geom_point(x = plotdoor$A_x, y = plotdoor$A_y,  size = 3, colour = "black") +
    geom_point(x = plotdoor$B_x, y = plotdoor$B_y,  size = 3, colour = "black") +
    geom_point(x = plotdoor$C_x, y = plotdoor$C_y,  size = 3, colour = "black") +
    geom_point(x = plotdoor$D_x, y = plotdoor$D_y,  size = 3, colour = "black") +
    geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 8, colour = "darkgreen", alpha = 1/20) +
    geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 3, colour = "green") +
    geom_path() #+ 
    #scale_colour_distiller(palette = "Reds") +
    #scale_color_grey(start = 0.8, end = 0.2) +
    #scale_color_viridis(option = "D") +
    
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  print(plots)
  #  }
})

path_data <- read.csv("~/data/cue_learning/distance/distance.csv") %>% 
  mutate(trial = str_remove(trial, "^T"))
path_data$trial <- str_sort(path_data$trial, numeric = TRUE)
path_data$trial <- as.integer(path_data$trial)

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = food_door, color = season))+
  geom_smooth(mapping = aes(x = as.numeric(trial), y = walked_to, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  #geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  geom_point(mapping = aes(x = as.numeric(trial), y = straightness_to_food, color = season)) +
  scale_color_manual(values = c("green", "coral", "cornflowerblue")) + 
  facet_wrap(~ID)

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  #geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  geom_point(mapping = aes(x = as.numeric(trial), y = straightness_exploration, color = season)) +
  scale_color_manual(values = c("green", "coral", "cornflowerblue")) + 
  facet_wrap(~ID)

ggplot(data = path_data) + 
  geom_smooth(mapping = aes(x = as.numeric(trial), y = straightness_to_food, color = season))+
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season)) +
  scale_color_manual(values = c("green", "coral", "cornflowerblue")) #+ 
  #facet_wrap(~season)

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_back, color = season))+
  geom_smooth(mapping = aes(x = as.numeric(trial), y = straightness_back, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))
ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_back, color = season))+
  geom_smooth(mapping = aes(x = as.numeric(trial), y = straightness_exploration, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

#########
#recurse#
#########
#recurse need a dataframe with x,y, timestamp and ID
#timestamp should be POSIXct format. but mine is just seconds

#check leo the vulture r script and try to adapt it!