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
# STEP 1: open data for all trials ------------------------------------------
setwd("~/data/cue_learning")
doors <- read.csv("trial_door.csv") %>%
  mutate(Trial = paste0("T",trial_n))

coords <- read.csv("food_door_coordinates.csv", sep = ",",header = TRUE) %>%
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
#transform inf values in NA, then drop them
tracking[c('ANGLE')][sapply(tracking[c('ANGLE')], is.infinite)] <- NA
tracking <- tracking %>% 
  drop_na(ANGLE)

#empty dataframe for saving info on visit to other doors on trip back
other_door_visits <- data.frame(ID = NULL, door = NULL, Trial = NULL, Season = NULL) #we can also add the length of time spent at the door if needed

trial_ls <- split(tracking, tracking$unique_trial_ID)

#prepare cluster for parallel computation
#mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)
mycl <- makeCluster(3)

clusterExport(mycl, c("coords", "trial_ls", "doors", "other_door_visits")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
})

b <- Sys.time()
other_door_visits_ls <- lapply(trial_ls, function(x){
  #to call all blocks one at a time
  #x = trial_ls[[15]]
  
  #extract food coordinates for this trial AND convert to a sf object
  food_coords <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("FOOD_x", "FOOD_y")) %>%
    st_as_sf(coords = c("FOOD_x", "FOOD_y"))
  
  food_buffer <- food_coords %>%
    st_buffer(dist = 2.5) #half of the length of the largest possible shrew
  
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
    st_buffer(dist = 3.5)
  
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
  if (nrow(at_food) > 0){
    track_sf_2 <- track_sf %>%
      full_join(at_food[c("frame", "visit_seq")]) %>%
      arrange(frame) %>%
      mutate(food_journey = ifelse(frame == head(at_food$frame,1), "arrival", #the first point in the at_food data is the arrival
                                   ifelse(frame == tail(at_food[at_food$visit_seq == 1, "frame"],1), "departure", #the last point in the first visit to food
                                          ifelse(frame < head(at_food$frame,1), "trip_to", #points before point of arrival are trip to food
                                                 ifelse(between(frame, head(at_food[at_food$visit_seq == 1, "frame"],1), tail(at_food[at_food$visit_seq == 1,"frame"],1)), "at_food", #time spent at food, between arrival and departure from food
                                                        ifelse(frame %in% at_food[at_food$visit_seq != 1,"frame"],
                                                               paste("trip_back_revisit", visit_seq, sep = "_"),
                                                               "trip_back")))))) %>% #all other points are trip back from food
      rowwise() %>%
      mutate(visit_to_prev_door = ifelse(Trial == "T1", NA, #assign NA for the first trial
                                         ifelse(nrow(other_door_visits[other_door_visits$Trial == Trial,]) == 0, "No", #if this trial isn't in the other_door_visits, assign no
                                                ifelse(lag(doors,1) %in% as.character(other_door_visits[other_door_visits$Trial == Trial, "door_ID"]), "Yes", "No"))))
      #I need the x,y in the .csv sheet. so I drop the geometry before saving the .csv file. but track_sf2 has to keep the geometry, so I create a new variable to save the file
      track_save <- track_sf_2 %>% 
        extract(geometry, c('x', 'y'), '\\((.*), (.*)\\)', convert = TRUE) %>% 
        relocate(x, .after = frame) %>% 
        relocate(y, .after = x) %>% 
        relocate(unique_trial_ID, .before = ID)
    
    write.csv(track_save, file = paste0("/home/ceci/data/cue_learning/results/", unique(x$unique_trial_ID),".csv"))
    
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
    df = data.frame(x[2,2],x[1,4],x[2,3] ,  dist_doorfood, walk_to, walk_back, straight_index_to, straight_index_back)
    colnames(df)[1] = "ID"
    colnames(df)[2] = "season"
    colnames(df)[3] = "trial"
    colnames(df)[4] = "food_door"
    colnames(df)[5] = "walked_to"
    colnames(df)[6] = "walked_back"
    colnames(df)[7] = "straightness_to_food"
    colnames(df)[8] = "straightness_back"
    
    write.csv(df, file = paste0("/home/ceci/data/cue_learning/distance/", unique(x$unique_trial_ID),".csv"))
    # check for overlap between the return trip and other doors
    
  }
  #just put this whole block out of the previous loop
  other_doors <- track_sf_2 %>%
    filter(food_journey == "trip_back") %>%
    st_intersection(all_doors_buffer %>% filter(door_ID != trial_door_ID))
  
  #if there was a visit to another door, save info to other_door_visits dataframe
  if(nrow(other_doors) > 0){
    new_visits <- other_doors %>%
      group_by(door_ID) %>%  slice(1) %>%  #this will give you one row per other door visited
      dplyr::select(c("ID", "Season", "Trial", "door_ID")) %>%
      st_drop_geometry() #convert back to non-spatial object
    
    #append to other_door_visits
    other_door_visits <<- rbind(other_door_visits,new_visits) #double arrow assignment operator allows to modify the dataframe in the global environment  
    #return other door visits
    return(new_visits)
  
#    for x in trial_ls
      
  } else {
    print(paste0("trial ", unique(x$unique_trial_ID), " not done"))
    #write.csv(data.frame(NULL), file = paste0("/home/ceci/data/cue_learning/results/", unique(x$unique_trial_ID),"_empty.csv"))
  }
  
  print(paste0("trial ", unique(x$unique_trial_ID), " completed."))
  
})

write.csv(other_door_visits, "/home/ceci/Documents/data/cue_learning/other_door_visit.csv", sep = ",", row.names = FALSE)

staSys.time() - b
stopCluster(mycl)

#mapview(trial_door_coords) + mapview(food_buffer) + mapview(food_coords) + mapview(trial_door_buffer) + mapview(track_sf)

# other_door_visits <- other_door_visits_ls %>%
#   reduce(rbind)

#saveRDS(sp_prep, file = "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_sp.rds")
saveRDS(other_door_visits_ls, file = "/home/ceci/Documents/data/cue_learning/trials_sp.rds")
#saveRDS(other_door_visits, "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_other_doors.rds")
saveRDS(other_door_visits, "/home/ceci/Documents/data/trials_other_doors.rds")

# STEP 4: Calculate speed, distance, tortuosity ------------------------------------------

sp_prep <- readRDS("/home/ceci/Documents/data/cue_learning/trials_sp.rds")
sp_prep$`spring_T1_20201031-1`
trial_summaries <- sp_prep %>%
  reduce(rbind) %>% #merge the list elements
  group_by(unique_trial_id, food_journey) %>%  #calculate length and duration for each leg of the food journey within each trial for each individual
  arrange(frame) %>%
  summarise(do_union = FALSE,
            duration = tail(time,1) - head(time,1)) %>%  #calculate the number of seconds spent in each part of the food journey (the difference between the last timestamp and the first)
  st_cast("LINESTRING") %>%  #convert to line object to calculate length.
  ungroup() %>%
  mutate(length = st_length(.)) %>%  #in cm. make sure it makes sense.
  rowwise() %>%
  mutate(tortuosity = ifelse(food_journey %in% c("trip_to","trip_back"),
                             length/as.numeric(doors_pt[doors_pt$Door == Door,"shortest_path"])[1], NA)) %>% #calculate tortuosity as the ratio between length of track and the shortest path (only for way to food and back)
  st_drop_geometry() #remove the spatial portion of the dataset

saveRDS(trial_summaries, "/home/enourani/ownCloud/Work/Collaborations/Cecilia_2022/5_trials_summaries.rds")


###########
###PLOTS###
###########

library(viridis)
install.packages("ggsci")
library(ggsci)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

#I need x, y and food_journey
#all three can be found in the .csv saved from the lapply call
#make track_all a list, with split(track_all, track_all$unique_trial_ID), then for loop plots
track_all <- read.csv("/home/ceci/data/cue_learning/results/all_track.csv", header = TRUE)

#there are some trials with more than 20 revisits, I want to filter them
df_revisits <- dplyr::filter(track_all, grepl('revisit_40', food_journey))
df_revisits['unique_trial_ID']
#List based on unique_trial_ID
track_all <- split(track_all, track_all$unique_trial_ID)
select.list(track_all, dplyr::starts_with("spring_T10_20201106-1"))
a <-track_all$`spring_T10_20201106-1`
plotfood <- coords %>%
  as.data.frame(plotfood, row.names = NULL) %>% 
  filter(unique_trial_ID == unique(a$unique_trial_ID)) %>%
  dplyr::select(c("FOOD_x", "FOOD_y", "unique_trial_ID"))
plotdoor <- coords %>%
  filter(unique_trial_ID == unique(a$unique_trial_ID)) %>%
  select(4:11, 16)
plot <- a %>% ggplot(aes(x, y, colour = food_journey)) +
  ggtitle(a$unique_trial_ID) +
  geom_point(x = plotdoor$A_x, y = plotdoor$A_y,  size = 3, colour = "black") +
  geom_point(x = plotdoor$B_x, y = plotdoor$B_y,  size = 3, colour = "black") +
  geom_point(x = plotdoor$C_x, y = plotdoor$C_y,  size = 3, colour = "black") +
  geom_point(x = plotdoor$D_x, y = plotdoor$D_y,  size = 3, colour = "black") +
  geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 8, colour = "darkgreen", alpha = 1/20) +
  geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 3, colour = "green") +
  geom_path()
print(plot)
#List based on ID
# track_all <- split(track_all, track_all$ID)
track_all$food_journey <- as.factor(track_all$food_journey) 
#OR 
##track_all[, 'food_journey'] <- as.factor(track_all[, 'food_journey'])

sapply(track_all, levels)

testa <- head(track_all)
lapply(testa, function(i){
  #i = track_all[[c("spring_T10_20201106-1")]]
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
    #  scale_color_grey(start = 0.8, end = 0.2) +
    #  scale_color_viridis(option = "D") +
    
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  print(plots)
  #  }
})


#########
#recurse#
#########
#recurse need a dataframe with x,y, timestamp and ID
#timestamp should be POSIXct format. but mine is just seconds

#check leo the vulture r script and try to adapt it!

####################
#Plots STRAIGHTNESS#
####################

distance <- read_csv("/home/ceci/data/cue_learning/distance/all.csv")

ggplot(data = distance) + 
  geom_point(mapping = aes(x = trial, y = walked_to)) +
  facet_grid(~season)

ggplot(data = distance) + 
  geom_point(mapping = aes(x = trial, y = walked_back, color = season))+
  geom_smooth(mapping = aes(x = trial, y = walked_back, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))
  

ggplot(data = distance) + 
  geom_smooth(mapping = aes(x = trial, y = walked_to, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

ggplot(data = distance) + 
  geom_point(mapping = aes(x = trial, y = food_door, color = season))+
  geom_smooth(mapping = aes(x = trial, y = walked_to, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))  

ggplot(data = distance) + 
  #geom_point(mapping = aes(x = trial, y = food_door, color = season))+
  geom_smooth(mapping = aes(x = trial, y = walked_back, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

ggplot(data = distance) + 
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

ggplot(data = distance) + 
  geom_point(mapping = aes(x = trial, y = straightness_back, color = season))+
  geom_smooth(mapping = aes(x = trial, y = straightness_back, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))
