#script for spatial analysis of shrew movement
#Elham Nourani, PhD. May 19.2022. Konstanz, Germany.
#update to analyse all trials in all seasons: July 4. 2022
#update Oct. 12. 2022
# Libraries-------------------------------------------------------------------

library(tidyverse)
library(sf) #for spatial manipulation
library(mapview) #for interactive visualization of the spatial data
library(parallel)
library(ggplot2)
library(trajr)
library(brms)
library(png)
library(transformr)
library(gganimate)
library(RColorBrewer)
# display.brewer.pal(n = 9, name = "Set1")
library(scales)

# STEP 1: open data for all trials ------------------------------------------

# coords <- read.csv("~/data/cue_learning/food_door_coordinates.csv", sep = ",",header = TRUE) %>%
#   #mutate(Trial = paste0("T",TRIAL)) %>%
#   mutate(ID = as.factor(ID),
#          SEASON = as.factor(SEASON),
#          TRIAL = as.factor(TRIAL)) %>%
#   #filter(SEASON != "spring") %>%
#   droplevels()
# str(coords)
# str(coords$SEASON)
# #THE FUCK
# coords$SEASON <- fct_collapse(coords$SEASON,
#                               summer = "summer",
#                               winter = c("winter", "winter "))
# coords <- coords %>%
#   mutate(unique_trial_ID = as.factor(paste(SEASON, ID, TRIAL, sep = "_")))
# tracking <- read.csv("~/data/cue_learning/csv/all_track.csv") %>%
#   mutate(status = case_when(season == "summer" ~ "summer_wild",
#                             season == "winter" & grepl("^2021", ID) ~ "winter_wild",
#                             season == "winter" & grepl("^2020", ID) ~ "winter_captive",
#                             season == "winter" & grepl("^2024", ID) ~ "winter_wild",
#                             season == "spring" & grepl("^2024", ID) ~ "spring_wild",
#                             season == "spring" & grepl("^2020", ID) ~ "spring_captive"),
#          unique_trial_ID = paste(season, ID, trial, sep = "_")) %>%
#   select(-time) %>%
#   relocate(status, .before = frame)
# tracking <- tracking %>%
#   mutate(unique_trial_ID = str_trim(unique_trial_ID, side = "right"),
#          trial = str_trim(trial, side = "right"))
# write.csv(tracking, "~/data/cue_learning/tracking.csv", row.names = FALSE)
# write.csv(coords, "~/data/cue_learning/coords.csv", row.names = FALSE)

doors <- read.csv("~/data/cue_learning/trial_door.csv")  %>%
  mutate(Trial = paste0("T",trial_n))
coords <- read.csv("~/data/cue_learning/coords.csv", sep = ",", header = TRUE)

tracking <- read.csv("~/data/cue_learning/tracking.csv")
tracking[c('x')][sapply(tracking[c('x')], is.infinite)] <- NA #transform inf values in NA, then drop them
tracking <- tracking %>% 
  drop_na(x)
#empty dataframe for saving info on visit to other doors on trip back
other_door_visits <- data.frame(ID = NULL, door = NULL, Trial = NULL, Season = NULL) #we can also add the length of time spent at the door if needed
no_visits <- data.frame(unique_trial_ID = NULL, other_door_visits = NULL, ID = NULL, trial = NULL, season = NULL)

trial_ls <- split(tracking, tracking$unique_trial_ID)
trial_ls <- lapply(trial_ls, function(df) {
  df <- df %>%
    mutate(frame = seq(0, nrow(df) - 1)) %>% 
    # Create the time column based on the frame column
    mutate(time = frame / 30) %>% 
    
    # Move the time column before the frame column
    relocate(time, .before = frame)
  
  return(df)  # Return the modified dataframe
})

#prepare cluster for parallel computation
#mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)
mycl <- makeCluster(4)
clusterExport(mycl, c("coords", "trial_ls", "doors", "other_door_visits")) #define the variable that will be used within the ParLapply call
clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
})

library(plotly)
#food out of buffer:
#spring_20201103-5_T2
#spring_20240526-1_T3 this has a weird space after 3
#summer_20200731-2_T2
y = trial_ls[["winter_20201031-1_T8"]]
food_coords <- coords %>%
  filter(unique_trial_ID == unique(y$unique_trial_ID))
doors_x <-  coords %>%
  filter(unique_trial_ID == unique(y$unique_trial_ID)) %>%
  select(4:11) %>%
  pivot_longer(cols = contains("x"), names_to = "door", values_to = "x") %>%
  select(x)
doors_coords <-  coords %>%
  filter(unique_trial_ID == unique(y$unique_trial_ID)) %>%
  select(4:11) %>%
  pivot_longer(cols = contains("y"), names_to = "door", values_to = "y") %>%
  mutate(door_ID = substr(door,1,1)) %>%
  select(c("door_ID", "y")) %>%
  bind_cols(doors_x)
p <- ggplot() +
  geom_point(data = y, aes(x = x, y = y, text = paste("X:", x, "Y:", y))) +  
  geom_point(data = doors_coords, aes(x=x, y=y, text = paste("Door", door_ID)), colour = "green", size = 5) +
  geom_point(data = food_coords, aes(x = FOOD_x, y = FOOD_y, text = paste("FOOD_X:", FOOD_x, "FOOD_Y:", FOOD_y)), colour = "red", size = 5) +
  ggtitle(paste("Trial ID:", y$unique_trial_ID[1]))
ggplotly(p, tooltip = "text")

## START LOOP ####

results <- "~/data/cue_learning/results/master_results.csv"
distance <- "~/data/cue_learning/distance/master_distance.csv"
other_door_visits_ls <- lapply(trial_ls, function(x){
  #to call all blocks one at a time
  #x = trial_ls[["spring_20240520-4_T8"]]
  
  #extract food coordinates for this trial AND convert to a sf object
  food_coords <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("FOOD_x", "FOOD_y")) %>%
    st_as_sf(coords = c("FOOD_x", "FOOD_y"))
  food_buffer <- food_coords %>%
    st_buffer(dist = 4) #half of the length of the largest possible shrew
  season_filter <- if (any(str_starts(x$unique_trial_ID, "winter_2024"))) {
    "winter_2024"
  } else {
    "other"
  }
  #extract door coordinates for this trial AND convert to a sf object
  trial_door_ID <- doors %>%
    filter(Trial == unique(x$trial)) %>%
    filter(season == season_filter) %>%
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
  #find the points of overlap between track and food
  #WARNING MESSAGE
  at_food <- track_sf %>%
    st_intersection(food_buffer) %>%
    as.data.frame() %>%
    arrange(frame) %>% #arrange by time/frame
    mutate(timediff = frame - lag(frame)) %>%
    #mutate(new_timediff = ifelse(is.na(timediff) | timediff != 1, 1,0)) %>%
    mutate(new_timediff = ifelse(is.na(timediff) | timediff >= 30, 1, 0)) %>%
    mutate(visit_seq = cumsum(new_timediff))
  #add food journey info to x IF food was reached
  
  if (nrow(at_food) == 0) {
    print(paste0("Food out of buffer zone in ", unique(x$unique_trial_ID)))
    return(NULL)
  }
  
  track_sf_2 <- track_sf %>%
    full_join(at_food[c("frame", "visit_seq")]) %>%
    arrange(frame) %>%
    mutate(old_food_journey = ifelse(frame == head(at_food$frame,1), "arrival", #the first point in the at_food data is the arrival
                                     ifelse(frame == tail(at_food[at_food$visit_seq == 1, "frame"],1), "departure", #the last point in the first visit to food
                                            ifelse(frame < head(at_food$frame,1), "trip_to", #points before point of arrival are trip to food
                                                   ifelse(between(frame, head(at_food[at_food$visit_seq == 1, "frame"],1), tail(at_food[at_food$visit_seq == 1,"frame"],1)), "at_food", #time spent at food, between arrival and departure from food
                                                          ifelse(frame %in% at_food[at_food$visit_seq != 1, "frame"], 
                                                                 paste("trip_back_revisit", visit_seq, sep = "_"), "trip_back")))))) #all other points are trip back from food
  
  #"being at the exit" sequences, and use this information to label "exploration"
  #if "being at the exit" only happens once, keep it as trip_back. from the frame after the first "being at the exit" it's exploration.
  #find overlap exit and track
  at_exit <- track_sf_2 %>% 
    filter(old_food_journey == "trip_back") %>% 
    st_intersection(all_doors_buffer %>% filter(door_ID == trial_door_ID)) %>% #herefind frames with intersection with the exit buffer
    as.data.frame() %>% 
    mutate(timediff = frame - lag(frame)) %>%
    mutate(new_timediff = ifelse(is.na(timediff) | timediff != 1, 1,0)) %>%
    mutate(exit_seq = cumsum(new_timediff))
  
  if (nrow(at_exit) > 0) { #if at_exit exists
    track_sf_2 <- track_sf_2 %>% 
      full_join(at_exit[c("frame", "exit_seq")]) %>%
      arrange(frame) %>% 
      # Label exploration visits with visit_seq similar to trip_back_revisit
      mutate(food_journey = ifelse(frame > tail(at_exit[at_exit$exit_seq == 1, "frame"], 1), 
                                   "exploration", 
                                   as.character(old_food_journey))) %>%
      mutate(food_journey = ifelse(food_journey == "exploration" & !is.na(visit_seq), 
                                   paste0("exploration_revisit_", visit_seq), 
                                   food_journey)) %>%
      select(-exit_seq)
    
  } else { #if at_exit does not exist, keep the old column
    track_sf_2 <- track_sf_2 %>%
      mutate(food_journey = old_food_journey)
  }
  #I need the x,y in the .csv sheet. so I drop the geometry before saving the .csv file. but track_sf2 has to keep the geometry, so I create a new variable to save the file
  track_save <- track_sf_2 %>%
    mutate(coordinates = st_coordinates(geometry)) %>%
    mutate(x = coordinates[, "X"], y = coordinates[, "Y"]) %>%
    relocate(x, .after = frame) %>% 
    relocate(y, .after = x) %>% 
    relocate(unique_trial_ID, .before = season) %>% 
    select(-old_food_journey, -visit_seq) %>% 
    st_drop_geometry()
  
#  write.csv(track_save, file = paste0("~/data/cue_learning/results/", unique(x$unique_trial_ID),".csv"), row.names = FALSE)
  write.table(track_save, file = results, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(results))
  
  if (sum(track_sf_2$food_journey == "exploration") > 0) { #some problems here, exploration is sometimes negative
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
    at_food <- x %>%
      filter(food_journey == "at_food") %>%
      select(time)
    
    if(nrow(at_food) > 1) {
      time_at_food <- max(at_food$time) - min(at_food$time)
    } else if(nrow(at_food) == 1) {
      time_at_food <- 1 / 30  # Assuming one frame represents 1/30th of a second
    } else {
      time_at_food <- 0  # No time spent at food
    }
    trj1 <- TrajFromCoords(trip_to, fps = 30, spatialUnits = "cm")
    #if to use smoothed_to: change in TrajDistance(smoothed_to, startIndex = 1, endIndex = nrow(trj1))
    dist_doorfood <- TrajDistance(trj1, startIndex = 1, endIndex = nrow(trj1))
    trj2 <- TrajFromCoords(trip_back, fps = 30, spatialUnits = "cm")
    #smoothed_back <- TrajSmoothSG(trj2, p=3, n=19)
    trj3 <- TrajFromCoords(exploration, fps = 30, spatialUnits = "cm")
    walk_to <- TrajLength(trj1, startIndex = 1, endIndex = nrow(trj1))
    walk_back <- TrajLength(trj2, startIndex = 1, endIndex = nrow(trj2))
    walk_expl <- TrajLength(trj3, startIndex = 1, endIndex = nrow(trj3))
    straight_index_to <- TrajStraightness(trj1)
    straight_index_back <- TrajStraightness(trj2)
    straight_exploration <- TrajStraightness(trj3)
    
    smooth_straightness_to <- NA  #initialize empty rows
    smooth_straightness_back <- NA
    
    adjust_to_odd <- function(n) {
      if (n %% 2 == 0) n - 1 else n
    }
    
    if(nrow(trip_to) > 20) {
      n_value <- adjust_to_odd(20)
      smoothed_to <- TrajSmoothSG(trj1, p = 3, n = n_value)
      smoothed_back <- TrajSmoothSG(trj2, p = 3, n = n_value)
      smooth_straightness_to <- TrajStraightness(smoothed_to)
      smooth_straightness_back <- TrajStraightness(smoothed_back)
    } else if (11 < nrow(trip_to) && nrow(trip_to) <= 20) {
      n_value <- adjust_to_odd(10)
      smoothed_to <- TrajSmoothSG(trj1, p = 3, n = n_value)
      smoothed_back <- TrajSmoothSG(trj2, p = 3, n = n_value)
      smooth_straightness_to <- TrajStraightness(smoothed_to)
      smooth_straightness_back <- TrajStraightness(smoothed_back)
    } else if (1 < nrow(trip_to) && nrow(trip_to) < 10){
      smooth_straightness_to <- straight_index_to
      smooth_straightness_back <- straight_index_back
    }
    
    # df = data.frame(x[1,1],x[1,2],x[1,3], x[1,4],
    #                 dist_doorfood, walk_to, walk_back, walk_expl, 
    #                 straight_index_to, straight_index_back, straight_exploration,
    #                 smooth_straightness_to, smooth_straightness_back,
    #                 time_at_food)
    # colnames(df) <- c(
    #   "season", "ID", "trial", "status",
    #   "food_door", "walked_to", "walked_back", "walk_exploration",
    #   "straightness_to_food", "straightness_back", "straightness_exploration",
    #   "smooth_straightness_to", "smooth_straightness_back",
    #   "time_at_food")
    df <- data.frame(
      season = as.character(x[1, 1]),
      ID = as.character(x[1, 2]),
      trial = as.character(x[1, 3]),
      status = as.character(x[1, 4]),
      dist_doorfood, walk_to, walk_back, walk_expl, 
      straight_index_to, straight_index_back, straight_exploration,
      smooth_straightness_to, smooth_straightness_back,
      time_at_food,
      stringsAsFactors = FALSE
    )
    colnames(df) <- c(
      "season", "ID", "trial", "status",
      "food_door", "walked_to", "walked_back", "walk_exploration",
      "straightness_to_food", "straightness_back", "straightness_exploration",
      "smooth_straightness_to", "smooth_straightness_back",
      "time_at_food"
    )
    
   # write.csv(df, file = paste0("~/data/cue_learning/distance/", unique(x$unique_trial_ID),".csv"), row.names = FALSE) #CHANGE DIRECTORY
    write.table(df, file = distance, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(distance))
    
  } else {
    track_sf_2$food_journey -> x$food_journey
    trip_to <- x %>%
      select(x, y, time, food_journey) %>%
      filter(food_journey == "trip_to")
    trip_back <- x %>%
      select(x, y, time, food_journey) %>%
      filter(food_journey == "trip_back")
    at_food <- x %>%
      filter(food_journey == "at_food") %>%
      select(time)
    
    if(nrow(at_food) > 1) {
      time_at_food <- max(at_food$time) - min(at_food$time)
    } else if(nrow(at_food) == 1) {
      time_at_food <- 1 / 30
    } else {
      time_at_food <- 0
    }
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
    walk_expl <- NA
    straight_exploration <- NA
    
    smooth_straightness_to <- NA  #initialize empty rows
    smooth_straightness_back <- NA
    
    adjust_to_odd <- function(n) {
      if (n %% 2 == 0) n - 1 else n
    }
    
    if(nrow(trip_to) > 20) {
      n_value <- adjust_to_odd(20)
      smoothed_to <- TrajSmoothSG(trj1, p = 3, n = n_value)
      smoothed_back <- TrajSmoothSG(trj2, p = 3, n = n_value)
      smooth_straightness_to <- TrajStraightness(smoothed_to)
      smooth_straightness_back <- TrajStraightness(smoothed_back)
    } else if (11 < nrow(trip_to) && nrow(trip_to) <= 20) {
      n_value <- adjust_to_odd(10)
      smoothed_to <- TrajSmoothSG(trj1, p = 3, n = n_value)
      smoothed_back <- TrajSmoothSG(trj2, p = 3, n = n_value)
      smooth_straightness_to <- TrajStraightness(smoothed_to)
      smooth_straightness_back <- TrajStraightness(smoothed_back)
    } else if (1 < nrow(trip_to) && nrow(trip_to) < 10){
      smooth_straightness_to <- straight_index_to
      smooth_straightness_back <- straight_index_back
    }
    
    df <- data.frame(
      season = as.character(x[1, 1]),
      ID = as.character(x[1, 2]),
      trial = as.character(x[1, 3]),
      status = as.character(x[1, 4]),
      dist_doorfood, walk_to, walk_back, walk_expl, 
      straight_index_to, straight_index_back, straight_exploration,
      smooth_straightness_to, smooth_straightness_back,
      time_at_food,
      stringsAsFactors = FALSE
    )
    colnames(df) <- c(
      "season", "ID", "trial", "status",
      "food_door", "walked_to", "walked_back", "walk_exploration",
      "straightness_to_food", "straightness_back", "straightness_exploration",
      "smooth_straightness_to", "smooth_straightness_back",
      "time_at_food"
    )
    
    #write.csv(df, file = paste0("~/data/cue_learning/distance/", unique(x$unique_trial_ID),".csv"), row.names = FALSE) #CHANGE DIRECTORY
    write.table(df, file = distance, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(distance))
  }
  #check for overlap between the return trip and other doors
  other_doors <- track_sf_2 %>%
    filter(food_journey %in% c("trip_back", "exploration")) %>%
    st_intersection(all_doors_buffer %>% filter(door_ID != trial_door_ID))
  
  #if there was a visit to another door, save info to other_door_visits dataframe
  
  if(nrow(other_doors) > 0){
    new_visits <- other_doors %>%
      group_by(door_ID) %>%  slice(1) %>%  #this will give you one row per other door visited
      dplyr::select(c("ID", "season", "trial", "door_ID", "food_journey")) %>%
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
  
#write.csv(no_visits, "/home/ceci/data/cue_learning/no_visits.csv", row.names = FALSE)
write.csv(other_door_visits, "~/data/cue_learning/other_door_visitcsv", row.names = FALSE)

stopCluster(mycl)

# saveRDS(other_door_visits_ls, file = "/home/ceci/Documents/data/cue_learning/trials_sp.rds")
# saveRDS(other_door_visits, "/home/ceci/Documents/data/trials_other_doors.rds")

tracking_results <- read.csv("~/data/cue_learning/results/master_results.csv") %>%
  mutate(season = as.factor(season),
         ID = as.factor(ID),
         status = as.factor(status),
         food_journey = as.factor(food_journey),
         trial = as.factor(trial),
         unique_trial_ID =as.factor(unique_trial_ID)) %>% 
  #filter(!is.na(season) & season != "1") %>% 
  droplevels()
str(tracking_results)
levels(tracking_results$food_journey)
tracking_results %>%
  filter(food_journey == "exploration_revisit_2") %>%
  select(unique_trial_ID) %>%
  distinct()
which(tracking_results$food_journey == "")
list_new <- split(tracking_results, tracking_results$unique_trial_ID)


#SPEED###########

for (i in 1:length(list_new)) {
  data <- list_new[[i]]
  
  x_diff <- c(NA, diff(data$x))
  y_diff <- c(NA, diff(data$y))
  time_diff <- c(NA, diff(data$time))

  distance <- sqrt(x_diff^2 + y_diff^2)
  speed <- ifelse(time_diff > 0, distance / time_diff, NA)

  data$speed_new <- speed
  
  list_new[[i]] <- data %>%
    filter(!(food_journey %in% c("arrival", "departure")))
}
tracked_speed <- do.call(rbind, list_new)
max(tracked_speed$speed_new, na.rm = TRUE)
# df <- list_new[[1]]
# ggplot(df, aes(x = time, y = speed_new, colour = food_journey)) +
#   geom_line() +
#   labs(title = "Speed Over Time", x = "Time (s)", y = "Speed (units/s)")
# 
# ggplot(df, aes(x = food_journey, y = speed_new, fill = food_journey)) +
#   geom_boxplot() +
#   labs(title = "Speed Distribution by Food Journey State", x = "Food Journey State", y = "Speed (units/s)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve label readability
# 
# # speed_summary <- df %>%
# #   filter(!(food_journey %in% c("arrival", "departure", "at_food"))) %>% 
# #   group_by(food_journey) %>%
# #   summarise(mean_speed = mean(speed, na.rm = TRUE),
# #             sd_speed = sd(speed, na.rm = TRUE),
# #             se_speed = sd_speed / sqrt(n()))

# calculate_speed_summary <- function(df) {
#   df %>%
#     filter(food_journey %in% c("trip_to", "trip_back")) %>%
#     group_by(unique_trial_ID, food_journey) %>%
#     summarise(mean_speed = mean(speed_new, na.rm = TRUE), .groups = 'drop') %>%
#     pivot_wider(names_from = food_journey, values_from = mean_speed,
#                 names_prefix = "speed_") %>%
#     ungroup()
# }

calculate_speed_summary <- function(df) {
  df %>%
    filter(food_journey %in% c("trip_to", "trip_back", "exploration")) %>%
    group_by(unique_trial_ID, food_journey) %>%
    summarise(mean_speed = mean(speed_new, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(
      names_from = food_journey, 
      values_from = mean_speed,
      names_prefix = "speed_", 
      values_fill = list(mean_speed = NA)
    ) %>%
    ungroup() %>%
    # Add missing columns explicitly if they do not exist
    mutate(
      speed_trip_to = if("speed_trip_to" %in% names(.)) speed_trip_to else NA,
      speed_trip_back = if("speed_trip_back" %in% names(.)) speed_trip_back else NA,
      speed_exploration = if("speed_exploration" %in% names(.)) speed_exploration else NA
    )
}

speed_summaries <- lapply(list_new, calculate_speed_summary)
speed_summary <- bind_rows(speed_summaries) %>%
  select(unique_trial_ID, speed_trip_to, speed_trip_back, speed_exploration)
max(speed_summary$speed_trip_to, na.rm = TRUE)
mean(speed_summary$speed_trip_to, na.rm = TRUE)

ggplot(speed_summary %>%
         pivot_longer(cols = c(speed_trip_to, speed_trip_back, speed_exploration),
                      names_to = "speed_type", values_to = "speed"),
       aes(x = speed, fill = speed > 50)) +
  geom_histogram(binwidth = 5, color = "black") +
  facet_wrap(~speed_type, scales = "free") +  # One plot per speed type
  scale_fill_manual(values = c("lightblue", "red"), labels = c("Normal", "Outlier")) +
  labs(title = "Speed Distribution for Different Types", x = "Speed (cm/s)", y = "Count") +
  theme_minimal() +
  theme(legend.title = element_blank())

ggplot(speed_summary %>%
         filter(speed_trip_to < 50), aes(x = speed_trip_to)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
  labs(title = "Filtered Speed Trip To Distribution", x = "Speed (cm/s)", y = "Count") +
  theme_minimal()

speed_summary %>%
  filter(speed_trip_to > 100) %>%
  distinct(unique_trial_ID)
speed_summary %>%
  filter(speed_trip_to > 1000) %>%
  distinct(unique_trial_ID)

#exclude these two trials and plot again

speed_summary %>%
  filter(!unique_trial_ID %in% c("winter_20201031-1_T8", "winter_20240303-3_T6", "winter_20240225-1_T7", "winter_20240228-1_T10", "winter_20240228-2_T9", "winter_20240228-3_T4", "winter_20240303-3_T6", "winter_20240303-3_T5")) %>%
  pivot_longer(cols = c(speed_trip_to, speed_trip_back, speed_exploration),
               names_to = "speed_type", values_to = "speed") %>%
  ggplot(aes(x = speed, fill = speed > 100)) +
  geom_histogram(binwidth = 5, color = "black") +
  facet_wrap(~speed_type, scales = "free") +  # One plot per speed type
  scale_fill_manual(values = c("lightblue", "red"), labels = c("Normal", "Outlier")) +
  labs(title = "Speed Distribution for Different Types (Excluding Outliers)", 
       x = "Speed (cm/s)", y = "Count")


tracked_speed <- do.call(rbind, list_new)
rownames(tracked_speed) <- NULL

speed_summary <- tracked_speed %>%
  filter(food_journey %in% c("trip_to", "trip_back", "exploration")) %>% 
  group_by(food_journey, unique_trial_ID) %>%
  summarise(mean_speed = mean(speed_new, na.rm = TRUE),
            sd_speed = sd(speed_new, na.rm = TRUE),
            se_speed = sd_speed / sqrt(n()),
            count = n(),
            .groups = 'drop') %>% 
  filter(mean_speed < 40)

ggplot(speed_summary, aes(x = food_journey, y = mean_speed, fill = food_journey)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_speed - sd_speed, ymax = mean_speed + se_speed), width = 0.2) +
  labs(title = "Mean Speed by Food Journey State with Error Bars", x = "Food Journey State", y = "Mean Speed (units/s)")

ggplot(speed_summary, aes(x = mean_speed, fill = food_journey)) +
  geom_density(alpha = 0.4) +
  labs(x = "Speed (cm/s)", y = "Density") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~status)
ggplot(speed_summary, aes(x = status, y = mean_speed)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1")

#Previous door#################

other_door_visits <- read.csv("/home/cbaldoni/data/cue_learning/other_door_visit.csv", header = TRUE) %>% 
  mutate(unique_trial_ID = paste(season, ID, trial, sep = "_"),
         unique_trial_ID = as.factor(unique_trial_ID),
         trial_n = as.integer(str_remove(trial, "T")),
         door_ID = as.factor(door_ID),
         sequence = as.factor(if_else(str_starts(unique_trial_ID, "winter_2024"), "winter_2024", "other")))
doors <- read.csv("~/data/cue_learning/trial_door.csv") %>%
  mutate(door = as.factor(door),
         sequence = as.factor(season)) %>%
  relocate(sequence, .before = season) %>%
  select(-season)

doors <- doors %>%
  arrange(trial_n, sequence) %>%
  group_by(sequence) %>%
  mutate(door = as.factor(door),
         previous_door_ID = ifelse(trial_n == 1, as.character(door), as.character(lag(door)))) %>%
  ungroup()
doors <- doors %>%
  mutate(previous_door_ID = as.factor(previous_door_ID))

result <- other_door_visits %>%
  mutate(previous_trial_n = trial_n - 1)  %>% 
  left_join(doors, by = c("previous_trial_n" = "trial_n", "sequence")) %>%
  mutate(door_match = door_ID == previous_door_ID,
         season = as.factor(season),
         status = case_when(season == "summer" ~ "summer_wild",
                            season == "winter" & grepl("^2021", ID) ~ "winter_wild",
                            season == "winter" & grepl("^2020", ID) ~ "winter_captive",
                            season == "winter" & grepl("^2024", ID) ~ "winter_wild",
                            season == "spring" & grepl("^2024", ID) ~ "spring_wild",
                            season == "spring" & grepl("^2020", ID) ~ "spring_captive"),
         status = as.factor(status),
         food_journey = as.factor(food_journey)) %>% 
  filter(!is.na(door_match)) %>% 
  select(- sequence, -previous_trial_n)


#relevel the factors, t10 comes always after t1
#result$trial <- factor(result$trial, levels = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"), ordered = TRUE)
ggplot(result, aes(x = trial_n, fill = door_match)) +
  geom_bar(stat = "count") +  # Counts the occurrences
  scale_x_continuous(breaks = pretty_breaks()) +
  #scale_fill_brewer(palette = "Set1") +
  facet_grid(food_journey ~ status, 
             labeller = labeller(
               food_journey = c(trip_to = "Trip To", at_food = "At Food", trip_back = "Trip Back"),
               status = c(summer_wild = "Summer Wild", winter_captive = "Winter Captive", winter_wild = "Winter Wild", spring_captive = "Spring Captive", spring_wild = "Spring Wild")
             )) +
  labs(x = "Trial", y = "Count", fill = "Door Match Status") +
  theme_bw()


library(BayesFactor)

table_door_status <- table(result$door_match, result$status)
bf_door_status <- contingencyTableBF(table_door_status, sampleType = "indepMulti", fixedMargin = "rows")
print(bf_door_status)

table_door_food_journey <- table(result$door_match, result$food_journey)
bf_door_food_journey <- contingencyTableBF(table_door_food_journey, sampleType = "indepMulti", fixedMargin = "rows")
print(bf_door_food_journey)


model_door <- brm(formula = bf(door_match ~ status + (1 | ID)),
             data = result, family = bernoulli("logit"), 
             prior = c(prior(exponential(0.65), class = "sd"),
                       prior(normal(0,1), class ="b")),
             chains = 4, iter = 4000, cores = 4,
             control = list(adapt_delta = 0.95, max_treedepth = 15))

model_door <- brm(formula = bf(door_match ~ status*food_journey + (1 | ID)),
                  data = result, family = bernoulli("logit"), 
                  prior = c(prior(exponential(0.65), class = "sd"),
                            prior(normal(0,1), class ="b")),
                  chains = 4, iter = 4000, cores = 4,
                  control = list(adapt_delta = 0.95, max_treedepth = 15))


#####TURNING ANGLE##################
calculate_turns <- function(x, y, window_size = 3) {
  n <- length(x)
  # Calculate base angles using a sliding window approach
  angles <- atan2(y[window_size:n] - y[1:(n-window_size+1)], 
                  x[window_size:n] - x[1:(n-window_size+1)])
  # Initialize a vector of NAs for turning angles
  turning_angles <- rep(NA, n)
  # Fill in the calculated angles appropriately offset to align with the midpoint of the window
  turning_angles[ceiling(window_size / 2):(n - floor(window_size / 2))] <- angles
  # Calculate the difference in angles to get turning
  turning_angle_diffs <- c(NA, diff(turning_angles))
  # Scale to degrees
  turning_angle_diffs <- turning_angle_diffs * 180 / pi
  return(turning_angle_diffs)
}

for (i in seq_along(list_new)) {
  df <- list_new[[i]]
  turning_angles <- calculate_turns(df$x, df$y)
  df$turning_angle <- turning_angles
  # Define abrupt turns based on a threshold, e.g., 95th percentile
  threshold <- quantile(abs(turning_angles), 0.95, na.rm = TRUE)
  df$abrupt_turns <- abs(turning_angles) > threshold
  list_new[[i]] <- df
}




####PLOTS####

library(viridis)
library(ggsci)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

#there are some trials with more than 20 revisits, I want to filter them
df_revisits <- dplyr::filter(tracking_results, grepl('revisit_', food_journey))
df_revisits['unique_trial_ID']

#List based on unique_trial_ID
list_new
#select.list(tracking_results, dplyr::starts_with("spring_T10_20201106-1"))

#b <-list_new[["spring_20240520-4_T8"]]
#b$food_journey <- factor(b$food_journey, levels = c("trip_to", "arrival", "at_food", "departure", "trip_back", "exploration"), ordered = TRUE)
plotfood <- coords %>%
  as.data.frame(plotfood, row.names = NULL) %>% 
  filter(unique_trial_ID == unique(b$unique_trial_ID)) %>%
  dplyr::select(c("FOOD_x", "FOOD_y", "unique_trial_ID"))
plotdoor <- coords %>%
  filter(unique_trial_ID == unique(b$unique_trial_ID)) %>%
  select(4:11, 15)
# plot <- 
b %>% ggplot(aes(x, y, colour = food_journey)) +
  ggtitle(b$unique_trial_ID) +
  geom_point(x = plotdoor$A_x, y = plotdoor$A_y,  size = 4, colour = "black") +
  geom_point(x = plotdoor$B_x, y = plotdoor$B_y,  size = 4, colour = "black") +
  geom_point(x = plotdoor$C_x, y = plotdoor$C_y,  size = 4, colour = "black") +
  geom_point(x = plotdoor$D_x, y = plotdoor$D_y,  size = 4, colour = "black") +
  geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 16, colour = "darkgreen", alpha = 1/20) +
  geom_point(x = plotfood$FOOD_x, y = plotfood$FOOD_y,  size = 8, colour = "green") +
  geom_point(size = 2) +
  theme(legend.position="none")

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
    geom_path() + 
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

path_data <- read.csv("~/data/cue_learning/distance/master_distance.csv") %>% 
  mutate(trial = str_remove(trial, "^T"))
path_data$trial <- str_sort(path_data$trial, numeric = TRUE)
path_data$trial <- as.integer(path_data$trial)
str(path_data$trial)
ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = food_door, color = season))+
  geom_smooth(mapping = aes(x = trial, y = walked_to, color = status)) +
  facet_wrap(~status)
  #scale_color_manual(values = c("green", "coral", "cornflowerblue"))

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  #geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season)) +
  scale_color_manual(values = c("green", "coral", "cornflowerblue")) + 
  facet_wrap(~ID)

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  #geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  geom_point(mapping = aes(x = as.numeric(trial), y = straightness_exploration, color = season)) +
  scale_color_manual(values = c("green", "coral", "cornflowerblue")) + 
  facet_wrap(~ID)

ggplot(data = path_data) + 
  geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season), se = FALSE) +
  geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season), alpha = 0.3) +
  scale_color_manual(values = c("green", "coral", "cornflowerblue"), name = "Season") +
  scale_x_continuous(limits = c(1, 10), breaks = seq(1, 10, 1), name = "Trial Number") +
  scale_y_continuous(name = "Straightness (to food)") +
  labs(title = "Straightness to Food by Trial Number")

# + 
  #facet_wrap(~season)

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_back, color = season))+
  geom_smooth(mapping = aes(x = as.numeric(trial), y = straightness_back, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))
ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_back, color = season))+
  geom_smooth(mapping = aes(x = as.numeric(trial), y = straightness_exploration, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

