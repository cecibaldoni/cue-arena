# Libraries ----

library(tidyverse)
library(here)
library(sf) #for spatial manipulation
library(zoo) #performs linear interpolation on NA values
library(parallel)
library(ggplot2)
library(trajr)
library(plotly)
library(RColorBrewer)
library(scales)

# Load files ----
doors <- read.csv(here("trial_door.csv")) %>%
  mutate(Trial = paste0("T", trial_n))

coords <- read.csv(here("coords.csv"))
tracking <- read.csv(here("tracking.csv"))

#empty dataframe for saving info on visit to other doors on trip back
other_door_visits <- data.frame(ID = NULL, door = NULL, Trial = NULL, Season = NULL) #we can also add the length of time spent at the door if needed
no_visits <- data.frame(unique_trial_ID = NULL, other_door_visits = NULL, ID = NULL, trial = NULL, season = NULL)

trial_ls <- split(tracking, tracking$unique_trial_ID)

#change Inf into NA and interpolate x and y (x_approx and y_approx)
trial_ls <- lapply(trial_ls, function(df) {
  df <- df %>%
    mutate(x = ifelse(is.infinite(x) | x == "inf", NA, as.numeric(x)),
           y = ifelse(is.infinite(y) | y == "inf", NA, as.numeric(y))) %>%
    mutate(x = na.approx(x, na.rm = FALSE),
           y = na.approx(y, na.rm = FALSE)) %>%
    mutate(frame = seq(0, nrow(df) - 1),
           time = frame / 30) %>%
    mutate(distance = sqrt((x - lag(x))^2 + (y - lag(y))^2), speed = distance * 30) %>%
    mutate(speed = ifelse(is.na(speed), 0, speed)) %>%
    select(-distance) %>%
    relocate(time, .before = frame)
  return(df)
})

#prepare cluster for parallel computation
mycl <- makeCluster(4) #the number of CPUs to use (adjust this based on your machine)
clusterExport(mycl, c("coords", "trial_ls", "doors", "other_door_visits")) #define the variable that will be used within the ParLapply call
clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
})

# how to check single trials and plot them (e.g. if coordinate system is correct)
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
  geom_point(data = y 
             #%>% filter(frame <= 15) #You can filter frames of interest
             , aes(x = x, y = y, color = frame, text = paste("X:", x, "Y:", y))) +  
  geom_point(data = doors_coords, aes(x=x, y=y, text = paste("Door", door_ID)), colour = "green", size = 5) +
  geom_point(data = food_coords, aes(x = FOOD_x, y = FOOD_y, text = paste("FOOD_X:", FOOD_x, "FOOD_Y:", FOOD_y)), colour = "red", size = 5) +
  ggtitle(paste("Trial ID:", y$unique_trial_ID[1]))
ggplotly(p, tooltip = "text")

## Start Loop through all df in trial_ls ----

#start empty csv files to store the data

results <- here("processed", "master_results.csv")
distance <-  here("processed", "master_distance.csv")

#start lapply function
other_door_visits_ls <- lapply(trial_ls, function(x){
  ## to call all blocks one at a time
  #x = trial_ls[["spring_20240520-4_T8"]]
  ## Or:
  #x= trial_ls[[1]]
  
  #extract food coordinates for this trial AND convert to a sf object
  # Extract food coordinates and buffer, convert to sf object
  food_coords <- coords %>%
    filter(unique_trial_ID == unique(x$unique_trial_ID)) %>%
    dplyr::select(c("FOOD_x", "FOOD_y")) %>%
    st_as_sf(coords = c("FOOD_x", "FOOD_y"))
  food_buffer <- st_buffer(food_coords, dist = 4)
  
  # Determine season
  season_filter <- if (any(str_starts(x$unique_trial_ID, "winter_2024"))) "winter_2024" else "other"
  
  # Door data
  trial_door_ID <- doors %>%
    filter(Trial == unique(x$trial), season == season_filter) %>%
    pull(door)
  
  # Extract and convert all doors
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
  
  # Convert track into an sf object
  track_sf <- x %>%
    st_as_sf(coords = c("x", "y"))
  
  # Find food intersections
  #WARNING MESSAGE
  at_food <- track_sf %>%
    st_intersection(food_buffer) %>%
    as.data.frame() %>%
    arrange(frame) %>% #arrange by time/frame
    mutate(timediff = frame - lag(frame)) %>%
    mutate(new_timediff = ifelse(is.na(timediff) | timediff >= 30, 1, 0)) %>%
    mutate(visit_seq = cumsum(new_timediff))
  
  # Add food journey info to x IF food was reached
  if (nrow(at_food) == 0) {
    print(paste0("Food out of buffer zone in ", unique(x$unique_trial_ID)))
    return(NULL)
  }
  
  ## Label food journey ----
  track_sf_2 <- track_sf %>%
    full_join(at_food[c("frame", "visit_seq")]) %>%
    arrange(frame) %>%
    mutate(old_food_journey = case_when(
      frame == head(at_food$frame, 1) ~ "arrival",
      frame == tail(at_food[at_food$visit_seq == 1, "frame"], 1) ~ "departure",
      frame < head(at_food$frame, 1) ~ "trip_to",
      between(frame, head(at_food[at_food$visit_seq == 1, "frame"], 1), tail(at_food[at_food$visit_seq == 1, "frame"], 1)) ~ "at_food",
      frame %in% at_food[at_food$visit_seq != 1, "frame"] ~ paste("trip_back_revisit", visit_seq, sep = "_"),
      TRUE ~ "trip_back"
    ))
  
  # "being at the exit" sequences, and use this information to label "exploration"
  # Find overlap exit and track
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
      mutate(food_journey = ifelse(frame > tail(at_exit[at_exit$exit_seq == 1, "frame"], 1), 
                                   "exploration", as.character(old_food_journey)),
             food_journey = ifelse(food_journey == "exploration" & !is.na(visit_seq),
                                   paste0("exploration_revisit_", visit_seq), food_journey)) %>%
      select(-exit_seq)
    
  } else { #if at_exit does not exist, keep the old column
    track_sf_2 <- track_sf_2 %>%
      mutate(food_journey = old_food_journey)
  }
  
  # Save CSV-ready version
  track_save <- track_sf_2 %>%
    mutate(coordinates = st_coordinates(geometry)) %>%
    mutate(x = coordinates[, "X"], y = coordinates[, "Y"]) %>%
    relocate(x, .after = frame) %>% 
    relocate(y, .after = x) %>% 
    relocate(unique_trial_ID, .before = season) %>% 
    select(-old_food_journey, -visit_seq, -coordinates) %>% 
    st_drop_geometry()
  
  write.table(track_save, file = results, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(results))
  
  ## Summary metrics ----
  
  x$food_journey <- track_sf_2$food_journey
  
  get_traj_metrics <- function(trip_df) {
    traj <- TrajFromCoords(trip_df, fps = 30, spatialUnits = "cm")
    list(length = TrajLength(traj),
         straightness = {s <- TrajStraightness(traj); if (length(s) == 0) NA_real_ else s},
         traj = traj)
  }
  
  trip_to <- x %>% filter(food_journey == "trip_to") %>% select(x, y, time)
  trip_back <- x %>% filter(food_journey == "trip_back") %>% select(x, y, time)
  exploration <- x %>% filter(food_journey == "exploration") %>% select(x, y, time)
  
  to_metrics <- get_traj_metrics(trip_to)
  back_metrics <- get_traj_metrics(trip_back)
  exploration_metrics <- get_traj_metrics(exploration)
  
  dist_doorfood <- if (!is.null(to_metrics$traj)) {
    TrajDistance(to_metrics$traj, startIndex = 1, endIndex = nrow(to_metrics$traj))
  } else NA_real_
  
  time_summary <- x %>%
    group_by(unique_trial_ID) %>%
    summarize(
      time_to_food = if (any(food_journey == "trip_to")) max(time[food_journey == "trip_to"], na.rm = TRUE) else NA_real_,
      time_back = if (any(food_journey == "trip_back")) max(time[food_journey == "trip_back"], na.rm = TRUE) else NA_real_,
      time_exploration = if (any(food_journey == "exploration")) diff(range(time[food_journey == "exploration"], na.rm = TRUE)) else NA_real_,
      time_at_food = {
        at_food_times <- time[food_journey == "at_food"]
        if (length(at_food_times) > 1) max(at_food_times) - min(at_food_times)
        else if (length(at_food_times) == 1) 1 / 30 else 0
      },
      .groups = "drop"
    )
  
  df <- data.frame(
    unique_trial_ID = as.factor(x$unique_trial_ID[1]),
    season = as.factor(x$season[1]),
    ID = as.factor(x$ID[1]),
    trial = as.factor(x$trial[1]),
    status = as.factor(x$status[1]),
    food_door = dist_doorfood,
    walked_to = to_metrics$length,
    walked_back = back_metrics$length,
    walk_exploration = exploration_metrics$length,
    straightness_to_food = to_metrics$straightness,
    straightness_back = back_metrics$straightness,
    straightness_exploration = exploration_metrics$straightness,
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(time_summary %>% select(-unique_trial_ID))
  
  write.table(df, file = distance, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(distance))
  
  # Check for other door visits
  other_doors <- track_sf_2 %>%
    filter(food_journey %in% c("trip_back", "exploration")) %>%
    st_intersection(all_doors_buffer %>% filter(door_ID != trial_door_ID))
  
  if (nrow(other_doors) > 0) {
    new_visits <- other_doors %>%
      group_by(door_ID) %>% slice(1) %>%
      dplyr::select(c("ID", "season", "trial", "door_ID", "food_journey")) %>%
      st_drop_geometry()
    other_door_visits <<- rbind(other_door_visits, new_visits)
    return(new_visits)
  } else {
    empty_data <- data.frame(unique_trial_ID = unique(x$unique_trial_ID), other_door_visits = 0) %>%
      separate(unique_trial_ID, into = c("season", "trial", "ID"), sep = "_")
    no_visits <<- rbind(no_visits, empty_data)
  }
  print(paste0("trial ", unique(x$unique_trial_ID), " completed."))
})

  
write.csv(other_door_visits, here("processed", "other_door_visit.csv"), row.names = FALSE)

stopCluster(mycl)

saveRDS(other_door_visits_ls, file = here("trials_sp.rds"))
saveRDS(other_door_visits, file= here("other_door_visits.rds"))

tracking_results <- read.csv(here("processed", "master_results.csv")) %>%
  mutate(season = as.factor(season),
         ID = as.factor(ID),
         status = as.factor(status),
         food_journey = as.factor(food_journey),
         trial = as.factor(trial),
         unique_trial_ID =as.factor(unique_trial_ID)) %>% 
  #filter(!is.na(season) & season != "1") %>% 
  droplevels()

trial_list <- split(tracking_results, tracking_results$unique_trial_ID)


# Speed ----


# Previous door ----

other_door_visits <- read.csv(here("processed", "other_door_visit.csv"), header = TRUE) %>% 
  mutate(unique_trial_ID = paste(season, ID, trial, sep = "_"),
         unique_trial_ID = as.factor(unique_trial_ID),
         trial_n = as.integer(str_remove(trial, "T")),
         door_ID = as.factor(door_ID),
         sequence = as.factor(if_else(str_starts(unique_trial_ID, "winter_2024"), "winter_2024", "other")))
doors <- read.csv(here("trial_door.csv")) %>%
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


## Plots ----

library(viridis)
library(ggsci)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

#there are some trials with more than 20 revisits, I want to filter them
df_revisits <- dplyr::filter(tracking_results, grepl('revisit_', food_journey))
df_revisits['unique_trial_ID']

#List based on unique_trial_ID
trial_list
#select.list(tracking_results, dplyr::starts_with("spring_T10_20201106-1"))

b <-trial_list[["summer_20200623-1_T4"]]
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

