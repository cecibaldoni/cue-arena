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
  #x = trial_ls[[1]]
  
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

write.csv(other_door_visits, "/home/ceci/data/cue_learning/other_door_visit.csv", sep = ",", row.names = FALSE)

staSys.time() - b
stopCluster(mycl)

#mapview(trial_door_coords) + mapview(food_buffer) + mapview(food_coords) + mapview(trial_door_buffer) + mapview(track_sf)

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

# Define threshold value for turning angle
angle_threshold <- 85

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
  i <- trial_ls[[i]]
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

# Loop through each dataframe in the list
for (i in 1:length(trial_ls)) {
  # Extract x, y, and frame columns from the dataframe
  x <- trial_ls[[i]]$x
  y <- trial_ls[[i]]$y
  frame <- trial_ls[[i]]$frame
  
  # Calculate speed based on x, y, and frame
  x_diff <- c(0, diff(x))
  y_diff <- c(0, diff(y))
  frame_diff <- c(0, diff(frame))
  distance <- sqrt(x_diff^2 + y_diff^2)
  speed <- distance / frame_diff
  
  # Add speed as a new column to the dataframe
  trial_ls[[i]]$speed_new <- speed
  
  # Calculate expected speed based on frame difference
  expected_speed <- distance / frame_diff
  # Calculate the threshold as the 90th percentile of the actual speed
  #it returns the value below which 90% of the data fall
  speed_threshold <- quantile(speed, 0.95, na.rm = TRUE)
  
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
library(brms)

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

#relevel the factors, for some reason t10 comes always after t1
#result$trial <- factor(result$trial, levels = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"), ordered = TRUE)
result$season <- factor(result$season, levels = c("summer", "winter", "spring"))
ggplot(result, aes(x = trial, #fill = factor(previous_door_ID))
                   )) +
  geom_bar() +
  facet_wrap(~ season) +
  scale_fill_manual(values = c("red", "darkgreen")) +
  xlab("Trial") +
  ylab("Count") +
  ggtitle("Visit to previous door, by Season and Trial")


prior_prev <- get_prior(previous_door_ID ~ season + trial + (season + trial | ID), data = result, family = bernoulli(link="logit"))

model_prevdoor <- brm(
  previous_door_ID ~ season + trial + (season + trial | ID),  # Dependent variable, predictors, and random slopes
  data = result,  # Dataset
  family = bernoulli(link="logit"),  # Bernoulli distribution for binary outcomes
  prior = prior_prev,  # Specify priors
  iter = 2000,  # Number of MCMC iterations
  chains = 4,  # Number of MCMC chains
  warmup = 1000,  # Number of warmup iterations
  control = list(adapt_delta = 0.99)  # Control parameters for the Stan sampler
)
summary(model_prevdoor)
plot(model_prevdoor)
pp_check(model_prevdoor)
proportion <- function(x){mean(x == 1)}
pp_check(model_prevdoor,
         plotfun = "stat", stat = "proportion") +
  xlab("proportion")

model_prev2 <- brm(
  previous_door_ID ~ season + trial + (1|ID),  # Dependent variable, predictors, and random effect
  data = result,  # Dataset
  family = bernoulli(link ="logit"),  # Bernoulli distribution for binary outcomes
  prior = c(set_prior("normal(0, 1)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b")),  # Specify priors
  iter = 2000,  # Number of MCMC iterations
  chains = 4,  # Number of MCMC chains
  warmup = 1000,  # Number of warmup iterations
  control = list(adapt_delta = 0.99)  # Control parameters for the Stan sampler
)

summary(model_prev2)
plot(model_2)
pp_check(model_2)

result$previous_door_ID <- as.numeric(result$previous_door_ID)

formula3 <- bf(previous_door_ID ~ trial + season)

model3 <- brm(formula3, data = result, family = bernoulli(link = "logit"), 
             prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                       set_prior("normal(0, 1)", class = "b")), 
             sample_prior = TRUE, cores = 4, chains = 4,
             iter = 4000, warmup = 1000)
summary(model3)
plot(model3)
pp_check(model3)

loo(model_prevdoor,model_prev2, model3, compare = TRUE)


library(bayesplot)

# Get the predicted probabilities from the model
newdata <- expand.grid(
  season = levels(result$season),
  trial = unique(result$trial)
)
pred <- posterior_predict(model3, newdata = newdata)
pred
hdi <- apply(pred, 1, hdi, prob = 0.95)

# Plot the predicted probabilities and the raw data
ggplot(result, aes(x = trial, y = previous_door_ID, color = season)) +
  geom_jitter(position = position_jitter(0.1)) +
  geom_smooth(aes(y = mean(pred), group = season)) +
  geom_ribbon(
     aes(ymin = hdi[, 1], ymax = hdi[, 2], group = season),
     alpha = 0.2
  ) +
  labs(
    x = "Trial",
    y = "Probability of Previous Door ID being 1",
    color = "Season"
  ) +
  theme_bw()

#Extract the posterior samples
post_samples <- posterior_samples(model3)
#summary() Calculate the posterior probabilities of the effects of "season" and "trial" on "previous_door_ID"
#Plot the posterior distributions of the effects of "season" and "trial" on "previous_door_ID" 
plot(model3, pars = c("b_season", "b_trial"))
library(tidybayes)
library(bayesplot)
library(ggdist)
library(bayesplot)


result$season <- factor(result$season)
result$ID <- factor(result$ID)
# Create a plot of the observed data
ggplot(result, aes(x = trial, y = previous_door_ID, color = season)) +
  geom_jitter(position = position_jitter(width = 0.01), alpha = 0.5) +
  scale_color_manual(values = c("#33a02c", "#e31a1c", "#1f78b4"), name = "Season") +
  labs(x = "Trial", y = "Previous Door ID") +
  theme_bw()

# Create a model plot with posterior predictions
posterior_predict <- posterior_predict(model_prevdoor, draws = 1000)

model_plot <- pp_check(model_prevdoor, ndraws = 100)

proportion <- function(x){mean(x == 1)}
pp_check(model_prevdoor,
         plotfun = "stat", stat = "proportion") +
  xlab("proportion")
summary(model_prevdoor)

result$trial <- factor(result$trial, levels = levels(model_prevdoor$fixed[[1]]$trial))


result_pp <- add_predicted_draws(object = model_prevdoor, newdata = result)


# Create a scatter plot with predicted values
model_fit <- result_pp %>%
  ggplot(aes(x = trial, y = .prediction, color = season)) +
  geom_point(position = position_jitter(width = 0.01), alpha = 0.5) +
  scale_color_manual(values = c("#33a02c", "#e31a1c", "#1f78b4"), name = "Season") +
  labs(x = "Trial", y = "Previous Door ID") +
  theme_bw()
model_fit

# Add the 95% credible interval as error bars
model_fit + geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.1)

model_plot <- model_plot +
  geom_jitter(data = result, mapping = aes(x = trial, y = previous_door_ID, color = season),
              position = position_jitter(width = 0.1), alpha = 0.5)

# Display the model plot
model_plot

str(result)



##############
#STRAIGHTNESS#
##############
library(brms)
library(rstanarm)
library(bayesplot)

path_data <- read_csv("/home/ceci/data/cue_learning/distance/distance.csv") %>% 
  select(-X)

path_data$trial <- gsub("T", "", path_data$trial) %>% 
  as.numeric(path_data$trial)
path_data <- path_data[order(path_data$trial), ]

path_data$season <- factor(path_data$season, levels = c("summer", "winter", "spring"))

# Define the model formula
model_formula2 <- bf(straightness_to_food ~ trial + (trial | ID))
prior2 <- get_prior(model_formula2, data = path_data, family = gaussian())

model2 <-  brm(straightness_to_food ~ trial + (trial | ID), 
               data = path_data,
               prior = prior,
               chains = 4, iter = 2000, warmup = 1000)

#model straightness index as a function of trial and season, with random intercepts for trial and season
model3 <- brm(straightness_to_food ~ trial + season + (1|trial) + (1|season), data = path_data)

summary(model2) #check for convergence to make sure that the MCMC chains have stabilized
traceplot(model2) #should show stable chains without obvious patterns

loo(model2,model3, compare = TRUE)
#model2 is slightly better
#the difference will be negative if the predictive accuracy of the first model is higher

##PIECEWISE MODEL
#Create a binary variable indicating whether trial is greater than 6 or not.
path_data$trial_binary <- ifelse(path_data$trial > 6, 1, 0)
subset <- path_data[path_data$trial_binary==0, ]
#new model with breakpoint
#trial:trial_binary  is the interaction term between trial and trial_binary
model_formula4 <- bf(straightness_to_food ~ trial + trial_binary + season + trial:trial_binary + (1|trial) + (1|season))
prior4 <- get_prior(model_formula4, data = path_data, family = gaussian())
model4 <- brm(straightness_to_food ~ trial + trial_binary + season + trial:trial_binary + (1|trial) + (1|season), 
              data = path_data[path_data$trial_binary==0, ],
              iter = 10000,
              warmup = 2000,
              control = list(adapt_delta = 0.95),
              prior = prior4)
summary(model4)

model_formula5 <- bf(straightness_to_food ~ trial + season + (1|trial) + (1|season))
prior4 <- get_prior(model_formula4, data = subset, family = gaussian())
model4 <- brm(straightness_to_food ~ trial + trial_binary + season + trial:trial_binary + (1|trial) + (1|season), 
              data = path_data[path_data$trial_binary==0, ],
              iter = 10000,
              warmup = 2000,
              control = list(adapt_delta = 0.95),
              prior = prior4)
summary(model4)

plot(model4, pars = c("Intercept", "trial", "trial_binary", "season", "trial:trial_binary"))
#generate predicted values for the straightness index using the modified model
predictions <- posterior_predict(model4, newdata = expand.grid(trial = 1:10, trial_binary = c(0, 1), season = c("summer", "winter", "spring")))
#reshape prediction dataframe
predictions <- predictions %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "trial_season_binary", values_to = "predicted_straightness_index") %>%
  separate(trial_season_binary, into = c("trial", "trial_binary", "season"), sep = "_")

#Plot the predicted straightness index for each trial and season combination, separating the lines for trial_binary = 0 and trial_binary = 1
ggplot(predictions, aes(x = trial, y = predicted_straightness_index, color = season, linetype = factor(trial_binary))) +
  geom_line() +
  geom_point() +
  labs(x = "Trial", y = "Predicted straightness index", color = "Season", linetype = "Trial binary")

#You can use posterior predictive checks to assess whether 
#the model adequately captures the patterns in the data. This involves comparing 
#the observed data to the data generated by the model using the posterior distribution. 
#If the model fits the data well, the simulated data should look similar to the observed data.
pp_check(model2, ndraws= 100)
    

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = food_door, color = season))+
  geom_smooth(mapping = aes(x = as.numeric(trial), y = walked_to, color = season))+
  scale_color_manual(values = c("green", "coral", "cornflowerblue"))

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  geom_smooth(mapping = aes(x = trial, y = straightness_to_food, color = season))+
  scale_color_manual(values = c("coral", "cornflowerblue", "green"))

ggplot(data = path_data) + 
  #geom_point(mapping = aes(x = trial, y = straightness_back, color = season))+
  geom_smooth(mapping = aes(x = trial, y = straightness_back, color = season))+
  scale_color_manual(values = c("coral", "cornflowerblue", "green"))


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
# df_revisits['unique_trial_ID']
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
  i = track_all[[c("spring_T1_20201131-1")]]
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
