library(sf)
library(tidyverse)

doors_coords
#assuming doors coords is defined, in my scenario is just a tibble with one triaL:
# > str(doors_coords)
# tibble [4 × 3] (S3: tbl_df/tbl/data.frame)
# $ door_ID: chr [1:4] "A" "B" "C" "D"
# $ y      : num [1:4] 170 104 60 126
# $ x      : num [1:4] 169 126 190 236


edges <- lapply(trial_ls, function(x){
  x = trial_ls[[1]]
  
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
  
  track_sf <- x %>%
    st_as_sf(coords = c("x", "y"))
  
  side_length <- 114
  calculate_intersection <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
    denom <- (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if (denom == 0) {
      return(c(NA, NA)) #this is for parallel lines
    }
    ua <- ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denom
    x <- x1 + ua * (x2 - x1)
    y <- y1 + ua * (y2 - y1)
    return(c(x, y))
  }
  
  calculate_perpendicular_vector <- function(x1, y1, x2, y2) {
    dx <- x2 - x1
    dy <- y2 - y1
    mag <- sqrt(dx^2 + dy^2)
    if (mag == 0) {
      return(c(0, 0)) # avoid division by zero
    }
    unit_vector <- c(dx / mag, dy / mag)
    # Perpendicular components
    perp_dx <- -unit_vector[2]
    perp_dy <- unit_vector[1]
    return(c(perp_dx, perp_dy))
  }
  
  # intersection of diagonals A-C and B-D
  intersection <- calculate_intersection(doors_coords$x[1], doors_coords$y[1], doors_coords$x[3], doors_coords$y[3], # A-C
                                         doors_coords$x[2], doors_coords$y[2], doors_coords$x[4], doors_coords$y[4])  # B-D
  
  
  # find endpoints with perpendiculars and stop at intersections
  endpoints_list <- list()
  corner_points <- list()
  for (i in 1:nrow(doors_coords)) {
    next_i <- ifelse(i == nrow(doors_coords), 1, i + 1)
    
    center_x <- doors_coords$x[i]
    center_y <- doors_coords$y[i]
    next_center_x <- doors_coords$x[next_i]
    next_center_y <- doors_coords$y[next_i]
    
    # perpendicular vector for current door
    perp_vector <- calculate_perpendicular_vector(center_x, center_y, intersection[1], intersection[2])
    
    # potential endpoints for the current side
    shift <- side_length / 2
    x1 <- center_x - shift * perp_vector[1]
    y1 <- center_y - shift * perp_vector[2]
    x2 <- center_x + shift * perp_vector[1]
    y2 <- center_y + shift * perp_vector[2]
    
    # intersection with the next side perpendicular line
    next_perp_vector <- calculate_perpendicular_vector(next_center_x, next_center_y, intersection[1], intersection[2])
    x3 <- next_center_x - shift * next_perp_vector[1]
    y3 <- next_center_y - shift * next_perp_vector[2]
    x4 <- next_center_x + shift * next_perp_vector[1]
    y4 <- next_center_y + shift * next_perp_vector[2]
    
    intersection_point <- calculate_intersection(x1, y1, x2, y2, x3, y3, x4, y4)
    
    if (!any(is.na(intersection_point))) {
      corner_points[[i]] <- intersection_point
    } else {
      corner_points[[i]] <- c(x2, y2) # second endpoint as fallback
    }
  }
  
  corner_points_df <- do.call(rbind, corner_points)
  side_lines <- st_multilinestring(list(
    rbind(corner_points_df[1,], corner_points_df[2,]),
    rbind(corner_points_df[2,], corner_points_df[3,]),
    rbind(corner_points_df[3,], corner_points_df[4,]),
    rbind(corner_points_df[4,], corner_points_df[1,])
  )) %>% st_sfc() %>% st_sf()
  
  buffer_distance <- 4
  edges_buffer <- st_buffer(side_lines, dist = buffer_distance)
  
  at_edge <- track_sf %>%
    st_intersection(edges_buffer) %>%
    as.data.frame() %>%
    arrange(frame) #arrange by time/frame
  
  x <- x %>%
    mutate(at_edge = as.vector(st_intersects(track_sf, edges_buffer, sparse = FALSE)))
  
  time_at_edge <- x %>%
    filter(at_edge == TRUE) %>%
    summarize(total_time = sum(time))
  
  calculate_distance <- function(data) {
    data <- data %>%
      arrange(frame) %>%
      mutate(
        x_lag = lag(x),
        y_lag = lag(y),
        dist = sqrt((x - x_lag)^2 + (y - y_lag)^2)
      )
    total_distance <- sum(data$dist, na.rm = TRUE)
    return(total_distance)
  }
  
  at_edge_cm <- calculate_distance(x %>% filter(at_edge == TRUE))
  out_edge_cm <- calculate_distance(x %>% filter(at_edge == FALSE))

  df = data.frame(x[1,1],x[1,2],x[1,3], x[1,4],
                  time_at_edge, at_edge_cm, out_edge_cm)
  colnames(df) <- c(
    "season", "ID", "trial", "status",
    "time_at_edge", "at_edge_cm", "out_edge_cm"
    )
  
  write.csv(df, file = paste0("~/data/cue_learning/distance/", unique(x$unique_trial_ID),".csv")) #CHANGE DIRECTORY
  
  # 
  # ggplot() +
  #   geom_sf(data = track_sf, aes(geometry = geometry, color = "Track")) +
  #   geom_sf(data = edges_buffer, aes(geometry = geometry, fill = "Edges Buffer"), alpha = 0.3) +
  #   geom_sf(data = at_edge, aes(geometry = geometry, color = "At Edge")) +
  #   scale_color_manual(values = c("Track" = "blue", "At Edge" = "red")) +
  #   scale_fill_manual(values = c("Edges Buffer" = "gray"))
  
  })




plot(st_geometry(buffer), col = 'lightblue', border = NA, main = "Arena with Side Buffer")
plot(st_geometry(side_lines), col = 'transparent', border = 'black', add = TRUE)
points(closed_corner_points[,1], closed_corner_points[,2], pch=16, col="red", cex=1.5)
points(doors_coords$x, doors_coords$y, pch=16, col="blue", cex=1.5)
text(doors_coords$x, doors_coords$y, labels=doors_coords$door_ID, pos=3, col="blue", cex=1.2)

legend("topright", legend=c("Center Points", "Corner Points", "Buffer"), col=c("blue", "red", "lightblue"), pch=16, bg="white")

