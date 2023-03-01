doors <- read.csv("~/Documents/data/cue_learning/trial_door.csv") %>% 
  mutate(Trial = paste0("T",trial_n))

coords <- read.csv("~/Documents/data/cue_learning/food_door_coordinates.csv", sep = ";",header = TRUE) %>%
  mutate(Trial = paste0("T",TRIAL)) %>%
  mutate(unique_trial_ID = paste(SEASON, Trial, ID, sep = "_"))

tracking <- lapply(list.files("~/Documents/data/cue_learning/tracking", full.names = T), read.csv) %>%
  reduce(rbind) %>%
  drop_na(frame) %>%
  mutate(unique_trial_ID = paste(Season, Trial, ID, sep = "_"))

newdata=c()
track_dynamics_data=c()
trialIDs=unique(tracking$unique_trial_ID)

steplengths=c()
for(i in 1:length(trialIDs)){
  trial_needed=trialIDs[i]
  temp_track=tracking[which(tracking$unique_trial_ID==trial_needed),]
  temp_track$z=temp_track$x+1i*temp_track$y
  
  steps=Mod(diff(temp_track$z))
  steps=as.numeric(na.omit(steps))
  if(sum(is.infinite(steps))>0){
    steps=steps[-which(is.infinite(steps)==TRUE)]
  }
  if(sum(is.nan(steps))>0){
    steps=steps[-which(is.nan(steps)==TRUE)]
  }
  
  steplengths[i]=list(steps)

}
steplengths=unlist(steplengths)
meanstep=mean(steplengths,na.rm=TRUE)

b <- Sys.time()
for(i in 1:length(trialIDs)){
  trial_needed=trialIDs[i]
  temp_track=tracking[which(tracking$unique_trial_ID==trial_needed),]
  temp_track$z=temp_track$x+1i*temp_track$y
  
  if(length(which(coords$unique_trial_ID==trial_needed))==0){
    print(paste("Trial missing:",trial_needed, sep = " " ))
    next}
  doors_temp=coords[which(coords$unique_trial_ID==trial_needed),]
  doors_temp$A_x=doors_temp$A_x#*doors_temp$cmxpixel fixed in the excel file (?)
  doors_temp$A_y=doors_temp$A_y#*doors_temp$cmxpixel
  doors_temp$B_x=doors_temp$B_x#*doors_temp$cmxpixel
  doors_temp$B_y=doors_temp$B_y#*doors_temp$cmxpixel
  doors_temp$C_x=doors_temp$C_x#*doors_temp$cmxpixel
  doors_temp$C_y=doors_temp$C_y#*doors_temp$cmxpixel
  doors_temp$D_x=doors_temp$D_x#*doors_temp$cmxpixel
  doors_temp$D_y=doors_temp$D_y#*doors_temp$cmxpixel
  doors_temp$FOOD_x=doors_temp$FOOD_x#*doors_temp$cmxpixel
  doors_temp$FOOD_y=doors_temp$FOOD_y#*doors_temp$cmxpixel
  #some unique trial ids are in "tracking" but not in coords
  
  Food=doors_temp$FOOD_x+1i*doors_temp$FOOD_y
  temp_track$distance_to_food=temp_track$z-Food
  temp_track$distance_to_food=Mod(temp_track$distance_to_food)
  journey=temp_track[which(temp_track$distance_to_food<=3.5),]
  journey=journey[which(journey$frame==min(journey$frame,na.rm=TRUE)),]
  temp_track$journey=NA
  
  temp_track$journey[which(temp_track$frame<unique(journey$frame))]="trip_to"
  temp_track$journey[which(temp_track$distance_to_food<=3.5)]="at_food"
  temp_track$journey[which(temp_track$frame>=unique(journey$frame))]="trip_back"
  
  # temp_track$at_food=temp_track$journey
  # temp_track$at_food[which(temp_track$distance_to_food<=3.5)]="At food"
  # 
  #FIZ THIS MESS
  # rleinfo=rle(temp_track$at_food)
  # rleinfo=data.frame(cbind(rleinfo$lengths, rleinfo$values))
  # rleinfo$start=NA
  # rleinfo$end=NA
  # rleinfo$X1=as.numeric(rleinfo$X1)
  # rleinfo$start=as.numeric(rleinfo$start)
  # rleinfo$end=as.numeric(rleinfo$end)
  # 
  # rleinfo$start[1]=1
  # rleinfo$end[1]=rleinfo$X1[1]-1
  # 
  # for(j in 2:nrow(rleinfo)){
  #   rleinfo$start[j]=rleinfo$end[j-1]+1
  #   rleinfo$end[j]=rleinfo$start[j]+rleinfo$X1[j]
  # }
  # importantones=rleinfo[which(rleinfo$X2=="At food"),]
  # 
  # for(k in 1:nrow(importantones)){
  #   index=seq(from=importantones$start[k], to = importantones$end[k], by=1)
  #   temp_track$at_food[index]=paste("At food", k, sep = "_")
  # }

  newdata[i]=list(temp_track)
  to_plot=ggplot()+geom_path(data=temp_track, aes(x=x, y=y, color=frame))+
    ggtitle(paste(trial_needed))+
    # annotate("text", label = "A",
    #   x=doors_temp$A_x, y=doors_temp$A_y,  size = 8, colour = "red") +
    # annotate("text", label = "B",
    #   x=doors_temp$B_x, y=doors_temp$B_y, size = 8, colour = "red") +
    # annotate("text", label = "C",
    #   x=doors_temp$C_x, y=doors_temp$C_y, size = 8, colour = "red") +
    # annotate("text", label = "D",
    #   x= doors_temp$D_x, y=doors_temp$D_y, size = 8, colour = "red") +
    geom_point(data=doors_temp, aes(x=FOOD_x, y=FOOD_y, size=2), color=I("orange"))+
    theme_classic()+ coord_equal()
  print(to_plot)
  
  
  # before_track=temp_track[which(temp_track$arrival=="before"),]
  # after_track=temp_track[which(temp_track$arrival=="after"),]
  # coords_before=before_track[,c(4,6,7)]
  # coords_after=after_track[,c(4,6,7)]
  # before_track2=TrajFromCoords(coords_before)
  # after_track2=TrajFromCoords(coords_after)
  # before_track2=TrajRediscretize(before_track2,meanstep)
  # after_track2=TrajRediscretize(after_track2,meanstep)
  # 
  # sinuosity_before=TrajSinuosity2(before_track2)
  # sinuosity_after=TrajSinuosity2(after_track2)
  # 
  # length_before=sum(Mod(diff(before_track$z)),na.rm=TRUE)
  # length_after=sum(Mod(diff(after_track$z)),na.rm=TRUE)
  # 
  # results=data.frame(cbind(trial_needed,length_before, length_after,sinuosity_before,sinuosity_after))
  # track_dynamics_data=list(results)
}

Sys.time() - b 
newdata=do.call(rbind,newdata)
track_dynamics_data=do.call(rbind,track_dynamics_data)

write.csv(newdata, "all_data.csv")
write.csv(track_dynamics_data, "all_dynamic_data.csv")


