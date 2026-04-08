rm(list=ls())
setwd('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Analysis')

#####################################################################################################

d <- read.csv('~/PhD_project/project_PONDERFUL/ponderful_DATA/ponderful_DATA_updated/PhyChe_XY_landcover_nov.txt')

library(dplyr)
d <- d %>%
  mutate(TP.t = log10(TP),  
         TN.t = log10(TN))  

d2 <- d[,c("TN.t","TP.t","bio1","bio4","bio5","bio6","bio7","bio12","bio15","bio17",
                                 'Area','Depth','Hydeoperiod_length',"Pond_dries", "Aquatic_500",
                                 "Cropland_500","Forest_500","Pastures.and.open.nature_500" ,
                                 "Urban_500", "Natural_5",'Animals_cont')]
model_df <- model_df %>%
  mutate(TP.t = log10(TP),  
         TN.t = log10(TN))  

model_df_2 <- model_df[,c("TN.t","TP.t",'bio1.s','bio7.s','bio5.s','bio12.s','Hydeoperiod_length.s','Animals_cont.s','Area.s','Fish.s','Depth.s',"Aquatic_500.s",   
                     "Cropland_500.s","Forest_500.s","Pastures.and.open.nature_500.s" ,
                     "Urban_500.s", "Natural_5.s")]

colnames(d2) <-c("TN","TP","Annual Mean Temperature","Temperature Seasonality","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range",
"Annual Precipitation","Precipitation Seasonality","Precipitation of Driest Quarter","Area","Depth","Hydroperiod length","Temporality","Aquatic","Cropland",
"Forest","Pastures and open nature","Urban","Natural 5","Livestock index")


library(gbm)
library(dismo)
library(ggplot2)

####################################################
#TN
####################################################

str(model_df_2)
d2N <- na.omit(d2)
response <- model_df_2$TN.t            
predictors <- model_df_2[, 3:17]    
predictors <- d2N[, c('Cropland','Depth','Annual Precipitation')]    

data <- data.frame(response, predictors)

set.seed(123)  

brt_TN <- gbm.step(
  data = data,
  gbm.x = 2:ncol(data),  # Predictor columns
  gbm.y = 1,            # Response variable column (TN)
  family = "gaussian",    
  tree.complexity = 3,   # Number of nodes in each tree
  learning.rate = 0.001, 
  bag.fraction = 0.75    # Proportion of data used for training at each step
)

brt_simp <- gbm.simplify(brt_TN, n.drops = 16)


print(brt_TN)
print(brt_TN$n.trees)
summary(brt_TN)
dev.off()
gbm.plot(brt_TN, n.plots = 3, write.title=FALSE)  
gbm.plot.fits(brt_TN,1)
gbm.plot(brt_TN, i.var = 1, j.var = 2)  # Interaction between two predictors

plot(brt_TN$residuals)  # Check residuals

find.int <- gbm.interactions(brt_TN)
find.int$interactions



explained_deviance <- (brt_TN$self.statistics$mean.null - brt_TN$self.statistics$mean.resid) / brt_TN$self.statistics$mean.null * 100


var_importance <- summary(brt_TN)

library(ggplot2)
ggplot(var_importance, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance in BRT Model", 
       x = "Predictor Variables", 
       y = "Relative Importance (%)") +
  theme_minimal()


###################################################
#TP
####################################################

d2P <- d2[!is.na(d2$TP), ]
response <- d2P$TP            
predictors <- d2P[, 3:21]     
data <- data.frame(response, predictors)

set.seed(123)  

brt_TP <- gbm.step(
  data = data,
  gbm.x = 2:ncol(data),  # Predictor columns
  gbm.y = 1,            # Response variable column (TP)
  family = "gaussian",    
  tree.complexity = 3,   # Number of nodes in each tree
  learning.rate = 0.001, 
  bag.fraction = 0.75    # Proportion of data used for training at each step
)

print(brt_TP)
print(brt_TP$n.trees)
summary(brt_TP)

gbm.plot(brt_TP, n.plots = 19, smooth = TRUE)  

explained_deviance <- (brt_TP$self.statistics$mean.null - brt_TP$self.statistics$mean.resid) / brt_TP$self.statistics$mean.null * 100


var_importance <- summary(brt_TP)

library(ggplot2)
ggplot(var_importance, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance in BRT Model", 
       x = "Predictor Variables", 
       y = "Relative Importance (%)") +
  theme_minimal()


#Simplified model
simplified_model <- gbm.simplify(
  brt_TP, 
  n.folds = 10 # Use the same cross-validation folds as in your original model
)

# Inspect the reduced predictor set
simplified_model$predictor.list

library(rpart)
library(rpart.plot)

# Fit a single decision tree
tree_model <- rpart(response ~ ., data = data)

# Plot the tree
rpart.plot(tree_model, type = 3, extra = 101, under = TRUE)



#######################
#RANDOM FOREST
#######################

####TN
library(randomForest)

d2N <- d2[!is.na(d2$TN), ]
response_TN <- d2N$TN  
predictors <- d2N[, 3:21]  
data_TN <- data.frame(response_TN, predictors)

set.seed(123)  
rf_TN <- randomForest(response_TN ~ ., data = data_TN, ntree = 500, mtry = 4)
print(rf_TN)
importance(rf_TN)



####TP
d2P <- d2[!is.na(d2$TP), ]
response_TP <- d2P$TP  
predictors <- d2P[, 3:21]  
data_TP <- data.frame(response_TP, predictors)

set.seed(123) 
rf_TP <- randomForest(response_TP ~ ., data = data_TP, ntree = 500, mtry = 4)
print(rf_TP)
importance(rf_TP)

# Plot the variable importance 
varImpPlot(rf_TN)
varImpPlot(rf_TP)

# Plot predictions vs actual values 
plot(response_TN, predict(rf_TN), main = "Random Forest for TN", xlab = "Actual", ylab = "Predicted")
plot(response_TP, predict(rf_TP), main = "Random Forest for TP", xlab = "Actual", ylab = "Predicted")


### Tree RF TN
single_tree <- getTree(rf_TN, k = 1, labelVar = TRUE)
print(single_tree)

library(rpart)
library(rpart.plot)

# Fit a simple rpart tree for visualization (using TN as the response)
rpart_tree <- rpart(response_TN ~ ., data = data_TN, method = "anova")

# Plot the decision tree
rpart.plot(rpart_tree, type = 4, extra = 101, fallen.leaves = TRUE, main = "Decision Tree for TN")

# View the splits and thresholds in the tree
print(rpart_tree)

# Loop through all trees in the random forest model
for (i in 1:rf_TN$ntree) {
  tree <- getTree(rf_TN, k = i, labelVar = TRUE)
  print(tree)  # View the splits for the i-th tree
}

