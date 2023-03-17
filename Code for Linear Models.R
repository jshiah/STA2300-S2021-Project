
##--------------------------------------
## Install and load any needed libraries
##--------------------------------------

library(lubridate)


##--------------------------------------
## Load the data
##--------------------------------------

data <- read.csv("Mines_Park_pH_Fault.csv", stringsAsFactors = FALSE)
attr(data$dateTime, "tzone") <- "America/Denver"

head(data)

colnames(data)
dim(data)
head(data$dateTime)
tail(data$dateTime)


#April 24, 2010, 10:00 AM.

month.val <- month(data$dateTime)
day.val <- day(data$dateTime)
hour.val <- hour(data$dateTime)
fault.time <- which(month.val==4 & day.val==24 & hour.val==10)[1]
data$dateTime[fault.time]


######################
# Plot of pH
# The biological community cannot survive if
# the pH is less than 6.5
# It clearly dips below this value
######################

ts.plot(data$ras_ph, xlab="Index", ylab="pH")
abline(h=6.5,col=4,lwd=2)
abline(v=fault.time,col=2,lwd=2)

######################
# Plot of MBR Flow
######################
ts.plot(data$mbr_1_perm_flow, xlab="Index", ylab="MBR 1, Permeate Flow")
abline(v=fault.time,col=2,lwd=2)

ts.plot(data$mbr_2_perm_flow, xlab="Index", ylab="MBR 2, Permeate Flow")
abline(v=fault.time,col=2,lwd=2)
abline(v=fault.time,col=2,lwd=2)





######################
# Sorting the variables into
# groups
######################
## cyclic variables
cyclic_vars <- c("cos_daily", 
                     "sin_daily", 
                     "cos_2hour", 
                     "sin_2hour", 
                     "cos_hourly", 
                     "sin_hourly")
##control variables
control_vars <- c("bio_1_blow_flow", 
                  "bio_2_blow_flow", 
                  "mbr_1_tmp", 
                  "mbr_2_tmp",
                  "ambient_temp", 
                  "bio_1_phase_1", 
                  "bio_1_phase_2", 
                  "bio_2_phase_1",
                  "bio_2_phase_2", 
                  "mbr_1_mode_1", 
                  "mbr_1_mode_2", 
                  "mbr_1_mode_4",
                  "mbr_2_mode_1", 
                  "mbr_2_mode_2", 
                  "mbr_2_mode_4")

## response variables
response_vars <- c("mbr_1_perm_flow", 
                   "mbr_2_perm_flow", 
                   "ras_temp",
                   "bio_1_do", 
                   "bio_2_do", 
                   "mbr_1_level", 
                   "mbr_2_level", 
                   "perm_turb", 
                   "sewage_flow", 
                   "bio_1_level", 
                   "bio_2_level",
                   "bio_1_temp", 
                   "bio_2_temp", 
                   "bio_1_tss", 
                   "bio_2_tss", 
                   "perm_tank_level", 
                   "ras_do", 
                   "ras_ph", 
                   "perm_cond", 
                   "ras_tss")
                   
#Other predictors                  
scale_predictors <- c("bio_1_blow_flow",
                      "bio_2_blow_flow", 
                      "ambient_temp",
                      "mbr_1_tmp", 
                      "mbr_2_tmp")

x.cols<-NULL
for(i in 1:length(scale_predictors)){
	x.cols<-c(x.cols,which(colnames(data)==scale_predictors[i]))}
for(i in 1:length(control_vars)){
	x.cols<-c(x.cols,which(colnames(data)==control_vars[i]))}
for(i in 1:length(cyclic_vars)){
	x.cols<-c(x.cols,which(colnames(data)==cyclic_vars[i]))}

y.cols<-NULL
for(i in 1:length(response_vars)){
	y.cols<-c(y.cols,which(colnames(data)==response_vars[i]))}


XX<-data[,x.cols]
YY<-data[,y.cols]


obs_hour <- hour(data$dateTime)
obs_minute <- minute(data$dateTime)
data$cos_2hour <- cos((obs_hour + obs_minute/60)*(360/24)*pi/180)
data$sin_2hour <- sin((obs_hour + obs_minute/60)*(360/24)*pi/180)







proportion <- 0.8

range <- c(3084:10084)

n <- length(range)


index <- floor(n *  proportion)

train_inputs <- XX[range[1]: range[index],]
test_inputs <- XX[range[index+1]: range[n], ]

train_labels <- YY[range[1]: range[index], "ras_ph"]
test_labels <- YY[range[index+1]: range[n], "ras_ph"]


nrow(train_inputs)
nrow(test_inputs)




ts.plot(YY$ras_ph, xlab="Index", ylab="pH")
abline(v=range[1],col=2,lwd=2)
abline(v=range[index],col=2,lwd=2)
abline(v=range[n],col=2,lwd=2)

predictionpH <- lm(train_labels ~  mbr_1_mode_1 +mbr_1_mode_2 + mbr_1_mode_4 , data=train_inputs)




predictedvalues <- predict(predictionpH, test_inputs)
ts.plot(test_labels, xlab="Index", ylab="pH" )
lines(predictedvalues, col="red")



n <- nrow(train_inputs) 
BICbackwardmodel <- step(predictionpH, direction="backward", k=log(n), trace = 0)
BICpredict <- predict(BICbackwardmodel, test_inputs)
ts.plot(test_labels, xlab="Index", ylab="pH" )
lines(BICpredict, col="red")

summary(BICbackwardmodel)

n <- nrow(train_inputs) 
BICforwardmodel <- step(predictionpH, direction="forward", k=log(n), trace = 0)
BICfpredict <- predict(BICforwardmodel, test_inputs)
ts.plot(test_labels, xlab="Index", ylab="pH" )
lines(BICfpredict, col="red")

summary(BICforwardmodel)

backwardpredictions <- predict(BICbackwardmodel, test_inputs)
ts.plot(test$ras_ph, xlab="Index", ylab="pH" )
lines(backwardpredictions, col="red")



suppressMessages(library(glmnet))
leave_out <- which((colnames(train) == "ras_ph")| (colnames(train) == "dateTime"))

x <- as.matrix(train[,-leave_out])
y <- train$ras_ph

lasso_fit <- cv.glmnet(x, y)

coef(lasso_fit)

lassopredictions <- predict(lasso_fit, newx=x)
ts.plot(test$ras_ph, xlab="Index", ylab="pH" )
lines(lassopredictions, col="red")


BICresiduals <- residuals(BICbackwardmodel)


BIC_RMSE <- sqrt(mean(BICresiduals^2))
BIC_MAE <- mean(abs(BICresiduals))


BIC_RMSE
BIC_MAE


proportion <- 0.25

fault_range <- c(fault.time:(fault.time + 7000))
n_fault <- length(fault_range)
n_fault
end_index <- floor(proportion * n_fault)
end_index


fault_train_inputs <- XX[fault_range[1]:fault_range[end_index], ]
fault_test_inputs <- XX[fault_range[end_index + 1]:fault_range[n_fault],]

fault_train_labels <- YY[fault_range[1]:fault_range[end_index], "ras_ph"]
fault_test_labels <- YY[fault_range[end_index + 1]:fault_range[n_fault], "ras_ph"]

fault_model <- lm(fault_train_labels ~., data = fault_train_inputs)

fault_prediction <- predict(fault_model, fault_test_inputs)
ts.plot(fault_test_labels, xlab="Index", ylab="pH" )
lines(fault_prediction, col="red")


Back_fault_model <- step(fault_model, direction = "backward", k=log(n_fault), trace=0)

fault_prediction <- predict(Back_fault_model, fault_test_inputs)
ts.plot(fault_test_labels, xlab="Index", ylab="pH" )
lines(fault_prediction, col="red")

summary(Back_fault_model)


BICresiduals <- residuals(Back_fault_model)


BIC_RMSE <- sqrt(mean(BICresiduals^2))
BIC_MAE <- mean(abs(BICresiduals))


BIC_RMSE
BIC_MAE

