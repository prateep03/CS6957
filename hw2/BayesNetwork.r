source("~/CS6957/hw2/BayesianNetworks_SD.r");
###############Question 1###############################################################################
meta_root="~/CS6957/hw2"
df <- read.table(file.path(meta_root,"RiskFactorData.csv"),header=TRUE,sep=",",colClasses=c("integer","integer","integer","integer","integer","integer","integer","integer","integer","integer"))
income = createCPT.fromData(df,c("income"))
bmi.income_exercise=createCPT.fromData(df,c("bmi","income","exercise"))
exercise.income=createCPT.fromData(df,c("exercise","income"))
bp.exercise_income_smoke=createCPT.fromData(df,c("bp","exercise","income","smoke"))
smoke.income=createCPT.fromData(df,c("smoke","income"))
diabetes.bmi=createCPT.fromData(df,c("diabetes","bmi"))
stroke.bmi_bp_cholesterol=createCPT.fromData(df,c("stroke","bmi","bp","cholesterol"))
attack.bmi_bp_cholesterol=createCPT.fromData(df,c("attack","bmi","bp","cholesterol"))
angina.bmi_bp_cholesterol=createCPT.fromData(df,c("angina","bmi","bp","cholesterol"))
cholesterol.exercise_income_smoke=createCPT.fromData(df,c("cholesterol","exercise","income","smoke"))

bayesNet = list("income" = income, "bmi" = bmi.income_exercise,"exercise"= exercise.income,
                "bp"=bp.exercise_income_smoke,"smoke"=smoke.income,"diabetes"=diabetes.bmi,
                "stroke"=stroke.bmi_bp_cholesterol,"attack"=attack.bmi_bp_cholesterol,
                "angina"=angina.bmi_bp_cholesterol,"cholesterol"=cholesterol.exercise_income_smoke)

##############################################Question 2####################################################################
## Probability of health outcomes if I have bad habits 
print(infer(bayesNet, c("income","bmi","bp","stroke","attack","angina","cholesterol"), c("smoke","exercise"), c(1,2)))
print(infer(bayesNet, c("income","bmi","bp","diabetes","attack","angina","cholesterol"), c("smoke","exercise"), c(1,2)))
print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","angina","cholesterol"), c("smoke","exercise"), c(1,2)))
print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","attack","cholesterol"), c("smoke","exercise"), c(1,2)))
## Probability of health outcomes if I have good habits
print(infer(bayesNet, c("income","bmi","bp","stroke","attack","angina","cholesterol"), c("smoke","exercise"), c(2,1)))
print(infer(bayesNet, c("income","bmi","bp","diabetes","attack","angina","cholesterol"), c("smoke","exercise"), c(2,1)))
print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","angina","cholesterol"), c("smoke","exercise"), c(2,1)))
print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","attack","cholesterol"), c("smoke","exercise"), c(2,1)))
## Probability of health outcomes if I have poor health
print(infer(bayesNet, c("income","stroke","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
print(infer(bayesNet, c("income","diabetes","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
print(infer(bayesNet, c("income","diabetes","stroke","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
print(infer(bayesNet, c("income","diabetes","stroke","attack","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
## Probability of health outcomes if I have good health
print(infer(bayesNet, c("income","stroke","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
print(infer(bayesNet, c("income","diabetes","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
print(infer(bayesNet, c("income","diabetes","stroke","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
print(infer(bayesNet, c("income","diabetes","stroke","attack","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
###############################################Question 3#################################################################
non3<-function()
{
## diabetes
d=c();
for(i in 1:8)
{
  diabetes_res=infer(bayesNet,c("bmi","stroke","attack","exercise","smoke","bp","angina","cholesterol"),c("income"),c(i))  
  res<-diabetes_res$probs[1]
  d[i]<-res[[1]]
}

## stroke
s=c();
for(i in 1:8)
{
  stroke_res=infer(bayesNet,c("bmi","diabetes","attack","exercise","smoke","bp","angina","cholesterol"),c("income"),c(i))
  res<- stroke_res$probs[1]
  s[i]<-res[[1]]
}
## heart attack
h=c();
for(i in 1:8)
{
  attack_res=infer(bayesNet,c("bmi","diabetes","stroke","exercise","smoke","bp","angina","cholesterol"),c("income"),c(i))
  res<- attack_res$probs[1]
  h[i]<-res[[1]]
}
## angina
a=c();
for(i in 1:8)
{
  angina_res=infer(bayesNet,c("bmi","diabetes","stroke","exercise","smoke","bp","attack","cholesterol"),c("income"),c(i))
  res<- angina_res$probs[1]
  a[i]<-res[[1]]
}

x<-c(1,2,3,4,5,6,7,8)
plot(x,d,type='o', col='black', lwd=3, ylim=c(0.0,0.15),xlab ="status of income",ylab=expression(paste("p(y=", 1, "|income=",i,")")))
lines(s, type="o", col='red', lwd=3)
lines(h, type="o", col='blue', lwd=3)
lines(a, type="o", col="green", lwd=3)
legend('topright',c("diabetes","stroke","attack","angina"),
       col = c("black","red","blue","green"))
}
############################################Question 4######################################################################################
non1<-function()
{
income = createCPT.fromData(df,c("income"))
bmi.income_exercise=createCPT.fromData(df,c("bmi","income","exercise"))
exercise.income=createCPT.fromData(df,c("exercise","income"))
bp.exercise_income_smoke=createCPT.fromData(df,c("bp","exercise","income","smoke"))
smoke.income=createCPT.fromData(df,c("smoke","income"))
diabetes.bmi_smoke_exercise=createCPT.fromData(df,c("diabetes","bmi","smoke","exercise"))
stroke.bmi_bp_cholesterol_smoke_exercise=createCPT.fromData(df,c("stroke","bmi","bp","cholesterol","smoke","exercise"))
attack.bmi_bp_cholesterol_smoke_exercise=createCPT.fromData(df,c("attack","bmi","bp","cholesterol","smoke","exercise"))
angina.bmi_bp_cholesterol_smoke_exercise=createCPT.fromData(df,c("angina","bmi","bp","cholesterol","smoke","exercise"))
cholesterol.exercise_income_smoke=createCPT.fromData(df,c("cholesterol","exercise","income","smoke"))

bayesNet_old = list("income" = income, "bmi" = bmi.income_exercise,"exercise"= exercise.income,
                "bp"=bp.exercise_income_smoke,"smoke"=smoke.income,"diabetes"=diabetes.bmi_smoke_exercise,
                "stroke"=stroke.bmi_bp_cholesterol_smoke_exercise,"attack"=attack.bmi_bp_cholesterol_smoke_exercise,
                "angina"=angina.bmi_bp_cholesterol_smoke_exercise,"cholesterol"=cholesterol.exercise_income_smoke)
}
##########################Queries in 2#####################################
## Probability of health outcomes if I have bad habits 
#print(infer(bayesNet, c("income","bmi","bp","stroke","attack","angina","cholesterol"), c("smoke","exercise"), c(1,2)))
#print(infer(bayesNet, c("income","bmi","bp","diabetes","attack","angina","cholesterol"), c("smoke","exercise"), c(1,2)))
#print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","angina","cholesterol"), c("smoke","exercise"), c(1,2)))
#print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","attack","cholesterol"), c("smoke","exercise"), c(1,2)))
## Probability of health outcomes if I have good habits
#print(infer(bayesNet, c("income","bmi","bp","stroke","attack","angina","cholesterol"), c("smoke","exercise"), c(2,1)))
#print(infer(bayesNet, c("income","bmi","bp","diabetes","attack","angina","cholesterol"), c("smoke","exercise"), c(2,1)))
#print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","angina","cholesterol"), c("smoke","exercise"), c(2,1)))
#print(infer(bayesNet, c("income","bmi","bp","diabetes","stroke","attack","cholesterol"), c("smoke","exercise"), c(2,1)))
## Probability of health outcomes if I have poor health
#print(infer(bayesNet, c("income","stroke","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
#print(infer(bayesNet, c("income","diabetes","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
#print(infer(bayesNet, c("income","diabetes","stroke","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
#print(infer(bayesNet, c("income","diabetes","stroke","attack","exercise","smoke"), c("bp","cholesterol","bmi"), c(1,1,3)))
## Probability of health outcomes if I have good health
#print(infer(bayesNet, c("income","stroke","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
#print(infer(bayesNet, c("income","diabetes","attack","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
#print(infer(bayesNet, c("income","diabetes","stroke","angina","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
#print(infer(bayesNet, c("income","diabetes","stroke","attack","exercise","smoke"), c("bp","cholesterol","bmi"), c(3,2,2)))
############################################Question 5######################################################################################
none<-function()
{
income = createCPT.fromData(df,c("income"))
bmi.income_exercise=createCPT.fromData(df,c("bmi","income","exercise"))
exercise.income=createCPT.fromData(df,c("exercise","income"))
bp.exercise_income_smoke=createCPT.fromData(df,c("bp","exercise","income","smoke"))
smoke.income=createCPT.fromData(df,c("smoke","income"))
diabetes.bmi_smoke_exercise=createCPT.fromData(df,c("diabetes","bmi","smoke","exercise"))
stroke.bmi_bp_cholesterol_smoke_exercise_diabetes=createCPT.fromData(df,c("stroke","bmi","bp","cholesterol","smoke","exercise","diabetes"))
attack.bmi_bp_cholesterol_smoke_exercise=createCPT.fromData(df,c("attack","bmi","bp","cholesterol","smoke","exercise"))
angina.bmi_bp_cholesterol_smoke_exercise=createCPT.fromData(df,c("angina","bmi","bp","cholesterol","smoke","exercise"))
cholesterol.exercise_income_smoke=createCPT.fromData(df,c("cholesterol","exercise","income","smoke"))

bayesNet = list("income" = income, "bmi" = bmi.income_exercise,"exercise"= exercise.income,
                "bp"=bp.exercise_income_smoke,"smoke"=smoke.income,"diabetes"=diabetes.bmi_smoke_exercise,
                "stroke"=stroke.bmi_bp_cholesterol_smoke_exercise_diabetes,"attack"=attack.bmi_bp_cholesterol_smoke_exercise,
                "angina"=angina.bmi_bp_cholesterol_smoke_exercise,"cholesterol"=cholesterol.exercise_income_smoke)

result_new<-infer(bayesNet,c("bmi","income","attack","exercise","smoke","bp","angina","cholesterol"),c("diabetes"),c(1))
print(result_new)
result_old<-infer(bayesNet_old,c("bmi","income","attack","exercise","smoke","bp","angina","cholesterol"),c("diabetes"),c(1))
print(result_old)

result_new<-infer(bayesNet,c("bmi","income","attack","exercise","smoke","bp","angina","cholesterol"),c("diabetes"),c(3))
print(result_new)
result_old<-infer(bayesNet_old,c("bmi","income","attack","exercise","smoke","bp","angina","cholesterol"),c("diabetes"),c(3))
print(result_old)
}
##########################################################END##################################################################################