rm(list = ls())
source('~/CS6957/hw2/BayesianNetworks.r');

##########################################################################
## "RiskFactors.csv"
##########################################################################
data = read.csv("RiskFactorData.csv", header = TRUE)

dataNet = list(createCPT.fromData(data,"income"),
               createCPT.fromData(data,c("smoke","income")),
               createCPT.fromData(data,c("bmi","income","exercise")),
               createCPT.fromData(data,c("exercise","income")),
               createCPT.fromData(data,c("bp","exercise","income","smoke")),
               createCPT.fromData(data,c("cholesterol","smoke","income","exercise")),
               createCPT.fromData(data,c("diabetes","bmi")),
               createCPT.fromData(data,c("stroke","bmi","bp","cholesterol")),
               createCPT.fromData(data,c("attack","bmi","bp","cholesterol")),
               createCPT.fromData(data,c("angina","bmi","bp","cholesterol")))

##############################################
## Q 2

# question (2) diabets
margvar = c('stroke','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q2(a) diabets, bad habits, ")
print(result)

margvar = c('stroke','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q2(a) diabets, good habits, ")
print(result)

margvar = c('stroke','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q2(b) diabets, poor health, ")
print(result)

margvar = c('stroke','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q2(b) diabets, good health, ")
print(result)

# question (2) stroke
margvar = c('diabetes','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q2(a) stroke, bad habits, ")
print(result)

margvar = c('diabetes','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q2(a) stroke, good habits, ")
print(result)

margvar = c('diabetes','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q2(b) stroke, poor health, ")
print(result)

margvar = c('diabetes','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q2(b) stroke, good health, ")
print(result)

# question (2) heart attack
margvar = c('stroke','diabetes','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q2(a) attack, bad habits, ")
print(result)

margvar = c('stroke','diabetes','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q2(a) attack, good habits, ")
print(result)

margvar = c('stroke','diabetes','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q2(b) attack, poor health, ")
print(result)

margvar = c('stroke','diabetes','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q2(b) attack, good health, ")
print(result)

# question (2) angina
margvar = c('stroke','diabetes','attack', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q2(a) angina, bad habits, ")
print(result)

margvar = c('stroke','diabetes','attack', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q2(a) angina, good habits, ")
print(result)

margvar = c('stroke','diabetes','attack', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q2(b) angina, poor health, ")
print(result)

margvar = c('stroke','diabetes','attack', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q2(b) angina, good health, ")
print(result)

####===================================================================================
# question 3
var = c('angina', 'attack', 'diabetes', 'stroke', 'exercise', 'smoke', 'bmi', 'bp', 
        'cholesterol', 'income')
observar = c('income')

# diabetes
margvar = setdiff(var, c('diabetes', 'income'))
ret.diabetes = income.healthOutcome(dataNet, margvar, observar)

# stroke
margvar = setdiff(var, c('stroke', 'income'))
ret.stroke = income.healthOutcome(dataNet, margvar, observar)

# heat attack
margvar = setdiff(var, c('attack', 'income'))
ret.attack = income.healthOutcome(dataNet, margvar, observar)

# angina
margvar = setdiff(var, c('angina', 'income'))
ret.angina = income.healthOutcome(dataNet, margvar, observar)
# plot out the results
#pdf('3_plots.pdf')

plot(ret.diabetes, type='o', col='black', lwd=3, ylim=c(0.0, 0.15), 
     xlab = "level of income",
     ylab = expression(paste("p(", "y=", 1, "|", "income", "=", i, ")"))
     )
lines(ret.stroke, type='o', col='red', lwd=3)
lines(ret.attack, type='o', col='blue', lwd=3)
lines(ret.angina, type='o', col='green', lwd=3)
legend('topright', c("diabetes", "stroke", "attack", "angina"),
       col = c("black", "red", "blue", "green"), lwd = 3)
#dev.off()

####===================================================================================
# question 4
dataNet = list(createCPT.fromData(data, c("angina", "bmi", "bp", "cholesterol", "smoke", "exercise")),
               createCPT.fromData(data, c("attack", "bmi", "bp", "cholesterol", "smoke", "exercise")),
               createCPT.fromData(data, c("bmi", "exercise", "income")),
               createCPT.fromData(data, c("bp", "exercise", "income", "smoke")),
               createCPT.fromData(data, c("cholesterol", "exercise", "income", "smoke")),
               createCPT.fromData(data, c("diabetes", "bmi", "smoke", "exercise")),
               createCPT.fromData(data, c("exercise", "income")),
               createCPT.fromData(data, "income"),
               createCPT.fromData(data, c("smoke", "income")),
               createCPT.fromData(data, c("stroke", "bmi", "bp", "cholesterol", "smoke", "exercise")))

# question (4) diabets
margvar = c('stroke','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q4(a) diabets, bad habits, ")
print(result)

margvar = c('stroke','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q4(a) diabets, good habits, ")
print(result)

margvar = c('stroke','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q4(b) diabets, poor health, ")
print(result)

margvar = c('stroke','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q4(b) diabets, good health, ")
print(result)

# question (4) stroke
margvar = c('diabetes','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q4(a) stroke, bad habits, ")
print(result)

margvar = c('diabetes','attack','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q4(a) stroke, good habits, ")
print(result)

margvar = c('diabetes','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q4(b) stroke, poor health, ")
print(result)

margvar = c('diabetes','attack','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q4(b) stroke, good health, ")
print(result)

# question (4) heart attack
margvar = c('stroke','diabetes','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q4(a) attack, bad habits, ")
print(result)

margvar = c('stroke','diabetes','angina', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q4(a) attack, good habits, ")
print(result)

margvar = c('stroke','diabetes','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q4(b) attack, poor health, ")
print(result)

margvar = c('stroke','diabetes','angina', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q4(b) attack, good health, ")
print(result)

# question (4) angina
margvar = c('stroke','diabetes','attack', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(1, 2))
print("Q4(a) angina, bad habits, ")
print(result)

margvar = c('stroke','diabetes','attack', 'income', 'bmi', 'bp', 'cholesterol')
result = infer(dataNet, margvar, c('smoke', 'exercise'), c(2, 1))
print("Q4(a) angina, good habits, ")
print(result)

margvar = c('stroke','diabetes','attack', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(1, 1, 3))
print("Q4(b) angina, poor health, ")
print(result)

margvar = c('stroke','diabetes','attack', 'income', 'smoke', 'exercise')
result = infer(dataNet, margvar, c('bp', 'cholesterol', 'bmi'), c(3, 2, 2))
print("Q4(b) angina, good health, ")
print(result)

####===================================================================================
# question 5 part one
diffvar = c('diabetes', 'stroke')
margvar = setdiff(var, diffvar)
result = infer(dataNet, margvar, c('diabetes'), c(1))
print("Q5 | Q4")
print(result)

result = infer(dataNet, margvar, c('diabetes'), c(3))
print("Q5 | Q4")
print(result)

# question 5 part two
dataNet = list(createCPT.fromData(data, c("angina", "bmi", "bp", "cholesterol", "smoke", "exercise")),
               createCPT.fromData(data, c("attack", "bmi", "bp", "cholesterol", "smoke", "exercise")),
               createCPT.fromData(data, c("bmi", "exercise", "income")),
               createCPT.fromData(data, c("bp", "exercise", "income", "smoke")),
               createCPT.fromData(data, c("cholesterol", "exercise", "income", "smoke")),
               createCPT.fromData(data, c("diabetes", "bmi", "smoke", "exercise")),
               createCPT.fromData(data, c("exercise", "income")),
               createCPT.fromData(data, "income"),
               createCPT.fromData(data, c("smoke", "income")),
               createCPT.fromData(data, c("stroke", "diabetes", "bmi", "bp", "cholesterol", "smoke", "exercise")))
result = infer(dataNet, margvar, c('diabetes'), c(1))
print("Q5 | New network")
print(result)

result = infer(dataNet, margvar, c('diabetes'), c(3))
print("Q5 | New network, ")
print(result)