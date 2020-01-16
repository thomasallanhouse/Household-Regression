# Load magittr library for data piping and dplyr for tidy data
library("magrittr")
library("dplyr")
library("stats4")

# Load the data
plosdata <- readr::read_csv("plosData_300519.csv")

# Drop HH variable
select(plosdata,-HH) -> plosdataClean

# Build multivariate logistic regression model
#contactageP = glm(secondary_case_with_ILI ~ AgeContact0_4 + AgeContact18_49 + AgeContact5_17 + AgeIndex0_4 + AgeIndex18_49 + AgeIndex5_17 + Chronic + HHmember2_3 + HHmember4_5 + ResidenceUrban + Season2009 + Gender + Shared_bedroom_index + Vaccinated, data = plosdata, family = binomial(link = cloglog)) # Binomial implies logistic
contactageP = glm(secondary_case_with_ILI ~ ., data = plosdataClean, family = binomial(link = cloglog)) # Binomial implies logistic

ctable <- coef(summary(contactageP))
#ORcontactageP <- cbind(exp(cbind("OR"=coef(contactageP),confint(contactageP,level = 0.95))),"p-value"=ctable[,"Pr(>|z|)"])
CIs <- cbind(coef(contactageP),confint(contactageP,level = 0.95),"p-value"=ctable[,"Pr(>|z|)"])
View(CIs)
aic = AIC(contactageP)