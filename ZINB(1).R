require(ggplot2)
require(pscl)
require(MASS)
require(boot)
zinb <- read.csv("C:\\Users\\tanma\\Desktop\\ALL\\2nd same isi\\ISI project\\fish data.csv")

hist(zinb$count,
     breaks = 30,
     col = "skyblue",
     main = "Histogram of Fish Count",
     xlab = "Number of Fish Caught",
     ylab = "Frequency")
table(zinb$count == 0)
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})

summary(zinb)
#############
## histogram with x axis in log10 scale
ggplot(zinb, aes(count, fill = camper)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(camper ~ ., margins=TRUE, scales="free_y")
##########################
m1 <- zeroinfl(count ~ child + camper | persons,
               data = zinb, dist = "negbin")
summary(m1)
##########################################################################################
# Load necessary libraries
library(pscl)     # for zeroinfl
library(MASS)     # for glm.nb

######## Poisson Model
pois_model <- glm(count ~ child + camper, data = zinb, family = "poisson")
nb_model
###### Negative Binomial Model
nb_model <- glm.nb(count ~ child + camper, data = zinb)
nb_model
#####Zero-Inflated Poisson (ZIP)
zip_model <- zeroinfl(count ~ child + camper | persons, data = zinb, dist = "poisson")
zip_model
###### Zero-Inflated Negative Binomial (ZINB)
zinb_model <- zeroinfl(count ~ child + camper |persons, data = zinb, dist = "negbin")
zinb_model
# AIC comparison
AIC(pois_model, nb_model, zip_model, zinb_model)

# Vuong test: compare ZINB vs NB (non-nested models)
vuong(zinb_model, nb_model)
