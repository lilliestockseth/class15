# EBIO 338/538: Analysis and Visualization of Biological Data
# Class 15: Poisson & Binomial models

# Plotting data and functions
plot(cars$dist ~ cars$speed)


# Simple linear model
plot(cars$dist ~ cars$speed)
m <- lm(cars$dist ~ cars$speed)
k <- coef(m) # Save coefficients as object called k
curve(k[1] + k[2] * x, 
      from = min(cars$speed),
      to = max(cars$speed), 
      add = TRUE)


# Exponential model

library(bbmle)
plot(cars$dist ~ cars$speed)
m <- lm(cars$dist ~ cars$speed)
k <- coef(m)
curve(k[1] + k[2] * x, 
      from = min(cars$speed),
      to = max(cars$speed), 
      add = TRUE)
m2 <- mle2(cars$dist ~ dnorm(mean = exp(a + b * cars$speed), 
                             sd = exp(s)),
           start = list(a = 0, b = 0, s = 2), data = cars)
k <- coef(m2)
curve(exp(k["a"] + k["b"] * x),
      from = min(cars$speed),
      to = max(cars$speed), 
      col = "orange",
      lwd = 2,
      add = TRUE)



# Eagles example

# Two options for how NOT to analyze:
# 1) Convert to proportions

library(MASS)
data(eagles)
eagles
?eagles
hist(eagles$y/eagles$n,
     col = "orange", xlab = "proportion successful", 
     main = "")

# What are a couple reasons why a proprotion is a bad outcome variable
# to analyze with a lm?


#2) Model y as normally distributed outcome
hist(eagles$y, col = "orange", 
     xlab = "number of successes",
     main = "", breaks = 30)

# What are the problems with just analyzing the successes?



# Binomial distribution example figures matching the ones from the lecture
# 

# Make a simple utility function to help plotting (this is basically just a histogram)
prep.bar <- function( d , dmax=max(d) ) {
  x <- rep(0,dmax+1)
  for ( i in 1:(dmax+1) ) {
    x[i] <- length( d[d==(i-1)] )
  }
  x
}

# Show examples random samples from a binomial distribution
# First where there's just one draw and a probability of success of 0.1
y <- rbinom(1000, prob = 0.1, size = 1)
barplot(prep.bar(y), names.arg = 0:max(y), 
        xlab= "y", ylab = "Frequency", col = "orange")

# Second with a probability of success of 0.7
y <- rbinom(1000, prob = 0.7, size = 1)
barplot(prep.bar(y), names.arg = 0:max(y), 
        xlab= "y", ylab = "Frequency", col = "orange")


# Alter prep.bar so both plot on same scale
prep.bar2 <- function( d , dmax=max(d) ) {
  x <- rep(0,dmax+1)
  for ( i in 1:(dmax+4) ) {
    x[i] <- length( d[d==(i-1)] )
  }
  x
}

# Now using those same two probability values but with 10 trials each
y <- rbinom(10000, prob = 0.1, size = 10)
barplot(prep.bar2(y), names.arg = 0 : max(y+3), 
        xlab= "y", ylab = "Frequency", col = "orange")


y <- rbinom(10000, prob = 0.7, size = 10)
barplot(prep.bar(y), names.arg = 0 : max(y), 
        xlab= "y", ylab = "Frequency", col = "orange")


# Show likelihood imposed over barplot
ss <- 10
set.seed(4)
y <- rbinom(ss, prob = 0.1, size = 10)
barplot(prep.bar(y), names.arg = 0 : max(y),
        xlab = "y", ylab = "Frequency", space = 0, col = "orange")
p <- ss*dbinom(0 : max(y), prob = 0.1, size = 10)
points(c(0 : max(y) + 0.5), p, pch = 16, cex = 2)
# These points show the theoretical expectation for how many samples should
# have 0, 1, 2... successes out of 10 given the probability parameter of the
# binomial distribution equal to 0.1



# Now let's think about how we'd go about modeling probability using
# our typical Y = intercept + slope * X formula
curve(0.1 + 0.2 * x, from = 0, to = 3, 
      ylim = c(0, 1.2), ylab = "probability")
curve(0.1 + 1 * x, from = 0, to = 3, ylim = c(0, 1), 
      ylab = "probability", add = TRUE)
lines(c(-1, 4), c(1, 1), lty = 2)
rect(-1, 1, 4, 2, col = "orange", density = 10, angle = -45)
# Does this make sense? Nope! Can't have a probability greater than 1


# But it will work if we use a logit link. Note that inside this,
# we have the same Y = intercept + slope * X formula
curve(1 / (1 + exp(1 - 0.2 * x)), from = 0, to = 3, 
      ylim=c(0,1.2) , ylab="probability")
curve(1 / (1 + exp(1 - 1 * x)), from = 0 , to = 3, 
      ylim = c(0, 1), ylab = "probability", add = TRUE)
curve(1 / (1 + exp(1 - 2 * x)), from = 0, to = 3,
      ylim = c(0, 1), ylab = "probability", add = TRUE)
curve(1 / (1 + exp(1 - 6 * x)), from = 0, to = 3, 
      ylim = c(0, 1), ylab = "probability", add = TRUE)
rect(-1, 1, 4, 2, col = "orange", density = 10, angle = -45)
lines(c(-1, 4), c(1, 1), lty = 2)
# Using the logit link function, we never exceed 1 (great!) and 
# we still get to use the linear model structure internally.


# Another way to understand the relationship between the probability
# scale and the logit scale:

# First, get a uniform distribution of probability values between 0 and 1
rand_prob <- runif(1000, 0, 1)

# Then logit transform these to get them on the probability scale
logit_rand_prob <- log(rand_prob / (1-rand_prob))
# Note that these values are no longer just between 0 and 1. In fact, they
# can go from negative infinity to infinity (just like the normal distribution!)

# Plot the relationship between the values on the probability scale
# and the logit scale
plot(x = rand_prob, 
     y = logit_rand_prob,
     pch = 16,
     col = "orange",
     xlab = "Probability scale",
     ylab = "Logit scale")


# To go back from the logit scale to the probability scale, we use the
# inverse logit on the logit scaled values.
invlogit_logit_rand_prob <- 1 / (1 + exp(-logit_rand_prob))
plot(x = rand_prob, 
     y = invlogit_logit_rand_prob,
     pch = 16,
     col = "orange",
     xlab = "Probability scale",
     ylab = "Back to the orobability scale again!")





# Analyzing the eagle data (the right way)

library(MASS)
data(eagles)
d <-  eagles

#make dummy variable for pirate
d$pirate <- ifelse(d$P == "L", 1, 0)
d

# intercept only model
m0 <- mle2(d$y ~ dbinom(size = d$n, 
                        prob = 1/(1+exp(a))),
           start = list(a = 0), data = d)

#size of pirating eagle model
mP <- mle2(d$y ~ dbinom(size=d$n, prob = 1/(1+exp(a + bp * d$pirate))),
           start = list(a = 0, bp = 0), data = d)

AICtab(m0, mP, base = TRUE, weights = TRUE)

# Which is the better performing model?




# Note that we can use glm (similar to lm) for a version that's simpler to code
?glm
# 
# A couple things are different because we're using a binomial distribution
# for the outcome variable rather than the normal distribution.

# Here cbind(y, n - y) gives the number of successes (first column) and 
# number of failures (second column)
cbind(d$y, d$n - d$y)

mP_glm <- glm(cbind(d$y, d$n - d$y) ~ d$pirate,
              family = "binomial")

# How do the coefficient estimates compare?
summary(mP_glm)
summary(mP)






# Salamander data
d <- read.csv("salamanders.csv")

barplot(prep.bar(d$SALAMAN), names.arg = 0 : max(d$SALAMAN), 
        xlab = "salamanders", ylab = "Frequency", col = "orange")



# Poisson examples

y <- rpois(1000, lambda = 10)
barplot(prep.bar(y), names.arg = 0 : max(y),
        xlab = "y", ylab = "Frequency", col = "orange" )


y <- rpois(1000, lambda = 0.1)
barplot(prep.bar(y), names.arg = 0 : max(y),
        xlab = "y", ylab = "Frequency", col = "orange" )



#Getting started on salamander example

# intercept model
m <- mle2(d$SALAMAN ~ dpois(a), 
          start = list(a = mean(d$SALAMAN)), data = d)
# linear cover model
mc <- mle2(d$SALAMAN ~ dpois(a + bc * d$PCTCOVER), 
           start = list(a = mean(d$SALAMAN), bc = 0), data = d)

# Compare the models using AIC
AICtab(m, mc, base = TRUE, weights = TRUE)
# Which models perform better

# How would you fit the mc model using glm() rather than mle2()?
# Hint: family = "poisson"

