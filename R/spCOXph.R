# 
# library(simcausal)
# n <-  500
# D <- DAG.empty()
# 
# D <- D +
# node("W1", distr = "runif", min = -1, max = 1) +
#   node("W2", distr = "runif", min = -1, max = 1) +
#     node("A", distr = "rbinom", size = 1, prob = plogis(0.5*(W1 + W2  ) ))+
#   node("time", distr = "rexp",rate = exp(0.25*(-1+A )) )
#   
#   setD <- set.DAG(D ) 
# data <- sim(setD, n = n)
# data$status<-1
# data$id <- data$ID
# data <- data[,-1]
#  p<-poissonize(data, factors  = c("A", "W1", "W2"), compress = F)
# p
# X <- data.frame(  Y=out$data$`y..coxph`, bin ,  A= out$data$A,out$data[,paste0("X.",c("W1", "W2"))]  )
# X<- as.matrix(p[,c(3,4,5)])
# xt <- model.matrix(~ . -1, as.data.frame(p$interval))
# X <- cbind(xt,X)
# coef(cv.glmnet(X, p[,2],  offset = log(p[,6]),  family="poisson"), s = "lambda.min")
# coef(glm.fit(X, p[,2],  offset = log(p[,6]),  family=poisson()))
# coef(coxph(Surv(time, status) ~ ., data, ties='breslow'))
# 
# 
#  