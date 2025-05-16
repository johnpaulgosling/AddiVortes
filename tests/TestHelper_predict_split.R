library(AddiVortes)

set.seed(111333)
X <- matrix(rnorm(100),10,10)
Y <- rnorm(10)
X_test <- matrix(rnorm(100),10,10)
Y_test <- rnorm(10)

results <- AddiVortes(Y,X,10,
                      90,10,
                      6,0.85,3,0.8,3,25,
                      Y_test,X_test,
                      IntialSigma = "Linear")

#expect_equal(round(results[[2]][[1]],3), 0.724)
#expect_equal(round(results[[2]][[2]],3), 1.052)

preds <- PredictAddiVortes(results[[1]],
                           X_test,
                           Y_test)
preds[[1]]

plot(Y_test, preds[[2]])
