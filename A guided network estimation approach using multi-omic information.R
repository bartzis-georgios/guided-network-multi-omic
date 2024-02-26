# Empty environment ----
rm(list = ls())
# If necessary, change the working directory to the one containing the data ----
workingDir<-"XXX"
setwd(workingDir)

# Load relevant packages ----
library(glmgraph)
library(data.table)
library(huge)
library(LassoNet)

# Should the glmgraph run once or multiple times based on random subset of data? ----
runOnce <- T
# Should the LassoNet Laplacian be used or glmgraph transformation? ----
LassoNetTransformation <- T
# Is the guiding data network structrure given or estimated ----
givenGuidingStructure <- F
# Read target data ----
targetData <- fread("./metDataDT.csv")
# Read guiding data ---- 
guidingData <- fread("./genDataDT.csv")
if(givenGuidingStructure){
  # Read guiding data structure ----
  guidingDataStructure <- read.csv("./givenguidingstructure.csv", header = TRUE,row.names=1, dec=".",sep=",")
}else{
  guidingDataStructure <- NULL
}


# Keep only common samples between datasets ----
commonSamples <- intersect(targetData$sampleName,guidingData$sampleName)
targetData <- targetData[sampleName %chin% commonSamples]
guidingData <- guidingData[sampleName %chin% commonSamples]
sampleNameVariable <- "sampleName"
targetDatavars <- setdiff(names(targetData),sampleNameVariable)
guidingDataVars <- setdiff(names(guidingData),sampleNameVariable)


if(givenGuidingStructure == F){
  # Estimate the network structure of the guiding data using Graphical Lasso (GL) ----
  ## GL regularization parameter sequence ----
  regParSequence <- seq(from=0.99, to=0.01, by= -0.01)
  ## Sequence of GL ----
  GLsequence <- huge(
    as.matrix(
      guidingData[,-c(sampleNameVariable),with = F]
    ),
    method="glasso",
    lambda = regParSequence
  )
  ## Select optimal GL based on StARS ----
  GLselection <- huge.select(
    est = GLsequence,
    criterion = "stars",
    rep.num = 100,
    stars.thresh = 0.001
  )
  ## Optimal icov matrix ----
  optIcov <- GLselection$opt.icov
  colnames(optIcov) <- rownames(optIcov) <- guidingDataVars
}else{
  optIcov <- guidingDataStructure
}

## Create Laplacian of icov matrix ----
absOptIcov <- abs(optIcov)
if(LassoNetTransformation){
  ### Laplacian based on Weber 2023 (LassoNet) ----
  L <- mat.to.laplacian(absOptIcov)
}else{
  ### Laplacian based on Li & Li package suggestion (glmgraph) ----
  diag(absOptIcov) <- 0
  diagL <- apply(absOptIcov,1,sum)
  L <- -absOptIcov
  diag(L) <- diagL
}

## Function for running penalty selection in glmgraph using cv ----
fit.lm <- function(target, guiding, Laplacian){
  cv.obj <-
    cv.glmgraph(
      X = guiding,
      Y = target,
      L = Laplacian,
      nfolds = 5,
      family = "gaussian",
      type.measure = "mse",
      standardize = TRUE,
      penalty = "lasso",
    )
  result <- data.table(
    lambda1 = cv.obj$lambda1.min, 
    lambda2 = cv.obj$lambda2.min
  )
  return(result)
}

if(runOnce){
  ## Extract L1 and L2 for all variables in target data ----
  l1l2s <- lapply(targetDatavars,FUN = function(i){
    target <- as.matrix(targetData[,i,with = F])
    guiding <- as.matrix(guidingData[,-c(sampleNameVariable),with = F])
    Laplacian <- L
    fit.lm(target,guiding,Laplacian)
  })
  names(l1l2s) <- targetDatavars
  l1l2s <- rbindlist(l = l1l2s,use.names = T,idcol = "TargetData")
}else{
  ## Run glmgraph using cv multiple times on random subset of data and keep avg l1, l2
  times=5
  pTrain <- (10 * sqrt(N)/N)*N
  trainingSamples <- lapply(1:times, FUN= function(i){
    sample(N,pTrain)
  })
  trainingSamples <- do.call(cbind, trainingSamples)
  testSamples <- lapply(1:ncol(trainingSamples),FUN = function(i,x = trainingSamples){
    c(1:N)[-x[,i]]
  })
  testSamples <- do.call(cbind, testSamples)
  l1l2s <- lapply(1:times,FUN = function(j){
    l1l2sIteration <- lapply(targetDatavars,FUN = function(i,iteration = j){
      target <- as.matrix(targetData[trainingSamples[,iteration],i,with = F])
      guiding <- as.matrix(guidingData[trainingSamples[,iteration],-c(sampleNameVariable),with = F])
      Laplacian <- L
      fit.lm(target,guiding,Laplacian)
    })
    names(l1l2sIteration) <- targetDatavars
    l1l2sIteration <- as.matrix(rbindlist(l = l1l2sIteration,use.names = T))
    rownames(l1l2sIteration) <- targetDatavars
    return(l1l2sIteration)
  })
  l1l2s <- Reduce("+", l1l2s) / length(l1l2s)
  l1l2s <- as.data.table(l1l2s,keep.rownames = T)
  setnames(l1l2s,"rn","TargetData")
}

## Extract target data predicted values based on selected penalties ----
guiding <- as.matrix(guidingData[,-c(sampleNameVariable),with = F])
predictedValuesAndCoefficients <- lapply(targetDatavars,FUN = function(i){
  Y <- as.matrix(targetData[,i, with = F])
  fit <- glmgraph(
    X = guiding,
    Y = Y,
    L = L,
    family="gaussian",
    # type.measure = "mse",
    standardize = TRUE,
    penalty = "lasso",
    lambda1 = l1l2s[TargetData == i]$lambda1,
    lambda2 = l1l2s[TargetData == i]$lambda2
  )
  coef(fit)
  return(
    list(
      predValues = predict(fit,guiding)[[1]],
      coefficients = coef(fit)[[1]][-1,]
    )
  )
})
names(predictedValuesAndCoefficients) <- targetDatavars
predictedValues <- as.data.table(sapply(predictedValuesAndCoefficients,"[[",1))
modelCoefficients <- sapply(predictedValuesAndCoefficients,"[[",2)
  
# Build network of original target data ----
## Sequence of GL ----
regParSequence <- sequ<-seq(from=1,to=0.001, by=-0.0005)
GLSequenceTargetORIG <- huge(
  x = as.matrix(targetData[,-c(sampleNameVariable), with = F]),
  method="glasso",
  lambda = regParSequence
)
## Select optimal GL based on StARS for original target data ----
GLselectionTargetORIG <- huge.select(
  est = GLSequenceTargetORIG,
  criterion ="stars",
  rep.num = 100,
  stars.thresh = 0.05
)
huge.plot(GLselectionTargetORIG$opt.icov)


# Build network of target data conditional on guiding ----
## Sequence of GL ----
GLSequenceTargetNEW <- huge(
  x = as.matrix(predictedValues),
  method="glasso",
  lambda = regParSequence
)
sparsityOfOrig <- GLselectionTargetORIG$opt.sparsity*(ncol(predictedValues)-1)*ncol(predictedValues)/2
## Select network based on sparsity of original data (so that networks can be comparable) ----
optIndexNEW <- max(which(GLSequenceTargetNEW$sparsity*(ncol(predictedValues)-1)*ncol(predictedValues)/2 <= sparsityOfOrig))
huge.plot(GLSequenceTargetNEW$icov[[optIndexNEW]])
