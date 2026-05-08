benchmark_datasets<-list()

### Regression datasets
{
  # 1. Longley Economic Data (datasets) - 16 obs, 6 covariates (Regression)        x SVM 0.4187, flexBART 0.7534,AddiVortes 0.9
  # Stress Test: Extreme micro-sample with high collinearity.
  data(longley, package = "datasets")
  longley <- na.omit(longley)
  benchmark_datasets$longley <- list(
    X = model.matrix(Employed ~ . - 1, data = longley),
    Y = as.numeric(longley$Employed)
  )
  rm(longley)
  
  # 2. Freeny's Revenue Data (datasets) - 39 obs, 4 covariates (Regression)        x  SVM 0.01452, RF 0.01820, AddiVortes 0.0191
  # Stress Test: Very small N.
  data(freeny, package = "datasets")
  freeny <- na.omit(freeny)
  benchmark_datasets$freeny <- list(
    X = model.matrix(y ~ . - 1, data = freeny),
    Y = as.numeric(freeny$y)
  )
  rm(freeny)
  
  # 3. Cabbages (MASS) - 60 obs, 3 covariates (Regression)                          x  ALL
  data(cabbages, package = "MASS")
  cabbages <- na.omit(cabbages)
  benchmark_datasets$cabbages <- list(
    X = model.matrix(HeadWt ~ . - 1, data = cabbages),
    Y = as.numeric(cabbages$HeadWt)
  )
  rm(cabbages)
  
  # 4. Anatomical Data from Cats (MASS) - 144 obs, 2 covariates (Regression)       x SVM, flexBART, wBART
  # Target: Heart weight (Hwt)
  data(cats, package = "MASS")
  cats <- na.omit(cats)
  benchmark_datasets$cats <- list(
    X = model.matrix(Hwt ~ . - 1, data = cats),
    Y = as.numeric(cats$Hwt)
  )
  rm(cats)
  
  # 5. Chemical Manufacturing Process (AppliedPredictiveModeling) - 176 obs, 57 covariates     x
  data(ChemicalManufacturingProcess, package = "AppliedPredictiveModeling")
  ChemicalManufacturingProcess<-na.omit(ChemicalManufacturingProcess)
  benchmark_datasets$ChemManuf <- list(
    X = model.matrix(Yield ~ . - 1, data = ChemicalManufacturingProcess),
    Y = as.numeric(ChemicalManufacturingProcess$Yield)
  )
  rm(ChemicalManufacturingProcess)
  
  # 6. Road Casualties (datasets) - 192 obs, 7 covariates (Regression)
  data(Seatbelts, package = "datasets")
  Seatbelts_df <- na.omit(as.data.frame(Seatbelts))
  benchmark_datasets$Seatbelts <- list(
    X = model.matrix(DriversKilled ~ . - 1, data = Seatbelts_df),
    Y = as.numeric(Seatbelts_df$DriversKilled)
  )
  rm(Seatbelts, Seatbelts_df)
  
  # 7. Blood-Brain Barrier (caret) - 208 obs, 134 covariates
  # Using caret for more diverse benchmark options
  data(BloodBrain, package = "caret")
  benchmark_datasets$BloodBrain <- list(
    X = as.matrix(bbbDescr),
    Y = as.numeric(logBBB)
  )
  rm(bbbDescr, logBBB)
  
  # 8. Computer CPU Performance (MASS) - 209 obs, 7 covariates (Regression)
  # Target: Relative CPU performance (perf)
  data(cpus, package = "MASS")
  cpus <- na.omit(cpus)
  benchmark_datasets$cpus <- list(
    X = model.matrix(perf ~ . - 1, data = cpus),
    Y = as.numeric(cpus$perf)
  )
  rm(cpus)
  
  # 9. Tecator (caret) - 215 obs, 100 covariates
  data(tecator, package = "caret")
  benchmark_datasets$Tecator <- list(
    X = as.matrix(absorp),
    Y = as.numeric(endpoints[, 2])
  )
  rm(absorp, endpoints)
  
  # 10. Vehicle Fuel Economy (ggplot2) - 234 obs, 10 covariates (Regression)
  # Target: Highway miles per gallon (hwy)
  data(mpg, package = "ggplot2")
  mpg <- na.omit(mpg)
  benchmark_datasets$mpg <- list(
    X = model.matrix(hwy ~ . - 1, data = mpg),
    Y = as.numeric(mpg$hwy)
  )
  rm(mpg)
  
  # 11. Cox2 Activity (caret) - 266 obs, 255 covariates (Regression)
  # Stress Test: Ultra-High Dimensionality (P ≈ N). Perfect test for Dirichlet Sparsity.
  data(cox2, package = "caret")
  cox2_data <- na.omit(data.frame(cox2IC50, cox2Descr))
  benchmark_datasets$Cox2 <- list(
    X = model.matrix(cox2IC50 ~ . - 1, data = cox2_data),
    Y = as.numeric(cox2_data$cox2IC50)
  )
  rm(cox2Descr, cox2IC50, cox2Class, cox2_data)
  
  # 12. Baseball Hitters (ISLR2) - 322 obs, 19 covariates (Regression)
  # Target: Player Salary
  data(Hitters, package = "ISLR2")
  Hitters <- na.omit(Hitters)
  benchmark_datasets$Hitters <- list(
    X = model.matrix(Salary ~ . - 1, data = Hitters),
    Y = as.numeric(Hitters$Salary)
  )
  rm(Hitters)
  
  # 13. Ozone Readings (mlbench) - 366 obs, 12 covariates (Regression)
  # Target: Daily maximum one-hour-average ozone reading (V4)
  data(Ozone, package = "mlbench")
  Ozone <- na.omit(Ozone)
  benchmark_datasets$Ozone <- list(
    X = model.matrix(V4 ~ . - 1, data = Ozone),
    Y = as.numeric(Ozone$V4)
  )
  rm(Ozone)
  
  # 14. Auto MPG (ISLR2) - 392 obs, 8 covariates (Regression)
  # Target: Miles Per Gallon (mpg)
  data(Auto, package = "ISLR2")
  Auto <- na.omit(Auto)
  benchmark_datasets$Auto <- list(
    X = model.matrix(mpg ~ . - 1, data = Auto),
    Y = as.numeric(Auto$mpg)
  )
  rm(Auto)
  
  # 15. Carseats Sales (ISLR2) - 400 obs, 10 covariates (Regression)
  # Target: Unit Sales at each location
  data(Carseats, package = "ISLR2")
  Carseats <- na.omit(Carseats)
  benchmark_datasets$Carseats <- list(
    X = model.matrix(Sales ~ . - 1, data = Carseats),
    Y = as.numeric(Carseats$Sales)
  )
  rm(Carseats)
  
  # 16. Credit Card Balance (ISLR2) - 400 obs, 11 covariates (Regression)
  # Target: Average credit card balance
  data(Credit, package = "ISLR2")
  Credit <- na.omit(Credit)
  benchmark_datasets$Credit <- list(
    X = model.matrix(Balance ~ . - 1, data = Credit),
    Y = as.numeric(Credit$Balance)
  )
  rm(Credit)
  
  # 17. Chick Weight (datasets) - 578 obs, 3 covariates (Regression)
  data(ChickWeight, package = "datasets")
  ChickWeight <- na.omit(ChickWeight)
  benchmark_datasets$ChickWeight <- list(
    X = model.matrix(weight ~ . - 1, data = ChickWeight),
    Y = as.numeric(ChickWeight$weight)
  )
  rm(ChickWeight)
  
  # 18. SoyBean (mlbench) - 683 obs, 35 covariates ## ordinal classification ###
  data(Soybean, package = "mlbench")
  Soybean<-na.omit(Soybean)
  benchmark_datasets$Soybean <- list(
    X = model.matrix(Class ~ . - 1, data = Soybean),
    Y = as.numeric(Soybean$Class)
  )
  rm(Soybean)
  
  # 19. College (ISLR2) - 777 obs, 17 covariates (Regression)
  # Target: Number of applications received
  data(College, package = "ISLR2")
  College <- na.omit(College)
  benchmark_datasets$College <- list(
    X = model.matrix(Apps ~ . - 1, data = College),
    Y = as.numeric(College$Apps)
  )
  rm(College)
  
  # 20. Cars 2004 (caret) - 804 obs, 17 covariates (Regression)
  data(cars, package = "caret")
  cars <- na.omit(cars)
  benchmark_datasets$Cars2004 <- list(
    X = model.matrix(Price ~ . - 1, data = cars),
    Y = as.numeric(cars$Price)
  )
  rm(cars)
  
  # 21. Sacramento Real Estate (caret) - 932 obs, 8 covariates (Regression)
  data(Sacramento, package = "caret")
  Sacramento <- na.omit(Sacramento)
  benchmark_datasets$Sacramento <- list(
    X = model.matrix(price ~ . - 1, data = Sacramento),
    Y = as.numeric(Sacramento$price)
  )
  rm(Sacramento)
  
  # 22. Solubility (AppliedPredictiveModeling) - 951 obs, 228 covariates (Regression)
  # Excellent for testing massive P dimensionality
  data(solubility, package = "AppliedPredictiveModeling")
  solubility_df <- na.omit(data.frame(solubility = solTrainY, solTrainX))
  benchmark_datasets$Solubility <- list(
    X = model.matrix(solubility ~ . - 1, data = solubility_df),
    Y = as.numeric(solubility_df$solubility)
  )
  rm(solTrainX, solTrainY, solTestX, solTestY, solubility_df)
  
  # 23. Fiji Earthquakes (datasets) - 1000 obs, 4 covariates (Regression)
  # Target: Earthquake Magnitude (mag)
  data(quakes, package = "datasets")
  quakes <- na.omit(quakes)
  benchmark_datasets$quakes <- list(
    X = model.matrix(mag ~ . - 1, data = quakes),
    Y = as.numeric(quakes$mag)
  )
  rm(quakes)
  
  # 24. Communities and Crime (fairml) (observations 1968, covariates 102)
  # Note: Ensure install.packages("fairml") is run if not present
  data(communities.and.crime, package = "fairml")
  communities.and.crime<-communities.and.crime[,-2]
  communities.and.crime<-na.omit(communities.and.crime)
  benchmark_datasets$CommunitiesCrime <- list(
    X = model.matrix(ViolentCrimesPerPop ~ . - 1, data = communities.and.crime),
    Y = as.numeric(communities.and.crime$ViolentCrimesPerPop)
  )
  rm(communities.and.crime)
  
  # 25. Dutch Schools Data (MASS) - 2287 obs, 5 covariates (Regression)
  # Stress Test: Large N, Low P scalability test.
  data(nlschools, package = "MASS")
  nlschools <- na.omit(nlschools)
  benchmark_datasets$nlschools <- list(
    X = model.matrix(lang ~ . - 1, data = nlschools),
    Y = as.numeric(nlschools$lang)
  )
  rm(nlschools)
  
  # 26. Ames Housing (AmesHousing) - 2930 obs, 74 covariates (Regression)
  # Note: make_ames() cleans the raw data automatically for immediate ML use.
  # Ensure install.packages("AmesHousing") is run if not present.
  ames <- AmesHousing::make_ames()
  ames <- na.omit(ames)
  benchmark_datasets$AmesHousing <- list(
    X = model.matrix(Sale_Price ~ . - 1, data = ames),
    Y = as.numeric(ames$Sale_Price)
  )
  rm(ames)
  
  # 27. Wage (ISLR2) - 3000 obs, 10 covariates (Regression)
  data(Wage, package = "ISLR2")
  Wage <- na.omit(Wage)
  benchmark_datasets$Wage <- list(
    X = model.matrix(wage ~ . - 1, data = Wage),
    Y = as.numeric(Wage$wage)
  )
  rm(Wage)
  
  # 28. DNA (observations 3186, Covariates 181) #### Ordinal Classification ###
  data(DNA, package = "mlbench")
  benchmark_datasets$DNA <- list(
    X = model.matrix(Class ~ . - 1, data = DNA),
    Y = as.numeric(DNA$Class)
  )
  rm(DNA)
  
  # 29. Abalone (AppliedPredictiveModeling) - 4177 obs, 8 covariates (Regression) ## ordinal classification ###
  data(abalone, package = "AppliedPredictiveModeling")
  abalone <- na.omit(abalone)
  benchmark_datasets$Abalone <- list(
    X = model.matrix(Rings ~ . - 1, data = abalone),
    Y = as.numeric(abalone$Rings)
  )
  rm(abalone)
  
  # 30. Satellite (observations 6235, Covariates 36) #### Ordinal Classification ###
  data(Satellite, package = "mlbench")
  benchmark_datasets$Satellite <- list(
    X = model.matrix(classes ~ . - 1, data = Satellite),
    Y = as.numeric(Satellite$classes)
  )
  rm(Satellite)
  
  # 31. Texas Housing (ggplot2) - 8602 obs, 8 covariates (Regression)
  # Target: Median sale price
  data(txhousing, package = "ggplot2")
  txhousing <- na.omit(txhousing)
  benchmark_datasets$TexasHousing <- list(
    X = model.matrix(median ~ . - 1, data = txhousing),
    Y = as.numeric(txhousing$median)
  )
  rm(txhousing)
  
  # 32. Bike Sharing (ISLR2) - 8645 obs, 14 covariates (Regression)
  data(Bikeshare, package = "ISLR2")
  Bikeshare <- na.omit(Bikeshare)
  benchmark_datasets$Bikeshare <- list(
    X = model.matrix(bikers ~ . - 1, data = Bikeshare),
    Y = as.numeric(Bikeshare$bikers)
  )
  rm(Bikeshare)
  
  # 33. Diamonds (ggplot2) - 53940 obs, 9 covariates (Regression)
  data(diamonds, package = "ggplot2")
  diamonds <- na.omit(diamonds)
  benchmark_datasets$Diamonds <- list(
    X = model.matrix(price ~ . - 1, data = diamonds),
    Y = as.numeric(diamonds$price)
  )
  rm(diamonds)
  
  # Final memory reclamation
  gc()
  }

### Binary Datasets
{
  
  # 1. Urine Crystals (boot) - 79 obs (after NA drop), 6 covariates
  # Target: Presence of calcium oxalate crystals (r).
  data(urine, package = "boot")
  urine <- na.omit(urine)
  unique_vals <- sort(unique(urine$r))
  benchmark_datasets$urine <- list(
    X = model.matrix(r ~ . - 1, data = urine),
    Y = as.numeric(urine$r == unique_vals[2])
  )
  rm(urine, unique_vals)
  
  # 2. Prostate Cancer Progression (rpart) - 134 obs (after NA drop), 7 covariates
  # Target: Tumor progression status (pgstat).
  data(stagec, package = "rpart")
  stagec <- na.omit(stagec)
  unique_vals <- sort(unique(stagec$pgstat))
  benchmark_datasets$stagec <- list(
    X = model.matrix(pgstat ~ . - 1, data = stagec),
    Y = as.numeric(stagec$pgstat == unique_vals[2])
  )
  rm(stagec, unique_vals)
  
  # 3. Quine School Attendance (MASS) - 146 obs, 4 covariates
  # Target: Learner status (Lrn) - Average Learner (AL) vs Slow Learner (SL).
  data(quine, package = "MASS")
  quine <- na.omit(quine)
  unique_vals <- sort(unique(quine$Lrn))
  benchmark_datasets$quine <- list(
    X = model.matrix(Lrn ~ . - 1, data = quine),
    Y = as.numeric(quine$Lrn == unique_vals[2])
  )
  rm(quine, unique_vals)
  
  # 4. Student Survey (MASS) - 168 obs (after NA drop), 11 covariates
  # Target: Sex of the student (Male vs Female).
  data(survey, package = "MASS")
  survey <- na.omit(survey)
  unique_vals <- sort(unique(survey$Sex))
  benchmark_datasets$survey <- list(
    X = model.matrix(Sex ~ . - 1, data = survey),
    Y = as.numeric(survey$Sex == unique_vals[2])
  )
  rm(survey, unique_vals)
  
  # 5. Infant Birth Weight (MASS) - 189 obs, 8 covariates (Binary Classification)
  # Target: Low birth weight (1) vs Normal (0). (Note: 'bwt' col removed as it is the target in grams)
  data(birthwt, package = "MASS")
  birthwt <- na.omit(birthwt)
  birthwt <- birthwt[, names(birthwt) != "bwt"] # Drop the continuous version of the target
  benchmark_datasets$birthwt <- list(
    X = model.matrix(low ~ . - 1, data = birthwt),
    Y = as.numeric(birthwt$low)
  )
  rm(birthwt)
  
  # 6. Crabs Species (MASS) - 200 obs, 6 covariates
  # Target: Crab species (sp) - Blue (B) vs Orange (O).
  data(crabs, package = "MASS")
  crabs <- na.omit(crabs)
  unique_vals <- sort(unique(crabs$sp))
  benchmark_datasets$crabs <- list(
    X = model.matrix(sp ~ . - 1, data = crabs),
    Y = as.numeric(crabs$sp == unique_vals[2])
  )
  rm(crabs, unique_vals)
  
  # 7. Melanoma Ulceration (boot) - 205 obs, 6 covariates
  # Target: Presence of ulceration (ulcer).
  data(melanoma, package = "boot")
  melanoma <- na.omit(melanoma)
  unique_vals <- sort(unique(melanoma$ulcer))
  benchmark_datasets$melanoma <- list(
    X = model.matrix(ulcer ~ . - 1, data = melanoma),
    Y = as.numeric(melanoma$ulcer == unique_vals[2])
  )
  rm(melanoma, unique_vals)
  
  # 8. Sonar
  data(Sonar, package = "mlbench")
  benchmark_datasets$Sonar <- list(
    X = model.matrix(Class ~ . - 1, data = Sonar),
    Y = as.numeric(Sonar$Class)
  )
  rm(Sonar)
  
  # 9. Congressional Voting Records 1984 (mlbench) - 232 obs, 16 covariates (Binary Classification)
  # Target: Democrat (0) vs Republican (1). Dropped NAs reduce it from 435 to 232 obs.
  data(HouseVotes84, package = "mlbench")
  HouseVotes84 <- na.omit(HouseVotes84)
  benchmark_datasets$HouseVotes84 <- list(
    X = model.matrix(Class ~ . - 1, data = HouseVotes84),
    Y = as.numeric(HouseVotes84$Class == "republican")
  )
  rm(HouseVotes84)
  
  # 10. Presence of Bacteria (MASS) - 232 obs, 5 covariates (Binary Classification)
  # Target: Presence of H. influenzae (y/n).
  data(bacteria, package = "MASS")
  bacteria <- na.omit(bacteria)
  benchmark_datasets$bacteria <- list(
    X = model.matrix(y ~ . - 1, data = bacteria),
    Y = as.numeric(bacteria$y == "y")
  )
  rm(bacteria)
  
  # 11. Infertility (datasets) - 248 obs, 8 covariates
  # Stress Test: Highly stratified medical data.
  # Target: case indicator (0 vs 1).
  data(infert, package = "datasets")
  infert <- na.omit(infert)
  unique_vals <- sort(unique(infert$case))
  benchmark_datasets$infert <- list(
    X = model.matrix(case ~ . - 1, data = infert),
    Y = as.numeric(infert$case == unique_vals[2])
  )
  rm(infert, unique_vals)
  
  # 12. Space Shuttle Autolander (MASS) - 256 obs, 6 covariates
  # Stress Test: Strictly categorical predictor space (converts to dummies).
  # Target: landing type (use) - auto vs noauto.
  data(shuttle, package = "MASS")
  shuttle <- na.omit(shuttle)
  unique_vals <- sort(unique(shuttle$use))
  benchmark_datasets$shuttle <- list(
    X = model.matrix(use ~ . - 1, data = shuttle),
    Y = as.numeric(shuttle$use == unique_vals[2])
  )
  rm(shuttle, unique_vals)
  
  # 13. DHFR Inhibitors (caret) - 325 obs, 228 covariates
  # Stress Test: Ultra-high dimensionality where P is almost equal to N.
  # Target: Biological activity (Y).
  data(dhfr, package = "caret")
  dhfr <- na.omit(dhfr)
  unique_vals <- sort(unique(dhfr$Y))
  benchmark_datasets$DHFR <- list(
    X = model.matrix(Y ~ . - 1, data = dhfr),
    Y = as.numeric(dhfr$Y == unique_vals[2])
  )
  rm(dhfr, unique_vals)
  
  # 14. Boston River Bound (MASS) - 506 obs, 13 covariates
  # Stress Test: Repurposed standard regression task into binary classification.
  # Target: Charles River dummy variable (chas).
  data(Boston, package = "MASS")
  Boston <- na.omit(Boston)
  unique_vals <- sort(unique(Boston$chas))
  benchmark_datasets$Boston_chas <- list(
    X = model.matrix(chas ~ . - 1, data = Boston),
    Y = as.numeric(Boston$chas == unique_vals[2])
  )
  rm(Boston, unique_vals)
  
  # 15. MDRR Bioconcentration (caret) - 528 obs, 342 covariates
  # Stress Test: Massive covariate space (P=342). Excellent for testing adaptive variable selection.
  # Target: Molecule activity class (mdrrClass).
  data(mdrr, package = "caret")
  mdrr_data <- na.omit(data.frame(Class = mdrrClass, mdrrDescr))
  unique_vals <- sort(unique(mdrr_data$Class))
  benchmark_datasets$MDRR <- list(
    X = model.matrix(Class ~ . - 1, data = mdrr_data),
    Y = as.numeric(mdrr_data$Class == unique_vals[2])
  )
  rm(mdrrClass, mdrrDescr, mdrr_data, unique_vals)
  
  # 16. Breast Cancer Biopsy (MASS) - 683 obs (after NA drop), 9 covariates
  # Note: 'ID' column is removed as it holds no predictive value.
  # Target: Tumor class (benign vs malignant).
  data(biopsy, package = "MASS")
  biopsy <- na.omit(biopsy)
  biopsy <- biopsy[, names(biopsy) != "ID"]
  unique_vals <- sort(unique(biopsy$class))
  benchmark_datasets$biopsy <- list(
    X = model.matrix(class ~ . - 1, data = biopsy),
    Y = as.numeric(biopsy$class == unique_vals[2])
  )
  rm(biopsy, unique_vals)
  
  # 17. Breast Cancer (mlbench) - 699 obs, 9 covariates
  data(BreastCancer, package = "mlbench")
  BreastCancer<-na.omit(BreastCancer)
  benchmark_datasets$BreastCancer <- list(
    X = model.matrix(Class ~ . - 1, data = BreastCancer),
    Y = as.numeric(BreastCancer$Class)
  )
  rm(BreastCancer)
  
  # 18. Pima Indians Diabetes (mlbench) - 768 obs, 8 covariates
  data(PimaIndiansDiabetes, package = "mlbench")
  benchmark_datasets$PimaDiabetes <- list(
    X = model.matrix(diabetes ~ . - 1, data = PimaIndiansDiabetes),
    Y = as.numeric(PimaIndiansDiabetes$diabetes)
  )
  rm(PimaIndiansDiabetes)
  
  # 19. German Credit Data (caret) - 1000 obs, 61 covariates (Binary Classification)
  # Standard benchmark for financial risk classification.
  data(GermanCredit, package = "caret")
  GermanCredit <- na.omit(GermanCredit)
  benchmark_datasets$GermanCredit <- list(
    X = model.matrix(Class ~ . - 1, data = GermanCredit),
    Y = as.numeric(GermanCredit$Class == "Good")
  )
  rm(GermanCredit)
  
  # 20. Orange Juice Brand Purchase (ISLR2) - 1070 obs, 17 covariates (Binary Classification)
  # Target: Purchase of Minute Maid (MM) vs Citrus Hill (CH).
  data(OJ, package = "ISLR2")
  OJ <- na.omit(OJ)
  benchmark_datasets$OJ <- list(
    X = model.matrix(Purchase ~ . - 1, data = OJ),
    Y = as.numeric(OJ$Purchase == "MM")
  )
  rm(OJ)
  
  # 21. Weekly Market Direction (ISLR2) - 1089 obs, 7 covariates (Binary Classification)
  # Similar to Smarket but weekly data over 21 years.
  data(Weekly, package = "ISLR2")
  Weekly <- na.omit(Weekly)
  Weekly <- Weekly[, names(Weekly) != "Today"] # Drop the continuous version of the target
  benchmark_datasets$Weekly <- list(
    X = model.matrix(Direction ~ . - 1, data = Weekly),
    Y = as.numeric(Weekly$Direction == "Up")
  )
  rm(Weekly)
  
  # 22. Stock Market Direction (ISLR2) - 1250 obs, 7 covariates (Binary Classification)
  # Target: Did the market go Up (1) or Down (0)? (Note: 'Today' col removed as it is the target continuous return)
  data(Smarket, package = "ISLR2")
  Smarket <- na.omit(Smarket)
  Smarket <- Smarket[, names(Smarket) != "Today"] # Drop the continuous version of the target
  benchmark_datasets$Smarket <- list(
    X = model.matrix(Direction ~ . - 1, data = Smarket),
    Y = as.numeric(Smarket$Direction == "Up")
  )
  rm(Smarket)
  
  # 23. Cell Body Segmentation (caret) - 2019 obs, 58 covariates
  # Target: Cell class (PS vs WS).
  data(segmentationData, package = "caret")
  segmentationData <- na.omit(segmentationData)
  unique_vals <- sort(unique(segmentationData$Class))
  benchmark_datasets$Segmentation <- list(
    X = model.matrix(Class ~ . - 1, data = segmentationData),
    Y = as.numeric(segmentationData$Class == unique_vals[2])
  )
  rm(segmentationData, unique_vals)
  
  # 24. Titanic Survival (datasets) - 2201 obs, 3 covariates
  # Stress Test: Large N, purely categorical predictors.
  # Note: Converts base R frequency table into 2201 individual binary rows.
  data(Titanic, package = "datasets")
  Titanic_df <- as.data.frame(Titanic)
  # Expand the frequency table into individual rows
  Titanic_expanded <- Titanic_df[rep(1:nrow(Titanic_df), Titanic_df$Freq), names(Titanic_df) != "Freq"]
  Titanic_expanded <- na.omit(Titanic_expanded)
  unique_vals <- sort(unique(Titanic_expanded$Survived))
  
  benchmark_datasets$Titanic <- list(
    X = model.matrix(Survived ~ . - 1, data = Titanic_expanded),
    Y = as.numeric(Titanic_expanded$Survived == unique_vals[2])
  )
  rm(Titanic, Titanic_df, Titanic_expanded, unique_vals)
  
  # 25. AIDS Clinical Trial Group 320 (MASS) - 2843 obs, 4 covariates
  # Stress Test: Large N, very small P.
  # Target: Status indicator (A vs D).
  data(Aids2, package = "MASS")
  Aids2 <- na.omit(Aids2)
  unique_vals <- sort(unique(Aids2$status))
  benchmark_datasets$Aids2 <- list(
    X = model.matrix(status ~ . - 1, data = Aids2),
    Y = as.numeric(Aids2$status == unique_vals[2])
  )
  rm(Aids2, unique_vals)
  
  # 26. Caravan Insurance (ISLR2) - 5822 obs, 85 covariates (Binary Classification)
  # Stress Test: Large N combined with high dimensionality (P=85).
  data(Caravan, package = "ISLR2")
  Caravan <- na.omit(Caravan)
  benchmark_datasets$Caravan <- list(
    X = model.matrix(Purchase ~ . - 1, data = Caravan),
    Y = as.numeric(Caravan$Purchase == "Yes")
  )
  rm(Caravan)
  
  # 27. Credit Card Default (ISLR2) - 10,000 obs, 3 covariates (Binary Classification)
  # Stress Test: Large N scalability with heavy class imbalance.
  data(Default, package = "ISLR2")
  Default <- na.omit(Default)
  benchmark_datasets$Default <- list(
    X = model.matrix(default ~ . - 1, data = Default),
    Y = as.numeric(Default$default == "Yes")
  )
  rm(Default)
  
  # Final memory reclamation
  gc()
}
