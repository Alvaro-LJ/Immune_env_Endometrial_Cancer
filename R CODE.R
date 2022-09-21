#Dependencies
library(tidyverse)
library(readxl)
library(corrplot)
library(gplots)
library(mlr)
library(pROC)
library(rpart.plot)
library(parallel)
library(parallelMap)

#Code to perform Clustering analysis
CELL_DENSITIES <- read_excel("DATA.xlsx", sheet = "Cell_densities")

Immune_vars <- CELL_DENSITIES %>% select(-c(1:7))

Immune_vars_scaled <- scale(Immune_vars)

Dist_immune <- dist(Immune_vars_scaled, method = "manhattan")

Clust_immune <- hclust(Dist_immune, method="ward.D2")


#Lets define the optimal number of clusters 

#First create a function that will return a list of internal measure indexes
#given data, number of clusters and a distance matrix
cluster_metrics <- function(data, clusters, dist_matrix){
  list(db = clusterSim::index.DB(data, clusters)$DB,
       G1 = clusterSim::index.G1(data, clusters),
       clusters = length(unique(clusters))
  )
}

#we will do a function that create a list of bootstrapped tibles of our original tibble
ImmuneBoot <- map(1:100, ~ {Immune_vars_scaled %>% as_tibble() %>% sample_n(size=nrow(.), replace=T)})
str(ImmuneBoot)
#lets combine both functions

metricsTib <- map_df(ImmuneBoot, function(boot){
  d <- dist(boot, method="manhattan")
  cl <- hclust(d, method="ward.D2")
  map_df(2:20, function(k){
    cut <- cutree(cl, k=k)
    cluster_metrics(data=boot, clusters=cut, dist_matrix=d)
  })
})

metricsTib <- metricsTib %>%
  gather(key="Metric", value="Value", -clusters)
metricsTib


ggplot(metricsTib, aes(as.factor(clusters), Value)) +
  facet_wrap(~ Metric, scales = "free_y") +
  geom_line(stat = "summary", fun.y = "mean", aes(group = 1)) +
  stat_summary(fun.data="mean_cl_boot",  
               geom="crossbar", width = 0.5, fill = "white") +
  theme_bw()

#Perform dendrogram spliting
CELL_DENSITIES <- CELL_DENSITIES %>% mutate(Cluster = cutree(Clust_immune, k=5))


#Perform clinical outcome prediction analysis
CLINICAL <- read_excel("DATA.xlsx", sheet = "Clinical")

CLINICAL_filtered <- CLINICAL %>%  select(contains(c("Has")), RELAPSE_UPDATE, GRADE, LVI, FIGO)
CLINICAL_filtered <- CLINICAL_filtered %>% mutate_if(is.double, as.factor)

Prediction_task <- makeClassifTask(data = CLINICAL_filtered, target = "RELAPSE_UPDATE")

logistic_rg <- makeLearner("classif.logreg", predict.type = "prob")

featureSelection <- makeFeatSelControlExhaustive()

set.seed(21)
CrossVal <- makeResampleDesc(method = "RepCV", folds = 2, reps = 20, stratify = T)

parallelStartSocket(cpus = detectCores())

selFeatures <- selectFeatures(learner = logistic_rg, task = Prediction_task,
                              resampling = CrossVal, control = featureSelection)


model_alternative <- glm(data = CLINICAL_filtered, RELAPSE_UPDATE ~ Has_1 + Has_2 + Has_4 + Has_5, family = binomial)

model1 <- glm(data = CLINICAL_filtered, RELAPSE_UPDATE ~ FIGO + LVI + GRADE, family = binomial)

classic_model <- vector(mode = "double", length = 10000)
alternative_model <- vector(mode = "double", length = 10000)

for(i in (1:length(classic_model))){
  Patient_filtered1 <- CLINICAL_filtered %>% sample_n(size = nrow(CLINICAL_filtered), replace = T)
  classic_model[i] <- roc(Patient_filtered1$RELAPSE_UPDATE, predict(model1, newdata = Patient_filtered1, type="response"))$auc
  alternative_model[i] <- roc(Patient_filtered1$RELAPSE_UPDATE, predict(model_alternative, newdata = Patient_filtered1, type="response"))$auc
}

boot_results <- as_tibble(data.frame(CLASSIC = classic_model, Alternative = alternative_model, Boot = 1:10000))
boot_results$Improvement <- boot_results$CLASSIC < boot_results$Alternative

sum(boot_results$Improvement)
quantile(boot_results$Alternative, 0.05)
quantile(boot_results$Alternative, 0.95)

quantile(boot_results$CLASSIC, 0.05)
quantile(boot_results$CLASSIC, 0.95)


CLINICAL_filtered$Predictions_alternative <- predict(model_alternative, newdata = CLINICAL_filtered, type="response")

CLINICAL_filtered$Predictions_classic <- predict(model1, newdata = CLINICAL_filtered, type="response")

Standard_coords <- coords(roc(data = CLINICAL_filtered, RELAPSE_UPDATE, Predictions_classic), x="best", input="threshold", best.method = "youden")
Alternative_coords <- coords(roc(data = CLINICAL_filtered, RELAPSE_UPDATE, Predictions_alternative), x="best", input="threshold",
                             best.method = "youden")

CLINICAL_FILTERED <- CLINICAL_filtered %>% mutate(Standard_model = factor(case_when(Predictions_classic >= Standard_coords$threshold ~ "HIGH RISK",
                                                                          TRUE ~ "LOW RISK"), levels = c("LOW RISK", "HIGH RISK")),
                             Alternative_model = factor(case_when(Predictions_alternative >= Alternative_coords$threshold ~ "HIGH RISK",
                                                                             TRUE ~ "LOW RISK"), levels = c("LOW RISK", "HIGH RISK")))



#Code to fit decision tree to stratify clusters according to immune densities
CELL_DENSITIES <- read_excel("DATA.xlsx", sheet = "Cell_densities")

Algorithm_var <- CELL_DENSITIES %>% select(Cluster, CD68_Epithelium, CD8_Epithelium, FoxP3_Epithelium, CK_PDL1_Epithelium)

Algorithm_var[1] <- as.factor(Algorithm_var[[1]])

Cluster_task <- makeClassifTask(data = Algorithm_var, target = "Cluster")

tree_learner <- makeLearner("classif.rpart", par.vals = list("maxdepth" = 4, "cp" = 0.001))

TreeModel<-train(learner = tree_learner, task = Cluster_task)
TreeModel

treeModelData<-getLearnerModel(TreeModel)

rpart.plot(treeModelData, roundint=F,
           box.palette = list("#f24646", "#ebbb1e", "#018f25", "#0c5deb", "#f56ed5"),
           type=5)