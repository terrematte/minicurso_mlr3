A Preprocessing analysis of clinical data of TCGA-KIRC patients and
building a model with mrl3
================

This project contains a pipeline for analysis of The Cancer Genome Atlas
Kidney - Renal Clear Cell Carcinoma (TCGA-KIRC) clinical data, from
[Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/).

``` r
# Avoid duplicate label error of knitr::purl

options(knitr.duplicate.label = 'allow')
# Code to browse the markdown file with renderized images.
knitr::opts_chunk$set(
  fig.path = "figs/tutorial_"
)
```

# Intro

First of all we are going to load required packages and the data. The
data is part of the mlr3data package.

``` r
if(!require("mlr3")){install.packages("mlr3")}  # mlr3 base package

if(!require("mlr3learners")){install.packages("mlr3learners")}  # additional ML algorithms

if(!require("mlr3extralearners")){install.packages("mlr3extralearners")}  # extra ML algorithms

if(!require("mlr3pipelines")){install.packages("mlr3pipelines")} # create ML pipelines

if(!require("mlr3data")){install.packages("mlr3data")}  # another way to obtain data sets

if(!require("mlr3misc")){install.packages("mlr3misc")} # contains some helper functions

if(!require("mlr3tuning")){install.packages("mlr3tuning")} # tuning ML algorithms

if(!require("paradox")){install.packages("paradox")} # hyperparameter space

if(!require("mlr3viz")){install.packages("mlr3viz")}  # autoplot for benchmarks

if(!require("skimr")){install.packages("skimr")} # Compact and Flexible Summaries of Data

if(!require("finalfit")){install.packages("finalfit")} #  Quickly Create Elegant Regression Results Tables and Plots when Modelling

if(!require("tidyverse")){install.packages("tidyverse")} # R packages for data science

if(!require("bestNormalize")){install.packages("bestNormalize")} # Normalizing Transformation Functions 

if(!require("smotefamily")){install.packages("smotefamily")}  # SMOTE algorithm for imbalance correction

if(!require("VennDiagram")){install.packages("VennDiagram")}  # Generate High-Resolution Venn and Euler Plots
```

# Loading data

``` r
load("data/tcga_kirc.RData")
```

# Exploratory Data Analysis

We can use the skimr package in order to get a first overview of the
data:

``` r
## clinical data size : 
dim(kirc_cli)
```

    ## [1] 298  10

``` r
skimr::skim(kirc_cli)
```

|                                                  |           |
| :----------------------------------------------- | :-------- |
| Name                                             | kirc\_cli |
| Number of rows                                   | 298       |
| Number of columns                                | 10        |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |           |
| Column type frequency:                           |           |
| character                                        | 1         |
| factor                                           | 6         |
| numeric                                          | 3         |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |           |
| Group variables                                  | None      |

Data summary

**Variable type: character**

| skim\_variable | n\_missing | complete\_rate | min | max | empty | n\_unique | whitespace |
| :------------- | ---------: | -------------: | --: | --: | ----: | --------: | ---------: |
| patient\_id    |          0 |              1 |  12 |  12 |     0 |       298 |          0 |

**Variable type: factor**

| skim\_variable | n\_missing | complete\_rate | ordered | n\_unique | top\_counts                       |
| :------------- | ---------: | -------------: | :------ | --------: | :-------------------------------- |
| prior.dx       |          0 |           1.00 | FALSE   |         2 | no: 254, yes: 44                  |
| gender         |          0 |           1.00 | FALSE   |         2 | mal: 190, fem: 108                |
| race           |          0 |           1.00 | FALSE   |         4 | whi: 258, bla: 34, asi: 5, not: 1 |
| metastasis     |          2 |           0.99 | FALSE   |         3 | M0: 199, M1: 78, MX: 19           |
| neoplasm       |          0 |           1.00 | FALSE   |         3 | NX: 175, N0: 117, N1: 6           |
| ajcc.stage     |          0 |           1.00 | FALSE   |         4 | T1: 224, T3: 56, T2: 10, T4: 8    |

**Variable type: numeric**

| skim\_variable | n\_missing | complete\_rate |    mean |     sd | p0 |    p25 |    p50 |     p75 | p100 | hist  |
| :------------- | ---------: | -------------: | ------: | -----: | -: | -----: | -----: | ------: | ---: | :---- |
| age            |          0 |              1 |   58.54 |  11.77 | 26 |  50.25 |   59.0 |   66.75 |   86 | ▁▅▇▆▂ |
| status         |          0 |              1 |    0.21 |   0.41 |  0 |   0.00 |    0.0 |    0.00 |    1 | ▇▁▁▁▂ |
| obs.time       |          0 |              1 | 1338.14 | 998.08 |  3 | 476.25 | 1166.5 | 1954.25 | 4537 | ▇▆▃▂▁ |

Filtering rows only with M0 and M1.

``` r
kirc_cli <- kirc_cli %>% 
  dplyr::filter(metastasis %in% c("M0", "M1")) %>%
  droplevels()

skimr::skim(kirc_cli)
```

|                                                  |           |
| :----------------------------------------------- | :-------- |
| Name                                             | kirc\_cli |
| Number of rows                                   | 277       |
| Number of columns                                | 10        |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |           |
| Column type frequency:                           |           |
| character                                        | 1         |
| factor                                           | 6         |
| numeric                                          | 3         |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |           |
| Group variables                                  | None      |

Data summary

**Variable type: character**

| skim\_variable | n\_missing | complete\_rate | min | max | empty | n\_unique | whitespace |
| :------------- | ---------: | -------------: | --: | --: | ----: | --------: | ---------: |
| patient\_id    |          0 |              1 |  12 |  12 |     0 |       277 |          0 |

**Variable type: factor**

| skim\_variable | n\_missing | complete\_rate | ordered | n\_unique | top\_counts                       |
| :------------- | ---------: | -------------: | :------ | --------: | :-------------------------------- |
| prior.dx       |          0 |              1 | FALSE   |         2 | no: 235, yes: 42                  |
| gender         |          0 |              1 | FALSE   |         2 | mal: 182, fem: 95                 |
| race           |          0 |              1 | FALSE   |         4 | whi: 258, bla: 13, asi: 5, not: 1 |
| metastasis     |          0 |              1 | FALSE   |         2 | M0: 199, M1: 78                   |
| neoplasm       |          0 |              1 | FALSE   |         3 | NX: 155, N0: 116, N1: 6           |
| ajcc.stage     |          0 |              1 | FALSE   |         4 | T1: 203, T3: 56, T2: 10, T4: 8    |

**Variable type: numeric**

| skim\_variable | n\_missing | complete\_rate |    mean |     sd | p0 | p25 |  p50 |  p75 | p100 | hist  |
| :------------- | ---------: | -------------: | ------: | -----: | -: | --: | ---: | ---: | ---: | :---- |
| age            |          0 |              1 |   58.43 |  11.84 | 26 |  51 |   59 |   66 |   86 | ▁▅▇▆▂ |
| status         |          0 |              1 |    0.23 |   0.42 |  0 |   0 |    0 |    0 |    1 | ▇▁▁▁▂ |
| obs.time       |          0 |              1 | 1394.49 | 998.93 | 11 | 552 | 1257 | 1986 | 4537 | ▇▇▅▂▁ |

## Cleaning expression data and pre-selecting genes

``` r
dim(kirc_rna)
```

    ## [1]   301 58387

``` r
head(kirc_rna[, c(1:10)])
```

    ##                  TSPAN6 TNMD DPM1 SCYL3 C1orf112  FGR   CFH FUCA2 GCLC NFYA
    ## TCGA-A3-3387-01A   2261   17 1642   895      361 1903  6194  4448 2623 1952
    ## TCGA-BP-4769-01A   3101   56 1669   891      173  377  5319  3266 1531 1753
    ## TCGA-BP-4977-01A   5404  312 1471  1278      356 1041  3807  7086 3236 2131
    ## TCGA-B0-5080-01A   3808    8 1615   724      247  318 11994  3901 1095 1853
    ## TCGA-CZ-4862-01A   1955   12 1275   948      268 1255  3091  3829 2310 1529
    ## TCGA-BP-4758-01A   2099   43 1096   820      167 1244  2371  4092 1496 1630

``` r
# Check if there are duplicated gene symbols
colnames(kirc_rna)[duplicated(colnames(kirc_rna))]
```

    ## character(0)

We performed a differential expression analysis to select differentially
expressed genes, on script `job_differential_gene_expression.R`

We also selected a list with 252 genes of papers on genes signatures,
obtained from search of Pubmed with the keywords: `renal AND ‘gene
signature’ OR kidney AND ‘gene signature’`

We also selected all genes mapped from Kegg:

<https://www.genome.jp/dbget-bin/www_bget?path:map05211>

``` r
genes_DEA_M1 <- readLines("data/dea.M0.M1.lst")
genes_papers <- readLines("data/genes_papers.lst")
genes_kegg <- readLines("data/genes_kegg.lst")

genes <- union(genes_DEA_M1, union(names(genes_papers), names(genes_kegg)))
patients_id <-  rownames(kirc_rna) %in% rownames(kirc_cli)
  
kirc_rna <- kirc_rna[patients_id, genes] 

dim(kirc_rna)
```

    ## [1] 277 150

``` r
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
```

    ## NULL

``` r
venn.diagram(
  x = list(dea = genes_DEA_M1, papers = genes_papers, kegg = genes_kegg),
  cat.just=list(c(0.5,1) , c(2,-1) , c(-1,-22)),
  height = 1200, width = 1200,
  resolution = 200,
  filename = "figs/selected_features.png", 
  imagetype = "png",
  col = "black",
  fill = c("khaki1", "skyblue", "tomato3"),
  alpha = 0.50,
  lwd = 4,
  #cat.cex = 1.2,
  #cex = 1.5,
  cat.cex = 1,
  cex = 1,
  cat.fontface = "bold"
)
```

    ## [1] 1

``` r
knitr::include_graphics("figs/selected_features.png", dpi = NA)
```

<img src="figs/selected_features.png" width="200px" />

Assertion on ‘feature names’: Must have names according to R’s variable
naming conventions.

``` r
# Rename columns, removing "-" 
colnames(kirc_rna) <- gsub("-", "_", colnames(kirc_rna))
# Minimum count is set to 1 in order to prevent 0 division problem within classification models.
kirc_rna <- (kirc_rna +1)
```

Selecting the data to classify metastasis

``` r
kirc_data <- as.data.frame(kirc_rna)
kirc_data$metastasis <- kirc_cli$metastasis
```

# A first model

Setting up the task and learner `rpart`: Recursive Partitioning and
Regression Trees

List of learners:
<https://mlr3extralearners.mlr-org.com/articles/learners/list_learners.html>

``` r
head(kirc_data[,c(1:4)])
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["SLC4A1"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["HHATL"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["SLC38A5"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["ZIC2"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"9","2":"2","3":"87","4":"6","_rn_":"TCGA-A3-3387-01A"},{"1":"9","2":"5","3":"759","4":"1","_rn_":"TCGA-BP-4769-01A"},{"1":"1025","2":"32","3":"101","4":"9","_rn_":"TCGA-BP-4977-01A"},{"1":"7","2":"6","3":"2799","4":"4","_rn_":"TCGA-B0-5080-01A"},{"1":"9","2":"12","3":"358","4":"2","_rn_":"TCGA-CZ-4862-01A"},{"1":"30","2":"81","3":"152","4":"1","_rn_":"TCGA-BP-4758-01A"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
tsk_raw <- TaskClassif$new(id="kirc_raw", 
                           backend = kirc_data, 
                           target = "metastasis", 
                           positive = "M1")

p_bc = po("boxcox", 
          affect_columns = selector_type("numeric"))

kirc_norm = p_bc$train(list(tsk_raw))$output$data()

head(kirc_norm[,c(1:4)])
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["metastasis"],"name":[1],"type":["fct"],"align":["left"]},{"label":["AC003092.1"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["AC006262.4"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["AC006262.5"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"M0","2":"1.0007861","3":"0.9466137","4":"0.66346732"},{"1":"M0","2":"-0.4147400","3":"-1.1683487","4":"0.77201847"},{"1":"M0","2":"0.4485704","3":"-0.2441639","4":"-1.26501108"},{"1":"M1","2":"0.9271871","3":"1.8700639","4":"1.94150301"},{"1":"M0","2":"1.2527336","3":"0.8604927","4":"0.08017418"},{"1":"M0","2":"-0.4147400","3":"-1.1683487","4":"-1.26501108"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
tsk_cla <- TaskClassif$new(id="kirc_cla", 
                           backend = kirc_norm, 
                           target = "metastasis", 
                           positive = "M1")


tsk_cla
```

    ## <TaskClassif:kirc_cla> (277 x 151)
    ## * Target: metastasis
    ## * Properties: twoclass
    ## * Features (150):
    ##   - dbl (150): AC003092.1, AC006262.4, AC006262.5, AC007879.6,
    ##     AC104654.2, AC116614.1, AFM, AHSG, AMH, APCDD1L_AS1, AQP2, AQP6,
    ##     ATP6V0A4, ATP6V0D2, BARX1, BSND, C10orf99, C14orf180, CA1, CASP14,
    ##     CCNA1, CDC42P2, CHAT, CIDEC, CILP2, CITF22_24E5.1, CLCNKB, CLDN8,
    ##     CLMP, COL11A1, COL7A1, CPNE7, CTD_2008P7.9, CXCL13, CYP1A1, DMRT2,
    ##     DQX1, EN2, ESRP1, FAM83B, FDCSP, FGF5, FKBP9P1, FOXI2, FXYD4, GGT6,
    ##     GLB1L3, GOLGA6L2, GOLGA6L7P, GPR110, HEPACAM2, HHATL, HMGA2,
    ##     HS3ST3A1, IGF2BP3, IGFBP1, IGFL1P1, IGFN1, IGHV1_69, IGKV3_11,
    ##     IGLC7, IGLV3_19, INHBE, ITPKA, KCNJ1, KIRREL3, KLF17, KLK1, KLK15,
    ##     KNG1, KRT7, L1CAM, LECT1, LINC00890, LINC00942, LINC00973,
    ##     LINC01187, LINC01436, LINC01559, LRRTM1, MAGEC2, MAGEC3, MFI2,
    ##     MYH8, NFE4, NIPAL4, NKX2_2, NKX2_3, NMRK2, NUPR1L, OTX1, PADI3,
    ##     PAEP, PI3, PITX1, PLG, PRR15L, PSG9, PVALB, RAB25, [...]

## Train and Predict

Setting up the train/test splits of the data

``` r
set.seed(1)
train_set = sample(tsk_cla$nrow,  0.7 * tsk_cla$nrow)

test_set = setdiff(seq_len(tsk_cla$nrow), train_set)
```

The field `$model` stores the model that is produced in the training
step. Before the `$train()` method is called on a learner object, this
field is `NULL`:

``` r
learner = lrn("classif.rpart")
learner$model
```

    ## NULL

Next, the classification tree is trained using the train set of the task
by calling the $train() method of the Learner:

``` r
set.seed(1)

learner$train(tsk_cla, row_ids = train_set)

print(learner$model)
```

    ## n= 193 
    ## 
    ## node), split, n, loss, yval, (yprob)
    ##       * denotes terminal node
    ## 
    ##  1) root 193 52 M0 (0.26943005 0.73056995)  
    ##    2) TNNT1>=0.7619996 46 16 M1 (0.65217391 0.34782609)  
    ##      4) CLDN8< 0.4817347 30  4 M1 (0.86666667 0.13333333) *
    ##      5) CLDN8>=0.4817347 16  4 M0 (0.25000000 0.75000000) *
    ##    3) TNNT1< 0.7619996 147 22 M0 (0.14965986 0.85034014)  
    ##      6) CHAT>=1.561052 11  4 M1 (0.63636364 0.36363636) *
    ##      7) CHAT< 1.561052 136 15 M0 (0.11029412 0.88970588)  
    ##       14) C10orf99>=1.343281 9  4 M1 (0.55555556 0.44444444) *
    ##       15) C10orf99< 1.343281 127 10 M0 (0.07874016 0.92125984) *

## Predicting

``` r
prediction = learner$predict(tsk_cla, row_ids = test_set)
prediction$confusion
```

    ##         truth
    ## response M1 M0
    ##       M1  9 11
    ##       M0 17 47

``` r
prediction$score( msr("classif.acc"))
```

    ## classif.acc 
    ##   0.6666667

## Evaluating with distincs measures

`View(as.data.table(mlr_measures))`

``` r
measures = list(
  msr("classif.acc"), 
  msr("classif.bacc"),
  msr("classif.precision"),
  msr("classif.sensitivity"), 
  msr("classif.specificity")
  )

prediction$score(measures)
```

    ##         classif.acc        classif.bacc   classif.precision classif.sensitivity 
    ##           0.6666667           0.5782493           0.4500000           0.3461538 
    ## classif.specificity 
    ##           0.8103448

## Resampling

Setting up our resampling method

``` r
rsmp_cv = rsmp("cv", folds = 3L)$instantiate(tsk_cla)

res = resample(task = tsk_cla, 
               learner = learner, 
               resampling = rsmp_cv,
               store_models = TRUE)
```

    ## INFO  [15:59:34.582] [mlr3]  Applying learner 'classif.rpart' on task 'kirc_cla' (iter 3/3) 
    ## INFO  [15:59:34.944] [mlr3]  Applying learner 'classif.rpart' on task 'kirc_cla' (iter 1/3) 
    ## INFO  [15:59:35.004] [mlr3]  Applying learner 'classif.rpart' on task 'kirc_cla' (iter 2/3)

``` r
measures <- list(
  msr("classif.acc"),
  msr("classif.bacc"),
  msr("classif.precision"), 
  msr("classif.sensitivity"), 
  msr("classif.specificity")
)

agg <- res$aggregate(measures)

agg
```

    ##         classif.acc        classif.bacc   classif.precision classif.sensitivity 
    ##           0.7182484           0.6260972           0.5617284           0.4125067 
    ## classif.specificity 
    ##           0.8396876

# Filter Selection - Variable Importance Filters

``` r
tsk_filt <- TaskClassif$new(id="filt_rpart", 
                               backend = kirc_norm, 
                               target = "metastasis", 
                               positive = "M1")


lrn = lrn("classif.rpart")

library("mlr3filters")
filter = flt("importance", learner = lrn)

filter$calculate(tsk_filt)

head(as.data.table(filter), 20)
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["feature"],"name":[1],"type":["chr"],"align":["left"]},{"label":["score"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"ITPKA","2":"27.421568"},{"1":"CLMP","2":"14.195831"},{"1":"AC116614.1","2":"11.505887"},{"1":"TNNT1","2":"10.515757"},{"1":"ZIC2","2":"9.200185"},{"1":"INHBE","2":"8.974032"},{"1":"KLF17","2":"7.360148"},{"1":"COL7A1","2":"7.360148"},{"1":"SAA1","2":"7.279029"},{"1":"RP11_643A5.3","2":"6.944444"},{"1":"GLB1L3","2":"6.867568"},{"1":"SAA2_SAA4","2":"5.665860"},{"1":"CXCL13","2":"5.149545"},{"1":"PITX1","2":"4.495096"},{"1":"PI3","2":"3.901438"},{"1":"CPNE7","2":"3.576154"},{"1":"SLC22A8","2":"3.000820"},{"1":"CA1","2":"2.993203"},{"1":"RP11_586K2.1","2":"2.903704"},{"1":"RP11_440G9.1","2":"2.798729"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
cols <- head(as.data.table(filter), 20)$feature

tsk_filt$select(cols = cols)
```

# Feature selection - RFE

``` r
library(mlr3fselect)

tsk_rfe <- TaskClassif$new(id="rfe_part", 
                               backend = kirc_norm, 
                               target = "metastasis", 
                               positive = "M1")


terminator = trm("evals", n_evals = 100)

instance = FSelectInstanceSingleCrit$new(
  task = tsk_cla,
  learner = lrn("classif.rpart"),
  resampling = rsmp("cv", folds = 5),
  measure = msr("classif.bacc"),
  terminator = terminator,
  store_models = T
)

# Modifies the instance by reference ----
fselector = fs("rfe", min_features=10)

fselector$optimize(instance)
```

    ## INFO  [15:59:35.712] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRFE>' and '<TerminatorEvals> [n_evals=100]' 
    ## INFO  [15:59:35.715] [bbotk] Evaluating 1 configuration(s) 
    ## INFO  [15:59:35.813] [mlr3]  Running benchmark with 5 resampling iterations 
    ## INFO  [15:59:35.820] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/5) 
    ## INFO  [15:59:35.949] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 3/5) 
    ## INFO  [15:59:36.081] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 2/5) 
    ## INFO  [15:59:36.208] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 5/5) 
    ## INFO  [15:59:36.347] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 4/5) 
    ## INFO  [15:59:36.491] [mlr3]  Finished benchmark 
    ## INFO  [15:59:36.551] [bbotk] Result of batch 1: 
    ## INFO  [15:59:36.560] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [15:59:36.560] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE TRUE TRUE 
    ## INFO  [15:59:36.560] [bbotk]   AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [15:59:36.560] [bbotk]  TRUE        TRUE TRUE TRUE     TRUE     TRUE  TRUE TRUE     TRUE      TRUE 
    ## INFO  [15:59:36.560] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:36.560] [bbotk]  TRUE   TRUE  TRUE    TRUE TRUE  TRUE  TRUE          TRUE   TRUE  TRUE TRUE 
    ## INFO  [15:59:36.560] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [15:59:36.560] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE TRUE TRUE  TRUE   TRUE 
    ## INFO  [15:59:36.560] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [15:59:36.560] [bbotk]   TRUE TRUE    TRUE  TRUE  TRUE TRUE   TRUE     TRUE      TRUE   TRUE     TRUE 
    ## INFO  [15:59:36.560] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [15:59:36.560] [bbotk]   TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE  TRUE 
    ## INFO  [15:59:36.560] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [15:59:36.560] [bbotk]      TRUE  TRUE  TRUE  TRUE    TRUE  TRUE TRUE  TRUE TRUE TRUE  TRUE  TRUE 
    ## INFO  [15:59:36.560] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [15:59:36.560] [bbotk]       TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE   TRUE 
    ## INFO  [15:59:36.560] [bbotk]  MAGEC3 MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP  PI3 
    ## INFO  [15:59:36.560] [bbotk]    TRUE TRUE TRUE TRUE   TRUE   TRUE   TRUE  TRUE   TRUE TRUE  TRUE TRUE TRUE 
    ## INFO  [15:59:36.560] [bbotk]  PITX1  PLG PRR15L PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [15:59:36.560] [bbotk]   TRUE TRUE   TRUE TRUE  TRUE  TRUE TRUE         TRUE          TRUE 
    ## INFO  [15:59:36.560] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [15:59:36.560] [bbotk]           TRUE         TRUE          TRUE          TRUE          TRUE 
    ## INFO  [15:59:36.560] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [15:59:36.560] [bbotk]          TRUE        TRUE          TRUE         TRUE         TRUE          TRUE 
    ## INFO  [15:59:36.560] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [15:59:36.560] [bbotk]          TRUE         TRUE         TRUE        TRUE TRUE TRUE      TRUE TRUE 
    ## INFO  [15:59:36.560] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [15:59:36.560] [bbotk]   TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE   TRUE    TRUE 
    ## INFO  [15:59:36.560] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [15:59:36.560] [bbotk]     TRUE    TRUE  TRUE   TRUE TRUE    TRUE   TRUE      TRUE  TRUE TRUE   TRUE 
    ## INFO  [15:59:36.560] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5 ZIC2 ZIC5 classif.bacc 
    ## INFO  [15:59:36.560] [bbotk]    TRUE  TRUE    TRUE TRUE  TRUE TRUE TRUE    0.6936458 
    ## INFO  [15:59:36.560] [bbotk]                                 uhash 
    ## INFO  [15:59:36.560] [bbotk]  559853be-1e45-4976-851d-86318e8b2c21 
    ## INFO  [15:59:36.662] [bbotk] Evaluating 1 configuration(s) 
    ## INFO  [15:59:36.717] [mlr3]  Running benchmark with 5 resampling iterations 
    ## INFO  [15:59:36.724] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 3/5) 
    ## INFO  [15:59:36.829] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 2/5) 
    ## INFO  [15:59:36.933] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/5) 
    ## INFO  [15:59:37.038] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 5/5) 
    ## INFO  [15:59:37.145] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 4/5) 
    ## INFO  [15:59:37.258] [mlr3]  Finished benchmark 
    ## INFO  [15:59:37.340] [bbotk] Result of batch 2: 
    ## INFO  [15:59:37.348] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [15:59:37.348] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE       TRUE TRUE TRUE 
    ## INFO  [15:59:37.348] [bbotk]    AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:37.348] [bbotk]  FALSE        TRUE TRUE TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [15:59:37.348] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:37.348] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE FALSE  TRUE         FALSE   TRUE FALSE TRUE 
    ## INFO  [15:59:37.348] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [15:59:37.348] [bbotk]    FALSE   TRUE FALSE        FALSE   TRUE   TRUE  TRUE TRUE TRUE  TRUE  FALSE 
    ## INFO  [15:59:37.348] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [15:59:37.348] [bbotk]  FALSE TRUE   FALSE FALSE FALSE TRUE  FALSE    FALSE     FALSE   TRUE     TRUE 
    ## INFO  [15:59:37.348] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [15:59:37.348] [bbotk]  FALSE  TRUE     TRUE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE FALSE 
    ## INFO  [15:59:37.348] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [15:59:37.348] [bbotk]     FALSE  TRUE  TRUE FALSE    TRUE  TRUE FALSE  TRUE TRUE TRUE  TRUE FALSE 
    ## INFO  [15:59:37.348] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [15:59:37.348] [bbotk]      FALSE      TRUE     FALSE     FALSE      TRUE      TRUE  FALSE  FALSE 
    ## INFO  [15:59:37.348] [bbotk]  MAGEC3 MFI2 MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP  PI3 
    ## INFO  [15:59:37.348] [bbotk]   FALSE TRUE TRUE FALSE   TRUE  FALSE  FALSE FALSE  FALSE TRUE  TRUE TRUE TRUE 
    ## INFO  [15:59:37.348] [bbotk]  PITX1   PLG PRR15L  PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [15:59:37.348] [bbotk]   TRUE FALSE   TRUE FALSE FALSE FALSE TRUE        FALSE          TRUE 
    ## INFO  [15:59:37.348] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [15:59:37.348] [bbotk]           TRUE        FALSE         FALSE         FALSE         FALSE 
    ## INFO  [15:59:37.348] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [15:59:37.348] [bbotk]          TRUE        TRUE         FALSE        FALSE        FALSE         FALSE 
    ## INFO  [15:59:37.348] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1 SAA1 SAA2_SAA4  SBSN 
    ## INFO  [15:59:37.348] [bbotk]          TRUE        FALSE         TRUE        TRUE FALSE TRUE      TRUE FALSE 
    ## INFO  [15:59:37.348] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [15:59:37.348] [bbotk]  FALSE FALSE   FALSE   FALSE    TRUE   FALSE    TRUE    TRUE   TRUE   FALSE 
    ## INFO  [15:59:37.348] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [15:59:37.348] [bbotk]    FALSE   FALSE FALSE  FALSE TRUE    TRUE   TRUE     FALSE  TRUE TRUE   TRUE 
    ## INFO  [15:59:37.348] [bbotk]  TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2 ZIC5 classif.bacc 
    ## INFO  [15:59:37.348] [bbotk]   FALSE  TRUE    TRUE FALSE FALSE TRUE TRUE    0.6936458 
    ## INFO  [15:59:37.348] [bbotk]                                 uhash 
    ## INFO  [15:59:37.348] [bbotk]  b69c92c1-a849-4a6f-a2e3-216b91d2975b 
    ## INFO  [15:59:37.350] [bbotk] Evaluating 1 configuration(s) 
    ## INFO  [15:59:37.404] [mlr3]  Running benchmark with 5 resampling iterations 
    ## INFO  [15:59:37.411] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/5) 
    ## INFO  [15:59:37.472] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 3/5) 
    ## INFO  [15:59:37.537] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 2/5) 
    ## INFO  [15:59:37.602] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 4/5) 
    ## INFO  [15:59:37.667] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 5/5) 
    ## INFO  [15:59:37.739] [mlr3]  Finished benchmark 
    ## INFO  [15:59:37.792] [bbotk] Result of batch 3: 
    ## INFO  [15:59:37.799] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM  AHSG 
    ## INFO  [15:59:37.799] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE TRUE FALSE 
    ## INFO  [15:59:37.799] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:37.799] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [15:59:37.799] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:37.799] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE TRUE 
    ## INFO  [15:59:37.799] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [15:59:37.799] [bbotk]    FALSE   TRUE FALSE        FALSE   TRUE   TRUE FALSE TRUE TRUE  TRUE  FALSE 
    ## INFO  [15:59:37.799] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [15:59:37.799] [bbotk]  FALSE TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE   TRUE    FALSE 
    ## INFO  [15:59:37.799] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [15:59:37.799] [bbotk]  FALSE  TRUE     TRUE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE FALSE 
    ## INFO  [15:59:37.799] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [15:59:37.799] [bbotk]     FALSE  TRUE  TRUE FALSE    TRUE  TRUE FALSE FALSE TRUE TRUE FALSE FALSE 
    ## INFO  [15:59:37.799] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [15:59:37.799] [bbotk]      FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE  FALSE 
    ## INFO  [15:59:37.799] [bbotk]  MAGEC3 MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 PAEP 
    ## INFO  [15:59:37.799] [bbotk]   FALSE TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE TRUE 
    ## INFO  [15:59:37.799] [bbotk]    PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [15:59:37.799] [bbotk]  FALSE  TRUE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE         FALSE 
    ## INFO  [15:59:37.799] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [15:59:37.799] [bbotk]          FALSE        FALSE         FALSE         FALSE         FALSE 
    ## INFO  [15:59:37.799] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [15:59:37.799] [bbotk]         FALSE       FALSE         FALSE        FALSE        FALSE         FALSE 
    ## INFO  [15:59:37.799] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 SAA2_SAA4  SBSN 
    ## INFO  [15:59:37.799] [bbotk]          TRUE        FALSE        FALSE        TRUE FALSE FALSE      TRUE FALSE 
    ## INFO  [15:59:37.799] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [15:59:37.799] [bbotk]  FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE  FALSE   FALSE 
    ## INFO  [15:59:37.799] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [15:59:37.799] [bbotk]    FALSE   FALSE FALSE  FALSE TRUE   FALSE  FALSE     FALSE  TRUE TRUE  FALSE 
    ## INFO  [15:59:37.799] [bbotk]  TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2 ZIC5 classif.bacc 
    ## INFO  [15:59:37.799] [bbotk]   FALSE FALSE   FALSE FALSE FALSE TRUE TRUE    0.6828125 
    ## INFO  [15:59:37.799] [bbotk]                                 uhash 
    ## INFO  [15:59:37.799] [bbotk]  35f1008a-a9f6-4b28-9e6a-83b16811141c 
    ## INFO  [15:59:37.802] [bbotk] Evaluating 1 configuration(s) 
    ## INFO  [15:59:37.857] [mlr3]  Running benchmark with 5 resampling iterations 
    ## INFO  [15:59:37.863] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 2/5) 
    ## INFO  [15:59:37.924] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 5/5) 
    ## INFO  [15:59:37.986] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 4/5) 
    ## INFO  [15:59:38.056] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/5) 
    ## INFO  [15:59:38.114] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 3/5) 
    ## INFO  [15:59:38.174] [mlr3]  Finished benchmark 
    ## INFO  [15:59:38.229] [bbotk] Result of batch 4: 
    ## INFO  [15:59:38.237] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:38.237] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:38.237] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:38.237] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:38.237] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:38.237] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE TRUE 
    ## INFO  [15:59:38.237] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [15:59:38.237] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE FALSE FALSE TRUE  TRUE  FALSE 
    ## INFO  [15:59:38.237] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:38.237] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:38.237] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:38.237] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:38.237] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:38.237] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE    TRUE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:38.237] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:38.237] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:38.237] [bbotk]  MAGEC2 MAGEC3 MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:38.237] [bbotk]   FALSE  FALSE TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:38.237] [bbotk]  PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:38.237] [bbotk]  TRUE FALSE  TRUE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:38.237] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:38.237] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:38.237] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:38.237] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:38.237] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:38.237] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:38.237] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:38.237] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:38.237] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:38.237] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:38.237] [bbotk]  TNNT1  TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:38.237] [bbotk]   TRUE TRUE  FALSE  FALSE FALSE   FALSE FALSE FALSE TRUE FALSE    0.6877206 
    ## INFO  [15:59:38.237] [bbotk]                                 uhash 
    ## INFO  [15:59:38.237] [bbotk]  af5fd2b9-baaa-431b-8093-c90a96a23d5e 
    ## INFO  [15:59:38.245] [bbotk] Finished optimizing after 4 evaluation(s) 
    ## INFO  [15:59:38.245] [bbotk] Result: 
    ## INFO  [15:59:38.251] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [15:59:38.251] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE TRUE TRUE 
    ## INFO  [15:59:38.251] [bbotk]   AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [15:59:38.251] [bbotk]  TRUE        TRUE TRUE TRUE     TRUE     TRUE  TRUE TRUE     TRUE      TRUE 
    ## INFO  [15:59:38.251] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:38.251] [bbotk]  TRUE   TRUE  TRUE    TRUE TRUE  TRUE  TRUE          TRUE   TRUE  TRUE TRUE 
    ## INFO  [15:59:38.251] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [15:59:38.251] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE TRUE TRUE  TRUE   TRUE 
    ## INFO  [15:59:38.251] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [15:59:38.251] [bbotk]   TRUE TRUE    TRUE  TRUE  TRUE TRUE   TRUE     TRUE      TRUE   TRUE     TRUE 
    ## INFO  [15:59:38.251] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [15:59:38.251] [bbotk]   TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE  TRUE 
    ## INFO  [15:59:38.251] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [15:59:38.251] [bbotk]      TRUE  TRUE  TRUE  TRUE    TRUE  TRUE TRUE  TRUE TRUE TRUE  TRUE  TRUE 
    ## INFO  [15:59:38.251] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [15:59:38.251] [bbotk]       TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE   TRUE 
    ## INFO  [15:59:38.251] [bbotk]  MAGEC3 MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP  PI3 
    ## INFO  [15:59:38.251] [bbotk]    TRUE TRUE TRUE TRUE   TRUE   TRUE   TRUE  TRUE   TRUE TRUE  TRUE TRUE TRUE 
    ## INFO  [15:59:38.251] [bbotk]  PITX1  PLG PRR15L PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [15:59:38.251] [bbotk]   TRUE TRUE   TRUE TRUE  TRUE  TRUE TRUE         TRUE          TRUE 
    ## INFO  [15:59:38.251] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [15:59:38.251] [bbotk]           TRUE         TRUE          TRUE          TRUE          TRUE 
    ## INFO  [15:59:38.251] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [15:59:38.251] [bbotk]          TRUE        TRUE          TRUE         TRUE         TRUE          TRUE 
    ## INFO  [15:59:38.251] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [15:59:38.251] [bbotk]          TRUE         TRUE         TRUE        TRUE TRUE TRUE      TRUE TRUE 
    ## INFO  [15:59:38.251] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [15:59:38.251] [bbotk]   TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE   TRUE    TRUE 
    ## INFO  [15:59:38.251] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [15:59:38.251] [bbotk]     TRUE    TRUE  TRUE   TRUE TRUE    TRUE   TRUE      TRUE  TRUE TRUE   TRUE 
    ## INFO  [15:59:38.251] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [15:59:38.251] [bbotk]    TRUE  TRUE    TRUE TRUE  TRUE TRUE TRUE 
    ## INFO  [15:59:38.251] [bbotk]                                                               features 
    ## INFO  [15:59:38.251] [bbotk]  AC003092.1,AC006262.4,AC006262.5,AC007879.6,AC104654.2,AC116614.1,... 
    ## INFO  [15:59:38.251] [bbotk]  classif.bacc 
    ## INFO  [15:59:38.251] [bbotk]     0.6936458

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["AC003092.1"],"name":[1],"type":["lgl"],"align":["right"]},{"label":["AC006262.4"],"name":[2],"type":["lgl"],"align":["right"]},{"label":["AC006262.5"],"name":[3],"type":["lgl"],"align":["right"]},{"label":["AC007879.6"],"name":[4],"type":["lgl"],"align":["right"]},{"label":["AC104654.2"],"name":[5],"type":["lgl"],"align":["right"]},{"label":["AC116614.1"],"name":[6],"type":["lgl"],"align":["right"]},{"label":["AFM"],"name":[7],"type":["lgl"],"align":["right"]},{"label":["AHSG"],"name":[8],"type":["lgl"],"align":["right"]},{"label":["AMH"],"name":[9],"type":["lgl"],"align":["right"]},{"label":["APCDD1L_AS1"],"name":[10],"type":["lgl"],"align":["right"]},{"label":["AQP2"],"name":[11],"type":["lgl"],"align":["right"]},{"label":["AQP6"],"name":[12],"type":["lgl"],"align":["right"]},{"label":["ATP6V0A4"],"name":[13],"type":["lgl"],"align":["right"]},{"label":["ATP6V0D2"],"name":[14],"type":["lgl"],"align":["right"]},{"label":["BARX1"],"name":[15],"type":["lgl"],"align":["right"]},{"label":["BSND"],"name":[16],"type":["lgl"],"align":["right"]},{"label":["C10orf99"],"name":[17],"type":["lgl"],"align":["right"]},{"label":["C14orf180"],"name":[18],"type":["lgl"],"align":["right"]},{"label":["CA1"],"name":[19],"type":["lgl"],"align":["right"]},{"label":["CASP14"],"name":[20],"type":["lgl"],"align":["right"]},{"label":["CCNA1"],"name":[21],"type":["lgl"],"align":["right"]},{"label":["CDC42P2"],"name":[22],"type":["lgl"],"align":["right"]},{"label":["CHAT"],"name":[23],"type":["lgl"],"align":["right"]},{"label":["CIDEC"],"name":[24],"type":["lgl"],"align":["right"]},{"label":["CILP2"],"name":[25],"type":["lgl"],"align":["right"]},{"label":["CITF22_24E5.1"],"name":[26],"type":["lgl"],"align":["right"]},{"label":["CLCNKB"],"name":[27],"type":["lgl"],"align":["right"]},{"label":["CLDN8"],"name":[28],"type":["lgl"],"align":["right"]},{"label":["CLMP"],"name":[29],"type":["lgl"],"align":["right"]},{"label":["COL11A1"],"name":[30],"type":["lgl"],"align":["right"]},{"label":["COL7A1"],"name":[31],"type":["lgl"],"align":["right"]},{"label":["CPNE7"],"name":[32],"type":["lgl"],"align":["right"]},{"label":["CTD_2008P7.9"],"name":[33],"type":["lgl"],"align":["right"]},{"label":["CXCL13"],"name":[34],"type":["lgl"],"align":["right"]},{"label":["CYP1A1"],"name":[35],"type":["lgl"],"align":["right"]},{"label":["DMRT2"],"name":[36],"type":["lgl"],"align":["right"]},{"label":["DQX1"],"name":[37],"type":["lgl"],"align":["right"]},{"label":["EN2"],"name":[38],"type":["lgl"],"align":["right"]},{"label":["ESRP1"],"name":[39],"type":["lgl"],"align":["right"]},{"label":["FAM83B"],"name":[40],"type":["lgl"],"align":["right"]},{"label":["FDCSP"],"name":[41],"type":["lgl"],"align":["right"]},{"label":["FGF5"],"name":[42],"type":["lgl"],"align":["right"]},{"label":["FKBP9P1"],"name":[43],"type":["lgl"],"align":["right"]},{"label":["FOXI2"],"name":[44],"type":["lgl"],"align":["right"]},{"label":["FXYD4"],"name":[45],"type":["lgl"],"align":["right"]},{"label":["GGT6"],"name":[46],"type":["lgl"],"align":["right"]},{"label":["GLB1L3"],"name":[47],"type":["lgl"],"align":["right"]},{"label":["GOLGA6L2"],"name":[48],"type":["lgl"],"align":["right"]},{"label":["GOLGA6L7P"],"name":[49],"type":["lgl"],"align":["right"]},{"label":["GPR110"],"name":[50],"type":["lgl"],"align":["right"]},{"label":["HEPACAM2"],"name":[51],"type":["lgl"],"align":["right"]},{"label":["HHATL"],"name":[52],"type":["lgl"],"align":["right"]},{"label":["HMGA2"],"name":[53],"type":["lgl"],"align":["right"]},{"label":["HS3ST3A1"],"name":[54],"type":["lgl"],"align":["right"]},{"label":["IGF2BP3"],"name":[55],"type":["lgl"],"align":["right"]},{"label":["IGFBP1"],"name":[56],"type":["lgl"],"align":["right"]},{"label":["IGFL1P1"],"name":[57],"type":["lgl"],"align":["right"]},{"label":["IGFN1"],"name":[58],"type":["lgl"],"align":["right"]},{"label":["IGHV1_69"],"name":[59],"type":["lgl"],"align":["right"]},{"label":["IGKV3_11"],"name":[60],"type":["lgl"],"align":["right"]},{"label":["IGLC7"],"name":[61],"type":["lgl"],"align":["right"]},{"label":["IGLV3_19"],"name":[62],"type":["lgl"],"align":["right"]},{"label":["INHBE"],"name":[63],"type":["lgl"],"align":["right"]},{"label":["ITPKA"],"name":[64],"type":["lgl"],"align":["right"]},{"label":["KCNJ1"],"name":[65],"type":["lgl"],"align":["right"]},{"label":["KIRREL3"],"name":[66],"type":["lgl"],"align":["right"]},{"label":["KLF17"],"name":[67],"type":["lgl"],"align":["right"]},{"label":["KLK1"],"name":[68],"type":["lgl"],"align":["right"]},{"label":["KLK15"],"name":[69],"type":["lgl"],"align":["right"]},{"label":["KNG1"],"name":[70],"type":["lgl"],"align":["right"]},{"label":["KRT7"],"name":[71],"type":["lgl"],"align":["right"]},{"label":["L1CAM"],"name":[72],"type":["lgl"],"align":["right"]},{"label":["LECT1"],"name":[73],"type":["lgl"],"align":["right"]},{"label":["LINC00890"],"name":[74],"type":["lgl"],"align":["right"]},{"label":["LINC00942"],"name":[75],"type":["lgl"],"align":["right"]},{"label":["LINC00973"],"name":[76],"type":["lgl"],"align":["right"]},{"label":["LINC01187"],"name":[77],"type":["lgl"],"align":["right"]},{"label":["LINC01436"],"name":[78],"type":["lgl"],"align":["right"]},{"label":["LINC01559"],"name":[79],"type":["lgl"],"align":["right"]},{"label":["LRRTM1"],"name":[80],"type":["lgl"],"align":["right"]},{"label":["MAGEC2"],"name":[81],"type":["lgl"],"align":["right"]},{"label":["MAGEC3"],"name":[82],"type":["lgl"],"align":["right"]},{"label":["MFI2"],"name":[83],"type":["lgl"],"align":["right"]},{"label":["MYH8"],"name":[84],"type":["lgl"],"align":["right"]},{"label":["NFE4"],"name":[85],"type":["lgl"],"align":["right"]},{"label":["NIPAL4"],"name":[86],"type":["lgl"],"align":["right"]},{"label":["NKX2_2"],"name":[87],"type":["lgl"],"align":["right"]},{"label":["NKX2_3"],"name":[88],"type":["lgl"],"align":["right"]},{"label":["NMRK2"],"name":[89],"type":["lgl"],"align":["right"]},{"label":["NUPR1L"],"name":[90],"type":["lgl"],"align":["right"]},{"label":["OTX1"],"name":[91],"type":["lgl"],"align":["right"]},{"label":["PADI3"],"name":[92],"type":["lgl"],"align":["right"]},{"label":["PAEP"],"name":[93],"type":["lgl"],"align":["right"]},{"label":["PI3"],"name":[94],"type":["lgl"],"align":["right"]},{"label":["PITX1"],"name":[95],"type":["lgl"],"align":["right"]},{"label":["PLG"],"name":[96],"type":["lgl"],"align":["right"]},{"label":["PRR15L"],"name":[97],"type":["lgl"],"align":["right"]},{"label":["PSG9"],"name":[98],"type":["lgl"],"align":["right"]},{"label":["PVALB"],"name":[99],"type":["lgl"],"align":["right"]},{"label":["RAB25"],"name":[100],"type":["lgl"],"align":["right"]},{"label":["RHBG"],"name":[101],"type":["lgl"],"align":["right"]},{"label":["RP11_10O22.1"],"name":[102],"type":["lgl"],"align":["right"]},{"label":["RP11_150O12.1"],"name":[103],"type":["lgl"],"align":["right"]},{"label":["RP11_161D15.1"],"name":[104],"type":["lgl"],"align":["right"]},{"label":["RP11_310H4.6"],"name":[105],"type":["lgl"],"align":["right"]},{"label":["RP11_314M24.1"],"name":[106],"type":["lgl"],"align":["right"]},{"label":["RP11_400N13.3"],"name":[107],"type":["lgl"],"align":["right"]},{"label":["RP11_425D17.2"],"name":[108],"type":["lgl"],"align":["right"]},{"label":["RP11_440G9.1"],"name":[109],"type":["lgl"],"align":["right"]},{"label":["RP11_54H7.4"],"name":[110],"type":["lgl"],"align":["right"]},{"label":["RP11_554D15.1"],"name":[111],"type":["lgl"],"align":["right"]},{"label":["RP11_586K2.1"],"name":[112],"type":["lgl"],"align":["right"]},{"label":["RP11_643A5.3"],"name":[113],"type":["lgl"],"align":["right"]},{"label":["RP11_690G19.4"],"name":[114],"type":["lgl"],"align":["right"]},{"label":["RP11_95M15.2"],"name":[115],"type":["lgl"],"align":["right"]},{"label":["RP13_895J2.6"],"name":[116],"type":["lgl"],"align":["right"]},{"label":["RP4_568C11.4"],"name":[117],"type":["lgl"],"align":["right"]},{"label":["RP5_984P4.6"],"name":[118],"type":["lgl"],"align":["right"]},{"label":["RTL1"],"name":[119],"type":["lgl"],"align":["right"]},{"label":["SAA1"],"name":[120],"type":["lgl"],"align":["right"]},{"label":["SAA2_SAA4"],"name":[121],"type":["lgl"],"align":["right"]},{"label":["SBSN"],"name":[122],"type":["lgl"],"align":["right"]},{"label":["SFTPB"],"name":[123],"type":["lgl"],"align":["right"]},{"label":["SHOX2"],"name":[124],"type":["lgl"],"align":["right"]},{"label":["SLC12A3"],"name":[125],"type":["lgl"],"align":["right"]},{"label":["SLC18A3"],"name":[126],"type":["lgl"],"align":["right"]},{"label":["SLC22A8"],"name":[127],"type":["lgl"],"align":["right"]},{"label":["SLC30A8"],"name":[128],"type":["lgl"],"align":["right"]},{"label":["SLC34A1"],"name":[129],"type":["lgl"],"align":["right"]},{"label":["SLC38A5"],"name":[130],"type":["lgl"],"align":["right"]},{"label":["SLC4A1"],"name":[131],"type":["lgl"],"align":["right"]},{"label":["SLC6A15"],"name":[132],"type":["lgl"],"align":["right"]},{"label":["SLC6A18"],"name":[133],"type":["lgl"],"align":["right"]},{"label":["SLC7A13"],"name":[134],"type":["lgl"],"align":["right"]},{"label":["STAC2"],"name":[135],"type":["lgl"],"align":["right"]},{"label":["TCEAL2"],"name":[136],"type":["lgl"],"align":["right"]},{"label":["TCL6"],"name":[137],"type":["lgl"],"align":["right"]},{"label":["TMEM213"],"name":[138],"type":["lgl"],"align":["right"]},{"label":["TMEM61"],"name":[139],"type":["lgl"],"align":["right"]},{"label":["TMPRSS11E"],"name":[140],"type":["lgl"],"align":["right"]},{"label":["TNNT1"],"name":[141],"type":["lgl"],"align":["right"]},{"label":["TTR"],"name":[142],"type":["lgl"],"align":["right"]},{"label":["TUBA3E"],"name":[143],"type":["lgl"],"align":["right"]},{"label":["TUBBP6"],"name":[144],"type":["lgl"],"align":["right"]},{"label":["TUNAR"],"name":[145],"type":["lgl"],"align":["right"]},{"label":["UGT1A10"],"name":[146],"type":["lgl"],"align":["right"]},{"label":["UMOD"],"name":[147],"type":["lgl"],"align":["right"]},{"label":["WFDC5"],"name":[148],"type":["lgl"],"align":["right"]},{"label":["ZIC2"],"name":[149],"type":["lgl"],"align":["right"]},{"label":["ZIC5"],"name":[150],"type":["lgl"],"align":["right"]},{"label":["features"],"name":[151],"type":["list"],"align":["right"]},{"label":["classif.bacc"],"name":[152],"type":["dbl"],"align":["right"]}],"data":[{"1":"TRUE","2":"TRUE","3":"TRUE","4":"TRUE","5":"TRUE","6":"TRUE","7":"TRUE","8":"TRUE","9":"TRUE","10":"TRUE","11":"TRUE","12":"TRUE","13":"TRUE","14":"TRUE","15":"TRUE","16":"TRUE","17":"TRUE","18":"TRUE","19":"TRUE","20":"TRUE","21":"TRUE","22":"TRUE","23":"TRUE","24":"TRUE","25":"TRUE","26":"TRUE","27":"TRUE","28":"TRUE","29":"TRUE","30":"TRUE","31":"TRUE","32":"TRUE","33":"TRUE","34":"TRUE","35":"TRUE","36":"TRUE","37":"TRUE","38":"TRUE","39":"TRUE","40":"TRUE","41":"TRUE","42":"TRUE","43":"TRUE","44":"TRUE","45":"TRUE","46":"TRUE","47":"TRUE","48":"TRUE","49":"TRUE","50":"TRUE","51":"TRUE","52":"TRUE","53":"TRUE","54":"TRUE","55":"TRUE","56":"TRUE","57":"TRUE","58":"TRUE","59":"TRUE","60":"TRUE","61":"TRUE","62":"TRUE","63":"TRUE","64":"TRUE","65":"TRUE","66":"TRUE","67":"TRUE","68":"TRUE","69":"TRUE","70":"TRUE","71":"TRUE","72":"TRUE","73":"TRUE","74":"TRUE","75":"TRUE","76":"TRUE","77":"TRUE","78":"TRUE","79":"TRUE","80":"TRUE","81":"TRUE","82":"TRUE","83":"TRUE","84":"TRUE","85":"TRUE","86":"TRUE","87":"TRUE","88":"TRUE","89":"TRUE","90":"TRUE","91":"TRUE","92":"TRUE","93":"TRUE","94":"TRUE","95":"TRUE","96":"TRUE","97":"TRUE","98":"TRUE","99":"TRUE","100":"TRUE","101":"TRUE","102":"TRUE","103":"TRUE","104":"TRUE","105":"TRUE","106":"TRUE","107":"TRUE","108":"TRUE","109":"TRUE","110":"TRUE","111":"TRUE","112":"TRUE","113":"TRUE","114":"TRUE","115":"TRUE","116":"TRUE","117":"TRUE","118":"TRUE","119":"TRUE","120":"TRUE","121":"TRUE","122":"TRUE","123":"TRUE","124":"TRUE","125":"TRUE","126":"TRUE","127":"TRUE","128":"TRUE","129":"TRUE","130":"TRUE","131":"TRUE","132":"TRUE","133":"TRUE","134":"TRUE","135":"TRUE","136":"TRUE","137":"TRUE","138":"TRUE","139":"TRUE","140":"TRUE","141":"TRUE","142":"TRUE","143":"TRUE","144":"TRUE","145":"TRUE","146":"TRUE","147":"TRUE","148":"TRUE","149":"TRUE","150":"TRUE","151":"<chr [150]>","152":"0.6936458"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
print(instance$result_y)
```

    ## classif.bacc 
    ##    0.6936458

``` r
tsk_rfe <- TaskClassif$new(id="kirc_rfe", 
                               backend = kirc_norm, 
                               target = "metastasis", 
                               positive = "M1")


tsk_rfe$select(cols = instance$result_feature_set)
tsk_rfe
```

    ## <TaskClassif:kirc_rfe> (277 x 151)
    ## * Target: metastasis
    ## * Properties: twoclass
    ## * Features (150):
    ##   - dbl (150): AC003092.1, AC006262.4, AC006262.5, AC007879.6,
    ##     AC104654.2, AC116614.1, AFM, AHSG, AMH, APCDD1L_AS1, AQP2, AQP6,
    ##     ATP6V0A4, ATP6V0D2, BARX1, BSND, C10orf99, C14orf180, CA1, CASP14,
    ##     CCNA1, CDC42P2, CHAT, CIDEC, CILP2, CITF22_24E5.1, CLCNKB, CLDN8,
    ##     CLMP, COL11A1, COL7A1, CPNE7, CTD_2008P7.9, CXCL13, CYP1A1, DMRT2,
    ##     DQX1, EN2, ESRP1, FAM83B, FDCSP, FGF5, FKBP9P1, FOXI2, FXYD4, GGT6,
    ##     GLB1L3, GOLGA6L2, GOLGA6L7P, GPR110, HEPACAM2, HHATL, HMGA2,
    ##     HS3ST3A1, IGF2BP3, IGFBP1, IGFL1P1, IGFN1, IGHV1_69, IGKV3_11,
    ##     IGLC7, IGLV3_19, INHBE, ITPKA, KCNJ1, KIRREL3, KLF17, KLK1, KLK15,
    ##     KNG1, KRT7, L1CAM, LECT1, LINC00890, LINC00942, LINC00973,
    ##     LINC01187, LINC01436, LINC01559, LRRTM1, MAGEC2, MAGEC3, MFI2,
    ##     MYH8, NFE4, NIPAL4, NKX2_2, NKX2_3, NMRK2, NUPR1L, OTX1, PADI3,
    ##     PAEP, PI3, PITX1, PLG, PRR15L, PSG9, PVALB, RAB25, [...]

``` r
lrn.rpa = lrn("classif.rpart")
lrn.xgb = lrn("classif.xgboost")
lrn.rgn = lrn("classif.ranger")
lrn.svm = lrn("classif.ksvm")

grid = benchmark_grid(
  task = list(tsk_filt, tsk_rfe),
  learner = list(lrn.rpa, lrn.xgb, lrn.rgn, lrn.svm),
  resampling = rsmp("cv", folds = 3)
)
```

``` r
bmr = benchmark(grid, store_models = TRUE)
```

    ## INFO  [15:59:38.984] [mlr3]  Running benchmark with 24 resampling iterations 
    ## INFO  [15:59:38.992] [mlr3]  Applying learner 'classif.ranger' on task 'filt_rpart' (iter 1/3) 
    ## INFO  [15:59:39.094] [mlr3]  Applying learner 'classif.ksvm' on task 'kirc_rfe' (iter 3/3) 
    ## INFO  [15:59:39.387] [mlr3]  Applying learner 'classif.xgboost' on task 'kirc_rfe' (iter 1/3) 
    ## [15:59:39] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:39.528] [mlr3]  Applying learner 'classif.rpart' on task 'filt_rpart' (iter 2/3) 
    ## INFO  [15:59:39.546] [mlr3]  Applying learner 'classif.xgboost' on task 'filt_rpart' (iter 3/3) 
    ## [15:59:39] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:39.570] [mlr3]  Applying learner 'classif.ranger' on task 'kirc_rfe' (iter 2/3) 
    ## INFO  [15:59:39.709] [mlr3]  Applying learner 'classif.xgboost' on task 'kirc_rfe' (iter 2/3) 
    ## [15:59:39] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:39.732] [mlr3]  Applying learner 'classif.ranger' on task 'kirc_rfe' (iter 3/3) 
    ## INFO  [15:59:39.868] [mlr3]  Applying learner 'classif.rpart' on task 'filt_rpart' (iter 1/3) 
    ## INFO  [15:59:39.884] [mlr3]  Applying learner 'classif.xgboost' on task 'filt_rpart' (iter 2/3) 
    ## [15:59:39] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:39.900] [mlr3]  Applying learner 'classif.ksvm' on task 'kirc_rfe' (iter 1/3) 
    ## INFO  [15:59:40.486] [mlr3]  Applying learner 'classif.rpart' on task 'filt_rpart' (iter 3/3) 
    ## INFO  [15:59:40.512] [mlr3]  Applying learner 'classif.ranger' on task 'kirc_rfe' (iter 1/3) 
    ## INFO  [15:59:40.669] [mlr3]  Applying learner 'classif.xgboost' on task 'filt_rpart' (iter 1/3) 
    ## [15:59:40] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:40.686] [mlr3]  Applying learner 'classif.rpart' on task 'kirc_rfe' (iter 3/3) 
    ## INFO  [15:59:40.744] [mlr3]  Applying learner 'classif.ksvm' on task 'filt_rpart' (iter 2/3) 
    ## INFO  [15:59:40.786] [mlr3]  Applying learner 'classif.ranger' on task 'filt_rpart' (iter 3/3) 
    ## INFO  [15:59:40.908] [mlr3]  Applying learner 'classif.rpart' on task 'kirc_rfe' (iter 1/3) 
    ## INFO  [15:59:40.968] [mlr3]  Applying learner 'classif.rpart' on task 'kirc_rfe' (iter 2/3) 
    ## INFO  [15:59:41.020] [mlr3]  Applying learner 'classif.ksvm' on task 'kirc_rfe' (iter 2/3) 
    ## INFO  [15:59:41.096] [mlr3]  Applying learner 'classif.ksvm' on task 'filt_rpart' (iter 3/3) 
    ## INFO  [15:59:41.127] [mlr3]  Applying learner 'classif.ranger' on task 'filt_rpart' (iter 2/3) 
    ## INFO  [15:59:41.248] [mlr3]  Applying learner 'classif.xgboost' on task 'kirc_rfe' (iter 3/3) 
    ## [15:59:41] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:41.272] [mlr3]  Applying learner 'classif.ksvm' on task 'filt_rpart' (iter 1/3) 
    ## INFO  [15:59:41.320] [mlr3]  Finished benchmark

``` r
cols <- c("task_id", "learner_id",  "classif.bacc", "classif.sensitivity", "classif.specificity")

measures <- list(
  msr("classif.bacc"),
  msr("classif.sensitivity"), 
  msr("classif.specificity")
)

bmr_df <- bmr$aggregate(measures) %>%
  dplyr::select(cols) %>%
  dplyr::arrange(desc(classif.bacc))
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(cols)` instead of `cols` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

``` r
bmr_df
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["task_id"],"name":[1],"type":["chr"],"align":["left"]},{"label":["learner_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["classif.bacc"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["classif.sensitivity"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["classif.specificity"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"filt_rpart","2":"classif.ksvm","3":"0.7555397","4":"0.5714450","5":"0.9396344"},{"1":"filt_rpart","2":"classif.ranger","3":"0.7489372","4":"0.5681006","5":"0.9297739"},{"1":"kirc_rfe","2":"classif.ranger","3":"0.7168554","4":"0.5124459","5":"0.9212648"},{"1":"kirc_rfe","2":"classif.ksvm","3":"0.7083111","4":"0.4763709","5":"0.9402513"},{"1":"filt_rpart","2":"classif.xgboost","3":"0.6606630","4":"0.5113213","5":"0.8100048"},{"1":"kirc_rfe","2":"classif.xgboost","3":"0.6576482","4":"0.4983766","5":"0.8169198"},{"1":"filt_rpart","2":"classif.rpart","3":"0.6513989","4":"0.4623458","5":"0.8404521"},{"1":"kirc_rfe","2":"classif.rpart","3":"0.6093991","4":"0.4148629","5":"0.8039354"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

# Plotting Benchmark Results

``` r
library("mlr3viz")
library("ggplot2")

autoplot(bmr, measure= msr("classif.bacc")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](figs/tutorial_unnamed-chunk-24-1.png)<!-- -->

# Automating the Feature Selection

``` r
library("paradox")
library("mlr3fselect")
library("mlr3extralearners")

terminator = trm("evals", n_evals = 20)
fselector = fs("random_search")

lrn.rpa = lrn("classif.rpart")
lrn.xgb = lrn("classif.xgboost")
lrn.ranger = lrn("classif.ranger")
lrn.ksvm = lrn("classif.ksvm")
lrn.svm = lrn("classif.svm")


at.rpa = AutoFSelector$new(
  learner = lrn.rpa,
  resampling = rsmp("holdout"),
  measure = msr("classif.bacc"),
  terminator = terminator,
  fselector = fselector
)

at.xgb = AutoFSelector$new(
  learner = lrn.xgb,
  resampling = rsmp("holdout"),
  measure = msr("classif.bacc"),
  terminator = terminator,
  fselector = fselector
)

at.ranger = AutoFSelector$new(
  learner = lrn.ranger,
  resampling = rsmp("holdout"),
  measure = msr("classif.bacc"),
  terminator = terminator,
  fselector = fselector
)

at.svm = AutoFSelector$new(
  learner = lrn.svm,
  resampling = rsmp("holdout"),
  measure = msr("classif.bacc"),
  terminator = terminator,
  fselector = fselector
)

at.ksvm = AutoFSelector$new(
  learner = lrn.ksvm,
  resampling = rsmp("holdout"),
  measure = msr("classif.bacc"),
  terminator = terminator,
  fselector = fselector
)
```

``` r
tsk_cla <- TaskClassif$new(id="kirc_cla", 
                               backend = kirc_norm, 
                               target = "metastasis", 
                               positive = "M1")


grid = benchmark_grid(
  task = tsk_cla,
  learner = list(at.rpa, at.xgb, at.ranger, at.svm, at.ksvm),
  resampling = rsmp("cv", folds = 3)
)
```

``` r
bmr = benchmark(grid, store_models = TRUE)
```

    ## INFO  [15:59:43.541] [mlr3]  Running benchmark with 15 resampling iterations 
    ## INFO  [15:59:43.547] [mlr3]  Applying learner 'classif.xgboost.fselector' on task 'kirc_cla' (iter 2/3) 
    ## INFO  [15:59:43.741] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [15:59:43.745] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [15:59:44.975] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [15:59:44.982] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.045] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.112] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.175] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.239] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.314] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.375] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.433] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.496] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.566] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:45] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:45.630] [mlr3]  Finished benchmark 
    ## INFO  [15:59:46.113] [bbotk] Result of batch 1: 
    ## INFO  [15:59:46.130] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]        TRUE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE       TRUE       TRUE      FALSE       TRUE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:46.130] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE       FALSE FALSE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE    FALSE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE    FALSE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE        TRUE  TRUE FALSE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE FALSE  TRUE         FALSE   TRUE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE   TRUE FALSE    TRUE  TRUE FALSE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE  TRUE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  FALSE FALSE    TRUE FALSE FALSE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:46.130] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE   TRUE FALSE         TRUE  FALSE  FALSE FALSE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE   TRUE FALSE        FALSE   TRUE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE    TRUE  TRUE FALSE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:46.130] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]     FALSE FALSE FALSE    FALSE    TRUE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE  FALSE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:46.130] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [15:59:46.130] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:46.130] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:46.130] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE     TRUE FALSE FALSE FALSE   FALSE  TRUE  TRUE FALSE FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     TRUE FALSE FALSE FALSE    TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE     TRUE FALSE FALSE  TRUE   FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE  TRUE FALSE FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:46.130] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:46.130] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE  FALSE   TRUE   TRUE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE   TRUE  FALSE  FALSE FALSE  FALSE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE   TRUE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE FALSE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE FALSE FALSE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE  TRUE  TRUE  FALSE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE FALSE FALSE   TRUE FALSE  TRUE FALSE FALSE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:46.130] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:46.130] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:46.130] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE        FALSE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE        FALSE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [15:59:46.130] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:46.130] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE        FALSE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [15:59:46.130] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:46.130] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:46.130] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:46.130] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:46.130] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [15:59:46.130] [bbotk]      FALSE FALSE  TRUE  TRUE   FALSE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [15:59:46.130] [bbotk]      FALSE  TRUE  TRUE FALSE   FALSE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [15:59:46.130] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:46.130] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:46.130] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [15:59:46.130] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:46.130] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:46.130] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:46.130] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE   FALSE    TRUE    TRUE FALSE  FALSE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE   FALSE   FALSE    TRUE  TRUE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE   FALSE   FALSE   FALSE  TRUE  FALSE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE    TRUE   FALSE   FALSE FALSE  FALSE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE   FALSE   FALSE   FALSE  TRUE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [15:59:46.130] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5342625 
    ## INFO  [15:59:46.130] [bbotk]  FALSE  TRUE  FALSE   TRUE FALSE    TRUE FALSE  TRUE FALSE  TRUE    0.6649245 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE    TRUE FALSE FALSE FALSE FALSE    0.5807201 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE FALSE  TRUE FALSE    0.6887340 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE FALSE  TRUE    0.6277584 
    ## INFO  [15:59:46.130] [bbotk]   TRUE FALSE  FALSE   TRUE FALSE   FALSE FALSE  TRUE  TRUE  TRUE    0.6411150 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6399535 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE FALSE  TRUE FALSE    0.5336818 
    ## INFO  [15:59:46.130] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5923345 
    ## INFO  [15:59:46.130] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6649245 
    ## INFO  [15:59:46.130] [bbotk]                                 uhash 
    ## INFO  [15:59:46.130] [bbotk]  6a9b72fb-0390-44e4-88c4-c567782eb930 
    ## INFO  [15:59:46.130] [bbotk]  73fb9f41-9805-4818-8026-4b6495301a30 
    ## INFO  [15:59:46.130] [bbotk]  d4a2d5ba-d689-4b6d-807e-7d0f227da504 
    ## INFO  [15:59:46.130] [bbotk]  b8f17571-f2a1-404f-b498-11fdcac58878 
    ## INFO  [15:59:46.130] [bbotk]  c7e18ec5-ffc1-450d-b1d7-9b2ac30db675 
    ## INFO  [15:59:46.130] [bbotk]  2af49c0b-d761-4e86-8271-d22eea29b545 
    ## INFO  [15:59:46.130] [bbotk]  cc2b1aa1-99ec-4bc8-87c8-033ca8334bbd 
    ## INFO  [15:59:46.130] [bbotk]  b848cd5b-cc5f-40ec-a74b-d835be41e10f 
    ## INFO  [15:59:46.130] [bbotk]  d34a2999-44fd-4910-b5d7-79dd9a3669f8 
    ## INFO  [15:59:46.130] [bbotk]  a3c00ba6-adb1-443c-8d74-a35fdfa1de35 
    ## INFO  [15:59:46.133] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [15:59:47.347] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [15:59:47.358] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.425] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.490] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.550] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.612] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.695] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.760] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.822] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:47.898] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:47] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:48.003] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [15:59:48] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:48.074] [mlr3]  Finished benchmark 
    ## INFO  [15:59:48.603] [bbotk] Result of batch 2: 
    ## INFO  [15:59:48.614] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:48.614] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]       FALSE       TRUE      FALSE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]        TRUE      FALSE      FALSE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE       TRUE  TRUE FALSE 
    ## INFO  [15:59:48.614] [bbotk]        TRUE       TRUE      FALSE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:48.614] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE        TRUE FALSE FALSE    FALSE     TRUE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE        TRUE FALSE  TRUE    FALSE    FALSE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE       FALSE FALSE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:48.614] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE   TRUE FALSE    TRUE  TRUE FALSE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE FALSE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE FALSE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  FALSE FALSE   FALSE FALSE  TRUE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE  TRUE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:48.614] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE FALSE FALSE FALSE  TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE   TRUE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]     TRUE  FALSE FALSE        FALSE   TRUE  FALSE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE FALSE FALSE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE   TRUE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE    TRUE FALSE FALSE  TRUE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE   FALSE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE    TRUE FALSE FALSE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:48.614] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [15:59:48.614] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:48.614] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:48.614] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [15:59:48.614] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:48.614] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:48.614] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [15:59:48.614] [bbotk]      TRUE FALSE FALSE    FALSE    TRUE   TRUE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [15:59:48.614] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:48.614] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:48.614] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE    FALSE FALSE  TRUE  TRUE    TRUE FALSE  TRUE FALSE FALSE  TRUE FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE    FALSE  TRUE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE    FALSE FALSE  TRUE  TRUE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:48.614] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE      TRUE     FALSE      TRUE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     FALSE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE     FALSE      TRUE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:48.614] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:48.614] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE   TRUE FALSE FALSE  TRUE  FALSE  FALSE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE  FALSE  TRUE FALSE  TRUE  FALSE  FALSE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE  FALSE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE   TRUE FALSE FALSE  TRUE  FALSE  FALSE   TRUE  TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE  TRUE FALSE FALSE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE FALSE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:48.614] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:48.614] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:48.614] [bbotk]           TRUE        FALSE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE        FALSE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:48.614] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE        FALSE         TRUE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE        FALSE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE        FALSE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [15:59:48.614] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:48.614] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:48.614] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE   FALSE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [15:59:48.614] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [15:59:48.614] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE    TRUE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:48.614] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [15:59:48.614] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:48.614] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE    TRUE    TRUE   FALSE FALSE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE    TRUE   FALSE    TRUE FALSE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE   FALSE   FALSE   FALSE FALSE   TRUE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]    TRUE    TRUE   FALSE   FALSE  TRUE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:48.614] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE FALSE    0.6411150 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6771196 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE FALSE    0.7253194 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE    0.6765389 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6277584 
    ## INFO  [15:59:48.614] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE FALSE FALSE  TRUE FALSE    0.6416957 
    ## INFO  [15:59:48.614] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE  TRUE    0.6167247 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.5923345 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE    TRUE FALSE FALSE FALSE FALSE    0.5917538 
    ## INFO  [15:59:48.614] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5447154 
    ## INFO  [15:59:48.614] [bbotk]                                 uhash 
    ## INFO  [15:59:48.614] [bbotk]  eb964b6e-1f71-4a55-be07-62baec0fb777 
    ## INFO  [15:59:48.614] [bbotk]  f9082f78-9213-4f6e-a4db-fa8247f5269b 
    ## INFO  [15:59:48.614] [bbotk]  395fcae1-f072-4bfe-9961-b3d618902836 
    ## INFO  [15:59:48.614] [bbotk]  fa1c12f6-9ff7-456c-bf20-1b47408ba1ef 
    ## INFO  [15:59:48.614] [bbotk]  604c5b78-59ae-4d9e-94b0-d9ad90ba29ed 
    ## INFO  [15:59:48.614] [bbotk]  062e24eb-ef84-45e7-be71-09886a70e7ad 
    ## INFO  [15:59:48.614] [bbotk]  3e64da8a-1b18-4ea6-b7f5-5cd8f4ac60a6 
    ## INFO  [15:59:48.614] [bbotk]  179d9b3a-80c2-4ea7-ae48-2a02cc18699b 
    ## INFO  [15:59:48.614] [bbotk]  4f6d54d4-1631-447a-b9b5-52b70410aa7a 
    ## INFO  [15:59:48.614] [bbotk]  f95d17bc-e3e2-49ff-a4ed-e850670cdd29 
    ## INFO  [15:59:48.621] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [15:59:48.622] [bbotk] Result: 
    ## INFO  [15:59:48.628] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:48.628] [bbotk]       FALSE       TRUE      FALSE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [15:59:48.628] [bbotk]   AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:48.628] [bbotk]  TRUE        TRUE FALSE FALSE    FALSE     TRUE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [15:59:48.628] [bbotk]    CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:48.628] [bbotk]  FALSE   TRUE FALSE    TRUE TRUE FALSE FALSE          TRUE   TRUE FALSE TRUE 
    ## INFO  [15:59:48.628] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:48.628] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE   TRUE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [15:59:48.628] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:48.628] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [15:59:48.628] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:48.628] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:48.628] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15  KNG1 KRT7 L1CAM 
    ## INFO  [15:59:48.628] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE    TRUE  TRUE TRUE  TRUE FALSE TRUE FALSE 
    ## INFO  [15:59:48.628] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:48.628] [bbotk]   TRUE      TRUE     FALSE      TRUE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:48.628] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 
    ## INFO  [15:59:48.628] [bbotk]    TRUE   TRUE FALSE FALSE TRUE  FALSE  FALSE   TRUE FALSE   TRUE TRUE  TRUE 
    ## INFO  [15:59:48.628] [bbotk]  PAEP   PI3 PITX1  PLG PRR15L PSG9 PVALB RAB25  RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [15:59:48.628] [bbotk]  TRUE FALSE FALSE TRUE   TRUE TRUE  TRUE FALSE FALSE         TRUE         FALSE 
    ## INFO  [15:59:48.628] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [15:59:48.628] [bbotk]           TRUE         TRUE         FALSE          TRUE          TRUE 
    ## INFO  [15:59:48.628] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [15:59:48.628] [bbotk]          TRUE       FALSE          TRUE        FALSE         TRUE          TRUE 
    ## INFO  [15:59:48.628] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [15:59:48.628] [bbotk]         FALSE        FALSE        FALSE        TRUE FALSE TRUE     FALSE TRUE 
    ## INFO  [15:59:48.628] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [15:59:48.628] [bbotk]  FALSE FALSE   FALSE   FALSE   FALSE    TRUE    TRUE    TRUE  FALSE    TRUE 
    ## INFO  [15:59:48.628] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1   TTR TUBA3E 
    ## INFO  [15:59:48.628] [bbotk]     TRUE   FALSE FALSE  FALSE TRUE   FALSE   TRUE     FALSE  TRUE FALSE  FALSE 
    ## INFO  [15:59:48.628] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5  ZIC2  ZIC5 
    ## INFO  [15:59:48.628] [bbotk]    TRUE  TRUE   FALSE TRUE  TRUE FALSE FALSE 
    ## INFO  [15:59:48.628] [bbotk]                                                         features classif.bacc 
    ## INFO  [15:59:48.628] [bbotk]  AC006262.4,AC007879.6,AC104654.2,AC116614.1,AMH,APCDD1L_AS1,...    0.7253194 
    ## [15:59:48] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [15:59:48.706] [mlr3]  Applying learner 'classif.ranger.fselector' on task 'kirc_cla' (iter 3/3) 
    ## INFO  [15:59:48.904] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [15:59:48.906] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [15:59:49.806] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [15:59:49.814] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:49.943] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.064] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.190] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.360] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.490] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.636] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.768] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:50.931] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:51.066] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:51.210] [mlr3]  Finished benchmark 
    ## INFO  [15:59:51.634] [bbotk] Result of batch 1: 
    ## INFO  [15:59:51.645] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:51.645] [bbotk]        TRUE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]        TRUE       TRUE       TRUE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]       FALSE       TRUE       TRUE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:51.645] [bbotk]   TRUE        TRUE  TRUE FALSE    FALSE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE        TRUE FALSE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE       FALSE  TRUE  TRUE    FALSE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE        TRUE  TRUE  TRUE    FALSE    FALSE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE FALSE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE         FALSE   TRUE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:51.645] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]     TRUE   TRUE FALSE        FALSE   TRUE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]     TRUE   TRUE FALSE        FALSE   TRUE   TRUE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE   TRUE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:51.645] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE   FALSE FALSE  TRUE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE    TRUE FALSE FALSE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:51.645] [bbotk]     FALSE FALSE  TRUE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:51.645] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:51.645] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:51.645] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:51.645] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [15:59:51.645] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:51.645] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:51.645] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:51.645] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [15:59:51.645] [bbotk]     FALSE  TRUE  TRUE    FALSE   FALSE  FALSE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [15:59:51.645] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:51.645] [bbotk]  FALSE    FALSE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE    FALSE FALSE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE FALSE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE     TRUE FALSE FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE     TRUE  TRUE  TRUE FALSE    TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE   FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:51.645] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE     FALSE      TRUE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE      TRUE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [15:59:51.645] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:51.645] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE   TRUE FALSE FALSE FALSE   TRUE   TRUE   TRUE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE  FALSE FALSE FALSE FALSE  FALSE   TRUE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE  FALSE  FALSE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE  FALSE  TRUE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE  FALSE FALSE  TRUE  TRUE  FALSE   TRUE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE FALSE  TRUE  FALSE FALSE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE  TRUE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE FALSE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE FALSE FALSE  FALSE  TRUE  TRUE FALSE FALSE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:51.645] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [15:59:51.645] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [15:59:51.645] [bbotk]          FALSE         TRUE         TRUE        FALSE       FALSE  TRUE  TRUE 
    ## INFO  [15:59:51.645] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [15:59:51.645] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:51.645] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:51.645] [bbotk]      FALSE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:51.645] [bbotk]       TRUE FALSE FALSE  TRUE   FALSE   FALSE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [15:59:51.645] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:51.645] [bbotk]       TRUE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:51.645] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:51.645] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:51.645] [bbotk]       TRUE  TRUE FALSE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:51.645] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:51.645] [bbotk]      FALSE FALSE FALSE FALSE   FALSE    TRUE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:51.645] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:51.645] [bbotk]    TRUE   FALSE    TRUE   FALSE FALSE  FALSE FALSE    TRUE   TRUE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:51.645] [bbotk]    TRUE   FALSE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE  FALSE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]   FALSE    TRUE   FALSE    TRUE FALSE  FALSE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:51.645] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7130719 
    ## INFO  [15:59:51.645] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE FALSE FALSE  TRUE FALSE    0.6908497 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE  FALSE   TRUE FALSE   FALSE FALSE FALSE  TRUE  TRUE    0.6725490 
    ## INFO  [15:59:51.645] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE  TRUE    0.6725490 
    ## INFO  [15:59:51.645] [bbotk]   TRUE FALSE  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE  TRUE    0.6098039 
    ## INFO  [15:59:51.645] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.7313725 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5176471 
    ## INFO  [15:59:51.645] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7019608 
    ## INFO  [15:59:51.645] [bbotk]  FALSE  TRUE  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE  TRUE    0.6980392 
    ## INFO  [15:59:51.645] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6725490 
    ## INFO  [15:59:51.645] [bbotk]                                 uhash 
    ## INFO  [15:59:51.645] [bbotk]  d14664ef-a75f-43a2-9a0a-7ac5d9babbcd 
    ## INFO  [15:59:51.645] [bbotk]  71e3819f-aafd-428d-b4ba-e913e0158024 
    ## INFO  [15:59:51.645] [bbotk]  9551c531-e056-41da-986f-e253a7668cbb 
    ## INFO  [15:59:51.645] [bbotk]  a39947fa-15a0-4437-8051-1e1989efff42 
    ## INFO  [15:59:51.645] [bbotk]  3d30be39-151b-4e4e-8a0d-7b6a76150717 
    ## INFO  [15:59:51.645] [bbotk]  80856b54-94d0-4d84-89c0-49ddcbac60ba 
    ## INFO  [15:59:51.645] [bbotk]  c67f8726-59ee-479b-a132-06ff034a3c43 
    ## INFO  [15:59:51.645] [bbotk]  552d6485-1cc3-41b7-ad23-c48642eb46f0 
    ## INFO  [15:59:51.645] [bbotk]  19b32241-4d44-469c-9d9d-fba6be33834a 
    ## INFO  [15:59:51.645] [bbotk]  3b6ab7df-e102-4182-b209-29e36362c513 
    ## INFO  [15:59:51.648] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [15:59:52.498] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [15:59:52.511] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:52.689] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:52.814] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:52.941] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:53.089] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:53.303] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:53.456] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:53.597] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:53.737] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:53.898] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:54.061] [mlr3]  Finished benchmark 
    ## INFO  [15:59:54.527] [bbotk] Result of batch 2: 
    ## INFO  [15:59:54.546] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:54.546] [bbotk]        TRUE      FALSE      FALSE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]       FALSE       TRUE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:54.546] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE    FALSE FALSE FALSE     TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE    FALSE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:54.546] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE   TRUE FALSE    TRUE FALSE  TRUE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  FALSE FALSE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:54.546] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE  FALSE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  TRUE   FALSE FALSE  TRUE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:54.546] [bbotk]      TRUE FALSE  TRUE    FALSE   FALSE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [15:59:54.546] [bbotk]     FALSE FALSE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [15:59:54.546] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:54.546] [bbotk]      TRUE  TRUE FALSE     TRUE   FALSE   TRUE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [15:59:54.546] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:54.546] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE     TRUE FALSE  TRUE FALSE    TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE     TRUE  TRUE FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE     TRUE FALSE FALSE  TRUE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:54.546] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE     FALSE     FALSE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE      TRUE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:54.546] [bbotk]  MAGEC2 MAGEC3 MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:54.546] [bbotk]   FALSE   TRUE TRUE  TRUE  TRUE   TRUE  FALSE   TRUE FALSE  FALSE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE   TRUE TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE   TRUE TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE  FALSE TRUE FALSE  TRUE  FALSE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE   TRUE TRUE FALSE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE  FALSE TRUE  TRUE  TRUE  FALSE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE   TRUE TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE  FALSE TRUE FALSE  TRUE  FALSE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE   TRUE TRUE  TRUE FALSE   TRUE  FALSE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE  FALSE TRUE  TRUE FALSE  FALSE   TRUE   TRUE FALSE   TRUE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE FALSE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE  TRUE FALSE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]   TRUE FALSE  TRUE FALSE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:54.546] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [15:59:54.546] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:54.546] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:54.546] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:54.546] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE         TRUE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE         TRUE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:54.546] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE  TRUE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [15:59:54.546] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:54.546] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE FALSE FALSE    TRUE   FALSE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [15:59:54.546] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:54.546] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:54.546] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE   TRUE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]   FALSE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [15:59:54.546] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:54.546] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.6908497 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6908497 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE   FALSE FALSE  TRUE  TRUE FALSE    0.6614379 
    ## INFO  [15:59:54.546] [bbotk]  FALSE  TRUE   TRUE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5915033 
    ## INFO  [15:59:54.546] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE FALSE    0.6503268 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE FALSE FALSE  TRUE FALSE    0.6725490 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE FALSE  TRUE    0.6908497 
    ## INFO  [15:59:54.546] [bbotk]   TRUE FALSE  FALSE  FALSE  TRUE    TRUE FALSE FALSE FALSE FALSE    0.6797386 
    ## INFO  [15:59:54.546] [bbotk]  FALSE FALSE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.7202614 
    ## INFO  [15:59:54.546] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE    0.6908497 
    ## INFO  [15:59:54.546] [bbotk]                                 uhash 
    ## INFO  [15:59:54.546] [bbotk]  a6d73ced-6163-4d8f-addd-0afa65fce0a4 
    ## INFO  [15:59:54.546] [bbotk]  f96bd115-1518-4bb6-8d16-34e467d60d1e 
    ## INFO  [15:59:54.546] [bbotk]  603f9e32-989b-4a43-8202-b8ba41894bc1 
    ## INFO  [15:59:54.546] [bbotk]  8c374bc8-8075-4105-be9a-4cb1a531c08b 
    ## INFO  [15:59:54.546] [bbotk]  f598ab9b-8f80-4f32-87c1-da0999f7fcda 
    ## INFO  [15:59:54.546] [bbotk]  dc3ceb59-def9-452c-a73d-f95ebef05105 
    ## INFO  [15:59:54.546] [bbotk]  d4e42a0f-89ea-4207-b041-31cb1bf283b6 
    ## INFO  [15:59:54.546] [bbotk]  133b5703-fc2b-421a-a216-fe14a86eddef 
    ## INFO  [15:59:54.546] [bbotk]  2c20149a-7da5-47c9-975a-9dc0c57ee56b 
    ## INFO  [15:59:54.546] [bbotk]  97564b4c-79fd-4685-8920-b48259f98747 
    ## INFO  [15:59:54.557] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [15:59:54.558] [bbotk] Result: 
    ## INFO  [15:59:54.568] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:54.568] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:54.568] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:54.568] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [15:59:54.568] [bbotk]    CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:54.568] [bbotk]  FALSE  FALSE  TRUE   FALSE TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:54.568] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:54.568] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:54.568] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:54.568] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:54.568] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:54.568] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:54.568] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1 KRT7 L1CAM 
    ## INFO  [15:59:54.568] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE TRUE  TRUE 
    ## INFO  [15:59:54.568] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:54.568] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [15:59:54.568] [bbotk]  MAGEC2 MAGEC3  MFI2 MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:54.568] [bbotk]   FALSE   TRUE FALSE TRUE FALSE  FALSE  FALSE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [15:59:54.568] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:54.568] [bbotk]  FALSE FALSE  TRUE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:54.568] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:54.568] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:54.568] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:54.568] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:54.568] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:54.568] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:54.568] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:54.568] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:54.568] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:54.568] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:54.568] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2  ZIC5 
    ## INFO  [15:59:54.568] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE TRUE FALSE 
    ## INFO  [15:59:54.568] [bbotk]                                           features classif.bacc 
    ## INFO  [15:59:54.568] [bbotk]  APCDD1L_AS1,ATP6V0D2,C10orf99,CCNA1,CHAT,DQX1,...    0.7313725 
    ## INFO  [15:59:54.701] [mlr3]  Applying learner 'classif.rpart.fselector' on task 'kirc_cla' (iter 2/3) 
    ## INFO  [15:59:54.909] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [15:59:54.912] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [15:59:55.483] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [15:59:55.492] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:55.577] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:55.658] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:55.733] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:55.800] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:55.888] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:55.986] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:56.061] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:56.171] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:56.326] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:56.401] [mlr3]  Finished benchmark 
    ## INFO  [15:59:56.772] [bbotk] Result of batch 1: 
    ## INFO  [15:59:56.784] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:56.784] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]        TRUE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]       FALSE      FALSE      FALSE       TRUE       TRUE      FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]       FALSE      FALSE      FALSE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:56.784] [bbotk]  FALSE        TRUE FALSE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE    FALSE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE       FALSE  TRUE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE  TRUE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  FALSE  TRUE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:56.784] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE  FALSE FALSE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE   TRUE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE  TRUE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]     TRUE  FALSE  TRUE         TRUE  FALSE   TRUE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE FALSE    TRUE FALSE FALSE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE FALSE   FALSE FALSE  TRUE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE FALSE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:56.784] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE  FALSE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [15:59:56.784] [bbotk]     FALSE  TRUE FALSE     TRUE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]      TRUE FALSE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE   TRUE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [15:59:56.784] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:56.784] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:56.784] [bbotk]  FALSE     TRUE  TRUE  TRUE FALSE    TRUE  TRUE FALSE  TRUE FALSE FALSE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE     TRUE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE    FALSE FALSE  TRUE FALSE    TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE   FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE     TRUE FALSE FALSE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE     TRUE FALSE  TRUE FALSE    TRUE  TRUE FALSE FALSE  TRUE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:56.784] [bbotk]   TRUE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE      TRUE      TRUE      TRUE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE      TRUE      TRUE     FALSE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE      TRUE     FALSE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:56.784] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:56.784] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE  FALSE   TRUE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE  FALSE FALSE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE  FALSE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE  FALSE FALSE FALSE FALSE  FALSE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE  FALSE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE   TRUE  TRUE  TRUE FALSE   TRUE   TRUE  FALSE FALSE   TRUE FALSE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE   TRUE FALSE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE  TRUE FALSE FALSE FALSE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE FALSE FALSE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  TRUE FALSE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:56.784] [bbotk]          FALSE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:56.784] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:56.784] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE         TRUE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:56.784] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:56.784] [bbotk]           TRUE        FALSE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE         TRUE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [15:59:56.784] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:56.784] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:56.784] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [15:59:56.784] [bbotk]      FALSE FALSE FALSE FALSE   FALSE    TRUE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:56.784] [bbotk]      FALSE  TRUE FALSE FALSE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:56.784] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [15:59:56.784] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:56.784] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:56.784] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:56.784] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [15:59:56.784] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:56.784] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:56.784] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:56.784] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE    TRUE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE    TRUE   FALSE    TRUE FALSE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [15:59:56.784] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE   TRUE   TRUE FALSE    TRUE  TRUE FALSE FALSE  TRUE    0.6413690 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE    TRUE  TRUE FALSE FALSE FALSE    0.7172619 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.7127976 
    ## INFO  [15:59:56.784] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6056548 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE    0.5833333 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6666667 
    ## INFO  [15:59:56.784] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.5967262 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7276786 
    ## INFO  [15:59:56.784] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE FALSE    0.6056548 
    ## INFO  [15:59:56.784] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE    0.6056548 
    ## INFO  [15:59:56.784] [bbotk]                                 uhash 
    ## INFO  [15:59:56.784] [bbotk]  6eddce89-97b3-4e65-99d2-8e05a6be4663 
    ## INFO  [15:59:56.784] [bbotk]  ba86ecac-bdca-45b8-bf27-73b1733594a3 
    ## INFO  [15:59:56.784] [bbotk]  09edf4c7-5cc5-4438-b540-9645e46f878b 
    ## INFO  [15:59:56.784] [bbotk]  648908c4-2950-41de-93b5-7e385629addd 
    ## INFO  [15:59:56.784] [bbotk]  0b4d20ca-bba7-4940-b55f-521eba6c66b1 
    ## INFO  [15:59:56.784] [bbotk]  d1db8cae-12fa-44aa-95fa-87f4452b59a4 
    ## INFO  [15:59:56.784] [bbotk]  d3bd4d4e-1a95-4515-856d-bce21827f358 
    ## INFO  [15:59:56.784] [bbotk]  b3c5bcef-f835-40e4-818f-269b381d8a6e 
    ## INFO  [15:59:56.784] [bbotk]  fbc98f61-728a-4c56-8d8e-e2206e8d682e 
    ## INFO  [15:59:56.784] [bbotk]  fd32a2f6-53e6-4524-a6ec-5bc5e04bfa53 
    ## INFO  [15:59:56.787] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [15:59:57.497] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [15:59:57.505] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:57.581] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:57.752] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:57.880] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.015] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.147] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.235] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.349] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.503] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.604] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [15:59:58.709] [mlr3]  Finished benchmark 
    ## INFO  [15:59:59.103] [bbotk] Result of batch 2: 
    ## INFO  [15:59:59.129] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [15:59:59.129] [bbotk]       FALSE      FALSE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]       FALSE      FALSE       TRUE      FALSE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]        TRUE       TRUE      FALSE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]        TRUE      FALSE      FALSE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]       FALSE       TRUE       TRUE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [15:59:59.129] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE    FALSE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE     TRUE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [15:59:59.129] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE FALSE  TRUE         FALSE   TRUE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:59.129] [bbotk]     TRUE  FALSE  TRUE        FALSE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE  TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]     TRUE  FALSE  TRUE        FALSE   TRUE  FALSE FALSE  TRUE FALSE FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]     TRUE  FALSE  TRUE        FALSE  FALSE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [15:59:59.129] [bbotk]   TRUE FALSE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE    TRUE FALSE FALSE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [15:59:59.129] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [15:59:59.129] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [15:59:59.129] [bbotk]     FALSE FALSE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]     FALSE FALSE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]      TRUE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [15:59:59.129] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [15:59:59.129] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE     TRUE FALSE  TRUE  TRUE    TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE     TRUE FALSE  TRUE FALSE    TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE    FALSE  TRUE  TRUE  TRUE   FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [15:59:59.129] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE     FALSE     FALSE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE     FALSE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE     FALSE     FALSE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [15:59:59.129] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE  FALSE  FALSE  FALSE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE   TRUE FALSE FALSE  TRUE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE  FALSE FALSE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE FALSE FALSE  FALSE  TRUE  TRUE FALSE FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE FALSE FALSE   TRUE FALSE  TRUE FALSE FALSE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE  TRUE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [15:59:59.129] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [15:59:59.129] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [15:59:59.129] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE         TRUE         TRUE        FALSE       FALSE  TRUE  TRUE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [15:59:59.129] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [15:59:59.129] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [15:59:59.129] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [15:59:59.129] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:59.129] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [15:59:59.129] [bbotk]      FALSE FALSE  TRUE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [15:59:59.129] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:59.129] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE    TRUE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [15:59:59.129] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [15:59:59.129] [bbotk]      FALSE  TRUE  TRUE FALSE   FALSE   FALSE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [15:59:59.129] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [15:59:59.129] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [15:59:59.129] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [15:59:59.129] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   FALSE   FALSE    TRUE  TRUE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE    TRUE    TRUE    TRUE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   FALSE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]    TRUE   FALSE   FALSE   FALSE  TRUE  FALSE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [15:59:59.129] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [15:59:59.129] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.5803571 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE  TRUE    0.5997024 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.5907738 
    ## INFO  [15:59:59.129] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE   FALSE FALSE FALSE  TRUE FALSE    0.7187500 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE FALSE    0.6011905 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE    TRUE FALSE  TRUE  TRUE FALSE    0.7187500 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE    0.6711310 
    ## INFO  [15:59:59.129] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.6011905 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE    0.6056548 
    ## INFO  [15:59:59.129] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE FALSE    0.6086310 
    ## INFO  [15:59:59.129] [bbotk]                                 uhash 
    ## INFO  [15:59:59.129] [bbotk]  d7b5f9b9-efca-4ca8-ba17-8cdf4da518f1 
    ## INFO  [15:59:59.129] [bbotk]  e329c51e-4e77-40d2-ac39-7b79841ad589 
    ## INFO  [15:59:59.129] [bbotk]  cb6da8eb-67e6-43c7-bee0-40537fe7a8dd 
    ## INFO  [15:59:59.129] [bbotk]  3ca5ba3c-9a58-4f67-a9d2-e4a2fd4a71d7 
    ## INFO  [15:59:59.129] [bbotk]  62d1a138-5e83-4684-9ae5-9b1819743665 
    ## INFO  [15:59:59.129] [bbotk]  00f8b284-7f06-4d8c-969f-8e35dd33f1ff 
    ## INFO  [15:59:59.129] [bbotk]  fdd9d432-2790-4fd4-9416-8881741090b0 
    ## INFO  [15:59:59.129] [bbotk]  fe9b1257-cb08-4408-a6fd-67bf2b378e0b 
    ## INFO  [15:59:59.129] [bbotk]  058e6c1f-b527-4267-a69e-107480bdd979 
    ## INFO  [15:59:59.129] [bbotk]  eb3c6a90-9d29-4160-a174-6201aab2a3b4 
    ## INFO  [15:59:59.144] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [15:59:59.145] [bbotk] Result: 
    ## INFO  [15:59:59.161] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [15:59:59.161] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE TRUE TRUE 
    ## INFO  [15:59:59.161] [bbotk]   AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [15:59:59.161] [bbotk]  TRUE        TRUE TRUE TRUE     TRUE     TRUE  TRUE TRUE     TRUE      TRUE 
    ## INFO  [15:59:59.161] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [15:59:59.161] [bbotk]  TRUE   TRUE FALSE    TRUE TRUE  TRUE  TRUE          TRUE   TRUE  TRUE TRUE 
    ## INFO  [15:59:59.161] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [15:59:59.161] [bbotk]     TRUE  FALSE  TRUE         TRUE  FALSE   TRUE  TRUE TRUE FALSE  TRUE   TRUE 
    ## INFO  [15:59:59.161] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [15:59:59.161] [bbotk]  FALSE TRUE    TRUE FALSE  TRUE FALSE  FALSE     TRUE      TRUE  FALSE     TRUE 
    ## INFO  [15:59:59.161] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [15:59:59.161] [bbotk]   TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE FALSE 
    ## INFO  [15:59:59.161] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [15:59:59.161] [bbotk]      TRUE FALSE FALSE  TRUE   FALSE FALSE TRUE  TRUE TRUE TRUE  TRUE FALSE 
    ## INFO  [15:59:59.161] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [15:59:59.161] [bbotk]       TRUE     FALSE     FALSE      TRUE      TRUE      TRUE   TRUE   TRUE 
    ## INFO  [15:59:59.161] [bbotk]  MAGEC3 MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP  PI3 
    ## INFO  [15:59:59.161] [bbotk]    TRUE TRUE TRUE TRUE   TRUE   TRUE   TRUE  TRUE  FALSE TRUE FALSE TRUE TRUE 
    ## INFO  [15:59:59.161] [bbotk]  PITX1   PLG PRR15L PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [15:59:59.161] [bbotk]  FALSE FALSE   TRUE TRUE  TRUE  TRUE TRUE         TRUE         FALSE 
    ## INFO  [15:59:59.161] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [15:59:59.161] [bbotk]           TRUE        FALSE         FALSE          TRUE          TRUE 
    ## INFO  [15:59:59.161] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [15:59:59.161] [bbotk]          TRUE       FALSE          TRUE        FALSE        FALSE         FALSE 
    ## INFO  [15:59:59.161] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [15:59:59.161] [bbotk]          TRUE         TRUE         TRUE        TRUE TRUE TRUE      TRUE TRUE 
    ## INFO  [15:59:59.161] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [15:59:59.161] [bbotk]   TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   FALSE   TRUE   FALSE 
    ## INFO  [15:59:59.161] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [15:59:59.161] [bbotk]     TRUE    TRUE  TRUE   TRUE TRUE    TRUE   TRUE     FALSE  TRUE TRUE   TRUE 
    ## INFO  [15:59:59.161] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [15:59:59.161] [bbotk]   FALSE  TRUE    TRUE TRUE  TRUE TRUE TRUE 
    ## INFO  [15:59:59.161] [bbotk]                                                        features classif.bacc 
    ## INFO  [15:59:59.161] [bbotk]  AC003092.1,AC006262.4,AC006262.5,AC007879.6,AC104654.2,AFM,...    0.7276786 
    ## INFO  [15:59:59.291] [mlr3]  Applying learner 'classif.ksvm.fselector' on task 'kirc_cla' (iter 2/3) 
    ## INFO  [15:59:59.517] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [15:59:59.520] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:00.290] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:00.298] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:00.428] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:00.600] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:00.771] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:00.896] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:01.036] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:01.172] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:01.278] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:01.378] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:01.542] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:01.663] [mlr3]  Finished benchmark 
    ## INFO  [16:00:02.099] [bbotk] Result of batch 1: 
    ## INFO  [16:00:02.117] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:02.117] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]       FALSE      FALSE       TRUE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]       FALSE      FALSE       TRUE      FALSE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]        TRUE       TRUE      FALSE      FALSE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]        TRUE      FALSE      FALSE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:02.117] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE        TRUE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE    FALSE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:02.117] [bbotk]  FALSE   TRUE FALSE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE FALSE  TRUE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE FALSE FALSE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE   TRUE FALSE   FALSE  TRUE FALSE FALSE         FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE FALSE FALSE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:02.117] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]     TRUE   TRUE FALSE         TRUE  FALSE  FALSE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE   TRUE FALSE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE   TRUE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE  FALSE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]    FALSE   TRUE  TRUE        FALSE   TRUE   TRUE FALSE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE FALSE    TRUE FALSE FALSE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE  TRUE    TRUE  TRUE FALSE  TRUE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:02.117] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:02.117] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:02.117] [bbotk]      TRUE FALSE FALSE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:02.117] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:02.117] [bbotk]      TRUE FALSE  TRUE    FALSE    TRUE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:02.117] [bbotk]     FALSE  TRUE FALSE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:02.117] [bbotk]     FALSE FALSE FALSE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:02.117] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:02.117] [bbotk]      TRUE  TRUE FALSE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:02.117] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:02.117] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:02.117] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE     TRUE  TRUE FALSE FALSE    TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE    FALSE  TRUE  TRUE FALSE   FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE     TRUE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:02.117] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE      TRUE      TRUE      TRUE     FALSE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE     FALSE      TRUE      TRUE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE      TRUE     FALSE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:02.117] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:02.117] [bbotk]    TRUE   TRUE  TRUE  TRUE FALSE  FALSE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE   TRUE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE  FALSE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:02.117] [bbotk]  FALSE  TRUE FALSE FALSE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE FALSE FALSE   TRUE FALSE  TRUE FALSE FALSE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE FALSE  TRUE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]  FALSE  TRUE  TRUE  TRUE  FALSE FALSE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE FALSE FALSE   TRUE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE  TRUE FALSE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:02.117] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:02.117] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         TRUE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE        FALSE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:02.117] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE         TRUE         TRUE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE         TRUE         TRUE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]          FALSE        FALSE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE        FALSE         TRUE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:02.117] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:02.117] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:02.117] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:02.117] [bbotk]       TRUE FALSE  TRUE FALSE   FALSE    TRUE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:02.117] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:02.117] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:02.117] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:02.117] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:02.117] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:02.117] [bbotk]       TRUE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:02.117] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:02.117] [bbotk]      FALSE  TRUE FALSE  TRUE    TRUE    TRUE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:02.117] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:02.117] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE   TRUE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE    TRUE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]    TRUE    TRUE   FALSE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE    TRUE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:02.117] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7391304 
    ## INFO  [16:00:02.117] [bbotk]  FALSE FALSE  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE FALSE    0.6766304 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE FALSE    0.6345109 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7187500 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.8016304 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE   FALSE  TRUE  TRUE  TRUE  TRUE    0.7391304 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE    0.7391304 
    ## INFO  [16:00:02.117] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE    0.6970109 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE   FALSE  TRUE  TRUE  TRUE FALSE    0.7391304 
    ## INFO  [16:00:02.117] [bbotk]   TRUE  TRUE   TRUE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.8016304 
    ## INFO  [16:00:02.117] [bbotk]                                 uhash 
    ## INFO  [16:00:02.117] [bbotk]  9adb7eef-9919-4af7-a644-f7cc1b451c87 
    ## INFO  [16:00:02.117] [bbotk]  04616788-2421-4429-a453-556e729fc9f8 
    ## INFO  [16:00:02.117] [bbotk]  be7f3ab9-caca-4152-a1a3-3e3193e13895 
    ## INFO  [16:00:02.117] [bbotk]  0ed83aa0-e680-4f20-a8bb-806481a7266e 
    ## INFO  [16:00:02.117] [bbotk]  a884992a-ce8a-49f7-aae7-5bf1de7f1ad8 
    ## INFO  [16:00:02.117] [bbotk]  beb0049b-b83c-41c8-a29b-8edeb4891747 
    ## INFO  [16:00:02.117] [bbotk]  7ac4829f-c420-4aca-ad83-18d74f0c2953 
    ## INFO  [16:00:02.117] [bbotk]  8c0a042c-73dc-48c7-afb6-f27b74755ed7 
    ## INFO  [16:00:02.117] [bbotk]  8b4ba709-d074-4d13-b62b-bf5f970004bd 
    ## INFO  [16:00:02.117] [bbotk]  0817e699-bb21-4992-a7bd-4e112811fca1 
    ## INFO  [16:00:02.120] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:02.874] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:02.882] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:02.988] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.174] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.279] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.372] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.459] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.606] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.710] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:03.833] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:04.008] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:04.119] [mlr3]  Finished benchmark 
    ## INFO  [16:00:04.535] [bbotk] Result of batch 2: 
    ## INFO  [16:00:04.552] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:04.552] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]        TRUE      FALSE       TRUE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]       FALSE       TRUE       TRUE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]       FALSE       TRUE      FALSE      FALSE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:04.552] [bbotk]   TRUE       FALSE  TRUE FALSE    FALSE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE       FALSE  TRUE  TRUE    FALSE    FALSE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE     TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE   TRUE FALSE    TRUE  TRUE FALSE FALSE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE FALSE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:04.552] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE  TRUE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE   TRUE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE FALSE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]     TRUE  FALSE  TRUE        FALSE   TRUE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE   FALSE FALSE  TRUE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE   FALSE  TRUE FALSE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE   FALSE FALSE FALSE  TRUE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE FALSE   FALSE FALSE  TRUE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE FALSE   FALSE  TRUE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE FALSE    TRUE FALSE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:04.552] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:04.552] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:04.552] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:04.552] [bbotk]      TRUE FALSE  TRUE    FALSE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:04.552] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:04.552] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:04.552] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:04.552] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:04.552] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:04.552] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:04.552] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:04.552] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE    FALSE  TRUE  TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE     TRUE  TRUE  TRUE FALSE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:04.552] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE     FALSE     FALSE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE      TRUE     FALSE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:04.552] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:04.552] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE   TRUE FALSE FALSE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE FALSE FALSE FALSE  FALSE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE FALSE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE  TRUE FALSE  FALSE  TRUE  TRUE FALSE FALSE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:04.552] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:04.552] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:04.552] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         TRUE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:04.552] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:04.552] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         TRUE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE        FALSE        FALSE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:04.552] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:04.552] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:04.552] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:04.552] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:04.552] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:04.552] [bbotk]       TRUE  TRUE FALSE FALSE   FALSE    TRUE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:04.552] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:04.552] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:04.552] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:04.552] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:04.552] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:04.552] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:04.552] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:04.552] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:04.552] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE    TRUE   FALSE    TRUE FALSE   TRUE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE   FALSE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE  FALSE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:04.552] [bbotk]    TRUE   FALSE    TRUE   FALSE  TRUE  FALSE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:04.552] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE   TRUE   TRUE FALSE    TRUE  TRUE FALSE  TRUE FALSE    0.7391304 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6345109 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE  TRUE    0.8016304 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE  TRUE FALSE FALSE    0.7391304 
    ## INFO  [16:00:04.552] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.7391304 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6875000 
    ## INFO  [16:00:04.552] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7391304 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6766304 
    ## INFO  [16:00:04.552] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7798913 
    ## INFO  [16:00:04.552] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE    0.7703804 
    ## INFO  [16:00:04.552] [bbotk]                                 uhash 
    ## INFO  [16:00:04.552] [bbotk]  d44ff558-b1e2-463c-9d02-02e9393fa9c0 
    ## INFO  [16:00:04.552] [bbotk]  8b1abe75-4d3f-44a4-a498-9b0099670b0a 
    ## INFO  [16:00:04.552] [bbotk]  a8a50011-b157-4be0-a1e7-39387745d2dd 
    ## INFO  [16:00:04.552] [bbotk]  3b496230-4c0d-4174-815b-729c12bcb764 
    ## INFO  [16:00:04.552] [bbotk]  2064accd-0e93-4d49-8fbd-895b14ce0889 
    ## INFO  [16:00:04.552] [bbotk]  2de957d3-d729-4df0-b130-03e2cd8cca5d 
    ## INFO  [16:00:04.552] [bbotk]  d8ccf75e-6f08-4200-83bc-4519cda600f7 
    ## INFO  [16:00:04.552] [bbotk]  dfa2348f-e063-42ac-93e5-a3fb59928f5e 
    ## INFO  [16:00:04.552] [bbotk]  5b794528-b3c7-4fda-bef9-b66f7a97e03c 
    ## INFO  [16:00:04.552] [bbotk]  7a7d9704-7698-4c51-a1e9-3d3729c20337 
    ## INFO  [16:00:04.564] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:04.565] [bbotk] Result: 
    ## INFO  [16:00:04.576] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:04.576] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:04.576] [bbotk]   AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:04.576] [bbotk]  TRUE        TRUE TRUE TRUE    FALSE     TRUE  TRUE TRUE     TRUE     FALSE 
    ## INFO  [16:00:04.576] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [16:00:04.576] [bbotk]  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE TRUE 
    ## INFO  [16:00:04.576] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:04.576] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:04.576] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:04.576] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE TRUE   TRUE    FALSE     FALSE   TRUE     TRUE 
    ## INFO  [16:00:04.576] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:04.576] [bbotk]  FALSE  TRUE    FALSE    TRUE  FALSE    TRUE FALSE    FALSE    FALSE  TRUE 
    ## INFO  [16:00:04.576] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15  KNG1  KRT7 L1CAM LECT1 
    ## INFO  [16:00:04.576] [bbotk]     FALSE  TRUE  TRUE FALSE   FALSE  TRUE TRUE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:04.576] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:04.576] [bbotk]       TRUE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE   TRUE 
    ## INFO  [16:00:04.576] [bbotk]  MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP 
    ## INFO  [16:00:04.576] [bbotk]   FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE TRUE  TRUE TRUE 
    ## INFO  [16:00:04.576] [bbotk]    PI3 PITX1  PLG PRR15L  PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:04.576] [bbotk]  FALSE FALSE TRUE  FALSE FALSE FALSE  TRUE TRUE        FALSE         FALSE 
    ## INFO  [16:00:04.576] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:04.576] [bbotk]          FALSE        FALSE         FALSE         FALSE          TRUE 
    ## INFO  [16:00:04.576] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:04.576] [bbotk]          TRUE       FALSE         FALSE         TRUE        FALSE         FALSE 
    ## INFO  [16:00:04.576] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1 SAA1 SAA2_SAA4  SBSN 
    ## INFO  [16:00:04.576] [bbotk]         FALSE        FALSE        FALSE        TRUE FALSE TRUE     FALSE FALSE 
    ## INFO  [16:00:04.576] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:04.576] [bbotk]  FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   TRUE   FALSE 
    ## INFO  [16:00:04.576] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1   TTR TUBA3E 
    ## INFO  [16:00:04.576] [bbotk]    FALSE    TRUE FALSE   TRUE FALSE   FALSE   TRUE     FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:04.576] [bbotk]  TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2  ZIC5 
    ## INFO  [16:00:04.576] [bbotk]   FALSE FALSE   FALSE FALSE FALSE TRUE FALSE 
    ## INFO  [16:00:04.576] [bbotk]                                                   features classif.bacc 
    ## INFO  [16:00:04.576] [bbotk]  AC003092.1,AC006262.5,AC007879.6,AMH,APCDD1L_AS1,AQP2,...    0.8016304 
    ## INFO  [16:00:04.637] [mlr3]  Applying learner 'classif.svm.fselector' on task 'kirc_cla' (iter 2/3) 
    ## INFO  [16:00:04.898] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:04.900] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:05.655] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:05.668] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:05.747] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:05.823] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:05.895] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:05.993] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:06.079] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:06.157] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:06.278] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:06.388] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:06.469] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:06.537] [mlr3]  Finished benchmark 
    ## INFO  [16:00:06.970] [bbotk] Result of batch 1: 
    ## INFO  [16:00:06.981] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:06.981] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]       FALSE       TRUE       TRUE      FALSE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]        TRUE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]        TRUE       TRUE       TRUE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:06.981] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE       FALSE FALSE  TRUE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:06.981] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  FALSE  TRUE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE  TRUE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:06.981] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE   TRUE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE  FALSE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]    FALSE   TRUE FALSE         TRUE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE  FALSE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE    TRUE FALSE  TRUE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  TRUE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:06.981] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:06.981] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE FALSE  TRUE     TRUE   FALSE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:06.981] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE FALSE FALSE    FALSE    TRUE  FALSE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:06.981] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:06.981] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:06.981] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE FALSE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE     TRUE FALSE  TRUE FALSE    TRUE  TRUE FALSE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE    TRUE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:06.981] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:06.981] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:06.981] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE   TRUE  TRUE  TRUE FALSE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE  FALSE FALSE  TRUE  TRUE  FALSE   TRUE  FALSE FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE   TRUE  TRUE FALSE FALSE  FALSE  FALSE  FALSE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE   TRUE FALSE FALSE FALSE   TRUE   TRUE   TRUE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE  TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE FALSE  TRUE  FALSE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:06.981] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:06.981] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         TRUE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:06.981] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE        FALSE         TRUE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE         TRUE        FALSE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]           TRUE        FALSE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:06.981] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:06.981] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:06.981] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:06.981] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:06.981] [bbotk]       TRUE  TRUE FALSE  TRUE    TRUE   FALSE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:06.981] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:06.981] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:06.981] [bbotk]      FALSE FALSE FALSE FALSE   FALSE    TRUE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:06.981] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:06.981] [bbotk]      FALSE  TRUE FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:06.981] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE   FALSE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:06.981] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:06.981] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:06.981] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE    TRUE    TRUE   FALSE FALSE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]    TRUE   FALSE   FALSE   FALSE FALSE  FALSE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE    TRUE    TRUE   FALSE FALSE  FALSE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:06.981] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7001570 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.6899529 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6899529 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6899529 
    ## INFO  [16:00:06.981] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE FALSE    0.6412873 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE   FALSE  TRUE  TRUE FALSE  TRUE    0.6899529 
    ## INFO  [16:00:06.981] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.7590267 
    ## INFO  [16:00:06.981] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6232339 
    ## INFO  [16:00:06.981] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.6899529 
    ## INFO  [16:00:06.981] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.6412873 
    ## INFO  [16:00:06.981] [bbotk]                                 uhash 
    ## INFO  [16:00:06.981] [bbotk]  3dda08d7-3c25-4b3c-9ea2-4589e44b8900 
    ## INFO  [16:00:06.981] [bbotk]  5e8cad40-8647-4356-b312-8e66f5679ae9 
    ## INFO  [16:00:06.981] [bbotk]  f828643a-87f5-42d0-a699-e699b3b8e8c6 
    ## INFO  [16:00:06.981] [bbotk]  c764356a-d480-4b46-941d-429430105f56 
    ## INFO  [16:00:06.981] [bbotk]  bcac21db-9c50-41f6-9bc5-d16f34eeb104 
    ## INFO  [16:00:06.981] [bbotk]  961870e8-60ff-4102-95e4-166c889d65c0 
    ## INFO  [16:00:06.981] [bbotk]  7c454037-733f-49ff-8b40-3fa2a4643eb0 
    ## INFO  [16:00:06.981] [bbotk]  0977bd63-2394-4ddb-8f72-a14213c5ae31 
    ## INFO  [16:00:06.981] [bbotk]  05e7fcd0-cf6e-43ba-832a-861ea757aea9 
    ## INFO  [16:00:06.981] [bbotk]  c8d86a71-771b-4590-9be5-8bbddb338d7b 
    ## INFO  [16:00:06.984] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:07.778] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:07.797] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:07.895] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:07.983] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.048] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.137] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.278] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.371] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.446] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.553] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.640] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:08.718] [mlr3]  Finished benchmark 
    ## INFO  [16:00:09.130] [bbotk] Result of batch 2: 
    ## INFO  [16:00:09.144] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:09.144] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE       TRUE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]       FALSE       TRUE       TRUE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:09.144] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE       FALSE FALSE  TRUE     TRUE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE    FALSE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE       FALSE  TRUE  TRUE     TRUE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE        TRUE  TRUE  TRUE    FALSE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:09.144] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE   TRUE FALSE    TRUE  TRUE FALSE  TRUE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE   TRUE FALSE   FALSE  TRUE FALSE  TRUE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE   TRUE FALSE    TRUE FALSE FALSE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE   TRUE  TRUE   FALSE FALSE  TRUE  TRUE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:09.144] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE   TRUE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE  FALSE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE   FALSE FALSE  TRUE FALSE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE   FALSE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]     FALSE  TRUE FALSE     TRUE    TRUE  FALSE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:09.144] [bbotk]      TRUE  TRUE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:09.144] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:09.144] [bbotk]  FALSE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE     TRUE FALSE FALSE FALSE    TRUE FALSE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE    FALSE  TRUE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:09.144] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE     FALSE      TRUE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE      TRUE      TRUE     FALSE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:09.144] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:09.144] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE  FALSE  TRUE FALSE  TRUE  FALSE  FALSE  FALSE FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE   TRUE  FALSE   TRUE FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE   TRUE  FALSE  FALSE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE  FALSE  FALSE  FALSE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE FALSE  TRUE  FALSE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:09.144] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         FALSE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:09.144] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE        FALSE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:09.144] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE         TRUE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:09.144] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:09.144] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:09.144] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:09.144] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:09.144] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:09.144] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE   FALSE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:09.144] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:09.144] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE    TRUE    TRUE   FALSE FALSE   TRUE FALSE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE    TRUE   FALSE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE  FALSE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:09.144] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]   FALSE   FALSE    TRUE   FALSE  TRUE  FALSE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:09.144] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE FALSE  TRUE    0.7001570 
    ## INFO  [16:00:09.144] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE    TRUE  TRUE  TRUE FALSE FALSE    0.6514914 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6797488 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE    0.7001570 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE   FALSE FALSE  TRUE  TRUE FALSE    0.7488226 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE    0.6899529 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7001570 
    ## INFO  [16:00:09.144] [bbotk]   TRUE FALSE  FALSE  FALSE  TRUE   FALSE  TRUE FALSE FALSE FALSE    0.6334380 
    ## INFO  [16:00:09.144] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6899529 
    ## INFO  [16:00:09.144] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.6899529 
    ## INFO  [16:00:09.144] [bbotk]                                 uhash 
    ## INFO  [16:00:09.144] [bbotk]  b67d59a6-229c-497c-8bc0-ef0540b4305f 
    ## INFO  [16:00:09.144] [bbotk]  25db8820-d4e6-4dd3-94ac-08327ccb3d43 
    ## INFO  [16:00:09.144] [bbotk]  2c493bea-5af3-49bb-bdb8-6b940fcab8d9 
    ## INFO  [16:00:09.144] [bbotk]  d774ee22-5546-4bba-ad8f-bfcf93215316 
    ## INFO  [16:00:09.144] [bbotk]  a0596870-d6f2-4609-8902-a025c5196ccf 
    ## INFO  [16:00:09.144] [bbotk]  90b698e7-4f1b-4779-9c01-70586d05a378 
    ## INFO  [16:00:09.144] [bbotk]  76850202-c48b-4760-9088-c2571b42702f 
    ## INFO  [16:00:09.144] [bbotk]  f4067aa4-81af-4da2-8f45-1287b3b6d8b0 
    ## INFO  [16:00:09.144] [bbotk]  3e89ddd9-49a0-4eb4-b26a-85e8d5ddde89 
    ## INFO  [16:00:09.144] [bbotk]  36fefcd2-5726-4857-93b2-de4d42fb8fcb 
    ## INFO  [16:00:09.151] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:09.152] [bbotk] Result: 
    ## INFO  [16:00:09.158] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM AHSG 
    ## INFO  [16:00:09.158] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE TRUE 
    ## INFO  [16:00:09.158] [bbotk]    AMH APCDD1L_AS1 AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:09.158] [bbotk]  FALSE        TRUE TRUE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:09.158] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:09.158] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:09.158] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:09.158] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:09.158] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:09.158] [bbotk]  FALSE FALSE    TRUE FALSE  TRUE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:09.158] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:09.158] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:09.158] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15 KNG1  KRT7 L1CAM 
    ## INFO  [16:00:09.158] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE TRUE FALSE TRUE FALSE FALSE 
    ## INFO  [16:00:09.158] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:09.158] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:09.158] [bbotk]  MAGEC2 MAGEC3 MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:09.158] [bbotk]   FALSE  FALSE TRUE FALSE FALSE  FALSE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:09.158] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25 RHBG RP11_10O22.1 
    ## INFO  [16:00:09.158] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE TRUE         TRUE 
    ## INFO  [16:00:09.158] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:09.158] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:09.158] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:09.158] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:09.158] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:09.158] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:09.158] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:09.158] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:09.158] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:09.158] [bbotk]    TRUE   FALSE   FALSE   FALSE FALSE  FALSE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:09.158] [bbotk]  TNNT1  TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2  ZIC5 
    ## INFO  [16:00:09.158] [bbotk]  FALSE TRUE  FALSE  FALSE FALSE   FALSE FALSE FALSE TRUE FALSE 
    ## INFO  [16:00:09.158] [bbotk]                                             features classif.bacc 
    ## INFO  [16:00:09.158] [bbotk]  AC006262.5,AHSG,APCDD1L_AS1,AQP2,CDC42P2,CXCL13,...    0.7590267 
    ## INFO  [16:00:09.194] [mlr3]  Applying learner 'classif.ranger.fselector' on task 'kirc_cla' (iter 2/3) 
    ## INFO  [16:00:09.414] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:09.416] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:10.381] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:10.389] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:10.535] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:10.685] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:10.814] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:10.935] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:11.091] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:11.242] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:11.372] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:11.543] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:11.700] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:11.845] [mlr3]  Finished benchmark 
    ## INFO  [16:00:12.347] [bbotk] Result of batch 1: 
    ## INFO  [16:00:12.358] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:12.358] [bbotk]        TRUE       TRUE      FALSE      FALSE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:12.358] [bbotk]   TRUE        TRUE FALSE FALSE     TRUE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE       FALSE FALSE  TRUE    FALSE    FALSE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE    FALSE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE        TRUE FALSE  TRUE    FALSE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE     TRUE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE       FALSE  TRUE  TRUE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE FALSE  TRUE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE FALSE  TRUE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE FALSE  TRUE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:12.358] [bbotk]     TRUE  FALSE FALSE         TRUE  FALSE   TRUE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]    FALSE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE  FALSE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]     TRUE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:12.358] [bbotk]   TRUE  TRUE    TRUE FALSE FALSE  TRUE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:12.358] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE FALSE FALSE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE FALSE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]      TRUE FALSE FALSE    FALSE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:12.358] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:12.358] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:12.358] [bbotk]   TRUE    FALSE  TRUE  TRUE FALSE   FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE    FALSE FALSE  TRUE  TRUE   FALSE FALSE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE    TRUE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE     TRUE  TRUE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE    FALSE FALSE  TRUE  TRUE    TRUE FALSE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:12.358] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE     FALSE      TRUE     FALSE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:12.358] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:12.358] [bbotk]   FALSE   TRUE FALSE FALSE FALSE  FALSE  FALSE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE  FALSE  TRUE FALSE  TRUE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE   TRUE FALSE FALSE  TRUE   TRUE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE  FALSE FALSE  TRUE  TRUE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE FALSE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  TRUE FALSE  TRUE  FALSE FALSE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  TRUE FALSE   TRUE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE  TRUE FALSE FALSE         TRUE 
    ## INFO  [16:00:12.358] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE  TRUE FALSE FALSE FALSE         TRUE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  TRUE FALSE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  TRUE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:12.358] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:12.358] [bbotk]           TRUE        FALSE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         TRUE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]           TRUE         TRUE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         TRUE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:12.358] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         TRUE        FALSE         TRUE       FALSE FALSE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]           TRUE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]           TRUE        FALSE         TRUE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:12.358] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:12.358] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE    TRUE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE    TRUE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:12.358] [bbotk]       TRUE  TRUE  TRUE FALSE   FALSE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:12.358] [bbotk]       TRUE FALSE FALSE  TRUE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE  TRUE  TRUE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:12.358] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:12.358] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:12.358] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE   TRUE FALSE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE   FALSE    TRUE   FALSE  TRUE  FALSE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE    TRUE   FALSE   FALSE  TRUE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE    TRUE    TRUE    TRUE FALSE  FALSE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]    TRUE    TRUE   FALSE   FALSE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE    TRUE   FALSE   FALSE  TRUE  FALSE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:12.358] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:12.358] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE    0.6553030 
    ## INFO  [16:00:12.358] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE  TRUE    0.6439394 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5883838 
    ## INFO  [16:00:12.358] [bbotk]   TRUE  TRUE   TRUE  FALSE FALSE   FALSE  TRUE  TRUE  TRUE  TRUE    0.6439394 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5997475 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.6439394 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE  TRUE  TRUE    0.6717172 
    ## INFO  [16:00:12.358] [bbotk]   TRUE FALSE   TRUE   TRUE FALSE    TRUE  TRUE FALSE FALSE  TRUE    0.6717172 
    ## INFO  [16:00:12.358] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE    TRUE  TRUE FALSE FALSE  TRUE    0.6439394 
    ## INFO  [16:00:12.358] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6439394 
    ## INFO  [16:00:12.358] [bbotk]                                 uhash 
    ## INFO  [16:00:12.358] [bbotk]  34f3bf41-4f22-4ceb-9018-d639cc9faaa3 
    ## INFO  [16:00:12.358] [bbotk]  5c030ef3-d154-46dd-b209-0ff2b17722e6 
    ## INFO  [16:00:12.358] [bbotk]  98e62896-0b40-4667-b800-f9a050446cc8 
    ## INFO  [16:00:12.358] [bbotk]  ea4ef1cf-7942-4dc0-bfcc-608c4a81243e 
    ## INFO  [16:00:12.358] [bbotk]  7db96542-1520-46ca-82b3-718fab8497e4 
    ## INFO  [16:00:12.358] [bbotk]  183d3518-9ab2-46b2-9e21-943f8e60d028 
    ## INFO  [16:00:12.358] [bbotk]  7705996b-bf15-416d-bc3d-6a61493d997a 
    ## INFO  [16:00:12.358] [bbotk]  c9157d25-a00f-4a34-8636-a1cf2f49dfd5 
    ## INFO  [16:00:12.358] [bbotk]  6b378b89-045d-4d13-90aa-e6e2e9d25559 
    ## INFO  [16:00:12.358] [bbotk]  0eb67aeb-cf4b-44b5-a4a2-bc6e1efe0d15 
    ## INFO  [16:00:12.361] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:13.411] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:13.418] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:13.558] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:13.739] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:13.902] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.028] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.157] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.314] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.470] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.609] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.783] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:14.937] [mlr3]  Finished benchmark 
    ## INFO  [16:00:15.417] [bbotk] Result of batch 2: 
    ## INFO  [16:00:15.428] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:15.428] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:15.428] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]       FALSE       TRUE      FALSE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:15.428] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:15.428] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE       FALSE FALSE  TRUE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE  TRUE FALSE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE   TRUE FALSE   FALSE FALSE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:15.428] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]     TRUE   TRUE FALSE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE  FALSE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  TRUE   FALSE  TRUE FALSE FALSE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE  TRUE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE  TRUE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:15.428] [bbotk]     FALSE FALSE FALSE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]     FALSE FALSE  TRUE     TRUE    TRUE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:15.428] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:15.428] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:15.428] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]     FALSE FALSE  TRUE    FALSE    TRUE  FALSE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:15.428] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     TRUE FALSE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE     TRUE  TRUE  TRUE FALSE   FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE    TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE   FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE   FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:15.428] [bbotk]   TRUE      TRUE      TRUE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     FALSE     FALSE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     FALSE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE      TRUE     FALSE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.428] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:15.428] [bbotk]   FALSE   TRUE FALSE FALSE FALSE  FALSE  FALSE   TRUE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE   TRUE FALSE FALSE FALSE  FALSE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE  FALSE FALSE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE   TRUE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE  FALSE  FALSE  FALSE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:15.428] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE FALSE FALSE  FALSE  TRUE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE FALSE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE FALSE  TRUE FALSE FALSE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE FALSE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:15.428] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:15.428] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:15.428] [bbotk]           TRUE        FALSE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:15.428] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:15.428] [bbotk]          FALSE         TRUE        FALSE         TRUE       FALSE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE        FALSE         TRUE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:15.428] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:15.428] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:15.428] [bbotk]      FALSE  TRUE  TRUE  TRUE   FALSE    TRUE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:15.428] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:15.428] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:15.428] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:15.428] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:15.428] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:15.428] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:15.428] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:15.428] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:15.428] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:15.428] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:15.428] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE    TRUE   FALSE    TRUE FALSE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:15.428] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:15.428] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:15.428] [bbotk]   TRUE FALSE  FALSE  FALSE FALSE    TRUE  TRUE FALSE  TRUE  TRUE    0.6325758 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE   TRUE  FALSE  TRUE   FALSE  TRUE FALSE FALSE FALSE    0.6439394 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6553030 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6439394 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE    0.6275253 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6717172 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.6553030 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.4709596 
    ## INFO  [16:00:15.428] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE    0.6717172 
    ## INFO  [16:00:15.428] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE  TRUE FALSE    0.6553030 
    ## INFO  [16:00:15.428] [bbotk]                                 uhash 
    ## INFO  [16:00:15.428] [bbotk]  40dca05e-d81b-4199-bab5-e4b773d14328 
    ## INFO  [16:00:15.428] [bbotk]  5733ba35-ab70-4449-bce4-6cda6e17d804 
    ## INFO  [16:00:15.428] [bbotk]  c3b3747f-9a39-4688-8910-895c7e3e8cc3 
    ## INFO  [16:00:15.428] [bbotk]  d32581a3-ef68-48e2-a7a4-88d19eca9590 
    ## INFO  [16:00:15.428] [bbotk]  fa23b24f-2f61-4ee1-8ec9-655dea5ee82e 
    ## INFO  [16:00:15.428] [bbotk]  9cce594d-e9f3-4be6-8631-66f6e9c84766 
    ## INFO  [16:00:15.428] [bbotk]  a9f68169-0fd1-4d5c-9d2e-a7d8b397e94b 
    ## INFO  [16:00:15.428] [bbotk]  e69915d6-fa0b-4822-b8a9-02c8c82ea6e6 
    ## INFO  [16:00:15.428] [bbotk]  0a07e868-4e14-4348-87f0-011d06c6b6fe 
    ## INFO  [16:00:15.428] [bbotk]  194fc6df-74b0-4e34-8506-afee93b48657 
    ## INFO  [16:00:15.435] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:15.436] [bbotk] Result: 
    ## INFO  [16:00:15.442] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:15.442] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:15.442] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:15.442] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:15.442] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:15.442] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:15.442] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:15.442] [bbotk]     TRUE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE TRUE FALSE  FALSE 
    ## INFO  [16:00:15.442] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:15.442] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:15.442] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:15.442] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:15.442] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:15.442] [bbotk]  FALSE    FALSE FALSE  TRUE  TRUE    TRUE FALSE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:15.442] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:15.442] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:15.442] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:15.442] [bbotk]   FALSE   TRUE FALSE FALSE TRUE   TRUE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:15.442] [bbotk]   PAEP   PI3 PITX1  PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:15.442] [bbotk]  FALSE FALSE FALSE TRUE  FALSE FALSE  TRUE FALSE FALSE         TRUE 
    ## INFO  [16:00:15.442] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:15.442] [bbotk]          FALSE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:15.442] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:15.442] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:15.442] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:15.442] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:15.442] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:15.442] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:15.442] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:15.442] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE   TRUE TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:15.442] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [16:00:15.442] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE TRUE TRUE 
    ## INFO  [16:00:15.442] [bbotk]                                           features classif.bacc 
    ## INFO  [16:00:15.442] [bbotk]  AC007879.6,AC116614.1,CILP2,COL11A1,CPNE7,EN2,...    0.6717172 
    ## INFO  [16:00:15.569] [mlr3]  Applying learner 'classif.svm.fselector' on task 'kirc_cla' (iter 1/3) 
    ## INFO  [16:00:15.815] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:15.818] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:16.614] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:16.626] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:16.717] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:16.794] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:16.883] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:16.976] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:17.048] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:17.152] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:17.245] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:17.324] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:17.437] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:17.526] [mlr3]  Finished benchmark 
    ## INFO  [16:00:18.301] [bbotk] Result of batch 1: 
    ## INFO  [16:00:18.311] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:18.311] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]        TRUE       TRUE      FALSE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]       FALSE      FALSE       TRUE       TRUE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:18.311] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE       FALSE  TRUE FALSE    FALSE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE    FALSE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:18.311] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE   TRUE  TRUE   FALSE FALSE  TRUE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE   TRUE  TRUE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE  TRUE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  FALSE FALSE    TRUE FALSE  TRUE  TRUE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:18.311] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]     TRUE   TRUE FALSE         TRUE  FALSE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]     TRUE  FALSE FALSE        FALSE   TRUE   TRUE FALSE  TRUE  TRUE FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]    FALSE  FALSE  TRUE         TRUE  FALSE   TRUE FALSE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE  FALSE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]    FALSE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]    FALSE  FALSE FALSE         TRUE   TRUE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]    FALSE   TRUE FALSE        FALSE   TRUE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE    TRUE FALSE  TRUE FALSE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE  TRUE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:18.311] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:18.311] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:18.311] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:18.311] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:18.311] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:18.311] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:18.311] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:18.311] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:18.311] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:18.311] [bbotk]     FALSE FALSE FALSE     TRUE    TRUE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:18.311] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:18.311] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE     TRUE  TRUE FALSE FALSE   FALSE FALSE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:18.311] [bbotk]  FALSE      TRUE      TRUE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:18.311] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:18.311] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE   TRUE FALSE FALSE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE  FALSE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE  FALSE  TRUE FALSE FALSE  FALSE  FALSE   TRUE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:18.311] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE  TRUE FALSE  FALSE  TRUE FALSE  TRUE FALSE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE FALSE FALSE  FALSE  TRUE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:18.311] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:18.311] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:18.311] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE         TRUE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        FALSE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE        FALSE         TRUE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE         TRUE        FALSE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE        FALSE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:18.311] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:18.311] [bbotk]       TRUE  TRUE FALSE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:18.311] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:18.311] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]      FALSE FALSE FALSE FALSE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]      FALSE FALSE  TRUE  TRUE    TRUE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:18.311] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:18.311] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE   FALSE    TRUE    TRUE FALSE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE    TRUE   FALSE    TRUE FALSE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]   FALSE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE  FALSE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:18.311] [bbotk]   FALSE    TRUE   FALSE   FALSE  TRUE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:18.311] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6756098 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE  TRUE    0.6756098 
    ## INFO  [16:00:18.311] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6756098 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7006098 
    ## INFO  [16:00:18.311] [bbotk]  FALSE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE    0.6628049 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE    0.7134146 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6884146 
    ## INFO  [16:00:18.311] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE FALSE    0.6756098 
    ## INFO  [16:00:18.311] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE FALSE  TRUE  TRUE  TRUE    0.7006098 
    ## INFO  [16:00:18.311] [bbotk]   TRUE FALSE  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE  TRUE    0.6756098 
    ## INFO  [16:00:18.311] [bbotk]                                 uhash 
    ## INFO  [16:00:18.311] [bbotk]  f8b94add-1269-43df-b92b-84f17c88a146 
    ## INFO  [16:00:18.311] [bbotk]  641d53e6-8ff2-4cc3-9b71-ca93b2f74e4b 
    ## INFO  [16:00:18.311] [bbotk]  6a332cb2-f137-4d48-9ad6-2cc2dd86d5d1 
    ## INFO  [16:00:18.311] [bbotk]  c3258abe-dcda-421f-9483-5630ecf3c74d 
    ## INFO  [16:00:18.311] [bbotk]  0507adc3-11b3-4b97-9f14-fd62729c68e7 
    ## INFO  [16:00:18.311] [bbotk]  ec373a99-37a8-4f6b-adc6-23df23b10dcc 
    ## INFO  [16:00:18.311] [bbotk]  bf4c7996-222d-4487-9dad-2a8603afd9bc 
    ## INFO  [16:00:18.311] [bbotk]  cf6beb2a-e15c-4e24-b6d6-fd035e3bb1b9 
    ## INFO  [16:00:18.311] [bbotk]  32ef0ec7-86ba-4304-bcde-60c8763a2cca 
    ## INFO  [16:00:18.311] [bbotk]  655fc2c3-71a2-4efe-9479-15dddcd35e46 
    ## INFO  [16:00:18.313] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:18.871] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:18.878] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:18.941] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.021] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.084] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.156] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.215] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.292] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.354] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.419] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.492] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:19.566] [mlr3]  Finished benchmark 
    ## INFO  [16:00:19.879] [bbotk] Result of batch 2: 
    ## INFO  [16:00:19.889] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:19.889] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]       FALSE      FALSE      FALSE       TRUE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]        TRUE       TRUE      FALSE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:19.889] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE        TRUE FALSE FALSE     TRUE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE       FALSE  TRUE FALSE    FALSE    FALSE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:19.889] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  FALSE FALSE    TRUE FALSE  TRUE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE   TRUE  TRUE   FALSE FALSE FALSE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE   TRUE FALSE   FALSE  TRUE  TRUE  TRUE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:19.889] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE   TRUE  TRUE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE   TRUE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]     TRUE   TRUE FALSE        FALSE   TRUE   TRUE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]    FALSE  FALSE  TRUE         TRUE  FALSE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE  FALSE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE FALSE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE FALSE   FALSE  TRUE  TRUE FALSE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE FALSE    TRUE FALSE FALSE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:19.889] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:19.889] [bbotk]      TRUE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE  FALSE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:19.889] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:19.889] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]      TRUE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:19.889] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:19.889] [bbotk]   TRUE     TRUE FALSE  TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE     TRUE FALSE FALSE FALSE    TRUE FALSE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE     TRUE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:19.889] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE      TRUE     FALSE     FALSE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:19.889] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:19.889] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE   TRUE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE   TRUE FALSE  TRUE FALSE  FALSE  FALSE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE   TRUE  FALSE   TRUE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   TRUE FALSE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE   TRUE   TRUE  FALSE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:19.889] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE FALSE FALSE  FALSE  TRUE FALSE FALSE FALSE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  TRUE FALSE   TRUE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE FALSE FALSE FALSE   TRUE FALSE FALSE  TRUE FALSE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  TRUE  TRUE FALSE  FALSE  TRUE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:19.889] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         FALSE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:19.889] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE        FALSE        TRUE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:19.889] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:19.889] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:19.889] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:19.889] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:19.889] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:19.889] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:19.889] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:19.889] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:19.889] [bbotk]       TRUE FALSE FALSE  TRUE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:19.889] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:19.889] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:19.889] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:19.889] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:19.889] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:19.889] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE   FALSE   FALSE    TRUE  TRUE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:19.889] [bbotk]    TRUE   FALSE    TRUE   FALSE  TRUE   TRUE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:19.889] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:19.889] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6884146 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  FALSE   TRUE  TRUE   FALSE  TRUE FALSE FALSE  TRUE    0.6634146 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE   FALSE  TRUE  TRUE FALSE  TRUE    0.6878049 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.6756098 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.5128049 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE    0.7006098 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE FALSE  TRUE  TRUE    0.6634146 
    ## INFO  [16:00:19.889] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6756098 
    ## INFO  [16:00:19.889] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6256098 
    ## INFO  [16:00:19.889] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE  TRUE    0.7128049 
    ## INFO  [16:00:19.889] [bbotk]                                 uhash 
    ## INFO  [16:00:19.889] [bbotk]  dd965303-12f8-4603-a1fe-a86fcdaa9d0f 
    ## INFO  [16:00:19.889] [bbotk]  4e2f4b98-7b9b-4e65-8226-0209201a9ea5 
    ## INFO  [16:00:19.889] [bbotk]  d427dbc9-ce81-4510-9636-63f54522717c 
    ## INFO  [16:00:19.889] [bbotk]  4ce40f24-496c-481e-b3d4-18d33d15bd9c 
    ## INFO  [16:00:19.889] [bbotk]  dee10f2c-0ce1-4e7a-acaf-fd7e2638f1e0 
    ## INFO  [16:00:19.889] [bbotk]  3dd5bccb-395e-4505-a4a8-c59a262d5227 
    ## INFO  [16:00:19.889] [bbotk]  6f27b1ac-7679-4dd3-8253-9ae04372584b 
    ## INFO  [16:00:19.889] [bbotk]  1be253a7-0a54-403d-ac92-ffa1f3e4959f 
    ## INFO  [16:00:19.889] [bbotk]  1a0561a7-c88a-49bb-b45c-b4ce0cdbf76a 
    ## INFO  [16:00:19.889] [bbotk]  4b86e715-2520-42e9-b091-c72e7c31e270 
    ## INFO  [16:00:19.896] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:19.897] [bbotk] Result: 
    ## INFO  [16:00:19.903] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [16:00:19.903] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE TRUE TRUE 
    ## INFO  [16:00:19.903] [bbotk]   AMH APCDD1L_AS1  AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:19.903] [bbotk]  TRUE        TRUE FALSE TRUE     TRUE     TRUE FALSE TRUE     TRUE      TRUE 
    ## INFO  [16:00:19.903] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [16:00:19.903] [bbotk]  TRUE   TRUE  TRUE    TRUE TRUE  TRUE  TRUE         FALSE  FALSE FALSE TRUE 
    ## INFO  [16:00:19.903] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2 DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:19.903] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE  FALSE  TRUE TRUE TRUE  TRUE   TRUE 
    ## INFO  [16:00:19.903] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:19.903] [bbotk]   TRUE TRUE    TRUE  TRUE  TRUE TRUE   TRUE     TRUE      TRUE   TRUE     TRUE 
    ## INFO  [16:00:19.903] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:19.903] [bbotk]  FALSE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE    FALSE FALSE 
    ## INFO  [16:00:19.903] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15  KNG1 KRT7 L1CAM LECT1 
    ## INFO  [16:00:19.903] [bbotk]      TRUE  TRUE FALSE  TRUE    TRUE  TRUE TRUE  TRUE FALSE TRUE FALSE  TRUE 
    ## INFO  [16:00:19.903] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:19.903] [bbotk]       TRUE      TRUE      TRUE     FALSE      TRUE     FALSE   TRUE   TRUE 
    ## INFO  [16:00:19.903] [bbotk]  MAGEC3 MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3  PAEP  PI3 
    ## INFO  [16:00:19.903] [bbotk]    TRUE TRUE TRUE TRUE  FALSE   TRUE   TRUE FALSE   TRUE FALSE FALSE FALSE TRUE 
    ## INFO  [16:00:19.903] [bbotk]  PITX1  PLG PRR15L  PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:19.903] [bbotk]   TRUE TRUE   TRUE FALSE  TRUE  TRUE TRUE         TRUE          TRUE 
    ## INFO  [16:00:19.903] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:19.903] [bbotk]          FALSE         TRUE          TRUE          TRUE          TRUE 
    ## INFO  [16:00:19.903] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:19.903] [bbotk]          TRUE        TRUE          TRUE         TRUE         TRUE         FALSE 
    ## INFO  [16:00:19.903] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1 SAA1 SAA2_SAA4  SBSN 
    ## INFO  [16:00:19.903] [bbotk]          TRUE        FALSE         TRUE       FALSE TRUE TRUE     FALSE FALSE 
    ## INFO  [16:00:19.903] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:19.903] [bbotk]  FALSE FALSE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE  FALSE    TRUE 
    ## INFO  [16:00:19.903] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [16:00:19.903] [bbotk]     TRUE    TRUE  TRUE   TRUE TRUE    TRUE   TRUE     FALSE  TRUE TRUE  FALSE 
    ## INFO  [16:00:19.903] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5  ZIC2 ZIC5 
    ## INFO  [16:00:19.903] [bbotk]    TRUE  TRUE    TRUE TRUE  TRUE FALSE TRUE 
    ## INFO  [16:00:19.903] [bbotk]                                                        features classif.bacc 
    ## INFO  [16:00:19.903] [bbotk]  AC003092.1,AC006262.4,AC007879.6,AC104654.2,AC116614.1,AFM,...    0.7134146 
    ## INFO  [16:00:19.955] [mlr3]  Applying learner 'classif.xgboost.fselector' on task 'kirc_cla' (iter 3/3) 
    ## INFO  [16:00:20.182] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:20.184] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:21.336] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:21.343] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.405] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.464] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.529] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.590] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.651] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.717] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.792] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.851] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.915] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:21] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:21.981] [mlr3]  Finished benchmark 
    ## INFO  [16:00:22.472] [bbotk] Result of batch 1: 
    ## INFO  [16:00:22.483] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:22.483] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]        TRUE       TRUE       TRUE      FALSE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]        TRUE      FALSE      FALSE       TRUE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]       FALSE      FALSE       TRUE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:22.483] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE       FALSE FALSE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE       FALSE  TRUE  TRUE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:22.483] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE   TRUE  TRUE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE FALSE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  FALSE FALSE   FALSE  TRUE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:22.483] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]     TRUE  FALSE  TRUE         TRUE  FALSE   TRUE FALSE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]     TRUE  FALSE  TRUE        FALSE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE  TRUE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  TRUE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE FALSE   FALSE FALSE  TRUE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:22.483] [bbotk]      TRUE  TRUE  TRUE    FALSE   FALSE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:22.483] [bbotk]     FALSE FALSE  TRUE     TRUE   FALSE   TRUE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:22.483] [bbotk]      TRUE  TRUE FALSE     TRUE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:22.483] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:22.483] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:22.483] [bbotk]      TRUE FALSE  TRUE     TRUE   FALSE   TRUE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:22.483] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [16:00:22.483] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:22.483] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:22.483] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:22.483] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:22.483] [bbotk]   TRUE    FALSE FALSE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE   FALSE FALSE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE   FALSE FALSE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE    FALSE  TRUE  TRUE FALSE   FALSE FALSE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:22.483] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE     FALSE      TRUE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:22.483] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:22.483] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE  FALSE  TRUE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE   TRUE  TRUE FALSE FALSE  FALSE  FALSE   TRUE FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE   TRUE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE   TRUE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE  FALSE  FALSE   TRUE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  TRUE  TRUE FALSE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE FALSE FALSE FALSE  FALSE  TRUE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE FALSE FALSE   TRUE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE FALSE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:22.483] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:22.483] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE        TRUE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:22.483] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:22.483] [bbotk]          FALSE         TRUE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:22.483] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:22.483] [bbotk]       TRUE  TRUE FALSE FALSE    TRUE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:22.483] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:22.483] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:22.483] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:22.483] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:22.483] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:22.483] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:22.483] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:22.483] [bbotk]       TRUE FALSE FALSE FALSE    TRUE    TRUE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:22.483] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:22.483] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:22.483] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE    TRUE   FALSE   FALSE FALSE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE   FALSE   FALSE   FALSE  TRUE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:22.483] [bbotk]    TRUE   FALSE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:22.483] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:22.483] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE    0.5788690 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE    0.5997024 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE FALSE    0.6547619 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6413690 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE    TRUE  TRUE  TRUE FALSE  TRUE    0.6101190 
    ## INFO  [16:00:22.483] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE   FALSE  TRUE  TRUE  TRUE FALSE    0.6562500 
    ## INFO  [16:00:22.483] [bbotk]   TRUE FALSE   TRUE  FALSE  TRUE    TRUE  TRUE FALSE FALSE FALSE    0.6666667 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6190476 
    ## INFO  [16:00:22.483] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE    TRUE FALSE FALSE  TRUE  TRUE    0.7842262 
    ## INFO  [16:00:22.483] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6562500 
    ## INFO  [16:00:22.483] [bbotk]                                 uhash 
    ## INFO  [16:00:22.483] [bbotk]  aa0acbc6-a426-4b6a-abbd-12711023de0f 
    ## INFO  [16:00:22.483] [bbotk]  b0e426ad-bd65-4cd8-aca9-f1b3d7a613a2 
    ## INFO  [16:00:22.483] [bbotk]  494b187d-ab99-4c36-bcf4-3545401f6cbb 
    ## INFO  [16:00:22.483] [bbotk]  0aa4b70c-7944-410f-a7f3-6dc7368ae5f6 
    ## INFO  [16:00:22.483] [bbotk]  fd8ee2c3-e047-4d04-bc4b-13000da4bdda 
    ## INFO  [16:00:22.483] [bbotk]  74d0f13b-3fc0-4ae7-b109-3fcb96ab3e4d 
    ## INFO  [16:00:22.483] [bbotk]  88203f15-9a9a-4586-9d14-2501a03bcc08 
    ## INFO  [16:00:22.483] [bbotk]  44643c37-8e08-4bb9-a71f-f54216e39f8b 
    ## INFO  [16:00:22.483] [bbotk]  0384b43a-c512-4b2b-b60b-22288ac4e902 
    ## INFO  [16:00:22.483] [bbotk]  0e88981b-704d-4244-9466-2b380f9012fd 
    ## INFO  [16:00:22.485] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:23.632] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:23.640] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:23] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:23.722] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:23] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:23.786] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:23] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:23.847] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:23] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:23.914] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:23] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:23.981] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:24] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:24.045] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:24] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:24.113] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:24] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:24.180] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:24] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:24.257] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:24] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:24.322] [mlr3]  Finished benchmark 
    ## INFO  [16:00:24.818] [bbotk] Result of batch 2: 
    ## INFO  [16:00:24.828] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:24.828] [bbotk]        TRUE       TRUE       TRUE      FALSE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]        TRUE       TRUE       TRUE      FALSE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]       FALSE      FALSE       TRUE       TRUE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]        TRUE       TRUE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:24.828] [bbotk]  FALSE        TRUE  TRUE  TRUE    FALSE    FALSE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE    FALSE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  FALSE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE  TRUE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:24.828] [bbotk]    FALSE   TRUE FALSE         TRUE  FALSE   TRUE  TRUE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]     TRUE  FALSE  TRUE         TRUE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE   TRUE  TRUE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE    TRUE FALSE FALSE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE    TRUE FALSE FALSE FALSE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:24.828] [bbotk]     FALSE FALSE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:24.828] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:24.828] [bbotk]      TRUE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:24.828] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:24.828] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:24.828] [bbotk]  FALSE    FALSE FALSE  TRUE  TRUE    TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE   FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE   FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE     TRUE  TRUE FALSE  TRUE   FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:24.828] [bbotk]  FALSE      TRUE      TRUE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE     FALSE      TRUE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE     FALSE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:24.828] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:24.828] [bbotk]    TRUE  FALSE FALSE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   TRUE  TRUE FALSE FALSE   TRUE  FALSE   TRUE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE   TRUE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   TRUE FALSE  TRUE FALSE  FALSE  FALSE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE  TRUE FALSE   TRUE  TRUE FALSE FALSE  TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE FALSE FALSE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE  TRUE FALSE   TRUE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE  TRUE FALSE  FALSE FALSE FALSE  TRUE FALSE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:24.828] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE        FALSE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         TRUE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:24.828] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE         TRUE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE        FALSE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:24.828] [bbotk]           TRUE         TRUE        FALSE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [16:00:24.828] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:24.828] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:24.828] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:24.828] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:24.828] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:24.828] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:24.828] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:24.828] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:24.828] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:24.828] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:24.828] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:24.828] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE  FALSE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE    TRUE    TRUE   FALSE FALSE  FALSE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:24.828] [bbotk]   FALSE   FALSE    TRUE    TRUE  TRUE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:24.828] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:24.828] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE   FALSE  TRUE FALSE FALSE  TRUE    0.6041667 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE  TRUE    0.6354167 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE  TRUE FALSE  TRUE  TRUE    0.6354167 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE  FALSE   TRUE FALSE    TRUE FALSE FALSE  TRUE FALSE    0.5282738 
    ## INFO  [16:00:24.828] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.5639881 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.4925595 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE FALSE  TRUE FALSE  TRUE    0.6041667 
    ## INFO  [16:00:24.828] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7380952 
    ## INFO  [16:00:24.828] [bbotk]   TRUE FALSE   TRUE   TRUE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6413690 
    ## INFO  [16:00:24.828] [bbotk]  FALSE FALSE   TRUE   TRUE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE    0.6711310 
    ## INFO  [16:00:24.828] [bbotk]                                 uhash 
    ## INFO  [16:00:24.828] [bbotk]  7a1b0752-95fb-44da-8b14-4dcb10673fb3 
    ## INFO  [16:00:24.828] [bbotk]  c3d339ee-90f2-4db4-a35c-9c5978b0e028 
    ## INFO  [16:00:24.828] [bbotk]  d2c6025f-22cc-4aeb-af26-239a1c15aa82 
    ## INFO  [16:00:24.828] [bbotk]  3fd5cc3a-8a5a-49b4-bc29-72f1b30137e7 
    ## INFO  [16:00:24.828] [bbotk]  e33156d5-cdeb-4c37-aebb-cd9c41d94c10 
    ## INFO  [16:00:24.828] [bbotk]  ba6b4e27-cb61-4be0-8cc7-add5299a5f5c 
    ## INFO  [16:00:24.828] [bbotk]  140d4df2-a004-4536-9073-adfee280309d 
    ## INFO  [16:00:24.828] [bbotk]  62a14005-c635-45a0-8a10-688eaf618bb8 
    ## INFO  [16:00:24.828] [bbotk]  d7ed5537-7a53-49a5-b0ae-03757882e740 
    ## INFO  [16:00:24.828] [bbotk]  19e621ba-dc8d-4289-8f14-adc134551b47 
    ## INFO  [16:00:24.836] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:24.837] [bbotk] Result: 
    ## INFO  [16:00:24.843] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [16:00:24.843] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE TRUE TRUE 
    ## INFO  [16:00:24.843] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:24.843] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:24.843] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:24.843] [bbotk]  TRUE  FALSE FALSE   FALSE TRUE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:24.843] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:24.843] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:24.843] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:24.843] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:24.843] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:24.843] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:24.843] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15 KNG1 KRT7 L1CAM 
    ## INFO  [16:00:24.843] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE    TRUE  TRUE TRUE FALSE TRUE TRUE  TRUE 
    ## INFO  [16:00:24.843] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:24.843] [bbotk]  FALSE     FALSE      TRUE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:24.843] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:24.843] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:24.843] [bbotk]   PAEP   PI3 PITX1  PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:24.843] [bbotk]  FALSE FALSE  TRUE TRUE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:24.843] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:24.843] [bbotk]           TRUE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:24.843] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:24.843] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:24.843] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:24.843] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:24.843] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:24.843] [bbotk]       TRUE FALSE FALSE FALSE    TRUE    TRUE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:24.843] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:24.843] [bbotk]    TRUE   FALSE   FALSE   FALSE  TRUE  FALSE TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:24.843] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [16:00:24.843] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE    TRUE FALSE FALSE TRUE TRUE 
    ## INFO  [16:00:24.843] [bbotk]                                features classif.bacc 
    ## INFO  [16:00:24.843] [bbotk]  AC116614.1,AFM,AHSG,CA1,CHAT,CIDEC,...    0.7842262 
    ## [16:00:24] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:24.903] [mlr3]  Applying learner 'classif.xgboost.fselector' on task 'kirc_cla' (iter 1/3) 
    ## INFO  [16:00:25.080] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:25.082] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:26.358] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:26.366] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.429] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.495] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.574] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.659] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.723] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.782] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.844] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.911] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:26] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:26.979] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:27] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:27.045] [mlr3]  Finished benchmark 
    ## INFO  [16:00:27.599] [bbotk] Result of batch 1: 
    ## INFO  [16:00:27.610] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:27.610] [bbotk]        TRUE       TRUE       TRUE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]       FALSE      FALSE      FALSE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]       FALSE       TRUE      FALSE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:27.610] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE       FALSE  TRUE  TRUE    FALSE    FALSE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE    FALSE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  FALSE  TRUE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:27.610] [bbotk]     TRUE  FALSE  TRUE        FALSE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]    FALSE   TRUE FALSE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]    FALSE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  TRUE   FALSE FALSE  TRUE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  TRUE   FALSE  TRUE FALSE  TRUE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE   FALSE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:27.610] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:27.610] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:27.610] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:27.610] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:27.610] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:27.610] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:27.610] [bbotk]      TRUE  TRUE  TRUE    FALSE   FALSE   TRUE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:27.610] [bbotk]     FALSE  TRUE FALSE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:27.610] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:27.610] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:27.610] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE     TRUE FALSE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE    TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE   FALSE  TRUE FALSE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     TRUE FALSE  TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:27.610] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     FALSE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     FALSE     FALSE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:27.610] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:27.610] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE  FALSE  FALSE   TRUE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE   TRUE FALSE  TRUE  TRUE   TRUE  FALSE   TRUE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE  FALSE  TRUE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE   TRUE   TRUE  FALSE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE  FALSE   TRUE  FALSE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE   TRUE FALSE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:27.610] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE FALSE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE  TRUE FALSE  FALSE  TRUE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:27.610] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:27.610] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE        FALSE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:27.610] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:27.610] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE        FALSE         TRUE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:27.610] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:27.610] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:27.610] [bbotk]      FALSE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:27.610] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:27.610] [bbotk]      FALSE FALSE  TRUE  TRUE    TRUE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:27.610] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:27.610] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:27.610] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:27.610] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:27.610] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:27.610] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:27.610] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:27.610] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:27.610] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE  FALSE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]   FALSE   FALSE    TRUE    TRUE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:27.610] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:27.610] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6858974 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE    0.6931090 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE   FALSE FALSE FALSE  TRUE  TRUE    0.7315705 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6266026 
    ## INFO  [16:00:27.610] [bbotk]   TRUE FALSE  FALSE   TRUE FALSE    TRUE FALSE FALSE FALSE FALSE    0.6650641 
    ## INFO  [16:00:27.610] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5633013 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE    0.5745192 
    ## INFO  [16:00:27.610] [bbotk]  FALSE  TRUE  FALSE  FALSE  TRUE   FALSE  TRUE FALSE FALSE FALSE    0.4975962 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE FALSE    0.6474359 
    ## INFO  [16:00:27.610] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6266026 
    ## INFO  [16:00:27.610] [bbotk]                                 uhash 
    ## INFO  [16:00:27.610] [bbotk]  8d05dcf7-a27f-4877-b245-f45d4cc7a3b8 
    ## INFO  [16:00:27.610] [bbotk]  1bf08963-ed24-404f-b104-1977dafda9e6 
    ## INFO  [16:00:27.610] [bbotk]  7aeb50e6-874b-437e-9ba2-98d6c39ebb89 
    ## INFO  [16:00:27.610] [bbotk]  3d051331-b9ee-4847-b5f9-4f212f4cd240 
    ## INFO  [16:00:27.610] [bbotk]  4a98e29a-da36-4184-8fa6-b5bedc31df41 
    ## INFO  [16:00:27.610] [bbotk]  433c4041-6d13-4ea2-928d-8a2976239050 
    ## INFO  [16:00:27.610] [bbotk]  7664ba86-6db2-4ebd-8f68-e60978b8a95f 
    ## INFO  [16:00:27.610] [bbotk]  23f6cb7f-175e-4770-a862-ce2f86fc0188 
    ## INFO  [16:00:27.610] [bbotk]  aabb34ca-c787-440a-bec7-d83f68ca89cc 
    ## INFO  [16:00:27.610] [bbotk]  e68c7978-82ff-4c3d-9ab9-88d5ec28e9c7 
    ## INFO  [16:00:27.612] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:28.906] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:28.915] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:28] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:28.983] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.047] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.113] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.176] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.237] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.319] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.394] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.460] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.527] [mlr3]  Applying learner 'select.classif.xgboost' on task 'kirc_cla' (iter 1/1) 
    ## [16:00:29] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:29.596] [mlr3]  Finished benchmark 
    ## INFO  [16:00:30.154] [bbotk] Result of batch 2: 
    ## INFO  [16:00:30.165] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:30.165] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE      FALSE      FALSE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE       TRUE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:30.165] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE       FALSE FALSE  TRUE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE       FALSE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:30.165] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE  TRUE  TRUE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE  TRUE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:30.165] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE FALSE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE FALSE   FALSE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE   FALSE FALSE  TRUE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:30.165] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:30.165] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:30.165] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:30.165] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:30.165] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:30.165] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:30.165] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:30.165] [bbotk]      TRUE FALSE FALSE    FALSE    TRUE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:30.165] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:30.165] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:30.165] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:30.165] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE     TRUE FALSE  TRUE  TRUE   FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE   FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:30.165] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE     FALSE      TRUE      TRUE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE      TRUE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:30.165] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:30.165] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE  TRUE  TRUE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE FALSE  TRUE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:30.165] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:30.165] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:30.165] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:30.165] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE        FALSE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE        FALSE         TRUE        FALSE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:30.165] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:30.165] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:30.165] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:30.165] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:30.165] [bbotk]       TRUE  TRUE FALSE FALSE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:30.165] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE   FALSE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:30.165] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:30.165] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]   FALSE    TRUE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:30.165] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:30.165] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.6370192 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE  TRUE    0.5881410 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE FALSE    0.6826923 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE  TRUE    0.6442308 
    ## INFO  [16:00:30.165] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.5849359 
    ## INFO  [16:00:30.165] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6474359 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7139423 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE FALSE FALSE  TRUE    0.6754808 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6057692 
    ## INFO  [16:00:30.165] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE    0.6858974 
    ## INFO  [16:00:30.165] [bbotk]                                 uhash 
    ## INFO  [16:00:30.165] [bbotk]  e0f29548-1c56-4590-ac95-a760bfd2154f 
    ## INFO  [16:00:30.165] [bbotk]  ea3d452e-6ff2-4ce5-a48d-44a1dd22e48c 
    ## INFO  [16:00:30.165] [bbotk]  3c5d646c-27d4-41d6-96b8-8ddca4cdfa03 
    ## INFO  [16:00:30.165] [bbotk]  3d75b328-92b9-4e29-a2f5-4d276a5a3e4b 
    ## INFO  [16:00:30.165] [bbotk]  20e8ea8b-0cb6-4093-99b3-b72d25272837 
    ## INFO  [16:00:30.165] [bbotk]  2d7eed56-e120-4761-b485-be947bc38573 
    ## INFO  [16:00:30.165] [bbotk]  46ea533d-e49c-458a-a5d8-81d25f8bac16 
    ## INFO  [16:00:30.165] [bbotk]  73437002-b2da-450f-ab39-321c4fba341f 
    ## INFO  [16:00:30.165] [bbotk]  a05b491c-4efd-4a1a-96e4-d35342aa414c 
    ## INFO  [16:00:30.165] [bbotk]  c863ee0f-c72b-4fa7-a068-f08f305de7a3 
    ## INFO  [16:00:30.172] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:30.173] [bbotk] Result: 
    ## INFO  [16:00:30.179] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM AHSG 
    ## INFO  [16:00:30.179] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE TRUE 
    ## INFO  [16:00:30.179] [bbotk]    AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:30.179] [bbotk]  FALSE       FALSE TRUE TRUE    FALSE    FALSE FALSE TRUE    FALSE     FALSE 
    ## INFO  [16:00:30.179] [bbotk]   CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [16:00:30.179] [bbotk]  TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE FALSE TRUE 
    ## INFO  [16:00:30.179] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:30.179] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE FALSE TRUE FALSE   TRUE 
    ## INFO  [16:00:30.179] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:30.179] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE TRUE   TRUE    FALSE      TRUE   TRUE    FALSE 
    ## INFO  [16:00:30.179] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:30.179] [bbotk]   TRUE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE FALSE 
    ## INFO  [16:00:30.179] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1 KRT7 L1CAM LECT1 
    ## INFO  [16:00:30.179] [bbotk]      TRUE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE  TRUE FALSE TRUE FALSE  TRUE 
    ## INFO  [16:00:30.179] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:30.179] [bbotk]      FALSE     FALSE      TRUE      TRUE     FALSE     FALSE   TRUE   TRUE 
    ## INFO  [16:00:30.179] [bbotk]  MAGEC3  MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3  PAEP 
    ## INFO  [16:00:30.179] [bbotk]    TRUE FALSE TRUE TRUE   TRUE  FALSE   TRUE FALSE   TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:30.179] [bbotk]   PI3 PITX1   PLG PRR15L PSG9 PVALB RAB25  RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:30.179] [bbotk]  TRUE FALSE FALSE  FALSE TRUE FALSE  TRUE FALSE        FALSE         FALSE 
    ## INFO  [16:00:30.179] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:30.179] [bbotk]          FALSE         TRUE          TRUE         FALSE          TRUE 
    ## INFO  [16:00:30.179] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:30.179] [bbotk]          TRUE       FALSE          TRUE         TRUE         TRUE          TRUE 
    ## INFO  [16:00:30.179] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1  SAA1 SAA2_SAA4  SBSN 
    ## INFO  [16:00:30.179] [bbotk]         FALSE        FALSE        FALSE        TRUE TRUE FALSE     FALSE FALSE 
    ## INFO  [16:00:30.179] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:30.179] [bbotk]   TRUE  TRUE    TRUE   FALSE   FALSE   FALSE    TRUE   FALSE  FALSE    TRUE 
    ## INFO  [16:00:30.179] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [16:00:30.179] [bbotk]     TRUE   FALSE  TRUE  FALSE TRUE   FALSE   TRUE     FALSE  TRUE TRUE   TRUE 
    ## INFO  [16:00:30.179] [bbotk]  TUBBP6 TUNAR UGT1A10  UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [16:00:30.179] [bbotk]   FALSE  TRUE   FALSE FALSE FALSE TRUE TRUE 
    ## INFO  [16:00:30.179] [bbotk]                                                         features classif.bacc 
    ## INFO  [16:00:30.179] [bbotk]  AC003092.1,AC006262.4,AC006262.5,AC007879.6,AC104654.2,AHSG,...    0.7315705 
    ## [16:00:30] WARNING: amalgamation/../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'binary:logistic' was changed from 'error' to 'logloss'. Explicitly set eval_metric if you'd like to restore the old behavior.
    ## INFO  [16:00:30.232] [mlr3]  Applying learner 'classif.rpart.fselector' on task 'kirc_cla' (iter 3/3) 
    ## INFO  [16:00:30.394] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:30.396] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:30.926] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:30.934] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.033] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.098] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.180] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.267] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.365] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.440] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.539] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.624] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.694] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:31.817] [mlr3]  Finished benchmark 
    ## INFO  [16:00:32.185] [bbotk] Result of batch 1: 
    ## INFO  [16:00:32.197] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:32.197] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:32.197] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:32.197] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE        TRUE  TRUE  TRUE    FALSE    FALSE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE        TRUE FALSE  TRUE     TRUE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:32.197] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE   TRUE FALSE    TRUE FALSE  TRUE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:32.197] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:32.197] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE  FALSE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE  FALSE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:32.197] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:32.197] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]     FALSE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:32.197] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:32.197] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:32.197] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:32.197] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:32.197] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE   FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:32.197] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:32.197] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE      TRUE     FALSE      TRUE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:32.197] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:32.197] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE   TRUE FALSE FALSE  TRUE   TRUE  FALSE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   TRUE FALSE  TRUE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE FALSE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE  TRUE FALSE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:32.197] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:32.197] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE        FALSE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:32.197] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:32.197] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE        FALSE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:32.197] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE        FALSE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:32.197] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:32.197] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:32.197] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:32.197] [bbotk]       TRUE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:32.197] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:32.197] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:32.197] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:32.197] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:32.197] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:32.197] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:32.197] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:32.197] [bbotk]       TRUE  TRUE FALSE FALSE    TRUE    TRUE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:32.197] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:32.197] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE    TRUE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:32.197] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7166463 
    ## INFO  [16:00:32.197] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE    TRUE  TRUE  TRUE FALSE FALSE    0.6175031 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6848225 
    ## INFO  [16:00:32.197] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6493268 
    ## INFO  [16:00:32.197] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7050184 
    ## INFO  [16:00:32.197] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.5587515 
    ## INFO  [16:00:32.197] [bbotk]  FALSE FALSE  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE FALSE    0.6701346 
    ## INFO  [16:00:32.197] [bbotk]                                 uhash 
    ## INFO  [16:00:32.197] [bbotk]  66bd1341-f8ef-4067-9e22-3298fd02143c 
    ## INFO  [16:00:32.197] [bbotk]  61b300db-4ccb-4f40-b9b3-0f2fabf79ea0 
    ## INFO  [16:00:32.197] [bbotk]  439937c8-80b5-4e12-901b-beee54ccec68 
    ## INFO  [16:00:32.197] [bbotk]  092a2b32-5296-4e93-a419-9a62236b7ee5 
    ## INFO  [16:00:32.197] [bbotk]  5c31d222-94d5-4e34-b184-456c39da5239 
    ## INFO  [16:00:32.197] [bbotk]  01f490e0-ff4e-4206-b376-9172f3c08425 
    ## INFO  [16:00:32.197] [bbotk]  37c34c34-25a4-4238-9c09-07d18507dcf7 
    ## INFO  [16:00:32.197] [bbotk]  726ea0d5-8914-4915-99a3-cc0f6c43d082 
    ## INFO  [16:00:32.197] [bbotk]  e42f2aef-b08a-4222-b1d2-701c0ccc6da5 
    ## INFO  [16:00:32.197] [bbotk]  ace7cdf0-83b1-429e-9a67-1e2bd4aad415 
    ## INFO  [16:00:32.200] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:32.846] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:32.855] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:32.968] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.054] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.165] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.246] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.364] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.475] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.554] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.621] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.730] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:33.822] [mlr3]  Finished benchmark 
    ## INFO  [16:00:34.190] [bbotk] Result of batch 2: 
    ## INFO  [16:00:34.202] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:34.202] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]        TRUE      FALSE      FALSE       TRUE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:34.202] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE    FALSE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE       FALSE  TRUE FALSE    FALSE     TRUE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE       FALSE  TRUE  TRUE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:34.202] [bbotk]   TRUE   TRUE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE  TRUE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:34.202] [bbotk]     TRUE  FALSE  TRUE        FALSE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]    FALSE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE   TRUE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]    FALSE   TRUE FALSE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE FALSE   FALSE FALSE FALSE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:34.202] [bbotk]     FALSE  TRUE FALSE    FALSE    TRUE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:34.202] [bbotk]     FALSE FALSE FALSE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:34.202] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:34.202] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:34.202] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE FALSE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE    TRUE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE     TRUE  TRUE FALSE  TRUE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE FALSE  TRUE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE   FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:34.202] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE     FALSE     FALSE      TRUE     FALSE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE      TRUE      TRUE      TRUE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE     FALSE     FALSE      TRUE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:34.202] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:34.202] [bbotk]    TRUE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   TRUE  TRUE FALSE FALSE  FALSE   TRUE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE   TRUE FALSE FALSE  TRUE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:34.202] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE FALSE FALSE FALSE  FALSE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE FALSE  TRUE FALSE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE FALSE  TRUE FALSE  FALSE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE FALSE  TRUE   TRUE  TRUE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE FALSE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:34.202] [bbotk]          FALSE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:34.202] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:34.202] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE FALSE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE         TRUE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:34.202] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:34.202] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:34.202] [bbotk]       TRUE FALSE FALSE FALSE   FALSE    TRUE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:34.202] [bbotk]       TRUE  TRUE FALSE FALSE   FALSE    TRUE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:34.202] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE    TRUE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:34.202] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:34.202] [bbotk]       TRUE FALSE FALSE FALSE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:34.202] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:34.202] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:34.202] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE    TRUE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:34.202] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:34.202] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:34.202] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:34.202] [bbotk]   FALSE    TRUE    TRUE   FALSE FALSE   TRUE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE   FALSE    TRUE   FALSE FALSE   TRUE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   FALSE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:34.202] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:34.202] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:34.202] [bbotk]  FALSE  TRUE  FALSE   TRUE FALSE   FALSE FALSE  TRUE FALSE FALSE    0.6413709 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE   TRUE  FALSE  TRUE    TRUE FALSE FALSE  TRUE FALSE    0.6903305 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7050184 
    ## INFO  [16:00:34.202] [bbotk]   TRUE FALSE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.6787026 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6407589 
    ## INFO  [16:00:34.202] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7050184 
    ## INFO  [16:00:34.202] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6701346 
    ## INFO  [16:00:34.202] [bbotk]                                 uhash 
    ## INFO  [16:00:34.202] [bbotk]  c6abe1c1-221d-4db5-9f2d-0fe6e035003b 
    ## INFO  [16:00:34.202] [bbotk]  c1d0521a-634d-41c4-8f00-cf3e967519dd 
    ## INFO  [16:00:34.202] [bbotk]  5bacde87-1421-470b-b425-fd77eab456f8 
    ## INFO  [16:00:34.202] [bbotk]  96e8dd1b-ad79-404b-b1c4-4c57923c890e 
    ## INFO  [16:00:34.202] [bbotk]  f17e8655-e291-4b83-bf9d-d8bc7344b8f5 
    ## INFO  [16:00:34.202] [bbotk]  6ac45076-a206-47c0-af81-d4e06c1e787d 
    ## INFO  [16:00:34.202] [bbotk]  b71903ad-0ede-467d-af53-2db70bda55ca 
    ## INFO  [16:00:34.202] [bbotk]  e21f99e1-b58a-44d1-a5c3-87decff55590 
    ## INFO  [16:00:34.202] [bbotk]  aa034670-49af-45a0-9357-7ce4743acf12 
    ## INFO  [16:00:34.202] [bbotk]  99a8ac19-d2df-484f-9236-2aa09e2c6f65 
    ## INFO  [16:00:34.210] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:34.211] [bbotk] Result: 
    ## INFO  [16:00:34.218] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:34.218] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:34.218] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:34.218] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE FALSE TRUE    FALSE     FALSE 
    ## INFO  [16:00:34.218] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:34.218] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.218] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:34.218] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:34.218] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:34.218] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.218] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:34.218] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:34.218] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:34.218] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:34.218] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:34.218] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:34.218] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:34.218] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:34.218] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25 RHBG RP11_10O22.1 
    ## INFO  [16:00:34.218] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE  TRUE TRUE        FALSE 
    ## INFO  [16:00:34.218] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:34.218] [bbotk]           TRUE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:34.218] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:34.218] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:34.218] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1 SAA1 
    ## INFO  [16:00:34.218] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE TRUE 
    ## INFO  [16:00:34.218] [bbotk]  SAA2_SAA4 SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:34.218] [bbotk]      FALSE TRUE FALSE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:34.218] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:34.218] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:34.218] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 
    ## INFO  [16:00:34.218] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:34.218] [bbotk]                                         features classif.bacc 
    ## INFO  [16:00:34.218] [bbotk]  ATP6V0A4,ATP6V0D2,BSND,FDCSP,IGFN1,IGKV3_11,...    0.7166463 
    ## INFO  [16:00:34.251] [mlr3]  Applying learner 'classif.ksvm.fselector' on task 'kirc_cla' (iter 1/3) 
    ## INFO  [16:00:34.489] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:34.492] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:35.399] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:35.414] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:35.540] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:35.672] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:35.808] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:35.933] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:36.055] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:36.261] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:36.415] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:36.523] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:36.643] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:36.746] [mlr3]  Finished benchmark 
    ## INFO  [16:00:37.209] [bbotk] Result of batch 1: 
    ## INFO  [16:00:37.220] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:37.220] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]       FALSE       TRUE       TRUE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE       TRUE      FALSE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:37.220] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE    FALSE FALSE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE    FALSE  TRUE FALSE    FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE       FALSE FALSE  TRUE     TRUE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE        TRUE  TRUE  TRUE    FALSE    FALSE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE       FALSE  TRUE  TRUE    FALSE    FALSE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:37.220] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  FALSE  TRUE   FALSE  TRUE FALSE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:37.220] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]    FALSE  FALSE FALSE         TRUE   TRUE  FALSE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE  FALSE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]     TRUE   TRUE  TRUE        FALSE  FALSE  FALSE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  TRUE   FALSE  TRUE FALSE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE FALSE    TRUE FALSE FALSE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE FALSE   FALSE FALSE  TRUE  TRUE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:37.220] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:37.220] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:37.220] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:37.220] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:37.220] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:37.220] [bbotk]     FALSE  TRUE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:37.220] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:37.220] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:37.220] [bbotk]     FALSE  TRUE  TRUE    FALSE   FALSE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:37.220] [bbotk]      TRUE  TRUE FALSE     TRUE   FALSE   TRUE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [16:00:37.220] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:37.220] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE    TRUE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE     TRUE FALSE  TRUE FALSE   FALSE  TRUE FALSE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:37.220] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE      TRUE      TRUE     FALSE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE     FALSE     FALSE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:37.220] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:37.220] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE   TRUE FALSE FALSE FALSE   TRUE  FALSE  FALSE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE  FALSE  FALSE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE   TRUE  TRUE FALSE FALSE   TRUE  FALSE   TRUE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]   PAEP   PI3 PITX1  PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE  TRUE TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE FALSE  TRUE TRUE   TRUE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  TRUE  TRUE TRUE  FALSE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE FALSE FALSE TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  TRUE  TRUE TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE FALSE TRUE  FALSE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]  FALSE  TRUE FALSE TRUE   TRUE  TRUE  TRUE FALSE FALSE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE  TRUE TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE FALSE FALSE TRUE  FALSE  TRUE FALSE  TRUE FALSE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE  TRUE TRUE   TRUE  TRUE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:37.220] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE         FALSE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:37.220] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE        FALSE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE        TRUE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE        FALSE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:37.220] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:37.220] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]          FALSE        FALSE         TRUE         TRUE       FALSE FALSE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:37.220] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:37.220] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:37.220] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:37.220] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:37.220] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:37.220] [bbotk]      FALSE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:37.220] [bbotk]      FALSE  TRUE FALSE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:37.220] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE    TRUE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:37.220] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:37.220] [bbotk]       TRUE FALSE FALSE FALSE   FALSE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:37.220] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE    TRUE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:37.220] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:37.220] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:37.220] [bbotk]    TRUE   FALSE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:37.220] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7076389 
    ## INFO  [16:00:37.220] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.7701389 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE FALSE    0.6965278 
    ## INFO  [16:00:37.220] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6340278 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7076389 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE   FALSE  TRUE FALSE  TRUE  TRUE    0.7076389 
    ## INFO  [16:00:37.220] [bbotk]   TRUE FALSE   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE FALSE    0.7388889 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.7076389 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE    TRUE FALSE  TRUE FALSE FALSE    0.6652778 
    ## INFO  [16:00:37.220] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE FALSE  TRUE  TRUE    0.6763889 
    ## INFO  [16:00:37.220] [bbotk]                                 uhash 
    ## INFO  [16:00:37.220] [bbotk]  185bf3d5-d4a5-440c-8ee2-df98fbbd13cc 
    ## INFO  [16:00:37.220] [bbotk]  90e013ad-ab45-4de1-b5be-7efbab64fd3b 
    ## INFO  [16:00:37.220] [bbotk]  745c33ba-8366-4c75-80c6-a756e148742e 
    ## INFO  [16:00:37.220] [bbotk]  1eaa7f7c-570c-49bc-9541-17d2612d49bb 
    ## INFO  [16:00:37.220] [bbotk]  c8a1874a-3b8c-4b64-9474-4f26c6d6fb53 
    ## INFO  [16:00:37.220] [bbotk]  a63a300d-e95e-4b33-b4be-47ba78afc04e 
    ## INFO  [16:00:37.220] [bbotk]  7702bcb1-2d52-484c-b4ef-edfef5f8693f 
    ## INFO  [16:00:37.220] [bbotk]  3320f916-add0-4c60-aaf7-5a9fe0f9c63d 
    ## INFO  [16:00:37.220] [bbotk]  d449ef4a-c5ec-4682-bdef-386c0788b82b 
    ## INFO  [16:00:37.220] [bbotk]  c5dcc702-e941-430e-b115-954c29a4fcd4 
    ## INFO  [16:00:37.223] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:38.063] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:38.072] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.170] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.289] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.399] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.530] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.703] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.844] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:38.959] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:39.077] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:39.181] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:39.300] [mlr3]  Finished benchmark 
    ## INFO  [16:00:39.777] [bbotk] Result of batch 2: 
    ## INFO  [16:00:39.789] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:39.789] [bbotk]        TRUE       TRUE       TRUE      FALSE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]       FALSE       TRUE       TRUE      FALSE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]       FALSE      FALSE       TRUE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE       TRUE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]       FALSE      FALSE      FALSE       TRUE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:39.789] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE    FALSE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE        TRUE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE        TRUE  TRUE FALSE    FALSE     TRUE  TRUE FALSE    FALSE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE     TRUE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE   TRUE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE  TRUE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE  TRUE  TRUE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE FALSE  TRUE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:39.789] [bbotk]    FALSE  FALSE FALSE         TRUE   TRUE  FALSE  TRUE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE  FALSE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]     TRUE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]    FALSE   TRUE  TRUE        FALSE  FALSE  FALSE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]     TRUE  FALSE  TRUE        FALSE   TRUE  FALSE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE   TRUE FALSE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE   TRUE FALSE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE FALSE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE   FALSE  TRUE  TRUE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE    TRUE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  TRUE   FALSE  TRUE  TRUE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE   FALSE FALSE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE   FALSE  TRUE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:39.789] [bbotk]      TRUE FALSE  TRUE    FALSE    TRUE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:39.789] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]      TRUE FALSE  TRUE    FALSE    TRUE  FALSE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]      TRUE  TRUE  TRUE    FALSE   FALSE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]      TRUE  TRUE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]     FALSE  TRUE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:39.789] [bbotk]     FALSE  TRUE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:39.789] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:39.789] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE    TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE    TRUE FALSE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE     TRUE FALSE FALSE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:39.789] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE      TRUE      TRUE     FALSE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE      TRUE     FALSE     FALSE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:39.789] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:39.789] [bbotk]   FALSE  FALSE FALSE FALSE FALSE   TRUE   TRUE  FALSE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE  FALSE  TRUE FALSE FALSE   TRUE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE   TRUE  FALSE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE   TRUE FALSE FALSE FALSE   TRUE   TRUE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE  FALSE   TRUE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE   TRUE FALSE FALSE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE   TRUE FALSE FALSE FALSE   TRUE  FALSE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE FALSE FALSE FALSE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE  TRUE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE  TRUE FALSE  FALSE FALSE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE FALSE FALSE   TRUE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE FALSE FALSE  TRUE FALSE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE  TRUE  TRUE  FALSE FALSE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE FALSE FALSE  FALSE  TRUE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:39.789] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE          TRUE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:39.789] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:39.789] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE        FALSE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:39.789] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:39.789] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE         TRUE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         TRUE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:39.789] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:39.789] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:39.789] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:39.789] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]       TRUE FALSE FALSE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:39.789] [bbotk]       TRUE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]      FALSE FALSE  TRUE  TRUE   FALSE   FALSE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]      FALSE  TRUE  TRUE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:39.789] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:39.789] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:39.789] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:39.789] [bbotk]    TRUE    TRUE    TRUE   FALSE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]    TRUE   FALSE    TRUE   FALSE  TRUE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE    TRUE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:39.789] [bbotk]   FALSE    TRUE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:39.789] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE    0.7388889 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7187500 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE    TRUE FALSE  TRUE FALSE  TRUE    0.6763889 
    ## INFO  [16:00:39.789] [bbotk]  FALSE  TRUE  FALSE  FALSE  TRUE    TRUE FALSE FALSE FALSE FALSE    0.6854167 
    ## INFO  [16:00:39.789] [bbotk]   TRUE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE  TRUE FALSE FALSE    0.7388889 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE  FALSE   TRUE  TRUE    TRUE FALSE FALSE  TRUE FALSE    0.6763889 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE  TRUE  TRUE FALSE  TRUE    0.6562500 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE FALSE FALSE  TRUE FALSE    0.6854167 
    ## INFO  [16:00:39.789] [bbotk]  FALSE FALSE   TRUE   TRUE FALSE    TRUE FALSE FALSE  TRUE  TRUE    0.7500000 
    ## INFO  [16:00:39.789] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.6875000 
    ## INFO  [16:00:39.789] [bbotk]                                 uhash 
    ## INFO  [16:00:39.789] [bbotk]  73856ea7-f727-40a0-b380-c38baaf7e582 
    ## INFO  [16:00:39.789] [bbotk]  0aea4bce-ff41-43d7-87c0-f6884f8b3591 
    ## INFO  [16:00:39.789] [bbotk]  044b4e29-c0d4-4ece-a92a-7e7a6ec2f0da 
    ## INFO  [16:00:39.789] [bbotk]  ce17f2bd-1242-44ff-bdef-eddd5fff507c 
    ## INFO  [16:00:39.789] [bbotk]  c844e706-33c4-4dd1-a4df-2eef298a346c 
    ## INFO  [16:00:39.789] [bbotk]  cd0d8ef2-9157-4353-a570-6534ed19fa9d 
    ## INFO  [16:00:39.789] [bbotk]  9f240d33-0740-4731-8696-d102a1f104af 
    ## INFO  [16:00:39.789] [bbotk]  9e24b6cc-e682-4c57-8061-529bd8df732b 
    ## INFO  [16:00:39.789] [bbotk]  77f7cfeb-4652-4fa0-8a90-d76ce3bc4726 
    ## INFO  [16:00:39.789] [bbotk]  3161a162-67f0-44d2-a328-ce9f58e55c1d 
    ## INFO  [16:00:39.797] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:39.798] [bbotk] Result: 
    ## INFO  [16:00:39.805] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM AHSG 
    ## INFO  [16:00:39.805] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE TRUE 
    ## INFO  [16:00:39.805] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:39.805] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:39.805] [bbotk]    CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:39.805] [bbotk]  FALSE  FALSE FALSE   FALSE TRUE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:39.805] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:39.805] [bbotk]    FALSE  FALSE FALSE         TRUE   TRUE  FALSE  TRUE FALSE TRUE  TRUE   TRUE 
    ## INFO  [16:00:39.805] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:39.805] [bbotk]  FALSE TRUE   FALSE  TRUE FALSE FALSE   TRUE     TRUE     FALSE   TRUE    FALSE 
    ## INFO  [16:00:39.805] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:39.805] [bbotk]   TRUE  TRUE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE     TRUE FALSE 
    ## INFO  [16:00:39.805] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1 KRT7 L1CAM LECT1 
    ## INFO  [16:00:39.805] [bbotk]     FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE TRUE FALSE FALSE 
    ## INFO  [16:00:39.805] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:39.805] [bbotk]      FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE  FALSE 
    ## INFO  [16:00:39.805] [bbotk]  MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3  PAEP 
    ## INFO  [16:00:39.805] [bbotk]    TRUE FALSE FALSE FALSE   TRUE  FALSE  FALSE FALSE  FALSE TRUE FALSE FALSE 
    ## INFO  [16:00:39.805] [bbotk]    PI3 PITX1  PLG PRR15L  PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:39.805] [bbotk]  FALSE  TRUE TRUE   TRUE FALSE FALSE  TRUE TRUE        FALSE         FALSE 
    ## INFO  [16:00:39.805] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:39.805] [bbotk]          FALSE        FALSE         FALSE         FALSE         FALSE 
    ## INFO  [16:00:39.805] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:39.805] [bbotk]         FALSE        TRUE          TRUE        FALSE        FALSE          TRUE 
    ## INFO  [16:00:39.805] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1  SAA1 SAA2_SAA4  SBSN 
    ## INFO  [16:00:39.805] [bbotk]         FALSE        FALSE        FALSE       FALSE TRUE FALSE     FALSE FALSE 
    ## INFO  [16:00:39.805] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:39.805] [bbotk]  FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE  FALSE   FALSE 
    ## INFO  [16:00:39.805] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1   TTR TUBA3E 
    ## INFO  [16:00:39.805] [bbotk]    FALSE   FALSE FALSE   TRUE TRUE   FALSE  FALSE     FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:39.805] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5  ZIC2  ZIC5 
    ## INFO  [16:00:39.805] [bbotk]   FALSE FALSE   FALSE TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:39.805] [bbotk]                                               features classif.bacc 
    ## INFO  [16:00:39.805] [bbotk]  AC003092.1,AC116614.1,AHSG,APCDD1L_AS1,CHAT,CILP2,...    0.7701389 
    ## INFO  [16:00:39.861] [mlr3]  Applying learner 'classif.rpart.fselector' on task 'kirc_cla' (iter 1/3) 
    ## INFO  [16:00:40.170] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:40.172] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:40.814] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:40.823] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:40.925] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.028] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.113] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.219] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.330] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.453] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.574] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.758] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.853] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:41.964] [mlr3]  Finished benchmark 
    ## INFO  [16:00:42.371] [bbotk] Result of batch 1: 
    ## INFO  [16:00:42.392] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:42.392] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]        TRUE      FALSE      FALSE       TRUE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]        TRUE      FALSE      FALSE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]       FALSE      FALSE       TRUE      FALSE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:42.392] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE    FALSE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE    FALSE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE    FALSE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE       FALSE FALSE  TRUE    FALSE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE        TRUE FALSE FALSE    FALSE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:42.392] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE   TRUE FALSE    TRUE FALSE FALSE  TRUE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:42.392] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE   TRUE FALSE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]    FALSE   TRUE FALSE        FALSE   TRUE   TRUE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE  FALSE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:42.392] [bbotk]     FALSE  TRUE  TRUE     TRUE   FALSE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE  FALSE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]     FALSE FALSE  TRUE    FALSE    TRUE  FALSE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:42.392] [bbotk]     FALSE FALSE FALSE    FALSE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]     FALSE  TRUE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:42.392] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:42.392] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:42.392] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:42.392] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE     TRUE FALSE FALSE  TRUE   FALSE FALSE FALSE  TRUE FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE    TRUE FALSE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:42.392] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE      TRUE     FALSE     FALSE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:42.392] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:42.392] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE   TRUE  TRUE FALSE FALSE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE   TRUE FALSE FALSE FALSE   TRUE   TRUE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE   TRUE  TRUE  TRUE FALSE  FALSE   TRUE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE   TRUE  TRUE FALSE FALSE   TRUE  FALSE   TRUE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE FALSE FALSE   TRUE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE FALSE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:42.392] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         FALSE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:42.392] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE        TRUE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         TRUE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:42.392] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:42.392] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE         TRUE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:42.392] [bbotk]          FALSE        FALSE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:42.392] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:42.392] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:42.392] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE   FALSE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE FALSE FALSE FALSE    TRUE    TRUE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:42.392] [bbotk]      FALSE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:42.392] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:42.392] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:42.392] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE    TRUE    TRUE   FALSE FALSE  FALSE  TRUE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:42.392] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:42.392] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE FALSE    0.6640827 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE    TRUE FALSE  TRUE  TRUE FALSE    0.6453488 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6317829 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE FALSE FALSE    0.7635659 
    ## INFO  [16:00:42.392] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE   FALSE  TRUE  TRUE FALSE FALSE    0.7222222 
    ## INFO  [16:00:42.392] [bbotk]  FALSE  TRUE  FALSE   TRUE  TRUE   FALSE  TRUE FALSE FALSE FALSE    0.6989664 
    ## INFO  [16:00:42.392] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE    0.6873385 
    ## INFO  [16:00:42.392] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:42.392] [bbotk]                                 uhash 
    ## INFO  [16:00:42.392] [bbotk]  3cf01a7f-4215-4035-8014-6cfcc62ae36f 
    ## INFO  [16:00:42.392] [bbotk]  1979a832-d397-49da-91b8-d0730897ab40 
    ## INFO  [16:00:42.392] [bbotk]  3ebcfd81-b5a9-493f-84f3-7addd7f99c1e 
    ## INFO  [16:00:42.392] [bbotk]  a3f0f808-f4bc-4db8-ab94-f368b4b5c950 
    ## INFO  [16:00:42.392] [bbotk]  c29aa28c-17cd-4708-b5aa-41487ce162a6 
    ## INFO  [16:00:42.392] [bbotk]  9d374596-d2a8-4299-b17a-b03ad156c369 
    ## INFO  [16:00:42.392] [bbotk]  7dfbe5d8-9999-4a20-ac73-8e5366f405e2 
    ## INFO  [16:00:42.392] [bbotk]  7c3f38c0-2881-4290-98e4-0713c16edfc6 
    ## INFO  [16:00:42.392] [bbotk]  40cb727e-8fba-4b5e-9d24-8a7ddb884c3c 
    ## INFO  [16:00:42.392] [bbotk]  223dc748-dec9-45d1-a690-16855d46616e 
    ## INFO  [16:00:42.395] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:43.030] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:43.040] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.158] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.319] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.466] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.545] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.656] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.788] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:43.881] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:44.007] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:44.175] [mlr3]  Applying learner 'select.classif.rpart' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:44.286] [mlr3]  Finished benchmark 
    ## INFO  [16:00:44.690] [bbotk] Result of batch 2: 
    ## INFO  [16:00:44.712] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:44.712] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]        TRUE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]       FALSE       TRUE       TRUE      FALSE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]       FALSE      FALSE       TRUE      FALSE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]       FALSE      FALSE       TRUE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]        TRUE      FALSE       TRUE      FALSE       TRUE       TRUE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]        TRUE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:44.712] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE    FALSE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE       FALSE FALSE  TRUE    FALSE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE        TRUE FALSE FALSE     TRUE     TRUE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE       FALSE  TRUE  TRUE    FALSE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE        TRUE  TRUE FALSE     TRUE    FALSE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:44.712] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:44.712] [bbotk]   TRUE   TRUE FALSE   FALSE  TRUE  TRUE  TRUE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE         FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE   TRUE FALSE   FALSE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  FALSE  TRUE    TRUE FALSE FALSE FALSE         FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE   TRUE FALSE    TRUE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:44.712] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]     TRUE  FALSE  TRUE        FALSE  FALSE  FALSE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]     TRUE  FALSE  TRUE         TRUE  FALSE  FALSE FALSE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE  FALSE FALSE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE  FALSE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE  FALSE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]    FALSE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE   FALSE  TRUE FALSE  TRUE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:44.712] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:44.712] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]     FALSE FALSE FALSE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:44.712] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:44.712] [bbotk]     FALSE  TRUE FALSE     TRUE   FALSE   TRUE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:44.712] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:44.712] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE    TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE   FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE FALSE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE    TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     FALSE      TRUE      TRUE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE      TRUE      TRUE      TRUE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE     FALSE      TRUE      TRUE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:44.712] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:44.712] [bbotk]    TRUE  FALSE  TRUE  TRUE FALSE   TRUE  FALSE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE  FALSE  FALSE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE   TRUE FALSE FALSE  TRUE  FALSE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE   TRUE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE   TRUE FALSE  TRUE FALSE   TRUE  FALSE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  TRUE FALSE FALSE   TRUE FALSE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE  TRUE FALSE FALSE  TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE  TRUE FALSE FALSE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  TRUE FALSE FALSE   TRUE FALSE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]   TRUE FALSE  TRUE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         FALSE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:44.712] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE        FALSE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:44.712] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE        FALSE        FALSE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:44.712] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:44.712] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:44.712] [bbotk]       TRUE FALSE FALSE FALSE    TRUE   FALSE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:44.712] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:44.712] [bbotk]       TRUE FALSE  TRUE  TRUE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:44.712] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:44.712] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:44.712] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:44.712] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:44.712] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:44.712] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:44.712] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE   FALSE    TRUE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:44.712] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:44.712] [bbotk]   FALSE    TRUE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE    TRUE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE    TRUE    TRUE    TRUE FALSE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:44.712] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE    TRUE   FALSE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:44.712] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:44.712] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE   TRUE  FALSE  TRUE    TRUE FALSE FALSE  TRUE  TRUE    0.7054264 
    ## INFO  [16:00:44.712] [bbotk]   TRUE FALSE  FALSE  FALSE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE    0.7403101 
    ## INFO  [16:00:44.712] [bbotk]  FALSE  TRUE  FALSE   TRUE FALSE   FALSE  TRUE  TRUE  TRUE  TRUE    0.7564599 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE   TRUE   TRUE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE    0.6518088 
    ## INFO  [16:00:44.712] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6408269 
    ## INFO  [16:00:44.712] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE FALSE  TRUE  TRUE  TRUE    0.7635659 
    ## INFO  [16:00:44.712] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7009044 
    ## INFO  [16:00:44.712] [bbotk]                                 uhash 
    ## INFO  [16:00:44.712] [bbotk]  a21a25bd-eb5e-47d8-ae0a-b94125f9c8a2 
    ## INFO  [16:00:44.712] [bbotk]  d548c5a7-f4b5-4e9a-b5d0-4f98ee19b7ac 
    ## INFO  [16:00:44.712] [bbotk]  f11d0033-edda-424b-912f-1ff3f4750207 
    ## INFO  [16:00:44.712] [bbotk]  17adc298-5096-4b2a-aacc-24797df15016 
    ## INFO  [16:00:44.712] [bbotk]  cc1d59b6-704b-4eeb-8b5f-407e2794970b 
    ## INFO  [16:00:44.712] [bbotk]  69551cf8-56ed-4768-98f6-2dd62a3e6794 
    ## INFO  [16:00:44.712] [bbotk]  e47110ce-9697-4a19-b295-7e5f4d56cb31 
    ## INFO  [16:00:44.712] [bbotk]  7a68ab6a-72eb-44bd-928d-42ae5dfb00d2 
    ## INFO  [16:00:44.712] [bbotk]  05e39896-d948-4bd2-adbd-a376d416a3db 
    ## INFO  [16:00:44.712] [bbotk]  f36b4a9f-b77c-4c2b-ac0f-49b3c6a5e914 
    ## INFO  [16:00:44.726] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:44.727] [bbotk] Result: 
    ## INFO  [16:00:44.740] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM AHSG 
    ## INFO  [16:00:44.740] [bbotk]        TRUE      FALSE       TRUE       TRUE      FALSE       TRUE TRUE TRUE 
    ## INFO  [16:00:44.740] [bbotk]   AMH APCDD1L_AS1 AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:44.740] [bbotk]  TRUE        TRUE TRUE TRUE    FALSE     TRUE  TRUE TRUE     TRUE      TRUE 
    ## INFO  [16:00:44.740] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:44.740] [bbotk]  TRUE   TRUE  TRUE    TRUE TRUE  TRUE  TRUE          TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:44.740] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:44.740] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE TRUE  TRUE  FALSE 
    ## INFO  [16:00:44.740] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:44.740] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE TRUE   TRUE    FALSE      TRUE   TRUE     TRUE 
    ## INFO  [16:00:44.740] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:44.740] [bbotk]  FALSE FALSE     TRUE    TRUE  FALSE    TRUE FALSE     TRUE    FALSE  TRUE 
    ## INFO  [16:00:44.740] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17 KLK1 KLK15  KNG1 KRT7 L1CAM LECT1 
    ## INFO  [16:00:44.740] [bbotk]     FALSE  TRUE  TRUE  TRUE    TRUE FALSE TRUE  TRUE FALSE TRUE  TRUE FALSE 
    ## INFO  [16:00:44.740] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:44.740] [bbotk]       TRUE      TRUE     FALSE      TRUE      TRUE      TRUE  FALSE   TRUE 
    ## INFO  [16:00:44.740] [bbotk]  MAGEC3 MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP  PI3 
    ## INFO  [16:00:44.740] [bbotk]    TRUE TRUE FALSE FALSE  FALSE   TRUE   TRUE  TRUE   TRUE TRUE FALSE TRUE TRUE 
    ## INFO  [16:00:44.740] [bbotk]  PITX1  PLG PRR15L PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:44.740] [bbotk]   TRUE TRUE   TRUE TRUE  TRUE  TRUE TRUE        FALSE          TRUE 
    ## INFO  [16:00:44.740] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:44.740] [bbotk]          FALSE         TRUE         FALSE          TRUE          TRUE 
    ## INFO  [16:00:44.740] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:44.740] [bbotk]         FALSE        TRUE          TRUE         TRUE         TRUE          TRUE 
    ## INFO  [16:00:44.740] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6 RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [16:00:44.740] [bbotk]          TRUE         TRUE         TRUE        TRUE TRUE TRUE     FALSE TRUE 
    ## INFO  [16:00:44.740] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:44.740] [bbotk]   TRUE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   TRUE   FALSE 
    ## INFO  [16:00:44.740] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1   TTR TUBA3E 
    ## INFO  [16:00:44.740] [bbotk]     TRUE    TRUE  TRUE  FALSE TRUE    TRUE   TRUE      TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:44.740] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [16:00:44.740] [bbotk]   FALSE  TRUE    TRUE TRUE  TRUE TRUE TRUE 
    ## INFO  [16:00:44.740] [bbotk]                                                  features classif.bacc 
    ## INFO  [16:00:44.740] [bbotk]  AC003092.1,AC006262.5,AC007879.6,AC116614.1,AFM,AHSG,...    0.7635659 
    ## INFO  [16:00:44.805] [mlr3]  Applying learner 'classif.ranger.fselector' on task 'kirc_cla' (iter 1/3) 
    ## INFO  [16:00:45.030] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:45.033] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:46.162] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:46.171] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:46.312] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:46.447] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:46.599] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:46.763] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:46.917] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:47.071] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:47.216] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:47.348] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:47.479] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:47.634] [mlr3]  Finished benchmark 
    ## INFO  [16:00:48.188] [bbotk] Result of batch 1: 
    ## INFO  [16:00:48.202] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:48.202] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]       FALSE       TRUE       TRUE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]        TRUE       TRUE       TRUE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]        TRUE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]        TRUE      FALSE      FALSE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:48.202] [bbotk]   TRUE        TRUE FALSE FALSE    FALSE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE       FALSE  TRUE  TRUE     TRUE     TRUE FALSE FALSE     TRUE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE        TRUE FALSE  TRUE     TRUE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE        TRUE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE       FALSE  TRUE FALSE     TRUE     TRUE  TRUE FALSE    FALSE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE  TRUE         FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE   TRUE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE  TRUE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:48.202] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]    FALSE   TRUE  TRUE         TRUE  FALSE   TRUE FALSE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE  FALSE  TRUE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE  FALSE  TRUE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE FALSE   FALSE FALSE FALSE  TRUE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:48.202] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:48.202] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:48.202] [bbotk]      TRUE FALSE  TRUE    FALSE   FALSE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:48.202] [bbotk]      TRUE  TRUE FALSE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:48.202] [bbotk]      TRUE FALSE FALSE    FALSE    TRUE  FALSE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [16:00:48.202] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:48.202] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:48.202] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:48.202] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE     TRUE 
    ## INFO  [16:00:48.202] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:48.202] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:48.202] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE    FALSE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE     TRUE FALSE  TRUE FALSE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE    TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:48.202] [bbotk]  FALSE      TRUE      TRUE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE     FALSE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE      TRUE     FALSE     FALSE      TRUE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:48.202] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:48.202] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE  FALSE  FALSE  FALSE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]    TRUE  FALSE  TRUE FALSE  TRUE  FALSE   TRUE  FALSE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   TRUE  TRUE FALSE  TRUE   TRUE  FALSE  FALSE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   TRUE  TRUE FALSE  TRUE  FALSE  FALSE  FALSE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]    TRUE   TRUE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:48.202] [bbotk]   TRUE FALSE FALSE FALSE  FALSE  TRUE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE  TRUE FALSE  FALSE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  TRUE FALSE  TRUE   TRUE FALSE FALSE FALSE  TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE FALSE  TRUE   TRUE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:48.202] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE         FALSE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:48.202] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:48.202] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE        FALSE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:48.202] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE         TRUE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE        FALSE         TRUE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE        FALSE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE        FALSE         TRUE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:48.202] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:48.202] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:48.202] [bbotk]       TRUE  TRUE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:48.202] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:48.202] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE   FALSE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:48.202] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:48.202] [bbotk]       TRUE  TRUE  TRUE FALSE   FALSE    TRUE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:48.202] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:48.202] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:48.202] [bbotk]       TRUE  TRUE FALSE  TRUE    TRUE   FALSE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:48.202] [bbotk]       TRUE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:48.202] [bbotk]      FALSE  TRUE FALSE FALSE    TRUE   FALSE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:48.202] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:48.202] [bbotk]    TRUE   FALSE    TRUE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE    TRUE    TRUE   FALSE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   FALSE   FALSE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE    TRUE    TRUE    TRUE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:48.202] [bbotk]   FALSE    TRUE   FALSE   FALSE FALSE  FALSE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:48.202] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:48.202] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.6152882 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6271930 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE  FALSE  FALSE FALSE    TRUE  TRUE FALSE  TRUE  TRUE    0.5889724 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6271930 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE FALSE  TRUE FALSE    0.5983709 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6008772 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE  FALSE   TRUE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.5933584 
    ## INFO  [16:00:48.202] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6008772 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.5983709 
    ## INFO  [16:00:48.202] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE   FALSE  TRUE FALSE FALSE  TRUE    0.6510025 
    ## INFO  [16:00:48.202] [bbotk]                                 uhash 
    ## INFO  [16:00:48.202] [bbotk]  8a98cd1e-936b-48be-ad2b-f296dd2abb28 
    ## INFO  [16:00:48.202] [bbotk]  fed0e012-322f-4a67-a1c9-ab4d820abdab 
    ## INFO  [16:00:48.202] [bbotk]  ca45e595-9732-40dc-9701-2b53bfc0f0a7 
    ## INFO  [16:00:48.202] [bbotk]  d3c898df-b974-4404-ba59-4c1035f655bc 
    ## INFO  [16:00:48.202] [bbotk]  41f2b18a-ab8c-4b92-99f8-639f65ec7c20 
    ## INFO  [16:00:48.202] [bbotk]  4a910597-7e67-49cb-b83d-591b9a708df3 
    ## INFO  [16:00:48.202] [bbotk]  215196e0-515f-4338-8b98-8f942ab66bf6 
    ## INFO  [16:00:48.202] [bbotk]  bbc569b7-0397-47ac-a5b5-8764a2e56de4 
    ## INFO  [16:00:48.202] [bbotk]  704c1d99-1adc-4486-a47f-93497afa7001 
    ## INFO  [16:00:48.202] [bbotk]  1a87c7ad-a564-49dd-8cc0-b34824875b7d 
    ## INFO  [16:00:48.205] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:49.195] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:49.203] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:49.350] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:49.483] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:49.625] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:49.751] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:49.894] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:50.045] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:50.196] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:50.350] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:50.487] [mlr3]  Applying learner 'select.classif.ranger' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:50.634] [mlr3]  Finished benchmark 
    ## INFO  [16:00:51.109] [bbotk] Result of batch 2: 
    ## INFO  [16:00:51.120] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:51.120] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]        TRUE      FALSE       TRUE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:51.120] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE    FALSE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE        TRUE  TRUE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE        TRUE  TRUE FALSE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE       FALSE FALSE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE        TRUE FALSE  TRUE    FALSE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE   TRUE FALSE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  FALSE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE   TRUE  TRUE   FALSE  TRUE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE   TRUE  TRUE   FALSE FALSE  TRUE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:51.120] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]     TRUE  FALSE FALSE         TRUE  FALSE   TRUE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]     TRUE  FALSE FALSE         TRUE  FALSE  FALSE FALSE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]    FALSE  FALSE  TRUE        FALSE   TRUE  FALSE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]    FALSE  FALSE  TRUE         TRUE   TRUE  FALSE FALSE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]     TRUE   TRUE FALSE         TRUE   TRUE   TRUE FALSE  TRUE FALSE  TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE   TRUE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  TRUE    TRUE FALSE FALSE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE   FALSE  TRUE FALSE FALSE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:51.120] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:51.120] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:51.120] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:51.120] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:51.120] [bbotk]      TRUE FALSE FALSE     TRUE   FALSE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:51.120] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:51.120] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE  FALSE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:51.120] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:51.120] [bbotk]     FALSE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:51.120] [bbotk]     FALSE FALSE FALSE     TRUE   FALSE  FALSE   FALSE FALSE     TRUE     TRUE 
    ## INFO  [16:00:51.120] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     TRUE FALSE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE    FALSE  TRUE  TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     TRUE  TRUE FALSE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE      TRUE     FALSE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE      TRUE      TRUE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE     FALSE     FALSE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     FALSE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:51.120] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:51.120] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE  FALSE FALSE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE   TRUE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE   TRUE FALSE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE  FALSE  TRUE FALSE FALSE   TRUE   TRUE   TRUE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE  FALSE  FALSE  TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE  FALSE  TRUE FALSE FALSE   TRUE   TRUE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  TRUE FALSE  TRUE  FALSE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  TRUE  TRUE FALSE   TRUE  TRUE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  TRUE  TRUE FALSE   TRUE  TRUE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE  TRUE  TRUE FALSE  FALSE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE FALSE FALSE   TRUE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE FALSE  TRUE  TRUE FALSE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE          TRUE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:51.120] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:51.120] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE       FALSE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         TRUE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:51.120] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:51.120] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE         TRUE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE        FALSE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE        FALSE         TRUE        FALSE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE FALSE 
    ## INFO  [16:00:51.120] [bbotk]          FALSE         TRUE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:51.120] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:51.120] [bbotk]      FALSE  TRUE FALSE  TRUE   FALSE    TRUE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:51.120] [bbotk]       TRUE  TRUE FALSE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]       TRUE  TRUE  TRUE FALSE   FALSE    TRUE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]      FALSE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]      FALSE  TRUE FALSE  TRUE    TRUE   FALSE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE    TRUE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:51.120] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:51.120] [bbotk]       TRUE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:51.120] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:51.120] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   FALSE    TRUE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE    TRUE   FALSE   FALSE FALSE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   FALSE   FALSE   FALSE  TRUE  FALSE  TRUE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE   FALSE   FALSE    TRUE  TRUE  FALSE FALSE    TRUE   TRUE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE   TRUE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:51.120] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE   TRUE FALSE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:51.120] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE  TRUE  TRUE FALSE    0.5745614 
    ## INFO  [16:00:51.120] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.5889724 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE   FALSE FALSE FALSE  TRUE  TRUE    0.6535088 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE FALSE  TRUE    0.6152882 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE   FALSE  TRUE  TRUE FALSE  TRUE    0.6127820 
    ## INFO  [16:00:51.120] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6390977 
    ## INFO  [16:00:51.120] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE FALSE    0.6654135 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE  FALSE   TRUE  TRUE   FALSE FALSE FALSE  TRUE FALSE    0.6271930 
    ## INFO  [16:00:51.120] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.6390977 
    ## INFO  [16:00:51.120] [bbotk]  FALSE FALSE   TRUE   TRUE FALSE   FALSE  TRUE  TRUE  TRUE FALSE    0.6365915 
    ## INFO  [16:00:51.120] [bbotk]                                 uhash 
    ## INFO  [16:00:51.120] [bbotk]  94650b7e-b334-4fbd-b81b-1f56f87f2b10 
    ## INFO  [16:00:51.120] [bbotk]  bf9de174-df18-4fc0-8630-3cc008213c01 
    ## INFO  [16:00:51.120] [bbotk]  8022efa4-d295-468c-9fb3-38a9a66313c2 
    ## INFO  [16:00:51.120] [bbotk]  cc88e2c3-ffde-46fa-8408-d9c9af7a5bbc 
    ## INFO  [16:00:51.120] [bbotk]  0e07fdd0-0e1a-417a-8b77-1ff260e13c29 
    ## INFO  [16:00:51.120] [bbotk]  ba40f114-65df-40e1-874c-260e22e638ef 
    ## INFO  [16:00:51.120] [bbotk]  f97c13ce-0c4f-4c69-b97b-4004e90a8f99 
    ## INFO  [16:00:51.120] [bbotk]  ebbca09e-e553-41f4-8da9-c55e43d14950 
    ## INFO  [16:00:51.120] [bbotk]  42d1cb38-29bc-4ac1-995b-ae1fbfc7eabe 
    ## INFO  [16:00:51.120] [bbotk]  fabb2fee-191f-42f6-8907-c82cd97c4ee0 
    ## INFO  [16:00:51.128] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:51.128] [bbotk] Result: 
    ## INFO  [16:00:51.134] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM AHSG 
    ## INFO  [16:00:51.134] [bbotk]        TRUE      FALSE       TRUE      FALSE      FALSE       TRUE FALSE TRUE 
    ## INFO  [16:00:51.134] [bbotk]    AMH APCDD1L_AS1  AQP2 AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:51.134] [bbotk]  FALSE       FALSE FALSE TRUE     TRUE    FALSE  TRUE TRUE     TRUE      TRUE 
    ## INFO  [16:00:51.134] [bbotk]   CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:51.134] [bbotk]  TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:51.134] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:51.134] [bbotk]    FALSE  FALSE  TRUE         TRUE   TRUE  FALSE FALSE FALSE TRUE FALSE   TRUE 
    ## INFO  [16:00:51.134] [bbotk]  FDCSP FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:51.134] [bbotk]   TRUE TRUE   FALSE  TRUE  TRUE TRUE   TRUE    FALSE      TRUE   TRUE    FALSE 
    ## INFO  [16:00:51.134] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:51.134] [bbotk]   TRUE  TRUE     TRUE    TRUE  FALSE   FALSE FALSE     TRUE     TRUE FALSE 
    ## INFO  [16:00:51.134] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [16:00:51.134] [bbotk]      TRUE  TRUE FALSE  TRUE   FALSE  TRUE FALSE  TRUE TRUE TRUE FALSE FALSE 
    ## INFO  [16:00:51.134] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:51.134] [bbotk]       TRUE     FALSE      TRUE     FALSE     FALSE      TRUE   TRUE   TRUE 
    ## INFO  [16:00:51.134] [bbotk]  MAGEC3 MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3  PAEP  PI3 
    ## INFO  [16:00:51.134] [bbotk]    TRUE TRUE TRUE TRUE   TRUE  FALSE  FALSE  TRUE  FALSE TRUE  TRUE FALSE TRUE 
    ## INFO  [16:00:51.134] [bbotk]  PITX1   PLG PRR15L PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:51.134] [bbotk]   TRUE FALSE  FALSE TRUE  TRUE  TRUE TRUE        FALSE          TRUE 
    ## INFO  [16:00:51.134] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:51.134] [bbotk]           TRUE        FALSE         FALSE          TRUE         FALSE 
    ## INFO  [16:00:51.134] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:51.134] [bbotk]          TRUE        TRUE         FALSE         TRUE        FALSE         FALSE 
    ## INFO  [16:00:51.134] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [16:00:51.134] [bbotk]          TRUE        FALSE        FALSE       FALSE FALSE TRUE     FALSE TRUE 
    ## INFO  [16:00:51.134] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:51.134] [bbotk]  FALSE  TRUE    TRUE   FALSE    TRUE   FALSE    TRUE    TRUE   TRUE   FALSE 
    ## INFO  [16:00:51.134] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1   TTR TUBA3E 
    ## INFO  [16:00:51.134] [bbotk]    FALSE   FALSE  TRUE  FALSE TRUE    TRUE  FALSE      TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:51.134] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5  ZIC2  ZIC5 
    ## INFO  [16:00:51.134] [bbotk]    TRUE  TRUE    TRUE TRUE  TRUE FALSE FALSE 
    ## INFO  [16:00:51.134] [bbotk]                                                 features classif.bacc 
    ## INFO  [16:00:51.134] [bbotk]  AC003092.1,AC006262.5,AC116614.1,AHSG,AQP6,ATP6V0A4,...    0.6654135 
    ## INFO  [16:00:51.282] [mlr3]  Applying learner 'classif.ksvm.fselector' on task 'kirc_cla' (iter 3/3) 
    ## INFO  [16:00:51.518] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:51.520] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:52.261] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:52.274] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:52.385] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:52.519] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:52.664] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:52.785] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:52.972] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:53.096] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:53.176] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:53.262] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:53.390] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:53.518] [mlr3]  Finished benchmark 
    ## INFO  [16:00:53.927] [bbotk] Result of batch 1: 
    ## INFO  [16:00:53.938] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:53.938] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE       TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]        TRUE       TRUE      FALSE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE      FALSE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]        TRUE       TRUE      FALSE       TRUE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE      FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]       FALSE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]       FALSE       TRUE      FALSE      FALSE      FALSE       TRUE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:53.938] [bbotk]  FALSE        TRUE FALSE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE       FALSE  TRUE  TRUE    FALSE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE       FALSE FALSE  TRUE    FALSE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE    FALSE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE       FALSE  TRUE FALSE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE       FALSE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE        TRUE  TRUE  TRUE    FALSE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE   TRUE  TRUE    TRUE FALSE FALSE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:53.938] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE  TRUE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE   TRUE FALSE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE   TRUE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE  FALSE FALSE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE  FALSE FALSE        FALSE   TRUE  FALSE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE FALSE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE   FALSE FALSE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE FALSE    TRUE  TRUE FALSE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE FALSE   TRUE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:53.938] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:53.938] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:00:53.938] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [16:00:53.938] [bbotk]      TRUE  TRUE  TRUE     TRUE   FALSE  FALSE    TRUE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:53.938] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:00:53.938] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:53.938] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE    TRUE FALSE     TRUE    FALSE 
    ## INFO  [16:00:53.938] [bbotk]     FALSE FALSE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:53.938] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:53.938] [bbotk]      TRUE  TRUE  TRUE    FALSE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:53.938] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:53.938] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE   FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     TRUE  TRUE FALSE  TRUE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE    FALSE FALSE FALSE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE    FALSE FALSE  TRUE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:53.938] [bbotk]   TRUE      TRUE      TRUE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE      TRUE     FALSE     FALSE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     FALSE     FALSE     FALSE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE     FALSE      TRUE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     FALSE      TRUE     FALSE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE     FALSE     FALSE     FALSE     FALSE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:53.938] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:53.938] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE  FALSE  FALSE   TRUE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE  FALSE FALSE  TRUE  TRUE   TRUE  FALSE  FALSE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE  FALSE   TRUE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE  FALSE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE  FALSE   TRUE  TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE  FALSE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE FALSE  TRUE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE FALSE FALSE  TRUE   TRUE  TRUE FALSE  TRUE FALSE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE  TRUE FALSE   TRUE  TRUE FALSE  TRUE FALSE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE  TRUE  TRUE  FALSE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:53.938] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         TRUE        TRUE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         TRUE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         TRUE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:00:53.938] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:53.938] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE        FALSE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         TRUE        FALSE         TRUE       FALSE FALSE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:53.938] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:53.938] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:53.938] [bbotk]      FALSE  TRUE FALSE  TRUE    TRUE   FALSE    TRUE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:53.938] [bbotk]       TRUE  TRUE FALSE FALSE   FALSE    TRUE   FALSE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:53.938] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE   FALSE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:53.938] [bbotk]      FALSE FALSE FALSE  TRUE   FALSE   FALSE   FALSE    TRUE   FALSE   FALSE 
    ## INFO  [16:00:53.938] [bbotk]      FALSE FALSE  TRUE  TRUE   FALSE    TRUE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:53.938] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:53.938] [bbotk]       TRUE  TRUE  TRUE  TRUE   FALSE    TRUE    TRUE    TRUE    TRUE   FALSE 
    ## INFO  [16:00:53.938] [bbotk]      FALSE FALSE FALSE  TRUE    TRUE   FALSE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:53.938] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:53.938] [bbotk]       TRUE FALSE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:53.938] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:53.938] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE  FALSE  TRUE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE    TRUE   FALSE    TRUE FALSE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE    TRUE   FALSE    TRUE FALSE   TRUE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE    TRUE   FALSE   FALSE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE    TRUE   FALSE   FALSE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:53.938] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:53.938] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE    TRUE FALSE FALSE FALSE FALSE    0.7664141 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.6931818 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.7828283 
    ## INFO  [16:00:53.938] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7272727 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE   TRUE   TRUE FALSE    TRUE  TRUE  TRUE FALSE FALSE    0.7272727 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE  FALSE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7941919 
    ## INFO  [16:00:53.938] [bbotk]  FALSE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE FALSE  TRUE    0.7436869 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7550505 
    ## INFO  [16:00:53.938] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7386364 
    ## INFO  [16:00:53.938] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7941919 
    ## INFO  [16:00:53.938] [bbotk]                                 uhash 
    ## INFO  [16:00:53.938] [bbotk]  01aa268f-bb88-41fa-9a27-6249dca03953 
    ## INFO  [16:00:53.938] [bbotk]  6ee24602-12ac-4c0b-8460-2bfa5e8131f1 
    ## INFO  [16:00:53.938] [bbotk]  09bbff3b-7f1d-4160-97a9-986a7108b890 
    ## INFO  [16:00:53.938] [bbotk]  1407c9af-6baf-4702-abf4-3cf418f7f459 
    ## INFO  [16:00:53.938] [bbotk]  af07e0ed-a3ae-4d5f-95f3-7ff63619f85b 
    ## INFO  [16:00:53.938] [bbotk]  57706429-7b57-440d-8214-e9b00805ee7c 
    ## INFO  [16:00:53.938] [bbotk]  08bcefbd-2acc-4e4e-ab4a-20e832ffce1d 
    ## INFO  [16:00:53.938] [bbotk]  b6cf1ea7-a72d-456f-a4f3-f13d303f92ef 
    ## INFO  [16:00:53.938] [bbotk]  461b8516-509f-410a-b2e2-137185e8d8fd 
    ## INFO  [16:00:53.938] [bbotk]  70a52ddf-ff3a-4a4e-a140-62e6b04ee575 
    ## INFO  [16:00:53.942] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:54.642] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:54.651] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:54.754] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:54.943] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.081] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.182] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.284] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.420] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.602] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.750] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.881] [mlr3]  Applying learner 'select.classif.ksvm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:55.982] [mlr3]  Finished benchmark 
    ## INFO  [16:00:56.385] [bbotk] Result of batch 2: 
    ## INFO  [16:00:56.395] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:56.395] [bbotk]       FALSE      FALSE      FALSE      FALSE       TRUE      FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:56.395] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]        TRUE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]       FALSE       TRUE       TRUE      FALSE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:56.395] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE     TRUE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE       FALSE FALSE  TRUE     TRUE    FALSE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE  TRUE FALSE    FALSE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE       FALSE FALSE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:56.395] [bbotk]  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  FALSE  TRUE   FALSE  TRUE FALSE FALSE          TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE   TRUE FALSE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE   TRUE  TRUE    TRUE  TRUE FALSE FALSE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:56.395] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]    FALSE  FALSE FALSE        FALSE   TRUE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE FALSE FALSE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE   TRUE FALSE        FALSE   TRUE  FALSE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE  FALSE FALSE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE    TRUE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE    TRUE  TRUE  TRUE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE   FALSE  TRUE FALSE FALSE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE FALSE    TRUE FALSE  TRUE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:56.395] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE    TRUE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]     FALSE  TRUE  TRUE    FALSE    TRUE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]     FALSE FALSE  TRUE    FALSE   FALSE  FALSE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:00:56.395] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:56.395] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:00:56.395] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE   FALSE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:56.395] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:56.395] [bbotk]   TRUE    FALSE FALSE FALSE FALSE    TRUE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE   FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     TRUE FALSE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE   FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:56.395] [bbotk]  FALSE     FALSE      TRUE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE     FALSE      TRUE      TRUE     FALSE      TRUE     FALSE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:00:56.395] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:56.395] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE   TRUE  FALSE FALSE  FALSE  TRUE FALSE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE  FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE  FALSE FALSE  TRUE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   FALSE   TRUE  TRUE FALSE  TRUE   TRUE   TRUE  FALSE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:56.395] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE FALSE  TRUE  TRUE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE FALSE  TRUE  TRUE  FALSE  TRUE FALSE FALSE  TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE FALSE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE FALSE  TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:56.395] [bbotk]          FALSE         FALSE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         FALSE        FALSE          TRUE          TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:56.395] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:56.395] [bbotk]          FALSE          TRUE        FALSE          TRUE         FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:56.395] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:56.395] [bbotk]          FALSE        FALSE       FALSE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE        TRUE         FALSE        FALSE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]          FALSE        FALSE       FALSE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:00:56.395] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:56.395] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE  TRUE FALSE 
    ## INFO  [16:00:56.395] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE FALSE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE         TRUE         TRUE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:56.395] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:56.395] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:56.395] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]      FALSE  TRUE  TRUE  TRUE   FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]      FALSE FALSE FALSE FALSE    TRUE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:56.395] [bbotk]       TRUE  TRUE FALSE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]      FALSE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE   FALSE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]       TRUE  TRUE FALSE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:00:56.395] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:56.395] [bbotk]   FALSE    TRUE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE    TRUE    TRUE   FALSE  TRUE   TRUE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE    TRUE   FALSE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE   FALSE    TRUE   FALSE  TRUE  FALSE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:56.395] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.7500000 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE  TRUE FALSE FALSE    0.7272727 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7828283 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7272727 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7828283 
    ## INFO  [16:00:56.395] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.7386364 
    ## INFO  [16:00:56.395] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE FALSE    0.7828283 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE FALSE  TRUE  TRUE    0.7828283 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.7828283 
    ## INFO  [16:00:56.395] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE    0.7828283 
    ## INFO  [16:00:56.395] [bbotk]                                 uhash 
    ## INFO  [16:00:56.395] [bbotk]  f3d7e784-4ef3-4e71-aeba-7ec33b8f20fa 
    ## INFO  [16:00:56.395] [bbotk]  c326aa1f-6970-43d2-8c9d-25710ac125b2 
    ## INFO  [16:00:56.395] [bbotk]  1b57df12-60b3-4460-9cd6-930a097ed396 
    ## INFO  [16:00:56.395] [bbotk]  517b6608-f8bd-42cc-adb2-3f88ebe7a329 
    ## INFO  [16:00:56.395] [bbotk]  24052558-9dda-4be0-8363-4d7580583974 
    ## INFO  [16:00:56.395] [bbotk]  cadf49fb-29b0-427d-97af-d7c1b8d458d3 
    ## INFO  [16:00:56.395] [bbotk]  4eb69909-033a-4be7-86fa-fc601495f3fa 
    ## INFO  [16:00:56.395] [bbotk]  b2e6b61f-36f1-4542-9f08-0938966c4c42 
    ## INFO  [16:00:56.395] [bbotk]  b14c3735-e7b2-4254-8094-0160f1d45ac1 
    ## INFO  [16:00:56.395] [bbotk]  e539c01a-eb1b-4785-b2b7-74c91232fd9f 
    ## INFO  [16:00:56.403] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:00:56.403] [bbotk] Result: 
    ## INFO  [16:00:56.410] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1  AFM  AHSG 
    ## INFO  [16:00:56.410] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE      FALSE TRUE FALSE 
    ## INFO  [16:00:56.410] [bbotk]   AMH APCDD1L_AS1 AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1 BSND C10orf99 C14orf180 
    ## INFO  [16:00:56.410] [bbotk]  TRUE       FALSE TRUE FALSE     TRUE     TRUE  TRUE TRUE     TRUE      TRUE 
    ## INFO  [16:00:56.410] [bbotk]   CA1 CASP14 CCNA1 CDC42P2 CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8 CLMP 
    ## INFO  [16:00:56.410] [bbotk]  TRUE   TRUE  TRUE    TRUE TRUE FALSE  TRUE          TRUE  FALSE  TRUE TRUE 
    ## INFO  [16:00:56.410] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1  EN2 ESRP1 FAM83B 
    ## INFO  [16:00:56.410] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE FALSE TRUE  TRUE   TRUE 
    ## INFO  [16:00:56.410] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4 GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 HEPACAM2 
    ## INFO  [16:00:56.410] [bbotk]  FALSE FALSE    TRUE  TRUE FALSE TRUE   TRUE     TRUE     FALSE  FALSE     TRUE 
    ## INFO  [16:00:56.410] [bbotk]  HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 IGLC7 
    ## INFO  [16:00:56.410] [bbotk]   TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE  TRUE 
    ## INFO  [16:00:56.410] [bbotk]  IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15 KNG1 KRT7 L1CAM LECT1 
    ## INFO  [16:00:56.410] [bbotk]      TRUE  TRUE FALSE  TRUE    TRUE FALSE FALSE  TRUE TRUE TRUE  TRUE FALSE 
    ## INFO  [16:00:56.410] [bbotk]  LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 MAGEC2 
    ## INFO  [16:00:56.410] [bbotk]      FALSE      TRUE      TRUE      TRUE     FALSE     FALSE   TRUE   TRUE 
    ## INFO  [16:00:56.410] [bbotk]  MAGEC3 MFI2 MYH8 NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L OTX1 PADI3 PAEP  PI3 
    ## INFO  [16:00:56.410] [bbotk]   FALSE TRUE TRUE TRUE   TRUE  FALSE   TRUE  TRUE  FALSE TRUE FALSE TRUE TRUE 
    ## INFO  [16:00:56.410] [bbotk]  PITX1  PLG PRR15L PSG9 PVALB RAB25 RHBG RP11_10O22.1 RP11_150O12.1 
    ## INFO  [16:00:56.410] [bbotk]   TRUE TRUE   TRUE TRUE FALSE  TRUE TRUE         TRUE         FALSE 
    ## INFO  [16:00:56.410] [bbotk]  RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 RP11_425D17.2 
    ## INFO  [16:00:56.410] [bbotk]           TRUE         TRUE          TRUE          TRUE         FALSE 
    ## INFO  [16:00:56.410] [bbotk]  RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 RP11_690G19.4 
    ## INFO  [16:00:56.410] [bbotk]          TRUE        TRUE          TRUE         TRUE         TRUE          TRUE 
    ## INFO  [16:00:56.410] [bbotk]  RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1 SAA1 SAA2_SAA4 SBSN 
    ## INFO  [16:00:56.410] [bbotk]          TRUE        FALSE         TRUE       FALSE FALSE TRUE      TRUE TRUE 
    ## INFO  [16:00:56.410] [bbotk]  SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 SLC4A1 SLC6A15 
    ## INFO  [16:00:56.410] [bbotk]   TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE   FALSE   TRUE   FALSE 
    ## INFO  [16:00:56.410] [bbotk]  SLC6A18 SLC7A13 STAC2 TCEAL2 TCL6 TMEM213 TMEM61 TMPRSS11E TNNT1  TTR TUBA3E 
    ## INFO  [16:00:56.410] [bbotk]     TRUE    TRUE  TRUE  FALSE TRUE    TRUE   TRUE      TRUE  TRUE TRUE  FALSE 
    ## INFO  [16:00:56.410] [bbotk]  TUBBP6 TUNAR UGT1A10 UMOD WFDC5 ZIC2 ZIC5 
    ## INFO  [16:00:56.410] [bbotk]    TRUE FALSE    TRUE TRUE  TRUE TRUE TRUE 
    ## INFO  [16:00:56.410] [bbotk]                                                 features classif.bacc 
    ## INFO  [16:00:56.410] [bbotk]  AC003092.1,AC006262.5,AC007879.6,AC104654.2,AFM,AMH,...    0.7941919 
    ## INFO  [16:00:56.473] [mlr3]  Applying learner 'classif.svm.fselector' on task 'kirc_cla' (iter 3/3) 
    ## INFO  [16:00:56.748] [bbotk] Starting to optimize 150 parameter(s) with '<FSelectorRandomSearch>' and '<TerminatorEvals> [n_evals=20]' 
    ## INFO  [16:00:56.751] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:57.451] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:57.459] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.530] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.590] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.652] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.722] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.818] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.907] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:57.981] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:58.043] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:58.111] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:58.185] [mlr3]  Finished benchmark 
    ## INFO  [16:00:58.563] [bbotk] Result of batch 1: 
    ## INFO  [16:00:58.574] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:00:58.574] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]        TRUE      FALSE       TRUE      FALSE      FALSE      FALSE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]       FALSE       TRUE      FALSE      FALSE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]       FALSE       TRUE      FALSE      FALSE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]        TRUE       TRUE       TRUE       TRUE      FALSE       TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]        TRUE       TRUE      FALSE       TRUE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:00:58.574] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE       FALSE FALSE  TRUE    FALSE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE       FALSE FALSE FALSE    FALSE    FALSE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE       FALSE FALSE FALSE     TRUE     TRUE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE        TRUE FALSE  TRUE     TRUE     TRUE FALSE  TRUE    FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE        TRUE  TRUE FALSE    FALSE     TRUE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE        TRUE FALSE FALSE     TRUE     TRUE FALSE  TRUE     TRUE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE  FALSE FALSE   FALSE  TRUE  TRUE  TRUE          TRUE  FALSE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE  FALSE FALSE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE   TRUE FALSE   FALSE FALSE  TRUE FALSE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE  FALSE  TRUE   FALSE FALSE  TRUE  TRUE         FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  FALSE FALSE    TRUE FALSE FALSE FALSE         FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:00:58.574] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]    FALSE   TRUE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]    FALSE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE  TRUE FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]    FALSE  FALSE  TRUE        FALSE  FALSE  FALSE FALSE FALSE  TRUE FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]    FALSE  FALSE  TRUE         TRUE  FALSE   TRUE FALSE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]     TRUE   TRUE FALSE        FALSE  FALSE  FALSE FALSE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]     TRUE   TRUE  TRUE         TRUE  FALSE  FALSE  TRUE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]    FALSE   TRUE FALSE         TRUE   TRUE  FALSE  TRUE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]     TRUE  FALSE FALSE         TRUE   TRUE  FALSE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE   FALSE FALSE FALSE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  TRUE   FALSE  TRUE  TRUE  TRUE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE   FALSE FALSE  TRUE  TRUE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  TRUE    TRUE  TRUE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE   FALSE FALSE FALSE FALSE   TRUE     TRUE     FALSE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE    TRUE  TRUE FALSE FALSE   TRUE    FALSE      TRUE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:00:58.574] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]     FALSE FALSE FALSE     TRUE    TRUE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]      TRUE FALSE FALSE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]     FALSE  TRUE FALSE    FALSE   FALSE   TRUE   FALSE  TRUE    FALSE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]     FALSE  TRUE FALSE     TRUE   FALSE  FALSE    TRUE  TRUE     TRUE    FALSE 
    ## INFO  [16:00:58.574] [bbotk]      TRUE FALSE FALSE    FALSE   FALSE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:00:58.574] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:00:58.574] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE    FALSE  TRUE FALSE  TRUE    TRUE  TRUE FALSE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE    FALSE FALSE  TRUE FALSE   FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE    FALSE  TRUE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE    FALSE FALSE FALSE  TRUE   FALSE FALSE  TRUE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE     TRUE  TRUE  TRUE FALSE   FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE    FALSE  TRUE  TRUE  TRUE   FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE    TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE     TRUE  TRUE  TRUE FALSE    TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:00:58.574] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE      TRUE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE      TRUE     FALSE     FALSE     FALSE     FALSE     FALSE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE     FALSE      TRUE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE      TRUE     FALSE      TRUE     FALSE      TRUE     FALSE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE     FALSE     FALSE     FALSE      TRUE      TRUE      TRUE  FALSE 
    ## INFO  [16:00:58.574] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:00:58.574] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   TRUE FALSE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE  FALSE FALSE FALSE  TRUE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]    TRUE  FALSE FALSE  TRUE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]    TRUE   TRUE FALSE FALSE  TRUE  FALSE   TRUE   TRUE  TRUE  FALSE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE  FALSE FALSE  TRUE FALSE   TRUE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   TRUE FALSE FALSE  TRUE   TRUE   TRUE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE  TRUE FALSE FALSE  FALSE  TRUE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE  TRUE  TRUE  FALSE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE  TRUE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE FALSE FALSE   TRUE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  TRUE FALSE FALSE  FALSE  TRUE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE FALSE  TRUE  FALSE FALSE FALSE  TRUE  TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  TRUE  TRUE  TRUE   TRUE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE          TRUE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:00:58.574] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:00:58.574] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         TRUE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         TRUE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:00:58.574] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:00:58.574] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE        FALSE        FALSE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE         TRUE        FALSE       FALSE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE        FALSE        FALSE         TRUE       FALSE FALSE FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE         TRUE        FALSE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]           TRUE        FALSE        FALSE        FALSE        TRUE FALSE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:00:58.574] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:00:58.574] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:58.574] [bbotk]       TRUE  TRUE  TRUE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE  TRUE FALSE FALSE   FALSE    TRUE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE  TRUE  TRUE FALSE   FALSE   FALSE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE FALSE  TRUE  TRUE    TRUE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE    TRUE   FALSE 
    ## INFO  [16:00:58.574] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE   FALSE   FALSE   FALSE   FALSE    TRUE 
    ## INFO  [16:00:58.574] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE FALSE    TRUE   TRUE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE    TRUE   FALSE  TRUE  FALSE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE   TRUE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE    TRUE    TRUE FALSE  FALSE FALSE   FALSE  FALSE      TRUE 
    ## INFO  [16:00:58.574] [bbotk]   FALSE   FALSE   FALSE   FALSE  TRUE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:00:58.574] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6903305 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE  FALSE   TRUE FALSE    TRUE FALSE FALSE FALSE FALSE    0.6028152 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE FALSE    0.5734394 
    ## INFO  [16:00:58.574] [bbotk]   TRUE  TRUE   TRUE  FALSE  TRUE    TRUE  TRUE FALSE FALSE FALSE    0.6670747 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6028152 
    ## INFO  [16:00:58.574] [bbotk]  FALSE  TRUE  FALSE  FALSE FALSE    TRUE FALSE FALSE  TRUE FALSE    0.6376989 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE  FALSE  FALSE  TRUE    TRUE FALSE FALSE  TRUE FALSE    0.6523868 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE   TRUE   TRUE  TRUE   FALSE FALSE  TRUE FALSE FALSE    0.6028152 
    ## INFO  [16:00:58.574] [bbotk]   TRUE FALSE  FALSE   TRUE  TRUE   FALSE FALSE FALSE  TRUE FALSE    0.6376989 
    ## INFO  [16:00:58.574] [bbotk]  FALSE FALSE   TRUE  FALSE FALSE    TRUE FALSE  TRUE  TRUE FALSE    0.6260710 
    ## INFO  [16:00:58.574] [bbotk]                                 uhash 
    ## INFO  [16:00:58.574] [bbotk]  3803c438-7005-443c-8f3c-0540e82d92d6 
    ## INFO  [16:00:58.574] [bbotk]  3ccbaeae-b181-4c19-9ddc-6578bbf39567 
    ## INFO  [16:00:58.574] [bbotk]  7cfaf6e5-d71e-4fd8-baba-a3e9f1ceaa23 
    ## INFO  [16:00:58.574] [bbotk]  21a57db1-b9da-49c2-8142-50e018181d75 
    ## INFO  [16:00:58.574] [bbotk]  2178b5f4-709f-4e3b-9bc1-5bbd495ee9fb 
    ## INFO  [16:00:58.574] [bbotk]  90fc6a22-215f-440b-a40d-b6d6e8bfc5a5 
    ## INFO  [16:00:58.574] [bbotk]  ba327e0a-4e09-4ec9-887e-9d89d2fb9be2 
    ## INFO  [16:00:58.574] [bbotk]  2cbefef0-78d5-4c5c-916f-6673ba42f714 
    ## INFO  [16:00:58.574] [bbotk]  bb2fa684-9597-4153-ae2e-e2dba2a78951 
    ## INFO  [16:00:58.574] [bbotk]  c19193b0-e1ec-48c2-aee7-0ec198055c8f 
    ## INFO  [16:00:58.577] [bbotk] Evaluating 10 configuration(s) 
    ## INFO  [16:00:59.369] [mlr3]  Running benchmark with 10 resampling iterations 
    ## INFO  [16:00:59.376] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.438] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.516] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.615] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.717] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.811] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.876] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:00:59.936] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:01:00.006] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:01:00.074] [mlr3]  Applying learner 'select.classif.svm' on task 'kirc_cla' (iter 1/1) 
    ## INFO  [16:01:00.169] [mlr3]  Finished benchmark 
    ## INFO  [16:01:00.552] [bbotk] Result of batch 2: 
    ## INFO  [16:01:00.563] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:01:00.563] [bbotk]       FALSE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]       FALSE       TRUE      FALSE      FALSE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]        TRUE       TRUE      FALSE       TRUE       TRUE       TRUE  TRUE FALSE 
    ## INFO  [16:01:00.563] [bbotk]       FALSE      FALSE      FALSE       TRUE      FALSE      FALSE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]       FALSE       TRUE      FALSE      FALSE       TRUE       TRUE FALSE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]        TRUE       TRUE       TRUE       TRUE       TRUE       TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]        TRUE      FALSE       TRUE       TRUE       TRUE      FALSE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:01:00.563] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE FALSE FALSE     TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE        TRUE FALSE  TRUE    FALSE    FALSE  TRUE  TRUE    FALSE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE     TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE FALSE    FALSE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE        TRUE FALSE  TRUE     TRUE     TRUE FALSE  TRUE     TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE       FALSE  TRUE FALSE    FALSE    FALSE FALSE FALSE    FALSE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE        TRUE  TRUE  TRUE     TRUE     TRUE  TRUE  TRUE     TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:01:00.563] [bbotk]   TRUE   TRUE FALSE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  FALSE FALSE    TRUE  TRUE FALSE FALSE          TRUE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE FALSE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE         FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE   TRUE  TRUE   FALSE FALSE  TRUE FALSE          TRUE   TRUE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  FALSE  TRUE    TRUE  TRUE  TRUE  TRUE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE   TRUE  TRUE   FALSE  TRUE  TRUE FALSE          TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  FALSE  TRUE   FALSE  TRUE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE   TRUE FALSE    TRUE FALSE  TRUE  TRUE          TRUE  FALSE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:01:00.563] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]    FALSE  FALSE  TRUE         TRUE   TRUE   TRUE  TRUE FALSE FALSE  TRUE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]    FALSE  FALSE FALSE         TRUE   TRUE   TRUE  TRUE FALSE FALSE  TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]     TRUE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE  TRUE  TRUE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]    FALSE   TRUE  TRUE        FALSE   TRUE   TRUE  TRUE  TRUE FALSE FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]     TRUE   TRUE  TRUE         TRUE   TRUE   TRUE  TRUE  TRUE  TRUE  TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE  TRUE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]     TRUE   TRUE  TRUE        FALSE   TRUE  FALSE  TRUE  TRUE  TRUE FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE    TRUE  TRUE FALSE FALSE  FALSE     TRUE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  TRUE   FALSE FALSE FALSE  TRUE  FALSE     TRUE      TRUE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE   FALSE  TRUE  TRUE  TRUE   TRUE    FALSE     FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  TRUE   FALSE FALSE  TRUE  TRUE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  TRUE   FALSE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE   TRUE     TRUE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  TRUE    TRUE FALSE  TRUE  TRUE  FALSE    FALSE     FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:01:00.563] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE   TRUE    TRUE FALSE    FALSE     TRUE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE  FALSE    TRUE FALSE    FALSE    FALSE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE FALSE  TRUE     TRUE    TRUE  FALSE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE FALSE FALSE    FALSE    TRUE   TRUE   FALSE FALSE     TRUE    FALSE 
    ## INFO  [16:01:00.563] [bbotk]     FALSE  TRUE FALSE    FALSE    TRUE  FALSE   FALSE  TRUE    FALSE     TRUE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE  TRUE FALSE    FALSE   FALSE   TRUE    TRUE FALSE     TRUE     TRUE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE  TRUE FALSE     TRUE    TRUE   TRUE    TRUE  TRUE     TRUE     TRUE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE FALSE FALSE    FALSE   FALSE   TRUE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:01:00.563] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:01:00.563] [bbotk]      TRUE  TRUE  TRUE     TRUE    TRUE  FALSE   FALSE  TRUE     TRUE     TRUE 
    ## INFO  [16:01:00.563] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:01:00.563] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE    FALSE  TRUE  TRUE  TRUE   FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE    FALSE  TRUE  TRUE FALSE   FALSE  TRUE FALSE  TRUE FALSE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE     TRUE  TRUE FALSE  TRUE    TRUE  TRUE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE     TRUE  TRUE FALSE FALSE   FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE     TRUE  TRUE  TRUE  TRUE    TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE    FALSE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE     TRUE FALSE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE  TRUE FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE    FALSE  TRUE  TRUE FALSE    TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:01:00.563] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE     FALSE     FALSE      TRUE      TRUE     FALSE     FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE      TRUE     FALSE      TRUE      TRUE      TRUE     FALSE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE     FALSE     FALSE      TRUE      TRUE     FALSE      TRUE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE      TRUE      TRUE     FALSE      TRUE     FALSE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE     FALSE     FALSE      TRUE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE     FALSE      TRUE      TRUE      TRUE      TRUE      TRUE   TRUE 
    ## INFO  [16:01:00.563] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:01:00.563] [bbotk]    TRUE   TRUE  TRUE FALSE FALSE   TRUE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE  FALSE  TRUE  TRUE FALSE   TRUE  FALSE   TRUE FALSE   TRUE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]    TRUE   TRUE  TRUE  TRUE  TRUE  FALSE   TRUE   TRUE  TRUE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE   TRUE  FALSE   TRUE FALSE   TRUE FALSE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]    TRUE   TRUE FALSE FALSE FALSE   TRUE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE  FALSE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   TRUE  TRUE  TRUE  TRUE   TRUE   TRUE   TRUE FALSE   TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE   TRUE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE  FALSE  TRUE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE  FALSE  TRUE  TRUE  TRUE   TRUE  FALSE   TRUE  TRUE   TRUE FALSE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE  TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE  TRUE FALSE  FALSE  TRUE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE FALSE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE FALSE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE FALSE  TRUE FALSE   TRUE FALSE FALSE FALSE  TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE  TRUE FALSE   TRUE  TRUE  TRUE  TRUE FALSE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE  TRUE  TRUE   TRUE  TRUE  TRUE  TRUE  TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]   TRUE FALSE FALSE  TRUE  FALSE  TRUE  TRUE FALSE  TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:01:00.563] [bbotk]           TRUE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         FALSE         TRUE          TRUE         FALSE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE          TRUE         TRUE          TRUE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE          TRUE        FALSE          TRUE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE          TRUE         TRUE         FALSE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         FALSE         TRUE          TRUE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         FALSE        FALSE         FALSE          TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE         FALSE         TRUE         FALSE         FALSE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE          TRUE         TRUE          TRUE         FALSE 
    ## INFO  [16:01:00.563] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:01:00.563] [bbotk]           TRUE        FALSE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE        FALSE        TRUE          TRUE         TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE        TRUE          TRUE         TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE        TRUE         FALSE         TRUE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE        FALSE       FALSE         FALSE         TRUE        FALSE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE        TRUE          TRUE        FALSE         TRUE 
    ## INFO  [16:01:00.563] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:01:00.563] [bbotk]           TRUE        FALSE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE  TRUE FALSE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE         TRUE        FALSE       FALSE FALSE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE         TRUE        FALSE        FALSE        TRUE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE         TRUE         TRUE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]           TRUE         TRUE         TRUE         TRUE        TRUE  TRUE  TRUE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE        FALSE         TRUE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]          FALSE         TRUE        FALSE         TRUE        TRUE FALSE FALSE 
    ## INFO  [16:01:00.563] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:01:00.563] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]      FALSE  TRUE  TRUE FALSE    TRUE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]       TRUE  TRUE  TRUE  TRUE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]      FALSE FALSE  TRUE FALSE    TRUE   FALSE    TRUE   FALSE    TRUE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]       TRUE  TRUE  TRUE FALSE   FALSE   FALSE    TRUE    TRUE   FALSE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]       TRUE  TRUE  TRUE FALSE    TRUE    TRUE    TRUE    TRUE    TRUE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE   FALSE 
    ## INFO  [16:01:00.563] [bbotk]      FALSE FALSE  TRUE FALSE   FALSE   FALSE    TRUE   FALSE   FALSE    TRUE 
    ## INFO  [16:01:00.563] [bbotk]       TRUE FALSE  TRUE  TRUE    TRUE   FALSE   FALSE    TRUE    TRUE   FALSE 
    ## INFO  [16:01:00.563] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:01:00.563] [bbotk]    TRUE   FALSE    TRUE    TRUE  TRUE  FALSE  TRUE    TRUE   TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   FALSE    TRUE   FALSE  TRUE  FALSE  TRUE   FALSE  FALSE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   FALSE    TRUE    TRUE FALSE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   FALSE   FALSE    TRUE FALSE  FALSE FALSE   FALSE   TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   FALSE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]    TRUE    TRUE    TRUE    TRUE  TRUE   TRUE  TRUE   FALSE   TRUE      TRUE 
    ## INFO  [16:01:00.563] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]    TRUE   FALSE    TRUE   FALSE FALSE  FALSE FALSE    TRUE  FALSE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]    TRUE   FALSE   FALSE    TRUE FALSE  FALSE  TRUE   FALSE   TRUE     FALSE 
    ## INFO  [16:01:00.563] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 classif.bacc 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE   FALSE  TRUE  TRUE  TRUE FALSE    0.6291310 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  TRUE  FALSE  FALSE  TRUE   FALSE FALSE FALSE FALSE FALSE    0.6028152 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE  TRUE  TRUE    0.6028152 
    ## INFO  [16:01:00.563] [bbotk]   TRUE FALSE   TRUE  FALSE FALSE   FALSE FALSE  TRUE FALSE  TRUE    0.6144431 
    ## INFO  [16:01:00.563] [bbotk]  FALSE  TRUE   TRUE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE    0.6468788 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE  FALSE   TRUE  TRUE    TRUE  TRUE  TRUE  TRUE  TRUE    0.6670747 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE  TRUE  TRUE FALSE  TRUE    0.6291310 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE  TRUE FALSE FALSE FALSE    0.5881273 
    ## INFO  [16:01:00.563] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE  TRUE FALSE    0.6493268 
    ## INFO  [16:01:00.563] [bbotk]   TRUE  TRUE   TRUE   TRUE  TRUE    TRUE FALSE  TRUE FALSE  TRUE    0.6407589 
    ## INFO  [16:01:00.563] [bbotk]                                 uhash 
    ## INFO  [16:01:00.563] [bbotk]  971a9a93-b8af-49ee-bd54-f37ea7265bc9 
    ## INFO  [16:01:00.563] [bbotk]  0be8537a-ea3a-46bd-8bd3-260fc72b0fb5 
    ## INFO  [16:01:00.563] [bbotk]  67f103a9-2700-4003-9d42-f247c5df0eff 
    ## INFO  [16:01:00.563] [bbotk]  9a19f344-2651-4550-a07a-fea9ef91416e 
    ## INFO  [16:01:00.563] [bbotk]  669cb7a8-e370-467b-a1a8-277d0cbebae9 
    ## INFO  [16:01:00.563] [bbotk]  587c3a88-4469-43ed-be60-d81f715e9232 
    ## INFO  [16:01:00.563] [bbotk]  74efd6a2-bfa2-474c-b443-e848132474a7 
    ## INFO  [16:01:00.563] [bbotk]  d53c7650-f69f-497c-b3c3-294c5313a04a 
    ## INFO  [16:01:00.563] [bbotk]  d965ec82-1461-414c-a484-4e5a393bda7a 
    ## INFO  [16:01:00.563] [bbotk]  dbc99287-10fa-4068-a944-22b3cc6fb275 
    ## INFO  [16:01:00.571] [bbotk] Finished optimizing after 20 evaluation(s) 
    ## INFO  [16:01:00.572] [bbotk] Result: 
    ## INFO  [16:01:00.579] [bbotk]  AC003092.1 AC006262.4 AC006262.5 AC007879.6 AC104654.2 AC116614.1   AFM  AHSG 
    ## INFO  [16:01:00.579] [bbotk]       FALSE      FALSE      FALSE      FALSE      FALSE      FALSE FALSE FALSE 
    ## INFO  [16:01:00.579] [bbotk]    AMH APCDD1L_AS1  AQP2  AQP6 ATP6V0A4 ATP6V0D2 BARX1  BSND C10orf99 C14orf180 
    ## INFO  [16:01:00.579] [bbotk]  FALSE       FALSE FALSE FALSE    FALSE    FALSE FALSE FALSE    FALSE     FALSE 
    ## INFO  [16:01:00.579] [bbotk]    CA1 CASP14 CCNA1 CDC42P2  CHAT CIDEC CILP2 CITF22_24E5.1 CLCNKB CLDN8  CLMP 
    ## INFO  [16:01:00.579] [bbotk]  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE         FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.579] [bbotk]  COL11A1 COL7A1 CPNE7 CTD_2008P7.9 CXCL13 CYP1A1 DMRT2  DQX1   EN2 ESRP1 FAM83B 
    ## INFO  [16:01:00.579] [bbotk]    FALSE  FALSE FALSE        FALSE  FALSE  FALSE FALSE FALSE FALSE FALSE  FALSE 
    ## INFO  [16:01:00.579] [bbotk]  FDCSP  FGF5 FKBP9P1 FOXI2 FXYD4  GGT6 GLB1L3 GOLGA6L2 GOLGA6L7P GPR110 
    ## INFO  [16:01:00.579] [bbotk]  FALSE FALSE   FALSE FALSE FALSE FALSE  FALSE    FALSE      TRUE  FALSE 
    ## INFO  [16:01:00.579] [bbotk]  HEPACAM2 HHATL HMGA2 HS3ST3A1 IGF2BP3 IGFBP1 IGFL1P1 IGFN1 IGHV1_69 IGKV3_11 
    ## INFO  [16:01:00.579] [bbotk]     FALSE FALSE FALSE    FALSE   FALSE  FALSE   FALSE FALSE    FALSE    FALSE 
    ## INFO  [16:01:00.579] [bbotk]  IGLC7 IGLV3_19 INHBE ITPKA KCNJ1 KIRREL3 KLF17  KLK1 KLK15  KNG1  KRT7 L1CAM 
    ## INFO  [16:01:00.579] [bbotk]  FALSE    FALSE  TRUE FALSE FALSE   FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:01:00.579] [bbotk]  LECT1 LINC00890 LINC00942 LINC00973 LINC01187 LINC01436 LINC01559 LRRTM1 
    ## INFO  [16:01:00.579] [bbotk]  FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE  FALSE 
    ## INFO  [16:01:00.579] [bbotk]  MAGEC2 MAGEC3  MFI2  MYH8  NFE4 NIPAL4 NKX2_2 NKX2_3 NMRK2 NUPR1L  OTX1 PADI3 
    ## INFO  [16:01:00.579] [bbotk]   FALSE  FALSE FALSE FALSE FALSE  FALSE  FALSE  FALSE FALSE  FALSE FALSE FALSE 
    ## INFO  [16:01:00.579] [bbotk]   PAEP   PI3 PITX1   PLG PRR15L  PSG9 PVALB RAB25  RHBG RP11_10O22.1 
    ## INFO  [16:01:00.579] [bbotk]  FALSE FALSE FALSE FALSE  FALSE FALSE FALSE FALSE FALSE        FALSE 
    ## INFO  [16:01:00.579] [bbotk]  RP11_150O12.1 RP11_161D15.1 RP11_310H4.6 RP11_314M24.1 RP11_400N13.3 
    ## INFO  [16:01:00.579] [bbotk]          FALSE         FALSE        FALSE         FALSE         FALSE 
    ## INFO  [16:01:00.579] [bbotk]  RP11_425D17.2 RP11_440G9.1 RP11_54H7.4 RP11_554D15.1 RP11_586K2.1 RP11_643A5.3 
    ## INFO  [16:01:00.579] [bbotk]          FALSE        FALSE       FALSE         FALSE        FALSE        FALSE 
    ## INFO  [16:01:00.579] [bbotk]  RP11_690G19.4 RP11_95M15.2 RP13_895J2.6 RP4_568C11.4 RP5_984P4.6  RTL1  SAA1 
    ## INFO  [16:01:00.579] [bbotk]          FALSE        FALSE        FALSE        FALSE       FALSE FALSE FALSE 
    ## INFO  [16:01:00.579] [bbotk]  SAA2_SAA4  SBSN SFTPB SHOX2 SLC12A3 SLC18A3 SLC22A8 SLC30A8 SLC34A1 SLC38A5 
    ## INFO  [16:01:00.579] [bbotk]      FALSE FALSE FALSE FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
    ## INFO  [16:01:00.579] [bbotk]  SLC4A1 SLC6A15 SLC6A18 SLC7A13 STAC2 TCEAL2  TCL6 TMEM213 TMEM61 TMPRSS11E 
    ## INFO  [16:01:00.579] [bbotk]   FALSE   FALSE   FALSE   FALSE FALSE  FALSE FALSE   FALSE  FALSE     FALSE 
    ## INFO  [16:01:00.579] [bbotk]  TNNT1   TTR TUBA3E TUBBP6 TUNAR UGT1A10  UMOD WFDC5  ZIC2  ZIC5 
    ## INFO  [16:01:00.579] [bbotk]  FALSE FALSE  FALSE  FALSE FALSE   FALSE FALSE FALSE FALSE FALSE 
    ## INFO  [16:01:00.579] [bbotk]         features classif.bacc 
    ## INFO  [16:01:00.579] [bbotk]  GOLGA6L7P,INHBE    0.6903305 
    ## INFO  [16:01:00.762] [mlr3]  Finished benchmark

``` r
cols <- c("task_id", "learner_id",  "classif.bacc", "classif.sensitivity", "classif.specificity")

measures <- list(
  msr("classif.bacc"),
  msr("classif.sensitivity"), 
  msr("classif.specificity")
)

bmr_df <- bmr$aggregate(measures) %>%
  dplyr::select(cols) %>%
  dplyr::arrange(desc(classif.bacc))

bmr_df
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["task_id"],"name":[1],"type":["chr"],"align":["left"]},{"label":["learner_id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["classif.bacc"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["classif.sensitivity"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["classif.specificity"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"kirc_cla","2":"classif.ksvm.fselector","3":"0.7292395","4":"0.5034948","5":"0.9549842"},{"1":"kirc_cla","2":"classif.ranger.fselector","3":"0.7287103","4":"0.5025071","5":"0.9549134"},{"1":"kirc_cla","2":"classif.svm.fselector","3":"0.7067860","4":"0.4788034","5":"0.9347685"},{"1":"kirc_cla","2":"classif.xgboost.fselector","3":"0.6589721","4":"0.4783286","5":"0.8396156"},{"1":"kirc_cla","2":"classif.rpart.fselector","3":"0.6570917","4":"0.4392783","5":"0.8749052"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

# Plotting Benchmark Results

``` r
library("mlr3viz")
library("ggplot2")

autoplot(bmr, measure= msr("classif.bacc")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](figs/tutorial_unnamed-chunk-29-1.png)<!-- -->
