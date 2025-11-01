Loading required package: ggplot2
Loading required package: mlr3
Loading required package: mlr3learners
Loaded glmnet 4.1-8
Loading required package: iterators
Loading required package: parallel
Loading required package: future
randomForest 4.7-1.2
Type rfNews() to see new features/changes/bug fixes.

Attaching package: ‘randomForest’

The following object is masked from ‘package:ggplot2’:

    margin

── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.4     ✔ readr     2.1.5
✔ forcats   1.0.0     ✔ stringr   1.5.1
✔ lubridate 1.9.4     ✔ tibble    3.2.1
✔ purrr     1.0.4     ✔ tidyr     1.3.1
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ purrr::accumulate()    masks foreach::accumulate()
✖ dplyr::combine()       masks randomForest::combine()
✖ tidyr::expand()        masks Matrix::expand()
✖ dplyr::filter()        masks stats::filter()
✖ dplyr::lag()           masks stats::lag()
✖ randomForest::margin() masks ggplot2::margin()
✖ tidyr::pack()          masks Matrix::pack()
✖ tidyr::unpack()        masks Matrix::unpack()
✖ purrr::when()          masks foreach::when()
ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

Attaching package: ‘kernlab’

The following object is masked from ‘package:purrr’:

    cross

The following object is masked from ‘package:ggplot2’:

    alpha

Rows: 741 Columns: 1259
── Column specification ────────────────────────────────────────────────────────
Delimiter: ","
chr   (89): block, mainjob, i1, hhi1, ii4other, ii8other, ii10other, ii14oth...
dbl (1090): id, id2, district, camp, parish, village, household, person, no_...
lgl   (80): IDENTIFYING_VARIABLES, died, a_id, EMPLOYMENT_VARIABLES, EDUCATI...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
Warning message:
In summary.glm(object) :
  observations with zero weight not used for calculating dispersion
Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998991 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998938 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998983 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998958 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 482 first singular values only accounts for 0.998955 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998936 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998996 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998988 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998934 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, numdims = nc, printprogress = FALSE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998982 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  Truncated SVD with 481 first singular values only accounts for 0.998998 of the variance of "K". The biasbound optimization may not perform as expected. You many want to increase "maxnumdims" to capture more of the variance of "K" 

Warning in kbal::kbal(X, treatment = data$Z, printprogress = FALSE, mixed_data = TRUE,  :
  When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.

