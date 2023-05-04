#Sourced from https://cran.r-project.org/web/packages/pscore/vignettes/metsss.html


## install the pscore package
# install.packages("pscore")

library(pscore)
library(dplyr)
library(ggplot2)
library("mice")


?MetSSS
mydata <- data.frame(
  sbp = c(122, 102.5),
  dbp = c(76.5, 64),
  trigs = c(1.47, 1.27),
  hdl = c(2.22, 1.59),
  waist = c(71, 91),
  glucose = c(5.16, 5.82),
  sex = c("Female", "Male"))
summary(mydata)
paste0("c('", paste(colnames(mydata), collapse = "', '"), "')")





# ---- manual_checking ----

head(d)

thr_wst_f <- 80
thr_wst_m <- 94
thr_tri <- 1.7
thr_hdl_f <- 1.3
thr_hdl_m <- 1
thr_sbp <- 130
thr_dbp <- 85
thr_bgl <- 5.6


# create copy
d_check <- d


# ---- step 1 ----


d_check <-
  d_check %>%
  mutate(
    waist =
      case_when(
        (sex == "Female") & (waist > thr_wst_f) ~ thr_wst_f,
        (sex ==   "Male") & (waist > thr_wst_m) ~ thr_wst_m,
        TRUE                                    ~ waist
      )
  )

inner_join(
  d %>% select(record_id, sex, waist),
  d_check %>% select(record_id, waist),
  "record_id"
) %>%
  with(., table(waist.x == waist.y, sex))


d_check <-
  d_check %>%
  mutate(
    hdl  =
      case_when(
        (sex == "Female") & (hdl  < thr_hdl_f) ~ thr_hdl_f,
        (sex ==   "Male") & (hdl  < thr_hdl_m) ~ thr_hdl_m,
        TRUE                                    ~ hdl 
      )
  )

inner_join(
  d %>% select(record_id, sex, hdl),
  d_check %>% select(record_id, hdl),
  "record_id"
) %>%
  with(., table(hdl.x == hdl.y, sex))




d_check <-
  d_check %>%
  mutate(
    trigs =
      case_when(
        (trigs > thr_tri) ~ thr_tri,
        TRUE              ~ trigs
      )
  )

inner_join(
  d %>% select(record_id, trigs),
  d_check %>% select(record_id, trigs),
  "record_id"
) %>%
  with(., table(trigs.x == trigs.y))

d_check <-
  d_check %>%
  mutate(
    sbp =
      case_when(
        (sbp > thr_sbp) ~ thr_sbp,
        TRUE              ~ sbp
      )
  )

inner_join(
  d %>% select(record_id, sbp),
  d_check %>% select(record_id, sbp),
  "record_id"
) %>%
  with(., table(sbp.x == sbp.y))


d_check <-
  d_check %>%
  mutate(
    dbp =
      case_when(
        (dbp > thr_dbp) ~ thr_dbp,
        TRUE              ~ dbp
      )
  )

inner_join(
  d %>% select(record_id, dbp),
  d_check %>% select(record_id, dbp),
  "record_id"
) %>%
  with(., table(dbp.x == dbp.y))



d_check <-
  d_check %>%
  mutate(
    glucose =
      case_when(
        (glucose > thr_bgl) ~ thr_bgl,
        TRUE              ~ glucose
      )
  )

inner_join(
  d %>% select(record_id, glucose),
  d_check %>% select(record_id, glucose),
  "record_id"
) %>%
  with(., table(glucose.x == glucose.y))

cbind(d, d_check)


# ---- step 2 ----

# ---- step 2 ----

d_check <-
  d_check %>%
  mutate(
    waist =   if_else(sex == "Female", waist - thr_wst_f, waist - thr_wst_m),
    hdl =     if_else(sex == "Female", thr_hdl_f - hdl, thr_hdl_m - hdl),
    trigs =   thr_tri - trigs,
    sbp =     thr_sbp - sbp,
    dbp =     thr_dbp - dbp,
    glucose = thr_bgl - glucose
  )


# ---- step 3 ----

# print fuinction
MetSSS
# print object used in function
MetSSS_model

str(MetSSS_model)

sig_std <- MetSSS_model@CompositeReady@sigma
# 13.2809156  3.1395191  0.2966108  0.1659097  8.5998314  0.1261516
sig_std[1]


d_check <-
  d_check %>%
  mutate(
    waist =   waist / sig_std[5],
    hdl =     hdl /  sig_std[4],
    trigs =   trigs / sig_std[3],
    sbp =     sbp / sig_std[1],
    dbp =     dbp / sig_std[2] ,
    glucose = glucose / sig_std[6]
  )



# ---- step 4 ----

(orthog_post <- MetSSS_model@pca$loadings[1:6, 1:6])
(col_order <- attr(MetSSS_model@pca$loadings, "dimnames")[[1]])
class(orthog_post)


orthog_d_check <- d_check[, col_order]
head(orthog_d_check)
orthog_d_check <- as.matrix(orthog_d_check)
class(orthog_d_check)

orthog_d_check <- orthog_d_check %*% orthog_post

# ---- step 5 ----


pca_comp_sds <- MetSSS_model@pca$sdev
pca_comp_sds_mat <- diag(pca_comp_sds)
inv_pca_comp_sds_mat <- pca_comp_sds_mat
diag(inv_pca_comp_sds_mat) <- 1 / diag(pca_comp_sds_mat)
# check :
solve(pca_comp_sds_mat) == inv_pca_comp_sds_mat

std_orthog_d_check <- orthog_d_check %*% inv_pca_comp_sds_mat
std_orthog_d_check

# ---- step 6 ----

sq_std_orthog_d_check <- std_orthog_d_check ^ 2

# ---- step 7 ----

sum_sq_std_orthog_d_check <- rowSums(sq_std_orthog_d_check)
length(sum_sq_std_orthog_d_check)

# ---- step 8 ----

metsss_manual_calc <- sqrt(sum_sq_std_orthog_d_check)
names(metsss_manual_calc) <- NULL

### Are these the same as what was calculated using the pscore package??
# (if not what is missing, what has been done wrong? is it the thresholds etc?)
metsss_manual_calc

hist(metsss_manual_calc)
hist(metsss)


write.csv(metsss_manual_calc, file = "metsss_manual_calc.csv", row.names = FALSE)



