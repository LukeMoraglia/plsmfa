set.seed(42)
Z <- matrix(rnorm(100*7), 100, 7)


X.block1 <- matrix(Z[,1], ncol = 5, nrow = length(Z[,1]))
X.block1.withnoise <- X.block1 + matrix(rnorm(NROW(X.block1)*NCOL(X.block1), 0, 0.2), ncol = NCOL(X.block1), nrow = NROW(X.block1))


X.block2 <- cbind(Z[,2], Z[,3], Z[,4], Z[,5])


X <- cbind(X.block1.withnoise, X.block2)

rm(X.block1, X.block1.withnoise, X.block2)

Y.block1 <- matrix(Z[,1], ncol = 2, nrow = length(Z[,1]))
Y.block1.withnoise <- Y.block1 +
   matrix(rnorm(NROW(Y.block1)*NCOL(Y.block1), 0, 1),
          ncol = NCOL(Y.block1), nrow = NROW(Y.block1))

Y.block2 <- cbind(Z[,3],Z[,3],Z[,4],Z[,6])
Y.block2.withnoise <- Y.block2 +
   matrix(rnorm(NROW(Y.block2)*NCOL(Y.block2), 0, 0.6),
          ncol = NCOL(Y.block2), nrow = NROW(Y.block2))

Y <- cbind(Y.block1.withnoise, Y.block2.withnoise)

rm(Y.block1, Y.block1.withnoise, Y.block2, Y.block2.withnoise)
rm(Z)

colnames(X) <- paste0("X", 1:NCOL(X))
colnames(Y) <- paste0("Y", 1:NCOL(Y))

X.design <- matrix(c(rep("SubX1", times = 5), rep("SubX2", times = 4)), ncol = 1)
Y.design <- matrix(c(rep("SubY1", times = 2), rep("SubY2", times = 4)), ncol = 1)

qualitative <- c("#9932cc","#1b998b","#ff9b71","#544b3d","#f9dc5c","#4c1c00","#113969","#6969dd")

X.var.color <- list(
   oc = dplyr::recode(X.design, "SubX1" = qualitative[1], "SubX2" = qualitative[2]),
   gc = c("SubX1" = qualitative[1], "SubX2" = qualitative[2])
)

Y.var.color <- list(
   oc = dplyr::recode(Y.design, "SubY1" = qualitative[3], "SubY2" = qualitative[6]),
   gc = c("SubY1" = qualitative[3], "SubY2" = qualitative[6])
)


res_plsmfa_sim <- plsmfa(X, Y, X.design, Y.design)

kmeans_sim <- kmeans(cbind(res_plsmfa_sim$pls$lx[,1], res_plsmfa_sim$pls$ly[,1]), 3)
obs_design <- dplyr::recode(kmeans_sim$cluster, '1' = 'A', '2' = 'B', '3' = 'C')
col_obs_sim <- dplyr::recode(obs_design, 'A' = 'tan4', 'B' = 'darkgreen', 'C' = 'darkorchid')
names(col_obs_sim) <- obs_design
col_group_sim <- unique(col_obs_sim)
names(col_group_sim) <- unique(obs_design)

sim_data <- list(X = X,
                 Y = Y,
                 X_design = X.design,
                 Y_design = Y.design,
                 X_var_color = X.var.color,
                 Y_var_color = Y.var.color,
                 obs_design = obs_design,
                 obs_color = col_obs_sim,
                 group_color = col_group_sim)

usethis::use_data(sim_data, overwrite = TRUE)
