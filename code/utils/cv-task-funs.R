run_cv_hgp <- function(params, .seed = 2024) {
  ki <- params$ki
  mod <- params$mod  
  ## Ensure libraries are loaded in the forked process
  library(rstan)
  ##--- model-specific data ----
  dists <- if(mod <= 2) dists_ct else dists_h
  ln_val <- as.integer(mod <= 6)
  homosc <- mod %% 2 != 0
  types_aux <- ifelse(W[, 2] == 1, "point", "poly")
  current_W <- if(homosc) W[, 1, drop = FALSE] else W
  ## For priors
  hd_lim <- range(dists[upper.tri(dists)])
  tau_0  <- sd(xy[, 1]) * .5
  nu <- .7
  ##--- k-fold split ----
  N <- NROW(dists)
  set.seed(.seed)
  ids_split <- sample(x = rep(seq_len(kf), ceiling(N / kf)), size = N)
  ids_pred <- which(ids_split == ki)
  ids_obs  <- which(ids_split != ki)
  N_pred <- length(ids_pred)
  N_obs  <- length(ids_obs)
  ##--- data partitions ----
  pred_xy <- xy[ids_pred, ]
  obs_xy <- xy[ids_obs, ]
  pred_W <- current_W[ids_pred, , drop = FALSE]
  obs_W <- current_W[ids_obs, , drop = FALSE]
  pred_dists <- dists[ids_pred, ids_pred]
  c_dists <- dists[ids_obs, ids_pred]
  obs_dists <- dists[ids_obs, ids_obs]
  types_aux <- types_aux[ids_pred]
  ##--- fit model ----
  hgp_dat <- list(
      N = N_obs,
      k = NCOL(obs_xy) - 1,
      k2 = ncol(obs_W),
      dists = obs_dists,
      rho_0 = hd_lim[2] * .8,
      tau_0 = tau_0,
      prob_rho = .05,
      pt_u = .05,
      y = obs_xy[, 1],
      X = obs_xy[, -1, drop = FALSE],
      A = obs_W,
      nu = nu,
      ln = ln_val
  )
  ## select the models (ugly code :()
  fit_model <- if(mod <= 4)
                 stan_models$hgp_gaus else
                                        stan_models$hgp_ln
  
  fit_hgp <- sampling(
      object = fit_model,
      data = hgp_dat,
      iter = samples + warmup,
      warmup = warmup,
      chains = chains,
      cores = 1,
      save_warmup = FALSE,
      seed = 2024
  )
  ##--- predictions ---
  pred_dat <- c(hgp_dat,
                list(Xp = pred_xy[, -1, drop = FALSE],
                     Ap = pred_W,
                     Np = N_pred,
                     distsp = pred_dists,
                     cdists = c_dists))
  
  pars_hgp <- c("beta", "alpha", "rho", "tau")
  if (mod > 4)
    pars_hgp <- c(pars_hgp, "z")
  pred_model <- if(mod <= 4)
                  stan_models$pred_gaus else
                                          stan_models$pred_ln
  pred_samples <- gqs(
      pred_model,
      data = pred_dat,
      draws = as.matrix(fit_hgp, pars = pars_hgp)
  )
  ##--- output ----
  pred_hgp <- as.matrix(pred_samples, pars = "y_rep")
  predictions <- data.frame(
      model_id = mod,
      k = ki,
      type = types_aux,
      obs = pred_xy[, 1],
      median_prd = apply(pred_hgp, 2, median),
      sd_prd = apply(pred_hgp, 2, sd),
      lower = apply(pred_hgp, 2, quantile, probs = .025),
      upper = apply(pred_hgp, 2, quantile, probs = .975)
  )
  rownames(predictions) <- NULL
  return(predictions)
}

run_cv_inla <- function(params) {
  ki <- params$ki
  mod <- params$mod
  library(INLA)
  library(sp)
  ##--- k-fold split ----
  N <- NROW(my_dt)
  set.seed(2024) ## Consistent seed for reproducible splits
  ids_split <- sample(x = rep(seq_len(kf), ceiling(N / kf)), size = N)
  ids_pred <- which(ids_split == ki)
  ids_obs  <- which(ids_split != ki)
  pred_dt <- my_dt[ids_pred, ]
  obs_dt  <- my_dt[ids_obs, ]
  ##--- data-proc for inla ----
  is_point_pred <- grepl("^st_", pred_dt$id)
  types_aux <- ifelse(is_point_pred, "point", "poly")
  ##--- * splitting stations and satellite ----
  stations <- obs_dt[grepl("^st_", obs_dt$id), ]
  dt_estim <- obs_dt[!grepl("^st_", obs_dt$id), ]
  coop <- sf::st_coordinates(stations)
  yp <- stations$pm25
  spol <- dt_estim |>
    sf::st_geometry() |>
    sf::as_Spatial()
  spol@proj4string <- CRS()
  ya <- dt_estim$pm25
  ##--- * test data has to be split into point and poly as well ----
  pred_dt_points <- pred_dt[is_point_pred, ]
  pred_dt_polys  <- pred_dt[!is_point_pred, ]
  coopred_points <- if (nrow(pred_dt_points) > 0)
                     sf::st_coordinates(pred_dt_points) else NULL
  coopred_polys  <- if (nrow(pred_dt_polys) > 0)
                      sf::st_coordinates(sf::st_centroid(pred_dt_polys)) else NULL
  coopred <- rbind(coopred_points, coopred_polys)
  pred_dt_ordered <- rbind(pred_dt_points, pred_dt_polys)
  ypred <- rep(NA_real_, nrow(coopred))
  ##--- creating mesh ----
  grid <- ifelse(mod == 9, "sparse", "fine")
  max_edge_ct <- ifelse(grid == "fine", 5000, 10000)
  max_edge <- c(max_edge_ct,
                10 * max_edge_ct)
  cutoff <- max_edge_ct / 100
  offset <- c(max_edge_ct,
              5 * max_edge_ct)
  my_mesh <- inla.mesh.2d(rbind(coopred, coop),
                          max.edge =
                            max_edge,
                          cutoff = cutoff,
                          offset = offset)
  my_spde <- inla.spde2.matern(mesh = my_mesh,
                               alpha = 2)
  ##--- data stack ----
  Ap <- inla.spde.make.A(mesh = my_mesh, loc = coop)
  stk_p <- inla.stack(tag = "point",
                      data = list(y = yp),
                      A = list(Ap, 1),
                      effects = list(s = 1:my_spde$n.spde,
                                     data.frame(b0 = rep(1, nrow(coop)))))
  ##--- * for numerical integration ----
  locin <- my_mesh$loc[as.vector(which(!is.na(over(SpatialPoints(my_mesh$loc),
                                                   spol)))), ]
  block <- rep(0, nrow(locin))
  for (i in seq_len(length(spol))) {
    block[as.vector(which(!is.na(over(SpatialPoints(locin), spol[i]))))] <- i
  }
  Aa <- inla.spde.make.A(mesh = my_mesh,
                         loc = locin,
                         block = block,
                         block.rescale = "sum")
  stk_a <- inla.stack(tag = "areal",
                      data = list(y = ya),
                      A = list(Aa, 1),
                      effects = list(s = 1:my_spde$n.spde,
                                     data.frame(b0 = rep(1, NROW(ya)))))
  ##--- * for predictions ----
  ##--- ** points ----
  ##--- *** there are significantly less points than polys ----
  if (NROW(pred_dt_points) > 0) {
    coopred_points <- sf::st_coordinates(pred_dt_points)
    ypred_p <- rep(NA_real_, nrow(pred_dt_points))
    Apred_p <- inla.spde.make.A(mesh = my_mesh, loc = coopred_points)
    stk_pred_p <- inla.stack(tag = "pred_p", data = list(y = ypred_p), A = list(Apred_p, 1),
                             effects = list(s = 1:my_spde$n.spde,
                                            data.frame(b0 = rep(1, NROW(ypred_p)))))
  }
  ##--- ** polygons ----
  sppr_polys <- sf::as_Spatial(sf::st_geometry(pred_dt_polys))
  sppr_polys@proj4string <- CRS() # Assign empty CRS
  ypred_a <- rep(NA_real_, nrow(pred_dt_polys))
  locin_pred <- my_mesh$loc[as.vector(which(!is.na(over(SpatialPoints(my_mesh$loc),
                                                        sppr_polys)))), ]
  block_pred <- rep(0, nrow(locin_pred))
  for (i in seq_len(length(sppr_polys))) {
    block_pred[as.vector(which(!is.na(over(SpatialPoints(locin_pred),
                                           sppr_polys[i]))))] <- i
  }
  Apred_a <- inla.spde.make.A(mesh = my_mesh,
                              loc = locin_pred,
                              block = block_pred,
                              block.rescale = "sum")
  stk_pred_a <- inla.stack(tag = "pred_a",
                           data = list(y = ypred_a),
                           A = list(Apred_a, 1),
                           effects = list(s = 1:my_spde$n.spde,
                                        data.frame(b0 = rep(1, NROW(ypred_a)))))
  ##--- ** combining data ----
  if (NROW(pred_dt_points) > 0) {
    full_dt <- inla.stack(stk_p, stk_a, stk_pred_p, stk_pred_a)
  } else {
    full_dt <- inla.stack(stk_p, stk_a, stk_pred_a)
  }
  ##--- fit ----
  formula <- y ~ -1 + b0 + f(s, model = my_spde)
  my_fit <- inla(formula,
                 data = inla.stack.data(full_dt),
                 control.predictor = list(compute = TRUE,
                                          A = inla.stack.A(full_dt)),
                 control.compute = list(return.marginals.predictor = TRUE),
                 verbose = FALSE)
  ##--- preds ids ----
  id_pred_p <- inla.stack.index(full_dt, "pred_p")$data
  id_pred_a <- inla.stack.index(full_dt, "pred_a")$data
  if (NROW(pred_dt_points) > 0) {
    id_pred <- c(id_pred_p, id_pred_a)
  } else {
    id_pred <- id_pred_a
  }
  ##--- output ----
  predictions <- data.frame(
      model_id = mod,
      k = ki,
      type = types_aux,
      obs = pred_dt_ordered$pm25,
      median_prd = my_fit$summary.fitted.values[id_pred, "0.5quant"],
      sd_prd = my_fit$summary.fitted.values$sd[id_pred],
      lower = my_fit$summary.fitted.values[id_pred, "0.025quant"],
      upper = my_fit$summary.fitted.values[id_pred, "0.975quant"]
  )
  rownames(predictions) <- NULL
  return(predictions)
}
