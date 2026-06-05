

# Shift Southern hemisphere julian dates six months 
shift_s_date <- function(d){
  d <- d - 183
  if(d <= 0){
    d <- d + 365
  }
  return(d)
}

## To create confidence intervals ----
add_ci <- function(model, df){
  rhs <- formula(model, fixed.only = TRUE)[-2]
  X   <- model.matrix(rhs, df)
  beta <- fixef(model)
  V    <- vcov(model)
  
  fit <- X %*% beta
  se <- sqrt(rowSums((X %*% V) * X))
  
  df$lwr <- fit - 1.96 * se
  df$upr <- fit + 1.96 * se
  
  df
}

forward_selection <- function(start_model, candidates, aic_criterion, dat){
  
  best_model <- start_model
  best_aic <- AIC(best_model)
  
  result <- list(selected_vars = c(), 
                 aic_log = list(), 
                 r2_log = list()) 
  
  for(v in names(candidates)){
      
      test_model <- update(best_model,
                           paste(". ~ . + ", candidates[[v]]),
                           data = dat)
      test_aic <- AIC(test_model)
      
      # Keep variable if AIC improves 
      if((best_aic - test_aic) >= aic_criterion){
        
        best_aic <- test_aic
        best_model <- test_model
        
        result$selected_vars <- c(result$selected_vars, v)
        result$aic_log[[v]] <- best_aic
        
        if(class(best_model) != "gls")  result$r2_log[[v]] <- r.squaredGLMM(best_model)
        
      }
  cat("\r", which(names(candidates) == v), "/", length(candidates))    
  }
  cat("\n")
  result[["model"]] <- best_model
  return(result)
}

build_neighbour <- function(sf_m, dist, clustcol, clustval){
  
  coords <- st_coordinates(sf_m)
  ind <- sf_m[[clustcol]] == clustval
  coord_clust <- coords[ind, ]
  nb_clust <- dnearneigh(coord_clust, d1 = 0, d2 = dist) 
  lw_clust <- nb2listw(nb_clust, style = "W", zero.policy = TRUE)
  
  return(lw_clust)
}

forward_selection_gam <- function(start_model, candidates, aic_criterion, dat){
  
  best_model <- start_model
  best_aic <- AIC(best_model)
  best_formula <- formula(best_model)
  
  result <- list(selected_vars = c(), 
                 aic_log = list(), 
                 r2_log = list())
  
  result$aic_log$base <- AIC(best_model)
  result$r2_log$base <- summary(best_model)$r.sq 
  
  for(v in names(candidates)){
    
    test_formula <- update.formula(best_formula, paste(". ~ . +", candidates[v]))
    test_model <- mgcv::gam(
      test_formula,
      family = gaussian(),
      data = dat,
      method = "ML"
    )
    test_aic <- AIC(test_model)
    
    if(best_aic - test_aic >= aic_criterion){
      best_model <- test_model
      best_formula <- formula(best_model)
      best_aic <- test_aic
      result$aic_log[[candidates[[v]]]] <- test_aic
      result$r2_log[[candidates[[v]]]] <- summary(best_model)$r.sq
      result$selected_vars <- c(result$selected_vars, v)
    }
    cat("\r", which(names(candidates) == v), "/", length(candidates)) 
  }
  cat("\n")
  result[["model"]] <- best_model
  return(result)
}
