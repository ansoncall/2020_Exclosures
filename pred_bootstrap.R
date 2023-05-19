# Bootstrap and refit best predator models to assess stability of coefficients
# Must be sourced from data_analysis.R

# best models
best_mod_list %>% names

# define functions ####
# Subsample a dataframe, removing data from one field at a time
bootstrap_data <- function(site, field, data) {
  newdat <- data %>% filter(Site != site | Field != field)
  return(newdat)
}
# try to update a model with new data
# upon error, return "error"
rev_update <- function(new_data, mod) {
  message(mod$call)
  result <- tryCatch({
    update(mod, data = new_data)
  },
  warning = function(w) {
    message(mod$call)
    message(paste0(w, "\n"))
    return(update(mod, data = new_data))
  },
  error = function(e) {
    message(mod$call)
    message(paste0(e, "\n"))
    return(NA)
  },
  finally = {
    # cleanup-code
  })
  return(result)
}
bootstrapGLM <- function(mod, field_list, site_list, dFrame) {
  # determine whether mod contains veg survey data
  # if so, filter data exclude rows where veg data is missing
  if (sum(c("totalCover", "shan", "rich") %in% names(mod$frame)) > 0) {
    dFrame %<>% filter(!is.na(shan))
  }
  # build data frame of all combinations of field and site
  df <- expand.grid(field_list, site_list)
  # create list of names of dropped fields
  dfNames <- df %>% mutate(names = paste0(Var2, Var1)) %>% pull(names)
  # define data to be used in bootstrap_data
  bootstrap_data_def <- function(x, y){
    return(bootstrap_data(x, y, data = dFrame))
  }
  # create list of bootstrapped data frames
  df_list <- map2(df$Var2, df$Var1, bootstrap_data_def)
  # for each data frame, refit the model to the bootstrapped data
  list <- lapply(df_list, rev_update, mod)
  # name the list of results
  names(list) <- dfNames
  return(list)
}

# bootstrap and re-fit ####
# build list of names for bootstrapping
sites_sp <- c('Minden',
              'Lovelock',
              'Fallon',
              'Yerington')
sites_fa <- c('Lovelock',
              'Fallon',
              'Yerington')
fields <- c(1:3)
# separate spring and fall mods (bootstrapped data is different for each)
bsMods_sp <- lapply(best_mod_list[1:5],
                    bootstrapGLM,
                    field_list = fields,
                    site_list = sites_sp,
                    dFrame = df_sp)
bsMods_fa <- lapply(best_mod_list[6:10],
                    bootstrapGLM,
                    field_list = fields,
                    site_list = sites_fa,
                    dFrame = df_fa %>% mutate(Site = fct_drop(Site)))
# build names for tibble
spDF <- expand.grid(fields, sites_sp)
faDF <- expand.grid(fields, sites_fa)
dfNames_sp <- spDF %>% mutate(names = paste0(Var2, Var1)) %>% pull(names)
dfNames_fa <- faDF %>% mutate(names = paste0(Var2, Var1)) %>% pull(names)
# collect betas
# build empty tibble to hold stats
bs_stats_df_sp <- tibble(dropped = dfNames_sp,
                         MarginalR2 = c(0),
                         ConditionalR2 = c(0),
                         effects1 = c("none"),
                         effects2 = c("none"),
                         coefs1 = c(0),
                         coefs2 = c(0),
                         coefsmin1 = c(0),
                         coefsmax1 = c(0),
                         coefsmin2 = c(0),
                         coefsmax2 = c(0),
                         coefsse1 = c(0),
                         coefsse2 = c(0),
)
bs_stats_df_fa <- tibble(dropped = dfNames_fa,
                         MarginalR2 = c(0),
                         ConditionalR2 = c(0),
                         effects1 = c("none"),
                         effects2 = c("none"),
                         coefs1 = c(0),
                         coefs2 = c(0),
                         coefsmin1 = c(0),
                         coefsmax1 = c(0),
                         coefsmin2 = c(0),
                         coefsmax2 = c(0),
                         coefsse1 = c(0),
                         coefsse2 = c(0),
)
bsDFlist_sp <- setNames(replicate(length(bsMods_sp),
                                  bs_stats_df_sp,simplify=FALSE),
                        names(bsMods_sp))
bsDFlist_fa <- setNames(replicate(length(bsMods_fa),
                                  bs_stats_df_fa,simplify=FALSE),
                        names(bsMods_fa))
for (j in seq_along(bsMods_sp)){
  # fill tibble with stats
  for (i in seq_along(bsMods_sp[[j]])){
    tryCatch({
      bsDFlist_sp[[j]]$MarginalR2[[i]] <- tryCatch({
        r2(bsMods_sp[[j]][[i]])[[2]]
      },
      error = function(e){return(NA_real__)})
    }, error = function(e){
      bsDFlist_sp[[j]]$MarginalR2[[i]] <- NA_real__
    })
    tryCatch({
      bsDFlist_sp[[j]]$ConditionalR2[[i]] <- tryCatch({
        tryCatch({
          r2(bsMods_sp[[j]][[i]])[[1]]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$ConditionalR2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$effects1[[i]] <- tryCatch({
        tryCatch({
          names(bsMods_sp[[j]][[i]]$frame)[2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$effects1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$effects2[[i]] <- tryCatch({
        tryCatch({
          names(bsMods_sp[[j]][[i]]$frame)[3]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$effects2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefs1[[i]] <- tryCatch({
        tryCatch({
          fixef(bsMods_sp[[j]][[i]])$cond[2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefs1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefs2[[i]] <- tryCatch({
        tryCatch({
          fixef(bsMods_sp[[j]][[i]])$cond[3]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefs2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefsmin1[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_sp[[j]][[i]])[2,1]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefsmin1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefsmin2[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_sp[[j]][[i]])[3,1]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefsmin2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefsmax1[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_sp[[j]][[i]])[2,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefsmax1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefsmax2[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_sp[[j]][[i]])[3,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefsmax2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefsse1[[i]] <- tryCatch({
        tryCatch({
          summary(bsMods_sp[[j]][[i]])$coefficients$cond[2,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefsse1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_sp[[j]]$coefsse2[[i]] <- tryCatch({
        tryCatch({
          summary(bsMods_sp[[j]][[i]])$coefficients$cond[3,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_sp[[j]]$coefsse2[[i]] <- NA_real_
      })
  }
}
for (j in seq_along(bsMods_fa)){
  # fill tibble with stats
  for (i in seq_along(bsMods_fa[[j]])){
    tryCatch({
      bsDFlist_fa[[j]]$MarginalR2[[i]] <- tryCatch({
        r2(bsMods_fa[[j]][[i]])[[2]]
      },
      error = function(e){return(NA_real__)})
    }, error = function(e){
      bsDFlist_fa[[j]]$MarginalR2[[i]] <- NA_real__
    })
    tryCatch({
      bsDFlist_fa[[j]]$ConditionalR2[[i]] <- tryCatch({
        tryCatch({
          r2(bsMods_fa[[j]][[i]])[[1]]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$ConditionalR2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$effects1[[i]] <- tryCatch({
        tryCatch({
          names(bsMods_fa[[j]][[i]]$frame)[2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$effects1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$effects2[[i]] <- tryCatch({
        tryCatch({
          names(bsMods_fa[[j]][[i]]$frame)[3]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$effects2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefs1[[i]] <- tryCatch({
        tryCatch({
          fixef(bsMods_fa[[j]][[i]])$cond[2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefs1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefs2[[i]] <- tryCatch({
        tryCatch({
          fixef(bsMods_fa[[j]][[i]])$cond[3]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefs2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefsmin1[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_fa[[j]][[i]])[2,1]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefsmin1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefsmin2[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_fa[[j]][[i]])[3,1]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefsmin2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefsmax1[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_fa[[j]][[i]])[2,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefsmax1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefsmax2[[i]] <- tryCatch({
        tryCatch({
          confint(bsMods_fa[[j]][[i]])[3,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefsmax2[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefsse1[[i]] <- tryCatch({
        tryCatch({
          summary(bsMods_fa[[j]][[i]])$coefficients$cond[2,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefsse1[[i]] <- NA_real_
      })
    tryCatch({
      bsDFlist_fa[[j]]$coefsse2[[i]] <- tryCatch({
        tryCatch({
          summary(bsMods_fa[[j]][[i]])$coefficients$cond[3,2]},
          error = function(e){return(NA_real__)}
        )},
        error = function(e){return(NA_real__)}
      )},
      error = function(e){bsDFlist_fa[[j]]$coefsse2[[i]] <- NA_real_
      })
  }
}

# create mod stats tibble ####
# row bind tables for each predator into two "season" tables
spring_boot <- bind_rows(bsDFlist_sp, .id = "model")
fall_boot <- bind_rows(bsDFlist_fa, .id = "model")
# row bind tables for each season to make the final output
all_boot <- bind_rows(list("Spring" = spring_boot, "Fall" = fall_boot),
                      .id = "Season") %>%
  # replace NaN with NA
  mutate_all(~ifelse(is.nan(.), NA_real__, .)) %>%
  # make Taxon column
  mutate(Taxon = case_when(
    grepl(".ant.", model) ~"Anthocoridae",
    grepl(".ara.", model) ~"Arachnida",
    grepl(".coc.", model) ~"Coccinellidae",
    grepl(".geo.", model) ~"Geocoris",
    grepl(".ich.", model) ~"Ichneumonoidea",
  ), .after = Season)

