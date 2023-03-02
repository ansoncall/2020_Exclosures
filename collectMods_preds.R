### Spring ####
# Anthocoridae
best_nb_scaled <- get.models(nb_scaled$tab_nb_anth_sp_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_anth_sp_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_anth_sp_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_anth_sp_ranked, 1)[[1]]
anth_fams_sp <- model.sel(best_nb_scaled,
                         best_nb_ranked,
                         best_pois_scaled,
                         best_pois_ranked)
# Arachnida
best_nb_scaled <- get.models(nb_scaled$tab_nb_ara_sp_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_ara_sp_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_ara_sp_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_ara_sp_ranked, 1)[[1]]
ara_fams_sp <- model.sel(best_nb_scaled,
                        best_nb_ranked,
                        best_pois_scaled,
                        best_pois_ranked)
# Coccinellidae
best_nb_scaled <- get.models(nb_scaled$tab_nb_cocc_sp_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_cocc_sp_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_cocc_sp_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_cocc_sp_ranked, 1)[[1]]
cocc_fams_sp <- model.sel(best_nb_scaled,
                         best_nb_ranked,
                         best_pois_scaled,
                         best_pois_ranked)
# Geocoris
best_nb_scaled <- get.models(nb_scaled$tab_nb_geo_sp_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_geo_sp_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_geo_sp_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_geo_sp_ranked, 1)[[1]]
geo_fams_sp <- model.sel(best_nb_scaled,
                        best_nb_ranked,
                        best_pois_scaled,
                        best_pois_ranked)
# Ichneumonidae
best_nb_scaled <- get.models(nb_scaled$tab_nb_ich_sp_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_ich_sp_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_ich_sp_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_ich_sp_ranked, 1)[[1]]
ich_fams_sp <- model.sel(best_nb_scaled,
                        best_nb_ranked,
                        best_pois_scaled,
                        best_pois_ranked)

### Fall ####
# Anthocoridae
best_nb_scaled <- get.models(nb_scaled$tab_nb_anth_fa_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_anth_fa_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_anth_fa_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_anth_fa_ranked, 1)[[1]]
anth_fams_fa <- model.sel(best_nb_scaled,
                         best_nb_ranked,
                         best_pois_scaled,
                         best_pois_ranked)
# Arachnida
best_nb_scaled <- get.models(nb_scaled$tab_nb_ara_fa_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_ara_fa_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_ara_fa_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_ara_fa_ranked, 1)[[1]]
ara_fams_fa <- model.sel(best_nb_scaled,
                        best_nb_ranked,
                        best_pois_scaled,
                        best_pois_ranked)
# Coccinellidae
best_nb_scaled <- get.models(nb_scaled$tab_nb_cocc_fa_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_cocc_fa_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_cocc_fa_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_cocc_fa_ranked, 1)[[1]]
cocc_fams_fa <- model.sel(best_nb_scaled,
                         best_nb_ranked,
                         best_pois_scaled,
                         best_pois_ranked)
# Geocoris
best_nb_scaled <- get.models(nb_scaled$tab_nb_geo_fa_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_geo_fa_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_geo_fa_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_geo_fa_ranked, 1)[[1]]
geo_fams_fa <- model.sel(best_nb_scaled,
                        best_nb_ranked,
                        best_pois_scaled,
                        best_pois_ranked)
# Ichneumonidae
best_nb_scaled <- get.models(nb_scaled$tab_nb_ich_fa_scaled, 1)[[1]]
best_nb_ranked <- get.models(nb_ranked$tab_nb_ich_fa_ranked, 1)[[1]]
best_pois_scaled <- get.models(pois_scaled$tab_pois_ich_fa_scaled, 1)[[1]]
best_pois_ranked <- get.models(pois_ranked$tab_pois_ich_fa_ranked, 1)[[1]]
ich_fams_fa <- model.sel(best_nb_scaled,
                        best_nb_ranked,
                        best_pois_scaled,
                        best_pois_ranked)

## clean env ####
rm(best_nb_scaled, best_nb_ranked, best_pois_scaled, best_pois_ranked)
