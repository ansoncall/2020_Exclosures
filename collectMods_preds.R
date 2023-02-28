### Spring ####
# Anthocoridae
best.nb.scaled <- get.models(nb.scaled$tab.nb.anth.sp.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.anth.sp.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.anth.sp.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.anth.sp.ranked, 1)[[1]]
anthFams.sp <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Arachnida
best.nb.scaled <- get.models(nb.scaled$tab.nb.ara.sp.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.ara.sp.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.ara.sp.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.ara.sp.ranked, 1)[[1]]
araFams.sp <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Coccinellidae
best.nb.scaled <- get.models(nb.scaled$tab.nb.cocc.sp.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.cocc.sp.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.cocc.sp.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.cocc.sp.ranked, 1)[[1]]
coccFams.sp <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Geocoris
best.nb.scaled <- get.models(nb.scaled$tab.nb.geo.sp.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.geo.sp.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.geo.sp.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.geo.sp.ranked, 1)[[1]]
geoFams.sp <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Ichneumonidae
best.nb.scaled <- get.models(nb.scaled$tab.nb.ich.sp.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.ich.sp.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.ich.sp.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.ich.sp.ranked, 1)[[1]]
ichFams.sp <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)

### Fall ####
# Anthocoridae
best.nb.scaled <- get.models(nb.scaled$tab.nb.anth.fa.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.anth.fa.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.anth.fa.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.anth.fa.ranked, 1)[[1]]
anthFams.fa <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Arachnida
best.nb.scaled <- get.models(nb.scaled$tab.nb.ara.fa.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.ara.fa.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.ara.fa.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.ara.fa.ranked, 1)[[1]]
araFams.fa <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Coccinellidae
best.nb.scaled <- get.models(nb.scaled$tab.nb.cocc.fa.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.cocc.fa.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.cocc.fa.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.cocc.fa.ranked, 1)[[1]]
coccFams.fa <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Geocoris
best.nb.scaled <- get.models(nb.scaled$tab.nb.geo.fa.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.geo.fa.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.geo.fa.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.geo.fa.ranked, 1)[[1]]
geoFams.fa <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)
# Ichneumonidae
best.nb.scaled <- get.models(nb.scaled$tab.nb.ich.fa.scaled, 1)[[1]]
best.nb.ranked <- get.models(nb.ranked$tab.nb.ich.fa.ranked, 1)[[1]]
best.pois.scaled <- get.models(pois.scaled$tab.pois.ich.fa.scaled, 1)[[1]]
best.pois.ranked <- get.models(pois.ranked$tab.pois.ich.fa.ranked, 1)[[1]]
ichFams.fa <- model.sel(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)

## clean env ####
rm(best.nb.scaled, best.nb.ranked, best.pois.scaled, best.pois.ranked)