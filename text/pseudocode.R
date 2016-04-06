# Fit models to known data:
cmsy_status <- cmsy(...)
comsir_status <- comsir(...)
# Combine and reshape output into wide format (not shown)
head(data_known)
# > cmsy_status comsir_status known_status
# > 1.4         1.3           1.2
# > ...         ...           ...
# Build super ensemble model:
ens_model <- lm(known_status ~ cmsy_status + comsir_status,
  data = data_known)

# Fit models to data of interest:
cmsy_status <- cmsy(...)
comsir_status <- comsir(...)
# Combine and reshape output into wide format (not shown)
head(data_unknown)
# > cmsy_status comsir_status ...
# > 0.9         1.1
# > ...         ...
# Predict status using super ensemble model:
predict(model, newdata = data_unknown)
# > 1.06  0.98  0.85 ...
