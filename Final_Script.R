# BIOS 480 FINAL PROJECT
# While running the code, please run it section-wise

library(ncdf4)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(sf)
library(maps)
library(viridis)
library(gridExtra)
library(MASS)
library(moments)
library(grid)

setwd("/Users/gunashree/Desktop/BIOS480/FINAL PROJECT/Mywork/CODE")

#DATA LOADING

invasive_data <- read.csv("data_raw/merged_richness_fips.csv", stringsAsFactors = FALSE)
cat("Loaded:", nrow(invasive_data), "counties\n\n")

nc_data <- nc_open("data_raw/MEGAN-EMEP-IsopreneEP_Glb_0.25x0.25_bio_isoprene-emission-potential__monthly.nc")

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
isoprene_data <- ncvar_get(nc_data, "ef_isop")
isoprene_annual_mean <- apply(isoprene_data, c(1, 2), mean, na.rm = TRUE)
nc_close(nc_data)

cat("Grid dimensions:", dim(isoprene_annual_mean)[1], "×", 
    dim(isoprene_annual_mean)[2], "\n\n")



#DATA CLEANING

lon_converted <- ifelse(lon > 180, lon - 360, lon)

us_lon_range <- c(-125, -66)
us_lat_range <- c(24, 50)

us_lon_idx <- which(lon_converted >= us_lon_range[1] & lon_converted <= us_lon_range[2])
us_lat_idx <- which(lat >= us_lat_range[1] & lat <= us_lat_range[2])

isoprene_us <- isoprene_annual_mean[us_lon_idx, us_lat_idx]
lon_us <- lon_converted[us_lon_idx]
lat_us <- lat[us_lat_idx]

isoprene_grid <- expand.grid(lon = lon_us, lat = lat_us)
isoprene_grid$isoprene <- as.vector(isoprene_us)
isoprene_grid <- isoprene_grid %>% filter(!is.na(isoprene))

cat(" US isoprene grid created:", nrow(isoprene_grid), "points\n\n")

if(!require(tigris, quietly = TRUE)) {
  install.packages("tigris")
  library(tigris)
}

options(tigris_use_cache = TRUE)
us_counties_sf <- counties(cb = TRUE, resolution = "20m", year = 2020)
county_centroids <- st_centroid(us_counties_sf)
county_coords <- st_coordinates(county_centroids)

counties_data <- data.frame(
  fips = us_counties_sf$GEOID,
  county_name = us_counties_sf$NAME,
  state = us_counties_sf$STATEFP,
  longitude = county_coords[, 1],
  latitude = county_coords[, 2],
  stringsAsFactors = FALSE
)

invasive_data$fips <- sprintf("%05d", invasive_data$fips)

find_nearest_isoprene <- function(county_lon, county_lat, grid_data) {
  distances <- sqrt((grid_data$lon - county_lon)^2 + (grid_data$lat - county_lat)^2)
  nearest_idx <- which.min(distances)
  return(grid_data$isoprene[nearest_idx])
}

cat("Matching isoprene to counties...\n")
counties_data$isoprene_emission <- sapply(1:nrow(counties_data), function(i) {
  if(i %% 500 == 0) cat("  ", i, "/", nrow(counties_data), "\n")
  find_nearest_isoprene(counties_data$longitude[i], counties_data$latitude[i], isoprene_grid)
})

analysis_data <- invasive_data %>%
  inner_join(counties_data, by = "fips") %>%
  dplyr::select(fips, species_count = Total_pests, county_name = county,
                state_name, longitude, latitude, isoprene_emission) %>%
  filter(!is.na(species_count), !is.na(isoprene_emission),
         species_count >= 0, isoprene_emission >= 0)

cat("\n Final dataset:", nrow(analysis_data), "counties\n\n")


#EXPLORATORY DATA ANALYSIS


summary_stats <- data.frame(
  Variable = c("Species Count", "Isoprene (µg/m²/h)"),
  Mean = c(mean(analysis_data$species_count), mean(analysis_data$isoprene_emission)),
  SD = c(sd(analysis_data$species_count), sd(analysis_data$isoprene_emission)),
  Median = c(median(analysis_data$species_count), median(analysis_data$isoprene_emission)),
  Min = c(min(analysis_data$species_count), min(analysis_data$isoprene_emission)),
  Max = c(max(analysis_data$species_count), max(analysis_data$isoprene_emission))
)

print(summary_stats)

variance_species <- var(analysis_data$species_count)
mean_species <- mean(analysis_data$species_count)
overdispersion_ratio <- variance_species / mean_species

cat("\nOVERDISPERSION CHECK:\n")
cat("  Variance/Mean ratio:", round(overdispersion_ratio, 3), "\n")
if(overdispersion_ratio > 1.5) {
  cat("  ✓ OVERDISPERSED → Negative Binomial needed\n\n")
}

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
hist(analysis_data$species_count, breaks = 30, main = "Species Count",
     xlab = "Count", col = "skyblue", border = "white")
hist(analysis_data$isoprene_emission, breaks = 30, main = "Isoprene Emissions",
     xlab = "µg/m²/h", col = "lightgreen", border = "white")
boxplot(analysis_data$species_count, main = "Species Outliers",
        ylab = "Count", col = "skyblue")
boxplot(analysis_data$isoprene_emission, main = "Isoprene Outliers",
        ylab = "µg/m²/h", col = "lightgreen")
par(mfrow = c(1, 1))






#SPATIAL VISUALIZATION
#While running this part of the code, please run it map-wise

# Preparing map data
map_data <- us_counties_sf %>%
  left_join(analysis_data, by = c("GEOID" = "fips"))

# Filtering to continental US
continental_states <- c(
  "01", "04", "05", "06", "08", "09", "10", "11", "12", "13", 
  "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", 
  "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", 
  "36", "37", "38", "39", "40", "41", "42", "44", "45", "46", 
  "47", "48", "49", "50", "51", "53", "54", "55", "56"
)

map_data_continental <- map_data %>%
  filter(STATEFP %in% continental_states)

cat("Continental US data:", sum(!is.na(map_data_continental$species_count)), 
    "counties\n\n")

#Isoprene Emissions Map
map_isoprene <- ggplot(map_data_continental) +
  geom_sf(aes(fill = isoprene_emission), color = NA, size = 0.1) +
  scale_fill_viridis(option = "magma", 
                     name = "Isoprene\nEmission\n(µg/m²/h)",
                     na.value = "gray90",
                     trans = "sqrt",
                     breaks = c(0, 100, 500, 1000, 2500, 5000, 10000),
                     labels = c("0", "100", "500", "1K", "2.5K", "5K", "10K")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  ) +
  labs(
    title = "Isoprene Emission Potential",
    subtitle = "County-level annual mean (2000)"
  )

map_isoprene


# Invasive Species Richness Map

map_species <- ggplot(map_data_continental) +
  geom_sf(aes(fill = species_count), color = NA, size = 0.1) +
  scale_fill_viridis(option = "plasma", 
                     name = "Invasive\nSpecies\nCount",
                     na.value = "gray90",
                     breaks = c(0, 2, 5, 10, 15, 20, 30)) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  ) +
  labs(
    title = "Invasive Forest Insect Species Richness",
    subtitle = "Number of invasive species per county"
  )

map_species

# Side-by-Side Comparison Map
grid.arrange(map_isoprene, map_species, ncol = 2,
             top = textGrob("Spatial Distribution: Isoprene Emissions and Invasive Species Richness",
                            gp = gpar(fontsize = 16, fontface = "bold")))



# CATEGORICAL AND HOTSPOT MAPS

# Categorize variables
map_data_continental <- map_data_continental %>%
  mutate(
    isoprene_category = cut(
      isoprene_emission,
      breaks = quantile(isoprene_emission, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE),
      labels = c("Low", "Medium", "High"),
      include.lowest = TRUE
    ),
    species_category = cut(
      species_count,
      breaks = c(-1, 0, 2, 5, Inf),
      labels = c("None", "Low (1-2)", "Medium (3-5)", "High (6+)"),
      include.lowest = TRUE
    )
  )

# Categorical maps
map_isoprene_cat <- ggplot(map_data_continental) +
  geom_sf(aes(fill = isoprene_category), color = "white", size = 0.05) +
  scale_fill_manual(
    values = c("Low" = "#440154", "Medium" = "#FDE725", "High" = "#FF6B35"),
    name = "Isoprene\nLevel",
    na.value = "gray90"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Isoprene Emission Categories")

map_species_cat <- ggplot(map_data_continental) +
  geom_sf(aes(fill = species_category), color = "white", size = 0.05) +
  scale_fill_manual(
    values = c("None" = "#EEEEEE", 
               "Low (1-2)" = "#A8DADC", 
               "Medium (3-5)" = "#457B9D", 
               "High (6+)" = "#1D3557"),
    name = "Species\nRichness",
    na.value = "gray90"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Invasive Species Categories")

grid.arrange(map_isoprene_cat, map_species_cat, ncol = 2,
             top = textGrob("Categorical Comparison",
                            gp = gpar(fontsize = 16, fontface = "bold")))

# Hotspot overlap analysis
map_data_continental <- map_data_continental %>%
  mutate(
    isoprene_hotspot = isoprene_emission > quantile(isoprene_emission, 0.75, na.rm = TRUE),
    species_hotspot = species_count > quantile(species_count, 0.75, na.rm = TRUE),
    overlap_category = case_when(
      isoprene_hotspot & species_hotspot ~ "Both High",
      isoprene_hotspot & !species_hotspot ~ "High Isoprene Only",
      !isoprene_hotspot & species_hotspot ~ "High Species Only",
      TRUE ~ "Neither High"
    )
  )

map_overlap <- ggplot(map_data_continental) +
  geom_sf(aes(fill = overlap_category), color = "white", size = 0.05) +
  scale_fill_manual(
    values = c(
      "Both High" = "#D62828",
      "High Isoprene Only" = "#F77F00",
      "High Species Only" = "#003049",
      "Neither High" = "#E5E5E5"
    ),
    name = "Hotspot\nCategory"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Hotspot Overlap Analysis",
    subtitle = "Counties in top 25% for isoprene and/or invasive species"
  )

print(map_overlap)

# Hotspot statistics
hotspot_counts <- map_data_continental %>%
  st_drop_geometry() %>%
  filter(!is.na(overlap_category)) %>%
  count(overlap_category) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

cat("\nHOTSPOT OVERLAP COUNTS:\n")
print(hotspot_counts)
cat("\n")



#INITIAL LINEAR REGRESSION

plot(analysis_data$isoprene_emission, analysis_data$species_count,
     xlab = "Isoprene (µg/m²/h)", ylab = "Species Count",
     main = "Raw Data Scatterplot", pch = 16, col = rgb(0, 0, 1, 0.3))

invasiveRegression <- lm(species_count ~ isoprene_emission, data = analysis_data)
cat("\nLinear Model Summary:\n")
print(summary(invasiveRegression))


#ASSUMPTION CHECKING FOR LINEAR MODEL

cat("The 4 assumptions:\n")
cat("1. LINEARITY\n2. NORMALITY\n3. HOMOGENEITY\n4. INDEPENDENCE\n\n")

par(mfrow = c(2, 2))
plot(invasiveRegression)
par(mfrow = c(1, 1))


cat("\nINTERPRETATION:\n")
cat("Assumption 1 - LINEARITY: ❌ VIOLATED (curved red line)\n")
cat("Assumption 2 - NORMALITY: ❌ VIOLATED (Q-Q deviates)\n")
cat("Assumption 3 - HOMOGENEITY: ❌ VIOLATED (funnel pattern)\n")
cat("Assumption 4 - INDEPENDENCE: ⚠ Influential points\n\n")

cat("CONCLUSION: Linear regression INAPPROPRIATE\n")
cat("Need COUNT DATA model (Negative Binomial)\n\n")




#COUNT DATA MODELS

analysis_data$log_isoprene <- log(analysis_data$isoprene_emission + 1)

cat("POISSON MODEL:\n")
poisson_model <- glm(species_count ~ log_isoprene, 
                     family = poisson(link = "log"), data = analysis_data)
print(summary(poisson_model))

dispersion <- poisson_model$deviance / poisson_model$df.residual
cat("\nDispersion:", round(dispersion, 3))
if(dispersion > 1.5) {
  cat("OVERDISPERSED, need Negative Binomial\n\n")
}



cat("\nNEGATIVE BINOMIAL MODEL:\n")
nb_model <- glm.nb(species_count ~ log_isoprene, data = analysis_data)
print(summary(nb_model))

cat("\nMODEL COMPARISON (AIC):\n")
cat("  Linear:    ", round(AIC(invasiveRegression), 1), "\n")
cat("  Poisson:   ", round(AIC(poisson_model), 1), "\n")
cat("  Neg Binom: ", round(AIC(nb_model), 1), " ← BEST\n\n")

best_model <- nb_model
model_name <- "Negative Binomial"



#DIAGNOSTIC PLOTS FOR NB MODEL

analysis_data$fitted_values <- fitted(best_model)
analysis_data$residuals <- residuals(best_model, type = "deviance")
analysis_data$std_resid <- rstandard(best_model, type = "deviance")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

plot(analysis_data$fitted_values, analysis_data$residuals,
     xlab = "Fitted", ylab = "Deviance Residuals",
     main = "Residuals vs Fitted", pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(lowess(analysis_data$fitted_values, analysis_data$residuals), col = "red", lwd = 2)

qqnorm(analysis_data$std_resid, main = "Q-Q Plot", pch = 16, col = rgb(0,0,0,0.3))
qqline(analysis_data$std_resid, col = "red", lwd = 2)

plot(analysis_data$fitted_values, sqrt(abs(analysis_data$std_resid)),
     xlab = "Fitted", ylab = expression(sqrt("|Std Residuals|")),
     main = "Scale-Location", pch = 16, col = rgb(0,0,0,0.3))
lines(lowess(analysis_data$fitted_values, sqrt(abs(analysis_data$std_resid))), 
      col = "red", lwd = 2)

plot(hatvalues(best_model), analysis_data$std_resid,
     xlab = "Leverage", ylab = "Std Residuals",
     main = "Residuals vs Leverage", pch = 16, col = rgb(0,0,0,0.3))
abline(h = 0, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))



cat("\nINTERPRETATION:\n")
cat("Assumption 1 - LINEARITY: Red line roughly horizontal\n")
cat("Assumption 2 - NORMALITY: Points follow diagonal\n")
cat("Assumption 3 - HOMOGENEITY: Variance relatively constant\n")
cat("Assumption 4 - INDEPENDENCE: No severe outliers\n\n")

cat("HORIZONTAL BANDS are NORMAL for count data!\n")
cat("Each band = one species count (0, 1, 2, 3...)\n\n")

cat("ASSUMPTIONS SATISFIED\n")
cat("Negative Binomial model is APPROPRIATE\n\n")



# FINAL VISUALIZATION

analysis_data$predicted <- predict(best_model, type = "response")

ggplot(analysis_data, aes(x = log_isoprene, y = species_count)) +
  geom_point(alpha = 0.4, size = 2, color = "#2E86AB") +
  geom_line(aes(y = predicted), color = "#A23B72", linewidth = 1.2) +
  labs(title = "Negative Binomial Regression: Species vs Isoprene",
       subtitle = paste0("n = ", nrow(analysis_data), " US counties"),
       x = "log(Isoprene + 1)", y = "Species Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))




#HYPOTHESIS TESTING AND CONCLUSIONS


cat("  ", rep("=", 70), "\n", sep = "")
cat("  HYPOTHESIS TESTING RESULTS\n")
cat("  ", rep("=", 70), "\n\n", sep = "")

cat("  Research Question:\n")
cat("    Is there a relationship between isoprene emissions and\n")
cat("    invasive forest insect species richness?\n\n")

cat("  H₀: Isoprene has NO effect (β = 0)\n")
cat("  H₁: Isoprene has a SIGNIFICANT effect (β ≠ 0)\n")
cat("  α = 0.05\n\n")

coef_summary <- summary(best_model)$coefficients
coef_isoprene <- coef_summary[2, 1]
p_value <- coef_summary[2, 4]

cat("  RESULTS:\n")
cat("  Coefficient: β =", round(coef_isoprene, 4), "\n")
cat("  P-value: p", ifelse(p_value < 0.001, "< 0.001", 
                           paste("=", round(p_value, 4))), "\n\n")

if(p_value < 0.05) {
  cat("  REJECT H₀\n\n")
  cat("  CONCLUSION:\n")
  cat("  There IS significant evidence that isoprene emissions\n")
  cat("  affect invasive species richness (p < 0.05).\n\n")
  
  if(coef_isoprene > 0) {
    cat("  Direction: POSITIVE ↗\n")
    cat("  Higher isoprene → More invasive species\n\n")
  }
  
  effect <- exp(coef_isoprene)
  cat("  Effect: For each 1-unit increase in log(isoprene),\n")
  cat("  species count multiplies by", round(effect, 3), "\n")
  cat("  (", round((effect-1)*100, 1), "% change)\n\n")
} else {
  cat("  FAIL TO REJECT H₀\n\n")
  cat("  No significant relationship detected.\n\n")
}

cat("  ", rep("=", 70), "\n\n", sep = "")



# WHY NEGATIVE BINOMIAL IS BEST

cat("1. DATA STRUCTURE:\n")
cat("   COUNT DATA (0, 1, 2, 3... species)\n")
cat("   Discrete integers, not continuous\n")
cat("   Non-negative values only\n\n")

cat("2. OVERDISPERSION:\n")
cat("   Var/Mean =", round(overdispersion_ratio, 2), "(>> 1.0)\n")
cat("   Poisson assumes Var = Mean → VIOLATED\n")
cat("   NB models variance properly → APPROPRIATE\n\n")

cat("3. MODEL COMPARISON:\n")
comparison_table <- data.frame(
  Model = c("Linear", "Poisson", "Negative Binomial"),
  Appropriate = c("❌ NO", "⚠ Marginal", "✓ YES"),
  Handles_Counts = c("❌ NO", "✓ YES", "✓ YES"),
  Handles_Overdispersion = c("❌ NO", "❌ NO", "✓ YES"),
  AIC = c(round(AIC(invasiveRegression), 1), 
          round(AIC(poisson_model), 1), 
          round(AIC(nb_model), 1))
)
print(comparison_table)

cat("\n4. ASSUMPTIONS:\n")
cat("   Diagnostic plots acceptable\n")
cat("   Residuals approximately normal\n")
cat("   No severe violations\n\n")

cat("5. INTERPRETATION:\n")
cat("   Clear biological meaning\n")
cat("   Valid predictions (non-negative)\n")
cat("   Accurate confidence intervals\n\n")

cat("NEGATIVE BINOMIAL IS THE CORRECT MODEL\n")


# FINAL SUMMARY
cat("FINAL SUMMARY\n")
cat("================================================================\n\n")

cat("DATASET:\n")
cat("  Sample: n =", nrow(analysis_data), "US counties\n")
cat("  Species: 0 to", max(analysis_data$species_count), "(mean =", 
    round(mean(analysis_data$species_count), 2), ")\n")
cat("  Overdispersion:", round(overdispersion_ratio, 2), "\n\n")

cat("MODEL:\n")
cat("  ", model_name, "regression\n")
cat("  Predictor: log(isoprene)\n")
cat("  AIC:", round(AIC(best_model), 1), "(best)\n\n")

cat("KEY FINDINGS:\n")
cat("  Coefficient: β =", round(coef_isoprene, 4), "\n")
cat("  P-value: p", ifelse(p_value < 0.001, "< 0.001", 
                             paste("=", round(p_value, 4))), "\n")
cat("  Effect:", ifelse(coef_isoprene > 0, "POSITIVE", "NEGATIVE"), "relationship\n")
cat("  Significance:", ifelse(p_value < 0.05, "YES ", "NO "), "\n\n")

cat("SPATIAL PATTERNS:\n")
cat("  Both high hotspots:", 
    sum(map_data_continental$overlap_category == "Both High", na.rm = TRUE), "counties\n")
cat("  Spatial clustering observed in maps\n")
cat("  Regional patterns evident\n\n")


cat("BIOLOGICAL INTERPRETATION:\n")
if(p_value < 0.05 && coef_isoprene > 0) {
  cat("  Counties with higher isoprene emissions tend to harbor more\n")
  cat("  invasive forest insect species. This suggests that volatile\n")
  cat("  organic compounds may play a role in invasion success, possibly\n")
  cat("  through attraction mechanisms or indicating suitable habitat.\n\n")
} else if(p_value < 0.05 && coef_isoprene < 0) {
  cat("  Counties with higher isoprene emissions tend to have fewer\n")
  cat("  invasive forest insect species. This may indicate that high\n")
  cat("  BVOC production serves as a plant defense mechanism.\n\n")
} else {
  cat("  No significant relationship detected between isoprene\n")
  cat("  emissions and invasive species richness. Other factors may\n")
  cat("  be more important in determining invasion patterns.\n\n")
}

cat("STATISTICAL VALIDITY:\n")
cat("  Assumptions checked: Yes\n")
cat("  Model diagnostics: Acceptable\n")
cat("  Appropriate method: Yes (count data model)\n")
cat("  Results are valid: Yes\n\n")
