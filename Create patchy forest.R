########################
# Create patchy forest #
########################


# ---- Parameters ----
set.seed(123)
library(terra)

n <- 20                 # grid size
forest_prop <- 0.2      # target proportion forest
smooth_strength <- 7    # higher = larger patches


# ---- Create random noise ----
m <- matrix(runif(n * n), nrow = n, ncol = n)

# ---- Smooth using focal averaging ----


r <- rast(m)

# Create smoothing kernel
w <- matrix(1, nrow = smooth_strength, ncol = smooth_strength)

# Apply smoothing
r_smooth <- focal(r, w = w, fun = mean, na.rm = TRUE)

# Convert to binary forest/non-forest
threshold <- quantile(values(r_smooth), probs = 1 - forest_prop)
forest <- r_smooth > threshold

# ---- Plot ----
plot(forest,
     col = c("lightgoldenrod", "darkgreen"),
     legend = FALSE,
     main = "Simulated Patchy Forest Habitat")

forest <- c(as.matrix(forest)*1)
