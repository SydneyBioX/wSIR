library(dplyr)
g1x <- rnorm(20, mean = 1, sd = 0.1)
g1y <- rnorm(20, mean = 1, sd = 0.1)
g2x <- rnorm(20, mean = 1, sd = 0.1)
g2y <- rnorm(20, mean = -1, sd = 0.1)
g3x <- rnorm(20, mean = -1, sd = 0.1)
g3y <- rnorm(20, mean = 1, sd = 0.1)

df <- matrix(NA, nrow = 60, ncol = 3) %>% as.data.frame()
colnames(df) <- c("x", "y", "colour")
df$x <- c(g1x, g2x, g3x)
df$y <- c(g1y, g2y, g3y)
df$colour <- as.character(c(rep("blue", 20), rep("red", 20), rep("purple", 20)))

plot(df$x, df$y, col = df$colour)

pc <- prcomp(df[,1:2])
pc$rotation

col_lda <- lda(colour~., data = df)
loadings_lda <- col_lda$scaling

#sir <- sir_univariate(X = df[,1:2], Y = df[,3, drop = FALSE], categorical = TRUE)

sir_package <- dr(colour~., data = df, method = "sir")
sir_package$evectors
#### direction
plot(df$x, df$y, col = df$colour, asp = 1)
points(x = c(0, loadings_lda[1, 1]),
       y = c(0, loadings_lda[2, 1]), type = "l")
abline(a = 0, b = loadings_lda[2, 1]/loadings_lda[1, 1], col = "green")
abline(a = 0, b = loadings_lda[2, 2]/loadings_lda[1, 2], col = "darkgreen")
abline(a = 0, b = pc$rotation[2, 1]/pc$rotation[1, 1])
abline(a = 0, b = pc$rotation[2, 2]/pc$rotation[1, 2], col = "grey")
abline(a = 0, b = sir_package$evectors[2, 1]/sir_package$evectors[1, 1], col = "lightblue")
abline(a = 0, b = sir_package$evectors[2, 2]/sir_package$evectors[1, 2], col = "blue")

sliced_data <- slicer(X = df[,1:2], Y = df[,3, drop = FALSE], categorical = TRUE)
sliced_data_centered <- sweep(sliced_data, 2, colMeans(df[, 1:2]), "-")
sliced_data_centered
A <- t(as.matrix(sliced_data_centered)) %*% as.matrix(sliced_data_centered)
A
eigen(A)
covX <- cov(df[, 1:2])
covX
solve(covX, eigen(A)$vectors)
sir_package <- dr(colour~., data = df, method = "sir")
sir_package$evectors

### re-code what is described in the package
X <- df[, 1:2]
n <- nrow(X)
Xc <- scale(X, center = TRUE, scale = FALSE)
qr.Xc <- qr(Xc)
R <- qr.R(qr.Xc)
Z <- qr.Q(qr.Xc) * sqrt(n) # identity sample covarinance
sliced_mean_data <- 
  as.matrix(slicer(Z, df[, 3, drop = FALSE], slices = 3, categorical = TRUE))
W <- diag(c(table(df[, 3, drop = FALSE]))/n)
M <- t(sliced_mean_data) %*% W %*% sliced_mean_data
evalues
evectors <- backsolve(R, eigen(M)$vectors)
evectors <- apply(evectors, 2, function(x) x/sqrt(sum(x^2))) # normalized

sir_package$evectors







