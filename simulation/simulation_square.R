# Load necessary libraries
library(tidyverse)
library(MASS)  # For normal distribution
library(NBZIMM)  # For negative binomial distribution
library(umap) # For UMAP plot


#### Generate simulation data from concentric circles
simulation_square <- function(num_rows = 2, # number of squares per row
                              num_cols = 3, # number of squares per col
                              len_width = 1, # width of each square
                              len_height = 1, # height of each square
                              genes_signal = 200, # number of signal genes
                              genes_noise = 200, # number of noise genes
                              mean_expression = NULL, # mean of the signal gene, e.g. c(1,2,3,4). if NULL, use equal mean
                              dispersion = NULL, # dispersion of the signal gene, e.g. c(1,2,1,1). if NULL use equal dispersion
                              grid = T, # whether use the grid data. if False, use the randomly distribution
                              sd = NULL, # sd of noise genes, e.g. c(1,1,1,2). if NULL use equal sd
                              distribution = "NB", # type of distribution for signal genes, from c("NB", "ZINB", "P", "ZIP")
                              grid_density = 50, # number of points along one dimension
                              num_spots_base = 400, # number of spots per square
                              zero_inflation_prob = 0.2, # probability of zero inflation
                              seed = 1234 # random seed
                              )
{
  set.seed(seed)
  #### Parameters
  num_clusters <- num_rows * num_cols
  if(is.null(mean_expression))
    mean_expression <- rep(1, num_clusters)
  if(is.null(dispersion))
    dispersion <- rep(1, num_clusters)
  if(is.null(sd))
    sd <- rep(0.8, num_clusters)
  
  #### Generate data
  spots <- data.frame(x = numeric(), y = numeric(), cluster = factor())
  expressions <- matrix(nrow = 0, ncol = genes_signal + genes_noise)
  
  if(grid)
  {
    max_grid_x <- num_cols * len_width
    max_grid_y <- num_rows * len_height
    x_grid <- seq(0, max_grid_x, length.out = grid_density)
    y_grid <- seq(0, max_grid_y, length.out = grid_density)
    grid_data <- expand.grid(x = x_grid, y = y_grid) 
    
    for(row in 1:num_rows) 
    {
      for(col in 1:num_cols)
      {
        square <- (row-1) * num_cols + col
        x_min <- (col-1) * len_width
        x_max <- col * len_width
        y_min <- (row-1) * len_height
        y_max <- row * len_height
        
        #### Filter grid points for the current square
        current_spots <- subset(grid_data, (between(x, x_min, x_max)) & (between(y, y_min, y_max)))
        signal_index <- 1:genes_signal
        noise_index <- (genes_signal+1):(genes_signal+genes_noise)
        
        for(i in 1:nrow(current_spots)) 
        {
          #### Location data
          x <- current_spots$x[i]
          y <- current_spots$y[i]
          spots <- rbind(spots, data.frame(x = x, y = y, cluster = as.factor(square)))
          
          #### Gene expression data 
          if(distribution == "NB")
            signal <- rnbinom(n = genes_signal, size = mean_expression[square] / dispersion[square], mu = mean_expression[square])
          else if(distribution == "ZINB")
            signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                             rnbinom(n = genes_signal, size = mean_expression[square] / dispersion[square], mu = mean_expression[square]))
          else if(distribution == "P")
            signal <- rpois(genes_signal, lambda = mean_expression[square])
          else if(distribution == "ZIP")
            signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                             rpois(genes_signal, lambda = mean_expression[square]))
          noise <- round(abs(rnorm(genes_noise, mean = 0, sd = sd[square])))
          expression_tmp <- rep(NA, genes_noise+genes_signal)
          expression_tmp[signal_index] <- signal
          expression_tmp[noise_index] <- noise
          expressions <- rbind(expressions, expression_tmp)
        }
      }
    }
  }
  else
  {
    for(row in 1:num_rows) 
    {
      for(col in 1:num_cols)
      {
        signal_index <- 1:genes_signal
        noise_index <- (genes_signal+1):(genes_signal+genes_noise)
        
        square <- (row-1) * num_cols + col
        x_min <- (col-1) * len_width
        x_max <- col * len_width
        y_min <- (row-1) * len_height
        y_max <- row * len_height
        
        for (i in 1:num_spots_base) 
        {
          #### Location data
          x <- runif(1, x_min, x_max)
          y <- runif(1, y_min, y_max)
          spots <- rbind(spots, data.frame(x = x, y = y, cluster = as.factor(square)))
          
          # Gene expression data 
          if(distribution == "NB")
            signal <- rnbinom(n = genes_signal, size = mean_expression[square] / dispersion[square], mu = mean_expression[square])
          else if(distribution == "ZINB")
            signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                             rnbinom(n = genes_signal, size = mean_expression[square] / dispersion[square], mu = mean_expression[square]))
          else if(distribution == "P")
            signal <- rpois(genes_signal, lambda = mean_expression[square])
          else if(distribution == "ZIP")
            signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                             rpois(genes_signal, lambda = mean_expression[square]))
          noise <- round(abs(rnorm(genes_noise, mean = 0, sd = sd[square])))
          expression_tmp <- rep(NA, genes_noise+genes_signal)
          expression_tmp[signal_index] <- signal
          expression_tmp[noise_index] <- noise
          expressions <- rbind(expressions, expression_tmp)
        } 
      }
    }
  }
  
  return(list(spatial = spots, expression = expressions))
}

disp_list <- c(2) * 0.1
mean_expression <- c(1,2.5,1.5,3.5,2,3)
for(disp in disp_list)
{
  work_path <- "/Users/james/Documents/lab/2024/WEST/simulation/"
  
  folder_path <- paste0(work_path, "square_grid_disp_", disp, "/")
  if(! exists(folder_path))
    dir.create(folder_path)
  data = simulation_square(grid = T, mean_expression = mean_expression, dispersion = rep(disp, 6))
  write.csv(data$spatial, paste0(folder_path, "S.csv"), row.names = F)
  write.csv(data$expression, paste0(folder_path, "X.csv"), row.names = F)
  
  # Plotting
  clusters <- unique(data$spatial$cluster)
  colors <- RColorBrewer::brewer.pal(length(clusters), "Set1")  
  color_palette <- setNames(colors, clusters)
  
  p1 <- ggplot(data$spatial, aes(x = x, y = y, color = cluster)) + 
          geom_point(alpha = 0.6) +
          theme_minimal() +
          scale_color_manual(values = color_palette)
  ggsave(paste0(folder_path, "distribution.pdf"), p1)
  
  umap_result <- umap(data$expression)
  umap_data <- data.frame(UMAP1 = umap_result$layout[, 1], 
                          UMAP2 = umap_result$layout[, 2], 
                          Cluster = data$spatial$cluster)
  p2 <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
          geom_point(alpha = 0.8) +
          theme_minimal() +
          scale_color_manual(values = color_palette)
  ggsave(paste0(folder_path, "umap.pdf"), p2)
}

p1
p2



