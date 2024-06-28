# Load necessary libraries
library(tidyverse)
library(MASS)  # For normal distribution
library(NBZIMM)  # For negative binomial distribution
library(umap) # For UMAP plot


#### Generate simulation data from concentric circles
simulation_circle <- function(num_clusters = 6, # number of clusters
                              genes_signal = 200, # number of signal genes
                              genes_noise = 200, # number of noise genes
                              radius = NULL, # radius of each circle, e.g. c(1,2,3,5). if NULL, use equally increased redius
                              base_radius = 1, # radius of the most inner circle
                              mean_expression = NULL, # mean of the signal gene, e.g. c(1,2,3,4). if NULL, use equal mean
                              dispersion = NULL, # dispersion of the signal gene, e.g. c(1,2,1,1). if NULL use equal dispersion
                              grid = T, # whether use the grid data. if False, use the randomly distribution
                              sd = NULL, # sd of noise genes, e.g. c(1,1,1,2). if NULL use equal sd
                              distribution = "NB", # type of distribution for signal genes, from c("NB", "ZINB", "P", "ZIP")
                              grid_density = 50, # number of points along one dimension
                              num_spots_base = 20, # base number for random distributed spots
                              zero_inflation_prob = 0.2, # probability of zero inflation
                              seed = 1234 # random seed
                              )
{
  set.seed(seed)
  #### Parameters
  if(is.null(radius))
    radius <- seq(1, num_clusters)
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
    max_radius <- max(radius)
    x_grid <- seq(-max_radius, max_radius, length.out = grid_density)
    y_grid <- seq(-max_radius, max_radius, length.out = grid_density)
    grid_data <- expand.grid(x = x_grid, y = y_grid) 
    radius <- c(0, radius)
    
    for (circle in 1:num_clusters) 
    {
      inner_radius <- radius[circle]
      outer_radius <- radius[circle+1]
      
      #### Filter grid points for the current circle
      current_spots <- subset(grid_data, sqrt(x^2 + y^2) >= inner_radius & sqrt(x^2 + y^2) < outer_radius)
      signal_index <- 1:genes_signal
      noise_index <- (genes_signal+1):(genes_signal+genes_noise)
      
      for (i in 1:nrow(current_spots)) 
      {
        #### Location data
        x <- current_spots$x[i]
        y <- current_spots$y[i]
        spots <- rbind(spots, data.frame(x = x, y = y, cluster = as.factor(circle)))
        
        #### Gene expression data 
        if(distribution == "NB")
          signal <- rnbinom(n = genes_signal, size = mean_expression[circle] / dispersion[circle], mu = mean_expression[circle])
        else if(distribution == "ZINB")
          signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                           rnbinom(n = genes_signal, size = mean_expression[circle] / dispersion[circle], mu = mean_expression[circle]))
        else if(distribution == "P")
          signal <- rpois(genes_signal, lambda = mean_expression[circle])
        else if(distribution == "ZIP")
          signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                           rpois(genes_signal, lambda = mean_expression[circle]))
        noise <- round(abs(rnorm(genes_noise, mean = 0, sd = sd[circle])))
        expression_tmp <- rep(NA, genes_noise+genes_signal)
        expression_tmp[signal_index] <- signal
        expression_tmp[noise_index] <- noise
        expressions <- rbind(expressions, expression_tmp)
      }
    }
  }
  else
  {
    for (circle in 1:num_clusters) 
    {
      signal_index <- 1:genes_signal
      noise_index <- (genes_signal+1):(genes_signal+genes_noise)
      if (circle == 1) 
      {
        #### For the most inner circle, generate a solid circle
        num_spots <- ceiling(pi * (base_radius^2) * num_spots_base)  # Area proportional
        inner_radius <- 0
        outer_radius <- base_radius
      } 
      else 
      {
        #### For the remaining circles, generate as rings
        inner_radius <- radius[circle-1]
        outer_radius <- radius[circle]
        num_spots <- ceiling(pi * (outer_radius^2 - inner_radius^2) * num_spots_base)  # Area proportional
      }
      
      for (i in 1:num_spots) 
      {
        #### Location data
        r <- sqrt(runif(1, inner_radius^2, outer_radius^2))
        theta <- runif(1, 0, 2 * pi)
        x <- r * cos(theta)
        y <- r * sin(theta)
        spots <- rbind(spots, data.frame(x = x, y = y, cluster = as.factor(circle)))
        
        # Gene expression data 
        if(distribution == "NB")
          signal <- rnbinom(n = genes_signal, size = mean_expression[circle] / dispersion[circle], mu = mean_expression[circle])
        else if(distribution == "ZINB")
          signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                           rnbinom(n = genes_signal, size = mean_expression[circle] / dispersion[circle], mu = mean_expression[circle]))
        else if(distribution == "P")
          signal <- rpois(genes_signal, lambda = mean_expression[circle])
        else if(distribution == "ZIP")
          signal <- ifelse(runif(genes_signal) < zero_inflation_prob, 0, 
                           rpois(genes_signal, lambda = mean_expression[circle]))
        noise <- round(abs(rnorm(genes_noise, mean = 0, sd = sd[circle])))
        expression_tmp <- rep(NA, genes_noise+genes_signal)
        expression_tmp[signal_index] <- signal
        expression_tmp[noise_index] <- noise
        expressions <- rbind(expressions, expression_tmp)
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
  
  folder_path <- paste0(work_path, "circle_grid_disp_", disp, "/")
  if(! exists(folder_path))
    dir.create(folder_path)
  data = simulation_circle(grid = T, mean_expression = mean_expression, dispersion = rep(disp, 6))
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
