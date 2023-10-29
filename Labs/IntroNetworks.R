# Lab 1: Basic network stuff

# Ideas for questions
  # Find data, is the included edge definition meaningful? If not, make your own

# To do
  # X Add the node attribute volume_from to the graph
  # X Transform weights to be unit magnitude, not extremely wide magnitude before plotting
  # X PCA, eigen
  # X Heatmap
  # Network summaries
  # X Add Nimble logistic regression with edge attributes
  # Clean up and convert to .Rmd

# Bonus: 
  # Overlay network on world map: https://www.r-bloggers.com/2018/05/three-ways-of-visualizing-a-graph-on-a-map/
  # Find a second data set

## Set directories
data_path <- 'Data/'

# Load packages
library(tidyverse)
library(here)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(pheatmap)
library(nimble)

#-------------------------- PART 1: Reading and transforming the data ----------------------------------#

# Load the files
file_path1 <- here::here(paste0(data_path, "2021Trade/trade2021.txt"))
trade <- read.table(file_path1, sep = '\t',header = TRUE, skip = 0 ,quote='', comment='')
file_path2 <- here::here(paste0(data_path, "2021Trade/tradecountries.txt"))
trade.countries <- read.csv(file_path2, sep = '\t',header = TRUE, skip = 0 ,quote='', comment='')

# Look at the data: note that it's a directed, weighted edge list
head(trade)

# Note we have more to countries than from and only from countries are in our country list
unique(trade$from)
unique(trade$to)
trade.countries

# Now: we still have a lot of edges
nC.from <- length(unique(trade$from))
nC.to <- length(unique(trade$to))
nC.from*(nC.to - 1)
nrow(trade)

g <- graph_from_data_frame(trade, directed = TRUE)
plot(g) # This is kind of a mess

### How to turn this into a reasonably straight forward example
# The real question here is: how do we meaningfully define and edge
# Maybe we only want to consider countries which are in both the "to" and "from" categories
trade.plus <- trade[trade$to %in% trade$from, ] %>% # restrict countries
  group_by(from) %>% 
  mutate(amount_total_from = sum(amount)) %>% 
  mutate(prop_from = amount/amount_total_from ) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  mutate(amount_total_to = sum(amount)) %>% 
  mutate(prop_to = amount/amount_total_to)

# Add a new useful variable
trade.plus$country1 <- apply(trade.plus, 1, function(x) sort(unlist(c(x[1], x[2])))[1])  
trade.plus$country2 <- apply(trade.plus, 1, function(x) sort(unlist(c(x[1], x[2])))[2])  

# Let's come up with a meaninful definintion of an edge
# And the weight

trade.subset <- trade.plus %>% 
                filter(prop_from > 0.01 & prop_to > 0.01) %>%
                group_by(country1, country2) %>% 
                mutate(trade_volume = sum(amount)) %>% 
                select(country1, country2, trade_volume) %>% 
                unique()

nrow(trade.subset) # much more sparse!

# Preserve the node-specific vars above in the vertex df
trade.countries.subset <- trade.countries %>% 
        filter(code %in% union(unique(trade.subset$country1 ), unique(trade.subset$country2))) %>%
        inner_join(., trade.plus, by = c("code" = "from")) %>%
        select(code, name, amount_total_from) %>%
        unique() %>% 
        mutate(region = ifelse(code == "AUS" | code == "NZL", "Oceania", 
                        ifelse(code == "CHL" , "South America", 
                        ifelse(code == "ISR" | code == "TUR", "Middle East",
                        ifelse(code == "USA" | code == "MEX" | code == "CAN", "North America",
                        ifelse(code == "JPN", "Asia", "Europe"))))))

### A. Edge list to adjacency matrix
# In the real world networks are usually large and sparse and stored as edge lists

## Scenario 1: Directed, weighted, with annotations/meta data
g <- graph_from_data_frame(trade.subset, directed = FALSE, vertices = trade.countries.subset)

par(mfrow = c(2,2))
plot(g, main = "Default") # this looks terrible
plot(g, layout = layout.circle, main = "Circle")
plot(g, layout = layout.sphere, main = "Sphere")
plot(g, layout = layout.fruchterman.reingold, main = "Fruchterman Reingold")
par(mfrow = c(1,1))

# Nice to add edge weight = total volume of trade, and node size = total amount_from
hist(trade.subset$trade_volume)
hist(log(trade.subset$trade_volume))

# Add some fancy stuff
E(g)$weight <- log(trade.subset$trade_volume) # scale volume, as.numeric removes the scaling attribute
E(g)$weight <- as.numeric(scale(trade.subset$trade_volume, center = FALSE)) # scale volume, as.numeric removes the scaling attribute
#V(g)$size <- as.numeric(scale(trade.countries.subset$amount_total_from, center = FALSE))
V(g)$size <- log(trade.countries.subset$amount_total_from)# scale exports, as.numeric removes the scaling attribute
V(g)$name <- trade.countries.subset$code
V(g)$color <- as.factor(trade.countries.subset$region)
plot(g, edge.width = E(g)$weight/5, vertex.size = V(g)$size/5) # nodes too small, too much text, overlapping

# igraph objects can be passed to ggraph calls, which have better control of aesthetics 
ggraph(g, layout = 'stress') + 
  geom_edge_link(aes(width = weight), alpha = 0.25) + # Set edge weight with attributes of g
  scale_edge_width(range = c(0.5, 2.5)) +  # Constrain edge width
  geom_node_point(aes(size = size, color = color), alpha = 0.5) + # Set node color with attributes of g
  scale_size(range = c(2,10)) + # Constrain node sizes
  geom_node_text(aes(label = name), repel = TRUE, point.padding = unit(0.5, "lines")) + 
  theme_void() + 
  theme(legend.position = "none")

#----------------------------- TOPOLOGICAL SUMMARIES ----------------------------#

# See good details here: https://fukamilab.github.io/BIO202/09-B-networks.html

# We'll work from the adjacency matrix
# There are built in commands for a lot of these, but let's work from scratch
A.dgC <- as_adjacency_matrix(g) # Gives us a dgC matrix
A <- as.matrix(A.dgC)
nC <- nrow(trade.countries.subset)

## Vertex metrics

# Degree: most basic property, number of edges adjacent to a vertex
rowSums(A)

# Centrality
# Explain why eigen centrality (v) is proportional to the eigen centrality of its neighbors
eigen_centrality(g) # uses edge weights
eigen(A) # just uses the binary adjacency matrix
plot(eigen(A)$values)

# Specialization


## Network metrics
# Connectance: number of observed interactions divided by the number of possible interactions
  # Measure of complexity

sum(A[lower.tri(A)])/(nC*(nC -1)/2)

# Nestedness, modularity, specialization
  
## Use null models to determine if observed network structure is significantly different from random network structure


## More visualizations
df.an <- trade.countries.subset %>% # annotation df must have same rownames as A
         remove_rownames() %>%
         column_to_rownames(var = "code") %>%
         select(region)
  
pheatmap(A)
pheatmap(A, cluster_rows = T, cluster_cols = T,
         color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30), 
         main = "Observed Adjacency Matrix", 
         fontsize = 8, show_rownames = T, show_colnames = F, 
         legend = F, annotation_col = df.an)


# Look at edge weights: log will work better for this
E(g)$weight <- log(trade.subset$trade_volume) # scale volume, as.numeric removes the scaling attribute
A.weighted <- as_adjacency_matrix(g, attr = "weight")
pheatmap(A.weighted, cluster_rows = T, cluster_cols = T,
         color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30), 
         main = "Observed Adjacency Matrix", 
         fontsize = 8, show_rownames = T, show_colnames = F, 
         legend = T, annotation_col = df.an)

#Plot continents vs eigenvectors
# This looks kind of like a ring, evecs give rotation...
evecs <- eigen(A)$vectors
data.frame(evec1 = evecs[,1], evec2 = evecs[,2], region = as.factor(trade.countries.subset$region)) %>% 
  ggplot(aes(x = evec1, y = evec2, color = region), alpha = 0.25) + 
  geom_point() + 
  theme_minimal()

#---------------------   LOGISTIC REGRESSION -----------------------------------#


# Basic example
## define the model
glmmCode <- nimbleCode({
  beta0 ~ dnorm(0, sd = 10000)
  beta1 ~ dnorm(0, sd = 10000)
  sigma_RE ~ dunif(0, 1000)
  for (i in 1:N) {
    beta2[i] ~ dnorm(0, sd = sigma_RE)
    logit(p[i]) <- beta0 + beta1 * x[i] + beta2[i]
    r[i] ~ dbin(p[i], n[i])
  }
})

## constants, data, and initial values
glmmConsts <- list(N = 10)

glmmData <- list(
  r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46),
  n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79),
  x = c(0,  0,  0,  0,  0,  1, 1,  1,  1,  1)
)

glmmInits <- list(beta0 = 0, beta1 = 0, sigma_RE = 1)

## create the model object
glmmModel <- nimbleModel(code = code, constants = constants, data = data, 
                         inits = inits, check = FALSE)

# A. Easiest way to do this: 
niter <- 1000
nchains <- 2
mcmc.out <- nimbleMCMC(code = glmmCode, constants = glmmConsts,
                       data = glmmData, inits = glmmInits,
                       nchains = nchains, niter = niter,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('beta0', 'beta1'))


### Our data
# What's a reasonable model? Edge covariates: 1(same region), maybe random effects for each country


## define the model
glmCode <- nimbleCode({
  # Priors
  beta0 ~ dnorm(0, sd = 10000)
  beta1 ~ dnorm(0, sd = 10000)
  sigma_RE ~ dunif(0, 1000)
  
  # Model
  for(i in 1:N){ # Random effects for each vertex
    beta2[i] ~ dnorm(0, sd = sigma_RE)
  }
  for (i in 2:N) { # Likelihood 
    for (j in 1:(i-1)){
      logit(p[i,j]) <- beta0 + beta1 * x[i,j] + beta2[i] + beta2[j]
      y[i,j] ~ dbin(size = 1, prob = p[i,j])
    }
  }
})


## constants, data, and initial values
glmConsts <- list(N = nC)


# Prepare data: easier to do in R than in nimble model definition
reg <- as.factor(trade.countries.subset$region)
glmData <- list(
  y = A,
  x = outer(reg,reg, FUN = "==")*1 # x_ij = 1(region_i == region_j)
)

glmInits <- list(beta0 = 0, beta1 = 0, beta2 = rep(0,nC))

## create the model object
glmModel <- nimbleModel(code = glmCode, constants = glmConsts, data = glmData, 
                         inits = glmInits)

# A. Easiest way to do this: 
niter <- 1000
nchains <- 2
mcmc.out <- nimbleMCMC(code = glmCode, constants = glmConsts,
                       data = glmData, inits = glmInits,
                       nchains = nchains, niter = niter,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('beta0', 'beta1', 'beta2'))

mcmc.out$summary

# Questions: interpret coefficients, WAIC
# Add random effects for each country
