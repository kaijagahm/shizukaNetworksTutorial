# Dai Shizuka's networks tutorial: following along.
# The code is taken directly (in most cases) from https://dshizuka.github.io/networkanalysis/tutorials.html.

# Part 1: Introduction to Networks ----------------------------------------
library(igraph) # load the igraph package

g <- make_graph(~A-B-C-A, D-E-F-D, A-F) # won't use this function very often because usually we will be making graphs from data
plot(g)
class(g) # "igraph"
g # UN = "undirected" network with "names" of vertices. 6 vertices, 7 edges

# Look up the vertices
V(g) # listing the vertices of the graph

# Look up the edges
E(g) # listing the edges of the graph. Why is it shown as two dashes instead of just one?

#In addition, we can look up vertex attributes. We currently only have one vertex attribute, called ‘name’:
  
V(g)$name

# We can also create new vertex attributes using this syntax. Certain attribute names can be directly interpreted by igraph for example, the vertex attribute ‘color’ will automatically be interpreted for plotting the network. Let’s try this out:
  
V(g)$color <- c("white", "red", "green", "blue", "orange", "yellow") # a random set of colors 
plot(g) # that's cool! "color" always works for colors.

#We can also add edge attributes. Let’s try adding two edge attributes, weight and color.
E(g)$weight <- 1:7
E(g)$color <- rainbow(7) #rainbow() function chooses a specified number of colors 
plot(g) # okay wow now this is super complicated and impractical

# Part 2: Data Formats for Networks ---------------------------------------
# Adjacency Matrix
## extract the adjacency matrix:
am <- as_adjacency_matrix(g, sparse = F) # if we set sparse to T, then it would replace the 0's with periods

# Edge List
el <- as_edgelist(g) # for a directed network, the edge goes from the first column to the second column. If it's weighted, you'd add another column with the weights.

# Affiliation matrix (for when the network is constructed based on individuals' membership in groups)
A <- c(1,1,0,0) 
B <- c(1,0,1,0) 
C <- c(1,0,1,0) 
D <- c(0,1,0,1) 
E <- c(0,0,1,1) 
aff <- matrix(c(A,B,C,D,E), nrow = 5, byrow = TRUE) 
dimnames(aff) <- list(c("A","B","C","D","E"),
                      c("Group1","Group2","Group3","Group4"))
aff #The individual-by-group matrix. 1's and 0's represent whether or not each individual (A, B, C, D, or E) is a member of each group.

# To make this into a network, we have to figure out _in how many groups each pair of individuals co-occurs_. To do this, we multiply the matrix by its transpose.

adj <- aff %*% t(aff) # %*% means matrix multiplication in R
# that gives us an adjacency matrix
# diagonal is how many groups the individual co-occurs with itself in (i.e. how many groups it is in) and off-diagonals are how many groups that pair of individuals co-occurs in.
adj

# Use that adjacency matrix to make a network:
g2 <- graph_from_adjacency_matrix(adj, "undirected", weighted = T, diag = F)
plot(g2, edge.width = E(g2)$weight)

# Adjacency List
# This is sort of a wide-format data structure. Probably really awkward to deal with. Each row has one node in the first column, and then subsequent columns contain identities of the other nodes it's connected to.
# We can make it easier to deal with by converting it to a list in R, where each list element is a focal node and that element contains a vector of all the other nodes it's connected to.
as_adj_list(g) # but I can already see how this will be difficult to work with, since it stores information in the list names, and that's often not conducive to working with e.g. the tidyverse.

# Data formats for directed graphs
# There are some additional considerations when formatting data for directed graphs.

# Example graph:
dir.g <- make_graph(~A-+B-+C-+A, D-+E-+F-+D, A+-+F)
plot(dir.g)
as_adjacency_matrix(dir.g, sparse = F)
# in the adjacency matrix, we have a 1 if there is an edge going from the *row* to the *column*. So it's not symmetrical anymore. 

as_edgelist(dir.g) # for double-headed arrows, each edge will be listed separately, meaning that the edgelist will have more rows than it would in an undirected network (which would exclude one of the directions)

# Data formats for weighted networks
as_adjacency_matrix(g, sparse = F, attr = "weight")

# igraph allows us to display the graph as an edgelist instad.
as_data_frame(g)
# I like edgelists. Long format is good. They're tidy. Strongly prefer them to adjacency matrices.

# Going from Data to Networks
## importing an edge list
edge.dat <- read.csv("https://dshizuka.github.io/network2018/NetworkWorkshop_SampleData/sample_edgelist.csv") 
edge.dat

# can graph the edge list
set.seed(2)
eg <- graph_from_data_frame(edge.dat, directed = FALSE) 
eg
plot(eg, edge.width = E(eg)$weight)

## importing an adjacency matrix
am <- as.matrix(read.csv("https://dshizuka.github.io/network2018/NetworkWorkshop_SampleData/sample_adjmatrix.csv", 
                         header = T, row.names = 1)) # important to set the header and row names so that R knows what it's looking at.
am

# create a graph from the adjacency matrix. Have to specify that it's weighted.
g <- graph_from_adjacency_matrix(am, mode = "undirected", weighted = T)
plot(g, edge.width = E(g)$weight)

# Part 3: Plotting Basics -------------------------------------------------
# Let's start by making a graph out of an adjacency matrix
am <- as.matrix(read.csv("https://dshizuka.github.io/networkanalysis/SampleData/sample_adjmatrix.csv", header = T, row.names = 1))
am

g <- graph_from_adjacency_matrix(am, mode = "undirected", weighted = T)
plot(g, edge.width = E(g)$weight)

# option 1 (specify the layout in the plot call)
plot(g, layout = layout_in_circle(g))

# option 2 (specify the layout first, and then add it to the plot)
l<- layout_in_circle(g) 
plot(g, layout = l)
# the above two options give the same result.

# let's examine the layout object
class(l) # the layout is a matrix/array
l # coordinates for where to put the nodes
