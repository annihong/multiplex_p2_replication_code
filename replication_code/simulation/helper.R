# given a vector of dyad-wise results d in {1,2,3,4}, and the number of actors n, reconstruct the adj matrix
dyads_to_matrix_1d <- function(y, n) {
  res <- matrix(rep(0, n*n), nrow = n)
  counter = 1
  for (i in 1:n) {
    for ( j in i:n) {
      if (i != j) {
        d = y[counter]
        if (d == 2) {
          res[i,j] = 1
        } else if (d == 3) {
          res[j,i] = 1
        } else if (d == 4) {
          res[i,j] = 1
          res[j,i] = 1
        }
        counter = counter + 1
      }
    }
  }
  return(res)
}

#given y_nd_ij and the total number of networks
#return a list of t y_1d_ij_t \in {1,2,3,4}
y_nd_ij_to_y_1d_ij <- function(y_nd_ij, total){
  one_t <- function(t){
    res = ceiling(y_nd_ij/4^(t - 1)) %% 4
    res = ifelse(res==0, 4,res) #replace 0 with 4
    return(res) 
  } 
  return(sapply(1:total, one_t))
}

#given a y \in {1,...,2^(2*t)}^N and t denoting total number of networks
#return a matrix of y = t x N,  y[i,] \in {1,2,3,4}^N
dyads_nd_to_1d <- function(y_nd, t){
  return(sapply(y_nd, y_nd_ij_to_y_1d_ij, total=t))
}

#given a y \in {1,...,2^(2*t)}^N, total number of networks t, and number of actors n
#return a list of t networks in adj form
dyads_to_matrix_nd <- function(y_nd, t, n){
  networks <- list()
  if (t < 2) {
    networks[[1]] <- dyads_to_matrix_1d(y_nd,n)
  } else {
    y_1d <- dyads_nd_to_1d(y_nd,t)
    for (i in 1:t) {
      networks[[i]] <- dyads_to_matrix_1d(y_1d[i,],n)
    }
  }
  return(networks)
}

#given a list of t network matrices, return a vector y of size N, such that y \in {1,...,2^(2*t)}^N representing the outcomes on a pair of dyads.
matrix_to_dyads_nd <- function(Ms) {
  #given a vector of dyad-wise outcome for dyads i
  #return the multiplex outcome y \in {1,...,2^(2*t)}
  helper <- function(ns){
    res <- sapply(2:t, function(i) (ns[i] - 1)*4^(i - 1))
    res = sum(res) + ns[1]
    return(res)
  }
  t = length(Ms)
  y_1d <- do.call(cbind, lapply(Ms, matrix_to_dyads_1d))
  y = apply(y_1d, 1, helper)
  return(y)
}


#given a matrix M, return the dyad-wise outcome y \in {1,2,3,4}^N
matrix_to_dyads_1d <- function(M){
  n = nrow(M)
  y <- c()
  for (i in 1:n){
    for (j in i:n){
      i_to_j = M[i,j]
      j_to_i = M[j,i]
      if (i != j) {
        outcome = 1
        if (M[i,j] & ! M[j,i]) {
          outcome = 2
        } else if (!M[i,j] & M[j,i]) {
          outcome = 3
        } else if (M[i,j] & M[j,i]) {
          outcome =4
        }
        y <- c(y, outcome)
      }
    }
  }
  return(y)
}



#helper function to zero out the diag of any squarematrix
zero_out_diag <- function(M){
  n = nrow(M)
  for (i in 1:n){
    M[i,i] = 0
  }
  return(M)
}

#given a stan object and the list of parameters of interests (in regex form)
# return the model summary relevant only to the parameters given
#CI = T: give the 95% credible interval, F: only the mean estimates
extract_model_estimates <- function(stanfit, params=c("^mu", "^rho","Sigma", "cross_mu", "cross_rho"), CI=F){
  sum <- summary(stanfit)$summary
  res <- lapply(params, function(x) sum[grep(x,rownames(sum)),1])
  names(res) <- params
  if (CI) {
    res <- lapply(params, function(x) sum[grep(x,rownames(sum)),c(1,4,8)])
    res <- data.frame(do.call(rbind, res))
  }
  return(res)
}

# given a B x N dataframe where each row is a network in dyad form (1 x N) with B total such networks, t representing the number of layers in the networks, and n number of actors per network
# return a list of B list of t adj matrices representing the B multiplex networks
dyads_to_matrix_list <- function(dyad_df, n, t, network_type = "adj") {
  B <- nrow(dyad_df)
  matrices <- list()
  if (t > 1) {
    for (i in 1:B) {
      res <- dyads_to_matrix_nd(dyad_df[i,], t, n)
      if (network_type == "igraph") {
        matrices[[i]] <- lapply(res, igraph::graph_from_adjacency_matrix)
      } else if (network_type == "network") {
        matrices[[i]] <- lapply(res, network, directed=T)
      } else if(network_type == "adj") {
        matrices[[i]] <- res      
      }
    }
  } else {
    #back to the uniplex case
    for (i in 1:B) {
      res <- dyads_to_matrix_1d(dyad_df[i,], n)
      if (network_type == "igraph") {
        matrices[[i]] <- igraph::graph_from_adjacency_matrix(res)
      } else if (network_type == "network") {
        matrices[[i]] <- network(res, directed=T)
      } else if(network_type == "adj") {
        matrices[[i]] <- res      
      }
    }
  }
  return(matrices)
}

get_dyads <- function(n){
  res = list()
  counter = 1
  for (i in 1:n) {
    for (j in i:n) {
      if (i != j) {
        res[[counter]] = c(i,j)
        counter = counter + 1
      }
    }
  }
  if (n < 2) {
    res = list(c(0,0))
  }
  return(res)
}

#missing value checks


