# Polynomial constraints and rank constraints to test

#' Constructs matrix M
#' 
#' @param p Integer, number of nodes.
#' @param orders List of integers, orders of cumulants to include.
#' 
#' @return Matrix M containing the cumulants of the lowest order in the first and the reamaining cumulants in 
#' the remaining rows. 
#' 
#' @export
construct_M <- function(p, orders) {
  first_order <- orders[1]
  remaining_orders <- orders[2:length(orders)]
  # Columns are indexed by  moments of the first_order
  column_indices <- gtools::combinations(p, first_order, 1:p, repeats.allowed = TRUE)
  column_indices <- lapply(split(column_indices, seq(nrow(column_indices))), function(ind) unlist(ind, use.names = FALSE))
  # First row consists of moments of the first_order
  row_indices <- list(numeric(0))
  # Remaining orders are arranged underneath
  for (order in remaining_orders) {
    additional_row_indices <- gtools::combinations(p, order - first_order, 1:p, repeats.allowed = TRUE)
    additional_row_indices <- lapply(split(additional_row_indices, seq(nrow(additional_row_indices))), function(ind) unlist(ind, use.names = FALSE))
    row_indices <- append(row_indices, additional_row_indices)
  }
  
  M <- matrix("", nrow = length(row_indices), ncol = length(column_indices))
  for (i in 1:length(row_indices)) {
    for (j in 1:length(column_indices)) {
      indices <- sort(c(row_indices[[i]], column_indices[[j]]))
      M[[i, j]] <- get_cumulant_formula(indices)
    } 
  }
  
  return(M)
}

#' Constructs (k/2)th flattening matrix of a tensor in (\mathbb{R}^p)^{\otimes k}.
#' 
#' @param p Integer, number of nodes.
#' @param k Integer, order of tensor.
#' 
#' @return Matrix fl. 
#' 
#' @export
construct_flattening_matrix <- function(p, k) {
  if ((p != 2 || k != 6) && (p != 3 || k != 4)) stop("Not implemented for values other than p = 2 and k = 6, or p = 3 and k = 4") 
  
  s <- floor(k/2)
  indices <- gtools::combinations(p, s, 1:p, repeats.allowed = TRUE)
  indices <- lapply(split(indices, seq(nrow(indices))), function(ind) unlist(ind, use.names = FALSE))
  fl <- matrix(NA, length(indices), length(indices))
  
  for (i in 1:length(indices)) {
    for (j in 1:length(indices)) {
      fl[i, j] <- get_cumulant_formula(sort(c(indices[[i]], indices[[j]])))
    }
  }
  return(fl)
}

#' Calculate all minors of a specific size of a matrix
#' 
#' @param A A matrix, whose minors to compute. The format of the entries of the matrix needs to be "m1_2_2"
#' @param size size of the minors
#' 
#' @return List \code{all_minors} where each entry corresponds to one polynomial
#' 
#' @export
calculate_minors <- function(A, size) {
  # Cast list of lists into matrix
  col_subsets <- apply(gtools::combinations(ncol(A), size, v = 1:ncol(A), repeats.allowed=F), 1, unlist, simplify = FALSE)
  row_subsets <- apply(gtools::combinations(nrow(A), size, v = 1:nrow(A), repeats.allowed=F), 1, unlist, simplify = FALSE)
  all_minors <- list()
  
  for (row_subset in row_subsets) {
    for (col_subset in col_subsets) {
      all_minors[[length(all_minors) + 1]] <- calculate_determinant(A[row_subset, col_subset])
    }
  }
  
  return(all_minors)
}

calculate_determinant <- function(A) {
  det <- list()
  n <- nrow(A)
  perms <- permutations::allperms(n) 
  for (i_perm in 1:length(perms)) { 
    perm <- perms[i_perm]
    moms <- list()
    for (i in 1:n) {
      curr_moment <- A[[i, as.function(perm)(i)]]
      # A[[i, as.function(perm)(i)]] contains moment in format m1_2_2. Here, a vector is needed. 
      # Remove m
      curr_moment <- substr(curr_moment, 2, nchar((curr_moment)))
      # Remove commata
      mom_indices <- strsplit(curr_moment, "_")[[1]]
      mom_indices <- as.integer(mom_indices)
      moms[[i]] <- mom_indices
    }
    det[[length(det) + 1]] <- list(coef=permutations::sgn(perm), moms=moms)
  }
  return(det)
}


# Strassen inequality
# - Strassen inequality since U-stat tests for ineq <= 0.
# A polynomial is represented as list of lists, where each inner list represents one summand of the ploynomial.
# The entry "moms" is a list where each entry corresponds to one factor and contains all the indices of the 
# respective moment.
strassen_ineq <- list(
  list(coef=-1, moms = list(c(1, 1, 1), c(1, 1, 1), c(2, 2, 2), c(2, 2, 2))), 
  list(coef=-4, moms = list(c(1, 1, 1), c(1, 2, 2), c(1, 2, 2), c(1, 2, 2))), 
  list(coef=-4, moms = list(c(2, 2, 2), c(1, 1, 2), c(1, 1, 2), c(1, 1, 2))), 
  list(coef=3, moms = list(c(1, 1, 2), c(1, 1, 2), c(1, 2, 2), c(1, 2, 2))), 
  list(coef=6, moms = list(c(1, 1, 1), c(1, 1, 2), c(1, 2, 2), c(2, 2, 2)))
)

# Aronhold invariant
Aronhold_invariant <- list(
  list(coef=-1, moms = list(c(1, 2, 3), c(1, 2, 3), c(1, 2, 3), c(1, 2, 3))), 
  list(coef=-1, moms = list(c(1, 1, 2), c(1, 1, 2), c(2, 3, 3), c(2, 3, 3))), 
  list(coef=-1, moms = list(c(1, 1, 3), c(1, 1, 3), c(2, 2, 3), c(2, 2, 3))), 
  list(coef=-1, moms = list(c(1, 2, 2), c(1, 2, 2), c(1, 3, 3), c(1, 3, 3))), 
  list(coef=1, moms = list(c(1, 1, 1), c(1, 2, 2), c(2, 3, 3), c(2, 3, 3))), 
  list(coef=1, moms = list(c(1, 1, 1), c(1, 3, 3), c(2, 2, 3), c(2, 2, 3))), 
  list(coef=1, moms = list(c(1, 1, 2), c(2, 2, 2), c(1, 3, 3), c(1, 3, 3))), 
  list(coef=1, moms = list(c(1, 1, 3), c(3, 3, 3), c(1, 2, 2), c(1, 2, 2))), 
  list(coef=1, moms = list(c(2, 2, 2), c(2, 3, 3), c(1, 1, 3), c(1, 1, 3))), 
  list(coef=1, moms = list(c(2, 2, 3), c(3, 3, 3), c(1, 1, 2), c(1, 1, 2))), 
  list(coef=2, moms = list(c(1, 1, 2), c(2, 3, 3), c(1, 2, 3), c(1, 2, 3))), 
  list(coef=2, moms = list(c(1, 1, 3), c(2, 2, 3), c(1, 2, 3), c(1, 2, 3))), 
  list(coef=2, moms = list(c(1, 2, 2), c(1, 3, 3), c(1, 2, 3), c(1, 2, 3))), 
  list(coef=1, moms = list(c(1, 1, 1), c(1, 2, 3), c(2, 2, 2), c(3, 3, 3))), 
  list(coef=1, moms = list(c(1, 1, 2), c(1, 1, 3), c(2, 2, 3), c(2, 3, 3))), 
  list(coef=1, moms = list(c(1, 1, 2), c(1, 2, 2), c(1, 3, 3), c(2, 3, 3))), 
  list(coef=1, moms = list(c(1, 1, 3), c(1, 2, 2), c(1, 3, 3), c(2, 2, 3))), 
  list(coef=-1, moms = list(c(1, 1, 1), c(1, 2, 2), c(2, 2, 3), c(3, 3, 3))), 
  list(coef=-1, moms = list(c(1, 1, 1), c(1, 2, 3), c(2, 2, 3), c(2, 3, 3))), 
  list(coef=-1, moms = list(c(1, 1, 1), c(1, 3, 3), c(2, 2, 2), c(2, 3, 3))), 
  list(coef=-1, moms = list(c(1, 1, 2), c(1, 1, 3), c(2, 2, 2), c(3, 3, 3))), 
  list(coef=-1, moms = list(c(1, 1, 2), c(1, 2, 2), c(1, 2, 3), c(3, 3, 3))), 
  list(coef=-1, moms = list(c(1, 1, 3), c(1, 2, 3), c(1, 3, 3), c(2, 2, 2))), 
  list(coef=-3, moms = list(c(1, 1, 2), c(1, 2, 3), c(1, 3, 3), c(2, 2, 3))), 
  list(coef=-3, moms = list(c(1, 1, 3), c(1, 2, 2), c(1, 2, 3), c(2, 3, 3)))
)

#' Constructs Young flattening Y_T
#' 
#' @param p Integer, number of nodes.
#' @param k Integer, order of the tensor.
#' 
#' @return Young flattening Y_T as matrix.
#' 
#' @export
construct_Y_T <- function(p, k) {
  # Return A_T in basis representation for image space basis e_1 x e_1, e_1 x e_2, ...
  # and def space basis e_1 x (e_1 ^ e_2), ..., e_1 x (e_1 ^ e_a), e_1 x (e_2 ^ e_3), ...
  
  calculate_entry_Y_T <- function(e_im, e_def, k) {
    cutoff <- if (k==3) 1 else 2 
    i_0 <- e_im[1:cutoff]
    j_0 <- e_def[1:cutoff]
    i_1_to_end <- e_im[-(1:cutoff)]
    j_1_to_end <- e_def[-(1:cutoff)]
    if (all(i_1_to_end %in% j_1_to_end)) {
      differing_elem <- setdiff(j_1_to_end, i_1_to_end)
      
      # sign in this case depends on first position where both vectors differ
      # as they are both sorted. E.g. if it's the penultimate elem -> sign = -1,
      # third last -> sign = 1 and so on
      sgn <- (-1) ^ (length(j_1_to_end) - which(j_1_to_end==differing_elem))
      return(paste0(sgn, " * (", get_cumulant_formula(sort(c(j_0, i_0, differing_elem))), ")"))
    } else {
      return("0")
    }
  }
  
  a <- floor((p-1)/2)
  
  basis_Lambda_a <- apply(gtools::combinations(p, a, repeats.allowed=FALSE), 1, unlist, simplify = FALSE)
  basis_Lambda_a_plus_1 <- apply(gtools::combinations(p, a+1, repeats.allowed=FALSE), 1, unlist, simplify = FALSE)
  
  if (k==3) {
    basis_im_space <- apply(expand.grid(1:p, basis_Lambda_a), 1, function(row) unlist(row, use.names = FALSE), simplify = FALSE)
    basis_def_space <- apply(expand.grid(1:p, basis_Lambda_a_plus_1), 1, function(row) unlist(row, use.names = FALSE), simplify = FALSE)
  } else if (k==5) {
    basis_S_2 <- apply(gtools::combinations(p, 2, 1:p, repeats.allowed = TRUE), 1, unlist, simplify = FALSE)
    basis_im_space <- apply(expand.grid(basis_S_2, basis_Lambda_a), 1, function(row) unlist(row, use.names = FALSE), simplify = FALSE)
    basis_def_space <- apply(expand.grid(basis_S_2, basis_Lambda_a_plus_1), 1, function(row) unlist(row, use.names = FALSE), simplify = FALSE)
  } else {
    stop("Only implemented for k=3, 5")
  }

  Y_T <- matrix("", length(basis_im_space), length(basis_def_space))
  for (i in seq_along(basis_im_space)) {
    for (j in seq_along(basis_def_space)) {
      Y_T[[i, j]] <- calculate_entry_Y_T(basis_im_space[[i]], basis_def_space[[j]], k)
    }
  }
  
  return(Y_T) 
}

# Returns the formula of a cumulant with indices \code{ind} in terms of moments. E.g. for the input c(1,2,3,4), it returns "m1_2_3_4 - m1_2 * m3_4 - m1_3 * m2_4 - m1_4 * m2_3".
get_cumulant_formula <- function(ind) {
  k <- length(ind)
  if (k <= 3) {
    result <- paste0("m", paste0(ind, collapse="_"))
  }
  else if (k == 4) {
    result <- paste0("m", paste0(ind, collapse="_"), " - m", paste0(ind[c(1,2)], collapse="_"), " * m", paste0(ind[c(3,4)], collapse="_"), " - m", paste0(ind[c(1,3)], collapse="_"), " * m", paste0(ind[c(2,4)], collapse="_"), " - m", paste0(ind[c(1,4)], collapse="_"), " * m", paste0(ind[c(2,3)], collapse="_"))
  } else {
    # Need to sum over partitions 
    result <- "" 
    for (partition in partitions::listParts(k)) {
      # Only partitions where each set has at least two elements
      if (all(sapply(partition, length)>1)) {
        l <- length(partition)
        coeff <- (-1)^(l-1) * factorial(l-1)
        coeff <- if (coeff > 0) paste("+", coeff) else paste0("+ (", coeff, ")")
        moments <- sapply(partition, function(set) paste0(ind[set], collapse = "_"))
        result <- paste0(result, " ", coeff, paste0(" * m", moments, collapse = ""))
      }
    }
    # Remove leading " + "
    result <- substring(result, 4, nchar(result))
  }
  return(result)
} 