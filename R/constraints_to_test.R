# Polynomial constraints and rank constraints to test

#' Constructs matrix M
#' 
#' @param p Integer, number of nodes.
#' @param k Integer, order of highest cumulants to include.
#' 
#' @return Matrix M containing the second moments in the first and the third moments in 
#' the remaining rows. Matrix is represented as a list of lists such that it can passed to the C++ functions.
#' Each outer list corresponds to one column.
#' 
#' @export
construct_M <- function(p, k = c(3, 4)) {
  column_indices <- gtools::combinations(p, 2, 1:p, repeats.allowed = TRUE)
  column_indices <- lapply(split(column_indices, seq(nrow(column_indices))), function(ind) unlist(ind, use.names = FALSE))
  
  M <- matrix(list(), nrow = 1 + p, ncol = length(column_indices))
  for (i in 1:nrow(M)) {
    for (j in 1:length((column_indices))) {
      if (i == 1) {
        M[[i, j]] <- list(coeff=1, mom=column_indices[[j]])
      } else {
        M[[i, j]] <- list(coeff=1, mom=c(column_indices[[j]], i - 1))
      }
    }
  }
  
  # For order four, one row for each tuple i_1<=i_2 needs to be appended to M
  if (k==4) {
    order_four_part <- matrix(list(), nrow = length(column_indices), ncol = length(column_indices))
    for (row in 1:length(column_indices)) {
      for (col in 1:length(column_indices)) {
        order_four_part[[row, col]] <- list(coeff=1, mom=c(column_indices[[row]], column_indices[[col]]))
      }
    }
    
    M <- rbind(M, order_four_part)
  }
  
  # Cast into a list of lists
  return(apply(M, 2, as.list))
}

# Generate all minors of a specific size of a matrix
#' 
#' @param A List of list representing a matrix, whose minors to compute.
#' @param size size of the minors
#' 
#' @return List \code{all_minors} where each entry corresponds to one polynomial
#' 
#' @export
calculate_minors <- function(A, size) {
  # Cast list of lists into matrix
  A <- as.matrix(do.call(cbind, A)) 
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
  for (i_perm in seq_len(length(perms))) { 
    perm <- perms[i_perm]
    moms <- matrix(0, nrow=n, ncol=3)
    for (i in seq_len(n)) {
      curr_moment <- A[[i, as.function(perm)(i)]]$mom
      moms[i, 1:length(curr_moment)] <- curr_moment
    }
    det[[length(det) + 1]] <- list(coef=permutations::sgn(perm), moms=moms) 
  }
  return(det)
}

# Constraints stemming from lower symmetric tensor rank:
# p = 2
# - Strassen inequality since U-stat tests for ineq <= 0.
# A polynomial is represented as list of lists, where each inner list represents one summand of the ploynomial.
# The entry "mom" is a matrix where each row corresponds to one factor and contains all the indices of the 
# respective moment.
strassen_ineq <- list(
  list(coef=-4, moms=matrix(c(1,1,2, 1,1,2, 1,1,2, 2,2,2), ncol=3, byrow=TRUE)),
  list(coef=6, moms=matrix(c(1,1,1, 1,1,2, 1,2,2, 2,2,2), ncol=3, byrow=TRUE)),
  list(coef=3, moms=matrix(c(1,1,2, 1,1,2, 1,2,2, 1,2,2), ncol=3, byrow=TRUE)),
  list(coef=-4, moms=matrix(c(1,1,1, 1,2,2, 1,2,2, 1,2,2), ncol=3, byrow=TRUE)),
  list(coef=-1, moms=matrix(c(1,1,1, 1,1,1, 2,2,2, 2,2,2), ncol=3, byrow=TRUE))
)


# p = 3
Aronhold_invariant <- list(
  list(coef=1, moms = matrix(c(1,1,1, 1,2,3, 2,2,2, 3,3,3), ncol = 3, byrow = TRUE)),
  list(coef=-1, moms = matrix(c(1,1,2, 1,1,3, 2,2,2, 3,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(c(1,1,1, 1,2,2, 2,2,3, 3,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(c(1,1,1, 1,3,3, 2,2,2, 2,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(c(1,1,1, 1,2,3, 2,2,3, 2,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(c(1,1,3, 1,2,3, 1,3,3, 2,2,2), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(c(1,1,2, 1,2,3, 1,2,2, 3,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(1,1,1, 1,2,2, rep(c(2,3,3), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(1,1,1, 1,3,3, rep(c(2,2,3), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(1,1,2, 2,2,2, rep(c(1,3,3), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(2,2,2, 2,3,3, rep(c(1,1,3), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(2,2,3, 3,3,3, rep(c(1,1,2), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(1,1,3, 3,3,3, rep(c(1,2,2), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(rep(c(1,2,3), 4), ncol = 3, byrow = TRUE)), 
  list(coef=-2, moms = matrix(c(rep(c(1,2,3), 2), 1,2,2, 1,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-2, moms = matrix(c(rep(c(1,2,3), 2), 1,1,2, 2,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-2, moms = matrix(c(rep(c(1,2,3), 2), 1,1,3, 2,2,3), ncol = 3, byrow = TRUE)), 
  list(coef=-3, moms = matrix(c(1,2,3, 1,1,2, 2,2,3, 1,3,3), ncol = 3, byrow = TRUE)), 
  list(coef=-3, moms = matrix(c(1,2,3, 1,1,3, 1,2,2, 2,3,3), ncol = 3, byrow = TRUE)),
  list(coef=-1, moms = matrix(c(rep(c(1,2,2), 2), rep(c(1,3,3), 2)), ncol = 3, byrow = TRUE)),
  list(coef=-1, moms = matrix(c(rep(c(2,3,3), 2), rep(c(1,1,2), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=-1, moms = matrix(c(rep(c(1,1,3), 2), rep(c(2,2,3), 2)), ncol = 3, byrow = TRUE)), 
  list(coef=1, moms = matrix(c(2,3,3, 1,1,2, 1,1,3, 2,2,3), ncol = 3, byrow = TRUE)),
  list(coef=1, moms = matrix(c(1,1,3, 2,2,3, 1,2,2, 1,3,3), ncol = 3, byrow = TRUE)),  
  list(coef=1, moms = matrix(c(1,2,2, 1,3,3, 2,3,3, 1,1,2), ncol = 3, byrow = TRUE))
)

# larger p
construct_A_T <- function(p) {
  # Return A_T in basis representation for image space basis e_1 x e_1, e_1 x e_2, ...
  # and def space basis e_1 x (e_1 ^ e_2), ..., e_1 x (e_1 ^ e_a), e_1 x (e_2 ^ e_3), ...
  
  calculate_entry_A_T <- function(e_im, e_def) {
    i_0 <- e_im[1]
    j_0 <- e_def[1]
    i_1_to_end <- e_im[-1]
    j_1_to_end <- e_def[-1]
    if (all(i_1_to_end %in% j_1_to_end)) {
      differing_elem <- setdiff(j_1_to_end, i_1_to_end)
      
      # sign in this case depends on first position where both vectors differ
      # as they are both sorted. E.g. if it's the penultimate elem -> sign = -1,
      # third last -> sign = 1 and so on
      sgn <- (-1) ^ (length(j_1_to_end) - which(j_1_to_end==differing_elem))
      return(list(coeff=sgn, mom=c(j_0, i_0, differing_elem)))
    } else {
      return(list(coeff=1, mom=numeric(0)))
    }
  }
  
  a <- floor((p-1)/2)
  
  basis_Lambda_a <- gtools::combinations(p, a, repeats.allowed=F)
  basis_im_space <- cbind(rep(1:p, each = nrow(basis_Lambda_a)), do.call(rbind, replicate(p, basis_Lambda_a, simplify=FALSE)))
  basis_im_space <- apply(basis_im_space, 1, unlist, simplify = FALSE)
  
  basis_Lambda_a_plus_1 <- gtools::combinations(p, a+1, repeats.allowed=F)
  basis_def_space <- cbind(rep(1:p, each = nrow(basis_Lambda_a_plus_1)), do.call(rbind, replicate(p, basis_Lambda_a_plus_1, simplify=FALSE)))
  basis_def_space <- apply(basis_def_space, 1, unlist, simplify = FALSE)
  
  A_T <- matrix(list(), length(basis_im_space), length(basis_def_space))
  for (i in seq_along(basis_im_space)) {
    for (j in seq_along(basis_def_space)) {
      A_T[[i, j]] <- calculate_entry_A_T(basis_im_space[[i]], basis_def_space[[j]])
    }
  }
  
  return(apply(A_T, 2, as.list)) 
}