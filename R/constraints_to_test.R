# Polynomial constraints and rank constraints to test

#' Constructs matrix M
#' 
#' @param p Integer, number of nodes.
#' @param k Integer, order of highest cumulants to include.
#' 
#' @return Matrix M containing the second moments in the first and the third moments in 
#' the remaining rows. 
#' 
#' @export
construct_M <- function(p, k = c(3, 4)) {
  column_indices <- gtools::combinations(p, 2, 1:p, repeats.allowed = TRUE)
  column_indices <- lapply(split(column_indices, seq(nrow(column_indices))), function(ind) unlist(ind, use.names = FALSE))
  
  M <- matrix("", nrow = 1 + p, ncol = length(column_indices))
  for (i in 1:nrow(M)) {
    for (j in 1:length((column_indices))) {
      if (i == 1) {
        M[[i, j]] <- paste(c("m", column_indices[[j]]), collapse = "")
      } else {
        M[[i, j]] <- paste0(c("m", sort(c(column_indices[[j]], i-1))), collapse = "")
      }
    } 
  }
  
  # For order four, one row for each tuple i_1<=i_2 needs to be appended to M
  if (k==4) {
    order_four_part <- matrix("", nrow = length(column_indices), ncol = length(column_indices))
    for (row in 1:length(column_indices)) {
      for (col in 1:length(column_indices)) {
        order_four_part[[row, col]] <- paste0(c("c", sort(c(column_indices[[row]], column_indices[[col]]))), collapse = "")
      }
    }
    
    M <- rbind(M, order_four_part)
  }
  
  return(M)
}

# Generate all minors of a specific size of a matrix
#' 
#' @param A A matrix, whose minors to compute. The format of the entries of the matrix needs to be "m122"
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
      # A[[i, as.function(perm)(i)]] contains moment in format m122. Here, a vector is needed
      mom_indices <- strsplit(curr_moment, "")[[1]]
      # Remove "m"
      mom_indices <- mom_indices[2:length(mom_indices)]
      mom_indices <- as.integer(mom_indices)
      moms[[i]] <- mom_indices
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
  list(coef=-1, moms = list(c(1, 1, 1), c(1, 1, 1), c(2, 2, 2), c(2, 2, 2))), 
  list(coef=-4, moms = list(c(1, 1, 1), c(1, 2, 2), c(1, 2, 2), c(1, 2, 2))), 
  list(coef=-4, moms = list(c(2, 2, 2), c(1, 1, 2), c(1, 1, 2), c(1, 1, 2))), 
  list(coef=3, moms = list(c(1, 1, 2), c(1, 1, 2), c(1, 2, 2), c(1, 2, 2))), 
  list(coef=6, moms = list(c(1, 1, 1), c(1, 1, 2), c(1, 2, 2), c(2, 2, 2)))
)


# p = 3
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



# p = 2, l = 1
disc <- list(
  list(coef=1, moms = list(c(1, 1, 1, 1), c(1, 1, 1, 1), c(1, 1, 1, 1), c(2, 2, 2, 2), c(2, 2, 2, 2), c(2, 2, 2, 2))), 
  list(coef=-64, moms = list(c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 2, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=-27, moms = list(c(1, 1, 1, 1), c(1, 1, 1, 1), c(1, 2, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=-27, moms = list(c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(2, 2, 2, 2), c(2, 2, 2, 2))), 
  list(coef=-54, moms = list(c(1, 1, 1, 1), c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=-54, moms = list(c(2, 2, 2, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 1, 2, 2))), 
  list(coef=-18, moms = list(c(1, 1, 1, 1), c(1, 1, 1, 1), c(1, 1, 2, 2), c(1, 1, 2, 2), c(2, 2, 2, 2), c(2, 2, 2, 2))), 
  list(coef=36, moms = list(c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=81, moms = list(c(1, 1, 1, 1), c(2, 2, 2, 2), c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 1, 2, 2))), 
  list(coef=-12, moms = list(c(1, 1, 1, 2), c(1, 2, 2, 2), c(1, 1, 1, 1), c(1, 1, 1, 1), c(2, 2, 2, 2), c(2, 2, 2, 2))), 
  list(coef=-6, moms = list(c(1, 1, 1, 1), c(2, 2, 2, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=54, moms = list(c(1, 1, 1, 1), c(1, 1, 2, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(2, 2, 2, 2), c(2, 2, 2, 2))), 
  list(coef=54, moms = list(c(1, 1, 2, 2), c(2, 2, 2, 2), c(1, 1, 1, 1), c(1, 1, 1, 1), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=108, moms = list(c(1, 1, 1, 1), c(1, 1, 1, 2), c(1, 1, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2), c(1, 2, 2, 2))), 
  list(coef=108, moms = list(c(1, 1, 2, 2), c(1, 2, 2, 2), c(2, 2, 2, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 1, 2))), 
  list(coef=-180, moms = list(c(1, 1, 1, 1), c(1, 1, 1, 2), c(1, 2, 2, 2), c(2, 2, 2, 2), c(1, 1, 2, 2), c(1, 1, 2, 2)))
)

P <- list(
  list(coef=1, moms = list(c(1, 1, 1, 2), c(1, 1, 1, 2))), 
  list(coef=-1, moms = list(c(1, 1, 1, 1), c(1, 1, 2, 2)))
)

D <- list(
  list(coef=12, moms = list(c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 1, 2), c(1, 1, 1, 2))), 
  list(coef=-1, moms = list(c(2, 2, 2, 2), c(1, 1, 1, 1), c(1, 1, 1, 1), c(1, 1, 1, 1))), 
  list(coef=9, moms = list(c(1, 1, 1, 1), c(1, 1, 1, 1), c(1, 1, 2, 2), c(1, 1, 2, 2))), 
  list(coef=-24, moms = list(c(1, 1, 1, 1), c(1, 1, 2, 2), c(1, 1, 1, 2), c(1, 1, 1, 2))), 
  list(coef=4, moms = list(c(1, 1, 1, 2), c(1, 2, 2, 2), c(1, 1, 1, 1), c(1, 1, 1, 1)))
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
      return(paste0(c(sgn, " * m", sort(c(j_0, i_0, differing_elem))), collapse = ""))
    } else {
      return("0")
    }
  }
  
  a <- floor((p-1)/2)
  
  basis_Lambda_a <- gtools::combinations(p, a, repeats.allowed=F)
  basis_im_space <- cbind(rep(1:p, each = nrow(basis_Lambda_a)), do.call(rbind, replicate(p, basis_Lambda_a, simplify=FALSE)))
  basis_im_space <- apply(basis_im_space, 1, unlist, simplify = FALSE)
  
  basis_Lambda_a_plus_1 <- gtools::combinations(p, a+1, repeats.allowed=F)
  basis_def_space <- cbind(rep(1:p, each = nrow(basis_Lambda_a_plus_1)), do.call(rbind, replicate(p, basis_Lambda_a_plus_1, simplify=FALSE)))
  basis_def_space <- apply(basis_def_space, 1, unlist, simplify = FALSE)
  
  A_T <- matrix("", length(basis_im_space), length(basis_def_space))
  for (i in seq_along(basis_im_space)) {
    for (j in seq_along(basis_def_space)) {
      A_T[[i, j]] <- calculate_entry_A_T(basis_im_space[[i]], basis_def_space[[j]])
    }
  }
  
  return(A_T) 
}

# Returns the formula of a cumulant in terms of moments. E.g. for the input "c1134", it returns "m1134 - m11 * m34 - m13 * m14 - m14 * m13".
get_cumulant_formula <- function(cum_expr) {
  indices <- unlist(strsplit(substring(cum_expr, 2, 5), split=""))
  paste0(c("m", indices, " - m", indices[c(1,2)], " * m", indices[c(3,4)], " - m", indices[c(1,3)], " * m", indices[c(2,4)], " - m", indices[c(1,4)], " * m", indices[c(2,3)]), collapse = "")
}