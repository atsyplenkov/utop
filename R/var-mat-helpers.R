# Helpers for assembling semivariance matrices.
#
# These extract the two pieces of arithmetic that were copy-pasted across
# compute_var_mat_utop() (s7-var-mat.R) and compute_var_mat_default()
# (varMat.R): the nugget term added when params@nugget is set, and the
# regularisation that subtracts half the sum of diagonal semivariances from
# each off-diagonal entry.

# Nugget contribution for a pair of supports. `a1` and `a2` are area vectors
# (length n_obs and n_pred, or both n_obs for the obs-obs case); `overlap` is
# the n_obs x n_pred (or n_obs x n_obs) overlap matrix. Returns a matrix of
# the same shape as `overlap`, to be added to the semivariance matrix.
#' @noRd
nugget_matrix <- function(a1, a2, overlap, variogramModel) {
  f <- matrix(rep(a1, length(a2)), ncol = length(a2))
  s <- t(matrix(rep(a2, length(a1)), ncol = length(a1)))
  ared <- (1 / f + 1 / s - 2 * overlap / (f * s)) / 2
  matrix(
    mapply(FUN = nuggEx, ared, MoreArgs = list(variogramModel = variogramModel)),
    ncol = length(a2)
  )
}

# Regularise a square semivariance matrix: subtract half the sum of the two
# diagonal entries from each off-diagonal entry, mirroring across the
# diagonal. Used for the obs-obs matrix (equal supports).
#' @noRd
regularize_symmetric <- function(varMatrix) {
  ndim <- nrow(varMatrix)
  v_diag <- diag(varMatrix)
  for (ia in seq_len(ndim - 1L)) {
    for (ib in (ia + 1L):ndim) {
      varMatrix[ia, ib] <- varMatrix[ia, ib] - 0.5 * (v_diag[ia] + v_diag[ib])
      varMatrix[ib, ia] <- varMatrix[ia, ib]
    }
  }
  varMatrix
}

# Regularise a (possibly rectangular) semivariance matrix against two
# diagonal vectors: subtract 0.5*(sub1[i] + sub2[j]) from each entry [i,j].
# Used for cross (obs-pred) matrices and for the diag-exclusion case in
# compute_var_mat_matrix().
#' @noRd
regularize_cross <- function(varMatrix, sub1, sub2) {
  for (ia in seq_len(nrow(varMatrix))) {
    for (ib in seq_len(ncol(varMatrix))) {
      varMatrix[ia, ib] <- varMatrix[ia, ib] - 0.5 * (sub1[ia] + sub2[ib])
    }
  }
  varMatrix
}

# Element-wise nugget contribution for parallel area vectors, used by
# objfunc() during variogram fitting. Returns a vector of the same length as
# the inputs, to be added to the modelled gamma.
#' @noRd
nugget_vector <- function(a1, a2, overlap, variogramModel) {
  ared <- (1 / a1 + 1 / a2 - 2 * overlap / (a1 * a2)) / 2
  mapply(FUN = nuggEx, ared, MoreArgs = list(variogramModel = variogramModel))
}
