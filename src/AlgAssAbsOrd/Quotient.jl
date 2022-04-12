################################################################################
#
#  Quotients of orders as orders
#
################################################################################

function quotient_order(O::AlgAssAbsOrd, I::AlgAssAbsOrdIdl)
  M = basis_matrix_wrt(I, O)
  @assert isone(denominator(M))
  S, U, V = snf_with_transform(numerator(M))
  adjusted_basis_matrix = inv(V)
  @show adjusted_basis_matrix
  adjusted_basis_all = [elem_from_mat_row(O, adjusted_basis_matrix, i) for i in 1:degree(O)]
  k = findfirst(iszero, diagonal(S))
  adjusted_basis = adjusted_basis_all[k:end]
  l = length(adjusted_basis)
  mt = Array{fmpq, 3}(undef, l, l, l)
  for i in 1:l
    for j in 1:l
      mt[i, j, :] = (coordinates(adjusted_basis[i] * adjusted_basis[j]) * V)[k:end]
    end
  end
  quoAlg = AlgAss(QQ, mt)
end
