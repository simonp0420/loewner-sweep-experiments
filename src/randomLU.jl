function random_nonsingular_matrix(n::Integer, T::DataType = ComplexF64)
    mat = rand(T, n, n)
    mat = LowerTriangular(mat) * UpperTriangular(mat)
    return mat
end

function loewner_random_model(L, σL, V, W)
    n = size(L, 1)
    Wa = random_nonsingular_matrix(n)
    Wb = random_nonsingular_matrix(n)
    Wb .*= 2.0
    E = -(Wb * L  * Wa)
    A = -(Wb * σL * Wa)
    B =   Wb * V
    C =   W  * Wa
    return E, A, B, C
end

eval_random_model!(H, M, X, E, A, B, C, s) = eval_reduced_model!(H, M, X, E, A, B, C, s)
