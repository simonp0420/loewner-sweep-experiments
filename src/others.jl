# Eq. 16. Alternative
@inline function svd_matrix_pencil_robust(L, σL)
    return svd!([L; σL]), svd!([L σL])
end

# Eq 17. Alternative
function loewner_reduce_robust(L, σL, B, C, F1, F2, r)
    X_r = @view F1.V[:, 1:r]   # na*p × r
    Y_r = @view F2.U[:, 1:r]   # nb*p × r
    Er = -(Y_r' * L * X_r)
    Ar = -(Y_r' * σL * X_r)
    Br = Y_r' * B
    Cr = C * X_r
    return Er, Ar, Br, Cr
end