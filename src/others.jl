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

#=
function loewner_pseudo_errors!(errors, f, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, Δf)
    p  = size(C_r1, 1)
    r1 = size(E_r1, 1)
    r2 = size(E_r2, 1)

    sE_r1    = Matrix{ComplexF64}(undef, r1, r1)
    sE_r2    = Matrix{ComplexF64}(undef, r2, r2)
    sEB_r1   = Matrix{ComplexF64}(undef, r1, p)
    sEB_r2   = Matrix{ComplexF64}(undef, r2, p)
    H_r1_buf = Matrix{ComplexF64}(undef, p, p)
    H_r2_buf = Matrix{ComplexF64}(undef, p, p)
    diff_buf = Matrix{ComplexF64}(undef, p, p)

    #f₀_set = Set(f₀)

    for (idx, fᵢ) in enumerate(f)
        #if fᵢ ∈ f₀_set
        #    errors[idx] = -Inf
        #    continue
        #end

        if any(isapprox(fᵢ, fs; rtol=0, atol=1e-12) for fs in f₀)
            errors[idx] = -Inf
            continue
        end

        s = im * 2*π * fᵢ
        s_shift = im * 2*π * (fᵢ + Δf)

        @. sE_r1 = s       * E_r1 - A_r1
        @. sE_r2 = s_shift * E_r2 - A_r2

        ldiv!(sEB_r1, lu!(sE_r1), B_r1)   # lu! factoriza sE_r1 in-place
        ldiv!(sEB_r2, lu!(sE_r2), B_r2)   # ldiv! usa la factorización

        mul!(H_r1_buf, C_r1, sEB_r1)
        mul!(H_r2_buf, C_r2, sEB_r2)

        @. diff_buf = H_r2_buf - H_r1_buf

        errors[idx] = maximum(svdvals(diff_buf)) / maximum(svdvals(H_r1_buf))
    end
    return errors
end
=#