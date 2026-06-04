function loewner_random_model(L, σL, V, W)
    n = size(L, 1)
    Wa = rand(ComplexF64, n, n)
    Wb = 2 .* rand(ComplexF64, n, n)
    E = -(Wb * L  * Wa)
    A = -(Wb * σL * Wa)
    B =   Wb * V
    C =   W  * Wa
    return E, A, B, C
end

# qr! instead of lu!
function eval_random_model!(H, M, X, E, A, B, C, s)
    M .= @. s * E - A
    F = qr!(M)
    ldiv!(X, F, B)
    mul!(H, C, X)
    return H
end

function loewner_pseudo_errors_random!(errors, fgrid, f₀, 
                                        E1, A1, B1, C1,
                                        E2, A2, B2, C2)
    p  = size(C1, 1)
    r1 = size(E1, 1)
    r2 = size(E2, 1)
    f₀_set = Set(f₀)
    nthreads = Threads.nthreads()
    n = length(fgrid)
    chunk_size = ceil(Int, n / nthreads)

    tasks = map(1:nthreads) do t
        Threads.@spawn begin
            M1 = Matrix{ComplexF64}(undef, r1, r1)
            M2 = Matrix{ComplexF64}(undef, r2, r2)
            X1 = Matrix{ComplexF64}(undef, r1, p)
            X2 = Matrix{ComplexF64}(undef, r2, p)
            H1 = Matrix{ComplexF64}(undef, p, p)
            H2 = Matrix{ComplexF64}(undef, p, p)
            ΔH = Matrix{ComplexF64}(undef, p, p)

            chunk = (t-1)*chunk_size+1 : min(t*chunk_size, n)
            for idx in chunk
                fi = fgrid[idx]
                if fi ∈ f₀_set
                    errors[idx] = -Inf
                    continue
                end
                s = 2π * im * fi
                eval_random_model!(H1, M1, X1, E1, A1, B1, C1, s)
                eval_random_model!(H2, M2, X2, E2, A2, B2, C2, s)
                @. ΔH = H2 - H1
                denom = maximum(svdvals(H1))
                errors[idx] = denom == 0 ? Inf : maximum(svdvals(ΔH)) / denom
            end
        end
    end
    foreach(wait, tasks)
    return errors
end

function sweep_random(solve::Func, f, sweep_options) where Func

    d_flag = sweep_options.use_D
    p_flag = sweep_options.data_partition
    D = ones(sweep_options.p, sweep_options.p)

    fmin = minimum(f)
    fmax = maximum(f)
    n₀ = semi_adaptive_n(sweep_options, fmax)
    f₀ = Float64[]
    if sweep_options.adaptive
        append!(f₀, [fmin, fmax])
    else # semi-adaptive
        append!(f₀, semi_adaptive_freqs(fmin, fmax, max(2, n₀)))
    end

    H₀ = [solve(fᵢ) for fᵢ in f₀]
    if d_flag
        for Hᵢ in H₀
            Hᵢ .-= D
        end
    end

    sa, sb, Ha, Hb = loewner_data_partition(f₀, H₀, p_flag)
    L, σL = loewner_mfti(Ha, Hb, sa, sb)
    V, W = vw_mfti(Ha, Hb)

    errors = Vector{Float64}(undef, length(f))
    memory = 0
    while memory < sweep_options.memory
        E1, A1, B1, C1 = loewner_random_model(L, σL, V, W)
        E2, A2, B2, C2 = loewner_random_model(L, σL, V, W)
        loewner_pseudo_errors_random!(errors, f, f₀, E1, A1, B1, C1, E2, A2, B2, C2)
        f_new = f[argmax(errors)]
        
        H_new = solve(f_new)
        if d_flag
            H_new .-= D
        end

        s_new = im * 2 * π * f_new
        Hr_new = C1 * ((s_new * E1 - A1) \ B1) #+ D
        error_act = maximum(svdvals(Hr_new - H_new)) / maximum(svdvals(H_new))
        if error_act <= sweep_options.tol
            memory += 1
        else
            memory = 0
        end
        push!(f₀, f_new)
        push!(H₀, H_new)
        sort_idx = sortperm(f₀)
        permute!(f₀, sort_idx)
        permute!(H₀, sort_idx)
        sa, sb, Ha, Hb = loewner_data_partition(f₀, H₀, p_flag)
        L, σL = loewner_mfti(Ha, Hb, sa, sb)
        V, W = vw_mfti(Ha, Hb)
    end

    E1, A1, B1, C1 = loewner_random_model(L, σL, V, W)

    return LoewnerStateModel(E1, A1, B1, C1, f₀, H₀)
end