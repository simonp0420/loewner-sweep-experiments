using LinearAlgebra

@kwdef struct SweepOptions
    Δf::Float64 = 1e-5
    q1::Int = 8
    q2::Int = 12
    adaptive::Bool = true
    tol::Float64 = 10^(-60 / 20)
    l::Float64 = 0.1
    p::Int = 2
    memory::Int = 3
    use_D::Bool = true
    data_partition::Bool = true
    parallel::Bool = false
    nthreads::Int = Threads.nthreads()
end

struct LoewnerStateModel{T}
    E::Matrix{Complex{T}}
    A::Matrix{Complex{T}}
    B::Matrix{Complex{T}}
    C::Matrix{Complex{T}}
    f₀::Vector{T}
    H₀::Vector{Matrix{Complex{T}}}
end

function (interp::LoewnerStateModel)(f::Real)
    s = im * 2π * f
    return interp.C * ((s * interp.E - interp.A) \ interp.B)
end

function (interp::LoewnerStateModel)(f::AbstractVector)
    return map(fᵢ -> interp(fᵢ), f)
end

function loewner_partition_posneg(f, S)
    n = length(f)
    n == length(S) || throw(ArgumentError("f and S must have same length"))

    any(iszero, f) && throw(ArgumentError("f contains 0. Use fmin > 0 to avoid s == conj(s)."))

    sa = Vector{ComplexF64}(undef, n)
    sb = Vector{ComplexF64}(undef, n)
    Ha = Vector{Matrix{ComplexF64}}(undef, n)
    Hb = Vector{Matrix{ComplexF64}}(undef, n)

    for k in 1:n
        s = 2π * im * f[k]
        sb[k] = s
        Hb[k] = Matrix{ComplexF64}(S[k])
        sa[k] = conj(s)
        Ha[k] = conj(S[k])
    end

    return sa, sb, Ha, Hb
end


function loewner_partition_evenodd(f, S)
    n = length(f)
    n == length(S) || throw(ArgumentError("f and S must have same length"))
    any(iszero, f) && throw(ArgumentError("f contains 0. Use fmin > 0 to avoid s == conj(s)."))

    sa = Vector{ComplexF64}(undef, n)
    sb = Vector{ComplexF64}(undef, n)
    Ha = Vector{Matrix{ComplexF64}}(undef, n)
    Hb = Vector{Matrix{ComplexF64}}(undef, n)

    n_even = iseven(n) ? n : n - 1  

    for i in 1:2:n_even
        si = 2π * im * f[i]
        k = (i + 1) ÷ 2
        sb[2*k-1] = si
        sb[2*k]   = conj(si)
        Hb[2*k-1] = S[i]
        Hb[2*k]   = conj(S[i])
    end
    for i in 2:2:n_even
        si = 2π * im * f[i]
        k = i ÷ 2
        sa[2*k-1] = si
        sa[2*k]   = conj(si)
        Ha[2*k-1] = S[i]
        Ha[2*k]   = conj(S[i])
    end

    # último punto si n es impar: positivo-negativo
    if isodd(n)
        sn = 2π * im * f[n]
        sb[n] = sn
        Hb[n] = S[n]
        sa[n] = conj(sn)
        Ha[n] = conj(S[n])
    end

    return sa, sb, Ha, Hb
end

function loewner_data_partition(f, S, flag)
    if flag == true
        return loewner_partition_evenodd(f, S)
    else
        return loewner_partition_posneg(f, S)
    end
end

function loewner_mfti(Ha, Hb, sa::AbstractVector{Complex{T}}, sb::AbstractVector{Complex{T}}) where T
    p = size(Ha[1], 1)
    n = length(sa)
    L = Matrix{Complex{T}}(undef, n * p, n * p)
    σL = Matrix{Complex{T}}(undef, n * p, n * p)
    for i in eachindex(Hb)
        rows = (i-1)*p+1:i*p
        Hbi = Hb[i]
        sbi = sb[i]
        for j in eachindex(Ha)
            cols = (j-1)*p+1:j*p
            Haj = Ha[j]
            saj = sa[j]
            L[rows, cols] .= (Haj - Hbi) / (saj - sbi)
            σL[rows, cols] .= (saj * Haj - sbi * Hbi) / (saj - sbi)
        end
    end
    return L, σL
end


function vw_mfti(Ha, Hb)
    n = length(Ha)
    p = size(Ha[1], 1)
    W = Matrix{ComplexF64}(undef, p, n * p)
    V = Matrix{ComplexF64}(undef, n * p, p)
    for j in 1:n
        cols = (j-1)*p+1:j*p
        rows = (j-1)*p+1:j*p
        W[:, cols] .= Ha[j]
        V[rows, :] .= Hb[j]
    end
    return V, W
end

semi_adaptive_n(options, fmax) = ceil(Int, 15 * options.l * fmax / (options.p * 3E8))
semi_adaptive_freqs(fmin, fmax, n₀) = 2 * fmax .+ fmin .- 10.0 .^ (range(log10(2 * fmax), log10(fmax + fmin), n₀))

# Equation 19.
function loewner_order(σ::AbstractVector, q::Int=8)
    6 <= q <= 12 || throw(ArgumentError("loewner_order: q must satisfy 6 ≤ q ≤ 12"))
    total = sum(σ)
    total < eps() && return 1
    threshold = total * (1 - 10.0^(-q))
    accum = zero(eltype(σ))
    for (r, σv) in enumerate(σ)
        accum += σv
        accum >= threshold && return r
    end
    return length(σ)
end

# Eq. 16
@inline function svd_matrix_pencil(x, L, σL)
    return svd(x * L - σL)
end

# Eq. 17
function loewner_reduce(L, σL, B, C, F, r)
    Xr = @view F.V[:, 1:r] # n*p × r
    Yr = @view F.U[:, 1:r] # n*p × r
    Er = -(Yr' * L * Xr)
    Ar = -(Yr' * σL * Xr)
    Br = Yr' * B
    Cr = C * Xr
    return Er, Ar, Br, Cr
end


function eval_reduced_model!(H, M, X, E, A, B, C, s)
    M .= @. s * E - A
    F = lu!(M)
    ldiv!(X, F, B)
    mul!(H, C, X)
    return H
end


function loewner_pseudo_errors!(errors, f, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, Δf, sweep_options)
    if sweep_options.parallel == false
        return loewner_pseudo_errors_sequential!(errors, f, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, Δf)
    else
        return loewner_pseudo_errors_parallel!(errors, f, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, Δf, sweep_options.nthreads)
    end
end

function loewner_pseudo_errors_sequential!(errors, fgrid, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, Δf)
    p  = size(C_r1, 1)
    r1 = size(E_r1, 1)
    r2 = size(E_r2, 1)

    M1 = Matrix{ComplexF64}(undef, r1, r1)
    M2 = Matrix{ComplexF64}(undef, r2, r2)

    X1 = Matrix{ComplexF64}(undef, r1, p)
    X2 = Matrix{ComplexF64}(undef, r2, p)

    H1 = Matrix{ComplexF64}(undef, p, p)
    H2 = Matrix{ComplexF64}(undef, p, p)
    ΔH = Matrix{ComplexF64}(undef, p, p)

    f₀_set = Set(f₀)
    for idx in eachindex(fgrid)
        fi = fgrid[idx]
        if fi in f₀_set
            errors[idx] = -Inf
            continue
        end
        s1 = 2π * im * fi
        s2 = 2π * im * (fi + Δf)

        eval_reduced_model!(H1, M1, X1, E_r1, A_r1, B_r1, C_r1, s1)
        eval_reduced_model!(H2, M2, X2, E_r2, A_r2, B_r2, C_r2, s2)
        @. ΔH = H2 - H1
        denom = maximum(svdvals!(H1))
        errors[idx] = denom == 0 ? Inf : maximum(svdvals!(ΔH)) / denom
    end

    return errors
end

function loewner_pseudo_errors_parallel!(errors, fgrid, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, Δf, nthreads)
    p  = size(C_r1, 1)
    r1 = size(E_r1, 1)
    r2 = size(E_r2, 1)
    f₀_set = Set(f₀)
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
                s1 = 2*π * im * fi
                s2 = 2*π * im * (fi + Δf)
                eval_reduced_model!(H1, M1, X1, E_r1, A_r1, B_r1, C_r1, s1)
                eval_reduced_model!(H2, M2, X2, E_r2, A_r2, B_r2, C_r2, s2)
                @. ΔH = H2 - H1
                denom = maximum(svdvals(H1))
                errors[idx] = denom == 0 ? Inf : maximum(svdvals(ΔH)) / denom
            end
        end
    end
    foreach(wait, tasks)
    return errors
end


function sweep(solve::Func, f, sweep_options) where Func

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
        x = 2π * (f₀[1] + f₀[end]) / 3
        F = svd_matrix_pencil(x, L, σL)
        r1 = loewner_order(F.S, sweep_options.q1)
        E_r1, A_r1, B_r1, C_r1 = loewner_reduce(L, σL, V, W, F, r1)
        r2 = loewner_order(F.S, min(sweep_options.q1 + 4, 12))
        E_r2, A_r2, B_r2, C_r2 = loewner_reduce(L, σL, V, W, F, r2)

        loewner_pseudo_errors!(errors, f, f₀, E_r1, A_r1, B_r1, C_r1, E_r2, A_r2, B_r2, C_r2, sweep_options.Δf, sweep_options)
        f_new = f[argmax(errors)]
        
        H_new = solve(f_new)
        if d_flag
            H_new .-= D
        end

        s_new = im * 2 * π * f_new
        Hr_new = C_r1 * ((s_new * E_r1 - A_r1) \ B_r1) #+ D
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

    x = im * 2π * (f₀[1] + f₀[end]) / 3
    F = svd_matrix_pencil(x, L, σL)
    r1 = loewner_order(F.S, sweep_options.q1)
    E_r1, A_r1, B_r1, C_r1 = loewner_reduce(L, σL, V, W, F, r1)

    return LoewnerStateModel(E_r1, A_r1, B_r1, C_r1, f₀, H₀)
end
