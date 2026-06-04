using GLMakie

@inline function sparam(solver::F, i::Integer, j::Integer, f::Real) where F
    solver(f)[i, j]
end

function sparam(solver::F, i::Integer, j::Integer, f::AbstractVector) where F
    map(fᵢ -> sparam(solver, i, j, fᵢ), f)
end

let
    opts = SweepOptions(
        Δf=1e-5,
        q1=8,
        q2=12,
        adaptive=true,
        tol=1e-3,      # 10^(-60/20)
        l=0.25,        # quarter-wave transformer
        p=1,           # ZQT es SISO
        memory=3,
        use_D=false,
        data_partition = true
    )

    function Z_qwt(f; k=1 / 4)
        s = 2π * im * f
        70.7 * (170.7 * exp(s * k) + 29.3 * exp(-s * k)) /
        (170.7 * exp(s * k) - 29.3 * exp(-s * k))
    end

    S11_from_Zin(Zin; Z0=50.0) = (Zin - Z0) / (Zin + Z0)

    function solver_qwt_s11(f; Z0=50.0)
        Zin = Z_qwt(f)
        S11 = S11_from_Zin(Zin; Z0)
        return ComplexF64[S11;;]  # matriz 1×1
    end

    solver_qwt(f) = ComplexF64[Z_qwt(f);;]  # matriz 1×1
    fmin = 0.1
    fmax = 2.0
    fgrid = range(fmin, fmax, length=201)
    ITP = sweep(solver_qwt_s11, fgrid, opts);
    fgrid = range(fmin, fmax, length=501)
    ITP = sweep(solver_qwt_s11, fgrid, opts);
    fint = range(fmin, fmax, length=1201)
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, fint, 20*log10.(abs.(sparam(ITP, 1, 1, fint))))
    lines!(ax, fint, 20*log10.(abs.(sparam(solver_qwt_s11, 1, 1, fint))), linestyle = :dash)
    scatter!(ax, ITP.f₀, 20*log10.(abs.(sparam(ITP, 1, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")
    fig
end

let
    function solver_butterworth(f; n_order=5, fc=1.0)
        s = 2π * im * f / (2π * fc)  # frecuencia normalizada
        # polos de Butterworth orden n en el semiplano izquierdo
        poles = [exp(im * π * (2k + n_order - 1) / (2n_order)) 
                for k in 1:n_order]
        # H(s) = 1 / prod(s - pk)
        H = 1.0 / prod(s - p for p in poles)
        return ComplexF64[H;;]
    end

    opts = SweepOptions(
        Δf=1e-5,
        q1=8,
        q2=12,
        adaptive=true,
        tol=1e-3,
        l=0.25,
        memory=3,
        use_D=false,
        data_partition = true
    )

    fmin = 0.1
    fmax = 3.0
    fgrid = range(fmin, fmax, length=501)
    ITP = sweep(solver_butterworth, fgrid, opts);
    fint = range(fmin, fmax, length=1201)
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, fint, 20*log10.(abs.(sparam(ITP, 1, 1, fint))))
    lines!(ax, fint, 20*log10.(abs.(sparam(solver_butterworth, 1, 1, fint))), linestyle = :dash)
    scatter!(ax, ITP.f₀, 20*log10.(abs.(sparam(ITP, 1, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")
    fig
end

let
    opts = SweepOptions(
        Δf=1e-5,
        q1=8,
        q2=12,
        adaptive=true,
        tol=1e-3,
        l=0.25,
        memory=3,
        use_D=false,
        data_partition = true
    )

    function solver_bandpass(f; n_order=3, f0=1.0, BW=0.5)
        ω  = 2π * f
        ω0 = 2π * f0
        B  = 2π * BW
        # transformación paso-bajo → paso-banda
        s  = 2π * im * f
        sb = (s^2 + ω0^2) / (B * s)   # frecuencia normalizada paso-banda
        # polos de Butterworth orden n
        poles = [exp(im * π * (2*k + n_order - 1) / (2*n_order)) 
                 for k in 1:n_order]
        H = 1.0 / prod(sb - p for p in poles)
        return ComplexF64[H;;]
    end

    fmin = 0.1
    fmax = 3.0
    fgrid = range(fmin, fmax, length=101)
    ITP = sweep(solver_bandpass, fgrid, opts);
    fint = range(fmin, fmax, length=1201)
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, fint, 20*log10.(abs.(sparam(ITP, 1, 1, fint))))
    lines!(ax, fint, 20*log10.(abs.(sparam(solver_bandpass, 1, 1, fint))), linestyle = :dash)
    scatter!(ax, ITP.f₀, 20*log10.(abs.(sparam(ITP, 1, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")
    fig
end