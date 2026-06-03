using GLMakie

@inline function sparam(solver::F, i::Integer, j::Integer, f::Real) where F
    solver(f)[i, j]
end

function sparam(solver::F, i::Integer, j::Integer, f::AbstractVector) where F
    map(fᵢ -> sparam(solver, i, j, fᵢ), f)
end

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

function solver_nearresonance(f; f0=1.0, Q=30.0)
    s  = 2π * im * f
    s0 = 2π * im * f0
    H  = s0 / (s - s0 + s0/Q)
    return ComplexF64[H;;]
end

fmin = 0.1
fmax = 2.0
fgrid = range(fmin, fmax, length=501)
ITP = sweep(solver_nearresonance, fgrid, opts);
fint = range(fmin, fmax, length=1201)
lines(fint, 20*log10.(abs.(sparam(ITP, 1, 1, fint))))
scatter!(ITP.f₀, 20*log10.(abs.(sparam(ITP, 1, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")


function solver_butterworth(f; n_order=5, fc=1.0)
    s = 2π * im * f / (2π * fc)  # frecuencia normalizada
    # polos de Butterworth orden n en el semiplano izquierdo
    poles = [exp(im * π * (2k + n_order - 1) / (2n_order)) 
             for k in 1:n_order]
    # H(s) = 1 / prod(s - pk)
    H = 1.0 / prod(s - p for p in poles)
    return ComplexF64[H;;]
end

fmin = 0.1
fmax = 3.0
fgrid = range(fmin, fmax, length=501)
ITP = sweep(solver_butterworth, fgrid, opts);
fint = range(fmin, fmax, length=1201)
lines(fint, 20*log10.(abs.(sparam(ITP, 1, 1, fint))))
scatter!(ITP.f₀, 20*log10.(abs.(sparam(ITP, 1, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")


function solver_bandpass(f; n_order=3, f0=1.0, BW=0.5)
    ω  = 2π * f
    ω0 = 2π * f0
    B  = 2π * BW
    # transformación paso-bajo → paso-banda
    s  = 2π * im * f
    sb = (s^2 + ω0^2) / (B * s)   # frecuencia normalizada paso-banda
    # polos de Butterworth orden n
    poles = [exp(im * π * (2k + n_order - 1) / (2n_order)) 
             for k in 1:n_order]
    H = 1.0 / prod(sb - p for p in poles)
    return ComplexF64[H;;]
end

fmin = 0.1
fmax = 3.0
fgrid = range(fmin, fmax, length=101)
ITP = sweep(solver_bandpass, fgrid, opts);
fint = range(fmin, fmax, length=1201)
lines(fint, 20*log10.(abs.(sparam(ITP, 1, 1, fint))))
lines!(fint, 20*log10.(abs.(sparam(solver_bandpass, 1, 1, fint))))

scatter!(ITP.f₀, 20*log10.(abs.(sparam(ITP, 1, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")


function solver_coupler(f; f0=1.0, BW=0.5, Z0=50.0)
    ω  = 2π * f
    ω0 = 2π * f0
    B  = 2π * BW
    s  = 2π * im * f
    sb = (s^2 + ω0^2) / (B * s)

    # paso-bajo Butterworth orden 3 transformado a paso-banda
    poles = [exp(im * π * (2k + 2) / 6) for k in 1:3]
    H = 1.0 / prod(sb - p for p in poles)

    # matriz S de acoplador: S11≈0 en banda, S21≈H, S12=S21, S22=S11
    S21 = H / (1 + abs(H))
    S11 = sqrt(max(0.0, 1.0 - abs2(S21))) * exp(im * π * f / f0)
    return ComplexF64[S11 S21; S21 S11]
end

begin
opts = SweepOptions(
    Δf=1e-5,
    q1=8,
    q2=12,
    adaptive=true,
    p=2,
    memory=3,
    use_D = false,
    data_partition = true,
    parallel = true
)

fmin = 0.1
fmax = 2.5
fgrid = range(fmin, fmax, length=201)
ITP = sweep(solver_coupler, fgrid, opts);
fint = range(fmin, fmax, length=1201)

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, fint, 20*log10.(abs.(sparam(ITP, 2, 1, fint))))
lines!(ax, fint, 20*log10.(abs.(sparam(solver_coupler, 2, 1, fint))))
scatter!(ax, ITP.f₀, 20*log10.(abs.(sparam(ITP, 2, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")
fig
end


ITP = sweep(solver_coupler, fgrid, opts);
ITP2 = sweep_random(solver_coupler, fgrid, opts);

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, fint, 20*log10.(abs.(sparam(ITP, 2, 1, fint))))
lines!(ax, fint, 20*log10.(abs.(sparam(solver_coupler, 2, 1, fint))))
scatter!(ax, ITP.f₀, 20*log10.(abs.(sparam(ITP, 2, 1, (ITP.f₀)))), markersize=8, color=:red, label="muestras")
fig

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, fint, 20*log10.(abs.(sparam(ITP2, 2, 1, fint))))
lines!(ax, fint, 20*log10.(abs.(sparam(solver_coupler, 2, 1, fint))))
scatter!(ax, ITP2.f₀, 20*log10.(abs.(sparam(ITP2, 2, 1, (ITP2.f₀)))), markersize=8, color=:red, label="muestras")
fig


@btime sweep($solver_coupler, $fgrid, $opts);
@btime sweep_random($solver_coupler, $fgrid, $opts);