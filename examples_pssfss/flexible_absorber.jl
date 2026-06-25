include(joinpath(@__DIR__, "common.jl"))

function flexible_absorber_case(; nfreq = 196, square_n = 21, disk_ntri = 800)
    p = 20
    a, Rsquare = 0.5, 120.0
    d, Rcircle = 9.0, 400.0

    squares = rectstrip(;
        Px = p,
        Py = p,
        Lx = p - 2a,
        Ly = p - 2a,
        Nx = square_n,
        Ny = square_n,
        units = mm,
        Zsheet = Rsquare,
    )
    disks = polyring(;
        units = mm,
        sides = 40,
        a = [0.0],
        b = [d],
        s1 = [p, 0],
        s2 = [0, p],
        ntri = disk_ntri,
        Zsheet = Rcircle,
    )

    strata = [
        Layer()
        Layer(width = 4mm, ϵᵣ = 2.9, tanδ = 0.1)
        disks
        Layer(width = 4mm, ϵᵣ = 2.9, tanδ = 0.1)
        squares
        Layer(width = 6mm, epsr = 1.05, tandel = 0.02)
        pecsheet()
        Layer()
    ]

    return strata, collect(range(0.5, 20.0, length = nfreq)), (θ = 0, ϕ = 0)
end

semi_adaptive_n(options, fmax) = ceil(Int, 15 * options.l * fmax*10e9 / (options.p * 3E8))


function run(; nfreq = 196, square_n = 21, disk_ntri = 800, seed = 1, run_direct = false)

    options = SweepOptions(
      use_D = false,
      data_partition = true,
      adaptive = false,
      tol = 1e-3,
      memory = 3,
      q1 = 8,
      q2 = 12,
      parallel = true,
      l = 2e-3
  )
    options = DEFAULT_SWEEP_OPTIONS

    strata, freqs, steering = flexible_absorber_case(; nfreq, square_n, disk_ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Flexible absorber",
        options,
        seed,
        run_direct,
    )
end

results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
#results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "flexible_absorber.png"))
