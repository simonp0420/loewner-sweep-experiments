include(joinpath(@__DIR__, "common.jl"))

function resistive_square_patch_case(; R = 30.0, nfreq = 119)
    patch = rectstrip(;
        Nx = 10,
        Ny = 10,
        Px = 1,
        Py = 1,
        Lx = 0.5,
        Ly = 0.5,
        units = cm,
        Zsheet = R,
    )
    return [Layer(), patch, Layer()], collect(range(1.0, 60.0, length = nfreq)), (ϕ = 0, θ = 0)
end

function run(; R = 30.0, nfreq = 119, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    options = SweepOptions(
      p = 2,
      use_D = false,
      data_partition = true,
      adaptive = false,
      tol = 1e-4,
      memory = 4,
      q1 = 8,
      q2 = 12,
      parallel = true,
      l = 3e-4
  )
    strata, freqs, steering = resistive_square_patch_case(; R, nfreq)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Resistive square patch R=$(R) ohm",
        options,
        seed,
        run_direct,
    )
end

results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "resistive_square_patch.png"))

