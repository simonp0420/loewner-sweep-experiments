include(joinpath(@__DIR__, "common.jl"))

function symmetric_strip_case(; nfreq = 49, ny = 60)
    c = 11.802852677165355 # inch*GHz
    period = c
    py = period
    ly = period / 2
    px = lx = ly / 10
    nx = round(Int, ny * lx / ly)

    sheet = rectstrip(; Px = px, Py = py, Lx = lx, Ly = ly, Nx = nx, Ny = ny, units = inch)
    return [Layer(), sheet, Layer()], collect(range(0.02, 0.98, length = nfreq)), (θ = 0, ϕ = 0)
end

function run(; nfreq = 49, ny = 60, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = symmetric_strip_case(; nfreq, ny)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Symmetric strip grating",
        options,
        seed,
        run_direct,
    )
end

results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
#results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "symmetric_strip.png"))
