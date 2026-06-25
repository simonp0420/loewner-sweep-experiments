include(joinpath(@__DIR__, "common.jl"))

function cross_substrate_case(;
    epsr = 4.0,
    fmin = 1.0,
    fmax = epsr == 1 ? 30.0 : epsr == 2 ? 26.0 : 20.0,
    nfreq = 301,
    ntri = 600,
)
    sheet = loadedcross(;
        w = 1.0,
        L1 = 0.6875,
        L2 = 0.0625,
        s1 = [1.0, 0.0],
        s2 = [0.0, 1.0],
        ntri,
        units = cm,
    )

    strata = [
        Layer()
        sheet
        Layer(ϵᵣ = epsr, width = 3mm)
        Layer()
    ]

    return strata, collect(range(fmin, fmax, length = nfreq)), (ϕ = 0, θ = 0)
end

function run(; epsr = 4.0, nfreq = 301, ntri = 600, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = cross_substrate_case(; epsr, nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Cross substrate epsr=$(epsr)",
        options,
        seed,
        run_direct,
    )
end

results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
#results = run();
plot_comparison(results; savepath = joinpath(@__DIR__, "cross_substrate.png"))
