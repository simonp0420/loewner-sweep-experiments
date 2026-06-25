include(joinpath(@__DIR__, "common.jl"))

function band_pass_filter_case(; nfreq = 191, ntri = 1200)
    sheet = loadedcross(;
        class = 'M',
        w = 0.023,
        L1 = 0.8,
        L2 = 0.14,
        s1 = [0.861, 0.0],
        s2 = [0.0, 0.861],
        ntri,
        units = cm,
    )

    strata = [
        Layer()
        Layer(ϵᵣ = 1.3, width = 1.1cm)
        sheet
        Layer(ϵᵣ = 1.9, width = 0.6cm)
        sheet
        Layer(ϵᵣ = 1.3, width = 1.1cm)
        Layer()
    ]

    return strata, collect(range(1.0, 20.0, length = nfreq)), (ϕ = 0, θ = 0)
end

function run(; nfreq = 191, ntri = 1200, seed = 1, run_direct = false, options = DEFAULT_SWEEP_OPTIONS)
    strata, freqs, steering = band_pass_filter_case(; nfreq, ntri)
    return run_sweep_comparison(
        strata,
        freqs,
        steering;
        label = "Loaded cross band-pass filter",
        options,
        seed,
        run_direct,
    )
end


results = run(run_direct = haskey(ENV, "RUN_DIRECT"));
plot_comparison(results; i = 3, j = 1, savepath = joinpath(@__DIR__, "band_pass_filter.png"))
